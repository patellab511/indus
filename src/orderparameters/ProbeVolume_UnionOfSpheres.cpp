/* ProbeVolume_UnionOfSpheres.cpp
 *
 * ABOUT: see header
 */

#include "ProbeVolume_UnionOfSpheres.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

// Register this probe volume under two names
namespace ProbeVolumeRegistry {
static const Register<ProbeVolume_UnionOfSpheres>
	register_ProbeVolume_UnionOfSpheres("union_of_spheres");
static const Register<ProbeVolume_UnionOfSpheres>
	register_ProbeVolume_UnionOfSphericalShells("union_of_spherical_shells");
}


ProbeVolume_UnionOfSpheres::ProbeVolume_UnionOfSpheres(
	ProbeVolume::ProbeVolumeInputPack& input_pack
):
	ProbeVolume(input_pack),
	r_min_(-1.0),
	r_max_(1.0)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Radii (shared by all spheres)
	input_parameter_pack.readNumber("r_max", KeyType::Required, r_max_);
	input_parameter_pack.readNumber("r_min", KeyType::Optional, r_min_);
	if ( r_max_ <= 0.0 ) {
		throw std::runtime_error("union_of_spheres: r_max must be positive");
	}
	if ( r_min_ >= r_max_ ) {
		throw std::runtime_error("union_of_spheres: r_min must be less than r_max");
	}

	// Sphere centers
	input_parameter_pack.readString("reference_positions_file", KeyType::Required,
	                                reference_positions_file_);
	sphere_centers_ = readSphereCenters(reference_positions_file_);
	if ( sphere_centers_.empty() ) {
		throw std::runtime_error("union_of_spheres: no sphere centers were found in file \""
		                         + reference_positions_file_ + "\"");
	}

	// Build subvolumes
	// - The sphere constructor only uses 'input_pack' to set base-class options
	//   (sigma, alpha_c), which setGeometry() overwrites anyway
	const int num_subvolumes = sphere_centers_.size();
	subvolume_ptrs_.reserve(num_subvolumes);
	for ( int k=0; k<num_subvolumes; ++k ) {
		subvolume_ptrs_.push_back(
			SubvolumePtr( new ProbeVolume_Sphere(sphere_centers_[k], r_min_, r_max_, input_pack) )
		);
	}

	setGeometry();

	// One scratch buffer per OpenMP thread
	subvolumes_buffers_.resize( OpenMP::get_max_threads() );
}


void ProbeVolume_UnionOfSpheres::setGeometry()
{
	const int num_subvolumes = subvolume_ptrs_.size();
	for ( int k=0; k<num_subvolumes; ++k ) {
		ProbeVolume_Sphere& sphere = *subvolume_ptrs_[k];
		sphere.setShellWidths(width_shell_1_, width_shell_2_);
		sphere.setCoarseGrainingParameters(sigma_, alpha_c_);
		sphere.setGeometry(sphere_centers_[k], r_min_, r_max_);
	}

	// All spheres share the same radii, so one is representative
	rtilde_max_ = ( num_subvolumes > 0 ) ? subvolume_ptrs_[0]->get_rtilde_max() : 0.0;
}


std::vector<ProbeVolume_UnionOfSpheres::Real3> ProbeVolume_UnionOfSpheres::readSphereCenters(
	const std::string& file
)
{
	std::ifstream ifs(file);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("union_of_spheres: unable to open reference_positions_file \""
		                         + file + "\"");
	}

	std::vector<Real3> centers;
	Real3 position;
	int dim = 0;
	int line_number = 0;

	std::string line, token;
	while ( std::getline(ifs, line) ) {
		++line_number;

		// Strip comments
		const auto comment_start = line.find('#');
		if ( comment_start != std::string::npos ) {
			line.erase(comment_start);
		}

		// Values may be split across lines arbitrarily, but every 3 form one center
		std::stringstream ss(line);
		while ( ss >> token ) {
			try {
				position[dim] = std::stod(token);
			}
			catch ( const std::exception& ) {
				std::stringstream err_ss;
				err_ss << "union_of_spheres: could not parse \"" << token << "\" as a number"
				       << " (file " << file << ", line " << line_number << ")";
				throw std::runtime_error( err_ss.str() );
			}
			++dim;
			if ( dim == DIM_ ) {
				centers.push_back(position);
				dim = 0;
			}
		}
	}

	if ( dim != 0 ) {
		std::stringstream err_ss;
		err_ss << "union_of_spheres: reference_positions_file \"" << file << "\" must contain"
		       << " 3 values (x y z) per sphere center; found " << dim << " trailing value(s)";
		throw std::runtime_error( err_ss.str() );
	}

	return centers;
}


BoundingBox ProbeVolume_UnionOfSpheres::constructBoundingBox() const
{
	// *** Assumes orthorhombic box *** //
	// BoundingBox extends itself to the full box along any axis where the
	// corners cross the periodic boundaries

	Real3 x_lower, x_upper;
	x_lower.fill( std::numeric_limits<double>::max() );
	x_upper.fill( std::numeric_limits<double>::lowest() );
	for ( const Real3& center : sphere_centers_ ) {
		for ( int d=0; d<DIM_; ++d ) {
			x_lower[d] = std::min(x_lower[d], center[d]);
			x_upper[d] = std::max(x_upper[d], center[d]);
		}
	}

	const double dr_buffer = rtilde_max_ + bounding_box_tol_;
	for ( int d=0; d<DIM_; ++d ) {
		x_lower[d] -= dr_buffer;
		x_upper[d] += dr_buffer;
	}

	return BoundingBox(x_lower, x_upper, simulation_box_);
}


void ProbeVolume_UnionOfSpheres::calculateIndicator(
	const Real3& x,
	double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
) const
{
	// Defaults: outside everything
	h_v      = 0.0;
	htilde_v = 0.0;
	dhtilde_v_dx.fill(0.0);
	region = RegionEnum::Unimportant;

	// Thread-local scratch space
	const int thread_id = OpenMP::get_thread_num();
	if ( thread_id >= static_cast<int>(subvolumes_buffers_.size()) ) {
		throw std::runtime_error("union_of_spheres: more OpenMP threads than buffers");
	}
	SubvolumesBuffer& buffer = subvolumes_buffers_[thread_id];
	buffer.clear();

	// Scan subvolumes, accumulating prod_k (1 - htilde_k)
	double product = 1.0;
	double h_k, htilde_k;
	Real3  dhtilde_k_dx;
	RegionEnum region_k;
	const int num_subvolumes = subvolume_ptrs_.size();
	for ( int k=0; k<num_subvolumes; ++k ) {
		subvolume_ptrs_[k]->calculateIndicator(x, h_k, htilde_k, dhtilde_k_dx, region_k);

		// Any sphere that puts the atom in a stronger region wins
		region = strongerRegion(region, region_k);

		if ( htilde_k <= 0.0 ) {
			continue;
		}

		if ( h_k == 1.0 ) {
			h_v = 1.0;
		}

		if ( htilde_k >= 1.0 ) {
			// Fully inside this sphere: htilde_v = 1 and every derivative vanishes,
			// so the remaining spheres cannot change anything
			product = 0.0;
			break;
		}

		product *= (1.0 - htilde_k);
		buffer.htilde.push_back(htilde_k);
		if ( need_derivatives_ ) {
			buffer.dhtilde_dx.push_back(dhtilde_k_dx);
		}
	}

	// Finish the indicator (guard against round-off from many near-one factors)
	if ( product <= 0.0 ) {
		htilde_v = 1.0;
	}
	else if ( product >= 1.0 ) {
		htilde_v = 0.0;
	}
	else {
		htilde_v = 1.0 - product;
	}

	if ( htilde_v > 0.0 ) {
		region = RegionEnum::Vtilde;
	}

	// Derivatives: nonzero only in the buffer region of the union
	if ( need_derivatives_ and htilde_v > 0.0 and htilde_v < 1.0 ) {
		// d htilde_v / dx = sum_k [ prod_{l != k} (1 - htilde_l) ] d htilde_k / dx
		// - The number of active spheres is small, so the O(n^2) product is cheap
		//   and avoids dividing by (1 - htilde_k), which can be tiny
		const int num_active = buffer.htilde.size();
		for ( int k=0; k<num_active; ++k ) {
			double product_others = 1.0;
			for ( int l=0; l<num_active; ++l ) {
				if ( l != k ) {
					product_others *= (1.0 - buffer.htilde[l]);
				}
			}
			for ( int d=0; d<DIM_; ++d ) {
				dhtilde_v_dx[d] += product_others * buffer.dhtilde_dx[k][d];
			}
		}
	}
}


std::string ProbeVolume_UnionOfSpheres::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Probe_volume union_of_spheres\n"
	   << prepend_string << "  reference_positions_file = " << reference_positions_file_ << "\n"
	   << prepend_string << "  num_spheres = " << subvolume_ptrs_.size() << "\n"
	   << prepend_string << "  r_min = " << r_min_ << " [nm]\n"
	   << prepend_string << "  r_max = " << r_max_ << " [nm]\n";

	// List the centers, but keep the header short for large unions
	const int max_centers_to_print = 10;
	const int num_centers = sphere_centers_.size();
	ss << prepend_string << "  centers [nm]\n";
	for ( int k=0; k<std::min(num_centers, max_centers_to_print); ++k ) {
		const Real3& c = sphere_centers_[k];
		ss << prepend_string << "    " << k+1 << ": {" << c[0] << ", " << c[1] << ", " << c[2] << "}\n";
	}
	if ( num_centers > max_centers_to_print ) {
		ss << prepend_string << "    ... (" << num_centers - max_centers_to_print << " more)\n";
	}

	ss << getSharedAttributesSummary(prepend_string + "  ");

	return ss.str();
}
