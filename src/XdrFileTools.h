/* XdrFileTools.h
 *
 * ABOUT: Object with routines for working with .xtc (and .gro) files
 * DEVELOPMENT:
 */

#pragma once

#ifndef XDR_FILE_TOOLS_H
#define XDR_FILE_TOOLS_H

// Check whether PLUMED is defined using one of its preprocessor variables
#ifdef __PLUMED_HAS_MPI
#define XDR_FILE_TOOLS_PLUMED_MODE
#endif

// Standard headers
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

// Library headers
#ifndef XDR_FILE_TOOLS_PLUMED_MODE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment" /* ignore useless warnings from this header */
#include "xdrfile.h"
#include "xdrfile_xtc.h" /* read xtc files */
#include "xdrfile_trr.h" /* read trr files */
#pragma GCC diagnostic pop
#endif

// Project headers
#include "CommonTypes.h"
#include "StringTools.h"

class XdrFileTools
{
 public:
	static const int DIM_ = CommonTypes::DIM_; // Dimensionality of simulation	

	using Real3 = CommonTypes::Real3;
	using Matrix = CommonTypes::Matrix;
	using Rvec = CommonTypes::Rvec;
	using RvecArray = CommonTypes::RvecArray;

	XdrFileTools();
	~XdrFileTools();

#ifndef XDR_FILE_TOOLS_PLUMED_MODE
	//----- Gro Files -----//

	// Read a single frame from a .gro file
	void readFrameFromGroFile(
		const std::string&        groFile,        // path to .gro file
		// Output
		RvecArray& coords,                    // atom coordinates: (x,y,z) [nm]
		std::vector<std::string>& atomTypes, // Atom types (includes index in molecule)
		std::vector<int>& atomSerialNumbers, // Atom serial no. (start at 1) read directly from file
		Rvec boxL                            // Simulation box lengths [nm]
	);

	// Print files with atom coordinates and order parameter derivatives
	void printFrameToGroFile(
		const std::string&              gro_file_name,
		const std::string&              header,
		const std::vector<int>&         residue_numbers, 
		const std::vector<std::string>& residue_names, 
		const std::vector<std::string>& atom_names,
		const std::vector<int>&         atom_numbers,  // indexed from 1
		const std::vector<Real3>&       atom_positions, 
		const Matrix&                   box_matrix
	) const;


	// Read a single frame from a .gro file
	void readFrameFromXyzFile(
		const std::string&        xyzFile, 
		// Output
		RvecArray& coords,        // Pointer to atom coordinates: (x,y,z) [nm]
		std::vector<std::string>& atomTypes, // Atom types (includes index in molecule)
		std::vector<int>& atomSerialNumbers, // Atom serial no. (start at 1) read directly from file
		Rvec boxL                            // Simulation box lengths [nm]
	);

	// Compare two trr files
	// TODO currently only checks forces and positions
	//      find a use for tol
	void compare_trr_files(
		const std::string& file_1, 
		const std::string& file_2, 
		const double tol,    // absolute tolerance on norm(delta_f) (kJ/mol/nm)
		const double rel_tol // relative tolerance norm(delta_f) [rel. to norm(f1)]
	);

#endif /* XDR_FILE_TOOLS_PLUMED_MODE */

 protected:
	// Static constants
	static constexpr double PI = 3.14159265358979323846;
	static constexpr double TWO_PI = 2.0*PI;

	// Norm of a C-style array
	template<typename T, std::size_t N>
	T norm(const T vec[N]) const 
	{
		T norm_sq = 0;
		for ( unsigned i=0; i<N; ++i ) {
			norm_sq += vec[i]*vec[i];
		}
		return sqrt( norm_sq );
	}

	// Convert array from type T to type U
	template<typename T, typename U, std::size_t N>
	void convert_array_type(const T vec[N], U vec_out[N]) const 
	{
		for ( unsigned i=0; i<N; ++i ) {
			vec_out[i] = static_cast<U>( vec[i] );
		}
	}

	// Distance between two points, excluding one dimension
	template<typename T, std::size_t N>
	T dist_dminus1(const T x_1[N], const T x_2[N], const unsigned exclude_dim) const
	{
		T norm_sq = 0, dx;
		for ( unsigned i=0; i<N; ++i ) {
			if ( i != exclude_dim ) {
				dx = x_2[i] - x_1[i];
				norm_sq += dx*dx;
			}
		}
		return sqrt( norm_sq );
	}

	template<typename T>
	T h_v_cylinder(
		const T x[DIM], 
		const T x_base[DIM], const T x_top[DIM], const T r, 
		const unsigned axis  // dimensional index for cylinder axis
	) const
	{
		// Neglecting possiblity of needing PBCs
		T radial_distance = dist_dminus1<T,DIM>( x, x_base, axis );

		T h_v;
		if ( radial_distance <= r and
		     x[axis] >= x_base[axis] and x[axis] <= x_top[axis]
		) {
			h_v = 1.0;
		}
		else {
			h_v = 0.0;
		}

		return h_v;
	}

 private:
		StringTools stringTools_;
}; // XdrFileTools

#endif // XDR_FILE_TOOLS_H
