/* XdrFileTools.cpp
 *
 */

// Standard headers

// Project headers
#include "XdrFileTools.h"

// Constructor
XdrFileTools::XdrFileTools() 
 : stringTools_()
{ }


// Destructor
XdrFileTools::~XdrFileTools() { }


#ifndef XDR_FILE_TOOLS_PLUMED_MODE
// Read a single frame from a .gro file 
void XdrFileTools::readFrameFromGroFile(
		const std::string& groFile,
		// Output
		RvecArray& coords, std::vector<std::string>& atomTypes, 
		std::vector<int>& atomSerials, XdrFileTools::Rvec boxL)
{
	// Get the number of atoms for a consistency check vs. the xtc file
	std::ifstream      ifs(groFile);
	if ( ! ifs.is_open() ) {
		std::cerr << "XdrFileTools::readFrameFromGroFile - Unable to open .gro file.\n"
							<< "(input: " << groFile << ").\n";
		exit(1);
	}

	std::string        line;
	std::istringstream ss;

	int numAtoms = -1;
	getline(ifs, line);		// Ignore first line
	getline(ifs, line);		// Second line contains the number of atoms 
	ss.str(line);
	ss >> numAtoms;

	// Checks
	if ( numAtoms < 0 ) {
		std::cout << "  XdrFileTools::readFrameFromGroFile - Invalid number of atoms "
                  << "read from .gro file." << "\n";
		exit(1);
	}

	// Allocate memory
	coords.resize(numAtoms);
	atomTypes.resize(numAtoms);
	atomSerials.resize(numAtoms);

	// Continue reading
	int               atomCounter = 0, moleculeIndex;
	std::string       moleculeName, atomName;
	std::stringstream parsingBuffer;

	while ( (atomCounter<numAtoms) && getline(ifs, line) ) {
		// Parse line
		ss.str(line);

		// GROMACS FORMAT: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
		// - Number of characters per field: 5, 5, 5, 5, 8, 8, 8, 8, 8, 8
		// - Use a stringstream to neatly handle the different variable types

		// Clear the parsing buffer
		parsingBuffer.str( std::string() );
		parsingBuffer.clear();

		// Field 1: molecule index (0-4, 5 char)
		parsingBuffer << line.substr(0, 5) << "\t";

		// Field 2: molecule name (5-9, 5 char)
		parsingBuffer << line.substr(5, 5) << "\t";

		// Field 3: atom name (10-14, 5 char)
		parsingBuffer << line.substr(10, 5) << "\t";

		// Field 4: atom index (15-19, 5 char)
		parsingBuffer << line.substr(15, 5) << "\t";

		// Field 5: x-position (20-27, 8 char)
		parsingBuffer << line.substr(20, 8) << "\t";

		// Field 6: y-position (28-35, 8 char)
		parsingBuffer << line.substr(28, 8) << "\t";

		// Field 7: z-position (36-43, 8 char)
		parsingBuffer << line.substr(36, 8) << "\t";

		// Store variables
		parsingBuffer >> moleculeIndex;
		parsingBuffer >> moleculeName;
		parsingBuffer >> atomName;
		parsingBuffer >> atomSerials[atomCounter];
		for ( int d=0; d<DIM_; ++d ) {
			parsingBuffer >> coords[atomCounter][d];
		}
		atomTypes[atomCounter] = atomName;

		atomCounter++;
	}

	//----- Last line has the box lengths -----//
	getline(ifs, line);

	// Clear flags before using with a new line
	ss.clear();
	ss.str(line);

	// Clear the parsing buffer
	parsingBuffer.str("");
	parsingBuffer.clear();

	for ( int i=0; i<3; ++i ) {
		ss >> boxL[i];
	}

	// Close .gro file
	ifs.close();

	// Checks 
	for ( int i=0; i<3; ++i ) {
		if ( boxL[i] <= 0.0 ) {
			std::cerr << "XdrFileTools::readFrameFromGroFile - Invalid box length for dimension " 
                      << i+1 << " of 3 (boxL = " << boxL[i] << ")." << "\n";
			exit(1);
		}
	}

	return;
}


// Print files with atom coordinates and order parameter derivatives
void XdrFileTools::printFrameToGroFile(
		const std::string& gro_file_name, const std::string& header,
		const std::vector<int>& residue_numbers, const std::vector<std::string>& residue_names, 
		const std::vector<std::string>& atom_names, const std::vector<int>& atom_numbers,
		const std::vector<Real3>& atom_positions, const Matrix& box_matrix
) const 
{
	// File
	std::ofstream ofs( gro_file_name );

	// Header
	int num_atoms = atom_positions.size();
	ofs << header << "\n"
	    << num_atoms << "\n";

	// Notes
	// - GRO format: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
	// - All fixed-point numbers
	ofs << std::fixed;

	for ( int i=0; i<num_atoms; ++i ) {
		// Residue number (%5d)
		ofs << std::setw(5) << residue_numbers[i];

		// Residue name (%-5s, left-justified)
		ofs << std::left << std::setw(5) << residue_names[i]
				<< std::right;

		// Atom name (%5s)
		ofs << std::setw(5) << atom_names[i];
	
		// Atom number (%5d)
		ofs << std::setw(5) << atom_numbers[i];

		// Atom position (%8.3f for each component)
		for ( int d=0; d<DIM_; ++d ) {
			ofs << std::setw(8) << std::setprecision(3) << atom_positions[i][d];
		}

		ofs << "\n";
	}

	// End with box lengths (assumes orthogonal box)
	for ( int d=0; d<DIM_; ++d ) {
		ofs << "  " << box_matrix[d][d];
		for ( int e=0; e<DIM_; ++e ) {
			if ( d != e and box_matrix[d][e] != 0.0 ) {
				std::stringstream err_ss;
				err_ss << "Error in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")"
				       << " non-orthorhombic boxes are not currently supported";
				throw std::runtime_error( err_ss.str() );
			}
		}
	}
	ofs << "\n";

	ofs.close();
}





// Read a single frame from a .xyz file
void XdrFileTools::readFrameFromXyzFile(
		const std::string& xyzFile,
		// Output
		RvecArray& coords, std::vector<std::string>& atomTypes, 
		std::vector<int>& atomSerials, XdrFileTools::Rvec boxL)
{
	// Get the number of atoms for a consistency check vs. the xtc file
	std::ifstream      ifs(xyzFile);
	if ( ! ifs.is_open() ) {
		std::cerr << "XdrFileTools::readFrameFromXyzFile - Unable to open .xyz file.\n"
							<< "(input: " << xyzFile << ").\n";
		exit(1);
	}

	std::string        line;
	std::istringstream ss;

	int numAtoms = -1;
	getline(ifs, line);		// Ignore first line
	getline(ifs, line);		// Second line contains the number of atoms 
	ss.str(line);
	ss >> numAtoms;

	// Checks
	if ( numAtoms < 0 ) {
		std::cout << "  XdrFileTools::readFrameFromXyzFile - Invalid number of atoms "
                  << "read from .xyz file." << "\n";
		exit(1);
	}

	// Allocate memory
	coords.resize(numAtoms);
	atomTypes.resize(numAtoms);
	atomSerials.resize(numAtoms);

	// Continue reading
	int               atomCounter = 0;
	std::string       moleculeName, atomName;
	std::stringstream parsingBuffer;

	while ( (atomCounter<numAtoms) && getline(ifs, line) ) {
		// Parse line
		ss.clear();
		ss.str(line);

		ss >> atomTypes[atomCounter] >>  atomSerials[atomCounter];

		for ( int d=0; d<DIM_; ++d ) {
			ss >> coords[atomCounter][d];
		}

		atomCounter++;
	}

	// Last line has the box lengths 

	getline(ifs, line);
	ss.clear();
	ss.str(line);

	for ( int i=0; i<3; ++i ) {
		ss >> boxL[i];
	}

	// Cleanup
	ifs.close();

	// Checks 
	for ( int i=0; i<3; ++i ) {
		if ( boxL[i] <= 0.0 ) {
			std::cerr << "XdrFileTools::readFrameFromXyzFile - Invalid box length for dimension " 
                      << i+1 << " of 3 (boxL = " << boxL[i] << ")." << "\n";
			exit(1);
		}
	}

	return;
}

void XdrFileTools::compare_trr_files(
		const std::string& file_1, const std::string& file_2, 
		const double tol, const double rel_tol) 
{
	/*
	const std::map<int, std::string> axis_index_to_axis_label = {
		{0, "x"},
		{1, "y"},
		{2, "z"},
	};
	*/

	// Get number of atoms
	int num_atoms_1, num_atoms_2;
	std::string file_1_name_copy = file_1;
	std::string file_2_name_copy = file_2;
	if ( exdrOK != read_trr_natoms(const_cast<char*>(file_1_name_copy.c_str()), &num_atoms_1) ) {
		throw std::runtime_error("error reading natoms from file 1");
	}
	if ( exdrOK != read_trr_natoms(const_cast<char*>(file_2_name_copy.c_str()), &num_atoms_2) ) {
		throw std::runtime_error("error reading natoms from file 1");
	}
	if ( num_atoms_1 != num_atoms_2 ) {
		throw std::runtime_error("num atoms mismatch");
	}
	int num_atoms = num_atoms_1;
	if ( num_atoms < 1 ) {
		throw std::runtime_error("no atoms in trr files");
	}

	// User feedback
	std::cout << "MDAnalysis++: Comparing trr files\n"
	          << "  file_1 = " << file_1 << "\n"
	          << "  file_2 = " << file_2 << "\n"
	          << "  tol    = " << tol    << "\n"
	          << "  natoms = " << num_atoms << "\n";
	
	// Prepare output file
	std::string   output_file_name = "force_mismatch.out";
	std::ofstream ofs(output_file_name);
	ofs << "# Force mismatch log\n"
	    << "#   tol     = " << tol << " [kJ/(mol*nm)]\n"
	    << "#   rel_tol = " << rel_tol << "\n"
	    << "\n"
			<< "# atom_serial  x(1)    x(2)    f(1)   f(2)   "
	        << "norm(df)  norm(df)/norm[f(1)]  norm(df)_exceeds_tol?\n";
			//<< "# atom_serial   x(1)    x(2)   f(1)   f(2)   norm(df)/df_1\n"; // FIXME DEBUG

	// Allocate memory
	XdrFileTools::RvecArray x_1(num_atoms), v_1(num_atoms), f_1(num_atoms);
	XdrFileTools::RvecArray x_2(num_atoms), v_2(num_atoms), f_2(num_atoms);

	// Working variables
	float  lambda_1, lambda_2, time_1, time_2;
	int    step_1, step_2;
	matrix box_1, box_2;
	int    default_field_width = 8;

	// Open files
	XDRFILE* file_1_ptr = xdrfile_open(file_1.c_str(), "r");
	XDRFILE* file_2_ptr = xdrfile_open(file_2.c_str(), "r");
	if ( file_1_ptr == nullptr or file_2_ptr == nullptr ) {
		throw std::runtime_error("unable to open file");
	}

	// Read through frames
	while ( ( exdrOK == read_trr(file_1_ptr, num_atoms, 
	                             &step_1, &time_1, &lambda_1, box_1,
	                             x_1.data(), v_1.data(), f_1.data()) )
	        and
	        ( exdrOK == read_trr(file_2_ptr, num_atoms, 
	                             &step_2, &time_2, &lambda_2, box_2,
	                             x_2.data(), v_2.data(), f_2.data()) )
	) {
		if ( time_1 != time_2 ) {
			throw std::runtime_error("time mismatch");
		}

		// Record the time
		ofs << "# t= " << time_1 << " [ps]\n";

		for ( int i=0; i<num_atoms; ++i ) {
			// Compare positions and forces
			bool are_different = false;
			for ( int d=0; d<DIM_; ++d ) {
				if ( /* x_1[i][d] != x_2[i][d] or */ f_1[i][d] != f_2[i][d] ) { // FIXME DEBUG
					are_different = true;
					break;
				}
			}

			if ( are_different ) {
				// Norm of the difference between forces
				double df, norm_df = 0.0, norm_f_1 = 0.0;
				for ( int d=0; d<DIM_; ++d ) {
					df = f_2[i][d] - f_1[i][d];
					norm_df  += df*df;
					norm_f_1 += f_1[i][d]*f_1[i][d];
				}
				norm_df  = sqrt(norm_df);
				norm_f_1 = sqrt(norm_f_1);
				double rel_df = norm_df/norm_f_1;

				int num_digits = 2;

				// Print indices, positions, and forces for the offending atoms
				ofs << std::fixed << std::setw(default_field_width) << i+1;
				for ( int d=0; d<DIM_; ++d ) {
					// Position(1)
					ofs << " " << std::fixed << std::setprecision(num_digits) << std::setw(default_field_width) << x_1[i][d];
				}
				for ( int d=0; d<DIM_; ++d ) {
					// Position(2)
					ofs << " " << std::fixed << std::setprecision(num_digits) << std::setw(default_field_width) << x_2[i][d];
				}
				for ( int d=0; d<DIM_; ++d ) {
					// Force(1)
					ofs << " " << std::scientific << std::setprecision(num_digits) << std::setw(10) << f_1[i][d];
				}
				for ( int d=0; d<DIM_; ++d ) {
					// Force(2)
					ofs << " " << std::scientific << std::setprecision(num_digits) << std::setw(10) << f_2[i][d];
				}
				ofs << " " << std::scientific << std::setprecision(num_digits) << std::setw(10) << norm_df 
				    << " " << std::scientific << std::setprecision(num_digits) << std::setw(10) << rel_df
				    << " " << std::setw(4) << (norm_df > tol)
				    << " " << std::setw(4) << (rel_df > rel_tol)
				    << "\n";
			} // end if ( are_different )

			/*
			// FIXME DEBUG: Examine frames which disagree in close detail
			if ( i % 4 == 0 and (time_1 == 192.0 or time_1 == 4020.0) ) {
				// Cylinder (double)
				double r_max_d = 1.0;
				double x_base_d[DIM], x_top_d[DIM];
				x_base_d[0] = 2.0;   x_top_d[0] = 3.0;
				x_base_d[1] = 2.5;   x_top_d[1] = 2.5;
				x_base_d[2] = 2.5;   x_top_d[2] = 2.5;
				const unsigned axis = 0;  // cylinder along x-axis

				// Cylinder (float)
				float r_max_f = static_cast<float>( r_max_d );
				float x_base_f[DIM], x_top_f[DIM];
				for ( int d=0; d<DIM; ++d ) {
					x_base_f[d] = static_cast<float>( x_base_d[d] );
					x_top_f[d]  = static_cast<float>( x_top_d[d] );
				}

				// GMX INDUS (file 1)
				float h_v_1 = h_v_cylinder<float>(x_1[i], x_base_f, x_top_f, r_max_f, axis );

				// PLMD MDA++ (file 2)
				double x_1_d[DIM];
				convert_array_type<float,double,DIM>( x_1[i], x_1_d );  // switch to double
				double h_v_2 = h_v_cylinder<double>(x_1_d, x_base_d, x_top_d, r_max_d, axis );

				// Output
				if ( h_v_1 != 0.0 or h_v_2 != 0.0 ) {
					std::cout << "  " << i+1 << "  " << h_v_1 << "  " << h_v_2 << "\n";
				}
			}
			*/
		} // end loop over atoms
	} // end loop over frames

	ofs.close();
}

#endif /* XDR_FILE_TOOLS_PLUMED_MODE */
