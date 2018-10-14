/* Topology.cpp
 *
 */

#include "Topology.h"


Topology::Topology()
{
	this->clear(); // probably redundant
}


Topology::Topology(const std::string& input_file)
 : top_file_(""), gro_file_(""),
   num_atoms_(0), num_molecules_(0)
{
	// Deduce whether the file is .top or .gro from its file extension
	setFromFile(input_file);
}


Topology::Topology(const std::string& top_file, const std::string& gro_file)
 : top_file_(top_file), gro_file_(gro_file),
   num_atoms_(0), num_molecules_(0)
{
	setFromFile(top_file_, gro_file_);
}


void Topology::setFromFile(const std::string& input_file)
{
	this->clear();

	// Infer the type of file provided based on its extension
	std::string extension = string_tools_.getFileExtension(input_file);

	if ( extension == "top" ) {
		top_file_ = input_file;
		readGromacsTopology(top_file_);
	}
	else if ( extension == "gro" ) {
		gro_file_ = input_file;

		std::vector<std::string> residue_names;
		readGroFile(gro_file_, 
		            residue_names, atom_names_, 
		            reference_positions_, reference_box_matrix_);

		num_atoms_ = reference_positions_.size();
	}
	else if ( extension.empty() ) {
		std::stringstream err_ss;
		err_ss << "Error in Topology constructor: Could not determine type of topology file.\n";
		throw std::runtime_error( err_ss.str() );
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in Topology constructor: File type \"" << extension << "\" is invalid.\n";
		throw std::runtime_error( err_ss.str() );
	}
}


void Topology::setFromFile(const std::string& top_file, const std::string& gro_file)
{
	this->clear();

	top_file_ = top_file;
	gro_file_ = gro_file;

	// Read top file
	readGromacsTopology(top_file_);

	// Read gro file, but only store the positions and box info
	readReferenceStateFromGroFile(gro_file_);
}


void Topology::readGromacsTopology(const std::string& top_file)
{
	this->clear();

	// Check input
	std::string extension = string_tools_.getFileExtension(top_file);
	if ( extension != "top" ) {
		std::stringstream err_ss;
		err_ss << "Error in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")\n"
		       << "   Input file \"" << top_file << "\" does not appear to be a GROMACS .top file.\n";
		throw std::runtime_error( err_ss.str() );
	}

	std::ifstream ifs;
	try {
		ifs.open(top_file);
	}
	catch (std::ios_base::failure failed_to_open) {
		std::cerr << "  Topology::readGromacsTopology: Failed to open top file (exception: "
		          << failed_to_open.what() << ")." << "\n";
		throw;
	}

	// Working variables
	int line_number = 0;
	int num_molecules_total = 0;
	std::string line, trimmed_line, comment, comment_chars = ";";
	std::string record_name = "", current_molecule;
	std::vector<std::string> tokens;

	MoleculeType* molecule_type_ptr = nullptr;

	while ( getline(ifs, line) ) {
		++line_number;

		// Yank off any comment
		string_tools_.removeTrailingCommentFromLine(
			line, comment_chars, trimmed_line, comment);

		// Tokenize based on whitespace
		string_tools_.split(trimmed_line, tokens);

		// Skip blank lines and comments
		int numTokens = tokens.size();
		if ( numTokens == 0 ) {
			continue;
		}

		// Check for a new record
		if ( numTokens >= 3 and tokens[0] == "[" and tokens[2] == "]" ) {
			record_name = tokens[1];
			continue; // nothing more to do with this line
		}

		//----- [ atomtypes ] -----//

		if ( record_name == "atomtypes" and numTokens >= 7 ) {
			// Create a new entry
			std::string atom_type = tokens[0];
			atom_types_[atom_type] = AtomType();

			// Store definition
			AtomType& atom_typeRef = atom_types_[atom_type];
			atom_typeRef.name         = atom_type;
			atom_typeRef.atomic_number = std::stoi( tokens[1] );
			atom_typeRef.mass         = std::stod( tokens[2] );
			atom_typeRef.charge       = std::stod( tokens[3] );
			atom_typeRef.ptype        = tokens[4];
			atom_typeRef.sigma        = std::stod( tokens[5] );
			atom_typeRef.epsilon      = std::stod( tokens[6] );
		}

		//----- [ moleculetype ] -----//

		else if ( record_name == "moleculetype" and numTokens >= 2 ) {
			// Add the new molecule to the map and save a pointer to it so that
			// the hash table doesn't need to be repeatedly consulted
			current_molecule = tokens[0];
			molecule_types[current_molecule] = MoleculeType();
			molecule_type_ptr = &( molecule_types[current_molecule] );
			molecule_type_ptr->name  = current_molecule;
			molecule_type_ptr->nrexl = std::stoi(tokens[1]);
		}

		//----- [ atoms ] -----//

		else if ( record_name == "atoms" and numTokens >= 6 ) {
			if ( molecule_type_ptr != nullptr ) {
				std::string atom_type = tokens[1];
				molecule_type_ptr->atom_types.push_back( atom_type );
				molecule_type_ptr->residue_numbers.push_back( std::stoi(tokens[2]) );
				molecule_type_ptr->residue_types.push_back( tokens[3] );
				molecule_type_ptr->atom_names.push_back( tokens[4] );

				// Mass and charge may or may not be inline.
				// If not inline, take values from the atomtypes section.

				double charge;
				if ( numTokens >= 7 ) { charge = std::stod(tokens[6]);        }
				else                  { charge = atom_types_[atom_type].charge; }
				molecule_type_ptr->charges.push_back( charge );

				double mass;
				if ( numTokens >= 8 ) { mass = std::stod(tokens[7]);      }
				else                  { mass = atom_types_[atom_type].mass; }
				molecule_type_ptr->masses.push_back( mass );
			}
			else {
				std::stringstream err_ss;
				err_ss << "Error in Topology::parseGromacsTopology: \n"
				       << "  An [ atoms ] directive was found, but no [ moleculetype ] has been found.\n";
				throw std::runtime_error( err_ss.str() );
			}
		}

		//----- [ molecules ] -----//

		else if ( record_name == "molecules" and numTokens >= 2 ) {
			// Add new entry
			std::string molecule_name = tokens[0];
			molecules_directive_entries_.push_back( MoleculesEntry() );

			int index = static_cast<int>(molecules_directive_entries_.size()) - 1;
			molecules_directive_entries_[index].molecule_name = molecule_name;
			molecules_directive_entries_[index].count        = std::stoi(tokens[1]);

			// Track the total number of molecules
			num_molecules_total += molecules_directive_entries_[index].count;
		}
	}

	// Store the number of atoms in each molecule type explicitly, for convenience
	for ( auto it=molecule_types.begin(); it!=molecule_types.end(); ++it ) {
		const std::string& name = it->first;        // molecule name
		MoleculeType& type = molecule_types[name];  // reference to MoleculeType itself
		type.num_atoms = type.atom_types.size();      // use atom_types vector to infer # atoms
	}

	// Unwrap the [ molecules ] entries to get the atom-wise, residue-wise, 
	// and molecule-wise topologies
	molecules_.resize(num_molecules_total);
	int global_index = 0, molecule_counter = 0;
	int num_molecules_entries = molecules_directive_entries_.size();
	for ( int i=0; i<num_molecules_entries; ++i ) {
		// Find the corresponding molecule type
		std::string name   = molecules_directive_entries_[i].molecule_name;
		MoleculeType& type = molecule_types[name];

		// Enumerate all molecules of this type
		int num_molecules_of_type  = molecules_directive_entries_[i].count;
		int num_atoms_per_molecule = type.num_atoms;

		/*
		// DEBUG
		std::cout << "Molecule " << name << "\n"
		          << "  num. of type    = " << num_molecules_of_type << "\n"
		          << "  num. atoms each = " << num_atoms_per_molecule << "\n";
		*/

		for ( int j=0; j<num_molecules_of_type;
		      ++j, ++molecule_counter, global_index += num_atoms_per_molecule ) {
			// Molecule
			molecules_[molecule_counter].type        = name;
			molecules_[molecule_counter].start_index = global_index;
			molecules_[molecule_counter].num_atoms   = num_atoms_per_molecule;

			/*
			std::cout << "Molecule " << molecule_counter+1 << "\n" // DEBUG
			          << "  name  "  << name << "\n"
			          << "  start "  << global_index+1 << "\n";
			*/

			// Atoms in this molecule
			atoms_.insert(atoms_.end(), type.atom_types.begin(), type.atom_types.end());
			atom_masses_.insert(atom_masses_.end(), type.masses.begin(), type.masses.end());

			// Residues
			int last_residue_number_in_molecule = 0;
			for ( int k=0; k<num_atoms_per_molecule; ++k ) {
				if ( type.residue_numbers[k] != last_residue_number_in_molecule ) {
					// New residue encountered
					Residue new_residue;
					new_residue.type        = type.residue_types[k];
					new_residue.start_index = global_index + k;
					new_residue.num_atoms   = 0;
					new_residue.number      = type.residue_numbers[k];
					new_residue.atom_type_names.clear();
					residues_.push_back(new_residue);

					last_residue_number_in_molecule = type.residue_numbers[k];
				}

				Residue& last_residue_in_vec = residues_.back();
				++( last_residue_in_vec.num_atoms );
				last_residue_in_vec.atom_type_names.push_back( type.atom_types[k] );
			}
		}
	}

	num_atoms_     = atoms_.size();
	num_molecules_ = molecules_.size();

	/*
	// DEBUG
	int num_residues = residues_.size();
	std::cout << "RESIDUES (" << num_residues << " total)\n";
	for ( int i=0; i<num_residues; ++i ) {
		int num_atoms_in_residue = residues_[i].num_atoms;
		int first_serial = residues_[i].start_index + 1;
		int last_serial  = first_serial + num_atoms_in_residue;

		int num_heavy_atoms = 0;

		std::cout << "  " << residues_[i].number << ": " << residues_[i].type << "\n"
		          << "    starts    " << first_serial << "\n"
		          << "    ends      " << last_serial  << "\n"
		          << "    num_atoms " << num_atoms_in_residue << "\n";
	}
	*/

	ifs.close();
}


void Topology::readGroFile(
	const std::string& gro_file,
	std::vector<std::string>& residue_names, std::vector<std::string>& atom_names,
	RvecArray& positions, Matrix& box_matrix
) const
{
	// Check input
	std::string extension = string_tools_.getFileExtension(gro_file);
	if ( extension != "gro" ) {
		std::stringstream err_ss;
		err_ss << "Error in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")\n"
		       << "   Input file \"" << gro_file << "\" does not appear to be a GROMACS .gro file.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Try to open the file
	std::ifstream ifs;
	ifs.open(gro_file);
	if ( ! ifs.is_open() ) {
		std::stringstream err_ss;
		err_ss << "Topology::readGroFile - Unable to open .gro file.\n"
		       << "  Input: \"" + gro_file + "\"\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Working variables
	std::string line;
	std::istringstream ss;

	// Read header
	int num_atoms_from_gro_file = -1; // sentinel value
	getline(ifs, line);		// Ignore first line
	getline(ifs, line);		// Second line contains the number of atoms 
	ss.str(line);
	ss >> num_atoms_from_gro_file;

	// Checks
	if ( num_atoms_from_gro_file <= 0 ) {
		std::stringstream err_ss;
		err_ss << "  Topology::readGroFile - Invalid number of atoms ("
           << num_atoms_from_gro_file << ") reported by .gro file." << "\n"
		       << "  Input: " << gro_file << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Allocate memory
	residue_names.resize(num_atoms_from_gro_file);
	atom_names.resize(num_atoms_from_gro_file);
	positions.resize(num_atoms_from_gro_file);

	// Read information about atoms
	int atomCounter = 0, molecule_serial, atom_serial;
	std::stringstream parsing_buffer;
	while ( (atomCounter<num_atoms_from_gro_file) and getline(ifs, line) ) {
		// Parse line
		ss.str(line);

		// GROMACS FORMAT: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
		// - Number of characters per field: 5, 5, 5, 5, 8, 8, 8, 8, 8, 8
		// - Use a stringstream to neatly handle the different variable types

		// Clear the parsing buffer
		parsing_buffer.str( std::string() );
		parsing_buffer.clear();

		// Field 1: molecule index (0-4, 5 char)
		parsing_buffer << line.substr(0, 5) << "\t";

		// Field 2: molecule name (5-9, 5 char)
		parsing_buffer << line.substr(5, 5) << "\t";

		// Field 3: atom name (10-14, 5 char)
		parsing_buffer << line.substr(10, 5) << "\t";

		// Field 4: atom index (15-19, 5 char)
		parsing_buffer << line.substr(15, 5) << "\t";

		// Field 5: x-position (20-27, 8 char)
		parsing_buffer << line.substr(20, 8) << "\t";

		// Field 6: y-position (28-35, 8 char)
		parsing_buffer << line.substr(28, 8) << "\t";

		// Field 7: z-position (36-43, 8 char)
		parsing_buffer << line.substr(36, 8) << "\t";

		// Store variables as needed
		parsing_buffer >> molecule_serial;
		parsing_buffer >> residue_names[atomCounter];
		parsing_buffer >> atom_names[atomCounter];
		parsing_buffer >> atom_serial;
		for ( int d=0; d<DIM_; ++d ) {
			parsing_buffer >> positions[atomCounter][d];
		}

		atomCounter++;
	}

	if ( atomCounter != num_atoms_from_gro_file ) {
		ifs.close();
		throw std::runtime_error(
			"Error in Topology::readGroFile: file ended before all expected atoms were read."
		);
	}

	// Last line has the box lengths
	// (FIXME: box matrix, if triclinic? how to check if extra, unexpected lines?)
	getline(ifs, line);
	ss.clear();   // Clear flags before using with a new line
	ss.str(line);
	for ( int d=0; d<3; ++d ) {
		ss >> box_matrix[d][d];
	}

	// Wrap up
	ifs.close();

	// Checks
	for ( int d=0; d<DIM_; ++d ) {
		if ( box_matrix[d][d] <= 0.0 ) {
			std::stringstream err_ss;
			err_ss << "Topology::readGroFile - Invalid box length for dimension d=" 
						 << d+1 << " of 3 (reference_box_matrix[d][d] = " << box_matrix[d][d]
						 << ")." << "\n";
			throw std::runtime_error( err_ss.str() );
		}
	}
}


void Topology::readReferenceStateFromGroFile(const std::string& gro_file)
{
	// Read gro file and save its positions and box length as reference
	std::vector<std::string> residue_names, atom_names_from_gro_file;
	readGroFile(gro_file, 
	            residue_names, atom_names_from_gro_file, 
							reference_positions_, reference_box_matrix_);

	// If Topology information is present, validate what was read from the gro
	// file using it
	if ( not top_file_.empty() and num_atoms_ > 0 ) {
		// Check number of atoms
		int num_atoms_from_gro_file = reference_positions_.size();
		if ( num_atoms_ != 0 and num_atoms_from_gro_file != num_atoms_ ) {
			std::stringstream err_ss;
			err_ss << "  Topology::readGroFile - Stored number of atoms (" 
						 << num_atoms_ << ") and number of atoms in gro file (" 
						 << num_atoms_from_gro_file << ") do not match.\n";
			throw std::runtime_error( err_ss.str() );
		}

		// Check atom names
		if ( atom_names_.size() == num_atoms_ ) {
			for ( int i=0; i<num_atoms_; ++i ) {
				if ( atom_names_from_gro_file[i] != atom_names_[i] ) {
					std::stringstream err_ss;
					err_ss << "  Topology::readGroFile - Name of atom " << i+1 
					       << " (" << atom_names_from_gro_file[i] << ") in gro file does not match"
					       << " expected name (" << atom_names_[i] << "\n"
					       << "  Input gro file: " << gro_file << "\n";
					throw std::runtime_error( err_ss.str() );
				}
			}
		}
	}
}


// Get indices of target atoms based on string from input file
// - A "target atom" can be a virtual site, such as the center of mass of a molecule
void Topology::getTargets(
		const std::vector<std::string>& target_tokens,
		// Output
		std::vector<int>& target_indices, TargetFlags& target_flags
) const
{
	// Check input
	int num_tokens = target_tokens.size();
	if ( num_tokens < 2 ) {
		std::stringstream err_ss;
		err_ss << "Topology::getTargets - Target atoms contains fewer than 2 tokens.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// For reference, get a string recording the input
	std::string target_string;
	for ( int i=0; i<num_tokens; ++i ) {
		if ( i > 0 ) { 
			target_string += " "; 
		}
		target_string += target_tokens[i];
	}

	// Reset all flags to false
	target_flags = TargetFlags();


	//----- Figure out the parsing mode (by atom type, atom indices, etc.) -----//

	std::string token;
	std::string targetAtomType, targetAtomName, target_index_ranges;
	std::string targetMoleculeType;
	std::string fileName;

	// Ensure parsing mode token is lowercase
	std::string parsing_mode = target_tokens[0];
	std::transform(parsing_mode.begin(), parsing_mode.end(), parsing_mode.begin(), ::tolower);

	if ( parsing_mode == "atom_type" ) { 
		targetAtomType = target_tokens[1];
		target_flags.is_atom = true;
		target_flags.isType = true;
	}
	else if ( parsing_mode == "atom_name" ) { 
		targetAtomName = target_tokens[1];
		target_flags.is_atom = true;
		target_flags.isName = true;
	}
	else if ( parsing_mode == "atom_index" ) { 
		target_index_ranges = target_tokens[1];
		target_flags.is_atom  = true;
		target_flags.isIndex = true;
	}
	else if (	parsing_mode == "atom_index_com" ) {
		// Remove keyword from beginning of target string
		target_index_ranges = target_tokens[1];
		target_flags.is_atom         = true;
		target_flags.isIndex        = true;
		target_flags.isCenterOfMass = true;
	}
	else if ( parsing_mode == "atom_index_file" ) {
		fileName = target_tokens[1];
		target_flags.is_atom = true;
		target_flags.isIndex = true;
		target_flags.isFile  = true;
	}
	else if ( parsing_mode == "molecule_type_com" ) {
		// Centers of mass of molecules of the specified type
		targetMoleculeType = target_tokens[1];
		target_flags.isMolecule     = true;
		target_flags.isType         = true;
		target_flags.isCenterOfMass = true;
	}
	else if ( parsing_mode == "molecule_index_com" ) {
		// Remove keyword from beginning of target string
		target_index_ranges = target_tokens[1];
		target_flags.isMolecule     = true;
		target_flags.isIndex        = true;
		target_flags.isCenterOfMass = true;
	}
	else if ( parsing_mode == "fixed_positions_list" ) {
		target_flags.is_positions_list = true;
		if ( target_tokens[1] != "file" ) {
			std::stringstream err_ss;
			err_ss << "Error in Topology::getTargets\n"
						 << "  For selection \"" << parsing_mode << "\","
			         << " the only supported format is \"file\"\n";
			throw std::runtime_error( err_ss.str() );
		}
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in Topology::getTargets\n"
		       << "  Selection mode \"" << parsing_mode << "\" not recognized.\n"
		       << "  Note that selection modes are case-insensitive.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Determine target indices
	target_indices.clear();
	if ( target_flags.isType or target_flags.isName ) {

		//----- Target is a type: search the stored topology -----//

		if ( target_flags.is_atom ) {
			// Sanity checks
			if ( num_atoms_ == 0 ) {
				std::stringstream err_ss;
				err_ss << "Error in Topology::getTargets\n"
				       << "  Requested target atoms (input: \"" << target_string
				       << "\"), but this Topology\n  doesn't have that information.\n";
				throw std::runtime_error( err_ss.str() );
			}
			else if ( target_flags.isType and atoms_.size() != num_atoms_ ) {
				std::stringstream err_ss;
				err_ss << "Error in Topology::getTargets\n"
				       << "  Requested target atoms by type (input: \"" << target_string
				       << "\"), but this Topology\n  doesn't have that information.\n";
				throw std::runtime_error( err_ss.str() );
			}
			else if ( target_flags.isName and atom_names_.size() != num_atoms_ ) {
				std::stringstream err_ss;
				err_ss << "Error in Topology::getTargets\n"
				       << "  Requested target atoms by name (input: \"" << target_string
				       << "\"), but this Topology\n  doesn't have that information.\n";
				throw std::runtime_error( err_ss.str() );
			}
			

			for ( int i=0; i<num_atoms_; ++i ) {
				if (    (target_flags.isType and atoms_[i]     == targetAtomType)
				     or (target_flags.isName and atom_names_[i] == targetAtomName) ) {
					target_indices.push_back(i);
				}
			}
		}
		else if ( target_flags.isMolecule ) {
			int num_molecules = molecules_.size();
			if ( num_molecules == 0 ) {
				std::stringstream err_ss;
				err_ss << "Error in Topology::getTargets\n"
				       << "  Requested target molecules of type \"" << targetMoleculeType
				       << "\", but this Topology\n"
				       << "  doesn't have this information.\n";
				throw std::runtime_error( err_ss.str() );
			}

			for ( int j=0; j<num_molecules; ++j ) {
				if ( molecules_[j].type == targetMoleculeType ) {
					target_indices.push_back(j);
				}
			}
		}
	}
	else if ( target_flags.isIndex and (not target_flags.isFile) ) {
		//----- Target is a set of indices -----//

		// Format: comma-separated ranges of the form "{start}-{stop}:{stride}" (PLUMED syntax)
		//  - {stride} is optional (default = 1)
		//  - Can have multiple ranges, separated by commas

		int stride, start, stop;
		size_t pos_colon, pos_dash, next_comma = 0;
		size_t current_pos = 0; // functions as a pointer to the first character of the current range
		std::string range_string;

		bool done = false;
		while ( ! done ) {
			// Look for a comma, which indicates the end of the current range
			next_comma = target_index_ranges.find_first_of(",", current_pos);
			if ( next_comma != std::string::npos ) {
				range_string = target_index_ranges.substr(current_pos, next_comma - current_pos);
			}
			else {
				// Didn't find a comma: must be parsing the last range
				range_string = target_index_ranges.substr(current_pos);
				next_comma  = target_index_ranges.size();
			}

			// Stride: look for a colon, which separates the range from the stride
			pos_colon = range_string.find_first_of(":");
			if ( pos_colon != std::string::npos ) {
				token  = range_string.substr(pos_colon+1, range_string.size() - pos_colon);
				stride = atoi(token.c_str());
			}
			else {
				stride    = 1;
				pos_colon = next_comma;
			}

			// Range: look for a dash, which separates the start and stop indices
			pos_dash = range_string.find_first_of("-");
			if ( pos_dash != std::string::npos ) {
				token = range_string.substr(0, pos_dash);
				start = atoi( token.c_str() );
				
				token = range_string.substr(pos_dash+1, pos_colon);
				stop = atoi( token.c_str() );

				// Record indices
				int serial = start;
				while ( serial <= stop ) {
					target_indices.push_back(serial-1);
					serial += stride;
				}
			}
			else {
				// Single index
				int serial = atoi( range_string.c_str() );
				target_indices.push_back(serial-1);
			}

			// Move position to one past the current comma
			current_pos = next_comma + 1;
			if ( current_pos >= target_index_ranges.size() ) {
				done = true;
			}
		}
	}
	else if ( target_flags.isFile ) {
		//----- Target indices are in a file -----//

		// Format: whitespace-separated integers (can be over multiple lines)
		//  - Commented lines begin with "#"

		// Open file
		std::ifstream ifs(fileName);
		if ( ! ifs.is_open() ) {
			std::stringstream err_ss;
			err_ss << "Topology::getTargets - Unable to open file with target indices\n"
			       << "(input: " << fileName << ").\n";
			throw std::runtime_error( err_ss.str() );
		}

		int serial;
		std::string line;
		while ( std::getline(ifs, line) ) {
			std::stringstream ss(line);

			while ( ss >> token ) {
				// Skip commments
				if ( token[0] == '#' ) {
					break;
				}
				// Store index = serial - 1
				serial = std::atoi( token.c_str() );
				target_indices.push_back(serial - 1);
			}
		}
	}


	/*
	// DEBUG
	std::cout << "target_indices = \n";
	for ( int k=0; k<target_indices.size(); ++k ) {
		std::cout << "  " << target_indices[k] << "\n";
	}
	*/

	return;
}


// Compute the centers of mass of the molecules of the indicated indices
// - "positions_com" should point to molecule_indices.size() number of rvec arrays
void Topology::calculateMoleculesCentersOfMass(
	const std::vector<int>& molecule_indices, const RvecArray& coords,
	RvecArray& positions_com) const
{
	//----- Check input -----//

	const int num_atoms = coords.size();

	int numTargetMolecules = molecule_indices.size();
	bool isError = false;
	std::stringstream err_ss;
	if ( numTargetMolecules > static_cast<int>(molecules_.size()) ) {
		// Error: asking for too many molecules
		err_ss << "Error in Topology::calculateMoleculesCentersOfMass\n"
		       << "  The number of target molecules is greater than the number\n"
		       << "  of molecules in the topology.\n";
		isError = true;
	}
	else if ( num_atoms != num_atoms_ ) {
		err_ss << "Error in Topology::calculateMoleculesCentersOfMass\n"
		       << "  The number of atoms in the input does not match the number of\n"
		       << "  atoms in the topology.\n";
		isError = true;
	}
	if ( isError ) {
		/* 
		if ( coords     != nullptr ) {
			// About to crash, so it's fine to const_cast
			rvec* nonconst_coords = const_cast<rvec*>(coords);
			free(nonconst_coords);
			nonconst_coords = nullptr; 
		}
		*/
		//if ( positions_com != nullptr ) { free(positions_com); positions_com = nullptr; }
		throw std::runtime_error( err_ss.str() );

	}

	//----- Compute centers of mass -----//

	int molecule_index;
	std::string type_string = "", this_type_string = "";

	const MoleculeType* molecule_type_ptr = nullptr;

	for ( int j=0; j<numTargetMolecules; ++j ) {
		// Index in the molecules_ array
		molecule_index = molecule_indices[j];

		// Only update the molecule_type_ptr when a new type is needed.
		// - Since molecules of the same type are generally listed consecutively, this
		//   avoids excessively searching the std::map of strings to types
		// new type is found. 
		this_type_string = molecules_[molecule_index].type;
		if ( this_type_string != type_string ) {
			type_string = this_type_string;
			molecule_type_ptr = &( molecule_types.find(type_string)->second );

			// For whatever reason, the operator[] method for getting a const ptr to the
			// underlying object fails (because it uses "this->" syntax internally?)
			//molecule_type_ptr = &( molecule_types[type_string] );
		}

		// Weighted sum over atoms in the molecule
		const int num_atoms_in_molecule = molecule_type_ptr->num_atoms;
		const int global_start_index = molecules_[molecule_index].start_index;
		const std::vector<double>& masses = molecule_type_ptr->masses;
		this->calculateCenterOfMass(coords, global_start_index, num_atoms_in_molecule, masses,
		                            positions_com[j]);

		/*
		std::cout << "COM of molecule " << molecule_index+1 << "\n" // DEBUG
		          << "  type  "       << molecule_type_ptr->name << "\n"
		          << "  start_index " << molecules_[molecule_index].start_index+1 << "\n"
		          << "  #atoms      " << molecules_[molecule_index].num_atoms << "\n"
		          << "  m_tot       " << m_tot << "\n"
		          << "  COM(nm)     ";
		for ( int d=0; d<DIM_; ++d ) { std::cout << positions_com[j][d] << " "; }
		std::cout << "\n";
		*/
	}
}


// Compute the center of mass of a set of atoms in contiguous memory
void Topology::calculateCenterOfMass(
	const RvecArray& coords, const int global_start_index, const int num_atoms_in_set,
	const std::vector<double>& masses,
	Rvec com) const
{
	// Compute the center of mass in double precision to avoid roundoff issues
	// for large sets of atoms
	std::array<double, DIM_> com_d = {{ 0.0, 0.0, 0.0 }};

	double m_tot = 0.0, m_k;
	for ( int k=0; k<num_atoms_in_set; ++k ) {
		m_k = masses[k];
		for ( int d=0; d<DIM_; ++d ) {
			com_d[d] += m_k*coords[global_start_index+k][d];
		}
		m_tot += m_k;
	}

	// Normalize and copy to output vector
	for ( int d=0; d<DIM_; ++d ) {
		com_d[d] /= m_tot;
		com[d]   = com_d[d];
	}
}

// Compute the center of mass of a set of atoms determined by target indices
// - Masses are inferred from internal storage
void Topology::calculateCenterOfMass(
	const RvecArray& coords, const std::vector<int>& target_indices,
	Rvec com) const
{
	// Compute the center of mass in double precision to avoid roundoff issues
	// for large sets of atoms
	std::array<double, DIM_> com_d = {{ 0.0, 0.0, 0.0 }};

	int    i,   numTargets = target_indices.size();
	double m_i, m_tot = 0.0;
	for ( int k=0; k<numTargets; ++k ) {
		i   = target_indices[k]; // global index in coords array
		m_i = atom_masses_[i];
		for ( int d=0; d<DIM_; ++d ) {
			com_d[d] += m_i*coords[i][d];
		}
		m_tot += m_i;
	}

	// Normalize and copy to output vector
	for ( int d=0; d<DIM_; ++d ) {
		com_d[d] /= m_tot;
		com[d]   = com_d[d];
	}
}


//----- Utility Functions -----//

void Topology::clear()
{
	atom_types_.clear();
	molecule_types.clear();
	molecules_directive_entries_.clear();

	atoms_.clear();
	atom_names_.clear();
	atom_masses_.clear();
	molecules_.clear();
	residues_.clear();

	reference_positions_.clear();

	num_atoms_     = 0;
	num_molecules_ = 0;
}


void Topology::print() const
{
	// Print atom types
	std::cout << "[ atomtypes ]\n";
	for ( auto it=atom_types_.begin(); it!=atom_types_.end(); ++it ) {
		const AtomType& type = it->second;
		std::cout << type.name << "  " << type.atomic_number << "  "
		          << type.mass << "  " << type.charge << "  "
		          << type.ptype << "  " << type.sigma << "  " << type.epsilon << "\n";
	}

	// Print molecule types
	for ( auto it=molecule_types.begin(); it!=molecule_types.end(); ++it ) {
		const MoleculeType& type = it->second;

		int num_atoms = type.num_atoms;
		std::cout << "[ moleculetype ]\n"
		          << "  Name     " << type.name  << "\n"
		          << "  nrexl    " << type.nrexl << "\n"
		          << "  num_atoms " << num_atoms   << "\n"
		          << "\n";

		std::cout << "[ atoms ]\n";
		for ( int i=0; i<num_atoms; ++i ) {
			std::cout <<  i+1 << "  " << type.atom_types[i] << "  " 
			          << type.residue_numbers[i] << "  " << type.residue_types[i] << "  "
			          << type.atom_names[i] << "  " << type.charges[i] << "  " << type.masses[i] << "\n";
		}
		std::cout << "\n";
	}

	// Print list of molecules
	std::cout << "[ molecules ]  (unwrapped)\n";
	for ( int i=0; i<num_molecules_; ++i ) {
		std::cout << "  " << i+1 << ": type " << molecules_[i].type << " starting at index "
		          << molecules_[i].start_index+1 << " (" << molecules_[i].num_atoms << " atoms)\n";
	}
}
