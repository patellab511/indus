/* Topology.h
 *
 * ABOUT: Object for reading and storing processed GROMACS topologies
 * NOTES:
 *  - Parse the top file before the gro file; the gro file will not 
 * DEVELOPMENT:
 *  - Probably not a good idea to use num_atoms_ to check whether things
 *    have been initialized...
 */

#pragma once

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <algorithm> // for std::transform
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project headers
#include "CommonTypes.h"
#include "StringTools.h"
#include "utils.h"

// Library headers
#ifdef INDUS_HAS_XDRFILE
#include "xdrfile_xtc.h"
#endif

class Topology
{
 public:
	//---------------------------------------//
	//----- Type and Struct Definitions -----//
	//---------------------------------------//

	static const int DIM_  = CommonTypes::DIM_;
	using Rvec      = CommonTypes::Rvec;
	using RvecArray = CommonTypes::RvecArray;
	using Matrix    = CommonTypes::Matrix;

	//using Real3 = rvec;

	struct AtomType {
		std::string name, ptype;
		int         atomic_number;
		double      mass, charge;
		double      sigma, epsilon;
	};

	// [ moleculetype] directive
	// - All vectors are num_atoms atoms long, and describe each atom
	//   - For example residue_types[i] is the type of residue to which atom i
	//     belongs (indexing *within* the molecule)
	struct MoleculeType {
		std::string              name;
		int                      num_atoms, nrexl;
		std::vector<std::string> atom_types, atom_names;
		std::vector<std::string> residue_types;
		std::vector<int>         residue_numbers;
		std::vector<double>      charges;
		std::vector<double>      masses;
	};

	// Entry under [ molecules ] directive in a GROMACS topology
	struct MoleculesEntry {
		std::string   molecule_name;
		int           count;
	};

	// Defines the type and location of a molecule in the global positions array
	struct Molecule {
		std::string type;
		int         start_index;
		int         num_atoms;
	};

	// Defines the type of a residue in the unrolled topology, and its location in
	// the global positions array
	struct Residue {
		std::string type;
		int start_index; // index of 1st atom of this residue in the global positions array
		int num_atoms;
		int number;      // from PDB file
		std::vector<std::string> atom_type_names; // for checking for heavy atoms
	};


	// Bundle of flags describing a set of target atoms or molecules
	struct TargetFlags {
		bool is_atom     = false;
		bool isMolecule = false;
		bool is_positions_list = false; // set of fixed positions

		bool isCenterOfMass = false;
		
		bool isType  = false;
		bool isName  = false;
		bool isIndex = false; // list of serial numbers (inline)
		bool isFile  = false; // file with serial numbers
	};


	//------------------------//
	//----- Constructors -----//
	//------------------------//

	Topology();

	// Type of input file is inferred from extension (e.g. ".gro" or ".top")
	Topology(const std::string& input_file);

	// (1) Read the GROMACS topology file "top_file" and set the simulation's topology.
	// (2) Read the .gro coordinate file "gro_file" and store the positions and box lengths
	//     as a reference state (reference_positions_, atom_names_)
	Topology(const std::string& top_file, const std::string& gro_file);


	//-----------------------------------//
	//----- Public Member Functions -----//
	//-----------------------------------//

	// Get indices of target atoms based on string 
	// - A "target atom" can be a virtual site, such as the center of mass of a molecule
	void getTargets(
		const std::vector<std::string>& target_tokens,
		// Output
		std::vector<int>& target_indices,
		TargetFlags&      target_flags
	) const;

	// Takes a string with serial numbers expressed in the form used by PLUMED,
	//    <first0>-<last0>:<stride0>,<first1>-<last1>:<stride1>,...
	// and appends the corresponding indices to 'target_indices'
	// - <stride> is optional and defaults to 1 
	void parseTargetIndices(
		const std::string& target_serial_ranges,
		std::vector<int>& target_indices
	) const;

	// Read a gro file with a *single frame* and return vectors with its contents
	// TODO velocities? more than 1 set of positions?
	// - Note that a .gro file has limited information about the system topology
	//   - HAS:
	//     - Atom *names* (NOT types)
	//     - Number of atoms
	//     - Residue names
	//   - DOESN'T HAVE:
	//     - Atom *types*
	//     - Any information about molecules: constituent atoms, connectivity, etc.
	//     - Force field parameters: masses, charges, bonded/non-bonded potential parameters, etc.
	void readGroFile(
		const std::string& gro_file,
		// Output
		std::vector<std::string>& residue_names, // resname corresponding to each atom
		std::vector<std::string>& atom_names,    // atom *types* would be in a .top file
		RvecArray& positions,
		Matrix&    box_matrix
	) const;

	//----- Compute the center(s) of mass of a set of atoms/molecules -----//

	// Compute the centers of mass of the molecules of the indicated indices
	// - "reference_positionscom" should point to molecule_indices.size() number of rvec arrays
	void calculateMoleculesCentersOfMass(
		const std::vector<int>& molecule_indices,
		const RvecArray& coords,
		// Output: centers of mass of indicates molecules
		RvecArray& reference_positionscom
	) const;

	// Compute the center of mass of a set of atoms in contiguous memory
	void calculateCenterOfMass(
		const RvecArray&           coords,              // All coords in topology
		const int                  global_start_index,
		const int                  num_atoms_in_set,
		const std::vector<double>& masses,
		Rvec                       com
	) const;

	// Compute the center of mass of a set of atoms determined by target indices
	// - Masses are inferred from internal storage
	void calculateCenterOfMass(
		const RvecArray&        coords,        // Set of all positions
		const std::vector<int>& target_indices, // Global atom indices
		Rvec                    com
	) const;

	// Clear the stored topology
	void clear();

	// Print topology to stdout
	void print() const;


	//----- Get Functions -----//

	int get_num_atoms() const { return num_atoms_; }

	const RvecArray& get_reference_positions()  const { return reference_positions_; }
	const Matrix&    get_reference_box_matrix() const { return reference_box_matrix_; }


	//----- Set Functions -----//

	// Same logic as constructor with the same input arguments (see above)
	// - Type of input file is inferred from extension (e.g. ".gro" or ".top")
	void setFromFile(const std::string& input_file);

	// Same logic as the constructor with the same input arguments (see above)
	// (1) Read the GROMACS topology file "top_file" and set the simulation's topology.
	// (2) Read the .gro coordinate file "gro_file" and store the positions and box lengths
	//     as a reference state (reference_positions_, atom_names_)
	void setFromFile(const std::string& top_file, const std::string& gro_file);

	// Read atom numbers (serials) from file, and return the list of indices
	// TODO allow for selecting atoms by GROMACS-style group name
	// - Format: whitespace-separated integers (can be over multiple lines)
	//   - Commented lines begin with "#" or ";"
	//   - Lines beginning with "[" are also ignored
	void readIndicesFromFile(
		const std::string& index_file,  // file with serial numbers
		std::vector<int>&  indices
	) const;

 private:
	// Input file
	std::string top_file_, gro_file_;

	// Type definitions
	std::map<std::string, AtomType>     atom_types_;
	std::map<std::string, MoleculeType> molecule_types;

	// Summary of the information under the [ molecules ] directive of a 
	// GROMACS topology file (# molecules of each type)
	// - Molecules are listed in the order they appear in the global positions array
	std::vector<MoleculesEntry> molecules_directive_entries_;

	// Global list of atom information (length num_atoms_): types, names, and masses
	std::vector<std::string> atoms_, atom_names_;
	std::vector<double>      atom_masses_;

	// Contains the sequence of molecules in the global positions array and
	// defines their types and the portion of the positions array they occupy
	// - Length num_molecules_
	std::vector<Molecule> molecules_;

	// Global list of residues, in the order they appear in the processed 
	// topology file (or gro file)
	std::vector<Residue> residues_;

	// Explicitly store the numbers of atoms, molecules, and residues
	// - This is important if the topology file is e.g. a gro file, and the atoms_
	//   array can't be used to infer the number of atoms since a gro file only has
	//   atom NAMES (not types)
	int num_atoms_, num_molecules_;

	StringTools string_tools_;

	// If the type of file use to generate the topology contains a sep of positions,
	// store it here.
	RvecArray reference_positions_;
	Matrix    reference_box_matrix_;

	// Read the given *processed* GROMACS topology and store its information
	// - Clears any stored topology
	void readGromacsTopology(const std::string& top_file);

	// Reads a .gro coordinate file and stores the positions and the box lengths 
	// it contains as a reference state
	// - If a Topology file has been provided, the Gro file is checked against the Topology.
	//   - If they differ, an exception is thrown
	void readReferenceStateFromGroFile(const std::string& gro_file);
};

#endif /* ifndef TOPOLOGY_H */
