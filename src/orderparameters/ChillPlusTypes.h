
#ifndef CHILL_PLUS_TYPES_H
#define CHILL_PLUS_TYPES_H

#include <array>
#include <string>
#include <vector>

//--------------------------------------------//
//----- CHILL+ Types and Data Structures -----//
//--------------------------------------------//

namespace ChillPlusTypes {

	//----- CHILL+ -----//

	// ChillPlusType: categories for classifying different solid-like waters
	// - Encompasses original CHILL+ categories, and CHILL* extension for
	//   heterogeneous interfaces
	enum class ChillPlusType : int {
		Liquid,
		HexagonalIce, 
		CubicIce,
		InterfacialIce,
		Clathrate, 
		InterfacialClathrate,
		SurfaceIce,           // CHILL*
		SurfaceClathrate,     // CHILL*
		COUNT   // Quick way to get the number of types
	};
	static const int NUM_CHILL_PLUS_TYPES = static_cast<int>(ChillPlusType::COUNT);

	static const std::array<std::array<std::string, 2>, NUM_CHILL_PLUS_TYPES>
		chill_plus_type_names = {{
			{{ "Liquid",                  "Liq"       }},
			{{ "Hexagonal ice",           "Hex"       }},
			{{ "Cubic ice",               "Cubic"     }},
			{{ "Interfacial ice",         "IntIce"    }},
			{{ "Clathrate",               "Clath"     }},
			{{ "Interfacial clathrate",   "IntClath"  }},
			{{ "Surface ice",             "SurfIce"   }}, // CHILL*
			{{ "Surface clathrate",       "SurfClath" }}  // CHILL*
		}};

	// Lumped CHILL+ categories
	enum class ChillPlusAggregateType : int {
		IceLike,
		ClathrateLike,
		SolidLike,
		IceLikeChillStar,       // CHILL*
		ClathrateLikeChillStar, // CHILL*
		SolidLikeChillStar,     // CHILL*
		Any,                    // Catch-all category (given name "Total" in output)
		COUNT
	};
	static const int NUM_CHILL_PLUS_AGGREGATE_TYPES 
			= static_cast<int>(ChillPlusAggregateType::COUNT);

	static const std::array<std::array<std::string, 2>, NUM_CHILL_PLUS_AGGREGATE_TYPES>
		chill_plus_aggregate_type_names = {{
			{{ "Ice-like",                "IceTot"    }},
			{{ "Clathrate-like",          "ClathTot"  }},
			{{ "Solid-like",              "SolidTot"  }},
			{{ "Ice-like (CHILL*)",       "IceTot*"   }},
			{{ "Clathrate-like (CHILL*)", "ClathTot*" }},
			{{ "Solid-like (CHILL*)",     "SolidTot*" }},
			{{ "Total (all categories)",  "Total"     }}
		}};

	// Ice-like types (standard CHILL+ only)
	const std::vector<ChillPlusType> ice_like_types = {
		ChillPlusType::HexagonalIce, 
		ChillPlusType::CubicIce, 
		ChillPlusType::InterfacialIce
	};

	// Clathrate-like types (standard CHILL+ only)
	const std::vector<ChillPlusType> clathrate_like_types = {
		ChillPlusType::Clathrate, 
		ChillPlusType::InterfacialClathrate 
	};

	const std::vector<ChillPlusType> chill_star_types = {
		ChillPlusType::SurfaceIce,
		ChillPlusType::SurfaceClathrate
	};

	// CHILL+ order parameter ranges for eclipsed (E) and staggered (S) bonds
	// - "Shifted" eclipsed bond range for TIP4P/Ice:  { -0.45  0.18 }
	const std::array<double, 2> staggered_bond_range = { {-1.0,  -0.8} };
	const std::array<double, 2> eclipsed_bond_range  = { {-0.35,  0.25} };

	// ChillPlusAtom: Organize CHILL+ information about a particular atom of interest
	// - Original CHILL+ settings
	//   - ChillPlusAtoms can have at most 4 nearest neighbors
	//   - If an atom has more than 4 neighbors according to the NeighborSphere,
	//     its CHILL+ neighbors are the 4 closest ones
	struct ChillPlusAtom {
		// Classification
		ChillPlusType chill_plus_type;  // CHILL+

		// Lumped category
		ChillPlusAggregateType chill_plus_aggregate_type;

		// Index ...
		int group_index;  // ... in target atom group
		int global_index; // ... globally

		// Lumped categories
		bool is_ice_like,            is_clathrate_like,            is_solid_like;            // CHILL+
		bool is_ice_like_chill_star, is_clathrate_like_chill_star, is_solid_like_chill_star; // CHILL*

		// Neighbors and pseudo-bonds
		std::vector<int>    nn_chill_plus; // Steinhardt atom indices of neighbors
		std::vector<double> cij_l_vec;     // c_l(i,j) for CHILL+ neighbors j
		int num_staggered_bonds;
		int num_eclipsed_bonds;
	};

	// Organize output from CHILL+ algorithm
	struct ChillPlusOutput {
		std::array<double, NUM_CHILL_PLUS_TYPES> chill_plus_counts;
		std::array<double, NUM_CHILL_PLUS_AGGREGATE_TYPES> chill_plus_aggregate_counts;
	};

} /* namespace ChillPlusTypes */

#endif /* CHILL_PLUS_TYPES_H */
