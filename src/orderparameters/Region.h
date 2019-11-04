#ifndef REGION_H
#define REGION_H

namespace Region {

// Non-overlapping regions of the simulation box that require different
// treatments
enum class RegionEnum {
	Vtilde,      // in the coarse-grained probe volume
	Shell_1,     // in shell 1 (outside the coarse-grained probe volume)
	Shell_2,     // in shell 2 (outside shell 1)
	Unimportant  // in a region that does not need to be considered
};

inline
static bool is_in_vtilde(const RegionEnum& region) {
	return ( region == RegionEnum::Vtilde );
}

inline
static bool is_in_shell_1(const RegionEnum& region) {
	return ( region == RegionEnum::Shell_1 );
}

inline
static bool is_in_shell_2(const RegionEnum& region) {
	return ( region == RegionEnum::Shell_2 );
}

inline
static bool is_in_shells(const RegionEnum& region) {
	return ( region == RegionEnum::Shell_1 or region == RegionEnum::Shell_2 );
}

}
#endif /* REGION_H */
