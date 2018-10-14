/* main.h
 *
 * ABOUT: Header with useful post-processing routines
 * DEVELOPMENT:
 *
 */

#pragma once
#ifndef MAIN_H
#define MAIN_H

// Standard headers
#include <cassert>  // for quick debugging
#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Library headers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment" /* ignore useless warnings from this header */
#include "xdrfile_xtc.h" /* read XTC files */
#pragma GCC diagnostic pop

#ifdef MPI_ENABLED
#include <mpi.h>
#endif /* MPI_ENABLED */

// Project headers
#include "CommonTypes.h"
#include "Indus.h"
#include "SimulationBox.h"
#include "Topology.h"
#include "XdrFileTools.h"

template<typename T, std::size_t dim>
double array_norm(const std::array<T,dim>& a) {
	double norm_sq = 0.0;
	for ( unsigned d=0; d<dim; ++d ) {
		norm_sq += a[d]*a[d];
	}
	return sqrt(norm_sq);
}

#endif // MAIN_H
