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
#include "xdrfile/xdrfile_xtc.h" /* read XTC files */
#pragma GCC diagnostic pop

#ifdef MPI_ENABLED
#include <mpi.h>
#endif /* MPI_ENABLED */

// Project headers
#include "orderparameters/CommonTypes.h"
#include "orderparameters/GenericFactory.h"
#include "orderparameters/GptlWrappers.h"
#include "orderparameters/MpiCommunicator.h"
#include "orderparameters/OpenMP.h"
#include "orderparameters/OrderParametersDriver.h"
#include "orderparameters/ProbeVolume.h"
#include "orderparameters/SimulationBox.h"
#include "orderparameters/Statistics.h"
#include "orderparameters/Topology.h"
#include "orderparameters/XdrFileTools.h"

#endif // MAIN_H
