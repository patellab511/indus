// Helper function to check for static initialization of GenericFactory objects

#pragma once
#ifndef GENERIC_FACTORY_INITIALIZATION_H
#define GENERIC_FACTORY_INITIALIZATION_H

#include <exception>
#include <stdexcept>

#include "GenericFactory.h"

// Base classes that use registries
#include "ProbeVolume.h"

bool checkGenericFactoryInitialization();

#endif //ifndef GENERIC_FACTORY_INITIALIZATION_H
