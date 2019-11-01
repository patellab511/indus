#include "GenericFactoryInitialization.h"

bool checkGenericFactoryInitialization()
{
	// Ensure that each registry is populated
	const auto& probe_volume_registry = ProbeVolumeRegistry::Factory::factory().get_registry();
	if ( probe_volume_registry.size() < 1 ) {
		throw std::runtime_error("Error: ProbeVolume registry is empty");
	}

	return true;
}
