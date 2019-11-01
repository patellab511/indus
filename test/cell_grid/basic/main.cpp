// Test CellGrid

#include "orderparameters/CellGrid.h"
#include "orderparameters/CommonTypes.h"
#include "orderparameters/MpiCommunicator.h"

int main(int argc, char* argv[])
{
	using namespace CommonTypes;

	// Enable SimulationState with debug mode on
	bool is_debug_mode = true;
	SimulationState simulation_state;
	simulation_state.set_debug_mode(is_debug_mode);

	MpiEnvironment::init();
	{
		// Wrap the communicator in a separate namespace so that it's cleaned up
		// before MPI_Finalize is called
		MpiCommunicator comm;

		CellGrid cell_grid(simulation_state, comm);

		// Test different grid sizes
		int max_num_cells = 5;
		std::array<int,DIM_> grid_dimensions;
		for ( grid_dimensions[X_DIM] = 1; grid_dimensions[X_DIM] < max_num_cells; ++grid_dimensions[X_DIM] ) {
			for ( grid_dimensions[Y_DIM] = 1; grid_dimensions[Y_DIM] < max_num_cells; ++grid_dimensions[Y_DIM] ) {
				for ( grid_dimensions[Z_DIM] = 1; grid_dimensions[Z_DIM] < max_num_cells; ++grid_dimensions[Z_DIM] ) {
					// This will invoke makeCells(), which determines the cell neighbor list
					// - Since debug mode is set, extra internal checks will be run
					cell_grid.set_grid_dimensions( grid_dimensions );
				}
			}
		}
	}
	MpiEnvironment::finalize();
}
