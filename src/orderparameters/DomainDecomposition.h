// TODO Make this into a more sophisticated class
// - update()
//   - Updates CellGrid_DomainDecomposition based on new state
//     - Make sure the bounding box in its local cell object gets properly updated!

#pragma once
#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "BoundingBox.h"
#include "CellGrid.h"
#include "CellGrid_DomainDecomposition.h"
#include "MpiCommunicator.h"
#include "SimulationBox.h"
#include "SimulationState.h"

class DomainDecomposition
{
 public:
	using Box   = SimulationBox::Box;
	using Real3 = SimulationBox::Real3;

	// TODO: Constructor that takes a raw MPI_Comm
	DomainDecomposition(
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator
	):
		simulation_state_(simulation_state),
  	mpi_communicator_( mpi_communicator ),
		cell_grid_( std::array<int,CellGrid::DIM_>({{1, 1, 1}}), 0.0, 0.0,
		            simulation_state_, mpi_communicator_ )
	{}

	MpiCommunicator& access_mpi_communicator() { 
		return mpi_communicator_; 
	}
	const MpiCommunicator& get_mpi_communicator() const { 
		return mpi_communicator_;
	}

	CellGrid_DomainDecomposition& access_cell_grid() { 
		return cell_grid_; 
	}
	const CellGrid_DomainDecomposition& get_cell_grid() const { 
		return cell_grid_;
	}

	const BoundingBox& getLocalCellBoundingBox() const {
		return cell_grid_.getLocalCellBoundingBox();
	}

 private:
	const SimulationState& simulation_state_;
	MpiCommunicator& mpi_communicator_;
	CellGrid_DomainDecomposition cell_grid_;
};

#endif // ifndef DOMAIN_DECOMPOSITION_H
