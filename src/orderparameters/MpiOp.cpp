#include "MpiOp.h"

const std::unordered_map<MpiOp::StandardOp, std::string, MpiOp::EnumClassHash> 
	MpiOp::standard_mpi_op_names_ =  {
	{ StandardOp::Null,    "Null"    },
	{ StandardOp::Max,     "Max"     },
	{ StandardOp::Min,     "Min"     },
	{ StandardOp::Sum,     "Sum"     },
	{ StandardOp::Product, "Product" },
	{ StandardOp::Land,    "Land"    },
	{ StandardOp::Band,    "Band"    },
	{ StandardOp::Lor,     "Lor"     },
	{ StandardOp::Bor,     "Bor"     },
	{ StandardOp::Lxor,    "Lxor"    },
	{ StandardOp::Bxor,    "Bxor"    },
	{ StandardOp::Minloc,  "Minloc"  },
	{ StandardOp::Maxloc,  "Maxloc"  },
	{ StandardOp::Replace, "Replace" }
};


const std::string& MpiOp::get_name(const MpiOp::StandardOp& op) {
	const auto it = standard_mpi_op_names_.find(op);
	if ( it != standard_mpi_op_names_.end() ) {
		return it->second;
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  A standard MPI operation does not have a registered name."
		          << " This should never happen.\n";
		throw std::runtime_error( err_ss.str() );
	}
}
