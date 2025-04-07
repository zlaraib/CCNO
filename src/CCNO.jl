module CCNO

__precompile__()

# note, the order of these statements matter. Those without dependencies need to be first.

include("physics/Parameters.jl")
include("physics/SimulationState.jl")

include("initializations/initial_cond.jl")

include("utilities/store_data.jl")
include("utilities/save_plots.jl")
include("utilities/checkpoint_simulation_hdf5.jl")

include("physics/constants.jl")
include("physics/evolution.jl")
include("physics/gates_function.jl")
include("physics/heaviside.jl")
include("physics/momentum.jl")
include("physics/perturb.jl")
include("physics/shape_func.jl")


end # module
