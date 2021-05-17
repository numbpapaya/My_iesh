module My_iesh

using LinearAlgebra
using Distances
using Random
using StaticArrays
using Plots
include("read_au_data.jl")
include("setup_parameters.jl")
include("types.jl")
include("misc_func.jl")
include("constructor.jl")
include("get_neighbors.jl")
include("get_energies.jl")
include("get_forces.jl")
include("propagation.jl")


s = simulation_init();
simulation_constructor_x_v!(s);
simulation_constructor_nn!(s);
simulation_constructor_energy(s);
simulation_constructor_force(s);
propagate_init!(s);
simulate!(s);
end
