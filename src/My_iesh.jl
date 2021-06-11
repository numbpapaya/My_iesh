module My_iesh


curr_vers = "v.0.8"
using LinearAlgebra
using Random
using StaticArrays
using Plots
using HDF5
using Statistics
using Dates

include("read_au_data.jl");
include("setup_parameters.jl");
include("types.jl");
include("misc_func.jl");
include("constructor.jl");
include("get_neighbors.jl");
include("get_energies.jl");
include("get_forces.jl");
include("propagation.jl");
include("run_trajectory.jl");


multiple_trajectory()

end
