module My_iesh


using LinearAlgebra
using Distances
using NearestNeighbors
using PyCall


#define parameters
include("constants.jl")
include("fit_parameters.jl")

#define functions
include("neighbors.jl") #needs testing!!!!
#include("matrix_elements.jl")


#begin reading data
include("read_data.jl")


end
