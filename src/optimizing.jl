using Profile
using ProfileView
using BenchmarkTools




# nn_lattice!(s)
#
# @code_warntype nn_lattice!(s)
#
@btime nn_lattice!(s)


#get_V_au_au(s)

#@code_warntype get_V_au_au(s)

#@btime get_V_au_au(s)

# get_E_au_n_neutral(s)
#
# @code_warntype get_E_au_n_neutral(s)
#
# @btime get_E_au_n_neutral(s)

# get_E_au_n_coup(s)
#
# @code_warntype get_E_ion(s)
#
# @btime get_E_all(s)
#
# hp = get_E_all(s)
#
# @profile get_E_all(s)
F_temp = zeros(Float64, N+2, 3)

get_F_all(s)

@benchmark get_F_au_au(s)

@code_warntype get_F_all(s)
#
@profiler  get_F_all(s)

@benchmark get_dm(s)

@code_warntype get_blk_akl!(s)

@benchmark get_blk_akl!(s)

get_dm!(s)

get_occupied_unoccpied_states!(s)

get_bkl!(s);

get_blk_akl!(s)
