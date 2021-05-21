using Profile
using ProfileView
using BenchmarkTools




# nn_lattice!(s)
#
# @code_warntype nn_lattice!(s)
#
@btime nn_lattice!(s)
rq = eigen(s.H)

#get_V_au_au(s)
mi = deepcopy(s.H)
aha = LAPACK.syev!('V', 'L', mi)
s.λ .= aha[1]
s.Γ .= aha[2]
#@code_warntype get_V_au_au(s)

#@btime get_V_au_au(s)

# get_e_au_n_neutral(s)
#
# @code_warntype get_e_au_n_neutral(s)
#
# @btime get_e_au_n_neutral(s)

# get_e_au_n_coup(s)
#
# @code_warntype get_e_ion(s)
#
# @btime get_e_all(s)
#
# hp = get_e_all(s)
#
# @profile get_e_all(s)
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



get_F_au_n_neutral(s)

get_F_au_o_neutral(s)

get_F_n_o_neutral(s)

get_F_au_n_ion(s)
get_F_au_o_ion(s)
get_F_n_o_ion(s)
F_au_n_c .= get_F_au_n_coup(s)
F_au_o_c .= get_F_au_o_coup(s)
F_coupling .= @. F_au_n_c + F_au_o_c


F_au .= get_F_au_au(s)
F_au_n_n .= get_F_au_n_neutral(s)
F_au_o_n .= get_F_au_o_neutral(s)
F_n_o_n_temp .= get_F_n_o_neutral(s)
F_n_o_n .= [F_n_o_n_temp; zeros(Float64, N, 3)]
F_neutral .= @. F_au_n_n + F_au_o_n + F_n_o_n - F_au




psi2 = collect([1.2358616421754511e-2,
  -1.3211748286157800e-3,
    0.99992134240773300 ,
           1.3936110728146633e-3,
              6.3955059627936444e-4,
                 3.3549186782848723e-4 ,
                   2.9694158126455239e-4,
      3.4448005730076524e-4,
     2.9591588060183586e-4,
    2.2393790827156967e-4,
    1.4074068340587945e-4])

vh1 = collect([
8.05440079e-03,
  0.999966562,
         1.22047227e-3,
            5.40201319e-4,
               3.10849166e-4 ,
                 1.74339104e-4,
                    1.57985167e-4,
                       1.89679558e-4,
                          1.68481318e-4,
                             1.30383385e-4,
                                8.28867342e-5
])
