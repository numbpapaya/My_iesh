export simulate!
export propagate_init!

@doc """
    propagate_init!(s)

Part of the initialization.
Initializes:    wavefunction s.ψ,
                eigenvalues s.λ,
                eigenvectors s.Γ,
                s.surfp,
                s.surfh,
                s.surfpinit,
                s.occnum,
                coupling matrices s.dhdv and s.dhdea
                storage_arrays...
"""
function propagate_init!(s::Simulation) #check seems okay #everything good.
    compute_eigen_H!(s)
    H_derivatives!(s)
    generate_therm_surf!(s)

    #set initial wavefunction
    s.ψ .= convert(Array{ComplexF64, 2}, view(s.Γ, :, s.surfp))

    #calculate initial forces
    update_forces!(s)
    #fix position of the last layer of atoms
    s.F[399:530, :] .= 0.0
    s.v[399:530, :] .= 0.0

    #convert to Forces
    s.F .= -s.F

    storage_init!(s)

end


#stores data from initialization into first columns of storage arrays
@inline function storage_init!(s::Simulation)
    get_temperature!(s, 1)
    s.storage_K[1] = 0.5 * sum(s.v.^2 .* m_spread)
    s.storage_P[1] = s.hp[1] + sum(s.λ[s.surfp])
    rtemp = norm(s.Δ_no)
    s.storage_e[1, 1] = 0.5 * (mass_arr[1] + mass_arr[2]) *
     sum(((mass_arr[1].*s.v[1, :] + mass_arr[2].*s.v[2, :])/(mass_arr[1] + mass_arr[2])).^2)
    @views T_vib = 0.50 * μ * dot((s.v[1, :] .- s.v[2, :]), s.Δ_no/rtemp)
    U_vib = F_n *(1.0 - exp(-δ_n * (rtemp -  r_0_NO)))^2
    E_vib = U_vib + T_vib
    s.storage_e[2, 1] = E_vib
    @views T_tot = 0.5 * mass_arr[1] * sum(s.v[1, :].^2) + 0.5*mass_arr[2]*sum(s.v[2, :].^2)
    s.storage_e[3, 1] = T_tot - T_vib - s.storage_e[1, 1]
    s.storage_e[4, 1] = sum(s.λ[s.surfp]) - sum(s.λ[s.surfpinit])       #das hier macht keinen Sinn
    if extensive_logging == 1
        s.storage_xno[1:3, 1] = view(s.x, 1, :)
        s.storage_xno[4:6, 1] = view(s.x, 2, :)
        s.storage_vno[1:3, 1] = view(s.v, 1, :)
        s.storage_vno[4:6, 1] = view(s.v, 2, :)#check seems okay
        s.storage_psi[:, 1] = vec(s.ψ)
        s.storage_xau[:, 1] = vec(s.x[3:398, :])
        s.storage_vau[:, 1] = vec(s.v[3:398, :])
    end
end

@inline function compute_eigen_H!(s::Simulation)
    s.H .= Symmetric(h0 .+ vm .* s.hp[3] ./ sqrt_de)
    s.H[1, 1] = s.H[1, 1] + s.hp[2] - s.hp[1]
    s.λ .= eigvals(s.H)
    s.Γ .= eigvecs(s.H)#check seems okay
end


@inbounds @inline function H_derivatives!(s::Simulation) #seems ok, fixed error
    for i in 1:Ms+1
        s.dhdea[:, i] .= view(s.Γ, 1, :) .* s.Γ[1, i]
    end
    mul!(s.temp_vm_gamma, vm, s.Γ)
    mul!(s.dhdv, transpose(s.Γ), s.temp_vm_gamma)
    s.dhdv .= s.dhdv ./ sqrt_de
end


@inline @inbounds function generate_therm_surf!(s::Simulation) #check seems okay
    #generate thermalized populations for electronic states
    for j in 1:Ne
        s.surfp[j] = j
    end

    for j in 1:Ne
        s.occnum[s.surfp[j]] = 1
    end

    k = 1
    for j in 1:Ms+1
        if s.occnum[j] == 0
            s.surfh[k] = j
            k = k + 1
        end
    end

    for t in 1:Ne
        for k in 1:Ne*(Ms + 1 - Ne)
            #hoprand = rand(Float64, 3)
            hoprand = collect([0.2, 0.4, 0.01])
            jp = Int(floor(hoprand[1]*Ne)) + 1
            jh = Int(floor(hoprand[2] * (Ms + 1 - Ne))) + 1
            rtemp = exp(-(s.λ[s.surfh[jh]] - s.dhp_neutral[s.surfp[jp]])) / (kb * tsurf)
            if hoprand[3] < rtemp
                m = s.surfh[jh]
                s.surfh[jh] = s.surfp[jp]
                s.surfp[jp] = m
            end
        end
    end

    s.surfpinit .= deepcopy(s.surfp)
end



@inline @inbounds function update_forces!(s::Simulation) #check seems okay
    for j in 1:Ne
        s.F .= s.F .+ (s.dhdea[s.surfp[j], s.surfp[j]] .* (s.dhp_ion .- s.dhp_neutral) .+
        s.dhdv[s.surfp[j], s.surfp[j]] .* s.dhp_coup)
    end
    s.F .= s.F .+ s.dhp_neutral
end

#---------------------simulation main time loop for n=2:tsteps-----------------
function simulate!(s::Simulation)
    for n in 2:tsteps
        do_we_stop = v_z_condition_check(s)
        if do_we_stop == -1
            break
        end
        # Record minimial distance of z approach and orientational angle
        @views min_z_no = minimum(s.x[1:2, 3])
        if min_z_no <= s.trajzmin[1]
            s.trajzmin .= min_z_no
            s.trajtheta .= acos(s.Δ_no[3]/norm(s.Δ_no))
        end
        #t_prog = round(n/tsteps * 100.0; digits=3)
        #str_t_prog = string(t_prog)*"%"
        #traj_prog = round(min_z_no/Å; digits=3)
        #str_trajzmin_prog = string(traj_prog)
        #println(str_t_prog*", "*str_trajzmin_prog)


        # get state corresponding to current surface
        @inbounds for j in 1:Ne
            s.ϕ[:, j] = view(s.Γ, :, s.surfp[j])
        end

        #Calculate nonadiabatic coupling matrix DM between adiabatic orbitals
        get_dm!(s::Simulation)
        # Get occupied (particle), unoccupied (hole) states corresponding to current surface specified by surfp
        # surfp is Ne x 1 array of particle states
        # surfh is (Ms+1-Ne) x 1 array of hole states
        get_occupied_unoccpied_states!(s)
        # Get b_{kl} where k is current many-electron adiabatic state and l is possible new
        # 	many-electron adiabatic state
        # k corresponds to occupation of states specified by surfp
        # l corresponds to occupation of states specified by a one electron-hole-pair excitation out of surfp
        get_blk2_pbmaxest!(s)
        # Calculate upper bound of maximum hopping probability, Pbmaxest

        #beginn hopping
        #Generate random number for surface hopping
        hoprand = rand()
        #! Only attempt hop if hoprand < Pbmaxest to avoid expensive calculation
        #! of real hopping elements blk
        if hoprand < s.Pbmaxest[1]
            hopping!(s, hoprand, n)
        end
        # Propagate electronic Hamiltonian

        propagate_hamiltonian!(s)
        #storage
        storage_sim1!(s, n)

        #propagate nuclear motion
        propagate_nuclear!(s)
        #get adiabatic eigenvalues, eigenvectors
        get_eigen_sim!(s)
        # Get dH / d Ep matrices where Ep is parameter in Newns-Anderson Hamiltonian
        # (i.e. Ep = V or Ep = E_a = (EI-EN) )
        get_dhdea_dhdv_loop!(s)

        # Calculate forces
        s.F .= 0
        update_forces!(s)
        s.F[399:530, :] .= 0
        s.v[399:530, :] .= 0
        s.F .= -s.F

        #propagate velocities
        s.v .= s.vtemp .+ 0.5 ./m_spread .* s.F * dt


        s.storage_xno[1:3, n] = view(s.x, 1, :)
        s.storage_xno[4:6, n] = view(s.x, 2, :)
        s.storage_vno[1:3, n] =view(s.x, 1, :)
        s.storage_vno[4:6, n] = view(s.v, 2, :)#check seems okay
        if extensive_logging == 1
            s.storage_psi[:, n] = vec(s.ψ)
            s.storage_xau[:, n] = vec(s.x[3:398, :])
            s.storage_vau[:, n] = vec(s.v[3:398, :])
            s.storage_phi[:, n] = vec(s.ϕ)
        end
        
        s.nf .= s.nf .+ one(Int64)
    end
    s.storage_e[5, :] = sum(s.storage_e[1:4, :], dims=1) .- sum(s.storage_e[1:4, 1])

end
#-------------------------------------------------------------------------------

#check if the projectile has reached simulation boundaries
function v_z_condition_check(s::Simulation) #should be ok
    @views v_cond = (s.v[1, 3]*m_N + s.v[2, 3]*m_O)/(m_N + m_O)
    @views z_cond = (s.x[1, 3] + s.x[2, 3])/2.0
    if z_cond > z_end && v_cond > 0.0
        #println("Boundary reached:End Simulation")
        return -1
    end
    return 1
end


@inbounds @inline function get_dm!(s::Simulation) #should be ok
    # Calculate nonadiabatic coupling elements DM between adiabatic orbitals s.dm
    # phipsi = <phi_k|psi> is overlap between "current" adiabatic state |phi_k> and electronic state |psi>
    mul!(s.temp_phi_psi,transpose(s.ϕ) , s.ψ)
    s.phipsi[1] = det(s.temp_phi_psi)
    # Calculate akk as defined in Tully JCP 1990
    s.akk[1] = abs(s.phipsi[1])^2 #corrected error

    # Calculate nonadiabatic coupling elements DM between adiabatic orbitals

    temp_2 = s.dhp_ion .- s.dhp_neutral #preallocation
    temp_2 .= temp_2 .* s.v
    temp_2_sum = sum(temp_2)
    s.dm .= temp_2_sum .* s.dhdea

    temp_3 = s.dhp_coup .* s.v
    temp_3_sum = sum(temp_3)
    s.dm .= s.dm .+ temp_3_sum .* s.dhdv

    temp_rep1 = repeat(s.λ, 1, Ms+1)
    temp_rep2 = temp_rep1' .- temp_rep1
    s.dm .= s.dm ./ temp_rep2
    for j in 1:Ms+1
        s.dm[j,j] = zero(eltype(s.dm))
    end #should be ok
end


@inline @inbounds function get_occupied_unoccpied_states!(s::Simulation) #should be ok
    s.occnum .= 0
    for j in 1:Ne
        s.occnum[s.surfp[j]] = 1
    end
    k = 1
    for j in 1:Ms+1
        if s.occnum[j] == 0
            s.surfh[k] = j
            k = k+1
        end
    end #should be ok
end

@inline @inbounds function get_blk2_pbmaxest!(s::Simulation) #should be ok
    rtemp = abs(s.phipsi[1] * s.phipsi[1])
    if geaq(rtemp, 1.0)
        rtemp = 1.0
    end
    for jp in 1:Ne
        for jh in 1:Ms+1-Ne
            s.blk2[jp, jh] = 2.0 * abs(s.phipsi[1] * s.dm[s.surfp[jp], s.surfh[jh]])
        end
    end

    s.Pbmaxest[1] = sqrt(sum(s.blk2.^2)) * sqrt(1.0 - rtemp)/s.akk[1] * dt * Float64(thop) #equation not found in any papers, but cool approach.
end


#---------------------start hopping! hoprand < pbmaxest -----------------------
function hopping!(s::Simulation, hoprand::Float64, n::Int64)
    hopping!_pbmaxest!(s, hoprand, n)
end

function hopping!_pbmaxest!(s::Simulation, hoprand::Float64, n::Int64)
    get_blk_akl!(s)
    get_pb!(s)

    s.storage_phop[n] = s.Pb[Ne*(Ms+1-Ne)]
    #attempt hop conditioned on random number
    if hoprand < s.Pb[Ne*(Ms+1-Ne)]
        #print("Hopping! ")
        hopping!_pbmaxest!_attempt!(s, hoprand, n)
    end
end

@inline function get_blk_akl!(s::Simulation) #seems ok, but numerically unstable, cancellation of random phase needs to be proven
    for jp in 1:Ne
        for jh in 1:Ms+1-Ne
            s.surfpnew .= s.surfp
            s.surfpnew[jp] = s.surfh[jh]
            for j in 1:Ne
                s.ϕ[:, j] = view(s.Γ, :, s.surfpnew[j])
            end
            #Ctemp = <phi_l|psi> is overlap between "new" adiabatic state
            #|phi_l> and electronic state |psi>
            mul!(s.temp_akl, transpose(s.ϕ) , s.ψ)
            ctemp = det(s.temp_akl)
            s.akl[1] = s.phipsi[1] * conj(ctemp)
            s.blk[jp, jh] = 2.0 * real(s.akl[1] * s.dm[s.surfp[jp], s.surfh[jh]])
        end
    end
end

@inline function get_pb!(s::Simulation) #seems ok, but probabilities can be greater than one.
    rtemp = 0.0
    for jh in 1:Ms+1-Ne
        for jp in 1:Ne
            s.Pb[(jh-1) * Ne + jp] = rtemp
            s.Pb[(jh-1) * Ne + jp] = s.Pb[(jh-1)*Ne + jp] + s.blk[jp, jh]*(s.blk[jp, jh] > 0.0)
            rtemp = s.Pb[(jh - 1) * Ne + jp]
        end
    end
    s.Pb .= s.Pb ./ s.akk[1] * dt
end
#-----------------start hopping! hoprand < pbmaxest && hoprand < Pbmax ---------
function hopping!_pbmaxest!_attempt!(s::Simulation, hoprand::Float64, n::Int64)
    hopstate = which_state(s, hoprand) #ok
    jh= which_orbital_jh(hopstate) #ok
    jp= which_orbital_jp(hopstate) #ok
    rhop = which_direction(s, jp, jh) #ok
    e_el_old, e_el_new, k_rhop, k_temp = check_energy1(s, jp, jh, rhop)
    bool_e = check_energy2(e_el_old, e_el_new, k_rhop, k_temp)
    s.attnum .= s.attnum .+ one(Int64)

    if bool_e == true
        update_vel_populations!(s, jp, jh, rhop, k_temp, k_rhop, e_el_old, e_el_new, n) #needs more storage
    end
end

@inline function which_state(s::Simulation, hoprand::Float64)::Int64 #ok
    hopstate = one(Int64)
    while s.Pb[hopstate] <= hoprand
        hopstate = hopstate + one(Int64)
    end
    return hopstate
end

@inline function which_orbital_jh(hopstate::Int64)::Int64 #ok
    jh = Int(floor((hopstate - 1) / Ne)) + one(Int64)
    return jh
end
@inline function which_orbital_jp(hopstate::Int64)::Int64 #ok
    jp = ((hopstate - one(Int64)) % Ne) + one(Int64)
    return jp
end

@inline function which_direction(s::Simulation, jp::Int64, jh::Int64)::float_array #not sure what happenes here, literature unclear
    @views temp = s.dhdea[s.surfp[jp], s.surfh[jh]] .* (s.dhp_ion .- s.dhp_neutral)
    @views temp .= temp .+ s.dhdv[s.surfp[jp], s.surfh[jh]] .* s.dhp_coup
    temp2 = sqrt(sum(temp.^2)/sum(s.v .^2))
    temp .= temp ./ temp2
    temp3 = sum(temp .* s.v)
    temp .= temp .* copysign(1.0, temp3)
    return temp
end


@inline function check_energy1(s::Simulation, jp::Int64, jh::Int64, rhop::float_array) # seems ok
    s.surfpnew .= s.surfp
    s.surfpnew[jp] = s.surfh[jh]
    e_el_old = sum(s.λ[s.surfp])
    e_el_new = sum(s.λ[s.surfpnew])
    k_rhop = 0.5 * sum(m_spread .* rhop .^2)
    k_temp = sum(m_spread .* s.v .* rhop)

    return e_el_old, e_el_new, k_rhop, k_temp
end


@inline function check_energy2(e_el_old::Float64, e_el_new::Float64, k_rhop::Float64, k_temp::Float64)::Bool
    bool_e = k_temp^2 - 4.0*k_rhop * (e_el_new - e_el_old) > 0.0
    return bool_e
end

@inline function update_vel_populations!(s::Simulation, jp::Int64, jh::Int64, rhop::float_array,
    k_temp::Float64, k_rhop::Float64, e_el_old::Float64, e_el_new::Float64, n::Int64)

    s.vscale[1] = (k_temp - sqrt(k_temp^2 - 4.0 * k_rhop * (e_el_new - e_el_old)))/2.0/k_rhop
    k_no = 0.50 * sum(view(m_spread, 1:2, :).* view(s.v, 1:2, :).^2)
    k_au = 0.50 * sum(view(m_spread, 3:N+2, :) .* view(s.v, 3:N+2, :).^2)
    s.v .= s.v .- s.vscale[1] .* rhop
    s.storage_deltaKNO[s.exnum[1]+1] = k_no - 0.5*sum(view(m_spread, 1:2, :) .* view(s.v, 1:2, :) .^2)
    s.storage_deltaKAu[s.exnum[1]+1] = k_au - 0.5*sum(view(m_spread, 3:N+2, :) .* view(s.v, 3:N+2, :) .^2)
    s.storage_state[2*s.exnum[1] + 1] = s.surfp[jp]
    s.storage_state[2*s.exnum[1] + 2] = s.surfh[jh]
    s.surfp .= s.surfpnew
    s.storage_hoptimes[s.exnum[1]+1] = n
    s.exnum .= s.exnum .+ one(Int64)
    #endif
end
#-----------------END hopping! hoprand < pbmaxest && hoprand < Pbmax ---------



@inline function propagate_hamiltonian!(s::Simulation) # seems ok
    s.ψ .= transpose(s.Γ) * s.ψ
    s.uu .= exp.(-1im * s.λ * dt/hbar)
    s.uuu .= repeat(s.uu, 1, Ne)
    s.ψ .= s.uuu .* s.ψ
    s.ψ .= s.Γ * s.ψ
end

#storage_op_e
@inline function storage_sim1!(s::Simulation, n::Int64) #seems ok
    get_temperature!(s, n)
    s.storage_K[n] = 0.5 * sum(s.v .^2 .* m_spread)
    s.storage_P[n] = s.hp[1] + sum(s.λ[s.surfp])
    @views s.storage_aop[n] = sum(abs.(s.ψ[1, :]).^2) #rhoa
    store_temp = sum(abs.(transpose(s.Γ) * s.ψ).^2, dims=2)
    s.storage_op[:, n] .= store_temp[:, 1] #psip

    rtemp = norm(s.Δ_no)
    s.storage_e[1, n] = 0.5 * (mass_arr[1] + mass_arr[2]) *
     sum(((mass_arr[1].*s.v[1, :] + mass_arr[2].*s.v[2, :])/(mass_arr[1] + mass_arr[2])).^2)
    @views T_vib = 0.50 * μ * dot((s.v[1, :] .- s.v[2, :]), s.Δ_no/rtemp)
    U_vib = F_n *(1.0 - exp(-δ_n * (rtemp -  r_0_NO)))^2
    E_vib = U_vib + T_vib
    s.storage_e[2, n] = E_vib
    @views T_tot = 0.5 * mass_arr[1] * sum(s.v[1, :].^2) + 0.5*mass_arr[2]*sum(s.v[2, :].^2)
    s.storage_e[3, n] = T_tot - T_vib - s.storage_e[1, n]
    s.storage_e[4, n] = sum(s.λ[s.surfp]) - sum(s.λ[s.surfpinit])

end
@inline function propagate_nuclear!(s::Simulation)
    s.vdot .= inv_m_spread .* s.F
    s.x .= s.x .+ s.v * dt .+ 0.50*s.vdot * dt^2
    s.Δ_no .=  s.x[1, :] - s.x[2, :]
    s.vtemp .= s.v .+ 0.50*s.vdot*dt
end

@inline function get_eigen_sim!(s::Simulation)
    s.hp .= get_E_all(s)    #update energy surface
    get_F_all(s)            #update partial derivatives dhp_i
    s.H .= h0 .+ vm .* s.hp[3] ./sqrt_de
    s.H[1, 1] = s.H[1, 1] + s.hp[2] - s.hp[1]
    # s.Γ_hold = deepcopy(s.Γ) #ARRAY ALLOCATION !!!! solve later
    s.λ .= eigvals(s.H)
    s.Γ .= eigvecs(s.H)
    # @inbounds for j in 1:Ms+1
    #     s.temp_prop_gamma .= abs.(transpose(Γ_hold) * view(s.Γ, :, j))
    #     Γmaxloc = argmax(temp)
    #     temp2 = dot(view(Γ_hold,:, Γmaxloc), view(s.Γ, :, j))
    #     s.Γ[:, j] = s.Γ[:, j]/copysign(1.0, temp2)
    # end
end

function get_dhdea_dhdv_loop!(s::Simulation)
    s.blk .= 0.0
    @inbounds for j in 1:Ms+1
        s.dhdea[:, j] = view(s.Γ, 1, :) .* s.Γ[1, j]
    end
    mul!(s.temp_vm_gamma, vm , s.Γ)
    mul!(s.dhdv, transpose(s.Γ), s.temp_vm_gamma)
    s.dhdv .= s.dhdv ./sqrt_de
end

function get_temperature!(s::Simulation, n::Int)
    speed = [norm(s.v[i, :]) for i in 3:398]
    s.storage_temp[n] = 2*m_au*sum(speed.^2)/(3*(N*0.75-3)*kb) #is this correct?
end
