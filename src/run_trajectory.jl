export multiple_trajectories

@doc """
this function can be used to print in a threaded loop
see https://discourse.julialang.org/t/print-functions-in-a-threaded-loop/12112
"""
const print_lock = SpinLock()
const prints_pending = Vector{String}()
function tprintln(str)
	tid= Threads.threadid()
	str = "[Thread $tid]: " * string(str)
	lock(print_lock) do
		push!(prints_pending, str)
		if tid == 1 # Only first thread is allows to print
			println.(prints_pending)
			empty!(prints_pending)
		end
	end
end

@doc """
creates parent folder which is encoded with version and time info.
It then runs a parallel loop over the trajectories.
The results of each loop are logged and saved in a hdf5 file.
"""
function multiple_trajectory()
	current_time = Dates.now()
    str_time = string(Dates.format(current_time, "yyyy-dd-mm_HH-MM-SS"))
    #filepath_run = filepath_parent*curr_vers*"\\"*str_time *"\\"
    filepath_run = filepath_parent*curr_vers*"/"*str_time *"/"

    mkpath(filepath_run)
    Threads.@threads for traj in 1:numtraj
        tprintln(traj)

		#construct simulation
		s = simulation_init()
		simulation_constructor_x_v!(s)
		simulation_constructor_nn!(s)
		simulation_constructor_x_300K!(s)
		simulation_constructor_energy(s)
		simulation_constructor_force(s)
		propagate_init!(s)
		#run simulation
		simulate!(s)

		#logging
        hdf5filename = "traj_"*string(traj)*".h5"
        hdf5filepath = filepath_run*hdf5filename
		h5open(hdf5filepath, "w") do fid2
			write(fid2, "x_no", s.storage_xno)
			write(fid2, "v_no", s.storage_vno)
			write(fid2, "adsorbate orbital population", s.storage_aop)
			write(fid2, "orbital population", s.storage_op)
			write(fid2, "electronic excitation energies", s.storage_e)
			write(fid2, "inst temp", s.storage_temp)
			write(fid2, "max hopping probability", s.storage_phop)
			write(fid2, "kinetic energy", s.storage_K)
			write(fid2, "potential energy", s.storage_P)
			write(fid2, "hopping states", s.storage_state)
			write(fid2, "deltaK_hop NO", s.storage_deltaKNO)
			write(fid2, "deltaK_hop Au", s.storage_deltaKAu)
			write(fid2, "hoptimes", s.storage_hoptimes)
			write(fid2, "nf", s.nf)
			write(fid2, "attnum", s.attnum)
			write(fid2, "exnum", s.exnum)
			write(fid2, "temperature", tsurf)
			write(fid2, "e_trans_no",  e_trans_i)
			write(fid2, "e_vib_no", e_vib_i)
			write(fid2, "e_rot_no", e_rot_i)
			write(fid2, "incident angle", phi_inc)
			write(fid2, "trajectories", numtraj)
			write(fid2, "time stepsize", dt)
			write(fid2, "time steps", tsteps)
			write(fid2, "number electrons", Ne)
			write(fid2, "number of states", Ms)
			if extensive_logging == true
				write(fid2, "wavefunction", s.storage_psi)
				write(fid2, "adiabatic states", s.storage_phi)
				write(fid2, "x_au", s.storage_xau)
				write(fid2, "v_au", s.storage_vau)
			end
		end
	end
end
