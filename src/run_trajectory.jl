export multiple_trajectories
using Base.Threads

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


function one_trajectory()::Simulation
    s = simulation_init();
    simulation_constructor_x_v!(s);
    simulation_constructor_nn!(s);
    simulation_constructor_energy(s);
    simulation_constructor_force(s);
    propagate_init!(s);
    simulate!(s);
    return s
end

function multiple_trajectory()
    #filepath_parent = "W:\\ALL\\Theory Group\\iesh\\data\\"
    filepath_parent = "/home/braza2/data/"
    current_time = Dates.now()
    str_time = string(Dates.format(current_time, "yyyy-dd-mm_HH-MM-SS"))
    #filepath_run = filepath_parent*curr_vers*"\\"*str_time *"\\"
    filepath_run = filepath_parent*curr_vers*"/"*str_time *"/"

    mkpath(filepath_run)
    Threads.@threads for traj in 1:numtraj
        tprintln(traj)
        s = one_trajectory()
        hdf5filename = "traj_"*string(traj)*".h5"
        hdf5filepath = filepath_run*hdf5filename
        fid = h5open(hdf5filepath, "w")
        create_dataset(fid, "x_no", s.storage_xno)
        create_dataset(fid, "v", s.storage_vno)
        create_dataset(fid, "adsorbate orbital population", s.storage_aop)
        create_dataset(fid, "orbital population", s.storage_op)
        create_dataset(fid, "electronic excitation energies", s.storage_e)
        create_dataset(fid, "inst temp", s.storage_temp)
        create_dataset(fid, "max hopping probability", s.storage_phop)
        create_dataset(fid, "kinetic energy", s.storage_K)
        create_dataset(fid, "potential energy", s.storage_P)
        create_dataset(fid, "hopping states", s.storage_state)
        create_dataset(fid, "deltaK_hop NO", s.storage_deltaKNO)
        create_dataset(fid, "deltaK_hop Au", s.storage_deltaKAu)
        create_dataset(fid, "hoptimes", s.storage_hoptimes)
        create_dataset(fid, "nf", s.nf)
        create_dataset(fid, "attnum", s.attnum)
        create_dataset(fid, "exnum", s.exnum)
        create_dataset(fid, "temperature", tsurf)
        create_dataset(fid, "e_trans_no",  e_trans_i)
        create_dataset(fid, "e_vib_no", e_vib_i)
        create_dataset(fid, "e_rot_no", e_rot_i)
        create_dataset(fid, "incident angle", phi_inc)
        create_dataset(fid, "trajectories", numtraj)
        create_dataset(fid, "time stepsize", dt)
        create_dataset(fid, "time steps", tsteps)
        create_dataset(fid, "number electrons", Ne)
        create_dataset(fid, "Number of states", Ms)
        close(fid)
    end
end
