#path to data

#read velocity file
#filename = "C:\\Users\\Belal\\.julia\\dev\\My_iesh.jl\\data\\surface_Au111.dat"
#filename = "/home/razamaza/dev/My_iesh.jl/data/surface_Au111.dat"
#filename = "C:\\Users\\braza2\\Documents\\GitHub\\My_iesh.jl\\data\\surface_Au111.dat"
#filename = "/home/braza2/My_iesh.jl/data/surface_Au111.dat"
filename = datapath*"surface_Au111.dat"

f = open(filename)
data = readlines(f)
close(f)

const global N = Int64(parse(Float64, data[1])) #number of grid atoms
const global g = parse(Float64, data[2])
const global a = sqrt(2) * g


cell_temp = zeros(Float64, 3)
cell_temp[1] = parse(Float64, data[3])
cell_temp[2] = parse(Float64, data[4])
cell_temp[3] = sqrt(6) * parse(Float64, data[2])* Float64(100)
const global cell = SVector{3, Float64}(copy(cell_temp)) #super cell


#get equilibrium position vectors of surface (T=0K)
x_au_temp = zeros(Float64, (Integer(N), 3)) #create temporary array to hold positions of grid atoms
for i in 5:length(data)
    data_strip = strip(data[i])
    data_split = split(data_strip)
    str_to_float = parse.(Float64, data_split)
    x_au_temp[i-4, :] = str_to_float
end

const global x_au0 = copy(x_au_temp) #position vectors of grid atoms


#get velocities for thermalized surface
#filename ="C:\\Users\\braza2\\Documents\\GitHub\\My_iesh.jl\\data\\Latticev300.dat"
#filename = "/home/braza2/My_iesh.jl/data/Latticev300.dat"
filename = datapath*"Latticev300.dat"

f = open(filename)
data = readlines(f)
close(f)


v_au_temp = zeros(Float64, (Integer(N), 3))

for i in 1:length(data)
    data_strip = strip(data[i])
    data_split = split(data_strip)
    str_to_float = parse.(Float64, data_split)
    v_au_temp[i, :] = str_to_float
end
#x_no_temp = vcat(x_no, x_temp)
const global v_au0 = copy(v_au_temp)


#get position vectors of thermalized surface
#filename ="C:\\Users\\braza2\\Documents\\GitHub\\My_iesh.jl\\data\\Latticex300.dat"
#filename = "/home/braza2/My_iesh.jl/data/Latticex300.dat"
filename = datapath*"Latticex300.dat"

f = open(filename)
data = readlines(f)
close(f)

x_300K_temp = zeros(Float64, N, 3)
for i in 1:length(data)
    data_strip = strip(data[i])
    data_split = split(data_strip)
    str_to_float = parse.(Float64, data_split)
    x_300K_temp[i, :] = str_to_float
end

global const x_300K0 = copy(x_300K_temp)
