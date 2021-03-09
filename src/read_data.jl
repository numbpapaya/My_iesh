#path to data
filename = "C:\\Users\\Belal\\.julia\\dev\\My_iesh.jl\\data\\surface_Au111.dat"

f = open(filename)
data = readlines(f)
close(f)
const global N = Integer(parse(Float64, data[1]))
const global a = sqrt(2) * parse(Float64, data[2]) * Å
const global g = parse(Float64, data[2])
cell_temp = zeros(Float64, 3)
cell_temp[1] = parse(Float64, data[3])
cell_temp[2] = parse(Float64, data[4])
cell_temp[3] = sqrt(6) * parse(Float64, data[2])* Float64(100)
const global cell = copy(cell_temp)*Å

x_temp = zeros(Float64, (Integer(N), 3))

for i in 5:length(data)
    data_strip = strip(data[i])
    data_split = split(data_strip)
    str_to_float = parse.(Float64, data_split)
    x_temp[i-4, :] = str_to_float
end

const global x0 = copy(x_temp)
