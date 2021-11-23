using DelimitedFiles

include("read_trajec.jl")
include("calcs.jl")

path = "csh2po4_444_standard-pos-1.xyz"
pbc_path = "cdp_444"

pbc = readdlm(pbc_path)
@time coord, atom = read_trajectory(path, true)
if coord == nothing || atom == nothing
    println("Nothing is fine.")
else
    println("Everything is fine.")
end

#@time time_msd, msd = msd_direct(coord, 0.5, atom, "H", 100, 0.3)
#println(msd)
#println(time_msd)
#open("msd_julia.tsv", "w") do outfile
#    writedlm(outfile, [time_msd msd .* 1E4])
#end
