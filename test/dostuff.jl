using DelimitedFiles
using TrajTools

path = "small.xyz"
pbc_path = "cdp_444"

@time traj = TrajTools.read_trajectory(path, pbc_path, true, false)
if traj == nothing
    println("Nothing is fine.")
else
    println("Everything is fine.")
end
