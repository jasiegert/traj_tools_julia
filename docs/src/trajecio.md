# Reading in xyz-files

Trajectories can be read from xyz-files and a text file containing the periodic boundary conditions.

```@docs
TrajTools.read_trajectory
```

This will return a `Trajectory` instance containing all the relevant information. More specifically it has the following elements:

    - coords: coordinates (array of shape (atom_no, frame_no) of StaticVectors, each representing a point in space)
    - atomlabels: labels given in the first frame (vector of length atom_no)
    - mdbox: periodic boundary conditions
    - unwrapped: whether or not the trajectory has been unwrapped
    - timestep_in_fs: time step used during the simulation in fs

For a small trajectory 'mini.xyz' containing two frames of a water molecule, this might look as follows:
```@example ; continued = true
    using TrajTools
    traj = read_trajectory("mini.xyz", "mini.pbc")
```
```@setup traj
    using TrajTools
    traj = read_trajectory("../mini.xyz", "../mini.pbc")
```

Coordinates and atom types can be accessed via the aforementioned elements 'coords' and `atomlabels`:
```@example traj
    traj.coords
```
```@example traj
    traj.atomlabels
```

Specific coordinates can be sliced from the coords as needed. For example to get the hydrogen positions in the second frame:
```@example traj
    traj.coords[traj.atomlabels .== "H", 2]
```
