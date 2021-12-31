# TrajTools

This is a small julia package used to read xyz files. It also contains functions to calculate mean-square displacement (MSD) and radial-distribution function (RDF) from the trajectory.

# Installation

Assuming you have this package in a folder /path/to/TrajTools/, you can add them to your julia environment using pkg:

```julia
]
add /path/to/TrajTools/
```

Afterwards simply include it like any other package via `Using TrajTools`. If you want to modify this code, you might want to use `dev` instead of `add`:

```julia
]
dev /path/to/TrajTools/
```

# Usage

In order to read an xyz-file /path/to/trajectory.xyz, you can run the following commands in the REPL or include them in your script:

```julia
Using TrajTools
coords, atoms = read_trajectory("/path/to/trajectory.xyz", com = True)
```

This will return a 3-dimensional array in coords, which contains all coordinates in the trajectory and a 1-dimensional array in atoms, which contains all labels given to the atoms in the first trajectory frame.
This will also remove the center of mass movement from the trajectory. If you don't want that, pass `False` as a second argument.

Once the trajectory has been read the most versatile function is `pbc_dist`. It calculates the distance between two points in real space.

...Further documentation pending...
