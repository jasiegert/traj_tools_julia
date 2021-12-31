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
