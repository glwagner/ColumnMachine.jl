# Standalone solver test with CLiMA driver
#
# Sets up a vertical temperature profile in multiple vertical (Z) columns of a horizontal (XY) 
# domain and then invokes a simple DG model to solve a 1-d vertical diffusion equation of the
# sort associated with ocean convective mixing parameterizations with a background diffusivity
# in stably stratified locations and an elevated diffusivity in unstabke regions.
#

## One time setup in clean directory, dedicated project tree and with local MPI install in /usr/local
## Assume ClimateMachine.jl (https://github.com/clima/climatemachine.jl) is in clone two levels up
```
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_MPI_PATH"/usr/local"
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=../../ClimateMachine.jl mysetup.jl
```

## Run
```
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_MPI_PATH="/usr/local"
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=../../ClimateMachine.jl driver.jl
```

## Using Package
```
Pkg.add(PackageSpec(url="https://github.com/clima/climatemachine.jl",rev="master"))
```

