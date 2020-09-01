using MPI
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes

using Dates
using ClimateMachine.Diagnostics
using CLIMAParameters

macro oldnew()
 return eval(newstyle)
end

import  ClimateMachine.SystemSolvers:
    BatchedGeneralizedMinimalResidual, linearsolve!

# Needed to control diagnostics dispatch
using ClimateMachine.ConfigTypes
struct IVDCConfigType <: ClimateMachineConfigType end

include("MYDIAGSModel.jl")
include("IVDCModel.jl")
include("ivdc_diagnostics.jl")

 #####
 # Basic initialization
 #####
 ClimateMachine.init()
 ArrayType = ClimateMachine.array_type()
 mpicomm = MPI.COMM_WORLD
 FT = Float64

 #####
 # Create mesh
 #####
 const N = 4
 # const Nˣ = 20
 # const Nʸ = 20
 const Nˣ = 2
 const Nʸ = 2
 const Nᶻ = 20
 const Lˣ = 4e6  # m
 const Lʸ = 4e6  # m
 const H = 1000  # m

 xrange = range(FT(0); length = Nˣ + 1, stop = Lˣ)
 yrange = range(FT(0); length = Nʸ + 1, stop = Lʸ)
 zrange = range(FT(-H); length = Nᶻ + 1, stop = 0)

 gconf  = (N=N,Nˣ=Nˣ,Nʸ=Nʸ,Lˣ=Lˣ,Lʸ=Lʸ,H=H,xrange=xrange,yrange=yrange,zrange=zrange)

 brickrange_2D = (xrange, yrange)
 topl_2D       =
   BrickTopology( mpicomm, brickrange_2D, periodicity = (false, false) );
 grid_2D = DiscontinuousSpectralElementGrid(
   topl_2D,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
 )
 brickrange_3D = (xrange, yrange, zrange)
 topl_3D = StackedBrickTopology(
   mpicomm,
   brickrange_3D;
   periodicity = (false, false, false),
   boundary = ((1, 1), (1, 1), (2, 3)),
 )
 grid_3D = DiscontinuousSpectralElementGrid(
   topl_3D,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
 )

 # Wave speeds for numerical flux
 const cʰ = 1  # typical of ocean internal-wave speed
 const cᶻ = 0

 # Set a timestep for implicit solve
 dt = 5400/5

 # Create a container balance law for working with multiple DG models in Diagnostics

 # Create balance law and RHS arrays for diffusion equation
 ivdc_dg = IVDCDGModel(
  IVDCModel{FT}(;dt=dt,cʰ=cʰ,cᶻ=cᶻ,gconf=gconf),
  grid_3D,
  RusanovNumericalFlux(),
  CentralNumericalFluxSecondOrder(),
  CentralNumericalFluxGradient();
  direction=VerticalDirection(),
 );

 ivdc_Q   = init_ode_state(ivdc_dg, FT(0); init_on_cpu=true)
 ivdc_RHS = init_ode_state(ivdc_dg, FT(0); init_on_cpu=true)

 # Instantiate a batched GM-RES solver uaing balance law as its operator
 ivdc_bgm_solver=BatchedGeneralizedMinimalResidual(
   ivdc_dg,
   ivdc_Q;
   max_subspace_size=10);

 # Save initial state
 θ₀=deepcopy( ivdc_Q.θ )

 # Setup some callbacks and diagnostics
 diagsconf=(myivdcdg=ivdc_dg,myivdcq=ivdc_Q)
 mydiags_dg=MYDIAGSModel(
  MYDIAGSModel{FT}(;diagsconf=diagsconf),
  grid_3D,
  RusanovNumericalFlux(),
  CentralNumericalFluxSecondOrder(),
  CentralNumericalFluxGradient();
  direction=VerticalDirection(),
 );

 callbacks = ()
 struct EarthParameterSet <: AbstractEarthParameterSet end
 pset = EarthParameterSet()
 dgn_starttime = replace(string(now()), ":" => ".")
 od=ClimateMachine.ClimateMachine_Settings().output_dir
 Diagnostics.init(
   mpicomm,
   pset,
   mydiags_dg,
   ivdc_Q,
   dgn_starttime,
   od
  )

 # dgn_parms=(interval="1steps",type=IVDCConfigType(),name="default")
 dgngrp=setup_single_column_default_diagnostics( IVDCConfigType(), "1steps", "default_output_group")
 # dgngrp=setup_single_column_default_diagnostics(dgn_parms.type,dgn_parms.interval,dgn_parms.name)
 # dgn_config=ClimateMachine.DiagnosticsConfiguration([dgngrp])
 

 # Set up right hand side
 # for i=1:200000
 for i=1:50
  ivdc_RHS.θ   .= ivdc_Q.θ/dt
  #### ivdc_RHS.θ   .= ivdc_Q.θ
  ivdc_dg.state_auxiliary.θ_init .= ivdc_Q.θ

  println("Before maximum ", maximum(ivdc_Q.θ) )
  println("Before minimum ", minimum(ivdc_Q.θ) )
 
  # Evaluate operator
  #### ivdc_dg(ivdc_Q,ivdc_RHS,nothing,0;increment=false);
  #### println( maximum(ivdc_Q) )

  # Now try applying batched GM res solver
  lm!(y,x)=ivdc_dg(y,x,nothing,0;increment=false)
  solve_time = @elapsed iters = linearsolve!(lm!, ivdc_bgm_solver, ivdc_Q, ivdc_RHS);
  println("solver iters, time: ",iters, ", ", solve_time)

  println("After maximum ", maximum(ivdc_Q.θ) )
  println("After minimum ", minimum(ivdc_Q.θ) )
 end

 ## using Plots
 ## ni=2;ei=1
 ## nj=3;ej=1

 ## θᶠ=reshape( ivdc_Q.θ,(N+1,N+1,N+1,1,Nᶻ,Nˣ,Nʸ) )[ni,nj,:,1,:,ei,ej]
 ## θ⁰=reshape( θ₀,(N+1,N+1,N+1,1,Nᶻ,Nˣ,Nʸ) )[ni,nj,:,1,:,ei,ej]
 ## zc=reshape( grid_3D.vgeo,(N+1,N+1,N+1,16,Nᶻ,Nˣ,Nʸ) )[ni,nj,:,15,:,ei,ej]
 ## plot(θᶠ[:],zc[:],label="")
 ## plot!(θ⁰[:],zc[:],label="")
