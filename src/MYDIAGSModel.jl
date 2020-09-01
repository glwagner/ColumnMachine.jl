# Linear model equations, for ocean model implicit vertical diffusion
# convective adjustment step.
#
using ClimateMachine.DGMethods.NumericalFluxes:
    NumericalFluxFirstOrder, NumericalFluxSecondOrder, NumericalFluxGradient

using ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient, CentralNumericalFluxSecondOrder

using LinearAlgebra: I, dot, Diagonal

using ClimateMachine.BalanceLaws:
   BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux

using ClimateMachine.VariableTemplates

import ClimateMachine.DGMethods:
     init_state_auxiliary!, update_auxiliary_state!, update_auxiliary_state_gradient!, vars_state, VerticalDirection, boundary_state!, compute_gradient_flux!, init_state_prognostic!, flux_first_order!, flux_second_order!, source!, wavespeed, compute_gradient_argument!

using StaticArrays
using Random

"""
 MYDIAGSModel{M} <: BalanceLaw

 DG model container for my diagnostics

 # Usage

"""

# Create a new linear model instance
abstract type AbstractMYDIAGSModel <: BalanceLaw end
struct MYDIAGSModel{FT} <: AbstractMYDIAGSModel 
 diagsconf
 function MYDIAGSModel{FT}(
  ;
  diagsconf=()
 ) where {FT<: AbstractFloat}
    return new{FT}(
      diagsconf
    )
 end
end

function MYDIAGSModel(
   bl::MYDIAGSModel,
   grid,
   nfnondiff,
   nfdiff,
   gnf;
   kwargs...,
   )

   modeldata=(diagsconf=bl.diagsconf)

   return DGModel(bl,grid,nfnondiff,nfdiff,gnf;kwargs...,modeldata=modeldata,)
end
