# Linear model equations, for ocean model implicit vertical diffusion
# convective adjustment step.
#
using ClimateMachine.DGMethods.NumericalFluxes:
    NumericalFluxFirstOrder,
    NumericalFluxSecondOrder,
    NumericalFluxGradient

using ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient,
    CentralNumericalFluxSecondOrder

using LinearAlgebra: I, dot, Diagonal

using ClimateMachine.BalanceLaws:
   BalanceLaw,
   Prognostic,
   Auxiliary,
   Gradient,
   GradientFlux

using ClimateMachine.VariableTemplates

import ClimateMachine.DGMethods:
     init_state_auxiliary!,
     update_auxiliary_state!,
     update_auxiliary_state_gradient!,
     vars_state,
     VerticalDirection,
     boundary_state!,
     compute_gradient_flux!,
     init_state_prognostic!,
     flux_first_order!,
     flux_second_order!,
     source!,
     wavespeed,
     compute_gradient_argument!

using StaticArrays, Random

"""
    IVDCModel{M} <: BalanceLaw

    This code defines DG `BalanceLaw` terms for an operator, L, that is evaluated from iterative
    implicit solver to solve an equation of the form

    (L + 1/Δt) ϕ^{n+1} = ϕ^{n}/Δt

    where L is a vertical diffusion operator with a spatially varying diffusion
    coefficient.
"""
struct IVDCModel{FT, G} <: BalanceLaw
      dt :: FT
      κᶻ :: FT
      κᶜ :: FT
      cʰ :: FT
      cᶻ :: FT
   gconf :: G

   function IVDCModel{FT}(; dt = FT(5400),
                            κᶻ = FT(1e-4),
                            κᶜ = FT(0.1),
                            cʰ = FT(0),
                            cᶻ = FT(0),
                            gconf = ()
                          ) where FT

       return new{FT}(dt, κᶻ, κᶜ, cʰ, cᶻ, gconf)
   end
end

function IVDCDGModel(balance_law, grid, nfnondiff, nfdiff, gnf; kwargs...)

   model_data = (
                 dt = balance_law.dt,
                 κᶻ = balance_law.κᶻ,
                 κᶜ = balance_law.κᶜ,
                 cʰ = balance_law.cʰ,
                 cᶻ = balance_law.cᶻ,
                 gconf = balance_law.gconf
                )

   return DGModel(balance_law, grid, nfnondiff, nfdiff, gnf; modeldata=model_data, kwargs...)
end

"""
 Set model state variables and operators
"""

# State variable and initial value, just one for now, θ
##
vars_state(m::IVDCModel, ::Prognostic, FT) = @vars(θ::FT)

init_state_prognostic!( args...) = ( init_state_cp!( args... ) )
function init_state_cp!( m::IVDCModel, Q::Vars, A::Vars, coords, t,)
  @inbounds begin
    x = coords[1]
    y = coords[2]
    z = coords[3]
    # Wavy with stable vertical stratification
    # Q.θ = (5 + 4 * cos(y * π / p.Lʸ)) * (1 + z / H)

    # Constant everywhere
    # Q.θ = 9

    # Big number
    # Q.θ = 50

    # Unstable column
    #
    # Define function with inversion near surface and stable stratification below
    #
    # $ \theta(z) = Ae^{\frac{-(z+L)}{L}} - Ce^{\frac{-D(z+L)}{L}} + B + Ez$
    #
    H₀=m.gconf.H/4
    H₁=m.gconf.H*2
    ΔH=H₁-H₀
    sx=x/m.gconf.Lˣ
    sy=y/m.gconf.Lʸ
    D₀=2.
    D₁=3.
    ΔD=D₁-D₀
    # rz=rand(length(z)).-0.5
    rz=rand()-0.5
    H=H₀+sx*ΔH
    D=D₀+sy*ΔD
    L=H/10.;A=20.;C=50.;B=8.;E=5.e-4;
    ft(xx,L)=exp(-xx/L);
    th1(zz)=A*ft.(-zz .+ L, L);
    th2(zz)=C*ft(D*(-zz .+ L), L);
    phi1=th1.(z)
    phi2=th2.(z)
    val=phi1 .- phi2 .+ B .+ E .* z;
    # Q.θ = val .+ 0.005.*( val.*rz )
    Q.θ = val + 0.005*val*rz
  end
  return nothing
end
##

##
vars_state(m::IVDCModel, ::Auxiliary, FT) = @vars(θ_init::FT)

function init_state_auxiliary!(m::IVDCModel,A::Vars, _...)
  @inbounds begin
    A.θ_init = -0
  end
  return nothing
end
##
#

# Variables and operations used in differentiating first derivatives
##
vars_state(m::IVDCModel, ::Gradient, FT) = @vars(∇θ::FT, ∇θ_init::FT,)
@inline function compute_gradient_argument!(
    m::IVDCModel,
    G::Vars,
    Q::Vars,
    A,
    t,
)
    G.∇θ = Q.θ
    G.∇θ_init = A.θ_init

    return nothing
end
##
#

# Variables and operations used in differentiating second derivatives
##
vars_state(m::IVDCModel, ::GradientFlux, FT) = @vars(κ∇θ::SVector{3, FT})
@inline function compute_gradient_flux!(
    m::IVDCModel,
    D::Vars,
    G::Grad,
    Q::Vars,
    A::Vars,
    t,
)
    κ = diffusivity_tensor(m, G.∇θ_init[3])
    D.κ∇θ = -κ * G.∇θ 
    return nothing
end
##

##
## Set vertical diffusivity profile based on vertical hydrography profile
@inline function diffusivity_tensor(m::IVDCModel, ∂θ∂z)
    κᶻ = m.κᶻ
    κᶜ = m.κᶜ
    ∂θ∂z < 0 ? κ = (@SVector [0, 0, κᶜ]) : κ = (@SVector [0, 0, κᶻ])
    return Diagonal(-κ)
end
##
#

# Function to apply I to state variable

##
@inline function source!(
    m::IVDCModel,
    S::Vars,
    Q::Vars,
    D::Vars,
    A::Vars,
    t,
    direction,
)
    ivdc_dt = m.dt
    @inbounds begin
     S.θ = Q.θ/ivdc_dt
    end

    return nothing
end
##
#

# Numerical fluxes and boundaries

##
function flux_first_order!(::IVDCModel, _...) end

function flux_second_order!(
    m::IVDCModel,
    F::Grad,
    Q::Vars,
    D::Vars,
    H::Vars,
    A::Vars,
    t,
)
    F.θ += D.κ∇θ
end

function wavespeed(m::IVDCModel, n⁻, _...)
    C = abs(SVector(m.cʰ, m.cʰ, m.cᶻ)' * n⁻)
    return C
end
##

function boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient, CentralNumericalFluxGradient},
    m::IVDCModel,
    Q⁺,
    A⁺,
    n,
    Q⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    return nothing
end

###    From -  function numerical_boundary_flux_gradient! , DGMethods/NumericalFluxes.jl
###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        aux1⁻,
###    )

function boundary_state!(
    nf::Union{NumericalFluxSecondOrder,CentralNumericalFluxSecondOrder},
    m::IVDCModel,
    Q⁺,
    D⁺,
    A⁺,
    n⁻,
    Q⁻,
    D⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    D⁺.κ∇θ = n⁻ * -0
    # D⁺.κ∇θ = Q/(rho*cp)   ( at top boundary for heat flux )
    return nothing
end

###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

###    boundary_flux_second_order!(
###        numerical_flux,
###        balance_law,
###        Grad{S}(flux),
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_hyperdiffusive⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_hyperdiffusive⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

