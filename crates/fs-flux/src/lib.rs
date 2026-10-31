//! fs-flux — incompressible Navier–Stokes, FEEC-native (plan §8.3
//! [F], bead tfz.17): H(div)-conforming BDM1 velocities with P0
//! pressures — EXACTLY divergence-free discrete velocities, so
//! velocity errors are independent of the pressure (the de Rham
//! exactness cashing out as PRESSURE-ROBUSTNESS, the correctness
//! property most production codes lack). Interior-penalty viscosity
//! (jumps are purely tangential by conformity), upwinded DG
//! convection on the single-valued face flux w·n, Picard steady
//! solves, IMEX and fully implicit rheological BDF1 transients, and discrete
//! adjoints. 2D
//! triangle-mesh instantiation; 3D, BDM2+, projection time stepping,
//! and LES closures are recorded successors with honesty labels —
//! no turbulence model ships here, and nothing pretends otherwise.

#[cfg(feature = "continuum")]
pub mod ale;
#[cfg(feature = "continuum")]
pub mod bdm;
/// Newtonian circular-tube flow and gravity-aware capillary meniscus motion,
/// with explicit quasisteady and optional mean-flow inertial applicability.
pub mod capillary;
/// Material-driven spherical bubble dynamics with bounded adaptive stepping.
#[cfg(feature = "cavitation")]
pub mod cavitation;
/// Source-owned, conservative two-dimensional shock propagation and reflection.
#[cfg(feature = "experimental-compressible")]
pub mod compressible;
#[cfg(feature = "continuum")]
pub mod gas_film;
#[cfg(feature = "quarter-wave")]
pub mod lc;
/// Material-bound FEEC flow with coherent SI loads, pressure and kinetic energy.
#[cfg(feature = "continuum")]
#[path = "ns/material.rs"]
pub mod material;
#[cfg(feature = "continuum")]
pub mod ns;
/// Conservative isothermal solution transport through selective membranes.
#[cfg(feature = "continuum")]
pub mod osmosis;
#[cfg(feature = "continuum")]
pub mod reduced_aero;
/// Reciprocal quasistatic Reynolds pressure for deforming thin-gap interfaces.
#[path = "gas_film/resistive.rs"]
pub mod resistive_film;
#[cfg(feature = "continuum")]
pub mod trimesh;

pub use capillary::inertia::{
    InertialCapillaryEnd, InertialCapillaryError, InertialCapillaryInput, InertialCapillarySolve,
    InertialCapillaryStep, InertialContactLine, step_inertial_circular_meniscus,
};
pub use capillary::{
    CircularCapillaryError, CircularCapillaryInput, CircularCapillaryStep, CircularMeniscusEnd,
    CircularMeniscusInput, CircularMeniscusStep, ContactLineCircularMeniscusInput,
    ContactLineCircularMeniscusStep, ContactLineMotion, InclinedCircularMeniscusInput,
    InclinedCircularMeniscusStep, InclinedTwoFluidCircularMeniscusInput,
    InclinedTwoFluidCircularMeniscusStep, TwoFluidCircularMeniscusInput,
    TwoFluidCircularMeniscusStep, step_newtonian_circular_capillary,
    step_newtonian_circular_meniscus, step_newtonian_contact_line_circular_meniscus,
    step_newtonian_inclined_circular_meniscus, step_newtonian_inclined_two_fluid_circular_meniscus,
    step_newtonian_two_fluid_circular_meniscus,
};
#[cfg(feature = "continuum")]
pub use gas_film::{
    ContactExclusionMask, GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget,
    GasFilmCheckpoint, GasFilmError, GasFilmGrid1d, GasFilmIdentity, GasFilmInput,
    GasFilmInputAuthority, GasFilmReceipt, GasFilmStep, GasFilmUncertainty, IsothermalIdealGas,
    MovingWallInput, RoughnessPolicy, SlipPolicy, isothermal_compressible_reynolds_model_id,
    solve_isothermal_gas_film_1d,
};
#[cfg(feature = "continuum")]
pub use material::{
    MaterialBrinkmanSolution, MaterialFluxError, MaterialFluxNumerics, MaterialFluxRheologyPoint,
    MaterialFluxRheologySample, MaterialFluxRheologyStep, MaterialFluxSolution, MaterialFluxSystem,
    MaterialStokesBoundary, MaterialStokesRheologySolution,
};
#[cfg(feature = "continuum")]
pub use ns::adjoint::{FluxAdjointError, FluxAdjointSolution};
#[cfg(feature = "continuum")]
pub use ns::boundary::{
    StokesAssembly, StokesBoundary, StokesBoundaryError, StokesPressureReference,
    StokesRheologyControls, StokesRheologySolution, StokesSolution,
};
#[cfg(feature = "continuum")]
pub use ns::linearization::{
    FluxLinearizationError, FluxRheologyLinearization, FluxViscosityTangent,
};
#[cfg(feature = "continuum")]
pub use ns::newton::{FluxNewtonControls, FluxNewtonError, FluxNewtonStep};
#[cfg(feature = "continuum")]
pub use ns::porous::{BrinkmanAssembly, BrinkmanCell};
#[cfg(feature = "continuum")]
pub use ns::{FluxParams, FluxRheologyError, FluxRheologyStep, FluxSolution, FluxSystem};
#[cfg(feature = "continuum")]
pub use reduced_aero::{
    AlternativeWrenchSet, ApplicabilityEnvelope, BodyKinematics, CandidateWrench, ClosedRange,
    ComponentWrenches, ContributionFamily, CorrelationIdentity, CorrelationUncertainty,
    DiscGeometry, DiscPose, EdgeFlow, EstimateAuthority, FormDrag, GasProperties, GasPropertyCard,
    OrientationRateDamping, ReducedAeroComponents, ReducedAeroError, ReducedAeroInput,
    ReducedAeroModel, RotationalSkinFriction, SurfaceRoughness, Vec3, WorkReceipt, WorkWindow,
};
#[cfg(feature = "continuum")]
pub use trimesh::TriMesh;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
