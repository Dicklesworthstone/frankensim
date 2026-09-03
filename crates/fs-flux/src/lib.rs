//! fs-flux — incompressible Navier–Stokes, FEEC-native (plan §8.3
//! [F], bead tfz.17): H(div)-conforming BDM1 velocities with P0
//! pressures — EXACTLY divergence-free discrete velocities, so
//! velocity errors are independent of the pressure (the de Rham
//! exactness cashing out as PRESSURE-ROBUSTNESS, the correctness
//! property most production codes lack). Interior-penalty viscosity
//! (jumps are purely tangential by conformity), upwinded DG
//! convection on the single-valued face flux w·n, Picard steady
//! solves, IMEX BDF1 transients, and discrete adjoints. 2D
//! triangle-mesh instantiation; 3D, BDM2+, projection time stepping,
//! and LES closures are recorded successors with honesty labels —
//! no turbulence model ships here, and nothing pretends otherwise.

#[cfg(feature = "continuum")]
pub mod ale;
#[cfg(feature = "continuum")]
pub mod bdm;
/// Reduced circular-capillary flow screening with an explicit Newtonian,
/// fully-developed applicability boundary.
pub mod capillary;
#[cfg(feature = "continuum")]
pub mod gas_film;
#[cfg(feature = "quarter-wave")]
pub mod lc;
#[cfg(feature = "continuum")]
pub mod ns;
#[cfg(feature = "continuum")]
pub mod reduced_aero;
#[cfg(feature = "continuum")]
pub mod trimesh;

pub use capillary::{
    CircularCapillaryError, CircularCapillaryInput, CircularCapillaryStep,
    step_newtonian_circular_capillary,
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
pub use ns::{FluxParams, FluxSolution, FluxSystem};
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
