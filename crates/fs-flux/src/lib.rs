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

pub mod ale;
pub mod bdm;
pub mod ns;
pub mod reduced_aero;
pub mod trimesh;

pub use ns::{FluxParams, FluxSolution, FluxSystem};
pub use reduced_aero::{
    AlternativeWrenchSet, ApplicabilityEnvelope, BodyKinematics, CandidateWrench, ClosedRange,
    ComponentWrenches, ContributionFamily, CorrelationIdentity, CorrelationUncertainty,
    DiscGeometry, DiscPose, EdgeFlow, EstimateAuthority, FormDrag, GasProperties, GasPropertyCard,
    OrientationRateDamping, ReducedAeroComponents, ReducedAeroError, ReducedAeroInput,
    ReducedAeroModel, RotationalSkinFriction, SurfaceRoughness, Vec3, WorkReceipt, WorkWindow,
};
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
