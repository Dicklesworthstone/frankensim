//! fs-solid — the elasticity core (plan §8.2, bead tfz.13): linear
//! elasticity through finite-strain hyperelasticity with a measured
//! locking-free path, on BOTH the body-fitted and the CutFEM-on-SDF
//! frontends (topology optimization solves THIS physics on octrees).
//!
//! Layer: L3. Constitutive laws live in fs-material (energies, exact
//! AD stresses, consistent tangents); this crate owns KINEMATICS,
//! WEAK FORMS, and NEWTON LOOPS:
//! - [`mesh2`]: structured 2D body-fitted meshes — P1 triangles, Q1
//!   quads, and mapped quadrilateral panels (Cook's membrane).
//! - [`linear`]: plane-strain/plane-stress small-strain elasticity;
//!   standard displacement elements plus the B-BAR dilatation
//!   projection (the measured locking-free formulation).
//! - [`hyper2d`]: plane-strain finite-strain hyperelasticity through
//!   fs-material cards (Neo-Hookean, Mooney–Rivlin): exact residuals
//!   and consistent tangents from the 3D deformation gradient,
//!   Newton with backtracking line search and load stepping.
//! - [`cutfront`]: the CutFEM frontend — vector Q1 on fs-cutfem
//!   background quadtrees, symmetric Nitsche displacement conditions,
//!   componentwise ghost penalty.
//!
//! Element-selection guidance (fs-regime's structural indicators feed
//! these thresholds): [`select_formulation`] returns B-bar whenever
//! near-incompressibility (ν ≥ 0.45) or bending domination
//! (slenderness ≥ 5) threatens locking — the conformance battery
//! MEASURES the standard element's failure on exactly those regimes.
//! TDNNS-proper (normal-normal-stress continuity) awaits the
//! simplicial H(div)-family bead and is a recorded no-claim.

pub mod continuation;
pub mod cutfront;
pub mod hyper2d;
#[cfg(feature = "koiter-asymptotics")]
pub mod koiter;
pub mod linear;
pub mod mesh2;
pub mod stability;

pub use continuation::{ArcSettings, PathEvent, PathResidual, PathState, advance, switch_branch};
pub use cutfront::{CutElasticity, CutSolution};
pub use hyper2d::{HyperProblem, NewtonReport, NewtonSettings};
pub use linear::{Formulation, LinearProblem, PlaneKind};
pub use mesh2::{Mesh2, Patch};
pub use stability::{
    BucklingResult, buckling_loads, eigenvalue_derivative, expand_mode, group_stiffness,
    ks_aggregate, ks_aggregate_derivative, lambda_indicator, reduced_pencil,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Teaching errors (P10).
#[derive(Debug, Clone, PartialEq)]
pub enum SolidError {
    /// The linear solve missed its residual gate. Repair: refine the
    /// load step, enable a stabilized formulation, or raise the
    /// iteration cap.
    SolveFailed {
        /// Iterations performed.
        iters: usize,
        /// Relative residual reached.
        rel_residual: f64,
    },
    /// Newton did not converge within the iteration budget. The
    /// carried history names the last residual norms; the repair is
    /// more load steps (the fixture battery demonstrates the recipe).
    NewtonStalled {
        /// Residual norms per iteration.
        history: Vec<f64>,
    },
    /// The material card refused a state (typically det F ≤ 0 under a
    /// too-aggressive step). Line search normally absorbs this; if it
    /// reaches you, the load step is too large.
    MaterialRefused {
        /// The material's own message.
        what: String,
    },
    /// A boundary condition named a patch the mesh does not carry.
    UnknownPatch {
        /// The requested patch.
        patch: Patch,
    },
    /// An fs-scenario boundary condition is outside this crate's
    /// elasticity surface (wrong physics or an unsupported kind).
    UnsupportedBc {
        /// What was requested.
        what: String,
    },
}

impl core::fmt::Display for SolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SolidError::SolveFailed {
                iters,
                rel_residual,
            } => write!(
                f,
                "linear solve stalled at {rel_residual:.3e} after {iters} \
                 iterations; refine load steps or switch formulation"
            ),
            SolidError::NewtonStalled { history } => write!(
                f,
                "Newton stalled after {} iterations (last residual \
                 {:.3e}); use more load steps",
                history.len(),
                history.last().copied().unwrap_or(f64::NAN)
            ),
            SolidError::MaterialRefused { what } => write!(
                f,
                "material card refused the state ({what}); reduce the \
                 load step"
            ),
            SolidError::UnknownPatch { patch } => {
                write!(f, "mesh carries no patch {patch:?}")
            }
            SolidError::UnsupportedBc { what } => write!(
                f,
                "fs-scenario condition outside the elasticity surface: {what}"
            ),
        }
    }
}

impl std::error::Error for SolidError {}

/// Structural-regime indicators (fs-regime supplies these upstream;
/// the thresholds here are the documented element-selection policy).
#[derive(Debug, Clone, Copy)]
pub struct RegimeIndicators {
    /// Poisson ratio of the dominant material.
    pub poisson: f64,
    /// Characteristic length over thickness (bending domination).
    pub slenderness: f64,
}

/// Documented element-selection guidance: standard displacement
/// elements until near-incompressibility (ν ≥ 0.45) or bending
/// domination (slenderness ≥ 5) — then the B-bar projection. The
/// locking battery measures exactly this boundary.
#[must_use]
pub fn select_formulation(r: RegimeIndicators) -> Formulation {
    if r.poisson >= 0.45 || r.slenderness >= 5.0 {
        Formulation::BBar
    } else {
        Formulation::Standard
    }
}

/// Map an fs-scenario boundary condition onto this crate's elasticity
/// surface: `Dirichlet` and `Traction` under `Physics::Elasticity`
/// are accepted (the caller resolves the dimensioned value to
/// numbers — units plumbing is the scenario consumer's job).
///
/// # Errors
/// [`SolidError::UnsupportedBc`] for any other (physics, kind) pair.
pub fn accept_scenario_bc(
    bc: &fs_scenario::bc::BoundaryCondition,
) -> Result<fs_scenario::bc::BcKind, SolidError> {
    use fs_scenario::bc::{BcKind, Physics};
    if bc.physics != Physics::Elasticity {
        return Err(SolidError::UnsupportedBc {
            what: format!("physics {:?} (this crate is Elasticity)", bc.physics),
        });
    }
    match bc.kind {
        BcKind::Dirichlet | BcKind::Traction => Ok(bc.kind),
        other => Err(SolidError::UnsupportedBc {
            what: format!("kind {other:?} (supported: Dirichlet, Traction)"),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn selection_policy_matches_doc() {
        let std_ = select_formulation(RegimeIndicators {
            poisson: 0.3,
            slenderness: 1.0,
        });
        let inc = select_formulation(RegimeIndicators {
            poisson: 0.4999,
            slenderness: 1.0,
        });
        let thin = select_formulation(RegimeIndicators {
            poisson: 0.3,
            slenderness: 10.0,
        });
        assert_eq!(std_, Formulation::Standard);
        assert_eq!(inc, Formulation::BBar);
        assert_eq!(thin, Formulation::BBar);
    }
}
