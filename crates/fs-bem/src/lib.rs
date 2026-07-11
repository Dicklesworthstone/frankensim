//! fs-bem — Laplace BEM panel methods (plan §8.3 [F], bead tfz.20):
//! potential-flow screening for exterior aerodynamics, O(N)-class via
//! fs-fmm. Inviscid honesty labels apply everywhere: this is the
//! ornithoid flagship's WIDE-SEARCH stage, not a viscous truth source.
//!
//! Layer: L3.
//! - [`panel3d`]: 3D exterior flow — constant source panels on
//!   fs-rep-mesh surfaces, collocation Neumann conditions, GMRES with
//!   the FMM-accelerated gradient matvec (three Chebyshev passes for
//!   the vector kernel, dotted with target normals) and the dense
//!   direct path as the oracle; the sphere's analytic surface speed is
//!   the G2 gate, single-layer reciprocity the G0 identity.
//! - [`panel2d`]: 2D Hess–Smith airfoils — constant sources per panel
//!   plus one bound vortex, the KUTTA condition closing the system;
//!   thin-airfoil lift slope as the reference band; the ADJOINT
//!   (one transposed dense solve) is the committed gradient path,
//!   FD-gated.
//! - [`wake2d`]: unsteady free wakes — point-vortex sheets shed at the
//!   trailing edge with Kutta-determined strength (Kelvin circulation
//!   conservation), convected by the regularized induced flow; the
//!   impulsive-start fixture shows the Wagner-like circulation
//!   transient and stable roll-up.

use fs_fmm::FmmError;
use fs_la::factor::FactorError;

/// Invalid geometry, an inadmissible work request, or a numerical refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum BemError {
    /// A panel count is outside the admitted range or violates a constructor's
    /// parity requirement.
    InvalidPanelCount {
        /// Requested count.
        count: usize,
        /// Minimum count.
        min: usize,
        /// Maximum count.
        max: usize,
        /// Whether this constructor requires an even count.
        even_required: bool,
    },
    /// A scalar input was outside its documented physical/numerical domain.
    InvalidScalar {
        /// Input name.
        name: &'static str,
        /// Rejected value.
        value: f64,
        /// Concise admissibility rule.
        requirement: &'static str,
    },
    /// An airfoil node coordinate was invalid.
    InvalidNode {
        /// Node index.
        index: usize,
        /// Coordinate axis.
        axis: usize,
        /// Rejected value.
        value: f64,
    },
    /// A consecutive airfoil edge was too short relative to the section scale.
    DegeneratePanel {
        /// Panel index.
        index: usize,
    },
    /// The airfoil polygon enclosed no numerically meaningful area.
    DegenerateAirfoil,
    /// Airfoil nodes were not in the required clockwise order.
    WrongAirfoilOrientation,
    /// Two non-adjacent airfoil panels intersected.
    SelfIntersectingAirfoil {
        /// First panel index.
        first: usize,
        /// Second panel index.
        second: usize,
    },
    /// Parallel public panel vectors had inconsistent lengths.
    PanelDataLength {
        /// Vector name.
        field: &'static str,
        /// Expected length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },
    /// One 3D panel datum violated the surface contract.
    InvalidSurfacePanel {
        /// Panel index.
        index: usize,
        /// Invalid datum.
        field: &'static str,
    },
    /// A vector argument did not match the operator dimension.
    VectorLength {
        /// Argument name.
        field: &'static str,
        /// Expected length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },
    /// An explicitly bounded dense or transient work request was too large.
    WorkEnvelopeExceeded {
        /// Operation name.
        operation: &'static str,
        /// Requested logical elements/steps.
        requested: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// A bounded allocation could not be reserved.
    AllocationFailed {
        /// Operation being prepared.
        operation: &'static str,
    },
    /// Dense factorization refused a singular or invalid system.
    LinearSolve {
        /// Primal or adjoint stage.
        stage: &'static str,
        /// Factorization refusal.
        source: FactorError,
    },
    /// The FMM layer refused its inputs or work envelope.
    Fmm(FmmError),
    /// A supposedly successful numerical path produced a non-finite value.
    NonFiniteResult {
        /// Operation that produced it.
        operation: &'static str,
    },
    /// Trace downsampling requires a nonzero stride.
    InvalidTraceStride,
}

impl core::fmt::Display for BemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidPanelCount {
                count,
                min,
                max,
                even_required,
            } => write!(
                f,
                "panel count {count} is invalid; expected {min}..={max}{}",
                if *even_required { " and even" } else { "" }
            ),
            Self::InvalidScalar {
                name,
                value,
                requirement,
            } => write!(f, "{name}={value} is invalid; expected {requirement}"),
            Self::InvalidNode { index, axis, value } => {
                write!(
                    f,
                    "airfoil node {index} coordinate {axis} is invalid ({value})"
                )
            }
            Self::DegeneratePanel { index } => {
                write!(
                    f,
                    "airfoil panel {index} is degenerate at the section scale"
                )
            }
            Self::DegenerateAirfoil => f.write_str("airfoil polygon encloses no valid area"),
            Self::WrongAirfoilOrientation => f.write_str("airfoil nodes must be ordered clockwise"),
            Self::SelfIntersectingAirfoil { first, second } => {
                write!(f, "airfoil panels {first} and {second} intersect")
            }
            Self::PanelDataLength {
                field,
                expected,
                actual,
            }
            | Self::VectorLength {
                field,
                expected,
                actual,
            } => write!(
                f,
                "{field} length {actual} does not match expected length {expected}"
            ),
            Self::InvalidSurfacePanel { index, field } => {
                write!(f, "surface panel {index} has invalid {field}")
            }
            Self::WorkEnvelopeExceeded {
                operation,
                requested,
                max,
            } => write!(
                f,
                "{operation} work request {requested} exceeds the admitted maximum {max}"
            ),
            Self::AllocationFailed { operation } => {
                write!(f, "failed to reserve bounded storage for {operation}")
            }
            Self::LinearSolve { stage, source } => write!(f, "{stage} solve failed: {source}"),
            Self::Fmm(source) => write!(f, "FMM refused BEM operation: {source}"),
            Self::NonFiniteResult { operation } => {
                write!(f, "{operation} produced a non-finite result")
            }
            Self::InvalidTraceStride => f.write_str("wake trace stride must be nonzero"),
        }
    }
}

impl std::error::Error for BemError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fmm(source) => Some(source),
            _ => None,
        }
    }
}

impl From<FmmError> for BemError {
    fn from(source: FmmError) -> Self {
        Self::Fmm(source)
    }
}

pub mod panel2d;
pub mod panel3d;
pub mod wake2d;

pub use panel2d::{Airfoil2d, PanelSolution2d, naca4_symmetric};
pub use panel3d::{ExteriorSolution, ExteriorSolveError, SpherePanels, solve_exterior};
pub use wake2d::{WakeSim, WakeStep};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
