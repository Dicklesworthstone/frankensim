//! fs-ga: the geometric-algebra layer (plan §7.7, Bet 2). PGA Cl(3,0,1)
//! as the kinematics substrate — motors/screws replace the quaternion +
//! matrix + Plücker zoo and kill gimbal-class bugs BY CONSTRUCTION — and
//! CGA Cl(4,1) for sphere/tangency-rich construction. All multiplication
//! tables are CONST-EVALUATED from the metric signatures (P2-deterministic
//! fixed-order products, no runtime blade bookkeeping); conventional
//! Vec3/quaternion/matrix façades sit at the API boundary so no caller
//! pays a formalism tax.
//!
//! Layer: L2 (MORPH). Runtime deps: `std`, fs-math (deterministic trig
//! for exp/log so motors are bit-identical across platforms).

pub mod cga;
pub mod facade;
pub mod lie;
pub mod mv;
pub mod pga;
pub mod table;

pub use facade::{Mat34, Quat, Vec3};
pub use lie::{
    Mat3, Mat6, Se3, Se3Jacobian, So3, So3Tangent, Twist, Wrench, so3_left_jacobian,
    so3_left_jacobian_inverse, so3_right_jacobian, so3_right_jacobian_inverse,
};
pub use mv::{Cga, Pga};
pub use pga::{Line, Motor, Plane, Point, exp_bivector, motor_log};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured geometric-algebra failures (Decalogue P10).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GaError {
    /// A trivector with zero e123 weight — an ideal (at-infinity) point
    /// with no Cartesian form.
    IdealPoint,
    /// A conformal element with no finite representative.
    ZeroWeight {
        /// Which operation refused.
        context: &'static str,
    },
    /// One coordinate was NaN or infinite.
    NonFinite {
        /// Which operation refused the value.
        context: &'static str,
        /// Coordinate index in that operation's documented ordering.
        index: usize,
    },
    /// A norm required to be nonzero collapsed to a degenerate value.
    DegenerateNorm {
        /// Which representation was degenerate.
        context: &'static str,
        /// Observed squared norm.
        norm_squared: f64,
    },
    /// A group representative failed its unit-versor invariant.
    NotUnit {
        /// Which representation was invalid.
        context: &'static str,
        /// Observed absolute unit defect.
        defect: f64,
        /// Largest accepted defect.
        tolerance: f64,
    },
    /// A representation contains coefficients forbidden by its grade/type.
    InvalidRepresentation {
        /// Which representation was invalid.
        context: &'static str,
        /// Largest forbidden coefficient.
        defect: f64,
        /// Largest accepted coefficient.
        tolerance: f64,
    },
    /// An otherwise valid calculation is too ill-conditioned to certify.
    IllConditioned {
        /// Which calculation refused.
        context: &'static str,
        /// Deterministic conditioning indicator (smaller or larger is
        /// documented by the producing operation).
        measure: f64,
        /// Refusal threshold for that indicator.
        limit: f64,
    },
    /// A bounded deterministic series could not certify its requested tail.
    SeriesDidNotConverge {
        /// Which series refused.
        context: &'static str,
        /// Number of terms accumulated.
        terms: usize,
        /// Last analytic tail bound.
        tail_bound: f64,
    },
}

impl fmt::Display for GaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GaError::IdealPoint => {
                write!(f, "ideal (at-infinity) point has no Cartesian coordinates")
            }
            GaError::ZeroWeight { context } => write!(f, "zero-weight element: {context}"),
            GaError::NonFinite { context, index } => {
                write!(f, "non-finite coordinate {index} in {context}")
            }
            GaError::DegenerateNorm {
                context,
                norm_squared,
            } => write!(
                f,
                "degenerate norm in {context}: squared norm {norm_squared}"
            ),
            GaError::NotUnit {
                context,
                defect,
                tolerance,
            } => write!(
                f,
                "non-unit {context}: defect {defect} exceeds tolerance {tolerance}"
            ),
            GaError::InvalidRepresentation {
                context,
                defect,
                tolerance,
            } => write!(
                f,
                "invalid {context}: forbidden-component defect {defect} exceeds tolerance {tolerance}"
            ),
            GaError::IllConditioned {
                context,
                measure,
                limit,
            } => write!(
                f,
                "ill-conditioned {context}: measure {measure}, refusal limit {limit}"
            ),
            GaError::SeriesDidNotConverge {
                context,
                terms,
                tail_bound,
            } => write!(
                f,
                "{context} did not converge after {terms} terms (tail bound {tail_bound})"
            ),
        }
    }
}

impl std::error::Error for GaError {}
