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
pub mod mv;
pub mod pga;
pub mod table;

pub use facade::{Mat34, Quat, Vec3};
pub use mv::{Cga, Pga};
pub use pga::{Line, Motor, Plane, Point, exp_bivector, motor_log};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured geometric-algebra failures (Decalogue P10).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaError {
    /// A trivector with zero e123 weight — an ideal (at-infinity) point
    /// with no Cartesian form.
    IdealPoint,
    /// A conformal element with no finite representative.
    ZeroWeight {
        /// Which operation refused.
        context: &'static str,
    },
}

impl fmt::Display for GaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GaError::IdealPoint => {
                write!(f, "ideal (at-infinity) point has no Cartesian coordinates")
            }
            GaError::ZeroWeight { context } => write!(f, "zero-weight element: {context}"),
        }
    }
}

impl std::error::Error for GaError {}
