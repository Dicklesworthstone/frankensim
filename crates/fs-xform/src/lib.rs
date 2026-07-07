//! fs-xform: design parameterizations — the levers ASCENT pulls
//! (plan §7.6). A parameterization is a differentiable map θ → space
//! warp, shipped with its Jacobian ACTION (`δT(x) = (∂T/∂θ)·δθ`, the
//! boundary velocity field shape-gradient assembly consumes, §8.7) and
//! its spatial Jacobian `∂T/∂x` (fold-over detection, composition chain
//! rule).
//!
//! v1 levers (all LINEAR in θ, so their Jacobian actions are exact basis
//! contractions — verified against dual-number JVPs and finite
//! differences in the conformance battery):
//! - [`ffd::FfdLattice`] — trivariate Bernstein free-form deformation.
//! - [`rbf::RbfMorph`] — compactly supported Wendland-C2 handle morphs.
//! - [`levelset::VelocityBand`] — narrow-band velocity DOFs + a working
//!   upwind SDF advection step (the topology-optimization workhorse;
//!   Appendix C's `xform.level-set-velocity`).
//! - [`density::DensityField`] — raw SIMP densities with clamp
//!   diagnostics (Helmholtz filtering is topo-simp's, downstream).
//! - [`Composed`] — parameterizations compose (global + local levers);
//!   Jacobian actions compose by the chain rule through the spatial
//!   Jacobian.
//!
//! Layer: L2 (MORPH). Runtime deps: `std`, fs-geom.

pub mod density;
pub mod ffd;
pub mod levelset;
pub mod rbf;

pub use density::DensityField;
pub use ffd::FfdLattice;
pub use fs_geom::{Point3, Vec3};
pub use levelset::{VelocityBand, advect_sdf};
pub use rbf::RbfMorph;

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured parameterization failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum XformError {
    /// The θ (or δθ) slice has the wrong length.
    DofMismatch {
        /// Expected DOF count.
        expected: usize,
        /// Supplied length.
        got: usize,
    },
    /// A parameter violates its validity bounds.
    OutOfBounds {
        /// Which component.
        index: usize,
        /// Its value.
        value: f64,
        /// Human-readable bound description.
        bound: &'static str,
    },
    /// The warp folds space (det ∂T/∂x ≤ 0) at a probed point.
    FoldOver {
        /// Where the fold was detected.
        at: Point3,
        /// The determinant found.
        det: f64,
    },
}

impl fmt::Display for XformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            XformError::DofMismatch { expected, got } => write!(
                f,
                "theta length {got} does not match the parameterization's {expected} DOFs"
            ),
            XformError::OutOfBounds {
                index,
                value,
                bound,
            } => {
                write!(f, "theta[{index}] = {value} violates {bound}")
            }
            XformError::FoldOver { at, det } => write!(
                f,
                "fold-over: det(dT/dx) = {det} <= 0 at ({}, {}, {}) — reduce the step or \
                 the lever amplitude",
                at.x, at.y, at.z
            ),
        }
    }
}

impl std::error::Error for XformError {}

/// The parameterization contract: a differentiable space warp `T_θ`.
///
/// Implementations guarantee: `apply` and `jacobian_action` agree in the
/// sense `T(θ + ε·δθ)(x) − T(θ)(x) = ε·jacobian_action(θ, δθ, x) + o(ε)`
/// (for the linear levers here, exactly, with no `o(ε)` term) — the G0
/// consistency law of the conformance battery.
pub trait Parameterization {
    /// Number of design DOFs.
    fn dof(&self) -> usize;

    /// The warped position `T_θ(x)`.
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on a wrong-length θ.
    fn apply(&self, theta: &[f64], x: Point3) -> Result<Point3, XformError>;

    /// The design velocity `δT(x) = (∂T/∂θ)(x) · δθ`.
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on wrong-length θ or δθ.
    fn jacobian_action(&self, theta: &[f64], dtheta: &[f64], x: Point3)
    -> Result<Vec3, XformError>;

    /// The spatial Jacobian `∂T/∂x` (rows = output components).
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on a wrong-length θ.
    fn spatial_jacobian(&self, theta: &[f64], x: Point3) -> Result<[[f64; 3]; 3], XformError>;
}

/// 3×3 determinant.
#[must_use]
pub fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Probe a warp for fold-over at sample points: the first probe with
/// `det(∂T/∂x) ≤ 0` is a structured refusal (invertibility monitoring).
///
/// # Errors
/// [`XformError::FoldOver`] naming the offending point and determinant.
pub fn detect_foldover(
    p: &dyn Parameterization,
    theta: &[f64],
    samples: &[Point3],
) -> Result<(), XformError> {
    for &x in samples {
        let j = p.spatial_jacobian(theta, x)?;
        let det = det3(&j);
        if det <= 0.0 {
            return Err(XformError::FoldOver { at: x, det });
        }
    }
    Ok(())
}

/// Composition `then ∘ first` (e.g. global manifold-harmonic lever + local
/// FFD near an inlet — the ornithoid setup). θ is the concatenation
/// `[θ_first | θ_then]`; the Jacobian action composes by the chain rule:
///
/// `δT(x) = δT_then(y) + (∂T_then/∂y)(y) · δT_first(x)`, `y = T_first(x)`.
pub struct Composed<'a> {
    /// Applied first.
    pub first: &'a dyn Parameterization,
    /// Applied to the output of `first`.
    pub then: &'a dyn Parameterization,
}

impl Composed<'_> {
    fn split<'t>(&self, theta: &'t [f64]) -> Result<(&'t [f64], &'t [f64]), XformError> {
        let (a, b) = (self.first.dof(), self.then.dof());
        if theta.len() != a + b {
            return Err(XformError::DofMismatch {
                expected: a + b,
                got: theta.len(),
            });
        }
        Ok(theta.split_at(a))
    }
}

impl Parameterization for Composed<'_> {
    fn dof(&self) -> usize {
        self.first.dof() + self.then.dof()
    }

    fn apply(&self, theta: &[f64], x: Point3) -> Result<Point3, XformError> {
        let (ta, tb) = self.split(theta)?;
        self.then.apply(tb, self.first.apply(ta, x)?)
    }

    fn jacobian_action(
        &self,
        theta: &[f64],
        dtheta: &[f64],
        x: Point3,
    ) -> Result<Vec3, XformError> {
        let (ta, tb) = self.split(theta)?;
        let (da, db) = self.split(dtheta)?;
        let y = self.first.apply(ta, x)?;
        let inner = self.first.jacobian_action(ta, da, x)?;
        let outer = self.then.jacobian_action(tb, db, y)?;
        let jb = self.then.spatial_jacobian(tb, y)?;
        let carried = Vec3::new(
            jb[0][0] * inner.x + jb[0][1] * inner.y + jb[0][2] * inner.z,
            jb[1][0] * inner.x + jb[1][1] * inner.y + jb[1][2] * inner.z,
            jb[2][0] * inner.x + jb[2][1] * inner.y + jb[2][2] * inner.z,
        );
        Ok(Vec3::new(
            outer.x + carried.x,
            outer.y + carried.y,
            outer.z + carried.z,
        ))
    }

    fn spatial_jacobian(&self, theta: &[f64], x: Point3) -> Result<[[f64; 3]; 3], XformError> {
        let (ta, tb) = self.split(theta)?;
        let y = self.first.apply(ta, x)?;
        let ja = self.first.spatial_jacobian(ta, x)?;
        let jb = self.then.spatial_jacobian(tb, y)?;
        let mut out = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    out[r][c] += jb[r][k] * ja[k][c];
                }
            }
        }
        Ok(out)
    }
}
