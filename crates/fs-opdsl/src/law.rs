//! Pointwise differentiable laws: the opaque-but-differentiable nodes
//! constitutive maps enter through. A law supplies its value AND its
//! derivative from one definition; the plan's chain rule consumes the
//! derivative for both JVP and VJP (a pointwise diagonal is its own
//! transpose), and the battery verifies it mechanically against
//! fs-ad's dual numbers — the "certify the generated derivative"
//! gate.

use fs_qty::Dims;

/// Registry index of a law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LawId(pub(crate) usize);

/// A scalar pointwise map with its exact derivative.
pub trait PointwiseLaw {
    /// Display name.
    fn name(&self) -> &'static str;
    /// N(u) at one dof.
    fn value(&self, u: f64) -> f64;
    /// dN/du at one dof.
    fn derivative(&self, u: f64) -> f64;
    /// Output dimensions given input dimensions; `None` when the law's
    /// exponent scaling overflows the i8 domain (bead sj31i.11 — a
    /// clamped exponent could alias false physics, so the type checker
    /// refuses instead).
    fn out_dims(&self, in_dims: Dims) -> Option<Dims>;
}

/// The reference nonlinearity: N(u) = α·u³ (a reaction term), with
/// dN/du = 3α·u². Cubic-degree dims scaling.
pub struct CubicReaction {
    /// Coefficient α.
    pub alpha: f64,
}

impl PointwiseLaw for CubicReaction {
    fn name(&self) -> &'static str {
        "cubic-reaction"
    }

    fn value(&self, u: f64) -> f64 {
        self.alpha * u * u * u
    }

    fn derivative(&self, u: f64) -> f64 {
        3.0 * self.alpha * u * u
    }

    fn out_dims(&self, in_dims: Dims) -> Option<Dims> {
        in_dims.checked_times(3)
    }
}
