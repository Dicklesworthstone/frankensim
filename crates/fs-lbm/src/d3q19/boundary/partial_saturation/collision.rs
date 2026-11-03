//! Fluid relaxation removed by the original immersed solid collision.
//!
//! Apply relaxation to physical g = f + S/2, never to stored forced f. Both
//! parity modes conserve momentum there, so replacing their collision leaves
//! the original Guo source outside the solid blend and counts it exactly once.

use super::{CollisionModel3, OPP3, Q3};

#[derive(Clone, Copy)]
pub(super) struct ImmersedRelaxation3 {
    /// Even/shear relaxation controls the original saturation weight.
    pub(super) tau: f64,
    odd_rate: Option<f64>,
}

impl ImmersedRelaxation3 {
    /// Grid construction already validates rates. Nonlinear moment collisions
    /// still need their own immersed-boundary derivation and are not admitted.
    pub(super) fn new(model: CollisionModel3) -> Option<Self> {
        match model {
            CollisionModel3::Bgk { tau } => Some(Self {
                tau,
                odd_rate: None,
            }),
            CollisionModel3::Trt { tau, magic } => Some(Self {
                tau,
                odd_rate: Some(1.0 / (0.5 + magic / (tau - 0.5))),
            }),
            _ => None,
        }
    }

    /// One stack array per wet cell, shared by every body in that cell. The
    /// BGK branch retains its original operations. TRT evaluates the same
    /// nine opposite pairs as the native collision, with no extra collision,
    /// moment transform, allocation, or division inside the direction loop.
    pub(super) fn increments(self, physical: &[f64; Q3], equilibrium: &[f64; Q3]) -> [f64; Q3] {
        let Some(odd_rate) = self.odd_rate else {
            return core::array::from_fn(|q| (equilibrium[q] - physical[q]) / self.tau);
        };
        let even_rate = 1.0 / self.tau;
        let mut increments = [0.0; Q3];
        increments[0] = -even_rate * (physical[0] - equilibrium[0]);
        for q in (1..Q3).step_by(2) {
            let opposite = OPP3[q];
            let neq = physical[q] - equilibrium[q];
            let neq_opposite = physical[opposite] - equilibrium[opposite];
            let even = -even_rate * (0.5 * (neq + neq_opposite));
            let odd = -odd_rate * (0.5 * (neq - neq_opposite));
            increments[q] = even + odd;
            increments[opposite] = even - odd;
        }
        increments
    }
}
