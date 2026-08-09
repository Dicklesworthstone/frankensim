//! Absorbing layer (sponge) for the D2Q9 [`crate::core2::Grid`] —
//! the outlet-reflection treatment bead
//! frankensim-fsim-aeroacoustic-sources-9ok02 requires BEFORE any
//! aeroacoustic source spectrum can be trusted (the jet pilot's 6%
//! flux imbalance at Re 200 was outlet-reflection sensitivity).
//!
//! Mechanism: after each stream/collide step, cells inside the layer
//! blend toward the equilibrium of an authored far-field target
//! state, `f <- f + sigma(x) (f_eq(rho_t, u_t) - f)`, with a
//! QUADRATIC ramp `sigma(x) = sigma_max ((d + 1) / width)^2` growing
//! into the layer — a gentle impedance gradient, so the layer itself
//! reflects little (an abrupt sigma is its own reflector).
//!
//! The layer's ACTUAL reflection coefficient is measured, not
//! assumed: the conformance battery launches acoustic density pulses
//! (D2Q9 is weakly compressible; sound at `c_s = 1/sqrt(3)` lattice
//! units is exactly what re-enters the domain as spurious feedback)
//! and gates the reflected/incident amplitude ratio, against a
//! bounce-back wall CONTROL that catches the disabled-sponge
//! mutation.

use crate::core2::{Cell, Grid};
use crate::equilibrium;

/// Which x-side of the domain the layer occupies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpongeSide {
    /// Layer at low x (inlet side), strongest at x = 0.
    LeftX,
    /// Layer at high x (outlet side), strongest at x = nx - 1.
    RightX,
}

/// An absorbing layer over a contiguous x-range of the grid.
#[derive(Debug, Clone)]
pub struct Sponge2 {
    side: SpongeSide,
    width: usize,
    sigma_max: f64,
    target: [f64; crate::Q],
}

impl Sponge2 {
    /// Construct a checked layer: `width` cells on `side`, blending
    /// strength ramping quadratically from ~0 at the inner edge to
    /// `sigma_max` at the boundary, toward the equilibrium of
    /// `(target_rho, target_u)`.
    ///
    /// # Panics
    /// On non-finite or non-positive targets, `sigma_max` outside
    /// (0, 1], zero width, or a target speed outside the low-Mach
    /// envelope (crate boundary convention: construction-time
    /// asserts, matching `VelocityPressureX2`).
    #[must_use]
    pub fn new(
        side: SpongeSide,
        width: usize,
        sigma_max: f64,
        target_rho: f64,
        target_u: [f64; 2],
    ) -> Self {
        assert!(width > 0, "sponge width must be positive");
        assert!(
            sigma_max.is_finite() && sigma_max > 0.0 && sigma_max <= 1.0,
            "sponge sigma_max must lie in (0, 1]"
        );
        assert!(
            target_rho.is_finite() && target_rho > 0.0,
            "sponge target density must be positive and finite"
        );
        assert!(
            target_u.iter().all(|c| c.is_finite()),
            "sponge target velocity must be finite"
        );
        let speed_sq = target_u[0] * target_u[0] + target_u[1] * target_u[1];
        assert!(
            speed_sq < 0.03,
            "sponge target velocity exceeds the low-Mach boundary envelope"
        );
        Sponge2 {
            side,
            width,
            sigma_max,
            target: equilibrium(target_rho, target_u[0], target_u[1]),
        }
    }

    /// Blend strength at depth `d` into the layer (0 = inner edge,
    /// `width - 1` = domain boundary).
    #[must_use]
    fn sigma(&self, d: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let t = (d as f64 + 1.0) / self.width as f64;
        self.sigma_max * t * t
    }

    /// Apply one blending pass (call once per LBM step, after
    /// [`Grid::step`]). Only `Cell::Fluid` cells are touched.
    ///
    /// # Panics
    /// If the layer is wider than the grid.
    pub fn apply(&self, grid: &mut Grid) {
        assert!(self.width <= grid.nx, "sponge layer wider than the grid");
        for d in 0..self.width {
            let x = match self.side {
                SpongeSide::RightX => grid.nx - self.width + d,
                SpongeSide::LeftX => self.width - 1 - d,
            };
            let s = self.sigma(d);
            for y in 0..grid.ny {
                let i = grid.idx(x, y);
                if grid.flags[i] != Cell::Fluid {
                    continue;
                }
                for (fq, tq) in grid.f[i].iter_mut().zip(&self.target) {
                    *fq += s * (tq - *fq);
                }
            }
        }
    }
}
