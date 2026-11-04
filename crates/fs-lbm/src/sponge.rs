//! Absorbing layer (sponge) for the D2Q9 [`crate::core2::Grid`] —
//! the outlet-reflection treatment bead
//! frankensim-fsim-aeroacoustic-sources-9ok02 requires BEFORE any
//! aeroacoustic source spectrum can be trusted (the jet pilot's 6%
//! flux imbalance at Re 200 was outlet-reflection sensitivity).
//!
//! Mechanism: after each stream/collide step, cells inside the layer
//! blend toward the equilibrium of an authored far-field target
//! state, `h <- h + sigma(x) (f_eq(rho_t, u_t) - h)`, where
//! `h = f + S(u, F)/2` is the physical population in the Guo force
//! convention. Restoring the endpoint half-source retains the grid's
//! acceleration and external force without a spurious half-force impulse or
//! stress. This blends mass, physical momentum and physical second moments
//! with the same strength, including density-changing targets. With no force,
//! `h = f` and the original arithmetic is retained. The
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

use crate::core2::{Cell, Grid, raw_population_moments};
use crate::equilibrium;
use crate::thermal::forcing::reframe2;

/// Which x-side of the domain the layer occupies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpongeSide {
    /// Layer at low x (inlet side), strongest at x = 0.
    LeftX,
    /// Layer at high x (outlet side), strongest at x = nx - 1.
    RightX,
}

/// An absorbing layer over a contiguous x-range of the grid. The
/// target may be uniform ([`Sponge2::new`]) or a per-row PROFILE
/// ([`Sponge2::with_profile`]) — the latter turns the layer into a
/// FRINGE region (spectral-DNS style): in a periodic-x domain it
/// re-conditions outflow toward an authored inflow profile, so the
/// wrap delivers fresh inflow without open boundaries.
#[derive(Debug, Clone)]
pub struct Sponge2 {
    side: SpongeSide,
    width: usize,
    sigma_max: f64,
    /// One equilibrium target per y-row (length 1 = uniform).
    targets: Vec<[f64; crate::Q]>,
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
            targets: vec![equilibrium(target_rho, target_u[0], target_u[1])],
        }
    }

    /// Construct a fringe layer whose target varies per y-row:
    /// `profile[y] = (rho, [ux, uy])`. Same ramp and checks as
    /// [`Sponge2::new`], applied to every row of the profile.
    ///
    /// # Panics
    /// As [`Sponge2::new`], plus an empty profile.
    #[must_use]
    pub fn with_profile(
        side: SpongeSide,
        width: usize,
        sigma_max: f64,
        profile: &[(f64, [f64; 2])],
    ) -> Self {
        assert!(width > 0, "sponge width must be positive");
        assert!(
            sigma_max.is_finite() && sigma_max > 0.0 && sigma_max <= 1.0,
            "sponge sigma_max must lie in (0, 1]"
        );
        assert!(!profile.is_empty(), "sponge profile must be non-empty");
        let targets = profile
            .iter()
            .map(|&(rho, u)| {
                assert!(
                    rho.is_finite() && rho > 0.0,
                    "sponge target density must be positive and finite"
                );
                assert!(
                    u.iter().all(|c| c.is_finite()),
                    "sponge target velocity must be finite"
                );
                assert!(
                    u[0] * u[0] + u[1] * u[1] < 0.03,
                    "sponge target velocity exceeds the low-Mach boundary envelope"
                );
                equilibrium(rho, u[0], u[1])
            })
            .collect();
        Sponge2 {
            side,
            width,
            sigma_max,
            targets,
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
    /// The authored target velocity is physical, including when gravity or a
    /// spatial external force is present. The layer exchanges mass and momentum
    /// with its imposed far field; it is not a closed-system conservative step.
    ///
    /// # Panics
    /// If the layer is wider than the grid, or a per-row profile's
    /// length does not match `grid.ny`.
    pub fn apply(&self, grid: &mut Grid) {
        assert!(self.width <= grid.nx, "sponge layer wider than the grid");
        assert!(
            self.targets.len() == 1 || self.targets.len() == grid.ny,
            "sponge profile length must equal the grid height"
        );
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
                let target = if self.targets.len() == 1 {
                    &self.targets[0]
                } else {
                    &self.targets[y]
                };
                if grid.g == [0.0; 2] && grid.fext[i] == [0.0; 2] {
                    // Keep the original unforced acoustic/fringe path exactly.
                    for (fq, tq) in grid.f[i].iter_mut().zip(target) {
                        *fq += s * (tq - *fq);
                    }
                } else {
                    let old = grid.moments(i);
                    let force =
                        std::array::from_fn(|a| grid.g[a].mul_add(old.rho, grid.fext[i][a]));
                    let populations = &mut grid.f[i];
                    reframe2(populations, old.u, force, [0.0; 2]);
                    for (fq, tq) in populations.iter_mut().zip(target) {
                        *fq += s * (tq - *fq);
                    }
                    // The target can change density, so both rho*g and the
                    // velocity in S must use the actual blended physical state.
                    // Blending raw force-shifted populations, even toward a
                    // shifted target, would not preserve the second moment.
                    let (rho, momentum) = raw_population_moments(populations);
                    let velocity = momentum.map(|p| p / rho);
                    let force = std::array::from_fn(|a| grid.g[a].mul_add(rho, grid.fext[i][a]));
                    reframe2(populations, velocity, [0.0; 2], force);
                }
            }
        }
    }
}
