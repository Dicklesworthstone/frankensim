//! Non-Newtonian rheology (bead tfz.19): LOCAL relaxation-rate
//! adaptation. Each cell's shear rate comes from its own
//! non-equilibrium stress (no gradients, no neighbors), the apparent
//! viscosity from the constitutive law, and τ = 3ν + ½ — with a
//! DOCUMENTED stability floor: shear-thinning laws drive τ → ½ at
//! high shear, so τ is clamped to `TAU_FLOOR` and the clamp events
//! are counted (the scaling-assistant interplay the plan calls out:
//! a floored cell is a cell whose viscosity the lattice can no longer
//! represent at this resolution — the ledger says HOW MANY, honestly,
//! instead of silently going unstable).

use crate::core2::{Cell, Grid, shear_rate};
use crate::{CS2, Q};

/// The τ stability floor (ν > 0 with margin; BGK folklore-free zone
/// starts ≈ 0.505 at fixture scale).
pub const TAU_FLOOR: f64 = 0.505;

/// The τ ceiling: shear-thinning laws send ν → ∞ in plug regions
/// (γ̇ → 0), where BGK truncation and bounce-back slip errors grow
/// like (τ − ½)² — measured as a 12% global profile error with τ
/// reaching 20 before this cap existed. In a plug the stress is
/// ν·γ̇ ≈ 0 regardless of ν, so capping ν there leaves the momentum
/// balance intact; capped cells are counted in the ledger.
pub const TAU_CAP: f64 = 3.0;

fn require_positive_finite(name: &str, value: f64) {
    assert!(
        value.is_finite() && value > 0.0,
        "{name} must be positive and finite"
    );
}

fn require_nonnegative_finite(name: &str, value: f64) {
    assert!(
        value.is_finite() && value >= 0.0,
        "{name} must be nonnegative and finite"
    );
}

/// Constitutive law for the apparent kinematic viscosity ν(γ̇).
#[derive(Debug, Clone, Copy)]
pub enum Rheology {
    /// Fixed viscosity.
    Newtonian {
        /// Kinematic viscosity (lattice units).
        nu: f64,
    },
    /// Power law ν = k·γ̇^(n−1) (n < 1 shear-thinning).
    PowerLaw {
        /// Consistency index (lattice units).
        k: f64,
        /// Flow index.
        n: f64,
    },
    /// Carreau: ν = ν∞ + (ν₀ − ν∞)(1 + (λγ̇)²)^((n−1)/2).
    Carreau {
        /// Zero-shear viscosity.
        nu0: f64,
        /// Infinite-shear viscosity.
        nu_inf: f64,
        /// Relaxation time λ.
        lambda: f64,
        /// Flow index.
        n: f64,
    },
}

impl Rheology {
    /// Apparent kinematic viscosity at shear rate `gdot`.
    #[must_use]
    pub fn viscosity(&self, gdot: f64) -> f64 {
        require_nonnegative_finite("shear rate", gdot);
        match *self {
            Rheology::Newtonian { nu } => {
                require_positive_finite("Newtonian viscosity", nu);
                nu
            }
            Rheology::PowerLaw { k, n } => {
                require_positive_finite("power-law consistency k", k);
                require_positive_finite("power-law index n", n);
                // Guard γ̇ = 0 (unbounded for n < 1): the floor below
                // bounds τ anyway; use a tiny reference rate.
                let g = gdot.max(1e-12);
                k * g.powf(n - 1.0)
            }
            Rheology::Carreau {
                nu0,
                nu_inf,
                lambda,
                n,
            } => {
                require_positive_finite("Carreau zero-shear viscosity", nu0);
                require_positive_finite("Carreau infinite-shear viscosity", nu_inf);
                require_nonnegative_finite("Carreau relaxation time", lambda);
                require_positive_finite("Carreau flow index", n);
                let x = lambda * gdot;
                nu_inf + (nu0 - nu_inf) * x.mul_add(x, 1.0).powf((n - 1.0) / 2.0)
            }
        }
    }
}

/// Per-step rheology ledger.
#[derive(Debug, Clone, Copy, Default)]
pub struct RheologyStats {
    /// Cells whose τ hit the floor this update.
    pub floored: u32,
    /// Cells whose τ hit the ceiling this update (plug regions).
    pub capped: u32,
    /// Max shear rate seen.
    pub max_shear: f64,
    /// Min/max τ after update.
    pub tau_range: (f64, f64),
}

/// Update every fluid/interface cell's τ from its local shear rate
/// under `law`, flooring at [`TAU_FLOOR`]. Call between collide-stream
/// steps (one-step lag — the standard explicit local scheme).
#[must_use]
pub fn update_tau(grid: &mut Grid, law: Rheology) -> RheologyStats {
    let mut stats = RheologyStats {
        floored: 0,
        capped: 0,
        max_shear: 0.0,
        tau_range: (f64::INFINITY, 0.0),
    };
    for i in 0..grid.nx * grid.ny {
        if !matches!(grid.flags[i], Cell::Fluid | Cell::Interface) {
            continue;
        }
        let mm = grid.moments(i);
        let feq = crate::equilibrium(mm.rho, mm.u[0], mm.u[1]);
        let gdot = shear_rate(&grid.f[i], &feq, mm.rho, grid.tau[i]);
        let nu = law.viscosity(gdot);
        require_positive_finite("apparent viscosity", nu);
        let mut tau = nu / CS2 + 0.5;
        if tau < TAU_FLOOR {
            tau = TAU_FLOOR;
            stats.floored += 1;
        }
        if tau > TAU_CAP {
            tau = TAU_CAP;
            stats.capped += 1;
        }
        grid.tau[i] = tau;
        stats.max_shear = stats.max_shear.max(gdot);
        stats.tau_range.0 = stats.tau_range.0.min(tau);
        stats.tau_range.1 = stats.tau_range.1.max(tau);
    }
    if stats.tau_range.0 == f64::INFINITY {
        stats.tau_range = (0.0, 0.0);
    }
    stats
}

/// Analytic power-law Poiseuille x-velocity for a force-driven channel
/// with walls at y = −½ and y = ny − ½ (halfway bounce-back planes):
/// u(y) = (n/(n+1))·(gx/k)^(1/n)·(H^(1+1/n) − |y − y_c|^(1+1/n)) with
/// H the half-height. Peak at the centerline; n = 1 recovers the
/// parabola.
#[must_use]
pub fn powerlaw_poiseuille_analytic(gx: f64, k: f64, n: f64, ny: usize, y: usize) -> f64 {
    assert!(gx.is_finite(), "body force must be finite");
    require_positive_finite("power-law consistency k", k);
    require_positive_finite("power-law index n", n);
    assert!(ny > 0, "channel height must be positive");
    let h = ny as f64 / 2.0; // half-height between the halfway planes
    let yc = (ny as f64 - 1.0) / 2.0;
    let d = (y as f64 - yc).abs();
    let e = 1.0 + 1.0 / n;
    let forcing = gx / k;
    let signed_scale = forcing.signum() * forcing.abs().powf(1.0 / n);
    (n / (n + 1.0)) * signed_scale * (h.powf(e) - d.powf(e))
}

/// Convenience: run a force-driven periodic channel (walls top and
/// bottom) with local rheology updates until the centerline velocity
/// stabilizes; returns (grid, steps, final stats).
///
/// # Panics
/// Never at fixture scale (bounded loop).
#[must_use]
pub fn channel_flow(
    nx: usize,
    ny: usize,
    gx: f64,
    law: Rheology,
    max_steps: usize,
) -> (Grid, usize, RheologyStats) {
    assert!(nx > 0 && ny > 0, "channel dimensions must be positive");
    assert!(gx.is_finite(), "body force must be finite");
    assert!(
        ny < usize::MAX - 1,
        "channel height leaves no room for wall rows"
    );
    let grid_ny = ny + 2;
    let mut grid = Grid::uniform(nx, grid_ny, 1.0);
    grid.periodic_y = false;
    grid.g = [gx, 0.0];
    // Wall rows top and bottom.
    for x in 0..nx {
        let bottom = grid.idx(x, 0);
        grid.flags[bottom] = Cell::Wall;
        let top = grid.idx(x, grid_ny - 1);
        grid.flags[top] = Cell::Wall;
    }
    let mut scratch: Vec<[f64; Q]> = Vec::new();
    let mut stats = RheologyStats::default();
    let mut last_peak = 0.0f64;
    let mut steps = 0;
    for s in 0..max_steps {
        stats = update_tau(&mut grid, law);
        grid.step(&mut scratch);
        steps = s + 1;
        if s % 200 == 199 {
            let mid = grid.idx(0, usize::midpoint(0, grid_ny));
            let peak = grid.moments(mid).u[0];
            if (peak - last_peak).abs() < 1e-9 * peak.abs().max(1e-30) {
                break;
            }
            last_peak = peak;
        }
    }
    (grid, steps, stats)
}
