//! Non-Newtonian rheology (bead tfz.19): LOCAL relaxation-rate
//! adaptation. Each cell's shear rate comes from its own
//! non-equilibrium stress (no gradients, no neighbors), the apparent
//! viscosity from the constitutive law, and τ = 3ν + ½. The local
//! stress and viscosity are solved together, rather than interpreting
//! the current populations through the previous relaxation time. The
//! documented stability floor handles shear-thinning laws driving τ → ½ at
//! high shear, so τ is clamped to `TAU_FLOOR` and the clamp events
//! are counted (the scaling-assistant interplay the plan calls out:
//! a floored cell is a cell whose viscosity the lattice can no longer
//! represent at this resolution — the ledger says HOW MANY, honestly,
//! instead of silently going unstable).

use crate::core2::{
    Cell, CollisionModel2, Grid, force_corrected_shear_rate, force_corrected_stress,
};
use crate::{CS2, Q};

/// The τ stability floor (ν > 0 with margin; BGK folklore-free zone
/// starts ≈ 0.505 at fixture scale).
pub const TAU_FLOOR: f64 = 0.505;

/// The τ ceiling: shear-thinning laws send ν → ∞ in plug regions
/// (γ̇ → 0), where BGK truncation and bounce-back slip errors grow
/// like (τ − ½)² — measured as a 12% global profile error with τ
/// reaching 20 before this cap existed. This is a numerical viscosity
/// limit, not an equivalent constitutive law: capped cells are counted,
/// and physical rheology fidelity requires rescaling when it is active.
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
    /// Regularized Herschel–Bulkley yield-stress fluid:
    /// `nu = k max(gamma, gamma_min)^(n-1)
    ///       + Y (1 - exp(-m gamma)) / gamma`.
    ///
    /// `n=1` is the Bingham–Papanastasiou law. The yield term has the exact
    /// zero-rate limit `Y*m`; the explicit power-law plateau also prevents
    /// divergence for `n<1`. Finite regularization admits creep below yield,
    /// so this is not an exact rigid-plug or cessation model.
    ///
    /// Papanastasiou, J. Rheol. 31 (1987), DOI:10.1122/1.549926;
    /// Jeon & Hodges (2018), equation (6),
    /// <https://doi.org/10.5194/nhess-18-303-2018>.
    HerschelBulkley {
        /// Kinematic consistency index in lattice units; positive.
        k: f64,
        /// Positive flow index; below one thins, above one thickens.
        n: f64,
        /// Yield stress divided by density, in lattice velocity squared.
        /// For SI input use `(stress_pa / rho_kg_m3) * dt_s^2 / dx_m^2`.
        yield_stress_over_density: f64,
        /// Papanastasiou regularization time `m` in lattice time units.
        /// SI time converts as `m_si / dt_s`; zero disables the yield term.
        regularization_time: f64,
        /// Positive low-shear cutoff `gamma_min` for the power-law part,
        /// in inverse lattice time. This is a declared constitutive
        /// regularization, separate from the lattice tau floor/cap.
        minimum_shear_rate: f64,
    },
}

impl Rheology {
    /// Apparent kinematic viscosity at shear rate `gdot`.
    #[must_use]
    pub fn viscosity(&self, gdot: f64) -> f64 {
        require_nonnegative_finite("shear rate", gdot);
        self.validate();
        self.viscosity_and_slope(gdot).0
    }

    fn validate(&self) {
        match *self {
            Rheology::Newtonian { nu } => {
                require_positive_finite("Newtonian viscosity", nu);
            }
            Rheology::PowerLaw { k, n } => {
                require_positive_finite("power-law consistency k", k);
                require_positive_finite("power-law index n", n);
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
                assert!(
                    n <= 1.0 || nu0 >= nu_inf || lambda == 0.0,
                    "thickening Carreau law requires nu0 >= nu_inf to stay positive"
                );
            }
            Rheology::HerschelBulkley {
                k,
                n,
                yield_stress_over_density,
                regularization_time,
                minimum_shear_rate,
            } => {
                require_positive_finite("Herschel-Bulkley consistency k", k);
                require_positive_finite("Herschel-Bulkley index n", n);
                require_nonnegative_finite("kinematic yield stress", yield_stress_over_density);
                require_nonnegative_finite("yield regularization time", regularization_time);
                require_positive_finite("power-law minimum shear rate", minimum_shear_rate);
            }
        }
    }

    // Return (nu, gamma*dnu/dgamma); the scaled derivative avoids division
    // by gamma at zero strain. Parameters are validated once per grid update.
    fn viscosity_and_slope(&self, gdot: f64) -> (f64, f64) {
        match *self {
            Rheology::Newtonian { nu } => (nu, 0.0),
            Rheology::PowerLaw { k, n } => {
                // Preserve the explicitly regularized power law below 1e-12
                // in lattice shear units; it is constant in that small band.
                let g = gdot.max(1e-12);
                let nu = k * fs_math::det::pow(g, n - 1.0);
                (nu, if gdot > 1e-12 { (n - 1.0) * nu } else { 0.0 })
            }
            Rheology::Carreau {
                nu0,
                nu_inf,
                lambda,
                n,
            } => {
                if n == 1.0 || lambda == 0.0 || nu0 == nu_inf {
                    return (nu0, 0.0);
                }
                let x = lambda * gdot;
                let x2 = x * x;
                let shape = fs_math::det::pow(x.mul_add(x, 1.0), (n - 1.0) / 2.0);
                let fraction = if x2.is_finite() { x2 / (1.0 + x2) } else { 1.0 };
                let nu = if n < 1.0 {
                    // Positive interpolation avoids cancellation between two
                    // nearly equal large viscosities near zero strain.
                    nu0 * shape + nu_inf * (1.0 - shape)
                } else {
                    nu_inf + (nu0 - nu_inf) * shape
                };
                (nu, (nu0 - nu_inf) * shape * (n - 1.0) * fraction)
            }
            Rheology::HerschelBulkley {
                k,
                n,
                yield_stress_over_density,
                regularization_time,
                minimum_shear_rate,
            } => {
                let base = if n == 1.0 {
                    k
                } else {
                    k * fs_math::det::pow(gdot.max(minimum_shear_rate), n - 1.0)
                };
                let base_slope = if gdot > minimum_shear_rate {
                    (n - 1.0) * base
                } else {
                    0.0
                };
                let (yield_viscosity, yield_slope) =
                    regularized_yield(gdot, yield_stress_over_density, regularization_time);
                (base + yield_viscosity, base_slope + yield_slope)
            }
        }
    }

    pub(crate) fn prepare(self) -> PreparedRheology {
        self.validate();
        let constant = match self {
            Self::Newtonian { nu } => Some(nu),
            Self::PowerLaw { k, n } if n == 1.0 => Some(k),
            Self::Carreau {
                nu0,
                nu_inf,
                lambda,
                n,
            } if lambda == 0.0 || n == 1.0 || nu0 == nu_inf => Some(nu0),
            Self::HerschelBulkley {
                k,
                n,
                yield_stress_over_density,
                regularization_time,
                ..
            } if n == 1.0 && (yield_stress_over_density == 0.0 || regularization_time == 0.0) => {
                Some(k)
            }
            _ => None,
        };
        PreparedRheology {
            law: self,
            constant_tau: constant.map(|nu| 0.5 + nu / CS2),
        }
    }
}

// Yield contribution and gamma*dnu/dgamma. expm1 preserves the zero-rate
// plateau; the small-x series also keeps the derivative from subtracting
// nearly equal numbers. No platform libm, fit or rate-dependent iteration.
fn regularized_yield(rate: f64, yield_stress: f64, time: f64) -> (f64, f64) {
    if yield_stress == 0.0 || time == 0.0 {
        return (0.0, 0.0);
    }
    let x = time * rate;
    if x == 0.0 {
        return (yield_stress * time, 0.0);
    }
    if x < 1e-4 {
        // phi=(1-exp(-x))/x; x*phi'=-x/2+x²/3-x³/8+x⁴/30+O(x⁵).
        let zero_rate = yield_stress * time;
        let phi = 1.0 + x * (-0.5 + x * (1.0 / 6.0 + x * (-1.0 / 24.0 + x / 120.0)));
        let slope = x * (-0.5 + x * (1.0 / 3.0 + x * (-1.0 / 8.0 + x / 30.0)));
        return if zero_rate.is_finite() {
            (zero_rate * phi, zero_rate * slope)
        } else {
            // The declared zero-rate plateau may overflow even when its
            // scaled value at this nonzero rate remains representable.
            ((yield_stress * phi) * time, (yield_stress * slope) * time)
        };
    }
    // One strict exponential evaluation per query. Near zero expm1 retains
    // the onset accurately; outside its own small-argument branch it would
    // call exp again, so share that exp with the analytic derivative instead.
    let (onset, decay) = if x < 0.5 * core::f64::consts::LN_2 {
        let onset = -fs_math::det::expm1(-x);
        (onset, 1.0 - onset)
    } else {
        let decay = fs_math::det::exp(-x);
        (1.0 - decay, decay)
    };
    let viscosity = if x <= 1.0 {
        let plateau = yield_stress * time;
        let factor = onset / x;
        if plateau.is_finite() {
            plateau * factor
        } else {
            (yield_stress * factor) * time
        }
    } else {
        let unscaled = yield_stress / rate;
        if unscaled.is_finite() {
            unscaled * onset
        } else {
            (yield_stress * onset) / rate
        }
    };
    let ratio = if x.is_finite() {
        x * decay / onset
    } else {
        0.0
    };
    (viscosity, viscosity * (ratio - 1.0))
}

/// One validated law shared by the two lattice dimensions. Constant laws
/// evaluate their constitutive relation once, without per-cell iteration.
pub(crate) struct PreparedRheology {
    law: Rheology,
    constant_tau: Option<f64>,
}

pub(crate) struct LocalRelaxation {
    pub tau: f64,
    pub shear_rate: f64,
    pub floored: bool,
    pub capped: bool,
}

impl LocalRelaxation {
    fn from_target(target: f64, numerator: f64, fixed_rate: f64) -> Self {
        let tau = target.clamp(TAU_FLOOR, TAU_CAP);
        let shear_rate = combined_strain_rate(numerator / tau, fixed_rate);
        require_nonnegative_finite("applied shear rate", shear_rate);
        Self {
            tau,
            shear_rate,
            floored: target < TAU_FLOOR,
            capped: target > TAU_CAP,
        }
    }
}

// Preserve the original BGK arithmetic when every second moment has the same
// relaxation time. The fixed part is the strain from independently relaxed
// trace moments, orthogonal to the deviatoric part in the Frobenius norm.
fn combined_strain_rate(variable_rate: f64, fixed_rate: f64) -> f64 {
    if fixed_rate == 0.0 {
        variable_rate
    } else {
        // The fixed IEEE operation tree avoids platform-libm hypot while
        // scaling before squaring keeps small/large admitted rates usable.
        let scale = variable_rate.max(fixed_rate);
        let a = variable_rate / scale;
        let b = fixed_rate / scale;
        scale * a.mul_add(a, b * b).sqrt()
    }
}

impl PreparedRheology {
    /// Solve tau = 1/2 + nu(A/tau)/cs², with A=gamma*tau from the
    /// force-corrected population stress. For the admitted laws, gamma*nu
    /// is increasing, giving a unique positive root. The previous tau is
    /// only a warm start; it never determines the returned constitutive rate.
    pub(crate) fn solve(&self, numerator: f64, previous_tau: f64) -> LocalRelaxation {
        self.solve_with_fixed_strain(numerator, 0.0, previous_tau)
    }

    /// Independently relaxed trace modes contribute a fixed orthogonal strain
    /// invariant B, so gamma = sqrt((A/tau)^2 + B^2). Only the deviatoric
    /// contribution varies with shear tau. Treating B as A/tau would apply
    /// the shear relaxation rate to the collision's independent bulk mode.
    fn solve_with_fixed_strain(
        &self,
        numerator: f64,
        fixed_rate: f64,
        previous_tau: f64,
    ) -> LocalRelaxation {
        require_nonnegative_finite("non-equilibrium strain numerator", numerator);
        require_nonnegative_finite("independent strain rate", fixed_rate);
        if let Some(tau) = self.constant_tau {
            return LocalRelaxation::from_target(tau, numerator, fixed_rate);
        }
        if numerator == 0.0 {
            return LocalRelaxation::from_target(
                0.5 + self.law.viscosity_and_slope(fixed_rate).0 / CS2,
                0.0,
                fixed_rate,
            );
        }
        let target = |tau: f64| {
            let variable_rate = numerator / tau;
            let rate = combined_strain_rate(variable_rate, fixed_rate);
            let (nu, slope) = self.law.viscosity_and_slope(rate);
            // Zero/overflow may result from an extreme power law. Their
            // positive limiting responses still classify a floor/cap, while
            // NaN or a negative apparent viscosity cannot be published.
            assert!(!nu.is_nan() && nu >= 0.0, "invalid apparent viscosity");
            let fraction = if fixed_rate == 0.0 {
                1.0
            } else {
                (variable_rate / rate).powi(2)
            };
            (0.5 + nu / CS2, 1.0 + slope * fraction / (CS2 * tau))
        };
        let (mut lo, mut hi) = (TAU_FLOOR, TAU_CAP);
        let mut tau = if previous_tau.is_finite() {
            previous_tau.clamp(lo, hi)
        } else {
            0.5 * (lo + hi)
        };
        let (mut next, mut derivative) = target(tau);
        if (tau == lo && next <= lo) || (tau == hi && next >= hi) {
            return LocalRelaxation::from_target(next, numerator, fixed_rate);
        }
        let initial_residual = tau - next;
        if initial_residual.abs() <= 1e-12 * tau {
            return LocalRelaxation::from_target(tau, numerator, fixed_rate);
        }
        // A settled cell needs only the warm-start evaluation above. If it
        // changes, only the endpoint in the root's direction is evaluated.
        if initial_residual < 0.0 {
            let cap_target = if tau == hi { next } else { target(hi).0 };
            if cap_target >= hi {
                return LocalRelaxation::from_target(cap_target, numerator, fixed_rate);
            }
            lo = tau;
        } else {
            let floor_target = if tau == lo { next } else { target(lo).0 };
            if floor_target <= lo {
                return LocalRelaxation::from_target(floor_target, numerator, fixed_rate);
            }
            hi = tau;
        }
        let mut previous_step = hi - lo;
        for _ in 0..96 {
            let residual = tau - next;
            if residual.abs() <= 1e-12 * tau {
                return LocalRelaxation::from_target(tau, numerator, fixed_rate);
            }
            if residual < 0.0 {
                lo = tau;
            } else {
                hi = tau;
            }
            let correction = residual / derivative;
            let newton = tau - correction;
            let candidate = if derivative > 0.0
                && newton > lo
                && newton < hi
                && correction.abs() <= 0.5 * previous_step
            {
                newton
            } else {
                0.5 * (lo + hi)
            };
            previous_step = (candidate - tau).abs();
            tau = candidate;
            (next, derivative) = target(tau);
        }
        panic!("local rheology did not converge within its bounded solve");
    }
}

/// Per-step rheology ledger.
#[derive(Debug, Clone, Copy, Default)]
pub struct RheologyStats {
    /// Cells whose τ hit the floor this update.
    pub floored: u32,
    /// Cells whose τ hit the ceiling this update (plug regions).
    pub capped: u32,
    /// Max shear rate consistent with the applied relaxation time.
    pub max_shear: f64,
    /// Min/max τ after update.
    pub tau_range: (f64, f64),
}

/// Update every fluid/interface cell's τ from its local shear rate
/// under `law`, flooring at [`TAU_FLOOR`]. Call between collide-stream
/// steps. The local stress/viscosity relation is solved self-consistently;
/// this does not make the time integration implicit. The strain
/// includes the Guo half-source correction for both body acceleration and
/// local external force; uniform acceleration is not material shear.
/// Central-moment collision retains its unit bulk relaxation time: its trace
/// strain and the shear-tau-dependent deviatoric strain are resolved separately.
#[must_use]
pub fn update_tau(grid: &mut Grid, law: Rheology) -> RheologyStats {
    let law = law.prepare();
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
        let force = std::array::from_fn(|d| grid.g[d].mul_add(mm.rho, grid.fext[i][d]));
        let local = if matches!(grid.collision, CollisionModel2::CentralMoment) {
            let [sxx, sxy, syy] = force_corrected_stress(&grid.f[i], mm.rho, mm.u, force);
            // Pi_dev = -2 rho cs² tau D_dev and tr(Pi) =
            // -2 rho cs² tau_bulk tr(D), with tau_bulk=1 in this model.
            // In 2-D, 2 D:D = 4(D_dev,xx² + D_xy²) + tr(D)².
            let scale = 3.0 / mm.rho;
            let numerator = scale * combined_strain_rate((0.5 * (sxx - syy)).abs(), sxy.abs());
            let dilation = 0.5 * scale * (sxx + syy).abs();
            law.solve_with_fixed_strain(numerator, dilation, grid.tau[i])
        } else {
            let numerator = force_corrected_shear_rate(&grid.f[i], mm.rho, mm.u, force, 1.0);
            law.solve(numerator, grid.tau[i])
        };
        grid.tau[i] = local.tau;
        stats.floored += u32::from(local.floored);
        stats.capped += u32::from(local.capped);
        stats.max_shear = stats.max_shear.max(local.shear_rate);
        stats.tau_range.0 = stats.tau_range.0.min(local.tau);
        stats.tau_range.1 = stats.tau_range.1.max(local.tau);
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
    let signed_scale = forcing.signum() * fs_math::det::pow(forcing.abs(), 1.0 / n);
    (n / (n + 1.0)) * signed_scale * (fs_math::det::pow(h, e) - fs_math::det::pow(d, e))
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

#[cfg(test)]
mod tests {
    use super::{
        Rheology, TAU_CAP, TAU_FLOOR, channel_flow, powerlaw_poiseuille_analytic, update_tau,
    };
    use crate::{CS2, E, Grid, W, core2::CollisionModel2, equilibrium};

    /// G0/G3: a mixed volumetric/deviatoric deformation has the same physical
    /// invariant for BGK and central moments, although their trace relaxation
    /// times differ. Manufacture each stress tensor from its own collision law
    /// and independently reconstruct its populations, including the Guo frame.
    #[test]
    fn central_moment_rheology_uses_its_independent_bulk_relaxation() {
        for law in [
            Rheology::Newtonian { nu: 0.02 },
            Rheology::Carreau {
                nu0: 0.2,
                nu_inf: 0.01,
                lambda: 1000.0,
                n: 0.5,
            },
            Rheology::PowerLaw { k: 0.5, n: 1.5 },
        ] {
            for strain in [
                [[0.003, 0.0], [0.0, 0.003]],
                [[0.003, -0.004], [-0.004, -0.001]],
                [[-0.004, 0.002], [0.002, -0.001]],
            ] {
                let rate = (2.0_f64 * strain.iter().flatten().map(|s| s * s).sum::<f64>()).sqrt();
                let expected_tau = 0.5 + law.viscosity(rate) / CS2;
                assert!(expected_tau > TAU_FLOOR && expected_tau < TAU_CAP);
                for collision in [CollisionModel2::Bgk, CollisionModel2::CentralMoment] {
                    let bulk_tau = if collision == CollisionModel2::CentralMoment {
                        1.0
                    } else {
                        expected_tau
                    };
                    for rho in [0.4, 2.5] {
                        for force in [[0.0; 2], [0.003, -0.002]] {
                            let u = [0.04, -0.03];
                            let isotropic = 0.5 * (strain[0][0] + strain[1][1]);
                            let stress: [[f64; 2]; 2] = std::array::from_fn(|a| {
                                std::array::from_fn(|b| {
                                    let trace = if a == b { isotropic } else { 0.0 };
                                    -2.0 * rho
                                        * CS2
                                        * (expected_tau * (strain[a][b] - trace) + bulk_tau * trace)
                                })
                            });
                            let mut populations = equilibrium(rho, u[0], u[1]);
                            for (q, direction) in E.iter().enumerate() {
                                let e = [f64::from(direction.0), f64::from(direction.1)];
                                let mut contraction = 0.0;
                                for a in 0..2 {
                                    for b in 0..2 {
                                        let delta = if a == b { CS2 } else { 0.0 };
                                        contraction += (e[a] * e[b] - delta) * stress[a][b];
                                    }
                                }
                                let eu = e[0] * u[0] + e[1] * u[1];
                                let source: f64 = (0..2)
                                    .map(|a| (3.0 * (e[a] - u[a]) + 9.0 * eu * e[a]) * force[a])
                                    .sum();
                                populations[q] +=
                                    W[q] * (contraction / (2.0 * CS2 * CS2) - 0.5 * source);
                            }
                            for old_tau in [0.51, 1.0, 2.9] {
                                let mut grid = Grid::uniform(2, 2, old_tau);
                                grid.collision = collision;
                                grid.f.fill(populations);
                                grid.fext.fill(force);
                                let stats = update_tau(&mut grid, law);
                                assert!(
                                    (stats.max_shear - rate).abs() < 5e-14,
                                    "{collision:?}, {law:?}, strain={strain:?}, rho={rho}, force={force:?}: {} != {rate}",
                                    stats.max_shear,
                                );
                                assert!(grid.tau.iter().all(|t| (t - expected_tau).abs() < 4e-12));
                                assert_eq!((stats.floored, stats.capped), (0, 0));
                            }
                        }
                    }
                }
            }
        }
    }

    /// G0/G3: changing the Guo force frame must not change the physical
    /// strain or the apparent viscosity. The equilibrium case represents
    /// uniform accelerated flow, whose material strain is exactly zero.
    #[test]
    fn forced_rheology_recovers_physical_strain_and_viscosity() {
        let law = Rheology::Carreau {
            nu0: 0.2,
            nu_inf: 0.01,
            lambda: 1000.0,
            n: 0.5,
        };
        let u = [0.06, -0.03];
        for rho in [0.4, 2.5] {
            for tau in [0.7, 1.8] {
                for strain in [
                    [[0.0, 0.0], [0.0, 0.0]],
                    [[0.003, -0.004], [-0.004, -0.003]],
                ] {
                    let rate =
                        (2.0_f64 * strain.iter().flatten().map(|s| s * s).sum::<f64>()).sqrt();
                    let expected_tau = 0.5 + law.viscosity(rate) / CS2;
                    for (acceleration, local_force) in [
                        ([0.0, 0.0], [0.0, 0.0]),
                        ([0.004, -0.002], [0.003, 0.005]),
                        ([-0.003, 0.005], [-0.002, -0.004]),
                    ] {
                        let mut grid = Grid::uniform(2, 2, tau);
                        grid.g = acceleration;
                        grid.fext.fill(local_force);
                        let force: [f64; 2] =
                            std::array::from_fn(|d| acceleration[d].mul_add(rho, local_force[d]));
                        let mut populations = equilibrium(rho, u[0], u[1]);
                        for q in 0..E.len() {
                            let e = [f64::from(E[q].0), f64::from(E[q].1)];
                            let eu = e[0] * u[0] + e[1] * u[1];
                            let trace = strain[0][0] + strain[1][1];
                            let quadratic = e[0] * (strain[0][0] * e[0] + strain[0][1] * e[1])
                                + e[1] * (strain[1][0] * e[0] + strain[1][1] * e[1]);
                            // Hermite reconstruction of Π = -2 rho cs² tau D,
                            // followed by the independent population-space
                            // half-source transform f_bar = f - S/2.
                            let source: f64 = (0..2)
                                .map(|d| (3.0 * (e[d] - u[d]) + 9.0 * eu * e[d]) * force[d])
                                .sum();
                            populations[q] -=
                                3.0 * W[q] * rho * expected_tau * (quadratic - CS2 * trace)
                                    + 0.5 * W[q] * source;
                        }
                        grid.f.fill(populations);
                        let stats = update_tau(&mut grid, law);
                        assert!(
                            (stats.max_shear - rate).abs() < 2e-14,
                            "rho={rho}, tau={tau}, force={force:?}: {} != {rate}",
                            stats.max_shear,
                        );
                        assert!(grid.tau.iter().all(|t| (t - expected_tau).abs() < 2e-12));
                        assert_eq!((stats.floored, stats.capped), (0, 0));
                        let applied = grid.tau.clone();
                        let repeated = update_tau(&mut grid, law);
                        assert_eq!(
                            grid.tau, applied,
                            "unchanged populations must not relax the law again"
                        );
                        assert_eq!(repeated.max_shear, stats.max_shear);
                    }
                }
            }
        }
    }

    /// Invert manufactured stresses over both thinning and thickening
    /// regimes. Large flow indices make undamped fixed-point iteration
    /// diverge; the safeguarded constitutive root remains well defined.
    #[test]
    fn local_rheology_is_self_consistent_for_thinning_and_thickening() {
        let cases = [
            (Rheology::Newtonian { nu: 0.1 }, [0.01, 0.1, 1.0]),
            (Rheology::PowerLaw { k: 0.01, n: 0.2 }, [0.01, 0.1, 1.0]),
            (Rheology::PowerLaw { k: 0.3, n: 5.0 }, [0.4, 0.8, 1.1]),
            (
                Rheology::Carreau {
                    nu0: 0.3,
                    nu_inf: 0.01,
                    lambda: 50.0,
                    n: 0.3,
                },
                [0.01, 0.1, 1.0],
            ),
            (
                Rheology::Carreau {
                    nu0: 0.01,
                    nu_inf: 0.2,
                    lambda: 10.0,
                    n: 0.5,
                },
                [0.01, 0.1, 1.0],
            ),
            (
                Rheology::Carreau {
                    nu0: 0.1,
                    nu_inf: 0.01,
                    lambda: 2.0,
                    n: 3.0,
                },
                [0.01, 0.1, 1.0],
            ),
            (
                Rheology::HerschelBulkley {
                    k: 0.03,
                    n: 0.5,
                    yield_stress_over_density: 0.003,
                    regularization_time: 10.0,
                    minimum_shear_rate: 0.001,
                },
                [0.01, 0.1, 1.0],
            ),
            (
                Rheology::HerschelBulkley {
                    k: 0.1,
                    n: 2.0,
                    yield_stress_over_density: 0.02,
                    regularization_time: 20.0,
                    minimum_shear_rate: 1e-6,
                },
                [0.01, 0.1, 1.0],
            ),
        ];
        for (law, rates) in cases {
            let prepared = law.prepare();
            for rate in rates {
                let expected = 0.5 + law.viscosity(rate) / CS2;
                assert!(expected > TAU_FLOOR && expected < TAU_CAP);
                for (numerator, fixed_rate) in [
                    (rate * expected, 0.0),
                    (0.6 * rate * expected, 0.8 * rate),
                    (0.0, rate),
                ] {
                    for previous in [TAU_FLOOR, 0.8, 2.9, TAU_CAP, f64::NAN] {
                        let solved =
                            prepared.solve_with_fixed_strain(numerator, fixed_rate, previous);
                        assert!(
                            (solved.tau - expected).abs() < 4e-12,
                            "{law:?}, rate={rate}, fixed={fixed_rate}, previous={previous}"
                        );
                        let target = 0.5 + law.viscosity(solved.shear_rate) / CS2;
                        assert!((solved.tau - target).abs() <= 1.1e-12 * solved.tau);
                        assert!(!solved.floored && !solved.capped);
                    }
                }
            }
        }
    }

    #[test]
    fn local_rheology_clipping_uses_the_applied_strain_rate() {
        for law in [
            Rheology::Newtonian { nu: 1e-5 },
            Rheology::PowerLaw { k: 1e-8, n: 0.5 },
            Rheology::Newtonian { nu: 10.0 },
            Rheology::PowerLaw { k: 10.0, n: 0.5 },
        ] {
            let prepared = law.prepare();
            let numerator = 0.1;
            for previous in [TAU_FLOOR, 0.9, TAU_CAP] {
                let local = prepared.solve(numerator, previous);
                assert_eq!(local.shear_rate, numerator / local.tau);
                let target = 0.5 + law.viscosity(local.shear_rate) / CS2;
                assert_eq!(local.tau, target.clamp(TAU_FLOOR, TAU_CAP));
                assert_eq!(local.floored, target < TAU_FLOOR);
                assert_eq!(local.capped, target > TAU_CAP);
                assert!(local.floored || local.capped);
            }
        }
        // The inherited zero-shear power-law regularization is explicit.
        let zero = Rheology::PowerLaw { k: 0.1, n: 0.5 }
            .prepare()
            .solve(0.0, 0.8);
        assert_eq!(zero.tau, TAU_CAP);
        assert_eq!(zero.shear_rate, 0.0);
        assert!(zero.capped);
    }

    #[test]
    fn yield_regularization_has_declared_plateaus_and_limits() {
        let law = Rheology::HerschelBulkley {
            k: 0.05,
            n: 1.0,
            yield_stress_over_density: 0.00005,
            regularization_time: 15_000.0,
            minimum_shear_rate: 1e-12,
        };
        assert_eq!(law.viscosity(0.0), 0.05 + 0.00005 * 15_000.0);
        // Neither 0/0 nor 1-exp(-x) cancellation may erase the yield plateau.
        assert_eq!(law.viscosity(1e-25), law.viscosity(0.0));
        assert!((law.viscosity(1.0) - 0.05005).abs() < 1e-16);
        // A finite constitutive result must survive an overflowing
        // intermediate plateau or unregularized Y/gamma. The independent
        // reference multiplies the bounded onset before division.
        for (time, rate) in [(2.0, 0.25), (3.0, 0.5), (1.7977, 1e-5 / 1.7977)] {
            let extreme = Rheology::HerschelBulkley {
                k: 0.1,
                n: 1.0,
                yield_stress_over_density: 1e308,
                regularization_time: time,
                minimum_shear_rate: 1e-12,
            };
            let expected = (1e308 * (-(-time * rate).exp_m1())) / rate;
            let (actual, slope) = extreme.viscosity_and_slope(rate);
            assert!(actual.is_finite() && expected.is_finite());
            assert!((actual / expected - 1.0).abs() < 4e-15);
            assert!(slope.is_finite() && slope < 0.0);
        }
        let mut previous = 0.0;
        for rate in [0.0, 1e-25, 1e-10, 1e-6, 1e-4, 0.01, 1.0] {
            let stress = rate * law.viscosity(rate);
            assert!(stress >= previous, "yield stress response must be monotone");
            previous = stress;
        }
        for n in [0.5, 1.0, 1.5] {
            let yielded = Rheology::HerschelBulkley {
                k: 0.1,
                n,
                yield_stress_over_density: 0.0,
                regularization_time: 1000.0,
                minimum_shear_rate: 1e-12,
            };
            let power = Rheology::PowerLaw { k: 0.1, n };
            for rate in [0.0, 1e-25, 0.01, 1.0] {
                assert_eq!(yielded.viscosity(rate), power.viscosity(rate));
            }
        }
        // Declared regularization and numerical resolution limits are distinct:
        // a large zero-rate viscosity is reported as capped, not as exact yield.
        let capped = Rheology::HerschelBulkley {
            k: 0.1,
            n: 1.0,
            yield_stress_over_density: 0.01,
            regularization_time: 1000.0,
            minimum_shear_rate: 1e-8,
        }
        .prepare()
        .solve(0.0, 0.8);
        assert_eq!(capped.tau, TAU_CAP);
        assert!(capped.capped);
        assert!(!capped.floored);
        assert_eq!(capped.shear_rate, 0.0);
    }

    #[test]
    fn yield_law_rejects_undeclared_or_nonphysical_regularization() {
        for (yield_stress, time, minimum_rate) in [
            (-1.0, 1.0, 1e-8),
            (1.0, f64::INFINITY, 1e-8),
            (1.0, -1.0, 1e-8),
            (1.0, 1.0, 0.0),
            (1.0, 1.0, f64::NAN),
        ] {
            assert!(
                std::panic::catch_unwind(|| {
                    Rheology::HerschelBulkley {
                        k: 0.1,
                        n: 1.0,
                        yield_stress_over_density: yield_stress,
                        regularization_time: time,
                        minimum_shear_rate: minimum_rate,
                    }
                    .viscosity(0.1)
                })
                .is_err()
            );
        }
    }

    // Independent steady momentum balance: stress/rho = |g*y|. Invert the
    // scalar Bingham-Papanastasiou stress, then integrate gamma*d(stress)/g
    // analytically. This reference never invokes the production viscosity or
    // local kinetic constitutive solver, and uses platform exp as its oracle.
    fn bingham_velocity(
        y: f64,
        half_height: f64,
        g: f64,
        k: f64,
        yield_stress: f64,
        m: f64,
    ) -> f64 {
        let shear = |stress: f64| {
            let (mut lo, mut hi) = (0.0, stress / k);
            for _ in 0..80 {
                let gamma = 0.5 * (lo + hi);
                let value = k * gamma + yield_stress * (1.0 - (-m * gamma).exp());
                if value > stress {
                    hi = gamma;
                } else {
                    lo = gamma;
                }
            }
            0.5 * (lo + hi)
        };
        let a = shear(g * y.abs());
        let b = shear(g * half_height);
        (0.5 * k * (b * b - a * a)
            + yield_stress * ((a + 1.0 / m) * (-m * a).exp() - (b + 1.0 / m) * (-m * b).exp()))
            / g
    }

    /// G1: a real yielded-channel flow develops a nearly flat creeping core.
    /// Both native lattice dimensions share the same law and must agree with
    /// the continuum stress integral without any numerical viscosity cap.
    #[test]
    fn regularized_bingham_channel_matches_stress_integral_in_both_dimensions() {
        let (ny, g, k, yield_stress, m) = (33, 1e-5, 0.05, 0.00005, 15_000.0);
        let law = Rheology::HerschelBulkley {
            k,
            n: 1.0,
            yield_stress_over_density: yield_stress,
            regularization_time: m,
            minimum_shear_rate: 1e-12,
        };
        let max_steps = 60_000;
        let (grid2, steps2, stats2) = channel_flow(2, ny, g, law, max_steps);
        let (grid3, steps3, stats3) =
            crate::d3q19::plate_channel_flow3(2, 1, ny, g, law, max_steps);
        assert!(steps2 < max_steps && steps3 < max_steps);
        assert_eq!(
            (stats2.floored, stats2.capped, stats3.floored, stats3.capped),
            (0, 0, 0, 0)
        );
        let center = (ny - 1) as f64 / 2.0;
        let half_height = ny as f64 / 2.0;
        let peak = bingham_velocity(0.0, half_height, g, k, yield_stress, m);
        let mut errors = [0.0_f64; 2];
        for y in 0..ny {
            let exact = bingham_velocity(y as f64 - center, half_height, g, k, yield_stress, m);
            let u2 = grid2.moments(grid2.idx(0, y + 1)).u[0];
            let u3 = grid3.moments3(grid3.idx(0, 0, y + 1)).1[0];
            errors[0] = errors[0].max((u2 - exact).abs());
            errors[1] = errors[1].max((u3 - exact).abs());
            assert!((u2 - u3).abs() < 1e-9, "2D/3D constitutive drift at y={y}");
        }
        for error in errors {
            assert!(
                error / peak < 0.03,
                "Bingham profile error={}",
                error / peak
            );
        }
        let mid = ny / 2 + 1;
        let center_u = grid2.moments(grid2.idx(0, mid)).u[0];
        let core_edge_u = grid2.moments(grid2.idx(0, mid + 2)).u[0];
        assert!((center_u - core_edge_u) / center_u < 0.005);
        assert!(center_u < 0.6 * g * half_height * half_height / (2.0 * k));
        assert!((grid2.total_mass() - (2 * (ny + 2)) as f64).abs() < 1e-9);
        assert!((grid3.total_mass() - (2 * (ny + 2)) as f64).abs() < 1e-9);
        println!(
            "regularized Bingham: 2D={steps2} steps, 3D={steps3} steps, relative profile errors={:?}",
            errors.map(|e| e / peak),
        );
    }

    /// G1: n=5 is well outside the range where undamped local fixed-point
    /// iteration contracts near the wall. The actual flow must still settle
    /// to the independently integrated power-law momentum balance.
    #[test]
    fn strongly_thickening_channel_matches_analytic_profile() {
        let ny = 33;
        let force = 2.0 * 0.1 * 0.0017 / ny as f64;
        let n = 5.0;
        // Wall gamma=0.0017 and nu=0.1: low Mach, resolved viscous scale,
        // and tau=0.8 where the naive tau fixed-point slope has magnitude 1.5.
        let k = 0.1 / 0.0017_f64.powi(4);
        let (grid, steps, _) = channel_flow(2, ny, force, Rheology::PowerLaw { k, n }, 80_000);
        let mut error = 0.0_f64;
        let mut peak = 0.0_f64;
        for y in 1..=ny {
            let actual = grid.moments(grid.idx(0, y)).u[0];
            let exact = powerlaw_poiseuille_analytic(force, k, n, ny, y - 1);
            error = error.max((actual - exact).abs());
            peak = peak.max(exact);
        }
        assert!(
            steps < 80_000,
            "strongly thickening channel did not reach steady state"
        );
        assert!(
            error / peak < 0.03,
            "n=5 channel relative error={}",
            error / peak
        );
        println!(
            "n=5 channel: {} steps, relative profile error={:.6}%",
            steps,
            100.0 * error / peak
        );
    }
}
