//! Benchmark tangential laws (IRP M4.2-W4): rate-and-state friction
//! (Dieterich 1979; Ruina 1983) with the aging and slip evolution laws, and
//! LuGre (Canudas de Wit, Olsson, Astrom & Lischinsky 1995). They exist only
//! under `cfg(test)` or the non-default `benchmark-laws` feature: named
//! comparisons for experiments such as X-Bw2, never render laws.
//!
//! **Held-effort discretization.** Over one interval the relative speed
//! `V = |du|/dt` is constant, so both state equations integrate exactly:
//! - aging law `d theta/dt = 1 - V theta/Dc`:
//!   `theta_1 = Dc/V + (theta_0 - Dc/V) e^{-V dt/Dc}`;
//! - slip law `d theta/dt = -(V theta/Dc) ln(V theta/Dc)`:
//!   `ln(V theta_1/Dc) = ln(V theta_0/Dc) e^{-V dt/Dc}`;
//! - LuGre `dz/dt = v - sigma0 |v| z/g(v)`:
//!   `z_1 = z_ss + (z_0 - z_ss) e^{-sigma0 |v| dt/g(v)}`, `z_ss = g(v) sgn(v)/sigma0`.
//!
//! The decays use `expm1`, so tiny steps keep full relative accuracy.
//!
//! Rate-and-state friction over an interval uses the state at its start:
//! `mu = mu0 + a ln(V/V0) + b ln(V0 theta_0/Dc)` is then strictly increasing
//! in `V` for `a > 0`, so the compliant root is unique, and a prescribed-speed
//! history sampled at interval starts is exact. LuGre carries
//! `F = sigma0 z_1 + sigma1 (z_1 - z_0)/dt + sigma2 v` with the bristle
//! energy `sigma0 z^2/2` as reversible storage. A step whose port
//! dissipation would be negative refuses, because passivity of the full
//! force map depends on `sigma1` and is not assumed here.
//!
//! **Not implemented.** The Smith & Woodhouse (2000) lumped thermal model.
//! Its constitutive form and fitted parameters must come from that paper (or
//! Woodhouse 2003), which was not accessible when this module was written; no
//! form is guessed.

use crate::tangential::{TangentialLaw, TangentialStep, TangentialStepInput};
use crate::{
    FrictionRegime, TriboError, checked_add, checked_mul, finite, nonnegative_finite,
    positive_finite,
};
use fs_math::det;

fn checked(value: f64, field: &'static str) -> Result<f64, TriboError> {
    finite(value, field)?;
    Ok(value)
}
fn outside(field: &'static str, value: f64, range: [f64; 2]) -> TriboError {
    TriboError::OutsideApplicability {
        field,
        value,
        minimum: range[0],
        maximum: range[1],
    }
}

/// Validated common input; returns the released slip distance [m].
fn admit_input(x: &TangentialStepInput) -> Result<f64, TriboError> {
    x.solve.validate()?;
    nonnegative_finite(x.normal_n, "normal_force")?;
    nonnegative_finite(x.compliance_m_per_n, "compliance_m_per_n")?;
    positive_finite(x.dt_s, "dt_s")?;
    for value in x.slip_increment_m {
        finite(value, "slip_increment_m")?;
    }
    checked(
        if x.slip_increment_m[1] == 0.0 {
            x.slip_increment_m[0].abs()
        } else {
            det::hypot(x.slip_increment_m[0], x.slip_increment_m[1])
        },
        "free_slip_distance_m",
    )
}

/// Bisection for the unique root of a residual that increases with the
/// traction, inside `[lo, hi]`, with the caller's budget and tolerance.
fn bisect(
    x: &TangentialStepInput,
    mut lo: f64,
    mut hi: f64,
    residual: impl Fn(f64) -> Result<(f64, f64), TriboError>,
) -> Result<(f64, usize), TriboError> {
    let (r_lo, t_lo) = residual(lo)?;
    if r_lo.abs() <= t_lo {
        return Ok((lo, 0));
    }
    let (r_hi, t_hi) = residual(hi)?;
    if r_hi.abs() <= t_hi {
        return Ok((hi, 0));
    }
    if r_lo > 0.0 || r_hi < 0.0 {
        return Err(TriboError::InvalidInput {
            field: "benchmark traction root outside the admitted bracket",
        });
    }
    let (mut r, mut tol) = (r_hi, t_hi);
    let mut used = 0;
    for iteration in 1..=x.solve.max_iterations {
        used = iteration;
        let mid = f64::midpoint(lo, hi);
        (r, tol) = residual(mid)?;
        if r.abs() <= tol {
            return Ok((mid, iteration));
        }
        // No representable traction lies strictly between the endpoints.
        if !(lo < mid && mid < hi) {
            break;
        }
        if r < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(TriboError::TangentialSolve {
        residual_n: r,
        tolerance_n: tol,
        iterations: used,
    })
}

/// State evolution law of rate-and-state friction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Evolution {
    /// Dieterich aging law `d theta/dt = 1 - V theta/Dc`.
    Aging,
    /// Ruina slip law `d theta/dt = -(V theta/Dc) ln(V theta/Dc)`.
    Slip,
}

/// Dieterich-Ruina rate-and-state friction, `mu = mu0 + a ln(V/V0) + b
/// ln(V0 theta/Dc)`. Its steady state is `mu0 + (a - b) ln(V/V0)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RateState {
    /// Reference coefficient at `V0` in steady state.
    pub mu0: f64,
    /// Direct-effect parameter, `a >= 0`.
    pub a: f64,
    /// Evolution-effect parameter, `b >= 0`.
    pub b: f64,
    /// Reference slip speed [m/s], positive.
    pub v0_mps: f64,
    /// Characteristic slip distance [m], positive.
    pub dc_m: f64,
    /// State evolution law.
    pub evolution: Evolution,
    /// Inclusive admitted slip-speed range [m/s]; the log law is undefined at
    /// zero speed, so `0 < min < max` and slip outside it refuses.
    pub speed_range_mps: [f64; 2],
}

impl RateState {
    /// Admit the parameters: `a, b >= 0`, `V0, Dc > 0`, `0 < min < max`.
    ///
    /// # Errors
    /// [`TriboError::InvalidInput`] naming the parameter.
    pub fn validate(&self) -> Result<(), TriboError> {
        finite(self.mu0, "rate_state_mu0")?;
        nonnegative_finite(self.a, "rate_state_a")?;
        nonnegative_finite(self.b, "rate_state_b")?;
        positive_finite(self.v0_mps, "rate_state_v0")?;
        positive_finite(self.dc_m, "rate_state_dc")?;
        let [lo, hi] = self.speed_range_mps;
        positive_finite(lo, "rate_state_speed_min")?;
        positive_finite(hi, "rate_state_speed_max")?;
        if hi <= lo {
            return Err(TriboError::InvalidInput {
                field: "rate_state_speed_range",
            });
        }
        Ok(())
    }

    fn admit_speed(&self, v: f64) -> Result<(), TriboError> {
        let [lo, hi] = self.speed_range_mps;
        if !(lo..=hi).contains(&v) {
            return Err(outside("rate_state_slip_speed", v, self.speed_range_mps));
        }
        Ok(())
    }

    /// Friction coefficient at speed `v` [m/s] and state `theta` [s].
    ///
    /// # Errors
    /// A speed outside the admitted range or a non-positive state.
    pub fn coefficient(&self, v: f64, theta: f64) -> Result<f64, TriboError> {
        self.validate()?;
        self.admit_speed(v)?;
        positive_finite(theta, "rate_state_theta")?;
        checked(
            self.mu0
                + self.a * det::ln(v / self.v0_mps)
                + self.b * det::ln(self.v0_mps * theta / self.dc_m),
            "rate_state_mu",
        )
    }

    /// Steady-state coefficient `mu0 + (a - b) ln(V/V0)`.
    ///
    /// # Errors
    /// A speed outside the admitted range.
    pub fn steady_coefficient(&self, v: f64) -> Result<f64, TriboError> {
        self.validate()?;
        self.admit_speed(v)?;
        checked(
            self.mu0 + (self.a - self.b) * det::ln(v / self.v0_mps),
            "rate_state_mu_ss",
        )
    }

    /// Exact state after `dt` [s] at constant speed `v` [m/s]. At `v = 0` the
    /// aging law heals linearly and the slip law holds its state.
    ///
    /// # Errors
    /// Non-finite or non-positive inputs.
    pub fn evolve(&self, theta: f64, v: f64, dt: f64) -> Result<f64, TriboError> {
        self.validate()?;
        positive_finite(theta, "rate_state_theta")?;
        nonnegative_finite(v, "rate_state_slip_speed")?;
        positive_finite(dt, "dt_s")?;
        if v == 0.0 {
            return checked(
                match self.evolution {
                    Evolution::Aging => theta + dt,
                    Evolution::Slip => theta,
                },
                "rate_state_theta",
            );
        }
        let decay = det::expm1(-v * dt / self.dc_m);
        let next = match self.evolution {
            Evolution::Aging => {
                let steady = self.dc_m / v;
                theta - (steady - theta) * decay
            }
            Evolution::Slip => {
                let phi = det::ln(v * theta / self.dc_m);
                self.dc_m / v * det::exp(phi + phi * decay)
            }
        };
        let next = checked(next, "rate_state_theta")?;
        positive_finite(next, "rate_state_theta")?;
        Ok(next)
    }

    /// Ruina (1983) quasi-static spring-slider critical stiffness
    /// `k_c = normal (b - a)/Dc` [N/m]; steady sliding is unstable below it
    /// for a velocity-weakening interface (`b > a`).
    #[must_use]
    pub fn critical_stiffness_n_per_m(&self, normal_n: f64) -> f64 {
        normal_n * (self.b - self.a) / self.dc_m
    }
}

/// Rate-and-state checkpoint: the state variable and accumulated slip.
#[derive(Debug, Clone, PartialEq)]
pub struct RateStateState {
    theta_s: f64,
    slip_m: [f64; 2],
}

impl RateStateState {
    /// A state with `theta` [s] and zero accumulated slip.
    ///
    /// # Errors
    /// A non-positive or non-finite `theta`.
    pub fn new(theta_s: f64) -> Result<Self, TriboError> {
        positive_finite(theta_s, "rate_state_theta")?;
        Ok(Self {
            theta_s,
            slip_m: [0.0; 2],
        })
    }
    /// State variable [s].
    #[must_use]
    pub const fn theta_s(&self) -> f64 {
        self.theta_s
    }
    /// Accumulated solved slip [m].
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }
}

impl TangentialLaw for RateState {
    type State = RateStateState;

    fn model_id(&self) -> &str {
        match self.evolution {
            Evolution::Aging => "benchmark-rate-state-aging",
            Evolution::Slip => "benchmark-rate-state-slip",
        }
    }

    /// Steady sliding at the reference speed, `theta = Dc/V0`.
    fn initial_state(&self) -> Self::State {
        RateStateState {
            theta_s: self.dc_m / self.v0_mps,
            slip_m: [0.0; 2],
        }
    }

    fn stored_energy_j(&self, _s: &Self::State) -> f64 {
        0.0
    }

    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.validate()?;
        let distance = admit_input(x)?;
        let c = x.compliance_m_per_n;
        let (force, speed, iterations) = if x.normal_n == 0.0 {
            (0.0, distance / x.dt_s, 0)
        } else if c == 0.0 {
            let speed = distance / x.dt_s;
            let force = checked_mul(
                self.coefficient(speed, s.theta_s)?,
                x.normal_n,
                "traction_n",
            )?;
            (force, speed, 0)
        } else {
            let [vmin, vmax] = self.speed_range_mps;
            let lo = ((distance - vmax * x.dt_s) / c).max(0.0);
            let hi = ((distance - vmin * x.dt_s) / c).min(x.solve.maximum_force_n);
            if !(hi >= lo) {
                return Err(outside(
                    "rate_state_slip_speed",
                    distance / x.dt_s,
                    self.speed_range_mps,
                ));
            }
            let speed_at = |f: f64| ((distance - c * f) / x.dt_s).clamp(vmin, vmax);
            let residual = |f: f64| -> Result<(f64, f64), TriboError> {
                let expected = self.coefficient(speed_at(f), s.theta_s)? * x.normal_n;
                Ok((
                    checked(f - expected, "constitutive_residual_n")?,
                    x.solve.tolerance(f, expected)?,
                ))
            };
            // A root beyond either bracket end needs a speed outside the
            // admitted range (a slip event or arrest), refused by name.
            let (r_lo, t_lo) = residual(lo)?;
            let (r_hi, t_hi) = residual(hi)?;
            if r_lo > t_lo && lo > 0.0 {
                return Err(outside("rate_state_slip_speed", vmax, self.speed_range_mps));
            }
            if r_hi < -t_hi && hi < x.solve.maximum_force_n {
                return Err(outside("rate_state_slip_speed", vmin, self.speed_range_mps));
            }
            let (force, iterations) = bisect(x, lo, hi, residual)?;
            (force, speed_at(force), iterations)
        };
        if force < 0.0 {
            return Err(TriboError::InvalidInput {
                field: "negative rate-state friction coefficient",
            });
        }
        let theta_s = self.evolve(s.theta_s, speed, x.dt_s)?;
        let direction = if distance == 0.0 {
            [0.0; 2]
        } else {
            [
                x.slip_increment_m[0] / distance,
                x.slip_increment_m[1] / distance,
            ]
        };
        let mut traction_n = [0.0; 2];
        let mut slip_increment_m = [0.0; 2];
        let mut compatibility_residual_m = [0.0; 2];
        let mut slip_m = s.slip_m;
        for i in 0..2 {
            traction_n[i] = checked_mul(force, direction[i], "traction_n")?;
            let actual = checked(x.slip_increment_m[i] - c * traction_n[i], "resolved_slip_m")?;
            slip_increment_m[i] = speed * x.dt_s * direction[i];
            compatibility_residual_m[i] =
                checked(slip_increment_m[i] - actual, "compatibility_residual_m")?;
            slip_m[i] = checked_add(s.slip_m[i], slip_increment_m[i], "accumulated_slip_m")?;
        }
        let work = checked_add(
            traction_n[0] * slip_increment_m[0],
            traction_n[1] * slip_increment_m[1],
            "tangent_work_j",
        )?;
        nonnegative_finite(work, "tangent_work_j")?;
        Ok(TangentialStep {
            traction_n,
            slip_increment_m,
            state: RateStateState { theta_s, slip_m },
            regime: FrictionRegime::Sliding,
            stored_change_j: 0.0,
            dissipated_j: work,
            closure_residual_j: 0.0,
            compatibility_residual_m,
            iterations,
        })
    }
}

/// LuGre friction (Canudas de Wit et al. 1995) in its classic fixed-load
/// force form: `g(v) = Fc + (Fs - Fc) e^{-(v/vs)^2}`,
/// `F = sigma0 z + sigma1 dz/dt + sigma2 v`. One-dimensional: slip along the
/// second tangent axis refuses. `normal_n` does not enter the classic form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LuGre {
    /// Bristle stiffness [N/m], positive.
    pub sigma0_n_per_m: f64,
    /// Bristle damping [N s/m], non-negative.
    pub sigma1_ns_per_m: f64,
    /// Viscous coefficient [N s/m], non-negative.
    pub sigma2_ns_per_m: f64,
    /// Coulomb level [N], positive.
    pub coulomb_n: f64,
    /// Static level [N], at least the Coulomb level.
    pub static_n: f64,
    /// Stribeck speed [m/s], positive.
    pub stribeck_mps: f64,
}

impl LuGre {
    /// Admit the parameters.
    ///
    /// # Errors
    /// [`TriboError::InvalidInput`] naming the parameter.
    pub fn validate(&self) -> Result<(), TriboError> {
        positive_finite(self.sigma0_n_per_m, "lugre_sigma0")?;
        nonnegative_finite(self.sigma1_ns_per_m, "lugre_sigma1")?;
        nonnegative_finite(self.sigma2_ns_per_m, "lugre_sigma2")?;
        positive_finite(self.coulomb_n, "lugre_coulomb")?;
        positive_finite(self.static_n, "lugre_static")?;
        positive_finite(self.stribeck_mps, "lugre_stribeck_speed")?;
        if self.static_n < self.coulomb_n {
            return Err(TriboError::InvalidInput {
                field: "lugre_static_below_coulomb",
            });
        }
        Ok(())
    }

    /// Stribeck curve `g(v)` [N].
    #[must_use]
    pub fn stribeck_n(&self, v: f64) -> f64 {
        let r = v / self.stribeck_mps;
        self.coulomb_n + (self.static_n - self.coulomb_n) * det::exp(-r * r)
    }

    /// Steady-state force `g(v) sgn(v) + sigma2 v` [N].
    #[must_use]
    pub fn steady_force_n(&self, v: f64) -> f64 {
        let sign = if v > 0.0 {
            1.0
        } else if v < 0.0 {
            -1.0
        } else {
            0.0
        };
        self.stribeck_n(v) * sign + self.sigma2_ns_per_m * v
    }

    /// Exact bristle deflection after `dt` at constant velocity `v`.
    #[must_use]
    pub fn bristle_after(&self, z: f64, v: f64, dt: f64) -> f64 {
        if v == 0.0 {
            return z;
        }
        let g = self.stribeck_n(v);
        let steady = g * v.signum() / self.sigma0_n_per_m;
        z - (steady - z) * det::expm1(-self.sigma0_n_per_m * v.abs() * dt / g)
    }

    /// Exact `int z v dt` over one constant-velocity interval (the z-map work
    /// in Canudas de Wit's passivity property).
    #[must_use]
    pub fn bristle_work(&self, z: f64, v: f64, dt: f64) -> f64 {
        if v == 0.0 {
            return 0.0;
        }
        let g = self.stribeck_n(v);
        let steady = g * v.signum() / self.sigma0_n_per_m;
        let rate = self.sigma0_n_per_m * v.abs() / g;
        v * (steady * dt - (z - steady) * det::expm1(-rate * dt) / rate)
    }

    fn force(&self, z0: f64, v: f64, dt: f64) -> (f64, f64) {
        let z1 = self.bristle_after(z0, v, dt);
        (
            self.sigma0_n_per_m * z1
                + self.sigma1_ns_per_m * (z1 - z0) / dt
                + self.sigma2_ns_per_m * v,
            z1,
        )
    }
}

/// LuGre checkpoint: bristle deflection and accumulated slip.
#[derive(Debug, Clone, PartialEq)]
pub struct LuGreState {
    z_m: f64,
    slip_m: [f64; 2],
}

impl LuGreState {
    /// A state with bristle deflection `z` [m] and zero accumulated slip.
    ///
    /// # Errors
    /// A non-finite `z`.
    pub fn new(z_m: f64) -> Result<Self, TriboError> {
        finite(z_m, "lugre_bristle_m")?;
        Ok(Self {
            z_m,
            slip_m: [0.0; 2],
        })
    }
    /// Bristle deflection [m].
    #[must_use]
    pub const fn z_m(&self) -> f64 {
        self.z_m
    }
    /// Accumulated slip [m].
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }
}

impl TangentialLaw for LuGre {
    type State = LuGreState;

    fn model_id(&self) -> &'static str {
        "benchmark-lugre"
    }

    fn initial_state(&self) -> Self::State {
        LuGreState {
            z_m: 0.0,
            slip_m: [0.0; 2],
        }
    }

    fn stored_energy_j(&self, s: &Self::State) -> f64 {
        0.5 * self.sigma0_n_per_m * s.z_m * s.z_m
    }

    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.validate()?;
        admit_input(x)?;
        if x.slip_increment_m[1] != 0.0 {
            return Err(TriboError::InvalidInput {
                field: "lugre benchmark is one-dimensional",
            });
        }
        let (free, c, dt) = (x.slip_increment_m[0], x.compliance_m_per_n, x.dt_s);
        let velocity_at = |f: f64| (free - c * f) / dt;
        let (force, z1, iterations) = if c == 0.0 {
            let (force, z1) = self.force(s.z_m, free / dt, dt);
            (force, z1, 0)
        } else {
            let cap = x.solve.maximum_force_n;
            let (force, iterations) = bisect(x, -cap, cap, |f| {
                let (expected, _) = self.force(s.z_m, velocity_at(f), dt);
                Ok((
                    checked(f - expected, "constitutive_residual_n")?,
                    x.solve.tolerance(f, expected)?,
                ))
            })?;
            let (_, z1) = self.force(s.z_m, velocity_at(force), dt);
            (force, z1, iterations)
        };
        let force = checked(force, "traction_n")?;
        let z1 = checked(z1, "lugre_bristle_m")?;
        let v = velocity_at(force);
        let slip = checked(v * dt, "resolved_slip_m")?;
        let actual = checked(free - c * force, "resolved_slip_m")?;
        let state = LuGreState {
            z_m: z1,
            slip_m: [
                checked_add(s.slip_m[0], slip, "accumulated_slip_m")?,
                s.slip_m[1],
            ],
        };
        let stored_change_j = self.stored_energy_j(&state) - self.stored_energy_j(s);
        let work = checked(force * slip, "tangent_work_j")?;
        let dissipated_j = work - stored_change_j;
        let scale = work.abs().max(stored_change_j.abs());
        if dissipated_j < -1e-12 * scale {
            return Err(TriboError::InvalidInput {
                field: "lugre step is not dissipative (sigma1 passivity not assumed)",
            });
        }
        Ok(TangentialStep {
            traction_n: [force, 0.0],
            slip_increment_m: [slip, 0.0],
            state,
            regime: if v == 0.0 {
                FrictionRegime::Sticking
            } else {
                FrictionRegime::Sliding
            },
            stored_change_j,
            dissipated_j: dissipated_j.max(0.0),
            closure_residual_j: stored_change_j + dissipated_j.max(0.0) - work,
            compatibility_residual_m: [slip - actual, 0.0],
            iterations,
        })
    }
}
