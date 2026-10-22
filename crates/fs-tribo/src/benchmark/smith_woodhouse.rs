//! Smith & Woodhouse (2000) "plastic" thermal friction model of rosin.
//!
//! Source: J.H. Smith and J. Woodhouse, "The tribology of rosin", J. Mech.
//! Phys. Solids 48:1633-1681 (2000). Equation, table and section numbers
//! below are that paper's.
//!
//! **Constitutive law.** Eq. (21): `mu N = A k_y(T) sgn(v)`. The friction
//! force is set by a temperature-dependent shear yield stress `k_y` over
//! the contact area `A = pi a^2`. Sticking is set-valued with `|F| <= A k_y`.
//! `A` is proportional to `N` (section 6.4), so `a` scales as `N^{1/2}`.
//!
//! **Heat balance.** Eq. (6), here divided by `A`. With `T` the
//! temperature above ambient of the actively deforming layer, of thickness
//! `delta`:
//! `k_y |v| = |v| delta rho_r c_r T/(pi a) + rho_r delta c_r dT/dt + q_1 + q_2`.
//! - `q_1` and `q_2` are the conduction fluxes into the stationary and the
//!   moving bounding solids.
//! - They use the discrete one-dimensional half-space Green's function of
//!   Eqs. (10)-(14), per unit area:
//!   `g_0 = 2 beta h^{-1/2}`,
//!   `g_n = 2 beta h^{-1/2} (sqrt(n+1) - 2 sqrt(n) + sqrt(n-1))`,
//!   with `beta = K/sqrt(pi D)`.
//! - The stationary solid sees the whole temperature history.
//! - The moving solid sees only the time a point has been inside the contact
//!   (Eq. 16): `8a/(3 pi v)` at steady speed.
//!
//! **Calibration.** Section 5.3: `k_y(T)` is not measured. In steady sliding
//! the heat balance reduces to Eq. (20),
//! `mu N v = v a delta rho_r c_r T + T A K_2 sqrt(3v/(8 a D_2))`.
//! - With the steady-sliding fit `mu(v)` of Eq. (2), Eq. (20) gives `T(v)`,
//!   and Eq. (21) gives `k_y(v) = N mu(v)/A`.
//! - Together they define `k_y(T)` parametrically.
//! - That inversion needs `T(v)` to be increasing. Construction certifies
//!   this with an interval bound over the declared calibration range, and
//!   refuses otherwise.
//!
//! **Choices that are this implementation's, not the paper's.**
//! - *Time-varying speed.* The paper says only that the moving-surface limit
//!   "can be based on" the average time a point has spent in the contact.
//!   Here the window is the history within slip distance `8a/(3 pi)` of the
//!   step midpoint, and the pulse cut by that edge enters with its exact
//!   partial weight. At constant speed this is Eq. (16) without its integer
//!   rounding.
//! - *Time stepping.* The paper uses a forward difference (section 5.2).
//!   Here the losses are backward Euler and the friction uses the temperature
//!   at the start of the step, so a step is one linear solve. The scheme is
//!   first order and unconditionally stable. A forward difference with the
//!   `g_0` term is only stable for `h` below about `(rho_r delta c_r/beta)^2/10`.
//! - *Smoothing.* Section 6.4 smooths `sgn(v)` with `(2/pi) atan(v/v0)` for
//!   numerical stability. The held-effort solve here is set-valued, so no
//!   smoothing is applied.
//!
//! The time step is declared and uniform, as in Eqs. (10)-(13). A step of
//! any other length refuses.

use super::{admit_input, checked, outside};
use crate::tangential::{TangentialLaw, TangentialStep, TangentialStepInput};
use crate::{
    FrictionRegime, TriboError, checked_add, checked_mul, nonnegative_finite, positive_finite,
};
use core::f64::consts::PI;
use fs_math::det;
use std::sync::Arc;

/// Largest admitted history capacity, in steps.
pub const MAX_HISTORY_STEPS: usize = 1 << 22;
const CHUNK: usize = 1024;
const CERTIFICATE_PIECES: usize = 4096;

/// Dimensionless weight `w_n` of the discrete half-space Green's function,
/// Eq. (13) divided by `2 lambda h^{-1/2}`: `w_0 = 1` and
/// `w_n = sqrt(n+1) - 2 sqrt(n) + sqrt(n-1)`, evaluated without
/// cancellation as
/// `-2/((sqrt(n+1)+sqrt(n-1))(sqrt(n+1)+sqrt(n))(sqrt(n)+sqrt(n-1)))`.
#[must_use]
pub fn half_space_unit_weight(n: usize) -> f64 {
    if n == 0 {
        return 1.0;
    }
    #[allow(clippy::cast_precision_loss)]
    let x = n as f64;
    let (up, mid, down) = (det::sqrt(x + 1.0), det::sqrt(x), det::sqrt(x - 1.0));
    -2.0 / ((up + down) * (up + mid) * (mid + down))
}

/// Weight of a pulse of age `n` steps whose start is cut at the fraction
/// `phi` in `[0, 1)` of its step: the exact step-averaged half-space flux,
/// in units of `2 lambda h^{-1/2}`. At `phi = 0` it is `w_n`.
fn partial_unit_weight(n: usize, phi: f64) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let x = n as f64;
    if n == 0 {
        return det::sqrt(1.0 - phi);
    }
    1.0 / (det::sqrt(x + 1.0 - phi) + det::sqrt(x - phi))
        - 1.0 / (det::sqrt(x) + det::sqrt(x - 1.0))
}

/// Thermal properties of one solid bounding the rosin layer (Table 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalSolid {
    /// Thermal conductivity `K` [W/m/K], positive.
    pub conductivity_w_per_m_k: f64,
    /// Density [kg/m^3], positive.
    pub density_kg_per_m3: f64,
    /// Specific heat capacity [J/kg/K], positive.
    pub specific_heat_j_per_kg_k: f64,
}

impl ThermalSolid {
    fn validate(&self) -> Result<(), TriboError> {
        positive_finite(self.conductivity_w_per_m_k, "smith_woodhouse_conductivity")?;
        positive_finite(self.density_kg_per_m3, "smith_woodhouse_solid_density")?;
        positive_finite(
            self.specific_heat_j_per_kg_k,
            "smith_woodhouse_solid_specific_heat",
        )?;
        checked(self.half_space_coefficient(), "smith_woodhouse_half_space")?;
        Ok(())
    }

    /// Thermal diffusivity `D = K/(rho c)` [m^2/s] (Eq. 7).
    #[must_use]
    pub fn diffusivity_m2_per_s(&self) -> f64 {
        self.conductivity_w_per_m_k / (self.density_kg_per_m3 * self.specific_heat_j_per_kg_k)
    }

    /// `lambda/A = K/sqrt(pi D)` [W/m^2/K/s^{1/2}] (Eq. 14 per unit area).
    #[must_use]
    pub fn half_space_coefficient(&self) -> f64 {
        self.conductivity_w_per_m_k / det::sqrt(PI * self.diffusivity_m2_per_s())
    }
}

/// Contact of one rig: radius at a reference load, the actively deforming
/// layer, and the two bounding solids (section 4.5).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RosinContact {
    /// Contact radius `a` [m] at `reference_normal_n`, positive.
    pub radius_m: f64,
    /// Load [N] at which `radius_m` holds; `a` scales as `N^{1/2}`.
    pub reference_normal_n: f64,
    /// Thickness `delta` [m] of the actively deforming rosin layer, positive.
    pub layer_thickness_m: f64,
    /// Surface 1, stationary relative to the contact.
    pub stationary: ThermalSolid,
    /// Surface 2, moving through the contact at the sliding speed.
    pub moving: ThermalSolid,
}

impl RosinContact {
    fn validate(&self) -> Result<(), TriboError> {
        positive_finite(self.radius_m, "smith_woodhouse_radius")?;
        positive_finite(self.reference_normal_n, "smith_woodhouse_reference_normal")?;
        positive_finite(self.layer_thickness_m, "smith_woodhouse_layer_thickness")?;
        self.stationary.validate()?;
        self.moving.validate()
    }

    /// Contact radius [m] at `normal_n`, `a sqrt(N/N_ref)`.
    #[must_use]
    pub fn radius_at_m(&self, normal_n: f64) -> f64 {
        self.radius_m * det::sqrt(normal_n / self.reference_normal_n)
    }
}

/// Steady-sliding friction fit in the form of Eq. (2):
/// `mu(V) = c_1 e^{-V/V_1} + c_2 e^{-V/V_2} + c_inf`, with `V` the
/// sliding speed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteadySlidingFit {
    /// `c_1, c_2`, non-negative.
    pub amplitudes: [f64; 2],
    /// `V_1, V_2` [m/s], positive.
    pub decay_speeds_mps: [f64; 2],
    /// `c_inf`, positive.
    pub asymptote: f64,
}

impl SteadySlidingFit {
    fn validate(&self) -> Result<(), TriboError> {
        for c in self.amplitudes {
            nonnegative_finite(c, "smith_woodhouse_fit_amplitude")?;
        }
        for v in self.decay_speeds_mps {
            positive_finite(v, "smith_woodhouse_fit_speed")?;
        }
        positive_finite(self.asymptote, "smith_woodhouse_fit_asymptote")
    }

    /// Steady-sliding coefficient at speed `|v|` [m/s].
    #[must_use]
    pub fn coefficient(&self, v: f64) -> f64 {
        let v = v.abs();
        self.amplitudes[0] * det::exp(-v / self.decay_speeds_mps[0])
            + self.amplitudes[1] * det::exp(-v / self.decay_speeds_mps[1])
            + self.asymptote
    }
}

/// Every input of [`SmithWoodhouse`]. Nothing is defaulted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmithWoodhouseParams {
    /// Rosin density `rho_r` [kg/m^3], positive.
    pub rosin_density_kg_per_m3: f64,
    /// Rosin effective specific heat `c_r` [J/kg/K], positive.
    pub rosin_specific_heat_j_per_kg_k: f64,
    /// Steady-sliding curve used for calibration (Eq. 2).
    pub steady_fit: SteadySlidingFit,
    /// Steady-sliding rig; its `reference_normal_n` is the rig load.
    pub calibration: RosinContact,
    /// Upper end [m/s] of the calibration speed range. `k_y(T)` is defined
    /// for `0 <= T <= T_ss(max)`.
    pub calibration_speed_max_mps: f64,
    /// Contact being simulated.
    pub contact: RosinContact,
    /// Declared uniform time step `h` [s].
    pub step_s: f64,
    /// Steps of temperature history kept, at most [`MAX_HISTORY_STEPS`].
    pub history_capacity: usize,
}

/// Validated Smith & Woodhouse plastic thermal law.
#[derive(Debug, Clone, PartialEq)]
pub struct SmithWoodhouse {
    params: SmithWoodhouseParams,
    weights: Vec<f64>,
    temperature_limit_k: f64,
    // Calibration coefficients of `T_ss(w) = k(w) sqrt(w)/(alpha sqrt(w) + gamma)`.
    alpha: f64,
    gamma: f64,
}

impl SmithWoodhouse {
    /// Validate the parameters, certify that the calibration temperature
    /// increases with speed, and tabulate the Green's function weights.
    ///
    /// # Errors
    /// [`TriboError::InvalidInput`] naming the parameter, a capacity outside
    /// `1..=MAX_HISTORY_STEPS`, or an uncertified calibration.
    pub fn new(params: SmithWoodhouseParams) -> Result<Self, TriboError> {
        positive_finite(
            params.rosin_density_kg_per_m3,
            "smith_woodhouse_rosin_density",
        )?;
        positive_finite(
            params.rosin_specific_heat_j_per_kg_k,
            "smith_woodhouse_rosin_specific_heat",
        )?;
        params.steady_fit.validate()?;
        params.calibration.validate()?;
        params.contact.validate()?;
        positive_finite(
            params.calibration_speed_max_mps,
            "smith_woodhouse_calibration_speed_max",
        )?;
        positive_finite(params.step_s, "dt_s")?;
        if !(1..=MAX_HISTORY_STEPS).contains(&params.history_capacity) {
            return Err(TriboError::InvalidInput {
                field: "smith_woodhouse_history_capacity",
            });
        }
        let c = &params.calibration;
        let alpha = checked(
            c.layer_thickness_m
                * params.rosin_density_kg_per_m3
                * params.rosin_specific_heat_j_per_kg_k
                / (PI * c.radius_m),
            "smith_woodhouse_calibration_advection",
        )?;
        let gamma = c.moving.conductivity_w_per_m_k
            * det::sqrt(3.0 / (8.0 * c.radius_m * c.moving.diffusivity_m2_per_s()));
        positive_finite(gamma, "smith_woodhouse_calibration_conduction")?;
        let mut law = Self {
            params,
            weights: Vec::new(),
            temperature_limit_k: 0.0,
            alpha,
            gamma,
        };
        law.certify_calibration()?;
        law.temperature_limit_k = checked(
            law.steady_temperature(params.calibration_speed_max_mps),
            "smith_woodhouse_temperature_limit",
        )?;
        law.weights = (0..=params.history_capacity)
            .map(half_space_unit_weight)
            .collect();
        Ok(law)
    }

    /// The validated parameters.
    #[must_use]
    pub const fn params(&self) -> &SmithWoodhouseParams {
        &self.params
    }

    /// `T_ss(calibration_speed_max)` [K above ambient], the largest
    /// temperature at which `k_y` is defined.
    #[must_use]
    pub const fn temperature_limit_k(&self) -> f64 {
        self.temperature_limit_k
    }

    /// Calibration yield stress `k_y(w) = N mu(w)/(pi a^2)` [Pa] at speed `w`.
    fn calibration_yield(&self, w: f64) -> f64 {
        let c = &self.params.calibration;
        c.reference_normal_n * self.params.steady_fit.coefficient(w)
            / (PI * c.radius_m * c.radius_m)
    }

    /// Eq. (20) solved for `T` at speed `w >= 0`; zero at `w = 0`. With
    /// `alpha = delta rho_r c_r/(pi a)` and `gamma = K_2 sqrt(3/(8 a D_2))`
    /// it reads `T = k_y(w) sqrt(w)/(alpha sqrt(w) + gamma)`.
    fn steady_temperature(&self, w: f64) -> f64 {
        let root = det::sqrt(w);
        self.calibration_yield(w) * root / (self.alpha * root + self.gamma)
    }

    /// Steady-sliding temperature [K above ambient] of the calibration rig
    /// at speed `|v|`, from Eq. (20) with the fit of Eq. (2).
    ///
    /// # Errors
    /// `|v|` above the calibration range.
    pub fn calibration_steady_temperature_k(&self, v: f64) -> Result<f64, TriboError> {
        let w = checked(v, "slip_speed_mps")?.abs();
        let max = self.params.calibration_speed_max_mps;
        if w > max {
            return Err(outside("smith_woodhouse_calibration_speed", w, [0.0, max]));
        }
        checked(
            self.steady_temperature(w),
            "smith_woodhouse_steady_temperature",
        )
    }

    /// `d T_ss/dw > 0` on `(0, w_max]` is equivalent to
    /// `2 (-w mu')(1 + (alpha/gamma) sqrt(w)) < mu`. On each piece of a
    /// uniform partition the left side is bounded above (each `x e^{-x}` by
    /// its interval maximum) and `mu` below (by its value at the right end).
    fn certify_calibration(&self) -> Result<(), TriboError> {
        let fit = &self.params.steady_fit;
        let max = self.params.calibration_speed_max_mps;
        let x_exp_max = |lo: f64, hi: f64| -> f64 {
            if lo <= 1.0 && 1.0 <= hi {
                det::exp(-1.0)
            } else if hi < 1.0 {
                hi * det::exp(-hi)
            } else {
                lo * det::exp(-lo)
            }
        };
        #[allow(clippy::cast_precision_loss)]
        let pieces = CERTIFICATE_PIECES as f64;
        for i in 0..CERTIFICATE_PIECES {
            #[allow(clippy::cast_precision_loss)]
            let (lo, hi) = (max * i as f64 / pieces, max * (i + 1) as f64 / pieces);
            let slope: f64 = fit
                .amplitudes
                .iter()
                .zip(fit.decay_speeds_mps)
                .map(|(c, v)| c * x_exp_max(lo / v, hi / v))
                .sum();
            let bound = 2.0 * slope * (1.0 + self.alpha / self.gamma * det::sqrt(hi));
            if !(bound < fit.coefficient(hi)) {
                return Err(TriboError::InvalidInput {
                    field: "smith_woodhouse_calibration_temperature_not_certified_increasing",
                });
            }
        }
        Ok(())
    }

    /// Calibration speed `w` with `T_ss(w) = t`, by bisection to adjacent
    /// floats; `T_ss` is certified increasing.
    fn calibration_speed_for(&self, t: f64) -> Result<f64, TriboError> {
        nonnegative_finite(t, "smith_woodhouse_temperature_rise")?;
        if t > self.temperature_limit_k {
            return Err(outside(
                "smith_woodhouse_temperature_rise",
                t,
                [0.0, self.temperature_limit_k],
            ));
        }
        if t == 0.0 {
            return Ok(0.0);
        }
        let (mut lo, mut hi) = (0.0, self.params.calibration_speed_max_mps);
        loop {
            let mid = f64::midpoint(lo, hi);
            if !(lo < mid && mid < hi) {
                return Ok(hi);
            }
            if self.steady_temperature(mid) < t {
                lo = mid;
            } else {
                hi = mid;
            }
        }
    }

    /// Shear yield stress `k_y(T)` [Pa] at a temperature `t` [K above
    /// ambient], calibrated through Eqs. (2), (20) and (21).
    ///
    /// # Errors
    /// A negative or non-finite `t`, or `t` above [`Self::temperature_limit_k`].
    pub fn yield_stress_pa(&self, t: f64) -> Result<f64, TriboError> {
        let w = self.calibration_speed_for(t)?;
        checked(self.calibration_yield(w), "smith_woodhouse_yield_stress")
    }

    /// Weighted history sum `sum_{k=from}^{m-1} w_{m-k} T_k`.
    fn history_sum(&self, history: &History, from: usize) -> f64 {
        let m = history.len();
        (from..m)
            .map(|k| self.weights[m - k] * history.get(k).temperature_rise_k)
            .sum()
    }
}

/// One step of history: the temperature over the step and the slip in it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Pulse {
    temperature_rise_k: f64,
    slip_m: f64,
}

/// Append-only history in shared sealed chunks, so a candidate state clones
/// in `O(len/CHUNK + CHUNK)`.
#[derive(Clone, PartialEq, Default)]
struct History {
    sealed: Vec<Arc<[Pulse]>>,
    open: Vec<Pulse>,
}

impl History {
    fn len(&self) -> usize {
        self.sealed.len() * CHUNK + self.open.len()
    }
    fn get(&self, k: usize) -> Pulse {
        let sealed = self.sealed.len() * CHUNK;
        if k < sealed {
            self.sealed[k / CHUNK][k % CHUNK]
        } else {
            self.open[k - sealed]
        }
    }
    fn push(&mut self, pulse: Pulse) {
        self.open.push(pulse);
        if self.open.len() == CHUNK {
            self.sealed.push(Arc::from(core::mem::take(&mut self.open)));
        }
    }
}

impl core::fmt::Debug for History {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let last = self.len().checked_sub(1).map(|k| self.get(k));
        f.debug_struct("History")
            .field("len", &self.len())
            .field("last", &last)
            .finish()
    }
}

/// Smith & Woodhouse checkpoint: the temperature and slip history from
/// the ambient start, and the accumulated slip.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SmithWoodhouseState {
    history: History,
    slip_m: [f64; 2],
}

impl SmithWoodhouseState {
    /// Contact temperature [K above ambient] over the last step; zero at
    /// the ambient start.
    #[must_use]
    pub fn temperature_rise_k(&self) -> f64 {
        self.history
            .len()
            .checked_sub(1)
            .map_or(0.0, |k| self.history.get(k).temperature_rise_k)
    }
    /// Steps taken.
    #[must_use]
    pub fn steps(&self) -> usize {
        self.history.len()
    }
    /// Accumulated solved slip [m].
    #[must_use]
    pub const fn slip_m(&self) -> [f64; 2] {
        self.slip_m
    }
}

impl TangentialLaw for SmithWoodhouse {
    type State = SmithWoodhouseState;

    fn model_id(&self) -> &'static str {
        "benchmark-smith-woodhouse-plastic-thermal"
    }

    /// The paper's initial condition: rosin at ambient temperature.
    fn initial_state(&self) -> Self::State {
        SmithWoodhouseState::default()
    }

    fn stored_energy_j(&self, _s: &Self::State) -> f64 {
        0.0
    }

    #[allow(clippy::too_many_lines)]
    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let p = &self.params;
        let distance = admit_input(x)?;
        positive_finite(x.normal_n, "smith_woodhouse_normal_force")?;
        if x.dt_s.to_bits() != p.step_s.to_bits() {
            return Err(TriboError::InvalidInput {
                field: "smith_woodhouse_step_differs_from_declared_step",
            });
        }
        let m = s.history.len();
        if m >= p.history_capacity {
            #[allow(clippy::cast_precision_loss)]
            return Err(outside(
                "smith_woodhouse_history_steps",
                (m + 1) as f64,
                [1.0, p.history_capacity as f64],
            ));
        }
        let h = p.step_s;
        let t_prev = s.temperature_rise_k();
        let yield_pa = self.yield_stress_pa(t_prev)?;
        let radius = checked(p.contact.radius_at_m(x.normal_n), "smith_woodhouse_radius")?;
        let area = checked(PI * radius * radius, "smith_woodhouse_area")?;
        let limit = checked_mul(area, yield_pa, "static_limit")?;
        let c = x.compliance_m_per_n;
        let (force, slip, regime) = if distance == 0.0 {
            (0.0, 0.0, FrictionRegime::Sticking)
        } else if c == 0.0 {
            (limit, distance, FrictionRegime::Sliding)
        } else {
            let required = checked(distance / c, "sticking_force_n")?;
            if required <= limit {
                (required, 0.0, FrictionRegime::Sticking)
            } else {
                let slip = (distance - c * limit).max(0.0);
                let regime = if slip > 0.0 {
                    FrictionRegime::Sliding
                } else {
                    FrictionRegime::Sticking
                };
                (limit, slip, regime)
            }
        };
        if force > x.solve.maximum_force_n {
            return Err(outside(
                "maximum_force_n",
                force,
                [0.0, x.solve.maximum_force_n],
            ));
        }

        // Moving-surface window: the history within slip `8a/(3 pi)` of the
        // step midpoint. `window` is (oldest pulse index, cut fraction).
        let reach = 8.0 * radius / (3.0 * PI);
        let mut window = None;
        let mut remaining = reach - 0.5 * slip;
        if remaining <= 0.0 {
            window = Some((m, 0.5 - reach / slip));
        } else {
            for k in (0..m).rev() {
                let s_k = s.history.get(k).slip_m;
                if s_k > 0.0 && s_k >= remaining {
                    window = Some((k, 1.0 - remaining / s_k));
                    break;
                }
                remaining -= s_k;
            }
        }
        let stationary = p.contact.stationary.half_space_coefficient();
        let moving = p.contact.moving.half_space_coefficient();
        let scale = 2.0 / det::sqrt(h);
        let all_history = self.history_sum(&s.history, 0);
        let (moving_self, moving_known) = match window {
            None => (1.0, all_history),
            Some((k, phi)) if k == m => (partial_unit_weight(0, phi), 0.0),
            Some((k, phi)) => (
                1.0,
                self.history_sum(&s.history, k + 1)
                    + partial_unit_weight(m - k, phi) * s.history.get(k).temperature_rise_k,
            ),
        };
        let capacity = p.rosin_density_kg_per_m3
            * p.contact.layer_thickness_m
            * p.rosin_specific_heat_j_per_kg_k;
        let speed = slip / h;
        let advection = speed * capacity / (PI * radius);
        let heating = if slip > 0.0 { yield_pa * speed } else { 0.0 };
        let diagonal = capacity / h + advection + scale * (stationary + moving * moving_self);
        let known = scale * (stationary * all_history + moving * moving_known);
        let t_next = checked(
            (capacity * t_prev / h + heating - known) / diagonal,
            "smith_woodhouse_temperature_rise",
        )?;
        nonnegative_finite(t_next, "smith_woodhouse_temperature_rise")?;

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
            slip_increment_m[i] = slip * direction[i];
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
        let mut history = s.history.clone();
        history.push(Pulse {
            temperature_rise_k: t_next,
            slip_m: slip,
        });
        Ok(TangentialStep {
            traction_n,
            slip_increment_m,
            state: SmithWoodhouseState { history, slip_m },
            regime,
            stored_change_j: 0.0,
            dissipated_j: work,
            closure_residual_j: 0.0,
            compatibility_residual_m,
            iterations: 0,
        })
    }
}
