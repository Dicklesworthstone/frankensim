//! Smith & Woodhouse (2000) "plastic" thermal friction model of rosin.
//!
//! **A render law.** It is the one law promoted out of `benchmark` (2026-11-06),
//! because it is the bowed-string thermal rosin model, for the fs-couple
//! bowed-string rosin row. `benchmark` re-exports it under the same
//! names, so the X-Bw2 bakeoff and the g1 tests are unchanged. Every input is
//! the caller's: nothing here admits a material or carries bow geometry.
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
//! - With a steady-sliding curve `mu(v)`, Eq. (20) gives `T(v)`, and Eq. (21)
//!   gives `k_y(v) = N mu(v)/A`. The curve is the paper's Eq. (2) fit or the
//!   shifted power law of [`FrictionLaw::ShiftedPower`], the form of Woodhouse
//!   & Galluzzo (2025) eq. (13); [`RosinCalibration`] names that one with its rig.
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

use crate::tangential::{TangentialLaw, TangentialStep, TangentialStepInput};
use crate::{
    FrictionLaw, FrictionRegime, TriboError, checked_add, checked_mul, finite, nonnegative_finite,
    positive_finite,
};
use core::f64::consts::PI;
use fs_math::det;
use std::sync::Arc;

/// Largest admitted history capacity, in steps.
pub const MAX_HISTORY_STEPS: usize = 1 << 22;
const CHUNK: usize = 1024;
const CERTIFICATE_PIECES: usize = 4096;

pub(crate) fn checked(value: f64, field: &'static str) -> Result<f64, TriboError> {
    finite(value, field)?;
    Ok(value)
}
pub(crate) fn outside(field: &'static str, value: f64, range: [f64; 2]) -> TriboError {
    TriboError::OutsideApplicability {
        field,
        value,
        minimum: range[0],
        maximum: range[1],
    }
}

/// Validated common input; returns the released slip distance [m].
pub(crate) fn admit_input(x: &TangentialStepInput) -> Result<f64, TriboError> {
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

/// Steady-sliding curve `mu(V)` that the calibration inverts through Eq. (20).
/// Both forms are non-increasing in speed, so `mu(w)` at the right end of an
/// interval bounds `mu` below on it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SteadySlidingCurve {
    /// The paper's Eq. (2).
    Exponential(SteadySlidingFit),
    /// `mu(V) = kinetic_mu + (static_mu - kinetic_mu) (1 + V/v0)^(-exponent)`,
    /// the curve of [`FrictionLaw::ShiftedPower`] with the same operations,
    /// so it equals that law bit for bit. Woodhouse & Galluzzo (2025)
    /// eq. (13) is this form.
    ShiftedPower {
        /// Coefficient at zero speed, positive.
        static_mu: f64,
        /// High-speed asymptote, `0 <= kinetic_mu <= static_mu`.
        kinetic_mu: f64,
        /// `v0` [m/s], positive.
        characteristic_speed_mps: f64,
        /// `n`, positive.
        exponent: f64,
    },
}

impl SteadySlidingCurve {
    /// The steady-sliding curve of a memoryless [`FrictionLaw`]:
    /// `ShiftedPower` as itself and `Hyperbolic` as `exponent = 1` (equal to
    /// rounding, since that law has its own quotient form).
    ///
    /// # Errors
    /// An invalid law, or any other law, refused by name: `Coulomb` would make
    /// `k_y` constant, `VelocityDependent` and viscous `Stribeck` curves can
    /// rise with speed, and `Stribeck`'s Gaussian decay is not bounded here.
    pub fn from_friction_law(law: &FrictionLaw) -> Result<Self, TriboError> {
        law.validate()?;
        let curve = match *law {
            FrictionLaw::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                exponent,
            } => Self::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed_mps: characteristic_speed,
                exponent,
            },
            FrictionLaw::Hyperbolic {
                static_mu,
                kinetic_mu,
                characteristic_speed,
            } => Self::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed_mps: characteristic_speed,
                exponent: 1.0,
            },
            FrictionLaw::Coulomb { .. }
            | FrictionLaw::VelocityDependent { .. }
            | FrictionLaw::Stribeck { .. } => {
                return Err(TriboError::InvalidInput {
                    field: "smith_woodhouse_calibration_curve_unsupported_law",
                });
            }
        };
        curve.validate()?;
        Ok(curve)
    }

    fn validate(&self) -> Result<(), TriboError> {
        match *self {
            Self::Exponential(fit) => fit.validate(),
            Self::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed_mps,
                exponent,
            } => {
                positive_finite(static_mu, "smith_woodhouse_curve_static_mu")?;
                nonnegative_finite(kinetic_mu, "smith_woodhouse_curve_kinetic_mu")?;
                if kinetic_mu > static_mu {
                    return Err(TriboError::InvalidInput {
                        field: "smith_woodhouse_curve_kinetic_mu_above_static",
                    });
                }
                positive_finite(characteristic_speed_mps, "smith_woodhouse_curve_speed")?;
                positive_finite(exponent, "smith_woodhouse_curve_exponent")
            }
        }
    }

    /// Steady-sliding coefficient at speed `|v|` [m/s].
    #[must_use]
    pub fn coefficient(&self, v: f64) -> f64 {
        match *self {
            Self::Exponential(fit) => fit.coefficient(v),
            Self::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed_mps,
                exponent,
            } => {
                let decay = det::pow(1.0 + v.abs() / characteristic_speed_mps, -exponent);
                kinetic_mu + (static_mu - kinetic_mu) * decay
            }
        }
    }

    /// An upper bound of `-w mu'(w)` on `[lo, hi]`. Each term is `c g(w/V)`
    /// with `g` unimodal (`x e^{-x}`, peak at 1; `x (1+x)^{-n-1}`, peak at
    /// `1/n`), bounded by its value at the peak or at the nearer end.
    fn weakening_bound(&self, lo: f64, hi: f64) -> f64 {
        let unimodal_max = |g: &dyn Fn(f64) -> f64, a: f64, b: f64, peak: f64| -> f64 {
            if a <= peak && peak <= b {
                g(peak)
            } else if b < peak {
                g(b)
            } else {
                g(a)
            }
        };
        match *self {
            Self::Exponential(fit) => fit
                .amplitudes
                .iter()
                .zip(fit.decay_speeds_mps)
                .map(|(c, v)| c * unimodal_max(&|x| x * det::exp(-x), lo / v, hi / v, 1.0))
                .sum(),
            Self::ShiftedPower {
                static_mu,
                kinetic_mu,
                characteristic_speed_mps: v0,
                exponent: n,
            } => {
                let g = |x: f64| x * det::pow(1.0 + x, -n - 1.0);
                n * (static_mu - kinetic_mu) * unimodal_max(&g, lo / v0, hi / v0, 1.0 / n)
            }
        }
    }
}

/// A named steady-sliding calibration: the curve, the rig it was measured on,
/// and its speed range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RosinCalibration {
    /// Steady-sliding curve.
    pub curve: SteadySlidingCurve,
    /// The rig; its `reference_normal_n` is the rig load.
    pub rig: RosinContact,
    /// Upper end [m/s] of the speed range: the measured range, or, from
    /// [`Self::woodhouse_galluzzo_2025_eq13_authors_extended`], the range the
    /// authors print the fit over, which is not measured above 0.3 m/s.
    pub speed_max_mps: f64,
}

impl RosinCalibration {
    /// Woodhouse & Galluzzo (2025) eq. (13), `mu = 0.2 (|v| + 0.011)^-0.4`
    /// ([`FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13`]): their revised fit
    /// to the steady-sliding measurements of their Fig. 1, which span about
    /// 0.3 m/s. Those are this paper's steady-sliding data: WG2025's eq. (1),
    /// static 1.2 and asymptote 0.35, is this paper's Eq. (2) fit to them.
    ///
    /// So the rig is this paper's steady-sliding rig: contact radius 500 um
    /// at 3 N and layer thickness 10 um (section 5.3; the "steady sliding"
    /// column of Table 2 of Galluzzo, Woodhouse & Mansour 2017, Acta Acustica
    /// 103:1080-1099). Both bounding solids are the caller's `rig_solid`. That
    /// rig's rod and wedge were perspex (Table 1), and Eq. (20) conducts into
    /// them. WG2025's eq. (5) writes the rosin's own `K_r` and `D_r` there
    /// instead; that is not the rig the data were measured on.
    ///
    /// # Errors
    /// Non-positive `rig_solid` properties.
    pub fn woodhouse_galluzzo_2025_eq13(rig_solid: ThermalSolid) -> Result<Self, TriboError> {
        const RIG_RADIUS_M: f64 = 500e-6;
        const RIG_LOAD_N: f64 = 3.0;
        const RIG_LAYER_M: f64 = 10e-6;
        const MEASURED_SPEED_MAX_MPS: f64 = 0.3;
        rig_solid.validate()?;
        Ok(Self {
            curve: SteadySlidingCurve::from_friction_law(
                &FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13(),
            )?,
            rig: RosinContact {
                radius_m: RIG_RADIUS_M,
                reference_normal_n: RIG_LOAD_N,
                layer_thickness_m: RIG_LAYER_M,
                stationary: rig_solid,
                moving: rig_solid,
            },
            speed_max_mps: MEASURED_SPEED_MAX_MPS,
        })
    }

    /// [`Self::woodhouse_galluzzo_2025_eq13`] with its speed range widened
    /// from the measured 0.3 m/s to 5 m/s, the range over which Woodhouse &
    /// Galluzzo (2025) print the curve's contact temperature (their Fig. 7,
    /// sliding speed 0 to 5 m/s). Their simulated contact temperatures go
    /// above the measured range too (their Fig. 13 peaks near 71 C). An
    /// opt-in only.
    ///
    /// Above 0.3 m/s the curve is the authors' extrapolated fit, not a
    /// measurement: no friction or rosin shear strength is printed there.
    /// Smith & Woodhouse (2000, p. 1666) keep their own runs below 40 K above
    /// ambient, "probably low enough for plausibility of this model".
    ///
    /// The curve and the rig are unchanged, so below the measured limit
    /// `k_y(T)` is the same function. On this rig Eq. (20) maps 5 m/s to
    /// about 45.27 K. Here `T_ss` peaks near 6.5 m/s at about 45.33 K: with
    /// `mu ~ w^-0.4`, it falls once the layer's advection reaches about a
    /// quarter of the conduction term. So no speed range reaches the about
    /// 53 K that the authors' Fig. 7 reads at 5 m/s; their eq. (5) contact
    /// conducts into rosin, not into the rig's solid.
    ///
    /// # Errors
    /// Non-positive `rig_solid` properties.
    pub fn woodhouse_galluzzo_2025_eq13_authors_extended(
        rig_solid: ThermalSolid,
    ) -> Result<Self, TriboError> {
        const FIG7_SPEED_MAX_MPS: f64 = 5.0;
        Ok(Self {
            speed_max_mps: FIG7_SPEED_MAX_MPS,
            ..Self::woodhouse_galluzzo_2025_eq13(rig_solid)?
        })
    }
}

/// Every input of [`SmithWoodhouse`]. Nothing is defaulted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmithWoodhouseParams {
    /// Rosin density `rho_r` [kg/m^3], positive.
    pub rosin_density_kg_per_m3: f64,
    /// Rosin effective specific heat `c_r` [J/kg/K], positive.
    pub rosin_specific_heat_j_per_kg_k: f64,
    /// Steady-sliding curve used for calibration: Eq. (2), or a
    /// [`FrictionLaw`] curve such as Woodhouse & Galluzzo (2025) eq. (13).
    pub steady_fit: SteadySlidingCurve,
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
    /// uniform partition the left side is bounded above (each unimodal
    /// weakening term by its interval maximum) and `mu` below (by its value at
    /// the right end).
    fn certify_calibration(&self) -> Result<(), TriboError> {
        let curve = &self.params.steady_fit;
        let max = self.params.calibration_speed_max_mps;
        #[allow(clippy::cast_precision_loss)]
        let pieces = CERTIFICATE_PIECES as f64;
        for i in 0..CERTIFICATE_PIECES {
            #[allow(clippy::cast_precision_loss)]
            let (lo, hi) = (max * i as f64 / pieces, max * (i + 1) as f64 / pieces);
            let slope = curve.weakening_bound(lo, hi);
            let bound = 2.0 * slope * (1.0 + self.alpha / self.gamma * det::sqrt(hi));
            if !(bound < curve.coefficient(hi)) {
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
