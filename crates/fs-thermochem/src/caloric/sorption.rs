//! Moisture sorption isotherms (GAB and Oswin) and Clausius-Clapeyron isosteric heat.
//!
//! Supplies equilibrium moisture content `X_eq(a_w, T)` [kg water / kg dry solid],
//! exact analytical inverse `a_w(X_w, T)`, exact analytical derivatives
//! `(da_w/dX_w)_T` and `(da_w/dT)_{X_w}` for implicit porous-media Newton solves,
//! and Clausius-Clapeyron net and total isosteric heats of sorption `q_st(X_w, T)`.
//!
//! - **GAB (Guggenheim-Anderson-de Boer):**
//!   `X_w(a_w, T) = X_m * C(T) * K(T) * a_w / ((1 - K(T)*a_w) * (1 - K(T)*a_w + C(T)*K(T)*a_w))`
//!   with Arrhenius `C(T) = C_0 * exp(Delta H_c / (R T))`, `K(T) = K_0 * exp(Delta H_k / (R T))`.
//! - **Oswin:**
//!   `X_w(a_w, T) = k(T) * (a_w / (1 - a_w))^n`
//!   with either linear (`k(T) = k_ref + dk_dT * (T - T_ref)`) or Arrhenius
//!   (`k(T) = k_0 * exp(Delta H_k / (R T))`) temperature dependence.
//!
//! All transcendentals route strictly through `fs_math::det` for deterministic
//! evaluation with zero external runtime dependencies.

use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use crate::water::{CRITICAL_TEMPERATURE_K, SaturatedWaterState, TRIPLE_POINT_TEMPERATURE_K};
use core::fmt;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy};
use fs_math::det::{exp, pow, sqrt};
use fs_qty::{Dims, Temperature};

/// Molar energy dimensions `[L^2 M T^-2 N^-1]` (J/mol).
const MOLAR_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, -1]);
/// Inverse-temperature dimensions `[Theta^-1]` (1/K).
const INV_TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, -1, 0, 0]);

/// Molar mass of ordinary water [kg/mol], IUPAC standard atomic weights.
pub const WATER_MOLAR_MASS_KG_PER_MOL_V1: f64 = 0.018_015_28;
/// Exact `fs-matdb` law id for a GAB moisture sorption isotherm card.
pub const GAB_ISOTHERM_LAW_ID_V1: &str = "gab-sorption-isotherm";
/// Exact `fs-matdb` law id for an Oswin moisture sorption isotherm card.
pub const OSWIN_ISOTHERM_LAW_ID_V1: &str = "oswin-sorption-isotherm";
/// Supported semantic version for sorption isotherm `fs-matdb` cards.
pub const SORPTION_LAW_VERSION_V1: u32 = 1;

/// Dedicated typed refusal for moisture sorption isotherm construction and evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum SorptionErrorV1 {
    /// Invalid isotherm parameters (non-finite, non-positive, or out of physical bounds).
    InvalidParameter(&'static str),
    /// Temperature is non-finite, non-positive, or outside the declared validity domain.
    TemperatureOutOfDomain {
        /// Rejected temperature [K].
        temperature_k: f64,
        /// Inclusive lower validity bound [K].
        min_k: f64,
        /// Inclusive upper validity bound [K].
        max_k: f64,
    },
    /// Water activity is non-finite or outside the valid open interval `(0, 1)`.
    WaterActivityOutOfDomain {
        /// Rejected water activity.
        water_activity: f64,
    },
    /// Dry-basis moisture content is non-finite, negative, or above the isotherm's
    /// finite bound-water ceiling at `a_w -> 1^-`.
    MoistureOutOfDomain {
        /// Rejected dry-basis moisture content [kg water / kg dry solid].
        moisture_dry_basis: f64,
    },
    /// GAB denominator `1 - K(T) * a_w` (or second factor) is non-positive or singular.
    SingularGabDenominator {
        /// Evaluated `K(T)`.
        k_value: f64,
        /// Requested water activity `a_w`.
        water_activity: f64,
        /// Evaluated `1 - K(T) * a_w`.
        denominator: f64,
    },
    /// An intermediate or derived quantity overflowed, underflowed, or became non-finite.
    Unrepresentable(&'static str),
}

impl fmt::Display for SorptionErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter(msg) => write!(f, "invalid sorption isotherm parameter: {msg}"),
            Self::TemperatureOutOfDomain {
                temperature_k,
                min_k,
                max_k,
            } => write!(
                f,
                "sorption temperature {temperature_k} K outside domain [{min_k}, {max_k}] K"
            ),
            Self::WaterActivityOutOfDomain { water_activity } => {
                write!(
                    f,
                    "water activity {water_activity} outside admissible domain (0, 1)"
                )
            }
            Self::MoistureOutOfDomain { moisture_dry_basis } => {
                write!(
                    f,
                    "dry-basis moisture {moisture_dry_basis} outside admissible sorption domain"
                )
            }
            Self::SingularGabDenominator {
                k_value,
                water_activity,
                denominator,
            } => write!(
                f,
                "singular GAB denominator 1 - K*a_w = {denominator} (K={k_value}, a_w={water_activity})"
            ),
            Self::Unrepresentable(msg) => write!(f, "unrepresentable sorption quantity: {msg}"),
        }
    }
}

impl core::error::Error for SorptionErrorV1 {}

/// Complete thermodynamic evaluation of a moisture sorption isotherm at `(X_w, a_w, T)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SorptionEvaluationV1 {
    temperature: Temperature,
    water_activity: f64,
    moisture_dry_basis: f64,
    moisture_wet_basis: f64,
    dx_daw: f64,
    daw_dx: f64,
    daw_dt: f64,
    net_isosteric_heat_j_per_mol: f64,
    net_isosteric_heat_j_per_kg: f64,
}

impl SorptionEvaluationV1 {
    /// Absolute temperature [K].
    #[must_use]
    pub const fn temperature(self) -> Temperature {
        self.temperature
    }

    /// Equilibrium water activity `a_w` in `[0, 1)`.
    #[must_use]
    pub const fn water_activity(self) -> f64 {
        self.water_activity
    }

    /// Dry-basis equilibrium moisture content `X_w` [kg water / kg dry solid].
    #[must_use]
    pub const fn moisture_dry_basis(self) -> f64 {
        self.moisture_dry_basis
    }

    /// Wet-basis moisture mass fraction `w_w = X_w / (1 + X_w)` [kg water / kg total].
    #[must_use]
    pub const fn moisture_wet_basis(self) -> f64 {
        self.moisture_wet_basis
    }

    /// Analytical isothermal slope `(dX_w / da_w)_T` [kg water / kg dry solid].
    #[must_use]
    pub const fn dx_daw(self) -> f64 {
        self.dx_daw
    }

    /// Analytical inverse isothermal slope `(da_w / dX_w)_T` [kg dry solid / kg water].
    #[must_use]
    pub const fn daw_dx(self) -> f64 {
        self.daw_dx
    }

    /// Analytical isosteric temperature sensitivity `(da_w / dT)_{X_w}` [1/K].
    #[must_use]
    pub const fn daw_dt(self) -> f64 {
        self.daw_dt
    }

    /// Net isosteric heat of sorption `q_st^net = -R (d ln a_w / d(1/T))_{X_w}` [J/mol].
    #[must_use]
    pub const fn net_isosteric_heat_j_per_mol(self) -> f64 {
        self.net_isosteric_heat_j_per_mol
    }

    /// Net isosteric heat of sorption per unit mass of water [J/kg].
    #[must_use]
    pub const fn net_isosteric_heat_j_per_kg(self) -> f64 {
        self.net_isosteric_heat_j_per_kg
    }
}

fn validate_sorption_temp_domain(
    min_t: Temperature,
    max_t: Temperature,
) -> Result<(f64, f64), SorptionErrorV1> {
    let lo = min_t.value();
    let hi = max_t.value();
    if !lo.is_finite() || lo <= 0.0 || !hi.is_finite() || hi < lo {
        return Err(SorptionErrorV1::InvalidParameter(
            "temperature domain must be finite, strictly positive, and ordered",
        ));
    }
    Ok((lo, hi))
}

fn check_sorption_temp(
    temperature: Temperature,
    min_k: f64,
    max_k: f64,
) -> Result<f64, SorptionErrorV1> {
    let t = temperature.value();
    if !t.is_finite() || t <= 0.0 || t < min_k || t > max_k {
        return Err(SorptionErrorV1::TemperatureOutOfDomain {
            temperature_k: t,
            min_k,
            max_k,
        });
    }
    Ok(t)
}

fn check_open_water_activity(water_activity: f64) -> Result<f64, SorptionErrorV1> {
    if !water_activity.is_finite() || water_activity <= 0.0 || water_activity >= 1.0 {
        return Err(SorptionErrorV1::WaterActivityOutOfDomain { water_activity });
    }
    Ok(water_activity)
}

fn water_latent_heat_j_per_kg(temperature: Temperature) -> Result<f64, SorptionErrorV1> {
    let t = temperature.value();
    let sat = SaturatedWaterState::at_temperature(temperature).map_err(|_| {
        SorptionErrorV1::TemperatureOutOfDomain {
            temperature_k: t,
            min_k: TRIPLE_POINT_TEMPERATURE_K,
            max_k: CRITICAL_TEMPERATURE_K,
        }
    })?;
    let h_fg = sat.latent_enthalpy().value();
    if !h_fg.is_finite() {
        return Err(SorptionErrorV1::Unrepresentable(
            "water latent heat of vaporization is non-finite",
        ));
    }
    Ok(h_fg)
}

/// Guggenheim-Anderson-de Boer (GAB) moisture sorption isotherm with Arrhenius
/// temperature dependence on `C(T)` and `K(T)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GabIsothermV1 {
    monolayer_moisture_xm: f64,
    c0: f64,
    delta_h_c_j_per_mol: f64,
    k0: f64,
    delta_h_k_j_per_mol: f64,
    temperature_min: Temperature,
    temperature_max: Temperature,
}

impl GabIsothermV1 {
    /// Construct a validated GAB moisture sorption isotherm.
    ///
    /// - `monolayer_moisture_xm`: `X_m > 0` [kg water / kg dry solid]
    /// - `c0`: Arrhenius pre-exponential `C_0 > 0` (dimensionless)
    /// - `delta_h_c_j_per_mol`: `Delta H_c = H_1 - H_m` [J/mol] (finite)
    /// - `k0`: Arrhenius pre-exponential `K_0 > 0` (dimensionless)
    /// - `delta_h_k_j_per_mol`: `Delta H_k = H_m - H_l` [J/mol] (finite)
    /// - `minimum_temperature..=maximum_temperature`: inclusive positive domain [K]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        monolayer_moisture_xm: f64,
        c0: f64,
        delta_h_c_j_per_mol: f64,
        k0: f64,
        delta_h_k_j_per_mol: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, SorptionErrorV1> {
        if !monolayer_moisture_xm.is_finite() || monolayer_moisture_xm <= 0.0 {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB monolayer moisture X_m must be finite and strictly positive",
            ));
        }
        if !c0.is_finite() || c0 <= 0.0 {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB C_0 must be finite and strictly positive",
            ));
        }
        if !k0.is_finite() || k0 <= 0.0 {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB K_0 must be finite and strictly positive",
            ));
        }
        if !delta_h_c_j_per_mol.is_finite() || !delta_h_k_j_per_mol.is_finite() {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB enthalpy differences must be finite",
            ));
        }
        validate_sorption_temp_domain(minimum_temperature, maximum_temperature)?;
        Ok(Self {
            monolayer_moisture_xm,
            c0,
            delta_h_c_j_per_mol,
            k0,
            delta_h_k_j_per_mol,
            temperature_min: minimum_temperature,
            temperature_max: maximum_temperature,
        })
    }

    /// Load a validated GAB isotherm from an immutable `fs-matdb` constitutive card.
    ///
    /// Expected parameters: `x_m`, `c0`, `delta_h_c`, `k0`, `delta_h_k`, and validity `"T"`.
    pub fn from_card(card: &ConstitutiveModelCard) -> Result<Self, SorptionErrorV1> {
        card.validate()
            .map_err(|_| SorptionErrorV1::InvalidParameter("GAB matdb card failed validation"))?;
        if card.law.0 != GAB_ISOTHERM_LAW_ID_V1 {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB matdb card law id mismatch",
            ));
        }
        if card.law_version != SORPTION_LAW_VERSION_V1 {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB matdb card law version mismatch",
            ));
        }
        if card.state_schema_version != 0
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(SorptionErrorV1::InvalidParameter(
                "GAB matdb card must be stateless",
            ));
        }
        let (t_min, t_max) = card
            .validity
            .bound("T")
            .ok_or(SorptionErrorV1::InvalidParameter(
                "GAB matdb card missing T validity axis",
            ))?;
        let param = |name: &'static str, expected_dims: Dims| -> Result<f64, SorptionErrorV1> {
            let p = card
                .parameters
                .get(name)
                .ok_or(SorptionErrorV1::InvalidParameter(
                    "GAB matdb card missing required parameter",
                ))?;
            if p.dims != expected_dims {
                return Err(SorptionErrorV1::InvalidParameter(
                    "GAB matdb card parameter dimension mismatch",
                ));
            }
            Ok(p.value)
        };
        Self::new(
            param("x_m", Dims::NONE)?,
            param("c0", Dims::NONE)?,
            param("delta_h_c", MOLAR_ENERGY_DIMS)?,
            param("k0", Dims::NONE)?,
            param("delta_h_k", MOLAR_ENERGY_DIMS)?,
            Temperature::new(t_min),
            Temperature::new(t_max),
        )
    }

    /// Monolayer moisture content `X_m` [kg water / kg dry solid].
    #[must_use]
    pub const fn monolayer_moisture_xm(self) -> f64 {
        self.monolayer_moisture_xm
    }

    /// Inclusive valid temperature domain `[T_min, T_max]` [K].
    #[must_use]
    pub const fn temperature_domain(self) -> [Temperature; 2] {
        [self.temperature_min, self.temperature_max]
    }

    /// Evaluate Guggenheim constant `C(T) = C_0 * exp(Delta H_c / (R T))`.
    pub fn c_at(self, temperature: Temperature) -> Result<f64, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let c = self.c0 * exp(self.delta_h_c_j_per_mol / (R * t));
        if !c.is_finite() || c <= 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB C(T) overflowed or underflowed",
            ));
        }
        Ok(c)
    }

    /// Evaluate multilayer factor `K(T) = K_0 * exp(Delta H_k / (R T))`.
    pub fn k_at(self, temperature: Temperature) -> Result<f64, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let k = self.k0 * exp(self.delta_h_k_j_per_mol / (R * t));
        if !k.is_finite() || k <= 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB K(T) overflowed or underflowed",
            ));
        }
        Ok(k)
    }

    /// Forward GAB isotherm `X_eq(a_w, T)` [kg water / kg dry solid].
    ///
    /// Refuses `a_w <= 0`, `a_w >= 1`, out-of-domain temperatures, and singular
    /// GAB denominators `1 - K(T) * a_w <= 0`.
    pub fn equilibrium_moisture(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_activity(water_activity, temperature)?
            .moisture_dry_basis)
    }

    /// Exact analytical inverse GAB isotherm `a_w(X_w, T)` on the unique monotonic
    /// physical branch `a_w in [0, 1)`.
    ///
    /// Refuses negative/non-finite `X_w` and moisture contents exceeding the finite
    /// bound-water sorption ceiling when `K(T) < 1`.
    pub fn water_activity(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .water_activity)
    }

    /// Full evaluation of the GAB isotherm and analytical derivatives starting from `(a_w, T)`.
    pub fn evaluate_from_activity(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let aw = check_open_water_activity(water_activity)?;
        let c = self.c_at(temperature)?;
        let k = self.k_at(temperature)?;

        let u = k * aw;
        let d1 = 1.0 - u;
        let d2 = d1 + c * u;
        if !d1.is_finite() || d1 <= 0.0 || !d2.is_finite() || d2 <= 0.0 {
            return Err(SorptionErrorV1::SingularGabDenominator {
                k_value: k,
                water_activity: aw,
                denominator: d1,
            });
        }

        let x_w = (self.monolayer_moisture_xm * c * u) / (d1 * d2);
        if !x_w.is_finite() || x_w < 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB equilibrium moisture is non-finite",
            ));
        }

        self.finish_evaluation(temperature, t, aw, x_w, c, k, u, d1, d2)
    }

    /// Full evaluation of the GAB isotherm and analytical derivatives starting from `(X_w, T)`.
    pub fn evaluate_from_moisture(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        if !moisture_dry_basis.is_finite() || moisture_dry_basis < 0.0 {
            return Err(SorptionErrorV1::MoistureOutOfDomain { moisture_dry_basis });
        }
        let c = self.c_at(temperature)?;
        let k = self.k_at(temperature)?;

        if moisture_dry_basis == 0.0 {
            let dx_daw = self.monolayer_moisture_xm * c * k;
            if !dx_daw.is_finite() || dx_daw <= 0.0 {
                return Err(SorptionErrorV1::Unrepresentable(
                    "GAB zero-moisture slope is unrepresentable",
                ));
            }
            let q_st_net = self.delta_h_k_j_per_mol + self.delta_h_c_j_per_mol;
            let q_st_kg = q_st_net / WATER_MOLAR_MASS_KG_PER_MOL_V1;
            if !q_st_net.is_finite() || !q_st_kg.is_finite() {
                return Err(SorptionErrorV1::Unrepresentable(
                    "GAB zero-moisture isosteric heat is non-finite",
                ));
            }
            return Ok(SorptionEvaluationV1 {
                temperature,
                water_activity: 0.0,
                moisture_dry_basis: 0.0,
                moisture_wet_basis: 0.0,
                dx_daw,
                daw_dx: 1.0 / dx_daw,
                daw_dt: 0.0,
                net_isosteric_heat_j_per_mol: q_st_net,
                net_isosteric_heat_j_per_kg: q_st_kg,
            });
        }

        // Solve X_w*(C - 1)*u^2 + [X_m*C - X_w*(C - 2)]*u - X_w = 0 for unique u in (0, 1).
        // Exact nonnegative discriminant: Delta = C^2*(X_m - X_w)^2 + 4*C*X_m*X_w.
        let xm = self.monolayer_moisture_xm;
        let xw = moisture_dry_basis;
        let diff = xm - xw;
        let disc = (c * diff) * (c * diff) + 4.0 * c * xm * xw;
        if !disc.is_finite() || disc < 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB inverse discriminant overflowed",
            ));
        }
        let sqrt_disc = sqrt(disc);
        let b = xm * c - xw * (c - 2.0);
        let u = if b >= 0.0 {
            // Cancellation-free Citardauq form (also handles C = 1 without division by 0).
            (2.0 * xw) / (b + sqrt_disc)
        } else {
            // b < 0 implies C > 2, so a = xw * (C - 1) > xw > 0 and -b > 0.
            let a = xw * (c - 1.0);
            (-b + sqrt_disc) / (2.0 * a)
        };
        if !u.is_finite() || u <= 0.0 || u >= 1.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB inverse root u outside (0, 1)",
            ));
        }
        let aw = u / k;
        if !aw.is_finite() || aw <= 0.0 || aw >= 1.0 {
            // Moisture exceeds the finite GAB saturation ceiling X_w(a_w -> 1^-, T) when K(T) < 1.
            return Err(SorptionErrorV1::MoistureOutOfDomain { moisture_dry_basis });
        }
        let d1 = 1.0 - u;
        let d2 = d1 + c * u;
        if !d1.is_finite() || d1 <= 0.0 || !d2.is_finite() || d2 <= 0.0 {
            return Err(SorptionErrorV1::SingularGabDenominator {
                k_value: k,
                water_activity: aw,
                denominator: d1,
            });
        }

        self.finish_evaluation(temperature, t, aw, xw, c, k, u, d1, d2)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_evaluation(
        self,
        temperature: Temperature,
        t: f64,
        aw: f64,
        x_w: f64,
        c: f64,
        k: f64,
        u: f64,
        d1: f64,
        d2: f64,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        // Numerator factor 1 + (C - 1)*u^2 = (1 - u^2) + C*u^2 > 0 strictly for u in (0, 1), C > 0.
        let curv_factor = (1.0 - u * u) + c * (u * u);
        let denom_sq = (d1 * d2) * (d1 * d2);
        if !curv_factor.is_finite()
            || curv_factor <= 0.0
            || !denom_sq.is_finite()
            || denom_sq <= 0.0
        {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB derivative denominator overflowed or underflowed",
            ));
        }
        let dx_daw = (self.monolayer_moisture_xm * c * k * curv_factor) / denom_sq;
        let daw_dx = denom_sq / (self.monolayer_moisture_xm * c * k * curv_factor);

        // Exact Clausius-Clapeyron net isosteric heat of sorption:
        // q_st^net = -R * (d ln a_w / d(1/T))_{X_w}
        //          = Delta H_k + [(1 - u)^2 / (1 + (C - 1)*u^2)] * Delta H_c
        let q_st_net =
            self.delta_h_k_j_per_mol + ((d1 * d1) / curv_factor) * self.delta_h_c_j_per_mol;
        let q_st_kg = q_st_net / WATER_MOLAR_MASS_KG_PER_MOL_V1;
        let daw_dt = aw * q_st_net / (R * t * t);
        let wet_basis = x_w / (1.0 + x_w);

        if ![dx_daw, daw_dx, q_st_net, q_st_kg, daw_dt, wet_basis]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(SorptionErrorV1::Unrepresentable(
                "GAB derived sorption quantities became non-finite",
            ));
        }

        Ok(SorptionEvaluationV1 {
            temperature,
            water_activity: aw,
            moisture_dry_basis: x_w,
            moisture_wet_basis: wet_basis,
            dx_daw,
            daw_dx,
            daw_dt,
            net_isosteric_heat_j_per_mol: q_st_net,
            net_isosteric_heat_j_per_kg: q_st_kg,
        })
    }

    /// Analytical derivative `(da_w / dX_w)_T` [kg dry solid / kg water].
    pub fn daw_dx(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .daw_dx)
    }

    /// Analytical derivative `(da_w / dT)_{X_w}` [1/K].
    pub fn daw_dt(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .daw_dt)
    }

    /// Clausius-Clapeyron net isosteric heat of sorption `q_st^net(X_w, T)` [J/mol].
    pub fn net_isosteric_heat_j_per_mol(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .net_isosteric_heat_j_per_mol)
    }

    /// Total isosteric heat of sorption `Q_st = q_st^net + h_fg(T)` [J/kg water],
    /// combining the analytical GAB Clausius-Clapeyron excess heat with IAPWS-SR1-86
    /// pure-water latent heat of vaporization.
    pub fn total_isosteric_heat_j_per_kg(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        let eval = self.evaluate_from_moisture(moisture_dry_basis, temperature)?;
        let h_fg = water_latent_heat_j_per_kg(temperature)?;
        let total = eval.net_isosteric_heat_j_per_kg + h_fg;
        if !total.is_finite() {
            return Err(SorptionErrorV1::Unrepresentable(
                "total isosteric heat is non-finite",
            ));
        }
        Ok(total)
    }

    /// Representative culinary preset for native potato starch (van den Berg 1984 /
    /// Al-Muhtaseb et al. 2004 ranges: `X_m = 0.095`, `C(25 C) ≈ 17.5`, `K(25 C) ≈ 0.78`).
    #[must_use]
    pub fn potato_starch_v1() -> Self {
        Self::new(
            0.095,
            0.015,
            17_500.0,
            0.60,
            650.0,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("potato starch GAB preset is valid")
    }

    /// Representative culinary preset for white bread crumb (Lomauro et al. 1985 /
    /// Bell & Labuza 2000 ranges: `X_m = 0.065`, `C(25 C) ≈ 12.0`, `K(25 C) ≈ 0.85`).
    #[must_use]
    pub fn bread_crumb_v1() -> Self {
        Self::new(
            0.065,
            0.025,
            15_300.0,
            0.65,
            660.0,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("bread crumb GAB preset is valid")
    }

    /// Representative culinary preset for lean beef muscle (Iglesias & Chirife /
    /// Trujillo et al. 2003 ranges: `X_m = 0.068`, `C(25 C) ≈ 10.0`, `K(25 C) ≈ 0.87`).
    #[must_use]
    pub fn lean_beef_v1() -> Self {
        Self::new(
            0.068,
            0.030,
            14_400.0,
            0.70,
            540.0,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("lean beef GAB preset is valid")
    }

    /// Representative culinary preset for potato tissue (Datta/Farkas drying ranges:
    /// `X_m = 0.078`, `C(25 C) ≈ 14.5`, `K(25 C) ≈ 0.82`).
    #[must_use]
    pub fn potato_tissue_v1() -> Self {
        Self::new(
            0.078,
            0.020,
            16_300.0,
            0.63,
            650.0,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("potato tissue GAB preset is valid")
    }
}

/// Temperature dependence of the Oswin scale parameter `k(T)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OswinTemperatureLawV1 {
    /// Modified Oswin linear form (ASAE D245.5): `k(T) = k_ref + dk_dt_per_k * (T - T_ref)`.
    Linear {
        /// Scale parameter `k(T_ref) > 0` [kg water / kg dry solid].
        k_ref: f64,
        /// Linear temperature slope `dk/dT` [1/K] (typically negative for foods).
        dk_dt_per_k: f64,
        /// Reference temperature `T_ref` [K].
        reference_temperature_k: f64,
    },
    /// Arrhenius form: `k(T) = k0 * exp(enthalpy_j_per_mol / (R T))`.
    Arrhenius {
        /// Pre-exponential factor `k_0 > 0` [kg water / kg dry solid].
        k0: f64,
        /// Effective sorption enthalpy `Delta H_k` [J/mol] (positive => `k` decreases with `T`).
        enthalpy_j_per_mol: f64,
    },
}

/// Oswin moisture sorption isotherm `X_w(a_w, T) = k(T) * (a_w / (1 - a_w))^n`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OswinIsothermV1 {
    temperature_law: OswinTemperatureLawV1,
    exponent_n: f64,
    temperature_min: Temperature,
    temperature_max: Temperature,
}

impl OswinIsothermV1 {
    /// Construct a modified Oswin isotherm with linear `k(T) = k_ref + dk_dt * (T - T_ref)`.
    ///
    /// Validates that `k(T) > 0` at both endpoints of `[minimum_temperature, maximum_temperature]`
    /// (and therefore everywhere on the interval by linearity).
    pub fn new_linear(
        k_ref: f64,
        dk_dt_per_k: f64,
        reference_temperature: Temperature,
        exponent_n: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, SorptionErrorV1> {
        let (t_min, t_max) =
            validate_sorption_temp_domain(minimum_temperature, maximum_temperature)?;
        let t_ref = reference_temperature.value();
        if !k_ref.is_finite()
            || k_ref <= 0.0
            || !dk_dt_per_k.is_finite()
            || !t_ref.is_finite()
            || t_ref <= 0.0
            || !exponent_n.is_finite()
            || exponent_n <= 0.0
        {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin linear parameters and exponent n must be finite and valid",
            ));
        }
        let k_lo = k_ref + dk_dt_per_k * (t_min - t_ref);
        let k_hi = k_ref + dk_dt_per_k * (t_max - t_ref);
        if !k_lo.is_finite() || k_lo <= 0.0 || !k_hi.is_finite() || k_hi <= 0.0 {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin linear k(T) must remain strictly positive across the temperature domain",
            ));
        }
        Ok(Self {
            temperature_law: OswinTemperatureLawV1::Linear {
                k_ref,
                dk_dt_per_k,
                reference_temperature_k: t_ref,
            },
            exponent_n,
            temperature_min: minimum_temperature,
            temperature_max: maximum_temperature,
        })
    }

    /// Construct an Oswin isotherm with Arrhenius `k(T) = k0 * exp(enthalpy / (R T))`.
    pub fn new_arrhenius(
        k0: f64,
        enthalpy_j_per_mol: f64,
        exponent_n: f64,
        minimum_temperature: Temperature,
        maximum_temperature: Temperature,
    ) -> Result<Self, SorptionErrorV1> {
        validate_sorption_temp_domain(minimum_temperature, maximum_temperature)?;
        if !k0.is_finite()
            || k0 <= 0.0
            || !enthalpy_j_per_mol.is_finite()
            || !exponent_n.is_finite()
            || exponent_n <= 0.0
        {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin Arrhenius parameters and exponent n must be finite and positive",
            ));
        }
        Ok(Self {
            temperature_law: OswinTemperatureLawV1::Arrhenius {
                k0,
                enthalpy_j_per_mol,
            },
            exponent_n,
            temperature_min: minimum_temperature,
            temperature_max: maximum_temperature,
        })
    }

    /// Load a validated Oswin isotherm from an immutable `fs-matdb` constitutive card.
    ///
    /// Supports either linear parameters (`k_ref`, `dk_dt`, `t_ref`, `n`) or Arrhenius
    /// parameters (`k0`, `enthalpy`, `n`), plus validity axis `"T"`.
    pub fn from_card(card: &ConstitutiveModelCard) -> Result<Self, SorptionErrorV1> {
        card.validate()
            .map_err(|_| SorptionErrorV1::InvalidParameter("Oswin matdb card failed validation"))?;
        if card.law.0 != OSWIN_ISOTHERM_LAW_ID_V1 {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card law id mismatch",
            ));
        }
        if card.law_version != SORPTION_LAW_VERSION_V1 {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card law version mismatch",
            ));
        }
        if card.state_schema_version != 0
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card must be stateless",
            ));
        }
        let (t_min, t_max) = card
            .validity
            .bound("T")
            .ok_or(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card missing T validity axis",
            ))?;
        let n_param = card
            .parameters
            .get("n")
            .ok_or(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card missing exponent n",
            ))?;
        if n_param.dims != Dims::NONE {
            return Err(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card exponent n must be dimensionless",
            ));
        }
        let n = n_param.value;
        if let (Some(k_ref), Some(dk_dt), Some(t_ref)) = (
            card.parameters.get("k_ref"),
            card.parameters.get("dk_dt"),
            card.parameters.get("t_ref"),
        ) {
            if k_ref.dims != Dims::NONE
                || dk_dt.dims != INV_TEMPERATURE_DIMS
                || t_ref.dims != Temperature::DIMS
            {
                return Err(SorptionErrorV1::InvalidParameter(
                    "Oswin linear matdb card parameter dimension mismatch",
                ));
            }
            Self::new_linear(
                k_ref.value,
                dk_dt.value,
                Temperature::new(t_ref.value),
                n,
                Temperature::new(t_min),
                Temperature::new(t_max),
            )
        } else if let (Some(k0), Some(enthalpy)) =
            (card.parameters.get("k0"), card.parameters.get("enthalpy"))
        {
            if k0.dims != Dims::NONE || enthalpy.dims != MOLAR_ENERGY_DIMS {
                return Err(SorptionErrorV1::InvalidParameter(
                    "Oswin Arrhenius matdb card parameter dimension mismatch",
                ));
            }
            Self::new_arrhenius(
                k0.value,
                enthalpy.value,
                n,
                Temperature::new(t_min),
                Temperature::new(t_max),
            )
        } else {
            Err(SorptionErrorV1::InvalidParameter(
                "Oswin matdb card requires either (k_ref, dk_dt, t_ref) or (k0, enthalpy)",
            ))
        }
    }

    /// Oswin shape exponent `n > 0`.
    #[must_use]
    pub const fn exponent_n(self) -> f64 {
        self.exponent_n
    }

    /// Temperature law for `k(T)`.
    #[must_use]
    pub const fn temperature_law(self) -> OswinTemperatureLawV1 {
        self.temperature_law
    }

    /// Inclusive valid temperature domain `[T_min, T_max]` [K].
    #[must_use]
    pub const fn temperature_domain(self) -> [Temperature; 2] {
        [self.temperature_min, self.temperature_max]
    }

    /// Evaluate `(k(T), dk/dT)` at an admitted temperature.
    pub fn k_and_derivative_at(
        self,
        temperature: Temperature,
    ) -> Result<(f64, f64), SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        match self.temperature_law {
            OswinTemperatureLawV1::Linear {
                k_ref,
                dk_dt_per_k,
                reference_temperature_k,
            } => {
                let k = k_ref + dk_dt_per_k * (t - reference_temperature_k);
                if !k.is_finite() || k <= 0.0 {
                    return Err(SorptionErrorV1::Unrepresentable(
                        "Oswin k(T) became non-positive or non-finite",
                    ));
                }
                Ok((k, dk_dt_per_k))
            }
            OswinTemperatureLawV1::Arrhenius {
                k0,
                enthalpy_j_per_mol,
            } => {
                let k = k0 * exp(enthalpy_j_per_mol / (R * t));
                let dk_dt = -k * enthalpy_j_per_mol / (R * t * t);
                if !k.is_finite() || k <= 0.0 || !dk_dt.is_finite() {
                    return Err(SorptionErrorV1::Unrepresentable(
                        "Oswin Arrhenius k(T) overflowed or underflowed",
                    ));
                }
                Ok((k, dk_dt))
            }
        }
    }

    /// Forward Oswin isotherm `X_eq(a_w, T) = k(T) * (a_w / (1 - a_w))^n`.
    pub fn equilibrium_moisture(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_activity(water_activity, temperature)?
            .moisture_dry_basis)
    }

    /// Analytical inverse Oswin isotherm `a_w(X_w, T)`.
    ///
    /// Returns the continuous dry limit `0.0` at `X_w = 0.0` for all `n > 0`.
    pub fn water_activity(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        if moisture_dry_basis == 0.0 {
            let _ = self.k_and_derivative_at(temperature)?;
            return Ok(0.0);
        }
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .water_activity)
    }

    /// Full evaluation starting from `(a_w, T)`.
    pub fn evaluate_from_activity(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        let aw = check_open_water_activity(water_activity)?;
        let (k, dk_dt) = self.k_and_derivative_at(temperature)?;
        let odds = aw / (1.0 - aw);
        let x_w = k * pow(odds, self.exponent_n);
        if !x_w.is_finite() || x_w <= 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "Oswin equilibrium moisture overflowed or underflowed",
            ));
        }
        self.finish_evaluation(temperature, t, aw, x_w, k, dk_dt)
    }

    /// Full evaluation starting from `(X_w, T)`.
    pub fn evaluate_from_moisture(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        let t = check_sorption_temp(
            temperature,
            self.temperature_min.value(),
            self.temperature_max.value(),
        )?;
        if !moisture_dry_basis.is_finite() || moisture_dry_basis < 0.0 {
            return Err(SorptionErrorV1::MoistureOutOfDomain { moisture_dry_basis });
        }
        let (k, dk_dt) = self.k_and_derivative_at(temperature)?;

        if moisture_dry_basis == 0.0 {
            let (dx_daw, daw_dx) = if (self.exponent_n - 1.0).abs() <= f64::EPSILON {
                (k, 1.0 / k)
            } else if self.exponent_n < 1.0 {
                // X ~ k * a_w^n with n < 1 => dX/da_w -> +inf, da_w/dX = 0 at X = 0.
                return Err(SorptionErrorV1::Unrepresentable(
                    "Oswin slope dX/da_w diverges at X_w = 0 for n < 1",
                ));
            } else {
                // n > 1 => dX/da_w = 0, da_w/dX -> +inf at X = 0.
                return Err(SorptionErrorV1::Unrepresentable(
                    "Oswin inverse slope da_w/dX diverges at X_w = 0 for n > 1",
                ));
            };
            let q_st_net = -(R * t * t / self.exponent_n) * (dk_dt / k);
            let q_st_kg = q_st_net / WATER_MOLAR_MASS_KG_PER_MOL_V1;
            if !q_st_net.is_finite() || !q_st_kg.is_finite() {
                return Err(SorptionErrorV1::Unrepresentable(
                    "Oswin zero-moisture isosteric heat is non-finite",
                ));
            }
            return Ok(SorptionEvaluationV1 {
                temperature,
                water_activity: 0.0,
                moisture_dry_basis: 0.0,
                moisture_wet_basis: 0.0,
                dx_daw,
                daw_dx,
                daw_dt: 0.0,
                net_isosteric_heat_j_per_mol: q_st_net,
                net_isosteric_heat_j_per_kg: q_st_kg,
            });
        }

        let ratio = moisture_dry_basis / k;
        let odds = pow(ratio, 1.0 / self.exponent_n);
        if !odds.is_finite() || odds <= 0.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "Oswin inverse odds ratio overflowed or underflowed",
            ));
        }
        let aw = odds / (1.0 + odds);
        if !aw.is_finite() || aw <= 0.0 || aw >= 1.0 {
            return Err(SorptionErrorV1::Unrepresentable(
                "Oswin inverse water activity rounded to boundary 0 or 1",
            ));
        }
        self.finish_evaluation(temperature, t, aw, moisture_dry_basis, k, dk_dt)
    }

    fn finish_evaluation(
        self,
        temperature: Temperature,
        t: f64,
        aw: f64,
        x_w: f64,
        k: f64,
        dk_dt: f64,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        let one_minus_aw = 1.0 - aw;
        let dx_daw = (self.exponent_n * x_w) / (aw * one_minus_aw);
        let daw_dx = (aw * one_minus_aw) / (self.exponent_n * x_w);
        let dln_k_dt = dk_dt / k;
        let q_st_net = -(R * t * t * one_minus_aw / self.exponent_n) * dln_k_dt;
        let q_st_kg = q_st_net / WATER_MOLAR_MASS_KG_PER_MOL_V1;
        let daw_dt = -(aw * one_minus_aw / self.exponent_n) * dln_k_dt;
        let wet_basis = x_w / (1.0 + x_w);

        if ![dx_daw, daw_dx, q_st_net, q_st_kg, daw_dt, wet_basis]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(SorptionErrorV1::Unrepresentable(
                "Oswin derived sorption quantities became non-finite",
            ));
        }

        Ok(SorptionEvaluationV1 {
            temperature,
            water_activity: aw,
            moisture_dry_basis: x_w,
            moisture_wet_basis: wet_basis,
            dx_daw,
            daw_dx,
            daw_dt,
            net_isosteric_heat_j_per_mol: q_st_net,
            net_isosteric_heat_j_per_kg: q_st_kg,
        })
    }

    /// Analytical derivative `(da_w / dX_w)_T` [kg dry solid / kg water].
    pub fn daw_dx(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .daw_dx)
    }

    /// Analytical derivative `(da_w / dT)_{X_w}` [1/K].
    pub fn daw_dt(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .daw_dt)
    }

    /// Clausius-Clapeyron net isosteric heat of sorption `q_st^net(X_w, T)` [J/mol].
    pub fn net_isosteric_heat_j_per_mol(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        Ok(self
            .evaluate_from_moisture(moisture_dry_basis, temperature)?
            .net_isosteric_heat_j_per_mol)
    }

    /// Total isosteric heat of sorption `Q_st = q_st^net + h_fg(T)` [J/kg water].
    pub fn total_isosteric_heat_j_per_kg(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        let eval = self.evaluate_from_moisture(moisture_dry_basis, temperature)?;
        let h_fg = water_latent_heat_j_per_kg(temperature)?;
        let total = eval.net_isosteric_heat_j_per_kg + h_fg;
        if !total.is_finite() {
            return Err(SorptionErrorV1::Unrepresentable(
                "total isosteric heat is non-finite",
            ));
        }
        Ok(total)
    }

    /// Representative culinary modified-Oswin preset for bread (`k(25 C) = 0.105`, `n = 0.43`).
    #[must_use]
    pub fn bread_oswin_v1() -> Self {
        Self::new_linear(
            0.105,
            -4.5e-4,
            Temperature::new(298.15),
            0.43,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("bread Oswin preset is valid")
    }

    /// Representative culinary modified-Oswin preset for lean meat (`k(25 C) = 0.110`, `n = 0.46`).
    #[must_use]
    pub fn meat_oswin_v1() -> Self {
        Self::new_linear(
            0.110,
            -4.2e-4,
            Temperature::new(298.15),
            0.46,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("meat Oswin preset is valid")
    }

    /// Representative culinary Arrhenius-Oswin preset for starch (`k(25 C) ≈ 0.145`, `n = 0.38`).
    #[must_use]
    pub fn starch_oswin_v1() -> Self {
        Self::new_arrhenius(
            0.032,
            3_750.0,
            0.38,
            Temperature::new(273.15),
            Temperature::new(423.15),
        )
        .expect("starch Oswin preset is valid")
    }
}

/// Unified sorption isotherm dispatch for downstream porous-media solvers (`fs-porous`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SorptionIsothermV1 {
    /// Guggenheim-Anderson-de Boer (GAB) isotherm.
    Gab(GabIsothermV1),
    /// Oswin (linear or Arrhenius) isotherm.
    Oswin(OswinIsothermV1),
}

impl SorptionIsothermV1 {
    /// Forward equilibrium moisture `X_eq(a_w, T)` [kg water / kg dry solid].
    pub fn equilibrium_moisture(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        match self {
            Self::Gab(g) => g.equilibrium_moisture(water_activity, temperature),
            Self::Oswin(o) => o.equilibrium_moisture(water_activity, temperature),
        }
    }

    /// Inverse water activity `a_w(X_w, T)` in `[0, 1)`.
    pub fn water_activity(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        match self {
            Self::Gab(g) => g.water_activity(moisture_dry_basis, temperature),
            Self::Oswin(o) => o.water_activity(moisture_dry_basis, temperature),
        }
    }

    /// Full evaluation from `(a_w, T)`.
    pub fn evaluate_from_activity(
        self,
        water_activity: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        match self {
            Self::Gab(g) => g.evaluate_from_activity(water_activity, temperature),
            Self::Oswin(o) => o.evaluate_from_activity(water_activity, temperature),
        }
    }

    /// Full evaluation from `(X_w, T)`.
    pub fn evaluate_from_moisture(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<SorptionEvaluationV1, SorptionErrorV1> {
        match self {
            Self::Gab(g) => g.evaluate_from_moisture(moisture_dry_basis, temperature),
            Self::Oswin(o) => o.evaluate_from_moisture(moisture_dry_basis, temperature),
        }
    }

    /// Total isosteric heat of sorption `Q_st = q_st^net + h_fg(T)` [J/kg water].
    pub fn total_isosteric_heat_j_per_kg(
        self,
        moisture_dry_basis: f64,
        temperature: Temperature,
    ) -> Result<f64, SorptionErrorV1> {
        match self {
            Self::Gab(g) => g.total_isosteric_heat_j_per_kg(moisture_dry_basis, temperature),
            Self::Oswin(o) => o.total_isosteric_heat_j_per_kg(moisture_dry_basis, temperature),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::det::ln;

    #[test]
    fn gab_roundtrip_and_analytical_derivatives_hold_across_all_c_regimes() {
        // Test C > 1 (normal food), C == 1 (degenerate linear quadratic), and 0 < C < 1.
        let models = [
            GabIsothermV1::potato_starch_v1(),
            GabIsothermV1::new(
                0.08,
                1.0,
                0.0,
                0.85,
                0.0,
                Temperature::new(273.15),
                Temperature::new(400.0),
            )
            .unwrap(),
            GabIsothermV1::new(
                0.08,
                0.45,
                2_000.0,
                0.75,
                300.0,
                Temperature::new(273.15),
                Temperature::new(400.0),
            )
            .unwrap(),
        ];

        for model in models {
            let t = Temperature::new(310.0);
            for aw in [1e-4, 0.05, 0.25, 0.55, 0.85, 0.995] {
                let fwd = model.evaluate_from_activity(aw, t).unwrap();
                let inv = model
                    .evaluate_from_moisture(fwd.moisture_dry_basis(), t)
                    .unwrap();
                assert!(
                    (inv.water_activity() - aw).abs() < 1e-12,
                    "roundtrip failed for aw={aw}: got {}",
                    inv.water_activity()
                );

                // Central finite difference for da_w/dX_w.
                let x = fwd.moisture_dry_basis();
                let dx = (x * 1e-6).max(1e-10);
                let aw_p = model.water_activity(x + dx, t).unwrap();
                let aw_m = model.water_activity(x - dx, t).unwrap();
                let fd_dx = (aw_p - aw_m) / (2.0 * dx);
                assert!(
                    (inv.daw_dx() - fd_dx).abs() / inv.daw_dx().abs().max(1.0) < 1e-8,
                    "da_w/dX_w mismatch at aw={aw}: analytic={}, fd={fd_dx}",
                    inv.daw_dx()
                );

                // Central finite difference for Clausius-Clapeyron q_st^net.
                let dt = 1e-3;
                let aw_tp = model
                    .water_activity(x, Temperature::new(t.value() + dt))
                    .unwrap();
                let aw_tm = model
                    .water_activity(x, Temperature::new(t.value() - dt))
                    .unwrap();
                let q_fd = R * t.value() * t.value() * (ln(aw_tp) - ln(aw_tm)) / (2.0 * dt);
                assert!(
                    (inv.net_isosteric_heat_j_per_mol() - q_fd).abs()
                        / inv.net_isosteric_heat_j_per_mol().abs().max(1.0)
                        < 1e-8,
                    "q_st mismatch at aw={aw}: analytic={}, fd={q_fd}",
                    inv.net_isosteric_heat_j_per_mol()
                );
            }
        }
    }

    #[test]
    fn oswin_roundtrip_and_derivatives_match_finite_differences() {
        for model in [
            OswinIsothermV1::bread_oswin_v1(),
            OswinIsothermV1::meat_oswin_v1(),
            OswinIsothermV1::starch_oswin_v1(),
        ] {
            let t = Temperature::new(325.0);
            for aw in [1e-4, 0.1, 0.4, 0.75, 0.995] {
                let fwd = model.evaluate_from_activity(aw, t).unwrap();
                let inv = model
                    .evaluate_from_moisture(fwd.moisture_dry_basis(), t)
                    .unwrap();
                assert!(
                    (inv.water_activity() - aw).abs() < 1e-12,
                    "Oswin roundtrip aw={aw}: got {}",
                    inv.water_activity()
                );

                let x = fwd.moisture_dry_basis();
                let dt = 1e-3;
                let aw_tp = model
                    .water_activity(x, Temperature::new(t.value() + dt))
                    .unwrap();
                let aw_tm = model
                    .water_activity(x, Temperature::new(t.value() - dt))
                    .unwrap();
                let q_fd = R * t.value() * t.value() * (ln(aw_tp) - ln(aw_tm)) / (2.0 * dt);
                assert!(
                    (inv.net_isosteric_heat_j_per_mol() - q_fd).abs()
                        / inv.net_isosteric_heat_j_per_mol().abs().max(1.0)
                        < 1e-8,
                    "Oswin q_st mismatch: analytic={}, fd={q_fd}",
                    inv.net_isosteric_heat_j_per_mol()
                );
            }
        }
    }

    #[test]
    fn gab_refuses_singular_denominator_and_out_of_domain_inputs() {
        // K(T) = 1.25 > 1 => singular at a_w >= 0.80.
        let singular_model = GabIsothermV1::new(
            0.08,
            10.0,
            0.0,
            1.25,
            0.0,
            Temperature::new(273.15),
            Temperature::new(400.0),
        )
        .unwrap();
        let t = Temperature::new(300.0);
        assert!(matches!(
            singular_model.equilibrium_moisture(0.80, t),
            Err(SorptionErrorV1::SingularGabDenominator { .. })
        ));
        assert!(matches!(
            singular_model.equilibrium_moisture(0.90, t),
            Err(SorptionErrorV1::SingularGabDenominator { .. })
        ));
        assert!(matches!(
            singular_model.equilibrium_moisture(0.0, t),
            Err(SorptionErrorV1::WaterActivityOutOfDomain { .. })
        ));
        assert!(matches!(
            singular_model.equilibrium_moisture(1.0, t),
            Err(SorptionErrorV1::WaterActivityOutOfDomain { .. })
        ));
    }
}
