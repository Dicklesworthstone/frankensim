//! Choi–Okos (1986) / ASHRAE Ch. 19 proximate thermophysical evaluator and
//! Kirchhoff-consistent regularized freezing-point enthalpy closure.
//!
//! Synthesizes density `rho(T)`, sensible and apparent specific heat
//! `c_p(T)`, `c_{p,app}(T) = dh/dT`, thermal conductivity `k(T)` (parallel,
//! series, Krischer, geometric, Kopelman isotropic, and Maxwell–Eucken),
//! dynamic ice fraction `X_i(T)`, and specific enthalpy `h(T)` on both the
//! local `0 C` Choi–Okos datum and the shared `janaf-298.15K-formation`
//! absolute enthalpy datum from a 6-component USDA proximate mass-fraction
//! vector `(X_w, X_p, X_f, X_c, X_fb, X_a)`.
//!
//! # Physics & Round-2 Invariants
//! - **Empirical scope (no Gibbs-potential overclaim):** Choi–Okos quadratic
//!   polynomials are empirical correlations valid over `[-40 C, +150 C]`.
//!   Above `+150 C` (up to `+250 C`) evaluations return
//!   [`ValidityStatus::AboveDesiccationWarning`] so callers transition dry
//!   crusts to charred/porous models. Outside `[-50 C, +250 C]` evaluation
//!   refuses with [`ProximateErrorV1::TemperatureOutOfRange`].
//! - **Piecewise `C^1`-regularized Chen–Heldman ice fraction:** Strictly
//!   `X_i(T) = 0` for `T >= T_f`, exact Celsius-offset Chen–Heldman
//!   `(X_w - X_b) * (1 - (T_f - 273.15) / (T - 273.15))` for
//!   `T <= T_f - delta_T`, and a monotone, overshoot-free cubic Hermite mushy
//!   transition across `[T_f - delta_T, T_f]` (`delta_T` defaulting to
//!   `0.05 K`). Because `T_f <= 273.15 K` and `delta_T > 0`, pure water
//!   (`T_f = 273.15 K`) never divides `0 / 0 = NaN` or explodes `dX_i/dT`.
//! - **Automatic Kirchhoff consistency:** Total specific enthalpy is computed
//!   directly from phase enthalpies `sum_alpha X_alpha(T) h_alpha(T)` with
//!   `h_ice(0 C) = -L_{f,0}`, guaranteeing `dL_f/dT = c_{p,w}(T) - c_{p,ice}(T)`
//!   and exact `integral_{T_1}^{T_2} c_{p,app}(T) dT = h(T_2) - h(T_1)` across
//!   freezing without constant-`L_f` bias at `-40 C`.
//! - **Unified `janaf-298.15K-formation` stream datum:** Both food enthalpy
//!   and IAPWS-IF97 water/steam enthalpy expose exact formation-shifted
//!   accessors anchored to liquid `H2O` standard formation enthalpy
//!   `h_f0(298.15 K, 1 bar) = -285.830 kJ/mol = -15.866 MJ/kg`.
//!
//! # Primary Sources
//! - Y. Choi and M. R. Okos (1986), *Effects of Temperature and Composition on
//!   the Thermal Properties of Foods*, in *Food Engineering and Process
//!   Applications*, Vol. 1, pp. 93–101, Elsevier Applied Science.
//! - *ASHRAE Handbook — Refrigeration*, Chapter 19 (Thermal Properties of
//!   Foods, Table 1, including the ASHRAE erratum correction `-2.7604e-4` for
//!   the linear fat thermal conductivity term).
//! - J. K. Carson, S. J. Lovatt, D. J. Tanner, and A. C. Cleland (2005/2006),
//!   effective thermal conductivity bounds and Maxwell–Eucken / Kopelman /
//!   Krischer structural models for porous foods.
//! - CODATA / NIST-JANAF 298.15 K liquid water formation enthalpy
//!   (`-285.830 kJ/mol`) and IAPWS R7-97(2012) molar mass (`18.015268 g/mol`).

use crate::water::conductivity::ThermalConductivity;
use crate::water::if97::{If97Region, If97State};
use crate::{
    MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use core::fmt;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawParameter};
use fs_math::det;
use fs_qty::{Density, Dims, Pressure, Temperature};
use std::collections::BTreeMap;

/// Deterministic evaluator version for Choi–Okos proximate thermophysics.
pub const CHOI_OKOS_EVALUATOR_VERSION_V1: u32 = 1;
/// Exact `fs-matdb` law id for a Choi–Okos proximate coefficient table card.
pub const CHOI_OKOS_LAW_ID_V1: &str = "choi-okos-proximate";
/// Supported semantic law version for [`CHOI_OKOS_LAW_ID_V1`] cards.
pub const CHOI_OKOS_LAW_VERSION_V1: u32 = 1;
/// Maximum allowed `|sum(X_j) - 1.0|` tolerance when admitting proximate fractions.
pub const PROXIMATE_COMPOSITION_SUM_TOLERANCE_V1: f64 = 1e-4;
/// Normal freezing point of pure water at 1 atm [K] (`0.0 C`).
pub const ICE_POINT_K: f64 = 273.15;
/// Standard JANAF reference temperature [K] (`25.0 C`).
pub const JANAF_REFERENCE_TEMPERATURE_K: f64 = 298.15;
/// Standard JANAF reference pressure [Pa] (`1 bar`).
pub const JANAF_REFERENCE_PRESSURE_PA: f64 = 100_000.0;
/// Canonical identifier for the shared coupled-stream absolute enthalpy datum.
pub const JANAF_FORMATION_DATUM_ID_V1: &str = "janaf-298.15K-formation";
/// IAPWS-IF97 / IAPWS-95 molar mass of ordinary water `H2O` [kg/mol].
pub const H2O_MOLAR_MASS_KG_PER_MOL: f64 = 0.018_015_268;
/// Standard molar enthalpy of formation of liquid `H2O(l)` at `298.15 K, 1 bar` [J/mol]
/// (CODATA 1989 / NIST-JANAF).
pub const H2O_LIQUID_MOLAR_FORMATION_ENTHALPY_298K_J_PER_MOL: f64 = -285_830.0;
/// Mass-specific standard formation enthalpy of liquid `H2O(l)` at `298.15 K, 1 bar` [J/kg]
/// (`-15.865987 MJ/kg`).
pub const H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG: f64 =
    H2O_LIQUID_MOLAR_FORMATION_ENTHALPY_298K_J_PER_MOL / H2O_MOLAR_MASS_KG_PER_MOL;
/// Standard latent heat of fusion of pure water/ice at `0 C` [J/kg] (ASHRAE Ch. 19).
pub const ICE_FUSION_ENTHALPY_0C_J_KG: f64 = 333_600.0;
/// Default `C^1` mushy-zone regularization width `delta_T` [K] across `[T_f - delta_T, T_f]`.
pub const DEFAULT_MUSHY_ZONE_WIDTH_K: f64 = 0.05;
/// Lower inclusive Choi–Okos calibration bound [K] (`-40 C`).
pub const MINIMUM_VALID_TEMPERATURE_K: f64 = 233.15;
/// Upper inclusive Choi–Okos calibration bound [K] (`+150 C`, desiccation onset above).
pub const MAXIMUM_VALID_TEMPERATURE_K: f64 = 423.15;
/// Lower hard admitted temperature bound [K] (`-50 C`); temperatures below refuse.
pub const MINIMUM_ADMITTED_TEMPERATURE_K: f64 = 223.15;
/// Upper hard admitted temperature bound [K] (`+250 C`); temperatures above refuse.
pub const MAXIMUM_ADMITTED_TEMPERATURE_K: f64 = 523.15;
/// ASHRAE Ch. 19 Table 7 enthalpy zero-reference temperature [K] (`-40 C`).
pub const ASHRAE_ENTHALPY_REFERENCE_TEMPERATURE_K: f64 = 233.15;
/// Dry-air specific gas constant [J/(kg K)] used for optional gas-pore void density.
pub const DRY_AIR_GAS_CONSTANT_J_KG_K: f64 = 287.058;
/// Standard sea-level atmospheric pressure [Pa] used for default gas-pore void density.
pub const STANDARD_ATMOSPHERE_PA: f64 = 101_325.0;
/// Provenance classification tag for the built-in ASHRAE / Choi–Okos empirical table.
pub const ASHRAE_CHOI_OKOS_SOURCE_TAG: &str = "empirical";
/// Bibliographic provenance citation for the built-in ASHRAE / Choi–Okos table.
pub const ASHRAE_CHOI_OKOS_CITATION: &str = "Choi & Okos (1986), Food Engineering and Process Applications Vol. 1, pp. 93-101; ASHRAE Handbook - Refrigeration Ch. 19 (Thermal Properties of Foods, Table 1, with ASHRAE fat thermal conductivity erratum -2.7604e-4 W/(m K^2))";

const CP_C0_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CP_C1_DIMS: Dims = Dims([2, 0, -2, -2, 0, 0]);
const CP_C2_DIMS: Dims = Dims([2, 0, -2, -3, 0, 0]);
const RHO_C0_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const RHO_C1_DIMS: Dims = Dims([-3, 1, 0, -1, 0, 0]);
const RHO_C2_DIMS: Dims = Dims([-3, 1, 0, -2, 0, 0]);
const K_C0_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const K_C1_DIMS: Dims = Dims([1, 1, -3, -2, 0, 0]);
const K_C2_DIMS: Dims = Dims([1, 1, -3, -3, 0, 0]);
const SPECIFIC_ENERGY_DIMS: Dims = Dims([2, 0, -2, 0, 0, 0]);

/// Dedicated typed refusal enum for Choi–Okos proximate thermophysical construction
/// and evaluation (isolated from [`crate::ThermochemErrorV1`]).
#[derive(Debug, Clone, PartialEq)]
pub enum ProximateErrorV1 {
    /// A proximate mass fraction or bound-water fraction is non-finite or negative.
    InvalidFraction {
        /// Offending field name.
        field: &'static str,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
        /// Stable diagnostic reason.
        reason: &'static str,
    },
    /// The six proximate mass fractions do not sum to 1 within the admitted tolerance.
    CompositionSumOutOfTolerance {
        /// Observed raw sum of the six fractions.
        sum: f64,
        /// Maximum admitted `|sum - 1.0|` tolerance.
        tolerance: f64,
    },
    /// Unfrozen bound-water mass fraction `x_b` exceeds total water mass fraction `water`.
    BoundWaterExceedsWater {
        /// Offered bound-water fraction.
        x_b: f64,
        /// Offered total water mass fraction.
        water: f64,
    },
    /// Initial freezing point `t_f_k` is non-finite, `<= 0 K`, or above `273.15 K`.
    InvalidFreezingPoint {
        /// Exact rejected IEEE-754 bits.
        t_f_k_bits: u64,
        /// Stable diagnostic reason.
        reason: &'static str,
    },
    /// Gas-pore void fraction is non-finite or outside `[0.0, 1.0)`.
    InvalidPorosity {
        /// Exact rejected IEEE-754 bits.
        porosity_bits: u64,
    },
    /// Mushy-zone smoothing width `delta_t_k` is non-finite or not strictly positive.
    InvalidMushyZoneWidth {
        /// Exact rejected IEEE-754 bits.
        delta_t_k_bits: u64,
    },
    /// Conductivity mixing model parameter (such as Krischer `anisotropy_factor`) is invalid.
    InvalidConductivityParameter {
        /// Parameter name.
        parameter: &'static str,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
        /// Stable diagnostic reason.
        reason: &'static str,
    },
    /// Evaluation temperature is non-finite or `<= 0 K`.
    InvalidTemperature {
        /// Exact rejected IEEE-754 bits.
        temperature_k_bits: u64,
    },
    /// Positive finite temperature lies outside the admitted evaluator envelope.
    TemperatureOutOfRange {
        /// Requested temperature [K].
        temperature_k: f64,
        /// Inclusive lower admitted temperature [K].
        minimum_k: f64,
        /// Inclusive upper admitted temperature [K].
        maximum_k: f64,
    },
    /// Coefficient table metadata, polynomial coefficients, or evaluated pure-component
    /// properties violated finiteness or positivity invariants.
    InvalidCoefficientTable {
        /// Component name.
        component: &'static str,
        /// Property or field name.
        property: &'static str,
        /// Stable diagnostic reason.
        reason: &'static str,
    },
    /// An `fs-matdb` constitutive model card failed validation, law, dimension, or axis checks.
    MaterialCard {
        /// Human-readable refusal diagnostic.
        reason: String,
    },
    /// An IAPWS-IF97 reference state or enthalpy argument was invalid.
    InvalidIf97Reference {
        /// Stable diagnostic reason.
        reason: &'static str,
    },
    /// Fixed-order evaluation produced a non-finite or physically inadmissible value.
    NonFiniteEvaluation {
        /// Derived property name.
        property: &'static str,
        /// Exact non-finite or invalid IEEE-754 bits.
        bits: u64,
    },
}

impl fmt::Display for ProximateErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidFraction {
                field,
                bits,
                reason,
            } => write!(
                f,
                "proximate fraction {field:?} refused (bits {bits:#018x}): {reason}"
            ),
            Self::CompositionSumOutOfTolerance { sum, tolerance } => write!(
                f,
                "proximate mass fractions sum to {sum:.8}, outside 1.0 +/- {tolerance:.1e}"
            ),
            Self::BoundWaterExceedsWater { x_b, water } => write!(
                f,
                "bound water fraction x_b={x_b} exceeds total water fraction water={water}"
            ),
            Self::InvalidFreezingPoint { t_f_k_bits, reason } => write!(
                f,
                "initial freezing point refused (bits {t_f_k_bits:#018x}): {reason}"
            ),
            Self::InvalidPorosity { porosity_bits } => write!(
                f,
                "gas porosity must be finite and in [0, 1) (bits {porosity_bits:#018x})"
            ),
            Self::InvalidMushyZoneWidth { delta_t_k_bits } => write!(
                f,
                "mushy-zone width delta_T must be finite and strictly positive (bits {delta_t_k_bits:#018x})"
            ),
            Self::InvalidConductivityParameter {
                parameter,
                bits,
                reason,
            } => write!(
                f,
                "conductivity parameter {parameter:?} refused (bits {bits:#018x}): {reason}"
            ),
            Self::InvalidTemperature { temperature_k_bits } => write!(
                f,
                "temperature must be finite and strictly positive (bits {temperature_k_bits:#018x})"
            ),
            Self::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            } => write!(
                f,
                "temperature {temperature_k:.3} K is outside admitted Choi-Okos range [{minimum_k:.2}, {maximum_k:.2}] K"
            ),
            Self::InvalidCoefficientTable {
                component,
                property,
                reason,
            } => write!(
                f,
                "Choi-Okos table invalid for {component}.{property}: {reason}"
            ),
            Self::MaterialCard { reason } => {
                write!(f, "Choi-Okos constitutive model card refused: {reason}")
            }
            Self::InvalidIf97Reference { reason } => {
                write!(f, "IAPWS-IF97 JANAF datum reference refused: {reason}")
            }
            Self::NonFiniteEvaluation { property, bits } => write!(
                f,
                "Choi-Okos {property:?} evaluation became non-finite or non-positive (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for ProximateErrorV1 {}

/// Explicit physical validity classification for a Choi–Okos temperature point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValidityStatus {
    /// Inside the published Choi–Okos / ASHRAE calibration range `[-40 C, +150 C]`.
    Valid,
    /// Between `-50 C` and `-40 C`: guarded cryogenic extension below published calibration.
    BelowCalibrationWarning,
    /// Strictly above `+150 C` up to `+250 C`: surface crust desiccation warning zone;
    /// callers should transition dry crusts to charred/porous conductivity laws.
    AboveDesiccationWarning,
    /// Outside `[-50 C, +250 C]` or non-positive/non-finite: refused by [`ChoiOkosEvaluator::evaluate`].
    RefusedOutOfRange,
}

/// Classify a temperature against the standard ASHRAE / Choi–Okos validity bands.
#[must_use]
pub fn classify_temperature_validity(temperature: Temperature) -> ValidityStatus {
    classify_temperature_with_bounds(
        temperature.value(),
        MINIMUM_ADMITTED_TEMPERATURE_K,
        MINIMUM_VALID_TEMPERATURE_K,
        MAXIMUM_VALID_TEMPERATURE_K,
        MAXIMUM_ADMITTED_TEMPERATURE_K,
    )
}

fn classify_temperature_with_bounds(
    t_k: f64,
    admitted_min_k: f64,
    cal_min_k: f64,
    cal_max_k: f64,
    admitted_max_k: f64,
) -> ValidityStatus {
    if !t_k.is_finite() || t_k < admitted_min_k || t_k > admitted_max_k {
        ValidityStatus::RefusedOutOfRange
    } else if t_k < cal_min_k {
        ValidityStatus::BelowCalibrationWarning
    } else if t_k > cal_max_k {
        ValidityStatus::AboveDesiccationWarning
    } else {
        ValidityStatus::Valid
    }
}

/// Raw 6-component proximate mass fractions on a condensed food mass basis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProximateFractions {
    /// Total water mass fraction `X_w` (unfrozen liquid + ice).
    pub water: f64,
    /// Protein mass fraction `X_p`.
    pub protein: f64,
    /// Fat / lipid mass fraction `X_f`.
    pub fat: f64,
    /// Carbohydrate mass fraction `X_c`.
    pub carbohydrate: f64,
    /// Dietary fiber mass fraction `X_fb`.
    pub fiber: f64,
    /// Ash / mineral mass fraction `X_a`.
    pub ash: f64,
}

/// Validated 6-component proximate food composition with bound water `x_b`,
/// initial freezing point `t_f_k`, and optional gas-pore void fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProximateComposition {
    water: f64,
    protein: f64,
    fat: f64,
    carbohydrate: f64,
    fiber: f64,
    ash: f64,
    x_b: f64,
    t_f_k: f64,
    porosity_gas: f64,
}

impl ProximateComposition {
    /// Construct and validate a proximate food composition with explicit bound
    /// water `x_b` and initial freezing point `t_f_k` [K].
    ///
    /// All six mass fractions must be finite and non-negative, with
    /// `|sum - 1.0| <= PROXIMATE_COMPOSITION_SUM_TOLERANCE_V1` (`1e-4`).
    /// Small USDA table rounding residuals within `1e-4` are normalized to
    /// exact unit sum after validating `0 <= x_b <= water`.
    pub fn new(
        fractions: ProximateFractions,
        x_b: f64,
        t_f_k: f64,
    ) -> Result<Self, ProximateErrorV1> {
        let checks = [
            ("water", fractions.water),
            ("protein", fractions.protein),
            ("fat", fractions.fat),
            ("carbohydrate", fractions.carbohydrate),
            ("fiber", fractions.fiber),
            ("ash", fractions.ash),
        ];
        for (field, value) in checks {
            if !value.is_finite() {
                return Err(ProximateErrorV1::InvalidFraction {
                    field,
                    bits: value.to_bits(),
                    reason: "fraction must be finite",
                });
            }
            if value < 0.0 {
                return Err(ProximateErrorV1::InvalidFraction {
                    field,
                    bits: value.to_bits(),
                    reason: "fraction must be non-negative",
                });
            }
        }
        let sum = fractions.water
            + fractions.protein
            + fractions.fat
            + fractions.carbohydrate
            + fractions.fiber
            + fractions.ash;
        if !sum.is_finite() || (sum - 1.0).abs() > PROXIMATE_COMPOSITION_SUM_TOLERANCE_V1 {
            return Err(ProximateErrorV1::CompositionSumOutOfTolerance {
                sum,
                tolerance: PROXIMATE_COMPOSITION_SUM_TOLERANCE_V1,
            });
        }
        if !x_b.is_finite() || x_b < 0.0 {
            return Err(ProximateErrorV1::InvalidFraction {
                field: "x_b",
                bits: x_b.to_bits(),
                reason: "bound-water fraction must be finite and non-negative",
            });
        }
        if x_b > fractions.water {
            return Err(ProximateErrorV1::BoundWaterExceedsWater {
                x_b,
                water: fractions.water,
            });
        }
        validate_freezing_point(t_f_k)?;

        let (water, protein, fat, carbohydrate, fiber, ash, x_b_norm) =
            if sum.to_bits() == 1.0_f64.to_bits() {
                (
                    fractions.water,
                    fractions.protein,
                    fractions.fat,
                    fractions.carbohydrate,
                    fractions.fiber,
                    fractions.ash,
                    x_b,
                )
            } else {
                let inv = 1.0 / sum;
                let w = fractions.water * inv;
                (
                    w,
                    fractions.protein * inv,
                    fractions.fat * inv,
                    fractions.carbohydrate * inv,
                    fractions.fiber * inv,
                    fractions.ash * inv,
                    (x_b * inv).min(w),
                )
            };

        Ok(Self {
            water,
            protein,
            fat,
            carbohydrate,
            fiber,
            ash,
            x_b: x_b_norm,
            t_f_k,
            porosity_gas: 0.0,
        })
    }

    /// Pure liquid/solid water (`water = 1.0`, `x_b = 0.0`, `t_f_k = 273.15 K`).
    #[must_use]
    pub const fn pure_water() -> Self {
        Self {
            water: 1.0,
            protein: 0.0,
            fat: 0.0,
            carbohydrate: 0.0,
            fiber: 0.0,
            ash: 0.0,
            x_b: 0.0,
            t_f_k: ICE_POINT_K,
            porosity_gas: 0.0,
        }
    }

    /// Construct from proximate fractions using ASHRAE Ch. 19 / Schwartzberg
    /// bound-water estimation (`x_b = min(b_protein * protein, water)`) and
    /// Chen–Heldman / Raoult freezing-point depression from soluble solids.
    ///
    /// Typical values: `bound_water_per_protein = 0.4` (Schwartzberg 1976),
    /// `effective_solute_molar_mass_kg_mol` in `[0.06, 0.35]` kg/mol.
    pub fn from_fractions_with_raoult_freezing(
        fractions: ProximateFractions,
        bound_water_per_protein: f64,
        effective_solute_molar_mass_kg_mol: f64,
    ) -> Result<Self, ProximateErrorV1> {
        if !bound_water_per_protein.is_finite() || bound_water_per_protein < 0.0 {
            return Err(ProximateErrorV1::InvalidFraction {
                field: "bound_water_per_protein",
                bits: bound_water_per_protein.to_bits(),
                reason: "bound-water coefficient must be finite and non-negative",
            });
        }
        // Validate fractions first by admitting with provisional freezing point.
        let base = Self::new(fractions, 0.0, ICE_POINT_K)?;
        let x_b = (bound_water_per_protein * base.protein).min(base.water);
        // Soluble solids: carbohydrates + ash (proteins and fats contribute negligibly to molarity).
        let x_solute = base.carbohydrate + base.ash;
        let t_f = estimate_freezing_point_raoult(
            base.water,
            x_b,
            x_solute,
            effective_solute_molar_mass_kg_mol,
        )?;
        Self::new(
            ProximateFractions {
                water: base.water,
                protein: base.protein,
                fat: base.fat,
                carbohydrate: base.carbohydrate,
                fiber: base.fiber,
                ash: base.ash,
            },
            x_b,
            t_f.value(),
        )
    }

    /// Attach an explicit gas-pore void fraction `porosity_gas` in `[0.0, 1.0)`
    /// for aerated dough, bread crumb, fried crust, or overrun ice cream.
    pub fn with_porosity_gas(mut self, porosity_gas: f64) -> Result<Self, ProximateErrorV1> {
        if !porosity_gas.is_finite() || !(0.0..1.0).contains(&porosity_gas) {
            return Err(ProximateErrorV1::InvalidPorosity {
                porosity_bits: porosity_gas.to_bits(),
            });
        }
        self.porosity_gas = porosity_gas;
        Ok(self)
    }

    /// Override the unfrozen bound-water mass fraction `x_b` in `[0.0, water]`.
    pub fn with_bound_water(mut self, x_b: f64) -> Result<Self, ProximateErrorV1> {
        if !x_b.is_finite() || x_b < 0.0 {
            return Err(ProximateErrorV1::InvalidFraction {
                field: "x_b",
                bits: x_b.to_bits(),
                reason: "bound-water fraction must be finite and non-negative",
            });
        }
        if x_b > self.water {
            return Err(ProximateErrorV1::BoundWaterExceedsWater {
                x_b,
                water: self.water,
            });
        }
        self.x_b = x_b;
        Ok(self)
    }

    /// Override the initial freezing point `t_f` [K], requiring `0 < t_f <= 273.15 K`.
    pub fn with_freezing_point(mut self, t_f: Temperature) -> Result<Self, ProximateErrorV1> {
        validate_freezing_point(t_f.value())?;
        self.t_f_k = t_f.value();
        Ok(self)
    }

    /// Total water mass fraction `X_w`.
    #[must_use]
    pub const fn water(&self) -> f64 {
        self.water
    }

    /// Protein mass fraction `X_p`.
    #[must_use]
    pub const fn protein(&self) -> f64 {
        self.protein
    }

    /// Fat / lipid mass fraction `X_f`.
    #[must_use]
    pub const fn fat(&self) -> f64 {
        self.fat
    }

    /// Carbohydrate mass fraction `X_c`.
    #[must_use]
    pub const fn carbohydrate(&self) -> f64 {
        self.carbohydrate
    }

    /// Dietary fiber mass fraction `X_fb`.
    #[must_use]
    pub const fn fiber(&self) -> f64 {
        self.fiber
    }

    /// Ash / mineral mass fraction `X_a`.
    #[must_use]
    pub const fn ash(&self) -> f64 {
        self.ash
    }

    /// Unfrozen bound-water mass fraction `x_b`.
    #[must_use]
    pub const fn bound_water(&self) -> f64 {
        self.x_b
    }

    /// Maximum freezable water mass fraction `X_w - x_b`.
    #[must_use]
    pub fn freezable_water(&self) -> f64 {
        (self.water - self.x_b).max(0.0)
    }

    /// Initial freezing point [K].
    #[must_use]
    pub const fn initial_freezing_point(&self) -> Temperature {
        Temperature::new(self.t_f_k)
    }

    /// Raw initial freezing point scalar [K].
    #[must_use]
    pub const fn initial_freezing_point_k(&self) -> f64 {
        self.t_f_k
    }

    /// Gas-pore void fraction `phi_g` in `[0, 1)`.
    #[must_use]
    pub const fn porosity_gas(&self) -> f64 {
        self.porosity_gas
    }

    /// Total dry solids mass fraction `1 - X_w`.
    #[must_use]
    pub fn solids(&self) -> f64 {
        self.protein + self.fat + self.carbohydrate + self.fiber + self.ash
    }
}

fn validate_freezing_point(t_f_k: f64) -> Result<(), ProximateErrorV1> {
    if !t_f_k.is_finite() {
        return Err(ProximateErrorV1::InvalidFreezingPoint {
            t_f_k_bits: t_f_k.to_bits(),
            reason: "initial freezing point must be finite",
        });
    }
    if t_f_k <= 0.0 {
        return Err(ProximateErrorV1::InvalidFreezingPoint {
            t_f_k_bits: t_f_k.to_bits(),
            reason: "initial freezing point must be strictly positive in kelvin",
        });
    }
    if t_f_k > ICE_POINT_K {
        return Err(ProximateErrorV1::InvalidFreezingPoint {
            t_f_k_bits: t_f_k.to_bits(),
            reason: "initial freezing point cannot exceed pure water freezing point (273.15 K)",
        });
    }
    Ok(())
}

/// Estimate initial freezing point `T_f` [K] via Chen–Heldman / Raoult
/// freezing-point depression:
/// `1 / T_f = 1 / T_0 - (R / (L_{f,0} M_w)) * ln(x_w)`.
///
/// When `x_solute == 0.0` or no freezable water is present, returns `273.15 K`.
pub fn estimate_freezing_point_raoult(
    water: f64,
    x_b: f64,
    solute_mass_fraction: f64,
    effective_solute_molar_mass_kg_mol: f64,
) -> Result<Temperature, ProximateErrorV1> {
    if !water.is_finite() || !(0.0..=1.0).contains(&water) {
        return Err(ProximateErrorV1::InvalidFraction {
            field: "water",
            bits: water.to_bits(),
            reason: "water mass fraction must be finite and in [0, 1]",
        });
    }
    if !x_b.is_finite() || x_b < 0.0 {
        return Err(ProximateErrorV1::InvalidFraction {
            field: "x_b",
            bits: x_b.to_bits(),
            reason: "bound-water fraction must be finite and non-negative",
        });
    }
    if x_b > water {
        return Err(ProximateErrorV1::BoundWaterExceedsWater { x_b, water });
    }
    if !solute_mass_fraction.is_finite() || !(0.0..=1.0).contains(&solute_mass_fraction) {
        return Err(ProximateErrorV1::InvalidFraction {
            field: "solute_mass_fraction",
            bits: solute_mass_fraction.to_bits(),
            reason: "solute mass fraction must be finite and in [0, 1]",
        });
    }
    if !effective_solute_molar_mass_kg_mol.is_finite() || effective_solute_molar_mass_kg_mol <= 0.0
    {
        return Err(ProximateErrorV1::InvalidFraction {
            field: "effective_solute_molar_mass_kg_mol",
            bits: effective_solute_molar_mass_kg_mol.to_bits(),
            reason: "effective solute molar mass must be finite and strictly positive",
        });
    }
    let free_water = water - x_b;
    if free_water == 0.0 || solute_mass_fraction == 0.0 {
        return Ok(Temperature::new(ICE_POINT_K));
    }
    let n_water = free_water / H2O_MOLAR_MASS_KG_PER_MOL;
    let n_solute = solute_mass_fraction / effective_solute_molar_mass_kg_mol;
    let x_water_mol = n_water / (n_water + n_solute);
    let inv_t_f = (1.0 / ICE_POINT_K)
        - (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            / (ICE_FUSION_ENTHALPY_0C_J_KG * H2O_MOLAR_MASS_KG_PER_MOL))
            * det::ln(x_water_mol);
    let t_f_k = 1.0 / inv_t_f;
    if !t_f_k.is_finite() || t_f_k <= 0.0 || t_f_k > ICE_POINT_K {
        return Err(ProximateErrorV1::InvalidFreezingPoint {
            t_f_k_bits: t_f_k.to_bits(),
            reason: "Raoult freezing-point estimate became non-physical",
        });
    }
    Ok(Temperature::new(t_f_k))
}

/// Quadratic Celsius-temperature polynomial `c0 + c1 * t_c + c2 * t_c^2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticCelsiusPoly {
    /// Constant term at `0 C`.
    pub c0: f64,
    /// Linear temperature coefficient (per `K`).
    pub c1: f64,
    /// Quadratic temperature coefficient (per `K^2`).
    pub c2: f64,
}

impl QuadraticCelsiusPoly {
    /// Construct a quadratic polynomial `c0 + c1 * t_c + c2 * t_c^2`.
    #[must_use]
    pub const fn new(c0: f64, c1: f64, c2: f64) -> Self {
        Self { c0, c1, c2 }
    }

    /// Evaluate `c0 + t_c * (c1 + t_c * c2)` using deterministic IEEE-754 arithmetic.
    #[must_use]
    pub fn eval(self, t_c: f64) -> f64 {
        self.c0 + t_c * (self.c1 + t_c * self.c2)
    }

    /// Exact analytical derivative `c1 + 2 * c2 * t_c`.
    #[must_use]
    pub fn derivative(self, t_c: f64) -> f64 {
        self.c1 + 2.0 * self.c2 * t_c
    }

    /// Exact analytical integral from `0 C` to `t_c`:
    /// `t_c * (c0 + t_c * (c1 / 2 + t_c * c2 / 3))`.
    #[must_use]
    pub fn integrate_from_zero(self, t_c: f64) -> f64 {
        t_c * (self.c0 + t_c * (0.5 * self.c1 + (self.c2 / 3.0) * t_c))
    }

    fn is_finite(self) -> bool {
        self.c0.is_finite() && self.c1.is_finite() && self.c2.is_finite()
    }
}

/// Pure-component quadratic property curves for one dry food constituent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComponentPropertyCurves {
    /// Specific heat `c_p(t_c)` [J/(kg K)].
    pub cp: QuadraticCelsiusPoly,
    /// Density `rho(t_c)` [kg/m^3].
    pub density: QuadraticCelsiusPoly,
    /// Thermal conductivity `k(t_c)` [W/(m K)].
    pub conductivity: QuadraticCelsiusPoly,
}

impl ComponentPropertyCurves {
    fn validate(self, component: &'static str) -> Result<(), ProximateErrorV1> {
        if !self.cp.is_finite() {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component,
                property: "cp",
                reason: "polynomial coefficients must be finite",
            });
        }
        if !self.density.is_finite() {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component,
                property: "density",
                reason: "polynomial coefficients must be finite",
            });
        }
        if !self.conductivity.is_finite() {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component,
                property: "conductivity",
                reason: "polynomial coefficients must be finite",
            });
        }
        Ok(())
    }
}

/// Complete Choi–Okos (1986) / ASHRAE Ch. 19 proximate polynomial table.
///
/// Built-in default [`ChoiOkosTable::ashrae_1986`] carries `"empirical"`
/// provenance and bibliographic citation; custom tables may also be loaded
/// from immutable `fs-matdb` [`ConstitutiveModelCard`] records via
/// [`ChoiOkosTable::from_card`].
#[derive(Debug, Clone, PartialEq)]
pub struct ChoiOkosTable {
    /// Provenance classification tag (e.g. `"empirical"`).
    pub provenance_tag: String,
    /// Bibliographic source citation.
    pub citation: String,
    /// Liquid water specific heat for `t_c >= 0 C` [J/(kg K)].
    pub water_cp_above_zero: QuadraticCelsiusPoly,
    /// Supercooled / unfrozen water specific heat for `t_c < 0 C` [J/(kg K)].
    pub water_cp_below_zero: QuadraticCelsiusPoly,
    /// Liquid water density [kg/m^3].
    pub water_density: QuadraticCelsiusPoly,
    /// Liquid water thermal conductivity [W/(m K)].
    pub water_conductivity: QuadraticCelsiusPoly,
    /// Pure ice specific heat [J/(kg K)].
    pub ice_cp: QuadraticCelsiusPoly,
    /// Pure ice density [kg/m^3].
    pub ice_density: QuadraticCelsiusPoly,
    /// Pure ice thermal conductivity [W/(m K)].
    pub ice_conductivity: QuadraticCelsiusPoly,
    /// Protein property curves.
    pub protein: ComponentPropertyCurves,
    /// Fat / lipid property curves.
    pub fat: ComponentPropertyCurves,
    /// Carbohydrate property curves.
    pub carbohydrate: ComponentPropertyCurves,
    /// Dietary fiber property curves.
    pub fiber: ComponentPropertyCurves,
    /// Ash / mineral property curves.
    pub ash: ComponentPropertyCurves,
    /// Pore-gas thermal conductivity [W/(m K)] (dry air default).
    pub gas_conductivity: QuadraticCelsiusPoly,
    /// Pore-gas specific gas constant [J/(kg K)] (dry air default `287.058`).
    pub gas_specific_constant_j_kg_k: f64,
    /// Reference pressure for pore-gas ideal density [Pa] (`101_325.0`).
    pub gas_reference_pressure_pa: f64,
    /// Pure water/ice latent heat of fusion at `0 C` [J/kg] (`333_600.0`).
    pub fusion_enthalpy_0c_j_kg: f64,
    /// Inclusive lower nominal calibration temperature [K] (`233.15 K = -40 C`).
    pub calibration_min_k: f64,
    /// Inclusive upper nominal calibration temperature [K] (`423.15 K = +150 C`).
    pub calibration_max_k: f64,
    /// Inclusive lower hard admitted temperature [K] (`223.15 K = -50 C`).
    pub admitted_min_k: f64,
    /// Inclusive upper hard admitted temperature [K] (`523.15 K = +250 C`).
    pub admitted_max_k: f64,
}

impl ChoiOkosTable {
    /// Canonical Choi–Okos (1986) / ASHRAE Handbook — Refrigeration Ch. 19 Table 1
    /// parameter set in coherent SI units, tagged `"empirical"`.
    #[must_use]
    pub fn ashrae_1986() -> Self {
        Self {
            provenance_tag: ASHRAE_CHOI_OKOS_SOURCE_TAG.to_string(),
            citation: ASHRAE_CHOI_OKOS_CITATION.to_string(),
            // Liquid water (0..150 C): 4.1762 - 9.0864e-5 t + 5.4731e-6 t^2 kJ/(kg K)
            water_cp_above_zero: QuadraticCelsiusPoly::new(4176.2, -9.0864e-2, 5.4731e-3),
            // Supercooled/unfrozen water (-40..0 C): 4.0817 - 5.3062e-3 t + 9.9516e-4 t^2 kJ/(kg K)
            water_cp_below_zero: QuadraticCelsiusPoly::new(4081.7, -5.3062, 9.9516e-1),
            // Water density: 9.9718e2 + 3.1439e-3 t - 3.7574e-3 t^2 kg/m^3
            water_density: QuadraticCelsiusPoly::new(997.18, 3.1439e-3, -3.7574e-3),
            // Water k: 5.7109e-1 + 1.7625e-3 t - 6.7036e-6 t^2 W/(m K)
            water_conductivity: QuadraticCelsiusPoly::new(5.7109e-1, 1.7625e-3, -6.7036e-6),
            // Ice cp: 2.0623 + 6.0769e-3 t kJ/(kg K)
            ice_cp: QuadraticCelsiusPoly::new(2062.3, 6.0769, 0.0),
            // Ice density: 9.1689e2 - 1.3071e-1 t kg/m^3
            ice_density: QuadraticCelsiusPoly::new(916.89, -1.3071e-1, 0.0),
            // Ice k: 2.2196 - 6.2489e-3 t + 1.0154e-4 t^2 W/(m K)
            ice_conductivity: QuadraticCelsiusPoly::new(2.2196, -6.2489e-3, 1.0154e-4),
            protein: ComponentPropertyCurves {
                // 2.0082 + 1.2089e-3 t - 1.3129e-6 t^2 kJ/(kg K)
                cp: QuadraticCelsiusPoly::new(2008.2, 1.2089, -1.3129e-3),
                // 1.3299e3 - 5.1840e-1 t kg/m^3
                density: QuadraticCelsiusPoly::new(1329.9, -5.1840e-1, 0.0),
                // 1.7881e-1 + 1.1958e-3 t - 2.7178e-6 t^2 W/(m K)
                conductivity: QuadraticCelsiusPoly::new(1.7881e-1, 1.1958e-3, -2.7178e-6),
            },
            fat: ComponentPropertyCurves {
                // 1.9842 + 1.4733e-3 t - 4.8008e-6 t^2 kJ/(kg K)
                cp: QuadraticCelsiusPoly::new(1984.2, 1.4733, -4.8008e-3),
                // 9.2559e2 - 4.1757e-1 t kg/m^3
                density: QuadraticCelsiusPoly::new(925.59, -4.1757e-1, 0.0),
                // ASHRAE erratum-corrected: 1.8071e-1 - 2.7604e-4 t - 1.7749e-7 t^2 W/(m K)
                conductivity: QuadraticCelsiusPoly::new(1.8071e-1, -2.7604e-4, -1.7749e-7),
            },
            carbohydrate: ComponentPropertyCurves {
                // 1.5488 + 1.9625e-3 t - 5.9399e-6 t^2 kJ/(kg K)
                cp: QuadraticCelsiusPoly::new(1548.8, 1.9625, -5.9399e-3),
                // 1.5991e3 - 3.1046e-1 t kg/m^3
                density: QuadraticCelsiusPoly::new(1599.1, -3.1046e-1, 0.0),
                // 2.0141e-1 + 1.3874e-3 t - 4.3312e-6 t^2 W/(m K)
                conductivity: QuadraticCelsiusPoly::new(2.0141e-1, 1.3874e-3, -4.3312e-6),
            },
            fiber: ComponentPropertyCurves {
                // 1.8459 + 1.8306e-3 t - 4.6509e-6 t^2 kJ/(kg K)
                cp: QuadraticCelsiusPoly::new(1845.9, 1.8306, -4.6509e-3),
                // 1.3115e3 - 3.6589e-1 t kg/m^3
                density: QuadraticCelsiusPoly::new(1311.5, -3.6589e-1, 0.0),
                // 1.8331e-1 + 1.2497e-3 t - 3.1683e-6 t^2 W/(m K)
                conductivity: QuadraticCelsiusPoly::new(1.8331e-1, 1.2497e-3, -3.1683e-6),
            },
            ash: ComponentPropertyCurves {
                // 1.0926 + 1.8896e-3 t - 3.6817e-6 t^2 kJ/(kg K)
                cp: QuadraticCelsiusPoly::new(1092.6, 1.8896, -3.6817e-3),
                // 2.4238e3 - 2.8063e-1 t kg/m^3
                density: QuadraticCelsiusPoly::new(2423.8, -2.8063e-1, 0.0),
                // 3.2962e-1 + 1.4011e-3 t - 2.9069e-6 t^2 W/(m K)
                conductivity: QuadraticCelsiusPoly::new(3.2962e-1, 1.4011e-3, -2.9069e-6),
            },
            // Dry air k(t_c) over -50..250 C (NIST reference match < 0.5%).
            gas_conductivity: QuadraticCelsiusPoly::new(2.414e-2, 7.57e-5, -2.6e-8),
            gas_specific_constant_j_kg_k: DRY_AIR_GAS_CONSTANT_J_KG_K,
            gas_reference_pressure_pa: STANDARD_ATMOSPHERE_PA,
            fusion_enthalpy_0c_j_kg: ICE_FUSION_ENTHALPY_0C_J_KG,
            calibration_min_k: MINIMUM_VALID_TEMPERATURE_K,
            calibration_max_k: MAXIMUM_VALID_TEMPERATURE_K,
            admitted_min_k: MINIMUM_ADMITTED_TEMPERATURE_K,
            admitted_max_k: MAXIMUM_ADMITTED_TEMPERATURE_K,
        }
    }

    /// Validate structural finiteness, positivity, and temperature-bound ordering.
    pub fn validate(&self) -> Result<(), ProximateErrorV1> {
        if self.provenance_tag.trim().is_empty() {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "table",
                property: "provenance_tag",
                reason: "provenance tag must not be empty",
            });
        }
        if self.citation.trim().is_empty() {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "table",
                property: "citation",
                reason: "bibliographic citation must not be empty",
            });
        }
        for (component, property, poly) in [
            ("water", "cp_above_zero", self.water_cp_above_zero),
            ("water", "cp_below_zero", self.water_cp_below_zero),
            ("water", "density", self.water_density),
            ("water", "conductivity", self.water_conductivity),
            ("ice", "cp", self.ice_cp),
            ("ice", "density", self.ice_density),
            ("ice", "conductivity", self.ice_conductivity),
            ("gas", "conductivity", self.gas_conductivity),
        ] {
            if !poly.is_finite() {
                return Err(ProximateErrorV1::InvalidCoefficientTable {
                    component,
                    property,
                    reason: "polynomial coefficients must be finite",
                });
            }
        }
        self.protein.validate("protein")?;
        self.fat.validate("fat")?;
        self.carbohydrate.validate("carbohydrate")?;
        self.fiber.validate("fiber")?;
        self.ash.validate("ash")?;

        if !self.gas_specific_constant_j_kg_k.is_finite()
            || self.gas_specific_constant_j_kg_k <= 0.0
        {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "gas",
                property: "specific_constant",
                reason: "gas specific constant must be finite and strictly positive",
            });
        }
        if !self.gas_reference_pressure_pa.is_finite() || self.gas_reference_pressure_pa <= 0.0 {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "gas",
                property: "reference_pressure",
                reason: "gas reference pressure must be finite and strictly positive",
            });
        }
        if !self.fusion_enthalpy_0c_j_kg.is_finite() || self.fusion_enthalpy_0c_j_kg <= 0.0 {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "ice",
                property: "fusion_enthalpy_0c",
                reason: "fusion enthalpy at 0 C must be finite and strictly positive",
            });
        }
        let bounds = [
            self.admitted_min_k,
            self.calibration_min_k,
            self.calibration_max_k,
            self.admitted_max_k,
        ];
        if !bounds.iter().all(|v| v.is_finite() && *v > 0.0)
            || self.admitted_min_k > self.calibration_min_k
            || self.calibration_min_k >= self.calibration_max_k
            || self.calibration_max_k > self.admitted_max_k
        {
            return Err(ProximateErrorV1::InvalidCoefficientTable {
                component: "table",
                property: "temperature_bounds",
                reason: "temperature bounds must satisfy 0 < admitted_min <= cal_min < cal_max <= admitted_max",
            });
        }
        Ok(())
    }

    /// Export all table parameters as canonical dimensioned `fs-matdb`
    /// [`LawParameter`] entries suitable for a `constitutive-v1` card with
    /// `law = "choi-okos-proximate"`, `law_version = 1`.
    #[must_use]
    pub fn to_card_parameters(&self) -> BTreeMap<String, LawParameter> {
        let mut params = BTreeMap::new();
        insert_poly(
            &mut params,
            "water_cp_pos",
            self.water_cp_above_zero,
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "water_cp_neg",
            self.water_cp_below_zero,
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "water_rho",
            self.water_density,
            RHO_C0_DIMS,
            RHO_C1_DIMS,
            RHO_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "water_k",
            self.water_conductivity,
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "ice_cp",
            self.ice_cp,
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "ice_rho",
            self.ice_density,
            RHO_C0_DIMS,
            RHO_C1_DIMS,
            RHO_C2_DIMS,
        );
        insert_poly(
            &mut params,
            "ice_k",
            self.ice_conductivity,
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
        );
        for (prefix, curves) in [
            ("protein", self.protein),
            ("fat", self.fat),
            ("carbohydrate", self.carbohydrate),
            ("fiber", self.fiber),
            ("ash", self.ash),
        ] {
            insert_poly(
                &mut params,
                &format!("{prefix}_cp"),
                curves.cp,
                CP_C0_DIMS,
                CP_C1_DIMS,
                CP_C2_DIMS,
            );
            insert_poly(
                &mut params,
                &format!("{prefix}_rho"),
                curves.density,
                RHO_C0_DIMS,
                RHO_C1_DIMS,
                RHO_C2_DIMS,
            );
            insert_poly(
                &mut params,
                &format!("{prefix}_k"),
                curves.conductivity,
                K_C0_DIMS,
                K_C1_DIMS,
                K_C2_DIMS,
            );
        }
        insert_poly(
            &mut params,
            "gas_k",
            self.gas_conductivity,
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
        );
        params.insert(
            "gas_specific_constant".to_string(),
            LawParameter {
                value: self.gas_specific_constant_j_kg_k,
                dims: CP_C0_DIMS,
            },
        );
        params.insert(
            "gas_reference_pressure".to_string(),
            LawParameter {
                value: self.gas_reference_pressure_pa,
                dims: Pressure::DIMS,
            },
        );
        params.insert(
            "fusion_enthalpy_0c".to_string(),
            LawParameter {
                value: self.fusion_enthalpy_0c_j_kg,
                dims: SPECIFIC_ENERGY_DIMS,
            },
        );
        params
    }

    /// Load and strictly validate a [`ChoiOkosTable`] from an immutable
    /// `fs-matdb` [`ConstitutiveModelCard`].
    ///
    /// Requires `law == "choi-okos-proximate"`, `law_version == 1`, stateless
    /// policy (`state_schema_version == 0`, `ZeroInternalState`), exact
    /// dimensioned parameter block with zero foreign parameters, and a valid
    /// positive finite `T` validity axis.
    pub fn from_card(card: &ConstitutiveModelCard) -> Result<Self, ProximateErrorV1> {
        let (cal_min_k, cal_max_k, admitted_min_k, admitted_max_k) = validate_card_header(card)?;

        let mut consumed = 0usize;
        let water_cp_above_zero = read_poly(
            card,
            "water_cp_pos",
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
            &mut consumed,
        )?;
        let water_cp_below_zero = read_poly(
            card,
            "water_cp_neg",
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
            &mut consumed,
        )?;
        let water_density = read_poly(
            card,
            "water_rho",
            RHO_C0_DIMS,
            RHO_C1_DIMS,
            RHO_C2_DIMS,
            &mut consumed,
        )?;
        let water_conductivity = read_poly(
            card,
            "water_k",
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
            &mut consumed,
        )?;
        let ice_cp = read_poly(
            card,
            "ice_cp",
            CP_C0_DIMS,
            CP_C1_DIMS,
            CP_C2_DIMS,
            &mut consumed,
        )?;
        let ice_density = read_poly(
            card,
            "ice_rho",
            RHO_C0_DIMS,
            RHO_C1_DIMS,
            RHO_C2_DIMS,
            &mut consumed,
        )?;
        let ice_conductivity = read_poly(
            card,
            "ice_k",
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
            &mut consumed,
        )?;
        let protein = read_component_curves(card, "protein", &mut consumed)?;
        let fat = read_component_curves(card, "fat", &mut consumed)?;
        let carbohydrate = read_component_curves(card, "carbohydrate", &mut consumed)?;
        let fiber = read_component_curves(card, "fiber", &mut consumed)?;
        let ash = read_component_curves(card, "ash", &mut consumed)?;
        let gas_conductivity = read_poly(
            card,
            "gas_k",
            K_C0_DIMS,
            K_C1_DIMS,
            K_C2_DIMS,
            &mut consumed,
        )?;
        let gas_specific_constant_j_kg_k =
            read_scalar(card, "gas_specific_constant", CP_C0_DIMS, &mut consumed)?;
        let gas_reference_pressure_pa = read_scalar(
            card,
            "gas_reference_pressure",
            Pressure::DIMS,
            &mut consumed,
        )?;
        let fusion_enthalpy_0c_j_kg = read_scalar(
            card,
            "fusion_enthalpy_0c",
            SPECIFIC_ENERGY_DIMS,
            &mut consumed,
        )?;

        if card.parameters.len() != consumed {
            return Err(ProximateErrorV1::MaterialCard {
                reason: format!(
                    "card contains unexpected foreign parameters (expected {consumed}, found {})",
                    card.parameters.len()
                ),
            });
        }

        let table = Self {
            provenance_tag: ASHRAE_CHOI_OKOS_SOURCE_TAG.to_string(),
            citation: card.provenance.source.clone(),
            water_cp_above_zero,
            water_cp_below_zero,
            water_density,
            water_conductivity,
            ice_cp,
            ice_density,
            ice_conductivity,
            protein,
            fat,
            carbohydrate,
            fiber,
            ash,
            gas_conductivity,
            gas_specific_constant_j_kg_k,
            gas_reference_pressure_pa,
            fusion_enthalpy_0c_j_kg,
            calibration_min_k: cal_min_k,
            calibration_max_k: cal_max_k,
            admitted_min_k,
            admitted_max_k,
        };
        table.validate()?;
        Ok(table)
    }

    /// Pure liquid water specific heat `c_{p,w}(t_c)` [J/(kg K)] (piecewise at `0 C`).
    #[must_use]
    pub fn water_cp(&self, t_c: f64) -> f64 {
        if t_c >= 0.0 {
            self.water_cp_above_zero.eval(t_c)
        } else {
            self.water_cp_below_zero.eval(t_c)
        }
    }

    /// Pure liquid water specific enthalpy `h_w(t_c)` [J/kg] relative to `0 C` liquid water.
    #[must_use]
    pub fn water_enthalpy_0c_ref(&self, t_c: f64) -> f64 {
        if t_c >= 0.0 {
            self.water_cp_above_zero.integrate_from_zero(t_c)
        } else {
            self.water_cp_below_zero.integrate_from_zero(t_c)
        }
    }

    /// Pure ice specific enthalpy `h_i(t_c)` [J/kg] relative to `0 C` liquid water
    /// (`h_i(0 C) = -L_{f,0}`), guaranteeing automatic Kirchhoff consistency
    /// `dL_f/dT = c_{p,w}(T) - c_{p,ice}(T)`.
    #[must_use]
    pub fn ice_enthalpy_0c_ref(&self, t_c: f64) -> f64 {
        -self.fusion_enthalpy_0c_j_kg + self.ice_cp.integrate_from_zero(t_c)
    }

    /// Kirchhoff-consistent latent heat of fusion `L_f(t_c) = h_w(t_c) - h_i(t_c)` [J/kg].
    #[must_use]
    pub fn latent_heat_fusion(&self, t_c: f64) -> f64 {
        self.fusion_enthalpy_0c_j_kg
            + (self.water_enthalpy_0c_ref(t_c) - self.ice_cp.integrate_from_zero(t_c))
    }

    /// Constant shift [J/kg] mapping the table's `0 C` liquid-water enthalpy datum
    /// to the shared `janaf-298.15K-formation` datum (`h_f0(298.15 K) - h_w(25 C)`).
    #[must_use]
    pub fn water_janaf_datum_offset_j_kg(&self) -> f64 {
        let t_298_c = JANAF_REFERENCE_TEMPERATURE_K - ICE_POINT_K;
        H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG - self.water_enthalpy_0c_ref(t_298_c)
    }
}

fn validate_card_header(
    card: &ConstitutiveModelCard,
) -> Result<(f64, f64, f64, f64), ProximateErrorV1> {
    card.validate()
        .map_err(|err| ProximateErrorV1::MaterialCard {
            reason: format!("card failed fs-matdb validation: {err}"),
        })?;
    if card.law.0 != CHOI_OKOS_LAW_ID_V1 {
        return Err(ProximateErrorV1::MaterialCard {
            reason: format!(
                "expected law id {CHOI_OKOS_LAW_ID_V1:?}, found {:?}",
                card.law.0
            ),
        });
    }
    if card.law_version != CHOI_OKOS_LAW_VERSION_V1 {
        return Err(ProximateErrorV1::MaterialCard {
            reason: format!(
                "expected law version {CHOI_OKOS_LAW_VERSION_V1}, found {}",
                card.law_version
            ),
        });
    }
    if card.state_schema_version != 0 || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(ProximateErrorV1::MaterialCard {
            reason: "Choi-Okos card must be stateless (state_schema_version=0, ZeroInternalState)"
                .to_string(),
        });
    }
    if let Some(axis) = card
        .validity
        .bounds()
        .keys()
        .find(|axis| axis.as_str() != "T" && axis.as_str() != "T_admitted")
    {
        return Err(ProximateErrorV1::MaterialCard {
            reason: format!("unexpected validity axis {axis:?} on Choi-Okos card"),
        });
    }
    let (cal_min_k, cal_max_k) =
        card.validity
            .bound("T")
            .ok_or_else(|| ProximateErrorV1::MaterialCard {
                reason: "card validity is missing the required \"T\" temperature axis".to_string(),
            })?;
    if !cal_min_k.is_finite()
        || !cal_max_k.is_finite()
        || cal_min_k <= 0.0
        || cal_min_k >= cal_max_k
    {
        return Err(ProximateErrorV1::MaterialCard {
            reason: "card \"T\" validity bounds must be positive, finite, and strictly increasing"
                .to_string(),
        });
    }
    let (admitted_min_k, admitted_max_k) = card.validity.bound("T_admitted").unwrap_or((
        cal_min_k.min(MINIMUM_ADMITTED_TEMPERATURE_K),
        cal_max_k.max(MAXIMUM_ADMITTED_TEMPERATURE_K),
    ));
    if !admitted_min_k.is_finite()
        || !admitted_max_k.is_finite()
        || admitted_min_k <= 0.0
        || admitted_min_k > cal_min_k
        || admitted_max_k < cal_max_k
    {
        return Err(ProximateErrorV1::MaterialCard {
            reason: "card \"T_admitted\" bounds must be positive, finite, and enclose \"T\""
                .to_string(),
        });
    }
    Ok((cal_min_k, cal_max_k, admitted_min_k, admitted_max_k))
}

fn insert_poly(
    params: &mut BTreeMap<String, LawParameter>,
    prefix: &str,
    poly: QuadraticCelsiusPoly,
    d0: Dims,
    d1: Dims,
    d2: Dims,
) {
    params.insert(
        format!("{prefix}_c0"),
        LawParameter {
            value: poly.c0,
            dims: d0,
        },
    );
    params.insert(
        format!("{prefix}_c1"),
        LawParameter {
            value: poly.c1,
            dims: d1,
        },
    );
    params.insert(
        format!("{prefix}_c2"),
        LawParameter {
            value: poly.c2,
            dims: d2,
        },
    );
}

fn read_scalar(
    card: &ConstitutiveModelCard,
    name: &str,
    expected_dims: Dims,
    consumed: &mut usize,
) -> Result<f64, ProximateErrorV1> {
    let param = card
        .parameters
        .get(name)
        .ok_or_else(|| ProximateErrorV1::MaterialCard {
            reason: format!("missing required parameter {name:?}"),
        })?;
    if param.dims != expected_dims {
        return Err(ProximateErrorV1::MaterialCard {
            reason: format!(
                "parameter {name:?} has dims {:?}, expected {expected_dims:?}",
                param.dims
            ),
        });
    }
    *consumed += 1;
    Ok(param.value)
}

fn read_poly(
    card: &ConstitutiveModelCard,
    prefix: &str,
    d0: Dims,
    d1: Dims,
    d2: Dims,
    consumed: &mut usize,
) -> Result<QuadraticCelsiusPoly, ProximateErrorV1> {
    let c0 = read_scalar(card, &format!("{prefix}_c0"), d0, consumed)?;
    let c1 = read_scalar(card, &format!("{prefix}_c1"), d1, consumed)?;
    let c2 = read_scalar(card, &format!("{prefix}_c2"), d2, consumed)?;
    Ok(QuadraticCelsiusPoly::new(c0, c1, c2))
}

fn read_component_curves(
    card: &ConstitutiveModelCard,
    prefix: &str,
    consumed: &mut usize,
) -> Result<ComponentPropertyCurves, ProximateErrorV1> {
    let cp = read_poly(
        card,
        &format!("{prefix}_cp"),
        CP_C0_DIMS,
        CP_C1_DIMS,
        CP_C2_DIMS,
        consumed,
    )?;
    let density = read_poly(
        card,
        &format!("{prefix}_rho"),
        RHO_C0_DIMS,
        RHO_C1_DIMS,
        RHO_C2_DIMS,
        consumed,
    )?;
    let conductivity = read_poly(
        card,
        &format!("{prefix}_k"),
        K_C0_DIMS,
        K_C1_DIMS,
        K_C2_DIMS,
        consumed,
    )?;
    Ok(ComponentPropertyCurves {
        cp,
        density,
        conductivity,
    })
}

/// Structural thermal-conductivity mixing law for multiphase proximate foods.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConductivityMixingModel {
    /// Parallel (Voigt) upper Wiener bound: `k_par = sum(eps_j * k_j)`.
    Parallel,
    /// Series (Reuss) lower Wiener bound: `k_ser = (sum(eps_j / k_j))^-1`.
    Series,
    /// Krischer harmonic blend between parallel and series:
    /// `k = ((1 - f) / k_par + f / k_ser)^-1`, with `anisotropy_factor` (`f` in `[0, 1]`)
    /// equal to `0.0` for parallel-to-fiber heat flow and `1.0` for perpendicular-to-fiber flow.
    Krischer {
        /// Series distribution factor in `[0.0, 1.0]`.
        anisotropy_factor: f64,
    },
    /// Volume-weighted geometric mean: `k_geom = exp(sum(eps_j * ln(k_j)))`.
    Geometric,
    /// Kopelman (1966) isotropic dispersed-phase structural model.
    KopelmanIsotropic,
    /// Multi-component Maxwell–Eucken dispersion model (two-stage for porous foods
    /// with continuous condensed food matrix and dispersed gas pores).
    MaxwellEucken,
}

/// Evaluation configuration controlling mushy-zone regularization and effective
/// thermal-conductivity mixing selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProximateEvalConfig {
    /// `C^1` mushy-zone transition width `delta_T` [K] across `[T_f - delta_T, T_f]`.
    /// Default: `0.05 K` (Round-2 Physics #7).
    pub mushy_zone_width_k: f64,
    /// Selected effective thermal-conductivity mixing model.
    pub conductivity_model: ConductivityMixingModel,
}

impl Default for ProximateEvalConfig {
    fn default() -> Self {
        Self {
            mushy_zone_width_k: DEFAULT_MUSHY_ZONE_WIDTH_K,
            conductivity_model: ConductivityMixingModel::MaxwellEucken,
        }
    }
}

impl ProximateEvalConfig {
    /// Validate configuration parameters.
    pub fn validate(&self) -> Result<(), ProximateErrorV1> {
        if !self.mushy_zone_width_k.is_finite() || self.mushy_zone_width_k <= 0.0 {
            return Err(ProximateErrorV1::InvalidMushyZoneWidth {
                delta_t_k_bits: self.mushy_zone_width_k.to_bits(),
            });
        }
        if let ConductivityMixingModel::Krischer { anisotropy_factor } = self.conductivity_model
            && (!anisotropy_factor.is_finite() || !(0.0..=1.0).contains(&anisotropy_factor))
        {
            return Err(ProximateErrorV1::InvalidConductivityParameter {
                parameter: "anisotropy_factor",
                bits: anisotropy_factor.to_bits(),
                reason: "Krischer anisotropy factor must be finite and in [0.0, 1.0]",
            });
        }
        Ok(())
    }
}

/// Complete immutable Choi–Okos proximate evaluation result at one temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProximateEvalResult {
    temperature_k: f64,
    validity: ValidityStatus,
    ice_mass_fraction: f64,
    unfrozen_water_mass_fraction: f64,
    d_ice_fraction_d_t: f64,
    density_condensed_kg_m3: f64,
    density_bulk_kg_m3: f64,
    cp_sensible_j_kg_k: f64,
    cp_apparent_j_kg_k: f64,
    latent_heat_fusion_j_kg: f64,
    enthalpy_0c_ref_j_kg: f64,
    enthalpy_janaf_j_kg: f64,
    k_parallel_w_m_k: f64,
    k_series_w_m_k: f64,
    k_maxwell_eucken_w_m_k: f64,
    k_kopelman_w_m_k: f64,
    k_effective_w_m_k: f64,
    thermal_diffusivity_sensible_m2_s: f64,
    thermal_diffusivity_apparent_m2_s: f64,
}

impl ProximateEvalResult {
    /// Evaluated temperature [K].
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        Temperature::new(self.temperature_k)
    }

    /// Physical calibration / desiccation validity classification at this temperature.
    #[must_use]
    pub const fn validity(&self) -> ValidityStatus {
        self.validity
    }

    /// Dynamic ice mass fraction `X_i(T)` on the condensed food mass basis.
    #[must_use]
    pub const fn ice_mass_fraction(&self) -> f64 {
        self.ice_mass_fraction
    }

    /// Unfrozen liquid water mass fraction `X_{uw}(T) = X_w - X_i(T)`.
    #[must_use]
    pub const fn unfrozen_water_mass_fraction(&self) -> f64 {
        self.unfrozen_water_mass_fraction
    }

    /// Analytical temperature derivative `dX_i/dT` [1/K] (`<= 0.0`).
    #[must_use]
    pub const fn d_ice_fraction_d_t(&self) -> f64 {
        self.d_ice_fraction_d_t
    }

    /// Effective bulk density [kg/m^3] (equals condensed density when `porosity_gas == 0`).
    #[must_use]
    pub const fn density(&self) -> Density {
        Density::new(self.density_bulk_kg_m3)
    }

    /// Pore-free condensed matrix density [kg/m^3].
    #[must_use]
    pub const fn density_condensed(&self) -> Density {
        Density::new(self.density_condensed_kg_m3)
    }

    /// Frozen-phase sensible specific heat `sum_alpha X_alpha(T) c_{p,alpha}(T)` [J/(kg K)].
    #[must_use]
    pub const fn cp_sensible(&self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cp_sensible_j_kg_k)
    }

    /// Exact analytical apparent specific heat `c_{p,app}(T) = dh/dT` [J/(kg K)],
    /// incorporating Kirchhoff-consistent latent heat release/absorption across freezing.
    #[must_use]
    pub const fn cp_apparent(&self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cp_apparent_j_kg_k)
    }

    /// Kirchhoff-consistent latent heat of fusion `L_f(T) = h_w(T) - h_i(T)` [J/kg].
    #[must_use]
    pub const fn latent_heat_fusion(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.latent_heat_fusion_j_kg)
    }

    /// Total specific enthalpy `sum_alpha X_alpha(T) h_alpha(T)` [J/kg] anchored at
    /// `0 C` liquid water and `0 C` dry solids (`h_ice(0 C) = -L_{f,0}`).
    #[must_use]
    pub const fn specific_enthalpy_0c_datum(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_0c_ref_j_kg)
    }

    /// Total specific enthalpy [J/kg] shifted onto the shared `janaf-298.15K-formation`
    /// datum for water mass fraction `X_w`, matching NASA-9 / IF97 kitchen streams.
    #[must_use]
    pub const fn specific_enthalpy_janaf(&self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_janaf_j_kg)
    }

    /// Parallel (Voigt) upper Wiener thermal-conductivity bound [W/(m K)].
    #[must_use]
    pub const fn k_parallel(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.k_parallel_w_m_k)
    }

    /// Series (Reuss) lower Wiener thermal-conductivity bound [W/(m K)].
    #[must_use]
    pub const fn k_series(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.k_series_w_m_k)
    }

    /// Maxwell–Eucken effective thermal conductivity [W/(m K)].
    #[must_use]
    pub const fn k_maxwell_eucken(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.k_maxwell_eucken_w_m_k)
    }

    /// Kopelman (1966) isotropic effective thermal conductivity [W/(m K)].
    #[must_use]
    pub const fn k_kopelman(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.k_kopelman_w_m_k)
    }

    /// Selected effective thermal conductivity per [`ProximateEvalConfig::conductivity_model`] [W/(m K)].
    #[must_use]
    pub const fn k_effective(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.k_effective_w_m_k)
    }

    /// Sensible thermal diffusivity `k_effective / (rho_bulk * cp_sensible)` [m^2/s].
    #[must_use]
    pub const fn thermal_diffusivity_sensible_m2_s(&self) -> f64 {
        self.thermal_diffusivity_sensible_m2_s
    }

    /// Apparent thermal diffusivity `k_effective / (rho_bulk * cp_apparent)` [m^2/s].
    #[must_use]
    pub const fn thermal_diffusivity_apparent_m2_s(&self) -> f64 {
        self.thermal_diffusivity_apparent_m2_s
    }
}

/// Validated, allocation-free Choi–Okos proximate thermophysical evaluator.
#[derive(Debug, Clone, PartialEq)]
pub struct ChoiOkosEvaluator {
    table: ChoiOkosTable,
    config: ProximateEvalConfig,
}

impl ChoiOkosEvaluator {
    /// Construct an evaluator using the canonical ASHRAE Ch. 19 / Choi–Okos (1986)
    /// empirical table and default `0.05 K` mushy-zone regularization + Maxwell–Eucken mixing.
    #[must_use]
    pub fn ashrae_default() -> Self {
        Self {
            table: ChoiOkosTable::ashrae_1986(),
            config: ProximateEvalConfig::default(),
        }
    }

    /// Construct an evaluator from an explicit coefficient table and evaluation config.
    pub fn new(
        table: ChoiOkosTable,
        config: ProximateEvalConfig,
    ) -> Result<Self, ProximateErrorV1> {
        table.validate()?;
        config.validate()?;
        Ok(Self { table, config })
    }

    /// Return a new evaluator with an updated validated [`ProximateEvalConfig`].
    pub fn with_config(mut self, config: ProximateEvalConfig) -> Result<Self, ProximateErrorV1> {
        config.validate()?;
        self.config = config;
        Ok(self)
    }

    /// Retained Choi–Okos coefficient table.
    #[must_use]
    pub const fn table(&self) -> &ChoiOkosTable {
        &self.table
    }

    /// Retained evaluation configuration.
    #[must_use]
    pub const fn config(&self) -> ProximateEvalConfig {
        self.config
    }

    /// Evaluate the piecewise `C^1`-regularized Chen–Heldman ice mass fraction
    /// `X_i(T)` and its exact analytical temperature derivative `dX_i/dT` [1/K].
    ///
    /// Guarantees:
    /// - `X_i(T) == 0.0` and `dX_i/dT == 0.0` strictly for `T >= T_f`.
    /// - Exact Chen–Heldman `(X_w - x_b) * (1 - (T_f - 273.15) / (T - 273.15))`
    ///   for `T <= T_f - delta_T`.
    /// - Monotone, overshoot-free `C^1` cubic Hermite mushy zone over
    ///   `[T_f - delta_T, T_f]`, preventing `0 / 0 = NaN` or infinite slope for
    ///   pure water (`T_f = 273.15 K`) or small freezing depressions
    ///   (`273.15 - T_f < delta_T`).
    pub fn ice_fraction_and_derivative(
        &self,
        composition: &ProximateComposition,
        temperature: Temperature,
    ) -> Result<(f64, f64), ProximateErrorV1> {
        let t_k = temperature.value();
        self.admit_temperature(t_k)?;
        Ok(regularized_chen_heldman_ice_fraction(
            composition.freezable_water(),
            composition.t_f_k,
            t_k,
            self.config.mushy_zone_width_k,
        ))
    }

    /// Evaluate all proximate thermophysical properties at `temperature` without heap allocation.
    pub fn evaluate(
        &self,
        composition: &ProximateComposition,
        temperature: Temperature,
    ) -> Result<ProximateEvalResult, ProximateErrorV1> {
        let t_k = temperature.value();
        let validity = self.admit_temperature(t_k)?;
        let t_c = t_k - ICE_POINT_K;

        let (x_ice, dx_ice_dt) = regularized_chen_heldman_ice_fraction(
            composition.freezable_water(),
            composition.t_f_k,
            t_k,
            self.config.mushy_zone_width_k,
        );
        let x_uw = (composition.water - x_ice).max(0.0);
        let mass_fractions = [
            x_uw,
            x_ice,
            composition.protein,
            composition.fat,
            composition.carbohydrate,
            composition.fiber,
            composition.ash,
        ];

        let (cps, rhos, ks, enthalpies, latent_fusion) =
            eval_pure_component_arrays(&self.table, t_c, &mass_fractions, dx_ice_dt)?;

        let (density_condensed, eps_cond) = condensed_density_and_volumes(&mass_fractions, &rhos)?;
        let (density_bulk, k_gas) = bulk_density_and_gas_k(
            &self.table,
            density_condensed,
            composition.porosity_gas,
            t_k,
            t_c,
        )?;

        // Direct phase-enthalpy summation (automatic Kirchhoff consistency):
        // h(T) = sum_alpha X_alpha(T) * h_alpha(T)
        let mut cp_sensible = 0.0;
        let mut enthalpy_0c_ref = 0.0;
        for idx in 0..7 {
            if mass_fractions[idx] > 0.0 {
                cp_sensible += mass_fractions[idx] * cps[idx];
                enthalpy_0c_ref += mass_fractions[idx] * enthalpies[idx];
            }
        }
        // Analytical derivative of sum X_alpha(T) h_alpha(T):
        // dh/dT = cp_sensible - dX_i/dT * (h_w - h_i) = cp_sensible - dX_i/dT * L_f(T).
        let cp_apparent = cp_sensible - dx_ice_dt * latent_fusion;
        let enthalpy_janaf =
            enthalpy_0c_ref + composition.water * self.table.water_janaf_datum_offset_j_kg();

        let (k_parallel, k_series, k_me, k_kop, k_geom) =
            compute_conductivities(&eps_cond, &ks, composition.porosity_gas, k_gas);
        let k_effective = select_effective_conductivity(
            self.config.conductivity_model,
            k_parallel,
            k_series,
            k_me,
            k_kop,
            k_geom,
        );
        let alpha_sens = k_effective / (density_bulk * cp_sensible);
        let alpha_app = k_effective / (density_bulk * cp_apparent);

        let result = ProximateEvalResult {
            temperature_k: t_k,
            validity,
            ice_mass_fraction: x_ice,
            unfrozen_water_mass_fraction: x_uw,
            d_ice_fraction_d_t: dx_ice_dt,
            density_condensed_kg_m3: density_condensed,
            density_bulk_kg_m3: density_bulk,
            cp_sensible_j_kg_k: cp_sensible,
            cp_apparent_j_kg_k: cp_apparent,
            latent_heat_fusion_j_kg: latent_fusion,
            enthalpy_0c_ref_j_kg: enthalpy_0c_ref,
            enthalpy_janaf_j_kg: enthalpy_janaf,
            k_parallel_w_m_k: k_parallel,
            k_series_w_m_k: k_series,
            k_maxwell_eucken_w_m_k: k_me,
            k_kopelman_w_m_k: k_kop,
            k_effective_w_m_k: k_effective,
            thermal_diffusivity_sensible_m2_s: alpha_sens,
            thermal_diffusivity_apparent_m2_s: alpha_app,
        };
        validate_eval_result(&result)?;
        Ok(result)
    }

    /// Evaluate specific enthalpy relative to ASHRAE Table 7's `-40 C` (`233.15 K`)
    /// reference datum: `h(T) - h(233.15 K)` [J/kg].
    pub fn enthalpy_ashrae_minus40c(
        &self,
        composition: &ProximateComposition,
        temperature: Temperature,
    ) -> Result<MassSpecificEnergyQuantityV1, ProximateErrorV1> {
        let state = self.evaluate(composition, temperature)?;
        let ref_state = self.evaluate(
            composition,
            Temperature::new(ASHRAE_ENTHALPY_REFERENCE_TEMPERATURE_K),
        )?;
        let dh = state.enthalpy_0c_ref_j_kg - ref_state.enthalpy_0c_ref_j_kg;
        if !dh.is_finite() {
            return Err(ProximateErrorV1::NonFiniteEvaluation {
                property: "enthalpy_ashrae_minus40c",
                bits: dh.to_bits(),
            });
        }
        Ok(MassSpecificEnergyQuantityV1::new(dh))
    }

    fn admit_temperature(&self, t_k: f64) -> Result<ValidityStatus, ProximateErrorV1> {
        if !t_k.is_finite() || t_k <= 0.0 {
            return Err(ProximateErrorV1::InvalidTemperature {
                temperature_k_bits: t_k.to_bits(),
            });
        }
        let status = classify_temperature_with_bounds(
            t_k,
            self.table.admitted_min_k,
            self.table.calibration_min_k,
            self.table.calibration_max_k,
            self.table.admitted_max_k,
        );
        if status == ValidityStatus::RefusedOutOfRange {
            return Err(ProximateErrorV1::TemperatureOutOfRange {
                temperature_k: t_k,
                minimum_k: self.table.admitted_min_k,
                maximum_k: self.table.admitted_max_k,
            });
        }
        Ok(status)
    }
}

/// Piecewise `C^1`-regularized Chen–Heldman ice mass fraction `X_i(T)` and
/// analytical derivative `dX_i/dT` [1/K].
///
/// Regions:
/// 1. `T >= T_f`: `X_i = 0`, `dX_i/dT = 0` (strictly unfrozen).
/// 2. `T <= T_f - delta_T`: exact Celsius-offset Chen–Heldman
///    `f = (T_f - T) / (273.15 - T) = 1 - (T_f - 273.15) / (T - 273.15)`,
///    with `f = 1, df/dT = 0` when `T_f == 273.15 K` (pure water).
/// 3. `T_f - delta_T < T < T_f`: monotone cubic Hermite transition matching
///    both value and derivative at `T_f - delta_T` and `(0, 0)` at `T_f`.
///    With `r = (273.15 - T_f) / ((273.15 - T_f) + delta_T) in [0, 1)`,
///    `f_m = delta_T / ((273.15 - T_f) + delta_T)`, and
///    `s = (T_f - T) / delta_T in (0, 1)`:
///    `f(s) = f_m * s^2 * ((3 - r) - (2 - r) * s)` and
///    `df/dT = -(f_m / delta_T) * s * (2 * (3 - r) - 3 * (2 - r) * s) <= 0`.
fn regularized_chen_heldman_ice_fraction(
    freezable_water: f64,
    t_f_k: f64,
    t_k: f64,
    delta_t_k: f64,
) -> (f64, f64) {
    if freezable_water <= 0.0 || t_k >= t_f_k {
        return (0.0, 0.0);
    }
    let delta_below_tf = t_f_k - t_k;
    let abs_tf = ICE_POINT_K - t_f_k;
    if delta_below_tf < delta_t_k {
        let denom = abs_tf + delta_t_k;
        let f_m = delta_t_k / denom;
        let r = abs_tf / denom;
        let s = delta_below_tf / delta_t_k;
        let a = 3.0 - r;
        let b = 2.0 - r;
        let f = (f_m * (s * s) * (a - b * s)).clamp(0.0, 1.0);
        let slope_inner = (2.0 * a - 3.0 * b * s).max(0.0);
        let df_dt = -(f_m / delta_t_k) * s * slope_inner;
        (freezable_water * f, freezable_water * df_dt)
    } else if abs_tf == 0.0 {
        (freezable_water, 0.0)
    } else {
        // Cancellation-free exact Celsius-offset Chen-Heldman:
        // 1 - t_f / t = (t_f - t) / (273.15 - t) = delta_below_tf / abs_t.
        let abs_t = ICE_POINT_K - t_k;
        let f = (delta_below_tf / abs_t).clamp(0.0, 1.0);
        let df_dt = -abs_tf / (abs_t * abs_t);
        (freezable_water * f, freezable_water * df_dt)
    }
}

type PureComponentEvaluation = ([f64; 7], [f64; 7], [f64; 7], [f64; 7], f64);

fn eval_pure_component_arrays(
    table: &ChoiOkosTable,
    t_c: f64,
    mass_fractions: &[f64; 7],
    dx_ice_dt: f64,
) -> Result<PureComponentEvaluation, ProximateErrorV1> {
    const NAMES: [&str; 7] = [
        "unfrozen_water",
        "ice",
        "protein",
        "fat",
        "carbohydrate",
        "fiber",
        "ash",
    ];
    let cps = [
        table.water_cp(t_c),
        table.ice_cp.eval(t_c),
        table.protein.cp.eval(t_c),
        table.fat.cp.eval(t_c),
        table.carbohydrate.cp.eval(t_c),
        table.fiber.cp.eval(t_c),
        table.ash.cp.eval(t_c),
    ];
    let rhos = [
        table.water_density.eval(t_c),
        table.ice_density.eval(t_c),
        table.protein.density.eval(t_c),
        table.fat.density.eval(t_c),
        table.carbohydrate.density.eval(t_c),
        table.fiber.density.eval(t_c),
        table.ash.density.eval(t_c),
    ];
    let ks = [
        table.water_conductivity.eval(t_c),
        table.ice_conductivity.eval(t_c),
        table.protein.conductivity.eval(t_c),
        table.fat.conductivity.eval(t_c),
        table.carbohydrate.conductivity.eval(t_c),
        table.fiber.conductivity.eval(t_c),
        table.ash.conductivity.eval(t_c),
    ];
    for idx in 0..7 {
        if mass_fractions[idx] > 0.0 {
            ensure_positive_component(NAMES[idx], "cp", cps[idx])?;
            ensure_positive_component(NAMES[idx], "density", rhos[idx])?;
            ensure_positive_component(NAMES[idx], "conductivity", ks[idx])?;
        }
    }
    let latent_fusion = table.latent_heat_fusion(t_c);
    if dx_ice_dt != 0.0 {
        ensure_positive_component("unfrozen_water", "cp", cps[0])?;
        ensure_positive_component("ice", "cp", cps[1])?;
        ensure_positive_component("ice", "latent_heat_fusion", latent_fusion)?;
    } else if mass_fractions[1] > 0.0 {
        ensure_positive_component("ice", "latent_heat_fusion", latent_fusion)?;
    }
    let enthalpies = [
        table.water_enthalpy_0c_ref(t_c),
        table.ice_enthalpy_0c_ref(t_c),
        table.protein.cp.integrate_from_zero(t_c),
        table.fat.cp.integrate_from_zero(t_c),
        table.carbohydrate.cp.integrate_from_zero(t_c),
        table.fiber.cp.integrate_from_zero(t_c),
        table.ash.cp.integrate_from_zero(t_c),
    ];
    Ok((cps, rhos, ks, enthalpies, latent_fusion))
}

fn condensed_density_and_volumes(
    mass_fractions: &[f64; 7],
    rhos: &[f64; 7],
) -> Result<(f64, [f64; 7]), ProximateErrorV1> {
    let mut specific_volume_cond = 0.0;
    for idx in 0..7 {
        if mass_fractions[idx] > 0.0 {
            specific_volume_cond += mass_fractions[idx] / rhos[idx];
        }
    }
    if !specific_volume_cond.is_finite() || specific_volume_cond <= 0.0 {
        return Err(ProximateErrorV1::NonFiniteEvaluation {
            property: "density_condensed",
            bits: specific_volume_cond.to_bits(),
        });
    }
    let density_condensed = 1.0 / specific_volume_cond;
    let mut eps_cond = [0.0; 7];
    for idx in 0..7 {
        if mass_fractions[idx] > 0.0 {
            eps_cond[idx] = (mass_fractions[idx] / rhos[idx]) * density_condensed;
        }
    }
    Ok((density_condensed, eps_cond))
}

fn bulk_density_and_gas_k(
    table: &ChoiOkosTable,
    density_condensed: f64,
    phi_g: f64,
    t_k: f64,
    t_c: f64,
) -> Result<(f64, f64), ProximateErrorV1> {
    if phi_g > 0.0 {
        let rho_gas = table.gas_reference_pressure_pa / (table.gas_specific_constant_j_kg_k * t_k);
        let k_g = table.gas_conductivity.eval(t_c);
        ensure_positive_component("gas", "density", rho_gas)?;
        ensure_positive_component("gas", "conductivity", k_g)?;
        Ok(((1.0 - phi_g) * density_condensed + phi_g * rho_gas, k_g))
    } else {
        Ok((density_condensed, 0.0))
    }
}

fn select_effective_conductivity(
    model: ConductivityMixingModel,
    k_parallel: f64,
    k_series: f64,
    k_me: f64,
    k_kop: f64,
    k_geom: f64,
) -> f64 {
    match model {
        ConductivityMixingModel::Parallel => k_parallel,
        ConductivityMixingModel::Series => k_series,
        ConductivityMixingModel::Krischer { anisotropy_factor } => {
            1.0 / ((1.0 - anisotropy_factor) / k_parallel + anisotropy_factor / k_series)
        }
        ConductivityMixingModel::Geometric => k_geom,
        ConductivityMixingModel::KopelmanIsotropic => k_kop,
        ConductivityMixingModel::MaxwellEucken => k_me,
    }
}

fn validate_eval_result(result: &ProximateEvalResult) -> Result<(), ProximateErrorV1> {
    for (property, value, require_positive) in [
        ("density_condensed", result.density_condensed_kg_m3, true),
        ("density_bulk", result.density_bulk_kg_m3, true),
        ("cp_sensible", result.cp_sensible_j_kg_k, true),
        ("cp_apparent", result.cp_apparent_j_kg_k, true),
        ("latent_heat_fusion", result.latent_heat_fusion_j_kg, false),
        ("enthalpy_0c_ref", result.enthalpy_0c_ref_j_kg, false),
        ("enthalpy_janaf", result.enthalpy_janaf_j_kg, false),
        ("k_parallel", result.k_parallel_w_m_k, true),
        ("k_series", result.k_series_w_m_k, true),
        ("k_maxwell_eucken", result.k_maxwell_eucken_w_m_k, true),
        ("k_kopelman", result.k_kopelman_w_m_k, true),
        ("k_effective", result.k_effective_w_m_k, true),
        (
            "thermal_diffusivity_sensible",
            result.thermal_diffusivity_sensible_m2_s,
            true,
        ),
        (
            "thermal_diffusivity_apparent",
            result.thermal_diffusivity_apparent_m2_s,
            true,
        ),
    ] {
        if !value.is_finite() || (require_positive && value <= 0.0) {
            return Err(ProximateErrorV1::NonFiniteEvaluation {
                property,
                bits: value.to_bits(),
            });
        }
    }
    Ok(())
}

fn ensure_positive_component(
    component: &'static str,
    property: &'static str,
    value: f64,
) -> Result<(), ProximateErrorV1> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ProximateErrorV1::InvalidCoefficientTable {
            component,
            property,
            reason: "evaluated component property must be finite and strictly positive",
        });
    }
    Ok(())
}

/// Deterministic cube root on `x >= 0` using `fs_math::det::pow` refined by
/// two IEEE-754 Newton steps (`x_{n+1} = (2 x_n + a / x_n^2) / 3`).
fn deterministic_cbrt_nonnegative(a: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    if a.to_bits() == 1.0_f64.to_bits() {
        return 1.0;
    }
    let mut x = det::pow(a, 1.0 / 3.0);
    if x > 0.0 && x.is_finite() {
        x = (2.0 * x + a / (x * x)) / 3.0;
        x = (2.0 * x + a / (x * x)) / 3.0;
    }
    x
}

fn maxwell_eucken_two_phase(k_continuous: f64, k_dispersed: f64, eps_dispersed: f64) -> f64 {
    if eps_dispersed <= 0.0 {
        return k_continuous;
    }
    let diff = k_continuous - k_dispersed;
    let base = 2.0 * k_continuous + k_dispersed;
    k_continuous * (base - 2.0 * diff * eps_dispersed) / (base + diff * eps_dispersed)
}

fn kopelman_two_phase(k_continuous: f64, k_dispersed: f64, eps_dispersed: f64) -> f64 {
    if eps_dispersed <= 0.0 {
        return k_continuous;
    }
    if eps_dispersed >= 1.0 {
        return k_dispersed;
    }
    let m = deterministic_cbrt_nonnegative(eps_dispersed);
    let q = (m * m) * (1.0 - k_dispersed / k_continuous);
    k_continuous * (1.0 - q) / (1.0 - q * (1.0 - m))
}

fn compute_conductivities(
    eps_cond: &[f64; 7],
    ks: &[f64; 7],
    phi_g: f64,
    k_gas: f64,
) -> (f64, f64, f64, f64, f64) {
    let solid_scale = 1.0 - phi_g;

    let mut k_par_cond = 0.0;
    let mut inv_k_ser_cond = 0.0;
    let mut ln_k_geom_cond = 0.0;
    let mut dominant_idx = 0usize;
    let mut dominant_eps = -1.0;

    for idx in 0..7 {
        let eps = eps_cond[idx];
        if eps > 0.0 {
            let k = ks[idx];
            k_par_cond += eps * k;
            inv_k_ser_cond += eps / k;
            ln_k_geom_cond += eps * det::ln(k);
            if eps > dominant_eps {
                dominant_eps = eps;
                dominant_idx = idx;
            }
        }
    }

    // Stage 1: Condensed-matrix Maxwell-Eucken (dominant condensed phase continuous).
    let k_c = ks[dominant_idx];
    let mut me_num = 0.0;
    let mut me_den = 0.0;
    let mut eps_disp = 0.0;
    let mut k_disp_weighted = 0.0;
    for idx in 0..7 {
        let eps = eps_cond[idx];
        if eps > 0.0 {
            let k = ks[idx];
            let weight = eps / (2.0 * k_c + k);
            me_num += weight * k;
            me_den += weight;
            if idx != dominant_idx {
                eps_disp += eps;
                k_disp_weighted += eps * k;
            }
        }
    }
    let k_me_cond = me_num / me_den;
    let k_kop_cond = if eps_disp > 0.0 {
        let k_disp = k_disp_weighted / eps_disp;
        kopelman_two_phase(k_c, k_disp, eps_disp)
    } else {
        k_c
    };

    // Stage 2: Incorporate optional gas-pore void fraction phi_g (dispersed inside condensed matrix).
    let (k_parallel, k_series, k_me, k_kop, k_geom) = if phi_g > 0.0 {
        let k_par = solid_scale * k_par_cond + phi_g * k_gas;
        let k_ser = 1.0 / (solid_scale * inv_k_ser_cond + phi_g / k_gas);
        let k_g_mean = det::exp(solid_scale * ln_k_geom_cond + phi_g * det::ln(k_gas));
        let k_me_porous = maxwell_eucken_two_phase(k_me_cond, k_gas, phi_g);
        let k_kop_porous = kopelman_two_phase(k_kop_cond, k_gas, phi_g);
        (k_par, k_ser, k_me_porous, k_kop_porous, k_g_mean)
    } else {
        (
            k_par_cond,
            1.0 / inv_k_ser_cond,
            k_me_cond,
            k_kop_cond,
            det::exp(ln_k_geom_cond),
        )
    };

    (k_parallel, k_series, k_me, k_kop, k_geom)
}

/// Evaluate IAPWS-IF97 liquid water enthalpy at the JANAF standard reference state
/// `(298.15 K, 1 bar = 100_000 Pa)` [J/kg].
pub fn water_if97_reference_enthalpy_298k_1bar()
-> Result<MassSpecificEnergyQuantityV1, ProximateErrorV1> {
    let state = If97State::at_pressure_temperature(
        Pressure::new(JANAF_REFERENCE_PRESSURE_PA),
        Temperature::new(JANAF_REFERENCE_TEMPERATURE_K),
        If97Region::Region1,
    )
    .map_err(|_| ProximateErrorV1::InvalidIf97Reference {
        reason: "failed to evaluate IF97 Region 1 liquid water at (298.15 K, 1 bar)",
    })?;
    Ok(state.specific_enthalpy())
}

/// Constant datum offset `h_f0(298.15 K) - h_IF97(298.15 K, 1 bar)` [J/kg]
/// shifting IAPWS-IF97 water/steam enthalpy onto the shared
/// `janaf-298.15K-formation` absolute enthalpy datum.
pub fn water_if97_janaf_offset() -> Result<MassSpecificEnergyQuantityV1, ProximateErrorV1> {
    let h_ref = water_if97_reference_enthalpy_298k_1bar()?.value();
    let offset = H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG - h_ref;
    if !offset.is_finite() {
        return Err(ProximateErrorV1::InvalidIf97Reference {
            reason: "IF97-to-JANAF datum offset is not finite",
        });
    }
    Ok(MassSpecificEnergyQuantityV1::new(offset))
}

/// Convert an IAPWS-IF97 water or steam specific enthalpy `h_IF97` [J/kg] onto
/// the shared `janaf-298.15K-formation` datum:
/// `h_H2O_JANAF = h_IF97 + h_f0(298.15 K) - h_IF97(298.15 K, 1 bar)`.
pub fn water_enthalpy_janaf_from_if97(
    if97_enthalpy: MassSpecificEnergyQuantityV1,
) -> Result<MassSpecificEnergyQuantityV1, ProximateErrorV1> {
    let h_if97 = if97_enthalpy.value();
    if !h_if97.is_finite() {
        return Err(ProximateErrorV1::InvalidIf97Reference {
            reason: "IF97 specific enthalpy input must be finite",
        });
    }
    let offset = water_if97_janaf_offset()?.value();
    let shifted = h_if97 + offset;
    if !shifted.is_finite() {
        return Err(ProximateErrorV1::NonFiniteEvaluation {
            property: "water_enthalpy_janaf",
            bits: shifted.to_bits(),
        });
    }
    Ok(MassSpecificEnergyQuantityV1::new(shifted))
}

/// Convert an admitted [`If97State`] specific enthalpy onto the shared
/// `janaf-298.15K-formation` datum.
pub fn water_state_enthalpy_janaf(
    state: &If97State,
) -> Result<MassSpecificEnergyQuantityV1, ProximateErrorV1> {
    water_enthalpy_janaf_from_if97(state.specific_enthalpy())
}
