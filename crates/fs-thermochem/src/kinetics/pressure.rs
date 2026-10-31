//! Authored collider and falloff laws. Their common multiplier preserves the
//! original thermodynamic equilibrium constant in both reaction directions.
//! Equations follow Cantera 3.2's third-body and Lindemann/Troe documentation:
//! <https://cantera.org/3.2/reference/kinetics/rate-constants.html>.

use super::{
    ArrheniusRateV1, BoundReaction, ElementaryReactionV1, KineticsErrorV1, SpeciesId, Sum, filled,
    positive, rate_from_log,
};
use crate::equilibrium::composition::MAX_EQUILIBRIUM_SPECIES_V1;
use fs_qty::Temperature;

/// Immutable, temperature-independent collider efficiencies.
///
/// The effective concentration is `M/c* = sum(epsilon_i c_i/c*)`, with the same
/// `c*=1 mol/m3` normalization as the elementary rate. Unlisted species use the
/// explicit default. Zero efficiencies exclude a species without removing it
/// from the thermodynamic bank or from elemental accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct ThirdBodyEfficienciesV1 {
    default_efficiency: f64,
    overrides: Vec<(SpeciesId, f64)>,
}

impl ThirdBodyEfficienciesV1 {
    /// Admit finite nonnegative efficiencies, at most 128 unique overrides.
    /// The mechanism subsequently checks every named species against its bank.
    pub fn new(
        default_efficiency: f64,
        mut overrides: Vec<(SpeciesId, f64)>,
    ) -> Result<Self, KineticsErrorV1> {
        if !default_efficiency.is_finite()
            || default_efficiency < 0.0
            || overrides.len() > MAX_EQUILIBRIUM_SPECIES_V1
            || overrides
                .iter()
                .any(|(_, value)| !value.is_finite() || *value < 0.0)
        {
            return Err(KineticsErrorV1::Invalid("invalid third-body efficiencies"));
        }
        overrides.sort_by(|a, b| a.0.cmp(&b.0));
        if overrides.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(KineticsErrorV1::Invalid("duplicate third-body species"));
        }
        Ok(Self {
            default_efficiency,
            overrides,
        })
    }
    /// Authored efficiency for every species not explicitly overridden.
    pub const fn default_efficiency(&self) -> f64 {
        self.default_efficiency
    }
    /// Canonical original species-specific overrides, including exact zeros.
    pub fn overrides(&self) -> &[(SpeciesId, f64)] {
        &self.overrides
    }

    fn bind(&self, species: &[SpeciesId]) -> Result<Vec<f64>, KineticsErrorV1> {
        let log = |value| {
            if value == 0.0 {
                f64::NEG_INFINITY
            } else {
                fs_math::det::ln(value)
            }
        };
        let mut logs = filled(species.len(), log(self.default_efficiency))?;
        for (id, efficiency) in &self.overrides {
            let index = species.binary_search(id).map_err(|_| {
                KineticsErrorV1::Invalid("third-body species is absent from the source bank")
            })?;
            logs[index] = log(*efficiency);
        }
        Ok(logs)
    }
}

impl Default for ThirdBodyEfficienciesV1 {
    /// The explicitly selected conventional unit efficiency for every species.
    fn default() -> Self {
        Self {
            default_efficiency: 1.0,
            overrides: Vec::new(),
        }
    }
}

/// Three- or four-parameter Troe broadening with source temperatures in kelvin.
///
/// `Fcent=(1-alpha) exp(-T/T3)+alpha exp(-T/T1)+exp(-T2/T)`.
/// The final term is absent when `temperature_2` is `None`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TroeFalloffV1 {
    alpha: f64,
    temperature_3: Temperature,
    temperature_1: Temperature,
    temperature_2: Option<Temperature>,
}

impl TroeFalloffV1 {
    /// Admit `alpha` in [0,1] and strictly positive finite source temperatures.
    /// These parameters are caller-authored kinetic data, not a fitted mechanism.
    pub fn new(
        alpha: f64,
        temperature_3: Temperature,
        temperature_1: Temperature,
        temperature_2: Option<Temperature>,
    ) -> Result<Self, KineticsErrorV1> {
        if !alpha.is_finite()
            || !(0.0..=1.0).contains(&alpha)
            || !positive(temperature_3.value())
            || !positive(temperature_1.value())
            || temperature_2.is_some_and(|t| !positive(t.value()))
        {
            return Err(KineticsErrorV1::Invalid("invalid Troe falloff parameters"));
        }
        Ok(Self {
            alpha,
            temperature_3,
            temperature_1,
            temperature_2,
        })
    }
    /// Original dimensionless central weighting parameter.
    pub const fn alpha(&self) -> f64 {
        self.alpha
    }
    /// Original T3 [K].
    pub const fn temperature_3(&self) -> Temperature {
        self.temperature_3
    }
    /// Original T1 [K].
    pub const fn temperature_1(&self) -> Temperature {
        self.temperature_1
    }
    /// Original optional T2 [K].
    pub const fn temperature_2(&self) -> Option<Temperature> {
        self.temperature_2
    }

    fn center(self, t: f64) -> Result<(f64, f64), KineticsErrorV1> {
        let logs = [
            if self.alpha == 1.0 {
                f64::NEG_INFINITY
            } else {
                fs_math::det::ln(1.0 - self.alpha) - t / self.temperature_3.value()
            },
            if self.alpha == 0.0 {
                f64::NEG_INFINITY
            } else {
                fs_math::det::ln(self.alpha) - t / self.temperature_1.value()
            },
            self.temperature_2
                .map_or(f64::NEG_INFINITY, |t2| -t2.value() / t),
        ];
        let maximum = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !maximum.is_finite() {
            return Err(KineticsErrorV1::Arithmetic("Troe central logarithm"));
        }
        let slopes = [
            -1.0 / self.temperature_3.value(),
            -1.0 / self.temperature_1.value(),
            self.temperature_2.map_or(0.0, |t2| (t2.value() / t) / t),
        ];
        let mut sum = Sum::default();
        let mut derivative = Sum::default();
        for i in 0..3 {
            let weight = fs_math::det::exp(logs[i] - maximum);
            sum.add(weight);
            // A negligible exponential must not create 0*infinity from an
            // irrelevant, extremely small authored temperature parameter.
            if weight != 0.0 {
                derivative.add(weight * slopes[i]);
            }
        }
        let log = maximum + fs_math::det::ln(sum.value());
        let dlog = derivative.value() / sum.value();
        if !log.is_finite() || !dlog.is_finite() {
            return Err(KineticsErrorV1::Arithmetic("Troe central derivative"));
        }
        Ok((log, dlog))
    }

    /// ln F, d ln F/d ln Pr, and d ln F/dT at fixed Pr.
    fn evaluate(self, t: f64, log_pr: f64) -> Result<(f64, f64, f64), KineticsErrorV1> {
        let (center, center_t) = self.center(t)?;
        const LN10: f64 = core::f64::consts::LN_10;
        if log_pr == f64::NEG_INFINITY {
            // The exact Troe correlation has this finite limiting broadening,
            // not F=1. Its value determines the derivative at exactly M=0.
            const WEIGHT: f64 = (0.14 * 0.14) / (1.0 + 0.14 * 0.14);
            return Ok((center * WEIGHT, 0.0, center_t * WEIGHT));
        }
        let c = -0.4 - 0.67 * (center / LN10);
        let n = 0.75 - 1.27 * (center / LN10);
        let y = log_pr / LN10 + c;
        let d = n - 0.14 * y;
        let scale = y.abs().max(d.abs()).max(n.abs());
        if !positive(scale) {
            return Err(KineticsErrorV1::Arithmetic("Troe reduced-pressure shape"));
        }
        let ys = y / scale;
        let ds = d / scale;
        let ns = n / scale;
        let norm = ys * ys + ds * ds;
        let weight = ds * ds / norm;
        let factor = 2.0 * ds * ys / (norm * norm);
        // Evaluate d²/(d²+y²), not 1/(1+(y/d)²): the usual rational form has
        // an artificial pole at d=0 even though F and both derivatives are smooth.
        let log_f = center * weight;
        let derivative_pr = -(center / scale) * factor * ns / LN10;
        let derivative_t =
            center_t * (weight + (center / scale) * factor * (-1.27 * ys + 0.67 * ns) / LN10);
        if !log_f.is_finite() || !derivative_pr.is_finite() || !derivative_t.is_finite() {
            return Err(KineticsErrorV1::Arithmetic("Troe broadening derivative"));
        }
        Ok((log_f, derivative_pr, derivative_t))
    }
}

/// Authored broadening of the normal Lindemann falloff expression.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum FalloffBroadeningV1 {
    /// F=1: k=k0 M/(1+k0 M/kinfinity).
    #[default]
    Lindemann,
    /// Original three-/four-parameter Troe correlation.
    Troe(TroeFalloffV1),
}

/// Optional pressure law for one reversible elementary stoichiometric column.
///
/// Both forward and reverse rates receive the same positive kinetic multiplier.
/// The collider does not alter stoichiometry, Kc, formation energy or affinity.
/// Chemically activated, P-log, SRI and temperature-dependent efficiencies are
/// separate unsupported laws; none is silently approximated by these variants.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum PressureDependenceV1 {
    /// The existing elementary, pressure-independent Arrhenius law.
    #[default]
    Independent,
    /// Multiply the existing forward law by the effective collider concentration.
    ThirdBody(ThirdBodyEfficienciesV1),
    /// Normal falloff. The elementary forward law supplies `kinfinity`.
    /// `low_pressure_rate` supplies k0 normalized by one additional c*=1mol/m³.
    /// Either exact zero rate scale disables both directions and all derivatives.
    Falloff {
        /// Original low-pressure modified Arrhenius scale.
        low_pressure_rate: ArrheniusRateV1,
        /// Original collider efficiencies.
        colliders: ThirdBodyEfficienciesV1,
        /// F(T,Pr), shared by both reaction directions.
        broadening: FalloffBroadeningV1,
    },
}

impl PressureDependenceV1 {
    pub(super) fn bind(&self, species: &[SpeciesId]) -> Result<Vec<f64>, KineticsErrorV1> {
        match self {
            Self::Independent => Ok(Vec::new()),
            Self::ThirdBody(efficiencies)
            | Self::Falloff {
                colliders: efficiencies,
                ..
            } => efficiencies.bind(species),
        }
    }
}

/// A derivative may be negative; logarithmic magnitude avoids losing a tiny
/// slope before a large source concentration or efficiency rescales it.
#[derive(Clone, Copy)]
pub(super) struct SignedLog {
    sign: f64,
    log_abs: f64,
}

impl SignedLog {
    const ZERO: Self = Self {
        sign: 0.0,
        log_abs: f64::NEG_INFINITY,
    };
    fn new(sign: f64, log_abs: f64) -> Result<Self, KineticsErrorV1> {
        if !log_abs.is_finite() {
            return Err(KineticsErrorV1::Arithmetic(
                "pressure-rate derivative logarithm",
            ));
        }
        Ok(Self { sign, log_abs })
    }
    fn add(self, other: Self) -> Result<Self, KineticsErrorV1> {
        if self.sign == 0.0 {
            return Ok(other);
        }
        if other.sign == 0.0 {
            return Ok(self);
        }
        let (larger, smaller) = if self.log_abs >= other.log_abs {
            (self, other)
        } else {
            (other, self)
        };
        let difference = smaller.log_abs - larger.log_abs;
        if larger.sign == smaller.sign {
            Self::new(
                larger.sign,
                larger.log_abs + log_one_plus(fs_math::det::exp(difference)),
            )
        } else if difference == 0.0 {
            Ok(Self::ZERO)
        } else {
            Self::new(
                larger.sign,
                larger.log_abs + fs_math::det::ln(-fs_math::det::expm1(difference)),
            )
        }
    }
    pub(super) fn scaled(self, log_efficiency: f64) -> Result<f64, KineticsErrorV1> {
        if self.sign == 0.0 || log_efficiency == f64::NEG_INFINITY {
            return Ok(0.0);
        }
        Ok(self.sign * rate_from_log(self.log_abs + log_efficiency)?)
    }
}

pub(super) struct RateScale {
    pub(super) log_forward: f64,
    pub(super) dlog_temperature: f64,
    /// Signed d k_effective / d(M/c*), before species mass-action products.
    collider: SignedLog,
}

impl RateScale {
    fn zero() -> Self {
        Self {
            log_forward: f64::NEG_INFINITY,
            dlog_temperature: 0.0,
            collider: SignedLog::ZERO,
        }
    }
    pub(super) fn inactive(&self) -> bool {
        self.log_forward == f64::NEG_INFINITY && self.collider.sign == 0.0
    }
    /// Collider contribution to net progress, before each species efficiency.
    /// Using the thermodynamic affinity preserves near-equilibrium differences.
    pub(super) fn collider_response(
        &self,
        forward_product: f64,
        reverse_product: f64,
        log_kc: f64,
        affinity: f64,
    ) -> Result<SignedLog, KineticsErrorV1> {
        if self.collider.sign == 0.0
            || (forward_product == f64::NEG_INFINITY && reverse_product == f64::NEG_INFINITY)
        {
            return Ok(SignedLog::ZERO);
        }
        let (sign, log_product) = if forward_product == f64::NEG_INFINITY {
            (-1.0, reverse_product - log_kc)
        } else if reverse_product == f64::NEG_INFINITY {
            (1.0, forward_product)
        } else if affinity > 0.0 {
            (
                1.0,
                forward_product + fs_math::det::ln(-fs_math::det::expm1(-affinity)),
            )
        } else if affinity < 0.0 {
            (
                -1.0,
                reverse_product - log_kc + fs_math::det::ln(-fs_math::det::expm1(affinity)),
            )
        } else {
            return Ok(SignedLog::ZERO);
        };
        SignedLog::new(
            sign * self.collider.sign,
            self.collider.log_abs + log_product,
        )
    }
}

fn collider_log(log_c: &[f64], log_efficiencies: &[f64]) -> f64 {
    let maximum = log_c
        .iter()
        .zip(log_efficiencies)
        .map(|(&c, &efficiency)| c + efficiency)
        .fold(f64::NEG_INFINITY, f64::max);
    if maximum == f64::NEG_INFINITY {
        return maximum;
    }
    let mut sum = Sum::default();
    for (&c, &efficiency) in log_c.iter().zip(log_efficiencies) {
        sum.add(fs_math::det::exp(c + efficiency - maximum));
    }
    maximum + fs_math::det::ln(sum.value())
}

/// Accurate ln(1+x) for 0<=x<=1, using the deterministic logarithm owner.
fn log_one_plus(x: f64) -> f64 {
    let sum = 1.0 + x;
    if sum == 1.0 {
        x
    } else {
        fs_math::det::ln(sum) * (x / (sum - 1.0))
    }
}

pub(super) fn effective_rate(
    source: &ElementaryReactionV1,
    reaction: &BoundReaction,
    temperature: f64,
    log_c: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RateScale, KineticsErrorV1> {
    // Exact authored shutdown takes precedence over unused limiting-rate
    // arithmetic. In particular, k0=0 must not evaluate an irrelevant kinfinity
    // law that cannot be represented at this temperature.
    if source.forward_rate.disabled()
        || matches!(&source.pressure_dependence,
            PressureDependenceV1::Falloff { low_pressure_rate, .. }
                if low_pressure_rate.reference_rate().value() == 0.0)
        || (!reaction.collider_logs.is_empty()
            && reaction
                .collider_logs
                .iter()
                .all(|log| *log == f64::NEG_INFINITY))
    {
        return Ok(RateScale::zero());
    }
    let (high, high_t) = source.forward_rate.log_rate(temperature, cancelled)?;
    if high == f64::NEG_INFINITY {
        return Ok(RateScale::zero());
    }
    if source.pressure_dependence == PressureDependenceV1::Independent {
        return Ok(RateScale {
            log_forward: high,
            dlog_temperature: high_t,
            collider: SignedLog::ZERO,
        });
    }
    let log_m = collider_log(log_c, &reaction.collider_logs);
    match &source.pressure_dependence {
        PressureDependenceV1::Independent => unreachable!(),
        PressureDependenceV1::ThirdBody(_) => {
            let log_forward = high + log_m;
            if log_forward != f64::NEG_INFINITY && !log_forward.is_finite() {
                return Err(KineticsErrorV1::Arithmetic("third-body rate logarithm"));
            }
            Ok(RateScale {
                log_forward,
                dlog_temperature: high_t,
                collider: SignedLog::new(1.0, high)?,
            })
        }
        PressureDependenceV1::Falloff {
            low_pressure_rate,
            broadening,
            ..
        } => {
            let (low, low_t) = low_pressure_rate.log_rate(temperature)?;
            if low == f64::NEG_INFINITY {
                return Ok(RateScale::zero());
            }
            let log_pr = low + log_m - high;
            if log_pr != f64::NEG_INFINITY && !log_pr.is_finite() {
                return Err(KineticsErrorV1::Arithmetic("falloff reduced pressure"));
            }
            let (log_f, f_pr, f_t) = match broadening {
                FalloffBroadeningV1::Lindemann => (0.0, 0.0, 0.0),
                FalloffBroadeningV1::Troe(troe) => troe.evaluate(temperature, log_pr)?,
            };
            if log_m == f64::NEG_INFINITY {
                return Ok(RateScale {
                    log_forward: f64::NEG_INFINITY,
                    dlog_temperature: 0.0,
                    collider: SignedLog::new(1.0, low + log_f)?,
                });
            }
            // Neither Pr nor k0*M is exponentiated. Extreme low/high-pressure
            // limits therefore remain usable when those intermediates overflow.
            let small = fs_math::det::exp(-log_pr.abs());
            let soft = log_one_plus(small);
            let log_one_plus_pr = log_pr.max(0.0) + soft;
            let low_weight = if log_pr >= 0.0 {
                small / (1.0 + small)
            } else {
                1.0 / (1.0 + small)
            };
            let high_weight = if log_pr >= 0.0 {
                1.0 / (1.0 + small)
            } else {
                small / (1.0 + small)
            };
            let log_forward = if log_pr >= 0.0 {
                high - soft + log_f
            } else {
                low + log_m - soft + log_f
            };
            let dlog_temperature =
                low_t * low_weight + high_t * high_weight + f_pr * (low_t - high_t) + f_t;
            let lindemann_derivative = SignedLog::new(1.0, low + log_f - 2.0 * log_one_plus_pr)?;
            let broadening_derivative = if f_pr == 0.0 {
                SignedLog::ZERO
            } else {
                SignedLog::new(
                    f_pr.signum(),
                    log_forward - log_m + fs_math::det::ln(f_pr.abs()),
                )?
            };
            if !log_forward.is_finite() || !dlog_temperature.is_finite() {
                return Err(KineticsErrorV1::Arithmetic(
                    "falloff rate or temperature derivative",
                ));
            }
            Ok(RateScale {
                log_forward,
                dlog_temperature,
                collider: lindemann_derivative.add(broadening_derivative)?,
            })
        }
    }
}
