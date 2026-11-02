//! Source-declared batch/contact-time recovery with explicit rate populations.
//!
//! These empirical curves retain the original constitutive card, source hashes
//! and typed applicability coordinates. Arithmetic verification does not turn
//! an authored or fitted curve into a calibrated plant forecast. Time is actual
//! contact time; a tank's hydraulic residence-time distribution is separate.
//!
//! The single and uniform-rate first-order forms are listed in *Minerals*
//! 10 (2020), 1027, Table 1, <https://doi.org/10.3390/min10111027>. This is a
//! model-family reference; it does not supply a rate or recovery for an ore.

#[path = "recovery/residence.rs"]
pub mod residence;

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{Dimensionless, Dims, Frequency, QuantitySpec, Time};

/// Constant-rate recoverable population: R(t)=Rinf*(1-exp(-k*t)).
pub const SINGLE_RATE_LAW: &str = "single-rate-recovery";
/// Uniform rate population on [0,kmax]; R(t)=Rinf*(1-(1-exp(-kmax*t))/(kmax*t)).
pub const UNIFORM_RATE_LAW: &str = "uniform-rate-recovery";
/// Typed seconds axis supplied by [`RecoveryKinetics::evaluate`].
pub const ELAPSED_TIME_AXIS: &str = "elapsed-time";
/// Dimensionless fraction of feed metal in the recoverable population.
pub const ULTIMATE_RECOVERY_PARAMETER: &str = "ultimate-recovery";
/// Single rate, in inverse seconds.
pub const RATE_PARAMETER: &str = "rate-constant";
/// Uniform distribution's upper endpoint, in inverse seconds.
pub const MAXIMUM_RATE_PARAMETER: &str = "maximum-rate-constant";

/// Unsupported source declaration, applicability, or numerical range.
#[derive(Clone, Debug, PartialEq)]
pub enum RecoveryError {
    /// Original material-card validation or typed query admission failed.
    Source(MatDbError),
    /// Invalid law, dimensions, source support, or parameter bounds.
    Invalid(&'static str),
    /// At least one required coordinate is missing, mistyped or outside support.
    OutsideValidity,
    /// A positive response cannot be represented without numerical range loss.
    Arithmetic(&'static str),
}
impl core::fmt::Display for RecoveryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "recovery kinetics: {self:?}")
    }
}
impl std::error::Error for RecoveryError {}
impl From<MatDbError> for RecoveryError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}

/// Explicit choice of rate population; neither choice is inferred from ore names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecoveryRateDistribution {
    /// All recoverable material has the supplied rate k.
    Single,
    /// Recoverable material has rates uniformly distributed from zero to kmax.
    UniformFromZero,
}

/// Immutable, source-bound empirical contact-time curve.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveryKinetics {
    source: ConstitutiveModelCard,
    distribution: RecoveryRateDistribution,
    ultimate: Dimensionless,
    rate: Frequency,
}

impl RecoveryKinetics {
    /// Admit exactly two finite dimensioned parameters, law version 1, no hidden
    /// state, nonempty original source hashes, and a finite typed time window.
    /// Additional typed validity axes (e.g. grind, pH or temperature) are retained.
    pub fn from_card(source: ConstitutiveModelCard) -> Result<Self, RecoveryError> {
        source.validate()?;
        let (distribution, rate_name) = match source.law.0.as_str() {
            SINGLE_RATE_LAW => (RecoveryRateDistribution::Single, RATE_PARAMETER),
            UNIFORM_RATE_LAW => (
                RecoveryRateDistribution::UniformFromZero,
                MAXIMUM_RATE_PARAMETER,
            ),
            _ => return Err(RecoveryError::Invalid("unsupported recovery law")),
        };
        if source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 2
            || source.sources.is_empty()
        {
            return Err(RecoveryError::Invalid(
                "requires version 1, zero state, two parameters and original source hashes",
            ));
        }
        let ultimate = source
            .parameters
            .get(ULTIMATE_RECOVERY_PARAMETER)
            .ok_or(RecoveryError::Invalid("missing ultimate-recovery"))?;
        let rate = source
            .parameters
            .get(rate_name)
            .ok_or(RecoveryError::Invalid(
                "missing rate parameter for declared distribution",
            ))?;
        if ultimate.dims != Dims::NONE
            || !(0.0..=1.0).contains(&ultimate.value)
            || rate.dims != Frequency::DIMS
            || rate.value < 0.0
        {
            return Err(RecoveryError::Invalid(
                "ultimate recovery must be a fraction and rate nonnegative inverse seconds",
            ));
        }
        let Some((lo, hi)) = source.validity.bound(ELAPSED_TIME_AXIS) else {
            return Err(RecoveryError::Invalid(
                "declare a finite elapsed-time validity window",
            ));
        };
        if !lo.is_finite()
            || !hi.is_finite()
            || lo < 0.0
            || hi <= lo
            || source.validity.axis_quantities().get(ELAPSED_TIME_AXIS)
                != Some(&QuantitySpec::dimensional(Time::DIMS))
            || source.validity.bounds().iter().any(|(axis, (lo, hi))| {
                !lo.is_finite()
                    || !hi.is_finite()
                    || lo > hi
                    || !source.validity.axis_quantities().contains_key(axis)
            })
        {
            return Err(RecoveryError::Invalid(
                "all applicability axes need finite bounds and explicit quantities; time must be nonnegative seconds",
            ));
        }
        let ultimate = Dimensionless::new(ultimate.value);
        let rate = Frequency::new(rate.value);
        Ok(Self {
            source,
            distribution,
            ultimate,
            rate,
        })
    }

    /// Complete original source, including covariance/fit evidence if its source
    /// artifact supplies it; this evaluator does not manufacture uncertainty.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }

    /// Declared population of rates.
    pub const fn distribution(&self) -> RecoveryRateDistribution {
        self.distribution
    }

    /// Maximum recoverable fraction, not an assertion of operational recovery.
    pub const fn ultimate_recovery(&self) -> Dimensionless {
        self.ultimate
    }

    /// k for a single population, or kmax for a uniform population.
    pub const fn rate(&self) -> Frequency {
        self.rate
    }

    /// Evaluate at actual contact time and all source conditions. The separate
    /// typed time argument supplies `elapsed-time`; conflicting duplicate time
    /// coordinates refuse. Unsupported conditions never silently extrapolate.
    pub fn evaluate(
        &self,
        elapsed: Time,
        conditions: &QueryPoint,
    ) -> Result<Dimensionless, RecoveryError> {
        let t = elapsed.value();
        if !t.is_finite() || t < 0.0 {
            return Err(RecoveryError::Invalid(
                "contact time must be finite and nonnegative",
            ));
        }
        if let Some(value) = conditions.axes().get(ELAPSED_TIME_AXIS) {
            if *value != t
                || conditions.axis_quantities().get(ELAPSED_TIME_AXIS)
                    != Some(&QuantitySpec::dimensional(Time::DIMS))
            {
                return Err(RecoveryError::Invalid(
                    "conflicting duplicate contact-time coordinate",
                ));
            }
        }
        let point = conditions.clone().with_quantity(
            ELAPSED_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            t,
        )?;
        if !self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            return Err(RecoveryError::OutsideValidity);
        }
        if t == 0.0 || self.rate.value() == 0.0 || self.ultimate.value() == 0.0 {
            return Ok(Dimensionless::new(0.0));
        }
        let x = self.rate.value() * t;
        if !x.is_finite() || x == 0.0 {
            return Err(RecoveryError::Arithmetic("rate-time product"));
        }
        let fraction = match self.distribution {
            RecoveryRateDistribution::Single => -det::expm1(-x),
            RecoveryRateDistribution::UniformFromZero if x < 0.125 => {
                // Integral of (1-exp(-k*t)) over k/kmax in [0,1]. Positive
                // leading term avoids cancellation of 1-(1-exp(-x))/x.
                let mut term = x / 2.0;
                let mut sum = term;
                for n in 2..=12 {
                    term *= -x / f64::from(n + 1);
                    sum += term;
                }
                sum
            }
            RecoveryRateDistribution::UniformFromZero => 1.0 + det::expm1(-x) / x,
        };
        let value = self.ultimate.value() * fraction;
        if !value.is_finite() || value <= 0.0 || value > self.ultimate.value() {
            return Err(RecoveryError::Arithmetic("recovery fraction"));
        }
        Ok(Dimensionless::new(value))
    }
}

#[cfg(test)]
#[path = "recovery/tests.rs"]
mod tests;
