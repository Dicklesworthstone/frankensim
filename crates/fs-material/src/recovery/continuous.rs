//! Source-declared steady recovery in equal, continuously stirred tanks.
//!
//! A stationary local first-order rate is a different claim from a batch curve
//! tested over a finite contact-time window. These distinct laws admit their
//! own original cards and operating domains; no batch-card conversion exists.
//! In particular, they never evaluate batch recovery at mean residence time or
//! silently extend a finite elapsed-time validity window over an infinite RTD.
//!
//! Every tank is full and perfectly mixed, with the same actual carrier flow,
//! active volume and constant rate law. The recoverable feed population has
//! one rate, uniform rates on [0,kmax], or gamma-distributed rates. Each rate class
//! retains its rate as it passes through the train; the surviving population is
//! not reset to the feed distribution after each tank. The remainder of feed
//! metal is inert under this law. Separation does not alter the carrier flow.
//!
//! For one rate, steady mass balances give survival `(1+k*V/(N*Q))^(-N)`.
//! See Fogler, *Elements of Chemical Reaction Engineering*, §17.4.2, eq.5-15:
//! <https://websites.umich.edu/~elements/5e/17chap/Ch17_Web-Additional%20Material.pdf>.
//! The uniform case analytically averages that survival over the declared rate
//! population. Source support and arithmetic verification do not establish plant
//! calibration. Transients, adsorption, reagent limits and volume changes are
//! outside this steady, dilute-constituent model.
//!
//! Gamma responses integrate the rate-distribution Laplace transform against
//! the actual Erlang residence distribution, with bounded work, explicit tail
//! bounds and quadrature estimates. They do not evaluate batch recovery at the
//! mean residence. Shape, inverse-second scale and original local-law support
//! are mandatory; no batch-to-continuous calibration is inferred.

use super::{ELAPSED_TIME_AXIS, RecoveryError, RecoveryRateDistribution, residence::MAX_TANKS};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_math::det;
use fs_qty::{
    Dimensionless, Dims, Frequency, MassFlowRate, QuantitySpec, Time, Volume, VolumetricFlowRate,
};

pub use super::gamma::GammaIntegrationReport;

/// Stationary local single-rate recovery; not the finite-time batch law.
pub const CONTINUOUS_SINGLE_RATE_LAW: &str = "cstr-single-rate-recovery";
/// Stationary local uniform rate population; not the finite-time batch law.
pub const CONTINUOUS_UNIFORM_RATE_LAW: &str = "cstr-uniform-rate-recovery";
/// Stationary local gamma rate population; not the finite-time batch law.
pub const CONTINUOUS_GAMMA_RATE_LAW: &str = "cstr-gamma-rate-recovery";
/// Operating-domain coordinate V/Q, in seconds; not an individual particle age.
pub const MEAN_RESIDENCE_TIME_AXIS: &str = "mean-residence-time";
/// Dimensionless integer count of equal ideal tanks in the declared model.
pub const TANK_COUNT_AXIS: &str = "tank-count";

/// A separate source declaration for stationary local first-order treatment.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousRecoveryCard {
    source: ConstitutiveModelCard,
    distribution: RecoveryRateDistribution,
    ultimate: f64,
    rate: f64,
    shape: Option<f64>,
}

/// A steady constituent balance with the complete original source and conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousTreatment {
    /// Complete original card; no uncertainty or calibration evidence is added.
    pub source_card: ConstitutiveModelCard,
    /// Actual typed coordinates, including derived mean residence and tank count.
    pub conditions: QueryPoint,
    /// Equal full tanks represented by the model.
    pub tank_count: usize,
    /// Total active mixed volume of the train, in cubic metres.
    pub active_volume: Volume,
    /// Actual flow through the train, in cubic metres per second.
    pub actual_train_flow: VolumetricFlowRate,
    /// Hydraulic mean residence V/Q, in seconds.
    pub mean_residence_time: Time,
    /// Fraction of entering constituent recovered by this steady model.
    pub recovery: Dimensionless,
    /// Independently evaluated inert plus surviving recoverable fraction.
    pub unrecovered_fraction: Dimensionless,
    /// Entering contained constituent mass, in kg/s; not bulk ore throughput.
    pub feed_metal_flow: MassFlowRate,
    /// Recovered contained constituent mass, in kg/s.
    pub recovered_metal_flow: MassFlowRate,
    /// Unrecovered contained constituent mass leaving the train, in kg/s.
    pub tail_metal_flow: MassFlowRate,
    /// Recovered plus tail minus feed, in kg/s; no residual balancing plug.
    pub conservation_residual: MassFlowRate,
    /// Gamma numerical work/error diagnostics; analytic single/uniform laws
    /// and explicitly zero-rate/zero-recoverable responses need no quadrature.
    pub gamma_integration: Option<GammaIntegrationReport>,
}

impl ContinuousRecoveryCard {
    /// Admit a distinct stationary law with fraction Rinf and inverse-second k
    /// (or kmax, or gamma scale plus dimensionless shape), version 1 and no
    /// hidden state. The source must explicitly bound
    /// positive hydraulic mean residence and integer tank count, and type every
    /// additional applicability axis. Any batch `elapsed-time` axis refuses.
    pub fn from_card(source: ConstitutiveModelCard) -> Result<Self, RecoveryError> {
        source.validate()?;
        let (distribution, rate_name) = match source.law.0.as_str() {
            CONTINUOUS_SINGLE_RATE_LAW => (RecoveryRateDistribution::Single, super::RATE_PARAMETER),
            CONTINUOUS_UNIFORM_RATE_LAW => (
                RecoveryRateDistribution::UniformFromZero,
                super::MAXIMUM_RATE_PARAMETER,
            ),
            CONTINUOUS_GAMMA_RATE_LAW => (
                RecoveryRateDistribution::Gamma,
                super::GAMMA_SCALE_PARAMETER,
            ),
            _ => {
                return Err(RecoveryError::Invalid(
                    "requires a separately sourced continuous CSTR law",
                ));
            }
        };
        if source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.sources.is_empty()
            || source.parameters.len()
                != if distribution == RecoveryRateDistribution::Gamma {
                    3
                } else {
                    2
                }
            || source.validity.bound(ELAPSED_TIME_AXIS).is_some()
        {
            return Err(RecoveryError::Invalid(
                "requires version 1, zero state, the exact sourced parameters and no batch elapsed-time domain",
            ));
        }
        let ultimate = source
            .parameters
            .get(super::ULTIMATE_RECOVERY_PARAMETER)
            .ok_or(RecoveryError::Invalid(
                "missing recoverable population fraction",
            ))?;
        let rate = source
            .parameters
            .get(rate_name)
            .ok_or(RecoveryError::Invalid("missing declared local rate"))?;
        if ultimate.dims != Dims::NONE
            || !(0.0..=1.0).contains(&ultimate.value)
            || rate.dims != Frequency::DIMS
            || rate.value < 0.0
        {
            return Err(RecoveryError::Invalid(
                "requires fraction Rinf in [0,1] and nonnegative local rate in inverse seconds",
            ));
        }
        let shape = super::gamma_shape(&source, distribution)?.map(Dimensionless::value);
        let validity = &source.validity;
        let mean = validity.bound(MEAN_RESIDENCE_TIME_AXIS);
        let count = validity.bound(TANK_COUNT_AXIS);
        if !mean.is_some_and(|(lo, hi)| lo > 0.0 && hi >= lo)
            || !count.is_some_and(|(lo, hi)| {
                lo >= 1.0
                    && hi <= MAX_TANKS as f64
                    && hi >= lo
                    && lo.fract() == 0.0
                    && hi.fract() == 0.0
            })
            || validity.axis_quantities().get(MEAN_RESIDENCE_TIME_AXIS)
                != Some(&QuantitySpec::dimensional(Time::DIMS))
            || validity.axis_quantities().get(TANK_COUNT_AXIS)
                != Some(&QuantitySpec::dimensional(Dims::NONE))
            || validity.bounds().iter().any(|(axis, (lo, hi))| {
                !lo.is_finite()
                    || !hi.is_finite()
                    || lo > hi
                    || !validity.axis_quantities().contains_key(axis)
            })
        {
            return Err(RecoveryError::Invalid(
                "declare finite typed operating domains, including positive mean residence and integer tank count within 1..=64",
            ));
        }
        let (ultimate, rate) = (ultimate.value, rate.value);
        Ok(Self {
            source,
            distribution,
            ultimate,
            rate,
            shape,
        })
    }

    /// Complete original evidence for this stationary local-rate declaration.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Declared feed population; no distribution is inferred from ore identity.
    pub const fn distribution(&self) -> RecoveryRateDistribution {
        self.distribution
    }
    /// Recoverable fraction of entering contained constituent.
    pub fn ultimate_recovery(&self) -> Dimensionless {
        Dimensionless::new(self.ultimate)
    }
    /// Local single k or upper endpoint kmax, in inverse seconds.
    pub fn rate(&self) -> Frequency {
        Frequency::new(self.rate)
    }

    /// Evaluate a stationary train using actual active volume and actual flow.
    /// Zero flow refuses: there is no steady through-flow treatment during a
    /// stop. Zero contained-metal feed produces zero output rates. Complete
    /// conditions must match the card; conflicting duplicate derived axes refuse.
    pub fn evaluate_train(
        &self,
        tank_count: usize,
        active_volume: Volume,
        actual_train_flow: VolumetricFlowRate,
        feed_metal_flow: MassFlowRate,
        conditions: &QueryPoint,
    ) -> Result<ContinuousTreatment, RecoveryError> {
        let (volume, flow, feed) = (
            active_volume.value(),
            actual_train_flow.value(),
            feed_metal_flow.value(),
        );
        if !(1..=MAX_TANKS).contains(&tank_count)
            || !volume.is_finite()
            || volume <= 0.0
            || !flow.is_finite()
            || flow <= 0.0
            || !feed.is_finite()
            || feed < 0.0
        {
            return Err(RecoveryError::Invalid(
                "requires 1..=64 tanks, positive finite active volume/actual flow and nonnegative finite contained-metal flow",
            ));
        }
        let mean = volume / flow;
        if !mean.is_finite() || mean == 0.0 {
            return Err(RecoveryError::Arithmetic("hydraulic mean residence"));
        }
        let mut point = conditions.clone();
        for (name, quantity, value) in [
            (
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                mean,
            ),
            (
                TANK_COUNT_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                tank_count as f64,
            ),
        ] {
            if point.axes().get(name).is_some_and(|v| {
                *v != value || point.axis_quantities().get(name) != Some(&quantity)
            }) {
                return Err(RecoveryError::Invalid(
                    "conflicting derived steady operating coordinate",
                ));
            }
            point = point.with_quantity(name, quantity, value)?;
        }
        if !self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            return Err(RecoveryError::OutsideValidity);
        }
        let mut gamma_integration = None;
        let (converted, survival) = if self.rate == 0.0 || self.ultimate == 0.0 {
            (0.0, 1.0)
        } else if let Some(shape) = self.shape {
            let result = super::gamma::continuous(shape, self.rate, mean, tank_count)?;
            gamma_integration = Some(result.report);
            (result.converted, result.surviving)
        } else {
            let mut a = self.rate * (mean / tank_count as f64);
            if a == 0.0 {
                a = (self.rate * mean) / tank_count as f64;
            }
            if !a.is_finite() || a == 0.0 {
                return Err(RecoveryError::Arithmetic(
                    "local rate times individual tank residence",
                ));
            }
            fractions(self.distribution, tank_count, a)?
        };
        let recovery = checked_product(self.ultimate, converted, "recovery fraction")?;
        let tail = (1.0 - self.ultimate) + self.ultimate * survival;
        if !tail.is_finite() || tail <= 0.0 || tail > 1.0 || recovery > self.ultimate {
            return Err(RecoveryError::Arithmetic("steady output fractions"));
        }
        let recovered = checked_product(feed, recovery, "recovered constituent flow")?;
        let unrecovered = checked_product(feed, tail, "tail constituent flow")?;
        let residual = (recovered + unrecovered) - feed;
        if !residual.is_finite() || residual.abs() > 2e-13 * feed {
            return Err(RecoveryError::Arithmetic("steady constituent conservation"));
        }
        Ok(ContinuousTreatment {
            source_card: self.source.clone(),
            conditions: point,
            tank_count,
            active_volume,
            actual_train_flow,
            mean_residence_time: Time::new(mean),
            recovery: Dimensionless::new(recovery),
            unrecovered_fraction: Dimensionless::new(tail),
            feed_metal_flow,
            recovered_metal_flow: MassFlowRate::new(recovered),
            tail_metal_flow: MassFlowRate::new(unrecovered),
            conservation_residual: MassFlowRate::new(residual),
            gamma_integration,
        })
    }
}

fn checked_product(a: f64, b: f64, context: &'static str) -> Result<f64, RecoveryError> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(RecoveryError::Arithmetic(context))
    } else {
        Ok(value)
    }
}

// Compensate the rounded addition in ln(1+a); for sub-ULP positive a the
// leading term a is already correctly rounded. Uses the deterministic log.
fn log_one_plus(a: f64) -> f64 {
    let sum = 1.0 + a;
    if sum == 1.0 {
        a
    } else {
        det::ln(sum) * (a / (sum - 1.0))
    }
}

fn fractions(
    distribution: RecoveryRateDistribution,
    n: usize,
    a: f64,
) -> Result<(f64, f64), RecoveryError> {
    let logarithm = log_one_plus(a);
    let result = match distribution {
        RecoveryRateDistribution::Single => {
            let exponent = -(n as f64) * logarithm;
            (-det::expm1(exponent), det::exp(exponent))
        }
        RecoveryRateDistribution::UniformFromZero if n == 1 && a < 0.125 => {
            // 1-ln(1+a)/a = a/2-a²/3+a³/4-...; avoids cancellation.
            let mut power = a;
            let mut converted = 0.0;
            for j in 1..=20 {
                converted += power / f64::from(j + 1);
                power *= -a;
            }
            (converted, 1.0 - converted)
        }
        RecoveryRateDistribution::UniformFromZero if n == 1 => {
            let survival = logarithm / a;
            (1.0 - survival, survival)
        }
        RecoveryRateDistribution::UniformFromZero => {
            // ∫₀¹(1+a*u)^(-n)du = Σ_(j=1)^(n-1)(1+a)^(-j)/(n-1).
            // Average converted and surviving fractions independently, so tiny
            // conversion or surviving tails are not lost through subtraction.
            let (mut converted, mut survival) = (0.0, 0.0);
            for j in 1..n {
                let exponent = -(j as f64) * logarithm;
                converted += -det::expm1(exponent);
                survival += det::exp(exponent);
            }
            (converted / (n - 1) as f64, survival / (n - 1) as f64)
        }
        RecoveryRateDistribution::Gamma => {
            return Err(RecoveryError::Invalid(
                "gamma closure requires the original shape",
            ));
        }
    };
    if !result.0.is_finite()
        || !result.1.is_finite()
        || result.0 <= 0.0
        || result.0 > 1.0
        || result.1 <= 0.0
        || result.1 > 1.0
    {
        Err(RecoveryError::Arithmetic(
            "positive continuous conversion or survival is unrepresentable",
        ))
    } else {
        Ok(result)
    }
}

#[cfg(test)]
#[path = "continuous_tests.rs"]
mod tests;
