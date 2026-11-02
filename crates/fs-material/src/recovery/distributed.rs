//! Source-bound first-order recovery of one declared mineral/metal population.
//!
//! Batch and ideal plug flow use `R(t) = R_inf * E[1-exp(-k*t)]`.
//! For a uniform rate on `[0,k_max]` this is the rectangular/Klimpel curve;
//! gamma rates use shape `a` and rate-scale `theta` [1/s], with mean `a*theta`.
//! These expressions are tabulated in *Representation of Kinetics Models in
//! Batch Flotation as Distributed First-Order Reactions* (2020), Table 1:
//! <https://doi.org/10.3390/min10100913>.
//!
//! A steady equal-CSTR cascade instead retains a fixed-rate population by
//! `(1+k*t/N)^(-N)`, from each tank's balance `Q(C_in-C_out)=k*V*C_out`.
//! Discrete populations are evaluated separately before averaging. The need
//! to combine kinetics with the actual residence-time distribution is treated
//! in *Flotation rate distribution in the collection zone of industrial cells*
//! (2010): <https://doi.org/10.1016/j.mineng.2010.05.008>.
//!
//! Reactor choice and source conditions are explicit. A batch fit is not a
//! plant scale-up; no rate, ultimate recovery, entrainment, froth loss, reagent
//! response, covariance, or operational correction is inferred. Parameters
//! remain constant over the declared residence and conditions. CSTR responses
//! for continuous rate distributions refuse rather than substituting their
//! means or a hidden quadrature. No transient hold-up or recycle is modeled.
//!
//! [`super::RecoveryKinetics`] retains the original source-card single/uniform
//! contact-time interface. This additive owner supplies gamma and explicit
//! discrete populations, a declared reactor, and independently retained residue.
//! [`super::residence`] owns transient hydraulic hold-up without reactions;
//! these steady kinetic fractions do not advance its retained inventories.

use fs_evidence::ValidityDomain;
use fs_matdb::{Provenance, QueryPoint};
use fs_math::det;
use fs_qty::{Dimensionless, Frequency, QuantitySpec, Time};

/// Required validity axis; evaluation supplies it from the typed time argument.
pub const RESIDENCE_TIME_AXIS: &str = "residence_time_s";
/// Canonical seconds, not minutes or a dimensionless residence coordinate.
pub const RESIDENCE_TIME_QUANTITY: QuantitySpec = QuantitySpec::dimensional(Time::DIMS);
/// Bounded work for explicitly supplied discrete rate populations.
pub const MAX_RECOVERY_POPULATIONS: usize = 256;

/// Immutable source context; numerical examples do not establish calibration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecoveryTestContext {
    /// Sample and recovered component identity, including size class if relevant.
    pub specimen: String,
    /// Preparation, test procedure, and justification of reactor applicability.
    pub method: String,
    /// Original test/fitting source and usage terms.
    pub provenance: Provenance,
}

/// A positive relative share of the recoverable component, not of whole ore.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RatePopulation {
    /// Relative weights are explicitly normalized for averaging; units cancel.
    pub weight: Dimensionless,
    /// Positive first-order rate in reciprocal seconds.
    pub rate: Frequency,
}

/// The source's distribution choice remains inspectable; no fitted defaults.
#[derive(Clone, Debug, PartialEq)]
pub enum RateDistribution {
    /// One fixed positive rate [1/s].
    Single {
        /// Positive rate shared by the entire recoverable population [1/s].
        rate: Frequency,
    },
    /// Uniform rates from zero to a positive maximum [1/s], not the mean.
    Rectangular {
        /// Positive upper endpoint of the uniform rate distribution [1/s].
        maximum_rate: Frequency,
    },
    /// Gamma shape is dimensionless; scale is a rate [1/s], not a time.
    Gamma {
        /// Positive dimensionless shape parameter.
        shape: Dimensionless,
        /// Positive rate-scale parameter [1/s]; mean rate is shape times scale.
        scale: Frequency,
    },
    /// Positive relative recoverable-mass weights and positive rate constants.
    /// Permanently unrecoverable material belongs in `1 - R_inf` instead.
    Discrete {
        /// One through `MAX_RECOVERY_POPULATIONS` explicitly weighted populations.
        populations: Vec<RatePopulation>,
    },
}

/// Residence semantics are part of the card, never inferred from elapsed time.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecoveryReactor {
    /// Closed batch; time since the specified initially unrecovered state.
    Batch,
    /// Ideal plug flow; every feed element has the same residence time.
    PlugFlow,
    /// Steady equal-volume perfect mixers with unchanged flow and rates.
    /// The supplied time is the TOTAL mean residence, not the time per tank.
    EqualCstrCascade {
        /// Positive number of equal-volume, perfectly mixed tanks in series.
        tanks: u32,
    },
}

/// Fractions of the original component amount, not concentrate grade or yield.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecoveryResponse {
    /// Recovered fraction, between zero and the source's ultimate recovery.
    pub recovered: Dimensionless,
    /// Recoverable component still in the residue/tailings.
    pub recoverable_remaining: Dimensionless,
    /// The source's permanently unrecoverable component, `1 - R_inf`.
    pub intrinsically_unrecoverable: Dimensionless,
}

/// Refusal preserves the source and never returns a clamped/extrapolated value.
#[derive(Clone, Debug, PartialEq)]
pub enum RecoveryError {
    /// Invalid parameters, context, residence time, or numerical result.
    Invalid(
        /// The failed physical, source, or numerical admission rule.
        &'static str,
    ),
    /// A declared source condition was not supplied.
    MissingCondition {
        /// Required source-coordinate name absent from the supplied query.
        axis: String,
    },
    /// A condition used a different or missing quantity convention.
    QuantityMismatch {
        /// Coordinate whose quantity descriptor differs from the source.
        axis: String,
    },
    /// A condition left the source's inclusive support.
    OutsideValidity {
        /// Coordinate outside its declared source support.
        axis: String,
        /// Supplied value in the source coordinate's canonical units.
        value: f64,
        /// Inclusive lower source bound in the same units.
        minimum: f64,
        /// Inclusive upper source bound in the same units.
        maximum: f64,
    },
    /// Continuous rate distributions have no implemented mixed-tank closure.
    UnsupportedReactorDistribution,
}

impl core::fmt::Display for RecoveryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "recovery: {self:?}")
    }
}
impl std::error::Error for RecoveryError {}

/// Validated immutable rates, reactor choice, test context, and source support.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveryCard {
    ultimate_recovery: Dimensionless,
    rates: RateDistribution,
    reactor: RecoveryReactor,
    context: RecoveryTestContext,
    validity: ValidityDomain,
}

impl RecoveryCard {
    /// Admit explicitly authored parameters. The validity domain MUST include
    /// `RESIDENCE_TIME_AXIS` with `RESIDENCE_TIME_QUANTITY` and finite bounds.
    /// Additional numeric conditions (grind, pH, dosage, temperature, etc.) use
    /// the shared domain's exact quantity descriptors; no missing axis defaults.
    pub fn new(
        ultimate_recovery: Dimensionless,
        rates: RateDistribution,
        reactor: RecoveryReactor,
        context: RecoveryTestContext,
        validity: ValidityDomain,
    ) -> Result<Self, RecoveryError> {
        let ultimate = ultimate_recovery.value();
        require(
            ultimate.is_finite() && (0.0..=1.0).contains(&ultimate),
            "ultimate recovery must be in [0,1]",
        )?;
        require(
            !context.specimen.trim().is_empty()
                && !context.method.trim().is_empty()
                && !context.provenance.source.trim().is_empty()
                && !context.provenance.license.trim().is_empty(),
            "test context and source/license are required",
        )?;
        for (axis, &(lo, hi)) in validity.bounds() {
            require(
                !axis.trim().is_empty() && lo.is_finite() && hi.is_finite() && lo <= hi,
                "source bounds must be finite and ordered",
            )?;
        }
        require(
            validity
                .bound(RESIDENCE_TIME_AXIS)
                .is_some_and(|(lo, _)| lo >= 0.0)
                && validity.axis_quantities().get(RESIDENCE_TIME_AXIS)
                    == Some(&RESIDENCE_TIME_QUANTITY),
            "a typed nonnegative residence-time domain is required",
        )?;
        match &rates {
            RateDistribution::Single { rate } => positive(rate.value())?,
            RateDistribution::Rectangular { maximum_rate } => positive(maximum_rate.value())?,
            RateDistribution::Gamma { shape, scale } => {
                positive(shape.value())?;
                positive(scale.value())?;
            }
            RateDistribution::Discrete { populations } => {
                require(
                    !populations.is_empty() && populations.len() <= MAX_RECOVERY_POPULATIONS,
                    "discrete population count must be in 1..=256",
                )?;
                let scale = populations
                    .iter()
                    .map(|p| p.weight.value())
                    .fold(0.0, f64::max);
                for population in populations {
                    positive(population.weight.value())?;
                    positive(population.rate.value())?;
                    require(
                        population.weight.value() / scale > 0.0,
                        "positive relative population weight is unrepresentable",
                    )?;
                }
            }
        }
        if let RecoveryReactor::EqualCstrCascade { tanks } = reactor {
            require(tanks > 0, "a CSTR cascade requires at least one tank")?;
            if !matches!(
                rates,
                RateDistribution::Single { .. } | RateDistribution::Discrete { .. }
            ) {
                return Err(RecoveryError::UnsupportedReactorDistribution);
            }
        }
        Ok(Self {
            ultimate_recovery,
            rates,
            reactor,
            context,
            validity,
        })
    }

    /// Original recoverable fraction of the named component.
    pub fn ultimate_recovery(&self) -> Dimensionless {
        self.ultimate_recovery
    }
    /// Original rates and original discrete weights, without rewriting them.
    pub fn rates(&self) -> &RateDistribution {
        &self.rates
    }
    /// Explicit residence/mixing assumption.
    pub fn reactor(&self) -> RecoveryReactor {
        self.reactor
    }
    /// Original testwork and applicability context.
    pub fn test_context(&self) -> &RecoveryTestContext {
        &self.context
    }
    /// Original support, never widened by evaluation.
    pub fn validity(&self) -> &ValidityDomain {
        &self.validity
    }

    /// Evaluate within source support. Omit the reserved residence axis from
    /// `conditions`; it is inserted from `residence_time`, and duplicates refuse.
    /// Zero residence returns zero recovery only when the source admits zero.
    pub fn evaluate(
        &self,
        residence_time: Time,
        conditions: &QueryPoint,
    ) -> Result<RecoveryResponse, RecoveryError> {
        let t = residence_time.value();
        require(
            t.is_finite() && t >= 0.0,
            "residence time must be finite and nonnegative",
        )?;
        require(
            !conditions.axes().contains_key(RESIDENCE_TIME_AXIS),
            "residence time must be supplied only by the typed argument",
        )?;
        let point = conditions
            .clone()
            .with_quantity(RESIDENCE_TIME_AXIS, RESIDENCE_TIME_QUANTITY, t)
            .map_err(|_| RecoveryError::Invalid("invalid typed residence coordinate"))?;
        for (axis, &(minimum, maximum)) in self.validity.bounds() {
            let &value = point
                .axes()
                .get(axis)
                .ok_or_else(|| RecoveryError::MissingCondition { axis: axis.clone() })?;
            if point.axis_quantities().get(axis) != self.validity.axis_quantities().get(axis) {
                return Err(RecoveryError::QuantityMismatch { axis: axis.clone() });
            }
            if value < minimum || value > maximum {
                return Err(RecoveryError::OutsideValidity {
                    axis: axis.clone(),
                    value,
                    minimum,
                    maximum,
                });
            }
        }
        let (recovered_share, remaining_share) = if t == 0.0 {
            (0.0, 1.0)
        } else {
            match &self.rates {
                RateDistribution::Single { rate } => fixed_rate(rate.value(), t, self.reactor),
                RateDistribution::Rectangular { maximum_rate } => {
                    rectangular(maximum_rate.value(), t)
                }
                RateDistribution::Gamma { shape, scale } => gamma(shape.value(), scale.value(), t),
                RateDistribution::Discrete { populations } => {
                    let scale = populations
                        .iter()
                        .map(|p| p.weight.value())
                        .fold(0.0, f64::max);
                    let mut sums = [(0.0, 0.0); 3];
                    for population in populations {
                        let weight = population.weight.value() / scale;
                        let (recovered_share, remaining_share) =
                            fixed_rate(population.rate.value(), t, self.reactor);
                        accumulate(&mut sums[0], weight);
                        accumulate(&mut sums[1], weight * recovered_share);
                        accumulate(&mut sums[2], weight * remaining_share);
                    }
                    let [total, recovered, remaining] =
                        sums.map(|(sum, correction)| sum + correction);
                    (recovered / total, remaining / total)
                }
            }
        };
        require(
            recovered_share.is_finite() && (0.0..=1.0).contains(&recovered_share),
            "recovery response is outside its physical range",
        )?;
        require(
            remaining_share.is_finite() && (0.0..=1.0).contains(&remaining_share),
            "remaining response is outside its physical range",
        )?;
        let ultimate = self.ultimate_recovery.value();
        let recovered = ultimate * recovered_share;
        Ok(RecoveryResponse {
            recovered: Dimensionless::new(recovered),
            recoverable_remaining: Dimensionless::new(ultimate * remaining_share),
            intrinsically_unrecoverable: Dimensionless::new(1.0 - ultimate),
        })
    }
}

fn require(ok: bool, message: &'static str) -> Result<(), RecoveryError> {
    if ok {
        Ok(())
    } else {
        Err(RecoveryError::Invalid(message))
    }
}
fn positive(value: f64) -> Result<(), RecoveryError> {
    require(
        value.is_finite() && value > 0.0,
        "rates, shape, and population weights must be finite and positive",
    )
}

// Neumaier accumulation retains small populations in all three accounts.
fn accumulate(sum: &mut (f64, f64), value: f64) {
    let next = sum.0 + value;
    sum.1 += if sum.0 >= value {
        (sum.0 - next) + value
    } else {
        (value - next) + sum.0
    };
    sum.0 = next;
}

// Deterministic log(1+x), preserving an increment smaller than 1 ULP at one.
fn log_one_plus(x: f64) -> f64 {
    let sum = 1.0 + x;
    if sum == 1.0 {
        x
    } else {
        det::ln(sum) * (x / (sum - 1.0))
    }
}

fn fixed_rate(rate: f64, time: f64, reactor: RecoveryReactor) -> (f64, f64) {
    let x = rate * time;
    let exponent = match reactor {
        RecoveryReactor::Batch | RecoveryReactor::PlugFlow => x,
        RecoveryReactor::EqualCstrCascade { tanks } => {
            let n = f64::from(tanks);
            if x.is_infinite() {
                n * (det::ln(rate) + det::ln(time) - det::ln(n))
            } else {
                let per_tank = x / n;
                if per_tank == 0.0 {
                    x
                } else {
                    x * (log_one_plus(per_tank) / per_tank)
                }
            }
        }
    };
    (-det::expm1(-exponent), det::exp(-exponent))
}

fn rectangular(rate: f64, time: f64) -> (f64, f64) {
    let x = rate * time;
    if x <= 0.125 {
        // 1-(1-exp(-x))/x: finite Taylor sum avoids cancellation at x=0.
        let (mut term, mut result) = (x * 0.5, x * 0.5);
        for n in 2..=12 {
            term *= -x / f64::from(n + 1);
            result += term;
        }
        (result, 1.0 - result)
    } else {
        let remaining = if x.is_infinite() {
            (1.0 / rate) / time
        } else {
            -det::expm1(-x) / x
        };
        (1.0 - remaining, remaining)
    }
}

fn gamma(shape: f64, scale: f64, time: f64) -> (f64, f64) {
    let x = scale * time;
    let exponent = if x < f64::MIN_POSITIVE {
        // Here log(1+x)=x to binary64 precision. Reorder multiplication so
        // a large shape can retain a representable shape*scale*time.
        if shape > 1.0 {
            (shape * scale.min(time)) * scale.max(time)
        } else {
            shape * x
        }
    } else if x.is_infinite() {
        // log(1+scale*time) remains representable even when the product does not.
        shape * (det::ln(scale) + det::ln(time))
    } else {
        shape * log_one_plus(x)
    };
    (-det::expm1(-exponent), det::exp(-exponent))
}

#[cfg(test)]
#[path = "distributed/tests.rs"]
mod tests;
