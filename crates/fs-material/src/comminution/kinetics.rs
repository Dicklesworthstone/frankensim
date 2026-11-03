//! Source-bound, constant-rate size-population evolution during grinding.
//!
//! The discrete mass balance is `dm_i/dt = -S_i m_i + sum_j b_ij S_j m_j`.
//! Here S is a measured/fitted selection rate and b is the mass distribution
//! of daughters after one breakage. See Nava, Llorens and Menéndez-Aguado,
//! *Metals* 10 (2020), 1687, §2, <https://doi.org/10.3390/met10121687>.
//! Active columns have no self daughters; zero-rate classes are sinks.
//! A size class and its assayed elements move congruently, without liberation.
//!
//! With q=max(S), P=I+Q/q is nonnegative. Poisson uniformization evaluates
//! exp(Qt) as sum_k exp(-qt)(qt)^k/k! P^k; see van Dijk et al. (2018),
//! <https://doi.org/10.1016/j.peva.2017.09.008>. The declared approximation
//! caps the number of uniformization events at K: all remaining probability
//! stays at P^K, rather than disappearing or renormalizing source columns.
//! The nominal real-arithmetic L1 error is at most twice the Poisson tail.
//! Floating-point/source-column rounding is separate, and no interval
//! certificate or ore calibration is claimed. Final coefficients must pass
//! the original redistribution owner's conservation admission unchanged.
//!
//! Rates and operating conditions remain constant over the contact interval.
//! The result can serve as the existing circuit's per-pass matrix, at the
//! exact prepared conditions. It is not a residence-time distribution,
//! power-draw law, DEM trajectory, changing-rate or mineral-liberation model.

use super::breakage::{
    BreakageColumn, BreakageError, BreakageTestContext, BreakageTransfer, SizeRedistributionCard,
};
use crate::separation::ParticleStream;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{Dimensionless, Frequency, QuantitySpec, Time};

/// Version-1 stateless selection-rate model in an ordinary constitutive card.
pub const SELECTION_BREAKAGE_LAW: &str = "constant-selection-breakage";
/// Owned coordinate for the complete grinding contact interval, in seconds.
pub const CONTACT_TIME_AXIS: &str = "grinding-contact-time";
/// Model parameters are `selection-rate:<class-id>`, in inverse seconds.
pub const RATE_PREFIX: &str = "selection-rate:";
const MAX_POISSON_MEAN: f64 = 64.0;
const MAX_TERMS: usize = 4096;

/// Source admission, unsupported operation or bounded numerical-work failure.
#[derive(Debug, Clone, PartialEq)]
pub enum GrindingKineticsError {
    /// The original material-card/query owner refused the input.
    Source(MatDbError),
    /// The original conservative redistribution owner refused the result/feed.
    Breakage(BreakageError),
    /// Invalid source schema, rate, time or requested numerical budget.
    Invalid(&'static str),
    /// One original source does not cover the complete interval and conditions.
    OutsideValidity,
    /// The declared numerical work is insufficient; no partial pass is returned.
    Budget(&'static str),
    /// Caller cancellation was observed; original cards/streams are unchanged.
    Cancelled,
    /// A required positive number or finite coefficient cannot be represented.
    Arithmetic(&'static str),
}
impl core::fmt::Display for GrindingKineticsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "selection-breakage kinetics: {self:?}")
    }
}
impl std::error::Error for GrindingKineticsError {}
impl From<MatDbError> for GrindingKineticsError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<BreakageError> for GrindingKineticsError {
    fn from(error: BreakageError) -> Self {
        Self::Breakage(error)
    }
}
type Result<T> = core::result::Result<T, GrindingKineticsError>;

/// Explicit truncation and scalar-work limits for cold pass preparation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformizationBudget {
    /// Requested nominal L1 error from the capped Poisson tail, in (0,1).
    /// This excludes binary64 rounding and source-coefficient uncertainty.
    pub maximum_l1_tail: Dimensionless,
    /// Total Poisson terms, including the positive tail summation; at most 4096.
    pub max_poisson_terms: usize,
    /// Maximum scalar products in all sparse transition applications.
    pub max_coefficient_products: usize,
}

/// Immutable original daughter distribution and source-specific selection rates.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectionBreakageModel {
    daughters: SizeRedistributionCard,
    source: ConstitutiveModelCard,
    rates: Vec<Frequency>,
}

/// A prepared nonnegative pass with its original sources and numerical result.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticGrindingPass {
    model: SelectionBreakageModel,
    transition: SizeRedistributionCard,
    point: QueryPoint,
    elapsed: Time,
    /// Largest retained number of uniformization events; later events freeze here.
    pub capped_events: usize,
    /// Nominal analytical bound on probability assigned to that frozen tail.
    /// Calculated in binary64, not an outward-rounded probability enclosure.
    pub poisson_tail_bound: Dimensionless,
    /// Actual scalar products used to prepare all original input columns.
    pub coefficient_products: usize,
}
impl KineticGrindingPass {
    /// Original selection-rate card, daughter fractions and source conditions.
    pub fn model(&self) -> &SelectionBreakageModel {
        &self.model
    }
    /// Derived pass accepted by the existing conservative redistribution owner.
    /// Its validity is pinned to the exact prepared coordinates, including time.
    pub fn transition(&self) -> &SizeRedistributionCard {
        &self.transition
    }
    /// Complete typed operating point that prepared this pass.
    pub fn point(&self) -> &QueryPoint {
        &self.point
    }
    /// Contact duration used by the population balance.
    pub const fn elapsed(&self) -> Time {
        self.elapsed
    }
    /// Congruently transport all dry solids/elements; carrier water is unchanged.
    /// Input stream layouts and numerical ranges remain the original owner's rules.
    pub fn transform(&self, feed: &ParticleStream) -> Result<ParticleStream> {
        Ok(self.transition.transform(feed, &self.point)?)
    }
}

impl SelectionBreakageModel {
    /// Bind one finite nonnegative inverse-second rate to every original class.
    /// Positive-rate classes must have no self daughter; zero-rate classes
    /// remain unchanged regardless of their dormant daughters. The smallest class
    /// cannot lose mass into an unrepresented fines bin.
    pub fn from_card(
        daughters: SizeRedistributionCard,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        if source.law.0 != SELECTION_BREAKAGE_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != daughters.classes().len()
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(GrindingKineticsError::Invalid(
                "requires sourced version 1, zero state and one rate per class",
            ));
        }
        let time_quantity = QuantitySpec::dimensional(Time::DIMS);
        if !matches!(source.validity.bound(CONTACT_TIME_AXIS), Some((lo, hi)) if lo == 0.0 && hi.is_finite() && hi > 0.0)
            || source.validity.axis_quantities().get(CONTACT_TIME_AXIS) != Some(&time_quantity)
            || source.validity.bounds().iter().any(|(axis, (lo, hi))| {
                !lo.is_finite()
                    || !hi.is_finite()
                    || lo > hi
                    || !source.validity.axis_quantities().contains_key(axis)
            })
        {
            return Err(GrindingKineticsError::Invalid(
                "requires typed finite conditions and contact-time support from zero",
            ));
        }
        let mut rates = Vec::with_capacity(daughters.classes().len());
        for class in daughters.classes() {
            let rate = source
                .parameters
                .get(&format!("{RATE_PREFIX}{}", class.id))
                .ok_or(GrindingKineticsError::Invalid(
                    "missing original class selection rate",
                ))?;
            if rate.dims != Frequency::DIMS || !rate.value.is_finite() || rate.value < 0.0 {
                return Err(GrindingKineticsError::Invalid(
                    "selection rates must be finite nonnegative inverse seconds",
                ));
            }
            let self_fraction = daughters.retained_fraction(&class.id)?.value();
            if rate.value > 0.0 && self_fraction != 0.0 {
                return Err(GrindingKineticsError::Invalid(
                    "positive-rate classes require strictly finer daughters",
                ));
            }
            rates.push(Frequency::new(rate.value));
        }
        Ok(Self {
            daughters,
            source,
            rates,
        })
    }

    /// Original source model; selection rates never come from mineral names or Wi.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Original single-breakage daughter distribution and its specimen evidence.
    pub fn daughters(&self) -> &SizeRedistributionCard {
        &self.daughters
    }
    /// Rates in the original daughter's canonical class-identity order.
    pub fn rates(&self) -> &[Frequency] {
        &self.rates
    }

    /// Prepare one fixed-contact-time pass. Poisson means `q*t > 64` refuse
    /// rather than silently selecting equilibrium or a different approximation.
    /// Cancellation is polled per Poisson term and per sparse column application.
    /// Failed preparation never publishes a partial pass or mutates source state.
    pub fn prepare(
        &self,
        elapsed: Time,
        conditions: &QueryPoint,
        budget: UniformizationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticGrindingPass> {
        let t = elapsed.value();
        if !t.is_finite()
            || t < 0.0
            || !budget.maximum_l1_tail.value().is_finite()
            || !(0.0 < budget.maximum_l1_tail.value() && budget.maximum_l1_tail.value() < 1.0)
            || budget.max_poisson_terms < 2
            || budget.max_poisson_terms > MAX_TERMS
        {
            return Err(GrindingKineticsError::Invalid(
                "invalid time or uniformization budget",
            ));
        }
        poll(&mut cancelled)?;
        if let Some(value) = conditions.axes().get(CONTACT_TIME_AXIS) {
            if *value != t
                || conditions.axis_quantities().get(CONTACT_TIME_AXIS)
                    != Some(&QuantitySpec::dimensional(Time::DIMS))
            {
                return Err(GrindingKineticsError::Invalid(
                    "conflicting contact-time coordinate",
                ));
            }
        }
        let point = conditions.clone().with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            t,
        )?;
        let start = conditions.clone().with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
        )?;
        for domain in [&self.source.validity, self.daughters.validity()] {
            if !domain.contains_typed(point.axes(), point.axis_quantities())
                || !domain.contains_typed(start.axes(), start.axis_quantities())
            {
                return Err(GrindingKineticsError::OutsideValidity);
            }
        }
        let q = self.rates.iter().map(|s| s.value()).fold(0.0_f64, f64::max);
        let mut products = 0;
        let (weights, tail, tail_bound) = if t == 0.0 || q == 0.0 {
            (vec![1.0], 0.0, 0.0)
        } else {
            let x = product(q, t)?;
            if x > MAX_POISSON_MEAN {
                return Err(GrindingKineticsError::Budget(
                    "maximum selection-rate times contact time exceeds 64",
                ));
            }
            poisson_weights(x, budget, &mut cancelled)?
        };
        let n = self.rates.len();
        let mut transition_columns = Vec::with_capacity(n);
        // Resolve source identities once; each application is sparse and fixed-order.
        let mut steps = vec![Vec::new(); n];
        for (j, (rate, column)) in self.rates.iter().zip(self.daughters.columns()).enumerate() {
            let selected = if weights.len() == 1 {
                0.0
            } else {
                rate.value() / q
            };
            if weights.len() > 1 && rate.value() > 0.0 && selected == 0.0 {
                return Err(GrindingKineticsError::Arithmetic(
                    "selection-rate ratio underflows",
                ));
            }
            if selected < 1.0 {
                steps[j].push((j, 1.0 - selected));
            }
            if selected > 0.0 {
                for transfer in &column.transfers {
                    if transfer.mass_fraction.value() > 0.0 {
                        let i = self
                            .daughters
                            .classes()
                            .binary_search_by(|c| c.id.cmp(&transfer.output))
                            .map_err(|_| {
                                GrindingKineticsError::Invalid("unknown original daughter")
                            })?;
                        steps[j].push((i, product(selected, transfer.mass_fraction.value())?));
                    }
                }
            }
        }
        for input in 0..n {
            poll(&mut cancelled)?;
            // A zero-rate class is exactly absorbing. Re-evaluating sum p_k
            // would add numerical error to an unchanged physical inventory.
            if self.rates[input].value() == 0.0 || weights.len() == 1 {
                let id = &self.daughters.classes()[input].id;
                transition_columns.push(BreakageColumn {
                    input: id.clone(),
                    transfers: vec![BreakageTransfer {
                        output: id.clone(),
                        mass_fraction: Dimensionless::new(1.0),
                    }],
                });
                continue;
            }
            let mut state = vec![0.0; n];
            state[input] = 1.0;
            let mut sum = vec![Sum::default(); n];
            for (k, &weight) in weights.iter().enumerate() {
                if k != 0 {
                    state = apply(
                        &steps,
                        &state,
                        &mut products,
                        budget.max_coefficient_products,
                        &mut cancelled,
                    )?;
                }
                for (total, &value) in sum.iter_mut().zip(&state) {
                    total.add(product(weight, value)?)?;
                }
            }
            // Tail is evaluated as a positive series, never as 1-rounded_sum.
            // Capping events is explicit; no coefficient is renormalized.
            for (total, &value) in sum.iter_mut().zip(&state) {
                total.add(product(tail, value)?)?;
            }
            let transfers = self
                .daughters
                .classes()
                .iter()
                .zip(sum)
                .map(|(class, value)| {
                    Ok(BreakageTransfer {
                        output: class.id.clone(),
                        mass_fraction: Dimensionless::new(value.value()?),
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            transition_columns.push(BreakageColumn {
                input: self.daughters.classes()[input].id.clone(),
                transfers,
            });
        }
        let mut validity = ValidityDomain::unconstrained();
        for (axis, &value) in point.axes() {
            let quantity =
                point
                    .axis_quantities()
                    .get(axis)
                    .ok_or(GrindingKineticsError::Invalid(
                        "all prepared coordinates require explicit quantities",
                    ))?;
            validity = validity.with_quantity(axis, quantity.clone(), value, value);
        }
        let context = BreakageTestContext {
            specimen: self.daughters.test_context().specimen.clone(),
            method: format!(
                "computed constant-selection population balance at {t} s; daughter method: {}; rate model: {:?}",
                self.daughters.test_context().method,
                self.source.content_hash()
            ),
            provenance: self.source.provenance.clone(),
        };
        let transition = SizeRedistributionCard::new(
            self.daughters.classes().to_vec(),
            transition_columns,
            context,
            validity,
        )?;
        Ok(KineticGrindingPass {
            model: self.clone(),
            transition,
            point,
            elapsed,
            capped_events: weights.len() - 1,
            poisson_tail_bound: Dimensionless::new(tail_bound),
            coefficient_products: products,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(GrindingKineticsError::Cancelled)
    } else {
        Ok(())
    }
}
fn product(a: f64, b: f64) -> Result<f64> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(GrindingKineticsError::Arithmetic(
            "positive coefficient product loses range",
        ))
    } else {
        Ok(value)
    }
}

// A geometric majorant for sum_{j=k+1}^infinity p_j once x/(k+2)<1.
fn tail_bound(x: f64, k: usize, next: f64) -> f64 {
    let ratio = x / (k + 2) as f64;
    if ratio < 1.0 {
        next / (1.0 - ratio)
    } else {
        f64::INFINITY
    }
}
fn poisson_weights(
    x: f64,
    budget: UniformizationBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64, f64)> {
    let mut weights = vec![det::exp(-x)];
    let bound;
    let mut next;
    loop {
        poll(cancelled)?;
        let k = weights.len() - 1;
        next = product(weights[k], x / (k + 1) as f64)?;
        let candidate = tail_bound(x, k, next);
        // At least one positive-event term keeps tiny observable fragmentation.
        if k > 0 && 2.0 * candidate <= budget.maximum_l1_tail.value() {
            bound = candidate;
            break;
        }
        if weights.len() >= budget.max_poisson_terms {
            return Err(GrindingKineticsError::Budget(
                "Poisson truncation tolerance not reached",
            ));
        }
        weights.push(next);
    }
    let mut tail = Sum::default();
    let mut k = weights.len();
    loop {
        poll(cancelled)?;
        if k >= budget.max_poisson_terms {
            return Err(GrindingKineticsError::Budget(
                "positive Poisson-tail summation budget exhausted",
            ));
        }
        tail.add(next)?;
        let following = product(next, x / (k + 1) as f64)?;
        if tail_bound(x, k, following) <= f64::EPSILON * tail.value()? {
            break;
        }
        next = following;
        k += 1;
    }
    Ok((weights, tail.value()?, bound))
}

fn apply(
    steps: &[Vec<(usize, f64)>],
    state: &[f64],
    products: &mut usize,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let mut out = vec![Sum::default(); state.len()];
    for (column, &value) in steps.iter().zip(state) {
        poll(cancelled)?;
        if value == 0.0 {
            continue;
        }
        for &(i, fraction) in column {
            if *products >= maximum {
                return Err(GrindingKineticsError::Budget(
                    "transition product budget exhausted",
                ));
            }
            *products += 1;
            out[i].add(product(value, fraction)?)?;
        }
    }
    out.into_iter().map(Sum::value).collect()
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = self.sum + value;
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.value().map(|_| ())
    }
    fn value(self) -> Result<f64> {
        let value = self.sum + self.correction;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(GrindingKineticsError::Arithmetic(
                "invalid accumulated coefficient",
            ))
        }
    }
}

#[cfg(test)]
#[path = "kinetics/tests.rs"]
mod tests;
