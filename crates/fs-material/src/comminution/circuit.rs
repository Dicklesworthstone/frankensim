//! Steady closed grinding/classification with original, fixed process cards.
//!
//! For fresh feed f, mill feed x, per-pass redistribution B and classifier
//! underflow fractions r, solve `x = f + diag(r) B x`. Strictly smaller daughter
//! classes make this triangular in decreasing diameter. Each material account
//! is solved once; no fixed-point iteration, guessed circulating load, inferred
//! liberation or population-balance kinetics is introduced. Water obeys its
//! own balance. An unforced trapped class is assigned zero initial circulation;
//! a fed trapped class has no finite steady solution and refuses.
//!
//! Coefficients must apply to the actual solved circuit. The original mill and
//! classifier evaluate that state before any result is published. Fixed source
//! coefficients are an explicit approximation, not a prediction of grinding
//! environment effects (see Petrakis 2023, doi:10.3390/ma16247687).
//! No power law is inferred. Bond energy, if used by a consumer, is on fresh
//! circuit throughput and circuit F80/P80; multiplying it by mill circulating
//! tonnage would change its basis (GMG01-MP-2021, sections 4–5).

use super::breakage::{BreakageError, SizeRedistributionCard};
use crate::separation::{
    ClassifierCard, ParticleClass, ParticleStream, SOLIDS_FRACTION_AXIS, SeparationError,
};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, MassFlowRate, QuantitySpec, chemistry::ElementId};
use std::collections::BTreeMap;

/// Optional source axis bound by the solver to actual dry mill-feed throughput.
pub const MILL_DRY_FLOW_AXIS: &str = "mill-dry-solids-flow";

/// Caller-declared numerical closure allowance and physical circulating-load limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircuitSettings {
    /// Absolute allowance for every dry-solid, water and element balance, kg/s.
    pub absolute_balance_tolerance: MassFlowRate,
    /// Additional allowance times the larger of the two compared flows.
    pub relative_balance_tolerance: Dimensionless,
    /// Maximum recycle dry solids / fresh dry solids, as a fraction, not percent.
    pub maximum_circulating_load: Dimensionless,
}

/// Refusal returns no partially accepted or clipped circuit result.
#[derive(Clone, Debug, PartialEq)]
pub enum CircuitError {
    /// Original mill-source or per-pass redistribution refusal.
    Breakage(BreakageError),
    /// Original classifier or material-stream refusal.
    Separation(SeparationError),
    /// Invalid settings, source coordinates or unsupported feed layout.
    Invalid(&'static str),
    /// A positively fed size class has neither a product nor a smaller-size exit.
    TrappedClass {
        /// Original size-class identity.
        class: String,
    },
    /// Positive water input has no water outlet.
    TrappedWater,
    /// Positive or finite physical flows cannot be represented in binary64.
    Arithmetic(&'static str),
    /// Solved dry circulating load exceeds the caller's declared limit.
    CirculatingLoadExceeded {
        /// Recycle / fresh dry solids.
        actual: Dimensionless,
        /// Caller-declared maximum on that same basis.
        maximum: Dimensionless,
    },
    /// An independently recomputed local or overall balance exceeds its allowance.
    BalanceExceeded {
        /// Absolute discrepancy in kg/s.
        absolute_error: MassFlowRate,
        /// Discrepancy divided by the larger compared flow.
        relative_error: Dimensionless,
    },
}
impl core::fmt::Display for CircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed grinding circuit: {self:?}")
    }
}
impl std::error::Error for CircuitError {}
impl From<BreakageError> for CircuitError {
    fn from(value: BreakageError) -> Self {
        Self::Breakage(value)
    }
}
impl From<SeparationError> for CircuitError {
    fn from(value: SeparationError) -> Self {
        Self::Separation(value)
    }
}

/// Numerical discrepancies from recomputed physical stream accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct CircuitClosure {
    /// Largest absolute discrepancy across local and overall accounts, kg/s.
    pub maximum_absolute_error: MassFlowRate,
    /// Largest relative discrepancy; zero for two zero compared flows.
    pub maximum_relative_error: Dimensionless,
    /// Product minus fresh-feed dry solids, kg/s.
    pub product_minus_fresh_solids: MassFlowRate,
    /// Product minus fresh-feed water, kg/s.
    pub product_minus_fresh_water: MassFlowRate,
    /// Product minus fresh-feed explicitly assayed elements, kg/s.
    pub product_minus_fresh_elements: BTreeMap<ElementId, MassFlowRate>,
}

/// Closed steady streams; classifier overflow is the circuit's external product.
#[derive(Clone, Debug, PartialEq)]
pub struct CircuitSolution {
    /// Original external input, including all makeup water.
    pub fresh_feed: ParticleStream,
    /// Fresh feed plus classifier underflow, before the mill pass.
    pub mill_feed: ParticleStream,
    /// Conservatively redistributed mill discharge, before classification.
    pub mill_discharge: ParticleStream,
    /// Classifier underflow returned to the mill, including its water.
    pub recycle: ParticleStream,
    /// Classifier overflow leaving the circuit.
    pub product: ParticleStream,
    /// Recycle dry solids divided by fresh dry solids; 2.5 means 250%.
    pub circulating_load: Dimensionless,
    /// Independently checked closure, with no balancing plug or repaired assays.
    pub closure: CircuitClosure,
}

/// Explicit fixed mill-pass and classifier cards; neither is fitted by this owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedGrindingCircuit {
    mill: SizeRedistributionCard,
    classifier: ClassifierCard,
}

impl ClosedGrindingCircuit {
    /// Retain both original cards. Actual source admission happens at solved state.
    pub const fn new(mill: SizeRedistributionCard, classifier: ClassifierCard) -> Self {
        Self { mill, classifier }
    }
    /// Original per-pass redistribution, source context and validity.
    pub const fn mill(&self) -> &SizeRedistributionCard {
        &self.mill
    }
    /// Original classifier parameters, evidence and validity.
    pub const fn classifier(&self) -> &ClassifierCard {
        &self.classifier
    }

    /// Solve a positive fresh-solids feed with all makeup water already included.
    ///
    /// The mill may declare `solids-mass-fraction` and/or `mill-dry-solids-flow`;
    /// these axes are supplied from actual mill feed and cannot be overridden.
    /// All other conditions remain explicit caller inputs. The classifier sees
    /// actual discharge composition. Work is bounded by 128 size classes and
    /// quadratic in that declared count, linear in the explicit element count.
    pub fn solve(
        &self,
        fresh: &ParticleStream,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        settings: CircuitSettings,
    ) -> Result<CircuitSolution, CircuitError> {
        settings.validate()?;
        let dry_fresh = fresh.dry_solids()?.value();
        require(
            dry_fresh > 0.0,
            "circuit requires positive fresh dry-solids flow",
        )?;
        let classes = self.mill.classes();
        require(
            fresh.classes().len() == classes.len()
                && fresh.classes().iter().zip(classes).all(|(input, class)| {
                    input.id() == class.id && input.diameter() == class.diameter
                }),
            "fresh feed must match the mill card's complete size-class layout",
        )?;
        for axis in [SOLIDS_FRACTION_AXIS, MILL_DRY_FLOW_AXIS] {
            require(
                !mill_conditions.axes().contains_key(axis),
                "caller overrides an owned mill coordinate",
            )?;
        }
        let elements: Vec<_> = fresh.classes()[0].elements().keys().cloned().collect();
        let channels = elements.len() + 1;
        let n = classes.len();
        let mut order: Vec<_> = (0..n).collect();
        order.sort_by(|&a, &b| {
            classes[b]
                .diameter
                .value()
                .total_cmp(&classes[a].diameter.value())
                .then(a.cmp(&b))
        });
        let partitions = classes
            .iter()
            .map(|class| self.classifier.curve_partition(class.diameter))
            .collect::<Result<Vec<_>, _>>()?;
        let mut flows = vec![vec![0.0; channels]; n];
        for (position, &i) in order.iter().enumerate() {
            let r = partitions[i].product().value();
            let q = partitions[i].reject().value();
            // Positive escape channels avoid cancellation in 1-r*Bii near a
            // closed loop. The card's independently retained off-diagonal
            // fractions define its breakage exit; columns are never normalized.
            let fragmented = self.mill.fragmented_fraction(&classes[i].id)?.value();
            // This is a coefficient sum, not a physical flow. Fusing avoids
            // rejecting a subnormal r*fragmented when q dominates the result;
            // the actual transported masses are still checked independently.
            let escape = r.mul_add(fragmented, q);
            if !escape.is_finite() || (escape == 0.0 && (q > 0.0 || (r > 0.0 && fragmented > 0.0)))
            {
                return Err(CircuitError::Arithmetic(
                    "escape coefficient is unrepresentable",
                ));
            }
            let mut incoming = vec![Sum::default(); channels];
            for &j in &order[..position] {
                let b = self.mill.fraction(&classes[i].id, &classes[j].id)?.value();
                for (account, value) in incoming.iter_mut().zip(&flows[j]) {
                    account.add(product(b, *value)?)?;
                }
            }
            for (channel, account) in incoming.into_iter().enumerate() {
                let original = if channel == 0 {
                    fresh.classes()[i].dry_solids().value()
                } else {
                    fresh.classes()[i].elements()[&elements[channel - 1]].value()
                };
                let numerator = sum([original, product(r, account.value()?)?])?;
                if escape == 0.0 && numerator > 0.0 {
                    return Err(CircuitError::TrappedClass {
                        class: classes[i].id.clone(),
                    });
                }
                flows[i][channel] = if numerator == 0.0 {
                    0.0
                } else {
                    quotient(numerator, escape)?
                };
            }
        }
        let water_exit = self.classifier.water_partition().reject().value();
        let water = match (fresh.water().value(), water_exit) {
            (0.0, _) => 0.0,
            (_, 0.0) => return Err(CircuitError::TrappedWater),
            (value, exit) => quotient(value, exit)?,
        };
        let feed_classes = classes
            .iter()
            .zip(&flows)
            .map(|(class, row)| {
                ParticleClass::new(
                    class.id.clone(),
                    class.diameter,
                    MassFlowRate::new(row[0]),
                    elements
                        .iter()
                        .cloned()
                        .zip(row[1..].iter().copied().map(MassFlowRate::new))
                        .collect(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mill_feed = ParticleStream::new(feed_classes, MassFlowRate::new(water))?;
        let point = self.actual_mill_conditions(&mill_feed, mill_conditions)?;
        let mill_discharge = self.mill.transform(&mill_feed, &point)?;
        let separated = self
            .classifier
            .classify(&mill_discharge, classifier_conditions)?;
        let circulating_load = Dimensionless::new(quotient(
            separated.product.dry_solids()?.value(),
            dry_fresh,
        )?);
        if circulating_load.value() > settings.maximum_circulating_load.value() {
            return Err(CircuitError::CirculatingLoadExceeded {
                actual: circulating_load,
                maximum: settings.maximum_circulating_load,
            });
        }
        let reconstructed = ParticleStream::mix(&[fresh, &separated.product])?;
        let closure = check_closure(
            fresh,
            &mill_feed,
            &reconstructed,
            &separated.reject,
            &elements,
            settings,
        )?;
        Ok(CircuitSolution {
            fresh_feed: fresh.clone(),
            mill_feed,
            mill_discharge,
            recycle: separated.product,
            product: separated.reject,
            circulating_load,
            closure,
        })
    }

    fn actual_mill_conditions(
        &self,
        feed: &ParticleStream,
        conditions: &QueryPoint,
    ) -> Result<QueryPoint, CircuitError> {
        let mut point = conditions.clone();
        for (axis, quantity, value) in [
            (
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dimensionless::DIMS),
                feed.solids_mass_fraction()?.value(),
            ),
            (
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                feed.dry_solids()?.value(),
            ),
        ] {
            if self.mill.validity().bound(axis).is_some() {
                point = if let Some(declared) = self.mill.validity().axis_quantities().get(axis) {
                    require(
                        *declared == quantity,
                        "owned mill axis has incompatible quantity convention",
                    )?;
                    point.with_quantity(axis, quantity, value)
                } else {
                    point.with(axis, value)
                }
                .map_err(|_| CircuitError::Invalid("invalid actual mill coordinate"))?;
            }
        }
        Ok(point)
    }
}

impl CircuitSettings {
    fn validate(self) -> Result<(), CircuitError> {
        require(
            self.absolute_balance_tolerance.value().is_finite()
                && self.absolute_balance_tolerance.value() >= 0.0,
            "absolute balance tolerance must be finite and nonnegative",
        )?;
        require(
            self.relative_balance_tolerance.value().is_finite()
                && (0.0..1.0).contains(&self.relative_balance_tolerance.value()),
            "relative balance tolerance must be finite and in [0,1)",
        )?;
        require(
            self.maximum_circulating_load.value().is_finite()
                && self.maximum_circulating_load.value() >= 0.0,
            "circulating-load limit must be finite and nonnegative",
        )
    }
}

fn check_closure(
    fresh: &ParticleStream,
    feed: &ParticleStream,
    reconstructed: &ParticleStream,
    outlet: &ParticleStream,
    elements: &[ElementId],
    settings: CircuitSettings,
) -> Result<CircuitClosure, CircuitError> {
    let mut maximum_absolute = 0.0_f64;
    let mut maximum_relative = 0.0_f64;
    let mut check = |left: f64, right: f64| -> Result<f64, CircuitError> {
        let residual = left - right;
        let absolute = residual.abs();
        let scale = left.max(right);
        let relative = if scale > 0.0 { absolute / scale } else { 0.0 };
        if absolute > settings.absolute_balance_tolerance.value()
            && relative
                > settings.relative_balance_tolerance.value()
                    + settings.absolute_balance_tolerance.value() / scale
        {
            return Err(CircuitError::BalanceExceeded {
                absolute_error: MassFlowRate::new(absolute),
                relative_error: Dimensionless::new(relative),
            });
        }
        maximum_absolute = maximum_absolute.max(absolute);
        maximum_relative = maximum_relative.max(relative);
        Ok(residual)
    };
    check(feed.water().value(), reconstructed.water().value())?;
    for (left, right) in feed.classes().iter().zip(reconstructed.classes()) {
        check(left.dry_solids().value(), right.dry_solids().value())?;
        for element in elements {
            check(
                left.elements()[element].value(),
                right.elements()[element].value(),
            )?;
        }
    }
    let dry = check(outlet.dry_solids()?.value(), fresh.dry_solids()?.value())?;
    let water = check(outlet.water().value(), fresh.water().value())?;
    let element_residuals = elements
        .iter()
        .map(|element| {
            Ok((
                element.clone(),
                MassFlowRate::new(check(
                    outlet.element_flow(element)?.value(),
                    fresh.element_flow(element)?.value(),
                )?),
            ))
        })
        .collect::<Result<BTreeMap<_, _>, CircuitError>>()?;
    Ok(CircuitClosure {
        maximum_absolute_error: MassFlowRate::new(maximum_absolute),
        maximum_relative_error: Dimensionless::new(maximum_relative),
        product_minus_fresh_solids: MassFlowRate::new(dry),
        product_minus_fresh_water: MassFlowRate::new(water),
        product_minus_fresh_elements: element_residuals,
    })
}

fn require(value: bool, reason: &'static str) -> Result<(), CircuitError> {
    if value {
        Ok(())
    } else {
        Err(CircuitError::Invalid(reason))
    }
}
fn product(left: f64, right: f64) -> Result<f64, CircuitError> {
    let value = left * right;
    if !value.is_finite() || (left > 0.0 && right > 0.0 && value == 0.0) {
        Err(CircuitError::Arithmetic(
            "positive transfer is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}
fn quotient(numerator: f64, denominator: f64) -> Result<f64, CircuitError> {
    let value = numerator / denominator;
    if !value.is_finite() || (numerator > 0.0 && value == 0.0) {
        Err(CircuitError::Arithmetic(
            "circulating flow or ratio is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}
#[derive(Clone, Copy, Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), CircuitError> {
        let next = self.total + value;
        self.correction += if self.total >= value {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
        self.value().map(|_| ())
    }
    fn value(self) -> Result<f64, CircuitError> {
        let value = self.total + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(CircuitError::Arithmetic("circulating sum overflow"))
        }
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, CircuitError> {
    let mut result = Sum::default();
    for value in values {
        result.add(value)?;
    }
    result.value()
}

#[cfg(test)]
#[path = "circuit/tests.rs"]
mod tests;
