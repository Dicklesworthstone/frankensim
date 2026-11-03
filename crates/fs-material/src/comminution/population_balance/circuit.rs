//! Continuous population-balance grinding closed around a classifier and purge.
//!
//! With fresh feed f, mill discharge p, returned underflow fraction s, classifier
//! underflow r_i, selection S_i and residence tau, each class satisfies
//! `(1 - s*r_i + tau*S_i)*p_i = f_i + tau*sum_j b_ij*S_j*p_j`.
//! Strictly coarser j contribute, so dry mass and every original elemental
//! account have one deterministic triangular solve. Carrier water obeys its own
//! classifier/recycle balance. The continuous perfect-mixer law remains the
//! original [`ContinuousPopulationBalanceMill`], not a batch survival matrix.
//!
//! The positive escape `purged + returned*overflow_i` retains independently
//! supplied small partition tails instead of subtracting a rounded return from
//! one. A fed trapped class or water loop refuses; unforced trapped circulation
//! is set to zero, matching the existing fixed-pass circuit convention.
//!
//! The original mill evaluates the reconstructed circulating feed, including
//! its actual residence, dry throughput and solids fraction. The original
//! classifier evaluates that mill's discharge. Independently recomputed local
//! and external balances must satisfy the caller's original circuit tolerances.
//! No rates, daughter fractions, assays or partitions are normalized to close
//! a balance. Source rounding and binary64 residuals are not certificates.
//!
//! Residence is explicitly controlled. The returned hold-up is its consequence;
//! no residence, power draw, equipment capacity, liberation, or plant calibration
//! is inferred from a fresh-feed rate, Bond index or equipment name.

use super::{
    ContinuousMillSolution, ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS,
    MILL_DRY_FLOW_AXIS, PopulationBalanceError, Sum, product, quotient, require,
};
use crate::comminution::circuit::{CircuitError, CircuitSettings, CircuitSolution, check_closure};
use crate::separation::{
    ClassifierCard, ParticleClass, ParticleStream, Partition, SOLIDS_FRACTION_AXIS, SeparationError,
};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, MassFlowRate, Time};

/// Original model or stream refusal, with no partially accepted circuit.
#[derive(Clone, Debug, PartialEq)]
pub enum ContinuousCircuitError {
    /// The continuous mill refused its inputs, source support or arithmetic.
    Population(PopulationBalanceError),
    /// The original classifier/circuit owner refused flow, balance or settings.
    Circuit(CircuitError),
}
impl core::fmt::Display for ContinuousCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous grinding circuit: {self:?}")
    }
}
impl std::error::Error for ContinuousCircuitError {}
impl From<PopulationBalanceError> for ContinuousCircuitError {
    fn from(error: PopulationBalanceError) -> Self {
        Self::Population(error)
    }
}
impl From<CircuitError> for ContinuousCircuitError {
    fn from(error: CircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<SeparationError> for ContinuousCircuitError {
    fn from(error: SeparationError) -> Self {
        Self::Circuit(CircuitError::Separation(error))
    }
}
type Result<T> = core::result::Result<T, ContinuousCircuitError>;

/// Original continuous mill and classifier, without fitted or inferred settings.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousGrindingCircuit {
    mill: ContinuousPopulationBalanceMill,
    classifier: ClassifierCard,
}

/// Closed steady streams together with original source-admitted mill hold-up.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousCircuitSolution {
    model: ContinuousGrindingCircuit,
    /// Existing circuit accounts: feed, discharge, overflow, recycle and purge.
    pub circuit: CircuitSolution,
    /// Recomputed original mill response, retained source query and class hold-up.
    pub mill: ContinuousMillSolution,
}
impl ContinuousCircuitSolution {
    /// Original selection/daughter data and complete classifier source card.
    pub fn model(&self) -> &ContinuousGrindingCircuit {
        &self.model
    }
}

impl ContinuousGrindingCircuit {
    /// Retain both original laws; their actual solved operating point is admitted
    /// before returning any circuit result.
    pub const fn new(mill: ContinuousPopulationBalanceMill, classifier: ClassifierCard) -> Self {
        Self { mill, classifier }
    }
    /// Original continuous selection and daughter source.
    pub const fn mill(&self) -> &ContinuousPopulationBalanceMill {
        &self.mill
    }
    /// Original classifier, including its water partition and source conditions.
    pub const fn classifier(&self) -> &ClassifierCard {
        &self.classifier
    }

    /// Solve one continuous closed circuit at a prescribed common residence.
    /// Fresh feed includes makeup water. `underflow_recycle.product()` returns
    /// to the mill; its independently retained reject fraction leaves as purge.
    /// Every class, element and water uses that same nonselective split.
    ///
    /// Sources apply at the solved circulating state. Caller conditions cannot
    /// replace actual residence, dry flow or solids fraction. The returned
    /// overflow alone need not equal fresh feed when purge leaves; their combined
    /// external accounts must close. Inputs are not mutated on any refusal.
    pub fn at_residence_time(
        &self,
        fresh: &ParticleStream,
        residence: Time,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        settings: CircuitSettings,
        underflow_recycle: Partition,
    ) -> Result<ContinuousCircuitSolution> {
        settings.validate()?;
        require(
            residence.value().is_finite() && residence.value() > 0.0,
            "continuous circuit residence must be finite and positive",
        )?;
        let dry_fresh = fresh.dry_solids()?.value();
        require(
            dry_fresh > 0.0,
            "continuous circuit requires positive fresh dry flow",
        )?;
        let classes = self.mill.classes();
        require(
            fresh.classes().len() == classes.len()
                && fresh
                    .classes()
                    .iter()
                    .zip(classes)
                    .all(|(a, b)| a.id() == b.id && a.diameter() == b.diameter),
            "fresh feed must match the original mill's complete class layout",
        )?;
        for axis in [
            MEAN_RESIDENCE_TIME_AXIS,
            MILL_DRY_FLOW_AXIS,
            SOLIDS_FRACTION_AXIS,
        ] {
            require(
                !mill_conditions.axes().contains_key(axis),
                "caller overrides an owned continuous-mill coordinate",
            )?;
        }
        let returned = underflow_recycle.product().value();
        let purged = underflow_recycle.reject().value();
        let partitions = classes
            .iter()
            .map(|class| self.classifier.curve_partition(class.diameter))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let water_exit = escape(
            returned,
            purged,
            self.classifier.water_partition().reject().value(),
        )?;
        let water = match (fresh.water().value(), water_exit) {
            (0.0, _) => 0.0,
            (_, 0.0) => return Err(CircuitError::TrappedWater.into()),
            (amount, exit) => quotient(amount, exit)?,
        };
        let elements: Vec<_> = fresh.classes()[0].elements().keys().cloned().collect();
        let channels = elements.len() + 1;
        let mut incoming = vec![vec![Sum::default(); channels]; classes.len()];
        let mut discharge = vec![vec![0.0; channels]; classes.len()];
        for (i, class) in fresh.classes().iter().enumerate() {
            incoming[i][0].add(class.dry_solids().value())?;
            for (j, element) in elements.iter().enumerate() {
                incoming[i][j + 1].add(class.elements()[element].value())?;
            }
        }
        let mut order: Vec<_> = (0..classes.len()).collect();
        order.sort_by(|&a, &b| {
            classes[b]
                .diameter
                .value()
                .total_cmp(&classes[a].diameter.value())
                .then_with(|| classes[a].id.cmp(&classes[b].id))
        });
        for i in order {
            let column = &self.mill.kinetics()[i];
            let selected = product(residence.value(), column.selection_rate.value())?;
            let escape = escape(returned, purged, partitions[i].reject().value())?;
            let denominator = escape + selected;
            if !denominator.is_finite() {
                return Err(
                    CircuitError::Arithmetic("continuous class denominator overflows").into(),
                );
            }
            for channel in 0..channels {
                let numerator = incoming[i][channel].value()?;
                if denominator == 0.0 && numerator > 0.0 {
                    return Err(CircuitError::TrappedClass {
                        class: classes[i].id.clone(),
                    }
                    .into());
                }
                let out = if numerator == 0.0 {
                    0.0
                } else {
                    quotient(numerator, denominator)?
                };
                discharge[i][channel] = out;
                let broken = product(selected, out)?;
                for daughter in &column.daughters {
                    let destination = classes
                        .binary_search_by(|class| class.id.cmp(&daughter.output))
                        .expect("original mill admitted every daughter identity");
                    incoming[destination][channel]
                        .add(product(broken, daughter.mass_fraction.value())?)?;
                }
            }
        }
        let predicted_classes = classes
            .iter()
            .zip(discharge)
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
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let predicted = ParticleStream::new(predicted_classes, MassFlowRate::new(water))?;
        let classifier_splits = classes
            .iter()
            .zip(&partitions)
            .map(|(class, partition)| (class.id.clone(), *partition))
            .collect::<Vec<_>>();
        let recycle_splits = classes
            .iter()
            .map(|class| (class.id.clone(), underflow_recycle))
            .collect::<Vec<_>>();
        // Candidate arithmetic uses original partitions without minting source
        // authority. Both original physical owners admit the actual state below.
        let candidate = predicted.split(&classifier_splits, self.classifier.water_partition())?;
        let candidate_recycle = candidate
            .product
            .split(&recycle_splits, underflow_recycle)?;
        let mill_feed = ParticleStream::mix(&[fresh, &candidate_recycle.product])?;
        let mill = self
            .mill
            .at_residence_time(&mill_feed, residence, mill_conditions)?;
        let separated = self
            .classifier
            .classify(&mill.outlet, classifier_conditions)?;
        let split = separated
            .product
            .split(&recycle_splits, underflow_recycle)?;
        let circulating_load =
            Dimensionless::new(quotient(split.product.dry_solids()?.value(), dry_fresh)?);
        if circulating_load.value() > settings.maximum_circulating_load.value() {
            return Err(CircuitError::CirculatingLoadExceeded {
                actual: circulating_load,
                maximum: settings.maximum_circulating_load,
            }
            .into());
        }
        let reconstructed_feed = ParticleStream::mix(&[fresh, &split.product])?;
        let reconstructed_underflow = ParticleStream::mix(&[&split.product, &split.reject])?;
        let external = ParticleStream::mix(&[&separated.reject, &split.reject])?;
        let closure = check_closure(
            fresh,
            &[
                (&predicted, &mill.outlet),
                (&mill_feed, &reconstructed_feed),
                (&separated.product, &reconstructed_underflow),
            ],
            &separated.reject,
            &external,
            &elements,
            settings,
        )?;
        let circuit = CircuitSolution {
            fresh_feed: fresh.clone(),
            mill_feed,
            mill_discharge: mill.outlet.clone(),
            classifier_underflow: separated.product,
            underflow_recycle,
            recycle: split.product,
            purge: split.reject,
            product: separated.reject,
            circulating_load,
            closure,
        };
        Ok(ContinuousCircuitSolution {
            model: self.clone(),
            circuit,
            mill,
        })
    }
}

fn escape(returned: f64, purged: f64, overflow: f64) -> Result<f64> {
    let coefficient = returned.mul_add(overflow, purged);
    if !coefficient.is_finite()
        || (coefficient == 0.0 && (purged > 0.0 || (returned > 0.0 && overflow > 0.0)))
    {
        return Err(CircuitError::Arithmetic("positive circuit escape is unrepresentable").into());
    }
    Ok(coefficient)
}

#[cfg(test)]
#[path = "circuit/tests.rs"]
mod tests;
