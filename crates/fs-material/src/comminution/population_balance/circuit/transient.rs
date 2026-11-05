//! Finite-inventory response of the original closed grinding/classifier circuit.
//!
//! With instantaneous return and no classifier or return-line inventory,
//! `dH_i/dt = fresh_i - (S_i + escape_i/tau) H_i + sum_j b_ij S_j H_j`.
//! `escape_i = purge + returned*overflow_i` preserves the original small tails.
//! This is accumulation added to the parent circuit's admitted mass balance;
//! it uses the same local perfect-mixer law, fixed selection/daughter data and
//! original classifier. Gross mill discharge is H_i/tau, not external product.
//!
//! The shared positive backward-Euler step belongs to population_balance::transient.
//! Every accepted substep re-evaluates the original classifier on its actual
//! discharge and the original mill's source support on fresh plus actual return.
//! Inventory and integrated mill feed/discharge, return, overflow and purge are
//! separately retained. Internal return is excluded from whole-circuit output.
//! No balance remainder repairs a stream. Conservation is independently checked
//! both around the mill and around the full circuit.
//!
//! Residence, fresh feed, source parameters and return fraction are fixed during
//! each interval. Feed changes use separate calls carrying inventory forward.
//! The source checks cover the initial point and each numerical endpoint, not
//! a continuous-time enclosure. This first-order trajectory supplies no error
//! bound, dynamic sump/pipe delay, shutdown model, power/rate calibration,
//! liberation or transfer of laboratory batch parameters to a plant.

use super::{ContinuousCircuitError, ContinuousGrindingCircuit, escape};
use crate::comminution::circuit::CircuitError;
use crate::comminution::population_balance::{
    ContinuousPopulationBalanceMill, PopulationBalanceError, Sum, product, quotient, require,
    transient::{
        MillInventory, TransientMillBalance, TransientMillError, TransientMillSettings,
        admit_settings, close_accounts, implicit_step, inventory_from_channels, poll,
    },
};
use crate::separation::{ParticleStream, Partition, SeparationError};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Time, chemistry::ElementId};

/// Original source/stream refusal, exhausted work, cancellation or failed balance.
#[derive(Clone, Debug, PartialEq)]
pub enum TransientCircuitError {
    /// The original finite-inventory numerical/material owner refused.
    Mill(TransientMillError),
    /// The original classifier or circulating-load constraint refused.
    Circuit(CircuitError),
}
impl core::fmt::Display for TransientCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient closed grinding: {self:?}")
    }
}
impl std::error::Error for TransientCircuitError {}
impl From<TransientMillError> for TransientCircuitError {
    fn from(value: TransientMillError) -> Self {
        Self::Mill(value)
    }
}
impl From<PopulationBalanceError> for TransientCircuitError {
    fn from(value: PopulationBalanceError) -> Self {
        Self::Mill(value.into())
    }
}
impl From<CircuitError> for TransientCircuitError {
    fn from(value: CircuitError) -> Self {
        Self::Circuit(value)
    }
}
impl From<SeparationError> for TransientCircuitError {
    fn from(value: SeparationError) -> Self {
        Self::Circuit(value.into())
    }
}
impl From<ContinuousCircuitError> for TransientCircuitError {
    fn from(value: ContinuousCircuitError) -> Self {
        match value {
            ContinuousCircuitError::Population(e) => e.into(),
            ContinuousCircuitError::Circuit(e) => e.into(),
        }
    }
}
type Result<T> = core::result::Result<T, TransientCircuitError>;

/// Explicit time resolution, work/conservation budget and dry circulation ceiling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientCircuitSettings {
    /// Existing finite-inventory step, work and extensive kg balance tolerances.
    pub time: TransientMillSettings,
    /// Maximum returned / fresh dry kg/s, admitted at every numerical point.
    pub maximum_circulating_load: Dimensionless,
}

/// Actual instantaneous streams, not a steady-state circuit solution.
#[derive(Clone, Debug, PartialEq)]
pub struct TransientCircuitPoint {
    /// Fresh feed plus the original classifier's returned underflow, kg/s.
    pub mill_feed: ParticleStream,
    /// Gross size-resolved discharge from retained inventory divided by residence.
    pub mill_discharge: ParticleStream,
    /// Original classifier underflow before return/purge splitting.
    pub classifier_underflow: ParticleStream,
    /// Returned underflow, internal to the circuit boundary.
    pub recycle: ParticleStream,
    /// Nonselective underflow bleed, external to the circuit boundary.
    pub purge: ParticleStream,
    /// Original classifier overflow, external to the circuit boundary.
    pub product: ParticleStream,
    /// Actual returned / fresh dry flow at this point.
    pub circulating_load: Dimensionless,
    /// Original mill query with actual circulating feed, residence and solids fraction.
    pub mill_conditions: QueryPoint,
}

/// Admitted finite interval with the original model, inputs and both boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct TransientCircuitStep {
    model: ContinuousGrindingCircuit,
    initial: MillInventory,
    fresh: ParticleStream,
    /// Original classifier query; size and solids fraction are computed at each point.
    pub classifier_conditions: QueryPoint,
    /// Unmodified return/purge split shared by every class and carrier water.
    pub underflow_recycle: Partition,
    /// Fixed common mill residence, distinct from interval elapsed time.
    pub residence: Time,
    /// Actual finite interval advanced.
    pub elapsed: Time,
    /// Fully admitted initial streams with this interval's fresh feed.
    pub initial_point: TransientCircuitPoint,
    /// Retained final size, dry, elemental and water inventory.
    pub inventory: MillInventory,
    /// Actual final instantaneous streams; output need not equal fresh input.
    pub final_point: TransientCircuitPoint,
    /// Original fresh feed integrated over the interval, kg.
    pub integrated_fresh: MillInventory,
    /// Actual fresh plus returned flow integrated over accepted steps, kg.
    pub integrated_mill_feed: MillInventory,
    /// Gross mill discharge integrated over accepted steps, kg.
    pub integrated_mill_discharge: MillInventory,
    /// Classifier underflow integrated before return/purge splitting, kg.
    pub integrated_underflow: MillInventory,
    /// Internal return integrated over accepted steps, kg.
    pub integrated_recycle: MillInventory,
    /// External overflow product integrated over accepted steps, kg.
    pub integrated_product: MillInventory,
    /// External purge integrated over accepted steps, kg.
    pub integrated_purge: MillInventory,
    /// Initial + actual mill feed - final - gross mill discharge.
    pub mill_balance: TransientMillBalance,
    /// Initial + fresh - final - overflow - purge. Recycle is internal.
    pub boundary_balance: TransientMillBalance,
    /// Largest admitted dry circulating-load ratio among the initial/step points.
    pub peak_circulating_load: Dimensionless,
    /// Original time, work, conservation and physical circulation limits.
    pub settings: TransientCircuitSettings,
    /// Number of completed implicit steps.
    pub substeps: usize,
}
impl TransientCircuitStep {
    /// Original mill kinetics and classifier source, unchanged by this response.
    pub fn model(&self) -> &ContinuousGrindingCircuit {
        &self.model
    }
    /// Original immutable initial inventory.
    pub fn initial(&self) -> &MillInventory {
        &self.initial
    }
    /// Original constant fresh feed, including makeup water.
    pub fn fresh_feed(&self) -> &ParticleStream {
        &self.fresh
    }
}

impl ContinuousGrindingCircuit {
    /// Advance the closed circuit at fixed local source parameters and residence.
    ///
    /// Fresh dry flow must be positive so the declared circulating-load limit
    /// remains defined. Initial wet inventory must give an admitted classifier
    /// solids fraction. This scope covers retained-inventory feed changes; it
    /// does not invent a start-from-empty classifier point or a shutdown policy.
    /// A fed trapped loop may accumulate finite inventory, unlike a steady solve,
    /// but every numerical point must still satisfy original support and limits.
    /// Failure or cancellation returns no partial result and mutates no input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_inventory(
        &self,
        initial: &MillInventory,
        fresh: &ParticleStream,
        residence: Time,
        elapsed: Time,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        underflow_recycle: Partition,
        settings: TransientCircuitSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransientCircuitStep> {
        poll(&mut cancelled)?;
        let substeps = admit_settings(elapsed, settings.time)?;
        require(
            settings.maximum_circulating_load.value().is_finite()
                && settings.maximum_circulating_load.value() >= 0.0,
            "dry circulating-load ceiling must be finite and nonnegative",
        )?;
        require(
            fresh.dry_solids()?.value() > 0.0,
            "transient closed circuit requires positive fresh dry flow",
        )?;
        require(
            initial.classes().len() == fresh.classes().len()
                && initial.classes().iter().zip(fresh.classes()).all(|(a, b)| {
                    a.id == b.id()
                        && a.diameter == b.diameter()
                        && a.elements.keys().eq(b.elements().keys())
                }),
            "inventory and fresh feed require identical classes and explicit assays",
        )?;
        let recycle_splits: Vec<_> = self
            .mill
            .classes()
            .iter()
            .map(|c| (c.id.clone(), underflow_recycle))
            .collect();
        let point = |inventory: &MillInventory| {
            self.transient_point(
                inventory,
                fresh,
                residence,
                mill_conditions,
                classifier_conditions,
                &recycle_splits,
                underflow_recycle,
                settings,
            )
        };
        // Validate the real initial circulating point before any numerical work.
        let initial_point = point(initial)?;
        let mut final_point = initial_point.clone();
        let mut peak = initial_point.circulating_load.value();
        let returned = underflow_recycle.product().value();
        let purged = underflow_recycle.reject().value();
        let exit = self
            .mill
            .classes()
            .iter()
            .map(|c| {
                escape(
                    returned,
                    purged,
                    self.classifier
                        .curve_partition(c.diameter)?
                        .reject()
                        .value(),
                )
            })
            .collect::<core::result::Result<Vec<_>, ContinuousCircuitError>>()?;
        let water_exit = escape(
            returned,
            purged,
            self.classifier.water_partition().reject().value(),
        )?;
        let elements: Vec<_> = initial.classes()[0].elements.keys().cloned().collect();
        let channels = elements.len() + 1;
        let mut masses: Vec<Vec<f64>> = initial
            .classes()
            .iter()
            .map(|c| {
                core::iter::once(c.dry_solids.value())
                    .chain(elements.iter().map(|e| c.elements[e].value()))
                    .collect()
            })
            .collect();
        let rates: Vec<Vec<f64>> = fresh
            .classes()
            .iter()
            .map(|c| {
                core::iter::once(c.dry_solids().value())
                    .chain(elements.iter().map(|e| c.elements()[e].value()))
                    .collect()
            })
            .collect();
        let mut water = initial.water().value();
        let mut order: Vec<_> = (0..self.mill.classes().len()).collect();
        order.sort_by(|&a, &b| {
            self.mill.classes()[b]
                .diameter
                .value()
                .total_cmp(&self.mill.classes()[a].diameter.value())
                .then_with(|| self.mill.classes()[a].id.cmp(&self.mill.classes()[b].id))
        });
        let mut accounts: Vec<_> = (0..7)
            .map(|_| Accumulator::new(masses.len(), channels))
            .collect();
        let integrated_fresh = MillInventory::integrate(fresh, elapsed)?;
        let nominal_step = quotient(elapsed.value(), substeps as f64)?;
        let mut inventory = initial.clone();
        for step in 0..substeps {
            poll(&mut cancelled)?;
            let dt = if step + 1 == substeps {
                elapsed.value() - nominal_step * step as f64
            } else {
                nominal_step
            };
            require(
                dt.is_finite() && dt > 0.0 && dt <= settings.time.maximum_step.value(),
                "requested time partition is not representable within the step limit",
            )?;
            implicit_step(
                &self.mill,
                &mut masses,
                &rates,
                &mut water,
                fresh.water().value(),
                dt,
                residence,
                &exit,
                water_exit,
                &order,
                &mut cancelled,
            )?;
            inventory = inventory_from_channels(&self.mill, &elements, masses.clone(), water)?;
            final_point = point(&inventory)?;
            peak = peak.max(final_point.circulating_load.value());
            let external = ParticleStream::mix(&[&final_point.product, &final_point.purge])?;
            for (account, stream) in accounts.iter_mut().zip([
                &final_point.mill_feed,
                &final_point.mill_discharge,
                &final_point.classifier_underflow,
                &final_point.recycle,
                &final_point.product,
                &final_point.purge,
                &external,
            ]) {
                account.add(stream, &elements, dt)?;
            }
        }
        let mut totals = accounts
            .into_iter()
            .map(|a| a.finish(&self.mill, &elements));
        let mut next = || totals.next().expect("seven declared integral accounts");
        let integrated_mill_feed = next()?;
        let integrated_mill_discharge = next()?;
        let integrated_underflow = next()?;
        let integrated_recycle = next()?;
        let integrated_product = next()?;
        let integrated_purge = next()?;
        let integrated_external = next()?;
        let mill_balance = close_accounts(
            initial,
            &integrated_mill_feed,
            &inventory,
            &integrated_mill_discharge,
            settings.time,
        )?;
        let boundary_balance = close_accounts(
            initial,
            &integrated_fresh,
            &inventory,
            &integrated_external,
            settings.time,
        )?;
        poll(&mut cancelled)?;
        Ok(TransientCircuitStep {
            model: self.clone(),
            initial: initial.clone(),
            fresh: fresh.clone(),
            classifier_conditions: classifier_conditions.clone(),
            underflow_recycle,
            residence,
            elapsed,
            initial_point,
            inventory,
            final_point,
            integrated_fresh,
            integrated_mill_feed,
            integrated_mill_discharge,
            integrated_underflow,
            integrated_recycle,
            integrated_product,
            integrated_purge,
            mill_balance,
            boundary_balance,
            peak_circulating_load: Dimensionless::new(peak),
            settings,
            substeps,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn transient_point(
        &self,
        inventory: &MillInventory,
        fresh: &ParticleStream,
        residence: Time,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        recycle_splits: &[(String, Partition)],
        underflow_recycle: Partition,
        settings: TransientCircuitSettings,
    ) -> Result<TransientCircuitPoint> {
        let mill_discharge = inventory.discharge_at(residence)?;
        let classified = self
            .classifier
            .classify(&mill_discharge, classifier_conditions)?;
        let split = classified
            .product
            .split(recycle_splits, underflow_recycle)?;
        let circulating_load = Dimensionless::new(quotient(
            split.product.dry_solids()?.value(),
            fresh.dry_solids()?.value(),
        )?);
        if circulating_load.value() > settings.maximum_circulating_load.value() {
            return Err(CircuitError::CirculatingLoadExceeded {
                actual: circulating_load,
                maximum: settings.maximum_circulating_load,
            }
            .into());
        }
        let mill_feed = ParticleStream::mix(&[fresh, &split.product])?;
        let actual = self
            .mill
            .admitted_feed(&mill_feed, residence, mill_conditions, false)?;
        Ok(TransientCircuitPoint {
            mill_feed,
            mill_discharge,
            classifier_underflow: classified.product,
            recycle: split.product,
            purge: split.reject,
            product: classified.reject,
            circulating_load,
            mill_conditions: actual,
        })
    }
}

struct Accumulator {
    masses: Vec<Vec<Sum>>,
    water: Sum,
}
impl Accumulator {
    fn new(classes: usize, channels: usize) -> Self {
        Self {
            masses: vec![vec![Sum::default(); channels]; classes],
            water: Sum::default(),
        }
    }
    fn add(&mut self, stream: &ParticleStream, elements: &[ElementId], dt: f64) -> Result<()> {
        for (row, class) in self.masses.iter_mut().zip(stream.classes()) {
            row[0].add(product(dt, class.dry_solids().value())?)?;
            for (cell, element) in row[1..].iter_mut().zip(elements) {
                cell.add(product(dt, class.elements()[element].value())?)?;
            }
        }
        self.water.add(product(dt, stream.water().value())?)?;
        Ok(())
    }
    fn finish(
        self,
        model: &ContinuousPopulationBalanceMill,
        elements: &[ElementId],
    ) -> Result<MillInventory> {
        Ok(inventory_from_channels(
            model,
            elements,
            self.masses
                .into_iter()
                .map(|row| row.into_iter().map(Sum::value).collect())
                .collect::<core::result::Result<_, PopulationBalanceError>>()?,
            self.water.value()?,
        )?)
    }
}

#[cfg(test)]
#[path = "transient/tests.rs"]
mod tests;
