//! Finite perfectly mixed classifier-feed sump in the original grinding circuit.
//!
//! DESIGN: mill mass M and sump mass Z obey, for each dry/elemental size class,
//! `M' = fresh + return*Z/tau_s - (selection+1/tau_m)*M + daughters` and
//! `Z' = M/tau_m - Z/tau_s`. Water has no breakage and a separately declared
//! makeup-water rate enters Z. Classifier and return pipe have no other stock.
//! Both drain laws and all interval parameters are explicit local assumptions;
//! no tank volume, density, pump curve, slip or plant calibration is inferred.
//!
//! The positive two-volume backward-Euler system has a scalar Schur complement
//! in each size class. With q=tau_s/(tau_s+dt), t=dt/(tau_s+dt), and the original
//! external escape e, the mill's effective escape is q+t*e and its extra source
//! is return*Z_old/(tau_s+dt). This delegates the complete breakage solve to the
//! original coarse-to-fine implicit step. Z_new=q*(Z_old+dt*M_new/tau_m).
//! Makeup water enters both the water Schur source and the independent Z update.
//! There is no explicit delayed-feed approximation or new dense linear solver.
//!
//! Actual classifier feed is Z/tau_s, independently of mill discharge M/tau_m.
//! The original classifier and circulating-feed source checks run at initial
//! and accepted numerical endpoints. Independently integrated flows must close
//! mill, sump and whole-circuit balances. A failure publishes no partial state.
//! First-order trajectories and endpoint support checks are not error enclosures.

use super::{
    Accumulator, ContinuousGrindingCircuit, Result, TransientCircuitPoint,
    TransientCircuitSettings, escape,
};
use crate::comminution::population_balance::{
    ContinuousPopulationBalanceMill, Sum, product, quotient, require,
    transient::{
        MillInventory, TransientMillBalance, admit_settings, close_accounts, implicit_step,
        inventory_from_channels, poll,
    },
};
use crate::separation::{ParticleStream, Partition};
use fs_matdb::{Provenance, QueryPoint};
use fs_qty::{Dimensionless, Mass, MassFlowRate, Time, chemistry::ElementId};

/// An explicitly declared common perfectly mixed sump drain time.
/// Every size class, its elements and water share this time; no two-phase slip.
#[derive(Clone, Debug, PartialEq)]
pub struct ClassifierFeedSump {
    residence: Time,
    provenance: Provenance,
}
impl ClassifierFeedSump {
    /// Keep the positive drain time and its original source/assumption together.
    /// Zero residence is a different model: use the existing instantaneous API.
    pub fn new(residence: Time, provenance: Provenance) -> Result<Self> {
        require(
            residence.value().is_finite()
                && residence.value() > 0.0
                && !provenance.source.trim().is_empty()
                && !provenance.license.trim().is_empty(),
            "sump requires positive finite residence and explicit provenance",
        )?;
        Ok(Self {
            residence,
            provenance,
        })
    }
    pub const fn residence(&self) -> Time {
        self.residence
    }
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }
}

/// Separate size-resolved physical inventories. Carry both into the next call.
#[derive(Clone, Debug, PartialEq)]
pub struct SumpCircuitInventory {
    pub mill: MillInventory,
    pub sump: MillInventory,
}

/// Admitted finite response with source owners and all three physical boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct SumpCircuitStep {
    model: ContinuousGrindingCircuit,
    sump_model: ClassifierFeedSump,
    initial: SumpCircuitInventory,
    fresh: ParticleStream,
    pub makeup_water: MassFlowRate,
    pub mill_residence: Time,
    pub elapsed: Time,
    pub mill_conditions: QueryPoint,
    pub classifier_conditions: QueryPoint,
    pub underflow_recycle: Partition,
    pub initial_point: TransientCircuitPoint,
    pub final_point: TransientCircuitPoint,
    pub inventory: SumpCircuitInventory,
    pub integrated_fresh: MillInventory,
    /// Water entering the sump directly; absent from fresh mill feed.
    pub integrated_makeup_water: Mass,
    pub integrated_mill_feed: MillInventory,
    pub integrated_mill_discharge: MillInventory,
    /// Actual mill discharge plus independently integrated makeup water.
    pub integrated_sump_feed: MillInventory,
    /// Actual sump discharge, also the original classifier's inlet.
    pub integrated_classifier_feed: MillInventory,
    pub integrated_underflow: MillInventory,
    pub integrated_recycle: MillInventory,
    pub integrated_product: MillInventory,
    pub integrated_purge: MillInventory,
    pub mill_balance: TransientMillBalance,
    pub sump_balance: TransientMillBalance,
    /// Both initial stocks + fresh + makeup - both final stocks - overflow - purge.
    pub boundary_balance: TransientMillBalance,
    pub peak_circulating_load: Dimensionless,
    pub settings: TransientCircuitSettings,
    pub substeps: usize,
}
impl SumpCircuitStep {
    pub fn model(&self) -> &ContinuousGrindingCircuit {
        &self.model
    }
    pub fn sump_model(&self) -> &ClassifierFeedSump {
        &self.sump_model
    }
    pub fn initial(&self) -> &SumpCircuitInventory {
        &self.initial
    }
    pub fn fresh_feed(&self) -> &ParticleStream {
        &self.fresh
    }
}

impl ContinuousGrindingCircuit {
    /// Advance both finite inventories with constant interval controls.
    ///
    /// Fresh dry flow must be positive. The initial sump must supply an admitted
    /// classifier feed; an undefined empty-sump solids fraction refuses. Sump
    /// makeup water does not dilute the mill immediately: only its returned share
    /// reaches that volume. Common residence implies no preferential settling.
    /// Complete inputs are immutable on domain, arithmetic, work or cancellation
    /// refusal. No source support between numerical endpoints is asserted.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_sump(
        &self,
        sump: &ClassifierFeedSump,
        initial: &SumpCircuitInventory,
        fresh: &ParticleStream,
        makeup_water: MassFlowRate,
        mill_residence: Time,
        elapsed: Time,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        underflow_recycle: Partition,
        settings: TransientCircuitSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SumpCircuitStep> {
        poll(&mut cancelled)?;
        let substeps = admit_settings(elapsed, settings.time)?;
        require(
            makeup_water.value().is_finite()
                && makeup_water.value() >= 0.0
                && mill_residence.value().is_finite()
                && mill_residence.value() > 0.0
                && settings.maximum_circulating_load.value().is_finite()
                && settings.maximum_circulating_load.value() >= 0.0
                && fresh.dry_solids()?.value() > 0.0,
            "require nonnegative makeup/circulation, positive residence and fresh dry flow",
        )?;
        for inventory in [&initial.mill, &initial.sump] {
            require(
                inventory.classes().len() == fresh.classes().len()
                    && inventory
                        .classes()
                        .iter()
                        .zip(fresh.classes())
                        .all(|(a, b)| {
                            a.id == b.id()
                                && a.diameter == b.diameter()
                                && a.elements.keys().eq(b.elements().keys())
                        }),
                "mill, sump and fresh require identical class and explicit element layouts",
            )?;
        }
        let splits: Vec<_> = self
            .mill
            .classes()
            .iter()
            .map(|c| (c.id.clone(), underflow_recycle))
            .collect();
        let point = |state: &SumpCircuitInventory| {
            self.transient_point_from_streams(
                state.mill.discharge_at(mill_residence)?,
                state.sump.discharge_at(sump.residence)?,
                fresh,
                mill_residence,
                mill_conditions,
                classifier_conditions,
                &splits,
                underflow_recycle,
                settings,
            )
        };
        let initial_point = point(initial)?;
        let mut final_point = initial_point.clone();
        let mut peak = initial_point.circulating_load.value();
        let returned = underflow_recycle.product().value();
        let purged = underflow_recycle.reject().value();
        let partitions = self
            .mill
            .classes()
            .iter()
            .map(|c| self.classifier.curve_partition(c.diameter))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let returns = partitions
            .iter()
            .map(|p| product(returned, p.product().value()))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let exits = partitions
            .iter()
            .map(|p| escape(returned, purged, p.reject().value()))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let water_return = product(
            returned,
            self.classifier.water_partition().product().value(),
        )?;
        let water_exit = escape(
            returned,
            purged,
            self.classifier.water_partition().reject().value(),
        )?;
        let elements: Vec<_> = initial.mill.classes()[0].elements.keys().cloned().collect();
        let channels = elements.len() + 1;
        let values = |inventory: &MillInventory| {
            inventory
                .classes()
                .iter()
                .map(|c| {
                    core::iter::once(c.dry_solids.value())
                        .chain(elements.iter().map(|e| c.elements[e].value()))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let mut masses = values(&initial.mill);
        let mut stored = values(&initial.sump);
        let fresh_rates: Vec<Vec<f64>> = fresh
            .classes()
            .iter()
            .map(|c| {
                core::iter::once(c.dry_solids().value())
                    .chain(elements.iter().map(|e| c.elements()[e].value()))
                    .collect()
            })
            .collect();
        let mut rates = fresh_rates.clone();
        let mut effective_exit = vec![0.0; masses.len()];
        let mut water = initial.mill.water().value();
        let mut sump_water = initial.sump.water().value();
        let mut order: Vec<_> = (0..masses.len()).collect();
        order.sort_by(|&a, &b| {
            self.mill.classes()[b]
                .diameter
                .value()
                .total_cmp(&self.mill.classes()[a].diameter.value())
                .then_with(|| self.mill.classes()[a].id.cmp(&self.mill.classes()[b].id))
        });
        let mut accounts: Vec<_> = (0..8)
            .map(|_| Accumulator::new(masses.len(), channels))
            .collect();
        let integrated_fresh = MillInventory::integrate(fresh, elapsed)?;
        let integrated_makeup_water = Mass::new(product(makeup_water.value(), elapsed.value())?);
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
            let denominator = sum([sump.residence.value(), dt])?;
            let retained = quotient(sump.residence.value(), denominator)?;
            let passed = quotient(dt, denominator)?;
            for i in 0..masses.len() {
                effective_exit[i] = if returned == 0.0 {
                    1.0
                } else {
                    sum([retained, product(passed, exits[i])?])?
                };
                let from_sump = quotient(returns[i], denominator)?;
                for c in 0..channels {
                    rates[i][c] = sum([fresh_rates[i][c], product(from_sump, stored[i][c])?])?;
                }
            }
            let sump_water_input = sum([sump_water, product(dt, makeup_water.value())?])?;
            let water_feed = sum([
                fresh.water().value(),
                product(quotient(water_return, denominator)?, sump_water_input)?,
            ])?;
            let effective_water_exit = if returned == 0.0 {
                1.0
            } else {
                sum([retained, product(passed, water_exit)?])?
            };
            implicit_step(
                &self.mill,
                &mut masses,
                &rates,
                &mut water,
                water_feed,
                dt,
                mill_residence,
                &effective_exit,
                effective_water_exit,
                &order,
                &mut cancelled,
            )?;
            let drained = quotient(dt, mill_residence.value())?;
            for (row, mill) in stored.iter_mut().zip(&masses) {
                for (value, mass) in row.iter_mut().zip(mill) {
                    *value = product(retained, sum([*value, product(drained, *mass)?])?)?;
                }
            }
            sump_water = product(retained, sum([sump_water_input, product(drained, water)?])?)?;
            inventory = SumpCircuitInventory {
                mill: inventory_from_channels(&self.mill, &elements, masses.clone(), water)?,
                sump: inventory_from_channels(&self.mill, &elements, stored.clone(), sump_water)?,
            };
            final_point = point(&inventory)?;
            peak = peak.max(final_point.circulating_load.value());
            let external = ParticleStream::mix(&[&final_point.product, &final_point.purge])?;
            for (account, stream) in accounts.iter_mut().zip([
                &final_point.mill_feed,
                &final_point.mill_discharge,
                &final_point.classifier_feed,
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
        let mut next = || totals.next().expect("eight declared integral accounts");
        let integrated_mill_feed = next()?;
        let integrated_mill_discharge = next()?;
        let integrated_classifier_feed = next()?;
        let integrated_underflow = next()?;
        let integrated_recycle = next()?;
        let integrated_product = next()?;
        let integrated_purge = next()?;
        let integrated_external = next()?;
        let integrated_sump_feed = add_water(&integrated_mill_discharge, integrated_makeup_water)?;
        let integrated_external_feed = add_water(&integrated_fresh, integrated_makeup_water)?;
        let mill_balance = close_accounts(
            &initial.mill,
            &integrated_mill_feed,
            &inventory.mill,
            &integrated_mill_discharge,
            settings.time,
        )?;
        let sump_balance = close_accounts(
            &initial.sump,
            &integrated_sump_feed,
            &inventory.sump,
            &integrated_classifier_feed,
            settings.time,
        )?;
        let boundary_balance = close_accounts(
            &combined(&self.mill, &elements, initial)?,
            &integrated_external_feed,
            &combined(&self.mill, &elements, &inventory)?,
            &integrated_external,
            settings.time,
        )?;
        poll(&mut cancelled)?;
        Ok(SumpCircuitStep {
            model: self.clone(),
            sump_model: sump.clone(),
            initial: initial.clone(),
            fresh: fresh.clone(),
            makeup_water,
            mill_residence,
            elapsed,
            mill_conditions: mill_conditions.clone(),
            classifier_conditions: classifier_conditions.clone(),
            underflow_recycle,
            initial_point,
            final_point,
            inventory,
            integrated_fresh,
            integrated_makeup_water,
            integrated_mill_feed,
            integrated_mill_discharge,
            integrated_sump_feed,
            integrated_classifier_feed,
            integrated_underflow,
            integrated_recycle,
            integrated_product,
            integrated_purge,
            mill_balance,
            sump_balance,
            boundary_balance,
            peak_circulating_load: Dimensionless::new(peak),
            settings,
            substeps,
        })
    }
}

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut total = Sum::default();
    for value in values {
        total.add(value)?;
    }
    Ok(total.value()?)
}
fn add_water(inventory: &MillInventory, water: Mass) -> Result<MillInventory> {
    Ok(MillInventory::new(
        inventory.classes().to_vec(),
        Mass::new(sum([inventory.water().value(), water.value()])?),
    )?)
}
fn combined(
    model: &ContinuousPopulationBalanceMill,
    elements: &[ElementId],
    state: &SumpCircuitInventory,
) -> Result<MillInventory> {
    let values = state
        .mill
        .classes()
        .iter()
        .zip(state.sump.classes())
        .map(|(a, b)| {
            core::iter::once(sum([a.dry_solids.value(), b.dry_solids.value()]))
                .chain(
                    elements
                        .iter()
                        .map(|e| sum([a.elements[e].value(), b.elements[e].value()])),
                )
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(inventory_from_channels(
        model,
        elements,
        values,
        sum([state.mill.water().value(), state.sump.water().value()])?,
    )?)
}
