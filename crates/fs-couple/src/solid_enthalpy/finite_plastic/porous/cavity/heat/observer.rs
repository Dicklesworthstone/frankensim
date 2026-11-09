//! Explicit gas-wall links on the original observer/thermal timeline.
use super::{
    CavityHeatError, CavityHeatObserverError, CavityWallHeatLink,
    CavityWallHeatReceipt, CavityWallSourceReceipt, ThermalCavityGas, number,
};
use super::source::{SourceFlux, SourceSample, observer_error};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
    collect_source_rates,
};
use std::collections::BTreeMap;

fn finite<E>(value: f64) -> Result<f64, CavityHeatObserverError<E>> {
    if value.is_finite() { Ok(value) }
    else { Err(CavityHeatObserverError::Invalid("unrepresentable gas-wall heat flux")) }
}

fn integral<E>(power: f64, dt: f64) -> Result<Dd, CavityHeatObserverError<E>> {
    let energy = Dd::from_f64(power) * Dd::from_f64(dt);
    let rounded = finite(energy.to_f64())?;
    if power != 0.0 && dt != 0.0 && rounded == 0.0 {
        return Err(CavityHeatObserverError::Invalid("gas-wall heat integral underflows"));
    }
    Ok(energy)
}

#[derive(Clone)]
struct LinkSample {
    heat_j: Dd,
    temperature_k: f64,
    power_w: f64,
}

#[derive(Clone)]
pub(super) struct WallObserver<'s, 'm, O, G> {
    links: &'s [CavityWallHeatLink],
    gas_temperatures: &'s [f64],
    source: Option<SourceFlux<'s, 'm, G>>,
    pub(super) inner: O,
    samples: Vec<LinkSample>,
    source_samples: Vec<SourceSample>,
    source_heat_j: Vec<Dd>,
    applied_w: BTreeMap<CellId, Dd>,
    conductances_w_k: BTreeMap<CellId, f64>,
    pub(super) gas_heat_j: Vec<Dd>,
    pub(super) applied_heat_j: Dd,
    pub(super) interface_defect_j: Dd,
}

impl<'s, 'm, O, G> WallObserver<'s, 'm, O, G> {
    pub(super) fn new(
        links: &'s [CavityWallHeatLink], gas_temperatures: &'s [f64],
        source: Option<SourceFlux<'s, 'm, G>>, inner: O,
    ) -> Self {
        let source_heat_j = vec![Dd::ZERO; source.as_ref().map_or(0, SourceFlux::len)];
        Self {
            links, gas_temperatures, inner,
            source, source_samples: Vec::new(), source_heat_j,
            samples: links.iter().map(|_| LinkSample {
                heat_j: Dd::ZERO, temperature_k: 0.0, power_w: 0.0,
            }).collect(),
            applied_w: BTreeMap::new(), conductances_w_k: BTreeMap::new(),
            gas_heat_j: vec![Dd::ZERO; gas_temperatures.len()],
            applied_heat_j: Dd::ZERO, interface_defect_j: Dd::ZERO,
        }
    }

    pub(super) fn receipts<E>(&self) -> Result<Vec<CavityWallHeatReceipt>, CavityHeatError<E>> {
        self.links.iter().zip(&self.samples).map(|(link, sample)| {
            Ok(CavityWallHeatReceipt {
                id: link.id, cavity: link.cavity, cell: link.cell,
                heat_into_solid_j: number(sample.heat_j.to_f64())?,
                last_cell_temperature_k: sample.temperature_k,
                last_power_into_solid_w: sample.power_w,
            })
        }).collect()
    }

    pub(super) fn source_receipts<E>(&self) -> Result<Vec<CavityWallSourceReceipt>, CavityHeatError<E>> {
        self.source_samples.iter().zip(&self.source_heat_j).map(|(sample, heat)| {
            let mut receipt = sample.receipt.clone();
            receipt.heat_into_solid_j = number(heat.to_f64())?;
            Ok(receipt)
        }).collect()
    }
}

impl<O: DepositionObserver, G: ThermalCavityGas> DepositionObserver for WallObserver<'_, '_, O, G> {
    type Error = CavityHeatObserverError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        for sample in &mut self.samples {
            sample.heat_j = Dd::ZERO;
            sample.temperature_k = 0.0;
            sample.power_w = 0.0;
        }
        self.applied_w.clear();
        self.source_samples.clear();
        self.source_heat_j.fill(Dd::ZERO);
        self.conductances_w_k.clear();
        self.gas_heat_j.fill(Dd::ZERO);
        self.applied_heat_j = Dd::ZERO;
        self.interface_defect_j = Dd::ZERO;
        self.inner.begin().map_err(CavityHeatObserverError::Observer)
    }

    fn sources(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, checkpoint))
    }
    fn visit_sources(
        &mut self, thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner.visit_sources(thermal, emit, checkpoint)
            .map_err(CavityHeatObserverError::Observer)?;
        let mut powers = BTreeMap::<CellId, Dd>::new();
        let mut conductances = BTreeMap::<CellId, Dd>::new();
        for (link, sample) in self.links.iter().zip(&mut self.samples) {
            checkpoint()?;
            let temperature = thermal.temperature_k(link.cell)?;
            let contrast = finite(self.gas_temperatures[link.cavity] - temperature)?;
            let power = finite(link.conductance_w_k * contrast)?;
            if link.conductance_w_k > 0.0 && contrast != 0.0 && power == 0.0 {
                return Err(CavityHeatObserverError::Invalid("gas-wall heat power underflows"));
            }
            emit(link.cell, power)?;
            let total = powers.entry(link.cell).or_insert(Dd::ZERO);
            *total = *total + Dd::from_f64(power);
            let decay = conductances.entry(link.cell).or_insert(Dd::ZERO);
            *decay = *decay + Dd::from_f64(link.conductance_w_k);
            sample.temperature_k = temperature;
            sample.power_w = power;
        }
        if let Some(source) = &self.source {
            let samples = source.samples(thermal, &mut || checkpoint()).map_err(observer_error)?;
            for sample in &samples {
                checkpoint()?;
                let row = &sample.receipt;
                emit(row.cell, row.last_power_into_solid_w)?;
                let total = powers.entry(row.cell).or_insert(Dd::ZERO);
                *total = *total + Dd::from_f64(row.last_power_into_solid_w);
                let decay = conductances.entry(row.cell).or_insert(Dd::ZERO);
                *decay = *decay + Dd::from_f64(row.last_explicit_decay_conductance_w_k);
            }
            self.source_samples = samples;
        }
        for sum in powers.values() {
            checkpoint()?;
            finite(sum.to_f64())?;
        }
        // Per-cell reporting retains the same contributions sent individually
        // to the thermal owner; none are flattened beside an unrelated source.
        self.applied_w = powers;
        self.conductances_w_k = conductances.into_iter().map(|(cell, sum)| {
            Ok((cell, finite(sum.to_f64())?))
        }).collect::<Result<_, Self::Error>>()?;
        Ok(())
    }

    fn conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner.conductances(thermal, checkpoint).map_err(CavityHeatObserverError::Observer)
    }
    fn reservoir_conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.reservoir_conductances(thermal, checkpoint).map_err(CavityHeatObserverError::Observer)
    }
    fn ambient_conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.ambient_conductances(thermal, checkpoint).map_err(CavityHeatObserverError::Observer)
    }
    fn heat_transfers(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner.heat_transfers(thermal, checkpoint).map_err(CavityHeatObserverError::Observer)
    }
    fn limit_timestep_s(
        &mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        let proposed = self.inner.limit_timestep_s(thermal, proposed_s, checkpoint)
            .map_err(CavityHeatObserverError::Observer)?;
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed, &self.conductances_w_k, checkpoint,
        )?)
    }
    fn requires_step_interval(&self) -> bool { self.inner.requires_step_interval() }
    fn next_transition_time_s(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner.next_transition_time_s(thermal, checkpoint).map_err(CavityHeatObserverError::Observer)
    }
    fn step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start: f64, dt: f64, end: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start, dt, end, emit, checkpoint)
        })
    }
    fn visit_step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start: f64, dt: f64, end: f64,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner.visit_step_sources(thermal, start, dt, end, emit, checkpoint)
            .map_err(CavityHeatObserverError::Observer)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.inner.event(event, thermal).map_err(CavityHeatObserverError::Observer)
    }
    fn step(
        &mut self, start: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner.step(start, dt, old, thermal, checkpoint)
            .map_err(CavityHeatObserverError::Observer)?;
        let mut applied = Dd::ZERO;
        for &power in self.applied_w.values() {
            checkpoint()?;
            applied = applied + power * Dd::from_f64(dt);
            finite(applied.to_f64())?;
        }
        let mut exchanged = Dd::ZERO;
        for (link, sample) in self.links.iter().zip(&mut self.samples) {
            checkpoint()?;
            let heat = integral(sample.power_w, dt)?;
            exchanged = exchanged + heat;
            sample.heat_j = sample.heat_j + heat;
            self.gas_heat_j[link.cavity] = self.gas_heat_j[link.cavity] + heat;
            finite(sample.heat_j.to_f64())?;
            finite(self.gas_heat_j[link.cavity].to_f64())?;
        }
        for (sample, accumulated) in self.source_samples.iter().zip(&mut self.source_heat_j) {
            checkpoint()?;
            let row = &sample.receipt;
            let heat = integral(row.last_power_into_solid_w, dt)?;
            exchanged = exchanged + heat;
            *accumulated = *accumulated + heat;
            self.gas_heat_j[row.cavity] = self.gas_heat_j[row.cavity] + heat;
            finite(accumulated.to_f64())?;
            finite(self.gas_heat_j[row.cavity].to_f64())?;
        }
        self.applied_heat_j = self.applied_heat_j + applied;
        // Compare the delivered per-cell integral with independent per-patch
        // gas accounting. Final source/storage rounding belongs to thermal.
        self.interface_defect_j = self.interface_defect_j + Dd::from_f64(
            finite((applied - exchanged).to_f64())?.abs(),
        );
        finite(self.applied_heat_j.to_f64())?;
        finite(self.interface_defect_j.to_f64())?;
        Ok(())
    }
}
