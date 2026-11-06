//! Finite radioactive material moving between the original thermal cells.
//!
//! One complete compartment cohort evolves reactions, daughter ingrowth,
//! movement and external production together. Local event energy heats the cell
//! in which that event occurs. Movement of the tracked nuclei does not imply
//! movement of bulk carrier mass or sensible enthalpy.

use std::{collections::BTreeSet, sync::Arc};

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NeutronFluxSchedule, NuclearHeating};
use crate::burnup::production::BurnupProductionSchedule;
use crate::burnup::transport::{BurnupTransportBudget, BurnupTransportState, BurnupTransportStep};
use crate::burnup::{BurnupBudget, BurnupStep};
use crate::{EnergyChannels, FissionError};

mod schedule;
pub use schedule::{NuclideTransportChange, NuclideTransportSchedule};

#[derive(Clone, Debug)]
struct Forcing {
    fluxes: Vec<NeutronFluxSchedule>,
    movement: Arc<NuclideTransportSchedule>,
    production: Vec<BurnupProductionSchedule>,
}

/// One handle into a complete immutable compartment history. Construct the
/// collection with bind or bind_with_schedules, then advance it together.
#[derive(Clone, Debug)]
pub struct BurnupTransportCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    compartment_index: usize,
    state: Arc<BurnupTransportState>,
    forcing: Arc<Forcing>,
}

/// Product-resolved radioactive movement and local heat on the native clock.
pub type BurnupTransportHeating = NuclearHeating<BurnupTransportCellSource>;

/// Local reaction report plus the one shared, fully accounted movement step.
#[derive(Clone, Debug)]
pub struct BurnupTransportCellStep {
    compartment_index: usize,
    transport: Arc<BurnupTransportStep>,
}

impl BurnupTransportCellStep {
    #[must_use]
    pub fn burnup(&self) -> &BurnupStep {
        &self.transport.burnup[self.compartment_index]
    }
    #[must_use]
    pub fn transport(&self) -> &BurnupTransportStep {
        &self.transport
    }
    #[must_use]
    pub fn compartment_index(&self) -> usize {
        self.compartment_index
    }
}

impl BurnupTransportCellSource {
    /// Bind all compartments to distinct original cells, using the state's
    /// current movement rates and zero external nuclide production. Each flux
    /// history remains explicit, including an authored zero-flux history.
    pub fn bind(
        state: Arc<BurnupTransportState>,
        cells: Vec<(CellId, MaterialId)>,
        fluxes: Vec<NeutronFluxSchedule>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        let movement = Arc::new(NuclideTransportSchedule::constant(&state));
        let production = state
            .inventories()
            .iter()
            .map(|inventory| {
                BurnupProductionSchedule::new(
                    Arc::clone(inventory.model()),
                    vec![0.0; inventory.populations().len()],
                    Vec::new(),
                    0,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::bind_with_schedules(state, cells, fluxes, movement, production)
    }

    /// Bind independent absolute-time flux, movement and production histories.
    /// Their union clips the existing thermal step exactly. A binding created
    /// after time zero uses its supplied current inventories without replaying
    /// earlier schedules. Compartment volumes must match the thermal geometry.
    pub fn bind_with_schedules(
        state: Arc<BurnupTransportState>,
        cells: Vec<(CellId, MaterialId)>,
        fluxes: Vec<NeutronFluxSchedule>,
        movement: Arc<NuclideTransportSchedule>,
        production: Vec<BurnupProductionSchedule>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        let count = state.inventories().len();
        if cells.len() != count || fluxes.len() != count || production.len() != count {
            return Err(FissionThermalError::Invalid(
                "complete movement-to-thermal cells, spectra and sources required",
            ));
        }
        movement.validate_binding(&state)?;
        for ((inventory, flux), source) in state.inventories().iter().zip(&fluxes).zip(&production)
        {
            if flux.flux_at(0.0)?.groups() != inventory.model().groups()
                || !Arc::ptr_eq(source.model(), inventory.model())
            {
                return Err(FissionThermalError::Invalid(
                    "movement histories must retain original nuclides and neutron groups",
                ));
            }
        }
        let forcing = Arc::new(Forcing {
            fluxes,
            movement,
            production,
        });
        let mut used = BTreeSet::new();
        let mut sources = Vec::with_capacity(count);
        for (compartment_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !used.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct movement thermal cells required",
                ));
            }
            sources.push(Self {
                cell_id,
                material_id,
                compartment_index,
                state: Arc::clone(&state),
                forcing: Arc::clone(&forcing),
            });
        }
        Ok(sources)
    }

    #[must_use]
    pub fn state(&self) -> &Arc<BurnupTransportState> {
        &self.state
    }
    #[must_use]
    pub fn compartment_index(&self) -> usize {
        self.compartment_index
    }
    #[must_use]
    pub fn flux_schedule(&self) -> &NeutronFluxSchedule {
        &self.forcing.fluxes[self.compartment_index]
    }
    #[must_use]
    pub fn production_schedule(&self) -> &BurnupProductionSchedule {
        &self.forcing.production[self.compartment_index]
    }
    #[must_use]
    pub fn transport_schedule(&self) -> &Arc<NuclideTransportSchedule> {
        &self.forcing.movement
    }

    pub(super) fn shares_forcing(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.forcing, &other.forcing)
    }
}

impl sealed::Sealed for BurnupTransportCellSource {}
impl NuclearHeatSource for BurnupTransportCellSource {
    type Budget = BurnupTransportBudget;
    type Step = BurnupTransportCellStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn expected_volume_m3(&self) -> Option<f64> {
        Some(self.state.volumes_m3()[self.compartment_index])
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }
    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(budget.nuclei_tolerance)
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        [
            self.flux_schedule().next_after(time_s),
            self.production_schedule().next_after(time_s),
            self.transport_schedule().next_after(time_s),
        ]
        .into_iter()
        .flatten()
        .min_by(f64::total_cmp)
    }
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        if self.compartment_index == 0 {
            let production = self
                .forcing
                .production
                .iter()
                .map(|schedule| schedule.rates_at(time_s).map(<[f64]>::to_vec))
                .collect::<Result<Vec<_>, _>>()?;
            self.state.validate_inputs(
                self.forcing.movement.transport_at(time_s)?,
                &production,
                budget.limits,
                || Ok(()),
            )?;
        }
        self.state.inventories()[self.compartment_index].trial_with_production(
            0.0,
            self.flux_schedule().flux_at(time_s)?,
            temperature_k,
            self.production_schedule().rates_at(time_s)?,
            BurnupBudget {
                limits: budget.limits.burnup,
                nuclei_tolerance: budget.nuclei_tolerance,
                energy_tolerance_j: budget.energy_tolerance_j,
            },
            || Ok(()),
        )?;
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.state.model().validate_temperature(temperature_k)?;
        Ok(())
    }
    fn validate_sources<'source>(
        mut sources: impl ExactSizeIterator<Item = &'source Self>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionThermalError>,
    ) -> Result<(), FissionThermalError>
    where
        Self: 'source,
    {
        checkpoint()?;
        let count = sources.len();
        let first = sources.next().ok_or(FissionThermalError::Invalid(
            "movement collection must be nonempty",
        ))?;
        if count != first.state.inventories().len() {
            return Err(FissionThermalError::Invalid(
                "all movement compartments must advance together",
            ));
        }
        let mut seen = vec![false; count];
        for source in std::iter::once(first).chain(sources) {
            checkpoint()?;
            if !Arc::ptr_eq(&first.state, &source.state)
                || !first.shares_forcing(source)
                || source.compartment_index >= count
                || seen[source.compartment_index]
            {
                return Err(FissionThermalError::Invalid(
                    "movement cells must share one complete accepted history and schedules",
                ));
            }
            seen[source.compartment_index] = true;
        }
        Ok(())
    }
    fn trial_source(
        &self,
        _: f64,
        _: f64,
        _: f64,
        _: Self::Budget,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        Err(FissionThermalError::Invalid(
            "nuclide movement requires the complete compartment collection",
        ))
    }
    fn trial_sources(
        sources: &[&Self],
        time_s: f64,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(Self, Self::Step)>, FissionThermalError> {
        Self::validate_sources(sources.iter().copied(), &mut || {
            checkpoint().map_err(Into::into)
        })?;
        if temperatures_k.len() != sources.len() {
            return Err(FissionThermalError::Invalid(
                "one temperature per movement compartment required",
            ));
        }
        let first = sources[0];
        let mut temperatures = vec![0.0; sources.len()];
        for (source, &temperature) in sources.iter().zip(temperatures_k) {
            checkpoint()?;
            if source
                .next_transition(time_s)
                .is_some_and(|next| duration_s > next - time_s)
            {
                return Err(FissionThermalError::Invalid(
                    "movement trial crosses a scheduled flux, flow or source change",
                ));
            }
            temperatures[source.compartment_index] = temperature;
        }
        let mut fluxes = Vec::with_capacity(sources.len());
        let mut production = Vec::with_capacity(sources.len());
        for (flux, source) in first.forcing.fluxes.iter().zip(&first.forcing.production) {
            checkpoint()?;
            fluxes.push(flux.flux_at(time_s)?.clone());
            production.push(source.rates_at(time_s)?.to_vec());
        }
        let (state, report) = first.state.trial_with_inputs(
            duration_s,
            &fluxes,
            &temperatures,
            first.forcing.movement.transport_at(time_s)?,
            &production,
            budget,
            &mut *checkpoint,
        )?;
        let state = Arc::new(state);
        let report = Arc::new(report);
        let mut candidates = Vec::with_capacity(sources.len());
        for source in sources {
            checkpoint()?;
            let mut candidate = (*source).clone();
            candidate.state = Arc::clone(&state);
            candidates.push((
                candidate,
                BurnupTransportCellStep {
                    compartment_index: source.compartment_index,
                    transport: Arc::clone(&report),
                },
            ));
        }
        Ok(candidates)
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.burnup().released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.burnup().energy_roundoff_j
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        step.burnup().population_roundoff_nuclei
    }
    fn decays(step: &Self::Step) -> f64 {
        step.burnup().decays
    }
    fn produced_nuclei(step: &Self::Step) -> f64 {
        step.transport.produced_nuclei_by_compartment[step.compartment_index]
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.burnup().fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        step.burnup().captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.inventories()[self.compartment_index].cumulative_released())
    }
}
