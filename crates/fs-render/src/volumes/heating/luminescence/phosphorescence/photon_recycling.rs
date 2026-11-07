//! Photon-number re-excitation through the original directed optical snapshot.
//!
//! Only explicitly mapped receivers acquire spin-state stores. All other
//! absorption stays passive; every activated photon pays its destination-state
//! energy before any excess becomes receiver heat. The native positive network
//! integrates crossing and arbitrarily many recycled visits in one interval.

use super::super::absorption::AbsorbingLuminescenceEnvironment;
use super::super::{Cell, LuminescenceSource, LuminescentOpticalHeating};
use super::annihilation::NetworkAnnihilation;
use super::{
    CoupledStep, HC_J_NM, SingletTripletFeedback, SingletTripletLimits, finite, invalid,
    preserve_checkpoint, sum,
};
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_thermal::luminescence::transfer::SingletTripletDecay;
use fs_thermal::luminescence::transfer::recycling::{
    EmissionChannel, PhotonAbsorption, SingletTripletRecyclingPlan,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Explicit excitation branching conditional on absorption of one photon.
/// Nodes follow the original common singlet/triplet emission wavelength order.
/// Each pair contains [singlet, triplet] probabilities and sums to at most one.
/// The unactivated fraction, and any photon energy above the selected state,
/// become heat at this receiver. These probabilities are authored inputs;
/// neither a long lifetime nor optical absorption alone establishes them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhosphorescenceReceiver<const M: usize> {
    pub grid_cell: usize,
    pub excitation_fractions: [[f64; 2]; M],
}

#[derive(Debug)]
struct Topology {
    channel_bins: Vec<usize>,
    absorption_receivers: Vec<CellId>,
    fraction_roundoff: f64,
    limits: SingletTripletLimits,
    identity: ContentHash,
}

#[derive(Debug)]
struct Prepared {
    duration_s: f64,
    plan: SingletTripletRecyclingPlan,
}

/// Shared immutable topology and propagator only; physical stores stay in the
/// existing accepted/tentative optical owner and are never cached here.
#[derive(Clone, Debug)]
pub(super) struct PhotonRecyclingFeedback {
    topology: Arc<Topology>,
    prepared: Arc<Prepared>,
}

#[derive(Debug)]
pub(in super::super) struct CoupledRecyclingStep<const M: usize> {
    pub cells: Vec<CoupledStep<M>>,
    pub recycled_excitation_j: Vec<f64>,
    pub secondary_heat_j: BTreeMap<CellId, f64>,
    pub escaped_j: [f64; M],
    pub unresolved_j: [f64; M],
    /// Transport and receiver-partition defects only. The parent separately
    /// charges each cell's kinetic balance, so the native global residual is
    /// not charged again as an additional local defect.
    pub absolute_roundoff_j: f64,
}

impl<const N: usize, const M: usize> LuminescentOpticalHeating<N, M> {
    /// Recycle emitted photons through explicitly admitted spin-state channels.
    /// Bind after the singlet/triplet source and absorbing optical environment,
    /// before optional annihilation binding and the first native advance.
    /// Missing receiver mappings stay passive.
    /// A mapped, initially unilluminated receiver acquires the same sourced
    /// singlet/triplet material law, refreshed at its own native temperature.
    ///
    /// Every positive excitation branch requires photon energy at least equal
    /// to that destination state's energy. No missing energy is silently drawn
    /// from the bath, and no branch or transport fraction is renormalized.
    /// The existing native owner admits at most 64 states (32 active cells);
    /// passive absorbers do not consume excitation states.
    pub fn with_singlet_triplet_recycling(
        mut self,
        receivers: &[PhosphorescenceReceiver<M>],
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self.recycling.is_some()
            || self.singlet_triplet.as_ref().is_none_or(|coupling| {
                coupling.photon_recycling.is_some() || coupling.annihilation.is_some()
            })
        {
            return Err(invalid(
                "phosphorescent recycling requires unadvanced coupling before annihilation binding",
            ));
        }
        self.check_state(thermal, true, &mut checkpoint)?;
        let environment = self.emission_transport.as_ref().ok_or_else(|| {
            invalid("phosphorescent recycling requires original photon transport")
        })?;
        let bindings: BTreeMap<_, _> = environment.bindings().collect();
        let mut authored = BTreeMap::new();
        for receiver in receivers {
            checkpoint()?;
            let thermal_cell = *bindings
                .get(&receiver.grid_cell)
                .ok_or(OpticalHeatingError::UnmappedCell(receiver.grid_cell))?;
            environment.check_emitter(receiver.grid_cell, thermal_cell)?;
            for branches in receiver.excitation_fractions {
                checkpoint()?;
                if branches
                    .iter()
                    .any(|p| !p.is_finite() || !(0.0..=1.0).contains(p))
                    || sum(branches)? > 1.0
                {
                    return Err(invalid(
                        "one absorbed photon may create at most one sourced excitation",
                    ));
                }
            }
            if authored
                .insert(receiver.grid_cell, receiver.excitation_fractions)
                .is_some()
            {
                return Err(OpticalHeatingError::DuplicateCell(receiver.grid_cell));
            }
        }
        let mut active: BTreeMap<_, _> = self
            .cells
            .iter()
            .cloned()
            .map(|cell| (cell.grid_cell, cell))
            .collect();
        for &grid_cell in authored.keys() {
            checkpoint()?;
            active.entry(grid_cell).or_insert_with(|| Cell {
                grid_cell,
                thermal_cell: bindings[&grid_cell],
                reservoir_pump_w: [0.0; N],
                prompt_heat_w: [0.0; N],
                absorbed_w: [0.0; N],
                stored_j: 0.0,
                state_stored_j: [0.0; 2],
            });
        }
        let coupling = self
            .singlet_triplet
            .as_ref()
            .expect("admitted singlet/triplet source");
        if active.is_empty()
            || active
                .len()
                .checked_mul(2)
                .is_none_or(|states| states > coupling.source.limits.maximum_nodes.min(64))
        {
            return Err(invalid(
                "phosphorescent recycling exceeds the native excitation-state budget",
            ));
        }
        let cells: Vec<_> = active.into_values().collect();
        let mut models = Vec::with_capacity(cells.len());
        for cell in &cells {
            checkpoint()?;
            environment.check_emitter(cell.grid_cell, cell.thermal_cell)?;
            let model = coupling
                .model_at_temperature(&self.source, thermal.temperature_k(cell.thermal_cell)?)?;
            if let Some(branches) = authored.get(&cell.grid_cell) {
                let energies = model.excitation_energies_j();
                for (lane, probabilities) in branches.iter().enumerate() {
                    checkpoint()?;
                    let photon_j = HC_J_NM / self.source.fluorescence().wavelengths_nm()[lane];
                    for state in 0..2 {
                        if probabilities[state] > 0.0 && photon_j < energies[state] {
                            return Err(invalid(
                                "recycled photon cannot pay for its selected spin state",
                            ));
                        }
                    }
                }
            }
            models.push(model);
        }
        let feedback = PhotonRecyclingFeedback::prepare(
            &self.source,
            coupling,
            &cells,
            &models,
            environment,
            &authored,
            &mut checkpoint,
        )?;
        checkpoint()?;
        self.cells = cells;
        self.cell_emission_w = self
            .cells
            .iter()
            .map(|cell| (cell.grid_cell, [0.0; M]))
            .collect();
        self.mean_cell_emission_w = self.cell_emission_w.clone();
        let coupling = self
            .singlet_triplet
            .as_mut()
            .expect("admitted singlet/triplet source");
        coupling.photon_recycling = Some(feedback);
        #[cfg(feature = "material-state")]
        {
            coupling.maximum_cached_temperatures = self.cells.len().saturating_mul(2).max(1);
        }
        Ok(self)
    }

    /// Whether emitted photons can explicitly repopulate the spin-state stores.
    pub fn has_singlet_triplet_recycling(&self) -> bool {
        self.singlet_triplet
            .as_ref()
            .is_some_and(|coupling| coupling.photon_recycling.is_some())
    }
}

impl<const M: usize> SingletTripletFeedback<M> {
    /// Evolve the complete optical feedback network once, before local native
    /// heat sources are reduced. Preparation is keyed by every actual model,
    /// including temperature-refreshed material rates, and the interval length.
    pub(in super::super) fn advance_recycling<const N: usize>(
        &mut self,
        singlet: &LuminescenceSource<M>,
        cells: &[Cell<N>],
        scale: &[f64; N],
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<CoupledRecyclingStep<M>>, OpticalHeatingError> {
        if self.photon_recycling.is_none() {
            return Ok(None);
        }
        let mut models = Vec::with_capacity(cells.len());
        let mut populations = Vec::with_capacity(cells.len());
        let mut pump_rates = Vec::with_capacity(cells.len());
        for cell in cells {
            checkpoint()?;
            let model =
                self.model_at_temperature(singlet, thermal.temperature_k(cell.thermal_cell)?)?;
            let energies = model.excitation_energies_j();
            populations.push([
                finite(cell.state_stored_j[0] / energies[0])?,
                finite(cell.state_stored_j[1] / energies[1])?,
            ]);
            let mut pump = CheckedSum::default();
            for (lane, factor) in scale.iter().enumerate() {
                pump.add(finite(cell.reservoir_pump_w[lane] * factor)?)?;
            }
            pump_rates.push([finite(pump.value()? / energies[0])?, 0.0]);
            models.push(model);
        }
        let annihilation = self
            .annihilation
            .as_ref()
            .and_then(|feedback| feedback.network_parameters(cells));
        self.photon_recycling
            .as_mut()
            .expect("bound photon recycling")
            .advance::<M>(
                &models,
                &populations,
                &pump_rates,
                duration_s,
                annihilation.as_ref(),
                checkpoint,
            )
            .map(Some)
    }
}

impl PhotonRecyclingFeedback {
    #[allow(clippy::too_many_arguments)]
    fn prepare<const N: usize, const M: usize>(
        singlet: &LuminescenceSource<M>,
        coupling: &SingletTripletFeedback<M>,
        cells: &[Cell<N>],
        models: &[SingletTripletDecay],
        environment: &AbsorbingLuminescenceEnvironment<M>,
        authored: &BTreeMap<usize, [[f64; 2]; M]>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        let operator = environment.operator();
        let indices: BTreeMap<_, _> = cells
            .iter()
            .enumerate()
            .map(|(index, cell)| (cell.grid_cell, index))
            .collect();
        let bindings: BTreeMap<_, _> = environment.bindings().collect();
        let mut channels = Vec::new();
        let mut absorption = Vec::new();
        let mut channel_bins = Vec::new();
        let mut absorption_receivers = Vec::new();
        let maximum_links = coupling.source.limits.maximum_edges;
        for (index, cell) in cells.iter().enumerate() {
            checkpoint()?;
            let row = &operator.rows[operator
                .rows
                .binary_search_by_key(&cell.grid_cell, |row| row.emitter_cell)
                .map_err(|_| invalid("missing original recycling emitter row"))?];
            for (state, source) in [singlet, coupling.source.triplet()].into_iter().enumerate() {
                for (lane, &weight) in source.fluorescence().photon_fractions().iter().enumerate() {
                    checkpoint()?;
                    if weight == 0.0 {
                        // A shared-grid node unused by this state is not an
                        // emitted photon and cannot violate its state energy.
                        continue;
                    }
                    if channels.len() >= maximum_links {
                        return Err(invalid("phosphorescent photon-channel budget exceeded"));
                    }
                    let channel = channels.len();
                    channels.push(EmissionChannel {
                        source_cell: index,
                        source_state: state,
                        photon_energy_j: finite(HC_J_NM / operator.wavelengths_nm[lane])?,
                        photon_fraction: weight,
                        escape_fraction: row.escape.mean[lane],
                        unresolved_fraction: row.unresolved.mean[lane],
                    });
                    channel_bins.push(lane);
                    for link in &row.absorption {
                        checkpoint()?;
                        let fraction = link.fraction.mean[lane];
                        if fraction == 0.0 {
                            continue;
                        }
                        if absorption.len() >= maximum_links {
                            return Err(invalid("phosphorescent absorption-link budget exceeded"));
                        }
                        let (receiver_cell, excitation_fractions) =
                            match authored.get(&link.cell_index) {
                                Some(branches) => (Some(indices[&link.cell_index]), branches[lane]),
                                None => (None, [0.0; 2]),
                            };
                        absorption.push(PhotonAbsorption {
                            channel,
                            receiver_cell,
                            absorbed_fraction: fraction,
                            excitation_fractions,
                        });
                        absorption_receivers.push(
                            *bindings
                                .get(&link.cell_index)
                                .ok_or(OpticalHeatingError::UnmappedCell(link.cell_index))?,
                        );
                    }
                }
            }
        }
        let fraction_roundoff = environment.tolerance().fraction_roundoff;
        let limits = coupling.source.limits;
        let plan = preserve_checkpoint(checkpoint, |cancelled| {
            SingletTripletRecyclingPlan::prepare(
                models,
                &channels,
                &absorption,
                0.0,
                fraction_roundoff,
                limits,
                cancelled,
            )
        })?;
        let identity = snapshot_identity(environment, authored, checkpoint)?;
        Ok(Self {
            topology: Arc::new(Topology {
                channel_bins,
                absorption_receivers,
                fraction_roundoff,
                limits,
                identity,
            }),
            prepared: Arc::new(Prepared {
                duration_s: 0.0,
                plan,
            }),
        })
    }

    pub(super) fn identity(&self) -> ContentHash {
        self.topology.identity
    }

    fn advance<const M: usize>(
        &mut self,
        models: &[SingletTripletDecay],
        populations: &[[f64; 2]],
        pump_rates: &[[f64; 2]],
        duration_s: f64,
        annihilation: Option<&NetworkAnnihilation>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledRecyclingStep<M>, OpticalHeatingError> {
        let same_models = self.prepared.plan.models() == models;
        if !same_models || self.prepared.duration_s != duration_s {
            let plan = preserve_checkpoint(checkpoint, |cancelled| {
                if same_models {
                    self.prepared.plan.with_duration(duration_s, cancelled)
                } else {
                    SingletTripletRecyclingPlan::prepare(
                        models,
                        self.prepared.plan.channels(),
                        self.prepared.plan.absorption_edges(),
                        duration_s,
                        self.topology.fraction_roundoff,
                        self.topology.limits,
                        cancelled,
                    )
                }
            })?;
            self.prepared = Arc::new(Prepared { duration_s, plan });
        }
        let (native, annihilation_heat_j) = match annihilation {
            Some(annihilation) => {
                let step = preserve_checkpoint(checkpoint, |cancelled| {
                    self.prepared.plan.advance_with_annihilation(
                        populations,
                        pump_rates,
                        &annihilation.laws,
                        annihilation.limits,
                        cancelled,
                    )
                })?;
                (step.linear, Some(step.annihilation_heat_j))
            }
            None => (
                preserve_checkpoint(checkpoint, |cancelled| {
                    self.prepared
                        .plan
                        .advance_with_cancel(populations, pump_rates, cancelled)
                })?,
                None,
            ),
        };
        let mut cells = Vec::with_capacity(models.len());
        let mut recycled_excitation_j = Vec::with_capacity(models.len());
        let mut receiver_partition = vec![CheckedSum::default(); models.len()];
        for index in 0..models.len() {
            checkpoint()?;
            cells.push(CoupledStep {
                stored_energy_j: sum(native.stored_energy_j[index])?,
                state_stored_energy_j: native.stored_energy_j[index],
                emitted_energy_j: [0.0; M],
                nonradiative_heat_j: sum(native.nonradiative_heat_j[index])?,
                relaxation_heat_j: sum(native.relaxation_heat_j[index])?,
                crossing_heat_j: native.crossing_heat_j[index],
                reverse_crossing_heat_absorbed_j: native.reverse_crossing_heat_absorbed_j[index],
                annihilation_heat_j: annihilation_heat_j.as_ref().map_or(0.0, |heat| heat[index]),
            });
            let recycled = sum(native.recycled_energy_j[index])?;
            recycled_excitation_j.push(recycled);
            receiver_partition[index].add(-recycled)?;
        }
        let channels = self.prepared.plan.channels();
        let edges = self.prepared.plan.absorption_edges();
        let mut absorbed_by_channel = vec![CheckedSum::default(); channels.len()];
        let mut heat = BTreeMap::<CellId, CheckedSum>::new();
        let mut roundoff = CheckedSum::default();
        for (index, edge) in edges.iter().enumerate() {
            checkpoint()?;
            let absorbed = native.absorbed_energy_j[index];
            let secondary = native.absorption_heat_j[index];
            absorbed_by_channel[edge.channel].add(absorbed)?;
            if let Some(receiver) = edge.receiver_cell {
                receiver_partition[receiver].add(absorbed)?;
                receiver_partition[receiver].add(-secondary)?;
            } else {
                let mut passive_balance = CheckedSum::default();
                passive_balance.add(absorbed)?;
                passive_balance.add(-secondary)?;
                roundoff.add(passive_balance.value()?.abs())?;
            }
            if secondary != 0.0 {
                heat.entry(self.topology.absorption_receivers[index])
                    .or_default()
                    .add(secondary)?;
            }
        }
        let mut escaped = [CheckedSum::default(); M];
        let mut unresolved = [CheckedSum::default(); M];
        for balance in receiver_partition {
            roundoff.add(balance.value()?.abs())?;
        }
        for (index, channel) in channels.iter().enumerate() {
            checkpoint()?;
            let lane = self.topology.channel_bins[index];
            cells[channel.source_cell].emitted_energy_j[lane] = finite(
                cells[channel.source_cell].emitted_energy_j[lane] + native.emitted_energy_j[index],
            )?;
            escaped[lane].add(native.escaped_energy_j[index])?;
            unresolved[lane].add(native.unresolved_energy_j[index])?;
            let mut balance = CheckedSum::default();
            for energy in [
                native.emitted_energy_j[index],
                -native.escaped_energy_j[index],
                -native.unresolved_energy_j[index],
                -absorbed_by_channel[index].value()?,
            ] {
                balance.add(energy)?;
            }
            roundoff.add(balance.value()?.abs())?;
        }
        let mut escaped_j = [0.0; M];
        let mut unresolved_j = [0.0; M];
        for lane in 0..M {
            checkpoint()?;
            escaped_j[lane] = escaped[lane].value()?;
            unresolved_j[lane] = unresolved[lane].value()?;
        }
        let mut secondary_heat_j = BTreeMap::new();
        for (receiver, energy) in heat {
            checkpoint()?;
            secondary_heat_j.insert(receiver, energy.value()?);
        }
        Ok(CoupledRecyclingStep {
            cells,
            recycled_excitation_j,
            secondary_heat_j,
            escaped_j,
            unresolved_j,
            absolute_roundoff_j: roundoff.value()?,
        })
    }
}

/// Bind the original report, sampling diagnostics and explicit receiver
/// branching without relying on a synthetic source label or row normalization.
fn snapshot_identity<const M: usize>(
    environment: &AbsorbingLuminescenceEnvironment<M>,
    authored: &BTreeMap<usize, [[f64; 2]; M]>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<ContentHash, OpticalHeatingError> {
    let operator = environment.operator();
    let mut hash = DomainHasher::new("org.frankensim.render.phosphorescent-photon-transport.v1");
    hash.update(operator.medium_identity.as_bytes());
    hash.update(operator.absorption_identity.as_bytes());
    for integer in [
        M as u64,
        operator.grid_cell_count as u64,
        operator.material_cells.len() as u64,
        operator.rows.len() as u64,
        operator.traced_segments,
        u64::from(operator.config.histories_per_emitter_lane),
        u64::from(operator.config.max_collisions_per_history),
        operator.config.max_total_segments,
        operator.config.max_material_cells as u64,
        operator.config.max_emitter_cells as u64,
        operator.config.max_absorption_links as u64,
        operator.config.seed,
        operator.config.source_id,
    ] {
        hash.update(&integer.to_le_bytes());
    }
    for value in [
        operator.wavelength_bounds_nm.0,
        operator.wavelength_bounds_nm.1,
        environment.tolerance().fraction_roundoff,
        environment.tolerance().max_fraction_standard_error,
        environment.tolerance().max_unresolved_energy_j,
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    for wavelength in operator.wavelengths_nm {
        hash.update(&wavelength.to_bits().to_le_bytes());
    }
    for &cell in &operator.material_cells {
        checkpoint()?;
        hash.update(&(cell as u64).to_le_bytes());
    }
    hash.update(&(operator.material_cells.len() as u64).to_le_bytes());
    for (grid_cell, thermal_cell) in environment.bindings() {
        checkpoint()?;
        hash.update(&(grid_cell as u64).to_le_bytes());
        hash.update(&thermal_cell.0.to_le_bytes());
    }
    for row in &operator.rows {
        checkpoint()?;
        hash.update(&(row.emitter_cell as u64).to_le_bytes());
        hash.update(&(row.absorption.len() as u64).to_le_bytes());
        for estimate in [&row.absorbed, &row.escape, &row.unresolved] {
            for value in estimate.mean.iter().chain(&estimate.standard_error) {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        for &count in &row.precision_limited_histories {
            hash.update(&count.to_le_bytes());
        }
        for link in &row.absorption {
            checkpoint()?;
            hash.update(&(link.cell_index as u64).to_le_bytes());
            for value in link
                .fraction
                .mean
                .iter()
                .chain(&link.fraction.standard_error)
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    hash.update(&(authored.len() as u64).to_le_bytes());
    for (&cell, probabilities) in authored {
        checkpoint()?;
        hash.update(&(cell as u64).to_le_bytes());
        for branches in probabilities {
            for probability in branches {
                hash.update(&probability.to_bits().to_le_bytes());
            }
        }
    }
    Ok(hash.finalize())
}
