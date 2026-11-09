//! Finite emitter inventories on the original thermal and observation clock.
//!
//! A prescribed local photon field supplies a per-ground-emitter hazard
//! `W = cross_section * photon_flux`. Actual absorption is `W * N_ground`,
//! integrated by the native three-state kinetic owner. An optically thin or
//! externally maintained pump field must be justified by the caller for the
//! default prescribed-field mode. The optional [`beam`] adapter instead solves
//! pump depletion against the finite ground inventory and separately accounts
//! incident/transmitted photons. Only actual molecular pump events and passive
//! host absorption supply external energy. Emitted light may use the existing
//! passive absorption/scattering environment; no active photon recycling is inferred.
//!
//! Original scalar or material-bound singlet/triplet rates refresh at each
//! thermal interval. Rates and local pump flux are frozen over that interval;
//! thermal-step refinement assesses feedback error. Intrinsic lifetimes do
//! not force lifetime-sized thermal steps. Stimulated emission, excited-state
//! absorption and irreversible photobleaching remain outside this model.
//! Explicit triplet-pair laws include finite-inventory recovery and delayed
//! fluorescence; reversible ground-state bleaching is part of the beam solve.

/// Authored or original-material triplet-pair kinetics on the shared inventory.
pub mod annihilation;
/// Self-consistent finite-inventory depletion of a directed spectral pump.
pub mod beam;

use super::phosphorescence::preserve_checkpoint;
use super::{
    AbsorbingLuminescenceEnvironment, Accounts, CheckedSum, LuminescenceEnergy, LuminescenceSource,
    OpticalEnergyBudget, OpticalHeatingError, SingletTripletSource, finite, invalid, nonnegative,
    partition,
};
use crate::volumes::transport::AbsorbingGrid;
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use fs_thermal::luminescence::transfer::SingletTripletDecay;
use fs_thermal::luminescence::transfer::saturation::SaturableSingletTripletPlan;
use std::{collections::BTreeMap, sync::Arc};

const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;

/// Two independently authored optical pump channels, ordered `[G -> S, G -> T]`.
/// Flux is photons/(m² s), cross section is m² per ground-state emitter, and
/// wavelength is vacuum nm. A zero cross section or flux disables that channel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedPhotonPump {
    source_identity: ContentHash,
    wavelengths_nm: [f64; 2],
    cross_sections_m2: [f64; 2],
    photon_flux_m2_s: [f64; 2],
    rates_s: [f64; 2],
    energies_j: [f64; 2],
}

impl PrescribedPhotonPump {
    /// Bind explicit cross sections and the local undepleted photon field.
    /// A source identity is attribution, not a measured-material certificate.
    pub fn try_new(
        source_identity: ContentHash,
        wavelengths_nm: [f64; 2],
        cross_sections_m2: [f64; 2],
        photon_flux_m2_s: [f64; 2],
    ) -> Result<Self, OpticalHeatingError> {
        if source_identity == ContentHash([0; 32])
            || !wavelengths_nm.iter().all(|x| x.is_finite() && *x > 0.0)
            || !cross_sections_m2.into_iter().all(nonnegative)
            || !photon_flux_m2_s.into_iter().all(nonnegative)
        {
            return Err(invalid("invalid prescribed photon pump"));
        }
        let mut rates_s = [0.0; 2];
        let mut energies_j = [0.0; 2];
        for i in 0..2 {
            rates_s[i] = finite(cross_sections_m2[i] * photon_flux_m2_s[i])?;
            energies_j[i] = finite(HC_J_NM / wavelengths_nm[i])?;
            if energies_j[i] <= 0.0
                || (cross_sections_m2[i] > 0.0 && photon_flux_m2_s[i] > 0.0 && rates_s[i] == 0.0)
            {
                return Err(invalid("prescribed pump exceeds numerical representation"));
            }
        }
        Ok(Self {
            source_identity,
            wavelengths_nm,
            cross_sections_m2,
            photon_flux_m2_s,
            rates_s,
            energies_j,
        })
    }

    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    pub const fn wavelengths_nm(&self) -> [f64; 2] {
        self.wavelengths_nm
    }
    pub const fn cross_sections_m2(&self) -> [f64; 2] {
        self.cross_sections_m2
    }
    pub const fn photon_flux_m2_s(&self) -> [f64; 2] {
        self.photon_flux_m2_s
    }
    /// Hazards per available ground-state emitter, inverse seconds.
    pub const fn rates_s(&self) -> [f64; 2] {
        self.rates_s
    }
}

/// One explicit inventory mapped to its original optical and thermal cells.
/// Populations are mean counts `[ground, singlet, triplet]`, not concentrations.
/// Initial excitation is included in the energy account, never treated as heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturableEmitterCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub initial_populations: [f64; 3],
    pub pump: PrescribedPhotonPump,
}

#[derive(Clone, Debug)]
struct Emitter {
    binding: SaturableEmitterCell,
    populations: [f64; 3],
    prepared: Option<Arc<SaturableSingletTripletPlan>>,
}

#[derive(Clone, Debug)]
struct Pending<const M: usize> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    populations: Vec<[f64; 3]>,
    mean_emission_w: Vec<(usize, [f64; M])>,
    accounts: Accounts<2, M>,
    beam_accounts: Option<beam::BeamAccounts>,
}

struct PreparedState {
    beam: Option<Arc<beam::PreparedPumps>>,
    annihilation: Option<Arc<annihilation::PreparedAnnihilation>>,
}

/// Saturation, triplet shelving, heat and observable afterglow on one timeline.
/// Clones share only immutable prepared propagators and optical environment;
/// populations and accepted accounts remain transactional. Preparation and
/// stepping are bounded by the supplied cell and native kinetic limits.
#[derive(Clone, Debug)]
pub struct SaturableLuminescentHeating<const M: usize> {
    singlet: LuminescenceSource<M>,
    triplet: SingletTripletSource<M>,
    emitters: Vec<Emitter>,
    environment: Option<AbsorbingLuminescenceEnvironment<M>>,
    pump_window_s: Option<(f64, f64)>,
    medium_identity: ContentHash,
    beam: Option<beam::DepletingBeam>,
    annihilation: Option<annihilation::AnnihilationFeedback>,
    source_identity: ContentHash,
    time_s: f64,
    energy_tolerance_j: f64,
    accounts: Accounts<2, M>,
    pump_account: OpticalEnergyBudget<2>,
    energy: LuminescenceEnergy<M>,
    emission_w: Vec<(usize, [f64; M])>,
    mean_emission_w: Vec<(usize, [f64; M])>,
    pending: Option<Pending<M>>,
}

impl<const M: usize> SaturableLuminescentHeating<M> {
    /// Bind finite populations, original spectra and a prescribed pump field.
    /// The complete medium must be transparent at emission nodes unless an
    /// explicit passive transport environment supplies absorption and escape.
    /// The medium describes emission transport; its pump-band Beer coefficient
    /// is not used as a substitute for the supplied molecular cross section.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        medium: &AbsorbingGrid<'_>,
        singlet: LuminescenceSource<M>,
        triplet: SingletTripletSource<M>,
        bindings: &[SaturableEmitterCell],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        energy_tolerance_j: f64,
        environment: Option<AbsorbingLuminescenceEnvironment<M>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if !nonnegative(energy_tolerance_j)
            || medium.cell_count() > max_grid_cells
            || bindings.len() > max_grid_cells
            || singlet.fluorescence().wavelengths_nm()
                != triplet.triplet().fluorescence().wavelengths_nm()
        {
            return Err(invalid(
                "saturable source spectrum, grid budget or energy allowance",
            ));
        }
        let wavelengths = singlet.fluorescence().wavelengths_nm();
        if let Some(environment) = &environment {
            if environment.operator().absorption_identity != medium.content_identity()
                || environment.operator().wavelengths_nm != *wavelengths
            {
                return Err(invalid(
                    "saturable emission transport differs from original medium",
                ));
            }
            environment.check_state(thermal, &mut checkpoint)?;
        } else {
            for cell in 0..medium.cell_count() {
                checkpoint()?;
                let coefficients = medium
                    .coefficients_at_cell(cell, wavelengths)
                    .map_err(|_| invalid("unsupported saturable emission wavelength"))?;
                if medium.cell_density_kg_m3(cell).is_some_and(|rho| rho > 0.0)
                    && coefficients
                        .into_iter()
                        .any(|coefficient| coefficient != 0.0)
                {
                    return Err(invalid(
                        "absorbing emission requires explicit passive transport",
                    ));
                }
            }
        }
        let mut ordered = BTreeMap::new();
        for binding in bindings {
            checkpoint()?;
            if binding.grid_cell >= medium.cell_count()
                || !binding.initial_populations.into_iter().all(nonnegative)
            {
                return Err(invalid("invalid finite emitter inventory or optical cell"));
            }
            if ordered.insert(binding.grid_cell, *binding).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.saturable-luminescence.v1");
        hash.update(medium.content_identity().as_bytes());
        hash.update(triplet.rate_source_identity().as_bytes());
        for source in [&singlet, triplet.triplet()] {
            hash_source(&mut hash, source);
        }
        hash.update(&(ordered.len() as u64).to_le_bytes());
        let mut initial = CheckedSum::default();
        let mut emitters = Vec::with_capacity(ordered.len());
        let mut emission_w = Vec::with_capacity(ordered.len());
        let mut cached_model = None;
        for binding in ordered.into_values() {
            checkpoint()?;
            let temperature = thermal.temperature_k(binding.thermal_cell)?;
            let model = resolve_model(&singlet, &triplet, temperature, &mut cached_model)?;
            for (i, spectrum) in [singlet.fluorescence(), triplet.triplet().fluorescence()]
                .into_iter()
                .enumerate()
            {
                let (lo, hi) = spectrum.excitation_bounds_nm();
                if binding.pump.rates_s[i] > 0.0
                    && !(lo..=hi).contains(&binding.pump.wavelengths_nm[i])
                {
                    return Err(invalid(
                        "prescribed pump exceeds original excitation support",
                    ));
                }
            }
            let plan = preserve_checkpoint(&mut checkpoint, |cancelled| {
                SaturableSingletTripletPlan::prepare(
                    model,
                    binding.pump.rates_s,
                    binding.pump.energies_j,
                    0.0,
                    triplet.evolution_limits(),
                    cancelled,
                )
            })?;
            let state = preserve_checkpoint(&mut checkpoint, |cancelled| {
                plan.advance_with_cancel(binding.initial_populations, cancelled)
            })?;
            for value in state.stored_energy_j {
                initial.add(value)?;
            }
            if let Some(environment) = &environment {
                environment.check_emitter(binding.grid_cell, binding.thermal_cell)?;
            }
            hash.update(&(binding.grid_cell as u64).to_le_bytes());
            hash.update(&binding.thermal_cell.0.to_le_bytes());
            hash.update(binding.pump.source_identity.as_bytes());
            for value in binding
                .initial_populations
                .into_iter()
                .chain(binding.pump.wavelengths_nm)
                .chain(binding.pump.cross_sections_m2)
                .chain(binding.pump.photon_flux_m2_s)
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
            emission_w.push((
                binding.grid_cell,
                endpoint_emission(
                    &singlet,
                    triplet.triplet(),
                    model,
                    binding.initial_populations,
                )?,
            ));
            emitters.push(Emitter {
                populations: binding.initial_populations,
                binding,
                prepared: Some(Arc::new(plan)),
            });
        }
        let accounts = Accounts {
            initial_excitation_j: initial.value()?,
            ..Accounts::default()
        };
        let (pump_account, energy) = accounts.report(initial.value()?)?;
        checkpoint()?;
        Ok(Self {
            singlet,
            triplet,
            emitters,
            environment,
            pump_window_s: None,
            medium_identity: medium.content_identity(),
            beam: None,
            annihilation: None,
            source_identity: hash.finalize(),
            time_s: thermal.time_s(),
            energy_tolerance_j,
            accounts,
            pump_account,
            energy,
            mean_emission_w: Vec::new(),
            emission_w,
            pending: None,
        })
    }

    /// Restrict all pumps to one exact global `[start, end)` pulse.
    /// Individual beam windows, when bound, intersect this global gate.
    /// Both edges split the native thermal interval; afterglow continues in the
    /// dark. Bind before advancing. Without a window the pump stays on.
    pub fn with_pump_window(
        mut self,
        start_s: f64,
        end_s: f64,
    ) -> Result<Self, OpticalHeatingError> {
        if self.energy.duration_s != 0.0
            || self.pending.is_some()
            || self.pump_window_s.is_some()
            || !start_s.is_finite()
            || !end_s.is_finite()
            || start_s < self.time_s
            || end_s <= start_s
        {
            return Err(invalid("invalid finite-population pump window"));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.saturable-pump-window.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(&start_s.to_bits().to_le_bytes());
        hash.update(&end_s.to_bits().to_le_bytes());
        self.source_identity = hash.finalize();
        self.pump_window_s = Some((start_s, end_s));
        Ok(self)
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn energy(&self) -> &LuminescenceEnergy<M> {
        &self.energy
    }
    /// Actual external pump consumption, ordered `[G -> S, G -> T]`, joules.
    pub const fn absorbed_pump_energy_j(&self) -> [f64; 2] {
        self.pump_account.absorbed_j
    }
    /// Net heat already supplied to the original thermal cells, joules.
    pub const fn thermal_input_j(&self) -> f64 {
        self.pump_account.thermal_input_j
    }
    pub const fn emission_source_identity(&self) -> ContentHash {
        self.source_identity
    }
    /// Accepted endpoint counts `[ground, singlet, triplet]` for an original cell.
    pub fn populations(&self, grid_cell: usize) -> Option<[f64; 3]> {
        self.emitters
            .binary_search_by_key(&grid_cell, |e| e.binding.grid_cell)
            .ok()
            .map(|i| self.emitters[i].populations)
    }
    /// Endpoint bin-integrated emitted watts, before passive propagation.
    pub fn cell_emission_w(&self) -> &[(usize, [f64; M])] {
        &self.emission_w
    }
    /// Last accepted interval's bin-integrated mean watts, suitable for a fixed
    /// geometry/optics exposure. This is distinct from the endpoint afterglow.
    pub fn last_mean_cell_emission_w(&self) -> &[(usize, [f64; M])] {
        &self.mean_emission_w
    }

    fn check_state(
        &self,
        thermal: &DepositionSimulation<'_>,
        match_clock: bool,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<PreparedState, OpticalHeatingError> {
        checkpoint()?;
        if match_clock && self.time_s != thermal.time_s() {
            return Err(invalid("finite-population and thermal clocks differ"));
        }
        let mut cached_model = None;
        for emitter in &self.emitters {
            checkpoint()?;
            resolve_model(
                &self.singlet,
                &self.triplet,
                thermal.temperature_k(emitter.binding.thermal_cell)?,
                &mut cached_model,
            )?;
        }
        if let Some(environment) = &self.environment {
            environment.check_state(thermal, checkpoint)?;
        }
        let annihilation = self
            .annihilation
            .as_ref()
            .map(|pairs| pairs.prepare(self, thermal, checkpoint))
            .transpose()?;
        let beam = self
            .beam
            .as_ref()
            .map(|beam| beam.check_state(self, thermal, checkpoint))
            .transpose()?;
        Ok(PreparedState { beam, annihilation })
    }

    fn next_edge(&self, after_s: f64) -> Option<f64> {
        let global = self.pump_window_s.and_then(|(start, end)| {
            if after_s < start {
                Some(start)
            } else if after_s < end {
                Some(end)
            } else {
                None
            }
        });
        let beam = self
            .beam
            .as_ref()
            .and_then(|beams| beams.next_edge(after_s));
        match (global, beam) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (edge, None) | (None, edge) => edge,
        }
    }
}

impl<const M: usize> DepositionObserver for SaturableLuminescentHeating<M> {
    type Error = OpticalHeatingError;
    fn requires_step_interval(&self) -> bool {
        true
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.check_state(thermal, true, checkpoint)?;
        Ok(self.next_edge(self.time_s))
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_state(thermal, true, checkpoint)?;
        Ok(BTreeMap::new())
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let prepared = self.check_state(thermal, true, checkpoint)?;
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !end_s.is_finite()
            || end_s <= start_s
            || (end_s != start_s + duration_s && duration_s != end_s - start_s)
            || self.next_edge(start_s).is_some_and(|edge| edge < end_s)
        {
            return Err(invalid(
                "invalid finite-population native interval or pulse edge",
            ));
        }
        let pumped = self
            .pump_window_s
            .is_none_or(|(start, end)| start_s >= start && start_s < end);
        let mut accounts = self.accounts.clone();
        accounts.duration.add(duration_s)?;
        let mut populations = Vec::with_capacity(self.emitters.len());
        let mut mean_emission_w = Vec::with_capacity(self.emitters.len());
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut beam_candidate = self
            .beam
            .as_ref()
            .map(|beam| {
                beam.advance(
                    self,
                    thermal,
                    duration_s,
                    pumped,
                    prepared.annihilation.as_deref(),
                    checkpoint,
                )
            })
            .transpose()?;
        if let Some(candidate) = &beam_candidate {
            for &(cell, heat_j) in &candidate.passive_heat_j {
                checkpoint()?;
                sources
                    .entry(cell)
                    .or_default()
                    .add(finite(heat_j / duration_s)?)?;
                accounts.passive_pump.add(heat_j)?;
            }
            accounts.roundoff.add(candidate.absolute_roundoff_j)?;
        }
        let mut cached_model = None;
        let mut previous_plan: Option<Arc<SaturableSingletTripletPlan>> = None;
        for (emitter_index, emitter) in self.emitters.iter_mut().enumerate() {
            checkpoint()?;
            let binding = emitter.binding;
            let step = if let Some(step) = beam_candidate
                .as_mut()
                .and_then(|candidate| candidate.steps[emitter_index].take())
            {
                step
            } else {
                let model = resolve_model(
                    &self.singlet,
                    &self.triplet,
                    thermal.temperature_k(binding.thermal_cell)?,
                    &mut cached_model,
                )?;
                let rates = if pumped && self.beam.is_none() {
                    binding.pump.rates_s
                } else {
                    [0.0; 2]
                };
                let matches = |plan: &SaturableSingletTripletPlan| {
                    *plan.model() == model
                        && plan.duration_s() == duration_s
                        && plan.pump_rates_s() == rates
                        && plan.pump_energies_j() == binding.pump.energies_j
                };
                if emitter.prepared.as_ref().is_none_or(|plan| !matches(plan)) {
                    if let Some(plan) = previous_plan.as_ref().filter(|plan| matches(plan)) {
                        emitter.prepared = Some(Arc::clone(plan));
                    } else {
                        emitter.prepared =
                            Some(Arc::new(preserve_checkpoint(checkpoint, |cancelled| {
                                SaturableSingletTripletPlan::prepare(
                                    model,
                                    rates,
                                    binding.pump.energies_j,
                                    duration_s,
                                    self.triplet.evolution_limits(),
                                    cancelled,
                                )
                            })?));
                    }
                }
                let plan = emitter
                    .prepared
                    .as_ref()
                    .expect("prepared finite emitter interval");
                previous_plan = Some(Arc::clone(plan));
                if let Some(pairs) = &prepared.annihilation {
                    preserve_checkpoint(checkpoint, |cancelled| {
                        plan.advance_with_annihilation(
                            emitter.populations,
                            pairs.laws[emitter_index],
                            pairs.limits,
                            cancelled,
                        )
                    })?
                    .kinetics
                } else {
                    preserve_checkpoint(checkpoint, |cancelled| {
                        plan.advance_with_cancel(emitter.populations, cancelled)
                    })?
                }
            };
            let (singlet_light, singlet_heat) =
                partition(&self.singlet, step.radiative_decay_energy_j[0])?;
            let (triplet_light, triplet_heat) =
                partition(self.triplet.triplet(), step.radiative_decay_energy_j[1])?;
            let mut light = [0.0; M];
            let mut mean = [0.0; M];
            let mut local = CheckedSum::default();
            local.add(step.residual_j)?;
            let relaxation = finite(singlet_heat + triplet_heat)?;
            let mut heat = CheckedSum::default();
            heat.add(relaxation)?;
            heat.add(step.annihilation_heat_j)?;
            accounts.relaxation.add(relaxation)?;
            accounts.annihilation.add(step.annihilation_heat_j)?;
            for i in 0..2 {
                accounts.absorbed[i].add(step.absorbed_pump_energy_j[i])?;
                accounts.prompt.add(step.pump_relaxation_heat_j[i])?;
                accounts.dark.add(step.nonradiative_heat_j[i])?;
                heat.add(step.pump_relaxation_heat_j[i])?;
                heat.add(step.nonradiative_heat_j[i])?;
                local.add(step.radiative_decay_energy_j[i])?;
            }
            heat.add(step.crossing_heat_j)?;
            heat.add(-step.reverse_crossing_heat_absorbed_j)?;
            accounts.crossing.add(step.crossing_heat_j)?;
            accounts
                .reverse_crossing
                .add(step.reverse_crossing_heat_absorbed_j)?;
            sources
                .entry(binding.thermal_cell)
                .or_default()
                .add(finite(heat.value()? / duration_s)?)?;
            local.add(-relaxation)?;
            for bin in 0..M {
                light[bin] = finite(singlet_light[bin] + triplet_light[bin])?;
                accounts.emitted[bin].add(light[bin])?;
                local.add(-light[bin])?;
                mean[bin] = finite(light[bin] / duration_s)?;
            }
            accounts.roundoff.add(local.value()?.abs())?;
            if let Some(environment) = &self.environment {
                let routed = environment.route(
                    binding.grid_cell,
                    light,
                    duration_s,
                    &mut sources,
                    checkpoint,
                )?;
                for bin in 0..M {
                    accounts.escaped[bin].add(routed.escaped_j[bin])?;
                    accounts.unresolved[bin].add(routed.unresolved_j[bin])?;
                }
                accounts.secondary.add(routed.heat_j)?;
                accounts.roundoff.add(routed.absolute_roundoff_j)?;
            } else {
                for (bin, value) in light.into_iter().enumerate() {
                    accounts.escaped[bin].add(value)?;
                }
            }
            populations.push(step.populations);
            mean_emission_w.push((binding.grid_cell, mean));
        }
        if accounts.roundoff.value()? > self.energy_tolerance_j {
            return Err(invalid(
                "finite-population local energy defects exceed allowance",
            ));
        }
        if let Some(environment) = &self.environment {
            let mut unresolved = CheckedSum::default();
            for value in &accounts.unresolved {
                unresolved.add(value.value()?)?;
            }
            if unresolved.value()? > environment.tolerance().max_unresolved_energy_j {
                return Err(invalid(
                    "finite-population unresolved light exceeds allowance",
                ));
            }
        }
        let sources = sources
            .into_iter()
            .map(|(id, power)| Ok((id, power.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        checkpoint()?;
        self.pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            populations,
            mean_emission_w,
            accounts,
            beam_accounts: beam_candidate.map(|candidate| candidate.accounts),
        });
        Ok(sources)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_state(thermal, true, &mut || Ok(())).map(|_| ())
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let prepared = self.check_state(thermal, false, checkpoint)?;
        let pending = self
            .pending
            .as_ref()
            .ok_or_else(|| invalid("finite emitter step lacks candidate"))?;
        if pending.start_s != start_s
            || pending.duration_s != duration_s
            || pending.end_s != thermal.time_s()
            || start_s != self.time_s
        {
            return Err(invalid(
                "finite emitter candidate differs from thermal step",
            ));
        }
        let mut stored = CheckedSum::default();
        let mut emission_w = Vec::with_capacity(self.emitters.len());
        let mut cached_model = None;
        for (emitter, &populations) in self.emitters.iter().zip(&pending.populations) {
            checkpoint()?;
            let model = resolve_model(
                &self.singlet,
                &self.triplet,
                thermal.temperature_k(emitter.binding.thermal_cell)?,
                &mut cached_model,
            )?;
            for (i, energy) in model.excitation_energies_j().into_iter().enumerate() {
                stored.add(finite(populations[i + 1] * energy)?)?;
            }
            emission_w.push((
                emitter.binding.grid_cell,
                endpoint_emission(&self.singlet, self.triplet.triplet(), model, populations)?,
            ));
        }
        let (pump_account, energy) = pending.accounts.report(stored.value()?)?;
        if energy.residual_j.abs() > self.energy_tolerance_j {
            return Err(invalid(
                "finite emitter pump, heat, light and excitation do not balance",
            ));
        }
        checkpoint()?;
        let pending = self
            .pending
            .take()
            .expect("admitted finite emitter candidate");
        for (emitter, populations) in self.emitters.iter_mut().zip(pending.populations) {
            emitter.populations = populations;
        }
        self.accounts = pending.accounts;
        if let (Some(beam), Some(accounts), Some(prepared)) =
            (&mut self.beam, pending.beam_accounts, prepared.beam)
        {
            beam.accept(accounts, prepared);
        }
        if let (Some(pairs), Some(prepared)) = (&mut self.annihilation, prepared.annihilation) {
            pairs.accept(prepared);
        }
        self.mean_emission_w = pending.mean_emission_w;
        self.emission_w = emission_w;
        self.pump_account = pump_account;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

// Consecutive cells commonly share a physical temperature. Keep this cache
// local to one pass: its memory is constant, no stale thermal state survives,
// and every original cell still performs its own cancellation/state lookup.
fn resolve_model<const M: usize>(
    singlet: &LuminescenceSource<M>,
    triplet: &SingletTripletSource<M>,
    temperature_k: f64,
    cached: &mut Option<(u64, SingletTripletDecay)>,
) -> Result<SingletTripletDecay, OpticalHeatingError> {
    let key = temperature_k.to_bits();
    if let Some((previous, model)) = cached
        && *previous == key
    {
        return Ok(*model);
    }
    let model = triplet.model_at_temperature(singlet, temperature_k)?;
    *cached = Some((key, model));
    Ok(model)
}

fn endpoint_emission<const M: usize>(
    singlet: &LuminescenceSource<M>,
    triplet: &LuminescenceSource<M>,
    model: SingletTripletDecay,
    populations: [f64; 3],
) -> Result<[f64; M], OpticalHeatingError> {
    let energies = model.excitation_energies_j();
    let rates = model.radiative_rates_s();
    let (s, _) = partition(
        singlet,
        finite(finite(populations[1] * energies[0])? * rates[0])?,
    )?;
    let (t, _) = partition(
        triplet,
        finite(finite(populations[2] * energies[1])? * rates[1])?,
    )?;
    let mut result = [0.0; M];
    for i in 0..M {
        result[i] = finite(s[i] + t[i])?;
    }
    Ok(result)
}

fn hash_source<const M: usize>(hash: &mut DomainHasher, source: &LuminescenceSource<M>) {
    hash.update(source.fluorescence().source_identity().as_bytes());
    let (lo, hi) = source.fluorescence().temperature_bounds_k();
    let (a, b) = source.fluorescence().excitation_bounds_nm();
    for value in [source.excitation_wavelength_nm(), lo, hi, a, b] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    for value in source
        .fluorescence()
        .wavelengths_nm()
        .iter()
        .chain(source.fluorescence().photon_fractions())
    {
        hash.update(&value.to_bits().to_le_bytes());
    }
    match source.decay {
        Some(decay) => {
            hash.update(&[1]);
            for value in [
                decay.lifetime_s(),
                decay.quantum_yield(),
                decay.nonradiative_yield(),
            ] {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        None => hash.update(&[0]),
    }
    match source.quenching {
        Some(law) => {
            hash.update(&[1]);
            hash.update(&law.prefactor_rate_s().to_bits().to_le_bytes());
            hash.update(&law.activation_energy_j().to_bits().to_le_bytes());
        }
        None => hash.update(&[0]),
    }
}
