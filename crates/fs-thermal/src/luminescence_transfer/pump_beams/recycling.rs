//! Quasistatic emitted-photon transport with finite ground-state absorption.
//!
//! External pumps and emitted rays visit the SAME original emitter inventories.
//! Every optical stage solves the mean ground residence and its own mean
//! radiative source together. A source event creates one conditional spectral
//! photon, with explicit source Stokes loss. Reabsorption uses the actual
//! receiver ground population; it does not inject excitations into a full cell.
//!
//! Paths describe authored ballistic transport, not scattering or coherent
//! interference. Photon flight time is neglected relative to the kinetic clock.
//! Spectrum/geometry quadrature error belongs to the path producer. Temporal
//! refinement and bounded accelerated source closure are numerical estimates,
//! not certified error enclosures. No coefficient or active volume is inferred.

use super::super::saturation::annihilation::{
    self as pairs, AnnihilationLimits, PairHistory, TripletAnnihilation,
};
use super::evolution::SourceAttribution;
use super::*;

/// Disposition of photons reaching the end of an authored emitted path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmissionTerminal {
    Escape,
    /// A declared unresolved transport remainder, never converted to heat.
    Unresolved,
}

/// One directional/spectral quadrature contribution from an original state.
/// Photon fractions include the directional weight and sum to one over all
/// paths and lanes belonging to each emitting source state. No normalization
/// is performed. Empty paths describe directly outgoing photons.
#[derive(Clone, Debug, PartialEq)]
pub struct SaturableEmissionPath<const B: usize> {
    pub source_cell: usize,
    pub source_state: usize,
    pub photon_fractions: [f64; B],
    pub path: SaturablePumpPath<B>,
    pub terminal: EmissionTerminal,
}

/// The original optical/kinetic controls plus bounded emitted-source closure.
/// Nonlinear temporal controls apply even when every pair law is absent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturableRecyclingLimits {
    pub pump: PumpBeamsLimits,
    pub annihilation: AnnihilationLimits,
    /// Per physical stage, including safeguarded trial source evaluations.
    /// All stages additionally share `pump.column.maximum_cell_evaluations`.
    pub maximum_feedback_iterations: usize,
    /// Maximum retained inverse-secant corrections; storage and application
    /// scale with this limit times the number of radiating source states.
    /// Zero uses safeguarded physical source iteration without acceleration.
    pub maximum_secant_history: usize,
}

impl Default for SaturableRecyclingLimits {
    fn default() -> Self {
        Self {
            pump: PumpBeamsLimits::default(),
            annihilation: AnnihilationLimits::default(),
            maximum_feedback_iterations: 96,
            maximum_secant_history: 8,
        }
    }
}

/// Immutable shared geometry and original spectra. Counts and accepted
/// histories are caller-owned. Source closure uses a bounded rank-one history,
/// without imposing a global reaction-node limit on the radiating states.
/// Native per-cell kinetic limits and aggregate work limits remain in force.
#[derive(Clone, Debug)]
pub struct SaturablePhotonRecycling<const B: usize> {
    transport: SaturablePumpBeams<B>,
    external_paths: usize,
    emission: Vec<SaturableEmissionPath<B>>,
    fraction_roundoff: f64,
}

/// Accepted candidate in original cell, external-path and emitted-path order.
/// Cell pump channels contain ALL absorbed photons. Their explicit attribution
/// below separates true external supply from recycled internal traffic.
#[derive(Clone, Debug, PartialEq)]
pub struct SaturableRecyclingStep<const B: usize> {
    pub cells: Vec<SaturableSingletTripletStep>,
    pub passive_heat_j: Vec<f64>,
    pub external_beams: Vec<PumpBeamStep<B>>,
    pub emission_paths: Vec<PumpBeamStep<B>>,
    /// Actual local-field supply, separate from every transported photon path.
    /// Original cell order when active; empty when all local hazards are zero.
    pub local_pumps: Vec<LocalPumpStep>,
    pub external_absorbed_pump_energy_j: Vec<[f64; 2]>,
    pub recycled_absorbed_pump_energy_j: Vec<[f64; 2]>,
    pub external_pump_events: Vec<[f64; 2]>,
    pub recycled_pump_events: Vec<[f64; 2]>,
    pub external_pump_relaxation_heat_j: Vec<[f64; 2]>,
    pub recycled_pump_relaxation_heat_j: Vec<[f64; 2]>,
    pub external_passive_heat_j: Vec<f64>,
    pub recycled_passive_heat_j: Vec<f64>,
    /// Per source cell/state, before any receiving path absorption.
    pub source_relaxation_heat_j: Vec<[f64; 2]>,
    pub emitted_energy_j: Vec<[f64; 2]>,
    pub escaped_energy_j: Vec<[f64; 2]>,
    pub unresolved_energy_j: Vec<[f64; 2]>,
    /// Initial excitation + external incident + actual local absorption +
    /// uphill bath energy - final excitation - external transmitted - escaped -
    /// unresolved - physical heat.
    /// Internal photon traffic is excluded from both residual and error scale.
    pub residual_j: f64,
    pub substeps: usize,
    pub cell_evaluations: usize,
    pub coupling_iterations: usize,
    pub feedback_iterations: usize,
    pub max_energy_discrepancy_j: f64,
}

struct History<const B: usize> {
    beams: PumpBeamsStep<B>,
    external_absorbed: Vec<[f64; 2]>,
    recycled_absorbed: Vec<[f64; 2]>,
    external_events: Vec<[f64; 2]>,
    recycled_events: Vec<[f64; 2]>,
    external_relaxation: Vec<[f64; 2]>,
    recycled_relaxation: Vec<[f64; 2]>,
    external_passive: Vec<f64>,
    recycled_passive: Vec<f64>,
    /// Composed source-closure bound this history carries, in joules: for each
    /// accepted feedback stage, the allowance it closed to plus that stage's
    /// rounding and fraction terms, summed over the corrector stages of a
    /// midpoint history. Zero for a stationary or zero-duration history.
    /// admit_composed_source() states the bound; finish() adds it to the
    /// single-step allowance of the total balance, which holds the signed
    /// source mismatch.
    source_allowance_j: f64,
}

impl<const B: usize> SaturablePhotonRecycling<B> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        models: &[SingletTripletDecay],
        external_paths: &[SaturablePumpPath<B>],
        emission_paths: &[SaturableEmissionPath<B>],
        fraction_roundoff: f64,
        maximum_cells: usize,
        maximum_paths: usize,
        maximum_segments: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        poll(&mut cancelled)?;
        if models.is_empty() || !nonnegative(fraction_roundoff) || fraction_roundoff >= 1.0 {
            return Err(SingletTripletError::Invalid(
                "finite recycling cells or fraction allowance",
            ));
        }
        // Admit the complete original geometry before allocating spectra or
        // cloning any path. The common transport repeats these checks for its
        // own owner; a rejected wrapper must not allocate that owner first.
        if B == 0 || B > 64 || maximum_cells == 0 || maximum_paths == 0 || maximum_segments == 0 {
            return Err(SingletTripletError::Invalid(
                "pump beam lanes or admission limits",
            ));
        }
        cap("pump beam cells", models.len(), maximum_cells)?;
        let path_count = external_paths
            .len()
            .checked_add(emission_paths.len())
            .ok_or(SingletTripletError::Unrepresentable(
                "finite recycling path count",
            ))?;
        cap("pump beams", path_count, maximum_paths)?;
        let mut segment_count = 0usize;
        for path in external_paths
            .iter()
            .chain(emission_paths.iter().map(|emitted| &emitted.path))
        {
            poll(&mut cancelled)?;
            segment_count = segment_count.checked_add(path.segments.len()).ok_or(
                SingletTripletError::Unrepresentable("pump beam segment count"),
            )?;
            cap("pump beam segments", segment_count, maximum_segments)?;
        }
        segment_count
            .checked_mul(B)
            .ok_or(SingletTripletError::Unrepresentable(
                "pump beam channel count",
            ))?;
        admit_storage::<[Sum; 2]>(models.len(), "finite recycling spectral weights")?;
        admit_storage::<SingletTripletDecay>(models.len(), "finite recycling models")?;
        admit_storage::<SaturablePumpPath<B>>(path_count, "finite recycling paths")?;
        admit_storage::<SaturableEmissionPath<B>>(
            emission_paths.len(),
            "finite recycling emitted paths",
        )?;
        admit_storage::<SaturablePumpSegment<B>>(segment_count, "finite recycling path segments")?;
        let mut weights = vec![[Sum::default(); 2]; models.len()];
        let mut paths = Vec::with_capacity(path_count);
        paths.extend_from_slice(external_paths);
        for emitted in emission_paths {
            poll(&mut cancelled)?;
            if emitted.source_cell >= models.len() || emitted.source_state >= 2 {
                return Err(SingletTripletError::Invalid(
                    "finite recycling source cell/state",
                ));
            }
            let energy = models[emitted.source_cell].excitation_energies_j()[emitted.source_state];
            for lane in 0..B {
                let weight = emitted.photon_fractions[lane];
                if !nonnegative(weight)
                    || weight > 1.0
                    || (weight > 0.0 && emitted.path.photon_energies_j[lane] > energy)
                {
                    return Err(SingletTripletError::Invalid(
                        "finite recycling photon weight or source energy",
                    ));
                }
                weights[emitted.source_cell][emitted.source_state].add(weight)?;
            }
            paths.push(emitted.path.clone());
        }
        for (index, model) in models.iter().enumerate() {
            poll(&mut cancelled)?;
            for state in 0..2 {
                let sum = weights[index][state].value()?;
                if (model.radiative_rates_s()[state] > 0.0 || sum > 0.0)
                    && (sum - 1.0).abs() > fraction_roundoff
                {
                    return Err(SingletTripletError::Invalid(
                        "finite recycling conditional spectra must sum to one",
                    ));
                }
            }
        }
        let transport = SaturablePumpBeams::try_new(
            models,
            &paths,
            maximum_cells,
            maximum_paths,
            maximum_segments,
            &mut cancelled,
        )?;
        Ok(Self {
            transport,
            external_paths: external_paths.len(),
            emission: emission_paths.to_vec(),
            fraction_roundoff,
        })
    }

    pub fn models(&self) -> &[SingletTripletDecay] {
        self.transport.models()
    }
    pub fn external_paths(&self) -> &[SaturablePumpPath<B>] {
        &self.transport.paths[..self.external_paths]
    }
    pub fn emission_paths(&self) -> &[SaturableEmissionPath<B>] {
        &self.emission
    }
    pub const fn fraction_roundoff(&self) -> f64 {
        self.fraction_roundoff
    }

    /// Advance all original populations once under constant external beams.
    /// Emitted photon rates are solved from this SAME stage's radiative events.
    /// Laws follow original unique cell order; absent laws retain ordinary
    /// finite kinetics. Every refusal leaves caller state untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        populations: &[[f64; 3]],
        incident_external_s: &[[f64; B]],
        duration_s: f64,
        laws: &[Option<TripletAnnihilation>],
        limits: SaturableRecyclingLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SaturableRecyclingStep<B>, SingletTripletError> {
        self.advance_with_local_pumps(
            populations,
            incident_external_s,
            &[],
            duration_s,
            laws,
            limits,
            cancelled,
        )
    }

    /// Advance local photon hazards, external beams, radiative return and pair
    /// encounters on one shared finite inventory. Local entries are in original
    /// cell order. Empty/all-dark local supply preserves the ordinary path;
    /// nonzero local absorption is reported as separate true external input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_local_pumps(
        &self,
        populations: &[[f64; 3]],
        incident_external_s: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        duration_s: f64,
        laws: &[Option<TripletAnnihilation>],
        limits: SaturableRecyclingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SaturableRecyclingStep<B>, SingletTripletError> {
        let local_pumps = if local::validate(self.models(), local_pumps, &mut cancelled)? {
            local_pumps
        } else {
            &[]
        };
        pairs::validate(
            self.models(),
            laws,
            populations,
            duration_s,
            limits.annihilation,
            &mut cancelled,
        )?;
        if incident_external_s.len() != self.external_paths
            || limits.maximum_feedback_iterations == 0
        {
            return Err(SingletTripletError::Invalid(
                "finite recycling external sources or feedback limit",
            ));
        }
        let mut incident = incident_external_s.to_vec();
        incident.extend((0..self.emission.len()).map(|_| [0.0; B]));
        self.transport.validate_advance(
            populations,
            &incident,
            duration_s,
            limits.pump,
            &mut cancelled,
        )?;
        let column = limits.pump.column;
        let mut cells = 0usize;
        let mut coupling = 0usize;
        let mut feedback = 0usize;
        // Local mutable models; immutable geometry is copied only once per call.
        let mut frozen = self.transport.clone();
        if duration_s == 0.0 {
            let history = self
                .frozen_history(
                    &mut frozen,
                    populations,
                    &incident,
                    local_pumps,
                    duration_s,
                    limits,
                    &mut cells,
                    &mut coupling,
                    &mut cancelled,
                )?
                .ok_or(SingletTripletError::Invalid("zero-duration photon closure"))?;
            return self.finish(
                populations,
                history,
                1,
                cells,
                coupling,
                feedback,
                0.0,
                column.kinetic.relative_tolerance,
                &mut cancelled,
            );
        }
        let stationary_rates = pairs::frozen_rates(populations, laws)?;
        for (index, model) in frozen.models.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            *model =
                pairs::frozen_model(self.models()[index], laws[index], stationary_rates[index])?;
        }
        feedback += 1;
        if let Some(mut history) = self.stationary_history(
            &frozen,
            populations,
            incident_external_s,
            local_pumps,
            duration_s,
            limits,
            &mut cells,
            &mut cancelled,
        )? {
            for (index, cell) in history.beams.cells.iter_mut().enumerate() {
                pairs::restore(
                    self.models()[index],
                    frozen.models[index],
                    populations[index],
                    laws[index],
                    stationary_rates[index],
                    cell,
                    column.kinetic.relative_tolerance,
                )?;
            }
            return self.finish(
                populations,
                history,
                1,
                cells,
                coupling,
                feedback,
                0.0,
                column.kinetic.relative_tolerance,
                &mut cancelled,
            );
        }
        let mut temporal = limits.annihilation;
        temporal.maximum_substeps = temporal.maximum_substeps.min(column.maximum_substeps);
        let fastest = self.fastest_rate(
            populations,
            incident_external_s,
            local_pumps,
            &mut cancelled,
        )?;
        let result = pairs::midpoint(
            populations,
            laws,
            duration_s,
            fastest,
            temporal,
            Some((
                column.relative_tolerance,
                column.absolute_energy_tolerance_j,
            )),
            |initial, rates, duration, cancelled| {
                for (index, model) in frozen.models.iter_mut().enumerate() {
                    poll(cancelled)?;
                    *model = pairs::frozen_model(self.models()[index], laws[index], rates[index])?;
                }
                let Some(mut history) = self.feedback(
                    &mut frozen,
                    initial,
                    incident_external_s,
                    local_pumps,
                    duration,
                    limits,
                    &mut cells,
                    &mut coupling,
                    &mut feedback,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                for (index, cell) in history.beams.cells.iter_mut().enumerate() {
                    pairs::restore(
                        self.models()[index],
                        frozen.models[index],
                        initial[index],
                        laws[index],
                        rates[index],
                        cell,
                        column.kinetic.relative_tolerance,
                    )?;
                }
                Ok(Some(history))
            },
            |history, cancelled| {
                self.audit_history(
                    populations,
                    history,
                    column.kinetic.relative_tolerance,
                    cancelled,
                )
            },
            |history| self.components(history, duration_s),
            &mut cancelled,
        )?;
        self.finish(
            populations,
            result.step,
            result.substeps,
            cells,
            coupling,
            feedback,
            result.discrepancy_j,
            column.kinetic.relative_tolerance,
            &mut cancelled,
        )
    }

    fn fastest_rate(
        &self,
        initial: &[[f64; 3]],
        external: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, SingletTripletError> {
        let mut upper = external.to_vec();
        for emitted in &self.emission {
            let total = inventory(initial[emitted.source_cell])?;
            let rate = positive_product(
                total,
                self.models()[emitted.source_cell].radiative_rates_s()[emitted.source_state],
                "upper radiative rate",
            )?;
            let mut lanes = [0.0; B];
            for lane in 0..B {
                lanes[lane] =
                    positive_product(rate, emitted.photon_fractions[lane], "upper emitted lane")?;
            }
            upper.push(lanes);
        }
        let mut hazards = vec![Sum::default(); initial.len()];
        for (index, pump) in local_pumps.iter().enumerate() {
            for rate in pump.rates_s {
                hazards[index].add(rate)?;
            }
        }
        for (path, rates) in self.transport.paths.iter().zip(upper) {
            for segment in &path.segments {
                poll(cancelled)?;
                for lane in 0..B {
                    for target in 0..2 {
                        hazards[segment.cell_index].add(positive_product(
                            rates[lane],
                            segment.optical_depth_per_ground[lane][target],
                            "upper recycled hazard",
                        )?)?;
                    }
                }
            }
        }
        let mut fastest = pairs::fastest_rate(self.models());
        for hazard in hazards {
            fastest = fastest.max(hazard.value()?);
        }
        Ok(fastest)
    }
}

fn inventory(populations: [f64; 3]) -> Result<f64, SingletTripletError> {
    let mut total = Sum::default();
    for value in populations {
        total.add(value)?;
    }
    total.value()
}

fn admit_storage<T>(length: usize, quantity: &'static str) -> Result<(), SingletTripletError> {
    let bytes = length
        .checked_mul(std::mem::size_of::<T>())
        .ok_or(SingletTripletError::Unrepresentable(quantity))?;
    if bytes > isize::MAX as usize {
        return Err(SingletTripletError::Unrepresentable(quantity));
    }
    Ok(())
}

#[path = "recycling/accounts.rs"]
mod accounts;
#[path = "recycling/feedback.rs"]
mod feedback;
#[cfg(test)]
#[path = "recycling/local_oracles.rs"]
mod local_oracles;
#[path = "recycling/secant.rs"]
mod secant;
#[path = "recycling/stationary.rs"]
mod stationary;
#[cfg(test)]
#[path = "recycling/tests.rs"]
mod tests;
