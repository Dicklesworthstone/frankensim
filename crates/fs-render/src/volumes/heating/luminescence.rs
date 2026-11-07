//! Finite-lifetime luminescence on the original conservative thermal clock.
//!
//! Absorption first relaxes to an explicitly authored emitting-state energy.
//! That prompt energy difference heats the absorbing cell. The remaining
//! excitation is a separate store, advanced by `fs_thermal::luminescence`.
//! Radiative decay produces the original photon-weighted spectrum; dark decay
//! and final vibrational relaxation heat the same native cell exactly once.
//!
//! The default is a linear, unsaturated emitting-state model. Explicit
//! [`phosphorescence`] adds a separate triplet, intersystem crossing, distinct
//! emission spectra and the forward/reverse crossing heat. Neither path infers
//! a trap distribution, saturation, bleaching or stimulated emission. The fast path
//! requires transparent emission-band escape. An explicitly supplied directed
//! transport operator also supports passive secondary absorption and elastic
//! multiple scattering in an open, index-one medium. Unresolved photons remain
//! separate from heat and escape. Explicit single-line resonant recycling
//! returns absorbed photons to equal-energy stores and evolves the coupled
//! network within the same interval. Supplied rates/spectra require their own
//! source support; a long lifetime alone does not establish phosphorescence.

pub mod absorption;
pub mod phosphorescence;
pub use absorption::{AbsorbingLuminescenceEnvironment, LuminescenceAbsorptionTolerance};
pub use phosphorescence::{SingletTripletLimits, SingletTripletSource};
mod recycling;
pub use fs_thermal::luminescence::recycling::RecyclingLimits;
use phosphorescence::SingletTripletFeedback;
use recycling::ResonantFeedback;

use super::feedback::OpticalFeedbackError;
use super::fluorescence::{FluorescenceSource, ratio_power, validate_fluorescent_medium};
use super::illumination::schedule::SpectralPowerSchedule;
use super::{
    CheckedSum, OpticalEnergyBudget, OpticalHeating, OpticalHeatingError, OpticalPowerBudget,
    finite, invalid, nonnegative,
};
use crate::volumes::transport::AbsorbingGrid;
use crate::volumes::transport::deposition::BeamPowerDeposition;
use crate::volumes::transport::scattering::ScatteringGrid;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use fs_thermal::luminescence::{ArrheniusQuenching, ExcitationDecay, ExcitationStep};
use std::collections::BTreeMap;

/// One sourced emitting state, with a discrete photon spectrum.
/// The source's quantum yield and `lifetime_s` specify the unquenched state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuminescenceSource<const M: usize> {
    fluorescence: FluorescenceSource<M>,
    excitation_wavelength_nm: f64,
    decay: ExcitationDecay,
    quenching: Option<ArrheniusQuenching>,
}

impl<const M: usize> LuminescenceSource<M> {
    /// `excitation_wavelength_nm` explicitly defines stored energy per excited
    /// emitter as hc/lambda. It must be at or beyond the entire pump band and
    /// no longer than the shortest positive-weight emission wavelength. No level energy is
    /// inferred from a color, temperature, material name or spectral centroid.
    pub fn try_new(
        fluorescence: FluorescenceSource<M>,
        excitation_wavelength_nm: f64,
        lifetime_s: f64,
    ) -> Result<Self, OpticalHeatingError> {
        if !excitation_wavelength_nm.is_finite()
            || excitation_wavelength_nm < fluorescence.excitation_bounds_nm().1
            || fluorescence
                .wavelengths_nm()
                .iter()
                .zip(fluorescence.photon_fractions())
                .any(|(&wavelength, &weight)| weight > 0.0 && excitation_wavelength_nm > wavelength)
        {
            return Err(invalid(
                "stored excitation must lie energetically between pump and every emission node",
            ));
        }
        let decay = ExcitationDecay::try_new(lifetime_s, fluorescence.quantum_yield())
            .map_err(|_| invalid("invalid luminescence lifetime or quantum yield"))?;
        Ok(Self {
            fluorescence,
            excitation_wavelength_nm,
            decay,
            quenching: None,
        })
    }

    /// Add an explicitly supplied activated nonradiative channel. Lifetime and
    /// quantum yield change together; the radiative decay rate stays fixed.
    /// Temperature validity remains the original fluorescence source interval.
    pub fn with_quenching(mut self, quenching: ArrheniusQuenching) -> Self {
        self.quenching = Some(quenching);
        self
    }

    /// Original photon spectrum, unquenched yield and source validity domains.
    pub const fn fluorescence(&self) -> &FluorescenceSource<M> {
        &self.fluorescence
    }
    /// Authored emitting-state energy expressed as its equivalent vacuum nm.
    pub const fn excitation_wavelength_nm(&self) -> f64 {
        self.excitation_wavelength_nm
    }

    /// Effective decay at an admitted native material temperature. Rates are
    /// frozen at the start of each accepted thermal substep; refine that step
    /// to assess temperature feedback, independently of the analytic decay.
    pub fn decay_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ExcitationDecay, OpticalHeatingError> {
        let (lo, hi) = self.fluorescence.temperature_bounds_k();
        if !temperature_k.is_finite() || !(lo..=hi).contains(&temperature_k) {
            return Err(invalid(
                "temperature exceeds the luminescence source domain",
            ));
        }
        match self.quenching {
            None => Ok(self.decay),
            Some(law) => {
                let rate = law
                    .rate_at_temperature(temperature_k)
                    .map_err(|_| invalid("unrepresentable luminescence quenching rate"))?;
                self.decay
                    .with_nonradiative_rate(rate)
                    .map_err(|_| invalid("unrepresentable quenched luminescence lifetime"))
            }
        }
    }
}

/// Accepted radiation, thermal loss and excitation storage in joules.
/// All heat fields are already included in the native external-heat
/// audit. `stored_excitation_j` is separate from the thermal enthalpy store.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuminescenceEnergy<const M: usize> {
    /// Sum of accepted native interval durations, seconds.
    pub duration_s: f64,
    /// Escaped joules in the supplied discrete emission-node order.
    pub escaped_j: [f64; M],
    /// Total radiation produced before propagation, including subsequently
    /// absorbed light and repeated recycling. Internal emission can exceed
    /// external pump input; this diagnostic is not an additional energy sink.
    pub emitted_j: [f64; M],
    /// Unassigned radiant energy from admitted transport cutoffs, joules.
    /// It is not a physical excitation store or escaped radiation.
    pub unresolved_radiation_j: [f64; M],
    /// Immediate pump-to-emitting-state relaxation heat, joules.
    pub prompt_heat_j: f64,
    /// Stored excitation released by nonradiative transitions, joules.
    pub nonradiative_heat_j: f64,
    /// Heat supplied by downhill singlet-to-triplet intersystem crossing.
    pub crossing_heat_j: f64,
    /// Heat drawn from the thermal bath by uphill triplet-to-singlet crossing.
    pub reverse_crossing_heat_absorbed_j: f64,
    /// Radiative-transition energy remaining after emitted photon energy, joules.
    pub relaxation_heat_j: f64,
    /// Emitted radiation absorbed by explicitly passive receivers, joules.
    pub secondary_absorption_heat_j: f64,
    /// Emitted radiation returned to resonant excitation stores, joules.
    /// This internal circulation is neither extra pump input nor heat.
    pub recycled_excitation_j: f64,
    /// Current total nonthermal excitation store, including both spin states
    /// when singlet/triplet coupling is explicitly bound, joules.
    pub stored_excitation_j: f64,
    /// Absorbed pump minus escape, unresolved light, all heat and excitation.
    pub residual_j: f64,
    /// Sum of absolute local energy defects over the accepted history.
    pub absolute_roundoff_j: f64,
}

#[derive(Clone, Debug)]
struct Cell<const N: usize> {
    grid_cell: usize,
    thermal_cell: CellId,
    reservoir_pump_w: [f64; N],
    prompt_heat_w: [f64; N],
    absorbed_w: [f64; N],
    stored_j: f64,
    state_stored_j: [f64; 2],
}

#[derive(Clone, Debug)]
struct Accounts<const N: usize, const M: usize> {
    duration: CheckedSum,
    incident: [CheckedSum; N],
    transmitted: [CheckedSum; N],
    absorbed: [CheckedSum; N],
    pump_residual: [CheckedSum; N],
    escaped: [CheckedSum; M],
    emitted: [CheckedSum; M],
    unresolved: [CheckedSum; M],
    prompt: CheckedSum,
    dark: CheckedSum,
    crossing: CheckedSum,
    reverse_crossing: CheckedSum,
    relaxation: CheckedSum,
    secondary: CheckedSum,
    recycled: CheckedSum,
    roundoff: CheckedSum,
}

impl<const N: usize, const M: usize> Default for Accounts<N, M> {
    fn default() -> Self {
        Self {
            duration: CheckedSum::default(),
            incident: [CheckedSum::default(); N],
            transmitted: [CheckedSum::default(); N],
            absorbed: [CheckedSum::default(); N],
            pump_residual: [CheckedSum::default(); N],
            escaped: [CheckedSum::default(); M],
            emitted: [CheckedSum::default(); M],
            unresolved: [CheckedSum::default(); M],
            prompt: CheckedSum::default(),
            dark: CheckedSum::default(),
            crossing: CheckedSum::default(),
            reverse_crossing: CheckedSum::default(),
            relaxation: CheckedSum::default(),
            secondary: CheckedSum::default(),
            recycled: CheckedSum::default(),
            roundoff: CheckedSum::default(),
        }
    }
}

#[derive(Clone, Debug)]
struct Pending<const N: usize, const M: usize> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    stored_j: Vec<f64>,
    state_stored_j: Vec<[f64; 2]>,
    mean_cell_emission_w: Vec<(usize, [f64; M])>,
    accounts: Accounts<N, M>,
}

/// Optical pumping, delayed emission and local heat on one native timeline.
/// Cloning retains excitation; a refused advance publishes no partial state.
/// The prepared Beer absorption remains linear and fixed. Independent decay
/// takes constant work per cell/node. Explicit resonant recycling uses the
/// existing bounded positive network propagator, cached for unchanged rates
/// and duration, without lifetime-sized thermal steps or photon generations.
#[derive(Clone, Debug)]
pub struct LuminescentOpticalHeating<const N: usize, const M: usize> {
    source: LuminescenceSource<M>,
    reference: OpticalPowerBudget<N>,
    cells: Vec<Cell<N>>,
    emission_transport: Option<AbsorbingLuminescenceEnvironment<M>>,
    recycling: Option<ResonantFeedback>,
    singlet_triplet: Option<SingletTripletFeedback<M>>,
    schedule: Option<SpectralPowerSchedule<N>>,
    time_s: f64,
    energy_tolerance_j: f64,
    accounts: Accounts<N, M>,
    pump_energy: OpticalEnergyBudget<N>,
    energy: LuminescenceEnergy<M>,
    cell_emission_w: Vec<(usize, [f64; M])>,
    mean_cell_emission_w: Vec<(usize, [f64; M])>,
    pending: Option<Pending<N, M>>,
}

impl<const N: usize, const M: usize> LuminescentOpticalHeating<N, M> {
    /// Prepare one initially unexcited source from the original beam report.
    /// Prompt fluorescence admission verifies report conservation, complete
    /// bindings, source domains and transparent escape in the complete medium.
    /// `energy_tolerance_j` bounds cumulative local and global roundoff; it is
    /// not a physical model or temperature-integration error estimate.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        medium: &AbsorbingGrid<'_>,
        source: LuminescenceSource<M>,
        deposition: BeamPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        power_tolerance_w: f64,
        energy_tolerance_j: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        Self::prepare(
            medium,
            source,
            deposition,
            cell_bindings,
            thermal,
            max_grid_cells,
            power_tolerance_w,
            energy_tolerance_j,
            None,
            checkpoint,
        )
    }

    /// Transport emitted light through an explicitly passive absorbing and
    /// possibly scattering environment. The original pump remains the bound
    /// non-scattering beam. Every pumped emitter needs an operator row and the
    /// same native binding in the environment. Reexcitation is not inferred
    /// from absorption: the declared emission-band absorption channels relax
    /// directly to heat. The prepared environment is shared across clones.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_emission_transport(
        medium: &ScatteringGrid<'_>,
        source: LuminescenceSource<M>,
        deposition: BeamPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        power_tolerance_w: f64,
        energy_tolerance_j: f64,
        environment: AbsorbingLuminescenceEnvironment<M>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        let absorption = medium.absorption_grid();
        if environment.operator().medium_identity != medium.content_identity()
            || absorption.cell_count() > max_grid_cells
        {
            return Err(invalid(
                "luminescence pump and emission transport snapshots differ",
            ));
        }
        let mut previous = None;
        for cell in 0..absorption.cell_count() {
            checkpoint()?;
            if absorption.cell_density_kg_m3(cell) == Some(0.0) {
                continue;
            }
            let law = medium
                .scattering_at_cell(cell)
                .map_err(|_| invalid("invalid original pump scattering law"))?;
            if previous.is_some_and(|old| core::ptr::eq(old, law)) {
                continue;
            }
            for &wavelength in &deposition.wavelengths_nm {
                if law
                    .at(wavelength)
                    .map_err(|_| invalid("pump exceeds scattering support"))?
                    != 0.0
                {
                    return Err(invalid(
                        "luminescence pump scattering requires a spatial scattering deposition owner",
                    ));
                }
            }
            previous = Some(law);
        }
        Self::prepare(
            absorption,
            source,
            deposition,
            cell_bindings,
            thermal,
            max_grid_cells,
            power_tolerance_w,
            energy_tolerance_j,
            Some(environment),
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        medium: &AbsorbingGrid<'_>,
        source: LuminescenceSource<M>,
        deposition: BeamPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        power_tolerance_w: f64,
        energy_tolerance_j: f64,
        emission_transport: Option<AbsorbingLuminescenceEnvironment<M>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        if !nonnegative(energy_tolerance_j) {
            return Err(invalid(
                "finite nonnegative luminescence energy allowance required",
            ));
        }
        validate_fluorescent_medium(
            medium,
            &source.fluorescence,
            &deposition,
            cell_bindings.len(),
            max_grid_cells,
            power_tolerance_w,
            emission_transport.is_none(),
            &mut checkpoint,
        )?;
        if let Some(environment) = &emission_transport {
            if environment.operator().absorption_identity != medium.content_identity()
                || environment.operator().wavelengths_nm != *source.fluorescence.wavelengths_nm()
            {
                return Err(invalid(
                    "luminescence source and secondary transport differ",
                ));
            }
            environment.check_state(thermal, &mut checkpoint)?;
        }
        let pump = OpticalHeating::try_new(
            deposition.clone(),
            cell_bindings,
            thermal,
            power_tolerance_w,
            &mut checkpoint,
        )?;
        let bindings: BTreeMap<_, _> = cell_bindings.iter().copied().collect();
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(deposition.cells.len())
            .map_err(|_| invalid("luminescence cell allocation refused"))?;
        for deposited in &deposition.cells {
            checkpoint()?;
            if deposited.absorbed_power_w.iter().all(|p| *p == 0.0) {
                continue;
            }
            let thermal_cell = bindings[&deposited.cell_index];
            if let Some(environment) = &emission_transport {
                environment.check_emitter(deposited.cell_index, thermal_cell)?;
            }
            source.decay_at_temperature(thermal.temperature_k(thermal_cell)?)?;
            let mut reservoir_pump_w = [0.0; N];
            let mut prompt_heat_w = [0.0; N];
            for lane in 0..N {
                let absorbed = deposited.absorbed_power_w[lane];
                let pump_nm = deposition.wavelengths_nm[lane];
                let stored_nm = source.excitation_wavelength_nm;
                reservoir_pump_w[lane] = ratio_power(absorbed, pump_nm, stored_nm)?;
                // Retain a tiny positive Stokes loss without subtracting watts.
                prompt_heat_w[lane] = ratio_power(absorbed, stored_nm - pump_nm, stored_nm)?;
            }
            cells.push(Cell {
                grid_cell: deposited.cell_index,
                thermal_cell,
                reservoir_pump_w,
                prompt_heat_w,
                absorbed_w: deposited.absorbed_power_w,
                stored_j: 0.0,
                state_stored_j: [0.0; 2],
            });
        }
        let cell_emission_w: Vec<_> = cells.iter().map(|c| (c.grid_cell, [0.0; M])).collect();
        let accounts = Accounts::default();
        let (pump_energy, energy) = accounts.report(0.0)?;
        checkpoint()?;
        Ok(Self {
            source,
            reference: *pump.power(),
            cells,
            emission_transport,
            recycling: None,
            singlet_triplet: None,
            schedule: None,
            time_s: thermal.time_s(),
            energy_tolerance_j,
            accounts,
            pump_energy,
            energy,
            mean_cell_emission_w: cell_emission_w.clone(),
            cell_emission_w,
            pending: None,
        })
    }

    /// Supply rectangular spectral pulses using the existing waveform owner.
    /// Every segment must be constant and no lane may exceed its prepared beam
    /// power. Gaps and times outside support are dark. Linear ramps refuse:
    /// using their mean would preserve fluence but give incorrect decay timing.
    /// The native transition hook splits every pulse edge exactly.
    pub fn with_pump_schedule(
        mut self,
        schedule: SpectralPowerSchedule<N>,
    ) -> Result<Self, OpticalHeatingError> {
        if self.pending.is_some()
            || schedule
                .segments()
                .iter()
                .any(|s| s.start_power_w != s.end_power_w)
            || schedule
                .peak_power_w()
                .iter()
                .zip(self.reference.incident_w)
                .any(|(peak, reference)| *peak > reference)
        {
            return Err(invalid(
                "luminescence pumping requires rectangular pulses within the prepared beam envelope",
            ));
        }
        self.schedule = Some(schedule);
        Ok(self)
    }

    /// Immutable supplied emitting-state law and photon spectrum.
    pub const fn source(&self) -> &LuminescenceSource<M> {
        &self.source
    }
    /// Accepted native thermal clock, seconds.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Cumulative accepted light/heat accounts and current excitation store.
    pub const fn energy(&self) -> &LuminescenceEnergy<M> {
        &self.energy
    }
    /// Original optical pump accounts and actual heat delivered to native cells.
    pub const fn pump_energy(&self) -> &OpticalEnergyBudget<N> {
        &self.pump_energy
    }

    /// Shared source-to-receiver transport and its explicit sampling/tail
    /// allowances. `None` denotes the validated transparent-escape fast path.
    pub fn emission_transport(&self) -> Option<&AbsorbingLuminescenceEnvironment<M>> {
        self.emission_transport.as_ref()
    }

    /// Whether emitted-photon absorption creates equal-energy excitation
    /// instead of following the default passive absorption-to-heat channel.
    pub const fn has_resonant_recycling(&self) -> bool {
        self.recycling.is_some()
    }

    /// Emitted integrated watts per original grid cell at the accepted endpoint.
    /// These remain positive during dark afterglow. No radiance/temperature or
    /// arbitrary spectral-bin width is inferred from these powers.
    pub fn cell_emission_w(&self) -> &[(usize, [f64; M])] {
        &self.cell_emission_w
    }

    /// Mean emitted watts during the last accepted interval, suitable for an
    /// exposure over THAT interval. Empty-history values are zero. They must
    /// not be confused with instantaneous endpoint powers or per-nm spectra.
    pub fn last_mean_cell_emission_w(&self) -> &[(usize, [f64; M])] {
        &self.mean_cell_emission_w
    }

    /// Stored nonthermal excitation of an original optical cell, in joules.
    pub fn stored_energy_j(&self, grid_cell: usize) -> Option<f64> {
        self.cells
            .iter()
            .find(|c| c.grid_cell == grid_cell)
            .map(|c| c.stored_j)
    }

    /// Accepted separate `[singlet, triplet]` nonthermal stores, joules.
    /// The original scalar model uses only the first entry.
    pub fn state_stored_energy_j(&self, grid_cell: usize) -> Option<[f64; 2]> {
        self.cells
            .iter()
            .find(|c| c.grid_cell == grid_cell)
            .map(|c| c.state_stored_j)
    }

    fn check_state(
        &self,
        thermal: &DepositionSimulation<'_>,
        check_clock: bool,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalHeatingError> {
        checkpoint()?;
        if check_clock && thermal.time_s() != self.time_s {
            return Err(invalid("luminescence and thermal clocks differ"));
        }
        for cell in &self.cells {
            checkpoint()?;
            if !thermal.cells().contains_key(&cell.thermal_cell) {
                return Err(OpticalHeatingError::InactiveCell(cell.thermal_cell));
            }
            self.source
                .decay_at_temperature(thermal.temperature_k(cell.thermal_cell)?)?;
            if let Some(coupling) = &self.singlet_triplet {
                coupling.model_at_temperature(
                    &self.source,
                    thermal.temperature_k(cell.thermal_cell)?,
                )?;
            }
        }
        if let Some(environment) = &self.emission_transport {
            environment.check_state(thermal, checkpoint)?;
        }
        Ok(())
    }
}

impl<const N: usize, const M: usize> DepositionObserver for LuminescentOpticalHeating<N, M> {
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
        self.schedule.as_ref().map_or(Ok(None), |s| {
            s.next_transition_time_s(self.time_s)
                .map_err(|_| invalid("invalid luminescence pulse clock"))
        })
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
        self.check_state(thermal, true, checkpoint)?;
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !end_s.is_finite()
            || end_s <= start_s
            || (end_s != start_s + duration_s && duration_s != end_s - start_s)
        {
            return Err(invalid(
                "matching final native luminescence interval required",
            ));
        }
        if let Some(schedule) = &self.schedule
            && schedule
                .next_transition_time_s(start_s)
                .map_err(|_| invalid("invalid luminescence pulse interval"))?
                .is_some_and(|edge| edge < end_s)
        {
            return Err(invalid("luminescence interval crosses a pump transition"));
        }
        let actual = match &self.schedule {
            Some(schedule) => schedule
                .mean_power_w(start_s, duration_s, end_s, &mut *checkpoint)
                .map_err(schedule_error)?,
            None => self.reference.incident_w,
        };
        let mut scale = [0.0; N];
        let mut accounts = self.accounts.clone();
        accounts.duration.add(duration_s)?;
        for lane in 0..N {
            scale[lane] = if self.reference.incident_w[lane] == 0.0 {
                0.0
            } else {
                actual[lane] / self.reference.incident_w[lane]
            };
            if actual[lane] > 0.0 && scale[lane] < f64::MIN_POSITIVE {
                return Err(invalid(
                    "scheduled pump scaling is below supported representation",
                ));
            }
            accounts.incident[lane].add(finite(actual[lane] * duration_s)?)?;
            accounts.transmitted[lane].add(finite(
                self.reference.transmitted_w[lane] * scale[lane] * duration_s,
            )?)?;
            accounts.absorbed[lane].add(finite(
                self.reference.absorbed_w[lane] * scale[lane] * duration_s,
            )?)?;
            accounts.pump_residual[lane].add(finite(
                self.reference.residual_w[lane] * scale[lane] * duration_s,
            )?)?;
        }
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let mut stored_j = Vec::with_capacity(self.cells.len());
        let mut state_stored_j = Vec::with_capacity(self.cells.len());
        let mut mean_cell_emission_w = Vec::with_capacity(self.cells.len());
        let recycled = match &mut self.recycling {
            Some(feedback) => Some(feedback.advance(
                &self.source,
                &self.cells,
                &scale,
                duration_s,
                thermal,
                checkpoint,
            )?),
            None => None,
        };
        for (index, cell) in self.cells.iter().enumerate() {
            checkpoint()?;
            let mut pump = CheckedSum::default();
            let mut prompt = CheckedSum::default();
            let mut absorbed = CheckedSum::default();
            for (lane, factor) in scale.iter().enumerate() {
                pump.add(finite(cell.reservoir_pump_w[lane] * factor)?)?;
                prompt.add(finite(cell.prompt_heat_w[lane] * factor * duration_s)?)?;
                absorbed.add(finite(cell.absorbed_w[lane] * factor * duration_s)?)?;
            }
            let coupled = match &mut self.singlet_triplet {
                Some(coupling) => Some(coupling.advance(
                    &self.source,
                    cell.state_stored_j,
                    pump.value()?,
                    duration_s,
                    thermal.temperature_k(cell.thermal_cell)?,
                    checkpoint,
                )?),
                None => None,
            };
            let step = match (&coupled, &recycled) {
                (Some(step), _) => ExcitationStep {
                    stored_energy_j: step.stored_energy_j,
                    pumped_energy_j: finite(pump.value()? * duration_s)?,
                    radiative_decay_energy_j: 0.0,
                    nonradiative_heat_j: step.nonradiative_heat_j,
                    residual_j: 0.0,
                },
                (None, Some(step)) => ExcitationStep {
                    stored_energy_j: step.stored_energy_j[index],
                    pumped_energy_j: step.pumped_energy_j[index],
                    radiative_decay_energy_j: step.emitted_energy_j[index],
                    nonradiative_heat_j: step.nonradiative_heat_j[index],
                    // Local physical closure is checked below, including the
                    // independently retained incoming recycled photon energy.
                    residual_j: 0.0,
                },
                (None, None) => self
                    .source
                    .decay_at_temperature(thermal.temperature_k(cell.thermal_cell)?)?
                    .advance(cell.stored_j, pump.value()?, duration_s)
                    .map_err(|_| invalid("luminescence reservoir advance refused"))?,
            };
            let prompt = prompt.value()?;
            let (emitted, relaxation, crossing, reverse_crossing) = match &coupled {
                Some(step) => (
                    step.emitted_energy_j,
                    step.relaxation_heat_j,
                    step.crossing_heat_j,
                    step.reverse_crossing_heat_absorbed_j,
                ),
                None => {
                    let (emitted, relaxation) =
                        partition(&self.source, step.radiative_decay_energy_j)?;
                    (emitted, relaxation, 0.0, 0.0)
                }
            };
            let mut heat = CheckedSum::default();
            heat.add(prompt)?;
            heat.add(step.nonradiative_heat_j)?;
            heat.add(relaxation)?;
            heat.add(crossing)?;
            heat.add(-reverse_crossing)?;
            let heat = heat.value()?;
            sources
                .entry(cell.thermal_cell)
                .or_default()
                .add(finite(heat / duration_s)?)?;
            let mut balance = CheckedSum::default();
            for value in [
                cell.stored_j,
                absorbed.value()?,
                -step.stored_energy_j,
                -heat,
            ] {
                balance.add(value)?;
            }
            if let Some(step) = &recycled {
                balance.add(step.reabsorbed_energy_j[index])?;
                accounts.recycled.add(step.reabsorbed_energy_j[index])?;
            }
            let mut mean = [0.0; M];
            for bin in 0..M {
                accounts.emitted[bin].add(emitted[bin])?;
                balance.add(-emitted[bin])?;
                mean[bin] = finite(emitted[bin] / duration_s)?;
            }
            if recycled.is_some() {
                // All resolved absorption already entered the coupled stores.
                // External escape/tail totals are reduced once after this loop.
            } else if let Some(environment) = &self.emission_transport {
                let routed = environment.route(
                    cell.grid_cell,
                    emitted,
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
                for bin in 0..M {
                    accounts.escaped[bin].add(emitted[bin])?;
                }
            }
            accounts.roundoff.add(balance.value()?.abs())?;
            accounts.prompt.add(prompt)?;
            accounts.dark.add(step.nonradiative_heat_j)?;
            accounts.crossing.add(crossing)?;
            accounts.reverse_crossing.add(reverse_crossing)?;
            accounts.relaxation.add(relaxation)?;
            stored_j.push(step.stored_energy_j);
            state_stored_j.push(
                coupled
                    .as_ref()
                    .map_or([step.stored_energy_j, 0.0], |step| {
                        step.state_stored_energy_j
                    }),
            );
            mean_cell_emission_w.push((cell.grid_cell, mean));
        }
        if let Some(step) = &recycled {
            // Recycling is admitted only by the single-line constructor.
            accounts.escaped[0].add(step.escaped_energy_j)?;
            accounts.unresolved[0].add(step.unresolved_energy_j)?;
        }
        if let Some(environment) = &self.emission_transport {
            let mut unresolved = CheckedSum::default();
            for lane in accounts.unresolved {
                unresolved.add(lane.value()?)?;
            }
            if unresolved.value()? > environment.tolerance().max_unresolved_energy_j {
                return Err(invalid(
                    "luminescence unresolved radiation exceeds the supplied energy allowance",
                ));
            }
        }
        if accounts.roundoff.value()? > self.energy_tolerance_j {
            return Err(invalid(
                "luminescence local energy defects exceed the supplied allowance",
            ));
        }
        let sources = sources
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        checkpoint()?;
        // Implicit thermal retries replace tentative work at the SAME accepted
        // start. Every candidate above was rebuilt from the accepted accounts
        // and stores, never from a preceding rejected interval.
        self.pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            stored_j,
            state_stored_j,
            mean_cell_emission_w,
            accounts,
        });
        Ok(sources)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_state(thermal, true, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_state(thermal, false, checkpoint)?;
        let pending = self
            .pending
            .as_ref()
            .ok_or_else(|| invalid("luminescence step lacks prepared interval"))?;
        if pending.start_s != start_s
            || pending.duration_s != duration_s
            || pending.end_s != thermal.time_s()
            || start_s != self.time_s
        {
            return Err(invalid(
                "luminescence interval does not match accepted thermal step",
            ));
        }
        let mut stored = CheckedSum::default();
        let mut cell_emission_w = Vec::with_capacity(self.cells.len());
        for ((cell, &energy), &state_energy) in self
            .cells
            .iter()
            .zip(&pending.stored_j)
            .zip(&pending.state_stored_j)
        {
            checkpoint()?;
            stored.add(energy)?;
            let temperature = thermal.temperature_k(cell.thermal_cell)?;
            let emission = match &self.singlet_triplet {
                Some(coupling) => {
                    coupling.emission_power_w(&self.source, state_energy, temperature)?
                }
                None => {
                    let decay = self.source.decay_at_temperature(temperature)?;
                    let radiative_power =
                        rate_power(energy, decay.quantum_yield(), decay.lifetime_s())?;
                    partition(&self.source, radiative_power)?.0
                }
            };
            cell_emission_w.push((cell.grid_cell, emission));
        }
        let (pump_energy, energy) = pending.accounts.report(stored.value()?)?;
        if energy.residual_j.abs() > self.energy_tolerance_j {
            return Err(invalid(
                "luminescence heat, light and storage do not conserve the pump energy",
            ));
        }
        checkpoint()?;
        let pending = self
            .pending
            .take()
            .expect("validated prepared luminescence interval");
        for ((cell, energy), states) in self
            .cells
            .iter_mut()
            .zip(pending.stored_j)
            .zip(pending.state_stored_j)
        {
            cell.stored_j = energy;
            cell.state_stored_j = states;
        }
        self.accounts = pending.accounts;
        self.pump_energy = pump_energy;
        self.energy = energy;
        self.cell_emission_w = cell_emission_w;
        self.mean_cell_emission_w = pending.mean_cell_emission_w;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

impl<const N: usize, const M: usize> Accounts<N, M> {
    fn report(
        &self,
        stored_j: f64,
    ) -> Result<(OpticalEnergyBudget<N>, LuminescenceEnergy<M>), OpticalHeatingError> {
        let mut energy = LuminescenceEnergy {
            duration_s: self.duration.value()?,
            escaped_j: [0.0; M],
            emitted_j: [0.0; M],
            unresolved_radiation_j: [0.0; M],
            prompt_heat_j: self.prompt.value()?,
            nonradiative_heat_j: self.dark.value()?,
            crossing_heat_j: self.crossing.value()?,
            reverse_crossing_heat_absorbed_j: self.reverse_crossing.value()?,
            relaxation_heat_j: self.relaxation.value()?,
            secondary_absorption_heat_j: self.secondary.value()?,
            recycled_excitation_j: self.recycled.value()?,
            stored_excitation_j: stored_j,
            residual_j: 0.0,
            absolute_roundoff_j: self.roundoff.value()?,
        };
        let mut heat = CheckedSum::default();
        for value in [
            energy.prompt_heat_j,
            energy.nonradiative_heat_j,
            energy.crossing_heat_j,
            -energy.reverse_crossing_heat_absorbed_j,
            energy.relaxation_heat_j,
            energy.secondary_absorption_heat_j,
        ] {
            heat.add(value)?;
        }
        let mut pump = OpticalEnergyBudget {
            duration_s: energy.duration_s,
            incident_j: [0.0; N],
            transmitted_j: [0.0; N],
            absorbed_j: [0.0; N],
            scattered_j: [0.0; N],
            residual_j: [0.0; N],
            thermal_input_j: heat.value()?,
        };
        let mut balance = CheckedSum::default();
        for lane in 0..N {
            pump.incident_j[lane] = self.incident[lane].value()?;
            pump.transmitted_j[lane] = self.transmitted[lane].value()?;
            pump.absorbed_j[lane] = self.absorbed[lane].value()?;
            pump.residual_j[lane] = self.pump_residual[lane].value()?;
            balance.add(pump.absorbed_j[lane])?;
        }
        balance.add(-stored_j)?;
        balance.add(-pump.thermal_input_j)?;
        for bin in 0..M {
            energy.emitted_j[bin] = self.emitted[bin].value()?;
            energy.escaped_j[bin] = self.escaped[bin].value()?;
            energy.unresolved_radiation_j[bin] = self.unresolved[bin].value()?;
            balance.add(-energy.escaped_j[bin])?;
            balance.add(-energy.unresolved_radiation_j[bin])?;
        }
        energy.residual_j = balance.value()?;
        Ok((pump, energy))
    }
}

fn partition<const M: usize>(
    source: &LuminescenceSource<M>,
    released: f64,
) -> Result<([f64; M], f64), OpticalHeatingError> {
    let mut emitted = [0.0; M];
    let mut relaxation = CheckedSum::default();
    for (bin, value) in emitted.iter_mut().enumerate() {
        if source.fluorescence.photon_fractions()[bin] == 0.0 {
            continue;
        }
        let weight = finite(released * source.fluorescence.photon_fractions()[bin])?;
        let wavelength = source.fluorescence.wavelengths_nm()[bin];
        *value = ratio_power(weight, source.excitation_wavelength_nm, wavelength)?;
        relaxation.add(ratio_power(
            weight,
            wavelength - source.excitation_wavelength_nm,
            wavelength,
        )?)?;
    }
    Ok((emitted, relaxation.value()?))
}

fn rate_power(energy: f64, fraction: f64, lifetime: f64) -> Result<f64, OpticalHeatingError> {
    if energy == 0.0 || fraction == 0.0 {
        return Ok(0.0);
    }
    let rate = (energy * fraction) / lifetime;
    if rate.is_finite() && rate > 0.0 {
        return Ok(rate);
    }
    finite(fs_math::det::exp(
        fs_math::det::ln(energy) + fs_math::det::ln(fraction) - fs_math::det::ln(lifetime),
    ))
}

fn schedule_error(error: OpticalFeedbackError) -> OpticalHeatingError {
    match error {
        OpticalFeedbackError::Thermal(error) => error.into(),
        _ => invalid("luminescence pulse integration refused"),
    }
}
