//! Photon-conservative pumping of an ordered finite-population optical column.
//!
//! Each cell is spatially homogeneous. A lane's active optical depth is its
//! current ground inventory times the sum of its two authored depth-per-emitter
//! coefficients. For a ray segment these are `cross_section * length / volume`;
//! the caller supplies the original cell volume and advances its full inventory
//! only once. A narrow beam therefore uses a cell-average population closure,
//! whose spatial accuracy requires cell refinement. No beam area is invented.
//!
//! A bounded scalar solve makes each interval's Beer field consistent with the
//! EXACT integrated ground residence of the original stiff three-state kinetic
//! owner. Actual pump events, independent transmitted photons, and passive host
//! absorption close each lane. Full-column time refinement resolves changes in
//! opacity and the time correlation of upstream bleaching with downstream
//! exposure. Conservation alone is not an accuracy certificate. Only the finer
//! accepted trajectory contributes physical accounts; no extrapolation,
//! clipping, normalization, or residual-to-heat correction is applied.
//!
//! Pump photon flight is quasistatic. The column has one directed pass, fixed
//! coefficients during this call, and no pump scattering, reflection, stimulated
//! emission, excited-state absorption, or active recycling of emitted light.
//! Thermal owners refresh models and optical coefficients on their own clock.

use super::saturation::{SaturableSingletTripletPlan, SaturableSingletTripletStep};
use super::{SingletTripletDecay, SingletTripletError, finite, nonnegative, positive_product};
use fs_la::decay::{DecayError, DecayLimits};
use fs_math::decay::{attenuate, decayed_fraction};

/// One homogeneous original cell in entry-to-exit order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturablePumpCell<const B: usize> {
    pub model: SingletTripletDecay,
    /// Lane then destination `[singlet, triplet]`, dimensionless per emitter.
    /// Optical geometry compiles this as molecular cross section times ray
    /// length divided by the full original population's cell volume.
    pub optical_depth_per_ground: [[f64; 2]; B],
    /// Fixed passive host absorption depth, separate from finite emitters.
    pub passive_optical_depth: [f64; B],
}

/// Bounded scalar and temporal work; physical balances use `kinetic` accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpColumnLimits {
    pub kinetic: DecayLimits,
    pub maximum_root_iterations: usize,
    /// Largest accepted/trial uniform temporal grid. Nonlinear comparison
    /// needs two grids; their sizes double deterministically.
    pub maximum_substeps: usize,
    /// Aggregate local kinetic preparations over ALL scalar trials and grids.
    pub maximum_cell_evaluations: usize,
    /// Componentwise temporal convergence indicator, not a rigorous enclosure.
    pub relative_tolerance: f64,
    pub absolute_energy_tolerance_j: f64,
}

impl Default for PumpColumnLimits {
    fn default() -> Self {
        Self {
            kinetic: DecayLimits::default(),
            maximum_root_iterations: 64,
            maximum_substeps: 1024,
            maximum_cell_evaluations: 262_144,
            relative_tolerance: 1e-6,
            absolute_energy_tolerance_j: 0.0,
        }
    }
}

/// Integrated response of one original cell, with spectral and target accounts
/// kept distinct. `kinetics.absorbed_pump_energy_j` is ordered by destination;
/// `absorbed_energy_j` below is ordered by original pump lane.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpColumnCellStep<const B: usize> {
    pub kinetics: SaturableSingletTripletStep,
    pub incident_photons: [f64; B],
    pub transmitted_photons: [f64; B],
    pub absorbed_photons: [f64; B],
    pub passive_absorbed_photons: [f64; B],
    pub incident_energy_j: [f64; B],
    pub transmitted_energy_j: [f64; B],
    pub absorbed_energy_j: [f64; B],
    pub passive_absorbed_energy_j: [f64; B],
    pub photon_residual: [f64; B],
    pub passive_heat_j: f64,
}

/// Independent whole-column candidate. Energy closure counts external incident
/// and transmitted pump photons; internal cell-to-cell transport is not an
/// additional source. Emitted state energy still precedes spectral conversion.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpColumnStep<const B: usize> {
    pub cells: Vec<PumpColumnCellStep<B>>,
    pub incident_photons: [f64; B],
    pub transmitted_photons: [f64; B],
    pub absorbed_photons: [f64; B],
    pub passive_absorbed_photons: [f64; B],
    pub incident_energy_j: [f64; B],
    pub transmitted_energy_j: [f64; B],
    pub absorbed_energy_j: [f64; B],
    pub passive_absorbed_energy_j: [f64; B],
    pub photon_residual: [f64; B],
    pub residual_j: f64,
    pub substeps: usize,
    pub cell_evaluations: usize,
    pub max_energy_discrepancy_j: f64,
}

/// Immutable authored column; no population or incident field is cached.
#[derive(Clone, Debug)]
pub struct SaturablePumpColumn<const B: usize> {
    cells: Vec<SaturablePumpCell<B>>,
    photon_energies_j: [f64; B],
}

impl<const B: usize> SaturablePumpColumn<B> {
    pub fn try_new(
        cells: &[SaturablePumpCell<B>],
        photon_energies_j: [f64; B],
        maximum_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        poll(&mut cancelled)?;
        if B == 0
            || B > 64
            || maximum_cells == 0
            || !photon_energies_j
                .into_iter()
                .all(|x| x.is_finite() && x > 0.0)
        {
            return Err(SingletTripletError::Invalid(
                "pump column lanes, energies or cell limit",
            ));
        }
        cap("pump column cells", cells.len(), maximum_cells)?;
        for cell in cells {
            poll(&mut cancelled)?;
            if !cell
                .optical_depth_per_ground
                .iter()
                .flatten()
                .copied()
                .all(nonnegative)
                || !cell.passive_optical_depth.into_iter().all(nonnegative)
            {
                return Err(SingletTripletError::Invalid("pump column optical depths"));
            }
            for (lane, energy) in photon_energies_j.iter().copied().enumerate() {
                finite(
                    cell.optical_depth_per_ground[lane].iter().sum(),
                    "summed active depth",
                )?;
                for target in 0..2 {
                    if cell.optical_depth_per_ground[lane][target] > 0.0
                        && energy < cell.model.excitation_energies_j()[target]
                    {
                        return Err(SingletTripletError::Invalid(
                            "pump photon cannot populate its destination",
                        ));
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            cells: cells.to_vec(),
            photon_energies_j,
        })
    }

    pub fn cells(&self) -> &[SaturablePumpCell<B>] {
        &self.cells
    }
    pub const fn photon_energies_j(&self) -> &[f64; B] {
        &self.photon_energies_j
    }

    /// Constant incident photon rates, per lane, and fixed constitutive models
    /// during this interval. Every candidate is independent of the inputs.
    pub fn advance(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: [f64; B],
        duration_s: f64,
        limits: PumpColumnLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PumpColumnStep<B>, SingletTripletError> {
        self.validate_advance(
            populations,
            incident_photon_rates_s,
            duration_s,
            limits,
            &mut cancelled,
        )?;
        let nonlinear = duration_s > 0.0
            && self.cells.iter().any(|cell| {
                (0..B).any(|lane| {
                    incident_photon_rates_s[lane] > 0.0
                        && cell.optical_depth_per_ground[lane].iter().any(|&x| x > 0.0)
                })
            });
        let mut evaluations = 0usize;
        let mut segments = 1usize;
        let mut previous: Option<Vec<f64>> = None;
        loop {
            poll(&mut cancelled)?;
            cap(
                "pump column temporal substeps",
                segments,
                limits.maximum_substeps,
            )?;
            let mut candidate = self.uniform(
                populations,
                incident_photon_rates_s,
                duration_s,
                segments,
                limits,
                &mut evaluations,
                &mut cancelled,
            )?;
            candidate.cell_evaluations = evaluations;
            if !nonlinear {
                return Ok(candidate);
            }
            let components = self.components(&candidate, duration_s)?;
            if let Some(old) = &previous {
                let mut converged = true;
                let mut discrepancy = 0.0_f64;
                for (&before, &after) in old.iter().zip(&components) {
                    let difference = (after - before).abs();
                    discrepancy = discrepancy.max(difference);
                    converged &= difference
                        <= limits.absolute_energy_tolerance_j
                            + limits.relative_tolerance * before.abs().max(after.abs());
                }
                if converged {
                    candidate.max_energy_discrepancy_j = discrepancy;
                    poll(&mut cancelled)?;
                    return Ok(candidate);
                }
            }
            previous = Some(components);
            segments = segments
                .checked_mul(2)
                .ok_or(SingletTripletError::Unrepresentable(
                    "pump column refinement counter",
                ))?;
        }
    }

    fn validate_advance(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: [f64; B],
        duration_s: f64,
        limits: PumpColumnLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        poll(&mut cancelled)?;
        if populations.len() != self.cells.len()
            || !populations.iter().flatten().copied().all(nonnegative)
            || !incident_photon_rates_s.into_iter().all(nonnegative)
            || !nonnegative(duration_s)
            || limits.maximum_root_iterations < 2
            || limits.maximum_substeps == 0
            || limits.maximum_cell_evaluations == 0
            || !nonnegative(limits.relative_tolerance)
            || limits.relative_tolerance >= 1.0
            || !nonnegative(limits.absolute_energy_tolerance_j)
            || (limits.relative_tolerance == 0.0 && limits.absolute_energy_tolerance_j == 0.0)
        {
            return Err(SingletTripletError::Invalid(
                "pump column populations, illumination or refinement limits",
            ));
        }
        // Admit the original kinetic limits even for an empty/transparent path.
        SaturableSingletTripletPlan::prepare(
            SingletTripletDecay::try_new([1.0; 2], [0.0; 2], [0.0; 2], [0.0; 2])?,
            [0.0; 2],
            [0.0; 2],
            0.0,
            limits.kinetic,
            &mut cancelled,
        )?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn uniform(
        &self,
        initial: &[[f64; 3]],
        incident_rate: [f64; B],
        duration: f64,
        segments: usize,
        limits: PumpColumnLimits,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PumpColumnStep<B>, SingletTripletError> {
        let mut populations = initial.to_vec();
        let mut totals: Vec<_> = self.cells.iter().map(|_| CellTotals::<B>::new()).collect();
        let mut transmitted = [Sum::default(); B];
        let mut start = 0.0;
        for segment in 0..segments {
            poll(cancelled)?;
            let end = duration * ((segment + 1) as f64 / segments as f64);
            let h = end - start;
            if duration > 0.0 && h <= 0.0 {
                return Err(SingletTripletError::Unrepresentable(
                    "pump column temporal resolution",
                ));
            }
            let mut incoming = incident_rate;
            for (index, cell) in self.cells.iter().enumerate() {
                poll(cancelled)?;
                let trial = solve_cell(
                    cell,
                    populations[index],
                    incoming,
                    self.photon_energies_j,
                    h,
                    limits,
                    evaluations,
                    cancelled,
                )?;
                populations[index] = trial.step.kinetics.populations;
                incoming = trial.transmitted_rates;
                totals[index].add(&trial.step)?;
            }
            for lane in 0..B {
                transmitted[lane].add(positive_product(
                    incoming[lane],
                    h,
                    "transmitted photons",
                )?)?;
            }
            start = end;
        }
        self.finish_history(
            initial,
            &populations,
            totals,
            transmitted,
            incident_rate,
            duration,
            segments,
            limits,
            *evaluations,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_history(
        &self,
        initial: &[[f64; 3]],
        populations: &[[f64; 3]],
        totals: Vec<CellTotals<B>>,
        transmitted: [Sum; B],
        incident_rate: [f64; B],
        duration: f64,
        segments: usize,
        limits: PumpColumnLimits,
        evaluations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PumpColumnStep<B>, SingletTripletError> {
        let mut cells = Vec::with_capacity(self.cells.len());
        for (index, total) in totals.into_iter().enumerate() {
            poll(cancelled)?;
            cells.push(total.finish(
                self.cells[index].model,
                initial[index],
                populations[index],
                self.photon_energies_j,
                limits.kinetic.relative_tolerance,
            )?);
        }
        let mut incident = [0.0; B];
        let mut outgoing = [0.0; B];
        let mut absorbed = [Sum::default(); B];
        let mut passive = [Sum::default(); B];
        for lane in 0..B {
            incident[lane] = positive_product(incident_rate[lane], duration, "incident photons")?;
            outgoing[lane] = transmitted[lane].value()?;
        }
        for cell in &cells {
            for lane in 0..B {
                absorbed[lane].add(cell.absorbed_photons[lane])?;
                passive[lane].add(cell.passive_absorbed_photons[lane])?;
            }
        }
        let absorbed = values(absorbed)?;
        let passive = values(passive)?;
        let photon_residual = photon_balance(incident, outgoing, absorbed, passive)?;
        admit_photons(
            incident,
            outgoing,
            absorbed,
            passive,
            photon_residual,
            limits.kinetic.relative_tolerance,
        )?;
        let incident_energy_j = energies(incident, self.photon_energies_j)?;
        let transmitted_energy_j = energies(outgoing, self.photon_energies_j)?;
        let absorbed_energy_j = energies(absorbed, self.photon_energies_j)?;
        let passive_absorbed_energy_j = energies(passive, self.photon_energies_j)?;
        let mut residual = Sum::default();
        let mut scale = Sum::default();
        for lane in 0..B {
            residual.add(incident_energy_j[lane])?;
            residual.add(-transmitted_energy_j[lane])?;
            scale.add(incident_energy_j[lane])?;
            scale.add(transmitted_energy_j[lane])?;
        }
        let mut heat = Sum::default();
        for (index, cell) in cells.iter().enumerate() {
            let state_energies = self.cells[index].model.excitation_energies_j();
            for state in 0..2 {
                let old = positive_product(
                    initial[index][state + 1],
                    state_energies[state],
                    "initial excitation",
                )?;
                residual.add(old)?;
                scale.add(old)?;
                for value in [
                    cell.kinetics.stored_energy_j[state],
                    cell.kinetics.radiative_decay_energy_j[state],
                ] {
                    residual.add(-value)?;
                    scale.add(value)?;
                }
                heat.add(cell.kinetics.pump_relaxation_heat_j[state])?;
                heat.add(cell.kinetics.nonradiative_heat_j[state])?;
            }
            heat.add(cell.kinetics.crossing_heat_j)?;
            heat.add(-cell.kinetics.reverse_crossing_heat_absorbed_j)?;
            heat.add(cell.kinetics.annihilation_heat_j)?;
            heat.add(cell.passive_heat_j)?;
        }
        let heat = heat.value()?;
        residual.add(-heat)?;
        scale.add(heat.abs())?;
        let residual_j = residual.value()?;
        admit_energy(
            residual_j,
            scale.value()?,
            limits.kinetic.relative_tolerance,
        )?;
        Ok(PumpColumnStep {
            cells,
            incident_photons: incident,
            transmitted_photons: outgoing,
            absorbed_photons: absorbed,
            passive_absorbed_photons: passive,
            incident_energy_j,
            transmitted_energy_j,
            absorbed_energy_j,
            passive_absorbed_energy_j,
            photon_residual,
            residual_j,
            substeps: segments,
            cell_evaluations: evaluations,
            max_energy_discrepancy_j: 0.0,
        })
    }

    fn components(
        &self,
        step: &PumpColumnStep<B>,
        duration: f64,
    ) -> Result<Vec<f64>, SingletTripletError> {
        let mut result = Vec::new();
        for (model, cell) in self.cells.iter().zip(&step.cells) {
            let k = &cell.kinetics;
            result.extend(cell.transmitted_energy_j);
            result.extend(cell.absorbed_energy_j);
            result.extend(cell.passive_absorbed_energy_j);
            result.extend(k.stored_energy_j);
            result.extend(k.radiative_decay_energy_j);
            result.extend(k.nonradiative_heat_j);
            result.extend(k.pump_relaxation_heat_j);
            result.push(k.crossing_heat_j);
            result.push(k.reverse_crossing_heat_absorbed_j);
            result.push(k.annihilation_heat_j);
            let energies = model.model.excitation_energies_j();
            result.push(positive_product(
                k.annihilation_events,
                energies[1],
                "pair convergence weight",
            )?);
            result.push(positive_product(
                k.annihilation_singlets,
                energies[0],
                "pair singlet convergence weight",
            )?);
            result.push(positive_product(
                k.populations[0],
                energies[0],
                "ground convergence weight",
            )?);
            for target in 0..2 {
                // The visible endpoint power cannot hide behind a tiny fast
                // state's energy. Weight it by the original whole interval.
                result.push(product3(
                    k.stored_energy_j[target],
                    model.model.radiative_rates_s()[target],
                    duration,
                )?);
                for events in [
                    k.radiative_events[target],
                    k.nonradiative_events[target],
                    k.crossing_events[target],
                ] {
                    result.push(positive_product(
                        events,
                        energies[target],
                        "event convergence weight",
                    )?);
                }
            }
        }
        result.extend(step.transmitted_energy_j);
        Ok(result)
    }
}

#[path = "pump_column/cell.rs"]
pub(super) mod cell;

#[path = "pump_column/annihilation.rs"]
mod annihilation;

struct Trial<const B: usize> {
    step: PumpColumnCellStep<B>,
    transmitted_rates: [f64; B],
}

#[allow(clippy::too_many_arguments)]
fn solve_cell<const B: usize>(
    cell: &SaturablePumpCell<B>,
    populations: [f64; 3],
    incoming: [f64; B],
    photon_energies: [f64; B],
    duration: f64,
    limits: PumpColumnLimits,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial<B>, SingletTripletError> {
    let channels: [_; B] = std::array::from_fn(|lane| cell::Channel {
        incident_rate: incoming[lane],
        energy: photon_energies[lane],
        depth: cell.optical_depth_per_ground[lane],
        passive_depth: cell.passive_optical_depth[lane],
    });
    let trial = cell::solve(
        cell.model,
        populations,
        &channels,
        duration,
        limits,
        evaluations,
        cancelled,
    )?;
    let incident = std::array::from_fn(|lane| trial.channels[lane].incident);
    let transmitted = std::array::from_fn(|lane| trial.channels[lane].transmitted);
    let absorbed = std::array::from_fn(|lane| trial.channels[lane].absorbed);
    let passive = std::array::from_fn(|lane| trial.channels[lane].passive);
    Ok(Trial {
        step: PumpColumnCellStep {
            kinetics: trial.kinetics,
            incident_photons: incident,
            transmitted_photons: transmitted,
            absorbed_photons: absorbed,
            passive_absorbed_photons: passive,
            incident_energy_j: energies(incident, photon_energies)?,
            transmitted_energy_j: energies(transmitted, photon_energies)?,
            absorbed_energy_j: energies(absorbed, photon_energies)?,
            passive_absorbed_energy_j: energies(passive, photon_energies)?,
            photon_residual: std::array::from_fn(|lane| trial.channels[lane].residual),
            passive_heat_j: trial.passive_heat_j,
        },
        transmitted_rates: std::array::from_fn(|lane| trial.channels[lane].transmitted_rate),
    })
}

#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    pub(super) fn add(&mut self, value: f64) -> Result<(), SingletTripletError> {
        let (sum, error) = fs_math::eft::two_sum(self.value, value);
        self.value = finite(sum, "pump column partial sum")?;
        self.correction = finite(self.correction + error, "pump column sum correction")?;
        Ok(())
    }
    pub(super) fn value(self) -> Result<f64, SingletTripletError> {
        finite(self.value + self.correction, "pump column sum")
    }
}

struct CellTotals<const B: usize> {
    incident: [Sum; B],
    transmitted: [Sum; B],
    absorbed: [Sum; B],
    passive: [Sum; B],
    kinetics: KineticTotals,
}

impl<const B: usize> CellTotals<B> {
    fn new() -> Self {
        Self {
            incident: [Sum::default(); B],
            transmitted: [Sum::default(); B],
            absorbed: [Sum::default(); B],
            passive: [Sum::default(); B],
            kinetics: KineticTotals::new(),
        }
    }

    fn add(&mut self, step: &PumpColumnCellStep<B>) -> Result<(), SingletTripletError> {
        add_values(&mut self.incident, step.incident_photons)?;
        add_values(&mut self.transmitted, step.transmitted_photons)?;
        add_values(&mut self.absorbed, step.absorbed_photons)?;
        add_values(&mut self.passive, step.passive_absorbed_photons)?;
        self.kinetics.add(&step.kinetics)?;
        Ok(())
    }

    fn finish(
        self,
        model: SingletTripletDecay,
        initial: [f64; 3],
        populations: [f64; 3],
        photon_energies: [f64; B],
        tolerance: f64,
    ) -> Result<PumpColumnCellStep<B>, SingletTripletError> {
        let kinetics = self
            .kinetics
            .finish(model, initial, populations, tolerance)?;
        let incident = values(self.incident)?;
        let transmitted = values(self.transmitted)?;
        let absorbed = values(self.absorbed)?;
        let passive = values(self.passive)?;
        let photon_residual = photon_balance(incident, transmitted, absorbed, passive)?;
        admit_photons(
            incident,
            transmitted,
            absorbed,
            passive,
            photon_residual,
            tolerance,
        )?;
        let absorbed_energy_j = energies(absorbed, photon_energies)?;
        let passive_absorbed_energy_j = energies(passive, photon_energies)?;
        let mut passive_heat = Sum::default();
        let mut target_minus_spectral = Sum::default();
        let mut pump_scale = Sum::default();
        for value in kinetics.absorbed_pump_energy_j {
            target_minus_spectral.add(value)?;
            pump_scale.add(value)?;
        }
        for lane in 0..B {
            passive_heat.add(passive_absorbed_energy_j[lane])?;
            target_minus_spectral.add(-absorbed_energy_j[lane])?;
            pump_scale.add(absorbed_energy_j[lane])?;
        }
        admit_energy(
            target_minus_spectral.value()?,
            pump_scale.value()?,
            tolerance,
        )?;
        Ok(PumpColumnCellStep {
            kinetics,
            incident_photons: incident,
            transmitted_photons: transmitted,
            absorbed_photons: absorbed,
            passive_absorbed_photons: passive,
            incident_energy_j: energies(incident, photon_energies)?,
            transmitted_energy_j: energies(transmitted, photon_energies)?,
            absorbed_energy_j,
            passive_absorbed_energy_j,
            photon_residual,
            passive_heat_j: passive_heat.value()?,
        })
    }
}

/// Shared compensated physical history; one entry per original emitter cell.
pub(super) struct KineticTotals {
    pump_events: [Sum; 2],
    pump_energy: [Sum; 2],
    pump_heat: [Sum; 2],
    radiative_events: [Sum; 2],
    dark_events: [Sum; 2],
    crossing_events: [Sum; 2],
    radiative_energy: [Sum; 2],
    dark_heat: [Sum; 2],
    crossing_heat: Sum,
    reverse_heat: Sum,
    annihilation_events: Sum,
    annihilation_singlets: Sum,
    annihilation_heat: Sum,
}

impl KineticTotals {
    pub(super) fn new() -> Self {
        Self {
            pump_events: [Sum::default(); 2],
            pump_energy: [Sum::default(); 2],
            pump_heat: [Sum::default(); 2],
            radiative_events: [Sum::default(); 2],
            dark_events: [Sum::default(); 2],
            crossing_events: [Sum::default(); 2],
            radiative_energy: [Sum::default(); 2],
            dark_heat: [Sum::default(); 2],
            crossing_heat: Sum::default(),
            reverse_heat: Sum::default(),
            annihilation_events: Sum::default(),
            annihilation_singlets: Sum::default(),
            annihilation_heat: Sum::default(),
        }
    }

    pub(super) fn add(
        &mut self,
        k: &SaturableSingletTripletStep,
    ) -> Result<(), SingletTripletError> {
        add_values(&mut self.pump_events, k.pump_events)?;
        add_values(&mut self.pump_energy, k.absorbed_pump_energy_j)?;
        add_values(&mut self.pump_heat, k.pump_relaxation_heat_j)?;
        add_values(&mut self.radiative_events, k.radiative_events)?;
        add_values(&mut self.dark_events, k.nonradiative_events)?;
        add_values(&mut self.crossing_events, k.crossing_events)?;
        add_values(&mut self.radiative_energy, k.radiative_decay_energy_j)?;
        add_values(&mut self.dark_heat, k.nonradiative_heat_j)?;
        self.crossing_heat.add(k.crossing_heat_j)?;
        self.reverse_heat.add(k.reverse_crossing_heat_absorbed_j)?;
        self.annihilation_events.add(k.annihilation_events)?;
        self.annihilation_singlets.add(k.annihilation_singlets)?;
        self.annihilation_heat.add(k.annihilation_heat_j)?;
        Ok(())
    }

    pub(super) fn finish(
        self,
        model: SingletTripletDecay,
        initial: [f64; 3],
        populations: [f64; 3],
        tolerance: f64,
    ) -> Result<SaturableSingletTripletStep, SingletTripletError> {
        let mut kinetics = SaturableSingletTripletStep {
            populations,
            pump_events: values(self.pump_events)?,
            absorbed_pump_energy_j: values(self.pump_energy)?,
            pump_relaxation_heat_j: values(self.pump_heat)?,
            radiative_events: values(self.radiative_events)?,
            nonradiative_events: values(self.dark_events)?,
            crossing_events: values(self.crossing_events)?,
            stored_energy_j: energies(
                [populations[1], populations[2]],
                model.excitation_energies_j(),
            )?,
            radiative_decay_energy_j: values(self.radiative_energy)?,
            nonradiative_heat_j: values(self.dark_heat)?,
            crossing_heat_j: self.crossing_heat.value()?,
            reverse_crossing_heat_absorbed_j: self.reverse_heat.value()?,
            annihilation_events: self.annihilation_events.value()?,
            annihilation_singlets: self.annihilation_singlets.value()?,
            annihilation_heat_j: self.annihilation_heat.value()?,
            population_residual: 0.0,
            residual_j: 0.0,
        };
        super::saturation::annihilation::admit_excited(initial, &kinetics, tolerance)?;
        let mut population_residual = Sum::default();
        let mut population_scale = Sum::default();
        for i in 0..3 {
            population_residual.add(populations[i])?;
            population_residual.add(-initial[i])?;
            population_scale.add(populations[i])?;
            population_scale.add(initial[i])?;
        }
        kinetics.population_residual = population_residual.value()?;
        let allowance = tolerance * population_scale.value()?;
        if kinetics.population_residual.abs() > allowance {
            return Err(SingletTripletError::PopulationBalance {
                residual: kinetics.population_residual,
                allowance,
            });
        }
        kinetics.residual_j = kinetic_energy_residual(model, initial, &kinetics, tolerance)?;
        Ok(kinetics)
    }
}

pub(super) fn values<const N: usize>(sums: [Sum; N]) -> Result<[f64; N], SingletTripletError> {
    let mut result = [0.0; N];
    for (destination, sum) in result.iter_mut().zip(sums) {
        *destination = sum.value()?;
    }
    Ok(result)
}

pub(super) fn add_values<const N: usize>(
    sums: &mut [Sum; N],
    values: [f64; N],
) -> Result<(), SingletTripletError> {
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value)?;
    }
    Ok(())
}

pub(super) fn energies<const N: usize>(
    photons: [f64; N],
    quantum: [f64; N],
) -> Result<[f64; N], SingletTripletError> {
    let mut result = [0.0; N];
    for lane in 0..N {
        result[lane] = positive_product(photons[lane], quantum[lane], "spectral photon energy")?;
    }
    Ok(result)
}

pub(super) fn photon_balance<const N: usize>(
    incident: [f64; N],
    transmitted: [f64; N],
    absorbed: [f64; N],
    passive: [f64; N],
) -> Result<[f64; N], SingletTripletError> {
    let mut result = [0.0; N];
    for lane in 0..N {
        let mut sum = Sum::default();
        for value in [
            incident[lane],
            -transmitted[lane],
            -absorbed[lane],
            -passive[lane],
        ] {
            sum.add(value)?;
        }
        result[lane] = sum.value()?;
    }
    Ok(result)
}

pub(super) fn admit_photons<const N: usize>(
    incident: [f64; N],
    transmitted: [f64; N],
    absorbed: [f64; N],
    passive: [f64; N],
    residual: [f64; N],
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    for lane in 0..N {
        let mut scale = Sum::default();
        for value in [
            incident[lane],
            transmitted[lane],
            absorbed[lane],
            passive[lane],
        ] {
            scale.add(value)?;
        }
        let allowance = tolerance * scale.value()?;
        if residual[lane].abs() > allowance {
            return Err(SingletTripletError::PopulationBalance {
                residual: residual[lane],
                allowance,
            });
        }
    }
    Ok(())
}

pub(super) fn admit_energy(
    residual_j: f64,
    scale: f64,
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    let allowance_j = tolerance * scale;
    if residual_j.abs() > allowance_j {
        return Err(SingletTripletError::EnergyBalance {
            residual_j,
            allowance_j,
        });
    }
    Ok(())
}

pub(super) fn kinetic_energy_residual(
    model: SingletTripletDecay,
    initial: [f64; 3],
    kinetics: &SaturableSingletTripletStep,
    tolerance: f64,
) -> Result<f64, SingletTripletError> {
    let mut residual = Sum::default();
    let mut scale = Sum::default();
    let mut heat = Sum::default();
    for state in 0..2 {
        for value in [
            positive_product(
                initial[state + 1],
                model.excitation_energies_j()[state],
                "initial excitation",
            )?,
            kinetics.absorbed_pump_energy_j[state],
            -kinetics.stored_energy_j[state],
            -kinetics.radiative_decay_energy_j[state],
        ] {
            residual.add(value)?;
            scale.add(value.abs())?;
        }
        heat.add(kinetics.pump_relaxation_heat_j[state])?;
        heat.add(kinetics.nonradiative_heat_j[state])?;
    }
    heat.add(kinetics.crossing_heat_j)?;
    heat.add(-kinetics.reverse_crossing_heat_absorbed_j)?;
    heat.add(kinetics.annihilation_heat_j)?;
    let heat = heat.value()?;
    residual.add(-heat)?;
    scale.add(heat.abs())?;
    let residual = residual.value()?;
    admit_energy(residual, scale.value()?, tolerance)?;
    Ok(residual)
}

// Avoid needlessly overflowing an intermediate when the complete positive
// product is representable. Required positive products may not vanish.
pub(super) fn product3(a: f64, b: f64, c: f64) -> Result<f64, SingletTripletError> {
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return Ok(0.0);
    }
    for (x, y, z) in [(a, b, c), (a, c, b), (b, c, a)] {
        let product = (x * y) * z;
        if product.is_finite() && product > 0.0 {
            return Ok(product);
        }
    }
    Err(SingletTripletError::Unrepresentable(
        "pump column positive product",
    ))
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SingletTripletError> {
    if cancelled() {
        Err(DecayError::Cancelled.into())
    } else {
        Ok(())
    }
}

pub(super) fn cap(
    quantity: &'static str,
    required: usize,
    maximum: usize,
) -> Result<(), SingletTripletError> {
    if required > maximum {
        Err(DecayError::Budget {
            quantity,
            required: required as u64,
            maximum: maximum as u64,
        }
        .into())
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "pump_column/tests.rs"]
mod tests;
