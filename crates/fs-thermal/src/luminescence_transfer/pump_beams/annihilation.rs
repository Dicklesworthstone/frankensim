//! Pair encounters and all simultaneous pumps share one positive stage.

use super::super::saturation::annihilation::{
    self as pairs, AnnihilationLimits, PairHistory, TripletAnnihilation,
};
use super::accounts::BeamTotals;
use super::*;

impl<const B: usize> SaturablePumpBeams<B> {
    /// Couple one optional pair law per UNIQUE original cell to all incident
    /// beams. Original photon order and per-beam ledgers are retained. Optical
    /// and nonlinear accuracy allowances both constrain the accepted history;
    /// neither weakens the native count/energy audit.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_annihilation(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: &[[f64; B]],
        duration_s: f64,
        laws: &[Option<TripletAnnihilation>],
        limits: PumpBeamsLimits,
        annihilation: AnnihilationLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PumpBeamsStep<B>, SingletTripletError> {
        self.advance_with_local_pumps_and_annihilation(
            populations,
            incident_photon_rates_s,
            &[],
            duration_s,
            laws,
            limits,
            annihilation,
            cancelled,
        )
    }

    /// Local photon fields, depleted beams and optional pair encounters share
    /// every original finite-population stage. Local fields are external supply
    /// only through their actual ground-state residence and absorbed photons.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_local_pumps_and_annihilation(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        duration_s: f64,
        laws: &[Option<TripletAnnihilation>],
        limits: PumpBeamsLimits,
        mut annihilation: AnnihilationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PumpBeamsStep<B>, SingletTripletError> {
        let local_pumps = if local::validate(&self.models, local_pumps, &mut cancelled)? {
            local_pumps
        } else {
            &[]
        };
        let active = pairs::validate(
            &self.models,
            laws,
            populations,
            duration_s,
            annihilation,
            &mut cancelled,
        )?;
        if !active || duration_s == 0.0 {
            return self.advance_with_local_pumps(
                populations,
                incident_photon_rates_s,
                local_pumps,
                duration_s,
                limits,
                cancelled,
            );
        }
        self.validate_advance(
            populations,
            incident_photon_rates_s,
            duration_s,
            limits,
            &mut cancelled,
        )?;
        let column = limits.column;
        if local_pumps.is_empty()
            && incident_photon_rates_s
                .iter()
                .flatten()
                .all(|&rate| rate == 0.0)
            && self.models.iter().copied().all(pairs::is_stable)
        {
            let mut step = self.advance_with_local_pumps(
                populations,
                incident_photon_rates_s,
                local_pumps,
                duration_s,
                limits,
                &mut cancelled,
            )?;
            for (index, cell) in step.cells.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                pairs::apply_isolated(
                    self.models[index],
                    populations[index],
                    laws[index],
                    duration_s,
                    cell,
                    column.kinetic.relative_tolerance,
                )?;
            }
            self.finish_pair_history(
                populations,
                incident_photon_rates_s,
                duration_s,
                column.kinetic.relative_tolerance,
                &mut step,
                &mut cancelled,
            )?;
            return Ok(step);
        }
        annihilation.maximum_substeps = annihilation.maximum_substeps.min(column.maximum_substeps);
        let nonlinear = self
            .paths
            .iter()
            .zip(incident_photon_rates_s)
            .any(|(path, rates)| {
                path.segments.iter().any(|segment| {
                    (0..B).any(|lane| {
                        rates[lane] > 0.0
                            && segment.optical_depth_per_ground[lane]
                                .iter()
                                .any(|&depth| depth > 0.0)
                    })
                })
            });
        let mut upper_hazards = vec![Sum::default(); self.models.len()];
        for (index, pump) in local_pumps.iter().enumerate() {
            for rate in pump.rates_s {
                upper_hazards[index].add(rate)?;
            }
        }
        for (path, rates) in self.paths.iter().zip(incident_photon_rates_s) {
            for segment in &path.segments {
                poll(&mut cancelled)?;
                for lane in 0..B {
                    for target in 0..2 {
                        upper_hazards[segment.cell_index].add(positive_product(
                            rates[lane],
                            segment.optical_depth_per_ground[lane][target],
                            "upper simultaneous pump hazard",
                        )?)?;
                    }
                }
            }
        }
        let mut fastest = pairs::fastest_rate(&self.models);
        for rate in upper_hazards {
            fastest = fastest.max(rate.value()?);
        }
        // Geometry is copied ONCE per call, not once per predictor/corrector.
        // Only this local model vector changes during nonlinear stages.
        let mut frozen = self.clone();
        let mut evaluations = 0usize;
        let mut coupling_iterations = 0usize;
        let result = pairs::midpoint(
            populations,
            laws,
            duration_s,
            fastest,
            annihilation,
            Some((
                column.relative_tolerance,
                column.absolute_energy_tolerance_j,
            )),
            |state, rates, duration, cancelled| {
                for (index, model) in frozen.models.iter_mut().enumerate() {
                    poll(cancelled)?;
                    *model = pairs::frozen_model(self.models[index], laws[index], rates[index])?;
                }
                let Some(mut step) = frozen.uniform_with_local(
                    state,
                    incident_photon_rates_s,
                    local_pumps,
                    duration,
                    1,
                    nonlinear,
                    limits,
                    &mut evaluations,
                    &mut coupling_iterations,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                for (index, cell) in step.cells.iter_mut().enumerate() {
                    pairs::restore(
                        self.models[index],
                        frozen.models[index],
                        state[index],
                        laws[index],
                        rates[index],
                        cell,
                        column.kinetic.relative_tolerance,
                    )?;
                }
                Ok(Some(step))
            },
            |step, cancelled| {
                self.finish_pair_history(
                    populations,
                    incident_photon_rates_s,
                    duration_s,
                    column.kinetic.relative_tolerance,
                    step,
                    cancelled,
                )
            },
            |step| self.components(step, duration_s),
            &mut cancelled,
        )?;
        let mut step = result.step;
        step.substeps = result.substeps;
        step.cell_evaluations = evaluations;
        step.coupling_iterations = coupling_iterations;
        step.max_energy_discrepancy_j = result.discrepancy_j;
        Ok(step)
    }

    fn finish_pair_history(
        &self,
        initial: &[[f64; 3]],
        incident: &[[f64; B]],
        duration: f64,
        tolerance: f64,
        step: &mut PumpBeamsStep<B>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        for (index, cell) in step.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            pairs::audit(self.models[index], initial[index], cell, tolerance)?;
        }
        for (index, beam) in step.beams.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut total = BeamTotals::new();
            total.add_summary(beam)?;
            *beam = total.finish(
                incident[index],
                duration,
                self.paths[index].photon_energies_j,
                tolerance,
            )?;
        }
        step.residual_j = self.energy_balance(
            initial,
            &step.cells,
            &step.passive_heat_j,
            &step.beams,
            &step.local_pumps,
            tolerance,
        )?;
        Ok(())
    }
}

impl<const B: usize> PairHistory for PumpBeamsStep<B> {
    fn populations(&self) -> Vec<[f64; 3]> {
        self.cells.iter().map(|cell| cell.populations).collect()
    }

    fn visit_integrated(
        &mut self,
        visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
    ) -> Result<(), SingletTripletError> {
        for cell in &mut self.cells {
            pairs::visit_kinetics(cell, visitor)?;
        }
        for value in &mut self.passive_heat_j {
            visitor(value)?;
        }
        for beam in &mut self.beams {
            for row in [
                &mut beam.incident_photons,
                &mut beam.transmitted_photons,
                &mut beam.absorbed_photons,
                &mut beam.passive_absorbed_photons,
                &mut beam.incident_energy_j,
                &mut beam.transmitted_energy_j,
                &mut beam.absorbed_energy_j,
                &mut beam.passive_absorbed_energy_j,
            ] {
                for value in row {
                    visitor(value)?;
                }
            }
        }
        for local in &mut self.local_pumps {
            for row in [
                &mut local.events,
                &mut local.absorbed_energy_j,
                &mut local.relaxation_heat_j,
            ] {
                for value in row {
                    visitor(value)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "annihilation/tests.rs"]
mod tests;
