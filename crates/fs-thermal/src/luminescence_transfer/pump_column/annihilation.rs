//! Joint pair kinetics and directed pump depletion on one accepted history.

use super::super::saturation::annihilation::{
    self as pairs, AnnihilationLimits, PairHistory, TripletAnnihilation,
};
use super::*;

impl<const B: usize> SaturablePumpColumn<B> {
    /// Resolve local triplet-pair encounters together with changing pump
    /// opacity. Laws are ordered by original column cell. Absent/zero laws use
    /// the unchanged column path. Predictors, correctors and rejected optical
    /// trials share the original aggregate cell-evaluation budget.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_annihilation(
        &self,
        populations: &[[f64; 3]],
        incident_photon_rates_s: [f64; B],
        duration_s: f64,
        laws: &[Option<TripletAnnihilation>],
        limits: PumpColumnLimits,
        mut annihilation: AnnihilationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PumpColumnStep<B>, SingletTripletError> {
        let models: Vec<_> = self.cells.iter().map(|cell| cell.model).collect();
        let active = pairs::validate(
            &models,
            laws,
            populations,
            duration_s,
            annihilation,
            &mut cancelled,
        )?;
        if !active || duration_s == 0.0 {
            return self.advance(
                populations,
                incident_photon_rates_s,
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
        if incident_photon_rates_s.iter().all(|&rate| rate == 0.0)
            && models.iter().copied().all(pairs::is_stable)
        {
            let mut step = self.advance(
                populations,
                incident_photon_rates_s,
                duration_s,
                limits,
                &mut cancelled,
            )?;
            for (index, cell) in step.cells.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                pairs::apply_isolated(
                    models[index],
                    populations[index],
                    laws[index],
                    duration_s,
                    &mut cell.kinetics,
                    limits.kinetic.relative_tolerance,
                )?;
            }
            self.finish_pair_history(
                populations,
                incident_photon_rates_s,
                duration_s,
                limits,
                &mut step,
                &mut cancelled,
            )?;
            return Ok(step);
        }
        annihilation.maximum_substeps = annihilation.maximum_substeps.min(limits.maximum_substeps);
        let mut fastest = pairs::fastest_rate(&models);
        for cell in &self.cells {
            let mut ground_rate = Sum::default();
            for lane in 0..B {
                for target in 0..2 {
                    ground_rate.add(positive_product(
                        incident_photon_rates_s[lane],
                        cell.optical_depth_per_ground[lane][target],
                        "upper pump hazard",
                    )?)?;
                }
            }
            fastest = fastest.max(ground_rate.value()?);
        }
        let mut frozen = self.clone();
        let mut evaluations = 0usize;
        let result = pairs::midpoint(
            populations,
            laws,
            duration_s,
            fastest,
            annihilation,
            Some((
                limits.relative_tolerance,
                limits.absolute_energy_tolerance_j,
            )),
            |state, rates, duration, cancelled| {
                for (index, cell) in frozen.cells.iter_mut().enumerate() {
                    poll(cancelled)?;
                    cell.model = pairs::frozen_model(models[index], laws[index], rates[index])?;
                }
                // One residence/transport stage. The outer physical schedule
                // controls BOTH density and opacity time error; no nested
                // independent pair kick or per-stage renewed work budget.
                let mut step = frozen.uniform(
                    state,
                    incident_photon_rates_s,
                    duration,
                    1,
                    limits,
                    &mut evaluations,
                    cancelled,
                )?;
                for (index, cell) in step.cells.iter_mut().enumerate() {
                    pairs::restore(
                        models[index],
                        frozen.cells[index].model,
                        state[index],
                        laws[index],
                        rates[index],
                        &mut cell.kinetics,
                        limits.kinetic.relative_tolerance,
                    )?;
                }
                Ok(Some(step))
            },
            |step, cancelled| {
                self.finish_pair_history(
                    populations,
                    incident_photon_rates_s,
                    duration_s,
                    limits,
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
        step.max_energy_discrepancy_j = result.discrepancy_j;
        Ok(step)
    }

    fn finish_pair_history(
        &self,
        initial: &[[f64; 3]],
        incident: [f64; B],
        duration: f64,
        limits: PumpColumnLimits,
        step: &mut PumpColumnStep<B>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        let mut totals = Vec::with_capacity(step.cells.len());
        let mut populations = Vec::with_capacity(step.cells.len());
        for cell in &step.cells {
            poll(cancelled)?;
            let mut total = CellTotals::new();
            total.add(cell)?;
            totals.push(total);
            populations.push(cell.kinetics.populations);
        }
        let mut transmitted = [Sum::default(); B];
        add_values(&mut transmitted, step.transmitted_photons)?;
        *step = self.finish_history(
            initial,
            &populations,
            totals,
            transmitted,
            incident,
            duration,
            step.substeps,
            limits,
            step.cell_evaluations,
            cancelled,
        )?;
        Ok(())
    }
}

impl<const B: usize> PairHistory for PumpColumnStep<B> {
    fn populations(&self) -> Vec<[f64; 3]> {
        self.cells
            .iter()
            .map(|cell| cell.kinetics.populations)
            .collect()
    }

    fn visit_integrated(
        &mut self,
        visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
    ) -> Result<(), SingletTripletError> {
        for cell in &mut self.cells {
            pairs::visit_kinetics(&mut cell.kinetics, visitor)?;
            for row in [
                &mut cell.incident_photons,
                &mut cell.transmitted_photons,
                &mut cell.absorbed_photons,
                &mut cell.passive_absorbed_photons,
                &mut cell.incident_energy_j,
                &mut cell.transmitted_energy_j,
                &mut cell.absorbed_energy_j,
                &mut cell.passive_absorbed_energy_j,
            ] {
                for value in row {
                    visitor(value)?;
                }
            }
            visitor(&mut cell.passive_heat_j)?;
        }
        for row in [
            &mut self.incident_photons,
            &mut self.transmitted_photons,
            &mut self.absorbed_photons,
            &mut self.passive_absorbed_photons,
            &mut self.incident_energy_j,
            &mut self.transmitted_energy_j,
            &mut self.absorbed_energy_j,
            &mut self.passive_absorbed_energy_j,
        ] {
            for value in row {
                visitor(value)?;
            }
        }
        Ok(())
    }
}
