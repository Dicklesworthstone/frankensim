//! Physical external accounts and compensated original-history accumulation.

use super::*;

impl<const B: usize> PairHistory for History<B> {
    fn populations(&self) -> Vec<[f64; 3]> {
        self.beams.populations()
    }

    fn visit_integrated(
        &mut self,
        visitor: &mut impl FnMut(&mut f64) -> Result<(), SingletTripletError>,
    ) -> Result<(), SingletTripletError> {
        self.beams.visit_integrated(visitor)?;
        for rows in [
            &mut self.external_absorbed,
            &mut self.recycled_absorbed,
            &mut self.external_events,
            &mut self.recycled_events,
            &mut self.external_relaxation,
            &mut self.recycled_relaxation,
        ] {
            for row in rows {
                for value in row {
                    visitor(value)?;
                }
            }
        }
        for row in [&mut self.external_passive, &mut self.recycled_passive] {
            for value in row {
                visitor(value)?;
            }
        }
        Ok(())
    }
}

impl<const B: usize> SaturablePhotonRecycling<B> {
    pub(super) fn external_scale(
        &self,
        initial: &[[f64; 3]],
        history: &PumpBeamsStep<B>,
    ) -> Result<f64, SingletTripletError> {
        let mut scale = Sum::default();
        for (index, model) in self.models().iter().enumerate() {
            for state in 0..2 {
                scale.add(positive_product(
                    initial[index][state + 1],
                    model.excitation_energies_j()[state],
                    "initial recycling store",
                )?)?;
            }
        }
        for beam in &history.beams[..self.external_paths] {
            for value in beam.incident_energy_j {
                scale.add(value)?;
            }
        }
        for local in &history.local_pumps {
            for energy in local.absorbed_energy_j {
                scale.add(energy)?;
            }
        }
        // Deliberately omit gross photon circulation and temporary frozen-pair
        // reaction coordinates. Neither is external physical energy supply.
        scale.value()
    }

    pub(super) fn audit_history(
        &self,
        initial: &[[f64; 3]],
        history: &mut History<B>,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        for (index, cell) in history.beams.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            pairs::audit(self.models()[index], initial[index], cell, tolerance)?;
            for state in 0..2 {
                let mut residual = Sum::default();
                residual.add(history.external_absorbed[index][state])?;
                residual.add(history.recycled_absorbed[index][state])?;
                if let Some(local) = history.beams.local_pumps.get(index) {
                    residual.add(local.absorbed_energy_j[state])?;
                }
                residual.add(-cell.absorbed_pump_energy_j[state])?;
                admit_energy(
                    residual.value()?,
                    cell.absorbed_pump_energy_j[state],
                    tolerance,
                )?;
            }
            let mut passive = Sum::default();
            passive.add(history.external_passive[index])?;
            passive.add(history.recycled_passive[index])?;
            passive.add(-history.beams.passive_heat_j[index])?;
            admit_energy(
                passive.value()?,
                history.beams.passive_heat_j[index],
                tolerance,
            )?;
        }
        for (beam, path) in history.beams.beams.iter_mut().zip(&self.transport.paths) {
            poll(cancelled)?;
            beam.photon_residual = photon_balance(
                beam.incident_photons,
                beam.transmitted_photons,
                beam.absorbed_photons,
                beam.passive_absorbed_photons,
            )?;
            admit_photons(
                beam.incident_photons,
                beam.transmitted_photons,
                beam.absorbed_photons,
                beam.passive_absorbed_photons,
                beam.photon_residual,
                tolerance,
            )?;
            beam.incident_energy_j = energies(beam.incident_photons, path.photon_energies_j)?;
            beam.transmitted_energy_j = energies(beam.transmitted_photons, path.photon_energies_j)?;
            beam.absorbed_energy_j = energies(beam.absorbed_photons, path.photon_energies_j)?;
            beam.passive_absorbed_energy_j =
                energies(beam.passive_absorbed_photons, path.photon_energies_j)?;
        }
        Ok(())
    }

    pub(super) fn components(
        &self,
        history: &History<B>,
        duration: f64,
    ) -> Result<Vec<f64>, SingletTripletError> {
        let mut components = self.transport.components(&history.beams, duration)?;
        for rows in [
            &history.external_absorbed,
            &history.recycled_absorbed,
            &history.external_relaxation,
            &history.recycled_relaxation,
        ] {
            for row in rows {
                components.extend(row);
            }
        }
        for rows in [&history.external_events, &history.recycled_events] {
            for (index, row) in rows.iter().enumerate() {
                for state in 0..2 {
                    components.push(positive_product(
                        row[state],
                        self.models()[index].excitation_energies_j()[state],
                        "attributed event convergence energy",
                    )?);
                }
            }
        }
        components.extend(&history.external_passive);
        components.extend(&history.recycled_passive);
        Ok(components)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        &self,
        initial: &[[f64; 3]],
        mut history: History<B>,
        substeps: usize,
        cell_evaluations: usize,
        coupling_iterations: usize,
        feedback_iterations: usize,
        max_energy_discrepancy_j: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SaturableRecyclingStep<B>, SingletTripletError> {
        self.audit_history(initial, &mut history, tolerance, cancelled)?;
        let cells = self.models().len();
        let mut source_heat = vec![[Sum::default(); 2]; cells];
        let mut emitted = vec![[Sum::default(); 2]; cells];
        let mut escaped = vec![[Sum::default(); 2]; cells];
        let mut unresolved = vec![[Sum::default(); 2]; cells];
        let mut source_discrepancy = Sum::default();
        for (path, photons) in self
            .emission
            .iter()
            .zip(&history.beams.beams[self.external_paths..])
        {
            poll(cancelled)?;
            let c = path.source_cell;
            let s = path.source_state;
            let state_energy = self.models()[c].excitation_energies_j()[s];
            for lane in 0..B {
                let count = positive_product(
                    history.beams.cells[c].radiative_events[s],
                    path.photon_fractions[lane],
                    "emitted spectral photons",
                )?;
                let photon_energy = path.path.photon_energies_j[lane];
                emitted[c][s].add(positive_product(
                    count,
                    photon_energy,
                    "emitted spectral energy",
                )?)?;
                if path.photon_fractions[lane] > 0.0 {
                    source_heat[c][s].add(positive_product(
                        count,
                        state_energy - photon_energy,
                        "source Stokes heat",
                    )?)?;
                }
                source_discrepancy.add(positive_product(
                    (count - photons.incident_photons[lane]).abs(),
                    photon_energy,
                    "integrated emitted source mismatch",
                )?)?;
                let target = match path.terminal {
                    EmissionTerminal::Escape => &mut escaped[c][s],
                    EmissionTerminal::Unresolved => &mut unresolved[c][s],
                };
                target.add(photons.transmitted_energy_j[lane])?;
            }
        }
        admit_energy(
            source_discrepancy.value()?,
            self.external_scale(initial, &history.beams)?,
            tolerance,
        )?;
        let source_relaxation_heat_j = source_heat
            .into_iter()
            .map(values)
            .collect::<Result<Vec<_>, _>>()?;
        let emitted_energy_j = emitted
            .into_iter()
            .map(values)
            .collect::<Result<Vec<_>, _>>()?;
        let escaped_energy_j = escaped
            .into_iter()
            .map(values)
            .collect::<Result<Vec<_>, _>>()?;
        let unresolved_energy_j = unresolved
            .into_iter()
            .map(values)
            .collect::<Result<Vec<_>, _>>()?;
        let mut residual = Sum::default();
        let mut scale = Sum::default();
        for local in &history.beams.local_pumps {
            for energy in local.absorbed_energy_j {
                residual.add(energy)?;
                scale.add(energy)?;
            }
        }
        for beam in &history.beams.beams[..self.external_paths] {
            for lane in 0..B {
                residual.add(beam.incident_energy_j[lane])?;
                residual.add(-beam.transmitted_energy_j[lane])?;
                scale.add(beam.incident_energy_j[lane])?;
                scale.add(beam.transmitted_energy_j[lane])?;
            }
        }
        for (index, cell) in history.beams.cells.iter().enumerate() {
            poll(cancelled)?;
            for state in 0..2 {
                let old = positive_product(
                    initial[index][state + 1],
                    self.models()[index].excitation_energies_j()[state],
                    "initial recycling excitation",
                )?;
                residual.add(old)?;
                scale.add(old)?;
                for value in [
                    cell.stored_energy_j[state],
                    cell.nonradiative_heat_j[state],
                    cell.pump_relaxation_heat_j[state],
                    source_relaxation_heat_j[index][state],
                    escaped_energy_j[index][state],
                    unresolved_energy_j[index][state],
                ] {
                    residual.add(-value)?;
                    scale.add(value)?;
                }
            }
            residual.add(cell.reverse_crossing_heat_absorbed_j)?;
            scale.add(cell.reverse_crossing_heat_absorbed_j)?;
            for value in [
                cell.crossing_heat_j,
                cell.annihilation_heat_j,
                history.beams.passive_heat_j[index],
            ] {
                residual.add(-value)?;
                scale.add(value)?;
            }
        }
        let residual_j = residual.value()?;
        admit_energy(residual_j, scale.value()?, tolerance)?;
        let emission_paths = history.beams.beams.split_off(self.external_paths);
        poll(cancelled)?;
        Ok(SaturableRecyclingStep {
            cells: history.beams.cells,
            passive_heat_j: history.beams.passive_heat_j,
            external_beams: history.beams.beams,
            emission_paths,
            local_pumps: history.beams.local_pumps,
            external_absorbed_pump_energy_j: history.external_absorbed,
            recycled_absorbed_pump_energy_j: history.recycled_absorbed,
            external_pump_events: history.external_events,
            recycled_pump_events: history.recycled_events,
            external_pump_relaxation_heat_j: history.external_relaxation,
            recycled_pump_relaxation_heat_j: history.recycled_relaxation,
            external_passive_heat_j: history.external_passive,
            recycled_passive_heat_j: history.recycled_passive,
            source_relaxation_heat_j,
            emitted_energy_j,
            escaped_energy_j,
            unresolved_energy_j,
            residual_j,
            substeps,
            cell_evaluations,
            coupling_iterations,
            feedback_iterations,
            max_energy_discrepancy_j,
        })
    }
}
