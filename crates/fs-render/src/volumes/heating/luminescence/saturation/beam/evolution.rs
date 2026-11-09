//! One shared physical state, one accepted thermal timeline.

use super::*;

impl DepletingBeam {
    pub(super) fn prepare<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedPumps>, OpticalHeatingError> {
        if let Some(prepared) = &self.prepared {
            let mut matches = true;
            for (&index, expected) in self.model.cells.iter().zip(&prepared.temperatures) {
                checkpoint()?;
                if thermal
                    .temperature_k(owner.emitters[index].binding.thermal_cell)?
                    .to_bits()
                    != *expected
                {
                    matches = false;
                    break;
                }
            }
            if matches {
                return Ok(Arc::clone(prepared));
            }
        }
        let mut temperatures = Vec::with_capacity(self.model.cells.len());
        let mut models = Vec::with_capacity(self.model.cells.len());
        let mut cache = None;
        for &index in &self.model.cells {
            checkpoint()?;
            let temperature = thermal.temperature_k(owner.emitters[index].binding.thermal_cell)?;
            temperatures.push(temperature.to_bits());
            models.push(resolve_model(
                &owner.singlet,
                &owner.triplet,
                temperature,
                &mut cache,
            )?);
        }
        let mut paths = Vec::with_capacity(self.model.paths.len());
        let mut segment_count = 0usize;
        for path in &self.model.paths {
            let mut segments = Vec::with_capacity(path.cells.len());
            for cell in &path.cells {
                checkpoint()?;
                let sections = cell.cross_sections.evaluate(
                    path.wavelengths_nm,
                    f64::from_bits(temperatures[cell.cell_index]),
                )?;
                let mut depth = [[0.0; 2]; 2];
                for lane in 0..2 {
                    for target in 0..2 {
                        let source = if target == 0 {
                            &owner.singlet
                        } else {
                            owner.triplet.triplet()
                        };
                        let (lo, hi) = source.fluorescence().excitation_bounds_nm();
                        if sections[lane][target] > 0.0
                            && !(lo..=hi).contains(&path.wavelengths_nm[lane])
                        {
                            return Err(invalid(
                                "molecular pump exceeds original excitation support",
                            ));
                        }
                        depth[lane][target] = super::super::super::rate_power(
                            sections[lane][target],
                            cell.length_m,
                            cell.volume_m3,
                        )?;
                        if sections[lane][target] > 0.0 && depth[lane][target] == 0.0 {
                            return Err(invalid(
                                "positive molecular pump optical depth is below representation",
                            ));
                        }
                    }
                }
                segments.push(SaturablePumpSegment {
                    cell_index: cell.cell_index,
                    optical_depth_per_ground: depth,
                    passive_optical_depth: cell.passive_optical_depth,
                });
            }
            segment_count = segment_count
                .checked_add(segments.len())
                .ok_or_else(|| invalid("pump segment count overflow"))?;
            paths.push(SaturablePumpPath {
                photon_energies_j: path.photon_energies_j,
                segments,
            });
        }
        let solver = if paths.len() == 1 {
            let path = &paths[0];
            let cells: Vec<_> = path
                .segments
                .iter()
                .map(|segment| SaturablePumpCell {
                    model: models[segment.cell_index],
                    optical_depth_per_ground: segment.optical_depth_per_ground,
                    passive_optical_depth: segment.passive_optical_depth,
                })
                .collect();
            PreparedSolver::Column(preserve_checkpoint(checkpoint, |cancelled| {
                SaturablePumpColumn::try_new(
                    &cells,
                    path.photon_energies_j,
                    cells.len().max(1),
                    cancelled,
                )
            })?)
        } else {
            PreparedSolver::Beams(preserve_checkpoint(checkpoint, |cancelled| {
                SaturablePumpBeams::try_new(
                    &models,
                    &paths,
                    models.len().max(1),
                    paths.len(),
                    segment_count.max(1),
                    cancelled,
                )
            })?)
        };
        Ok(Arc::new(PreparedPumps {
            temperatures,
            solver,
            paths,
        }))
    }

    /// Remap the already prepared optical paths to the complete original
    /// emitter inventory, including cells outside every external pump.
    pub(in crate::volumes::heating::luminescence::saturation) fn recycling_paths<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<SaturablePumpPath<2>>, OpticalHeatingError> {
        let prepared = self.prepare(owner, thermal, checkpoint)?;
        let mut paths = prepared.paths.clone();
        for path in &mut paths {
            for segment in &mut path.segments {
                checkpoint()?;
                segment.cell_index = self.model.cells[segment.cell_index];
            }
        }
        Ok(paths)
    }

    pub(in crate::volumes::heating::luminescence::saturation) fn incident_rates(
        &self,
        time_s: f64,
        globally_pumped: bool,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<[f64; 2]>, OpticalHeatingError> {
        self.model
            .paths
            .iter()
            .map(|path| {
                checkpoint()?;
                Ok(
                    if globally_pumped
                        && path
                            .pump_window_s
                            .is_none_or(|(start, end)| time_s >= start && time_s < end)
                    {
                        path.incident_photon_rates_s
                    } else {
                        [0.0; 2]
                    },
                )
            })
            .collect()
    }

    pub(in crate::volumes::heating::luminescence::saturation) fn recycling_accounts(
        &self,
        step: &fs_thermal::luminescence::transfer::pump_beams::recycling::SaturableRecyclingStep<2>,
        duration_s: f64,
        tolerance_j: f64,
    ) -> Result<BeamAccounts, OpticalHeatingError> {
        let mut accounts = self.accounts.clone();
        accounts.add_paths(
            &step.external_beams,
            duration_s,
            step.substeps,
            step.cell_evaluations,
            step.coupling_iterations,
            step.max_energy_discrepancy_j,
            tolerance_j,
        )?;
        Ok(accounts)
    }

    pub(in crate::volumes::heating::luminescence::saturation) fn check_state<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedPumps>, OpticalHeatingError> {
        self.prepare(owner, thermal, checkpoint)
    }

    pub(in crate::volumes::heating::luminescence::saturation) fn next_edge(
        &self,
        after_s: f64,
    ) -> Option<f64> {
        self.model
            .paths
            .iter()
            .filter_map(|path| path.pump_window_s)
            .flat_map(|(start, end)| [start, end])
            .filter(|&time| time > after_s)
            .min_by(f64::total_cmp)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::volumes::heating::luminescence::saturation) fn advance<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        duration_s: f64,
        globally_pumped: bool,
        annihilation: Option<&super::super::annihilation::PreparedAnnihilation>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BeamCandidate, OpticalHeatingError> {
        let prepared = self.prepare(owner, thermal, checkpoint)?;
        let rates = self.incident_rates(owner.time_s, globally_pumped, checkpoint)?;
        let mut accounts = self.accounts.clone();
        let mut roundoff = CheckedSum::default();
        let mut steps = vec![None; owner.emitters.len()];
        let mut passive_heat_j = Vec::with_capacity(self.model.cells.len());
        match &prepared.solver {
            PreparedSolver::Column(column) => {
                let path = &self.model.paths[0];
                let populations: Vec<_> = path
                    .cells
                    .iter()
                    .map(|cell| owner.emitters[self.model.cells[cell.cell_index]].populations)
                    .collect();
                let step = if let Some(pairs) = annihilation {
                    let laws: Vec<_> = path
                        .cells
                        .iter()
                        .map(|cell| pairs.laws[self.model.cells[cell.cell_index]])
                        .collect();
                    preserve_checkpoint(checkpoint, |cancelled| {
                        column.advance_with_annihilation(
                            &populations,
                            rates[0],
                            duration_s,
                            &laws,
                            self.model.limits.column,
                            pairs.limits,
                            cancelled,
                        )
                    })?
                } else {
                    preserve_checkpoint(checkpoint, |cancelled| {
                        column.advance(
                            &populations,
                            rates[0],
                            duration_s,
                            self.model.limits.column,
                            cancelled,
                        )
                    })?
                };
                accounts.add_column(&step, duration_s, owner.energy_tolerance_j)?;
                roundoff.add(step.residual_j.abs())?;
                for lane in 0..2 {
                    roundoff.add(finite(
                        step.photon_residual[lane].abs() * path.photon_energies_j[lane],
                    )?)?;
                }
                for (cell, result) in path.cells.iter().zip(step.cells) {
                    checkpoint()?;
                    let index = self.model.cells[cell.cell_index];
                    passive_heat_j.push((
                        owner.emitters[index].binding.thermal_cell,
                        result.passive_heat_j,
                    ));
                    steps[index] = Some(result.kinetics);
                }
            }
            PreparedSolver::Beams(beams) => {
                let populations: Vec<_> = self
                    .model
                    .cells
                    .iter()
                    .map(|&index| owner.emitters[index].populations)
                    .collect();
                let step = if let Some(pairs) = annihilation {
                    let laws: Vec<_> = self
                        .model
                        .cells
                        .iter()
                        .map(|&index| pairs.laws[index])
                        .collect();
                    preserve_checkpoint(checkpoint, |cancelled| {
                        beams.advance_with_annihilation(
                            &populations,
                            &rates,
                            duration_s,
                            &laws,
                            self.model.limits,
                            pairs.limits,
                            cancelled,
                        )
                    })?
                } else {
                    preserve_checkpoint(checkpoint, |cancelled| {
                        beams.advance(
                            &populations,
                            &rates,
                            duration_s,
                            self.model.limits,
                            cancelled,
                        )
                    })?
                };
                accounts.add_beams(&step, duration_s, owner.energy_tolerance_j)?;
                roundoff.add(step.residual_j.abs())?;
                for (path, result) in self.model.paths.iter().zip(&step.beams) {
                    for lane in 0..2 {
                        roundoff.add(finite(
                            result.photon_residual[lane].abs() * path.photon_energies_j[lane],
                        )?)?;
                    }
                }
                if step.cells.len() != self.model.cells.len()
                    || step.passive_heat_j.len() != self.model.cells.len()
                {
                    return Err(invalid(
                        "coupled pump result differs from original cell inventory",
                    ));
                }
                for ((&index, kinetics), heat) in self
                    .model
                    .cells
                    .iter()
                    .zip(step.cells)
                    .zip(step.passive_heat_j)
                {
                    checkpoint()?;
                    passive_heat_j.push((owner.emitters[index].binding.thermal_cell, heat));
                    steps[index] = Some(kinetics);
                }
            }
        }
        checkpoint()?;
        Ok(BeamCandidate {
            steps,
            passive_heat_j,
            absolute_roundoff_j: roundoff.value()?,
            accounts,
        })
    }

    pub(in crate::volumes::heating::luminescence::saturation) fn accept(
        &mut self,
        accounts: BeamAccounts,
        prepared: Arc<PreparedPumps>,
    ) {
        self.accounts = accounts;
        self.prepared = Some(prepared);
    }
}
