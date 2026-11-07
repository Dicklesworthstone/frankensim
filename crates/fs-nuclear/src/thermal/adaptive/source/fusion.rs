//! Fusion coordinates on the existing nuclear/material time-refinement driver.
//!
//! Compare finite populations and current-interval releases, including separate
//! kinetic and nuclear channels. Large lifetime exports cannot hide a small
//! late heat signal. Shared neutron/nuclide transfers are visited only once.

use super::*;
use crate::diffusion::DiffusionState;
use crate::fusion_neutrons::NeutronKineticModel;

fn same_cell<S: NuclearHeatSource>(
    initial: &S,
    coarse: &S,
    fine: &S,
    indices: [usize; 3],
    forcing: [&FusionNeutronForcing; 3],
) -> Result<usize, FissionThermalError> {
    same_binding(initial, coarse, fine)?;
    if indices[0] != indices[1]
        || indices[0] != indices[2]
        || !forcing[0].same_history(forcing[1])
        || !forcing[0].same_history(forcing[2])
    {
        return Err(FissionThermalError::Invalid(
            "fusion refinement requires original cell ordering and source history",
        ));
    }
    Ok(indices[0])
}

fn neutron_state(
    metric: &mut Comparison<'_>,
    cell: usize,
    models: [&Arc<NeutronKineticModel>; 3],
    states: [&DiffusionState; 3],
) -> Result<(), FissionThermalError> {
    same_model(models[0].as_ref(), models[1].as_ref(), models[2].as_ref())?;
    let groups = models[0].group_energies_j().len();
    let offset = cell * groups;
    let mut energy = [Sum::default(); 3];
    for group in 0..groups {
        let counts = states.map(|state| state.neutron_counts()[offset + group]);
        metric.neutron_count(counts[0], counts[1], counts[2])?;
        for (total, count) in energy.iter_mut().zip(counts) {
            total.add(finite(count * models[0].group_energies_j()[group])?)?;
        }
    }
    metric.energy(energy[0].value()?, energy[1].value()?, energy[2].value()?)
}

fn energy_steps(
    metric: &mut Comparison<'_>,
    a: EnergyChannels,
    b: EnergyChannels,
    c: EnergyChannels,
) -> Result<(), FissionThermalError> {
    for ((a, b), c) in channels(a).into_iter().zip(channels(b)).zip(channels(c)) {
        metric.energy(0.0, a, finite(b + c)?)?;
    }
    Ok(())
}

impl AdaptiveNuclearHeatSource for NeutronKineticCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let cell = same_cell(
            self,
            coarse,
            fine,
            [self.cell_index(), coarse.cell_index(), fine.cell_index()],
            [self.forcing(), coarse.forcing(), fine.forcing()],
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        neutron_state(
            &mut metric,
            cell,
            [
                self.state().model(),
                coarse.state().model(),
                fine.state().model(),
            ],
            [
                self.state().neutrons(),
                coarse.state().neutrons(),
                fine.state().neutrons(),
            ],
        )?;
        Ok(metric.error)
    }

    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        if coarse.cell_index() != first.cell_index()
            || coarse.cell_index() != second.cell_index()
            || !coarse.forcing().same_history(first.forcing())
            || !coarse.forcing().same_history(second.forcing())
        {
            return Err(FissionThermalError::Invalid(
                "fusion refinement changed interval binding",
            ));
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        let cell = coarse.cell_index();
        metric.energy(
            0.0,
            coarse.spatial().cell_nuclear_reaction_energy_j[cell],
            finite(
                first.spatial().cell_nuclear_reaction_energy_j[cell]
                    + second.spatial().cell_nuclear_reaction_energy_j[cell],
            )?,
        )?;
        if coarse.cell_index() == 0 {
            diffusion::neutron_steps(
                &mut metric,
                &coarse.spatial().diffusion,
                &first.spatial().diffusion,
                &second.spatial().diffusion,
            )?;
        }
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for FusionNeutronActivationCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let cell = same_cell(
            self,
            coarse,
            fine,
            [self.cell_index(), coarse.cell_index(), fine.cell_index()],
            [self.forcing(), coarse.forcing(), fine.forcing()],
        )?;
        for candidate in [coarse, fine] {
            match (
                self.state().spatial().transport(),
                candidate.state().spatial().transport(),
            ) {
                (None, None) => {}
                (Some(a), Some(b))
                    if a.volumes_m3() == b.volumes_m3() && a.transport() == b.transport() => {}
                _ => {
                    return Err(FissionThermalError::Invalid(
                        "fusion refinement changed material transport",
                    ));
                }
            }
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        neutron_state(
            &mut metric,
            cell,
            [
                self.state().model(),
                coarse.state().model(),
                fine.state().model(),
            ],
            [
                self.state().neutrons(),
                coarse.state().neutrons(),
                fine.state().neutrons(),
            ],
        )?;
        burnup(
            &mut metric,
            &self.state().inventories()[cell],
            &coarse.state().inventories()[cell],
            &fine.state().inventories()[cell],
        )?;
        Ok(metric.error)
    }

    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let cell = coarse.cell_index();
        if cell != first.cell_index()
            || cell != second.cell_index()
            || !coarse.forcing().same_history(first.forcing())
            || !coarse.forcing().same_history(second.forcing())
        {
            return Err(FissionThermalError::Invalid(
                "fusion refinement changed activation binding",
            ));
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        burnup_steps(
            &mut metric,
            coarse.burnup(),
            first.burnup(),
            second.burnup(),
        )?;
        let [a, b, c] = [coarse.activation(), first.activation(), second.activation()];
        // Opposing changes in kinetic deposition and nuclear Q must not cancel
        // in the combined heat channel and falsely imply temporal convergence.
        energy_steps(
            &mut metric,
            a.kinetic.cell_released[cell],
            b.kinetic.cell_released[cell],
            c.kinetic.cell_released[cell],
        )?;
        metric.energy(
            0.0,
            a.kinetic.cell_nuclear_reaction_energy_j[cell],
            finite(
                b.kinetic.cell_nuclear_reaction_energy_j[cell]
                    + c.kinetic.cell_nuclear_reaction_energy_j[cell],
            )?,
        )?;
        energy_steps(
            &mut metric,
            a.spatial.burnup[cell].released,
            b.spatial.burnup[cell].released,
            c.spatial.burnup[cell].released,
        )?;
        metric.reaction_counts(
            &a.spatial.produced_nuclei[cell],
            &b.spatial.produced_nuclei[cell],
            &c.spatial.produced_nuclei[cell],
        )?;
        if cell == 0 {
            diffusion::neutron_steps(
                &mut metric,
                &a.spatial.neutrons,
                &b.spatial.neutrons,
                &c.spatial.neutrons,
            )?;
            match (
                &a.spatial.transport,
                &b.spatial.transport,
                &c.spatial.transport,
            ) {
                (None, None, None) => {}
                (Some(a), Some(b), Some(c)) => {
                    for (a, b, c) in [
                        (&a.transfer_counts, &b.transfer_counts, &c.transfer_counts),
                        (&a.outflow_counts, &b.outflow_counts, &c.outflow_counts),
                    ] {
                        if a.len() != b.len() || a.len() != c.len() {
                            return Err(FissionThermalError::Invalid(
                                "fusion refinement changed material transfer reports",
                            ));
                        }
                        for ((a, b), c) in a.iter().zip(b).zip(c) {
                            metric.reaction_counts(a, b, c)?;
                        }
                    }
                }
                _ => {
                    return Err(FissionThermalError::Invalid(
                        "fusion refinement changed material transport reports",
                    ));
                }
            }
        }
        Ok(metric.error)
    }
}
