//! Cell-local coordinates and once-per-cohort spatial reaction/transfer accounts.
use super::*;

impl AdaptiveNuclearHeatSource for DiffusionCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        let index = self.cell_index();
        if coarse.cell_index() != index || fine.cell_index() != index {
            return Err(FissionThermalError::Invalid(
                "refinement changed spatial cell ordering",
            ));
        }
        let initial = self.state().neutrons();
        let a = coarse.state().neutrons();
        let b = fine.state().neutrons();
        let model = initial.model();
        let groups = model.groups().len();
        let delayed = model.delayed_groups().len();
        for candidate in [coarse, fine] {
            match (self.state().transport(), candidate.state().transport()) {
                (None, None) => {}
                (Some(initial), Some(candidate))
                    if initial.volumes_m3() == candidate.volumes_m3()
                        && initial.transport() == candidate.transport() => {}
                _ => {
                    return Err(FissionThermalError::Invalid(
                        "refinement changed spatial material transport",
                    ));
                }
            }
        }
        for state in [a, b] {
            let candidate = state.model();
            // Burnup deliberately rebinds macroscopic reaction coefficients.
            // Compare the retained coordinate layout, not that operator's Arc.
            if !Arc::ptr_eq(model.groups(), candidate.groups())
                || model.cell_count() != candidate.cell_count()
                || model.cell_volumes_m3()[index] != candidate.cell_volumes_m3()[index]
                || model.speeds_m_s() != candidate.speeds_m_s()
                || model.delayed_groups() != candidate.delayed_groups()
            {
                return Err(FissionThermalError::Invalid(
                    "refinement changed spatial neutron coordinates",
                ));
            }
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        burnup(
            &mut metric,
            &self.state().inventories()[index],
            &coarse.state().inventories()[index],
            &fine.state().inventories()[index],
        )?;
        for group in 0..groups {
            let at = index * groups + group;
            metric.neutron_count(
                initial.neutron_counts()[at],
                a.neutron_counts()[at],
                b.neutron_counts()[at],
            )?;
        }
        for group in 0..delayed {
            let at = index * delayed + group;
            metric.neutron_count(
                initial.precursor_counts()[at],
                a.precursor_counts()[at],
                b.precursor_counts()[at],
            )?;
        }
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
        if coarse.cell_index() != first.cell_index() || coarse.cell_index() != second.cell_index() {
            return Err(FissionThermalError::Invalid(
                "refinement changed spatial report ordering",
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
        let cell = coarse.cell_index();
        metric.reaction_counts(
            &coarse.spatial().produced_nuclei[cell],
            &first.spatial().produced_nuclei[cell],
            &second.spatial().produced_nuclei[cell],
        )?;
        if coarse.cell_index() == 0 {
            match (
                &coarse.spatial().transport,
                &first.spatial().transport,
                &second.spatial().transport,
            ) {
                (None, None, None) => {}
                (Some(a), Some(b), Some(c)) => {
                    // Nuclide movement is a separate extensive account from
                    // neutron/precursor movement and nuclear reaction counts.
                    for (a, b, c) in [
                        (&a.transfer_counts, &b.transfer_counts, &c.transfer_counts),
                        (&a.outflow_counts, &b.outflow_counts, &c.outflow_counts),
                    ] {
                        if a.len() != b.len() || a.len() != c.len() {
                            return Err(FissionThermalError::Invalid(
                                "refinement changed spatial material path ordering",
                            ));
                        }
                        for ((a, b), c) in a.iter().zip(b).zip(c) {
                            metric.reaction_counts(a, b, c)?;
                        }
                    }
                }
                _ => {
                    return Err(FissionThermalError::Invalid(
                        "refinement changed spatial material transport reports",
                    ));
                }
            }
            // Shared transport is solved once and its full edge/group accounts
            // are compared once, not once for every thermal cell handle.
            let a = &coarse.spatial().neutrons;
            let b = &first.spatial().neutrons;
            let c = &second.spatial().neutrons;
            neutron_steps(&mut metric, a, b, c)?;
        }
        Ok(metric.error)
    }
}

/// Original extensive neutron accounts, compared once per complete cohort.
pub(super) fn neutron_steps(
    metric: &mut Comparison<'_>,
    a: &crate::diffusion::DiffusionStep,
    b: &crate::diffusion::DiffusionStep,
    c: &crate::diffusion::DiffusionStep,
) -> Result<(), FissionThermalError> {
    for (a, b, c) in [
        (&a.reaction_counts, &b.reaction_counts, &c.reaction_counts),
        (
            &a.link_a_to_b_counts,
            &b.link_a_to_b_counts,
            &c.link_a_to_b_counts,
        ),
        (
            &a.link_b_to_a_counts,
            &b.link_b_to_a_counts,
            &c.link_b_to_a_counts,
        ),
        (&a.leakage_counts, &b.leakage_counts, &c.leakage_counts),
        (&a.delayed_decays, &b.delayed_decays, &c.delayed_decays),
        (
            &a.injected_neutron_counts,
            &b.injected_neutron_counts,
            &c.injected_neutron_counts,
        ),
        (
            &a.injected_precursor_counts,
            &b.injected_precursor_counts,
            &c.injected_precursor_counts,
        ),
        (
            &a.precursor_transfer_counts,
            &b.precursor_transfer_counts,
            &c.precursor_transfer_counts,
        ),
        (
            &a.precursor_outflow_counts,
            &b.precursor_outflow_counts,
            &c.precursor_outflow_counts,
        ),
    ] {
        metric.particle_counts(a, b, c)?;
    }
    Ok(())
}
