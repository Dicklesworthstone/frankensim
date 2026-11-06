//! Local nuclides and once-per-cohort finite-material movement comparisons.

use super::*;

impl AdaptiveNuclearHeatSource for BurnupTransportCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        let index = self.compartment_index();
        for candidate in [coarse, fine] {
            if candidate.compartment_index() != index
                || !self.shares_forcing(candidate)
                || self.state().volumes_m3() != candidate.state().volumes_m3()
            {
                return Err(FissionThermalError::Invalid(
                    "refinement changed movement coordinates or prescribed histories",
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
        let index = coarse.compartment_index();
        if first.compartment_index() != index || second.compartment_index() != index {
            return Err(FissionThermalError::Invalid(
                "refinement changed movement report ordering",
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
        metric.reaction_counts(
            &coarse.transport().produced_nuclei[index],
            &first.transport().produced_nuclei[index],
            &second.transport().produced_nuclei[index],
        )?;
        if index == 0 {
            // These extensive expected-nuclei counts belong to the complete
            // cohort. Compare every directed path and isotope exactly once,
            // rather than allowing equal total transfers to hide redistribution.
            for (a, b, c) in [
                (
                    &coarse.transport().transfer_counts,
                    &first.transport().transfer_counts,
                    &second.transport().transfer_counts,
                ),
                (
                    &coarse.transport().outflow_counts,
                    &first.transport().outflow_counts,
                    &second.transport().outflow_counts,
                ),
            ] {
                if a.len() != b.len() || a.len() != c.len() {
                    return Err(FissionThermalError::Invalid(
                        "refinement changed movement path ordering",
                    ));
                }
                for ((a, b), c) in a.iter().zip(b).zip(c) {
                    metric.reaction_counts(a, b, c)?;
                }
            }
        }
        Ok(metric.error)
    }
}
