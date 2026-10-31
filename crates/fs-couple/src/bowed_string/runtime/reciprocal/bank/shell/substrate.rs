//! Explicitly selected substrate histories in the existing source-memory block.
//! They precede coating histories, preserving the coating-only heat observers.
use super::*;
use crate::shear_shell::constant_loss::relaxation::ShellSubstrateRelaxation;

impl MaterialStringBankSource {
    /// Complete causal substrate storage/loss on the original shell basis.
    /// This is selected on the shell specimen, then attached by its bank
    /// constructor. Cold gas/topology rebuilds retain the same physical rows.
    pub fn material_shell_substrate_relaxation(&self) -> Option<&ShellSubstrateRelaxation> {
        self.material_shell_image.as_ref().and_then(|image| image.substrate_relaxation())
            .or_else(|| self.material_shell_assembly_image.as_ref()
                .and_then(|image| image.substrate_relaxation()))
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn shell_substrate_memory_branches(
        &self, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<Vec<fs_phs::RelaxationBranch>, ReciprocalStringError> {
        poll(gate)?;
        let Some(image) = self.material_shell_substrate_relaxation() else {
            return Ok(Vec::new());
        };
        let body = self.strings.len();
        let component = system.components().get(body)
            .ok_or_else(|| invalid("substrate memory lost its original body component"))?;
        if component.modes().len() != image.mode_count()
            || component.modes().iter().any(|mode| mode.damping_ratio != 0.0) {
            return Err(invalid("substrate memory needs its unchanged complete equilibrium body basis"));
        }
        let start: usize = system.components()[..body].iter()
            .map(|component| component.modes().len()).sum();
        let dimension = system.mode_count().checked_mul(2)
            .ok_or_else(|| invalid("substrate memory embedding extent overflow"))?;
        if image.coupling_rows_per_s().len() != image.relaxation_times_s().len()
            || image.coupling_rows_per_s().len().checked_mul(dimension)
                .is_none_or(|work| work > self.options.network.max_setup_terms) {
            return Err(invalid("complete substrate histories exceed the original setup allowance"));
        }
        let mut branches = Vec::with_capacity(image.coupling_rows_per_s().len());
        for (row, &time) in image.coupling_rows_per_s().iter().zip(image.relaxation_times_s()) {
            poll(gate)?;
            if row.len() != image.mode_count() {
                return Err(invalid("substrate memory lost its complete modal projection"));
            }
            let mut projection = vec![0.0; dimension];
            for (mode, &value) in row.iter().enumerate() {
                projection[2*(start+mode)] = value;
            }
            branches.push(fs_phs::RelaxationBranch {
                projection, stiffness: 1.0, relaxation_time_s: time,
            });
        }
        Ok(branches)
    }
}

impl ReciprocalStringBank {
    // The specimen method already explicitly selected the passive continuation.
    // Construction must attach it atomically; it cannot publish an equilibrium
    // body whose frequency-dependent stiffness and loss are silently missing.
    pub(super) fn finish_substrate_relaxation(mut self, gate: &CancelGate)
        -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let Some(image) = self.source.material_shell_substrate_relaxation() else {
            return Ok(self);
        };
        let count = image.relaxation_times_s().len();
        let source_count = self.source.source_memory_count(&self.cards, gate)?;
        if self.system.samples_rendered() != 0
            || self.system.material_relaxation_branches().len().checked_add(count) != Some(source_count) {
            return Err(invalid("substrate source histories require one cold complete attachment"));
        }
        if count != 0 {
            let branches = self.source.shell_substrate_memory_branches(&self.system, gate)?;
            self.system = self.system.with_additional_material_relaxation(branches, None, gate)
                .map_err(ReciprocalStringError::Friction)?;
            self.source.validate_body_memory_frequencies(&self.system, gate)?;
            self.source.validate_material_bridge_frequencies(&self.system, gate)?;
        }
        poll(gate)?;
        Ok(self)
    }
}
