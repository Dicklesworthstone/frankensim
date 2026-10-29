//! Projected film histories in the original reciprocal string/body owner.
//! The source fit is finite-band. Selecting this time-domain image explicitly
//! accepts the positive Prony model's causal continuation outside that band;
//! it does not turn the continuation into measured material evidence.
use super::*;
use crate::shear_shell::coating::relaxation::ShellCoatingRelaxation;

impl MaterialStringBankSource {
    /// Bound film spectrum, original fit band and numerical fit certificates.
    /// Its rows use the same equilibrium body basis as wood loss, bridge work
    /// and acoustic volume. The source survives cold gas/torsion rebuilds.
    pub fn material_shell_coating_relaxation(&self) -> Option<&ShellCoatingRelaxation> {
        self.material_shell_coating_relaxation.as_ref()
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn shell_coating_memory_branches(
        &self,
        system: &FrictionModalSystem,
        gate: &CancelGate,
    ) -> Result<Vec<fs_phs::RelaxationBranch>, ReciprocalStringError> {
        poll(gate)?;
        let Some(image) = &self.material_shell_coating_relaxation else {
            return Ok(Vec::new());
        };
        let body = self.strings.len();
        let component = system
            .components()
            .get(body)
            .ok_or_else(|| invalid("coating memory lost its original shell component"))?;
        if component.modes().len() != image.mode_count()
            || component
                .modes()
                .iter()
                .any(|mode| mode.damping_ratio != 0.0)
        {
            return Err(invalid(
                "coating memory requires its unchanged complete equilibrium shell basis",
            ));
        }
        let start: usize = system.components()[..body]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let dimension = system
            .mode_count()
            .checked_mul(2)
            .ok_or_else(|| invalid("coating memory embedding extent overflow"))?;
        if image.coupling_rows_per_s().len() != image.relaxation_times_s().len()
            || image
                .coupling_rows_per_s()
                .len()
                .checked_mul(dimension)
                .is_none_or(|n| n > self.options.network.max_setup_terms)
        {
            return Err(invalid(
                "complete coating history embedding exceeds the source setup allowance",
            ));
        }
        let mut branches = Vec::with_capacity(image.coupling_rows_per_s().len());
        for (row, &time) in image
            .coupling_rows_per_s()
            .iter()
            .zip(image.relaxation_times_s())
        {
            poll(gate)?;
            if row.len() != image.mode_count() {
                return Err(invalid("coating memory lost its complete modal projection"));
            }
            let mut projection = vec![0.0; dimension];
            for (mode, &value) in row.iter().enumerate() {
                projection[2 * (start + mode)] = value;
            }
            branches.push(fs_phs::RelaxationBranch {
                projection,
                stiffness: 1.0,
                relaxation_time_s: time,
            });
        }
        Ok(branches)
    }
}

impl ReciprocalStringBank {
    /// Bind every projected coating branch to the SAME coupled mechanics,
    /// pressure acceleration and material heat account before playback.
    ///
    /// The shell must already have selected `with_coating_relaxation`: its
    /// pencil contains film equilibrium stiffness and its damping contains
    /// wood/interface viscosity only. The positive film branches add storage
    /// and loss once. Full modal cross terms survive the source projection.
    ///
    /// This explicitly selects a passive causal continuation outside each
    /// source's measured fit band. Retained shell modes and fit errors were
    /// admitted by the shell producer; the continuation does not claim source
    /// accuracy at DC, outside that band or for finite-amplitude excitation.
    /// Existing string/plate source-band guards remain in force. The complete
    /// instantaneous stiffness still must fit the original Nyquist and energy
    /// limits, and every history counts against the chosen propagation budget.
    ///
    /// `None` initializes ONLY the new histories relaxed at the current body
    /// displacement. `Some` supplies one energy-normalized value per row.
    /// Existing string histories remain an unchanged prefix; later explicit
    /// string and thermal histories follow this body block.
    pub fn with_material_shell_coating_relaxation(
        mut self,
        initial_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let source_count = self.source.source_memory_count(&self.cards, gate)?;
        if self.system.samples_rendered() != 0
            || self.source.material_shell_coating_relaxation.is_some()
            || self.source.material_chart_relaxation.is_some()
            || self.source.material_chart_tensor_relaxation.is_some()
            || self.source.string_bending_relaxation.is_some()
            || self.system.material_relaxation_branches().len() != source_count
        {
            return Err(invalid(
                "coating memory requires one cold binding before additional string or thermal histories",
            ));
        }
        let image = self
            .source
            .material_shell_image
            .as_ref()
            .and_then(|image| image.coating_relaxation())
            .or_else(|| {
                self.source
                    .material_shell_assembly_image
                    .as_ref()
                    .and_then(|image| image.coating_relaxation())
            })
            .ok_or_else(|| {
                invalid("coating memory requires the original equilibrium shell projection")
            })?
            .clone();
        self.source.material_shell_coating_relaxation = Some(image);
        let branches = self
            .source
            .shell_coating_memory_branches(&self.system, gate)?;
        self.system = self
            .system
            .with_additional_material_relaxation(branches, initial_memory, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source
            .validate_body_memory_frequencies(&self.system, gate)?;
        self.source
            .validate_material_bridge_frequencies(&self.system, gate)?;
        poll(gate)?;
        Ok(self)
    }
}
