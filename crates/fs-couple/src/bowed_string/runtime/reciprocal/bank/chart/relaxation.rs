//! Material history reacts through the original shared body and physical bridges.
use super::*;
use crate::plate_viscosity::relaxation::{
    InitialPlateMemory, ModalPlateRelaxation, PlateBendingRelaxation, PlateRelaxationLimits,
};

impl MaterialStringBankSource {
    /// All original regional spectra and their full same-basis strain factors.
    /// This image survives gas/torque rebuilding and multibow/PCM ownership.
    #[must_use]
    pub const fn material_chart_relaxation(&self) -> Option<&ModalPlateRelaxation> {
        self.material_chart_relaxation.as_ref()
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_body_memory_frequencies(
        &self, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        poll(gate)?;
        if let Some(image) = &self.material_chart_relaxation {
            // A coupled body receives forcing and relaxation below its bare
            // resonances as well as at the instantaneous coupled frequencies.
            // Validate all source coefficients, even explicitly zero arms.
            image.source().validate_band(0.0, system.linear_frequency_upper_bound_rad_s(),
                1.0 / f64::from(self.options.sample_rate_hz))
                .map_err(ReciprocalStringError::Material)?;
        }
        Ok(())
    }

    // Reconstruct only the shared owner's coordinate embedding after a cold
    // topology change. Body modes/material factors are NEVER resolved again.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn restore_body_memory(
        &self, mut system: FrictionModalSystem, memory: Option<Vec<f64>>, gate: &CancelGate,
    ) -> Result<FrictionModalSystem, ReciprocalStringError> {
        poll(gate)?;
        if let Some(image) = &self.material_chart_relaxation {
            let body = self.strings.len();
            let component = system.components().get(body)
                .ok_or_else(|| invalid("regional memory lost its original body component"))?;
            if component.modes().len() != image.mode_count()
                || component.modes().iter().any(|m| m.damping_ratio != 0.0) {
                return Err(invalid("regional body memory requires the unchanged complete elastic body basis"));
            }
            let start: usize = system.components()[..body].iter().map(|c| c.modes().len()).sum();
            let dimension = system.mode_count().checked_mul(2)
                .ok_or_else(|| invalid("body memory state embedding overflow"))?;
            let branches = image.coupling_rows_per_s().iter().zip(image.relaxation_times_s())
                .map(|(row, &time)| {
                    let mut projection = vec![0.0; dimension];
                    for (i, &value) in row.iter().enumerate() { projection[2*(start+i)] = value; }
                    fs_phs::RelaxationBranch { projection, stiffness: 1.0, relaxation_time_s: time }
                }).collect();
            system = system.with_material_relaxation(branches, memory, gate)
                .map_err(ReciprocalStringError::Friction)?;
            self.validate_body_memory_frequencies(&system, gate)?;
            self.validate_material_bridge_frequencies(&system, gate)?;
        }
        poll(gate)?;
        Ok(system)
    }
}

impl ReciprocalStringBank {
    /// Resolve independent regional material spectra into the ORIGINAL body's
    /// retained DKT basis and attach their internal strains to the SAME coupled
    /// string/body/cavity/friction owner. Material changes therefore modify both
    /// physical pressure and the bridge reactions back onto every string.
    ///
    /// Each region's selected excess modulus scales its actual rotated bending
    /// tensor; all modal cross terms survive. The equilibrium elastic pencil,
    /// conservative prestress, nodal bridge footprints, gas and observer remain
    /// in that one basis. This is not a general tensor-Prony plane-stress law.
    ///
    /// The existing dense exact capacity and setup allowance include every real
    /// first-order memory state. Prepared mechanical-only propagation refuses.
    /// Source validity conservatively covers DC through the complete coupled
    /// instantaneous stiffness bound. This does not qualify finite-amplitude
    /// harmonic content or reproduce a measured violin/cello/bass by itself.
    ///
    /// # Errors
    /// Requires one cold attachment to an original chart with no selected scalar
    /// or tensor viscosity or authored body decay. Missing source properties,
    /// unresolved positive forms, excessive state/projection work, source/clock
    /// refusal and cancellation return no partially assembled performer.
    pub fn with_material_plate_relaxation(
        mut self, properties: &[(String, String)], initial: InitialPlateMemory,
        limits: PlateRelaxationLimits, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if self.system.samples_rendered() != 0 || self.source.material_chart_relaxation.is_some()
            || self.source.material_chart_viscosity.is_some()
            || !self.system.material_relaxation_branches().is_empty() {
            return Err(invalid("body material memory requires one cold attachment without parallel intrinsic damping"));
        }
        let body = self.source.material_chart.as_ref()
            .ok_or_else(|| invalid("body material memory requires an original regional plate chart"))?;
        let basis = self.source.material_chart_basis.as_ref()
            .ok_or_else(|| invalid("body material memory requires the original full DKT basis"))?;
        if body.reduction.damping_ratio != 0.0 || body.specimen.bending_viscosity().is_some() {
            return Err(invalid("body material memory cannot duplicate authored or selected intrinsic body loss"));
        }
        let law = PlateBendingRelaxation::from_chart(&body.specimen, properties, gate)
            .map_err(ReciprocalStringError::Material)?;
        let image = law.project(&body.specimen, basis, limits, gate)
            .map_err(ReciprocalStringError::Material)?;
        let memory = match initial {
            InitialPlateMemory::Relaxed => None,
            InitialPlateMemory::Unrelaxed => Some(vec![0.0; image.coupling_rows_per_s().len()]),
        };
        self.source.material_chart_relaxation = Some(image);
        self.system = self.source.restore_body_memory(self.system, memory, gate)?;
        Ok(self)
    }
}
