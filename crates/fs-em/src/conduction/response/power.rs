//! Constant-power closure through the original nonlinear material owner.

use super::*;

impl PlanarConductionResponse {
    /// Close original bulk/contact laws under a held terminal power set point.
    /// Each constitutive trial has its own native resistance and actual E/J;
    /// the accepted nonlinear field is never rescaled after source closure.
    /// One native basis solve is charged per constitutive trial.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_power_candidate(
        &self,
        positive: &[usize],
        negative: &[usize],
        source: ConductionPowerSource,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_power_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            positive,
            negative,
            source,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same closure at an explicitly refreshed thermal/phase context.
    /// Supply compliance applies to the converged physical operating point,
    /// not an unconverged coefficient iterate. All original material domains,
    /// current/power residuals, contact closure, budgets and cancellation gates
    /// still apply to every trial. Unsupported or nonconvergent sources refuse;
    /// this bounded Picard path promises neither uniqueness nor convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_power_with_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        positive: &[usize],
        negative: &[usize],
        source: ConductionPowerSource,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        poll(cancelled)?;
        source
            .validate()
            .map_err(ConductionResponseError::Invalid)?;
        let trial_source = ConductionPowerSource {
            maximum_voltage_v: f64::MAX,
            maximum_current_a: f64::MAX,
            ..source
        };
        let (candidate, resistance_ohm) = self
            .solve_coupled_with_materials_and_contacts_candidate(
                seeds,
                contacts,
                assembly,
                field_controls,
                controls,
                cancelled,
                |native, cancelled| {
                    let response = native.resistance_at_reference_voltage(
                        positive,
                        negative,
                        1.0,
                        field_controls,
                        &mut || cancelled(),
                    )?;
                    let resistance_ohm = response.resistance_ohm();
                    let field = response.field_at_power(
                        trial_source,
                        field_controls.current_tolerance_a,
                        field_controls.power_tolerance_w,
                        &mut || cancelled(),
                    )?;
                    Ok::<_, ConductionResponseError>((field, resistance_ohm))
                },
            )?;
        source
            .operating_point(resistance_ohm)
            .map_err(ConductionResponseError::Invalid)?;
        poll(cancelled)?;
        Ok(candidate)
    }
}
