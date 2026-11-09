//! Multi-terminal source controls share the original electrothermal residual.

use super::*;
use fs_em::conduction::response::terminals::ConductionTerminalDrive;

#[derive(Clone, Copy)]
pub(in crate::circuit_enthalpy::conduction) enum Drive<'a> {
    Field(ConductionDrive<'a>),
    Terminals(ConductionTerminalDrive<'a>),
}

impl Drive<'_> {
    pub(in crate::circuit_enthalpy::conduction) fn supplied_power_w(
        self,
        field: &ConductionResponseReport,
        terminals: Option<&ConductionTerminalReport>,
    ) -> Result<f64, Error> {
        match self {
            Self::Field(drive) => drive.supplied_power_w(
                &field.voltage_v,
                field.power_w,
                field.power_balance_residual_w,
            ),
            Self::Terminals(_) => terminals
                .map(|report| report.supplied_power_w)
                .ok_or_else(|| invalid("nonlinear terminal operating point is absent")),
        }
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance ideal voltage/current/mixed electrodes and original nonlinear
    /// bulk/contact sources with the same finite enthalpy endpoint. Floating
    /// electrodes retain zero net current while transferring current locally.
    /// Each terminal field closes actual E/J and declared interface coordinates
    /// after refreshing the original thermal/phase state.
    ///
    /// Supply one original contact heat binding per physical interface, or an
    /// empty slice when there are none. The existing source identity, caloric,
    /// Fourier and geometric admission applies. Every native field solve counts
    /// against maximum_circuit_trials, including discarded constitutive iterates
    /// and final fresh admission. Source work comes from authored current
    /// controls and native voltage reactions at retained physical drops.
    ///
    /// Accepted reports retain the ordered terminal operating point. Both
    /// electrical and thermal owners continue together; refusal or cancellation
    /// leaves both supplied owners reusable. This adds no electrical storage,
    /// field-dependent caloric/Fourier law or uniqueness/convergence claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_terminal_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionTerminalDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_response_driven_resolved_candidate(
            thermal,
            target_s,
            problem,
            Drive::Terminals(drive),
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }
}
