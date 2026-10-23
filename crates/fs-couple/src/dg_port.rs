//! Source-bound exact mean-port loads shared by nonlinear acoustic couplers.
//!
//! The original DG owner retains geometry, material gas, wall memory, storage
//! and propagation. This module only admits one cold inlet and makes a complete
//! held-flow step transactional. It never reconstructs or fits a second bore.

/// Authored geometry, material gas and numerical declarations for this load.
pub mod input;

use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgPortResponse, DgPropagator, DgState, DgStepOutput, PortId, PortInput,
};

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// One declared inlet, one original propagator and one accepted physical history.
/// Read its exact mean impedance before constructing a nonlinear source.
#[derive(Debug, Clone)]
pub struct DgPortLoad {
    duct: DgDuct,
    propagator: DgPropagator,
    state: DgState,
    inlet: PortId,
    impedance_pa_s_m3: f64,
}

impl DgPortLoad {
    /// Prepare an initially quiescent duct with one inlet port and a closed,
    /// ideal-open, or matched outlet. The DOF allowance covers acoustic AND wall
    /// states before construction of the original dense exponential.
    ///
    /// The complete material wall band must lie below the mechanical Nyquist
    /// frequency. Geometry, material, pole and exponential admissions remain
    /// owned by the supplied duct. A zero return-period impedance is refused.
    ///
    /// # Errors
    /// Missing inlet, a second unknown port, state-budget excess, unresolved
    /// exponential, inadmissible wall band or nonpositive mean-port impedance.
    pub fn new(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(inlet) = duct.inlet() else {
            return Err(invalid("DG load requires one declared inlet port"));
        };
        if matches!(duct.outlet(), DgBoundary::Port(_)) {
            return Err(invalid("DG load cannot independently solve a second unknown port"));
        }
        if duct.total_dofs() > maximum_dofs {
            return Err(invalid("DG acoustic plus memory states exceed the declared DOF budget"));
        }
        if let Some(band) = duct.viscothermal_band_rad_s() {
            if band[1] * time_step_s > core::f64::consts::PI {
                return Err(invalid("DG viscothermal band exceeds the mechanical Nyquist frequency"));
            }
        }
        let propagator = duct
            .propagator(time_step_s, exponential_tolerance)
            .map_err(AcousticRealizeError::Duct)?;
        let state = duct.zero_state();
        let impedance_pa_s_m3 = propagator
            .port_response(&state, inlet)
            .map_err(AcousticRealizeError::Duct)?
            .impedance_pa_s_m3;
        if impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("DG load requires strictly positive exact mean-port impedance"));
        }
        Ok(Self { duct, propagator, state, inlet, impedance_pa_s_m3 })
    }

    /// Exact-ZOH mean pressure per held inlet volume flow [Pa s/m³].
    #[must_use]
    pub const fn impedance_pa_s_m3(&self) -> f64 { self.impedance_pa_s_m3 }

    /// Original geometry and material source.
    #[must_use]
    pub const fn duct(&self) -> &DgDuct { &self.duct }

    /// Original exponential, admitted clock and certificate.
    #[must_use]
    pub const fn propagator(&self) -> &DgPropagator { &self.propagator }

    /// Complete accepted state, including current physical nodal fields.
    #[must_use]
    pub const fn state(&self) -> &DgState { &self.state }

    /// The unique inlet port receiving positive flow into the duct.
    #[must_use]
    pub const fn inlet(&self) -> PortId { self.inlet }

    /// Affine step-mean pressure from the original accepted physical history.
    ///
    /// # Errors
    /// The original propagator's source-binding or finite-response refusal.
    pub fn port_response(&self) -> Result<DgPortResponse, AcousticRealizeError> {
        self.propagator.port_response(&self.state, self.inlet)
            .map_err(AcousticRealizeError::Duct)
    }

    /// Retain a private candidate until the nonlinear participant has validated
    /// its own pressure, work and state. The exclusive borrow prevents a trial
    /// from being accepted into another source or an already advanced history.
    pub(crate) fn preview_step(
        &mut self, volume_flow_m3_s: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        let mut state = self.state.clone();
        let output = self.propagator.step(
            &mut state, &[PortInput { id: self.inlet, effort: volume_flow_m3_s }],
        ).map_err(AcousticRealizeError::Duct)?;
        self.propagator.sync_nodal(&mut state).map_err(AcousticRealizeError::Duct)?;
        let ledger = &output.ledger;
        if output.port_outputs.len() != 1
            || !output.port_outputs[0].conjugate_rate.is_finite()
            || ![
                state.stored_energy_j(), ledger.stored_energy_j,
                ledger.port_work_step_j, ledger.cumulative_port_work_j,
                ledger.cumulative_dissipated_j, ledger.energy_balance_residual_j,
                ledger.relative_energy_drift,
            ].iter().all(|value| value.is_finite())
            || state.pressure().iter().chain(state.volume_flow()).any(|x| !x.is_finite())
        {
            return Err(invalid("DG port candidate left the finite pressure/state/work domain"));
        }
        Ok(DgPortTrial { owner: self, state, output })
    }
}

/// A trial borrows the only accepted history. Dropping it publishes nothing.
pub(crate) struct DgPortTrial<'a> {
    owner: &'a mut DgPortLoad,
    state: DgState,
    output: DgStepOutput,
}

impl DgPortTrial<'_> {
    pub(crate) const fn previous_state(&self) -> &DgState { &self.owner.state }
    pub(crate) const fn state(&self) -> &DgState { &self.state }
    pub(crate) const fn output(&self) -> &DgStepOutput { &self.output }

    /// Publish the already validated candidate; all fallible checks precede this.
    pub(crate) fn accept(self) -> DgStepOutput {
        self.owner.state = self.state;
        self.output
    }
}

#[cfg(test)]
mod tests;
