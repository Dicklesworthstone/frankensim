//! A driven exterior channel in the SAME rational terminal owner.
//! There is no second DG propagator, rational integrator or delayed feedback.

use super::{AcousticRealizeError, DgPortLoad, DgPortResponse, TerminalCoupling, invalid};
use crate::bernoulli_aperture::radiation::environment::multiport::bake::mixed::DrivenMouthSurface;
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError;
use fs_exec::CancelGate;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};

mod independent;
pub use independent::{DgIndependentMouthExchange, DgIndependentMouthFrame,
    DgIndependentMouthResponse, DgIndependentMouthsFrame, DgIndependentMouthsResponse};

/// Accepted exterior mouth pressure and signed work from the original joint
/// model. This is energy exchange, not a heat or rational-storage certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgDrivenMouthFrame {
    /// Outward into the exterior: exactly the negative accepted DG inlet flow.
    pub outward_flow_m3_s: f64,
    pub mean_pressure_pa: f64,
    /// Original joint MIMO mean pressure times outward flow and duration [J].
    pub load_work_j: f64,
    /// Actual bore inlet mean pressure minus exterior mouth mean pressure [Pa].
    pub source_pressure_pa: f64,
    /// Original DG inlet work plus original exterior-mouth load work [J].
    pub source_work_j: f64,
    /// Independently evaluated source work minus its affine mean effort work.
    pub pressure_work_residual_j: f64,
}

impl DgDrivenMouthFrame {
    pub(super) fn new(
        inlet_flow: f64, mouth_pressure: f64, inlet_pressure: f64, inlet_work: f64,
        source: &DgPortResponse, dt: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let outward_flow_m3_s = -inlet_flow;
        let load_work_j = mouth_pressure * (outward_flow_m3_s * dt);
        let source_pressure_pa = inlet_pressure - mouth_pressure;
        let source_work_j = inlet_work + load_work_j;
        let expected = source.history_pressure_pa + source.impedance_pa_s_m3 * inlet_flow;
        let expected_work = expected * (inlet_flow * dt);
        let pressure_work_residual_j = source_work_j - expected_work;
        let pressure_scale = source.history_pressure_pa.abs()
            + (source.impedance_pa_s_m3 * inlet_flow).abs()
            + inlet_pressure.abs() + mouth_pressure.abs();
        let work_scale = inlet_work.abs() + load_work_j.abs() + expected_work.abs()
            + pressure_scale * inlet_flow.abs() * dt;
        if ![outward_flow_m3_s, load_work_j, source_pressure_pa, source_work_j,
            pressure_work_residual_j, pressure_scale, work_scale].iter().all(|value| value.is_finite())
            || (source_pressure_pa - expected).abs() > 1e-10 * pressure_scale.max(1e-300)
            || pressure_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
        {
            return Err(invalid("DG inlet and original joint mouth load disagree on source effort or work"));
        }
        Ok(Self { outward_flow_m3_s, mean_pressure_pa: mouth_pressure, load_work_j,
            source_pressure_pa, source_work_j, pressure_work_residual_j })
    }
}

impl DgPortLoad {
    /// Whether this load has any leading exterior-mouth coordinates. Scalar
    /// nonlinear adapters must explicitly support the singleton pressure jump;
    /// a larger group requires the complete independent-flow interface.
    pub fn has_driven_mouth(&self) -> bool {
        self.terminal_coupling.as_ref().is_some_and(|coupling| coupling.driven_mouth)
    }

    pub const fn last_driven_mouth_frame(&self) -> Option<&DgDrivenMouthFrame> {
        self.last_driven_mouth_frame.as_ref()
    }

    fn admit_cold_mouth(&self) -> Result<(), AcousticRealizeError> {
        if self.state.step_count() != 0
            || self.branches.iter().any(|branch| branch.state.step_count() != 0)
            || self.terminal.as_ref().is_none_or(|load|
                load.pole_states().iter().any(|state| state.re != 0.0 || state.im != 0.0))
            || self.terminal_coupling.as_ref().is_none_or(|coupling| coupling.driven_mouth)
        {
            return Err(invalid("joint mouth requires one cold DG/branch/terminal source"));
        }
        Ok(())
    }

    // Private to the shared owner and its tests. Public nonlinear adapters do
    // not silently acquire a different source-effort meaning.
    fn install_driven_mouth(&mut self, terminal: PreparedMimoLoad)
        -> Result<(), AcousticRealizeError>
    {
        self.install_independent_mouths(terminal, 1)
    }

    fn install_independent_mouths(&mut self, terminal: PreparedMimoLoad, mouths_count: usize)
        -> Result<(), AcousticRealizeError>
    {
        self.admit_cold_mouth()?;
        if mouths_count == 0
            || self.terminal_ports().len().checked_add(mouths_count) != Some(terminal.ports())
            || terminal.time_step_s().to_bits() != self.propagator.dt().to_bits()
            || terminal.pole_states().iter().any(|state| state.re != 0.0 || state.im != 0.0)
        {
            return Err(invalid("joint mouth requires cold rational memory on the original DG clock"));
        }
        let response = self.multiport_response()?;
        let coupling = TerminalCoupling::with_driven_mouth(&response, self.inlet,
            self.terminal_ports().to_vec(), &terminal, &self.terminal_channels)?;
        let impedance = coupling.relation(response, &terminal, &self.terminal_channels)?.inlet.impedance_pa_s_m3;
        // All construction and physical checks precede publication. Original
        // bore and branch owners/states remain exactly where they were.
        self.terminal = Some(terminal);
        self.terminal_coupling = Some(coupling);
        self.impedance_pa_s_m3 = impedance;
        Ok(())
    }

    /// Cold, source-specific replacement used by the flue's pressure-jump loop.
    /// Existing bore/branch operators and histories are preserved. The complete
    /// current physical source record replaces its former circular-only image.
    pub(crate) fn with_driven_mouth(
        self, mouth: DrivenMouthSurface, mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        self.with_driven_mouths(&[mouth], mechanics_rate_hz, gate)
    }

    pub(crate) fn with_driven_mouths(
        mut self, mouths: &[DrivenMouthSurface], mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        self.admit_cold_mouth().map_err(MimoRadiationInputError::Physics)?;
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let record = self.radiation_record.as_ref().ok_or_else(||
            MimoRadiationInputError::Physics(invalid("joint mouth requires the original computed exterior source")))?;
        let (model, record) = record.with_driven_mouths(mouths, mechanics_rate_hz, gate)?;
        let limits = record.limits();
        let terminal = PreparedMimoLoad::new(&model, self.propagator.dt(), MimoLoadLimits {
            band_rad_s: [core::f64::consts::TAU * record.band_hz()[0],
                core::f64::consts::PI / self.propagator.dt()],
            sampled_passivity_points: limits.sampled_passivity_points,
            max_states: limits.max_states, max_memory_bytes: limits.max_memory_bytes,
        }).map_err(|error| MimoRadiationInputError::Physics(invalid(error.0)))?;
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        if mouths.len() == 1 {
            self.install_driven_mouth(terminal)
        } else {
            self.install_independent_mouths(terminal, mouths.len())
        }.map_err(MimoRadiationInputError::Physics)?;
        self.radiation_record = Some(record);
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        Ok(self)
    }
}

#[cfg(test)]
mod tests;

