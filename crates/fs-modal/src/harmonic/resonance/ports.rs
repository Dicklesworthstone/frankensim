//! Actual complex load-port observations for resonance design.
//!
//! Port flow and effort are the original power-conjugate variables supplied by
//! the physical load owner. A coherent receiver is a transfer row, without
//! conjugating its weights. Observation changes no force, dissipation or work.

use super::super::{
    HarmonicFrequencyParameterResponse, HarmonicFrequencyResponse, HarmonicResponseError, invalid,
};
use super::{
    HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet, HarmonicTransferSample,
};
use fs_math::c64::C64;

impl HarmonicFrequencyResponse {
    /// Observe the original port flow through fixed real receiver weights.
    ///
    /// Flow is B^T*v, in the physical owner's original port order and units.
    /// The receiver may combine ports coherently with either sign.
    pub fn port_flow_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, true)?;
        HarmonicTransferSample::project_real(receiver, values, derivative)
    }

    /// Observe original port effort through fixed real receiver weights.
    ///
    /// Effort is Z*flow. Its frequency tangent includes the complete Z_omega
    /// action, so transmitted-force and motion resonances may differ.
    pub fn port_effort_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, false)?;
        HarmonicTransferSample::project_real(receiver, values, derivative)
    }

    /// Observe port flow through fixed complex coherent transfer weights.
    pub fn port_flow_transfer_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, true)?;
        HarmonicTransferSample::project_complex(receiver, values, derivative)
    }

    /// Observe port effort through fixed complex coherent transfer weights.
    pub fn port_effort_transfer_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, false)?;
        HarmonicTransferSample::project_complex(receiver, values, derivative)
    }

    /// Observe port flow through an authored frequency-dependent receiver.
    ///
    /// Retains both r_omega*flow and r*flow_omega. Receiver propagation or
    /// filter dispersion is separate from the original physical load law.
    pub fn port_flow_transfer_complex_with_derivative(
        &self,
        receiver: &[C64],
        receiver_angular_derivative: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, true)?;
        HarmonicTransferSample::project_complex_with_derivative(
            receiver,
            receiver_angular_derivative,
            values,
            derivative,
        )
    }

    /// Observe port effort through an authored frequency-dependent receiver.
    pub fn port_effort_transfer_complex_with_derivative(
        &self,
        receiver: &[C64],
        receiver_angular_derivative: &[C64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        let (values, derivative) = frequency_port_fields(self, false)?;
        HarmonicTransferSample::project_complex_with_derivative(
            receiver,
            receiver_angular_derivative,
            values,
            derivative,
        )
    }
}

impl HarmonicFrequencyParameterResponse {
    /// Borrow the complete physical port-flow jet without another solve.
    ///
    /// The parameter and mixed terms include both the motion change and the
    /// authored attachment change: u_theta=B^T*v_theta+B_theta^T*v and
    /// u_omega_theta=B^T*v_omega_theta+B_theta^T*v_omega.
    /// Use the shared receiver-jet projection if receiver weights also change.
    pub fn port_flow_field_parameter_jet(
        &self,
    ) -> Result<HarmonicTransferFieldParameterJet<'_>, HarmonicResponseError> {
        parameter_port_fields(self, true)
    }

    /// Borrow the complete physical effort jet without another solve.
    ///
    /// This includes all impedance, flow, frequency and parameter product
    /// terms. It does not infer a receiver or a force-map derivative from the
    /// published effort. Those observation derivatives remain caller-owned.
    pub fn port_effort_field_parameter_jet(
        &self,
    ) -> Result<HarmonicTransferFieldParameterJet<'_>, HarmonicResponseError> {
        parameter_port_fields(self, false)
    }

    /// Observe the port-flow parameter jet through a fixed real receiver.
    pub fn port_flow_transfer_parameter(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_real(receiver, self.port_flow_field_parameter_jet()?)
    }

    /// Observe the port-effort parameter jet through a fixed real receiver.
    pub fn port_effort_transfer_parameter(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_real(
            receiver,
            self.port_effort_field_parameter_jet()?,
        )
    }

    /// Observe the port-flow parameter jet through fixed complex transfer weights.
    pub fn port_flow_transfer_parameter_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_complex(
            receiver,
            self.port_flow_field_parameter_jet()?,
        )
    }

    /// Observe the port-effort parameter jet through fixed complex transfer weights.
    pub fn port_effort_transfer_parameter_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_complex(
            receiver,
            self.port_effort_field_parameter_jet()?,
        )
    }
}

fn frequency_port_fields(
    response: &HarmonicFrequencyResponse,
    flow: bool,
) -> Result<(&[C64], &[C64]), HarmonicResponseError> {
    let primal = response.response.port_load.as_ref().ok_or_else(|| {
        invalid("a physical port load is required for a port transfer observation")
    })?;
    let first = response
        .port_load
        .as_ref()
        .ok_or_else(|| invalid("complete physical port frequency derivatives are required"))?;
    Ok(if flow {
        (&primal.flow, &first.flow_derivative)
    } else {
        (&primal.effort, &first.effort_derivative)
    })
}

fn parameter_port_fields(
    response: &HarmonicFrequencyParameterResponse,
    flow: bool,
) -> Result<HarmonicTransferFieldParameterJet<'_>, HarmonicResponseError> {
    let frequency = &response.frequency_response.frequency_response;
    let (values, angular_derivative) = frequency_port_fields(frequency, flow)?;
    let curvature = response
        .frequency_response
        .port_load
        .as_ref()
        .ok_or_else(|| invalid("complete physical port frequency curvatures are required"))?;
    let parameter = response
        .parameter_response
        .port_load
        .as_ref()
        .ok_or_else(|| invalid("complete physical port parameter derivatives are required"))?;
    let mixed = response
        .port_load
        .as_ref()
        .ok_or_else(|| invalid("complete physical port mixed derivatives are required"))?;
    Ok(if flow {
        HarmonicTransferFieldParameterJet {
            values,
            angular_derivative,
            angular_second_derivative: &curvature.flow_second_derivative,
            parameter_derivative: &parameter.flow_derivative,
            angular_parameter_derivative: &mixed.flow_mixed_derivative,
        }
    } else {
        HarmonicTransferFieldParameterJet {
            values,
            angular_derivative,
            angular_second_derivative: &curvature.effort_second_derivative,
            parameter_derivative: &parameter.effort_derivative,
            angular_parameter_derivative: &mixed.effort_mixed_derivative,
        }
    })
}
