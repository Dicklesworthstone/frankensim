//! Original mechanical response jets observed without another physical solve.

use super::super::{HarmonicFrequencyParameterResponse, HarmonicResponseError};
use super::{HarmonicTransferFieldParameterJet, HarmonicTransferParameterJet};
use fs_math::c64::C64;

impl HarmonicFrequencyParameterResponse {
    /// Borrow the complete displacement field jet in original coordinates.
    ///
    /// Feed this carrier to the shared complex receiver-jet projection when
    /// the physical observation itself depends on frequency or parameter.
    pub fn displacement_field_parameter_jet(&self) -> HarmonicTransferFieldParameterJet<'_> {
        let frequency = &self.frequency_response.frequency_response;
        HarmonicTransferFieldParameterJet {
            values: &frequency.response.displacement,
            angular_derivative: &frequency.displacement_derivative,
            angular_second_derivative: &self.frequency_response.displacement_second_derivative,
            parameter_derivative: &self.parameter_response.displacement_derivative,
            angular_parameter_derivative: &self.displacement_mixed_derivative,
        }
    }

    /// Borrow the complete velocity field jet, including all kinematic terms.
    ///
    /// Displacement and velocity have distinct driven extrema; all derivatives
    /// here include differentiation of v = -i*omega*q, including at DC.
    pub fn velocity_field_parameter_jet(&self) -> HarmonicTransferFieldParameterJet<'_> {
        let frequency = &self.frequency_response.frequency_response;
        HarmonicTransferFieldParameterJet {
            values: &frequency.response.velocity,
            angular_derivative: &frequency.velocity_derivative,
            angular_second_derivative: &self.frequency_response.velocity_second_derivative,
            parameter_derivative: &self.parameter_response.velocity_derivative,
            angular_parameter_derivative: &self.velocity_mixed_derivative,
        }
    }

    /// Observe the displacement parameter jet through a fixed real receiver.
    pub fn displacement_transfer_parameter(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_real(
            receiver,
            self.displacement_field_parameter_jet(),
        )
    }

    /// Observe the velocity parameter jet through a fixed real receiver.
    pub fn velocity_transfer_parameter(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_real(receiver, self.velocity_field_parameter_jet())
    }

    /// Observe displacement through fixed complex coherent transfer weights.
    pub fn displacement_transfer_parameter_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_complex(
            receiver,
            self.displacement_field_parameter_jet(),
        )
    }

    /// Observe velocity through fixed complex coherent transfer weights.
    pub fn velocity_transfer_parameter_complex(
        &self,
        receiver: &[C64],
    ) -> Result<HarmonicTransferParameterJet, HarmonicResponseError> {
        HarmonicTransferParameterJet::project_complex(receiver, self.velocity_field_parameter_jet())
    }
}
