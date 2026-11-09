//! Objective adjoints of isolated cycles with the original phase border intact.

use super::{AutonomousSensitivityPoint, PreparedAutonomousSensitivity, Result, direction};
use crate::sensitivity::{
    adjoint::{ForcedObjectiveGradient, PreparedForcedAdjoint},
    invalid, poll,
};

/// Total objective derivative through the periodic waveform AND solved frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AutonomousObjectiveGradient {
    /// `q_p - lambda_R^T R_p` in physical objective/parameter units. The fixed
    /// reference phase has zero direct parameter partial; the solved frequency
    /// contribution is already included by the bordered adjoint.
    pub derivative: f64,
    /// Complete balance-direction sampling change in normalized equation units
    /// per declared parameter increment, including directions orthogonal to lambda.
    pub direction_sampling_defect: f64,
    /// Contracted objective derivative's sampling change, normalized by the
    /// declared objective unit and parameter increment.
    pub gradient_sampling_defect: f64,
}

impl From<ForcedObjectiveGradient> for AutonomousObjectiveGradient {
    fn from(value: ForcedObjectiveGradient) -> Self {
        Self {
            derivative: value.derivative,
            direction_sampling_defect: value.direction_sampling_defect,
            gradient_sampling_defect: value.gradient_sampling_defect,
        }
    }
}

/// One objective's transpose solve of `[R_X R_omega; phase_X 0]`.
///
/// The complete original balance and fixed-reference phase row are retained.
/// No equation is discarded, frequency frozen, phase direction regularized, or
/// extra factorization performed. Phase-sensitive objectives refer to the SAME
/// fixed physical phase chart as the prepared cycle. Singular/nonisolated cycles
/// retain the original preparation refusal; no Floquet or harmonic-tail claim
/// follows from an accepted dual residual.
pub struct PreparedAutonomousAdjoint<'a> {
    adjoint: PreparedForcedAdjoint<'a>,
}

impl PreparedAutonomousSensitivity {
    /// Prepare the total gradient of one real objective of a complete cycle.
    ///
    /// `objective_gradient` contains EVERY ordinary real coefficient partial
    /// at fixed omega: DC states, then real and imaginary state blocks for each
    /// positive harmonic. `objective_frequency_derivative` is the independent
    /// partial q_omega at fixed coefficients, not a prescribed omega_p. Include
    /// explicit frequency and time-chart dependence here. The original bordered
    /// solve supplies omega_p; do not count that chain a second time in q_p.
    ///
    /// Positive objective_scale controls numerical acceptance only. The original
    /// sensitivity resource limits, relative residual limit and two-grid physical
    /// dual-equation limit apply, including the frequency equation and phase row.
    /// No source, state tangent, frequency-column or phase callback is repeated.
    pub fn prepare_adjoint(
        &mut self,
        objective_gradient: &[f64],
        objective_frequency_derivative: f64,
        objective_scale: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAutonomousAdjoint<'_>> {
        poll(&mut cancelled)?;
        if objective_gradient.len() != self.prepared.pack.coeff_len()
            || !objective_frequency_derivative.is_finite()
        {
            return Err(invalid(
                "autonomous adjoint needs all coefficient partials and a finite independent frequency partial",
            ));
        }
        let mut complete = Vec::with_capacity(objective_gradient.len() + 1);
        complete.extend_from_slice(objective_gradient);
        complete.push(objective_frequency_derivative);
        Ok(PreparedAutonomousAdjoint {
            adjoint: self
                .prepared
                .prepare_square_adjoint(&complete, objective_scale, cancelled)?,
        })
    }
}

impl PreparedAutonomousAdjoint<'_> {
    /// Original physical balance weights in canonical real coefficient order.
    pub fn residual_weights(&self) -> &[f64] {
        let weights = self.adjoint.residual_weights();
        &weights[..weights.len() - 1]
    }

    /// Weight of the original dimensionless fixed-reference phase row, in
    /// objective units. Its value depends on the authored phase chart.
    pub fn phase_weight(&self) -> f64 {
        let weights = self.adjoint.residual_weights();
        weights[weights.len() - 1]
    }

    /// Componentwise backward error of the complete bordered transpose solve.
    pub fn relative_residual(&self) -> f64 {
        self.adjoint.relative_residual
    }

    /// Maximum scaled physical dual-equation norm on the two original grids.
    pub fn dual_residual(&self) -> f64 {
        self.adjoint.dual_residual
    }

    /// Contract a complete physical parameter partial of the original balances
    /// at fixed coefficients AND omega. Write only the original balance rows;
    /// the unchanged reference phase has a zero direct parameter partial.
    /// `explicit_derivative` is q_p at fixed coefficients and omega, in physical
    /// objective/parameter units. Frequency-dependent ports and parameterized
    /// time/material laws belong in the complete callback as in the forward solve.
    pub fn gradient(
        &self,
        parameter_scale: f64,
        explicit_derivative: f64,
        mut residual_direction: impl FnMut(AutonomousSensitivityPoint<'_>, &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        self.adjoint
            .gradient(
                direction(parameter_scale),
                explicit_derivative,
                |point, out| {
                    residual_direction(
                        AutonomousSensitivityPoint {
                            angular_frequency: point.angular_frequency,
                            harmonics: point.harmonics,
                            aft_samples: point.aft_samples,
                            coefficients: point.coefficients,
                        },
                        out,
                    )
                },
                cancelled,
            )
            .map(Into::into)
    }

    /// Differentiate a complete time-local island with its linear port fixed.
    /// The callback supplies the direct parameter partial at unchanged state
    /// and omega, not J*x_p or R_omega*omega_p. The cached original waveforms and
    /// centered AFT owner are reused; no baseline source or per-direction solve
    /// occurs. The bordered adjoint still includes the unknown frequency response.
    pub fn gradient_island_direction(
        &mut self,
        parameter_scale: f64,
        explicit_derivative: f64,
        island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        self.adjoint
            .gradient_island_direction(
                direction(parameter_scale),
                explicit_derivative,
                island_direction,
                cancelled,
            )
            .map(Into::into)
    }
}
