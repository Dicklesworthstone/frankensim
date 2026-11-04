//! Driven port-pressure peaks of the complete original DG spatial operator.
//!
//! The shared search owns frequency refinement; the duct owns every physical
//! solve, analytic tangent, material band and residual. No pole fit or reduced
//! model replaces the original gas, boundary or wall-memory equations.

use super::{C64, DgDuct, DgMedium, DuctError, PortId, invalid, poll};
use fs_phs::harmonic::resonance::{ResonanceSearchError, survey_resonances};
pub use fs_phs::harmonic::resonance::{
    ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
};

impl<M: DgMedium> DgDuct<M> {
    /// Discover and refine driven pressure peaks on the caller's frequency grid.
    ///
    /// `inflows` are fixed complex peak volume flows INTO the original ports
    /// [m3/s], with the native `exp(-i omega t)` convention. Every returned
    /// observation is pressure [Pa] at `observed_port`; its tangent is Pa per
    /// rad/s. All input phases and their interference remain in the joint solve.
    /// This is a pressure response, not an impedance normalized by one input.
    ///
    /// The grid must be finite, positive and strictly increasing. Only peaks
    /// resolved by its sampled slopes are refined; an empty result does not
    /// prove the absence of narrower or overlapping features between queries.
    /// Half-power bandwidths are absent: neighboring peaks can share a flank,
    /// and a band survey does not establish an isolated linewidth or decay Q.
    ///
    /// Each frequency tangent shares one original LU with its primary response.
    /// `maximum_work_per_response` bounds each complete native solve, while the
    /// survey controls bound their total count and the number of returned peaks.
    /// Physical material bands, original residual checks and cancellation apply
    /// to every trial. A singular lossless response or any other refusal returns
    /// no partial survey. No accepted state, time or material changes.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_port_pressure_resonances(
        &self,
        frequencies_rad_s: &[f64],
        inflows: &[(PortId, C64)],
        observed_port: PortId,
        controls: ResonanceSurveyControls,
        maximum_work_per_response: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<DuctError>> {
        poll(cancelled).map_err(ResonanceSearchError::Source)?;
        let index = self
            .continuous_operator()
            .ports()
            .iter()
            .position(|(id, _)| *id == observed_port)
            .ok_or_else(|| {
                ResonanceSearchError::Source(invalid(
                    "DG resonance observation names an undeclared port",
                ))
            })?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, derivative, cancelled| {
                if derivative {
                    let result = self.harmonic_response_with_frequency_derivative(
                        omega,
                        inflows,
                        None,
                        maximum_work_per_response,
                        relative_tolerance,
                        cancelled,
                    )?;
                    Ok((
                        result.response.ports[index].pressure_pa,
                        result.derivative.ports[index].pressure_pa,
                        result
                            .response
                            .solve_relative_residual
                            .max(result.derivative.solve_relative_residual),
                    ))
                } else {
                    let response = self.harmonic_response(
                        omega,
                        inflows,
                        maximum_work_per_response,
                        relative_tolerance,
                        cancelled,
                    )?;
                    Ok((
                        response.ports[index].pressure_pa,
                        C64::ZERO,
                        response.solve_relative_residual,
                    ))
                }
            },
        )
    }
}

#[cfg(test)]
mod tests;
