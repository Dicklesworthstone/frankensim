//! Driven pressure peaks of the original full-space heterogeneous operator.

use super::{
    BoundaryAdmittance, HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe,
    VelocityTrace, poll,
};
use fs_phs::harmonic::resonance::{ResonanceSearchError, survey_resonances};
pub use fs_phs::harmonic::resonance::{
    ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
};

impl HeterogeneousHelmholtzFe<'_, '_> {
    /// Discover driven pressure peaks resolved by a supplied angular-frequency grid.
    ///
    /// `observation` names an ORIGINAL tetrahedron and its barycentric point.
    /// Complex peak pressure [Pa] uses `exp(+i omega t)`; all constant velocity
    /// traces interfere before observation. The tangent is per rad/s, with no
    /// pressure/velocity normalization, conjugation or finite differencing.
    /// The original hierarchical pressure basis supplies both interpolations.
    ///
    /// Material and geometry are frozen in this original operator. The uniform
    /// boundary admittances and source phasors are explicitly frequency
    /// independent across the supplied grid. A sampled dispersive wall law
    /// requires its own complete frequency derivative and is outside this API.
    /// Ideal pressure-release faces retain exact zero pressure and tangent.
    ///
    /// The existing shared search owns all frequency refinement and query/peak
    /// budgets. Each query factors the original full pressure space once, then
    /// reuses that LU for the analytic tangent with its original residual gate.
    /// No medium is resampled, no reduced model is fitted and no loss is added.
    /// `limits` applies separately to every native preparation/solve.
    ///
    /// This is a GRID-RESOLVED driven survey, not a free-mode census. Unresolved
    /// or unobservable modes may be absent; an empty result is not proof of no
    /// resonances. The shared result supplies no half-power width or decay Q.
    /// Native source failures and cancellation publish no partial survey.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_pressure_resonances(
        &self,
        frequencies_rad_s: &[f64],
        admittances: &[BoundaryAdmittance],
        pressure_release_faces: &[[u32; 3]],
        inward_velocity: &[VelocityTrace],
        observation: (usize, [f64; 4]),
        limits: HarmonicLimits,
        controls: ResonanceSurveyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<HeterogeneousError>> {
        poll(cancelled).map_err(ResonanceSearchError::Source)?;
        self.validate_sample(self.ndof(), observation.0, observation.1)
            .map_err(HeterogeneousError::from)
            .map_err(ResonanceSearchError::Source)?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, _, cancelled| {
                let prepared = self.prepare_with_pressure_release(
                    omega / core::f64::consts::TAU,
                    admittances,
                    pressure_release_faces,
                    limits,
                    cancelled,
                )?;
                let result =
                    prepared.solve_with_frequency_derivative(inward_velocity, None, cancelled)?;
                Ok((
                    self.pressure_at(&result.response.pressure_pa, observation.0, observation.1)?,
                    self.pressure_at(
                        &result.pressure_derivative_pa_per_rad_s,
                        observation.0,
                        observation.1,
                    )?,
                    result
                        .response
                        .backward_error
                        .max(result.derivative_backward_error),
                ))
            },
        )
    }
}

#[cfg(test)]
mod tests;
