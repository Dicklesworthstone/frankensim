//! Explicit analytic boundary dispersion in the original pressure equations.

use super::super::time::mechanical::MechanicalBoundary;
use super::{
    BoundaryAdmittance, C64, HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe,
    PreparedHeterogeneous, ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
    VelocityTrace, finite, frequency, invalid, poll,
};
use fs_phs::harmonic::resonance::{ResonanceSearchError, survey_resonances};

/// One uniform wall law and its analytic slope at the same declared frequency.
///
/// `Y` has units [m/(Pa s)] in `exp(+i omega t)`. Its slope is per angular
/// frequency [rad/s], even though the sample's frequency is recorded in Hz.
/// `Re(Y) >= 0` is required; either component of `dY/domega` may have either
/// sign. Supplied samples do not themselves certify causality or bandwidth.
#[derive(Clone, Debug)]
pub struct BoundaryAdmittanceFrequencySample {
    pub frequency_hz: f64,
    pub boundary: BoundaryAdmittance,
    pub derivative_m_pa_s_per_rad_s: C64,
}

impl BoundaryAdmittanceFrequencySample {
    /// Evaluate the original causal surface mass/spring/resistance and every
    /// original Maxwell arm of an existing mechanical boundary.
    ///
    /// This is exact local harmonic condensation of that declared wall law,
    /// not a shell model or complete coupled-field resonance qualification.
    /// A lossless zero of the wall impedance is refused: the reciprocal law
    /// is undefined there even when the uncondensed coupled system is regular.
    /// Material validity, spatial resolution and wall-law adequacy remain
    /// caller obligations. No fitted dispersion, artificial damping or pole
    /// clipping is introduced. Cancellation publishes no partial sample.
    pub fn mechanical(
        frequency_hz: f64,
        boundary: &MechanicalBoundary,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, HeterogeneousError> {
        let (admittance, derivative) = frequency::mechanical_admittance(
            core::f64::consts::TAU * frequency_hz,
            boundary.wall.surface_density,
            boundary.wall.stiffness_per_area,
            boundary.wall.resistance,
            boundary
                .relaxation
                .iter()
                .map(|arm| (arm.stiffness_per_area, arm.relaxation_time_s)),
            cancelled,
        )?;
        if boundary.faces.is_empty() {
            return Err(invalid("nonempty mechanical boundary faces required").into());
        }
        Ok(Self {
            frequency_hz,
            boundary: BoundaryAdmittance {
                faces: boundary.faces.clone(),
                admittance_m_pa_s: admittance,
            },
            derivative_m_pa_s_per_rad_s: derivative,
        })
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare the full pressure operator together with explicit wall slopes.
    ///
    /// Samples must match this exact frequency and occupy disjoint original
    /// faces, disjoint also from pressure-release faces. A fixed admittance
    /// can be supplied with zero slope. The same original face quadrature
    /// assembles `C_Y` and `C_Y'`; the slope is never treated as a passive heat
    /// form. Ordinary source solves retain their physical pressure and work.
    /// Reported stored energies remain the fluid's energies; condensed wall
    /// storage is not relabeled as acoustic kinetic or compression energy.
    /// Both frequency-derivative solve methods automatically include `C_Y'`
    /// with the retained LU, original residuals and coherent source slopes.
    /// No wall-law callback runs during subsequent solves.
    pub fn prepare_with_dispersive_admittances<'m>(
        &'m self,
        frequency_hz: f64,
        samples: &[BoundaryAdmittanceFrequencySample],
        pressure_release_faces: &[[u32; 3]],
        limits: HarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHeterogeneous<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let mut values = Vec::with_capacity(samples.len());
        for sample in samples {
            poll(cancelled)?;
            if sample.frequency_hz != frequency_hz || !finite(sample.derivative_m_pa_s_per_rad_s) {
                return Err(
                    invalid("matching frequency and finite boundary slope required").into(),
                );
            }
            values.push(sample.boundary.clone());
        }
        let mut prepared = self.prepare_with_pressure_release(
            frequency_hz,
            &values,
            pressure_release_faces,
            limits,
            cancelled,
        )?;
        if samples
            .iter()
            .all(|sample| sample.derivative_m_pa_s_per_rad_s == C64::ZERO)
        {
            return Ok(prepared);
        }
        let derivative = self.boundary_coefficient_matrices(
            samples.iter().map(|sample| {
                (
                    sample.boundary.faces.as_slice(),
                    sample.derivative_m_pa_s_per_rad_s,
                )
            }),
            false,
            cancelled,
        )?;
        for row in 0..self.ndof() {
            poll(cancelled)?;
            if derivative
                .0
                .row(row)
                .1
                .iter()
                .chain(derivative.1.row(row).1)
                .any(|v| !v.is_finite())
            {
                return Err(invalid("boundary frequency derivative assembly overflow").into());
            }
        }
        prepared.boundary_frequency_derivative = Some(derivative);
        Ok(prepared)
    }

    /// Survey driven pressure peaks with the original causal mechanical walls.
    ///
    /// Fixed admittances and coherent inward-velocity phasors remain constant
    /// across the grid. Every mechanical wall is resampled analytically with
    /// its complete Maxwell spectrum at each query; its derivative is consumed
    /// by the same full-pressure solve. The observation is an original
    /// tetrahedron/barycentric point in the hierarchical pressure field.
    ///
    /// The existing shared survey owns refinement and all query/peak budgets.
    /// This discovers grid-resolved driven peaks, not a free-mode census,
    /// guaranteed completeness, half-power width or decay Q. Every local wall
    /// impedance must remain finite and nonzero at every queried frequency.
    /// A local condensation failure or cancellation returns no partial survey;
    /// no pole is clipped and no artificial loss is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_pressure_resonances_with_mechanical_walls(
        &self,
        frequencies_rad_s: &[f64],
        admittances: &[BoundaryAdmittance],
        walls: &[MechanicalBoundary],
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
                let frequency_hz = omega / core::f64::consts::TAU;
                let mut samples: Vec<_> = admittances
                    .iter()
                    .map(|boundary| BoundaryAdmittanceFrequencySample {
                        frequency_hz,
                        boundary: boundary.clone(),
                        derivative_m_pa_s_per_rad_s: C64::ZERO,
                    })
                    .collect();
                for wall in walls {
                    samples.push(BoundaryAdmittanceFrequencySample::mechanical(
                        frequency_hz,
                        wall,
                        cancelled,
                    )?);
                }
                let result = self
                    .prepare_with_dispersive_admittances(
                        frequency_hz,
                        &samples,
                        pressure_release_faces,
                        limits,
                        cancelled,
                    )?
                    .solve_with_frequency_derivative(inward_velocity, None, cancelled)?;
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
