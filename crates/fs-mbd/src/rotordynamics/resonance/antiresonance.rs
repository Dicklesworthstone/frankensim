//! Actual coherent cancellation under the original instantaneous or hereditary rotor.

use super::*;
use crate::rotordynamics::{MAX_ROTOR_DOFS, relaxation::RotorMaxwellOperatingPoint};
use fs_phs::harmonic::resonance::search_antiresonance;

impl RotorOperatingPoint<'_> {
    /// Find one observed displacement/velocity minimum at fixed shaft speed
    /// and fixed original complex force/torque. Exact destructive zeros retain
    /// the complete coupled response; no coordinate is divided by its amplitude.
    ///
    /// The interval must decrease into and increase out of one isolated notch.
    /// Optional recovery_power is an explicit level of the selected observation
    /// squared, with both crossings required inside the interval. Its bandwidth
    /// is not a modal Q. A singular operator or silent experiment refuses.
    pub fn driven_antiresonance(
        &self,
        interval_rad_s: (f64, f64),
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.model.n)
            .map_err(ResonanceSearchError::Source)?;
        let zero = [C64::ZERO; MAX_ROTOR_DOFS];
        search_antiresonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, tangent, cancel| {
                if !tangent {
                    let response =
                        self.harmonic_response(omega, force, relative_tolerance, cancel)?;
                    let values = match observation {
                        RotorObservation::Displacement(_) => &response.displacement,
                        RotorObservation::Velocity(_) => &response.velocity,
                    };
                    return observation.observe_phasors(
                        values,
                        &zero[..self.model.n],
                        response.relative_residual,
                    );
                }
                let result = self.harmonic_response_with_derivative(
                    omega,
                    force,
                    RotorHarmonicDirection {
                        frequency_derivative: 1.0,
                        spin_derivative: 0.0,
                        force_derivative: &zero[..self.model.n],
                    },
                    relative_tolerance,
                    cancel,
                )?;
                observation.observe(&result)
            },
        )
    }
}

impl RotorMaxwellOperatingPoint<'_> {
    /// Find a fixed-spin coherent notch while retaining every original bearing
    /// memory, material-dispersion derivative and coupled mechanical coordinate.
    /// The source and receiver stay fixed over the band. A zero observed motion
    /// does not imply zero stress, internal motion, stored energy or branch heat.
    ///
    /// Recovery levels and bandwidth use the same original observation units as
    /// the instantaneous-bearing method. Existing source admission, full-row
    /// physical residuals, cancellation and shared query budget remain active.
    pub fn driven_antiresonance(
        &self,
        interval_rad_s: (f64, f64),
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        if force.is_empty() || force.len() > MAX_ROTOR_DOFS {
            return Err(ResonanceSearchError::Source(RotorError::Invalid(
                "bounded original rotor force coordinates required",
            )));
        }
        observation
            .validate(force.len())
            .map_err(ResonanceSearchError::Source)?;
        let zero = [C64::ZERO; MAX_ROTOR_DOFS];
        search_antiresonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, tangent, cancel| {
                if !tangent {
                    let response =
                        self.harmonic_response(omega, force, relative_tolerance, cancel)?;
                    let values = match observation {
                        RotorObservation::Displacement(_) => &response.motion.displacement,
                        RotorObservation::Velocity(_) => &response.motion.velocity,
                    };
                    return observation.observe_phasors(
                        values,
                        &zero[..force.len()],
                        response
                            .motion
                            .relative_residual
                            .max(response.relaxation_relative_residual),
                    );
                }
                let result = self.harmonic_response_with_derivative(
                    omega,
                    force,
                    RotorHarmonicDirection {
                        frequency_derivative: 1.0,
                        spin_derivative: 0.0,
                        force_derivative: &zero[..force.len()],
                    },
                    relative_tolerance,
                    cancel,
                )?;
                let (values, slopes) = match observation {
                    RotorObservation::Displacement(_) => (
                        &result.response.motion.displacement,
                        &result.displacement_derivative,
                    ),
                    RotorObservation::Velocity(_) => (
                        &result.response.motion.velocity,
                        &result.velocity_derivative,
                    ),
                };
                observation.observe_phasors(
                    values,
                    slopes,
                    result
                        .response
                        .motion
                        .relative_residual
                        .max(result.response.relaxation_relative_residual)
                        .max(result.derivative_relative_residual)
                        .max(result.relaxation_derivative_relative_residual),
                )
            },
        )
    }
}
