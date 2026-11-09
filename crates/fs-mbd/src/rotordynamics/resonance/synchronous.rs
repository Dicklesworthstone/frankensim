//! Synchronous speed response with the original gyroscopic and unbalance laws.

use super::super::{
    RotorBearingModel, RotorError, finite, poll,
    rigid::{RigidRotor, RotorUnbalance},
};
use super::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, RotorHarmonicDerivative,
    RotorHarmonicDirection, RotorObservation,
};
use fs_math::c64::C64;
use fs_phs::harmonic::resonance::search_resonance;

impl RotorBearingModel {
    /// Original response along `omega = excitation_order * Omega`, and its
    /// TOTAL derivative per shaft speed [rad/s]. The order is finite and
    /// positive; signed shaft speed retains reverse rotation. The supplied
    /// force derivative is also per shaft speed in the original coordinates.
    ///
    /// Both the forcing frequency and `C + Omega G` change in this query:
    /// `D' = -2*order*omega*M + i*order*(C+Omega*G) + i*omega*G`.
    /// M/C/K/G stay fixed. Speed-dependent bearings, temperature, preload or
    /// constitutive laws require their actual coefficient derivatives and are
    /// not inferred. This compares frozen operating points; it does not
    /// integrate startup, braking or any angular-acceleration force.
    pub fn synchronous_response_with_derivative(
        &self,
        spin_rad_s: f64,
        excitation_order: f64,
        force: &[C64],
        force_spin_derivative: &[C64],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorHarmonicDerivative, RotorError> {
        poll(&mut cancelled)?;
        if !excitation_order.is_finite() || excitation_order <= 0.0 {
            return Err(RotorError::Invalid(
                "finite positive synchronous excitation order required",
            ));
        }
        let omega = finite(excitation_order * spin_rad_s)?;
        if omega == 0.0 && spin_rad_s != 0.0 {
            return Err(RotorError::Invalid(
                "nonzero synchronous excitation frequency is unrepresentable",
            ));
        }
        self.at_speed(spin_rad_s, &mut cancelled)?
            .harmonic_response_with_derivative(
                omega,
                force,
                RotorHarmonicDirection {
                    frequency_derivative: excitation_order,
                    spin_derivative: 1.0,
                    force_derivative: force_spin_derivative,
                },
                relative_tolerance,
                cancelled,
            )
    }
}

impl RigidRotor {
    /// Original 1x unbalance response and its total shaft-speed derivative.
    /// Each original station retains its phase and signed bending moment.
    /// The existing source supplies `F = U*Omega^2`; the tangent includes
    /// `F' = 2*F/Omega` as well as both frequency and gyroscopic derivatives.
    /// At zero speed F and F' are both zero, without division by zero.
    ///
    /// The nominal inertia, bearing coefficients and source geometry remain
    /// unchanged. The steady perturbation does not add eccentric mass again
    /// or stand in for transient spin acceleration. Singular original response,
    /// unrepresentable source derivatives and cancellation refuse atomically.
    pub fn unbalance_response_with_derivative(
        &self,
        spin_rad_s: f64,
        unbalances: &[RotorUnbalance],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorHarmonicDerivative, RotorError> {
        let force = self.unbalance_force(spin_rad_s, unbalances, &mut cancelled)?;
        let mut force_derivative = [C64::ZERO; 4];
        if spin_rad_s != 0.0 {
            let slope = |component: f64| -> Result<f64, RotorError> {
                // Do not form 2/Omega: that can overflow even when the physical
                // derivative is finite. Source cancellation is retained first.
                let value = finite((component / spin_rad_s) * 2.0)?;
                if component != 0.0 && value == 0.0 {
                    return Err(RotorError::Invalid(
                        "nonzero rotor unbalance derivative is unrepresentable",
                    ));
                }
                Ok(value)
            };
            for (derivative, value) in force_derivative.iter_mut().zip(force) {
                *derivative = C64::new(slope(value.re)?, slope(value.im)?);
            }
        }
        self.model().synchronous_response_with_derivative(
            spin_rad_s,
            1.0,
            &force,
            &force_derivative,
            relative_tolerance,
            cancelled,
        )
    }

    /// Resolve an observed unbalance resonance over POSITIVE shaft speeds.
    /// The source is reevaluated at every original operating point, retaining
    /// its actual Omega^2 magnitude and phases. A fixed-force frequency peak or
    /// a Campbell-line crossing generally occurs at a different speed.
    ///
    /// The shared search retains one isolated peak and optional half-height
    /// crossings of the chosen displacement/velocity observation. Its returned
    /// angular frequencies and bandwidth refer to shaft speed, which equals
    /// the 1x excitation frequency. Its Q is an observed response-width ratio,
    /// not a free-mode damping Q. The caller supplies a valid isolated interval
    /// and evaluation limits; this is not a startup trajectory or a guarantee
    /// that every resonance is found. Fixed-bearing and linearized-rotor limits
    /// are the same as `unbalance_response_with_derivative`.
    pub fn unbalance_resonance(
        &self,
        interval_rad_s: (f64, f64),
        unbalances: &[RotorUnbalance],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.model().dofs())
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(interval_rad_s, controls, cancelled, |spin, _, cancelled| {
            let result = self.unbalance_response_with_derivative(
                spin,
                unbalances,
                relative_tolerance,
                cancelled,
            )?;
            observation.observe(&result)
        })
    }
}
