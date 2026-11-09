//! Synchronous shaft-speed response with the original hereditary bearings.

use super::{
    DrivenResonance, MaxwellDampingBranch, ResonanceControls, ResonanceSearchError,
    RotorMaxwellHarmonicDerivative, RotorMaxwellOperatingPoint, finite, poll, strain,
};
use crate::rotordynamics::{
    RotorError,
    resonance::{RotorHarmonicDirection, RotorObservation},
    rigid::{RigidRotor, RotorUnbalance},
};
use fs_math::c64::C64;
use fs_phs::harmonic::resonance::search_resonance;

impl RotorMaxwellOperatingPoint<'_> {
    /// Response at omega = order*Omega and its TOTAL shaft-speed derivative.
    ///
    /// The finite nonzero order and shaft speed may each be signed. The source
    /// phase is unchanged: a negative order reverses the driving frequency,
    /// without an implicit conjugation. The supplied force derivative is per
    /// shaft speed, in the same original generalized force/torque coordinates.
    ///
    /// The native Maxwell tangent includes d(omega)/dOmega=order, the separate
    /// gyroscopic d(Omega)/dOmega=1, and every branch's frequency dispersion.
    /// Original stresses, stress derivatives, heat and source residuals remain
    /// available. This compares steady operating states; no angular-acceleration
    /// force or speed/temperature dependence of bearing coefficients is inferred.
    pub fn synchronous_response_with_derivative(
        &self,
        excitation_order: f64,
        force: &[C64],
        force_spin_derivative: &[C64],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellHarmonicDerivative, RotorError> {
        poll(&mut cancelled)?;
        validate_order(excitation_order)?;
        let omega = finite(excitation_order * self.point.spin_rad_s)?;
        if omega == 0.0 && self.point.spin_rad_s != 0.0 {
            return Err(RotorError::Invalid(
                "nonzero synchronous excitation frequency is unrepresentable",
            ));
        }
        self.harmonic_response_with_derivative(
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

    /// Resolve a fixed-source synchronous peak over POSITIVE SHAFT SPEEDS.
    ///
    /// Every sample rebuilds C+Omega*G and evaluates the original branch laws at
    /// omega=order*Omega. The prepared point supplies the immutable physical
    /// model and branches; its stored speed does not fix the search samples.
    ///
    /// The returned peak, brackets and half-height width are in SHAFT rad/s.
    /// For order != 1 these are not excitation angular frequencies. The
    /// observation derivative is likewise per shaft speed. The caller brackets
    /// one isolated maximum and supplies the shared total evaluation limit.
    /// This is not a Campbell crossing, free-pole Q or transient spin-up result.
    pub fn synchronous_resonance(
        &self,
        interval_spin_rad_s: (f64, f64),
        excitation_order: f64,
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        validate_order(excitation_order).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.point.model.n)
            .map_err(ResonanceSearchError::Source)?;
        let zero = [C64::ZERO; super::MAX_ROTOR_DOFS];
        search_resonance(
            interval_spin_rad_s,
            controls,
            cancelled,
            |spin, _, cancel| {
                let point = self.at_spin(spin, &mut *cancel)?;
                let result = point.synchronous_response_with_derivative(
                    excitation_order,
                    force,
                    &zero[..self.point.model.n],
                    relative_tolerance,
                    cancel,
                )?;
                observe(observation, &result)
            },
        )
    }

    // Branch admission is reusable: these original borrowed arrays are immutable.
    // Only the signed gyroscopic operating-point assembly changes between samples.
    fn at_spin(
        &self,
        spin_rad_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellOperatingPoint<'_>, RotorError> {
        Ok(RotorMaxwellOperatingPoint {
            point: self.point.model.at_speed(spin_rad_s, cancelled)?,
            branches: self.branches,
        })
    }
}

impl RigidRotor {
    /// Original 1x unbalance response with every declared Maxwell branch retained.
    ///
    /// The native rigid source supplies F=U*Omega^2 with its original station
    /// moments and phases. Its analytic derivative 2*F/Omega and both speed
    /// directions enter the same hereditary tangent solve. At zero speed the
    /// source and its derivative are zero, while a singular static model still
    /// refuses. The existing rotor stiffness is the RELAXED bearing stiffness;
    /// branch geometry and constitutive parameters are supplied by the caller.
    pub fn unbalance_maxwell_response_with_derivative(
        &self,
        spin_rad_s: f64,
        unbalances: &[RotorUnbalance],
        branches: &[MaxwellDampingBranch],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellHarmonicDerivative, RotorError> {
        let force = self.unbalance_force(spin_rad_s, unbalances, &mut cancelled)?;
        let derivative = unbalance_slope(&force, spin_rad_s)?;
        self.model()
            .at_speed_with_maxwell_memory(spin_rad_s, branches, &mut cancelled)?
            .synchronous_response_with_derivative(
                1.0,
                &force,
                &derivative,
                relative_tolerance,
                cancelled,
            )
    }

    /// Actual unbalance amplitude peak and optional half-height width in SHAFT
    /// rad/s over a caller-supplied POSITIVE speed interval.
    ///
    /// The original Omega^2 source, signed station moments, gyroscopic operator
    /// and each hereditary branch are reevaluated at every sampled operating
    /// speed. One source law and branch admission are retained across the search;
    /// no frozen-force approximation or single-frequency damping fit is used.
    /// Output Q belongs to the selected response, not a material or modal loss.
    pub fn unbalance_maxwell_resonance(
        &self,
        interval_spin_rad_s: (f64, f64),
        unbalances: &[RotorUnbalance],
        branches: &[MaxwellDampingBranch],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.model().dofs())
            .map_err(ResonanceSearchError::Source)?;
        let source = self
            .model()
            .at_speed_with_maxwell_memory(0.0, branches, &mut cancelled)
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(
            interval_spin_rad_s,
            controls,
            cancelled,
            |spin, _, cancel| {
                let force = self.unbalance_force(spin, unbalances, &mut *cancel)?;
                let derivative = unbalance_slope(&force, spin)?;
                let point = source.at_spin(spin, &mut *cancel)?;
                let result = point.synchronous_response_with_derivative(
                    1.0,
                    &force,
                    &derivative,
                    relative_tolerance,
                    cancel,
                )?;
                observe(observation, &result)
            },
        )
    }
}

fn validate_order(order: f64) -> Result<(), RotorError> {
    if !order.is_finite() || order == 0.0 {
        Err(RotorError::Invalid(
            "finite nonzero signed synchronous order required",
        ))
    } else {
        Ok(())
    }
}

fn unbalance_slope(force: &[C64; 4], spin: f64) -> Result<[C64; 4], RotorError> {
    let mut derivative = [C64::ZERO; 4];
    if spin != 0.0 {
        let slope = |value: f64| -> Result<f64, RotorError> {
            // Retain original source cancellation first; never form 2/Omega.
            let result = finite((value / spin) * 2.0)?;
            if value != 0.0 && result == 0.0 {
                return Err(RotorError::Invalid(
                    "nonzero rotor unbalance derivative is unrepresentable",
                ));
            }
            Ok(result)
        };
        for (target, value) in derivative.iter_mut().zip(force) {
            *target = C64::new(slope(value.re)?, slope(value.im)?);
        }
    }
    Ok(derivative)
}

fn observe(
    observation: RotorObservation<'_>,
    result: &RotorMaxwellHarmonicDerivative,
) -> Result<(C64, C64, f64), RotorError> {
    let (weights, values, slopes) = match observation {
        RotorObservation::Displacement(weights) => (
            weights,
            &result.response.motion.displacement,
            &result.displacement_derivative,
        ),
        RotorObservation::Velocity(weights) => (
            weights,
            &result.response.motion.velocity,
            &result.velocity_derivative,
        ),
    };
    Ok((
        strain(weights, values)?,
        strain(weights, slopes)?,
        result
            .response
            .motion
            .relative_residual
            .max(result.response.relaxation_relative_residual)
            .max(result.derivative_relative_residual)
            .max(result.relaxation_derivative_relative_residual),
    ))
}
