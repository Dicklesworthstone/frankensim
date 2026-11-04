//! Driven peaks of original modal motion, including connections, causal material
//! memory and reciprocal air loading. The shared physical-response resolver owns
//! all bracket refinement; these adapters only supply native observations.

use super::*;
use fs_math::dd::Dd;
use fs_phs::harmonic::resonance::search_resonance;
pub use fs_phs::harmonic::resonance::{
    DrivenResonance, HalfPowerBandwidth, ResonanceControls, ResonanceFrequency,
    ResonanceSearchError,
};

/// A fixed real observation of the original mass-normalized mechanical modes.
/// Weights retain their signs and physical units without normalization. For a
/// physical point, use its original displacement/velocity participation column.
#[derive(Clone, Copy, Debug)]
pub enum ModalHarmonicObservation<'a> {
    /// Weighted displacement; weights in 1/sqrt(kg) give a result in metres.
    Displacement(&'a [f64]),
    /// Weighted velocity; the same weights give a result in metres/second.
    Velocity(&'a [f64]),
    /// Weighted physical acceleration; the same weights give metres/second².
    /// This observes `-i omega v`, retaining the inertial acceleration of the
    /// complete coupled response. Its driven peak generally differs from both
    /// displacement and velocity peaks, even for one damped mode.
    Acceleration(&'a [f64]),
}

/// Search refusals preserve the original coupled admission/solve error.
pub type ModalDrivenResonanceError = ResonanceSearchError<ModalCouplingError>;

impl ModalHarmonicObservation<'_> {
    pub(super) fn validate(self, modes: usize) -> Result<(), ModalCouplingError> {
        let (Self::Displacement(weights) | Self::Velocity(weights) | Self::Acceleration(weights)) =
            self;
        if weights.len() != modes
            || weights.iter().any(|weight| !weight.is_finite())
            || !weights.iter().any(|weight| *weight != 0.0)
        {
            return Err(invalid(
                "a finite nonzero observation must retain every original mode",
            ));
        }
        Ok(())
    }

    pub(super) fn observe(
        self,
        response: &ModalHarmonicResponse,
        tangent: Option<&ModalHarmonicFrequencyResponse>,
        gate: &CancelGate,
    ) -> Result<(C64, C64), ModalCouplingError> {
        let (weights, values, derivatives) = match self {
            Self::Displacement(weights) => (
                weights,
                &response.displacement_m_sqrt_kg,
                tangent.map(|value| &value.displacement_derivative),
            ),
            Self::Velocity(weights) | Self::Acceleration(weights) => (
                weights,
                &response.velocity_m_sqrt_kg_per_s,
                tangent.map(|value| &value.velocity_derivative),
            ),
        };
        let value = project(weights, values, gate)?;
        let derivative = derivatives
            .map(|values| project(weights, values, gate))
            .transpose()?
            .unwrap_or(C64::ZERO);
        if matches!(self, Self::Acceleration(_)) {
            // Differentiate a = -i omega v analytically. Multiplying velocity
            // once avoids forming omega², and projection precedes that scale
            // so signed modal cancellation remains in the physical observer.
            let s = C64::new(0.0, -response.omega_rad_s);
            let acceleration = s * value;
            let slope = if tangent.is_some() {
                s * derivative + C64::new(0.0, -1.0) * value
            } else {
                C64::ZERO
            };
            if !complex_finite(acceleration) || !complex_finite(slope) {
                return Err(invalid(
                    "physical acceleration observation is not representable",
                ));
            }
            Ok((acceleration, slope))
        } else {
            Ok((value, derivative))
        }
    }
}

// Signed physical projections can be much quieter than their individual modal
// contributions. Accumulate the products in the existing double-double owner;
// never replace a cancellation or a true zero by an arbitrary amplitude floor.
fn project(weights: &[f64], values: &[C64], gate: &CancelGate) -> Result<C64, ModalCouplingError> {
    if weights.len() != values.len() {
        return Err(invalid(
            "resonance observation lost an original modal coordinate",
        ));
    }
    let (mut real, mut imaginary) = (Dd::ZERO, Dd::ZERO);
    for (i, (&weight, value)) in weights.iter().zip(values).enumerate() {
        if i % 64 == 0 {
            poll(Some(gate))?;
        }
        real = real + Dd::from_f64(weight) * Dd::from_f64(value.re);
        imaginary = imaginary + Dd::from_f64(weight) * Dd::from_f64(value.im);
    }
    let result = C64::new(real.to_f64(), imaginary.to_f64());
    if !complex_finite(result) {
        return Err(invalid(
            "physical resonance observation is not representable",
        ));
    }
    Ok(result)
}

impl CoupledModalSystem {
    /// Resolve a local driven peak and optional half-power width of a FIXED
    /// physical displacement, velocity or acceleration observation under fixed force phasors.
    /// All original modes, connection directions, Maxwell arms and reciprocal
    /// air poles remain in every solve. Displacement, velocity and acceleration
    /// generally have different peak frequencies; none is replaced by a free pole.
    ///
    /// The caller must isolate one peak with increasing amplitude at the lower
    /// endpoint and decreasing amplitude at the upper endpoint. Both endpoints
    /// must also lie below its squared-amplitude half-height when width is
    /// requested. The shared resolver returns numerical brackets, not a proof of
    /// uniqueness, an interval certificate or a material decay quality factor.
    /// Silent observations, exact singularities and missing brackets refuse.
    ///
    /// Exact frequency tangents reuse each native LU. Half-height refinement
    /// only requests the primary response. `maximum_work_per_frequency` retains
    /// native complete-state admission; `controls.maximum_evaluations` bounds
    /// the number of frequency queries. Material/source validity stays with its
    /// original owner. No accepted state, history or sample clock changes.
    pub fn driven_resonance(
        &self,
        force: &[C64],
        observation: ModalHarmonicObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        maximum_work_per_frequency: usize,
        gate: &CancelGate,
    ) -> Result<DrivenResonance, ModalDrivenResonanceError> {
        poll(Some(gate)).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.mode_count())
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(
            interval_rad_s,
            controls,
            || gate.is_requested(),
            |omega, derivative, _| {
                if derivative {
                    let result = self.harmonic_response_with_frequency_derivative(
                        omega,
                        force,
                        None,
                        maximum_work_per_frequency,
                        gate,
                    )?;
                    let (value, slope) =
                        observation.observe(&result.response, Some(&result), gate)?;
                    Ok((
                        value,
                        slope,
                        result
                            .relative_derivative_residual
                            .max(result.response.solve_relative_residual),
                    ))
                } else {
                    let result =
                        self.harmonic_response(omega, force, maximum_work_per_frequency, gate)?;
                    let (value, slope) = observation.observe(&result, None, gate)?;
                    Ok((value, slope, result.solve_relative_residual))
                }
            },
        )
    }
}

#[cfg(test)]
mod tests;
