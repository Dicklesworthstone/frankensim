//! Actual driven rotor peaks, with complete gyroscopic and source tangents.
//!
//! For `D = K - omega² M + i omega (C + Omega G)`, differentiating the
//! ORIGINAL equilibrium gives `D q' = f' - D' q`. The primal and tangent
//! share one equilibrated LU; residuals retain the original physical rows.
//! Neither a Campbell-line crossing nor a free pole substitutes for the
//! selected displacement/velocity peak. Phasors use `exp(+i omega t)`, as
//! does the shared `fs-phs` resonance query owner.

use super::{
    RotorError, RotorHarmonicResponse, RotorOperatingPoint, action, finite, finite_complex, norm,
    poll,
};
use fs_math::{c64::C64, dd::Dd};
pub use fs_phs::harmonic::resonance::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, ResonanceSurvey,
    ResonanceSurveyControls, ResonanceSurveyError,
};
use fs_phs::harmonic::resonance::{search_resonance, survey_resonances};

mod synchronous;

/// Original physical modes, force residues and shaft-speed sensitivities.
pub mod poles;

/// Direction through forcing frequency, shaft speed and physical force phasors.
/// All derivatives are with respect to the SAME caller parameter. M/C/K/G are
/// fixed; speed-dependent bearing, prestress or material coefficients require
/// their own physical model and are not inferred from these three entries.
#[derive(Clone, Copy, Debug)]
pub struct RotorHarmonicDirection<'a> {
    /// d(omega)/d(parameter), including its sign under `exp(+i omega t)`.
    pub frequency_derivative: f64,
    /// d(Omega)/d(parameter) in the original `C + Omega G` operator.
    pub spin_derivative: f64,
    /// d(f_hat)/d(parameter) in the original generalized force/torque coordinates.
    pub force_derivative: &'a [C64],
}

/// Physical response and its analytic derivative without changing coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorHarmonicDerivative {
    /// The same response, original residual and power accounting as the native solve.
    pub response: RotorHarmonicResponse,
    /// d(q_hat)/d(parameter), with no mode normalization or modal truncation.
    pub displacement_derivative: Vec<C64>,
    /// d(i omega q_hat)/d(parameter), retaining both product-rule terms.
    pub velocity_derivative: Vec<C64>,
    /// Original differentiated equilibrium residual divided by its separate
    /// force-component norms. This is not a forward-error bound near a pole.
    pub derivative_relative_residual: f64,
}

/// Fixed real physical receiver; signed weights retain coherent interference.
/// Receiver weights and their output units belong to the caller. No norm or
/// absolute value is applied to individual coordinates before observation.
#[derive(Clone, Copy, Debug)]
pub enum RotorObservation<'a> {
    /// Weighted original displacements/angles.
    Displacement(&'a [f64]),
    /// Weighted original velocities/angular rates.
    Velocity(&'a [f64]),
}

impl RotorObservation<'_> {
    pub(super) fn validate(self, dofs: usize) -> Result<(), RotorError> {
        let (Self::Displacement(weights) | Self::Velocity(weights)) = self;
        if weights.len() != dofs || weights.iter().any(|x| !x.is_finite()) {
            return Err(RotorError::Invalid(
                "one finite rotor receiver weight per coordinate required",
            ));
        }
        Ok(())
    }

    pub(super) fn observe(
        self,
        result: &RotorHarmonicDerivative,
    ) -> Result<(C64, C64, f64), RotorError> {
        let (weights, values, derivatives) = match self {
            Self::Displacement(weights) => (
                weights,
                &result.response.displacement,
                &result.displacement_derivative,
            ),
            Self::Velocity(weights) => (
                weights,
                &result.response.velocity,
                &result.velocity_derivative,
            ),
        };
        self.validate(values.len())?;
        let mut sums = [Dd::from_f64(0.0); 4];
        for ((&weight, value), derivative) in weights.iter().zip(values).zip(derivatives) {
            for (sum, part) in
                sums.iter_mut()
                    .zip([value.re, value.im, derivative.re, derivative.im])
            {
                *sum = *sum + Dd::from_f64(weight) * Dd::from_f64(part);
            }
        }
        let values = sums.map(Dd::to_f64);
        for value in values {
            finite(value)?;
        }
        Ok((
            C64::new(values[0], values[1]),
            C64::new(values[2], values[3]),
            result
                .response
                .relative_residual
                .max(result.derivative_relative_residual),
        ))
    }
}

impl RotorOperatingPoint<'_> {
    /// Differentiate the original full equilibrium with one additional solve
    /// against its ALREADY factored operator. For a fixed-spin frequency
    /// derivative use frequency_derivative=1 and spin_derivative=0. A sweep
    /// where both change must retain both derivatives, as must a varying source.
    ///
    /// `D' = -2 omega omega' M + i omega' (C + Omega G) + i omega Omega' G`.
    /// Nonsymmetric/active bearing C/K and gyroscopic work retain their original
    /// meaning. This is a particular-solution derivative, not a stability claim.
    /// Native singularity, residual, arithmetic and cancellation failures publish
    /// neither a partial derivative nor an artificial regularized response.
    pub fn harmonic_response_with_derivative(
        &self,
        frequency_rad_s: f64,
        force: &[C64],
        direction: RotorHarmonicDirection<'_>,
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorHarmonicDerivative, RotorError> {
        poll(&mut cancelled)?;
        let m = self.model;
        if direction.force_derivative.len() != m.n
            || !finite_complex(direction.force_derivative)
            || !direction.frequency_derivative.is_finite()
            || !direction.spin_derivative.is_finite()
        {
            return Err(RotorError::Invalid(
                "finite matching rotor source and frequency/spin derivatives required",
            ));
        }
        let (response, factor, scales) = self.harmonic_response_factored(
            frequency_rad_s,
            force,
            relative_tolerance,
            &mut cancelled,
        )?;
        let w = frequency_rad_s;
        let dw = direction.frequency_derivative;
        let ds = direction.spin_derivative;
        let mut mass_slope = action(&m.mass, &response.displacement, &mut cancelled)?;
        let mut velocity_slope = action(
            &self.velocity_operator,
            &response.displacement,
            &mut cancelled,
        )?;
        let mut spin_slope = action(&m.gyro, &response.displacement, &mut cancelled)?;
        for i in 0..m.n {
            mass_slope[i] = mass_slope[i].scale(-w).scale(dw).scale(2.0);
            velocity_slope[i] = times_i(velocity_slope[i].scale(dw));
            spin_slope[i] = times_i(spin_slope[i].scale(w).scale(ds));
        }
        let mut displacement_derivative = Vec::with_capacity(m.n);
        for i in 0..m.n {
            poll(&mut cancelled)?;
            let rhs = sum([
                direction.force_derivative[i],
                -mass_slope[i],
                -velocity_slope[i],
                -spin_slope[i],
            ])?;
            displacement_derivative.push(C64::new(rhs.re / scales[i], rhs.im / scales[i]));
        }
        if !finite_complex(&displacement_derivative) {
            return Err(RotorError::Invalid(
                "scaled rotor tangent forcing is unrepresentable",
            ));
        }
        factor.solve_with_cancel(&mut displacement_derivative, &mut cancelled)?;
        let mut velocity_derivative = Vec::with_capacity(m.n);
        for (&dq, &q) in displacement_derivative.iter().zip(&response.displacement) {
            velocity_derivative.push(times_i(sum([dq.scale(w), q.scale(dw)])?));
        }
        // Audit D q' + D' q - f' using the ORIGINAL separate physical rows.
        // The velocity used in D q' is i omega q', not the total velocity jet.
        let rates: Vec<_> = displacement_derivative
            .iter()
            .map(|q| times_i(q.scale(w)))
            .collect();
        let mut inertia = action(&m.mass, &displacement_derivative, &mut cancelled)?;
        for value in &mut inertia {
            *value = value.scale(-w).scale(w);
        }
        let damping = action(&m.damping, &rates, &mut cancelled)?;
        let mut gyro = action(&m.gyro, &rates, &mut cancelled)?;
        for value in &mut gyro {
            *value = value.scale(self.spin_rad_s);
        }
        let stiffness = action(&m.stiffness, &displacement_derivative, &mut cancelled)?;
        let mut residual = Vec::with_capacity(m.n);
        for i in 0..m.n {
            residual.push(sum([
                inertia[i],
                damping[i],
                gyro[i],
                stiffness[i],
                mass_slope[i],
                velocity_slope[i],
                spin_slope[i],
                -direction.force_derivative[i],
            ])?);
        }
        let terms = [
            &inertia[..],
            &damping,
            &gyro,
            &stiffness,
            &mass_slope,
            &velocity_slope,
            &spin_slope,
            direction.force_derivative,
        ];
        let norms = terms.map(norm).into_iter().collect::<Result<Vec<_>, _>>()?;
        let scale = norms.iter().copied().fold(0.0_f64, f64::max);
        let derivative_relative_residual = if scale == 0.0 {
            norm(&residual)?
        } else {
            finite((norm(&residual)? / scale) / norms.iter().map(|x| x / scale).sum::<f64>())?
        };
        if derivative_relative_residual > relative_tolerance {
            return Err(RotorError::Residual {
                relative: derivative_relative_residual,
                limit: relative_tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(RotorHarmonicDerivative {
            response,
            displacement_derivative,
            velocity_derivative,
            derivative_relative_residual,
        })
    }

    /// Locate the selected physical peak at FIXED shaft speed and fixed complex
    /// force/torque, optionally measuring that observation's half-power width.
    /// The caller brackets one isolated maximum; this does not count free modes
    /// or equate response linewidth with decay Q. Original losses are unchanged.
    pub fn driven_resonance(
        &self,
        interval_rad_s: (f64, f64),
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.model.n)
            .map_err(ResonanceSearchError::Source)?;
        let zero = vec![C64::ZERO; self.model.n];
        search_resonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, _, cancelled| {
                let result = self.harmonic_response_with_derivative(
                    omega,
                    force,
                    RotorHarmonicDirection {
                        frequency_derivative: 1.0,
                        spin_derivative: 0.0,
                        force_derivative: &zero,
                    },
                    relative_tolerance,
                    cancelled,
                )?;
                observation.observe(&result)
            },
        )
    }

    /// Discover observable driven peaks resolved by the caller's frequency grid.
    /// Fixed spin, physical force and receiver match `driven_resonance`; a
    /// global shared query budget includes every original solve and refinement.
    /// An empty result is not proof that the rotor has no free resonances.
    pub fn survey_driven_resonances(
        &self,
        frequencies_rad_s: &[f64],
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceSurveyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.model.n)
            .map_err(ResonanceSearchError::Source)?;
        let zero = vec![C64::ZERO; self.model.n];
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, _, cancelled| {
                let result = self.harmonic_response_with_derivative(
                    omega,
                    force,
                    RotorHarmonicDirection {
                        frequency_derivative: 1.0,
                        spin_derivative: 0.0,
                        force_derivative: &zero,
                    },
                    relative_tolerance,
                    cancelled,
                )?;
                observation.observe(&result)
            },
        )
    }
}

fn times_i(z: C64) -> C64 {
    C64::new(-z.im, z.re)
}

fn sum<const N: usize>(terms: [C64; N]) -> Result<C64, RotorError> {
    let mut re = Dd::from_f64(0.0);
    let mut im = Dd::from_f64(0.0);
    for term in terms {
        finite(term.re)?;
        finite(term.im)?;
        re = re + Dd::from_f64(term.re);
        im = im + Dd::from_f64(term.im);
    }
    Ok(C64::new(finite(re.to_f64())?, finite(im.to_f64())?))
}
