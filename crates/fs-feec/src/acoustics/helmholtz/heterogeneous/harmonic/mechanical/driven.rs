//! Physical observed resonances under the complete coherent acoustic/wall drive.

use super::{
    BoundaryVelocitySource, C64, HeterogeneousError, HeterogeneousHelmholtzFe, MechanicalBoundary,
    MechanicalHarmonicField, MechanicalHarmonicFrequencyResponse, MechanicalHarmonicLimits,
    ResistiveBoundary, VelocityTrace, VolumeVelocitySource, WallTractionProfile, invalid, poll,
};
use fs_math::dd::Dd;
pub use fs_phs::harmonic::resonance::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, ResonanceSurvey,
    ResonanceSurveyControls, ResonanceSurveyError,
};
use fs_phs::harmonic::resonance::{search_resonance, survey_resonances};

/// A FIXED coherent peak-phasor experiment across the requested frequency band.
/// Velocity/flow sources retain their original SI units and FE owner binding.
/// Applied traction amplitudes are Pa, positive outward like cavity pressure.
/// Values interfere in the original equations before observation; they are not
/// independent powers. Every source amplitude has zero derivative per rad/s.
#[derive(Clone, Copy, Default)]
pub struct MechanicalResonanceDrive<'s, 'm, 'f, 'c> {
    pub inward_velocity: &'s [VelocityTrace],
    pub boundary_sources: &'s [(&'s BoundaryVelocitySource<'m, 'f, 'c>, C64)],
    pub volume_sources: &'s [(&'s VolumeVelocitySource<'m, 'f, 'c>, C64)],
    /// Original profiles are projected and bound to each actual-frequency owner.
    /// A traction prepared for a different LU is never reused as a foreign owner.
    pub wall_tractions: &'s [(&'s WallTractionProfile, C64)],
}

/// One physical pressure sample in a coherent linear receiver observation.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalPressureObservation {
    pub tetrahedron: usize,
    pub barycentric: [f64; 4],
    /// Finite, signed and dimensionless; no automatic normalization occurs.
    pub weight: f64,
}

/// A frequency-independent real observation of ORIGINAL physical coordinates.
/// Negative weights retain interference and antiresonances. Wall weights follow
/// the native sorted wall_dofs order, including pressure-release-adjacent DOFs.
#[derive(Clone, Copy, Debug)]
pub enum MechanicalResonanceObservation<'a> {
    PressurePoint {
        tetrahedron: usize,
        barycentric: [f64; 4],
    },
    /// At most 128 original physical pressure samples; result units are Pa.
    CoherentPressure(&'a [MechanicalPressureObservation]),
    /// Dimensionless weights on original displacement coefficients; result m.
    WallDisplacement { wall: usize, weights: &'a [f64] },
    /// Dimensionless weights on original normal velocity coefficients; result m/s.
    WallVelocity { wall: usize, weights: &'a [f64] },
}

impl HeterogeneousHelmholtzFe<'_, '_> {
    /// Resolve one actual driven pressure or wall-motion maximum and optionally
    /// BOTH selected-amplitude half-height crossings in the supplied interval.
    ///
    /// The fixed drive includes original boundary/volume injection and applied
    /// mechanical traction as well as velocity traces. Every trial prepares the
    /// complete pressure, displacement, velocity and original Maxwell memory
    /// equations once and reuses that LU for the analytic frequency tangent.
    /// Fluid RHS differentiation includes i*b even for fixed velocity/flow;
    /// traction remains a direct momentum force, never divided by frequency.
    ///
    /// Bulk/geometry and source wall coefficients remain fixed. No wall impedance
    /// is inverted, no memory is fitted or discarded, and no damping is added.
    /// Original spatial source admission, owner identity, release-face refusal
    /// and physical field accounting remain active at every queried frequency.
    ///
    /// The interval must isolate an increasing-then-decreasing observed maximum.
    /// Missing half-height crossings, exact singularities, silent observations,
    /// antiresonances, exhausted budgets or cancellation refuse. The result is a
    /// numerical observed-amplitude peak/linewidth, not a pole-decay Q, radiated
    /// sound power, mode count or continuum completeness certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_resonance_with_coupled_mechanical_walls(
        &self,
        interval_rad_s: (f64, f64),
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        drive: MechanicalResonanceDrive<'_, '_, '_, '_>,
        observation: MechanicalResonanceObservation<'_>,
        limits: MechanicalHarmonicLimits,
        controls: ResonanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<HeterogeneousError>> {
        self.admit_mechanical_experiment(drive, observation, boundaries.len(), cancelled)
            .map_err(ResonanceSearchError::Source)?;
        search_resonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, tangent, cancelled| {
                self.mechanical_resonance_sample(
                    omega,
                    tangent,
                    resistive,
                    boundaries,
                    pressure_release_faces,
                    drive,
                    observation,
                    limits,
                    cancelled,
                )
            },
        )
    }

    /// Discover grid-resolved observed peaks with the SAME complete coherent
    /// experiment and physical observations as the local peak method above.
    /// Each actual frequency factors once and shares the LU with its tangent.
    /// The shared survey owns total evaluations/peak limits and cancellation.
    /// A grid can miss narrow or silent modes; no complete census or bandwidth
    /// is inferred. Use the local method on an isolated returned peak's interval
    /// when both half-height crossings can be supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_resonances_with_coupled_mechanical_walls(
        &self,
        frequencies_rad_s: &[f64],
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        drive: MechanicalResonanceDrive<'_, '_, '_, '_>,
        observation: MechanicalResonanceObservation<'_>,
        limits: MechanicalHarmonicLimits,
        controls: ResonanceSurveyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<HeterogeneousError>> {
        self.admit_mechanical_experiment(drive, observation, boundaries.len(), cancelled)
            .map_err(ResonanceSearchError::Source)?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, tangent, cancelled| {
                self.mechanical_resonance_sample(
                    omega,
                    tangent,
                    resistive,
                    boundaries,
                    pressure_release_faces,
                    drive,
                    observation,
                    limits,
                    cancelled,
                )
            },
        )
    }

    fn admit_mechanical_experiment(
        &self,
        drive: MechanicalResonanceDrive<'_, '_, '_, '_>,
        observation: MechanicalResonanceObservation<'_>,
        wall_count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        poll(cancelled)?;
        super::super::super::source::admit_count(drive.volume_sources.len())?;
        super::super::super::source::admit_count(drive.wall_tractions.len())?;
        super::super::super::boundary_source::admit_count(drive.boundary_sources.len())?;
        match observation {
            MechanicalResonanceObservation::PressurePoint {
                tetrahedron,
                barycentric,
            } => {
                self.validate_sample(self.ndof(), tetrahedron, barycentric)?;
            }
            MechanicalResonanceObservation::CoherentPressure(points) => {
                if points.is_empty()
                    || points.len() > 128
                    || points.iter().any(|p| !p.weight.is_finite())
                    || !points.iter().any(|p| p.weight != 0.0)
                {
                    return Err(invalid(
                        "one to 128 finite nonzero coherent pressure weights required",
                    )
                    .into());
                }
                for point in points {
                    poll(cancelled)?;
                    self.validate_sample(self.ndof(), point.tetrahedron, point.barycentric)?;
                }
            }
            MechanicalResonanceObservation::WallDisplacement { wall, weights }
            | MechanicalResonanceObservation::WallVelocity { wall, weights } => {
                if wall >= wall_count
                    || weights.is_empty()
                    || weights.len() > super::super::MAX_HARMONIC_DOFS
                    || weights.iter().any(|w| !w.is_finite())
                    || !weights.iter().any(|w| *w != 0.0)
                {
                    return Err(invalid(
                        "original wall and finite nonzero physical observation required",
                    )
                    .into());
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn mechanical_resonance_sample(
        &self,
        omega: f64,
        tangent: bool,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        drive: MechanicalResonanceDrive<'_, '_, '_, '_>,
        observation: MechanicalResonanceObservation<'_>,
        limits: MechanicalHarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(C64, C64, f64), HeterogeneousError> {
        let prepared = self.prepare_with_coupled_mechanical_walls(
            omega / core::f64::consts::TAU,
            resistive,
            boundaries,
            pressure_release_faces,
            limits,
            cancelled,
        )?;
        let mut tractions = Vec::with_capacity(drive.wall_tractions.len());
        for &(profile, amplitude) in drive.wall_tractions {
            poll(cancelled)?;
            tractions.push((
                prepared.prepare_wall_traction(profile, cancelled)?,
                amplitude,
            ));
        }
        if tangent {
            let boundary: Vec<_> = drive
                .boundary_sources
                .iter()
                .map(|&(source, value)| (source, value, C64::ZERO))
                .collect();
            let volume: Vec<_> = drive
                .volume_sources
                .iter()
                .map(|&(source, value)| (source, value, C64::ZERO))
                .collect();
            let traction: Vec<_> = tractions
                .iter()
                .map(|(source, value)| (source, *value, C64::ZERO))
                .collect();
            let result = prepared.solve_with_wall_tractions_frequency_derivative(
                drive.inward_velocity,
                None,
                &boundary,
                &volume,
                &traction,
                cancelled,
            )?;
            let (value, slope) = self.mechanical_observation(
                observation,
                &result.response,
                Some(&result),
                cancelled,
            )?;
            Ok((
                value,
                slope,
                result
                    .response
                    .acoustic
                    .backward_error
                    .max(result.derivative_backward_error),
            ))
        } else {
            let traction: Vec<_> = tractions
                .iter()
                .map(|(source, value)| (source, *value))
                .collect();
            let result = prepared.solve_with_wall_tractions(
                drive.inward_velocity,
                drive.boundary_sources,
                drive.volume_sources,
                &traction,
                cancelled,
            )?;
            let (value, _) = self.mechanical_observation(observation, &result, None, cancelled)?;
            Ok((value, C64::ZERO, result.acoustic.backward_error))
        }
    }

    fn mechanical_observation(
        &self,
        observation: MechanicalResonanceObservation<'_>,
        field: &MechanicalHarmonicField,
        tangent: Option<&MechanicalHarmonicFrequencyResponse>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(C64, C64), HeterogeneousError> {
        poll(cancelled)?;
        let mut sum = [[Dd::ZERO; 2]; 2];
        let mut pressure = |tetrahedron, barycentric, weight| -> Result<(), HeterogeneousError> {
            accumulate(
                &mut sum[0],
                weight,
                self.pressure_at(&field.acoustic.pressure_pa, tetrahedron, barycentric)?,
            );
            if let Some(tangent) = tangent {
                accumulate(
                    &mut sum[1],
                    weight,
                    self.pressure_at(
                        &tangent.pressure_derivative_pa_per_rad_s,
                        tetrahedron,
                        barycentric,
                    )?,
                );
            }
            Ok(())
        };
        match observation {
            MechanicalResonanceObservation::PressurePoint {
                tetrahedron,
                barycentric,
            } => {
                pressure(tetrahedron, barycentric, 1.0)?;
            }
            MechanicalResonanceObservation::CoherentPressure(points) => {
                for point in points {
                    poll(cancelled)?;
                    pressure(point.tetrahedron, point.barycentric, point.weight)?;
                }
            }
            MechanicalResonanceObservation::WallDisplacement { wall, weights }
            | MechanicalResonanceObservation::WallVelocity { wall, weights } => {
                let original = field
                    .walls
                    .get(wall)
                    .ok_or_else(|| invalid("original mechanical observation wall required"))?;
                let displacement = matches!(
                    observation,
                    MechanicalResonanceObservation::WallDisplacement { .. }
                );
                let values = if displacement {
                    &original.displacement_m
                } else {
                    &original.normal_velocity_m_s
                };
                if weights.len() != values.len() {
                    return Err(invalid(
                        "one observation weight per original wall coordinate required",
                    )
                    .into());
                }
                for (i, (&weight, &value)) in weights.iter().zip(values).enumerate() {
                    poll(cancelled)?;
                    accumulate(&mut sum[0], weight, value);
                    if let Some(tangent) = tangent {
                        let slope = &tangent.wall_derivatives[wall];
                        let value = if displacement {
                            slope.displacement_m_per_rad_s[i]
                        } else {
                            slope.normal_velocity_m_s_per_rad_s[i]
                        };
                        accumulate(&mut sum[1], weight, value);
                    }
                }
            }
        }
        let result = sum.map(|value| C64::new(value[0].to_f64(), value[1].to_f64()));
        if result
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(
                invalid("physical coherent resonance observation is unrepresentable").into(),
            );
        }
        poll(cancelled)?;
        Ok((result[0], result[1]))
    }
}

fn accumulate(sum: &mut [Dd; 2], weight: f64, value: C64) {
    sum[0] = sum[0] + Dd::from_f64(weight) * Dd::from_f64(value.re);
    sum[1] = sum[1] + Dd::from_f64(weight) * Dd::from_f64(value.im);
}

#[cfg(test)]
mod tests;
