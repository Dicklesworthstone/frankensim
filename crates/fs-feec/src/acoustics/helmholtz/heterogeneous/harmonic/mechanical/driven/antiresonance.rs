//! Coherent notches in the complete original cavity/wall experiment.

use super::*;
use fs_phs::harmonic::resonance::search_antiresonance;

impl HeterogeneousHelmholtzFe<'_, '_> {
    /// Resolve one physical pressure or wall-motion minimum, including an exact
    /// destructive zero, under the same complete coherent drive as the peak API.
    /// Integrated fluid sources and wall tractions retain their original units,
    /// spatial projection and owner binding. Every wall/Maxwell state remains.
    ///
    /// The same native full-system factor supplies the response and its exact
    /// frequency tangent. No wall impedance or observed response is inverted,
    /// so a quiet receiver does not delete finite internal motion or memory.
    ///
    /// The interval must isolate a decreasing-then-increasing minimum. Optional
    /// recovery_power is an explicit level of this observation's squared
    /// amplitude (Pa², m² or m²/s²), with both crossings inside the interval.
    /// The recovered width is an observed suppression band, not pole-decay Q.
    /// Silent experiments, native singularities, missing crossings, exhausted
    /// budgets and cancellation retain explicit refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_antiresonance_with_coupled_mechanical_walls(
        &self,
        interval_rad_s: (f64, f64),
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        drive: MechanicalResonanceDrive<'_, '_, '_, '_>,
        observation: MechanicalResonanceObservation<'_>,
        limits: MechanicalHarmonicLimits,
        controls: AntiresonanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, ResonanceSearchError<HeterogeneousError>> {
        self.admit_mechanical_experiment(drive, observation, boundaries.len(), cancelled)
            .map_err(ResonanceSearchError::Source)?;
        search_antiresonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, tangent, cancel| {
                self.mechanical_resonance_sample(
                    omega,
                    tangent,
                    resistive,
                    boundaries,
                    pressure_release_faces,
                    drive,
                    observation,
                    limits,
                    cancel,
                )
            },
        )
    }
}
