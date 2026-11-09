//! Coherent optical cancellation with the same original material experiment.

use super::*;
use fs_phs::harmonic::resonance::search_antiresonance;
pub use fs_phs::harmonic::resonance::{
    AntiresonanceControls, AntiresonanceRecovery, DrivenAntiresonance,
};

impl PlanarDielectricLayers<'_> {
    /// Resolve one coherent reflected/transmitted FIELD-amplitude minimum,
    /// including an exact reflection zero without dividing by that amplitude.
    /// Every sample uses the same original materials, geometry, unit incident
    /// tangential E and declared angle/transverse-wave-number constraint as the
    /// peak API. Dispersion and all coherent finite-film interference remain.
    ///
    /// Select ReflectedElectricAmplitude, TransmittedElectricAmplitude or
    /// TransmittedMagneticAmplitude. Their raw complex fields remain analytic
    /// at a destructive zero. Power square-root observations are not admitted
    /// here because their phase-free amplitude has a cusp at zero.
    ///
    /// For reflected E, observation.norm_sq() is physical reflectance. The
    /// transmitted field choices retain their named E/H ratios; their squared
    /// amplitudes are not generally transmitted energy-flux fractions.
    ///
    /// The caller brackets an isolated decreasing-then-increasing minimum.
    /// Optional recovery_power declares a field-ratio-squared level and requires
    /// both crossings inside the interval. No phase delay or modal Q is assigned
    /// to a zero. Source-domain/underflow refusals, cancellation and the shared
    /// total evaluation budget remain active. This is a real-frequency driven
    /// cancellation, not a complex scattering-zero or free-mode computation.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_antiresonance(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
        observation: LayeredResonanceObservation,
        interval_rad_s: (f64, f64),
        controls: AntiresonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, ResonanceSearchError<ReflectionError>> {
        if !matches!(
            observation,
            LayeredResonanceObservation::ReflectedElectricAmplitude
                | LayeredResonanceObservation::TransmittedElectricAmplitude
                | LayeredResonanceObservation::TransmittedMagneticAmplitude
        ) {
            return Err(ResonanceSearchError::Source(bad(
                "coherent electric or magnetic field required for a film antiresonance",
            )));
        }
        let experiment = experiment::FrequencyExperiment::new(
            self,
            incidence_cosine,
            polarization,
            constraint,
            observation,
        )
        .map_err(ResonanceSearchError::Source)?;
        search_antiresonance(interval_rad_s, controls, cancelled, |omega, tangent, _| {
            experiment.sample(omega, tangent)
        })
    }
}
