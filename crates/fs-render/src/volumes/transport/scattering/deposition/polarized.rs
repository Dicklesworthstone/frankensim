//! Full-Stokes illumination on the original conservative heat/force histories.
//!
//! Extinction is isotropic in this admitted independent-particle medium. Only
//! its angular redistribution depends on polarization. Each history retains
//! that state through every collision, with a matched power-preserving sampler.
//! This does not infer dichroic absorption, aligned particles or optical spin.

use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, PolarizedBeam};
use crate::volumes::transport::scattering::polarized::{PolarizedScatteringError, RayleighGrid};

/// Original spatial/history limits plus a bounded physical angular sampler.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PolarizedDepositionConfig {
    /// Native free-flight, collision, spatial-output and reproducibility limits.
    pub transport: ScatteringDepositionConfig,
    /// Maximum scalar angular proposals per accepted polarized collision.
    /// At most two proposals are needed on average. Exhaustion refuses the
    /// whole query, never publishes a conditioned or partly deposited history.
    pub max_phase_proposals: u32,
}

impl ScatteringGrid<'_> {
    /// Spatial heating from a collimated full-Stokes source. Stokes components
    /// are already spectral-quadrature-weighted watts in the common `frame`;
    /// each wavelength keeps its own independent polarization history.
    ///
    /// Molecular, depolarized-molecular and retained spherical/Mie tables and
    /// their physical spectral mixtures are admitted before tracing. Scalar
    /// HG/isotropic/intensity-only tables cannot supply a polarized response.
    /// Pure absorption still uses the original exact Beer deposition.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_polarized_multiple_scattering_power<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        frame: PolarizationFrame,
        wavelengths_nm: &[f64; N],
        incident_stokes_w: &[Stokes; N],
        config: PolarizedDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, PolarizedScatteringError> {
        self.deposit_polarized_observed(
            cx,
            origin,
            frame,
            wavelengths_nm,
            incident_stokes_w,
            config,
            None,
        )
    }

    /// Heating and per-cell radiation force from the SAME polarized histories.
    /// Actual accepted outgoing directions determine recoil; their average
    /// cosine is never substituted. Each history closes energy and momentum.
    /// These stationary index-one forces do not include spin torque or work
    /// done on moving matter. A refusal returns neither heat nor force.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_polarized_multiple_scattering_momentum<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        frame: PolarizationFrame,
        wavelengths_nm: &[f64; N],
        incident_stokes_w: &[Stokes; N],
        config: PolarizedDepositionConfig,
    ) -> Result<ScatteringMomentumDeposition<N>, PolarizedScatteringError> {
        let mut tallies = MomentumTallies::new(config.transport.max_deposition_cells);
        let power = self.deposit_polarized_observed(
            cx,
            origin,
            frame,
            wavelengths_nm,
            incident_stokes_w,
            config,
            Some(&mut tallies),
        )?;
        let momentum = tallies.finish(&mut || checkpoint(cx), &power, frame.direction())?;
        Ok(ScatteringMomentumDeposition { power, momentum })
    }

    #[allow(clippy::too_many_arguments)]
    fn deposit_polarized_observed<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        frame: PolarizationFrame,
        wavelengths_nm: &[f64; N],
        incident_stokes_w: &[Stokes; N],
        config: PolarizedDepositionConfig,
        momentum: Option<&mut MomentumTallies<N>>,
    ) -> Result<ScatteringPowerDeposition<N>, PolarizedScatteringError> {
        checkpoint(cx)?;
        if N == 0 {
            return Err(VolumeTransportError::Sample.into());
        }
        if config.max_phase_proposals == 0 {
            return Err(VolumeTransportError::Budget.into());
        }
        validate_ray(origin, frame.direction(), f64::INFINITY)?;
        RayleighGrid::try_new(cx, self)?;
        let carrier = DielectricGlass::new(
            CauchyIor::try_constant(1.0).map_err(BeamError::from)?,
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        );
        let initial = PolarizedPhase {
            beam: PolarizedBeam::try_new(
                origin,
                frame,
                wavelengths_nm[0],
                carrier,
                Stokes::unpolarized(1.0),
            )?,
            max_proposals: config.max_phase_proposals,
        };
        let mut phases = [initial; N];
        for lane in 0..N {
            let s = incident_stokes_w[lane];
            if Stokes::try_new(s.i, s.q, s.u, s.v).is_none() {
                return Err(BeamError::InvalidStokes.into());
            }
            let normalized = if s.i == 0.0 {
                Stokes::unpolarized(1.0)
            } else {
                Stokes {
                    i: 1.0,
                    q: s.q / s.i,
                    u: s.u / s.i,
                    v: s.v / s.i,
                }
            };
            phases[lane].beam =
                PolarizedBeam::try_new(origin, frame, wavelengths_nm[lane], carrier, normalized)?;
        }
        let powers = incident_stokes_w.map(|s| s.i);
        if self.absorption.absorption_palette().is_some() || self.scattering_palette().is_some() {
            let medium = PaletteOptics::try_new(self, *wavelengths_nm, &mut || checkpoint(cx))?;
            return Ok(deposit_prepared_sampled(
                &medium,
                &mut || checkpoint(cx),
                origin,
                frame.direction(),
                &powers,
                config.transport,
                momentum,
                &phases,
                AbsorptionEstimator::Collision,
            )?);
        }
        let (absorption, scattering, extinction) = self.coefficients(wavelengths_nm)?;
        let medium = UniformOptics {
            medium: self,
            wavelengths: *wavelengths_nm,
            absorption,
            scattering,
            extinction,
        };
        Ok(deposit_prepared_sampled(
            &medium,
            &mut || checkpoint(cx),
            origin,
            frame.direction(),
            &powers,
            config.transport,
            momentum,
            &phases,
            AbsorptionEstimator::Collision,
        )?)
    }
}

// A fresh copy starts every independent history. Isotropic absorption changes
// watts in the native owner, not this unit-intensity coherency state.
#[derive(Clone, Copy)]
struct PolarizedPhase {
    beam: PolarizedBeam,
    max_proposals: u32,
}

impl DepositionPhase for PolarizedPhase {
    fn sample(
        &mut self,
        phase: &PhaseFunction,
        wavelength_nm: f64,
        direction: Vec3,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        random: &mut impl FnMut() -> f64,
    ) -> Result<Vec3, VolumeTransportError> {
        debug_assert_eq!(wavelength_nm, self.beam.wavelength_nm());
        debug_assert_eq!(direction, self.beam.frame().direction());
        let sample = phase
            .sample_polarized_power_with_checkpoint(check, self.beam, self.max_proposals, random)
            .map_err(|error| match error {
                PolarizedScatteringError::Volume(error) => error,
                PolarizedScatteringError::UnsupportedPhase => VolumeTransportError::Phase,
                _ => VolumeTransportError::Unrepresentable,
            })?;
        self.beam = sample.beam;
        Ok(self.beam.frame().direction())
    }
}

#[cfg(test)]
mod tests;
