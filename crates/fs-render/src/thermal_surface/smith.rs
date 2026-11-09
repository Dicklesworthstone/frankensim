//! Thermal emission through the complete reciprocal Smith conductor walk.
//!
//! Each sampled facet emits its original positive film/substrate Kirchhoff
//! radiance into the reversed camera ray. The current Mueller prefix contracts
//! that source, then the original native reflection continues. For isothermal
//! optics the local terms telescope to `B_lambda * (1 - escaped_power)`, but
//! implementation sums positive absorption-derived sources instead of forming
//! that cancelling residual. Separate film temperatures remain separate Planck
//! sources. Geometry and every draw are shared with the native Smith sampler.
//!
//! This combines the reciprocal Smith process of Heitz et al. (2016),
//! <https://eheitzresearch.wordpress.com/240-2/>, with the local Kirchhoff source
//! already owned by [`CoatingThermalEmission`], Greffet et al. (2018),
//! <https://doi.org/10.1103/PhysRevX.8.021008>. Facets retain coherent coatings;
//! different random paths and thermal sources add incoherently. Statistical
//! height introduces no metric distance, temperature profile, delay, internal
//! heat transfer, new wavelength support or bolometric cooling law.

use super::{CoatingEmissionError, CoatingThermalEmission};
use crate::conductor::microsurface::{SmithGgx, unit};
use crate::conductor::{ConductorOptics, ConductorSurface};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::connection::{
    normalized_camera, reciprocal_contraction,
};
use crate::polarization::dielectric::transport::microsurface::{
    SmithBeamError, SmithBeamSample, SmithWalkSettings,
};
use crate::polarization::dielectric::transport::{
    BeamError, BeamMedium, PolarizationFrame, PolarizedBeam,
};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// A thermal state bound to exactly its original native opaque optical law.
/// Currently admitted native laws are bare substrates and one original
/// lossless or absorbing coating; unmatched/periodic/transmitting stacks refuse.
/// The copied optical law is small fixed state; thermal films share storage.
#[derive(Clone, Debug)]
pub struct SmithCoatingThermalEmission {
    source: CoatingThermalEmission,
    optics: ConductorOptics,
}

/// One stochastic outgoing radiance observation and the same completed
/// reciprocal reflection path. This has no outgoing marginal direction PDF.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithThermalSample {
    radiance: f64,
    reflected: SmithBeamSample,
}

impl SmithThermalSample {
    /// W/(m² sr nm), already contracted with the supplied UNIT camera Stokes
    /// importance. The original camera prefix and index conversion belong to
    /// the caller and must be applied once. No source cosine multiplies radiance.
    pub const fn radiance_w_m2_sr_nm(self) -> f64 {
        self.radiance
    }

    /// Original native reflection, with exactly the same draws, Mueller
    /// history, frames and hidden-path density as a source-free Smith walk.
    pub const fn reflected(self) -> SmithBeamSample {
        self.reflected
    }
}

/// A refused complete observation; no partial emitted sample is returned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithThermalError {
    /// Source layer and substrate laws differ from the reflected optical law.
    OpticalMismatch,
    /// The existing source emits into vacuum/exterior air only.
    VacuumRequired,
    /// The physical camera-prefix intensity must remain with its caller.
    UnitCameraRequired,
    /// Original finite-band layer absorption or Planck-source refusal.
    Emission(CoatingEmissionError),
    /// Original geometry, polarized reflection, work limit or cancellation.
    Walk(SmithBeamError),
}

impl From<CoatingEmissionError> for SmithThermalError {
    fn from(error: CoatingEmissionError) -> Self {
        Self::Emission(error)
    }
}
impl From<SmithBeamError> for SmithThermalError {
    fn from(error: SmithBeamError) -> Self {
        Self::Walk(error)
    }
}
impl From<BeamError> for SmithThermalError {
    fn from(error: BeamError) -> Self {
        Self::Walk(error.into())
    }
}
impl core::fmt::Display for SmithThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Smith conductor thermal emission: {self:?}")
    }
}
impl core::error::Error for SmithThermalError {}

impl SmithCoatingThermalEmission {
    /// Bind the original thermal source to matching native conductor optics.
    /// Matching compares complete laws and provenance, not sampled colours.
    pub fn try_new(
        source: CoatingThermalEmission,
        optics: ConductorOptics,
    ) -> Result<Self, SmithThermalError> {
        if !source.matches_smooth_conductor(&optics) {
            return Err(SmithThermalError::OpticalMismatch);
        }
        Ok(Self { source, optics })
    }

    /// Original optical films, individual temperatures and source identities.
    pub const fn source(&self) -> &CoatingThermalEmission {
        &self.source
    }

    /// The exact optical law applied at every microscopic reflection.
    pub const fn optics(&self) -> &ConductorOptics {
        &self.optics
    }

    /// Observe complete directional emission with the native isotropic Smith
    /// continuation. Smooth surfaces reduce to one exact planar observation and
    /// consume no random draws. The original macro normal is face-forwarded.
    ///
    /// `camera` travels toward the surface with intensity exactly one. The
    /// callback has the original zero-based `[flight, facet_u, facet_v]` domain.
    /// Limits/cancellation refuse the whole observation, never partial radiance.
    #[allow(clippy::too_many_arguments)]
    pub fn sample(
        &self,
        cx: &Cx<'_>,
        camera: PolarizedBeam,
        normal: Vec3,
        surface: ConductorSurface,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithThermalSample, SmithThermalError> {
        cx.checkpoint().map_err(SmithBeamError::from)?;
        self.admit(camera)?;
        let normal = incident_normal(camera, normal)?;
        if surface.is_delta() {
            let radiance = self.facet_radiance(camera, normal)?;
            let reflected = camera.sample_smith_conductor(
                cx,
                normal,
                &self.optics,
                surface,
                settings,
                samples,
            )?;
            return Ok(SmithThermalSample {
                radiance,
                reflected,
            });
        }
        let reference = if normal.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        self.sample_anisotropic(
            cx,
            camera,
            normal,
            reference,
            surface.roughness_alpha(),
            surface.roughness_alpha(),
            settings,
            samples,
        )
    }

    /// Use original GGX slope parameters in a caller-supplied physical tangent
    /// frame. `tangent` is projected into the macro plane, just as for native
    /// Stokes frames. This is geometric anisotropy, without an empirical
    /// tilted shading-normal correction or an inferred roughness measurement.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_anisotropic(
        &self,
        cx: &Cx<'_>,
        camera: PolarizedBeam,
        normal: Vec3,
        tangent: Vec3,
        alpha_x: f64,
        alpha_y: f64,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithThermalSample, SmithThermalError> {
        cx.checkpoint().map_err(SmithBeamError::from)?;
        self.admit(camera)?;
        let normal = incident_normal(camera, normal)?;
        let frame = PolarizationFrame::try_new(normal, tangent)?;
        let model = SmithGgx::try_new(alpha_x, alpha_y).map_err(SmithBeamError::from)?;
        let mut radiance = PositiveSum::default();
        let reflected = camera.sample_smith_conductor_in_frame_observed(
            cx,
            &self.optics,
            model,
            [frame.horizontal(), frame.vertical(), normal],
            settings,
            samples,
            |arriving, facet| {
                radiance.add(self.facet_radiance(arriving, facet)?)?;
                Ok::<_, SmithThermalError>(())
            },
        )?;
        Ok(SmithThermalSample {
            radiance: radiance.value()?,
            reflected,
        })
    }

    fn admit(&self, camera: PolarizedBeam) -> Result<(), SmithThermalError> {
        if camera.power().i != 1.0 {
            return Err(SmithThermalError::UnitCameraRequired);
        }
        let index = match camera.medium() {
            BeamMedium::Homogeneous(medium) => medium
                .ior()
                .eval(camera.wavelength_nm())
                .map_err(BeamError::from)?,
            BeamMedium::GradientIndex(local) => {
                if local.position() != camera.ray().origin
                    || local.wavelength_nm() != camera.wavelength_nm()
                {
                    return Err(BeamError::MediumMismatch.into());
                }
                local.sample().phase_index()
            }
        };
        if index != 1.0 {
            return Err(SmithThermalError::VacuumRequired);
        }
        Ok(())
    }

    fn facet_radiance(
        &self,
        camera: PolarizedBeam,
        normal: Vec3,
    ) -> Result<f64, SmithThermalError> {
        if camera.power().i == 0.0 {
            return Ok(0.0);
        }
        let normal = unit(normal).map_err(SmithBeamError::from)?;
        let outgoing = camera.frame().direction().scale(-1.0);
        let cosine = normal.dot(outgoing);
        let cosine = if cosine > 1.0 && cosine <= 1.0 + 16.0 * f64::EPSILON {
            1.0
        } else {
            cosine
        };
        let radiance = self
            .source
            .evaluate(camera.wavelength_nm(), cosine)?
            .polarized_radiance_w_m2_sr_nm();
        let s = Vec3::new(
            outgoing.y * normal.z - outgoing.z * normal.y,
            outgoing.z * normal.x - outgoing.x * normal.z,
            outgoing.x * normal.y - outgoing.y * normal.x,
        );
        let s = if s.x == 0.0 && s.y == 0.0 && s.z == 0.0 {
            camera.frame().horizontal()
        } else {
            s
        };
        let frame = PolarizationFrame::try_new(outgoing, s)?;
        let normalized = normalized_camera(camera, Stokes::unpolarized(1.0))?;
        let mut total = PositiveSum::default();
        for (emitted, q) in radiance.into_iter().zip([1.0, -1.0]) {
            if emitted == 0.0 {
                continue;
            }
            let overlap = reciprocal_contraction(
                normalized,
                frame,
                Stokes {
                    i: 1.0,
                    q,
                    u: 0.0,
                    v: 0.0,
                },
            )?;
            total.add(positive_product([emitted, camera.power().i, overlap])?)?;
        }
        total.value().map_err(SmithThermalError::from)
    }
}

fn incident_normal(camera: PolarizedBeam, normal: Vec3) -> Result<Vec3, SmithThermalError> {
    let mut normal = unit(normal).map_err(SmithBeamError::from)?;
    if normal.dot(camera.frame().direction()) > 0.0 {
        normal = normal.scale(-1.0);
    }
    if normal.dot(camera.frame().direction()) >= 0.0 {
        return Err(BeamError::InvalidFrame.into());
    }
    Ok(normal)
}

fn positive_product(factors: [f64; 3]) -> Result<f64, BeamError> {
    if factors
        .iter()
        .any(|factor| !factor.is_finite() || *factor < 0.0)
    {
        return Err(BeamError::NumericalRange);
    }
    if factors.contains(&0.0) {
        return Ok(0.0);
    }
    let direct = factors.into_iter().product::<f64>();
    if direct > 0.0 && direct.is_finite() {
        return Ok(direct);
    }
    let recovered = det::exp(factors.into_iter().map(det::ln).sum());
    if recovered <= 0.0 || !recovered.is_finite() {
        return Err(BeamError::NumericalRange);
    }
    Ok(recovered)
}

#[derive(Default)]
struct PositiveSum {
    sum: f64,
    correction: f64,
}
impl PositiveSum {
    fn add(&mut self, value: f64) -> Result<(), BeamError> {
        let next = self.sum + value;
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.value().map(|_| ())
    }
    fn value(&self) -> Result<f64, BeamError> {
        let value = self.sum + self.correction;
        if !value.is_finite() || value < 0.0 {
            return Err(BeamError::NumericalRange);
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
