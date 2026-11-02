//! Scalar camera radiance through native smooth, reciprocal refractive paths.
//!
//! The camera launches a unit unpolarized reciprocal beam. Each native smooth
//! boundary retains its full Fresnel Mueller response divided by its actual
//! scalar branch proposal; bulk propagation retains Beer absorption and frame
//! transport. For reciprocal paths, the reverse Mueller matrix is a transpose
//! up to polarization-convention sign matrices D. With unpolarized emission
//! and intensity-only detection, e0^T D M^T D e0 = e0^T M e0. This does NOT
//! provide a general polarized camera or an adjoint for a rough/shading normal.
//!
//! Radiance obeys L/n^2 along smooth isotropic lossless propagation. Consequently
//! the sampled reciprocal power is multiplied ONCE by (n_camera/n_endpoint)^2,
//! using the actual local phase indices. This includes continuous GRIN and
//! discrete refraction together. No bundle-area Jacobian, source inverse-square
//! factor or additional Fresnel probability belongs in this radiance estimator.
//!
//! References: Pediredla et al. (2020), "Path Tracing Estimators for Refractive Radiative Transfer",
//! Eq. 3, <https://imaging.cs.cmu.edu/rrte/assets/Renderer_AuthorCopy.pdf>;
//! Sekera (1966), polarized reciprocity, <https://doi.org/10.1364/JOSA.56.001732>.
//!
//! Finite emitters are uniform, unpolarized radiance on original opaque geometry.
//! An environment is evaluated only after proven homogeneous escape. Nonzero
//! homogeneous absorption over the remaining infinite path gives exactly zero;
//! no volume emission is assumed. Numerical distance/event limits are errors,
//! never background. Native zero power cannot distinguish exact extinction from
//! floating-point underflow and is unresolved before a positive finite source.
//! Diffuse/rough surface scattering, point-light connections, moving
//! fields, gain and nonreciprocal media remain outside this consumer.
//! With `volumes`, `participating` adds the original reciprocal molecular
//! scattering walk and retains this same endpoint radiometry.

use super::{CurvedPathResult, CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings};
use crate::charts::Ray;
use crate::dielectric::DielectricError;
use crate::lighting::EnvironmentMap;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, PathTermination};
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame, PolarizedBeam};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_exec::Cx;
use fs_geom::Vec3;

/// Physical camera and spectral pixel integration over these native paths.
pub mod raster;

/// Molecular multiple-scattering camera paths on the original curved owner.
#[cfg(feature = "volumes")]
pub mod participating;

/// Maximum knots in one immutable source spectrum, admitted before storage.
pub const MAX_RADIANCE_KNOTS: usize = 1024;

/// Authored local emitted radiance in W/(m² sr nm), with linear interpolation
/// between retained vacuum-wavelength knots. No extrapolation, RGB lift, source
/// normalization or material accuracy is inferred from these samples.
#[derive(Clone, Debug, PartialEq)]
pub struct RadianceSpectrum {
    knots: Vec<(f64, f64)>,
}

impl RadianceSpectrum {
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CameraRadianceError> {
        if !(2..=MAX_RADIANCE_KNOTS).contains(&knots.len())
            || knots.first().map(|x| x.0) != Some(LAMBDA_MIN)
            || knots.last().map(|x| x.0) != Some(LAMBDA_MAX)
            || knots
                .iter()
                .any(|&(w, r)| !w.is_finite() || !r.is_finite() || r < 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(CameraRadianceError::InvalidSpectrum);
        }
        Ok(Self { knots })
    }

    pub fn try_constant(radiance_w_m2_sr_nm: f64) -> Result<Self, CameraRadianceError> {
        Self::try_new(vec![
            (LAMBDA_MIN, radiance_w_m2_sr_nm),
            (LAMBDA_MAX, radiance_w_m2_sr_nm),
        ])
    }

    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    pub fn evaluate(&self, wavelength_nm: f64) -> Result<f64, CameraRadianceError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(CameraRadianceError::Optics(
                DielectricError::InvalidWavelength,
            ));
        }
        let upper = self.knots.partition_point(|&(w, _)| w < wavelength_nm);
        if upper == 0 || self.knots[upper].0 == wavelength_nm {
            return Ok(self.knots[upper].1);
        }
        let (lo, a) = self.knots[upper - 1];
        let (hi, b) = self.knots[upper];
        let fraction = (wavelength_nm - lo) / (hi - lo);
        finite_nonnegative((1.0 - fraction) * a + fraction * b)
    }
}

/// Per-side radiance is unchanged when both sides emit. This is not a total
/// lumen allocation: a two-sided surface radiates twice a one-sided surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmissionSides {
    /// Emit toward the original geometric winding normal.
    Front,
    Both,
}

/// Uniform unpolarized source on an actual native opaque object. Supplied
/// radiance applies in that object's surrounding local medium. It is not a
/// vacuum thermal spectrum silently promoted into a dielectric medium.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceEmitter {
    object_id: u64,
    spectrum: RadianceSpectrum,
    sides: EmissionSides,
}

impl SurfaceEmitter {
    pub const fn new(object_id: u64, spectrum: RadianceSpectrum, sides: EmissionSides) -> Self {
        Self {
            object_id,
            spectrum,
            sides,
        }
    }
    pub const fn object_id(&self) -> u64 {
        self.object_id
    }
    pub const fn spectrum(&self) -> &RadianceSpectrum {
        &self.spectrum
    }
    pub const fn sides(&self) -> EmissionSides {
        self.sides
    }
}

/// Explicit radiometric calibration of a native scene-linear RGB environment.
/// The original map's nonnegative RGB-to-spectrum reconstruction is multiplied
/// by `scale_w_m2_sr_nm`. This is a caller-authored spectrum and calibration;
/// the RGB samples do not determine a measured spectral power distribution.
/// Calibration applies in the final homogeneous surrounding medium.
#[derive(Clone, Copy, Debug)]
pub struct EnvironmentRadiance<'a> {
    map: &'a EnvironmentMap,
    scale_w_m2_sr_nm: f64,
}

impl<'a> EnvironmentRadiance<'a> {
    pub fn try_new(
        map: &'a EnvironmentMap,
        scale_w_m2_sr_nm: f64,
    ) -> Result<Self, CameraRadianceError> {
        if !scale_w_m2_sr_nm.is_finite() || scale_w_m2_sr_nm < 0.0 {
            return Err(CameraRadianceError::InvalidEnvironmentScale);
        }
        Ok(Self {
            map,
            scale_w_m2_sr_nm,
        })
    }

    pub const fn map(&self) -> &'a EnvironmentMap {
        self.map
    }

    pub const fn scale_w_m2_sr_nm(&self) -> f64 {
        self.scale_w_m2_sr_nm
    }

    fn evaluate(&self, direction: Vec3, wavelength_nm: f64) -> Result<f64, CameraRadianceError> {
        if self.scale_w_m2_sr_nm == 0.0 {
            return Ok(0.0);
        }
        let sample = self
            .map
            .evaluate(direction)
            .ok_or(CameraRadianceError::NumericalRange)?;
        finite_nonnegative(
            sample.emission.0.eval(wavelength_nm) * sample.emission.1 * self.scale_w_m2_sr_nm,
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CameraRadianceError {
    InvalidSpectrum,
    InvalidEnvironmentScale,
    InvalidEmitter {
        object_id: u64,
    },
    UnsupportedSurface {
        object_id: u64,
    },
    /// Surviving light with unresolved continuation is not a black sample.
    Unresolved(PathTermination),
    AllocationRefused,
    NumericalRange,
    Optics(DielectricError),
    Scene(CurvedSceneError),
}
impl From<CurvedSceneError> for CameraRadianceError {
    fn from(error: CurvedSceneError) -> Self {
        Self::Scene(error)
    }
}
impl From<DielectricError> for CameraRadianceError {
    fn from(error: DielectricError) -> Self {
        Self::Optics(error)
    }
}
impl From<fs_exec::Cancelled> for CameraRadianceError {
    fn from(error: fs_exec::Cancelled) -> Self {
        Self::Scene(error.into())
    }
}
impl core::fmt::Display for CameraRadianceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved camera radiance: {self:?}")
    }
}
impl core::error::Error for CameraRadianceError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadianceEndpoint {
    Surface { object_id: u64 },
    Environment,
}

/// One spectral Monte Carlo estimate, before camera exposure or spectral/pixel
/// integration. Surface branch probabilities are already included in its
/// transport weight. The retained native path carries contact resolution,
/// optical lengths, medium identity and the original termination.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CameraRadianceSample {
    radiance_w_m2_sr_nm: f64,
    transport_weight: f64,
    endpoint: RadianceEndpoint,
    path: CurvedPathResult,
}
impl CameraRadianceSample {
    pub const fn radiance_w_m2_sr_nm(&self) -> f64 {
        self.radiance_w_m2_sr_nm
    }
    pub const fn transport_weight(&self) -> f64 {
        self.transport_weight
    }
    pub const fn endpoint(&self) -> RadianceEndpoint {
        self.endpoint
    }
    pub const fn path(&self) -> &CurvedPathResult {
        &self.path
    }
}

/// Prepared emission bindings over the original curved traversal. Admission
/// rejects unsupported interactions even when one particular ray would miss
/// them. Each trace allocates no camera-owned workspace and preserves the
/// native cancellation, integration-step and geometry-work budgets.
pub struct CurvedCameraScene<'s, 'm> {
    scene: &'s CurvedPolarizedScene<'m>,
    emitters: Vec<Option<&'s SurfaceEmitter>>,
    environment: Option<EnvironmentRadiance<'s>>,
}

impl<'s, 'm> CurvedCameraScene<'s, 'm> {
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'s CurvedPolarizedScene<'m>,
        emitters: &'s [SurfaceEmitter],
        environment: Option<EnvironmentRadiance<'s>>,
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        if emitters.len() > scene.bindings.len() {
            return Err(CameraRadianceError::InvalidEmitter {
                object_id: emitters[scene.bindings.len()].object_id,
            });
        }
        let mut slots = Vec::new();
        slots
            .try_reserve_exact(scene.bindings.len())
            .map_err(|_| CameraRadianceError::AllocationRefused)?;
        for binding in &scene.bindings {
            cx.checkpoint()?;
            let supported = match binding.interaction {
                Interaction::Dielectric { surface, .. } => surface.is_delta(),
                Interaction::CoatedDielectric { .. }
                | Interaction::SmoothConductor { .. }
                | Interaction::Opaque => true,
                Interaction::RoughConductor { .. }
                | Interaction::BrushedConductor { .. }
                | Interaction::Detector { .. } => false,
            };
            if !supported {
                return Err(CameraRadianceError::UnsupportedSurface {
                    object_id: binding.object_id,
                });
            }
            slots.push(None);
        }
        for emitter in emitters {
            cx.checkpoint()?;
            let index = scene
                .bindings
                .binary_search_by_key(&emitter.object_id, |b| b.object_id)
                .map_err(|_| CameraRadianceError::InvalidEmitter {
                    object_id: emitter.object_id,
                })?;
            if slots[index].is_some() || scene.bindings[index].interaction != Interaction::Opaque {
                return Err(CameraRadianceError::InvalidEmitter {
                    object_id: emitter.object_id,
                });
            }
            slots[index] = Some(emitter);
        }
        cx.checkpoint()?;
        Ok(Self {
            scene,
            emitters: slots,
            environment,
        })
    }

    pub const fn source_scene(&self) -> &'s CurvedPolarizedScene<'m> {
        self.scene
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let path = self
            .scene
            .trace(cx, beam, initial_medium, settings, seed, path_id)?;
        self.finish(cx, path, initial_index)
    }

    /// Caller variates have the native [microfacet_u, microfacet_v, branch]
    /// convention. In this smooth-only consumer only the third selects a
    /// dielectric branch. The supplied probabilities are not resampled or
    /// conditioned on reaching an emitter.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let path = self
            .scene
            .trace_with_samples(cx, beam, initial_medium, settings, samples)?;
        self.finish(cx, path, initial_index)
    }

    fn launch(
        &self,
        ray: Ray,
        wavelength_nm: f64,
        medium: usize,
    ) -> Result<PolarizedBeam, CameraRadianceError> {
        // The launch is unpolarized, so any stable transverse basis gives the
        // same scalar response. Choosing the least-aligned Cartesian axis
        // avoids a singular reference for vertical/axial camera rays.
        let [x, y, z] = [ray.dir.x.abs(), ray.dir.y.abs(), ray.dir.z.abs()];
        let axis = if x <= y && x <= z {
            Vec3::new(1.0, 0.0, 0.0)
        } else if y <= z {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let frame = PolarizationFrame::try_new(ray.dir, axis).map_err(CurvedSceneError::from)?;
        Ok(self.scene.launch_beam(
            ray.origin,
            frame,
            wavelength_nm,
            medium,
            Stokes::unpolarized(1.0),
        )?)
    }

    fn finish(
        &self,
        cx: &Cx<'_>,
        path: CurvedPathResult,
        initial_index: f64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = path.beam();
        let wavelength = beam.wavelength_nm();
        let mut tail_transmission = 1.0;
        let (endpoint, emitted) = match path.termination() {
            PathTermination::Opaque => {
                let contact = path
                    .last_contact()
                    .ok_or(CameraRadianceError::NumericalRange)?;
                let emitter = self
                    .emitters
                    .get(contact.witness.boundary_id)
                    .ok_or(CameraRadianceError::NumericalRange)?;
                let emitted = if let Some(emitter) = emitter {
                    if emitter.object_id != contact.object_id {
                        return Err(CameraRadianceError::NumericalRange);
                    }
                    if emitter.sides == EmissionSides::Both
                        || contact.witness.normal.dot(beam.ray().dir) < 0.0
                    {
                        emitter.spectrum.evaluate(wavelength)?
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                (
                    RadianceEndpoint::Surface {
                        object_id: contact.object_id,
                    },
                    emitted,
                )
            }
            PathTermination::Escaped => {
                let BeamMedium::Homogeneous(material) = beam.medium() else {
                    return Err(CameraRadianceError::Unresolved(PathTermination::Escaped));
                };
                if material.absorption().extinction_per_m(wavelength)? > 0.0 {
                    tail_transmission = 0.0;
                }
                let emitted = if tail_transmission == 0.0 {
                    // Infinite absorbing propagation has exactly zero response,
                    // including arbitrarily bright finite calibrated sources.
                    0.0
                } else if let Some(environment) = self.environment {
                    environment.evaluate(beam.ray().dir, wavelength)?
                } else {
                    0.0
                };
                (RadianceEndpoint::Environment, emitted)
            }
            other => return Err(CameraRadianceError::Unresolved(other)),
        };
        if emitted > 0.0 && beam.power().i == 0.0 {
            return Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished,
            ));
        }
        let endpoint_index = phase_index(beam.medium(), wavelength)?;
        let ratio = initial_index / endpoint_index;
        let transport_weight = if tail_transmission == 0.0 {
            0.0
        } else {
            finite_nonnegative(beam.power().i * (ratio * ratio))?
        };
        let radiance_w_m2_sr_nm = finite_nonnegative(emitted * transport_weight)?;
        if emitted > 0.0 && radiance_w_m2_sr_nm == 0.0 {
            return Err(CameraRadianceError::NumericalRange);
        }
        cx.checkpoint()?;
        Ok(CameraRadianceSample {
            radiance_w_m2_sr_nm,
            transport_weight,
            endpoint,
            path,
        })
    }
}

fn phase_index(medium: BeamMedium, wavelength: f64) -> Result<f64, CameraRadianceError> {
    Ok(match medium {
        BeamMedium::Homogeneous(material) => material.ior().eval(wavelength)?,
        BeamMedium::GradientIndex(local) => local.sample().phase_index(),
    })
}

fn finite_nonnegative(value: f64) -> Result<f64, CameraRadianceError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(CameraRadianceError::NumericalRange)
    }
}

#[cfg(test)]
mod tests;
