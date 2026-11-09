//! Scalar camera radiance through native reciprocal refractive paths.
//!
//! The camera launches a unit unpolarized reciprocal beam. Each native
//! boundary retains its complete Mueller response divided by its actual
//! branch/directional proposal; bulk propagation retains Beer absorption and frame
//! transport. For reciprocal paths, the reverse Mueller matrix is a transpose
//! up to polarization-convention sign matrices D. With unpolarized emission
//! and intensity-only detection, e0^T D M^T D e0 = e0^T M e0.
//! [`CameraAnalyzer`] supplies an explicit passive detector row; optional
//! [`SourcePolarization`] retains the emitter's full spectral Stokes state.
//! Their original reciprocal contraction aligns the physical source frame at
//! each encounter, retaining Q/U/V and absolute throughput through this walk.
//!
//! Radiance obeys L/n^2 along smooth isotropic lossless propagation. Consequently
//! each source contribution is multiplied ONCE by (n_camera/n_source)^2,
//! using the actual local phase indices at that encounter. This includes continuous GRIN and
//! discrete refraction together. No bundle-area Jacobian, source inverse-square
//! factor or additional Fresnel probability belongs in this radiance estimator.
//! Native GGX dielectric/conductor weights already include the solid-angle
//! Jacobian, geometric cosine and masking. The default path retains their
//! single-scattering approximation, including exact-zero rejected microfacets.
//! [`SmithCameraScene`] explicitly selects complete internal conductor
//! reflections in a separate path-only estimator with its own facet budget.
//! Body-locked brushing uses the native effective radiance law, including
//! authored normal tilt: F_mueller*G2/G1(wo) under its native VNDF proposal. The
//! forward-power owner instead uses the importance adjoint with its shading-
//! normal projection correction. The camera selects its law at each facet,
//! before forming that correction; directions, PDFs and both horizons agree.
//! Decorative normals remain a shading approximation, not resolved grooves.
//!
//! References: Pediredla et al. (2020), "Path Tracing Estimators for Refractive Radiative Transfer",
//! Eq. 3, <https://imaging.cs.cmu.edu/rrte/assets/Renderer_AuthorCopy.pdf>;
//! Sekera (1966), polarized reciprocity, <https://doi.org/10.1364/JOSA.56.001732>.
//! Rough transport and shading-normal conventions: PBRT 4e, section 9.7,
//! <https://pbr-book.org/4ed/Reflection_Models/Rough_Dielectric_BSDF>, and PBRT 3e,
//! section 16.1.3, <https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation>.
//!
//! Finite emitters supply uniform radiance on original native surfaces, with
//! optional authored spectral polarization in an original body-local frame.
//! Emission is added before the same path reflects or transmits; it does not
//! replace the surface's optical interaction. Supplied radiance is independently
//! calibrated: no thermal emissivity or material/source energy coupling is inferred.
//! Explicit diffuse thermal bindings instead query the original finite-band
//! Planck/Kirchhoff owner and verify the native surface's reciprocal reflectance.
//! [`CoatingThermalEmitter`] binds the original directional s/p thermal law to
//! a matching smooth conductor face, retaining each film's own temperature.
//! An environment is evaluated only after proven homogeneous escape. Nonzero
//! homogeneous absorption over the remaining infinite path gives exactly zero;
//! explicitly authored environment Stokes fractions follow that escaped
//! direction, the original map texel/rotation, and the reciprocal camera frame.
//! optional isothermal LTE volume sources supply their own finite thermal tail.
//! Numerical distance/event limits are errors,
//! never background. Native zero power cannot distinguish exact extinction from
//! floating-point underflow and is unresolved before a positive finite source.
//! Explicit depolarizing diffuse surfaces use their native cosine proposal and
//! spectral reflectance while preserving the same curved path and medium.
//! With `photometry`, `point_lights` adds selected calibrated fixture branches
//! at native diffuse contacts. Moving fields, gain and nonreciprocal media
//! remain outside this consumer.
//! With `volumes`, `participating` adds the original reciprocal molecular
//! scattering walk and retains this same endpoint radiometry.

use super::{CurvedPathResult, CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings};
use crate::charts::Ray;
use crate::dielectric::DielectricError;
use crate::lighting::EnvironmentMap;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, PathTermination};
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame, PolarizedBeam};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum};
use fs_exec::Cx;
use fs_geom::Vec3;

mod emission;
/// Authored polarized environment observations and direct-source sampling.
pub mod environment;
mod polarized_source;
pub use polarized_source::SourcePolarization;
/// Finite-source emission histories on native group-delay paths.
pub mod transient;
use emission::SurfaceEmission;

mod thermal;
pub use thermal::{
    CoatingThermalEmitter, DiffuseThermalEmitter, IsothermalVolumeEmitter,
    SpatialVolumeEmitter, TemperatureField,
};
#[cfg(feature = "material-state")]
pub use thermal::AxialGasTemperatureField;

/// Passive polarization analysis at the original camera launch plane.
pub mod analyzer;
pub use analyzer::{CameraAnalyzer, SpectralCameraAnalyzer};

// Complete reciprocal lens state on the original surface and molecular paths.
mod sequential;

/// Complete Smith conductor walks through original reciprocal camera paths.
pub mod microsurface;
pub use microsurface::SmithCameraScene;

/// Physical camera and spectral pixel integration over these native paths.
pub mod raster;

/// Finite-area incident irradiance through original native optical images.
pub mod area_lights;

/// Molecular multiple-scattering camera paths on the original curved owner.
#[cfg(feature = "volumes")]
pub mod participating;

/// Selected calibrated point-fixture branches at native diffuse contacts.
#[cfg(feature = "photometry")]
pub mod point_lights;

/// Maximum knots in one immutable source spectrum, admitted before storage.
pub const MAX_RADIANCE_KNOTS: usize = 1024;

/// Authored local emitted radiance in W/(m² sr nm). Tabulated sources linearly
/// interpolate retained vacuum-wavelength knots; lifted sources evaluate their
/// original visible polynomial and explicit SI scale without fitting knots.
/// No extrapolation, source normalization or measured spectral accuracy is
/// inferred from either representation.
#[derive(Clone, Debug, PartialEq)]
pub struct RadianceSpectrum {
    source: RadianceSource,
}

#[derive(Clone, Debug, PartialEq)]
enum RadianceSource {
    Tabulated(Vec<(f64, f64)>),
    Lifted {
        spectrum: LiftedSpectrum,
        scale_w_m2_sr_nm: f64,
    },
}

impl RadianceSpectrum {
    /// Retain a visible source spanning exactly 380–780 nm. Nonvisible source
    /// data must use [`Self::try_new_in_band`] with its own explicit support.
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CameraRadianceError> {
        if knots.first().map(|x| x.0) != Some(LAMBDA_MIN)
            || knots.last().map(|x| x.0) != Some(LAMBDA_MAX)
        {
            return Err(CameraRadianceError::InvalidSpectrum);
        }
        Self::try_new_in_band(knots)
    }

    /// Retain original emitted radiance over its first/last positive finite
    /// vacuum wavelengths. Linear interpolation never extends that source band.
    /// This enables monochromatic radiometric paths; it does not extend the
    /// visible CIE observer or infer a spectrum from RGB data.
    pub fn try_new_in_band(knots: Vec<(f64, f64)>) -> Result<Self, CameraRadianceError> {
        if !(2..=MAX_RADIANCE_KNOTS).contains(&knots.len())
            || knots
                .iter()
                .any(|&(w, r)| !w.is_finite() || w <= 0.0 || !r.is_finite() || r < 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(CameraRadianceError::InvalidSpectrum);
        }
        Ok(Self { source: RadianceSource::Tabulated(knots) })
    }

    pub fn try_constant(radiance_w_m2_sr_nm: f64) -> Result<Self, CameraRadianceError> {
        Self::try_new(vec![
            (LAMBDA_MIN, radiance_w_m2_sr_nm),
            (LAMBDA_MAX, radiance_w_m2_sr_nm),
        ])
    }

    /// Retain an original lifted visible spectrum with explicit local radiance
    /// calibration. This constructor allocates nothing. The original polynomial
    /// is evaluated at every actual wavelength; there is no sampled surrogate.
    /// Its support is exactly the native 380–780 nm band, including for black
    /// sources. Finite coefficients are admitted here; unrepresentable evaluated
    /// radiance refuses when queried. A lifted RGB source is authored, not a
    /// recovered measured spectrum or a thermal source.
    pub fn try_lifted(
        spectrum: LiftedSpectrum,
        scale_w_m2_sr_nm: f64,
    ) -> Result<Self, CameraRadianceError> {
        if spectrum.c.iter().any(|value| !value.is_finite())
            || !scale_w_m2_sr_nm.is_finite() || scale_w_m2_sr_nm < 0.0
        {
            return Err(CameraRadianceError::InvalidSpectrum);
        }
        Ok(Self { source: RadianceSource::Lifted { spectrum, scale_w_m2_sr_nm } })
    }

    /// Original tabulated knots. An analytic lifted source has no knots and
    /// returns an empty slice; inspect [`Self::lifted_source`] for that source.
    pub fn knots(&self) -> &[(f64, f64)] {
        match &self.source {
            RadianceSource::Tabulated(knots) => knots,
            RadianceSource::Lifted { .. } => &[],
        }
    }

    /// Original lifted coefficients and SI calibration, if this is analytic.
    pub const fn lifted_source(&self) -> Option<(LiftedSpectrum, f64)> {
        match &self.source {
            RadianceSource::Lifted { spectrum, scale_w_m2_sr_nm } => Some((*spectrum, *scale_w_m2_sr_nm)),
            RadianceSource::Tabulated(_) => None,
        }
    }

    /// Inclusive source support in vacuum nanometres.
    #[must_use]
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        match &self.source {
            RadianceSource::Tabulated(knots) => (knots[0].0, knots[knots.len() - 1].0),
            RadianceSource::Lifted { .. } => (LAMBDA_MIN, LAMBDA_MAX),
        }
    }

    pub fn evaluate(&self, wavelength_nm: f64) -> Result<f64, CameraRadianceError> {
        let (minimum, maximum) = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || !(minimum..=maximum).contains(&wavelength_nm) {
            return Err(CameraRadianceError::Optics(
                DielectricError::InvalidWavelength,
            ));
        }
        match &self.source {
            RadianceSource::Tabulated(knots) => {
                let upper = knots.partition_point(|&(w, _)| w < wavelength_nm);
                if upper == 0 || knots[upper].0 == wavelength_nm {
                    return Ok(knots[upper].1);
                }
                let (lo, a) = knots[upper - 1];
                let (hi, b) = knots[upper];
                let fraction = (wavelength_nm - lo) / (hi - lo);
                finite_nonnegative((1.0 - fraction) * a + fraction * b)
            }
            RadianceSource::Lifted { spectrum, scale_w_m2_sr_nm } => {
                if *scale_w_m2_sr_nm == 0.0 { return Ok(0.0); }
                finite_nonnegative(spectrum.eval(wavelength_nm) * scale_w_m2_sr_nm)
            }
        }
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

/// Uniform source on an actual native surface, initially unpolarized. Reflection or
/// transmission continues under the original material after emission is added.
/// No thermal emissivity is inferred from the independently calibrated spectrum.
/// Supplied
/// radiance applies in that object's surrounding local medium. It is not a
/// vacuum thermal spectrum silently promoted into a dielectric medium.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceEmitter {
    object_id: u64,
    spectrum: RadianceSpectrum,
    sides: EmissionSides,
    polarization: Option<SourcePolarization>,
}

impl SurfaceEmitter {
    pub const fn new(object_id: u64, spectrum: RadianceSpectrum, sides: EmissionSides) -> Self {
        Self {
            object_id,
            spectrum,
            sides,
            polarization: None,
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

    /// Retain a physical spectral polarization in the emitting body's local
    /// coordinates. Its declared support must cover the complete radiance
    /// spectrum; neither source is extrapolated. Intensity and sidedness remain
    /// those of this emitter, with no automatic polarizer transmission loss.
    pub fn with_polarization(
        mut self,
        polarization: SourcePolarization,
    ) -> Result<Self, CameraRadianceError> {
        if let Some((low, high)) = polarization.wavelength_bounds_nm() {
            let (source_low, source_high) = self.spectrum.wavelength_bounds_nm();
            if low > source_low || high < source_high {
                return Err(CameraRadianceError::InvalidPolarization);
            }
        }
        self.polarization = Some(polarization);
        Ok(self)
    }

    pub const fn polarization(&self) -> Option<&SourcePolarization> {
        self.polarization.as_ref()
    }
}

/// Explicit radiometric calibration of a native scene-linear RGB environment.
/// The original map's nonnegative RGB-to-spectrum reconstruction is multiplied
/// by `scale_w_m2_sr_nm`. This is a caller-authored spectrum and calibration;
/// the RGB samples do not determine a measured spectral power distribution.
/// Calibration applies in the final homogeneous surrounding medium.
/// Optional explicit polarization is declared in the environment map's local
/// coordinates and follows its original rotation and piecewise-constant texels.
#[derive(Clone, Copy, Debug)]
pub struct EnvironmentRadiance<'a> {
    map: &'a EnvironmentMap,
    scale_w_m2_sr_nm: f64,
    polarization: Option<environment::EnvironmentPolarization<'a>>,
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
            polarization: None,
        })
    }

    pub const fn map(&self) -> &'a EnvironmentMap {
        self.map
    }

    pub const fn scale_w_m2_sr_nm(&self) -> f64 {
        self.scale_w_m2_sr_nm
    }

    fn evaluate(&self, direction: Vec3, wavelength_nm: f64) -> Result<f64, CameraRadianceError> {
        self.evaluate_texel(direction, wavelength_nm).map(|(radiance, _)| radiance)
    }

    fn evaluate_texel(
        &self,
        direction: Vec3,
        wavelength_nm: f64,
    ) -> Result<(f64, usize), CameraRadianceError> {
        // RGB reconstruction only declares a visible spectrum. Admitting an
        // infrared ray through native glass must not extrapolate that source.
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(CameraRadianceError::Optics(
                DielectricError::InvalidWavelength,
            ));
        }
        if self.scale_w_m2_sr_nm == 0.0 {
            return Ok((0.0, 0));
        }
        let sample = self
            .map
            .evaluate(direction)
            .ok_or(CameraRadianceError::NumericalRange)?;
        let radiance = finite_nonnegative(
            sample.emission.0.eval(wavelength_nm) * sample.emission.1 * self.scale_w_m2_sr_nm,
        )?;
        Ok((radiance, sample.texel_index))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CameraRadianceError {
    InvalidSpectrum,
    InvalidPolarization,
    InvalidEnvironmentScale,
    /// The original scene must resolve the lens's authored exterior phase
    /// index at the actual outgoing ray origin and wavelength.
    IncompatibleLensExterior,
    InvalidEmitter {
        object_id: u64,
    },
    InvalidVolumeEmitter {
        medium_index: usize,
    },
    UnsupportedSurface {
        object_id: u64,
    },
    /// Surviving light with unresolved continuation is not a black sample.
    Unresolved(PathTermination),
    AllocationRefused,
    NumericalRange,
    Optics(DielectricError),
    Thermal(crate::thermal_surface::ThermalSurfaceError),
    CoatingThermal(crate::thermal_surface::CoatingEmissionError),
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
    Surface {
        object_id: u64,
    },
    Environment,
    /// Exact zero continuation from the native single-scattering microfacet
    /// proposal. Already encountered surface emission remains in the sample.
    /// The retained path ends at this surface with its incoming beam state.
    RejectedMicrofacet {
        object_id: u64,
    },
}

/// One spectral Monte Carlo estimate, before camera exposure or spectral/pixel
/// integration. Radiance sums emission at all encountered sources. Surface
/// branch probabilities are already included in the terminal transport weight;
/// that weight does not multiply earlier source contributions. The retained native path carries contact resolution,
/// optical lengths, medium identity and the original termination. For tilted
/// brushing, the retained reciprocal beam carries the reduced-radiance
/// estimator's weight, not the public forward-light beam's importance weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CameraRadianceSample {
    radiance_w_m2_sr_nm: f64,
    launch_phase_index: f64,
    transport_weight: f64,
    endpoint: RadianceEndpoint,
    path: CurvedPathResult,
}
impl CameraRadianceSample {
    pub const fn radiance_w_m2_sr_nm(&self) -> f64 {
        self.radiance_w_m2_sr_nm
    }
    /// Actual local phase index at the camera launch, before curved transport.
    /// Optical camera consumers use this to admit their exterior medium.
    pub const fn launch_phase_index(&self) -> f64 {
        self.launch_phase_index
    }
    /// Native scalar prefix including the local phase-index conversion.
    /// An environment's source polarization is contracted separately into
    /// the reported radiance, as with source contributions along the path.
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
    // Source axes are bound once to the immutable native placement. Ordinary
    // unpolarized scenes keep this empty and allocate no additional table.
    emitter_polarization_axes: Vec<Option<Vec3>>,
    diffuse_thermal_emitters: Vec<Option<DiffuseThermalEmitter>>,
    // Independently authored positive/negative thermal faces, allocated only
    // when a directional coating source is bound to the original geometry.
    coating_thermal_emitters: Vec<[Option<CoatingThermalEmitter>; 2]>,
    environment: Option<EnvironmentRadiance<'s>>,
    volume_emitters: Vec<Option<IsothermalVolumeEmitter>>,
    spatial_volume_emitters: Vec<Option<SpatialVolumeEmitter<'s>>>,
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
                Interaction::Dielectric { .. }
                | Interaction::CoatedDielectric { .. }
                | Interaction::SmoothConductor { .. }
                | Interaction::RoughConductor { .. }
                | Interaction::BrushedConductor { .. }
                | Interaction::DepolarizingDiffuse { .. }
                | Interaction::RoughDiffuse { .. }
                | Interaction::DepolarizingLiftedDiffuse { .. }
                | Interaction::CoatedDiffuse { .. }
                | Interaction::Opaque => true,
                Interaction::Detector { .. } => false,
            };
            if !supported {
                return Err(CameraRadianceError::UnsupportedSurface {
                    object_id: binding.object_id,
                });
            }
            slots.push(None);
        }
        let mut emitter_polarization_axes = Vec::new();
        if emitters.iter().any(|emitter| emitter.polarization.is_some()) {
            emitter_polarization_axes
                .try_reserve_exact(scene.bindings.len())
                .map_err(|_| CameraRadianceError::AllocationRefused)?;
            emitter_polarization_axes.resize(scene.bindings.len(), None);
        }
        for emitter in emitters {
            cx.checkpoint()?;
            let index = scene
                .bindings
                .binary_search_by_key(&emitter.object_id, |b| b.object_id)
                .map_err(|_| CameraRadianceError::InvalidEmitter {
                    object_id: emitter.object_id,
                })?;
            if slots[index].is_some() {
                return Err(CameraRadianceError::InvalidEmitter {
                    object_id: emitter.object_id,
                });
            }
            if let Some(polarization) = &emitter.polarization {
                let instances = scene.geometry.instances();
                let instance = instances
                    .binary_search_by_key(&emitter.object_id, |instance| instance.object_id())
                    .map_err(|_| CameraRadianceError::InvalidEmitter {
                        object_id: emitter.object_id,
                    })?;
                emitter_polarization_axes[index] = Some(
                    instances[instance]
                        .transform()
                        .transform_vector(polarization.reference_axis()),
                );
            }
            slots[index] = Some(emitter);
        }
        cx.checkpoint()?;
        Ok(Self {
            scene,
            emitters: slots,
            emitter_polarization_axes,
            diffuse_thermal_emitters: Vec::new(),
            coating_thermal_emitters: Vec::new(),
            environment,
            volume_emitters: Vec::new(),
            spatial_volume_emitters: Vec::new(),
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
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    fn sample_launched_ray(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self, initial_index, beam);
        let result = self.scene.trace_with_seed_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            true,
            &self.spatial_volume_emitters,
            |event| emission.observe_event(event),
        );
        let (path, emitted) = emission.finish(result)?;
        self.finish(cx, path, initial_index, emitted)
    }

    /// Caller variates have the native [microfacet_u, microfacet_v, branch]
    /// convention. Smooth dielectrics use only the third; rough dielectrics
    /// use all three; rough/brushed conductors and diffuse surfaces use the
    /// first two. Exact-black diffuse surfaces require no draw. Proposals are
    /// not resampled or conditioned on reaching an emitter.
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
        self.sample_launched_ray_with_samples(cx, beam, initial_medium, settings, samples)
    }

    fn sample_launched_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self, initial_index, beam);
        let result = self.scene.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            true,
            &self.spatial_volume_emitters,
            |event| emission.observe_event(event),
        );
        let (path, emitted) = emission.finish(result)?;
        self.finish(cx, path, initial_index, emitted)
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
        surface_emission: f64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        if path.termination() == PathTermination::RejectedMicrofacet {
            let contact = path
                .last_contact()
                .ok_or(CameraRadianceError::NumericalRange)?;
            // A rejected native VNDF draw is a known zero of the authored
            // BSDF, not unknown continuation or underflow. Preserve the draw
            // in the Monte Carlo denominator instead of trying another facet.
            return Ok(CameraRadianceSample {
                radiance_w_m2_sr_nm: surface_emission,
                launch_phase_index: initial_index,
                transport_weight: 0.0,
                endpoint: RadianceEndpoint::RejectedMicrofacet {
                    object_id: contact.object_id,
                },
                path,
            });
        }
        let beam = path.beam();
        let wavelength = beam.wavelength_nm();
        let mut tail_transmission = 1.0;
        let (endpoint, emitted, polarization_weight) = match path.termination() {
            PathTermination::Opaque => {
                let contact = path
                    .last_contact()
                    .ok_or(CameraRadianceError::NumericalRange)?;
                // The arriving-surface callback already counted this terminal
                // emitter, including an exact-black diffuse material. Count it
                // once regardless of the interaction that stopped the walk.
                (
                    RadianceEndpoint::Surface {
                        object_id: contact.object_id,
                    },
                    0.0,
                    1.0,
                )
            }
            PathTermination::Escaped => {
                let BeamMedium::Homogeneous(material) = beam.medium() else {
                    return Err(CameraRadianceError::Unresolved(PathTermination::Escaped));
                };
                if material.absorption().extinction_per_m(wavelength)? > 0.0 {
                    // Geometry proves this ray escapes, but a point query
                    // cannot determine a varying temperature's infinite tail.
                    if self.spatial_volume_emitters.get(path.medium_index()).is_some_and(Option::is_some) {
                        return Err(CameraRadianceError::Unresolved(PathTermination::Escaped));
                    }
                    tail_transmission = 0.0;
                }
                let (emitted, polarization_weight) = if tail_transmission == 0.0 {
                    // Infinite absorbing propagation has exactly zero response,
                    // including arbitrarily bright finite calibrated sources.
                    (0.0, 1.0)
                } else if let Some(environment) = self.environment {
                    environment.evaluate_for_camera(*beam)?
                } else {
                    (0.0, 1.0)
                };
                (RadianceEndpoint::Environment, emitted, polarization_weight)
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
        let endpoint_emission = environment::endpoint_radiance(
            emitted, transport_weight, polarization_weight,
        )?;
        let radiance_w_m2_sr_nm = finite_nonnegative(surface_emission + endpoint_emission)?;
        cx.checkpoint()?;
        Ok(CameraRadianceSample {
            radiance_w_m2_sr_nm,
            launch_phase_index: initial_index,
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

#[cfg(test)]
mod rough_tests;
