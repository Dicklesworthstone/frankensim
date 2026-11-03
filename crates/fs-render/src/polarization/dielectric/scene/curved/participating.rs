//! Molecular multiple scattering along actual curved physical paths.
//!
//! Each declared medium supplies either a mass-scattering spectrum with density
//! or an explicitly bounded local spectral molecular law. Its original phase-
//! index field continues to own absorption and group delay. A scattering-only
//! flight is measured in arc length, not endpoint separation or a tangent ray.
//! Native geometry decides whether a surface interrupts that flight. At a
//! collision the original Rayleigh Mueller/PDF sampler acts on the locally
//! transported polarization frame and current physical refractive index.
//!
//! Sampling scattering separately from continuous Beer absorption is the
//! survival-biased power estimator: scattering survival is already represented
//! by flight probability, so no extra scattering extinction or albedo weight
//! multiplies the packet. Source power receives no radiance/index-squared
//! factor. A volume collision may redirect a ray toward a previously missed
//! detector; an active scattering medium cannot use a straight escape proof.
//!
//! Heterogeneous laws use majorant thinning on that same curved arc, never grid
//! chords. A local spectral law can vary both scattering rate and molecular
//! depolarization with composition and wavelength. Only molecular Rayleigh/
//! depolarized-Rayleigh laws have the requisite native Stokes response; scalar
//! particle laws refuse. All source snapshots and majorants remain immutable
//! for the complete path; local vacuum cannot prove escape from a positive
//! global majorant.

use super::{
    CurvedPathResult, CurvedPolarizedScene, CurvedSceneError, CurvedTraceEvent, CurvedTraceSettings,
    CurvedVolumeFlight, MAX_POLARIZED_PATH_EVENTS, PolarizedBeam, SceneError, SceneMedium,
};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::polarized::{
    PolarizedPhaseSample, PolarizedScatteringError,
};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use core::cell::{Cell, RefCell};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_math::det;
use fs_rand::philox::philox4x32_10;

/// Scattering survival along a solved native source-to-receiver branch.
pub mod transmittance;

/// Immutable density snapshot queried only at delta-tracking candidates.
/// Implementations must return kg/m³ for one frozen physical state.
pub trait CurvedDensityField: Sync {
    fn density_kg_m3(&self, position: Point3) -> Result<f64, VolumeTransportError>;
}

impl CurvedDensityField for VolumeGrid<'_> {
    fn density_kg_m3(&self, position: Point3) -> Result<f64, VolumeTransportError> {
        if !position.x.is_finite() || !position.y.is_finite() || !position.z.is_finite() {
            return Err(VolumeTransportError::Density { index: 0 });
        }
        let density = self.sigma_at([position.x, position.y, position.z]);
        if !density.is_finite() || density < 0.0 {
            return Err(VolumeTransportError::Density { index: 0 });
        }
        Ok(density)
    }
}

/// Local molecular scattering, with a rate per metre of physical ray arc.
/// The depolarization ratio is the unpolarized 90-degree parallel/perpendicular
/// intensity ratio used by the native Rayleigh Mueller and phase samplers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedRayleighSample {
    pub scattering_per_m: f64,
    pub depolarization_ratio: f64,
}

/// Immutable molecular optics of the SAME source-bound refractive snapshot.
///
/// The majorant must bound the local scattering rate over the entire field at
/// each admitted vacuum wavelength. A zero majorant certifies zero scattering
/// everywhere at that wavelength. Implementations must refuse unrepresentable
/// positive rates rather than rounding them to vacuum. The local phase may vary
/// with composition and wavelength; no independent density or phase inference
/// is made by this owner. Absorption remains in the original refractive field.
pub trait CurvedRayleighField: Sync {
    fn source_field_identity(&self) -> ContentHash;

    fn scattering_majorant_per_m(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError>;

    fn sample_scattering(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError>;
}

/// Borrow either the original density-scaled law or a local spectral molecular
/// field. Mixed tables retain the original arithmetic and random stream for
/// every density-scaled entry.
#[derive(Clone, Copy)]
pub enum CurvedRayleighBinding<'a> {
    Density(CurvedRayleighMedium<'a>),
    Field(&'a dyn CurvedRayleighField),
}

impl CurvedRayleighBinding<'_> {
    fn scattering_per_m(&self, wavelength_nm: f64) -> Result<f64, CurvedSceneError> {
        match self {
            Self::Density(law) => law.scattering_per_m(wavelength_nm),
            Self::Field(field) => {
                let rate = field
                    .scattering_majorant_per_m(wavelength_nm)
                    .map_err(volume_error)?;
                if !rate.is_finite() || rate < 0.0 {
                    return Err(volume_error(VolumeTransportError::ScatteringSpectrum));
                }
                Ok(rate)
            }
        }
    }

    fn scatter_candidate(
        &self,
        beam: PolarizedBeam,
        samples: [f64; 2],
        acceptance_sample: f64,
        sampled_majorant: Option<f64>,
        mut observe_accepted: impl FnMut(&PhaseFunction) -> Result<(), CurvedSceneError>,
        mut observe_departure: impl FnMut(
            &PhaseFunction,
            PolarizedPhaseSample,
        ) -> Result<(), CurvedSceneError>,
    ) -> Result<Option<PolarizedBeam>, CurvedSceneError> {
        match self {
            Self::Density(law) => {
                if !law.accepts_candidate(beam.ray().origin, acceptance_sample)? {
                    return Ok(None);
                }
                observe_accepted(law.phase)?;
                let next = law
                    .phase
                    .sample_polarized(beam, samples)
                    .map_err(scattering_error)?;
                observe_departure(law.phase, next)?;
                Ok(Some(next.beam))
            }
            Self::Field(field) => {
                if !acceptance_sample.is_finite() || !(0.0..1.0).contains(&acceptance_sample) {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let majorant = sampled_majorant.ok_or(CurvedSceneError::NumericalRange)?;
                let local = field
                    .sample_scattering(beam.ray().origin, beam.wavelength_nm())
                    .map_err(volume_error)?;
                if !local.scattering_per_m.is_finite()
                    || local.scattering_per_m < 0.0
                    || local.scattering_per_m > majorant
                {
                    return Err(volume_error(VolumeTransportError::ScatteringSpectrum));
                }
                // Validate the declared phase even at a null candidate: invalid
                // molecular optics are not an absent collision.
                let phase = PhaseFunction::rayleigh_depolarized(local.depolarization_ratio)
                    .map_err(volume_error)?;
                if local.scattering_per_m == 0.0 {
                    return Ok(None);
                }
                let probability = local.scattering_per_m / majorant;
                if !probability.is_finite() || probability <= 0.0 {
                    return Err(volume_error(VolumeTransportError::Unrepresentable));
                }
                if acceptance_sample >= probability {
                    return Ok(None);
                }
                observe_accepted(&phase)?;
                let next = phase
                    .sample_polarized(beam, samples)
                    .map_err(scattering_error)?;
                observe_departure(&phase, next)?;
                Ok(Some(next.beam))
            }
        }
    }
}

#[derive(Clone, Copy)]
enum CurvedRayleighBindings<'a, 'optics> {
    Density(&'a [Option<CurvedRayleighMedium<'optics>>]),
    Fields(&'a [Option<CurvedRayleighBinding<'optics>>]),
}

impl<'optics> CurvedRayleighBindings<'_, 'optics> {
    fn get(&self, index: usize) -> Option<CurvedRayleighBinding<'optics>> {
        match self {
            Self::Density(media) => media[index].map(CurvedRayleighBinding::Density),
            Self::Fields(media) => media[index],
        }
    }
}

#[derive(Clone, Copy)]
enum CurvedDensityProfile<'a> {
    Uniform(f64),
    Field {
        field: &'a (dyn CurvedDensityField + Sync),
        majorant_kg_m3: f64,
        identity: ContentHash,
    },
}

/// Immutable scattering law for one explicit native medium. Absorption remains
/// owned by the original refractive field. Heterogeneous density uses unbiased
/// delta tracking against an explicit global density majorant.
#[derive(Clone, Copy)]
pub struct CurvedRayleighMedium<'a> {
    scattering: &'a MassScatteringSpectrum,
    density: CurvedDensityProfile<'a>,
    phase: &'a PhaseFunction,
}

impl<'a> CurvedRayleighMedium<'a> {
    /// Admit kg/m³ and m²/kg without inferring either from refractive index.
    /// Every declared phase is checked, including a zero-density medium.
    pub fn try_new(
        scattering: &'a MassScatteringSpectrum,
        density_kg_m3: f64,
        phase: &'a PhaseFunction,
    ) -> Result<Self, CurvedSceneError> {
        if !density_kg_m3.is_finite() || density_kg_m3 < 0.0 {
            return Err(volume_error(VolumeTransportError::Density { index: 0 }));
        }
        phase.rayleigh_mueller(0.0).map_err(scattering_error)?;
        Ok(Self {
            scattering,
            density: CurvedDensityProfile::Uniform(density_kg_m3),
            phase,
        })
    }

    /// Bind the native zero-copy piecewise-constant volume grid directly.
    /// Admission validates every physical density cell, derives the exact
    /// cell majorant, and hashes the complete immutable grid snapshot. Points
    /// outside the grid are vacuum. No resampling or copied density buffer is
    /// introduced. The explicit cell budget bounds the one-time scan.
    pub fn try_new_volume_grid(
        scattering: &'a MassScatteringSpectrum,
        density: &'a VolumeGrid<'_>,
        phase: &'a PhaseFunction,
        max_grid_cells: usize,
    ) -> Result<Self, CurvedSceneError> {
        let cells = density
            .dims
            .iter()
            .try_fold(1usize, |count, dim| count.checked_mul(*dim))
            .ok_or_else(|| volume_error(VolumeTransportError::Grid))?;
        if max_grid_cells == 0 || cells != density.data.len() || cells > max_grid_cells {
            return Err(volume_error(if cells > max_grid_cells {
                VolumeTransportError::Budget
            } else {
                VolumeTransportError::Grid
            }));
        }
        if density.dims.contains(&0)
            || density
                .origin
                .iter()
                .chain(density.cell.iter())
                .any(|value| !value.is_finite())
            || density.cell.iter().any(|value| *value <= 0.0)
        {
            return Err(volume_error(VolumeTransportError::Grid));
        }
        let mut majorant = 0.0_f64;
        let mut hash = DomainHasher::new("org.frankensim.render.curved-density-grid.v1");
        for axis in 0..3 {
            hash.update(&(density.dims[axis] as u64).to_le_bytes());
            hash.update(&density.origin[axis].to_bits().to_le_bytes());
            hash.update(&density.cell[axis].to_bits().to_le_bytes());
        }
        for (index, value) in density.data.iter().copied().enumerate() {
            if !value.is_finite() || value < 0.0 {
                return Err(volume_error(VolumeTransportError::Density { index }));
            }
            majorant = majorant.max(value);
            hash.update(&value.to_bits().to_le_bytes());
        }
        if majorant == 0.0 {
            return Self::try_new(scattering, 0.0, phase);
        }
        Self::try_new_heterogeneous(scattering, density, majorant, hash.finalize(), phase)
    }

    /// Bind a frozen heterogeneous density field. The positive finite majorant
    /// must bound the field everywhere this medium can be reached.
    pub fn try_new_heterogeneous(
        scattering: &'a MassScatteringSpectrum,
        density: &'a (dyn CurvedDensityField + Sync),
        density_majorant_kg_m3: f64,
        density_identity: ContentHash,
        phase: &'a PhaseFunction,
    ) -> Result<Self, CurvedSceneError> {
        if !density_majorant_kg_m3.is_finite()
            || density_majorant_kg_m3 <= 0.0
            || density_identity == ContentHash([0; 32])
        {
            return Err(volume_error(VolumeTransportError::Density { index: 0 }));
        }
        phase.rayleigh_mueller(0.0).map_err(scattering_error)?;
        Ok(Self {
            scattering,
            density: CurvedDensityProfile::Field {
                field: density,
                majorant_kg_m3: density_majorant_kg_m3,
                identity: density_identity,
            },
            phase,
        })
    }

    pub const fn scattering_spectrum(&self) -> &'a MassScatteringSpectrum {
        self.scattering
    }
    /// Uniform density, or the global majorant for a heterogeneous field.
    pub const fn density_kg_m3(&self) -> f64 {
        match self.density {
            CurvedDensityProfile::Uniform(value) => value,
            CurvedDensityProfile::Field { majorant_kg_m3, .. } => majorant_kg_m3,
        }
    }
    pub const fn is_heterogeneous(&self) -> bool {
        matches!(self.density, CurvedDensityProfile::Field { .. })
    }
    pub const fn density_identity(&self) -> Option<ContentHash> {
        match self.density {
            CurvedDensityProfile::Uniform(_) => None,
            CurvedDensityProfile::Field { identity, .. } => Some(identity),
        }
    }
    pub const fn phase(&self) -> &'a PhaseFunction {
        self.phase
    }

    /// Majorant scattering rate per metre of ray arc length.
    pub fn scattering_per_m(&self, wavelength_nm: f64) -> Result<f64, CurvedSceneError> {
        let mass = self.scattering.at(wavelength_nm).map_err(volume_error)?;
        let density = self.density_kg_m3();
        let rate = mass * density;
        if !rate.is_finite() || (mass > 0.0 && density > 0.0 && rate == 0.0) {
            return Err(volume_error(VolumeTransportError::Unrepresentable));
        }
        Ok(rate)
    }

    fn accepts_candidate(
        &self,
        position: Point3,
        acceptance_sample: f64,
    ) -> Result<bool, CurvedSceneError> {
        match self.density {
            CurvedDensityProfile::Uniform(_) => Ok(true),
            CurvedDensityProfile::Field {
                field,
                majorant_kg_m3,
                ..
            } => {
                if !acceptance_sample.is_finite() || !(0.0..1.0).contains(&acceptance_sample) {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let density = field.density_kg_m3(position).map_err(volume_error)?;
                if !density.is_finite() || density < 0.0 || density > majorant_kg_m3 {
                    return Err(volume_error(VolumeTransportError::Density { index: 0 }));
                }
                if density == 0.0 {
                    return Ok(false);
                }
                let probability = density / majorant_kg_m3;
                if !probability.is_finite() || probability <= 0.0 {
                    return Err(volume_error(VolumeTransportError::Unrepresentable));
                }
                Ok(acceptance_sample < probability)
            }
        }
    }
}

/// Independent physical-distance, surface-event and volume-collision limits.
#[derive(Clone, Copy, Debug)]
pub struct CurvedRayleighTraceSettings {
    pub surfaces: CurvedTraceSettings,
    pub max_volume_events: u32,
}

/// The original detector-consumable path plus actual volume-work counts.
/// Surviving power at a volume cap has `PathTermination::VolumeLimit`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedRayleighPath {
    path: CurvedPathResult,
    volume_events: u32,
    free_flight_draws: u32,
}

/// Original camera arrival, before a surface response or an accepted molecular
/// phase sample. Null candidates never produce a volume vertex. The local
/// phase is the exact admitted law used for the unchanged continuation.
pub(super) enum CurvedRayleighVertex<'a> {
    Surface(&'a CurvedPathResult),
    Volume {
        path: &'a CurvedPathResult,
        phase: &'a PhaseFunction,
    },
}

/// Original surface and molecular proposals in physical path order. A real
/// collision has one arrival followed by its successful sampled departure;
/// rejected tracking candidates produce neither event. The phase PDF is the
/// native solid-angle proposal, separate from the Mueller-weighted beam power.
/// No inferred no-scatter density or extra survival weight is introduced.
pub(super) enum CurvedRayleighEvent<'a> {
    Surface(CurvedTraceEvent<'a>),
    VolumeArrival {
        path: &'a CurvedPathResult,
        phase: &'a PhaseFunction,
    },
    VolumeDeparture {
        /// Original path counters with the native outgoing beam. Its medium
        /// and position are unchanged by this local molecular response.
        path: &'a CurvedPathResult,
        incoming_beam: PolarizedBeam,
        phase: &'a PhaseFunction,
        pdf: f64,
    },
}
impl CurvedRayleighPath {
    pub const fn path(&self) -> &CurvedPathResult {
        &self.path
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
}

/// Borrowed original geometry and one optional scattering law per medium.
/// Binding and tracing allocate no heap workspace or copied optical tables.
pub struct CurvedRayleighScene<'a, 'scene, 'optics> {
    scene: &'a CurvedPolarizedScene<'scene>,
    media: CurvedRayleighBindings<'a, 'optics>,
    radiance_brushing: bool,
}
impl<'a, 'scene, 'optics> CurvedRayleighScene<'a, 'scene, 'optics> {
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'a CurvedPolarizedScene<'scene>,
        media: &'a [Option<CurvedRayleighMedium<'optics>>],
    ) -> Result<Self, CurvedSceneError> {
        cx.checkpoint()?;
        if media.len() != scene.media.len() {
            return Err(SceneError::InvalidBindings.into());
        }
        Ok(Self {
            scene,
            media: CurvedRayleighBindings::Density(media),
            radiance_brushing: false,
        })
    }

    /// Bind local spectral molecular fields to their exact original GRIN
    /// source/snapshot identities. An equal index at one point cannot substitute
    /// for that identity; homogeneous media cannot stand in for source fields.
    /// Spectral rate and phase admission occurs at the queried wavelength.
    pub fn try_new_with_fields(
        cx: &Cx<'_>,
        scene: &'a CurvedPolarizedScene<'scene>,
        media: &'a [Option<CurvedRayleighBinding<'optics>>],
    ) -> Result<Self, CurvedSceneError> {
        cx.checkpoint()?;
        if media.len() != scene.media.len() {
            return Err(SceneError::InvalidBindings.into());
        }
        for (index, (law, medium)) in media.iter().zip(scene.media).enumerate() {
            if index % 1024 == 0 {
                cx.checkpoint()?;
            }
            if let Some(CurvedRayleighBinding::Field(field)) = law {
                let source_identity = field.source_field_identity();
                if source_identity == ContentHash([0; 32]) {
                    return Err(CurvedSceneError::MissingFieldIdentity);
                }
                match medium {
                    SceneMedium::GradientIndex { identity, .. }
                    | SceneMedium::Differentiable { identity, .. }
                        if *identity == source_identity => {}
                    _ => return Err(SceneError::InvalidBindings.into()),
                }
            }
        }
        Ok(Self {
            scene,
            media: CurvedRayleighBindings::Fields(media),
            radiance_brushing: false,
        })
    }

    pub const fn surface_scene(&self) -> &'a CurvedPolarizedScene<'scene> {
        self.scene
    }

    // A stack-only borrowed view; original power paths and material/field
    // identities remain unchanged. All three sampling APIs share this mode.
    pub(super) fn for_camera(&self) -> Self {
        Self {
            scene: self.scene,
            media: self.media,
            radiance_brushing: true,
        }
    }

    /// Surface draws keep the native domain and surface ordinal. Flights use a
    /// disjoint Philox domain and their own ordinal, including interrupted
    /// flights. Sharing path IDs across wavelengths intentionally correlates
    /// spectral estimates. Zero scattering consumes no random flight draw.
    #[allow(clippy::too_many_arguments)]
    pub fn trace(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_observed(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            |_| Ok(()),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
        mut observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_vertices_observed(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            |vertex| match vertex {
                CurvedRayleighVertex::Surface(path) => observe_surface(path),
                CurvedRayleighVertex::Volume { .. } => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_vertices_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
        mut observe_vertex: impl FnMut(CurvedRayleighVertex<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_events_observed(
            cx, beam, initial_medium, settings, seed, path_id,
            &[],
            |event| match event {
                CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                    observe_vertex(CurvedRayleighVertex::Surface(path))
                }
                CurvedRayleighEvent::VolumeArrival { path, phase } => {
                    observe_vertex(CurvedRayleighVertex::Volume { path, phase })
                }
                _ => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_events_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
        thermal_sources: &super::SpatialSources<'_>,
        observe_event: impl FnMut(CurvedRayleighEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        let draws = |ordinal, domain| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, ordinal, domain],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        };
        self.trace_with_tracking_samples_events_impl(
            cx,
            beam,
            initial_medium,
            settings,
            |event| draws(event, 0x706f_6c32),
            |flight| draws(flight, 0x706f_6c34),
            true,
            thermal_sources,
            observe_event,
        )
    }

    /// Deterministic custom sampler for the original uniform-density path.
    /// Heterogeneous density and positive local-field majorants need an
    /// independent fourth acceptance variate and refuse here. A local field's
    /// exact zero majorant consumes no flight sample and needs no acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_samples_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            |_| Ok(()),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_samples_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
        mut observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_samples_vertices_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            |vertex| match vertex {
                CurvedRayleighVertex::Surface(path) => observe_surface(path),
                CurvedRayleighVertex::Volume { .. } => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_samples_vertices_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
        mut observe_vertex: impl FnMut(CurvedRayleighVertex<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_samples_events_observed(
            cx, beam, initial_medium, settings, surface_samples, flight_samples,
            &[],
            |event| match event {
                CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                    observe_vertex(CurvedRayleighVertex::Surface(path))
                }
                CurvedRayleighEvent::VolumeArrival { path, phase } => {
                    observe_vertex(CurvedRayleighVertex::Volume { path, phase })
                }
                _ => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_samples_events_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut flight_samples: impl FnMut(u32) -> [f64; 3],
        thermal_sources: &super::SpatialSources<'_>,
        observe_event: impl FnMut(CurvedRayleighEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_events_impl(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            |flight| {
                let draw = flight_samples(flight);
                [draw[0], draw[1], draw[2], 0.0]
            },
            false,
            thermal_sources,
            observe_event,
        )
    }

    /// Deterministic custom sampler for uniform or heterogeneous density.
    /// Draws are [majorant CDF, phase cosine, phase azimuth, delta acceptance].
    #[allow(clippy::too_many_arguments)]
    pub fn trace_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            |_| Ok(()),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_tracking_samples_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        mut observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_vertices_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            |vertex| match vertex {
                CurvedRayleighVertex::Surface(path) => observe_surface(path),
                CurvedRayleighVertex::Volume { .. } => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_tracking_samples_vertices_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        observe_vertex: impl FnMut(CurvedRayleighVertex<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_impl(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            true,
            observe_vertex,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_with_tracking_samples_impl(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        heterogeneous_samples: bool,
        mut observe_vertex: impl FnMut(CurvedRayleighVertex<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_events_impl(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            heterogeneous_samples,
            &[],
            |event| match event {
                CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                    observe_vertex(CurvedRayleighVertex::Surface(path))
                }
                CurvedRayleighEvent::VolumeArrival { path, phase } => {
                    observe_vertex(CurvedRayleighVertex::Volume { path, phase })
                }
                CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure { .. })
                | CurvedRayleighEvent::VolumeDeparture { .. } => Ok(()),
            },
        )
    }

    /// Observe the original sampled directions and proposal densities without
    /// drawing a second phase or surface sample. Null candidates remain absent
    /// from the physical-vertex stream; all original flight limits still apply.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_with_tracking_samples_events_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        thermal_sources: &super::SpatialSources<'_>,
        observe_event: impl FnMut(CurvedRayleighEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_events_impl(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            true,
            thermal_sources,
            observe_event,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_with_tracking_samples_events_impl(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut flight_samples: impl FnMut(u32) -> [f64; 4],
        heterogeneous_samples: bool,
        thermal_sources: &super::SpatialSources<'_>,
        observe_event: impl FnMut(CurvedRayleighEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        cx.checkpoint()?;
        if settings.max_volume_events == 0 || settings.max_volume_events > MAX_POLARIZED_PATH_EVENTS
        {
            return Err(SceneError::InvalidSettings.into());
        }
        let mut free_flight_draws = 0_u32;
        // A path keeps its vacuum wavelength. Reuse the last spectral field's
        // majorant, including across null candidates, and thin against the
        // exact rate used to sample that flight. Density-scaled laws retain
        // their original evaluation and arithmetic. No per-medium workspace.
        let field_majorant: Cell<Option<(usize, f64)>> = Cell::new(None);
        // Native surface and molecular callbacks execute sequentially on this
        // stack. One observer retains cumulative source work in actual path
        // order without allocating a vertex history or retracing the camera.
        let observe_event = RefCell::new(observe_event);
        let (path, volume_events) = self.scene.trace_with_flights_events(
            cx,
            beam,
            initial_medium,
            settings.surfaces,
            surface_samples,
            self.radiance_brushing,
            thermal_sources,
            |event| (observe_event.borrow_mut())(CurvedRayleighEvent::Surface(event)),
            settings.max_volume_events,
            |medium, wavelength| {
                let Some(law) = self.media.get(medium) else {
                    return Ok(None);
                };
                if matches!(law, CurvedRayleighBinding::Density(law) if law.is_heterogeneous())
                    && !heterogeneous_samples
                {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let rate = if matches!(law, CurvedRayleighBinding::Field(_)) {
                    if let Some((cached_medium, rate)) = field_majorant.get()
                        && cached_medium == medium
                    {
                        rate
                    } else {
                        let rate = law.scattering_per_m(wavelength)?;
                        field_majorant.set(Some((medium, rate)));
                        rate
                    }
                } else {
                    law.scattering_per_m(wavelength)?
                };
                if rate == 0.0 {
                    return Ok(None);
                }
                if matches!(law, CurvedRayleighBinding::Field(_)) && !heterogeneous_samples {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let draws = flight_samples(free_flight_draws);
                free_flight_draws = free_flight_draws
                    .checked_add(1)
                    .ok_or(CurvedSceneError::NumericalRange)?;
                if !draws
                    .into_iter()
                    .all(|u| u.is_finite() && (0.0..1.0).contains(&u))
                {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                // Compensate the rounded subtraction in 1-u. If it rounds to
                // one, u itself is the correctly resolved first-order depth.
                let complement = 1.0 - draws[0];
                let depth = if complement == 1.0 {
                    draws[0]
                } else {
                    -det::ln(complement) * (draws[0] / (1.0 - complement))
                };
                let distance_m = depth / rate;
                if distance_m.is_nan() || distance_m < 0.0 || (depth > 0.0 && distance_m == 0.0) {
                    return Err(volume_error(VolumeTransportError::Unrepresentable));
                }
                Ok(Some(CurvedVolumeFlight {
                    distance_m,
                    phase_samples: [draws[1], draws[2]],
                    acceptance_sample: draws[3],
                }))
            },
            |medium, path, samples, acceptance_sample| {
                let law = self
                    .media
                    .get(medium)
                    .ok_or(CurvedSceneError::NumericalRange)?;
                let sampled_majorant = field_majorant
                    .get()
                    .filter(|(cached_medium, _)| *cached_medium == medium)
                    .map(|(_, rate)| rate);
                law.scatter_candidate(
                    *path.beam(),
                    samples,
                    acceptance_sample,
                    sampled_majorant,
                    |phase| {
                        (observe_event.borrow_mut())(CurvedRayleighEvent::VolumeArrival {
                            path,
                            phase,
                        })
                    },
                    |phase, next| {
                        let mut departure = *path;
                        departure.beam = next.beam;
                        (observe_event.borrow_mut())(CurvedRayleighEvent::VolumeDeparture {
                            path: &departure,
                            incoming_beam: *path.beam(),
                            phase,
                            pdf: next.pdf,
                        })
                    },
                )
            },
        )?;
        cx.checkpoint()?;
        Ok(CurvedRayleighPath {
            path,
            volume_events,
            free_flight_draws,
        })
    }
}

fn scattering_error(error: PolarizedScatteringError) -> CurvedSceneError {
    SceneError::from(error).into()
}
fn volume_error(error: VolumeTransportError) -> CurvedSceneError {
    scattering_error(PolarizedScatteringError::Volume(error))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod spectral_field_tests;

#[cfg(test)]
mod grid_binding_tests {
    use super::*;

    #[test]
    fn zero_copy_grid_derives_majorant_identity_and_null_collision_probability() {
        let cells = [0.5, 1.5];
        let grid = VolumeGrid::new([1, 1, 2], &cells, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]);
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let medium =
            CurvedRayleighMedium::try_new_volume_grid(&spectrum, &grid, &phase, 2).unwrap();
        assert_eq!(medium.density_kg_m3(), 1.5);
        assert!(medium.is_heterogeneous());
        assert_ne!(medium.density_identity(), Some(ContentHash([0; 32])));
        assert!(
            !medium
                .accepts_candidate(Point3::new(0.0, 0.0, 0.25), 0.8)
                .unwrap()
        );
        assert!(
            medium
                .accepts_candidate(Point3::new(0.0, 0.0, 0.75), 0.9)
                .unwrap()
        );
        assert!(
            !medium
                .accepts_candidate(Point3::new(0.0, 0.0, 2.0), 0.0)
                .unwrap()
        );

        let changed = [0.5, 1.4];
        let changed_grid = VolumeGrid::new([1, 1, 2], &changed, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]);
        let changed_medium =
            CurvedRayleighMedium::try_new_volume_grid(&spectrum, &changed_grid, &phase, 2).unwrap();
        assert_ne!(medium.density_identity(), changed_medium.density_identity());

        assert!(CurvedRayleighMedium::try_new_volume_grid(&spectrum, &grid, &phase, 1).is_err());
    }
}
