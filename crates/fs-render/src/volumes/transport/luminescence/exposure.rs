//! Positive temporal quadrature of complete accepted luminescence observations.
//!
//! Every time node transports its own powers and opacity before light is
//! averaged. This owner estimates a declared quadrature; it neither evolves
//! the physical state nor certifies the caller's temporal resolution.

use super::super::camera::{CameraSample, Compensated, approximate_xyz, render_camera_radiance};
use super::super::{VolumeTransportError, checkpoint, positive_product};
use super::{
    LuminescenceCameraError, LuminescenceCameraSettings, LuminescenceGrid,
    LuminescenceScatteringBudget, ScatteringLuminescenceGrid,
};
use crate::camera::PhysicalCamera;
use crate::volumes::DvrError;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// One accepted field, with its original native transport owner.
#[derive(Clone, Copy)]
pub enum LuminescenceExposureField<'a, 's, 'm, 'g, const M: usize> {
    /// Exact cell source/absorption integrals without elastic scattering.
    Direct(&'a LuminescenceGrid<'m, 'g, M>),
    /// Ballistic and multiply scattered emission on the native path walk.
    Scattered {
        /// The complete accepted emission, opacity and scattering snapshot.
        field: &'a ScatteringLuminescenceGrid<'s, 'm, 'g, M>,
        /// Independent limits for each completed temporal path history.
        budget: LuminescenceScatteringBudget,
    },
}

/// The caller binds an accepted physical snapshot to its actual time.
/// Positive dimensionless weights describe a normalized temporal quadrature;
/// they must sum to one within 64 machine epsilons. Only this rounding-sized
/// normalization is removed by the compensated weighted mean.
#[derive(Clone, Copy)]
pub struct LuminescenceExposureNode<'a, 's, 'm, 'g, const M: usize> {
    /// Caller-supplied observation time in seconds within the declared bounds.
    pub time_s: f64,
    /// Positive dimensionless quadrature weight, normalized across all nodes.
    pub weight: f64,
    /// Powers and opacity belonging to this accepted observation time.
    pub field: LuminescenceExposureField<'a, 's, 'm, 'g, M>,
}

/// Temporal work in addition to the original camera/image/native ray limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LuminescenceExposureBudget {
    /// Maximum supplied temporal nodes, including zero-emission snapshots.
    pub max_nodes: usize,
    /// Maximum pixels * primary samples per pixel * temporal nodes.
    pub max_observations: u64,
}

/// Native observation refusals and temporal admission diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LuminescenceExposureError {
    /// Original physical-camera, allocation, transport or observer refusal.
    Camera(LuminescenceCameraError),
    /// Nonfinite/collapsed bounds or unordered/out-of-window node times.
    Time,
    /// Empty nodes, nonpositive weights, or weights not summing to one.
    Quadrature,
    /// Node source identities or exact wavelength arrays differ.
    Correspondence,
    /// Temporal node/evaluation/allocation limits are exceeded.
    Budget,
}

impl From<LuminescenceCameraError> for LuminescenceExposureError {
    fn from(error: LuminescenceCameraError) -> Self { Self::Camera(error) }
}
impl From<VolumeTransportError> for LuminescenceExposureError {
    fn from(error: VolumeTransportError) -> Self { Self::Camera(error.into()) }
}
impl core::fmt::Display for LuminescenceExposureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "luminescence exposure refused: {self:?}")
    }
}
impl core::error::Error for LuminescenceExposureError {}

/// Complete normalized spectral observation over one declared exposure.
/// The camera is retained explicitly: matching pixel counts alone never
/// establishes correspondence between independently rendered images.
#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceCameraExposure<const M: usize> {
    camera: PhysicalCamera,
    settings: LuminescenceCameraSettings,
    interval_s: [f64; 2],
    wavelengths_nm: [f64; M],
    source_identity: ContentHash,
    field_history_identity: ContentHash,
    radiance: Vec<[f64; M]>,
    cell_visits: u64,
}

impl<const M: usize> LuminescenceCameraExposure<M> {
    /// One frozen camera shared by every temporal and spatial observation.
    #[must_use]
    pub const fn camera(&self) -> &PhysicalCamera { &self.camera }
    /// Original native raster/pupil sampling settings and whole-image limits.
    #[must_use]
    pub const fn settings(&self) -> LuminescenceCameraSettings { self.settings }
    /// Declared increasing finite exposure bounds in seconds.
    #[must_use]
    pub const fn interval_s(&self) -> [f64; 2] { self.interval_s }
    /// Finite positive duration used to convert mean radiance to its integral.
    #[must_use]
    pub fn duration_s(&self) -> f64 { self.interval_s[1] - self.interval_s[0] }
    /// Original discrete spectral-node/bin wavelengths in vacuum nanometres.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] { &self.wavelengths_nm }
    /// Shared declared photophysical source and condition identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash { self.source_identity }
    /// Ordered times, weights and complete accepted field identities. The
    /// camera and spatial sampling settings are retained separately above.
    #[must_use]
    pub const fn field_history_identity(&self) -> ContentHash { self.field_history_identity }
    /// Top-row-first row-major, bin-integrated mean radiance in W/(m² sr).
    #[must_use]
    pub fn mean_radiance_w_m2_sr(&self) -> &[[f64; M]] { &self.radiance }
    /// Actual DDA visits over every time node and native primary/path sample.
    #[must_use]
    pub const fn cell_visits(&self) -> u64 { self.cell_visits }

    /// Time-integrated radiance in J/(m² sr). This is not sensor irradiance,
    /// aperture gain or calibrated display brightness. A nonrepresentable
    /// integral refuses while the already computed finite mean remains usable.
    pub fn integrated_radiance_j_m2_sr(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; M]>, LuminescenceExposureError> {
        checkpoint(cx)?;
        let mut result = Vec::new();
        result.try_reserve_exact(self.radiance.len()).map_err(|_| {
            LuminescenceCameraError::from(DvrError::ImageAllocationRefused { pixels: self.radiance.len() })
        })?;
        for pixel in &self.radiance {
            checkpoint(cx)?;
            let integrated = pixel.map(|value| positive_product(value, self.duration_s(), 1.0));
            if integrated.iter().zip(pixel).any(|(value, source)| {
                !value.is_finite() || *value < 0.0 || (*source > 0.0 && *value == 0.0)
            }) {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
            result.push(integrated);
        }
        Ok(result)
    }

    /// The existing relative visible observer applied to mean bin-integrated
    /// radiance, without an extra wavelength-width or opacity factor.
    pub fn approximate_mean_xyz(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; 3]>, LuminescenceExposureError> {
        Ok(approximate_xyz(cx, &self.wavelengths_nm, &[1.0; M], self.radiance.len(),
            |pixel, lane| (self.radiance[pixel][lane], 1.0, 0.0, 0.0))?)
    }
}

impl<const M: usize> LuminescenceExposureField<'_, '_, '_, '_, M> {
    fn wavelengths(&self) -> &[f64; M] {
        match self { Self::Direct(field) => field.wavelengths_nm(), Self::Scattered { field, .. } => field.wavelengths_nm() }
    }
    fn source(&self) -> ContentHash {
        match self { Self::Direct(field) => field.source_identity(), Self::Scattered { field, .. } => field.source_identity() }
    }
    fn identity(&self) -> ContentHash {
        match self { Self::Direct(field) => field.content_identity(), Self::Scattered { field, .. } => field.content_identity() }
    }

    fn trace(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        remaining_cells: u64,
        seed: u64,
        path_id: u64,
    ) -> Result<([f64; M], u32), LuminescenceCameraError> {
        match self {
            Self::Direct(field) => {
                let sample = field.transport_counted(cx, origin, direction, f64::INFINITY, remaining_cells)?;
                Ok((sample.transport.emitted_radiance_w_m2_sr, sample.cell_visits))
            }
            Self::Scattered { field, budget } => {
                let sample = field.sample_ray(cx, origin, direction, LuminescenceScatteringBudget {
                    max_cell_visits: budget.max_cell_visits.min(u32::try_from(remaining_cells).unwrap_or(u32::MAX)),
                    ..*budget
                }, seed, path_id)?;
                Ok((sample.radiance_w_m2_sr, sample.cell_visits))
            }
        }
    }
}

/// Observe changing accepted powers AND opacity through one fixed physical
/// camera. Nodes must be strictly time-ordered within finite increasing bounds,
/// with one original source and exact spectral-node order. Each node's time and
/// accepted physical state are supplied by the caller; this function does not
/// construct or validate a kinetic/thermal history.
///
/// Complete native light transport precedes temporal averaging. The camera's
/// original pixel/pupil rays and path streams are shared across time nodes;
/// correlated scattering samples remain Monte Carlo estimates. Constant fields
/// therefore retain the original frozen observation. Positive quadrature gives
/// no temporal error certificate: resolve pump edges and changing populations,
/// compare temporal refinements, and retain each original material domain.
///
/// Every temporal field query is admitted before tracing, and every actual DDA
/// visit consumes the same whole-image cell budget. The native per-history
/// scattering and ray limits remain binding. No partial image is returned on
/// cancellation or refusal, and no physical field is modified. Backgrounds,
/// camera motion, flight delays and exposure/sensor gain remain separate.
pub fn render_luminescence_exposure<const M: usize>(
    cx: &Cx<'_>,
    nodes: &[LuminescenceExposureNode<'_, '_, '_, '_, M>],
    camera: &PhysicalCamera,
    settings: LuminescenceCameraSettings,
    interval_s: [f64; 2],
    budget: LuminescenceExposureBudget,
) -> Result<LuminescenceCameraExposure<M>, LuminescenceExposureError> {
    checkpoint(cx)?;
    let duration = interval_s[1] - interval_s[0];
    if !interval_s.iter().all(|value| value.is_finite()) || !duration.is_finite() || duration <= 0.0 {
        return Err(LuminescenceExposureError::Time);
    }
    if nodes.is_empty() { return Err(LuminescenceExposureError::Quadrature); }
    if nodes.len() > budget.max_nodes { return Err(LuminescenceExposureError::Budget); }
    let observations = settings.resolution[0].checked_mul(settings.resolution[1])
        .and_then(|pixels| u64::try_from(pixels).ok())
        .and_then(|pixels| pixels.checked_mul(u64::from(settings.samples_per_pixel)))
        .and_then(|samples| u64::try_from(nodes.len()).ok().and_then(|nodes| samples.checked_mul(nodes)))
        .ok_or(LuminescenceExposureError::Budget)?;
    if observations > budget.max_observations { return Err(LuminescenceExposureError::Budget); }
    let wavelengths_nm = *nodes[0].field.wavelengths();
    let source_identity = nodes[0].field.source();
    let mut hash = DomainHasher::new("org.frankensim.render.luminescence-exposure-fields.v1");
    hash.update(source_identity.as_bytes());
    for value in interval_s.into_iter().chain(wavelengths_nm) {
        hash.update(&(if value == 0.0 { 0.0 } else { value }).to_bits().to_le_bytes());
    }
    hash.update(&(nodes.len() as u64).to_le_bytes());
    let mut fractions = Vec::new();
    fractions.try_reserve_exact(nodes.len()).map_err(|_| LuminescenceExposureError::Budget)?;
    let mut weight = Compensated::default();
    let mut previous_time = f64::NEG_INFINITY;
    for node in nodes {
        checkpoint(cx)?;
        if !node.time_s.is_finite() || node.time_s <= previous_time
            || !(interval_s[0]..=interval_s[1]).contains(&node.time_s) {
            return Err(LuminescenceExposureError::Time);
        }
        if !node.weight.is_finite() || node.weight <= 0.0 || node.weight > 1.0 {
            return Err(LuminescenceExposureError::Quadrature);
        }
        if node.field.source() != source_identity || node.field.wavelengths() != &wavelengths_nm {
            return Err(LuminescenceExposureError::Correspondence);
        }
        let previous_weight = weight.value;
        weight.add(node.weight);
        if !weight.value.is_finite() || weight.value <= 0.0
            || weight.value > 1.0 + 64.0 * f64::EPSILON {
            return Err(LuminescenceExposureError::Quadrature);
        }
        fractions.push((previous_weight / weight.value, node.weight / weight.value));
        previous_time = node.time_s;
        hash.update(&(if node.time_s == 0.0 { 0.0 } else { node.time_s }).to_bits().to_le_bytes());
        hash.update(&node.weight.to_bits().to_le_bytes());
        hash.update(node.field.identity().as_bytes());
    }
    if !weight.value.is_finite() || (weight.value - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(LuminescenceExposureError::Quadrature);
    }
    let mut path_id = 0_u64;
    let mut cell_visits = 0_u64;
    let radiance = render_camera_radiance(cx, camera, settings, 0x6c75_6d63,
        |origin, direction, remaining_cells| {
            let mut means = [Compensated::default(); M];
            let mut visits = 0_u32;
            for (node, &(old_fraction, new_fraction)) in nodes.iter().zip(&fractions) {
                checkpoint(cx)?;
                let remaining = remaining_cells.checked_sub(u64::from(visits))
                    .ok_or(VolumeTransportError::Budget)?;
                let (radiance, used) = node.field.trace(cx, origin, direction, remaining, settings.seed, path_id)?;
                visits = visits.checked_add(used).ok_or(VolumeTransportError::Budget)?;
                if u64::from(visits) > remaining_cells { return Err(VolumeTransportError::Budget.into()); }
                for lane in 0..M {
                    means[lane].update_weighted_mean(radiance[lane], old_fraction, new_fraction);
                }
            }
            path_id = path_id.checked_add(1).ok_or(VolumeTransportError::Budget)?;
            cell_visits = cell_visits.checked_add(u64::from(visits)).ok_or(VolumeTransportError::Budget)?;
            Ok(CameraSample { emission: means.map(|mean| mean.value), transmittance: [1.0; M],
                optical_depth: [0.0; M], cell_visits: visits })
        })?;
    Ok(LuminescenceCameraExposure { camera: camera.clone(), settings, interval_s,
        wavelengths_nm, source_identity, field_history_identity: hash.finalize(), radiance, cell_visits })
}

#[cfg(test)]
mod tests;
