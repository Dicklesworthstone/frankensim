//! Original material laws in native cell transport and physical-camera images.
//!
//! Preparation evaluates each source once per requested wavelength and material
//! model. The immutable packet has no off-node query or interpolation API.
//! Geometry and density stay borrowed, and models share the original cell map.

use super::{MaterialAbsorptionBinding, MaterialAbsorptionError};
use crate::camera::PhysicalCamera;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use crate::volumes::MAX_OPACITY_WAVELENGTHS;
use crate::volumes::transport::camera::{
    AbsorbingCameraImage, AbsorbingCameraRadianceImage, VolumeCameraBackground, VolumeCameraError,
    VolumeCameraSettings, render_absorbing_field, render_absorbing_field_with_background,
};
use crate::volumes::transport::{
    AbsorbingGrid, CellPalette, GridSegmentTransport, SegmentTransport, VolumeTransportError,
};
use core::marker::PhantomData;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

mod deposition;
mod mixture;
pub use mixture::{
    MaterialAbsorptionComponent, MaterialAbsorptionMixture, MaterialAbsorptionMixtureField,
};

/// Explicit source choice at each original material condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialAbsorptionEmission {
    /// Absorption only. No material temperature is required or evaluated.
    None,
    /// Each absorbing material emits at its own accepted absolute temperature.
    /// No average temperature, independent emissivity or state evolution is
    /// inferred. Mixture constituents need not share a temperature; each native
    /// Planck source remains linked to that same original absorber.
    AcceptedLte,
}

/// Preparation storage and source-evaluation limits, checked before allocation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MaterialAbsorptionPreparationBudget {
    /// Maximum prepared models, including unused palette entries.
    pub max_models: usize,
    /// Maximum constituent/wavelength pairs, including zero-weight mixture
    /// entries. Each pair evaluates the original absorption law and, in
    /// accepted LTE mode, at most one Planck value.
    pub max_spectral_evaluations: u64,
}

#[derive(Clone, Copy)]
struct PreparedMaterial<const N: usize> {
    absorption: [f64; N],
    source: Option<[f64; N]>,
}

/// Exact source-law point packets on an original physical density grid.
///
/// This follows the native prepared-field convention: `geometry` supplies its
/// admitted metre grid, kg/m³ density, numerical identity and traversal ceiling.
/// Its existing absorption law and temperature sources are not used. All
/// absorption and optional LTE source state come from the material bindings.
/// Intrinsic material density already normalizes each mass coefficient; the
/// borrowed grid's density remains its declared reference mass per bulk m³:
/// the selected material for a single binding, or the recipe's explicit basis.
/// The caller owns that continuum model; no interface reflection, particle
/// scattering, partial-coverage law or density/temperature extrapolation follows.
///
/// The original bindings and immutable cell map remain retained. Preparation
/// storage is O(models * N), with no dense per-cell spectral or thermal copy,
/// and tracing evaluates neither material sources nor Planck functions. Every
/// transport and camera method uses exactly the stored nodes. A new wavelength
/// packet or accepted material state requires a new preparation.
pub struct MaterialAbsorptionField<
    'g,
    'field,
    'material,
    const N: usize,
    Model = MaterialAbsorptionBinding<'material>,
> {
    geometry: &'g AbsorbingGrid<'field>,
    materials: CellPalette<Model>,
    wavelengths_nm: [f64; N],
    source_bounds_nm: [f64; 2],
    models: Vec<PreparedMaterial<N>>,
    emission: MaterialAbsorptionEmission,
    vacuum: bool,
    has_emission: bool,
    identity: ContentHash,
    _material_lifetime: PhantomData<&'material ()>,
}

impl<'g, 'field, 'material, const N: usize> MaterialAbsorptionField<'g, 'field, 'material, N> {
    /// Prepare 1--16 increasing, positive wavelength nodes from each original
    /// material law, within every supplied source's original closed support.
    /// Different source bands may overlap; they are never widened or resampled.
    /// Unused models are validated but cannot create absorption or emission.
    pub fn try_new(
        cx: &Cx<'_>,
        geometry: &'g AbsorbingGrid<'field>,
        materials: CellPalette<MaterialAbsorptionBinding<'material>>,
        wavelengths_nm: [f64; N],
        emission: MaterialAbsorptionEmission,
        budget: MaterialAbsorptionPreparationBudget,
    ) -> Result<Self, MaterialAbsorptionFieldError> {
        Self::prepare(
            cx,
            geometry,
            materials,
            wavelengths_nm,
            emission,
            budget,
            "org.frankensim.render.material-absorption-field.v1",
            |_| 1,
            |material, index| {
                let model = prepare_material(material, &wavelengths_nm, emission, index)?;
                Ok((model, material.wavelength_bounds_nm(), material.identity()))
            },
        )
    }
}

impl<'g, 'field, 'material, const N: usize, Model>
    MaterialAbsorptionField<'g, 'field, 'material, N, Model>
{
    // One prepared packet, map and observer owner serves individual materials
    // and constituent mixtures. The source-specific closures run only here.
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        cx: &Cx<'_>,
        geometry: &'g AbsorbingGrid<'field>,
        materials: CellPalette<Model>,
        wavelengths_nm: [f64; N],
        emission: MaterialAbsorptionEmission,
        budget: MaterialAbsorptionPreparationBudget,
        identity_domain: &'static str,
        component_count: impl Fn(&Model) -> usize,
        mut prepare_model: impl FnMut(
            &Model,
            usize,
        ) -> Result<
            (PreparedMaterial<N>, [f64; 2], ContentHash),
            MaterialAbsorptionFieldError,
        >,
    ) -> Result<Self, MaterialAbsorptionFieldError> {
        checkpoint(cx)?;
        let model_count = materials.values().len();
        if model_count > budget.max_models {
            return Err(VolumeTransportError::Budget.into());
        }
        let mut evaluations = 0_u64;
        for model in materials.values() {
            checkpoint(cx)?;
            evaluations = u64::try_from(component_count(model))
                .ok()
                .and_then(|count| count.checked_mul(u64::try_from(N).ok()?))
                .and_then(|count| evaluations.checked_add(count))
                .ok_or(VolumeTransportError::Budget)?;
            if evaluations > budget.max_spectral_evaluations {
                return Err(VolumeTransportError::Budget.into());
            }
        }
        if materials.cell_indices().len() != geometry.cell_count() {
            return Err(VolumeTransportError::Palette.into());
        }
        if N == 0
            || N > MAX_OPACITY_WAVELENGTHS
            || wavelengths_nm.iter().any(|w| !w.is_finite() || *w <= 0.0)
            || wavelengths_nm.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(MaterialAbsorptionError::InvalidWavelengthSamples.into());
        }
        let mut models = Vec::new();
        models
            .try_reserve_exact(model_count)
            .map_err(|_| VolumeTransportError::Budget)?;
        let mut hash = DomainHasher::new(identity_domain);
        hash.update(geometry.content_identity().as_bytes());
        hash.update(&(model_count as u64).to_le_bytes());
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&[u8::from(
            emission == MaterialAbsorptionEmission::AcceptedLte,
        )]);
        for wavelength in wavelengths_nm {
            hash.update(&wavelength.to_bits().to_le_bytes());
        }
        let mut source_bounds_nm = [0.0_f64, f64::INFINITY];
        for (index, material) in materials.values().iter().enumerate() {
            checkpoint(cx)?;
            let (model, bounds, source_identity) = prepare_model(material, index)?;
            source_bounds_nm[0] = source_bounds_nm[0].max(bounds[0]);
            source_bounds_nm[1] = source_bounds_nm[1].min(bounds[1]);
            hash.update(source_identity.as_bytes());
            for coefficient in model.absorption {
                hash.update(&coefficient.to_bits().to_le_bytes());
            }
            hash.update(&[u8::from(model.source.is_some())]);
            if let Some(source) = model.source {
                for value in source {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            models.push(model);
        }
        let mut vacuum = true;
        let mut has_emission = false;
        for (cell, &index) in materials.cell_indices().iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            hash.update(&index.to_le_bytes());
            if geometry
                .cell_density_kg_m3(cell)
                .expect("admitted cell map")
                > 0.0
            {
                let model = &models[index as usize];
                if model.absorption.iter().any(|value| *value > 0.0) {
                    vacuum = false;
                    has_emission |= model.source.is_some();
                }
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            geometry,
            materials,
            wavelengths_nm,
            source_bounds_nm,
            models,
            emission,
            vacuum,
            has_emission,
            identity: hash.finalize(),
            _material_lifetime: PhantomData,
        })
    }

    /// Original immutable source models or recipes, with the native cell map.
    #[must_use]
    pub const fn materials(&self) -> &CellPalette<Model> {
        &self.materials
    }

    /// The only wavelengths observed by this prepared field, in vacuum nm.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; N] {
        &self.wavelengths_nm
    }

    /// Intersection of original source supports. This records source scope;
    /// the prepared field still permits only its immutable wavelength nodes.
    #[must_use]
    pub const fn source_wavelength_bounds_nm(&self) -> [f64; 2] {
        self.source_bounds_nm
    }

    /// Accepted geometry, density, material identities, cell map, node values
    /// and explicitly chosen emission mode.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Original explicit source choice; no thermal mode is inferred.
    #[must_use]
    pub const fn emission_mode(&self) -> MaterialAbsorptionEmission {
        self.emission
    }

    /// No mapped positive-density cell absorbs at any prepared node.
    #[must_use]
    pub const fn is_vacuum(&self) -> bool {
        self.vacuum
    }

    /// Exact Beer--Lambert transmission at the stored nodes. No source is read.
    pub fn transmittance(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
    ) -> Result<[f64; N], VolumeTransportError> {
        Ok(self
            .integrate::<false, false>(cx, origin, direction, max_distance_m, u32::MAX)?
            .transport
            .transmittance)
    }

    /// Exact native cell transport at the stored nodes, including the explicitly
    /// prepared LTE sources. Radiance is W/(m² sr nm). Positive infinity denotes
    /// a distant endpoint; a finite distance excludes everything behind it.
    pub fn transport(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
    ) -> Result<SegmentTransport<N>, VolumeTransportError> {
        Ok(self
            .integrate::<true, false>(cx, origin, direction, max_distance_m, u32::MAX)?
            .transport)
    }

    /// Physical pixel/pupil observation with the native camera's exact sampling,
    /// stable thin/opaque limits and whole-image cell budget. These normalized
    /// radiances include no exposure, scattering, interface loss or refraction.
    pub fn render_camera(
        &self,
        cx: &Cx<'_>,
        camera: &PhysicalCamera,
        settings: VolumeCameraSettings,
    ) -> Result<AbsorbingCameraImage<N>, VolumeCameraError> {
        render_absorbing_field(
            cx,
            camera,
            self.wavelengths_nm,
            self.identity,
            settings,
            |origin, direction, distance, cells| {
                self.integrate::<true, true>(cx, origin, direction, distance, cells)
            },
        )
    }

    /// Compose each actual camera ray's thermal source and attenuated background
    /// before pixel/pupil averaging, preserving illumination/opacity correlation.
    /// The callback's radiance is W/(m² sr nm) at exactly the supplied nodes.
    /// Finite source distances, cancellation and refusal retain native behavior.
    pub fn render_camera_with_background(
        &self,
        cx: &Cx<'_>,
        camera: &PhysicalCamera,
        settings: VolumeCameraSettings,
        background: impl FnMut(
            &Cx<'_>,
            Point3,
            Vec3,
            &[f64; N],
        ) -> Result<VolumeCameraBackground<N>, VolumeCameraError>,
    ) -> Result<AbsorbingCameraRadianceImage<N>, VolumeCameraError> {
        render_absorbing_field_with_background(
            cx,
            camera,
            self.wavelengths_nm,
            self.identity,
            settings,
            |origin, direction, distance, cells| {
                self.integrate::<true, true>(cx, origin, direction, distance, cells)
            },
            background,
        )
    }

    fn integrate<const EMIT: bool, const COUNTED: bool>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        distance: f64,
        cells: u32,
    ) -> Result<GridSegmentTransport<N>, VolumeTransportError> {
        let mut previous_model = None;
        self.geometry.integrate_prepared::<N, EMIT, COUNTED>(
            cx,
            origin,
            direction,
            distance,
            cells,
            self.vacuum,
            self.has_emission,
            |cell, coefficients| {
                let model = self.materials.cell_indices()[cell];
                if previous_model != Some(model) {
                    *coefficients = self.models[model as usize].absorption;
                    previous_model = Some(model);
                }
                Ok(())
            },
            |cell, lane| {
                Ok(self.models[self.materials.cell_indices()[cell] as usize]
                    .source
                    .as_ref()
                    .map(|source| source[lane]))
            },
        )
    }
}

fn prepare_material<const N: usize>(
    material: &MaterialAbsorptionBinding<'_>,
    wavelengths_nm: &[f64; N],
    emission: MaterialAbsorptionEmission,
    index: usize,
) -> Result<PreparedMaterial<N>, MaterialAbsorptionFieldError> {
    let samples = material.sample(wavelengths_nm)?;
    let mut absorption = [0.0; N];
    absorption.copy_from_slice(samples.mass_absorption_m2_kg());
    let source = if emission == MaterialAbsorptionEmission::AcceptedLte
        && absorption.iter().any(|value| *value > 0.0)
    {
        let temperature = material
            .accepted_temperature_k()
            .ok_or(MaterialAbsorptionFieldError::Temperature { model: index })?;
        let mut source = [0.0; N];
        for lane in 0..N {
            if absorption[lane] > 0.0 {
                source[lane] = planck_radiance_w_m2_sr_nm(temperature, wavelengths_nm[lane])
                    .map_err(|_| VolumeTransportError::Unrepresentable)?;
            }
        }
        Some(source)
    } else {
        None
    };
    Ok(PreparedMaterial { absorption, source })
}

/// Preparation publishes no partial field on any original source or native
/// admission refusal. Transport and camera calls retain their native errors.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialAbsorptionFieldError {
    /// Original material law, source support or wavelength-node refusal.
    Absorption(MaterialAbsorptionError),
    /// Geometry/map, budget, arithmetic or cancellation refusal.
    Transport(VolumeTransportError),
    /// An absorbing model lacks one positive accepted Kelvin temperature.
    Temperature { model: usize },
    /// An active mixture constituent lacks one positive accepted temperature.
    ComponentTemperature { model: usize, component: usize },
}

impl From<MaterialAbsorptionError> for MaterialAbsorptionFieldError {
    fn from(error: MaterialAbsorptionError) -> Self {
        Self::Absorption(error)
    }
}
impl From<VolumeTransportError> for MaterialAbsorptionFieldError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for MaterialAbsorptionFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native material absorption preparation refused: {self:?}"
        )
    }
}
impl core::error::Error for MaterialAbsorptionFieldError {}

fn checkpoint(cx: &Cx<'_>) -> Result<(), VolumeTransportError> {
    cx.checkpoint().map_err(|_| VolumeTransportError::Cancelled)
}

#[cfg(test)]
mod tests;
