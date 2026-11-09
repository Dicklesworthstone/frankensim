//! Photoelectron emission driven by the original film's absorbed optical power.
//!
//! The original coherent depth operator retains standing-wave interference and
//! longitudinal TM absorption. Its entire-film loss, multiplied by physical
//! area once, supplies an **internal** emitted-electron yield. That yield must
//! be calibrated for this exact film thickness and selected exposed face; it
//! includes electron transport and escape losses. No escape depth or external
//! quantum efficiency is inferred from a bulk material name.
//!
//! [`EscapeFace`] declares a vacuum escape boundary. The electromagnetic owner
//! proves only that the corresponding optical exterior is lossless and supports
//! its propagating mode; refractive index alone cannot identify a vacuum.
//! Sources retain a fixed monochromatic frequency and complete material state.
//! Electron collection, contact potentials, space charge, Schottky lowering,
//! dynamic charging and feedback on the optical field remain outside this law.
//! The returned residual heat replaces this film's total optical deposition;
//! a consumer must not also deposit the original absorbed power as heat.

use super::{IlluminationBudget, PreparedLayerDepthInterval, PreparedLayerIllumination};
use crate::plane_wave::ReflectionError;
use fs_material::photoelectric::source::{PhotoelectricMaterial, PhotoelectricSourceError};
use fs_material::photoelectric::{PhotoelectricError, PhotoelectricLaw, PhotoelectricResponse};
use fs_math::c64::C64;

/// Explicit declaration of which external face emits into vacuum. The selected
/// source yield describes escape through this face, not the sum of both faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EscapeFace {
    /// The first finite film's entrance face borders a caller-declared vacuum.
    EntranceVacuum,
    /// The last finite film's exit face borders a caller-declared vacuum.
    ExitVacuum,
}

#[derive(Debug)]
pub enum LayerPhotoelectricError {
    Invalid(&'static str),
    Optical(ReflectionError),
    Source(PhotoelectricSourceError),
    Emission(PhotoelectricError),
}

impl core::fmt::Display for LayerPhotoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photoelectric film: {self:?}")
    }
}
impl std::error::Error for LayerPhotoelectricError {}
impl From<ReflectionError> for LayerPhotoelectricError {
    fn from(value: ReflectionError) -> Self {
        Self::Optical(value)
    }
}
impl From<PhotoelectricSourceError> for LayerPhotoelectricError {
    fn from(value: PhotoelectricSourceError) -> Self {
        Self::Source(value)
    }
}
impl From<PhotoelectricError> for LayerPhotoelectricError {
    fn from(value: PhotoelectricError) -> Self {
        Self::Emission(value)
    }
}

/// Immutable whole-film absorption and original photoemission law. Preparation
/// checks the complete sources once and reuses the original layer loss modes.
/// Each subsequent coherent evaluation costs O(1) and allocates nothing.
/// The original optical stack has already prepared its bounded layer sweeps.
pub struct PreparedLayerPhotoelectric<'p, 's, 'm, 'd> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    absorption: PreparedLayerDepthInterval<'p, 's, 'm>,
    material: &'d PhotoelectricMaterial,
    law: PhotoelectricLaw,
    layer_index: usize,
    face: EscapeFace,
    area_m2: f64,
    vacuum_wavelength_m: f64,
}

/// Surface emission, energy partition and original physical support. Keeping
/// the original source references prevents an unattached watt value from
/// masquerading as a source- and geometry-bound thermal/electrical observation.
pub struct LayerPhotoelectricResponse<'p, 's, 'm, 'd> {
    illumination: &'p PreparedLayerIllumination<'s, 'm>,
    material: &'d PhotoelectricMaterial,
    layer_index: usize,
    face: EscapeFace,
    area_m2: f64,
    vacuum_wavelength_m: f64,
    absorbed_power_w_m2: f64,
    emission: PhotoelectricResponse,
}

impl<'p, 's, 'm, 'd> PreparedLayerPhotoelectric<'p, 's, 'm, 'd> {
    /// Attach the calibrated emitting film to its original complete optical
    /// material snapshot. Buried films cannot emit through another material.
    /// `area_m2` is physical film area, so no further incidence cosine applies
    /// to the original normal Poynting power density.
    pub fn try_new(
        illumination: &'p PreparedLayerIllumination<'s, 'm>,
        layer_index: usize,
        material: &'d PhotoelectricMaterial,
        face: EscapeFace,
        area_m2: f64,
    ) -> Result<Self, LayerPhotoelectricError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(LayerPhotoelectricError::Invalid(
                "positive finite physical film area required",
            ));
        }
        let stack = illumination.stack();
        let layer = stack
            .layer(layer_index)
            .ok_or(LayerPhotoelectricError::Invalid(
                "emitting film index is out of range",
            ))?;
        let exposed = match face {
            EscapeFace::EntranceVacuum => layer_index == 0,
            EscapeFace::ExitVacuum => layer_index + 1 == stack.layer_count(),
        };
        if !exposed {
            return Err(LayerPhotoelectricError::Invalid(
                "selected film has no exposed face on the declared vacuum boundary",
            ));
        }
        if material.resolved() != layer.response.material() {
            return Err(LayerPhotoelectricError::Invalid(
                "emission and optical loss require the same complete material snapshot",
            ));
        }
        if layer.thickness_m != material.thickness_m() {
            return Err(LayerPhotoelectricError::Invalid(
                "film thickness differs from the original photoemission calibration",
            ));
        }
        // Photon energy follows the vacuum frequency, regardless of either
        // exterior refractive index or the internal longitudinal wave number.
        let vacuum_wavelength_m =
            (core::f64::consts::TAU * 299_792_458.0) / stack.angular_frequency_rad_s();
        if !vacuum_wavelength_m.is_finite() || vacuum_wavelength_m <= 0.0 {
            return Err(LayerPhotoelectricError::Invalid(
                "vacuum photon wavelength is unrepresentable",
            ));
        }
        let law = material.at_wavelength(vacuum_wavelength_m * 1e9)?;
        let absorption =
            illumination.prepare_depth_interval(layer_index, 0.0, layer.thickness_m)?;
        Ok(Self {
            illumination,
            absorption,
            material,
            law,
            layer_index,
            face,
            area_m2,
            vacuum_wavelength_m,
        })
    }

    /// Simultaneous coherent incoming peak tangential E amplitudes [V/m] at
    /// the entrance/exit stack faces, on the original common phasor clock.
    pub fn evaluate(
        &self,
        input: [C64; 2],
    ) -> Result<LayerPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        self.from_absorption(self.absorption.evaluate(input)?)
    }

    /// Statistically independent face sources. Their absorbed powers add;
    /// phase-coherent interference is deliberately absent from this request.
    pub fn evaluate_independent_sources(
        &self,
        input: [C64; 2],
    ) -> Result<LayerPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        self.from_absorption(self.absorption.evaluate_independent_sources(input)?)
    }

    /// Original positive source-covariance factors at this one frequency.
    /// Each factor is internally coherent; powers add between independent
    /// factors. The original interval owner admits all work and cancellation
    /// before publication. Emission is evaluated once on the completed power;
    /// a refusal exposes neither a partial optical nor an emission response.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[[C64; 2]],
        budget: IlluminationBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LayerPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        self.from_absorption(
            self.absorption
                .evaluate_coherence_modes(factors, budget, cancelled)?,
        )
    }

    fn from_absorption(
        &self,
        absorbed_power_w_m2: f64,
    ) -> Result<LayerPhotoelectricResponse<'p, 's, 'm, 'd>, LayerPhotoelectricError> {
        let power = absorbed_power_w_m2 * self.area_m2;
        if !power.is_finite() || (absorbed_power_w_m2 > 0.0 && power == 0.0) {
            return Err(LayerPhotoelectricError::Invalid(
                "physical absorbed film power is unrepresentable",
            ));
        }
        let emission = self.law.evaluate(self.vacuum_wavelength_m, power)?;
        Ok(LayerPhotoelectricResponse {
            illumination: self.illumination,
            material: self.material,
            layer_index: self.layer_index,
            face: self.face,
            area_m2: self.area_m2,
            vacuum_wavelength_m: self.vacuum_wavelength_m,
            absorbed_power_w_m2,
            emission,
        })
    }
}

impl<'p, 's, 'm, 'd> LayerPhotoelectricResponse<'p, 's, 'm, 'd> {
    pub const fn illumination(&self) -> &'p PreparedLayerIllumination<'s, 'm> {
        self.illumination
    }
    pub const fn material(&self) -> &'d PhotoelectricMaterial {
        self.material
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn escape_face(&self) -> EscapeFace {
        self.face
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn vacuum_wavelength_m(&self) -> f64 {
        self.vacuum_wavelength_m
    }
    pub const fn absorbed_power_w_m2(&self) -> f64 {
        self.absorbed_power_w_m2
    }
    /// Residual heat, electron kinetic power and removal work partition the
    /// original absorbed watts. The current is surface charging, not collection.
    pub const fn photoelectric(&self) -> &PhotoelectricResponse {
        &self.emission
    }
}
