//! Original tabulated material phase, group and loss on frozen thermal fields.
//!
//! The shared material owner retains the selected wavelength/index/absorption
//! claims and their exact temperature interpolation. This adapter composes
//! n_T and n_TT with the original spatial temperature derivatives:
//! grad(n) = n_T grad(T), Hess(n) = n_T Hess(T) + n_TT grad(T) grad(T)^T.
//! No fitted thermo-optic correction, spectral resampling or finite difference
//! is added to the source. Preparation admits one smooth interpolation cell
//! over the entire declared temperature interval and caches its wavelength.
//!
//! The reference glass identifies the same native scene medium; it must equal
//! the complete original source spectrum, including absorption and gap policy.
//! Local optics are evaluated from the original thermal path, rather than
//! adding a temperature correction to an already temperature-dependent index.
//! An accepted snapshot remains frozen during transport. Material interfaces,
//! thermal contacts and genuine source cusps require separate native segments.

use super::thermo_optic::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureSecondOrderField,
};
use crate::dielectric::{DielectricError, DielectricGlass, DielectricIor};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use fs_blake3::ContentHash;
use fs_geom::Point3;
use fs_material::state_point::{
    PreparedSpectralDielectricTemperatureCell, SpectralDielectricTemperatureError,
    SpectralDielectricTemperaturePath, SpectralDielectricTemperatureSample,
};

/// One original material source and an immutable accepted temperature snapshot.
pub struct MaterialSpectralRefractiveField<'a, F: TemperatureGradientField + ?Sized> {
    reference: &'a DielectricGlass,
    source: &'a SpectralDielectricTemperaturePath,
    temperature: &'a F,
    snapshot_identity: ContentHash,
    wavelength_range_nm: [f64; 2],
}

impl<'a, F: TemperatureGradientField + ?Sized> MaterialSpectralRefractiveField<'a, F> {
    /// Retain the original table, its interpolation-gap limit, and the exact
    /// common support throughout the declared temperature interval. The caller
    /// associates the source and frozen thermal reconstruction with the same
    /// specimen and held physical coordinates; no material query occurs here.
    pub fn try_new(
        reference: &'a DielectricGlass,
        source: &'a SpectralDielectricTemperaturePath,
        temperature: &'a F,
        snapshot_identity: ContentHash,
    ) -> Result<Self, SpectralRefractiveFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(SpectralRefractiveFieldError::MissingSnapshotIdentity);
        }
        let DielectricIor::Tabulated(table) = reference.ior() else {
            return Err(SpectralRefractiveFieldError::ReferenceSpectrumMismatch);
        };
        let original = super::bulk_optics::spectral_tabulated(
            source.source(),
            table.max_knot_span_nm(),
        ).map_err(|_| SpectralRefractiveFieldError::ReferenceSpectrumMismatch)?;
        if &original != reference {
            return Err(SpectralRefractiveFieldError::ReferenceSpectrumMismatch);
        }
        // A moving source grid must retain the SAME native table policy at
        // every admitted temperature, not only at the reference snapshot.
        let maximum_gap_nm = source.maximum_wavelength_knot_span_m() * 1.0e9;
        if !maximum_gap_nm.is_finite() || maximum_gap_nm > table.max_knot_span_nm() {
            return Err(DielectricError::InvalidTable.into());
        }
        let (minimum, maximum) = source.common_wavelength_bounds_m()
            .ok_or(SpectralRefractiveFieldError::DisjointSpectralSupport)?;
        let wavelength_range_nm = [minimum * 1.0e9, maximum * 1.0e9];
        if wavelength_range_nm.iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(DielectricError::InvalidWavelength.into());
        }
        Ok(Self {
            reference,
            source,
            temperature,
            snapshot_identity,
            wavelength_range_nm,
        })
    }

    pub const fn reference_glass(&self) -> &'a DielectricGlass {
        self.reference
    }
    pub const fn source(&self) -> &'a SpectralDielectricTemperaturePath {
        self.source
    }
    pub const fn temperature_source(&self) -> &'a F {
        self.temperature
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.source.temperature_range_k()
    }

    // Native mesh traversal replaces only the cursor of this SAME snapshot.
    // Preserve the entire already admitted source, support and reference policy.
    pub(crate) fn with_temperature_source<'b, G: TemperatureGradientField + ?Sized>(
        &'b self,
        temperature: &'b G,
    ) -> Result<MaterialSpectralRefractiveField<'b, G>, SpectralRefractiveFieldError> {
        Ok(MaterialSpectralRefractiveField {
            reference: self.reference,
            source: self.source,
            temperature,
            snapshot_identity: self.snapshot_identity,
            wavelength_range_nm: self.wavelength_range_nm,
        })
    }

    /// Admit the original smooth source cell once. Hot spatial sampling borrows
    /// its exact rational interpolation and performs no allocation/card query.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialSpectralRefractiveField<'_, 'a, F>, SpectralRefractiveFieldError> {
        if !wavelength_nm.is_finite()
            || !(self.wavelength_range_nm[0]..=self.wavelength_range_nm[1]).contains(&wavelength_nm)
        {
            return Err(DielectricError::InvalidWavelength.into());
        }
        let (minimum, maximum) = self.source.common_wavelength_bounds_m()
            .ok_or(SpectralRefractiveFieldError::DisjointSpectralSupport)?;
        // Preserve exact original SI endpoints across the nm presentation.
        let wavelength_m = if wavelength_nm == self.wavelength_range_nm[0] {
            minimum
        } else if wavelength_nm == self.wavelength_range_nm[1] {
            maximum
        } else {
            wavelength_nm * 1.0e-9
        };
        Ok(PreparedMaterialSpectralRefractiveField {
            field: self,
            material: self.source.prepare(wavelength_m)?,
            wavelength_nm,
        })
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), GradientIndexError> {
        let [minimum, maximum] = self.temperature_range_k();
        if !temperature_k.is_finite() || !(minimum..=maximum).contains(&temperature_k) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(())
    }
}

/// Borrowed wavelength preparation retaining the original thermal source cell.
pub struct PreparedMaterialSpectralRefractiveField<
    'binding,
    'source,
    F: TemperatureGradientField + ?Sized,
> {
    field: &'binding MaterialSpectralRefractiveField<'source, F>,
    material: PreparedSpectralDielectricTemperatureCell<'binding>,
    wavelength_nm: f64,
}

impl<'binding, 'source, F: TemperatureGradientField + ?Sized>
    PreparedMaterialSpectralRefractiveField<'binding, 'source, F>
{
    pub const fn source_field(&self) -> &'binding MaterialSpectralRefractiveField<'source, F> {
        self.field
    }
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    pub fn snapshot_identity(&self) -> ContentHash {
        self.field.snapshot_identity
    }

    pub fn optical_sample(&self, p: Point3) -> Result<GradientIndexSample, SpectralRefractiveFieldError> {
        admit_position(p)?;
        let temperature = self.field.temperature.temperature_gradient(p)?;
        self.convert(temperature).map(|(optical, _)| optical)
    }

    fn convert(
        &self,
        temperature: TemperatureGradientSample,
    ) -> Result<(GradientIndexSample, SpectralDielectricTemperatureSample), SpectralRefractiveFieldError> {
        let sample = self.material.sample(temperature.temperature_k)?;
        let optical = GradientIndexSample::try_new(
            sample.phase_index,
            temperature.gradient_k_per_m.scale(sample.derivative_temperature_per_k),
            sample.absorption_per_m,
            sample.group_index,
        )?;
        Ok((optical, sample))
    }
}

impl<F: TemperatureSecondOrderField + ?Sized>
    PreparedMaterialSpectralRefractiveField<'_, '_, F>
{
    pub fn optical_jet(&self, p: Point3) -> Result<IndexJet, SpectralRefractiveFieldError> {
        admit_position(p)?;
        let temperature = self.field.temperature.temperature_jet(p)?;
        let (optical, sample) = self.convert(temperature.sample)?;
        let first = sample.derivative_temperature_per_k;
        let second = sample.second_derivative_temperature_per_k2;
        let g = temperature.sample.gradient_k_per_m;
        let diagonal = core::array::from_fn(|i| {
            let value = [g.x, g.y, g.z][i];
            first * temperature.diagonal_k_per_m2[i] + (second * value) * value
        });
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = [(g.x, g.y), (g.x, g.z), (g.y, g.z)][i];
            first * temperature.off_diagonal_k_per_m2[i] + (second * a) * b
        });
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| GradientIndexError::NumericalRange)?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<F: TemperatureGradientField + ?Sized> GradientIndexField
    for MaterialSpectralRefractiveField<'_, F>
{
    fn sample(&self, p: Point3, wavelength_nm: f64) -> Result<GradientIndexSample, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|prepared| prepared.optical_sample(p))
            .map_err(SpectralRefractiveFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialSpectralRefractiveField<'_, F>
{
    fn sample_second_order(&self, p: Point3, wavelength_nm: f64) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|prepared| prepared.optical_jet(p))
            .map_err(SpectralRefractiveFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> GradientIndexField
    for PreparedMaterialSpectralRefractiveField<'_, '_, F>
{
    fn sample(&self, p: Point3, wavelength_nm: f64) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(p).map_err(SpectralRefractiveFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.field.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for PreparedMaterialSpectralRefractiveField<'_, '_, F>
{
    fn sample_second_order(&self, p: Point3, wavelength_nm: f64) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(p).map_err(SpectralRefractiveFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField
    for MaterialSpectralRefractiveField<'_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        admit_position(p)?;
        let temperature = self.temperature.temperature_k(p)?;
        self.admit_temperature(temperature)?;
        Ok(temperature)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField
    for PreparedMaterialSpectralRefractiveField<'_, '_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.field.temperature_k(p)
    }
}

fn admit_position(p: Point3) -> Result<(), GradientIndexError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(GradientIndexError::OutsideFieldDomain)
    }
}

/// Refusal of original source smoothness, support or accepted-state binding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectralRefractiveFieldError {
    MissingSnapshotIdentity,
    ReferenceSpectrumMismatch,
    DisjointSpectralSupport,
    Source(SpectralDielectricTemperatureError),
    Optics(DielectricError),
    Field(GradientIndexError),
}

impl SpectralRefractiveFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Optics(error) => error.into(),
            Self::Source(SpectralDielectricTemperatureError::TemperatureOutsideRange) => {
                GradientIndexError::OutsideFieldDomain
            }
            Self::Source(SpectralDielectricTemperatureError::WavelengthOutsideRange) => {
                DielectricError::InvalidWavelength.into()
            }
            Self::Source(SpectralDielectricTemperatureError::NumericalRange) => {
                GradientIndexError::NumericalRange
            }
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<SpectralDielectricTemperatureError> for SpectralRefractiveFieldError {
    fn from(error: SpectralDielectricTemperatureError) -> Self { Self::Source(error) }
}
impl From<DielectricError> for SpectralRefractiveFieldError {
    fn from(error: DielectricError) -> Self { Self::Optics(error) }
}
impl From<GradientIndexError> for SpectralRefractiveFieldError {
    fn from(error: GradientIndexError) -> Self { Self::Field(error) }
}
impl core::fmt::Display for SpectralRefractiveFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material tabulated thermal refraction: {self:?}")
    }
}
impl core::error::Error for SpectralRefractiveFieldError {}

#[cfg(test)]
pub(super) mod tests;
