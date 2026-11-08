//! Original NASA CHARMS silica optics through native lens and curved-ray owners.
//!
//! The material's absolute phase/group law retains its 30--300 K, 400--2600 nm
//! source domain. Homogeneous lenses reuse the existing Sellmeier owner; smooth
//! temperature fields compose the original analytic n_T and n_TT with their
//! actual spatial derivatives. No fit, table resampling or finite differences
//! are introduced here. Absorption is a separate explicit input: the index
//! measurements do not establish transparency or thermal emissivity.
//!
//! The temperature derivatives are derivatives of the published global index
//! fit, not the source's independently processed thermo-optic measurements.
//! The source owner documents their discrepancy near temperature endpoints.
//! An accepted snapshot is frozen during transport. This adapter supplies no
//! stress optics, thermal expansion or discontinuous-interface reconstruction.

use super::thermo_optic::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureSecondOrderField,
};
use crate::dielectric::{
    BeerLambertAbsorption, DielectricError, DielectricGlass, GlassProvenance, SellmeierIor,
};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use fs_blake3::ContentHash;
use fs_geom::Point3;
use fs_material::silica_refraction::{
    Corning7980Silica, SilicaRefractionError, SilicaRefractionSample, WAVELENGTH_RANGE_NM,
};

/// Resolve the measured absolute-index fit into the original homogeneous glass
/// owner. Supplying `CLEAR` explicitly declares zero absorption; it is not a
/// claim made by NASA's index data. Source bands must have a nonempty overlap.
pub fn corning7980_glass(
    temperature_k: f64,
    absorption: BeerLambertAbsorption,
) -> Result<DielectricGlass, SilicaFieldError> {
    let source = Corning7980Silica::at_temperature(temperature_k)?;
    let (strengths, resonances_um2) = source.sellmeier_coefficients();
    let ior =
        SellmeierIor::try_new_in_band(strengths, resonances_um2, source.wavelength_bounds_nm())?;
    let glass = DielectricGlass::new(ior, absorption, GlassProvenance::CharmsCorning7980V1);
    glass.wavelength_bounds_nm()?;
    Ok(glass)
}

/// One immutable thermal reconstruction and an independent absorption spectrum.
pub struct MaterialSilicaRefractiveField<'a, F: TemperatureGradientField + ?Sized> {
    temperature: &'a F,
    absorption: &'a BeerLambertAbsorption,
    snapshot_identity: ContentHash,
    wavelength_range_nm: [f64; 2],
}

impl<'a, F: TemperatureGradientField + ?Sized> MaterialSilicaRefractiveField<'a, F> {
    pub fn try_new(
        temperature: &'a F,
        absorption: &'a BeerLambertAbsorption,
        snapshot_identity: ContentHash,
    ) -> Result<Self, SilicaFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(SilicaFieldError::MissingSnapshotIdentity);
        }
        let mut band = WAVELENGTH_RANGE_NM;
        if let Some((low, high)) = absorption.wavelength_bounds_nm() {
            band[0] = band[0].max(low);
            band[1] = band[1].min(high);
        }
        if band[0] > band[1] {
            return Err(SilicaFieldError::DisjointSpectralSupport);
        }
        Ok(Self {
            temperature,
            absorption,
            snapshot_identity,
            wavelength_range_nm: band,
        })
    }

    pub const fn temperature_source(&self) -> &'a F {
        self.temperature
    }
    pub const fn absorption_source(&self) -> &'a BeerLambertAbsorption {
        self.absorption
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }
    pub fn law_identity(&self) -> ContentHash {
        Corning7980Silica::source_identity()
    }

    /// Cache the wavelength and bulk absorption without freezing the spatial
    /// temperature. Each stage still reads the original accepted thermal state.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialSilicaRefractiveField<'_, 'a, F>, SilicaFieldError> {
        if !wavelength_nm.is_finite()
            || !(self.wavelength_range_nm[0]..=self.wavelength_range_nm[1]).contains(&wavelength_nm)
        {
            return Err(DielectricError::InvalidWavelength.into());
        }
        Ok(PreparedMaterialSilicaRefractiveField {
            field: self,
            wavelength_nm,
            absorption_per_m: self.absorption.extinction_per_m(wavelength_nm)?,
        })
    }
}

/// Allocation-free wavelength preparation borrowing the original field.
pub struct PreparedMaterialSilicaRefractiveField<
    'binding,
    'source,
    F: TemperatureGradientField + ?Sized,
> {
    field: &'binding MaterialSilicaRefractiveField<'source, F>,
    wavelength_nm: f64,
    absorption_per_m: f64,
}

impl<F: TemperatureGradientField + ?Sized> PreparedMaterialSilicaRefractiveField<'_, '_, F> {
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    pub fn snapshot_identity(&self) -> ContentHash {
        self.field.snapshot_identity
    }

    pub fn optical_sample(&self, p: Point3) -> Result<GradientIndexSample, SilicaFieldError> {
        admit_position(p)?;
        let temperature = self.field.temperature.temperature_gradient(p)?;
        self.convert(temperature).map(|(sample, _)| sample)
    }

    fn convert(
        &self,
        temperature: TemperatureGradientSample,
    ) -> Result<(GradientIndexSample, SilicaRefractionSample), SilicaFieldError> {
        let source = Corning7980Silica::at_temperature(temperature.temperature_k)?;
        let sample = source.evaluate(self.wavelength_nm)?;
        let gradient = temperature
            .gradient_k_per_m
            .scale(sample.derivative_temperature_per_k);
        let optical = GradientIndexSample::try_new(
            sample.phase_index,
            gradient,
            self.absorption_per_m,
            Some(sample.group_index),
        )?;
        Ok((optical, sample))
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> PreparedMaterialSilicaRefractiveField<'_, '_, F> {
    pub fn optical_jet(&self, p: Point3) -> Result<IndexJet, SilicaFieldError> {
        admit_position(p)?;
        let temperature = self.field.temperature.temperature_jet(p)?;
        let (optical, source) = self.convert(temperature.sample)?;
        let first = source.derivative_temperature_per_k;
        let second = source.second_derivative_temperature_per_k2;
        let g = temperature.sample.gradient_k_per_m;
        let diagonal = core::array::from_fn(|i| {
            let v = [g.x, g.y, g.z][i];
            first * temperature.diagonal_k_per_m2[i] + (second * v) * v
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
    for MaterialSilicaRefractiveField<'_, F>
{
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|field| field.optical_sample(p))
            .map_err(SilicaFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialSilicaRefractiveField<'_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|field| field.optical_jet(p))
            .map_err(SilicaFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> GradientIndexField
    for PreparedMaterialSilicaRefractiveField<'_, '_, F>
{
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(p)
            .map_err(SilicaFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.field.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for PreparedMaterialSilicaRefractiveField<'_, '_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(p)
            .map_err(SilicaFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField
    for MaterialSilicaRefractiveField<'_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        admit_position(p)?;
        let temperature = self.temperature.temperature_k(p)?;
        Corning7980Silica::at_temperature(temperature)
            .map_err(|error| SilicaFieldError::Material(error).into_gradient_error())?;
        Ok(temperature)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField
    for PreparedMaterialSilicaRefractiveField<'_, '_, F>
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

/// Refusal of a source, snapshot, or requested physical domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SilicaFieldError {
    MissingSnapshotIdentity,
    DisjointSpectralSupport,
    Material(SilicaRefractionError),
    Optical(DielectricError),
    Field(GradientIndexError),
}

impl SilicaFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Optical(error) => error.into(),
            Self::Material(SilicaRefractionError::TemperatureOutsideRange) => {
                GradientIndexError::OutsideFieldDomain
            }
            Self::Material(SilicaRefractionError::WavelengthOutsideRange) => {
                DielectricError::InvalidWavelength.into()
            }
            Self::Material(SilicaRefractionError::NumericalRange) => {
                GradientIndexError::NumericalRange
            }
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<SilicaRefractionError> for SilicaFieldError {
    fn from(error: SilicaRefractionError) -> Self {
        Self::Material(error)
    }
}
impl From<DielectricError> for SilicaFieldError {
    fn from(error: DielectricError) -> Self {
        Self::Optical(error)
    }
}
impl From<GradientIndexError> for SilicaFieldError {
    fn from(error: GradientIndexError) -> Self {
        Self::Field(error)
    }
}
impl core::fmt::Display for SilicaFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material silica refraction: {self:?}")
    }
}
impl core::error::Error for SilicaFieldError {}

#[cfg(test)]
mod tests;
