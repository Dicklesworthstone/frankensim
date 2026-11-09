//! Temperature-driven refraction on the native curved-ray and bundle owners.
//!
//! The material owner evaluates its declared thermo-optic law. This adapter
//! composes its exact temperature derivatives with an original smooth thermal
//! snapshot: grad(n) = n_T grad(T), and Hess(n) = n_T Hess(T)
//! + n_TT grad(T) grad(T)^T. The original reference dispersion supplies both
//! phase and group indices. Preparing a wavelength reuses that evaluation at
//! every subsequent optical stage, with no heap allocation or finite differences.
//!
//! Reference glass must describe absolute vacuum-wavelength indices at the
//! law's reference temperature. The snapshot supplies the actual material
//! temperature, not the reservoir temperature. Absorption retains the declared
//! reference-glass law unless an original material absorption path is supplied.
//! That path evaluates the same local temperature without material queries in
//! the ray stages. No thermal expansion, stress birefringence, moving-medium
//! correction or missing material data is inferred.
//! An optical path sees one frozen, isotropic state. At material interfaces or
//! thermal contact jumps, callers must use separate regions and native boundary
//! events rather than integrate a smooth segment across the discontinuity.

use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};
use fs_material::state_point::{
    AbsorptionTemperatureSampleError, SpectralDielectricAbsorptionTemperaturePath,
};
use fs_material::thermo_optic::{
    PreparedThermoOptic, SchottThermoOpticLaw, ThermoOpticError, ThermoOpticSample,
};

use crate::dielectric::{
    BeerLambertAbsorption, DielectricError, DielectricGlass, GlassProvenance, SchottThermoOpticIor,
};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexSample, HomogeneousIndexField,
};

/// Accepted material endpoints retaining one original reference optical source.
pub mod binding;
/// Original solved planar conduction profiles as thermal-lensing snapshots.
pub mod conduction;
/// Original temperature-dependent conductivity profiles, without source-knot smoothing.
pub mod nonlinear_conduction;
/// Original transient P1 temperature snapshots on their tetrahedral elements.
pub mod transient;

/// Bind a uniform accepted temperature to the same material law used by smooth
/// thermal fields. The result enters original lenses, surface/nested-medium
/// transport, films and bulk phase/group propagation without a spectral refit.
///
/// Reference absorption is retained verbatim and intersects the thermal law's
/// source band. It is not a temperature-dependent absorption measurement. The
/// result is a caller-authored composition; its native identity binds the exact
/// reference dispersion, material-law source and frozen temperature. Changes
/// of geometry, surface state or material-card claim selection remain explicit.
pub fn schott_glass_at_temperature(
    reference: &DielectricGlass,
    law: &SchottThermoOpticLaw,
    temperature_k: f64,
) -> Result<DielectricGlass, ThermoOpticFieldError> {
    let ior = SchottThermoOpticIor::try_new(*reference.ior(), *law, temperature_k)?;
    let glass = DielectricGlass::new(ior, *reference.absorption_ref(), GlassProvenance::Custom);
    glass.wavelength_bounds_nm()?;
    Ok(glass)
}

/// Temperature and the analytic gradient of that same immutable reconstruction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureGradientSample {
    pub temperature_k: f64,
    pub gradient_k_per_m: Vec3,
}

/// One temperature/gradient query plus its symmetric analytic spatial Hessian.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureJet {
    pub sample: TemperatureGradientSample,
    /// `(xx, yy, zz)`, in kelvin per square metre.
    pub diagonal_k_per_m2: [f64; 3],
    /// `(xy, xz, yz)`, in kelvin per square metre.
    pub off_diagonal_k_per_m2: [f64; 3],
}

/// An immutable thermal reconstruction with an original analytic gradient.
/// Values must agree with [`TemperatureField`] at the same point. Providers
/// must refuse outside their physical region and resolve their spatial features
/// within the optical step envelope; cell averages alone do not supply this law.
pub trait TemperatureGradientField: TemperatureField {
    fn temperature_gradient(
        &self,
        position_m: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError>;

    /// True only for a temperature constant over all finite positions. A bounded
    /// region must keep the default so a native uniform step cannot skip it.
    fn is_spatially_uniform(&self) -> bool {
        false
    }
}

/// Explicit second derivative capability; a locally zero gradient proves no
/// Hessian. The returned sample must agree with the first-order query.
pub trait TemperatureSecondOrderField: TemperatureGradientField {
    fn temperature_jet(&self, position_m: Point3) -> Result<TemperatureJet, GradientIndexError>;
}

/// Original reference optics and thermo-optic law bound to one thermal snapshot.
pub struct MaterialThermoOpticField<'a, F: TemperatureGradientField + ?Sized> {
    reference: &'a DielectricGlass,
    law: &'a SchottThermoOpticLaw,
    temperature: &'a F,
    snapshot_identity: ContentHash,
    wavelength_range_nm: [f64; 2],
    temperature_range_k: [f64; 2],
    absorption: Option<&'a SpectralDielectricAbsorptionTemperaturePath>,
}

impl<'a, F: TemperatureGradientField + ?Sized> MaterialThermoOpticField<'a, F> {
    /// Retain the exact original sources and their common inclusive wavelength
    /// band. A shared endpoint is usable for a monochromatic path. The supplied
    /// identity associates the immutable thermal snapshot; it authenticates no
    /// measurement. A changed accepted state requires a fresh binding.
    pub fn try_new(
        reference: &'a DielectricGlass,
        law: &'a SchottThermoOpticLaw,
        temperature: &'a F,
        snapshot_identity: ContentHash,
    ) -> Result<Self, ThermoOpticFieldError> {
        Self::from_sources(reference, law, temperature, snapshot_identity, None)
    }

    /// Combine original SCHOTT refraction with the original material's prepared
    /// temperature-dependent intensity absorption. Both read the same pointwise
    /// temperature, including inside native element-transition collars.
    ///
    /// The caller associates the phase reference, thermal law and absorption
    /// source with the same physical specimen and held environment. The path
    /// retains its original selected claims and all supplied requirements; its
    /// construction checks their continuous temperature support. A changed
    /// snapshot requires a fresh binding. No spatial absorption fit is made.
    ///
    /// Only wavelengths supported throughout the declared absorption path are
    /// available. The discarded reference loss does not restrict this band.
    /// Thermal support is the intersection of the phase and absorption laws.
    /// Prepared ray/bundle samples allocate nothing and perform no card queries.
    pub fn try_with_absorption_source(
        reference: &'a DielectricGlass,
        law: &'a SchottThermoOpticLaw,
        temperature: &'a F,
        absorption: &'a SpectralDielectricAbsorptionTemperaturePath,
        snapshot_identity: ContentHash,
    ) -> Result<Self, ThermoOpticFieldError> {
        Self::from_sources(
            reference,
            law,
            temperature,
            snapshot_identity,
            Some(absorption),
        )
    }

    fn from_sources(
        reference: &'a DielectricGlass,
        law: &'a SchottThermoOpticLaw,
        temperature: &'a F,
        snapshot_identity: ContentHash,
        absorption: Option<&'a SpectralDielectricAbsorptionTemperaturePath>,
    ) -> Result<Self, ThermoOpticFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(ThermoOpticFieldError::MissingSnapshotIdentity);
        }
        let (reference_min, reference_max) = if absorption.is_some() {
            reference.ior().wavelength_bounds_nm()
        } else {
            reference.wavelength_bounds_nm()?
        };
        let [law_min, law_max] = law.wavelength_range_nm();
        let mut wavelength_range_nm = [reference_min.max(law_min), reference_max.min(law_max)];
        let mut temperature_range_k = law.temperature_range_k();
        if let Some(source) = absorption {
            let (lo, hi) = source
                .common_wavelength_bounds_m()
                .ok_or(ThermoOpticFieldError::DisjointSpectralSupport)?;
            wavelength_range_nm[0] = wavelength_range_nm[0].max(lo * 1.0e9);
            wavelength_range_nm[1] = wavelength_range_nm[1].min(hi * 1.0e9);
            let [minimum, maximum] = source.temperature_range_k();
            temperature_range_k[0] = temperature_range_k[0].max(minimum);
            temperature_range_k[1] = temperature_range_k[1].min(maximum);
        }
        if wavelength_range_nm[0] > wavelength_range_nm[1] {
            return Err(ThermoOpticFieldError::DisjointSpectralSupport);
        }
        if temperature_range_k[0] > temperature_range_k[1] {
            return Err(ThermoOpticFieldError::DisjointTemperatureSupport);
        }
        Ok(Self {
            reference,
            law,
            temperature,
            snapshot_identity,
            wavelength_range_nm,
            temperature_range_k,
            absorption,
        })
    }

    pub const fn reference_glass(&self) -> &DielectricGlass {
        self.reference
    }
    pub const fn law(&self) -> &SchottThermoOpticLaw {
        self.law
    }
    pub const fn temperature_source(&self) -> &F {
        self.temperature
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    pub const fn absorption_source(&self) -> Option<&SpectralDielectricAbsorptionTemperaturePath> {
        self.absorption
    }

    // Only replace the geometry cursor of the SAME frozen snapshot. Internal
    // face traversal must retain the exact absorption policy and prepared source.
    fn with_temperature_source<'b, G: TemperatureGradientField + ?Sized>(
        &'b self,
        temperature: &'b G,
    ) -> MaterialThermoOpticField<'b, G> {
        MaterialThermoOpticField {
            reference: self.reference,
            law: self.law,
            temperature,
            snapshot_identity: self.snapshot_identity,
            wavelength_range_nm: self.wavelength_range_nm,
            temperature_range_k: self.temperature_range_k,
            absorption: self.absorption,
        }
    }

    /// Evaluate reference phase/group/absorption and thermo-optic spectral
    /// coefficients once. The prepared field still queries original temperature
    /// and derivatives at each stage; it never freezes a spatially varying T.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialThermoOpticField<'_, 'a, F>, ThermoOpticFieldError> {
        if !wavelength_nm.is_finite()
            || !(self.wavelength_range_nm[0]..=self.wavelength_range_nm[1]).contains(&wavelength_nm)
        {
            return Err(DielectricError::InvalidWavelength.into());
        }
        let (reference, absorption) = if let Some(source) = self.absorption {
            // Reuse native phase/group admission without evaluating discarded
            // reference alpha or accidentally retaining its narrower band.
            let phase = DielectricGlass::with_ior(
                *self.reference.ior(),
                BeerLambertAbsorption::CLEAR,
                self.reference.provenance(),
            );
            let sample = HomogeneousIndexField::new(&phase)
                .sample(Point3::new(0.0, 0.0, 0.0), wavelength_nm)?;
            let (lo, hi) = source
                .common_wavelength_bounds_m()
                .ok_or(ThermoOpticFieldError::DisjointSpectralSupport)?;
            // Preserve exact original SI endpoints across the nm presentation.
            let wavelength_m = if wavelength_nm == lo * 1.0e9 {
                lo
            } else if wavelength_nm == hi * 1.0e9 {
                hi
            } else {
                wavelength_nm * 1.0e-9
            };
            (
                sample,
                PreparedAbsorption::Current {
                    source,
                    wavelength_m,
                },
            )
        } else {
            let sample = HomogeneousIndexField::new(self.reference)
                .sample(Point3::new(0.0, 0.0, 0.0), wavelength_nm)?;
            (
                sample,
                PreparedAbsorption::Reference(sample.absorption_per_m()),
            )
        };
        let prepared = self.law.prepare(
            wavelength_nm,
            reference.phase_index(),
            reference.group_index(),
        )?;
        Ok(PreparedMaterialThermoOpticField {
            field: self,
            prepared,
            absorption,
        })
    }

    pub fn optical_sample(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, ThermoOpticFieldError> {
        self.prepare(wavelength_nm)?.optical_sample(position_m)
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), GradientIndexError> {
        let [minimum, maximum] = self.temperature_range_k;
        if !temperature_k.is_finite() || !(minimum..=maximum).contains(&temperature_k) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum PreparedAbsorption<'a> {
    Reference(f64),
    Current {
        source: &'a SpectralDielectricAbsorptionTemperaturePath,
        wavelength_m: f64,
    },
}

/// Stack-owned spectral preparation borrowing the original thermal snapshot.
pub struct PreparedMaterialThermoOpticField<'binding, 'source, F: TemperatureGradientField + ?Sized>
{
    field: &'binding MaterialThermoOpticField<'source, F>,
    prepared: PreparedThermoOptic,
    absorption: PreparedAbsorption<'source>,
}

impl<F: TemperatureGradientField + ?Sized> PreparedMaterialThermoOpticField<'_, '_, F> {
    pub fn wavelength_nm(&self) -> f64 {
        self.prepared.wavelength_nm()
    }

    pub fn snapshot_identity(&self) -> ContentHash {
        self.field.snapshot_identity
    }

    pub fn optical_sample(
        &self,
        position_m: Point3,
    ) -> Result<GradientIndexSample, ThermoOpticFieldError> {
        admit_position(position_m)?;
        let sample = self.field.temperature.temperature_gradient(position_m)?;
        self.convert_sample(sample).map(|(optical, _)| optical)
    }

    fn convert_sample(
        &self,
        temperature: TemperatureGradientSample,
    ) -> Result<(GradientIndexSample, ThermoOpticSample), ThermoOpticFieldError> {
        self.field.admit_temperature(temperature.temperature_k)?;
        let gradient = temperature.gradient_k_per_m;
        if ![gradient.x, gradient.y, gradient.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(GradientIndexError::InvalidSample.into());
        }
        let sample = self.prepared.sample(temperature.temperature_k)?;
        let absorption_per_m = match self.absorption {
            PreparedAbsorption::Reference(value) => value,
            PreparedAbsorption::Current {
                source,
                wavelength_m,
            } => source
                .absorption_per_m_at(temperature.temperature_k, wavelength_m)
                .map_err(ThermoOpticFieldError::Absorption)?,
        };
        let optical = GradientIndexSample::try_new(
            sample.phase_index,
            gradient.scale(sample.derivative_temperature_per_k),
            absorption_per_m,
            sample.group_index,
        )?;
        Ok((optical, sample))
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> PreparedMaterialThermoOpticField<'_, '_, F> {
    pub fn optical_jet(&self, position_m: Point3) -> Result<IndexJet, ThermoOpticFieldError> {
        admit_position(position_m)?;
        let temperature = self.field.temperature.temperature_jet(position_m)?;
        if !temperature
            .diagonal_k_per_m2
            .into_iter()
            .chain(temperature.off_diagonal_k_per_m2)
            .all(f64::is_finite)
        {
            return Err(GradientIndexError::InvalidSample.into());
        }
        let (optical, material) = self.convert_sample(temperature.sample)?;
        let g = temperature.sample.gradient_k_per_m;
        let first = material.derivative_temperature_per_k;
        let second = material.second_derivative_temperature_per_k2;
        let mut diagonal = temperature.diagonal_k_per_m2.map(|value| first * value);
        let mut off_diagonal = temperature.off_diagonal_k_per_m2.map(|value| first * value);
        for (entry, gradient) in diagonal.iter_mut().zip([g.x, g.y, g.z]) {
            // Multiply the temperature coefficient first: an intermediate
            // squared gradient can overflow although the physical jet is finite.
            *entry += (second * gradient) * gradient;
        }
        for (entry, (a, b)) in off_diagonal
            .iter_mut()
            .zip([(g.x, g.y), (g.x, g.z), (g.y, g.z)])
        {
            *entry += (second * a) * b;
        }
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| GradientIndexError::NumericalRange)?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<F: TemperatureGradientField + ?Sized> GradientIndexField for MaterialThermoOpticField<'_, F> {
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.optical_sample(p, wavelength_nm)
            .map_err(ThermoOpticFieldError::into_gradient_error)
    }

    fn is_spatially_uniform(&self) -> bool {
        self.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialThermoOpticField<'_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|prepared| prepared.optical_jet(p))
            .map_err(ThermoOpticFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> GradientIndexField
    for PreparedMaterialThermoOpticField<'_, '_, F>
{
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.prepared.wavelength_nm() {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(p)
            .map_err(ThermoOpticFieldError::into_gradient_error)
    }

    fn is_spatially_uniform(&self) -> bool {
        self.field.temperature.is_spatially_uniform()
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for PreparedMaterialThermoOpticField<'_, '_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.prepared.wavelength_nm() {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(p)
            .map_err(ThermoOpticFieldError::into_gradient_error)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField for MaterialThermoOpticField<'_, F> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        admit_position(p)?;
        let temperature = self.temperature.temperature_k(p)?;
        self.admit_temperature(temperature)?;
        Ok(temperature)
    }
}

impl<F: TemperatureGradientField + ?Sized> TemperatureField
    for PreparedMaterialThermoOpticField<'_, '_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.field.temperature_k(p)
    }
}

/// Source or physical-domain refusal before publishing an optical sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermoOpticFieldError {
    MissingSnapshotIdentity,
    DisjointSpectralSupport,
    DisjointTemperatureSupport,
    Absorption(AbsorptionTemperatureSampleError),
    Material(ThermoOpticError),
    Optics(DielectricError),
    Field(GradientIndexError),
}

impl ThermoOpticFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Optics(error) => error.into(),
            Self::Material(ThermoOpticError::NumericalRange) => GradientIndexError::NumericalRange,
            Self::Material(ThermoOpticError::WavelengthOutsideRange) => {
                DielectricError::InvalidWavelength.into()
            }
            Self::Material(ThermoOpticError::TemperatureOutsideRange) => {
                GradientIndexError::OutsideFieldDomain
            }
            Self::Absorption(AbsorptionTemperatureSampleError::TemperatureOutsideRange) => {
                GradientIndexError::OutsideFieldDomain
            }
            Self::Absorption(AbsorptionTemperatureSampleError::WavelengthOutsideRange) => {
                DielectricError::InvalidWavelength.into()
            }
            Self::Absorption(AbsorptionTemperatureSampleError::NumericalRange) => {
                GradientIndexError::NumericalRange
            }
            _ => GradientIndexError::InvalidSample,
        }
    }
}

impl From<ThermoOpticError> for ThermoOpticFieldError {
    fn from(error: ThermoOpticError) -> Self {
        Self::Material(error)
    }
}
impl From<DielectricError> for ThermoOpticFieldError {
    fn from(error: DielectricError) -> Self {
        Self::Optics(error)
    }
}
impl From<GradientIndexError> for ThermoOpticFieldError {
    fn from(error: GradientIndexError) -> Self {
        Self::Field(error)
    }
}
impl core::fmt::Display for ThermoOpticFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material thermo-optic field: {self:?}")
    }
}
impl core::error::Error for ThermoOpticFieldError {}

fn admit_position(p: Point3) -> Result<(), GradientIndexError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(GradientIndexError::OutsideFieldDomain)
    }
}

#[cfg(test)]
mod tests;
