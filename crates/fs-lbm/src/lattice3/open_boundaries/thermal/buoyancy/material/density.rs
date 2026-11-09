//! Physical material density on the retained native buoyancy force policy.
//!
//! The original immutable source is queried at each actual Kelvin temperature.
//! Density differences affect buoyancy only; numerical mass and sensible storage
//! retain the reference density and heat capacity. Nonlinear or nonmonotone
//! density is admitted only within the explicit small-departure approximation.
//! This source does not add variable-mass transport or pressure work.

use super::*;
use crate::material::MaterialLbmError;
use crate::thermal::open::buoyancy::{BoussinesqDensityLaw3, DensityDepartureLaw3};

/// Original source/configuration, SI scaling, or native force-policy refusal.
#[derive(Debug)]
pub enum MaterialDensitySourceError3 {
    Material(MaterialTransportError),
    Scaling(MaterialLbmError),
    Law(BoussinesqLawError3),
}
impl From<MaterialTransportError> for MaterialDensitySourceError3 {
    fn from(error: MaterialTransportError) -> Self {
        Self::Material(error)
    }
}
impl From<MaterialLbmError> for MaterialDensitySourceError3 {
    fn from(error: MaterialLbmError) -> Self {
        Self::Scaling(error)
    }
}
impl From<BoussinesqLawError3> for MaterialDensitySourceError3 {
    fn from(error: BoussinesqLawError3) -> Self {
        Self::Law(error)
    }
}
impl core::fmt::Display for MaterialDensitySourceError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native material density source: {self:?}")
    }
}
impl std::error::Error for MaterialDensitySourceError3 {}

/// Original reference source, frozen reference transport and immutable resolver.
///
/// The resolver must reproduce the same physical material/composition/context
/// at the requested temperature. Its original concrete source value supplies
/// density and heat capacity; neither is fitted to an expansion coefficient.
/// Callback effects are outside native state and must be repeatable.
///
/// The retained configuration supplies the SAME Kelvin mapping and neglected
/// rho/Cp bounds as the material transport workspace. Workspace/report byte
/// caps apply when that workspace is used, not to this allocation-free binding.
pub struct MaterialDensitySource3<'a, M: FluidTransport + ?Sized, R> {
    reference_source: &'a M,
    reference: FluidTransportState,
    config: BoussinesqOpenMaterialConfig3,
    resolve: R,
}

impl<M: FluidTransport + ?Sized, R> MaterialDensitySource3<'_, M, R> {
    /// Original reference owner, including its source identity and receipts.
    pub const fn reference_source(&self) -> &M {
        self.reference_source
    }

    /// The exact transport snapshot used by the original SI scaling.
    pub const fn reference_transport(&self) -> &FluidTransportState {
        &self.reference
    }

    pub const fn config(&self) -> BoussinesqOpenMaterialConfig3 {
        self.config
    }

    /// Immutable original resolver; no mutable law/source rebind is exposed.
    pub const fn resolver(&self) -> &R {
        &self.resolve
    }

    /// Query the original source in Kelvin without projecting its physical data.
    ///
    /// The material transport facade can use this same resolver/configuration
    /// for its local viscosity and conductivity queries. Native force sampling
    /// additionally applies the original Boussinesq rho/Cp admission below.
    pub fn resolve_material<F: FluidTransport>(
        &self,
        temperature_k: f64,
    ) -> Result<F, MaterialError>
    where
        R: Fn(f64, &FluidTransportState) -> Result<F, MaterialError>,
    {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(MaterialError::State {
                what: "material density requires a positive finite Kelvin temperature".into(),
            });
        }
        (self.resolve)(temperature_k, &self.reference)
    }
}

impl<M, R, F> DensityDepartureLaw3 for MaterialDensitySource3<'_, M, R>
where
    M: FluidTransport + ?Sized,
    R: Fn(f64, &FluidTransportState) -> Result<F, MaterialError>,
    F: FluidTransport,
{
    type Error = MaterialError;

    fn relative_density_departure(&self, theta: f64) -> Result<f64, Self::Error> {
        let temperature_k = self
            .config
            .temperature(theta)
            .ok_or_else(|| MaterialError::State {
                what: "native temperature cannot be represented in positive finite Kelvin".into(),
            })?;
        let material = BoussinesqMaterial::try_new(
            self.resolve_material(temperature_k)?,
            &self.reference,
            self.config.projection(),
        )?;
        let density = material.resolved_transport().density_kg_m3();
        let reference = self.reference.density_kg_m3();
        Ok((density - reference) / reference)
    }
}

/// Bind original temperature-dependent material density to native 3-D buoyancy.
///
/// Native reference density is one, as prescribed by MaterialScaling. Gravity
/// uses its original SI acceleration conversion; theta uses the explicitly
/// supplied affine Kelvin coordinates, including a negative nonzero scale.
/// Reference rho/Cp remain fixed storage. Source density is used only in
/// (rho(T)-rho_ref)/rho_ref times the native reference density and gravity.
///
/// No source is queried speculatively across an interval: the retained native
/// owner admits actual initialization, source-substep, wall/reservoir, preview
/// and accepted temperatures. Source/domain errors remain MaterialError inside
/// BoussinesqDensityLawError3. A law is immutable after native attachment.
///
/// Use the same reference, coordinate configuration and physical source when
/// composing material viscosity/conductivity or SGS transport. This factory
/// does not infer that two independently supplied source models are identical.
pub fn material_density_law3<'a, M, R, F>(
    reference: &'a MaterialScaling<'_, M>,
    config: BoussinesqOpenMaterialConfig3,
    gravity_m_s2: [f64; 3],
    resolve: R,
) -> Result<BoussinesqDensityLaw3<MaterialDensitySource3<'a, M, R>>, MaterialDensitySourceError3>
where
    M: FluidTransport + ?Sized,
    R: Fn(f64, &FluidTransportState) -> Result<F, MaterialError>,
    F: FluidTransport,
{
    config.projection().validate()?;
    if !config.temperature_at_zero_k.is_finite()
        || !config.temperature_scale_k.is_finite()
        || config.temperature_scale_k == 0.0
    {
        return Err(MaterialTransportError::Invalid(
            "material density requires finite nonzero Kelvin coordinates",
        )
        .into());
    }
    let gravity = [
        reference.acceleration_lattice(gravity_m_s2[0])?,
        reference.acceleration_lattice(gravity_m_s2[1])?,
        reference.acceleration_lattice(gravity_m_s2[2])?,
    ];
    let source = MaterialDensitySource3 {
        reference_source: reference.material(),
        reference: *reference.transport(),
        config,
        resolve,
    };
    Ok(BoussinesqDensityLaw3::try_new(
        1.0,
        gravity,
        config.maximum_relative_density_departure,
        source,
    )?)
}
