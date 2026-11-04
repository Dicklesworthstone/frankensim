//! Water optics driven by an original thermal field at declared fixed pressure.
//!
//! This adapter evaluates the existing IF97 Region 1 or 2 source at every
//! temperature query. It therefore retains thermal expansion, including its
//! sign change near water's density maximum, instead of assigning an unrelated
//! fixed density to a varying temperature. The same source supplies analytic
//! density gradients and Hessians for the native ray and bundle owners.
//!
//! The caller selects a pure ordinary-water branch, one constant SI pressure,
//! and an accepted immutable temperature reconstruction. Every query retains
//! that explicit branch; crossing saturation or entering another IF97 region
//! refuses rather than interpolating liquid/vapor or inventing a phase change.
//! Spatially varying pressure, advective/Doppler optics, free-surface geometry,
//! and thermal-mesh face transitions are not supplied by this adapter. IF97
//! thermodynamic and IAPWS optical source uncertainties remain independent.

use super::{
    MaterialWaterRefractiveField, WaterRefractiveFieldError, WaterStateGradientField,
    WaterStateGradientSample, WaterStateJet, WaterStateSecondOrderField,
};
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use crate::material_state::thermo_optic::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureJet,
    TemperatureSecondOrderField,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::Point3;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{
    EVALUATOR_VERSION, If97Error, If97Region, If97State, MAXIMUM_PRESSURE_PA, SOURCE_RELEASE,
};

/// One original thermal reconstruction with an explicit constant-pressure EOS.
pub struct IsobaricWaterField<'a, F: TemperatureGradientField + ?Sized> {
    temperature: &'a F,
    pressure: Pressure,
    region: If97Region,
    snapshot_identity: ContentHash,
    identity: ContentHash,
}

impl<'a, F: TemperatureGradientField + ?Sized> IsobaricWaterField<'a, F> {
    /// Bind the original field without copying, smoothing or fitting it.
    ///
    /// The snapshot identity names the supplied temperature field and its
    /// reconstruction. It is a declared binding, not a global source-domain or
    /// smoothness proof. The native EOS and optical law admit actual query states.
    pub fn try_new(
        temperature: &'a F,
        pressure: Pressure,
        region: If97Region,
        snapshot_identity: ContentHash,
    ) -> Result<Self, IsobaricWaterError> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 || p > MAXIMUM_PRESSURE_PA {
            return Err(IsobaricWaterError::InvalidPressure);
        }
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(IsobaricWaterError::MissingSnapshotIdentity);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.isobaric-if97-water.v1");
        hash.update(snapshot_identity.as_bytes());
        hash.update(SOURCE_RELEASE.as_bytes());
        hash.update(&EVALUATOR_VERSION.to_le_bytes());
        hash.update(&p.to_le_bytes());
        hash.update(&[match region {
            If97Region::Region1 => 1,
            If97Region::Region2 => 2,
        }]);
        Ok(Self {
            temperature,
            pressure,
            region,
            snapshot_identity,
            identity: hash.finalize(),
        })
    }

    pub const fn temperature_source(&self) -> &'a F {
        self.temperature
    }
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    pub const fn region(&self) -> If97Region {
        self.region
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    /// Original snapshot, EOS release/evaluator, pressure and explicit branch.
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Bind the existing IAPWS optical owner and independent absorption law.
    /// Wavelength preparation, phase/group transport and bundle propagation
    /// remain with that owner; no alternate optical engine is introduced.
    pub fn bind_refraction<'s>(
        &'s self,
        absorption: &'s BeerLambertAbsorption,
    ) -> Result<MaterialWaterRefractiveField<'s, Self>, WaterRefractiveFieldError> {
        MaterialWaterRefractiveField::try_new(self, absorption, self.identity)
    }

    /// Original admitted IF97 state at the supplied field's actual temperature.
    pub fn sample_state(&self, position_m: Point3) -> Result<If97State, IsobaricWaterError> {
        admit_position(position_m)?;
        let t = self
            .temperature
            .temperature_k(position_m)
            .map_err(IsobaricWaterError::Field)?;
        self.state_at_temperature(t)
    }

    fn state_at_temperature(&self, t: f64) -> Result<If97State, IsobaricWaterError> {
        If97State::at_pressure_temperature(self.pressure, Temperature::new(t), self.region)
            .map_err(IsobaricWaterError::Thermodynamics)
    }

    pub fn sample_gradient(
        &self,
        position_m: Point3,
    ) -> Result<WaterStateGradientSample, IsobaricWaterError> {
        admit_position(position_m)?;
        let temperature = self
            .temperature
            .temperature_gradient(position_m)
            .map_err(IsobaricWaterError::Field)?;
        admit_temperature_gradient(temperature)?;
        let state = self.state_at_temperature(temperature.temperature_k)?;
        let density = state.density().value();
        let slope = -density * state.thermal_expansion().value();
        let gradient = temperature.gradient_k_per_m.scale(slope);
        if ![gradient.x, gradient.y, gradient.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(IsobaricWaterError::Unrepresentable);
        }
        Ok(WaterStateGradientSample {
            temperature,
            density_kg_m3: density,
            density_gradient_kg_m4: gradient,
        })
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> IsobaricWaterField<'_, F> {
    pub fn sample_jet(&self, position_m: Point3) -> Result<WaterStateJet, IsobaricWaterError> {
        admit_position(position_m)?;
        let temperature: TemperatureJet = self
            .temperature
            .temperature_jet(position_m)
            .map_err(IsobaricWaterError::Field)?;
        admit_temperature_gradient(temperature.sample)?;
        if !temperature
            .diagonal_k_per_m2
            .into_iter()
            .chain(temperature.off_diagonal_k_per_m2)
            .all(f64::is_finite)
        {
            return Err(IsobaricWaterError::Field(GradientIndexError::InvalidSample));
        }
        let state = self.state_at_temperature(temperature.sample.temperature_k)?;
        let density = state
            .isobaric_density_jet()
            .map_err(IsobaricWaterError::Thermodynamics)?;
        let first = density.temperature_derivative.value();
        let second = density.temperature_second_derivative.value();
        let v = temperature.sample.gradient_k_per_m;
        let g = [v.x, v.y, v.z];
        let diagonal = core::array::from_fn(|i| {
            first.mul_add(temperature.diagonal_k_per_m2[i], (second * g[i]) * g[i])
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = pairs[i];
            first.mul_add(temperature.off_diagonal_k_per_m2[i], (second * g[a]) * g[b])
        });
        let gradient = v.scale(first);
        if ![gradient.x, gradient.y, gradient.z]
            .into_iter()
            .chain(diagonal)
            .chain(off_diagonal)
            .all(f64::is_finite)
        {
            return Err(IsobaricWaterError::Unrepresentable);
        }
        Ok(WaterStateJet {
            temperature,
            density_kg_m3: density.density.value(),
            density_gradient_kg_m4: gradient,
            density_diagonal_kg_m5: diagonal,
            density_off_diagonal_kg_m5: off_diagonal,
        })
    }
}

impl<F: TemperatureGradientField + ?Sized> WaterStateGradientField for IsobaricWaterField<'_, F> {
    fn water_state_gradient(
        &self,
        p: Point3,
    ) -> Result<WaterStateGradientSample, GradientIndexError> {
        self.sample_gradient(p)
            .map_err(IsobaricWaterError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.temperature.is_spatially_uniform()
    }
}
impl<F: TemperatureSecondOrderField + ?Sized> WaterStateSecondOrderField
    for IsobaricWaterField<'_, F>
{
    fn water_state_jet(&self, p: Point3) -> Result<WaterStateJet, GradientIndexError> {
        self.sample_jet(p)
            .map_err(IsobaricWaterError::into_gradient_error)
    }
}
impl<F: TemperatureGradientField + ?Sized> TemperatureField for IsobaricWaterField<'_, F> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_state(p)
            .map(|s| s.temperature().value())
            .map_err(IsobaricWaterError::into_gradient_error)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricWaterError {
    InvalidPressure,
    MissingSnapshotIdentity,
    Field(GradientIndexError),
    Thermodynamics(If97Error),
    Unrepresentable,
}
impl IsobaricWaterError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Unrepresentable | Self::Thermodynamics(If97Error::Unrepresentable(_)) => {
                GradientIndexError::NumericalRange
            }
            Self::Thermodynamics(_) => GradientIndexError::OutsideFieldDomain,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for IsobaricWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isobaric IF97 water field: {self:?}")
    }
}
impl core::error::Error for IsobaricWaterError {}

fn admit_position(p: Point3) -> Result<(), IsobaricWaterError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(IsobaricWaterError::Field(
            GradientIndexError::OutsideFieldDomain,
        ))
    }
}
fn admit_temperature_gradient(s: TemperatureGradientSample) -> Result<(), IsobaricWaterError> {
    let g = s.gradient_k_per_m;
    if [s.temperature_k, g.x, g.y, g.z]
        .into_iter()
        .all(f64::is_finite)
    {
        Ok(())
    } else {
        Err(IsobaricWaterError::Field(GradientIndexError::InvalidSample))
    }
}

#[cfg(test)]
mod tests;
