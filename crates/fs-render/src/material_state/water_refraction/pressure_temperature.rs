//! Water refraction from the original pressure and temperature reconstruction.
//!
//! Each query evaluates one explicitly selected IF97 liquid/vapor branch at
//! the actual local pressure and temperature. Compressibility, expansion and
//! their mixed response all contribute to the density gradient and Hessian.
//! The existing IAPWS optical owner then supplies phase/group dispersion and
//! the existing ray/bundle owners supply bending and focusing.
//!
//! The source must be one immutable, smooth, single-phase ordinary-water
//! reconstruction. Pressure means absolute thermodynamic pressure, not a
//! pressure perturbation alone. A frozen compressive/acoustic snapshot can
//! drive optics; this does not infer its thermal closure or model moving-fluid
//! Doppler effects. Phase crossings, other IF97 regions and mesh transitions
//! remain with their native owners. Source uncertainties are not reduced.

use super::{
    MaterialWaterRefractiveField, WaterRefractiveFieldError, WaterStateGradientField,
    WaterStateGradientSample, WaterStateJet, WaterStateSecondOrderField,
};
use crate::dielectric::BeerLambertAbsorption;
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{EVALUATOR_VERSION, If97Error, If97Region, If97State, SOURCE_RELEASE};

/// Pressure and temperature gradients from the same accepted physical state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureTemperatureGradientSample {
    pub temperature: TemperatureGradientSample,
    pub pressure_pa: f64,
    pub pressure_gradient_pa_per_m: Vec3,
}

/// Complete spatial Hessians of that same pressure/temperature reconstruction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureTemperatureJet {
    pub temperature: TemperatureJet,
    pub pressure_pa: f64,
    pub pressure_gradient_pa_per_m: Vec3,
    /// `(xx, yy, zz)` pressure Hessian, in Pa/m².
    pub pressure_diagonal_pa_per_m2: [f64; 3],
    /// `(xy, xz, yz)` pressure Hessian, in Pa/m².
    pub pressure_off_diagonal_pa_per_m2: [f64; 3],
}
impl PressureTemperatureJet {
    pub const fn sample(&self) -> PressureTemperatureGradientSample {
        PressureTemperatureGradientSample {
            temperature: self.temperature.sample,
            pressure_pa: self.pressure_pa,
            pressure_gradient_pa_per_m: self.pressure_gradient_pa_per_m,
        }
    }
}

/// One frozen reconstruction, queried jointly to retain state correspondence.
pub trait PressureTemperatureGradientField {
    fn pressure_temperature_gradient(
        &self,
        position_m: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError>;

    /// True only when BOTH fields are constant over ALL finite positions.
    /// Bounded fields keep the default even when their local gradients vanish.
    fn is_spatially_uniform(&self) -> bool {
        false
    }
}

/// Analytic second derivatives of the original joint reconstruction.
/// The value/gradient agrees with `pressure_temperature_gradient`.
pub trait PressureTemperatureSecondOrderField: PressureTemperatureGradientField {
    fn pressure_temperature_jet(
        &self,
        position_m: Point3,
    ) -> Result<PressureTemperatureJet, GradientIndexError>;
}

/// Original single-phase IF97 water state driven by varying pressure and temperature.
pub struct PressureTemperatureWaterField<'a, F: PressureTemperatureGradientField + ?Sized> {
    source: &'a F,
    region: If97Region,
    snapshot_identity: ContentHash,
    identity: ContentHash,
}
impl<'a, F: PressureTemperatureGradientField + ?Sized> PressureTemperatureWaterField<'a, F> {
    /// Bind the source without copying, fitting or smoothing its reconstruction.
    /// The snapshot names BOTH fields and their common reconstruction; it is
    /// not a proof of global source validity. Actual queries retain EOS admission.
    pub fn try_new(
        source: &'a F,
        region: If97Region,
        snapshot_identity: ContentHash,
    ) -> Result<Self, PressureTemperatureWaterError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(PressureTemperatureWaterError::MissingSnapshotIdentity);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.pressure-temperature-if97-water.v1");
        hash.update(snapshot_identity.as_bytes());
        hash.update(SOURCE_RELEASE.as_bytes());
        hash.update(&EVALUATOR_VERSION.to_le_bytes());
        hash.update(&[match region {
            If97Region::Region1 => 1,
            If97Region::Region2 => 2,
        }]);
        Ok(Self { source, region, snapshot_identity, identity: hash.finalize() })
    }

    pub const fn state_source(&self) -> &'a F {
        self.source
    }
    pub const fn region(&self) -> If97Region {
        self.region
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    /// Joint snapshot, original EOS release/evaluator and explicit phase branch.
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Bind the existing prepared IAPWS optical owner and independent absorption law.
    pub fn bind_refraction<'s>(
        &'s self,
        absorption: &'s BeerLambertAbsorption,
    ) -> Result<MaterialWaterRefractiveField<'s, Self>, WaterRefractiveFieldError> {
        MaterialWaterRefractiveField::try_new(self, absorption, self.identity)
    }

    /// Original admitted IF97 state at the reconstruction's actual local `(p,T)`.
    pub fn sample_state(&self, position_m: Point3) -> Result<If97State, PressureTemperatureWaterError> {
        self.state_at_sample(self.source_sample(position_m)?)
    }

    pub fn sample_gradient(
        &self,
        position_m: Point3,
    ) -> Result<WaterStateGradientSample, PressureTemperatureWaterError> {
        let sample = self.source_sample(position_m)?;
        let state = self.state_at_sample(sample)?;
        let rho = state.density().value();
        // Ordinary ray queries require only the already retained first EOS
        // responses, not the extra Gibbs third partials needed by bundles.
        let gradient = density_gradient(
            sample,
            -rho * state.thermal_expansion().value(),
            rho * state.isothermal_compressibility().value(),
        );
        if ![gradient.x, gradient.y, gradient.z].into_iter().all(f64::is_finite) {
            return Err(PressureTemperatureWaterError::Unrepresentable);
        }
        Ok(WaterStateGradientSample {
            temperature: sample.temperature,
            density_kg_m3: rho,
            density_gradient_kg_m4: gradient,
        })
    }

    fn source_sample(
        &self,
        position_m: Point3,
    ) -> Result<PressureTemperatureGradientSample, PressureTemperatureWaterError> {
        admit_position(position_m)?;
        let sample = self.source.pressure_temperature_gradient(position_m)
            .map_err(PressureTemperatureWaterError::Field)?;
        admit_sample(sample)?;
        Ok(sample)
    }

    fn state_at_sample(
        &self,
        sample: PressureTemperatureGradientSample,
    ) -> Result<If97State, PressureTemperatureWaterError> {
        If97State::at_pressure_temperature(
            Pressure::new(sample.pressure_pa),
            Temperature::new(sample.temperature.temperature_k),
            self.region,
        ).map_err(PressureTemperatureWaterError::Thermodynamics)
    }
}

impl<F: PressureTemperatureSecondOrderField + ?Sized> PressureTemperatureWaterField<'_, F> {
    pub fn sample_jet(&self, position_m: Point3) -> Result<WaterStateJet, PressureTemperatureWaterError> {
        admit_position(position_m)?;
        let jet = self.source.pressure_temperature_jet(position_m)
            .map_err(PressureTemperatureWaterError::Field)?;
        let sample = jet.sample();
        admit_sample(sample)?;
        if !jet.temperature.diagonal_k_per_m2.into_iter()
            .chain(jet.temperature.off_diagonal_k_per_m2)
            .chain(jet.pressure_diagonal_pa_per_m2)
            .chain(jet.pressure_off_diagonal_pa_per_m2).all(f64::is_finite)
        {
            return Err(PressureTemperatureWaterError::Field(GradientIndexError::InvalidSample));
        }
        let density = self.state_at_sample(sample)?.density_jet()
            .map_err(PressureTemperatureWaterError::Thermodynamics)?;
        let (dt, dp) = (density.temperature_derivative.value(), density.pressure_derivative.value());
        let dtt = density.temperature_second_derivative.value();
        let dpp = density.pressure_second_derivative.value();
        let dtp = density.pressure_temperature_derivative.value();
        let t = sample.temperature.gradient_k_per_m;
        let p = sample.pressure_gradient_pa_per_m;
        let (gt, gp) = ([t.x, t.y, t.z], [p.x, p.y, p.z]);
        let entry = |a: usize, b: usize, ht: f64, hp: f64| {
            let mixed = (dtp * gt[a]).mul_add(gp[b], (dtp * gp[a]) * gt[b]);
            dt.mul_add(ht, dp.mul_add(hp,
                (dtt * gt[a]).mul_add(gt[b], (dpp * gp[a]).mul_add(gp[b], mixed))))
        };
        let diagonal = core::array::from_fn(|i| {
            entry(i, i, jet.temperature.diagonal_k_per_m2[i], jet.pressure_diagonal_pa_per_m2[i])
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = pairs[i];
            entry(a, b, jet.temperature.off_diagonal_k_per_m2[i], jet.pressure_off_diagonal_pa_per_m2[i])
        });
        let gradient = density_gradient(sample, dt, dp);
        if ![gradient.x, gradient.y, gradient.z].into_iter()
            .chain(diagonal).chain(off_diagonal).all(f64::is_finite)
        {
            return Err(PressureTemperatureWaterError::Unrepresentable);
        }
        Ok(WaterStateJet {
            temperature: jet.temperature,
            density_kg_m3: density.density.value(),
            density_gradient_kg_m4: gradient,
            density_diagonal_kg_m5: diagonal,
            density_off_diagonal_kg_m5: off_diagonal,
        })
    }
}

impl<F: PressureTemperatureGradientField + ?Sized> WaterStateGradientField
    for PressureTemperatureWaterField<'_, F>
{
    fn water_state_gradient(&self, p: Point3) -> Result<WaterStateGradientSample, GradientIndexError> {
        self.sample_gradient(p).map_err(PressureTemperatureWaterError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.source.is_spatially_uniform()
    }
}
impl<F: PressureTemperatureSecondOrderField + ?Sized> WaterStateSecondOrderField
    for PressureTemperatureWaterField<'_, F>
{
    fn water_state_jet(&self, p: Point3) -> Result<WaterStateJet, GradientIndexError> {
        self.sample_jet(p).map_err(PressureTemperatureWaterError::into_gradient_error)
    }
}
impl<F: PressureTemperatureGradientField + ?Sized> TemperatureField
    for PressureTemperatureWaterField<'_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_state(p).map(|s| s.temperature().value())
            .map_err(PressureTemperatureWaterError::into_gradient_error)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PressureTemperatureWaterError {
    MissingSnapshotIdentity,
    Field(GradientIndexError),
    Thermodynamics(If97Error),
    Unrepresentable,
}
impl PressureTemperatureWaterError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Unrepresentable | Self::Thermodynamics(If97Error::Unrepresentable(_)) => {
                GradientIndexError::NumericalRange
            }
            Self::Thermodynamics(_) => GradientIndexError::OutsideFieldDomain,
            Self::MissingSnapshotIdentity => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for PressureTemperatureWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pressure/temperature IF97 water field: {self:?}")
    }
}
impl core::error::Error for PressureTemperatureWaterError {}

fn density_gradient(sample: PressureTemperatureGradientSample, dt: f64, dp: f64) -> Vec3 {
    let t = sample.temperature.gradient_k_per_m;
    let p = sample.pressure_gradient_pa_per_m;
    Vec3::new(dt.mul_add(t.x, dp * p.x), dt.mul_add(t.y, dp * p.y), dt.mul_add(t.z, dp * p.z))
}
fn admit_position(p: Point3) -> Result<(), PressureTemperatureWaterError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(PressureTemperatureWaterError::Field(GradientIndexError::OutsideFieldDomain))
    }
}
fn admit_sample(s: PressureTemperatureGradientSample) -> Result<(), PressureTemperatureWaterError> {
    let t = s.temperature.gradient_k_per_m;
    let p = s.pressure_gradient_pa_per_m;
    if [s.temperature.temperature_k, s.pressure_pa, t.x, t.y, t.z, p.x, p.y, p.z]
        .into_iter().all(f64::is_finite)
    {
        Ok(())
    } else {
        Err(PressureTemperatureWaterError::Field(GradientIndexError::InvalidSample))
    }
}

#[cfg(test)]
mod tests;
