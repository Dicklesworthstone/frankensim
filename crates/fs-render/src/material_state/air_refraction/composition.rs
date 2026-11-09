//! Atmospheric composition gradients and curvature on original Ciddor optics.
//!
//! Borrow the accepted p/T and composition reconstructions at the same world
//! position. Water is a mole fraction; CO2 is µmol/mol of the DRY component,
//! not relative humidity or total-mixture ppm. Their joint immutable snapshot
//! supplies every value, gradient and Hessian; no EOS, mixing or condensation
//! closure is inferred. The explicit nonabsorbing approximation matches the
//! fixed-composition atmospheric adapter.

use super::{AirFieldError, admit_position, field_error};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::material_state::water_refraction::pressure_temperature::{
    PressureTemperatureGradientField, PressureTemperatureGradientSample,
    PressureTemperatureSecondOrderField,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_material::air_refraction::{
    AirRefractiveStateJet, CiddorAir, PreparedCiddorAirSpectrum, WAVELENGTH_RANGE_NM,
    ciddor_model_identity,
};

/// Original admitted composition and gradients of its two independent coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AirCompositionGradientSample {
    pub composition: CiddorAir,
    pub water_gradient_per_m: Vec3,
    pub co2_gradient_micromol_per_mol_m: Vec3,
}

/// Analytic spatial Hessians; diagonals are `(xx,yy,zz)`, mixed entries `(xy,xz,yz)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AirCompositionJet {
    pub sample: AirCompositionGradientSample,
    pub water_diagonal_per_m2: [f64; 3],
    pub water_off_diagonal_per_m2: [f64; 3],
    pub co2_diagonal_micromol_per_mol_m2: [f64; 3],
    pub co2_off_diagonal_micromol_per_mol_m2: [f64; 3],
}

/// Smooth original composition reconstruction in world metres, at the p/T snapshot.
pub trait AirCompositionGradientField {
    fn air_composition_gradient(
        &self,
        position_m: Point3,
    ) -> Result<AirCompositionGradientSample, GradientIndexError>;

    /// True only when both composition coordinates are constant at ALL finite positions.
    /// A bounded domain or locally zero gradient must keep the default.
    fn is_spatially_uniform(&self) -> bool {
        false
    }
}

/// Original Hessians of the same composition returned by the first-order query.
pub trait AirCompositionSecondOrderField: AirCompositionGradientField {
    fn air_composition_jet(
        &self,
        position_m: Point3,
    ) -> Result<AirCompositionJet, GradientIndexError>;
}

/// Nonabsorbing atmospheric field retaining original p/T and composition providers.
pub struct MaterialAirCompositionRefractiveField<
    'a,
    P: PressureTemperatureGradientField + ?Sized,
    C: AirCompositionGradientField + ?Sized,
> {
    state: &'a P,
    composition: &'a C,
    snapshot_identity: ContentHash,
    content_identity: ContentHash,
}

impl<'a, P: PressureTemperatureGradientField + ?Sized, C: AirCompositionGradientField + ?Sized>
    MaterialAirCompositionRefractiveField<'a, P, C>
{
    /// The nonzero snapshot binds BOTH reconstructions and their common physical time.
    /// Source phase validity and spatial domain admission remain with their providers.
    pub fn try_nonabsorbing(
        state: &'a P,
        composition: &'a C,
        snapshot_identity: ContentHash,
    ) -> Result<Self, AirFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(AirFieldError::MissingSnapshotIdentity);
        }
        let mut hash = DomainHasher::new("fs-render:nonabsorbing-ciddor-composition-field:v1");
        hash.update(ciddor_model_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        Ok(Self {
            state,
            composition,
            snapshot_identity,
            content_identity: hash.finalize(),
        })
    }
    pub const fn state_source(&self) -> &'a P {
        self.state
    }
    pub const fn composition_source(&self) -> &'a C {
        self.composition
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.content_identity
    }
    pub fn model_identity(&self) -> ContentHash {
        ciddor_model_identity()
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        WAVELENGTH_RANGE_NM
    }

    /// Cache composition-independent spectral coefficients outside integration stages.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialAirCompositionRefractiveField<'_, 'a, P, C>, AirFieldError> {
        Ok(PreparedMaterialAirCompositionRefractiveField {
            source: self,
            law: PreparedCiddorAirSpectrum::try_new(wavelength_nm)?,
        })
    }
}

/// Stack-owned wavelength preparation; every stage still queries the original state.
pub struct PreparedMaterialAirCompositionRefractiveField<
    'f,
    'a,
    P: PressureTemperatureGradientField + ?Sized,
    C: AirCompositionGradientField + ?Sized,
> {
    source: &'f MaterialAirCompositionRefractiveField<'a, P, C>,
    law: PreparedCiddorAirSpectrum,
}

impl<'f, 'a, P: PressureTemperatureGradientField + ?Sized, C: AirCompositionGradientField + ?Sized>
    PreparedMaterialAirCompositionRefractiveField<'f, 'a, P, C>
{
    pub const fn source_field(&self) -> &'f MaterialAirCompositionRefractiveField<'a, P, C> {
        self.source
    }
    pub fn wavelength_nm(&self) -> f64 {
        self.law.wavelength_nm()
    }

    pub fn optical_sample(&self, point: Point3) -> Result<GradientIndexSample, AirFieldError> {
        admit_position(point)?;
        let state = self.source.state.pressure_temperature_gradient(point)?;
        let composition = self.source.composition.air_composition_gradient(point)?;
        self.evaluate(state, composition)
            .map(|(sample, _, _)| sample)
    }

    fn evaluate(
        &self,
        state: PressureTemperatureGradientSample,
        composition: AirCompositionGradientSample,
    ) -> Result<(GradientIndexSample, AirRefractiveStateJet, [[f64; 3]; 4]), AirFieldError> {
        let law = self.law.sample_jet(
            state.temperature.temperature_k,
            state.pressure_pa,
            composition.composition,
        )?;
        let gradients = [
            state.temperature.gradient_k_per_m,
            state.pressure_gradient_pa_per_m,
            composition.water_gradient_per_m,
            composition.co2_gradient_micromol_per_mol_m,
        ]
        .map(|g| [g.x, g.y, g.z]);
        let gradient: [f64; 3] = core::array::from_fn(|i| {
            (0..4)
                .map(|a| law.phase_gradient[a] * gradients[a][i])
                .sum()
        });
        let optical = GradientIndexSample::try_new(
            law.sample.phase_index,
            Vec3::new(gradient[0], gradient[1], gradient[2]),
            0.0,
            Some(law.sample.group_index),
        )?;
        Ok((optical, law, gradients))
    }
}

impl<P: PressureTemperatureSecondOrderField + ?Sized, C: AirCompositionSecondOrderField + ?Sized>
    PreparedMaterialAirCompositionRefractiveField<'_, '_, P, C>
{
    /// H(n)=sum(n_a H(q_a))+sum(n_ab grad(q_a) grad(q_b)^T), q=[T,P,xw,CO2].
    /// All six mixed material responses survive even for affine state fields.
    pub fn optical_jet(&self, point: Point3) -> Result<IndexJet, AirFieldError> {
        admit_position(point)?;
        let state = self.source.state.pressure_temperature_jet(point)?;
        let composition = self.source.composition.air_composition_jet(point)?;
        let (optical, law, gradients) = self.evaluate(state.sample(), composition.sample)?;
        let diagonal = [
            state.temperature.diagonal_k_per_m2,
            state.pressure_diagonal_pa_per_m2,
            composition.water_diagonal_per_m2,
            composition.co2_diagonal_micromol_per_mol_m2,
        ];
        let mixed = [
            state.temperature.off_diagonal_k_per_m2,
            state.pressure_off_diagonal_pa_per_m2,
            composition.water_off_diagonal_per_m2,
            composition.co2_off_diagonal_micromol_per_mol_m2,
        ];
        let entry = |i: usize, j: usize, state_hessian: [f64; 4]| {
            let mut sum = 0.0;
            for a in 0..4 {
                sum += law.phase_gradient[a] * state_hessian[a]
                    + law.phase_hessian_diagonal[a] * gradients[a][i] * gradients[a][j];
            }
            for (k, (a, b)) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
                .into_iter()
                .enumerate()
            {
                sum += law.phase_hessian_off_diagonal[k]
                    * (gradients[a][i] * gradients[b][j] + gradients[b][i] * gradients[a][j]);
            }
            sum
        };
        let hessian = IndexHessian::try_new(
            core::array::from_fn(|i| entry(i, i, diagonal.map(|h| h[i]))),
            core::array::from_fn(|k| {
                let (i, j) = [(0, 1), (0, 2), (1, 2)][k];
                entry(i, j, mixed.map(|h| h[k]))
            }),
        )
        .map_err(|_| GradientIndexError::NumericalRange)?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<P: PressureTemperatureGradientField + ?Sized, C: AirCompositionGradientField + ?Sized>
    GradientIndexField for MaterialAirCompositionRefractiveField<'_, P, C>
{
    fn sample(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|f| f.optical_sample(point))
            .map_err(field_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.state.is_spatially_uniform() && self.composition.is_spatially_uniform()
    }
}

impl<P: PressureTemperatureGradientField + ?Sized, C: AirCompositionGradientField + ?Sized>
    GradientIndexField for PreparedMaterialAirCompositionRefractiveField<'_, '_, P, C>
{
    fn sample(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(point).map_err(field_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.source.is_spatially_uniform()
    }
}

impl<P: PressureTemperatureSecondOrderField + ?Sized, C: AirCompositionSecondOrderField + ?Sized>
    SecondOrderIndexField for MaterialAirCompositionRefractiveField<'_, P, C>
{
    fn sample_second_order(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|f| f.optical_jet(point))
            .map_err(field_error)
    }
}

impl<P: PressureTemperatureSecondOrderField + ?Sized, C: AirCompositionSecondOrderField + ?Sized>
    SecondOrderIndexField for PreparedMaterialAirCompositionRefractiveField<'_, '_, P, C>
{
    fn sample_second_order(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(point).map_err(field_error)
    }
}

#[cfg(test)]
mod tests;
