//! Composition-driven bending and focusing from original molecular optics.
//!
//! Partial densities, gradients and Hessians refer to one immutable physical
//! reconstruction. Their source order and units are bound explicitly. The
//! additive polarizability is differentiated before the nonlinear material
//! phase law, retaining refraction even when total density is constant.

use super::*;
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use fs_material::gas_mixture::MAX_GAS_MIXTURE_COMPONENTS;
use fs_material::gas_scattering::composition::{
    CompositionRefractiveIndex, MolecularRefractionModel, MolecularRefractionWavelength,
    PartialDensityBasis,
};

/// One component's partial density and spatial gradient in the declared basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasPartialDensitySample {
    /// kg/m³ or molecules/m³; this is not a normalized fraction.
    pub density: f64,
    /// Spatial derivative of the same partial density, per metre.
    pub gradient_per_m: Vec3,
}

/// Value, gradient and symmetric Hessian of one actual partial-density field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasPartialDensityJet {
    pub sample: GasPartialDensitySample,
    /// `(xx, yy, zz)` second spatial derivatives, per square metre.
    pub hessian_diagonal_per_m2: [f64; 3],
    /// `(xy, xz, yz)` mixed spatial derivatives, per square metre.
    pub hessian_off_diagonal_per_m2: [f64; 3],
}

const UNWRITTEN: GasPartialDensitySample = GasPartialDensitySample {
    density: f64::NAN,
    gradient_per_m: Vec3::new(f64::NAN, f64::NAN, f64::NAN),
};
const UNWRITTEN_JET: GasPartialDensityJet = GasPartialDensityJet {
    sample: UNWRITTEN,
    hessian_diagonal_per_m2: [f64::NAN; 3],
    hessian_off_diagonal_per_m2: [f64::NAN; 3],
};

/// Smooth frozen reconstruction in the prepared molecular model's exact order.
/// Metadata must remain fixed for the binding's lifetime. Every output row must
/// be written, including zero-density species; missing rows are not vacuum.
/// Species absent from the source model require a newly sourced model.
pub trait GasCompositionGradientField {
    fn model_identity(&self) -> ContentHash;
    fn component_count(&self) -> usize;
    fn density_basis(&self) -> PartialDensityBasis;
    fn sample_partial_densities(
        &self,
        position_m: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError>;
}

/// Optional analytic curvature of the same partial densities as the first-order
/// query. No finite differences or smoothing of cell averages are inferred.
pub trait GasCompositionSecondOrderField: GasCompositionGradientField {
    fn sample_partial_density_jets(
        &self,
        position_m: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError>;
}

/// Allocation-free native GRIN field tied to original species and a snapshot.
/// The material's dilute visible-band model and explicit spectral sources apply.
/// Local composition is supplied, not inferred by another equation of state.
/// Phase controls bending; group index uses the full molecular spectral
/// derivative when dispersion is attached. True absorption is zero. Molecular
/// scattering remains separate and is not integrated along these curved paths.
pub struct MaterialGasCompositionRefractiveField<'a, F: GasCompositionGradientField + ?Sized> {
    model: &'a MolecularRefractionModel,
    field: &'a F,
    basis: PartialDensityBasis,
    identity: ContentHash,
}

impl<'a, F: GasCompositionGradientField + ?Sized> MaterialGasCompositionRefractiveField<'a, F> {
    /// Bind exact source/order identity and row count, then retain the caller's
    /// nonzero snapshot/reconstruction identity. An identity records a declared
    /// reconstruction; it does not certify its interpolation or thermodynamics.
    pub fn try_new(
        model: &'a MolecularRefractionModel,
        field: &'a F,
        snapshot_identity: ContentHash,
    ) -> Result<Self, GasRefractiveFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(GasRefractiveFieldError::MissingSnapshotIdentity);
        }
        if field.model_identity() != model.content_identity() {
            return Err(GasRefractiveFieldError::Gas(
                GasScatteringError::StateMismatch,
            ));
        }
        if field.component_count() != model.components().len() {
            return Err(GasRefractiveFieldError::Gas(GasScatteringError::Components));
        }
        let basis = field.density_basis();
        let mut hash = DomainHasher::new("org.frankensim.render.gas-composition-refraction.v1");
        hash.update(model.content_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        hash.update(&[match basis {
            PartialDensityBasis::MassKgPerM3 => 0,
            PartialDensityBasis::MoleculesPerM3 => 1,
        }]);
        Ok(Self {
            model,
            field,
            basis,
            identity: hash.finalize(),
        })
    }

    pub const fn model(&self) -> &MolecularRefractionModel {
        self.model
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Original provider/material refusals remain available for diagnostics.
    pub fn sample_optics(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GasRefractiveFieldError> {
        admit_position(position_m)?;
        let mut buffer = [UNWRITTEN; MAX_GAS_MIXTURE_COMPONENTS];
        let samples = &mut buffer[..self.model.components().len()];
        self.field
            .sample_partial_densities(position_m, samples)
            .map_err(GasRefractiveFieldError::Field)?;
        let spectral = self
            .model
            .at_wavelength(wavelength_nm, self.basis)
            .map_err(GasRefractiveFieldError::Gas)?;
        self.optics(samples.iter(), &spectral)
            .map(|(optical, _, _)| optical)
    }

    fn optics<'b>(
        &self,
        samples: impl Iterator<Item = &'b GasPartialDensitySample>,
        spectral: &MolecularRefractionWavelength<'_>,
    ) -> Result<(GradientIndexSample, CompositionRefractiveIndex, [f64; 3]), GasRefractiveFieldError>
    {
        let mut densities = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
        let mut gradient_l = [0.0; 3];
        for (i, (sample, coefficient)) in samples.zip(spectral.coefficients()).enumerate() {
            densities[i] = sample.density;
            let g = [
                sample.gradient_per_m.x,
                sample.gradient_per_m.y,
                sample.gradient_per_m.z,
            ];
            if g.into_iter().any(|x| !x.is_finite()) {
                return Err(GasRefractiveFieldError::Field(
                    GradientIndexError::InvalidSample,
                ));
            }
            for j in 0..3 {
                gradient_l[j] += coefficient.coefficient * g[j];
            }
        }
        let index = spectral
            .refractive_index(&densities[..self.model.components().len()])
            .map_err(GasRefractiveFieldError::Gas)?;
        let optical = GradientIndexSample::try_new(
            index.phase_index,
            Vec3::new(gradient_l[0], gradient_l[1], gradient_l[2])
                .scale(index.derivative_lorentz_lorenz),
            0.0,
            Some(index.group_index),
        )
        .map_err(GasRefractiveFieldError::Field)?;
        Ok((optical, index, gradient_l))
    }
}

impl<F: GasCompositionSecondOrderField + ?Sized> MaterialGasCompositionRefractiveField<'_, F> {
    /// `H_n = n_L sum(K_i H_i) + n_LL grad(L) grad(L)^T` in either density basis.
    /// One provider query supplies the complete jet with no heap allocation.
    pub fn sample_optical_jet(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GasRefractiveFieldError> {
        admit_position(position_m)?;
        let mut buffer = [UNWRITTEN_JET; MAX_GAS_MIXTURE_COMPONENTS];
        let samples = &mut buffer[..self.model.components().len()];
        self.field
            .sample_partial_density_jets(position_m, samples)
            .map_err(GasRefractiveFieldError::Field)?;
        let spectral = self
            .model
            .at_wavelength(wavelength_nm, self.basis)
            .map_err(GasRefractiveFieldError::Gas)?;
        let (optical, index, g) = self.optics(samples.iter().map(|jet| &jet.sample), &spectral)?;
        let mut diagonal = [0.0; 3];
        let mut off_diagonal = [0.0; 3];
        for (jet, coefficient) in samples.iter().zip(spectral.coefficients()) {
            if jet
                .hessian_diagonal_per_m2
                .into_iter()
                .chain(jet.hessian_off_diagonal_per_m2)
                .any(|x| !x.is_finite())
            {
                return Err(GasRefractiveFieldError::Field(
                    GradientIndexError::InvalidSample,
                ));
            }
            let k = coefficient.coefficient;
            for i in 0..3 {
                diagonal[i] += k * jet.hessian_diagonal_per_m2[i];
                off_diagonal[i] += k * jet.hessian_off_diagonal_per_m2[i];
            }
        }
        let first = index.derivative_lorentz_lorenz;
        let second = index.second_derivative_lorentz_lorenz;
        for (i, (a, b)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            diagonal[i] = first * diagonal[i] + (second * g[i]) * g[i];
            off_diagonal[i] = first * off_diagonal[i] + (second * g[a]) * g[b];
        }
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| GasRefractiveFieldError::Field(GradientIndexError::NumericalRange))?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<F: GasCompositionGradientField + ?Sized> GradientIndexField
    for MaterialGasCompositionRefractiveField<'_, F>
{
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.sample_optics(position, wavelength_nm)
            .map_err(GasRefractiveFieldError::into_gradient_error)
    }
}
impl<F: GasCompositionSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialGasCompositionRefractiveField<'_, F>
{
    fn sample_second_order(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.sample_optical_jet(position, wavelength_nm)
            .map_err(GasRefractiveFieldError::into_gradient_error)
    }
}

fn admit_position(p: Point3) -> Result<(), GasRefractiveFieldError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(GasRefractiveFieldError::Field(
            GradientIndexError::InvalidSample,
        ))
    }
}

#[cfg(test)]
mod tests;
