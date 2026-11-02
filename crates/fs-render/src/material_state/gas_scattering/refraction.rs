//! Material-derived gas focusing from the same immutable density reconstruction.
//!
//! The density Hessian and nonlinear Lorentz--Lorenz chain-rule term both enter
//! the native ray-bundle equations. One provider query returns the complete
//! spatial jet: no finite differences, independent optical fit or per-ray heap
//! allocation is introduced. Molecular scattering still has its separate owner.

use super::*;
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};

/// One consistent density value, gradient and symmetric spatial Hessian.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDensityJet {
    /// Density [kg/m³] and gradient [kg/m⁴] of the same smooth reconstruction.
    pub sample: GasDensityGradientSample,
    /// `(xx, yy, zz)` second spatial derivatives [kg/m⁵].
    pub hessian_diagonal_kg_m5: [f64; 3],
    /// `(xy, xz, yz)` mixed spatial derivatives [kg/m⁵]; symmetry is structural.
    pub hessian_off_diagonal_kg_m5: [f64; 3],
}

/// Optional analytic second derivatives of a frozen physical gas snapshot.
///
/// The returned value and gradient must match [`GasDensityGradientField`] at
/// the same point. Hessians are those of the actual density reconstruction,
/// not of a separately smoothed field. Cell averages alone do not supply this
/// capability. All positions and derivative units are physical SI quantities.
pub trait GasDensitySecondOrderField: GasDensityGradientField {
    fn sample_density_second_order(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityJet, GradientIndexError>;
}

impl<F: GasDensitySecondOrderField + ?Sized> MaterialGasRefractiveField<'_, F> {
    /// Material phase, gradient and Hessian for native focusing and detector
    /// transport, retaining the original gas/provider refusal for diagnostics.
    ///
    /// `H_n = n' H_rho + n'' grad(rho) grad(rho)^T`. The second term remains
    /// present even for an affine density; conversely a zero local gradient
    /// does not remove focusing by a nonzero density Hessian. The returned
    /// optical sample uses the same phase/group law and zero true absorption
    /// as ordinary propagation. This does not supply molecular dispersion,
    /// composition derivatives or finite-aperture caustic regularization.
    pub fn sample_optical_jet(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GasRefractiveFieldError> {
        if ![position_m.x, position_m.y, position_m.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(GasRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let density = self
            .density_field
            .sample_density_second_order(position_m)
            .map_err(GasRefractiveFieldError::Field)?;
        let gradient = density.sample.gradient_kg_m4;
        let g = [gradient.x, gradient.y, gradient.z];
        if g.into_iter()
            .chain(density.hessian_diagonal_kg_m5)
            .chain(density.hessian_off_diagonal_kg_m5)
            .any(|x| !x.is_finite())
        {
            return Err(GasRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let jet = self
            .state
            .refractive_index_jet_at_density(density.sample.density_kg_m3, wavelength_nm)
            .map_err(GasRefractiveFieldError::Gas)?;
        let first = jet.index.derivative_density_m3_kg;
        let second = jet.second_derivative_density_m6_kg2;
        let optical = GradientIndexSample::try_new(
            jet.index.phase_index,
            gradient.scale(first),
            0.0,
            Some(jet.index.phase_index),
        )
        .map_err(GasRefractiveFieldError::Field)?;
        let diagonal = core::array::from_fn(|i| {
            first * density.hessian_diagonal_kg_m5[i] + (second * g[i]) * g[i]
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = pairs[i];
            first * density.hessian_off_diagonal_kg_m5[i] + (second * g[a]) * g[b]
        });
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| GasRefractiveFieldError::Field(GradientIndexError::NumericalRange))?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<F: GasDensitySecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialGasRefractiveField<'_, F>
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

#[cfg(test)]
mod tests;
