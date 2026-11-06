//! Original nodal gas composition on the native tetrahedral mesh.
//!
//! Each borrowed component row is a partial density, in the original molecular
//! model's exact order and declared units. One native barycentric admission
//! supplies every species. This is an accepted continuous P1 reconstruction,
//! not an interpolation of cell averages, a fraction normalization, or another
//! equation of state. The temperature snapshot supplies existing geometry and
//! topology; its temperature is never used as a concentration surrogate.
//!
//! The original molecular composition owner supplies phase/group dispersion,
//! density-to-index gradients and the nonlinear optical Hessian. Piecewise
//! affine partial densities have zero in-element Hessians, but their optical
//! Hessian is generally nonzero. Shared original vertex IDs make concentration
//! and index continuous; native transport retains unsmoothed gradient jumps.

use crate::gradient_index::GradientIndexError;
use crate::material_state::gas_scattering::{
    GasCompositionGradientField, GasCompositionSecondOrderField, GasPartialDensityJet,
    GasPartialDensitySample, GasRefractiveFieldError, MaterialGasCompositionRefractiveField,
};
use crate::material_state::thermo_optic::transient::{
    P1ConductionElement, P1ConductionTemperatureSnapshot, TransientTemperatureFieldError,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_material::gas_mixture::MAX_GAS_MIXTURE_COMPONENTS;
use fs_material::gas_scattering::composition::{MolecularRefractionModel, PartialDensityBasis};

/// Native mixed-gas rays and angular bundles through original element faces.
pub mod transport;

/// Borrowed original component-major partial-density buffers and native geometry.
/// All rows have one value per original vertex; absent species have an explicit
/// zero row. The source identity names the accepted geometry/state/reconstruction.
#[derive(Clone, Copy)]
pub struct P1GasCompositionSnapshot<'a> {
    geometry: P1ConductionTemperatureSnapshot<'a>,
    model: &'a MolecularRefractionModel,
    partial_densities: &'a [&'a [f64]],
    basis: PartialDensityBasis,
    snapshot_identity: ContentHash,
    identity: ContentHash,
}

impl<'a> P1GasCompositionSnapshot<'a> {
    /// Validate the original buffers once, without copying or allocating.
    /// Equal lengths do not establish physical correspondence: the caller
    /// declares the original model ordering and nodal P1 reconstruction here.
    /// Partial densities need not sum to one; all-zero composition is vacuum.
    pub fn try_new(
        cx: &Cx<'_>,
        geometry: P1ConductionTemperatureSnapshot<'a>,
        model: &'a MolecularRefractionModel,
        partial_densities: &'a [&'a [f64]],
        basis: PartialDensityBasis,
        snapshot_identity: ContentHash,
    ) -> Result<Self, GasCompositionSnapshotError> {
        cx.checkpoint()
            .map_err(|_| GasCompositionSnapshotError::Cancelled)?;
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(GasCompositionSnapshotError::MissingSnapshotIdentity);
        }
        if partial_densities.len() != model.components().len() {
            return Err(GasCompositionSnapshotError::ComponentCount {
                expected: model.components().len(),
                found: partial_densities.len(),
            });
        }
        for (component, values) in partial_densities.iter().enumerate() {
            if values.len() != geometry.mesh().vertex_count() {
                return Err(GasCompositionSnapshotError::FieldLength {
                    component,
                    expected: geometry.mesh().vertex_count(),
                    found: values.len(),
                });
            }
            for (vertex, &value) in values.iter().enumerate() {
                if vertex % 512 == 0 {
                    cx.checkpoint()
                        .map_err(|_| GasCompositionSnapshotError::Cancelled)?;
                }
                if !value.is_finite() || value < 0.0 {
                    return Err(GasCompositionSnapshotError::InvalidDensity { component, vertex });
                }
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.p1-gas-composition.v1");
        hash.update(model.content_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        hash.update(&[match basis {
            PartialDensityBasis::MassKgPerM3 => 0,
            PartialDensityBasis::MoleculesPerM3 => 1,
        }]);
        Ok(Self {
            geometry,
            model,
            partial_densities,
            basis,
            snapshot_identity,
            identity: hash.finalize(),
        })
    }

    pub const fn geometry_snapshot(self) -> P1ConductionTemperatureSnapshot<'a> {
        self.geometry
    }
    pub const fn model(self) -> &'a MolecularRefractionModel {
        self.model
    }
    pub const fn nodal_partial_densities(self) -> &'a [&'a [f64]] {
        self.partial_densities
    }
    pub const fn density_basis(self) -> PartialDensityBasis {
        self.basis
    }
    pub const fn snapshot_identity(self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn content_identity(self) -> ContentHash {
        self.identity
    }

    /// Prepare one original element in O(component count), with fixed-capacity
    /// stack storage. Native geometry gradients are reused without inversion.
    pub fn element(
        self,
        element: usize,
    ) -> Result<P1GasCompositionElement<'a>, GasCompositionSnapshotError> {
        Ok(P1GasCompositionElement {
            snapshot: self,
            geometry: self.geometry.element(element)?,
            polynomials: self.polynomials(element)?,
        })
    }

    fn polynomials(
        self,
        element: usize,
    ) -> Result<[DensityPolynomial; MAX_GAS_MIXTURE_COMPONENTS], GasCompositionSnapshotError> {
        let vertices = self.geometry.mesh().complex().tets.get(element).ok_or(
            GasCompositionSnapshotError::Geometry(TransientTemperatureFieldError::InvalidElement {
                element,
            }),
        )?;
        let gradients = self.geometry.mesh().geometry().grads[element];
        let mut polynomials = [DensityPolynomial::ZERO; MAX_GAS_MIXTURE_COMPONENTS];
        for (out, values) in polynomials.iter_mut().zip(self.partial_densities) {
            let values = vertices.map(|vertex| values[vertex as usize]);
            let differences = [
                values[1] - values[0],
                values[2] - values[0],
                values[3] - values[0],
            ];
            // Original nodal differences preserve an exactly uniform species.
            let component = |axis: usize| {
                differences[0].mul_add(
                    gradients[1][axis],
                    differences[1].mul_add(gradients[2][axis], differences[2] * gradients[3][axis]),
                )
            };
            let gradient = Vec3::new(component(0), component(1), component(2));
            if ![gradient.x, gradient.y, gradient.z]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(GasCompositionSnapshotError::NumericalRange);
            }
            *out = DensityPolynomial { values, gradient };
        }
        Ok(polynomials)
    }
}

/// One original element; one geometric query serves every component and jet.
/// It also implements the existing composition field traits for ordinary native
/// refraction, bundle, scene or scattering bindings restricted to this element.
pub struct P1GasCompositionElement<'a> {
    snapshot: P1GasCompositionSnapshot<'a>,
    geometry: P1ConductionElement<'a>,
    polynomials: [DensityPolynomial; MAX_GAS_MIXTURE_COMPONENTS],
}

impl<'a> P1GasCompositionElement<'a> {
    pub const fn snapshot(&self) -> P1GasCompositionSnapshot<'a> {
        self.snapshot
    }
    pub const fn geometry_element(&self) -> P1ConductionElement<'a> {
        self.geometry
    }
    pub fn bind_refraction(
        &self,
    ) -> Result<MaterialGasCompositionRefractiveField<'_, Self>, GasRefractiveFieldError> {
        MaterialGasCompositionRefractiveField::try_new(
            self.snapshot.model,
            self,
            self.snapshot.content_identity(),
        )
    }
    pub fn sample(
        &self,
        position_m: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GasCompositionSnapshotError> {
        let lambda = self.geometry.barycentric(position_m)?;
        sample_rows(&self.polynomials[..self.component_count()], lambda, output)
    }
    pub fn sample_jets(
        &self,
        position_m: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GasCompositionSnapshotError> {
        let lambda = self.geometry.barycentric(position_m)?;
        sample_jet_rows(&self.polynomials[..self.component_count()], lambda, output)
    }
}

impl GasCompositionGradientField for P1GasCompositionElement<'_> {
    fn model_identity(&self) -> ContentHash {
        self.snapshot.model.content_identity()
    }
    fn component_count(&self) -> usize {
        self.snapshot.partial_densities.len()
    }
    fn density_basis(&self) -> PartialDensityBasis {
        self.snapshot.basis
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        self.sample(p, output)
            .map_err(GasCompositionSnapshotError::into_gradient_error)
    }
}
impl GasCompositionSecondOrderField for P1GasCompositionElement<'_> {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        self.sample_jets(p, output)
            .map_err(GasCompositionSnapshotError::into_gradient_error)
    }
}

#[derive(Clone, Copy)]
struct DensityPolynomial {
    values: [f64; 4],
    gradient: Vec3,
}
impl DensityPolynomial {
    const ZERO: Self = Self {
        values: [0.0; 4],
        gradient: Vec3::new(0.0, 0.0, 0.0),
    };

    fn sample(self, lambda: [f64; 4]) -> Result<f64, GasCompositionSnapshotError> {
        if let Some(vertex) = lambda.iter().enumerate().find_map(|(index, &value)| {
            (value == 1.0
                && lambda
                    .iter()
                    .enumerate()
                    .all(|(i, &v)| i == index || v == 0.0))
            .then_some(index)
        }) {
            return Ok(self.values[vertex]);
        }
        // Anchor at the smallest nonnegative node. All increments are
        // nonnegative, including beside a zero-density species boundary;
        // an exactly constant row stays exact without clamping the field.
        let base = self.values.into_iter().fold(f64::INFINITY, f64::min);
        let mut value = base;
        for (weight, node) in lambda.into_iter().zip(self.values) {
            value = weight.mul_add(node - base, value);
        }
        if !value.is_finite() || value < 0.0 {
            return Err(GasCompositionSnapshotError::NumericalRange);
        }
        Ok(value)
    }
}

fn sampled_values(
    polynomials: &[DensityPolynomial],
    lambda: [f64; 4],
    output_len: usize,
) -> Result<[f64; MAX_GAS_MIXTURE_COMPONENTS], GasCompositionSnapshotError> {
    if output_len != polynomials.len() {
        return Err(GasCompositionSnapshotError::OutputCount {
            expected: polynomials.len(),
            found: output_len,
        });
    }
    let mut values = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
    for (value, polynomial) in values.iter_mut().zip(polynomials) {
        *value = polynomial.sample(lambda)?;
    }
    Ok(values)
}
fn sample_rows(
    polynomials: &[DensityPolynomial],
    lambda: [f64; 4],
    output: &mut [GasPartialDensitySample],
) -> Result<(), GasCompositionSnapshotError> {
    let values = sampled_values(polynomials, lambda, output.len())?;
    for ((out, polynomial), density) in output.iter_mut().zip(polynomials).zip(values) {
        *out = GasPartialDensitySample {
            density,
            gradient_per_m: polynomial.gradient,
        };
    }
    Ok(())
}
fn sample_jet_rows(
    polynomials: &[DensityPolynomial],
    lambda: [f64; 4],
    output: &mut [GasPartialDensityJet],
) -> Result<(), GasCompositionSnapshotError> {
    let values = sampled_values(polynomials, lambda, output.len())?;
    for ((out, polynomial), density) in output.iter_mut().zip(polynomials).zip(values) {
        *out = GasPartialDensityJet {
            sample: GasPartialDensitySample {
                density,
                gradient_per_m: polynomial.gradient,
            },
            hessian_diagonal_per_m2: [0.0; 3],
            hessian_off_diagonal_per_m2: [0.0; 3],
        };
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasCompositionSnapshotError {
    MissingSnapshotIdentity,
    ComponentCount {
        expected: usize,
        found: usize,
    },
    FieldLength {
        component: usize,
        expected: usize,
        found: usize,
    },
    InvalidDensity {
        component: usize,
        vertex: usize,
    },
    OutputCount {
        expected: usize,
        found: usize,
    },
    Geometry(TransientTemperatureFieldError),
    NumericalRange,
    Cancelled,
}
impl GasCompositionSnapshotError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Geometry(error) => error.into_gradient_error(),
            Self::NumericalRange => GradientIndexError::NumericalRange,
            Self::Cancelled => GradientIndexError::Cancelled,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<TransientTemperatureFieldError> for GasCompositionSnapshotError {
    fn from(error: TransientTemperatureFieldError) -> Self {
        Self::Geometry(error)
    }
}
impl core::fmt::Display for GasCompositionSnapshotError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original P1 gas partial densities: {self:?}")
    }
}
impl core::error::Error for GasCompositionSnapshotError {}

#[cfg(test)]
mod tests;
