//! Native event and bundle continuation for original gas partial densities.
//!
//! A private two-element cache retains species polynomials. One original
//! geometry query selects all species on the same side, including inside the
//! admitted arrival/departure collar. The original molecular law is prepared
//! once per wavelength and supplies both phase/group paths and optical jets.
//! No stage allocates, smooths concentration gradients, renormalizes fractions,
//! or substitutes total density for the individual molecular polarizabilities.

use super::*;
use crate::charts::Ray;
use crate::material_state::thermo_optic::transient::mesh::bundle::{
    P1MeshBundleOptions, P1MeshBundlePropagation,
};
use crate::material_state::thermo_optic::transient::mesh::{
    MAX_INTERNAL_FACES, MeshTemperature, P1ConductionMesh, P1MeshError, P1MeshOptions,
    P1MeshPropagation,
};
use core::cell::{Cell, RefCell};

impl<'a> P1ConductionMesh<'a> {
    /// Trace one wavelength through the original continuous gas composition.
    /// The element must borrow this exact mesh and geometry snapshot. Actual
    /// species dispersion and its original spectral support remain authoritative.
    /// Exterior and duplicated contact-node faces retain native termination.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_gas_composition(
        &self,
        cx: &Cx<'_>,
        element: &P1GasCompositionElement<'a>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, GasCompositionMeshError> {
        admit_options(options, max_distance_m)?;
        let cursor = self.composition_cursor(element)?;
        let source = MeshComposition::new(&cursor, element);
        let original = MaterialGasCompositionRefractiveField::try_new(
            source.snapshot.model,
            &source,
            source.snapshot.content_identity(),
        )?;
        let prepared = original.at_wavelength(wavelength_nm)?;
        self.propagate_prepared(
            cx,
            &cursor,
            &prepared,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
        )
        .map_err(GasCompositionMeshError::Mesh)
    }

    /// Retain native launch-angle derivatives through concentration-gradient
    /// jumps and their crossing-time correction. Zero P1 species Hessians do
    /// not discard the nonlinear Lorentz--Lorenz optical Hessian. Continuous
    /// index at an internal face introduces no artificial Fresnel event.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_gas_composition_bundle(
        &self,
        cx: &Cx<'_>,
        element: &P1GasCompositionElement<'a>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, GasCompositionMeshError> {
        admit_options(options.mesh, max_distance_m)?;
        let cursor = self.composition_cursor(element)?;
        let source = MeshComposition::new(&cursor, element);
        let original = MaterialGasCompositionRefractiveField::try_new(
            source.snapshot.model,
            &source,
            source.snapshot.content_identity(),
        )?;
        let prepared = original.at_wavelength(wavelength_nm)?;
        self.propagate_bundle_prepared(
            cx,
            &cursor,
            &prepared,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
            |previous, next| source.gradient_jump(previous, next),
        )
        .map_err(GasCompositionMeshError::Mesh)
    }

    fn composition_cursor(
        &self,
        element: &P1GasCompositionElement<'a>,
    ) -> Result<MeshTemperature<'a>, GasCompositionMeshError> {
        let source = element.snapshot.geometry;
        if !core::ptr::eq(source.mesh(), self.snapshot.mesh())
            || !core::ptr::eq(
                source.nodal_temperature_k(),
                self.snapshot.nodal_temperature_k(),
            )
        {
            return Err(P1MeshError::ForeignSnapshot.into());
        }
        Ok(MeshTemperature {
            current: Cell::new(element.geometry),
            previous: Cell::new(None),
        })
    }
}

fn admit_options(options: P1MeshOptions, distance_m: f64) -> Result<(), GasCompositionMeshError> {
    options.event.validate().map_err(P1MeshError::from)?;
    if !distance_m.is_finite()
        || distance_m < 0.0
        || options.max_internal_faces > MAX_INTERNAL_FACES
    {
        return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
    }
    Ok(())
}

struct CachedPolynomials {
    element: usize,
    rows: [DensityPolynomial; MAX_GAS_MIXTURE_COMPONENTS],
}
struct PolynomialCache {
    current: CachedPolynomials,
    previous: Option<CachedPolynomials>,
}
struct MeshComposition<'cursor, 'source> {
    cursor: &'cursor MeshTemperature<'source>,
    snapshot: P1GasCompositionSnapshot<'source>,
    cache: RefCell<PolynomialCache>,
}

impl<'cursor, 'source> MeshComposition<'cursor, 'source> {
    fn new(
        cursor: &'cursor MeshTemperature<'source>,
        element: &P1GasCompositionElement<'source>,
    ) -> Self {
        Self {
            cursor,
            snapshot: element.snapshot,
            cache: RefCell::new(PolynomialCache {
                current: CachedPolynomials {
                    element: element.geometry.element_index(),
                    rows: element.polynomials,
                },
                previous: None,
            }),
        }
    }

    // Borrow the bounded cache instead of copying all 128 possible component
    // polynomials at every ODE stage. Cache changes never choose a geometric
    // side: that choice remains solely with the original native cursor.
    fn with_polynomials<R>(
        &self,
        element: usize,
        sample: impl FnOnce(&[DensityPolynomial]) -> Result<R, GasCompositionSnapshotError>,
    ) -> Result<R, GradientIndexError> {
        let mut cache = self.cache.borrow_mut();
        if cache.current.element != element {
            if cache
                .previous
                .as_ref()
                .is_some_and(|previous| previous.element == element)
            {
                let PolynomialCache { current, previous } = &mut *cache;
                if let Some(previous) = previous.as_mut() {
                    core::mem::swap(current, previous);
                }
            } else {
                let rows = self
                    .snapshot
                    .polynomials(element)
                    .map_err(GasCompositionSnapshotError::into_gradient_error)?;
                let previous =
                    core::mem::replace(&mut cache.current, CachedPolynomials { element, rows });
                cache.previous = Some(previous);
            }
        }
        sample(&cache.current.rows[..self.component_count()])
            .map_err(GasCompositionSnapshotError::into_gradient_error)
    }

    fn gradient_jump(
        &self,
        previous: P1ConductionElement<'source>,
        next: P1ConductionElement<'source>,
    ) -> Result<bool, GradientIndexError> {
        let gradients: [Vec3; MAX_GAS_MIXTURE_COMPONENTS] =
            self.with_polynomials(previous.element_index(), |rows| {
                Ok(core::array::from_fn(|i| {
                    rows.get(i)
                        .map_or(Vec3::new(0.0, 0.0, 0.0), |row| row.gradient)
                }))
            })?;
        self.with_polynomials(next.element_index(), |rows| {
            Ok(rows
                .iter()
                .zip(gradients)
                .any(|(next, previous)| next.gradient != previous))
        })
    }
}

impl GasCompositionGradientField for MeshComposition<'_, '_> {
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
        let (element, lambda) = self.cursor.sample_barycentric(p)?;
        self.with_polynomials(element.element_index(), |rows| {
            sample_rows(rows, lambda, output)
        })
    }
}
impl GasCompositionSecondOrderField for MeshComposition<'_, '_> {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        let (element, lambda) = self.cursor.sample_barycentric(p)?;
        self.with_polynomials(element.element_index(), |rows| {
            sample_jet_rows(rows, lambda, output)
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasCompositionMeshError {
    Mesh(P1MeshError),
    Optical(GasRefractiveFieldError),
}
impl From<P1MeshError> for GasCompositionMeshError {
    fn from(error: P1MeshError) -> Self {
        Self::Mesh(error)
    }
}
impl From<GasRefractiveFieldError> for GasCompositionMeshError {
    fn from(error: GasRefractiveFieldError) -> Self {
        Self::Optical(error)
    }
}
impl core::fmt::Display for GasCompositionMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native P1 gas composition transport: {self:?}")
    }
}
impl core::error::Error for GasCompositionMeshError {}
