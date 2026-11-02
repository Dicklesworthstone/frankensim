//! Original multiaxis phase conductivity integrated over P1 enthalpy fields.
//!
//! The caller explicitly declares the scalar source isotropic. Its original
//! caloric object, source paths, fixed coordinates and selected claim survive.
//! h_h is the P1 interpolant of nodal specific enthalpy. The existing scalar
//! moment leaf integrates k(h_h) and lambda_i dk/dh across every crossed source
//! knot; this is not interpolation of nodal k or a sample at mean phase.
//! The existing nodal T(h) interpolation, reference masses, sparse thermal
//! residual, Newton, radiation, adaptive stepping and adjoints remain owners.
//! This fixed-reference model introduces no phase mixture or latent account.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;

use super::super::{
    ASSEMBLY_TILE, EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError,
    adjoint::EnthalpyAdjointError, finite, poll,
};
use crate::{
    ConductionMesh, ConductionProblem, ConductivityModel, ScalarField, ThermalBoundary,
    material::integrate_scalar_path,
};

/// One original phase source interpreted explicitly as isotropic conductivity.
/// The independently declared rho0 is fixed reference mass per mesh volume.
/// A directional scalar source does not establish physical isotropy by itself.
#[derive(Debug)]
pub struct MaterialPhaseTransport<'a> {
    source: &'a PreparedPhaseConductivity<'a>,
    // The original source knots are converted once into the shared leaf's
    // coordinate/value view. Their coordinate remains h [J/kg], never T [K].
    knots: Vec<(f64, f64)>,
    reference_density_kg_m3: f64,
    numerical_basis: ConductivityModel,
    identity: ContentHash,
}

impl<'a> MaterialPhaseTransport<'a> {
    /// Retain one complete original scalar k(T,f) source on its caloric h path.
    /// Native FEM conduction requires strictly positive k; insulating zero-k
    /// charts remain available to the separate finite-volume conduction owner.
    /// No conductivity or density is inferred from a material name.
    pub fn try_new_isotropic(
        cx: &Cx<'_>,
        source: &'a PreparedPhaseConductivity<'a>,
        reference_density_kg_m3: f64,
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        if !reference_density_kg_m3.is_finite() || reference_density_kg_m3 <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "positive fixed reference density required for phase field transport",
            ));
        }
        let mut knots = Vec::with_capacity(source.knots().len());
        for (index, &[h, k]) in source.knots().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if !k.is_finite() || k <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "spatial phase conductivity must be strictly positive",
                ));
            }
            knots.push((h, k));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-conduction.phase-field-transport.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&reference_density_kg_m3.to_bits().to_le_bytes());
        poll(cx, knots.len())?;
        Ok(Self {
            source,
            knots,
            reference_density_kg_m3,
            numerical_basis: ConductivityModel::isotropic_declared(1.0)?,
            identity: hash.finalize(),
        })
    }

    /// Identity of the original source, isotropic interpretation and fixed rho0.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact original source, including its borrowed caloric chart and paths.
    #[must_use]
    pub const fn conductivity_source(&self) -> &'a PreparedPhaseConductivity<'a> {
        self.source
    }

    /// Fixed reference mass per mesh volume, independent of current bulk density.
    #[must_use]
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }

    /// Bind the SAME caloric chart to existing P1 reference-mass storage.
    /// The ordinary forward, endpoint, radiation, adaptive and adjoint APIs
    /// then evaluate the source-integrated conductivity through that owner.
    pub fn prepare_uniform<'m, 's>(
        &'s self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        budget: EnthalpyBudget,
    ) -> Result<EnthalpyBackwardEuler<'m, 's>, EnthalpyError> {
        let mut storage = EnthalpyBackwardEuler::uniform(
            cx,
            mesh,
            self.source.caloric_curve(),
            self.reference_density_kg_m3,
            budget,
        )?;
        storage.sourced_phase = Some(self);
        Ok(storage)
    }

    /// Retain this law's unit numerical basis. An unrelated material tensor,
    /// element override or second phase multiplier cannot multiply physical k.
    pub fn problem<'s>(
        &'s self,
        mesh: &'s ConductionMesh,
        boundary: &'s ThermalBoundary,
        source: &'s ScalarField,
    ) -> ConductionProblem<'s> {
        ConductionProblem {
            mesh,
            boundary,
            material: &self.numerical_basis,
            element_materials: None,
            source,
        }
    }

    /// Element-average physical heat flux, -<k(h_h)> grad(T_h), in W/m².
    /// T_h remains the native P1 interpolant of nodal chart temperatures.
    /// This uses the complete original conductivity, not the unit numerical
    /// basis retained by the problem. It is not a face-conservative recovery.
    pub fn element_heat_flux(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<Vec<[f64; 3]>, EnthalpyError> {
        let (conductivity, _) = self.assemble(cx, mesh, h, false)?;
        let mut temperature = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            temperature.push(
                self.source
                    .caloric_curve()
                    .state_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .temperature_k(),
            );
        }
        let mut flux = Vec::with_capacity(mesh.element_count());
        for (element, (vertices, &k)) in mesh.complex().tets.iter().zip(&conductivity).enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let mut gradient = [0.0; 3];
            for (&vertex, shape_gradient) in vertices.iter().zip(&mesh.geometry().grads[element]) {
                for (value, &shape) in gradient.iter_mut().zip(shape_gradient) {
                    *value = temperature[vertex as usize].mul_add(shape, *value);
                }
            }
            flux.push([
                finite(-k * gradient[0])?,
                finite(-k * gradient[1])?,
                finite(-k * gradient[2])?,
            ]);
        }
        poll(cx, mesh.element_count())?;
        Ok(flux)
    }

    pub(in crate::transient::enthalpy) fn matches_numerical_basis(
        &self,
        problem: ConductionProblem<'_>,
    ) -> bool {
        std::ptr::eq(problem.material, &self.numerical_basis) && problem.element_materials.is_none()
    }

    pub(in crate::transient::enthalpy) fn admit_enthalpies(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
    ) -> Result<(), EnthalpyError> {
        let [lo, hi] = self.source.enthalpy_bounds_j_kg();
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if !value.is_finite() {
                return Err(EnthalpyError::InvalidInput(
                    "finite nodal enthalpy required for phase conductivity",
                ));
            }
            if value < lo || value > hi {
                return Err(EnthalpyError::SourceEnthalpyOutsideDomain {
                    vertex,
                    specific_enthalpy_j_kg: value,
                    lower_j_kg: lo,
                    upper_j_kg: hi,
                });
            }
        }
        Ok(())
    }

    /// Integrate each original k(h_h) path and its consistent nodal derivatives.
    /// The leaf allocates nothing and polls crossed source knots. Element
    /// buffers are O(elements); a forward-only residual skips all derivatives.
    pub(in crate::transient::enthalpy) fn assemble(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
        with_derivatives: bool,
    ) -> Result<(Vec<f64>, Option<Vec<[f64; 4]>>), EnthalpyError> {
        poll(cx, 0)?;
        if h.len() != mesh.vertex_count() {
            return Err(EnthalpyError::InvalidInput(
                "one phase-source enthalpy per mesh vertex required",
            ));
        }
        self.admit_enthalpies(cx, h)?;
        let mut values = Vec::with_capacity(mesh.element_count());
        let mut derivatives = with_derivatives.then(|| Vec::with_capacity(mesh.element_count()));
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let (average, tangent) = integrate_scalar_path(
                Some(cx),
                &self.knots,
                vertices.map(|v| h[v as usize]),
                with_derivatives,
            )?;
            values.push(average);
            if let Some(derivatives) = &mut derivatives {
                derivatives.push(tangent);
            }
        }
        poll(cx, mesh.element_count())?;
        Ok((values, derivatives))
    }

    // A nonuniform P1 integral is differentiable across an interior source
    // knot. A uniform element at a true slope corner has only a semismooth
    // choice. Domain endpoints also prohibit a two-sided source perturbation.
    #[allow(clippy::float_cmp)] // Exact corners delimit classical differentiability.
    pub(in crate::transient::enthalpy) fn admit_adjoint(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<(), EnthalpyAdjointError> {
        self.admit_enthalpies(cx, h)?;
        let [lo, hi] = self.source.enthalpy_bounds_j_kg();
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if value == lo || value == hi {
                return Err(EnthalpyAdjointError::MaterialKink);
            }
        }
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let nodal = vertices.map(|v| h[v as usize]);
            if nodal.iter().any(|&value| value != nodal[0]) {
                continue;
            }
            if let Ok(index) = self.knots.binary_search_by(|k| k.0.total_cmp(&nodal[0])) {
                if index > 0 && index + 1 < self.knots.len() {
                    let a = self.knots[index - 1];
                    let b = self.knots[index];
                    let c = self.knots[index + 1];
                    if (b.1 - a.1) / (b.0 - a.0) != (c.1 - b.1) / (c.0 - b.0) {
                        return Err(EnthalpyAdjointError::MaterialKink);
                    }
                }
            }
        }
        Ok(())
    }
}
