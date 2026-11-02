//! One source-moment owner for uniform and assigned material phase fields.
use std::sync::Arc;

use fs_exec::Cx;

use super::super::super::{
    ASSEMBLY_TILE, EnthalpyError, adjoint::EnthalpyAdjointError, finite, poll,
};
use super::MaterialPhaseTransport;
use crate::{
    ConductionMesh, ConductionProblem, ConductivityModel, material::integrate_scalar_path,
};

/// Assignments share the heterogeneous storage owner's admitted immutable maps.
/// Every element's nodes have exactly its material identity; no source is ever
/// queried at a temperature, enthalpy or fraction averaged across materials.
#[derive(Debug)]
pub(in crate::transient::enthalpy) enum PhaseFieldTransport<'c> {
    Uniform(&'c MaterialPhaseTransport<'c>),
    Assigned {
        materials: &'c [MaterialPhaseTransport<'c>],
        element_material_ids: Arc<[usize]>,
        vertex_material_ids: Arc<[usize]>,
    },
}

impl<'c> PhaseFieldTransport<'c> {
    pub(in crate::transient::enthalpy) fn material(
        &self,
        id: usize,
    ) -> Option<&'c MaterialPhaseTransport<'c>> {
        match self {
            Self::Uniform(material) => (id == 0).then_some(*material),
            Self::Assigned { materials, .. } => (*materials).get(id),
        }
    }

    fn for_vertex(&self, vertex: usize) -> &'c MaterialPhaseTransport<'c> {
        match self {
            Self::Uniform(material) => *material,
            Self::Assigned {
                materials,
                vertex_material_ids,
                ..
            } => &materials[vertex_material_ids[vertex]],
        }
    }

    fn for_element(&self, element: usize) -> &'c MaterialPhaseTransport<'c> {
        match self {
            Self::Uniform(material) => *material,
            Self::Assigned {
                materials,
                element_material_ids,
                ..
            } => &materials[element_material_ids[element]],
        }
    }

    pub(in crate::transient::enthalpy) fn numerical_basis(&self) -> &'c ConductivityModel {
        match self {
            Self::Uniform(material) => &material.numerical_basis,
            // The heterogeneous constructor admits a nonempty material table.
            Self::Assigned { materials, .. } => &materials[0].numerical_basis,
        }
    }

    pub(in crate::transient::enthalpy) fn matches_numerical_basis(
        &self,
        problem: ConductionProblem<'_>,
    ) -> bool {
        std::ptr::eq(problem.material, self.numerical_basis())
            && problem.element_materials.is_none()
    }

    pub(in crate::transient::enthalpy) fn admit_enthalpies(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
    ) -> Result<(), EnthalpyError> {
        poll(cx, 0)?;
        if let Self::Assigned {
            vertex_material_ids,
            ..
        } = self
        {
            if h.len() != vertex_material_ids.len() {
                return Err(EnthalpyError::InvalidInput(
                    "phase field must retain the original vertex assignment",
                ));
            }
        }
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if !value.is_finite() {
                return Err(EnthalpyError::InvalidInput(
                    "finite nodal enthalpy required for phase conductivity",
                ));
            }
            let [lo, hi] = self.for_vertex(vertex).source.enthalpy_bounds_j_kg();
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

    /// Integrate original k(h_h) and its consistent derivatives in one pass.
    /// No material table is copied into per-element objects. Residual-only
    /// work skips derivative moments; allocated buffers remain O(elements).
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
        if let Self::Assigned {
            element_material_ids,
            ..
        } = self
        {
            if element_material_ids.len() != mesh.element_count() {
                return Err(EnthalpyError::InvalidInput(
                    "phase field must retain the original element assignment",
                ));
            }
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
                &self.for_element(element).knots,
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

    pub(in crate::transient::enthalpy) fn element_heat_flux(
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
                self.for_vertex(vertex)
                    .source
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

    // Interior source knots have zero volume in a nonuniform P1 field. At a
    // uniform true source corner, retain the existing classical-adjoint gate.
    // Each material's source-domain endpoints have their own two-sided gate.
    #[allow(clippy::float_cmp)]
    pub(in crate::transient::enthalpy) fn admit_adjoint(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<(), EnthalpyAdjointError> {
        self.admit_enthalpies(cx, h)?;
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let [lo, hi] = self.for_vertex(vertex).source.enthalpy_bounds_j_kg();
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
            let knots = &self.for_element(element).knots;
            if let Ok(index) = knots.binary_search_by(|k| k.0.total_cmp(&nodal[0])) {
                if index > 0 && index + 1 < knots.len() {
                    let a = knots[index - 1];
                    let b = knots[index];
                    let c = knots[index + 1];
                    if (b.1 - a.1) / (b.0 - a.0) != (c.1 - b.1) / (c.0 - b.0) {
                        return Err(EnthalpyAdjointError::MaterialKink);
                    }
                }
            }
        }
        Ok(())
    }
}
