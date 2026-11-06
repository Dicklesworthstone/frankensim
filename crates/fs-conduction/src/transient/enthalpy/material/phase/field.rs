//! One source-moment owner for uniform and assigned material phase fields.
use std::sync::Arc;

use fs_exec::Cx;

use super::super::super::{
    ASSEMBLY_TILE, EnthalpyError, adjoint::EnthalpyAdjointError, finite, poll,
};
use super::MaterialPhaseTransport;
use crate::{
    ConductionMesh, ConductionProblem, ConductivityModel,
    material::{ElementConductivity, integrate_caloric_path},
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
    /// The scalar fast path requires explicit isotropy throughout the material
    /// table, including unused entries and the first entry owning the basis.
    pub(in crate::transient::enthalpy) fn is_isotropic(&self) -> bool {
        match self {
            Self::Uniform(material) => material.principal_axes.is_none(),
            Self::Assigned { materials, .. } => materials
                .iter()
                .all(|material| material.principal_axes.is_none()),
        }
    }

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
            let [lo, hi] = self.for_vertex(vertex).enthalpy_bounds_j_kg();
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

    /// Preserve compact scalar buffers for explicitly isotropic materials.
    /// Integrate one original source path per element without staging full
    /// tensors or allocating an unused derivative array.
    pub(in crate::transient::enthalpy) fn assemble_isotropic(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
        with_derivatives: bool,
    ) -> Result<(Vec<f64>, Option<Vec<[f64; 4]>>), EnthalpyError> {
        self.admit_field(cx, mesh, h)?;
        if !self.is_isotropic() {
            return Err(EnthalpyError::InvalidInput(
                "scalar phase assembly requires every material table entry to be explicitly isotropic",
            ));
        }
        let mut values = Vec::with_capacity(mesh.element_count());
        let mut derivatives = with_derivatives.then(|| Vec::with_capacity(mesh.element_count()));
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let material = self.for_element(element);
            let (average, tangent) = integrate_caloric_path(
                Some(cx),
                &material.knots[0],
                &material.heat_capacities[0],
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

    /// Integrate original K(h_h) and its consistent enthalpy derivatives.
    /// No material table is copied into per-element objects. Residual-only
    /// consumers use [`Self::assemble_tensors`] to avoid derivative storage.
    pub(in crate::transient::enthalpy) fn assemble(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
        with_derivatives: bool,
    ) -> Result<Vec<ElementConductivity>, EnthalpyError> {
        self.map_elements(cx, mesh, h, with_derivatives, |moments| moments)
    }

    /// Complete physical tensor means with no allocated derivative buffer.
    pub(in crate::transient::enthalpy) fn assemble_tensors(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<Vec<[[f64; 3]; 3]>, EnthalpyError> {
        self.map_elements(cx, mesh, h, false, |moments| moments.tensor)
    }

    fn admit_field(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<(), EnthalpyError> {
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
        self.admit_enthalpies(cx, h)
    }

    fn map_elements<T>(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
        with_derivatives: bool,
        mut output: impl FnMut(ElementConductivity) -> T,
    ) -> Result<Vec<T>, EnthalpyError> {
        self.admit_field(cx, mesh, h)?;
        let mut values = Vec::with_capacity(mesh.element_count());
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            values.push(output(self.for_element(element).element_moments(
                cx,
                vertices.map(|v| h[v as usize]),
                with_derivatives,
            )?));
        }
        poll(cx, mesh.element_count())?;
        Ok(values)
    }

    pub(in crate::transient::enthalpy) fn element_heat_flux(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<Vec<[f64; 3]>, EnthalpyError> {
        self.admit_field(cx, mesh, h)?;
        let mut temperature = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            temperature.push(
                self.for_vertex(vertex)
                    .source
                    .caloric_curve()
                    .values_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .temperature_k,
            );
        }
        let mut flux = Vec::with_capacity(mesh.element_count());
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let k = self
                .for_element(element)
                .element_moments(cx, vertices.map(|v| h[v as usize]), false)?
                .tensor;
            let mut gradient = [0.0; 3];
            for (&vertex, shape_gradient) in vertices.iter().zip(&mesh.geometry().grads[element]) {
                for (value, &shape) in gradient.iter_mut().zip(shape_gradient) {
                    *value = temperature[vertex as usize].mul_add(shape, *value);
                }
            }
            flux.push([
                finite(-k[0][0].mul_add(
                    gradient[0],
                    k[0][1].mul_add(gradient[1], k[0][2] * gradient[2]),
                ))?,
                finite(-k[1][0].mul_add(
                    gradient[0],
                    k[1][1].mul_add(gradient[1], k[1][2] * gradient[2]),
                ))?,
                finite(-k[2][0].mul_add(
                    gradient[0],
                    k[2][1].mul_add(gradient[1], k[2][2] * gradient[2]),
                ))?,
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
            let [lo, hi] = self.for_vertex(vertex).enthalpy_bounds_j_kg();
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
            let material = self.for_element(element);
            for (knots, capacities) in material.knots[..material.conductivity_sources().len()]
                .iter().zip(&material.heat_capacities) {
                if let Ok(index) = knots.binary_search_by(|k| k.0.total_cmp(&nodal[0])) {
                    if index > 0 && index + 1 < knots.len() {
                        let a = knots[index - 1];
                        let b = knots[index];
                        let c = knots[index + 1];
                        let different = if let (Some(left), Some(right)) =
                            (capacities[index - 1], capacities[index]) {
                            // Compare the ORIGINAL k(T) slope and shared Cp.
                            // Reconstructing dk/dh from two h-chord products
                            // can manufacture a kink through rounding alone.
                            let temperature = |h| material.source.caloric_curve()
                                .temperature_at_specific_enthalpy(h)
                                .map_err(|source| EnthalpyError::Phase {
                                    vertex: vertices[0] as usize, source,
                                });
                            let (ta, tb, tc) = (temperature(a.0)?, temperature(b.0)?, temperature(c.0)?);
                            left[1] != right[0]
                                || (b.1 - a.1) / (tb - ta) != (c.1 - b.1) / (tc - tb)
                        } else {
                            let left_scale = capacities[index - 1].map_or(1.0, |[ca, cb]| {
                                0.5 * (ca / cb) + 0.5
                            });
                            let right_scale = capacities[index].map_or(1.0, |[ca, cb]| {
                                0.5 + 0.5 * (cb / ca)
                            });
                            ((b.1 - a.1) / (b.0 - a.0)) * left_scale
                                != ((c.1 - b.1) / (c.0 - b.0)) * right_scale
                        };
                        if different {
                            return Err(EnthalpyAdjointError::MaterialKink);
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

impl MaterialPhaseTransport<'_> {
    fn element_moments(
        &self,
        cx: &Cx<'_>,
        nodal: [f64; 4],
        with_derivatives: bool,
    ) -> Result<ElementConductivity, EnthalpyError> {
        let mut moments = ElementConductivity {
            tensor: [[0.0; 3]; 3],
            derivatives: [[[0.0; 3]; 3]; 4],
        };
        let Some(axes) = self.principal_axes else {
            // Preserve the single scalar integration and its exact diagonal
            // entries for the original isotropic API.
            let (average, tangent) =
                integrate_caloric_path(Some(cx), &self.knots[0], &self.heat_capacities[0], nodal, with_derivatives)?;
            for i in 0..3 {
                moments.tensor[i][i] = average;
                if with_derivatives {
                    for c in 0..4 {
                        moments.derivatives[c][i][i] = tangent[c];
                    }
                }
            }
            return Ok(moments);
        };
        for ((axis, knots), capacities) in axes.iter().zip(&self.knots).zip(&self.heat_capacities) {
            let (average, tangent) =
                integrate_caloric_path(Some(cx), knots, capacities, nodal, with_derivatives)?;
            for i in 0..3 {
                for j in i..3 {
                    let value = finite((average * axis[i]).mul_add(axis[j], moments.tensor[i][j]))?;
                    moments.tensor[i][j] = value;
                    moments.tensor[j][i] = value;
                    if with_derivatives {
                        for c in 0..4 {
                            let value = finite(
                                (tangent[c] * axis[i])
                                    .mul_add(axis[j], moments.derivatives[c][i][j]),
                            )?;
                            moments.derivatives[c][i][j] = value;
                            moments.derivatives[c][j][i] = value;
                        }
                    }
                }
            }
        }
        Ok(moments)
    }
}
