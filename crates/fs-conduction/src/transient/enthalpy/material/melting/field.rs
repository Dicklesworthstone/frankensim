//! Original sharp-melting sources on uniform or assigned native material nodes.
use std::sync::Arc;

use fs_exec::Cx;

use super::super::super::{ASSEMBLY_TILE, EnthalpyError, finite, poll};
use super::{MaterialMeltingTransport, MeltingConductivitySample};
use crate::{ConductionMesh, ConductionProblem, ConductivityModel};

/// The existing heterogeneous constructor guarantees that every element has
/// one material identity. Borrow its admitted vertex map; do not copy source
/// records into elements or average different caloric enthalpy datums.
#[derive(Debug)]
pub(in crate::transient::enthalpy) enum MeltingFieldTransport<'c> {
    Uniform(&'c MaterialMeltingTransport),
    Assigned {
        materials: &'c [MaterialMeltingTransport],
        vertex_material_ids: Arc<[usize]>,
    },
}

impl<'c> MeltingFieldTransport<'c> {
    pub(in crate::transient::enthalpy) fn material(
        &self,
        id: usize,
    ) -> Option<&'c MaterialMeltingTransport> {
        match self {
            Self::Uniform(material) => (id == 0).then_some(*material),
            Self::Assigned { materials, .. } => materials.get(id),
        }
    }

    fn for_vertex(&self, vertex: usize) -> &'c MaterialMeltingTransport {
        match self {
            Self::Uniform(material) => material,
            Self::Assigned {
                materials,
                vertex_material_ids,
            } => &materials[vertex_material_ids[vertex]],
        }
    }

    pub(in crate::transient::enthalpy) fn numerical_basis(&self) -> &'c ConductivityModel {
        // Assigned storage admits a nonempty original material table.
        &self
            .material(0)
            .expect("admitted nonempty melting table")
            .numerical_basis
    }

    pub(in crate::transient::enthalpy) fn matches_numerical_basis(
        &self,
        problem: ConductionProblem<'_>,
    ) -> bool {
        self.material(0)
            .expect("admitted nonempty melting table")
            .matches_numerical_basis(problem)
    }

    pub(in crate::transient::enthalpy) fn sample_for_vertex(
        &self,
        h: f64,
        vertex: usize,
    ) -> Result<MeltingConductivitySample, EnthalpyError> {
        self.for_vertex(vertex).sample_for_vertex(h, vertex)
    }

    pub(in crate::transient::enthalpy) fn conductivity_kink_at(
        &self,
        h: f64,
        vertex: usize,
    ) -> Result<bool, EnthalpyError> {
        self.for_vertex(vertex).conductivity_kink_at(h, vertex)
    }

    pub(in crate::transient::enthalpy) fn element_heat_flux(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<Vec<[f64; 3]>, EnthalpyError> {
        poll(cx, 0)?;
        if h.len() != mesh.vertex_count()
            || matches!(self,
            Self::Assigned { vertex_material_ids, .. } if vertex_material_ids.len() != h.len())
        {
            return Err(EnthalpyError::InvalidInput(
                "melting heat flux requires the original nodal material layout",
            ));
        }
        let mut conductivity = Vec::with_capacity(h.len());
        let mut temperature = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let source = self.for_vertex(vertex);
            conductivity.push(source.sample_for_vertex(value, vertex)?.conductivity_w_m_k);
            temperature.push(
                source
                    .caloric
                    .curve()
                    .temperature_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?,
            );
        }
        let mut flux = Vec::with_capacity(mesh.element_count());
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            // Same declared P1 nodal conductivity average as StepContext.
            let k: f64 = vertices
                .iter()
                .map(|&v| conductivity[v as usize] / 4.0)
                .sum();
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
}
