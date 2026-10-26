//! Heterogeneous heat storage on the history runner's unchanged physical mesh.
//!
//! Capacity is independent of conductivity: a material's conductivity does not
//! determine its rho*c_p. Values are declared or explicitly resolved from the
//! original material card. Preparation uses the existing lumped P1 capacity
//! rule and captures the matrix and any source domains/receipts; borrowed input
//! tables need not outlive it.

use super::poll;
use crate::{ConductionError, ConductionMesh, ElementMaterials, MaterialId};
use crate::transient::{MaterialHeatCapacity, VolumetricHeatCapacity, backward_euler::BackwardEuler};
use fs_exec::Cx;
use std::collections::BTreeMap;

/// Explicit constant heat-capacity distribution, in J/(m^3 K).
#[derive(Debug, Clone, Copy)]
pub enum ThermalHistoryCapacity<'a> {
    /// One declared value for the whole mesh, retaining the original path.
    Uniform(VolumetricHeatCapacity),
    /// Exactly one declared value per tetrahedron, in the original mesh order.
    PerElement(&'a [VolumetricHeatCapacity]),
    /// Exactly one value per material actually used by `element_materials`.
    /// Input ordering is irrelevant. No missing value inherits a fallback, and
    /// duplicate or unused IDs refuse instead of silently changing a model.
    ByMaterial(&'a [(MaterialId, VolumetricHeatCapacity)]),
    /// Uniform density and specific heat selected from one material card.
    /// Source validity is checked against initial and accepted nodal fields.
    MaterialUniform(&'a MaterialHeatCapacity),
    /// Source-derived rho*c_p using the SAME material IDs as conductivity.
    /// Each source remains distinct; no representative material is averaged in.
    MaterialById(&'a [(MaterialId, MaterialHeatCapacity)]),
}

impl ThermalHistoryCapacity<'_> {
    pub(super) fn prepare<'m>(
        self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: Option<&ElementMaterials>,
    ) -> Result<BackwardEuler<'m>, ConductionError> {
        poll(cx, 0)?;
        match self {
            Self::Uniform(value) => BackwardEuler::uniform(cx, mesh, value),
            Self::PerElement(values) => BackwardEuler::per_element(cx, mesh, values),
            Self::MaterialUniform(source) => BackwardEuler::uniform_material(cx, mesh, source),
            Self::MaterialById(rows) => {
                let materials = materials.ok_or_else(|| invalid(
                    "material-bound capacity requires an explicit element material assignment".into()))?;
                materials.validate_for(mesh)?;
                let mut indexed = BTreeMap::new();
                let mut sources = Vec::with_capacity(rows.len());
                for (slot, (id, source)) in rows.iter().enumerate() {
                    if slot % 512 == 0 { poll(cx, slot)?; }
                    if indexed.insert(*id, slot).is_some() {
                        return Err(invalid(format!("duplicate heat-capacity source for material {id:?}")));
                    }
                    sources.push(source.clone());
                }
                let mut assigned = Vec::with_capacity(mesh.element_count());
                for (element, id) in materials.of_element().iter().enumerate() {
                    if element % 512 == 0 { poll(cx, element)?; }
                    assigned.push(*indexed.get(id).ok_or_else(|| invalid(format!(
                        "element {element} uses material {id:?}, which has no heat-capacity source")))?);
                }
                BackwardEuler::material_assigned(cx, mesh, &sources, &assigned)
            }
            Self::ByMaterial(rows) => {
                let materials = materials.ok_or_else(|| invalid(
                    "material-bound capacity requires an explicit element material assignment".into()))?;
                materials.validate_for(mesh)?;
                let mut indexed = BTreeMap::new();
                for (slot, &(id, value)) in rows.iter().enumerate() {
                    if slot % 512 == 0 { poll(cx, slot)?; }
                    if indexed.insert(id, (value, false)).is_some() {
                        return Err(invalid(format!("duplicate heat-capacity declaration for material {id:?}")));
                    }
                }
                let mut values = Vec::with_capacity(mesh.element_count());
                for (element, id) in materials.of_element().iter().enumerate() {
                    if element % 512 == 0 { poll(cx, element)?; }
                    let (value, used) = indexed.get_mut(id).ok_or_else(|| invalid(format!(
                        "element {element} uses material {id:?}, which has no declared heat capacity")))?;
                    *used = true;
                    values.push(*value);
                }
                for (slot, (id, (_, used))) in indexed.iter().enumerate() {
                    if slot % 512 == 0 { poll(cx, slot)?; }
                    if !*used {
                        return Err(invalid(format!("heat-capacity declaration names unused material {id:?}")));
                    }
                }
                poll(cx, values.len())?;
                BackwardEuler::per_element(cx, mesh, &values)
            }
        }
    }
}

fn invalid(what: String) -> ConductionError {
    ConductionError::Config { parameter: "heat-capacity assignment", what }
}
