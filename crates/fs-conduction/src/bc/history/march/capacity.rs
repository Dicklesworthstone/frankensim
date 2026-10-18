//! Heterogeneous heat storage on the history runner's unchanged physical mesh.
//!
//! Capacity is independent of conductivity: a material's conductivity does not
//! determine its rho*c_p. All values here are explicitly declared, not inferred
//! material-card queries. Preparation uses the existing lumped P1 capacity rule
//! and captures an immutable matrix; borrowed input tables need not outlive it.

use super::poll;
use crate::{ConductionError, ConductionMesh, ElementMaterials, MaterialId};
use crate::transient::{VolumetricHeatCapacity, backward_euler::BackwardEuler};
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
