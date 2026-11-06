//! Original P1 reference inventory in the full tensor-history membrane.

use super::{
    MembraneError, MembraneFacetMaterial, Result, TensorMembrane, bad, source::MembraneMass,
};
use crate::tensor_material::TensorMaterialBudget;
use fs_exec::CancelGate;
use fs_material::visco_tensor::TensorProjectionBudget;
use fs_plate::shell::ShellMesh;
use fs_solid::membrane::inertia::{p1_reference_consistent_mass, p1_reference_row_masses};

impl TensorMembrane {
    /// Prepare the original hereditary membrane with an explicit P1 inventory.
    ///
    /// Each triangle supplies three positive reference areal densities [kg/m²]
    /// in its ORIGINAL stored corner order. The field is linear within that
    /// triangle and may jump across facets sharing a node. The positive
    /// geometric reference area is independent of facet winding. The same
    /// native triangle kernels as `fs_solid::membrane` integrate density with
    /// the original displacement shapes: exact degree-three consistent mass,
    /// or its positive row sums. Supports restrict both matrix axes; they never
    /// move constrained inventory to a free node.
    ///
    /// This field replaces only `MembraneFacetMaterial::density_kg_m3` for
    /// inertia. It is already integrated through thickness and is not multiplied
    /// by thickness again. The original thickness, frame, tensor law, material
    /// quadrature and complete hereditary memory are retained. The caller owns
    /// the physical consistency of those independent authored inputs; this
    /// entry does not infer material stiffness or source authority from density.
    ///
    /// Preparation retains immutable corner inventory and the existing mass
    /// operator. History initialization and all later native steps use the same
    /// canonical momentum, kinetic energy and force inverse. No deformed
    /// geometry is sampled for mass, and no heat or parameter work is invented
    /// by assigning the initial density. This is not a live mass-transfer API.
    /// The original dense exact propagator and all state/work caps still apply.
    ///
    /// # Errors
    /// Original geometry/material/dynamics refusal, incomplete corner coverage,
    /// invalid or unrepresentable reference inertia, or cancellation. Refusal
    /// returns no partial membrane or constitutive history.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_reference_areal_density(
        mesh: ShellMesh,
        materials: Vec<MembraneFacetMaterial>,
        density_kg_m2: Vec<[f64; 3]>,
        fixed_dofs: &[usize],
        mass: MembraneMass,
        projection: TensorProjectionBudget,
        dynamics: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        Self::assemble_with_reference_density(
            mesh,
            materials,
            fixed_dofs,
            projection,
            dynamics,
            gate,
            mass == MembraneMass::Consistent,
            Some(density_kg_m2),
        )
    }

    /// Explicit immutable corner inventory [kg/m²] in original facet order.
    /// `None` means the original scalar density and thickness own each facet's
    /// inventory. These are reference values, not current-geometry densities.
    pub fn reference_corner_areal_density_kg_m2(&self) -> Option<&[[f64; 3]]> {
        self.reference_areal_density.as_deref()
    }
}

pub(super) fn row_masses(
    area: f64,
    material: &MembraneFacetMaterial,
    density: Option<[f64; 3]>,
) -> Result<[f64; 3]> {
    let volume = material.thickness_m * area;
    if let Some(density) = density {
        if !material.thickness_m.is_finite()
            || material.thickness_m <= 0.0
            || !volume.is_finite()
            || volume <= 0.0
        {
            return Err(bad(
                "each facet needs finite positive thickness and material volume",
            ));
        }
        return p1_reference_row_masses(area, density).map_err(MembraneError::Inertia);
    }
    // Preserve the original scalar constructor's multiplication order and
    // refusal, independently of the additive explicit-field preparation.
    let mass = material.density_kg_m3 * volume / 3.0;
    if !material.thickness_m.is_finite()
        || material.thickness_m <= 0.0
        || !material.density_kg_m3.is_finite()
        || material.density_kg_m3 <= 0.0
        || !volume.is_finite()
        || volume <= 0.0
        || !mass.is_finite()
        || mass <= 0.0
    {
        return Err(bad(
            "each facet needs finite positive thickness, density, volume and nodal mass",
        ));
    }
    Ok([mass; 3])
}

pub(super) fn consistent_mass(
    area: f64,
    material: &MembraneFacetMaterial,
    density: Option<[f64; 3]>,
) -> Result<[[f64; 3]; 3]> {
    if let Some(density) = density {
        return p1_reference_consistent_mass(area, density).map_err(MembraneError::Inertia);
    }
    // This is exactly the original constant-density arithmetic, including its
    // multiplication by two after division by twelve on each diagonal.
    let weight = material.density_kg_m3 * (material.thickness_m * area) / 12.0;
    Ok(core::array::from_fn(|a| {
        core::array::from_fn(|b| if a == b { 2.0 * weight } else { weight })
    }))
}

#[cfg(test)]
mod tests;
