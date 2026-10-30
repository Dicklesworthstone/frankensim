//! Source material frames on the original full-memory membrane discretization.

use super::{MembraneError, MembraneFacetMaterial, Result, TensorMembrane, bad, poll};
use crate::tensor_material::TensorMaterialBudget;
use fs_exec::CancelGate;
use fs_material::visco_tensor::{TensorProjectionBudget, source::ResolvedTensorMaterial};
use fs_plate::shell::ShellMesh;

/// One facet's complete source, thickness and independently declared orientation.
#[derive(Clone, Debug)]
pub struct MembraneMaterialSource {
    /// Immutable same-card density, equilibrium tensor and all Maxwell branches.
    pub material: ResolvedTensorMaterial,
    /// Actual homogeneous thickness [m]. Neither density nor sound selects it.
    pub thickness_m: f64,
    /// Proper rotation: global_vector = matrix * material_vector.
    /// Its columns are the material axes expressed in global Cartesian axes.
    /// Facet winding does not select a material direction.
    pub material_to_global: [[f64; 3]; 3],
}

/// Explicit choice between the original two physical-inertia discretizations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MembraneMass {
    /// Exact P1 consistent mass, including cross-node velocity coupling.
    Consistent,
    /// Original row-sum lumped mass; constrained nodes do not donate mass.
    Lumped,
}

impl TensorMembrane {
    /// Assemble every original facet from the same source law and density.
    ///
    /// The full tensor is rotated material -> global -> facet BEFORE the native
    /// hereditary plane-stress elimination. Equilibrium and all nonproportional
    /// branch tensors therefore share the actual material axes. The native
    /// membrane owns geometry, quadrature, inertia, history and time integration.
    /// Missing facets, invalid rotations, exceeded budgets and cancellation
    /// refuse without publishing a partly bound body. This adds no prestress,
    /// bending, thermal evolution or acoustic-radiation approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_sources(
        mesh: ShellMesh,
        sources: Vec<MembraneMaterialSource>,
        fixed_dofs: &[usize],
        mass: MembraneMass,
        projection: TensorProjectionBudget,
        dynamics: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        if sources.is_empty()
            || sources.len() != mesh.tris.len()
            || sources.len() > dynamics.max_material_points
            || !(1..=4096).contains(&dynamics.max_material_points)
            || !(1..=super::super::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM)
                .contains(&dynamics.max_state_coordinates)
            || !(1..=256).contains(&projection.max_coordinates)
            || !(3..=4096).contains(&mesh.nodes.len())
        {
            return Err(bad(
                "source membrane needs a complete bounded per-facet material assignment",
            ));
        }
        let mut materials = Vec::with_capacity(sources.len());
        for (element, source) in sources.iter().enumerate() {
            poll(gate)?;
            if !source.thickness_m.is_finite() || source.thickness_m <= 0.0 {
                return Err(bad("source membrane thickness must be finite and positive"));
            }
            let facet = mesh.facet(element).map_err(MembraneError::Geometry)?;
            // Frame rows are the local basis vectors in global coordinates.
            // Thus F*Q maps material components directly into facet components.
            let rotation = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    (0..3).fold(0.0, |sum, k| {
                        sum + facet.frame[i][k] * source.material_to_global[k][j]
                    })
                })
            });
            let law = source
                .material
                .rotated_law(rotation)
                .map_err(MembraneError::MaterialSource)?;
            materials.push(MembraneFacetMaterial {
                law,
                thickness_m: source.thickness_m,
                density_kg_m3: source.material.density_kg_m3(),
            });
        }
        let mut membrane = Self::assemble(
            mesh,
            materials,
            fixed_dofs,
            projection,
            dynamics,
            gate,
            mass == MembraneMass::Consistent,
        )?;
        poll(gate)?;
        membrane.material_sources = Some(sources);
        Ok(membrane)
    }

    /// Complete original per-facet sources and global material orientations.
    /// Authored-law constructors return `None`; they never acquire provenance
    /// merely because their numerical coefficients match a material source.
    pub fn material_sources(&self) -> Option<&[MembraneMaterialSource]> {
        self.material_sources.as_deref()
    }
}
