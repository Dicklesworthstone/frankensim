//! Mesh-backed small-strain membrane dynamics with full tensor material memory.
//!
//! The original fs-plate facet owner supplies frames, areas and P1 derivatives.
//! Each triangle has one exact constant-strain quadrature point through a
//! homogeneous thickness, and the original plane-stress history supplies its
//! stress/force response. Translational inertia can retain the exact P1 triangle
//! mass or use explicit row-sum lumping, rho*h*A/3 per vertex. No elastic
//! stiffness is assembled a second time.
//!
//! This is a MEMBRANE discretization: no bending, drilling rotations, geometric
//! stiffness, prestress, finite deformation, or acoustic radiation is supplied.
//! A flat unprestressed membrane has no linear transverse restoring force.
//! The existing dense exact propagator's full-state cap is not bypassed.
use super::{InitialMemory, PlaneStressPoint, TensorMaterialBody, TensorMaterialBudget,
    TensorMaterialError, TensorMaterialHistory};
use fs_exec::CancelGate;
use fs_material::visco_tensor::{TensorProny, TensorProjectionBudget};
use fs_plate::shell::{FacetGeometry, ShellMesh};
use std::collections::BTreeSet;

mod loading;
pub use loading::MembraneEdgeTraction;

/// Explicit homogeneous-through-thickness material on ONE original triangle.
#[derive(Clone, Debug)]
pub struct MembraneFacetMaterial {
    /// Full 3D law expressed in that facet's ex/ey/ez frame, not global axes.
    /// Nonzero Kelvin viscosity retains the plane-stress owner's refusal.
    pub law: TensorProny,
    /// Positive thickness [m]; no geometry-dependent default or averaging.
    pub thickness_m: f64,
    /// Positive volumetric density [kg/m^3], supplied independently of stiffness.
    pub density_kg_m3: f64,
}

/// Geometry, loading, or inherited material/dynamics refusal.
#[derive(Debug)]
pub enum MembraneError {
    /// Invalid input, topology, coordinate map or resource envelope.
    Invalid(&'static str),
    /// Original facet geometry refusal, including its element/node diagnosis.
    Geometry(fs_plate::PlateError),
    /// Original material, exact propagation, or cancellation refusal.
    Dynamics(TensorMaterialError),
}
impl core::fmt::Display for MembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "tensor membrane: {what}"),
            Self::Geometry(e) => write!(f, "tensor membrane geometry: {e}"),
            Self::Dynamics(e) => write!(f, "tensor membrane dynamics: {e}"),
        }
    }
}
impl std::error::Error for MembraneError {}
type Result<T> = core::result::Result<T, MembraneError>;
fn bad(what: &'static str) -> MembraneError { MembraneError::Invalid(what) }
fn poll(gate: &CancelGate) -> Result<()> { super::poll(gate).map_err(MembraneError::Dynamics) }

/// Immutable geometry/material assembly and the original force-driven body.
/// Every free Cartesian displacement and six coordinates of EACH material
/// branch at EACH triangle survive. No mode selection or diagonal loss fit.
pub struct TensorMembrane {
    mesh: ShellMesh,
    materials: Vec<MembraneFacetMaterial>,
    facets: Vec<FacetGeometry>,
    dof_map: Vec<Option<usize>>,
    nodal_mass: Vec<f64>,
    body: TensorMaterialBody,
}
impl TensorMembrane {
    /// Assemble row-sum lumped inertia on existing triangle connectivity.
    /// `fixed_dofs` lists distinct
    /// 3*node+axis addresses (x=0,y=1,z=2), held at zero displacement/velocity.
    /// Partial supports and fully free surfaces are explicit; no rigid modes
    /// are pinned automatically. Duplicate faces and unused vertices refuse.
    /// Geometry intersections and global manifoldness are not certified here.
    ///
    /// The projection work budget is a TOTAL: it is divided evenly among all
    /// facets before any point preparation. Its coordinate ceiling includes
    /// the three local normal strains. Dense state admission precedes all
    /// material condensation. Exhausting either budget never drops a facet.
    pub fn new(mesh: ShellMesh, materials: Vec<MembraneFacetMaterial>, fixed_dofs: &[usize],
        projection: TensorProjectionBudget, dynamics: TensorMaterialBudget,
        gate: &CancelGate) -> Result<Self>
    {
        Self::assemble(mesh,materials,fixed_dofs,projection,dynamics,gate,false)
    }

    /// Assemble exact consistent P1 translational inertia. Each facet contributes
    /// rho*h*A/12 times [[2,1,1],[1,2,1],[1,1,2]] independently to each Cartesian
    /// component. Both matrix axes are restricted to the free DOFs; constrained
    /// neighbors do not donate their inertia to a free vertex. This retains the
    /// cross-node velocity coupling of the consistent finite-element vibration spectrum.
    /// Geometry, constraints, material history and budgets follow `new`.
    pub fn new_with_consistent_mass(mesh:ShellMesh,materials:Vec<MembraneFacetMaterial>,
        fixed_dofs:&[usize],projection:TensorProjectionBudget,dynamics:TensorMaterialBudget,
        gate:&CancelGate)->Result<Self>
    {
        Self::assemble(mesh,materials,fixed_dofs,projection,dynamics,gate,true)
    }

    fn assemble(mesh:ShellMesh,materials:Vec<MembraneFacetMaterial>,fixed_dofs:&[usize],
        projection:TensorProjectionBudget,dynamics:TensorMaterialBudget,gate:&CancelGate,
        consistent:bool)->Result<Self>
    {
        poll(gate)?;
        let nodes = mesh.nodes.len(); let elements = mesh.tris.len();
        if !(3..=4096).contains(&nodes) || elements == 0 || elements != materials.len()
            || !(1..=4096).contains(&dynamics.max_material_points)
            || elements > dynamics.max_material_points
            || !(1..=super::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM).contains(&dynamics.max_state_coordinates)
            || !(1..=256).contains(&projection.max_coordinates)
            || mesh.nodes.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("bounded nonempty mesh, complete per-facet materials and valid budgets required"));
        }
        let mut fixed = vec![false; 3*nodes];
        for &dof in fixed_dofs {
            let entry = fixed.get_mut(dof).ok_or_else(|| bad("constraint names an absent Cartesian DOF"))?;
            if *entry { return Err(bad("duplicate Cartesian constraint")); }
            *entry = true;
        }
        let mut free = 0;
        let dof_map: Vec<_> = fixed.iter().map(|&fixed| if fixed { None } else {
            let i = free; free += 1; Some(i)
        }).collect();
        let memory = materials.iter().try_fold(0_usize, |n, m|
            m.law.branches().len().checked_mul(6).and_then(|b| n.checked_add(b)))
            .ok_or_else(|| bad("membrane memory extent overflow"))?;
        let states = memory.checked_add(2*free).ok_or_else(|| bad("membrane state extent overflow"))?;
        let block = states.checked_add(free+1).and_then(|n| n.checked_mul(2))
            .ok_or_else(|| bad("membrane exponential extent overflow"))?;
        if free == 0 || free > 128 || free+3 > projection.max_coordinates
            || states > dynamics.max_state_coordinates || block > super::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(bad("complete membrane mesh and memory exceed capacity; no coordinates removed"));
        }
        let per_point = TensorProjectionBudget {
            max_scalar_terms: projection.max_scalar_terms/elements, ..projection
        };
        let mut nodal_mass = vec![0.0; nodes];
        let mut facets = Vec::with_capacity(elements);
        let mut unique = BTreeSet::new();
        // Admit all geometry and inertia before constructing any material point.
        for (e, (tri, material)) in mesh.tris.iter().zip(&materials).enumerate() {
            poll(gate)?;
            let facet = mesh.facet(e).map_err(MembraneError::Geometry)?;
            let mut key = *tri; key.sort_unstable();
            if !unique.insert(key) { return Err(bad("duplicate membrane triangle")); }
            let volume = material.thickness_m*facet.area_m2;
            let mass = material.density_kg_m3*volume/3.0;
            if !material.thickness_m.is_finite() || material.thickness_m <= 0.0
                || !material.density_kg_m3.is_finite() || material.density_kg_m3 <= 0.0
                || !volume.is_finite() || volume <= 0.0 || !mass.is_finite() || mass <= 0.0 {
                return Err(bad("each facet needs finite positive thickness, density, volume and nodal mass"));
            }
            for &node in tri { nodal_mass[node] += mass; }
            facets.push(facet);
        }
        if nodal_mass.iter().any(|m| !m.is_finite() || *m <= 0.0)
            || !nodal_mass.iter().sum::<f64>().is_finite() {
            return Err(bad("isolated membrane vertex or unrepresentable assembled mass"));
        }
        let mut points = Vec::with_capacity(elements);
        for ((tri, facet), material) in mesh.tris.iter().zip(&facets).zip(&materials) {
            poll(gate)?;
            let mut columns = vec![[0.0; 3]; free];
            for (local, &node) in tri.iter().enumerate() {
                let [dx, dy] = facet.gradient[local];
                for axis in 0..3 {
                    if let Some(i) = dof_map[3*node+axis] {
                        let ex = facet.frame[0][axis]; let ey = facet.frame[1][axis];
                        columns[i] = [dx*ex, dy*ey, dy*ex+dx*ey];
                    }
                }
            }
            points.push(PlaneStressPoint::new(&material.law, &columns,
                material.thickness_m*facet.area_m2, per_point, gate).map_err(MembraneError::Dynamics)?);
        }
        let body = if consistent {
            // The same constant density/thickness used by the material quadrature
            // makes this degree-two velocity integral exact on each P1 triangle.
            let mut mass=vec![0.0;free*free];
            for ((tri,facet),material) in mesh.tris.iter().zip(&facets).zip(&materials) {
                poll(gate)?;
                let weight=material.density_kg_m3*(material.thickness_m*facet.area_m2)/12.0;
                for (a,&node_a) in tri.iter().enumerate() {
                    for (b,&node_b) in tri.iter().enumerate() {
                        for axis in 0..3 {
                            if let (Some(i),Some(j))=(dof_map[3*node_a+axis],dof_map[3*node_b+axis]) {
                                mass[i*free+j]+=if a==b {2.0*weight} else {weight};
                            }
                        }
                    }
                }
            }
            TensorMaterialBody::from_plane_stress_mass_matrix(free,mass,points,dynamics,gate)
        } else {
            let masses:Vec<_>=dof_map.iter().enumerate()
                .filter_map(|(dof,index)|index.map(|_|nodal_mass[dof/3])).collect();
            TensorMaterialBody::from_plane_stress(masses,points,dynamics,gate)
        }.map_err(MembraneError::Dynamics)?;
        poll(gate)?;
        Ok(Self { mesh, materials, facets, dof_map, nodal_mass, body })
    }
    /// Original mesh; it cannot change after the strain maps have been prepared.
    pub fn mesh(&self) -> &ShellMesh { &self.mesh }
    /// Original per-facet constitutive inputs in facet order.
    pub fn materials(&self) -> &[MembraneFacetMaterial] { &self.materials }
    /// Original facet frames, areas and gradients, in element order.
    pub fn facet_geometry(&self) -> &[FacetGeometry] { &self.facets }
    /// Full 3*node+axis to retained displacement address (None is a fixed DOF).
    pub fn dof_map(&self) -> &[Option<usize>] { &self.dof_map }
    /// Full-mesh row-sum mass at each node, INCLUDING constrained nodes [kg].
    /// These weights are independent of the kinetic discretization; consistent
    /// dynamics use `body().mass_matrix_kg()`, including off-diagonal entries.
    pub fn nodal_masses_kg(&self) -> &[f64] { &self.nodal_mass }
    /// The original exact-ZOH body; it owns all mechanics and material history.
    pub fn body(&self) -> &TensorMaterialBody { &self.body }

    /// Work-conjugate restriction of nodal force [N] to free Cartesian DOFs.
    /// Forces on fixed DOFs do zero work and are not applied to other nodes.
    /// This is not a support-reaction calculation.
    pub fn project_nodal_forces(&self, force: &[[f64; 3]]) -> Result<Vec<f64>> {
        self.check_nodal(force)?;
        Ok(self.dof_map.iter().enumerate().filter_map(|(i, address)|
            address.map(|_| force[i/3][i%3])).collect())
    }
    fn check_nodal(&self, values: &[[f64; 3]]) -> Result<()> {
        if values.len() != self.mesh.nodes.len() || values.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("complete finite Cartesian nodal field required"));
        }
        Ok(())
    }
    fn restrict_state(&self, values: &[[f64; 3]]) -> Result<Vec<f64>> {
        self.check_nodal(values)?;
        for (i, address) in self.dof_map.iter().enumerate() {
            if address.is_none() && values[i/3][i%3] != 0.0 {
                return Err(bad("initial state violates a fixed Cartesian DOF"));
            }
        }
        self.project_nodal_forces(values)
    }
    /// Reconstruct Cartesian displacement or velocity from the retained vector.
    /// Fixed DOFs are exactly zero; no interpolation or nodal reordering occurs.
    pub fn expand_coordinates(&self, values: &[f64]) -> Result<Vec<[f64; 3]>> {
        if values.len() != self.body.coordinates() || values.iter().any(|v| !v.is_finite()) {
            return Err(bad("complete finite retained coordinate vector required"));
        }
        Ok((0..self.mesh.nodes.len()).map(|node| core::array::from_fn(|axis|
            self.dof_map[3*node+axis].map_or(0.0, |i| values[i]))).collect())
    }
    /// Initialize from physical nodal fields and an explicit constitutive prehistory.
    /// Subsequent force-driven steps use TensorMaterialHistory::advance unchanged.
    pub fn history(&self, displacement_m: &[[f64; 3]], velocity_m_s: &[[f64; 3]],
        memory: InitialMemory, gate: &CancelGate) -> Result<TensorMaterialHistory<'_>> {
        poll(gate)?;
        self.body.history(&self.restrict_state(displacement_m)?, &self.restrict_state(velocity_m_s)?,
            memory, gate).map_err(MembraneError::Dynamics)
    }
}

#[cfg(test)]
mod tests;
