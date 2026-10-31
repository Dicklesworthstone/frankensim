//! Consistent reference-geometry surface and volume loads for P1 tetrahedra.
//!
//! Boundary triangles are sorted by global node triple and admitted only with
//! one incident tet. Their normals point away from that tet's opposite vertex,
//! independently of connectivity winding. Interior faces never receive pressure.
//! Face incidence checks do not certify global injectivity or T-junction freedom:
//! the input still requires a conforming, non-overlapping P1 volume mesh.
//!
//! Uniform traction gives A*t/3 per face node; constant element force density
//! gives V*b/4 per tet node. These preserve resultant and first moment on the
//! represented affine elements, to rounding. Positive pressure acts inward on
//! the REFERENCE normal. No follower loading or finite-deformation law is added.

use crate::linear3::{TetAssemblyBudget, TetElasticError, TetLinearElasticProblem,
    TetStaticSolveConfig, TetThermalDisplacementSolution, TetThermalStrainField,
    TetThermalStressReport};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use std::collections::BTreeMap;

/// Typed time/frame histories feeding these same physical load integrals.
pub mod history;
/// Spatially varying P1 pressure, traction and volume fields.
pub mod fields;
/// Sparse inertial trajectories on the same reference mesh and load operators.
pub mod transient;

#[derive(Debug, Clone, PartialEq)]
pub enum LoadError3 {
    Invalid(&'static str),
    Budget,
    Face { nodes: [usize; 3], what: &'static str },
    Solver(TetElasticError),
    Cancelled,
}
impl std::fmt::Display for LoadError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D distributed load refused: {self:?}")
    }
}
impl std::error::Error for LoadError3 {}
impl From<TetElasticError> for LoadError3 {
    fn from(error: TetElasticError) -> Self {
        if error == TetElasticError::Cancelled { Self::Cancelled } else { Self::Solver(error) }
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), LoadError3> { cx.checkpoint().map_err(|_| LoadError3::Cancelled) }
fn finite(v: f64) -> Result<f64, LoadError3> {
    if v.is_finite() { Ok(v) } else { Err(LoadError3::Invalid("nonfinite load or geometry arithmetic")) }
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { std::array::from_fn(|i| a[i] - b[i]) }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn norm(a: [f64; 3]) -> f64 { a[0].hypot(a[1]).hypot(a[2]) }
fn vector(v: [f64; 3]) -> Result<[f64; 3], LoadError3> {
    for x in v { finite(x)?; }
    Ok(v)
}

/// Read-only exterior face. IDs are its indices in `TetLoadGeometry3::faces`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoadFace3 {
    pub nodes: [usize; 3],
    pub centroid_m: [f64; 3],
    pub area_m2: f64,
    pub outward_normal: [f64; 3],
}
#[derive(Debug)]
struct Incidence { face: LoadFace3, positive_side: bool, owners: u8 }

/// Borrowed immutable geometry plus reusable boundary and element-volume data.
#[derive(Debug)]
pub struct TetLoadGeometry3<'m> {
    nodes: &'m [[f64; 3]],
    tets: &'m [[usize; 4]],
    faces: Vec<LoadFace3>,
    volumes: Vec<f64>,
}
impl<'m> TetLoadGeometry3<'m> {
    /// Bound face-incidence work by `maximum_face_entries` before allocation.
    /// Node/tet/quality limits reuse the elastic assembly budget; no matrix is
    /// assembled. Cancellation is checked per node, element and face operation.
    /// Both tet windings are supported, like the existing P1 elasticity owner.
    pub fn new(nodes: &'m [[f64; 3]], tets: &'m [[usize; 4]], budget: TetAssemblyBudget,
        maximum_face_entries: usize, cx: &Cx<'_>) -> Result<Self, LoadError3> {
        poll(cx)?;
        if nodes.is_empty() || tets.is_empty() { return Err(LoadError3::Invalid("empty volume mesh")); }
        if nodes.len() > budget.maximum_nodes || tets.len() > budget.maximum_tetrahedra
            || tets.len().checked_mul(4).is_none_or(|n| n > maximum_face_entries)
            || nodes.len().checked_mul(3).is_none() { return Err(LoadError3::Budget); }
        if !budget.minimum_scaled_jacobian.is_finite() || budget.minimum_scaled_jacobian <= 0.0 {
            return Err(LoadError3::Invalid("positive finite scaled Jacobian floor required"));
        }
        for p in nodes { poll(cx)?; vector(*p)?; }
        let mut incident = vec![false; nodes.len()];
        let mut entries: BTreeMap<[usize; 3], Incidence> = BTreeMap::new();
        let mut volumes = Vec::with_capacity(tets.len());
        for tet in tets {
            poll(cx)?;
            let mut key = *tet;
            key.sort_unstable();
            if key[3] >= nodes.len() || key.windows(2).any(|p| p[0] == p[1]) {
                return Err(LoadError3::Invalid("invalid tetrahedron connectivity"));
            }
            for &i in tet { incident[i] = true; }
            let p = tet.map(|i| nodes[i]);
            let mut h = 0.0_f64;
            for i in 0..4 { for j in i+1..4 { h = h.max(finite(norm(sub(p[i], p[j])))?); } }
            if h == 0.0 { return Err(LoadError3::Invalid("zero-size tetrahedron")); }
            let edge = |i| sub(p[i], p[0]).map(|x| x / h);
            let quality = finite(dot(edge(1), cross(edge(2), edge(3))))?.abs();
            if quality < budget.minimum_scaled_jacobian { return Err(LoadError3::Invalid("tetrahedron below scaled Jacobian floor")); }
            let volume = finite(((quality / 6.0 * h) * h) * h)?;
            if volume <= 0.0 { return Err(LoadError3::Invalid("unrepresentable positive tetrahedron volume")); }
            volumes.push(volume);
            for opposite in 0..4 {
                poll(cx)?;
                let mut ids = [0; 3];
                let mut slot = 0;
                for (i, &node) in tet.iter().enumerate() { if i != opposite { ids[slot] = node; slot += 1; } }
                ids.sort_unstable();
                let [a, b, c] = ids.map(|i| nodes[i]);
                let normal = cross(sub(b, a).map(|x| x/h), sub(c, a).map(|x| x/h));
                let length = finite(norm(normal))?;
                let side = finite(dot(normal, sub(p[opposite], a).map(|x| x/h)))?;
                if length <= 0.0 || side == 0.0 { return Err(LoadError3::Invalid("degenerate face or opposite vertex")); }
                let positive_side = side > 0.0;
                let outward_normal = normal.map(|x| if positive_side { -x/length } else { x/length });
                let area_m2 = finite((0.5 * length * h) * h)?;
                if area_m2 <= 0.0 { return Err(LoadError3::Invalid("unrepresentable positive face area")); }
                let centroid_m = vector(std::array::from_fn(|i| a[i] + (b[i]-a[i])/3.0 + (c[i]-a[i])/3.0))?;
                match entries.entry(ids) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert(Incidence { face: LoadFace3 { nodes: ids, centroid_m, area_m2, outward_normal }, positive_side, owners: 1 });
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        let previous = entry.get_mut();
                        if previous.owners != 1 || previous.positive_side == positive_side {
                            return Err(LoadError3::Face { nodes: ids, what: "nonmanifold face or incident tetrahedra on the same side" });
                        }
                        previous.owners = 2;
                    }
                }
            }
        }
        if incident.iter().any(|v| !v) { return Err(LoadError3::Invalid("mesh has an unreferenced node")); }
        let mut faces = Vec::new();
        for entry in entries.into_values() { poll(cx)?; if entry.owners == 1 { faces.push(entry.face); } }
        poll(cx)?;
        Ok(Self { nodes, tets, faces, volumes })
    }
    pub fn faces(&self) -> &[LoadFace3] { &self.faces }
    pub fn element_volumes_m3(&self) -> &[f64] { &self.volumes }

    /// Integrate one explicitly superposed load case. Each exterior face may be
    /// assigned once; combine multiple laws explicitly before this call. Loads
    /// on fixed nodes are retained, not discarded. Nodal DOFs must be sorted and
    /// unique. Uniform and elementwise body values have units N/m^3, NOT N or m/s^2.
    pub fn assemble<'g>(&'g self, surfaces: &[SurfaceLoad3<'_>], body: VolumeLoad3<'_>,
        nodal: &[(usize, f64)], cx: &Cx<'_>) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
        poll(cx)?;
        let full_dofs = self.nodes.len() * 3;
        if surfaces.len() > self.faces.len() || nodal.len() > full_dofs { return Err(LoadError3::Budget); }
        let body_at = |i| match body {
            VolumeLoad3::None => [0.0; 3], VolumeLoad3::Uniform(v) => v,
            VolumeLoad3::PerElement(v) => v[i],
        };
        if let VolumeLoad3::PerElement(v) = body {
            if v.len() != self.tets.len() { return Err(LoadError3::Invalid("one force-density vector per tetrahedron required")); }
        }
        let mut forces = self.surface_forces(surfaces, cx)?;
        if !matches!(body, VolumeLoad3::None) {
            for (i, tet) in self.tets.iter().enumerate() {
                poll(cx)?;
                let density = vector(body_at(i))?;
                for &node in tet { for c in 0..3 {
                    let contribution = finite((self.volumes[i] / 4.0) * density[c])?;
                    forces[3*node+c] = finite(forces[3*node+c] + contribution)?;
                } }
            }
        }
        self.finish_forces(forces, nodal, cx)
    }

    // Reference surface integration is shared by force-density and
    // material-mass acceleration loading; neither duplicates a pressure law.
    fn surface_forces(&self, surfaces: &[SurfaceLoad3<'_>], cx: &Cx<'_>)
        -> Result<Vec<f64>, LoadError3> {
        poll(cx)?;
        if surfaces.len() > self.faces.len() { return Err(LoadError3::Budget); }
        let mut forces = vec![0.0; self.nodes.len() * 3];
        let mut used = vec![false; self.faces.len()];
        for row in surfaces {
            poll(cx)?;
            let ids = match *row { SurfaceLoad3::Traction { faces, .. } | SurfaceLoad3::Pressure { faces, .. } => faces };
            if ids.is_empty() || ids.len() > self.faces.len() { return Err(LoadError3::Invalid("surface selection must be nonempty and bounded")); }
            for &id in ids {
                poll(cx)?;
                let face = self.faces.get(id).ok_or(LoadError3::Invalid("unknown exterior face index"))?;
                if std::mem::replace(&mut used[id], true) { return Err(LoadError3::Invalid("duplicate surface load assignment")); }
                let traction = vector(match *row {
                    SurfaceLoad3::Traction { pascals, .. } => pascals,
                    SurfaceLoad3::Pressure { pascals, .. } => face.outward_normal.map(|n| -pascals*n),
                })?;
                for &node in &face.nodes { for c in 0..3 {
                    let contribution = finite((face.area_m2 / 3.0) * traction[c])?;
                    forces[3*node+c] = finite(forces[3*node+c] + contribution)?;
                } }
            }
        }
        Ok(forces)
    }

    // Shared nodal superposition and canonical publication for uniform/P1 fields.
    fn finish_forces<'g>(&'g self, mut forces: Vec<f64>, nodal: &[(usize, f64)],
        cx: &Cx<'_>) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
        let full_dofs = self.nodes.len() * 3;
        let mut previous = None;
        for &(dof, force) in nodal {
            poll(cx)?;
            if dof >= full_dofs || previous.is_some_and(|p| dof <= p) { return Err(LoadError3::Invalid("nodal force DOFs must be sorted, unique and in range")); }
            previous = Some(dof);
            forces[dof] = finite(forces[dof] + finite(force)?)?;
        }
        let mut nodal_forces_n = Vec::new();
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 { poll(cx)?; }
            if force != 0.0 { nodal_forces_n.push((dof, force)); }
        }
        poll(cx)?;
        Ok(TetLoadCase3 { geometry: self, nodal_forces_n })
    }
}

#[derive(Debug, Clone, Copy)]
pub enum SurfaceLoad3<'a> {
    Traction { faces: &'a [usize], pascals: [f64; 3] },
    /// Positive pressure compresses the solid on its reference exterior.
    Pressure { faces: &'a [usize], pascals: f64 },
}
#[derive(Debug, Clone, Copy)]
pub enum VolumeLoad3<'a> {
    None,
    Uniform([f64; 3]),
    PerElement(&'a [[f64; 3]]),
}

/// A load bound to the exact immutable geometry on which it was integrated.
#[derive(Debug)]
pub struct TetLoadCase3<'g, 'm> {
    geometry: &'g TetLoadGeometry3<'m>,
    nodal_forces_n: Vec<(usize, f64)>,
}
/// Native displacement/stress plus actual support reactions (not thermal forces).
#[derive(Debug)]
pub struct TetLoadSolution3 {
    pub displacement: TetThermalDisplacementSolution,
    pub stress: TetThermalStressReport,
    /// Kinematic-support force = integral B^T sigma - external force, on fixed DOFs.
    pub support_reactions_n: Vec<[f64; 3]>,
    /// Absolute free-node physical force residual, N; not a discretization bound.
    pub free_force_residual_n: f64,
}
impl TetLoadCase3<'_, '_> {
    pub fn nodal_forces_n(&self) -> &[(usize, f64)] { &self.nodal_forces_n }
    /// Resultant and moment about an explicitly chosen reference point.
    pub fn resultant(&self, origin_m: [f64; 3], cx: &Cx<'_>) -> Result<([f64; 3], [f64; 3]), LoadError3> {
        poll(cx)?;
        vector(origin_m)?;
        let mut force = [0.0; 3];
        let mut moment = [0.0; 3];
        for &(dof, value) in &self.nodal_forces_n {
            poll(cx)?;
            let mut f = [0.0; 3]; f[dof % 3] = value;
            let m = vector(cross(vector(sub(self.geometry.nodes[dof/3], origin_m))?, f))?;
            force[dof % 3] = finite(force[dof % 3] + value)?;
            for c in 0..3 { moment[c] = finite(moment[c] + m[c])?; }
        }
        Ok((force, moment))
    }

    /// Use the existing thermoelastic solve and input-bound stress recovery,
    /// with this exact mechanical load. Thermal strains/materials/supports stay
    /// explicit. A zero thermal state selects ordinary mechanical equilibrium.
    /// No stiffness, Krylov, stress law or thermal approximation is duplicated.
    /// Mesh mutation/reordering refuses; homogeneous component constraints are
    /// inherited unchanged. This is a quasi-static P1 endpoint, not dynamics.
    pub fn solve(&self, problem: &TetLinearElasticProblem<'_>, thermal: TetThermalStrainField<'_>,
        config: TetStaticSolveConfig, mesh_frame: ContentHash, maximum_strain: f64,
        cx: &Cx<'_>) -> Result<TetLoadSolution3, LoadError3> {
        poll(cx)?;
        if problem.nodes_m != self.geometry.nodes || problem.tetrahedra != self.geometry.tets {
            return Err(LoadError3::Invalid("loads belong to a different reference volume mesh"));
        }
        let displacement = problem.solve_thermal_displacement_with_nodal_forces(
            thermal, &self.nodal_forces_n, config, cx)?;
        let stress = problem.recover_thermal_stress_with_nodal_forces(&displacement,
            thermal, &self.nodal_forces_n, mesh_frame, maximum_strain, cx)?;
        self.finish_solution(problem, displacement, stress, cx)
    }

    /// Solve the same mechanical/thermal load with explicit absolute support
    /// motion. Each sorted `(full_xyz_dof, metres)` entry must already be fixed
    /// in the problem; omitted fixed components stay zero. Loads at moving
    /// supports remain in the reaction calculation. The returned displacement
    /// is the FULL field, not the zero-support correction used inside PCG.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_prescribed_dofs(
        &self, problem: &TetLinearElasticProblem<'_>, thermal: TetThermalStrainField<'_>,
        prescribed_m: &[(usize, f64)], maximum_prescribed_displacement_m: f64,
        config: TetStaticSolveConfig, mesh_frame: ContentHash, maximum_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<TetLoadSolution3, LoadError3> {
        poll(cx)?;
        if problem.nodes_m != self.geometry.nodes || problem.tetrahedra != self.geometry.tets {
            return Err(LoadError3::Invalid("loads belong to a different reference volume mesh"));
        }
        let displacement = problem.solve_thermal_displacement_with_prescribed_dofs(
            thermal, &self.nodal_forces_n, prescribed_m, maximum_prescribed_displacement_m,
            config, cx)?;
        let stress = problem.recover_thermal_stress_with_prescribed_dofs(&displacement,
            thermal, &self.nodal_forces_n, prescribed_m, maximum_prescribed_displacement_m,
            mesh_frame, maximum_strain, cx)?;
        self.finish_solution(problem, displacement, stress, cx)
    }

    fn finish_solution(&self, problem: &TetLinearElasticProblem<'_>,
        displacement: TetThermalDisplacementSolution, stress: TetThermalStressReport,
        cx: &Cx<'_>) -> Result<TetLoadSolution3, LoadError3> {
        let mut imbalance = stress.nodal_internal_force_n().to_vec();
        for &(dof, value) in &self.nodal_forces_n { poll(cx)?; imbalance[dof/3][dof%3] = finite(imbalance[dof/3][dof%3] - value)?; }
        let mut support_reactions_n = vec![[0.0; 3]; self.geometry.nodes.len()];
        for &dof in problem.fixed_dofs {
            poll(cx)?;
            support_reactions_n[dof/3][dof%3] = imbalance[dof/3][dof%3];
            imbalance[dof/3][dof%3] = 0.0;
        }
        let mut free_force_residual_n = 0.0_f64;
        for v in imbalance { poll(cx)?; for x in v { free_force_residual_n = finite(free_force_residual_n.hypot(x))?; } }
        poll(cx)?;
        Ok(TetLoadSolution3 { displacement, stress, support_reactions_n, free_force_residual_n })
    }
}
