//! Three-dimensional current-free magnetostatics on tetrahedral P1 elements.
//!
//! H = -grad(phi), B = mu H + Br, div(B)=0. The original tetrahedral Whitney
//! Hodge supplies K=d0^T M_mu d0. Remanence contributes its complete weak load,
//! including interface and exterior terms. Material sources, recoil directions
//! and temperature are fixed. See the scalar-potential formulation in
//! https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.51.html.
//!
//! This is a full 3D field on CURRENT-FREE domains with a single-valued magnetic
//! potential. Free currents, nonzero circulation around linked currents,
//! induction, nonlinear saturation and moving media require other formulations.
//! Exterior truncation is explicit: unlisted exterior faces have B.n=0, while
//! supplied potential or normal-flux traces define the actual finite problem.
//! A finite box is not an infinite-domain or far-field accuracy certificate.
//!
//! The supplied tetrahedra must form a conforming, non-overlapping volume mesh.
//! Incidence admission checks duplicate cells, nonmanifold faces and same-side
//! neighbors; it does not certify absence of intersections between remote cells.
//! P1 H/B are cellwise constant. Interface normal-flux jumps are reported rather
//! than presenting the raw field as an exactly equilibrated H(div) field.

use crate::planar::permanent::MagneticEnergyReference;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
    whitney::{ElementGeometry, ElementGeometryError, try_element_geometry},
};
use fs_material::{
    field_storage::{
        FieldStorageError, MAGNETIC_PERMEABILITY_DIMS, linear_magnetic_permeability,
        permanent::{RecoilOperatingPoint, ResolvedRecoilMaterial},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TetComplex;
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::Precond};
use std::collections::BTreeSet;

/// Material, geometry, boundary or numerical refusal; no partial field returns.
#[derive(Debug)]
pub enum SpatialMagneticError {
    Invalid(&'static str),
    Material(FieldStorageError),
    Geometry(ElementGeometryError),
    Assembly(WeightedError),
    Recoil {
        cell: usize,
        error: FieldStorageError,
    },
    Nonconvergence(SolveReport),
    Cancelled,
}
impl core::fmt::Display for SpatialMagneticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial magnetic field: {self:?}")
    }
}
impl std::error::Error for SpatialMagneticError {}
impl From<FieldStorageError> for SpatialMagneticError {
    fn from(e: FieldStorageError) -> Self {
        Self::Material(e)
    }
}
fn invalid(s: &'static str) -> SpatialMagneticError {
    SpatialMagneticError::Invalid(s)
}
fn finite(x: f64) -> Result<f64, SpatialMagneticError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("nonfinite magnetic arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpatialMagneticError> {
    if cancelled() {
        Err(SpatialMagneticError::Cancelled)
    } else {
        Ok(())
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|j| a[j] - b[j])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}

/// Original induced or recoil source, with its distinct constitutive meaning.
#[derive(Clone, Debug, PartialEq)]
pub enum SpatialMagneticMaterial {
    Linear(ResolvedMaterialStatePoint),
    Recoil(SpatialRecoilMaterial),
}

/// A source-resolved magnet with explicit unit direction in the mesh frame.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialRecoilMaterial {
    material: ResolvedRecoilMaterial,
    direction: [f64; 3],
    remanence_t: [f64; 3],
}
impl SpatialRecoilMaterial {
    pub fn new(
        material: ResolvedRecoilMaterial,
        direction: [f64; 3],
    ) -> Result<Self, SpatialMagneticError> {
        let length = norm(direction);
        if !length.is_finite() || (length - 1.0).abs() > 1e-12 {
            return Err(invalid(
                "3D remanence direction must be a finite unit vector",
            ));
        }
        let direction = direction.map(|x| x / length);
        let remanence_t = direction.map(|x| x * material.remanence_t());
        for x in remanence_t {
            finite(x)?;
            finite(x / material.permeability_h_m())?;
        }
        Ok(Self {
            material,
            direction,
            remanence_t,
        })
    }
    pub const fn material(&self) -> &ResolvedRecoilMaterial {
        &self.material
    }
    pub const fn direction(&self) -> [f64; 3] {
        self.direction
    }
    pub const fn remanence_t(&self) -> [f64; 3] {
        self.remanence_t
    }
}

/// Dirichlet trace on one EXTERIOR face; values follow its sorted vertex order.
/// Magnetic scalar potential has units A, not volts or vector-potential units.
#[derive(Clone, Copy, Debug)]
pub struct MagneticPotentialFace {
    pub face: usize,
    pub potential_a: [f64; 3],
}
/// Constant outward B.n on one exterior face [T]. Missing flux faces are zero.
#[derive(Clone, Copy, Debug)]
pub struct MagneticFluxFace {
    pub face: usize,
    pub outward_flux_density_t: f64,
}
#[derive(Clone, Debug, Default)]
pub struct SpatialMagneticBoundary {
    pub potential: Vec<MagneticPotentialFace>,
    pub flux: Vec<MagneticFluxFace>,
}

/// Shared Krylov work limits plus dimensional physical acceptance allowances.
#[derive(Clone, Copy, Debug)]
pub struct SpatialMagneticControls {
    pub cg: CheckedCgConfig,
    /// Total-flux compatibility and maximum free/gauge weak residual [Wb].
    pub flux_tolerance_wb: f64,
    /// Independent recovered-energy versus Hodge-energy allowance [J].
    pub energy_tolerance_j: f64,
}

struct Face {
    cells: Vec<usize>,
    area_m2: f64,
    outward_normal: [f64; 3], // outward from cells[0]
}
struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((out, &value), &scale) in z.iter_mut().zip(r).zip(&self.0) {
            *out = value * scale;
        }
    }
}

/// Prepared geometry, source-bound permeability/recoil and sparse scalar Hodge.
/// Matrix preparation is reused across boundary excitations at fixed materials.
pub struct SpatialMagnetostatics {
    mesh: TetComplex,
    positions_m: Vec<[f64; 3]>,
    geometry: ElementGeometry,
    materials: Vec<SpatialMagneticMaterial>,
    permeability_h_m: Vec<f64>,
    remanence_t: Vec<[f64; 3]>,
    faces: Vec<Face>,
    components: Vec<usize>,
    stiffness: Csr,
    remanence_load_wb: Vec<f64>,
}

fn root(parent: &mut [usize], mut v: usize) -> usize {
    while parent[v] != v {
        parent[v] = parent[parent[v]];
        v = parent[v];
    }
    v
}

impl SpatialMagnetostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<SpatialMagneticMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SpatialMagneticError> {
        poll(cancelled)?;
        if tetrahedra.is_empty()
            || tetrahedra.len() != materials.len()
            || tetrahedra.len() > limits.max_cells
            || positions_m.len() > limits.max_dofs
            || tetrahedra
                .len()
                .checked_mul(36)
                .is_none_or(|n| n > limits.max_triplets)
        {
            return Err(invalid(
                "nonempty matching tetrahedra/materials within assembly limits required",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut used = vec![false; positions_m.len()];
        for tet in &tetrahedra {
            poll(cancelled)?;
            let mut sorted = *tet;
            sorted.sort_unstable();
            if sorted.windows(2).any(|p| p[0] == p[1])
                || sorted.iter().any(|&v| v as usize >= positions_m.len())
                || !seen.insert(sorted)
            {
                return Err(invalid("invalid, repeated or duplicate tetrahedron"));
            }
            for &v in tet {
                used[v as usize] = true;
            }
        }
        if used.iter().any(|v| !v) {
            return Err(invalid("isolated vertices are not magnetic volume"));
        }
        let mesh = TetComplex::from_tets(positions_m.len(), tetrahedra);
        if mesh.edges.len() > limits.max_dofs {
            return Err(invalid("magnetic edge assembly exceeds dof limit"));
        }
        let geometry =
            try_element_geometry(&mesh, &positions_m).map_err(SpatialMagneticError::Geometry)?;
        poll(cancelled)?;
        let mut parent: Vec<_> = (0..positions_m.len()).collect();
        for &[a, b] in &mesh.edges {
            let ra = root(&mut parent, a as usize);
            let rb = root(&mut parent, b as usize);
            if ra != rb {
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
        let components = (0..parent.len())
            .map(|v| root(&mut parent, v))
            .collect::<Vec<_>>();
        let mut adjacency = vec![Vec::new(); mesh.faces.len()];
        for (cell, tet) in mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            for skip in 0..4 {
                let mut face = [0u32; 3];
                let mut n = 0;
                for (j, &v) in tet.iter().enumerate() {
                    if j != skip {
                        face[n] = v;
                        n += 1;
                    }
                }
                face.sort_unstable();
                let id = mesh
                    .faces
                    .binary_search(&face)
                    .map_err(|_| invalid("missing tetrahedral face"))?;
                adjacency[id].push(cell);
                if adjacency[id].len() > 2 {
                    return Err(invalid("nonmanifold magnetic face"));
                }
            }
        }
        let mut faces = Vec::with_capacity(mesh.faces.len());
        for (vertices, cells) in mesh.faces.iter().zip(adjacency) {
            poll(cancelled)?;
            let [a, b, c] = vertices.map(|v| positions_m[v as usize]);
            let raw = cross(sub(b, a), sub(c, a));
            let length = finite(norm(raw))?;
            if length <= 0.0 {
                return Err(invalid("degenerate magnetic face"));
            }
            let opposite = |cell: usize| {
                mesh.tets[cell]
                    .iter()
                    .copied()
                    .find(|v| !vertices.contains(v))
                    .unwrap()
            };
            let side = finite(dot(raw, sub(positions_m[opposite(cells[0]) as usize], a)))?;
            if side == 0.0 {
                return Err(invalid("unrepresentable tetrahedron side"));
            }
            if cells.len() == 2 {
                let other = finite(dot(raw, sub(positions_m[opposite(cells[1]) as usize], a)))?;
                if other == 0.0 || side.is_sign_positive() == other.is_sign_positive() {
                    return Err(invalid(
                        "neighbor tetrahedra must lie on opposite face sides",
                    ));
                }
            }
            let direction = if side > 0.0 { -1.0 } else { 1.0 };
            faces.push(Face {
                cells,
                area_m2: 0.5 * length,
                outward_normal: raw.map(|v| direction * v / length),
            });
        }
        let mut permeability_h_m = Vec::with_capacity(materials.len());
        let mut remanence_t = Vec::with_capacity(materials.len());
        for material in &materials {
            poll(cancelled)?;
            let (mu, br) = match material {
                SpatialMagneticMaterial::Linear(source) => {
                    (linear_magnetic_permeability(source)?, [0.0; 3])
                }
                SpatialMagneticMaterial::Recoil(source) => {
                    (source.material.permeability_h_m(), source.remanence_t)
                }
            };
            permeability_h_m.push(mu);
            remanence_t.push(br);
        }
        let weights: Vec<_> = permeability_h_m
            .iter()
            .copied()
            .map(CellWeight::Scalar)
            .collect();
        let mass = weighted_mass_matrix(
            &mesh,
            &geometry,
            1,
            &weights,
            MAGNETIC_PERMEABILITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|e| {
            if e == WeightedError::Cancelled {
                SpatialMagneticError::Cancelled
            } else {
                SpatialMagneticError::Assembly(e)
            }
        })?;
        let stiffness = stiffness(&incidence_to_csr(&mesh.d0()), mass.matrix());
        let mut remanence_load_wb = vec![0.0; positions_m.len()];
        for (cell, tet) in mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            let volume = geometry.vol_signed[cell].abs();
            let mut first = 0.0;
            for local in 1..4 {
                let term = finite(volume * dot(geometry.grads[cell][local], remanence_t[cell]))?;
                remanence_load_wb[tet[local] as usize] =
                    finite(remanence_load_wb[tet[local] as usize] + term)?;
                first = finite(first - term)?;
            }
            remanence_load_wb[tet[0] as usize] =
                finite(remanence_load_wb[tet[0] as usize] + first)?;
        }
        poll(cancelled)?;
        Ok(Self {
            mesh,
            positions_m,
            geometry,
            materials,
            permeability_h_m,
            remanence_t,
            faces,
            components,
            stiffness,
            remanence_load_wb,
        })
    }

    pub const fn mesh(&self) -> &TetComplex {
        &self.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }
    pub fn materials(&self) -> &[SpatialMagneticMaterial] {
        &self.materials
    }
    /// Exterior face IDs in the original canonical face table.
    pub fn exterior_faces(&self) -> Vec<usize> {
        self.faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| (f.cells.len() == 1).then_some(i))
            .collect()
    }
    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.faces
            .get(face)
            .filter(|f| f.cells.len() == 1)
            .map(|f| f.outward_normal)
    }

    /// Solve the declared finite current-free problem. Components with only
    /// flux boundaries receive one deterministic zero-potential gauge AFTER
    /// their total imposed flux is checked; a gauge cannot supply monopole flux.
    pub fn solve(
        &self,
        boundary: &SpatialMagneticBoundary,
        controls: SpatialMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMagneticField<'_>, SpatialMagneticError> {
        poll(cancelled)?;
        if !controls.flux_tolerance_wb.is_finite()
            || controls.flux_tolerance_wb <= 0.0
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || !controls.cg.tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.cg.tolerance)
            || controls.cg.tolerance == 0.0
        {
            return Err(invalid(
                "positive finite magnetic residual/energy controls required",
            ));
        }
        let n = self.positions_m.len();
        let mut prescribed = vec![None; n];
        let mut has_potential = vec![false; n];
        let mut used_faces = BTreeSet::new();
        for trace in &boundary.potential {
            poll(cancelled)?;
            if self
                .faces
                .get(trace.face)
                .is_none_or(|f| f.cells.len() != 1)
                || !used_faces.insert(trace.face)
                || trace.potential_a.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "potential trace needs unique exterior faces and finite values",
                ));
            }
            for (&v, &value) in self.mesh.faces[trace.face].iter().zip(&trace.potential_a) {
                let v = v as usize;
                if prescribed[v].is_some_and(|old| old != value) {
                    return Err(invalid("inconsistent potential at shared trace vertex"));
                }
                prescribed[v] = Some(value);
                has_potential[self.components[v]] = true;
            }
        }
        let physical_prescribed: Vec<_> = prescribed.iter().map(Option::is_some).collect();
        let mut rhs = self.remanence_load_wb.clone();
        let mut net_flux = vec![0.0; n];
        for trace in &boundary.flux {
            poll(cancelled)?;
            if self
                .faces
                .get(trace.face)
                .is_none_or(|f| f.cells.len() != 1)
                || !used_faces.insert(trace.face)
                || !trace.outward_flux_density_t.is_finite()
            {
                return Err(invalid(
                    "flux trace needs unique exterior faces disjoint from potential traces",
                ));
            }
            let amount = finite(trace.outward_flux_density_t * self.faces[trace.face].area_m2)?;
            let vertices = self.mesh.faces[trace.face];
            let component = self.components[vertices[0] as usize];
            net_flux[component] = finite(net_flux[component] + amount)?;
            for v in vertices {
                rhs[v as usize] = finite(rhs[v as usize] - amount / 3.0)?;
            }
        }
        let mut gauges = Vec::new();
        let mut references = vec![0.0; n];
        for v in 0..n {
            if self.components[v] != v {
                continue;
            }
            if !has_potential[v] {
                if net_flux[v].abs() > controls.flux_tolerance_wb {
                    return Err(invalid(
                        "net outward magnetic flux must vanish in a pure-flux component",
                    ));
                }
                prescribed[v] = Some(0.0);
                gauges.push(v);
            } else {
                references[v] = (0..n)
                    .find_map(|i| (self.components[i] == v).then(|| prescribed[i]).flatten())
                    .unwrap();
            }
        }
        // Work in component-relative coordinates so a common gauge offset does
        // not create magnetic fields or contaminate the Hodge energy.
        for v in 0..n {
            if let Some(value) = prescribed[v] {
                prescribed[v] = Some(finite(value - references[self.components[v]])?);
            }
        }
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; n];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        let mut loads = Vec::with_capacity(free.len());
        let mut inverse = Vec::with_capacity(free.len());
        for &v in &free {
            poll(cancelled)?;
            let mut load = rhs[v];
            let mut diagonal = 0.0;
            let (cols, vals) = self.stiffness.row(v);
            for (&w, &value) in cols.iter().zip(vals) {
                if w == v {
                    diagonal = value;
                }
                if let Some(p) = prescribed[w] {
                    load = finite((-value).mul_add(p, load))?;
                } else {
                    columns.push(index[w]);
                    values.push(value);
                }
            }
            if !diagonal.is_finite() || diagonal <= 0.0 || !diagonal.recip().is_finite() {
                return Err(invalid(
                    "magnetic Jacobi diagonal is not representably positive",
                ));
            }
            inverse.push(diagonal.recip());
            loads.push(load);
            row_ptr.push(columns.len());
        }
        let mut relative_potential_a: Vec<_> =
            prescribed.iter().map(|v| v.unwrap_or(0.0)).collect();
        let report = if free.is_empty() {
            None
        } else {
            let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
            let result = checked_cg(
                &CsrOp::symmetric(matrix),
                &Diagonal(inverse),
                &loads,
                controls.cg,
                |_| poll(cancelled),
            )
            .map_err(|e| match e {
                CheckedCgError::InvalidInput(s) => invalid(s),
                CheckedCgError::Interrupted(e) => e,
            })?;
            if !result.report.converged_euclidean() {
                return Err(SpatialMagneticError::Nonconvergence(result.report));
            }
            for (&v, &value) in free.iter().zip(&result.x) {
                relative_potential_a[v] = value;
            }
            Some(result.report)
        };
        let mut applied = vec![0.0; n];
        for v in 0..n {
            poll(cancelled)?;
            let (cols, vals) = self.stiffness.row(v);
            applied[v] = finite(
                cols.iter()
                    .zip(vals)
                    .map(|(&j, &a)| a * relative_potential_a[j])
                    .sum(),
            )?;
        }
        let mut maximum_weak_residual_wb = 0.0_f64;
        for v in 0..n {
            if !physical_prescribed[v] {
                maximum_weak_residual_wb =
                    maximum_weak_residual_wb.max((applied[v] - rhs[v]).abs());
            }
        }
        if maximum_weak_residual_wb > controls.flux_tolerance_wb {
            return Err(invalid(
                "magnetic weak flux residual exceeds its physical allowance",
            ));
        }
        let mut h = Vec::with_capacity(self.materials.len());
        let mut b = Vec::with_capacity(self.materials.len());
        let mut cell_energy_j = Vec::with_capacity(self.materials.len());
        let mut recoil_operating_points = Vec::with_capacity(self.materials.len());
        let mut energy_j = 0.0;
        for (cell, tet) in self.mesh.tets.iter().enumerate() {
            poll(cancelled)?;
            let mut field = [0.0; 3];
            for local in 1..4 {
                let delta = finite(
                    relative_potential_a[tet[local] as usize]
                        - relative_potential_a[tet[0] as usize],
                )?;
                for (j, out) in field.iter_mut().enumerate() {
                    *out = finite((-delta).mul_add(self.geometry.grads[cell][local][j], *out))?;
                }
            }
            let flux = std::array::from_fn(|j| {
                self.permeability_h_m[cell].mul_add(field[j], self.remanence_t[cell][j])
            });
            for value in flux {
                finite(value)?;
            }
            let operating = if let SpatialMagneticMaterial::Recoil(law) = &self.materials[cell] {
                let parallel = finite(dot(field, law.direction))?;
                let transverse = norm(std::array::from_fn(|j| {
                    field[j] - parallel * law.direction[j]
                }));
                Some(
                    law.material
                        .admit_operating_point(parallel, transverse)
                        .map_err(|error| SpatialMagneticError::Recoil { cell, error })?,
                )
            } else {
                None
            };
            let energy = finite(
                0.5 * self.permeability_h_m[cell]
                    * dot(field, field)
                    * self.geometry.vol_signed[cell].abs(),
            )?;
            energy_j = finite(energy_j + energy)?;
            cell_energy_j.push(energy);
            h.push(field);
            b.push(flux);
            recoil_operating_points.push(operating);
        }
        let hodge_energy = finite(
            0.5 * relative_potential_a
                .iter()
                .zip(&applied)
                .map(|(a, b)| a * b)
                .sum::<f64>(),
        )?;
        if (energy_j - hodge_energy).abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "recovered 3D magnetic energy disagrees with the original weighted Hodge",
            ));
        }
        let mut interface_flux_jumps_wb = Vec::new();
        for (face_id, face) in self.faces.iter().enumerate() {
            poll(cancelled)?;
            if face.cells.len() == 2 {
                interface_flux_jumps_wb.push((
                    face_id,
                    finite(
                        face.area_m2
                            * dot(sub(b[face.cells[0]], b[face.cells[1]]), face.outward_normal),
                    )?,
                ));
            }
        }
        let nodal_boundary_flux_wb = self
            .remanence_load_wb
            .iter()
            .zip(&applied)
            .map(|(source, value)| finite(source - value))
            .collect::<Result<Vec<_>, _>>()?;
        let potential_a = relative_potential_a
            .iter()
            .enumerate()
            .map(|(v, &p)| finite(p + references[self.components[v]]))
            .collect::<Result<Vec<_>, _>>()?;
        poll(cancelled)?;
        Ok(SpatialMagneticField {
            problem: self,
            potential_a,
            relative_potential_a,
            field_strength_a_m: h,
            flux_density_t: b,
            cell_energy_j,
            energy_j,
            recoil_operating_points,
            nodal_boundary_flux_wb,
            interface_flux_jumps_wb,
            maximum_weak_residual_wb,
            gauge_vertices: gauges,
            report,
        })
    }
}

/// Accepted fixed-state solution. B/H retain original tetrahedron ordering.
pub struct SpatialMagneticField<'m> {
    problem: &'m SpatialMagnetostatics,
    potential_a: Vec<f64>,
    relative_potential_a: Vec<f64>,
    field_strength_a_m: Vec<[f64; 3]>,
    flux_density_t: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    nodal_boundary_flux_wb: Vec<f64>,
    interface_flux_jumps_wb: Vec<(usize, f64)>,
    maximum_weak_residual_wb: f64,
    gauge_vertices: Vec<usize>,
    report: Option<SolveReport>,
}
impl<'m> SpatialMagneticField<'m> {
    pub const fn problem(&self) -> &'m SpatialMagnetostatics {
        self.problem
    }
    pub fn potential_a(&self) -> &[f64] {
        &self.potential_a
    }
    /// Component-relative potential actually used for stable field recovery.
    pub fn relative_potential_a(&self) -> &[f64] {
        &self.relative_potential_a
    }
    pub fn field_strength_a_m(&self) -> &[[f64; 3]] {
        &self.field_strength_a_m
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.energy_j
    }
    pub fn recoil_operating_points(&self) -> &[Option<RecoilOperatingPoint>] {
        &self.recoil_operating_points
    }
    /// Weak outward boundary flux paired with nodal test functions [Wb].
    pub fn nodal_boundary_flux_wb(&self) -> &[f64] {
        &self.nodal_boundary_flux_wb
    }
    /// Interior face normal-flux jumps of the raw P1 field [Wb]; not equilibrated.
    pub fn interface_flux_jumps_wb(&self) -> &[(usize, f64)] {
        &self.interface_flux_jumps_wb
    }
    pub const fn maximum_weak_residual_wb(&self) -> f64 {
        self.maximum_weak_residual_wb
    }
    pub fn gauge_vertices(&self) -> &[usize] {
        &self.gauge_vertices
    }
    pub const fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if self
            .problem
            .materials
            .iter()
            .any(|m| matches!(m, SpatialMagneticMaterial::Recoil(_)))
        {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
}
