//! One relative gauge, electric-charge space and Whitney conductor operator.
//!
//! Linear, nonlinear and thermal induction retain these original geometric
//! identities. Constitutive owners supply magnetic effort and its tangent;
//! this module supplies no replacement magnetic material or time integrator.

use super::super::super::{SpatialGeometry, cross, dot, invalid, norm, poll, sub};
use super::super::{CellVectors, VectorTopology, apply};
use super::{
    Result, SpatialInductionCell, SpatialInductionControls, SpatialInductionLimits, add, csr,
    finite, positive, product, row_map,
};
use fs_material::conductor::spatial::{SpatialConductivity, SpatialConductorMaterial};
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;
use std::collections::{BTreeMap, VecDeque};

#[derive(Clone, Copy)]
pub(super) struct Geometry<'a> {
    pub geometry: &'a SpatialGeometry,
    pub cells: &'a [CellVectors],
    pub curl: &'a Csr,
}
impl Geometry<'_> {
    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn flux(
        &self,
        potential: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>> {
        let mut out = Vec::with_capacity(self.cells.len());
        for cell in self.cells {
            poll(cancelled)?;
            let mut b = [0.0; 3];
            for e in 0..6 {
                for j in 0..3 {
                    b[j] = finite(potential[cell.edges[e]].mul_add(cell.curls[e][j], b[j]))?;
                }
            }
            out.push(b);
        }
        Ok(out)
    }
}

#[derive(Clone)]
pub(super) struct ConductionCell {
    pub cell: usize,
    pub component: usize,
    pub sigma: SpatialConductivity,
    pub volume: f64,
    pub indices: [Option<usize>; 10],
    pub charges: [usize; 4],
    /// Whitney edges followed by nodal gradients, at original vertices.
    pub basis: [[[f64; 3]; 10]; 4],
}
#[derive(Clone)]
pub(super) struct Coupling<C: SpatialConductorMaterial> {
    pub sources: Vec<Option<C>>,
    pub cells: Vec<ConductionCell>,
    pub edge_unknown: Vec<Option<usize>>,
    pub fixed_edges: Vec<bool>,
    pub tree_edges: Vec<usize>,
    pub charge_unknown: Vec<Option<usize>>,
    pub charge_vertices: Vec<(usize, usize)>,
    pub mass: Csr,
    pub physical_mass_diagonal: Vec<f64>,
}
pub(super) struct Current {
    pub cells: Vec<SpatialInductionCell>,
    pub edge_current_a: Vec<f64>,
    pub maximum_charge_residual_a: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
}
fn root(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}
fn join(parent: &mut [usize], a: usize, b: usize) -> bool {
    let (a, b) = (root(parent, a), root(parent, b));
    if a == b {
        false
    } else {
        parent[a.max(b)] = a.min(b);
        true
    }
}
impl<C: SpatialConductorMaterial> Coupling<C> {
    /// Prepare once from a source-admitted field's actual trace and cochain.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        g: Geometry<'_>,
        initial: &[f64],
        fixed_edges: Vec<bool>,
        sources: Vec<Option<C>>,
        magnetic_entry_bound: usize,
        limits: SpatialInductionLimits,
        flux_tolerance_wb: f64,
        mut admit: impl FnMut(usize, &C) -> Result<()>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Vec<f64>)> {
        poll(cancelled)?;
        let mesh = g.mesh();
        let nc = sources.iter().filter(|s| s.is_some()).count();
        if sources.len() != mesh.tets.len()
            || nc == 0
            || nc > limits.maximum_conducting_cells
            || !positive(flux_tolerance_wb)
            || initial.len() != mesh.edges.len()
            || fixed_edges.len() != mesh.edges.len()
            || initial.iter().any(|v| !v.is_finite())
            || nc
                .checked_mul(100)
                .and_then(|n| n.checked_add(magnetic_entry_bound))
                .is_none_or(|n| n > limits.maximum_assembly_terms)
        {
            return Err(invalid(
                "complete bounded conductor, trace and potential maps required",
            ));
        }
        let mut fixed_parent: Vec<_> = (0..mesh.vertex_count).collect();
        for (e, &[a, b]) in mesh.edges.iter().enumerate() {
            poll(cancelled)?;
            if fixed_edges[e] {
                join(&mut fixed_parent, a as usize, b as usize);
            }
        }
        let fixed_component: Vec<_> = (0..mesh.vertex_count)
            .map(|v| root(&mut fixed_parent, v))
            .collect();
        let mut forest = fixed_parent.clone();
        let mut tree_edges = vec![];
        let mut adjacency = vec![vec![]; mesh.vertex_count];
        for (e, &[a, b]) in mesh.edges.iter().enumerate() {
            poll(cancelled)?;
            if !fixed_edges[e] && join(&mut forest, a as usize, b as usize) {
                tree_edges.push(e);
                let (a, b) = (fixed_component[a as usize], fixed_component[b as usize]);
                adjacency[a].push((b, -initial[e]));
                adjacency[b].push((a, initial[e]));
            }
        }
        let mut shift = vec![None; mesh.vertex_count];
        for v in 0..mesh.vertex_count {
            if fixed_component[v] != v || shift[v].is_some() {
                continue;
            }
            shift[v] = Some(0.0);
            let mut q = VecDeque::from([v]);
            while let Some(a) = q.pop_front() {
                poll(cancelled)?;
                for &(b, delta) in &adjacency[a] {
                    if shift[b].is_none() {
                        shift[b] = Some(finite(shift[a].unwrap() + delta)?);
                        q.push_back(b);
                    }
                }
            }
        }
        let mut potential = initial.to_vec();
        for (e, &[a, b]) in mesh.edges.iter().enumerate() {
            potential[e] = finite(
                potential[e] + shift[fixed_component[b as usize]].unwrap()
                    - shift[fixed_component[a as usize]].unwrap(),
            )?;
        }
        let mut gauged = vec![false; mesh.edges.len()];
        for &e in &tree_edges {
            gauged[e] = true;
            potential[e] = 0.0;
        }
        let old_flux = apply(g.curl, initial, cancelled)?;
        let new_flux = apply(g.curl, &potential, cancelled)?;
        if old_flux
            .iter()
            .zip(&new_flux)
            .any(|(a, b)| (a - b).abs() > flux_tolerance_wb)
        {
            return Err(invalid(
                "relative tree gauge changed initial magnetic face flux",
            ));
        }
        let mut edge_unknown = vec![None; mesh.edges.len()];
        let mut unknowns = 0;
        for e in 0..mesh.edges.len() {
            if !fixed_edges[e] && !gauged[e] {
                edge_unknown[e] = Some(unknowns);
                unknowns += 1;
            }
        }
        let mut regions: Vec<_> = (0..mesh.tets.len()).collect();
        for face in &g.geometry.faces {
            poll(cancelled)?;
            if face.cells.len() == 2
                && sources[face.cells[0]].is_some()
                && sources[face.cells[1]].is_some()
            {
                join(&mut regions, face.cells[0], face.cells[1]);
            }
        }
        let mut charge_map = BTreeMap::new();
        let mut component = vec![usize::MAX; mesh.tets.len()];
        for (cell, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            let Some(source) = source else {
                continue;
            };
            admit(cell, source)?;
            if !source.has_admitted_static_sources()
                || !positive(source.conductivity().maximum_conductivity_s_m())
            {
                return Err(invalid("positive fixed-condition conductivity required"));
            }
            component[cell] = root(&mut regions, cell);
            for &v in &mesh.tets[cell] {
                charge_map.insert((component[cell], v as usize), 0);
            }
        }
        let mut charge_unknown = vec![];
        let mut charge_vertices = vec![];
        let mut last_region = None;
        for (&key, index) in &mut charge_map {
            *index = charge_vertices.len();
            charge_vertices.push(key);
            if last_region != Some(key.0) {
                charge_unknown.push(None);
                last_region = Some(key.0);
            } else {
                charge_unknown.push(Some(unknowns));
                unknowns += 1;
            }
        }
        if unknowns > limits.maximum_unknowns {
            return Err(invalid("coupled A-phi unknown budget exceeded"));
        }
        let mut cells = Vec::with_capacity(nc);
        for (cell, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            let Some(source) = source else {
                continue;
            };
            let tet = mesh.tets[cell];
            let grads = &g.geometry.geometry.grads[cell];
            let mut basis = [[[0.0; 3]; 10]; 4];
            let mut indices = [None; 10];
            let mut e = 0;
            for a in 0..4 {
                for b in a + 1..4 {
                    let (a, b) = if tet[a] < tet[b] { (a, b) } else { (b, a) };
                    indices[e] = edge_unknown[g.cells[cell].edges[e]];
                    basis[a][e] = grads[b];
                    basis[b][e] = grads[a].map(|v| -v);
                    e += 1;
                }
            }
            let charges = tet.map(|v| charge_map[&(component[cell], v as usize)]);
            for a in 0..4 {
                indices[a + 6] = charge_unknown[charges[a]];
                for row in &mut basis {
                    row[a + 6] = grads[a];
                }
            }
            cells.push(ConductionCell {
                cell,
                component: component[cell],
                sigma: source.conductivity(),
                volume: g.geometry.geometry.vol_signed[cell].abs(),
                indices,
                charges,
                basis,
            });
        }
        let (mass, physical_mass_diagonal) = assemble_mass(g, &cells, unknowns, cancelled)?;
        Ok((
            Self {
                sources,
                cells,
                edge_unknown,
                fixed_edges,
                tree_edges,
                charge_unknown,
                charge_vertices,
                mass,
                physical_mass_diagonal,
            },
            potential,
        ))
    }

    /// Reweight the SAME electric support and cached original local basis.
    /// The constitutive caller separately admits every retained source identity.
    pub fn refresh_sources(
        &self,
        g: Geometry<'_>,
        sources: Vec<Option<C>>,
        mut admit: impl FnMut(usize, &C) -> Result<()>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if sources.len() != self.sources.len()
            || sources
                .iter()
                .zip(&self.sources)
                .any(|(a, b)| a.is_some() != b.is_some())
        {
            return Err(invalid(
                "thermal induction must retain the original conducting support",
            ));
        }
        let mut cells = self.cells.clone();
        for cell in &mut cells {
            poll(cancelled)?;
            let source = sources[cell.cell]
                .as_ref()
                .expect("unchanged conducting support");
            admit(cell.cell, source)?;
            cell.sigma = source.conductivity();
            if !positive(cell.sigma.maximum_conductivity_s_m())
                || !source.has_admitted_static_sources()
            {
                return Err(invalid("positive fixed-condition conductivity required"));
            }
        }
        let (mass, physical_mass_diagonal) =
            assemble_mass(g, &cells, self.mass.nrows(), cancelled)?;
        poll(cancelled)?;
        Ok(Self {
            sources,
            cells,
            mass,
            physical_mass_diagonal,
            edge_unknown: self.edge_unknown.clone(),
            fixed_edges: self.fixed_edges.clone(),
            tree_edges: self.tree_edges.clone(),
            charge_unknown: self.charge_unknown.clone(),
            charge_vertices: self.charge_vertices.clone(),
        })
    }

    pub fn restrict(&self, matrix: &Csr, cancelled: &mut impl FnMut() -> bool) -> Result<Csr> {
        let mut rows = row_map(self.mass.nrows());
        for (e, index) in self.edge_unknown.iter().enumerate() {
            poll(cancelled)?;
            if let Some(i) = *index {
                let (cols, values) = matrix.row(e);
                for (&f, &value) in cols.iter().zip(values) {
                    if let Some(j) = self.edge_unknown[f] {
                        add(&mut rows, i, j, value)?;
                    }
                }
            }
        }
        Ok(csr(rows))
    }

    pub fn current_load(
        &self,
        g: Geometry<'_>,
        current: &[[f64; 3]],
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let topology = VectorTopology {
            geometry: g.geometry,
            cells: g.cells,
        };
        for (i, face) in g.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            let jump = if face.cells.len() == 2 {
                sub(current[face.cells[0]], current[face.cells[1]])
            } else if topology
                .face_edges(i)
                .unwrap()
                .iter()
                .all(|&e| self.fixed_edges[e])
            {
                continue;
            } else {
                current[face.cells[0]]
            };
            if finite(face.area_m2 * dot(jump, face.outward_normal))?.abs() > tolerance {
                return Err(invalid(
                    "impressed current violates interface or natural-boundary continuity",
                ));
            }
        }
        let mut load = vec![0.0; g.mesh().edges.len()];
        for (cell, vectors) in g.cells.iter().enumerate() {
            poll(cancelled)?;
            for e in 0..6 {
                let edge = vectors.edges[e];
                load[edge] = finite(load[edge] + dot(current[cell], vectors.integrals[e]))?;
            }
        }
        Ok(load)
    }

    /// Full sparse sigma*(v cross curl) and its action on the accepted field.
    pub fn motion(
        &self,
        g: Geometry<'_>,
        velocity: &[[f64; 3]],
        old_b: &[[f64; 3]],
        permeability_bound: &[f64],
        maximum_peclet: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Csr, Vec<f64>, f64)> {
        let n = self.mass.nrows();
        let mut rows = row_map(n);
        let mut drive = vec![0.0; n];
        let mut observed = 0.0_f64;
        for c in &self.cells {
            poll(cancelled)?;
            let tet = g.mesh().tets[c.cell];
            let v = tet.map(|i| velocity[i as usize]);
            let mut diameter = 0.0_f64;
            let mut speed = 0.0_f64;
            for a in 0..4 {
                speed = speed.max(norm(v[a]));
                for b in a + 1..4 {
                    diameter = diameter.max(norm(sub(
                        g.geometry.positions_m[tet[a] as usize],
                        g.geometry.positions_m[tet[b] as usize],
                    )));
                }
            }
            let peclet = finite(0.5 * permeability_bound[c.cell]
                * c.sigma.maximum_conductivity_s_m() * speed * diameter)?;
            observed = observed.max(peclet);
            if peclet > maximum_peclet {
                return Err(invalid(
                    "moving-conductor magnetic Peclet allowance exceeded",
                ));
            }
            let initial_motion = v.map(|v| cross(v, old_b[c.cell]));
            for a in 0..10 {
                let Some(i) = c.indices[a] else {
                    continue;
                };
                drive[i] = finite(
                    drive[i] + conductivity_product(c.sigma, c.basis.map(|q| q[a]), initial_motion, c.volume),
                )?;
                for b in 0..6 {
                    if let Some(j) = c.indices[b] {
                        let motion = v.map(|v| cross(v, g.cells[c.cell].curls[b]));
                        add(
                            &mut rows,
                            i,
                            j,
                            finite(conductivity_product(c.sigma, c.basis.map(|q| q[a]), motion, c.volume))?,
                        )?;
                    }
                }
            }
        }
        Ok((csr(rows), drive, observed))
    }

    /// Recover all physical transfers and all charge equations from ONE current.
    #[allow(clippy::too_many_arguments)]
    pub fn recover(
        &self,
        g: Geometry<'_>,
        delta: &[f64],
        solution: &[f64],
        stage_b: &[[f64; 3]],
        velocity: &[[f64; 3]],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Current> {
        let mut edge_current_a = vec![0.0; g.mesh().edges.len()];
        let mut charge = vec![0.0; self.charge_vertices.len()];
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut heat_j = 0.0;
        let mut mechanical_work_j = 0.0;
        for c in &self.cells {
            poll(cancelled)?;
            let vectors = &g.cells[c.cell];
            let velocity = g.mesh().tets[c.cell].map(|v| velocity[v as usize]);
            let potential_v = c
                .charges
                .map(|i| self.charge_unknown[i].map_or(0.0, |j| solution[j] / dt));
            let mut laboratory = [[0.0; 3]; 4];
            for k in 0..4 {
                for j in 0..3 {
                    let mut electric = 0.0;
                    for a in 0..6 {
                        electric -= delta[vectors.edges[a]] / dt * c.basis[k][a][j];
                    }
                    for a in 0..4 {
                        electric -= potential_v[a] * c.basis[k][a + 6][j];
                    }
                    laboratory[k][j] = finite(electric)?;
                }
            }
            let conductor = std::array::from_fn(|k| {
                let motion = cross(velocity[k], stage_b[c.cell]);
                std::array::from_fn(|j| laboratory[k][j] + motion[j])
            });
            let current = conductor.map(|e| c.sigma.current_density_a_m2(e));
            let force = current.map(|j| cross(j, stage_b[c.cell]));
            let heat = finite(match c.sigma {
                SpatialConductivity::Scalar(s) => dt * s * product(conductor, conductor, c.volume),
                SpatialConductivity::Tensor(_) => dt * product(conductor, current, c.volume),
            })?;
            let work = finite(dt * product(velocity, force, c.volume))?;
            if heat < 0.0 {
                return Err(invalid("negative Ohmic dissipation"));
            }
            heat_j = finite(heat_j + heat)?;
            mechanical_work_j = finite(mechanical_work_j + work)?;
            for a in 0..10 {
                let value = finite(product(c.basis.map(|q| q[a]), current, c.volume))?;
                if a < 6 {
                    let e = vectors.edges[a];
                    edge_current_a[e] = finite(edge_current_a[e] + value)?;
                } else {
                    let i = c.charges[a - 6];
                    charge[i] = finite(charge[i] + value)?;
                }
            }
            cells.push(SpatialInductionCell {
                cell: c.cell,
                component: c.component,
                volume_m3: c.volume,
                potential_v,
                laboratory_electric_field_v_m: laboratory,
                conductor_electric_field_v_m: conductor,
                current_density_a_m2: current,
                lorentz_force_density_n_m3: force,
                heat_j: heat,
                mechanical_work_j: work,
            });
        }
        Ok(Current {
            cells,
            edge_current_a,
            maximum_charge_residual_a: charge.iter().fold(0.0_f64, |m, v| m.max(v.abs())),
            heat_j,
            mechanical_work_j,
        })
    }
}

// Exact coefficient-weighted P1 integral; scalar sources retain their original
// arithmetic while the tensor acts inside the same Whitney/charge bilinear form.
fn conductivity_product(
    sigma: SpatialConductivity,
    a: [[f64; 3]; 4],
    b: [[f64; 3]; 4],
    volume: f64,
) -> f64 {
    match sigma {
        SpatialConductivity::Scalar(s) => s * product(a, b, volume),
        SpatialConductivity::Tensor(_) =>
            product(a, b.map(|v| sigma.current_density_a_m2(v)), volume),
    }
}

fn assemble_mass(
    g: Geometry<'_>,
    cells: &[ConductionCell],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Csr, Vec<f64>)> {
    let mut mass = row_map(n);
    let mut diagonal = vec![0.0; g.mesh().edges.len()];
    for c in cells {
        poll(cancelled)?;
        for a in 0..10 {
            for b in a..10 {
                let value = finite(
                    conductivity_product(c.sigma, c.basis.map(|q| q[a]), c.basis.map(|q| q[b]), c.volume),
                )?;
                if a == b && a < 6 {
                    let e = g.cells[c.cell].edges[a];
                    diagonal[e] = finite(diagonal[e] + value)?;
                }
                if let (Some(i), Some(j)) = (c.indices[a], c.indices[b]) {
                    add(&mut mass, i, j, value)?;
                    if i != j {
                        add(&mut mass, j, i, value)?;
                    }
                }
            }
        }
    }
    Ok((csr(mass), diagonal))
}

pub(super) fn admit_step<C: SpatialConductorMaterial>(
    g: Geometry<'_>,
    coupling: &Coupling<C>,
    old_time: f64,
    current: &[[f64; 3]],
    velocity: &[[f64; 3]],
    dt: f64,
    c: SpatialInductionControls,
) -> Result<(f64, usize)> {
    if !positive(dt)
        || !positive(c.relative_tolerance)
        || c.relative_tolerance >= 1.0
        || [
            c.absolute_residual_a,
            c.current_tolerance_a,
            c.flux_tolerance_wb,
            c.energy_tolerance_j,
            c.maximum_cell_magnetic_peclet,
        ]
        .into_iter()
        .any(|v| !positive(v))
        || c.restart == 0
        || c.maximum_cycles == 0
        || current.len() != g.mesh().tets.len()
        || current.iter().flatten().any(|v| !v.is_finite())
        || velocity.len() != g.mesh().vertex_count
        || velocity.iter().flatten().any(|v| !v.is_finite())
    {
        return Err(invalid(
            "complete finite motion/current and positive bounded induction controls required",
        ));
    }
    let time = finite(old_time + dt)?;
    if time <= old_time {
        return Err(invalid(
            "induction interval disappears in the retained clock",
        ));
    }
    let n = coupling.mass.nrows();
    let restart = c.restart.min(n.max(1));
    let memory = restart
        .checked_mul(2)
        .and_then(|v| v.checked_add(12))
        .and_then(|v| v.checked_mul(n))
        .and_then(|v| {
            restart
                .checked_add(8)
                .and_then(|r| r.checked_mul(restart))
                .and_then(|r| v.checked_add(r))
        })
        .and_then(|v| {
            c.maximum_cycles
                .checked_mul(4)
                .and_then(|r| v.checked_add(r))
        });
    if memory.is_none_or(|v| v > c.maximum_krylov_scalars)
        || restart.checked_mul(c.maximum_cycles).is_none()
    {
        return Err(invalid("induction Krylov work/storage budget exceeded"));
    }
    Ok((time, restart))
}
