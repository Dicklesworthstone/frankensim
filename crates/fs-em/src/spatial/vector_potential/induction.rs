//! Fixed-domain 3D moving-conductor magnetoquasistatics on the native Whitney owner.
//!
//! Midpoint Ampere and charge equations use the SAME current:
//! J = sigma*(-Adot-grad(phi)+v cross B_mid).
//! A deterministic relative tree gauge removes only gradient freedom. Harmonic
//! circulation remains a physical unknown: A is never projected independently
//! of the electrical potential. Conducting components have independent floating
//! potentials and insulating natural current boundaries, including internal
//! conductor/insulator interfaces. Tangential A traces stay fixed.
//!
//! Linear/recoil material, conductivity and Eulerian occupancy remain fixed.
//! This models prescribed nonrelativistic motion preserving those coefficient
//! fields, not moving interfaces/ALE, terminals, saturation, hysteresis or heat
//! feedback. Joule and Lorentz transfers require their receivers' own accounts.
//! Source pairing is current-controlled, with piecewise linear impressed current
//! and held velocity. Displacement current is outside the MQS approximation.
//! A-v formulation reference: https://getdp.info/doc/texinfo/getdp.html .

use super::super::{
    SpatialMagneticError, SpatialMagneticMaterial, cross, dot, finite, invalid, poll, sub,
};
use super::{SpatialVectorField, SpatialVectorMagnetostatics, apply};
use fs_material::conductor::ResolvedConductorState;
use fs_solver::{
    FgmresState, FlexiblePreconditioner, LinearOp, ResidualClaim, SolveReport, SolverRunProgress,
};
use fs_sparse::{Csr, ops};
use std::{
    cell::{Cell, RefCell},
    collections::{BTreeMap, VecDeque},
    sync::Arc,
};

type Result<T> = core::result::Result<T, SpatialMagneticError>;

#[derive(Clone, Copy, Debug)]
pub struct SpatialInductionLimits {
    pub maximum_unknowns: usize,
    pub maximum_conducting_cells: usize,
    /// Bounds the original K entries plus 100 local coupled mass entries/cell.
    pub maximum_assembly_terms: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct SpatialInductionControls {
    pub relative_tolerance: f64,
    pub absolute_residual_a: f64,
    pub restart: usize,
    pub maximum_cycles: usize,
    pub maximum_krylov_scalars: usize,
    pub current_tolerance_a: f64,
    pub flux_tolerance_wb: f64,
    pub energy_tolerance_j: f64,
    /// Caller-admitted mu*sigma*max_vertex_speed*diameter/2 for unstabilized FE.
    pub maximum_cell_magnetic_peclet: f64,
}
#[derive(Clone, Debug)]
pub struct SpatialInductionCell {
    pub cell: usize,
    pub component: usize,
    pub volume_m3: f64,
    /// Original local vertex order; phi is floating independently per component.
    pub potential_v: [f64; 4],
    pub laboratory_electric_field_v_m: [[f64; 3]; 4],
    pub conductor_electric_field_v_m: [[f64; 3]; 4],
    pub current_density_a_m2: [[f64; 3]; 4],
    pub lorentz_force_density_n_m3: [[f64; 3]; 4],
    pub heat_j: f64,
    /// Positive is energy delivered by the field to prescribed motion.
    pub mechanical_work_j: f64,
}
pub struct SpatialInductionStep {
    pub cells: Vec<SpatialInductionCell>,
    pub stage_flux_density_t: Vec<[f64; 3]>,
    pub endpoint_flux_density_t: Vec<[f64; 3]>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
    pub impressed_work_j: f64,
    pub energy_residual_j: f64,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub maximum_cell_magnetic_peclet: f64,
    pub report: Option<SolveReport>,
}
pub struct SpatialInductionCandidate<'m> {
    state: SpatialInduction<'m>,
    step: SpatialInductionStep,
}
impl<'m> SpatialInductionCandidate<'m> {
    pub fn state(&self) -> &SpatialInduction<'m> {
        &self.state
    }
    pub fn step(&self) -> &SpatialInductionStep {
        &self.step
    }
    pub fn into_parts(self) -> (SpatialInduction<'m>, SpatialInductionStep) {
        (self.state, self.step)
    }
}
#[derive(Clone)]
pub struct SpatialInduction<'m> {
    prepared: Arc<Prepared<'m>>,
    potential_wb: Vec<f64>,
    impressed_current_a_m2: Vec<[f64; 3]>,
    time_s: f64,
}
struct Prepared<'m> {
    problem: &'m SpatialVectorMagnetostatics,
    sources: Vec<Option<ResolvedConductorState>>,
    cells: Vec<ConductionCell>,
    edge_unknown: Vec<Option<usize>>,
    fixed_edges: Vec<bool>,
    tree_edges: Vec<usize>,
    charge_unknown: Vec<Option<usize>>,
    charge_vertices: Vec<(usize, usize)>,
    mass: Csr,
    magnetic: Csr,
    physical_mass_diagonal: Vec<f64>,
}
struct ConductionCell {
    cell: usize,
    component: usize,
    sigma: f64,
    volume: f64,
    indices: [Option<usize>; 10],
    charges: [usize; 4],
    /// Whitney edges followed by nodal gradients, evaluated at original vertices.
    basis: [[[f64; 3]; 10]; 4],
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
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
/// Exact tetrahedral integral of two P1 vector fields.
fn product(a: [[f64; 3]; 4], b: [[f64; 3]; 4], volume: f64) -> f64 {
    let mut sa = [0.0; 3];
    let mut sb = [0.0; 3];
    let mut diagonal = 0.0;
    for k in 0..4 {
        for j in 0..3 {
            sa[j] += a[k][j];
            sb[j] += b[k][j];
        }
        diagonal += dot(a[k], b[k]);
    }
    volume / 20.0 * (dot(sa, sb) + diagonal)
}
fn row_map(n: usize) -> Vec<BTreeMap<usize, f64>> {
    (0..n).map(|_| BTreeMap::new()).collect()
}
fn add(rows: &mut [BTreeMap<usize, f64>], r: usize, c: usize, value: f64) -> Result<()> {
    let v = rows[r].entry(c).or_default();
    *v = finite(*v + value)?;
    Ok(())
}
fn csr(rows: Vec<BTreeMap<usize, f64>>) -> Csr {
    let n = rows.len();
    let mut ptr = vec![0];
    let mut cols = vec![];
    let mut vals = vec![];
    for row in rows {
        for (c, v) in row {
            if v != 0.0 {
                cols.push(c);
                vals.push(v);
            }
        }
        ptr.push(cols.len());
    }
    Csr::from_parts(n, n, ptr, cols, vals)
}
fn diagonal(matrix: &Csr, r: usize) -> f64 {
    let (c, v) = matrix.row(r);
    c.binary_search(&r).ok().map(|i| v[i]).unwrap_or(0.0)
}
fn field(
    problem: &SpatialVectorMagnetostatics,
    a: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<[f64; 3]>, f64)> {
    let mut b = Vec::with_capacity(problem.cells.len());
    let mut energy = 0.0;
    for (cell, vectors) in problem.cells.iter().enumerate() {
        poll(cancelled)?;
        let mut flux = [0.0; 3];
        for e in 0..6 {
            for j in 0..3 {
                flux[j] = finite(a[vectors.edges[e]].mul_add(vectors.curls[e][j], flux[j]))?;
            }
        }
        let h = std::array::from_fn(|j| {
            (flux[j] - problem.remanence[cell][j]) / problem.permeability[cell]
        });
        for v in h {
            finite(v)?;
        }
        if let SpatialMagneticMaterial::Recoil(m) = &problem.materials[cell] {
            let parallel = finite(dot(h, m.direction()))?;
            let transverse =
                super::super::norm(std::array::from_fn(|j| h[j] - parallel * m.direction()[j]));
            m.material()
                .admit_operating_point(parallel, transverse)
                .map_err(|error| SpatialMagneticError::Recoil { cell, error })?;
        }
        let volume = problem.geometry.geometry.vol_signed[cell].abs();
        energy = finite(energy + 0.5 * problem.permeability[cell] * dot(h, h) * volume)?;
        b.push(flux);
    }
    Ok((b, energy))
}
impl<'m> SpatialInduction<'m> {
    /// Continue an admitted native static field with fixed tangential A traces.
    /// None is an insulator; Some retains the exact original multiphysics bundle.
    /// Conducting components are joined by faces, never by incidental point/edge
    /// contact. They have separate floating scalar potentials and no terminals.
    pub fn new(
        initial: &SpatialVectorField<'m>,
        sources: Vec<Option<ResolvedConductorState>>,
        limits: SpatialInductionLimits,
        flux_tolerance_wb: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let problem = initial.problem();
        let mesh = problem.mesh();
        let nc = sources.iter().filter(|s| s.is_some()).count();
        if sources.len() != mesh.tets.len()
            || nc == 0
            || nc > limits.maximum_conducting_cells
            || !positive(flux_tolerance_wb)
            || nc
                .checked_mul(100)
                .and_then(|n| n.checked_add(problem.stiffness.nnz()))
                .is_none_or(|n| n > limits.maximum_assembly_terms)
        {
            return Err(invalid(
                "complete nonempty bounded conducting/insulating map required",
            ));
        }
        // Relative tree gauge: contract fixed trace edges first, then select a
        // canonical forest among the remaining edge columns. No handle is cut.
        let fixed_edges = initial.prescribed_edges().to_vec();
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
                adjacency[a].push((b, -initial.potential_wb()[e]));
                adjacency[b].push((a, initial.potential_wb()[e]));
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
        let mut potential_wb = initial.potential_wb().to_vec();
        for (e, &[a, b]) in mesh.edges.iter().enumerate() {
            potential_wb[e] = finite(
                potential_wb[e] + shift[fixed_component[b as usize]].unwrap()
                    - shift[fixed_component[a as usize]].unwrap(),
            )?;
        }
        let mut gauged = vec![false; mesh.edges.len()];
        for &e in &tree_edges {
            gauged[e] = true;
            potential_wb[e] = 0.0;
        }
        let old_flux = apply(&problem.curl, initial.potential_wb(), cancelled)?;
        let new_flux = apply(&problem.curl, &potential_wb, cancelled)?;
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
        for face in &problem.geometry.faces {
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
            let Some(source) = source else { continue };
            let magnetic = match &problem.materials[cell] {
                SpatialMagneticMaterial::Linear(m) => m,
                SpatialMagneticMaterial::Recoil(m) => m.material().material(),
            };
            if source.material() != magnetic
                || !source
                    .source_property()
                    .is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "conductivity must retain the exact fixed original magnetic material bundle",
                ));
            }
            let sigma = source.resistivity_ohm_m().recip();
            if !positive(sigma) {
                return Err(invalid("positive representable conductivity required"));
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
        let mut mass = row_map(unknowns);
        let mut physical_mass_diagonal = vec![0.0; mesh.edges.len()];
        for (cell, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            let Some(source) = source else { continue };
            let tet = mesh.tets[cell];
            let grads = &problem.geometry.geometry.grads[cell];
            let mut basis = [[[0.0; 3]; 10]; 4];
            let mut indices = [None; 10];
            let mut e = 0;
            for a in 0..4 {
                for b in a + 1..4 {
                    let (a, b) = if tet[a] < tet[b] { (a, b) } else { (b, a) };
                    indices[e] = edge_unknown[problem.cells[cell].edges[e]];
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
            let volume = problem.geometry.geometry.vol_signed[cell].abs();
            let sigma = source.resistivity_ohm_m().recip();
            for a in 0..10 {
                for b in a..10 {
                    let value =
                        finite(sigma * product(basis.map(|q| q[a]), basis.map(|q| q[b]), volume))?;
                    if a == b && a < 6 {
                        physical_mass_diagonal[problem.cells[cell].edges[a]] += value;
                    }
                    if let (Some(i), Some(j)) = (indices[a], indices[b]) {
                        add(&mut mass, i, j, value)?;
                        if i != j {
                            add(&mut mass, j, i, value)?;
                        }
                    }
                }
            }
            cells.push(ConductionCell {
                cell,
                component: component[cell],
                sigma,
                volume,
                indices,
                charges,
                basis,
            });
        }
        let mut magnetic = row_map(unknowns);
        for (e, index) in edge_unknown.iter().enumerate() {
            poll(cancelled)?;
            if let Some(i) = *index {
                let (cols, vals) = problem.stiffness.row(e);
                for (&f, &value) in cols.iter().zip(vals) {
                    if let Some(j) = edge_unknown[f] {
                        add(&mut magnetic, i, j, value)?;
                    }
                }
            }
        }
        field(problem, &potential_wb, cancelled)?;
        Ok(Self {
            prepared: Arc::new(Prepared {
                problem,
                sources,
                cells,
                edge_unknown,
                fixed_edges,
                tree_edges,
                charge_unknown,
                charge_vertices,
                mass: csr(mass),
                magnetic: csr(magnetic),
                physical_mass_diagonal,
            }),
            potential_wb,
            impressed_current_a_m2: initial.impressed_current_density_a_m2().to_vec(),
            time_s: 0.0,
        })
    }
    pub fn potential_wb(&self) -> &[f64] {
        &self.potential_wb
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn conductors(&self) -> &[Option<ResolvedConductorState>] {
        &self.prepared.sources
    }
    pub fn tree_edges(&self) -> &[usize] {
        &self.prepared.tree_edges
    }
    pub fn problem(&self) -> &'m SpatialVectorMagnetostatics {
        self.prepared.problem
    }

    /// Piecewise linear impressed-current history, held nodal velocity and exact
    /// midpoint magnetic storage. Returned transfers use the same P1 stage current.
    /// No field/clock/source history is published on any admission failure.
    pub fn step_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionCandidate<'m>> {
        poll(cancelled)?;
        let p = &self.prepared;
        let problem = p.problem;
        let mesh = problem.mesh();
        let n = p.mass.nrows();
        if !positive(duration_s)
            || !positive(controls.relative_tolerance)
            || controls.relative_tolerance >= 1.0
            || [
                controls.absolute_residual_a,
                controls.current_tolerance_a,
                controls.flux_tolerance_wb,
                controls.energy_tolerance_j,
                controls.maximum_cell_magnetic_peclet,
            ]
            .into_iter()
            .any(|v| !positive(v))
            || controls.restart == 0
            || controls.maximum_cycles == 0
            || endpoint_current_a_m2.len() != mesh.tets.len()
            || endpoint_current_a_m2
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
            || velocity_m_s.len() != mesh.vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite motion/current and positive bounded induction controls required",
            ));
        }
        let time_s = finite(self.time_s + duration_s)?;
        if time_s <= self.time_s {
            return Err(invalid(
                "induction interval disappears in the retained clock",
            ));
        }
        let restart = controls.restart.min(n.max(1));
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
                controls
                    .maximum_cycles
                    .checked_mul(4)
                    .and_then(|r| v.checked_add(r))
            });
        if memory.is_none_or(|v| v > controls.maximum_krylov_scalars)
            || restart.checked_mul(controls.maximum_cycles).is_none()
        {
            return Err(invalid("induction Krylov work/storage budget exceeded"));
        }
        let endpoint_load = current_load(
            p,
            endpoint_current_a_m2,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let midpoint_current: Vec<_> = self
            .impressed_current_a_m2
            .iter()
            .zip(endpoint_current_a_m2)
            .map(|(a, b)| std::array::from_fn(|j| f64::midpoint(a[j], b[j])))
            .collect();
        let load = current_load(
            p,
            &midpoint_current,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let (b0, stored_before_j) = field(problem, &self.potential_wb, cancelled)?;
        let old_magnetic = apply(&problem.stiffness, &self.potential_wb, cancelled)?;
        let mut rhs = vec![0.0; n];
        for (e, index) in p.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                rhs[i] = finite(load[e] + problem.remanence_load[e] - old_magnetic[e])?;
            }
        }
        let mut motion = row_map(n);
        let mut maximum_cell_magnetic_peclet = 0.0_f64;
        for c in &p.cells {
            poll(cancelled)?;
            let tet = mesh.tets[c.cell];
            let v = tet.map(|i| velocity_m_s[i as usize]);
            let mut diameter = 0.0_f64;
            let mut speed = 0.0_f64;
            for a in 0..4 {
                speed = speed.max(super::super::norm(v[a]));
                for b in a + 1..4 {
                    diameter = diameter.max(super::super::norm(sub(
                        problem.positions_m()[tet[a] as usize],
                        problem.positions_m()[tet[b] as usize],
                    )));
                }
            }
            let peclet = finite(0.5 * problem.permeability[c.cell] * c.sigma * speed * diameter)?;
            maximum_cell_magnetic_peclet = maximum_cell_magnetic_peclet.max(peclet);
            if peclet > controls.maximum_cell_magnetic_peclet {
                return Err(invalid(
                    "moving-conductor magnetic Peclet allowance exceeded",
                ));
            }
            let drive = v.map(|v| cross(v, b0[c.cell]));
            for a in 0..10 {
                let Some(i) = c.indices[a] else { continue };
                rhs[i] =
                    finite(rhs[i] + c.sigma * product(c.basis.map(|q| q[a]), drive, c.volume))?;
                for b in 0..6 {
                    if let Some(j) = c.indices[b] {
                        let motional = v.map(|v| cross(v, problem.cells[c.cell].curls[b]));
                        add(
                            &mut motion,
                            i,
                            j,
                            finite(c.sigma * product(c.basis.map(|q| q[a]), motional, c.volume))?,
                        )?;
                    }
                }
            }
        }
        let motion = csr(motion);
        let transpose = ops::transpose(&motion);
        let inverse: Vec<_> = (0..n)
            .map(|i| {
                let d = finite(0.5 * diagonal(&p.magnetic, i) + diagonal(&p.mass, i) / duration_s)?;
                if !positive(d) || !positive(d.recip()) {
                    Err(invalid(
                        "positive representable induction Jacobi scale required",
                    ))
                } else {
                    Ok(d.recip())
                }
            })
            .collect::<Result<_>>()?;
        let rhs_norm = finite(rhs.iter().fold(0.0_f64, |s, v| s.hypot(*v)))?;
        let target = controls
            .absolute_residual_a
            .max(finite(controls.relative_tolerance * rhs_norm)?);
        let mut x = vec![0.0; n];
        let mut solver_report = None;
        {
            let cancel = RefCell::new(&mut *cancelled);
            let stopped = Cell::new(false);
            let check = || {
                if stopped.get() {
                    true
                } else {
                    let yes = (*cancel.borrow_mut())();
                    stopped.set(yes);
                    yes
                }
            };
            let operator = Increment {
                magnetic: &p.magnetic,
                mass: &p.mass,
                motion: &motion,
                transpose: &transpose,
                dt: duration_s,
                check: &check,
            };
            if rhs_norm > target {
                let mut solver = FgmresState::new(&rhs, restart);
                let tolerance = finite((0.25 * target / rhs_norm).min(0.5))?;
                if !positive(tolerance) {
                    return Err(invalid("unrepresentable induction residual target"));
                }
                let result = solver.run_controlled(
                    &operator,
                    &Jacobi(inverse),
                    &rhs,
                    tolerance,
                    controls.maximum_cycles,
                    &mut || check(),
                );
                if result.progress != SolverRunProgress::Complete || check() {
                    return Err(SpatialMagneticError::Cancelled);
                }
                x = solver.x;
                solver_report = Some(result.report);
            }
            // Re-evaluate the actually retained A increment, not a rounded-away
            // trial. phi has no prior history and remains the solved variable.
            for (e, index) in p.edge_unknown.iter().enumerate() {
                if let Some(i) = *index {
                    let next = finite(self.potential_wb[e] + x[i])?;
                    let actual = finite(next - self.potential_wb[e])?;
                    // Tiny equilibrium corrections may round to zero. The
                    // true residual below tests the retained increment, so any
                    // physically significant loss still refuses the candidate.
                    x[i] = actual;
                }
            }
            let mut ax = vec![0.0; n];
            operator.apply(&x, &mut ax);
            if check() {
                return Err(SpatialMagneticError::Cancelled);
            }
            let residual = finite(
                rhs.iter()
                    .zip(&ax)
                    .fold(0.0_f64, |r, (b, a)| r.hypot(b - a)),
            )?;
            let scale = if rhs_norm == 0.0 { 1.0 } else { rhs_norm };
            let report = SolveReport::from_claim(
                solver_report.as_ref().map_or(0, |r| r.iters),
                ResidualClaim::TrueEuclidean(residual / scale),
                target / scale,
                solver_report.as_ref().map_or_else(
                    || vec![residual / scale],
                    |r| {
                        let mut h = r.history.clone();
                        h.push(residual / scale);
                        h
                    },
                ),
            );
            if !report.converged_euclidean() {
                return Err(SpatialMagneticError::Nonconvergence(report));
            }
            solver_report = Some(report);
        }
        let mut next = self.potential_wb.clone();
        let mut delta = vec![0.0; next.len()];
        for (e, index) in p.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                next[e] = finite(next[e] + x[i])?;
                delta[e] = finite(next[e] - self.potential_wb[e])?;
            }
        }
        let mid: Vec<_> = self
            .potential_wb
            .iter()
            .zip(&next)
            .map(|(a, b)| f64::midpoint(*a, *b))
            .collect();
        let (stage_flux_density_t, _) = field(problem, &mid, cancelled)?;
        let (endpoint_flux_density_t, stored_after_j) = field(problem, &next, cancelled)?;
        let midpoint_magnetic = apply(&problem.stiffness, &mid, cancelled)?;
        let endpoint_magnetic = apply(&problem.stiffness, &next, cancelled)?;
        let mut current_edges = vec![0.0; mesh.edges.len()];
        let mut charge = vec![0.0; p.charge_vertices.len()];
        let mut cells = vec![];
        let mut heat_j = 0.0;
        let mut mechanical_work_j = 0.0;
        for c in &p.cells {
            poll(cancelled)?;
            let vectors = &problem.cells[c.cell];
            let tet = mesh.tets[c.cell];
            let velocity = tet.map(|v| velocity_m_s[v as usize]);
            let potential_v = c
                .charges
                .map(|i| p.charge_unknown[i].map_or(0.0, |j| x[j] / duration_s));
            let mut laboratory = [[0.0; 3]; 4];
            for k in 0..4 {
                for j in 0..3 {
                    let mut electric = 0.0;
                    for a in 0..6 {
                        electric -= delta[vectors.edges[a]] / duration_s * c.basis[k][a][j];
                    }
                    for a in 0..4 {
                        electric -= potential_v[a] * c.basis[k][a + 6][j];
                    }
                    laboratory[k][j] = finite(electric)?;
                }
            }
            let conductor = std::array::from_fn(|k| {
                let motional = cross(velocity[k], stage_flux_density_t[c.cell]);
                std::array::from_fn(|j| laboratory[k][j] + motional[j])
            });
            let current = conductor.map(|e| e.map(|v| v * c.sigma));
            let force = current.map(|j| cross(j, stage_flux_density_t[c.cell]));
            let heat = finite(duration_s * c.sigma * product(conductor, conductor, c.volume))?;
            let work = finite(duration_s * product(velocity, force, c.volume))?;
            if heat < 0.0 {
                return Err(invalid("negative Ohmic dissipation"));
            }
            heat_j = finite(heat_j + heat)?;
            mechanical_work_j = finite(mechanical_work_j + work)?;
            for a in 0..10 {
                let value = finite(product(c.basis.map(|q| q[a]), current, c.volume))?;
                if a < 6 {
                    let e = vectors.edges[a];
                    current_edges[e] = finite(current_edges[e] + value)?;
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
        let maximum_charge_residual_a = charge.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let mut maximum_ampere_residual_a = 0.0_f64;
        for e in 0..mesh.edges.len() {
            poll(cancelled)?;
            if !p.fixed_edges[e] {
                let residual = finite(
                    midpoint_magnetic[e] - problem.remanence_load[e] - load[e] - current_edges[e],
                )?;
                maximum_ampere_residual_a = maximum_ampere_residual_a.max(residual.abs());
                if p.physical_mass_diagonal[e] == 0.0
                    && (endpoint_magnetic[e] - problem.remanence_load[e] - endpoint_load[e]).abs()
                        > controls.current_tolerance_a
                {
                    return Err(invalid("endpoint insulating-domain Ampere equation failed"));
                }
            }
        }
        if maximum_ampere_residual_a > controls.current_tolerance_a
            || maximum_charge_residual_a > controls.current_tolerance_a
        {
            return Err(invalid(
                "full Ampere or floating-component charge equation failed",
            ));
        }
        let face_flux = apply(&problem.curl, &next, cancelled)?;
        for row in mesh.d2().rows {
            let flux = finite(row.iter().map(|&(f, s)| f64::from(s) * face_flux[f]).sum())?;
            if flux.abs() > controls.flux_tolerance_wb {
                return Err(invalid("induction violates discrete div-curl allowance"));
            }
        }
        let impressed_work_j = finite(load.iter().zip(&delta).map(|(j, a)| j * a).sum())?;
        let energy_residual_j = finite(
            stored_after_j - stored_before_j + heat_j + mechanical_work_j - impressed_work_j,
        )?;
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "magnetic, Ohmic, Lorentz and impressed work do not close",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialInductionCandidate {
            state: Self {
                prepared: Arc::clone(&self.prepared),
                potential_wb: next,
                impressed_current_a_m2: endpoint_current_a_m2.to_vec(),
                time_s,
            },
            step: SpatialInductionStep {
                cells,
                stage_flux_density_t,
                endpoint_flux_density_t,
                stored_before_j,
                stored_after_j,
                heat_j,
                mechanical_work_j,
                impressed_work_j,
                energy_residual_j,
                maximum_ampere_residual_a,
                maximum_charge_residual_a,
                maximum_cell_magnetic_peclet,
                report: solver_report,
            },
        })
    }
}
fn current_load(
    p: &Prepared<'_>,
    current: &[[f64; 3]],
    tol: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let problem = p.problem;
    for (i, face) in problem.geometry.faces.iter().enumerate() {
        poll(cancelled)?;
        let jump = if face.cells.len() == 2 {
            sub(current[face.cells[0]], current[face.cells[1]])
        } else if problem
            .face_edges(i)
            .unwrap()
            .iter()
            .all(|&e| p.fixed_edges[e])
        {
            continue;
        } else {
            current[face.cells[0]]
        };
        if finite(face.area_m2 * dot(jump, face.outward_normal))?.abs() > tol {
            return Err(invalid(
                "impressed current violates interface or natural-boundary continuity",
            ));
        }
    }
    let mut load = vec![0.0; problem.mesh().edges.len()];
    for (cell, vectors) in problem.cells.iter().enumerate() {
        poll(cancelled)?;
        for e in 0..6 {
            let edge = vectors.edges[e];
            load[edge] = finite(load[edge] + dot(current[cell], vectors.integrals[e]))?;
        }
    }
    Ok(load)
}
struct Increment<'a> {
    magnetic: &'a Csr,
    mass: &'a Csr,
    motion: &'a Csr,
    transpose: &'a Csr,
    dt: f64,
    check: &'a dyn Fn() -> bool,
}
impl Increment<'_> {
    fn action(&self, x: &[f64], y: &mut [f64], motion: &Csr) {
        for i in 0..y.len() {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let row = |a: &Csr| {
                let (c, v) = a.row(i);
                c.iter().zip(v).map(|(&j, &v)| v * x[j]).sum::<f64>()
            };
            y[i] = 0.5 * row(self.magnetic) + row(self.mass) / self.dt - 0.5 * row(motion);
        }
    }
}
impl LinearOp for Increment<'_> {
    fn n(&self) -> usize {
        self.mass.nrows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, self.motion)
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, self.transpose)
    }
}
struct Jacobi(Vec<f64>);
impl FlexiblePreconditioner for Jacobi {
    fn apply(&self, _: usize, r: &[f64], z: &mut [f64]) {
        for ((z, r), d) in z.iter_mut().zip(r).zip(&self.0) {
            *z = r * d;
        }
    }
}
