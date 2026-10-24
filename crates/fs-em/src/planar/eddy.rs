//! Stationary planar induction with floating, electrically isolated conductors.
//!
//! The SAME P1 vector potential and material FEEC stiffness obey
//! `K A_mid + C_sigma (A_next-A_old)/dt = load(J_impressed_mid)`.
//! Each connected conductor has an axial scalar-potential gradient chosen so
//! its integrated axial current is zero: `E_z = mean_sigma(dA/dt)-dA/dt`.
//! Thus `C_sigma` is the consistent P1 mass minus its regional weighted-mean
//! projection. It is positive semidefinite; the original anchored magnetic
//! stiffness makes the increment solve positive definite. No lumped or fitted
//! magnetic loss is added. Triangle integrals of sigma*E_z^2 supply actual heat.
//!
//! This is a fixed-geometry, linear-mu magnetoquasistatic section, invariant in z.
//! Opposite axial currents close outside the modeled section; their end fields
//! and end losses are excluded. Fixed A boundary traces and balanced impressed
//! current profiles drive the represented section. There is no terminal circuit,
//! displacement current, Hall effect, hysteresis, moving conductor or force law.
//! Implicit midpoint is second order for held coefficients; thermal endpoint
//! refresh in a composed solve generally reduces coupled accuracy to first order.

use super::{MagneticError, PlanarField, PlanarMagnetostatics};
use fs_material::{
    conductor::ResolvedConductorState,
    field_storage::{FieldStorageError, admit_fixed_magnetic_refresh},
};
use fs_solver::{
    CheckedCgConfig, CheckedCgError, LinearOp, ResidualClaim, SolveReport, checked_cg,
};
use fs_sparse::precond::IdentityPrecond;
use std::{
    cell::{Cell, RefCell},
    collections::BTreeMap,
};

/// Complete conductor source and explicit isolated axial-current region.
#[derive(Clone, Debug, PartialEq)]
pub struct EddyConductor {
    pub triangle: usize,
    pub region: u32,
    pub conductor: ResolvedConductorState,
}

/// Physical acceptance allowances and the existing checked solver's work budget.
#[derive(Clone, Copy, Debug)]
pub struct EddyControls {
    pub linear: CheckedCgConfig,
    pub energy_tolerance_j: f64,
    /// Net axial current of EACH conductor/impressed-current component and
    /// insulating/floating algebraic equilibrium residual allowance [A].
    pub current_balance_tolerance_a: f64,
}

#[derive(Debug)]
pub enum EddyError {
    Magnetic(MagneticError),
    Material(FieldStorageError),
    Invalid(&'static str),
    Balance(&'static str),
    Cancelled,
}
impl From<MagneticError> for EddyError {
    fn from(e: MagneticError) -> Self {
        match e {
            MagneticError::Cancelled => Self::Cancelled,
            other => Self::Magnetic(other),
        }
    }
}
impl From<FieldStorageError> for EddyError {
    fn from(e: FieldStorageError) -> Self {
        Self::Material(e)
    }
}
impl core::fmt::Display for EddyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar induction: {self:?}")
    }
}
impl std::error::Error for EddyError {}

/// Accepted original nodal field, complete conductor sources and physical time.
#[derive(Clone)]
pub struct PlanarEddy<'m> {
    problem: &'m PlanarMagnetostatics,
    conductors: Vec<EddyConductor>,
    insulators: Vec<usize>,
    free: Vec<usize>,
    components: Vec<usize>,
    algebraic: Vec<Vec<usize>>,
    impressed_jz_a_m2: Vec<f64>,
    potential_wb_m: Vec<f64>,
    energy_j: f64,
    time_s: f64,
}

/// Linear nodal E/J on one original triangle; heat integrates their product.
#[derive(Clone, Debug, PartialEq)]
pub struct EddyCellHeat {
    pub source: EddyConductor,
    pub electric_field_v_m: [f64; 3],
    pub current_density_a_m2: [f64; 3],
    pub volume_m3: f64,
    pub power_w: f64,
    pub heat_j: f64,
}
#[derive(Clone, Debug)]
pub struct EddyStep {
    pub cells: Vec<EddyCellHeat>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub heat_j: f64,
    /// Midpoint impressed J_z paired with the actual accepted delta A_z [J].
    pub supplied_work_j: f64,
    pub energy_residual_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
    /// Actual accepted increment's recomputed free-coordinate residual.
    pub report: Option<SolveReport>,
}
pub struct PlanarEddyCandidate<'m> {
    state: PlanarEddy<'m>,
    step: EddyStep,
}
impl<'m> PlanarEddyCandidate<'m> {
    pub fn state(&self) -> &PlanarEddy<'m> {
        &self.state
    }
    pub fn step(&self) -> &EddyStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarEddy<'m>, EddyStep) {
        (self.state, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EddyError> {
    if cancelled() {
        Err(EddyError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64, EddyError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(EddyError::Invalid("unrepresentable induction arithmetic"))
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn controls(c: EddyControls) -> Result<(), EddyError> {
    if !positive(c.energy_tolerance_j)
        || !positive(c.current_balance_tolerance_a)
        || !positive(c.linear.tolerance)
        || c.linear.tolerance >= 1.0
    {
        return Err(EddyError::Invalid(
            "positive finite energy/current allowances and solver tolerance in (0,1) required",
        ));
    }
    Ok(())
}

struct ConductionCell {
    vertices: [usize; 3],
    group: usize,
    sigma_volume: f64,
}
struct Conduction {
    cells: Vec<ConductionCell>,
    regions: Vec<u32>,
    weights: Vec<f64>,
    anchors: Vec<usize>,
}
impl Conduction {
    // Subtract a regional coordinate BEFORE summation; constant gauge values
    // cannot manufacture axial field through cancellation of large offsets.
    fn means(&self, x: &[f64], check: &dyn Fn() -> bool) -> Option<Vec<f64>> {
        let mut means = vec![0.0; self.regions.len()];
        for cell in &self.cells {
            if check() {
                return None;
            }
            let anchor = x[self.anchors[cell.group]];
            means[cell.group] += cell.sigma_volume
                * cell
                    .vertices
                    .iter()
                    .map(|&v| (x[v] - anchor) / 3.0)
                    .sum::<f64>();
        }
        for (v, &weight) in means.iter_mut().zip(&self.weights) {
            *v /= weight;
        }
        Some(means)
    }
    fn apply(&self, x: &[f64], y: &mut [f64], check: &dyn Fn() -> bool) -> bool {
        y.fill(0.0);
        let Some(means) = self.means(x, check) else {
            return false;
        };
        for cell in &self.cells {
            if check() {
                return false;
            }
            let anchor = x[self.anchors[cell.group]];
            let d = cell.vertices.map(|v| (x[v] - anchor) - means[cell.group]);
            let sum = d.iter().sum::<f64>();
            for i in 0..3 {
                y[cell.vertices[i]] += cell.sigma_volume / 12.0 * (d[i] + sum);
            }
        }
        true
    }
}

// An operator adapter only; the existing checked CG owns every recurrence.
struct Increment<'a> {
    problem: &'a PlanarMagnetostatics,
    free: &'a [usize],
    conduction: &'a Conduction,
    dt: f64,
    interrupted: &'a dyn Fn() -> bool,
}
impl LinearOp for Increment<'_> {
    fn n(&self) -> usize {
        self.free.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut full = vec![0.0; self.problem.mesh.vertices().len()];
        for (&v, &value) in self.free.iter().zip(x) {
            full[v] = value;
        }
        let mut c = vec![0.0; full.len()];
        if !self.conduction.apply(&full, &mut c, self.interrupted) {
            y.fill(f64::NAN);
            return;
        }
        for (row, &v) in self.free.iter().enumerate() {
            if (self.interrupted)() {
                y.fill(f64::NAN);
                return;
            }
            let (cols, vals) = self.problem.stiffness.row(v);
            let k = cols
                .iter()
                .zip(vals)
                .map(|(&w, &value)| value * full[w])
                .sum::<f64>();
            y[row] = 0.5 * k + c[v] / self.dt;
        }
    }
}

fn bind(
    problem: &PlanarMagnetostatics,
    mut sources: Vec<EddyConductor>,
    previous: Option<&[EddyConductor]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<EddyConductor>, Conduction), EddyError> {
    if sources.len() > problem.mesh.faces().len() {
        return Err(EddyError::Invalid("too many conductor triangles"));
    }
    sources.sort_by_key(|s| s.triangle);
    if previous.is_some_and(|p| p.len() != sources.len()) {
        return Err(EddyError::Invalid("conductor correspondence changed"));
    }
    let region_ids: BTreeMap<_, _> = sources.iter().map(|s| (s.region, ())).collect();
    let groups: BTreeMap<_, _> = region_ids
        .keys()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect();
    let mut conduction = Conduction {
        cells: Vec::with_capacity(sources.len()),
        regions: groups.keys().copied().collect(),
        weights: vec![0.0; groups.len()],
        anchors: vec![usize::MAX; groups.len()],
    };
    for (i, source) in sources.iter().enumerate() {
        poll(cancelled)?;
        if source.triangle >= problem.mesh.faces().len()
            || i > 0 && sources[i - 1].triangle == source.triangle
        {
            return Err(EddyError::Invalid(
                "duplicate or foreign conductor triangle",
            ));
        }
        if previous.is_some_and(|p| {
            p[i].triangle != source.triangle
                || p[i].region != source.region
                || p[i].conductor.source_basis() != source.conductor.source_basis()
        }) {
            return Err(EddyError::Invalid(
                "conductor region, triangle or independent source basis changed",
            ));
        }
        admit_fixed_magnetic_refresh(
            &problem.materials[source.triangle],
            source.conductor.material(),
        )?;
        if !source
            .conductor
            .source_property()
            .is_static_at_fixed_external_conditions()
        {
            return Err(EddyError::Invalid(
                "eddy currents require a static isotropic Ohmic scalar at fixed external conditions",
            ));
        }
        let sigma = source.conductor.resistivity_ohm_m().recip();
        let volume = problem
            .mesh
            .face_measure(source.triangle)
            .ok_or(EddyError::Invalid("missing conductor volume"))?;
        let weight = sigma * volume;
        if !positive(sigma) || !positive(weight) || !positive(weight / 12.0) {
            return Err(EddyError::Invalid(
                "conductivity pairing must be representably positive",
            ));
        }
        let group = groups[&source.region];
        let vertices = problem.mesh.faces()[source.triangle].map(|v| v as usize);
        conduction.weights[group] = finite(conduction.weights[group] + weight)?;
        if conduction.anchors[group] == usize::MAX {
            conduction.anchors[group] = vertices[0];
        }
        conduction.cells.push(ConductionCell {
            vertices,
            group,
            sigma_volume: weight,
        });
    }
    Ok((sources, conduction))
}

// Electrical contact across a whole edge cannot be silently insulated by a
// region label. A point touch does not create positive-area electrical contact.
fn regions(
    problem: &PlanarMagnetostatics,
    conductors: &[EddyConductor],
    insulators: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, EddyError> {
    let count = problem.mesh.faces().len();
    if conductors.len() + insulators.len() != count {
        return Err(EddyError::Invalid(
            "every triangle must explicitly be conducting or insulating",
        ));
    }
    let mut roles = vec![None; count];
    for (i, s) in conductors.iter().enumerate() {
        roles[s.triangle] = Some(Some(i));
    }
    for &cell in insulators {
        poll(cancelled)?;
        if cell >= count || roles[cell].replace(None).is_some() {
            return Err(EddyError::Invalid(
                "duplicate or foreign insulating triangle",
            ));
        }
    }
    if roles.iter().any(Option::is_none) {
        return Err(EddyError::Invalid("missing electrical triangle role"));
    }
    let mut edges = BTreeMap::<(u32, u32), Vec<usize>>::new();
    let mut neighbors = vec![Vec::new(); conductors.len()];
    let mut parent: Vec<_> = (0..problem.mesh.vertices().len()).collect();
    fn root(parent: &[usize], mut i: usize) -> usize {
        while parent[i] != i {
            i = parent[i];
        }
        i
    }
    for (cell, face) in problem.mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        for k in 0..3 {
            let a = face[k];
            let b = face[(k + 1) % 3];
            edges.entry((a.min(b), a.max(b))).or_default().push(cell);
            let a = root(&parent, a as usize);
            let b = root(&parent, b as usize);
            parent[a.max(b)] = a.min(b);
        }
    }
    for adjacent in edges.values() {
        poll(cancelled)?;
        if adjacent.len() > 2 {
            return Err(EddyError::Invalid("nonmanifold conductor interface"));
        }
        if let [a, b] = adjacent.as_slice() {
            if let (Some(Some(a)), Some(Some(b))) = (roles[*a], roles[*b]) {
                if conductors[a].region != conductors[b].region {
                    return Err(EddyError::Invalid(
                        "edge-connected conductors must share their floating current region",
                    ));
                }
                neighbors[a].push(b);
                neighbors[b].push(a);
            }
        }
    }
    let mut visited = vec![false; conductors.len()];
    let mut used_regions = BTreeMap::new();
    for start in 0..conductors.len() {
        if visited[start] {
            continue;
        }
        if used_regions.insert(conductors[start].region, ()).is_some() {
            return Err(EddyError::Invalid(
                "one floating region cannot join disconnected conductor sections",
            ));
        }
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(i) = stack.pop() {
            poll(cancelled)?;
            for &j in &neighbors[i] {
                if !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    Ok(problem
        .mesh
        .faces()
        .iter()
        .map(|f| root(&parent, f[0] as usize))
        .collect())
}

// Null(C_sigma) consists of air-only nodal hats and a constant on every
// conducting nodal component that does not touch the fixed boundary. These
// are algebraic coordinates, not magnetic time-history degrees of freedom.
fn algebraic_blocks(
    problem: &PlanarMagnetostatics,
    conductors: &[EddyConductor],
    prescribed: &[Option<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Vec<usize>>, EddyError> {
    let n = prescribed.len();
    let mut parent: Vec<_> = (0..n).collect();
    let mut active = vec![false; n];
    fn root(parent: &[usize], mut i: usize) -> usize {
        while parent[i] != i {
            i = parent[i];
        }
        i
    }
    for c in conductors {
        poll(cancelled)?;
        let face = problem.mesh.faces()[c.triangle].map(|v| v as usize);
        for &v in &face {
            active[v] = true;
        }
        for &v in &face[1..] {
            let a = root(&parent, face[0]);
            let b = root(&parent, v);
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    let mut blocks = Vec::new();
    for v in 0..n {
        poll(cancelled)?;
        if active[v] {
            groups.entry(root(&parent, v)).or_default().push(v);
        } else if prescribed[v].is_none() {
            blocks.push(vec![v]);
        }
    }
    for group in groups.into_values() {
        if group.iter().all(|&v| prescribed[v].is_none()) {
            blocks.push(group);
        }
    }
    Ok(blocks)
}

fn load(
    problem: &PlanarMagnetostatics,
    current: &[f64],
    components: &[usize],
    tolerance_a: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, EddyError> {
    if current.len() != components.len() || current.iter().any(|v| !v.is_finite()) {
        return Err(EddyError::Invalid(
            "complete finite impressed current required",
        ));
    }
    let mut result = vec![0.0; problem.mesh.vertices().len()];
    let mut balances = BTreeMap::<usize, f64>::new();
    for (cell, face) in problem.mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        let volume = problem
            .mesh
            .face_measure(cell)
            .ok_or(EddyError::Invalid("missing mesh volume"))?;
        let sum = balances.entry(components[cell]).or_default();
        *sum = finite(*sum + current[cell] * (volume / problem.thickness_m))?;
        for &v in face {
            result[v as usize] = finite(result[v as usize] + current[cell] * volume / 3.0)?;
        }
    }
    if balances.values().any(|v| v.abs() > tolerance_a) {
        return Err(EddyError::Invalid(
            "impressed axial current requires an explicit balanced return on every mesh component",
        ));
    }
    Ok(result)
}

fn magnetic_load(
    problem: &PlanarMagnetostatics,
    potential: &[f64],
    components: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, EddyError> {
    let mut centered = potential.to_vec();
    for (cell, face) in problem.mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        for &v in face {
            centered[v as usize] = finite(potential[v as usize] - potential[components[cell]])?;
        }
    }
    let mut result = Vec::with_capacity(potential.len());
    for v in 0..potential.len() {
        poll(cancelled)?;
        let (cols, vals) = problem.stiffness.row(v);
        result.push(finite(
            cols.iter().zip(vals).map(|(&w, &a)| a * centered[w]).sum(),
        )?);
    }
    Ok(result)
}

fn admit_algebraic(
    problem: &PlanarMagnetostatics,
    magnetic: &[f64],
    load: &[f64],
    blocks: &[Vec<usize>],
    tolerance_a: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), EddyError> {
    for block in blocks {
        let mut residual = 0.0;
        for &v in block {
            poll(cancelled)?;
            residual = finite(residual + (magnetic[v] - load[v]) / problem.thickness_m)?;
        }
        if residual.abs() > tolerance_a {
            return Err(EddyError::Invalid(
                "magnetic field and impressed source violate an algebraic air/floating-region constraint",
            ));
        }
    }
    Ok(())
}

impl<'m> PlanarEddy<'m> {
    /// Begin from an actual field of the SAME prepared owner. No state is reset.
    /// All triangles have an explicit electrical role. Initial conductivity may
    /// use resistivity or conductivity source data; later refresh retains it.
    /// Boundary values remain fixed for this history and must match the field.
    /// Initial impressed current is explicitly zero. Nonzero initial forcing
    /// uses `new_with_impressed_current`; inconsistent algebraic air fields refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &PlanarField<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        time_s: f64,
        limits: EddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EddyError> {
        Self::new_with_impressed_current(
            initial,
            conductors,
            insulating_triangles,
            boundary,
            &vec![0.0; initial.problem.mesh.faces().len()],
            time_s,
            limits,
            cancelled,
        )
    }

    /// Begin with an explicit accepted impressed-current field. This history
    /// makes changing input and the algebraic insulating field consistent at
    /// BOTH endpoints. No initialization impulse or hidden source work is used.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_impressed_current(
        initial: &PlanarField<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        impressed_jz_a_m2: &[f64],
        time_s: f64,
        limits: EddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EddyError> {
        poll(cancelled)?;
        controls(limits)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(EddyError::Invalid(
                "finite nonnegative initial time required",
            ));
        }
        let problem = initial.problem;
        let prescribed = problem.admit_boundary(boundary, cancelled)?;
        for (i, value) in prescribed.iter().enumerate() {
            if value.is_some_and(|a| a.to_bits() != initial.potential_wb_m[i].to_bits()) {
                return Err(EddyError::Invalid(
                    "fixed induction boundary does not match the original field",
                ));
            }
        }
        let (conductors, _) = bind(problem, conductors, None, cancelled)?;
        let components = regions(problem, &conductors, insulating_triangles, cancelled)?;
        let algebraic = algebraic_blocks(problem, &conductors, &prescribed, cancelled)?;
        let initial_load = load(
            problem,
            impressed_jz_a_m2,
            &components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let magnetic = magnetic_load(problem, &initial.potential_wb_m, &components, cancelled)?;
        admit_algebraic(
            problem,
            &magnetic,
            &initial_load,
            &algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            problem,
            conductors,
            insulators: insulating_triangles.to_vec(),
            free: (0..prescribed.len())
                .filter(|&i| prescribed[i].is_none())
                .collect(),
            components,
            algebraic,
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            potential_wb_m: initial.potential_wb_m.clone(),
            energy_j: initial.energy_j,
            time_s,
        })
    }
    pub fn problem(&self) -> &'m PlanarMagnetostatics {
        self.problem
    }
    pub fn conductors(&self) -> &[EddyConductor] {
        &self.conductors
    }
    pub fn insulating_triangles(&self) -> &[usize] {
        &self.insulators
    }
    pub fn potential_wb_m(&self) -> &[f64] {
        &self.potential_wb_m
    }
    pub fn impressed_current_density_a_m2(&self) -> &[f64] {
        &self.impressed_jz_a_m2
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Recover B/H and reuse the existing explicit position/frame transfer.
    pub fn field(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'m>, EddyError> {
        Ok(self
            .problem
            .recover(self.potential_wb_m.clone(), None, cancelled)?)
    }

    /// One candidate interval from unchanged accepted A/history. Impressed
    /// current is the endpoint value, varying linearly from the accepted source,
    /// with a return on each mesh component. The actual midpoint source does the
    /// accounted work. Retained source history prevents spurious alternating
    /// algebraic fields in insulating air when an excitation changes.
    /// Spatial Ohmic heat already accounts for all represented eddy loss; never
    /// add a core-loss map for it. Publish this candidate with its heat receiver.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        conductors: Vec<EddyConductor>,
        impressed_jz_a_m2: &[f64],
        limits: EddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarEddyCandidate<'m>, EddyError> {
        poll(cancelled)?;
        controls(limits)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite()
            || !positive(dt)
            || impressed_jz_a_m2.len() != self.components.len()
            || impressed_jz_a_m2.iter().any(|v| !v.is_finite())
        {
            return Err(EddyError::Invalid(
                "positive representable interval and complete finite impressed current required",
            ));
        }
        let (conductors, conduction) =
            bind(self.problem, conductors, Some(&self.conductors), cancelled)?;
        let endpoint_load = load(
            self.problem,
            impressed_jz_a_m2,
            &self.components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let midpoint_current: Vec<_> = self
            .impressed_jz_a_m2
            .iter()
            .zip(impressed_jz_a_m2)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let loads = load(
            self.problem,
            &midpoint_current,
            &self.components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let old_magnetic = magnetic_load(
            self.problem,
            &self.potential_wb_m,
            &self.components,
            cancelled,
        )?;
        let mut rhs = Vec::with_capacity(self.free.len());
        for &v in &self.free {
            poll(cancelled)?;
            rhs.push(finite(loads[v] - old_magnetic[v])?);
        }
        let mut potential = self.potential_wb_m.clone();
        let report = if self.free.is_empty() {
            None
        } else {
            let cancel = RefCell::new(&mut *cancelled);
            let interrupted = Cell::new(false);
            let check = || {
                if interrupted.get() || (*cancel.borrow_mut())() {
                    interrupted.set(true);
                    true
                } else {
                    false
                }
            };
            let op = Increment {
                problem: self.problem,
                free: &self.free,
                conduction: &conduction,
                dt,
                interrupted: &check,
            };
            let result = checked_cg(&op, &IdentityPrecond, &rhs, limits.linear, |_| {
                if check() {
                    Err(EddyError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .map_err(|e| match e {
                CheckedCgError::InvalidInput(s) => EddyError::Invalid(s),
                CheckedCgError::Interrupted(e) => e,
            })?;
            if !result.report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(result.report).into());
            }
            let mut actual = Vec::with_capacity(self.free.len());
            for (&v, &delta) in self.free.iter().zip(&result.x) {
                potential[v] = finite(potential[v] + delta)?;
                let change = finite(potential[v] - self.potential_wb_m[v])?;
                if delta != 0.0 && change == 0.0 {
                    return Err(EddyError::Invalid(
                        "field increment disappeared in accepted potential",
                    ));
                }
                actual.push(change);
            }
            let mut applied = vec![0.0; rhs.len()];
            op.apply(&actual, &mut applied);
            if check() {
                return Err(EddyError::Cancelled);
            }
            let scale = rhs
                .iter()
                .chain(&applied)
                .map(|v| v.abs())
                .fold(0.0, f64::max);
            let relative = if scale == 0.0 {
                0.0
            } else {
                let denominator = rhs.iter().map(|v| (v / scale).powi(2)).sum::<f64>().sqrt();
                let numerator = rhs
                    .iter()
                    .zip(&applied)
                    .map(|(b, a)| ((b - a) / scale).powi(2))
                    .sum::<f64>()
                    .sqrt();
                if denominator == 0.0 {
                    if numerator == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    numerator / denominator
                }
            };
            let mut history = result.report.history;
            history.push(relative);
            let report = SolveReport::from_claim(
                result.report.iters,
                ResidualClaim::TrueEuclidean(relative),
                limits.linear.tolerance,
                history,
            );
            if !report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(report).into());
            }
            Some(report)
        };
        let new_magnetic = magnetic_load(self.problem, &potential, &self.components, cancelled)?;
        admit_algebraic(
            self.problem,
            &new_magnetic,
            &endpoint_load,
            &self.algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let field = self.problem.recover(potential.clone(), None, cancelled)?;
        let increments: Vec<_> = potential
            .iter()
            .zip(&self.potential_wb_m)
            .map(|(a, b)| a - b)
            .collect();
        let means = {
            let cancel = RefCell::new(&mut *cancelled);
            conduction
                .means(&increments, &|| (*cancel.borrow_mut())())
                .ok_or(EddyError::Cancelled)?
        };
        let mut cells = Vec::with_capacity(conductors.len());
        let mut net_current = vec![0.0; conduction.regions.len()];
        let mut heat_j = 0.0;
        for (source, cell) in conductors.iter().zip(&conduction.cells) {
            poll(cancelled)?;
            let anchor = increments[conduction.anchors[cell.group]];
            let electric_field_v_m = cell
                .vertices
                .map(|v| -((increments[v] - anchor) - means[cell.group]) / dt);
            let sigma = source.conductor.resistivity_ohm_m().recip();
            let current_density_a_m2 = electric_field_v_m.map(|e| sigma * e);
            let volume_m3 = self
                .problem
                .mesh
                .face_measure(source.triangle)
                .ok_or(EddyError::Invalid("missing conductor volume"))?;
            let e_sum = electric_field_v_m.iter().sum::<f64>();
            let power_w = finite(
                cell.sigma_volume / 12.0
                    * (e_sum * e_sum + electric_field_v_m.iter().map(|e| e * e).sum::<f64>()),
            )?;
            let heat = finite(power_w * dt)?;
            if current_density_a_m2.iter().any(|v| !v.is_finite())
                || power_w < 0.0
                || electric_field_v_m != [0.0; 3] && power_w == 0.0
                || power_w != 0.0 && heat == 0.0
            {
                return Err(EddyError::Balance(
                    "represented eddy heat is not finite and nonnegative",
                ));
            }
            net_current[cell.group] = finite(
                net_current[cell.group]
                    + current_density_a_m2.iter().sum::<f64>() / 3.0 * volume_m3
                        / self.problem.thickness_m,
            )?;
            heat_j = finite(heat_j + heat)?;
            cells.push(EddyCellHeat {
                source: source.clone(),
                electric_field_v_m,
                current_density_a_m2,
                volume_m3,
                power_w,
                heat_j: heat,
            });
        }
        if net_current
            .iter()
            .any(|i| i.abs() > limits.current_balance_tolerance_a)
        {
            return Err(EddyError::Balance(
                "floating conductor has unresolved net axial current",
            ));
        }
        let supplied_work_j = finite(loads.iter().zip(&increments).map(|(f, d)| f * d).sum())?;
        let energy_residual_j = finite(field.energy_j - self.energy_j + heat_j - supplied_work_j)?;
        if energy_residual_j.abs() > limits.energy_tolerance_j {
            return Err(EddyError::Balance(
                "magnetic storage, spatial eddy heat and impressed-source work do not close",
            ));
        }
        let step = EddyStep {
            cells,
            stored_before_j: self.energy_j,
            stored_after_j: field.energy_j,
            heat_j,
            supplied_work_j,
            energy_residual_j,
            region_currents_a: conduction.regions.into_iter().zip(net_current).collect(),
            report,
        };
        let state = Self {
            problem: self.problem,
            conductors,
            insulators: self.insulators.clone(),
            free: self.free.clone(),
            components: self.components.clone(),
            algebraic: self.algebraic.clone(),
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            potential_wb_m: potential,
            energy_j: field.energy_j,
            time_s: target_s,
        };
        poll(cancelled)?;
        Ok(PlanarEddyCandidate { state, step })
    }
}
