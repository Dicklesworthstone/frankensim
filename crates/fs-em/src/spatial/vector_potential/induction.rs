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
//! The base step holds linear/recoil material, conductivity and Eulerian
//! occupancy fixed. The nonlinear child supplies reversible saturation; the
//! thermal child refreshes source-bound linear material and conductivity.
//! Motion is prescribed and nonrelativistic. Moving interfaces/ALE, terminals
//! and hysteresis are outside these owners. Joule and Lorentz transfers require
//! their receivers' own accounts.
//! Source pairing is current-controlled, with piecewise linear impressed current
//! and held velocity. Displacement current is outside the MQS approximation.
//! A-v formulation reference: https://getdp.info/doc/texinfo/getdp.html .

mod coupling;
pub mod nonlinear;
pub mod thermal;

use super::super::{SpatialMagneticError, SpatialMagneticMaterial, dot, finite, invalid, poll};
use super::{SpatialVectorField, SpatialVectorMagnetostatics, apply};
use coupling::{Coupling, Geometry, admit_step};
use fs_material::conductor::ResolvedConductorState;
use fs_solver::{
    FgmresState, FlexiblePreconditioner, LinearOp, ResidualClaim, SolveReport, SolverRunProgress,
};
use fs_sparse::{Csr, ops};
use std::{
    cell::{Cell, RefCell},
    collections::BTreeMap,
    ops::Deref,
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
    /// Original per-cell magnetic energies and source-parameter work.
    pub magnetic_cells: Vec<SpatialInductionMagneticCell>,
    /// Signed magnetic parameter work; not heat or impressed-current supply.
    pub parameter_work_j: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
    pub impressed_work_j: f64,
    pub energy_residual_j: f64,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub maximum_cell_magnetic_peclet: f64,
    pub report: Option<SolveReport>,
}
#[derive(Clone, Debug)]
pub struct SpatialInductionMagneticCell {
    pub cell: usize,
    pub volume_m3: f64,
    pub flux_density_t: [[f64; 3]; 2],
    pub free_energy_j: [f64; 2],
    pub parameter_work_j: f64,
}
/// Original constitutive recovery of an accepted induction cochain.
/// Private fields prevent construction from detached arbitrary B/H samples.
pub struct SpatialInductionField {
    flux_density_t: Vec<[f64; 3]>,
    field_strength_a_m: Vec<[f64; 3]>,
    stored_energy_j: f64,
}
impl SpatialInductionField {
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    pub fn field_strength_a_m(&self) -> &[[f64; 3]] {
        &self.field_strength_a_m
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
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
    problem: MagneticOwner<'m>,
    coupling: Coupling,
    magnetic: Csr,
}
impl Deref for Prepared<'_> {
    type Target = Coupling;
    fn deref(&self) -> &Self::Target {
        &self.coupling
    }
}
#[derive(Clone)]
enum MagneticOwner<'m> {
    Borrowed(&'m SpatialVectorMagnetostatics),
    Owned(Arc<SpatialVectorMagnetostatics>),
}
impl Deref for MagneticOwner<'_> {
    type Target = SpatialVectorMagnetostatics;
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Borrowed(p) => p,
            Self::Owned(p) => p,
        }
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
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
fn geometry(problem: &SpatialVectorMagnetostatics) -> Geometry<'_> {
    Geometry {
        geometry: &problem.geometry,
        cells: &problem.cells,
        curl: &problem.curl,
    }
}
fn field(
    problem: &SpatialVectorMagnetostatics,
    a: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialInductionField> {
    let mut b = Vec::with_capacity(problem.cells.len());
    let mut strengths = Vec::with_capacity(problem.cells.len());
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
        strengths.push(h);
    }
    Ok(SpatialInductionField {
        flux_density_t: b,
        field_strength_a_m: strengths,
        stored_energy_j: energy,
    })
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
        let (coupling, potential_wb) = Coupling::prepare(
            geometry(problem),
            initial.potential_wb(),
            initial.prescribed_edges().to_vec(),
            sources,
            problem.stiffness.nnz(),
            limits,
            flux_tolerance_wb,
            |cell, source| {
                let magnetic = match &problem.materials[cell] {
                    SpatialMagneticMaterial::Linear(m) => m,
                    SpatialMagneticMaterial::Recoil(m) => m.material().material(),
                };
                if source.material() != magnetic {
                    return Err(invalid(
                        "conductivity must retain the exact original magnetic material bundle",
                    ));
                }
                Ok(())
            },
            cancelled,
        )?;
        let magnetic = coupling.restrict(&problem.stiffness, cancelled)?;
        field(problem, &potential_wb, cancelled)?;
        Ok(Self {
            prepared: Arc::new(Prepared {
                problem: MagneticOwner::Borrowed(problem),
                coupling,
                magnetic,
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
    pub fn impressed_current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.impressed_current_a_m2
    }
    pub fn tree_edges(&self) -> &[usize] {
        &self.prepared.tree_edges
    }
    pub fn problem(&self) -> &SpatialVectorMagnetostatics {
        &self.prepared.problem
    }
    /// Recover the accepted physical B without a static re-equilibration.
    pub fn flux_density_t(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Vec<[f64; 3]>> {
        geometry(self.problem()).flux(&self.potential_wb, cancelled)
    }
    /// Recover B/H from accepted A and its original source, with no static solve.
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        field(self.problem(), &self.potential_wb, cancelled)
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
        self.step_prepared_candidate(
            Arc::clone(&self.prepared),
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared_candidate(
        &self,
        prepared: Arc<Prepared<'m>>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionCandidate<'m>> {
        poll(cancelled)?;
        let p = &prepared;
        let problem = &*p.problem;
        let original = self.problem();
        let mesh = problem.mesh();
        let n = p.mass.nrows();
        let g = geometry(problem);
        let (time_s, restart) = admit_step(
            g,
            &p.coupling,
            self.time_s,
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
        )?;
        let endpoint_load = p.current_load(
            g,
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
        let load = p.current_load(
            g,
            &midpoint_current,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let old_field = field(original, &self.potential_wb, cancelled)?;
        let stored_before_j = old_field.stored_energy_j;
        let b0 = old_field.flux_density_t;
        let old_magnetic = apply(&original.stiffness, &self.potential_wb, cancelled)?;
        let refreshed_magnetic = apply(&problem.stiffness, &self.potential_wb, cancelled)?;
        let mut rhs = vec![0.0; n];
        for (e, index) in p.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                rhs[i] = finite(
                    load[e] + f64::midpoint(original.remanence_load[e], problem.remanence_load[e])
                        - f64::midpoint(old_magnetic[e], refreshed_magnetic[e]),
                )?;
            }
        }
        let (motion, drive, maximum_cell_magnetic_peclet) = p.motion(
            g,
            velocity_m_s,
            &b0,
            &problem.permeability,
            controls.maximum_cell_magnetic_peclet,
            cancelled,
        )?;
        for (rhs, drive) in rhs.iter_mut().zip(drive) {
            *rhs = finite(*rhs + drive)?;
        }
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
                magnetic_scale: 0.5,
                motion_scale: 0.5,
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
        let stage_flux_density_t = field(problem, &mid, cancelled)?.flux_density_t;
        let endpoint_field = field(problem, &next, cancelled)?;
        let stored_after_j = endpoint_field.stored_energy_j;
        let endpoint_flux_density_t = endpoint_field.flux_density_t;
        let endpoint_magnetic = apply(&problem.stiffness, &next, cancelled)?;
        let current = p.recover(
            g,
            &delta,
            &x,
            &stage_flux_density_t,
            velocity_m_s,
            duration_s,
            cancelled,
        )?;
        let maximum_charge_residual_a = current.maximum_charge_residual_a;
        let mut maximum_ampere_residual_a = 0.0_f64;
        for e in 0..mesh.edges.len() {
            poll(cancelled)?;
            if !p.fixed_edges[e] {
                let residual = finite(
                    f64::midpoint(old_magnetic[e], endpoint_magnetic[e])
                        - f64::midpoint(original.remanence_load[e], problem.remanence_load[e])
                        - load[e]
                        - current.edge_current_a[e],
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
        let mut magnetic_cells = Vec::with_capacity(mesh.tets.len());
        let mut parameter_work_j = 0.0;
        for cell in 0..mesh.tets.len() {
            poll(cancelled)?;
            let volume = problem.geometry.geometry.vol_signed[cell].abs();
            let old_mu = original.permeability[cell];
            let new_mu = problem.permeability[cell];
            let old_h =
                std::array::from_fn(|j| (b0[cell][j] - original.remanence[cell][j]) / old_mu);
            let new_h = std::array::from_fn(|j| {
                (endpoint_flux_density_t[cell][j] - problem.remanence[cell][j]) / new_mu
            });
            let work = finite(
                0.5 * volume
                    * dot(b0[cell], endpoint_flux_density_t[cell])
                    * ((old_mu - new_mu) / old_mu / new_mu),
            )?;
            parameter_work_j = finite(parameter_work_j + work)?;
            magnetic_cells.push(SpatialInductionMagneticCell {
                cell,
                volume_m3: volume,
                flux_density_t: [b0[cell], endpoint_flux_density_t[cell]],
                free_energy_j: [
                    finite(0.5 * volume * old_mu * dot(old_h, old_h))?,
                    finite(0.5 * volume * new_mu * dot(new_h, new_h))?,
                ],
                parameter_work_j: work,
            });
        }
        let energy_residual_j = finite(
            stored_after_j - stored_before_j + current.heat_j + current.mechanical_work_j
                - impressed_work_j
                - parameter_work_j,
        )?;
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "magnetic, Ohmic, Lorentz and impressed work do not close",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialInductionCandidate {
            state: Self {
                prepared,
                potential_wb: next,
                impressed_current_a_m2: endpoint_current_a_m2.to_vec(),
                time_s,
            },
            step: SpatialInductionStep {
                cells: current.cells,
                stage_flux_density_t,
                endpoint_flux_density_t,
                stored_before_j,
                stored_after_j,
                magnetic_cells,
                parameter_work_j,
                heat_j: current.heat_j,
                mechanical_work_j: current.mechanical_work_j,
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
struct Increment<'a> {
    magnetic: &'a Csr,
    mass: &'a Csr,
    motion: &'a Csr,
    transpose: &'a Csr,
    dt: f64,
    magnetic_scale: f64,
    motion_scale: f64,
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
            y[i] = self.magnetic_scale * row(self.magnetic) + row(self.mass) / self.dt
                - self.motion_scale * row(motion);
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
