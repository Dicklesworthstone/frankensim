//! Eulerian planar moving-conductor induction on the original field owner.
//!
//! With B=(A_y,-A_x), (v cross B)_z=-v.grad(A). Midpoint evolution solves
//! `K A_mid + P_sigma(Adot + v.grad(A_mid)) = load(J_impressed_mid)`.
//! P_sigma is the SAME exact P1 mass and floating-region mean projection as
//! stationary induction. The generally nonsymmetric operator uses FGMRES.
//! Current, Joule heat and J cross B use that same stage field, closing
//! `delta W + Q + W_mechanical = W_impressed`.
//!
//! Velocity is prescribed at original mesh vertices and held over the step.
//! Geometry and Eulerian material coefficients remain fixed. This represents
//! flow through a stationary homogeneous conductor domain or motion leaving
//! its occupancy/material field invariant; it does not move a finite body,
//! material interfaces or magnetic sources. It is not an ALE/sliding-mesh,
//! terminal-current, nonlinear or permanent-magnet evolution model.
//! The nonrelativistic MQS/Lorentz convention and fixed-domain restriction are
//! described in COMSOL's AC/DC guide, "Velocity (Lorentz Term)":
//! https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html.

use super::{
    Conduction, EddyCellHeat, EddyConductor, EddyError, PlanarEddy, admit_algebraic, bind, finite,
    load, magnetic_load, poll, positive,
};
use crate::planar::{MagneticError, gradients};
use fs_rep_mesh::TriComplex2;
use fs_solver::{
    FgmresState, FlexiblePreconditioner, LinearOp, ResidualClaim, SolveReport, SolverRunProgress,
};
use fs_sparse::Csr;
use std::{
    cell::{Cell, RefCell},
    collections::BTreeMap,
};

/// Explicit Krylov, spatial-resolution and physical-accounting allowances.
#[derive(Clone, Copy, Debug)]
pub struct MovingEddyControls {
    /// True relative Euclidean residual of the accepted potential increment.
    pub relative_tolerance: f64,
    /// Absolute free-coordinate weak residual allowance [A m]. The actual
    /// target is max(this, relative_tolerance * norm(rhs)).
    pub absolute_residual_a_m: f64,
    /// Arnoldi vectors per restart, bounded additionally by free DOF count.
    pub restart: usize,
    /// Maximum restart cycles; no failed solve publishes a candidate.
    pub max_cycles: usize,
    /// Maximum retained Krylov scalars, including the Hessenberg matrix.
    pub max_krylov_scalars: usize,
    /// Caller-admitted upper bound on mu*sigma*max_vertex_speed*cell_diameter/2.
    /// Exceeding it requires mesh refinement or a different spatial formulation.
    pub max_cell_magnetic_peclet: f64,
    pub energy_tolerance_j: f64,
    pub current_balance_tolerance_a: f64,
}

pub(super) fn admit_controls(
    limits: MovingEddyControls,
    free_count: usize,
    dt: f64,
) -> Result<usize, EddyError> {
    if !positive(dt)
        || !positive(limits.relative_tolerance)
        || limits.relative_tolerance >= 1.0
        || !positive(limits.absolute_residual_a_m)
        || limits.restart == 0
        || limits.max_cycles == 0
        || !positive(limits.max_cell_magnetic_peclet)
        || !positive(limits.energy_tolerance_j)
        || !positive(limits.current_balance_tolerance_a)
    {
        return Err(EddyError::Invalid(
            "positive finite moving-induction interval and controls required",
        ));
    }
    let restart = limits.restart.min(free_count.max(1));
    // Includes staged checkpoints, Arnoldi work vectors, rotations, history,
    // and triangular-solve coefficients, for the existing FGMRES owner.
    let scalars = restart
        .checked_mul(2)
        .and_then(|v| v.checked_add(12))
        .and_then(|v| v.checked_mul(free_count))
        .and_then(|v| {
            restart
                .checked_add(8)
                .and_then(|m| m.checked_mul(restart))
                .and_then(|h| v.checked_add(h))
        })
        .and_then(|v| {
            limits
                .max_cycles
                .checked_mul(4)
                .and_then(|history| v.checked_add(history))
        });
    let scalars = if free_count == 0 { Some(0) } else { scalars };
    if scalars.is_none_or(|n| n > limits.max_krylov_scalars)
        || restart.checked_mul(limits.max_cycles).is_none()
    {
        return Err(EddyError::Invalid(
            "moving-induction Krylov storage budget exceeded",
        ));
    }
    Ok(restart)
}

/// Stage fields and exact P1 mechanical-power quadrature on one original cell.
#[derive(Clone, Debug)]
pub struct MovingEddyCell {
    /// `electric_field_v_m` here is E+v cross B, in the conductor frame.
    pub joule: EddyCellHeat,
    pub laboratory_electric_field_v_m: [f64; 3],
    pub motional_electric_field_v_m: [f64; 3],
    /// J_z e_z cross B at each local vertex, in original planar axes [N/m³].
    pub lorentz_force_density_n_m3: [[f64; 2]; 3],
    /// Positive work is electromagnetic work delivered to prescribed motion.
    pub mechanical_work_j: f64,
}

#[derive(Clone, Debug)]
pub struct MovingEddyStep {
    pub cells: Vec<MovingEddyCell>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
    pub supplied_work_j: f64,
    /// Delta stored + Joule heat + mechanical work - impressed source work.
    pub energy_residual_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
    pub max_cell_magnetic_peclet: f64,
    pub report: Option<SolveReport>,
}

pub struct MovingEddyCandidate<'m> {
    state: PlanarEddy<'m>,
    step: MovingEddyStep,
}
impl<'m> MovingEddyCandidate<'m> {
    pub fn state(&self) -> &PlanarEddy<'m> {
        &self.state
    }
    pub fn step(&self) -> &MovingEddyStep {
        &self.step
    }
    /// The receiving mechanical and thermal owners must close the reported work.
    pub fn into_parts(self) -> (PlanarEddy<'m>, MovingEddyStep) {
        (self.state, self.step)
    }
}

pub(super) struct MotionCell {
    gradients: [[f64; 2]; 3],
    velocity: [[f64; 2]; 3],
}
pub(super) struct Motion {
    cells: Vec<MotionCell>,
    pub(super) max_peclet: f64,
}
impl Motion {
    pub(super) fn new(
        mesh: &TriComplex2,
        conduction: &Conduction,
        velocity: &[[f64; 2]],
        maximum_permeability_h_m: &[f64],
        limit: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EddyError> {
        if velocity.len() != mesh.vertices().len()
            || velocity.iter().flatten().any(|v| !v.is_finite())
            || maximum_permeability_h_m.len() != conduction.cells.len()
            || maximum_permeability_h_m.iter().any(|v| !positive(*v))
        {
            return Err(EddyError::Invalid(
                "finite planar velocity is required at every original mesh vertex",
            ));
        }
        let mut cells = Vec::with_capacity(conduction.cells.len());
        let mut max_peclet: f64 = 0.0;
        for (&permeability, cell) in maximum_permeability_h_m.iter().zip(&conduction.cells) {
            poll(cancelled)?;
            let points = cell.vertices.map(|v| mesh.vertices()[v]);
            let velocity = cell.vertices.map(|v| velocity[v]);
            let speed = velocity
                .iter()
                .map(|v| v[0].hypot(v[1]))
                .fold(0.0, f64::max);
            let diameter = (0..3)
                .map(|i| {
                    let a = points[i];
                    let b = points[(i + 1) % 3];
                    (a[0] - b[0]).hypot(a[1] - b[1])
                })
                .fold(0.0, f64::max);
            let peclet = finite(permeability * cell.conductivity_s_m * speed * (0.5 * diameter))?;
            if peclet > limit {
                return Err(EddyError::Invalid(
                    "moving conductor exceeds the declared cell magnetic Peclet limit",
                ));
            }
            max_peclet = max_peclet.max(peclet);
            cells.push(MotionCell {
                gradients: gradients(points)?,
                velocity,
            });
        }
        Ok(Self { cells, max_peclet })
    }

    pub(super) fn offsets(
        &self,
        conduction: &Conduction,
        potential: &[f64],
        factor: f64,
        output: &mut [[f64; 3]],
        check: &dyn Fn() -> bool,
    ) -> bool {
        for ((cell, motion), output) in conduction.cells.iter().zip(&self.cells).zip(output) {
            if check() {
                return false;
            }
            let anchor = potential[cell.vertices[0]];
            let mut grad = [0.0; 2];
            for i in 1..3 {
                let delta = potential[cell.vertices[i]] - anchor;
                for (axis, value) in grad.iter_mut().enumerate() {
                    *value = delta.mul_add(motion.gradients[i][axis], *value);
                }
            }
            for (value, velocity) in output.iter_mut().zip(&motion.velocity) {
                *value = factor * velocity[0].mul_add(grad[0], velocity[1] * grad[1]);
            }
        }
        true
    }
}

pub(super) struct Scratch {
    pub(super) full: Vec<f64>,
    pub(super) rate: Vec<f64>,
    pub(super) conduction: Vec<f64>,
    pub(super) offsets: Vec<[f64; 3]>,
}
pub(super) struct MovingIncrement<'a> {
    pub(super) magnetic: &'a Csr,
    pub(super) stage_factor: f64,
    pub(super) free: &'a [usize],
    pub(super) conduction: &'a Conduction,
    pub(super) motion: &'a Motion,
    pub(super) dt: f64,
    pub(super) check: &'a dyn Fn() -> bool,
    pub(super) scratch: RefCell<Scratch>,
}

// Positive diffusion/storage diagonal, leaving the nonsymmetric velocity term
// to FGMRES. Regional subtraction is the exact floating-mass diagonal; omitting
// it would precondition the wrong electrical constraint.
pub(super) struct MovingJacobi(Vec<f64>);
impl MovingJacobi {
    pub(super) fn new(
        magnetic: &Csr,
        free: &[usize],
        conduction: &Conduction,
        dt: f64,
        stage_factor: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EddyError> {
        let mut diagonal = vec![0.0; magnetic.nrows()];
        let mut moments = BTreeMap::<(usize, usize), f64>::new();
        for cell in &conduction.cells {
            poll(cancelled)?;
            for &v in &cell.vertices {
                diagonal[v] += cell.sigma_volume / 6.0;
                *moments.entry((v, cell.group)).or_default() += cell.sigma_volume / 3.0;
            }
        }
        for ((v, group), moment) in moments {
            poll(cancelled)?;
            diagonal[v] -= moment * (moment / conduction.weights[group]);
        }
        let mut inverse = Vec::with_capacity(free.len());
        for &v in free {
            poll(cancelled)?;
            let (cols, vals) = magnetic.row(v);
            let stiffness = cols
                .iter()
                .position(|&col| col == v)
                .map_or(0.0, |i| vals[i]);
            let value = stage_factor * stiffness + diagonal[v] / dt;
            if !positive(value) || !positive(value.recip()) {
                return Err(EddyError::Invalid(
                    "moving-induction diffusion/storage diagonal is not representably positive",
                ));
            }
            inverse.push(value.recip());
        }
        Ok(Self(inverse))
    }
}
impl FlexiblePreconditioner for MovingJacobi {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        for ((out, &value), &scale) in output.iter_mut().zip(residual).zip(&self.0) {
            *out = value * scale;
        }
    }
}

impl LinearOp for MovingIncrement<'_> {
    fn n(&self) -> usize {
        self.free.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        // Reuse full-mesh buffers through every Arnoldi action; the sparse
        // stiffness and regional projection never become dense matrices.
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            full,
            rate,
            conduction,
            offsets,
        } = &mut *scratch;
        full.fill(0.0);
        rate.fill(0.0);
        for (&v, &value) in self.free.iter().zip(x) {
            full[v] = value;
            rate[v] = value / self.dt;
        }
        if !self.motion.offsets(
            self.conduction,
            full,
            self.stage_factor,
            offsets,
            self.check,
        ) || !self
            .conduction
            .apply_with_offsets(rate, offsets, conduction, self.check)
        {
            y.fill(f64::NAN);
            return;
        }
        for (row, &v) in self.free.iter().enumerate() {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let (cols, vals) = self.magnetic.row(v);
            let stiffness = cols
                .iter()
                .zip(vals)
                .map(|(&j, &a)| a * full[j])
                .sum::<f64>();
            y[row] = self.stage_factor * stiffness + conduction[v];
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            full,
            rate,
            conduction,
            ..
        } = &mut *scratch;
        full.fill(0.0);
        for (&v, &value) in self.free.iter().zip(x) {
            full[v] = value;
        }
        if !self.conduction.apply(full, conduction, self.check) {
            y.fill(f64::NAN);
            return;
        }
        let Some(means) = self.conduction.means(full, self.check) else {
            y.fill(f64::NAN);
            return;
        };
        rate.fill(0.0);
        for (cell, motion) in self.conduction.cells.iter().zip(&self.motion.cells) {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let anchor = full[self.conduction.anchors[cell.group]];
            let d = cell
                .vertices
                .map(|v| (full[v] - anchor) - means[cell.group]);
            let sum = d.iter().sum::<f64>();
            let mut weighted_velocity = [0.0; 2];
            for (i, &d) in d.iter().enumerate() {
                for (axis, value) in weighted_velocity.iter_mut().enumerate() {
                    *value += cell.sigma_volume / 12.0 * (d + sum) * motion.velocity[i][axis];
                }
            }
            for (i, &v) in cell.vertices.iter().enumerate() {
                rate[v] += motion.gradients[i][0].mul_add(
                    weighted_velocity[0],
                    motion.gradients[i][1] * weighted_velocity[1],
                );
            }
        }
        for (row, &v) in self.free.iter().enumerate() {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let (cols, vals) = self.magnetic.row(v);
            let stiffness = cols
                .iter()
                .zip(vals)
                .map(|(&j, &a)| a * full[j])
                .sum::<f64>();
            y[row] = self.stage_factor * (stiffness + rate[v]) + conduction[v] / self.dt;
        }
    }
}

pub(super) struct MotionHeat {
    pub(super) cells: Vec<MovingEddyCell>,
    pub(super) heat_j: f64,
    pub(super) mechanical_work_j: f64,
    pub(super) uniform_drive_work_j: f64,
    pub(super) motional_drive_work_j: f64,
    pub(super) region_motional_electric_field_v_m: Vec<(u32, f64)>,
    pub(super) region_currents_a: Vec<(u32, f64)>,
}

// One constitutive current and one exact P1 product quadrature serve midpoint
// linear induction and endpoint nonlinear induction. Only the magnetic stage
// and its own stored-energy/time-discretization account differ.
#[allow(clippy::too_many_arguments)]
pub(super) fn recover_motion(
    conduction: &Conduction,
    conductors: &[EddyConductor],
    motion: &Motion,
    stage_flux_density_t: &[[f64; 2]],
    rate: &[f64],
    dt: f64,
    thickness_m: f64,
    current_balance_tolerance_a: f64,
    drives: Option<(&[f64], &[f64])>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MotionHeat, EddyError> {
    let mut offsets = vec![[0.0; 3]; conduction.cells.len()];
    for ((offset, motion), source) in offsets.iter_mut().zip(&motion.cells).zip(conductors) {
        poll(cancelled)?;
        let b = stage_flux_density_t[source.triangle];
        let grad = [-b[1], b[0]];
        for (value, velocity) in offset.iter_mut().zip(&motion.velocity) {
            *value = finite(velocity[0].mul_add(grad[0], velocity[1] * grad[1]))?;
        }
    }
    let heat = {
        let cancel = RefCell::new(&mut *cancelled);
        conduction.heat_with_offsets(
            rate,
            &offsets,
            dt,
            thickness_m,
            current_balance_tolerance_a,
            drives,
            &|| (*cancel.borrow_mut())(),
        )?
    };
    let mut uniform_drive_work_j = 0.0;
    let mut motional_drive_work_j = 0.0;
    let mut region_motional_electric_field_v_m = Vec::new();
    if let Some((electric, currents)) = drives {
        let means = {
            let cancel = RefCell::new(&mut *cancelled);
            conduction.mean_offsets(&offsets, &|| (*cancel.borrow_mut())())?
        };
        region_motional_electric_field_v_m.reserve(conduction.regions.len());
        for (group, (&region, &mean)) in conduction.regions.iter().zip(&means).enumerate() {
            poll(cancelled)?;
            let dc =
                finite(dt * finite(thickness_m * finite(electric[group] * currents[group])?)?)?;
            if dc < 0.0 {
                return Err(EddyError::Balance(
                    "solid DC drive work must be nonnegative",
                ));
            }
            uniform_drive_work_j = finite(uniform_drive_work_j + dc)?;
            motional_drive_work_j = finite(
                motional_drive_work_j
                    + finite(dt * finite(thickness_m * finite(mean * currents[group])?)?)?,
            )?;
            // q uses +v.grad(A), while the physical (v cross B)_z is its
            // negative. The winding back-EMF is -depth*n times this mean.
            region_motional_electric_field_v_m.push((region, -mean));
        }
    }
    let mut mechanical_work_j = 0.0;
    let mut cells = Vec::with_capacity(heat.cells.len());
    for ((source, heat), offsets) in conductors.iter().zip(heat.cells).zip(&offsets) {
        poll(cancelled)?;
        let b = stage_flux_density_t[source.triangle];
        let grad = [-b[1], b[0]];
        let j = heat.current_density_a_m2;
        // Exact integral J*(v.grad A) = (J cross B).v, including the floating
        // conductor's potential gradient already present in this SAME current.
        let work = finite(
            dt * heat.volume_m3 / 12.0
                * (j.iter().sum::<f64>() * offsets.iter().sum::<f64>()
                    + j.iter().zip(offsets).map(|(a, b)| a * b).sum::<f64>()),
        )?;
        mechanical_work_j = finite(mechanical_work_j + work)?;
        let laboratory_electric_field_v_m =
            std::array::from_fn(|i| heat.electric_field_v_m[i] + offsets[i]);
        let lorentz_force_density_n_m3 = j.map(|j| grad.map(|g| j * g));
        if laboratory_electric_field_v_m
            .iter()
            .chain(lorentz_force_density_n_m3.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(EddyError::Invalid(
                "unrepresentable moving-conductor field or force",
            ));
        }
        cells.push(MovingEddyCell {
            joule: EddyCellHeat {
                source: source.clone(),
                electric_field_v_m: heat.electric_field_v_m,
                current_density_a_m2: j,
                volume_m3: heat.volume_m3,
                power_w: heat.power_w,
                heat_j: heat.heat_j,
            },
            laboratory_electric_field_v_m,
            motional_electric_field_v_m: offsets.map(|v| -v),
            lorentz_force_density_n_m3,
            mechanical_work_j: work,
        });
    }
    Ok(MotionHeat {
        cells,
        heat_j: heat.heat_j,
        mechanical_work_j,
        uniform_drive_work_j,
        motional_drive_work_j,
        region_motional_electric_field_v_m,
        region_currents_a: heat.region_currents_a,
    })
}

impl<'m> PlanarEddy<'m> {
    /// Continue the original history with prescribed Eulerian velocity in m/s.
    /// Magnetic material, conductor source/regions, mesh and fixed A boundary
    /// remain unchanged. Velocity is held at the midpoint; endpoint impressed
    /// current varies linearly from the accepted source. Floating electric
    /// potential gradients enforce zero net axial current in EACH conductor.
    ///
    /// Use only when the declared motion leaves material occupancy/coefficient
    /// fields fixed (including homogeneous through-flow with supplied boundary
    /// data). Translating interfaces, moving magnets and changing geometry need
    /// a moving-mesh/material owner and cannot be represented by this operation.
    pub fn step_with_eulerian_velocity_candidate(
        &self,
        target_s: f64,
        velocity_at_vertices_m_s: &[[f64; 2]],
        impressed_jz_a_m2: &[f64],
        limits: MovingEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingEddyCandidate<'m>, EddyError> {
        poll(cancelled)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite() {
            return Err(EddyError::Invalid(
                "positive finite moving-induction interval and controls required",
            ));
        }
        let restart = admit_controls(limits, self.free.len(), dt)?;
        let (_, conduction) = bind(
            &*self.problem,
            self.conductors.clone(),
            Some(&self.conductors),
            cancelled,
        )?;
        let motion = Motion::new(
            &self.problem.mesh,
            &conduction,
            velocity_at_vertices_m_s,
            &self
                .conductors
                .iter()
                .map(|c| self.problem.permeability_h_m[c.triangle])
                .collect::<Vec<_>>(),
            limits.max_cell_magnetic_peclet,
            cancelled,
        )?;
        let preconditioner = MovingJacobi::new(
            &self.problem.stiffness,
            &self.free,
            &conduction,
            dt,
            0.5,
            cancelled,
        )?;
        let endpoint_load = load(
            &*self.problem,
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
            &*self.problem,
            &midpoint_current,
            &self.components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let old_magnetic = magnetic_load(
            &self.problem,
            &self.potential_wb_m,
            &self.components,
            cancelled,
        )?;
        let n = self.potential_wb_m.len();
        let mut potential = self.potential_wb_m.clone();
        let mut offsets = vec![[0.0; 3]; conduction.cells.len()];
        let report;
        {
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
            let mut advected = vec![0.0; n];
            if !motion.offsets(&conduction, &self.potential_wb_m, 1.0, &mut offsets, &check)
                || !conduction.apply_with_offsets(&vec![0.0; n], &offsets, &mut advected, &check)
            {
                return Err(EddyError::Cancelled);
            }
            let rhs: Vec<_> = self
                .free
                .iter()
                .map(|&v| finite(loads[v] - old_magnetic[v] - advected[v]))
                .collect::<Result<_, _>>()?;
            let rhs_norm = rhs.iter().fold(0.0_f64, |norm, &value| norm.hypot(value));
            let tolerance = if rhs_norm == 0.0 {
                limits.relative_tolerance
            } else {
                limits
                    .relative_tolerance
                    .max((limits.absolute_residual_a_m / rhs_norm).min(f64::MAX))
            };
            if self.free.is_empty() {
                report = None;
            } else {
                let op = MovingIncrement {
                    magnetic: &self.problem.stiffness,
                    stage_factor: 0.5,
                    free: &self.free,
                    conduction: &conduction,
                    motion: &motion,
                    dt,
                    check: &check,
                    scratch: RefCell::new(Scratch {
                        full: vec![0.0; n],
                        rate: vec![0.0; n],
                        conduction: vec![0.0; n],
                        offsets: vec![[0.0; 3]; conduction.cells.len()],
                    }),
                };
                let mut solver = FgmresState::new(&rhs, restart);
                let result = solver.run_controlled(
                    &op,
                    &preconditioner,
                    &rhs,
                    tolerance,
                    limits.max_cycles,
                    &mut || check(),
                );
                if check() || result.progress == SolverRunProgress::Paused {
                    return Err(EddyError::Cancelled);
                }
                if result.progress != SolverRunProgress::Complete {
                    return Err(EddyError::Invalid(
                        "moving-induction solver callback failed",
                    ));
                }
                if !result.report.converged_euclidean() {
                    return Err(MagneticError::Nonconvergence(result.report).into());
                }
                let mut actual = Vec::with_capacity(self.free.len());
                for (&v, &delta) in self.free.iter().zip(&solver.x) {
                    potential[v] = finite(potential[v] + delta)?;
                    let change = finite(potential[v] - self.potential_wb_m[v])?;
                    if delta != 0.0 && change == 0.0 {
                        return Err(EddyError::Invalid(
                            "moving field increment disappeared in accepted potential",
                        ));
                    }
                    actual.push(change);
                }
                let mut applied = vec![0.0; rhs.len()];
                op.apply(&actual, &mut applied);
                if check() {
                    return Err(EddyError::Cancelled);
                }
                if applied.iter().any(|v| !v.is_finite()) {
                    return Err(EddyError::Invalid(
                        "moving induction operator is unrepresentable",
                    ));
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
                let accepted = SolveReport::from_claim(
                    result.report.iters,
                    ResidualClaim::TrueEuclidean(relative),
                    tolerance,
                    history,
                );
                if !accepted.converged_euclidean() {
                    return Err(MagneticError::Nonconvergence(accepted).into());
                }
                report = Some(accepted);
            }
        }
        let new_magnetic = magnetic_load(&self.problem, &potential, &self.components, cancelled)?;
        admit_algebraic(
            &*self.problem,
            &new_magnetic,
            &endpoint_load,
            &self.algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let field = self.problem.recover(potential.clone(), None, cancelled)?;
        let old_field = self.field(cancelled)?;
        let increments: Vec<_> = potential
            .iter()
            .zip(&self.potential_wb_m)
            .map(|(a, b)| a - b)
            .collect();
        let rate: Vec<_> = increments
            .iter()
            .map(|a| finite(a / dt))
            .collect::<Result<_, _>>()?;
        // Average endpoint B, avoiding addition of large common gauge offsets.
        let stage_b: Vec<_> = old_field
            .magnetic_flux_density_t()
            .iter()
            .zip(field.magnetic_flux_density_t())
            .map(|(a, b)| [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])])
            .collect();
        let heat = recover_motion(
            &conduction,
            &self.conductors,
            &motion,
            &stage_b,
            &rate,
            dt,
            self.problem.thickness_m,
            limits.current_balance_tolerance_a,
            None,
            cancelled,
        )?;
        let supplied_work_j = finite(loads.iter().zip(&increments).map(|(f, d)| f * d).sum())?;
        let energy_residual_j = finite(
            field.energy_j - self.energy_j + heat.heat_j + heat.mechanical_work_j - supplied_work_j,
        )?;
        if energy_residual_j.abs() > limits.energy_tolerance_j {
            return Err(EddyError::Balance(
                "moving magnetic storage, conductor-frame heat, Lorentz work and source work do not close",
            ));
        }
        let energy_j = field.energy_j;
        let step = MovingEddyStep {
            cells: heat.cells,
            stored_before_j: self.energy_j,
            stored_after_j: energy_j,
            heat_j: heat.heat_j,
            mechanical_work_j: heat.mechanical_work_j,
            supplied_work_j,
            energy_residual_j,
            region_currents_a: heat.region_currents_a,
            max_cell_magnetic_peclet: motion.max_peclet,
            report,
        };
        let mut state = self.clone();
        state.potential_wb_m = potential;
        state.energy_j = energy_j;
        state.impressed_jz_a_m2 = impressed_jz_a_m2.to_vec();
        state.time_s = target_s;
        poll(cancelled)?;
        Ok(MovingEddyCandidate { state, step })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planar::PlanarMagnetostatics;
    use crate::planar::eddy::conduction::ConductionCell;
    use fs_rep_mesh::{Metric2, TriComplex2};
    use fs_sparse::Coo;

    #[test]
    fn moving_operator_has_a_true_transpose_with_projected_variable_velocity() {
        let mesh = TriComplex2::from_indexed_triangles(
            "moving-adjoint-algebra",
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
            Metric2::planar(0.3).unwrap(),
        )
        .unwrap();
        // Only the adapter algebra is exercised here; the integration tests
        // separately use actual source-admitted FEEC magnetic stiffness.
        let mut diagonal = Coo::new(5, 5);
        for i in 0..5 {
            diagonal.push(i, i, 1.0 + i as f64);
        }
        let problem = PlanarMagnetostatics {
            mesh,
            materials: vec![],
            permeability_h_m: vec![],
            recoil_materials: vec![],
            stiffness: diagonal.assemble(),
            thickness_m: 0.3,
        };
        let conduction = Conduction {
            cells: problem
                .mesh
                .faces()
                .iter()
                .map(|f| ConductionCell {
                    vertices: f.map(|v| v as usize),
                    group: 0,
                    sigma_volume: 1.35,
                    conductivity_s_m: 18.0,
                    volume_m3: 0.075,
                })
                .collect(),
            regions: vec![1],
            weights: vec![5.4],
            anchors: vec![0],
        };
        let motion = Motion {
            cells: conduction
                .cells
                .iter()
                .map(|cell| {
                    let p = cell.vertices.map(|v| problem.mesh.vertices()[v]);
                    MotionCell {
                        gradients: gradients(p).unwrap(),
                        velocity: p.map(|p| [p[0] + 0.3 * p[1], -0.2 * p[0] + 1.1 * p[1]]),
                    }
                })
                .collect(),
            max_peclet: 0.0,
        };
        let check = || false;
        let op = MovingIncrement {
            magnetic: &problem.stiffness,
            stage_factor: 0.5,
            free: &[0, 1, 2, 3, 4],
            conduction: &conduction,
            motion: &motion,
            dt: 0.17,
            check: &check,
            scratch: RefCell::new(Scratch {
                full: vec![0.0; 5],
                rate: vec![0.0; 5],
                conduction: vec![0.0; 5],
                offsets: vec![[0.0; 3]; 4],
            }),
        };
        let x = [0.1, -0.7, 2.0, 1.5, -0.2];
        let y = [0.4, 1.7, -0.9, 0.2, 2.2];
        let mut ax = [0.0; 5];
        let mut aty = [0.0; 5];
        op.apply(&x, &mut ax);
        op.apply_transpose(&y, &mut aty);
        let lhs = ax.iter().zip(y).map(|(a, b)| a * b).sum::<f64>();
        let rhs = aty.iter().zip(x).map(|(a, b)| a * b).sum::<f64>();
        assert!((lhs - rhs).abs() < 1e-11, "{lhs} != {rhs}");
        op.apply_transpose(&x, &mut aty);
        assert!(ax.iter().zip(aty).any(|(a, b)| (a - b).abs() > 1e-5));
    }
}
