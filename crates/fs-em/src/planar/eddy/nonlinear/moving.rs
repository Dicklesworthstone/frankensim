//! Eulerian motion with the original reversible H(B) and fixed-recoil laws.
//!
//! Backward Euler solves K_H(A_next) + P_sigma(delta A/dt + v.grad A_next)
//! = load(J_next). Its actual constitutive Hessian plus the nonsymmetric
//! motion derivative is solved with FGMRES, with residual-decreasing Newton
//! backtracking. No tangent permeability replaces the source's stored energy.
//!
//! The SAME endpoint current supplies Joule heat and J cross B. Thus
//! delta psi + Q_Ohmic + W_mechanical + D_BE = W_impressed, where the convex
//! Bregman remainder D_BE is numerical dissipation, never physical heat.
//! The Lorentz/frame convention is that of the parent moving-induction owner.

use super::super::{
    Conduction, EddyError, admit_algebraic, bind, finite, load,
    moving::{
        Motion, MovingEddyCell, MovingEddyControls, MovingIncrement, MovingJacobi, Scratch,
        admit_controls, recover_motion,
    },
    poll, positive,
};
use super::{
    IsotropicMagneticMaterial, NonlinearEddyError, NonlinearEddyMagneticCell,
    NonlinearMagneticError, NonlinearMagneticReport, PlanarNonlinearEddy,
    linear_magnetic_permeability, numerical_loss_density,
};
use crate::planar::{MagneticError, PlanarMagnetostatics};
use fs_solver::{FgmresState, SolveReport, SolverRunProgress};
use fs_sparse::Csr;
use std::cell::{Cell, RefCell};

/// Separate true nonlinear residual, Newton work, and nonsymmetric Krylov limits.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearMovingEddyControls {
    /// FGMRES, spatial Peclet, current and physical energy allowances.
    pub moving: MovingEddyControls,
    /// Absolute free-coordinate physical residual [A m].
    pub absolute_residual_tolerance_a_m: f64,
    /// Relative to max(initial residual, impressed-current load).
    pub relative_residual_tolerance: f64,
    /// Zero admits only an already converged initial field.
    pub max_newton_iterations: usize,
    /// Maximum halvings after each full Newton trial.
    pub max_backtracks: usize,
}

#[derive(Clone, Debug)]
pub struct NonlinearMovingEddyStep {
    /// Conductor-frame current/heat, laboratory E, and endpoint Lorentz work.
    pub cells: Vec<MovingEddyCell>,
    /// Original constitutive energy and positive backward-Euler loss per cell.
    /// Parameter work is zero: all magnetic/material sources are fixed here.
    pub magnetic_cells: Vec<NonlinearEddyMagneticCell>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub heat_j: f64,
    pub mechanical_work_j: f64,
    pub supplied_work_j: f64,
    /// Numerical time-discretization loss; not a thermal credit.
    pub numerical_dissipation_j: f64,
    /// Delta psi + heat + mechanical work + numerical loss - source work.
    pub energy_residual_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
    /// Conservative Peclet bound over each complete admitted constitutive curve.
    pub max_cell_magnetic_peclet: f64,
    pub report: NonlinearMagneticReport,
}

pub struct NonlinearMovingEddyCandidate<'m> {
    state: PlanarNonlinearEddy<'m>,
    step: NonlinearMovingEddyStep,
}
impl<'m> NonlinearMovingEddyCandidate<'m> {
    pub fn state(&self) -> &PlanarNonlinearEddy<'m> {
        &self.state
    }
    pub fn step(&self) -> &NonlinearMovingEddyStep {
        &self.step
    }
    /// Publish the field alongside its distinct thermal and mechanical transfers.
    pub fn into_parts(self) -> (PlanarNonlinearEddy<'m>, NonlinearMovingEddyStep) {
        (self.state, self.step)
    }
}

fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0_f64, |sum, value| sum.hypot(value))
}

fn maximum_permeability(
    material: &IsotropicMagneticMaterial,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, EddyError> {
    let value = match material {
        IsotropicMagneticMaterial::Linear(source) => linear_magnetic_permeability(source)?,
        IsotropicMagneticMaterial::Recoil(law) => law.material().permeability_h_m(),
        IsotropicMagneticMaterial::Reversible(law) => {
            // H(0)=0 and the positive source slopes make H/B their positive
            // weighted average. The minimum segment slope therefore bounds
            // BOTH radial dH/dB and tangential H/B over the entire source range.
            let mut minimum = f64::INFINITY;
            for pair in law.source_path().knots.windows(2) {
                poll(cancelled)?;
                let slope = (pair[1].1.evidence.value.value - pair[0].1.evidence.value.value)
                    / (pair[1].0 - pair[0].0);
                minimum = minimum.min(slope * law.field_scale());
            }
            minimum.recip()
        }
    };
    if !positive(value) {
        return Err(EddyError::Invalid(
            "magnetic transport bound is unrepresentable",
        ));
    }
    Ok(value)
}

#[allow(clippy::too_many_arguments)]
fn residual(
    state: &PlanarNonlinearEddy<'_>,
    potential: &[f64],
    loads: &[f64],
    conduction: &Conduction,
    motion: &Motion,
    dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, NonlinearEddyError> {
    let mut result = state.problem.magnetic_load(potential, cancelled)?;
    let rate = potential
        .iter()
        .zip(&state.potential_wb_m)
        .map(|(a, b)| finite((a - b) / dt))
        .collect::<Result<Vec<_>, _>>()?;
    let mut offsets = vec![[0.0; 3]; conduction.cells.len()];
    let mut current_load = vec![0.0; potential.len()];
    {
        let cancel = RefCell::new(&mut *cancelled);
        let check = || (*cancel.borrow_mut())();
        if !motion.offsets(conduction, potential, 1.0, &mut offsets, &check)
            || !conduction.apply_with_offsets(&rate, &offsets, &mut current_load, &check)
        {
            return Err(NonlinearEddyError::Cancelled);
        }
    }
    for ((value, applied), load) in result.iter_mut().zip(current_load).zip(loads) {
        *value = finite(*value + applied - load)?;
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn direction(
    magnetic: &Csr,
    free: &[usize],
    conduction: &Conduction,
    motion: &Motion,
    rhs: &[f64],
    dt: f64,
    limits: MovingEddyControls,
    restart: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, SolveReport), NonlinearEddyError> {
    let preconditioner = MovingJacobi::new(magnetic, free, conduction, dt, 1.0, cancelled)?;
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
    let op = MovingIncrement {
        magnetic,
        stage_factor: 1.0,
        free,
        conduction,
        motion,
        dt,
        check: &check,
        scratch: RefCell::new(Scratch {
            full: vec![0.0; magnetic.nrows()],
            rate: vec![0.0; magnetic.nrows()],
            conduction: vec![0.0; magnetic.nrows()],
            offsets: vec![[0.0; 3]; conduction.cells.len()],
        }),
    };
    let rhs_norm = finite(norm(rhs.iter().copied()))?;
    let tolerance = limits
        .relative_tolerance
        .max((limits.absolute_residual_a_m / rhs_norm).min(f64::MAX));
    let mut solver = FgmresState::new(rhs, restart);
    let result = solver.run_controlled(
        &op,
        &preconditioner,
        rhs,
        tolerance,
        limits.max_cycles,
        &mut || check(),
    );
    if check() || result.progress == SolverRunProgress::Paused {
        return Err(NonlinearEddyError::Cancelled);
    }
    if result.progress != SolverRunProgress::Complete {
        return Err(EddyError::Invalid("moving-induction solver callback failed").into());
    }
    if !result.report.converged_euclidean() {
        return Err(MagneticError::Nonconvergence(result.report).into());
    }
    Ok((solver.x, result.report))
}

impl<'m> PlanarNonlinearEddy<'m> {
    /// Continue saturating-core or fixed-recoil induction with held Eulerian v.
    ///
    /// Velocity is prescribed at every original mesh vertex. The complete
    /// magnetic laws, remanence direction, conductor regions, mesh and A boundary
    /// stay fixed. Motion must preserve their Eulerian occupancy/coefficient
    /// fields; translating finite interfaces or rotating magnet directions need
    /// a different geometry/material transport owner. This method does not
    /// move magnets, infer hysteresis/temperature changes, or advance mechanics.
    ///
    /// First-order backward Euler uses endpoint impressed current and velocity.
    /// Its Newton linearization is nonsymmetric even for a convex H(B) law.
    /// Source-envelope, convergence, energy or cancellation failure returns no
    /// replacement of the accepted history. Numerical loss is reported separately.
    pub fn step_with_eulerian_velocity_candidate(
        &self,
        target_s: f64,
        velocity_at_vertices_m_s: &[[f64; 2]],
        impressed_jz_a_m2: &[f64],
        limits: NonlinearMovingEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearMovingEddyCandidate<'m>, NonlinearEddyError> {
        poll(cancelled)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite()
            || !positive(limits.absolute_residual_tolerance_a_m)
            || !limits.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&limits.relative_residual_tolerance)
        {
            return Err(
                EddyError::Invalid("finite nonlinear moving-induction controls required").into(),
            );
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_boundary(
            self.problem.mesh(),
            &self.boundary,
            cancelled,
        )?;
        let free: Vec<_> = (0..prescribed.len())
            .filter(|&v| prescribed[v].is_none())
            .collect();
        let restart = admit_controls(limits.moving, free.len(), dt)?;
        let (_, conduction) = bind(
            &*self.problem,
            self.conductors.clone(),
            Some(&self.conductors),
            cancelled,
        )?;
        let mut permeabilities = Vec::with_capacity(self.conductors.len());
        for cell in &self.conductors {
            poll(cancelled)?;
            permeabilities.push(maximum_permeability(
                &self.problem.materials()[cell.triangle],
                cancelled,
            )?);
        }
        let motion = Motion::new(
            self.problem.mesh(),
            &conduction,
            velocity_at_vertices_m_s,
            &permeabilities,
            limits.moving.max_cell_magnetic_peclet,
            cancelled,
        )?;
        // A fixed recoil source does not transport permanent magnetization.
        // Stationary magnets can bias adjacent moving conductors; silently
        // applying Lorentz motion inside a magnet would omit its source and
        // magnetization-boundary work, even if a fitted Br were held constant.
        for (triangle, material) in self.problem.materials().iter().enumerate() {
            poll(cancelled)?;
            if matches!(material, IsotropicMagneticMaterial::Recoil(_))
                && self.problem.mesh().faces()[triangle].iter().any(|&v| {
                    velocity_at_vertices_m_s[v as usize]
                        .iter()
                        .any(|&component| component != 0.0)
                })
            {
                return Err(EddyError::Invalid(
                    "fixed-recoil magnet cells require zero material velocity",
                )
                .into());
            }
        }
        let loads = load(
            &*self.problem,
            impressed_jz_a_m2,
            &self.components,
            limits.moving.current_balance_tolerance_a,
            cancelled,
        )?;
        let mut potential = self.potential_wb_m.clone();
        let mut r = residual(
            self,
            &potential,
            &loads,
            &conduction,
            &motion,
            dt,
            cancelled,
        )?;
        let initial_norm = finite(norm(free.iter().map(|&v| r[v])))?;
        let reference = initial_norm.max(finite(norm(free.iter().map(|&v| loads[v])))?);
        let target = limits
            .absolute_residual_tolerance_a_m
            .max(finite(limits.relative_residual_tolerance * reference)?);
        let mut report = NonlinearMagneticReport {
            newton_iterations: 0,
            linear_iterations: 0,
            backtracks: 0,
            residual_norm_a_m: initial_norm,
            residual_target_a_m: target,
        };
        loop {
            poll(cancelled)?;
            report.residual_norm_a_m = finite(norm(free.iter().map(|&v| r[v])))?;
            if report.residual_norm_a_m <= target {
                break;
            }
            if report.newton_iterations == limits.max_newton_iterations {
                return Err(NonlinearMagneticError::Nonconvergence(report).into());
            }
            let (_, tangent) = self.problem.magnetic_linearization(&potential, cancelled)?;
            let rhs: Vec<_> = free.iter().map(|&v| -r[v]).collect();
            let (delta, linear_report) = direction(
                &tangent,
                &free,
                &conduction,
                &motion,
                &rhs,
                dt,
                limits.moving,
                restart,
                cancelled,
            )?;
            report.linear_iterations = report
                .linear_iterations
                .checked_add(linear_report.iters)
                .ok_or(EddyError::Invalid("nonlinear linear-work count overflow"))?;
            let mut alpha = 1.0_f64;
            let mut accepted = None;
            for backtrack in 0..=limits.max_backtracks {
                poll(cancelled)?;
                if backtrack > 0 {
                    report.backtracks = report
                        .backtracks
                        .checked_add(1)
                        .ok_or(EddyError::Invalid("nonlinear backtrack count overflow"))?;
                }
                let mut trial = potential.clone();
                for (&v, &d) in free.iter().zip(&delta) {
                    trial[v] = finite(alpha.mul_add(d, potential[v]))?;
                }
                match residual(self, &trial, &loads, &conduction, &motion, dt, cancelled) {
                    Ok(next) => {
                        let error = finite(norm(free.iter().map(|&v| next[v])))?;
                        if error <= target
                            || error <= (1.0 - 1e-4 * alpha) * report.residual_norm_a_m
                        {
                            accepted = Some((trial, next));
                            break;
                        }
                    }
                    Err(NonlinearEddyError::Magnetic(
                        NonlinearMagneticError::FieldRange { .. }
                        | NonlinearMagneticError::RecoilRange { .. },
                    )) => {}
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some((next_potential, next_residual)) = accepted else {
                return Err(NonlinearMagneticError::Nonconvergence(report).into());
            };
            potential = next_potential;
            r = next_residual;
            report.newton_iterations += 1;
        }
        let magnetic = self.problem.magnetic_load(&potential, cancelled)?;
        admit_algebraic(
            &*self.problem,
            &magnetic,
            &loads,
            &self.algebraic,
            limits.moving.current_balance_tolerance_a,
            cancelled,
        )?;
        let field = self
            .problem
            .recover(potential.clone(), report.clone(), cancelled)?;
        let increments = potential
            .iter()
            .zip(&self.potential_wb_m)
            .map(|(a, b)| finite(a - b))
            .collect::<Result<Vec<_>, _>>()?;
        let rate = increments
            .iter()
            .map(|d| finite(d / dt))
            .collect::<Result<Vec<_>, _>>()?;
        let heat = recover_motion(
            &conduction,
            &self.conductors,
            &motion,
            field.magnetic_flux_density_t(),
            &rate,
            dt,
            self.problem.thickness_m(),
            limits.moving.current_balance_tolerance_a,
            cancelled,
        )?;
        let before = self.field(cancelled)?;
        let mut magnetic_cells = Vec::with_capacity(self.problem.materials().len());
        let mut numerical_dissipation_j = 0.0;
        for (triangle, material) in self.problem.materials().iter().enumerate() {
            poll(cancelled)?;
            let b0 = before.magnetic_flux_density_t()[triangle];
            let b1 = field.magnetic_flux_density_t()[triangle];
            let volume_m3 = self
                .problem
                .mesh()
                .face_measure(triangle)
                .ok_or(EddyError::Invalid("missing magnetic volume"))?;
            let numerical = finite(volume_m3 * numerical_loss_density(material, b0, b1)?)?;
            if numerical < 0.0 {
                return Err(EddyError::Balance(
                    "convex magnetic numerical loss must be nonnegative",
                )
                .into());
            }
            numerical_dissipation_j = finite(numerical_dissipation_j + numerical)?;
            magnetic_cells.push(NonlinearEddyMagneticCell {
                triangle,
                volume_m3,
                flux_density_t: [b0, b1],
                free_energy_j: [
                    before.cell_energy_j()[triangle],
                    field.cell_energy_j()[triangle],
                ],
                numerical_dissipation_j: numerical,
                parameter_work_j: 0.0,
            });
        }
        let supplied_work_j = finite(loads.iter().zip(&increments).map(|(a, b)| a * b).sum())?;
        let energy_residual_j = finite(
            field.energy_j() - self.energy_j
                + heat.heat_j
                + heat.mechanical_work_j
                + numerical_dissipation_j
                - supplied_work_j,
        )?;
        if energy_residual_j.abs() > limits.moving.energy_tolerance_j {
            return Err(EddyError::Balance("nonlinear moving storage, Ohmic heat, Lorentz work, numerical loss and source work do not close").into());
        }
        let step = NonlinearMovingEddyStep {
            cells: heat.cells,
            magnetic_cells,
            stored_before_j: self.energy_j,
            stored_after_j: field.energy_j(),
            heat_j: heat.heat_j,
            mechanical_work_j: heat.mechanical_work_j,
            supplied_work_j,
            numerical_dissipation_j,
            energy_residual_j,
            region_currents_a: heat.region_currents_a,
            max_cell_magnetic_peclet: motion.max_peclet,
            report: report.clone(),
        };
        let mut state = self.clone();
        state.potential_wb_m = potential;
        state.energy_j = field.energy_j();
        state.impressed_jz_a_m2 = impressed_jz_a_m2.to_vec();
        state.time_s = target_s;
        state.report = report;
        poll(cancelled)?;
        Ok(NonlinearMovingEddyCandidate { state, step })
    }
}
