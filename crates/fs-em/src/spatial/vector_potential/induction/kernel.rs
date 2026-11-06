//! Shared endpoint A–phi solve with an owned constitutive trial response.
//!
//! Constitutive callbacks retain their own accepted memory. This kernel owns
//! the original Newton/FGMRES iteration, exact tangent assembly, rounded
//! cochain increments, full charge/Ampere checks and div-curl admission.

use super::super::super::{SpatialMagneticError, finite, invalid, poll};
use super::super::apply;
use super::{
    Increment, Jacobi, SpatialInductionField,
    coupling::{Coupling, Current, Geometry},
    diagonal,
    nonlinear::{SpatialNonlinearInductionControls, SpatialNonlinearInductionReport},
    terminals::{SpatialInductionTerminalTransfer, TerminalDrive},
};
use fs_feec::{
    assembly::stiffness,
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_qty::Dims;
use fs_solver::{FgmresState, LinearOp, ResidualClaim, SolveReport, SolverRunProgress};
use fs_sparse::{Csr, ops};
use std::cell::{Cell, RefCell};

pub(super) struct Magnetic<S = ()> {
    pub field: SpatialInductionField,
    pub energy_j: Vec<f64>,
    pub effort_a: Vec<f64>,
    pub tangent: Vec<CellWeight>,
    /// Actual constitutive candidate; never a mutation of accepted history.
    #[cfg_attr(not(feature = "magnetic-hysteresis"), allow(dead_code))]
    pub state: S,
}

struct Evaluation<S> {
    magnetic: Magnetic<S>,
    potential: Vec<f64>,
    residual: Vec<f64>,
    norm: f64,
}

pub(super) struct Solution<S> {
    pub magnetic: Magnetic<S>,
    pub potential: Vec<f64>,
    pub unknowns: Vec<f64>,
    pub report: SpatialNonlinearInductionReport,
}

pub(super) struct Physical {
    pub delta: Vec<f64>,
    pub current: Current,
    pub maximum_ampere_residual_a: f64,
    pub maximum_charge_residual_a: f64,
    pub terminals: Vec<SpatialInductionTerminalTransfer>,
    pub terminal_work_j: f64,
}

pub(super) struct Kernel<'a, C: SpatialConductorMaterial> {
    pub geometry: Geometry<'a>,
    pub coupling: &'a Coupling<C>,
    pub accepted_potential: &'a [f64],
    pub duration_s: f64,
    pub motion: &'a Csr,
    pub drive: &'a [f64],
    pub load: &'a [f64],
    pub assembly: WeightedAssemblyLimits,
    pub controls: SpatialNonlinearInductionControls,
    pub restart: usize,
}

impl<C: SpatialConductorMaterial> Kernel<'_, C> {
    fn evaluate<S, E: From<SpatialMagneticError>>(
        &self,
        x: &mut [f64],
        response: &mut impl FnMut(&[f64], &mut dyn FnMut() -> bool) -> Result<Magnetic<S>, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation<S>, E> {
        let mut potential = self.accepted_potential.to_vec();
        for (edge, index) in self.coupling.edge_unknown.iter().enumerate() {
            poll(cancelled)?;
            if let Some(i) = *index {
                potential[edge] = finite(self.accepted_potential[edge] + x[i])?;
                // Every trial and physical gate sees the actually retained A.
                x[i] = finite(potential[edge] - self.accepted_potential[edge])?;
            }
        }
        let magnetic = response(&potential, &mut || cancelled())?;
        let mass = apply(&self.coupling.mass, x, cancelled)?;
        let movement = apply(self.motion, x, cancelled)?;
        let mut residual: Vec<_> = mass
            .iter()
            .zip(&movement)
            .zip(self.drive)
            .map(|((m, v), d)| finite(m / self.duration_s - v - d))
            .collect::<super::Result<_>>()?;
        for (edge, index) in self.coupling.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                residual[i] = finite(residual[i] + magnetic.effort_a[edge] - self.load[edge])?;
            }
        }
        let norm = finite(residual.iter().fold(0.0_f64, |n, v| n.hypot(*v)))?;
        Ok(Evaluation {
            magnetic,
            potential,
            residual,
            norm,
        })
    }

    pub fn solve<S, E: From<SpatialMagneticError>>(
        &self,
        mut response: impl FnMut(&[f64], &mut dyn FnMut() -> bool) -> Result<Magnetic<S>, E>,
        retryable_source_refusal: impl Fn(&E) -> bool,
        nonconvergence: impl Fn(SpatialNonlinearInductionReport) -> E,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Solution<S>, E> {
        poll(cancelled)?;
        let settings = self.controls;
        let controls = settings.induction;
        if settings
            .maximum_backtracks
            .checked_add(1)
            .and_then(|trials| settings.maximum_newton_iterations.checked_mul(trials))
            .is_none()
            || settings
                .maximum_newton_iterations
                .checked_mul(controls.maximum_cycles)
                .and_then(|n| n.checked_mul(self.restart))
                .is_none()
        {
            return Err(invalid("nonlinear induction work bound overflows").into());
        }
        let transpose = ops::transpose(self.motion);
        let mut x = vec![0.0; self.coupling.mass.nrows()];
        let mut value = self.evaluate(&mut x, &mut response, cancelled)?;
        let initial_norm = value.norm;
        let target = controls
            .absolute_residual_a
            .max(finite(controls.relative_tolerance * initial_norm)?);
        let mut report = SpatialNonlinearInductionReport {
            newton_iterations: 0,
            linear_iterations: 0,
            backtracks: 0,
            residual_norm_a: value.norm,
            residual_target_a: target,
        };
        while value.norm > target {
            poll(cancelled)?;
            if report.newton_iterations >= settings.maximum_newton_iterations {
                return Err(nonconvergence(report));
            }
            let star = weighted_mass_matrix(
                self.geometry.mesh(),
                &self.geometry.geometry.geometry,
                2,
                &value.magnetic.tangent,
                Dims([-1, -1, 2, 0, 2, 0]),
                self.assembly,
                cancelled,
            )
            .map_err(|error| {
                if error == WeightedError::Cancelled {
                    SpatialMagneticError::Cancelled
                } else {
                    SpatialMagneticError::Assembly(error)
                }
            })?;
            let magnetic = self
                .coupling
                .restrict(&stiffness(self.geometry.curl, star.matrix()), cancelled)?;
            let inverse = (0..x.len())
                .map(|i| {
                    let d = finite(
                        diagonal(&magnetic, i) + diagonal(&self.coupling.mass, i) / self.duration_s,
                    )?;
                    if !super::positive(d) || !super::positive(d.recip()) {
                        Err(invalid(
                            "positive representable nonlinear induction Jacobi scale required",
                        ))
                    } else {
                        Ok(d.recip())
                    }
                })
                .collect::<super::Result<Vec<_>>>()?;
            let rhs: Vec<_> = value.residual.iter().map(|r| -r).collect();
            let forcing = 0.1_f64.min((value.norm / initial_norm).sqrt());
            let linear_target = (0.1 * target).max(forcing * value.norm);
            let (direction, linear) = linear_solve(
                &magnetic,
                &self.coupling.mass,
                self.motion,
                &transpose,
                self.duration_s,
                &rhs,
                inverse,
                linear_target,
                self.restart,
                controls.maximum_cycles,
                cancelled,
            )?;
            report.linear_iterations = report
                .linear_iterations
                .checked_add(linear.iters)
                .ok_or_else(|| invalid("nonlinear induction linear work count overflows"))?;
            let mut alpha = 1.0_f64;
            let mut accepted = None;
            let mut source_refusal = None;
            for backtrack in 0..=settings.maximum_backtracks {
                poll(cancelled)?;
                let mut trial = x
                    .iter()
                    .zip(&direction)
                    .map(|(a, d)| finite(alpha.mul_add(*d, *a)))
                    .collect::<super::Result<Vec<_>>>()?;
                match self.evaluate(&mut trial, &mut response, cancelled) {
                    Ok(next)
                        if next.norm <= target
                            || next.norm <= (1.0 - 1e-4 * alpha) * value.norm =>
                    {
                        accepted = Some((trial, next, backtrack));
                        break;
                    }
                    Ok(_) => {}
                    Err(error) if retryable_source_refusal(&error) => source_refusal = Some(error),
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some((next_x, next, backtracks)) = accepted else {
                return Err(source_refusal.unwrap_or_else(|| nonconvergence(report)));
            };
            x = next_x;
            value = next;
            report.newton_iterations += 1;
            report.backtracks = report
                .backtracks
                .checked_add(backtracks)
                .ok_or_else(|| invalid("nonlinear induction backtrack count overflows"))?;
            report.residual_norm_a = value.norm;
        }
        poll(cancelled)?;
        Ok(Solution {
            magnetic: value.magnetic,
            potential: value.potential,
            unknowns: x,
            report,
        })
    }

    /// One physical recovery for both constitutive owners, including equations
    /// eliminated from the algebraic system by its relative tree/charge gauges.
    pub fn recover<S>(
        &self,
        solution: &Solution<S>,
        velocity_m_s: &[[f64; 3]],
        terminal_drive: Option<&TerminalDrive<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> super::Result<Physical> {
        let c = self.coupling;
        let controls = self.controls.induction;
        let delta = solution
            .potential
            .iter()
            .zip(self.accepted_potential)
            .map(|(a, b)| finite(a - b))
            .collect::<super::Result<Vec<_>>>()?;
        let current = if let Some(drive) = terminal_drive {
            c.recover_with_potentials(
                self.geometry,
                &delta,
                &solution.unknowns,
                &solution.magnetic.field.flux_density_t,
                velocity_m_s,
                self.duration_s,
                drive.prescribed_potentials(),
                cancelled,
            )?
        } else {
            c.recover(
                self.geometry,
                &delta,
                &solution.unknowns,
                &solution.magnetic.field.flux_density_t,
                velocity_m_s,
                self.duration_s,
                cancelled,
            )?
        };
        let (terminals, terminal_work_j, maximum_charge_residual_a) =
            if let Some(drive) = terminal_drive {
                drive.transfers(
                    &current.charge_current_a,
                    self.duration_s,
                    &solution.unknowns,
                    cancelled,
                )?
            } else {
                (vec![], 0.0, current.maximum_charge_residual_a)
            };
        let mut maximum_ampere_residual_a = 0.0_f64;
        for edge in 0..delta.len() {
            poll(cancelled)?;
            if !c.fixed_edges[edge] {
                let residual = finite(
                    solution.magnetic.effort_a[edge]
                        - self.load[edge]
                        - current.edge_current_a[edge],
                )?;
                maximum_ampere_residual_a = maximum_ampere_residual_a.max(residual.abs());
            }
        }
        if maximum_ampere_residual_a > controls.current_tolerance_a
            || maximum_charge_residual_a > controls.current_tolerance_a
        {
            return Err(invalid(
                "full nonlinear Ampere or physical conductor charge equation failed",
            ));
        }
        let face_flux = apply(self.geometry.curl, &solution.potential, cancelled)?;
        for row in self.geometry.mesh().d2().rows {
            poll(cancelled)?;
            let flux = finite(row.iter().map(|&(f, s)| f64::from(s) * face_flux[f]).sum())?;
            if flux.abs() > controls.flux_tolerance_wb {
                return Err(invalid(
                    "nonlinear induction violates discrete div-curl allowance",
                ));
            }
        }
        Ok(Physical {
            delta,
            current,
            maximum_ampere_residual_a,
            maximum_charge_residual_a,
            terminals,
            terminal_work_j,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn linear_solve(
    magnetic: &Csr,
    mass: &Csr,
    motion: &Csr,
    transpose: &Csr,
    dt: f64,
    rhs: &[f64],
    inverse: Vec<f64>,
    target: f64,
    restart: usize,
    maximum_cycles: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> super::Result<(Vec<f64>, SolveReport)> {
    let cancel = RefCell::new(cancelled);
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
        magnetic,
        mass,
        motion,
        transpose,
        dt,
        magnetic_scale: 1.0,
        motion_scale: 1.0,
        check: &check,
    };
    let norm = finite(rhs.iter().fold(0.0_f64, |n, v| n.hypot(*v)))?;
    let tolerance = finite((0.25 * target / norm).min(0.5))?;
    if !super::positive(tolerance) {
        return Err(invalid("unrepresentable Newton linear residual target"));
    }
    let mut solver = FgmresState::new(rhs, restart);
    let result = solver.run_controlled(
        &operator,
        &Jacobi(inverse),
        rhs,
        tolerance,
        maximum_cycles,
        &mut || check(),
    );
    if result.progress != SolverRunProgress::Complete || check() {
        return Err(SpatialMagneticError::Cancelled);
    }
    let mut ax = vec![0.0; rhs.len()];
    operator.apply(&solver.x, &mut ax);
    if check() {
        return Err(SpatialMagneticError::Cancelled);
    }
    let residual = finite(rhs.iter().zip(ax).fold(0.0_f64, |n, (b, a)| n.hypot(b - a)))?;
    let mut history = result.report.history;
    history.push(residual / norm);
    let report = SolveReport::from_claim(
        result.report.iters,
        ResidualClaim::TrueEuclidean(residual / norm),
        target / norm,
        history,
    );
    if !report.converged_euclidean() {
        return Err(SpatialMagneticError::Nonconvergence(report));
    }
    Ok((solver.x, report))
}
