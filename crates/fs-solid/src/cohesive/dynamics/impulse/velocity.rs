//! The mass-only event kernel; constitutive history stays with the caller.

use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::{LinearOp, SolveReport};
use fs_sparse::precond::Precond;

#[cfg(test)]
#[path = "velocity/tests.rs"]
mod tests;

/// Explicit solver, workspace and physical acceptance limits.
#[derive(Debug, Clone, Copy)]
pub struct ImpulseControls {
    /// Native checked CG, including a shared iteration/correction budget.
    pub solve: CheckedCgConfig,
    /// Conservative scalar workspace bound for the numerical event, including
    /// native Krylov histories. Excludes operator/preconditioner-owned storage
    /// and the existing accepted-state clone returned by the motion owner.
    pub max_workspace_entries: usize,
    /// Each component of `M (v_after - v_before) - J` must pass this absolute
    /// allowance plus the relative allowance on that component's two terms.
    pub absolute_momentum_tolerance_ns: f64,
    pub relative_momentum_tolerance: f64,
    /// Work is checked at the event scale, never against a large unchanged
    /// background kinetic energy that could conceal a lost impulse.
    pub absolute_work_tolerance_j: f64,
    pub relative_work_tolerance: f64,
}

/// Incremental event balances in the original reduced coordinates. A geometric
/// trace is required to interpret these as Cartesian linear/angular momentum.
#[derive(Debug, Clone)]
pub struct ImpulseReport {
    /// Native recomputed-residual report for the solved velocity increment.
    pub solve: SolveReport,
    pub corrections: usize,
    /// Rechecked on the actual representable returned velocity increment.
    pub max_momentum_residual_ns: f64,
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    /// Evaluated as the mean velocity dotted with the actual mass increment,
    /// avoiding subtraction of nearly equal large endpoint energies.
    pub kinetic_energy_change_j: f64,
    /// Signed mean velocity dotted with the supplied generalized impulse.
    pub applied_impulse_work_j: f64,
    /// Kinetic change minus supplied work, accumulated from the local defects.
    pub work_residual_j: f64,
    pub work_tolerance_j: f64,
}

#[derive(Debug, Clone)]
pub enum ImpulseError {
    Invalid(&'static str),
    Cancelled,
    Budget(&'static str),
    Solve(SolveReport),
    Velocity {
        component: usize,
        value_m_s: f64,
    },
    Momentum {
        component: usize,
        residual_ns: f64,
        tolerance_ns: f64,
    },
    Work {
        residual_j: f64,
        tolerance_j: f64,
    },
}

impl core::fmt::Display for ImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive impulse refused: {self:?}")
    }
}
impl std::error::Error for ImpulseError {}

struct BorrowedMass<'m>(&'m dyn LinearOp);
impl LinearOp for BorrowedMass<'_> {
    fn n(&self) -> usize {
        self.0.n()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply_transpose(x, y);
    }
}

pub(super) fn velocity_candidate<P: Precond, C: FnMut() -> bool>(
    mass: &dyn LinearOp,
    old_velocity: &[f64],
    impulse: &[f64],
    preconditioner: &P,
    maximum_velocity: f64,
    controls: ImpulseControls,
    cancelled: &mut C,
) -> Result<(Vec<f64>, ImpulseReport), ImpulseError> {
    poll(cancelled)?;
    let n = mass.n();
    if n == 0 || old_velocity.len() != n || impulse.len() != n {
        return Err(ImpulseError::Invalid(
            "complete, nonempty reduced momentum required",
        ));
    }
    if !maximum_velocity.is_finite() || maximum_velocity <= 0.0 {
        return Err(ImpulseError::Invalid(
            "positive finite velocity limit required",
        ));
    }
    for tolerance in [
        controls.absolute_momentum_tolerance_ns,
        controls.relative_momentum_tolerance,
        controls.absolute_work_tolerance_j,
        controls.relative_work_tolerance,
    ] {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(ImpulseError::Invalid(
                "finite nonnegative physical tolerances required",
            ));
        }
    }
    // Account for all live numerical vectors and geometrically grown/cloned
    // native CG histories. Do this before any event allocation/operator work.
    let histories = controls
        .solve
        .max_iterations
        .checked_mul(4)
        .and_then(|v| {
            controls
                .solve
                .max_corrections
                .min(controls.solve.max_iterations)
                .checked_mul(2)
                .and_then(|h| v.checked_add(h))
        })
        .and_then(|v| v.checked_add(32));
    let workspace = n
        .checked_mul(16)
        .and_then(|v| histories.and_then(|h| v.checked_add(h)))
        .ok_or(ImpulseError::Budget("workspace size overflow"))?;
    if workspace > controls.max_workspace_entries {
        return Err(ImpulseError::Budget("impulse workspace"));
    }
    for (i, (&v, &j)) in old_velocity.iter().zip(impulse).enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        if !v.is_finite() || !j.is_finite() {
            return Err(ImpulseError::Invalid(
                "finite velocity and impulse required",
            ));
        }
        if v.abs() > maximum_velocity {
            return Err(ImpulseError::Velocity {
                component: i,
                value_m_s: v,
            });
        }
    }
    let solved = checked_cg(
        &BorrowedMass(mass),
        preconditioner,
        impulse,
        controls.solve,
        |_| poll(cancelled),
    )
    .map_err(|error| match error {
        CheckedCgError::InvalidInput(message) => ImpulseError::Invalid(message),
        CheckedCgError::Interrupted(error) => error,
    })?;
    if !solved.report.converged_euclidean() {
        return Err(ImpulseError::Solve(solved.report));
    }
    let mut velocity = Vec::with_capacity(n);
    let mut delta = solved.x;
    for (i, (dv, &old)) in delta.iter_mut().zip(old_velocity).enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        let value = old + *dv;
        if !value.is_finite() || value.abs() > maximum_velocity {
            return Err(ImpulseError::Velocity {
                component: i,
                value_m_s: value,
            });
        }
        // A mathematically solved increment can disappear when added to a
        // large velocity. Gate the actual returned motion, not the CG vector.
        *dv = finite(value - old)?;
        velocity.push(value);
    }
    let mut mass_delta = vec![0.0; n];
    poll(cancelled)?;
    mass.apply(&delta, &mut mass_delta);
    poll(cancelled)?;

    let mut change = Sum::default();
    let mut work = Sum::default();
    let mut defect_work = Sum::default();
    let mut work_scale: f64 = 0.0;
    let mut max_residual: f64 = 0.0;
    for i in 0..n {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        let actual = finite(mass_delta[i])?;
        let defect = finite(actual - impulse[i])?;
        let tolerance = finite(controls.relative_momentum_tolerance.mul_add(
            actual.abs().max(impulse[i].abs()),
            controls.absolute_momentum_tolerance_ns,
        ))?;
        if defect.abs() > tolerance {
            return Err(ImpulseError::Momentum {
                component: i,
                residual_ns: defect,
                tolerance_ns: tolerance,
            });
        }
        max_residual = max_residual.max(defect.abs());
        let mean = 0.5 * old_velocity[i] + 0.5 * velocity[i];
        let dk = finite(mean * actual)?;
        let w = finite(mean * impulse[i])?;
        change.add(dk)?;
        work.add(w)?;
        defect_work.add(finite(mean * defect)?)?;
        work_scale = finite(work_scale + dk.abs().max(w.abs()))?;
    }
    let kinetic_energy_change_j = change.total()?;
    let applied_impulse_work_j = work.total()?;
    let work_residual_j = defect_work.total()?;
    let work_tolerance_j = finite(
        controls
            .relative_work_tolerance
            .mul_add(work_scale, controls.absolute_work_tolerance_j),
    )?;
    if work_residual_j.abs() > work_tolerance_j {
        return Err(ImpulseError::Work {
            residual_j: work_residual_j,
            tolerance_j: work_tolerance_j,
        });
    }
    // Endpoint diagnostics use the same mass. They are not subtracted to admit
    // the event, and their magnitude cannot loosen either acceptance gate.
    let kinetic_energy_before_j = kinetic(mass, old_velocity, &mut mass_delta, cancelled)?;
    let kinetic_energy_after_j = kinetic(mass, &velocity, &mut mass_delta, cancelled)?;
    poll(cancelled)?;
    Ok((
        velocity,
        ImpulseReport {
            solve: solved.report,
            corrections: solved.corrections,
            max_momentum_residual_ns: max_residual,
            kinetic_energy_before_j,
            kinetic_energy_after_j,
            kinetic_energy_change_j,
            applied_impulse_work_j,
            work_residual_j,
            work_tolerance_j,
        },
    ))
}

fn kinetic<C: FnMut() -> bool>(
    mass: &dyn LinearOp,
    velocity: &[f64],
    applied: &mut [f64],
    cancelled: &mut C,
) -> Result<f64, ImpulseError> {
    poll(cancelled)?;
    mass.apply(velocity, applied);
    poll(cancelled)?;
    let mut energy = Sum::default();
    for (i, (&v, &p)) in velocity.iter().zip(applied.iter()).enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        energy.add(finite((0.5 * v) * finite(p)?)?)?;
    }
    let value = energy.total()?;
    if value < 0.0 {
        return Err(ImpulseError::Invalid(
            "negative kinetic energy from supplied mass",
        ));
    }
    Ok(value)
}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), ImpulseError> {
    if cancelled() {
        Err(ImpulseError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, ImpulseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ImpulseError::Invalid(
            "non-finite momentum or work arithmetic",
        ))
    }
}

/// Neumaier accumulation retains small event work amid opposing contributions.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), ImpulseError> {
        let next = finite(self.value + term)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= term.abs() {
                    (self.value - next) + term
                } else {
                    (term - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, ImpulseError> {
        finite(self.value + self.correction)
    }
}
