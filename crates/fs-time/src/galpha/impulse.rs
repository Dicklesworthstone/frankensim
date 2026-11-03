//! Instantaneous loading of the original structural mass operator.
//!
//! At a fixed configuration, M (v+ - v-) = J. Finite elastic and viscous
//! forces supply no impulse over this zero-duration event. They do determine
//! the new acceleration: M a+ = f+ - R(q, v+). Reusing the old acceleration
//! after a velocity jump would seed the next generalized-alpha step with a
//! residual of the wrong physical state, especially for damped contact models.

use super::initialization::{
    InitialSolveConfig, InitialSolveError, initial_workspace_components, second_order_acceleration,
    solve,
};
use super::{SecondOrderProblem, SecondOrderState};
use fs_math::dd::Dd;
use fs_solver::{FlexiblePreconditioner, SolveReport};

/// Independent allowances for the mass solve and the actual stored transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpulseControls {
    /// The original bounded mass solver, used for velocity and acceleration.
    pub mass_solve: InitialSolveConfig,
    /// Absolute componentwise allowance in the supplied generalized impulse units.
    /// Mixed translational/rotational coordinates require compatible scaling.
    pub absolute_impulse_tolerance: f64,
    /// Relative componentwise allowance against max(|J_i|, |M delta-v_i|).
    /// Existing background momentum never enlarges this allowance.
    pub relative_impulse_tolerance: f64,
    /// Absolute allowance on kinetic increment minus supplied impulse work.
    /// This is joules when coordinates and generalized impulses use SI units.
    pub work_tolerance: f64,
}

/// A refused event leaves the complete trajectory checkpoint unchanged.
#[derive(Clone, Debug)]
pub enum ImpulseError {
    Invalid(&'static str),
    /// The existing initialization/mass-solve owner refused the event.
    Solve(InitialSolveError),
    /// The rounded stored velocity cannot deliver the requested impulse.
    Momentum {
        component: usize,
        residual: f64,
        allowance: f64,
    },
    /// Actual kinetic increment and applied impulse work disagree.
    Work {
        residual: f64,
        allowance: f64,
    },
    NonFinite,
    Cancelled,
}
impl core::fmt::Display for ImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "structural impulse refused: {self:?}")
    }
}
impl core::error::Error for ImpulseError {}
impl From<InitialSolveError> for ImpulseError {
    fn from(error: InitialSolveError) -> Self {
        match error {
            InitialSolveError::Cancelled => Self::Cancelled,
            other => Self::Solve(other),
        }
    }
}

/// Actual accepted changes at the unchanged configuration and physical time.
#[derive(Clone, Debug)]
pub struct SecondOrderImpulse {
    /// Difference of the stored endpoint velocities, not the ideal solve vector.
    pub velocity_change: Vec<f64>,
    /// Original mass applied to that actual velocity change.
    pub momentum_change: Vec<f64>,
    /// Componentwise M delta-v - J in generalized impulse units.
    pub momentum_residual: Vec<f64>,
    /// Polarized change mean(v) dot M delta-v, avoiding subtraction of energies.
    pub kinetic_energy_change: f64,
    /// Signed applied work mean(v) dot J; a braking impulse can extract energy.
    pub supplied_work: f64,
    /// Independently accumulated mean(v) dot (M delta-v - J).
    pub work_residual: f64,
    pub velocity_solve: SolveReport,
    pub acceleration_solve: SolveReport,
}

/// Conservative live scalar-storage bound for this operation. Callback-owned
/// storage, caller inputs and the pre-existing trajectory history are excluded.
#[must_use]
pub fn impulse_workspace_components(n: usize, config: InitialSolveConfig) -> Option<usize> {
    let required = initial_workspace_components(n, 0, config)?.checked_add(n.checked_mul(8)?)?;
    (required <= isize::MAX as usize / core::mem::size_of::<f64>()).then_some(required)
}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), ImpulseError> {
    if cancelled() {
        Err(ImpulseError::Cancelled)
    } else {
        Ok(())
    }
}

/// Apply a prescribed generalized impulse to any admitted structural model.
///
/// The full mass action and caller-supplied mass preconditioner are reused;
/// there is no lumping, modal truncation, new integrator, or division by a
/// timestep. `forcing_after` is the finite force at the event's physical time
/// immediately after the impulse, distinct from J. Displacement, time, accepted
/// step count and history remain unchanged. Velocity and acceleration publish
/// together only after both solves, transfer/work checks and a final poll.
///
/// The physical model must have a constant symmetric positive-definite mass
/// in independent generalized coordinates, as for supported structural systems.
/// The mass solve is not a proof of that assumption. Singular constrained DAEs,
/// moving supports, collision/restitution laws, event localization and event
/// adjoints have separate owners. A caller must split a trajectory at this
/// event rather than replay it as a smooth uninterrupted trajectory.
///
/// Cost is two existing bounded matrix-free mass solves plus O(n) work/storage.
/// Cancellation is polled between solver cycles, around model callbacks, every
/// 256 transfer components and before publication. A model callback must bound
/// its own work. No viscous heat or elastic change occurs at this zero duration.
#[allow(clippy::too_many_arguments)]
pub fn second_order_impulse<M, P, C>(
    state: &mut SecondOrderState,
    model: &M,
    impulse: &[f64],
    forcing_after: &[f64],
    mass_preconditioner: &P,
    controls: ImpulseControls,
    max_workspace_components: usize,
    cancelled: &mut C,
) -> Result<SecondOrderImpulse, ImpulseError>
where
    M: SecondOrderProblem + ?Sized,
    P: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    poll(cancelled)?;
    let n = model.dimension();
    let config = controls.mass_solve;
    if n == 0
        || !state.t.is_finite()
        || config.restart == 0
        || config.max_cycles == 0
        || !config.tolerance.is_finite()
        || config.tolerance <= 0.0
        || !controls.absolute_impulse_tolerance.is_finite()
        || controls.absolute_impulse_tolerance < 0.0
        || !controls.relative_impulse_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_impulse_tolerance)
        || !controls.work_tolerance.is_finite()
        || controls.work_tolerance < 0.0
    {
        return Err(ImpulseError::Invalid(
            "finite state and valid impulse controls required",
        ));
    }
    let required = impulse_workspace_components(n, config).ok_or(ImpulseError::Invalid(
        "impulse workspace dimension overflow",
    ))?;
    if required > max_workspace_components {
        return Err(InitialSolveError::WorkspaceLimit {
            required,
            limit: max_workspace_components,
        }
        .into());
    }
    if [&state.q[..], &state.v, &state.a, impulse, forcing_after]
        .iter()
        .any(|values| values.len() != n || values.iter().any(|x| !x.is_finite()))
    {
        return Err(ImpulseError::Invalid(
            "finite dimension-matched state, impulse and force required",
        ));
    }
    let admission = model.validate_state(&state.q, &state.v);
    poll(cancelled)?;
    admission.map_err(InitialSolveError::ModelRefused)?;
    let velocity_solve = solve(
        impulse,
        |input, output| {
            model.mass_apply(input, output);
            Ok(())
        },
        mass_preconditioner,
        config,
        cancelled,
    )?;
    let mut velocity = vec![0.0; n];
    let mut change = vec![0.0; n];
    for i in 0..n {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        velocity[i] = state.v[i] + velocity_solve.value[i];
        change[i] = velocity[i] - state.v[i];
        if !velocity[i].is_finite() || !change[i].is_finite() {
            return Err(ImpulseError::NonFinite);
        }
    }
    let mut momentum = vec![f64::NAN; n];
    poll(cancelled)?;
    model.mass_apply(&change, &mut momentum);
    poll(cancelled)?;
    let mut residual = vec![0.0; n];
    let mut kinetic = Dd::ZERO;
    let mut work = Dd::ZERO;
    let mut work_residual = Dd::ZERO;
    for i in 0..n {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        residual[i] = momentum[i] - impulse[i];
        let allowance = controls.absolute_impulse_tolerance
            + controls.relative_impulse_tolerance * impulse[i].abs().max(momentum[i].abs());
        if !momentum[i].is_finite() || !residual[i].is_finite() || !allowance.is_finite() {
            return Err(ImpulseError::NonFinite);
        }
        if residual[i].abs() > allowance {
            return Err(ImpulseError::Momentum {
                component: i,
                residual: residual[i],
                allowance,
            });
        }
        // Half before adding avoids overflow in a representable mean velocity.
        let mean = Dd::from_f64(0.5 * state.v[i]) + Dd::from_f64(0.5 * velocity[i]);
        kinetic = kinetic + mean * Dd::from_f64(momentum[i]);
        work = work + mean * Dd::from_f64(impulse[i]);
        work_residual =
            work_residual + mean * (Dd::from_f64(momentum[i]) - Dd::from_f64(impulse[i]));
    }
    let kinetic = kinetic.to_f64();
    let work = work.to_f64();
    let work_residual = work_residual.to_f64();
    if [kinetic, work, work_residual]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(ImpulseError::NonFinite);
    }
    if work_residual.abs() > controls.work_tolerance {
        return Err(ImpulseError::Work {
            residual: work_residual,
            allowance: controls.work_tolerance,
        });
    }
    let acceleration = second_order_acceleration(
        model,
        &state.q,
        &velocity,
        forcing_after,
        mass_preconditioner,
        config,
        max_workspace_components - 8 * n,
        cancelled,
    )?;
    poll(cancelled)?;
    state.v = velocity;
    state.a = acceleration.value;
    Ok(SecondOrderImpulse {
        velocity_change: change,
        momentum_change: momentum,
        momentum_residual: residual,
        kinetic_energy_change: kinetic,
        supplied_work: work,
        work_residual,
        velocity_solve: velocity_solve.report,
        acceleration_solve: acceleration.report,
    })
}
