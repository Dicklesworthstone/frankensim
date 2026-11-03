//! Shared scalar constraint step for two original finite rigid bodies.
//!
//! Geometry supplies the Jacobian direction and its actual material arms. The
//! numerical owner is the original distance step: initial physical reaction,
//! bounded residual-decreasing native free trials, then endpoint reaction.

use super::distance::{DistanceStepControls, DistanceStepError, DistanceStepReceipt};
use super::free_motion::{FreeFlow, SpatialFlow};
use crate::spatial_loads::{PreparedSpatialLoads, SpatialStepReceipt};
use crate::{
    ActionReactionImpulseReceipt, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3,
    Wrench, apply_equal_and_opposite_impulse_at_body_points,
};

#[derive(Clone, Copy)]
pub(super) struct Geometry {
    /// An impulse `j*normal` on body 1 increases the scalar gap.
    pub normal: Vec3,
    pub arms_body_m: [Vec3; 2],
    pub gap: f64,
    pub rate: f64,
}

pub(super) trait Constraint {
    fn properties(&self) -> [MassProperties; 2];
    fn tolerances(&self) -> [f64; 2];
    fn geometry(&self, states: [RigidBodyState; 2]) -> Result<Geometry, DistanceStepError>;
}

pub(super) fn advance(
    constraint: &impl Constraint,
    states: [RigidBodyState; 2],
    integrator: RigidBodyIntegrator,
    wrenches: [Wrench; 2],
    duration_s: f64,
    controls: DistanceStepControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DistanceStepReceipt, DistanceStepError> {
    advance_with(
        constraint,
        states,
        integrator,
        &wrenches[..],
        duration_s,
        controls,
        cancelled,
    )
}

pub(super) fn advance_spatial(
    constraint: &impl Constraint,
    states: [RigidBodyState; 2],
    integrator: RigidBodyIntegrator,
    loads: [PreparedSpatialLoads; 2],
    duration_s: f64,
    controls: DistanceStepControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DistanceStepReceipt<SpatialStepReceipt>, DistanceStepError> {
    poll(&mut cancelled)?;
    let bounds = [loads[0].torque_norm_bound()?, loads[1].torque_norm_bound()?];
    advance_with(
        constraint,
        states,
        integrator,
        SpatialFlow {
            loads: &loads,
            bounds: &bounds,
        },
        duration_s,
        controls,
        cancelled,
    )
}

fn advance_with<F: FreeFlow>(
    constraint: &impl Constraint,
    states: [RigidBodyState; 2],
    integrator: RigidBodyIntegrator,
    flow: F,
    duration_s: f64,
    controls: DistanceStepControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DistanceStepReceipt<F::Receipt>, DistanceStepError> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || !positive(controls.maximum_impulse_ns)
        || !positive(controls.maximum_rotation_increment_rad)
        || controls.maximum_rotation_increment_rad > 0.25
        || controls.maximum_position_iterations == 0
        || controls.maximum_backtracks == 0
    {
        return Err(DistanceStepError::Invalid(
            "positive finite local controls required",
        ));
    }
    let [position_tolerance, velocity_tolerance] = constraint.tolerances();
    let properties = constraint.properties();
    let initial = constraint.geometry(states)?;
    if initial.gap.abs() > position_tolerance || initial.rate.abs() > velocity_tolerance {
        return Err(DistanceStepError::InitialConstraint);
    }
    let response = duration_s * mobility(properties, states, initial)?;
    if !positive(response) {
        return Err(DistanceStepError::Invalid("distance impulse response"));
    }
    let problem = Problem {
        constraint,
        properties,
        states,
        integrator,
        flow,
        duration_s,
        controls,
        initial,
    };
    let mut impulse = 0.0;
    let mut trial = problem.evaluate(impulse, &mut cancelled)?;
    let mut iterations = 0;
    while trial.geometry.gap.abs() > position_tolerance {
        if iterations == controls.maximum_position_iterations {
            return Err(DistanceStepError::NotConverged);
        }
        let mut correction = -trial.geometry.gap / response;
        if !correction.is_finite() {
            return Err(DistanceStepError::LocalBranch);
        }
        let mut accepted = None;
        for _ in 0..controls.maximum_backtracks {
            poll(&mut cancelled)?;
            let next_impulse = impulse + correction;
            if next_impulse == impulse {
                break;
            }
            if next_impulse.abs() <= controls.maximum_impulse_ns {
                match problem.evaluate(next_impulse, &mut cancelled) {
                    Ok(candidate) if candidate.geometry.gap.abs() < trial.geometry.gap.abs() => {
                        accepted = Some((next_impulse, candidate));
                        break;
                    }
                    Ok(_) | Err(DistanceStepError::LocalBranch) => {}
                    Err(error) => return Err(error),
                }
            }
            correction *= 0.5;
        }
        let Some((next_impulse, candidate)) = accepted else {
            return Err(DistanceStepError::LocalBranch);
        };
        impulse = next_impulse;
        trial = candidate;
        iterations += 1;
    }
    poll(&mut cancelled)?;
    let free = trial.steps.map(F::state_after);
    let final_impulse = -trial.geometry.rate / mobility(properties, free, trial.geometry)?;
    if !final_impulse.is_finite() || final_impulse.abs() > controls.maximum_impulse_ns {
        return Err(DistanceStepError::LocalBranch);
    }
    let velocity_impulse = apply(properties, free, trial.geometry, final_impulse)?;
    let states_after = [
        velocity_impulse.first.state_after,
        velocity_impulse.second.state_after,
    ];
    let final_geometry = constraint.geometry(states_after)?;
    if final_geometry.gap.abs() > position_tolerance
        || final_geometry.rate.abs() > velocity_tolerance
    {
        return Err(DistanceStepError::NotConverged);
    }
    let constraint_work_j = trial.impulse.impulse_work_j + velocity_impulse.impulse_work_j;
    if !constraint_work_j.is_finite() {
        return Err(DistanceStepError::Invalid("distance constraint work"));
    }
    poll(&mut cancelled)?;
    Ok(DistanceStepReceipt {
        states_after,
        position_impulse: trial.impulse,
        free_steps: trial.steps,
        velocity_impulse,
        constraint_work_j,
        final_gap_m: final_geometry.gap,
        final_range_rate_mps: final_geometry.rate,
        position_iterations: iterations,
    })
}

fn mobility(
    properties: [MassProperties; 2],
    states: [RigidBodyState; 2],
    g: Geometry,
) -> Result<f64, DistanceStepError> {
    let a = states[0].directional_effective_mass(properties[0], g.arms_body_m[0], g.normal)?;
    let b = states[1].directional_effective_mass(properties[1], g.arms_body_m[1], g.normal)?;
    let result = a.inverse_mass_kg_inverse + b.inverse_mass_kg_inverse;
    if !positive(result) {
        return Err(DistanceStepError::Invalid("pair mobility"));
    }
    Ok(result)
}

fn apply(
    properties: [MassProperties; 2],
    states: [RigidBodyState; 2],
    g: Geometry,
    impulse: f64,
) -> Result<ActionReactionImpulseReceipt, DistanceStepError> {
    Ok(apply_equal_and_opposite_impulse_at_body_points(
        states[0],
        properties[0],
        g.arms_body_m[0],
        states[1],
        properties[1],
        g.arms_body_m[1],
        g.normal.scale(impulse),
    )?)
}

struct Trial<R> {
    impulse: ActionReactionImpulseReceipt,
    steps: [R; 2],
    geometry: Geometry,
}
struct Problem<'a, C, F> {
    constraint: &'a C,
    properties: [MassProperties; 2],
    states: [RigidBodyState; 2],
    integrator: RigidBodyIntegrator,
    flow: F,
    duration_s: f64,
    controls: DistanceStepControls,
    initial: Geometry,
}
impl<C: Constraint, F: FreeFlow> Problem<'_, C, F> {
    fn evaluate(
        &self,
        impulse_ns: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial<F::Receipt>, DistanceStepError> {
        poll(cancelled)?;
        let impulse = apply(self.properties, self.states, self.initial, impulse_ns)?;
        let states = [impulse.first.state_after, impulse.second.state_after];
        for i in 0..2 {
            let principal = self.properties[i].principal_inertia_body();
            let minimum = principal.x.min(principal.y).min(principal.z);
            let travel = (norm(states[i].angular_momentum_body())
                + self.duration_s * self.flow.torque_bound(i))
                * (self.duration_s / minimum);
            if !travel.is_finite() || travel > self.controls.maximum_rotation_increment_rad {
                return Err(DistanceStepError::LocalBranch);
            }
        }
        let first = self.flow.step(
            0,
            self.integrator,
            states[0],
            self.properties[0],
            self.duration_s,
        )?;
        poll(cancelled)?;
        let second = self.flow.step(
            1,
            self.integrator,
            states[1],
            self.properties[1],
            self.duration_s,
        )?;
        poll(cancelled)?;
        let geometry = self
            .constraint
            .geometry([F::state_after(first), F::state_after(second)])?;
        let max_chord = 2.0 * fs_math::det::sin(0.5 * self.controls.maximum_rotation_increment_rad);
        if norm(geometry.normal.sub(self.initial.normal)) > max_chord {
            return Err(DistanceStepError::LocalBranch);
        }
        Ok(Trial {
            impulse,
            steps: [first, second],
            geometry,
        })
    }
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DistanceStepError> {
    if cancelled() {
        Err(DistanceStepError::Cancelled)
    } else {
        Ok(())
    }
}
