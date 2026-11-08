//! Original finite material-point work inside the simultaneous actuator owner.

use super::{
    Dd, Gravity, HeldExternalLoads, InternalWrenchStepError, LoadChange, MassProperties,
    RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench, add, cross, dd, midpoint,
    native_loaded_step, sub, world_vector,
};
use crate::spatial_loads::{
    PreparedSpatialLoads, basis_world,
    midpoint::{SpatialMidpointConfig, SpatialMidpointReceipt},
};

/// One original prepared spatial load per body, plus a uniform gravity field.
///
/// Each load retains world-held material-point/COM forces and body/world free
/// couples. Preparation owns the bounded point reduction; finite motion uses
/// its constant-sized moments without allocation or copying the point bank.
/// The midpoint limits apply to each nonzero prepared load. A zero load keeps
/// the existing native free-body chart and applies gravity once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeldSpatialLoads<'a> {
    pub bodies: &'a [PreparedSpatialLoads],
    pub gravity: Gravity,
    pub midpoint: SpatialMidpointConfig,
}

#[derive(Clone, Copy)]
pub(super) enum BodyLoads<'a> {
    Unloaded,
    Held(HeldExternalLoads<'a>),
    Spatial(HeldSpatialLoads<'a>),
}

impl BodyLoads<'_> {
    pub(super) fn check_count(self, bodies: usize) -> Result<(), InternalWrenchStepError> {
        let count = match self {
            Self::Unloaded => bodies,
            Self::Held(loads) => loads.wrenches.len(),
            Self::Spatial(loads) => loads.bodies.len(),
        };
        if bodies != count {
            return Err(InternalWrenchStepError::ExternalLoadCount {
                bodies,
                wrenches: count,
            });
        }
        Ok(())
    }

    pub(super) fn step(
        self,
        body: usize,
        before: RigidBodyState,
        properties: MassProperties,
        duration_s: f64,
        origin: Vec3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(StepReceipt, LoadChange), InternalWrenchStepError> {
        let (wrench, gravity) = match self {
            Self::Unloaded => (Wrench::ZERO, Gravity::ZERO),
            Self::Held(loads) => (loads.wrenches[body], loads.gravity),
            Self::Spatial(loads) => {
                let applied = &loads.bodies[body];
                if *applied != PreparedSpatialLoads::ZERO {
                    let receipt = RigidBodyIntegrator::new(loads.gravity).step_spatial_midpoint(
                        before,
                        properties,
                        applied,
                        duration_s,
                        loads.midpoint,
                        cancelled,
                    )?;
                    let source =
                        spatial_change(&receipt, properties, applied, loads.gravity, origin);
                    return Ok((receipt.motion.step, source));
                }
                (Wrench::ZERO, loads.gravity)
            }
        };
        let step = native_loaded_step(before, properties, wrench, gravity, duration_s)?;
        let source = LoadChange::from_step(&step, properties, wrench, gravity, origin);
        Ok((step, source))
    }
}

fn spatial_change(
    receipt: &SpatialMidpointReceipt,
    properties: MassProperties,
    loads: &PreparedSpatialLoads,
    gravity: Gravity,
    origin: Vec3,
) -> LoadChange {
    let step = &receipt.motion.step;
    let (before, after) = (step.state_before, step.state_after);
    let (duration, mass, half) = (
        Dd::from_f64(step.duration_seconds),
        Dd::from_f64(properties.mass()),
        Dd::from_f64(0.5),
    );
    let force = dd(loads.force_world);
    let acceleration = dd(gravity.acceleration_world());
    let linear = core::array::from_fn(|i| (force[i] + mass * acceleration[i]) * duration);
    let body_couple = add(
        world_vector(before, loads.torque_body),
        world_vector(after, loads.torque_body),
    )
    .map(|v| v * half);
    let mut moment = add(dd(loads.torque_world), body_couple);
    let axes0 = basis_world(before.pose().orientation());
    let axes1 = basis_world(after.pose().orientation());
    for j in 0..3 {
        // Sum (Rbar r) cross F from the ORIGINAL compensated first moments.
        // No instantaneous frozen arm or measured angular-momentum difference
        // is substituted for the physical external source quadrature.
        moment = add(
            moment,
            cross(
                midpoint(axes0[j], axes1[j]),
                loads.point_force_moments_world[j].0,
            ),
        );
    }
    let radius = sub(
        midpoint(
            before.pose().position_world(),
            after.pose().position_world(),
        ),
        dd(origin),
    );
    let angular = add(moment.map(|v| v * duration), cross(radius, linear));
    let force_work = Dd::from_f64(receipt.motion.force_work_j);
    let torque_work = Dd::from_f64(receipt.motion.couple_work_j);
    LoadChange {
        linear,
        angular,
        force_work,
        torque_work,
        applied_work: force_work + torque_work,
        gravity_work: Dd::from_f64(step.work.gravity_work_j),
    }
}
