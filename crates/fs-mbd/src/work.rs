//! Work conjugate to the native rigid-body force and torque stages.

use super::{
    DynamicsError, Gravity, MassProperties, MomentumWork, RigidBodyState, Vec3, finite_derived,
};
use fs_math::dd::Dd;

/// Signed work and the mechanical-energy balance of a smooth rigid-body step.
///
/// Applied world force is paired with the actual accepted COM displacement.
/// Held-force work is therefore `F . (x1-x0)`; endpoint Verlet loads use the
/// trapezoidal displacement quadrature `mean(F0,F1) . (x1-x0)`. The latter is
/// a numerical quadrature, not an exact integral for an arbitrary force field.
/// Body torque uses the original integrator's conjugate angular increment:
/// `dt * mean(omega0,omega1)` for implicit midpoint, or the sum of the native
/// torque-kick impulse/mean-velocity pairings for split methods. Each kick's
/// torque and velocity use the same body frame, including endpoint loads.
///
/// Mass, body inertia and uniform gravity are fixed over this receipt's
/// interval. Changes of those parameters between receipts require separate
/// parameter-work or momentum-transfer accounting; this balance does not
/// include the potential change from changing gravity at a fixed position.
///
/// Storage increments use polarized compensated arithmetic, independently of
/// the source work. The residual consequently retains free-rotation splitting
/// error, force quadrature error and loss from rounding the returned state.
/// None of these is physical friction, heat or actuator work. In particular,
/// kick-drift's force work is measured over its displacement; it is not replaced
/// by its kinetic-energy increment to manufacture an exact balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepWorkDiagnostics {
    /// Work done by the prescribed force on the COM [J]; negative for braking.
    pub applied_force_work_j: f64,
    /// Work done by prescribed body torques [J]; negative for regeneration.
    pub applied_torque_work_j: f64,
    /// Signed applied force plus applied torque work [J], excluding gravity.
    pub applied_work_j: f64,
    /// Conservative gravitational force work over the accepted displacement [J].
    pub gravity_work_j: f64,
    /// Actual kinetic change between the represented endpoint momenta [J].
    pub kinetic_energy_change_j: f64,
    /// Gravity potential change [J], equal to minus gravitational work.
    pub gravitational_potential_energy_change_j: f64,
    /// Actual kinetic plus gravity potential change [J].
    pub mechanical_energy_change_j: f64,
    /// `delta(K+V) - applied_work` [J], without subtracting rounded totals.
    /// Positive means numerical energy creation; negative means numerical loss.
    /// No nonnegative clipping or tolerance based on stored energy is applied.
    pub work_energy_residual_j: f64,
}

pub(super) fn angular_torque_work(
    properties: MassProperties,
    before: [Dd; 3],
    after: [Dd; 3],
    torque: Vec3,
    duration: f64,
    weight: f64,
) -> MomentumWork {
    let mut work = MomentumWork::default();
    if torque != Vec3::ZERO {
        work.add_angular_work(
            components(properties.principal_inertia_body),
            properties.principal_to_body.components(),
            core::array::from_fn(|i| {
                if before[i] == after[i] {
                    before[i]
                } else {
                    Dd::from_f64(0.5) * before[i] + Dd::from_f64(0.5) * after[i]
                }
            }),
            components(torque).map(Dd::from_f64),
            duration,
            weight,
        );
    }
    work
}

pub(super) fn finish(
    before: RigidBodyState,
    after: RigidBodyState,
    properties: MassProperties,
    gravity: Gravity,
    force_before: Vec3,
    force_after: Vec3,
    torque_work: MomentumWork,
) -> Result<StepWorkDiagnostics, DynamicsError> {
    // The step owner has already admitted the states and properties. Reuse the
    // fixed-inertia storage form; work never comes from subtracting its totals.
    let mut kinetic_change = MomentumWork::default();
    kinetic_change.add_fixed_inertia(
        properties.mass,
        components(properties.principal_inertia_body),
        properties.principal_to_body.components(),
        [
            components(before.linear_momentum_world),
            components(before.angular_momentum_body),
        ],
        [
            components(after.linear_momentum_world),
            components(after.angular_momentum_body),
        ],
    );
    let mut force_work = MomentumWork::default();
    let mut gravity_work = MomentumWork::default();
    if force_before != Vec3::ZERO || force_after != Vec3::ZERO || gravity != Gravity::ZERO {
        let old = components(before.pose.position_world);
        let new = components(after.pose.position_world);
        let displacement = core::array::from_fn(|i| Dd::from_f64(new[i]) - Dd::from_f64(old[i]));
        let first = components(force_before);
        let last = components(force_after);
        for i in 0..3 {
            if first[i] == last[i] {
                force_work.add_product(Dd::from_f64(first[i]), displacement[i]);
            } else {
                force_work.add_scaled_product(Dd::from_f64(first[i]), displacement[i], 0.5);
                force_work.add_scaled_product(Dd::from_f64(last[i]), displacement[i], 0.5);
            }
        }
        gravity_work.add_compensated(
            components(gravity.acceleration_world)
                .map(|g| Dd::from_f64(properties.mass) * Dd::from_f64(g)),
            displacement,
        );
    }
    let mut applied = force_work;
    applied.merge(torque_work);
    let mut mechanical_change = kinetic_change;
    mechanical_change.subtract(gravity_work);
    let residual = mechanical_change.difference(applied);
    let gravity_work_j = finite_derived(gravity_work.total().to_f64(), "step_work.gravity_j")?;
    Ok(StepWorkDiagnostics {
        applied_force_work_j: finite_derived(force_work.total().to_f64(), "step_work.force_j")?,
        applied_torque_work_j: finite_derived(torque_work.total().to_f64(), "step_work.torque_j")?,
        applied_work_j: finite_derived(applied.total().to_f64(), "step_work.applied_j")?,
        gravity_work_j,
        kinetic_energy_change_j: finite_derived(
            kinetic_change.total().to_f64(),
            "step_work.kinetic_change_j",
        )?,
        gravitational_potential_energy_change_j: -gravity_work_j,
        mechanical_energy_change_j: finite_derived(
            mechanical_change.total().to_f64(),
            "step_work.mechanical_change_j",
        )?,
        work_energy_residual_j: finite_derived(residual.to_f64(), "step_work.residual_j")?,
    })
}

fn components(vector: Vec3) -> [f64; 3] {
    [vector.x, vector.y, vector.z]
}
