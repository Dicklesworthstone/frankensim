//! Oriented finite support and body-attached drag on the original native owners.
use super::*;
use fs_couple::porous_momentum::free_body::{
    FreePorousBody, PorousBodySupport, PrincipalPorousResistance,
};

const SEMIAXES: [f64; 3] = [0.64, 0.46, 0.32];
fn support() -> PorousBodySupport {
    PorousBodySupport::Ellipsoid {
        semiaxes_m: SEMIAXES,
    }
}
fn law() -> PrincipalPorousResistance {
    PrincipalPorousResistance {
        darcy_principal: [0.06, 0.14, 0.27],
        forchheimer: 2.0,
    }
}
fn controls() -> PorousRotorControls {
    let mut c = super::controls();
    c.flow.maximum_native_trials = 256;
    c.flow.maximum_trial_cells = 256 * CELLS;
    c
}
fn fixture(shape: PorousBodySupport) -> FreePorousBody {
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.47).unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, 1.0, -1.0), 0.31).unwrap();
    let properties =
        MassProperties::from_principal_frame(16.0, Vec3::ZERO, Vec3::new(1.2, 1.6, 2.0), principal)
            .unwrap();
    let body = RigidBodyState::new(
        Pose::new(COM, orientation).unwrap(),
        Vec3::new(-0.006, 0.008, -0.004).scale(properties.mass()),
        properties
            .angular_momentum_body_checked(Vec3::new(0.015, -0.02, 0.06))
            .unwrap(),
    )
    .unwrap();
    FreePorousBody::new_with_support(
        grid(true),
        body,
        properties,
        units(),
        shape,
        CELLS,
        &CancelGate::new(),
    )
    .unwrap()
}
fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

/// Independent native Cayley half angle, rather than the implementation's
/// normalized quaternion half sum.
fn native_midpoint(
    before: RigidBodyState,
    after: RigidBodyState,
    properties: MassProperties,
) -> UnitQuaternion {
    let omega = properties
        .angular_velocity_body_checked(
            before
                .angular_momentum_body()
                .add(after.angular_momentum_body())
                .scale(0.5),
        )
        .unwrap();
    let speed = norm(omega);
    if speed == 0.0 {
        before.pose().orientation()
    } else {
        before
            .pose()
            .orientation()
            .right_exp(omega.scale((0.5 * units().timestep_s * speed).atan() / speed))
            .unwrap()
    }
}
fn in_ellipsoid(arm: Vec3, frame: UnitQuaternion) -> bool {
    let local = frame.rotate_world_to_body(arm);
    a(local)
        .into_iter()
        .zip(SEMIAXES)
        .map(|(r, a)| (r / a).powi(2))
        .sum::<f64>()
        <= 1.0
}
fn drag_in_world(slip: Vec3, frame: UnitQuaternion, law: PrincipalPorousResistance) -> Vec3 {
    let local = frame.rotate_world_to_body(slip);
    let linear = frame.rotate_body_to_world(Vec3::new(
        law.darcy_principal[0] * local.x,
        law.darcy_principal[1] * local.y,
        law.darcy_principal[2] * local.z,
    ));
    linear.add(slip.scale(law.forchheimer * norm(slip)))
}

#[test]
fn moving_ellipsoid_closes_anisotropic_native_force_torque_and_work() {
    let before = fixture(support());
    let populations = snapshot(before.fluid());
    let body_before = before.body();
    let force = Vec3::new(0.009, -0.004, 0.006);
    let couple = Vec3::new(0.004, -0.003, 0.002);
    let result = before
        .step_candidate_with_principal_resistance(
            law(),
            forcing,
            force,
            couple,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let body_after = result.motion.body();
    let frame = native_midpoint(body_before, body_after, before.properties());
    let centre = body_before
        .pose()
        .position_world()
        .add(body_after.pose().position_world())
        .scale(0.5);
    let u = units();
    let speed_unit = u.cell_size_m / u.timestep_s;
    let impulse_unit = u.mass_per_density_cell_kg * speed_unit;
    let angular_unit = impulse_unit * u.cell_size_m;
    let mut impulse = Vec3::ZERO;
    let mut torque = Vec3::ZERO;
    let mut loss = 0.0;
    let mut active = 0;
    let mut bounding_sphere_only = 0;
    let mut frozen_axis_defect = Vec3::ZERO;
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let arm = position([x, y, z]).sub(centre);
                let inside = in_ellipsoid(arm, frame);
                if !inside {
                    bounding_sphere_only += usize::from(norm(arm) < before.radius_m());
                    continue;
                }
                active += 1;
                let imposed = result
                    .midpoint_velocity_m_s
                    .add(result.angular_velocity_world_rad_s.cross(arm))
                    .scale(1.0 / speed_unit);
                let (rho0, p0) = cell(before.fluid(), x, y, z);
                let (rho1, p1) = cell(result.motion.fluid(), x, y, z);
                let slip0 = p0.scale(1.0 / rho0).sub(imposed);
                let slip1 = p1.scale(1.0 / rho1).sub(imposed);
                let j = drag_in_world(slip0, frame, law())
                    .add(drag_in_world(slip1, frame, law()))
                    .scale(0.5 * impulse_unit);
                let frozen = drag_in_world(slip0, body_before.pose().orientation(), law())
                    .add(drag_in_world(
                        slip1,
                        body_before.pose().orientation(),
                        law(),
                    ))
                    .scale(0.5 * impulse_unit);
                frozen_axis_defect = frozen_axis_defect.add(frozen.sub(j));
                impulse = impulse.add(j);
                torque = torque.add(arm.cross(j));
                loss += j.dot(slip0.add(slip1).scale(0.5 * speed_unit));
            }
        }
    }
    assert!(active > 4 && bounding_sphere_only > active);
    assert!(norm(frozen_axis_defect) > 1e-7);
    close_v(
        v(result.exchange.skeleton_impulse).scale(impulse_unit),
        impulse,
        4e-12,
    );
    close_v(
        v(result.exchange.skeleton_angular_impulse).scale(angular_unit),
        torque,
        4e-12,
    );
    close_v(
        body_after
            .linear_momentum_world()
            .sub(body_before.linear_momentum_world()),
        impulse.add(force.scale(u.timestep_s)),
        controls().flow.linear_tolerance_kg_m_s,
    );
    close_v(
        world_spin(body_after).sub(world_spin(body_before)),
        torque.add(couple.scale(u.timestep_s)),
        controls().flow.angular_tolerance_kg_m2_s,
    );
    close(result.relative_dissipation_j, loss, 4e-12);
    close(
        result.body_kinetic_change_j + result.fluid_drag_work_j + result.relative_dissipation_j
            - result.body_external_work_j,
        0.0,
        controls().flow.work_tolerance_j,
    );
    assert!(result.relative_dissipation_j > 0.0);
    assert!(result.absolute_linear_residual_kg_m_s <= controls().flow.linear_tolerance_kg_m_s);
    assert!(result.absolute_angular_residual_kg_m2_s <= controls().flow.angular_tolerance_kg_m2_s);
    assert_eq!(result.motion.support(), support());
    assert_eq!(result.motion.properties(), before.properties());
    assert_eq!(snapshot(before.fluid()), populations);
    assert_eq!(before.body(), body_before);
    assert_eq!(before.accepted_steps(), 0);
    assert_eq!(result.motion.accepted_steps(), 1);
}

#[test]
fn scalar_sphere_entry_and_principal_specialization_publish_identical_native_states() {
    let before = fixture(PorousBodySupport::Sphere { radius_m: 0.55 });
    let original = PorousRotorResistance {
        darcy: 0.12,
        forchheimer: 2.0,
    };
    let scalar = before
        .step_candidate(
            original,
            forcing,
            Vec3::ZERO,
            Vec3::ZERO,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let principal = before
        .step_candidate_with_principal_resistance(
            PrincipalPorousResistance::isotropic(original.darcy, original.forchheimer),
            forcing,
            Vec3::ZERO,
            Vec3::ZERO,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(scalar.mechanics, principal.mechanics);
    assert_eq!(scalar.exchange, principal.exchange);
    assert_eq!(
        snapshot(scalar.motion.fluid()),
        snapshot(principal.motion.fluid())
    );
    assert_eq!(scalar.motion.time_s(), principal.motion.time_s());
    let old_name: &fs_couple::porous_momentum::free_spherical::FreeSphericalPorousBody = &before;
    assert_eq!(old_name.radius_m(), 0.55);
}

#[test]
fn anisotropic_refusal_and_cancel_preserve_both_original_owners() {
    let before = fixture(support());
    let population = snapshot(before.fluid());
    let body = before.body();
    let bad = PrincipalPorousResistance {
        darcy_principal: [0.1, -0.2, 0.3],
        forchheimer: 0.0,
    };
    assert!(matches!(
        before.step_candidate_with_principal_resistance(
            bad,
            |_| panic!("invalid resistance must refuse before source sampling"),
            Vec3::ZERO,
            Vec3::ZERO,
            controls(),
            &CancelGate::new(),
        ),
        Err(PorousMomentumError::Invalid(_))
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        before.step_candidate_with_principal_resistance(
            law(),
            forcing,
            Vec3::ZERO,
            Vec3::ZERO,
            controls(),
            &gate,
        ),
        Err(PorousMomentumError::Cancelled)
    ));
    let mut strict = controls();
    strict.flow.maximum_native_trials = 2;
    strict.flow.maximum_trial_cells = 2 * CELLS;
    assert!(
        before
            .step_candidate_with_principal_resistance(
                law(),
                forcing,
                Vec3::ZERO,
                Vec3::ZERO,
                strict,
                &CancelGate::new(),
            )
            .is_err()
    );
    assert_eq!(snapshot(before.fluid()), population);
    assert_eq!(before.body(), body);
    assert_eq!(before.accepted_steps(), 0);
    assert_eq!(before.time_s(), 0.0);
}
