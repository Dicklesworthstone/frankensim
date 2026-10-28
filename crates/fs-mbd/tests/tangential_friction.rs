use fs_contact::tangential::nonsmooth::{
    NonsmoothContact, NonsmoothControls, NonsmoothTangentialState,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn close_vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn controls() -> RigidFrictionControls {
    RigidFrictionControls {
        maximum_bodies: 8,
        maximum_contacts: 8,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 1e-8,
        contact: NonsmoothControls {
            maximum_evaluations: 4000,
            maximum_newton_iterations: 40,
            maximum_mechanical_dofs: 48,
            difference_step_mps: 1e-5,
            velocity_tolerance_mps: 1e-9,
            cone_tolerance: 1e-10,
            energy_tolerance_j: 1e-8,
            newton: Default::default(),
        },
    }
}
fn law() -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "synthetic/rigid-friction",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            tangential_stiffness_n_per_m: 10000.0,
            torsional_stiffness_nm_per_rad: 100.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .unwrap()
}
fn patch(id: &str, normal: f64) -> NormalPatchView {
    NormalPatchView::new(
        id,
        "synthetic/normal",
        "synthetic/rigid-contact",
        NormalPatchAuthority::SyntheticFixture,
        normal,
        0.02,
        0.01,
        1e-4,
    )
    .unwrap()
}
fn interface() -> PartialSlipInterface {
    PartialSlipInterface::new(
        "body-a->body-b",
        "history",
        "synthetic/interface",
        NormalPatchAuthority::SyntheticFixture,
    )
    .unwrap()
}
fn contact<'a>(
    law: &'a PartialSlipLaw,
    patch: &'a NormalPatchView,
    interface: &'a PartialSlipInterface,
) -> NonsmoothContact<'a> {
    NonsmoothContact {
        law,
        patch,
        interface,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        common_spin_rad_per_s: 0.0,
        torsional_length_m: 0.01,
        creepage_reference_speed_mps: 1.0,
        work: GeneralizedWorkOwnership::new(patch.patch_id(), "step-1", "x", "y", "spin").unwrap(),
    }
}
fn properties(mass: f64, inertia: Vec3) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, inertia).unwrap()
}
fn state(pose: Pose, momentum: Vec3, angular: Vec3) -> RigidBodyState {
    RigidBodyState::new(pose, momentum, angular).unwrap()
}

#[test]
fn rotated_anisotropic_bodies_use_point_arms_and_torsional_inertia() {
    let l = law();
    let p = patch("patch", 10000.0);
    let i = interface();
    let gate = CancelGate::new();
    let qa =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_2)
            .unwrap();
    let qb =
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), core::f64::consts::FRAC_PI_2)
            .unwrap();
    let poses = [
        Pose::new(Vec3::ZERO, qa).unwrap(),
        Pose::new(Vec3::new(0.0, 2.0, 0.0), qb).unwrap(),
    ];
    let props = [
        properties(2.0, Vec3::new(1.0, 1.5, 2.0)),
        properties(3.0, Vec3::new(0.7, 1.1, 1.4)),
    ];
    let mut bodies = [
        state(poses[0], Vec3::new(4.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)),
        state(
            poses[1],
            Vec3::new(-1.5, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.7),
        ),
    ];
    let before = bodies;
    let binding = RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: Vec3::new(1.0, 0.0, 0.0),
        arm_b_body: Vec3::new(0.0, 0.0, 1.0),
    };
    let mut c = contact(&l, &p, &i);
    c.frame = TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
    let mut history = NonsmoothTangentialState::new(1, 8).unwrap();
    let dt = 0.005;
    let result = apply_rigid_friction(
        &mut bodies,
        &props,
        &mut history,
        &[c],
        &[binding],
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    // Independent aligned-principal-axis solution. Force along world x has
    // unit COM arms and inertias I_Az=2, I_By=1.1. Spin about world y has
    // inertias I_Ax=1, I_Bz=1.4; pure torsion contributes no linear impulse.
    let inverse_mass = 0.5 + 1.0 / 3.0 + 1.0 / 2.0 + 1.0 / 1.1;
    let slip = 2.5 / (1.0 + dt * dt * 10000.0 * inverse_mass);
    let spin = 2.5 / (1.0 + dt * dt * 100.0 * (1.0 + 1.0 / 1.4));
    let impulse = -dt * dt * 10000.0 * slip;
    let couple = -dt * dt * 100.0 * spin;
    close(result.contacts[0].endpoint_velocity_mps[0], slip);
    close(result.contacts[0].endpoint_velocity_mps[2], 0.01 * spin);
    close(
        bodies[0].linear_momentum_world().x,
        before[0].linear_momentum_world().x + impulse,
    );
    close(
        bodies[1].linear_momentum_world().x,
        before[1].linear_momentum_world().x - impulse,
    );
    close(bodies[0].angular_momentum_body().x, 2.0 + couple);
    close(bodies[0].angular_momentum_body().z, -impulse);
    assert_eq!(bodies[0].pose(), before[0].pose());
    assert_eq!(bodies[1].pose(), before[1].pose());
    assert!(result.kinetic_after_j < result.kinetic_before_j);
    close(result.friction.physical_heat_j, 0.0);
    close(result.total_energy_residual_j, 0.0);
    close(result.work_energy_residual_j, 0.0);
    close_vector(result.linear_momentum_change_ns, Vec3::ZERO);
    close_vector(result.angular_momentum_change_nms, Vec3::ZERO);
    close(history.time_s(), dt);
}

#[test]
fn shared_body_contacts_publish_physical_recoil_independently_of_input_order() {
    let l = law();
    let pa = patch("a", 100.0);
    let pb = patch("b", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let props = [
        properties(2.0, Vec3::new(1.0, 1.0, 1.0)),
        properties(3.0, Vec3::new(1.0, 1.0, 1.0)),
        properties(4.0, Vec3::new(1.0, 1.0, 1.0)),
    ];
    let initial = [
        state(Pose::identity(), Vec3::new(100.0, 0.0, 0.0), Vec3::ZERO),
        state(Pose::identity(), Vec3::ZERO, Vec3::ZERO),
        state(Pose::identity(), Vec3::ZERO, Vec3::ZERO),
    ];
    let a = RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: Vec3::ZERO,
        arm_b_body: Vec3::ZERO,
    };
    let b = RigidContactBinding {
        body_a: 0,
        body_b: 2,
        ..a
    };
    let mut first = initial;
    let mut reversed = initial;
    let mut h1 = NonsmoothTangentialState::new(2, 8).unwrap();
    let mut h2 = h1.clone();
    let r1 = apply_rigid_friction(
        &mut first,
        &props,
        &mut h1,
        &[contact(&l, &pa, &i), contact(&l, &pb, &i)],
        &[a, b],
        0.001,
        controls(),
        &gate,
    )
    .unwrap();
    let r2 = apply_rigid_friction(
        &mut reversed,
        &props,
        &mut h2,
        &[contact(&l, &pb, &i), contact(&l, &pa, &i)],
        &[b, a],
        0.001,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(first, reversed);
    assert_eq!(h1, h2);
    assert_eq!(r1.contacts, r2.contacts);
    close(first[0].linear_momentum_world().x, 100.0 - 0.08);
    close(first[1].linear_momentum_world().x, 0.04);
    close(first[2].linear_momentum_world().x, 0.04);
    close(r1.contacts[0].endpoint_velocity_mps[0], 49.96 - 0.04 / 3.0);
    close(r1.contacts[1].endpoint_velocity_mps[0], 49.96 - 0.04 / 4.0);
    assert!(r1.friction.physical_heat_j > 0.0);
    close_vector(r1.linear_momentum_change_ns, Vec3::ZERO);
    close_vector(r1.angular_momentum_change_nms, Vec3::ZERO);
    close(r1.total_energy_residual_j, 0.0);
    let accepted = first;
    let accepted_history = h1.clone();
    assert!(
        apply_rigid_friction(
            &mut first,
            &props,
            &mut h1,
            &[contact(&l, &pa, &i), contact(&l, &pb, &i)],
            &[a, b],
            0.001,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(first, accepted);
    assert_eq!(h1, accepted_history);
}

#[test]
fn geometry_momentum_and_cancellation_refusals_are_atomic() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let props = [
        properties(2.0, Vec3::new(1.0, 1.0, 1.0)),
        properties(3.0, Vec3::new(1.0, 1.0, 1.0)),
    ];
    let initial = [
        state(Pose::identity(), Vec3::new(100.0, 0.0, 0.0), Vec3::ZERO),
        state(
            Pose::new(Vec3::new(0.0, 1e-6, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        ),
    ];
    let binding = RigidContactBinding {
        body_a: 0,
        body_b: 1,
        arm_a_body: Vec3::ZERO,
        arm_b_body: Vec3::ZERO,
    };
    let mut bodies = initial;
    let mut history = NonsmoothTangentialState::new(1, 8).unwrap();
    let original_history = history.clone();
    let gate = CancelGate::new();
    assert!(matches!(
        apply_rigid_friction(
            &mut bodies,
            &props,
            &mut history,
            &[contact(&l, &p, &i)],
            &[binding],
            0.001,
            controls(),
            &gate
        ),
        Err(RigidFrictionError::Invalid(_))
    ));
    let mut loose_geometry = controls();
    loose_geometry.geometry_tolerance_m = 1e-5;
    loose_geometry.angular_momentum_tolerance_nms = 1e-14;
    // Proximity admission cannot excuse an actual angular-momentum defect:
    // applying equal/opposite impulses at distinct points leaves a couple.
    assert!(matches!(
        apply_rigid_friction(
            &mut bodies,
            &props,
            &mut history,
            &[contact(&l, &p, &i)],
            &[binding],
            0.001,
            loose_geometry,
            &gate
        ),
        Err(RigidFrictionError::Audit(_))
    ));
    assert_eq!(bodies, initial);
    assert_eq!(history, original_history);
    gate.request();
    assert!(matches!(
        apply_rigid_friction(
            &mut bodies,
            &props,
            &mut history,
            &[contact(&l, &p, &i)],
            &[binding],
            0.001,
            loose_geometry,
            &gate
        ),
        Err(RigidFrictionError::Cancelled)
    ));
    assert_eq!(bodies, initial);
    assert_eq!(history, original_history);
}
