mod common;
use common::*;
use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::principal::PrincipalFrameError;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};

fn zero() -> Vec3 {
    Vec3::new(0.0, 0.0, 0.0)
}
fn ga(v: fs_mbd::Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}
fn mbd(v: Vec3) -> fs_mbd::Vec3 {
    fs_mbd::Vec3::new(v.x, v.y, v.z)
}
fn same(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn frame() -> Se3 {
    Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.4, -0.3, 0.7))).unwrap(),
        Vec3::new(3.0, 4.0, -2.0),
    )
    .unwrap()
}
fn body(density: f64) -> ResolvedBody {
    let card = density_card("source", density, density);
    let p = point(300.0);
    let placed = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.6, 0.2, -0.4))).unwrap(),
        Vec3::new(0.8, -0.3, 0.4),
    )
    .unwrap();
    let geometry = VolumeMoments::cuboid(Vec3::new(1.0, 2.0, 3.0))
        .unwrap()
        .transformed(placed)
        .unwrap();
    resolve_body(
        &[RegionBinding {
            key: "solid",
            geometry,
            material: &card,
            point: &p,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn compare_pose(a: Se3, b: Se3) {
    same(a.translation(), b.translation());
    for (a, b) in a.rotation().matrix().m.iter().zip(b.rotation().matrix().m) {
        close(*a, b);
    }
}

#[test]
fn full_material_tensor_reconstructs_in_the_actual_principal_pose() {
    let body = body(7.0);
    let principal = body.principal_body(1e-12, || false).unwrap();
    assert!(core::ptr::eq(principal.source(), &body));
    close(principal.properties().mass(), body.inertia().mass());
    let moments = principal.properties().principal_inertia_body();
    close(moments.x, 42.0 * 5.0 / 12.0);
    close(moments.y, 42.0 * 10.0 / 12.0);
    close(moments.z, 42.0 * 13.0 / 12.0);
    let r = principal.body_from_principal().rotation().matrix();
    let diag = fs_ga::Mat3 {
        m: [
            moments.x, 0.0, 0.0, 0.0, moments.y, 0.0, 0.0, 0.0, moments.z,
        ],
    };
    let reconstructed = r.compose(diag).compose(r.transpose());
    for (a, b) in reconstructed.m.iter().zip(body.inertia().inertia_com().m) {
        close(*a, b);
    }
    assert!(principal.relative_tensor_residual() <= 1e-12);
    same(
        principal.body_from_principal().translation(),
        body.inertia().center_of_mass(),
    );
}

#[test]
fn frame_conversion_preserves_geometry_point_velocity_momentum_and_kinetic_energy() {
    let body = body(7.0);
    let principal = body.principal_body(1e-12, || false).unwrap();
    let world = frame();
    let twist = Twist::new(Vec3::new(0.2, -0.4, 0.7), Vec3::new(1.2, 0.3, -0.2));
    let state = principal.state_from_kinematics(world, twist).unwrap();
    compare_pose(principal.world_from_body(state).unwrap(), world);
    let point = Vec3::new(-0.2, 0.6, 1.4);
    let arm = principal.arm_for_body_point(point).unwrap();
    let kinematics = state.point_kinematics(principal.properties(), arm).unwrap();
    same(
        ga(kinematics.point_world),
        world.rotation().rotate(point).unwrap() + world.translation(),
    );
    same(
        ga(kinematics.point_velocity_world),
        world
            .rotation()
            .rotate(twist.linear + twist.angular.cross(point))
            .unwrap(),
    );
    close(
        state.kinetic_energy(principal.properties()).unwrap(),
        body.inertia().kinetic_energy(twist).unwrap(),
    );
    let expected = body.inertia().momentum(twist).unwrap();
    let p_body = world
        .rotation()
        .inverse()
        .rotate(ga(state.linear_momentum_world()))
        .unwrap();
    let h_com = principal
        .body_from_principal()
        .rotation()
        .rotate(ga(state.angular_momentum_body()))
        .unwrap();
    // Spatial angular momentum is about the authored origin, not the centroid.
    same(p_body, expected.force);
    same(
        h_com + body.inertia().center_of_mass().cross(p_body),
        expected.torque,
    );
}

#[test]
fn material_swap_changes_actual_contact_effective_mass_and_impulse_response() {
    let point = Vec3::new(-0.2, 0.6, 1.4);
    let direction = Vec3::new(0.0, 1.0, 0.0);
    let impulse = direction.scale(2.0);
    let mut masses = Vec::new();
    let mut changes = Vec::new();
    for density in [7.0, 21.0] {
        let body = body(density);
        let principal = body.principal_body(1e-12, || false).unwrap();
        let state = principal
            .state_from_kinematics(frame(), Twist::zero())
            .unwrap();
        let arm = principal.arm_for_body_point(point).unwrap();
        let effective = state
            .directional_effective_mass(principal.properties(), arm, mbd(direction))
            .unwrap();
        let response = state
            .apply_impulse_at_body_point(principal.properties(), arm, mbd(impulse))
            .unwrap();
        let change = ga(response.point_kinematics_after.point_velocity_world)
            - ga(response.point_kinematics_before.point_velocity_world);
        close(
            direction.dot(change),
            2.0 * effective.inverse_mass_kg_inverse,
        );
        close(
            response.work.kinetic_energy_change_j,
            response.work.impulse_work_j,
        );
        assert!(response.work.work_energy_residual_j.abs() < 1e-12);
        compare_pose(
            principal.world_from_body(response.state_after).unwrap(),
            frame(),
        );
        masses.push(effective.effective_mass_kg);
        changes.push(change);
    }
    close(masses[1], 3.0 * masses[0]);
    same(changes[1].scale(3.0), changes[0]);
}

#[test]
fn moving_material_balance_does_not_move_authored_geometry_or_contact_points() {
    let a = density_card("a", 1.0, 1.0);
    let b = density_card("b", 3.0, 3.0);
    let p = point(300.0);
    let shape = VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let left = shape
        .transformed(Se3::from_parts(So3::identity(), Vec3::new(-2.0, -1.0, 0.0)).unwrap())
        .unwrap();
    let right = shape
        .transformed(Se3::from_parts(So3::identity(), Vec3::new(2.0, 1.0, 0.0)).unwrap())
        .unwrap();
    let mut centres = Vec::new();
    let point = Vec3::new(2.5, 1.5, 0.5);
    let twist = Twist::new(Vec3::new(0.1, 0.2, -0.3), Vec3::new(1.0, -2.0, 0.5));
    for (m_left, m_right) in [(&a, &b), (&b, &a)] {
        let body = resolve_body(
            &[
                RegionBinding {
                    key: "left",
                    geometry: left,
                    material: m_left,
                    point: &p,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
                RegionBinding {
                    key: "right",
                    geometry: right,
                    material: m_right,
                    point: &p,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
            ],
            2,
            || false,
        )
        .unwrap();
        let principal = body.principal_body(1e-12, || false).unwrap();
        let state = principal.state_from_kinematics(frame(), twist).unwrap();
        let point_state = state
            .point_kinematics(
                principal.properties(),
                principal.arm_for_body_point(point).unwrap(),
            )
            .unwrap();
        compare_pose(principal.world_from_body(state).unwrap(), frame());
        same(
            ga(point_state.point_world),
            frame().rotation().rotate(point).unwrap() + frame().translation(),
        );
        same(
            ga(point_state.point_velocity_world),
            frame()
                .rotation()
                .rotate(twist.linear + twist.angular.cross(point))
                .unwrap(),
        );
        centres.push(body.inertia().center_of_mass());
    }
    same(centres[0], centres[1].scale(-1.0));
    assert!(centres[0].x.abs() > 0.5);
}

#[test]
fn degenerate_axes_and_half_turns_preserve_the_physics_without_axis_uniqueness_claims() {
    let material = density_card("isotropic-inertia", 2.0, 2.0);
    let p = point(300.0);
    let shape = VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let body = resolve_body(
        &[RegionBinding {
            key: "cube",
            geometry: shape,
            material: &material,
            point: &p,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap();
    let principal = body.principal_body(1e-12, || false).unwrap();
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let world = Se3::from_parts(
            So3::exp(So3Tangent::new(axis.scale(core::f64::consts::PI))).unwrap(),
            zero(),
        )
        .unwrap();
        let state = principal
            .state_from_kinematics(world, Twist::new(axis, zero()))
            .unwrap();
        compare_pose(principal.world_from_body(state).unwrap(), world);
        close(
            state.kinetic_energy(principal.properties()).unwrap(),
            1.0 / 6.0,
        );
    }
}

#[test]
fn frame_budget_and_cancellation_refuse_without_modifying_the_source() {
    let body = body(7.0);
    let saved = body.clone();
    for tolerance in [0.0, -1.0, 1e-2, f64::NAN] {
        assert!(body.principal_body(tolerance, || false).is_err());
    }
    assert!(matches!(
        body.principal_body(1e-12, || true),
        Err(PrincipalFrameError::Cancelled)
    ));
    let mut count = 0;
    assert!(matches!(
        body.principal_body(1e-12, || {
            count += 1;
            count == 2
        }),
        Err(PrincipalFrameError::Cancelled)
    ));
    let principal = body.principal_body(1e-12, || false).unwrap();
    assert!(
        principal
            .state_from_kinematics(frame(), Twist::new(Vec3::new(f64::NAN, 0.0, 0.0), zero()))
            .is_err()
    );
    assert!(
        principal
            .arm_for_body_point(Vec3::new(0.0, f64::INFINITY, 0.0))
            .is_err()
    );
    assert_eq!(body, saved);
}

#[test]
fn configuration_transfer_preserves_world_momenta_through_new_centroid_and_axes() {
    use fs_mbd::articulated::SpatialInertia;
    use fs_mbd_material::principal::PrincipalFrame;
    use fs_mbd_material::principal::transfer::ConfigurationTransferControls;
    let material_body = body(7.0);
    let original = *material_body
        .principal_body(1e-12, || false)
        .unwrap()
        .frame();
    let rotated = So3::exp(So3Tangent::new(Vec3::new(-0.8, 0.3, 0.2)))
        .unwrap()
        .matrix();
    let tensor = rotated
        .compose(fs_ga::Mat3 {
            m: [90., 0., 0., 0., 70., 0., 0., 0., 60.],
        })
        .compose(rotated.transpose());
    let inertia = SpatialInertia::new(
        material_body.inertia().mass(),
        Vec3::new(-1., 0.8, 1.2),
        tensor,
    )
    .unwrap();
    let next = PrincipalFrame::from_inertia(inertia, 1e-12, || false).unwrap();
    let world = frame();
    let state = original
        .state_from_kinematics(
            world,
            Twist::new(Vec3::new(0.3, -0.5, 0.7), Vec3::new(1.2, -0.6, 0.4)),
        )
        .unwrap();
    let controls = ConfigurationTransferControls {
        relative_momentum_tolerance: 1e-12,
        maximum_parameter_work_j: 100.,
    };
    let transfer = next
        .transfer_state_from(&original, state, controls, || false)
        .unwrap();
    let after = transfer.state_after;
    assert_eq!(state.pose().position_world(), after.pose().position_world());
    assert_eq!(state.linear_momentum_world(), after.linear_momentum_world());
    let integrator = fs_mbd::RigidBodyIntegrator::new(fs_mbd::Gravity::ZERO);
    let old_d = integrator
        .diagnostics(state, original.properties())
        .unwrap();
    let new_d = integrator.diagnostics(after, next.properties()).unwrap();
    same(
        ga(old_d.angular_momentum_world),
        ga(new_d.angular_momentum_world),
    );
    let placed = next.world_from_body(after).unwrap();
    for (a, b) in placed
        .rotation()
        .matrix()
        .m
        .into_iter()
        .zip(world.rotation().matrix().m)
    {
        close(a, b);
    }
    same(
        placed.translation() + placed.rotation().rotate(inertia.center_of_mass()).unwrap(),
        ga(state.pose().position_world()),
    );
    close(
        transfer.parameter_work_j,
        new_d.mechanical_energy - old_d.mechanical_energy,
    );
    assert!(transfer.parameter_work_j.abs() > 0.01);
    assert!(transfer.relative_momentum_residual < 1e-12);
    let back = original
        .transfer_state_from(&next, after, controls, || false)
        .unwrap();
    compare_pose(original.world_from_body(back.state_after).unwrap(), world);
    same(
        ga(back.state_after.angular_momentum_body()),
        ga(state.angular_momentum_body()),
    );
    close(back.parameter_work_j, -transfer.parameter_work_j);
    assert_eq!(
        original
            .transfer_state_from(&original, state, controls, || false)
            .unwrap()
            .state_after,
        state
    );
    let zero_work = ConfigurationTransferControls {
        maximum_parameter_work_j: 0.,
        ..controls
    };
    assert!(
        next.transfer_state_from(&original, state, zero_work, || false)
            .is_err()
    );
    assert!(matches!(
        next.transfer_state_from(&original, state, controls, || true),
        Err(PrincipalFrameError::Cancelled)
    ));
    let other_body = body(14.0);
    assert!(
        other_body
            .principal_body(1e-12, || false)
            .unwrap()
            .frame()
            .transfer_state_from(&original, state, controls, || false)
            .is_err()
    );
}
