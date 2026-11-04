//! Independent exterior Stokes-law and actual native-motion regressions.
//! Synthetic solid cards and declared fluid values are implementation fixtures,
//! not measured properties of named materials.
mod common;

use common::{density_card, point};
use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3 as GaVec3};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::{Gravity, RigidBodyIntegrator, Vec3, Wrench};
use fs_mbd_material::creeping_flow::{CreepingFlowLimits, CreepingSphere, CreepingStepControls};
use fs_mbd_material::ellipsoid_flow::{
    CreepingEllipsoid, EllipsoidBoundary, EllipsoidFlowError, EllipsoidQuadrature,
    EllipsoidStepControls,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};

fn quadrature() -> EllipsoidQuadrature {
    EllipsoidQuadrature {
        initial_panels: 2,
        maximum_panels: 512,
        relative_tolerance: 1e-11,
    }
}
fn limits() -> CreepingFlowLimits {
    CreepingFlowLimits {
        maximum_reynolds: 0.05,
        maximum_diffusion_ratio: 0.05,
    }
}
fn controls(dt: f64) -> EllipsoidStepControls {
    EllipsoidStepControls {
        mechanical: CreepingStepControls {
            duration_s: dt,
            maximum_relaxation_fraction: 0.5,
            maximum_energy_defect_j: 1e-12,
        },
        maximum_coupling_iterations: 64,
        relative_wrench_tolerance: 1e-12,
        maximum_rotation_increment_rad: 0.2,
    }
}
fn shape(a: f64, b: f64, c: f64) -> EllipsoidBoundary {
    EllipsoidBoundary {
        semiaxes_m: GaVec3::new(a, b, c),
        body_from_ellipsoid: Se3::identity(),
    }
}
fn body(boundary: EllipsoidBoundary, density: f64) -> ResolvedBody {
    let card = density_card("synthetic ellipsoid", density, density);
    resolve_body(
        &[RegionBinding {
            key: "ellipsoid",
            geometry: boundary.volume_moments().unwrap(),
            material: &card,
            point: &point(300.0),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn fluid(rho: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap()
}
fn rotation(x: f64, y: f64, z: f64) -> So3 {
    So3::exp(So3Tangent::new(GaVec3::new(x, y, z))).unwrap()
}
fn mv(v: GaVec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}
fn gv(v: Vec3) -> GaVec3 {
    GaVec3::new(v.x, v.y, v.z)
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e} at relative tolerance {tolerance:e}"
    );
}
fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        magnitude(actual.sub(expected)) <= tolerance,
        "{actual:?} != {expected:?} at absolute tolerance {tolerance:e}"
    );
}

// Closed forms B6/B7, Anand & Narsimhan, JFM 983 (2024), A28,
// DOI 10.1017/jfm.2024.136. These do not use the product quadrature.
// The symmetry semiaxis a is x, with equal transverse semiaxes b.
fn spheroid(a: f64, b: f64, mu: f64) -> ([f64; 3], [f64; 3]) {
    let (xa, ya, xc, yc, length) = if a > b {
        let e = (1.0 - (b / a).powi(2)).sqrt();
        let e2 = e * e;
        let e3 = e2 * e;
        let logarithm = ((1.0 + e) / (1.0 - e)).ln();
        (
            (8.0 / 3.0) * e3 / (-2.0 * e + (1.0 + e2) * logarithm),
            (16.0 / 3.0) * e3 / (2.0 * e + (3.0 * e2 - 1.0) * logarithm),
            (4.0 / 3.0) * e3 * (1.0 - e2) / (2.0 * e - (1.0 - e2) * logarithm),
            (4.0 / 3.0) * e3 * (2.0 - e2) / (-2.0 * e + (1.0 + e2) * logarithm),
            a,
        )
    } else {
        let root = a / b;
        let e = (1.0 - root * root).sqrt();
        let e2 = e * e;
        let e3 = e2 * e;
        let angle = e.atan2(root);
        (
            (4.0 / 3.0) * e3 / ((2.0 * e2 - 1.0) * angle + e * root),
            (8.0 / 3.0) * e3 / ((2.0 * e2 + 1.0) * angle - e * root),
            (2.0 / 3.0) * e3 / (angle - e * root),
            (2.0 / 3.0) * e3 * (2.0 - e2) / (e * root - (1.0 - 2.0 * e2) * angle),
            b,
        )
    };
    let c = 6.0 * core::f64::consts::PI * mu * length;
    let d = 8.0 * core::f64::consts::PI * mu * length.powi(3);
    ([c * xa, c * ya, c * ya], [d * xc, d * yc, d * yc])
}

#[test]
fn sphere_reduction_matches_the_original_driven_native_owner() {
    let radius: f64 = 0.01;
    let boundary = shape(radius, radius, radius);
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = Vec3::new(0.01, -0.02, 0.015);
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::ZERO,
        ambient,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let sphere = CreepingSphere::new(&material, radius, Gravity::ZERO, ambient, limits()).unwrap();
    assert!(core::ptr::eq(medium.source(), &material));
    let coefficients = medium.resistance();
    assert_eq!(coefficients.quadrature_panels, 0);
    for axis in 0..3 {
        relative(
            coefficients.translation_kg_s[axis],
            6.0 * core::f64::consts::PI * radius,
            2e-15,
        );
        relative(
            coefficients.rotation_kg_m2_s[axis],
            8.0 * core::f64::consts::PI * radius.powi(3),
            2e-15,
        );
    }
    let state = native
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(
                GaVec3::new(0.03, -0.02, 0.01),
                GaVec3::new(0.02, -0.03, 0.04),
            ),
        )
        .unwrap();
    let applied = Wrench {
        force_world: Vec3::new(0.005, -0.002, 0.001),
        torque_body: Vec3::new(1e-5, -2e-6, 3e-6),
    };
    let candidate = medium
        .step_with_wrench(&native, state, applied, controls(0.01), || false)
        .unwrap()
        .flow;
    let reference = sphere
        .step_with_wrench(&native, state, applied, controls(0.01).mechanical, || false)
        .unwrap();
    vector_close(
        candidate.dynamics.state_after.linear_momentum_world(),
        reference.dynamics.state_after.linear_momentum_world(),
        1e-15,
    );
    vector_close(
        candidate.dynamics.state_after.angular_momentum_body(),
        reference.dynamics.state_after.angular_momentum_body(),
        1e-15,
    );
    vector_close(
        candidate.dynamics.state_after.pose().position_world(),
        reference.dynamics.state_after.pose().position_world(),
        1e-15,
    );
    relative(candidate.viscous_heat_j, reference.viscous_heat_j, 2e-12);
    relative(candidate.applied_work_j, reference.applied_work_j, 2e-12);
    assert!(candidate.energy_defect_j.abs() < 1e-13);
}

#[test]
fn all_six_coefficients_match_independent_prolate_and_oblate_closed_forms() {
    let material = fluid(1.0, 0.4);
    for (a, b) in [(0.03, 0.01), (0.005, 0.02), (0.11, 0.0011), (0.0002, 0.02)] {
        let medium = CreepingEllipsoid::new(
            &material,
            shape(a, b, b),
            Gravity::ZERO,
            Vec3::ZERO,
            limits(),
            quadrature(),
            || false,
        )
        .unwrap();
        let actual = medium.resistance();
        let (translation, rotation) = spheroid(a, b, 0.4);
        for axis in 0..3 {
            relative(actual.translation_kg_s[axis], translation[axis], 2e-9);
            relative(actual.rotation_kg_m2_s[axis], rotation[axis], 2e-9);
        }
        assert!(actual.quadrature_panels >= 4);
        assert!(actual.estimated_relative_change <= quadrature().relative_tolerance);
    }
}

#[test]
fn triaxial_force_torque_and_heat_use_authored_axes_after_native_principal_sorting() {
    let boundary = EllipsoidBoundary {
        semiaxes_m: GaVec3::new(0.01, 0.03, 0.02),
        body_from_ellipsoid: Se3::from_parts(rotation(0.4, -0.7, 0.2), GaVec3::new(0.2, -0.1, 0.3))
            .unwrap(),
    };
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = Vec3::new(0.01, 0.02, -0.01);
    let gravity = Gravity::new(Vec3::new(0.0, -0.001, 0.0)).unwrap();
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        gravity,
        ambient,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let world_from_body =
        Se3::from_parts(rotation(-0.3, 0.1, 0.6), GaVec3::new(1.0, 2.0, -1.0)).unwrap();
    let slip = GaVec3::new(0.03, -0.02, 0.01);
    let spin = GaVec3::new(0.06, 0.04, -0.03);
    let body_rotation = boundary.body_from_ellipsoid.rotation();
    let omega_body = body_rotation.rotate(spin).unwrap();
    let center_velocity_body = world_from_body
        .rotation()
        .inverse()
        .rotate(gv(ambient))
        .unwrap()
        + body_rotation.rotate(slip).unwrap();
    let origin_velocity_body =
        center_velocity_body - omega_body.cross(boundary.body_from_ellipsoid.translation());
    let state = native
        .state_from_kinematics(
            world_from_body,
            Twist::new(omega_body, origin_velocity_body),
        )
        .unwrap();
    let load = medium.load(&native, state, || false).unwrap();
    let c = medium.resistance().translation_kg_s;
    let d = medium.resistance().rotation_kg_m2_s;
    let drag_ellipsoid = GaVec3::new(-c[0] * slip.x, -c[1] * slip.y, -c[2] * slip.z);
    let torque_ellipsoid = GaVec3::new(-d[0] * spin.x, -d[1] * spin.y, -d[2] * spin.z);
    let force_world = world_from_body
        .rotation()
        .rotate(body_rotation.rotate(drag_ellipsoid).unwrap())
        .unwrap();
    let torque_principal = native
        .body_from_principal()
        .rotation()
        .inverse()
        .rotate(body_rotation.rotate(torque_ellipsoid).unwrap())
        .unwrap();
    vector_close(load.drag_force_world_n, mv(force_world), 2e-15);
    vector_close(load.rotary_torque_body_nm, mv(torque_principal), 2e-17);
    let power = -drag_ellipsoid.dot(slip) - torque_ellipsoid.dot(spin);
    relative(load.dissipated_power_w, power, 2e-13);
    vector_close(
        load.buoyancy_force_world_n,
        gravity
            .acceleration_world()
            .scale(-material.density_kg_m3() * boundary.volume_moments().unwrap().volume_m3()),
        1e-18,
    );
    let inertia = native.properties().principal_inertia_body();
    // Supplied semiaxes are deliberately unsorted, whereas the native inertia is sorted.
    assert!(inertia.x < inertia.y && inertia.y < inertia.z);
}

#[test]
fn driven_principal_axis_trajectory_refines_to_independent_translation_and_spin_solutions() {
    let (a, b, density) = (0.03, 0.01, 20000.0);
    let boundary = shape(a, b, b);
    let source = body(boundary, density);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let (c, d) = spheroid(a, b, 1.0);
    let mass = density * (4.0 / 3.0) * core::f64::consts::PI * a * b * b;
    let inertia_x = (2.0 / 5.0) * mass * b * b;
    let (terminal_v, terminal_w, time): (f64, f64, f64) = (0.08, 0.3, 0.8);
    let applied = Wrench {
        force_world: Vec3::new(c[0] * terminal_v, 0.0, 0.0),
        torque_body: mv(native
            .body_from_principal()
            .rotation()
            .inverse()
            .rotate(GaVec3::new(d[0] * terminal_w, 0.0, 0.0))
            .unwrap()),
    };
    let exact_v = -terminal_v * (-time * c[0] / mass).exp_m1();
    let exact_w = -terminal_w * (-time * d[0] / inertia_x).exp_m1();
    let exact_x = terminal_v * (time + mass / c[0] * (-time * c[0] / mass).exp_m1());
    let mut errors = Vec::new();
    for steps in [40_usize, 80, 160] {
        let mut state = native
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        let mut heat = 0.0;
        let mut work = 0.0;
        for _ in 0..steps {
            let result = medium
                .step_with_wrench(
                    &native,
                    state,
                    applied,
                    controls(time / steps as f64),
                    || false,
                )
                .unwrap()
                .flow;
            heat += result.viscous_heat_j;
            work += result.applied_work_j;
            state = result.dynamics.state_after;
        }
        let v = state
            .center_of_mass_velocity_world(native.properties())
            .unwrap();
        let omega = native
            .properties()
            .angular_velocity_body_checked(state.angular_momentum_body())
            .unwrap();
        let omega_world = state.pose().orientation().rotate_body_to_world(omega);
        errors.push(
            (v.x - exact_v).abs()
                + (omega_world.x - exact_w).abs()
                + (state.pose().position_world().x - exact_x).abs(),
        );
        assert!(v.y.abs() + v.z.abs() + omega_world.y.abs() + omega_world.z.abs() < 1e-13);
        assert!(heat > 0.0);
        assert!((state.kinetic_energy(native.properties()).unwrap() + heat - work).abs() < 2e-12);
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8, "errors: {errors:?}");
    }
}

#[test]
fn rotating_triaxial_body_is_world_covariant_and_closes_drive_flow_and_heat_work() {
    let boundary = EllipsoidBoundary {
        semiaxes_m: GaVec3::new(0.01, 0.025, 0.016),
        body_from_ellipsoid: Se3::from_parts(rotation(0.3, -0.2, 0.1), GaVec3::new(0.0, 0.0, 0.0))
            .unwrap(),
    };
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let rotate_world = rotation(-0.5, 0.4, 0.3);
    let map = |v: Vec3| mv(rotate_world.rotate(gv(v)).unwrap());
    let gravity = Gravity::new(Vec3::new(0.0, -0.001, 0.0)).unwrap();
    let rotated_gravity = Gravity::new(map(gravity.acceleration_world())).unwrap();
    let ambient = Vec3::new(0.01, -0.02, 0.015);
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        gravity,
        ambient,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let rotated_medium = CreepingEllipsoid::new(
        &material,
        boundary,
        rotated_gravity,
        map(ambient),
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let pose = Se3::from_parts(rotation(0.1, 0.2, -0.3), GaVec3::new(0.0, 0.0, 0.0)).unwrap();
    let rotated_pose = Se3::from_parts(rotate_world, GaVec3::new(0.0, 0.0, 0.0))
        .unwrap()
        .compose(pose)
        .unwrap();
    let twist = Twist::new(
        GaVec3::new(0.12, -0.07, 0.09),
        GaVec3::new(0.03, 0.02, -0.01),
    );
    let mut state = native.state_from_kinematics(pose, twist).unwrap();
    let initial_orientation = state.pose().orientation();
    let mut rotated_state = native.state_from_kinematics(rotated_pose, twist).unwrap();
    let applied = Wrench {
        force_world: Vec3::new(0.004, -0.003, 0.002),
        torque_body: Vec3::new(2e-5, -3e-5, 1e-5),
    };
    let rotated_applied = Wrench {
        force_world: map(applied.force_world),
        ..applied
    };
    let mut totals = [0.0; 4];
    for _ in 0..64 {
        let candidate = medium
            .step_with_wrench(&native, state, applied, controls(0.005), || false)
            .unwrap();
        let rotated_candidate = rotated_medium
            .step_with_wrench(
                &native,
                rotated_state,
                rotated_applied,
                controls(0.005),
                || false,
            )
            .unwrap();
        let a = candidate.flow;
        let b = rotated_candidate.flow;
        // This is the original dynamics owner, with the physically accepted midpoint load.
        let reference = RigidBodyIntegrator::new(gravity)
            .step(
                state,
                native.properties(),
                Wrench {
                    force_world: a.load_midpoint.wrench.force_world.add(applied.force_world),
                    torque_body: a.load_midpoint.wrench.torque_body.add(applied.torque_body),
                },
                0.005,
            )
            .unwrap();
        // Picard closure tolerance permits a tiny last-iterate torque difference.
        vector_close(
            a.dynamics.state_after.angular_momentum_body(),
            reference.state_after.angular_momentum_body(),
            2e-15,
        );
        state = a.dynamics.state_after;
        rotated_state = b.dynamics.state_after;
        vector_close(
            rotated_state.linear_momentum_world(),
            map(state.linear_momentum_world()),
            2e-13,
        );
        vector_close(
            rotated_state.pose().position_world(),
            map(state.pose().position_world()),
            2e-13,
        );
        vector_close(
            rotated_state.angular_momentum_body(),
            state.angular_momentum_body(),
            2e-15,
        );
        relative(a.viscous_heat_j, b.viscous_heat_j, 2e-11);
        totals[0] += a.stored_energy_change_j;
        totals[1] += a.viscous_heat_j;
        totals[2] += a.background_work_j;
        totals[3] += a.applied_work_j;
        assert!(candidate.rotational_iterations > 0);
    }
    assert!(totals[1] > 0.0);
    assert!((totals[0] + totals[1] - totals[2] - totals[3]).abs() < 2e-12);
    assert_ne!(state.pose().orientation(), initial_orientation);
}

#[test]
fn source_viscosity_changes_actual_anisotropic_response_at_the_same_solid_inertia() {
    let boundary = shape(0.01, 0.025, 0.016);
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let a = fluid(1.0, 1.0);
    let b = fluid(1.0, 3.0);
    let first_medium = CreepingEllipsoid::new(
        &a,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let second_medium = CreepingEllipsoid::new(
        &b,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let state = native
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(GaVec3::new(0.0, 0.0, 0.0), GaVec3::new(0.03, 0.02, -0.01)),
        )
        .unwrap();
    let first = first_medium
        .step(&native, state, controls(0.005), || false)
        .unwrap()
        .flow;
    let second = second_medium
        .step(&native, state, controls(0.005), || false)
        .unwrap()
        .flow;
    assert_eq!(
        first.dynamics.diagnostics_before,
        second.dynamics.diagnostics_before
    );
    assert!(
        magnitude(second.dynamics.state_after.linear_momentum_world())
            < magnitude(first.dynamics.state_after.linear_momentum_world())
    );
    assert!(second.viscous_heat_j > first.viscous_heat_j);
}

#[test]
fn geometry_regime_iteration_and_cancellation_refuse_without_advancing_state() {
    let boundary = shape(0.03, 0.01, 0.02);
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let state = native
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(
                GaVec3::new(0.1, 0.07, -0.05),
                GaVec3::new(0.03, 0.02, -0.01),
            ),
        )
        .unwrap();
    let saved = state;
    let mut one_iteration = controls(0.005);
    one_iteration.maximum_coupling_iterations = 1;
    assert!(matches!(
        medium.step(&native, state, one_iteration, || false),
        Err(EllipsoidFlowError::Coupling { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        medium.step(&native, state, controls(0.005), || {
            polls += 1;
            polls == 3
        }),
        Err(EllipsoidFlowError::Cancelled)
    ));
    assert_eq!(state, saved);
    let wrong_source = body(shape(0.01, 0.01, 0.01), 20000.0);
    let wrong = wrong_source.principal_body(1e-12, || false).unwrap();
    assert!(medium.load(&wrong, state, || false).is_err());
    let fast = native
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(GaVec3::new(0.0, 0.0, 0.0), GaVec3::new(10.0, 0.0, 0.0)),
        )
        .unwrap();
    assert!(medium.load(&native, fast, || false).is_err());
    let dense = fluid(1e6, 1.0);
    let invalid_unsteady = CreepingEllipsoid::new(
        &dense,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap();
    let rest = native
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    assert!(
        invalid_unsteady
            .step(&native, rest, controls(0.005), || false)
            .is_err()
    );
    assert!(
        medium
            .step_with_wrench(
                &native,
                state,
                Wrench {
                    force_world: Vec3::new(f64::NAN, 0.0, 0.0),
                    ..Wrench::ZERO
                },
                controls(0.005),
                || false,
            )
            .is_err()
    );
    assert!(
        CreepingEllipsoid::new(
            &material,
            shape(0.1, 0.0001, 0.0001),
            Gravity::ZERO,
            Vec3::ZERO,
            limits(),
            EllipsoidQuadrature {
                initial_panels: 1,
                maximum_panels: 2,
                relative_tolerance: 1e-13,
            },
            || false,
        )
        .is_err()
    );
    assert!(
        medium
            .step(&native, state, controls(0.005), || false)
            .is_ok()
    );
    assert!(VolumeMoments::ellipsoid(GaVec3::new(0.0, 0.01, 0.02)).is_err());
    assert_eq!(state, saved);
}
