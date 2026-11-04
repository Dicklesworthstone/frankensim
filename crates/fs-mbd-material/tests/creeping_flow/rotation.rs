//! G1/G3/G4: independent rotary impedance, continuum spin-up and native exchange.
use super::*;
use fs_ga::{So3, So3Tangent};
use fs_mbd::{Vec3 as NativeVec3, Wrench};
use fs_mbd_material::creeping_flow::unsteady::rotation::RotatingStepControls;
use fs_mbd_material::creeping_flow::unsteady::{BassetResolution, UnsteadyStepControls};

fn resolution(lo: f64, hi: f64, intervals: usize) -> BassetResolution {
    BassetResolution {
        minimum_wavenumber_radius: lo,
        maximum_wavenumber_radius: hi,
        intervals,
    }
}
fn controls(dt: f64) -> RotatingStepControls {
    RotatingStepControls {
        mechanical: UnsteadyStepControls {
            duration_s: dt,
            maximum_energy_defect_j: 1e-18,
            maximum_momentum_defect_ns: 1e-16,
        },
        maximum_angular_momentum_defect_nms: 1e-17,
        maximum_rotation_increment_rad: 0.2,
        relative_torque_tolerance: 1e-12,
        maximum_coupling_iterations: 32,
    }
}
fn magnitude(v: NativeVec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn vector_near(a: NativeVec3, b: NativeVec3, tolerance: f64) {
    assert!(magnitude(a.sub(b)) <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn rotary_memory_matches_the_exact_impedance_and_has_no_instantaneous_added_inertia() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let d = 8.0 * core::f64::consts::PI * RADIUS.powi(3);
    let tau = RADIUS.powi(2) / 0.001;
    for omega in [1.0, 10.0, 100.0] {
        let p = omega * tau;
        let root = (p / 2.0).sqrt();
        let denominator = (1.0 + root).powi(2) + root * root;
        let exact = [
            d * (1.0 + p * root / (3.0 * denominator)),
            d * p * (1.0 + root) / (3.0 * denominator),
        ];
        let mut previous = f64::INFINITY;
        for (lo, hi, n) in [
            (0.02, 40.0, 32),
            (0.001, 1000.0, 128),
            (0.0001, 10000.0, 256),
        ] {
            let motion = medium
                .unsteady_rotating_motion(&body, rest, resolution(lo, hi, n), || false)
                .unwrap();
            let actual = motion.rotational_impedance_kg_m2_s(omega).unwrap();
            let error = (actual[0] - exact[0]).hypot(actual[1] - exact[1]);
            assert!(
                error < 0.15 * previous,
                "rotary spectral ladder: {error:e} >= {previous:e}"
            );
            assert!(error <= 1.01 * motion.rotational_omitted_tail_bound_kg_m2_s(omega).unwrap());
            assert!(
                motion
                    .rotational_branches()
                    .iter()
                    .all(|b| b.resistance_kg_m2_s > 0.0 && b.inertia_kg_m2 > 0.0)
            );
            previous = error;
        }
        assert!(previous / exact[0].hypot(exact[1]) < 3e-4);
    }
    let motion = medium
        .unsteady_rotating_motion(&body, rest, resolution(1e-4, 1e4, 256), || false)
        .unwrap();
    let fluid_inertia: f64 = motion
        .rotational_branches()
        .iter()
        .map(|b| b.inertia_kg_m2)
        .sum();
    let continuum_band =
        2.0 * d * tau / (3.0 * core::f64::consts::PI) * (10000.0_f64.atan() - 0.0001_f64.atan());
    assert!((fluid_inertia / continuum_band - 1.0).abs() < 1e-6);
    // The finite low-frequency entrained inertia is carried by memory. It must
    // not reduce the instantaneous torque / original dry inertia acceleration.
    let applied = Wrench {
        torque_body: NativeVec3::new(d * 0.01, 0.0, 0.0),
        ..Wrench::ZERO
    };
    let dt = 1e-12;
    let (next, _) = motion
        .step_with_wrench(applied, controls(dt), || false)
        .unwrap();
    let omega = body
        .properties()
        .angular_velocity_body_checked(next.state().angular_momentum_body())
        .unwrap()
        .x;
    let expected = applied.torque_body.x / body.properties().principal_inertia_body().x;
    assert!((omega / dt / expected - 1.0).abs() < 1e-6);
    assert!(fluid_inertia > body.properties().principal_inertia_body().x);
}

// Independent branch-cut inversion of 1/(I*s+Z_rot(s)), using the analytic
// continuum impedance rather than the product's finite memory branches.
fn continuum_spin_fraction(time_s: f64, inertia: f64, d: f64, tau: f64) -> f64 {
    let panels = 32768;
    let h = 24.0 / panels as f64;
    let kappa = inertia / (d * tau);
    let mut total = 0.0;
    for i in 0..=panels {
        let x = (-12.0 + i as f64 * h).exp();
        let x2 = x * x;
        let real = 1.0 - x2 / (3.0 * (1.0 + x2)) - kappa * x2;
        let imaginary = x * x2 / (3.0 * (1.0 + x2));
        let integrand = (2.0 / core::f64::consts::PI) * imaginary
            / (real * real + imaginary * imaginary)
            * -(-x2 * time_s / tau).exp_m1();
        total += integrand
            * if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
    }
    total * h / 3.0
}

#[test]
fn native_spin_up_refines_to_the_independent_continuum_response() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let d = 8.0 * core::f64::consts::PI * RADIUS.powi(3);
    let expected =
        0.01 * continuum_spin_fraction(0.1, body.properties().principal_inertia_body().x, d, 0.1);
    let mut previous = f64::INFINITY;
    for (lo, hi, n, steps) in [
        (0.01, 100.0, 32, 100),
        (0.001, 1000.0, 128, 400),
        (0.0001, 10000.0, 256, 1600),
    ] {
        let mut motion = medium
            .unsteady_rotating_motion(&body, rest, resolution(lo, hi, n), || false)
            .unwrap();
        let mut energy = [0.0; 3];
        for _ in 0..steps {
            let (next, receipt) = motion
                .step_with_wrench(
                    Wrench {
                        torque_body: NativeVec3::new(d * 0.01, 0.0, 0.0),
                        ..Wrench::ZERO
                    },
                    controls(0.1 / steps as f64),
                    || false,
                )
                .unwrap();
            assert!(receipt.rotation.unwrap().viscous_heat_j >= 0.0);
            energy[0] += receipt.stored_energy_change_j;
            energy[1] += receipt.viscous_heat_j;
            energy[2] += receipt.applied_work_j;
            motion = next;
        }
        let omega = body
            .properties()
            .angular_velocity_body_checked(motion.state().angular_momentum_body())
            .unwrap()
            .x;
        let error = (omega / expected - 1.0).abs();
        assert!(
            error < previous * 0.15,
            "continuum spin-up ladder {error:e} >= {previous:e}"
        );
        assert!((energy[0] + energy[1] - energy[2]).abs() < 1e-16 * energy[2].abs().max(1.0));
        previous = error;
    }
    assert!(previous < 1e-5);
}

#[test]
fn changing_three_axis_drives_close_native_linear_and_angular_exchange_with_memory() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.4))).unwrap(),
        Vec3::new(1.0, -2.0, 3.0),
    )
    .unwrap();
    let state = body.state_from_kinematics(pose, Twist::zero()).unwrap();
    let mut motion = medium
        .unsteady_rotating_motion(&body, state, resolution(1e-4, 1e4, 192), || false)
        .unwrap();
    let mut totals = [0.0; 3];
    for i in 0..300 {
        let applied = if i < 100 {
            Wrench {
                force_world: NativeVec3::new(1e-4, -2e-4, 3e-4),
                torque_body: NativeVec3::new(1e-7, 3e-7, -2e-7),
            }
        } else if i < 200 {
            Wrench {
                force_world: NativeVec3::new(-2e-4, 1e-4, 0.0),
                torque_body: NativeVec3::new(-3e-7, 2e-7, 4e-7),
            }
        } else {
            Wrench::ZERO
        };
        let before = motion.state();
        let (next, audit) = motion
            .step_with_wrench(applied, controls(0.001), || false)
            .unwrap();
        let exact_native = RigidBodyIntegrator::new(Gravity::ZERO)
            .step_cayley(
                before,
                body.properties(),
                Wrench {
                    force_world: audit.fluid_wrench.force_world.add(applied.force_world),
                    torque_body: audit.fluid_wrench.torque_body.add(applied.torque_body),
                },
                0.001,
            )
            .unwrap();
        assert_eq!(audit.dynamics, exact_native);
        let spin = audit.rotation.unwrap();
        assert!(spin.viscous_heat_j >= 0.0);
        assert!(magnitude(spin.angular_momentum_residual_nms) < 1e-17);
        assert!(magnitude(audit.momentum_residual_ns) < 1e-16);
        assert!(audit.energy_defect_j.abs() < 1e-18);
        let omega_body = body
            .properties()
            .angular_velocity_body_checked(
                before
                    .angular_momentum_body()
                    .scale(0.5)
                    .add(next.state().angular_momentum_body().scale(0.5)),
            )
            .unwrap();
        let omega_world = before.pose().orientation().rotate_body_to_world(omega_body);
        // Independently reconstruct fluid storage and physical history torque.
        let mut memory_energy = 0.0;
        let mut history = NativeVec3::ZERO;
        for ((branch, u0), u1) in motion
            .rotational_branches()
            .iter()
            .zip(motion.fluid_angular_velocities_world_rad_s())
            .zip(next.fluid_angular_velocities_world_rad_s())
        {
            let midpoint = u0.scale(0.5).add(u1.scale(0.5));
            memory_energy += 0.5 * branch.inertia_kg_m2 * (u1.norm_squared() - u0.norm_squared());
            history = history.sub(omega_world.sub(midpoint).scale(branch.resistance_kg_m2_s));
        }
        near(memory_energy, spin.memory_energy_change_j, 1e-23);
        vector_near(history, spin.history_torque_world_nm, 1e-18);
        totals[0] += audit.stored_energy_change_j;
        totals[1] += audit.viscous_heat_j;
        totals[2] += audit.applied_work_j;
        motion = next;
    }
    assert!(
        motion
            .fluid_angular_velocities_world_rad_s()
            .iter()
            .any(|v| magnitude(*v) > 0.0)
    );
    assert!((totals[0] + totals[1] - totals[2]).abs() < 1e-18);
}

#[test]
fn world_rotation_transforms_both_fluid_histories_without_a_spurious_body_frame_memory() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(-0.5, 0.2, 0.7))).unwrap();
    let map = |v: NativeVec3| {
        let v = rotation.rotate(Vec3::new(v.x, v.y, v.z)).unwrap();
        NativeVec3::new(v.x, v.y, v.z)
    };
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let rotated = body
        .state_from_kinematics(
            Se3::from_parts(rotation, Vec3::new(0.0, 0.0, 0.0)).unwrap(),
            Twist::zero(),
        )
        .unwrap();
    let mut first = medium
        .unsteady_rotating_motion(&body, state, resolution(1e-4, 1e4, 128), || false)
        .unwrap();
    let mut second = medium
        .unsteady_rotating_motion(&body, rotated, resolution(1e-4, 1e4, 128), || false)
        .unwrap();
    for i in 0..100 {
        let applied = Wrench {
            force_world: NativeVec3::new(1e-4, -2e-4, 3e-4),
            torque_body: if i < 50 {
                NativeVec3::new(3e-7, -1e-7, 2e-7)
            } else {
                NativeVec3::new(-2e-7, 4e-7, 1e-7)
            },
        };
        let (a, ar) = first
            .step_with_wrench(applied, controls(0.001), || false)
            .unwrap();
        let (b, br) = second
            .step_with_wrench(
                Wrench {
                    force_world: map(applied.force_world),
                    ..applied
                },
                controls(0.001),
                || false,
            )
            .unwrap();
        vector_near(
            b.state().linear_momentum_world(),
            map(a.state().linear_momentum_world()),
            1e-18,
        );
        vector_near(
            b.state().angular_momentum_body(),
            a.state().angular_momentum_body(),
            1e-20,
        );
        for (x, y) in a
            .fluid_angular_velocities_world_rad_s()
            .iter()
            .zip(b.fluid_angular_velocities_world_rad_s())
        {
            vector_near(*y, map(*x), 1e-13);
        }
        near(ar.viscous_heat_j, br.viscous_heat_j, 1e-22);
        first = a;
        second = b;
    }
}

#[test]
fn cancelled_or_invalid_rotation_never_advances_either_fluid_history() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let motion = medium
        .unsteady_rotating_motion(&body, rest, resolution(1e-4, 1e4, 64), || false)
        .unwrap();
    let drive = Wrench {
        force_world: NativeVec3::new(1e-4, 0.0, 0.0),
        torque_body: NativeVec3::new(1e-7, -2e-7, 3e-7),
    };
    let (motion, _) = motion
        .step_with_wrench(drive, controls(0.001), || false)
        .unwrap();
    let saved_state = motion.state();
    let saved_time = motion.time_s();
    let saved_translation = motion.fluid_velocities_world_m_s().to_vec();
    let saved_rotation = motion.fluid_angular_velocities_world_rad_s().to_vec();
    let mut calls = 0;
    let (expected, expected_receipt) = motion
        .step_with_wrench(drive, controls(0.001), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls - 1, calls] {
        let mut seen = 0;
        assert_eq!(
            motion
                .step_with_wrench(drive, controls(0.001), || {
                    seen += 1;
                    seen == stop
                })
                .err(),
            Some(CreepingFlowError::Cancelled)
        );
    }
    let mut invalid = controls(0.001);
    invalid.maximum_coupling_iterations = 0;
    assert!(motion.step_with_wrench(drive, invalid, || false).is_err());
    let mut invalid = controls(0.001);
    invalid.maximum_rotation_increment_rad = 1e-15;
    assert!(motion.step_with_wrench(drive, invalid, || false).is_err());
    assert!(
        motion
            .step_with_wrench(
                Wrench {
                    torque_body: NativeVec3::new(f64::NAN, 0.0, 0.0),
                    ..drive
                },
                controls(0.001),
                || false
            )
            .is_err()
    );
    assert!(
        motion
            .step_with_wrench(
                Wrench {
                    torque_body: NativeVec3::new(1e6, 0.0, 0.0),
                    ..drive
                },
                controls(0.001),
                || false
            )
            .is_err()
    );
    assert_eq!(motion.state(), saved_state);
    assert_eq!(motion.time_s(), saved_time);
    assert_eq!(motion.fluid_velocities_world_m_s(), saved_translation);
    assert_eq!(
        motion.fluid_angular_velocities_world_rad_s(),
        saved_rotation
    );
    let (retry, receipt) = motion
        .step_with_wrench(drive, controls(0.001), || false)
        .unwrap();
    assert_eq!(receipt, expected_receipt);
    assert_eq!(retry.state(), expected.state());
    assert_eq!(
        retry.fluid_angular_velocities_world_rad_s(),
        expected.fluid_angular_velocities_world_rad_s()
    );
    let original = medium
        .unsteady_motion(&body, rest, resolution(1e-4, 1e4, 64), || false)
        .unwrap();
    assert!(
        original
            .step_with_wrench(drive, controls(0.001).mechanical, || false)
            .is_err()
    );
}
