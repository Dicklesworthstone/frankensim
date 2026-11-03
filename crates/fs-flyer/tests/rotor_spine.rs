//! G0/G1/G3 full carrier/rotor momentum, work and phase integration.
use fs_flyer::spine::{
    Loads, RigidBody, SixDofState,
    rotors::{
        AxialRotor, CarrierImpulse, RotorAssembly, RotorImpulses, RotorState, RotorTorques,
        apply_impulse, geared, step_held,
    },
};
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};

const ZERO: [f64; 3] = [0.0; 3];
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn rotation(s: &SixDofState) -> So3 {
    So3::try_from_quat(Quat {
        w: s.quat[0],
        x: s.quat[1],
        y: s.quat[2],
        z: s.quat[3],
    })
    .unwrap()
}
fn rigid() -> SixDofState {
    SixDofState {
        pos_m: [2.0, -3.0, 1.0],
        vel_mps: ZERO,
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: ZERO,
    }
}
fn body() -> RigidBody {
    RigidBody::new(3.0, [1.0, 1.5, 2.0]).unwrap()
}
fn rotor(axis_body: [f64; 3]) -> AxialRotor {
    AxialRotor {
        axis_body,
        axial_inertia_kgm2: 0.2,
    }
}
fn no_loads() -> Loads {
    Loads {
        force_n: ZERO,
        moment_nm: ZERO,
    }
}
fn total<const N: usize>(
    a: &RotorAssembly<N>,
    s: &RotorState<N>,
    origin: [f64; 3],
) -> (Vec3, Vec3) {
    let p = v(s.rigid.vel_mps).scale(a.reduced_body.mass_kg);
    let mut l = a
        .reduced_body
        .inertia_tensor()
        .matrix()
        .apply(v(s.rigid.omega_body));
    for i in 0..N {
        l = l + v(a.rotors[i].axis_body).scale(s.axial_momentum[i]);
    }
    (
        p,
        rotation(&s.rigid).rotate(l).unwrap() + (v(s.rigid.pos_m) - v(origin)).cross(p),
    )
}

#[test]
fn g1_motor_recoil_phase_and_energy_are_analytic() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let s = a.state_from_relative_spin(rigid(), [0.0], [0.0]).unwrap();
    let dt = 0.05;
    let torque = 0.7;
    let r = step_held(
        &a,
        &s,
        no_loads(),
        RotorTorques {
            motor_nm: [torque],
            aerodynamic_nm: [0.0],
        },
        ZERO,
        dt,
    )
    .unwrap();
    let relative = torque * dt * (1.0 / 0.2 + 1.0);
    assert!((r.state.axial_momentum[0] - torque * dt).abs() < 1e-15);
    assert!((r.state.rigid.omega_body[0] + torque * dt).abs() < 1e-15);
    assert!((r.relative_speed_rad_s[0] - relative).abs() < 1e-15);
    assert!((r.state.relative_phase_rad[0] - 0.5 * dt * relative).abs() < 1e-15);
    assert!((rotation(&r.state.rigid).log().angular.x + 0.5 * torque * dt * dt).abs() < 1e-15);
    assert!(v(r.angular_momentum_residual).norm() < 1e-15);
    assert!(r.kinetic_work_residual_j.abs() < 1e-15);
    assert!((r.motor_work_j - torque * r.state.relative_phase_rad[0]).abs() < 1e-15);
    assert_eq!(r.external_work_j, 0.0);
}

#[test]
fn g1_full_tensor_carrier_has_analytic_cross_axis_motor_recoil() {
    let a = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            fs_ga::Mat3 {
                m: [2.0, 0.6, 0.0, 0.6, 1.5, 0.0, 0.0, 0.0, 2.5],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let s = a.state_from_relative_spin(rigid(), [0.0], [0.0]).unwrap();
    let dt = 0.05;
    let torque = 0.7;
    let r = step_held(
        &a,
        &s,
        no_loads(),
        RotorTorques {
            motor_nm: [torque],
            aerodynamic_nm: [0.0],
        },
        ZERO,
        dt,
    )
    .unwrap();
    let p = torque * dt;
    let determinant = 2.0 * 1.5 - 0.6 * 0.6;
    // Initially zero total angular momentum stays zero: I omega = -p axis.
    // Independent inverse of the physical xy inertia block supplies recoil.
    let expected = Vec3::new(-1.5 * p / determinant, 0.6 * p / determinant, 0.0);
    assert!((v(r.state.rigid.omega_body) - expected).norm() < 2e-14);
    assert!((r.state.axial_momentum[0] - p).abs() < 1e-15);
    assert!(v(r.angular_momentum_residual).norm() < 2e-14);
    assert!(
        (r.kinetic_energy_after_j - 0.5 * p * p * (1.0 / 0.2 + 1.5 / determinant)).abs() < 2e-14
    );
    assert!(r.kinetic_work_residual_j.abs() < 2e-14);
}

#[test]
fn g1_aerodynamic_braking_is_external_without_spurious_motor_reaction() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([0.0, 0.0, 1.0])],
    };
    let s = a.state_from_relative_spin(rigid(), [0.3], [10.0]).unwrap();
    let dt = 0.05;
    let drag = -0.4;
    let r = step_held(
        &a,
        &s,
        no_loads(),
        RotorTorques {
            motor_nm: [0.0],
            aerodynamic_nm: [drag],
        },
        ZERO,
        dt,
    )
    .unwrap();
    assert!(v(r.state.rigid.omega_body).norm() < 1e-14);
    assert!(rotation(&r.state.rigid).log().angular.norm() < 1e-14);
    assert!((r.state.axial_momentum[0] - (2.0 + drag * dt)).abs() < 1e-15);
    assert!((r.external_angular_impulse_world[2] - drag * dt).abs() < 1e-15);
    assert_eq!(r.motor_work_j, 0.0);
    assert!(r.external_work_j < 0.0);
    assert!(r.kinetic_work_residual_j.abs() < 2e-15);
}

#[test]
fn g0_multiple_driven_rotors_close_spin_orbital_and_external_impulses() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [
            rotor([1.0, 0.0, 0.0]),
            rotor([0.0, 1.0, 0.0]),
            rotor([0.0, 0.0, 1.0]),
        ],
    };
    let mut initial = rigid();
    initial.omega_body = [0.4, -0.2, 0.3];
    initial.vel_mps = [0.7, 0.3, -0.2];
    let mut state = a
        .state_from_relative_spin(initial, [0.1, -0.3, 0.2], [12.0, -7.0, 5.0])
        .unwrap();
    let origin = [-3.0, 1.0, -2.0];
    let (p0, l0) = total(&a, &state, origin);
    let mut jp = Vec3::default();
    let mut jl = Vec3::default();
    for _ in 0..500 {
        let r = step_held(
            &a,
            &state,
            Loads {
                force_n: [0.2, -0.1, 0.3],
                moment_nm: [0.01, -0.02, 0.03],
            },
            RotorTorques {
                motor_nm: [0.04, -0.03, 0.02],
                aerodynamic_nm: [-0.05, 0.04, -0.02],
            },
            origin,
            0.002,
        )
        .unwrap();
        assert!(v(r.linear_momentum_residual).norm() < 3e-15);
        assert!(v(r.angular_momentum_residual).norm() < 3e-14);
        jp = jp + v(r.external_linear_impulse_world);
        jl = jl + v(r.external_angular_impulse_world);
        state = r.state;
    }
    let (p1, l1) = total(&a, &state, origin);
    assert!((p1 - p0 - jp).norm() < 2e-12);
    assert!((l1 - l0 - jl).norm() < 2e-11);
}

#[test]
fn g0_balanced_counter_rotors_reduce_to_native_free_carrier() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 2],
    };
    let mut initial = rigid();
    initial.omega_body = [0.3, -0.4, 0.2];
    let state = RotorState {
        rigid: initial,
        axial_momentum: [2.0, -2.0],
        relative_phase_rad: [0.0; 2],
    };
    let r = step_held(
        &a,
        &state,
        no_loads(),
        RotorTorques {
            motor_nm: [0.3, -0.3],
            aerodynamic_nm: [-0.2, 0.2],
        },
        ZERO,
        0.03,
    )
    .unwrap();
    let (exact_r, exact_w) = fs_time::lie::rigid_body_step(
        rotation(&initial),
        v(initial.omega_body),
        body().inertia_tensor().principal_moments(),
        0.03,
    )
    .unwrap();
    assert_eq!(rotation(&r.state.rigid), exact_r);
    assert_eq!(v(r.state.rigid.omega_body), exact_w);
    assert_eq!(r.external_angular_impulse_world, ZERO);
    let empty = RotorAssembly {
        reduced_body: body(),
        rotors: [],
    };
    let zero = empty.state_from_relative_spin(initial, [], []).unwrap();
    let r0 = step_held(
        &empty,
        &zero,
        no_loads(),
        RotorTorques {
            motor_nm: [],
            aerodynamic_nm: [],
        },
        ZERO,
        0.03,
    )
    .unwrap();
    assert_eq!(r.state.rigid, r0.state.rigid);
}

#[test]
fn g0_reduced_inertia_and_relative_spin_match_independent_locked_assembly() {
    let t = 0.5_f64.sqrt();
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([t, t, 0.0])],
    };
    let mut initial = rigid();
    initial.omega_body = [0.4, -0.2, 0.3];
    let relative = 3.0;
    let state = a
        .state_from_relative_spin(initial, [0.0], [relative])
        .unwrap();
    let projected = v(a.rotors[0].axis_body).dot(v(initial.omega_body));
    let j = a.rotors[0].axial_inertia_kgm2;
    assert_eq!(state.axial_momentum[0], j * (relative + projected));
    let carrier_e = 0.5
        * v(initial.omega_body).dot(
            body()
                .inertia_tensor()
                .matrix()
                .apply(v(initial.omega_body)),
        );
    let locked_e = carrier_e
        + 0.5 * j * projected * projected
        + j * relative * projected
        + 0.5 * j * relative * relative;
    let r = step_held(
        &a,
        &state,
        no_loads(),
        RotorTorques {
            motor_nm: [0.0],
            aerodynamic_nm: [0.0],
        },
        ZERO,
        0.01,
    )
    .unwrap();
    assert!((r.kinetic_energy_before_j - locked_e).abs() < 2e-15);
    assert!(
        (r.state.rigid.omega_body[2] - initial.omega_body[2]).abs() > 1e-4,
        "stored rotor momentum must actually drive gyroscopic precession"
    );
}

#[test]
fn g3_rotor_motion_phase_work_and_replay_converge_at_second_order() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]), rotor([0.0, 1.0, 0.0])],
    };
    let mut initial = rigid();
    initial.omega_body = [0.4, 0.7, 1.1];
    let state = a
        .state_from_relative_spin(initial, [0.2, -0.3], [4.0, -3.0])
        .unwrap();
    let run = |n: usize| {
        let mut s = state;
        let mut residual = 0.0;
        for _ in 0..n {
            let r = step_held(
                &a,
                &s,
                Loads {
                    force_n: [0.2, -0.1, 0.3],
                    moment_nm: [0.03, -0.02, 0.01],
                },
                RotorTorques {
                    motor_nm: [0.2, -0.1],
                    aerodynamic_nm: [-0.04, 0.03],
                },
                ZERO,
                1.0 / n as f64,
            )
            .unwrap();
            s = r.state;
            residual += r.kinetic_work_residual_j;
        }
        (s, residual)
    };
    let distance = |a: RotorState<2>, b: RotorState<2>| {
        rotation(&a.rigid)
            .body_minus(rotation(&b.rigid))
            .unwrap()
            .angular
            .norm()
            + (v(a.rigid.omega_body) - v(b.rigid.omega_body)).norm()
            + (a.relative_phase_rad[0] - b.relative_phase_rad[0]).abs()
            + (a.relative_phase_rad[1] - b.relative_phase_rad[1]).abs()
    };
    let coarse = run(40);
    let medium = run(80);
    let fine = run(160);
    assert_eq!(fine, run(160));
    let ratio = distance(coarse.0, medium.0) / distance(medium.0, fine.0);
    println!(
        "rotor spine order ratio {ratio}; work defects {}, {}, {}",
        coarse.1, medium.1, fine.1
    );
    assert!((3.7..4.3).contains(&ratio));
    assert!((3.7..4.3).contains(&(coarse.1 / medium.1)));
    assert!((3.7..4.3).contains(&(medium.1 / fine.1)));
}

#[test]
fn g0_admission_is_bounded_and_refusal_never_changes_state() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let s = a.state_from_relative_spin(rigid(), [0.0], [2.0]).unwrap();
    let original = s;
    let torques = RotorTorques {
        motor_nm: [0.0],
        aerodynamic_nm: [0.0],
    };
    for bad in [f64::NAN, f64::INFINITY] {
        let load = RotorTorques {
            motor_nm: [bad],
            ..torques
        };
        assert!(step_held(&a, &s, no_loads(), load, ZERO, 0.01).is_err());
    }
    let huge = RotorTorques {
        motor_nm: [f64::MAX],
        aerodynamic_nm: [f64::MAX],
    };
    assert!(step_held(&a, &s, no_loads(), huge, ZERO, 0.01).is_err());
    assert!(step_held(&a, &s, no_loads(), torques, ZERO, 0.0).is_err());
    let invalid = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([2.0, 0.0, 0.0])],
    };
    assert!(step_held(&invalid, &s, no_loads(), torques, ZERO, 0.01).is_err());
    let at_cap = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 64],
    };
    let at = at_cap
        .state_from_relative_spin(rigid(), [0.0; 64], [0.0; 64])
        .unwrap();
    assert!(
        step_held(
            &at_cap,
            &at,
            no_loads(),
            RotorTorques {
                motor_nm: [0.0; 64],
                aerodynamic_nm: [0.0; 64]
            },
            ZERO,
            0.01
        )
        .is_ok()
    );
    let excessive = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 65],
    };
    assert!(excessive.admit().is_err());
    assert_eq!(s, original);
    assert_eq!(
        step_held(&a, &s, no_loads(), torques, ZERO, 0.01).unwrap(),
        step_held(&a, &original, no_loads(), torques, ZERO, 0.01).unwrap()
    );
}

#[test]
fn g1_instant_motor_and_aero_impulses_have_distinct_analytic_recoil_and_work() {
    let a = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [2.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 4.0],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let s = a.state_from_relative_spin(rigid(), [0.37], [0.0]).unwrap();
    let jump = apply_impulse(
        &a,
        &s,
        RotorImpulses {
            carrier: CarrierImpulse::default(),
            motor_nms: [0.2],
            aerodynamic_nms: [-0.05],
        },
        ZERO,
    )
    .unwrap();
    // I^-1(-0.2,0,0) = (-0.12,0.04,0). Aerodynamic impulse
    // changes absolute rotor spin, with no duplicate carrier reaction.
    assert!((v(jump.state.rigid.omega_body) - v([-0.12, 0.04, 0.0])).norm() < 1e-14);
    assert!((jump.state.axial_momentum[0] - 0.15).abs() < 1e-15);
    assert!((jump.relative_speed_rad_s[0] - 0.87).abs() < 1e-14);
    assert!((jump.motor_work_j - 0.087).abs() < 1e-14);
    assert!((jump.external_work_j + 0.01875).abs() < 1e-14);
    assert!((jump.kinetic_energy_after_j - 0.06825).abs() < 1e-14);
    assert!(jump.kinetic_work_residual_j.abs() < 1e-14);
    assert!(
        (total(&a, &jump.state, ZERO).1 - total(&a, &s, ZERO).1 - v([-0.05, 0.0, 0.0])).norm()
            < 1e-14
    );
    assert_eq!(jump.state.rigid.pos_m, s.rigid.pos_m);
    assert_eq!(jump.state.rigid.quat, s.rigid.quat);
    assert_eq!(jump.state.relative_phase_rad, s.relative_phase_rad);
}

fn set_rotation(state: &mut SixDofState, r: So3) {
    let q = r.as_quat();
    state.quat = [q.w, q.x, q.y, q.z];
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}

#[test]
fn g3_off_center_impulse_closes_global_momentum_and_is_world_frame_covariant() {
    let a = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [3.0, 0.2, 0.1, 0.2, 4.0, 0.3, 0.1, 0.3, 5.0],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0]), rotor([0.0, 1.0, 0.0])],
    };
    let mut r = rigid();
    r.vel_mps = [2.0, -1.0, 3.0];
    r.omega_body = [0.2, -0.3, 0.4];
    set_rotation(
        &mut r,
        So3::exp(So3Tangent::new(v([0.3, -0.4, 0.2]))).unwrap(),
    );
    let s = a
        .state_from_relative_spin(r, [0.1, -0.2], [7.0, -4.0])
        .unwrap();
    let arm = [0.7, 0.2, -0.3];
    let impulse = [1.2, -0.8, 0.5];
    let origin = [4.0, 5.0, 6.0];
    let inputs = RotorImpulses {
        carrier: CarrierImpulse::at_body_point(&s.rigid, arm, impulse).unwrap(),
        motor_nms: [0.12, -0.07],
        aerodynamic_nms: [-0.03, 0.02],
    };
    let jump = apply_impulse(&a, &s, inputs, origin).unwrap();
    let world = rotation(&s.rigid);
    let point = v(s.rigid.pos_m) + world.rotate(v(arm)).unwrap();
    let expected_angular =
        (point - v(origin)).cross(v(impulse)) + world.rotate(v([-0.03, 0.02, 0.0])).unwrap();
    let (p0, h0) = total(&a, &s, origin);
    let (p1, h1) = total(&a, &jump.state, origin);
    assert!((p1 - p0 - v(impulse)).norm() < 1e-12);
    assert!((h1 - h0 - expected_angular).norm() < 1e-12);
    assert!((v(jump.external_angular_impulse_world) - expected_angular).norm() < 1e-12);
    assert!(v(jump.angular_momentum_residual).norm() < 1e-12);
    assert!(jump.kinetic_work_residual_j.abs() < 1e-12);
    assert!(
        (jump.kinetic_energy_after_j
            - jump.kinetic_energy_before_j
            - jump.motor_work_j
            - jump.external_work_j)
            .abs()
            < 1e-12
    );
    let q = So3::exp(So3Tangent::new(v([-0.2, 0.1, 0.3]))).unwrap();
    let mut transformed = s;
    transformed.rigid.pos_m = array(q.rotate(v(s.rigid.pos_m)).unwrap());
    transformed.rigid.vel_mps = array(q.rotate(v(s.rigid.vel_mps)).unwrap());
    set_rotation(&mut transformed.rigid, q.compose(world).unwrap());
    let rotated_inputs = RotorImpulses {
        carrier: CarrierImpulse::at_body_point(
            &transformed.rigid,
            arm,
            array(q.rotate(v(impulse)).unwrap()),
        )
        .unwrap(),
        ..inputs
    };
    let rotated = apply_impulse(
        &a,
        &transformed,
        rotated_inputs,
        array(q.rotate(v(origin)).unwrap()),
    )
    .unwrap();
    assert!((v(rotated.state.rigid.omega_body) - v(jump.state.rigid.omega_body)).norm() < 1e-12);
    assert!(
        (v(rotated.state.rigid.vel_mps) - q.rotate(v(jump.state.rigid.vel_mps)).unwrap()).norm()
            < 1e-12
    );
    assert!(
        (v(rotated.external_angular_impulse_world) - q.rotate(expected_angular).unwrap()).norm()
            < 1e-12
    );
    assert_eq!(rotated.state.relative_phase_rad, s.relative_phase_rad);
}

fn gear_momentum(a: &geared::GearedAssembly<2>, s: &geared::GearedState, origin: [f64; 3]) -> Vec3 {
    let mut spin = a.locked_inertia_kgm2().apply(v(s.rigid.omega_body));
    for rotor in a.rotors() {
        spin = spin
            + v(rotor.axis_body)
                .scale(rotor.axial_inertia_kgm2 * rotor.ratio * a.speed(s).unwrap());
    }
    rotation(&s.rigid).rotate(spin).unwrap()
        + (v(s.rigid.pos_m) - v(origin)).cross(v(s.rigid.vel_mps).scale(a.mass_kg()))
}

#[test]
fn g1_geared_impulse_retains_signed_constraints_and_closes_all_physical_ports() {
    let a = geared::GearedAssembly::new(
        3.0,
        Mat3 {
            m: [4.0, 0.3, 0.2, 0.3, 5.0, 0.1, 0.2, 0.1, 6.0],
        },
        [
            geared::GearedRotor {
                axis_body: [1.0, 0.0, 0.0],
                axial_inertia_kgm2: 0.3,
                ratio: 1.0,
            },
            geared::GearedRotor {
                axis_body: [0.0, 1.0, 0.0],
                axial_inertia_kgm2: 0.2,
                ratio: -2.0,
            },
        ],
    )
    .unwrap();
    let mut r = rigid();
    r.vel_mps = [2.0, -1.0, 3.0];
    r.omega_body = [0.2, -0.3, 0.4];
    let s = a.state_from_speed(r, 0.37, 7.0).unwrap();
    let internal = geared::apply_impulse(
        &a,
        &s,
        geared::GearedImpulses {
            carrier: CarrierImpulse::default(),
            motor_nms: 0.75,
            aerodynamic_nms: [0.0; 2],
        },
        ZERO,
    )
    .unwrap();
    assert!(
        (gear_momentum(&a, &internal.state, ZERO) - gear_momentum(&a, &s, ZERO)).norm() < 1e-12
    );
    assert!(internal.kinetic_work_residual_j.abs() < 1e-12);
    let input = geared::GearedImpulses {
        carrier: CarrierImpulse::at_body_point(&s.rigid, [0.2, -0.7, 0.4], [1.0, -0.5, 0.2])
            .unwrap(),
        motor_nms: 0.75,
        aerodynamic_nms: [-0.12, 0.08],
    };
    let jump = geared::apply_impulse(&a, &s, input, ZERO).unwrap();
    assert!((jump.state.momentum - s.momentum - 0.47).abs() < 1e-14);
    assert_eq!(
        jump.relative_speed_rad_s[1],
        -2.0 * jump.relative_speed_rad_s[0]
    );
    assert_eq!(jump.state.phase_rad, s.phase_rad);
    assert_eq!(jump.state.rigid.pos_m, s.rigid.pos_m);
    assert_eq!(jump.state.rigid.quat, s.rigid.quat);
    assert!(
        (gear_momentum(&a, &jump.state, ZERO)
            - gear_momentum(&a, &s, ZERO)
            - v(jump.external_angular_impulse_world))
        .norm()
            < 1e-12
    );
    let mean_w = (v(s.rigid.omega_body) + v(jump.state.rigid.omega_body)).scale(0.5);
    let mean_speed = 0.5 * (a.speed(&s).unwrap() + jump.speed_rad_s);
    let work = v(input.carrier.linear_world_ns)
        .dot((v(s.rigid.vel_mps) + v(jump.state.rigid.vel_mps)).scale(0.5))
        + v(input.carrier.angular_body_nms).dot(mean_w)
        + input.aerodynamic_nms[0] * (mean_speed + mean_w.x)
        + input.aerodynamic_nms[1] * (-2.0 * mean_speed + mean_w.y);
    assert!((jump.external_work_j - work).abs() < 1e-12);
    assert!((jump.motor_work_j - 0.75 * mean_speed).abs() < 1e-12);
    assert!(
        (jump.kinetic_energy_after_j
            - jump.kinetic_energy_before_j
            - jump.motor_work_j
            - jump.external_work_j)
            .abs()
            < 1e-12
    );
    assert!(jump.kinetic_work_residual_j.abs() < 1e-12);
    let mut unresolved = s;
    unresolved.momentum = 1e30;
    assert_eq!(
        geared::apply_impulse(&a, &unresolved, input, ZERO)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
}

#[test]
fn g0_impulses_refuse_dropped_coordinates_but_resolve_recoil_beside_large_other_rotors() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 2],
    };
    let s = a
        .state_from_relative_spin(rigid(), [0.0; 2], [0.0; 2])
        .unwrap();
    let none = RotorImpulses {
        carrier: CarrierImpulse::default(),
        motor_nms: [0.0; 2],
        aerodynamic_nms: [0.0; 2],
    };
    assert_eq!(apply_impulse(&a, &s, none, ZERO).unwrap().state, s);
    let mut enormous = s;
    enormous.axial_momentum[0] = 1e30;
    let lost = RotorImpulses {
        motor_nms: [1.0, 0.0],
        ..none
    };
    assert_eq!(
        apply_impulse(&a, &enormous, lost, ZERO).unwrap_err().code,
        "impulse-unrepresentable"
    );
    let resolved = apply_impulse(
        &a,
        &enormous,
        RotorImpulses {
            motor_nms: [0.0, 1.0],
            ..none
        },
        ZERO,
    )
    .unwrap();
    assert_eq!(resolved.state.axial_momentum[1], 1.0);
    assert_eq!(resolved.state.rigid.omega_body[0], -1.0);
    assert!(v(resolved.angular_momentum_residual).norm() < 1e-14);
    assert!(resolved.kinetic_work_residual_j.abs() < 1e-14);
    let mut moving = s;
    moving.rigid.vel_mps[0] = 1e30;
    let contact = RotorImpulses {
        carrier: CarrierImpulse {
            linear_world_ns: [1.0, 0.0, 0.0],
            ..CarrierImpulse::default()
        },
        ..none
    };
    assert_eq!(
        apply_impulse(&a, &moving, contact, ZERO).unwrap_err().code,
        "impulse-unrepresentable"
    );
    let mut spinning = s;
    spinning.rigid.omega_body[0] = 1e30;
    assert_eq!(
        apply_impulse(&a, &spinning, lost, ZERO).unwrap_err().code,
        "impulse-unrepresentable"
    );
    assert_eq!(enormous.axial_momentum, [1e30, 0.0]);
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(
            apply_impulse(
                &a,
                &s,
                RotorImpulses {
                    motor_nms: [bad, 0.0],
                    ..none
                },
                ZERO
            )
            .is_err()
        );
    }
}

#[test]
fn g0_finite_steps_refuse_erased_shaft_transfers_and_aggregate_recoil() {
    let a = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 2],
    };
    let s = a
        .state_from_relative_spin(rigid(), [0.0; 2], [0.0; 2])
        .unwrap();
    let mut enormous = s;
    enormous.axial_momentum[0] = 1e30;
    for motor in [[1.0, 0.0], [0.0, 1.0]] {
        // First: the shaft's own increment disappears. Second: shaft 1
        // resolves its kick, but aggregation with shaft 0 would erase recoil.
        let error = step_held(
            &a,
            &enormous,
            no_loads(),
            RotorTorques {
                motor_nm: motor,
                aerodynamic_nm: [0.0; 2],
            },
            ZERO,
            0.05,
        )
        .unwrap_err();
        assert_eq!(error.code, "impulse-unrepresentable");
    }
    let zero = RotorTorques {
        motor_nm: [0.0; 2],
        aerodynamic_nm: [0.0; 2],
    };
    assert!(step_held(&a, &enormous, no_loads(), zero, ZERO, 0.05).is_ok());
    let underflow = RotorTorques {
        motor_nm: [f64::from_bits(1), 0.0],
        ..zero
    };
    assert_eq!(
        step_held(&a, &s, no_loads(), underflow, ZERO, 1e-6)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    assert!(step_held(&a, &s, no_loads(), zero, ZERO, 1e-6).is_ok());
    assert!(step_held(&a, &s, no_loads(), zero, ZERO, 0.0).is_err());
    assert_eq!(enormous.axial_momentum, [1e30, 0.0]);

    let gear = geared::GearedAssembly::new(
        3.0,
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0],
        },
        [geared::GearedRotor {
            axis_body: [1.0, 0.0, 0.0],
            axial_inertia_kgm2: 0.2,
            ratio: 1.0,
        }],
    )
    .unwrap();
    let mut state = gear.state_from_speed(rigid(), 0.0, 1.0).unwrap();
    state.momentum = 1e30;
    let loads = geared::GearedLoads {
        carrier: no_loads(),
        motor_nm: 1.0,
        aerodynamic_nm: [0.0],
    };
    assert_eq!(
        geared::step_held(&gear, &state, loads, ZERO, 0.05)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    let mut calls = 0;
    let error = geared::step(&gear, &state, 0.0, 0.05, ZERO, |_, _| {
        calls += 1;
        Ok(loads)
    })
    .unwrap_err();
    assert_eq!(error.code, "impulse-unrepresentable");
    assert_eq!(
        calls, 1,
        "refused half-step must not evaluate a fictitious midpoint"
    );
    let ordinary = gear.state_from_speed(rigid(), 0.0, 0.0).unwrap();
    let tiny = geared::GearedLoads {
        motor_nm: f64::from_bits(1),
        ..loads
    };
    assert_eq!(
        geared::step_held(&gear, &ordinary, tiny, ZERO, 1e-6)
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    assert_eq!(state.momentum, 1e30);
}

#[test]
fn g1_near_balanced_finite_torques_retain_the_existing_quantized_trajectory() {
    let gear = geared::GearedAssembly::new(
        3.0,
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0],
        },
        [geared::GearedRotor {
            axis_body: [1.0, 0.0, 0.0],
            axial_inertia_kgm2: 0.2,
            ratio: 1.0,
        }],
    )
    .unwrap();
    let state = gear.state_from_speed(rigid(), 0.0, 500.0).unwrap();
    let loads = geared::GearedLoads {
        carrier: no_loads(),
        motor_nm: 1.0,
        aerodynamic_nm: [-1.0 + 1e-6],
    };
    let dt = 0.1;
    let held = geared::step_held(&gear, &state, loads, ZERO, dt).unwrap();
    let midpoint = geared::step(&gear, &state, 0.0, dt, ZERO, |_, _| Ok(loads)).unwrap();
    assert_eq!(state.momentum, 100.0);
    assert_eq!(
        held.state.momentum,
        state.momentum + dt * (loads.motor_nm + loads.aerodynamic_nm[0])
    );
    assert!(held.state.momentum > state.momentum);
    assert_eq!(held, midpoint);
    let actual = held.state.momentum - state.momentum;
    let requested = dt * (loads.motor_nm + loads.aerodynamic_nm[0]);
    assert!(
        (actual - requested).abs() > 4096.0 * f64::EPSILON * requested.abs(),
        "fixture must exercise normal addition quantization beyond the explicit-impulse budget"
    );
    let independent = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let s = independent
        .state_from_relative_spin(rigid(), [0.0], [500.0])
        .unwrap();
    let r = step_held(
        &independent,
        &s,
        no_loads(),
        RotorTorques {
            motor_nm: [loads.motor_nm],
            aerodynamic_nm: loads.aerodynamic_nm,
        },
        ZERO,
        dt,
    )
    .unwrap();
    assert_eq!(
        r.state.axial_momentum[0],
        dt.mul_add(
            loads.motor_nm + loads.aerodynamic_nm[0],
            s.axial_momentum[0]
        )
    );
}

#[path = "rotor_spine/passive.rs"]
mod passive;
