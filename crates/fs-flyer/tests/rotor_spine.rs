//! G0/G1/G3 full carrier/rotor momentum, work and phase integration.
use fs_flyer::spine::{
    Loads, RigidBody, SixDofState,
    rotors::{AxialRotor, RotorAssembly, RotorState, RotorTorques, step_held},
};
use fs_ga::{Quat, So3, Vec3};

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
    RigidBody {
        mass_kg: 3.0,
        inertia_kgm2: [1.0, 1.5, 2.0],
    }
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
    let mut l = v(core::array::from_fn(|i| {
        a.reduced_body.inertia_kgm2[i] * s.rigid.omega_body[i]
    }));
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
        v(body().inertia_kgm2),
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
        * v(initial.omega_body).dot(v(core::array::from_fn(|i| {
            body().inertia_kgm2[i] * initial.omega_body[i]
        })));
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
