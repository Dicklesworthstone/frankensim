//! G0/G1/G3: simultaneous fluid/body/hinge motion, independent physical oracles.
use fs_flyer::addedmass::hinged::{HingedStrip, assemble_hinged_strips};
use fs_flyer::addedmass::{Strip, assemble_analytic_strip};
use fs_flyer::partitioned::coupled::{CoupledLoads, CoupledStepConfig, step};
use fs_flyer::partitioned::{HingeParams, MemoryState, PartitionedState};
use fs_flyer::spine::{Loads, RigidBody, SixDofState};
use fs_ga::{Quat, So3, Vec3};

fn config() -> CoupledStepConfig {
    CoupledStepConfig {
        max_iterations: 64,
        linear_tolerance_mps: 1e-13,
        angular_tolerance_rad_s: 1e-13,
        hinge_tolerance_rad_s: 1e-13,
    }
}
fn rest() -> PartitionedState {
    PartitionedState {
        rigid: SixDofState {
            pos_m: [0.0; 3],
            vel_mps: [0.0; 3],
            quat: [1.0, 0.0, 0.0, 0.0],
            omega_body: [0.0; 3],
        },
        hinge_q: 0.0,
        hinge_v: 0.0,
        memory: vec![],
    }
}
fn plate(pivot: [f64; 3]) -> HingedStrip {
    HingedStrip::new(
        Strip {
            name: "control",
            chord_m: 2.0,
            span_m: 3.0,
            position_m: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            span_axis: [0.0, 1.0, 0.0],
            control_coord: None,
            control_gain: 0.0,
            control_angular_gain: [0.0; 3],
        },
        pivot,
        [0.0, 1.0, 0.0],
        0,
    )
    .unwrap()
}
fn zero_loads() -> Loads {
    Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn stiff_hinge_and_body_recoil_share_the_same_midpoint_mass() {
    let body = RigidBody::new(10.0, [8.0, 9.0, 10.0]).unwrap();
    let hinge = HingeParams {
        inertia: 2.0,
        damping: 3.0,
        stiffness: 200_000.0,
    };
    let plate = plate([0.0; 3]);
    let rho = 1.2;
    let j = rho * std::f64::consts::PI * 2.0_f64.powi(4) * 3.0 / 128.0;
    let recoil = -j / (9.0 + j);
    // Exact Schur complement of the independent two-angular-coordinate mass.
    let effective = hinge.inertia + j - j * j / (9.0 + j);
    let mut state = rest();
    state.hinge_q = 0.4;
    state.hinge_v = 0.3;
    state.rigid.omega_body[1] = recoil * state.hinge_v;
    let dt = 0.05;
    let h = 0.5 * dt;
    let denom = effective + hinge.damping * h + hinge.stiffness * h * h;
    let mut energy =
        0.5 * effective * state.hinge_v.powi(2) + 0.5 * hinge.stiffness * state.hinge_q.powi(2);
    for tick in 0..64 {
        let q = state.hinge_q;
        let v = state.hinge_v;
        // Closed-form Cayley response, independent of solve_joint or the step.
        let expected_q = ((effective + hinge.damping * h - hinge.stiffness * h * h) * q
            + 2.0 * effective * h * v)
            / denom;
        let expected_v = (-2.0 * hinge.stiffness * h * q
            + (effective - hinge.damping * h - hinge.stiffness * h * h) * v)
            / denom;
        let out = step(
            &body,
            &hinge,
            &state,
            &[],
            tick as f64 * dt,
            dt,
            config(),
            |_, mid, nu| {
                Ok(CoupledLoads {
                    fluid: assemble_hinged_strips(
                        rho,
                        &[],
                        std::slice::from_ref(&plate),
                        &[mid.hinge_q],
                        nu,
                        &[mid.hinge_v],
                    )?,
                    rigid: zero_loads(),
                    hinge_torque_nm: 0.0,
                })
            },
        )
        .unwrap();
        close(out.state.hinge_q, expected_q, 2e-12);
        close(out.state.hinge_v, expected_v, 2e-11);
        close(out.state.rigid.omega_body[1], recoil * expected_v, 2e-11);
        assert!(out.iterations <= 4);
        let next_energy = 0.5 * effective * out.state.hinge_v.powi(2)
            + 0.5 * hinge.stiffness * out.state.hinge_q.powi(2);
        let loss = dt * hinge.damping * (0.5 * (v + out.state.hinge_v)).powi(2);
        close(next_energy + loss, energy, 2e-8);
        energy = next_energy;
        state = out.state;
    }
}

#[test]
fn changing_body_coordinates_preserve_free_translation_to_second_order() {
    let body = RigidBody::new(3.0, [2.0; 3]).unwrap();
    let hinge = HingeParams {
        inertia: 1.0,
        damping: 0.0,
        stiffness: 0.0,
    };
    let run = |dt: f64| {
        let mut state = rest();
        state.rigid.vel_mps = [2.0, -1.0, 0.3];
        state.rigid.omega_body = [0.4, -0.6, 0.8];
        let n = (0.8 / dt).round() as usize;
        for tick in 0..n {
            state = step(
                &body,
                &hinge,
                &state,
                &[],
                tick as f64 * dt,
                dt,
                config(),
                |_, _, nu| {
                    Ok(CoupledLoads {
                        fluid: assemble_analytic_strip(1.2, &[], 1, nu, &[0.0])?,
                        rigid: zero_loads(),
                        hinge_torque_nm: 0.0,
                    })
                },
            )
            .unwrap()
            .state;
        }
        state
            .rigid
            .vel_mps
            .iter()
            .zip([2.0, -1.0, 0.3])
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let errors = [run(0.04), run(0.02), run(0.01)];
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

fn rotation(state: &PartitionedState) -> So3 {
    let q = state.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

#[test]
fn finite_hinge_motion_has_convergent_complete_kinetic_energy_balance() {
    let body = RigidBody::new(10.0, [8.0, 9.0, 10.0]).unwrap();
    let hinge = HingeParams {
        inertia: 2.0,
        damping: 0.0,
        stiffness: 3.0,
    };
    let pivot = [-0.4, 0.0, 0.1];
    let plate = plate(pivot);
    let rho = 0.05;
    let energy = |state: &PartitionedState| {
        let strip = plate.at_angle(state.hinge_q).unwrap();
        let omega = vec(state.rigid.omega_body);
        let velocity = rotation(state)
            .inverse()
            .rotate(vec(state.rigid.vel_mps))
            .unwrap();
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let relative = axis
            .cross(vec(strip.position_m) - vec(pivot))
            .scale(state.hinge_v);
        let normal_speed =
            (velocity + omega.cross(vec(strip.position_m)) + relative).dot(vec(strip.normal));
        let angular = omega + axis.scale(state.hinge_v);
        let m = rho * std::f64::consts::PI * strip.chord_m.powi(2) * strip.span_m / 4.0;
        let fluid = 0.5 * m * normal_speed.powi(2)
            + 0.5 * (m * strip.chord_m.powi(2) / 32.0) * angular.dot(vec(strip.span_axis)).powi(2)
            + 0.5
                * (m * strip.span_m.powi(2) / 12.0)
                * angular
                    .dot(vec(strip.span_axis).cross(vec(strip.normal)))
                    .powi(2);
        0.5 * body.mass_kg * velocity.dot(velocity)
            + 0.5 * omega.dot(body.inertia_tensor().matrix().apply(omega))
            + 0.5 * hinge.inertia * state.hinge_v.powi(2)
            + 0.5 * hinge.stiffness * state.hinge_q.powi(2)
            + fluid
    };
    let mut initial = rest();
    initial.rigid.vel_mps = [1.0, 0.4, -0.2];
    initial.rigid.omega_body = [0.25, 0.2, -0.1];
    initial.hinge_q = 0.35;
    initial.hinge_v = 0.4;
    let e0 = energy(&initial);
    let run = |dt: f64| {
        let mut state = initial.clone();
        for tick in 0..(0.4 / dt).round() as usize {
            state = step(
                &body,
                &hinge,
                &state,
                &[],
                tick as f64 * dt,
                dt,
                config(),
                |_, mid, nu| {
                    Ok(CoupledLoads {
                        fluid: assemble_hinged_strips(
                            rho,
                            &[],
                            std::slice::from_ref(&plate),
                            &[mid.hinge_q],
                            nu,
                            &[mid.hinge_v],
                        )?,
                        rigid: zero_loads(),
                        hinge_torque_nm: 0.0,
                    })
                },
            )
            .unwrap()
            .state;
        }
        (energy(&state) - e0).abs()
    };
    let errors = [run(0.02), run(0.01), run(0.005)];
    for pair in errors.windows(2) {
        assert!((3.5..4.5).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert!(errors[2] < 1e-6, "{errors:?}");
}

#[test]
fn exact_memory_is_sampled_at_midpoint_and_committed_once() {
    let body = RigidBody::new(3.0, [2.0; 3]).unwrap();
    let hinge = HingeParams {
        inertia: 1.0,
        damping: 0.0,
        stiffness: 0.0,
    };
    let mut state = rest();
    state.memory = vec![MemoryState {
        tau_s: 0.02,
        value: 0.4,
    }];
    let dt = 0.01;
    let next = step(
        &body,
        &hinge,
        &state,
        &[1.0],
        0.0,
        dt,
        config(),
        |time, mid, nu| {
            close(time, dt / 2.0, 0.0);
            close(mid.memory[0].value, 1.0 - 0.6 * (-time / 0.02).exp(), 1e-14);
            Ok(CoupledLoads {
                fluid: assemble_analytic_strip(1.2, &[], 1, nu, &[0.0])?,
                rigid: zero_loads(),
                hinge_torque_nm: 0.0,
            })
        },
    )
    .unwrap();
    close(
        next.state.memory[0].value,
        1.0 - 0.6 * (-dt / 0.02).exp(),
        1e-14,
    );
    assert_eq!(state.memory[0].value, 0.4);
}

#[test]
fn refused_coupling_does_not_publish_partial_state_or_accept_duplicate_loads() {
    let body = RigidBody::new(3.0, [2.0; 3]).unwrap();
    let hinge = HingeParams {
        inertia: 1.0,
        damping: 0.0,
        stiffness: 0.0,
    };
    let state = rest();
    let saved = state.clone();
    let mut budget = config();
    budget.max_iterations = 1;
    let error = step(&body, &hinge, &state, &[], 0.0, 0.01, budget, |_, _, nu| {
        Ok(CoupledLoads {
            fluid: assemble_analytic_strip(1.2, &[], 1, nu, &[0.0])?,
            rigid: Loads {
                force_n: [1.0, 0.0, 0.0],
                moment_nm: [0.0; 3],
            },
            hinge_torque_nm: 0.0,
        })
    })
    .unwrap_err();
    assert_eq!(error.code, "coupled-step-not-converged");
    assert_eq!(state, saved);
    for bad_shape in [false, true] {
        let error = step(
            &body,
            &hinge,
            &state,
            &[],
            0.0,
            0.01,
            config(),
            |_, _, nu| {
                let mut fluid = assemble_analytic_strip(1.2, &[], 1, nu, &[0.0])?;
                if bad_shape {
                    fluid.m_added_rc.clear();
                } else {
                    fluid.q_rigid_nonaccel[0] = 1.0;
                }
                Ok(CoupledLoads {
                    fluid,
                    rigid: zero_loads(),
                    hinge_torque_nm: 0.0,
                })
            },
        )
        .unwrap_err();
        assert_eq!(
            error.code,
            if bad_shape {
                "coupled-load-shape"
            } else {
                "coupled-load-ownership"
            }
        );
    }
}
