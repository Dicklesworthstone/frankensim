//! Actual powered native bodies under original local-speed Hertz traction.
use super::*;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::hertz::driven::HertzDrivenMode;
use fs_mbd::tangential_friction::hertz::rate_dependent::RigidHertzRateContact;
use fs_mbd::tangential_friction::hertz::rate_dependent::driven::*;
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

fn patch(slope: f64) -> HertzFrictionPatch {
    HertzFrictionPatch::new(
        10.0,
        R,
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.3,
            slope_per_speed: slope,
        },
    )
    .unwrap()
}
fn contact<'a>(f: &Fixture, law: &'a HertzFrictionPatch) -> RigidHertzRateContact<'a> {
    RigidHertzRateContact {
        binding: f.binding,
        frame: f.frame,
        law,
    }
}
fn quadrature() -> HertzFrictionQuadrature {
    HertzFrictionQuadrature {
        radial_panels: 8,
        angular_panels: 8,
    }
}
fn rate_controls() -> HertzRateDrivenControls {
    HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_step_s: 0.02,
            ..transient()
        },
        rest_scaled_rate_tolerance_m_s: 2e-11,
    }
}
fn native_controls() -> RigidHertzControls {
    let mut c = controls();
    c.contact.maximum_force_n = 1e4;
    c.contact.maximum_spin_moment_nm = 1e3;
    c
}
fn centre(rate: f64, spin: bool) -> Fixture {
    let mut f = fixture();
    f.binding.arm_a_body = Vec3::ZERO;
    f.binding.arm_b_body = Vec3::ZERO;
    let pose = Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap();
    f.bodies = [
        RigidBodyState::new(
            pose,
            Vec3::new(
                if spin {
                    0.0
                } else {
                    f.properties[0].mass() * rate
                },
                0.0,
                0.0,
            ),
            Vec3::new(
                0.0,
                0.0,
                if spin {
                    f.properties[0].principal_inertia_body().z * rate
                } else {
                    0.0
                },
            ),
        )
        .unwrap(),
        RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO).unwrap(),
    ];
    f
}
fn run(
    f: &Fixture,
    law: &HertzFrictionPatch,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    dt: f64,
) -> RigidHertzRateDrivenCandidate {
    prepare_rigid_hertz_rate_driven(
        &f.bodies,
        &f.properties,
        contact(f, law),
        external,
        mode,
        dt,
        native_controls(),
        quadrature(),
        rate_controls(),
        &CancelGate::new(),
    )
    .unwrap()
}

// Independent closed-form body-coordinate reference on one side of zero.
// The dry traction and the viscous pressure moment are integrated separately;
// heat is not recovered by subtracting endpoint kinetic energies.
fn affine_segment(
    initial: f64,
    direction: f64,
    acceleration: f64,
    mobility: f64,
    dry: f64,
    viscous: f64,
    dt: f64,
) -> [f64; 3] {
    let decay = mobility * viscous;
    let terminal = (acceleration - mobility * dry * direction) / decay;
    let offset = initial - terminal;
    let phi = -(-decay * dt).exp_m1() / decay;
    let travel = terminal * dt + offset * phi;
    let square = terminal * terminal * dt
        + 2.0 * terminal * offset * phi
        + offset * offset * (-(-2.0 * decay * dt).exp_m1()) / (2.0 * decay);
    [
        terminal + offset * (-decay * dt).exp(),
        travel,
        dry * direction * travel + viscous * square,
    ]
}

#[test]
fn speed_dependent_braking_holds_or_reverses_with_actual_body_work() {
    let law = patch(0.8);
    for spin in [false, true] {
        for reverse in [false, true] {
            let initial = if spin { 3.0 } else { 0.2 };
            let f = centre(initial, spin);
            let inertia = f.properties.map(|p| {
                if spin {
                    p.principal_inertia_body().z
                } else {
                    p.mass()
                }
            });
            let mobility = inertia[0].recip() + inertia[1].recip();
            let dry = 3.0
                * if spin {
                    3.0 * core::f64::consts::PI * R / 16.0
                } else {
                    1.0
                };
            let viscous = 8.0 * if spin { 2.0 * R * R / 5.0 } else { 1.0 };
            // A below-static drive retains common body acceleration after
            // stopping. An opposing over-static drive produces a real reversal.
            let supplied = inertia[0] * mobility * dry * if reverse { -4.0 } else { 0.5 };
            let acceleration = supplied / inertia[0];
            let decay = mobility * viscous;
            let stop = (decay * initial / (mobility * dry - acceleration)).ln_1p() / decay;
            let dt = stop + 0.3;
            let before_stop =
                affine_segment(initial, 1.0, acceleration, mobility, dry, viscous, stop);
            let after_stop = if reverse {
                affine_segment(0.0, -1.0, acceleration, mobility, dry, viscous, dt - stop)
            } else {
                [0.0; 3]
            };
            let endpoint = after_stop[0];
            let travel = before_stop[1] + after_stop[1];
            let heat = before_stop[2] + after_stop[2];
            let impulse = (endpoint - initial - acceleration * dt) / mobility;
            let external = [
                if spin {
                    Wrench {
                        torque_body: Vec3::new(0.0, 0.0, supplied),
                        ..Wrench::ZERO
                    }
                } else {
                    Wrench {
                        force_world: Vec3::new(supplied, 0.0, 0.0),
                        ..Wrench::ZERO
                    }
                },
                Wrench::ZERO,
            ];
            let result = run(&f, &law, external, HertzDrivenMode::Sliding, dt);
            let m = &result.report.mechanical;
            let path = &m.driven;
            assert!(result.report.exact_invariant);
            assert!(path.arrest_time_is_exact);
            let bracket = result.report.stopping_time_bracket_s.unwrap();
            assert_eq!(bracket[0], bracket[1]);
            close(bracket[0], stop, 2e-9);
            assert_eq!(
                path.mode,
                if reverse {
                    HertzDrivenMode::Sliding
                } else {
                    HertzDrivenMode::Resting
                }
            );
            close(path.first_arrest_time_s.unwrap(), stop, 2e-9);
            close(
                if spin {
                    path.spin_rate_after_rad_s
                } else {
                    path.velocity_after_m_s[0]
                },
                endpoint,
                2e-8,
            );
            close(
                if spin {
                    path.spin_increment_rad
                } else {
                    path.slip_increment_m[0]
                },
                travel,
                3e-9,
            );
            close(
                if spin {
                    path.spin_angular_impulse_nms
                } else {
                    path.force_impulse_ns[0]
                },
                impulse,
                3e-9,
            );
            close(m.dissipated_energy_j, heat, 3e-9);
            let total = inertia[0] + inertia[1];
            let travel_a = inertia[0] * initial / total * dt
                + 0.5 * supplied / total * dt * dt
                + inertia[1] / total * travel;
            close(m.external_work_j, supplied * travel_a, 4e-9);
            close(m.kinetic_change_j + heat, supplied * travel_a, 4e-9);
            for (i, expected) in [inertia[0] * initial + supplied * dt + impulse, -impulse]
                .into_iter()
                .enumerate()
            {
                close(
                    if spin {
                        result.bodies_after[i].angular_momentum_body().z
                    } else {
                        result.bodies_after[i].linear_momentum_world().x
                    },
                    expected,
                    3e-9,
                );
            }

            // Split before the event, retaining the actual two finite bodies.
            // The continuation must resolve the same physical stop/release.
            let first_dt = 0.4 * stop;
            let first = run(&f, &law, external, HertzDrivenMode::Sliding, first_dt);
            let continued = Fixture {
                bodies: first.bodies_after,
                ..centre(initial, spin)
            };
            let second = run(
                &continued,
                &law,
                external,
                first.report.mechanical.driven.mode,
                dt - first_dt,
            );
            for i in 0..2 {
                close_vec(
                    second.bodies_after[i].linear_momentum_world(),
                    result.bodies_after[i].linear_momentum_world(),
                    3e-9,
                );
                close_vec(
                    second.bodies_after[i].angular_momentum_body(),
                    result.bodies_after[i].angular_momentum_body(),
                    3e-10,
                );
            }
            close(
                first_dt + second.report.mechanical.driven.first_arrest_time_s.unwrap(),
                stop,
                2e-9,
            );
            close(
                first.report.mechanical.dissipated_energy_j
                    + second.report.mechanical.dissipated_energy_j,
                heat,
                4e-9,
            );
            close(
                first.report.mechanical.external_work_j + second.report.mechanical.external_work_j,
                m.external_work_j,
                5e-9,
            );

            if !reverse {
                let held = Fixture {
                    bodies: result.bodies_after,
                    ..centre(initial, spin)
                };
                let hold = run(&held, &law, external, HertzDrivenMode::Resting, 0.2);
                assert_eq!(hold.report.mechanical.driven.mode, HertzDrivenMode::Resting);
                assert_eq!(hold.report.mechanical.dissipated_energy_j, 0.0);
                assert!(hold.report.mechanical.external_work_j > 0.0);
                close(
                    hold.report.mechanical.kinetic_change_j,
                    hold.report.mechanical.external_work_j,
                    4e-9,
                );
            }
        }
    }
}

#[test]
fn driven_affine_slide_spin_and_release_match_exact_terminal_motion_and_body_power() {
    let law = patch(0.8);
    for spin in [false, true] {
        for initial in [0.0, if spin { 3.0 } else { 0.2 }] {
            let f = centre(initial, spin);
            let inertias = f.properties.map(|p| {
                if spin {
                    p.principal_inertia_body().z
                } else {
                    p.mass()
                }
            });
            let mobility = inertias[0].recip() + inertias[1].recip();
            let supplied = if spin { 0.3 } else { 18.0 };
            let external = [
                if spin {
                    Wrench {
                        torque_body: Vec3::new(0.0, 0.0, supplied),
                        ..Wrench::ZERO
                    }
                } else {
                    Wrench {
                        force_world: Vec3::new(supplied, 0.0, 0.0),
                        ..Wrench::ZERO
                    }
                },
                Wrench::ZERO,
            ];
            let c = 3.0
                * if spin {
                    3.0 * core::f64::consts::PI * R / 16.0
                } else {
                    1.0
                };
            let b = 8.0 * if spin { 2.0 * R * R / 5.0 } else { 1.0 };
            let k = mobility * b;
            let terminal = (supplied / inertias[0] / mobility - c) / b;
            let dt = 0.4;
            let offset = initial - terminal;
            let phi = -(-k * dt).exp_m1() / k;
            let endpoint = terminal + offset * (-k * dt).exp();
            let travel = terminal * dt + offset * phi;
            let square = terminal * terminal * dt
                + 2.0 * terminal * offset * phi
                + offset * offset * (-(-2.0 * k * dt).exp_m1()) / (2.0 * k);
            let heat = c * travel + b * square;
            let impulse = -c * dt - b * travel;
            // COM (or common angular) motion carries source work even when
            // the relative contact coordinates are approaching terminal speed.
            let total = inertias[0] + inertias[1];
            let travel_a = inertias[0] * initial / total * dt
                + 0.5 * supplied / total * dt * dt
                + inertias[1] / total * travel;
            let result = run(
                &f,
                &law,
                external,
                if initial == 0.0 {
                    HertzDrivenMode::Resting
                } else {
                    HertzDrivenMode::Sliding
                },
                dt,
            );
            let m = &result.report.mechanical;
            let path = &m.driven;
            assert_eq!(path.mode, HertzDrivenMode::Sliding);
            assert!(result.report.exact_invariant);
            assert!(!result.report.adaptive_error_estimates_available);
            close(
                if spin {
                    path.spin_rate_after_rad_s
                } else {
                    path.velocity_after_m_s[0]
                },
                endpoint,
                2e-8,
            );
            close(
                if spin {
                    path.spin_increment_rad
                } else {
                    path.slip_increment_m[0]
                },
                travel,
                2e-9,
            );
            close(
                if spin {
                    path.spin_angular_impulse_nms
                } else {
                    path.force_impulse_ns[0]
                },
                impulse,
                2e-9,
            );
            close(m.dissipated_energy_j, heat, 3e-9);
            close(m.external_work_j, supplied * travel_a, 3e-9);
            close(m.kinetic_change_j + heat, supplied * travel_a, 4e-9);
            for (i, expected) in [inertias[0] * initial + supplied * dt + impulse, -impulse]
                .into_iter()
                .enumerate()
            {
                close(
                    if spin {
                        result.bodies_after[i].angular_momentum_body().z
                    } else {
                        result.bodies_after[i].linear_momentum_world().x
                    },
                    expected,
                    2e-9,
                );
            }
        }
    }
}

// Independent original body-coordinate RK4, including physical heat and actual
// external power. It does not use the production reduced mobility/work formula.
fn body_reference(
    f: &Fixture,
    law: &HertzFrictionPatch,
    external: [Wrench; 2],
    dt: f64,
) -> [f64; 11] {
    let arms = [f.binding.arm_a_body, f.binding.arm_b_body];
    let rhs = |y: [f64; 11]| {
        let v = [Vec3::new(y[0], y[1], 0.0), Vec3::new(y[3], y[4], 0.0)];
        let omega = [Vec3::new(0.0, 0.0, y[2]), Vec3::new(0.0, 0.0, y[5])];
        let q = v[0]
            .add(omega[0].cross(arms[0]))
            .sub(v[1].add(omega[1].cross(arms[1])));
        let spin = y[2] - y[5];
        let w = law.evaluate([q.x, q.y], spin, quadrature()).unwrap();
        let force = Vec3::new(w.force_n[0], w.force_n[1], 0.0);
        let mut d = [0.0; 11];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let acceleration = external[i]
                .force_world
                .add(force.scale(sign))
                .scale(f.properties[i].mass().recip());
            d[3 * i] = acceleration.x;
            d[3 * i + 1] = acceleration.y;
            d[3 * i + 2] = (external[i].torque_body.z
                + sign * (arms[i].cross(force).z + w.spin_moment_nm))
                / f.properties[i].principal_inertia_body().z;
            d[7] += external[i].force_world.dot(v[i]) + external[i].torque_body.dot(omega[i]);
        }
        d[6] = w.dissipated_power_w;
        d[8] = q.x;
        d[9] = q.y;
        d[10] = spin;
        d
    };
    let mut y = [0.0; 11];
    for i in 0..2 {
        y[3 * i] = f.bodies[i].linear_momentum_world().x / f.properties[i].mass();
        y[3 * i + 1] = f.bodies[i].linear_momentum_world().y / f.properties[i].mass();
        y[3 * i + 2] =
            f.bodies[i].angular_momentum_body().z / f.properties[i].principal_inertia_body().z;
    }
    let h = dt / 4096.0;
    for _ in 0..4096 {
        let a = rhs(y);
        let b = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn mixed_weakening_drive_matches_original_body_ode_work_and_world_frame() {
    let f = fixture();
    let law = HertzFrictionPatch::new(
        10.0,
        R,
        FrictionLaw::VelocityDependent {
            static_mu: 1.0,
            mu_zero: 0.8,
            slope_per_speed: -1.3,
        },
    )
    .unwrap();
    let external = [
        Wrench {
            force_world: Vec3::new(1.2, 0.5, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.004),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.7, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.002),
        },
    ];
    let dt = 0.004;
    let result = run(&f, &law, external, HertzDrivenMode::Sliding, dt);
    let expected = body_reference(&f, &law, external, dt);
    let m = &result.report.mechanical;
    let w = f.mobility();
    assert!(w[0][1] != 0.0 && w[0][2] != 0.0 && w[1][2] != 0.0);
    for i in 0..2 {
        close(
            result.bodies_after[i].linear_momentum_world().x / f.properties[i].mass(),
            expected[3 * i],
            3e-9,
        );
        close(
            result.bodies_after[i].linear_momentum_world().y / f.properties[i].mass(),
            expected[3 * i + 1],
            3e-9,
        );
        close(
            R * result.bodies_after[i].angular_momentum_body().z
                / f.properties[i].principal_inertia_body().z,
            R * expected[3 * i + 2],
            3e-9,
        );
    }
    close(m.dissipated_energy_j, expected[6], 3e-10);
    close(m.external_work_j, expected[7], 3e-10);
    close(m.driven.slip_increment_m[0], expected[8], 3e-10);
    close(m.driven.slip_increment_m[1], expected[9], 3e-10);
    close(R * m.driven.spin_increment_rad, R * expected[10], 3e-10);
    close(m.audit.total_energy_residual_j, 0.0, 2e-9);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.7).unwrap();
    let a = |v: Vec3| [v.x, v.y, v.z];
    let mut transformed = fixture();
    transformed.bodies = f.bodies.map(|body| {
        RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(body.pose().position_world()),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(body.linear_momentum_world()),
            body.angular_momentum_body(),
        )
        .unwrap()
    });
    transformed.frame = TangentFrame::new(
        a(rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
        a(rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
    )
    .unwrap();
    let rotated = run(
        &transformed,
        &law,
        external.map(|w| Wrench {
            force_world: rotation.rotate_body_to_world(w.force_world),
            ..w
        }),
        HertzDrivenMode::Sliding,
        dt,
    );
    for i in 0..2 {
        close_vec(
            rotated.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(result.bodies_after[i].linear_momentum_world()),
            3e-10,
        );
        close_vec(
            rotated.bodies_after[i].angular_momentum_body(),
            result.bodies_after[i].angular_momentum_body(),
            3e-11,
        );
    }
    close(
        rotated.report.mechanical.external_work_j,
        m.external_work_j,
        3e-10,
    );
    close(
        rotated.report.mechanical.dissipated_energy_j,
        m.dissipated_energy_j,
        3e-10,
    );
}

#[test]
fn prepared_rate_trials_and_atomic_mode_publication_keep_original_geometry_and_spectators() {
    let f = centre(0.2, false);
    let law = patch(0.8);
    let external = [
        Wrench {
            force_world: Vec3::new(18.0, 0.0, 0.0),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    let gate = CancelGate::new();
    let prepared = PreparedRigidHertzRateDriven::new(
        &f.bodies,
        &f.properties,
        contact(&f, &law),
        external,
        HertzDrivenMode::Sliding,
        0.1,
        native_controls(),
        quadrature(),
        rate_controls(),
        &gate,
    )
    .unwrap();
    let first = prepared.advance(&law, &gate).unwrap();
    let changed = prepared.advance(&patch(1.2), &gate).unwrap();
    assert_ne!(changed.bodies_after, first.bodies_after);
    assert_eq!(prepared.advance(&law, &gate).unwrap(), first);
    let wrong = HertzFrictionPatch::new(11.0, R, law.friction_law().clone()).unwrap();
    assert!(prepared.advance(&wrong, &gate).is_err());
    let mut bodies = vec![f.bodies[0], f.bodies[1], f.bodies[1]];
    let properties = vec![f.properties[0], f.properties[1], f.properties[1]];
    let before = bodies.clone();
    let mut mode = HertzDrivenMode::Sliding;
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_hertz_rate_driven(
            &mut bodies,
            &properties,
            contact(&f, &law),
            external,
            &mut mode,
            0.1,
            native_controls(),
            quadrature(),
            rate_controls(),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(bodies, before);
    assert_eq!(mode, HertzDrivenMode::Sliding);
    let brake = [
        Wrench {
            force_world: Vec3::new(1.0, 0.0, 0.0),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    let mut limited = native_controls();
    limited.contact.maximum_force_n = 0.01;
    assert!(
        apply_rigid_hertz_rate_driven(
            &mut bodies,
            &properties,
            contact(&f, &law),
            brake,
            &mut mode,
            1.0,
            limited,
            quadrature(),
            rate_controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, before);
    assert_eq!(mode, HertzDrivenMode::Sliding);
    let report = apply_rigid_hertz_rate_driven(
        &mut bodies,
        &properties,
        contact(&f, &law),
        external,
        &mut mode,
        0.1,
        native_controls(),
        quadrature(),
        rate_controls(),
        &gate,
    )
    .unwrap();
    assert_eq!([bodies[0], bodies[1]], first.bodies_after);
    assert_eq!(bodies[2], before[2]);
    assert_eq!(mode, report.mechanical.driven.mode);
}
