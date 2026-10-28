//! G1/G3/G4: actual two-axis recoil, anisotropic mobility and retained shear.

use super::*;
use fs_mbd::articulated_contacts::tangential::planar::{
    ArticulatedPlanarTangentialConfig, ArticulatedPlanarTangentialContact,
    ArticulatedPlanarTangentialResult, advance_articulated_planar_tangential,
};
use fs_mbd::articulated_contacts::tangential::{
    ArticulatedTangentialConfig, ArticulatedTangentialDrive,
};
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::{
    FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};

fn config() -> ArticulatedPlanarTangentialConfig {
    ArticulatedPlanarTangentialConfig {
        mechanical: ArticulatedTangentialConfig {
            solve: TangentialSolveConfig {
                max_iterations: 96,
                maximum_force_n: 100.0,
                force_absolute_tolerance_n: 1e-12,
                force_relative_tolerance: 1e-12,
            },
            maximum_step_s: 0.1,
            maximum_slip_speed_m_s: 10.0,
            geometry_tolerance_m: 1e-10,
            displacement_tolerance_m: 1e-11,
            velocity_tolerance_m_s: 1e-10,
            energy_tolerance_j: 1e-11,
        },
        coupling: MatrixTangentialConfig {
            maximum_evaluations: 256,
            displacement_tolerance_m: 1e-13,
            energy_tolerance_j: 1e-12,
            matrix_relative_tolerance: 1e-12,
        },
    }
}

fn drive(dt_s: f64) -> ArticulatedTangentialDrive {
    ArticulatedTangentialDrive {
        dt_s,
        normal_force_n: 10.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
    }
}

fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}

fn port(arm: Vec3, support: Vec3) -> ArticulatedPlanarTangentialContact {
    ArticulatedPlanarTangentialContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: arm,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: arm,
            velocity_world: support,
        },
        frame_world: frame(),
    }
}

fn layer(yield_force_n: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(400.0, yield_force_n, 0.4, 300.0, 1.0).unwrap()
}

fn preload<L: TangentialLaw>(law: &L, displacement: [f64; 2]) -> L::State {
    law.advance(
        &law.initial_state(),
        &TangentialStepInput {
            normal_n: 10.0,
            slip_increment_m: displacement,
            dt_s: 0.01,
            contact_temperature_k: Some(300.0),
            nominal_pressure_pa: None,
            compliance_m_per_n: 0.0,
            solve: config().mechanical.solve,
        },
    )
    .unwrap()
    .state
}

fn close2(actual: [f64; 2], expected: [f64; 2]) {
    close(actual[0], expected[0]);
    close(actual[1], expected[1]);
}

fn audit<S>(result: &ArticulatedPlanarTangentialResult<S>) {
    close(result.energy_balance_residual_j, 0.0);
    close(result.contact_work_residual_j, 0.0);
    close2(result.compatibility_residual_m, [0.0; 2]);
    close(result.endpoint_velocity_residual_m_s, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());
    assert!(result.contact.dissipated_j >= 0.0);
    close(
        result.kinetic_energy_change_j + result.contact_energy_after_j
            - result.contact_energy_before_j
            + result.contact.dissipated_j,
        result.support_work_j,
    );
}

// Independent rigid-body point mobility for isotropic rotational inertia:
// W_ij = delta_ij/m + (|r|² delta_ij-r_i r_j)/I, projected onto world x,y.
fn point_mobility(mass: f64, moment: f64, r: Vec3) -> [[f64; 2]; 2] {
    [
        [
            1.0 / mass + (r.y * r.y + r.z * r.z) / moment,
            -r.x * r.y / moment,
        ],
        [
            -r.x * r.y / moment,
            1.0 / mass + (r.x * r.x + r.z * r.z) / moment,
        ],
    ]
}

fn solve2(a: [[f64; 2]; 2], b: [f64; 2]) -> [f64; 2] {
    let determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    [
        (a[1][1] * b[0] - a[0][1] * b[1]) / determinant,
        (a[0][0] * b[1] - a[1][0] * b[0]) / determinant,
    ]
}

fn elastic_force(w: [[f64; 2]; 2], dt: f64, e: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    let c = 100.0 * dt * dt; // (k/2)*(dt²/2), k=400 N/m
    solve2(
        [
            [1.0 + c * w[0][0], c * w[0][1]],
            [c * w[1][0], 1.0 + c * w[1][1]],
        ],
        [
            400.0 * e[0] - 200.0 * dt * v[0],
            400.0 * e[1] - 200.0 * dt * v[1],
        ],
    )
}

#[test]
fn off_centre_two_axis_recoil_uses_cross_mobility_and_actual_support_work() {
    let model = body(2.0, 1.0);
    let linear = Vec3::new(2.0, -1.0, 0.0);
    let angular = Vec3::new(0.2, 0.3, -0.4);
    let initial = system(&model, Vec3::default(), linear, angular);
    let arm = Vec3::new(0.6, -0.8, 0.2);
    let support = Vec3::new(0.4, 0.2, 0.0);
    let velocity = linear + angular.cross(arm) - support;
    let dt = 0.02;
    let law = layer(100.0);
    let result = advance_articulated_planar_tangential(
        &[initial],
        port(arm, support),
        &law,
        &law.initial_state(),
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    let w = point_mobility(2.0, 1.0, arm);
    let force = elastic_force(w, dt, [0.0; 2], [velocity.x, velocity.y]);
    for i in 0..2 {
        close2(result.mobility_kg_inverse[i], w[i]);
    }
    close2(result.contact.traction_n, force);
    assert!(w[0][1].abs() > 0.4);
    let diagonal = elastic_force(
        [[w[0][0], 0.0], [0.0, w[1][1]]],
        dt,
        [0.0; 2],
        [velocity.x, velocity.y],
    );
    assert!((force[0] - diagonal[0]).abs() > 1e-3);
    let impulse = Vec3::new(dt * force[0], dt * force[1], 0.0);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        linear + impulse.scale(0.5),
    );
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        angular + arm.cross(impulse),
    );
    close(result.support_work_j, support.dot(impulse));
    close(result.contact.dissipated_j, 0.0);
    assert_eq!(
        result.responses[0].base_after.world_from_base,
        initial.base.world_from_base
    );
    audit(&result);
}

#[test]
fn retained_two_axis_preload_releases_energy_into_stationary_bodies() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let law = layer(100.0);
    let elastic = [0.01, -0.02];
    let state = preload(&law, elastic);
    let saved = state.clone();
    let arm = Vec3::new(0.6, -0.8, 0.2);
    let dt = 0.03;
    let result = advance_articulated_planar_tangential(
        &[initial],
        port(arm, Vec3::default()),
        &law,
        &state,
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    close2(
        result.contact.traction_n,
        elastic_force(point_mobility(2.0, 1.0, arm), dt, elastic, [0.0; 2]),
    );
    assert!(result.kinetic_energy_change_j > 0.0);
    assert!(result.contact.stored_change_j < 0.0);
    assert!(result.responses[0].base_after.twist_body.linear.x > 0.0);
    assert!(result.responses[0].base_after.twist_body.linear.y < 0.0);
    close(result.contact.dissipated_j, 0.0);
    assert_eq!(state, saved);
    audit(&result);
}

#[test]
fn nonproportional_layer_paths_retain_elastic_and_plastic_memory() {
    let model = body(2.0, 0.8);
    let mut mechanics = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let law = layer(0.5);
    let mut history = law.initial_state();
    let dt = 0.01;
    let mut heat = 0.0;
    for (index, support) in [
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, 0.0],
        [0.0, -1.0],
        [0.5, 0.8],
        [-0.7, 0.4],
        [0.2, -0.9],
    ]
    .into_iter()
    .enumerate()
    {
        let before = mechanics.base.twist_body.linear;
        let input = TangentialStepInput {
            normal_n: 10.0,
            slip_increment_m: [dt * (support[0] - before.x), dt * (support[1] - before.y)],
            dt_s: dt,
            contact_temperature_k: Some(300.0),
            nominal_pressure_pa: None,
            compliance_m_per_n: 0.5 * dt * dt / 2.0,
            solve: config().mechanical.solve,
        };
        // Independent consumer limit: at the centre of an isotropic free body,
        // the ORIGINAL layer owner admits this scalar compliance exactly.
        let expected = law.advance(&history, &input).unwrap();
        let result = advance_articulated_planar_tangential(
            &[mechanics],
            port(Vec3::default(), Vec3::new(support[0], support[1], 0.0)),
            &law,
            &history,
            drive(dt),
            config(),
            || false,
        )
        .unwrap();
        close2(result.contact.traction_n, expected.traction_n);
        close2(
            result.contact.state.elastic_slip_m(),
            expected.state.elastic_slip_m(),
        );
        close2(
            result.contact.state.plastic_slip_m(),
            expected.state.plastic_slip_m(),
        );
        if index == 1 {
            let reset = law.advance(&law.initial_state(), &input).unwrap();
            assert!((result.contact.traction_n[0] - reset.traction_n[0]).abs() > 0.1);
            assert!(result.contact.state.plastic_slip_m()[0] > 0.0);
            assert!(result.contact.state.plastic_slip_m()[1] > 0.0);
        }
        heat += result.contact.dissipated_j;
        audit(&result);
        mechanics.base = result.responses[0].base_after;
        history = result.contact.state;
    }
    assert!(heat > 0.01);
    assert!(history.elastic_slip_m()[0] != 0.0 && history.elastic_slip_m()[1] != 0.0);
}

#[test]
fn dynamic_pairs_exchange_full_wrench_and_rotate_covariantly() {
    let models = [body(2.0, 0.8), body(3.0, 1.2)];
    let law = layer(100.0);
    let history = preload(&law, [0.003, -0.004]);
    let arm = Vec3::new(0.4, -0.3, 0.1);
    let dt = 0.01;
    let mut baseline: Option<Vec<FreeFloatingBaseState>> = None;
    for rotation in [
        So3::identity(),
        So3::try_from_quat(Quat {
            w: 0.5_f64.sqrt(),
            x: 0.0,
            y: 0.0,
            z: 0.5_f64.sqrt(),
        })
        .unwrap(),
    ] {
        let states: Vec<_> = models
            .iter()
            .zip([Vec3::new(1.2, -0.6, 0.0), Vec3::new(-0.2, 0.4, 0.0)])
            .map(|(model, velocity)| ArticulatedContactSystem {
                model,
                base: FreeFloatingBaseState::new(
                    Se3::from_parts(rotation, Vec3::default()).unwrap(),
                    Twist::new(Vec3::default(), velocity),
                ),
                joint_position: &[],
                joint_velocity: &[],
            })
            .collect();
        let normal = rotation.rotate(Vec3::new(0.0, 0.0, 1.0)).unwrap();
        let tangent = rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let contact = ArticulatedPlanarTangentialContact {
            other: ContactCounterpart::Body(LinkContactPoint {
                system_index: 1,
                link_index: 0,
                point_body: arm,
            }),
            frame_world: TangentFrame::new(
                [normal.x, normal.y, normal.z],
                [tangent.x, tangent.y, tangent.z],
            )
            .unwrap(),
            ..port(arm, Vec3::default())
        };
        let result = advance_articulated_planar_tangential(
            &states,
            contact,
            &law,
            &history,
            drive(dt),
            config(),
            || false,
        )
        .unwrap();
        let wa = point_mobility(2.0, 0.8, arm);
        let wb = point_mobility(3.0, 1.2, arm);
        let mut w = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                w[i][j] = wa[i][j] + wb[i][j];
            }
        }
        close2(
            result.contact.traction_n,
            elastic_force(w, dt, [0.003, -0.004], [1.4, -1.0]),
        );
        for i in 0..2 {
            close2(result.mobility_kg_inverse[i], w[i]);
        }
        close(result.support_work_j, 0.0);
        audit(&result);
        if let Some(reference) = &baseline {
            for (response, expected) in result.responses.iter().zip(reference) {
                close_vector(
                    response.base_after.twist_body.linear,
                    expected.twist_body.linear,
                );
                close_vector(
                    response.base_after.twist_body.angular,
                    expected.twist_body.angular,
                );
            }
        } else {
            baseline = Some(result.responses.iter().map(|r| r.base_after).collect());
        }
    }
}

#[test]
fn anisotropic_coulomb_static_and_sliding_match_independent_kkt() {
    let model = body(2.0, 1.0);
    let arm = Vec3::new(0.6, -0.8, 0.2);
    let dt = 0.02;
    let w = point_mobility(2.0, 1.0, arm);
    let c = [
        [0.5 * dt * dt * w[0][0], 0.5 * dt * dt * w[0][1]],
        [0.5 * dt * dt * w[1][0], 0.5 * dt * dt * w[1][1]],
    ];
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    for velocity in [
        [0.001, -0.002],
        // Along the strong-mobility eigenvector: global static effort is
        // 5.8 N, below 6 N. A small isotropic probe alone can wrongly settle
        // on the distinct 4 N kinetic root even though full-matrix stick fits.
        [-0.071456, -0.053592],
        [1.0, -0.4],
    ] {
        let initial = system(
            &model,
            Vec3::default(),
            Vec3::new(velocity[0], velocity[1], 0.0),
            Vec3::default(),
        );
        let result = advance_articulated_planar_tangential(
            &[initial],
            port(arm, Vec3::default()),
            &law,
            &law.initial_state(),
            drive(dt),
            config(),
            || false,
        )
        .unwrap();
        let released = [-dt * velocity[0], -dt * velocity[1]];
        let static_force = solve2(c, released);
        if static_force[0].hypot(static_force[1]) <= 6.0 {
            assert_eq!(result.contact.regime, FrictionRegime::Sticking);
            close2(result.contact.traction_n, static_force);
            close2(result.slip_increment_m, [0.0; 2]);
            close2(
                result.relative_speed_after_m_s,
                [-velocity[0], -velocity[1]],
            );
            close(result.contact.dissipated_j, 0.0);
        } else {
            // Coulomb KKT: (C+alpha I)F = released, |F|=mu_k N.
            // This scalar multiplier oracle shares no matrix-resolvent iteration.
            let mut lo = 0.0;
            let mut hi = released[0].hypot(released[1]) / 4.0;
            for _ in 0..96 {
                let alpha = 0.5 * (lo + hi);
                let force = solve2(
                    [[c[0][0] + alpha, c[0][1]], [c[1][0], c[1][1] + alpha]],
                    released,
                );
                if force[0].hypot(force[1]) > 4.0 {
                    lo = alpha;
                } else {
                    hi = alpha;
                }
            }
            let alpha = 0.5 * (lo + hi);
            let force = solve2(
                [[c[0][0] + alpha, c[0][1]], [c[1][0], c[1][1] + alpha]],
                released,
            );
            assert_eq!(result.contact.regime, FrictionRegime::Sliding);
            close2(result.contact.traction_n, force);
            close2(
                result.slip_increment_m,
                [alpha * force[0], alpha * force[1]],
            );
            close(result.contact.dissipated_j, alpha * 16.0);
        }
        audit(&result);
    }
}

#[test]
fn internal_prismatic_pair_keeps_rank_one_mobility_and_transverse_preload() {
    let model =
        super::friction::two_link(JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap());
    let initial = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.2, -0.1, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[0.7],
    };
    let contact = ArticulatedPlanarTangentialContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 1,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        }),
        frame_world: frame(),
    };
    let law = layer(100.0);
    let history = preload(&law, [0.01, 0.02]);
    let dt = 0.02;
    let result = advance_articulated_planar_tangential(
        &[initial],
        contact,
        &law,
        &history,
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    close2(result.mobility_kg_inverse[0], [1.5, 0.0]);
    close2(result.mobility_kg_inverse[1], [0.0, 0.0]);
    let force = elastic_force([[1.5, 0.0], [0.0, 0.0]], dt, [0.01, 0.02], [0.7, 0.0]);
    close2(result.contact.traction_n, force);
    close(result.contact.state.elastic_slip_m()[1], 0.02);
    close(result.slip_increment_m[1], 0.0);
    close(
        result.responses[0].joint_velocity_after[0],
        0.7 + dt * 1.5 * force[0],
    );
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(0.2 - dt * force[0] / 2.0, -0.1, 0.0),
    );
    audit(&result);
}

#[test]
fn frictionless_two_axis_motion_retains_original_free_body_state() {
    let model = body(2.0, 0.8);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.4, -0.2, 0.1),
        Vec3::new(0.1, 0.2, 0.3),
    );
    let arm = Vec3::new(0.3, -0.4, 0.2);
    let law = FrictionLaw::Coulomb {
        static_mu: 0.0,
        kinetic_mu: 0.0,
    };
    let dt = 0.01;
    let result = advance_articulated_planar_tangential(
        &[initial],
        port(arm, Vec3::default()),
        &law,
        &law.initial_state(),
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(result.contact.traction_n, [0.0; 2]);
    assert_eq!(result.responses[0].base_after, initial.base);
    let velocity = initial.base.twist_body.linear + initial.base.twist_body.angular.cross(arm);
    close2(
        result.contact.state.slip_m(),
        [-dt * velocity.x, -dt * velocity.y],
    );
    assert_eq!(result.contact.dissipated_j, 0.0);
    audit(&result);
}

#[test]
fn late_refusal_and_cancellation_preserve_both_axes_for_exact_retry() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(1.0, -0.5, 0.0),
        Vec3::default(),
    );
    let law = layer(100.0);
    let history = preload(&law, [0.003, -0.004]);
    let saved = history.clone();
    let contact = port(Vec3::new(0.6, -0.8, 0.2), Vec3::default());
    let dt = 0.02;
    let mut polls = 0;
    let baseline = advance_articulated_planar_tangential(
        &[initial],
        contact,
        &law,
        &history,
        drive(dt),
        config(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            advance_articulated_planar_tangential(
                &[initial],
                contact,
                &law,
                &history,
                drive(dt),
                config(),
                || {
                    seen += 1;
                    seen >= stop
                }
            ),
            Err(ArticulatedContactError::Cancelled)
        ));
    }
    let mut invalid = config();
    invalid.mechanical.maximum_slip_speed_m_s = 1e-8;
    assert!(
        advance_articulated_planar_tangential(
            &[initial],
            contact,
            &law,
            &history,
            drive(dt),
            invalid,
            || false
        )
        .is_err()
    );
    invalid = config();
    invalid.coupling.maximum_evaluations = 1;
    assert!(
        advance_articulated_planar_tangential(
            &[initial],
            contact,
            &law,
            &history,
            drive(dt),
            invalid,
            || false
        )
        .is_err()
    );
    let mut wrong = contact;
    wrong.other = ContactCounterpart::PrescribedPlane {
        point_world: Vec3::default(),
        velocity_world: Vec3::default(),
    };
    assert!(
        advance_articulated_planar_tangential(
            &[initial],
            wrong,
            &law,
            &history,
            drive(dt),
            config(),
            || false
        )
        .is_err()
    );
    // Ill-conditioned projection must refuse at construction or at the actual
    // consumer; a shear coordinate cannot acquire a hidden normal component.
    if let Ok(frame) = TangentFrame::new([1.0, 1.0, 0.0], [1.0, 1.0, 1e-16]) {
        wrong = contact;
        wrong.frame_world = frame;
        assert!(
            advance_articulated_planar_tangential(
                &[initial],
                wrong,
                &law,
                &history,
                drive(dt),
                config(),
                || false
            )
            .is_err()
        );
    }
    let foreign = ThirdBodyLayerContact::new(401.0, 100.0, 0.4, 300.0, 1.0).unwrap();
    assert!(
        advance_articulated_planar_tangential(
            &[initial],
            contact,
            &foreign,
            &history,
            drive(dt),
            config(),
            || false
        )
        .is_err()
    );
    let mut wrong_temperature = drive(dt);
    wrong_temperature.contact_temperature_k = Some(301.0);
    assert!(
        advance_articulated_planar_tangential(
            &[initial],
            contact,
            &law,
            &history,
            wrong_temperature,
            config(),
            || false
        )
        .is_err()
    );
    assert_eq!(history, saved);
    let retried = advance_articulated_planar_tangential(
        &[initial],
        contact,
        &law,
        &history,
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(retried.contact, baseline.contact);
    assert_eq!(
        retried.responses[0].base_after,
        baseline.responses[0].base_after
    );
    audit(&retried);
}
