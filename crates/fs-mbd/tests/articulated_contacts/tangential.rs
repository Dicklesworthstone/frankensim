//! G1/G3/G4: actual articulated recoil, retained shear and physical work.

use super::*;
use fs_mbd::articulated_contacts::tangential::{
    ArticulatedTangentialConfig, ArticulatedTangentialContact, ArticulatedTangentialDrive,
    ArticulatedTangentialResult, advance_articulated_tangential,
};
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn config() -> ArticulatedTangentialConfig {
    ArticulatedTangentialConfig {
        solve: TangentialSolveConfig {
            max_iterations: 80,
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

fn layer(yield_force: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(400.0, yield_force, 0.4, 300.0, 1.0).unwrap()
}

fn contact(arm: Vec3, support_speed: f64) -> ArticulatedTangentialContact {
    ArticulatedTangentialContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: arm,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: arm,
            velocity_world: Vec3::new(support_speed, 0.0, 0.0),
        },
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
    }
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
            solve: config().solve,
        },
    )
    .unwrap()
    .state
}

fn audit<S>(result: &ArticulatedTangentialResult<S>) {
    close(result.energy_balance_residual_j, 0.0);
    close(result.contact_work_residual_j, 0.0);
    close(result.compatibility_residual_m, 0.0);
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

#[test]
fn elastic_shear_uses_exact_off_centre_mobility_and_midpoint_work() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, 0.0, 0.0),
        Vec3::default(),
    );
    let law = layer(100.0);
    let dt = 0.02;
    let result = advance_articulated_tangential(
        &[initial],
        contact(Vec3::new(0.0, -1.0, 0.0), 0.0),
        &law,
        &law.initial_state(),
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    // Independent point mobility: 1/m + r^2/I = 1.5. The elastic midpoint
    // load line is F = (k/2)*(-v*dt-C*F), C=dt^2*W/2.
    let w = 1.5;
    let c = 0.5 * dt * dt * w;
    let force = -0.5 * 400.0 * 3.0 * dt / (1.0 + 0.5 * 400.0 * c);
    close(result.mobility_kg_inverse, w);
    close(result.compliance_m_per_n, c);
    close(result.contact.traction_n[0], force);
    close(result.relative_speed_after_m_s, 3.0 + w * dt * force);
    close(
        result.responses[0].base_after.twist_body.linear.x,
        3.0 + dt * force / 2.0,
    );
    close(
        result.responses[0].base_after.twist_body.angular.z,
        dt * force,
    );
    close(
        result.contact.state.elastic_slip_m()[0],
        -3.0 * dt - c * force,
    );
    close(result.contact.dissipated_j, 0.0);
    assert_eq!(
        result.responses[0].base_after.world_from_base,
        initial.base.world_from_base
    );
    audit(&result);
}

#[test]
fn retained_shear_releases_energy_into_initially_stationary_bodies() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let law = layer(100.0);
    let state = preload(&law, [0.01, 0.0]);
    let before = state.clone();
    let dt = 0.02;
    let result = advance_articulated_tangential(
        &[initial],
        contact(Vec3::default(), 0.0),
        &law,
        &state,
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    let c = 0.5 * dt * dt / 2.0;
    let force = 400.0 * 0.01 / (1.0 + 200.0 * c);
    close(result.contact.traction_n[0], force);
    close(
        result.responses[0].base_after.twist_body.linear.x,
        dt * force / 2.0,
    );
    assert!(result.kinetic_energy_change_j > 0.0);
    assert!(result.contact.stored_change_j < 0.0);
    close(result.contact.dissipated_j, 0.0);
    assert_eq!(state, before);
    audit(&result);
}

#[test]
fn prismatic_and_revolute_recoil_retain_the_actual_joint_mobility() {
    let law = layer(100.0);
    let dt = 0.02;
    let model =
        super::friction::two_link(JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap());
    let initial = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.3, 0.0, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[0.7],
    };
    let mut port = contact(Vec3::default(), 0.0);
    port.point.link_index = 1;
    let result = advance_articulated_tangential(
        &[initial],
        port,
        &law,
        &law.initial_state(),
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    let force = -200.0 * dt / (1.0 + 200.0 * 0.5 * dt * dt);
    close(result.mobility_kg_inverse, 1.0); // child mass, not total tree mass
    close(result.contact.traction_n[0], force);
    close(result.responses[0].base_after.twist_body.linear.x, 0.3);
    close(
        result.responses[0].joint_velocity_after[0],
        0.7 + dt * force,
    );
    audit(&result);

    let model =
        super::friction::two_link(JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap());
    let initial = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, 0.1, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[2.0],
    };
    let mut port = contact(Vec3::new(1.0, 0.0, 0.0), 0.0);
    port.point.link_index = 1;
    port.tangent_world = Vec3::new(0.0, 1.0, 0.0);
    let result = advance_articulated_tangential(
        &[initial],
        port,
        &law,
        &law.initial_state(),
        drive(dt),
        config(),
        || false,
    )
    .unwrap();
    let w = 1.0 / 3.0 + 1.0 / 0.5;
    let force = -200.0 * 2.1 * dt / (1.0 + 200.0 * 0.5 * dt * dt * w);
    close(result.mobility_kg_inverse, w);
    close(
        result.responses[0].base_after.twist_body.linear.y,
        0.1 + dt * force / 3.0,
    );
    close(result.responses[0].base_after.twist_body.angular.z, 0.0);
    close(
        result.responses[0].joint_velocity_after[0],
        2.0 + dt * force / 0.5,
    );
    audit(&result);
}

#[test]
fn dynamic_counterparts_conserve_momenta_and_rotate_covariantly() {
    let models = [body(2.0, 0.8), body(3.0, 1.2)];
    let law = layer(100.0);
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
            .zip([1.2, -0.4])
            .map(|(model, v)| ArticulatedContactSystem {
                model,
                base: FreeFloatingBaseState::new(
                    Se3::from_parts(rotation, Vec3::default()).unwrap(),
                    Twist::new(Vec3::default(), Vec3::new(v, 0.0, 0.0)),
                ),
                joint_position: &[],
                joint_velocity: &[],
            })
            .collect();
        let arm = Vec3::new(0.0, -0.3, 0.0);
        let port = ArticulatedTangentialContact {
            other: ContactCounterpart::Body(LinkContactPoint {
                system_index: 1,
                link_index: 0,
                point_body: arm,
            }),
            tangent_world: rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            ..contact(arm, 0.0)
        };
        let result = advance_articulated_tangential(
            &states,
            port,
            &law,
            &law.initial_state(),
            drive(dt),
            config(),
            || false,
        )
        .unwrap();
        let w = 1.0 / 2.0 + 1.0 / 3.0 + 0.09 * (1.0 / 0.8 + 1.0 / 1.2);
        let force = -200.0 * 1.6 * dt / (1.0 + 200.0 * 0.5 * dt * dt * w);
        close(result.mobility_kg_inverse, w);
        close(result.contact.traction_n[0], force);
        close(
            result.responses[0].base_after.twist_body.linear.x,
            1.2 + dt * force / 2.0,
        );
        close(
            result.responses[1].base_after.twist_body.linear.x,
            -0.4 - dt * force / 3.0,
        );
        if let Some(reference) = &baseline {
            for (actual, reference) in result.responses.iter().zip(reference) {
                close_vector(
                    actual.base_after.twist_body.linear,
                    reference.twist_body.linear,
                );
                close_vector(
                    actual.base_after.twist_body.angular,
                    reference.twist_body.angular,
                );
            }
        } else {
            baseline = Some(result.responses.iter().map(|r| r.base_after).collect());
        }
        close(result.support_work_j, 0.0);
        audit(&result);
    }
}

#[test]
fn actual_mindlin_backbone_and_reversal_memory_drive_the_body() {
    let law = MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: 0.01,
        shear_moduli_pa: [3e5, 2e5],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: 0.5,
        maximum_history_intervals: 32,
    })
    .unwrap();
    let model = body(2.0, 1.0);
    let mut base = FreeFloatingBaseState::new(Se3::identity(), Twist::zero());
    let mut state = law.initial_state();
    let dt = 0.001;
    let speed = 0.4 * law.gross_slip_displacement_m() / dt;
    let mut total_heat = 0.0;
    for (index, sign) in [1.0, -1.0, 1.0].into_iter().enumerate() {
        let body = ArticulatedContactSystem {
            model: &model,
            base,
            joint_position: &[],
            joint_velocity: &[],
        };
        let result = advance_articulated_tangential(
            &[body],
            contact(Vec3::default(), sign * speed),
            &law,
            &state,
            drive(dt),
            config(),
            || false,
        )
        .unwrap();
        if index == 0 {
            let u = result.slip_increment_m;
            let d = law.gross_slip_displacement_m();
            let endpoint_force = law.gross_slip_capacity_n() * (1.0 - (1.0 - u / d).powf(1.5));
            close(result.contact.state.endpoint_force_n(), endpoint_force);
            assert!(result.contact.traction_n[0] < endpoint_force);
        } else {
            assert!(result.contact.state.history_intervals() > 1);
        }
        assert!(result.contact.dissipated_j > 0.0);
        close(
            result.support_work_j,
            result.contact.traction_n[0] * sign * speed * dt,
        );
        audit(&result);
        total_heat += result.contact.dissipated_j;
        base = result.responses[0].base_after;
        state = result.contact.state;
    }
    assert!(total_heat > 0.0);
}

#[test]
fn dry_sliding_heat_and_late_refusal_preserve_all_original_states() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::default(),
    );
    let law = FrictionLaw::Coulomb {
        static_mu: 0.03,
        kinetic_mu: 0.03,
    };
    let state = law.initial_state();
    let saved = state.clone();
    let mut polls = 0;
    let result = advance_articulated_tangential(
        &[initial],
        contact(Vec3::default(), 0.0),
        &law,
        &state,
        drive(0.1),
        config(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    close(result.contact.traction_n[0], -0.3);
    close(result.relative_speed_after_m_s, 0.985);
    close(result.contact.dissipated_j, 0.3 * 0.1 * 0.9925);
    audit(&result);
    let mut count = 0;
    assert!(matches!(
        advance_articulated_tangential(
            &[initial],
            contact(Vec3::default(), 0.0),
            &law,
            &state,
            drive(0.1),
            config(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(ArticulatedContactError::Cancelled)
    ));
    let mut too_small = config();
    too_small.maximum_slip_speed_m_s = 0.1;
    assert!(
        advance_articulated_tangential(
            &[initial],
            contact(Vec3::default(), 0.0),
            &law,
            &state,
            drive(0.1),
            too_small,
            || false
        )
        .is_err()
    );
    assert_eq!(state, saved);
    let retry = advance_articulated_tangential(
        &[initial],
        contact(Vec3::default(), 0.0),
        &law,
        &state,
        drive(0.1),
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(result.responses, retry.responses);
    assert_eq!(result.contact, retry.contact);

    let elastic = layer(100.0);
    let transverse = preload(&elastic, [0.0, 0.01]);
    assert!(
        advance_articulated_tangential(
            &[initial],
            contact(Vec3::default(), 0.0),
            &elastic,
            &transverse,
            drive(0.01),
            config(),
            || false
        )
        .is_err()
    );
    let invalid = ArticulatedTangentialContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(0.1, 0.0, 0.0),
            velocity_world: Vec3::default(),
        },
        ..contact(Vec3::default(), 0.0)
    };
    assert!(
        advance_articulated_tangential(
            &[initial],
            invalid,
            &law,
            &state,
            drive(0.01),
            config(),
            || false
        )
        .is_err()
    );
}
