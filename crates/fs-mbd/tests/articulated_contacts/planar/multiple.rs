//! G1/G3/G4: simultaneous retained histories on shared native bodies.

#[path = "transport.rs"]
mod transport;

use super::*;
use fs_mbd::articulated_contacts::tangential::planar::multiple::{
    ArticulatedMultiTangentialConfig, ArticulatedMultiTangentialResult, ArticulatedTangentialSite,
    advance_articulated_tangential_contacts,
};
use fs_tribo::tangential::multiple::MultipleTangentialConfig;

fn controls() -> ArticulatedMultiTangentialConfig {
    let base = config();
    ArticulatedMultiTangentialConfig {
        mechanical: base.mechanical,
        coupling: MultipleTangentialConfig {
            maximum_contacts: 8,
            maximum_matrix_terms: 256,
            maximum_sweeps: 64,
            maximum_law_evaluations: 4096,
            local: MatrixTangentialConfig {
                displacement_tolerance_m: 1e-15,
                ..base.coupling
            },
            total_energy_tolerance_j: 1e-12,
        },
    }
}

fn layers() -> [ThirdBodyLayerContact; 2] {
    [
        ThirdBodyLayerContact::new(400.0, 100.0, 0.4, 300.0, 1.0).unwrap(),
        ThirdBodyLayerContact::new(700.0, 100.0, 0.7, 300.0, 1.0).unwrap(),
    ]
}

fn opposing() -> [ArticulatedPlanarTangentialContact; 2] {
    let lower = port(Vec3::new(0.0, 0.0, -0.5), Vec3::new(0.1, 0.05, 0.0));
    let mut upper = port(Vec3::new(0.0, 0.0, 0.5), Vec3::new(-0.05, 0.2, 0.0));
    upper.frame_world = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
    [lower, upper]
}

fn sites<'a, L: TangentialLaw>(
    contacts: &[ArticulatedPlanarTangentialContact],
    laws: &'a [L],
    states: &'a [L::State],
    dt: f64,
) -> Vec<ArticulatedTangentialSite<'a, L>> {
    contacts
        .iter()
        .enumerate()
        .map(|(i, &contact)| ArticulatedTangentialSite {
            contact,
            law: &laws[i],
            state: &states[i],
            drive: drive(dt),
        })
        .collect()
}

fn axes(contact: ArticulatedPlanarTangentialContact) -> [Vec3; 2] {
    [
        contact.frame_world.longitudinal(),
        contact.frame_world.lateral(),
    ]
    .map(|v| Vec3::new(v[0], v[1], v[2]))
}

// Independent rigid-body point Jacobian Gram matrix, without any ABA probe.
// Identically oriented dynamic counterparts add their inverse masses/inertias.
fn mobility(
    contacts: &[ArticulatedPlanarTangentialContact],
    inv_mass: f64,
    inv_moment: f64,
) -> Vec<f64> {
    let rows: Vec<_> = contacts
        .iter()
        .flat_map(|&contact| axes(contact).map(|axis| (axis, contact.point.point_body.cross(axis))))
        .collect();
    rows.iter()
        .flat_map(|(linear, angular)| {
            rows.iter()
                .map(move |(l, a)| inv_mass * linear.dot(*l) + inv_moment * angular.dot(*a))
        })
        .collect()
}

fn released(contacts: &[ArticulatedPlanarTangentialContact], relative: Twist, dt: f64) -> Vec<f64> {
    contacts
        .iter()
        .flat_map(|&contact| {
            let support = match contact.other {
                ContactCounterpart::PrescribedPlane { velocity_world, .. } => velocity_world,
                ContactCounterpart::Body(_) => Vec3::default(),
            };
            let speed =
                relative.linear + relative.angular.cross(contact.point.point_body) - support;
            axes(contact).map(|axis| -dt * axis.dot(speed))
        })
        .collect()
}

// Independent tiny dense elastic oracle; the production owner instead calls
// each original constitutive law through the full cross-contact coupling solve.
fn solve_dense(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        assert!(a[pivot * n + k].abs() > 1e-15);
        for j in 0..n {
            a.swap(k * n + j, pivot * n + j);
        }
        b.swap(k, pivot);
        for i in k + 1..n {
            let ratio = a[i * n + k] / a[k * n + k];
            for j in k..n {
                a[i * n + j] -= ratio * a[k * n + j];
            }
            b[i] -= ratio * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let rest: f64 = (i + 1..n).map(|j| a[i * n + j] * x[j]).sum();
        x[i] = (b[i] - rest) / a[i * n + i];
    }
    x
}

fn elastic_reference(
    w: &[f64],
    released: &[f64],
    stiffness: &[f64],
    old_elastic: &[[f64; 2]],
    dt: f64,
) -> Vec<f64> {
    let n = released.len();
    let mut a: Vec<_> = w.iter().map(|v| 0.5 * dt * dt * v).collect();
    let mut rhs = released.to_vec();
    for i in 0..n {
        a[i * n + i] += 2.0 / stiffness[i / 2];
        rhs[i] += 2.0 * old_elastic[i / 2][i % 2];
    }
    solve_dense(a, rhs)
}

fn resultant(
    contacts: &[ArticulatedPlanarTangentialContact],
    force: &[f64],
    dt: f64,
) -> (Vec3, Vec3, f64) {
    let mut impulse = Vec3::default();
    let mut torque = Vec3::default();
    let mut support_work = 0.0;
    for (i, &contact) in contacts.iter().enumerate() {
        let basis = axes(contact);
        let local = basis[0].scale(dt * force[2 * i]) + basis[1].scale(dt * force[2 * i + 1]);
        impulse = impulse + local;
        torque = torque + contact.point.point_body.cross(local);
        if let ContactCounterpart::PrescribedPlane { velocity_world, .. } = contact.other {
            support_work += local.dot(velocity_world);
        }
    }
    (impulse, torque, support_work)
}

fn audit_all<S>(result: &ArticulatedMultiTangentialResult<S>) {
    close(result.energy_balance_residual_j, 0.0);
    close(result.endpoint_velocity_residual_m_s, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());
    let mut change = 0.0;
    let mut heat = 0.0;
    for i in 0..result.contacts.len() {
        close2(result.compatibility_residual_m[i], [0.0; 2]);
        close(result.contact_work_residual_j[i], 0.0);
        assert!(result.contacts[i].dissipated_j >= 0.0);
        change += result.contact_energy_after_j[i] - result.contact_energy_before_j[i];
        heat += result.contacts[i].dissipated_j;
    }
    close(heat, result.dissipated_energy_j);
    close(
        result.kinetic_energy_change_j + change + heat,
        result.support_work_j,
    );
}

#[test]
fn opposing_contacts_share_full_body_response_and_match_a_dense_elastic_oracle() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.7, -0.4, 0.1),
        Vec3::new(0.3, -0.2, 0.15),
    );
    let contacts = opposing();
    let laws = layers();
    let old_elastic = [[0.004, -0.002], [-0.003, 0.005]];
    let history = [
        preload(&laws[0], old_elastic[0]),
        preload(&laws[1], old_elastic[1]),
    ];
    let dt = 0.02;
    let result = advance_articulated_tangential_contacts(
        &[initial],
        &sites(&contacts, &laws, &history, dt),
        controls(),
        || false,
    )
    .unwrap();
    let w = mobility(&contacts, 0.5, 1.0);
    let free = released(&contacts, initial.base.twist_body, dt);
    let expected = elastic_reference(&w, &free, &[400.0, 700.0], &old_elastic, dt);
    for (actual, target) in result.mobility_kg_inverse.iter().zip(&w) {
        close(*actual, *target);
    }
    close(w[2], 0.25); // site 0 x / site 1 x cross response
    close(w[7], -0.25); // opposing right-handed frames reverse the y cross sign
    for i in 0..2 {
        close2(
            result.contacts[i].traction_n,
            [expected[2 * i], expected[2 * i + 1]],
        );
        let k = laws[i].stiffness_n_m();
        let ignored = (k * old_elastic[i][0] + 0.5 * k * free[2 * i])
            / (1.0 + 0.25 * k * dt * dt * w[(2 * i) * 4 + 2 * i]);
        assert!((ignored - result.contacts[i].traction_n[0]).abs() > 1e-3);
        close(result.contacts[i].dissipated_j, 0.0);
    }
    let (impulse, torque, work) = resultant(&contacts, &expected, dt);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        initial.base.twist_body.linear + impulse.scale(0.5),
    );
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        initial.base.twist_body.angular + torque,
    );
    close(result.support_work_j, work);
    assert_eq!(
        result.responses[0].base_after.world_from_base,
        initial.base.world_from_base
    );
    audit_all(&result);
}

#[test]
fn three_orthogonal_interfaces_retain_cross_terms_and_permute_covariantly() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.3, 0.2, -0.1),
        Vec3::new(0.1, -0.2, 0.4),
    );
    let contacts = [
        port(Vec3::new(0.0, 0.0, 0.5), Vec3::default()),
        ArticulatedPlanarTangentialContact {
            frame_world: TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
            ..port(Vec3::new(0.5, 0.0, 0.0), Vec3::default())
        },
        ArticulatedPlanarTangentialContact {
            frame_world: TangentFrame::new([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap(),
            ..port(Vec3::new(0.0, 0.5, 0.0), Vec3::default())
        },
    ];
    let stiffness = [400.0, 600.0, 800.0];
    let laws = stiffness.map(|k| ThirdBodyLayerContact::new(k, 100.0, 0.4, 300.0, 1.0).unwrap());
    let elastic = [[0.004, -0.002], [-0.003, 0.005], [0.002, 0.003]];
    let histories = [
        preload(&laws[0], elastic[0]),
        preload(&laws[1], elastic[1]),
        preload(&laws[2], elastic[2]),
    ];
    let dt = 0.02;
    let mut baseline: Option<FreeFloatingBaseState> = None;
    let mut exact = controls();
    exact.coupling.maximum_contacts = 3;
    exact.coupling.maximum_matrix_terms = 36;
    for order in [[0, 1, 2], [2, 0, 1]] {
        let ordered = order.map(|i| contacts[i]);
        let input: Vec<_> = order
            .iter()
            .map(|&i| ArticulatedTangentialSite {
                contact: contacts[i],
                law: &laws[i],
                state: &histories[i],
                drive: drive(dt),
            })
            .collect();
        let result =
            advance_articulated_tangential_contacts(&[initial], &input, exact, || false).unwrap();
        let w = mobility(&ordered, 0.5, 1.0);
        let free = released(&ordered, initial.base.twist_body, dt);
        let expected = elastic_reference(
            &w,
            &free,
            &order.map(|i| stiffness[i]),
            &order.map(|i| elastic[i]),
            dt,
        );
        for (actual, target) in result.mobility_kg_inverse.iter().zip(w) {
            close(*actual, target);
        }
        for (i, step) in result.contacts.iter().enumerate() {
            close2(step.traction_n, [expected[2 * i], expected[2 * i + 1]]);
        }
        if let Some(reference) = baseline {
            close_vector(
                result.responses[0].base_after.twist_body.linear,
                reference.twist_body.linear,
            );
            close_vector(
                result.responses[0].base_after.twist_body.angular,
                reference.twist_body.angular,
            );
        } else {
            baseline = Some(result.responses[0].base_after);
        }
        audit_all(&result);
    }
}

#[test]
fn two_dynamic_counterparts_receive_one_coupled_impulse_and_conserve_momentum() {
    let models = [body(2.0, 1.0), body(3.0, 1.5)];
    let systems = [
        system(
            &models[0],
            Vec3::default(),
            Vec3::new(0.8, -0.3, 0.1),
            Vec3::new(0.2, -0.1, 0.0),
        ),
        system(
            &models[1],
            Vec3::default(),
            Vec3::new(-0.4, 0.2, -0.1),
            Vec3::new(-0.1, 0.15, 0.0),
        ),
    ];
    let mut contacts = opposing();
    for contact in &mut contacts {
        contact.other = ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: contact.point.point_body,
        });
    }
    let laws = layers();
    let elastic = [[0.004, -0.002], [-0.003, 0.005]];
    let histories = [preload(&laws[0], elastic[0]), preload(&laws[1], elastic[1])];
    let dt = 0.02;
    let result = advance_articulated_tangential_contacts(
        &systems,
        &sites(&contacts, &laws, &histories, dt),
        controls(),
        || false,
    )
    .unwrap();
    let relative = Twist::new(
        systems[0].base.twist_body.angular - systems[1].base.twist_body.angular,
        systems[0].base.twist_body.linear - systems[1].base.twist_body.linear,
    );
    let w = mobility(&contacts, 0.5 + 1.0 / 3.0, 1.0 + 1.0 / 1.5);
    let expected = elastic_reference(
        &w,
        &released(&contacts, relative, dt),
        &[400.0, 700.0],
        &elastic,
        dt,
    );
    for i in 0..2 {
        close2(
            result.contacts[i].traction_n,
            [expected[2 * i], expected[2 * i + 1]],
        );
    }
    let (impulse, torque, _) = resultant(&contacts, &expected, dt);
    for (i, sign, mass, moment) in [(0, 1.0, 2.0, 1.0), (1, -1.0, 3.0, 1.5)] {
        close_vector(
            result.responses[i].base_after.twist_body.linear,
            systems[i].base.twist_body.linear + impulse.scale(sign / mass),
        );
        close_vector(
            result.responses[i].base_after.twist_body.angular,
            systems[i].base.twist_body.angular + torque.scale(sign / moment),
        );
    }
    close(result.support_work_j, 0.0);
    audit_all(&result);
}

#[test]
fn nonproportional_plastic_histories_advance_together_without_reset() {
    let model = body(2.0, 1.0);
    let mut mechanics = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let laws = [
        ThirdBodyLayerContact::new(400.0, 0.5, 0.4, 300.0, 1.0).unwrap(),
        ThirdBodyLayerContact::new(700.0, 0.8, 0.7, 300.0, 1.0).unwrap(),
    ];
    let mut histories = [laws[0].initial_state(), laws[1].initial_state()];
    let mut contacts = opposing();
    let dt = 0.01;
    let mut total_heat = 0.0;
    for (iteration, support) in [
        [1.0, 0.2],
        [-0.3, 1.0],
        [-1.0, -0.4],
        [0.2, -1.0],
        [0.8, 0.7],
    ]
    .into_iter()
    .enumerate()
    {
        for (i, contact) in contacts.iter_mut().enumerate() {
            let velocity = if i == 0 {
                Vec3::new(support[0], support[1], 0.0)
            } else {
                Vec3::new(-0.4 * support[1], 0.7 * support[0], 0.0)
            };
            contact.other = ContactCounterpart::PrescribedPlane {
                point_world: contact.point.point_body,
                velocity_world: velocity,
            };
        }
        let result = advance_articulated_tangential_contacts(
            &[mechanics],
            &sites(&contacts, &laws, &histories, dt),
            controls(),
            || false,
        )
        .unwrap();
        for i in 0..2 {
            let step = &result.contacts[i];
            let old_e = histories[i].elastic_slip_m();
            let old_p = histories[i].plastic_slip_m();
            let new_e = step.state.elastic_slip_m();
            let new_p = step.state.plastic_slip_m();
            let norm = step.traction_n[0].hypot(step.traction_n[1]);
            let distance =
                dt / laws[i].plastic_resistance_n_s_m() * (norm - laws[i].yield_force_n()).max(0.0);
            for j in 0..2 {
                let plastic = if norm == 0.0 {
                    0.0
                } else {
                    distance * step.traction_n[j] / norm
                };
                close(new_p[j], old_p[j] + plastic);
                close(new_e[j], old_e[j] + result.slip_increment_m[i][j] - plastic);
                close(
                    step.traction_n[j],
                    0.5 * laws[i].stiffness_n_m() * (old_e[j] + new_e[j]),
                );
            }
            close(
                step.dissipated_j,
                laws[i].yield_force_n() * distance
                    + laws[i].plastic_resistance_n_s_m() / dt * distance * distance,
            );
            if iteration > 0 {
                assert_ne!(step.state, laws[i].initial_state());
            }
        }
        total_heat += result.dissipated_energy_j;
        audit_all(&result);
        mechanics.base = result.responses[0].base_after;
        histories = [
            result.contacts[0].state.clone(),
            result.contacts[1].state.clone(),
        ];
    }
    assert!(total_heat > 0.01);
    for state in histories {
        assert!(state.plastic_slip_m()[0] != 0.0 && state.plastic_slip_m()[1] != 0.0);
    }
}

#[test]
fn manufactured_mixed_stick_and_slide_closes_each_physical_work_port() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let mut contacts = opposing();
    let laws = [
        FrictionLaw::Coulomb {
            static_mu: 0.5,
            kinetic_mu: 0.3,
        },
        FrictionLaw::Coulomb {
            static_mu: 0.2,
            kinetic_mu: 0.1,
        },
    ];
    let histories = [laws[0].initial_state(), laws[1].initial_state()];
    let expected_force = [1.0, -0.5, 0.8, 0.6];
    let expected_slip = [0.0, 0.0, 0.002 * 0.8, 0.002 * 0.6];
    let dt = 0.02;
    let w = mobility(&contacts, 0.5, 1.0);
    let mut free = expected_slip;
    for i in 0..4 {
        for j in 0..4 {
            free[i] += 0.5 * dt * dt * w[4 * i + j] * expected_force[j];
        }
    }
    for (i, contact) in contacts.iter_mut().enumerate() {
        let basis = axes(*contact);
        contact.other = ContactCounterpart::PrescribedPlane {
            point_world: contact.point.point_body,
            velocity_world: basis[0].scale(free[2 * i] / dt) + basis[1].scale(free[2 * i + 1] / dt),
        };
    }
    let result = advance_articulated_tangential_contacts(
        &[initial],
        &sites(&contacts, &laws, &histories, dt),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.contacts[0].regime, FrictionRegime::Sticking);
    assert_eq!(result.contacts[1].regime, FrictionRegime::Sliding);
    for i in 0..2 {
        close2(
            result.contacts[i].traction_n,
            [expected_force[2 * i], expected_force[2 * i + 1]],
        );
        close2(
            result.slip_increment_m[i],
            [expected_slip[2 * i], expected_slip[2 * i + 1]],
        );
        close2(
            result.contacts[i].state.slip_m(),
            result.slip_increment_m[i],
        );
    }
    close(result.contacts[0].dissipated_j, 0.0);
    close(result.contacts[1].dissipated_j, 0.002);
    let (impulse, torque, work) = resultant(&contacts, &expected_force, dt);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        impulse.scale(0.5),
    );
    close_vector(result.responses[0].base_after.twist_body.angular, torque);
    close(result.support_work_j, work);
    audit_all(&result);
}

#[test]
fn whole_contact_set_refuses_limits_and_late_cancellation_without_partial_history() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.7, -0.4, 0.1),
        Vec3::new(0.3, -0.2, 0.15),
    );
    let contacts = opposing();
    let laws = layers();
    let histories = [
        preload(&laws[0], [0.004, -0.002]),
        preload(&laws[1], [-0.003, 0.005]),
    ];
    let saved = histories.clone();
    let dt = 0.02;
    let input = sites(&contacts, &laws, &histories, dt);
    let mut polls = 0;
    let baseline = advance_articulated_tangential_contacts(&[initial], &input, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            advance_articulated_tangential_contacts(&[initial], &input, controls(), || {
                seen += 1;
                seen >= stop
            }),
            Err(ArticulatedContactError::Cancelled)
        ));
    }
    let mut invalid = controls();
    invalid.coupling.maximum_contacts = 1;
    assert!(
        advance_articulated_tangential_contacts(&[initial], &input, invalid, || false).is_err()
    );
    invalid = controls();
    invalid.coupling.maximum_matrix_terms = 15;
    assert!(
        advance_articulated_tangential_contacts(&[initial], &input, invalid, || false).is_err()
    );
    invalid = controls();
    invalid.coupling.maximum_law_evaluations = 1;
    assert!(
        advance_articulated_tangential_contacts(&[initial], &input, invalid, || false).is_err()
    );
    invalid = controls();
    invalid.mechanical.maximum_slip_speed_m_s = 1e-8;
    assert!(
        advance_articulated_tangential_contacts(&[initial], &input, invalid, || false).is_err()
    );
    let mut wrong = sites(&contacts, &laws, &histories, dt);
    wrong[1].drive.dt_s = 0.03;
    assert!(
        advance_articulated_tangential_contacts(&[initial], &wrong, controls(), || false).is_err()
    );
    wrong = sites(&contacts, &laws, &histories, dt);
    wrong[1].drive.contact_temperature_k = Some(301.0);
    assert!(
        advance_articulated_tangential_contacts(&[initial], &wrong, controls(), || false).is_err()
    );
    let empty: [ArticulatedTangentialSite<'_, ThirdBodyLayerContact>; 0] = [];
    assert!(
        advance_articulated_tangential_contacts(&[initial], &empty, controls(), || false).is_err()
    );
    assert_eq!(histories, saved);
    let retried =
        advance_articulated_tangential_contacts(&[initial], &input, controls(), || false).unwrap();
    assert_eq!(retried.contacts, baseline.contacts);
    for (actual, expected) in retried.responses.iter().zip(&baseline.responses) {
        assert_eq!(actual.base_after, expected.base_after);
        assert_eq!(actual.joint_velocity_after, expected.joint_velocity_after);
    }
    audit_all(&retried);
}
