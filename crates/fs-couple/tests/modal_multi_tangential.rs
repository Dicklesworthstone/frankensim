//! G1/G3: simultaneous retained shear against full-body and analytic references.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::multiple::{
    ModalTangentialContact, MultiTangentialModalConfig, MultiTangentialModalSystem,
};
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::planar::PlanarTangentialModalDrive;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::tangential::multiple::MultipleTangentialConfig;
use fs_tribo::{FrictionLaw, FrictionRegime, TangentialLaw};

#[path = "modal_multi_tangential/heat.rs"]
mod heat;

const OMEGAS: [[f64; 2]; 3] = [[800.0, 1100.0], [1700.0, 2200.0], [1300.0, 1900.0]];
const LEFT: [[[f64; 2]; 2]; 2] = [[[8.0, 2.0], [3.0, 6.0]], [[-3.0, 7.0], [4.0, 1.0]]];
const RIGHT: [[[f64; 2]; 2]; 2] = [[[4.0, -1.0], [2.0, 5.0]], [[5.0, 3.0], [-2.0, 6.0]]];

fn limits() -> ModalFrictionConfig {
    ModalFrictionConfig {
        max_iterations: 128,
        maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-13,
        force_relative_tolerance: 1e-12,
        sticking_tolerance_m: 2e-12,
    }
}
fn config() -> MultiTangentialModalConfig {
    MultiTangentialModalConfig {
        coupling: MultipleTangentialConfig {
            maximum_contacts: 8,
            maximum_matrix_terms: 256,
            maximum_sweeps: 128,
            maximum_law_evaluations: 65536,
            local: MatrixTangentialConfig {
                maximum_evaluations: 512,
                displacement_tolerance_m: 2e-14,
                energy_tolerance_j: 2e-12,
                matrix_relative_tolerance: 1e-12,
            },
            total_energy_tolerance_j: 1e-11,
        },
        maximum_projection_terms: 4096,
    }
}
fn network_config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 4,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    }
}
fn model(omegas: &[f64], rate: u32) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        rate,
        omegas
            .iter()
            .map(|&omega| ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::new(1.0, 0.0),
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn network(connected: bool) -> CoupledModalSystem {
    let connections = if connected {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 1,
                shapes: vec![1.0, 0.7],
            },
            right: ModalAttachment {
                component: 2,
                shapes: vec![1.0, -0.5],
            },
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        OMEGAS.iter().map(|o| model(o, 48000)).collect(),
        connections,
        network_config(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn rotate(x: [f64; 2], angle: f64) -> [f64; 2] {
    [
        angle.cos() * x[0] - angle.sin() * x[1],
        angle.sin() * x[0] + angle.cos() * x[1],
    ]
}
fn contact<L>(site: usize, law: L, angle: f64) -> ModalTangentialContact<L> {
    let rows = |shapes: [[f64; 2]; 2], component| {
        let transformed: [[f64; 2]; 2] =
            std::array::from_fn(|mode| rotate([shapes[0][mode], shapes[1][mode]], angle));
        std::array::from_fn(|axis| ModalAttachment {
            component,
            shapes: (0..2).map(|mode| transformed[mode][axis]).collect(),
        })
    };
    ModalTangentialContact {
        left: rows(LEFT[site], 0),
        right: Some(rows(RIGHT[site], site + 1)),
        law,
        limits: limits(),
    }
}
fn layer(site: usize, yield_force: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new([2e6, 3e6][site], yield_force, 0.5, 300.0, 0.01).unwrap()
}
fn layered(connected: bool, yield_force: f64) -> MultiTangentialModalSystem<ThirdBodyLayerContact> {
    MultiTangentialModalSystem::new(
        network(connected),
        (0..2)
            .map(|i| contact(i, layer(i, yield_force), 0.0))
            .collect(),
        config(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn drive(speed: [f64; 2]) -> PlanarTangentialModalDrive {
    PlanarTangentialModalDrive {
        speed_m_s: speed,
        normal_force_n: 1.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
    }
}
fn states<L: TangentialLaw>(owner: &MultiTangentialModalSystem<L>) -> Vec<ModalAcousticState> {
    owner
        .network()
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}
// Independent small dense Gaussian elimination, only for the manufactured reference.
fn solve(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> [f64; 4] {
    for k in 0..4 {
        let pivot = (k..4)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..4 {
            let ratio = a[i][k] / a[k][k];
            for j in k..4 {
                a[i][j] -= ratio * a[k][j];
            }
            b[i] -= ratio * b[k];
        }
    }
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        x[i] = (b[i] - (i + 1..4).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
fn columns() -> [[f64; 6]; 4] {
    std::array::from_fn(|i| {
        let site = i / 2;
        let axis = i % 2;
        let mut row = [0.0; 6];
        row[..2].copy_from_slice(&LEFT[site][axis]);
        row[2 * (site + 1)..2 * (site + 2)].copy_from_slice(&RIGHT[site][axis]);
        for x in &mut row[2 * (site + 1)..2 * (site + 2)] {
            *x = -*x;
        }
        row
    })
}

#[test]
fn joint_elastic_forces_and_all_body_recoils_match_independent_four_axis_solution() {
    let mut owner = layered(false, 100.0);
    let dt = owner.network().sample_period_s();
    let b = columns();
    let omega: Vec<f64> = OMEGAS.into_iter().flatten().collect();
    let response: [f64; 6] =
        std::array::from_fn(|m| (1.0 - (omega[m] * dt).cos()) / omega[m].powi(2));
    let c: [[f64; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..6).map(|m| b[i][m] * b[j][m] * response[m]).sum())
    });
    for i in 0..4 {
        for j in 0..4 {
            close(owner.compliance_m_per_n()[4 * i + j], c[i][j], 2e-20);
        }
    }
    assert!(c[0][2].abs() + c[0][3].abs() > 0.2 * c[0][0]);
    let speeds = [0.02, -0.01, -0.015, 0.025];
    let q = [1e6, 1e6, 1.5e6, 1.5e6];
    let a = std::array::from_fn(|i| {
        std::array::from_fn(|j| q[i] * c[i][j] + if i == j { 1.0 } else { 0.0 })
    });
    let rhs = std::array::from_fn(|i| q[i] * speeds[i] * dt);
    let expected = solve(a, rhs);
    let frame = owner
        .step(
            &[0.0; 6],
            &[drive([speeds[0], speeds[1]]), drive([speeds[2], speeds[3]])],
        )
        .unwrap()
        .clone();
    for i in 0..4 {
        close(frame.contacts[i / 2].traction_n[i % 2], expected[i], 2e-8);
    }
    for (m, body) in states(&owner).iter().enumerate() {
        let force: f64 = (0..4).map(|i| b[i][m] * expected[i]).sum();
        close(body.displacement_m_sqrt_kg, force * response[m], 1e-13);
        close(
            body.velocity_m_sqrt_kg_per_s,
            force * (omega[m] * dt).sin() / omega[m],
            2e-10,
        );
    }
    let mut uncoupled = a;
    for i in 0..4 {
        for j in 0..4 {
            if i / 2 != j / 2 {
                uncoupled[i][j] = 0.0;
            }
        }
    }
    let wrong = solve(uncoupled, rhs);
    assert!((0..4).any(|i| (wrong[i] - expected[i]).abs() > 1e-4));
    assert_eq!(frame.friction_dissipation_j, 0.0);
    close(owner.total_energy_j().unwrap(), frame.drive_work_j, 2e-12);
}

#[test]
fn nonproportional_histories_share_body_motion_and_preserve_total_work() {
    let mut owner = layered(true, 0.04);
    let mut supplied = 0.0;
    let mut losses = 0.0;
    let mut released = [false; 2];
    for k in 0..400 {
        let speed = match k / 100 {
            0 => [0.12, 0.0],
            1 => [0.0, 0.1],
            2 => [-0.12, 0.0],
            _ => [0.0, -0.1],
        };
        let second = rotate(speed, 0.8);
        let frame = owner
            .step(
                &[0.001, -0.002, 0.003, -0.001, 0.002, 0.0],
                &[drive(speed), drive(second)],
            )
            .unwrap();
        supplied += frame.external_work_j + frame.drive_work_j;
        losses += frame.network_dissipation_j + frame.friction_dissipation_j;
        for (i, point) in frame.contacts.iter().enumerate() {
            released[i] |= point.stored_change_j < -1e-12;
            assert!(point.dissipated_j >= 0.0);
            assert!(point.work_residual_j.abs() < 1e-11);
        }
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert!(released.into_iter().all(|yes| yes));
    assert!(losses > 1e-6);
    assert!(
        owner
            .states()
            .iter()
            .all(|s| s.plastic_slip_m().iter().all(|u| u.abs() > 1e-8))
    );
    assert!(
        owner
            .network()
            .components()
            .iter()
            .all(|m| m.states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0))
    );
    close(owner.total_energy_j().unwrap() + losses, supplied, 2e-9);
}

#[test]
fn independent_tangent_rotations_and_contact_permutation_preserve_physical_solution() {
    let angles = [0.63, -0.41];
    let mut original = layered(true, 0.04);
    let mut transformed = MultiTangentialModalSystem::new(
        network(true),
        [1, 0]
            .into_iter()
            .map(|i| contact(i, layer(i, 0.04), angles[i]))
            .collect(),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    for k in 0..160 {
        let speeds = if k < 80 {
            [[0.08, -0.04], [-0.03, 0.09]]
        } else {
            [[-0.07, 0.05], [0.1, -0.02]]
        };
        let a = original
            .step(&[0.0; 6], &speeds.map(drive))
            .unwrap()
            .clone();
        let b = transformed
            .step(
                &[0.0; 6],
                &[
                    drive(rotate(speeds[1], angles[1])),
                    drive(rotate(speeds[0], angles[0])),
                ],
            )
            .unwrap();
        for i in 0..2 {
            let expected = rotate(a.contacts[i].traction_n, angles[i]);
            for axis in 0..2 {
                close(b.contacts[1 - i].traction_n[axis], expected[axis], 3e-7);
            }
        }
        close(a.friction_dissipation_j, b.friction_dissipation_j, 2e-12);
        close(a.contact_energy_j, b.contact_energy_j, 2e-12);
    }
    for (a, b) in states(&original).iter().zip(states(&transformed)) {
        close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 2e-11);
        close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s, 2e-8);
    }
}

#[test]
fn static_and_sliding_contacts_satisfy_the_same_coupled_motion() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    let mut owner = MultiTangentialModalSystem::new(
        network(false),
        (0..2).map(|i| contact(i, law.clone(), 0.0)).collect(),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let c = owner.compliance_m_per_n().to_vec();
    let dt = owner.network().sample_period_s();
    let static_force = [0.1, -0.15, -0.12, 0.08];
    let free: [f64; 4] =
        std::array::from_fn(|i| (0..4).map(|j| c[4 * i + j] * static_force[j]).sum());
    let first = owner
        .step(
            &[0.0; 6],
            &[
                drive([free[0] / dt, free[1] / dt]),
                drive([free[2] / dt, free[3] / dt]),
            ],
        )
        .unwrap();
    for i in 0..4 {
        close(
            first.contacts[i / 2].traction_n[i % 2],
            static_force[i],
            1e-7,
        );
        assert_eq!(first.contacts[i / 2].regime, FrictionRegime::Sticking);
    }
    let frame = owner
        .step(&[0.0; 6], &[drive([0.2, 0.15]), drive([0.0, 0.0])])
        .unwrap();
    assert_eq!(frame.contacts[0].regime, FrictionRegime::Sliding);
    assert_eq!(frame.contacts[1].regime, FrictionRegime::Sticking);
    let sliding = frame.contacts[0];
    close(
        sliding.traction_n[0].hypot(sliding.traction_n[1]),
        0.4,
        1e-10,
    );
    close(
        sliding.traction_n[0] * sliding.slip_distance_m[1]
            - sliding.traction_n[1] * sliding.slip_distance_m[0],
        0.0,
        2e-12,
    );
    close(
        sliding.dissipated_j,
        sliding.traction_n[0] * sliding.slip_distance_m[0]
            + sliding.traction_n[1] * sliding.slip_distance_m[1],
        2e-12,
    );
    assert!(
        frame.contacts[1]
            .slip_distance_m
            .iter()
            .all(|u| u.abs() < 2e-12)
    );
}

#[test]
fn redundant_contact_rows_retain_two_distinct_histories_without_regularization() {
    let mut owner = MultiTangentialModalSystem::new(
        network(false),
        vec![
            contact(0, layer(0, 100.0), 0.0),
            contact(0, layer(1, 100.0), 0.0),
        ],
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let c = owner.compliance_m_per_n();
    for i in 0..2 {
        for j in 0..4 {
            assert_eq!(c[i * 4 + j], c[(i + 2) * 4 + j]);
        }
    }
    let frame = owner
        .step(&[0.0; 6], &[drive([0.01, -0.02]), drive([-0.015, 0.025])])
        .unwrap()
        .clone();
    let dt = owner.network().sample_period_s();
    close(
        frame.contacts[0].slip_distance_m[0] - frame.contacts[1].slip_distance_m[0],
        0.025 * dt,
        1e-14,
    );
    close(
        frame.contacts[0].slip_distance_m[1] - frame.contacts[1].slip_distance_m[1],
        -0.045 * dt,
        1e-14,
    );
    assert_ne!(
        owner.states()[0].elastic_slip_m(),
        owner.states()[1].elastic_slip_m()
    );
    close(owner.total_energy_j().unwrap(), frame.drive_work_j, 2e-12);
}

#[test]
fn a_late_contact_refusal_or_cancellation_preserves_all_loaded_histories_and_clock() {
    let contacts = || {
        let first = contact(0, layer(0, 0.04), 0.0);
        let mut second = contact(1, layer(1, 0.04), 0.0);
        second.limits.maximum_slip_speed_m_s = 0.02;
        vec![first, second]
    };
    let mut owner =
        MultiTangentialModalSystem::new(network(true), contacts(), config(), &CancelGate::new())
            .unwrap();
    let mut control =
        MultiTangentialModalSystem::new(network(true), contacts(), config(), &CancelGate::new())
            .unwrap();
    let initial = [drive([0.01, -0.005]), drive([-0.003, 0.007])];
    owner.step(&[0.0; 6], &initial).unwrap();
    control.step(&[0.0; 6], &initial).unwrap();
    let histories = owner.states().to_vec();
    let bodies = states(&owner);
    let frame = owner.last_frame().unwrap().clone();
    let sample = owner.network().samples_rendered();
    // The shared constitutive solve succeeds; the original physical speed bound
    // on the second contact refuses only after the complete body update stages.
    assert!(
        owner
            .step(&[0.0; 6], &[drive([0.01, 0.0]), drive([0.2, 0.1])])
            .is_err()
    );
    let mut bad = initial;
    bad[1].contact_temperature_k = Some(301.0);
    assert!(owner.step(&[0.0; 6], &bad).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(owner.step_under_gate(&[0.0; 6], &initial, &gate).is_err());
    assert_eq!(owner.states(), histories);
    assert_eq!(states(&owner), bodies);
    assert_eq!(owner.last_frame(), Some(&frame));
    assert_eq!(owner.network().samples_rendered(), sample);
    let retried = owner.step(&[0.0; 6], &initial).unwrap().clone();
    let expected = control.step(&[0.0; 6], &initial).unwrap().clone();
    assert_eq!(retried, expected);
    assert_eq!(owner.states(), control.states());
    assert_eq!(states(&owner), states(&control));
}

#[test]
fn two_contact_elastic_frequency_converges_to_the_independent_continuum_solution() {
    let duration = 0.008;
    let mut errors = Vec::new();
    for rate in [8000, 16000, 32000] {
        let mut body = model(&[300.0], rate);
        body.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 0.1,
        }])
        .unwrap();
        let network =
            CoupledModalSystem::new(vec![body], vec![], network_config(), &CancelGate::new())
                .unwrap();
        let contacts = [40000.0, 90000.0]
            .into_iter()
            .map(|stiffness| ModalTangentialContact {
                left: [
                    ModalAttachment {
                        component: 0,
                        shapes: vec![1.0],
                    },
                    ModalAttachment {
                        component: 0,
                        shapes: vec![0.0],
                    },
                ],
                right: None,
                law: ThirdBodyLayerContact::new(stiffness, 100.0, 0.5, 300.0, 0.01).unwrap(),
                limits: limits(),
            })
            .collect();
        let mut owner =
            MultiTangentialModalSystem::new(network, contacts, config(), &CancelGate::new())
                .unwrap();
        for _ in 0..(duration * rate as f64) as usize {
            let f = owner
                .step(&[0.0], &[drive([0.0; 2]), drive([0.0; 2])])
                .unwrap();
            assert_eq!(f.friction_dissipation_j, 0.0);
        }
        let omega = (300.0_f64.powi(2) + 40000.0 + 90000.0).sqrt();
        let state = states(&owner)[0];
        let expected_q = 0.1 / omega * (omega * duration).sin();
        let expected_v = 0.1 * (omega * duration).cos();
        errors.push(
            ((state.displacement_m_sqrt_kg - expected_q) * omega)
                .hypot(state.velocity_m_sqrt_kg_per_s - expected_v),
        );
        close(owner.total_energy_j().unwrap(), 0.005, 1e-10);
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}
