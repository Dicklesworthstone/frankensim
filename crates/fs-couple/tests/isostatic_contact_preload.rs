//! G1/G3 physical support balance on the existing free-coordinate/contact owners.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::multiple::{
    MultiContactConfig, MultiContactModalSystem,
};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn budget() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 16,
        max_connections: 8,
        max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 1000.0,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-10,
        energy_relative_tolerance: 1e-8,
    }
}
fn contact_budget() -> ModalContactConfig {
    ModalContactConfig {
        max_iterations: 96,
        maximum_force_n: 1e4,
        maximum_penetration_m: 0.05,
        force_absolute_tolerance_n: 1e-10,
        force_relative_tolerance: 1e-11,
    }
}
fn set_budget() -> MultiContactConfig {
    MultiContactConfig {
        max_contacts: 8,
        max_sweeps: 64,
        max_setup_terms: 100_000,
    }
}
fn free_coordinates(count: usize) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_free_coordinates(
        48000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::new(0.0, 0.0),
            };
            count
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn mass(kg: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_free_mass(
        48000,
        kg,
        0.0,
        0.0,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn scalar(component: usize, shape: f64) -> ModalAttachment {
    ModalAttachment {
        component,
        shapes: vec![shape],
    }
}
fn point(
    left: ModalAttachment,
    right: ModalAttachment,
    gap: f64,
    k: f64,
    weight: f64,
    alpha: f64,
) -> ModalContact {
    ModalContact {
        left,
        right,
        law: Obstacle::new(
            vec![-1.0],
            1,
            1,
            vec![gap],
            vec![weight],
            k,
            alpha,
            "authored isostatic support".into(),
        )
        .unwrap()
        .with_internal_loss(0.2)
        .unwrap(),
    }
}
fn support(column: &[f64], gap: f64, k: f64, weight: f64, alpha: f64) -> ModalContact {
    point(
        ModalAttachment {
            component: 0,
            shapes: column.to_vec(),
        },
        ModalAttachment {
            component: 0,
            shapes: vec![0.0; column.len()],
        },
        gap,
        k,
        weight,
        alpha,
    )
}
fn states(s: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states())
        .copied()
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tolerance {tolerance:e})"
    );
}
fn single_component(count: usize, config: ModalCouplingConfig) -> CoupledModalSystem {
    CoupledModalSystem::new(
        vec![free_coordinates(count)],
        vec![],
        config,
        &CancelGate::new(),
    )
    .unwrap()
}

#[test]
fn g1_three_supports_balance_force_and_torque_and_continue_on_original_contact_dynamics() {
    // Fixed small-rotation normal maps, q=(sqrt(M) z,sqrt(Ix) theta_x,sqrt(Iy) theta_y).
    // The fixture does not model finite rigid-body rotation or contact discovery.
    let roots = [0.6_f64.sqrt(), 0.004_f64.sqrt(), 0.005_f64.sqrt()];
    let xy = [[-0.035, -0.022], [0.041, -0.017], [0.006, 0.046]];
    let required = [3.0_f64, 5.0, 7.0];
    let gaps = [0.0001, 0.0002, 0.00015];
    let laws = [(3e6, 0.8, 1.5), (9e6, 0.6, 2.0), (1e5, 1.2, 1.0)];
    let mut contacts = Vec::new();
    let mut external = [0.0; 3];
    let mut closures = [0.0; 3];
    let mut expected_energy = 0.0;
    for i in 0..3 {
        let b = [1.0 / roots[0], xy[i][1] / roots[1], -xy[i][0] / roots[2]];
        for j in 0..3 {
            external[j] += b[j] * required[i];
        }
        let (k, weight, alpha) = laws[i];
        let penetration = (required[i] / (k * weight)).powf(1.0 / alpha);
        closures[i] = gaps[i] + penetration;
        expected_energy += required[i] * penetration / (alpha + 1.0);
        contacts.push((support(&b, gaps[i], k, weight, alpha), contact_budget()));
    }
    // Independent intersection of the three displaced support planes.
    let dx = [xy[1][0] - xy[0][0], xy[2][0] - xy[0][0]];
    let dy = [xy[1][1] - xy[0][1], xy[2][1] - xy[0][1]];
    let dz = [closures[1] - closures[0], closures[2] - closures[0]];
    let determinant = dx[0] * dy[1] - dy[0] * dx[1];
    let theta_x = (dx[0] * dz[1] - dz[0] * dx[1]) / determinant;
    let theta_y = (dy[0] * dz[1] - dz[0] * dy[1]) / determinant;
    let z = closures[0] - xy[0][1] * theta_x + xy[0][0] * theta_y;
    let expected = [z * roots[0], theta_x * roots[1], theta_y * roots[2]];
    let mut first = None;
    for reverse in [false, true] {
        let basis = if reverse {
            vec![
                vec![0.6, 0.8, 0.0],
                vec![-0.8, 0.6, 0.0],
                vec![0.0, 0.0, -1.0],
            ]
        } else {
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ]
        };
        let mut ordered = contacts.clone();
        if reverse {
            ordered.reverse();
        }
        let mut s = single_component(3, budget());
        let energy = s
            .initialize_contact_equilibrium_with_gauges(
                &external,
                &ordered,
                &basis,
                set_budget(),
                &CancelGate::new(),
            )
            .unwrap();
        for (actual, wanted) in states(&s).iter().zip(expected) {
            close(actual.displacement_m_sqrt_kg, wanted, 3e-12);
            assert_eq!(actual.velocity_m_sqrt_kg_per_s, 0.0);
        }
        close(energy, expected_energy, 3e-12);
        if first.is_none() {
            first = Some((s, energy));
        }
    }
    let (s, energy) = first.unwrap();
    let mut dynamic =
        MultiContactModalSystem::new(s, contacts, set_budget(), &CancelGate::new()).unwrap();
    for _ in 0..32 {
        dynamic.step(&external).unwrap();
    }
    for model in dynamic.components() {
        for state in model.states() {
            assert!(state.velocity_m_sqrt_kg_per_s.abs() < 2e-8);
        }
    }
    let frame = dynamic.last_frame().unwrap();
    for (contact, reaction) in frame.contacts.iter().zip(required) {
        close(contact.normal_force_n, reaction, 2e-8);
    }
    for _ in 0..48 {
        dynamic.step(&[0.0; 3]).unwrap();
    }
    assert!(dynamic.components()[0].states()[0].velocity_m_sqrt_kg_per_s < 0.0);
    assert!(dynamic.total_energy_j().unwrap() <= energy + 1e-9);
}

#[test]
fn g1_stacked_unequal_assemblies_keep_rest_lengths_and_reciprocal_elastic_receiver() {
    let masses = [0.04_f64, 0.09, 0.07, 0.13];
    let roots = masses.map(f64::sqrt);
    let k = [1200.0, 1700.0];
    let rests = [0.0012, -0.0003];
    let loads = [0.2_f64, 0.7, 0.4, 0.6];
    let external = [
        loads[0] / roots[0],
        loads[1] / roots[1],
        loads[2] / roots[2],
        loads[3] / roots[3],
        0.3,
    ];
    let contacts = vec![
        (
            point(
                scalar(1, 1.0 / roots[1]),
                scalar(4, 2.0),
                0.0002,
                3e6,
                0.8,
                1.5,
            ),
            contact_budget(),
        ),
        (
            point(
                scalar(3, 1.0 / roots[3]),
                scalar(1, 1.0 / roots[1]),
                0.00015,
                8e6,
                0.7,
                2.0,
            ),
            contact_budget(),
        ),
    ];
    let norm_a = (masses[0] + masses[1]).sqrt();
    let norm_b = (masses[2] + masses[3]).sqrt();
    let a = [roots[0] / norm_a, roots[1] / norm_a, 0.0, 0.0, 0.0];
    let b = [0.0, 0.0, roots[2] / norm_b, roots[3] / norm_b, 0.0];
    let r = [1.9_f64, 1.0_f64];
    let penetration = [
        (r[0] / (3e6 * 0.8)).powf(2.0 / 3.0),
        (r[1] / (8e6 * 0.7)).sqrt(),
    ];
    let receiver = (external[4] + 2.0 * r[0]) / 800.0_f64.powi(2);
    let x1 = 2.0 * receiver + 0.0002 + penetration[0];
    let x3 = x1 + 0.00015 + penetration[1];
    let expected = [
        (x1 + rests[0] + loads[0] / k[0]) * roots[0],
        x1 * roots[1],
        (x3 + rests[1] + loads[2] / k[1]) * roots[2],
        x3 * roots[3],
        receiver,
    ];
    let expected_energy = 0.5 * loads[0] * loads[0] / k[0]
        + 0.5 * loads[2] * loads[2] / k[1]
        + 0.5 * 800.0_f64.powi(2) * receiver * receiver
        + r[0] * penetration[0] / 2.5
        + r[1] * penetration[1] / 3.0;
    for rotate in [false, true] {
        let receiver_model = ModalAcousticTimeModel::try_new(
            48000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 800.0,
                damping_ratio: 0.02,
                pressure_per_modal_velocity: C64::new(1.0, 0.0),
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        let mut models: Vec<_> = masses.into_iter().map(mass).collect();
        models.push(receiver_model);
        let links: Vec<_> = (0..2)
            .map(|i| ModalConnection {
                left: scalar(2 * i, 1.0 / roots[2 * i]),
                right: scalar(2 * i + 1, 1.0 / roots[2 * i + 1]),
                stiffness_n_m: k[i],
                damping_n_s_m: 2.0,
                rest_extension_m: rests[i],
            })
            .collect();
        let mut s = CoupledModalSystem::new(models, links, budget(), &CancelGate::new()).unwrap();
        let basis = if rotate {
            vec![
                a.iter().zip(b).map(|(a, b)| 0.6 * a + 0.8 * b).collect(),
                a.iter().zip(b).map(|(a, b)| -0.8 * a + 0.6 * b).collect(),
            ]
        } else {
            vec![a.to_vec(), b.to_vec()]
        };
        let energy = s
            .initialize_contact_equilibrium_with_gauges(
                &external,
                &contacts,
                &basis,
                set_budget(),
                &CancelGate::new(),
            )
            .unwrap();
        let actual = states(&s);
        for (actual, expected) in actual.iter().zip(expected) {
            close(actual.displacement_m_sqrt_kg, expected, 2e-12);
        }
        close(energy, expected_energy, 2e-12);
        // Direct original per-mass balances, not only the reduced support map.
        let x: Vec<_> = actual[..4]
            .iter()
            .zip(roots)
            .map(|(q, r)| q.displacement_m_sqrt_kg / r)
            .collect();
        let springs = [
            k[0] * (x[0] - x[1] - rests[0]),
            k[1] * (x[2] - x[3] - rests[1]),
        ];
        let reactions = [
            3e6 * 0.8 * (x[1] - 2.0 * actual[4].displacement_m_sqrt_kg - 0.0002).powf(1.5),
            8e6 * 0.7 * (x[3] - x[1] - 0.00015).powi(2),
        ];
        for residual in [
            loads[0] - springs[0],
            loads[1] + springs[0] - reactions[0] + reactions[1],
            loads[2] - springs[1],
            loads[3] + springs[1] - reactions[1],
            external[4] + 2.0 * reactions[0] - 800.0_f64.powi(2) * actual[4].displacement_m_sqrt_kg,
        ] {
            close(residual, 0.0, 2e-9);
        }
    }
}

#[test]
fn g1_zero_load_has_a_feasible_open_or_touching_pose_without_contact_storage() {
    let mut s = single_component(2, budget());
    let contacts = vec![
        (support(&[1.0, 1.0], 0.001, 1e6, 1.0, 1.5), contact_budget()),
        (
            support(&[1.0, -1.0], -0.0004, 2e6, 0.8, 2.0),
            contact_budget(),
        ),
    ];
    let energy = s
        .initialize_contact_equilibrium_with_gauges(
            &[0.0, 0.0],
            &contacts,
            &[vec![0.6, 0.8], vec![-0.8, 0.6]],
            set_budget(),
            &CancelGate::new(),
        )
        .unwrap();
    let q = states(&s);
    close(
        q[0].displacement_m_sqrt_kg + q[1].displacement_m_sqrt_kg,
        0.0,
        1e-17,
    );
    close(
        q[0].displacement_m_sqrt_kg - q[1].displacement_m_sqrt_kg,
        -0.0004,
        1e-17,
    );
    close(energy, 0.0, 1e-24);
    // One positive reaction and one zero reaction also retain the open closure.
    let mut s = single_component(2, budget());
    let contacts = vec![
        (
            support(&[1.0, 0.0], 0.0001, 1e6, 1.0, 2.0),
            contact_budget(),
        ),
        (support(&[0.0, 1.0], 0.001, 2e6, 1.0, 1.5), contact_budget()),
    ];
    s.initialize_contact_equilibrium_with_gauges(
        &[1.0, 0.0],
        &contacts,
        &[vec![1.0, 0.0], vec![0.0, 1.0]],
        set_budget(),
        &CancelGate::new(),
    )
    .unwrap();
    close(states(&s)[1].displacement_m_sqrt_kg, 0.0, 0.0);
}

#[test]
fn g3_one_support_uses_the_original_owner_bit_for_bit() {
    let stop = support(&[1.0], 0.0003, 1e6, 0.7, 1.5);
    let mut original = single_component(1, budget());
    let mut set = single_component(1, budget());
    let a = original
        .initialize_contact_equilibrium_with_gauge(
            &[1.3],
            &stop,
            contact_budget(),
            &[1.0],
            &CancelGate::new(),
        )
        .unwrap();
    let b = set
        .initialize_contact_equilibrium_with_gauges(
            &[1.3],
            &[(stop, contact_budget())],
            &[vec![1.0]],
            set_budget(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(states(&original), states(&set));
    assert_eq!(a, b);
}

#[test]
fn g0_tensile_redundant_incomplete_and_unresolved_support_sets_refuse_atomically() {
    let independent = vec![
        (support(&[1.0, 0.0], 0.0, 1e6, 1.0, 2.0), contact_budget()),
        (support(&[0.0, 1.0], 0.0, 1e6, 1.0, 2.0), contact_budget()),
    ];
    let basis = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let mut s = single_component(2, budget());
    let initial = states(&s);
    for load in [-1.0, -f64::MIN_POSITIVE, -f64::from_bits(1)] {
        assert!(
            s.initialize_contact_equilibrium_with_gauges(
                &[1.0, load],
                &independent,
                &basis,
                set_budget(),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(states(&s), initial);
    }
    for nullspace in [
        vec![],
        vec![vec![1.0, 0.0]],
        vec![vec![1.0, 0.0], vec![1.0, 0.0]],
        vec![vec![1.0, 0.0], vec![0.0, 2.0]],
        vec![vec![1.0, 0.0], vec![0.0, f64::NAN]],
    ] {
        assert!(
            s.initialize_contact_equilibrium_with_gauges(
                &[1.0, 1.0],
                &independent,
                &nullspace,
                set_budget(),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(states(&s), initial);
    }
    for column in [[1.0, 1.0], [1.0, 1.0 + 1e-14], [0.0, 0.0]] {
        let dependent = vec![
            (support(&[1.0, 1.0], 0.0, 1e6, 1.0, 2.0), contact_budget()),
            (support(&column, 0.0, 1e6, 1.0, 2.0), contact_budget()),
        ];
        assert!(
            s.initialize_contact_equilibrium_with_gauges(
                &[1.0, 1.0],
                &dependent,
                &basis,
                set_budget(),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(states(&s), initial);
    }
    let mut redundant = independent.clone();
    redundant.push(independent[0].clone());
    assert!(
        s.initialize_contact_equilibrium_with_gauges(
            &[1.0, 1.0],
            &redundant,
            &basis,
            set_budget(),
            &CancelGate::new()
        )
        .is_err()
    );
    s.initialize_contact_equilibrium_with_gauges(
        &[1.0, 1.0],
        &independent,
        &basis,
        set_budget(),
        &CancelGate::new(),
    )
    .unwrap();
}

#[test]
fn g0_original_force_penetration_energy_work_and_cancellation_caps_remain_transactional() {
    let basis = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let contacts = vec![
        (
            support(&[1.0, 0.0], 0.0003, 1e6, 1.0, 2.0),
            contact_budget(),
        ),
        (
            support(&[0.0, 1.0], 0.0002, 1e6, 1.0, 1.5),
            contact_budget(),
        ),
    ];
    let mut s = single_component(2, budget());
    let initial = states(&s);
    for config in [
        MultiContactConfig {
            max_contacts: 1,
            ..set_budget()
        },
        MultiContactConfig {
            max_setup_terms: 1,
            ..set_budget()
        },
        MultiContactConfig {
            max_sweeps: 0,
            ..set_budget()
        },
    ] {
        assert!(
            s.initialize_contact_equilibrium_with_gauges(
                &[1.0, 1.0],
                &contacts,
                &basis,
                config,
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(states(&s), initial);
    }
    for limits in [
        ModalContactConfig {
            maximum_force_n: 0.5,
            ..contact_budget()
        },
        ModalContactConfig {
            maximum_penetration_m: 0.0005,
            ..contact_budget()
        },
        ModalContactConfig {
            max_iterations: 1,
            ..contact_budget()
        },
    ] {
        let mut bounded = contacts.clone();
        bounded[0].1 = limits;
        assert!(
            s.initialize_contact_equilibrium_with_gauges(
                &[1.0, 1.0],
                &bounded,
                &basis,
                set_budget(),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(states(&s), initial);
    }
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        s.initialize_contact_equilibrium_with_gauges(
            &[1.0, 1.0],
            &contacts,
            &basis,
            set_budget(),
            &gate
        ),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(states(&s), initial);
    s.initialize_contact_equilibrium_with_gauges(
        &[1.0, 1.0],
        &contacts,
        &basis,
        set_budget(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut short = single_component(
        2,
        ModalCouplingConfig {
            maximum_total_energy_j: 1e-6,
            ..budget()
        },
    );
    assert!(
        short
            .initialize_contact_equilibrium_with_gauges(
                &[1.0, 1.0],
                &contacts,
                &basis,
                set_budget(),
                &CancelGate::new()
            )
            .is_err()
    );
    assert!(
        states(&short)
            .iter()
            .all(|q| q.displacement_m_sqrt_kg == 0.0)
    );
    short
        .initialize_contact_equilibrium_with_gauges(
            &[0.001, 0.001],
            &contacts,
            &basis,
            set_budget(),
            &CancelGate::new(),
        )
        .unwrap();
}

#[test]
fn g3_support_equilibration_preserves_distinct_coordinate_and_contact_scales() {
    let columns = [[1e-4, 2e-2], [2e3, -1e6]];
    let expected_reactions = [0.3_f64, 0.0001];
    let load = [
        columns[0][0] * expected_reactions[0] + columns[1][0] * expected_reactions[1],
        columns[0][1] * expected_reactions[0] + columns[1][1] * expected_reactions[1],
    ];
    let contacts = vec![
        (support(&columns[0], 0.0, 1e6, 1.0, 1.5), contact_budget()),
        (support(&columns[1], 0.0, 2e6, 0.7, 2.0), contact_budget()),
    ];
    let mut s = single_component(
        2,
        ModalCouplingConfig {
            maximum_total_energy_j: 1e3,
            ..budget()
        },
    );
    s.initialize_contact_equilibrium_with_gauges(
        &load,
        &contacts,
        &[vec![1.0, 0.0], vec![0.0, 1.0]],
        set_budget(),
        &CancelGate::new(),
    )
    .unwrap();
    let q = states(&s);
    let closures =
        columns.map(|b| b[0] * q[0].displacement_m_sqrt_kg + b[1] * q[1].displacement_m_sqrt_kg);
    let reactions = [1e6 * closures[0].powf(1.5), 2e6 * 0.7 * closures[1].powi(2)];
    for (actual, wanted) in reactions.into_iter().zip(expected_reactions) {
        close(actual, wanted, 2e-9);
    }
    for i in 0..2 {
        close(
            load[i] - columns[0][i] * reactions[0] - columns[1][i] * reactions[1],
            0.0,
            2e-4,
        );
    }
}
