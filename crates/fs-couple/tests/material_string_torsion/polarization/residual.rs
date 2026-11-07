//! The omitted bending tail follows both physical bridge axes during contact.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{
    MaterialStringStop, StringStopSide,
};
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;

#[path = "coupled_residual.rs"]
mod coupled;

#[path = "rotational_residual.rs"]
mod rotational;

fn multi(bank: ReciprocalStringBank, gate: &CancelGate) -> MultiBowStringBank {
    bank.with_additional_bows(
        vec![],
        MultiFrictionConfig {
            max_ports: 4,
            max_sweeps: 128,
            max_setup_terms: 65_536,
        },
        gate,
    )
    .unwrap()
}

fn stop(direction: TransverseDirection) -> MaterialStringStop {
    MaterialStringStop {
        location: StringBowLocation {
            string: 0,
            station_fraction: 0.43,
        },
        side: StringStopSide::Along(direction),
        body: None,
        law: Obstacle::new(
            vec![-1.0],
            1,
            1,
            vec![0.0],
            vec![1.0],
            5000.0,
            1.0,
            "synthetic two-axis residual contact".into(),
        )
        .unwrap()
        .with_internal_loss(0.2)
        .unwrap(),
        limits: ModalContactConfig {
            max_iterations: 96,
            maximum_force_n: 10.0,
            maximum_penetration_m: 0.001,
            force_absolute_tolerance_n: 1e-11,
            force_relative_tolerance: 1e-9,
        },
    }
}

fn stop_limits() -> FrictionObstacleConfig {
    FrictionObstacleConfig {
        max_contacts: 4,
        max_setup_terms: 65_536,
    }
}

// Independent convergent sine series. EI > 0 makes this omitted tail O(n^-4).
fn omitted(card: &fs_couple::bowed_string::BowedStringCard, x: f64, y: f64) -> f64 {
    (card.mode_count + 1..=20_000)
        .map(|n| {
            let a = n as f64 * core::f64::consts::PI;
            let k = a / card.length_m;
            2.0 * (a * x).sin() * (a * y).sin()
                / (card.length_m
                    * (card.tension_n * k * k + card.bending_stiffness_n_m2 * k.powi(4)))
        })
        .sum()
}

#[test]
fn g4_bowed_residual_second_axis_matches_its_scalar_mirror_and_preserves_work() {
    with_cx(|cx, gate| {
        let mirror = StringBodyTranslationConnection {
            string: 0,
            station_fraction: 0.8,
            direction: TransverseDirection::Y,
            body: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 0.25],
            },
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        };
        let make = |bow| {
            let bank = two_axis_bank(cx, gate, 7800.0, directions(bow, [1.0, 0.0]))
                .with_residual_flexibility(gate)
                .unwrap()
                .with_translational_body_connections(vec![mirror.clone()], gate)
                .unwrap();
            let original_torsion = bank.system().port().left.shapes[8..].to_vec();
            let bank = multi(bank, gate)
                .with_bow_residual_flexibility(gate)
                .unwrap();
            assert_eq!(
                bank.system().port(0).unwrap().left.shapes[8..],
                original_torsion
            );
            bank
        };
        let (mut x, mut y) = (make([1.0, 0.0]), make([0.0, 1.0]));
        assert_eq!(
            x.system().friction_residual_compliance(),
            y.system().friction_residual_compliance()
        );
        let zero = vec![0.0; x.system().mode_count()];
        let (mut peak, mut tail_energy) = (0.0_f64, 0.0_f64);
        for _ in 0..128 {
            let a = x.step(&zero, &[drive()]).unwrap().clone();
            let b = y.step(&zero, &[drive()]).unwrap().clone();
            near(a.observer_pressure_pa, b.observer_pressure_pa, 1e-12);
            near(a.drive_work_j, b.drive_work_j, 1e-17);
            near(
                a.friction_residual_energy_j,
                b.friction_residual_energy_j,
                1e-17,
            );
            near(
                a.friction_residual_loss_j,
                b.friction_residual_loss_j,
                1e-17,
            );
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            peak = peak.max(b.observer_pressure_pa.abs());
            tail_energy = tail_energy.max(b.friction_residual_energy_j);
        }
        assert!(
            peak > 1e-8 && tail_energy > 0.0,
            "second-axis bow did not load the body and omitted tail"
        );
    });
}

#[test]
fn g1_oblique_bows_condense_both_axes_in_the_joint_omitted_mode_matrix() {
    with_cx(|cx, gate| {
        let bank = two_axis_bank(cx, gate, 7800.0, directions([0.6, 0.8], [1.0, 0.0]))
            .with_translational_body_connections(vec![second_axis(TransverseDirection::Y)], gate)
            .unwrap()
            .with_residual_flexibility(gate)
            .unwrap()
            .with_additional_bows(
                vec![MaterialBowContact {
                    location: StringBowLocation {
                        string: 0,
                        station_fraction: 0.37,
                    },
                    law: law(),
                    friction: options().friction,
                }],
                MultiFrictionConfig {
                    max_ports: 4,
                    max_sweeps: 128,
                    max_setup_terms: 65_536,
                },
                gate,
            )
            .unwrap();
        let original: Vec<_> = (0..2)
            .map(|i| bank.system().port(i).unwrap().left.shapes.clone())
            .collect();
        let bank = bank.with_bow_residual_flexibility(gate).unwrap();
        let card = &bank.cards()[0];
        let gbb = omitted(card, 0.8, 0.8);
        let (kx, ky) = (3000.0 / (1.0 + 3000.0 * gbb), 1200.0 / (1.0 + 1200.0 * gbb));
        let stations = [0.2, 0.37];
        let cross = stations.map(|s| omitted(card, s, 0.8));
        let matrix = bank.system().friction_residual_compliance().unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let expected = omitted(card, stations[i], stations[j])
                    - (0.36 * kx + 0.64 * ky) * cross[i] * cross[j];
                assert!(
                    (matrix[2 * i + j] - expected).abs() < 1e-13,
                    "joint residual {i},{j}"
                );
            }
            let shapes = &bank.system().port(i).unwrap().left.shapes;
            let norm = (0.5 * card.length_m * card.linear_density_kg_m).sqrt();
            for n in 0..4 {
                let bridge = ((n + 1) as f64 * core::f64::consts::PI * 0.8).sin() / norm;
                near(
                    shapes[n],
                    original[i][n] - 0.6 * kx * cross[i] * bridge,
                    1e-9,
                );
                near(
                    shapes[4 + n],
                    original[i][4 + n] - 0.8 * ky * cross[i] * bridge,
                    1e-9,
                );
            }
            assert_eq!(shapes[8..], original[i][8..]);
        }
        assert_eq!(matrix[1].to_bits(), matrix[2].to_bits());
    });
}

#[test]
fn g3_oblique_stops_keep_two_axis_residuals_under_a_common_basis_rotation() {
    with_cx(|cx, gate| {
        let make = |axes: TransverseDirections, second| {
            let bank = two_axis_bank(cx, gate, 7800.0, axes)
                .with_translational_body_connections(vec![second_axis(second)], gate)
                .unwrap()
                .with_residual_flexibility(gate)
                .unwrap();
            multi(bank, gate)
                .with_string_stops(
                    vec![stop(TransverseDirection::try_new(axes.bow()).unwrap())],
                    stop_limits(),
                    gate,
                )
                .unwrap()
        };
        let mut a = make(directions([0.6, 0.8], [1.0, 0.0]), TransverseDirection::Y);
        let mut b = make(
            directions([-0.28, 0.96], [0.6, 0.8]),
            TransverseDirection::try_new([-0.8, 0.6]).unwrap(),
        );
        let tail = a.system().obstacle_residual_compliance().unwrap()[0];
        near(
            tail,
            b.system().obstacle_residual_compliance().unwrap()[0],
            1e-17,
        );
        let card = &a.cards()[0];
        let gbb = omitted(card, 0.8, 0.8);
        let cross = omitted(card, 0.43, 0.8);
        let expected = omitted(card, 0.43, 0.43)
            - (0.36 * 3000.0 / (1.0 + 3000.0 * gbb) + 0.64 * 1200.0 / (1.0 + 1200.0 * gbb))
                * cross
                * cross;
        assert!((tail - expected).abs() < 1e-13);
        let force = |bank: &MultiBowStringBank| {
            let port = StringForcePort {
                location: StringBowLocation {
                    string: 0,
                    station_fraction: 0.43,
                },
                coupling: StringForceCoupling::Centreline,
            };
            bank.string_force_columns(&[port], 4096, gate).unwrap()[0]
                .iter()
                .map(|v| 0.001 * v)
                .collect::<Vec<_>>()
        };
        let (fa, fb) = (force(&a), force(&b));
        let zero = vec![0.0; fa.len()];
        let released = [FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        }];
        let mut engaged = false;
        for n in 0..96 {
            let x = a
                .step(if n < 32 { &fa } else { &zero }, &released)
                .unwrap()
                .clone();
            let y = b
                .step(if n < 32 { &fb } else { &zero }, &released)
                .unwrap()
                .clone();
            near(x.observer_pressure_pa, y.observer_pressure_pa, 1e-12);
            near(
                x.obstacle_residual_energy_j,
                y.obstacle_residual_energy_j,
                1e-17,
            );
            assert!(x.energy_residual_j.abs() <= x.energy_tolerance_j);
            assert!(y.energy_residual_j.abs() <= y.energy_tolerance_j);
            engaged |= a.system().last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
        }
        assert!(engaged, "oblique two-axis residual stop never engaged");
        assert!(
            a.system().components()[0].states()[8..]
                .iter()
                .all(|s| s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0)
        );
    });
}

// Identical physical modes/ports, stored as one four-mode component or two
// independent two-mode bodies. Splitting storage cannot change mechanics.
fn partitioned_bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    split: bool,
    bow: [f64; 2],
) -> MultiBowStringBank {
    let base = two_axis_bank(cx, gate, 7800.0, directions(bow, [1.0, 0.0]));
    let body = base.system().components()[1].clone();
    let mut strings = base.source().strings().to_vec();
    let mut other = second_axis(TransverseDirection::Y);
    let bodies = if split {
        other.body.component = 1;
        vec![body.clone(), body]
    } else {
        strings[0].connections[0].body.shapes.resize(4, 0.0);
        let second = other.body.shapes.clone();
        other.body.shapes = vec![0.0; 2];
        other.body.shapes.extend(second);
        vec![
            ModalAcousticTimeModel::try_new(
                48_000,
                body.modes().iter().chain(body.modes()).copied().collect(),
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        ]
    };
    let bank = ReciprocalStringBank::new(cx, strings, gas(), bodies, law(), options(), gate)
        .unwrap()
        .with_translational_body_connections(vec![other], gate)
        .unwrap()
        .with_residual_flexibility(gate)
        .unwrap();
    multi(bank, gate)
}

#[test]
fn g3_oblique_residual_bow_keeps_independent_bridge_bodies_and_their_work() {
    with_cx(|cx, gate| {
        let make = |split, bow| {
            partitioned_bank(cx, gate, split, bow)
                .with_bow_residual_flexibility(gate)
                .unwrap()
        };
        // Zero-axis corrections remain admitted and add no artificial receiver.
        let _ = make(true, [1.0, 0.0]);
        let _ = make(true, [0.0, 1.0]);
        let (mut combined, mut split) = (make(false, [0.6, 0.8]), make(true, [0.6, 0.8]));
        assert_eq!(
            combined.system().friction_residual_compliance(),
            split.system().friction_residual_compliance()
        );
        let zero = vec![0.0; split.system().mode_count()];
        for _ in 0..128 {
            let a = combined.step(&zero, &[drive()]).unwrap().clone();
            let b = split.step(&zero, &[drive()]).unwrap().clone();
            near(a.observer_pressure_pa, b.observer_pressure_pa, 1e-12);
            near(a.drive_work_j, b.drive_work_j, 1e-17);
            near(
                a.friction_residual_energy_j,
                b.friction_residual_energy_j,
                1e-17,
            );
            near(
                a.friction_residual_loss_j,
                b.friction_residual_loss_j,
                1e-17,
            );
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        }
        for body in &split.system().components()[1..] {
            assert!(
                body.states()
                    .iter()
                    .any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-12),
                "an independent bridge body received no reciprocal reaction"
            );
        }
        for (a, b) in combined.system().components()[1].states().iter().zip(
            split.system().components()[1..]
                .iter()
                .flat_map(|body| body.states()),
        ) {
            near(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 1e-14);
            near(
                a.velocity_m_sqrt_kg_per_s,
                b.velocity_m_sqrt_kg_per_s,
                1e-12,
            );
        }
    });
}

#[test]
fn g3_oblique_residual_stop_keeps_independent_bridge_bodies_and_their_work() {
    with_cx(|cx, gate| {
        let normal = TransverseDirection::try_new([0.6, 0.8]).unwrap();
        let make = |split| {
            partitioned_bank(cx, gate, split, normal.components())
                .with_string_stops(vec![stop(normal)], stop_limits(), gate)
                .unwrap()
        };
        let (mut combined, mut split) = (make(false), make(true));
        let column = split.system().obstacle_normal_column(0).unwrap();
        assert!(
            column[11..13].iter().any(|x| *x != 0.0) && column[13..15].iter().any(|x| *x != 0.0)
        );
        assert_eq!(column, combined.system().obstacle_normal_column(0).unwrap());
        let port = StringForcePort {
            location: StringBowLocation {
                string: 0,
                station_fraction: 0.43,
            },
            coupling: StringForceCoupling::Centreline,
        };
        let force: Vec<_> = split.string_force_columns(&[port], 4096, gate).unwrap()[0]
            .iter()
            .map(|x| 0.001 * x)
            .collect();
        let zero = vec![0.0; force.len()];
        let released = [FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        }];
        let mut engaged = false;
        for n in 0..96 {
            let external = if n < 32 { &force } else { &zero };
            let a = combined.step(external, &released).unwrap().clone();
            let b = split.step(external, &released).unwrap().clone();
            near(a.observer_pressure_pa, b.observer_pressure_pa, 1e-12);
            near(a.external_work_j, b.external_work_j, 1e-17);
            near(
                a.obstacle_residual_energy_j,
                b.obstacle_residual_energy_j,
                1e-17,
            );
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            engaged |= split.system().last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
        }
        assert!(
            engaged,
            "independently supported oblique residual stop never engaged"
        );
    });
}
