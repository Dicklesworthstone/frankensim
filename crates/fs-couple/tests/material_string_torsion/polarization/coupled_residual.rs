//! Full two-support condensation on the original bowed/obstacle owner.
use super::*;

fn coupled_bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    split: bool,
    rotated: bool,
    rest: f64,
    max_connections: usize,
) -> Result<ReciprocalStringBank, ReciprocalStringError> {
    let (first, second) = if rotated {
        ([0.6, 0.8], [-0.28, 0.96])
    } else {
        ([1.0, 0.0], [0.6, 0.8])
    };
    let base = two_axis_bank(cx, gate, 7800.0, directions(second, first));
    let body = base.system().components()[1].clone();
    let mut strings = base.source().strings().to_vec();
    strings[0].connections[0].rest_extension_m = rest;
    let mut other = second_axis(TransverseDirection::try_new(second).unwrap());
    other.station_fraction = 0.67;
    other.rest_extension_m = -0.7 * rest;
    let bodies = if split {
        other.body.component = 1;
        vec![body.clone(), body]
    } else {
        strings[0].connections[0].body.shapes.resize(4, 0.0);
        other.body.shapes = vec![0.0, 0.0, 0.35, 0.8];
        vec![
            ModalAcousticTimeModel::try_new(
                48_000,
                body.modes().iter().chain(body.modes()).copied().collect(),
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        ]
    };
    let mut config = options();
    config.network.max_connections = max_connections;
    ReciprocalStringBank::new(cx, strings, gas(), bodies, law(), config, gate)?
        .with_translational_body_connections(vec![other], gate)?
        .with_residual_flexibility(gate)
}

// Independent 2x2 inverse from a converged omitted sine sum, not the production
// sequential factors. Signed cross entries follow the physical support axes.
fn stiffness(card: &fs_couple::bowed_string::BowedStringCard) -> [f64; 3] {
    let a = omitted(card, 0.8, 0.8);
    let b = 0.6 * omitted(card, 0.8, 0.67);
    let d = omitted(card, 0.67, 0.67);
    let (k, l) = (3000.0, 1200.0);
    let determinant = 1.0 + k * a + l * d + k * l * (a * d - b * b);
    [
        k * (1.0 + l * d) / determinant,
        -k * l * b / determinant,
        l * (1.0 + k * a) / determinant,
    ]
}

#[test]
fn g1_nonorthogonal_bows_and_stops_use_the_complete_reciprocal_schur_map() {
    with_cx(|cx, gate| {
        let bank = coupled_bank(cx, gate, true, false, 0.0, 8)
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
            .map(|i| bank.system().friction_port_column(i).unwrap().to_vec())
            .collect();
        let bank = bank.with_bow_residual_flexibility(gate).unwrap();
        let card = &bank.cards()[0];
        let k = stiffness(card);
        let stations = [0.2, 0.37];
        let toward = stations.map(|s| [0.6 * omitted(card, s, 0.8), omitted(card, s, 0.67)]);
        let matrix = bank.system().friction_residual_compliance().unwrap();
        let norm = (0.5 * card.length_m * card.linear_density_kg_m).sqrt();
        for i in 0..2 {
            let g = toward[i];
            let factor = [k[0] * g[0] + k[1] * g[1], k[1] * g[0] + k[2] * g[1]];
            let mut expected = original[i].clone();
            for n in 0..4 {
                let order = (n + 1) as f64 * core::f64::consts::PI;
                let (a, b) = ((order * 0.8).sin() / norm, (order * 0.67).sin() / norm);
                expected[n] -= factor[0] * a + factor[1] * 0.6 * b;
                expected[4 + n] -= factor[1] * 0.8 * b;
            }
            for (offset, shapes, coefficient) in
                [(11, [1.0, 0.25], factor[0]), (13, [0.35, 0.8], factor[1])]
            {
                for j in 0..2 {
                    expected[offset + j] += coefficient * shapes[j];
                }
            }
            for (a, b) in bank
                .system()
                .friction_port_column(i)
                .unwrap()
                .iter()
                .zip(expected)
            {
                near(*a, b, 1e-9);
            }
            for j in 0..2 {
                near(
                    matrix[2 * i + j],
                    omitted(card, stations[i], stations[j])
                        - factor[0] * toward[j][0]
                        - factor[1] * toward[j][1],
                    1e-13,
                );
            }
        }
        assert_eq!(matrix[1].to_bits(), matrix[2].to_bits());
        let mut normal = stop(TransverseDirection::try_new([0.6, 0.8]).unwrap());
        normal.location.station_fraction = 0.2;
        let stopped = multi(coupled_bank(cx, gate, true, false, 0.0, 8).unwrap(), gate)
            .with_string_stops(vec![normal], stop_limits(), gate)
            .unwrap();
        let bow = bank.system().friction_port_column(0).unwrap();
        let stop = stopped.system().obstacle_normal_column(0).unwrap();
        assert_eq!(&bow[..8], &stop[..8]);
        assert_eq!(&bow[11..], &stop[11..]);
        assert!(stop[8..11].iter().all(|x| *x == 0.0));
        assert_eq!(
            matrix[0],
            stopped.system().obstacle_residual_compliance().unwrap()[0]
        );
    });
}

#[test]
fn g3_coupled_bridges_preserve_body_partition_rotation_and_contact_work() {
    with_cx(|cx, gate| {
        let make = |split, rotated| {
            multi(
                coupled_bank(cx, gate, split, rotated, 0.0, 8).unwrap(),
                gate,
            )
            .with_bow_residual_flexibility(gate)
            .unwrap()
        };
        let (mut combined, mut split, mut rotated) =
            (make(false, false), make(true, false), make(true, true));
        assert_eq!(
            combined.system().friction_port_column(0),
            split.system().friction_port_column(0)
        );
        let zero = vec![0.0; split.system().mode_count()];
        let mut pressure = 0.0_f64;
        for _ in 0..128 {
            let a = combined.step(&zero, &[drive()]).unwrap().clone();
            let b = split.step(&zero, &[drive()]).unwrap().clone();
            let c = rotated.step(&zero, &[drive()]).unwrap().clone();
            for frame in [&b, &c] {
                near(a.observer_pressure_pa, frame.observer_pressure_pa, 1e-12);
                near(a.drive_work_j, frame.drive_work_j, 1e-17);
                near(a.network_dissipation_j, frame.network_dissipation_j, 1e-17);
                near(
                    a.friction_residual_energy_j,
                    frame.friction_residual_energy_j,
                    1e-17,
                );
            }
            for frame in [&a, &b, &c] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                assert!(frame.network_dissipation_j >= 0.0);
            }
            pressure = pressure.max(b.observer_pressure_pa.abs());
            let q = split.system().components()[0].states();
            let r = rotated.system().components()[0].states();
            for i in 0..4 {
                near(
                    r[i].displacement_m_sqrt_kg,
                    0.6 * q[i].displacement_m_sqrt_kg - 0.8 * q[4 + i].displacement_m_sqrt_kg,
                    1e-14,
                );
                near(
                    r[4 + i].displacement_m_sqrt_kg,
                    0.8 * q[i].displacement_m_sqrt_kg + 0.6 * q[4 + i].displacement_m_sqrt_kg,
                    1e-14,
                );
            }
        }
        assert!(pressure > 1e-8);
        for body in &split.system().components()[1..] {
            assert!(
                body.states()
                    .iter()
                    .any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1e-12)
            );
        }
    });
}

#[test]
fn g0_affine_network_energy_survives_but_unrepresented_contact_offsets_and_budgets_refuse() {
    with_cx(|cx, gate| {
        let rest = 2e-6;
        let bank = coupled_bank(cx, gate, true, false, rest, 8).unwrap();
        let k = stiffness(&bank.cards()[0]);
        let other = -0.7 * rest;
        let expected =
            0.5 * (k[0] * rest * rest + 2.0 * k[1] * rest * other + k[2] * other * other);
        near(bank.system().total_energy_j().unwrap(), expected, 1e-17);
        let error = multi(bank, gate)
            .with_bow_residual_flexibility(gate)
            .err()
            .unwrap();
        assert!(error.to_string().contains("zero bridge rest extensions"));
        assert!(
            multi(coupled_bank(cx, gate, true, false, rest, 8).unwrap(), gate)
                .with_string_stops(vec![stop(TransverseDirection::Y)], stop_limits(), gate)
                .is_err()
        );
        assert!(coupled_bank(cx, gate, true, false, 0.0, 2).is_err());
    });
}
