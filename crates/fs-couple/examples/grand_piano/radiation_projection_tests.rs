//! G0/G3: acoustic port reduction preserves the complete mechanical power pair.
use super::*;

fn basis(board_ports: usize) -> PortBasis {
    let mut rows = vec![vec![0.; board_ports]; 2];
    rows[0][0] = 0.6;
    rows[0][board_ports - 1] = 0.8;
    rows[1][board_ports / 2] = 1.;
    PortBasis::new(board_ports, rows).unwrap()
}

fn model(zeta: f64) -> Model {
    Model {
        ports: 2,
        poles: vec![
            Pole {
                omega: 1200.,
                zeta,
                coupling: vec![180., -90.],
            },
            Pole {
                omega: 3200.,
                zeta,
                coupling: vec![50., 220.],
            },
        ],
    }
}

fn kinetic(v: &[f64]) -> f64 {
    0.5 * v.iter().map(|x| x * x).sum::<f64>()
}

#[test]
fn projected_air_closes_full_board_energy_and_leaves_the_nullspace_alone() {
    let basis = basis(128);
    assert_eq!(basis.board_ports(), 128);
    assert_eq!(basis.ports(), 2);
    for zeta in [0., 0.2] {
        let mut air = Prepared::new_projected(&model(zeta), 192_000, &basis).unwrap();
        let mut v = vec![0.; 128];
        v[0] = 0.03;
        v[127] = -0.02;
        v[64] = 0.05;
        v[89] = 0.011; // a mechanical coordinate outside both acoustic rows
        let initial = kinetic(&v);
        let null_velocity = 0.8 * v[0] - 0.6 * v[127];
        let mut loss = 0.;
        for _ in 0..2400 {
            loss += air.before(&mut v).unwrap();
            loss += air.after(&mut v).unwrap();
        }
        let defect = kinetic(&v) + air.energy() + loss - initial;
        assert!(defect.abs() < 2e-13, "full-board work defect {defect:e}");
        assert!((0.8 * v[0] - 0.6 * v[127] - null_velocity).abs() < 2e-13);
        assert_eq!(v[89].to_bits(), 0.011_f64.to_bits());
        assert!(air.energy() > 1e-10);
        if zeta > 0. {
            assert!(loss > 1e-8);
        } else {
            assert!(loss.abs() < 2e-13);
        }
    }

    // A velocity in ker(Q^T) cannot create air energy or acoustic drag.
    let mut air = Prepared::new_projected(&model(0.2), 192_000, &basis).unwrap();
    let mut v = vec![0.; 128];
    v[0] = 0.8;
    v[127] = -0.6;
    v[89] = 0.03;
    let original = v.clone();
    for _ in 0..64 {
        assert_eq!(air.before(&mut v).unwrap(), 0.);
        assert_eq!(air.after(&mut v).unwrap(), 0.);
    }
    assert_eq!(v, original);
    assert_eq!(air.energy(), 0.);
}

#[test]
fn lifted_impedance_preserves_complex_force_velocity_work() {
    let basis = basis(40);
    let v: Vec<_> = (0..40)
        .map(|j| C64::new((j % 5) as f64 - 2., j as f64 / 41.))
        .collect();
    let projected: Vec<_> = basis
        .vectors()
        .iter()
        .map(|row| {
            row.iter()
                .zip(&v)
                .fold(C64::ZERO, |sum, (&q, &v)| sum + v.scale(q))
        })
        .collect();
    for omega in [50., 1200., 3200., 10000.] {
        let z = model(0.2).impedance(omega).unwrap();
        let lifted = basis.lift_impedance(&z).unwrap();
        assert_eq!(lifted.len(), 40 * 40);
        let recovered = basis.project_impedance(&lifted).unwrap();
        for (a, b) in z.iter().zip(recovered) {
            assert!((*a - b).abs() < 1e-12 * (1. + a.abs()));
        }
        let work = |velocity: &[C64], impedance: &[C64]| {
            let n = velocity.len();
            (0..n).fold(C64::ZERO, |sum, i| {
                sum + velocity[i].conj()
                    * (0..n).fold(C64::ZERO, |force, j| {
                        force + impedance[i * n + j] * velocity[j]
                    })
            })
        };
        let full = work(&v, &lifted);
        let reduced = work(&projected, &z);
        assert!((full - reduced).abs() < 1e-12 * (1. + reduced.abs()));
        assert!(full.re >= -1e-12);
    }
    assert!(basis.project_impedance(&[C64::ZERO; 4]).is_err());
    assert!(basis.lift_impedance(&[C64::ZERO; 3]).is_err());
}

#[test]
fn invalid_port_spaces_are_refused_and_failed_flows_preserve_air_history() {
    for (count, rows) in [
        (0, vec![]),
        (129, vec![vec![0.; 129]]),
        (4, vec![]),
        (4, vec![vec![1., 0., 0.]]),
        (4, vec![vec![0.; 4]]),
        (4, vec![vec![2., 0., 0., 0.]]),
        (4, vec![vec![f64::NAN, 0., 0., 0.]]),
        (4, vec![vec![1., 0., 0., 0.], vec![1., 0., 0., 0.]]),
        (
            40,
            (0..33)
                .map(|j| (0..40).map(|i| if i == j { 1. } else { 0. }).collect())
                .collect(),
        ),
    ] {
        assert!(PortBasis::new(count, rows).is_err());
    }
    let basis = basis(40);
    let mut bad = model(0.2);
    bad.ports = 1;
    assert!(Prepared::new_projected(&bad, 192_000, &basis).is_err());
    assert!(Prepared::new_projected(&model(0.2), 1000, &basis).is_err());

    let mut air = Prepared::new_projected(&model(0.2), 192_000, &basis).unwrap();
    let mut v = vec![0.01; 40];
    air.before(&mut v).unwrap();
    air.after(&mut v).unwrap();
    air.checkpoint();
    let original = v.clone();
    let states = air.state.clone();
    for after in [false, true] {
        let mut invalid = original.clone();
        invalid[17] = f64::NAN; // outside the acoustic rows, still invalid mechanics
        let bits: Vec<_> = invalid.iter().map(|x| x.to_bits()).collect();
        assert!(
            if after {
                air.after(&mut invalid)
            } else {
                air.before(&mut invalid)
            }
            .is_err()
        );
        assert_eq!(air.state, states);
        assert_eq!(
            invalid.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            bits
        );
        assert!(
            if after {
                air.after(&mut [0.; 39])
            } else {
                air.before(&mut [0.; 39])
            }
            .is_err()
        );
        assert_eq!(air.state, states);
    }
    let loss = air.before(&mut v).unwrap() + air.after(&mut v).unwrap();
    let result = v.clone();
    let candidate = air.state.clone();
    air.restore();
    v.copy_from_slice(&original);
    let replay = air.before(&mut v).unwrap() + air.after(&mut v).unwrap();
    assert_eq!(v, result);
    assert_eq!(air.state, candidate);
    assert_eq!(loss.to_bits(), replay.to_bits());
}

fn piano() -> super::super::Instrument {
    use super::super::super::{geometry, linear::BoardMode};
    let course = geometry::demonstration_scale().unwrap()[48];
    let board: Vec<_> = (0..40)
        .map(|j| {
            let mut bridge = [0.; 88];
            bridge[48] = 0.005 * (j + 1) as f64;
            BoardMode {
                frequency_hz: 90. + 7. * j as f64,
                damping_ratio: 0.015,
                bridge,
                volume: 0.001 / (j + 1) as f64,
            }
        })
        .collect();
    super::super::Instrument::new(vec![course], &board, 48_000, 4, 12, true).unwrap()
}

#[test]
fn projected_engine_admission_is_cold_atomic_and_zero_coupling_preserves_mechanics() {
    let mut loaded = piano();
    let mut bare = piano();
    let basis = basis(40);
    assert!(
        loaded
            .configure_projected_radiation(
                &model(0.2),
                &super::PortBasis::new(
                    39,
                    vec![(0..39).map(|j| if j == 0 { 1. } else { 0. }).collect()]
                )
                .unwrap()
            )
            .is_err()
    );
    assert!(!loaded.has_radiation());
    assert_eq!(loaded.bank.q, bare.bank.q);
    let mut zero = model(0.2);
    for pole in &mut zero.poles {
        pole.coupling.fill(0.);
    }
    loaded.configure_projected_radiation(&zero, &basis).unwrap();
    assert_eq!(loaded.bank.board_count, 40);
    assert!(loaded.configure_projected_radiation(&zero, &basis).is_err());
    loaded.note_on(69, 0.5).unwrap();
    bare.note_on(69, 0.5).unwrap();
    assert!(bare.configure_projected_radiation(&zero, &basis).is_err());
    assert!(!bare.has_radiation());
    for _ in 0..480 {
        loaded.step().unwrap();
        bare.step().unwrap();
    }
    assert_eq!(loaded.bank.q, bare.bank.q);
    assert_eq!(loaded.bank.v, bare.bank.v);
    assert_eq!(loaded.radiation_energy_j(), 0.);
    assert_eq!(loaded.accounting.radiation_loss_j, 0.);
}

#[test]
fn failed_piano_samples_restore_projected_air_and_every_board_coordinate() {
    let mut a = piano();
    let mut b = piano();
    let basis = basis(40);
    a.configure_projected_radiation(&model(0.2), &basis)
        .unwrap();
    b.configure_projected_radiation(&model(0.2), &basis)
        .unwrap();
    a.note_on(69, 0.5).unwrap();
    b.note_on(69, 0.5).unwrap();
    let mut trace = vec![0.; a.board_trace_len()];
    for _ in 0..480 {
        a.step_with_board_trace(&mut trace).unwrap();
        for sub in 0..b.substeps {
            b.mechanics_step().unwrap();
            assert_eq!(
                &trace[sub * 40..(sub + 1) * 40],
                &b.bank.v[b.bank.modes.len()..]
            );
        }
    }
    assert!(a.radiation_energy_j() > 0.);
    let energy = a.energy_j();
    let loss = a.accounting.radiation_loss_j;
    a.damper_drag_ns_m = f64::NAN;
    assert!(a.step_with_board_trace(&mut trace).is_err());
    assert_eq!(a.energy_j(), energy);
    assert_eq!(a.accounting.radiation_loss_j, loss);
    assert_eq!(a.bank.q, b.bank.q);
    assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
    a.damper_drag_ns_m = 0.4;
    a.step().unwrap();
    b.step().unwrap();
    assert_eq!(a.bank.q, b.bank.q);
    assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
}
