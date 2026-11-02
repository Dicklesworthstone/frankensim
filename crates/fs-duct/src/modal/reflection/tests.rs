use super::*;
use fs_material::gas::GasSpec;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn multiply(matrix: &[C64], vector: &[C64]) -> Vec<C64> {
    let n = vector.len();
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |sum, j| sum + matrix[i * n + j] * vector[j]))
        .collect()
}

fn mixed_nodes() -> ModalReflectionResponse {
    // An exact lossless pressure/flow transformer: P0 = 2 P1, U1 = -2 U0.
    // Native impedances 4 and 1 make R = [[0, 2], [1/2, 0]], whose eigenvalues
    // are exactly +1 and -1. Neither full impedance nor full admittance exists.
    ModalReflectionResponse {
        omega: 1000.0,
        n_modes: 2,
        reflection_matrix: vec![C64::ZERO, C64::from_re(2.0), C64::from_re(0.5), C64::ZERO],
        characteristic_impedances: vec![C64::from_re(4.0), C64::ONE],
        min_shear_number: f64::INFINITY,
        mouth_ka: 0.1,
        cutoffs_at_throat_hz: vec![100.0],
        steps_used: 1,
    }
}

#[test]
fn exact_mixed_nodes_convert_without_impedance_and_preserve_physical_work() {
    let response = mixed_nodes();
    assert!(matches!(
        impedance_from_reflection(&response.reflection_matrix, &response.characteristic_impedances, 2),
        Err(DuctError::Singular)
    ));
    for references in [[4.0, 1.0], [3.0, 7.0], [1.0e-4, 1.0e4]] {
        let reflection = response.at_references(&references).unwrap();
        let outgoing = [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)];
        let incoming = multiply(&reflection, &outgoing);
        let pressure: Vec<_> = outgoing.iter().zip(&incoming).map(|(a, b)| *a + *b).collect();
        let flow: Vec<_> = outgoing.iter().zip(&incoming).zip(references)
            .map(|((a, b), z)| (*a - *b).scale(1.0 / z)).collect();
        assert!((pressure[0] - pressure[1].scale(2.0)).abs() < 2.0e-12);
        assert!((flow[1] + flow[0].scale(2.0)).abs() < 2.0e-11);
        let incident_power: f64 = outgoing.iter().zip(references).map(|(a, z)| a.abs().powi(2) / z).sum();
        let reflected_power: f64 = incoming.iter().zip(references).map(|(b, z)| b.abs().powi(2) / z).sum();
        assert!((incident_power - reflected_power).abs() < 2.0e-12 * incident_power);
        assert!((reflection[1].scale(1.0 / references[0])
            - reflection[2].scale(1.0 / references[1])).abs() < 2.0e-12);
    }
}

#[test]
fn plane_source_preserves_zero_higher_flow_instead_of_matching_it() {
    let response = mixed_nodes();
    for reference in [0.5, 3.0, 20.0] {
        let plane = response.plane_reflection(reference).unwrap();
        // U1 = 0 forces U0 = 0 for this physical transformer, so a plane
        // source sees a rigid boundary. Matching port 1 permits power transfer.
        assert!((plane - C64::ONE).abs() < 2.0e-14);
        let matched = response.at_references(&[reference; 2]).unwrap();
        assert!((matched[0] - plane).abs() > 0.1);
    }
}

#[test]
fn prepared_reflection_matches_impedance_away_from_poles_for_every_load() {
    let gas = air();
    let exterior = GasState::try_new_moist_air(303.15, gas.pressure, 0.7).unwrap();
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.006, length: 0.025 },
        Segment::Cylinder { radius: 0.012, length: 0.035 },
        Segment::Cylinder { radius: 0.008, length: 0.020 },
    ] };
    let table = crate::TabulatedLoad::try_new(
        vec![500.0, 20_000.0],
        vec![C64::new(2.0e6, -3.0e5), C64::new(3.0e6, -5.0e5)],
        "manufactured finite passive mouth for reflection equivalence", 0.0,
    ).unwrap();
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.006_f64.powi(2));
    for n in [1, 3] {
        let geometry = PreparedModalDuct::new(&duct, n, 1).unwrap();
        let references: Vec<_> = (0..n).map(|i| z0 * (1.0 + 0.7 * i as f64)).collect();
        let reference_complex: Vec<_> = references.iter().map(|z| C64::from_re(*z)).collect();
        for omega in [1000.0, 8000.0, 16_000.0] {
            for loss in [LossModel::Lossless, LossModel::Bessel] {
                let mut cases = Vec::new();
                for termination in [Termination::Closed, Termination::IdealOpen, Termination::FlangedOpen] {
                    cases.push((
                        geometry.input_reflection(&gas, omega, loss, termination).unwrap(),
                        geometry.input_impedance(&gas, omega, loss, termination).unwrap(),
                    ));
                }
                cases.push((
                    geometry.input_reflection_tabulated(&gas, omega, loss, &table).unwrap(),
                    geometry.input_impedance_tabulated(&gas, omega, loss, &table).unwrap(),
                ));
                cases.push((
                    geometry.input_reflection_load(&gas, omega, loss, C64::from_re(2.0e6)).unwrap(),
                    geometry.input_impedance_load(&gas, omega, loss, C64::from_re(2.0e6)).unwrap(),
                ));
                cases.push((
                    geometry.input_reflection_flanged(&gas, &exterior, omega, loss).unwrap(),
                    geometry.input_impedance_flanged(&gas, &exterior, omega, loss).unwrap(),
                ));
                for (waves, impedance) in cases {
                    let expected = modal_reflection_from_impedance(&impedance.impedance_matrix, &reference_complex, n).unwrap();
                    let actual = waves.at_references(&references).unwrap();
                    for (a, b) in actual.iter().zip(expected) {
                        assert!((*a - b).abs() < 5.0e-10, "modes={n} omega={omega}: {a:?} vs {b:?}");
                    }
                    let z = impedance.plane_impedance;
                    let expected_plane = (z - C64::from_re(z0)) / (z + C64::from_re(z0));
                    assert!((waves.plane_reflection(z0).unwrap() - expected_plane).abs() < 5.0e-10);
                    assert_eq!(waves.steps_used, impedance.steps_used);
                    assert_eq!(waves.min_shear_number, impedance.min_shear_number);
                    assert_eq!(waves.mouth_ka, impedance.mouth_ka);
                    assert_eq!(waves.cutoffs_at_throat_hz, impedance.cutoffs_at_throat_hz);
                }
            }
        }
    }
}

#[test]
fn physical_closed_and_open_cylinders_keep_reflection_at_input_nodes() {
    let gas = air();
    let (radius, length) = (0.006, 0.18);
    let duct = Duct { segments: vec![Segment::Cylinder { radius, length }] };
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    for n in [1, 4] {
        let geometry = PreparedModalDuct::new(&duct, n, 1).unwrap();
        for termination in [Termination::Closed, Termination::IdealOpen] {
            let sign = if matches!(termination, Termination::Closed) { 1.0 } else { -1.0 };
            for phase in [core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
                for detuning in [-1.0e-12, 0.0, 1.0e-12] {
                    let phase = phase * (1.0 + detuning);
                    let omega = phase * gas.sound_speed / length;
                    let waves = geometry.input_reflection(&gas, omega, LossModel::Lossless, termination).unwrap();
                    // Uniform cylinder: plane reflection is the terminal sign
                    // times its physical round-trip phase, including both nodes.
                    let expected = C64::new((2.0 * phase).cos(), (2.0 * phase).sin()).scale(sign);
                    let plane = waves.plane_reflection(reference).unwrap();
                    let matrix = waves.at_references(&vec![reference; n]).unwrap();
                    assert!((plane - expected).abs() < 5.0e-13);
                    assert!((matrix[0] - expected).abs() < 5.0e-13);
                    assert!((plane.abs() - 1.0).abs() < 5.0e-14);
                }
            }
        }
    }
}

#[test]
fn reflection_refuses_invalid_references_and_keeps_load_admission() {
    let mut response = mixed_nodes();
    for references in [vec![], vec![1.0], vec![1.0, 0.0], vec![1.0, -1.0], vec![1.0, f64::NAN]] {
        assert!(response.at_references(&references).is_err());
    }
    for reference in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(response.plane_reflection(reference).is_err());
    }
    response.n_modes = usize::MAX;
    assert!(response.at_references(&[1.0, 1.0]).is_err());
    assert!(response.plane_reflection(1.0).is_err());
    let gas = air();
    let duct = Duct { segments: vec![Segment::Cylinder { radius: 0.006, length: 0.18 }] };
    let geometry = PreparedModalDuct::new(&duct, 2, 1).unwrap();
    assert!(geometry.input_reflection(&gas, 0.0, LossModel::Lossless, Termination::Closed).is_err());
    assert!(geometry.input_reflection_load(&gas, 1000.0, LossModel::Lossless, C64::from_re(-1.0)).is_err());
    let other_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
    assert!(geometry.input_reflection_flanged(&gas, &other_pressure, 1000.0, LossModel::Lossless).is_err());
}
