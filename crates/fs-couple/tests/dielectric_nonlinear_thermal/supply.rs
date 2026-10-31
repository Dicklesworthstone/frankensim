//! Independent source work, reversible heat and retained thermal authority.
use super::*;

fn field_state(q: f64, temperature: f64, alpha: f64) -> (f64, f64) {
    let scale = 1.0 + alpha * (temperature - 300.0);
    let z = (q / scale).abs();
    (scale * density(z), alpha * (z * voltage(z) - density(z)))
}

#[test]
fn lossless_charge_drive_has_signed_reversible_heat_and_exact_original_source_work() {
    let capacity = 0.3;
    let dt = 0.1;
    for alpha in [-0.001, 0.001] {
        for sign in [-1.0, 1.0] {
            let q0 = sign * 1.9;
            let q1 = sign * 2.1;
            let current = (q1 - q0) / dt;
            let card = fixture_scaled(capacity, false, 0.0, [250.0, 500.0], alpha, 0.0);
            let curve = caloric(&card);
            let thermal = thermal(&[&curve], &[0.0]);
            let context = point(300.0);
            let cells = [cell(0, &card, &curve, &context)];
            let network = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let (seed, poles, law) = initial_scaled(&card, q0, None, false);
            let circuit = DebyeCircuit::new_with_thermal_reversible_instantaneous(
                CircuitGraph {
                    node_count: 2,
                    branches: vec![
                        (1, 0, Branch::Capacitor { farads: f64::NAN }),
                        (1, 0, Branch::CurrentSource { port: 0 }),
                    ],
                    transformers: vec![],
                },
                seed.source().clone(),
                seed.reversible_instantaneous().unwrap().clone(),
                seed.thermal_reversible_storage().unwrap().clone(),
                vec![],
                vec![q0],
                vec![],
                &[current],
                || false,
            )
            .unwrap();
            let initial = field_state(q0, 300.0, alpha);
            // Same native endpoint-temperature work split, solved separately
            // using only hand inverse/integral formulas and a scalar bracket.
            let equation = |t| {
                let bridge = field_state(q0, t, alpha);
                let after = field_state(q1, t, alpha);
                capacity * (t - 300.0) + bridge.0 - initial.0 + t * after.1 - 300.0 * initial.1
            };
            let (mut lo, mut hi) = (290.0, 310.0);
            assert!(equation(lo) < 0.0 && equation(hi) > 0.0);
            for _ in 0..90 {
                let mid = f64::midpoint(lo, hi);
                if equation(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let expected_t = f64::midpoint(lo, hi);
            let bridge = field_state(q0, expected_t, alpha);
            let after = field_state(q1, expected_t, alpha);
            let mut bound = binding(&poles, &[]);
            bound.maximum_instantaneous_charge_c = 3.0;
            let next = dielectric::step_candidate_with_reversible_storage(
                &network,
                &thermal,
                dt,
                &circuit,
                &[current],
                bound,
                &law,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
            let report = next.report();
            let source = &report.coupled.source;
            near(next.circuit().instantaneous_charge_c(), q1, 2e-11);
            near(
                next.thermal().temperature_k(DIELECTRIC).unwrap(),
                expected_t,
                2e-8,
            );
            near(source.parameter_work_j, bridge.0 - initial.0, 2e-10);
            near(
                source.thermal_offset_change_j,
                expected_t * after.1 - 300.0 * initial.1,
                2e-10,
            );
            near(
                source.reversible_heat_j,
                capacity * (expected_t - 300.0),
                2e-10,
            );
            assert!(alpha * source.reversible_heat_j < 0.0);
            assert_eq!(report.coupled.electrical.polarization_heat_j, 0.0);
            assert_eq!(report.coupled.electrical.other_resistor_heat_j, 0.0);
            near(
                report.coupled.electrical.native.supplied,
                after.0 - bridge.0,
                2e-10,
            );
            near(report.total_energy_residual_j, 0.0, 4e-10);
            assert!(q0.abs() < 2.0 && q1.abs() > 2.0 * (1.0 + alpha * (expected_t - 300.0)));
            assert_eq!(circuit.instantaneous_charge_c(), q0);
            assert_eq!(thermal.temperature_k(DIELECTRIC).unwrap(), 300.0);
        }
    }
}

#[test]
fn thermal_authority_survives_native_rebind_and_legacy_or_inconsistent_paths_refuse() {
    let card = fixture_scaled(0.3, true, 0.0, [250.0, 500.0], 0.001, 0.0006);
    let (circuit, poles, law) = initial_scaled(&card, 1.9, Some(0.4), false);
    let saved = circuit.state().to_vec();
    let source = poles.at_temperature(310.0).unwrap().source().clone();
    let instantaneous = circuit
        .reversible_instantaneous()
        .unwrap()
        .at_resolved_temperature(&law.at_temperature(310.0).unwrap())
        .unwrap();
    assert!(
        circuit
            .rebind_source_candidate(source.clone(), &[], || false)
            .is_err()
    );
    assert!(
        circuit
            .rebind_source_and_resistors_candidate(source.clone(), vec![], &[], || false)
            .is_err()
    );
    assert!(
        circuit
            .rebind_reversible_sources_candidate(source.clone(), instantaneous.clone(), &[], || {
                false
            })
            .is_err()
    );
    assert!(
        circuit
            .rebind_reversible_sources_and_resistors_candidate(
                source.clone(),
                instantaneous.clone(),
                vec![],
                &[],
                || false
            )
            .is_err()
    );
    let refreshed = circuit
        .rebind_thermal_reversible_sources_candidate(source, instantaneous, &[], || false)
        .unwrap();
    for &slot in &circuit.dae().charge_index {
        assert_eq!(
            refreshed.circuit().state()[slot].to_bits(),
            saved[slot].to_bits()
        );
    }
    assert_eq!(
        refreshed.circuit().thermal_reversible_storage(),
        circuit.thermal_reversible_storage()
    );
    assert!(std::ptr::eq(
        refreshed.circuit().thermal_reversible_storage().unwrap(),
        circuit.thermal_reversible_storage().unwrap(),
    ));
    let before = field_state(1.9, 300.0, 0.001).0 + 0.4_f64.powi(2) / 2.0;
    let after = field_state(1.9, 310.0, 0.001).0 + 0.4_f64.powi(2) / (2.0 * (1.0 + 0.0006 * 10.0));
    near(refreshed.parameter_work_j(), after - before, 2e-12);

    // This same card agrees at Tref but supplies a different epsilon slope.
    let inconsistent = ResolvedDebyeCapacitor::resolve(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inconsistent",
        circuit.source().pole_properties(),
        1.0,
        1.0,
    )
    .unwrap();
    let inconsistent =
        ThermalDebyeCapacitor::new(&card, &inconsistent, "T", [250.0, 500.0]).unwrap();
    assert!(ThermalDebyeStorage::new(&inconsistent, &law).is_err());

    let curve = caloric(&card);
    let thermal = thermal(&[&curve], &[0.1]);
    let context = point(300.0);
    let cells = [cell(0, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let mut bound = binding(&poles, &[1.0]);
    bound.maximum_instantaneous_charge_c = 3.0;
    assert!(
        dielectric::step_candidate(
            &network,
            &thermal,
            0.1,
            &circuit,
            &[],
            bound,
            controls(),
            &[0.0],
            || false
        )
        .is_err()
    );
    let mut bound = binding(&poles, &[1.0]);
    bound.maximum_instantaneous_charge_c = 3.0;
    let mut polls = 0;
    assert!(matches!(
        dielectric::step_candidate_with_reversible_storage(
            &network,
            &thermal,
            0.1,
            &circuit,
            &[],
            bound,
            &law,
            controls(),
            &[0.0],
            || {
                polls += 1;
                polls > 20
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(polls > 20);
    assert_eq!(circuit.state(), saved);
    assert_eq!(thermal.time_s(), 0.0);
}
