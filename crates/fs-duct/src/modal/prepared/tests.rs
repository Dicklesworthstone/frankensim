use super::*;
use fs_material::gas::GasSpec;

fn duct() -> Duct {
    Duct { segments: vec![
        Segment::Cylinder { radius: 0.006, length: 0.15 },
        Segment::Cone { inlet_radius: 0.006, outlet_radius: 0.010, length: 0.12 },
        Segment::Cylinder { radius: 0.010, length: 0.05 },
        Segment::Cone { inlet_radius: 0.010, outlet_radius: 0.007, length: 0.09 },
    ] }
}

fn gas(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}

fn bits(values: &[C64]) -> Vec<(u64, u64)> {
    values.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
}

fn same_response(a: &ModalResponse, b: &ModalResponse) {
    assert_eq!(a.omega.to_bits(), b.omega.to_bits());
    assert_eq!(a.n_modes, b.n_modes);
    assert_eq!(bits(&a.impedance_matrix), bits(&b.impedance_matrix));
    assert_eq!(bits(&[a.plane_impedance]), bits(&[b.plane_impedance]));
    assert_eq!(a.min_shear_number.to_bits(), b.min_shear_number.to_bits());
    assert_eq!(a.mouth_ka.to_bits(), b.mouth_ka.to_bits());
    assert_eq!(a.cutoffs_at_throat_hz, b.cutoffs_at_throat_hz);
    assert_eq!(a.steps_used, b.steps_used);
    assert_eq!(a.tail_bound_estimate, b.tail_bound_estimate);
}

#[test]
fn prepared_geometry_reuses_the_exact_original_roots_and_junctions() {
    let source = duct();
    for modes in [1, 3] {
        for slices in [1, 2] {
            let prepared = PreparedModalDuct::new(&source, modes, slices).unwrap();
            assert!(core::ptr::eq(prepared.source(), &source));
            assert_eq!(prepared.n_modes(), modes);
            let stations = stations_of(&source, slices).unwrap();
            let gammas = j1_roots(modes - 1).unwrap();
            assert_eq!(prepared.gammas, gammas);
            assert_eq!(prepared.stations.len(), stations.len());
            for (cached, original) in prepared.stations.iter().zip(&stations) {
                assert_eq!(cached.radius.to_bits(), original.radius.to_bits());
                assert_eq!(cached.length.to_bits(), original.length.to_bits());
            }
            let mut expansions = 0;
            let mut contractions = 0;
            for (pair, cached) in stations.windows(2).zip(&prepared.junctions) {
                let a = pair[0].radius;
                let b = pair[1].radius;
                if (a - b).abs() <= 1e-12 * a.max(1e-300) {
                    assert!(cached.is_none());
                    continue;
                }
                let cached = cached.as_ref().unwrap();
                let expected = if a < b {
                    expansions += 1;
                    junction_projection(a, b, &gammas, modes)
                } else {
                    contractions += 1;
                    junction_projection(b, a, &gammas, modes)
                };
                assert_eq!(cached.expansion, a < b);
                assert_eq!(bits(&cached.projection), bits(&expected));
                assert_eq!(bits(&cached.transpose), bits(&transpose(&expected, modes)));
            }
            assert!(expansions > 0 && contractions > 0);
        }
    }
}

#[test]
fn prepared_sweep_preserves_frequency_loss_and_live_material_state() {
    let source = duct();
    let start = std::time::Instant::now();
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let preparation_time = start.elapsed();
    let original = gas(293.15, 101325.0);
    let changed_pressure = gas(293.15, 95000.0);
    let changed_temperature = gas(313.15, 101325.0);
    let changed_composition = GasState::try_new_moist_air(293.15, 101325.0, 0.75).unwrap();
    let load = C64::new(2.0e6, -3.0e5);
    let mut responses = Vec::new();
    let mut reused_time = std::time::Duration::ZERO;
    let mut repeated_time = std::time::Duration::ZERO;
    for material in [&original, &changed_pressure, &changed_temperature, &changed_composition] {
        let mut sweep = Vec::new();
        for frequency in [173.0, 620.0, 1310.0] {
            let omega = core::f64::consts::TAU * frequency;
            for loss in [LossModel::Lossless, LossModel::Bessel] {
                let start = std::time::Instant::now();
                let cached = geometry.input_impedance_load(material, omega, loss, load).unwrap();
                reused_time += start.elapsed();
                let start = std::time::Instant::now();
                let fresh = mm_input_impedance_load(&source, material, omega, loss, load, 3, 1).unwrap();
                repeated_time += start.elapsed();
                same_response(&cached, &fresh);
                sweep.extend(bits(&cached.impedance_matrix));
            }
        }
        responses.push(sweep);
    }
    assert_ne!(responses[0], responses[1], "pressure must not be cached with geometry");
    assert_ne!(responses[0], responses[2], "temperature must not be cached with geometry");
    assert_ne!(responses[0], responses[3], "mixture composition must not be cached with geometry");
    // Report measured work without making a machine-dependent timing assertion.
    println!("24 modal queries: geometry preparation {preparation_time:?}, reused {reused_time:?}, repeated geometry {repeated_time:?}; exact response bits preserved");
}

#[test]
fn prepared_analytic_tabulated_and_full_modal_loads_match_one_frequency_queries() {
    let source = duct();
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let material = gas(293.15, 101325.0);
    let exterior = GasState::try_new_moist_air(303.15, 101325.0, 0.75).unwrap();
    let omega = 2400.0;
    for termination in [Termination::Closed, Termination::IdealOpen, Termination::FlangedOpen] {
        let cached = geometry.input_impedance(&material, omega, LossModel::Bessel, termination).unwrap();
        let fresh = mm_input_impedance(&source, &material, omega, LossModel::Bessel, termination, 3, 1).unwrap();
        same_response(&cached, &fresh);
    }
    let table = crate::TabulatedLoad::try_new(
        vec![1000.0, 4000.0],
        vec![C64::new(2.0e6, -3.0e5), C64::new(3.0e6, -5.0e5)],
        "manufactured passive load, not measured data", 0.0,
    ).unwrap();
    let cached = geometry.input_impedance_tabulated(&material, omega, LossModel::Bessel, &table).unwrap();
    let fresh = mm_input_impedance_tabulated(&source, &material, omega, LossModel::Bessel, &table, 3, 1).unwrap();
    same_response(&cached, &fresh);
    assert!(geometry.input_impedance_tabulated(&material, 500.0, LossModel::Bessel, &table).is_err());
    let same_gas = geometry.input_impedance_flanged(&material, &material, omega, LossModel::Bessel).unwrap();
    let changed_gas = geometry.input_impedance_flanged(&material, &exterior, omega, LossModel::Bessel).unwrap();
    let fresh = mm_input_impedance_flanged(&source, &material, &exterior, omega, LossModel::Bessel, 3, 1).unwrap();
    same_response(&changed_gas, &fresh);
    assert_ne!(bits(&same_gas.impedance_matrix), bits(&changed_gas.impedance_matrix));
    assert!(geometry.input_impedance_flanged(&material, &gas(293.15, 95000.0), omega, LossModel::Bessel).is_err());
}

#[test]
fn prepared_construction_and_queries_keep_refusals_and_cancellation() {
    let source = duct();
    for modes in [0, MAX_MODES + 1] {
        assert!(PreparedModalDuct::new(&source, modes, 1).is_err());
    }
    assert!(PreparedModalDuct::new(&Duct { segments: vec![] }, 1, 1).is_err());
    let invalid = Duct { segments: vec![Segment::Cylinder { radius: -0.01, length: 1.0 }] };
    assert!(PreparedModalDuct::new(&invalid, 1, 1).is_err());
    let unsupported = Duct { segments: vec![Segment::ToneHole {
        hole_radius: 0.001, chimney_height: 0.002, bore_radius: 0.006,
        state: crate::HoleState::Open,
    }] };
    assert!(PreparedModalDuct::new(&unsupported, 1, 1).is_err());
    assert!(matches!(PreparedModalDuct::new_controlled(&source, 3, 1, || true),
        Err(DuctError::Cancelled { completed_steps: 0 })));
    let mut calls = 0;
    assert!(matches!(PreparedModalDuct::new_controlled(&source, 3, 1, || {
        calls += 1;
        calls == 3
    }), Err(DuctError::Cancelled { completed_steps: 1 })));
    let prepared = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let material = gas(293.15, 101325.0);
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.input_impedance(&material, omega, LossModel::Bessel, Termination::Closed).is_err());
    }
    for load in [C64::new(-1.0, 0.0), C64::new(1.0, f64::NAN)] {
        assert!(prepared.input_impedance_load(&material, 1000.0, LossModel::Bessel, load).is_err());
    }
    let mut incoherent = material;
    incoherent.density *= 2.0;
    assert!(prepared.input_impedance(&incoherent, 1000.0, LossModel::Bessel, Termination::Closed).is_err());
}
