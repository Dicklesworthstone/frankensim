use super::*;

fn gas() -> GasState {
    GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap()
}
fn apertures() -> [CircularOpening; 2] {
    [
        CircularOpening { center_m: [0.0, 0.0], radius_m: 0.002 },
        CircularOpening { center_m: [0.015, 0.0], radius_m: 0.0015 },
    ]
}
fn pistons(openings: &[CircularOpening]) -> Vec<BaffledPiston> {
    openings.iter().map(|a| BaffledPiston {
        center_m: a.center_m, radius_m: a.radius_m,
    }).collect()
}
fn options() -> BaffledMimoBakeOptions {
    BaffledMimoBakeOptions {
        frequency_samples: 65, orders: [4, 8, 12], fit_iterations: 20,
        maximum_relative_error: 0.01, max_fit_memory_bytes: 64 * 1024 * 1024,
        causality: CausalityOptions { subtractions: 2, tail_bound_ratio: 100.0,
            relative_data_tolerance: 1e-8 },
    }
}
fn limits() -> MimoRadiationLimits {
    MimoRadiationLimits { max_states: 64, max_memory_bytes: 1024 * 1024,
        sampled_passivity_points: 257 }
}

#[test]
fn source_uses_acoustic_units_causal_phase_and_actual_mutual_geometry() {
    let gas = gas();
    let opening = apertures();
    let frequencies = [100.0, 1100.0, 3800.0];
    let gate = CancelGate::new_clock_free();
    let rows = sample_rows(&pistons(&opening), &gas, &frequencies, &gate).unwrap();
    for (&frequency, row) in frequencies.iter().zip(&rows) {
        for (port, a) in opening.iter().enumerate() {
            // This independent acoustic owner already returns p/Q. An extra
            // area factor or missing phasor conjugation cannot pass this check.
            let (resistance, reactance) = fs_phs::baffled_piston_impedance(
                gas.density, gas.sound_speed, a.radius_m,
                core::f64::consts::TAU * frequency).unwrap();
            let expected = C64::new(resistance, -reactance);
            assert!((row[port * 2 + port] - expected).abs() <= 4e-13 * expected.abs());
            assert!(row[port * 2 + port].re > 0.0 && row[port * 2 + port].im > 0.0);
        }
        assert!(row[1].abs() > 0.01 * row[0].abs());
        assert!((row[1] - row[2]).abs() < 1e-12 * row[1].abs());
    }
    let mut farther = opening;
    farther[1].center_m[0] *= 2.0;
    let distant = sample_rows(&pistons(&farther), &gas, &frequencies, &gate).unwrap();
    assert_eq!(rows[0][0], distant[0][0]);
    assert_eq!(rows[0][3], distant[0][3]);
    assert!(distant[0][1].abs() < 0.6 * rows[0][1].abs());
    assert!((distant[2][1] - rows[2][1]).abs() > 0.1 * rows[2][1].abs());
}

#[test]
fn source_recomputes_all_entries_from_pressure_and_mixture_data() {
    let ambient = gas();
    let reduced = GasState::try_new_moist_air(293.15, 81_060.0, 0.0).unwrap();
    let humid = GasState::try_new_moist_air(293.15, 101_325.0, 0.8).unwrap();
    let opening = apertures();
    let pistons = pistons(&opening);
    let gate = CancelGate::new_clock_free();
    let row = sample_rows(&pistons, &ambient, &[3800.0], &gate).unwrap().remove(0);
    let low = sample_rows(&pistons, &reduced, &[3800.0], &gate).unwrap().remove(0);
    let wet = sample_rows(&pistons, &humid, &[3800.0], &gate).unwrap().remove(0);
    assert_eq!(ambient.sound_speed, reduced.sound_speed);
    let ratio = reduced.density / ambient.density;
    for (a, b) in row.iter().zip(&low) {
        assert!((*b - a.scale(ratio)).abs() < 4e-12 * a.abs());
    }
    // Composition also changes c and therefore the mutual propagation phase;
    // a density-only update cannot reproduce this original physical sample.
    let density_only = row[1].scale(humid.density / ambient.density);
    assert!((wet[1] - density_only).abs() > 1e-5 * wet[1].abs());
    assert_ne!(declaration(&ambient, &opening, 8000, [20.0, 4000.0]),
        declaration(&humid, &opening, 8000, [20.0, 4000.0]));
}

#[test]
fn physical_bake_retains_three_orders_original_source_and_held_out_accuracy() {
    let exterior = gas();
    let openings = apertures();
    let opts = options();
    let gate = CancelGate::new_clock_free();
    let (source, report) = bake_baffled_mimo(exterior, &openings, 8000,
        [20.0, 4000.0], opts, limits(), &gate).unwrap();
    assert_eq!(source.exterior, exterior);
    assert_eq!(source.apertures, openings);
    assert_eq!(source.band_hz, [20.0, 4000.0]);
    assert!(source.impedance.is_stable());
    assert!(source.impedance.e.iter().all(|e| *e == 0.0));
    assert_eq!(report.rungs.iter().map(|r| r.order).collect::<Vec<_>>(), opts.orders);
    assert_eq!(report.held_out_samples, opts.frequency_samples - 1);
    assert!(report.maximum_fit_relative_error <= opts.maximum_relative_error);
    assert!(report.maximum_held_out_relative_error <= opts.maximum_relative_error);
    assert!(report.maximum_causality_ratio <= 1.0);
    assert!(report.continuous.passive && report.sampled.passive);
    let selected = report.rungs.iter().find(|r| r.order == report.selected_order).unwrap();
    assert!(selected.refusal.is_none() && selected.continuous_passive && selected.sampled_passive);

    // These points belong to neither identification nor its geometric midpoint
    // validation grid. Compare with new physical samples, including mutuals.
    let frequencies = [93.0, 713.0, 2671.0, 3983.0];
    let data = sample_rows(&pistons(&openings), &exterior, &frequencies, &gate).unwrap();
    for (&frequency, expected) in frequencies.iter().zip(&data) {
        let predicted = source.impedance.eval_iw(core::f64::consts::TAU * frequency);
        let reference_norm = det::sqrt(expected.iter().map(|z| z.norm_sq()).sum());
        let error = det::sqrt(predicted.iter().zip(expected).map(|(a, b)| (*a - *b).norm_sq()).sum());
        assert!(error <= opts.maximum_relative_error * reference_norm);
        assert!(predicted[1].abs() > 0.005 * predicted[0].abs());
    }
}

#[test]
fn source_refuses_short_band_over_budget_geometry_and_cancellation_before_fitting() {
    let gas = gas();
    let openings = apertures();
    let opts = options();
    let runtime = limits();
    let admission = |openings: &[CircularOpening], band, options, limits| {
        admit(&gas, openings, 8000, band, options, limits)
    };
    assert!(matches!(admission(&openings, [20.0, 3999.0], opts, runtime),
        Err(BaffledMimoBakeError::Invalid(_))));
    let mut overlap = openings;
    overlap[1].center_m = [0.003, 0.0];
    assert!(admission(&overlap, [20.0, 4000.0], opts, runtime).is_err());
    let mut too_small = opts;
    too_small.max_fit_memory_bytes = 1024;
    assert!(matches!(admission(&openings, [20.0, 4000.0], too_small, runtime),
        Err(BaffledMimoBakeError::Memory { kind: "fit scratch", .. })));
    let mut too_small = runtime;
    too_small.max_memory_bytes = 1024;
    assert!(matches!(admission(&openings, [20.0, 4000.0], opts, too_small),
        Err(BaffledMimoBakeError::Memory { kind: "runtime payload", .. })));
    too_small = runtime;
    too_small.max_states = 23;
    assert!(admission(&openings, [20.0, 4000.0], opts, too_small).is_err());
    let mut invalid_ladder = opts;
    invalid_ladder.orders = [4, 4, 12];
    assert!(admission(&openings, [20.0, 4000.0], invalid_ladder, runtime).is_err());
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(bake_baffled_mimo(gas, &openings, 8000, [20.0, 4000.0],
        opts, runtime, &gate), Err(BaffledMimoBakeError::Cancelled)));
}

#[test]
fn full_seventeen_opening_radiation_fits_the_existing_memory_and_state_envelopes() {
    let openings: Vec<_> = (0..17).map(|i| CircularOpening {
        center_m: [0.015 * f64::from(i), 0.0], radius_m: 0.001,
    }).collect();
    let options = options();
    let mut runtime = limits();
    runtime.max_states = 204;
    admit(&gas(), &openings, 8000, [20.0, 4000.0], options, runtime).unwrap();
    let required = fit_scratch_bytes(openings.len(), options);
    assert!(required < 32 * 1024 * 1024, "{required}");
    let mut short = options;
    short.max_fit_memory_bytes = required - 1;
    assert!(matches!(admit(&gas(), &openings, 8000, [20.0, 4000.0], short, runtime),
        Err(BaffledMimoBakeError::Memory { kind: "fit scratch", .. })));
    runtime.max_states = 203;
    assert!(admit(&gas(), &openings, 8000, [20.0, 4000.0], options, runtime).is_err());

    // Source rows still come from every original physical self/mutual pair.
    let rows = sample_rows(&pistons(&openings), &gas(), &[1000.0],
        &CancelGate::new_clock_free()).unwrap();
    let matrix = &rows[0];
    assert_eq!(matrix.len(), 17 * 17);
    for i in 0..17 {
        assert!(matrix[i * 17 + i].re > 0.0);
        for j in 0..17 {
            assert!(matrix[i * 17 + j].abs() > 0.0);
            assert!((matrix[i * 17 + j] - matrix[j * 17 + i]).abs()
                < 1e-12 * matrix[i * 17 + j].abs());
        }
    }
}
