//! Source-aperture declarations preserve original fluid and acoustic owners.
use super::*;

fn input(upstream: bool, background: f64) -> String {
    let gas = reference(100000.0);
    let mut text = fixture(gas, "dry-air 300 100000")
        .replace("characteristic 0 high 0.05", &format!(
            "characteristic 0 high 0.05\nsource_patch 1 2 1 2 pressure {background}"))
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 2");
    if upstream {
        let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
        let radius = (dx * dx / core::f64::consts::PI).sqrt();
        // Both nested roots spend this ONE allowance. The one-cell source
        // aperture needs more trials than the original full-face upstream case.
        text = text.replace("coupling 64 6912", "coupling 256 27648")
            .replace("observation inlet-pressure", &format!(
                "upstream_dg dg_tube 0.03 {radius:.17e} closed 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
            .replace("boundary pressure 1", "boundary upstream-pressure 1");
    }
    text
}

#[test]
fn native_source_aperture_drives_only_selected_source_cells_across_callback_partitions() {
    {
        let text = input(false, -0.25);
        let file = read(&text, 8);
        let system = file.renderer().system();
        assert!(system.coupled_patch().is_none());
        assert_eq!(system.source_patch(), Some(NsfCharacteristicPatch { start: [1, 1], end: [2, 2] }));
        assert_eq!(system.area_m2(), 9.0 * system.zone().law().dx_m.powi(2));
        let (pressure, complete) = render(file, 8);
        let (partitioned, parts) = render(read(&text, 3), 3);
        assert_eq!(pressure, partitioned);
        assert!(pressure.iter().any(|p| p.abs() > 1e-8));
        assert_eq!(complete.system().zone().state_hash(), parts.system().zone().state_hash());
        assert_eq!(complete.system().load().state(), parts.system().load().state());
        let frame = complete.last_frame().unwrap();
        let patch = frame.source_patch.unwrap();
        assert!(frame.coupled_patch.is_none());
        assert_eq!(frame.inlet_flow_m3_s, frame.faces[1].normal_volume_flow_m3_s);
        for trace in complete.system().zone().open_boundary_traces(false).unwrap() {
            if trace.cell[1] != 1 || trace.cell[2] != 1 {
                assert!((trace.pressure_pa - 99999.75).abs() < 1e-7);
            }
        }
        assert!(frame.upstream.is_none());
        assert!((patch.mean_acoustic_pressure_pa - 2.0).abs() < 1e-7);
        let (ambient, _) = render(read(&input(false, 0.0), 8), 8);
        assert_ne!(pressure, ambient, "source remainder pressure must reach actual motion");
    }
}

#[test]
fn reflected_source_aperture_completes_both_duct_histories_across_callback_partitions() {
    // G0/G3: retain the complete 32-sample, 2 Pa, closed 0.03 m upstream
    // history, including the returning wave that previously refused interval
    // 18. Linear pressure/characteristic controls must retain their reference
    // advective donor without weakening the original root or work allowances.
    let text = input(true, -0.25);
    let mut renderer = read(&text, 1).into_renderer();
    assert!(renderer.system().coupled_patch().is_none());
    assert_eq!(renderer.system().source_patch(),
        Some(NsfCharacteristicPatch { start: [1, 1], end: [2, 2] }));
    assert_eq!(renderer.system().area_m2(), 9.0 * renderer.system().zone().law().dx_m.powi(2));
    let snapshot = |renderer: &NsfDgPerformance<D3q27>| {
        let system = renderer.system();
        let zone = system.zone();
        (
            (zone.state_hash(), zone.steps(), zone.cumulative_ledger(),
                zone.boundary_drives(), zone.characteristic_patch_drive(),
                zone.last_boundary_exchange().copied(), zone.last_patch_exchange().copied()),
            (system.load().state().clone(), system.upstream_load().unwrap().state().clone(),
                system.accepted_steps()),
            (renderer.samples_rendered(), renderer.remaining_samples(),
                renderer.last_frame().map(|frame| (frame.step, frame.inlet_pressure_pa, frame.source_patch,
                    frame.upstream.as_ref().map(|source| (source.inlet_pressure_pa, source.inlet_flow_m3_s))))),
        )
    };
    let mut accepted = Vec::new();
    let mut history = Vec::new();
    for step in 1..=SAMPLES {
        let mut sample = [f64::NAN];
        renderer.block(&mut sample).unwrap();
        accepted.push(sample[0]);
        let frame = renderer.last_frame().unwrap();
        assert_eq!(frame.step, step as u64);
        assert_eq!(renderer.system().accepted_steps(), step as u64);
        assert_eq!(renderer.system().zone().steps(), step);
        assert_eq!(renderer.samples_rendered(), step as u64);
        let patch = frame.source_patch.unwrap();
        let source = frame.upstream.as_ref().unwrap();
        assert!(frame.coupled_patch.is_none());
        assert_eq!(frame.inlet_flow_m3_s, frame.faces[1].normal_volume_flow_m3_s);
        assert_eq!(source.inlet_flow_m3_s, patch.normal_volume_flow_m3_s);
        assert_ne!(source.inlet_flow_m3_s, frame.faces[0].normal_volume_flow_m3_s);
        assert!((patch.mean_acoustic_pressure_pa - source.inlet_pressure_pa - 2.0).abs() < 1e-6);
        for trace in renderer.system().zone().open_boundary_traces(false).unwrap() {
            if trace.cell[1] != 1 || trace.cell[2] != 1 {
                assert!((trace.pressure_pa - 99999.75).abs() < 1e-7);
            }
        }
        history.push(snapshot(&renderer));
    }
    assert_eq!(accepted.len(), SAMPLES);
    assert!(accepted.iter().all(|p| p.is_finite()));
    assert!(accepted.iter().any(|p| p.abs() > 1e-8));
    assert_eq!(renderer.remaining_samples(), 0);

    for block in [3, 8] {
        let mut partitioned = read(&text, block).into_renderer();
        let mut pressure = [f64::NAN; SAMPLES];
        let mut rendered = 0;
        for chunk in pressure.chunks_mut(block) {
            partitioned.block(chunk).unwrap();
            rendered += chunk.len();
            // Compare both complete original duct states and the fluid state
            // at every callback endpoint, not only their final pressure sample.
            assert_eq!(snapshot(&partitioned), history[rendered - 1]);
        }
        assert_eq!(pressure.as_slice(), accepted.as_slice());
        assert_eq!(partitioned.remaining_samples(), 0);
    }
}

#[test]
fn native_source_aperture_refuses_duplicate_roles_bad_controls_and_wrong_upstream_area() {
    let text = input(true, 0.0);
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (dx * dx / core::f64::consts::PI).sqrt();
    for invalid in [
        text.replace("source_patch 1 2 1 2 pressure 0", "source_patch 1 2 1 2 pressure NaN"),
        text.replace("source_patch 1 2 1 2 pressure 0", "source_patch 2 1 1 2 pressure 0"),
        text.replace("source_patch 1 2 1 2 pressure 0", "source_patch 1 4 1 2 pressure 0"),
        text.replace("source_patch 1 2 1 2 pressure 0", "source_patch 1 2 1 2 pressure 0\ncoupled_patch 1 2 1 2 pressure 0"),
        text.replace("source_patch 1 2 1 2 pressure 0", "coupled_patch 1 2 1 2 pressure 0\nsource_patch 1 2 1 2 pressure 0"),
        text.replace(&format!("upstream_dg dg_tube 0.03 {radius:.17e}"),
            &format!("upstream_dg dg_tube 0.03 {:.17e}", 3.0 * radius)),
    ] {
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err(), "{invalid}");
    }
    let incoming = input(false, 0.0).replace("source_patch 1 2 1 2 pressure 0",
        "source_patch 1 2 1 2 characteristic 0.001");
    assert_eq!(read(&incoming, 8).renderer().system().patch_background_drive(),
        Some(NsfBoundaryDrive::IncomingCharacteristic(0.001)));
}

