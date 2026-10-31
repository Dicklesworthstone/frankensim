//! G0/G3: physical breath ramps retain source time, units and native owners.
use super::*;

fn ramp_input() -> String {
    fixture(reference(100000.0), "dry-air 300 100000").replace(
        "boundary 1\nhold 32 -0.1",
        "boundary pressure-temperature 2\nramp 16 10 300 14 300.5\nhold 32 12 300.5",
    )
}

#[test]
fn g3_native_warm_ramp_replays_independently_authored_typed_controls() {
    let text = ramp_input();
    let file = read(&text, 8);
    assert_eq!(file.info().boundary_segments, 2);
    assert_eq!(file.info().boundary_ramps, 1);
    // Independent quarter-Pa and 1/32-K increments. The distinct held pressure
    // proves that the next segment owns the exclusive matching-clock endpoint.
    let controls: Vec<_> = (0..SAMPLES).map(|i| {
        NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: if i < 16 { 10.0 + i as f64 / 4.0 } else { 12.0 },
            inflow_temperature_k: 300.0 + i.min(16) as f64 / 32.0,
        })
    }).collect();
    assert_eq!(file.renderer().drive_history(), controls.as_slice());
    let system = read(&fixture(reference(100000.0), "dry-air 300 100000"), 8)
        .into_renderer().system;
    let mut typed = NsfDgPerformance::new_with_drives(system,
        ApertureObservation::Inlet, RATE, controls, 3,
        Arc::new(CancelGate::new_clock_free())).unwrap();
    let mut expected = vec![0.0; SAMPLES];
    for chunk in expected.chunks_mut(3) { typed.block(chunk).unwrap(); }
    let (observed, native) = render(file, 8);
    let (partitioned, replay) = render(read(&text, 3), 3);
    assert_eq!(observed, expected);
    assert_eq!(observed, partitioned);
    assert_eq!(native.system().zone().state_hash(), typed.system().zone().state_hash());
    assert_eq!(native.system().zone().state_hash(), replay.system().zone().state_hash());
    assert_eq!(native.system().load().state(), typed.system().load().state());
}

#[test]
fn g0_ramps_sample_original_interval_starts_on_a_nonmatching_observation_clock() {
    let text = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("clock 96000 1", "clock 72000 1")
        .replace("boundary 1\nhold 32 -0.1",
            "boundary pressure 3\nramp 8 10 18\nhold 16 9\nramp 32 9 13");
    let file = read(&text, 8);
    let history = file.renderer().drive_history();
    assert_eq!(history.len(), 43);
    assert_eq!(file.info().boundary_ramps, 2);
    let dt = file.renderer().system().load().propagator().dt();
    assert!(!same_clock(dt, 1.0 / 72000.0));
    for (i, control) in history.iter().enumerate() {
        // The events deliberately fall between source intervals, independently
        // of floating-point ties. The oracle uses physical time and Pa/s.
        let time = i as f64 * dt;
        let expected = if time < 8.0 / 72000.0 {
            10.0 + time * 72000.0
        } else if time < 16.0 / 72000.0 { 9.0 }
        else { 9.0 + (time - 16.0 / 72000.0) * 18000.0 };
        let NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(value)) = control else {
            panic!("pressure units lost");
        };
        assert!((value - expected).abs() < 2e-13);
    }
    assert_eq!(history[11], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(9.0)));
}

#[test]
fn g0_ramp_variants_keep_physical_location_and_legacy_holds_keep_bits() {
    for (kind, a, b, aperture, expected) in [
        (0, "-0.2", "0.2", None, NsfDgDrive::Boundary(NsfBoundaryDrive::IncomingCharacteristic(0.0))),
        (1, "10", "14", None, NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(12.0))),
        (2, "10", "14", None, NsfDgDrive::UpstreamPressureRise(12.0)),
        (3, "100000 300", "100004 302", Some(0), NsfDgDrive::Reservoir { pressure_pa: 100002.0, temperature_k: 301.0 }),
        (4, "10 300", "14 302", None, NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 12.0, inflow_temperature_k: 301.0 })),
        (5, "100000 300", "100004 302", Some(1), NsfDgDrive::ReservoirAperture {
            aperture: 1, pressure_pa: 100002.0, temperature_k: 301.0 }),
    ] {
        let mut first = Row { fields: a.split_ascii_whitespace(), line: 1 };
        let mut last = Row { fields: b.split_ascii_whitespace(), line: 1 };
        let segment = BoundarySegment { start: 4, end: 12,
            first: read_drive(&mut first, kind, aperture, [250.0, 600.0], None).unwrap(),
            last: Some(read_drive(&mut last, kind, aperture, [250.0, 600.0], None).unwrap()) };
        first.end().unwrap(); last.end().unwrap();
        assert_eq!(segment.at(8.0).unwrap(), expected);
    }
    let held = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("hold 32 -0.1", "hold 32 -0");
    let file = read(&held, 8);
    assert_eq!(file.info().boundary_ramps, 0);
    for control in file.renderer().drive_history() {
        let NsfDgDrive::Boundary(NsfBoundaryDrive::IncomingCharacteristic(value)) = control else {
            panic!("held units lost");
        };
        assert_eq!(value.to_bits(), (-0.0_f64).to_bits());
    }
    let indexed = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("initial rest",
            "reservoir_aperture 1 low reference 1\nreservoir_cell 4 0 1\nreservoir_aperture 1 high reference 1\nreservoir_cell 7 2 1\ninitial rest")
        .replace("boundary 1\nhold 32 -0.1",
            "boundary aperture 1 1\nramp 32 100000 300 100004 302");
    let file = read(&indexed, 8);
    assert_eq!(file.renderer().drive_history()[16],
        NsfDgDrive::ReservoirAperture { aperture: 1, pressure_pa: 100002.0, temperature_k: 301.0 });
}

#[test]
fn g0_unsampled_ramp_endpoints_and_malformed_histories_refuse() {
    let text = ramp_input();
    for replacement in ["ramp 16 10 300 14 601", "ramp 16 10 249 14 300",
        "ramp 16 10 300 14", "ramp 16 10 300 14 300 extra",
        "ramp 16 10 300 NaN 300", "ramp 0 10 300 14 300",
        "ramp 33 10 300 14 300", "step 16 10 300 14 300"] {
        let invalid = text.replace("ramp 16 10 300 14 300.5", replacement);
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err(), "accepted {replacement}");
    }
    let reservoir = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("initial rest", "reservoir_aperture 1 low reference 1\nreservoir_cell 6 0 1\ninitial rest")
        .replace("boundary 1\nhold 32 -0.1", "boundary reservoir 1\nramp 32 100000 300 100004 601");
    assert!(ResolvedDgFilePerformance::from_bytes(reservoir.as_bytes(), 8,
        Arc::new(CancelGate::new_clock_free())).is_err());
}

#[test]
fn g0_warm_ramp_example_uses_original_radiation_and_gas_observation() {
    let file = read(include_str!("../../../../../examples/resolved-dg-warm-ramp.performance"), 8);
    assert_eq!(file.info().boundary_ramps, 3);
    assert_eq!(file.renderer().drive_history()[0],
        NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 0.0, inflow_temperature_k: 300.0 }));
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
}
