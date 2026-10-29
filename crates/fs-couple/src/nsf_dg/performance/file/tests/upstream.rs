//! G0/G3: source declarations reach both original DG participants and pressure.
use super::*;

fn input(termination: &str, length: f64) -> String {
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    fixture(gas,"dry-air 300 100000")
        .replace("coupling 64 6912","coupling 128 13824")
        .replace("observation inlet-pressure",&format!(
            "upstream_dg dg_tube {length} {radius:.17e} {termination} 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
        .replace("boundary 1","boundary upstream-pressure 1")
        .replace("hold 32 -0.1","hold 32 2")
}

#[test]
fn g3_native_upstream_data_changes_the_same_pressure_and_keeps_block_identity() {
    let text = input("closed", 0.03);
    let file = read(&text, 8);
    assert!(file.renderer().system().upstream_load().is_some());
    assert!(file.upstream_wall_layout().is_some());
    assert!(
        file.renderer()
            .drive_history()
            .iter()
            .all(|&x| x == NsfDgDrive::UpstreamPressureRise(2.0))
    );
    let (a, owner) = render(file, 8);
    let (b, parts) = render(read(&text, 3), 3);
    assert_eq!(a, b);
    assert_eq!(
        owner.system().zone().state_hash(),
        parts.system().zone().state_hash()
    );
    assert_eq!(
        owner.system().upstream_load().unwrap().state(),
        parts.system().upstream_load().unwrap().state()
    );
    let (different, other) = render(read(&input("matched", 0.03), 8), 8);
    assert_ne!(a, different);
    assert_ne!(owner.system().load().state(), other.system().load().state());
    let (length, _) = render(read(&input("closed", 0.08), 8), 8);
    assert_ne!(a, length);
    let final_frame = owner.last_frame().unwrap();
    let upstream = final_frame.upstream.as_ref().unwrap();
    assert!(
        (final_frame.faces[upstream.face].mean_acoustic_pressure_pa
            - upstream.inlet_pressure_pa
            - 2.0)
            .abs()
            < 1e-6
    );
}

#[test]
fn g0_native_drive_location_and_original_source_domains_are_explicit() {
    let text = input("closed", 0.03);
    for invalid in [
        text.replace("boundary upstream-pressure", "boundary pressure"),
        text.replace("upstream_dg dg_tube", "upstream_dg guessed-tube"),
        text.replace("upstream_dg dg_tube 0.03", "upstream_dg dg_tube -0.03"),
        fixture(reference(100000.0), "dry-air 300 100000")
            .replace("boundary 1", "boundary upstream-pressure 1"),
    ] {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                invalid.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err()
        );
    }
}

#[test]
fn canonical_two_sided_fixture_reaches_original_external_microphone() {
    let file = read(
        include_str!("../../../../../examples/resolved-dg-upstream.performance"),
        8,
    );
    assert!(file.renderer().system().upstream_load().is_some());
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
    assert_eq!(
        file.renderer().drive_history()[0],
        NsfDgDrive::UpstreamPressureRise(20.0)
    );
    let mut renderer = file.into_renderer();
    renderer.block(&mut [0.0; 8]).unwrap();
    assert_eq!(
        renderer
            .system()
            .upstream_load()
            .unwrap()
            .state()
            .step_count(),
        8
    );
    assert!(renderer.last_frame().unwrap().upstream.is_some());
}
