//! The card row `outlet_mutual free-field` (delayed free-field mutual radiation
//! between compact outlets): where it is admitted, what it binds, that it
//! changes the rendered pressure, and that a card without it is untouched; and
//! `outlet_mutual pritchard`, the same term weighted by the openings' own
//! radiation resistance.
use fs_couple::flue_loop::network::outlet_mutual::OutletMutualForm;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const RECORDER_A1: &str =
    include_str!("../../../data/instruments/descriptions/recorder/recorder-a1.performance");
const LAST_OUTLET: &str = "outlet 21 -0.0863000 0\n";

fn gate() -> CancelGate { CancelGate::new_clock_free() }
/// The alto recorder's a1 card cut to a 1024-sample window with its onset only.
fn recorder_a1() -> String {
    assert!(RECORDER_A1.contains("\noutput 120000 5 "));
    let events = "events\t2\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(RECORDER_A1.contains(events), "the a1 gesture changed; update this test");
    assert!(RECORDER_A1.contains(LAST_OUTLET), "the a1 outlet rows changed; update this test");
    RECORDER_A1
        .replace("\noutput 120000 5 ", "\noutput 1024 1 ")
        .replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2")
}
fn with_row(text: &str, row: &str) -> String {
    text.replace(LAST_OUTLET, &format!("{LAST_OUTLET}{row}\n"))
}
fn parse(text: &str) -> Result<FlueFilePerformance, String> {
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate()).map_err(|e| e.to_string())
}
fn render(text: &str) -> Vec<f64> {
    let mut renderer = parse(text).unwrap().into_renderer();
    let mut pressure = vec![0.0_f64; 1024];
    for block in pressure.chunks_mut(256) { renderer.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    pressure
}

#[test]
fn the_row_binds_every_observed_outlet_in_order_with_its_admission_record() {
    let base = recorder_a1();
    let plain = parse(&base).unwrap();
    assert!(plain.renderer().voice().network().unwrap().outlet_mutual().is_none());
    let file = parse(&with_row(&base, "outlet_mutual free-field")).unwrap();
    let network = file.renderer().voice().network().unwrap();
    let mutual = network.outlet_mutual().expect("attached");
    assert_eq!(mutual.nodes(), [24, 19, 21]);
    let pairs = mutual.pairs();
    assert_eq!(pairs.len(), 3);
    // Primary at the origin; holes at their card centers.
    let distances: Vec<f64> = pairs.iter().map(|p| p.distance_m).collect();
    for (got, want) in distances.iter().zip([0.1199, 0.0863, 0.1199 - 0.0863]) {
        assert!((got - want).abs() < 1e-12, "{distances:?}");
    }
    let record = mutual.record();
    let binding = file.radiation_bindings()[0].load();
    assert_eq!(record.density_kg_m3.to_bits(), binding.density_kg_m3().to_bits());
    assert_eq!(record.sound_speed_m_s.to_bits(), binding.sound_speed_m_s().to_bits());
    assert!((record.time_step_s * 288_000.0 - 1.0).abs() < 1e-15);
    assert_eq!(record.band_hz, 2800.0);
    assert!(record.max_band_error < 0.05 && record.min_pivot_ratio > 0.0, "{record:?}");
    // The plain row carries no weight filter.
    assert_eq!(record.form, OutletMutualForm::FreeField);
    assert_eq!(record.min_weight_half_taps, 0);
    assert_eq!(record.max_free_field_deviation.to_bits(), record.max_band_error.to_bits());
}

#[test]
fn the_pritchard_row_binds_the_same_outlets_and_renders_apart_from_free_field() {
    let base = recorder_a1();
    let file = parse(&with_row(&base, "outlet_mutual pritchard")).unwrap();
    let mutual = file.renderer().voice().network().unwrap().outlet_mutual().expect("attached");
    assert_eq!(mutual.nodes(), [24, 19, 21]);
    assert_eq!(mutual.pairs().len(), 3);
    let record = mutual.record();
    assert_eq!(record.form, OutletMutualForm::Pritchard);
    assert!(record.min_weight_half_taps >= 8, "{record:?}");
    assert!(record.max_band_error < 0.05 && record.min_pivot_ratio > 0.0, "{record:?}");
    assert!(mutual.min_read_lag() >= 1);
    let (plain, free, weighted) = (render(&base), render(&with_row(&base, "outlet_mutual free-field")),
        render(&with_row(&base, "outlet_mutual pritchard")));
    assert_eq!(weighted, render(&with_row(&base, "outlet_mutual pritchard")), "the same bits twice");
    assert_ne!(weighted, plain, "planted negative: the weighted coupling must reach the microphone");
    assert_ne!(weighted, free, "planted negative: the weight must change the coupling");
}

#[test]
fn the_row_changes_the_rendered_pressure_and_its_absence_changes_nothing() {
    let base = recorder_a1();
    let plain = render(&base);
    assert_eq!(plain, render(&base), "a card without the row renders the same bits twice");
    let coupled = render(&with_row(&base, "outlet_mutual free-field"));
    assert_ne!(plain, coupled, "planted negative: the coupling must reach the microphone");
}

#[test]
fn the_row_is_refused_outside_its_scope() {
    let base = recorder_a1();
    assert!(parse(&with_row(&base, "outlet_mutual free-field")).is_ok());
    assert!(parse(&with_row(&base, "outlet_mutual pritchard")).is_ok());
    for (row, reason) in [("outlet_mutual", "free-field"), ("outlet_mutual half-space", "free-field"),
        ("outlet_mutual free-field extra", ""), ("outlet_mutual Pritchard", "pritchard"), ("outlet_mutual pritchard extra", "")] {
        let error = parse(&with_row(&base, row)).err().unwrap_or_else(|| panic!("{row} was admitted"));
        assert!(error.contains(reason), "{row}: {error}");
    }
    let prescribed = with_row(&base, "outlet_mutual free-field").replace("radiation compact-baffled", "radiation prescribed");
    let error = parse(&prescribed).err().expect("refused without compact loads");
    assert!(error.contains("compact-baffled"), "{error}");
    // Only directly after the outlet rows.
    let late = base.replace("mouth_reaction window-cylinder", "outlet_mutual free-field\nmouth_reaction window-cylinder");
    assert!(late != base && parse(&late).is_err());
}
