//! Authored observation policy reaches playback without changing the source.
use super::*;

const FILTER: &str = "observation_filter bandlimited 2000 8000 128 2000000 65536 65536\n";

fn input() -> String {
    fixture(reference(100000.0), "dry-air 300 100000")
        .replace("output 32 20", "output 256 20")
        .replace("hold 32 -0.1", "hold 256 -0.1")
}

#[test]
fn authored_filter_changes_observation_but_preserves_physical_history() {
    let original = input();
    let text = original.replace("gas dry-air", &format!("{FILTER}gas dry-air"));
    let source = read(&text, 19);
    let filter = source.renderer().bandlimited_observation().unwrap();
    let report = filter.report();
    assert_eq!(filter.policy().passband_hz, 2000.0);
    assert_eq!(filter.policy().stopband_hz, 8000.0);
    assert_eq!(report.channels, 1);
    assert!(report.latency_s > 0.0);
    assert_eq!(
        report.input_period_s,
        source.renderer().system().load().propagator().dt()
    );
    let (filtered, state) = render(source, 19);
    let (partitioned, other) = render(read(&text, 7), 7);
    let (unfiltered, original_state) = render(read(&original, 19), 19);
    assert_eq!(filtered, partitioned);
    assert_ne!(filtered, unfiltered);
    assert_eq!(
        state.system().zone().state_hash(),
        original_state.system().zone().state_hash()
    );
    assert_eq!(
        state.system().zone().state_hash(),
        other.system().zone().state_hash()
    );
    assert_eq!(
        state.system().load().state(),
        original_state.system().load().state()
    );
    assert_eq!(state.drive_history(), original_state.drive_history());
    assert_eq!(
        state.system().accepted_steps(),
        original_state.system().accepted_steps()
    );
}

#[test]
fn unavailable_observation_band_and_incomplete_budget_refuse_the_file() {
    let original = input();
    for policy in [
        "bandlimited 2000 60000 128 2000000 65536 65536",
        "bandlimited 8000 2000 128 2000000 65536 65536",
        "bandlimited 2000 8000 128 1 65536 65536",
        "bandlimited 2000 8000 128 2000000 1 65536",
        "bandlimited 2000 8000 128 2000000 65536 1",
    ] {
        let text = original.replace(
            "gas dry-air",
            &format!("observation_filter {policy}\ngas dry-air"),
        );
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                text.as_bytes(),
                19,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err(),
            "{policy}"
        );
    }
    let repeated = original.replace("gas dry-air", &format!("{FILTER}{FILTER}gas dry-air"));
    assert!(
        ResolvedDgFilePerformance::from_bytes(
            repeated.as_bytes(),
            19,
            Arc::new(CancelGate::new_clock_free())
        )
        .is_err()
    );
}
