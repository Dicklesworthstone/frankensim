use super::*;

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 4800 20000\n\
    ambient 293.15 101325 0\n";
const TAIL: &str = "listener_m 1\ntermination flanged\nsegments 1\ncylinder 0.0073 0.5\n\
    compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t2.8e3\n";
const SIX: &str = "reed 0.0004 0.013 6000 3.8e-6 1400 0.077";

fn source(reed_row: &str) -> String {
    format!("{HEAD}{reed_row}\n{TAIL}")
}

fn refusal(reed_row: &str) -> &'static str {
    match ReedPerformance::from_bytes(source(reed_row).as_bytes(), 64) {
        Err(ReedPerformanceError::Input { what, .. }) => what,
        other => panic!(
            "expected an input refusal for {reed_row:?}, got {:?}",
            other.err()
        ),
    }
}

fn pressure(reed_row: &str) -> Vec<f64> {
    let mut performance = ReedPerformance::from_bytes(source(reed_row).as_bytes(), 64)
        .expect("a valid reed row must load");
    let mut out = vec![0.0; 4800];
    for chunk in out.chunks_mut(64) {
        performance
            .block(chunk)
            .expect("the short window must render");
    }
    out
}

#[test]
fn a_four_field_lay_loads_and_the_six_field_row_is_untouched() {
    let plain = ReedPerformance::from_bytes(source(SIX).as_bytes(), 64).expect("six fields");
    let with_lay =
        ReedPerformance::from_bytes(source(&format!("{SIX} 7.9e6 2 2e-4 6.47")).as_bytes(), 64)
            .expect("six fields and a lay");
    assert!(plain.info().massive_reed && with_lay.info().massive_reed);
    assert_ne!(plain.info().input_hash, with_lay.info().input_hash);
}

#[test]
fn a_malformed_lay_is_refused_at_its_field() {
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4")),
        "missing field",
        "a lay needs all four fields"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4 6.47 1")),
        "unexpected extra field"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 0.5 2e-4 6.47")),
        "lay exponent must be at least 1 and its loss nonnegative"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4 -1")),
        "lay exponent must be at least 1 and its loss nonnegative"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 4e-4 6.47")),
        "lay onset must lie between zero and the reed's rest opening"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 -1e-4 6.47")),
        "lay onset must lie between zero and the reed's rest opening"
    );
    assert_eq!(
        refusal(&format!("{SIX} 0 2 2e-4 6.47")),
        "physical dimension/scale must be positive"
    );
    assert_eq!(
        refusal("reed 0.0004 0.013 6000 0 0 0.35 7.9e6 2 2e-4 6.47"),
        "a measured lay needs a massive reed"
    );
}

#[test]
fn a_lay_that_meets_the_reed_before_the_slit_shuts_changes_the_sound() {
    // At 2.8 kPa on a 6 kPa reed the static opening is about 0.21 mm, so an onset
    // of 0.3 mm keeps the reed in contact for the whole window while the default
    // lay (onset 0) is untouched: the two renders cannot agree.
    let plain = pressure(SIX);
    let curved = pressure(&format!("{SIX} 7.9e6 2 3e-4 6.47"));
    assert!(plain.iter().chain(&curved).all(|p| p.is_finite()));
    let largest = plain
        .iter()
        .zip(&curved)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        largest > 1e-6,
        "the lay left the sound unchanged: {largest}"
    );
    // The default render is deterministic, so the difference above is the lay.
    assert_eq!(pressure(SIX), plain);
}
