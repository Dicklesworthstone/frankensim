use super::*;

const CAP: usize = 256;

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 481 20000\n\
    ambient 293.15 101325 0\nreed 0.0004 0.013 6000 0 0 0.35\nlistener_m 1\n\
    termination flanged\n";
const TAIL: &str = "compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t2.8e3\n";

fn source(records: usize, body: &str) -> String {
    format!("{HEAD}segments {records}\n{body}{TAIL}")
}

fn cylinders(count: usize) -> String {
    "cylinder 0.0073 0.002\n".repeat(count)
}

#[test]
fn the_authored_and_regional_caps_agree_at_256() {
    assert_eq!(MAX_SEGMENTS, CAP);
    assert_eq!(crate::driving_point::regional::MAX_REGIONAL_SEGMENTS, CAP);
}

#[test]
fn a_source_at_the_cap_loads_and_one_more_record_refuses_at_the_count() {
    let loaded = ReedPerformance::from_bytes(source(CAP, &cylinders(CAP)).as_bytes(), 64)
        .expect("a source at the cap must load");
    assert_eq!(loaded.info().segments, CAP);
    match ReedPerformance::from_bytes(source(CAP + 1, &cylinders(CAP + 1)).as_bytes(), 64) {
        Err(ReedPerformanceError::Input { what, .. }) => {
            assert_eq!(what, "count exceeds the input budget");
        }
        other => panic!(
            "cap + 1 records must refuse at the count, got {:?}",
            other.err()
        ),
    }
}

#[test]
fn a_split_bore_with_closed_chimneys_keeps_every_authored_record() {
    // 100 axial pieces separated by 99 closed chimneys is 199 records: more than
    // the former 64-record budget and the shape a measured lattice needs.
    let mut body = String::new();
    for piece in 0..100 {
        body.push_str("cylinder 0.0073 0.005\n");
        if piece != 99 {
            body.push_str("hole 0.003 0.006 0.0073 0\n");
        }
    }
    let loaded = ReedPerformance::from_bytes(source(199, &body).as_bytes(), 64)
        .expect("199 authored records must load");
    assert_eq!(loaded.info().segments, 199);
    assert_eq!(loaded.info().tone_holes, 99);
}
