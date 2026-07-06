//! Conformance suite for fs-qty (plan §13.3): any reimplementation must pass
//! these cases. Each case logs a JSON-lines verdict so failures are
//! diagnosable from output alone (observability standard; fs-obs schema
//! adoption pending that crate's landing).

use fs_qty::parse::parse_qty;
use fs_qty::{Dims, DynViscosity, Length, Pressure, QtyAny, Time};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-qty/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

/// qty-001: the FrankenScript literal battery — exact values and dimensions
/// for every literal form the plan's example studies use (Appendix C).
#[test]
fn qty_001_appendix_c_literal_battery() {
    let cases: &[(&str, f64, [i8; 5])] = &[
        ("0.12Pa*s", 0.12, [-1, 1, -1, 0, 0]),
        ("0.061N/m", 0.061, [0, 1, -2, 0, 0]),
        ("0.5L/s", 5e-4, [3, 0, -1, 0, 0]),
        ("12mm", 0.012, [1, 0, 0, 0, 0]),
        (
            "65deg",
            65.0 * core::f64::consts::PI / 180.0,
            [0, 0, 0, 0, 0],
        ),
        ("2h", 7200.0, [0, 0, 1, 0, 0]),
        ("0.03m2/s3", 0.03, [2, 0, -3, 0, 0]),
        ("15rad/s", 15.0, [0, 0, -1, 0, 0]),
        ("2e-2", 0.02, [0, 0, 0, 0, 0]),
    ];
    for (text, value, dims) in cases {
        let q = parse_qty(text).unwrap_or_else(|e| panic!("{text}: {e}"));
        let ok = (q.value - value).abs() <= 1e-12 * value.abs().max(1.0) && q.dims == Dims(*dims);
        verdict(
            &format!("qty-001/{text}"),
            ok,
            &format!("value={} dims={:?}", q.value, q.dims),
        );
    }
}

/// qty-002: typed and erased algebra agree bit-for-bit.
#[test]
fn qty_002_typed_erased_agreement() {
    let typed = (Length::new(0.37) / Time::new(1.61)).value();
    let erased = (Length::new(0.37).erase() / Time::new(1.61).erase()).value;
    verdict(
        "qty-002/bit-agreement",
        typed.to_bits() == erased.to_bits(),
        &format!("typed={typed:?} erased={erased:?}"),
    );
}

/// qty-003: JSON round-trip is bit-exact and shape-canonical.
#[test]
fn qty_003_json_round_trip() {
    let q = DynViscosity::new(0.12).erase();
    let text = fs_qty::json::to_json(q).expect("finite");
    let canonical = text == r#"{"value":0.12,"dims":[-1,1,-1,0,0]}"#;
    let back = fs_qty::json::from_json(&text).expect("parses");
    verdict(
        "qty-003/round-trip",
        canonical && back.value.to_bits() == q.value.to_bits() && back.dims == q.dims,
        &text,
    );
}

/// qty-004: dimension safety — runtime mismatches produce structured,
/// teaching errors; downcasts check dimensions.
#[test]
fn qty_004_dimension_safety() {
    let e = Pressure::new(1.0)
        .erase()
        .try_add(Time::new(1.0).erase())
        .unwrap_err();
    verdict(
        "qty-004/mismatch-error",
        e.to_string().contains("dimension mismatch"),
        &e.to_string(),
    );
    let bad: Result<Pressure, _> = QtyAny::dimensionless(1.0).to_typed();
    verdict(
        "qty-004/downcast-checked",
        bad.is_err(),
        "dimensionless -> Pressure must fail",
    );
}

/// qty-005: the parser never panics (a compressed re-run of the hardening
/// battery at conformance level — reimplementations must hold this too).
#[test]
fn qty_005_parser_total_over_garbage() {
    let mut state: u64 = 0x00C0_FFEE;
    for _ in 0..2_000 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bytes: Vec<u8> = (0..(state % 16))
            .map(|i| ((state >> (i % 56)) & 0x7F) as u8)
            .collect();
        let s = String::from_utf8_lossy(&bytes);
        let _ = parse_qty(&s);
    }
    verdict("qty-005/no-panic", true, "2000 garbage inputs, no panic");
}
