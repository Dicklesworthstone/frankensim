//! Nominal-versus-as-built projection tests (bead
//! `frankensim-extreal-program-f85xj.12.2` item (d)).
//!
//! G0 for the row algebra (shift, half-width extraction, resolution tri-state)
//! and G3-flavoured checks on the rendered projection: the correlation sentence
//! must follow the terms actually supplied, not the story the caller wants.
//!
//! What is deliberately NOT tested here: that the upstream propagation is
//! correct. This crate is presentation only; `fs-asbuilt`'s own suites own the
//! covariance mathematics.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::uncertainty::{
    DistributionTerm, EngineeringUncertaintyKind, EngineeringUncertaintyTerm, TermValue,
    UncertaintyArtifactRef,
};
use fs_report::{AsBuiltDeltaError, AsBuiltQoiDelta, nominal_vs_as_built_markdown};

const ROLE: &str = "as-built-pose-propagation";

fn digest(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-report.test.as-built.v1",
        label.as_bytes(),
    )
}

fn artifact(label: &str) -> UncertaintyArtifactRef {
    UncertaintyArtifactRef::new(ROLE, digest(label)).expect("valid artifact fixture")
}

/// A geometry term standing in for one produced by
/// `fs_asbuilt::propagate::GeometryPropagation::geometry_term`.
fn geometry_term(record: &str, half_width: f64) -> EngineeringUncertaintyTerm {
    let replay = artifact(record);
    let value = TermValue::Distribution(DistributionTerm {
        mean: 0.0,
        standard_deviation: half_width / 2.0,
        conservative_half_width: half_width,
        level: 0.95,
        replay: replay.clone(),
    });
    EngineeringUncertaintyTerm::try_new(EngineeringUncertaintyKind::Geometry, value, replay)
        .expect("valid geometry term")
}

fn unknown_geometry_term(record: &str) -> EngineeringUncertaintyTerm {
    let value = TermValue::unknown(
        "linearization rejected: max relative gap 0.42 exceeds tolerance 0.05 at qoi t_junction",
    )
    .expect("named unknown");
    EngineeringUncertaintyTerm::try_new(
        EngineeringUncertaintyKind::Geometry,
        value,
        artifact(record),
    )
    .expect("valid geometry term")
}

#[test]
fn ab_001_row_construction_is_fail_closed() {
    let geometry = geometry_term("record-a", 0.5);

    // The row must cite the geometry slot. Citing any other term would let a
    // discretization or model-form budget masquerade as pose attribution.
    let parameters = EngineeringUncertaintyTerm::try_new(
        EngineeringUncertaintyKind::Parameters,
        TermValue::negligible("exact in this fixture").expect("named negligible"),
        artifact("record-a"),
    )
    .expect("valid parameters term");
    let refusal = AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.0, &parameters)
        .expect_err("a non-geometry term must be refused");
    assert_eq!(
        refusal,
        AsBuiltDeltaError::NotGeometryTerm {
            kind: EngineeringUncertaintyKind::Parameters
        }
    );
    assert!(refusal.to_string().contains("geometry term"));

    for (nominal, as_built, field) in [
        (f64::NAN, 351.0, "nominal"),
        (350.0, f64::INFINITY, "as_built"),
    ] {
        let refusal =
            AsBuiltQoiDelta::try_new("t_junction", "kelvin", nominal, as_built, &geometry)
                .expect_err("non-finite values must be refused");
        assert_eq!(refusal, AsBuiltDeltaError::NonFiniteValue { field });
    }

    let row = AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &geometry)
        .expect("valid row");
    assert_eq!(row.qoi(), "t_junction");
    assert_eq!(row.unit(), "kelvin");
    assert_eq!(row.nominal(), 350.0);
    assert_eq!(row.as_built(), 351.25);
    assert_eq!(row.shift(), 1.25);
    assert_eq!(row.geometry().kind(), EngineeringUncertaintyKind::Geometry);
    assert!(row.describe().contains("shift +1.25 kelvin"));
}

#[test]
fn ab_002_resolution_is_a_tri_state_not_a_boolean() {
    // Shift outside the pose-propagated band: the measurement resolves it.
    let wide = geometry_term("record-a", 0.5);
    let resolved =
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &wide).expect("row");
    assert_eq!(resolved.geometry_half_width(), Some(0.5));
    assert_eq!(resolved.shift_resolved(), Some(true));

    // Shift inside the band: the two solves differ, but this measurement
    // cannot tell the difference from its own pose noise.
    let narrow_shift =
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 350.25, &wide).expect("row");
    assert_eq!(narrow_shift.shift_resolved(), Some(false));

    // Exactly at the half-width is NOT resolved: the comparison is strict, so
    // a boundary case never reads as evidence of a real shift.
    let boundary =
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 350.5, &wide).expect("row");
    assert_eq!(boundary.shift().abs(), 0.5);
    assert_eq!(boundary.shift_resolved(), Some(false));

    // A rejected linearization upstream produces an Unknown term. The question
    // then has no answer — which must not collapse to "not resolved".
    let unknown = unknown_geometry_term("record-a");
    let undecidable =
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 380.0, &unknown).expect("row");
    assert_eq!(undecidable.geometry_half_width(), None);
    assert_eq!(undecidable.shift_resolved(), None);

    // An interval-bound geometry term publishes a half-width too. This variant
    // brackets the HALF-WIDTH itself in [lower, upper] — it is not a [min, max]
    // band around the value — so the conservative reading is the upper
    // endpoint, matching what fs-evidence uses when aggregating the marginal.
    // Reading it as a midpoint would report less uncertainty than the budget.
    let interval = EngineeringUncertaintyTerm::try_new(
        EngineeringUncertaintyKind::Geometry,
        TermValue::interval(0.25, 0.5).expect("ordered non-negative half-width bracket"),
        artifact("record-a"),
    )
    .expect("valid interval term");
    let banded =
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 350.375, &interval).expect("row");
    assert_eq!(banded.geometry_half_width(), Some(0.5));
    // 0.375 is inside the conservative 0.5 half-width but outside the optimistic
    // 0.25 one: this case fails if the upper endpoint is ever swapped for a
    // midpoint or a lower bound.
    assert_eq!(banded.shift(), 0.375);
    assert_eq!(banded.shift_resolved(), Some(false));
}

#[test]
fn ab_003_shared_record_renders_as_one_correlated_block() {
    // Both QoIs cite the SAME propagation record — the upstream contract for
    // one pose covariance moving everything at once.
    let junction = geometry_term("record-a", 0.5);
    let case = geometry_term("record-a", 0.2);
    assert_eq!(
        junction.provenance().digest(),
        case.provenance().digest(),
        "fixture must model one shared record"
    );

    let rows = vec![
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &junction).expect("row"),
        AsBuiltQoiDelta::try_new("t_case", "kelvin", 320.0, 320.125, &case).expect("row"),
    ];
    let rendered = nominal_vs_as_built_markdown(&rows);

    assert!(rendered.starts_with("## Nominal versus as-built\n\n"));
    assert!(
        rendered.contains("**`t_junction`:** `350 kelvin` → `351.25 kelvin`; shift `+1.25 kelvin`")
    );
    // Fixture values are exactly representable in binary on purpose: the
    // renderer prints the exact f64, so a fixture like 320.1 would assert on
    // decimal-conversion noise instead of on the projection.
    assert!(
        rendered.contains("**`t_case`:** `320 kelvin` → `320.125 kelvin`; shift `+0.125 kelvin`")
    );
    assert!(rendered.contains("resolved: the shift exceeds the pose-propagated half-width"));
    assert!(
        rendered.contains("NOT resolved: the shift sits inside the pose-propagated half-width")
    );
    assert!(rendered.contains("all 2 geometry terms cite one propagation record"));
    assert!(rendered.contains("single correlated block"));
    assert!(rendered.contains(&digest("record-a").to_hex()));
    assert!(rendered.contains(ROLE));
    assert!(rendered.contains("**No-claim boundary:**"));
    assert!(rendered.contains("not a measurement of the physical part"));

    // Presentation is deterministic: same rows, same bytes.
    assert_eq!(rendered, nominal_vs_as_built_markdown(&rows));

    // Row order is caller order, not a sort. The report must follow the study.
    let reversed = vec![
        AsBuiltQoiDelta::try_new("t_case", "kelvin", 320.0, 320.125, &case).expect("row"),
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &junction).expect("row"),
    ];
    let reversed_render = nominal_vs_as_built_markdown(&reversed);
    assert!(
        reversed_render.find("t_case").expect("present")
            < reversed_render.find("t_junction").expect("present")
    );
}

#[test]
fn ab_004_distinct_records_refuse_the_correlated_block_claim() {
    // Two QoIs attributed to DIFFERENT propagation records. Nothing here
    // establishes that one pose moved both, so the renderer must say so
    // instead of implying a block that was never computed.
    let junction = geometry_term("record-a", 0.5);
    let case = geometry_term("record-b", 0.2);
    assert_ne!(junction.provenance().digest(), case.provenance().digest());

    let rows = vec![
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &junction).expect("row"),
        AsBuiltQoiDelta::try_new("t_case", "kelvin", 320.0, 320.125, &case).expect("row"),
    ];
    let rendered = nominal_vs_as_built_markdown(&rows);
    assert!(rendered.contains("cite 2 distinct propagation records"));
    assert!(rendered.contains("cross-QoI correlation is not established here"));
    // The affirmative sentence must be absent. Assert on the affirmative
    // phrasing itself: the refusal text legitimately contains the words "single
    // correlated block" inside "NO single correlated block", so a bare
    // substring check would pass for the wrong reason.
    assert!(!rendered.contains("cite one propagation record"));
    assert!(rendered.contains("NO single correlated block"));
}

#[test]
fn ab_005_single_row_and_empty_table_make_no_correlation_claim() {
    let junction = geometry_term("record-a", 0.5);
    let rows = vec![
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 351.25, &junction).expect("row"),
    ];
    let rendered = nominal_vs_as_built_markdown(&rows);
    assert!(rendered.contains("one QoI, one propagation record"));
    assert!(rendered.contains("no cross-QoI correlation is exercised"));
    assert!(!rendered.contains("single correlated block"));

    let empty = nominal_vs_as_built_markdown(&[]);
    assert!(empty.contains("no QoI was compared"));
    assert!(empty.contains("makes no as-built claim"));
    assert!(!empty.contains("**Correlation:**"));
}

#[test]
fn ab_006_unknown_geometry_term_refuses_to_attribute_the_shift() {
    // A 30 K shift is enormous, but the upstream linearization was rejected.
    // The renderer must not let the size of the number imply attribution.
    let unknown = unknown_geometry_term("record-a");
    let rows = vec![
        AsBuiltQoiDelta::try_new("t_junction", "kelvin", 350.0, 380.0, &unknown).expect("row"),
    ];
    let rendered = nominal_vs_as_built_markdown(&rows);
    assert!(rendered.contains("shift `+30 kelvin`"));
    assert!(rendered.contains("geometry term: `unknown`"));
    assert!(rendered.contains("linearization rejected"));
    assert!(rendered.contains("the shift is NOT attributed"));
    assert!(rendered.contains("not decidable: the geometry term publishes no half-width"));
    assert!(!rendered.contains("half-width `"));
}
