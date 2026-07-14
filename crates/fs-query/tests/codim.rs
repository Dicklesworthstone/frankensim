//! Codimensional thickness battery (bead rjnd, part 5).
//!
//! - gd-001 G0: effective-gap brackets contain the analytic values and
//!   the three verdicts fire on their exact sides.
//! - gd-002 G0: composition with a real convex separation matches the
//!   analytic offset-body geometry (shelled spheres), and hull-only
//!   evidence downgrades contact claims instead of asserting them.
//! - gd-003 G3: thickening is monotone — growing either radius never
//!   turns a contact verdict into a clear one, and the bracket only
//!   moves down.
//! - gd-004 G0: thickness/distance refusals fail closed typed.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_query::{
    CodimThickness, CodimVerdict, ConvexSphere, QueryError, codim_gap, codim_gap_from_separation,
    convex_separation,
};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-query/codim\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC0D1,
                kernel_id: 18,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn t(radius: f64) -> CodimThickness {
    CodimThickness::new(radius).expect("valid thickness")
}

#[test]
fn gd_001_effective_gap_brackets_and_verdicts() {
    // Clear: distance 1.0, combined thickness 0.7 → effective 0.3.
    let clear = codim_gap(1.0, 1.0 + 1e-12, t(0.3), t(0.4)).expect("clear bracket");
    assert!(
        clear.lo <= 0.3 && 0.3 <= clear.hi,
        "clear bracket [{}, {}] must contain 0.3",
        clear.lo,
        clear.hi
    );
    assert_eq!(clear.verdict, CodimVerdict::ProvenClear);

    // Contact: distance at most 0.12, combined thickness 0.7.
    let contact = codim_gap(0.1, 0.12, t(0.3), t(0.4)).expect("contact bracket");
    assert!(
        contact.hi < 0.0,
        "contact bracket [{}, {}] proves interpenetration",
        contact.lo,
        contact.hi
    );
    assert_eq!(contact.verdict, CodimVerdict::ProvenContact);

    // Straddling: distance [0.5, 0.9] vs 0.7 → no claim.
    let unresolved = codim_gap(0.5, 0.9, t(0.3), t(0.4)).expect("straddling bracket");
    assert_eq!(unresolved.verdict, CodimVerdict::Unresolved);

    // Zero thickness: the effective gap IS the midsurface distance.
    let bare = codim_gap(0.5, 0.9, t(0.0), t(0.0)).expect("bare bracket");
    assert!(bare.lo <= 0.5 && bare.hi >= 0.9);
    assert_eq!(bare.verdict, CodimVerdict::ProvenClear);
    verdict(
        "gd-001",
        true,
        &format!(
            "clear [{:.3}, {:.3}], contact [{:.3}, {:.3}], straddle unresolved",
            clear.lo, clear.hi, contact.lo, contact.hi
        ),
    );
}

#[test]
fn gd_002_composition_matches_offset_body_geometry() {
    // Sphere midsurfaces (radius 0.5) with centers 4 apart: midsurface
    // distance is 3.0 exactly.
    let a = ConvexSphere::new(Point3::new(-2.0, 0.0, 0.0), 0.5).expect("a");
    let b = ConvexSphere::new(Point3::new(2.0, 0.0, 0.0), 0.5).expect("b");
    let sep = with_cx(|cx| convex_separation(&a, &b, 256, cx)).expect("separation");

    // Shells 0.6 thick each: offset radii 1.1, gap 4 - 2.2 = 1.8 > 0.
    let clear = codim_gap_from_separation(&sep, t(0.6), t(0.6), true).expect("clear");
    assert_eq!(clear.verdict, CodimVerdict::ProvenClear);
    assert!(
        clear.lo <= 1.8 && 1.8 <= clear.hi,
        "shelled clear bracket [{}, {}] must contain 1.8",
        clear.lo,
        clear.hi
    );

    // Shells 1.6 thick each: offset radii 2.1, sum 4.2 > 4 → contact.
    let contact = codim_gap_from_separation(&sep, t(1.6), t(1.6), true).expect("contact");
    assert_eq!(contact.verdict, CodimVerdict::ProvenContact);
    assert!(
        contact.lo <= -0.2 && -0.2 <= contact.hi,
        "shelled contact bracket [{}, {}] must contain -0.2",
        contact.lo,
        contact.hi
    );

    // The same numbers with hull-only evidence must NOT claim contact.
    let hull_only = codim_gap_from_separation(&sep, t(1.6), t(1.6), false).expect("hull only");
    assert_eq!(
        hull_only.verdict,
        CodimVerdict::Unresolved,
        "hull-only witnesses cannot prove contact"
    );
    // ...but clear verdicts survive hull-only evidence (lower bounds
    // only shrink under hulls).
    let hull_clear = codim_gap_from_separation(&sep, t(0.6), t(0.6), false).expect("hull clear");
    assert_eq!(hull_clear.verdict, CodimVerdict::ProvenClear);
    verdict(
        "gd-002",
        true,
        &format!(
            "shelled spheres: clear [{:.3}, {:.3}] ∋ 1.8, contact [{:.3}, {:.3}] ∋ -0.2, \
             hull-only downgraded",
            clear.lo, clear.hi, contact.lo, contact.hi
        ),
    );
}

#[test]
fn gd_003_thickening_is_monotone() {
    let radii = [0.0, 0.1, 0.25, 0.4, 0.7, 1.3];
    let mut previous_hi = f64::INFINITY;
    let mut seen_contact = false;
    for r in radii {
        let gap = codim_gap(1.0, 1.001, t(r), t(0.2)).expect("bracket");
        assert!(
            gap.hi <= previous_hi,
            "thickening must move the bracket down: hi {} after {}",
            gap.hi,
            previous_hi
        );
        if seen_contact {
            assert_eq!(
                gap.verdict,
                CodimVerdict::ProvenContact,
                "a proven contact cannot revert to clear under thickening"
            );
        }
        seen_contact |= gap.verdict == CodimVerdict::ProvenContact;
        previous_hi = gap.hi;
    }
    assert!(seen_contact, "the largest radius crosses into contact");
    verdict(
        "gd-003",
        true,
        "brackets move monotonically down; contact never reverts",
    );
}

#[test]
fn gd_004_refusals_fail_closed() {
    assert!(matches!(
        CodimThickness::new(-0.1),
        Err(QueryError::CodimInvalidThickness { .. })
    ));
    assert!(matches!(
        CodimThickness::new(f64::NAN),
        Err(QueryError::CodimInvalidThickness { .. })
    ));
    assert!(matches!(
        CodimThickness::new(f64::INFINITY),
        Err(QueryError::CodimInvalidThickness { .. })
    ));
    for (lo, hi) in [(f64::NAN, 1.0), (1.0, f64::NAN), (2.0, 1.0), (-0.5, 1.0)] {
        assert!(
            matches!(
                codim_gap(lo, hi, t(0.1), t(0.1)),
                Err(QueryError::CodimInvalidDistance { .. })
            ),
            "distance enclosure ({lo}, {hi}) must refuse"
        );
    }
    verdict(
        "gd-004",
        true,
        "thickness and distance-enclosure refusals all typed",
    );
}
