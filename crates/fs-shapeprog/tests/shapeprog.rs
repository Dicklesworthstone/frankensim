//! Battery for geometry program synthesis (fs-shapeprog). Covers DSL
//! round-trip, SDF semantics, the load-bearing rewrite-safety property (SDF
//! preserved within the declared certificate, checked by sampling),
//! canonicalization/dedup, and seeded shape-grammar derivation.

use fs_shapeprog::{
    Certificate, Geom, ParseError, linear_repeat, max_sdf_discrepancy, parse, simplify,
    stochastic_repeat,
};

fn grid() -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            for k in -3..=3 {
                pts.push([f64::from(i), f64::from(j), f64::from(k)]);
            }
        }
    }
    pts
}

#[test]
fn the_dsl_round_trips() {
    let g = Geom::sphere(2.0)
        .offset(1.0)
        .translate([1.0, 2.0, 3.0])
        .union(Geom::cube(3.0));
    assert_eq!(parse(&g.to_sexpr()), Ok(g.clone()));
    // print/parse/print is stable.
    assert_eq!(parse(&g.to_sexpr()).unwrap().to_sexpr(), g.to_sexpr());
}

#[test]
fn sdf_semantics_are_correct() {
    let s = Geom::sphere(2.0);
    assert!((s.sdf([0.0, 0.0, 0.0]) + 2.0).abs() < 1e-12); // centre: -radius
    assert!((s.sdf([3.0, 0.0, 0.0]) - 1.0).abs() < 1e-12); // 3-2 = 1
    // union is the min of the two SDFs.
    let u = Geom::sphere(1.0).union(Geom::sphere(1.0).translate([5.0, 0.0, 0.0]));
    assert!((u.sdf([5.0, 0.0, 0.0]) + 1.0).abs() < 1e-12);
}

#[test]
fn exact_rewrites_preserve_the_sdf_exactly() {
    // offset(offset(a, 1), 2) == offset(a, 3).
    let pre = Geom::sphere(2.0).offset(1.0).offset(2.0);
    let out = simplify(&pre, 1e-6);
    assert_eq!(out.program, Geom::sphere(2.0).offset(3.0));
    assert!(out.rewrites.iter().any(|r| r.rule == "offset-compose"));
    assert!((out.max_error - 0.0).abs() < 1e-15);
    // the SAFETY property: the SDF is preserved everywhere.
    assert!(max_sdf_discrepancy(&pre, &out.program, &grid()) < 1e-9);
}

#[test]
fn identity_and_distribution_rewrites_hold() {
    // union(a, empty) -> a.
    let pre = Geom::sphere(2.0).union(Geom::Empty);
    let out = simplify(&pre, 1e-6);
    assert_eq!(out.program, Geom::sphere(2.0));
    assert!(max_sdf_discrepancy(&pre, &out.program, &grid()) < 1e-9);
    // translate(union(a,b), t) -> union(translate a, translate b).
    let pre2 = Geom::sphere(1.0)
        .union(Geom::cube(1.0))
        .translate([2.0, 0.0, 0.0]);
    let out2 = simplify(&pre2, 1e-6);
    assert!(
        out2.rewrites
            .iter()
            .any(|r| r.rule == "translate-distributes")
    );
    assert!(max_sdf_discrepancy(&pre2, &out2.program, &grid()) < 1e-9);
}

#[test]
fn a_certified_approximate_rewrite_stays_within_its_bound() {
    // an offset below tolerance is dropped, certified within |r|.
    let pre = Geom::sphere(2.0).offset(0.001);
    let out = simplify(&pre, 0.01);
    assert_eq!(out.program, Geom::sphere(2.0));
    let bound = out
        .rewrites
        .iter()
        .find(|r| r.rule == "drop-tiny-offset")
        .map(|r| r.certificate);
    assert!(matches!(bound, Some(Certificate::Approximate { .. })));
    assert!((out.max_error - 0.001).abs() < 1e-12);
    // and the ACTUAL SDF discrepancy respects the certified bound.
    assert!(max_sdf_discrepancy(&pre, &out.program, &grid()) <= out.max_error + 1e-12);
}

#[test]
fn canonicalization_deduplicates_commutative_programs() {
    let a = Geom::sphere(1.0);
    let b = Geom::cube(2.0);
    let ab = a.clone().union(b.clone());
    let ba = b.union(a);
    // union is commutative -> the same content hash (archive/ledger dedup).
    assert_eq!(ab.canonical_hash(), ba.canonical_hash());
    assert_eq!(ab.canonical(), ba.canonical());
}

#[test]
fn grammar_derivations_are_reproducible_from_seeds() {
    let unit = Geom::sphere(0.4);
    let three = linear_repeat(&unit, 3, [1.0, 0.0, 0.0]);
    // three unit spheres -> each centre is inside.
    for i in 0..3 {
        assert!(three.sdf([f64::from(i), 0.0, 0.0]) < 0.0);
    }
    // the same seed yields the same derivation; different seeds may differ.
    let s1 = stochastic_repeat(&unit, 5, [1.0, 0.0, 0.0], 42);
    let s1b = stochastic_repeat(&unit, 5, [1.0, 0.0, 0.0], 42);
    assert_eq!(s1, s1b);
    assert_eq!(s1.canonical_hash(), s1b.canonical_hash());
}

#[test]
fn malformed_programs_are_rejected() {
    assert!(matches!(parse("("), Err(ParseError::UnexpectedEnd)));
    assert!(matches!(
        parse("(banana 3)"),
        Err(ParseError::Unexpected(_))
    ));
    assert!(matches!(parse("(sphere x)"), Err(ParseError::BadNumber(_))));
    // a truncated primitive reads no number -> unexpected end.
    assert!(matches!(parse("(sphere"), Err(ParseError::UnexpectedEnd)));
}
