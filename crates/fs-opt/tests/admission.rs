//! fs-opt sealed-IR admission battery (bead frankensim-sj31i.48).
//!
//! G0 leaf-policy tables (manifolds, finiteness, weights, tags, caps,
//! checked dimension arithmetic), builder-rollback and
//! builder/admission agreement laws, checked ID accessors, binding
//! validation, domain-separated identity (semantic vs wire vs
//! quarantined legacy), deterministic admission-report ordering, and
//! the v3 wire round trip with legacy bilevel quarantine. Direct
//! `Problem` construction/mutation is unavailable by privacy (sealed
//! fields); a trybuild-style compile-fail harness is tracked follow-up
//! work — the seal itself is enforced by the compiler.

use fs_opt::{
    AdmissionCaps, AdmissionViolation, BilevelRef, Manifold, NodeId, OptError, ProblemBuilder,
    ProblemSemanticId, ProblemTag, Sense, VarId, WireVersion, canonical_v2_migration_target, eval,
    parse, parse_with_version, problem_hash, serialize, serialize_with_id,
};
use fs_qty::Dims;

const METER: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn simple_problem(constant: f64) -> fs_opt::Problem {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("var");
    let r = b.var_ref(v).expect("ref");
    let n = b.norm_sq(r).expect("norm_sq");
    let c = b.konst(constant, Dims::NONE).expect("konst");
    let obj = b.add(n, c).expect("add");
    b.objective(obj, Sense::Minimize, 1.0).expect("objective");
    b.finish()
}

/// adm-001 — manifold leaf policies: Rn zero-dim, Sphere ambient
/// boundaries, Stiefel frame bounds, and CHECKED point/tangent
/// formulas (no u32 wrap can reach a Problem).
#[test]
fn adm_001_manifold_leaf_policies() {
    let mut b = ProblemBuilder::new();
    assert!(matches!(
        b.var("z", Manifold::Rn { dim: 0 }, Dims::NONE),
        Err(OptError::ManifoldInvalid { .. })
    ));
    for ambient in [0u32, 1] {
        assert!(matches!(
            b.var("s", Manifold::Sphere { ambient }, Dims::NONE),
            Err(OptError::ManifoldInvalid { .. })
        ));
    }
    b.var("s2", Manifold::Sphere { ambient: 2 }, Dims::NONE)
        .expect("the smallest admitted sphere");
    // A representable but enormous sphere trips the per-variable cap,
    // not an overflow wrap.
    assert!(matches!(
        b.var("sx", Manifold::Sphere { ambient: u32::MAX }, Dims::NONE),
        Err(OptError::CapExceeded {
            what: "variable point storage",
            ..
        })
    ));
    assert!(matches!(
        b.var("st0", Manifold::Stiefel { n: 4, p: 0 }, Dims::NONE),
        Err(OptError::ManifoldInvalid { .. })
    ));
    assert!(matches!(
        b.var("stp", Manifold::Stiefel { n: 2, p: 3 }, Dims::NONE),
        Err(OptError::ManifoldInvalid { .. })
    ));
    // n * p overflows u32: the CHECKED formula refuses instead of
    // wrapping to a small bogus storage length.
    let huge = Manifold::Stiefel {
        n: 1 << 20,
        p: 1 << 20,
    };
    assert_eq!(huge.point_dim(), None, "checked formula reports overflow");
    assert!(matches!(
        b.var("sth", huge, Dims::NONE),
        Err(OptError::ManifoldInvalid { .. })
    ));
    b.var("st", Manifold::Stiefel { n: 4, p: 2 }, Dims::NONE)
        .expect("Stiefel(4,2) admits");
    assert_eq!(Manifold::Stiefel { n: 4, p: 2 }.point_dim(), Some(8));
    assert_eq!(Manifold::Stiefel { n: 4, p: 2 }.tangent_dim(), Some(5));
}

/// adm-002 — finite-payload policies: non-finite constants refuse with
/// the exact bit pattern retained; objective weights are finite
/// nonnegative with `-0.0` refused (one meaning, one wire identity).
#[test]
fn adm_002_nonfinite_and_weight_policies() {
    let mut b = ProblemBuilder::new();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let err = b.konst(bad, Dims::NONE).expect_err("non-finite constant");
        match err {
            OptError::NonFinite { what, bits } => {
                assert_eq!(what, "constant value");
                assert_eq!(bits, bad.to_bits(), "bit pattern is retained");
            }
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }
    let c = b.konst(1.0, Dims::NONE).expect("finite konst");
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            b.objective(c, Sense::Minimize, bad),
            Err(OptError::NonFinite { .. })
        ));
    }
    assert!(matches!(
        b.objective(c, Sense::Minimize, -1.0),
        Err(OptError::BadParam { .. })
    ));
    assert!(
        matches!(
            b.objective(c, Sense::Minimize, -0.0),
            Err(OptError::BadParam { .. })
        ),
        "-0.0 would give one meaning two bit-pattern wire identities"
    );
    b.objective(c, Sense::Minimize, 0.0).expect("+0.0 admits");
    b.objective(c, Sense::Maximize, 2.5)
        .expect("positive admits");
}

/// adm-003 — tag policies: nonzero capped fidelity levels, finite open
/// (0, 1) chance probabilities, and typed bilevel references.
#[test]
fn adm_003_tag_policies() {
    let mut b = ProblemBuilder::new();
    assert!(matches!(
        b.tag(ProblemTag::MultiFidelity { levels: 0 }),
        Err(OptError::BadParam { .. })
    ));
    let cap = AdmissionCaps::default().max_fidelity_levels;
    assert!(matches!(
        b.tag(ProblemTag::MultiFidelity { levels: cap + 1 }),
        Err(OptError::CapExceeded { .. })
    ));
    b.tag(ProblemTag::MultiFidelity { levels: 1 })
        .expect("one level");
    for bad in [0.0, 1.0, -0.5, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                b.tag(ProblemTag::ChanceConstrained { prob: bad }),
                Err(OptError::BadParam { .. })
            ),
            "chance prob {bad} must refuse"
        );
    }
    b.tag(ProblemTag::ChanceConstrained { prob: 0.5 })
        .expect("open interval");
    let semantic = ProblemSemanticId::from_hex(
        "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
    )
    .expect("64 hex chars");
    b.tag(ProblemTag::Bilevel {
        inner: BilevelRef::Semantic(semantic),
    })
    .expect("semantic bilevel reference");
}

/// adm-004 — rejection is rollback-free: a build interleaved with
/// refused operations finishes IDENTICAL (ids, bytes, hashes) to the
/// same build without them.
#[test]
fn adm_004_rejection_leaves_builder_unchanged() {
    let build = |with_failures: bool| {
        let mut b = ProblemBuilder::new();
        let v = b.var("x", Manifold::Rn { dim: 3 }, METER).expect("var");
        let r = b.var_ref(v).expect("ref");
        let n = b.norm_sq(r).expect("norm_sq");
        if with_failures {
            let _ = b.var("bad", Manifold::Rn { dim: 0 }, Dims::NONE);
            let _ = b.konst(f64::NAN, Dims::NONE);
            let meter = b.konst(1.0, METER).expect("legal konst is shared");
            let _ = b.add(n, meter); // m² + m refuses
            let _ = b.ln(meter); // dimensioned ln refuses
            let _ = b.tag(ProblemTag::MultiFidelity { levels: 0 });
            let _ = b.objective(n, Sense::Minimize, -1.0);
            let _ = b.component(r, 99);
        } else {
            let _ = b.konst(1.0, METER).expect("legal konst is shared");
        }
        let m2 = b.konst(4.0, Dims([2, 0, 0, 0, 0, 0])).expect("konst");
        let excess = b.sub(n, m2).expect("sub");
        b.objective(excess, Sense::Minimize, 1.0)
            .expect("objective");
        b.finish()
    };
    let clean = build(false);
    let stormy = build(true);
    assert_eq!(clean, stormy, "rejections must not disturb ids or state");
    assert_eq!(problem_hash(&clean), problem_hash(&stormy));
    assert_eq!(serialize(&clean), serialize(&stormy));
}

/// adm-005 — builder/admission agreement + identity replay: every
/// builder-finished problem re-admits; identical builds mint identical
/// semantic ids; any semantic edit (constant bits, weight, budget,
/// name, sense, tag) changes the id.
#[test]
fn adm_005_builder_admission_agreement_and_identity() {
    let p = simple_problem(2.0);
    let a1 = p.admit().expect("builder output admits");
    let a2 = simple_problem(2.0).admit().expect("identical build admits");
    assert_eq!(a1.semantic_id(), a2.semantic_id(), "replay-stable identity");
    assert_eq!(a1.schema_version(), fs_opt::ADMISSION_SCHEMA_VERSION);
    assert_eq!(a1.var_count(), 1);
    assert!(a1.quarantined_legacy_identities().is_empty());

    let base = simple_problem(2.0).admit().expect("admit").semantic_id();
    let mutants: Vec<fs_opt::Problem> = vec![
        simple_problem(2.0000000000000004), // one ulp: bit-pattern sensitivity
        {
            let mut b = ProblemBuilder::new();
            let v = b
                .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
                .expect("var");
            let r = b.var_ref(v).expect("ref");
            let n = b.norm_sq(r).expect("norm_sq");
            let c = b.konst(2.0, Dims::NONE).expect("konst");
            let obj = b.add(n, c).expect("add");
            b.objective(obj, Sense::Minimize, 2.0).expect("weight edit");
            b.finish()
        },
        {
            let mut b = ProblemBuilder::new();
            let v = b
                .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
                .expect("var");
            let r = b.var_ref(v).expect("ref");
            let n = b.norm_sq(r).expect("norm_sq");
            let c = b.konst(2.0, Dims::NONE).expect("konst");
            let obj = b.add(n, c).expect("add");
            b.objective(obj, Sense::Maximize, 1.0).expect("sense edit");
            b.finish()
        },
        {
            let mut b = ProblemBuilder::new();
            let v = b
                .var("y", Manifold::Rn { dim: 2 }, Dims::NONE)
                .expect("var");
            let r = b.var_ref(v).expect("ref");
            let n = b.norm_sq(r).expect("norm_sq");
            let c = b.konst(2.0, Dims::NONE).expect("konst");
            let obj = b.add(n, c).expect("add");
            b.objective(obj, Sense::Minimize, 1.0).expect("objective");
            b.finish()
        },
        {
            let mut b = ProblemBuilder::new();
            let v = b
                .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
                .expect("var");
            let r = b.var_ref(v).expect("ref");
            let n = b.norm_sq(r).expect("norm_sq");
            let c = b.konst(2.0, Dims::NONE).expect("konst");
            let obj = b.add(n, c).expect("add");
            b.objective(obj, Sense::Minimize, 1.0).expect("objective");
            b.set_budget(7);
            b.finish()
        },
        {
            let mut b = ProblemBuilder::new();
            let v = b
                .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
                .expect("var");
            let r = b.var_ref(v).expect("ref");
            let n = b.norm_sq(r).expect("norm_sq");
            let c = b.konst(2.0, Dims::NONE).expect("konst");
            let obj = b.add(n, c).expect("add");
            b.objective(obj, Sense::Minimize, 1.0).expect("objective");
            b.tag(ProblemTag::MultiFidelity { levels: 2 }).expect("tag");
            b.finish()
        },
    ];
    for (i, m) in mutants.iter().enumerate() {
        let id = m.admit().expect("mutant admits").semantic_id();
        assert_ne!(id, base, "semantic field {i} must be mutation-sensitive");
    }
}

/// adm-006 — admission under explicit caps is a COMPLETE report in
/// deterministic order, and re-running it is bit-identical.
#[test]
fn adm_006_admission_report_is_complete_and_deterministic() {
    let mut b = ProblemBuilder::new();
    for i in 0..3 {
        b.var(&format!("v{i}"), Manifold::Rn { dim: 2 }, Dims::NONE)
            .expect("var");
    }
    let r0 = b.var_ref(VarId(0)).expect("ref");
    let n = b.norm_sq(r0).expect("norm_sq");
    b.objective(n, Sense::Minimize, 1.0).expect("objective");
    b.constraint(n, fs_opt::ConstraintKind::LeZero, "c")
        .expect("constraint");
    let p = b.finish();

    let mut caps = AdmissionCaps::default();
    caps.max_vars = 1;
    caps.max_constraints = 0;
    let report1 = p.admit_with_caps(&caps).expect_err("caps must refuse");
    let report2 = p
        .admit_with_caps(&caps)
        .expect_err("caps must refuse again");
    assert_eq!(report1, report2, "deterministic report replay");
    assert!(
        report1.violations().len() >= 2,
        "COMPLETE report, not first-error"
    );
    assert!(matches!(
        report1.violations()[0],
        AdmissionViolation::Aggregate {
            what: "variables",
            ..
        }
    ));
    assert!(
        report1.violations().iter().any(|v| matches!(
            v,
            AdmissionViolation::Aggregate {
                what: "constraints",
                ..
            }
        )),
        "constraint cap violation is also reported"
    );
    let text = report1.to_string();
    assert!(text.contains("schema v1") && text.contains("variables"));

    p.admit().expect("default caps admit the same problem");
}

/// adm-007 — identity domain separation: semantic, wire, and legacy
/// identities are distinct types with distinct values; the wire id
/// exists only through serialization/parsing and round-trips exactly.
#[test]
fn adm_007_identity_domain_separation() {
    let p = simple_problem(3.0);
    let semantic = p.admit().expect("admit").semantic_id();
    let (text, wire) = serialize_with_id(&p);
    assert_ne!(
        semantic.as_hash(),
        wire.as_hash(),
        "domain separation: same problem, different id preimages"
    );
    let reparsed = parse_with_version(&text).expect("canonical v3 parses");
    assert_eq!(reparsed.wire_content_id(), wire);
    assert_eq!(reparsed.source_version(), WireVersion::V3);
    assert_eq!(
        reparsed
            .problem()
            .admit()
            .expect("reparse admits")
            .semantic_id(),
        semantic,
        "wire round trip preserves semantic identity"
    );
    // The legacy hash is Display/UpperHex 16-hex and stays quarantined.
    let legacy = problem_hash(&p);
    assert_eq!(format!("{legacy}").len(), 16);
    assert_eq!(format!("{legacy:016X}"), format!("{:016X}", legacy.get()));
}

/// adm-008 — v3 bilevel identities: semantic references round-trip
/// full-width; historical v2 artifacts parse with their FNV identities
/// QUARANTINED (typed, listed by admission, re-encoded as the explicit
/// `bilevel_legacy` v3 spelling).
#[test]
fn adm_008_v3_bilevel_roundtrip_and_legacy_quarantine() {
    let semantic = ProblemSemanticId::from_hex(
        "00112233445566778899aabbccddeeff00112233445566778899aabbccddeeff",
    )
    .expect("hex");
    let mut b = ProblemBuilder::new();
    let c = b.konst(1.0, Dims::NONE).expect("konst");
    b.objective(c, Sense::Minimize, 1.0).expect("objective");
    b.tag(ProblemTag::Bilevel {
        inner: BilevelRef::Semantic(semantic),
    })
    .expect("semantic tag");
    let p = b.finish();
    let text = serialize(&p);
    assert!(text.starts_with("fsopt v3\n"));
    assert!(text.contains(&format!("tag bilevel {}", semantic.to_hex())));
    let back = parse(&text).expect("v3 round trip");
    assert_eq!(back, p);
    assert!(
        back.admit()
            .expect("admits")
            .quarantined_legacy_identities()
            .is_empty()
    );

    // Historical v2 spelling: quarantined on read, explicit on re-write.
    let mut lb = ProblemBuilder::new();
    let lc = lb.konst(1.0, Dims::NONE).expect("konst");
    lb.objective(lc, Sense::Minimize, 1.0).expect("objective");
    lb.tag(ProblemTag::Bilevel {
        inner: BilevelRef::LegacyFnv(fs_opt::LegacyProblemHash::new(0xDEAD_BEEF_0123_4567)),
    })
    .expect("legacy tag");
    let lp = lb.finish();
    let v2 = canonical_v2_migration_target(&lp);
    assert!(v2.contains("tag bilevel DEADBEEF01234567"));
    let decoded = parse_with_version(&v2).expect("exact v2 bytes decode");
    assert_eq!(decoded.source_version(), WireVersion::V2);
    let admission = decoded.problem().admit().expect("admits with quarantine");
    assert_eq!(
        admission.quarantined_legacy_identities(),
        &[fs_opt::LegacyProblemHash::new(0xDEAD_BEEF_0123_4567)],
        "the legacy identity is listed, not laundered"
    );
    assert!(
        serialize(decoded.problem()).contains("tag bilevel_legacy DEADBEEF01234567"),
        "v3 re-encoding keeps the quarantine explicit"
    );
}

/// adm-009 — checked ID accessors: unknown node/var ids refuse with
/// typed errors everywhere a panic used to be reachable.
#[test]
fn adm_009_checked_id_accessors() {
    let p = simple_problem(1.0);
    let bogus = NodeId(9999);
    assert!(matches!(
        p.shape(bogus),
        Err(OptError::UnknownNode { id: 9999 })
    ));
    assert!(matches!(
        p.node_dims(bogus),
        Err(OptError::UnknownNode { .. })
    ));
    assert!(matches!(p.class(bogus), Err(OptError::UnknownNode { .. })));
    assert!(matches!(p.expr(bogus), Err(OptError::UnknownNode { .. })));
    assert!(matches!(
        p.reachable(bogus),
        Err(OptError::UnknownNode { .. })
    ));
    assert!(matches!(
        p.variable(VarId(77)),
        Err(OptError::UnknownVar { id: 77 })
    ));
    let root = p.objectives()[0].node;
    assert!(p.reachable(root).expect("valid root").contains(&root));
    assert!(matches!(
        eval(&p, bogus, &[vec![0.0, 0.0]]),
        Err(OptError::UnknownNode { .. })
    ));
}

/// adm-010 — binding validation: surplus bindings and wrong-length
/// bindings refuse with typed errors; a missing prefix binding teaches
/// through `UnknownVar` when actually referenced.
#[test]
fn adm_010_binding_validation() {
    let p = simple_problem(1.0);
    let root = p.objectives()[0].node;
    assert!(matches!(
        eval(&p, root, &[vec![0.0, 0.0], vec![1.0]]),
        Err(OptError::BindingCount { vars: 1, got: 2 })
    ));
    assert!(matches!(
        eval(&p, root, &[vec![0.0, 0.0, 0.0]]),
        Err(OptError::BindingLen {
            var: 0,
            expected: 2,
            got: 3
        })
    ));
    assert!(matches!(
        eval(&p, root, &[]),
        Err(OptError::UnknownVar { id: 0 })
    ));
    eval(&p, root, &[vec![0.5, -0.5]]).expect("exact binding evaluates");
}

/// adm-011 — caps bind incrementally: node/var caps refuse BEFORE
/// mutation, while hash-consed duplicates still return their existing
/// id under a full arena.
#[test]
fn adm_011_caps_bind_incrementally() {
    let mut caps = AdmissionCaps::default();
    caps.max_nodes = 2;
    caps.max_vars = 1;
    caps.max_name_bytes = 8;
    let mut b = ProblemBuilder::with_caps(caps);
    assert!(matches!(
        b.var("waaaaay-too-long-name", Manifold::Rn { dim: 1 }, Dims::NONE),
        Err(OptError::CapExceeded {
            what: "variable name",
            ..
        })
    ));
    let v = b
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("var");
    assert!(matches!(
        b.var("y", Manifold::Rn { dim: 1 }, Dims::NONE),
        Err(OptError::CapExceeded {
            what: "variables",
            ..
        })
    ));
    let r = b.var_ref(v).expect("node 1 of 2");
    let n = b.norm_sq(r).expect("node 2 of 2");
    assert_eq!(
        b.var_ref(v).expect("interned duplicate"),
        r,
        "CSE under a full arena"
    );
    assert!(matches!(
        b.abs(n),
        Err(OptError::CapExceeded {
            what: "expression nodes",
            ..
        })
    ));
}

/// adm-012 — dimension arithmetic is CHECKED end to end: exponent
/// combinations that used to saturate silently now refuse.
#[test]
fn adm_012_checked_dimension_arithmetic() {
    let mut b = ProblemBuilder::new();
    let big = b.konst(1.0, Dims([100, 0, 0, 0, 0, 0])).expect("konst");
    let err = b.mul(big, big).expect_err("m^200 cannot fit i8 dims");
    assert!(matches!(err, OptError::DimSumOverflow { op: "mul", .. }));
    let neg = b.konst(1.0, Dims([-100, 0, 0, 0, 0, 0])).expect("konst");
    assert!(matches!(
        b.div(neg, big),
        Err(OptError::DimSumOverflow { op: "div", .. })
    ));
    let fine = b.mul(big, neg).expect("m^0 is representable");
    assert_eq!(b.finish().node_dims(fine).expect("known node"), Dims::NONE);
}

/// adm-013 — arbitrary component indices cannot panic: asking for
/// component u32::MAX of a scalar is a REPORTED shape mismatch with
/// saturated diagnostic arithmetic, never a wrap or debug panic.
#[test]
fn adm_013_max_component_index_reports_instead_of_panicking() {
    let mut b = ProblemBuilder::new();
    let scalar = b.konst(1.0, Dims::NONE).expect("konst");
    let err = b
        .component(scalar, u32::MAX)
        .expect_err("a scalar has no components");
    match err {
        OptError::ShapeMismatch { op, right, .. } => {
            assert_eq!(op, "component");
            assert_eq!(
                format!("{right:?}"),
                format!("{:?}", fs_opt::Shape::Vector(u32::MAX)),
                "the required-shape diagnostic saturates at the u32 boundary"
            );
        }
        other => panic!("expected a shape mismatch, got {other:?}"),
    }
}

/// adm-014 — Min/Max are scalar-only END TO END: the shared leaf rule
/// refuses vector operands (both evaluators implement only the scalar
/// case), so no admitted program can reach their panic arms.
#[test]
fn adm_014_vector_min_max_refused_at_the_shared_rule() {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 3 }, Dims::NONE)
        .expect("var");
    let r = b.var_ref(v).expect("ref");
    let err = b.min_of(r, r).expect_err("vector min must refuse");
    assert!(
        matches!(
            err,
            OptError::ShapeMismatch {
                op: "min",
                right: fs_opt::Shape::Scalar,
                ..
            }
        ),
        "vector min refuses with the scalar requirement, got {err:?}"
    );
    let err = b.max_of(r, r).expect_err("vector max must refuse");
    assert!(matches!(err, OptError::ShapeMismatch { op: "max", .. }));
    // The scalar case is untouched.
    let s1 = b.konst(1.0, Dims::NONE).expect("konst");
    let s2 = b.konst(2.0, Dims::NONE).expect("konst");
    let m = b.min_of(s1, s2).expect("scalar min admits");
    b.objective(m, Sense::Minimize, 1.0).expect("objective");
    b.finish().admit().expect("re-admits");
}
