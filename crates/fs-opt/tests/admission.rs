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
    assert!(
        text.contains(&format!("schema v{}", fs_opt::ADMISSION_SCHEMA_VERSION))
            && text.contains("variables")
    );

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
    let v2 = canonical_v2_migration_target(&lp)
        .expect("a legacy bilevel reference is representable in v2");
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

/// adm-010 — binding validation: surplus, wrong-length, and non-finite
/// bindings refuse with typed errors; a missing prefix binding teaches
/// through `UnknownVar` when actually referenced. Runtime points are
/// rejected before their components can acquire graph authority.
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
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            eval(&p, root, &[vec![0.0, bad]]),
            Err(OptError::BindingNonFinite {
                var: 0,
                component: 1,
                bits: bad.to_bits(),
            }),
            "the exact malformed binding coordinate and bits must be retained"
        );
    }
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

/// adm-015 — recursive work and retained allocation are bounded before
/// mutation. A depth overflow and an aggregate-byte overflow leave the
/// builder in an admissible state, while oversized external identifiers
/// refuse before they are cloned into an expression.
#[test]
fn adm_015_depth_and_retained_bytes_fail_closed() {
    let mut depth_caps = AdmissionCaps::default();
    depth_caps.max_graph_depth = 3;
    let mut depth = ProblemBuilder::with_caps(depth_caps.clone());
    let n1 = depth.konst(1.0, Dims::NONE).expect("depth 1");
    let n2 = depth.neg(n1).expect("depth 2");
    let n3 = depth.neg(n2).expect("depth 3");
    assert!(matches!(
        depth.neg(n3),
        Err(OptError::CapExceeded {
            what: "graph depth",
            count: 4,
            cap: 3
        })
    ));
    depth
        .objective(n3, Sense::Minimize, 1.0)
        .expect("rollback left a valid root");
    depth
        .finish()
        .admit_with_caps(&depth_caps)
        .expect("depth rejection did not mutate the graph");

    let mut retained_caps = AdmissionCaps::default();
    retained_caps.max_total_retained_bytes = 258;
    let mut retained = ProblemBuilder::with_caps(retained_caps.clone());
    let v = retained
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("name charge is bounded");
    retained.var_ref(v).expect("exact aggregate byte cap");
    assert!(matches!(
        retained.norm_sq(NodeId(0)),
        Err(OptError::CapExceeded {
            what: "total retained/canonical bytes",
            ..
        })
    ));
    let admitted = retained
        .finish()
        .admit_with_caps(&retained_caps)
        .expect("aggregate rejection did not retain the candidate");
    assert_eq!(admitted.total_retained_bytes(), 258);

    let mut string_caps = AdmissionCaps::default();
    string_caps.max_string_bytes = 4;
    let mut strings = ProblemBuilder::with_caps(string_caps);
    let v = strings
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("var");
    assert!(matches!(
        strings.pde_residual("12345", v, false, Dims::NONE),
        Err(OptError::CapExceeded {
            what: "PDE study",
            count: 5,
            cap: 4
        })
    ));
    assert_eq!(
        strings
            .pde_residual("ok", v, false, Dims::NONE)
            .expect("rejection happened before mutation"),
        NodeId(0)
    );
}

/// G4 admission envelope — aggregate work is charged incrementally and
/// cap+1 refuses before retaining a new item. The admitted receipt
/// re-derives the exact same charge.
#[test]
fn aggregate_work_cap_is_incremental_and_receipted() {
    let mut caps = AdmissionCaps::default();
    caps.max_total_work = 2;
    let mut builder = ProblemBuilder::with_caps(caps.clone());
    let variable = builder
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("work unit 1");
    builder.var_ref(variable).expect("work unit 2");
    assert!(matches!(
        builder.tag(ProblemTag::MultiFidelity { levels: 1 }),
        Err(OptError::CapExceeded {
            what: "total admission work",
            count: 3,
            cap: 2,
        })
    ));
    let problem = builder.finish();
    assert!(problem.tags().is_empty(), "cap+1 tag was not retained");
    assert_eq!(problem.total_admission_work(), 2);
    assert_eq!(
        problem
            .admit_with_caps(&caps)
            .expect("bounded prefix re-admits")
            .total_work(),
        2
    );
}

/// G4/G5 depth boundary — a chain exactly at the default depth cap
/// builds, parses, admits, and evaluates deterministically. The next
/// node refuses in the builder; a problem deliberately built under a
/// looser cap is rejected by default admission, parsing, and evaluation
/// before recursive work begins.
#[test]
fn graph_depth_boundary_is_closed_everywhere() {
    fn chain(caps: AdmissionCaps, depth: u32) -> (fs_opt::Problem, NodeId) {
        let mut builder = ProblemBuilder::with_caps(caps);
        let mut root = builder.konst(1.0, Dims::NONE).expect("depth 1");
        for _ in 1..depth {
            root = builder.neg(root).expect("next admitted depth");
        }
        builder
            .objective(root, Sense::Minimize, 1.0)
            .expect("scalar objective");
        (builder.finish(), root)
    }

    let default_caps = AdmissionCaps::default();
    let limit = default_caps.max_graph_depth;
    let (at_limit, root) = chain(default_caps.clone(), limit);
    let receipt = at_limit.admit().expect("exact depth cap admits");
    assert_eq!(receipt.graph_depth(), limit);
    let value = eval(&at_limit, root, &[])
        .expect("depth-gated evaluator accepts the boundary")
        .scalar()
        .expect("scalar chain");
    let expected: f64 = if limit % 2 == 0 { -1.0 } else { 1.0 };
    assert_eq!(value.to_bits(), expected.to_bits());
    parse(&serialize(&at_limit)).expect("wire parser accepts the exact boundary");

    let mut builder = ProblemBuilder::with_caps(default_caps.clone());
    let mut capped = builder.konst(1.0, Dims::NONE).expect("depth 1");
    for _ in 1..limit {
        capped = builder.neg(capped).expect("within depth cap");
    }
    assert!(matches!(
        builder.neg(capped),
        Err(OptError::CapExceeded {
            what: "graph depth",
            count,
            cap,
        }) if count == u64::from(limit) + 1 && cap == u64::from(limit)
    ));

    let mut relaxed = default_caps;
    relaxed.max_graph_depth = limit + 1;
    let (over_limit, over_root) = chain(relaxed, limit + 1);
    assert!(
        over_limit.admit().is_err(),
        "default admission refuses max+1"
    );
    assert!(matches!(
        eval(&over_limit, over_root, &[]),
        Err(OptError::CapExceeded {
            what: "graph depth",
            count,
            cap,
        }) if count == u64::from(limit) + 1 && cap == u64::from(limit)
    ));
    assert!(matches!(
        parse(&serialize(&over_limit)),
        Err(OptError::Parse { what, .. }) if what.contains("graph depth")
    ));
}

/// adm-013 (bead j3vb5, review High #6) — descent leaf gating: a
/// non-finite start component or degenerate step policy refuses TYPED
/// through the public API before any descent arithmetic, for both the
/// closure-driven and IR-driven entry points. (Invalid manifolds and
/// wrong-length starts are pinned by opt-005 in the conformance
/// suite.)
#[test]
fn adm_013_descent_leaf_gating() {
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x13,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let quadratic = |x: &[f64]| x[0] * x[0];
        let opts = fs_opt::DescentOptions::default();

        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let refusal =
                fs_opt::descend_fn(Manifold::Rn { dim: 1 }, &quadratic, &[bad], opts, 0, &cx)
                    .expect_err("non-finite start must refuse");
            match refusal {
                OptError::NonFinite { what, bits } => {
                    assert_eq!(what, "descent initial point component");
                    assert_eq!(bits, bad.to_bits(), "offending bit pattern retained");
                }
                other => panic!("expected NonFinite, got {other:?}"),
            }
        }

        for bad_fd in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
            let mut o = fs_opt::DescentOptions::default();
            o.fd_h = bad_fd;
            assert!(
                matches!(
                    fs_opt::descend_fn(Manifold::Rn { dim: 1 }, &quadratic, &[1.0], o, 0, &cx),
                    Err(OptError::BadParam { .. })
                ),
                "fd_h {bad_fd} must refuse before any FD division"
            );
        }
        for bad_lr in [0.0, -0.5, f64::NAN] {
            let mut o = fs_opt::DescentOptions::default();
            o.lr = bad_lr;
            assert!(
                matches!(
                    fs_opt::descend_fn(Manifold::Rn { dim: 1 }, &quadratic, &[1.0], o, 0, &cx),
                    Err(OptError::BadParam { .. })
                ),
                "lr {bad_lr} must refuse (descent, not ascent or a no-op)"
            );
        }

        // The IR entry point must share the same ordering: invalid
        // options refuse before even an unevaluable objective is asked
        // to run.
        let mut unevaluable_builder = ProblemBuilder::new();
        let design = unevaluable_builder
            .var("design", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("design variable");
        let pde = unevaluable_builder
            .pde_residual("leaf-gate-order", design, true, Dims::NONE)
            .expect("PDE node");
        unevaluable_builder
            .objective(pde, Sense::Minimize, 1.0)
            .expect("objective");
        let unevaluable = unevaluable_builder.finish();
        let mut invalid_options = opts;
        invalid_options.fd_h = 0.0;
        assert!(matches!(
            fs_opt::descend_ir(&unevaluable, &[0.0], invalid_options, &cx),
            Err(OptError::BadParam { .. })
        ));

        // The IR-driven variant shares the seam: a NaN start refuses
        // typed instead of descending on garbage.
        let p = simple_problem(1.0);
        assert!(matches!(
            fs_opt::descend_ir(&p, &[f64::NAN, 0.0], opts, &cx),
            Err(OptError::NonFinite { .. })
        ));

        // Valid descent is unchanged: the guarded path still converges.
        let rep = fs_opt::descend_fn(Manifold::Rn { dim: 1 }, &quadratic, &[1.0], opts, 0, &cx)
            .expect("valid descent unchanged");
        assert!(rep.f_final < rep.f0, "still actually descends");
    });
}

/// adm-016 — the public raw retraction boundary is total over arbitrary
/// slice lengths. Every live manifold refuses both short and long point
/// or step storage before zip truncation or direct indexing can occur.
#[test]
fn adm_016_retraction_storage_lengths_are_exact() {
    let cases = [
        (Manifold::Rn { dim: 2 }, vec![0.0, 0.0], vec![0.0, 0.0]),
        (
            Manifold::Sphere { ambient: 3 },
            vec![1.0, 0.0, 0.0],
            vec![0.0; 3],
        ),
        (Manifold::So3, vec![1.0, 0.0, 0.0, 0.0], vec![0.0; 3]),
        (
            Manifold::Stiefel { n: 3, p: 2 },
            vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            vec![0.0; 6],
        ),
    ];

    for (manifold, point, step) in cases {
        let point_dim = manifold.point_dim().expect("valid fixture manifold");
        let step_dim = manifold.param_dim().expect("valid fixture manifold");

        let mut long_point = point.clone();
        long_point.push(0.0);
        for malformed in [&point[..point.len() - 1], long_point.as_slice()] {
            assert!(matches!(
                manifold.retract(malformed, &step),
                Err(OptError::RetractionLen {
                    input: "retraction point",
                    expected,
                    got,
                }) if expected == point_dim && got == malformed.len() as u64
            ));
        }

        let mut long_step = step.clone();
        long_step.push(0.0);
        for malformed in [&step[..step.len() - 1], long_step.as_slice()] {
            assert!(matches!(
                manifold.retract(&point, malformed),
                Err(OptError::RetractionLen {
                    input: "retraction step",
                    expected,
                    got,
                }) if expected == step_dim && got == malformed.len() as u64
            ));
        }

        assert_eq!(
            manifold
                .retract(&point, &step)
                .expect("exact storage retracts")
                .len(),
            point.len()
        );
    }
}

/// adm-017 — finite runtime bindings do not authorize NaN/Inf graph
/// results. Evaluation identifies the exact node/component, and IR
/// descent propagates a probe-domain refusal instead of panicking or
/// publishing a non-finite report.
#[test]
fn adm_017_non_finite_runtime_results_refuse() {
    let mut builder = ProblemBuilder::new();
    let variable = builder
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("variable");
    let point = builder.var_ref(variable).expect("point");
    let scalar = builder.component(point, 0).expect("component");
    let logarithm = builder.ln(scalar).expect("dimensionless logarithm");
    builder
        .objective(logarithm, Sense::Minimize, 1.0)
        .expect("objective");
    let problem = builder.finish();

    assert!(matches!(
        eval(&problem, logarithm, &[vec![-1.0]]),
        Err(OptError::EvalNonFinite {
            node,
            component: None,
            bits,
        }) if node == logarithm.0 && f64::from_bits(bits).is_nan()
    ));

    let mut vector_builder = ProblemBuilder::new();
    let vector_var = vector_builder
        .var("v", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("vector variable");
    let vector = vector_builder.var_ref(vector_var).expect("vector point");
    let doubled = vector_builder.add(vector, vector).expect("vector sum");
    let vector_problem = vector_builder.finish();
    assert!(matches!(
        eval(&vector_problem, doubled, &[vec![f64::MAX, 0.0]]),
        Err(OptError::EvalNonFinite {
            node,
            component: Some(0),
            bits,
        }) if node == doubled.0 && f64::from_bits(bits).is_infinite()
    ));

    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x17,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let opts = fs_opt::DescentOptions {
            steps: 1,
            lr: 0.1,
            fd_h: 1e-6,
        };
        assert!(matches!(
            fs_opt::descend_ir(&problem, &[5e-7], opts, &cx),
            Err(OptError::EvalNonFinite {
                node,
                component: None,
                bits,
            }) if node == logarithm.0 && f64::from_bits(bits).is_nan()
        ));

        let raw = |x: &[f64]| if x[0] < 0.0 { f64::NAN } else { x[0] };
        assert!(matches!(
            fs_opt::descend_fn(Manifold::Rn { dim: 1 }, &raw, &[5e-7], opts, 0, &cx),
            Err(OptError::NonFinite {
                what: "descent objective result",
                bits,
            }) if f64::from_bits(bits).is_nan()
        ));
    });
}

/// adm-018 — retractions refuse non-finite storage, off-manifold bases,
/// and singular candidates instead of normalizing them into fabricated
/// points. Descent applies the same membership gate before f0.
#[test]
#[allow(clippy::too_many_lines)] // one explicit refusal matrix across all manifold families
fn adm_018_retraction_domain_is_fail_closed() {
    let sphere = Manifold::Sphere { ambient: 3 };
    assert!(matches!(
        sphere.retract(&[0.0, 0.0, 0.0], &[0.0; 3]),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            ..
        })
    ));
    assert!(matches!(
        sphere.retract(&[1.0, 0.0, 0.0], &[-1.0, 0.0, 0.0]),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            ..
        })
    ));

    assert!(matches!(
        Manifold::So3.retract(&[0.0; 4], &[0.0; 3]),
        Err(OptError::RetractionDomain {
            manifold: "SO(3)",
            ..
        })
    ));
    assert!(matches!(
        Manifold::So3.retract(&[1.0, 0.0, 0.0, 0.0], &[f64::MAX; 3]),
        Err(OptError::RetractionDomain {
            manifold: "SO(3)",
            ..
        })
    ));

    let stiefel = Manifold::Stiefel { n: 3, p: 2 };
    let orthonormal = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let duplicate_columns = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    assert!(matches!(
        stiefel.retract(&duplicate_columns, &[0.0; 6]),
        Err(OptError::RetractionDomain {
            manifold: "Stiefel",
            ..
        })
    ));
    assert!(matches!(
        stiefel.retract(&orthonormal, &[0.0, 0.0, 0.0, 1.0, -1.0, 0.0]),
        Err(OptError::RetractionDomain {
            manifold: "Stiefel",
            ..
        })
    ));

    assert!(matches!(
        Manifold::Rn { dim: 2 }.retract(&[0.0, f64::INFINITY], &[0.0; 2]),
        Err(OptError::RetractionNonFinite {
            input: "retraction point",
            component: 1,
            ..
        })
    ));
    assert!(matches!(
        Manifold::Rn { dim: 1 }.retract(&[f64::MAX], &[f64::MAX]),
        Err(OptError::RetractionNonFinite {
            input: "retraction output",
            component: 0,
            ..
        })
    ));

    let sphere_out = sphere
        .retract(&[1.0, 0.0, 0.0], &[0.0, 0.25, 0.0])
        .expect("regular sphere candidate");
    let sphere_norm_sq = sphere_out.iter().map(|value| value * value).sum::<f64>();
    assert!((sphere_norm_sq - 1.0).abs() <= 1e-10);
    let so3_out = Manifold::So3
        .retract(&[1.0, 0.0, 0.0, 0.0], &[0.25, 0.0, 0.0])
        .expect("regular SO(3) candidate");
    let so3_norm_sq = so3_out.iter().map(|value| value * value).sum::<f64>();
    assert!((so3_norm_sq - 1.0).abs() <= 1e-10);
    let stiefel_out = stiefel
        .retract(&orthonormal, &[0.0; 6])
        .expect("regular Stiefel candidate");
    assert!(
        stiefel_out
            .iter()
            .zip(orthonormal)
            .all(|(got, expected)| (*got - expected).abs() <= 1e-12)
    );

    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x18,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let calls = std::cell::Cell::new(0u32);
        let guarded = |_: &[f64]| {
            calls.set(calls.get() + 1);
            0.0
        };
        assert!(matches!(
            fs_opt::descend_fn(
                sphere,
                &guarded,
                &[2.0, 0.0, 0.0],
                fs_opt::DescentOptions {
                    steps: 0,
                    ..fs_opt::DescentOptions::default()
                },
                0,
                &cx,
            ),
            Err(OptError::RetractionDomain { .. })
        ));
        assert_eq!(calls.get(), 0, "invalid base must refuse before f0");

        let mut builder = ProblemBuilder::new();
        let variable = builder
            .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("variable");
        let point = builder.var_ref(variable).expect("point");
        let objective = builder.norm_sq(point).expect("objective node");
        builder
            .objective(objective, Sense::Minimize, 1.0)
            .expect("objective");
        let problem = builder.finish();
        assert!(matches!(
            fs_opt::descend_ir(
                &problem,
                &[1.0],
                fs_opt::DescentOptions {
                    steps: 1,
                    lr: f64::MAX,
                    fd_h: 1e-6,
                },
                &cx,
            ),
            Err(OptError::RetractionNonFinite {
                input: "retraction step",
                ..
            })
        ));
    });
}
