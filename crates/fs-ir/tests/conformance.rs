//! fs-ir conformance suite (plan §13.3; the gp3.4 bead): any
//! reimplementation must pass. Acceptance criteria: the Appendix C studies
//! parse as canonical fixtures; s-expr ↔ JSON ↔ AST round-trips are
//! shape-stable (property-tested on generated programs); verb lowering is
//! explicit and inspectable; version pinning round-trips; both parsers are
//! total (fuzz: structured rejections with in-bounds spans, never panics).

use fs_ir::{IR_VERSION, Node, NodeKind, Study, VersionedProgram, json, lower, sexpr};
use fs_qty::Dims;

fn count_qty(n: &Node, hits: &mut usize) {
    match &n.kind {
        NodeKind::Qty { .. } => *hits += 1,
        NodeKind::List(items) => items.iter().for_each(|i| count_qty(i, hits)),
        _ => {}
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ir/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// Appendix C, the spout study (plan §11.1) — verbatim.
const SPOUT: &str = r#"(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let vessel (frep (revolve (cheb-profile "body.chb")) (fillet :edge lip :r 3mm)))
  (let lever  (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (let pour   (flux.free-surface-lbm vessel
                (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
                (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth pour :at lip :modes (1 .. 8))))
  (ascent.optimize J :over lever :method (lbfgs :m 17)
    :until (any (grad-norm 1e-5) (e-value 20) (budget-exhausted))
    :emit (pareto ledger report)))"#;

/// Appendix C, the seismic frame study — verbatim.
const FRAME: &str = r#"(study "frame-seismic-cvar-v9"
  (seed 0xF00D0002) (versions (constellation :lock "2026-07"))
  (capability :cores 96 :mem 384GiB :wall 36h :ops (flux.* ascent.* topo.* uq.*))
  (budget (qoi "P(drift>2e-2)" :rel-error 0.15 :confidence 0.95))
  (let site   (uq.ground-motion (kanai-tajimi :S0 0.03m2/s3 :wg 15rad/s :zg 0.6)
                                (records "PEER-set-A") (mlmc :levels 4)))
  (let ground (topo.ground-structure (grid 8 x 5 x 24m) :knn 14 :rules "AISC-cat.json"))
  (let layout (ascent.solve-lp (min (member-volume ground)) :method pdhg
                               :oracle (michell :tol 0.08)))
  (let frame  (topo.size layout :method tr-newton-krylov
                :constraints ((buckling :code "AISC-E3") (drift-elastic 5e-3))))
  (let resp   (flux.fiber-frame frame site :integrator variational :dt-adapt true))
  (let frag   (uq.probability (exceeds (peak-drift resp) 2e-2)
                :stop (e-process :alpha 0.05)))
  (ascent.optimize (min (mass frame)) :over (sections frame)
    :subject-to ((cvar frag :beta 0.9 :le 0.02) (constructable :catalog "AISC"))
    :method augmented-lagrangian
    :emit (frame frag report ledger)))"#;

#[test]
fn ir_001_appendix_c_studies_parse_as_fixtures() {
    let spout = sexpr::parse(SPOUT).expect("spout study parses");
    let s = Study::from_node(&spout).expect("spout is a study");
    assert_eq!(s.name, "spout-laminar-v3");
    assert_eq!(s.seed, Some(0x5EED_0001));
    assert_eq!(s.constellation_lock(), Some("2026-07"));
    assert_eq!(s.lets.len(), 4);
    assert_eq!(s.lets[0].0, "vessel");
    assert!(s.budget.is_some() && s.versions.is_some());
    assert_eq!(s.body.len(), 1, "one optimize clause");

    let frame = sexpr::parse(FRAME).expect("frame study parses");
    let f = Study::from_node(&frame).expect("frame is a study");
    assert_eq!(f.seed, Some(0xF00D_0002));
    assert!(f.capability.is_some());
    assert_eq!(f.lets.len(), 6);
    // The typed nouns landed as typed nouns: find 0.12Pa*s in the spout.
    let mut hits = 0;
    count_qty(&spout, &mut hits);
    assert!(
        hits >= 7,
        "spout carries at least 7 dimensioned quantities, found {hits}"
    );
    verdict(
        "ir-001",
        "Appendix C spout + frame parse; seeds/locks/lets/nouns extracted",
    );
}

#[test]
fn ir_002_sexpr_json_ast_isomorphism_property() {
    let mut seed = 0x5EED_18AA_0000_0002u64;
    for round in 0..200 {
        let node = gen_node(&mut seed, 0);
        // s-expr round trip.
        let s = sexpr::print(&node).expect("generated AST is valid");
        let back = sexpr::parse(&s)
            .unwrap_or_else(|e| panic!("round {round}: sexpr reparse failed: {e}\nsrc: {s}"));
        assert!(
            back.same_shape(&node),
            "round {round}: sexpr shape drift\nsrc: {s}"
        );
        // JSON round trip.
        let j = json::print(&node).expect("generated AST is valid");
        let back = json::parse(&j)
            .unwrap_or_else(|e| panic!("round {round}: json reparse failed: {e}\nsrc: {j}"));
        assert!(
            back.same_shape(&node),
            "round {round}: json shape drift\nsrc: {j}"
        );
        // Cross-syntax: sexpr → AST → json → AST → sexpr → AST.
        let json_bytes = json::print(&node).expect("generated AST is valid");
        let json_node = json::parse(&json_bytes).expect("JSON reparse");
        let sexpr_bytes = sexpr::print(&json_node).expect("JSON AST is valid");
        let cross = sexpr::parse(&sexpr_bytes).expect("cross-syntax reparse");
        assert!(cross.same_shape(&node), "round {round}: cross-syntax drift");
    }
    // The fixtures survive the full cross-syntax cycle too.
    for src in [SPOUT, FRAME] {
        let ast = sexpr::parse(src).unwrap();
        let cycled = json::parse(&json::print(&ast).expect("fixture AST is valid")).unwrap();
        assert!(cycled.same_shape(&ast), "fixture cross-syntax drift");
    }
    verdict(
        "ir-002",
        "200 generated programs + fixtures: sexpr/json/AST round-trips shape-stable",
    );
}

#[test]
fn ir_003_parsers_are_total_with_in_bounds_spans() {
    let mut seed = 0x5EED_F022_0000_0003u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    let alphabet: Vec<char> = "()\"\\:;.0123456789abcxyzPa*s/-+eE h∞é{}[],#"
        .chars()
        .collect();
    let mut rejects = 0usize;
    for _ in 0..4000 {
        let len = (lcg() % 40) as usize;
        let s: String = (0..len)
            .map(|_| alphabet[(lcg() % alphabet.len() as u64) as usize])
            .collect();
        for result in [sexpr::parse(&s), json::parse(&s)] {
            match result {
                Ok(_) => {}
                Err(e) => {
                    rejects += 1;
                    assert!(
                        e.span.start <= e.span.end && e.span.end <= s.len() + 1,
                        "span out of bounds: {:?} for input len {}",
                        e.span,
                        s.len()
                    );
                    assert!(!e.hint.is_empty(), "refusals must teach: {e}");
                }
            }
        }
    }
    assert!(
        rejects > 1000,
        "garbage battery too easy: only {rejects} rejections"
    );
    // Adversarial nesting: structured rejection, not a stack overflow.
    let deep = "(".repeat(100_000);
    let e = sexpr::parse(&deep).unwrap_err();
    assert_eq!(e.kind.code(), "IrTooDeep");
    let deep_json = "[".repeat(100_000);
    let e = json::parse(&deep_json).unwrap_err();
    assert_eq!(e.kind.code(), "IrTooDeep");
    verdict(
        "ir-003",
        "8000 garbage parses total; structured rejections; depth cap enforced",
    );
}

#[test]
fn ir_004_error_spans_point_at_the_offense() {
    let src = "(study \"x\" (seed 0xZZZ))";
    let e = sexpr::parse(src).unwrap_err();
    assert_eq!(
        &src[e.span.start..e.span.end],
        "0xZZZ",
        "span must isolate the bad seed"
    );
    assert_eq!(e.kind.code(), "IrBadSeed");
    let src = "(budget (wall 2hh))";
    let e = sexpr::parse(src).unwrap_err();
    assert_eq!(&src[e.span.start..e.span.end], "2hh");
    assert_eq!(e.kind.code(), "IrBadQuantity");
    verdict(
        "ir-004",
        "error spans isolate offending tokens; kinds are stable codes",
    );
}

#[test]
fn ir_005_verb_lowering_is_explicit_and_inspectable() {
    let src = "(study \"s\" (seed 0x1) \
               (optimize-shape :min (drag wing) :over (levers wing)) \
               (simulate-pour vessel fluid schedule))";
    let ast = sexpr::parse(src).unwrap();
    let lowered = lower(&ast).expect("lowering succeeds");
    let printed = sexpr::print(&lowered.node).expect("lowered AST is valid");
    // The shorthand is gone; the explicit ops and injected defaults are in.
    assert!(!printed.contains("optimize-shape"));
    assert!(printed.contains("ascent.optimize"));
    assert!(
        printed.contains("(lbfgs :m 17)"),
        "default method must be explicit: {printed}"
    );
    assert!(
        printed.contains("grad-norm"),
        "default stop must be explicit"
    );
    assert!(printed.contains("flux.free-surface-lbm"));
    // The trace names every injection (nothing hidden).
    assert_eq!(lowered.trace.len(), 2);
    let opt = &lowered.trace[0];
    assert_eq!(opt.verb, "optimize-shape");
    assert!(opt.injected.iter().any(|i| i.contains("lbfgs")));
    assert!(!opt.expansion.is_empty());
    // Lowered IR reparses and re-lowers to a fixed point.
    let reparsed = sexpr::parse(&printed).unwrap();
    let again = lower(&reparsed).unwrap();
    assert!(
        again.node.same_shape(&lowered.node),
        "lowering must be idempotent"
    );
    assert!(again.trace.is_empty(), "no verbs remain after lowering");
    // Malformed verb usage refuses with a span.
    let bad = sexpr::parse("(optimize-shape :over x)").unwrap();
    let e = lower(&bad).unwrap_err();
    assert_eq!(e.kind.code(), "IrMalformedClause");
    for (source, expected) in [
        (
            "(optimize-shape :min j :min k :over x)",
            "duplicate optimize-shape argument :min",
        ),
        (
            "(optimize-shape :min j :over x :gpu true)",
            "unknown optimize-shape argument :gpu",
        ),
        (
            "(optimize-shape :min j :over x :method)",
            "dangling :method",
        ),
        (
            "(optimize-shape :min j :over x trailing)",
            "trailing argument",
        ),
        (
            "(simulate-pour vessel fluid schedule trailing)",
            "takes 3 arguments",
        ),
    ] {
        let node = sexpr::parse(source).expect("malformed shorthand still parses");
        let error = lower(&node).expect_err("ambiguous shorthand must not lower");
        assert_eq!(error.kind.code(), "IrMalformedClause");
        assert!(
            error.detail.contains(expected),
            "{source}: expected {expected:?}, got {:?}",
            error.detail
        );
    }
    verdict(
        "ir-005",
        "verbs lower to explicit IR; defaults named; malformed fields refuse; idempotent",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One exact-count authority and identity matrix.
fn ir_0xx_exact_count_literals_survive_admission_and_identity() {
    use fs_ir::ast::CountValue;
    // gp3.20: adjacent large integers stay distinguishable end to end.
    let lo = sexpr::parse("9007199254740992B").expect("2^53 parses");
    let hi = sexpr::parse("9007199254740993B").expect("2^53+1 parses");
    assert!(
        !lo.same_shape(&hi),
        "2^53 and 2^53+1 bytes must be distinct identities"
    );
    // Canonical print preserves the exact digits (identity binding).
    assert_eq!(
        sexpr::print(&hi).expect("count AST is valid"),
        "9007199254740993B"
    );
    assert_eq!(
        json::print(&hi).expect("count AST is valid"),
        "{\"c\":\"9007199254740993B\"}",
        "exact digits bind into the JSON identity too"
    );
    // Round trips are lossless in both syntaxes.
    assert!(
        sexpr::parse(&sexpr::print(&hi).expect("count AST is valid"))
            .expect("reparse")
            .same_shape(&hi)
    );
    assert!(
        json::parse(&json::print(&hi).expect("count AST is valid"))
            .expect("reparse")
            .same_shape(&hi)
    );
    // Exact scaling is checked in wide arithmetic: u64::MAX bytes is
    // admissible; one more is refused (overflow BEFORE rounding).
    let max = sexpr::parse("18446744073709551615B").expect("u64::MAX parses");
    let NodeKind::Count { value, unit } = &max.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), Some(u64::MAX));
    let over = sexpr::parse("18446744073709551616B").expect("2^64 parses as exact u128");
    let NodeKind::Count { value, unit } = &over.kind else {
        panic!("count expected");
    };
    assert_eq!(
        value.integral_bytes(*unit),
        None,
        "2^64 B refuses, never rounds"
    );
    // GiB scaling near the boundary: 2^34 GiB = 2^64 bytes refuses;
    // (2^34 - 1) GiB fits.
    let fits = sexpr::parse("17179869183GiB").expect("parses");
    let NodeKind::Count { value, unit } = &fits.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), Some((17_179_869_183u64) << 30));
    let over_gib = sexpr::parse("17179869184GiB").expect("parses");
    let NodeKind::Count { value, unit } = &over_gib.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), None);
    // Beyond u128: a structured refusal at parse, not saturation.
    let err = sexpr::parse("340282366920938463463374607431768211456B")
        .expect_err("2^128 cannot be an exact count");
    assert!(err.to_string().contains("exact count range"));
    let hostile = format!("{}B", "9".repeat(100_000));
    let err = sexpr::parse(&hostile).expect_err("oversized exact count must refuse");
    assert!(err.to_string().contains("exact count range"));
    assert!(err.span.end <= hostile.len());
    assert!(
        err.detail.len() < 256,
        "hostile token must not be copied wholesale into diagnostics"
    );
    // Decimal/exponent semantics are exact too: useful fractional-unit
    // spellings work without binary-float authority decisions.
    let frac = sexpr::parse("1.5GiB").expect("parses");
    let NodeKind::Count { value, unit } = &frac.kind else {
        panic!("count expected");
    };
    assert!(matches!(value, CountValue::Fractional(_)));
    assert_eq!(value.integral_bytes(*unit), Some(3 << 29));
    assert_eq!(sexpr::print(&frac).expect("count AST is valid"), "15e-1GiB");
    assert!(
        sexpr::parse(&sexpr::print(&frac).expect("count AST is valid"))
            .expect("decimal canonical reparse")
            .same_shape(&frac)
    );
    let exponent = sexpr::parse("1e3B").expect("exponent count parses");
    let NodeKind::Count { value, unit } = &exponent.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), Some(1_000));
    assert!(
        json::parse(&json::print(&exponent).expect("count AST is valid"))
            .expect("exponent JSON reparse")
            .same_shape(&exponent)
    );
    let half_kib = sexpr::parse("0.5KiB").expect("fractional unit parses");
    let NodeKind::Count { value, unit } = &half_kib.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), Some(512));
    let reduced =
        sexpr::parse("316912650130689144134521484375e-30GiB").expect("large exact decimal parses");
    let NodeKind::Count { value, unit } = &reduced.kind else {
        panic!("count expected");
    };
    assert_eq!(
        value.integral_bytes(*unit),
        Some(340_282_367),
        "denominator factors must cancel before checked scaling"
    );

    // Exact decimal storage distinguishes an integer-valued decimal above
    // 2^53 from its neighbors without guessing through f64.
    let big_frac = sexpr::parse("9007199254740993.0B").expect("parses");
    let NodeKind::Count { value, unit } = &big_frac.kind else {
        panic!("count expected");
    };
    assert_eq!(value.integral_bytes(*unit), Some(9_007_199_254_740_993));
    for fractional_byte in [
        "0.1B",
        "0.99999999999999999B",
        "1.00000000000000001B",
        "1e-1B",
    ] {
        let parsed = sexpr::parse(fractional_byte).expect("bounded decimal parses");
        let NodeKind::Count { value, unit } = parsed.kind else {
            panic!("count expected");
        };
        assert_eq!(
            value.integral_bytes(unit),
            None,
            "{fractional_byte} must not round into a whole-byte claim"
        );
    }
    // Mixed written forms are distinct claims.
    let int_form = sexpr::parse("2B").expect("parses");
    let frac_form = sexpr::parse("2.0B").expect("parses");
    assert!(!int_form.same_shape(&frac_form));
}

#[test]
fn ir_006_version_pinning_round_trips() {
    assert_eq!(IR_VERSION, 3);
    let src = "(study \"v\" (seed 0x2) (versions (constellation :lock \"2026-07\")))";
    let ast = sexpr::parse(src).unwrap();
    // Through BOTH syntaxes, the pin survives verbatim.
    let via_json = json::parse(&json::print(&ast).expect("fixture AST is valid")).unwrap();
    let via_sexpr = sexpr::parse(&sexpr::print(&via_json).expect("JSON AST is valid")).unwrap();
    let study = Study::from_node(&via_sexpr).unwrap();
    assert_eq!(study.constellation_lock(), Some("2026-07"));

    let artifact = VersionedProgram::current(ast);
    let canonical_sexpr = artifact.print_sexpr();
    let canonical_json = artifact.print_json();
    let from_sexpr = VersionedProgram::parse_sexpr(&canonical_sexpr).expect("v3 envelope");
    let from_json = VersionedProgram::parse_json(&canonical_json).expect("v3 JSON envelope");
    assert_eq!(from_sexpr.version(), IR_VERSION);
    assert!(from_sexpr.program().same_shape(from_json.program()));
    assert_eq!(from_sexpr.print_sexpr(), canonical_sexpr);
    assert_eq!(from_json.print_json(), canonical_json);

    for unsupported in [1, 2, IR_VERSION + 1] {
        let source = canonical_sexpr.replacen(
            &format!(":version {IR_VERSION}"),
            &format!(":version {unsupported}"),
            1,
        );
        let error = VersionedProgram::parse_sexpr(&source)
            .expect_err("unsupported language semantics must refuse");
        assert_eq!(error.kind.code(), "IrUnsupportedVersion");
        assert_eq!(
            &source[error.span.start..error.span.end],
            unsupported.to_string()
        );
    }
    let newer_json = canonical_json.replacen(
        &format!("{{\"i\":{IR_VERSION}}}"),
        &format!("{{\"i\":{}}}", IR_VERSION + 1),
        1,
    );
    assert_eq!(
        VersionedProgram::parse_json(&newer_json)
            .expect_err("newer JSON envelope must refuse")
            .kind
            .code(),
        "IrUnsupportedVersion"
    );

    let legacy_count = sexpr::parse("384.0GiB").expect("legacy bare syntax parses syntax-only");
    let migrated = VersionedProgram::current(legacy_count);
    assert!(
        migrated.print_sexpr().contains("384e0GiB"),
        "migration is explicit re-emission under the v3 envelope"
    );
    let mol = sexpr::parse("2mol").expect("sixth-base quantity parses");
    let NodeKind::Qty { dims, .. } = mol.kind else {
        panic!("quantity expected");
    };
    assert_eq!(dims, Dims([0, 0, 0, 0, 0, 1]));
    verdict(
        "ir-006",
        "constellation lock and six-base IR v3 envelope survive both syntaxes; stale versions refuse",
    );
}

#[test]
fn ir_006b_quantity_budget_refusal_is_bounded_and_deterministic() {
    let source = format!("1{}", "x".repeat(5_000));
    let first = sexpr::parse(&source).expect_err("oversized quantity token must refuse");
    let second = sexpr::parse(&source).expect_err("repeat must refuse identically");
    assert_eq!(first, second);
    assert_eq!(first.kind.code(), "IrBadQuantity");
    assert_eq!(first.span.start, 0);
    assert_eq!(first.span.end, source.len());
    assert!(
        first.detail.contains("InputBytes")
            && first.detail.contains("unavailable-before-byte-admission")
    );
    assert!(
        first.detail.len() < 1_024,
        "bounded fs-qty diagnostics must not retain the 5 KiB source: {} bytes",
        first.detail.len()
    );
    let normal = sexpr::parse("2mol").expect("ordinary quantity remains admitted");
    assert!(matches!(normal.kind, NodeKind::Qty { .. }));
    verdict(
        "ir-006b",
        "fs-ir explicitly uses the bounded fs-qty entry point and retains bounded deterministic diagnostics",
    );
}

#[test]
fn ir_006c_json_numbers_and_strings_are_strict_rfc_8259() {
    for source in [
        r#"{"f":0}"#,
        r#"{"f":-0}"#,
        r#"{"f":1}"#,
        r#"{"f":0.1}"#,
        r#"{"f":1e2}"#,
        r#"{"f":1E+2}"#,
        r#"{"f":1e-2}"#,
    ] {
        let node = json::parse(source).unwrap_or_else(|error| panic!("{source}: {error}"));
        let first = json::print_checked(&node).expect("canonical JSON");
        let second = json::print_checked(&json::parse(&first).expect("canonical reparse"))
            .expect("canonical JSON");
        assert_eq!(first, second, "canonical bytes must be idempotent");
    }

    for source in [
        r#"{"f":+1}"#,
        r#"{"f":01}"#,
        r#"{"f":-01}"#,
        r#"{"f":.1}"#,
        r#"{"f":1.}"#,
        r#"{"f":1e}"#,
        r#"{"f":1e+}"#,
        r#"{"f":--1}"#,
        r#"{"f":1e9999}"#,
    ] {
        let error = json::parse(source).expect_err("non-RFC/nonfinite number must refuse");
        assert_eq!(error.kind.code(), "IrBadNumber", "{source}: {error}");
    }
    let leading_zero = r#"{"f":01}"#;
    let error = json::parse(leading_zero).expect_err("leading zero must refuse");
    assert_eq!(&leading_zero[error.span.start..error.span.end], "01");

    let escaped = json::parse(r#"{"s":"line\n\t\r\b\f"}"#).expect("escaped controls parse");
    assert!(matches!(escaped.kind, NodeKind::Str(_)));
    let rocket = json::parse(r#"{"s":"\uD83D\uDE80"}"#).expect("surrogate pair parses");
    assert!(matches!(&rocket.kind, NodeKind::Str(text) if text == "🚀"));
    for source in [
        r#"{"s":"\uD83D"}"#,
        r#"{"s":"\uDE80"}"#,
        r#"{"s":"\uD83D\u0041"}"#,
    ] {
        assert_eq!(
            json::parse(source)
                .expect_err("unpaired surrogate must refuse")
                .kind
                .code(),
            "IrBadEscape"
        );
    }
    let raw_control = "{\"s\":\"line\nbreak\"}";
    let error = json::parse(raw_control).expect_err("raw JSON control must refuse");
    assert_eq!(error.kind.code(), "IrJsonSyntax");
    assert_eq!(&raw_control[error.span.start..error.span.end], "\n");

    for source in [r#"{"sym":"1"}"#, r#"{"kw":""}"#] {
        let error = json::parse(source).expect_err("cross-syntax-invalid atom must refuse");
        assert!(error.detail.contains("invalid AST at $"));
    }
}

// ---------------------------------------------------------------------------
// Random AST generator (seeded LCG; atoms drawn from the real noun pool)
// ---------------------------------------------------------------------------

fn gen_node(seed: &mut u64, depth: usize) -> Node {
    let mut next = |m: u64| {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *seed % m
    };
    let leaf = depth >= 4 || next(3) > 0;
    if leaf {
        let kind = match next(8) {
            0 => NodeKind::Int(next(2_000_000).cast_signed() - 1_000_000),
            1 => {
                let v = (next(1_000_000) as f64 / 997.0) - 500.0;
                NodeKind::Float(v)
            }
            2 => {
                let qty = [
                    "0.12Pa*s",
                    "65deg",
                    "3s",
                    "0.061N/m",
                    "15rad/s",
                    "0.03m2/s3",
                    "2e-2m",
                    "5mm",
                    "0.5L/s",
                    "36h",
                    "2mol",
                ][next(11) as usize];
                match sexpr::parse(qty).map(|n| n.kind) {
                    Ok(k @ NodeKind::Qty { .. }) => k,
                    _ => unreachable!("qty pool entries always parse"),
                }
            }
            3 => {
                let counts = [
                    "384GiB", "96cores", "512MiB", "7KiB", "42B", "1.5GiB", "0.5KiB", "1e3B",
                    "1e-21B",
                ];
                let c = counts[next(counts.len() as u64) as usize];
                match sexpr::parse(c).map(|n| n.kind) {
                    Ok(k @ NodeKind::Count { .. }) => k,
                    _ => unreachable!("count pool entries always parse"),
                }
            }
            4 => NodeKind::Seed(next(u64::MAX)),
            5 => NodeKind::Str(gen_text(&mut next, true)),
            6 => NodeKind::Keyword(gen_ident(&mut next)),
            _ => NodeKind::Symbol(gen_ident(&mut next)),
        };
        return Node::synthetic(kind);
    }
    let n = 1 + next(5) as usize;
    let items = (0..n).map(|_| gen_node(seed, depth + 1)).collect();
    Node::synthetic(NodeKind::List(items))
}

fn gen_ident(next: &mut dyn FnMut(u64) -> u64) -> String {
    let pool = [
        "ascent.optimize",
        "frep",
        "vessel",
        "grad-norm",
        "e-value",
        "mlmc",
        "flux.fiber-frame",
        "budget-exhausted",
        "x",
        "wall",
    ];
    pool[next(pool.len() as u64) as usize].to_string()
}

fn gen_text(next: &mut dyn FnMut(u64) -> u64, escapes: bool) -> String {
    let pool: &[&str] = if escapes {
        &[
            "body.chb",
            "PEER-set-A",
            "with \"quotes\"",
            "tab\there",
            "line\nbreak",
            "π≈3",
            "",
        ]
    } else {
        &["plain"]
    };
    pool[next(pool.len() as u64) as usize].to_string()
}
