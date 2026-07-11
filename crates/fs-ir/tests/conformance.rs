//! fs-ir conformance suite (plan §13.3; the gp3.4 bead): any
//! reimplementation must pass. Acceptance criteria: the Appendix C studies
//! parse as canonical fixtures; s-expr ↔ JSON ↔ AST round-trips are
//! shape-stable (property-tested on generated programs); verb lowering is
//! explicit and inspectable; version pinning round-trips; both parsers are
//! total (fuzz: structured rejections with in-bounds spans, never panics).

use fs_ir::{IR_VERSION, Node, NodeKind, Study, json, lower, sexpr};

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
        let s = sexpr::print(&node);
        let back = sexpr::parse(&s)
            .unwrap_or_else(|e| panic!("round {round}: sexpr reparse failed: {e}\nsrc: {s}"));
        assert!(
            back.same_shape(&node),
            "round {round}: sexpr shape drift\nsrc: {s}"
        );
        // JSON round trip.
        let j = json::print(&node);
        let back = json::parse(&j)
            .unwrap_or_else(|e| panic!("round {round}: json reparse failed: {e}\nsrc: {j}"));
        assert!(
            back.same_shape(&node),
            "round {round}: json shape drift\nsrc: {j}"
        );
        // Cross-syntax: sexpr → AST → json → AST → sexpr → AST.
        let cross = sexpr::parse(&sexpr::print(&json::parse(&json::print(&node)).unwrap()))
            .expect("cross-syntax reparse");
        assert!(cross.same_shape(&node), "round {round}: cross-syntax drift");
    }
    // The fixtures survive the full cross-syntax cycle too.
    for src in [SPOUT, FRAME] {
        let ast = sexpr::parse(src).unwrap();
        let cycled = json::parse(&json::print(&ast)).unwrap();
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
    let printed = sexpr::print(&lowered.node);
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
    verdict(
        "ir-005",
        "verbs lower to explicit IR; defaults named in trace; idempotent",
    );
}

#[test]
fn ir_006_version_pinning_round_trips() {
    assert_eq!(IR_VERSION, 1);
    let src = "(study \"v\" (seed 0x2) (versions (constellation :lock \"2026-07\")))";
    let ast = sexpr::parse(src).unwrap();
    // Through BOTH syntaxes, the pin survives verbatim.
    let via_json = json::parse(&json::print(&ast)).unwrap();
    let via_sexpr = sexpr::parse(&sexpr::print(&via_json)).unwrap();
    let study = Study::from_node(&via_sexpr).unwrap();
    assert_eq!(study.constellation_lock(), Some("2026-07"));
    verdict(
        "ir-006",
        "constellation lock pin survives sexpr->json->sexpr verbatim",
    );
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
                ][next(10) as usize];
                match sexpr::parse(qty).map(|n| n.kind) {
                    Ok(k @ NodeKind::Qty { .. }) => k,
                    _ => unreachable!("qty pool entries always parse"),
                }
            }
            3 => {
                let c = ["384GiB", "96cores", "512MiB", "7KiB", "42B"][next(5) as usize];
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
