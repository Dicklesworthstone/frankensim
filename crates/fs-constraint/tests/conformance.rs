//! fs-constraint conformance suite (CONTRACT.md: any reimplementation
//! must pass). Kind taxonomy + treatments + serialization, evidence
//! statuses/roles/penalties, chance validity machinery (the BOUND
//! decides, not the raw rate), certification refusals + real interval
//! proofs, minimal unsat cores vs enumeration, and the worked repair
//! example with calibrated feasibility estimates. Completed aggregate
//! cases emit canonical fs-obs verdicts. Randomized input campaigns carry
//! their literal base seeds; fixed-input cases use zero, and the fixed Cx
//! stream remains separate execution provenance. Assertions and expectations
//! reached before a verdict remain ordinary Rust test diagnostics.

use asupersync::types::Budget;
use fs_constraint::{
    ChanceEstimator, ConError, ConstraintKind, ConstraintSpec, Diagnosis, DomainBox, DomainError,
    DomainRangeError, ElasticReport, PenaltyLaw, ProofKind, RepairAction, RepairKind, Status,
    Treatment, diagnose_infeasibility, elastic_solve, evaluate, interval_eval, parse_specs,
    prove_interval, serialize_specs,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{AdmissionCaps, Manifold, NodeId, Problem, ProblemBuilder};
use fs_qty::Dims;
use std::cell::Cell;

const FIXED_INPUT_SEED: u64 = 0;
const EXECUTION_SEED: u64 = 0xC0C0;
const FSCON_003_INPUT_SEED: u64 = 0x1001_2026_0707_0003;
const FSCON_003_STREAM_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;
const FSCON_004_INPUT_SEED: u64 = 0x1001_2026_0707_0004;
const FSCON_005_INPUT_SEED: u64 = 0x1001_2026_0707_0005;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-constraint/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-constraint/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("constraint verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("constraint verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Host problem: one Rn(2) variable; linear forms `a·x − b` as nodes.
struct Host {
    problem: Problem,
    nodes: Vec<NodeId>,
}

/// Build `a0·x0 + a1·x1 − b ≤ 0` constraint nodes.
fn linear_host(rows: &[(f64, f64, f64)]) -> Host {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("ref");
    let x0 = b.component(vr, 0).expect("x0");
    let x1 = b.component(vr, 1).expect("x1");
    let mut nodes = Vec::new();
    for &(a0, a1, rhs) in rows {
        let c0 = b.konst(a0, Dims::NONE).expect("finite konst");
        let c1 = b.konst(a1, Dims::NONE).expect("finite konst");
        let t0 = b.mul(c0, x0).expect("t0");
        let t1 = b.mul(c1, x1).expect("t1");
        let s = b.add(t0, t1).expect("s");
        let rb = b.konst(rhs, Dims::NONE).expect("finite konst");
        nodes.push(b.sub(s, rb).expect("g"));
    }
    // Anchor an objective so the problem is well-formed.
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    Host {
        problem: b.finish(),
        nodes,
    }
}

fn hard(name: &str, node: NodeId) -> ConstraintSpec {
    ConstraintSpec {
        name: name.to_string(),
        node,
        kind: ConstraintKind::Hard,
        active_tol: 1e-9,
    }
}

/// fscon-001 — taxonomy: every kind maps to its optimizer treatment;
/// the spec set round-trips through canonical serialization; ledger
/// rows validate through fs-obs.
#[test]
fn fscon_001_taxonomy_and_roundtrip() {
    let host = linear_host(&[(1.0, 0.0, 1.0), (0.0, 1.0, 1.0)]);
    let specs = vec![
        hard("wall", host.nodes[0]),
        ConstraintSpec {
            name: "pretty please".to_string(),
            node: host.nodes[1],
            kind: ConstraintKind::Soft(PenaltyLaw::Quadratic { weight: 3.5 }),
            active_tol: 1e-9,
        },
        ConstraintSpec {
            name: "yield-prob".to_string(),
            node: host.nodes[0],
            kind: ConstraintKind::Chance {
                level: 0.9,
                estimator: ChanceEstimator::MonteCarlo {
                    samples: 256,
                    delta: 0.05,
                },
            },
            active_tol: 1e-9,
        },
        ConstraintSpec {
            name: "load-envelope".to_string(),
            node: host.nodes[1],
            kind: ConstraintKind::Robust {
                half_widths: vec![0.1, 0.05],
            },
            active_tol: 1e-9,
        },
        ConstraintSpec {
            name: "stress-proof".to_string(),
            node: host.nodes[0],
            kind: ConstraintKind::Certification {
                proof: ProofKind::Interval,
            },
            active_tol: 1e-9,
        },
        ConstraintSpec {
            name: "min-wall".to_string(),
            node: host.nodes[1],
            kind: ConstraintKind::Fabrication {
                process: "cnc 3axis".to_string(),
            },
            active_tol: 1e-9,
        },
        ConstraintSpec {
            name: "member-slenderness".to_string(),
            node: host.nodes[0],
            kind: ConstraintKind::Code {
                standard: "aisc-360".to_string(),
            },
            active_tol: 1e-9,
        },
    ];
    let treatments_ok = specs[0].kind.treatment() == Treatment::FeasibilityRestoration
        && specs[1].kind.treatment() == Treatment::PenaltyTerm
        && specs[2].kind.treatment() == Treatment::EstimateThenBound
        && specs[3].kind.treatment() == Treatment::ProveOrEscalate
        && specs[4].kind.treatment() == Treatment::ProveOrEscalate
        && specs[5].kind.treatment() == Treatment::DomainCheck
        && specs[6].kind.treatment() == Treatment::DomainCheck;
    let text = serialize_specs(&specs);
    let back = parse_specs(&text).expect("round-trip");
    let roundtrip = back == specs;
    let garbage = parse_specs("fscon v2\nwat\n");
    let teaches = matches!(garbage, Err(ConError::Parse { line: 2, .. }));
    // Ledger row through fs-obs.
    let ev = evaluate(&host.problem, &specs[0], &[0.2, 0.2], None).expect("eval");
    let mut em = fs_obs::Emitter::new("fs-constraint/conformance", "fscon-001/ledger");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "constraint-ledger-row".to_string(),
                json: ev.to_ledger_row(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("ledger row validates");
    println!("{line}");
    verdict(
        "fscon-001",
        treatments_ok && roundtrip && teaches,
        "all seven kinds map to their optimizer treatments, the spec set \
         round-trips through canonical serialization IDENTICALLY, garbage refuses \
         with its line number, and ledger rows validate through fs-obs",
        FIXED_INPUT_SEED,
    );
}

/// fscon-002 — evaluation evidence: statuses, active-set roles, exact
/// violation certificates, and penalty laws.
#[test]
fn fscon_002_evidence() {
    let host = linear_host(&[(1.0, 1.0, 1.0)]); // x0 + x1 − 1 ≤ 0
    let spec_hard = hard("sum-cap", host.nodes[0]);
    let sat = evaluate(&host.problem, &spec_hard, &[0.2, 0.3], None).expect("sat");
    let act = evaluate(&host.problem, &spec_hard, &[0.5, 0.5], None).expect("act");
    let vio = evaluate(&host.problem, &spec_hard, &[0.9, 0.4], None).expect("vio");
    let statuses = sat.status == Status::Satisfied
        && act.status == Status::Active
        && vio.status == Status::Violated;
    let roles = sat.role == fs_constraint::ActiveRole::Inactive
        && act.role == fs_constraint::ActiveRole::Active
        && vio.role == fs_constraint::ActiveRole::Violating;
    let exact_mag = (vio.violation - 0.3).abs() < 1e-12
        && vio.certificate.kind == fs_evidence::NumericalKind::Exact;
    let soft = ConstraintSpec {
        kind: ConstraintKind::Soft(PenaltyLaw::Quadratic { weight: 2.0 }),
        ..spec_hard.clone()
    };
    let pen = evaluate(&host.problem, &soft, &[0.9, 0.4], None).expect("pen");
    let hinge = ConstraintSpec {
        kind: ConstraintKind::Soft(PenaltyLaw::Hinge { weight: 2.0 }),
        ..spec_hard
    };
    let pen2 = evaluate(&host.problem, &hinge, &[0.9, 0.4], None).expect("pen2");
    let laws =
        (pen.penalty - 2.0 * 0.3 * 0.3).abs() < 1e-12 && (pen2.penalty - 2.0 * 0.3).abs() < 1e-12;
    verdict(
        "fscon-002",
        statuses && roles && exact_mag && laws,
        "statuses and active-set roles classify correctly, violation magnitudes \
         carry EXACT certificates, and both penalty laws price violations as \
         declared",
        FIXED_INPUT_SEED,
    );
}

/// fscon-003 — chance validity machinery: an analytic uniform-noise
/// case where the raw empirical rate clears the level but the Hoeffding
/// BOUND does not — and the constraint REFUSES to claim satisfied.
#[test]
fn fscon_003_chance_bound_decides() {
    let host = linear_host(&[(1.0, 0.0, 1.0)]); // x0 − 1 ≤ 0
    let spec = |level: f64| ConstraintSpec {
        name: "chance-cap".to_string(),
        node: host.nodes[0],
        kind: ConstraintKind::Chance {
            level,
            estimator: ChanceEstimator::MonteCarlo {
                samples: 400,
                delta: 0.05,
            },
        },
        active_tol: 1e-9,
    };
    // Noise: u ~ U(0,1) on x0 (deterministic per-sample stream).
    // At x0 = 0.2: P(x0 + u ≤ 1) = P(u ≤ 0.8) = 0.8 exactly.
    let noise = |s: u64| -> Vec<f64> {
        let mut r = Lcg(FSCON_003_INPUT_SEED ^ (s.wrapping_mul(FSCON_003_STREAM_STRIDE)));
        vec![r.unit(), 0.0]
    };
    let x = [0.2, 0.0];
    let comfortable = evaluate(&host.problem, &spec(0.70), &x, Some(&noise)).expect("comfortable");
    let squeezed = evaluate(&host.problem, &spec(0.78), &x, Some(&noise)).expect("squeezed");
    let hopeless = evaluate(&host.problem, &spec(0.95), &x, Some(&noise)).expect("hopeless");
    // Half-width at n=400, delta=0.05: sqrt(ln 20 / 800) ≈ 0.0612.
    let comfortable_ok = comfortable.status == Status::Satisfied;
    let squeezed_ok = matches!(
        squeezed.status,
        Status::BoundNotCleared { empirical, lower_bound }
            if empirical >= 0.78 && lower_bound < 0.78
    );
    let hopeless_ok = hopeless.status == Status::Violated;
    let stat_carried = matches!(
        comfortable.statistical,
        fs_evidence::StatisticalCertificate::HalfWidth { confidence, .. }
            if (confidence - 0.95).abs() < 1e-12
    );
    verdict(
        "fscon-003",
        comfortable_ok && squeezed_ok && hopeless_ok && stat_carried,
        &format!(
            "the BOUND decides: level 0.70 satisfied, level 0.78 refused as \
             BoundNotCleared even though the raw rate clears it ({squeezed:?} \
             status), level 0.95 violated; the Hoeffding half-width travels as a \
             StatisticalCertificate; input seed 0x1001_2026_0707_0003, with \
             sample stream s derived as seed ^ \
             s.wrapping_mul(0x9E37_79B9_7F4A_7C15)",
            squeezed = squeezed.status
        ),
        FSCON_003_INPUT_SEED,
    );
}

/// G0: the public evaluator's v1 uncertainty model is Euclidean. Unsupported
/// host shapes must refuse before either the graph or noise model observes
/// work; manifold-aware chance perturbations are deliberately not inferred.
#[test]
fn evaluate_admits_the_single_rn_host_before_graph_or_noise_work() {
    fn problem_with(manifolds: &[Manifold]) -> Problem {
        let mut builder = ProblemBuilder::new();
        for (index, manifold) in manifolds.iter().copied().enumerate() {
            builder
                .var(&format!("x{index}"), manifold, Dims::NONE)
                .expect("valid host manifold fixture");
        }
        let objective = builder
            .konst(0.0, Dims::NONE)
            .expect("finite scalar objective");
        builder
            .objective(objective, fs_opt::Sense::Minimize, 1.0)
            .expect("objective entry");
        builder.finish()
    }

    // The forged node is a tripwire: reaching fs-opt graph evaluation would
    // return Eval instead of the host-specific InvalidDomain below.
    let spec = ConstraintSpec {
        name: "host-admission".to_string(),
        node: NodeId(u32::MAX),
        kind: ConstraintKind::Chance {
            level: 0.5,
            estimator: ChanceEstimator::MonteCarlo {
                samples: 1,
                delta: 0.05,
            },
        },
        active_tol: 0.0,
    };
    let noise_calls = Cell::new(0_usize);
    let noise = |_sample: u64| {
        noise_calls.set(noise_calls.get() + 1);
        Vec::new()
    };

    for (manifold, point) in [
        (Manifold::Sphere { ambient: 2 }, vec![1.0, 0.0]),
        (Manifold::So3, vec![1.0, 0.0, 0.0, 0.0]),
        (Manifold::Stiefel { n: 2, p: 1 }, vec![1.0, 0.0]),
    ] {
        let problem = problem_with(&[manifold]);
        assert_eq!(
            evaluate(&problem, &spec, &point, Some(&noise)).unwrap_err(),
            ConError::InvalidDomain(DomainError::HostVariableManifold { got: manifold })
        );
        assert_eq!(noise_calls.get(), 0, "{manifold:?} invoked noise");
    }

    let multi = problem_with(&[Manifold::Rn { dim: 1 }, Manifold::Rn { dim: 1 }]);
    assert_eq!(
        evaluate(&multi, &spec, &[0.0, 0.0], Some(&noise)).unwrap_err(),
        ConError::InvalidDomain(DomainError::HostVariableCount { got: 2 })
    );
    assert_eq!(noise_calls.get(), 0, "multi-variable host invoked noise");

    let rn = problem_with(&[Manifold::Rn { dim: 2 }]);
    assert_eq!(
        evaluate(&rn, &spec, &[0.0], Some(&noise)).unwrap_err(),
        ConError::InvalidDomain(DomainError::DimensionMismatch {
            expected: 2,
            got: 1,
        })
    );
    assert_eq!(noise_calls.get(), 0, "short Rn point invoked noise");
}

/// fscon-004 — certification refusals + REAL interval proofs (and the
/// G0 containment law for the interval engine), plus robust kinds
/// proven over uncertainty boxes.
#[test]
fn fscon_004_certification_and_robust() {
    let host = linear_host(&[(1.0, 1.0, 3.0)]); // x0 + x1 − 3 ≤ 0
    let cert = ConstraintSpec {
        name: "stress-proof".to_string(),
        node: host.nodes[0],
        kind: ConstraintKind::Certification {
            proof: ProofKind::Interval,
        },
        active_tol: 1e-9,
    };
    // Pointwise deep inside — still NeedsProof without an artifact.
    let ev = evaluate(&host.problem, &cert, &[0.1, 0.1], None).expect("eval");
    let refuses = matches!(
        ev.status,
        Status::NeedsProof {
            proof: ProofKind::Interval
        }
    );
    // Interval proof over [0,1]²: hi = 2 − 3 = −1 ≤ 0 → PROVEN.
    let (proven, artifact) =
        prove_interval(&host.problem, &cert, &[(0.0, 1.0), (0.0, 1.0)]).expect("proof");
    let proven_ok = proven.status == Status::Proven
        && artifact.interval_bound().hi <= 0.0
        && artifact.is_bound_to(&host.problem, cert.node, &[(0.0, 1.0), (0.0, 1.0)]);
    // Unprovable domain: [0,2]² has hi = 1 > 0 → honest refusal.
    let refused = prove_interval(&host.problem, &cert, &[(0.0, 2.0), (0.0, 2.0)]);
    let honest = matches!(refused, Err(ConError::NotProvable { ref why }) if why.contains("upper"));
    // G0 containment: random nonlinear graph, random boxes, samples in.
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("r");
    let n = b.norm_sq(vr).expect("n");
    let x0 = b.component(vr, 0).expect("x0");
    let t = b.tanh(x0).expect("t");
    let m = b.mul(n, n).expect("m");
    let s = b.add(m, t).expect("s");
    let a = b.abs(s).expect("a");
    let one = b.konst(1.0, Dims::NONE).expect("finite konst");
    let g = b.min_of(a, one).expect("g");
    b.objective(g, fs_opt::Sense::Minimize, 1.0).expect("o");
    let nl = b.finish();
    let mut rng = Lcg(FSCON_004_INPUT_SEED);
    let mut contained = true;
    for _ in 0..300 {
        let c = [rng.range(-1.5, 1.5), rng.range(-1.5, 1.5)];
        let h = [rng.range(0.05, 0.8), rng.range(0.05, 0.8)];
        let boxes = [(c[0] - h[0], c[0] + h[0]), (c[1] - h[1], c[1] + h[1])];
        let iv = interval_eval(&nl, g, &boxes).expect("interval");
        for _ in 0..10 {
            let p = vec![
                rng.range(boxes[0].0, boxes[0].1),
                rng.range(boxes[1].0, boxes[1].1),
            ];
            let val = fs_opt::eval(&nl, g, std::slice::from_ref(&p))
                .expect("eval")
                .scalar()
                .expect("s");
            contained &= val >= iv.lo - 1e-9 && val <= iv.hi + 1e-9;
        }
    }
    // Robust kind: proven at a safe point, violated near the edge.
    let robust = ConstraintSpec {
        name: "load-envelope".to_string(),
        node: host.nodes[0],
        kind: ConstraintKind::Robust {
            half_widths: vec![0.2, 0.2],
        },
        active_tol: 1e-9,
    };
    let safe = evaluate(&host.problem, &robust, &[1.0, 1.0], None).expect("safe");
    let edgy = evaluate(&host.problem, &robust, &[1.4, 1.4], None).expect("edgy");
    let robust_ok = safe.status == Status::Proven
        && edgy.status == Status::Violated
        && safe.certificate.kind == fs_evidence::NumericalKind::Enclosure;
    verdict(
        "fscon-004",
        refuses && proven_ok && honest && contained && robust_ok,
        "certification kinds refuse satisfied without an artifact (pointwise \
         goodness is not a proof), the interval engine PROVES over provable \
         domains and refuses honestly otherwise, containment holds over 300 \
         random nonlinear boxes x 10 samples, and robust kinds carry enclosure \
         certificates; seed 0x1001_2026_0707_0004",
        FSCON_004_INPUT_SEED,
    );
}

/// Grid-enumeration feasibility of a subset (ground truth).
fn grid_feasible(
    problem: &Problem,
    specs: &[ConstraintSpec],
    subset: &[usize],
    domain: &DomainBox,
) -> bool {
    let n = 80;
    for i in 0..=n {
        for j in 0..=n {
            let x = vec![
                domain.ranges[0].0
                    + (domain.ranges[0].1 - domain.ranges[0].0) * f64::from(i) / f64::from(n),
                domain.ranges[1].0
                    + (domain.ranges[1].1 - domain.ranges[1].0) * f64::from(j) / f64::from(n),
            ];
            let ok = subset.iter().all(|&k| {
                fs_opt::eval(problem, specs[k].node, std::slice::from_ref(&x))
                    .expect("eval")
                    .scalar()
                    .expect("s")
                    <= 1e-9
            });
            if ok {
                return true;
            }
        }
    }
    false
}

/// fscon-005 — minimal unsat cores: seeded fixtures verified against
/// enumeration — the FULL set is infeasible, the core is infeasible,
/// and dropping ANY core member restores feasibility (G0 minimality).
#[test]
fn fscon_005_unsat_cores() {
    with_cx(|cx| {
        // Triangle infeasibility: x+y ≥ 3, x ≤ 1, y ≤ 1 (+ a bystander).
        let host = linear_host(&[
            (-1.0, -1.0, -3.0), // 3 − x − y ≤ 0  ⇔  x+y ≥ 3
            (1.0, 0.0, 1.0),    // x ≤ 1
            (0.0, 1.0, 1.0),    // y ≤ 1
            (-1.0, 0.0, 5.0),   // x ≥ −5 (satisfiable bystander)
        ]);
        let specs: Vec<ConstraintSpec> = ["need-sum", "cap-x", "cap-y", "floor-x"]
            .iter()
            .zip(&host.nodes)
            .map(|(n, &node)| hard(n, node))
            .collect();
        let domain = DomainBox {
            ranges: vec![(-5.0, 5.0), (-5.0, 5.0)],
        };
        let diag = diagnose_infeasibility(&host.problem, &specs, &domain, cx).expect("diag");
        let core_right = !diag.feasible && diag.core == vec![0, 1, 2];
        // G0 minimality vs enumeration: the core is infeasible; every
        // deletion is feasible.
        let full_infeasible = !grid_feasible(&host.problem, &specs, &diag.core, &domain);
        let mut deletions_feasible = true;
        for &drop in &diag.core {
            let rest: Vec<usize> = diag.core.iter().copied().filter(|&i| i != drop).collect();
            deletions_feasible &= grid_feasible(&host.problem, &specs, &rest, &domain);
        }
        // A feasible system reports a witness and no core.
        let feasible_host = linear_host(&[(1.0, 0.0, 1.0), (0.0, 1.0, 1.0)]);
        let fspecs: Vec<ConstraintSpec> = ["a", "b"]
            .iter()
            .zip(&feasible_host.nodes)
            .map(|(n, &node)| hard(n, node))
            .collect();
        let fd = diagnose_infeasibility(&feasible_host.problem, &fspecs, &domain, cx)
            .expect("feasible diag");
        let feasible_ok = fd.feasible && fd.core.is_empty() && fd.witness.is_some();
        // Seeded random family: elastic verdict matches enumeration.
        let mut rng = Lcg(FSCON_005_INPUT_SEED);
        let mut agree = 0;
        for _ in 0..12 {
            let rows: Vec<(f64, f64, f64)> = (0..4)
                .map(|_| {
                    (
                        rng.range(-1.0, 1.0),
                        rng.range(-1.0, 1.0),
                        rng.range(-1.5, 1.5),
                    )
                })
                .collect();
            let h = linear_host(&rows);
            let ss: Vec<ConstraintSpec> =
                (0..4).map(|i| hard(&format!("c{i}"), h.nodes[i])).collect();
            let d = diagnose_infeasibility(&h.problem, &ss, &domain, cx).expect("d");
            let truth = grid_feasible(&h.problem, &ss, &[0, 1, 2, 3], &domain);
            if d.feasible == truth {
                agree += 1;
            }
        }
        verdict(
            "fscon-005",
            core_right && full_infeasible && deletions_feasible && feasible_ok && agree == 12,
            &format!(
                "the triangle fixture yields the minimal core {{need-sum, cap-x, \
                 cap-y}} (bystander excluded); enumeration confirms the core is \
                 infeasible and EVERY single deletion restores feasibility; feasible \
                 systems return witnesses; elastic feasibility verdicts match \
                 enumeration on 12/12 seeded random fixtures ({agree}/12); \
                 input seed 0x1001_2026_0707_0005; fixed Cx execution seed \
                 0xC0C0"
            ),
            FSCON_005_INPUT_SEED,
        );
    });
}

/// Regression for frankensim-js9b: the constraints violated at the elastic
/// sum-optimum can be jointly feasible. The diagnosis must expand that seed
/// before deletion filtering rather than mislabeling the feasible support as
/// an unsat core.
#[test]
fn feasible_elastic_support_is_expanded_before_deletion_filtering() {
    with_cx(|cx| {
        // A: x >= 1, B: y >= 1, C: x + y <= 1. A and B are jointly
        // feasible at (1, 1). Scaling C by two makes the elastic sum attain
        // its minimum at (0.5, 0.5), where A and B are violated while C is
        // satisfied, so the raw support is the feasible set {A, B}.
        let host = linear_host(&[(-1.0, 0.0, -1.0), (0.0, -1.0, -1.0), (2.0, 2.0, 2.0)]);
        let specs: Vec<ConstraintSpec> = ["floor-x", "floor-y", "sum-cap"]
            .iter()
            .zip(&host.nodes)
            .map(|(name, &node)| hard(name, node))
            .collect();
        let domain = DomainBox {
            ranges: vec![(0.0, 1.0), (0.0, 1.0)],
        };
        let first = diagnose_infeasibility(&host.problem, &specs, &domain, cx).expect("first");
        let support: Vec<usize> = first
            .elastic
            .violations
            .iter()
            .enumerate()
            .filter(|&(_, &violation)| violation > 1e-6)
            .map(|(index, _)| index)
            .collect();
        assert_eq!(
            support,
            vec![0, 1],
            "the fixture must exercise the feasible elastic-support path"
        );
        assert!(
            grid_feasible(&host.problem, &specs, &support, &domain),
            "the deliberately feasible elastic support is the regression precondition"
        );
        let replay = diagnose_infeasibility(&host.problem, &specs, &domain, cx).expect("replay");
        let jointly_infeasible = !grid_feasible(&host.problem, &specs, &first.core, &domain);
        let deletions_feasible = first.core.iter().copied().all(|drop| {
            let rest: Vec<usize> = first
                .core
                .iter()
                .copied()
                .filter(|&index| index != drop)
                .collect();
            grid_feasible(&host.problem, &specs, &rest, &domain)
        });
        let deterministic = first.core == replay.core;

        verdict(
            "fscon-005-feasible-support",
            !first.feasible
                && first.core == vec![0, 1, 2]
                && jointly_infeasible
                && deletions_feasible
                && deterministic,
            &format!(
                "the feasible elastic support {{floor-x, floor-y}} expands to the deterministic minimal jointly-infeasible core {:?}; every deletion is feasible; fixed Cx execution seed 0xC0C0",
                first.core
            ),
            FIXED_INPUT_SEED,
        );
    });
}

/// fscon-006 — the worked repair example: ranked repairs whose
/// feasibility estimates are CALIBRATED against enumeration, and the
/// full diagnosis payload ships through fs-obs.
#[test]
fn fscon_006_repairs_calibrated() {
    with_cx(|cx| {
        // Mass budget vs strength floors: infeasible by 0.15 kg.
        let host = linear_host(&[
            (1.0, 1.0, 1.2),    // mass: x0 + x1 ≤ 1.2
            (-1.0, 0.0, -0.8),  // strength: x0 ≥ 0.8
            (0.0, -1.0, -0.55), // stiffness: x1 ≥ 0.55
        ]);
        let specs = vec![
            hard("mass-budget", host.nodes[0]),
            hard("strength-floor", host.nodes[1]),
            ConstraintSpec {
                name: "stiffness-pref".to_string(),
                node: host.nodes[2],
                kind: ConstraintKind::Soft(PenaltyLaw::Hinge { weight: 1.0 }),
                active_tol: 1e-9,
            },
        ];
        let domain = DomainBox {
            ranges: vec![(0.0, 2.0), (0.0, 2.0)],
        };
        let diag = diagnose_infeasibility(&host.problem, &specs, &domain, cx).expect("diag");
        let core_full = diag.core.len() == 3;
        let has_repairs = !diag.repairs.is_empty();
        // Ranking is by estimate, descending.
        let ranked = diag
            .repairs
            .windows(2)
            .all(|w| w[0].feasibility_estimate >= w[1].feasibility_estimate);
        // Soft members offer a drop action.
        let offers_drop = diag
            .repairs
            .iter()
            .any(|r| matches!(r.kind, fs_constraint::RepairKind::DropSoft { index: 2 }));
        // CALIBRATION: each estimate vs exact grid volume fraction.
        let mut worst_gap = 0.0f64;
        for r in &diag.repairs {
            let (relax, drop) = match r.kind {
                fs_constraint::RepairKind::RelaxBound { index, slack } => {
                    (vec![(index, slack)], None)
                }
                fs_constraint::RepairKind::DropSoft { index } => (vec![], Some(index)),
            };
            let n = 100;
            let mut hits = 0u32;
            for i in 0..=n {
                for j in 0..=n {
                    let x = vec![
                        2.0 * f64::from(i) / f64::from(n),
                        2.0 * f64::from(j) / f64::from(n),
                    ];
                    let ok = specs.iter().enumerate().all(|(k, s)| {
                        if Some(k) == drop {
                            return true;
                        }
                        let slack = relax
                            .iter()
                            .find(|(idx, _)| *idx == k)
                            .map_or(0.0, |(_, sl)| *sl);
                        fs_opt::eval(&host.problem, s.node, std::slice::from_ref(&x))
                            .expect("eval")
                            .scalar()
                            .expect("s")
                            <= slack
                    });
                    if ok {
                        hits += 1;
                    }
                }
            }
            let actual = f64::from(hits) / f64::from((n + 1) * (n + 1));
            worst_gap = worst_gap.max((r.feasibility_estimate - actual).abs());
        }
        let calibrated = worst_gap < 0.05;
        let payload = diag.to_json(&specs);
        assert!(
            !payload.contains("\"valid\":false"),
            "solver-produced diagnosis must satisfy its own evidence invariants: {payload}"
        );
        let mut em = fs_obs::Emitter::new("fs-constraint/conformance", "fscon-006/diagnosis");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "constraint-diagnosis".to_string(),
                    json: payload.clone(),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("diagnosis payload validates");
        println!("{line}");
        verdict(
            "fscon-006",
            core_full && has_repairs && ranked && offers_drop && calibrated,
            &format!(
                "the mass/strength/stiffness fixture diagnoses a 3-member core with \
                 RANKED repairs (drop-soft offered for the soft member); feasibility \
                 estimates are calibrated against exact enumeration (worst gap \
                 {worst_gap:.3} < 0.05); the full diagnosis payload ships through \
                 fs-obs under fixed Cx execution seed 0xC0C0: {payload}"
            ),
            FIXED_INPUT_SEED,
        );
    });
}

#[test]
fn elastic_domain_admission_refuses_malformed_ranges_before_solving() {
    with_cx(|cx| {
        let host = linear_host(&[(1.0, 0.0, 1.0)]);
        let specs = [hard("cap", host.nodes[0])];

        for (domain, expected_axis, expected_reason) in [
            (
                DomainBox {
                    ranges: vec![(f64::NAN, 1.0), (0.0, 1.0)],
                },
                0,
                DomainRangeError::NonFiniteEndpoint,
            ),
            (
                DomainBox {
                    ranges: vec![(0.0, 1.0), (0.0, f64::INFINITY)],
                },
                1,
                DomainRangeError::NonFiniteEndpoint,
            ),
            (
                DomainBox {
                    ranges: vec![(f64::NEG_INFINITY, 1.0), (0.0, 1.0)],
                },
                0,
                DomainRangeError::NonFiniteEndpoint,
            ),
            (
                DomainBox {
                    ranges: vec![(1.0, -1.0), (0.0, 1.0)],
                },
                0,
                DomainRangeError::Reversed,
            ),
            (
                DomainBox {
                    ranges: vec![(-f64::MAX, f64::MAX), (0.0, 1.0)],
                },
                0,
                DomainRangeError::UnrepresentableSpan,
            ),
        ] {
            for active_specs in [&specs[..], &[]] {
                assert!(matches!(
                    elastic_solve(&host.problem, active_specs, &domain, &[], cx),
                    Err(ConError::InvalidDomain(DomainError::InvalidRange {
                        axis,
                        reason,
                        ..
                    })) if axis == expected_axis && reason == expected_reason
                ));
            }
            assert!(matches!(
                elastic_solve(&host.problem, &specs, &domain, &[0], cx),
                Err(ConError::InvalidDomain(DomainError::InvalidRange { .. }))
            ));
        }

        for got in [1, 3] {
            let domain = DomainBox {
                ranges: vec![(0.0, 1.0); got],
            };
            assert_eq!(
                elastic_solve(&host.problem, &specs, &domain, &[], cx).unwrap_err(),
                ConError::InvalidDomain(DomainError::DimensionMismatch { expected: 2, got })
            );
        }

        let mut multi_builder = ProblemBuilder::new();
        let x = multi_builder
            .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("x");
        multi_builder
            .var("y", Manifold::Rn { dim: 1 }, Dims::NONE)
            .expect("y");
        let x_ref = multi_builder.var_ref(x).expect("x ref");
        let objective = multi_builder.norm_sq(x_ref).expect("objective");
        multi_builder
            .objective(objective, fs_opt::Sense::Minimize, 1.0)
            .expect("objective entry");
        let multi_problem = multi_builder.finish();
        assert_eq!(
            elastic_solve(
                &multi_problem,
                &[],
                &DomainBox {
                    ranges: vec![(0.0, 1.0)],
                },
                &[],
                cx,
            )
            .unwrap_err(),
            ConError::InvalidDomain(DomainError::HostVariableCount { got: 2 })
        );

        let mut sphere_builder = ProblemBuilder::new();
        let sphere = sphere_builder
            .var("sphere", Manifold::Sphere { ambient: 2 }, Dims::NONE)
            .expect("sphere");
        let sphere_ref = sphere_builder.var_ref(sphere).expect("sphere ref");
        let objective = sphere_builder.norm_sq(sphere_ref).expect("objective");
        sphere_builder
            .objective(objective, fs_opt::Sense::Minimize, 1.0)
            .expect("objective entry");
        let sphere_problem = sphere_builder.finish();
        assert_eq!(
            elastic_solve(
                &sphere_problem,
                &[],
                &DomainBox {
                    ranges: vec![(0.0, 1.0), (0.0, 1.0)],
                },
                &[],
                cx,
            )
            .unwrap_err(),
            ConError::InvalidDomain(DomainError::HostVariableManifold {
                got: Manifold::Sphere { ambient: 2 },
            })
        );

        let forged = [hard("forged", NodeId(u32::MAX))];
        let invalid = DomainBox {
            ranges: vec![(1.0, -1.0), (0.0, 1.0)],
        };
        assert!(matches!(
            elastic_solve(&host.problem, &forged, &invalid, &[], cx),
            Err(ConError::InvalidDomain(DomainError::InvalidRange { .. }))
        ));

        let fixed = DomainBox {
            ranges: vec![(0.0, 0.0), (1.0, 1.0)],
        };
        let report = elastic_solve(&host.problem, &specs, &fixed, &[], cx)
            .expect("zero-width axes are valid fixed coordinates");
        assert_eq!(report.x, vec![0.0, 1.0]);
        assert_eq!(
            report.total_violation,
            report.violations.iter().sum::<f64>(),
            "published components are the authority for the published total"
        );
        assert!(report.total_violation <= 1e-6);
    });
}

#[test]
fn json_surfaces_escape_untrusted_text_and_nonfinite_numbers() {
    let hostile = "name\"\\\n\r\t\u{0008}\u{000c}\u{0001}";
    let host = linear_host(&[(1.0, 0.0, 1.0)]);
    let spec = hard(hostile, host.nodes[0]);
    let mut evidence = evaluate(&host.problem, &spec, &[0.0, 0.0], None).expect("evidence");
    evidence.violation = f64::NAN;
    evidence.penalty = f64::INFINITY;
    let row = evidence.to_ledger_row();
    assert!(row.contains("\\\""));
    assert!(row.contains("\\\\"));
    assert!(row.contains("\\n"));
    assert!(row.contains("\\r"));
    assert!(row.contains("\\t"));
    assert!(row.contains("\\b"));
    assert!(row.contains("\\f"));
    assert!(row.contains("\\u0001"));
    assert!(row.contains("\"valid\":false"));
    assert!(row.contains("\"reason\":\"nonfinite-violation\""));
    assert!(row.contains("\"status\":\"no-claim\""));
    assert!(row.contains("\"violation\":null,\"penalty\":null"));
    assert!(!row.chars().any(|ch| ch <= '\u{001f}'));
    assert_eq!(row, evidence.to_ledger_row());

    let mut bad_penalty = evaluate(&host.problem, &spec, &[0.0, 0.0], None).expect("evidence");
    bad_penalty.penalty = f64::INFINITY;
    let bad_penalty_row = bad_penalty.to_ledger_row();
    assert!(bad_penalty_row.contains("\"status\":\"no-claim\""));
    assert!(bad_penalty_row.contains("\"reason\":\"nonfinite-penalty\""));

    let mut bad_certificate = evaluate(&host.problem, &spec, &[0.0, 0.0], None).expect("evidence");
    bad_certificate.certificate.lo = f64::NAN;
    let bad_certificate_row = bad_certificate.to_ledger_row();
    assert!(bad_certificate_row.contains("\"status\":\"no-claim\""));
    assert!(bad_certificate_row.contains("\"reason\":\"nonfinite-numerical-certificate\""));

    let mut wrong_certificate_kind =
        evaluate(&host.problem, &spec, &[0.0, 0.0], None).expect("evidence");
    wrong_certificate_kind.certificate.kind = fs_evidence::NumericalKind::Estimate;
    let wrong_certificate_row = wrong_certificate_kind.to_ledger_row();
    assert!(wrong_certificate_row.contains("\"status\":\"no-claim\""));
    assert!(wrong_certificate_row.contains("\"reason\":\"numerical-certificate-kind-mismatch\""));

    let diagnosis = Diagnosis {
        feasible: true,
        witness: Some(vec![0.0, 0.0]),
        core: vec![0, usize::MAX],
        repairs: vec![RepairAction {
            description: hostile.to_string(),
            kind: RepairKind::RelaxBound {
                index: 0,
                slack: 0.0,
            },
            feasibility_estimate: f64::NEG_INFINITY,
        }],
        elastic: ElasticReport {
            x: Vec::new(),
            total_violation: f64::NAN,
            violations: Vec::new(),
            evals: 0,
        },
    };
    let payload = diagnosis.to_json(std::slice::from_ref(&spec));
    assert!(payload.contains("\"valid\":false"));
    assert!(payload.contains("\"reason\":\"nonfinite-total-violation\""));
    assert!(payload.contains("\"feasible\":false"));
    assert!(payload.contains("\"total_violation\":null"));
    assert!(payload.contains("\"core\":[],\"repairs\":[]"));
    assert!(!payload.chars().any(|ch| ch <= '\u{001f}'));
    assert_eq!(payload, diagnosis.to_json(std::slice::from_ref(&spec)));

    let mut inconsistent_claim = diagnosis.clone();
    inconsistent_claim.elastic.total_violation = 0.0;
    inconsistent_claim.elastic.violations = vec![0.0];
    let inconsistent_payload = inconsistent_claim.to_json(std::slice::from_ref(&spec));
    assert!(inconsistent_payload.contains("\"valid\":false"));
    assert!(inconsistent_payload.contains("\"feasible\":false"));
    assert!(inconsistent_payload.contains("\"reason\":\"unknown-core-constraint\""));

    let valid_hostile = Diagnosis {
        feasible: false,
        witness: None,
        core: vec![0],
        repairs: vec![RepairAction {
            description: hostile.to_string(),
            kind: RepairKind::RelaxBound {
                index: 0,
                slack: 1.0,
            },
            feasibility_estimate: 0.5,
        }],
        elastic: ElasticReport {
            x: vec![0.0, 0.0],
            total_violation: 1.0,
            violations: vec![1.0],
            evals: 1,
        },
    };
    let valid_payload = valid_hostile.to_json(std::slice::from_ref(&spec));
    assert!(valid_payload.contains("\\\""));
    assert!(valid_payload.contains("\\\\"));
    assert!(valid_payload.contains("\\n"));
    assert!(!valid_payload.chars().any(|ch| ch <= '\u{001f}'));

    for (name, json) in [
        ("hostile-ledger", row),
        ("invalid-diagnosis", payload),
        ("hostile-diagnosis", valid_payload),
    ] {
        let mut emitter = fs_obs::Emitter::new("fs-constraint/conformance", name);
        let line = emitter
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: name.to_string(),
                    json,
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("hostile JSON payload remains valid");
    }
}

#[test]
fn malformed_policy_values_fail_closed_with_admitted_boundary_neighbors() {
    let host = linear_host(&[(1.0, 0.0, 1.0)]);
    let x = [0.0, 0.0];
    let min_subnormal = f64::from_bits(1);

    for active_tol in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY, -min_subnormal] {
        let mut spec = hard("tol", host.nodes[0]);
        spec.active_tol = active_tol;
        assert!(matches!(
            evaluate(&host.problem, &spec, &x, None),
            Err(ConError::BadParam {
                what: "active-set tolerance",
                ..
            })
        ));
    }
    for active_tol in [-0.0, 0.0, min_subnormal] {
        let mut spec = hard("tol-neighbor", host.nodes[0]);
        spec.active_tol = active_tol;
        assert!(evaluate(&host.problem, &spec, &x, None).is_ok());
    }

    for weight in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY, -min_subnormal] {
        let spec = ConstraintSpec {
            name: "soft".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Soft(PenaltyLaw::Hinge { weight }),
            active_tol: 0.0,
        };
        assert!(matches!(
            evaluate(&host.problem, &spec, &[2.0, 0.0], None),
            Err(ConError::BadParam {
                what: "soft-constraint weight",
                ..
            })
        ));
    }
    for weight in [-0.0, 0.0, min_subnormal, 1.0] {
        let spec = ConstraintSpec {
            name: "soft-neighbor".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Soft(PenaltyLaw::Quadratic { weight }),
            active_tol: 0.0,
        };
        let evidence = evaluate(&host.problem, &spec, &[2.0, 0.0], None).expect("valid weight");
        assert!(evidence.penalty.is_finite() && evidence.penalty >= 0.0);
    }
    let overflowing_penalty = ConstraintSpec {
        name: "overflow".into(),
        node: host.nodes[0],
        kind: ConstraintKind::Soft(PenaltyLaw::Quadratic { weight: f64::MAX }),
        active_tol: 0.0,
    };
    assert!(matches!(
        evaluate(&host.problem, &overflowing_penalty, &[3.0, 0.0], None),
        Err(ConError::BadParam {
            what: "soft-constraint penalty result",
            ..
        })
    ));

    for half_widths in [
        vec![-min_subnormal, 0.0],
        vec![f64::NAN, 0.0],
        vec![f64::INFINITY, 0.0],
        vec![0.0],
        Vec::new(),
    ] {
        let spec = ConstraintSpec {
            name: "robust".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Robust { half_widths },
            active_tol: 0.0,
        };
        assert!(matches!(
            evaluate(&host.problem, &spec, &x, None),
            Err(ConError::BadParam { .. })
        ));
    }
    let span_overflow = ConstraintSpec {
        name: "robust-overflow".into(),
        node: host.nodes[0],
        kind: ConstraintKind::Robust {
            half_widths: vec![f64::MAX, 0.0],
        },
        active_tol: 0.0,
    };
    assert!(matches!(
        evaluate(&host.problem, &span_overflow, &x, None),
        Err(ConError::BadParam {
            what: "robust interval range",
            ..
        })
    ));
    for half_width in [-0.0, 0.0, min_subnormal, f64::MAX / 2.0] {
        let spec = ConstraintSpec {
            name: "robust-neighbor".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Robust {
                half_widths: vec![half_width, 0.0],
            },
            active_tol: 0.0,
        };
        assert!(evaluate(&host.problem, &spec, &x, None).is_ok());
    }

    let confidence_rounding_tie = f64::EPSILON / 4.0;
    for delta in [min_subnormal, confidence_rounding_tie] {
        let unrepresentable_delta = ConstraintSpec {
            name: "unrepresentable-delta".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Chance {
                level: 0.5,
                estimator: ChanceEstimator::MonteCarlo { samples: 1, delta },
            },
            active_tol: 0.0,
        };
        assert!(matches!(
            evaluate(&host.problem, &unrepresentable_delta, &x, None),
            Err(ConError::BadParam {
                what: "chance confidence representation",
                ..
            })
        ));
    }
    let first_delta_above_tie = f64::from_bits(confidence_rounding_tie.to_bits() + 1);
    let representable_delta = ConstraintSpec {
        name: "representable-delta".into(),
        node: host.nodes[0],
        kind: ConstraintKind::Chance {
            level: 0.5,
            estimator: ChanceEstimator::MonteCarlo {
                samples: 1,
                delta: first_delta_above_tie,
            },
        },
        active_tol: 0.0,
    };
    let zero_noise = |_sample: u64| vec![0.0, 0.0];
    let representable_delta_evidence =
        evaluate(&host.problem, &representable_delta, &x, Some(&zero_noise))
            .expect("a confidence-distinguishable delta is admitted");
    assert!(matches!(
        representable_delta_evidence.statistical,
        fs_evidence::StatisticalCertificate::HalfWidth { half_width, .. }
            if half_width.is_finite()
    ));
    let mut wrong_statistical_kind = representable_delta_evidence.clone();
    wrong_statistical_kind.statistical = fs_evidence::StatisticalCertificate::EValue {
        e: 1.0,
        alpha: 0.05,
    };
    assert!(
        wrong_statistical_kind
            .to_ledger_row()
            .contains("\"reason\":\"chance-statistical-certificate-kind-mismatch\"")
    );
    let short_noise = |_sample: u64| vec![0.0];
    let long_noise = |_sample: u64| vec![0.0, 0.0, 0.0];
    for bad_noise in [
        &short_noise as &dyn Fn(u64) -> Vec<f64>,
        &long_noise as &dyn Fn(u64) -> Vec<f64>,
    ] {
        assert!(matches!(
            evaluate(&host.problem, &representable_delta, &x, Some(bad_noise)),
            Err(ConError::BadParam {
                what: "chance noise draw dimension",
                ..
            })
        ));
    }
}

#[test]
fn interval_proofs_refuse_bad_boxes_and_wrong_proof_kinds_and_bind_subjects() {
    let host = linear_host(&[(1.0, 0.0, 1.0), (0.0, 1.0, 1.0)]);
    let interval_spec = ConstraintSpec {
        name: "interval".into(),
        node: host.nodes[0],
        kind: ConstraintKind::Certification {
            proof: ProofKind::Interval,
        },
        active_tol: 0.0,
    };
    for domain in [
        [(1.0, 0.0), (0.0, 1.0)],
        [(f64::NAN, 1.0), (0.0, 1.0)],
        [(0.0, f64::INFINITY), (0.0, 1.0)],
        [(-f64::MAX, f64::MAX), (0.0, 1.0)],
    ] {
        assert!(matches!(
            interval_eval(&host.problem, interval_spec.node, &domain),
            Err(fs_constraint::IvalError::BadBindings)
        ));
        assert!(matches!(
            prove_interval(&host.problem, &interval_spec, &domain),
            Err(ConError::NotProvable { .. })
        ));
    }
    for domain in [[(0.0, 0.0), (1.0, 1.0)], [(-f64::MAX, 0.0), (0.0, 1.0)]] {
        assert!(interval_eval(&host.problem, interval_spec.node, &domain).is_ok());
    }
    for domain in [vec![(0.0, 1.0)], vec![(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]] {
        assert!(matches!(
            interval_eval(&host.problem, interval_spec.node, &domain),
            Err(fs_constraint::IvalError::BadBindings)
        ));
        assert!(matches!(
            prove_interval(&host.problem, &interval_spec, &domain),
            Err(ConError::NotProvable { .. })
        ));
    }

    let mut multi_builder = ProblemBuilder::new();
    let multi_x = multi_builder
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("x");
    multi_builder
        .var("y", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("y");
    let multi_ref = multi_builder.var_ref(multi_x).expect("x ref");
    let multi_node = multi_builder.norm_sq(multi_ref).expect("objective");
    multi_builder
        .objective(multi_node, fs_opt::Sense::Minimize, 1.0)
        .expect("objective");
    let multi_problem = multi_builder.finish();
    assert!(matches!(
        interval_eval(&multi_problem, multi_node, &[(0.0, 1.0)]),
        Err(fs_constraint::IvalError::BadBindings)
    ));

    let mut constant_builder = ProblemBuilder::new();
    let constant_node = constant_builder.konst(-1.0, Dims::NONE).expect("constant");
    constant_builder
        .objective(constant_node, fs_opt::Sense::Minimize, 1.0)
        .expect("objective");
    let constant_problem = constant_builder.finish();
    let constant_spec = ConstraintSpec {
        name: "constant-proof".into(),
        node: constant_node,
        kind: ConstraintKind::Certification {
            proof: ProofKind::Interval,
        },
        active_tol: 0.0,
    };
    let (constant_evidence, constant_artifact) =
        prove_interval(&constant_problem, &constant_spec, &[]).expect("constant proof");
    assert!(constant_artifact.verifies_evidence(&constant_evidence));
    assert!(constant_artifact.is_bound_to(&constant_problem, constant_node, &[]));
    assert!(matches!(
        interval_eval(&constant_problem, constant_node, &[(0.0, 0.0)]),
        Err(fs_constraint::IvalError::BadBindings)
    ));

    let mut sphere_builder = ProblemBuilder::new();
    let sphere = sphere_builder
        .var("sphere", Manifold::Sphere { ambient: 2 }, Dims::NONE)
        .expect("sphere");
    let sphere_ref = sphere_builder.var_ref(sphere).expect("sphere ref");
    let sphere_node = sphere_builder.norm_sq(sphere_ref).expect("objective");
    sphere_builder
        .objective(sphere_node, fs_opt::Sense::Minimize, 1.0)
        .expect("objective");
    let sphere_problem = sphere_builder.finish();
    assert!(matches!(
        interval_eval(&sphere_problem, sphere_node, &[(0.0, 1.0), (0.0, 1.0)]),
        Err(fs_constraint::IvalError::BadBindings)
    ));

    let sos = ConstraintSpec {
        kind: ConstraintKind::Certification {
            proof: ProofKind::Sos,
        },
        ..interval_spec.clone()
    };
    let sos_intent = evaluate(&host.problem, &sos, &[0.0, 0.0], None)
        .expect("SOS intent remains representable without a forged artifact");
    assert_eq!(
        sos_intent.status,
        Status::NeedsProof {
            proof: ProofKind::Sos
        }
    );
    assert!(
        sos_intent
            .to_ledger_row()
            .contains("\"status\":\"needs-proof:Sos\"")
    );
    assert_eq!(
        prove_interval(&host.problem, &sos, &[(0.0, 1.0), (0.0, 1.0)]).unwrap_err(),
        ConError::ProofKindMismatch {
            requested: Some(ProofKind::Sos),
            attempted: ProofKind::Interval,
        }
    );
    let hard_spec = hard("not-certification", host.nodes[0]);
    assert_eq!(
        prove_interval(&host.problem, &hard_spec, &[(0.0, 1.0), (0.0, 1.0)],).unwrap_err(),
        ConError::ProofKindMismatch {
            requested: None,
            attempted: ProofKind::Interval,
        }
    );

    let domain = [(0.0, 0.5), (0.0, 0.5)];
    let (mut proven_evidence, artifact) =
        prove_interval(&host.problem, &interval_spec, &domain).expect("admitted proof");
    assert!(artifact.verifies_evidence(&proven_evidence));
    assert_eq!(
        proven_evidence.proof_bound(),
        Some(artifact.interval_bound())
    );
    let proven_row = proven_evidence.to_ledger_row();
    assert!(proven_row.contains("\"status\":\"proven\""));
    assert!(proven_row.contains("\"proof_subject\":"));
    assert!(proven_row.contains("\"proof_bound_bits\":"));
    assert!(proven_row.contains(&artifact.subject().problem().to_hex()));
    let neighbor_domain = [(0.0, 0.25), (0.0, 0.5)];
    let (neighbor_evidence, neighbor_artifact) =
        prove_interval(&host.problem, &interval_spec, &neighbor_domain).expect("neighbor proof");
    assert!(neighbor_artifact.verifies_evidence(&neighbor_evidence));
    assert!(!artifact.verifies_evidence(&neighbor_evidence));
    assert!(!neighbor_artifact.verifies_evidence(&proven_evidence));
    proven_evidence.status = Status::Satisfied;
    let downgraded = proven_evidence.to_ledger_row();
    assert!(downgraded.contains("\"status\":\"no-claim\""));
    assert!(downgraded.contains("\"reason\":\"unexpected-proof-evidence\""));
    assert!(!artifact.verifies_evidence(&proven_evidence));
    assert!(artifact.is_bound_to(&host.problem, interval_spec.node, &domain));
    assert!(!artifact.is_bound_to(
        &host.problem,
        interval_spec.node,
        &[(0.0, 0.5), (0.0, 0.75)],
    ));
    assert!(!artifact.is_bound_to(&host.problem, host.nodes[1], &domain));
    let signed_zero_domain = [(-0.0, 0.5), (0.0, 0.5)];
    assert!(!artifact.is_bound_to(&host.problem, interval_spec.node, &signed_zero_domain,));
    assert_eq!(artifact.subject().node(), interval_spec.node);
    assert_eq!(artifact.subject().domain_bits()[0].0, 0.0f64.to_bits());
}

#[test]
fn fscon_v2_wire_encoding_is_injective_and_raw_policies_are_admitted() {
    let host = linear_host(&[(1.0, 0.0, 1.0)]);
    let min_subnormal = f64::from_bits(1);
    let confidence_rounding_tie = f64::EPSILON / 4.0;
    let first_delta_above_tie = f64::from_bits(confidence_rounding_tie.to_bits() + 1);
    let specs = vec![
        ConstraintSpec {
            name: "literal%20 versus space \n tab\t and 🧬".into(),
            node: host.nodes[0],
            kind: ConstraintKind::Fabrication {
                process: "%20 process name".into(),
            },
            active_tol: -0.0,
        },
        ConstraintSpec {
            name: String::new(),
            node: host.nodes[0],
            kind: ConstraintKind::Code {
                standard: String::new(),
            },
            active_tol: f64::MIN_POSITIVE,
        },
    ];
    let encoded = serialize_specs(&specs);
    assert!(encoded.starts_with("fscon v2\n"));
    assert!(encoded.contains("literal%2520%20versus%20space%20%0A%20tab%09"));
    assert_eq!(parse_specs(&encoded).expect("canonical round trip"), specs);
    assert_eq!(serialize_specs(&parse_specs(&encoded).unwrap()), encoded);

    let mut literal = hard("%20", host.nodes[0]);
    literal.active_tol = 0.0;
    let mut space = literal.clone();
    space.name = " ".into();
    assert_ne!(serialize_specs(&[literal]), serialize_specs(&[space]));
    let empty = hard("", host.nodes[0]);
    let percent = hard("%", host.nodes[0]);
    assert_ne!(serialize_specs(&[empty]), serialize_specs(&[percent]));

    let raw = |tol: f64, kind: &str| {
        format!("fscon v2\nconstraint bad 0 {:016X} {kind}\n", tol.to_bits())
    };
    for text in [
        raw(f64::NAN, "hard"),
        raw(f64::INFINITY, "hard"),
        raw(-min_subnormal, "hard"),
        raw(0.0, &format!("soft hinge {:016X}", f64::NAN.to_bits())),
        raw(0.0, &format!("soft quadratic {:016X}", (-1.0f64).to_bits())),
        raw(
            0.0,
            &format!("robust {:016X},{:016X}", f64::INFINITY.to_bits(), 0),
        ),
        raw(
            0.0,
            &format!(
                "chance {:016X} mc 4 {:016X}",
                0.9f64.to_bits(),
                f64::NAN.to_bits()
            ),
        ),
        raw(
            0.0,
            &format!(
                "chance {:016X} mc 4 {:016X}",
                f64::INFINITY.to_bits(),
                0.05f64.to_bits()
            ),
        ),
        raw(
            0.0,
            &format!(
                "chance {:016X} mc 1 {:016X}",
                0.5f64.to_bits(),
                min_subnormal.to_bits()
            ),
        ),
        raw(
            0.0,
            &format!(
                "chance {:016X} mc 1 {:016X}",
                0.5f64.to_bits(),
                confidence_rounding_tie.to_bits()
            ),
        ),
    ] {
        assert!(matches!(parse_specs(&text), Err(ConError::Parse { .. })));
    }
    for text in [
        raw(-0.0, "hard"),
        raw(0.0, &format!("soft hinge {:016X}", (-0.0f64).to_bits())),
        raw(
            0.0,
            &format!("robust {:016X},{:016X}", (-0.0f64).to_bits(), 0),
        ),
        raw(
            0.0,
            &format!(
                "chance {:016X} mc 1 {:016X}",
                0.5f64.to_bits(),
                first_delta_above_tie.to_bits()
            ),
        ),
    ] {
        assert!(
            parse_specs(&text).is_ok(),
            "boundary neighbor must admit: {text}"
        );
    }
    assert!(parse_specs("fscon v2\nconstraint %2f 0 0000000000000000 hard\n").is_err());
    assert!(parse_specs("fscon v2\nconstraint %41 0 0000000000000000 hard\n").is_err());
    assert!(parse_specs("fscon v1\n").is_err());
}

#[test]
fn a_pre_cancelled_elastic_solve_refuses_before_expression_evaluation() {
    let host = linear_host(&[(1.0, 0.0, 1.0)]);
    let forged = [hard("must-not-evaluate", NodeId(u32::MAX))];
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(matches!(
            elastic_solve(
                &host.problem,
                &forged,
                &DomainBox {
                    ranges: vec![(0.0, 1.0), (0.0, 1.0)],
                },
                &[],
                &cx,
            ),
            Err(ConError::Eval(fs_opt::OptError::Cancelled))
        ));
    });
}

const _: () = {
    // Compile-time reminder that Diagnosis is the agent-facing artifact.
    fn _assert_payload(d: &Diagnosis, s: &[ConstraintSpec]) -> String {
        d.to_json(s)
    }
};

#[test]
fn a_non_finite_constraint_value_is_never_certified_feasible() {
    // Regression: a design point OUTSIDE a constraint's domain (here `sqrt` of a
    // negative argument -> NaN) must be maximally VIOLATED, never Satisfied.
    // Every IEEE comparison with NaN is false, so the old status ladder fell
    // through to `Satisfied` and `NaN.max(0.0) == 0.0` attached an EXACT
    // zero-violation certificate -- certifying an undefined constraint as
    // strictly feasible (a false certificate).
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 1 }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("ref");
    let x0 = b.component(vr, 0).expect("x0");
    let g = b.sqrt(x0).expect("sqrt"); // g = sqrt(x0): NaN for x0 < 0
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    let problem = b.finish();
    let spec = hard("domain", g);

    // Out-of-domain point (x0 = -0.5 -> sqrt = NaN): the explicit-stack
    // evaluator (bead xf8v7) refuses with a typed non-finite error naming
    // the exact node — strictly stronger than the previous
    // Violated-with-infinite-violation contract, and still never
    // certified feasible.
    let nan_err = evaluate(&problem, &spec, &[-0.5], None)
        .expect_err("a non-finite constraint value must refuse evaluation");
    assert!(
        matches!(
            nan_err,
            ConError::Eval(fs_opt::OptError::EvalNonFinite { .. })
        ),
        "expected a typed EvalNonFinite refusal, got {nan_err:?}"
    );
    // A finite in-domain point still classifies with a FINITE violation.
    let finite_ev = evaluate(&problem, &spec, &[4.0], None).expect("evaluate returns Ok");
    assert!(
        finite_ev.violation.is_finite(),
        "a finite constraint value must yield a finite violation, got {}",
        finite_ev.violation
    );
    verdict(
        "fscon-nonfinite",
        matches!(
            nan_err,
            ConError::Eval(fs_opt::OptError::EvalNonFinite { .. })
        ) && finite_ev.violation.is_finite(),
        "a NaN (out-of-domain) constraint refuses with a typed EvalNonFinite, never a feasible \
         exact-0 certificate",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn a_chance_constraint_with_a_bad_delta_or_zero_samples_is_refused() {
    // Regression: the Hoeffding half-width sqrt(ln(1/delta)/(2n)) is NaN for
    // delta >= 1 (ln(1/delta) <= 0) and +inf for n = 0, and confidence = 1-delta
    // falls outside [0,1] for delta outside (0,1). Unvalidated, these produced a
    // garbage certificate; they must be refused as BadParam, like the level.
    let host = linear_host(&[(1.0, 0.0, 1.0)]);
    let noise = |_s: u64| -> Vec<f64> { vec![0.5, 0.0] };
    let x = [0.2, 0.0];
    let chance = |samples: u32, delta: f64| ConstraintSpec {
        name: "chance".to_string(),
        node: host.nodes[0],
        kind: ConstraintKind::Chance {
            level: 0.9,
            estimator: ChanceEstimator::MonteCarlo { samples, delta },
        },
        active_tol: 1e-9,
    };
    for (s, d, why) in [
        (400u32, 1.5f64, "delta >= 1"),
        (400, 0.0, "delta = 0"),
        (0, 0.05, "zero samples"),
    ] {
        assert!(
            matches!(
                evaluate(&host.problem, &chance(s, d), &x, Some(&noise)),
                Err(ConError::BadParam { .. })
            ),
            "{why} must be refused as BadParam"
        );
    }
    // A valid (delta, samples) still evaluates.
    assert!(evaluate(&host.problem, &chance(400, 0.05), &x, Some(&noise)).is_ok());
    verdict(
        "fscon-chance-params",
        true,
        "invalid chance delta / zero samples are refused, not turned into a NaN certificate",
        FIXED_INPUT_SEED,
    );
}

/// Forged/stale NodeIds are typed refusals, never index panics
/// (batch-verify High #2): interval evaluation checks the arena
/// boundary before touching any expression.
#[test]
fn forged_node_ids_refuse_instead_of_panicking() {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 1 }, fs_qty::Dims::NONE)
        .expect("var");
    let r = b.var_ref(v).expect("ref");
    // Objectives are scalar-only under the sealed-IR rules; project the
    // 1-dim vector ref down before declaring it.
    let s = b.component(r, 0).expect("scalar component");
    b.objective(s, fs_opt::Sense::Minimize, 1.0).expect("o");
    let small = b.finish();
    let boxes = [(0.0, 1.0)];
    for forged in [NodeId(u32::MAX), NodeId(1_000)] {
        let err =
            interval_eval(&small, forged, &boxes).expect_err("out-of-arena node id must refuse");
        assert!(
            matches!(
                err,
                fs_constraint::IvalError::UnknownNode { node } if node == forged.0
            ),
            "typed UnknownNode refusal expected, got {err:?}"
        );
    }
}

/// Interval work follows the admitted DAG rather than its exponentially
/// large tree expansion, and integer powers take logarithmic work even at
/// the public i32 boundary.
#[test]
fn interval_eval_bounds_shared_dag_and_powi_work() {
    let mut b = ProblemBuilder::new();
    let one = b.konst(1.0, Dims::NONE).expect("one");
    let mut doubled = one;
    for _ in 0..40 {
        doubled = b.add(doubled, doubled).expect("shared DAG level");
    }
    let huge_power = b
        .powi(one, i32::MAX)
        .expect("positive exponent is admitted");
    let root = b.add(doubled, huge_power).expect("root");
    b.objective(root, fs_opt::Sense::Minimize, 1.0)
        .expect("objective");
    let problem = b.finish();

    let enclosure = interval_eval(&problem, root, &[]).expect("bounded interval work");
    let expected = 2.0f64.powi(40) + 1.0;
    assert_eq!(enclosure.lo.to_bits(), expected.to_bits());
    assert_eq!(enclosure.hi.to_bits(), expected.to_bits());
}

/// G4 depth boundary: interval recursion is tied to fs-opt's admission
/// schedule. The exact cap evaluates without overflowing the stack; a
/// graph built under a deliberately looser policy refuses before
/// recursive interval work begins.
#[test]
fn interval_eval_respects_the_admitted_depth_boundary() {
    fn chain(caps: AdmissionCaps, depth: u32) -> (Problem, NodeId) {
        let mut builder = ProblemBuilder::with_caps(caps);
        let mut root = builder.konst(1.0, Dims::NONE).expect("depth 1");
        for _ in 1..depth {
            root = builder.neg(root).expect("next depth");
        }
        builder
            .objective(root, fs_opt::Sense::Minimize, 1.0)
            .expect("objective");
        (builder.finish(), root)
    }

    let defaults = AdmissionCaps::default();
    let limit = defaults.max_graph_depth;
    let (at_limit, root) = chain(defaults.clone(), limit);
    let interval = interval_eval(&at_limit, root, &[]).expect("exact depth cap evaluates");
    let expected: f64 = if limit % 2 == 0 { -1.0 } else { 1.0 };
    assert_eq!(interval.lo.to_bits(), expected.to_bits());
    assert_eq!(interval.hi.to_bits(), expected.to_bits());

    let mut relaxed = defaults;
    relaxed.max_graph_depth = limit + 1;
    let (over_limit, over_root) = chain(relaxed, limit + 1);
    assert!(matches!(
        interval_eval(&over_limit, over_root, &[]),
        Err(fs_constraint::IvalError::CapExceeded {
            what: "graph depth",
            count,
            cap,
        }) if count == u64::from(limit) + 1 && cap == u64::from(limit)
    ));
}
