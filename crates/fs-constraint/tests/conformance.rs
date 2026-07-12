//! fs-constraint conformance suite (CONTRACT.md: any reimplementation
//! must pass). Kind taxonomy + treatments + serialization, evidence
//! statuses/roles/penalties, chance validity machinery (the BOUND
//! decides, not the raw rate), certification refusals + real interval
//! proofs, minimal unsat cores vs enumeration, and the worked repair
//! example with calibrated feasibility estimates. JSON-line verdicts;
//! seeded cases carry seeds.

use asupersync::types::Budget;
use fs_constraint::{
    ChanceEstimator, ConError, ConstraintKind, ConstraintSpec, Diagnosis, DomainBox, PenaltyLaw,
    ProofKind, Status, Treatment, diagnose_infeasibility, evaluate, interval_eval, parse_specs,
    prove_interval, serialize_specs,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder};
use fs_qty::Dims;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-constraint/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
                seed: 0xC0C0,
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
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE);
    let vr = b.var_ref(v).expect("ref");
    let x0 = b.component(vr, 0).expect("x0");
    let x1 = b.component(vr, 1).expect("x1");
    let mut nodes = Vec::new();
    for &(a0, a1, rhs) in rows {
        let c0 = b.konst(a0, Dims::NONE);
        let c1 = b.konst(a1, Dims::NONE);
        let t0 = b.mul(c0, x0).expect("t0");
        let t1 = b.mul(c1, x1).expect("t1");
        let s = b.add(t0, t1).expect("s");
        let rb = b.konst(rhs, Dims::NONE);
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
    let garbage = parse_specs("fscon v1\nwat\n");
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
        let mut r = Lcg(0x1001_2026_0707_0003 ^ (s.wrapping_mul(0x9E37_79B9_7F4A_7C15)));
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
             StatisticalCertificate; seed 0x1001_2026_0707_0003",
            squeezed = squeezed.status
        ),
    );
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
        && matches!(artifact, fs_constraint::ProofArtifact::IntervalBound { hi } if hi <= 0.0);
    // Unprovable domain: [0,2]² has hi = 1 > 0 → honest refusal.
    let refused = prove_interval(&host.problem, &cert, &[(0.0, 2.0), (0.0, 2.0)]);
    let honest = matches!(refused, Err(ConError::NotProvable { ref why }) if why.contains("upper"));
    // G0 containment: random nonlinear graph, random boxes, samples in.
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE);
    let vr = b.var_ref(v).expect("r");
    let n = b.norm_sq(vr).expect("n");
    let x0 = b.component(vr, 0).expect("x0");
    let t = b.tanh(x0).expect("t");
    let m = b.mul(n, n).expect("m");
    let s = b.add(m, t).expect("s");
    let a = b.abs(s).expect("a");
    let one = b.konst(1.0, Dims::NONE);
    let g = b.min_of(a, one).expect("g");
    b.objective(g, fs_opt::Sense::Minimize, 1.0).expect("o");
    let nl = b.finish();
    let mut rng = Lcg(0x1001_2026_0707_0004);
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
        let mut rng = Lcg(0x1001_2026_0707_0005);
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
                 seed 0x1001_2026_0707_0005"
            ),
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
                 fs-obs: {payload}"
            ),
        );
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
    let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE);
    let vr = b.var_ref(v).expect("ref");
    let x0 = b.component(vr, 0).expect("x0");
    let g = b.sqrt(x0).expect("sqrt"); // g = sqrt(x0): NaN for x0 < 0
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    let problem = b.finish();
    let spec = hard("domain", g);

    // Out-of-domain point (x0 = -0.5 -> sqrt = NaN): must be Violated, never
    // Satisfied, with a positive (infinite) violation.
    let nan_ev = evaluate(&problem, &spec, &[-0.5], None).expect("evaluate returns Ok");
    assert!(
        matches!(nan_ev.status, Status::Violated),
        "a non-finite constraint must be Violated, got {:?}",
        nan_ev.status
    );
    assert!(
        nan_ev.violation > 0.0,
        "a non-finite constraint must report a positive violation, got {}",
        nan_ev.violation
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
        matches!(nan_ev.status, Status::Violated) && finite_ev.violation.is_finite(),
        "a NaN (out-of-domain) constraint is Violated, not a feasible exact-0 certificate",
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
    );
}
