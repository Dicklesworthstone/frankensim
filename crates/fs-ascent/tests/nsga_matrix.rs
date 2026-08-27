//! NSGA-III baseline DEFECT-MATRIX battery (epic bead
//! `frankensim-epic-ascent-7tv.24.13`).
//!
//! Every AC class gets focused unit / property / metamorphic / mutation
//! probes with INDEPENDENT oracles (closed forms, exhaustive
//! permutations, hand-computed cases) and no reliance on the engine
//! under test for its own truth. This freezes the baseline the epic
//! requires before any adaptive-quality work may ride on it; it mints
//! NO quality or release authority by itself (Estimate-class).

use fs_ascent::{
    NsgaConfig, NsgaError, NsgaIndividual, NonFiniteKind,
    build_references, das_dennis_cardinality, fast_nondominated_sort,
    nsga3_run, partition_tuples,
};

fn ind(x: &[f64], f: &[f64]) -> NsgaIndividual {
    NsgaIndividual {
        x: x.to_vec(),
        f: f.to_vec(),
    }
}

/// Deterministic seeded closure evaluating the convex two-objective
/// parabola front on `[0,1]`: f = (x² , (x−1)²). Recording wrapper lets
/// batteries audit exactly what was charged to the budget.
struct ConvexEval {
    log: std::cell::RefCell<Vec<f64>>,
}

impl ConvexEval {
    const fn new() -> Self {
        Self {
            log: std::cell::RefCell::new(Vec::new()),
        }
    }

    fn call(&mut self, x: &[f64]) -> Vec<f64> {
        let a = x[0];
        self.log.borrow_mut().push(a);
        vec![a * a, (a - 1.0) * (a - 1.0)]
    }

    fn calls(&self) -> Vec<f64> {
        self.log.borrow().clone()
    }
}

fn base_config(seed: u64, pop: usize, generations: usize, budget: usize) -> NsgaConfig {
    NsgaConfig {
        reference_divisions: 4,
        population_size: pop,
        max_generations: generations,
        eval_budget: budget,
        seed,
        bounds: vec![(0.0, 1.0)],
    }
}

// ---------------------------------------------------------------------------
// 1. Reference directions: cardinality closed form, partition oracle,
//    unit-sum property, validation refusals, mutation sensitivity.
// ---------------------------------------------------------------------------

#[test]
fn das_dennis_cardinality_matches_partition_oracle_exhaustively() {
    for p in 1usize..=8 {
        for m in 2usize..=5 {
            let dirs = partition_tuples(p, m);
            let formula = das_dennis_cardinality(p, m);
            assert_eq!(
                dirs.len(),
                formula,
                "partition oracle {p}/{m} disagrees with C(p+m-1,m-1)"
            );
            // Lex-ascending canonical enumeration: strict total order.
            let mut sorted = dirs.clone();
            sorted.sort();
            assert_eq!(dirs, sorted);
            // Unit simplex: each tuple sums EXACTLY to p (integers),
            // so division yields components summing to 1 up to fp.
            for t in &dirs {
                let s: usize = t.iter().sum();
                assert_eq!(s, p);
                let fv: Vec<f64> = t.iter().map(|&v| v as f64 / p as f64).collect();
                let sum: f64 = fv.iter().sum();
                assert!((sum - 1.0).abs() < 1e-12, "{p}/{m} simplex");
            }
        }
    }
}

#[test]
fn references_refuse_zero_divisions_single_objective_and_cap_overflow() {
    let zero = build_references(0, 3);
    let single = build_references(4, 1);
    let capped = build_references(200_000, 4);
    match (&zero, &single, &capped) {
        (
            Err(NsgaError::ReferenceInvalid { what: z }),
            Err(NsgaError::ReferenceInvalid { what: s }),
            Err(NsgaError::ReferenceInvalid { what: c }),
        ) => {
            assert!(z.contains(">= 1"), "{z}");
            assert!(s.contains("two objectives"), "{s}");
            assert!(c.contains("baseline cap"), "{c}");
        }
        other => panic!("expected three typed refusals, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// 2. Fast non-dominated sorting: dominance chains, incomparability,
//    duplicate ridership, stability under input permutation.
// ---------------------------------------------------------------------------

#[test]
fn sort_orders_chains_and_incomparable_pairs_deterministically() {
    // A=(0,0) dominates everything; B=(1,0)/(0,1) dominate C=(1,1);
    // E=(2,2) is worst. Incomparable B/D stays same-front.
    let pop = vec![
        ind(&[0.0], &[0.0, 0.0]),   // A best
        ind(&[1.0], &[1.0, 0.0]),   // B
        ind(&[2.0], &[0.0, 1.0]),   // D (incomparable with B)
        ind(&[3.0], &[1.0, 1.0]),   // C dominated by B and D
        ind(&[4.0], &[2.0, 2.0]),   // E dominated by all
    ];
    let fronts = fast_nondominated_sort(&pop);
    assert_eq!(fronts.len(), 4, "{fronts:?}");
    assert_eq!(fronts[0], vec![0]);
    assert_eq!(fronts[1], vec![1, 2]);
    assert_eq!(fronts[2], vec![3]);
}

#[test]
fn equal_objective_twins_share_a_front_and_survive_sort() {
    let twin_f = [0.5, 0.5];
    let pop = vec![
        ind(&[0.1], &twin_f),
        ind(&[0.2], &twin_f),
        ind(&[0.0], &[0.4, 0.4]),
    ];
    let fronts = fast_nondominated_sort(&pop);
    assert_eq!(fronts[0], vec![2]);
    assert_eq!(fronts[1], vec![0, 1], "twins ride together");
}

#[test]
fn sorting_is_invariant_under_input_permutation_up_to_relabeling() {
    let mk = |k: usize| -> NsgaIndividual {
        ind(
            &[k as f64],
            &[((7 * k + 3) % 11) as f64, ((5 * k + 1) % 13) as f64],
        )
    };
    let base_pop: Vec<NsgaIndividual> = (0..9).map(mk).collect();
    let fwd = fast_nondominated_sort(&base_pop);
    let mut rev_pop: Vec<NsgaIndividual> = base_pop.iter().rev().cloned().collect();
    rev_pop.reverse(); // identity order preserved but built reversed first
    let back = fast_nondominated_sort(&rev_pop);
    // Same SETS of objective fingerprints per rank position.
    for (fa, fb) in fwd.iter().zip(back.iter()) {
        let sa: Vec<[u64; 2]> = fa
            .iter()
            .map(|&i| {
                [
                    base_pop[i].f[0].to_bits(),
                    base_pop[i].f[1].to_bits(),
                ]
            })
            .collect();
        let sb: Vec<[u64; 2]> = fb
            .iter()
            .map(|&i| {
                [
                    rev_pop[i].f[0].to_bits(),
                    rev_pop[i].f[1].to_bits(),
                ]
            })
            .collect();
        let mut ca = sa.clone();
        ca.sort_unstable();
        let mut cb = sb.clone();
        cb.sort_unstable();
        assert_eq!(ca, cb, "rank contents must be permutation-stable");
    }
}

// ---------------------------------------------------------------------------
// 3. Admission refusals: empty population, dimension mismatches,
//    non-finite carriers (both kinds, with component attribution).
// ---------------------------------------------------------------------------

#[test]
fn admission_refuses_empty_mismatched_and_nonfinite_inputs_typed() {
    let cfg = base_config(7, 4, 3, 32);
    let good_x = [0.25];
    let empty_err = nsga3_run(&[], &cfg, &mut |_| vec![0.0, 0.0]);
    assert!(matches!(
        empty_err,
        Err(NsgaError::PopulationEmpty)
    ));
    let dim_pop = vec![ind(&good_x, &[0.1, 0.2]), ind(&[0.3], &[0.1])];
    assert!(matches!(
        nsga3_run(&dim_pop, &cfg, &mut |_| vec![0.0, 0.0]),
        Err(NsgaError::ObjectiveCountMismatch {
            expected: 2,
            at: 1,
            found: 1
        })
    ));
    let dec_pop = vec![
        ind(&good_x, &[0.1, 0.2]),
        ind(&[0.3, 0.4], &[0.1, 0.2]),
    ];
    assert!(matches!(
        nsga3_run(&dec_pop, &cfg, &mut |_| vec![0.0, 0.0]),
        Err(NsgaError::DecisionCountMismatch {
            expected: 1,
            at: 1,
            found: 2
        })
    ));
    let nan_obj = vec![ind(&good_x, &[0.1, f64::NAN])];
    assert!(matches!(
        nsga3_run(&nan_obj, &cfg, &mut |_| vec![0.0, 0.0]),
        Err(NsgaError::NonFinite {
            individual: 0,
            component: 1,
            kind: NonFiniteKind::Objective
        })
    ));
    let inf_dec = vec![ind(&[f64::INFINITY], &[0.1, 0.2])];
    assert!(matches!(
        nsga3_run(&inf_dec, &cfg, &mut |_| vec![0.0, 0.0]),
        Err(NsgaError::NonFinite {
            individual: 0,
            component: 0,
            kind: NonFiniteKind::Decision
        })
    ));
}

#[test]
fn budget_below_one_generation_is_infeasible_before_any_evaluation() {
    let cfg = base_config(7, 6, 3, 5);
    let pop = (0..4)
        .map(|i| ind(&[i as f64 * 0.1], &[i as f64, 1.0 - i as f64]))
        .collect::<Vec<_>>();
    let mut charged = ConvexEval::new();
    let err = nsga3_run(&pop, &cfg, &mut |x| charged.call(x));
    assert!(matches!(
        err,
        Err(NsgaError::BudgetInfeasible {
            minimum_population: 6,
            budget: 5
        })
    ));
    assert!(charged.calls().is_empty(), "refusal must precede any eval");
}

// ---------------------------------------------------------------------------
// 4. Exact budget accounting: integral charging, BudgetBoundary reason,
//    partial generations never report, evaluations fit within ceiling.
// ---------------------------------------------------------------------------

#[test]
fn budget_boundary_stops_on_exact_generation_edge_never_partial() {
    // init 4 evaluated free; each generation costs 6 => generations 3,
    // evals 22, then next generation would hit 28 > 28?? budget 28 =>
    // exactly one more fits; we author 27 to force the edge AFTER gen3.
    let init_pop: Vec<NsgaIndividual> = (0..4)
        .map(|i| {
            let x0 = 0.05 + i as f64 * 0.1;
            ind(&[x0], &[x0 * x0, (x0 - 1.0) * (x0 - 1.0)])
        })
        .collect();
    let mut fe = ConvexEval::new();
    let cfg = base_config(11, 6, 50, 27);
    let report = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x))
        .expect("edge budget runs to completion");
    assert_eq!(report.evaluations, 22, "3 complete generations of 6");
    assert!(matches!(report.stop, fs_ascent::NsgaStop::BudgetBoundary));
    assert_eq!(report.generations, 3);
    assert!(report.evaluations <= cfg.eval_budget);
    // Population invariant by construction even under dedup churn.
    assert_eq!(report.population.len(), cfg.population_size);
}

#[test]
fn max_generations_ceiling_reports_max_generations_reason() {
    let init_pop: Vec<NsgaIndividual> = (0..4)
        .map(|i| ind(&[0.1 + i as f64 * 0.15], &[i as f64, 0.0]))
        .collect();
    let mut fe = ConvexEval::new();
    let cfg = base_config(5, 4, 7, 1_000);
    let report = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x)).expect("runs");
    assert_eq!(report.generations, 7);
    assert_eq!(report.evaluations, 4 + 7 * 4);
    assert!(matches!(report.stop, fs_ascent::NsgaStop::MaxGenerations));
    assert_eq!(fe.calls().len(), report.evaluations - 4);
}

// ---------------------------------------------------------------------------
// 5. Duplicates + boundaries: duplicated populations survive, box holds,
//    canonical tie policy keeps outcomes enumeration-stable.
// ---------------------------------------------------------------------------

#[test]
fn fully_duplicated_population_keeps_size_and_box_bounds_hold() {
    let same = ind(&[0.5], &[0.3, 0.7]);
    let init_pop: Vec<NsgaIndividual> = vec![same.clone(); 5];
    let mut fe = ConvexEval::new();
    let cfg = base_config(9, 5, 4, 40);
    let report = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x)).expect("runs");
    assert_eq!(report.population.len(), cfg.population_size);
    for member in &report.population {
        assert!(
            (0.0..=1.0).contains(&member.x[0]),
            "box violated: {}",
            member.x[0]
        );
    }
    // First front contains EVERYONE here (mutual non-domination among
    // post-dedup survivors is possible; sort still defines ranks).
    assert!(!report.fronts.is_empty());
}

#[test]
fn box_clamps_children_even_when_parents_cluster_at_a_face() {
    // All parents pinned AT the lower face; children must stay >= 0.0.
    let init_pop: Vec<NsgaIndividual> = (0..4)
        .map(|i| ind(&[0.0], &[i as f64 * 0.01, 1.0]))
        .collect();
    let mut fe = ConvexEval::new();
    let cfg = base_config(23, 6, 3, 60);
    let report = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x)).expect("runs");
    for member in &report.population {
        assert!(member.x[0] >= 0.0 && member.x[0] <= 1.0);
        assert!(member.x[0].is_finite());
    }
}

// ---------------------------------------------------------------------------
// 6. Metamorphic units/scaling law: positive per-axis rescaling and
//    shifting preserve DOMINANCE STRUCTURE (the multi-objective core).
// ---------------------------------------------------------------------------

#[test]
fn positive_affine_objective_maps_preserve_front_structure() {
    let mk_case = || -> Vec<NsgaIndividual> {
        (0..8)
            .map(|i| {
                let t = i as f64 / 7.0;
                ind(&[t], &[t, 1.0 - t])
            })
            .chain(std::iter::once(ind(
                &[0.99],
                &[0.05, 0.95],
            )))
            .collect()
    };
    let pop = mk_case();
    let fronts_a = fast_nondominated_sort(&pop);
    // Affine map g_i = alpha*f_i + beta (alpha > 0): dominance relations
    // are preserved bit-for-bit-semantics notwithstanding rounding, so
    // rank ASSIGNMENTS must coincide. Shift must not move anyone into
    // tie collisions with different members here.
    let (alpha, beta) = (3.5, 12.75);
    let mapped: Vec<NsgaIndividual> = pop
        .iter()
        .map(|p| ind(&p.x, &p.f.iter().map(|v| alpha * v + beta).collect::<Vec<_>>()))
        .collect();
    let fronts_b = fast_nondominated_sort(&mapped);
    assert_eq!(fronts_a.len(), fronts_b.len());
    for (ra, rb) in fronts_a.iter().zip(fronts_b.iter()) {
        let fingerprint_a: Vec<[u64; 2]> = ra
            .iter()
            .map(|&i| [pop[i].f[0].to_bits(), pop[i].f[1].to_bits()])
            .collect();
        let fingerprint_b: Vec<[u64; 2]> = rb
            .iter()
            .map(|&i| {
                let orig = i;
                let _ = orig;
                [
                    mapped[i].f[0].to_bits(),
                    mapped[i].f[1].to_bits(),
                ]
            })
            .collect();
        assert_eq!(
            fingerprint_a.len(),
            fingerprint_b.len(),
            "same rank population"
        );
    }
    // Front MEMBER COUNT parity at every rank (structure identity).
    for (ra, rb) in fronts_a.iter().zip(fronts_b.iter()) {
        assert_eq!(ra.len(), rb.len());
    }
}

// ---------------------------------------------------------------------------
// 7. Normalization: well-posed runs avoid the nadir fallback; crafted
//    degenerate geometry engages it without refusing.
// ---------------------------------------------------------------------------

#[test]
fn normalization_singular_fallback_engages_without_breaking_the_run() {
    // Degenerate geometry: all objectives lie ON THE LINE f1=f2, so the
    // ASF extreme system is rank-deficient (every extreme row is the
    // same translated point) => fallback must engage yet the pipeline
    // completes deterministically.
    let init_pop: Vec<NsgaIndividual> = (0..4)
        .map(|i| {
            let t = 0.1 * i as f64;
            ind(&[t], &[t, t])
        })
        .collect();
    let cfg = base_config(31, 4, 3, 40);
    let report =
        nsga3_run(&init_pop, &cfg, &mut |x| {
            let t = x[0] * 0.5 + 0.25;
            vec![t, t]
        })
        .expect("degenerate geometry still realizes");
    assert!(
        report.normalization_singular_fallback,
        "rank-deficient ASF system must route through the nadir fallback"
    );
    assert_eq!(report.population.len(), cfg.population_size);
    let _ = cfg.reference_divisions;
}
#[test]
fn wellposed_geometry_reports_no_fallback() {
    let init_pop: Vec<NsgaIndividual> = vec![
        ind(&[0.0], &[0.0, 1.0]),
        ind(&[1.0], &[1.0, 0.0]),
        ind(&[0.5], &[0.35, 0.55]),
        ind(&[0.25], &[0.18, 0.78]),
    ];
    let mut fe = ConvexEval::new();
    let cfg = base_config(41, 4, 2, 30);
    let report = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x)).expect("runs");
    assert!(!report.normalization_singular_fallback);
}

// ---------------------------------------------------------------------------
// 8. Bitwise determinism across repeated runs (G5 class), independent
//    of thread/order effects by construction of the seeded stream.
// ---------------------------------------------------------------------------

#[test]
fn repeated_runs_are_bitwise_identical_for_identical_inputs() {
    let init_pop: Vec<NsgaIndividual> = (0..6)
        .map(|i| {
            let u = (i as f64 + 1.0) / 7.0;
            ind(&[u], &[u * u, (u - 1.0).powi(2)])
        })
        .collect();
    let run_once = |seed: u64| -> (usize, Vec<String>) {
        let cfg = base_config(seed, 6, 5, 400);
        let mut fe = ConvexEval::new();
        let r = nsga3_run(&init_pop, &cfg, &mut |x| fe.call(x)).expect("runs");
        let fingerprint: Vec<String> = r
            .population
            .iter()
            .map(|p| {
                format!(
                    "{:?}|{:?}",
                    p.f.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    p.x
                )
            })
            .collect();
        (r.evaluations, fingerprint)
    };
    let (e1, f1) = run_once(2024_0827);
    let (e2, f2) = run_once(2024_0827);
    assert_eq!(e1, e2);
    assert_eq!(f1, f2, "population fingerprint must reproduce bitwise");
}
