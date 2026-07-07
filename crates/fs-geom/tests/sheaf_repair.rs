//! Sheaf-repair conformance (the wqd.14 bead; runs under the
//! `sheaf-repair` feature). Acceptance: exact-component defects
//! auto-repair to certificate-passing without exceeding chart budgets;
//! coexact seeding (flipped orientation) diagnoses as converter-side;
//! harmonic seeding declares unrepairable-locally with the correct
//! cut-set; predicted post-repair norms match actuals; the decomposition
//! matches a dense oracle; repair is idempotent and budget-safe.
#![cfg(feature = "sheaf-repair")]

use fs_geom::router::{ConverterSpec, ErrorModel, MemoryCostOracle, RouteRequest, Router};
use fs_geom::sheaf_repair::{SheafSkeleton, apply_gauge, hodge_decompose, plan_repair};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf-repair\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A 3-patch triangle complex (one triple junction).
fn triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// A 4-patch ring (cycle, NO triangles): H¹ is nontrivial by design.
fn ring() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 2), (2, 3), (0, 3)],
        triangles: vec![],
    }
}

fn norm_inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
}

/// Dense least-squares oracle: minimize ‖m − A x‖² by normal equations
/// solved with Gaussian elimination (partial pivot) — an independent
/// code path from the module's Gauss–Seidel.
fn dense_projection(m: &[f64], columns: &[Vec<f64>]) -> Vec<f64> {
    let n = columns.len();
    let mut ata = vec![vec![0.0f64; n]; n];
    let mut atb = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            ata[i][j] = columns[i].iter().zip(&columns[j]).map(|(a, b)| a * b).sum();
        }
        atb[i] = columns[i].iter().zip(m).map(|(a, b)| a * b).sum();
    }
    // Ridge the (rank-deficient) gauge direction minimally.
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += 1e-12;
    }
    // Gaussian elimination.
    let mut aug: Vec<Vec<f64>> = ata
        .iter()
        .zip(&atb)
        .map(|(row, &b)| {
            let mut r = row.clone();
            r.push(b);
            r
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&a, &b| aug[a][col].abs().total_cmp(&aug[b][col].abs()))
            .expect("rows");
        aug.swap(col, pivot);
        let p = aug[col][col];
        if p.abs() < 1e-300 {
            continue;
        }
        for r in 0..n {
            if r != col {
                let f = aug[r][col] / p;
                for k in col..=n {
                    aug[r][k] -= f * aug[col][k];
                }
            }
        }
    }
    (0..n)
        .map(|i| {
            if aug[i][i].abs() < 1e-300 {
                0.0
            } else {
                aug[i][n] / aug[i][i]
            }
        })
        .collect()
}

#[test]
fn sr_001_decomposition_matches_dense_oracle() {
    let sk = triangle();
    // A mixed cochain: gauge part + circulation part.
    let gauge_part = sk.d0(&[0.0, 0.7, -0.3]);
    let circ_part = sk.d1t(&[0.4]);
    let m: Vec<f64> = gauge_part
        .iter()
        .zip(&circ_part)
        .map(|(a, b)| a + b)
        .collect();
    let split = hodge_decompose(&sk, &m);
    // Oracle: dense projection onto im δ⁰ (columns = δ⁰ of unit vertex
    // vectors, vertex 0 pinned) and im δ¹ᵀ.
    let d0_cols: Vec<Vec<f64>> = (1..sk.n_patches)
        .map(|i| {
            let mut e = vec![0.0; sk.n_patches];
            e[i] = 1.0;
            sk.d0(&e)
        })
        .collect();
    let c_oracle = dense_projection(&m, &d0_cols);
    let exact_oracle = {
        let mut full = vec![0.0; sk.n_patches];
        full[1..].copy_from_slice(&c_oracle);
        sk.d0(&full)
    };
    for (got, want) in split.exact.iter().zip(&exact_oracle) {
        assert!(
            (got - want).abs() < 1e-8,
            "exact component vs dense oracle: {got} vs {want}"
        );
    }
    // On a triangle (contractible), harmonic must vanish.
    assert!(
        norm_inf(&split.harmonic) < 1e-8,
        "contractible complex has no harmonic part: {:?}",
        split.harmonic
    );
    // Orthogonality residuals: δ⁰ᵀh ≈ 0 and δ¹h ≈ 0.
    assert!(norm_inf(&sk.d0t(&split.harmonic)) < 1e-8);
    assert!(norm_inf(&sk.d1(&split.harmonic)) < 1e-8);
    // Energy fractions sum to ~1 on an orthogonal split.
    let (fe, fc, fh) = split.fractions;
    assert!(
        (fe + fc + fh - 1.0).abs() < 1e-6,
        "fractions partition energy: {fe} + {fc} + {fh}"
    );
    verdict(
        "sr-001",
        "exact component matches the dense-oracle projection; contractible harmonic \
         vanishes; energy partitions",
    );
}

#[test]
fn sr_002_exact_defect_auto_repairs_within_budget() {
    let sk = triangle();
    // Seed a pure gauge defect: patch 2 drifted by +0.012.
    let mismatch = sk.d0(&[0.0, 0.0, 0.012]);
    let budgets = [0.02, 0.02, 0.02];
    let plan = plan_repair(&sk, &mismatch, &budgets, None);
    assert!(plan.auto_repairable, "within budgets: auto-repairable");
    assert!(plan.split.fractions.0 > 0.999, "pure exact defect");
    assert!(plan.obstruction_cutset.is_empty(), "no obstruction");
    // Predicted-vs-actual: apply the gauge, re-measure.
    let predicted = plan.proposals[0].expected_post_norm;
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    let actual = norm_inf(&repaired);
    assert!(
        (predicted - actual).abs() < 1e-9,
        "prediction {predicted} vs actual {actual}"
    );
    assert!(actual < 1e-9, "certificate-passing after repair");
    // Repair SAFETY: offsets stay within each chart's declared budget.
    for (off, b) in plan.gauge.iter().zip(&budgets) {
        assert!(off.abs() <= *b, "repair never exceeds a budget");
    }
    // Repair IDEMPOTENCE: repairing the repaired model is a no-op.
    let plan2 = plan_repair(&sk, &repaired, &budgets, None);
    assert!(
        norm_inf(&plan2.gauge) < 1e-9,
        "no residual gauge on a passing model: {:?}",
        plan2.gauge
    );
    let repaired2 = apply_gauge(&sk, &repaired, &plan2.gauge);
    assert!(
        (norm_inf(&repaired2) - actual).abs() < 1e-12,
        "no-op repair"
    );
    // Over-budget variant: the SAME defect with a tight budget must NOT
    // auto-apply (needs explicit acceptance).
    let tight = [0.001, 0.001, 0.001];
    let gated = plan_repair(&sk, &mismatch, &tight, None);
    assert!(
        !gated.auto_repairable,
        "budget gate blocks silent distortion"
    );
    assert!(
        gated.proposals[0].action.contains("EXCEEDS"),
        "the proposal says so: {}",
        gated.proposals[0].action
    );
    verdict(
        "sr-002",
        "gauge defect repaired to ~0 with exact prediction; idempotent; budget gate \
         blocks over-budget auto-apply",
    );
}

#[test]
fn sr_003_coexact_seeding_diagnoses_converter_side() {
    let sk = triangle();
    // Seed a pure circulation (the flipped-orientation signature): the
    // image of δ¹ᵀ.
    let mismatch = sk.d1t(&[0.05]);
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], None);
    assert!(
        plan.split.fractions.1 > 0.999,
        "pure coexact defect: {:?}",
        plan.split.fractions
    );
    let converter_proposal = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("CONVERTER"))
        .expect("converter-side diagnosis present");
    assert!(
        converter_proposal.action.contains("(0, 1, 2)"),
        "localized to the triple junction: {}",
        converter_proposal.action
    );
    // Gauge repair CANNOT fix circulation: applying it leaves the norm.
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    assert!(
        norm_inf(&repaired) > 0.9 * norm_inf(&mismatch),
        "circulation is not gauge-repairable"
    );
    verdict(
        "sr-003",
        "circulation seeding is >99.9% coexact, diagnosed converter-side at the right \
         junction, and provably not gauge-fixable",
    );
}

#[test]
fn sr_004_harmonic_seeding_declares_unrepairable_with_cutset() {
    let sk = ring();
    // A circulation around the 4-cycle: with no 2-cells, nothing coexact
    // exists and no gauge kills a loop sum — genuinely harmonic.
    // Orientation: edges (0,1),(1,2),(2,3) run low→high along the loop;
    // (0,3) runs AGAINST it, so the loop cochain is (ε, ε, ε, −ε).
    let eps = 0.03;
    let mismatch = vec![eps, eps, eps, -eps];
    let plan = plan_repair(&sk, &mismatch, &[1.0; 4], None);
    assert!(
        plan.split.fractions.2 > 0.999,
        "pure harmonic: {:?}",
        plan.split.fractions
    );
    assert!(!plan.auto_repairable || plan.split.fractions.0 < 1e-9);
    assert_eq!(
        plan.obstruction_cutset.len(),
        4,
        "the whole cycle is the cut-set"
    );
    let unrepairable = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("NO local fix"))
        .expect("honest unrepairable proposal");
    assert!(unrepairable.cost_s.is_infinite(), "no local cost claim");
    // And indeed gauge repair achieves nothing.
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    assert!(norm_inf(&repaired) > 0.9 * eps, "harmonic survives gauge");
    verdict(
        "sr-004",
        "cycle circulation is >99.9% harmonic, declared unrepairable-locally with the \
         full-cycle cut-set",
    );
}

#[test]
fn sr_005_router_reroute_proposal_ranks_by_expected_norm() {
    let sk = triangle();
    let mismatch = sk.d0(&[0.0, 0.0, 0.012]);
    // A router with one certified conversion available for the worst
    // patch's chart kind.
    let mut router = Router::new();
    router
        .register(ConverterSpec {
            name: "sdf->mesh/dc-interval".to_string(),
            from: "sdf".to_string(),
            to: "mesh".to_string(),
            base_cost_s: 2.0,
            error: ErrorModel::AdditiveAbs(5e-7),
            certified: true,
        })
        .expect("register");
    let oracle = MemoryCostOracle::new();
    let req = RouteRequest {
        from: "sdf".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 100.0,
    };
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], Some((&router, &oracle, &req)));
    let reroute = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("reroute"))
        .expect("router proposal present");
    assert!(reroute.action.contains("dc-interval"), "{}", reroute.action);
    assert!((reroute.cost_s - 2.0).abs() < 1e-9, "router-modeled cost");
    // Ranking: proposals sorted by expected post-repair norm — the gauge
    // repair (→ ~0) outranks the reroute (→ 5e-7) only if smaller; both
    // must be ordered non-decreasingly.
    for pair in plan.proposals.windows(2) {
        assert!(
            pair[0].expected_post_norm <= pair[1].expected_post_norm + 1e-12,
            "proposals ranked by expected norm"
        );
    }
    verdict(
        "sr-005",
        "router reroute proposal carries the planned route + modeled cost; ranking is \
         by expected post-repair norm",
    );
}
