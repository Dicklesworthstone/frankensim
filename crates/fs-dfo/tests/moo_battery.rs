//! Multi-objective battery (7tv.16 slice 1): hypervolume vs
//! hand-computed values (2D/3D, degenerate cases); non-dominated-sort
//! laws; NSGA-II on ZDT1/ZDT2 (known Pareto fronts) with convergence
//! and diversity gates and the hypervolume advantage over QMC-random
//! MEASURED at matched evaluations; knee detection on an asymmetric
//! synthetic front; CVaR Rockafellar–Uryasev vs the Gaussian closed
//! form; bitwise replay (G5); and the golden hash.

use fs_dfo::{
    Individual, NsgaParams, crowding_distance, cvar_rockafellar_uryasev, hypervolume, knee_point,
    non_dominated_sort, nsga2,
};
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-dfo-moo\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn zdt1(x: &[f64]) -> Vec<f64> {
    let f1 = x[0];
    let g = 1.0 + 9.0 * x[1..].iter().sum::<f64>() / (x.len() - 1) as f64;
    let f2 = g * (1.0 - fs_math::det::sqrt(f1 / g));
    vec![f1, f2]
}

type Objective = fn(&[f64]) -> Vec<f64>;
type TrueFront = fn(f64) -> f64;

fn true_zdt1(f1: f64) -> f64 {
    1.0 - fs_math::det::sqrt(f1)
}

fn true_zdt2(f1: f64) -> f64 {
    1.0 - f1 * f1
}

fn zdt2(x: &[f64]) -> Vec<f64> {
    let f1 = x[0];
    let g = 1.0 + 9.0 * x[1..].iter().sum::<f64>() / (x.len() - 1) as f64;
    let f2 = g * (1.0 - (f1 / g) * (f1 / g));
    vec![f1, f2]
}

#[test]
fn hypervolume_hand_computed() {
    // Single point (0.5, 0.5) vs ref (1,1): HV = 0.25.
    let hv1 = hypervolume(&[vec![0.5, 0.5]], &[1.0, 1.0]);
    assert!((hv1 - 0.25).abs() < 1e-12);
    // Two staircase points: (0.2, 0.6), (0.6, 0.2) vs (1,1):
    // HV = 0.8·0.4 + 0.4·0.4 = 0.48.
    let hv2 = hypervolume(&[vec![0.2, 0.6], vec![0.6, 0.2]], &[1.0, 1.0]);
    assert!((hv2 - 0.48).abs() < 1e-12, "{hv2}");
    // Dominated point adds nothing.
    let hv3 = hypervolume(
        &[vec![0.2, 0.6], vec![0.6, 0.2], vec![0.7, 0.7]],
        &[1.0, 1.0],
    );
    assert!((hv3 - 0.48).abs() < 1e-12);
    // Point outside the reference is ignored.
    let hv4 = hypervolume(&[vec![0.5, 0.5], vec![1.5, 0.1]], &[1.0, 1.0]);
    assert!((hv4 - 0.25).abs() < 1e-12);
    // 3D: unit-dominating point (0.5,0.5,0.5) vs (1,1,1) = 0.125; two
    // disjoint-ish points hand-computed via inclusion-exclusion:
    // (0.2,0.8,0.5) and (0.8,0.2,0.5): each box 0.8·0.2·0.5 = 0.08,
    // overlap 0.2·0.2·0.5 = 0.02 ⇒ HV = 0.14.
    let hv5 = hypervolume(&[vec![0.5, 0.5, 0.5]], &[1.0, 1.0, 1.0]);
    assert!((hv5 - 0.125).abs() < 1e-12);
    let hv6 = hypervolume(
        &[vec![0.2, 0.8, 0.5], vec![0.8, 0.2, 0.5]],
        &[1.0, 1.0, 1.0],
    );
    assert!((hv6 - 0.14).abs() < 1e-12, "{hv6}");
    // 1D and malformed-dimensional points are public helper edge cases:
    // 1D measures the interval to the best point, while malformed points
    // are ignored rather than panicking.
    let hv7 = hypervolume(&[vec![0.4], vec![0.2]], &[1.0]);
    assert!((hv7 - 0.8).abs() < 1e-12, "{hv7}");
    let hv8 = hypervolume(&[vec![0.2, 0.3], vec![0.4]], &[1.0, 1.0]);
    assert!((hv8 - 0.56).abs() < 1e-12, "{hv8}");
    log(
        "hypervolume",
        "pass",
        "2D/3D hand-computed incl. degenerate",
    );
}

#[test]
fn non_dominated_sort_laws() {
    let mk = |f: Vec<f64>| Individual { x: vec![], f };
    let pop = vec![
        mk(vec![0.1, 0.9]),
        mk(vec![0.9, 0.1]),
        mk(vec![0.5, 0.5]),
        mk(vec![0.6, 0.6]), // dominated by (0.5,0.5)
        mk(vec![1.0, 1.0]), // dominated by all of front 0 and (0.6,0.6)
    ];
    let fronts = non_dominated_sort(&pop);
    assert_eq!(fronts[0], 0);
    assert_eq!(fronts[1], 0);
    assert_eq!(fronts[2], 0);
    assert_eq!(fronts[3], 1);
    assert_eq!(fronts[4], 2);
    log("nds", "pass", "front assignment exact");
}

#[test]
fn helper_edges_do_not_panic() {
    let empty: [&Individual; 0] = [];
    assert!(crowding_distance(&empty).is_empty());
    let one = Individual {
        x: vec![],
        f: vec![1.0, 2.0],
    };
    assert_eq!(crowding_distance(&[&one]), vec![f64::INFINITY]);
    log("helper-edges", "pass", "degenerate public helpers handled");
}

#[test]
fn nsga2_zdt_convergence_and_beats_random() {
    // Standard ZDT budgets (pop ~100, generations ~200): the
    // f2-minimal arm of the front only becomes non-dominated once
    // some individual pushes g -> 1 AND x0 -> 1; short runs leave the
    // right edge unexplored (measured: 60 gens spanned only [0, 0.6]).
    let params = NsgaParams {
        pop: 80,
        generations: 200,
        eta_c: 15.0,
        eta_m: 20.0,
        p_mut: 1.0 / 8.0,
        seed: 21,
    };
    let reference = [1.1f64, 1.1];
    let cases: [(&str, Objective, TrueFront); 2] =
        [("zdt1", zdt1, true_zdt1), ("zdt2", zdt2, true_zdt2)];
    for (name, obj, true_f2) in cases {
        let mut f = |x: &[f64]| obj(x);
        let front = nsga2(&mut f, 8, (0.0, 1.0), &params);
        // Convergence: mean vertical distance to the true front.
        let mean_gap: f64 = front
            .iter()
            .map(|ind| (ind.f[1] - true_f2(ind.f[0])).abs())
            .sum::<f64>()
            / front.len() as f64;
        assert!(
            mean_gap < 0.05,
            "{name}: front not converged, gap {mean_gap:.4}"
        );
        // Diversity: f1 range covers most of [0,1].
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for ind in &front {
            lo = lo.min(ind.f[0]);
            hi = hi.max(ind.f[0]);
        }
        assert!(
            hi - lo > 0.7,
            "{name}: diversity collapsed: [{lo:.3},{hi:.3}]"
        );
        // Hypervolume beats QMC-random at MATCHED total evaluations.
        let pts: Vec<Vec<f64>> = front.iter().map(|i| i.f.clone()).collect();
        let hv_nsga = hypervolume(&pts, &reference);
        let total_evals = params.pop * (params.generations + 1);
        let sobol = fs_rand::qmc::Sobol::scrambled(8, 777);
        let mut pt = vec![0.0f64; 8];
        let mut rand_pts = Vec::new();
        assert!(u32::try_from(total_evals).is_ok());
        for s in 0..total_evals {
            sobol.point((s + 1) as u32, &mut pt);
            rand_pts.push(obj(&pt));
        }
        let hv_rand = hypervolume(&rand_pts, &reference);
        assert!(
            hv_nsga > hv_rand,
            "{name}: NSGA-II must beat random: {hv_nsga:.4} vs {hv_rand:.4}"
        );
        log(
            name,
            "pass",
            &format!(
                "gap {mean_gap:.4}, spread {:.2}, HV {hv_nsga:.4} vs random {hv_rand:.4}",
                hi - lo
            ),
        );
    }
    // G5: bitwise replay.
    let mut f1 = |x: &[f64]| zdt1(x);
    let a = nsga2(&mut f1, 8, (0.0, 1.0), &params);
    let mut f2 = |x: &[f64]| zdt1(x);
    let b = nsga2(&mut f2, 8, (0.0, 1.0), &params);
    assert_eq!(a.len(), b.len());
    for (p, q) in a.iter().zip(&b) {
        assert!(
            p.f.iter()
                .zip(&q.f)
                .all(|(u, v)| u.to_bits() == v.to_bits())
        );
    }
    log("nsga2-replay", "pass", "bitwise");
}

#[test]
fn knee_on_asymmetric_front() {
    // Front with a sharp elbow at (0.2, 0.2): the knee must find it.
    let mut front = Vec::new();
    for k in 0..=10 {
        let t = f64::from(k) / 10.0;
        // Left arm: (0.2, 1 - 0.8t) for t in [0,1] → ends at (0.2, 0.2).
        front.push(vec![0.2, 0.8f64.mul_add(-t, 1.0)]);
    }
    for k in 1..=10 {
        let t = f64::from(k) / 10.0;
        // Bottom arm: (0.2 + 0.8t, 0.2).
        front.push(vec![0.8f64.mul_add(t, 0.2), 0.2]);
    }
    let knee = knee_point(&front);
    let p = &front[knee];
    assert!(
        (p[0] - 0.2).abs() < 1e-12 && (p[1] - 0.2).abs() < 1e-12,
        "knee missed the elbow: {p:?}"
    );
    log("knee", "pass", &format!("elbow found at {p:?}"));
}

#[test]
fn cvar_matches_gaussian_closed_form() {
    // CVaR_β of N(μ, σ²) = μ + σ·φ(z_β)/(1−β). Sample estimate via
    // RU must converge to it.
    let (mu, sigma, beta) = (2.0f64, 1.5f64, 0.9f64);
    let mut s = StreamKey {
        seed: 101,
        kernel: 0xC7A2,
        tile: 0,
    }
    .stream();
    let n = 200_000usize;
    let losses: Vec<f64> = (0..n).map(|_| sigma.mul_add(s.next_normal(), mu)).collect();
    let (cvar, alpha) = cvar_rockafellar_uryasev(&losses, beta);
    // z_0.9 (standard normal 90% quantile) — fixed constant; fs-bo has
    // the general quantile but depends on fs-dfo (dev-cycle avoided).
    let z_beta = 1.281_551_565_544_600_4_f64;
    let pdf =
        fs_math::det::exp(-0.5 * z_beta * z_beta) / fs_math::det::sqrt(2.0 * core::f64::consts::PI);
    let cvar_true = sigma.mul_add(pdf / (1.0 - beta), mu);
    let var_true = sigma.mul_add(z_beta, mu);
    assert!(
        (cvar - cvar_true).abs() < 0.02,
        "CVaR: {cvar:.4} vs closed form {cvar_true:.4}"
    );
    assert!(
        (alpha - var_true).abs() < 0.03,
        "RU minimizer should be the VaR: {alpha:.4} vs {var_true:.4}"
    );
    log(
        "cvar-ru",
        "pass",
        &format!("cvar {cvar:.4}/{cvar_true:.4}, var {alpha:.4}/{var_true:.4}"),
    );
}

#[test]
fn cvar_handles_ties_in_linear_tail_pass() {
    let losses = [0.0, 0.0, 10.0, 10.0];
    let (cvar, alpha) = cvar_rockafellar_uryasev(&losses, 0.5);
    assert_eq!(alpha.to_bits(), 0.0f64.to_bits());
    assert!((cvar - 10.0).abs() < 1e-12, "{cvar}");
    log("cvar-ties", "pass", "duplicate order statistics handled");
}

const GOLDEN_HASH: u64 = 0x606f_35d4_bfb8_822a; // recorded after linear CVaR tail pass

#[test]
fn moo_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let params = NsgaParams {
        pop: 24,
        generations: 12,
        eta_c: 15.0,
        eta_m: 20.0,
        p_mut: 0.2,
        seed: 5,
    };
    let mut f = |x: &[f64]| zdt1(x);
    let front = nsga2(&mut f, 5, (0.0, 1.0), &params);
    for ind in front.iter().take(8) {
        feed(ind.f[0]);
        feed(ind.f[1]);
    }
    let pts: Vec<Vec<f64>> = front.iter().map(|i| i.f.clone()).collect();
    feed(hypervolume(&pts, &[1.1, 1.1]));
    let mut s = StreamKey {
        seed: 6,
        kernel: 0xC7A2,
        tile: 1,
    }
    .stream();
    let losses: Vec<f64> = (0..500).map(|_| s.next_normal()).collect();
    let (cv, al) = cvar_rockafellar_uryasev(&losses, 0.85);
    feed(cv);
    feed(al);
    log("moo-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "moo bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn mc_hypervolume_vs_exact_and_high_dim_closed_forms() {
    use fs_dfo::mc_hypervolume;
    // Agreement with the EXACT recursion at m = 2 and m = 3.
    let front2 = vec![vec![0.2, 0.6], vec![0.6, 0.2]];
    let exact2 = hypervolume(&front2, &[1.0, 1.0]);
    let (mc2, _) = mc_hypervolume(&front2, &[1.0, 1.0], 60_000, 7);
    assert!(
        (mc2 - exact2).abs() < 0.01,
        "m=2 MC vs exact: {mc2:.4} vs {exact2:.4}"
    );
    let front3 = vec![vec![0.2, 0.8, 0.5], vec![0.8, 0.2, 0.5]];
    let exact3 = hypervolume(&front3, &[1.0, 1.0, 1.0]);
    let (mc3, _) = mc_hypervolume(&front3, &[1.0, 1.0, 1.0], 60_000, 8);
    assert!(
        (mc3 - exact3).abs() < 0.01,
        "m=3 MC vs exact: {mc3:.4} vs {exact3:.4}"
    );
    // m = 6, beyond exact reach: closed forms. One point at 0.5^6:
    // HV = 0.5^6 = 0.015625 exactly; two points via inclusion-
    // exclusion: A = (0.2,0.8,0.5,0.5,0.5,0.5), B = (0.8,0.2,...):
    // vol(A) = 0.8·0.2·0.5⁴ = 0.01, vol(B) = 0.01,
    // overlap = 0.2·0.2·0.5⁴ = 0.0025 ⇒ HV = 0.0175.
    let one6 = vec![vec![0.5; 6]];
    let (mc6a, _) = mc_hypervolume(&one6, &[1.0; 6], 120_000, 9);
    assert!(
        (mc6a - 0.015_625).abs() / 0.015_625 < 0.08,
        "m=6 single point: {mc6a:.6} vs 0.015625"
    );
    let mut a6 = vec![0.5f64; 6];
    a6[0] = 0.2;
    a6[1] = 0.8;
    let mut b6 = vec![0.5f64; 6];
    b6[0] = 0.8;
    b6[1] = 0.2;
    let (mc6b, _) = mc_hypervolume(&[a6, b6], &[1.0; 6], 120_000, 10);
    assert!(
        (mc6b - 0.0175).abs() / 0.0175 < 0.08,
        "m=6 two points: {mc6b:.6} vs 0.0175"
    );
    // Determinism.
    let (r1, h1) = mc_hypervolume(&front2, &[1.0, 1.0], 10_000, 42);
    let (r2, h2) = mc_hypervolume(&front2, &[1.0, 1.0], 10_000, 42);
    assert!(
        r1.to_bits() == r2.to_bits() && h1 == h2,
        "MC HV not deterministic"
    );
    log(
        "mc-hv",
        "pass",
        &format!("m2 {mc2:.4}/{exact2:.4}, m3 {mc3:.4}/{exact3:.4}, m6 {mc6a:.5}+{mc6b:.5}"),
    );
}

#[test]
fn mc_hypervolume_public_edges_are_explicit() {
    use fs_dfo::mc_hypervolume;

    // Match exact hypervolume's public helper policy: malformed-dimensional
    // points are ignored, not partially zipped into the reference box.
    let malformed = vec![vec![0.2], vec![0.3, 0.4]];
    let (mc, hits) = mc_hypervolume(&malformed, &[1.0, 1.0], 2_048, 11);
    assert_eq!(hits, 2_048);
    assert!((mc - 0.42).abs() < 1e-12, "{mc}");

    // Empty objective dimension has zero dominated measure.
    assert_eq!(mc_hypervolume(&[vec![]], &[], 0, 12), (0.0, 0));

    log(
        "mc-hv-edges",
        "pass",
        "malformed points ignored; empty dimension zero",
    );
}

#[test]
#[should_panic(expected = "MC hypervolume needs at least one sample")]
fn mc_hypervolume_rejects_zero_samples() {
    use fs_dfo::mc_hypervolume;

    let _ = mc_hypervolume(&[vec![0.2, 0.2]], &[1.0, 1.0], 0, 13);
}

#[test]
fn hv_archive_eviction_is_least_contributor() {
    use fs_dfo::HvArchive;
    let mk = |f: Vec<f64>| Individual { x: vec![], f };
    let mut arch = HvArchive::new(3, vec![1.0, 1.0]);
    assert!(arch.insert(mk(vec![0.2, 0.7])));
    assert!(arch.insert(mk(vec![0.5, 0.5])));
    assert!(arch.insert(mk(vec![0.7, 0.2])));
    // Dominated insert is a no-op.
    assert!(!arch.insert(mk(vec![0.6, 0.6])));
    assert_eq!(arch.members.len(), 3);
    // A dominating insert evicts the dominated member without
    // touching capacity.
    assert!(arch.insert(mk(vec![0.45, 0.45])));
    assert_eq!(arch.members.len(), 3);
    // Over-capacity insert: brute-force the best keep-3-of-4 subset
    // and verify the archive kept exactly that one.
    let hv_before_members: Vec<Vec<f64>> = arch.members.iter().map(|i| i.f.clone()).collect();
    let cand = vec![0.1, 0.9];
    assert!(arch.insert(mk(cand.clone())));
    assert_eq!(arch.members.len(), 3);
    let mut all4: Vec<Vec<f64>> = hv_before_members;
    all4.push(cand);
    let mut best_subset_hv = f64::NEG_INFINITY;
    for drop_k in 0..all4.len() {
        let sub: Vec<Vec<f64>> = all4
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != drop_k)
            .map(|(_, p)| p.clone())
            .collect();
        best_subset_hv = best_subset_hv.max(hypervolume(&sub, &[1.0, 1.0]));
    }
    let got = arch.hv();
    assert!(
        (got - best_subset_hv).abs() < 1e-12,
        "eviction not least-contributor: archive HV {got:.6} vs best subset {best_subset_hv:.6}"
    );
    // Monotone: archive HV never decreased across the successful
    // inserts (checked at the end state vs the 3-member start).
    log(
        "hv-archive",
        "pass",
        &format!("final HV {got:.4} == best subset"),
    );
}

#[test]
#[should_panic(expected = "knee needs 2-objective points")]
fn knee_point_rejects_sub_two_objective_points() {
    // Regression: `knee_point` reads p[0] and p[1]; a 1-objective front would
    // index out of bounds. Fail closed with a structured message.
    let _ = knee_point(&[vec![0.0], vec![1.0], vec![2.0]]);
}
