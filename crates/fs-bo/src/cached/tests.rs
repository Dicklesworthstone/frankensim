use super::*;
use crate::{Kernel, Matern};
use fs_la::factor::cholesky;

fn config() -> NoisyBoConfig {
    NoisyBoConfig { bounds: (-2.0, 2.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.4] },
        prior_mean: 0.0, q: 2, mc_samples: 32, acq_starts: 2, acq_evals: 24, seed: 73 }
}
fn gp() -> Gp { Gp::fit(&[vec![0.0]], &[0.2], config().kernel, 0.3) }

// One-training-point posterior evaluated from scalar Gaussian equations, with
// a full joined factorization rather than a prepared baseline/extension.
fn dense_score(gp: &Gp, base: &[Vec<f64>], candidates: &[Vec<f64>], bank: &[f64], jitter: f64) -> f64 {
    let p: Vec<_> = base.iter().chain(candidates).cloned().collect();
    let n = p.len(); let mut sigma = vec![0.0; n * n];
    let k: Vec<_> = p.iter().map(|p| gp.kernel.eval(p, &[0.0])).collect();
    let mu: Vec<_> = k.iter().map(|k| k * 0.2 / 1.3).collect();
    for i in 0..n {
        for j in 0..=i {
            let v = gp.kernel.eval(&p[i], &p[j]) - k[i] * k[j] / 1.3;
            sigma[i * n + j] = if i == j { v.max(0.0) + jitter } else { v };
        }
    }
    let l = cholesky(&sigma, n).unwrap();
    let mut sum = 0.0;
    for row in bank.chunks_exact(n) {
        let mut old = f64::INFINITY; let mut new = f64::INFINITY;
        for i in 0..n {
            let mut value = mu[i];
            for j in 0..=i { value = l.l(i, j).mul_add(row[j], value); }
            if i < base.len() { old = old.min(value); } else { new = new.min(value); }
        }
        sum += (old - new).max(0.0);
    }
    sum / (bank.len() / n) as f64
}

#[test]
fn g1_cached_scores_match_full_joint_draws_under_the_same_jitter() {
    let gp = gp();
    let base = vec![vec![0.0], vec![0.7], vec![-0.6]]; // observed plus pending
    for candidates in [vec![vec![1.3]], vec![vec![1.3], vec![-1.1]]] {
        let bank = joint_normal_bank(2048, base.len() + candidates.len(), 37);
        let cache = CachedNoisyImprovement::new(&gp, &base, candidates.len(), &bank, 1e-7).unwrap();
        let got = cache.score(&candidates).unwrap();
        let expected = dense_score(&gp, &base, &candidates, &bank, 1e-7);
        assert!(got > 0.01);
        assert!((got - expected).abs() < 2e-11, "{got} != {expected}");
        assert_eq!(got.to_bits(), cache.score(&candidates).unwrap().to_bits());
    }
}

#[test]
fn g3_duplicate_columns_and_signed_zero_share_the_first_latent_draw() {
    let gp = gp();
    let bank = joint_normal_bank(256, 4, 41);
    let mut padded = Vec::new();
    for r in bank.chunks_exact(4) { padded.extend_from_slice(&[r[0], 1e9, r[1], r[2], -1e9, r[3]]); }
    let a = CachedNoisyImprovement::new(&gp, &[vec![0.0], vec![0.8]], 2, &bank, 1e-7).unwrap();
    let b = CachedNoisyImprovement::new(&gp, &[vec![0.0], vec![-0.0], vec![0.8]], 3, &padded, 1e-7).unwrap();
    assert_eq!(b.baseline_size(), 2);
    assert_eq!(a.score(&[vec![1.3], vec![-0.7]]).unwrap().to_bits(),
        b.score(&[vec![1.3], vec![1.3], vec![-0.7]]).unwrap().to_bits());
    assert_eq!(b.score(&[vec![-0.0], vec![0.8]]).unwrap(), 0.0);
}

#[test]
fn g3_fixed_sample_columns_make_batch_growth_monotone() {
    let gp = gp();
    let bank = joint_normal_bank(512, 3, 43);
    let cache = CachedNoisyImprovement::new(&gp, &[vec![0.0]], 2, &bank, 1e-7).unwrap();
    let a = cache.score(&[vec![0.8]]).unwrap();
    let b = cache.score(&[vec![0.8], vec![-0.9]]).unwrap();
    assert!(a > 0.01 && b >= a);
    assert_eq!(a.to_bits(), cache.score(&[vec![0.8], vec![0.8]]).unwrap().to_bits());
}

#[test]
fn g0_bad_banks_candidates_and_cancellation_do_not_poison_reuse() {
    let gp = gp(); let baseline = vec![vec![0.0]];
    for bank in [vec![], vec![0.0], vec![0.0, f64::NAN]] {
        assert!(matches!(CachedNoisyImprovement::new(&gp, &baseline, 1, &bank, 1e-7),
            Err(PreparedError::InvalidInput)));
    }
    let bank = joint_normal_bank(64, 2, 47);
    let cache = CachedNoisyImprovement::new(&gp, &baseline, 1, &bank, 1e-7).unwrap();
    let x = vec![vec![0.7]];
    let expected = cache.score(&x).unwrap();
    let mut checks = 0;
    assert_eq!(cache.score_controlled(&x, &mut || { checks += 1; checks < 7 }), Err(PreparedError::Cancelled));
    assert_eq!(cache.score(&x).unwrap().to_bits(), expected.to_bits());
    assert_eq!(cache.score(&[vec![0.7], vec![0.8]]), Err(PreparedError::InvalidInput));
    assert_eq!(cache.score(&[vec![f64::NAN]]), Err(PreparedError::InvalidInput));
}

#[test]
fn g5_cached_search_replays_and_excludes_all_reserved_designs() {
    let gp = gp(); let c = config();
    let base = vec![vec![0.0]]; let pending = vec![vec![0.4]];
    let a = c.propose_pending_cached(&gp, &base, &pending, 2, 5, 1e-7, &mut || true).unwrap();
    let b = c.propose_pending_cached(&gp, &base, &pending, 2, 5, 1e-7, &mut || true).unwrap();
    assert_eq!(a, b); assert_eq!(a.len(), 2); assert_ne!(a[0], a[1]);
    assert!(a.iter().all(|p| !base.contains(p) && !pending.contains(p)));
    assert!(a.iter().flatten().all(|v| (-2.0..=2.0).contains(v)));
    assert_eq!(c.propose_pending_cached(&gp, &base, &pending, 2, 5, 1e-7, &mut || false),
        Err(CachedProposalError::Cancelled));
    assert_eq!(c.propose_pending_cached(&gp, &base, &pending, 3, 5, 1e-7, &mut || true),
        Err(CachedProposalError::InvalidRequest));
    let mut checks = 0;
    // More calls than preparation needs: pause inside acquisition scoring.
    assert_eq!(c.propose_pending_cached(&gp, &base, &pending, 2, 5, 1e-7,
        &mut || { checks += 1; checks < 65 }), Err(CachedProposalError::Cancelled));
    assert_eq!(checks, 65, "no more posterior callbacks after cancellation");
}
