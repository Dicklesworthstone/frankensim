//! Wasserstein-DRO oracle battery (complements dro_battery.rs, which
//! covers kink recovery + saturation + contract guards): closed-form
//! endpoints, monotonicity in ρ, STRONG DUALITY machine-checked
//! against an exact tiny-LP enumeration oracle, the robust-decision
//! shift demo, and the golden hash.

use fs_dfo::wasserstein_worst_case;
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-dfo-dro-oracle\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

/// Exact primal LP oracle by basic-solution enumeration: maximize
/// Σ πᵢⱼ ℓⱼ s.t. Σⱼ πᵢⱼ = 1/n, π ≥ 0, ⟨π, C⟩ ≤ ρ. Basic solutions
/// assign each row to one destination, except at most ONE row that
/// splits between two destinations to bind the budget.
fn primal_oracle(losses: &[f64], costs: &[f64], n: usize, rho: f64) -> f64 {
    let m = losses.len();
    let wn = 1.0 / n as f64;
    let mut best = f64::NEG_INFINITY;
    let mut total = 1usize;
    for _ in 0..n {
        total *= m;
    }
    for code in 0..total {
        let mut cc = code;
        let mut cost = 0.0f64;
        let mut val = 0.0f64;
        let mut assign = vec![0usize; n];
        for (i, a) in assign.iter_mut().enumerate() {
            *a = cc % m;
            cc /= m;
            cost += wn * costs[i * m + *a];
            val += wn * losses[*a];
        }
        if cost <= rho + 1e-12 {
            best = best.max(val);
        }
        for i in 0..n {
            let j0 = assign[i];
            for k in 0..m {
                if k == j0 {
                    continue;
                }
                let dc = wn * (costs[i * m + k] - costs[i * m + j0]);
                if dc.abs() < 1e-15 {
                    continue;
                }
                let t = (rho - cost) / dc;
                if (0.0..=1.0).contains(&t) {
                    best = best.max(t.mul_add(wn * (losses[k] - losses[j0]), val));
                }
            }
        }
    }
    best
}

fn q_expectation(losses: &[f64], q: &[f64]) -> f64 {
    losses.iter().zip(q).map(|(loss, mass)| loss * mass).sum()
}

#[test]
fn closed_form_endpoints_and_monotonicity() {
    let mut s = StreamKey {
        seed: 131,
        kernel: 0x0D20,
        tile: 1,
    }
    .stream();
    let n = 4usize;
    let m = 5usize;
    let losses: Vec<f64> = (0..m).map(|_| s.next_f64() * 3.0).collect();
    let mut costs = vec![0.0f64; n * m];
    for i in 0..n {
        for j in 0..m {
            costs[i * m + j] = if i == j { 0.0 } else { 0.2 + s.next_f64() };
        }
    }
    let rep0 = wasserstein_worst_case(&losses, &costs, n, 0.0);
    let mean0: f64 = (0..n).map(|i| losses[i]).sum::<f64>() / n as f64;
    assert!(
        (rep0.worst_case - mean0).abs() < 1e-9,
        "rho=0 must equal the empirical mean: {} vs {mean0}",
        rep0.worst_case
    );
    let lmax = losses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let diam = costs.iter().copied().fold(0.0f64, f64::max);
    let rep_big = wasserstein_worst_case(&losses, &costs, n, diam * 2.0);
    assert!(
        (rep_big.worst_case - lmax).abs() < 1e-9,
        "large rho must reach the max loss: {} vs {lmax}",
        rep_big.worst_case
    );
    let mut prev = rep0.worst_case;
    for k in 1..=10 {
        let rho = diam * 0.04 * f64::from(k);
        let r = wasserstein_worst_case(&losses, &costs, n, rho);
        assert!(r.worst_case >= prev - 1e-9, "must be monotone in rho");
        prev = r.worst_case;
    }
    log(
        "endpoints",
        "pass",
        &format!("mean {mean0:.4}, max {lmax:.4}, monotone"),
    );
}

#[test]
fn strong_duality_vs_lp_oracle() {
    let mut worst_gap = 0.0f64;
    for inst in 0..20u32 {
        let mut s = StreamKey {
            seed: 132,
            kernel: 0x0D21,
            tile: inst,
        }
        .stream();
        let n = 3usize;
        let m = 3usize;
        let losses: Vec<f64> = (0..m).map(|_| s.next_f64() * 2.0).collect();
        let mut costs = vec![0.0f64; n * m];
        for i in 0..n {
            for j in 0..m {
                costs[i * m + j] = if i == j { 0.0 } else { 0.1 + s.next_f64() };
            }
        }
        let rho = 0.15 * (1.0 + s.next_f64());
        let report = wasserstein_worst_case(&losses, &costs, n, rho);
        let dual = report.worst_case;
        let primal = primal_oracle(&losses, &costs, n, rho);
        worst_gap = worst_gap.max((dual - primal).abs());
        let q_value = q_expectation(&losses, &report.q);
        assert!(
            (q_value - dual).abs() < 1e-8,
            "recovered q must realize the dual value for inst {inst}: q_value {q_value:.17e}, \
             dual {dual:.17e}, primal {primal:.17e}, rho {rho:.17e}, lambda {:.17e}, losses \
             {losses:?}, costs {costs:?}, q {:?}",
            report.lambda,
            report.q
        );
        assert!(
            (report.q.iter().sum::<f64>() - 1.0).abs() < 1e-12
                && report.q.iter().all(|mass| *mass >= -1e-12),
            "recovered q must be a probability distribution: {:?}",
            report.q
        );
    }
    assert!(
        worst_gap < 1e-8,
        "duality gap vs the LP oracle: {worst_gap:.3e}"
    );
    log(
        "strong-duality",
        "pass",
        &format!("worst gap {worst_gap:.1e} over 20 instances"),
    );
}

#[test]
fn robust_decision_shifts_conservatively() {
    let support = [0.0f64, 1.0, 4.0];
    let samples = [0.0f64, 0.0, 1.0, 1.0];
    let n = samples.len();
    let m = support.len();
    let mut costs = vec![0.0f64; n * m];
    for (i, &xi) in samples.iter().enumerate() {
        for (j, &zj) in support.iter().enumerate() {
            costs[i * m + j] = (xi - zj).abs();
        }
    }
    let objective = |x: f64, rho: f64| -> f64 {
        let losses: Vec<f64> = support.iter().map(|&z| (x - z) * (x - z)).collect();
        wasserstein_worst_case(&losses, &costs, n, rho).worst_case
    };
    let minimize_1d = |rho: f64| -> f64 {
        let (mut a, mut b) = (-1.0f64, 5.0f64);
        let phi = 0.618_033_988_749_894_9f64;
        let mut c = b - phi * (b - a);
        let mut d = a + phi * (b - a);
        let mut fc = objective(c, rho);
        let mut fd = objective(d, rho);
        for _ in 0..120 {
            if fc < fd {
                b = d;
                d = c;
                fd = fc;
                c = b - phi * (b - a);
                fc = objective(c, rho);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + phi * (b - a);
                fd = objective(d, rho);
            }
        }
        f64::midpoint(a, b)
    };
    let x_emp = minimize_1d(0.0);
    let x_dro = minimize_1d(0.8);
    assert!(
        (x_emp - 0.5).abs() < 0.05,
        "empirical minimizer should sit at the sample mean: {x_emp:.3}"
    );
    assert!(
        x_dro > x_emp + 0.2,
        "DRO minimizer must shift toward the adversarial tail: {x_dro:.3} vs {x_emp:.3}"
    );
    log(
        "robust-shift",
        "pass",
        &format!("x* {x_emp:.3} -> {x_dro:.3}"),
    );
}

const GOLDEN_HASH: u64 = 0xd21c_d092_b4a5_ba98;

#[test]
fn dro_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut s = StreamKey {
        seed: 133,
        kernel: 0x0D22,
        tile: 0,
    }
    .stream();
    let (n, m) = (4usize, 4usize);
    let losses: Vec<f64> = (0..m).map(|_| s.next_f64() * 2.0).collect();
    let mut costs = vec![0.0f64; n * m];
    for i in 0..n {
        for j in 0..m {
            costs[i * m + j] = if i == j { 0.0 } else { 0.1 + s.next_f64() };
        }
    }
    for &rho in &[0.0f64, 0.1, 0.3, 0.9] {
        let rep = wasserstein_worst_case(&losses, &costs, n, rho);
        feed(rep.worst_case);
        feed(rep.lambda);
        for v in &rep.q {
            feed(*v);
        }
    }
    log("dro-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "dro bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
