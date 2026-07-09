//! NSGA-III battery (vcia many-objective lane): Das–Dennis direction
//! counts + simplex membership; DTLZ2(3) convergence to the known
//! unit-sphere-octant front with reference-direction COVERAGE
//! (diversity measured by association); the many-objective claim at
//! m = 5 — NSGA-III beats NSGA-II on MC-estimated hypervolume at
//! matched budget; bitwise replay; golden.

use fs_dfo::{NsgaParams, das_dennis, hypervolume, mc_hypervolume, nsga2, nsga3};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-dfo-nsga3\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

/// DTLZ2 with m objectives, n = m − 1 + k variables in [0,1].
fn dtlz2(x: &[f64], m: usize) -> Vec<f64> {
    let k = x.len() - (m - 1);
    let g: f64 = x[m - 1..]
        .iter()
        .map(|v| (v - 0.5) * (v - 0.5))
        .sum::<f64>();
    let _ = k;
    let half_pi = core::f64::consts::FRAC_PI_2;
    (0..m)
        .map(|i| {
            let mut f = 1.0 + g;
            for &xj in &x[..m - 1 - i] {
                f *= fs_math::det::cos(xj * half_pi);
            }
            if i > 0 {
                f *= fs_math::det::sin(x[m - 1 - i] * half_pi);
            }
            f
        })
        .collect()
}

#[test]
fn das_dennis_counts_and_simplex() {
    // C(p+m−1, m−1): m=3, p=12 → C(14,2) = 91; m=5, p=4 → C(8,4) = 70.
    let d3 = das_dennis(3, 12);
    assert_eq!(d3.len(), 91);
    let d5 = das_dennis(5, 4);
    assert_eq!(d5.len(), 70);
    for dir in d3.iter().chain(&d5) {
        let s: f64 = dir.iter().sum();
        assert!(
            (s - 1.0).abs() < 1e-12,
            "direction off the simplex: {dir:?}"
        );
        assert!(dir.iter().all(|&v| v >= 0.0));
    }
    log("das-dennis", "pass", "91 @ (3,12), 70 @ (5,4), on-simplex");
}

#[test]
fn dtlz2_m3_convergence_and_coverage() {
    let m = 3usize;
    let dirs = das_dennis(m, 12);
    // Standard DTLZ2 budgets (~250 generations); 150 left the worst
    // straggler at 0.0515 against the 0.05 gate.
    let params = NsgaParams {
        pop: 92,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: 17,
    };
    let mut f = |x: &[f64]| dtlz2(x, m);
    let front = nsga3(&mut f, 7, (0.0, 1.0), &dirs, &params);
    // Convergence: the true front is ‖f‖₂ = 1.
    let mut worst_norm = 0.0f64;
    for ind in &front {
        let n2: f64 = ind.f.iter().map(|v| v * v).sum();
        worst_norm = worst_norm.max((fs_math::det::sqrt(n2) - 1.0).abs());
    }
    assert!(
        worst_norm < 0.05,
        "DTLZ2 front not converged: worst | ||f||-1 | = {worst_norm:.4}"
    );
    // Coverage: fraction of reference directions holding an associate.
    let covered = {
        let mut hit = vec![false; dirs.len()];
        for ind in &front {
            let mut best = (0usize, f64::INFINITY);
            for (k, dir) in dirs.iter().enumerate() {
                let dd: f64 = dir.iter().map(|d| d * d).sum();
                let t: f64 = ind.f.iter().zip(dir).map(|(a, b)| a * b).sum::<f64>() / dd;
                let d2: f64 = ind
                    .f
                    .iter()
                    .zip(dir)
                    .map(|(a, b)| {
                        let r = t.mul_add(-b, *a);
                        r * r
                    })
                    .sum();
                if d2 < best.1 {
                    best = (k, d2);
                }
            }
            hit[best.0] = true;
        }
        hit.iter().filter(|&&h| h).count() as f64 / dirs.len() as f64
    };
    assert!(
        covered > 0.6,
        "reference-direction coverage too low: {covered:.2}"
    );
    log(
        "dtlz2-m3",
        "pass",
        &format!(
            "worst norm dev {worst_norm:.4}, coverage {covered:.2}, front {}",
            front.len()
        ),
    );
}

#[test]
fn many_objective_m5_beats_nsga2_on_hv() {
    let m = 5usize;
    let dirs = das_dennis(m, 4);
    let params = NsgaParams {
        pop: 70,
        generations: 120,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 9.0,
        seed: 23,
    };
    let mut f3 = |x: &[f64]| dtlz2(x, m);
    let front3 = nsga3(&mut f3, 9, (0.0, 1.0), &dirs, &params);
    let mut f2 = |x: &[f64]| dtlz2(x, m);
    let front2 = nsga2(&mut f2, 9, (0.0, 1.0), &params);
    let reference = vec![1.5f64; m];
    let pts3: Vec<Vec<f64>> = front3.iter().map(|i| i.f.clone()).collect();
    let pts2: Vec<Vec<f64>> = front2.iter().map(|i| i.f.clone()).collect();
    let (hv3, _) = mc_hypervolume(&pts3, &reference, 200_000, 99);
    let (hv2, _) = mc_hypervolume(&pts2, &reference, 200_000, 99);
    assert!(
        hv3 > hv2,
        "NSGA-III should beat NSGA-II at m=5: {hv3:.4} vs {hv2:.4}"
    );
    // Bitwise replay of NSGA-III.
    let mut fr = |x: &[f64]| dtlz2(x, m);
    let ra = nsga3(&mut fr, 9, (0.0, 1.0), &dirs, &params);
    let mut fr2 = |x: &[f64]| dtlz2(x, m);
    let rb = nsga3(&mut fr2, 9, (0.0, 1.0), &dirs, &params);
    assert_eq!(ra.len(), rb.len());
    for (p, q) in ra.iter().zip(&rb) {
        assert!(
            p.f.iter()
                .zip(&q.f)
                .all(|(u, v)| u.to_bits() == v.to_bits())
        );
    }
    log(
        "m5-vs-nsga2",
        "pass",
        &format!("HV nsga3 {hv3:.4} vs nsga2 {hv2:.4} at matched budget, replay bitwise"),
    );
}

const GOLDEN_HASH: u64 = 0xd912_6c49_f1b1_6897; // recorded at vcia NSGA-III lane, frozen

#[test]
fn nsga3_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let dirs = das_dennis(3, 6);
    for d in dirs.iter().step_by(5) {
        for v in d {
            feed(*v);
        }
    }
    let params = NsgaParams {
        pop: 28,
        generations: 30,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 0.2,
        seed: 3,
    };
    let mut f = |x: &[f64]| dtlz2(x, 3);
    let front = nsga3(&mut f, 5, (0.0, 1.0), &dirs, &params);
    for ind in front.iter().take(10) {
        for v in &ind.f {
            feed(*v);
        }
    }
    log("nsga3-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "nsga3 bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn moead_zdt1_and_dtlz2_competitive() {
    use fs_dfo::{MoeadParams, moead};
    // ZDT1 convergence + spread (the NSGA-II gates, decomposition path).
    fn zdt1(x: &[f64]) -> Vec<f64> {
        let f1 = x[0];
        let g = 1.0 + 9.0 * x[1..].iter().sum::<f64>() / (x.len() - 1) as f64;
        vec![f1, g * (1.0 - fs_math::det::sqrt(f1 / g))]
    }
    let weights2 = das_dennis(2, 79); // 80 subproblems
    let params = MoeadParams {
        neighbors: 12,
        max_replace: 2,
        generations: 220,
        eta_c: 20.0,
        eta_m: 20.0,
        p_mut: 1.0 / 8.0,
        seed: 29,
    };
    let mut f = |x: &[f64]| zdt1(x);
    let front = moead(&mut f, 8, (0.0, 1.0), &weights2, &params);
    let mean_gap: f64 = front
        .iter()
        .map(|ind| (ind.f[1] - (1.0 - fs_math::det::sqrt(ind.f[0]))).abs())
        .sum::<f64>()
        / front.len() as f64;
    assert!(mean_gap < 0.05, "MOEA/D ZDT1 not converged: {mean_gap:.4}");
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for ind in &front {
        lo = lo.min(ind.f[0]);
        hi = hi.max(ind.f[0]);
    }
    assert!(hi - lo > 0.7, "MOEA/D diversity collapsed: [{lo:.3},{hi:.3}]");
    // DTLZ2(m=3): competitive with NSGA-III on hypervolume at matched
    // budget (within 10% — both are legitimate; numbers ledgered).
    let m = 3usize;
    let dirs = das_dennis(m, 12);
    let params3 = MoeadParams {
        neighbors: 15,
        max_replace: 2,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: 31,
    };
    let mut fd = |x: &[f64]| dtlz2(x, m);
    let front_md = moead(&mut fd, 7, (0.0, 1.0), &dirs, &params3);
    let nsga_params = NsgaParams {
        pop: 92,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: 31,
    };
    let mut fn3 = |x: &[f64]| dtlz2(x, m);
    let front_n3 = nsga3(&mut fn3, 7, (0.0, 1.0), &dirs, &nsga_params);
    let reference = vec![1.5f64; m];
    let pts_md: Vec<Vec<f64>> = front_md.iter().map(|i| i.f.clone()).collect();
    let pts_n3: Vec<Vec<f64>> = front_n3.iter().map(|i| i.f.clone()).collect();
    let hv_md = hypervolume(&pts_md, &reference);
    let hv_n3 = hypervolume(&pts_n3, &reference);
    assert!(
        hv_md > 0.9 * hv_n3,
        "MOEA/D should be competitive with NSGA-III: {hv_md:.4} vs {hv_n3:.4}"
    );
    // Bitwise replay.
    let mut fr = |x: &[f64]| zdt1(x);
    let ra = moead(&mut fr, 8, (0.0, 1.0), &weights2, &params);
    let mut fr2 = |x: &[f64]| zdt1(x);
    let rb = moead(&mut fr2, 8, (0.0, 1.0), &weights2, &params);
    assert_eq!(ra.len(), rb.len());
    for (p, q) in ra.iter().zip(&rb) {
        assert!(p.f.iter().zip(&q.f).all(|(u, v)| u.to_bits() == v.to_bits()));
    }
    log(
        "moead",
        "pass",
        &format!(
            "ZDT1 gap {mean_gap:.4} spread {:.2}; DTLZ2 HV {hv_md:.4} vs NSGA-III {hv_n3:.4}",
            hi - lo
        ),
    );
}
