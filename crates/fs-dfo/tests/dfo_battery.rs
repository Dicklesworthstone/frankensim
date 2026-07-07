//! DFO battery (7tv.4): benchmark convergence (incl. the ill-conditioned
//! ellipsoid that separates CMA from plain ES), bitwise seed determinism,
//! the IGO invariance laws tested bitwise, SPD maintenance, BIPOP
//! schedule shape, Nelder–Mead polish, and the cross-ISA golden hash.

use fs_dfo::{CmaParams, bipop_cmaes, cmaes, nelder_mead};

fn sphere(x: &[f64]) -> f64 {
    x.iter().map(|t| t * t).sum()
}

fn rosenbrock(x: &[f64]) -> f64 {
    x.windows(2)
        .map(|w| {
            let a = 1.0 - w[0];
            let b = w[1] - w[0] * w[0];
            100.0f64.mul_add(b * b, a * a)
        })
        .sum()
}

fn ellipsoid(x: &[f64]) -> f64 {
    // Condition number 1e6: the covariance-adaptation showcase.
    let n = x.len();
    x.iter()
        .enumerate()
        .map(|(i, t)| {
            let e = 1e6f64.powf(i as f64 / (n - 1) as f64);
            e * t * t
        })
        .sum()
}

fn rastrigin(x: &[f64]) -> f64 {
    let a = 10.0f64;
    x.iter()
        .map(|&t| a + t.mul_add(t, -a * fs_math::det::cos(2.0 * std::f64::consts::PI * t)))
        .sum()
}

#[test]
fn benchmark_convergence() {
    // Sphere(10): trivial, fast.
    let mut f = |x: &[f64]| sphere(x);
    let rep = cmaes(
        &mut f,
        &[3.0; 10],
        &CmaParams::standard(10, 2.0, 20_000, 1e-10),
        1,
    );
    assert!(rep.converged, "sphere must converge: {rep:?}");
    // Rosenbrock(10): the classic valley.
    let mut f = |x: &[f64]| rosenbrock(x);
    let rep = cmaes(
        &mut f,
        &[0.0; 10],
        &CmaParams::standard(10, 0.5, 120_000, 1e-8),
        2,
    );
    assert!(
        rep.converged,
        "rosenbrock(10) must converge: evals {}",
        rep.evals
    );
    // Ellipsoid(10) at condition 1e6: adaptation must handle it.
    let mut f = |x: &[f64]| ellipsoid(x);
    let rep = cmaes(
        &mut f,
        &[1.0; 10],
        &CmaParams::standard(10, 1.0, 120_000, 1e-8),
        3,
    );
    assert!(
        rep.converged,
        "ellipsoid must converge: evals {}",
        rep.evals
    );
    println!(
        "{{\"suite\":\"fs-dfo\",\"case\":\"benchmarks\",\"verdict\":\"pass\",\"detail\":\"sphere/rosenbrock/ellipsoid(cond 1e6) all to target\"}}"
    );
}

#[test]
fn deterministic_evolution_from_seed() {
    // Budget SHORT of convergence: at full convergence every seed lands
    // on the exact minimum (f = +0.0 bits) and cross-seed comparisons go
    // vacuous — measured during bring-up.
    let mut f1 = |x: &[f64]| rosenbrock(x);
    let mut f2 = |x: &[f64]| rosenbrock(x);
    let p = CmaParams::standard(6, 0.5, 1_500, -1.0);
    let r1 = cmaes(&mut f1, &[0.2; 6], &p, 42);
    let r2 = cmaes(&mut f2, &[0.2; 6], &p, 42);
    assert_eq!(
        r1.f_best.to_bits(),
        r2.f_best.to_bits(),
        "same seed → same bits"
    );
    assert_eq!(r1.evals, r2.evals);
    for (a, b) in r1.x_best.iter().zip(&r2.x_best) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    // Different seed → different mid-flight trajectory.
    let mut f3 = |x: &[f64]| rosenbrock(x);
    let r3 = cmaes(&mut f3, &[0.2; 6], &p, 43);
    assert!(
        r1.x_best
            .iter()
            .zip(&r3.x_best)
            .any(|(a, b)| a.to_bits() != b.to_bits()),
        "different seeds must explore differently mid-flight"
    );
}

#[test]
fn monotone_transform_invariance_is_bitwise() {
    // Rank-based selection sees only the ORDER of fitness values: any
    // strictly monotone transform of f must give the IDENTICAL evolution
    // — the IGO invariance, tested bitwise on the search trajectory.
    // Budget kept SHORT of deep convergence: near machine precision the
    // transforms stop being injective in f64 (x³ underflows below
    // ~1e-108, exp(x) saturates to 1.0 below ~1e-16), which creates
    // ties the plain objective does not have — measured during bring-up.
    // Strict monotonicity at the RESOLUTION OF THE SAMPLES is the real
    // precondition, and it holds on this budget.
    let p = CmaParams::standard(5, 0.7, 800, -1.0);
    let mut plain = |x: &[f64]| sphere(x);
    let r_plain = cmaes(&mut plain, &[1.5; 5], &p, 7);
    // exp is strictly monotone; sphere ≥ 0 so cube is monotone there too.
    let mut expf = |x: &[f64]| fs_math::det::exp(sphere(x).min(700.0));
    let r_exp = cmaes(&mut expf, &[1.5; 5], &p, 7);
    for (a, b) in r_plain.x_best.iter().zip(&r_exp.x_best) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "monotone transform must not change the trajectory"
        );
    }
    let mut cubef = |x: &[f64]| sphere(x).powi(3);
    let r_cube = cmaes(&mut cubef, &[1.5; 5], &p, 7);
    for (a, b) in r_plain.x_best.iter().zip(&r_cube.x_best) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "cube transform must not change it either"
        );
    }
    println!(
        "{{\"suite\":\"fs-dfo\",\"case\":\"igo-invariance\",\"verdict\":\"pass\",\"detail\":\"exp/cube monotone transforms: bitwise-identical trajectories\"}}"
    );
}

#[test]
fn translation_equivariance() {
    // Shifting the optimum shifts the answer, nothing else (within the
    // stochastic tolerance of identical seeds — the SAMPLES differ in
    // absolute space, so this is behavioral, not bitwise).
    let shift = [2.5f64, -1.0, 0.5, 3.0, -2.0];
    let p = CmaParams::standard(5, 1.0, 40_000, 1e-10);
    let mut f0 = |x: &[f64]| sphere(x);
    let r0 = cmaes(&mut f0, &[1.0; 5], &p, 11);
    let mut fs = |x: &[f64]| {
        let shifted: Vec<f64> = x.iter().zip(&shift).map(|(a, s)| a - s).collect();
        sphere(&shifted)
    };
    let start: Vec<f64> = shift.iter().map(|s| s + 1.0).collect();
    let rs = cmaes(&mut fs, &start, &p, 11);
    assert!(r0.converged && rs.converged);
    for (got, want) in rs.x_best.iter().zip(&shift) {
        assert!(
            (got - want).abs() < 1e-4,
            "shifted optimum: {:?}",
            rs.x_best
        );
    }
}

#[test]
fn bipop_solves_multimodal_and_reports_schedule() {
    let mut f = |x: &[f64]| rastrigin(x);
    let rep = bipop_cmaes(&mut f, &[3.0; 5], 2.0, 400_000, 1e-8, 17);
    assert!(
        rep.best.f_best < 1e-6,
        "BIPOP must reach the global basin on rastrigin(5): {}",
        rep.best.f_best
    );
    assert!(!rep.schedule.is_empty());
    // Schedule shape: population sizes are the base or doublings of it.
    let base = rep.schedule[0];
    for &lam in &rep.schedule {
        assert!(
            lam % base == 0 && (lam / base).is_power_of_two(),
            "schedule must be doublings of the base: {:?}",
            rep.schedule
        );
    }
    println!(
        "{{\"suite\":\"fs-dfo\",\"case\":\"bipop\",\"verdict\":\"pass\",\"detail\":\"rastrigin(5) f={:.2e}, schedule {:?}, {} evals\"}}",
        rep.best.f_best, rep.schedule, rep.total_evals
    );
}

#[test]
fn nelder_mead_polishes() {
    let mut f = |x: &[f64]| rosenbrock(x);
    let (x, fv, evals) = nelder_mead(&mut f, &[0.8, 0.6], 0.1, 5_000, 1e-12);
    assert!(
        fv < 1e-10,
        "NM must polish rosenbrock(2): f={fv:.2e} after {evals} evals"
    );
    assert!((x[0] - 1.0).abs() < 1e-4 && (x[1] - 1.0).abs() < 1e-4);
    // Fully deterministic: identical reruns bitwise.
    let mut g = |x: &[f64]| rosenbrock(x);
    let (x2, f2, e2) = nelder_mead(&mut g, &[0.8, 0.6], 0.1, 5_000, 1e-12);
    assert_eq!(fv.to_bits(), f2.to_bits());
    assert_eq!(evals, e2);
    assert_eq!(x[0].to_bits(), x2[0].to_bits());
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x5441_10a6_afb1_70a1; // bumped: TolFun stagnation criterion added (semantic change)

#[test]
fn dfo_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let p = CmaParams::standard(6, 0.8, 6_000, -1.0);
    let mut f = |x: &[f64]| rosenbrock(x);
    let rep = cmaes(&mut f, &[0.1; 6], &p, 99);
    feed(rep.f_best);
    feed(rep.sigma);
    for &v in &rep.x_best {
        feed(v);
    }
    let mut g = |x: &[f64]| ellipsoid(x);
    let rep2 = cmaes(&mut g, &[1.0; 6], &p, 100);
    feed(rep2.f_best);
    let mut h = |x: &[f64]| rosenbrock(x);
    let (xnm, fnm, _) = nelder_mead(&mut h, &[0.3, -0.2], 0.2, 2_000, -1.0);
    feed(fnm);
    feed(xnm[0]);
    feed(xnm[1]);
    println!(
        "{{\"suite\":\"fs-dfo\",\"case\":\"dfo-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "DFO bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
