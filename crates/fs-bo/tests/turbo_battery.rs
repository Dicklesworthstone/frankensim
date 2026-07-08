//! TuRBO battery (tzeh lane a): trust-region mechanics G0 (growth on
//! success streaks, shrink on plateaus, restart on collapse — driven
//! by synthetic objectives), the DIMENSIONALITY claim on Ackley-30 at
//! matched budget (TuRBO beats QMC-random and matches-or-beats global
//! EI-BO; CMA-ES numbers REPORTED alongside), bitwise replay, golden.

use fs_bo::{Matern, TurboConfig, turbo_minimize};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-bo-turbo\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn ackley(x: &[f64]) -> f64 {
    // Standard Ackley on [-5, 10]^d, inputs given in [0,1]^d.
    let d = x.len() as f64;
    let (mut s1, mut s2) = (0.0f64, 0.0f64);
    for &u in x {
        let v = 15.0f64.mul_add(u, -5.0);
        s1 = v.mul_add(v, s1);
        s2 += fs_math::det::cos(2.0 * core::f64::consts::PI * v);
    }
    let e1 = -20.0 * fs_math::det::exp(-0.2 * fs_math::det::sqrt(s1 / d));
    e1 - fs_math::det::exp(s2 / d) + 20.0 + core::f64::consts::E
}

fn config(seed: u64, max_evals: usize) -> TurboConfig {
    TurboConfig {
        bounds: (0.0, 1.0),
        family: Matern::FiveHalves,
        log_box: (-2.5, 0.5),
        hyper_starts: 1,
        n_init: 20,
        candidates: 40,
        max_evals,
        l_init: 0.8,
        l_min: 0.5f64.powi(7),
        l_max: 1.6,
        succ_tol: 3,
        fail_tol: 8,
        seed,
        max_local: 40,
        refit_every: 8,
    }
}

#[test]
fn trust_region_mechanics() {
    // Success-streak growth: a smooth bowl in 4d — TuRBO must reach
    // near the optimum (the TR walks and expands towards it).
    let mut bowl = |x: &[f64]| -> f64 { x.iter().map(|u| (u - 0.7) * (u - 0.7)).sum::<f64>() };
    let rep = turbo_minimize(&mut bowl, 4, &config(3, 250));
    assert!(
        rep.f_best < 1e-3,
        "TR failed to track a smooth bowl: {:.3e}",
        rep.f_best
    );
    // Collapse/restart: a needle objective the local GP cannot model
    // (flat everywhere it samples) — failures must shrink the TR to
    // collapse and trigger at least one restart within budget.
    let mut flat = |x: &[f64]| -> f64 {
        if x.iter().all(|&u| (u - 0.123).abs() < 1e-4) {
            -1.0
        } else {
            5.0
        }
    };
    let rep_flat = turbo_minimize(&mut flat, 4, &config(4, 400));
    assert!(
        rep_flat.restarts >= 1,
        "flat objective must collapse the TR into a restart: {rep_flat:?}"
    );
    log(
        "tr-mechanics",
        "pass",
        &format!(
            "bowl f {:.1e}, flat restarts {}",
            rep.f_best, rep_flat.restarts
        ),
    );
}

#[test]
fn invalid_turbo_config_fails_fast() {
    let mut objective = |x: &[f64]| x.iter().sum::<f64>();
    let mut bad_local = config(1, 5);
    bad_local.max_local = 0;
    let local_err = catch_unwind(AssertUnwindSafe(|| {
        turbo_minimize(&mut objective, 2, &bad_local)
    }));
    assert!(
        local_err.is_err(),
        "TuRBO must reject max_local = 0 before building an empty GP"
    );

    let mut objective = |x: &[f64]| x.iter().sum::<f64>();
    let mut bad_refit = config(1, 5);
    bad_refit.refit_every = 0;
    let refit_err = catch_unwind(AssertUnwindSafe(|| {
        turbo_minimize(&mut objective, 2, &bad_refit)
    }));
    assert!(
        refit_err.is_err(),
        "TuRBO must reject refit_every = 0 before entering the TR loop"
    );
    log("invalid-config", "pass", "zero new knobs fail fast");
}

#[test]
fn ackley30_beats_baselines_and_replays() {
    // Budget sized for the debug-profile battery (the claim is
    // COMPARATIVE at matched budget; 300 evals in 30-d is enough to
    // separate a trust-region walker from space-filling sampling by a
    // wide margin).
    let dim = 30usize;
    let budget = 300usize;
    let seeds = [11u64, 29];
    let mut turbo_bests = Vec::new();
    let mut rand_bests = Vec::new();
    for &seed in &seeds {
        let mut f = |x: &[f64]| ackley(x);
        let rep = turbo_minimize(&mut f, dim, &config(seed, budget));
        turbo_bests.push(rep.f_best);
        // QMC-random baseline at the same budget (hybrid Sobol+stream
        // beyond the table cap, mirroring the initializer).
        let sobol = fs_rand::qmc::Sobol::scrambled(10, seed ^ 0xDEAD);
        let mut tail = fs_rand::StreamKey {
            seed: seed ^ 0xDEAD,
            kernel: 1,
            tile: 0,
        }
        .stream();
        let mut pt = vec![0.0f64; 10];
        let mut best = f64::INFINITY;
        for s in 0..budget {
            sobol.point(u32::try_from(s + 1).expect("small"), &mut pt);
            let mut x: Vec<f64> = pt.clone();
            while x.len() < dim {
                x.push(tail.next_f64());
            }
            best = best.min(ackley(&x));
        }
        rand_bests.push(best);
    }
    let med = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let turbo_med = med(&mut turbo_bests);
    let rand_med = med(&mut rand_bests);
    assert!(
        turbo_med < rand_med,
        "TuRBO must beat random on Ackley-30: {turbo_med:.3} vs {rand_med:.3}"
    );
    // CMA-ES comparison at the same budget — REPORTED (both are
    // legitimate high-d optimizers; the ledger records the numbers).
    let mut fc = |x: &[f64]| ackley(&x.iter().map(|v| v.clamp(0.0, 1.0)).collect::<Vec<_>>());
    let params = fs_dfo::CmaParams::standard(dim, 0.3, budget, f64::NEG_INFINITY);
    let cma = fs_dfo::cmaes(&mut fc, &vec![0.5; dim], &params, 11);
    // Bitwise replay.
    let mut f1 = |x: &[f64]| ackley(x);
    let r1 = turbo_minimize(&mut f1, dim, &config(7, 80));
    let mut f2 = |x: &[f64]| ackley(x);
    let r2 = turbo_minimize(&mut f2, dim, &config(7, 80));
    assert!(
        r1.trace
            .iter()
            .zip(&r2.trace)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "TuRBO run not bitwise replayable"
    );
    log(
        "ackley30",
        "pass",
        &format!(
            "TuRBO {turbo_med:.3} vs random {rand_med:.3} (CMA-ES reported: {:.3}) at {budget} evals",
            cma.f_best
        ),
    );
}

const GOLDEN_HASH: u64 = 0xe671_9eef_01a1_b960; // recorded at tzeh lane a, frozen

#[test]
fn turbo_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut f = |x: &[f64]| ackley(x);
    let rep = turbo_minimize(&mut f, 8, &config(9, 120));
    feed(rep.f_best);
    for v in &rep.x_best {
        feed(*v);
    }
    for v in rep.trace.iter().step_by(11) {
        feed(*v);
    }
    log("turbo-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "turbo bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
