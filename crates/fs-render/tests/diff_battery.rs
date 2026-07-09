//! fs-render differentiable battery (bead qfx.5, smoke tier):
//! edge-aware gradients vs FD of the render, the naive-autodiff
//! negative control, quadrature-bias characterization, inverse
//! rendering, the combined appearance+physics fixture, and bitwise
//! replay.

use fs_render::diff::{loss_and_grad, render, render_grad, NPARAMS, RenderCfg};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// The two-sphere fixture: overlapping blended spheres, off-center.
fn theta0() -> [f64; NPARAMS] {
    [0.38, 0.45, 0.0, 0.22, 0.62, 0.58, 0.1, 0.17, 0.08]
}

/// dr-001: the edge-aware gradient of the full-image L2 loss matches
/// central FD of the RENDERED loss for every parameter — silhouette
/// (centers, radii) and shading/blend derivatives together.
#[test]
fn dr_001_gradient_matches_fd() {
    let cfg = RenderCfg::default();
    let th = theta0();
    // Target: render at a shifted θ so the loss is not at a minimum.
    let mut tht = th;
    tht[0] += 0.06;
    tht[3] -= 0.03;
    tht[8] += 0.02;
    let target = render(&tht, cfg);
    let (_, grad) = loss_and_grad(&th, &target, cfg);
    let h = 1e-5;
    let mut worst = 0.0f64;
    let mut details = String::new();
    for k in 0..NPARAMS {
        let (mut tp, mut tm) = (th, th);
        tp[k] += h;
        tm[k] -= h;
        let lp: f64 = render(&tp, cfg)
            .iter()
            .zip(&target)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            / (cfg.res * cfg.res) as f64;
        let lm: f64 = render(&tm, cfg)
            .iter()
            .zip(&target)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            / (cfg.res * cfg.res) as f64;
        let fd = (lp - lm) / (2.0 * h);
        let scale = grad[k].abs().max(fd.abs()).max(1e-8);
        let rel = (grad[k] - fd).abs() / scale;
        worst = worst.max(rel);
        use std::fmt::Write as _;
        let _ = write!(details, "p{k}:{rel:.1e} ");
    }
    verdict(
        "dr-001-gradient-vs-fd",
        worst < 1e-4,
        &format!("edge-aware grad vs central FD, worst rel {worst:.2e} ({details})"),
    );
}

/// dr-002: the NEGATIVE CONTROL — freezing the crossings (what naive
/// pointwise autodiff computes) must be measurably WRONG on the
/// silhouette-dominated parameters. This is the bead's core claim:
/// visibility discontinuities carry a boundary term autodiff misses.
#[test]
fn dr_002_naive_autodiff_is_silently_wrong() {
    let cfg = RenderCfg::default();
    let th = theta0();
    let mut tht = th;
    tht[0] += 0.06;
    let target = render(&tht, cfg);
    let scale = 1.0 / (cfg.res * cfg.res) as f64;
    // Edge-aware and naive gradients of the same loss (d/d cx).
    let full = render_grad(&th, cfg, true);
    let naive = render_grad(&th, cfg, false);
    let (mut g_full, mut g_naive) = (0.0f64, 0.0f64);
    for i in 0..full.len() {
        let r = full[i].re - target[i];
        g_full += 2.0 * r * full[i].eps[0] * scale;
        g_naive += 2.0 * r * naive[i].eps[0] * scale;
    }
    // FD reference for the x-center parameter (silhouette-dominated).
    let h = 1e-5;
    let (mut tp, mut tm) = (th, th);
    tp[0] += h;
    tm[0] -= h;
    let lp: f64 = render(&tp, cfg)
        .iter()
        .zip(&target)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        * scale;
    let lm: f64 = render(&tm, cfg)
        .iter()
        .zip(&target)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        * scale;
    let fd = (lp - lm) / (2.0 * h);
    let err_full = (g_full - fd).abs() / fd.abs().max(1e-12);
    let err_naive = (g_naive - fd).abs() / fd.abs().max(1e-12);
    verdict(
        "dr-002-naive-negative-control",
        err_full < 1e-4 && err_naive > 100.0 * err_full,
        &format!(
            "d(loss)/d(cx): FD {fd:.6e}, edge-aware rel err {err_full:.2e}, NAIVE rel err {err_naive:.2e} — the boundary term is {:.0}x of the honest error budget",
            err_naive / err_full.max(1e-300)
        ),
    );
}

/// dr-003: bias discipline — the deterministic quadrature's
/// discretization error against a much finer reference shrinks with
/// the coarse sampling knobs (no variance; bias measured, not vibes).
#[test]
fn dr_003_quadrature_bias_shrinks() {
    let th = theta0();
    let fine = render(
        &th,
        RenderCfg {
            res: 32,
            subrows: 16,
            xsamples: 16,
        },
    );
    let err_at = |subrows: usize| -> f64 {
        let img = render(
            &th,
            RenderCfg {
                res: 32,
                subrows,
                xsamples: 4,
            },
        );
        img.iter()
            .zip(&fine)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max)
    };
    let e1 = err_at(1);
    let e2 = err_at(2);
    let e4 = err_at(4);
    let e8 = err_at(8);
    verdict(
        "dr-003-bias-shrinks",
        e2 < e1 && e4 < e2 && e8 < e4 && e8 < 0.02,
        &format!("max-pixel bias vs fine reference: subrows 1/2/4/8 -> {e1:.2e}/{e2:.2e}/{e4:.2e}/{e8:.2e} (monotone shrink)"),
    );
}

/// dr-004: inverse rendering — recover sphere-1 position and radius
/// from a target image by gradient descent with backtracking, using
/// ONLY the edge-aware gradient.
#[test]
fn dr_004_inverse_rendering_recovers_shape() {
    let cfg = RenderCfg::default();
    let truth = theta0();
    let target = render(&truth, cfg);
    let mut th = truth;
    th[0] += 0.05;
    th[1] -= 0.04;
    th[3] += 0.03;
    let free = [0usize, 1, 3];
    let mut step = 0.5f64;
    let (mut loss, mut grad) = loss_and_grad(&th, &target, cfg);
    for _ in 0..60 {
        // Backtracking line search along the masked gradient.
        let mut trial = th;
        loop {
            for &k in &free {
                trial[k] = th[k] - step * grad[k];
            }
            let (lt, gt) = loss_and_grad(&trial, &target, cfg);
            if lt < loss {
                th = trial;
                loss = lt;
                grad = gt;
                step *= 1.3;
                break;
            }
            step *= 0.5;
            assert!(step > 1e-12, "line search collapsed at loss {loss:.3e}");
        }
        if loss < 1e-12 {
            break;
        }
    }
    let dev = free
        .iter()
        .map(|&k| (th[k] - truth[k]).abs())
        .fold(0.0f64, f64::max);
    verdict(
        "dr-004-inverse-rendering",
        loss < 1e-10 && dev < 1e-4,
        &format!("recovered (cx, cy, r1): worst |dev| {dev:.2e}, final loss {loss:.2e}"),
    );
}

/// dr-005: aesthetic + physics — one combined objective (image L2 +
/// a volume-budget penalty) optimized end-to-end through the SHARED
/// gradient path; both terms must drop.
#[test]
fn dr_005_combined_appearance_physics() {
    let cfg = RenderCfg::default();
    let truth = theta0();
    let target = render(&truth, cfg);
    // Physics term: (sum of sphere volumes − budget)², budget set at
    // the true radii so the optimum is compatible with the image term.
    let vol = |th: &[f64]| -> (f64, [f64; NPARAMS]) {
        let c = 4.0 / 3.0 * std::f64::consts::PI;
        let v = c * (th[3].powi(3) + th[7].powi(3));
        let mut g = [0.0f64; NPARAMS];
        g[3] = 3.0 * c * th[3] * th[3];
        g[7] = 3.0 * c * th[7] * th[7];
        (v, g)
    };
    let (vstar, _) = vol(&truth);
    let w = 0.5f64;
    let objective = |th: &[f64], target: &[f64]| -> (f64, [f64; NPARAMS]) {
        let (li, gi) = loss_and_grad(th, target, cfg);
        let (v, gv) = vol(th);
        let mut g = gi;
        for k in 0..NPARAMS {
            g[k] += w * 2.0 * (v - vstar) * gv[k];
        }
        (li + w * (v - vstar) * (v - vstar), g)
    };
    let mut th = truth;
    th[0] += 0.04;
    th[3] += 0.04;
    th[7] -= 0.02;
    let free = [0usize, 3, 7];
    let (l0, _) = objective(&th, &target);
    let (mut loss, mut grad) = objective(&th, &target);
    let mut step = 0.2f64;
    for _ in 0..60 {
        let mut trial = th;
        loop {
            for &k in &free {
                trial[k] = th[k] - step * grad[k];
            }
            let (lt, gt) = objective(&trial, &target);
            if lt < loss {
                th = trial;
                loss = lt;
                grad = gt;
                step *= 1.3;
                break;
            }
            step *= 0.5;
            assert!(step > 1e-12, "combined line search collapsed at {loss:.3e}");
        }
        if loss < 1e-11 {
            break;
        }
    }
    let (vfinal, _) = vol(&th);
    verdict(
        "dr-005-combined-objective",
        loss < 1e-8 && loss < 1e-3 * l0 && (vfinal - vstar).abs() < 1e-3,
        &format!(
            "combined loss {l0:.3e} -> {loss:.3e}; volume {vfinal:.5} vs budget {vstar:.5} — appearance and physics share one gradient path"
        ),
    );
}

/// dr-006: determinism — bitwise replay of both the render and the
/// gradient.
#[test]
fn dr_006_bitwise_replay() {
    let cfg = RenderCfg::default();
    let th = theta0();
    let a = render(&th, cfg);
    let b = render(&th, cfg);
    let bit_img = a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits());
    let ga = render_grad(&th, cfg, true);
    let gb = render_grad(&th, cfg, true);
    let bit_grad = ga.iter().zip(&gb).all(|(x, y)| {
        x.re.to_bits() == y.re.to_bits()
            && x.eps
                .iter()
                .zip(&y.eps)
                .all(|(u, v)| u.to_bits() == v.to_bits())
    });
    verdict(
        "dr-006-bitwise-replay",
        bit_img && bit_grad,
        "render and edge-aware gradient replay bitwise",
    );
}
