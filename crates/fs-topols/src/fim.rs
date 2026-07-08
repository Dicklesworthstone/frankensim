//! Fast-iterative-method redistancing: restore the signed-distance
//! property |∇φ| = 1 without moving the zero level set (the audit
//! quantifies how well "without" holds).
//!
//! Interface cells are FROZEN at values reconstructed from linear
//! zero-crossings of the input field (the interface is data, not an
//! unknown); everything else relaxes through Godunov upwind eikonal
//! updates on an active list (FIM — the parallel-friendly cousin of
//! fast marching; the list converges in O(band) sweeps here and the
//! update order does not affect the fixed point, which is what makes
//! the method many-core-safe).
//!
//! Audits: [`RedistanceAudit`] carries the sampled zero-level-set
//! Hausdorff drift and |∇φ|−1 statistics — the redistancing-frequency
//! policy consumes exactly these numbers.

#![allow(clippy::cast_possible_wrap)] // lattice indices ≤ 2^17 — i64 stencil arithmetic is exact

use crate::gridsdf::GridSdf;
use std::fmt::Write as _;

/// Redistancing evidence.
#[derive(Debug, Clone)]
pub struct RedistanceAudit {
    /// Sampled one-sided Hausdorff drift of the zero level set
    /// (before → after), in units of h.
    pub interface_drift_h: f64,
    /// Mean |(|∇φ|) − 1| within the band after redistancing.
    pub grad_dev_mean: f64,
    /// Max |(|∇φ|) − 1| within the band after redistancing.
    pub grad_dev_max: f64,
    /// FIM sweeps to convergence.
    pub sweeps: usize,
}

impl RedistanceAudit {
    /// Ledger-style JSON row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"interface_drift_h\":{:.3e},\"grad_dev_mean\":{:.3e},\
             \"grad_dev_max\":{:.3e},\"sweeps\":{}}}",
            self.interface_drift_h, self.grad_dev_mean, self.grad_dev_max, self.sweeps
        );
        s
    }
}

/// Zero crossings of the nodal field along grid lines (linear
/// interpolation) — the sampled interface.
#[must_use]
pub fn zero_crossings(phi: &GridSdf) -> Vec<[f64; 2]> {
    let n = phi.n();
    let mut pts = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let v = phi.node(i, j);
            if i < n {
                let w = phi.node(i + 1, j);
                if v == 0.0 || v * w < 0.0 {
                    let t = if v == 0.0 { 0.0 } else { v / (v - w) };
                    let p = phi.pos(i, j);
                    pts.push([p[0] + t * phi.h(), p[1]]);
                }
            }
            if j < n {
                let w = phi.node(i, j + 1);
                if v * w < 0.0 {
                    let t = v / (v - w);
                    let p = phi.pos(i, j);
                    pts.push([p[0], p[1] + t * phi.h()]);
                }
            }
        }
    }
    pts
}

/// One-sided sampled Hausdorff distance max_a min_b |a − b|.
#[must_use]
pub fn hausdorff(a: &[[f64; 2]], b: &[[f64; 2]]) -> f64 {
    let mut worst = 0.0f64;
    for p in a {
        let mut best = f64::INFINITY;
        for q in b {
            let d = (p[0] - q[0]).hypot(p[1] - q[1]);
            if d < best {
                best = d;
            }
        }
        worst = worst.max(best);
    }
    worst
}

/// Godunov eikonal update from upwind neighbor values.
fn eikonal_update(ax: f64, ay: f64, h: f64) -> f64 {
    let (a, b) = if ax <= ay { (ax, ay) } else { (ay, ax) };
    if (b - a) >= h {
        a + h
    } else {
        0.5 * (a + b + (2.0 * h * h - (b - a) * (b - a)).sqrt())
    }
}

/// Redistance in place; returns the audit.
#[allow(clippy::too_many_lines)] // freeze + FIM + audit is one narrative
pub fn redistance(phi: &mut GridSdf, band_cells: f64) -> RedistanceAudit {
    let n = phi.n();
    let h = phi.h();
    let before = zero_crossings(phi);
    let stride = n + 1;
    let big = 1e9f64;
    // Freeze interface-adjacent nodes at reconstructed distances:
    // d = |φ| / |∇φ|_local (one-term Taylor), sign preserved.
    let mut dist = vec![big; stride * stride];
    let mut sign = vec![0.0f64; stride * stride];
    let mut frozen = vec![false; stride * stride];
    for j in 0..=n {
        for i in 0..=n {
            let k = i + j * stride;
            sign[k] = if phi.node(i, j) >= 0.0 { 1.0 } else { -1.0 };
            let v = phi.node(i, j);
            let mut near = false;
            for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (ni, nj) = (i as i64 + di, j as i64 + dj);
                if ni < 0 || nj < 0 || ni > n as i64 || nj > n as i64 {
                    continue;
                }
                #[allow(clippy::cast_sign_loss)]
                let w = phi.node(ni as usize, nj as usize);
                if v == 0.0 || v * w < 0.0 {
                    near = true;
                }
            }
            if near {
                let g = phi.gradient_at(phi.pos(i, j));
                let gn = g[0].hypot(g[1]).max(1e-9);
                dist[k] = (v.abs() / gn).min(h);
                frozen[k] = true;
            }
        }
    }
    // FIM relaxation: Gauss–Seidel sweeps in four alternating orders
    // until no update exceeds the tolerance (order-independent fixed
    // point; alternating orders just reach it faster).
    let tol = 1e-12;
    let mut sweeps = 0usize;
    loop {
        sweeps += 1;
        let mut changed = false;
        let dirs = [(false, false), (true, false), (false, true), (true, true)];
        let (rev_i, rev_j) = dirs[(sweeps - 1) % 4];
        for jj in 0..=n {
            let j = if rev_j { n - jj } else { jj };
            for ii in 0..=n {
                let i = if rev_i { n - ii } else { ii };
                let k = i + j * stride;
                if frozen[k] {
                    continue;
                }
                let get = |i: i64, j: i64| -> f64 {
                    if i < 0 || j < 0 || i > n as i64 || j > n as i64 {
                        big
                    } else {
                        #[allow(clippy::cast_sign_loss)]
                        {
                            dist[i as usize + j as usize * stride]
                        }
                    }
                };
                let (i64i, i64j) = (i as i64, j as i64);
                let ax = get(i64i - 1, i64j).min(get(i64i + 1, i64j));
                let ay = get(i64i, i64j - 1).min(get(i64i, i64j + 1));
                let candidate = eikonal_update(ax.min(big), ay.min(big), h);
                if candidate < dist[k] - tol {
                    dist[k] = candidate;
                    changed = true;
                }
            }
        }
        if !changed || sweeps > 4 * (n + 1) {
            break;
        }
    }
    for j in 0..=n {
        for i in 0..=n {
            let k = i + j * stride;
            *phi.node_mut(i, j) = sign[k] * dist[k].min(big);
        }
    }
    // Audits.
    let after = zero_crossings(phi);
    let drift = hausdorff(&before, &after).max(hausdorff(&after, &before)) / h;
    let limit = band_cells * h;
    let mut dev_sum = 0.0f64;
    let mut dev_max = 0.0f64;
    let mut count = 0usize;
    for j in 1..n {
        for i in 1..n {
            if phi.node(i, j).abs() > limit {
                continue;
            }
            let gx = (phi.node(i + 1, j) - phi.node(i - 1, j)) / (2.0 * h);
            let gy = (phi.node(i, j + 1) - phi.node(i, j - 1)) / (2.0 * h);
            let dev = (gx.hypot(gy) - 1.0).abs();
            dev_sum += dev;
            dev_max = dev_max.max(dev);
            count += 1;
        }
    }
    RedistanceAudit {
        interface_drift_h: drift,
        #[allow(clippy::cast_precision_loss)]
        grad_dev_mean: dev_sum / count.max(1) as f64,
        grad_dev_max: dev_max,
        sweeps,
    }
}
