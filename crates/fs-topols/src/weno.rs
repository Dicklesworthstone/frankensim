//! WENO5 narrow-band level-set advection with TVD-RK3 time stepping.
//!
//! Two Hamiltonians, one stencil machinery:
//! - LINEAR advection `φ_t + u·∇φ = 0` (upwinded per component) — the
//!   G1 order battery runs a rotating field one revolution and
//!   measures the design order on smooth data.
//! - NORMAL flow `φ_t + v_n|∇φ| = 0` (Godunov/Rouy–Tourin flux) — the
//!   optimizer's evolution law, with v_n extended off the interface.
//!
//! The band: nodes with |φ| ≤ band·h advance; the frozen far field is
//! refreshed by redistancing (fim). CFL is tied to the band: the
//! interface must not cross more than `cfl` cells per step.

#![allow(clippy::cast_possible_wrap)] // lattice indices ≤ 2^17 — i64 stencil arithmetic is exact

use crate::gridsdf::GridSdf;

/// One-sided WENO5 derivative from the five upwind differences
/// (Jiang–Shu weights, ε = 1e-6).
fn weno5(v1: f64, v2: f64, v3: f64, v4: f64, v5: f64) -> f64 {
    let eps = 1e-6;
    let p1 = v1 / 3.0 - 7.0 * v2 / 6.0 + 11.0 * v3 / 6.0;
    let p2 = -v2 / 6.0 + 5.0 * v3 / 6.0 + v4 / 3.0;
    let p3 = v3 / 3.0 + 5.0 * v4 / 6.0 - v5 / 6.0;
    let s1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).powi(2) + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3).powi(2);
    let s2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).powi(2) + 0.25 * (v2 - v4).powi(2);
    let s3 = 13.0 / 12.0 * (v3 - 2.0 * v4 + v5).powi(2) + 0.25 * (3.0 * v3 - 4.0 * v4 + v5).powi(2);
    let a1 = 0.1 / (eps + s1).powi(2);
    let a2 = 0.6 / (eps + s2).powi(2);
    let a3 = 0.3 / (eps + s3).powi(2);
    (a1 * p1 + a2 * p2 + a3 * p3) / (a1 + a2 + a3)
}

/// The clamped-index nodal difference D⁻φ along an axis.
fn diff(phi: &GridSdf, i: i64, j: i64, axis: usize) -> f64 {
    let n = phi.n() as i64;
    let cl = |k: i64| k.clamp(0, n) as usize;
    let (a, b) = if axis == 0 {
        ((cl(i), cl(j)), (cl(i - 1), cl(j)))
    } else {
        ((cl(i), cl(j)), (cl(i), cl(j - 1)))
    };
    (phi.node(a.0, a.1) - phi.node(b.0, b.1)) / phi.h()
}

/// WENO5 left/right derivative pair at a node along an axis.
fn weno_pair(phi: &GridSdf, i: usize, j: usize, axis: usize) -> (f64, f64) {
    let (i, j) = (i as i64, j as i64);
    let d = |k: i64| {
        if axis == 0 {
            diff(phi, i + k, j, 0)
        } else {
            diff(phi, i, j + k, 1)
        }
    };
    let minus = weno5(d(-2), d(-1), d(0), d(1), d(2));
    let plus = weno5(d(3), d(2), d(1), d(0), d(-1));
    (minus, plus)
}

/// One Hamiltonian evaluation at every band node: returns −H (the
/// time derivative), zero off the band.
fn rhs(phi: &GridSdf, band: &[bool], velocity: &Velocity<'_>) -> Vec<f64> {
    let n = phi.n();
    let mut out = vec![0.0f64; (n + 1) * (n + 1)];
    for j in 0..=n {
        for i in 0..=n {
            let k = i + j * (n + 1);
            if !band[k] {
                continue;
            }
            let p = phi.pos(i, j);
            let (xm, xp) = weno_pair(phi, i, j, 0);
            let (ym, yp) = weno_pair(phi, i, j, 1);
            let h = match velocity {
                Velocity::Linear(u) => {
                    let uv = u(p[0], p[1]);
                    uv[0] * if uv[0] >= 0.0 { xm } else { xp }
                        + uv[1] * if uv[1] >= 0.0 { ym } else { yp }
                }
                Velocity::Normal(vn) => {
                    let v = vn[k];
                    // Godunov (Rouy–Tourin) |∇φ| flux.
                    let grad = if v >= 0.0 {
                        (xm.max(0.0).powi(2).max(xp.min(0.0).powi(2))
                            + ym.max(0.0).powi(2).max(yp.min(0.0).powi(2)))
                        .sqrt()
                    } else {
                        (xp.max(0.0).powi(2).max(xm.min(0.0).powi(2))
                            + yp.max(0.0).powi(2).max(ym.min(0.0).powi(2)))
                        .sqrt()
                    };
                    v * grad
                }
            };
            out[k] = -h;
        }
    }
    out
}

/// The advection velocity: an external field or a nodal normal speed.
pub enum Velocity<'a> {
    /// External linear advection u(x, y).
    Linear(&'a dyn Fn(f64, f64) -> [f64; 2]),
    /// Nodal normal speed (one value per node).
    Normal(&'a [f64]),
}

/// Nodes within `band_cells` cells of the interface.
#[must_use]
pub fn build_band(phi: &GridSdf, band_cells: f64) -> Vec<bool> {
    let limit = band_cells * phi.h();
    phi.nodes().iter().map(|&v| v.abs() <= limit).collect()
}

/// Advance the level set by `total_time` with TVD-RK3 + WENO5 under a
/// CFL number; returns the number of steps taken.
pub fn advect(
    phi: &mut GridSdf,
    band: &[bool],
    velocity: &Velocity<'_>,
    total_time: f64,
    cfl: f64,
) -> usize {
    let n = phi.n();
    let vmax = match velocity {
        Velocity::Linear(u) => {
            let mut m = 0.0f64;
            for j in 0..=n {
                for i in 0..=n {
                    let p = phi.pos(i, j);
                    let uv = u(p[0], p[1]);
                    m = m.max(uv[0].abs().max(uv[1].abs()));
                }
            }
            m
        }
        Velocity::Normal(vn) => vn.iter().fold(0.0f64, |m, v| m.max(v.abs())),
    };
    if vmax <= 0.0 || total_time <= 0.0 {
        return 0;
    }
    let dt_max = cfl * phi.h() / vmax;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let steps = (total_time / dt_max).ceil() as usize;
    #[allow(clippy::cast_precision_loss)]
    let dt = total_time / steps as f64;
    for _ in 0..steps {
        let phi0 = phi.clone();
        // Stage 1.
        let k1 = rhs(phi, band, velocity);
        axpy(phi.nodes_mut(), phi0.nodes(), &k1, 1.0, dt);
        // Stage 2: φ = 3/4 φ0 + 1/4 (φ1 + dt L(φ1)).
        let k2 = rhs(phi, band, velocity);
        let phi1 = phi.nodes().to_vec();
        for (idx, v) in phi.nodes_mut().iter_mut().enumerate() {
            *v = 0.75 * phi0.nodes()[idx] + 0.25 * (phi1[idx] + dt * k2[idx]);
        }
        // Stage 3: φ = 1/3 φ0 + 2/3 (φ2 + dt L(φ2)).
        let k3 = rhs(phi, band, velocity);
        let phi2 = phi.nodes().to_vec();
        for (idx, v) in phi.nodes_mut().iter_mut().enumerate() {
            *v = phi0.nodes()[idx] / 3.0 + 2.0 / 3.0 * (phi2[idx] + dt * k3[idx]);
        }
    }
    steps
}

fn axpy(out: &mut [f64], base: &[f64], k: &[f64], a: f64, dt: f64) {
    for ((o, b), kv) in out.iter_mut().zip(base).zip(k) {
        *o = a * b + dt * kv;
    }
}
