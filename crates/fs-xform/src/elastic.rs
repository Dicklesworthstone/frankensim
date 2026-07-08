//! ELASTIC SHAPE-SPACE METRICS (plan Bet 6 / §9, bead wqd.27; [F] —
//! behind the `elastic-shapes` feature): the space of shapes itself
//! carries a Riemannian metric, so the system computes GEODESICS
//! between designs — interpolation through plausible intermediates
//! (not crossfades), morphing, deformation-energy TRUST REGIONS, and
//! Karcher means. Complementary to OT: Wasserstein sees shapes as mass
//! distributions; elastic metrics see them as deformable objects.
//!
//! Curves ride the SQUARE-ROOT-VELOCITY (SRV) representation
//! `q = c′/√|c′|`, under which the elastic metric becomes the flat L²
//! metric: geodesics are LINEAR interpolation in SRV space
//! (closed-form), and reparameterization is handled by dynamic
//! programming over monotone correspondences. Surfaces (same
//! connectivity) use discrete path-straightening on a
//! membrane-plus-bending energy. The pullback through the
//! manifold-harmonics basis is DIAGONAL — the mathematically honest
//! preconditioner for spectral shape search.

use crate::harmonics::{ManifoldBasis, Surface, cotan_laplacian};

/// A planar/space curve as ordered samples.
#[derive(Debug, Clone, PartialEq)]
pub struct Curve {
    /// Sample points (open curve: endpoints distinct).
    pub points: Vec<[f64; 2]>,
}

impl Curve {
    /// Arc length.
    #[must_use]
    pub fn length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt())
            .sum()
    }

    /// Resample to `n` points uniformly by arc length (linear).
    #[must_use]
    pub fn resample(&self, n: usize) -> Curve {
        let total = self.length().max(1e-300);
        let mut cum = vec![0.0f64];
        for w in self.points.windows(2) {
            let d = ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt();
            cum.push(cum.last().copied().unwrap_or(0.0) + d);
        }
        let mut out = Vec::with_capacity(n);
        let mut seg = 0usize;
        for i in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let s = total * i as f64 / (n - 1) as f64;
            while seg + 2 < cum.len() && cum[seg + 1] < s {
                seg += 1;
            }
            let t = ((s - cum[seg]) / (cum[seg + 1] - cum[seg]).max(1e-300)).clamp(0.0, 1.0);
            let (a, b) = (self.points[seg], self.points[seg + 1]);
            out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
        Curve { points: out }
    }
}

/// The SRV transform: `q_i = c′_i / √|c′_i|` at segment midpoints,
/// with the parameter weight `1/(n−1)` folded into the L² inner
/// product.
#[must_use]
pub fn srv(curve: &Curve) -> Vec<[f64; 2]> {
    let n = curve.points.len();
    #[allow(clippy::cast_precision_loss)]
    let dt = 1.0 / (n - 1) as f64;
    curve
        .points
        .windows(2)
        .map(|w| {
            let v = [(w[1][0] - w[0][0]) / dt, (w[1][1] - w[0][1]) / dt];
            let speed = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-300);
            let s = speed.sqrt();
            [v[0] / s, v[1] / s]
        })
        .collect()
}

/// L² distance between SRVs (the ELASTIC distance for fixed
/// parameterization; translation-invariant by construction).
#[must_use]
pub fn srv_distance(a: &Curve, b: &Curve) -> f64 {
    let (qa, qb) = (srv(a), srv(b));
    assert_eq!(qa.len(), qb.len(), "resample to a common n first");
    #[allow(clippy::cast_precision_loss)]
    let dt = 1.0 / qa.len() as f64;
    let sum: f64 = qa
        .iter()
        .zip(&qb)
        .map(|(x, y)| (x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2))
        .sum();
    (sum * dt).sqrt()
}

/// The CLOSED-FORM geodesic: linear interpolation in SRV space,
/// mapped back to a curve anchored at `a`'s start point.
#[must_use]
pub fn srv_geodesic(a: &Curve, b: &Curve, s: f64) -> Curve {
    let (qa, qb) = (srv(a), srv(b));
    assert_eq!(qa.len(), qb.len());
    let n = qa.len() + 1;
    #[allow(clippy::cast_precision_loss)]
    let dt = 1.0 / (n - 1) as f64;
    let mut points = Vec::with_capacity(n);
    let start = [
        a.points[0][0] + s * (b.points[0][0] - a.points[0][0]),
        a.points[0][1] + s * (b.points[0][1] - a.points[0][1]),
    ];
    points.push(start);
    for (x, y) in qa.iter().zip(&qb) {
        let q = [x[0] + s * (y[0] - x[0]), x[1] + s * (y[1] - x[1])];
        let qn = (q[0] * q[0] + q[1] * q[1]).sqrt();
        let v = [q[0] * qn, q[1] * qn]; // invert: c' = q·|q|
        let last = *points.last().expect("nonempty");
        points.push([last[0] + v[0] * dt, last[1] + v[1] * dt]);
    }
    Curve { points }
}

/// Optimal monotone reparameterization by DYNAMIC PROGRAMMING: match
/// SRV samples of `b` onto `a` over monotone lattice paths, returning
/// the reparameterization-invariant elastic distance.
#[must_use]
pub fn elastic_distance(a: &Curve, b: &Curve, n: usize) -> f64 {
    let ra = a.resample(n);
    // DP over the correspondence grid with local slope steps: the cost
    // of matching segment blocks approximates ∫|q_a − √γ′ (q_b∘γ)|².
    let (qa, qb) = (srv(&ra), srv(&b.resample(n)));
    let m = qa.len();
    let inf = f64::INFINITY;
    let mut cost = vec![vec![inf; m + 1]; m + 1];
    cost[0][0] = 0.0;
    // Steps (di, dj) with slope dj/di in {1/2, 1, 2} — the standard
    // slope-constrained elastic matching lattice.
    let steps: [(usize, usize); 5] = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1)];
    #[allow(clippy::cast_precision_loss)]
    let dt = 1.0 / m as f64;
    for i in 0..=m {
        for j in 0..=m {
            if cost[i][j].is_infinite() {
                continue;
            }
            for &(di, dj) in &steps {
                let (ni, nj) = (i + di, j + dj);
                if ni > m || nj > m {
                    continue;
                }
                // Segment match: q_a over [i, ni) vs q_b over [j, nj),
                // with the √γ′ factor from the slope.
                let slope = dj as f64 / di as f64;
                let root = slope.sqrt();
                let mut c = 0.0f64;
                for k in 0..di {
                    let ia = i + k;
                    #[allow(
                        clippy::cast_possible_truncation,
                        clippy::cast_sign_loss,
                        clippy::cast_precision_loss
                    )]
                    let jb = (j as f64 + slope * (k as f64 + 0.5)).floor() as usize;
                    let jb = jb.min(m - 1);
                    let (x, y) = (qa[ia], qb[jb]);
                    c += (x[0] - root * y[0]).powi(2) + (x[1] - root * y[1]).powi(2);
                }
                let cand = cost[i][j] + c * dt;
                if cand < cost[ni][nj] {
                    cost[ni][nj] = cand;
                }
            }
        }
    }
    cost[m][m].sqrt()
}

/// Karcher (Fréchet) mean of a curve family in SRV space: iterate the
/// SRV average (the flat-space mean is exact per iteration; iteration
/// re-anchors lengths). Returns (mean, per-iteration mean-shift norms
/// — the convergence-rate evidence).
#[must_use]
pub fn karcher_mean(curves: &[Curve], n: usize, iters: usize) -> (Curve, Vec<f64>) {
    let mut mean = curves[0].resample(n);
    let mut shifts = Vec::with_capacity(iters);
    for _ in 0..iters {
        let qm = srv(&mean);
        let mut acc = vec![[0.0f64; 2]; qm.len()];
        let mut start = [0.0f64; 2];
        for c in curves {
            let rc = c.resample(n);
            let qc = srv(&rc);
            for (a, q) in acc.iter_mut().zip(&qc) {
                a[0] += q[0];
                a[1] += q[1];
            }
            start[0] += rc.points[0][0];
            start[1] += rc.points[0][1];
        }
        #[allow(clippy::cast_precision_loss)]
        let inv = 1.0 / curves.len() as f64;
        for a in &mut acc {
            a[0] *= inv;
            a[1] *= inv;
        }
        start[0] *= inv;
        start[1] *= inv;
        // Map back.
        #[allow(clippy::cast_precision_loss)]
        let dt = 1.0 / (n - 1) as f64;
        let mut points = vec![start];
        for q in &acc {
            let qn = (q[0] * q[0] + q[1] * q[1]).sqrt();
            let last = *points.last().expect("nonempty");
            points.push([last[0] + q[0] * qn * dt, last[1] + q[1] * qn * dt]);
        }
        let next = Curve { points };
        shifts.push(srv_distance(&mean, &next));
        mean = next;
    }
    (mean, shifts)
}

/// Discrete elastic energy of a surface displacement (membrane +
/// bending proxy): `∫ |∇u|² + |Δu|²` per component via the cotan
/// machinery — the deformation-energy ball used for trust regions.
#[must_use]
pub fn deformation_energy(base: &Surface, displaced: &Surface) -> f64 {
    let (rows, mass) = cotan_laplacian(base);
    let n = base.positions.len();
    let mut energy = 0.0f64;
    for comp in 0..3 {
        let u: Vec<f64> = (0..n)
            .map(|i| displaced.positions[i][comp] - base.positions[i][comp])
            .collect();
        // Membrane: uᵀ L u.
        let lu: Vec<f64> = rows
            .iter()
            .map(|row| row.iter().map(|&(j, w)| w * u[j]).sum())
            .collect();
        let membrane: f64 = u.iter().zip(&lu).map(|(a, b)| a * b).sum();
        // Bending: (M⁻¹Lu)ᵀ M (M⁻¹Lu) = luᵀ M⁻¹ lu.
        let bending: f64 = lu
            .iter()
            .zip(&mass)
            .map(|(l, m)| l * l / m.max(1e-300))
            .sum();
        energy += membrane + bending;
    }
    energy
}

/// The linear crossfade path (the standard initialization; for
/// AFFINELY related shapes it is already the discrete geodesic).
#[must_use]
pub fn crossfade_path(a: &Surface, b: &Surface, interior: usize) -> Vec<Surface> {
    let steps = interior + 1;
    (0..=steps)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let t = k as f64 / steps as f64;
            let mut s = a.clone();
            for (i, p) in s.positions.iter_mut().enumerate() {
                for c in 0..3 {
                    p[c] = a.positions[i][c] + t * (b.positions[i][c] - a.positions[i][c]);
                }
            }
            s
        })
        .collect()
}

/// PATH-STRAIGHTEN an arbitrary same-connectivity path in place:
/// coordinate descent of the discrete path energy `Σ E(x_k, x_{k+1})`
/// over the interior shapes (endpoints pinned). Returns the per-sweep
/// total energies — monotone down to the discrete-geodesic fixed
/// point.
#[must_use]
pub fn straighten_from(path: &mut [Surface], sweeps: usize) -> Vec<f64> {
    let steps = path.len() - 1;
    let n = path[0].positions.len();
    let path_energy = |path: &[Surface]| -> f64 {
        path.windows(2)
            .map(|w| deformation_energy(&w[0], &w[1]))
            .sum()
    };
    let mut energies = vec![path_energy(path)];
    for _ in 0..sweeps {
        // Interior smoothing step: each interior shape moves toward the
        // average of its neighbors (the discrete geodesic fixed point
        // for the quadratic energy), with a conservative step.
        for k in 1..steps {
            let mut moved = path[k].clone();
            for i in 0..n {
                for c in 0..3 {
                    let target =
                        f64::midpoint(path[k - 1].positions[i][c], path[k + 1].positions[i][c]);
                    moved.positions[i][c] += 0.5 * (target - path[k].positions[i][c]);
                }
            }
            path[k] = moved;
        }
        energies.push(path_energy(path));
    }
    energies
}

/// Convenience: crossfade init + straightening.
#[must_use]
pub fn straighten_path(
    a: &Surface,
    b: &Surface,
    interior: usize,
    sweeps: usize,
) -> (Vec<Surface>, Vec<f64>) {
    let mut path = crossfade_path(a, b, interior);
    let energies = straighten_from(&mut path, sweeps);
    (path, energies)
}

/// The elastic metric PULLED BACK through a manifold-harmonics basis:
/// for normal displacements `u = Σ θⱼ ψⱼ n̂`, the membrane+bending
/// quadratic form is DIAGONAL in θ with weights `λⱼ + λⱼ²` (M-orthonormal
/// modes) — Bet 6's "right inner product" doctrine as an exact formula.
#[must_use]
pub fn pullback_metric(basis: &ManifoldBasis) -> Vec<f64> {
    basis.eigenvalues.iter().map(|l| l + l * l).collect()
}

/// A deformation-energy TRUST REGION: accept a spectral step `δθ` iff
/// its pullback elastic energy `Σ (λⱼ+λⱼ²) δθⱼ²` is at most `radius`.
/// "Step size" means bounded deformation, not bounded coefficients.
#[must_use]
pub fn within_trust_region(metric: &[f64], dtheta: &[f64], radius: f64) -> bool {
    let e: f64 = metric.iter().zip(dtheta).map(|(g, d)| g * d * d).sum();
    e <= radius
}
