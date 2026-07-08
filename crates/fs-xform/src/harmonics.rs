//! MANIFOLD HARMONICS (plan §7.6, bead wqd.20; [F] — behind the
//! `manifold-harmonics` feature): eigenfunctions of the Laplace–
//! Beltrami operator ON THE CURRENT SHAPE — an automatically smooth,
//! hierarchical, low-dimensional shape spectrum for CMA-ES-scale
//! global search. Low eigenvalues = low frequencies: exploring the
//! first coefficients first is FREE preconditioning.
//!
//! Cotan-weight assembly with a robust clamp for bad triangles, lumped
//! mass, generalized eigenpairs via the symmetric similarity
//! `M^{-1/2} L M^{-1/2}` fed to fs-la's matrix-free LOBPCG, a
//! DETERMINISTIC sign/ordering convention (P2: CMA-ES coordinates mean
//! the same thing across runs), normal-displacement parameterization,
//! and basis REFRESH with coefficient transfer so optimizers survive
//! recomputation.

use fs_la::eigen::{EigenPair, LobpcgState, lobpcg_run};

/// A triangle surface (positions + triangle indices) — the input is
/// deliberately plain so no body-fitted volumetric mesh is required.
#[derive(Debug, Clone, PartialEq)]
pub struct Surface {
    /// Vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Triangles (CCW outward).
    pub triangles: Vec<[u32; 3]>,
}

impl Surface {
    /// Area-weighted vertex normals (unit).
    #[must_use]
    pub fn vertex_normals(&self) -> Vec<[f64; 3]> {
        let mut n = vec![[0.0f64; 3]; self.positions.len()];
        for t in &self.triangles {
            let [a, b, c] = [
                self.positions[t[0] as usize],
                self.positions[t[1] as usize],
                self.positions[t[2] as usize],
            ];
            let e1 = sub(b, a);
            let e2 = sub(c, a);
            let fx = cross(e1, e2); // area-weighted
            for &vi in t {
                for k in 0..3 {
                    n[vi as usize][k] += fx[k];
                }
            }
        }
        for v in &mut n {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-300);
            for k in 0..3 {
                v[k] /= l;
            }
        }
        n
    }

    /// Total surface area.
    #[must_use]
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let e1 = sub(
                    self.positions[t[1] as usize],
                    self.positions[t[0] as usize],
                );
                let e2 = sub(
                    self.positions[t[2] as usize],
                    self.positions[t[0] as usize],
                );
                0.5 * norm(cross(e1, e2))
            })
            .sum()
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Sparse symmetric rows of the (positive semidefinite) cotan
/// Laplacian, plus the LUMPED mass vector (barycentric thirds).
///
/// Robust variant: cotangents are clamped to `>= 0` (obtuse triangles
/// would inject negative off-diagonal weights that break the maximum
/// principle on bad meshes — the clamp trades a little consistency for
/// unconditional stability; documented no-claim).
#[must_use]
pub fn cotan_laplacian(surface: &Surface) -> (Vec<Vec<(usize, f64)>>, Vec<f64>) {
    let nv = surface.positions.len();
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nv];
    let mut mass = vec![0.0f64; nv];
    let mut push = |i: usize, j: usize, w: f64| {
        if let Some(e) = rows[i].iter_mut().find(|(c, _)| *c == j) {
            e.1 += w;
        } else {
            rows[i].push((j, w));
        }
    };
    for t in &surface.triangles {
        let idx = [t[0] as usize, t[1] as usize, t[2] as usize];
        let p = [
            surface.positions[idx[0]],
            surface.positions[idx[1]],
            surface.positions[idx[2]],
        ];
        let area = 0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        for corner in 0..3 {
            mass[idx[corner]] += area / 3.0;
            // The edge OPPOSITE this corner gets cot(angle at corner)/2.
            let (i, j, k) = (idx[(corner + 1) % 3], idx[(corner + 2) % 3], corner);
            let u = sub(p[(k + 1) % 3], p[k]);
            let v = sub(p[(k + 2) % 3], p[k]);
            let cs = dot(u, v);
            let sn = norm(cross(u, v)).max(1e-300);
            let w = (0.5 * cs / sn).max(0.0); // robust clamp
            push(i, i, w);
            push(j, j, w);
            push(i, j, -w);
            push(j, i, -w);
        }
    }
    (rows, mass)
}

/// One manifold-harmonic basis: k modes of `L ψ = λ M ψ`,
/// M-orthonormal, deterministically signed and ordered. Mode 0 is the
/// EXPLICIT constant (λ = 0): θ₀ is uniform normal offset — inflation
/// is a legitimate design direction, so the eigensolver's kernel
/// deflation re-admits it analytically rather than numerically.
#[derive(Debug, Clone)]
pub struct ManifoldBasis {
    /// The surface the basis was computed on.
    pub surface: Surface,
    /// Eigenvalues ascending (first NON-trivial mode first).
    pub eigenvalues: Vec<f64>,
    /// M-orthonormal eigenfunctions, `modes[j][vertex]`.
    pub modes: Vec<Vec<f64>>,
    /// Unit vertex normals frozen at computation time.
    pub normals: Vec<[f64; 3]>,
    /// LOBPCG true residuals per mode (G0 evidence).
    pub residuals: Vec<f64>,
}

impl ManifoldBasis {
    /// Compute the first `k` non-trivial manifold harmonics.
    ///
    /// The generalized problem is symmetrized: `B = M^{-1/2} L M^{-1/2}`
    /// (lumped mass), smallest-k via LOBPCG on `sI − B` (spectral
    /// shift-and-negate keeps the solver's LARGEST-mode path exact),
    /// then `ψ = M^{-1/2} y`, M-normalized.
    ///
    /// DETERMINISM CONVENTION (P2): modes ordered by ascending
    /// eigenvalue with vertex-0-entry tie-break; each mode's sign fixed
    /// so its largest-|entry| coordinate (lowest index on ties) is
    /// positive. CMA-ES coordinate j means the same thing every run.
    #[must_use]
    pub fn compute(surface: &Surface, k: usize, iters: usize) -> ManifoldBasis {
        let (rows, mass) = cotan_laplacian(surface);
        let n = mass.len();
        let sqm: Vec<f64> = mass.iter().map(|m| m.max(1e-300).sqrt()).collect();
        // Gershgorin upper bound for B = M^{-1/2} L M^{-1/2}.
        let mut smax = 0.0f64;
        for (i, row) in rows.iter().enumerate() {
            let r: f64 = row
                .iter()
                .map(|&(j, w)| (w / (sqm[i] * sqm[j])).abs())
                .sum();
            smax = smax.max(r);
        }
        let shift = smax * 1.01 + 1.0;
        // Operator: y <- (shift·I − B)·x (largest of this = smallest of B).
        let op = |x: &[f64], y: &mut [f64]| {
            for (i, row) in rows.iter().enumerate() {
                let mut acc = 0.0f64;
                for &(j, w) in row {
                    acc += w / (sqm[i] * sqm[j]) * x[j];
                }
                y[i] = shift * x[i] - acc;
            }
        };
        let ident = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
        // k requested + 1 trivial + slack for convergence.
        let b = (k + 3).min(n / 3);
        let mut state = LobpcgState::new(n, b);
        let pairs: Vec<EigenPair> = lobpcg_run(&op, &mut state, iters, true, &ident);
        // Convert back: eigenvalue of B = shift − ritz; ψ = M^{-1/2} y.
        let mut converted: Vec<(f64, Vec<f64>, f64)> = pairs
            .into_iter()
            .map(|p| {
                let lam = shift - p.value;
                let mut psi: Vec<f64> = p.vector.iter().zip(&sqm).map(|(v, s)| v / s).collect();
                // M-normalize.
                let mnorm: f64 = psi
                    .iter()
                    .zip(&mass)
                    .map(|(v, m)| v * v * m)
                    .sum::<f64>()
                    .sqrt()
                    .max(1e-300);
                for v in &mut psi {
                    *v /= mnorm;
                }
                (lam, psi, p.residual)
            })
            .collect();
        converted.sort_by(|a, b| a.0.total_cmp(&b.0));
        // Deflate the numeric kernel, then re-admit the constant
        // ANALYTICALLY as mode 0 (θ₀ = uniform inflation).
        let total_mass: f64 = mass.iter().sum();
        let mut eigenvalues = vec![0.0f64];
        let mut modes = vec![vec![1.0 / total_mass.sqrt(); n]];
        let mut residuals = vec![0.0f64];
        for (lam, mut psi, res) in converted {
            if lam < 1e-8 * shift {
                continue; // constant / rigid kernel
            }
            if eigenvalues.len() == k {
                break;
            }

            // Deterministic sign: largest-|entry| coordinate positive.
            let mut arg = 0usize;
            for (i, v) in psi.iter().enumerate() {
                if v.abs() > psi[arg].abs() + 1e-14 {
                    arg = i;
                }
            }
            if psi[arg] < 0.0 {
                for v in &mut psi {
                    *v = -*v;
                }
            }
            eigenvalues.push(lam);
            modes.push(psi);
            residuals.push(res);
        }
        ManifoldBasis {
            surface: surface.clone(),
            eigenvalues,
            modes,
            normals: surface.vertex_normals(),
            residuals,
        }
    }

    /// Number of coefficients.
    #[must_use]
    pub fn dof(&self) -> usize {
        self.modes.len()
    }

    /// The displaced surface for spectral coefficients θ:
    /// `x_i + Σ_j θ_j ψ_j(i) n̂_i` (normal displacement).
    #[must_use]
    pub fn displace(&self, theta: &[f64]) -> Surface {
        let mut out = self.surface.clone();
        for (i, p) in out.positions.iter_mut().enumerate() {
            let mut amp = 0.0f64;
            for (t, mode) in theta.iter().zip(&self.modes) {
                amp += t * mode[i];
            }
            for k in 0..3 {
                p[k] += amp * self.normals[i][k];
            }
        }
        out
    }

    /// M-weighted projection of an arbitrary normal-displacement field
    /// onto this basis (the coefficient-transfer primitive).
    #[must_use]
    pub fn project(&self, normal_field: &[f64]) -> Vec<f64> {
        let (_, mass) = cotan_laplacian(&self.surface);
        self.modes
            .iter()
            .map(|mode| {
                mode.iter()
                    .zip(normal_field)
                    .zip(&mass)
                    .map(|((psi, f), m)| psi * f * m)
                    .sum()
            })
            .collect()
    }

    /// Dirichlet energy of mode j (== λ_j for M-orthonormal modes —
    /// the smoothness-ordering witness).
    #[must_use]
    pub fn dirichlet_energy(&self, j: usize) -> f64 {
        let (rows, _) = cotan_laplacian(&self.surface);
        let psi = &self.modes[j];
        rows.iter()
            .enumerate()
            .map(|(i, row)| {
                let li: f64 = row.iter().map(|&(c, w)| w * psi[c]).sum();
                psi[i] * li
            })
            .sum()
    }
}

/// Basis-refresh drift criterion: refresh when the mean vertex
/// displacement exceeds `frac` of the bounding-box diagonal.
#[must_use]
pub fn needs_refresh(old: &Surface, new: &Surface, frac: f64) -> bool {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &old.positions {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let diag = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    let mean: f64 = old
        .positions
        .iter()
        .zip(&new.positions)
        .map(|(a, b)| norm(sub(*b, *a)))
        .sum::<f64>()
        / old.positions.len().max(1) as f64;
    mean > frac * diag
}

/// COEFFICIENT TRANSFER across a basis refresh: express the OLD basis's
/// displacement in the NEW basis (M-weighted projection of the old
/// normal field onto the new modes), so optimizers survive refreshes.
/// Returns the new coefficients.
#[must_use]
pub fn transfer(old: &ManifoldBasis, new: &ManifoldBasis, theta_old: &[f64]) -> Vec<f64> {
    // The old displacement's normal amplitude per vertex (old and new
    // bases share the vertex set across a refresh-in-place).
    let n = old.surface.positions.len();
    let mut field = vec![0.0f64; n];
    for (t, mode) in theta_old.iter().zip(&old.modes) {
        for (fi, m) in field.iter_mut().zip(mode) {
            *fi += t * m;
        }
    }
    // Re-express along the NEW normals: scale by n̂_old · n̂_new.
    for ((fi, no), nn) in field.iter_mut().zip(&old.normals).zip(&new.normals) {
        *fi *= dot(*no, *nn);
    }
    new.project(&field)
}
