//! MANIFOLD HARMONICS (plan §7.6, bead wqd.20; [F] — behind the
//! `manifold-harmonics` feature): eigenfunctions of the Laplace–
//! Beltrami operator ON THE CURRENT SHAPE — an automatically smooth,
//! hierarchical, low-dimensional shape spectrum for CMA-ES-scale
//! global search. Low eigenvalues = low frequencies: exploring the
//! first coefficients first is FREE preconditioning.
//!
//! Signed finite-element cotan weights, lumped mass, generalized
//! eigenpairs via a dimensionless symmetric similarity
//! `M^{-1/2} L M^{-1/2}` fed to fs-la's matrix-free LOBPCG, a
//! DETERMINISTIC sign/ordering convention (P2: CMA-ES coordinates mean
//! the same thing across runs), normal-displacement parameterization,
//! and basis REFRESH with coefficient transfer so optimizers survive
//! recomputation.

use fs_la::eigen::{LobpcgState, jacobi_eigh, lobpcg_run};

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
                for (nk, fk) in n[vi as usize].iter_mut().zip(&fx) {
                    *nk += fk;
                }
            }
        }
        for v in &mut n {
            let l = norm(*v);
            if l > 0.0 {
                for vk in v.iter_mut() {
                    *vk /= l;
                }
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
                let e1 = sub(self.positions[t[1] as usize], self.positions[t[0] as usize]);
                let e2 = sub(self.positions[t[2] as usize], self.positions[t[0] as usize]);
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
    let scale = a.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let scaled = a.map(|x| x / scale);
    scale * dot(scaled, scaled).sqrt()
}

/// Sparse symmetric rows of the (positive semidefinite) cotan
/// Laplacian, plus the LUMPED mass vector (barycentric thirds).
///
/// Retains negative cotangents on obtuse triangles: the local stiffness
/// form is exactly `area * |grad(u)|²` for a piecewise-linear scalar field,
/// so its assembly is positive semidefinite without an edge-weight clamp.
/// Clamping changes that energy and the harmonic spectrum. A discrete
/// maximum principle additionally requires a suitable mesh; it is not
/// implied by positive semidefiniteness.
///
/// # Panics
/// On non-finite coordinates, invalid indices, degenerate triangles, or
/// unrepresentable triangle areas or cotangent weights.
#[must_use]
pub fn cotan_laplacian(surface: &Surface) -> (Vec<Vec<(usize, f64)>>, Vec<f64>) {
    assert!(
        surface.positions.iter().flatten().all(|x| x.is_finite()),
        "surface coordinates must be finite"
    );
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
        let twice_area = norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        assert!(
            twice_area.is_finite() && twice_area > 0.0,
            "surface triangles must have finite positive area"
        );
        let area = 0.5 * twice_area;
        for corner in 0..3 {
            mass[idx[corner]] += area / 3.0;
            // The edge OPPOSITE this corner gets cot(angle at corner)/2.
            let (i, j, k) = (idx[(corner + 1) % 3], idx[(corner + 2) % 3], corner);
            let u = sub(p[(k + 1) % 3], p[k]);
            let v = sub(p[(k + 2) % 3], p[k]);
            let cs = dot(u, v);
            let w = 0.5 * (cs / twice_area);
            assert!(w.is_finite(), "surface cotangents must be finite");
            push(i, i, w);
            push(j, j, w);
            push(i, j, -w);
            push(j, i, -w);
        }
    }
    (rows, mass)
}

/// One manifold-harmonic basis: k modes of `L ψ = λ M ψ`,
/// M-orthonormal, deterministically signed and ordered. The first modes
/// are EXPLICIT constants (λ = 0), one per connected component ordered
/// by its lowest vertex index. For a connected surface, θ₀ is uniform
/// normal offset. Independent inflation of disconnected components is
/// retained rather than collapsed into a single global constant.
#[derive(Debug, Clone)]
pub struct ManifoldBasis {
    /// The surface the basis was computed on.
    pub surface: Surface,
    /// Eigenvalues ascending, with component constants first.
    pub eigenvalues: Vec<f64>,
    /// M-orthonormal eigenfunctions, `modes[j][vertex]`.
    pub modes: Vec<Vec<f64>>,
    /// Unit vertex normals frozen at computation time.
    pub normals: Vec<[f64; 3]>,
    /// Recomputed whitened residuals `||M^-1/2 L ψ - λ M^1/2 ψ||₂`,
    /// in inverse-square geometry units (G0 numerical evidence).
    pub residuals: Vec<f64>,
}

impl ManifoldBasis {
    /// Compute `k` manifold harmonics, including component constants.
    ///
    /// The generalized problem is symmetrized: `B = M^{-1/2} L M^{-1/2}`
    /// (lumped mass), then divided by a Gershgorin bound before solving.
    /// The LOBPCG operator is `2I − B/scale − 2QQᵀ`, where Q contains
    /// the known component constants in mass-whitened coordinates.
    /// Its largest eigenvectors are the smallest positive modes of B;
    /// component constants are returned analytically. No absolute
    /// eigenvalue cutoff can erase a mode after changing geometry units.
    /// Surfaces with at most 64 vertices use the existing dense Jacobi
    /// kernel, which also permits requesting their complete spectrum.
    /// It is also the bounded fallback up to 128 vertices when LOBPCG's
    /// search-block or practical-work admission refuses the requested count.
    ///
    /// DETERMINISM CONVENTION (P2): modes ordered by ascending
    /// eigenvalue; each mode's sign fixed
    /// so its largest-|entry| coordinate (lowest index on ties) is
    /// positive. CMA-ES coordinate j means the same thing every run.
    /// Repeated eigenspaces have no unique individual eigenvectors.
    ///
    /// # Panics
    /// On invalid surface geometry, zero-mass vertices, unrepresentable
    /// mass/operator scales, `k > vertex_count`, or a request exceeding
    /// LOBPCG's block/work limits on a surface larger than 128 vertices.
    /// A search block
    /// that loses rank after nullspace projection refuses instead of
    /// publishing fewer modes than requested. Returned Ritz estimates may
    /// still be unconverged; inspect their residuals.
    #[must_use]
    pub fn compute(surface: &Surface, k: usize, iters: usize) -> ManifoldBasis {
        let (mut rows, mass) = cotan_laplacian(surface);
        let n = mass.len();
        assert!(k <= n, "harmonic count must not exceed the vertex count");
        if k == 0 {
            return ManifoldBasis {
                surface: surface.clone(),
                eigenvalues: Vec::new(),
                modes: Vec::new(),
                normals: surface.vertex_normals(),
                residuals: Vec::new(),
            };
        }
        assert!(
            mass.iter().all(|m| m.is_finite() && *m > 0.0),
            "every harmonic vertex must have finite positive mass"
        );
        let sqm: Vec<f64> = mass.iter().map(|m| m.sqrt()).collect();
        let components = connected_components(&rows);
        let component_mass: Vec<f64> = components
            .iter()
            .map(|vertices| vertices.iter().map(|&i| mass[i]).sum())
            .collect();
        assert!(
            component_mass.iter().all(|m| m.is_finite()),
            "component mass must be representable"
        );
        let mut null_weight = vec![0.0; n];
        for (vertices, m) in components.iter().zip(&component_mass) {
            let root_mass = m.sqrt();
            for &i in vertices {
                null_weight[i] = sqm[i] / root_mass;
            }
        }
        // Cache the whitened coefficients once, rather than divide in
        // every iterative matvec. Solve only after nondimensionalizing:
        // a dimensional additive shift previously erased large-scale modes.
        let mut smax = 0.0f64;
        for (i, row) in rows.iter_mut().enumerate() {
            for (j, w) in row.iter_mut() {
                *w /= sqm[i] * sqm[*j];
                assert!(w.is_finite(), "whitened harmonic weights must be finite");
            }
            smax = smax.max(row.iter().map(|(_, w)| w.abs()).sum());
        }
        assert!(
            smax.is_finite() && smax > 0.0,
            "harmonic operator scale must be finite and positive"
        );
        for row in &mut rows {
            for (_, w) in row {
                *w /= smax;
            }
        }
        let apply = |x: &[f64], y: &mut [f64]| {
            for (i, row) in rows.iter().enumerate() {
                y[i] = row.iter().map(|&(j, w)| w * x[j]).sum();
            }
        };
        let op = |x: &[f64], y: &mut [f64]| {
            apply(x, y);
            for (yi, xi) in y.iter_mut().zip(x) {
                *yi = 2.0 * xi - *yi;
            }
            for vertices in &components {
                let projection: f64 = vertices.iter().map(|&i| null_weight[i] * x[i]).sum();
                for &i in vertices {
                    y[i] -= 2.0 * null_weight[i] * projection;
                }
            }
        };
        let residual = |psi: &[f64], lambda: f64| {
            let y: Vec<f64> = psi.iter().zip(&sqm).map(|(p, m)| p * m).collect();
            let mut by = vec![0.0; n];
            apply(&y, &mut by);
            by.iter()
                .zip(&y)
                .map(|(b, v)| (b - lambda * v).powi(2))
                .sum::<f64>()
                .sqrt()
                * smax
        };
        let mut eigenvalues = Vec::with_capacity(k);
        let mut modes = Vec::with_capacity(k);
        let mut residuals = Vec::with_capacity(k);
        for (vertices, m) in components.iter().zip(&component_mass).take(k) {
            let mut psi = vec![0.0; n];
            let constant = 1.0 / m.sqrt();
            for &i in vertices {
                psi[i] = constant;
            }
            eigenvalues.push(0.0);
            residuals.push(residual(&psi, 0.0));
            modes.push(psi);
        }
        let positive_count = k - modes.len();
        let mut block_width = positive_count
            .saturating_add(3)
            .min(n / 3)
            .min(n - components.len());
        // Oversampling is optional; the numerical kernel's actual memory
        // and work admission is stronger than the geometric 3b <= n rule.
        while block_width > positive_count && !LobpcgState::shape_is_admitted(n, block_width) {
            block_width -= 1;
        }
        let iterative_admitted =
            block_width >= positive_count && LobpcgState::shape_is_admitted(n, block_width);
        let candidates: Vec<Vec<f64>> = if positive_count == 0 {
            Vec::new()
        } else if n <= 64 || (n <= 128 && !iterative_admitted) {
            let mut dense = vec![0.0; n * n];
            for (i, row) in rows.iter().enumerate() {
                for &(j, w) in row {
                    dense[i * n + j] = -w;
                }
                dense[i * n + i] += 2.0;
            }
            for vertices in &components {
                for &i in vertices {
                    for &j in vertices {
                        dense[i * n + j] -= 2.0 * (null_weight[i] * null_weight[j]);
                    }
                }
            }
            let (_, vectors) = jacobi_eigh(&dense, n);
            (n - positive_count..n)
                .rev()
                .map(|j| (0..n).map(|i| vectors[i * n + j]).collect())
                .collect()
        } else {
            assert!(
                iterative_admitted,
                "requested positive harmonics exceed LOBPCG block/work limits"
            );
            let mut state = LobpcgState::new(n, block_width);
            let ident = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
            lobpcg_run(&op, &mut state, iters, true, &ident)
                .into_iter()
                .map(|p| p.vector)
                .collect()
        };
        let mut converted = Vec::new();
        for mut y in candidates {
            // Remove only the analytically known nullspace, never modes
            // selected by an absolute or relative eigenvalue threshold.
            for _ in 0..2 {
                for vertices in &components {
                    let projection: f64 = vertices.iter().map(|&i| null_weight[i] * y[i]).sum();
                    for &i in vertices {
                        y[i] -= null_weight[i] * projection;
                    }
                }
            }
            let ynorm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
            if ynorm <= 128.0 * f64::EPSILON {
                continue;
            }
            for v in &mut y {
                *v /= ynorm;
            }
            let mut by = vec![0.0; n];
            apply(&y, &mut by);
            let lambda: f64 = y.iter().zip(&by).map(|(v, b)| v * b).sum();
            let psi: Vec<f64> = y.iter().zip(&sqm).map(|(v, m)| v / m).collect();
            let res = residual(&psi, lambda);
            converted.push((lambda * smax, psi, res));
        }
        converted.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (lam, mut psi, res) in converted {
            if eigenvalues.len() == k {
                break;
            }

            // Deterministic sign: largest-|entry| coordinate positive.
            let mut arg = 0usize;
            for (i, v) in psi.iter().enumerate() {
                if v.abs() > psi[arg].abs() {
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
        assert_eq!(
            modes.len(),
            k,
            "harmonic search block lost rank after nullspace projection; increase iterations"
        );
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
            for (pk, nk) in p.iter_mut().zip(&self.normals[i]) {
                *pk += amp * nk;
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

// Triangle incidence is retained in the rows even for a zero cotangent.
// Thus connectivity does not change when a right angle makes one weight zero.
fn connected_components(rows: &[Vec<(usize, f64)>]) -> Vec<Vec<usize>> {
    let mut visited = vec![false; rows.len()];
    let mut components = Vec::new();
    for seed in 0..rows.len() {
        if visited[seed] {
            continue;
        }
        let mut vertices = Vec::new();
        let mut stack = vec![seed];
        visited[seed] = true;
        while let Some(i) = stack.pop() {
            vertices.push(i);
            for &(j, _) in &rows[i] {
                if !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
        vertices.sort_unstable();
        components.push(vertices);
    }
    components
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
