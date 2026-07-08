//! 3D exterior potential flow: constant SOURCE panels with collocation
//! Neumann conditions. The influence matrix row is
//! `A_ij = n_i · ∇_x G(x_i, y_j) · area_j` with the flat-panel
//! centroid approximation for well-separated pairs and the exact
//! self-term σ/2 on the diagonal — a documented screening-grade
//! discretization whose measured convergence on the sphere IS the
//! gate. The GMRES matvec runs three fs-fmm passes (the gradient's
//! components are each a smooth kernel) dotted with target normals;
//! the dense assembly is the oracle it must match.

use fs_fmm::{Fmm, Kernel};
use fs_geom::Point3;
use fs_rep_mesh::shapes::icosphere;
use fs_solver::krylov::GmresState;
use fs_solver::op::LinearOp;

/// One TRUE gradient component of the Laplace kernel:
/// `∂G/∂x_c = −(x_c − y_c)/(4π|x−y|³)` — sign conventions verified by
/// the uniform-sphere Gauss identity in the battery (row action on
/// ones ≈ −1).
struct GradKernel {
    c: usize,
}

impl Kernel for GradKernel {
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> f64 {
        let d = [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 < 1e-300 {
            return 0.0;
        }
        let r = r2.sqrt();
        -d[self.c] / (4.0 * std::f64::consts::PI * r2 * r)
    }
}

/// A panelized closed surface: centroids, outward normals, areas.
pub struct SpherePanels {
    /// Panel centroids.
    pub centroids: Vec<[f64; 3]>,
    /// Outward unit normals.
    pub normals: Vec<[f64; 3]>,
    /// Panel areas.
    pub areas: Vec<f64>,
}

impl SpherePanels {
    /// Panelize an icosphere (fs-rep-mesh) of given radius/subdivisions.
    #[must_use]
    pub fn icosphere(radius: f64, subdivisions: u32) -> SpherePanels {
        let soup = icosphere(Point3::new(0.0, 0.0, 0.0), radius, subdivisions);
        let mut centroids = Vec::new();
        let mut normals = Vec::new();
        let mut areas = Vec::new();
        for t in 0..soup.triangles.len() {
            let [a, b, c] = soup.tri(t);
            let cx = [
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            ];
            let u = [b.x - a.x, b.y - a.y, b.z - a.z];
            let v = [c.x - a.x, c.y - a.y, c.z - a.z];
            let mut n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let area = 0.5 * len;
            for x in &mut n {
                *x /= len;
            }
            // Outward (sphere centered at origin): flip if pointing in.
            if n[0] * cx[0] + n[1] * cx[1] + n[2] * cx[2] < 0.0 {
                for x in &mut n {
                    *x = -*x;
                }
            }
            centroids.push(cx);
            normals.push(n);
            areas.push(area);
        }
        SpherePanels {
            centroids,
            normals,
            areas,
        }
    }

    /// Dense influence matrix (row-major n×n): normal velocity at
    /// centroid i induced by unit source density on panel j.
    #[must_use]
    pub fn dense_matrix(&self) -> Vec<f64> {
        let n = self.centroids.len();
        let gk = [
            GradKernel { c: 0 },
            GradKernel { c: 1 },
            GradKernel { c: 2 },
        ];
        let mut a = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    a[i * n + j] = -0.5; // outside-limit jump −σ/2
                    continue;
                }
                let mut v = 0.0;
                for (c, k) in gk.iter().enumerate() {
                    v += self.normals[i][c] * k.eval(self.centroids[i], self.centroids[j]);
                }
                a[i * n + j] = v * self.areas[j];
            }
        }
        a
    }

    /// FMM-accelerated matvec of the SAME operator (three gradient
    /// passes dotted with target normals; the P2P near field inside
    /// fs-fmm keeps close pairs exact at centroid resolution).
    #[must_use]
    pub fn fmm_matvec(&self, sigma: &[f64], order: usize) -> Vec<f64> {
        let n = self.centroids.len();
        let weighted: Vec<f64> = sigma.iter().zip(&self.areas).map(|(s, a)| s * a).collect();
        let mut out = vec![0.0f64; n];
        for c in 0..3 {
            let k = GradKernel { c };
            let fmm = Fmm::new(&k, self.centroids.clone(), order, 40);
            let comp = fmm.potentials(&weighted);
            for i in 0..n {
                out[i] += self.normals[i][c] * comp[i];
            }
        }
        for i in 0..n {
            out[i] += -0.5 * sigma[i];
        }
        out
    }
}

/// The GMRES operator wrapping the FMM matvec.
struct FmmOp<'a> {
    panels: &'a SpherePanels,
    order: usize,
}

impl LinearOp for FmmOp<'_> {
    fn n(&self) -> usize {
        self.panels.centroids.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let v = self.panels.fmm_matvec(x, self.order);
        y.copy_from_slice(&v);
    }
    fn apply_transpose(&self, _x: &[f64], _y: &mut [f64]) {
        unimplemented!("transpose path unused by GMRES here")
    }
}

/// Solve the exterior Neumann problem for uniform onset flow `u_inf`:
/// source densities σ with `A·σ = −u_inf·n`. Returns (σ, iterations,
/// relative residual).
#[must_use]
pub fn solve_exterior(
    panels: &SpherePanels,
    u_inf: [f64; 3],
    order: usize,
    tol: f64,
) -> (Vec<f64>, usize, f64) {
    let n = panels.centroids.len();
    let rhs: Vec<f64> = (0..n)
        .map(|i| {
            -(u_inf[0] * panels.normals[i][0]
                + u_inf[1] * panels.normals[i][1]
                + u_inf[2] * panels.normals[i][2])
        })
        .collect();
    let op = FmmOp { panels, order };
    let mut st = GmresState::new(&rhs, 60);
    let _ = st.run(&op, &rhs, tol, 8, false);
    (st.x.clone(), st.iters, st.rel_residual())
}

/// Surface velocity at panel centroids for a solved σ (onset +
/// induced; the tangential projection is the physical speed on the
/// body since the normal component vanishes by construction).
#[must_use]
pub fn surface_velocity(
    panels: &SpherePanels,
    sigma: &[f64],
    u_inf: [f64; 3],
    order: usize,
) -> Vec<[f64; 3]> {
    let n = panels.centroids.len();
    let weighted: Vec<f64> = sigma
        .iter()
        .zip(&panels.areas)
        .map(|(s, a)| s * a)
        .collect();
    let mut out = vec![u_inf; n];
    for c in 0..3 {
        let k = GradKernel { c };
        let fmm = Fmm::new(&k, panels.centroids.clone(), order, 40);
        let comp = fmm.potentials(&weighted);
        for (o, ci) in out.iter_mut().zip(&comp) {
            o[c] += ci;
        }
    }
    // Project out the (small residual) normal component.
    for (o, nrm) in out.iter_mut().zip(&panels.normals) {
        let vn = o[0] * nrm[0] + o[1] * nrm[1] + o[2] * nrm[2];
        for (oc, nc) in o.iter_mut().zip(nrm) {
            *oc -= vn * nc;
        }
    }
    out
}
