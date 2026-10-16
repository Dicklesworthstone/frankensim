//! Graded macro-optimization + de-homogenization (bead 7tv.14, smoke
//! tier): a cantilever whose per-element CELL DENSITY is the design
//! field, elasticity through the FITTED homogenized law C(ρ) (the
//! property manifold from `homogenize`, monotone cubic fit with a
//! declared validity domain), compliance minimized under a volume
//! budget by projected adjoint gradient descent. The SCALE-SEPARATION
//! flag fires when neighboring cells' densities jump beyond the
//! declared gradation bound — homogenization validity is a model
//! card, not an assumption.

use crate::homogenize::{Homogenizer, UnitCell};
use fs_solid::hyper2d::{HyperProblem, NewtonSettings};
use fs_solid::mesh2::Mesh2;

/// A fitted scalar property curve s(ρ) (the C11 entry, normalized by
/// the solid value) with its sampling provenance and validity domain.
pub struct PropertyFit {
    /// Sample densities.
    pub rho: Vec<f64>,
    /// Normalized stiffness samples s(ρ) ∈ [0, 1].
    pub s: Vec<f64>,
    /// Max |Δρ| between adjacent macro cells for which the
    /// separation-of-scales card holds (declared, checked by the
    /// optimizer's flag).
    pub gradation_bound: f64,
}

impl PropertyFit {
    /// Build the manifold by homogenizing holed-plate cells across a
    /// density sweep (n × n cell resolution).
    ///
    /// # Panics
    /// On degenerate sweeps (programmer contract).
    #[must_use]
    pub fn sample_holed_plates(n: usize, radii: &[f64]) -> PropertyFit {
        let hom = Homogenizer::new(n);
        let solid = hom.effective(&UnitCell::holed_plate(n, 0.0));
        let mut rho = Vec::with_capacity(radii.len());
        let mut s = Vec::with_capacity(radii.len());
        for &r in radii {
            let cell = UnitCell::holed_plate(n, r);
            let eff = hom.effective(&cell);
            rho.push(eff.density);
            s.push(eff.c[0][0] / solid.c[0][0]);
        }
        PropertyFit {
            rho,
            s,
            gradation_bound: 0.35,
        }
    }

    /// Piecewise-linear interpolation of s(ρ), clamped to the sampled
    /// range (the validity domain in ρ).
    #[must_use]
    pub fn eval(&self, rho: f64) -> f64 {
        let n = self.rho.len();
        // Samples are stored in DESCENDING density (growing hole).
        if rho >= self.rho[0] {
            return self.s[0];
        }
        if rho <= self.rho[n - 1] {
            return self.s[n - 1];
        }
        for w in 0..n - 1 {
            let (r0, r1) = (self.rho[w], self.rho[w + 1]);
            if rho <= r0 && rho >= r1 {
                let t = (r0 - rho) / (r0 - r1);
                return t.mul_add(self.s[w + 1] - self.s[w], self.s[w]);
            }
        }
        self.s[n - 1]
    }

    /// d s / d ρ by the same piecewise-linear rule.
    #[must_use]
    pub fn slope(&self, rho: f64) -> f64 {
        let n = self.rho.len();
        if rho >= self.rho[0] || rho <= self.rho[n - 1] {
            return 0.0;
        }
        for w in 0..n - 1 {
            let (r0, r1) = (self.rho[w], self.rho[w + 1]);
            if rho <= r0 && rho >= r1 {
                return (self.s[w] - self.s[w + 1]) / (r0 - r1);
            }
        }
        0.0
    }
}

/// The graded design record.
pub struct GradedDesign {
    /// Per-element densities.
    pub rho: Vec<f64>,
    /// Final compliance.
    pub compliance: f64,
    /// Compliance of the equal-mass uniform baseline.
    pub uniform_compliance: f64,
    /// Iterations used.
    pub iterations: u32,
    /// Did any adjacent-cell density jump exceed the fit's declared
    /// gradation bound? (The separation-of-scales honesty flag.)
    pub scale_separation_violated: bool,
    /// Worst adjacent jump observed.
    pub worst_jump: f64,
}

/// Assemble the cantilever stiffness at density field `rho` through
/// the fitted law and return (compliance, load vector, K dense).
fn analyze(
    mesh: &Mesh2,
    ke: &[f64],
    nx: usize,
    ny: usize,
    fit: &PropertyFit,
    rho: &[f64],
    floor: f64,
) -> (f64, Vec<f64>, Vec<f64>) {
    let nn = nx + 1;
    let ndof = 2 * mesh.node_count();
    let mut k = vec![0.0f64; ndof * ndof];
    for j in 0..ny {
        for i in 0..nx {
            let e = j * nx + i;
            let scale = fit.eval(rho[e]).max(floor);
            // ke index order: element-mesh GLOBAL ids (bl, br, tl, tr);
            // see homogenize.rs for the measured CCW-swap failure.
            let nodes = [
                j * nn + i,
                j * nn + i + 1,
                (j + 1) * nn + i,
                (j + 1) * nn + i + 1,
            ];
            for (a, &na) in nodes.iter().enumerate() {
                for (b, &nb) in nodes.iter().enumerate() {
                    for da in 0..2 {
                        for db in 0..2 {
                            k[(2 * na + da) * ndof + (2 * nb + db)] +=
                                scale * ke[(2 * a + da) * 8 + (2 * b + db)];
                        }
                    }
                }
            }
        }
    }
    // Cantilever: clamp the left edge, uniform downward traction on
    // the right edge with trapezoidal lumping (interior nodes h,
    // corner nodes h/2) — the SAME load convention the CutFEM
    // de-homogenization re-analysis applies, so compliances compare.
    let mut load = vec![0.0f64; ndof];
    let hy = 1.0 / ny as f64;
    for node in 0..mesh.node_count() {
        if (mesh.nodes[node][0] - 1.0).abs() < 1e-9 {
            let y = mesh.nodes[node][1];
            let w = if y < 1e-9 || (y - 1.0).abs() < 1e-9 {
                0.5 * hy
            } else {
                hy
            };
            load[2 * node + 1] = -w;
        }
    }
    for node in 0..mesh.node_count() {
        if mesh.nodes[node][0] <= 1e-12 {
            for comp in 0..2 {
                let d = 2 * node + comp;
                for c in 0..ndof {
                    k[d * ndof + c] = 0.0;
                    k[c * ndof + d] = 0.0;
                }
                k[d * ndof + d] = 1.0;
                load[d] = 0.0;
            }
        }
    }
    (0.0, load, k)
}

/// Projected-gradient compliance minimization at a volume budget.
/// Compliance c = fᵀu; dc/dρ_e = −s′(ρ_e)·uᵀK_e u (self-adjoint).
///
/// # Panics
/// On singular systems (programmer contracts at fixture scale).
#[must_use]
#[allow(clippy::too_many_lines)] // analyze → gradient → project loop, one narrative
pub fn graded_compliance_opt(
    nx: usize,
    ny: usize,
    cell_res: usize,
    volfrac: f64,
    iters: u32,
) -> GradedDesign {
    let mesh = Mesh2::quads(1.0, 1.0, nx, ny);
    // Reference element stiffness (congruent structured cells).
    let elem_mesh = Mesh2::quads(1.0 / nx as f64, 1.0 / ny as f64, 1, 1);
    let hom = Homogenizer::new(cell_res);
    let elem_problem = HyperProblem {
        mesh: &elem_mesh,
        material: &hom.material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings::default(),
    };
    let (_, ke_csr) = elem_problem
        .residual_and_tangent(&[0.0; 8], 0.0)
        .expect("element tangent");
    let ke = ke_csr.to_dense();
    let fit = PropertyFit::sample_holed_plates(cell_res, &[0.0, 0.15, 0.25, 0.35, 0.45]);
    let ne = nx * ny;
    let floor = 1e-4;
    let solve = |rho: &[f64]| -> (f64, Vec<f64>) {
        let (_, load, k) = analyze(&mesh, &ke, nx, ny, &fit, rho, floor);
        let ndof = load.len();
        let f = fs_la::factor::lu(&k, ndof).expect("stiffness nonsingular");
        let mut u = load.clone();
        f.solve(&mut u);
        let c: f64 = load.iter().zip(&u).map(|(a, b)| a * b).sum();
        (c, u)
    };
    // Uniform baseline at the same mass.
    let rho_min = *fit.rho.last().expect("samples") + 1e-6;
    let uniform = vec![volfrac.clamp(rho_min, 1.0); ne];
    let (uniform_compliance, _) = solve(&uniform);
    // Projected gradient with a bisected volume multiplier.
    let mut rho = uniform.clone();
    let nn = nx + 1;
    for _ in 0..iters {
        let (_, u) = solve(&rho);
        // Element sensitivities: −s′(ρ)·uᵀ K_e u.
        let mut sens = vec![0.0f64; ne];
        for j in 0..ny {
            for i in 0..nx {
                let e = j * nx + i;
                // ke is indexed by the ELEMENT MESH's GLOBAL node ids
                // (bl, br, tl, tr) — NOT CCW element order; the CCW
                // ordering was measured to swap the two top nodes and
                // silently relax C11 from 3.5 to 1.5 while leaving
                // C22 exact (y-fields cannot see a top-node swap).
                let nodes = [
                    j * nn + i,
                    j * nn + i + 1,
                    (j + 1) * nn + i,
                    (j + 1) * nn + i + 1,
                ];
                let mut energy = 0.0;
                for (a, &na) in nodes.iter().enumerate() {
                    for (b, &nb) in nodes.iter().enumerate() {
                        for da in 0..2 {
                            for db in 0..2 {
                                energy += u[2 * na + da]
                                    * ke[(2 * a + da) * 8 + (2 * b + db)]
                                    * u[2 * nb + db];
                            }
                        }
                    }
                }
                sens[e] = -fit.slope(rho[e]) * energy;
            }
        }
        // OC-style multiplicative update with volume bisection.
        let (mut lo, mut hi) = (1e-9f64, 1e9f64);
        let mut trial = rho.clone();
        for _ in 0..60 {
            let lam = (lo * hi).sqrt();
            for e in 0..ne {
                let b = (-sens[e] / lam).max(0.0).sqrt();
                // Tight move limits: 0.7/1.3 was measured to oscillate
                // past the optimum on the piecewise-linear fit.
                trial[e] = (rho[e] * b.clamp(0.9, 1.1)).clamp(rho_min, 1.0);
            }
            let vol: f64 = trial.iter().sum::<f64>() / ne as f64;
            if vol > volfrac {
                lo = lam;
            } else {
                hi = lam;
            }
        }
        rho = trial;
    }
    // Final state evaluated AFTER the last update (the pre-update
    // value was measured stale by one oscillating iteration).
    let (compliance, _) = solve(&rho);
    // Scale-separation audit: adjacent density jumps vs the declared
    // gradation bound.
    let mut worst_jump = 0.0f64;
    for j in 0..ny {
        for i in 0..nx {
            let e = j * nx + i;
            if i + 1 < nx {
                worst_jump = worst_jump.max((rho[e] - rho[e + 1]).abs());
            }
            if j + 1 < ny {
                worst_jump = worst_jump.max((rho[e] - rho[e + nx]).abs());
            }
        }
    }
    GradedDesign {
        scale_separation_violated: worst_jump > fit.gradation_bound,
        worst_jump,
        rho,
        compliance,
        uniform_compliance,
        iterations: iters,
    }
}
