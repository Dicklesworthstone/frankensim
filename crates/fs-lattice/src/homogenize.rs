//! Periodic unit-cell homogenization (bead 7tv.14): the effective
//! elasticity of a microstructured cell from three cell problems.
//! Displacements split as u = E·x + u_per with u_per PERIODIC; the
//! cell stiffness comes from the fs-solid hyper2d tangent at u = 0
//! (exact linearization — no separate linear-elastic assembly to
//! drift), the periodic constraint is a representative-node reduction on the
//! structured mesh's exact opposite-edge node correspondence, and the
//! effective tensor is the energy average
//! C_hom[i][j] = (1/|Y|)·(E_i·x + χ_i)ᵀ K (E_j·x + χ_j).
//! Microstructure enters as a per-element density multiplier on the
//! stiffness (void ≈ `void_eps` — the classic contrast approach; the
//! CutFEM-exact variant is the recorded successor at this tier).

use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper2d::{HyperProblem, NewtonSettings};
use fs_solid::mesh2::Mesh2;

/// A 2D macroscopic strain (Voigt: xx, yy, xy with engineering shear).
pub const UNIT_STRAINS: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// The homogenized 2D elasticity tensor in Voigt form (3×3), plus its
/// provenance.
#[derive(Debug, Clone, Copy)]
pub struct EffectiveTensor {
    /// Voigt matrix (xx, yy, xy)².
    pub c: [[f64; 3]; 3],
    /// Solid volume fraction of the cell.
    pub density: f64,
}

/// A microstructured square unit cell: per-element solid indicator on
/// an `n × n` quad grid.
pub struct UnitCell {
    /// Elements per edge.
    pub n: usize,
    /// Per-element density in [0, 1] (1 = solid).
    pub rho: Vec<f64>,
}

impl UnitCell {
    /// A plate with a centered circular hole of radius `r` (cell is
    /// the unit square; element counted void when its center is
    /// inside the hole).
    #[must_use]
    pub fn holed_plate(n: usize, r: f64) -> UnitCell {
        let mut rho = vec![1.0f64; n * n];
        for j in 0..n {
            for i in 0..n {
                let cx = (i as f64 + 0.5) / n as f64 - 0.5;
                let cy = (j as f64 + 0.5) / n as f64 - 0.5;
                if cx.hypot(cy) < r {
                    rho[j * n + i] = 0.0;
                }
            }
        }
        UnitCell { n, rho }
    }

    /// A cross (plus-shaped) strut cell: arms of half-width `w` along
    /// both axes, void elsewhere.
    #[must_use]
    pub fn cross_strut(n: usize, w: f64) -> UnitCell {
        let mut rho = vec![0.0f64; n * n];
        for j in 0..n {
            for i in 0..n {
                let cx = (i as f64 + 0.5) / n as f64 - 0.5;
                let cy = (j as f64 + 0.5) / n as f64 - 0.5;
                if cx.abs() < w || cy.abs() < w {
                    rho[j * n + i] = 1.0;
                }
            }
        }
        UnitCell { n, rho }
    }

    /// Solid volume fraction.
    #[must_use]
    pub fn density(&self) -> f64 {
        self.rho.iter().sum::<f64>() / self.rho.len() as f64
    }
}

/// Homogenization workspace for one cell resolution: the base-cell
/// stiffness rows are assembled once per call (fixture scale — dense).
pub struct Homogenizer {
    /// The unit-square mesh.
    pub mesh: Mesh2,
    /// Base material.
    pub material: Hyperelastic,
    /// Void stiffness contrast (relative to solid).
    pub void_eps: f64,
}

impl Homogenizer {
    /// Standard plane-stress-ish base material (μ = 1, λ = 1.5 gives
    /// E ≈ 2.6, ν ≈ 0.3 in plane strain — the ORACLE gates use the
    /// measured solid-cell tensor as their reference, so the exact
    /// convention cancels).
    #[must_use]
    pub fn new(n: usize) -> Homogenizer {
        Homogenizer {
            mesh: Mesh2::quads(1.0, 1.0, n, n),
            material: Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 1.0,
                    lambda: 1.5,
                },
                1.5,
            )
            .expect("card"),
            void_eps: 1e-6,
        }
    }

    /// Effective tensor of `cell` by three periodic cell problems.
    ///
    /// # Panics
    /// On mesh/cell size mismatch or singular reduced systems
    /// (programmer contracts at fixture scale).
    #[must_use]
    #[allow(clippy::too_many_lines)] // assemble → reduce → solve → average, one narrative
    pub fn effective(&self, cell: &UnitCell) -> EffectiveTensor {
        let n = cell.n;
        let nn = n + 1;
        assert_eq!(self.mesh.node_count(), nn * nn, "mesh/cell mismatch");
        let ndof = 2 * nn * nn;
        // Element-wise density-scaled stiffness: assemble the SOLID
        // tangent once per element block via a per-element trick —
        // fixture-simple: assemble the full solid K, then rescale each
        // element's contribution. hyper2d has no per-element hook, so
        // assemble per-element by zeroing all other densities: O(n²)
        // assemblies would be wasteful — instead use the contrast
        // field DIRECTLY: assemble K(ρ) by scaling the MATERIAL per
        // element. hyper2d supports one material; the honest
        // fixture-scale route: assemble the solid K on the full mesh
        // and the per-element K_e via a single-element mesh, then sum
        // scaled scatters. The single reference element is exact
        // because the structured cells are congruent.
        let elem_mesh = Mesh2::quads(1.0 / n as f64, 1.0 / n as f64, 1, 1);
        let elem_problem = HyperProblem {
            mesh: &elem_mesh,
            material: &self.material,
            dirichlet: vec![],
            traction: vec![],
            settings: NewtonSettings::default(),
        };
        let (_, ke_csr) = elem_problem
            .residual_and_tangent(&[0.0; 8], 0.0)
            .expect("element tangent");
        let ke = ke_csr.to_dense(); // 8×8
        // Scatter: global K(ρ) dense.
        let mut k = vec![0.0f64; ndof * ndof];
        for j in 0..n {
            for i in 0..n {
                let e = j * n + i;
                let scale = cell.rho[e].max(self.void_eps);
                // Element nodes (Q1, CCW): matches Mesh2::quads local
                // ordering (id = i + j·(n+1)).
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
        // Periodic representative reduction: right edge → left edge, top →
        // bottom; corners all collapse to node 0. Map each node to its
        // retained representative.
        let mut representative: Vec<usize> = (0..nn * nn).collect();
        for j in 0..nn {
            representative[j * nn + n] = j * nn; // right → left
        }
        for i in 0..nn {
            representative[n * nn + i] = representative[i]; // top → bottom (post left-fold)
        }
        // Reduced dof indexing over representative nodes, minus rigid modes:
        // pin node 0 fully (translations; periodic + energy average is
        // rotation-free for the cell problems).
        let mut red_index = vec![usize::MAX; ndof];
        let mut nred = 0usize;
        for node in 0..nn * nn {
            if representative[node] == node && node != 0 {
                red_index[2 * node] = nred;
                red_index[2 * node + 1] = nred + 1;
                nred += 2;
            }
        }
        let red_of = |dof: usize| -> usize {
            let node = dof / 2;
            let comp = dof % 2;
            let m = representative[node];
            if m == 0 {
                usize::MAX
            } else {
                red_index[2 * m + comp]
            }
        };
        // Reduced stiffness.
        let mut kr = vec![0.0f64; nred * nred];
        for r in 0..ndof {
            let rr = red_of(r);
            if rr == usize::MAX {
                continue;
            }
            for c in 0..ndof {
                let cc = red_of(c);
                if cc == usize::MAX {
                    continue;
                }
                kr[rr * nred + cc] += k[r * ndof + c];
            }
        }
        let f = fs_la::factor::lu(&kr, nred).expect("reduced cell stiffness nonsingular");
        // Affine fields and cell-problem solutions.
        let affine = |voigt: [f64; 3]| -> Vec<f64> {
            let mut u = vec![0.0f64; ndof];
            for (node, p) in self.mesh.nodes.iter().enumerate() {
                u[2 * node] = voigt[0].mul_add(p[0], 0.5 * voigt[2] * p[1]);
                u[2 * node + 1] = voigt[1].mul_add(p[1], 0.5 * voigt[2] * p[0]);
            }
            u
        };
        let mut totals: Vec<Vec<f64>> = Vec::with_capacity(3);
        for voigt in UNIT_STRAINS {
            let ua = affine(voigt);
            // rhs_red = −(K ua) reduced.
            let mut kua = vec![0.0f64; ndof];
            for (r, out) in kua.iter_mut().enumerate() {
                let mut acc = 0.0;
                for c in 0..ndof {
                    acc += k[r * ndof + c] * ua[c];
                }
                *out = acc;
            }
            let mut rhs = vec![0.0f64; nred];
            for (r, &kr) in kua.iter().enumerate() {
                let rr = red_of(r);
                if rr != usize::MAX {
                    rhs[rr] -= kr;
                }
            }
            f.solve(&mut rhs);
            // Total field u = ua + P·u_per.
            let mut u = ua.clone();
            for (dof, val) in u.iter_mut().enumerate() {
                let rr = red_of(dof);
                if rr != usize::MAX {
                    *val += rhs[rr];
                }
            }
            totals.push(u);
        }
        // Energy averages: C[i][j] = u_iᵀ K u_j / |Y| (|Y| = 1).
        let mut c = [[0.0f64; 3]; 3];
        for i in 0..3 {
            let mut ku = vec![0.0f64; ndof];
            for r in 0..ndof {
                let mut acc = 0.0;
                for cc in 0..ndof {
                    acc += k[r * ndof + cc] * totals[i][cc];
                }
                ku[r] = acc;
            }
            for j in 0..3 {
                c[i][j] = totals[j].iter().zip(&ku).map(|(a, b)| a * b).sum();
            }
        }
        EffectiveTensor {
            c,
            density: cell.density(),
        }
    }
}

/// Voigt and Reuss bounds on the (xx, xx) stiffness entry for a
/// solid/void mixture at volume fraction `f_solid` — void makes Reuss
/// collapse to ~0; Voigt is the linear mixture (the hard upper bound
/// every homogenized entry must respect).
#[must_use]
pub fn voigt_bound(c_solid: f64, f_solid: f64, void_eps: f64) -> f64 {
    f_solid.mul_add(c_solid, (1.0 - f_solid) * void_eps * c_solid)
}
