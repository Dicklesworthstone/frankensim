//! Stage 1: PARAMETERIZE. The smoke-tier ornithoid candidate is a
//! sectional design: wing thickness, trim angle of attack, inlet
//! position along the chord, and a flapping gait (amplitude ×
//! reduced frequency). Every lever exposes a JACOBIAN ACTION — the
//! landed fs-bem adjoint gives a solve-based ∂cl/∂α with finite-difference
//! output partials; the remaining levers
//! carry deterministic central differences (documented — the F-rep
//! manifold-harmonics global form with analytic chart Jacobians is
//! the recorded successor).

use fs_bem::panel2d::{Airfoil2d, dcl_dalpha_adjoint, naca4_symmetric, solve};

/// Gene vector length (thickness, alpha, inlet_x, flap_amp, flap_freq).
pub const GENE_DIM: usize = 5;

/// A smoke-tier ornithoid candidate.
#[derive(Debug, Clone, Copy)]
pub struct OrnithCandidate {
    /// Wing section thickness ratio (NACA-4 symmetric family).
    pub thickness: f64,
    /// Trim angle of attack (radians).
    pub alpha: f64,
    /// Inlet position along the chord ∈ (0, 1).
    pub inlet_x: f64,
    /// Flapping amplitude (circulation modulation fraction).
    pub flap_amp: f64,
    /// Reduced flapping frequency.
    pub flap_freq: f64,
}

impl OrnithCandidate {
    /// Decode a gene vector in [0, 1]^5 into the bounded design box.
    #[must_use]
    pub fn from_genes(g: &[f64]) -> OrnithCandidate {
        assert_eq!(g.len(), GENE_DIM, "gene vector has {GENE_DIM} entries");
        let clamp01 = |v: f64| v.clamp(0.0, 1.0);
        OrnithCandidate {
            thickness: 0.06 + 0.12 * clamp01(g[0]),
            alpha: 0.02 + 0.12 * clamp01(g[1]),
            inlet_x: 0.15 + 0.6 * clamp01(g[2]),
            flap_amp: 0.05 + 0.4 * clamp01(g[3]),
            flap_freq: 0.2 + 1.4 * clamp01(g[4]),
        }
    }

    /// The wing section (unit chord, closed trailing edge).
    #[must_use]
    pub fn section(&self, panels: usize) -> Airfoil2d {
        naca4_symmetric(self.thickness, panels)
            .expect("ornith candidate section must satisfy the bounded NACA contract")
    }

    /// Inlet mass-flow satisfaction proxy: the panel tangential speed
    /// at the inlet station (suction feeds the inlet); target band
    /// [0.9, 1.6] of freestream.
    #[must_use]
    pub fn inlet_mass_flow(&self, panels: usize) -> f64 {
        let foil = self.section(panels);
        let sol = solve(&foil, self.alpha)
            .expect("ornith candidate must satisfy the bounded panel-solve contract");
        // Panel midpoints run TE → LE (lower) → TE (upper); find the
        // upper-surface panel nearest the inlet station.
        let n = foil.nodes().len();
        let mut best = (f64::INFINITY, 0usize);
        for i in 0..n {
            let a = foil.nodes()[i];
            let b = foil.nodes()[(i + 1) % n];
            let mid = [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])];
            if mid[1] > 0.0 {
                let d = (mid[0] - self.inlet_x).abs();
                if d < best.0 {
                    best = (d, i);
                }
            }
        }
        sol.vt[best.1].abs()
    }

    /// The Jacobian ACTION of cl w.r.t. the design levers: ∂cl/∂α by
    /// the landed BEM adjoint with solve-free finite-difference output
    /// partials, ∂cl/∂thickness by deterministic central difference.
    #[must_use]
    pub fn cl_gradient(&self, panels: usize) -> [f64; 2] {
        let foil = self.section(panels);
        let dcl_da = dcl_dalpha_adjoint(&foil, self.alpha)
            .expect("ornith candidate must satisfy the bounded adjoint contract");
        let h = 1e-5;
        let up = naca4_symmetric(self.thickness + h, panels)
            .expect("bounded positive thickness perturbation");
        let dn = naca4_symmetric(self.thickness - h, panels)
            .expect("bounded positive thickness perturbation");
        let dcl_dt = (solve(&up, self.alpha).expect("bounded panel solve").cl
            - solve(&dn, self.alpha).expect("bounded panel solve").cl)
            / (2.0 * h);
        [dcl_da, dcl_dt]
    }
}
