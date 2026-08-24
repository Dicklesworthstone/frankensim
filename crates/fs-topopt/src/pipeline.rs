//! The design pipeline ρ → ρ̃ (Helmholtz filter) → ρ̄ (Heaviside
//! projection) → E(ρ̄) (SIMP), with the EXACT reverse chain
//! dc/dρ = filterᵀ(projection′ ⊙ SIMP′ ⊙ dc/dE-contraction) —
//! FD-verified at multiple continuation stages in the battery, per
//! the acceptance.

use crate::elasticity::DensityElasticity;
use crate::filter::{DensityFilter, heaviside, heaviside_derivative};

/// SIMP + continuation parameters.
#[derive(Debug, Clone, Copy)]
pub struct SimpParams {
    /// Void modulus floor (relative to E₀ = 1).
    pub e_min: f64,
    /// Penalization exponent p.
    pub penal: f64,
    /// Heaviside sharpness β.
    pub beta: f64,
    /// Heaviside threshold η.
    pub eta: f64,
}

impl Default for SimpParams {
    fn default() -> SimpParams {
        SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
        }
    }
}

impl SimpParams {
    pub(crate) fn assert_valid(self) {
        assert!(
            self.e_min.is_finite() && self.e_min > 0.0 && self.e_min <= 1.0,
            "SIMP void-modulus floor must be finite and lie in (0, 1]"
        );
        assert!(
            self.penal.is_finite() && self.penal >= 1.0,
            "SIMP penalization exponent must be finite and at least one"
        );
        assert!(
            self.beta.is_finite() && self.beta >= 0.0,
            "Heaviside sharpness beta must be finite and nonnegative"
        );
        assert!(
            self.eta.is_finite() && (0.0..=1.0).contains(&self.eta),
            "Heaviside threshold eta must be finite and lie in [0, 1]"
        );
    }
}

/// The chained design-to-physics pipeline.
pub struct DesignPipeline {
    /// The filter stage.
    pub filter: DensityFilter,
    /// SIMP/projection parameters (continuation mutates these).
    pub params: SimpParams,
}

impl DesignPipeline {
    /// Forward: raw design ρ → (ρ̃ filtered, ρ̄ projected, E moduli).
    #[must_use]
    pub fn forward(&self, rho: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let p = &self.params;
        p.assert_valid();
        assert!(
            rho.iter()
                .all(|value| value.is_finite() && (0.0..=1.0).contains(value)),
            "raw design densities must be finite and lie in [0, 1]"
        );
        let rho_tilde = self.filter.apply(rho);
        let rho_bar: Vec<f64> = rho_tilde
            .iter()
            .map(|&r| heaviside(r, p.beta, p.eta))
            .collect();
        let moduli: Vec<f64> = rho_bar
            .iter()
            .map(|&r| {
                let rc = r.clamp(0.0, 1.0);
                p.e_min + (1.0 - p.e_min) * fs_math::det::pow(rc.max(1e-12), p.penal)
            })
            .collect();
        (rho_tilde, rho_bar, moduli)
    }

    /// Reverse chain: given dc/dE per cell (the physics-level
    /// sensitivity), pull back to dc/dρ through SIMP′, projection′,
    /// and the transposed filter.
    #[must_use]
    pub fn pullback(&self, rho_tilde: &[f64], dc_de: &[f64]) -> Vec<f64> {
        let p = &self.params;
        p.assert_valid();
        assert_eq!(
            rho_tilde.len(),
            dc_de.len(),
            "SIMP pullback requires one physics sensitivity per filtered density"
        );
        assert!(
            dc_de.iter().all(|value| value.is_finite()),
            "SIMP pullback sensitivities must be finite"
        );
        let chained: Vec<f64> = rho_tilde
            .iter()
            .zip(dc_de)
            .map(|(&rt, &de)| {
                let rb = heaviside(rt, p.beta, p.eta).clamp(0.0, 1.0);
                let dsimp =
                    (1.0 - p.e_min) * p.penal * fs_math::det::pow(rb.max(1e-12), p.penal - 1.0);
                let dproj = heaviside_derivative(rt, p.beta, p.eta);
                de * dsimp * dproj
            })
            .collect();
        self.filter.apply_transpose(&chained)
    }

    /// Compliance objective and its EXACT design gradient for the
    /// elasticity problem: c = fᵀu (self-adjoint: λ = u, so
    /// dc/dE_c = −u_cᵀK_cu_c — no extra solve), then the reverse
    /// chain. Returns (compliance, u, dc/dρ).
    pub fn compliance_and_gradient(
        &self,
        elasticity: &mut DensityElasticity,
        rho: &[f64],
        force: &[f64],
    ) -> (f64, Vec<f64>, Vec<f64>) {
        assert_eq!(
            elasticity.cells(),
            rho.len(),
            "elasticity and design must have the same cell count"
        );
        assert_eq!(
            force.len(),
            elasticity.n(),
            "force must contain one finite value per elasticity dof"
        );
        assert!(
            force.iter().all(|value| value.is_finite()),
            "force must contain only finite values"
        );
        let (rho_tilde, _rho_bar, moduli) = self.forward(rho);
        elasticity.moduli = moduli;
        let u = solve(elasticity, force);
        let compliance: f64 = force.iter().zip(&u).map(|(f, ui)| f * ui).sum();
        // dc/dE_c = −u K_c u (unit-modulus cell energies).
        let energies = elasticity.cell_energies(&u);
        let dc_de: Vec<f64> = energies.iter().map(|e| -e).collect();
        let grad = self.pullback(&rho_tilde, &dc_de);
        (compliance, u, grad)
    }
}

fn solve(op: &DensityElasticity, b: &[f64]) -> Vec<f64> {
    let mut st = fs_solver::CgState::new(op, &fs_sparse::precond::IdentityPrecond, b);
    let rep = st.run(op, &fs_sparse::precond::IdentityPrecond, 1e-11, 50_000);
    assert!(rep.converged, "elasticity solve failed: {rep:?}");
    st.x
}
