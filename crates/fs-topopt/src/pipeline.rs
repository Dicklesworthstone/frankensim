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
    /// Refuse parameter sets that do not define the documented SIMP model.
    pub(crate) fn assert_valid(self) {
        assert!(
            self.e_min.is_finite() && self.e_min > 0.0 && self.e_min < 1.0,
            "SIMP void-modulus floor must be finite and lie in (0, 1)"
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

/// One independently applied load, not a term in a simultaneous force sum.
#[derive(Debug, Clone, Copy)]
pub struct LoadCase<'a> {
    /// One finite force per vector dof. Loads on fixed dofs are reactions
    /// and do not do displacement work under homogeneous Dirichlet conditions.
    pub force: &'a [f64],
    /// Nonnegative objective weight. Weights are not implicitly normalized.
    pub weight: f64,
}

/// Weighted compliance and exact design sensitivity of independent load cases.
#[derive(Debug, Clone)]
pub struct MultiLoadCompliance {
    /// Sum of weight times compliance, with no cross-load terms.
    pub compliance: f64,
    /// Unweighted compliance for each load, in input order.
    pub case_compliances: Vec<f64>,
    /// Displacement for each load, including exact zeros at fixed dofs.
    pub displacements: Vec<Vec<f64>>,
    /// Gradient of the weighted objective with respect to raw densities.
    pub gradient: Vec<f64>,
}

pub(crate) fn assert_valid_load_cases(elasticity: &DensityElasticity, loads: &[LoadCase<'_>]) {
    assert!(!loads.is_empty(), "at least one load case is required");
    assert!(
        loads.iter().any(|load| load.weight > 0.0),
        "at least one load weight must be positive"
    );
    for load in loads {
        assert!(
            load.weight.is_finite() && load.weight >= 0.0,
            "load weights must be finite and nonnegative"
        );
        assert_eq!(
            load.force.len(),
            elasticity.n(),
            "each load must contain one force per elasticity dof"
        );
        assert!(
            load.force.iter().all(|value| value.is_finite()),
            "load forces must be finite"
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

    /// Projected physical volume fraction and its raw-design gradient.
    /// The derivative is Fᵀ(H′ ⊙ V / ΣV), NOT the raw cell-volume
    /// vector and NOT the SIMP stiffness pullback.
    #[must_use]
    pub fn volume_and_gradient(&self, rho: &[f64], cell_vol: &[f64]) -> (f64, Vec<f64>) {
        assert_eq!(rho.len(), cell_vol.len(), "one volume per design cell is required");
        assert!(
            !cell_vol.is_empty()
                && cell_vol.iter().all(|volume| volume.is_finite() && *volume > 0.0),
            "cell volumes must be finite, positive, and nonempty"
        );
        let total: f64 = cell_vol.iter().sum();
        assert!(total.is_finite(), "total cell volume must be finite");
        let (filtered, projected, _) = self.forward(rho);
        let weights: Vec<f64> = cell_vol.iter().map(|volume| volume / total).collect();
        let volume = projected.iter().zip(&weights).map(|(r, w)| r * w).sum();
        let local: Vec<f64> = filtered
            .iter()
            .zip(&weights)
            .map(|(&r, &w)| w * heaviside_derivative(r, self.params.beta, self.params.eta))
            .collect();
        (volume, self.filter.apply_transpose(&local))
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
        let energies = elasticity.cell_energies(&u);
        let dc_de: Vec<f64> = energies.iter().map(|e| -e).collect();
        let grad = self.pullback(&rho_tilde, &dc_de);
        (compliance, u, grad)
    }

    /// Evaluate independent loading scenarios against ONE design/operator.
    /// Each load gets its own equilibrium solve. Weighted element energies
    /// are accumulated before the single filter-transpose solve. Summing
    /// forces instead would introduce cross terms and can erase opposite loads.
    pub fn multi_load_compliance_and_gradient(
        &self,
        elasticity: &mut DensityElasticity,
        rho: &[f64],
        loads: &[LoadCase<'_>],
    ) -> MultiLoadCompliance {
        assert_valid_load_cases(elasticity, loads);
        assert_eq!(elasticity.cells(), rho.len(), "one density per cell is required");
        let (filtered, _, moduli) = self.forward(rho);
        elasticity.moduli = moduli;
        let mut compliance = 0.0;
        let mut dc_de = vec![0.0; rho.len()];
        let mut case_compliances = Vec::with_capacity(loads.len());
        let mut displacements = Vec::with_capacity(loads.len());
        for load in loads {
            let u = solve(elasticity, load.force);
            let c: f64 = load.force.iter().zip(&u).map(|(f, value)| f * value).sum();
            assert!(c.is_finite(), "load compliance must remain finite");
            compliance += load.weight * c;
            for (sensitivity, energy) in dc_de.iter_mut().zip(elasticity.cell_energies(&u)) {
                *sensitivity -= load.weight * energy;
            }
            case_compliances.push(c);
            displacements.push(u);
        }
        assert!(compliance.is_finite(), "weighted compliance must remain finite");
        let gradient = self.pullback(&filtered, &dc_de);
        MultiLoadCompliance { compliance, case_compliances, displacements, gradient }
    }
}

fn solve(op: &DensityElasticity, b: &[f64]) -> Vec<f64> {
    // The operator uses identity rows at fixed dofs to remain SPD. Those
    // rows enforce u=0, not u=f: support loads belong to the reaction balance.
    let rhs: Vec<f64> = b.iter().zip(op.free()).map(|(&f, &free)| if free { f } else { 0.0 }).collect();
    if rhs.iter().all(|value| value.abs() <= 0.0) {
        return vec![0.0; op.n()];
    }
    let mut st = fs_solver::CgState::new(op, &fs_sparse::precond::IdentityPrecond, &rhs);
    let rep = st.run(op, &fs_sparse::precond::IdentityPrecond, 1e-11, 50_000);
    assert!(rep.converged, "elasticity solve failed: {rep:?}");
    st.x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>, Vec<f64>) {
        let (complex, positions) = fs_feec::kuhn_cube(2);
        let elasticity = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p| p[0] < 1e-12);
        let mut force = vec![0.0; elasticity.n()];
        for (v, p) in positions.iter().enumerate() {
            if p[0] > 1.0 - 1e-12 {
                force[3 * v + 2] = -1.0;
            }
        }
        let volumes = fs_feec::element_geometry(&complex, &positions)
            .vol_signed.iter().map(|v| v.abs()).collect();
        let rho = (0..elasticity.cells()).map(|i| 0.35 + 0.02 * (i % 11) as f64).collect();
        let pipeline = DesignPipeline {
            filter: DensityFilter::new(&complex, &positions, 0.15),
            params: SimpParams::default(),
        };
        (pipeline, elasticity, rho, force, volumes)
    }

    #[test]
    fn opposite_load_cases_do_not_cancel() {
        let (pipeline, mut elasticity, rho, force, _) = fixture();
        let opposite: Vec<f64> = force.iter().map(|f| -f).collect();
        let (single, _, grad) = pipeline.compliance_and_gradient(&mut elasticity, &rho, &force);
        let result = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &rho, &[
            LoadCase { force: &force, weight: 0.25 },
            LoadCase { force: &opposite, weight: 0.75 },
        ]);
        assert!(result.compliance > 0.0);
        assert!((result.compliance - single).abs() < 1e-10 * single);
        for (actual, expected) in result.gradient.iter().zip(&grad) {
            assert!((actual - expected).abs() < 1e-8 * expected.abs().max(1.0));
        }
        for (a, b) in result.displacements[0].iter().zip(&result.displacements[1]) {
            assert!((a + b).abs() < 1e-9 * a.abs().max(1.0));
        }
    }

    #[test]
    fn support_loads_do_not_create_displacement_or_compliance() {
        let (pipeline, mut elasticity, rho, force, _) = fixture();
        let (expected, _, _) = pipeline.compliance_and_gradient(&mut elasticity, &rho, &force);
        let mut loaded = force;
        for (f, free) in loaded.iter_mut().zip(elasticity.free()) {
            if !free { *f = 1000.0; }
        }
        let (actual, u, _) = pipeline.compliance_and_gradient(&mut elasticity, &rho, &loaded);
        assert!((actual - expected).abs() < 1e-10 * expected);
        for (value, free) in u.iter().zip(elasticity.free()) {
            if !free { assert!(value.abs() <= 0.0); }
        }
        for (f, free) in loaded.iter_mut().zip(elasticity.free()) {
            if *free { *f = 0.0; }
        }
        let (c, u, g) = pipeline.compliance_and_gradient(&mut elasticity, &rho, &loaded);
        assert!(c.abs() <= 0.0);
        assert!(u.iter().chain(&g).all(|v| v.abs() <= 0.0));
    }

    #[test]
    fn independent_load_gradient_matches_finite_differences() {
        let (mut pipeline, mut elasticity, rho, force, _) = fixture();
        let mut lateral = vec![0.0; force.len()];
        for (src, dst) in force.chunks_exact(3).zip(lateral.chunks_exact_mut(3)) {
            dst[1] = src[2];
        }
        let loads = [LoadCase { force: &force, weight: 0.3 }, LoadCase { force: &lateral, weight: 0.7 }];
        for beta in [0.0, 2.0, 8.0] {
            pipeline.params.beta = beta;
            let result = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &rho, &loads);
            for i in [0, 5, 12, 43] {
                let mut plus = rho.clone();
                let mut minus = rho.clone();
                plus[i] += 1e-5;
                minus[i] -= 1e-5;
                let cp = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &plus, &loads).compliance;
                let cm = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &minus, &loads).compliance;
                let fd = (cp - cm) / 2e-5;
                assert!((fd - result.gradient[i]).abs() < 2e-5 * fd.abs().max(1.0), "beta={beta}, cell={i}, fd={fd}, analytic={}", result.gradient[i]);
            }
        }
    }

    #[test]
    fn physical_volume_gradient_includes_filter_and_projection() {
        let (mut pipeline, _, rho, _, volumes) = fixture();
        pipeline.params.beta = 8.0;
        let (_, gradient) = pipeline.volume_and_gradient(&rho, &volumes);
        for i in [0, 5, 12, 43] {
            let mut plus = rho.clone();
            let mut minus = rho.clone();
            plus[i] += 1e-5;
            minus[i] -= 1e-5;
            let vp = pipeline.volume_and_gradient(&plus, &volumes).0;
            let vm = pipeline.volume_and_gradient(&minus, &volumes).0;
            assert!(((vp - vm) / 2e-5 - gradient[i]).abs() < 1e-7);
        }
    }
}
