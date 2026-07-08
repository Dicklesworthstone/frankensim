//! Robust erode/dilate three-field formulation (the
//! manufacturable-by-construction trick): project the SAME filtered
//! field at three thresholds — eroded (η_e > ½), nominal (½), dilated
//! (η_d < ½) — minimize the ERODED compliance (it dominates the
//! worst case for stiffness problems) subject to volume on the
//! DILATED design (the standard practical form, Wang–Lazarov–Sigmund
//! class). Because the projection is pointwise-decreasing in η, the
//! three realizations are POINTWISE ORDERED (ρ̄_e ≤ ρ̄_n ≤ ρ̄_d —
//! tested, not assumed), and surviving erosion is precisely the
//! minimum-length-scale guarantee: features thinner than the
//! erode/dilate band get annihilated by the eroded projection and
//! therefore cannot carry the optimized load path.

use crate::elasticity::DensityElasticity;
use crate::filter::{DensityFilter, heaviside, heaviside_derivative};
use crate::pipeline::SimpParams;

/// The three-field design state at one evaluation.
#[derive(Debug, Clone)]
pub struct ThreeField {
    /// Filtered densities ρ̃.
    pub rho_tilde: Vec<f64>,
    /// Eroded projection (η_e).
    pub eroded: Vec<f64>,
    /// Nominal projection (η = ½).
    pub nominal: Vec<f64>,
    /// Dilated projection (η_d).
    pub dilated: Vec<f64>,
}

/// Robust three-field pipeline: one filter, three projections.
pub struct RobustPipeline {
    /// The shared filter.
    pub filter: DensityFilter,
    /// SIMP parameters (β/penal continuation mutates these; η is the
    /// NOMINAL threshold, with `eta_offset` giving η_e/η_d = η ± δ).
    pub params: SimpParams,
    /// Threshold offset δ (η_e = η + δ, η_d = η − δ).
    pub eta_offset: f64,
}

impl RobustPipeline {
    /// All three projected fields from a raw design.
    #[must_use]
    pub fn three_fields(&self, rho: &[f64]) -> ThreeField {
        let p = &self.params;
        let rho_tilde = self.filter.apply(rho);
        let project = |eta: f64| -> Vec<f64> {
            rho_tilde.iter().map(|&r| heaviside(r, p.beta, eta)).collect()
        };
        ThreeField {
            eroded: project(p.eta + self.eta_offset),
            nominal: project(p.eta),
            dilated: project(p.eta - self.eta_offset),
            rho_tilde,
        }
    }

    /// SIMP moduli from a projected field.
    fn moduli(&self, rho_bar: &[f64]) -> Vec<f64> {
        let p = &self.params;
        rho_bar
            .iter()
            .map(|&r| {
                let rc = r.clamp(0.0, 1.0);
                p.e_min + (1.0 - p.e_min) * fs_math::det::pow(rc.max(1e-12), p.penal)
            })
            .collect()
    }

    /// ERODED compliance and its exact design gradient (the robust
    /// objective): same self-adjoint structure as the nominal path,
    /// chained through the eroded projection's slope.
    pub fn eroded_compliance_and_gradient(
        &self,
        elasticity: &mut DensityElasticity,
        rho: &[f64],
        force: &[f64],
    ) -> (f64, Vec<f64>) {
        let p = &self.params;
        let eta_e = p.eta + self.eta_offset;
        let tf = self.three_fields(rho);
        elasticity.moduli = self.moduli(&tf.eroded);
        let u = solve(elasticity, force);
        let compliance: f64 = force.iter().zip(&u).map(|(f, ui)| f * ui).sum();
        let energies = elasticity.cell_energies(&u);
        let chained: Vec<f64> = tf
            .rho_tilde
            .iter()
            .zip(&energies)
            .map(|(&rt, &e)| {
                let rb = heaviside(rt, p.beta, eta_e).clamp(0.0, 1.0);
                let dsimp = (1.0 - p.e_min)
                    * p.penal
                    * fs_math::det::pow(rb.max(1e-12), p.penal - 1.0);
                let dproj = heaviside_derivative(rt, p.beta, eta_e);
                -e * dsimp * dproj
            })
            .collect();
        (compliance, self.filter.apply_transpose(&chained))
    }
}

/// Outcome of a robust OC run.
#[derive(Debug, Clone)]
pub struct RobustReport {
    /// Final raw design.
    pub rho: Vec<f64>,
    /// Eroded-compliance trace.
    pub compliance_eroded: Vec<f64>,
    /// Final volume fractions (eroded, nominal, dilated).
    pub volumes: (f64, f64, f64),
    /// Erosion retention vol(eroded)/vol(nominal) — the measured
    /// minimum-length-scale signal (features thinner than the band
    /// die under erosion and drag this ratio down).
    pub erosion_retention: f64,
}

fn volume_fraction(field: &[f64], cell_vol: &[f64]) -> f64 {
    let total: f64 = cell_vol.iter().sum();
    field.iter().zip(cell_vol).map(|(r, v)| r * v).sum::<f64>() / total
}

/// Robust OC: minimize eroded compliance with the volume constraint
/// on the DILATED field, whose target is ADAPTED each iteration so
/// the NOMINAL design meets `vol_frac` (the standard
/// Wang–Lazarov–Sigmund practice — without adaptation the nominal
/// budget drifts far below target and the robust design is starved
/// relative to a non-robust baseline; the first draft measured
/// exactly that: retention 0.545 vs 0.622 AGAINST the robust run
/// purely from the budget mismatch).
#[allow(clippy::too_many_arguments)]
pub fn robust_optimality_criteria(
    pipeline: &RobustPipeline,
    elasticity: &mut DensityElasticity,
    force: &[f64],
    rho0: &[f64],
    cell_vol: &[f64],
    vol_frac: f64,
    move_limit: f64,
    iters: usize,
) -> RobustReport {
    let nc = rho0.len();
    let mut rho = rho0.to_vec();
    let mut trace = Vec::with_capacity(iters);
    let mut dilated_target = vol_frac;
    for _ in 0..iters {
        // Adapt the dilated target so the NOMINAL field hits vol_frac —
        // DAMPED (0.3 blend toward the instantaneous estimate). The
        // undamped update measured a period-2 limit cycle on the
        // uniform cold start (0.4 ↔ 0.6 bouncing, targets 1.0 ↔ 0.5
        // exactly out of phase with the move-limit-saturated steps);
        // damping contracts the cycle to the fixed point vn = vol_frac.
        let tf_now = pipeline.three_fields(&rho);
        let vn_now = volume_fraction(&tf_now.nominal, cell_vol);
        let vd_now = volume_fraction(&tf_now.dilated, cell_vol);
        if vn_now > 1e-12 {
            let instant = (vol_frac * vd_now / vn_now).clamp(vol_frac, 1.0);
            dilated_target = 0.3f64.mul_add(instant - dilated_target, dilated_target);
        }
        let (c, grad) = pipeline.eroded_compliance_and_gradient(elasticity, &rho, force);
        trace.push(c);
        let sensitivity: Vec<f64> = grad.iter().map(|g| (-g).max(1e-30)).collect();
        let mut lo = 1e-12f64;
        let mut hi = 1e12f64;
        let mut candidate = rho.clone();
        for _ in 0..80 {
            let lambda = fs_math::det::sqrt(lo * hi);
            for i in 0..nc {
                let scale = fs_math::det::sqrt(sensitivity[i] / (lambda * cell_vol[i]));
                candidate[i] = (rho[i] * scale)
                    .clamp(rho[i] - move_limit, rho[i] + move_limit)
                    .clamp(1e-3, 1.0);
            }
            let dilated = pipeline.three_fields(&candidate).dilated;
            if volume_fraction(&dilated, cell_vol) > dilated_target {
                lo = lambda;
            } else {
                hi = lambda;
            }
        }
        rho = candidate;
    }
    let tf = pipeline.three_fields(&rho);
    let ve = volume_fraction(&tf.eroded, cell_vol);
    let vn = volume_fraction(&tf.nominal, cell_vol);
    let vd = volume_fraction(&tf.dilated, cell_vol);
    RobustReport {
        rho,
        compliance_eroded: trace,
        volumes: (ve, vn, vd),
        erosion_retention: ve / vn.max(1e-30),
    }
}

fn solve(op: &DensityElasticity, b: &[f64]) -> Vec<f64> {
    let mut st = fs_solver::CgState::new(op, &fs_sparse::precond::IdentityPrecond, b);
    let rep = st.run(op, &fs_sparse::precond::IdentityPrecond, 1e-11, 50_000);
    assert!(rep.converged, "elasticity solve failed: {rep:?}");
    st.x
}
