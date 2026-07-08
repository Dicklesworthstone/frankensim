//! Stress-constrained SIMP: relaxed per-cell von Mises measures with
//! qp-RELAXATION (stress exponent q < penal, so void-cell stresses
//! vanish in the MEASURE while E_min keeps u finite — the singularity
//! trap's standard treatment, gated by regression), p-norm
//! aggregation with ADAPTIVE normalization (c_k = σ_max/PN tracked
//! per iteration so the aggregate follows the true maximum), and the
//! NON-self-adjoint design gradient: one extra adjoint solve
//! K·λ = (∂PN/∂u)ᵀ, then explicit + implicit chain through the
//! landed filter/projection pipeline.

use crate::elasticity::DensityElasticity;
use crate::pipeline::DesignPipeline;

/// Von Mises scalar from a Voigt stress vector.
#[must_use]
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let sq = (s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2);
    let shear = s[3].powi(2) + s[4].powi(2) + s[5].powi(2);
    fs_math::det::sqrt(3.0f64.mul_add(shear, 0.5 * sq))
}

/// d(von Mises)/dσ (Voigt 6-vector).
#[must_use]
pub fn von_mises_derivative(s: &[f64; 6]) -> [f64; 6] {
    let vm = von_mises(s).max(1e-12);
    [
        (2.0 * s[0] - s[1] - s[2]) / (2.0 * vm) * 1.0,
        (2.0 * s[1] - s[0] - s[2]) / (2.0 * vm) * 1.0,
        (2.0 * s[2] - s[0] - s[1]) / (2.0 * vm) * 1.0,
        3.0 * s[3] / vm,
        3.0 * s[4] / vm,
        3.0 * s[5] / vm,
    ]
}

/// Aggregated-stress report.
#[derive(Debug, Clone)]
pub struct StressReport {
    /// The p-norm aggregate PN.
    pub pnorm: f64,
    /// True maximum relaxed stress.
    pub sigma_max: f64,
    /// The adaptive normalization c = σ_max/PN (multiply PN by this
    /// to track the max — the honesty knob, reported per evaluation).
    pub adaptive_c: f64,
    /// Design gradient of PN w.r.t. the RAW design.
    pub gradient: Vec<f64>,
    /// Per-cell relaxed measures σ̃ (the aggregate's inputs — lets
    /// callers audit WHICH cells drive the constraint).
    pub relaxed: Vec<f64>,
}

/// Evaluate the relaxed p-norm stress aggregate and its EXACT design
/// gradient. `q_relax` < penal is the qp-relaxation exponent
/// (typical 0.5); `p_agg` the aggregation power (typical 8).
pub fn stress_aggregate(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    rho: &[f64],
    force: &[f64],
    q_relax: f64,
    p_agg: f64,
) -> StressReport {
    let p = &pipeline.params;
    let (rho_tilde, rho_bar, moduli) = pipeline.forward(rho);
    elasticity.moduli = moduli;
    let u = solve(elasticity, force);
    let sig = elasticity.cell_stress(&u);
    let nc = sig.len();
    // Relaxed measures σ̃_c = ρ̄^q·E(ρ̄)·vm(σ_unit) — the PHYSICAL
    // stress carries the modulus (σ = E·CB·u); relaxation multiplies
    // by ρ̄^{q−p}-ish. We use the direct form: σ̃ = ρ̄^q·vm(E·σ_unit)
    // = ρ̄^q·E·vm(σ_unit).
    let vm: Vec<f64> = sig.iter().map(von_mises).collect();
    let e_of = |rb: f64| -> f64 {
        let rc = rb.clamp(0.0, 1.0);
        p.e_min + (1.0 - p.e_min) * fs_math::det::pow(rc.max(1e-12), p.penal)
    };
    let relax = |rb: f64| -> f64 {
        fs_math::det::pow(rb.clamp(0.0, 1.0).max(1e-12), q_relax)
    };
    let tilde: Vec<f64> = rho_bar
        .iter()
        .zip(&vm)
        .map(|(&rb, &v)| relax(rb) * e_of(rb) * v)
        .collect();
    let sigma_max = tilde.iter().copied().fold(0.0f64, f64::max);
    // PN = (1/N · Σ σ̃^P)^{1/P}.
    let nfin = nc as f64;
    let sum_p: f64 = tilde
        .iter()
        .map(|t| fs_math::det::pow(t.max(1e-30), p_agg))
        .sum::<f64>()
        / nfin;
    let pnorm = fs_math::det::pow(sum_p, 1.0 / p_agg);
    // ∂PN/∂σ̃_c = (1/N)·σ̃^{P−1}·PN^{1−P}.
    let dpn_dtilde: Vec<f64> = tilde
        .iter()
        .map(|t| {
            fs_math::det::pow(t.max(1e-30), p_agg - 1.0)
                * fs_math::det::pow(pnorm, 1.0 - p_agg)
                / nfin
        })
        .collect();
    // IMPLICIT part: ∂PN/∂u = Σ_c ∂PN/∂σ̃_c·ρ̄^q·E·d(vm)/dσ·CB_c.
    let dj_dsigma: Vec<[f64; 6]> = sig
        .iter()
        .zip(rho_bar.iter().zip(&dpn_dtilde))
        .map(|(s, (&rb, &dp))| {
            let dvm = von_mises_derivative(s);
            let scale = dp * relax(rb) * e_of(rb);
            [
                scale * dvm[0],
                scale * dvm[1],
                scale * dvm[2],
                scale * dvm[3],
                scale * dvm[4],
                scale * dvm[5],
            ]
        })
        .collect();
    let rhs = elasticity.stress_pullback(&dj_dsigma);
    let lambda = solve(elasticity, &rhs);
    // Implicit design term: −λᵀ·E′·K_c·u per cell (K's density chain).
    let mut cross = vec![0.0f64; nc];
    // uᵀK_cλ via the polarization identity on cell_energies:
    // e(u+λ) − e(u) − e(λ) = 2·uᵀK_cλ.
    let sum_vec: Vec<f64> = u.iter().zip(&lambda).map(|(a, c)| a + c).collect();
    let e_sum = elasticity.cell_energies(&sum_vec);
    let e_u = elasticity.cell_energies(&u);
    let e_l = elasticity.cell_energies(&lambda);
    for c in 0..nc {
        cross[c] = 0.5 * (e_sum[c] - e_u[c] - e_l[c]);
    }
    // Total dPN/dρ̄: explicit (relaxation + modulus factors on σ̃) −
    // implicit (through u).
    let dlam_drhobar: Vec<f64> = rho_bar
        .iter()
        .zip(vm.iter().zip(dpn_dtilde.iter().zip(&cross)))
        .map(|(&rb, (&v, (&dp, &cr)))| {
            let rc = rb.clamp(0.0, 1.0).max(1e-12);
            let dsimp = (1.0 - p.e_min) * p.penal * fs_math::det::pow(rc, p.penal - 1.0);
            let drelax = q_relax * fs_math::det::pow(rc, q_relax - 1.0);
            let explicit = dp * v * drelax.mul_add(e_of(rb), relax(rb) * dsimp);
            let implicit = dsimp * cr;
            explicit - implicit
        })
        .collect();
    let chained: Vec<f64> = rho_tilde
        .iter()
        .zip(&dlam_drhobar)
        .map(|(&rt, &d)| d * crate::filter::heaviside_derivative(rt, p.beta, p.eta))
        .collect();
    let gradient = pipeline.filter.apply_transpose(&chained);
    StressReport {
        pnorm,
        sigma_max,
        adaptive_c: sigma_max / pnorm.max(1e-30),
        gradient,
        relaxed: tilde,
    }
}

fn solve(op: &DensityElasticity, b: &[f64]) -> Vec<f64> {
    let mut st = fs_solver::CgState::new(op, &fs_sparse::precond::IdentityPrecond, b);
    let rep = st.run(op, &fs_sparse::precond::IdentityPrecond, 1e-11, 50_000);
    assert!(rep.converged, "stress solve failed: {rep:?}");
    st.x
}
