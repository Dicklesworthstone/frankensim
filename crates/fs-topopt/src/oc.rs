//! The optimality-criteria update — the classical driver for
//! compliance + volume (documented choice: OC is THE standard for
//! this problem class; fs-ascent's augmented Lagrangian is the
//! general constrained path). Deterministic: fixed bisection on the
//! volume multiplier, fixed move limits, fixed iteration schedule —
//! a whole run replays bitwise.

use crate::elasticity::DensityElasticity;
use crate::pipeline::DesignPipeline;

/// Outcome of an OC run.
#[derive(Debug, Clone)]
pub struct OcReport {
    /// Final raw design.
    pub rho: Vec<f64>,
    /// Compliance trace per iteration.
    pub compliance: Vec<f64>,
    /// Volume-fraction trace (of the PROJECTED design).
    pub volume: Vec<f64>,
    /// Iterations run.
    pub iters: usize,
    /// Max design change in the final iteration.
    pub final_change: f64,
}

/// Volume fraction of a projected design under cell volumes.
fn volume_fraction(pipeline: &DesignPipeline, rho: &[f64], cell_vol: &[f64]) -> f64 {
    let (_, rho_bar, _) = pipeline.forward(rho);
    let total: f64 = cell_vol.iter().sum();
    rho_bar
        .iter()
        .zip(cell_vol)
        .map(|(r, v)| r * v)
        .sum::<f64>()
        / total
}

/// Validate the common nominal/robust OC modeling boundary before iteration.
pub(crate) fn assert_valid_oc_inputs(
    elasticity: &DensityElasticity,
    force: &[f64],
    rho0: &[f64],
    cell_vol: &[f64],
    vol_frac: f64,
    move_limit: f64,
) {
    assert_eq!(
        rho0.len(),
        elasticity.cells(),
        "initial design must contain one density per elasticity cell"
    );
    assert_eq!(
        cell_vol.len(),
        rho0.len(),
        "OC requires one cell volume per design density"
    );
    assert_eq!(
        force.len(),
        elasticity.n(),
        "force must contain one value per elasticity dof"
    );
    assert!(
        rho0.iter()
            .all(|value| value.is_finite() && (0.0..=1.0).contains(value)),
        "initial design densities must be finite and lie in [0, 1]"
    );
    assert!(
        force.iter().all(|value| value.is_finite()),
        "force must contain only finite values"
    );
    assert!(
        cell_vol
            .iter()
            .all(|volume| volume.is_finite() && *volume > 0.0),
        "cell volumes must be finite and positive"
    );
    assert!(
        cell_vol.iter().sum::<f64>().is_finite(),
        "total cell volume must be representable as a finite value"
    );
    assert!(
        vol_frac.is_finite() && vol_frac > 0.0 && vol_frac <= 1.0,
        "target volume fraction must be finite and lie in (0, 1]"
    );
    assert!(
        move_limit.is_finite() && move_limit > 0.0 && move_limit <= 1.0,
        "OC move limit must be finite and lie in (0, 1]"
    );
}

/// Run OC iterations at FIXED continuation parameters (drivers wrap
/// this with β/p schedules). Volume constraint is on the projected
/// design; the multiplier is found by deterministic bisection.
#[allow(clippy::too_many_arguments)]
pub fn optimality_criteria(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    force: &[f64],
    rho0: &[f64],
    cell_vol: &[f64],
    vol_frac: f64,
    move_limit: f64,
    iters: usize,
) -> OcReport {
    assert_valid_oc_inputs(elasticity, force, rho0, cell_vol, vol_frac, move_limit);
    let nc = rho0.len();
    let mut rho = rho0.to_vec();
    let mut compliance_trace = Vec::with_capacity(iters);
    let mut volume_trace = Vec::with_capacity(iters);
    let mut final_change = 0.0f64;
    for _ in 0..iters {
        let (c, _u, grad) = pipeline.compliance_and_gradient(elasticity, &rho, force);
        compliance_trace.push(c);
        // OC multiplicative update with move limits: the compliance
        // gradient is ≤ 0 (energies), so −grad ≥ 0 drives growth.
        let sensitivity: Vec<f64> = grad.iter().map(|g| (-g).max(1e-30)).collect();
        let mut lo = 1e-12f64;
        let mut hi = 1e12f64;
        let mut candidate = rho.clone();
        for _ in 0..80 {
            let lambda = fs_math::det::sqrt(lo * hi);
            for i in 0..nc {
                let scale = fs_math::det::sqrt(sensitivity[i] / (lambda * cell_vol[i]));
                let stepped = rho[i] * scale;
                candidate[i] = stepped
                    .clamp(rho[i] - move_limit, rho[i] + move_limit)
                    .clamp(1e-3, 1.0);
            }
            if volume_fraction(pipeline, &candidate, cell_vol) > vol_frac {
                lo = lambda;
            } else {
                hi = lambda;
            }
        }
        final_change = rho
            .iter()
            .zip(&candidate)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        rho = candidate;
        volume_trace.push(volume_fraction(pipeline, &rho, cell_vol));
    }
    OcReport {
        rho,
        compliance: compliance_trace,
        volume: volume_trace,
        iters,
        final_change,
    }
}
