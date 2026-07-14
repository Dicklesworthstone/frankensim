//! Stage 4: ANYTIME-VALID fragility — P(peak drift ratio > limit)
//! over a Kanai–Tajimi ensemble, estimated with an fs-eproc
//! confidence sequence that is valid AT the data-dependent stopping
//! time by construction: the study stops itself once the confidence-sequence
//! radius is below the requested margin. An fs-uq MLMC report over dt-refinement
//! levels rides along for the level-design evidence. Exceedance
//! indicators are ½-sub-Gaussian, so σ = ½ is a hard bound, not a
//! plug-in estimate.

use crate::assert_ground_motion_ensemble;
use crate::history::{StoryFrame, StoryParams, peak_drift};
use fs_eproc::GaussianMixtureCs;
use fs_scenario::ensemble::StochasticEnsemble;
use fs_uq::mlmc::{MlmcReport, mlmc_estimate};

/// The fragility record.
pub struct FragilityReport {
    /// Members actually consumed (the e-stop decides).
    pub members_used: u32,
    /// Exceedance estimate (CS center).
    pub p_hat: f64,
    /// CS radius at stop.
    pub radius: f64,
    /// Did the study stop before exhausting the budget?
    pub stopped_early: bool,
    /// Exceedances observed.
    pub exceedances: u32,
    /// The MLMC dt-level report on the peak-drift QoI.
    pub mlmc: MlmcReport,
}

/// Run the e-stopped fragility study: consume ensemble members one at
/// a time, feed exceedance indicators to the confidence sequence, and
/// stop when decision-grade (`radius ≤ margin`). `alpha` is the anytime
/// validity level. This API has no decision threshold distinct from the
/// physical `drift_limit`, so it makes no one-sided-decision stopping claim.
///
/// # Panics
/// On ensemble realization errors (spec defects — programmer
/// contracts at smoke scale).
#[must_use]
pub fn e_stopped_fragility(
    ensemble: &StochasticEnsemble,
    params: StoryParams,
    drift_limit: f64,
    alpha: f64,
    margin: f64,
) -> FragilityReport {
    assert_ground_motion_ensemble(ensemble);
    let dt = ensemble.dt.value;
    let mut cs = GaussianMixtureCs::new(0.5, 8.0, alpha);
    let mut exceedances = 0u32;
    let mut used = 0u32;
    let mut stopped_early = false;
    for member in 0..ensemble.members {
        let real = ensemble.realize(member).expect("ensemble realizes");
        let mut frame = StoryFrame::new(params);
        let drifts = frame.run(&real.values, dt);
        let pd = peak_drift(&drifts, params.h);
        let x = if pd > drift_limit { 1.0 } else { 0.0 };
        if x > 0.5 {
            exceedances += 1;
        }
        cs.observe(x);
        used = member + 1;
        if let Some((_, radius)) = cs.interval() {
            // Decision-grade: the interval is tight enough, and we've
            // seen enough members for the asymptotics to mean anything.
            if used >= 8 && radius <= margin {
                stopped_early = used < ensemble.members;
                break;
            }
        }
    }
    let (p_hat, radius) = cs.interval().expect("observed at least one member");
    // MLMC over dt levels (coarse 2dt vs dt) on the CONTINUOUS QoI
    // (peak drift) — the level-design evidence for the full tier.
    let mut sampler = |level: usize, g: u64| -> f64 {
        let member = u32::try_from(g % u64::from(ensemble.members)).expect("small");
        let real = ensemble.realize(member).expect("ensemble realizes");
        let run_at = |factor: usize| -> f64 {
            let mut frame = StoryFrame::new(params);
            let coarse: Vec<f64> = real.values.iter().copied().step_by(factor).collect();
            let drifts = frame.run(&coarse, dt * factor as f64);
            peak_drift(&drifts, params.h)
        };
        if level == 0 {
            run_at(2)
        } else {
            run_at(1) - run_at(2)
        }
    };
    let mlmc = mlmc_estimate(&mut sampler, &[1.0, 3.0], 6, 1e-6);
    FragilityReport {
        members_used: used,
        p_hat,
        radius,
        stopped_early,
        exceedances,
        mlmc,
    }
}
