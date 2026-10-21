//! Checked feasible seeds for the entire combined-load set, with fixed caps.
use super::*;
pub use crate::sdf3::robust::envelope::limits::start::{
    ComplianceStartOptions3 as LoadSpaceStartOptions3, ComplianceStartStop3 as LoadSpaceStartStop3,
};

#[derive(Debug, Clone)]
pub struct LoadSpaceStartReport3 {
    /// Fully solved, cap-feasible point; retained even after a late stop.
    /// This is a restoration seed, not an optimized minimum-material design.
    pub point: Option<LoadSpaceDesignPoint3>,
    /// raw = incoming + blend*(1-incoming), matching the retained point.
    pub blend: Option<f64>,
    pub stop: LoadSpaceStartStop3,
    pub error: Option<EvaluationStop>,
    pub full_density_smooth_maximum: Option<f64>,
    pub full_density_worst_compliance: Option<f64>,
    pub evaluations: usize,
    pub options: LoadSpaceStartOptions3,
    pub work: SolveWork,
}

/// Bounded search toward raw one, re-solving every basis and cross-compliance
/// spectrum at each probe. Never freeze a previously selected worst direction,
/// infer a cap from a baseline, loosen a requirement or change protected voids.
///
/// Both the physical worst-work cap and conservative smooth constraint must
/// pass, with options.margin reserved BELOW the smooth boundary. The margin is
/// not permission to exceed the physical cap. A refused full endpoint is not a
/// global infeasibility proof. Every returned point is actually checked; no
/// monotonicity or root-location certificate is claimed by this seed search.
///
/// All probes share the original SolveControl. Returned failures preserve the
/// incoming nominal model/scales and any previously checked feasible endpoint.
/// Invalid policies return before the first solve; at most 2+max_bisections
/// complete load-space evaluations are attempted. Optimization is separate.
#[allow(clippy::too_many_arguments)]
pub fn controlled_load_space_feasible_start<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, basis: &[&[f64]], incoming: &[f64], cap: &LoadSpaceCap3,
    design: LoadSpaceLimitedOptions3, options: LoadSpaceStartOptions3, control: &mut SolveControl<'_>,
) -> Result<LoadSpaceStartReport3, ProjectedAlError<EvaluationStop>> {
    super::super::admit(study, basis, incoming, cap.policy(), engine_options(design))?;
    if options.max_bisections > 64 || !options.margin.is_finite() || options.margin < 0.0
        || options.margin >= 1.0-cap.policy().bias_bound() {
        return Err(ProjectedAlError::Invalid("invalid load-space cap restoration policy"));
    }
    let mut report = LoadSpaceStartReport3 { point: None, blend: None,
        stop: LoadSpaceStartStop3::IncomingSatisfied, error: None,
        full_density_smooth_maximum: None, full_density_worst_compliance: None,
        evaluations: 0, options, work: control.work() };
    let mut evaluations = 0;
    let mut evaluate = |rho: Vec<f64>, control: &mut SolveControl<'_>| {
        control.checkpoint("sdf3-load-space-cap-start-probe")?;
        evaluations += 1;
        let evaluation = study.evaluate_robust_load_space(&rho, basis, design.robustness, cap.policy(), control)?;
        Ok::<_, EvaluationStop>(LoadSpaceDesignPoint3 { rho, evaluation })
    };
    let accepts = |p: &LoadSpaceDesignPoint3| cap.admits(&p.evaluation, 0.0)
        && p.evaluation.spectrum.smooth_maximum <= 1.0-options.margin;
    let attempt = (|| -> Result<(), EvaluationStop> {
        let first = evaluate(incoming.to_vec(), control)?;
        if accepts(&first) { report.point = Some(first); report.blend = Some(0.0); return Ok(()); }
        let full = evaluate(vec![1.0; incoming.len()], control)?;
        report.full_density_smooth_maximum = Some(full.evaluation.spectrum.smooth_maximum);
        report.full_density_worst_compliance = Some(full.evaluation.worst_compliance);
        if !accepts(&full) { report.stop = LoadSpaceStartStop3::FullDensityRejected; return Ok(()); }
        report.point = Some(full); report.blend = Some(1.0);
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for _ in 0..options.max_bisections {
            control.checkpoint("sdf3-load-space-cap-start-bisect")?;
            let mid = f64::midpoint(lo, hi);
            if mid == lo || mid == hi { break; }
            let rho = incoming.iter().map(|r| (r+mid*(1.0-r)).clamp(*r, 1.0)).collect();
            let trial = evaluate(rho, control)?;
            if accepts(&trial) { hi = mid; report.point = Some(trial); report.blend = Some(mid); }
            else { lo = mid; }
        }
        report.stop = LoadSpaceStartStop3::Restored;
        Ok(())
    })();
    if let Err(error) = attempt { report.stop = LoadSpaceStartStop3::EvaluationStopped; report.error = Some(error); }
    report.evaluations = evaluations;
    report.work = control.work();
    Ok(report)
}

#[cfg(test)]
#[path = "start_tests.rs"]
mod tests;
