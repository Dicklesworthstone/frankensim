//! Bounded, freshly solved feasible starts for compliance-limited design.
use super::{ComplianceCaps3, ComplianceLimitedOptions3, EnvelopeDesignPoint3,
    CutDensityStudy3, Sdf3Elasticity, EvaluationStop, ProjectedAlError, SolveControl, SolveWork, admit};

#[derive(Debug, Clone, Copy)]
pub struct ComplianceStartOptions3 {
    /// Bisection count after the incoming/full-density endpoint solves. Zero
    /// returns the checked full endpoint. At most 64 bisections are admitted.
    pub max_bisections: usize,
    /// Nonnegative reserve BELOW the smooth constraint boundary, not allowance
    /// above a physical cap. Additional to conservative log-sum-exp smoothing.
    pub margin: f64,
}
impl Default for ComplianceStartOptions3 {
    fn default() -> Self { Self { max_bisections: 16, margin: 1e-3 } }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplianceStartStop3 {
    IncomingSatisfied,
    Restored,
    /// The tested full raw endpoint failed the requested smooth margin.
    /// This is not a proof that every possible design is infeasible.
    FullDensityRejected,
    EvaluationStopped,
}
#[derive(Debug, Clone)]
pub struct ComplianceStartReport3 {
    /// A complete point meeting all physical caps AND the requested smooth
    /// margin. Retained after late interruption; never an unfunded trial.
    pub point: Option<EnvelopeDesignPoint3>,
    /// x = incoming + blend * (1 - incoming), matching point when present.
    pub blend: Option<f64>,
    pub stop: ComplianceStartStop3,
    pub error: Option<EvaluationStop>,
    pub full_density_envelope: Option<f64>,
    pub full_density_maximum_ratio: Option<f64>,
    pub evaluations: usize,
    pub options: ComplianceStartOptions3,
    pub work: SolveWork,
}

/// Restore a numerical feasible start without silently changing any cap,
/// material parameter, physical region or reference load. Each probe solves
/// the actual eroded family; nominal/dilated compliance is never substituted.
///
/// Search the segment toward raw one, retaining only actually checked feasible
/// endpoints. The returned point is NOT a minimum-material result; optimization
/// and its gradient gate are separate. No global monotonicity/infeasibility or
/// sharp root-location claim is needed. Prescribed voids remain void at raw one.
/// The source's nominal state and input controls are unchanged on all outcomes.
///
/// All probes share control. A stopped report may contain a previously solved
/// feasible endpoint, separately labeled from a completed restoration. Policy
/// errors return before any physical solve. Bisection is a seed search only,
/// never finite-difference optimization of the material objective.
#[allow(clippy::too_many_arguments)]
pub fn controlled_compliance_feasible_start<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, forces: &[&[f64]], incoming: &[f64], caps: &ComplianceCaps3,
    design: ComplianceLimitedOptions3, options: ComplianceStartOptions3, control: &mut SolveControl<'_>,
) -> Result<ComplianceStartReport3, ProjectedAlError<EvaluationStop>> {
    admit(study, forces, incoming, caps, design)?;
    if options.max_bisections > 64 || !options.margin.is_finite() || options.margin < 0.0
        || options.margin >= 1.0 - caps.envelope().bias_bound() {
        return Err(ProjectedAlError::Invalid("invalid compliance-start margin or bisection budget"));
    }
    let mut report = ComplianceStartReport3 {
        point: None, blend: None, stop: ComplianceStartStop3::IncomingSatisfied, error: None,
        full_density_envelope: None, full_density_maximum_ratio: None,
        evaluations: 0, options, work: control.work(),
    };
    let mut evaluations = 0;
    let mut evaluate = |rho: Vec<f64>, control: &mut SolveControl<'_>| {
        control.checkpoint("sdf3-compliance-start-probe")?;
        evaluations += 1;
        let evaluation = study.evaluate_robust_envelope(&rho, forces, design.robustness, caps.envelope(), control)?;
        Ok::<_, EvaluationStop>(EnvelopeDesignPoint3 { rho, evaluation })
    };
    let accepts = |p: &EnvelopeDesignPoint3| {
        p.evaluation.envelope.value <= 1.0 - options.margin
            && caps.satisfied_by(&p.evaluation.physical.eroded.objective.case_compliances)
    };
    let attempt = (|| -> Result<(), EvaluationStop> {
        let first = evaluate(incoming.to_vec(), control)?;
        if accepts(&first) { report.point = Some(first); report.blend = Some(0.0); return Ok(()); }
        let full = evaluate(vec![1.0; incoming.len()], control)?;
        report.full_density_envelope = Some(full.evaluation.envelope.value);
        report.full_density_maximum_ratio = Some(full.evaluation.envelope.maximum);
        if !accepts(&full) { report.stop = ComplianceStartStop3::FullDensityRejected; return Ok(()); }
        report.point = Some(full); report.blend = Some(1.0);
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for _ in 0..options.max_bisections {
            control.checkpoint("sdf3-compliance-start-bisect")?;
            let mid = f64::midpoint(lo, hi);
            if mid == lo || mid == hi { break; }
            let rho = incoming.iter().map(|r| (r + mid*(1.0-r)).clamp(*r, 1.0)).collect();
            let trial = evaluate(rho, control)?;
            if accepts(&trial) { hi = mid; report.point = Some(trial); report.blend = Some(mid); }
            else { lo = mid; }
        }
        report.stop = ComplianceStartStop3::Restored;
        Ok(())
    })();
    if let Err(error) = attempt {
        report.stop = ComplianceStartStop3::EvaluationStopped; report.error = Some(error);
    }
    report.evaluations = evaluations;
    report.work = control.work();
    Ok(report)
}

#[cfg(test)]
#[path = "start_tests.rs"]
mod tests;
