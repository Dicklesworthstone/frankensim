//! Bounded enumeration of physical binary endpoints, not a monotone root search.
//! Upper/lower reaction windows need not be monotone in threshold or volume.
use super::*;

/// Endpoint enumeration limits, kept separate from cumulative physics work.
#[derive(Debug, Clone, Copy)]
pub struct BinaryScanOptions3 {
    pub limits: ResponseLimitOptions3,
    pub void_stiffness: f64,
    /// Maximum length of the caller's explicit threshold schedule.
    pub max_candidates: usize,
    /// Candidate's numerical solid volume must not exceed this cap. There is
    /// no automatic assumption that it equals the source's gray material volume.
    pub volume_cap: f64,
    /// Normalized numerical tolerance applied to ALL original limits and G.
    pub feasibility_tolerance: f64,
}
impl Default for BinaryScanOptions3 {
    fn default() -> Self {
        Self { limits: ResponseLimitOptions3::default(), void_stiffness: 1e-3,
            max_candidates: 64, volume_cap: 1.0, feasibility_tolerance: 0.0 }
    }
}
/// One fully solved trial. Retain summaries, not every candidate displacement.
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryCandidate3 {
    pub threshold: f64,
    pub volume_fraction: f64,
    pub constraint: f64,
    pub worst_residual: f64,
    pub response_feasible: bool,
    pub volume_feasible: bool,
}
/// Accepted prefix of an explicit threshold sweep. `best` is the least-volume
/// feasible COMPLETED trial only, not a global binary/topology optimum.
#[derive(Debug, Clone)]
pub struct BinaryScan3 {
    pub candidates: Vec<BinaryCandidate3>,
    pub best: Option<BinaryLimitEvaluation3>,
    /// Index in the original caller schedule. Equal-volume ties keep the first.
    pub best_index: Option<usize>,
    /// Errors after source filtering retain completed evidence and spent work.
    /// An Ok outcome with no best means this schedule found no feasible endpoint.
    pub outcome: Result<(), ResponseError3>,
    pub work: SolveWork,
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Solve every threshold in caller order under ONE shared SolveControl.
    /// Source filtering occurs once, but each candidate prepares its own current
    /// binary stiffness and re-solves every independent experiment. No interpolation
    /// of feasibility, bisection on force windows, cached gray fields, or implicit
    /// threshold/volume repair. Duplicate masks are deliberately still solved.
    ///
    /// Bad schedule/family/design inputs refuse before a reference callback or
    /// physical solve. A later stop retains only fully completed candidates and
    /// their best feasible field. Incoming study scales are never replaced, even
    /// on success. Per-candidate work is cumulative at its solve; report.work also
    /// includes discarded trial work. Memory keeps only the incumbent and current
    /// field families plus scalar trial summaries, not all candidate fields.
    pub fn scan_binary_response_limits(&mut self, rho: &[f64], cases: &[ResponseLimitCase3<'_>],
        thresholds: &[f64], options: BinaryScanOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryScan3, ResponseError3> {
        control.checkpoint("response-binary-scan-start")?;
        if thresholds.is_empty() || thresholds.len() > options.max_candidates
            || !options.volume_cap.is_finite() || !(0.0..=1.0).contains(&options.volume_cap)
            || !options.feasibility_tolerance.is_finite() || options.feasibility_tolerance < 0.0 {
            return Err(ResponseError3::Invalid("invalid binary threshold schedule or feasibility policy"));
        }
        for &threshold in thresholds {
            BinaryEndpointOptions3 { threshold, void_stiffness: options.void_stiffness, limits: options.limits }.validate()?;
        }
        admit_limits(self, cases, options.limits)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(ResponseError3::Invalid("invalid binary-scan source controls"));
        }
        let source = self.design(rho, control)?;
        let mut report = BinaryScan3 { candidates: Vec::new(), best: None,
            best_index: None, outcome: Ok(()), work: control.work() };
        report.outcome = (|| {
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("response-binary-scan-candidate")?;
                let endpoint = self.binary_from_projected(&source.projected, source.volume, cases,
                    BinaryEndpointOptions3 { threshold, void_stiffness: options.void_stiffness, limits: options.limits }, control)?;
                let response_feasible = endpoint.feasible(options.feasibility_tolerance);
                let volume_feasible = endpoint.volume_fraction <= options.volume_cap;
                control.checkpoint("response-binary-scan-retain")?;
                report.candidates.push(BinaryCandidate3 { threshold, volume_fraction: endpoint.volume_fraction,
                    constraint: endpoint.constraint, worst_residual: endpoint.worst_residual,
                    response_feasible, volume_feasible });
                if response_feasible && volume_feasible
                    && report.best.as_ref().is_none_or(|best| endpoint.volume_fraction < best.volume_fraction) {
                    report.best = Some(endpoint);
                    report.best_index = Some(index);
                }
                // Stop after this atomic acceptance retains the completed field.
                control.checkpoint("response-binary-scan-accepted")?;
            }
            control.checkpoint("response-binary-scan-publish")?;
            Ok(())
        })();
        report.work = control.work();
        Ok(report)
    }
}
