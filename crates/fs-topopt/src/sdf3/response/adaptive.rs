//! Fit -> selectively rebuild -> refit using immutable reference experiments.
//! Only raw densities transfer. Loads and observations are integrated on the
//! candidate's own quadrature, motion liftings are rebuilt by the original
//! evaluator, and a new optimizer starts at a freshly solved baseline. A coarse
//! optimizer's multipliers/gradients are never relabeled as fine-grid state.
use super::*;
use super::refinement::ReferenceResponseCase3;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;

type FitError = ProjectedAlError<ResponseError3>;

/// Complete accepted fitting prefix on ONE grid. Initializing the grid can
/// fail before this exists. Once initialization succeeds, an interrupted run
/// returns its error in `outcome` alongside the last accepted physical state.
#[derive(Debug)]
pub struct ReferenceResponseFit3 {
    /// New-grid baseline, not prolonged coarse displacements or old loss values.
    pub initial: ResponseEvaluation3,
    pub accepted: ResponseEvaluation3,
    pub history: Vec<ProjectedResponseIteration3>,
    /// Iteration/evaluation limits and stalls are not convergence. A successful
    /// optimizer return does not itself assert nonlinear volume feasibility.
    pub outcome: Result<ProjectedAlReport, FitError>,
    /// Cumulative physical work, including rejected/interrupted candidate work.
    pub work: SolveWork,
}

/// Unpublished refined candidate. The source study is borrowed immutably and
/// remains untouched even when initialization, restoration or fitting fails.
/// Inspect the fit's stop and feasibility BEFORE replacing the source study.
pub struct ResponseRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub fit: ReferenceResponseFit3,
    /// Exactly parent-inherited raw values, before new-grid volume restoration.
    pub inherited_rho: Vec<f64>,
    /// Actual feasible starting raw values for the new-grid baseline solve.
    pub starting_rho: Vec<f64>,
}

fn physics(error: impl Into<ResponseError3>) -> FitError {
    ProjectedAlError::Evaluation(error.into())
}

// Admission before any reference callback, dense nodal observation allocation,
// or physical solve. Reuse the existing optimizer's own dimension/work policy.
fn admit_reference<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>,
    cases: &[ReferenceResponseCase3<'_>], options: ProjectedResponseOptions3) -> Result<(), FitError> {
    options.optimizer.validate(study.cells())?;
    if !options.volume_cap.is_finite() || options.volume_cap <= 0.0 || options.volume_cap > 1.0
        || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
        || !options.objective_scale.is_finite() || options.objective_scale <= 0.0
        || !options.response.volume_weight.is_finite() || options.response.volume_weight < 0.0
        || cases.is_empty() || cases.len() > options.response.max_cases {
        return Err(ProjectedAlError::Invalid("invalid reference response family or fitting policy"));
    }
    let mut count = 0usize;
    let mut positive = options.response.volume_weight > 0.0;
    for case in cases {
        count = count.checked_add(case.targets.len()).filter(|n| *n <= options.response.max_observations)
            .ok_or(ProjectedAlError::Invalid("reference response observation allowance"))?;
        if case.targets.is_empty() || (case.prescribed.is_some()
            && study.operator.adaptive().embedded_dirichlet_penalty().is_none()) {
            return Err(ProjectedAlError::Invalid("empty targets or absent embedded motion support"));
        }
        for target in case.targets {
            if !target.target.is_finite() || !target.scale.is_finite() || target.scale <= 0.0
                || !target.weight.is_finite() || target.weight < 0.0 {
                return Err(ProjectedAlError::Invalid("invalid reference response target"));
            }
            positive |= target.weight > 0.0;
        }
    }
    if !positive { return Err(ProjectedAlError::Invalid("zero reference response objective")); }
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Bind reference force/observation laws to this grid and run the EXISTING
    /// projected response optimizer. The same pure laws can then be used on a
    /// refined grid; nodal q/f arrays are never interpolated or nearest-matched.
    /// Each law is integrated once for this fit, not on every optimizer trial.
    /// Density-dependent prescribed motion remains in evaluate_responses.
    ///
    /// A fresh optimizer is intentionally created here. For continuation on the
    /// SAME grid with retained multiplier/spectral state, use
    /// ProjectedResponseStudy3 directly. A failure before its first accepted
    /// sample returns Err without changing this study's material. Later errors
    /// remain in `outcome` with an aligned accepted state and spent work.
    pub fn fit_reference_responses(&mut self, cases: &[ReferenceResponseCase3<'_>], rho: &[f64],
        options: ProjectedResponseOptions3, steps: usize, control: &mut SolveControl<'_>)
        -> Result<ReferenceResponseFit3, FitError> {
        control.checkpoint("response-reference-start").map_err(physics)?;
        admit_reference(self, cases, options)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("invalid reference fitting start"));
        }
        let mut forces = Vec::with_capacity(cases.len());
        let mut observations = Vec::with_capacity(cases.len());
        for case in cases {
            control.checkpoint("response-reference-assemble").map_err(physics)?;
            let op = self.operator.adaptive();
            forces.push(op.reference_load(case.load, || poll(control)).map_err(physics)?);
            let mut rows = Vec::with_capacity(case.targets.len());
            for target in case.targets {
                control.checkpoint("response-reference-assemble").map_err(physics)?;
                rows.push(op.reference_load(target.observation, || poll(control)).map_err(physics)?);
            }
            observations.push(rows);
        }
        // These borrowed views live only for the synchronous fit. Reports own
        // their fields, so no self-referential experiment or dangling view exists.
        let targets: Vec<Vec<ResponseTarget3<'_>>> = cases.iter().zip(&observations).map(|(case, rows)| {
            case.targets.iter().zip(rows).map(|(target, q)| ResponseTarget3 {
                q, target: target.target, scale: target.scale, weight: target.weight,
            }).collect()
        }).collect();
        let nodal: Vec<ResponseCase3<'_>> = cases.iter().zip(&forces).zip(&targets)
            .map(|((case, force), targets)| ResponseCase3 { force, prescribed: case.prescribed, targets }).collect();
        let mut session = ProjectedResponseStudy3::new(self, &nodal, rho, options, control)?;
        let initial = session.accepted().clone();
        let outcome = session.run(steps);
        Ok(ReferenceResponseFit3 {
            initial, accepted: session.accepted().clone(), history: session.history().to_vec(),
            outcome, work: session.work(),
        })
    }

    /// Transfer an accepted raw design into an explicitly constructed refined
    /// pipeline, restore its actual projected volume, and fit fresh reference
    /// experiments. The candidate is consumed, never installed into `self`.
    /// No old primal/adjoint field, gradient, multiplier or spectral step is used
    /// by its optimizer. The original study remains usable after EVERY result.
    ///
    /// The existing AdaptiveTransfer3 checks nested active support, geometry
    /// bounds, material/support method and clamp compatibility. Equality of the
    /// SDF and patch callbacks remains a caller obligation. Pipeline parameters
    /// and filter radius are explicit in `candidate`; changing them deliberately
    /// produces a fresh baseline, never a cross-grid descent claim.
    ///
    /// Only complete raw-design provenance is required from `accepted`; its
    /// loss and displacement caches do NOT enter the fine fit. The same pure
    /// experiment laws/targets should be supplied to both fits and the estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_reference_responses<F: AdaptiveSdf3Elasticity>(&self, mut candidate: CutDensityStudy3<F>,
        accepted: &ResponseEvaluation3, cases: &[ReferenceResponseCase3<'_>], options: ProjectedResponseOptions3,
        steps: usize, max_transfer_terms: usize, control: &mut SolveControl<'_>)
        -> Result<ResponseRefit3<F>, FitError> {
        control.checkpoint("response-refit-start").map_err(physics)?;
        admit_reference(&candidate, cases, options)?;
        if accepted.rho.len() != self.cells() || accepted.rho.iter()
            .any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("source raw design violates refit shape or bounds"));
        }
        let source_design = self.design(&accepted.rho, control).map_err(physics)?;
        if source_design.scales != accepted.scales || source_design.projected != accepted.projected_rho {
            return Err(ProjectedAlError::Invalid("source evaluation belongs to another density pipeline"));
        }
        let inherited_rho = {
            let coarse = self.operator.adaptive();
            let fine = candidate.operator.adaptive();
            let transfer = AdaptiveTransfer3::new(coarse, fine, max_transfer_terms, || poll(control)).map_err(physics)?;
            if !fine.leaves().iter().zip(transfer.parents()).any(|(leaf, &parent)| leaf.level() > coarse.leaves()[parent].level()) {
                return Err(ProjectedAlError::Invalid("refit requires an enriched active space"));
            }
            let mut raw = Vec::with_capacity(candidate.cells());
            for &parent in transfer.parents() {
                control.checkpoint("response-refit-transfer").map_err(physics)?;
                raw.push(accepted.rho[parent]);
            }
            raw
        };
        // Use the declared optimizer floor, not the unrelated fixed OC floor.
        // Keep an actually evaluated feasible endpoint, with no cap tolerance
        // silently consuming the optimizer's later feasibility allowance.
        let starting_rho = candidate.feasible_start_with_floor(&inherited_rho, options.volume_cap,
            0.0, options.density_floor, control).map_err(physics)?;
        control.checkpoint("response-refit-baseline").map_err(physics)?;
        let fit = candidate.fit_reference_responses(cases, &starting_rho, options, steps, control)?;
        Ok(ResponseRefit3 { study: candidate, fit, inherited_rho, starting_rho })
    }
}
