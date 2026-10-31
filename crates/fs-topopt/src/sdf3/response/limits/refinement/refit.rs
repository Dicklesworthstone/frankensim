//! Bind unchanged reference inequalities to a new grid and restart the existing
//! minimum-material optimizer. Unlike a volume-capped compliance design there
//! is NO general monotone density scaling that restores force-window feasibility.
use super::*;
use crate::sdf3::regions::inherit_physical_regions;
use super::super::design::{MinimumMaterialStudy3, MinimumMaterialOptions3, MinimumMaterialIteration3};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport};

type FitError = ProjectedAlError<ResponseError3>;
fn physics(e: impl Into<ResponseError3>) -> FitError { ProjectedAlError::Evaluation(e.into()) }

/// ONE grid's complete fitting prefix. Initial field failure returns Err before
/// this exists; later failures are retained beside the last accepted state.
#[derive(Debug)]
pub struct ReferenceMinimumMaterialFit3 {
    pub initial: ResponseLimitEvaluation3,
    pub accepted: ResponseLimitEvaluation3,
    /// Least-objective accepted point satisfying BOTH smooth and original limits.
    /// None is a genuine absence of a feasible result, not an empty placeholder.
    pub best_feasible: Option<ResponseLimitEvaluation3>,
    pub history: Vec<MinimumMaterialIteration3>,
    /// Stop/KKT belong to `accepted`, not necessarily `best_feasible`.
    pub outcome: Result<ProjectedAlReport, FitError>,
    pub work: SolveWork,
}
/// Candidate is never installed into its immutably borrowed source. Its scales
/// match `fit.accepted`; consuming `into_best_feasible` selects an earlier
/// feasible point explicitly, without transferring its predecessor's KKT claim.
pub struct MinimumMaterialRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub fit: ReferenceMinimumMaterialFit3,
    pub inherited_rho: Vec<f64>,
}
impl<O: AdaptiveSdf3Elasticity> MinimumMaterialRefit3<O> {
    pub fn into_best_feasible(mut self, control: &mut SolveControl<'_>)
        -> Result<(CutDensityStudy3<O>, ResponseLimitEvaluation3), ResponseError3> {
        control.checkpoint("minimum-material-refit-select")?;
        let mut best = self.fit.best_feasible.take()
            .ok_or(ResponseError3::Invalid("refitted grid has no accepted response-limit-feasible design"))?;
        self.study.operator.set_scales(&best.scales)?;
        best.work = control.work();
        Ok((self.study, best))
    }
}
fn admit<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>, cases: &[ReferenceResponseLimitCase3<'_>],
    rho: &[f64], options: MinimumMaterialOptions3) -> Result<(), FitError> {
    options.optimizer.validate(study.cells())?;
    if rho.len() != study.cells() || !options.density_floor.is_finite() || options.density_floor <= 0.0
        || options.density_floor >= 1.0 || !options.volume_scale.is_finite() || options.volume_scale <= 0.0
        || !(1.0/options.volume_scale).is_finite()
        || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
        return Err(ProjectedAlError::Invalid("invalid reference minimum-material starting design or policy"));
    }
    super::super::design::admit_objective(study, options)?;
    admit_reference_limits(study, cases, options.limits).map_err(physics)
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Integrate external and displacement-observation laws once on this grid,
    /// then use the ORIGINAL minimum-material session. Reaction modes remain
    /// current-density observations during every trial, not frozen nodal q's.
    /// The whole family and options are admitted before reference callbacks.
    /// A later optimization error retains fields/work rather than refunding or
    /// substituting a coarse solution. Nonlinear infeasible starts are allowed.
    pub fn fit_reference_minimum_material(&mut self, cases: &[ReferenceResponseLimitCase3<'_>],
        rho: &[f64], options: MinimumMaterialOptions3, steps: usize, control: &mut SolveControl<'_>)
        -> Result<ReferenceMinimumMaterialFit3, FitError> {
        control.checkpoint("minimum-material-reference-start").map_err(physics)?;
        admit(self, cases, rho, options)?;
        let mut forces = Vec::with_capacity(cases.len());
        let mut observations = Vec::with_capacity(cases.len());
        for case in cases {
            control.checkpoint("minimum-material-reference-assemble").map_err(physics)?;
            let op = self.operator.adaptive();
            forces.push(op.reference_load(case.load, || poll(control)).map_err(physics)?);
            let mut row = Vec::with_capacity(case.limits.len());
            for limit in case.limits {
                row.push(match limit.response {
                    ReferenceBoundedResponse3::Displacement(law) =>
                        Some(op.reference_load(law, || poll(control)).map_err(physics)?),
                    ReferenceBoundedResponse3::Reaction(_) => None,
                });
            }
            observations.push(row);
        }
        let limits: Vec<Vec<ResponseLimit3<'_>>> = cases.iter().zip(&observations).map(|(case, row)| {
            case.limits.iter().zip(row).map(|(limit, q)| ResponseLimit3 {
                sense: limit.sense, bound: limit.bound, scale: limit.scale,
                response: match limit.response {
                    ReferenceBoundedResponse3::Displacement(_) => BoundedResponse3::Displacement(q.as_ref().expect("assembled observation")),
                    ReferenceBoundedResponse3::Reaction(mode) => BoundedResponse3::Reaction(mode),
                },
            }).collect()
        }).collect();
        let nodal: Vec<_> = cases.iter().zip(&forces).zip(&limits).map(|((case, force), limits)|
            ResponseLimitCase3 { force, prescribed: case.prescribed, limits }).collect();
        let mut session = MinimumMaterialStudy3::new(self, &nodal, rho, options, control)?;
        let initial = session.accepted().clone();
        let outcome = session.run(steps);
        Ok(ReferenceMinimumMaterialFit3 { initial, accepted: session.accepted().clone(),
            best_feasible: session.best_feasible().cloned(), history: session.history().to_vec(),
            outcome, work: session.work() })
    }

    /// Refit raw controls on an explicitly constructed refined pipeline. No
    /// coarse displacement, aggregate adjoint, gradient, multiplier or spectral
    /// step is transferred. The first fine response/feasibility is a new solve.
    /// Prescribed physical regions inherit their exact checked parent labels.
    ///
    /// A coarse feasible design may start INFEASIBLE here: changing density
    /// toward zero/one does not generally restore upper/lower force windows.
    /// The original AL handles that start; best_feasible is absent until an
    /// actual accepted fine point passes every declared numerical limit. No
    /// volume restoration or invented response-feasibility certificate occurs.
    ///
    /// Candidate filter/radius/SIMP choices are explicit. Geometry, reference
    /// material/support laws and response limits must be held consistent by the
    /// caller. Source fields provide raw-design provenance, never a fine result.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_reference_minimum_material<F: AdaptiveSdf3Elasticity>(&self, mut candidate: CutDensityStudy3<F>,
        source: &ResponseLimitEvaluation3, cases: &[ReferenceResponseLimitCase3<'_>], options: MinimumMaterialOptions3,
        steps: usize, max_transfer_terms: usize, control: &mut SolveControl<'_>) -> Result<MinimumMaterialRefit3<F>, FitError> {
        control.checkpoint("minimum-material-refit-start").map_err(physics)?;
        if matches!(options.objective, super::super::design::MinimumMaterialObjective3::MassKg { .. })
            != source.mass_kg.is_some() {
            return Err(ProjectedAlError::Invalid("refit objective disagrees with the source material inventory"));
        }
        admit(self, cases, &source.rho, options)?;
        options.optimizer.validate(candidate.cells())?;
        super::super::design::admit_objective(&candidate, options)?;
        if max_transfer_terms == 0 { return Err(ProjectedAlError::Invalid("empty material-refit transfer allowance")); }
        admit_reference_limits(&candidate, cases, options.limits).map_err(physics)?;
        let design = self.design(&source.rho, control).map_err(physics)?;
        if design.scales != source.scales || design.projected != source.projected_rho {
            return Err(ProjectedAlError::Invalid("source material design does not match its pipeline"));
        }
        let (inherited_rho, labels) = {
            let coarse = self.operator.adaptive(); let fine = candidate.operator.adaptive();
            let transfer = AdaptiveTransfer3::new(coarse, fine, max_transfer_terms, || poll(control)).map_err(physics)?;
            if !fine.leaves().iter().zip(transfer.parents()).any(|(leaf, &p)| leaf.level() > coarse.leaves()[p].level()) {
                return Err(ProjectedAlError::Invalid("minimum-material refit requires an enriched active space"));
            }
            let mut rho = Vec::with_capacity(candidate.cells());
            for &p in transfer.parents() {
                control.checkpoint("minimum-material-refit-transfer").map_err(physics)?;
                rho.push(source.rho[p]);
            }
            let labels = self.physical_regions().map(|regions|
                inherit_physical_regions(regions, transfer.parents(), control)).transpose().map_err(physics)?;
            (rho, labels)
        };
        if let Some(labels) = labels { candidate = candidate.with_physical_regions(labels).map_err(physics)?; }
        let fit = candidate.fit_reference_minimum_material(cases, &inherited_rho, options, steps, control)?;
        Ok(MinimumMaterialRefit3 { study: candidate, fit, inherited_rho })
    }
}
