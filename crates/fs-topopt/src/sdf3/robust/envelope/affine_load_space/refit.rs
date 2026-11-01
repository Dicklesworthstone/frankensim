//! Rebuild a complete affine load model and restart the existing design owner.
//! Raw controls are inherited; physical assessment and reoptimization remain
//! different experiments. Adding material is not a self-weight restoration.
use super::AffineLoadSpace3;
use super::design::{
    AffineLoadSpaceDesign3, AffineLoadSpaceDesignGoal3, AffineLoadSpaceDesignIteration3,
    AffineLoadSpaceDesignOptions3, AffineLoadSpaceDesignPoint3, AffineLoadSpaceGradientCheck3,
    admit_model, physical_feasible,
};
use super::response::AffineResponseLimits3;
use super::response::reference::{self, ReferenceAffineResponseLimits3};
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlStop};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;

// The same numerical/refinement error owners apply; no alternative optimizer
// or physics error mapping is introduced for affine loading.
pub use crate::sdf3::robust::envelope::limits::refit::ComplianceCapRefitError3 as AffineLoadSpaceRefitError3;

#[derive(Debug, Clone, Copy)]
pub struct AffineLoadSpaceRefitOptions3 {
    pub design: AffineLoadSpaceDesignOptions3,
    /// Accepted optimizer updates allowed on this grid, independently of the
    /// cumulative physical solve budget shared by the complete campaign.
    pub accepted_steps: usize,
    pub max_transfer_terms: usize,
}
impl Default for AffineLoadSpaceRefitOptions3 {
    fn default() -> Self {
        Self {
            design: AffineLoadSpaceDesignOptions3::default(),
            accepted_steps: 100,
            max_transfer_terms: 2_000_000,
        }
    }
}

/// An uninstalled grid and its complete admitted optimization prefix.
/// Initial, latest accepted, best feasible and freshly rechecked fields are
/// distinct: the optimizer's final KKT report belongs to `accepted` only.
/// A failed baseline/gradient gate retains the grid with an error outcome and
/// no admitted initial point; later failures retain every accepted field below.
pub struct AffineLoadSpaceRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub params: SimpParams,
    pub policy: AffineLoadSpace3,
    pub goal: AffineLoadSpaceDesignGoal3,
    /// Supplied initial controls, or literal parent copies on an actual refit.
    pub inherited_rho: Vec<f64>,
    pub gradient_check: Option<AffineLoadSpaceGradientCheck3>,
    /// Candidate-grid baseline, suitable for comparison to its own endpoint.
    /// No cross-grid objective descent is implied by a successful refit.
    pub initial: Option<AffineLoadSpaceDesignPoint3>,
    pub accepted: Option<AffineLoadSpaceDesignPoint3>,
    pub best: Option<AffineLoadSpaceDesignPoint3>,
    pub history: Vec<AffineLoadSpaceDesignIteration3>,
    pub outcome: Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>>,
    /// Fresh COMPLETE affine-envelope solve after a completed AL allowance.
    /// A failure leaves `best` intact but does not authorize grid promotion.
    pub recheck: Option<Result<AffineLoadSpaceDesignPoint3, EvaluationStop>>,
    pub options: AffineLoadSpaceRefitOptions3,
    pub work: SolveWork,
}
impl<O: AdaptiveSdf3Elasticity> AffineLoadSpaceRefit3<O> {
    /// A bounded optimizer allowance completed and a fresh full-envelope
    /// recheck passes the actual physical constraint. AL stopping tolerance
    /// never excuses a cap excess. This is not stationarity, cross-grid
    /// descent, global design optimality or continuum feasibility.
    #[must_use]
    pub fn completed_and_rechecked(&self) -> bool {
        matches!(&self.outcome, Ok(report)
            if matches!(report.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit))
            && self.recheck.as_ref().is_some_and(|result| {
                result.as_ref().is_ok_and(|point| {
                    self.gradient_check.as_ref().is_some_and(|audit| {
                        audit.passed()
                            && audit.response_options
                                == point
                                    .evaluation
                                    .responses
                                    .as_ref()
                                    .map(|response| response.options)
                    }) && self.initial.as_ref().is_some_and(|initial| {
                        match (&initial.evaluation.responses, &point.evaluation.responses) {
                            (None, None) => true,
                            (Some(a), Some(b)) => {
                                a.options == b.options
                                    && a.witnesses.len() == b.witnesses.len()
                                    && a.witnesses.iter().zip(&b.witnesses).all(|(a, b)| {
                                        a.bound == b.bound
                                            && a.sense == b.sense
                                            && a.scale == b.scale
                                    })
                            }
                            _ => false,
                        }
                    }) && physical_feasible(&point.evaluation, self.goal, &self.policy)
                })
            })
    }
}

fn admit_reference<O: AdaptiveSdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    rho: &[f64],
    basis_count: usize,
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceRefitOptions3,
) -> Result<(), AffineLoadSpaceRefitError3> {
    admit_model(study, rho, basis_count, policy, goal, options.design)?;
    if options.max_transfer_terms == 0 {
        return Err(GoalRefinementError3::Invalid("empty affine-refit transfer allowance").into());
    }
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Parent-inherit RAW controls and protected physical labels onto the
    /// supplied finer operator. Preserve the source filter radius and SIMP/
    /// projection parameters; rebuild the common material acceleration on the
    /// candidate's native quadrature/constraints/material owner. Reference dead
    /// and basis laws are independently reintegrated, never prolonged nodally.
    ///
    /// The caller preserves the implicit domain, supports, law order/amplitudes
    /// and correlations, policy, goal and geometry budget. Native transfer
    /// additionally checks complete source material identity and real active-
    /// space enrichment before acceleration or load callbacks. The source is
    /// immutable and this operation never installs the candidate.
    ///
    /// The fine start may be infeasible even when its source was feasible.
    /// Self-weight prevents a general monotone density-restoration argument;
    /// the SAME projected-AL owner handles the new fine constraint directly.
    /// No coarse displacement, adjoint, gradient or optimizer state is reused.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_reference_affine_load_space<F: AdaptiveSdf3Elasticity>(
        &self,
        operator: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<F>, AffineLoadSpaceRefitError3> {
        self.refit_reference_affine_load_space_optional(
            operator, rho, dead, basis, policy, goal, options, None, control,
        )
    }

    /// Refit the same compliance/material goal AND every original response
    /// inequality over the complete load ball. Displacement observations are
    /// reintegrated on the candidate grid; support-force modes stay live during
    /// every material trial. No coarse nodal observation or active direction is
    /// prolonged, and a failed candidate never replaces its source.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_reference_affine_load_space_with_responses<F: AdaptiveSdf3Elasticity>(
        &self,
        operator: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        responses: ReferenceAffineResponseLimits3<'_>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<F>, AffineLoadSpaceRefitError3> {
        self.refit_reference_affine_load_space_optional(
            operator,
            rho,
            dead,
            basis,
            policy,
            goal,
            options,
            Some(responses),
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn refit_reference_affine_load_space_optional<F: AdaptiveSdf3Elasticity>(
        &self,
        operator: F,
        rho: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        responses: Option<ReferenceAffineResponseLimits3<'_>>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<F>, AffineLoadSpaceRefitError3> {
        control.checkpoint("sdf3-affine-refit-start")?;
        self.require_material_body_recipe("sdf3-affine-refit-material-recipe")?;
        admit_reference(self, rho, basis.len(), policy, goal, options)?;
        if let Some(family) = responses {
            reference::admit(self.operator(), family)?;
            reference::admit(&operator, family)?;
        }
        let (inherited_rho, regions) =
            {
                let transfer = controlled_goal3(control, "sdf3-affine-refit-transfer", |poll| {
                    AdaptiveTransfer3::new(
                        self.operator().adaptive(),
                        operator.adaptive(),
                        options.max_transfer_terms,
                        poll,
                    )
                    .map_err(Into::into)
                })?;
                if !transfer.fine().leaves().iter().zip(transfer.parents()).any(
                    |(leaf, &parent)| leaf.level() > transfer.coarse().leaves()[parent].level(),
                ) {
                    return Err(GoalRefinementError3::Invalid(
                        "affine refit requires active-space enrichment",
                    )
                    .into());
                }
                let mut inherited = Vec::with_capacity(transfer.parents().len());
                for (i, &parent) in transfer.parents().iter().enumerate() {
                    if i % 64 == 0 {
                        control.checkpoint("sdf3-affine-refit-transfer")?;
                    }
                    inherited.push(rho[parent]);
                }
                let regions = self
                    .physical_regions()
                    .map(|regions| {
                        crate::sdf3::regions::inherit_physical_regions(
                            regions,
                            transfer.parents(),
                            control,
                        )
                    })
                    .transpose()?;
                (inherited, regions)
            };
        let mut candidate = CutDensityStudy3::new(operator, self.filter_radius, self.params());
        if let Some(regions) = regions {
            candidate = candidate.with_physical_regions(regions)?;
        }
        // Fine quadrature may alter the protected-solid fraction or inventory.
        // Admit that actual candidate before rebuilding any acceleration law.
        admit_reference(
            &candidate,
            &inherited_rho,
            basis.len(),
            policy,
            goal,
            options,
        )?;
        self.inherit_material_acceleration(&mut candidate, control)?;
        candidate.solve_reference_affine_load_space_optional(
            &inherited_rho,
            dead,
            basis,
            policy,
            goal,
            options,
            responses,
            control,
        )
    }

    /// Initial-grid and refitted-grid fitting share this exact owner. Assemble
    /// the authored dead baseline and each pure basis once, then use the native
    /// affine evaluator for B*rho, all worst-load solves, full-chain gradient
    /// checks, both design goals, and the existing bounded AL line search.
    ///
    /// Input and load-assembly failures precede an admitted design. Once the
    /// optimization attempt begins, failures are retained in `outcome` beside
    /// the grid and any completed prefix; no fallible publication step discards
    /// that evidence. A fresh full-set recheck is required for promotion. No
    /// monotone restoration or coarse worst-direction substitution occurs.
    /// Every evaluation restores incoming operator scales, so the returned
    /// points carry their own eroded scales instead of silently changing them.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_reference_affine_load_space(
        self,
        incoming: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<O>, AffineLoadSpaceRefitError3> {
        self.solve_reference_affine_load_space_optional(
            incoming, dead, basis, policy, goal, options, None, control,
        )
    }

    /// Fit the original affine material design with pure reference response
    /// laws. All policies, bounds and native support requirements are admitted
    /// before load, observation or reaction callbacks. The same design owner
    /// checks the full joint gradient and archives only points satisfying the
    /// original physical compliance/volume and signed response bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_reference_affine_load_space_with_responses(
        self,
        incoming: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        responses: ReferenceAffineResponseLimits3<'_>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<O>, AffineLoadSpaceRefitError3> {
        self.solve_reference_affine_load_space_optional(
            incoming,
            dead,
            basis,
            policy,
            goal,
            options,
            Some(responses),
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_reference_affine_load_space_optional(
        self,
        incoming: &[f64],
        dead: ReferenceLoad3<'_>,
        basis: &[ReferenceLoad3<'_>],
        policy: &AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceRefitOptions3,
        responses: Option<ReferenceAffineResponseLimits3<'_>>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceRefit3<O>, AffineLoadSpaceRefitError3> {
        control.checkpoint("sdf3-affine-reference-start")?;
        admit_reference(&self, incoming, basis.len(), policy, goal, options)?;
        if let Some(family) = responses {
            reference::admit(self.operator(), family)?;
        }
        let params = self.params();
        let mut study = self;
        let dead_force = controlled_goal3(control, "sdf3-affine-reference-load", |poll| {
            study
                .operator()
                .adaptive()
                .reference_load(dead, poll)
                .map_err(Into::into)
        })?;
        let mut forces = Vec::with_capacity(basis.len());
        for &law in basis {
            forces.push(controlled_goal3(
                control,
                "sdf3-affine-reference-load",
                |poll| {
                    study
                        .operator()
                        .adaptive()
                        .reference_load(law, poll)
                        .map_err(Into::into)
                },
            )?);
        }
        let borrowed: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        let observations = responses
            .map(|family| reference::assemble(study.operator(), family, control))
            .transpose()?;
        let response_limits = observations.as_ref().map(|prepared| prepared.limits());
        let nodal_responses = responses
            .zip(response_limits.as_deref())
            .map(|(family, limits)| AffineResponseLimits3 {
                limits,
                options: family.options,
            });
        let (mut gradient_check, mut initial, mut accepted, mut best, mut recheck) =
            (None, None, None, None, None);
        let mut history = Vec::new();
        let outcome = {
            let attempt = match nodal_responses {
                Some(family) => AffineLoadSpaceDesign3::new_with_response_limits(
                    &mut study,
                    &dead_force,
                    &borrowed,
                    incoming,
                    policy.clone(),
                    goal,
                    options.design,
                    family,
                    control,
                ),
                None => AffineLoadSpaceDesign3::new(
                    &mut study,
                    &dead_force,
                    &borrowed,
                    incoming,
                    policy.clone(),
                    goal,
                    options.design,
                    control,
                ),
            };
            match attempt {
                Err(error) => Err(error),
                Ok(mut session) => {
                    gradient_check = Some(session.gradient_check().clone());
                    initial = Some(session.accepted().clone());
                    let outcome = session.run(options.accepted_steps);
                    accepted = Some(session.accepted().clone());
                    best = session.best_feasible().cloned();
                    history = session.history().to_vec();
                    if matches!(&outcome, Ok(report)
                    if matches!(report.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit))
                    {
                        recheck = Some(session.recheck_best().and_then(|point| {
                            point.ok_or(EvaluationStop::Breakdown {
                                stage: "sdf3-affine-refit-no-feasible-point",
                            })
                        }));
                    }
                    outcome
                }
            }
        };
        Ok(AffineLoadSpaceRefit3 {
            study,
            params,
            policy: policy.clone(),
            goal,
            inherited_rho: incoming.to_vec(),
            gradient_check,
            initial,
            accepted,
            best,
            history,
            outcome,
            recheck,
            options,
            work: control.work(),
        })
    }
}
