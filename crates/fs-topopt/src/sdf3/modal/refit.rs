//! Native material-mass reconstruction and source-preserving frequency refits.
//! Every grid uses the existing fixed-grid design owner. Only raw controls and
//! protected physical labels transfer; modes, derivatives and optimizer state
//! belong to the new grid's own complete solves.

use super::design::{
    Sdf3EigenfrequencyDesign3, Sdf3EigenfrequencyDesignIteration3,
    Sdf3EigenfrequencyDesignOptions3, Sdf3EigenfrequencyDesignPoint3,
    Sdf3EigenfrequencyGoal3, Sdf3EigenfrequencyGradientCheck3,
    admit, physical_feasible, same_goal,
};
use super::refinement::Sdf3EigenfrequencyAdaptiveError3;
use crate::modal::MatrixFreeEigenError;
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3, Sdf3Elasticity};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3};
use crate::{SimpParams, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlStop};
use fs_cutfem::elastic3::MaterialMass3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_spectral::generalized::GeneralizedError;

#[derive(Debug, Clone, Copy)]
pub struct Sdf3EigenfrequencyRefitOptions3 {
    pub design: Sdf3EigenfrequencyDesignOptions3,
    /// Accepted optimizer updates on this grid. All rejected trials, gradient
    /// probes and fresh rechecks share the caller's cumulative solve budget.
    pub accepted_steps: usize,
    pub max_transfer_terms: usize,
    /// Native consistent-mass stencil bound, applied anew on every grid.
    pub max_mass_terms: usize,
}

impl Default for Sdf3EigenfrequencyRefitOptions3 {
    fn default() -> Self {
        Self {
            design: Sdf3EigenfrequencyDesignOptions3::default(),
            accepted_steps: 100,
            max_transfer_terms: 2_000_000,
            max_mass_terms: 2_000_000,
        }
    }
}

// Public diagnostics are convenient to inspect, but changing their model
// metadata must not reinterpret an already-computed physical recheck.
#[derive(Clone, Copy)]
struct RefitBinding3 {
    params: SimpParams,
    goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyRefitOptions3,
    work: SolveWork,
}

/// An uninstalled study, its OWN native consistent mass, and every admitted
/// optimization checkpoint. Moving this value preserves the native mass owner;
/// the cache must never be reused on a reassembled operator.
///
/// Initial, accepted, best feasible and freshly rechecked fields are distinct.
/// The optimizer's final KKT report belongs to `accepted`, not necessarily to
/// the best point. A failed baseline or gradient gate retains the candidate
/// with an error outcome; later failures retain the completed accepted prefix.
pub struct Sdf3EigenfrequencyRefit3<O: Sdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub mass: MaterialMass3,
    pub params: SimpParams,
    pub goal: Sdf3EigenfrequencyGoal3,
    /// Supplied initial controls, or literal parent copies on an actual refit.
    pub inherited_rho: Vec<f64>,
    pub gradient_check: Option<Sdf3EigenfrequencyGradientCheck3>,
    /// Baseline on THIS grid, for comparison with its own optimized endpoint.
    pub initial: Option<Sdf3EigenfrequencyDesignPoint3>,
    pub accepted: Option<Sdf3EigenfrequencyDesignPoint3>,
    pub best: Option<Sdf3EigenfrequencyDesignPoint3>,
    pub history: Vec<Sdf3EigenfrequencyDesignIteration3>,
    pub outcome: Result<ProjectedAlReport, ProjectedAlError<MatrixFreeEigenError>>,
    /// Independent complete selected-cluster solve at the retained best point.
    /// A failure preserves that point but never authorizes grid promotion.
    pub recheck: Option<Result<Sdf3EigenfrequencyDesignPoint3, MatrixFreeEigenError>>,
    pub options: Sdf3EigenfrequencyRefitOptions3,
    pub work: SolveWork,
    binding: RefitBinding3,
}

fn same_params(a: SimpParams, b: SimpParams) -> bool {
    [a.e_min, a.penal, a.beta, a.eta] == [b.e_min, b.penal, b.beta, b.eta]
}

fn same_options(a: Sdf3EigenfrequencyRefitOptions3, b: Sdf3EigenfrequencyRefitOptions3) -> bool {
    let x = a.design;
    let y = b.design;
    let xs = x.frequency.solver;
    let ys = y.frequency.solver;
    let xo = x.optimizer;
    let yo = y.optimizer;
    (a.accepted_steps, a.max_transfer_terms, a.max_mass_terms)
        == (b.accepted_steps, b.max_transfer_terms, b.max_mass_terms)
        && (xs.count, xs.oversampling, xs.max_iterations, xs.seed)
            == (ys.count, ys.oversampling, ys.max_iterations, ys.seed)
        && [xs.relative_tolerance, xs.linear_tolerance, x.frequency.aggregation_beta,
            x.eigenvalue_reference, x.density_floor,
            x.gradient_check.step, x.gradient_check.relative_tolerance]
            == [ys.relative_tolerance, ys.linear_tolerance, y.frequency.aggregation_beta,
                y.eigenvalue_reference, y.density_floor,
                y.gradient_check.step, y.gradient_check.relative_tolerance]
        && (xo.max_dimension, xo.max_evaluations, xo.max_multiplier_updates, xo.max_backtracks)
            == (yo.max_dimension, yo.max_evaluations, yo.max_multiplier_updates, yo.max_backtracks)
        && [xo.tolerance, xo.initial_penalty, xo.max_penalty, xo.inner_tolerance,
            xo.min_spectral_step, xo.max_spectral_step]
            == [yo.tolerance, yo.initial_penalty, yo.max_penalty, yo.inner_tolerance,
                yo.min_spectral_step, yo.max_spectral_step]
}

impl<O: Sdf3Elasticity> Sdf3EigenfrequencyRefit3<O> {
    /// A bounded optimizer allowance completed under the ORIGINAL model and
    /// policy, and a fresh solve meets the original physical requirement.
    /// AL stopping tolerance never excuses a cap excess or a frequency deficit.
    /// This does not assert stationarity, cross-grid descent, spectral
    /// completeness, global design optimality or continuum feasibility.
    #[must_use]
    pub fn completed_and_rechecked(&self) -> bool {
        if !same_params(self.params, self.binding.params)
            || !same_params(self.study.params(), self.binding.params)
            || !same_goal(self.goal, self.binding.goal)
            || !same_options(self.options, self.binding.options)
            || self.work != self.binding.work
            || self.study.operator().validate_material_mass(&self.mass).is_err()
            || !matches!(&self.outcome, Ok(report)
                if matches!(report.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit))
        {
            return false;
        }
        let (Some(audit), Some(initial), Some(accepted), Some(best), Some(Ok(point))) =
            (&self.gradient_check, &self.initial, &self.accepted, &self.best, &self.recheck)
        else {
            return false;
        };
        let options = self.binding.options.design;
        audit.passed()
            && same_goal(audit.goal, self.binding.goal)
            && audit.eigenvalue_reference == options.eigenvalue_reference
            && audit.options.step == options.gradient_check.step
            && audit.options.relative_tolerance == options.gradient_check.relative_tolerance
            && initial.rho == self.inherited_rho
            && best.rho == point.rho
            && admit(&self.study, &accepted.rho, self.goal, options).is_ok()
            && admit(&self.study, &point.rho, self.goal, options).is_ok()
            && point.evaluation.objective.modal.values.len() == options.frequency.solver.count
            && physical_feasible(&point.evaluation, self.goal, options)
    }
}

fn admit_refit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    rho: &[f64],
    goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyRefitOptions3,
) -> Result<(), Sdf3EigenfrequencyAdaptiveError3> {
    admit(study, rho, goal, options.design)?;
    if options.max_transfer_terms == 0 || options.max_mass_terms == 0 {
        return Err(GoalRefinementError3::Invalid(
            "empty modal refit transfer or native mass allowance",
        ).into());
    }
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Inherit RAW controls and protected labels onto a genuine native
    /// enrichment. Rebuild the graph filter at the original physical radius,
    /// preserve SIMP/projection and use the supplied goal and complete policy,
    /// and prepare consistent inertia from the candidate's own retained
    /// material density, quadrature and strong/hanging constraint map.
    ///
    /// The native transfer checks source-material identity and geometric
    /// ancestry before any material acceleration callback. The caller retains
    /// the same implicit domain, supports, goal, frequency policy and physical
    /// configuration. Optional acceleration recipes are rebuilt, even though
    /// these infinitesimal vibrations have no prestress or geometric stiffness.
    ///
    /// A coarse feasible design may be infeasible on the new grid. The existing
    /// AL owner solves that new problem directly: adding density is not a
    /// general frequency-restoration rule. No coarse modes, gradients, best
    /// fields or optimizer multipliers are reused. The source is immutable,
    /// and even a completed candidate is never installed by this method.
    pub fn refit_eigenfrequency_design<F: AdaptiveSdf3Elasticity>(
        &self,
        operator: F,
        rho: &[f64],
        goal: Sdf3EigenfrequencyGoal3,
        options: Sdf3EigenfrequencyRefitOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<Sdf3EigenfrequencyRefit3<F>, Sdf3EigenfrequencyAdaptiveError3> {
        control.checkpoint("sdf3-modal-refit-start")?;
        self.require_material_body_recipe("sdf3-modal-refit-material-recipe")?;
        admit_refit(self, rho, goal, options)?;
        let (inherited_rho, regions) = {
            let transfer = controlled_goal3(control, "sdf3-modal-refit-transfer", |poll| {
                AdaptiveTransfer3::new(
                    self.operator().adaptive(),
                    operator.adaptive(),
                    options.max_transfer_terms,
                    poll,
                ).map_err(Into::into)
            })?;
            if !transfer.fine().leaves().iter().zip(transfer.parents()).any(
                |(leaf, &parent)| leaf.level() > transfer.coarse().leaves()[parent].level(),
            ) {
                return Err(GoalRefinementError3::Invalid(
                    "modal refit requires active-space enrichment",
                ).into());
            }
            // Use the checked native parent map; the legacy raw-transfer
            // helper's fixed 1e-3 floor would change this design's own bounds.
            let mut inherited = Vec::with_capacity(transfer.parents().len());
            for (i, &parent) in transfer.parents().iter().enumerate() {
                if i % 64 == 0 {
                    control.checkpoint("sdf3-modal-refit-transfer")?;
                }
                inherited.push(rho[parent]);
            }
            let regions = self.physical_regions().map(|regions| {
                crate::sdf3::regions::inherit_physical_regions(
                    regions, transfer.parents(), control,
                )
            }).transpose()?;
            (inherited, regions)
        };
        let mut candidate = CutDensityStudy3::new(operator, self.filter_radius, self.params());
        if let Some(regions) = regions {
            candidate = candidate.with_physical_regions(regions)?;
        }
        // Fine quadrature can change the protected-solid fraction. Reject an
        // impossible volume cap before invoking the acceleration recipe.
        admit_refit(&candidate, &inherited_rho, goal, options)?;
        self.inherit_material_acceleration(&mut candidate, control)?;
        candidate.solve_eigenfrequency_design(&inherited_rho, goal, options, control)
    }
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Initial and refitted grids share this complete fitting path. Prepare
    /// native consistent mass once on the actual operator, then reuse the
    /// existing full-chain gradient gate and fixed-grid AL design owner for
    /// either a volume-constrained maximum frequency or frequency-constrained
    /// minimum volume/kg. Physical kg inventory remains linear occupancy.
    ///
    /// Input/transfer/native-mass setup failures are returned directly. Once
    /// optimization starts, failures retain the candidate and its full accepted
    /// prefix in the result. Recheck failure also retains the best feasible
    /// fields. No fallible publication step can discard those diagnostics.
    /// Every evaluation restores the incoming operator scales on success and
    /// returned errors, and all work remains charged to the original control.
    pub fn solve_eigenfrequency_design(
        self,
        incoming: &[f64],
        goal: Sdf3EigenfrequencyGoal3,
        options: Sdf3EigenfrequencyRefitOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<Sdf3EigenfrequencyRefit3<O>, Sdf3EigenfrequencyAdaptiveError3> {
        control.checkpoint("sdf3-modal-refit-start")?;
        admit_refit(&self, incoming, goal, options)?;
        let mass = controlled_goal3(control, "sdf3-modal-refit-mass", |poll| {
            self.operator().prepare_material_mass(options.max_mass_terms, poll).map_err(Into::into)
        })?;
        let params = self.params();
        let mut study = self;
        let (mut gradient_check, mut initial, mut accepted, mut best, mut recheck) =
            (None, None, None, None, None);
        let mut history = Vec::new();
        let outcome = {
            let attempt = Sdf3EigenfrequencyDesign3::new(
                &mut study, &mass, incoming, goal, options.design, control,
            );
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
                            point.ok_or(GeneralizedError::Breakdown(
                                "modal refit has no physically feasible point to recheck",
                            ))
                        }));
                    }
                    outcome
                }
            }
        };
        let work = control.work();
        Ok(Sdf3EigenfrequencyRefit3 {
            study, mass, params, goal, inherited_rho: incoming.to_vec(),
            gradient_check, initial, accepted, best, history, outcome, recheck,
            options, work, binding: RefitBinding3 { params, goal, options, work },
        })
    }
}
