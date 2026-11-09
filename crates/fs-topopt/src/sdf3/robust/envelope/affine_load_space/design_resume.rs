//! Restore the original affine optimizer from independently checked endpoints.
use super::*;
use crate::sdf3::response::limits::{BoundedResponse3, ResponseLimitSense3};
use fs_ascent::projected_al::ProjectedAlCheckpoint;
use fs_blake3::{ContentHash, DomainHasher};

/// Accepted optimizer state and the original full-chain gradient experiment.
/// No displacement, adjoint, factorization or load maximizer is imported.
///
/// `binding` covers the exact nodal loads, native material body map, inventory,
/// filter, regions and numerical/constraint declarations. The caller must ALSO
/// bind geometry, constitutive source, executable, solver setup and pure reaction
/// callbacks to their original source: endpoint equality cannot establish the
/// global identity of an arbitrary callback. Retained history is not replayed.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceDesignCheckpoint3 {
    pub optimizer: ProjectedAlCheckpoint,
    pub binding: ContentHash,
    pub best_feasible_density: Option<Vec<f64>>,
    pub history: Vec<AffineLoadSpaceDesignIteration3>,
    pub gradient_check: AffineLoadSpaceGradientCheck3,
    /// Endpoint reservations already included in `optimizer.work.evaluations`.
    /// A failed restore leaves these reservations charged on this checkpoint.
    pub restoration_evaluations: usize,
}
impl AffineLoadSpaceDesignCheckpoint3 {
    /// Re-solve the accepted point and, only when distinct, its feasible archive.
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(
            self.best_feasible_density
                .as_ref()
                .is_some_and(|rho| rho != &self.optimizer.point),
        )
    }
}

fn count(hash: &mut DomainHasher, value: usize) -> Result<(), EvaluationStop> {
    let value = u64::try_from(value).map_err(|_| failure("sdf3-affine-binding-count"))?;
    hash.update(&value.to_le_bytes());
    Ok(())
}
fn number(hash: &mut DomainHasher, value: f64) {
    hash.update(&value.to_bits().to_le_bytes());
}
fn vector(
    hash: &mut DomainHasher,
    values: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<(), EvaluationStop> {
    count(hash, values.len())?;
    for (i, &value) in values.iter().enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-binding")?;
        }
        number(hash, value);
    }
    Ok(())
}

// The prepared body-map fingerprint already binds native source attribution
// and every retained stencil coefficient. Do not substitute one sampled load.
#[allow(clippy::too_many_arguments)]
pub(super) fn binding<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    dead: &[f64],
    basis: &[&[f64]],
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
    constraints: AffineLoadSpaceConstraints3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<ContentHash, EvaluationStop> {
    control.checkpoint("sdf3-affine-binding")?;
    let mut hash = DomainHasher::new("frankensim.fs-topopt.affine-design-binding3.v1");
    count(&mut hash, study.cells())?;
    count(&mut hash, study.operator().n())?;
    number(&mut hash, study.filter_radius);
    vector(&mut hash, &study.mass, control)?;
    vector(&mut hash, &study.operator().volumes(), control)?;
    let fixed = study.operator().fixed();
    count(&mut hash, fixed.len())?;
    for (i, &fixed) in fixed.iter().enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-binding")?;
        }
        hash.update(&[u8::from(fixed)]);
    }
    let edges = study.operator().filter_edges();
    count(&mut hash, edges.len())?;
    for (i, &(a, b, distance)) in edges.iter().enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-binding")?;
        }
        count(&mut hash, a)?;
        count(&mut hash, b)?;
        number(&mut hash, distance);
    }
    hash.update(&[u8::from(study.physical_regions().is_some())]);
    if let Some(regions) = study.physical_regions() {
        count(&mut hash, regions.len())?;
        for region in regions {
            hash.update(&[match region {
                crate::sdf3::PhysicalRegion3::Design => 0,
                crate::sdf3::PhysicalRegion3::Solid => 1,
                crate::sdf3::PhysicalRegion3::Void => 2,
            }]);
        }
    }
    let masses = study.operator().reference_cell_masses_kg();
    hash.update(&[u8::from(masses.is_some())]);
    if let Some(masses) = masses {
        vector(&mut hash, masses, control)?;
    }
    let body = study.material_body_load_fingerprint();
    hash.update(&[u8::from(body.is_some())]);
    if let Some(body) = body {
        hash.update(body.as_bytes());
    }
    let penalty = study.operator().embedded_dirichlet_penalty();
    hash.update(&[u8::from(penalty.is_some())]);
    if let Some(penalty) = penalty {
        number(&mut hash, penalty);
    }
    vector(&mut hash, dead, control)?;
    count(&mut hash, basis.len())?;
    for force in basis {
        vector(&mut hash, force, control)?;
    }

    let params = study.params();
    let optimizer = options.optimizer;
    let oracle = policy.oracle_options();
    for value in [
        policy.dimension(),
        oracle.max_iterations,
        optimizer.max_dimension,
        optimizer.max_evaluations,
        optimizer.max_multiplier_updates,
        optimizer.max_backtracks,
    ] {
        count(&mut hash, value)?;
    }
    for value in [
        params.e_min,
        params.penal,
        params.beta,
        params.eta,
        policy.reference(),
        oracle.tolerance,
        oracle.separation_tolerance,
        options.robustness.threshold_offset,
        options.density_floor,
        options.gradient_check.step,
        options.gradient_check.relative_tolerance,
        optimizer.tolerance,
        optimizer.initial_penalty,
        optimizer.max_penalty,
        optimizer.inner_tolerance,
        optimizer.min_spectral_step,
        optimizer.max_spectral_step,
    ] {
        number(&mut hash, value);
    }
    match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => {
            hash.update(&[0]);
            number(&mut hash, volume_limit);
        }
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective,
            volume_scale,
        } => {
            hash.update(&[1]);
            number(&mut hash, volume_scale);
            match objective {
                MinimumMaterialObjective3::VolumeFraction => hash.update(&[0]),
                MinimumMaterialObjective3::MassKg { scale_kg } => {
                    hash.update(&[1]);
                    number(&mut hash, scale_kg);
                }
            }
        }
    }
    number(&mut hash, constraints.sharpness);
    hash.update(&[u8::from(constraints.responses.is_some())]);
    if let Some(responses) = constraints.responses {
        count(&mut hash, responses.options.max_limits)?;
        number(&mut hash, responses.options.sharpness);
        count(&mut hash, responses.limits.len())?;
        for limit in responses.limits {
            control.checkpoint("sdf3-affine-binding")?;
            hash.update(&[match limit.sense {
                ResponseLimitSense3::AtMost => 0,
                ResponseLimitSense3::AtLeast => 1,
            }]);
            number(&mut hash, limit.bound);
            number(&mut hash, limit.scale);
            match limit.response {
                BoundedResponse3::Displacement(q) => {
                    hash.update(&[0]);
                    vector(&mut hash, q, control)?;
                }
                // Function addresses are not portable source identities.
                BoundedResponse3::Reaction(_) => hash.update(&[1]),
            }
        }
    }
    hash.update(&[u8::from(constraints.stress.is_some())]);
    if let Some(stress) = constraints.stress {
        number(&mut hash, stress.limit);
        number(&mut hash, stress.options.relaxation_power);
        number(&mut hash, stress.options.aggregation_power);
        count(&mut hash, stress.options.max_cases)?;
        count(&mut hash, stress.options.max_points)?;
    }
    control.checkpoint("sdf3-affine-binding")?;
    Ok(hash.finalize())
}

fn archive_row(row: &AffineLoadSpaceDesignIteration3, goal: AffineLoadSpaceDesignGoal3) -> f64 {
    match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => row.maximum,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction,
            ..
        } => row.dilated_volume_fraction,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { .. },
            ..
        } => row.dilated_mass_kg.unwrap_or(f64::NAN),
    }
}

fn valid_row(
    row: &AffineLoadSpaceDesignIteration3,
    goal: AffineLoadSpaceDesignGoal3,
    policy: &AffineLoadSpace3,
    constrained: bool,
) -> bool {
    if !row.objective.is_finite()
        || row.objective < 0.0
        || !row.constraint.is_finite()
        || !row.maximum.is_finite()
        || row.maximum < 0.0
        || !row.worst_compliance.is_finite()
        || row.worst_compliance < 0.0
        || !(0.0..=1.0).contains(&row.dilated_volume_fraction)
        || row
            .dilated_mass_kg
            .is_some_and(|mass| !mass.is_finite() || mass < 0.0)
        || row.dilated_mass_kg.is_some() != goal.needs_mass()
    {
        return false;
    }
    let (objective, base, physical) = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            row.maximum,
            row.dilated_volume_fraction / volume_limit - 1.0,
            row.dilated_volume_fraction <= volume_limit,
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective,
            volume_scale,
        } => {
            let objective = match objective {
                MinimumMaterialObjective3::VolumeFraction => {
                    row.dilated_volume_fraction / volume_scale
                }
                MinimumMaterialObjective3::MassKg { scale_kg } => {
                    row.dilated_mass_kg.unwrap_or(f64::NAN) / scale_kg
                }
            };
            (
                objective,
                row.maximum - 1.0,
                row.maximum <= 1.0 && row.worst_compliance <= policy.reference(),
            )
        }
    };
    row.objective == objective
        && if constrained {
            row.constraint >= base
        } else {
            row.constraint == base
        }
        && (!row.feasible || (physical && row.constraint <= 0.0))
        && (constrained || row.feasible == (physical && row.constraint <= 0.0))
}

impl<'a, 'callback, O: Sdf3Elasticity> AffineLoadSpaceDesign3<'a, 'callback, O> {
    #[must_use]
    pub fn checkpoint(&self) -> AffineLoadSpaceDesignCheckpoint3 {
        AffineLoadSpaceDesignCheckpoint3 {
            optimizer: self.state.checkpoint(),
            binding: self.binding,
            best_feasible_density: self.best.as_ref().map(|point| point.rho.clone()),
            history: self.history.clone(),
            gradient_check: self.check.clone(),
            restoration_evaluations: self.restoration_evaluations,
        }
    }

    /// Restore a compliance/volume or compliance/inventory problem. Endpoint
    /// solves use the original control; initial finite differences and prior
    /// searches are not replayed. See `restore_with_constraints` for accounting.
    #[allow(clippy::too_many_arguments)]
    pub fn restore(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        checkpoint: &mut AffineLoadSpaceDesignCheckpoint3,
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::restore_with_constraints(
            study,
            dead_offset,
            basis,
            checkpoint,
            policy,
            goal,
            options,
            AffineLoadSpaceConstraints3::default(),
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn restore_with_response_limits(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        checkpoint: &mut AffineLoadSpaceDesignCheckpoint3,
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        responses: AffineResponseLimits3<'a>,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::restore_with_constraints(
            study,
            dead_offset,
            basis,
            checkpoint,
            policy,
            goal,
            options,
            AffineLoadSpaceConstraints3::from_responses(responses),
            control,
        )
    }

    /// Re-solve accepted and distinct feasible endpoints with every original
    /// response/strength inequality, checking complete objective/constraint
    /// derivatives and physical archive measurements exactly.
    ///
    /// After structural admission, reserve the full `restoration_cost()` in BOTH
    /// checkpoint counters before any endpoint work. Failure leaves only those
    /// counters changed; persist them and the control's actual physical work
    /// before retrying. An exhausted reservation or structural refusal changes
    /// no checkpoint. Original cumulative budgets remain the caller's duty.
    /// Every endpoint evaluation restores incoming operator scales on all exits.
    #[allow(clippy::too_many_arguments)]
    pub fn restore_with_constraints(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        checkpoint: &mut AffineLoadSpaceDesignCheckpoint3,
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        constraints: AffineLoadSpaceConstraints3<'a>,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        let invalid = ProjectedAlError::Invalid;
        control
            .checkpoint("sdf3-affine-restore-start")
            .map_err(ProjectedAlError::Evaluation)?;
        admit(
            study,
            dead_offset,
            basis,
            &checkpoint.optimizer.point,
            &policy,
            goal,
            options,
        )?;
        admit_constraints(study, options, constraints)?;
        let binding = binding(
            study,
            dead_offset,
            basis,
            &policy,
            goal,
            options,
            constraints,
            control,
        )
        .map_err(ProjectedAlError::Evaluation)?;
        if binding != checkpoint.binding {
            return Err(invalid(
                "affine checkpoint changed its loads, material map or model policy",
            ));
        }
        let n = study.cells();
        let check = &checkpoint.gradient_check;
        let p = study.params();
        let q = check.params;
        if check.goal != goal
            || check.constraints != constraints.policy()
            || [p.e_min, p.penal, p.beta, p.eta].map(f64::to_bits)
                != [q.e_min, q.penal, q.beta, q.eta].map(f64::to_bits)
            || check.options.step != options.gradient_check.step
            || check.options.relative_tolerance != options.gradient_check.relative_tolerance
            || !check.passed()
            || check.directions > 2
            || !check.baseline_objective.is_finite()
            || check.baseline_objective < 0.0
            || !check.baseline_constraint.is_finite()
            || check.work.linear_solves == 0
        {
            return Err(invalid(
                "affine checkpoint changed its original gradient experiment",
            ));
        }
        let lower = vec![options.density_floor; n];
        let upper = vec![1.0; n];
        let original = ProjectedAlState::try_restore::<EvaluationStop>(
            checkpoint.optimizer.clone(),
            &lower,
            &upper,
            options.optimizer,
        )?;
        let work = original.work();
        let minimum = work
            .iterations
            .checked_add(work.rejected_trials)
            .and_then(|v| v.checked_add(checkpoint.restoration_evaluations))
            .and_then(|v| v.checked_add(1))
            .ok_or(invalid("affine checkpoint evaluation count overflow"))?;
        if work.evaluations < minimum
            || work.iterations.checked_add(1) != Some(checkpoint.history.len())
            || checkpoint.history.iter().enumerate().any(|(i, row)| {
                row.iteration != i
                    || !valid_row(row, goal, &policy, constraints.policy().has_limits())
            })
            || checkpoint.history.first().is_none_or(|row| {
                row.objective != check.baseline_objective
                    || row.constraint != check.baseline_constraint
            })
            || checkpoint.history.last().is_none_or(|row| {
                row.objective != original.sample().objective
                    || row.constraint != original.sample().constraint
            })
        {
            return Err(invalid(
                "invalid affine checkpoint history, gradient evidence or work",
            ));
        }
        let best_row = checkpoint
            .history
            .iter()
            .filter(|row| row.feasible)
            .fold(
                None,
                |best: Option<&AffineLoadSpaceDesignIteration3>, row| {
                    if best.is_none_or(|best| archive_row(row, goal) < archive_row(best, goal)) {
                        Some(row)
                    } else {
                        best
                    }
                },
            )
            .cloned();
        if best_row.is_some() != checkpoint.best_feasible_density.is_some() {
            return Err(invalid(
                "affine checkpoint lost its physically feasible incumbent",
            ));
        }
        if let Some(rho) = &checkpoint.best_feasible_density {
            admit(study, dead_offset, basis, rho, &policy, goal, options)?;
        }
        let cost = checkpoint.restoration_cost();
        let mut reserved = checkpoint.optimizer.clone();
        reserved.work.evaluations = reserved
            .work
            .evaluations
            .checked_add(cost)
            .ok_or(invalid("affine restoration evaluation overflow"))?;
        let restoration_evaluations = checkpoint
            .restoration_evaluations
            .checked_add(cost)
            .ok_or(invalid("affine restoration counter overflow"))?;
        let state = ProjectedAlState::try_restore::<EvaluationStop>(
            reserved.clone(),
            &lower,
            &upper,
            options.optimizer,
        )?;
        checkpoint.optimizer = reserved;
        checkpoint.restoration_evaluations = restoration_evaluations;
        control
            .checkpoint("sdf3-affine-restore-reserved")
            .map_err(ProjectedAlError::Evaluation)?;
        let mut accepted = None;
        let evaluated = sample(
            study,
            state.point(),
            dead_offset,
            basis,
            &policy,
            goal,
            options,
            constraints,
            &mut accepted,
            control,
        )
        .map_err(ProjectedAlError::Evaluation)?;
        let accepted = accepted.ok_or(invalid("missing restored affine accepted endpoint"))?;
        if evaluated.as_ref() != Some(state.sample())
            || checkpoint.history.last()
                != Some(&row(
                    &accepted,
                    work.iterations,
                    state.sample(),
                    goal,
                    &policy,
                ))
        {
            return Err(invalid(
                "restored affine endpoint changed its physical response or derivatives",
            ));
        }
        let best = match &checkpoint.best_feasible_density {
            None => None,
            Some(rho) if rho == &accepted.rho => Some(accepted.clone()),
            Some(rho) => {
                control
                    .checkpoint("sdf3-affine-restore-incumbent")
                    .map_err(ProjectedAlError::Evaluation)?;
                let mut best = None;
                let evaluated = sample(
                    study,
                    rho,
                    dead_offset,
                    basis,
                    &policy,
                    goal,
                    options,
                    constraints,
                    &mut best,
                    control,
                )
                .map_err(ProjectedAlError::Evaluation)?;
                if evaluated.is_none() {
                    return Err(invalid(
                        "restored affine incumbent lost its smooth gradient",
                    ));
                }
                best
            }
        };
        if let (Some(best), Some(expected)) = (&best, &best_row) {
            let sample = optimization_sample(&best.evaluation, goal, true)
                .map_err(ProjectedAlError::Evaluation)?
                .ok_or(invalid(
                    "restored affine incumbent lost its smooth gradient",
                ))?;
            if row(best, expected.iteration, &sample, goal, &policy) != *expected {
                return Err(invalid(
                    "restored affine incumbent changed its physical bounds or inventory",
                ));
            }
        }
        control
            .checkpoint("sdf3-affine-restored")
            .map_err(ProjectedAlError::Evaluation)?;
        Ok(Self {
            study,
            dead: dead_offset,
            basis,
            constraints,
            control,
            policy,
            goal,
            options,
            check: checkpoint.gradient_check.clone(),
            state,
            accepted,
            best,
            history: checkpoint.history.clone(),
            binding,
            restoration_evaluations,
        })
    }
}

#[cfg(test)]
#[path = "design_resume/tests.rs"]
mod tests;
