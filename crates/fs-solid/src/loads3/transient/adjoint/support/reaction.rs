//! Derivatives of physical support force ON the solid in world components.
//!
//! R_c=(M*a+(eta*M+zeta*K)*v+K*q-fT-f_external)_c. A reaction objective
//! supplies weights w on fixed rows; its motion seed is (K' w,C' w,M' w).
//! Material and fixed-force terms are direct dependencies and must be added
//! to the implicit motion/trajectory pullback, including at initialization.
use super::super::material::{MaterialGradientBudget3, MaterialScaleGradient3};
use super::super::{TransientAdjointError3, check, number};
use super::{NodalMotionCotangent3, PrescribedTransientAdjoint3};
use crate::linear3::{TetLinearElasticProblem, TetThermalStrainField};
use crate::loads3::transient::{TetTransientRun3, TransientError3, response};
use fs_exec::Cx;

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionPullback3 {
    /// Seed the matching endpoint's motion/trajectory reverse with this.
    pub motion: NodalMotionCotangent3,
    /// Direct fixed-node force derivatives [objective/N], in input DOF order.
    /// Add these to the corresponding endpoint's mechanical-force gradient.
    /// The same -w applies independently to the full frozen thermal equivalent
    /// force; stiffness-induced thermal changes are already in `material`.
    pub force: Vec<(usize, f64)>,
    /// Direct material/damping terms at fixed full q/v/a. Add the existing
    /// motion reverse's implicit material terms for a complete design gradient.
    pub material: MaterialScaleGradient3,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Reverse a scalar objective's support-reaction observation at the EXACT
    /// accepted endpoint, including the solved initial acceleration. Supply
    /// sorted unique fixed DOFs and dJ/dR [objective/N], e.g. weighted residuals
    /// from `recover_response`. An omitted reaction row has zero weight.
    ///
    /// Original canonical K/M tiles provide the actual transpose. Full inertia,
    /// both Rayleigh terms and frozen thermal eigenstrain are retained without
    /// a new equilibrium solve, global full matrix, or finite difference.
    /// Works for stationary and driven supports, including fully prescribed
    /// bodies. Geometry, axes and free thermal strain remain fixed; this does
    /// not differentiate strain/acceptance gates or physical work reports.
    /// Refusal/cancellation returns no partial result and never changes motion.
    pub fn reaction_pullback(
        &self,
        reaction_cotangent: &[(usize, f64)],
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<ReactionPullback3, TransientAdjointError3> {
        let full = admit(self.problem, budget, cx)?;
        validate_rows(full, self.free_dofs(), reaction_cotangent, cx)?;
        let motion = self.nodal_motion(cx)?;
        self.reaction_pullback_nodal_at(&motion, reaction_cotangent, budget, cx)
    }

    pub(in crate::loads3::transient::adjoint) fn reaction_pullback_nodal_at(
        &self,
        motion: &[Vec<[f64; 3]>; 3],
        rows: &[(usize, f64)],
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<ReactionPullback3, TransientAdjointError3> {
        let full = admit(self.problem, budget, cx)?;
        let seed = reaction_seed(full, self.free_dofs(), rows, cx)?;
        let motion_seed = motion_pullback(
            self.problem,
            &seed,
            self.config.mass_damping_per_s,
            self.config.stiffness_damping_s,
            cx,
        )?;
        let material = direct_material(self, motion, &seed, cx)?;
        let mut force = Vec::with_capacity(rows.len());
        for &(dof, weight) in rows {
            check(cx)?;
            force.push((dof, -weight));
        }
        check(cx)?;
        Ok(ReactionPullback3 {
            motion: motion_seed,
            force,
            material,
        })
    }

    /// Internal recorded-endpoint hook; the caller owns association of the
    /// exact full q/v/a with this immutable run's material and thermal inputs.
    pub(in crate::loads3::transient::adjoint) fn reaction_material_pullback_nodal_at(
        &self,
        motion: &[Vec<[f64; 3]>; 3],
        rows: &[(usize, f64)],
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<MaterialScaleGradient3, TransientAdjointError3> {
        let full = admit(self.problem, budget, cx)?;
        let seed = reaction_seed(full, self.free_dofs(), rows, cx)?;
        direct_material(self, motion, &seed, cx)
    }
}

impl PrescribedTransientAdjoint3<'_, '_, '_, '_> {
    /// Full motion seed for support-reaction observations with fixed operators.
    /// This needs no primal motion tape. Rows are sorted unique fixed DOFs with
    /// finite dJ/dR weights; prepare-time full-DOF and element budgets apply.
    /// Chain the returned q/v/a through the matching motion reverse, and add
    /// direct force -w and material terms from the reaction's other dependencies.
    pub fn reaction_motion_pullback(
        &self,
        rows: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<NodalMotionCotangent3, TransientAdjointError3> {
        let seed = reaction_seed(self.full_dofs, self.free_dofs(), rows, cx)?;
        motion_pullback(
            self.problem,
            &seed,
            self.reduced.primal_config.mass_damping_per_s,
            self.reduced.primal_config.stiffness_damping_s,
            cx,
        )
    }
}

fn admit(
    problem: &TetLinearElasticProblem<'_>,
    budget: MaterialGradientBudget3,
    cx: &Cx<'_>,
) -> Result<usize, TransientAdjointError3> {
    check(cx)?;
    let full = problem
        .nodes_m
        .len()
        .checked_mul(3)
        .ok_or(TransientAdjointError3::Budget)?;
    if full > budget.maximum_full_dofs || problem.tetrahedra.len() > budget.maximum_elements {
        return Err(TransientAdjointError3::Budget);
    }
    Ok(full)
}

fn validate_rows(
    full: usize,
    free: &[usize],
    rows: &[(usize, f64)],
    cx: &Cx<'_>,
) -> Result<(), TransientAdjointError3> {
    check(cx)?;
    let mut previous = None;
    for &(dof, weight) in rows {
        check(cx)?;
        if dof >= full || previous.is_some_and(|p| p >= dof) || free.binary_search(&dof).is_ok() {
            return Err(TransientAdjointError3::Invalid(
                "reaction weights must name sorted unique fixed world-component DOFs",
            ));
        }
        number(weight)?;
        previous = Some(dof);
    }
    Ok(())
}

fn reaction_seed(
    full: usize,
    free: &[usize],
    rows: &[(usize, f64)],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, TransientAdjointError3> {
    validate_rows(full, free, rows, cx)?;
    let mut seed = vec![0.0; full];
    for &(dof, weight) in rows {
        check(cx)?;
        seed[dof] = weight;
    }
    Ok(seed)
}

fn motion_pullback(
    problem: &TetLinearElasticProblem<'_>,
    seed: &[f64],
    eta: f64,
    zeta: f64,
    cx: &Cx<'_>,
) -> Result<NodalMotionCotangent3, TransientAdjointError3> {
    let (mass, stiffness) = problem
        .transpose_full_kinematic_actions(seed, cx)
        .map_err(TransientError3::from)?;
    let mut motion = NodalMotionCotangent3 {
        q: vec![[0.0; 3]; problem.nodes_m.len()],
        v: vec![[0.0; 3]; problem.nodes_m.len()],
        a: vec![[0.0; 3]; problem.nodes_m.len()],
    };
    for dof in 0..seed.len() {
        if dof % 256 == 0 {
            check(cx)?;
        }
        motion.q[dof / 3][dof % 3] = stiffness[dof];
        motion.v[dof / 3][dof % 3] = number(eta * mass[dof] + zeta * stiffness[dof])?;
        motion.a[dof / 3][dof % 3] = mass[dof];
    }
    check(cx)?;
    Ok(motion)
}

fn direct_material(
    run: &TetTransientRun3<'_, '_, '_, '_>,
    motion: &[Vec<[f64; 3]>; 3],
    seed: &[f64],
    cx: &Cx<'_>,
) -> Result<MaterialScaleGradient3, TransientAdjointError3> {
    // The canonical contraction computes -lambda' dR, so lambda=-w gives the
    // DIRECT +w' dR contribution, including constrained-constrained entries.
    let mut lambda = Vec::with_capacity(run.problem.nodes_m.len());
    for (i, row) in seed.chunks_exact(3).enumerate() {
        if i % 256 == 0 {
            check(cx)?;
        }
        lambda.push([-row[0], -row[1], -row[2]]);
    }
    let thermal = run.thermal_strains.as_ref().map(|field| match field {
        response::FrozenThermal3::Uniform(value) => TetThermalStrainField::Uniform(value),
        response::FrozenThermal3::PerElement(values) => TetThermalStrainField::PerElement(values),
    });
    let (element_scales, rayleigh) = run
        .problem
        .transient_material_contraction(
            motion,
            &lambda,
            thermal,
            run.config.mass_damping_per_s,
            run.config.stiffness_damping_s,
            cx,
        )
        .map_err(TransientError3::from)?;
    check(cx)?;
    Ok(MaterialScaleGradient3 {
        element_scales,
        rayleigh,
    })
}
