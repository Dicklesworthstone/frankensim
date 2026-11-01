//! Frequency design assessed on two actual native adaptive grids.
//!
//! Enrichment inherits projected PHYSICAL occupancy and reference-material
//! stiffness from source parents. It never refilters those fields. Independent
//! low-mode solves and a fine-mass projector comparison provide numerical grid
//! evidence. Local modal residuals guide refinement; they are not certified
//! eigenvalue errors, continuum bounds, or a proof of spectral completeness.

use super::design::{Sdf3EigenfrequencyDesignOptions3, Sdf3EigenfrequencyGoal3};
use super::{Sdf3EigenfrequencyOptions3, selected_modes};
use crate::eigenfreq::{mass_interp, smooth_min};
use crate::modal::{MatrixFreeEigenError, MatrixFreeEigenReport};
use crate::pipeline::SimpParams;
use crate::sdf3::inventory::{mass_weights, volume_fraction, weighted_occupancy};
use crate::sdf3::operator::native_controlled;
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3, Sdf3Elasticity};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3};
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_ascent::projected_al::ProjectedAlError;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::{ElasticityError3, MaterialMass3};
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::{GoalError3, GoalMarking3, dorfler3};
use fs_la::factor::cholesky;
use std::{collections::BTreeMap, ops::ControlFlow};

/// Preserve the native setup, spectral and optimization refusal separately.
#[derive(Debug)]
pub enum Sdf3EigenfrequencyAdaptiveError3 {
    Refinement(GoalRefinementError3),
    Spectral(MatrixFreeEigenError),
    Optimization(ProjectedAlError<MatrixFreeEigenError>),
}
impl From<GoalRefinementError3> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Refinement(e) }
}
impl From<MatrixFreeEigenError> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: MatrixFreeEigenError) -> Self { Self::Spectral(e) }
}
impl From<ProjectedAlError<MatrixFreeEigenError>> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: ProjectedAlError<MatrixFreeEigenError>) -> Self { Self::Optimization(e) }
}
impl From<EvaluationStop> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: EvaluationStop) -> Self { Self::Refinement(e.into()) }
}
impl From<ElasticityError3> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: ElasticityError3) -> Self { Self::Refinement(e.into()) }
}
impl From<GoalError3> for Sdf3EigenfrequencyAdaptiveError3 {
    fn from(e: GoalError3) -> Self { Self::Refinement(e.into()) }
}
impl std::fmt::Display for Sdf3EigenfrequencyAdaptiveError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D eigenfrequency adaptivity stopped: {self:?}")
    }
}
impl std::error::Error for Sdf3EigenfrequencyAdaptiveError3 {}

/// Resource bounds and admission of numerical mass-space identities.
#[derive(Debug, Clone, Copy)]
pub struct Sdf3EigenfrequencyRefinementOptions3 {
    pub max_transfer_terms: usize,
    pub max_mass_terms: usize,
    /// Maximum absolute defect of the fine mass-orthonormal modal Gram matrix.
    pub orthogonality_tolerance: f64,
}
impl Default for Sdf3EigenfrequencyRefinementOptions3 {
    fn default() -> Self {
        Self {
            max_transfer_terms: 2_000_000,
            max_mass_terms: 4_000_000,
            orthogonality_tolerance: 1e-8,
        }
    }
}

/// Values on one actual native operator, without a fictitious design gradient.
#[derive(Debug, Clone)]
pub struct Sdf3EigenfrequencyGrid3 {
    pub modal: MatrixFreeEigenReport,
    pub aggregate: f64,
    pub physical_rho: Vec<f64>,
    pub scales: Vec<f64>,
    pub volume_fraction: f64,
    pub physical_mass_kg: f64,
    pub next_eigenvalue: Option<f64>,
}

/// Numerical two-grid observation bound to its original design and requirement.
#[derive(Debug, Clone)]
pub struct Sdf3EigenfrequencyRefinement3 {
    pub goal: Sdf3EigenfrequencyGoal3,
    pub design_options: Sdf3EigenfrequencyDesignOptions3,
    pub params: SimpParams,
    pub rho: Vec<f64>,
    pub coarse: Sdf3EigenfrequencyGrid3,
    pub fine: Sdf3EigenfrequencyGrid3,
    /// Absolute aggregate change divided by the original eigenvalue reference.
    pub aggregate_change: f64,
    /// Largest ordered eigenvalue change divided by that reference. This
    /// compares spectral values, without pairing individual eigenvectors.
    pub eigenvalue_change: f64,
    pub volume_change: f64,
    pub mass_change_kg: f64,
    /// ||P_coarse-P_fine||_F/sqrt(2*m) in the fine native mass inner product.
    /// The coarse basis is prolonged and normalized by its full Gram matrix.
    /// Invariant to nonsingular coarse basis changes and fine modal rotations.
    pub subspace_change: f64,
    /// Absolute goal-weighted localized modal residual, aggregated to source
    /// leaves. A diagnostic marking signal, NOT a DWR or eigenvalue error bound.
    pub marking_mass: BTreeMap<Octant3, f64>,
    pub work: SolveWork,
}

fn ratio(value: f64, scale: f64) -> Option<f64> {
    let result = value / scale;
    (value.is_finite() && value >= 0.0 && scale.is_finite() && scale > 0.0
        && result.is_finite() && (value == 0.0 || result > 0.0)).then_some(result)
}
fn goal_scale(goal: Sdf3EigenfrequencyGoal3) -> f64 {
    match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } => volume_limit,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction, volume_scale,
        } => volume_scale,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg }, ..
        } => scale_kg,
    }
}
impl Sdf3EigenfrequencyRefinement3 {
    /// Check the ORIGINAL physical inequality on both independent grid solves.
    /// Numerical optimizer tolerances never relax the cap or frequency floor.
    #[must_use]
    pub fn physically_feasible(&self) -> bool {
        let reference = self.design_options.eigenvalue_reference;
        if ratio(1.0, reference).is_none() || ratio(1.0, goal_scale(self.goal)).is_none() {
            return false;
        }
        [&self.coarse, &self.fine].iter().all(|grid| {
            let valid = grid.aggregate.is_finite()
                && grid.volume_fraction.is_finite()
                && (0.0..=1.0).contains(&grid.volume_fraction)
                && grid.physical_mass_kg.is_finite() && grid.physical_mass_kg >= 0.0
                && !grid.modal.values.is_empty()
                && grid.modal.values.iter().all(|v| v.is_finite() && *v > 0.0);
            valid && match self.goal {
                Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } =>
                    volume_limit <= 1.0 && grid.volume_fraction <= volume_limit,
                Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } =>
                    grid.aggregate >= reference && grid.modal.values[0] >= reference,
            }
        })
    }

    /// Change in the actual original objective, with its original scale.
    /// The unused volume scale does not enter a kg objective.
    #[must_use]
    pub fn objective_change(&self) -> Option<f64> {
        match self.goal {
            Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } =>
                ratio(self.aggregate_change, 1.0),
            Sdf3EigenfrequencyGoal3::MinimumMaterial {
                objective: MinimumMaterialObjective3::VolumeFraction, volume_scale,
            } => ratio(self.volume_change, volume_scale),
            Sdf3EigenfrequencyGoal3::MinimumMaterial {
                objective: MinimumMaterialObjective3::MassKg { scale_kg }, ..
            } => ratio(self.mass_change_kg, scale_kg),
        }
    }

    /// Convenience numerical gate; none of these observations certify the
    /// continuum spectrum. The campaign may impose separate per-metric limits.
    #[must_use]
    pub fn stable_and_feasible(&self, tolerance: f64, subspace_tolerance: f64) -> bool {
        let mass_scale = self.coarse.physical_mass_kg.max(self.fine.physical_mass_kg);
        let relative_mass = if mass_scale == 0.0 && self.mass_change_kg == 0.0 {
            Some(0.0)
        } else {
            ratio(self.mass_change_kg, mass_scale)
        };
        tolerance.is_finite() && tolerance >= 0.0
            && subspace_tolerance.is_finite() && subspace_tolerance >= 0.0
            && self.physically_feasible()
            && [self.aggregate_change, self.eigenvalue_change, self.volume_change]
                .iter().all(|v| v.is_finite() && *v >= 0.0 && *v <= tolerance)
            && self.objective_change().is_some_and(|v| v <= tolerance)
            && relative_mass.is_some_and(|v| v <= tolerance)
            && self.subspace_change.is_finite() && self.subspace_change >= 0.0
            && self.subspace_change <= subspace_tolerance
    }

    pub fn mark(
        &self, theta: f64, max_marks: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, theta, max_marks, checkpoint)
    }
}

fn invalid(message: &'static str) -> Sdf3EigenfrequencyAdaptiveError3 {
    GoalRefinementError3::Invalid(message).into()
}
fn finite(value: f64) -> Result<f64, Sdf3EigenfrequencyAdaptiveError3> {
    if value.is_finite() { Ok(value) } else {
        Err(invalid("nonfinite modal refinement contraction"))
    }
}
fn normalized(value: f64, scale: f64) -> Result<f64, Sdf3EigenfrequencyAdaptiveError3> {
    ratio(value, scale).ok_or_else(|| invalid("unrepresentable modal refinement change"))
}

fn grid<O: Sdf3Elasticity>(
    op: &O, mass: &MaterialMass3, physical_rho: Vec<f64>,
    options: Sdf3EigenfrequencyOptions3, control: &mut SolveControl<'_>,
) -> Result<Sdf3EigenfrequencyGrid3, Sdf3EigenfrequencyAdaptiveError3> {
    let occupancy: Vec<_> = physical_rho.iter().map(|&r| mass_interp(r)).collect();
    let (modal, next_eigenvalue) = selected_modes(op, mass, &occupancy, options, control)?;
    let (aggregate, _) = smooth_min(&modal.values, options.aggregation_beta);
    let aggregate = finite(aggregate)?;
    let volume_fraction = volume_fraction(op, &physical_rho, control)?;
    let physical_mass_kg = weighted_occupancy(mass_weights(op)?, &physical_rho, control)?;
    Ok(Sdf3EigenfrequencyGrid3 {
        modal, aggregate, physical_rho, scales: op.scales().to_vec(),
        volume_fraction, physical_mass_kg, next_eigenvalue,
    })
}

fn dot(
    left: &[f64], right: &[f64], control: &mut SolveControl<'_>,
) -> Result<f64, Sdf3EigenfrequencyAdaptiveError3> {
    if left.len() != right.len() {
        return Err(invalid("modal refinement vector shape mismatch"));
    }
    let mut value = 0.0;
    for (i, (&a, &b)) in left.iter().zip(right).enumerate() {
        if i % 1024 == 0 { control.checkpoint("sdf3-eigen-refinement-projector")?; }
        value = a.mul_add(b, value);
    }
    finite(value)
}

fn compare<F: Sdf3Elasticity>(
    transfer: &AdaptiveTransfer3<'_>, fine_op: &F, fine_mass: &MaterialMass3,
    coarse: &Sdf3EigenfrequencyGrid3, fine: &Sdf3EigenfrequencyGrid3,
    reference: f64, aggregation_beta: f64, tolerance: f64, control: &mut SolveControl<'_>,
) -> Result<(f64, BTreeMap<Octant3, f64>), Sdf3EigenfrequencyAdaptiveError3> {
    let m = coarse.modal.modes.len();
    if m == 0 || m > 64 || m != fine.modal.modes.len()
        || m != coarse.modal.values.len() || m != fine.modal.values.len()
    {
        return Err(invalid("modal refinement requires matching nonempty clusters"));
    }
    let occupancy: Vec<_> = fine.physical_rho.iter().map(|&r| mass_interp(r)).collect();
    let mut prolonged = Vec::with_capacity(m);
    for mode in &coarse.modal.modes {
        prolonged.push(controlled_goal3(control, "sdf3-eigen-refinement-transfer", |poll| {
            transfer.prolongate(mode, poll).map_err(Into::into)
        })?);
    }

    // Validate the fine orthonormal basis in the ACTUAL fine consistent mass.
    // The spectral residual gate alone is not an orthogonality test.
    for j in 0..m {
        let image = native_controlled(control, "sdf3-eigen-refinement-mass", |poll| {
            fine_mass.apply(&occupancy, &fine.modal.modes[j], poll)
        })?;
        for i in 0..m {
            let entry = dot(&fine.modal.modes[i], &image, control)?;
            let expected = if i == j { 1.0 } else { 0.0 };
            if (entry - expected).abs() > tolerance {
                return Err(invalid("fine modal basis is not mass orthonormal"));
            }
        }
    }

    // H = V^T M_f V and C = V^T M_f W. Only m-by-m matrices are dense.
    let mut gram = vec![0.0; m * m];
    let mut cross = vec![0.0; m * m];
    for j in 0..m {
        let image = native_controlled(control, "sdf3-eigen-refinement-mass", |poll| {
            fine_mass.apply(&occupancy, &prolonged[j], poll)
        })?;
        for i in 0..m {
            gram[i * m + j] = dot(&prolonged[i], &image, control)?;
            cross[j * m + i] = dot(&fine.modal.modes[i], &image, control)?;
        }
    }
    let mut norms = Vec::with_capacity(m);
    for i in 0..m {
        let diagonal = gram[i * m + i];
        if diagonal <= 0.0 {
            return Err(invalid("prolonged modal basis has zero fine mass"));
        }
        norms.push(diagonal.sqrt());
    }
    // Diagonal equilibration avoids a needless scale sensitivity in Cholesky;
    // it does not alter either subspace or the resulting projector distance.
    let mut overlap = cross.clone();
    for i in 0..m {
        control.checkpoint("sdf3-eigen-refinement-projector")?;
        for j in 0..m {
            gram[i * m + j] = finite((gram[i * m + j] / norms[i]) / norms[j])?;
            overlap[i * m + j] = finite(overlap[i * m + j] / norms[i])?;
        }
    }
    for i in 0..m {
        for j in 0..i {
            let a = gram[i * m + j];
            let b = gram[j * m + i];
            if (a - b).abs() > tolerance {
                return Err(invalid("prolonged fine-mass Gram matrix is not symmetric"));
            }
            let symmetric = 0.5 * a + 0.5 * b;
            gram[i * m + j] = symmetric;
            gram[j * m + i] = symmetric;
        }
    }
    control.checkpoint("sdf3-eigen-refinement-projector")?;
    let factor = cholesky(&gram, m)
        .map_err(|_| invalid("prolonged modal basis is not independent in fine mass"))?;
    control.checkpoint("sdf3-eigen-refinement-projector")?;
    // Existing Cholesky owns the factorization. Forward substitution alone
    // computes L^-1 C, whose squared Frobenius norm is projector overlap.
    let mut captured = 0.0;
    for j in 0..m {
        for i in 0..m {
            control.checkpoint("sdf3-eigen-refinement-projector")?;
            let mut entry = overlap[i * m + j];
            for k in 0..i {
                entry = (-factor.l(i, k)).mul_add(overlap[k * m + j], entry);
            }
            let value = finite(entry / factor.l(i, i))?;
            overlap[i * m + j] = value;
            captured = finite(value.mul_add(value, captured))?;
        }
    }
    let dimension = m as f64;
    if captured > dimension * (1.0 + tolerance) {
        return Err(invalid("modal projector overlap exceeds its admitted dimension"));
    }
    let subspace_change = ((1.0 - captured / dimension).max(0.0)).sqrt();

    // Project each prolonged coarse vector onto the WHOLE fine cluster.
    // No arbitrary matching of eigenvectors is meaningful near multiplicity.
    let (_, weights) = smooth_min(&coarse.modal.values, aggregation_beta);
    let mut signed = vec![0.0; fine_op.cells()];
    for k in 0..m {
        let mut delta = prolonged[k].iter().map(|v| -*v).collect::<Vec<_>>();
        for j in 0..m {
            for (dof, (&mode, value)) in fine.modal.modes[j].iter().zip(&mut delta).enumerate() {
                if dof % 1024 == 0 { control.checkpoint("sdf3-eigen-refinement-marking")?; }
                *value = cross[k * m + j].mul_add(mode, *value);
            }
        }
        if delta.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite modal subspace correction"));
        }
        let stiffness = native_controlled(control, "sdf3-eigen-refinement-stiffness", |poll| {
            fine_op.scale_bilinear_forms(&delta, &prolonged[k], poll)
        })?;
        let inertia = native_controlled(control, "sdf3-eigen-refinement-inertia", |poll| {
            fine_mass.occupancy_pullback(&delta, &prolonged[k], poll)
        })?;
        for cell in 0..fine_op.cells() {
            if cell % 64 == 0 { control.checkpoint("sdf3-eigen-refinement-marking")?; }
            let elastic = finite(fine.scales[cell] * stiffness[cell])?;
            let kinetic = finite(coarse.modal.values[k] * occupancy[cell] * inertia[cell])?;
            signed[cell] = finite(weights[k].mul_add(finite(elastic - kinetic)?, signed[cell]))?;
        }
    }
    let mut marking_mass = BTreeMap::new();
    for (cell, &parent) in transfer.parents().iter().enumerate() {
        if cell % 64 == 0 { control.checkpoint("sdf3-eigen-refinement-marking")?; }
        // Sum SIGNED modal contributions before taking absolute values. For
        // repeated eigenvalues their equal-weight trace is rotation invariant.
        let value = normalized(signed[cell].abs(), reference)?;
        let sum = marking_mass.entry(transfer.coarse().leaves()[parent]).or_insert(0.0);
        *sum = finite(*sum + value)?;
    }
    Ok((subspace_change, marking_mass))
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Assess one raw design on this grid and an actual compatible enrichment.
    ///
    /// The source filter/projection runs once. Both native pencils then use
    /// that SAME physical material field, with literal parent copies on fine
    /// cells. Fine density filtering belongs to the later refit, not this
    /// discretization experiment. Native source density supplies both inertia
    /// and independent linear-occupancy kg inventory on each geometry.
    ///
    /// The selected low-mode cluster is independently solved and admitted on
    /// each grid. Mode amplitudes are algebraic, so material strain ceilings
    /// are not applied. This is vibration about the undeformed reference:
    /// force, self-weight, prescribed support motion and prestress do not
    /// enter the pencil. Both incoming stiffness fields are restored after
    /// success and every returned error; all work uses the original control.
    pub fn estimate_eigenfrequency_enrichment<F: AdaptiveSdf3Elasticity>(
        &mut self,
        coarse_mass: &MaterialMass3,
        enriched: &mut F,
        rho: &[f64],
        goal: Sdf3EigenfrequencyGoal3,
        design_options: Sdf3EigenfrequencyDesignOptions3,
        options: Sdf3EigenfrequencyRefinementOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<Sdf3EigenfrequencyRefinement3, Sdf3EigenfrequencyAdaptiveError3> {
        control.checkpoint("sdf3-eigen-refinement-setup")?;
        let reference = design_options.eigenvalue_reference;
        let beta = design_options.frequency.aggregation_beta;
        if rho.len() != self.cells() || rho.is_empty()
            || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || options.max_transfer_terms == 0 || options.max_mass_terms == 0
            || !options.orthogonality_tolerance.is_finite()
            || options.orthogonality_tolerance <= 0.0
            || options.orthogonality_tolerance >= 1.0
            || !beta.is_finite() || beta <= 0.0
            || ratio(1.0, reference).is_none()
            || ratio(1.0, goal_scale(goal)).is_none()
            || matches!(goal, Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit }
                if volume_limit > 1.0)
        {
            return Err(invalid("invalid modal refinement design, goal or resource limits"));
        }
        self.operator.validate_material_mass(coarse_mass)?;
        let design = self.design(rho, control)?;
        let previous_coarse = self.operator.scales().to_vec();
        let previous_fine = enriched.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let result = (|| {
            let initial = controlled_goal3(control, "sdf3-eigen-refinement-setup", |poll| {
                AdaptiveTransfer3::new(
                    self.operator.adaptive(), enriched.adaptive(), options.max_transfer_terms, poll,
                ).map_err(Into::into)
            })?;
            if !initial.fine().leaves().iter().zip(initial.parents())
                .any(|(leaf, &parent)| leaf.level() > initial.coarse().leaves()[parent].level())
            {
                return Err(invalid("modal refinement requires actual active-cell enrichment"));
            }
            let mut fine_rho = Vec::with_capacity(initial.parents().len());
            for (cell, &parent) in initial.parents().iter().enumerate() {
                if cell % 1024 == 0 { control.checkpoint("sdf3-eigen-refinement-transfer")?; }
                fine_rho.push(design.projected[parent]);
            }
            let fine_scales = initial.inherited_scales();
            drop(initial);
            enriched.set_scales(&fine_scales)?;
            // This checkpoint is inside the two-operator restoration scope.
            control.checkpoint("sdf3-eigen-refinement-fields")?;
            let fine_mass = controlled_goal3(control, "sdf3-eigen-refinement-mass-setup", |poll| {
                enriched.prepare_material_mass(options.max_mass_terms, poll).map_err(Into::into)
            })?;
            let coarse = grid(
                &self.operator, coarse_mass, design.projected, design_options.frequency, control,
            )?;
            let fine = grid(
                enriched, &fine_mass, fine_rho, design_options.frequency, control,
            )?;
            let transfer = controlled_goal3(control, "sdf3-eigen-refinement-transfer", |poll| {
                AdaptiveTransfer3::new(
                    self.operator.adaptive(), enriched.adaptive(), options.max_transfer_terms, poll,
                ).map_err(Into::into)
            })?;
            let (subspace_change, marking_mass) = compare(
                &transfer, enriched, &fine_mass, &coarse, &fine, reference, beta,
                options.orthogonality_tolerance, control,
            )?;
            let aggregate_change = normalized((fine.aggregate - coarse.aggregate).abs(), reference)?;
            let mut eigenvalue_change = 0.0_f64;
            for (&a, &b) in coarse.modal.values.iter().zip(&fine.modal.values) {
                eigenvalue_change = eigenvalue_change.max(normalized((a - b).abs(), reference)?);
            }
            let volume_change = finite((fine.volume_fraction - coarse.volume_fraction).abs())?;
            let mass_change_kg = finite((fine.physical_mass_kg - coarse.physical_mass_kg).abs())?;
            let report = Sdf3EigenfrequencyRefinement3 {
                goal, design_options, params: self.params(), rho: rho.to_vec(),
                coarse, fine, aggregate_change, eigenvalue_change, volume_change,
                mass_change_kg, subspace_change, marking_mass, work: control.work(),
            };
            if report.objective_change().is_none() {
                return Err(invalid("unrepresentable original modal objective change"));
            }
            control.checkpoint("sdf3-eigen-refinement-publish")?;
            Ok(report)
        })();
        enriched.set_scales(&previous_fine).expect("previous enriched scales valid");
        self.operator.set_scales(&previous_coarse).expect("previous source scales valid");
        result
    }
}
