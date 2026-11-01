//! Whole-ball stress assessment on two native grids and matched fine sites.
//!
//! The actual coarse/fine peak aggregates keep their own retained quadrature
//! sets. A second comparison prolongates every coarse affine field into the
//! fine space and samples BOTH complete field families with the fine material
//! and fine rules. Native transfer admits identical inherited reference cards;
//! this owner also requires literal parent occupancy and stiffness copies.
//!
//! At a fine point, let d0,dj be the differences of the reference stress
//! columns, already divided by the original stress cap. Then
//! `D = VM(d0) + sqrt(sum_j VM(dj)^2)` bounds VM(d0+sum a_j*dj) for ||a||<=1
//! by the triangle and Cauchy inequalities. Multiplying D by rho^q or k(rho)
//! bounds the relaxed or physical stress change. The reverse triangle
//! inequality also bounds the change of the corresponding VM magnitudes.
//! This uses every coefficient direction, not either grid's maximizing vector.
//!
//! These computed bounds are finite-sample numerical evidence, not intervals
//! or continuum stress certificates. The p-norm aggregate changes with sample
//! count even for an unchanged field. Requiring its change to resolve is
//! deliberately conservative and does not promise adaptive convergence.
use super::{
    AffineStressEvaluation3, AffineStressLimit3, controlled, measure_prepared, power,
    stress_product, von_mises,
};
use crate::sdf3::failure;
use crate::sdf3::robust::envelope::affine_load_space::AffineLoadSpaceGrid3;
use crate::sdf3::robust::envelope::affine_load_space::grid::{finite, normalize};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};
use fs_ascent::quadratic_ball::QuadraticBallOptions;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;

/// Bound on complete affine stress-model change at one original FINE point.
/// All three entries are dimensionless and divided by the original stress cap.
#[derive(Debug, Clone)]
pub struct AffineStressChange3 {
    pub cell: usize,
    pub position: [f64; 3],
    pub weight: f64,
    pub normalized_reference_change: f64,
    pub normalized_relaxed_change: f64,
    pub normalized_physical_change: f64,
}

#[derive(Debug, Clone)]
pub struct AffineStressRefinement3 {
    /// Independently remaximized on the actual coarse retained quadrature.
    pub coarse: AffineStressEvaluation3,
    /// Independently remaximized on the actual fine retained quadrature.
    pub fine: AffineStressEvaluation3,
    /// Absolute change of the optimizer's actual p-norm aggregate / cap.
    /// This is distinct from a pointwise field change because sampling changes.
    pub normalized_aggregate_change: f64,
    pub changes: Vec<AffineStressChange3>,
    pub maximum_reference_change: f64,
    pub maximum_relaxed_change: f64,
    pub maximum_physical_change: f64,
    /// Maximum of the relaxed/physical forms actually constrained. Reference
    /// stress alone is retained as a diagnostic, not an extra accuracy target.
    pub maximum_model_change: f64,
    /// Cumulative shared work; this assessment adds no elasticity solves.
    pub work: SolveWork,
}
impl AffineStressRefinement3 {
    /// Resolve the actual constrained aggregate AND all matched whole-ball
    /// field changes, while enforcing both grids' original physical caps.
    /// Numerical resolution tolerance never relaxes the physical stress cap.
    #[must_use]
    pub fn stable_and_feasible(&self, tolerance: f64) -> bool {
        tolerance.is_finite()
            && tolerance >= 0.0
            && !self.changes.is_empty()
            && self.coarse.physically_feasible()
            && self.fine.physically_feasible()
            && [
                self.normalized_aggregate_change,
                self.maximum_relaxed_change,
                self.maximum_physical_change,
                self.maximum_model_change,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0 && *v <= tolerance)
            && self.maximum_reference_change.is_finite()
            && self.maximum_reference_change >= 0.0
            && self.changes.iter().all(|point| {
                point.normalized_reference_change.is_finite()
                    && point.normalized_reference_change >= 0.0
                    && [
                        point.normalized_relaxed_change,
                        point.normalized_physical_change,
                    ]
                    .iter()
                    .all(|v| v.is_finite() && *v >= 0.0 && *v <= tolerance)
            })
    }
}

// Prefer subtraction before normalization so a loose cap does not erase a
// small physical change. Opposite extreme finite values may overflow their
// difference even when the normalized difference is representable.
fn difference(a: f64, b: f64, cap: f64) -> Result<f64, EvaluationStop> {
    let delta = a - b;
    if delta.is_finite() {
        normalize(delta, cap, "sdf3-affine-stress-refinement-change")
    } else {
        finite(
            normalize(a, cap, "sdf3-affine-stress-refinement-change")?
                - normalize(b, cap, "sdf3-affine-stress-refinement-change")?,
            "sdf3-affine-stress-refinement-change",
        )
    }
}

fn stress_difference(a: &[f64; 6], b: &[f64; 6], cap: f64) -> Result<f64, EvaluationStop> {
    let mut delta = [0.0; 6];
    for ((value, &a), &b) in delta.iter_mut().zip(a).zip(b) {
        *value = difference(a, b, cap)?;
    }
    finite(von_mises(&delta).0, "sdf3-affine-stress-refinement-change")
}

// Incremental scaled sum of squares, with no overflowing intermediate square.
fn add_square(value: f64, scale: &mut f64, sum: &mut f64) {
    if value > *scale {
        *sum = if *scale == 0.0 {
            1.0
        } else {
            *sum * (*scale / value).powi(2) + 1.0
        };
        *scale = value;
    } else if value > 0.0 {
        *sum += (value / *scale).powi(2);
    }
}

/// Borrow the actual native transfer and complete freshly solved grid families
/// inside the caller's existing material transaction. Sampling the difference
/// field itself would wrongly treat it as a physical experiment for source
/// strain admission; each real family is instead admitted separately first.
#[allow(clippy::too_many_arguments)]
pub(in crate::sdf3::robust::envelope::affine_load_space) fn assess(
    transfer: &AdaptiveTransfer3<'_>,
    coarse: &AffineLoadSpaceGrid3,
    fine: &AffineLoadSpaceGrid3,
    coarse_projected: &[f64],
    fine_projected: &[f64],
    params: SimpParams,
    limit: AffineStressLimit3,
    oracle_options: QuadraticBallOptions,
    control: &mut SolveControl<'_>,
) -> Result<AffineStressRefinement3, EvaluationStop> {
    control.checkpoint("sdf3-affine-stress-refinement-start")?;
    let m = coarse.basis_displacements.len();
    if m != fine.basis_displacements.len()
        || coarse_projected.len() != transfer.coarse().cells()
        || fine_projected.len() != transfer.fine().cells()
        || transfer.parents().len() != fine_projected.len()
    {
        return Err(failure("sdf3-affine-stress-refinement-shape"));
    }
    for (child, &parent) in transfer.parents().iter().enumerate() {
        control.checkpoint("sdf3-affine-stress-refinement-material")?;
        if parent >= coarse_projected.len()
            || fine_projected[child] != coarse_projected[parent]
            || transfer.fine().scales()[child] != transfer.coarse().scales()[parent]
        {
            return Err(failure("sdf3-affine-stress-refinement-material"));
        }
    }
    let coarse_report = measure_prepared(
        transfer.coarse(),
        coarse_projected,
        params,
        coarse,
        limit,
        oracle_options,
        control,
    )?;
    let fine_report = measure_prepared(
        transfer.fine(),
        fine_projected,
        params,
        fine,
        limit,
        oracle_options,
        control,
    )?;
    let normalized_aggregate_change =
        difference(fine_report.aggregate, coarse_report.aggregate, limit.limit)?.abs();

    control.checkpoint("sdf3-affine-stress-refinement-transfer")?;
    let prolonged_offset = controlled(control, |poll| {
        transfer.prolongate(&coarse.offset_displacement, poll)
    })?;
    let mut prolonged_basis = Vec::with_capacity(m);
    for field in &coarse.basis_displacements {
        control.checkpoint("sdf3-affine-stress-refinement-transfer")?;
        prolonged_basis.push(controlled(control, |poll| {
            transfer.prolongate(field, poll)
        })?);
    }
    let coarse_basis: Vec<_> = prolonged_basis.iter().map(Vec::as_slice).collect();
    let fine_basis: Vec<_> = fine.basis_displacements.iter().map(Vec::as_slice).collect();
    // Both samplings use the fine reference card, scale and point order.
    // Source-domain admission includes all load directions of EACH family.
    control.checkpoint("sdf3-affine-stress-refinement-sample")?;
    let prolonged_samples = controlled(control, |poll| {
        transfer.fine().affine_bulk_stress(
            &prolonged_offset,
            &coarse_basis,
            limit.options.max_points,
            poll,
        )
    })?;
    let fine_samples = controlled(control, |poll| {
        transfer.fine().affine_bulk_stress(
            &fine.offset_displacement,
            &fine_basis,
            limit.options.max_points,
            poll,
        )
    })?;
    let points = fine_samples.baseline.len();
    if points != fine_report.point_count
        || prolonged_samples.baseline.len() != points
        || fine_samples.basis_reference_stresses.len() != m
        || prolonged_samples.basis_reference_stresses.len() != m
        || fine_samples
            .basis_reference_stresses
            .iter()
            .chain(&prolonged_samples.basis_reference_stresses)
            .any(|v| v.len() != points)
    {
        return Err(failure("sdf3-affine-stress-refinement-samples"));
    }
    let mut changes = Vec::with_capacity(points);
    let mut maximum_reference_change = 0.0_f64;
    let mut maximum_relaxed_change = 0.0_f64;
    let mut maximum_physical_change = 0.0_f64;
    for (index, (a, b)) in fine_samples
        .baseline
        .iter()
        .zip(&prolonged_samples.baseline)
        .enumerate()
    {
        control.checkpoint("sdf3-affine-stress-refinement-change")?;
        if a.cell != b.cell
            || a.position != b.position
            || a.weight != b.weight
            || a.cell >= fine_projected.len()
        {
            return Err(failure("sdf3-affine-stress-refinement-samples"));
        }
        let baseline = stress_difference(&a.reference_stress, &b.reference_stress, limit.limit)?;
        let mut scale = 0.0;
        let mut sum = 0.0;
        for (a, b) in fine_samples
            .basis_reference_stresses
            .iter()
            .zip(&prolonged_samples.basis_reference_stresses)
        {
            let value = stress_difference(&a[index], &b[index], limit.limit)?;
            add_square(value, &mut scale, &mut sum);
        }
        let radius = stress_product(scale, fs_math::det::sqrt(sum))?;
        let normalized_reference_change =
            finite(baseline + radius, "sdf3-affine-stress-refinement-change")?;
        let normalized_relaxed_change = stress_product(
            power(fine_projected[a.cell], limit.options.relaxation_power),
            normalized_reference_change,
        )?;
        let normalized_physical_change = stress_product(
            transfer.fine().scales()[a.cell],
            normalized_reference_change,
        )?;
        maximum_reference_change = maximum_reference_change.max(normalized_reference_change);
        maximum_relaxed_change = maximum_relaxed_change.max(normalized_relaxed_change);
        maximum_physical_change = maximum_physical_change.max(normalized_physical_change);
        changes.push(AffineStressChange3 {
            cell: a.cell,
            position: a.position,
            weight: a.weight,
            normalized_reference_change,
            normalized_relaxed_change,
            normalized_physical_change,
        });
    }
    let maximum_model_change = maximum_relaxed_change.max(maximum_physical_change);
    control.checkpoint("sdf3-affine-stress-refinement-publish")?;
    Ok(AffineStressRefinement3 {
        coarse: coarse_report,
        fine: fine_report,
        normalized_aggregate_change,
        changes,
        maximum_reference_change,
        maximum_relaxed_change,
        maximum_physical_change,
        maximum_model_change,
        work: control.work(),
    })
}
