//! One numerical owner for a complete affine load set on a fixed native grid.
use super::AffineLoadSpace3;
use crate::sdf3::{Sdf3Elasticity, checked_solve_preconditioned, failure};
use crate::{EvaluationStop, SolveControl};
use fs_ascent::quadratic_ball::{QuadraticBallMaximum, maximize_quadratic_ball};
use fs_sparse::precond::Precond;
use std::ops::ControlFlow;

const RECIPROCITY_TOLERANCE: f64 = 1e-9;

/// Affine load set on one fixed eroded discrete operator.
///
/// `gram`, `linear`, `offset` and `maximum` use the policy's single compliance
/// reference; forces, displacements and `worst_compliance` stay physical.
/// A hard/flat witness remains usable for values and frozen-goal DWR even when
/// the oracle refuses a differentiable envelope branch.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceGrid3 {
    /// Symmetric row-major G = F^T K^-1 F / reference.
    pub gram: Vec<f64>,
    /// h = F^T K^-1 b / reference.
    pub linear: Vec<f64>,
    /// c = b^T K^-1 b / reference, including native material self-weight.
    pub offset: f64,
    pub oracle: QuadraticBallMaximum,
    /// Complete normalized maximum, including `offset`.
    pub maximum: f64,
    pub reciprocity_defect: f64,
    pub witness_relative_defect: f64,
    pub witness_relative_residual: f64,
    pub offset_force: Vec<f64>,
    pub offset_displacement: Vec<f64>,
    pub basis_displacements: Vec<Vec<f64>>,
    pub worst_coefficients: Vec<f64>,
    pub worst_force: Vec<f64>,
    pub worst_displacement: Vec<f64>,
    /// Actual force/displacement work, before normalization.
    pub worst_compliance: f64,
}

pub(super) fn finite(value: f64, stage: &'static str) -> Result<f64, EvaluationStop> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(failure(stage))
    }
}

// Keep the physical value until division so normalization cannot silently
// erase a nonzero work or sensitivity. Overflow and underflow both refuse.
pub(super) fn normalize(
    value: f64,
    reference: f64,
    stage: &'static str,
) -> Result<f64, EvaluationStop> {
    let normalized = value / reference;
    if normalized.is_finite() && (value == 0.0 || normalized != 0.0) {
        Ok(normalized)
    } else {
        Err(failure(stage))
    }
}

// Raw basis forces are allowed to carry entries at strong clamps. Both the
// actual solves and every reported work/combination use the same masked space.
pub(super) fn force_work(
    force: &[f64],
    displacement: &[f64],
    fixed: &[bool],
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    let mut value = 0.0;
    for (d, (&f, &u)) in force.iter().zip(displacement).enumerate() {
        if d % 256 == 0 {
            control.checkpoint("sdf3-affine-load-space-work")?;
        }
        if !fixed[d / 3] {
            value += f * u;
        }
    }
    finite(value, "sdf3-affine-load-space-work")
}

pub(super) fn relative_difference(a: f64, b: f64, scale: f64) -> f64 {
    if scale == 0.0 {
        0.0
    } else {
        (a / scale - b / scale).abs()
    }
}

fn has_force(force: &[f64], fixed: &[bool]) -> bool {
    force
        .iter()
        .enumerate()
        .any(|(d, &f)| !fixed[d / 3] && f != 0.0)
}

fn witness_residual(
    op: &impl fs_solver::LinearOp,
    force: &[f64],
    displacement: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    control.checkpoint("sdf3-affine-load-space-witness-residual")?;
    let mut applied = vec![0.0; force.len()];
    op.apply(displacement, &mut applied);
    let scale = force.iter().map(|f| f.abs()).fold(0.0_f64, f64::max);
    let mut nr = 0.0;
    let mut nb = 0.0;
    for (i, (&f, &a)) in force.iter().zip(&applied).enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-load-space-witness-residual")?;
        }
        if !a.is_finite() || (scale == 0.0 && a != 0.0) {
            return Err(failure("sdf3-affine-load-space-witness-residual"));
        }
        if scale != 0.0 {
            let b = f / scale;
            let r = b - a / scale;
            nr += r * r;
            nb += b * b;
        }
    }
    let residual = if scale == 0.0 {
        0.0
    } else {
        fs_math::det::sqrt(nr / nb)
    };
    if !residual.is_finite() || residual > 1e-9 {
        return Err(failure("sdf3-affine-load-space-witness-residual"));
    }
    Ok(residual)
}

// No new equilibrium solves: frozen directions use the same admitted direct
// baseline/basis solutions and the same strong-clamp projection.
pub(super) fn combine_fields(
    fixed: &[bool],
    offset_force: &[f64],
    offset_displacement: &[f64],
    basis: &[&[f64]],
    basis_displacements: &[Vec<f64>],
    coefficients: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<(Vec<f64>, Vec<f64>), EvaluationStop> {
    let n = offset_force.len();
    let m = coefficients.len();
    let mut force = offset_force.to_vec();
    let mut displacement = offset_displacement.to_vec();
    for d in 0..n {
        if d % 256 == 0 {
            control.checkpoint("sdf3-affine-load-space-witness")?;
        }
        if fixed[d / 3] {
            force[d] = 0.0;
            displacement[d] = 0.0;
            continue;
        }
        for i in 0..m {
            force[d] += coefficients[i] * basis[i][d];
            displacement[d] += coefficients[i] * basis_displacements[i][d];
        }
        if !force[d].is_finite() || !displacement[d].is_finite() {
            return Err(failure("sdf3-affine-load-space-witness"));
        }
    }
    Ok((force, displacement))
}

// `offset_force` is the already assembled total physical baseline, including
// B*rho exactly once. The independent basis columns NEVER contain this offset.
pub(super) fn evaluate_grid<O: Sdf3Elasticity>(
    operator: &O,
    prepared: &impl Precond,
    offset_force: Vec<f64>,
    basis: &[&[f64]],
    policy: &AffineLoadSpace3,
    control: &mut SolveControl<'_>,
) -> Result<AffineLoadSpaceGrid3, EvaluationStop> {
    let m = policy.dimension;
    let offset_displacement = checked_solve_preconditioned(
        operator,
        prepared,
        &offset_force,
        1e-11,
        "sdf3-affine-load-space-offset",
        control,
    )?;
    let mut basis_displacements = Vec::with_capacity(m);
    for force in basis {
        let rhs = crate::sdf3::body::external_force(operator, force, None, control)?;
        basis_displacements.push(checked_solve_preconditioned(
            operator,
            prepared,
            &rhs,
            1e-11,
            "sdf3-affine-load-space-basis",
            control,
        )?);
    }
    let fixed = operator.fixed();
    let offset = normalize(
        force_work(&offset_force, &offset_displacement, fixed, control)?,
        policy.reference,
        "sdf3-affine-load-space-offset-work",
    )?;
    if offset < 0.0 || (offset == 0.0 && has_force(&offset_force, fixed)) {
        return Err(failure("sdf3-affine-load-space-offset-work"));
    }
    let mut gram = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..m {
            gram[i * m + j] = normalize(
                force_work(basis[i], &basis_displacements[j], fixed, control)?,
                policy.reference,
                "sdf3-affine-load-space-gram",
            )?;
        }
    }
    if (0..m)
        .any(|i| gram[i * m + i] < 0.0 || (gram[i * m + i] == 0.0 && has_force(basis[i], fixed)))
    {
        return Err(failure("sdf3-affine-load-space-gram"));
    }
    let scale = gram.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
    let mut reciprocity_defect = 0.0_f64;
    for i in 0..m {
        for j in i + 1..m {
            let (a, b) = (gram[i * m + j], gram[j * m + i]);
            reciprocity_defect = reciprocity_defect.max(relative_difference(a, b, scale));
            let mean = 0.5 * a + 0.5 * b;
            gram[i * m + j] = mean;
            gram[j * m + i] = mean;
        }
    }
    let mut linear = Vec::with_capacity(m);
    for i in 0..m {
        let a = normalize(
            force_work(basis[i], &offset_displacement, fixed, control)?,
            policy.reference,
            "sdf3-affine-load-space-linear",
        )?;
        let b = normalize(
            force_work(&offset_force, &basis_displacements[i], fixed, control)?,
            policy.reference,
            "sdf3-affine-load-space-linear",
        )?;
        // The energy Cauchy scale remains meaningful when the actual
        // cross term is zero or nearly cancels by orthogonality.
        let cross_scale = finite(
            (fs_math::det::sqrt(offset) * fs_math::det::sqrt(gram[i * m + i]))
                .max(a.abs())
                .max(b.abs()),
            "sdf3-affine-load-space-linear-scale",
        )?;
        reciprocity_defect = reciprocity_defect.max(relative_difference(a, b, cross_scale));
        linear.push(0.5 * a + 0.5 * b);
    }
    if !reciprocity_defect.is_finite() || reciprocity_defect > RECIPROCITY_TOLERANCE {
        return Err(failure("sdf3-affine-load-space-reciprocity"));
    }
    let mut stopped = None;
    let oracle = maximize_quadratic_ball(&gram, &linear, policy.oracle_options, || {
        match control.checkpoint("sdf3-affine-load-space-oracle") {
            Ok(()) => ControlFlow::Continue(()),
            Err(stop) => {
                stopped = Some(stop);
                ControlFlow::Break(())
            }
        }
    });
    if let Some(stop) = stopped {
        return Err(stop);
    }
    let oracle = oracle.map_err(|_| failure("sdf3-affine-load-space-oracle"))?;
    let maximum = finite(offset + oracle.value, "sdf3-affine-load-space-maximum")?;
    if maximum < 0.0 {
        return Err(failure("sdf3-affine-load-space-maximum"));
    }
    let worst_coefficients = oracle.coefficients.clone();
    let (worst_force, worst_displacement) = combine_fields(
        fixed,
        &offset_force,
        &offset_displacement,
        basis,
        &basis_displacements,
        &worst_coefficients,
        control,
    )?;
    let worst_compliance = force_work(&worst_force, &worst_displacement, fixed, control)?;
    let actual = normalize(
        worst_compliance,
        policy.reference,
        "sdf3-affine-load-space-witness-work",
    )?;
    let witness_scale = maximum
        .abs()
        .max(actual.abs())
        .max(offset.abs())
        .max(oracle.value.abs());
    let witness_relative_defect = relative_difference(maximum, actual, witness_scale);
    if worst_compliance < 0.0
        || (actual == 0.0 && has_force(&worst_force, fixed))
        || !witness_relative_defect.is_finite()
        || witness_relative_defect > (8.0 * m as f64 + 4.0) * RECIPROCITY_TOLERANCE
    {
        return Err(failure("sdf3-affine-load-space-witness-work"));
    }
    let witness_relative_residual =
        witness_residual(operator, &worst_force, &worst_displacement, control)?;
    Ok(AffineLoadSpaceGrid3 {
        gram,
        linear,
        offset,
        oracle,
        maximum,
        reciprocity_defect,
        witness_relative_defect,
        witness_relative_residual,
        offset_force,
        offset_displacement,
        basis_displacements,
        worst_coefficients,
        worst_force,
        worst_displacement,
        worst_compliance,
    })
}
