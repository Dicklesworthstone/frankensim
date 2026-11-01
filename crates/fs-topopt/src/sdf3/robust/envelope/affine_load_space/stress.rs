//! Material stress over the complete affine load ball on one native grid.
//!
//! Each retained point has its own adverse coefficient vector. Squared von
//! Mises stress is a positive-semidefinite quadratic in those coefficients, so
//! the existing bounded quadratic-ball oracle supplies the point envelope.
//! The unweighted p-norm of these envelopes bounds both relaxed and physical
//! stress at every retained point and load direction in exact arithmetic.
//! It is conservative relative to maximizing one global p-norm with a common
//! direction, and is neither a continuum bound nor a floating-point enclosure.
//!
//! All points share the original affine primal family. Aggregate cotangents
//! factor by [baseline, basis0, ...], requiring at most m+1 native adjoints.
//! Material self-weight contributes B^T z only for the baseline equation.
use super::AffineLoadSpaceGrid3;
use super::grid::{finite, normalize, relative_difference};
use crate::sdf3::stress::{self as stress_owner, StressOptions3, peak_scale, power, von_mises};
use crate::sdf3::{
    CutDensityStudy3, Design3, Sdf3Elasticity, checked_solve_preconditioned, failure,
};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};
use fs_ascent::quadratic_ball::{
    MAX_QUADRATIC_BALL_DIMENSION, QuadraticBallMaximum, QuadraticBallOptions,
    maximize_quadratic_ball,
};
use fs_cutfem::elastic3::ElasticityError3;
use fs_sparse::precond::Precond;
use std::ops::ControlFlow;

/// One explicit cap in the admitted reference material's stress units.
/// The constitutive tensor and material-domain ceiling remain native; this
/// declaration does not infer an allowable stress from an elastic modulus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineStressLimit3 {
    pub limit: f64,
    /// The existing qp relaxation and SampledPeakBound p-norm policy.
    /// max_cases bounds baseline+basis fields; max_points bounds their TOTAL
    /// point-field entries, not just the number of geometric sample locations.
    pub options: StressOptions3,
}
impl Default for AffineStressLimit3 {
    fn default() -> Self {
        Self {
            limit: 1.0,
            options: StressOptions3::default(),
        }
    }
}

/// Numerical worst-direction evidence at one original bulk quadrature point.
#[derive(Debug, Clone)]
pub struct AffineStressWitness3 {
    pub cell: usize,
    pub position: [f64; 3],
    pub weight: f64,
    /// Stress scale removed BEFORE constructing the dimensionless VM Gram.
    /// oracle.value omits the dimensionless baseline VM squared, as usual.
    pub normalization: f64,
    pub oracle: QuadraticBallMaximum,
    /// Native reference stress columns combined at oracle.coefficients.
    pub reference_stress: [f64; 6],
    /// Directly measured VM of the combined reference tensor.
    pub reference_von_mises: f64,
    pub relaxed_von_mises: f64,
    /// VM of the actually scaled combined tensor, independently evaluated.
    pub physical_von_mises: f64,
    /// Values from the complete quadratic maximum, including its constant.
    pub worst_reference_von_mises: f64,
    pub worst_relaxed_von_mises: f64,
    pub worst_physical_von_mises: f64,
    /// Relative identity defect between quadratic maximum and measured VM.
    pub relative_defect: f64,
}
impl AffineStressWitness3 {
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.oracle.coefficients
    }
}

#[derive(Debug, Clone)]
pub struct AffineStressEvaluation3 {
    pub limit: AffineStressLimit3,
    pub oracle_options: QuadraticBallOptions,
    /// Unweighted p-norm of BOTH relaxed and physical pointwise envelopes.
    pub aggregate: f64,
    pub constraint: f64,
    /// Raw-density derivative of aggregate/limit-1, not of aggregate alone.
    /// None retains valid values without claiming a smooth envelope branch.
    pub gradient: Option<Vec<f64>>,
    /// Actual combined-field sample maxima, separate from oracle values.
    pub sampled_relaxed_max: f64,
    pub sampled_physical_max: f64,
    pub witnesses: Vec<AffineStressWitness3>,
    pub point_count: usize,
    pub point_fields: usize,
    /// Aggregate adjoints in [baseline, basis0, ...] order; empty for values
    /// only or when any positive point contribution has an unresolved branch.
    pub adjoints: Vec<Vec<f64>>,
    pub work: SolveWork,
}
impl AffineStressEvaluation3 {
    /// Require the upper aggregate and BOTH physical stress measures at every
    /// retained point. No optimizer tolerance excuses an original cap excess.
    #[must_use]
    pub fn physically_feasible(&self) -> bool {
        let cap = self.limit.limit;
        !self.witnesses.is_empty()
            && cap.is_finite()
            && cap > 0.0
            && self.constraint.is_finite()
            && self.constraint <= 0.0
            && [
                self.aggregate,
                self.sampled_relaxed_max,
                self.sampled_physical_max,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0 && *v <= cap)
            && self.witnesses.iter().all(|w| {
                [
                    w.relaxed_von_mises,
                    w.physical_von_mises,
                    w.worst_relaxed_von_mises,
                    w.worst_physical_von_mises,
                ]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0 && *v <= cap)
            })
    }

    /// Numerical branch admission, independent of whether adjoints were
    /// requested. Positive hard/flat/near-hard contributions are not declared
    /// smooth. A zero point's p-th power has zero derivative for p>=2 when the
    /// whole aggregate is positive. The zero aggregate itself is not admitted.
    #[must_use]
    pub fn gradient_admissible(&self) -> bool {
        self.aggregate > 0.0
            && self.witnesses.iter().all(|w| {
                (w.worst_reference_von_mises == 0.0 && w.reference_von_mises == 0.0)
                    || w.oracle.gradient_admissible()
            })
    }
}

pub(super) fn admit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    limit: AffineStressLimit3,
) -> Result<(), EvaluationStop> {
    admit_policy(study.params, limit)
}

fn admit_policy(params: SimpParams, limit: AffineStressLimit3) -> Result<(), EvaluationStop> {
    stress_owner::admit_options(params, limit.options)
        .map_err(|_| failure("sdf3-affine-stress-policy"))?;
    if !limit.limit.is_finite() || limit.limit <= 0.0 || !(1.0 / limit.limit).is_finite() {
        return Err(failure("sdf3-affine-stress-cap"));
    }
    Ok(())
}

fn controlled<T>(
    control: &mut SolveControl<'_>,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, ElasticityError3>,
) -> Result<T, EvaluationStop> {
    control.checkpoint("sdf3-affine-stress-physics")?;
    let mut stopped = None;
    let result = action(
        &mut || match control.checkpoint("sdf3-affine-stress-physics") {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => {
                stopped = Some(error);
                ControlFlow::Break(())
            }
        },
    );
    if let Some(error) = stopped {
        return Err(error);
    }
    result.map_err(|_| failure("sdf3-affine-stress-physics"))
}

fn oracle(
    gram: &[f64],
    linear: &[f64],
    options: QuadraticBallOptions,
    control: &mut SolveControl<'_>,
) -> Result<QuadraticBallMaximum, EvaluationStop> {
    let mut stopped = None;
    let result = maximize_quadratic_ball(gram, linear, options, || {
        match control.checkpoint("sdf3-affine-stress-oracle") {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => {
                stopped = Some(error);
                ControlFlow::Break(())
            }
        }
    });
    if let Some(error) = stopped {
        return Err(error);
    }
    result.map_err(|_| failure("sdf3-affine-stress-oracle"))
}

// A factor of the EXISTING physical-Voigt von-Mises quadratic, not a second
// constitutive law. Inputs have already been normalized to absolute value <=1;
// dot products therefore cannot overflow merely by squaring physical stresses.
fn vm_coordinates(s: &[f64; 6]) -> [f64; 6] {
    let shear = fs_math::det::sqrt(3.0);
    [
        (s[0] - s[1]) * core::f64::consts::FRAC_1_SQRT_2,
        (s[1] - s[2]) * core::f64::consts::FRAC_1_SQRT_2,
        (s[2] - s[0]) * core::f64::consts::FRAC_1_SQRT_2,
        shear * s[3],
        shear * s[4],
        shear * s[5],
    ]
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn stress_product(a: f64, b: f64) -> Result<f64, EvaluationStop> {
    let value = a * b;
    if value.is_finite() && (a == 0.0 || b == 0.0 || value != 0.0) {
        Ok(value)
    } else {
        Err(failure("sdf3-affine-stress-range"))
    }
}

struct Measured {
    report: AffineStressEvaluation3,
    factors: Vec<(f64, f64)>,
    normalization: f64,
    root: f64,
}

#[allow(clippy::too_many_arguments)]
fn measure(
    operator: &impl Sdf3Elasticity,
    projected: &[f64],
    params: SimpParams,
    grid: &AffineLoadSpaceGrid3,
    limit: AffineStressLimit3,
    oracle_options: QuadraticBallOptions,
    control: &mut SolveControl<'_>,
) -> Result<Measured, EvaluationStop> {
    control.checkpoint("sdf3-affine-stress-assess")?;
    admit_policy(params, limit)?;
    if oracle_options.validate().is_err() {
        return Err(failure("sdf3-affine-stress-policy"));
    }
    let m = grid.basis_displacements.len();
    if m == 0
        || m > MAX_QUADRATIC_BALL_DIMENSION
        || m + 1 > limit.options.max_cases
        || projected.len() != operator.cells()
        || operator.scales().len() != projected.len()
    {
        return Err(failure("sdf3-affine-stress-field-shape"));
    }
    for field in std::iter::once(&grid.offset_displacement).chain(&grid.basis_displacements) {
        if field.len() != operator.n() {
            return Err(failure("sdf3-affine-stress-field-shape"));
        }
        for (i, value) in field.iter().enumerate() {
            if i % 256 == 0 {
                control.checkpoint("sdf3-affine-stress-field-shape")?;
            }
            if !value.is_finite() {
                return Err(failure("sdf3-affine-stress-field-shape"));
            }
        }
    }
    let mut relaxation = Vec::with_capacity(projected.len());
    let mut factors = Vec::with_capacity(projected.len());
    for (cell, &rho) in projected.iter().enumerate() {
        control.checkpoint("sdf3-affine-stress-coefficient")?;
        if !rho.is_finite() || !(0.0..=1.0).contains(&rho) {
            return Err(failure("sdf3-affine-stress-material"));
        }
        let stiffness = params.e_min + (1.0 - params.e_min) * power(rho, params.penal);
        if stiffness != operator.scales()[cell] {
            return Err(failure("sdf3-affine-stress-material"));
        }
        let relaxed = power(rho, limit.options.relaxation_power);
        relaxation.push(relaxed);
        factors.push(
            peak_scale(rho, relaxed, stiffness, params, limit.options)
                .map_err(|_| failure("sdf3-affine-stress-coefficient"))?,
        );
    }
    let basis: Vec<_> = grid.basis_displacements.iter().map(Vec::as_slice).collect();
    // Native sampling admits the WHOLE affine Mandel strain domain before
    // publishing any family. Direction columns are never treated as separate
    // physical experiments or allowed to bypass the retained material ceiling.
    let samples = controlled(control, |poll| {
        operator.affine_bulk_stress(
            &grid.offset_displacement,
            &basis,
            limit.options.max_points,
            poll,
        )
    })?;
    let point_count = samples.baseline.len();
    let point_fields = point_count
        .checked_mul(m + 1)
        .ok_or_else(|| failure("sdf3-affine-stress-point-fields"))?;
    if point_count == 0
        || point_fields > limit.options.max_points
        || samples.basis_reference_stresses.len() != m
        || samples
            .basis_reference_stresses
            .iter()
            .any(|v| v.len() != point_count)
    {
        return Err(failure("sdf3-affine-stress-point-fields"));
    }
    let mut witnesses = Vec::with_capacity(point_count);
    let mut sampled_relaxed_max = 0.0_f64;
    let mut sampled_physical_max = 0.0_f64;
    let mut normalization = 0.0_f64;
    let mut oracle_relaxed_max = 0.0_f64;
    let mut oracle_physical_max = 0.0_f64;
    // Reuse the small dense workspace across points. Only each oracle's own
    // bounded numerical report survives; material columns are not duplicated.
    let mut fields = Vec::with_capacity(m + 1);
    let mut coordinates = Vec::with_capacity(m + 1);
    let mut gram = vec![0.0; m * m];
    let mut linear = vec![0.0; m];
    for (index, point) in samples.baseline.iter().enumerate() {
        control.checkpoint("sdf3-affine-stress-witness")?;
        if point.cell >= projected.len()
            || !point.weight.is_finite()
            || point.weight <= 0.0
            || point.position.iter().any(|v| !v.is_finite())
        {
            return Err(failure("sdf3-affine-stress-sample"));
        }
        fields.clear();
        fields.push(point.reference_stress);
        fields.extend(samples.basis_reference_stresses.iter().map(|v| v[index]));
        if fields.iter().flatten().any(|v| !v.is_finite()) {
            return Err(failure("sdf3-affine-stress-sample"));
        }
        let stress_normalization = fields
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        if stress_normalization != 0.0 {
            for field in &mut fields {
                for value in field {
                    *value = normalize(*value, stress_normalization, "sdf3-affine-stress-range")?;
                }
            }
        }
        coordinates.clear();
        coordinates.extend(fields.iter().map(vm_coordinates));
        for i in 0..m {
            linear[i] = dot(&coordinates[0], &coordinates[i + 1]);
            for j in i..m {
                let value = dot(&coordinates[i + 1], &coordinates[j + 1]);
                gram[i * m + j] = value;
                gram[j * m + i] = value;
            }
        }
        let constant = dot(&coordinates[0], &coordinates[0]);
        let maximum = oracle(&gram, &linear, oracle_options, control)?;
        let complete = finite(constant + maximum.value, "sdf3-affine-stress-maximum")?;
        if complete < 0.0 {
            return Err(failure("sdf3-affine-stress-maximum"));
        }
        let worst_reference_von_mises =
            stress_product(stress_normalization, fs_math::det::sqrt(complete))?;
        let mut reference_stress = fields[0];
        for (field, &a) in fields[1..].iter().zip(&maximum.coefficients) {
            for (value, &direction) in reference_stress.iter_mut().zip(field) {
                *value = finite(*value + a * direction, "sdf3-affine-stress-witness")?;
            }
        }
        for value in &mut reference_stress {
            *value = stress_product(*value, stress_normalization)?;
        }
        let reference_von_mises =
            finite(von_mises(&reference_stress).0, "sdf3-affine-stress-witness")?;
        let relative_defect = relative_difference(
            worst_reference_von_mises,
            reference_von_mises,
            worst_reference_von_mises.max(reference_von_mises),
        );
        if !relative_defect.is_finite() || relative_defect > 8.0 * oracle_options.tolerance {
            return Err(failure("sdf3-affine-stress-witness-identity"));
        }
        let relaxed_von_mises = stress_product(relaxation[point.cell], reference_von_mises)?;
        let worst_relaxed_von_mises =
            stress_product(relaxation[point.cell], worst_reference_von_mises)?;
        let stiffness = operator.scales()[point.cell];
        let mut physical_stress = reference_stress;
        for value in &mut physical_stress {
            *value = stress_product(*value, stiffness)?;
        }
        let physical_von_mises =
            finite(von_mises(&physical_stress).0, "sdf3-affine-stress-witness")?;
        let worst_physical_von_mises = stress_product(stiffness, worst_reference_von_mises)?;
        sampled_relaxed_max = sampled_relaxed_max.max(relaxed_von_mises);
        sampled_physical_max = sampled_physical_max.max(physical_von_mises);
        oracle_relaxed_max = oracle_relaxed_max.max(worst_relaxed_von_mises);
        oracle_physical_max = oracle_physical_max.max(worst_physical_von_mises);
        normalization = normalization.max(stress_product(
            factors[point.cell].0,
            worst_reference_von_mises,
        )?);
        witnesses.push(AffineStressWitness3 {
            cell: point.cell,
            position: point.position,
            weight: point.weight,
            normalization: stress_normalization,
            oracle: maximum,
            reference_stress,
            reference_von_mises,
            relaxed_von_mises,
            physical_von_mises,
            worst_reference_von_mises,
            worst_relaxed_von_mises,
            worst_physical_von_mises,
            relative_defect,
        });
    }
    // Algebraic maximum rescaling stabilizes the same p-norm; its scale is not
    // frozen in differentiation. No volume or probability can dilute a point.
    let mut sum = 0.0;
    if normalization > 0.0 {
        for witness in &witnesses {
            control.checkpoint("sdf3-affine-stress-aggregate")?;
            let ratio = factors[witness.cell].0 * witness.worst_reference_von_mises / normalization;
            sum = finite(
                sum + power(ratio, limit.options.aggregation_power),
                "sdf3-affine-stress-aggregate",
            )?;
        }
    }
    let root = power(sum, 1.0 / limit.options.aggregation_power);
    let aggregate = stress_product(normalization, root)?;
    if [
        sampled_relaxed_max,
        sampled_physical_max,
        oracle_relaxed_max,
        oracle_physical_max,
    ]
    .iter()
    .any(|v| aggregate < *v)
    {
        return Err(failure("sdf3-affine-stress-upper-bound"));
    }
    let constraint = finite(
        normalize(aggregate, limit.limit, "sdf3-affine-stress-cap")? - 1.0,
        "sdf3-affine-stress-cap",
    )?;
    Ok(Measured {
        report: AffineStressEvaluation3 {
            limit,
            oracle_options,
            aggregate,
            constraint,
            gradient: None,
            sampled_relaxed_max,
            sampled_physical_max,
            witnesses,
            point_count,
            point_fields,
            adjoints: Vec::new(),
            work: control.work(),
        },
        factors,
        normalization,
        root,
    })
}

/// Same native point envelopes for coarse/fine assessment, without preparing,
/// filtering or solving any additional primal/adjoint. The caller owns the
/// current material transaction and the freshly admitted affine field family.
#[allow(clippy::too_many_arguments)]
pub(super) fn measure_prepared(
    operator: &impl Sdf3Elasticity,
    projected: &[f64],
    params: SimpParams,
    grid: &AffineLoadSpaceGrid3,
    limit: AffineStressLimit3,
    oracle_options: QuadraticBallOptions,
    control: &mut SolveControl<'_>,
) -> Result<AffineStressEvaluation3, EvaluationStop> {
    let mut report = measure(
        operator,
        projected,
        params,
        grid,
        limit,
        oracle_options,
        control,
    )?
    .report;
    control.checkpoint("sdf3-affine-stress-publish")?;
    report.work = control.work();
    Ok(report)
}

/// Differentiate the pointwise upper aggregate inside the parent's one eroded
/// material transaction/preparation. Every direct term uses the same retained
/// reference material and the same qp/physical peak coefficient as its value.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate_prepared<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    design: &Design3,
    prepared: &impl Precond,
    grid: &AffineLoadSpaceGrid3,
    limit: AffineStressLimit3,
    oracle_options: QuadraticBallOptions,
    control: &mut SolveControl<'_>,
) -> Result<AffineStressEvaluation3, EvaluationStop> {
    let Measured {
        mut report,
        factors,
        normalization,
        root,
    } = measure(
        &study.operator,
        &design.projected,
        study.params,
        grid,
        limit,
        oracle_options,
        control,
    )?;
    if report.gradient_admissible() {
        let mut occupancy = vec![0.0; study.cells()];
        let mut scale_gradient = vec![0.0; study.cells()];
        let mut cotangents = Vec::with_capacity(report.point_count);
        for witness in &report.witnesses {
            control.checkpoint("sdf3-affine-stress-derivative")?;
            let (coefficient, slope) = factors[witness.cell];
            let ratio = coefficient * witness.worst_reference_von_mises / normalization;
            let derivative = normalize(
                power(ratio / root, limit.options.aggregation_power - 1.0),
                limit.limit,
                "sdf3-affine-stress-derivative",
            )?;
            occupancy[witness.cell] = finite(
                occupancy[witness.cell] + derivative * witness.worst_reference_von_mises * slope,
                "sdf3-affine-stress-derivative",
            )?;
            let dvm = von_mises(&witness.reference_stress).1;
            let factor = finite(derivative * coefficient, "sdf3-affine-stress-derivative")?;
            let tangent = dvm.map(|v| factor * v);
            if tangent.iter().any(|v| !v.is_finite()) {
                return Err(failure("sdf3-affine-stress-derivative"));
            }
            cotangents.push(tangent);
        }
        let fields = grid.basis_displacements.len() + 1;
        let mut adjoints = Vec::with_capacity(fields);
        let mut basis_cotangents = Vec::with_capacity(cotangents.len());
        for field in 0..fields {
            control.checkpoint("sdf3-affine-stress-adjoint-rhs")?;
            let derivatives = if field == 0 {
                cotangents.as_slice()
            } else {
                basis_cotangents.clear();
                for (tangent, witness) in cotangents.iter().zip(&report.witnesses) {
                    control.checkpoint("sdf3-affine-stress-adjoint-rhs")?;
                    let a = witness.oracle.coefficients[field - 1];
                    let derivative = tangent.map(|v| a * v);
                    if derivative.iter().any(|v| !v.is_finite()) {
                        return Err(failure("sdf3-affine-stress-adjoint-rhs"));
                    }
                    basis_cotangents.push(derivative);
                }
                basis_cotangents.as_slice()
            };
            let rhs = controlled(control, |poll| {
                study
                    .operator
                    .reference_bulk_stress_pullback(derivatives, poll)
            })?;
            let z = checked_solve_preconditioned(
                &study.operator,
                prepared,
                &rhs,
                1e-11,
                "sdf3-affine-stress-adjoint",
                control,
            )?;
            let u = if field == 0 {
                &grid.offset_displacement
            } else {
                &grid.basis_displacements[field - 1]
            };
            let cross = controlled(control, |poll| {
                study.operator.scale_bilinear_forms(&z, u, poll)
            })?;
            if cross.len() != scale_gradient.len() {
                return Err(failure("sdf3-affine-stress-gradient"));
            }
            for (i, (value, contraction)) in scale_gradient.iter_mut().zip(cross).enumerate() {
                if i % 256 == 0 {
                    control.checkpoint("sdf3-affine-stress-gradient")?;
                }
                *value = finite(*value - contraction, "sdf3-affine-stress-gradient")?;
            }
            if field == 0 {
                study.add_material_body_gradient(&z, 1.0, &mut occupancy, control)?;
            }
            adjoints.push(z);
        }
        report.gradient =
            Some(study.physical_gradient(design, &scale_gradient, &occupancy, control)?);
        report.adjoints = adjoints;
    }
    control.checkpoint("sdf3-affine-stress-publish")?;
    report.work = control.work();
    Ok(report)
}
