//! Continuous affine uncertainty: f = b_dead + B*rho_eroded + F*a, ||a|| <= 1.
//!
//! The native bounded quadratic-ball owner supplies the global load witness.
//! This numerically solves the exact discrete quadratic formulation, with
//! numerical optimality evidence, rather than a smooth spectral surrogate.
//! A selected hard-case witness is useful evidence but is not a differentiable
//! objective. Only a numerically separated regular maximizer gets a gradient.
pub mod design;

use super::super::RobustOptions3;
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, checked_solve_preconditioned, failure};
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_ascent::quadratic_ball::{
    MAX_QUADRATIC_BALL_DIMENSION, QuadraticBallMaximum, QuadraticBallOptions,
    maximize_quadratic_ball,
};
use std::ops::ControlFlow;

const RECIPROCITY_TOLERANCE: f64 = 1e-9;

/// Fixed finite-dimensional uncertainty policy. Basis columns carry physical
/// amplitudes and correlations; they are not independent load probabilities.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineLoadSpace3 {
    dimension: usize,
    reference: f64,
    oracle_options: QuadraticBallOptions,
}
impl AffineLoadSpace3 {
    /// One positive compliance reference normalizes the whole affine set.
    /// The baseline, basis, reference and numerical policy stay fixed in a run.
    pub fn new(dimension: usize, reference: f64) -> Result<Self, &'static str> {
        if dimension == 0
            || dimension > MAX_QUADRATIC_BALL_DIMENSION
            || !reference.is_finite()
            || reference <= 0.0
            || !(1.0 / reference).is_finite()
            || 1.0 / reference == 0.0
        {
            return Err("invalid affine load-space dimension or compliance reference");
        }
        Ok(Self {
            dimension,
            reference,
            oracle_options: QuadraticBallOptions::default(),
        })
    }

    /// Replace the native bounded oracle policy after its own validation.
    pub fn with_oracle_options(
        mut self,
        options: QuadraticBallOptions,
    ) -> Result<Self, &'static str> {
        options
            .validate()
            .map_err(|_| "invalid affine load-space oracle policy")?;
        self.oracle_options = options;
        Ok(self)
    }
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }
    #[must_use]
    pub const fn reference(&self) -> f64 {
        self.reference
    }
    #[must_use]
    pub const fn oracle_options(&self) -> QuadraticBallOptions {
        self.oracle_options
    }
}

/// Three physical realizations from ONE forward filter at the reported raw
/// design. Only eroded material participates in the uncertainty solves.
#[derive(Debug, Clone)]
pub struct AffineLoadSpacePhysical3 {
    pub rho: Vec<f64>,
    pub eroded_rho: Vec<f64>,
    pub eroded_scales: Vec<f64>,
    pub eroded_volume_fraction: f64,
    pub nominal_rho: Vec<f64>,
    pub nominal_volume_fraction: f64,
    pub dilated_rho: Vec<f64>,
    pub dilated_volume_fraction: f64,
    pub dilated_volume_gradient: Vec<f64>,
    /// Requested by a mass-objective owner, using the SAME dilated design.
    pub dilated_mass_kg: Option<f64>,
    pub dilated_mass_gradient: Option<Vec<f64>>,
}

/// Complete affine uncertainty experiment on one eroded discrete operator.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceEvaluation3 {
    pub policy: AffineLoadSpace3,
    pub robustness: RobustOptions3,
    pub physical: AffineLoadSpacePhysical3,
    /// Symmetric row-major G = F^T K^-1 F / reference. Reciprocity is checked
    /// before averaging independently solved cross terms.
    pub gram: Vec<f64>,
    /// h = F^T K^-1 b / reference in a^T G a + 2*h^T a + c.
    pub linear: Vec<f64>,
    /// c = b^T K^-1 b / reference, including the native material body force.
    pub offset: f64,
    /// Native maximizing branch, witness and numerical optimality evidence.
    /// Its value/dual value exclude the constant `offset`.
    pub oracle: QuadraticBallMaximum,
    /// Complete normalized maximum, including the physical offset compliance.
    pub maximum: f64,
    /// Exact raw-design derivative only on a numerically separated regular
    /// branch. None does not invalidate the witness or assert a zero gradient.
    pub gradient: Option<Vec<f64>>,
    /// Largest admitted relative mismatch of G_ij/G_ji or F_i^T u0/b^T Ui.
    pub reciprocity_defect: f64,
    /// Relative discrepancy between the reduced quadratic value and actual
    /// force/displacement work of the returned combination.
    pub witness_relative_defect: f64,
    /// True equilibrium residual after combining all independently solved
    /// fields. Checked by one operator application, without another solve.
    pub witness_relative_residual: f64,
    /// Masked physical baseline b = b_dead + B*rho_eroded.
    pub offset_force: Vec<f64>,
    pub offset_displacement: Vec<f64>,
    /// Direct solutions Ui=K^-1 Fi; no subtraction of offset-loaded states.
    pub basis_displacements: Vec<Vec<f64>>,
    pub worst_coefficients: Vec<f64>,
    pub worst_force: Vec<f64>,
    pub worst_displacement: Vec<f64>,
    /// Actual physical work of the returned witness, before normalization.
    pub worst_compliance: f64,
    pub work: SolveWork,
}

fn finite(value: f64, stage: &'static str) -> Result<f64, EvaluationStop> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(failure(stage))
    }
}

// Keep the physical value until division so normalization cannot silently
// erase a nonzero work or sensitivity. Overflow and underflow both refuse.
fn normalize(value: f64, reference: f64, stage: &'static str) -> Result<f64, EvaluationStop> {
    let normalized = value / reference;
    if normalized.is_finite() && (value == 0.0 || normalized != 0.0) {
        Ok(normalized)
    } else {
        Err(failure(stage))
    }
}

// Raw basis forces are allowed to carry entries at strong clamps. Both the
// actual solves and every reported work/combination use the same masked space.
fn force_work(
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

fn relative_difference(a: f64, b: f64, scale: f64) -> f64 {
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

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate worst compliance for one fixed affine load ellipsoid, including
    /// an optional native material self-weight map. `dead_offset` must have
    /// full displacement-space length; supply zeros for no authored baseline.
    /// It is independent of design. Each basis column is another fixed force;
    /// amplitudes and correlations are encoded by F before this call.
    ///
    /// One eroded preparation serves the direct baseline and basis solves.
    /// Self-weight is assembled once at eroded physical occupancy and added
    /// only to the baseline. It is never replicated into the basis columns.
    /// The coefficient witness reconstructs the physical force/displacement
    /// without another elasticity solve. Regular-branch sensitivities include
    /// both stiffness and 2*B^T*u through the shared physical density pullback.
    ///
    /// Only homogeneous supports and fixed reference geometry/acceleration are
    /// covered. Eroded material is the declared manufacturing realization; with
    /// self-weight it is not asserted to be worst across threshold realizations.
    /// All outcomes restore incoming scales; nominal parameters never change.
    /// Numerical oracle evidence is not a floating-point or continuum bound.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_robust_affine_load_space(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        self.affine_load_space_with_mass(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            false,
            control,
        )
    }

    // The optimization owner requests kg from the already prepared dilated
    // field. No second forward filter, eta mutation or redundant equilibrium.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn affine_load_space_with_mass(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        include_mass: bool,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        let m = policy.dimension;
        let n = self.operator.n();
        if !robustness.valid(self.params)
            || policy.oracle_options.validate().is_err()
            || rho.len() != self.cells()
            || rho
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || dead_offset.len() != n
            || dead_offset.iter().any(|v| !v.is_finite())
            || basis.len() != m
            || basis
                .iter()
                .any(|f| f.len() != n || f.iter().any(|v| !v.is_finite()))
        {
            return Err(failure("sdf3-affine-load-space-input"));
        }
        if include_mass && self.operator.reference_cell_masses_kg().is_none() {
            return Err(failure("sdf3-material-mass-unavailable"));
        }
        control.checkpoint("sdf3-affine-load-space-start")?;
        let filtered = self.filter_apply(rho, control)?;
        let nominal = self.design_from_filtered(&filtered, self.params.eta, control)?;
        let dilated = self.design_from_filtered(
            &filtered,
            self.params.eta - robustness.threshold_offset,
            control,
        )?;
        let eroded = self.design_from_filtered(
            &filtered,
            self.params.eta + robustness.threshold_offset,
            control,
        )?;
        if eroded
            .projected
            .iter()
            .zip(&nominal.projected)
            .zip(&dilated.projected)
            .any(|((&e, &n), &d)| e > n + 32.0 * f64::EPSILON || n > d + 32.0 * f64::EPSILON)
            || [eroded.volume, nominal.volume, dilated.volume]
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(failure("sdf3-affine-load-space-ordering"));
        }
        let dilated_local: Vec<_> = dilated
            .slope
            .iter()
            .zip(&self.mass)
            .map(|(s, m)| s * m)
            .collect();
        let dilated_volume_gradient = self.pullback(&dilated_local, control)?;
        if dilated_volume_gradient.iter().any(|v| !v.is_finite()) {
            return Err(failure("sdf3-affine-load-space-volume-gradient"));
        }
        let mass = if include_mass {
            Some(self.mass_prepared(&dilated, control)?)
        } else {
            None
        };
        let previous = self.operator.scales().to_vec();
        self.operator
            .set_scales(&eroded.scales)
            .map_err(|_| failure("sdf3-affine-load-space-scales"))?;
        let result = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            let body = self.material_body_force(&eroded.projected, control)?;
            let offset_force = self.external_force(dead_offset, body.as_deref(), control)?;
            let offset_displacement = checked_solve_preconditioned(
                &self.operator,
                &prepared,
                &offset_force,
                1e-11,
                "sdf3-affine-load-space-offset",
                control,
            )?;
            let mut basis_displacements = Vec::with_capacity(m);
            for force in basis {
                let rhs = self.external_force(force, None, control)?;
                basis_displacements.push(checked_solve_preconditioned(
                    &self.operator,
                    &prepared,
                    &rhs,
                    1e-11,
                    "sdf3-affine-load-space-basis",
                    control,
                )?);
            }
            let fixed = self.operator.fixed();
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
            if (0..m).any(|i| {
                gram[i * m + i] < 0.0 || (gram[i * m + i] == 0.0 && has_force(basis[i], fixed))
            }) {
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
            let mut worst_force = offset_force.clone();
            let mut worst_displacement = offset_displacement.clone();
            for d in 0..n {
                if d % 256 == 0 {
                    control.checkpoint("sdf3-affine-load-space-witness")?;
                }
                if fixed[d / 3] {
                    worst_force[d] = 0.0;
                    worst_displacement[d] = 0.0;
                    continue;
                }
                for i in 0..m {
                    worst_force[d] += worst_coefficients[i] * basis[i][d];
                    worst_displacement[d] += worst_coefficients[i] * basis_displacements[i][d];
                }
                if !worst_force[d].is_finite() || !worst_displacement[d].is_finite() {
                    return Err(failure("sdf3-affine-load-space-witness"));
                }
            }
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
                witness_residual(&self.operator, &worst_force, &worst_displacement, control)?;
            let gradient = if oracle.gradient_admissible() {
                // Danskin's derivative holds for the resolved unique maximizer.
                // No derivative of its coefficients or the iterative solves.
                let mut local = crate::sdf3::body::controlled(control, |poll| {
                    self.operator.scale_bilinear_forms(
                        &worst_displacement,
                        &worst_displacement,
                        poll,
                    )
                })?;
                for value in &mut local {
                    *value = normalize(
                        -*value,
                        policy.reference,
                        "sdf3-affine-load-space-stiffness-gradient",
                    )?;
                }
                let occupancy = if let Some(load) = self.material_body_load() {
                    // Preserve the native unnormalized B^T*u until its own
                    // range check; multiplying by a reciprocal first can hide
                    // underflow before the factor of two is applied.
                    let mut work = crate::sdf3::body::controlled(control, |poll| {
                        load.occupancy_pullback(&worst_displacement, poll)
                    })?;
                    for value in &mut work {
                        *value = finite(
                            2.0 * normalize(
                                *value,
                                policy.reference,
                                "sdf3-affine-load-space-body-gradient",
                            )?,
                            "sdf3-affine-load-space-body-gradient",
                        )?;
                    }
                    work
                } else {
                    vec![0.0; self.cells()]
                };
                Some(self.physical_gradient(&eroded, &local, &occupancy, control)?)
            } else {
                None
            };
            let (dilated_mass_kg, dilated_mass_gradient) = match mass {
                Some((value, gradient)) => (Some(value), Some(gradient)),
                None => (None, None),
            };
            let physical = AffineLoadSpacePhysical3 {
                rho: rho.to_vec(),
                eroded_rho: eroded.projected,
                eroded_scales: eroded.scales,
                eroded_volume_fraction: eroded.volume,
                nominal_rho: nominal.projected,
                nominal_volume_fraction: nominal.volume,
                dilated_rho: dilated.projected,
                dilated_volume_fraction: dilated.volume,
                dilated_volume_gradient,
                dilated_mass_kg,
                dilated_mass_gradient,
            };
            control.checkpoint("sdf3-affine-load-space-publish")?;
            Ok(AffineLoadSpaceEvaluation3 {
                policy: policy.clone(),
                robustness,
                physical,
                gram,
                linear,
                offset,
                oracle,
                maximum,
                gradient,
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
                work: control.work(),
            })
        })();
        self.operator
            .set_scales(&previous)
            .expect("previous affine load-space scales remain valid");
        result
    }
}
