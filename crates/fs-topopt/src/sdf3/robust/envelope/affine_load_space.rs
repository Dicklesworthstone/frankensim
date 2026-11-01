//! Continuous affine uncertainty: f = b_dead + B*rho_eroded + F*a, ||a|| <= 1.
//!
//! The native bounded quadratic-ball owner supplies the global load witness.
//! This numerically solves the exact discrete quadratic formulation, with
//! numerical optimality evidence, rather than a smooth spectral surrogate.
//! A selected hard-case witness is useful evidence but is not a differentiable
//! objective. Only a numerically separated regular maximizer gets a gradient.
pub mod campaign;
pub mod constraints;
pub mod design;
mod grid;
pub mod refinement;
pub mod refit;
pub mod response;
pub mod stress;
pub use grid::AffineLoadSpaceGrid3;

use super::super::RobustOptions3;
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_ascent::quadratic_ball::{
    MAX_QUADRATIC_BALL_DIMENSION, QuadraticBallMaximum, QuadraticBallOptions,
};
use grid::{finite, normalize};

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
    /// Optional original response limits over the SAME complete affine ball.
    /// Each limit has its own adverse direction; compliance's witness is not
    /// substituted for a displacement or support-force extremum.
    pub responses: Option<response::AffineResponseEvaluation3>,
    /// Original optional constraint policy, shared with the audit and optimizer.
    pub constraints: constraints::AffineLoadSpaceConstraintPolicy3,
    /// Full-ball sampled strength, using the same native material and fields.
    pub stress: Option<stress::AffineStressEvaluation3>,
    pub work: SolveWork,
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
    /// The complete physical strain ball must satisfy every retained native
    /// material domain, even when no stress constraint is requested.
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
        self.affine_load_space_extended(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            include_mass,
            control,
            |_, _, _, _, _| Ok(()),
        )
        .map(|(evaluation, ())| evaluation)
    }

    // A component borrows the actual eroded design, native preparation and
    // freshly solved affine field family inside this one material transaction.
    // Nothing public admits cached fields, and the component cannot mutate the
    // operator or trigger a second filter/preconditioner preparation.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn affine_load_space_extended<T>(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        include_mass: bool,
        control: &mut SolveControl<'_>,
        extension: impl FnOnce(
            &Self,
            &crate::sdf3::Design3,
            &O::Prepared<'_>,
            &AffineLoadSpaceGrid3,
            &mut SolveControl<'_>,
        ) -> Result<T, EvaluationStop>,
    ) -> Result<(AffineLoadSpaceEvaluation3, T), EvaluationStop> {
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
            let grid = grid::evaluate_grid(
                &self.operator,
                &prepared,
                offset_force,
                basis,
                policy,
                control,
            )?;
            let extra = extension(self, &eroded, &prepared, &grid, control)?;
            let AffineLoadSpaceGrid3 {
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
            } = grid;
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
            Ok((
                AffineLoadSpaceEvaluation3 {
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
                    responses: None,
                    constraints: constraints::AffineLoadSpaceConstraintPolicy3::default(),
                    stress: None,
                    work: control.work(),
                },
                extra,
            ))
        })();
        self.operator
            .set_scales(&previous)
            .expect("previous affine load-space scales remain valid");
        result
    }
}
