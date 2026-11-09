//! Fixed free-strain experiments through the original CutFEM density pipeline.
//!
//! Each independent experiment solves `K(s) u = f + B rho_bar + b_free(s)`.
//! Its equivalent-load work is the COMPLETE right-hand-side work, not elastic
//! storage or externally supplied heat. Its derivative includes
//! `2 u^T db_free/ds - u^T dK/ds u`; reference-density body force retains its
//! separate occupancy derivative. Neither free strain nor temperature is a
//! density variable here. Source-temperature, changing material tensors, heat
//! feedback and continuum load balls require their own physical owners.
//!
//! Cases borrow preparations on this exact operator. They are passed per call,
//! so cold modal, refinement and load-ball consumers cannot silently inherit a
//! thermal state. Geometry, material tensors and the original source field stay
//! fixed while the existing filter/projection/SIMP chain changes its scales.

use super::operator::native_controlled;
use super::{CutDensityEvaluation3, CutDensityStudy3, Sdf3Elasticity, failure};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, SolveControl};
use fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3;

/// One independently solved fixed-field thermal/manufacturing experiment.
#[derive(Debug, Clone, Copy)]
pub struct CutEigenstrainCase3<'a> {
    /// External force in the original independent displacement coordinates.
    /// Configured source-density body force is added by the original study.
    pub force: &'a [f64],
    /// Nonnegative compliance weight. Zero-weight cases still undergo physical
    /// admission; they cannot hide a material-domain failure.
    pub weight: f64,
    /// Original physical free strain, with source/geometry binding retained.
    pub eigenstrain: &'a PreparedEigenstrain3,
}

/// Explicit limits on the physical state and retained thermal observations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutEigenstrainPolicy3 {
    /// Frobenius ceiling in (0,1) for actual and B-bar TOTAL strain.
    /// Original constitutive limits independently bound elastic strain.
    pub maximum_kinematic_strain: f64,
    /// Maximum original bulk points in EACH independent experiment.
    pub max_points_per_case: usize,
}

#[derive(Clone, Copy)]
pub(super) struct EigenstrainFamily3<'a, 'field> {
    pub(super) cases: &'a [CutEigenstrainCase3<'field>],
    pub(super) policy: CutEigenstrainPolicy3,
}

impl<'a, 'field> EigenstrainFamily3<'a, 'field> {
    pub(super) fn admit(
        operator: &impl Sdf3Elasticity,
        cases: &'a [CutEigenstrainCase3<'field>],
        policy: CutEigenstrainPolicy3,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, EvaluationStop> {
        control.checkpoint("sdf3-eigenstrain-admission")?;
        if !policy.maximum_kinematic_strain.is_finite()
            || policy.maximum_kinematic_strain <= 0.0
            || policy.maximum_kinematic_strain >= 1.0
            || policy.max_points_per_case == 0
            || cases.is_empty()
            || !cases.iter().any(|case| case.weight > 0.0)
        {
            return Err(failure("sdf3-eigenstrain-policy"));
        }
        for case in cases {
            control.checkpoint("sdf3-eigenstrain-admission")?;
            if !case.weight.is_finite()
                || case.weight < 0.0
                || case.force.len() != operator.n()
                || case.eigenstrain.bulk_points() > policy.max_points_per_case
            {
                return Err(failure("sdf3-eigenstrain-case"));
            }
            for chunk in case.force.chunks(192) {
                control.checkpoint("sdf3-eigenstrain-admission")?;
                if chunk.iter().any(|force| !force.is_finite()) {
                    return Err(failure("sdf3-eigenstrain-case"));
                }
            }
            operator
                .validate_eigenstrain(case.eigenstrain)
                .map_err(|_| failure("sdf3-eigenstrain-owner"))?;
        }
        Ok(Self { cases, policy })
    }

    pub(super) fn loads(self) -> Vec<LoadCase<'field>> {
        self.cases
            .iter()
            .map(|case| LoadCase {
                force: case.force,
                weight: case.weight,
            })
            .collect()
    }

    pub(super) fn add_load(
        self,
        operator: &impl Sdf3Elasticity,
        case: usize,
        rhs: &mut [f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), EvaluationStop> {
        let load = native_controlled(control, "sdf3-eigenstrain-load", |poll| {
            operator.eigenstrain_load(self.cases[case].eigenstrain, poll)
        })?;
        for (i, (rhs, free)) in rhs.iter_mut().zip(load).enumerate() {
            if i % 192 == 0 {
                control.checkpoint("sdf3-eigenstrain-load")?;
            }
            if free != 0.0 {
                *rhs += free;
                if !rhs.is_finite() {
                    return Err(failure("sdf3-eigenstrain-load-overflow"));
                }
            }
        }
        Ok(())
    }

    pub(super) fn admit_displacement(
        self,
        operator: &impl Sdf3Elasticity,
        case: usize,
        displacement: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), EvaluationStop> {
        native_controlled(control, "sdf3-eigenstrain-domain", |poll| {
            operator.admit_eigenstrain(
                self.cases[case].eigenstrain,
                displacement,
                self.policy.maximum_kinematic_strain,
                self.policy.max_points_per_case,
                poll,
            )
        })
    }

    pub(super) fn add_scale_gradient(
        self,
        operator: &impl Sdf3Elasticity,
        case: usize,
        dual: &[f64],
        multiplier: f64,
        gradient: &mut [f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), EvaluationStop> {
        let work = native_controlled(control, "sdf3-eigenstrain-pullback", |poll| {
            operator.eigenstrain_load_scale_pullback(self.cases[case].eigenstrain, dual, poll)
        })?;
        for (i, (value, work)) in gradient.iter_mut().zip(work).enumerate() {
            if i % 64 == 0 {
                control.checkpoint("sdf3-eigenstrain-pullback")?;
            }
            if work != 0.0 && multiplier != 0.0 {
                *value += multiplier * work;
                if !value.is_finite() {
                    return Err(failure("sdf3-eigenstrain-pullback-overflow"));
                }
            }
        }
        Ok(())
    }
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate fixed thermal experiments on the existing eroded realization
    /// with the unchanged nominal/dilated material-volume constraints.
    /// Thermal load scales follow the SAME eroded SIMP material as stiffness;
    /// material-density force follows its eroded physical occupancy. This does
    /// not assert that eroded thermal work is the worst across realizations.
    /// Success and failure both restore incoming scales, like the cold robust
    /// owner. No continuous load ball or automatic geometry transfer is implied.
    pub fn evaluate_eigenstrain_robust(
        &mut self,
        rho: &[f64],
        cases: &[CutEigenstrainCase3<'_>],
        policy: CutEigenstrainPolicy3,
        robustness: super::robust::RobustOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<super::robust::RobustEvaluation3, EvaluationStop> {
        let family = EigenstrainFamily3::admit(&self.operator, cases, policy, control)?;
        if !robustness.valid(self.params)
            || rho.len() != self.cells()
            || rho
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
        {
            return Err(failure("sdf3-eigenstrain-robust-policy"));
        }
        let loads = family.loads();
        self.evaluate_robust_family(rho, &loads, Some(family), robustness, control)
    }

    /// Differentiate the existing relaxed/reference or sampled-peak stress
    /// functional under independent fixed free-strain experiments.
    ///
    /// Reference stress is `C_reference:(epsilon - epsilon_free)`, with the
    /// native B-bar projection when selected. The existing stress owner retains
    /// its direct relaxation/peak derivative and reference-stress adjoint; the
    /// thermal load contributes `z^T db_free/dscale` through equilibrium.
    /// `case_compliances` retains that owner's EXTERNAL-force work convention,
    /// excluding the equivalent free-strain load. Use [`Self::evaluate_eigenstrain`]
    /// when complete equivalent-load work is the desired objective.
    ///
    /// Both success and failure restore incoming scales. Admission covers all
    /// retained physical bulk points, including zero-weight cases, without
    /// asserting a between-point bound for the original sampled free strain.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_eigenstrain_stress(
        &mut self,
        rho: &[f64],
        cases: &[CutEigenstrainCase3<'_>],
        policy: CutEigenstrainPolicy3,
        options: super::stress::StressOptions3,
        measure: super::stress::StressMeasure3,
        control: &mut SolveControl<'_>,
    ) -> Result<super::stress::StressEvaluation3, super::stress::StressError3> {
        let family = EigenstrainFamily3::admit(&self.operator, cases, policy, control)?;
        let loads = family.loads();
        super::stress::admit(self, rho, &loads, options)?;
        let points = cases
            .iter()
            .try_fold(0usize, |sum, case| {
                sum.checked_add(case.eigenstrain.bulk_points())
            })
            .ok_or(super::stress::StressError3::PointBudget)?;
        if points > options.max_points {
            return Err(super::stress::StressError3::PointBudget);
        }
        self.evaluate_stress_family(rho, &loads, None, Some(family), options, measure, control)
    }

    /// Solve independent thermal experiments and differentiate their complete
    /// equivalent-load work through the original density map.
    ///
    /// Both the external and thermal right-hand side participate in the work
    /// objective. This is not external-force-only compliance, physical elastic
    /// storage, or a heat balance. Positive density gradients are possible;
    /// the force-only OC update must not be used for this functional.
    ///
    /// Success installs this design's scales, like [`Self::evaluate`]. Every
    /// returned failure restores incoming scales. Existing body-load source
    /// identity, preconditioner reuse, work accounting and true residual gates
    /// remain active. A rebuilt/refined operator requires a new preparation.
    pub fn evaluate_eigenstrain(
        &mut self,
        rho: &[f64],
        cases: &[CutEigenstrainCase3<'_>],
        policy: CutEigenstrainPolicy3,
        control: &mut SolveControl<'_>,
    ) -> Result<CutDensityEvaluation3, EvaluationStop> {
        let family = EigenstrainFamily3::admit(&self.operator, cases, policy, control)?;
        if rho.len() != self.cells()
            || rho
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
        {
            return Err(failure("sdf3-eigenstrain-design"));
        }
        let loads = family.loads();
        let design = self.design(rho, control)?;
        self.evaluate_design_family(design, &loads, Some(family), control)
    }
}
