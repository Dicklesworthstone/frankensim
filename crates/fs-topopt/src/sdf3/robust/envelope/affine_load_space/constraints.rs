//! One constraint family over the shared material-aware affine primal fields.
use super::response::reference::{self as reference, ReferenceAffineResponseLimits3};
use super::response::{self, AffineResponseLimits3, AffineResponseOptions3};
use super::stress::{self, AffineStressLimit3};
use super::{AffineLoadSpace3, AffineLoadSpaceEvaluation3, RobustOptions3};
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{EvaluationStop, SolveControl};

/// Bound metadata retained by every audit, evaluation and refit. The response
/// report additionally retains the ordered original inequalities themselves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineLoadSpaceConstraintPolicy3 {
    pub responses: Option<AffineResponseOptions3>,
    pub stress: Option<AffineStressLimit3>,
    /// Inverse width of the unnormalized upper aggregate across the original
    /// compliance/volume constraint and the optional response/strength families.
    pub sharpness: f64,
}
impl Default for AffineLoadSpaceConstraintPolicy3 {
    fn default() -> Self {
        Self {
            responses: None,
            stress: None,
            sharpness: 40.0,
        }
    }
}
impl AffineLoadSpaceConstraintPolicy3 {
    #[must_use]
    pub fn has_limits(self) -> bool {
        self.responses.is_some() || self.stress.is_some()
    }

    pub(super) fn matches(self, evaluation: &AffineLoadSpaceEvaluation3) -> bool {
        self.matches_components(evaluation.responses.as_ref(), evaluation.stress.as_ref())
    }

    pub(super) fn matches_components(
        self,
        responses: Option<&response::AffineResponseEvaluation3>,
        stress: Option<&stress::AffineStressEvaluation3>,
    ) -> bool {
        self.responses == responses.map(|r| r.options) && self.stress == stress.map(|s| s.limit)
    }
}

/// Optional serviceability AND strength bounds, held fixed throughout a run.
/// The original design goal supplies the compliance or volume constraint.
/// Stress is the conservative sampled-peak measure over the complete load ball;
/// response laws retain their original signed physical inequalities.
#[derive(Clone, Copy)]
pub struct AffineLoadSpaceConstraints3<'a> {
    pub responses: Option<AffineResponseLimits3<'a>>,
    pub stress: Option<AffineStressLimit3>,
    pub sharpness: f64,
}
impl Default for AffineLoadSpaceConstraints3<'_> {
    fn default() -> Self {
        Self {
            responses: None,
            stress: None,
            sharpness: 40.0,
        }
    }
}
impl<'a> AffineLoadSpaceConstraints3<'a> {
    #[must_use]
    pub fn from_responses(responses: AffineResponseLimits3<'a>) -> Self {
        Self {
            responses: Some(responses),
            stress: None,
            sharpness: responses.options.sharpness,
        }
    }

    #[must_use]
    pub fn policy(self) -> AffineLoadSpaceConstraintPolicy3 {
        AffineLoadSpaceConstraintPolicy3 {
            responses: self.responses.map(|r| r.options),
            stress: self.stress,
            sharpness: self.sharpness,
        }
    }

    pub(super) fn admit(
        self,
        study: &CutDensityStudy3<impl Sdf3Elasticity>,
    ) -> Result<(), EvaluationStop> {
        admit_sharpness(self.sharpness)?;
        if let Some(responses) = self.responses {
            response::admit(study, responses)?;
        }
        if let Some(stress) = self.stress {
            stress::admit(study, stress)?;
        }
        Ok(())
    }
}

/// Grid-independent declarations for initial fitting and native adaptive
/// refits. Each grid binds its displacement observations afresh; stress always
/// samples that grid's retained quadrature and actual constitutive tensors.
#[derive(Clone, Copy)]
pub struct ReferenceAffineLoadSpaceConstraints3<'a> {
    pub responses: Option<ReferenceAffineResponseLimits3<'a>>,
    pub stress: Option<AffineStressLimit3>,
    pub sharpness: f64,
}
impl Default for ReferenceAffineLoadSpaceConstraints3<'_> {
    fn default() -> Self {
        Self {
            responses: None,
            stress: None,
            sharpness: 40.0,
        }
    }
}
impl<'a> ReferenceAffineLoadSpaceConstraints3<'a> {
    #[must_use]
    pub fn from_responses(responses: ReferenceAffineResponseLimits3<'a>) -> Self {
        Self {
            responses: Some(responses),
            stress: None,
            sharpness: responses.options.sharpness,
        }
    }

    #[must_use]
    pub fn policy(self) -> AffineLoadSpaceConstraintPolicy3 {
        AffineLoadSpaceConstraintPolicy3 {
            responses: self.responses.map(|r| r.options),
            stress: self.stress,
            sharpness: self.sharpness,
        }
    }

    pub(super) fn admit(
        self,
        study: &CutDensityStudy3<impl Sdf3Elasticity>,
    ) -> Result<(), EvaluationStop> {
        admit_sharpness(self.sharpness)?;
        if let Some(responses) = self.responses {
            reference::admit(study.operator(), responses)?;
        }
        if let Some(stress) = self.stress {
            stress::admit(study, stress)?;
        }
        Ok(())
    }
}

pub(super) fn admit_sharpness(sharpness: f64) -> Result<(), EvaluationStop> {
    if !sharpness.is_finite() || sharpness <= 0.0 || !(1.0 / sharpness).is_finite() {
        return Err(failure("sdf3-affine-constraint-policy"));
    }
    Ok(())
}

/// The SAME upper aggregate serves optimization, physical archives and grid
/// assessment. Weights are in [base, responses, stress] order, with zero for
/// absent components. No normalization by component count dilutes a violation.
pub(super) fn aggregate(
    base: f64,
    responses: Option<f64>,
    stress: Option<f64>,
    sharpness: f64,
) -> Result<(f64, [f64; 3]), EvaluationStop> {
    admit_sharpness(sharpness)?;
    if responses.is_none() && stress.is_none() {
        return super::grid::finite(base, "sdf3-affine-constraint-value")
            .map(|v| (v, [1.0, 0.0, 0.0]));
    }
    let mut values = vec![base];
    values.extend(responses);
    values.extend(stress);
    let (value, _, weights) = crate::sdf3::response::limits::aggregate(&values, sharpness)
        .map_err(|error| match error {
            crate::sdf3::response::ResponseError3::Evaluation(stop) => stop,
            _ => failure("sdf3-affine-design-joint-constraint"),
        })?;
    let mut result = [weights[0], 0.0, 0.0];
    let mut next = 1;
    if responses.is_some() {
        result[1] = weights[next];
        next += 1;
    }
    if stress.is_some() {
        result[2] = weights[next];
    }
    Ok((value, result))
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate compliance, original response limits and sampled strength over
    /// the SAME continuous load ball. A single physical map and stiffness
    /// preparation supply every component; each retains its own adverse loads.
    /// Native material-domain admission covers the complete strain ball.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_robust_affine_constraints(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        constraints: AffineLoadSpaceConstraints3<'_>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        self.affine_load_space_constraints_with_mass(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            constraints,
            false,
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn affine_load_space_constraints_with_mass(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        constraints: AffineLoadSpaceConstraints3<'_>,
        include_mass: bool,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        control.checkpoint("sdf3-affine-constraints-start")?;
        constraints.admit(self)?;
        if constraints.stress.is_some_and(|limit| {
            basis
                .len()
                .checked_add(1)
                .is_none_or(|n| n > limit.options.max_cases)
        }) {
            return Err(failure("sdf3-affine-stress-field-allowance"));
        }
        let (mut result, (responses, stress)) = self.affine_load_space_extended(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            include_mass,
            control,
            |study, design, prepared, grid, control| {
                let responses = constraints
                    .responses
                    .map(|family| {
                        response::evaluate(study, design, prepared, grid, family, control)
                    })
                    .transpose()?;
                let stress = constraints
                    .stress
                    .map(|limit| {
                        stress::evaluate_prepared(
                            study,
                            design,
                            prepared,
                            grid,
                            limit,
                            policy.oracle_options(),
                            control,
                        )
                    })
                    .transpose()?;
                Ok((responses, stress))
            },
        )?;
        control.checkpoint("sdf3-affine-constraints-publish")?;
        result.constraints = constraints.policy();
        result.responses = responses;
        result.stress = stress;
        result.work = control.work();
        if let Some(response) = &mut result.responses {
            response.work = control.work();
        }
        if let Some(stress) = &mut result.stress {
            stress.work = control.work();
        }
        Ok(result)
    }
}
