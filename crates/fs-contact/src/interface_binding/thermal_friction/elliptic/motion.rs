//! Actual-motion closure for the original directional coefficient sources.
//!
//! The matrix solve fixes slip for each queried distance. This owner closes
//! that scalar distance against the norm of the returned physical slip. Its
//! explicit bracket belongs to the caller's mechanical/source domain; general
//! matrix recoil does not bound slip by the norm of released displacement.

use core::cell::Cell;

use super::*;

/// Finite source-motion domain and accuracy/work budgets for one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionQuerySolve {
    /// Nonnegative queried slip lengths [m], with every endpoint inside the
    /// admitted source domain. Different brackets may select different roots.
    pub distance_bracket_m: [f64; 2],
    pub absolute_distance_tolerance_m: f64,
    pub relative_distance_tolerance: f64,
    /// Complete original-source/native-response evaluations, including endpoints.
    pub maximum_evaluations: usize,
}

/// One source and mechanical candidate; neither original owner is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct MotionQueryReport {
    pub step: TangentialStep<SourceThermalEllipticState>,
    /// Distance used for every owned speed and accumulated-path coordinate.
    /// Equals actual distance when the source declares no motion coordinates.
    pub queried_distance_m: f64,
    pub actual_distance_m: f64,
    /// Queried minus actual distance. Accepted motion-dependent roots are
    /// nonpositive, so the source path cannot run ahead of resolved motion.
    pub distance_residual_m: f64,
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MotionQueryError {
    Source(SourceThermalCoulombError),
    Invalid(&'static str),
    Cancelled,
    EvaluationBudget {
        evaluations: usize,
        maximum: usize,
    },
    DistanceOutsideBracket {
        distance_m: f64,
        bracket_m: [f64; 2],
    },
    Unbracketed {
        bracket_m: [f64; 2],
        residual_m: [f64; 2],
    },
    Unresolved {
        bracket_m: [f64; 2],
        residual_m: [f64; 2],
    },
}

impl core::fmt::Display for MotionQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional motion closure: {self:?}")
    }
}
impl std::error::Error for MotionQueryError {}
impl From<SourceThermalCoulombError> for MotionQueryError {
    fn from(error: SourceThermalCoulombError) -> Self {
        Self::Source(error)
    }
}

pub(super) fn finite(value: f64) -> Result<f64, MotionQueryError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MotionQueryError::Invalid(
            "unrepresentable source-motion coordinate",
        ))
    }
}

fn distance(value: [f64; 2]) -> Result<f64, MotionQueryError> {
    finite(fs_math::det::hypot(value[0], value[1]))
}

pub(super) fn validate_controls(
    controls: MotionQuerySolve,
    dt_s: f64,
) -> Result<(), MotionQueryError> {
    let [lower, upper] = controls.distance_bracket_m;
    if !lower.is_finite()
        || lower < 0.0
        || !upper.is_finite()
        || upper < lower
        || !controls.absolute_distance_tolerance_m.is_finite()
        || controls.absolute_distance_tolerance_m < 0.0
        || !controls.relative_distance_tolerance.is_finite()
        || controls.relative_distance_tolerance < 0.0
        || controls.absolute_distance_tolerance_m == 0.0
            && controls.relative_distance_tolerance == 0.0
        || controls.maximum_evaluations == 0
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(MotionQueryError::Invalid(
            "finite source-motion bracket, tolerances and duration required",
        ));
    }
    Ok(())
}

pub(super) fn query_for_distance(
    template: &QueryPoint,
    accepted: &QueryPoint,
    axes: &[FrictionQueryAxis],
    r: f64,
    dt_s: f64,
) -> Result<QueryPoint, MotionQueryError> {
    let mut point = template.clone();
    for axis in axes {
        let value = match axis.kind {
            FrictionQueryAxisKind::SlidingSpeed => {
                let speed = finite(r / dt_s)?;
                if r > 0.0 && speed == 0.0 {
                    return Err(MotionQueryError::Invalid(
                        "positive source sliding speed is unrepresentable",
                    ));
                }
                speed
            }
            FrictionQueryAxisKind::AccumulatedDistance => {
                let old = accepted.axes()[&axis.name];
                let next = finite(old + r)?;
                if r > 0.0 && next <= old {
                    return Err(MotionQueryError::Invalid(
                        "source path increment is unrepresentable at its accepted origin",
                    ));
                }
                next
            }
            _ => continue,
        };
        point = point
            .with(&axis.name, value)
            .map_err(SourceThermalCoulombError::from)?;
    }
    Ok(point)
}

pub(super) fn admissible_distance(
    controls: MotionQuerySolve,
    queried: f64,
    actual: f64,
) -> Result<bool, MotionQueryError> {
    let tolerance = finite(
        controls.absolute_distance_tolerance_m
            + controls.relative_distance_tolerance * queried.max(actual),
    )?;
    let residual = finite(queried - actual)?;
    Ok(residual <= 0.0 && -residual <= tolerance)
}

/// Shared lower-sided scalar closure. The physical owner supplies the actual
/// distance: a held-step norm or a continuous trajectory's unsigned path.
pub(super) fn solve_distance<T, E: From<MotionQueryError>>(
    controls: MotionQuerySolve,
    mut evaluate: impl FnMut(f64) -> Result<(T, f64), E>,
) -> Result<T, E> {
    let [mut lower, mut upper] = controls.distance_bracket_m;
    let (candidate, actual) = evaluate(lower)?;
    if admissible_distance(controls, lower, actual)? {
        return Ok(candidate);
    }
    let mut low_residual = finite(lower - actual)?;
    drop(candidate);
    if lower == upper {
        return Err(MotionQueryError::Unbracketed {
            bracket_m: [lower, upper],
            residual_m: [low_residual; 2],
        }
        .into());
    }
    let (candidate, actual) = evaluate(upper)?;
    if admissible_distance(controls, upper, actual)? {
        return Ok(candidate);
    }
    let mut high_residual = finite(upper - actual)?;
    drop(candidate);
    if low_residual.is_sign_positive() == high_residual.is_sign_positive() {
        return Err(MotionQueryError::Unbracketed {
            bracket_m: [lower, upper],
            residual_m: [low_residual, high_residual],
        }
        .into());
    }
    loop {
        let width = upper - lower;
        let fraction = low_residual / (low_residual - high_residual);
        let fraction = if fraction.is_finite() && fraction > 0.0625 && fraction < 0.9375 {
            fraction
        } else {
            0.5
        };
        let middle = lower + fraction * width;
        if middle <= lower || middle >= upper {
            return Err(MotionQueryError::Unresolved {
                bracket_m: [lower, upper],
                residual_m: [low_residual, high_residual],
            }
            .into());
        }
        let (candidate, actual) = evaluate(middle)?;
        if admissible_distance(controls, middle, actual)? {
            return Ok(candidate);
        }
        let residual = finite(middle - actual)?;
        if residual.is_sign_positive() == low_residual.is_sign_positive() {
            lower = middle;
            low_residual = residual;
        } else {
            upper = middle;
            high_residual = residual;
        }
    }
}

impl SourceThermalEllipticCoulomb {
    /// Resolve speed/path-dependent original claims against actual matrix slip.
    ///
    /// `thermal_query` supplies the complete coordinate schema, temperatures,
    /// pressure and fixed external conditions. This operation replaces every
    /// explicitly owned sliding-speed coordinate with `r / dt` and every owned
    /// accumulated-path coordinate with its own accepted source value plus `r`.
    /// Caller values in those motion coordinates are placeholders, not state.
    ///
    /// A safeguarded scalar solve requires a sign bracket and an actual distance
    /// residual; a small bracket alone never grants a candidate. Acceptance is
    /// from `r <= norm(slip)` within the declared absolute/relative tolerance.
    /// Discontinuous or unresolved responses refuse. Arbitrary velocity/path
    /// dependence does not establish existence, uniqueness or a lowest root.
    /// The caller still owns physical branch selection and arrest localization.
    ///
    /// Exactly prescribed slip and motion-independent sources need one query.
    /// Every source/native evaluation is budgeted and cancellation-polled; only
    /// this candidate carries refreshed receipts, history and source path.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_motion_consistent(
        &self,
        state: &SourceThermalEllipticState,
        input: &TangentialStepInput,
        thermal_query: &QueryPoint,
        mode: EllipticFrictionMode,
        compliance: SymmetricTangentialCompliance,
        controls: MotionQuerySolve,
        mut cancel: impl FnMut() -> bool,
    ) -> Result<MotionQueryReport, MotionQueryError> {
        if cancel() {
            return Err(MotionQueryError::Cancelled);
        }
        self.check(state)?;
        validate_controls(controls, input.dt_s)?;
        let [lower, upper] = controls.distance_bracket_m;
        if input.compliance_m_per_n != 0.0 {
            return Err(MotionQueryError::Invalid("sole matrix compliance required"));
        }
        let axes = self.principal[0].query_axes();
        let has_motion = axes.iter().any(|axis| {
            matches!(
                axis.kind,
                FrictionQueryAxisKind::SlidingSpeed | FrictionQueryAxisKind::AccumulatedDistance
            )
        });
        let accepted = state.sources[0].query();
        let evaluations = Cell::new(0usize);
        let mut evaluate = |r: f64| -> Result<MotionQueryReport, MotionQueryError> {
            if cancel() {
                return Err(MotionQueryError::Cancelled);
            }
            if evaluations.get() >= controls.maximum_evaluations {
                return Err(MotionQueryError::EvaluationBudget {
                    evaluations: evaluations.get(),
                    maximum: controls.maximum_evaluations,
                });
            }
            evaluations.set(evaluations.get() + 1);
            let point = query_for_distance(thermal_query, &accepted, axes, r, input.dt_s)?;
            let step =
                self.advance_at_query_with_compliance(state, input, &point, mode, compliance)?;
            if cancel() {
                return Err(MotionQueryError::Cancelled);
            }
            let actual = distance(step.slip_increment_m)?;
            Ok(MotionQueryReport {
                step,
                queried_distance_m: r,
                actual_distance_m: actual,
                distance_residual_m: finite(r - actual)?,
                evaluations: evaluations.get(),
            })
        };
        let in_bracket = |r: f64| -> Result<(), MotionQueryError> {
            if r < lower || r > upper {
                Err(MotionQueryError::DistanceOutsideBracket {
                    distance_m: r,
                    bracket_m: [lower, upper],
                })
            } else {
                Ok(())
            }
        };
        if !has_motion {
            let mut candidate = evaluate(0.0)?;
            in_bracket(candidate.actual_distance_m)?;
            candidate.queried_distance_m = candidate.actual_distance_m;
            candidate.distance_residual_m = 0.0;
            return Ok(candidate);
        }
        if compliance.matrix_m_per_n() == [[0.0; 2]; 2]
            || input.normal_n == 0.0
            || input.slip_increment_m == [0.0; 2]
        {
            let r = distance(input.slip_increment_m)?;
            in_bracket(r)?;
            let candidate = evaluate(r)?;
            if admissible_distance(controls, r, candidate.actual_distance_m)? {
                return Ok(candidate);
            }
            return Err(MotionQueryError::Unresolved {
                bracket_m: [r; 2],
                residual_m: [candidate.distance_residual_m; 2],
            });
        }
        solve_distance(controls, |r| {
            let candidate = evaluate(r)?;
            let actual = candidate.actual_distance_m;
            Ok((candidate, actual))
        })
    }
}
