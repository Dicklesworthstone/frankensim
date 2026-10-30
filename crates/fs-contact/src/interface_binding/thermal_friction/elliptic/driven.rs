//! Original source claims and accepted contact mode under external driving.
//!
//! The native driven owner supplies static hold, breakaway, kinetic motion and
//! explicit drive work. This adapter closes every declared motion coordinate
//! against that owner's unsigned travel, retaining original material pins and
//! the caller's accepted Resting/Sliding mode through every private source trial.

use super::motion::{
    MotionQueryError, MotionQuerySolve, finite, query_for_distance, solve_distance,
    validate_controls,
};
use super::*;
use fs_tribo::tangential::elliptic::driven::{
    EllipticDrivenControls, EllipticDrivenError, EllipticDrivenInput, EllipticDrivenReport,
    EllipticDrivenSegment, EllipticDrivenSegmentError,
};
use std::convert::Infallible;

#[derive(Clone, Debug)]
pub struct SourceDrivenReport {
    pub driven: EllipticDrivenReport,
    pub state: SourceThermalEllipticState,
    pub queried_distance_m: f64,
    pub actual_distance_m: f64,
    /// Query minus actual unsigned path; accepted motion roots are nonpositive.
    pub distance_residual_m: f64,
    pub evaluations: usize,
}

#[derive(Debug)]
pub enum SourceDrivenError {
    Query(MotionQueryError),
    Driven(EllipticDrivenError),
}
impl From<MotionQueryError> for SourceDrivenError {
    fn from(e: MotionQueryError) -> Self {
        Self::Query(e)
    }
}
impl From<SourceThermalCoulombError> for SourceDrivenError {
    fn from(e: SourceThermalCoulombError) -> Self {
        Self::Query(MotionQueryError::Source(e))
    }
}
impl core::fmt::Display for SourceDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source driven directional friction: {self:?}")
    }
}
impl std::error::Error for SourceDrivenError {}

/// Source trials and physical-segment replay retain distinct failure owners.
#[derive(Debug)]
pub enum SourceDrivenSegmentError<E> {
    Source(SourceDrivenError),
    Segments(EllipticDrivenSegmentError<E>),
    /// The accepted source trial and identical native replay must agree.
    ReplayMismatch,
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceDrivenSegmentError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source driven physical segments: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceDrivenSegmentError<E> {}

impl SourceThermalEllipticCoulomb {
    /// Resolve the original source pair against driven unsigned travel.
    ///
    /// The mechanical owner supplies the full tangent mobility, relative drive
    /// acceleration and accepted contact mode. Every trial begins at that same
    /// mode/history; a trial cannot turn numerical rest noise into breakaway.
    /// Temperatures, pressure and fixed conditions come from the complete query.
    /// Each declared speed becomes `r/dt`; each owned path starts at its own
    /// accepted coordinate and advances by `r`. Original claims never switch.
    ///
    /// The shared safeguarded distance root accepts only an actual residual
    /// from below. The caller supplies a bracket within the physical/source
    /// domain; driving and full-matrix recoil provide no free-slip norm bound.
    /// Source/native failures and cancellation preserve all borrowed owners.
    ///
    /// Coefficients are held within this source phase. The native curved driven
    /// trajectory has its own timestep accuracy controls; changing source
    /// temperature/speed/path requires separate outer phase refinement. Neither
    /// a source root nor an energy balance certifies a unique continuous solution.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_driven_motion_consistent(
        &self,
        state: &SourceThermalEllipticState,
        input: &EllipticDrivenInput,
        thermal_query: &QueryPoint,
        motion: MotionQuerySolve,
        driven: EllipticDrivenControls,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<SourceDrivenReport, SourceDrivenError> {
        self.resolve_driven_motion(state, input, thermal_query, motion, driven, false, cancel)
    }

    /// Resolve source motion and visit the accepted held-traction trajectory.
    ///
    /// Source-distance trials use the SAME segmented native method as replay,
    /// including its bounded zero-drive curve. They never call `visit`. After
    /// accepting a root, the selected coefficients are replayed from the original
    /// pre-phase history without re-querying material claims. The complete replay
    /// report must equal the accepted native trial before this method succeeds.
    /// The returned evaluation count counts source trials; replay performs one
    /// additional bounded native integration, with no additional source query.
    ///
    /// Each visitor receives actual accepted duration, force and motion. A pose
    /// owner must apply these segments individually; averaged phase force does
    /// not preserve their path or work. Visitors MUST mutate private staged state:
    /// a later segment refusal or cancellation cannot undo previous callbacks.
    /// Native/source histories remain borrowed and unchanged on every failure.
    ///
    /// Coefficients are held over this phase. Plane/body geometry and source
    /// sampling refresh require their own outer phase control. General curved
    /// segments are first-order even at zero drive; no exact curved flow is
    /// claimed by this segment-consumer entry point.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_driven_motion_consistent_with_segments<E>(
        &self,
        state: &SourceThermalEllipticState,
        input: &EllipticDrivenInput,
        thermal_query: &QueryPoint,
        motion: MotionQuerySolve,
        driven: EllipticDrivenControls,
        cancel: &mut impl FnMut() -> bool,
        visit: &mut impl FnMut(&EllipticDrivenSegment) -> Result<(), E>,
    ) -> Result<SourceDrivenReport, SourceDrivenSegmentError<E>> {
        let candidate = self
            .resolve_driven_motion(state, input, thermal_query, motion, driven, true, cancel)
            .map_err(SourceDrivenSegmentError::Source)?;
        let (law, original_history) = state
            .law
            .continue_coefficients(
                &state.constitutive,
                candidate.state.law.static_coefficients(),
                candidate.state.law.kinetic_coefficients(),
            )
            .map_err(|e| {
                SourceDrivenSegmentError::Source(SourceDrivenError::Driven(
                    EllipticDrivenError::Constitutive(e),
                ))
            })?;
        let replay = law
            .advance_driven_segments(&original_history, input, driven, cancel, visit)
            .map_err(SourceDrivenSegmentError::Segments)?;
        if replay != candidate.driven {
            return Err(SourceDrivenSegmentError::ReplayMismatch);
        }
        Ok(candidate)
    }

    /// Shared original source selection and safeguarded unsigned-distance root.
    #[allow(clippy::too_many_arguments)]
    fn resolve_driven_motion(
        &self,
        state: &SourceThermalEllipticState,
        input: &EllipticDrivenInput,
        thermal_query: &QueryPoint,
        motion: MotionQuerySolve,
        driven: EllipticDrivenControls,
        segmented: bool,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<SourceDrivenReport, SourceDrivenError> {
        if cancel() {
            return Err(MotionQueryError::Cancelled.into());
        }
        self.check(state)?;
        validate_controls(motion, input.dt_s)?;
        let axes = self.principal[0].query_axes();
        let has_motion = axes.iter().any(|axis| {
            matches!(
                axis.kind,
                FrictionQueryAxisKind::SlidingSpeed | FrictionQueryAxisKind::AccumulatedDistance
            )
        });
        let accepted = state.sources[0].query();
        let mut evaluations = 0;
        let mut evaluate = |r| -> Result<SourceDrivenReport, SourceDrivenError> {
            if cancel() {
                return Err(MotionQueryError::Cancelled.into());
            }
            if evaluations >= motion.maximum_evaluations {
                return Err(MotionQueryError::EvaluationBudget {
                    evaluations,
                    maximum: motion.maximum_evaluations,
                }
                .into());
            }
            evaluations += 1;
            let point = query_for_distance(thermal_query, &accepted, axes, r, input.dt_s)?;
            let sources = [0, 1].map(|i| {
                self.principal[i].source_at_query_conditions(
                    &state.sources[i],
                    input.contact_temperature_k,
                    input.nominal_pressure_pa,
                    &point,
                )
            });
            let [first, second] = sources;
            let sources = [first?, second?];
            let (static_mu, kinetic_mu) =
                principal_coefficients(&self.principal, [&sources[0], &sources[1]])?;
            let (law, history) = state
                .law
                .continue_coefficients(&state.constitutive, static_mu, kinetic_mu)
                .map_err(SourceThermalCoulombError::from)?;
            let result = if segmented {
                law.advance_driven_segments(&history, input, driven, cancel, &mut |_| {
                    Ok::<_, Infallible>(())
                })
                .map_err(|error| match error {
                    EllipticDrivenSegmentError::Driven(e) => SourceDrivenError::Driven(e),
                    EllipticDrivenSegmentError::Visitor(never) => match never {},
                })?
            } else {
                law.advance_driven(&history, input, driven, cancel)
                    .map_err(SourceDrivenError::Driven)?
            };
            if cancel() {
                return Err(MotionQueryError::Cancelled.into());
            }
            let actual = result.distance_m;
            Ok(SourceDrivenReport {
                state: SourceThermalEllipticState {
                    owner: self.identity,
                    law,
                    constitutive: result.state.clone(),
                    sources,
                },
                driven: result,
                queried_distance_m: r,
                actual_distance_m: actual,
                distance_residual_m: finite(r - actual)?,
                evaluations,
            })
        };
        if !has_motion {
            let mut candidate = evaluate(0.0)?;
            let [lower, upper] = motion.distance_bracket_m;
            if candidate.actual_distance_m < lower || candidate.actual_distance_m > upper {
                return Err(MotionQueryError::DistanceOutsideBracket {
                    distance_m: candidate.actual_distance_m,
                    bracket_m: [lower, upper],
                }
                .into());
            }
            candidate.queried_distance_m = candidate.actual_distance_m;
            candidate.distance_residual_m = 0.0;
            return Ok(candidate);
        }
        solve_distance(motion, |r| {
            let candidate = evaluate(r)?;
            let actual = candidate.actual_distance_m;
            Ok((candidate, actual))
        })
    }
}
