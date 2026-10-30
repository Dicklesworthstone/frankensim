//! Original coefficient queries closed against continuous sliding and arrest.
//!
//! A phase holds each queried coefficient fixed while the native law follows
//! the complete curved velocity path. The source speed is unsigned path / phase
//! duration; the source distance is its own accepted origin plus that path.
//! Refining phases resolves changing sources. A single phase is not an exact
//! integration of a continuously temperature-, speed- or path-dependent law.

use super::motion::{
    MotionQueryError, MotionQuerySolve, admissible_distance, finite, query_for_distance,
    solve_distance, validate_controls,
};
use super::*;
use fs_tribo::tangential::elliptic::relaxation::{
    EllipticRelaxationControls, EllipticRelaxationError, EllipticRelaxationInput,
    EllipticRelaxationReport,
};

/// One original-source and native trajectory candidate. Neither owner mutates.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceRelaxationReport {
    pub relaxation: EllipticRelaxationReport,
    pub state: SourceThermalEllipticState,
    pub queried_distance_m: f64,
    /// Unsigned path along the native curve, not the norm of net displacement.
    pub actual_distance_m: f64,
    /// Query minus resolved unsigned path; accepted motion roots are nonpositive.
    pub distance_residual_m: f64,
    /// Original-source/native trajectory evaluations, including bracket ends.
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SourceRelaxationError {
    Query(MotionQueryError),
    Relaxation(EllipticRelaxationError),
}
impl From<MotionQueryError> for SourceRelaxationError {
    fn from(e: MotionQueryError) -> Self {
        Self::Query(e)
    }
}
impl From<SourceThermalCoulombError> for SourceRelaxationError {
    fn from(e: SourceThermalCoulombError) -> Self {
        Self::Query(MotionQueryError::Source(e))
    }
}
impl core::fmt::Display for SourceRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional trajectory: {self:?}")
    }
}
impl std::error::Error for SourceRelaxationError {}

impl SourceThermalEllipticCoulomb {
    /// Resolve original source motion against the complete finite-arrest curve.
    ///
    /// `thermal_query` carries the complete original schema and physical T/P;
    /// this owner supplies every declared speed and accumulated-path axis.
    /// A source trial refreshes both pinned principal coefficient pairs, retains
    /// the accepted native texture/history and solves the full body mobility.
    /// The safeguarded scalar root is the same one used by the held-step owner.
    /// Each accepted query must agree with actual unsigned travel, from below,
    /// within the explicit distance tolerance. Nonmonotone source laws do not
    /// imply existence, uniqueness or selection of the smallest root.
    ///
    /// Native stopping-time and quadrature estimates remain separate from this
    /// source residual. Coefficients, load and geometry are constant inside this
    /// phase; source/thermal variation requires a convergent outer subdivision.
    #[allow(clippy::too_many_arguments)]
    pub fn relax_motion_consistent(
        &self,
        state: &SourceThermalEllipticState,
        input: &EllipticRelaxationInput,
        thermal_query: &QueryPoint,
        controls: MotionQuerySolve,
        relaxation: EllipticRelaxationControls,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<SourceRelaxationReport, SourceRelaxationError> {
        if cancel() {
            return Err(MotionQueryError::Cancelled.into());
        }
        self.check(state)?;
        validate_controls(controls, input.dt_s)?;
        let axes = self.principal[0].query_axes();
        let has_motion = axes.iter().any(|axis| {
            matches!(
                axis.kind,
                FrictionQueryAxisKind::SlidingSpeed | FrictionQueryAxisKind::AccumulatedDistance
            )
        });
        let accepted = state.sources[0].query();
        let mut evaluations = 0;
        let mut evaluate = |r| -> Result<SourceRelaxationReport, SourceRelaxationError> {
            if cancel() {
                return Err(MotionQueryError::Cancelled.into());
            }
            if evaluations >= controls.maximum_evaluations {
                return Err(MotionQueryError::EvaluationBudget {
                    evaluations,
                    maximum: controls.maximum_evaluations,
                }
                .into());
            }
            evaluations += 1;
            let point = query_for_distance(thermal_query, &accepted, axes, r, input.dt_s)?;
            let sources = [
                self.principal[0].source_at_query_conditions(
                    &state.sources[0],
                    input.contact_temperature_k,
                    input.nominal_pressure_pa,
                    &point,
                )?,
                self.principal[1].source_at_query_conditions(
                    &state.sources[1],
                    input.contact_temperature_k,
                    input.nominal_pressure_pa,
                    &point,
                )?,
            ];
            let (static_mu, kinetic_mu) =
                principal_coefficients(&self.principal, [&sources[0], &sources[1]])?;
            let (law, history) = state
                .law
                .continue_coefficients(&state.constitutive, static_mu, kinetic_mu)
                .map_err(SourceThermalCoulombError::from)?;
            let trajectory = law
                .relax_velocity(&history, input, relaxation, cancel)
                .map_err(SourceRelaxationError::Relaxation)?;
            if cancel() {
                return Err(MotionQueryError::Cancelled.into());
            }
            let actual = trajectory.distance_m;
            Ok(SourceRelaxationReport {
                state: SourceThermalEllipticState {
                    owner: self.identity,
                    law,
                    constitutive: trajectory.state.clone(),
                    sources,
                },
                relaxation: trajectory,
                queried_distance_m: r,
                actual_distance_m: actual,
                distance_residual_m: finite(r - actual)?,
                evaluations,
            })
        };
        let in_bracket = |r| -> Result<(), MotionQueryError> {
            let [lower, upper] = controls.distance_bracket_m;
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
        if input.normal_n == 0.0 || input.velocity_m_s == [0.0; 2] {
            let r = finite(
                fs_math::det::hypot(input.velocity_m_s[0], input.velocity_m_s[1]) * input.dt_s,
            )?;
            in_bracket(r)?;
            let candidate = evaluate(r)?;
            if admissible_distance(controls, r, candidate.actual_distance_m)? {
                return Ok(candidate);
            }
            return Err(MotionQueryError::Unresolved {
                bracket_m: [r; 2],
                residual_m: [candidate.distance_residual_m; 2],
            }
            .into());
        }
        solve_distance(controls, |r| {
            let candidate = evaluate(r)?;
            let actual = candidate.actual_distance_m;
            Ok((candidate, actual))
        })
    }
}
