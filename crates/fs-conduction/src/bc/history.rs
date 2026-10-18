//! Scenario time histories on an explicitly tagged thermal boundary.
//!
//! Signal interpolation, endpoint clamping and dimensional validation remain
//! owned by fs-scenario. Spatial profiles and typed payloads remain refused:
//! this consumes the scalar/time-signal vocabulary admitted for Thermal rows,
//! not a second scenario schema or an implicit spatial sampling rule.
//!
//! Sampling preserves the original face partition, boundary kinds, static nodal
//! data and lowest-declared-region precedence at shared Dirichlet vertices.
//! Robin reference temperature is an independently dimensioned history, or the
//! explicitly supplied scenario ambient. No ambient history is inferred from h.

use super::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use crate::{ConductionError, TEMPERATURE_DIMS};
use crate::mesh::ConductionMesh;
use fs_exec::Cx;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
use std::collections::BTreeSet;

impl ThermalBc {
    /// Lower a scenario row at an explicitly selected time in seconds.
    ///
    /// Uniform rows retain their steady semantics. A signal uses its declared
    /// interpolation and clamping; the steady entry point still refuses to
    /// choose an implicit instant. `t_ref` has the same kelvin/ambient semantics
    /// as [`Self::from_scenario_row`].
    ///
    /// # Errors
    /// Nonfinite time, invalid scenario row/signal, incompatible dimensions, or
    /// any steady lowering refusal. Spatial profiles are not time histories.
    pub fn from_scenario_row_at(
        row: &BoundaryCondition,
        environment: &Environment,
        t_ref: Option<f64>,
        time_s: f64,
    ) -> Result<Self, ConductionError> {
        crate::require_finite("boundary evaluation time", time_s)?;
        let mut violations = Vec::new();
        row.check(&mut violations);
        if let Some(first) = violations.into_iter().next() {
            return Err(ConductionError::ScenarioRow {
                region: row.region.clone(), what: first.what, fix: first.fix,
            });
        }
        let Some(BcValue::Signal(signal)) = &row.value else {
            return Self::from_scenario_row(row, environment, t_ref);
        };
        let value = signal.eval(time_s).map_err(|error| row_error(row, error.to_string()))?;
        // Do not clone the time table just to replace it with one scalar.
        let sampled = BoundaryCondition {
            region: row.region.clone(), physics: row.physics, kind: row.kind,
            value: Some(BcValue::Uniform(value)), compatibility: row.compatibility,
            frame: row.frame,
        };
        Self::from_scenario_row(&sampled, environment, t_ref)
    }
}

/// One scenario row and its independently declared Robin reference history.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalRegionHistory {
    /// The exact region name must already exist in the supplied partition.
    pub row: BoundaryCondition,
    /// Kelvin-valued history, permitted only for Robin. None explicitly selects
    /// the ambient temperature of the environment passed to the history owner.
    pub reference_temperature: Option<TimeSignal>,
}

impl ThermalRegionHistory {
    fn at(&self, environment: &Environment, time_s: f64) -> Result<ThermalBc, ConductionError> {
        let reference = match &self.reference_temperature {
            None => None,
            Some(signal) => {
                if self.row.kind != BcKind::Robin {
                    return Err(row_error(&self.row, "a reference-temperature history requires a Robin row".into()));
                }
                let value = signal.eval(time_s)
                    .map_err(|error| row_error(&self.row, error.to_string()))?;
                if value.dims != TEMPERATURE_DIMS {
                    return Err(ConductionError::Dimensions {
                        context: format!("reference-temperature history on {:?}", self.row.region),
                        expected: TEMPERATURE_DIMS.0, found: value.dims.0,
                    });
                }
                Some(value.value)
            }
        };
        ThermalBc::from_scenario_row_at(&self.row, environment, reference, time_s)
    }
}

/// Immutable histories bound to the exact mesh of a boundary builder.
///
/// Histories may replace a subset of named rows; all others remain unchanged.
/// No geometry selector is rerun at later times. The initial instant and every
/// requested endpoint are checked, not the continuum between them. This does
/// not establish substep extrema, continuous positivity or temporal accuracy.
#[derive(Debug)]
pub struct ThermalBoundaryHistory<'m> {
    mesh: &'m ConductionMesh,
    partition: ThermalBoundary,
    environment: Environment,
    rows: Vec<(usize, ThermalRegionHistory)>,
    initial_time_s: f64,
}

impl<'m> ThermalBoundaryHistory<'m> {
    /// Finish a mesh-bound partition and attach histories to its existing rows.
    ///
    /// Initial row values in the builder define the static partition/kinds;
    /// selected histories supply their values at `initial_time_s` and later
    /// samples. Unknown, repeated, empty, or kind-changing bindings refuse.
    /// Existing untagged-boundary and row-admission checks remain in force.
    pub fn new(
        cx: &Cx<'_>,
        builder: ThermalBoundaryBuilder<'m>,
        environment: &Environment,
        histories: Vec<ThermalRegionHistory>,
        initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        crate::require_finite("initial boundary time", initial_time_s)?;
        let mesh = builder.mesh;
        let partition = builder.finish()?;
        let mut seen = BTreeSet::new();
        let mut rows = Vec::with_capacity(histories.len());
        for (at, history) in histories.into_iter().enumerate() {
            poll(cx, at)?;
            let name = &history.row.region;
            let Some(index) = partition.names.iter().position(|candidate| candidate == name) else {
                return Err(row_error(&history.row, "history region is absent from the boundary partition".into()));
            };
            if !seen.insert(index) {
                return Err(ConductionError::DuplicateRegion { region: name.clone() });
            }
            if !partition.face_region.contains(&Some(index)) {
                return Err(row_error(&history.row, "history region contains no boundary faces".into()));
            }
            let sampled = history.at(environment, initial_time_s)?;
            if sampled.tag() != partition.conditions[index].tag() {
                return Err(row_error(&history.row, "history changes the declared boundary kind".into()));
            }
            rows.push((index, history));
        }
        poll(cx, rows.len())?;
        Ok(Self { mesh, partition, environment: environment.clone(), rows, initial_time_s })
    }

    /// The original mesh; a time-dependent row cannot retarget its faces.
    #[must_use]
    pub const fn mesh(&self) -> &'m ConductionMesh { self.mesh }

    /// Explicit initial instant, in seconds; also used by a history-driven run.
    #[must_use]
    pub const fn initial_time_s(&self) -> f64 { self.initial_time_s }

    /// Evaluate all selected rows together, rebuilding prescribed vertex values.
    ///
    /// The original partition is never mutated. A bad endpoint or cancellation
    /// returns no partially updated boundary. Cancellation is checked between
    /// existing signal evaluations and around partition reconstruction; this
    /// adds no stronger latency claim to those existing owners.
    pub fn at(&self, cx: &Cx<'_>, time_s: f64) -> Result<ThermalBoundary, ConductionError> {
        poll(cx, 0)?;
        crate::require_finite("boundary evaluation time", time_s)?;
        let mut conditions = self.partition.conditions.clone();
        for (at, (index, history)) in self.rows.iter().enumerate() {
            poll(cx, at)?;
            conditions[*index] = history.at(&self.environment, time_s)?;
        }
        poll(cx, self.rows.len())?;
        // The existing finish owns shared-vertex precedence. Reuse it rather
        // than updating only the condition array and leaving a stale lift.
        let boundary = ThermalBoundaryBuilder {
            mesh: self.mesh, names: self.partition.names.clone(), conditions,
            face_region: self.partition.face_region.clone(),
            adiabatic_remainder: self.partition.adiabatic_remainder != 0,
        }.finish()?;
        poll(cx, self.rows.len())?;
        Ok(boundary)
    }
}

fn row_error(row: &BoundaryCondition, what: String) -> ConductionError {
    ConductionError::ScenarioRow {
        region: row.region.clone(), what,
        fix: "supply a valid thermal scalar/time history on a matching named region, with kelvin reference data only for Robin".into(),
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "thermal-boundary-history", at })
}
