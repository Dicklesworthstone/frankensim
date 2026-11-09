//! Scenario time histories on an explicitly tagged thermal boundary.
//!
//! Signal and typed scalar/vector interpolation, declared tails and quantity/frame
//! validation remain owned by fs-scenario. Spatial profiles remain refused:
//! this consumes its existing thermal vocabulary without a second schema or
//! an implicit spatial sampling rule. Frames do not move the mesh partition.
//!
//! Sampling preserves the original face partition, boundary kinds, static nodal
//! data and lowest-declared-region precedence at shared Dirichlet vertices.
//! Robin reference temperature is an independently dimensioned history, or the
//! explicitly supplied scenario ambient. No ambient history is inferred from h.

/// Resumable, energy-checked backward Euler driven by these boundary histories.
pub mod march;

use super::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use crate::{ConductionError, TEMPERATURE_DIMS};
use crate::mesh::ConductionMesh;
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_scenario::FrameTree;
use fs_scenario::payload::Payload;
use fs_scenario::scalar::ScalarEvaluationError;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
use std::borrow::Cow;
use std::collections::BTreeSet;
use std::ops::ControlFlow;

impl ThermalBc {
    /// Lower a scenario row at an explicitly selected time in seconds.
    ///
    /// Uniform rows retain their steady semantics. A signal uses its declared
    /// interpolation and clamping; the steady entry point still refuses to
    /// choose an implicit instant. `t_ref` has the same kelvin/ambient semantics
    /// as [`Self::from_scenario_row`]. Typed values use the implicit world tree;
    /// use [`Self::from_scenario_row_at_in_frames`] for declared non-world frames.
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
        Self::from_scenario_row_at_in_frames(
            row, environment, t_ref, &FrameTree::new(), time_s,
            &mut || ControlFlow::Continue(()),
        )
    }

    /// Sample a thermal row using its original immutable frame tree.
    /// Source values can clamp or wrap; frames always use physical time.
    /// Typed temperatures preserve absolute-temperature semantics. Heat-flux
    /// vectors rotate to WORLD; the original native faces project their outward
    /// component without changing the partition. Scalar signs and the separate
    /// Robin ambient retain their meaning. Original sample/frame checkpoints apply.
    pub fn from_scenario_row_at_in_frames(
        row: &BoundaryCondition,
        environment: &Environment,
        t_ref: Option<f64>,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ConductionError> {
        check(checkpoint)?;
        crate::require_finite("boundary evaluation time", time_s)?;
        let mut violations = Vec::new();
        row.check(&mut violations);
        if let Some(first) = violations.into_iter().next() {
            return Err(ConductionError::ScenarioRow {
                region: row.region.clone(), what: first.what, fix: first.fix,
            });
        }
        if matches!(&row.value, Some(BcValue::Typed(Payload::Vector(_)))) {
            return typed_flux(row, frames, time_s, checkpoint);
        }
        let value = match &row.value {
            Some(BcValue::Signal(signal)) => signal.eval(time_s)
                .map_err(|error| row_error(row, error.to_string()))?,
            Some(BcValue::Typed(_)) => typed_value(row, frames, time_s, checkpoint)?,
            _ => super::uniform_value(row)?,
        };
        check(checkpoint)?;
        Self::from_scenario_value(row, environment, t_ref, value)
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
    fn at(&self, cx: &Cx<'_>, environment: &Environment, frames: &FrameTree,
        time_s: f64) -> Result<ThermalBc, ConductionError>
    {
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
        ThermalBc::from_scenario_row_at_in_frames(
            &self.row, environment, reference, frames, time_s, &mut || checkpoint(cx),
        )
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
    frames: Cow<'m, FrameTree>,
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
        let mesh = builder.mesh;
        Self::bind(cx, mesh, builder.finish()?, environment,
            Cow::Owned(FrameTree::new()), histories, initial_time_s)
    }

    /// Attach typed histories with an explicitly borrowed scenario frame tree.
    /// Its immutable borrow lasts for the trajectory; no motion table is copied.
    /// Frames rotate physical flux vectors but never retarget faces.
    pub fn new_in_frames(
        cx: &Cx<'_>, builder: ThermalBoundaryBuilder<'m>,
        environment: &Environment, frames: &'m FrameTree,
        histories: Vec<ThermalRegionHistory>, initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let mesh = builder.mesh;
        Self::bind(cx, mesh, builder.finish()?, environment,
            Cow::Borrowed(frames), histories, initial_time_s)
    }

    /// Attach histories to the original already-lowered face partition.
    /// The caller supplies the same mesh used to build that boundary; no
    /// selector is rerun and shared Dirichlet precedence stays unchanged.
    pub fn from_boundary(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, boundary: &ThermalBoundary,
        environment: &Environment, histories: Vec<ThermalRegionHistory>,
        initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        Self::bind_boundary(cx, mesh, boundary, environment,
            Cow::Owned(FrameTree::new()), histories, initial_time_s)
    }

    /// Bind the original partition and explicit frame tree to one trajectory.
    pub fn from_boundary_in_frames(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, boundary: &ThermalBoundary,
        environment: &Environment, frames: &'m FrameTree,
        histories: Vec<ThermalRegionHistory>, initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        Self::bind_boundary(cx, mesh, boundary, environment,
            Cow::Borrowed(frames), histories, initial_time_s)
    }

    fn bind_boundary(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, boundary: &ThermalBoundary,
        environment: &Environment, frames: Cow<'m, FrameTree>,
        histories: Vec<ThermalRegionHistory>, initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if boundary.face_region.len() != mesh.boundary().len() {
            return Err(ConductionError::ScenarioRow {
                region: "boundary-history".into(),
                what: "retained thermal partition differs from the supplied mesh".into(),
                fix: "keep the original mesh and its lowered boundary together".into(),
            });
        }
        let mut partition = ThermalBoundaryBuilder {
            mesh, names: boundary.names.clone(), conditions: boundary.conditions.clone(),
            face_region: boundary.face_region.clone(),
            adiabatic_remainder: boundary.adiabatic_remainder != 0,
        }.finish()?;
        partition.surface_conduction = boundary.surface_conduction.clone();
        Self::bind(cx, mesh, partition, environment, frames, histories, initial_time_s)
    }

    fn bind(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, partition: ThermalBoundary,
        environment: &Environment, frames: Cow<'m, FrameTree>,
        histories: Vec<ThermalRegionHistory>, initial_time_s: f64,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        crate::require_finite("initial boundary time", initial_time_s)?;
        if frames.frames.len() > 256 {
            return Err(ConductionError::Config {
                parameter: "thermal history frames",
                what: "thermal histories admit at most 256 frame rows".into(),
            });
        }
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
            let sampled = history.at(cx, environment, &frames, initial_time_s)?;
            if sampled.tag() != partition.conditions[index].tag() {
                return Err(row_error(&history.row, "history changes the declared boundary kind".into()));
            }
            rows.push((index, history));
        }
        poll(cx, rows.len())?;
        Ok(Self { mesh, partition, environment: environment.clone(), frames, rows, initial_time_s })
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
            conditions[*index] = history.at(cx, &self.environment, &self.frames, time_s)?;
        }
        poll(cx, self.rows.len())?;
        // The existing finish owns shared-vertex precedence. Reuse it rather
        // than updating only the condition array and leaving a stale lift.
        let mut boundary = ThermalBoundaryBuilder {
            mesh: self.mesh, names: self.partition.names.clone(), conditions,
            face_region: self.partition.face_region.clone(),
            adiabatic_remainder: self.partition.adiabatic_remainder != 0,
        }.finish()?;
        boundary.surface_conduction = self.partition.surface_conduction.clone();
        poll(cx, self.rows.len())?;
        Ok(boundary)
    }
}

fn typed_flux(
    row: &BoundaryCondition, frames: &FrameTree, time_s: f64,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<ThermalBc, ConductionError> {
    let flux = row.heat_flux_world_at(frames, time_s, checkpoint)
        .map_err(|error| sample_error(row, error))?;
    ThermalBc::neumann_vector(flux)
}

fn sample_error(row: &BoundaryCondition, error: ScalarEvaluationError) -> ConductionError {
    match error {
        ScalarEvaluationError::Cancelled => ConductionError::Cancelled {
            stage: "thermal-boundary-history", at: 0,
        },
        other => row_error(row, other.to_string()),
    }
}

fn typed_value(
    row: &BoundaryCondition, frames: &FrameTree, time_s: f64,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<QtyAny, ConductionError> {
    let sample = row.thermal_at(frames, time_s, checkpoint)
        .map_err(|error| sample_error(row, error))?;
    let dims = sample.contract.dims()
        .ok_or_else(|| row_error(row, "thermal scalar must have homogeneous dimensions".into()))?;
    Ok(QtyAny::new(sample.value, dims))
}

fn checkpoint(cx: &Cx<'_>) -> ControlFlow<()> {
    if cx.checkpoint().is_err() { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
}

fn check(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ConductionError> {
    if checkpoint().is_break() {
        Err(ConductionError::Cancelled { stage: "thermal-boundary-history", at: 0 })
    } else { Ok(()) }
}

fn row_error(row: &BoundaryCondition, what: String) -> ConductionError {
    ConductionError::ScenarioRow {
        region: row.region.clone(), what,
        fix: "supply a valid thermal history on a matching named region, with a scalar or three-component Neumann flux and kelvin reference data only for Robin".into(),
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "thermal-boundary-history", at })
}
