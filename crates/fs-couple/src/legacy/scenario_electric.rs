//! Canonical electrical scenario rows on the original 3-D Ohmic field owner.
//!
//! Geometry binding is explicit: a named patch is resolved to original exterior
//! faces by the caller. Potential and normal-current histories are sampled at
//! one physical time. The source-bound native operator still owns conductivity,
//! contact traces, gauges, charge continuity and local irreversible Joule heat.
//! These are fixed-geometry DC snapshots, not capacitive or inductive dynamics.

use fs_em::conduction::ConductionControls;
use fs_em::spatial::conduction::{
    ConductionCurrentFace, ConductionVoltageFace, SpatialConduction, SpatialConductionBoundary,
    SpatialConductionError, SpatialConductionField,
};
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_qty::Dims;
use fs_scenario::payload::{Payload, QuantityContract};
use fs_scenario::scalar::ScalarEvaluationError;
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics};
use std::collections::BTreeSet;
use std::ops::ControlFlow;

const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const CURRENT_DENSITY: Dims = Dims([-2, 0, 0, 0, 1, 0]);

/// The scenario row does not encode a normal direction. Bind it explicitly;
/// rotating component frames never reverses this geometric sign convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PositiveCurrent {
    IntoBody,
    OutOfBody,
}

/// Borrowed canonical row on an explicitly resolved exterior patch.
#[derive(Clone, Copy, Debug)]
pub enum ElectricBoundary3<'a> {
    /// Equipotential value on every selected face, in volts.
    Potential {
        faces: &'a [usize],
        row: &'a BoundaryCondition,
    },
    /// Uniform normal current density, in A/m², not a total current in amperes.
    NormalCurrent {
        faces: &'a [usize],
        row: &'a BoundaryCondition,
        positive: PositiveCurrent,
    },
}

#[derive(Debug)]
pub enum ScenarioElectricError {
    Invalid(&'static str),
    Scenario(String),
    Sample(ScalarEvaluationError),
    Native(SpatialConductionError),
    Cancelled,
}
impl core::fmt::Display for ScenarioElectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario electrical boundary refused: {self:?}")
    }
}
impl std::error::Error for ScenarioElectricError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ScenarioElectricError> {
    if cancelled() {
        Err(ScenarioElectricError::Cancelled)
    } else {
        Ok(())
    }
}

/// Lower every row together to the original native boundary representation.
///
/// The resulting boundary can also drive the existing spatial electrothermal
/// owner (`SpatialConductionDrive::Boundary`), so no electrical or heat law is
/// copied into scenario execution. Unselected faces retain the native insulating
/// default. A row is sampled once and its value broadcast over its bound patch;
/// no face area is used to reinterpret a density as a total current.
///
/// There are at most 1024 rows and at most one assignment per original exterior
/// face. Unsupported role, units, source/frame/reference policy, overlapping or
/// empty patches, invalid time and cancellation publish no partial boundary.
/// Shared-vertex voltage agreement, material-contact exclusions and component
/// current balance remain with the native solve.
pub fn conduction_boundary_at<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    rows: &[ElectricBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialConductionBoundary, ScenarioElectricError> {
    use ScenarioElectricError as E;
    poll(cancelled)?;
    if !time_s.is_finite() || rows.len() > 1024 {
        return Err(E::Invalid(
            "finite electrical sample time and at most 1024 rows required",
        ));
    }
    let mut count = 0usize;
    for row in rows {
        poll(cancelled)?;
        let faces = match row {
            ElectricBoundary3::Potential { faces, .. }
            | ElectricBoundary3::NormalCurrent { faces, .. } => faces,
        };
        count = count
            .checked_add(faces.len())
            .ok_or(E::Invalid("electrical face count overflow"))?;
        if faces.is_empty() || count > problem.mesh().faces.len() {
            return Err(E::Invalid(
                "nonempty electrical patches within the original face count required",
            ));
        }
    }
    let mut boundary = SpatialConductionBoundary::default();
    let mut assigned = BTreeSet::new();
    for input in rows {
        poll(cancelled)?;
        let (faces, row, kind, dims, sign) = match input {
            ElectricBoundary3::Potential { faces, row } => {
                (*faces, *row, BcKind::ElectricPotential, VOLTAGE, 1.0)
            }
            ElectricBoundary3::NormalCurrent {
                faces,
                row,
                positive,
            } => (
                *faces,
                *row,
                BcKind::NormalCurrentDensity,
                CURRENT_DENSITY,
                if *positive == PositiveCurrent::IntoBody {
                    1.0
                } else {
                    -1.0
                },
            ),
        };
        let mut violations = Vec::new();
        row.check(&mut violations);
        if let Some(first) = violations.into_iter().next() {
            return Err(E::Scenario(first.what));
        }
        if row.physics != Physics::Electrics || row.kind != kind {
            return Err(E::Invalid(
                "electrical patch role must match its canonical scenario row",
            ));
        }
        let Some(BcValue::Typed(Payload::Scalar(value))) = &row.value else {
            return Err(E::Invalid(
                "electrical boundary requires the canonical typed scalar carrier",
            ));
        };
        if value.meta().contract() != QuantityContract::Dimensions(dims) {
            return Err(E::Invalid(
                "electrical boundary requires its exact dimensions-only SI contract",
            ));
        }
        let sample = value
            .scalar_at(frames, time_s, &mut || {
                if cancelled() {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            })
            .map_err(|error| match error {
                ScalarEvaluationError::Cancelled => E::Cancelled,
                error => E::Sample(error),
            })?;
        for &face in faces {
            poll(cancelled)?;
            if problem.exterior_normal(face).is_none() || !assigned.insert(face) {
                return Err(E::Invalid(
                    "electrical patches need distinct original exterior faces",
                ));
            }
            if kind == BcKind::ElectricPotential {
                boundary.voltage.push(ConductionVoltageFace {
                    face,
                    voltage_v: [sample.value; 3],
                });
            } else {
                boundary.current.push(ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: sign * sample.value,
                });
            }
        }
    }
    poll(cancelled)?;
    Ok(boundary)
}

/// Native field, its actual controls, and the physical sample time stay together.
#[derive(Debug)]
pub struct ScenarioConduction3 {
    pub time_s: f64,
    pub boundary: SpatialConductionBoundary,
    pub field: SpatialConductionField,
}

/// Sample a complete scenario boundary and execute its original material solve.
///
/// Local current, bulk/contact heat and supplied power come directly from the
/// native field result. Source material and geometry are borrowed unchanged.
/// No clock is advanced and no energy is deposited merely by taking a snapshot.
pub fn solve_conduction3<C: SpatialConductorMaterial>(
    problem: &SpatialConduction<C>,
    rows: &[ElectricBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    controls: ConductionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScenarioConduction3, ScenarioElectricError> {
    let boundary = conduction_boundary_at(problem, rows, frames, time_s, cancelled)?;
    let field = problem
        .solve(&boundary, controls, cancelled)
        .map_err(|error| match error {
            SpatialConductionError::Cancelled => ScenarioElectricError::Cancelled,
            error => ScenarioElectricError::Native(error),
        })?;
    poll(cancelled)?;
    Ok(ScenarioConduction3 {
        time_s,
        boundary,
        field,
    })
}
