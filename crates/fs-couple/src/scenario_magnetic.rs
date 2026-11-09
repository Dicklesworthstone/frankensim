//! Canonical magnetic histories on the original three-dimensional field owners.
//!
//! Geometry binding is explicit and uses original exterior face identities.
//! The native operators retain material sources, compatibility, gauges, true
//! residuals and energy recovery. A time sample is a fixed-geometry static
//! solve: sampling a history does not infer induction or deposit heat.

use fs_em::spatial::{
    MagneticFluxFace, SpatialMagneticBoundary, SpatialMagneticControls, SpatialMagneticError,
    SpatialMagneticField, SpatialMagnetostatics,
};
use fs_qty::Dims;
use fs_scenario::payload::{OrientationParity, Payload, QuantityContract};
use fs_scenario::scalar::ScalarEvaluationError;
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics};
use std::collections::BTreeSet;
use std::ops::ControlFlow;

const FLUX_DENSITY: Dims = Dims([0, 1, -2, 0, -1, 0]);

/// The sign of the authored normal flux density relative to the original body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PositiveMagneticFlux {
    IntoBody,
    OutOfBody,
}

/// A canonical normal-flux row explicitly bound to original exterior faces.
#[derive(Clone, Copy, Debug)]
pub struct MagneticFluxBoundary3<'a> {
    pub faces: &'a [usize],
    pub row: &'a BoundaryCondition,
    pub positive: PositiveMagneticFlux,
}

/// A refused boundary or solve publishes no partial field.
#[derive(Debug)]
pub enum ScenarioMagneticError {
    Invalid(&'static str),
    Scenario(String),
    Sample(ScalarEvaluationError),
    Native(SpatialMagneticError),
    Cancelled,
}
impl core::fmt::Display for ScenarioMagneticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario magnetic boundary refused: {self:?}")
    }
}
impl std::error::Error for ScenarioMagneticError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ScenarioMagneticError> {
    if cancelled() {
        Err(ScenarioMagneticError::Cancelled)
    } else {
        Ok(())
    }
}

fn check_row(row: &BoundaryCondition, kind: BcKind) -> Result<(), ScenarioMagneticError> {
    if row.physics != Physics::Magnetics || row.kind != kind {
        return Err(ScenarioMagneticError::Invalid(
            "magnetic patch role must match its canonical scenario row",
        ));
    }
    let mut violations = Vec::new();
    row.check(&mut violations);
    if let Some(first) = violations.into_iter().next() {
        return Err(ScenarioMagneticError::Scenario(first.what));
    }
    Ok(())
}

fn sample_error(error: ScalarEvaluationError) -> ScenarioMagneticError {
    match error {
        ScalarEvaluationError::Cancelled => ScenarioMagneticError::Cancelled,
        error => ScenarioMagneticError::Sample(error),
    }
}

fn native_error(error: SpatialMagneticError) -> ScenarioMagneticError {
    match error {
        SpatialMagneticError::Cancelled => ScenarioMagneticError::Cancelled,
        error => ScenarioMagneticError::Native(error),
    }
}

/// Lower normal-flux histories together at one physical time.
///
/// The scalar is B.n [T], not total flux [Wb]. Because B is axial and the
/// geometric normal is polar, its canonical payload must declare odd parity.
/// Proper scenario frame rotations do not change this oriented scalar; only
/// the explicit `positive` convention changes its sign. A sample is uniform
/// over its named patch and the native owner integrates its actual face areas.
///
/// At most 1024 nonempty, nonoverlapping patches within the original face
/// count are admitted. Unselected faces retain zero outward flux. Net flux
/// compatibility is deliberately not repaired here: the native current-free
/// solver checks every connected component before installing its gauge.
pub fn magnetic_flux_boundary_at(
    problem: &SpatialMagnetostatics,
    rows: &[MagneticFluxBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialMagneticBoundary, ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    if !time_s.is_finite() || rows.len() > 1024 {
        return Err(E::Invalid(
            "finite magnetic sample time and at most 1024 rows required",
        ));
    }
    let mut count = 0usize;
    for input in rows {
        poll(cancelled)?;
        count = count
            .checked_add(input.faces.len())
            .ok_or(E::Invalid("magnetic face count overflow"))?;
        if input.faces.is_empty() || count > problem.mesh().faces.len() {
            return Err(E::Invalid(
                "nonempty magnetic patches within the original face count required",
            ));
        }
    }
    let mut boundary = SpatialMagneticBoundary::default();
    let mut assigned = BTreeSet::new();
    for input in rows {
        poll(cancelled)?;
        check_row(input.row, BcKind::NormalMagneticFluxDensity)?;
        let Some(BcValue::Typed(Payload::Scalar(value))) = &input.row.value else {
            return Err(E::Invalid(
                "normal magnetic flux requires a canonical typed scalar",
            ));
        };
        if value.meta().contract() != QuantityContract::Dimensions(FLUX_DENSITY)
            || value.meta().orientation() != OrientationParity::Odd
        {
            return Err(E::Invalid(
                "normal magnetic flux requires dimensions-only tesla and odd parity",
            ));
        }
        let sample = value
            .oriented_scalar_at(frames, time_s, &mut || {
                if cancelled() {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            })
            .map_err(sample_error)?;
        let sign = if input.positive == PositiveMagneticFlux::OutOfBody {
            1.0
        } else {
            -1.0
        };
        for &face in input.faces {
            poll(cancelled)?;
            if problem.exterior_normal(face).is_none() || !assigned.insert(face) {
                return Err(E::Invalid(
                    "magnetic patches need distinct original exterior faces",
                ));
            }
            boundary.flux.push(MagneticFluxFace {
                face,
                outward_flux_density_t: sign * sample.value,
            });
        }
    }
    poll(cancelled)?;
    Ok(boundary)
}

/// Native fixed-state field with the actual boundary and physical sample time.
pub struct ScenarioMagneticFlux3<'m> {
    pub time_s: f64,
    pub boundary: SpatialMagneticBoundary,
    pub field: SpatialMagneticField<'m>,
}

/// Execute an authored normal-flux snapshot on the original current-free mesh.
///
/// The native owner retains its single-valued scalar-potential restriction,
/// material laws, component flux balance, gauge and residual acceptance. The
/// returned B/H and stored energy can feed the existing field consumers.
pub fn solve_magnetic_flux3<'m>(
    problem: &'m SpatialMagnetostatics,
    rows: &[MagneticFluxBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    controls: SpatialMagneticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScenarioMagneticFlux3<'m>, ScenarioMagneticError> {
    let boundary = magnetic_flux_boundary_at(problem, rows, frames, time_s, cancelled)?;
    let field = problem
        .solve(&boundary, controls, cancelled)
        .map_err(native_error)?;
    poll(cancelled)?;
    Ok(ScenarioMagneticFlux3 {
        time_s,
        boundary,
        field,
    })
}
