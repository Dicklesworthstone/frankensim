//! Canonical magnetic histories on the original three-dimensional field owners.
//!
//! Geometry binding is explicit and uses original exterior face identities.
//! The native operators retain material sources, compatibility, gauges, true
//! residuals and energy recovery. A time sample is a fixed-geometry static
//! solve: sampling a history does not infer induction or deposit heat.

/// Affine scenario traces on the original retained Faraday-induction step.
pub mod induction;

use fs_em::spatial::vector_potential::{
    SpatialVectorField, SpatialVectorMagnetostatics, VectorMagneticBoundary,
    VectorMagneticControls, VectorPotentialFace,
};
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
const VECTOR_POTENTIAL: Dims = Dims([1, 1, -2, 0, -1, 0]);

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

/// Explicit spatial interpretation of canonical vector-potential histories.
#[derive(Clone, Copy, Debug)]
pub enum MagneticPotentialBoundary3<'a> {
    /// One spatially uniform A vector on each selected original face.
    /// A globally uniform vector is pure gauge, not an imposed nonzero B.
    Uniform {
        faces: &'a [usize],
        row: &'a BoundaryCondition,
    },
    /// Authored line-average A vectors along the original face's three edges.
    /// For sorted vertices [a,b,c], the order is [(a,b),(a,c),(b,c)].
    /// These are actual edge averages, not inferred midpoint point samples.
    EdgeAverages {
        face: usize,
        values: [&'a BoundaryCondition; 3],
    },
}

fn potential_sample(
    row: &BoundaryCondition,
    frames: &FrameTree,
    time_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    check_row(row, BcKind::MagneticVectorPotential)?;
    let Some(BcValue::Typed(Payload::Vector(value))) = &row.value else {
        return Err(E::Invalid(
            "magnetic vector potential requires a canonical typed vector",
        ));
    };
    if value.meta().contract() != QuantityContract::Dimensions(VECTOR_POTENTIAL)
        || value.meta().orientation() != OrientationParity::Even
    {
        return Err(E::Invalid(
            "magnetic vector potential requires dimensions-only Wb/m and even parity",
        ));
    }
    value
        .world_vector_at(frames, time_s, &mut || {
            if cancelled() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .map(|sample| sample.values)
        .map_err(sample_error)
}

fn potential_face(
    problem: &SpatialVectorMagnetostatics,
    face: usize,
    values: [[f64; 3]; 3],
    assigned: &mut BTreeSet<usize>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<VectorPotentialFace, ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    if problem.exterior_normal(face).is_none() || !assigned.insert(face) {
        return Err(E::Invalid(
            "vector-potential patches need distinct original exterior faces",
        ));
    }
    let edges = problem
        .face_edges(face)
        .ok_or(E::Invalid("original face edges are unavailable"))?;
    let mut edge_integrals_wb = [0.0; 3];
    for (index, edge) in edges.into_iter().enumerate() {
        poll(cancelled)?;
        let [a, b] = problem.mesh().edges[edge];
        let p = problem.positions_m()[a as usize];
        let q = problem.positions_m()[b as usize];
        let delta: [f64; 3] = std::array::from_fn(|axis| q[axis] - p[axis]);
        let a = values[index];
        let integral = a[0].mul_add(delta[0], a[1].mul_add(delta[1], a[2] * delta[2]));
        if !delta.iter().all(|value| value.is_finite()) || !integral.is_finite() {
            return Err(E::Invalid("unrepresentable magnetic edge integral"));
        }
        edge_integrals_wb[index] = integral;
    }
    Ok(VectorPotentialFace {
        face,
        edge_integrals_wb,
    })
}

/// Sample canonical A histories and integrate them on the original Whitney edges.
///
/// Mesh positions and sampled vectors share WORLD Cartesian axes. Translation
/// does not act on A; proper frame rotation acts once on its components.
/// The integral is A_average dot (x_b-x_a) [Wb], with the native canonical edge
/// orientation. Shared-edge agreement, current compatibility and coexact gauge
/// remain with the native field or induction solve; conflicting traces refuse
/// there rather than being silently averaged by this adapter.
///
/// Uniform face vectors have zero face circulation. Use explicit line averages
/// to prescribe nonzero flux; a general spatial profile is never guessed from
/// one vector. Unselected faces retain the native natural condition. At most
/// 1024 bindings and one assignment per original exterior face are admitted.
pub fn vector_potential_boundary_at(
    problem: &SpatialVectorMagnetostatics,
    rows: &[MagneticPotentialBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<VectorMagneticBoundary, ScenarioMagneticError> {
    use ScenarioMagneticError as E;
    poll(cancelled)?;
    if !time_s.is_finite() || rows.len() > 1024 {
        return Err(E::Invalid(
            "finite magnetic sample time and at most 1024 bindings required",
        ));
    }
    let mut count = 0usize;
    for input in rows {
        poll(cancelled)?;
        let added = match input {
            MagneticPotentialBoundary3::Uniform { faces, .. } => faces.len(),
            MagneticPotentialBoundary3::EdgeAverages { .. } => 1,
        };
        count = count
            .checked_add(added)
            .ok_or(E::Invalid("magnetic face count overflow"))?;
        if added == 0 || count > problem.mesh().faces.len() {
            return Err(E::Invalid(
                "nonempty magnetic patches within the original face count required",
            ));
        }
    }
    let mut boundary = VectorMagneticBoundary::default();
    let mut assigned = BTreeSet::new();
    for input in rows {
        poll(cancelled)?;
        match input {
            MagneticPotentialBoundary3::Uniform { faces, row } => {
                let sample = potential_sample(row, frames, time_s, cancelled)?;
                for &face in *faces {
                    boundary.potential.push(potential_face(
                        problem,
                        face,
                        [sample; 3],
                        &mut assigned,
                        cancelled,
                    )?);
                }
            }
            MagneticPotentialBoundary3::EdgeAverages { face, values } => {
                let mut samples = [[0.0; 3]; 3];
                for (sample, row) in samples.iter_mut().zip(values) {
                    *sample = potential_sample(row, frames, time_s, cancelled)?;
                }
                boundary.potential.push(potential_face(
                    problem,
                    *face,
                    samples,
                    &mut assigned,
                    cancelled,
                )?);
            }
        }
    }
    poll(cancelled)?;
    Ok(boundary)
}

/// Native current-carrying field with its actual boundary and sample time.
pub struct ScenarioVectorPotential3<'m> {
    pub time_s: f64,
    pub boundary: VectorMagneticBoundary,
    pub field: SpatialVectorField<'m>,
}

/// Execute the existing current-carrying 3D owner from canonical A histories.
///
/// Original impressed currents, material tensors, source identities, relative
/// topology and gauges remain native. This is a static snapshot, not a time
/// step; its field can initialize the original Faraday-induction owner.
#[allow(clippy::too_many_arguments)]
pub fn solve_vector_potential3<'m>(
    problem: &'m SpatialVectorMagnetostatics,
    rows: &[MagneticPotentialBoundary3<'_>],
    frames: &FrameTree,
    time_s: f64,
    current_a_m2: &[[f64; 3]],
    controls: VectorMagneticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScenarioVectorPotential3<'m>, ScenarioMagneticError> {
    let boundary = vector_potential_boundary_at(problem, rows, frames, time_s, cancelled)?;
    let field = problem
        .solve(current_a_m2, &boundary, controls, cancelled)
        .map_err(native_error)?;
    poll(cancelled)?;
    Ok(ScenarioVectorPotential3 {
        time_s,
        boundary,
        field,
    })
}
