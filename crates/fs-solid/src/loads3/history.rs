//! Scenario histories evaluated on one immutable reference-volume mesh.
//!
//! Traction is a world Cartesian polar vector [Pa]; pressure is a scalar [Pa]
//! acting on the retained reference normals. A declared acceleration uses the
//! CURRENT problem's uniform/per-element reference density, never a second
//! unrelated density table. This is external acceleration loading, not an
//! accelerating-coordinate formulation: no inertial/Coriolis terms are inferred.
//! Sampling yields one quasi-static load case with the original thermal/solver
//! interfaces. It does not integrate mechanical dynamics or update geometry.

use super::*;
use super::fields::AccelerationField3;
use fs_qty::Dims;
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics, TimeSignal};
use fs_scenario::payload::{OrientationParity, Payload, QuantityContract, VectorPayload};
use fs_scenario::vector::VectorEvaluationError;
use std::ops::ControlFlow;

pub const TRACTION_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
pub const FORCE_DENSITY_DIMS: Dims = Dims([-2, 1, -2, 0, 0, 0]);
pub const ACCELERATION_DIMS: Dims = Dims([1, 0, -2, 0, 0, 0]);

#[derive(Debug, Clone, Copy)]
pub enum SurfaceHistory3<'a> {
    Fixed(SurfaceLoad3<'a>),
    Traction { faces: &'a [usize], value: &'a VectorPayload },
    Pressure { faces: &'a [usize], value: &'a TimeSignal },
}
impl<'a> SurfaceHistory3<'a> {
    /// Consume the existing Elasticity/Neumann vector row, on explicitly
    /// resolved exterior faces. Scalar legacy Traction rows are not pressures
    /// or vectors by implication; they retain their existing IR semantics.
    pub fn scenario_traction(faces: &'a [usize], row: &'a BoundaryCondition)
        -> Result<Self, LoadHistoryError3> {
        let mut violations = Vec::new(); row.check(&mut violations);
        if let Some(first) = violations.first() { return Err(LoadHistoryError3::Scenario(first.what.clone())); }
        if row.physics != Physics::Elasticity || row.kind != BcKind::Neumann {
            return Err(LoadHistoryError3::Invalid("only an Elasticity/Neumann vector row is admitted"));
        }
        match &row.value {
            Some(BcValue::Typed(Payload::Vector(value))) => Ok(Self::Traction { faces, value }),
            _ => Err(LoadHistoryError3::Invalid("typed vector traction required")),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum VolumeHistory3<'a> {
    Fixed(VolumeLoad3<'a>),
    ForceDensity(&'a VectorPayload),
    /// Same acceleration at every element; density comes from the solve material.
    Acceleration(&'a VectorPayload),
    /// Held uniform or P1 acceleration in reference-mesh world axes. The same
    /// solve material supplies each element's reference density, including at
    /// material interfaces. This does not resample a spatial field during motion.
    FixedAcceleration(AccelerationField3<'a>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum LoadHistoryError3 {
    Invalid(&'static str),
    Scenario(String),
    Quantity { expected: Dims, found: QuantityContract },
    Vector(VectorEvaluationError),
    Loads(LoadError3),
    Cancelled,
}
impl std::fmt::Display for LoadHistoryError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D load history refused: {self:?}")
    }
}
impl std::error::Error for LoadHistoryError3 {}
impl From<LoadError3> for LoadHistoryError3 {
    fn from(error: LoadError3) -> Self {
        if error == LoadError3::Cancelled { Self::Cancelled } else { Self::Loads(error) }
    }
}
fn sample(value: &VectorPayload, dims: Dims, frames: &FrameTree, time_s: f64, cx: &Cx<'_>)
    -> Result<[f64; 3], LoadHistoryError3> {
    if value.meta().contract() != QuantityContract::Dimensions(dims) {
        return Err(LoadHistoryError3::Quantity { expected: dims, found: value.meta().contract() });
    }
    if value.meta().orientation() != OrientationParity::Even {
        return Err(LoadHistoryError3::Invalid("physical vector histories require ordinary even-parity vectors"));
    }
    let mut checkpoint = || if cx.checkpoint().is_err() { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    value.world_vector_at(frames, time_s, &mut checkpoint).map(|v| v.values).map_err(|error| {
        if error == VectorEvaluationError::Cancelled { LoadHistoryError3::Cancelled }
        else { LoadHistoryError3::Vector(error) }
    })
}

/// Borrowed histories plus explicit static nodal forces, superposed at an endpoint.
/// All time interpolation, outside policies and frame motion remain fs-scenario's.
/// Reusing this value neither advances nor resets a simulation clock.
#[derive(Debug, Clone, Copy)]
pub struct TetLoadHistory3<'a> {
    pub surfaces: &'a [SurfaceHistory3<'a>],
    pub volume: VolumeHistory3<'a>,
    pub nodal_forces_n: &'a [(usize, f64)],
    pub frames: &'a FrameTree,
}

/// An endpoint borrows its exact problem as well as its load geometry. Material
/// density cannot change between acceleration integration and the solve.
pub struct TetLoadSample3<'g, 'm, 'p, 's> {
    case: TetLoadCase3<'g, 'm>,
    problem: &'p TetLinearElasticProblem<'s>,
    time_s: f64,
}
impl TetLoadSample3<'_, '_, '_, '_> {
    pub fn time_s(&self) -> f64 { self.time_s }
    pub fn nodal_forces_n(&self) -> &[(usize, f64)] { self.case.nodal_forces_n() }
    pub fn resultant(&self, origin_m: [f64; 3], cx: &Cx<'_>)
        -> Result<([f64; 3], [f64; 3]), LoadError3> { self.case.resultant(origin_m, cx) }
    pub fn solve(&self, thermal: TetThermalStrainField<'_>, config: TetStaticSolveConfig,
        mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>)
        -> Result<TetLoadSolution3, LoadError3> {
        self.case.solve(self.problem, thermal, config, mesh_frame, maximum_strain, cx)
    }
}
impl TetLoadHistory3<'_> {
    /// Resolve every source at the same physical time before publishing a load.
    /// Vector frame time is never replaced by a table's clamped/wrapped time.
    /// Invalid quantities, sources, face bindings and cancellation publish no
    /// partial loads. The returned case may be passed straight to `solve` for
    /// input-bound stress/reactions; the original Cx remains active in that solve.
    /// Scalar TimeSignal evaluation retains its existing indivisible work unit;
    /// there is no new sub-kernel cancellation-latency or time-accuracy claim.
    pub fn at<'g, 'm, 'p, 's>(&self, geometry: &'g TetLoadGeometry3<'m>,
        problem: &'p TetLinearElasticProblem<'s>, time_s: f64, cx: &Cx<'_>)
        -> Result<TetLoadSample3<'g, 'm, 'p, 's>, LoadHistoryError3> {
        poll(cx)?;
        if !time_s.is_finite() || self.surfaces.len() > 1024 {
            return Err(LoadHistoryError3::Invalid("finite time and at most 1024 surface histories required"));
        }
        if problem.nodes_m != geometry.nodes || problem.tetrahedra != geometry.tets {
            return Err(LoadHistoryError3::Invalid("history and problem must share the same reference mesh"));
        }
        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        for row in self.surfaces {
            poll(cx)?;
            surfaces.push(match *row {
                SurfaceHistory3::Fixed(value) => value,
                SurfaceHistory3::Traction { faces, value } => SurfaceLoad3::Traction {
                    faces, pascals: sample(value, TRACTION_DIMS, self.frames, time_s, cx)?,
                },
                SurfaceHistory3::Pressure { faces, value } => {
                    let pressure = value.eval(time_s).map_err(|e| LoadHistoryError3::Scenario(e.to_string()))?;
                    if pressure.dims != TRACTION_DIMS {
                        return Err(LoadHistoryError3::Quantity { expected: TRACTION_DIMS, found: QuantityContract::Dimensions(pressure.dims) });
                    }
                    SurfaceLoad3::Pressure { faces, pascals: pressure.value }
                }
            });
        }
        let case = match self.volume {
            VolumeHistory3::Fixed(value) => geometry.assemble(
                &surfaces, value, self.nodal_forces_n, cx)?,
            VolumeHistory3::ForceDensity(value) => geometry.assemble(&surfaces,
                VolumeLoad3::Uniform(sample(value, FORCE_DENSITY_DIMS, self.frames, time_s, cx)?),
                self.nodal_forces_n, cx)?,
            VolumeHistory3::Acceleration(value) => {
                let acceleration = sample(value, ACCELERATION_DIMS, self.frames, time_s, cx)?;
                geometry.assemble_acceleration(&surfaces, problem.materials,
                    AccelerationField3::Uniform(acceleration), self.nodal_forces_n, cx)?
            }
            VolumeHistory3::FixedAcceleration(value) => geometry.assemble_acceleration(
                &surfaces, problem.materials, value, self.nodal_forces_n, cx)?,
        };
        poll(cx)?;
        Ok(TetLoadSample3 { case, problem, time_s })
    }
}

/// Explicit absolute support-motion histories sharing the load endpoint time.
pub mod supports;

/// Scheduled physical loads driving sparse inertial trajectories.
pub mod dynamics;
