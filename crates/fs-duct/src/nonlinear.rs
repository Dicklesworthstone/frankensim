//! Bidirectional finite-amplitude gas waves on the existing duct geometry.
//!
//! This image evolves mass, axial momentum and **total** gas energy through
//! [`fs_gas`]'s conservative quasi-one-dimensional Euler solver. Compression
//! steepens, shocks heat the gas, and reflected waves travel in both directions.
//! The caller supplies the original NASA-9 source through [`FrozenGas`], which
//! retains its temperature-dependent caloric law throughout propagation.
//!
//! This is an inviscid finite-volume engineering lane, also suitable as the
//! subcell building block for nonlinear DG. It does not inherit the linear DG
//! exponential, its convergence certificate, Zwikker--Kosten memory, or a
//! fully discrete entropy-stability theorem. Smooth-area plane-wave modeling
//! does not resolve transverse modes, abrupt-area junction losses, tone holes,
//! moving/compliant walls, or physical shock thickness.

use core::fmt;

use crate::{Duct, DuctError, RadiusProfile, Segment};
pub use fs_gas::dynamics::{
    Boundaries, Boundary, CaloricDomain, EndWave, FaceFlux, FrozenGas, PortStage, Primitive,
    StepControl, StepReport,
};
use fs_gas::dynamics::{Conserved, Error as GasError, Euler1d, Mesh1d};
use fs_qty::Time;

/// Source-owned molecular stress, heat and species transport in a nonlinear bore.
pub mod transport;

/// A sounding bore: radius steps, Zwikker-Kosten wall memory and a radiating end.
pub mod bore;

/// Explicit allocation bound for this geometric image; excess requests refuse.
pub const MAX_NONLINEAR_DUCT_CELLS: usize = fs_gas::dynamics::MAX_EULER_CELLS;
/// Bound on composite four-point Gaussian refinement of one smooth profile span.
pub const MAX_PROFILE_VOLUME_SUBDIVISIONS: usize = 256;

/// Geometric or gas-dynamics refusal, retaining the originating error.
#[derive(Debug, Clone, PartialEq)]
pub enum NonlinearDuctError {
    /// Existing duct primitive or profile admission failed.
    Geometry(DuctError),
    /// Gas state, solve budget, or cancellation refused.
    Gas(GasError),
    /// Geometry or requested discretization is outside this image's domain.
    UnsupportedGeometry(&'static str),
}

impl fmt::Display for NonlinearDuctError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Geometry(error) => write!(f, "FS-DUCT-NONLINEAR-GEOMETRY: {error}"),
            Self::Gas(error) => write!(f, "FS-DUCT-NONLINEAR-GAS: {error}"),
            Self::UnsupportedGeometry(reason) => {
                write!(f, "FS-DUCT-NONLINEAR-UNSUPPORTED: {reason}")
            }
        }
    }
}
impl std::error::Error for NonlinearDuctError {}
impl From<DuctError> for NonlinearDuctError {
    fn from(error: DuctError) -> Self {
        Self::Geometry(error)
    }
}
impl From<GasError> for NonlinearDuctError {
    fn from(error: GasError) -> Self {
        Self::Gas(error)
    }
}

/// Primitive initial data sampled at each axial cell midpoint.
///
/// These values initialize piecewise-constant cell averages. An initial shock
/// should coincide with a cell face when an exact discontinuity location matters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InitialGasState {
    /// Mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Signed axial velocity, positive from inlet to outlet [m/s].
    pub velocity_m_s: f64,
    /// Absolute pressure [Pa], not acoustic overpressure.
    pub pressure_pa: f64,
}

/// Cell-center gas observable, using the original area at its axial midpoint.
///
/// Volume velocity is a local diagnostic, not the accepted Riemann face flux.
/// Conservative exchange with another subsystem must use the solver's integrated
/// boundary exchange; a sampled pressure/flow pair is not an energy audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearDuctSample {
    /// Axial midpoint [m].
    pub axial_position_m: f64,
    /// Original geometric cross-sectional area at the midpoint [m²].
    pub area_m2: f64,
    /// Absolute gas pressure [Pa].
    pub pressure_pa: f64,
    /// Density [kg/m³].
    pub density_kg_m3: f64,
    /// Absolute temperature [K].
    pub temperature_k: f64,
    /// Signed axial velocity [m/s].
    pub velocity_m_s: f64,
    /// Frozen-composition sound speed from the evolving gas law [m/s].
    pub sound_speed_m_s: f64,
    /// Signed local volume velocity, A(x) u [m³/s].
    pub volume_velocity_m3_s: f64,
    /// Signed local mass flow, rho A(x) u [kg/s].
    pub mass_flow_kg_s: f64,
}

/// A nonlinear propagation image of an ordinary [`Duct`].
///
/// Geometry is prepared once; time stepping delegates the one shared gas solver.
/// The exact frozen composition, source cards and energy reference remain
/// available through [`Self::solver`].
pub struct NonlinearDuct {
    geometry: Duct,
    solver: Euler1d,
    centers_m: Vec<f64>,
    center_areas_m2: Vec<f64>,
}

impl NonlinearDuct {
    /// Build a conservative mesh, admitting one positive cell count per segment.
    ///
    /// Every segment boundary is a mesh face. Cylinder and cone volumes use
    /// closed forms; profile volumes use bounded converged Gaussian quadrature,
    /// split at original spline knots. Face and midpoint areas always evaluate
    /// the original geometry, never a radius interpolant over the solver mesh.
    ///
    /// Adjacent endpoint radii must agree to roundoff. Compact side branches and
    /// abrupt area jumps require a separately modeled junction and refuse here.
    /// No wall-loss or gas-composition default is selected implicitly.
    pub fn new(
        duct: &Duct,
        gas: FrozenGas,
        cells_per_segment: &[usize],
        boundaries: Boundaries,
        initial: impl FnMut(f64) -> Result<InitialGasState, GasError>,
    ) -> Result<Self, NonlinearDuctError> {
        Self::new_controlled(duct, gas, cells_per_segment, boundaries, initial, || false)
    }

    /// Cancellable construction, polling before each cell and each quadrature
    /// refinement. A failed construction returns no partially prepared image.
    pub fn new_controlled(
        duct: &Duct,
        gas: FrozenGas,
        cells_per_segment: &[usize],
        boundaries: Boundaries,
        mut initial: impl FnMut(f64) -> Result<InitialGasState, GasError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NonlinearDuctError> {
        poll(&mut cancelled)?;
        if duct.segments.is_empty() {
            return Err(DuctError::EmptyDuct.into());
        }
        if cells_per_segment.len() != duct.segments.len() || cells_per_segment.contains(&0) {
            return Err(unsupported(
                "one positive cell count is required per segment",
            ));
        }
        let count = cells_per_segment
            .iter()
            .try_fold(0usize, |total, &n| {
                total
                    .checked_add(n)
                    .filter(|&sum| sum <= MAX_NONLINEAR_DUCT_CELLS)
            })
            .ok_or_else(|| unsupported("cell count exceeds MAX_NONLINEAR_DUCT_CELLS"))?;
        // Admit topology before allocation and before invoking caller initial data.
        for (index, segment) in duct.segments.iter().enumerate() {
            segment.validate()?;
            if matches!(segment, Segment::ToneHole { .. }) {
                return Err(unsupported(
                    "tone holes require a conservative gas junction",
                ));
            }
            if index > 0 {
                let previous = duct.segments[index - 1].outlet_radius();
                let next = segment.inlet_radius();
                if (previous - next).abs() > 64.0 * f64::EPSILON * previous.max(next) {
                    return Err(unsupported(
                        "abrupt radius changes require a conservative gas junction",
                    ));
                }
            }
        }
        let mut positions = reserved(count + 1)?;
        let mut areas = reserved(count + 1)?;
        let mut volumes = reserved(count)?;
        let mut centers_m = reserved(count)?;
        let mut center_areas_m2 = reserved(count)?;
        let mut cells = reserved(count)?;
        let mut origin = 0.0;
        positions.push(origin);
        areas.push(area(duct.segments[0].inlet_radius())?);
        for (segment, &n) in duct.segments.iter().zip(cells_per_segment) {
            let length = segment.axial_length();
            for i in 0..n {
                poll(&mut cancelled)?;
                let left = length * (i as f64 / n as f64);
                let right = if i + 1 == n {
                    length
                } else {
                    length * ((i + 1) as f64 / n as f64)
                };
                let midpoint = left.midpoint(right);
                let global_right = origin + right;
                let global_midpoint = origin + midpoint;
                if !global_right.is_finite()
                    || global_right <= *positions.last().expect("initial face")
                    || !(global_midpoint > *positions.last().expect("initial face")
                        && global_midpoint < global_right)
                {
                    return Err(unsupported("axial mesh coordinates are not representable"));
                }
                let left_radius = radius(segment, left)?;
                let right_radius = radius(segment, right)?;
                let volume = match segment {
                    Segment::Cylinder { .. } => area(left_radius)? * (right - left),
                    Segment::Cone { .. } => {
                        core::f64::consts::PI
                            * ((right - left) / 3.0)
                            * (left_radius * left_radius
                                + left_radius * right_radius
                                + right_radius * right_radius)
                    }
                    Segment::Profile { profile, .. } => {
                        profile_volume(profile, left, right, &mut cancelled)?
                    }
                    Segment::ToneHole { .. } => unreachable!("topology admitted above"),
                };
                if !positive(volume) {
                    return Err(unsupported("cell volume is not positive and finite"));
                }
                let state = initial(global_midpoint)?;
                cells.push(gas.conserved(
                    state.density_kg_m3,
                    state.velocity_m_s,
                    state.pressure_pa,
                )?);
                positions.push(global_right);
                areas.push(area(right_radius)?);
                volumes.push(volume);
                centers_m.push(global_midpoint);
                center_areas_m2.push(area(radius(segment, midpoint)?)?);
            }
            origin += length;
        }
        poll(&mut cancelled)?;
        let mesh = Mesh1d::new(positions, areas, volumes)?;
        let solver = Euler1d::new_on_mesh(gas, mesh, cells, boundaries, &mut cancelled)?;
        Ok(Self {
            geometry: duct.clone(),
            solver,
            centers_m,
            center_areas_m2,
        })
    }

    /// Original immutable duct geometry.
    #[must_use]
    pub const fn geometry(&self) -> &Duct {
        &self.geometry
    }

    /// Shared conservative state, mesh, gas law, and time.
    #[must_use]
    pub const fn solver(&self) -> &Euler1d {
        &self.solver
    }

    /// Transfer the original mesh, source and accepted field into a gas network.
    /// No geometric resampling or thermodynamic reconstruction is performed.
    pub fn into_solver(self) -> Euler1d {
        self.solver
    }

    /// Accepted physical time.
    #[must_use]
    pub fn time(&self) -> Time {
        self.solver.time()
    }

    /// Accept one conservative step no longer than the requested duration.
    ///
    /// The source-owned [`Euler1d`] selects the admissible CFL-limited duration.
    /// Refusal or cancellation preserves this step's original state and time.
    /// Boundaries are physical gas boundaries selected at construction;
    /// frequency-domain radiation impedances are not substituted.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, NonlinearDuctError> {
        Ok(self.solver.step(maximum_duration, controls, cancelled)?)
    }

    /// Evaluate one cell without allocating a complete waveform snapshot.
    pub fn sample(&self, cell: usize) -> Result<NonlinearDuctSample, NonlinearDuctError> {
        let conserved: Conserved = *self
            .solver
            .states()
            .get(cell)
            .ok_or_else(|| unsupported("sample cell is outside the duct mesh"))?;
        let primitive = self.solver.gas().primitive(conserved, || false)?;
        let area_m2 = self.center_areas_m2[cell];
        let volume_velocity_m3_s = area_m2 * primitive.velocity;
        let mass_flow_kg_s = primitive.density * volume_velocity_m3_s;
        if !volume_velocity_m3_s.is_finite() || !mass_flow_kg_s.is_finite() {
            return Err(GasError::Invalid("duct volume or mass flow is unrepresentable").into());
        }
        Ok(NonlinearDuctSample {
            axial_position_m: self.centers_m[cell],
            area_m2,
            pressure_pa: primitive.pressure,
            density_kg_m3: primitive.density,
            temperature_k: primitive.temperature,
            velocity_m_s: primitive.velocity,
            sound_speed_m_s: primitive.sound_speed,
            volume_velocity_m3_s,
            mass_flow_kg_s,
        })
    }

    /// Cell-center pressure/flow snapshot in increasing axial order.
    pub fn samples(&self) -> Result<Vec<NonlinearDuctSample>, NonlinearDuctError> {
        let mut result = reserved(self.centers_m.len())?;
        for cell in 0..self.centers_m.len() {
            result.push(self.sample(cell)?);
        }
        Ok(result)
    }
}

fn reserved<T>(capacity: usize) -> Result<Vec<T>, NonlinearDuctError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| GasError::Allocation)?;
    Ok(result)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn unsupported(reason: &'static str) -> NonlinearDuctError {
    NonlinearDuctError::UnsupportedGeometry(reason)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), NonlinearDuctError> {
    if cancelled() {
        Err(GasError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn area(radius: f64) -> Result<f64, NonlinearDuctError> {
    let value = core::f64::consts::PI * radius * radius;
    if positive(value) {
        Ok(value)
    } else {
        Err(unsupported("bore area is not positive and finite"))
    }
}
fn radius(segment: &Segment, x: f64) -> Result<f64, NonlinearDuctError> {
    Ok(match segment {
        Segment::Cylinder { radius, .. } => *radius,
        Segment::Cone {
            inlet_radius,
            outlet_radius,
            length,
        } => {
            if x == 0.0 {
                *inlet_radius
            } else if x == *length {
                *outlet_radius
            } else {
                let t = x / length;
                (1.0 - t) * inlet_radius + t * outlet_radius
            }
        }
        Segment::Profile { profile, .. } => profile.radius(x)?,
        Segment::ToneHole { .. } => return Err(unsupported("tone hole has no axial cell")),
    })
}

// Four Gauss points integrate r(x)^2 exactly on every cubic spline span.
// Exponential and power profiles refine the same positive quadrature to a
// declared repeatability tolerance. This is a numerical estimate, not an enclosure.
fn profile_volume(
    profile: &RadiusProfile,
    left: f64,
    right: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, NonlinearDuctError> {
    let mut lower = left;
    let mut total = 0.0;
    if let Some(stations) = profile.stations() {
        let origin = stations[0].0;
        let first = stations.partition_point(|&(x, _)| x - origin <= left);
        for &(x, _) in &stations[first..] {
            let knot = x - origin;
            if knot >= right {
                break;
            }
            total += profile_span_volume(profile, lower, knot, cancelled)?;
            lower = knot;
        }
    }
    total += profile_span_volume(profile, lower, right, cancelled)?;
    if positive(total) {
        Ok(total)
    } else {
        Err(unsupported("profile volume is unrepresentable"))
    }
}
fn profile_span_volume(
    profile: &RadiusProfile,
    left: f64,
    right: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, NonlinearDuctError> {
    const NODES: [f64; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const WEIGHTS: [f64; 4] = [
        0.347_854_845_137_453_85,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_85,
    ];
    let mut previous: f64 = 0.0;
    let mut subdivisions = 1;
    loop {
        poll(cancelled)?;
        let mut integral = 0.0;
        for i in 0..subdivisions {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let a = left + (right - left) * (i as f64 / subdivisions as f64);
            let b = if i + 1 == subdivisions {
                right
            } else {
                left + (right - left) * ((i + 1) as f64 / subdivisions as f64)
            };
            if b <= a {
                return Err(unsupported(
                    "profile volume quadrature subinterval is unrepresentable",
                ));
            }
            let midpoint = a.midpoint(b);
            let half = (b - a) * 0.5;
            for (&node, &weight) in NODES.iter().zip(&WEIGHTS) {
                integral += half * weight * area(profile.radius(midpoint + half * node)?)?;
            }
        }
        if !positive(integral) {
            return Err(unsupported("profile volume quadrature is unrepresentable"));
        }
        if subdivisions > 1
            && (integral - previous).abs() <= 64.0 * f64::EPSILON * integral.max(previous)
        {
            return Ok(integral);
        }
        if subdivisions == MAX_PROFILE_VOLUME_SUBDIVISIONS {
            return Err(unsupported("profile volume quadrature did not converge"));
        }
        previous = integral;
        subdivisions *= 2;
    }
}

#[cfg(test)]
#[path = "nonlinear/tests.rs"]
mod tests;
