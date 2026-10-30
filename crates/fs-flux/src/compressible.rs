//! Two-dimensional compressible finite volumes over the original `fs-gas` EOS.
//!
//! Unsplit MUSCL / SSP RK2 shares one conservative flux per Cartesian face.
//! Normal and transverse momentum and both kinetic-energy components use the
//! same source-owned HLLC/HLL closure. Stationary walls use the thermally perfect
//! one-sided shock/rarefaction solution, not a linear acoustic pressure law.
//! This opt-in inviscid engineering lane has no entropy-stability, adjoint,
//! turbulence, physical shock-thickness or experimental-validation claim.

pub use fs_gas::dynamics::{Error, FrozenGas, Reconstruction, StepControl};
use fs_gas::dynamics::{MAX_EULER_CELLS, NormalConserved, NormalFaceFlux};
use fs_qty::Time;

/// Conservative density, Cartesian momentum and full referenced total energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct State {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Momentum density in the x/y axes [kg/(m² s)].
    pub momentum: [f64; 2],
    /// `rho*(e(T)-e(Tmin)+(u²+v²)/2)` [J/m³].
    pub energy: f64,
}

impl State {
    /// Construct through the original gas source, checking caloric precision.
    pub fn from_primitive(
        gas: &FrozenGas,
        density: f64,
        velocity: [f64; 2],
        pressure: f64,
    ) -> Result<Self, Error> {
        let q = gas.normal_conserved(density, velocity[0], velocity[1], pressure)?;
        Ok(Self {
            density: q.density,
            momentum: [q.normal_momentum, q.transverse_momentum],
            energy: q.energy,
        })
    }

    /// Decode pressure, temperature and sound speed from both kinetic components.
    pub fn primitive(
        self,
        gas: &FrozenGas,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let p = gas.normal_primitive(self.normal(0), cancelled)?;
        Ok(Primitive {
            density: p.density,
            velocity: [p.velocity, self.momentum[1] / self.density],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        })
    }

    fn normal(self, axis: usize) -> NormalConserved {
        NormalConserved {
            density: self.density,
            normal_momentum: self.momentum[axis],
            transverse_momentum: self.momentum[1 - axis],
            energy: self.energy,
        }
    }
    fn array(self) -> [f64; 4] {
        [
            self.density,
            self.momentum[0],
            self.momentum[1],
            self.energy,
        ]
    }
    fn from_array(a: [f64; 4]) -> Self {
        Self {
            density: a[0],
            momentum: [a[1], a[2]],
            energy: a[3],
        }
    }
}

/// Thermodynamic state in the grid frame, from the same frozen mixture.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Primitive {
    /// Mass density [kg/m³].
    pub density: f64,
    /// Cartesian velocity [m/s].
    pub velocity: [f64; 2],
    /// Absolute pressure [Pa].
    pub pressure: f64,
    /// Absolute temperature [K].
    pub temperature: f64,
    /// Frozen-composition sound speed [m/s].
    pub sound_speed: f64,
}
impl Primitive {
    fn array(self) -> [f64; 4] {
        [
            self.density,
            self.velocity[0],
            self.velocity[1],
            self.pressure,
        ]
    }
}

/// A physical side of the rectangular domain; face order increases tangentially.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundarySide {
    /// Minimum x coordinate, with outward normal (-1, 0).
    XMin,
    /// Maximum x coordinate, with outward normal (1, 0).
    XMax,
    /// Minimum y coordinate, with outward normal (0, -1).
    YMin,
    /// Maximum y coordinate, with outward normal (0, 1).
    YMax,
}
impl BoundarySide {
    fn axis(self) -> usize {
        match self {
            Self::XMin | Self::XMax => 0,
            Self::YMin | Self::YMax => 1,
        }
    }
    fn sign(self) -> f64 {
        match self {
            Self::XMin | Self::YMin => -1.0,
            Self::XMax | Self::YMax => 1.0,
        }
    }
}
const SIDES: [BoundarySide; 4] = [
    BoundarySide::XMin,
    BoundarySide::XMax,
    BoundarySide::YMin,
    BoundarySide::YMax,
];

/// Physical treatment at a stationary exterior face.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryCondition {
    /// Zero-gradient continuation, not a nonreflecting subsonic outlet.
    Transmissive,
    /// Impermeable, adiabatic slip wall with zero mass and energy flux.
    Reflecting,
}

/// Periodic axes share exactly the same face flux at both ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AxisBoundaries {
    /// Identify the two ends of this axis.
    Periodic,
    /// Independent exterior conditions.
    Ends {
        /// Treatment at the minimum coordinate.
        lower: BoundaryCondition,
        /// Treatment at the maximum coordinate.
        upper: BoundaryCondition,
    },
}

/// Independent x and y boundary pairs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Boundaries2d {
    /// Conditions normal to x.
    pub x: AxisBoundaries,
    /// Conditions normal to y.
    pub y: AxisBoundaries,
}
impl Boundaries2d {
    fn axis(self, axis: usize) -> AxisBoundaries {
        if axis == 0 { self.x } else { self.y }
    }
}

/// Uniform Cartesian cells with a physical extrusion depth.
///
/// States use row-major order `j*nx+i`; lengths are coherent SI metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CartesianGrid {
    shape: [usize; 2],
    spacing: [f64; 2],
    depth: f64,
}
impl CartesianGrid {
    /// Admit positive, representable geometry with at most MAX_EULER_CELLS cells.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, depth: f64) -> Result<Self, Error> {
        let n = nx
            .checked_mul(ny)
            .ok_or(Error::Invalid("Cartesian cell count overflow"))?;
        if nx == 0
            || ny == 0
            || n > MAX_EULER_CELLS
            || ![
                dx,
                dy,
                depth,
                dx * dy,
                dx * depth,
                dy * depth,
                dx * dy * depth,
                nx as f64 * dx,
                ny as f64 * dy,
                n as f64 * dx * dy * depth,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(Error::Invalid("Cartesian geometry or cell budget"));
        }
        Ok(Self {
            shape: [nx, ny],
            spacing: [dx, dy],
            depth,
        })
    }
    /// Number of cells along x and y.
    pub const fn shape(&self) -> [usize; 2] {
        self.shape
    }
    /// Cell widths [m].
    pub const fn spacing(&self) -> [f64; 2] {
        self.spacing
    }
    /// Physical out-of-plane depth [m].
    pub const fn depth(&self) -> f64 {
        self.depth
    }
    /// Cell volume [m³].
    pub fn cell_volume(&self) -> f64 {
        self.spacing[0] * self.spacing[1] * self.depth
    }
    /// One physical boundary-face area [m²].
    pub fn face_area(&self, side: BoundarySide) -> f64 {
        self.spacing[1 - side.axis()] * self.depth
    }
    /// Number of faces along a domain side.
    pub fn face_count(&self, side: BoundarySide) -> usize {
        self.shape[1 - side.axis()]
    }
    /// Centre of a cell [m], relative to the caller's origin.
    pub fn cell_center(&self, index: usize) -> Option<[f64; 2]> {
        (index < self.len()).then(|| {
            [
                (index % self.shape[0]) as f64 * self.spacing[0] + 0.5 * self.spacing[0],
                (index / self.shape[0]) as f64 * self.spacing[1] + 0.5 * self.spacing[1],
            ]
        })
    }
    fn len(&self) -> usize {
        self.shape[0] * self.shape[1]
    }
    fn boundary_range(&self, side: BoundarySide) -> std::ops::Range<usize> {
        let [nx, ny] = self.shape;
        let (start, len) = match side {
            BoundarySide::XMin => (0, ny),
            BoundarySide::XMax => (ny, ny),
            BoundarySide::YMin => (2 * ny, nx),
            BoundarySide::YMax => (2 * ny + nx, nx),
        };
        start..start + len
    }
}

/// Extensive mass [kg], momentum [kg m/s] and referenced energy [J].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Totals {
    /// Mass [kg].
    pub mass: f64,
    /// Cartesian momentum [kg m/s].
    pub momentum: [f64; 2],
    /// Referenced total energy [J].
    pub energy: f64,
}
impl Totals {
    fn array(self) -> [f64; 4] {
        [self.mass, self.momentum[0], self.momentum[1], self.energy]
    }
    fn from_array(a: [f64; 4]) -> Self {
        Self {
            mass: a[0],
            momentum: [a[1], a[2]],
            energy: a[3],
        }
    }
    /// Recover original physical total energy without changing the transport gauge.
    pub fn absolute_energy(self, gas: &FrozenGas) -> f64 {
        self.energy + gas.energy_reference() * self.mass
    }
}

/// Signed outward transfer across one face during the last accepted RK step.
/// At a reflecting wall momentum is the impulse exerted on its support.
pub type BoundaryImpulse = Totals;

/// Only accepted steps publish balances and time-integrated boundary loads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Actual accepted duration [s].
    pub elapsed: Time,
    /// Largest two-direction face CFL across both stages.
    pub maximum_cfl: f64,
    /// Number of rejected, halved attempts.
    pub halvings: usize,
    /// Number of accepted-stage faces using the shared HLL fallback.
    pub hll_fallback_faces: usize,
    /// Number of accepted-stage cells with source-limited reconstruction.
    pub limited_cells: usize,
    /// An admissibility refusal switched the whole trial to first order.
    pub first_order_retry: bool,
    /// Original extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Net transfer into the domain, integrated over both stages.
    pub boundary_change: Totals,
    /// `after - before - boundary_change`.
    pub conservation_residual: Totals,
}

struct Workspace {
    primitive: Vec<Primitive>,
    slopes: Vec<[[f64; 4]; 2]>,
    x_flux: Vec<NormalFaceFlux>,
    y_flux: Vec<NormalFaceFlux>,
    stage: Vec<State>,
    candidate: Vec<State>,
    boundary_trial: Vec<BoundaryImpulse>,
}

/// Source-owned, unsplit two-dimensional Euler field on a stationary rectangle.
///
/// Workspace is allocated once. Failed retries or cancellation preserve cells,
/// time and the previously accepted wall impulse. Source-domain violations
/// reduce the step; no mass, momentum, energy or temperature is clipped.
pub struct Euler2d {
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    states: Vec<State>,
    time: f64,
    last_step: Option<StepReport>,
    boundary_impulses: Vec<BoundaryImpulse>,
    work: Workspace,
}
impl Euler2d {
    /// Admit an original gas source, fixed grid and conservative cell averages.
    pub fn new(
        gas: FrozenGas,
        grid: CartesianGrid,
        states: Vec<State>,
        boundaries: Boundaries2d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if states.len() != grid.len() {
            return Err(Error::Invalid("Cartesian state count"));
        }
        let n = states.len();
        let [nx, ny] = grid.shape;
        let mut primitive = allocated(n, Primitive::default())?;
        for (p, q) in primitive.iter_mut().zip(&states) {
            *p = q.primitive(&gas, &mut cancelled)?;
        }
        if !checked_totals(&states, grid.cell_volume(), &mut cancelled)?
            .array()
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid("unrepresentable Cartesian extensive state"));
        }
        let result = Self {
            gas,
            grid,
            boundaries,
            states,
            time: 0.0,
            last_step: None,
            boundary_impulses: allocated(2 * (nx + ny), Totals::default())?,
            work: Workspace {
                primitive,
                slopes: allocated(n, [[0.0; 4]; 2])?,
                x_flux: allocated((nx + 1) * ny, NormalFaceFlux::default())?,
                y_flux: allocated(nx * (ny + 1), NormalFaceFlux::default())?,
                stage: allocated(n, State::default())?,
                candidate: allocated(n, State::default())?,
                boundary_trial: allocated(2 * (nx + ny), Totals::default())?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
    /// Original source-bearing thermodynamic owner.
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }
    /// Admitted geometry, including physical depth.
    pub const fn grid(&self) -> &CartesianGrid {
        &self.grid
    }
    /// Current exterior conditions.
    pub const fn boundaries(&self) -> Boundaries2d {
        self.boundaries
    }
    /// Current conservative averages in row-major order.
    pub fn states(&self) -> &[State] {
        &self.states
    }
    /// Accepted physical time.
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Most recent accepted step, absent before the first step.
    pub const fn last_step(&self) -> Option<&StepReport> {
        self.last_step.as_ref()
    }
    /// Outward per-face integrated transfer from that same accepted step.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.grid.boundary_range(side)]
    }
    /// Compensated extensive state in deterministic cell order.
    pub fn totals(&self) -> Totals {
        totals(&self.states, self.grid.cell_volume())
    }

    /// Advance at most the requested duration with a combined x/y CFL budget.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("Cartesian duration, CFL or retry controls"));
        }
        let mut rate: f64 = 0.0;
        for (p, q) in self.work.primitive.iter_mut().zip(&self.states) {
            *p = q.primitive(&self.gas, &mut cancelled)?;
            rate = rate.max(
                (p.velocity[0].abs() + p.sound_speed) / self.grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / self.grid.spacing[1],
            );
        }
        let mut dt = maximum_duration.value().min(controls.cfl / rate);
        let before = checked_totals(&self.states, self.grid.cell_volume(), &mut cancelled)?;
        let mut reconstruction = controls.reconstruction;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !self.time.is_finite()
                || self.time + dt <= self.time
                || !(self.time + dt).is_finite()
            {
                return Err(Error::Invalid("unrepresentable Cartesian time increment"));
            }
            let trial = attempt(
                &self.gas,
                &self.grid,
                self.boundaries,
                &self.states,
                &mut self.work,
                dt,
                controls.cfl,
                reconstruction,
                &mut cancelled,
            );
            match trial {
                Ok(stats) => {
                    let after = checked_totals(
                        &self.work.candidate,
                        self.grid.cell_volume(),
                        &mut cancelled,
                    )?;
                    let outward = sum_arrays(
                        self.work.boundary_trial.iter().map(|v| v.array()),
                        &mut cancelled,
                    )?;
                    let boundary_change = Totals::from_array(outward.map(|x| -x));
                    let a = after.array();
                    let b = before.array();
                    let c = boundary_change.array();
                    let report = StepReport {
                        elapsed: Time::new(dt),
                        maximum_cfl: stats.rate * dt,
                        halvings,
                        hll_fallback_faces: stats.fallback,
                        limited_cells: stats.limited,
                        first_order_retry: reconstruction != controls.reconstruction,
                        before,
                        after,
                        boundary_change,
                        conservation_residual: Totals::from_array(std::array::from_fn(|k| {
                            a[k] - b[k] - c[k]
                        })),
                    };
                    if !a.into_iter().chain(b).chain(c).all(f64::is_finite) {
                        return Err(Error::Invalid(
                            "unrepresentable Cartesian extensive balance",
                        ));
                    }
                    poll(&mut cancelled)?;
                    std::mem::swap(&mut self.states, &mut self.work.candidate);
                    self.boundary_impulses
                        .copy_from_slice(&self.work.boundary_trial);
                    self.time += dt;
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => {
                    if halvings == controls.maximum_halvings {
                        return Err(Error::StepRetries {
                            attempts: halvings + 1,
                            last: Box::new(error),
                        });
                    }
                    if !matches!(error, Error::CflExceeded { .. }) {
                        reconstruction = Reconstruction::FirstOrder;
                    }
                    dt *= 0.5;
                }
            }
        }
        unreachable!("bounded step attempts return on the last attempt")
    }
}

#[derive(Default)]
struct StageStats {
    rate: f64,
    fallback: usize,
    limited: usize,
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    states: &[State],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let first = spatial(
        gas,
        grid,
        boundaries,
        states,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.x_flux,
        &mut work.y_flux,
        reconstruction,
        cancelled,
    )?;
    check_cfl(first.rate * dt, cfl)?;
    update(
        gas,
        grid,
        states,
        states,
        &work.x_flux,
        &work.y_flux,
        &mut work.stage,
        dt,
        false,
        cancelled,
    )?;
    capture_boundaries(
        grid,
        &work.x_flux,
        &work.y_flux,
        0.5 * dt,
        &mut work.boundary_trial,
        false,
        cancelled,
    )?;
    let second = spatial(
        gas,
        grid,
        boundaries,
        &work.stage,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.x_flux,
        &mut work.y_flux,
        reconstruction,
        cancelled,
    )?;
    check_cfl(second.rate * dt, cfl)?;
    update(
        gas,
        grid,
        states,
        &work.stage,
        &work.x_flux,
        &work.y_flux,
        &mut work.candidate,
        dt,
        true,
        cancelled,
    )?;
    capture_boundaries(
        grid,
        &work.x_flux,
        &work.y_flux,
        0.5 * dt,
        &mut work.boundary_trial,
        true,
        cancelled,
    )?;
    Ok(StageStats {
        rate: first.rate.max(second.rate),
        fallback: first.fallback + second.fallback,
        limited: first.limited + second.limited,
    })
}

fn check_cfl(observed: f64, maximum: f64) -> Result<(), Error> {
    if !observed.is_finite() || observed > maximum * (1.0 + 16.0 * f64::EPSILON) {
        Err(Error::CflExceeded { observed, maximum })
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    states: &[State],
    primitive: &mut [Primitive],
    slopes: &mut [[[f64; 4]; 2]],
    x_flux: &mut [NormalFaceFlux],
    y_flux: &mut [NormalFaceFlux],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let mut stats = StageStats::default();
    for (p, q) in primitive.iter_mut().zip(states) {
        *p = q.primitive(gas, &mut *cancelled)?;
    }
    let [nx, ny] = grid.shape;
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if k % 64 == 0 {
                poll(cancelled)?;
            }
            slopes[k] = [[0.0; 4]; 2];
            if reconstruction == Reconstruction::FirstOrder {
                continue;
            }
            let center = primitive[k].array();
            for axis in 0..2 {
                let lower = neighbour(primitive, grid, boundaries, i, j, axis, false);
                let upper = neighbour(primitive, grid, boundaries, i, j, axis, true);
                slopes[k][axis] =
                    std::array::from_fn(|d| mc(center[d] - lower[d], upper[d] - center[d]));
            }
            let mut theta = 1.0;
            for _ in 0..32 {
                let valid = (0..2).all(|axis| {
                    [-0.5, 0.5].into_iter().all(|sign| {
                        let p: [f64; 4] =
                            std::array::from_fn(|d| center[d] + sign * theta * slopes[k][axis][d]);
                        let t = p[3] / (gas.gas_constant() * p[0]);
                        p.into_iter().all(f64::is_finite)
                            && positive(p[0])
                            && positive(p[3])
                            && t >= gas.domain().minimum_temperature
                            && t <= gas.domain().maximum_temperature
                    })
                });
                if valid {
                    break;
                }
                theta *= 0.5;
            }
            if theta < 1.0 {
                stats.limited += 1;
                if theta <= 2.0_f64.powi(-32) {
                    theta = 0.0;
                }
                for axis in &mut slopes[k] {
                    for s in axis {
                        *s *= theta;
                    }
                }
            }
        }
    }
    for axis in 0..2 {
        let count = grid.shape[axis];
        let lines = grid.shape[1 - axis];
        for line in 0..lines {
            for face in 0..=count {
                if face % 64 == 0 {
                    poll(cancelled)?;
                }
                let pair = boundaries.axis(axis);
                if face == count && pair == AxisBoundaries::Periodic {
                    let first = face_index(grid, axis, 0, line);
                    let last = face_index(grid, axis, count, line);
                    if axis == 0 {
                        x_flux[last] = x_flux[first];
                    } else {
                        y_flux[last] = y_flux[first];
                    }
                    continue;
                }
                let cell_index = |position| {
                    if axis == 0 {
                        line * nx + position
                    } else {
                        position * nx + line
                    }
                };
                let reconstruct = |k: usize, sign: f64| -> Result<NormalConserved, Error> {
                    if reconstruction == Reconstruction::FirstOrder {
                        return Ok(states[k].normal(axis));
                    }
                    let center = primitive[k].array();
                    let p: [f64; 4] =
                        std::array::from_fn(|d| center[d] + sign * slopes[k][axis][d]);
                    gas.normal_conserved(p[0], p[axis + 1], p[2 - axis], p[3])
                };
                let flux = if face > 0 && face < count {
                    gas.hllc_normal_flux(
                        reconstruct(cell_index(face - 1), 0.5)?,
                        reconstruct(cell_index(face), -0.5)?,
                        &mut *cancelled,
                    )?
                } else if pair == AxisBoundaries::Periodic {
                    gas.hllc_normal_flux(
                        reconstruct(cell_index(count - 1), 0.5)?,
                        reconstruct(cell_index(0), -0.5)?,
                        &mut *cancelled,
                    )?
                } else {
                    let AxisBoundaries::Ends { lower, upper } = pair else {
                        unreachable!()
                    };
                    let at_lower = face == 0;
                    let boundary = if at_lower { lower } else { upper };
                    let sign = if at_lower { -1.0 } else { 1.0 };
                    let q =
                        reconstruct(cell_index(if at_lower { 0 } else { count - 1 }), 0.5 * sign)?;
                    match boundary {
                        BoundaryCondition::Reflecting => {
                            gas.stationary_wall_normal_flux(q, sign, &mut *cancelled)?
                        }
                        BoundaryCondition::Transmissive => {
                            gas.hllc_normal_flux(q, q, &mut *cancelled)?
                        }
                    }
                };
                stats.fallback += usize::from(flux.used_hll_fallback);
                let index = face_index(grid, axis, face, line);
                if axis == 0 {
                    x_flux[index] = flux;
                } else {
                    y_flux[index] = flux;
                }
            }
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let x = j * (nx + 1) + i;
            let y = j * nx + i;
            let rate = x_flux[x].signal_speed.max(x_flux[x + 1].signal_speed) / grid.spacing[0]
                + y_flux[y].signal_speed.max(y_flux[y + nx].signal_speed) / grid.spacing[1];
            stats.rate = stats.rate.max(rate);
        }
    }
    Ok(stats)
}

fn face_index(grid: &CartesianGrid, axis: usize, face: usize, line: usize) -> usize {
    if axis == 0 {
        line * (grid.shape[0] + 1) + face
    } else {
        face * grid.shape[0] + line
    }
}

fn neighbour(
    p: &[Primitive],
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    i: usize,
    j: usize,
    axis: usize,
    upper: bool,
) -> [f64; 4] {
    let index = if axis == 0 { i } else { j };
    let length = grid.shape[axis];
    let edge = if upper {
        index + 1 == length
    } else {
        index == 0
    };
    let mut q = p[j * grid.shape[0] + i].array();
    let next = if !edge {
        Some(if upper { index + 1 } else { index - 1 })
    } else if boundaries.axis(axis) == AxisBoundaries::Periodic {
        Some(if upper { 0 } else { length - 1 })
    } else {
        None
    };
    if let Some(next) = next {
        q = p[if axis == 0 {
            j * grid.shape[0] + next
        } else {
            next * grid.shape[0] + i
        }]
        .array();
    } else if let AxisBoundaries::Ends { lower, upper: high } = boundaries.axis(axis) {
        if (if upper { high } else { lower }) == BoundaryCondition::Reflecting {
            q[axis + 1] = -q[axis + 1];
        }
    }
    q
}

fn mc(left: f64, right: f64) -> f64 {
    if left > 0.0 && right > 0.0 {
        (0.5 * left + 0.5 * right).min(2.0 * left).min(2.0 * right)
    } else if left < 0.0 && right < 0.0 {
        (0.5 * left + 0.5 * right).max(2.0 * left).max(2.0 * right)
    } else {
        0.0
    }
}

fn flux_array(f: NormalFaceFlux, axis: usize) -> [f64; 4] {
    if axis == 0 {
        [f.mass, f.normal_momentum, f.transverse_momentum, f.energy]
    } else {
        [f.mass, f.transverse_momentum, f.normal_momentum, f.energy]
    }
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    original: &[State],
    stage: &[State],
    x_flux: &[NormalFaceFlux],
    y_flux: &[NormalFaceFlux],
    output: &mut [State],
    dt: f64,
    final_stage: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let [nx, ny] = grid.shape;
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let x = j * (nx + 1) + i;
            let left = flux_array(x_flux[x], 0);
            let right = flux_array(x_flux[x + 1], 0);
            let lower = flux_array(y_flux[k], 1);
            let upper = flux_array(y_flux[k + nx], 1);
            let old = original[k].array();
            let current = stage[k].array();
            output[k] = State::from_array(std::array::from_fn(|d| {
                let delta = dt / grid.spacing[0] * (left[d] - right[d])
                    + dt / grid.spacing[1] * (lower[d] - upper[d]);
                if final_stage {
                    0.5 * old[d] + 0.5 * current[d] + 0.5 * delta
                } else {
                    current[d] + delta
                }
            }));
            output[k].primitive(gas, &mut *cancelled)?;
        }
    }
    Ok(())
}

fn capture_boundaries(
    grid: &CartesianGrid,
    x_flux: &[NormalFaceFlux],
    y_flux: &[NormalFaceFlux],
    duration: f64,
    output: &mut [BoundaryImpulse],
    add: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for side in SIDES {
        let axis = side.axis();
        let face = if side.sign() < 0.0 {
            0
        } else {
            grid.shape[axis]
        };
        let scale = side.sign() * duration * grid.face_area(side);
        for (line, target) in output[grid.boundary_range(side)].iter_mut().enumerate() {
            if line % 64 == 0 {
                poll(cancelled)?;
            }
            let index = face_index(grid, axis, face, line);
            let flux = flux_array(
                if axis == 0 {
                    x_flux[index]
                } else {
                    y_flux[index]
                },
                axis,
            );
            let previous = if add { target.array() } else { [0.0; 4] };
            *target = Totals::from_array(std::array::from_fn(|d| previous[d] + scale * flux[d]));
        }
    }
    Ok(())
}

fn sum_arrays(
    values: impl Iterator<Item = [f64; 4]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 4], Error> {
    let mut sum = [0.0; 4];
    let mut correction = [0.0; 4];
    for (index, a) in values.enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        for d in 0..4 {
            let y = a[d] - correction[d];
            let t = sum[d] + y;
            correction[d] = (t - sum[d]) - y;
            sum[d] = t;
        }
    }
    Ok(sum)
}
fn totals(states: &[State], volume: f64) -> Totals {
    checked_totals(states, volume, &mut || false).expect("infallible cancellation callback")
}
fn checked_totals(
    states: &[State],
    volume: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Totals, Error> {
    Ok(Totals::from_array(sum_arrays(
        states.iter().map(|q| q.array().map(|v| v * volume)),
        cancelled,
    )?))
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn allocated<T: Clone>(len: usize, value: T) -> Result<Vec<T>, Error> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| Error::Allocation)?;
    values.resize(len, value);
    Ok(values)
}

#[cfg(test)]
mod tests;
