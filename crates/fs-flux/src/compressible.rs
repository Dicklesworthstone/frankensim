//! Cartesian compressible finite volumes over the original `fs-gas` EOS.
//!
//! Unsplit MUSCL / SSP RK2 shares one conservative flux per Cartesian face.
//! Normal and transverse momentum and both kinetic-energy components use the
//! same source-owned HLLC/HLL closure. Stationary walls use the thermally perfect
//! one-sided shock/rarefaction solution, not a linear acoustic pressure law.
//! This opt-in inviscid engineering lane has no entropy-stability, adjoint,
//! turbulence, physical shock-thickness or experimental-validation claim.
//! [`spatial`] extends the same closure and conservative step semantics to
//! three-dimensional volumes with both transverse momentum components.

pub use fs_gas::dynamics::{Error, FrozenGas, Reconstruction, StepControl};
use fs_gas::dynamics::{MAX_EULER_CELLS, NormalConserved, NormalFaceFlux, NormalPrimitive};
use fs_qty::Time;

mod cartesian;
/// Evolving canonical species in the same conservative Cartesian discretization.
pub mod mixture;
/// Original-source Fourier heat transport and complete planar Newtonian stress.
pub mod molecular;
/// Volumetric shock propagation with three conservative momentum components.
pub mod spatial;
use cartesian::StageStats;

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
    fn component_mut(&mut self, component: usize) -> &mut f64 {
        match component {
            0 => &mut self.mass,
            1 => &mut self.momentum[0],
            2 => &mut self.momentum[1],
            3 => &mut self.energy,
            _ => unreachable!("four-component frozen state"),
        }
    }
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
    slopes: Vec<f64>,
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
                slopes: allocated(8 * n, 0.0)?,
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
        let before = checked_totals(&self.states, self.grid.cell_volume(), &mut cancelled)?;
        let accepted = cartesian::advance(
            maximum_duration.value(),
            self.time,
            rate,
            controls,
            &mut cancelled,
            |dt, reconstruction, cancelled| {
                attempt(
                    &self.gas,
                    &self.grid,
                    self.boundaries,
                    &self.states,
                    &mut self.work,
                    dt,
                    controls.cfl,
                    reconstruction,
                    cancelled,
                )
            },
        )?;
        let cartesian::AcceptedStep {
            dt,
            halvings,
            reconstruction,
            stats,
        } = accepted;
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
            conservation_residual: Totals::from_array(std::array::from_fn(|k| a[k] - b[k] - c[k])),
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
        Ok(report)
    }
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
    cartesian::rk2(cancelled, |final_stage, cancelled| {
        let stats = spatial(
            gas,
            grid,
            boundaries,
            if final_stage { &work.stage } else { states },
            &mut work.primitive,
            &mut work.slopes,
            &mut work.x_flux,
            &mut work.y_flux,
            reconstruction,
            cancelled,
        )?;
        cartesian::check_cfl::<Error>(stats.rate * dt, cfl)?;
        if final_stage {
            update(
                gas,
                grid,
                states,
                &work.stage,
                &work.x_flux,
                &work.y_flux,
                &mut work.candidate,
                dt,
                final_stage,
                cancelled,
            )?;
        } else {
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
        }
        capture_boundaries(
            grid,
            &work.x_flux,
            &work.y_flux,
            0.5 * dt,
            &mut work.boundary_trial,
            final_stage,
            cancelled,
        )?;
        Ok(stats)
    })
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    states: &[State],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    x_flux: &mut [NormalFaceFlux],
    y_flux: &mut [NormalFaceFlux],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let mut stats = StageStats::default();
    for (p, q) in primitive.iter_mut().zip(states) {
        *p = q.primitive(gas, &mut *cancelled)?;
    }
    stats.limited = cartesian::slopes::<Error>(
        grid,
        boundaries,
        4,
        [1, 2],
        reconstruction,
        |k, d| primitive[k].array()[d],
        slopes,
        |k, _, sign, slope| {
            let center = primitive[k].array();
            let p: [f64; 4] = std::array::from_fn(|d| center[d] + sign * slope[d]);
            let t = p[3] / (gas.gas_constant() * p[0]);
            p.into_iter().all(f64::is_finite)
                && positive(p[0])
                && positive(p[3])
                && t >= gas.domain().minimum_temperature
                && t <= gas.domain().maximum_temperature
        },
        cancelled,
    )?;
    cartesian::faces::<Error, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            let flux = if axis == 0 {
                &mut *x_flux
            } else {
                &mut *y_flux
            };
            flux[face.index] = flux[first];
            return Ok(());
        }
        let reconstruct = |k: usize, sign: f64| -> NormalPrimitive {
            let center = primitive[k].array();
            let p: [f64; 4] =
                std::array::from_fn(|d| center[d] + sign * slopes[(2 * k + axis) * 4 + d]);
            NormalPrimitive {
                density: p[0],
                normal_velocity: p[axis + 1],
                transverse_velocity: p[2 - axis],
                pressure: p[3],
            }
        };
        let first_order = reconstruction == Reconstruction::FirstOrder;
        let flux = if let cartesian::FaceKind::Interior { left, right } = face.kind {
            if first_order {
                gas.hllc_normal_flux(
                    states[left].normal(axis),
                    states[right].normal(axis),
                    &mut *cancelled,
                )?
            } else {
                gas.hllc_normal_flux_from_primitive(
                    reconstruct(left, 0.5),
                    reconstruct(right, -0.5),
                    &mut *cancelled,
                )?
            }
        } else {
            let cartesian::FaceKind::Boundary {
                cell: k,
                sign,
                condition: boundary,
            } = face.kind
            else {
                unreachable!()
            };
            match (first_order, boundary) {
                (true, BoundaryCondition::Reflecting) => {
                    gas.stationary_wall_normal_flux(states[k].normal(axis), sign, &mut *cancelled)?
                }
                (true, BoundaryCondition::Transmissive) => {
                    let q = states[k].normal(axis);
                    gas.hllc_normal_flux(q, q, &mut *cancelled)?
                }
                (false, BoundaryCondition::Reflecting) => gas
                    .stationary_wall_normal_flux_from_primitive(
                        reconstruct(k, 0.5 * sign),
                        sign,
                        &mut *cancelled,
                    )?,
                (false, BoundaryCondition::Transmissive) => {
                    let p = reconstruct(k, 0.5 * sign);
                    gas.hllc_normal_flux_from_primitive(p, p, &mut *cancelled)?
                }
            }
        };
        stats.fallback += usize::from(flux.used_hll_fallback);
        if axis == 0 {
            x_flux[face.index] = flux;
        } else {
            y_flux[face.index] = flux;
        }
        Ok(())
    })?;
    stats.rate = cartesian::rate::<Error>(
        grid,
        |axis, index| {
            if axis == 0 {
                x_flux[index].signal_speed
            } else {
                y_flux[index].signal_speed
            }
        },
        cancelled,
    )?;
    Ok(stats)
}

fn face_index(grid: &CartesianGrid, axis: usize, face: usize, line: usize) -> usize {
    cartesian::face_index(grid, axis, face, line)
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
    let (index, reflected) =
        cartesian::neighbour(grid, boundaries, j * grid.shape[0] + i, axis, upper);
    let mut q = p[index].array();
    if reflected {
        q[axis + 1] = -q[axis + 1];
    }
    q
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
    cartesian::update::<4, Error, _>(
        grid,
        4,
        dt,
        final_stage,
        |k, d| original[k].array()[d],
        |k, d| stage[k].array()[d],
        |axis, face, d| {
            flux_array(
                if axis == 0 {
                    x_flux[face]
                } else {
                    y_flux[face]
                },
                axis,
            )[d]
        },
        cancelled,
        |k, row, cancelled| {
            output[k] = State::from_array([row[0], row[1], row[2], row[3]]);
            output[k].primitive(gas, &mut *cancelled)?;
            Ok(())
        },
    )
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
    cartesian::capture::<Error>(
        grid,
        4,
        duration,
        add,
        |axis, face, d| {
            flux_array(
                if axis == 0 {
                    x_flux[face]
                } else {
                    y_flux[face]
                },
                axis,
            )[d]
        },
        |boundary, d, delta, add| {
            let target = output[boundary].component_mut(d);
            *target = (if add { *target } else { 0.0 }) + delta;
        },
        cancelled,
    )
}

fn sum_arrays<const N: usize>(
    values: impl Iterator<Item = [f64; N]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; N], Error> {
    let mut sum = [0.0; N];
    let mut correction = [0.0; N];
    for (index, a) in values.enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        for d in 0..N {
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
