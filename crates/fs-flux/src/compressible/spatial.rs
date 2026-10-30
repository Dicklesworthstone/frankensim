//! Three-dimensional, unsplit compressible finite volumes.
//!
//! The five conserved components share the original `fs-gas` face closure and
//! NASA9 source. MUSCL reconstruction and SSP RK2 use a combined three-axis
//! face CFL, including HLL fallback speeds. This stationary Cartesian Euler
//! lane is Estimated: it does not certify entropy stability, shock thickness,
//! turbulence, reactions, moving geometry or a shock adjoint.

use super::cartesian::{StageStats, check_cfl, mc};
pub use super::{AxisBoundaries, BoundaryCondition, Error, FrozenGas, Reconstruction, StepControl};
use super::{allocated, poll, positive, sum_arrays};
use fs_gas::dynamics::{MAX_EULER_CELLS, SpatialConserved, SpatialFaceFlux, SpatialPrimitive};
use fs_qty::Time;

mod evolve;
/// Source-owned molecular stress and heat transport for the volumetric field.
pub mod molecular;

/// Conservative averages in coherent SI units and the source's energy gauge.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct State {
    /// Density [kg/m³].
    pub density: f64,
    /// x/y/z momentum density [kg/(m² s)].
    pub momentum: [f64; 3],
    /// `rho*(e(T)-e(Tmin)+|velocity|²/2)` [J/m³].
    pub energy: f64,
}
impl State {
    /// Construct all three kinetic components through the original gas source.
    pub fn from_primitive(
        gas: &FrozenGas,
        density: f64,
        velocity: [f64; 3],
        pressure: f64,
    ) -> Result<Self, Error> {
        let q =
            gas.spatial_conserved(density, velocity[0], [velocity[1], velocity[2]], pressure)?;
        Ok(Self {
            density: q.density,
            momentum: [
                q.normal_momentum,
                q.transverse_momentum[0],
                q.transverse_momentum[1],
            ],
            energy: q.energy,
        })
    }
    /// Recover the exact admitted thermodynamics, retaining all kinetic energy.
    pub fn primitive(
        self,
        gas: &FrozenGas,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let p = gas.spatial_primitive(self.normal(0), cancelled)?;
        Ok(Primitive {
            density: p.density,
            velocity: [
                p.velocity,
                self.momentum[1] / self.density,
                self.momentum[2] / self.density,
            ],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        })
    }
    fn normal(self, axis: usize) -> SpatialConserved {
        let [a, b] = tangents(axis);
        SpatialConserved {
            density: self.density,
            normal_momentum: self.momentum[axis],
            transverse_momentum: [self.momentum[a], self.momentum[b]],
            energy: self.energy,
        }
    }
    fn array(self) -> [f64; 5] {
        [
            self.density,
            self.momentum[0],
            self.momentum[1],
            self.momentum[2],
            self.energy,
        ]
    }
    fn from_array(a: [f64; 5]) -> Self {
        Self {
            density: a[0],
            momentum: [a[1], a[2], a[3]],
            energy: a[4],
        }
    }
}

/// Decoded state; caloric quantities are always computed by the original EOS.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Primitive {
    /// Density [kg/m³].
    pub density: f64,
    /// x/y/z velocity [m/s].
    pub velocity: [f64; 3],
    /// Absolute pressure [Pa].
    pub pressure: f64,
    /// Temperature [K].
    pub temperature: f64,
    /// Frozen-composition sound speed [m/s].
    pub sound_speed: f64,
}
impl Primitive {
    fn array(self) -> [f64; 5] {
        [
            self.density,
            self.velocity[0],
            self.velocity[1],
            self.velocity[2],
            self.pressure,
        ]
    }
}

/// One exterior side. Face order uses the remaining axes in increasing order:
/// y then z on x sides, x then z on y sides, and x then y on z sides.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundarySide {
    /// Minimum x, outward normal (-1,0,0).
    XMin,
    /// Maximum x, outward normal (1,0,0).
    XMax,
    /// Minimum y, outward normal (0,-1,0).
    YMin,
    /// Maximum y, outward normal (0,1,0).
    YMax,
    /// Minimum z, outward normal (0,0,-1).
    ZMin,
    /// Maximum z, outward normal (0,0,1).
    ZMax,
}
impl BoundarySide {
    /// Stable order used for reductions and boundary storage.
    pub const ALL: [Self; 6] = [
        Self::XMin,
        Self::XMax,
        Self::YMin,
        Self::YMax,
        Self::ZMin,
        Self::ZMax,
    ];
    /// Normal Cartesian axis: x=0, y=1, z=2.
    pub const fn axis(self) -> usize {
        match self {
            Self::XMin | Self::XMax => 0,
            Self::YMin | Self::YMax => 1,
            Self::ZMin | Self::ZMax => 2,
        }
    }
    /// Sign of the outward normal along that axis.
    pub const fn sign(self) -> f64 {
        match self {
            Self::XMin | Self::YMin | Self::ZMin => -1.0,
            Self::XMax | Self::YMax | Self::ZMax => 1.0,
        }
    }
}

/// Independent boundary pair for each Cartesian direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Boundaries3d {
    /// Conditions normal to x.
    pub x: AxisBoundaries,
    /// Conditions normal to y.
    pub y: AxisBoundaries,
    /// Conditions normal to z.
    pub z: AxisBoundaries,
}
impl Boundaries3d {
    /// Conditions normal to a physical side.
    pub const fn side(self, side: BoundarySide) -> AxisBoundaries {
        match side.axis() {
            0 => self.x,
            1 => self.y,
            _ => self.z,
        }
    }
    /// Conditions for x=0, y=1 or z=2. Panics for any other index.
    pub fn axis(self, axis: usize) -> AxisBoundaries {
        [self.x, self.y, self.z][axis]
    }
}

/// Uniform physical cells. State order is `(k*ny+j)*nx+i` (x fastest).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CartesianGrid {
    shape: [usize; 3],
    spacing: [f64; 3],
}
impl CartesianGrid {
    /// Admit finite positive geometry, with a bounded, preallocated cell budget.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Result<Self, Error> {
        let n = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(nz))
            .ok_or(Error::Invalid("volumetric Cartesian cell count overflow"))?;
        if nx == 0
            || ny == 0
            || nz == 0
            || n > MAX_EULER_CELLS
            || ![
                dx,
                dy,
                dz,
                dx * dy,
                dx * dz,
                dy * dz,
                dx * dy * dz,
                nx as f64 * dx,
                ny as f64 * dy,
                nz as f64 * dz,
                n as f64 * dx * dy * dz,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(Error::Invalid(
                "volumetric Cartesian geometry or cell budget",
            ));
        }
        Ok(Self {
            shape: [nx, ny, nz],
            spacing: [dx, dy, dz],
        })
    }
    /// Cell counts in x/y/z.
    pub const fn shape(&self) -> [usize; 3] {
        self.shape
    }
    /// Cell widths [m].
    pub const fn spacing(&self) -> [f64; 3] {
        self.spacing
    }
    /// Physical cell volume [m³].
    pub fn cell_volume(&self) -> f64 {
        self.spacing[0] * self.spacing[1] * self.spacing[2]
    }
    /// Physical area of one boundary face [m²].
    pub fn face_area(&self, side: BoundarySide) -> f64 {
        let [a, b] = tangents(side.axis());
        self.spacing[a] * self.spacing[b]
    }
    /// Number of faces on this side.
    pub fn face_count(&self, side: BoundarySide) -> usize {
        self.line_count(side.axis())
    }
    /// Cell centre [m] relative to a caller-selected origin.
    pub fn cell_center(&self, index: usize) -> Option<[f64; 3]> {
        (index < self.len()).then(|| {
            let c = self.coordinates(index);
            std::array::from_fn(|a| c[a] as f64 * self.spacing[a] + 0.5 * self.spacing[a])
        })
    }
    /// Centre of a physical side's face [m], in the documented face order.
    pub fn face_center(&self, side: BoundarySide, face: usize) -> Option<[f64; 3]> {
        (face < self.face_count(side)).then(|| {
            let axis = side.axis();
            let [a, b] = tangents(axis);
            let mut p = [0.0; 3];
            p[axis] = if side.sign() < 0.0 {
                0.0
            } else {
                self.shape[axis] as f64 * self.spacing[axis]
            };
            p[a] = (face % self.shape[a]) as f64 * self.spacing[a] + 0.5 * self.spacing[a];
            p[b] = (face / self.shape[a]) as f64 * self.spacing[b] + 0.5 * self.spacing[b];
            p
        })
    }
    fn len(&self) -> usize {
        self.shape[0] * self.shape[1] * self.shape[2]
    }
    fn index(&self, c: [usize; 3]) -> usize {
        (c[2] * self.shape[1] + c[1]) * self.shape[0] + c[0]
    }
    fn coordinates(&self, index: usize) -> [usize; 3] {
        [
            index % self.shape[0],
            (index / self.shape[0]) % self.shape[1],
            index / (self.shape[0] * self.shape[1]),
        ]
    }
    fn line_count(&self, axis: usize) -> usize {
        let [a, b] = tangents(axis);
        self.shape[a] * self.shape[b]
    }
    fn face_len(&self, axis: usize) -> usize {
        (self.shape[axis] + 1) * self.line_count(axis)
    }
    fn face_strides(&self) -> [usize; 3] {
        [1, self.shape[0], self.shape[0] * self.shape[1]]
    }
    fn cell_face(&self, axis: usize, c: [usize; 3]) -> usize {
        let mut n = self.shape;
        n[axis] += 1;
        (c[2] * n[1] + c[1]) * n[0] + c[0]
    }
    fn boundary_face(&self, axis: usize, position: usize, line: usize) -> usize {
        let [a, b] = tangents(axis);
        let mut c = [0; 3];
        c[axis] = position;
        c[a] = line % self.shape[a];
        c[b] = line / self.shape[a];
        self.cell_face(axis, c)
    }
    fn boundary_len(&self) -> usize {
        2 * (self.line_count(0) + self.line_count(1) + self.line_count(2))
    }
    fn boundary_range(&self, side: BoundarySide) -> std::ops::Range<usize> {
        let axis = side.axis();
        let count = self.line_count(axis);
        let before: usize = (0..axis).map(|a| 2 * self.line_count(a)).sum();
        let start = before + if side.sign() < 0.0 { 0 } else { count };
        start..start + count
    }
}

fn tangents(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

/// Extensive mass, three-axis momentum and referenced total energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Totals {
    /// Mass [kg].
    pub mass: f64,
    /// Momentum [kg m/s].
    pub momentum: [f64; 3],
    /// Referenced energy [J].
    pub energy: f64,
}
impl Totals {
    fn array(self) -> [f64; 5] {
        [
            self.mass,
            self.momentum[0],
            self.momentum[1],
            self.momentum[2],
            self.energy,
        ]
    }
    fn from_array(a: [f64; 5]) -> Self {
        Self {
            mass: a[0],
            momentum: [a[1], a[2], a[3]],
            energy: a[4],
        }
    }
    /// Recover original energy including its fixed source reference [J].
    pub fn absolute_energy(self, gas: &FrozenGas) -> f64 {
        self.energy + gas.energy_reference() * self.mass
    }
}

/// Outward time-integrated face transfer. Wall momentum is support impulse.
pub type BoundaryImpulse = Totals;

/// A complete accepted RK step; failed attempts publish no field or load state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Accepted duration [s].
    pub elapsed: Time,
    /// Largest sum of face CFL contributions from x, y and z over both stages.
    pub maximum_cfl: f64,
    /// Rejected trials before acceptance.
    pub halvings: usize,
    /// Accepted-stage faces that used HLL fallback.
    pub hll_fallback_faces: usize,
    /// Accepted-stage cells with source-domain-limited reconstruction.
    pub limited_cells: usize,
    /// An admissibility refusal required a whole-trial first-order retry.
    pub first_order_retry: bool,
    /// Extensive state before this step.
    pub before: Totals,
    /// Extensive state after this step.
    pub after: Totals,
    /// Net inward transfer integrated with the same two RK stages.
    pub boundary_change: Totals,
    /// `after - before - boundary_change` in coherent SI units.
    pub conservation_residual: Totals,
}

struct Workspace {
    primitive: Vec<Primitive>,
    slopes: Vec<[[f64; 5]; 3]>,
    flux: [Vec<SpatialFaceFlux>; 3],
    stage: Vec<State>,
    candidate: Vec<State>,
    boundary_trial: Vec<BoundaryImpulse>,
}

/// Conservative volumetric shock field over the original immutable gas source.
///
/// All field and face buffers are allocated at construction. Retries operate on
/// private stage storage; cancellation preserves cells, time and accepted loads.
pub struct Euler3d {
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    states: Vec<State>,
    time: f64,
    last_step: Option<StepReport>,
    boundary_impulses: Vec<BoundaryImpulse>,
    work: Workspace,
}
impl Euler3d {
    /// Admit original gas thermodynamics and cell averages on a fixed box.
    pub fn new(
        gas: FrozenGas,
        grid: CartesianGrid,
        states: Vec<State>,
        boundaries: Boundaries3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if states.len() != grid.len() {
            return Err(Error::Invalid("volumetric Cartesian state count"));
        }
        let n = states.len();
        let mut primitive = allocated(n, Primitive::default())?;
        for (p, q) in primitive.iter_mut().zip(&states) {
            *p = q.primitive(&gas, &mut cancelled)?;
        }
        if !checked_totals(&states, grid.cell_volume(), &mut cancelled)?
            .array()
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid("unrepresentable volumetric extensive state"));
        }
        let result = Self {
            gas,
            grid,
            boundaries,
            states,
            time: 0.0,
            last_step: None,
            boundary_impulses: allocated(grid.boundary_len(), Totals::default())?,
            work: Workspace {
                primitive,
                slopes: allocated(n, [[0.0; 5]; 3])?,
                flux: [
                    allocated(grid.face_len(0), SpatialFaceFlux::default())?,
                    allocated(grid.face_len(1), SpatialFaceFlux::default())?,
                    allocated(grid.face_len(2), SpatialFaceFlux::default())?,
                ],
                stage: allocated(n, State::default())?,
                candidate: allocated(n, State::default())?,
                boundary_trial: allocated(grid.boundary_len(), Totals::default())?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
    /// Original source-bearing EOS owner.
    pub const fn gas(&self) -> &FrozenGas {
        &self.gas
    }
    /// Fixed physical grid.
    pub const fn grid(&self) -> &CartesianGrid {
        &self.grid
    }
    /// Boundary treatments.
    pub const fn boundaries(&self) -> Boundaries3d {
        self.boundaries
    }
    /// Current cell averages, x fastest.
    pub fn states(&self) -> &[State] {
        &self.states
    }
    /// Accepted time [s].
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Last accepted report; absent before the first step.
    pub const fn last_step(&self) -> Option<&StepReport> {
        self.last_step.as_ref()
    }
    /// Original per-face outward transfer integrated over the accepted step.
    pub fn boundary_impulses(&self, side: BoundarySide) -> &[BoundaryImpulse] {
        &self.boundary_impulses[self.grid.boundary_range(side)]
    }
    /// Compensated extensive totals in deterministic cell order.
    pub fn totals(&self) -> Totals {
        checked_totals(&self.states, self.grid.cell_volume(), &mut || false)
            .expect("infallible cancellation callback")
    }
    /// Advance at most this duration, respecting the combined three-axis CFL.
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
            return Err(Error::Invalid("volumetric duration, CFL or retry controls"));
        }
        let mut rate: f64 = 0.0;
        for (p, q) in self.work.primitive.iter_mut().zip(&self.states) {
            *p = q.primitive(&self.gas, &mut cancelled)?;
            let local: f64 = (0..3)
                .map(|a| (p.velocity[a].abs() + p.sound_speed) / self.grid.spacing[a])
                .sum();
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(controls.cfl / rate);
        let before = checked_totals(&self.states, self.grid.cell_volume(), &mut cancelled)?;
        let mut reconstruction = controls.reconstruction;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || self.time + dt <= self.time || !(self.time + dt).is_finite() {
                return Err(Error::Invalid("unrepresentable volumetric time increment"));
            }
            let trial = evolve::attempt(
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
                            "unrepresentable volumetric extensive balance",
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
        unreachable!("bounded attempts return at the last retry")
    }
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

#[cfg(test)]
mod tests;
