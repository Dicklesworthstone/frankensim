//! Native D3Q19/D3Q27 velocity and pressure reservoirs.
//!
//! The original collision and pull stencil fill the original transactional
//! buffers. Regularized open-face populations are completed before physical
//! admission and publication. Donor moments include the half force; target
//! populations use the native equilibrium and complete half-source correction.
//! This is the same second-moment boundary model as the specialized D3Q19
//! owner, extended to both native quadratures without selecting a collision law.
//!
//! Open faces occupy one Cartesian axis. At wall/open seams, only open-owned
//! incoming rows are replaced: prescribed full moments apply to face interiors,
//! not to the intersection of incompatible wall and reservoir conditions.
//! Pressure and velocity profiles are held on their native face planes;
//! geometry and the uniform body force remain fixed. Tangential wall motion uses the
//! native moving reflection and a separate wall impulse/work account.
//! Characteristic outlets remain separate.

use super::{
    EquilibriumLattice3, FirstError, Grid3, GridError3, Lattice3, StepData, coordinates, initial,
    macro_cell, raw_moments_cell, run_pass,
};
use crate::d3q19::{
    BoundaryForcingError3, BoundarySpec3, Face3, FaceBoundary3, OpenBoundaryMomentumTransport3,
};
use core::{marker::PhantomData, ops::ControlFlow};
use fs_exec::{CancelGate, Cancelled, Cx, KernelRunner, TileKernel, TilePlan};
use fs_math::dd::Dd;
use std::sync::Mutex;

/// Spatial force epochs on the same native reservoir and population owner.
pub mod forcing;
mod moving_walls;
mod pressure;
mod rheology;
/// Explicit local turbulent momentum and heat on original reservoir geometry.
pub mod sgs;
/// Conservative D3Q7 heat transport on the same native open-grid transaction.
pub mod thermal;
mod velocity;

use pressure::PressureProfile3;
use velocity::VelocityProfile3;

/// Original native grid with fixed-topology velocity/pressure face ownership.
/// Populations, physical time and measured reservoir exchange publish together.
pub struct OpenBoundaryGrid3<L: EquilibriumLattice3<Q>, const Q: usize> {
    grid: Grid3<L, Q>,
    plan: FacePlan,
    pressure_profiles: [Option<PressureProfile3>; 6],
    velocity_profiles: [Option<VelocityProfile3>; 6],
    wall_exchange: Option<super::moving_walls::MovingWallMomentum3>,
    pending_wall_exchange: Option<super::moving_walls::MovingWallMomentum3>,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Bind reservoirs without modifying populations or advancing time.
    /// Periodicity must match the original grid. The open axis needs at least
    /// three cells; every exposed fluid boundary cell needs a fluid donor one
    /// cell inward. Wall motion must be tangential, compatible at wall junctions
    /// and locally mass conserving on the actual native reflected links.
    pub fn try_new(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        boundaries: BoundarySpec3,
    ) -> Result<Self, GridError3<L::Error>> {
        Self::try_new_with_poll(grid, boundaries, &mut || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    pub(super) fn try_new_with_poll(
        grid: Grid3<L, Q>,
        boundaries: BoundarySpec3,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Self, GridError3<L::Error>> {
        let plan = FacePlan::admit(&grid, boundaries, poll)?;
        for (index, cell) in plan.cells().enumerate() {
            if index % 64 == 0 {
                poll()?;
            }
            if !grid.solid[cell]
                && let Some((_, donor)) = plan.locate(cell)
                && grid.solid[donor]
            {
                return Err(GridError3::MaterialField {
                    cell,
                    reason: "open face requires a fluid interior donor",
                });
            }
        }
        poll()?;
        Ok(Self {
            grid,
            plan,
            pressure_profiles: core::array::from_fn(|_| None),
            velocity_profiles: core::array::from_fn(|_| None),
            wall_exchange: None,
            pending_wall_exchange: None,
        })
    }

    /// Grid extents in x/y/z order; cell storage is x-fastest.
    pub fn dimensions(&self) -> [usize; 3] {
        self.grid.dimensions()
    }
    /// Accepted populations in their uniform half-source frame.
    pub fn populations(&self) -> &[[f64; Q]] {
        self.grid.populations()
    }
    /// Accepted density and physical, force-corrected velocity.
    pub fn macroscopic(&self, cell: usize) -> Result<(f64, [f64; 3]), GridError3<L::Error>> {
        self.grid.macroscopic(cell)
    }
    /// Fully accepted native steps.
    pub fn completed_steps(&self) -> usize {
        self.grid.completed_steps()
    }
    /// Held force density in lattice units, not acceleration.
    pub fn force_density(&self) -> [f64; 3] {
        self.grid.force
    }
    /// Current face targets with the original registry ordering.
    pub fn boundaries(&self) -> BoundarySpec3 {
        self.plan.boundaries
    }

    /// Atomically replace values on existing velocity/pressure/wall faces.
    /// Changing kinds or periodicity requires a different owner;
    /// this update performs no streaming, impulse or hidden population reset.
    /// Spatial pressure/velocity profiles remain authoritative: clear a profile
    /// before changing its held uniform fallback, or update the profile itself.
    /// Unrelated wall/reservoir updates preserve every attached profile.
    pub fn rebind_open_boundaries(
        &mut self,
        boundaries: BoundarySpec3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        let candidate = self.prepare_open_boundaries(boundaries, &mut || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.plan = candidate;
        Ok(())
    }

    fn prepare_open_boundaries(
        &self,
        boundaries: BoundarySpec3,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<FacePlan, GridError3<L::Error>> {
        poll()?;
        for face in Face3::ALL {
            match (self.plan.boundaries.face(face), boundaries.face(face)) {
                (FaceBoundary3::Velocity { .. }, FaceBoundary3::Velocity { .. })
                | (FaceBoundary3::Pressure { .. }, FaceBoundary3::Pressure { .. })
                | (FaceBoundary3::Wall { .. }, FaceBoundary3::Wall { .. }) => {}
                (old, new) if old == new => {}
                _ => {
                    return Err(GridError3::Input(
                        "open target update changes boundary kind",
                    ));
                }
            }
            if self.plan.boundaries.face(face) != boundaries.face(face)
                && self.pressure_profiles[face as usize].is_some()
            {
                return Err(GridError3::Input(
                    "clear or update the pressure profile before changing its uniform target",
                ));
            }
            if self.plan.boundaries.face(face) != boundaries.face(face)
                && self.velocity_profiles[face as usize].is_some()
            {
                return Err(GridError3::Input(
                    "clear or update the velocity profile before changing its uniform target",
                ));
            }
        }
        let candidate = FacePlan::admit(&self.grid, boundaries, poll)?;
        poll()?;
        Ok(candidate)
    }

    /// Advance the original collision/stream operator with held reservoirs.
    /// `collisions` optionally supplies the native per-cell material/SGS laws.
    /// `origin` is the lattice-coordinate origin for the measured angular
    /// exchange. No boundary pass or receipt failure publishes a partial step.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        self.step_serial_with_plan(cx, self.plan, collisions, origin)
    }

    /// Advance with new reservoir targets in the same native transaction.
    ///
    /// Boundary kinds, geometry and spatial-profile ownership retain the
    /// checks of [`Self::rebind_open_boundaries`]. The candidate targets drive
    /// this step's reconstruction and measured transport, but become held
    /// controls only when populations and the completed-step count publish.
    /// Collision, streaming, measurement or cancellation failure preserves
    /// the previous accepted targets as well as the previous fluid state.
    pub fn step_serial_with_open_boundaries(
        &mut self,
        cx: &Cx<'_>,
        boundaries: BoundarySpec3,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let plan = self.prepare_open_boundaries(boundaries, &mut || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.step_serial_with_plan(cx, plan, collisions, origin)
    }

    fn step_serial_with_plan(
        &mut self,
        cx: &Cx<'_>,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let (receipt, step) = self.stage_serial_with_plan(cx, plan, collisions, origin)?;
        core::mem::swap(&mut self.grid.f, &mut self.grid.next);
        self.grid.completed_steps = step;
        self.plan = plan;
        self.publish_wall_exchange();
        Ok(receipt)
    }

    // The thermal owner admits its final source kick and independent heat
    // receipt before either native population field becomes visible.
    fn stage_serial_with_plan(
        &mut self,
        cx: &Cx<'_>,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<(OpenBoundaryMomentumTransport3, usize), GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        self.stage_serial_with_forces(plan, collisions, None, None, origin, &mut poll)
    }

    // One original collision/stream operator serves both the uniform owner
    // and the source-history wrapper. End forces interpret the streamed donor
    // and target populations; they are not applied as a second collision.
    #[allow(clippy::too_many_arguments)]
    fn stage_serial_with_forces(
        &mut self,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        forces_at_start: Option<&[[f64; 3]]>,
        forces_at_end: Option<&[[f64; 3]]>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(OpenBoundaryMomentumTransport3, usize), GridError3<L::Error>> {
        self.stage_serial_with_source(
            plan,
            collisions,
            forces_at_start,
            forces_at_end,
            None,
            origin,
            poll,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn stage_serial_with_source(
        &mut self,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        forces_at_start: Option<&[[f64; 3]]>,
        forces_at_end: Option<&[[f64; 3]]>,
        acceleration: Option<(&[[f64; 3]], &mut [[f64; 3]])>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(OpenBoundaryMomentumTransport3, usize), GridError3<L::Error>> {
        let step = self.admit_step(collisions, origin, poll)?;
        for forces in [forces_at_start, forces_at_end].into_iter().flatten() {
            super::sources::admit_forces(&self.grid, forces, collisions, poll)?;
        }
        let (acceleration, mut resolved) =
            self.admit_acceleration(acceleration, collisions, poll)?;
        let grid = &mut self.grid;
        let data = StepData::<L, Q> {
            dims: grid.dims,
            periodic: grid.periodic,
            solid: &grid.solid,
            collision: grid.collision,
            collisions,
            force: grid.force,
            forces: forces_at_start,
            budget: grid.budget,
            marker: PhantomData,
        };
        for (tile, chunk) in grid.post.chunks_mut(64).enumerate() {
            poll()?;
            data.process(&grid.f, chunk, tile * 64, false)?;
        }
        let view = OpenStep {
            data: StepData {
                forces: forces_at_end,
                ..data
            },
            plan,
            pressure_profiles: &self.pressure_profiles,
            velocity_profiles: &self.velocity_profiles,
            acceleration,
        };
        for (tile, chunk) in grid.next.chunks_mut(64).enumerate() {
            poll()?;
            let field = resolved
                .as_deref_mut()
                .map(|field| &mut field[tile * 64..tile * 64 + chunk.len()]);
            view.process(&grid.post, chunk, tile * 64, field)?;
        }
        let receipt = view.measure(&grid.post, &grid.next, origin, poll)?;
        let walls = view.measure_walls(&grid.post, &grid.next, origin, poll)?;
        poll()?;
        self.pending_wall_exchange = walls;
        Ok((receipt, step))
    }

    /// Executor counterpart using the same fixed 64-cell collision and stream
    /// tiles. The receipt is reduced in canonical cell/direction order after
    /// the candidate pass, so worker count does not change its arithmetic.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        self.step_pooled_with_plan(pool, gate, self.plan, collisions, origin)
    }

    /// Executor counterpart of [`Self::step_serial_with_open_boundaries`].
    /// The same candidate targets feed the original tiled operator and
    /// canonical receipt reduction; all accepted state publishes together.
    pub fn step_pooled_with_open_boundaries<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        boundaries: BoundarySpec3,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let plan = self.prepare_open_boundaries(boundaries, &mut || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.step_pooled_with_plan(pool, gate, plan, collisions, origin)
    }

    fn step_pooled_with_plan<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let (receipt, step) =
            self.stage_pooled_with_source(pool, gate, plan, collisions, None, None, None, origin)?;
        core::mem::swap(&mut self.grid.f, &mut self.grid.next);
        self.grid.completed_steps = step;
        self.plan = plan;
        self.publish_wall_exchange();
        Ok(receipt)
    }

    #[allow(clippy::too_many_arguments)]
    fn stage_pooled_with_source<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        plan: FacePlan,
        collisions: Option<&[L::Collision]>,
        forces_at_start: Option<&[[f64; 3]]>,
        forces_at_end: Option<&[[f64; 3]]>,
        acceleration: Option<(&[[f64; 3]], &mut [[f64; 3]])>,
        origin: [f64; 3],
    ) -> Result<(OpenBoundaryMomentumTransport3, usize), GridError3<L::Error>> {
        let mut poll = || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        };
        let step = self.admit_step(collisions, origin, &mut poll)?;
        for forces in [forces_at_start, forces_at_end].into_iter().flatten() {
            super::sources::admit_forces(&self.grid, forces, collisions, &mut poll)?;
        }
        let (acceleration, resolved) =
            self.admit_acceleration(acceleration, collisions, &mut poll)?;
        let grid = &mut self.grid;
        let data = StepData::<L, Q> {
            dims: grid.dims,
            periodic: grid.periodic,
            solid: &grid.solid,
            collision: grid.collision,
            collisions,
            force: grid.force,
            forces: forces_at_start,
            budget: grid.budget,
            marker: PhantomData,
        };
        run_pass(&data, &grid.f, &mut grid.post, false, pool, gate)?;
        let view = OpenStep {
            data: StepData {
                forces: forces_at_end,
                ..data
            },
            plan,
            pressure_profiles: &self.pressure_profiles,
            velocity_profiles: &self.velocity_profiles,
            acceleration,
        };
        {
            let chunks = match resolved {
                Some(field) => grid
                    .next
                    .chunks_mut(64)
                    .zip(field.chunks_mut(64))
                    .map(|(pop, forces)| Mutex::new((pop, Some(forces))))
                    .collect(),
                None => grid
                    .next
                    .chunks_mut(64)
                    .map(|pop| Mutex::new((pop, None)))
                    .collect(),
            };
            let kernel = OpenStream {
                view: &view,
                src: &grid.post,
                chunks,
            };
            match pool.run_with_gate(&kernel, gate).0 {
                Ok(FirstError(None)) => {}
                Ok(FirstError(Some(error))) => return Err(error),
                Err(_) if gate.is_requested() => return Err(GridError3::Cancelled),
                Err(error) => return Err(GridError3::Pool(format!("{error:?}"))),
            }
        }
        let receipt = view.measure(&grid.post, &grid.next, origin, &mut poll)?;
        let walls = view.measure_walls(&grid.post, &grid.next, origin, &mut poll)?;
        poll()?;
        self.pending_wall_exchange = walls;
        Ok((receipt, step))
    }

    fn admit_acceleration<'a>(
        &self,
        field: Option<(&'a [[f64; 3]], &'a mut [[f64; 3]])>,
        collisions: Option<&[L::Collision]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(Option<&'a [[f64; 3]]>, Option<&'a mut [[f64; 3]]>), GridError3<L::Error>> {
        match field {
            None => Ok((None, None)),
            Some((acceleration, resolved)) => {
                super::sources::admit_forces(&self.grid, acceleration, collisions, poll)?;
                if resolved.len() != self.grid.f.len() {
                    return Err(GridError3::Input("endpoint force workspace length"));
                }
                Ok((Some(acceleration), Some(resolved)))
            }
        }
    }

    fn admit_step(
        &self,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<usize, GridError3<L::Error>> {
        poll()?;
        if origin.iter().any(|v| !v.is_finite()) {
            return Err(GridError3::Input("nonfinite open-reservoir moment origin"));
        }
        if let Some(field) = collisions {
            self.grid.admit_collisions(field, &mut *poll)?;
        }
        self.grid
            .completed_steps
            .checked_add(1)
            .ok_or(GridError3::Input("step count overflow"))
    }
}

#[derive(Clone, Copy)]
struct FacePlan {
    boundaries: BoundarySpec3,
    dims: [usize; 3],
    axis: Option<usize>,
    stride: usize,
    moving: bool,
}

fn is_open(boundary: FaceBoundary3) -> bool {
    matches!(
        boundary,
        FaceBoundary3::Velocity { .. } | FaceBoundary3::Pressure { .. }
    )
}

impl FacePlan {
    fn admit<L: EquilibriumLattice3<Q>, const Q: usize>(
        grid: &Grid3<L, Q>,
        boundaries: BoundarySpec3,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Self, GridError3<L::Error>> {
        boundaries.try_validate().map_err(|error| match error {
            BoundaryForcingError3::Invalid(reason) => GridError3::Input(reason),
            _ => GridError3::Input("invalid native open boundaries"),
        })?;
        let mut axis = None;
        let mut moving = false;
        for face in Face3::ALL {
            let a = face as usize / 2;
            let boundary = boundaries.face(face);
            if matches!(boundary, FaceBoundary3::Periodic) != grid.periodic[a] {
                return Err(GridError3::Input(
                    "boundary periodicity differs from grid topology",
                ));
            }
            match boundary {
                FaceBoundary3::Wall { velocity } => {
                    if velocity[a] != 0.0 {
                        return Err(GridError3::Input(
                            "fixed exterior wall velocity must be tangential",
                        ));
                    }
                    initial::admit_velocity::<L, Q>(
                        grid.collision,
                        velocity,
                        0,
                        grid.budget.maximum_mach,
                    )?;
                    moving |= velocity != [0.0; 3];
                }
                FaceBoundary3::Velocity { velocity } => {
                    initial::admit_velocity::<L, Q>(
                        grid.collision,
                        velocity,
                        0,
                        grid.budget.maximum_mach,
                    )?;
                }
                _ => {}
            }
            if is_open(boundary) {
                if grid.dims[a] < 3 {
                    return Err(GridError3::Input("open axis requires three or more cells"));
                }
                axis = Some(a);
            }
        }
        let plan = Self {
            boundaries,
            dims: grid.dims,
            axis,
            stride: axis.map_or(1, |a| grid.dims[..a].iter().product()),
            moving,
        };
        plan.admit_wall_motion(grid, poll)?;
        Ok(plan)
    }

    fn locate(&self, cell: usize) -> Option<(Face3, usize)> {
        let axis = self.axis?;
        let coordinate = cell / self.stride % self.dims[axis];
        let (index, donor) = if coordinate == 0 {
            (2 * axis, cell + self.stride)
        } else if coordinate + 1 == self.dims[axis] {
            (2 * axis + 1, cell - self.stride)
        } else {
            return None;
        };
        let face = Face3::ALL[index];
        is_open(self.boundaries.face(face)).then_some((face, donor))
    }

    // Remove the fixed normal coordinate, retaining x-fast plane ordering.
    fn plane_index(&self, cell: usize) -> usize {
        let axis = self.axis.expect("only open faces have profile coordinates");
        cell / (self.stride * self.dims[axis]) * self.stride + cell % self.stride
    }

    // Enumerate only exposed planes, still in the original x-fast cell order.
    // Admission and receipts scale with surface area and allocate no cell map.
    fn cells(self) -> impl Iterator<Item = usize> {
        self.axis.into_iter().flat_map(move |axis| {
            let dimension = self.dims[axis];
            let block = self.stride * dimension;
            let blocks = self.dims.iter().product::<usize>() / block;
            (0..blocks).flat_map(move |outer| {
                [(0, 2 * axis), (dimension - 1, 2 * axis + 1)]
                    .into_iter()
                    .filter(move |&(_, face)| is_open(self.boundaries.face(Face3::ALL[face])))
                    .flat_map(move |(coordinate, _)| {
                        let start = outer * block + coordinate * self.stride;
                        start..start + self.stride
                    })
            })
        })
    }

    // Match the native wall/open ownership rule: a link crossing any wall is
    // wall-owned even if it also crosses the open plane. Interior solid links
    // are walls; periodic coordinates wrap before their solid query.
    fn masks<L: Lattice3<Q>, const Q: usize>(
        &self,
        data: &StepData<'_, L, Q>,
        cell: usize,
    ) -> (u32, u32) {
        let p = coordinates(cell, data.dims);
        let mut walls = 0;
        let mut open = 0;
        for q in 1..Q {
            let mut from = [0; 3];
            let mut wall = false;
            let mut reservoir = false;
            for a in 0..3 {
                let v = p[a] as isize - L::E[q][a] as isize;
                if v < 0 || v >= data.dims[a] as isize {
                    let face = Face3::ALL[2 * a + usize::from(v >= 0)];
                    match self.boundaries.face(face) {
                        FaceBoundary3::Periodic => {
                            from[a] = v.rem_euclid(data.dims[a] as isize) as usize;
                        }
                        FaceBoundary3::Wall { .. } => wall = true,
                        _ => reservoir = true,
                    }
                } else {
                    from[a] = v as usize;
                }
            }
            if !wall && !reservoir {
                let source = (from[2] * data.dims[1] + from[1]) * data.dims[0] + from[0];
                wall = data.solid[source];
            }
            if wall {
                walls |= 1 << q;
            } else if reservoir {
                open |= 1 << q;
            }
        }
        (walls, open)
    }
}

struct OpenStep<'a, L: EquilibriumLattice3<Q>, const Q: usize> {
    data: StepData<'a, L, Q>,
    plan: FacePlan,
    pressure_profiles: &'a [Option<PressureProfile3>; 6],
    velocity_profiles: &'a [Option<VelocityProfile3>; 6],
    acceleration: Option<&'a [[f64; 3]]>,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenStep<'_, L, Q> {
    fn process(
        &self,
        src: &[[f64; Q]],
        dst: &mut [[f64; Q]],
        base: usize,
        mut resolved: Option<&mut [[f64; 3]]>,
    ) -> Result<(), GridError3<L::Error>> {
        self.data
            .process_stream_with_update::<false>(src, dst, base, |cell, output, _| {
                self.reflect_walls(src, output, cell)?;
                if let Some((face, donor)) = self.plan.locate(cell) {
                    let (walls, open) = self.plan.masks(&self.data, cell);
                    if open != 0 {
                        self.reconstruct(src, cell, donor, face, walls, open, output)?;
                    }
                }
                let (_, velocity, force) = self.endpoint_state(output, cell)?;
                let collision = self
                    .data
                    .collisions
                    .map_or(self.data.collision, |field| field[cell]);
                initial::admit_velocity::<L, Q>(
                    collision,
                    velocity,
                    cell,
                    self.data.budget.maximum_mach,
                )?;
                if let Some(field) = resolved.as_deref_mut() {
                    field[cell - base] = force;
                }
                Ok(())
            })
    }

    // Donor and final observations use their OWN actual density. Computing
    // moments once also supplies the resolved endpoint field inside its worker.
    fn endpoint_state(
        &self,
        populations: &[f64; Q],
        cell: usize,
    ) -> Result<(f64, [f64; 3], [f64; 3]), GridError3<L::Error>> {
        if let Some(field) = self.acceleration {
            let (rho, momentum) = raw_moments_cell::<L, Q>(populations, cell)?;
            let force = field[cell].map(|a| rho * a);
            let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
            Ok((rho, velocity, force))
        } else {
            let force = self
                .data
                .forces
                .map_or(self.data.force, |field| field[cell]);
            let (rho, velocity) =
                macro_cell::<L, Q>(populations, force, cell, self.data.budget.maximum_mach)?;
            Ok((rho, velocity, force))
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn reconstruct(
        &self,
        src: &[[f64; Q]],
        cell: usize,
        donor: usize,
        face: Face3,
        walls: u32,
        open: u32,
        output: &mut [f64; Q],
    ) -> Result<(), GridError3<L::Error>> {
        let mut source = [0.0; Q];
        self.data.pull_cell(src, &mut source, donor);
        self.reflect_walls(src, &mut source, donor)?;
        let (rho_source, velocity_source, force) = self.endpoint_state(&source, donor)?;
        let donor_law = self
            .data
            .collisions
            .map_or(self.data.collision, |field| field[donor]);
        initial::admit_velocity::<L, Q>(
            donor_law,
            velocity_source,
            donor,
            self.data.budget.maximum_mach,
        )?;
        let equilibrium = initial::equilibrium::<L, Q>(donor_law, rho_source, velocity_source);
        let mut stress = [[0.0; 3]; 3];
        for q in 0..Q {
            let difference = source[q] - equilibrium[q];
            for a in 0..3 {
                for b in 0..3 {
                    stress[a][b] += f64::from(L::E[q][a] * L::E[q][b]) * difference;
                }
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                stress[a][b] +=
                    0.5 * (velocity_source[a] * force[b] + force[a] * velocity_source[b]);
            }
        }
        let (rho, velocity) = match self.plan.boundaries.face(face) {
            FaceBoundary3::Velocity { velocity } => (
                rho_source,
                self.velocity_profiles[face as usize]
                    .as_ref()
                    .map_or(velocity, |profile| {
                        profile.values[self.plan.plane_index(cell)]
                    }),
            ),
            FaceBoundary3::Pressure { density } => (
                self.pressure_profiles[face as usize]
                    .as_ref()
                    .map_or(density, |profile| {
                        profile.values[self.plan.plane_index(cell)]
                    }),
                velocity_source,
            ),
            _ => unreachable!("only open faces have a reservoir reconstruction"),
        };
        let law = self
            .data
            .collisions
            .map_or(self.data.collision, |field| field[cell]);
        initial::admit_velocity::<L, Q>(law, velocity, cell, self.data.budget.maximum_mach)?;
        let mut reconstructed = initial::equilibrium::<L, Q>(law, rho, velocity);
        for q in 0..Q {
            let mut contraction = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let isotropic = if a == b { L::CS2 } else { 0.0 };
                    contraction += (f64::from(L::E[q][a] * L::E[q][b]) - isotropic) * stress[a][b];
                }
            }
            reconstructed[q] += L::W[q] * contraction / (2.0 * L::CS2 * L::CS2);
        }
        if let Some(field) = self.acceleration {
            return self.finish_accelerated_reconstruction(
                cell,
                law,
                velocity,
                field[cell],
                walls,
                open,
                reconstructed,
                output,
            );
        }
        let target_force = self
            .data
            .forces
            .map_or(self.data.force, |field| field[cell]);
        initial::attach_force::<L, Q>(&mut reconstructed, law, velocity, target_force);
        for q in 0..Q {
            if walls == 0 || open & (1 << q) != 0 {
                output[q] = reconstructed[q];
            }
        }
        Ok(())
    }

    // At a wall/open seam only reservoir-owned rows receive -S(u,rho*a)/2.
    // Their retained mass is affine in rho, so close it directly instead of
    // freezing the nominal pressure/donor density or iterating a source guess.
    #[allow(clippy::too_many_arguments)]
    fn finish_accelerated_reconstruction(
        &self,
        cell: usize,
        law: L::Collision,
        velocity: [f64; 3],
        acceleration: [f64; 3],
        walls: u32,
        open: u32,
        mut reconstructed: [f64; Q],
        output: &mut [f64; Q],
    ) -> Result<(), GridError3<L::Error>> {
        let owned = |q| walls == 0 || open & (1 << q) != 0;
        let mut density = reconstructed.iter().sum::<f64>();
        if walls != 0 {
            let mut response = [0.0; Q];
            initial::attach_force::<L, Q>(&mut response, law, velocity, acceleration);
            let mass: f64 = (0..Q)
                .map(|q| {
                    if owned(q) {
                        reconstructed[q]
                    } else {
                        output[q]
                    }
                })
                .sum();
            let denominator = 1.0
                - (0..Q)
                    .filter(|&q| owned(q))
                    .map(|q| response[q])
                    .sum::<f64>();
            if !mass.is_finite() || mass <= 0.0 || !denominator.is_finite() || denominator <= 0.0 {
                return Err(GridError3::MaterialField {
                    cell,
                    reason: "density-dependent reservoir source has no admissible root",
                });
            }
            density = mass / denominator;
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(GridError3::NonFiniteCell { cell });
        }
        initial::attach_force::<L, Q>(
            &mut reconstructed,
            law,
            velocity,
            acceleration.map(|a| density * a),
        );
        for q in 0..Q {
            if owned(q) {
                output[q] = reconstructed[q];
            }
        }
        let actual = output.iter().sum::<f64>();
        let scale = output.iter().map(|f| f.abs()).sum::<f64>().max(density);
        if !actual.is_finite()
            || !scale.is_finite()
            || (actual - density).abs() > 64.0 * f64::EPSILON * scale
        {
            return Err(GridError3::MaterialField {
                cell,
                reason: "density-dependent reservoir source lost its mass closure",
            });
        }
        Ok(())
    }

    fn measure(
        &self,
        post: &[[f64; Q]],
        next: &[[f64; Q]],
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let mut mass = Dd::ZERO;
        let mut impulse = [Dd::ZERO; 3];
        let mut angular = [Dd::ZERO; 3];
        let mut receipt = OpenBoundaryMomentumTransport3::default();
        for (index, cell) in self.plan.cells().enumerate() {
            if index % 64 == 0 {
                poll()?;
            }
            if self.data.solid[cell] {
                continue;
            }
            let (_, open) = self.plan.masks(&self.data, cell);
            if open == 0 {
                continue;
            }
            let mut before = [0.0; Q];
            self.data.pull_cell(post, &mut before, cell);
            self.reflect_walls(post, &mut before, cell)?;
            let p = coordinates(cell, self.data.dims);
            let arm: [Dd; 3] =
                core::array::from_fn(|a| Dd::from_f64(p[a] as f64 + 0.5 - origin[a]));
            for q in 0..Q {
                let difference = Dd::from_f64(next[cell][q]) - Dd::from_f64(before[q]);
                mass = mass + difference;
                let transfer = if open & (1 << q) != 0 {
                    receipt.open_links += 1;
                    difference + Dd::from_f64(2.0) * Dd::from_f64(post[cell][L::OPP[q]])
                } else {
                    difference
                };
                let momentum: [Dd; 3] =
                    core::array::from_fn(|a| transfer * Dd::from_f64(f64::from(L::E[q][a])));
                for a in 0..3 {
                    let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                    impulse[a] = impulse[a] + momentum[a];
                    angular[a] = angular[a] + arm[b] * momentum[c] - arm[c] * momentum[b];
                }
            }
            receipt.reconstructed_cells += 1;
        }
        receipt.mass_change = mass.to_f64();
        receipt.impulse = impulse.map(Dd::to_f64);
        receipt.angular_impulse = angular.map(Dd::to_f64);
        if [receipt.mass_change]
            .into_iter()
            .chain(receipt.impulse)
            .chain(receipt.angular_impulse)
            .any(|v| !v.is_finite())
        {
            return Err(GridError3::Input(
                "open-reservoir momentum account overflowed",
            ));
        }
        poll()?;
        Ok(receipt)
    }
}

struct OpenStream<'a, 'b, L: EquilibriumLattice3<Q>, const Q: usize> {
    view: &'a OpenStep<'b, L, Q>,
    src: &'a [[f64; Q]],
    chunks: Vec<Mutex<(&'a mut [[f64; Q]], Option<&'a mut [[f64; 3]]>)>>,
}
impl<L: EquilibriumLattice3<Q>, const Q: usize> TileKernel for OpenStream<'_, '_, L, Q> {
    type Out = FirstError<L::Error>;
    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm/grid3-open-stream", self.chunks.len() as u64)
    }
    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let outcome = match self.chunks[tile as usize].lock() {
            Ok(mut chunk) => {
                let (populations, forces) = &mut *chunk;
                self.view.process(
                    self.src,
                    populations,
                    tile as usize * 64,
                    forces.as_deref_mut(),
                )
            }
            Err(_) => Err(GridError3::Poisoned),
        };
        ControlFlow::Continue(FirstError(outcome.err()))
    }
}
