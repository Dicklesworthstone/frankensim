//! Implicit porous momentum on the original shared D3Q19/D3Q27 grid.
//!
//! The constitutive solve and force/work law are shared with BoundaryGrid3.
//! Collision families keep their own forcing implementation: no D3Q19 Guo
//! population transform is imposed on a D3Q27 product or cumulant model.

use super::{ForcedGrid3, ForcedLattice3};
use crate::d3q19::{BoundaryStepError3, WallMomentumExchange3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::lattice3::{
    FirstError, Grid3, GridError3, StepData, coordinates, macro_cell, raw_moments_cell, run_pass,
};
use core::{marker::PhantomData, ops::ControlFlow};
use fs_exec::{CancelGate, Cancelled, Cx, KernelRunner, TileKernel, TilePlan};
use std::sync::Mutex;

pub use crate::d3q19::{DarcyForchheimer3, LinearResistance3, ResistanceMomentumExchange3};

/// A frozen spatial porous medium on native lattice collision and transport.
///
/// Each fluid cell obeys `F = b + D(v_s-u) - beta |u-v_s| (u-v_s)` in lattice
/// units. Resistance, prescribed skeleton velocity and body force stay fixed
/// for this owner's lifetime. The current solved force drives collision; after
/// the native pull, the same law and `rho*u = m_raw + F/2` determine the next
/// force analytically or by the shared monotone scalar speed solve.
///
/// Stationary exterior/voxel walls and periodic axes retain the supplied grid's
/// topology. Every wall link is measured. Periodic angular transport is counted
/// directly in the native pull, including D3Q27 corner wraps. The receipt thus
/// supplies the complete momentum account without inferring any boundary term
/// from a global residual. Populations, force history, receipt and step count
/// publish only after all collision, endpoint, wall and cancellation checks.
///
/// Work is the shared midpoint force-work identity, not a transported fluid
/// energy balance; no heat or skeleton inertia evolves. Trapezoidal relaxation
/// is second order when resolved, but stiff steps can reverse relative velocity.
/// All original collision and Mach checks still apply. Moving/open boundaries
/// and changing resistance laws require their own source/frame coupling.
/// Endpoints whose drag or local force-work identity is lost to half-force
/// cancellation refuse; neither resistance nor velocity is clipped.
pub struct ResistanceGrid3<L: ForcedLattice3<Q>, const Q: usize> {
    flow: ForcedGrid3<L, Q>,
    laws: Vec<DarcyForchheimer3>,
    next_forces: Vec<[f64; 3]>,
    drag: Vec<[f64; 3]>,
    next_drag: Vec<[f64; 3]>,
    walls: WallMomentumProbe3<L, Q>,
    origin: [f64; 3],
}

impl<L: ForcedLattice3<Q>, const Q: usize> ResistanceGrid3<L, Q> {
    /// Attach a complete immutable law field to an UNSTEPPED zero-uniform-force
    /// grid. One law is supplied per cell in x-fastest order; solid entries are
    /// admitted but apply no force. All wall moments use the finite `origin`.
    ///
    /// This declares the INITIAL population frame, as for `ForcedGrid3::try_new`;
    /// it does not preserve a previously read unforced velocity. For each fluid
    /// cell the initial force is solved from the supplied raw populations and
    /// law. Read the initial physical velocity from `flow().macroscopic`.
    /// Both original raw-grid and final force-corrected Mach admission apply.
    /// No force rebinding, equilibrium reset or elapsed impulse is performed.
    pub fn try_new(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        laws: Vec<DarcyForchheimer3>,
        origin: [f64; 3],
    ) -> Result<Self, GridError3<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        if grid.completed_steps != 0 || grid.force.iter().any(|&f| f != 0.0) {
            return Err(GridError3::Input(
                "resistance requires an unstepped zero-force grid",
            ));
        }
        if laws.len() != grid.f.len() {
            return Err(GridError3::Input("resistance field length"));
        }
        let walls = grid.wall_momentum_probe(cx, origin, |_| true)?;
        let mut forces = Vec::new();
        forces
            .try_reserve_exact(laws.len())
            .map_err(|_| GridError3::Allocation)?;
        let mut drag = Vec::new();
        drag.try_reserve_exact(laws.len())
            .map_err(|_| GridError3::Allocation)?;
        for (cell, &law) in laws.iter().enumerate() {
            if cell.is_multiple_of(64) {
                cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
            }
            law.admit().map_err(|error| GridError3::MaterialField {
                cell,
                reason: match error {
                    BoundaryStepError3::EndpointForce { reason } => reason,
                    _ => "invalid resistance law",
                },
            })?;
            let (force, reaction) = if grid.solid[cell] {
                ([0.0; 3], [0.0; 3])
            } else {
                let (rho, momentum) = raw_moments_cell::<L, Q>(&grid.f[cell], cell)?;
                let (force, drag) = law
                    .endpoint_force(rho, momentum)
                    .ok_or(GridError3::NonFiniteCell { cell })?;
                let (_, velocity) =
                    macro_cell::<L, Q>(&grid.f[cell], force, cell, grid.budget.maximum_mach)?;
                admit_drag::<L, Q>(law, velocity, drag, cell)?;
                (force, drag)
            };
            forces.push(force);
            drag.push(reaction);
        }
        let mut next_forces = Vec::new();
        next_forces
            .try_reserve_exact(laws.len())
            .map_err(|_| GridError3::Allocation)?;
        next_forces.resize(laws.len(), [0.0; 3]);
        let mut next_drag = Vec::new();
        next_drag
            .try_reserve_exact(laws.len())
            .map_err(|_| GridError3::Allocation)?;
        next_drag.resize(laws.len(), [0.0; 3]);
        let flow = ForcedGrid3::try_new(grid, cx, forces)?;
        Ok(Self {
            flow,
            laws,
            next_forces,
            drag,
            next_drag,
            walls,
            origin,
        })
    }

    /// Read-only native populations, force history, macros and elapsed steps.
    pub const fn flow(&self) -> &ForcedGrid3<L, Q> {
        &self.flow
    }
    /// Held constitutive coefficients in the grid's x-fastest cell order.
    pub fn laws(&self) -> &[DarcyForchheimer3] {
        &self.laws
    }
    /// Fixed moment origin in lattice cell coordinates.
    pub const fn moment_origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Advance with the native collision and implicit porous endpoint. An
    /// optional collision field has the original material-viscosity meaning;
    /// None retains the constructor's uniform collision family and rates.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        self.step_inner::<fs_exec::TilePool>(None, collisions, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    /// Same fixed64-cell arithmetic and ordered receipts, independent of the
    /// executor's worker count. Reuses all native population and force buffers.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        self.step_inner(Some((pool, gate)), collisions, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })
    }

    fn step_inner<P: KernelRunner>(
        &mut self,
        pool: Option<(&P, &CancelGate)>,
        collisions: Option<&[L::Collision]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        checkpoint()?;
        let grid = &mut self.flow.grid;
        if let Some(field) = collisions {
            grid.admit_collisions(field, &mut checkpoint)?;
        }
        super::admit_forces(grid, &self.flow.forces, collisions, &mut checkpoint)?;
        let step = grid
            .completed_steps
            .checked_add(1)
            .ok_or(GridError3::Input("step count overflow"))?;
        let data = StepData::<L, Q> {
            dims: grid.dims,
            periodic: grid.periodic,
            solid: &grid.solid,
            collision: grid.collision,
            collisions,
            force: grid.force,
            forces: Some(&self.flow.forces),
            budget: grid.budget,
            marker: PhantomData,
        };
        if let Some((pool, gate)) = pool {
            run_pass(&data, &grid.f, &mut grid.post, false, pool, gate)?;
        } else {
            for (tile, chunk) in grid.post.chunks_mut(64).enumerate() {
                checkpoint()?;
                data.process(&grid.f, chunk, tile * 64, false)?;
            }
        }
        let mut receipt = ResistanceMomentumExchange3::default();
        {
            let kernel = EndpointPass {
                data: &data,
                initial: &grid.f,
                post: &grid.post,
                laws: &self.laws,
                initial_drag: &self.drag,
                origin: self.origin,
                chunks: grid
                    .next
                    .chunks_mut(64)
                    .zip(self.next_forces.chunks_mut(64))
                    .zip(self.next_drag.chunks_mut(64))
                    .map(|((populations, forces), drag)| {
                        Mutex::new(EndpointChunk {
                            populations,
                            forces,
                            drag,
                            receipt: ResistanceMomentumExchange3::default(),
                        })
                    })
                    .collect(),
            };
            if let Some((pool, gate)) = pool {
                match pool.run_with_gate(&kernel, gate).0 {
                    Ok(FirstError(None)) => {}
                    Ok(FirstError(Some(error))) => return Err(error),
                    Err(_) if gate.is_requested() => return Err(GridError3::Cancelled),
                    Err(error) => return Err(GridError3::Pool(format!("{error:?}"))),
                }
            } else {
                for tile in 0..kernel.chunks.len() {
                    checkpoint()?;
                    kernel.run_tile(tile)?;
                }
            }
            for chunk in &kernel.chunks {
                let contribution = chunk.lock().map_err(|_| GridError3::Poisoned)?.receipt;
                receipt.add(contribution);
                for a in 0..3 {
                    receipt.periodic.angular_momentum_change[a] +=
                        contribution.periodic.angular_momentum_change[a];
                }
                receipt.periodic.face_crossings += contribution.periodic.face_crossings;
            }
        }
        if !receipt.is_finite()
            || !receipt
                .periodic
                .angular_momentum_change
                .iter()
                .all(|v| v.is_finite())
        {
            return Err(GridError3::Input("resistance impulse or work overflow"));
        }
        let walls = self.walls.measure(grid, &mut checkpoint)?;
        receipt.walls = WallMomentumExchange3 {
            wall_impulse: walls.wall_impulse,
            fluid_population_impulse: walls.wall_impulse.map(|j| -j),
            wall_angular_impulse: walls.wall_angular_impulse,
            fluid_population_angular_impulse: walls.wall_angular_impulse.map(|j| -j),
            measured_links: walls.measured_links,
            ..WallMomentumExchange3::default()
        };
        checkpoint()?;
        core::mem::swap(&mut grid.f, &mut grid.next);
        core::mem::swap(&mut self.flow.forces, &mut self.next_forces);
        core::mem::swap(&mut self.drag, &mut self.next_drag);
        grid.completed_steps = step;
        Ok(receipt)
    }
}

struct EndpointChunk<'a, const Q: usize> {
    populations: &'a mut [[f64; Q]],
    forces: &'a mut [[f64; 3]],
    drag: &'a mut [[f64; 3]],
    receipt: ResistanceMomentumExchange3,
}

struct EndpointPass<'a, 'b, L: ForcedLattice3<Q>, const Q: usize> {
    data: &'a StepData<'b, L, Q>,
    initial: &'a [[f64; Q]],
    post: &'a [[f64; Q]],
    laws: &'a [DarcyForchheimer3],
    initial_drag: &'a [[f64; 3]],
    origin: [f64; 3],
    chunks: Vec<Mutex<EndpointChunk<'a, Q>>>,
}

impl<L: ForcedLattice3<Q>, const Q: usize> EndpointPass<'_, '_, L, Q> {
    fn run_tile(&self, tile: usize) -> Result<(), GridError3<L::Error>> {
        let mut chunk = self.chunks[tile].lock().map_err(|_| GridError3::Poisoned)?;
        let EndpointChunk {
            populations,
            forces,
            drag: stored_drag,
            receipt,
        } = &mut *chunk;
        let base = tile * 64;
        self.data.process_stream_with_endpoint::<true>(
            self.post,
            populations,
            base,
            |cell, pop, periodic| {
                let (rho, momentum) = raw_moments_cell::<L, Q>(pop, cell)?;
                let law = self.laws[cell];
                let (force, drag) = law
                    .endpoint_force(rho, momentum)
                    .ok_or(GridError3::NonFiniteCell { cell })?;
                let collision = self
                    .data
                    .collisions
                    .map_or(self.data.collision, |field| field[cell]);
                if force != [0.0; 3] && !L::supports_body_force(collision) {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "collision model does not admit body force",
                    });
                }
                let (_, velocity) =
                    macro_cell::<L, Q>(pop, force, cell, self.data.budget.maximum_mach)?;
                admit_drag::<L, Q>(law, velocity, drag, cell)?;
                let previous = self
                    .data
                    .forces
                    .map_or(self.data.force, |field| field[cell]);
                let (_, initial_velocity) = macro_cell::<L, Q>(
                    &self.initial[cell],
                    previous,
                    cell,
                    self.data.budget.maximum_mach,
                )?;
                let contribution = law
                    .momentum_exchange(
                        initial_velocity,
                        velocity,
                        self.initial_drag[cell],
                        drag,
                        coordinates(cell, self.data.dims),
                        self.origin,
                    )
                    .ok_or(GridError3::NonFiniteCell { cell })?;
                let defect = contribution.fluid_drag_work
                    + contribution.skeleton_work
                    + contribution.relative_dissipation;
                // Scale the check before subtracting opposing endpoint work.
                // Stiff trapezoidal reversal can make their net much smaller.
                let scale = contribution.relative_dissipation
                    + (0..3)
                        .map(|a| {
                            let midpoint = 0.5 * initial_velocity[a] + 0.5 * velocity[a];
                            (0.5 * self.initial_drag[cell][a].abs() + 0.5 * drag[a].abs())
                                * (midpoint.abs() + law.linear.skeleton_velocity[a].abs())
                        })
                        .sum::<f64>();
                if !scale.is_finite() || defect.abs() > 4096.0 * f64::EPSILON * scale {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "resistance work is unresolved in the represented velocity",
                    });
                }
                forces[cell - base] = force;
                stored_drag[cell - base] = drag;
                receipt.add(contribution);
                for a in 0..3 {
                    receipt.periodic.angular_momentum_change[a] +=
                        periodic.angular_momentum_change[a];
                }
                receipt.periodic.face_crossings += periodic.face_crossings;
                Ok(())
            },
        )
    }
}

// A stiff half-force cancellation can make the true endpoint velocity smaller
// than its stored representation resolves. Do not then substitute a different
// constitutive reaction for the force actually applied by collision. Accepted
// drag is retained verbatim, and an unresolved constitutive state refuses.
fn admit_drag<L: ForcedLattice3<Q>, const Q: usize>(
    law: DarcyForchheimer3,
    velocity: [f64; 3],
    drag: [f64; 3],
    cell: usize,
) -> Result<(), GridError3<L::Error>> {
    let represented = law.drag(velocity);
    let scale = drag
        .iter()
        .chain(&represented)
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    if !represented.iter().all(|v| v.is_finite())
        || (0..3).any(|a| (drag[a] - represented[a]).abs() > 4096.0 * f64::EPSILON * scale)
    {
        return Err(GridError3::MaterialField {
            cell,
            reason: "resistance force is unresolved in the represented velocity",
        });
    }
    Ok(())
}

impl<L: ForcedLattice3<Q>, const Q: usize> TileKernel for EndpointPass<'_, '_, L, Q> {
    type Out = FirstError<L::Error>;
    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm/grid3-resistance-stream", self.chunks.len() as u64)
    }
    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        ControlFlow::Continue(FirstError(self.run_tile(tile as usize).err()))
    }
}
