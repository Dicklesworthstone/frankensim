//! Implicit Darcy–Forchheimer flow through native velocity/pressure reservoirs.
//!
//! The original constitutive solve, collision, regularized reservoir and moving
//! wall laws are reused. Force/drag endpoints, populations and complete momentum
//! accounting publish together. Held steps and authored next-endpoint laws
//! retain their own force history and temporal work meaning.

use super::*;
use super::super::{EndpointSource3, EndpointWorkspace3, FaceBoundary3, StepData, coordinates,
    WallVelocityProfileReport3, wall_profile_values};
use crate::d3q19::BoundaryStepError3;
use crate::lattice3::moving_walls::WallExchange;
use crate::lattice3::sources::resistance::{SourceEndpoint3, admit_drag, held_exchange};
use crate::lattice3::sources::resistance::endpoint as source_endpoint;
use crate::lattice3::wall_loads::WallMomentumProbe3;
use core::marker::PhantomData;
use fs_exec::KernelRunner;

/// Original D3Q7 thermal transport with this implicit porous velocity owner.
pub mod thermal;

mod endpoint;
pub use endpoint::ResistanceOpenEndpointExchange3;

pub use crate::d3q19::{DarcyForchheimer3, LinearResistance3, ResistanceMomentumExchange3};

/// A porous medium with native D3Q19/D3Q27 through-flow and moving walls.
///
/// At each fluid cell `F = b - D (u-v_s) - beta |u-v_s| (u-v_s)`. Velocities
/// and material axes are resolved in the grid frame; `v_s` is prescribed
/// skeleton velocity in that SAME frame, not the wall velocity or fluid slip.
/// Coefficients are lattice force-density/velocity and force-density/velocity²,
/// already including any physical constitutive density factor. No extra rho is
/// applied. A free direction has zero resistance; infinite resistance refuses.
///
/// Accepted forces drive the original collision. The native pull and original
/// implicit law resolve endpoint force from `rho*u = m_raw + F/2`. Full face
/// traces use their prescribed physical velocity and law-derived half-source;
/// donors use their own law and density. Pressure/velocity profiles and optional
/// characteristic outlets retain their existing ownership and admission.
///
/// Nonzero drag on a mixed wall/open boundary cell refuses at attachment.
/// Such partially replaced populations do not prescribe a complete physical
/// velocity; this owner does not silently use that nominal trace for drag.
/// Zero-drag seams retain the original force-density path. Interior pores,
/// stationary voxel walls and admitted tangential exterior motion are supported.
///
/// Every interval returns the original constitutive force/work receipt plus
/// actual reservoir transport, all wall links and periodic angular transport.
/// `step_*` holds the accepted law; `advance_*` supplies its next physical
/// endpoint, including changing skeleton velocity, resistance and body force.
/// The held-law identity is `fluid_drag_work + skeleton_work = -relative_dissipation`;
/// it is not a transported fluid energy balance. The skeleton has prescribed
/// motion and fixed geometry; no porosity, skeleton inertia or drag heat evolves.
/// Original low-Mach/collision budgets and unresolved-force refusals apply.
pub struct ResistanceOpenGrid3<L: EquilibriumLattice3<Q>, const Q: usize> {
    flow: ForcedOpenBoundaryGrid3<L, Q>,
    laws: Vec<DarcyForchheimer3>,
    drag: Vec<[f64; 3]>,
    pending_drag: Vec<[f64; 3]>,
    walls: WallMomentumProbe3<L, Q>,
    origin: [f64; 3],
    last_exchange: Option<ResistanceMomentumExchange3>,
    last_endpoint_exchange: Option<ResistanceOpenEndpointExchange3>,
}

fn workspace<E>(cells: usize) -> Result<Vec<[f64; 3]>, GridError3<E>> {
    let mut field = Vec::new();
    field.try_reserve_exact(cells).map_err(|_| GridError3::Allocation)?;
    field.resize(cells, [0.0; 3]);
    Ok(field)
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ResistanceOpenGrid3<L, Q> {
    /// Attach one admitted law per x-fast cell without changing physical
    /// density, velocity, force-corrected stress, elapsed time or controls.
    /// Solid law entries are validated but apply no force or drag.
    ///
    /// This uses the native source-coordinate transform, not a momentum kick.
    /// The original background collision selects the initial source family.
    /// Scratch and the actual wall-link selection are allocated once here;
    /// serial steps invoke no user callbacks and allocate no new source storage.
    pub fn try_new_preserving_velocity(
        inner: OpenBoundaryGrid3<L, Q>,
        cx: &Cx<'_>,
        laws: Vec<DarcyForchheimer3>,
        origin: [f64; 3],
    ) -> Result<Self, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        poll()?;
        if laws.len() != inner.grid.f.len() {
            return Err(GridError3::Input("resistance field length"));
        }
        let plan = inner.plan;
        let walls = inner.grid.wall_momentum_probe(cx, origin, |link| {
            link.solid_cell.is_some()
                || link.exterior_faces.iter().enumerate().any(|(face, &crossed)| {
                    crossed && matches!(plan.boundaries.face(Face3::ALL[face]), FaceBoundary3::Wall { .. })
                })
        })?;
        let mut forces = workspace(laws.len())?;
        let mut drag = workspace(laws.len())?;
        let data = StepData::<L, Q> {
            dims: inner.grid.dims,
            periodic: inner.grid.periodic,
            solid: &inner.grid.solid,
            collision: inner.grid.collision,
            collisions: None,
            force: inner.grid.force,
            forces: None,
            budget: inner.grid.budget,
            marker: PhantomData,
        };
        for (cell, &law) in laws.iter().enumerate() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            law.admit().map_err(|error| GridError3::MaterialField {
                cell,
                reason: match error {
                    BoundaryStepError3::EndpointForce { reason } => reason,
                    _ => "invalid resistance law",
                },
            })?;
            if inner.grid.solid[cell] {
                continue;
            }
            if (law.linear.principal != [0.0; 3] || law.quadratic != 0.0)
                && plan.locate(cell).is_some()
            {
                let (wall, open) = plan.masks(&data, cell);
                if wall != 0 && open != 0 {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "active resistance requires a complete open-face trace; wall/open seams are unsupported",
                    });
                }
            }
            let (_, velocity) = inner.macroscopic(cell)?;
            forces[cell] = law.force(velocity);
            drag[cell] = law.drag(velocity);
        }
        let pending_drag = workspace(laws.len())?;
        let pending_forces = workspace(laws.len())?;
        let mut flow = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
            inner, forces, &mut poll,
        )?;
        for (cell, &law) in laws.iter().enumerate() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            if !flow.inner.grid.solid[cell] {
                let (_, velocity) = flow.macroscopic(cell)?;
                admit_drag::<L, Q>(law, velocity, drag[cell], cell)?;
            }
        }
        poll()?;
        flow.pending_forces = pending_forces;
        Ok(Self {
            flow, laws, drag, pending_drag, walls, origin, last_exchange: None,
            last_endpoint_exchange: None,
        })
    }

    /// Accepted native populations, forces, physical macros and boundary traces.
    /// A mutable force-only owner would permit stepping with stale drag history.
    pub const fn flow(&self) -> &ForcedOpenBoundaryGrid3<L, Q> {
        &self.flow
    }

    /// Constitutive field at the same accepted endpoint as force and velocity.
    pub fn laws(&self) -> &[DarcyForchheimer3] {
        &self.laws
    }

    /// Fixed moment origin in lattice cell coordinates.
    pub const fn moment_origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Held-law receipt from the last accepted interval, or None when that
    /// interval used `advance_*`. Failures retain both receipt selections.
    pub const fn last_exchange(&self) -> Option<ResistanceMomentumExchange3> {
        self.last_exchange
    }

    /// Original native collision/pull with an implicit porous endpoint.
    /// An optional spatial collision field keeps its original viscosity meaning.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        self.step_inner::<fs_exec::TilePool>(None, collisions, &mut || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    /// Same original tile arithmetic and canonical receipt order for any worker
    /// count. Source buffers are reused; scheduling remains owned by fs-exec.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        self.step_inner(Some((pool, gate)), collisions, &mut || {
            if gate.is_requested() { Err(GridError3::Cancelled) } else { Ok(()) }
        })
    }

    fn step_inner<P: KernelRunner>(
        &mut self,
        pool: Option<(&P, &CancelGate)>,
        collisions: Option<&[L::Collision]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ResistanceMomentumExchange3, GridError3<L::Error>> {
        let (receipt, step) = self.stage_inner(pool, collisions, poll)?;
        poll()?;
        self.publish_step(step, receipt);
        Ok(receipt)
    }

    // Composing heat/material owners admit their candidate before any physical
    // momentum, force, drag or receipt epoch is published.
    fn stage_inner<P: KernelRunner>(
        &mut self,
        pool: Option<(&P, &CancelGate)>,
        collisions: Option<&[L::Collision]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(ResistanceMomentumExchange3, usize), GridError3<L::Error>> {
        let (held, _, step) = self.stage_with_laws(pool, collisions, None, poll)?;
        Ok((held, step))
    }

    fn stage_with_laws<P: KernelRunner>(
        &mut self,
        pool: Option<(&P, &CancelGate)>,
        collisions: Option<&[L::Collision]>,
        next_laws: Option<(&[DarcyForchheimer3], f64)>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(ResistanceMomentumExchange3, ResistanceOpenEndpointExchange3, usize), GridError3<L::Error>> {
        poll()?;
        if let Some((laws, limit)) = next_laws {
            self.admit_next_laws(laws, limit, poll)?;
        }
        let endpoint = Some((
            EndpointSource3::Resistance(next_laws.map_or(self.laws.as_slice(), |(laws, _)| laws)),
            EndpointWorkspace3 {
                forces: &mut self.flow.pending_forces,
                drag: Some(&mut self.pending_drag),
            },
        ));
        let (open, step) = if let Some((pool, gate)) = pool {
            self.flow.inner.stage_pooled_with_source(
                pool, gate, self.flow.inner.plan, collisions,
                Some(&self.flow.forces), None, endpoint, self.origin,
            )?
        } else {
            self.flow.inner.stage_serial_with_source(
                self.flow.inner.plan, collisions, Some(&self.flow.forces),
                None, endpoint, self.origin, poll,
            )?
        };
        let (held, evolving) = self.measure(open, collisions, next_laws.map(|(laws, _)| laws), poll)?;
        if next_laws.is_some_and(|(_, limit)| {
            evolving.resistance.absolute_relative_quadrature_remainder > limit
        }) {
            return Err(GridError3::Input("resistance temporal work budget exceeded"));
        }
        Ok((held, evolving, step))
    }

    // The complete native candidate and every composing gate have succeeded.
    // Publication only swaps admitted buffers and cannot fail.
    fn publish_step(&mut self, step: usize, receipt: ResistanceMomentumExchange3) {
        self.publish_fields(step);
        self.last_exchange = Some(receipt);
        self.last_endpoint_exchange = None;
    }

    fn publish_fields(&mut self, step: usize) {
        core::mem::swap(&mut self.flow.inner.grid.f, &mut self.flow.inner.grid.next);
        core::mem::swap(&mut self.flow.forces, &mut self.flow.pending_forces);
        core::mem::swap(&mut self.drag, &mut self.pending_drag);
        self.flow.inner.grid.completed_steps = step;
        self.flow.inner.publish_wall_exchange();
    }

    fn measure(
        &self,
        open: OpenBoundaryMomentumTransport3,
        collisions: Option<&[L::Collision]>,
        next_laws: Option<&[DarcyForchheimer3]>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(ResistanceMomentumExchange3, ResistanceOpenEndpointExchange3), GridError3<L::Error>> {
        let grid = &self.flow.inner.grid;
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
        let mut receipt = ResistanceMomentumExchange3 {
            open, ..ResistanceMomentumExchange3::default()
        };
        let mut evolving = ResistanceOpenEndpointExchange3 {
            open, ..ResistanceOpenEndpointExchange3::default()
        };
        for (cell, &law) in self.laws.iter().enumerate() {
            if cell.is_multiple_of(64) {
                poll()?;
            }
            if grid.solid[cell] {
                continue;
            }
            let (_, initial_velocity) = self.flow.macroscopic(cell)?;
            let (_, velocity) = macro_cell::<L, Q>(
                &grid.next[cell], self.flow.pending_forces[cell], cell,
                grid.budget.maximum_mach,
            )?;
            let initial = SourceEndpoint3 { law, velocity: initial_velocity, drag: self.drag[cell] };
            let end = SourceEndpoint3 {
                law: next_laws.map_or(law, |field| field[cell]),
                velocity, drag: self.pending_drag[cell],
            };
            if next_laws.is_some() {
                evolving.resistance.add(source_endpoint::exchange::<L, Q>(
                    initial, end, cell, coordinates(cell, grid.dims), self.origin,
                )?);
            } else {
                receipt.add(held_exchange::<L, Q>(
                    initial, end, cell, coordinates(cell, grid.dims), self.origin,
                )?);
            }
            // Read the original pull for actual periodic coordinate crossings.
            // Reservoir replacement relative to this pull is already in open;
            // no second collision, source application or state evolution occurs.
            if grid.periodic.iter().any(|&periodic| periodic) {
                let mut pulled = [0.0; Q];
                let periodic = data.pull_cell_with_periodic::<true>(&grid.post, &mut pulled, cell);
                for a in 0..3 {
                    receipt.periodic.angular_momentum_change[a] += periodic.angular_momentum_change[a];
                }
                receipt.periodic.face_crossings += periodic.face_crossings;
            }
        }
        let mut walls = WallExchange::default();
        let wall_profiles = wall_profile_values(&self.flow.inner.wall_profiles);
        for (index, link) in self.walls.links.iter().enumerate() {
            if index.is_multiple_of(64) {
                poll()?;
            }
            let q = link.incoming;
            let velocity = if link.solid.is_some() {
                [0.0; 3]
            } else {
                self.flow.inner.plan.exterior_wall_velocity::<L::Error, Q>(
                    &wall_profiles,
                    coordinates(link.cell, grid.dims), L::E[q], q,
                )?.ok_or(GridError3::Input("resistance wall probe lost wall ownership"))?
            };
            walls.add(
                link.cell, grid.post[link.cell][L::OPP[q]], grid.next[link.cell][q],
                L::E[q].map(f64::from), velocity, link.arm,
            )?;
        }
        receipt.walls = walls.finish()?;
        evolving.resistance.walls = receipt.walls;
        evolving.resistance.periodic = receipt.periodic;
        if !receipt.is_finite()
            || (next_laws.is_some() && !evolving.resistance.is_finite())
            || receipt.periodic.angular_momentum_change.iter().any(|v| !v.is_finite())
        {
            return Err(GridError3::Input("resistance impulse or work overflow"));
        }
        poll()?;
        Ok((receipt, evolving))
    }

    /// Update prescribed reservoir or tangential wall values atomically.
    /// Face kinds, accepted fluid/source state and immutable laws remain held.
    pub fn rebind_open_boundaries(
        &mut self,
        boundaries: BoundarySpec3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.flow.rebind_open_boundaries(boundaries, gate)
    }

    /// Drive the same actual halfway-wall links used by porous streaming and
    /// its complete wall/skeleton/open receipt. Accepted force history is held.
    pub fn set_wall_velocity_profile(
        &mut self, face: Face3, maximum_storage_bytes: usize, gate: &CancelGate,
        sample: impl FnMut([f64; 3]) -> Result<[f64; 3], GridError3<L::Error>>,
    ) -> Result<WallVelocityProfileReport3, GridError3<L::Error>> {
        self.flow.set_wall_velocity_profile(face, maximum_storage_bytes, gate, sample)
    }

    pub fn clear_wall_velocity_profile(
        &mut self, face: Face3, gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.flow.clear_wall_velocity_profile(face, gate)
    }

    pub fn wall_velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.flow.wall_velocity_profile(face)
    }

    /// Resolve a spatial/time-updated inlet before the next native interval.
    pub fn set_velocity_profile(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        self.flow.set_velocity_profile(face, flux, gate, sample)
    }

    pub fn set_velocity_profile_with_speed_limit(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        maximum_speed_lu: f64,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        self.flow.set_velocity_profile_with_speed_limit(face, flux, maximum_speed_lu, gate, sample)
    }

    pub fn clear_velocity_profile(&mut self, face: Face3) {
        self.flow.clear_velocity_profile(face);
    }

    pub fn set_pressure_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<f64, BoundaryProfileError3>,
    ) -> Result<PressureProfileReport3, BoundaryProfileError3> {
        self.flow.set_pressure_profile(face, gate, sample)
    }

    pub fn clear_pressure_profile(&mut self, face: Face3) {
        self.flow.clear_pressure_profile(face);
    }

    pub fn set_characteristic_outlet(
        &mut self,
        face: Face3,
        control: CharacteristicOutlet3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.flow.set_characteristic_outlet(face, control, gate)
    }

    pub fn clear_characteristic_outlet(&mut self, face: Face3) {
        self.flow.clear_characteristic_outlet(face);
    }
}
