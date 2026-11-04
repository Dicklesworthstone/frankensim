//! Finite-inertia wall response on the original shared D3Q19/D3Q27 stream.
//!
//! A fixed planar coupling window admits two tangential velocity coordinates
//! and spin about its positive Cartesian normal axis. Its geometry and inertia
//! are frozen during this exchange; no translating footprint, pose integrator,
//! fresh-cell initialization or moving-volume remap is implied. Normal force
//! and transverse torque are carried by a zero-work support.
//!
//! Complete face-cell stencils give the actual wall generalized impulse
//! `j(v) = b - K v`, with `K = sum(2*w*rho/cs² * g*gᵀ)` positive semidefinite.
//! Here `g` projects an incoming direction onto the three wall modes, and
//! `b = -sum(2*f_out*g)` uses the native post-collision populations. Therefore
//! `(M + K/2) delta_v = b - K*v_before + j_external` closes the midpoint
//! response without repeating collision, iteration or a force extrapolation.
//!
//! Wall velocity is sampled once per face cell. For complete D3Q19/D3Q27
//! face stencils, moving mass sums to zero and their tangential second moment
//! is isotropic. The boundary-relative load therefore has this same linear
//! force/axial-torque response, including Wen's moving-mass terms. The final
//! measured populations and work must still close the declared budgets before
//! either fluid or wall state is accepted. Fluid kinetic-energy loss is not
//! identified as wall heat; this remains an isothermal LBM exchange.

use super::moving_walls::{MovingWallMomentum3, WallMotion3};
use super::wall_loads::{StationaryWallMomentum3, WallMomentumProbe3};
use super::{Grid3, GridError3, Lattice3, StepData, coordinates, run_pass_with_wall_motion};
use crate::d3q19::Face3;
use core::marker::PhantomData;
use fs_exec::{CancelGate, Cx, KernelRunner};
use fs_la::LuWorkspace;

/// Accepted planar mechanical state in lattice units (`dx = dt = 1`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarWallState3 {
    /// Positive mass for either tangential translation coordinate.
    pub mass: f64,
    /// Positive inertia about the fixed positive Cartesian normal axis.
    pub axial_inertia: f64,
    /// Two tangential speeds followed by normal-axis angular speed. For normal
    /// axis `n`, the tangential axes are `(n+1)%3` and `(n+2)%3`, in that order.
    pub velocity: [f64; 3],
}

/// Numerical acceptance, separate from the physical mass, inertia and drag.
#[derive(Clone, Copy, Debug)]
pub struct FiniteWallBudget3 {
    /// Each momentum defect divided by the sum of its inertia, fluid and drive
    /// impulse magnitudes. A zero denominator requires an exact zero defect.
    pub maximum_relative_impulse_defect: f64,
    /// Work defect divided by the sum of kinetic increment, fluid work and
    /// external work magnitudes, with the same exact-zero convention.
    pub maximum_relative_work_defect: f64,
}

/// Prepared fixed planar coupling window and its finite mechanical state.
///
/// Preparation binds all wall topology and allocates sparse link selections,
/// one reusable velocity field and the existing three-by-three LU workspace.
/// Stepping adds no population buffer or per-step wall-field/factor allocation.
/// The original executor retains its own ordinary tile scratch allocation.
///
/// Translation describes material moving tangentially through this fixed
/// coupling window. A translating finite footprint or changing body geometry
/// needs a separate conservative geometry update. The caller owns any pose,
/// gyroscopic response, external mechanism or evolution of inertia.
pub struct FinitePlanarWall3<L: Lattice3<Q>, const Q: usize> {
    face: Face3,
    center: [f64; 3],
    state: PlanarWallState3,
    solid: Vec<bool>,
    moving: WallMomentumProbe3<L, Q>,
    stationary: WallMomentumProbe3<L, Q>,
    velocities: Vec<[f64; 3]>,
    solver: LuWorkspace,
}

/// Accepted mechanical response and separate native fluid/support exchanges.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePlanarWallStep3 {
    pub before: PlanarWallState3,
    pub after: PlanarWallState3,
    pub velocity_midpoint: [f64; 3],
    /// Actual impulse, angular impulse and work ON the finite wall window.
    /// Angular impulses use its fixed center as the moment origin.
    pub hydrodynamic: MovingWallMomentum3,
    /// All other exterior and voxel walls, about that same center.
    pub stationary_walls: StationaryWallMomentum3,
    /// Reaction ON the finite wall, cancelling its normal hydrodynamic force.
    pub support_impulse: [f64; 3],
    /// Reaction ON the wall, cancelling transverse torque about its center.
    pub support_angular_impulse: [f64; 3],
    /// Two supplied tangential impulses and one axial angular impulse.
    pub external_impulse: [f64; 3],
    /// Held Guo fluid-force impulse, separate from wall and drive reactions.
    pub fluid_force_impulse: [f64; 3],
    pub kinetic_energy_change: f64,
    pub external_work: f64,
    pub momentum_residual: [f64; 3],
    /// Wall kinetic increment minus measured fluid work minus external work.
    pub work_residual: f64,
}

impl<L: Lattice3<Q>, const Q: usize> FinitePlanarWall3<L, Q> {
    /// Select complete cells on one nonperiodic exterior face. The predicate is
    /// called once per face cell, in native cell order. Side-wall corner links
    /// refuse; use periodic tangents or a window away from those corners.
    /// The center lies on the physical wall plane, at 0 or the grid dimension.
    /// GridBudget3 bounds traversal/storage; all allocation failure is explicit.
    pub fn new(
        grid: &Grid3<L, Q>,
        cx: &Cx<'_>,
        face: Face3,
        center: [f64; 3],
        state: PlanarWallState3,
        mut select_cell: impl FnMut([usize; 3]) -> bool,
    ) -> Result<Self, GridError3<L::Error>> {
        let axis = face as usize / 2;
        let plane = if (face as usize).is_multiple_of(2) {
            0.0
        } else {
            grid.dims[axis] as f64
        };
        if grid.periodic[axis]
            || !center.iter().all(|x| x.is_finite())
            || center[axis] != plane
            || !state.mass.is_finite()
            || state.mass <= 0.0
            || !state.axial_inertia.is_finite()
            || state.axial_inertia <= 0.0
            || !state.velocity.iter().all(|x| x.is_finite())
        {
            return Err(GridError3::Input("finite wall face, center or inertia"));
        }
        let mut last = None;
        let mut selected = false;
        let moving = grid.wall_momentum_probe(cx, center, |link| {
            if !link.exterior_faces[face as usize] {
                return false;
            }
            if last != Some(link.fluid_cell) {
                last = Some(link.fluid_cell);
                selected = select_cell(link.fluid_cell);
            }
            selected
        })?;
        if moving.links.is_empty() {
            return Err(GridError3::Input(
                "finite wall needs a nonempty complete planar window",
            ));
        }
        let sign = if (face as usize).is_multiple_of(2) {
            1
        } else {
            -1
        };
        let expected = L::E.iter().filter(|e| e[axis] == sign).count();
        let mut count = 0;
        for (i, (link, descriptor)) in moving.links.iter().zip(moving.wall_links()).enumerate() {
            if i % 64 == 0 {
                cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
            }
            if descriptor.solid_cell.is_some()
                || descriptor.exterior_faces.iter().filter(|&&x| x).count() != 1
                || L::E[link.incoming][axis] != sign
            {
                return Err(GridError3::Input(
                    "finite wall window crosses another surface",
                ));
            }
            count += 1;
            if i + 1 == moving.links.len() || moving.links[i + 1].cell != link.cell {
                if count != expected {
                    return Err(GridError3::Input(
                        "finite wall requires complete face-cell stencils",
                    ));
                }
                count = 0;
            }
        }
        let stationary = grid.wall_momentum_probe(cx, center, |link| {
            let p = link.fluid_cell;
            let cell = (p[2] * grid.dims[1] + p[1]) * grid.dims[0] + p[0];
            !link.exterior_faces[face as usize]
                || moving
                    .links
                    .binary_search_by_key(&cell, |link| link.cell)
                    .is_err()
        })?;
        let mut velocities = Vec::new();
        velocities
            .try_reserve_exact(moving.links.len())
            .map_err(|_| GridError3::Allocation)?;
        velocities.resize(moving.links.len(), [0.0; 3]);
        let mut solid = Vec::new();
        solid
            .try_reserve_exact(grid.solid.len())
            .map_err(|_| GridError3::Allocation)?;
        solid.extend_from_slice(&grid.solid);
        let mut out = Self {
            face,
            center,
            state,
            solid,
            moving,
            stationary,
            velocities,
            solver: LuWorkspace::new(3).map_err(|_| GridError3::Allocation)?,
        };
        out.fill_velocities(grid.dims, state.velocity);
        out.motion()
            .admit(grid, || cx.checkpoint().map_err(|_| GridError3::Cancelled))?;
        Ok(out)
    }

    pub const fn state(&self) -> PlanarWallState3 {
        self.state
    }
    pub const fn face(&self) -> Face3 {
        self.face
    }
    pub const fn center(&self) -> [f64; 3] {
        self.center
    }
    pub fn measured_links(&self) -> usize {
        self.moving.measured_links()
    }

    fn motion(&self) -> WallMotion3<'_, L, Q> {
        WallMotion3 {
            probe: &self.moving,
            velocities: &self.velocities,
        }
    }

    fn fill_velocities(&mut self, dims: [usize; 3], v: [f64; 3]) {
        let n = self.face as usize / 2;
        let a = (n + 1) % 3;
        let b = (n + 2) % 3;
        for (link, u) in self.moving.links.iter().zip(&mut self.velocities) {
            let p = coordinates(link.cell, dims);
            *u = [0.0; 3];
            u[a] = v[0] - v[2] * (p[b] as f64 + 0.5 - self.center[b]);
            u[b] = v[1] + v[2] * (p[a] as f64 + 0.5 - self.center[a]);
        }
    }

    fn solve(
        &mut self,
        post: &[[f64; Q]],
        dims: [usize; 3],
        external: [f64; 3],
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<([f64; 3], [f64; 3]), GridError3<L::Error>> {
        let n = self.face as usize / 2;
        let a = (n + 1) % 3;
        let b = (n + 2) % 3;
        // [b_0,b_1,b_2,K_00,...,K_22], in canonical link order.
        let mut sum = [0.0; 12];
        let mut correction = [0.0; 12];
        let mut last = usize::MAX;
        let mut rho = 0.0;
        for (index, link) in self.moving.links.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            let p = coordinates(link.cell, dims);
            if last != link.cell {
                rho = post[link.cell].iter().sum();
                last = link.cell;
            }
            let e = L::E[link.incoming].map(f64::from);
            let g = [
                e[a],
                e[b],
                (p[a] as f64 + 0.5 - self.center[a]) * e[b]
                    - (p[b] as f64 + 0.5 - self.center[b]) * e[a],
            ];
            let weight = 2.0 * L::W[link.incoming] * rho / L::CS2;
            let outgoing = post[link.cell][L::OPP[link.incoming]];
            for i in 0..12 {
                let value = if i < 3 {
                    -2.0 * outgoing * g[i]
                } else {
                    weight * g[(i - 3) / 3] * g[(i - 3) % 3]
                };
                let next = sum[i] + value;
                correction[i] += if sum[i].abs() >= value.abs() {
                    (sum[i] - next) + value
                } else {
                    (value - next) + sum[i]
                };
                sum[i] = next;
            }
        }
        for i in 0..12 {
            sum[i] += correction[i];
        }
        if !sum.iter().all(|x| x.is_finite()) {
            return Err(GridError3::Input("finite wall response overflow"));
        }
        let masses = [self.state.mass, self.state.mass, self.state.axial_inertia];
        let roots = masses.map(fs_math::det::sqrt);
        let matrix: [f64; 9] = core::array::from_fn(|i| {
            let identity = if i / 3 == i % 3 { 1.0 } else { 0.0 };
            identity + 0.5 * (sum[3 + i] / roots[i / 3]) / roots[i % 3]
        });
        let rhs = core::array::from_fn::<_, 3, _>(|i| {
            let resisting = (0..3)
                .map(|j| sum[3 + 3 * i + j] * self.state.velocity[j])
                .sum::<f64>();
            (external[i] + sum[i] - resisting) / roots[i]
        });
        let mut delta = [0.0; 3];
        self.solver
            .solve_into(&matrix, &rhs, &mut delta)
            .map_err(|_| GridError3::Input("finite wall inertia solve"))?;
        let after = core::array::from_fn(|i| self.state.velocity[i] + delta[i] / roots[i]);
        let midpoint = core::array::from_fn(|i| 0.5 * self.state.velocity[i] + 0.5 * after[i]);
        if !after.iter().chain(midpoint.iter()).all(|x| x.is_finite()) {
            return Err(GridError3::Input("nonfinite finite wall response"));
        }
        checkpoint()?;
        Ok((after, midpoint))
    }

    fn finish(
        &self,
        after: [f64; 3],
        midpoint: [f64; 3],
        hydrodynamic: MovingWallMomentum3,
        stationary_walls: StationaryWallMomentum3,
        external: [f64; 3],
        fluid_force_impulse: [f64; 3],
        budget: FiniteWallBudget3,
    ) -> Result<FinitePlanarWallStep3, GridError3<L::Error>> {
        let n = self.face as usize / 2;
        let received = [
            hydrodynamic.wall_impulse[(n + 1) % 3],
            hydrodynamic.wall_impulse[(n + 2) % 3],
            hydrodynamic.wall_angular_impulse[n],
        ];
        let masses = [self.state.mass, self.state.mass, self.state.axial_inertia];
        let changed: [f64; 3] =
            core::array::from_fn(|i| masses[i] * (after[i] - self.state.velocity[i]));
        let residual = core::array::from_fn(|i| changed[i] - received[i] - external[i]);
        for i in 0..3 {
            if !within_relative(
                residual[i],
                changed[i].abs() + received[i].abs() + external[i].abs(),
                budget.maximum_relative_impulse_defect,
            ) {
                return Err(GridError3::Input(
                    "finite wall impulse closure exceeds budget",
                ));
            }
        }
        let kinetic_energy_change = dot(changed, midpoint);
        let external_work = dot(external, midpoint);
        let work_residual = kinetic_energy_change - hydrodynamic.wall_work - external_work;
        if !within_relative(
            work_residual,
            kinetic_energy_change.abs() + hydrodynamic.wall_work.abs() + external_work.abs(),
            budget.maximum_relative_work_defect,
        ) || !fluid_force_impulse.iter().all(|x| x.is_finite())
        {
            return Err(GridError3::Input("finite wall work closure exceeds budget"));
        }
        let mut support_impulse = [0.0; 3];
        support_impulse[n] = -hydrodynamic.wall_impulse[n];
        let mut support_angular_impulse = hydrodynamic.wall_angular_impulse.map(|x| -x);
        support_angular_impulse[n] = 0.0;
        Ok(FinitePlanarWallStep3 {
            before: self.state,
            after: PlanarWallState3 {
                velocity: after,
                ..self.state
            },
            velocity_midpoint: midpoint,
            hydrodynamic,
            stationary_walls,
            support_impulse,
            support_angular_impulse,
            external_impulse: external,
            fluid_force_impulse,
            kinetic_energy_change,
            external_work,
            momentum_residual: residual,
            work_residual,
        })
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn within_relative(defect: f64, scale: f64, tolerance: f64) -> bool {
    defect.is_finite()
        && scale.is_finite()
        && if scale == 0.0 {
            defect == 0.0
        } else {
            (defect / scale).abs() <= tolerance
        }
}

impl<L: Lattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Solve actual drag and an external tangent impulse/axial couple into the
    /// wall's finite inertia, then stream at its solved midpoint velocity.
    /// Fluid and wall velocities publish together only after all admissions.
    /// Other walls remain stationary. Optional collisions retain the original
    /// material/SGS field meaning and act before the response is assembled.
    pub fn step_serial_with_finite_wall(
        &mut self,
        cx: &Cx<'_>,
        wall: &mut FinitePlanarWall3<L, Q>,
        external: [f64; 3],
        budget: FiniteWallBudget3,
        collisions: Option<&[L::Collision]>,
    ) -> Result<FinitePlanarWallStep3, GridError3<L::Error>> {
        self.finite_wall_step(
            wall,
            external,
            budget,
            collisions,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
            |data, src, dst, stream, motion| {
                for (tile, chunk) in dst.chunks_mut(64).enumerate() {
                    cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
                    match motion {
                        Some(motion) => data.process_moving_walls(src, chunk, tile * 64, motion)?,
                        None => data.process(src, chunk, tile * 64, stream)?,
                    }
                }
                Ok(())
            },
        )
    }

    /// Same response and reduction order with the original native tile executor.
    pub fn step_pooled_with_finite_wall<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        wall: &mut FinitePlanarWall3<L, Q>,
        external: [f64; 3],
        budget: FiniteWallBudget3,
        collisions: Option<&[L::Collision]>,
    ) -> Result<FinitePlanarWallStep3, GridError3<L::Error>> {
        self.finite_wall_step(
            wall,
            external,
            budget,
            collisions,
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
            |data, src, dst, stream, motion| {
                run_pass_with_wall_motion(data, src, dst, stream, pool, gate, motion)
            },
        )
    }

    fn finite_wall_step(
        &mut self,
        wall: &mut FinitePlanarWall3<L, Q>,
        external: [f64; 3],
        budget: FiniteWallBudget3,
        collisions: Option<&[L::Collision]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
        mut pass: impl FnMut(
            &StepData<'_, L, Q>,
            &[[f64; Q]],
            &mut [[f64; Q]],
            bool,
            Option<&WallMotion3<'_, L, Q>>,
        ) -> Result<(), GridError3<L::Error>>,
    ) -> Result<FinitePlanarWallStep3, GridError3<L::Error>> {
        checkpoint()?;
        if !external.iter().all(|x| x.is_finite())
            || [
                budget.maximum_relative_impulse_defect,
                budget.maximum_relative_work_defect,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x > 1e-3)
            || wall.solid != self.solid
        {
            return Err(GridError3::Input(
                "finite wall drive, budget or fixed topology",
            ));
        }
        let step = self
            .completed_steps
            .checked_add(1)
            .ok_or(GridError3::Input("step count overflow"))?;
        wall.stationary.validate(self, &mut checkpoint)?;
        wall.fill_velocities(self.dims, wall.state.velocity);
        wall.motion().admit(self, &mut checkpoint)?;
        if let Some(field) = collisions {
            self.admit_collisions(field, &mut checkpoint)?;
        }
        let data = StepData::<L, Q> {
            dims: self.dims,
            periodic: self.periodic,
            solid: &self.solid,
            collision: self.collision,
            collisions,
            force: self.force,
            forces: None,
            budget: self.budget,
            marker: PhantomData,
        };
        pass(&data, &self.f, &mut self.post, false, None)?;
        let (after, midpoint) = wall.solve(&self.post, self.dims, external, &mut checkpoint)?;
        wall.fill_velocities(self.dims, after);
        wall.motion().admit(self, &mut checkpoint)?;
        wall.fill_velocities(self.dims, midpoint);
        let motion = wall.motion();
        let moving = motion.admit(self, &mut checkpoint)?;
        pass(
            &data,
            &self.post,
            &mut self.next,
            true,
            moving.then_some(&motion),
        )?;
        let hydrodynamic = motion.measure(self, &mut checkpoint)?;
        let stationary = wall.stationary.measure(self, &mut checkpoint)?;
        let fluid_cells = self.solid.iter().filter(|&&solid| !solid).count() as f64;
        let receipt = wall.finish(
            after,
            midpoint,
            hydrodynamic,
            stationary,
            external,
            self.force.map(|f| f * fluid_cells),
            budget,
        )?;
        checkpoint()?;
        std::mem::swap(&mut self.f, &mut self.next);
        self.completed_steps = step;
        wall.state = receipt.after;
        Ok(receipt)
    }
}
