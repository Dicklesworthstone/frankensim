//! Spatial, time-dependent body forces on the existing lattice collision/stream owner.
//!
//! Force values are densities in lattice units, not accelerations. At accepted
//! time n, rho*u = sum_q(f_q*e_q) + F_n/2. Collision uses F_n; streamed-state
//! admission and readback use the supplied F_(n+1). Thus a periodic domain's
//! total momentum receives (sum(F_n) + sum(F_(n+1)))/2 per step. This is the
//! trapezoidal source convention, not a switch that adds a spurious half kick.
//! Source evolution is prescribed by the caller; this is not an implicit
//! velocity-dependent force solve, a heat solver, or a total-energy model.

use super::{Grid3, GridError3, Lattice3, StepData, macro_cell, run_pass};
use core::marker::PhantomData;
use fs_exec::{CancelGate, Cx, KernelRunner};

/// Source-preserving material viscosity and density-driven buoyancy.
pub mod material;

/// Native collision capability used before a force-bearing endpoint can be
/// published. This is separate from coefficient validity: D3Q19 and D3Q27
/// families retain their own explicit forcing implementations, including
/// physical-population forcing for the D3Q19 moment families.
pub trait ForcedLattice3<const Q: usize>: Lattice3<Q> {
    fn supports_body_force(collision: Self::Collision) -> bool;
}
impl ForcedLattice3<19> for super::D3q19 {
    fn supports_body_force(collision: Self::Collision) -> bool {
        collision.supports_body_force()
    }
}
impl ForcedLattice3<27> for crate::d3q27::D3q27 {
    fn supports_body_force(_collision: Self::Collision) -> bool {
        true
    }
}

/// A force-aware grid. Its inner grid is deliberately private: returning it
/// alone would erase the force needed to interpret its stored momentum.
/// Source history and populations are published together after both passes.
pub struct ForcedGrid3<L: ForcedLattice3<Q>, const Q: usize> {
    pub(super) grid: Grid3<L, Q>,
    forces: Vec<[f64; 3]>,
}

impl<L: ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Attach an explicit initial force field to an UNSTEPPED grid constructed
    /// with zero uniform force. This specifies the INITIAL population frame;
    /// it does not preserve a previously observed unforced velocity. Supply
    /// initial populations satisfying rho*u = sum(f*e) + F/2. Solid forces must
    /// be zero. Both the original grid admission and this corrected-velocity
    /// admission must pass. No field is attached to a previously evolved state.
    pub fn try_new(
        grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        forces: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        poll()?;
        if grid.completed_steps != 0 || grid.force.iter().any(|&v| v != 0.0) {
            return Err(GridError3::Input(
                "spatial forces require an unstepped zero-force grid",
            ));
        }
        admit_forces(&grid, &forces, None, &mut poll)?;
        for (cell, pop) in grid.f.iter().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            if !grid.solid[cell] {
                macro_cell::<L, Q>(pop, forces[cell], cell, grid.budget.maximum_mach)?;
            }
        }
        poll()?;
        Ok(Self { grid, forces })
    }

    /// Grid dimensions, in x/y/z order.
    pub fn dimensions(&self) -> [usize; 3] {
        self.grid.dimensions()
    }
    /// Accepted populations. Interpret momentum with `Self::macroscopic`.
    pub fn populations(&self) -> &[[f64; Q]] {
        self.grid.populations()
    }
    /// Number of fully accepted steps.
    pub fn completed_steps(&self) -> usize {
        self.grid.completed_steps()
    }
    /// Force densities at the SAME epoch as the accepted populations.
    pub fn forces(&self) -> &[[f64; 3]] {
        &self.forces
    }
    /// Accepted density and force-corrected velocity; no stale uniform force.
    pub fn macroscopic(&self, cell: usize) -> Result<(f64, [f64; 3]), GridError3<L::Error>> {
        let pop = self
            .grid
            .f
            .get(cell)
            .ok_or(GridError3::Input("cell index"))?;
        macro_cell::<L, Q>(pop, self.forces[cell], cell, self.grid.budget.maximum_mach)
    }

    /// Advance with current-epoch collision coefficients and a force field at
    /// the NEXT time. A failure leaves force history, populations and step count
    /// unchanged. The supplied endpoint vector is consumed, never partly kept.
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        serial(&mut self.grid, &self.forces, &forces_at_end, collisions, cx)?;
        self.forces = forces_at_end;
        Ok(())
    }
    /// Hold the accepted spatial force constant for one serial step.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: &[L::Collision],
    ) -> Result<(), GridError3<L::Error>> {
        serial(&mut self.grid, &self.forces, &self.forces, collisions, cx)
    }
    /// Executor counterpart with identical source timing and publication.
    pub fn advance_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(), GridError3<L::Error>> {
        pooled(
            &mut self.grid,
            &self.forces,
            &forces_at_end,
            collisions,
            pool,
            gate,
        )?;
        self.forces = forces_at_end;
        Ok(())
    }
    /// Hold the accepted spatial force constant for one executor step.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: &[L::Collision],
    ) -> Result<(), GridError3<L::Error>> {
        pooled(
            &mut self.grid,
            &self.forces,
            &self.forces,
            collisions,
            pool,
            gate,
        )
    }
}

fn admit_forces<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &Grid3<L, Q>,
    forces: &[[f64; 3]],
    collisions: Option<&[L::Collision]>,
    poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<(), GridError3<L::Error>> {
    poll()?;
    if forces.len() != grid.f.len() {
        return Err(GridError3::Input("force field length"));
    }
    for (cell, force) in forces.iter().enumerate() {
        if cell % 64 == 0 {
            poll()?;
        }
        if force.iter().any(|v| !v.is_finite()) {
            return Err(GridError3::MaterialField {
                cell,
                reason: "nonfinite force density",
            });
        }
        if grid.solid[cell] && force.iter().any(|&v| v != 0.0) {
            return Err(GridError3::MaterialField {
                cell,
                reason: "body force on solid cell",
            });
        }
        let collision = collisions.map_or(grid.collision, |field| field[cell]);
        if force.iter().any(|&v| v != 0.0) && !L::supports_body_force(collision) {
            return Err(GridError3::MaterialField {
                cell,
                reason: "collision model does not admit body force",
            });
        }
    }
    poll()
}

fn serial<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &mut Grid3<L, Q>,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    collisions: &[L::Collision],
    cx: &Cx<'_>,
) -> Result<(), GridError3<L::Error>> {
    let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
    grid.admit_collisions(collisions, &mut poll)?;
    admit_forces(grid, end, Some(collisions), &mut poll)?;
    let step = grid
        .completed_steps
        .checked_add(1)
        .ok_or(GridError3::Input("step count overflow"))?;
    let data = StepData::<L, Q> {
        dims: grid.dims,
        periodic: grid.periodic,
        solid: &grid.solid,
        collision: grid.collision,
        collisions: Some(collisions),
        force: grid.force,
        forces: Some(start),
        budget: grid.budget,
        marker: PhantomData,
    };
    for (tile, chunk) in grid.post.chunks_mut(64).enumerate() {
        poll()?;
        data.process(&grid.f, chunk, tile * 64, false)?;
    }
    let data = StepData {
        forces: Some(end),
        ..data
    };
    for (tile, chunk) in grid.next.chunks_mut(64).enumerate() {
        poll()?;
        data.process(&grid.post, chunk, tile * 64, true)?;
    }
    poll()?;
    core::mem::swap(&mut grid.f, &mut grid.next);
    grid.completed_steps = step;
    Ok(())
}

fn pooled<L: ForcedLattice3<Q>, const Q: usize, P: KernelRunner>(
    grid: &mut Grid3<L, Q>,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    collisions: &[L::Collision],
    pool: &P,
    gate: &CancelGate,
) -> Result<(), GridError3<L::Error>> {
    let mut poll = || {
        if gate.is_requested() {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    };
    grid.admit_collisions(collisions, &mut poll)?;
    admit_forces(grid, end, Some(collisions), &mut poll)?;
    let step = grid
        .completed_steps
        .checked_add(1)
        .ok_or(GridError3::Input("step count overflow"))?;
    let data = StepData::<L, Q> {
        dims: grid.dims,
        periodic: grid.periodic,
        solid: &grid.solid,
        collision: grid.collision,
        collisions: Some(collisions),
        force: grid.force,
        forces: Some(start),
        budget: grid.budget,
        marker: PhantomData,
    };
    run_pass(&data, &grid.f, &mut grid.post, false, pool, gate)?;
    let data = StepData {
        forces: Some(end),
        ..data
    };
    run_pass(&data, &grid.post, &mut grid.next, true, pool, gate)?;
    poll()?;
    core::mem::swap(&mut grid.f, &mut grid.next);
    grid.completed_steps = step;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::CollisionModel3;
    use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
    use crate::lattice3::{D3q19, GridBudget3};
    use fs_exec::{RunError, RunReport, TileKernel, TilePool};
    use std::sync::atomic::{AtomicUsize, Ordering};

    const N: usize = 130;
    const SEED: u64 = 0xF03C_E003;
    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        arenas.scope(|arena| {
            f(&Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: SEED,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ))
        })
    }
    fn bgk() -> CollisionConfig27 {
        CollisionConfig27 {
            model: CollisionModel27::Bgk {
                tau: 0.8,
                equilibrium: Equilibrium27::Quadratic,
            },
            maximum_mach: 0.2,
        }
    }
    fn grid<L: ForcedLattice3<Q>, const Q: usize>(
        collision: L::Collision,
        force: [f64; 3],
    ) -> Grid3<L, Q> {
        Grid3::try_new(
            [N, 1, 1],
            collision,
            force,
            [true; 3],
            GridBudget3 {
                maximum_cells: N,
                maximum_mach: 0.2,
            },
            |_| false,
            |_| L::W,
        )
        .unwrap()
    }
    fn bits<L: ForcedLattice3<Q>, const Q: usize>(pop: &[[f64; Q]]) -> Vec<u64> {
        pop.iter().flatten().map(|v| v.to_bits()).collect()
    }
    fn constant_parity<L: ForcedLattice3<Q>, const Q: usize>(collision: L::Collision) {
        let force = [1e-6, -2e-6, 0.0];
        let mut legacy = grid::<L, Q>(collision, force);
        let mut local = with_cx(|cx| {
            ForcedGrid3::try_new(grid::<L, Q>(collision, [0.0; 3]), cx, vec![force; N]).unwrap()
        });
        let laws = vec![collision; N];
        with_cx(|cx| {
            for _ in 0..5 {
                legacy.step_serial(cx).unwrap();
                local.step_serial(cx, &laws).unwrap();
                assert_eq!(
                    bits::<L, Q>(legacy.populations()),
                    bits::<L, Q>(local.populations())
                );
                assert_eq!(
                    legacy.macroscopic(7).unwrap(),
                    local.macroscopic(7).unwrap()
                );
            }
        });
    }
    #[test]
    fn constant_force_preserves_original_d19_and_d27_paths() {
        constant_parity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
        constant_parity::<D3q19, 19>(CollisionModel3::CentralMoment {
            second_order_rate: 1.25,
            higher_order_rate: 1.0,
        });
        constant_parity::<D3q19, 19>(CollisionModel3::ReducedCumulant {
            second_order_rate: 1.25,
            third_order_rate: 1.35,
            fourth_order_rate: 1.72,
        });
        constant_parity::<D3q27, 27>(bgk());
        let rates = Rates27 {
            omega: [1.25, 1.0, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
        };
        for model in [
            CollisionModel27::CentralMoment { rates },
            CollisionModel27::Cumulant { rates },
        ] {
            constant_parity::<D3q27, 27>(CollisionConfig27 {
                model,
                maximum_mach: 0.2,
            });
        }
    }
    #[test]
    fn force_reversal_has_trapezoidal_impulse_without_an_extra_half_kick() {
        let mut local = with_cx(|cx| {
            ForcedGrid3::try_new(grid::<D3q27, 27>(bgk(), [0.0; 3]), cx, vec![[0.0; 3]; N]).unwrap()
        });
        let laws = vec![bgk(); N];
        let mut previous = 0.0;
        let mut impulse = 0.0;
        for force in [1e-5, 1e-5, -1e-5, -1e-5, 0.0] {
            impulse += 0.5 * (previous + force);
            with_cx(|cx| {
                local
                    .advance_serial(cx, &laws, vec![[force, 0.0, 0.0]; N])
                    .unwrap()
            });
            let (rho, u) = local.macroscopic(0).unwrap();
            assert!((rho * u[0] - impulse).abs() < 2e-14);
            previous = force;
        }
    }
    #[test]
    fn spatial_forces_match_independent_push_streaming() {
        let forces: Vec<_> = (0..N).map(|i| [1e-7 * i as f64, -2e-6, 0.0]).collect();
        let mut local = with_cx(|cx| {
            ForcedGrid3::try_new(grid::<D3q27, 27>(bgk(), [0.0; 3]), cx, forces.clone()).unwrap()
        });
        let mut expected = vec![[0.0; 27]; N];
        for i in 0..N {
            let post = D3q27::collide(local.populations()[i], bgk(), forces[i]).unwrap();
            for q in 0..27 {
                let target = (i as i32 + D3q27::E[q][0]).rem_euclid(N as i32) as usize;
                expected[target][q] = post[q];
            }
        }
        with_cx(|cx| local.step_serial(cx, &vec![bgk(); N]).unwrap());
        assert_eq!(
            bits::<D3q27, 27>(&expected),
            bits::<D3q27, 27>(local.populations())
        );
        let mass: f64 = local.populations().iter().flatten().sum();
        assert!((mass - N as f64).abs() < 1e-10);
    }
    #[test]
    fn dynamic_source_history_is_worker_count_invariant() {
        let make = || {
            with_cx(|cx| {
                ForcedGrid3::try_new(grid::<D3q27, 27>(bgk(), [0.0; 3]), cx, vec![[0.0; 3]; N])
                    .unwrap()
            })
        };
        let field = |step| {
            (0..N)
                .map(|i| [0.0, (i % 5) as f64 * step as f64 * 1e-7, 0.0])
                .collect::<Vec<_>>()
        };
        let laws = vec![bgk(); N];
        let mut serial = make();
        for step in 1..=4 {
            with_cx(|cx| serial.advance_serial(cx, &laws, field(step)).unwrap());
        }
        for workers in [1, 2, 4] {
            let pool = TilePool::for_host(workers, SEED);
            let gate = CancelGate::new_clock_free();
            let mut parallel = make();
            for step in 1..=4 {
                parallel
                    .advance_pooled(&pool, &gate, &laws, field(step))
                    .unwrap();
            }
            assert_eq!(
                bits::<D3q27, 27>(serial.populations()),
                bits::<D3q27, 27>(parallel.populations())
            );
            assert_eq!(serial.forces(), parallel.forces());
        }
    }
    #[test]
    fn cancelled_source_transition_preserves_history_and_retries_exactly() {
        struct CancelPass {
            pool: TilePool,
            calls: AtomicUsize,
            after: usize,
        }
        impl KernelRunner for CancelPass {
            fn workers(&self) -> usize {
                self.pool.workers()
            }
            fn run_with_gate<K: TileKernel>(
                &self,
                kernel: &K,
                gate: &CancelGate,
            ) -> (Result<K::Out, RunError>, RunReport) {
                let result = self.pool.run_with_gate(kernel, gate);
                if self.calls.fetch_add(1, Ordering::AcqRel) == self.after {
                    gate.request();
                }
                result
            }
        }
        let make = || {
            with_cx(|cx| {
                ForcedGrid3::try_new(grid::<D3q27, 27>(bgk(), [0.0; 3]), cx, vec![[0.0; 3]; N])
                    .unwrap()
            })
        };
        let laws = vec![bgk(); N];
        let end = vec![[1e-6, 0.0, 0.0]; N];
        let mut expected = make();
        with_cx(|cx| expected.advance_serial(cx, &laws, end.clone()).unwrap());
        for after in [0, 1] {
            let runner = CancelPass {
                pool: TilePool::for_host(2, SEED),
                calls: AtomicUsize::new(0),
                after,
            };
            let mut actual = make();
            let before = bits::<D3q27, 27>(actual.populations());
            let gate = CancelGate::new_clock_free();
            assert!(matches!(
                actual.advance_pooled(&runner, &gate, &laws, end.clone()),
                Err(GridError3::Cancelled)
            ));
            assert_eq!(bits::<D3q27, 27>(actual.populations()), before);
            assert_eq!(actual.forces(), vec![[0.0; 3]; N]);
            assert_eq!(actual.completed_steps(), 0);
            actual
                .advance_pooled(
                    &runner.pool,
                    &CancelGate::new_clock_free(),
                    &laws,
                    end.clone(),
                )
                .unwrap();
            assert_eq!(
                bits::<D3q27, 27>(actual.populations()),
                bits::<D3q27, 27>(expected.populations())
            );
            assert_eq!(actual.forces(), expected.forces());
        }
    }
    #[test]
    fn malformed_or_late_mach_failure_never_publishes_new_forces() {
        let mut local = with_cx(|cx| {
            ForcedGrid3::try_new(grid::<D3q27, 27>(bgk(), [0.0; 3]), cx, vec![[0.0; 3]; N]).unwrap()
        });
        let laws = vec![bgk(); N];
        let before = bits::<D3q27, 27>(local.populations());
        with_cx(|cx| {
            for n in [0, N - 1, N + 1] {
                assert!(matches!(
                    local.advance_serial(cx, &laws, vec![[0.0; 3]; n]),
                    Err(GridError3::Input("force field length"))
                ));
            }
            for force in [f64::NAN, 1.0] {
                let mut end = vec![[0.0; 3]; N];
                end[N - 1][0] = force;
                assert!(local.advance_serial(cx, &laws, end).is_err());
            }
        });
        assert_eq!(bits::<D3q27, 27>(local.populations()), before);
        assert_eq!(local.forces(), vec![[0.0; 3]; N]);
        assert_eq!(local.completed_steps(), 0);
    }
    #[test]
    fn force_attachment_rejects_advanced_uniform_forced_and_solid_inputs() {
        with_cx(|cx| {
            let mut advanced = grid::<D3q27, 27>(bgk(), [0.0; 3]);
            advanced.step_serial(cx).unwrap();
            assert!(ForcedGrid3::try_new(advanced, cx, vec![[0.0; 3]; N]).is_err());
            assert!(
                ForcedGrid3::try_new(
                    grid::<D3q27, 27>(bgk(), [1e-6, 0.0, 0.0]),
                    cx,
                    vec![[0.0; 3]; N]
                )
                .is_err()
            );
            let mut solid = grid::<D3q27, 27>(bgk(), [0.0; 3]);
            solid.solid[N - 1] = true;
            let mut forces = vec![[0.0; 3]; N];
            forces[N - 1][0] = 1e-6;
            assert!(
                matches!(ForcedGrid3::try_new(solid, cx, forces), Err(GridError3::MaterialField { cell, .. }) if cell == N-1)
            );
        });
    }

    #[test]
    fn d19_moment_models_follow_trapezoidal_force_endpoints() {
        for model in [
            CollisionModel3::CentralMoment {
                second_order_rate: 1.25,
                higher_order_rate: 1.0,
            },
            CollisionModel3::ReducedCumulant {
                second_order_rate: 1.25,
                third_order_rate: 1.35,
                fourth_order_rate: 1.72,
            },
        ] {
            with_cx(|cx| {
                let mut local =
                    ForcedGrid3::try_new(grid::<D3q19, 19>(model, [0.0; 3]), cx, vec![[0.0; 3]; N])
                        .unwrap();
                let mut previous = 0.0;
                let mut impulse = 0.0;
                for force in [1e-5, 2e-5, -1e-5, 0.0] {
                    impulse += 0.5 * (previous + force);
                    local
                        .advance_serial(cx, &vec![model; N], vec![[force, 0.0, 0.0]; N])
                        .unwrap();
                    let (rho, u) = local.macroscopic(0).unwrap();
                    assert!((rho * u[0] - impulse).abs() < 5e-13, "{model:?}");
                    previous = force;
                }
                assert_eq!(local.completed_steps(), 4);
                assert_eq!(local.forces(), vec![[0.0; 3]; N]);
            });
        }
    }
}
