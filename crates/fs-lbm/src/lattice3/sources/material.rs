//! Same-sample material viscosity and buoyancy on the transient 3-D grid.
//!
//! This composes existing MaterialScaling, BoussinesqMaterial and ForcedGrid3.
//! The original reference source and dx/dt/rho0 remain pinned for the run.
//! Material samples may come from an external temperature/composition history;
//! this module advances MOMENTUM, not temperature, chemistry or phase geometry.

use core::fmt;
use fs_exec::{CancelGate, Cx, KernelRunner};
use fs_material::{MaterialError, fluid_transport::FluidTransport};
use crate::lattice3::{Grid3, GridError3, ViscousLattice3};
use crate::material::{MaterialLbmError, MaterialScaling, plan_scaling_from_material};
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use super::{ForcedGrid3, ForcedLattice3};

/// Keep material refusal and kinetic refusal distinct, with the original cell.
#[derive(Debug)]
pub enum MaterialForceError<E> {
    Grid(GridError3<E>),
    Source { cell: usize, source: MaterialError },
    Scaling { cell: usize, source: MaterialLbmError },
    Epoch { expected: usize, supplied: usize },
}
impl<E> From<GridError3<E>> for MaterialForceError<E> {
    fn from(error: GridError3<E>) -> Self { Self::Grid(error) }
}
impl<E: fmt::Debug> fmt::Display for MaterialForceError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-LBM-MATERIAL-FORCE: {self:?}")
    }
}
impl<E: fmt::Debug> std::error::Error for MaterialForceError<E> {}

/// One original material response and the exact applied numerical projection.
/// Source receipts are retained on `material.material()`, not fabricated here.
#[derive(Debug)]
pub struct MaterialForceSample3<F> {
    material: BoussinesqMaterial<F>,
    force_n_m3: [f64; 3],
    tau: f64,
}
impl<F> MaterialForceSample3<F> {
    pub const fn material(&self) -> &BoussinesqMaterial<F> { &self.material }
    pub const fn force_density_n_m3(&self) -> [f64; 3] { self.force_n_m3 }
    pub const fn shear_relaxation_time(&self) -> f64 { self.tau }
}

/// A live momentum owner retaining reference units, source epoch, collision
/// coefficients and force history. The underlying grid has no mutable escape.
///
/// Density departures produce `(rho_source-rho0)*gravity`; inertia keeps rho0.
/// Viscosity uses the SAME source response. Cp/k are retained and admitted by
/// the shared transport owners, but this route does not apply a heat equation.
/// Sources, pressure/composition validity and any prescribed temperature
/// history remain caller-owned. No separate beta*g load is added.
pub struct MaterialForceGrid3<'a, L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize, R: FluidTransport + ?Sized, F> {
    reference: MaterialScaling<'a, R>,
    config: BoussinesqConfig<3>,
    flow: ForcedGrid3<L, Q>,
    samples: Vec<Option<MaterialForceSample3<F>>>,
    collisions: Vec<L::Collision>,
    sampled_at_step: usize,
}

struct Prepared<L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize, F> {
    samples: Vec<Option<MaterialForceSample3<F>>>,
    collisions: Vec<L::Collision>,
    forces: Vec<[f64; 3]>,
}

impl<'a, L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize, R: FluidTransport + ?Sized, F: FluidTransport>
    MaterialForceGrid3<'a, L, Q, R, F>
{
    /// Bind an unstepped, zero-uniform-force grid to its original scaling and
    /// complete initial source field. Each fluid cell requires Some(source),
    /// each solid None. Initial populations follow ForcedGrid3's documented
    /// force-shifted convention. The scenario speed in `reference` also bounds
    /// runtime velocity; the grid's stricter Mach limit is never relaxed.
    pub fn try_new(
        mut grid: Grid3<L, Q>,
        cx: &Cx<'_>,
        reference: MaterialScaling<'a, R>,
        config: BoussinesqConfig<3>,
        sources: Vec<Option<F>>,
    ) -> Result<Self, MaterialForceError<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        if grid.completed_steps != 0 || grid.force.iter().any(|&v| v != 0.0) {
            return Err(GridError3::Input("material forcing needs an unstepped zero-force grid").into());
        }
        // Tighten admission only; this does not rescale any population or rate.
        grid.budget.maximum_mach = grid.budget.maximum_mach.min(reference.flow().u_lattice / L::CS2.sqrt());
        let prepared = prepare(&grid, &reference, config, sources, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        let flow = ForcedGrid3::try_new(grid, cx, prepared.forces)?;
        Ok(Self {
            reference, config, flow, samples: prepared.samples,
            collisions: prepared.collisions, sampled_at_step: 0,
        })
    }
    /// Original borrowed source, retained SI mapping and reference transport.
    pub const fn reference(&self) -> &MaterialScaling<'a, R> { &self.reference }
    /// Read-only accepted kinetic state, including its same-epoch force field.
    pub const fn flow(&self) -> &ForcedGrid3<L, Q> { &self.flow }
    /// Original responses used to prepare the current coefficients.
    pub fn samples(&self) -> &[Option<MaterialForceSample3<F>>] { &self.samples }
    pub fn collision_field(&self) -> &[L::Collision] { &self.collisions }
    /// Last explicitly supplied material endpoint; frozen steps do not pretend
    /// to re-evaluate the sources at a later time.
    pub const fn sampled_at_step(&self) -> usize { self.sampled_at_step }
    /// Numerical elapsed SI time from the pinned dt; every step refuses time
    /// overflow or a floating timestamp that can no longer advance.
    pub fn time_s(&self) -> f64 { self.flow.completed_steps() as f64 * self.reference.timestep_s() }

    /// Advance one step using current material viscosity and endpoint buoyancy.
    /// New viscosity becomes current only after success, so it drives the NEXT
    /// collision. Force impulse is trapezoidal; coefficient evolution is an
    /// explicit lagged update. Refine dt; this is not an implicit stiff-feedback
    /// solver. Returned failures roll back populations, forces AND source epoch.
    /// Sources must not mutate shared external simulation history when queried.
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        endpoint_step: usize,
        sources_at_end: Vec<Option<F>>,
    ) -> Result<(), MaterialForceError<L::Error>> {
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        self.admit_epoch(endpoint_step)?;
        let pending = prepare(&self.flow.grid, &self.reference, self.config, sources_at_end, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.flow.advance_serial(cx, &self.collisions, pending.forces)?;
        self.samples = pending.samples;
        self.collisions = pending.collisions;
        self.sampled_at_step = endpoint_step;
        Ok(())
    }
    /// Executor-backed endpoint update with the same atomic source timing.
    pub fn advance_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        endpoint_step: usize,
        sources_at_end: Vec<Option<F>>,
    ) -> Result<(), MaterialForceError<L::Error>> {
        if gate.is_requested() { return Err(GridError3::Cancelled.into()); }
        self.admit_epoch(endpoint_step)?;
        let pending = prepare(&self.flow.grid, &self.reference, self.config, sources_at_end, || {
            if gate.is_requested() { Err(GridError3::Cancelled) } else { Ok(()) }
        })?;
        self.flow.advance_pooled(pool, gate, &self.collisions, pending.forces)?;
        self.samples = pending.samples;
        self.collisions = pending.collisions;
        self.sampled_at_step = endpoint_step;
        Ok(())
    }
    /// Advance with the explicitly frozen material snapshot. No new source
    /// epoch is claimed. Caller must establish continued validity of that hold.
    pub fn step_serial(&mut self, cx: &Cx<'_>) -> Result<(), MaterialForceError<L::Error>> {
        self.admit_epoch(self.next_step()?)?;
        self.flow.step_serial(cx, &self.collisions)?;
        Ok(())
    }
    pub fn step_pooled<P: KernelRunner>(&mut self, pool: &P, gate: &CancelGate) -> Result<(), MaterialForceError<L::Error>> {
        self.admit_epoch(self.next_step()?)?;
        self.flow.step_pooled(pool, gate, &self.collisions)?;
        Ok(())
    }
    fn next_step(&self) -> Result<usize, MaterialForceError<L::Error>> {
        self.flow.completed_steps().checked_add(1).ok_or_else(|| GridError3::Input("step count overflow").into())
    }
    fn admit_epoch(&self, supplied: usize) -> Result<(), MaterialForceError<L::Error>> {
        let expected = self.next_step()?;
        if expected != supplied { return Err(MaterialForceError::Epoch { expected, supplied }); }
        let time = supplied as f64 * self.reference.timestep_s();
        if !time.is_finite() || time <= self.time_s() {
            return Err(GridError3::Input("material endpoint time is not representable").into());
        }
        Ok(())
    }
}

fn prepare<L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize, R: FluidTransport + ?Sized, F: FluidTransport>(
    grid: &Grid3<L, Q>,
    reference: &MaterialScaling<'_, R>,
    config: BoussinesqConfig<3>,
    sources: Vec<Option<F>>,
    mut poll: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<Prepared<L, Q, F>, MaterialForceError<L::Error>> {
    poll()?;
    config.validate().map_err(|_| GridError3::Input("invalid material buoyancy configuration"))?;
    if sources.len() != grid.f.len() { return Err(GridError3::Input("material source field length").into()); }
    // Check ownership before invoking any producer.
    for (cell, source) in sources.iter().enumerate() {
        if cell % 64 == 0 { poll()?; }
        if source.is_some() == grid.solid[cell] {
            return Err(GridError3::MaterialField { cell, reason: "material source/solid ownership mismatch" }.into());
        }
    }
    let mut samples = Vec::new();
    let mut collisions = Vec::new();
    let mut forces = Vec::new();
    samples.try_reserve_exact(sources.len()).map_err(|_| GridError3::Allocation)?;
    collisions.try_reserve_exact(sources.len()).map_err(|_| GridError3::Allocation)?;
    forces.try_reserve_exact(sources.len()).map_err(|_| GridError3::Allocation)?;
    let max_speed = reference.velocity_si(reference.flow().u_lattice)
        .map_err(|source| MaterialForceError::Scaling { cell: 0, source })?;
    for (cell, source) in sources.into_iter().enumerate() {
        poll()?;
        let Some(source) = source else {
            samples.push(None); collisions.push(grid.collision); forces.push([0.0; 3]);
            continue;
        };
        // Freeze the original source exactly once. Post-producer cancellation
        // is checked before interpreting either a value or a source refusal.
        let response = BoussinesqMaterial::try_new(source, reference.transport(), config);
        poll()?;
        let material = response.map_err(|source| MaterialForceError::Source { cell, source })?;
        let local = plan_scaling_from_material(&material, reference.cell_size_m(), reference.timestep_s(), max_speed)
            .map_err(|source| MaterialForceError::Scaling { cell, source })?;
        let tau = 0.5 + local.flow().viscosity / L::CS2;
        let collision = L::with_shear_relaxation(grid.collision, tau)
            .map_err(|source| GridError3::Collision { cell, source })?;
        let departure = material.resolved_transport().density_kg_m3() - reference.transport().density_kg_m3();
        let mut force_n_m3 = [0.0; 3];
        let mut force = [0.0; 3];
        for axis in 0..3 {
            let gravity = config.gravity_m_s2[axis];
            force_n_m3[axis] = departure * gravity;
            if departure != 0.0 && gravity != 0.0 && force_n_m3[axis] == 0.0 {
                return Err(GridError3::MaterialField { cell, reason: "buoyancy force underflow" }.into());
            }
            force[axis] = reference.force_density_lattice(force_n_m3[axis])
                .map_err(|source| MaterialForceError::Scaling { cell, source })?;
        }
        samples.push(Some(MaterialForceSample3 { material, force_n_m3, tau }));
        collisions.push(collision);
        forces.push(force);
    }
    poll()?;
    Ok(Prepared { samples, collisions, forces })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::CollisionModel3;
    use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
    use crate::lattice3::{D3q19, GridBudget3, Lattice3};
    use fs_material::fluid_transport::FluidTransportState;
    use fs_material::gas::{GasSpec, GasState};
    use fs_exec::{RunError, RunReport, TileKernel, TilePool};
    use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};

    const SEED: u64 = 0xB005_003;
    fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        arenas.scope(|arena| f(&Cx::new(gate, arena,
            fs_exec::StreamKey { seed: SEED, kernel_id: 3, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic)))
    }
    fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T { with_gate(&CancelGate::new_clock_free(), f) }
    fn config(gravity: [f64; 3]) -> BoussinesqConfig<3> {
        BoussinesqConfig { gravity_m_s2: gravity, max_relative_density_departure: 0.1, max_relative_specific_heat_departure: 0.0 }
    }
    fn fluid(rho: f64, mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(rho, 4000.0, 0.5, mu).unwrap()
    }
    fn bgk() -> CollisionConfig27 {
        CollisionConfig27 { model: CollisionModel27::Bgk { tau: 0.8, equilibrium: Equilibrium27::Quadratic }, maximum_mach: 0.2 }
    }
    fn grid<L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(n: usize, collision: L::Collision) -> Grid3<L, Q> {
        Grid3::try_new([n, 1, 1], collision, [0.0; 3], [true; 3],
            GridBudget3 { maximum_cells: n, maximum_mach: 0.2 }, |_| false, |_| L::W).unwrap()
    }
    fn bits<const Q: usize>(pop: &[[f64; Q]]) -> Vec<u64> {
        pop.iter().flatten().map(|v| v.to_bits()).collect()
    }

    #[test]
    fn gas_temperature_changes_viscosity_and_buoyancy_from_one_state() {
        let gas = |t| GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 100000.0).unwrap();
        let reference = gas(300.0);
        let scale = plan_scaling_from_material(&reference, 0.001, 0.001, 0.1).unwrap();
        let cfg = BoussinesqConfig { max_relative_specific_heat_departure: 0.1, ..config([0.0, -9.81, 0.0]) };
        let mut run = with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(8, bgk()), cx, scale, cfg, (0..8).map(|_| Some(gas(300.0))).collect()).unwrap());
        let initial_tau = run.samples()[0].as_ref().unwrap().shear_relaxation_time();
        with_cx(|cx| run.advance_serial(cx, 1, (0..8).map(|_| Some(gas(310.0))).collect()).unwrap());
        let sample = run.samples()[0].as_ref().unwrap();
        assert!(sample.material().resolved_transport().density_kg_m3() < run.reference().transport().density_kg_m3());
        assert_eq!(sample.material().reference_storage_transport().density_kg_m3(), run.reference().transport().density_kg_m3());
        assert!(sample.shear_relaxation_time() > initial_tau);
        assert!(sample.force_density_n_m3()[1] > 0.0);
        assert_eq!(sample.material().material().temperature, 310.0);
        assert!(core::ptr::eq(run.reference().material(), &reference));
        let expected_force = run.reference().force_density_lattice(sample.force_density_n_m3()[1]).unwrap();
        let (rho, velocity) = run.flow().macroscopic(0).unwrap();
        assert!((rho * velocity[1] - expected_force / 2.0).abs() < 2e-14);
        with_cx(|cx| run.step_serial(cx).unwrap());
        let (rho, velocity) = run.flow().macroscopic(0).unwrap();
        assert!((rho * velocity[1] - 1.5 * expected_force).abs() < 2e-14);
        assert_eq!(run.sampled_at_step(), 1); // holding is not re-sampling
        assert_eq!(run.flow().completed_steps(), 2);
        assert_eq!(run.time_s(), 0.002);
    }

    fn sinusoidal_buoyancy<L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(n: usize, model: L::Collision, mu: f64) -> f64 {
        let reference = fluid(1000.0, mu);
        let dx = 0.016 / n as f64;
        let dt = 0.01 * (16.0 / n as f64).powi(2);
        let scale = plan_scaling_from_material(&reference, dx, dt, 0.01).unwrap();
        let source = (0..n).map(|i| Some(fluid(
            1000.0 - 10.0 * fs_math::det::sin(core::f64::consts::TAU * i as f64 / n as f64), mu))).collect();
        let mut run = with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<L, Q>(n, model), cx, scale, config([0.0, -0.001, 0.0]), source).unwrap());
        with_cx(|cx| { for _ in 0..12*n*n { run.step_serial(cx).unwrap(); } });
        let k = core::f64::consts::TAU / 0.016;
        // Independent stationary Stokes solution: nu*u'' + b*sin(k*x)=0.
        let amplitude = (10.0 / 1000.0 * 0.001) / ((mu / 1000.0) * k * k);
        let mut error = 0.0_f64;
        for i in 0..n {
            let (rho, u) = run.flow().macroscopic(i).unwrap();
            assert!((rho - 1.0).abs() < 1e-9);
            let observed = run.reference().velocity_si(u[1]).unwrap();
            let exact = amplitude * fs_math::det::sin(core::f64::consts::TAU * i as f64 / n as f64);
            error = error.max((observed - exact).abs() / amplitude);
        }
        error
    }
    #[test]
    fn material_buoyancy_flow_refines_against_independent_stokes_solution() {
        for (coarse, fine) in [
            (sinusoidal_buoyancy::<D3q19, 19>(16, CollisionModel3::Bgk { tau: 0.8 }, 0.01),
             sinusoidal_buoyancy::<D3q19, 19>(32, CollisionModel3::Bgk { tau: 0.8 }, 0.01)),
            (sinusoidal_buoyancy::<D3q27, 27>(16, bgk(), 0.01), sinusoidal_buoyancy::<D3q27, 27>(32, bgk(), 0.01)),
            // Higher tau needs finer spatial resolution for the SAME error bounds.
            (sinusoidal_buoyancy::<D3q27, 27>(32, bgk(), 0.04), sinusoidal_buoyancy::<D3q27, 27>(64, bgk(), 0.04)),
        ] {
            assert!(coarse < 0.08, "coarse error {coarse}");
            assert!(fine < 0.025, "fine error {fine}");
            let order = fs_math::det::ln(coarse / fine) / fs_math::det::ln(2.0);
            assert!((1.8..=2.2).contains(&order), "order {order}: {coarse} -> {fine}");
        }
    }
    #[test]
    fn endpoint_viscosity_drives_the_next_collision_not_the_past_one() {
        let reference = fluid(1000.0, 0.01);
        let make = || with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(8, bgk()), cx,
            plan_scaling_from_material(&reference, 0.001, 0.01, 0.01).unwrap(),
            config([0.0, -0.001, 0.0]),
            (0..8).map(|i| Some(fluid(if i < 4 { 990.0 } else { 1000.0 }, 0.01))).collect()).unwrap());
        let mut changed = make();
        let mut fixed = make();
        let sources = |mu| (0..8).map(|i| Some(fluid(if i < 4 { 990.0 } else { 1000.0 }, mu))).collect();
        with_cx(|cx| {
            changed.advance_serial(cx, 1, sources(0.04)).unwrap();
            fixed.advance_serial(cx, 1, sources(0.01)).unwrap();
            assert_eq!(bits(changed.flow().populations()), bits(fixed.flow().populations()));
            assert_ne!(changed.collision_field(), fixed.collision_field());
            for _ in 0..4 { changed.step_serial(cx).unwrap(); fixed.step_serial(cx).unwrap(); }
        });
        assert_ne!(bits(changed.flow().populations()), bits(fixed.flow().populations()));
        assert_eq!(changed.flow().forces(), fixed.flow().forces());
    }
    #[test]
    fn stale_epoch_and_source_domain_failures_preserve_all_accepted_state() {
        let reference = fluid(1000.0, 0.01);
        let mut run = with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(8, bgk()), cx,
            plan_scaling_from_material(&reference, 0.001, 0.01, 0.01).unwrap(),
            config([0.0, -0.001, 0.0]), vec![Some(reference); 8]).unwrap());
        let before = bits(run.flow().populations());
        let laws = run.collision_field().to_vec();
        with_cx(|cx| {
            assert!(matches!(run.advance_serial(cx, 0, vec![Some(reference); 8]), Err(MaterialForceError::Epoch { expected: 1, supplied: 0 })));
            for bad in [fluid(800.0, 0.01), FluidTransportState::try_new(1000.0, 4100.0, 0.5, 0.01).unwrap()] {
                let mut next = vec![Some(reference); 8]; next[7] = Some(bad);
                assert!(matches!(run.advance_serial(cx, 1, next), Err(MaterialForceError::Source { cell: 7, .. })));
            }
            let mut missing = vec![Some(reference); 8]; missing[7] = None;
            assert!(run.advance_serial(cx, 1, missing).is_err());
            assert!(run.advance_serial(cx, 1, vec![Some(reference); 7]).is_err());
        });
        assert_eq!(bits(run.flow().populations()), before);
        assert_eq!(run.collision_field(), laws);
        assert_eq!(run.flow().forces(), vec![[0.0; 3]; 8]);
        assert_eq!(run.sampled_at_step(), 0);
        assert_eq!(run.flow().completed_steps(), 0);
    }
    #[test]
    fn source_producer_is_queried_once_and_its_cancellation_rolls_back() {
        struct Source { state: FluidTransportState, calls: Arc<AtomicUsize>, gate: Arc<CancelGate>, trip: bool }
        impl FluidTransport for Source {
            fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
                self.calls.fetch_add(1, Ordering::SeqCst);
                if self.trip { self.gate.request(); }
                Ok(self.state)
            }
        }
        let gate = Arc::new(CancelGate::new_clock_free());
        let calls = Arc::new(AtomicUsize::new(0));
        let reference = fluid(1000.0, 0.01);
        let sources = |trip| (0..8).map(|i| Some(Source { state: reference, calls: calls.clone(), gate: gate.clone(), trip: trip && i == 7 })).collect();
        let mut run = with_gate(&gate, |cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(8, bgk()), cx,
            plan_scaling_from_material(&reference, 0.001, 0.01, 0.01).unwrap(),
            config([0.0, -0.001, 0.0]), sources(false)).unwrap());
        assert_eq!(calls.load(Ordering::SeqCst), 8);
        let before = bits(run.flow().populations());
        with_gate(&gate, |cx| {
            assert!(matches!(run.advance_serial(cx, 1, sources(true)), Err(MaterialForceError::Grid(GridError3::Cancelled))));
        });
        assert_eq!(calls.load(Ordering::SeqCst), 16);
        assert_eq!(bits(run.flow().populations()), before);
        assert_eq!(run.flow().completed_steps(), 0);
        assert_eq!(run.sampled_at_step(), 0);
        assert!(run.samples().iter().all(|s| !s.as_ref().unwrap().material().material().trip));
    }
    #[test]
    fn material_source_publication_follows_pooled_collision_and_stream_success() {
        struct CancelAfter { pool: TilePool, calls: AtomicUsize }
        impl KernelRunner for CancelAfter {
            fn workers(&self) -> usize { self.pool.workers() }
            fn run_with_gate<K: TileKernel>(&self, k: &K, gate: &CancelGate) -> (Result<K::Out, RunError>, RunReport) {
                let result = self.pool.run_with_gate(k, gate);
                if self.calls.fetch_add(1, Ordering::SeqCst) == 1 { gate.request(); }
                result
            }
        }
        let reference = fluid(1000.0, 0.01);
        let make = || with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(130, bgk()), cx,
            plan_scaling_from_material(&reference, 0.001, 0.01, 0.01).unwrap(),
            config([0.0, -0.001, 0.0]), vec![Some(reference); 130]).unwrap());
        let next = || (0..130).map(|i| Some(fluid(1000.0 - (i % 5) as f64, 0.02))).collect();
        let mut run = make();
        let before = bits(run.flow().populations());
        let runner = CancelAfter { pool: TilePool::for_host(2, SEED), calls: AtomicUsize::new(0) };
        assert!(matches!(run.advance_pooled(&runner, &CancelGate::new_clock_free(), 1, next()), Err(MaterialForceError::Grid(GridError3::Cancelled))));
        assert_eq!(bits(run.flow().populations()), before);
        assert_eq!(run.sampled_at_step(), 0);
        assert_eq!(run.flow().completed_steps(), 0);
        assert!(run.samples().iter().all(|s| s.as_ref().unwrap().material().resolved_transport() == &reference));
        run.advance_pooled(&runner.pool, &CancelGate::new_clock_free(), 1, next()).unwrap();
        let mut serial = make();
        with_cx(|cx| serial.advance_serial(cx, 1, next()).unwrap());
        assert_eq!(bits(run.flow().populations()), bits(serial.flow().populations()));
        assert_eq!(run.flow().forces(), serial.flow().forces());
        assert_eq!(run.collision_field(), serial.collision_field());
    }
    #[test]
    fn retained_physical_speed_bound_refuses_a_material_force_endpoint() {
        let reference = fluid(1000.0, 0.01);
        let mut run = with_cx(|cx| MaterialForceGrid3::try_new(
            grid::<D3q27, 27>(8, bgk()), cx,
            plan_scaling_from_material(&reference, 0.001, 0.01, 1e-8).unwrap(),
            config([0.0, -1.0, 0.0]), vec![Some(reference); 8]).unwrap());
        let before = bits(run.flow().populations());
        with_cx(|cx| {
            assert!(matches!(run.advance_serial(cx, 1, vec![Some(fluid(990.0, 0.01)); 8]), Err(MaterialForceError::Grid(GridError3::Mach { .. }))));
        });
        assert_eq!(bits(run.flow().populations()), before);
        assert_eq!(run.sampled_at_step(), 0);
        assert_eq!(run.flow().completed_steps(), 0);
    }
}
