//! Generalized-Newtonian material viscosity on the existing 3-D flow owner.
//!
//! Each incoming velocity field supplies a wall-aware finite-difference strain
//! rate. The original material producer resolves mu at that physical rate; the
//! same response drives this step's shear relaxation. This is an explicit time
//! discretization, not a local nonequilibrium-stress fixed point, an implicit
//! flow solve, hereditary rheology or a heat/temperature integrator.

use super::sources::{ForcedGrid3, ForcedLattice3};
use super::{Grid3, GridError3, ViscousLattice3, coordinates, macro_cell};
use crate::material::{MaterialLbmError, MaterialScaling, plan_scaling_from_material};
use fs_exec::{CancelGate, Cx, KernelRunner};
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

/// Conditions at the incoming kinetic time level, with x-fastest cell indexing.
/// Temperature, composition and material validity remain with the producer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RheologyPoint3 {
    pub cell: usize,
    pub coordinates: [usize; 3],
    pub time_s: f64,
    /// sqrt(2 D:D) [1/s], D=sym(grad u), including every diagonal/shear term.
    pub shear_rate_per_s: f64,
}

/// The original source response and the exact numerical snapshot applied.
#[derive(Debug)]
pub struct RheologySample3<F> {
    point: RheologyPoint3,
    source: F,
    transport: FluidTransportState,
    tau: f64,
}

impl<F> RheologySample3<F> {
    pub const fn point(&self) -> RheologyPoint3 {
        self.point
    }
    pub const fn source(&self) -> &F {
        &self.source
    }
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    pub const fn shear_relaxation_time(&self) -> f64 {
        self.tau
    }
}

/// Returned only after the complete collision/stream step succeeds.
/// Samples and collision laws refer to the INCOMING time level; solid entries
/// have no material sample. Reusing these coefficients does not refresh rheology.
#[derive(Debug)]
pub struct RheologyStep3<C, F> {
    pub sampled_at_step: usize,
    pub samples: Vec<Option<RheologySample3<F>>>,
    pub collisions: Vec<C>,
}

#[derive(Debug)]
pub enum RheologyError3<E> {
    Grid(GridError3<E>),
    Source {
        cell: usize,
        source: MaterialError,
    },
    Scaling {
        cell: usize,
        source: MaterialLbmError,
    },
}

impl<E> From<GridError3<E>> for RheologyError3<E> {
    fn from(error: GridError3<E>) -> Self {
        Self::Grid(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RheologyError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-RHEOLOGY3: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RheologyError3<E> {}

// sqrt(2 D:D) in lattice units. Rigid rotation cancels before the norm; scaling
// the six independent entries avoids losing a representable rate by squaring it.
fn strain_rate(gradient: [[f64; 3]; 3]) -> f64 {
    let components = [
        core::f64::consts::SQRT_2 * gradient[0][0],
        core::f64::consts::SQRT_2 * gradient[1][1],
        core::f64::consts::SQRT_2 * gradient[2][2],
        gradient[0][1] + gradient[1][0],
        gradient[0][2] + gradient[2][0],
        gradient[1][2] + gradient[2][1],
    ];
    let scale = components.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    scale
        * components
            .iter()
            .map(|v| {
                let x = v / scale;
                x * x
            })
            .sum::<f64>()
            .sqrt()
}

fn prepare<L: ViscousLattice3<Q>, const Q: usize, R: FluidTransport + ?Sized, F: FluidTransport>(
    grid: &Grid3<L, Q>,
    reference: &MaterialScaling<'_, R>,
    forces: Option<&[[f64; 3]]>,
    mut resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    mut poll: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
    poll()?;
    let cells = grid.f.len();
    if forces.is_some_and(|f| f.len() != cells) {
        return Err(GridError3::Input("rheology force field length").into());
    }
    let dt = reference.timestep_s();
    let time = grid.completed_steps as f64 * dt;
    if !time.is_finite() {
        return Err(GridError3::Input("rheology time is unrepresentable").into());
    }
    let speed = reference
        .velocity_si(reference.flow().u_lattice)
        .map_err(|source| RheologyError3::Scaling { cell: 0, source })?;
    let mut velocity = Vec::new();
    velocity
        .try_reserve_exact(cells)
        .map_err(|_| GridError3::Allocation)?;
    for cell in 0..cells {
        if cell % 64 == 0 {
            poll()?;
        }
        let u = if grid.solid[cell] {
            [0.0; 3]
        } else {
            macro_cell::<L, Q>(
                &grid.f[cell],
                forces.map_or(grid.force, |f| f[cell]),
                cell,
                grid.budget.maximum_mach,
            )?
            .1
        };
        velocity.push(u);
    }
    let mut samples = Vec::new();
    let mut collisions = Vec::new();
    samples
        .try_reserve_exact(cells)
        .map_err(|_| GridError3::Allocation)?;
    collisions
        .try_reserve_exact(cells)
        .map_err(|_| GridError3::Allocation)?;
    for cell in 0..cells {
        poll()?;
        if grid.solid[cell] {
            samples.push(None);
            collisions.push(grid.collision);
            continue;
        }
        let lattice_rate = strain_rate(grid.velocity_gradient(cell, &velocity));
        let rate = lattice_rate / dt;
        if !rate.is_finite() || (lattice_rate != 0.0 && rate == 0.0) {
            return Err(GridError3::MaterialField {
                cell,
                reason: "physical strain rate is unrepresentable",
            }
            .into());
        }
        let point = RheologyPoint3 {
            cell,
            coordinates: coordinates(cell, grid.dims),
            time_s: time,
            shear_rate_per_s: rate,
        };
        let response = resolve(point, reference.transport());
        poll()?;
        let source = response.map_err(|source| RheologyError3::Source { cell, source })?;
        let response = source.transport_state();
        poll()?;
        let transport = response.map_err(|source| RheologyError3::Source { cell, source })?;
        if transport.density_kg_m3().to_bits() != reference.transport().density_kg_m3().to_bits() {
            return Err(GridError3::MaterialField {
                cell,
                reason: "rheology reference density differs",
            }
            .into());
        }
        // Reuse the original scaling owner with the frozen transport snapshot:
        // the producer is queried exactly once, and no source receipt is invented.
        let local = plan_scaling_from_material(&transport, reference.cell_size_m(), dt, speed)
            .map_err(|source| RheologyError3::Scaling { cell, source })?;
        let tau = 0.5 + local.flow().viscosity / L::CS2;
        let collision = L::with_shear_relaxation(grid.collision, tau)
            .map_err(|source| GridError3::Collision { cell, source })?;
        collisions.push(collision);
        samples.push(Some(RheologySample3 {
            point,
            source,
            transport,
            tau,
        }));
    }
    poll()?;
    Ok(RheologyStep3 {
        sampled_at_step: grid.completed_steps,
        samples,
        collisions,
    })
}

fn gate_poll<E>(gate: &CancelGate) -> Result<(), GridError3<E>> {
    if gate.is_requested() {
        Err(GridError3::Cancelled)
    } else {
        Ok(())
    }
}

impl<L: ViscousLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Query the original material at this incoming velocity's physical strain
    /// rate, then advance the real grid once with those molecular viscosities.
    ///
    /// Keep the original `reference` dx/dt/rho0 throughout the run. Sources must
    /// retain that density; only mu changes momentum transport. Cp/k remain in
    /// the original response but this path advances neither heat nor buoyancy.
    /// The callback may use an external temperature/composition history, and
    /// must not advance external history during these tentative queries.
    ///
    /// Gradients use physical half-force velocities, centered interior/periodic
    /// differences and the existing quadratic stationary-halfway-wall stencil.
    /// This is explicit incoming-time rheology, not an implicit constitutive or
    /// flow solve. Independent bulk/higher rates and equilibrium choices are
    /// preserved; D3Q19 moment models retain their existing bulk/shear tie.
    /// No viscosity floor, cap, fit, source-domain extension or warm-start is
    /// introduced. A source/collision/cancellation refusal publishes nothing.
    pub fn step_serial_with_material_rheology<R: FluidTransport + ?Sized, F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(self, reference, None, resolve, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.step_serial_with_collisions(cx, &report.collisions)?;
        Ok(report)
    }

    /// Same material samples and collision laws, with existing pooled stepping.
    /// Source queries remain serial and deterministic; no producer is cloned or
    /// called from workers. The returned receipt is published only on success.
    pub fn step_pooled_with_material_rheology<
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(self, reference, None, resolve, || gate_poll(gate))?;
        self.step_pooled_with_collisions(pool, gate, &report.collisions)?;
        Ok(report)
    }
}

impl<L: ViscousLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Incoming-time rheology using this owner's spatial force field in every
    /// macroscopic velocity. The existing force history is held for one step.
    /// See `Grid3::step_serial_with_material_rheology` for material/time semantics.
    pub fn step_serial_with_material_rheology<R: FluidTransport + ?Sized, F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(&self.grid, reference, Some(self.forces()), resolve, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.step_serial(cx, &report.collisions)?;
        Ok(report)
    }

    /// Use current physical strain for viscosity, and supply force at the NEXT
    /// time to the existing trapezoidal force owner. A new force endpoint must
    /// never be used to reinterpret the incoming velocity or material query.
    pub fn advance_serial_with_material_rheology<R: FluidTransport + ?Sized, F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        forces_at_end: Vec<[f64; 3]>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(&self.grid, reference, Some(self.forces()), resolve, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })?;
        self.advance_serial(cx, &report.collisions, forces_at_end)?;
        Ok(report)
    }

    /// Pooled counterpart retaining the current spatial force field.
    pub fn step_pooled_with_material_rheology<
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(&self.grid, reference, Some(self.forces()), resolve, || {
            gate_poll(gate)
        })?;
        self.step_pooled(pool, gate, &report.collisions)?;
        Ok(report)
    }

    /// Pooled counterpart with an explicit next-time force endpoint.
    pub fn advance_pooled_with_material_rheology<
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        forces_at_end: Vec<[f64; 3]>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<RheologyStep3<L::Collision, F>, RheologyError3<L::Error>> {
        let report = prepare(&self.grid, reference, Some(self.forces()), resolve, || {
            gate_poll(gate)
        })?;
        self.advance_pooled(pool, gate, &report.collisions, forces_at_end)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::CollisionModel3;
    use crate::d3q27::{
        CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
    };
    use crate::lattice3::{D3q19, GridBudget3, Lattice3};
    use std::{cell::Cell, rc::Rc};

    fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            f(&Cx::new(
                gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0x3D_4E0,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ))
        })
    }
    fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
        with_gate(&CancelGate::new_clock_free(), f)
    }
    fn fluid(mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(1000.0, 4000.0, 0.5, mu).unwrap()
    }
    fn collision(model: usize) -> CollisionConfig27 {
        let rates = Rates27 {
            omega: [1.0 / 0.9, 1.1, 0.8, 1.2, 0.9, 1.3, 0.7, 1.0, 1.4, 0.6],
        };
        CollisionConfig27 {
            model: match model {
                0 => CollisionModel27::Bgk {
                    tau: 0.9,
                    equilibrium: Equilibrium27::Quadratic,
                },
                1 => CollisionModel27::Bgk {
                    tau: 0.9,
                    equilibrium: Equilibrium27::Factorized,
                },
                2 => CollisionModel27::CentralMoment { rates },
                _ => CollisionModel27::Cumulant { rates },
            },
            maximum_mach: 0.2,
        }
    }
    fn grid27(
        dims: [usize; 3],
        model: usize,
        force: [f64; 3],
        periodic: [bool; 3],
        initial: impl Fn([usize; 3]) -> [f64; 27],
    ) -> Grid3<D3q27, 27> {
        Grid3::try_new(
            dims,
            collision(model),
            force,
            periodic,
            GridBudget3 {
                maximum_cells: dims.iter().product(),
                maximum_mach: 0.2,
            },
            |_| false,
            initial,
        )
        .unwrap()
    }
    fn bits<const Q: usize>(populations: &[[f64; Q]]) -> Vec<u64> {
        populations.iter().flatten().map(|v| v.to_bits()).collect()
    }

    fn constant_law<L: ViscousLattice3<Q>, const Q: usize>(model: L::Collision) -> L::Collision {
        struct Source {
            cell: usize,
            state: FluidTransportState,
            reads: Rc<Cell<usize>>,
        }
        impl FluidTransport for Source {
            fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
                self.reads.set(self.reads.get() + 1);
                Ok(self.state)
            }
        }
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let make = || {
            Grid3::<L, Q>::try_new(
                [3, 4, 2],
                model,
                [1e-6, -5e-7, 0.0],
                [true; 3],
                GridBudget3 {
                    maximum_cells: 24,
                    maximum_mach: 0.2,
                },
                |_| false,
                |_| L::W,
            )
            .unwrap()
        };
        let mut actual = make();
        let mut manual = make();
        let tau = 0.5 + reference.flow().viscosity / L::CS2;
        let expected = L::with_shear_relaxation(model, tau).unwrap();
        let reads = Rc::new(Cell::new(0));
        let queries = Cell::new(0);
        let report = with_cx(|cx| {
            let report = actual
                .step_serial_with_material_rheology(cx, &reference, |point, _| {
                    queries.set(queries.get() + 1);
                    Ok(Source {
                        cell: point.cell,
                        state: material,
                        reads: reads.clone(),
                    })
                })
                .unwrap();
            manual
                .step_serial_with_collisions(cx, &vec![expected; 24])
                .unwrap();
            report
        });
        assert_eq!(queries.get(), 24);
        assert_eq!(reads.get(), 24);
        assert_eq!(bits(actual.populations()), bits(manual.populations()));
        assert_eq!(actual.completed_steps(), 1);
        assert_eq!(report.sampled_at_step, 0);
        for (cell, sample) in report.samples.iter().enumerate() {
            let sample = sample.as_ref().unwrap();
            assert_eq!(sample.source().cell, cell);
            assert_eq!(sample.transport(), &material);
            assert_eq!(sample.point().time_s, 0.0);
        }
        report.collisions[0]
    }

    #[test]
    fn newtonian_material_reduces_to_existing_step_and_retains_original_sources() {
        for model in [
            CollisionModel3::Bgk { tau: 0.9 },
            CollisionModel3::CentralMoment {
                second_order_rate: 1.0 / 0.9,
                higher_order_rate: 1.2,
            },
            CollisionModel3::ReducedCumulant {
                second_order_rate: 1.0 / 0.9,
                third_order_rate: 1.2,
                fourth_order_rate: 0.8,
            },
        ] {
            constant_law::<D3q19, 19>(model);
        }
        for model in 0..4 {
            let updated = constant_law::<D3q27, 27>(collision(model));
            if let (
                CollisionModel27::CentralMoment { rates: before }
                | CollisionModel27::Cumulant { rates: before },
                CollisionModel27::CentralMoment { rates: after }
                | CollisionModel27::Cumulant { rates: after },
            ) = (collision(model).model, updated.model)
            {
                assert_eq!(before.omega[1..], after.omega[1..]);
            }
        }
    }

    #[test]
    fn carreau_queries_use_incoming_force_centered_strain_and_ignore_rotation() {
        let material = fluid(0.05);
        let dt = 0.02;
        let reference = plan_scaling_from_material(&material, 0.001, dt, 0.005).unwrap();
        for model in [2, 3] {
            for shear in [0.0, 0.006] {
                let force = |p: [usize; 3]| {
                    [
                        0.003 * (p[1] as f64 - 2.0),
                        0.002 * (p[0] as f64 - 2.0),
                        0.0,
                    ]
                };
                let grid = grid27([5, 5, 3], model, [0.0; 3], [true; 3], |p| {
                    let x = p[0] as f64 - 2.0;
                    let y = p[1] as f64 - 2.0;
                    let physical = [(shear - 0.004) * y, 0.004 * x, 0.0];
                    equilibrium27(
                        1.0,
                        core::array::from_fn(|a| physical[a] - force(p)[a] / 2.0),
                        Equilibrium27::Factorized,
                    )
                });
                let forces = (0..75)
                    .map(|cell| force(coordinates(cell, [5, 5, 3])))
                    .collect();
                let mut run = with_cx(|cx| ForcedGrid3::try_new(grid, cx, forces).unwrap());
                let carreau = |rate: f64| 0.005 + 0.045 * (1.0 + (5.0 * rate).powi(2)).powf(-0.3);
                let report = with_cx(|cx| {
                    run.advance_serial_with_material_rheology(
                        cx,
                        &reference,
                        vec![[0.0; 3]; 75],
                        |point, _| Ok(fluid(carreau(point.shear_rate_per_s))),
                    )
                })
                .unwrap();
                let center = 2 + 5 * (2 + 5);
                let sample = report.samples[center].as_ref().unwrap();
                let expected = shear / dt;
                assert!((sample.point().shear_rate_per_s - expected).abs() < 2e-14);
                assert!(
                    (sample.transport().dynamic_viscosity_pa_s() - carreau(expected)).abs() < 2e-14
                );
                assert_eq!(run.forces(), vec![[0.0; 3]; 75]);
            }
        }
    }

    #[test]
    fn wall_gradient_recovers_quadratic_physical_shear_at_every_cell() {
        let material = fluid(0.01);
        let dt = 0.02;
        let reference = plan_scaling_from_material(&material, 0.001, dt, 0.005).unwrap();
        let force = [8e-4, 0.0, 0.0];
        let mut grid = grid27([1, 9, 1], 3, force, [true, false, true], |p| {
            let y = p[1] as f64 + 0.5;
            equilibrium27(
                1.0,
                [0.001 * y * (9.0 - y) - force[0] / 2.0, 0.0, 0.0],
                Equilibrium27::Factorized,
            )
        });
        let report = with_cx(|cx| {
            grid.step_serial_with_material_rheology(cx, &reference, |_, _| Ok(material))
        })
        .unwrap();
        for (cell, sample) in report.samples.iter().enumerate() {
            let expected = (0.001 * (9.0 - 2.0 * (cell as f64 + 0.5))).abs() / dt;
            assert!(
                (sample.as_ref().unwrap().point().shear_rate_per_s - expected).abs() < 2e-14,
                "wall/interior shear at cell {cell}"
            );
        }
    }

    fn power_law_channel(n: usize) -> f64 {
        let exponent = 0.75;
        let height = 0.01;
        let rho = 1000.0;
        let consistency = 0.02;
        let wall_rate: f64 = 0.1;
        let mu_wall = consistency * wall_rate.powf(exponent - 1.0);
        let material = fluid(mu_wall);
        let dx = height / n as f64;
        let dt = 0.04 * dx * dx / (mu_wall / rho);
        let reference = plan_scaling_from_material(&material, dx, dt, 0.005).unwrap();
        let force_si = 2.0 * consistency * wall_rate.powf(exponent) / height;
        let force = reference.force_density_lattice(force_si).unwrap();
        // Independent power-law solution of (K |u'|^(n-1) u')' = -F, u(0)=u(H)=0.
        let peak = exponent / (exponent + 1.0) * wall_rate * height / 2.0;
        let exact = |y: f64| {
            peak * (1.0 - ((y - height / 2.0).abs() / (height / 2.0)).powf(1.0 + 1.0 / exponent))
        };
        // The analytic velocity is only an initial guess. Equilibrium populations
        // do not supply steady stress; actual collision/streaming must settle.
        let mut grid = grid27([1, n, 1], 2, [force, 0.0, 0.0], [true, false, true], |p| {
            let u = exact((p[1] as f64 + 0.5) * dx) * dt / dx;
            equilibrium27(1.0, [u - force / 2.0, 0.0, 0.0], Equilibrium27::Factorized)
        });
        let mut previous = grid.macroscopic(n / 2).unwrap().1[0];
        let mut converged = false;
        with_cx(|cx| {
            for step in 1..=100 * n * n {
                grid.step_serial_with_material_rheology(cx, &reference, |point, _| {
                    // This source's explicit low-rate regularization is never active:
                    // every actual query must remain on the analytic power-law branch.
                    assert!(point.shear_rate_per_s > 1e-8);
                    Ok(fluid(
                        consistency * point.shear_rate_per_s.max(1e-8).powf(exponent - 1.0),
                    ))
                })
                .unwrap();
                if step % (n * n) == 0 {
                    let current = grid.macroscopic(n / 2).unwrap().1[0];
                    let change = (current - previous).abs() * dx / dt / peak;
                    if step >= 4 * n * n && change < 2e-7 {
                        converged = true;
                        break;
                    }
                    previous = current;
                }
            }
        });
        assert!(converged, "steady-state budget n={n}");
        let mut error = 0.0_f64;
        for cell in 0..n {
            let (density, u) = grid.macroscopic(cell).unwrap();
            assert!((density - 1.0).abs() < 1e-9);
            let observed = reference.velocity_si(u[0]).unwrap();
            error = error.max((observed - exact((cell as f64 + 0.5) * dx)).abs() / peak);
        }
        error
    }

    #[test]
    fn central_moment_material_power_law_channel_refines_against_analytic_profile() {
        let coarse = power_law_channel(12);
        let fine = power_law_channel(24);
        eprintln!("material power-law CM n=0.75: {coarse} -> {fine}");
        assert!(
            coarse < 0.05 && fine < 0.025,
            "profile errors {coarse} -> {fine}"
        );
        assert!(
            fine < 0.65 * coarse,
            "non-refining physical profile {coarse} -> {fine}"
        );
    }

    #[test]
    fn source_refusal_and_callback_cancellation_leave_the_grid_unchanged() {
        let material = fluid(0.01);
        let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
        let make = || grid27([4, 1, 1], 3, [0.0; 3], [true; 3], |_| D3q27::W);
        let mut grid = make();
        let before = bits(grid.populations());
        let result = with_cx(|cx| {
            grid.step_serial_with_material_rheology(cx, &reference, |point, _| {
                if point.cell == 3 {
                    FluidTransportState::try_new(1000.0, 4000.0, 0.5, -1.0)
                } else {
                    Ok(material)
                }
            })
        });
        assert!(matches!(
            result,
            Err(RheologyError3::Source { cell: 3, .. })
        ));
        assert_eq!(bits(grid.populations()), before);
        assert_eq!(grid.completed_steps(), 0);
        let gate = CancelGate::new_clock_free();
        let queries = Cell::new(0);
        let result = with_gate(&gate, |cx| {
            grid.step_serial_with_material_rheology(cx, &reference, |point, _| {
                queries.set(queries.get() + 1);
                if point.cell == 2 {
                    gate.request();
                }
                Ok(material)
            })
        });
        assert!(matches!(
            result,
            Err(RheologyError3::Grid(GridError3::Cancelled))
        ));
        assert_eq!(queries.get(), 3);
        assert_eq!(bits(grid.populations()), before);
        assert_eq!(grid.completed_steps(), 0);
    }
}
