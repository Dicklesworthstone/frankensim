//! Source/flow integration controls, not calibrated turbulent rheology.

use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_exec::{CancelGate, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySgsError3, BoundarySpec3, BoundaryStepError3, CollisionModel3, E3,
    FaceBoundary3, Q3, collide_cell3, equilibrium3,
};
use fs_lbm::material::boundary::sgs::MaterialBoundarySgsError3;
use fs_lbm::material::immersed::MaterialImmersedPoint3;
use fs_lbm::material::{MaterialLbmError, plan_scaling_from_material};
use fs_lbm::sgs3::{SgsConfig, SgsError, SgsModel, SgsUse};
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const ORIGIN: [f64; 3] = [1.25, -2.5, 0.75];

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1.0, mu).unwrap()
}

fn settings(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn model(tau: f64, family: usize) -> CollisionModel3 {
    match family {
        0 => CollisionModel3::Bgk { tau },
        1 => CollisionModel3::Trt { tau, magic: 0.137 },
        2 => CollisionModel3::CentralMoment {
            second_order_rate: 1.0 / tau,
            higher_order_rate: 0.73,
        },
        _ => CollisionModel3::ReducedCumulant {
            second_order_rate: 1.0 / tau,
            third_order_rate: 0.87,
            fourth_order_rate: 1.63,
        },
    }
}

fn seed(grid: &mut BoundaryGrid3, velocity: impl Fn([usize; 3]) -> [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let u = velocity([x, y, z]);
                let force = grid.force_density(x, y, z);
                let mut f = equilibrium3(1.0, core::array::from_fn(|a| u[a] - 0.5 * force[a]));
                // Nonzero stress with exactly zero added mass and momentum.
                f[1] += 0.00012;
                f[2] += 0.00012;
                f[0] -= 0.00024;
                grid.set_populations(x, y, z, &f);
            }
        }
    }
}

fn fixture(
    dimensions: [usize; 3],
    boundaries: BoundarySpec3,
    law: CollisionModel3,
) -> BoundaryGrid3 {
    let [nx, ny, nz] = dimensions;
    let mut grid =
        BoundaryGrid3::with_collision_model(nx, ny, nz, law, [1e-6, -2e-6, 0.0], boundaries);
    grid.voxelize_sdf(|p| (p[0] - 4.0).powi(2) + (p[1] - 4.0).powi(2) + (p[2] - 4.0).powi(2) - 1.5);
    grid.set_force_field(|p| [1e-6 * (p[1] as f64 + 1.0), 0.0, -3e-7]);
    seed(&mut grid, |p| {
        [
            0.002 * p[1] as f64,
            -0.001 * p[2] as f64,
            0.0007 * p[0] as f64,
        ]
    });
    grid
}

fn bits(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    (0..nx * ny * nz)
        .map(|i| {
            grid.populations(i % nx, i / nx % ny, i / (nx * ny))
                .map(f64::to_bits)
        })
        .collect()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

struct SingleRead {
    point: MaterialImmersedPoint3,
    reads: Cell<usize>,
}
impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        assert_eq!(self.reads.replace(self.reads.get() + 1), 0);
        // Original heat transport must not impose a thermal lattice constraint.
        FluidTransportState::try_new(1000.0, 1000.0, 1e-250, 1.0)
    }
}

#[test]
fn constant_material_recovers_native_sgs_and_wall_loads_for_all_collision_families() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    for boundaries in [
        BoundarySpec3::lid_cavity([0.025, 0.0, 0.01]),
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.02], 1.0),
    ] {
        for family in 0..4 {
            for closure in [
                SgsModel::Wale { cw: 0.0 },
                SgsModel::Smagorinsky { cs: 0.18 },
                SgsModel::Vreman { cv: 0.07 },
            ] {
                let original = model(1.37, family);
                let mut actual = fixture([8; 3], boundaries, original);
                let mut expected = fixture([8; 3], boundaries, model(scaling.flow().tau, family));
                let links = actual.boundary_links();
                let mut work = scaling
                    .boundary_rheology_sgs_workspace(&actual, 512, &gate)
                    .unwrap();
                let bytes = work.workspace_bytes();
                let mut native = expected.sgs_workspace().unwrap();
                for _ in 0..2 {
                    let native_report = expected
                        .step_with_sgs(
                            &gate,
                            settings(closure),
                            &mut native,
                            Some((&links, ORIGIN)),
                        )
                        .unwrap();
                    let mut calls = 0;
                    let report = scaling
                        .step_boundary_with_material_rheology_and_sgs(
                            &mut actual,
                            &gate,
                            settings(closure),
                            &mut work,
                            Some((&links, ORIGIN)),
                            |point| {
                                calls += 1;
                                Ok(SingleRead {
                                    point,
                                    reads: Cell::new(0),
                                })
                            },
                        )
                        .unwrap();
                    assert_eq!(report, native_report);
                    assert_eq!(bits(&actual), bits(&expected));
                    assert_eq!(actual.collision_model(), original);
                    assert_eq!(calls, report.sgs.fluid_cells);
                    assert_eq!(work.applied().len(), calls);
                    assert_eq!(work.last_report(), Some(report));
                    assert_eq!(work.workspace_bytes(), bytes);
                    for sample in work.applied() {
                        assert_eq!(sample.point, sample.material.point);
                        assert_eq!(sample.material.reads.get(), 1);
                        assert_eq!(sample.shear_relaxation_time, scaling.flow().tau);
                        let [x, y, z] = sample.point.cell;
                        assert!(!actual.is_solid(x, y, z));
                    }
                }
            }
        }
    }
}

#[test]
fn local_molecular_and_eddy_viscosities_drive_the_actual_collision_and_periodic_transport() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let velocity = [0.01, 0.02, 0.01, 0.0];
    let derivative = [0.01, 0.0, -0.01, 0.0];
    let force = [2e-6, -1e-6, 5e-7];
    for family in 0..4 {
        for cs in [0.0, 0.3] {
            let mut grid = BoundaryGrid3::with_collision_model(
                4,
                4,
                4,
                model(1.37, family),
                force,
                BoundarySpec3::periodic(),
            );
            seed(&mut grid, |p| [velocity[p[1]], 0.0, 0.0]);
            let mut next = vec![[0.0; Q3]; 64];
            for z in 0..4 {
                for y in 0..4 {
                    for x in 0..4 {
                        let rate = f64::abs(derivative[y]) / DT;
                        let mu = 0.4 + 0.2 * x as f64 + 0.1 * rate;
                        let molecular = mu / 1000.0 * DT / (DX * DX);
                        let eddy = cs * cs * f64::abs(derivative[y]);
                        let post = collide_cell3(
                            grid.populations(x, y, z),
                            model(0.5 + 3.0 * (molecular + eddy), family),
                            force,
                        )
                        .unwrap();
                        // Independently push each post-collision population.
                        for (q, e) in E3.iter().enumerate() {
                            let to = [
                                (x as i32 + i32::from(e.0)).rem_euclid(4) as usize,
                                (y as i32 + i32::from(e.1)).rem_euclid(4) as usize,
                                (z as i32 + i32::from(e.2)).rem_euclid(4) as usize,
                            ];
                            next[(to[2] * 4 + to[1]) * 4 + to[0]][q] = post[q];
                        }
                    }
                }
            }
            let mut work = scaling
                .boundary_rheology_sgs_workspace(&grid, 64, &gate)
                .unwrap();
            scaling
                .step_boundary_with_material_rheology_and_sgs(
                    &mut grid,
                    &gate,
                    settings(SgsModel::Smagorinsky { cs }),
                    &mut work,
                    None,
                    |p| {
                        close(
                            p.shear_rate_per_s,
                            f64::abs(derivative[p.cell[1]]) / DT,
                            2e-13,
                        );
                        Ok(fluid(
                            0.4 + 0.2 * p.cell[0] as f64 + 0.1 * p.shear_rate_per_s,
                        ))
                    },
                )
                .unwrap();
            for (i, expected) in next.iter().enumerate() {
                for (a, b) in grid
                    .populations(i % 4, i / 4 % 4, i / 16)
                    .into_iter()
                    .zip(expected)
                {
                    close(a, *b, 3e-15);
                }
            }
        }
    }
}

#[test]
fn physical_couette_shear_drives_source_viscosity_and_local_sgs_dissipation() {
    let reference = fluid(1.0);
    let gate = CancelGate::new_clock_free();
    let periodic = FaceBoundary3::Periodic;
    let boundaries = BoundarySpec3::new([
        periodic,
        periodic,
        FaceBoundary3::Wall {
            velocity: [0.01, 0.0, 0.0],
        },
        FaceBoundary3::Wall {
            velocity: [0.05, 0.0, 0.0],
        },
        periodic,
        periodic,
    ]);
    let shear = 0.04 / 8.0;
    for dt in [DT, 2.0 * DT] {
        let scaling = plan_scaling_from_material(&reference, DX, dt, 0.05).unwrap();
        let mut grid = BoundaryGrid3::new(4, 8, 4, 1.37, [1e-5, 0.0, 0.0], boundaries);
        grid.set_force_field(|p| [2e-5 * (p[1] as f64 + 1.0), 0.0, 0.0]);
        seed(&mut grid, |p| {
            [0.01 + shear * (p[1] as f64 + 0.5), 0.0, 0.0]
        });
        let mut config = settings(SgsModel::Smagorinsky { cs: 0.2 });
        config.filter_width_lu = 1.5;
        let mut work = scaling
            .boundary_rheology_sgs_workspace(&grid, 128, &gate)
            .unwrap();
        let report = scaling
            .step_boundary_with_material_rheology_and_sgs(
                &mut grid,
                &gate,
                config,
                &mut work,
                None,
                |p| {
                    close(p.shear_rate_per_s, shear / dt, 2e-13);
                    Ok(fluid(
                        0.4 + 0.1 * p.shear_rate_per_s + 0.02 * p.cell[1] as f64,
                    ))
                },
            )
            .unwrap();
        let eddy = (0.2_f64 * 1.5).powi(2) * shear;
        let mut molecular_dissipation = 0.0;
        let mut maximum_ratio = 0.0_f64;
        for sample in work.applied() {
            let mu = 0.4 + 0.1 * shear / dt + 0.02 * sample.point.cell[1] as f64;
            let nu = mu / 1000.0 * dt / (DX * DX);
            close(sample.shear_relaxation_time, 0.5 + 3.0 * nu, 2e-14);
            close(sample.transport.dynamic_viscosity_pa_s(), mu, 2e-14);
            molecular_dissipation += nu * shear * shear;
            maximum_ratio = maximum_ratio.max(eddy / nu);
        }
        close(
            report.sgs.resolved_dissipation_lu,
            molecular_dissipation,
            3e-15,
        );
        close(
            report.sgs.sgs_dissipation_lu,
            128.0 * eddy * shear * shear,
            3e-15,
        );
        close(report.sgs.maximum_eddy_to_molecular, maximum_ratio, 3e-14);
    }
}

#[test]
fn source_ratio_and_cancellation_refusals_preserve_the_accepted_source_epoch() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut grid = fixture(
        [8; 3],
        BoundarySpec3::lid_cavity([0.02, 0.0, 0.0]),
        model(1.37, 0),
    );
    let config = settings(SgsModel::Smagorinsky { cs: 0.18 });
    let mut work = scaling
        .boundary_rheology_sgs_workspace(&grid, 512, &gate)
        .unwrap();
    scaling
        .step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            None,
            |_| Ok(fluid(0.8)),
        )
        .unwrap();
    let before = bits(&grid);
    let applied = work.applied().to_vec();
    let accepted = work.last_report();
    let source_error = MaterialError::State {
        what: "original source refused this state".into(),
    };
    let mut calls = 0;
    let failed = scaling.step_boundary_with_material_rheology_and_sgs(
        &mut grid,
        &gate,
        config,
        &mut work,
        None,
        |_| {
            calls += 1;
            if calls == 6 {
                Err(source_error.clone())
            } else {
                Ok(fluid(1.2))
            }
        },
    );
    assert!(
        matches!(failed, Err(MaterialBoundarySgsError3::Source { source: MaterialLbmError::Transport(e), .. }) if e == source_error)
    );
    assert_eq!(calls, 6);
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    assert!(matches!(
        scaling.step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            None,
            |_| FluidTransportState::try_new(1001.0, 1000.0, 1.0, 1.2)
        ),
        Err(MaterialBoundarySgsError3::Source {
            source: MaterialLbmError::Invalid(_),
            ..
        })
    ));
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    let mut capped = config;
    capped.maximum_eddy_to_molecular = 0.0;
    assert!(matches!(
        scaling.step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            capped,
            &mut work,
            None,
            |_| Ok(fluid(1.2))
        ),
        Err(MaterialBoundarySgsError3::Sgs(BoundarySgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        )))
    ));
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    let cancelled = CancelGate::new_clock_free();
    let mut calls = 0;
    assert!(matches!(
        scaling.step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &cancelled,
            config,
            &mut work,
            None,
            |_| {
                calls += 1;
                if calls == applied.len() {
                    cancelled.request();
                }
                Ok(fluid(1.2))
            }
        ),
        Err(MaterialBoundarySgsError3::Sgs(BoundarySgsError3::Boundary(
            BoundaryStepError3::Cancelled
        )))
    ));
    assert_eq!(calls, applied.len());
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    assert_eq!(work.last_report(), accepted);
    let other_scaling = plan_scaling_from_material(&reference, DX * 2.0, DT, 0.05).unwrap();
    assert!(matches!(
        other_scaling.step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            None,
            |_| panic!("incompatible units sampled")
        ),
        Err(MaterialBoundarySgsError3::Invalid(_))
    ));
    let links = grid.boundary_links();
    let duplicated = [links[0], links[0]];
    assert!(
        scaling
            .step_boundary_with_material_rheology_and_sgs(
                &mut grid,
                &gate,
                config,
                &mut work,
                Some((&duplicated, ORIGIN)),
                |_| panic!("invalid wall selection sampled")
            )
            .is_err()
    );
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    scaling
        .step_boundary_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            None,
            |_| Ok(fluid(1.2)),
        )
        .unwrap();
    assert!(work.applied().iter().all(|s| s.transport == fluid(1.2)));
}

struct CancelAfterStream {
    pool: TilePool,
    passes: AtomicUsize,
}
impl KernelRunner for CancelAfterStream {
    fn workers(&self) -> usize {
        self.pool.workers()
    }
    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let result = self.pool.run_with_gate(kernel, gate);
        if self.passes.fetch_add(1, Ordering::Relaxed) == 1 {
            gate.request();
        }
        result
    }
}

#[test]
fn pooled_wall_loads_replay_serial_and_late_native_cancellation_is_atomic() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let config = settings(SgsModel::Wale { cw: 0.5 });
    let resolve = |p: MaterialImmersedPoint3| {
        Ok(fluid(
            0.3 + 0.1 * p.shear_rate_per_s + 0.01 * p.cell[0] as f64,
        ))
    };
    for workers in [1, 3] {
        let make = || {
            fixture(
                [16, 8, 8],
                BoundarySpec3::lid_cavity([0.025, 0.0, 0.01]),
                model(1.37, 3),
            )
        };
        let mut serial = make();
        let mut pooled = make();
        let links = serial.boundary_links();
        let mut sw = scaling
            .boundary_rheology_sgs_workspace(&serial, 1024, &gate)
            .unwrap();
        let mut pw = scaling
            .boundary_rheology_sgs_workspace(&pooled, 1024, &gate)
            .unwrap();
        let pool = TilePool::for_host(workers, 853);
        for _ in 0..2 {
            let sr = scaling
                .step_boundary_with_material_rheology_and_sgs(
                    &mut serial,
                    &gate,
                    config,
                    &mut sw,
                    Some((&links, ORIGIN)),
                    resolve,
                )
                .unwrap();
            let pr = scaling
                .step_boundary_pooled_with_material_rheology_and_sgs(
                    &mut pooled,
                    &pool,
                    &gate,
                    config,
                    &mut pw,
                    Some((&links, ORIGIN)),
                    resolve,
                )
                .unwrap();
            assert_eq!(sr, pr);
            assert_eq!(bits(&serial), bits(&pooled));
            assert_eq!(sw.applied(), pw.applied());
        }
        let before = bits(&pooled);
        let applied = pw.applied().to_vec();
        let accepted = pw.last_report();
        let cancelled = CancelGate::new_clock_free();
        let stop = CancelAfterStream {
            pool,
            passes: AtomicUsize::new(0),
        };
        assert!(matches!(
            scaling.step_boundary_pooled_with_material_rheology_and_sgs(
                &mut pooled,
                &stop,
                &cancelled,
                config,
                &mut pw,
                Some((&links, ORIGIN)),
                resolve
            ),
            Err(MaterialBoundarySgsError3::Sgs(BoundarySgsError3::Boundary(
                BoundaryStepError3::Cancelled
            )))
        ));
        assert_eq!(bits(&pooled), before);
        assert_eq!(pw.applied(), applied);
        assert_eq!(pw.last_report(), accepted);
        scaling
            .step_boundary_with_material_rheology_and_sgs(
                &mut serial,
                &gate,
                config,
                &mut sw,
                Some((&links, ORIGIN)),
                resolve,
            )
            .unwrap();
        scaling
            .step_boundary_pooled_with_material_rheology_and_sgs(
                &mut pooled,
                &stop,
                &gate,
                config,
                &mut pw,
                Some((&links, ORIGIN)),
                resolve,
            )
            .unwrap();
        assert_eq!(bits(&serial), bits(&pooled));
        assert_eq!(sw.applied(), pw.applied());
        assert_eq!(sw.last_report(), pw.last_report());
    }
}
