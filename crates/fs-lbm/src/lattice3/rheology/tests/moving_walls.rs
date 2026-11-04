use super::*;
use crate::lattice3::wall_loads::WallMomentumProbe3;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law() -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(0.02, 0.003, 8.0, 0.4, 2.0, [0.0, 100.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn probe<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> WallMomentumProbe3<L, Q> {
    with_cx(|cx| grid.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap())
}

fn velocities<L: Lattice3<Q>, const Q: usize>(
    probe: &WallMomentumProbe3<L, Q>,
    walls: [[f64; 3]; 2],
) -> Vec<[f64; 3]> {
    probe
        .wall_links()
        .map(|link| walls[usize::from(link.incoming_direction[1] < 0)])
        .collect()
}

fn analytic_couette<L: ViscousLattice3<Q>, const Q: usize>(model: L::Collision) {
    let material = fluid(0.02);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let law = law();
    for height in [1, 6] {
        for boost in [0.0, 0.02] {
            let lower = [boost, 0.0, -0.5 * boost];
            let delta = [0.012, 0.0, -0.009];
            let upper = core::array::from_fn(|a| lower[a] + delta[a]);
            let force = [2e-5, 0.0, 0.0];
            let mut grid = Grid3::<L, Q>::try_new(
                [2, height, 2],
                model,
                force,
                [true, false, true],
                GridBudget3 {
                    maximum_cells: 4 * height,
                    maximum_mach: 0.2,
                },
                |_| false,
                |p| {
                    // Raw momentum must recover the same-time half-force
                    // physical velocity before the material rate is sampled.
                    let u: [f64; 3] = core::array::from_fn(|a| {
                        lower[a] + delta[a] * (p[1] as f64 + 0.5) / height as f64 - 0.5 * force[a]
                    });
                    let uu = u.iter().map(|v| v * v).sum::<f64>();
                    core::array::from_fn(|q| {
                        let eu = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum::<f64>();
                        L::W[q]
                            * (1.0 + eu / L::CS2 + 0.5 * eu * eu / L::CS2.powi(2)
                                - 0.5 * uu / L::CS2)
                    })
                },
            )
            .unwrap();
            let probe = probe(&grid);
            let drive = velocities(&probe, [lower, upper]);
            let mut workspace = grid.rheology_workspace().unwrap();
            let queries = Cell::new(0);
            with_cx(|cx| {
                let (report, load) = grid
                    .step_serial_with_channel_material_rheology_and_moving_walls(
                        cx,
                        &reference,
                        &probe,
                        &drive,
                        &mut workspace,
                        |point, _| {
                            queries.set(queries.get() + 1);
                            law.evaluate(300.0, point.shear_rate_per_s)?.bind(material)
                        },
                    )
                    .unwrap();
                // Independent analytic Couette + Carreau equation, including
                // both tangential shear components and physical time scaling.
                let rate = 0.015 / (height as f64 * reference.timestep_s());
                let mu = 0.003 + 0.017 * (1.0 + (8.0 * rate).powi(2)).powf(-0.3);
                let tau = 0.5 + (mu / 1000.0) * 0.01 / (L::CS2 * 0.001_f64.powi(2));
                for sample in report.samples.iter().flatten() {
                    assert!((sample.point().shear_rate_per_s - rate).abs() < 2e-12);
                    assert!((sample.transport().dynamic_viscosity_pa_s() - mu).abs() < 2e-14);
                    assert!((sample.shear_relaxation_time() - tau).abs() < 1e-12);
                    assert_eq!(sample.source().rheology().law(), &law);
                    assert_eq!(sample.source().original_transport(), &material);
                }
                assert!(load.fluid_mass_change.abs() < 2e-14);
                assert_eq!(load.measured_links, probe.measured_links());
            });
            assert_eq!(queries.get(), 4 * height);
            assert_eq!(grid.completed_steps(), 1);
        }
    }
}

#[test]
fn moving_wall_couette_uses_frame_independent_material_viscosity() {
    analytic_couette::<D3q19, 19>(CollisionModel3::CentralMoment {
        second_order_rate: 1.0 / 0.9,
        higher_order_rate: 1.2,
    });
    analytic_couette::<D3q27, 27>(collision(3));
}

#[test]
fn moving_material_channel_matches_pooled_native_transport_and_reuses_storage() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let law = law();
    let make = || grid27([3, 5, 2], 3, [0.0; 3], [true, false, true], |_| D3q27::W);
    let mut serial = make();
    let mut pooled = make();
    let mut native = make();
    let mut stationary = make();
    let probe = probe(&serial);
    let mut serial_workspace = serial.rheology_workspace().unwrap();
    let mut pooled_workspace = pooled.rheology_workspace().unwrap();
    let addresses = storage_addresses(&serial_workspace);
    let pool = fs_exec::TilePool::for_host(2, 0x3D_4E0);
    let gate = CancelGate::new_clock_free();
    let source = |point: RheologyPoint3, _: &FluidTransportState| {
        law.evaluate(300.0, point.shear_rate_per_s)?.bind(material)
    };
    for step in 0..4 {
        let drive = velocities(
            &probe,
            [[0.0; 3], [0.004 * step as f64, 0.0, -0.001 * step as f64]],
        );
        let serial_load = with_cx(|cx| {
            serial
                .step_serial_with_channel_material_rheology_and_moving_walls(
                    cx,
                    &reference,
                    &probe,
                    &drive,
                    &mut serial_workspace,
                    source,
                )
                .unwrap()
                .1
        });
        let pooled_load = pooled
            .step_pooled_with_channel_material_rheology_and_moving_walls(
                &pool,
                &gate,
                &reference,
                &probe,
                &drive,
                &mut pooled_workspace,
                source,
            )
            .unwrap()
            .1;
        let report = serial_workspace.report().unwrap();
        let native_load = with_cx(|cx| {
            native
                .step_serial_with_moving_wall_momentum(cx, &probe, &drive, Some(&report.collisions))
                .unwrap()
        });
        assert_eq!(serial_load, pooled_load);
        assert_eq!(serial_load, native_load);
        assert_eq!(bits(serial.populations()), bits(pooled.populations()));
        assert_eq!(bits(serial.populations()), bits(native.populations()));
        assert_eq!(storage_addresses(&serial_workspace), addresses);
        assert_eq!(report.sampled_at_step, step);
        if step == 0 {
            with_cx(|cx| {
                stationary
                    .step_serial_with_material_rheology(cx, &reference, source)
                    .unwrap()
            });
            assert_eq!(bits(serial.populations()), bits(stationary.populations()));
        } else {
            assert!(
                report.samples[12]
                    .as_ref()
                    .unwrap()
                    .point()
                    .shear_rate_per_s
                    > 0.0
            );
            assert!(
                report.samples[12]
                    .as_ref()
                    .unwrap()
                    .transport()
                    .dynamic_viscosity_pa_s()
                    < 0.02
            );
        }
        assert_eq!(
            report.collisions,
            pooled_workspace.report().unwrap().collisions
        );
        for local in &report.collisions {
            let CollisionModel27::Cumulant { rates } = local.model else {
                panic!("collision family changed")
            };
            let CollisionModel27::Cumulant { rates: original } = collision(3).model else {
                unreachable!()
            };
            assert_eq!(&rates.omega[1..], &original.omega[1..]);
        }
    }
}

#[test]
fn moving_material_channel_refusals_preserve_flow_and_source_epoch() {
    use fs_exec::{RunError, RunReport, TileKernel, TilePool};
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct CancelAfterStream {
        pool: TilePool,
        calls: AtomicUsize,
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
            if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
                gate.request();
            }
            result
        }
    }
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let make = || grid27([3, 5, 2], 3, [0.0; 3], [true, false, true], |_| D3q27::W);
    let mut grid = make();
    let mut control = make();
    let probe = probe(&grid);
    let drive = velocities(&probe, [[0.0; 3], [0.02, 0.0, 0.0]]);
    let mut workspace = grid.rheology_workspace().unwrap();
    let mut control_workspace = control.rheology_workspace().unwrap();
    let source = |_: RheologyPoint3, _: &FluidTransportState| Ok(material);
    with_cx(|cx| {
        let _ = grid
            .step_serial_with_channel_material_rheology_and_moving_walls(
                cx,
                &reference,
                &probe,
                &drive,
                &mut workspace,
                source,
            )
            .unwrap();
        let _ = control
            .step_serial_with_channel_material_rheology_and_moving_walls(
                cx,
                &reference,
                &probe,
                &drive,
                &mut control_workspace,
                source,
            )
            .unwrap();
    });
    let before = bits(grid.populations());
    let accepted_points: Vec<_> = workspace
        .report()
        .unwrap()
        .samples
        .iter()
        .map(|s| s.as_ref().unwrap().point())
        .collect();
    let addresses = storage_addresses(&workspace);
    let unchanged =
        |grid: &Grid3<D3q27, 27>,
         workspace: &RheologyWorkspace3<CollisionConfig27, FluidTransportState>| {
            assert_eq!(bits(grid.populations()), before);
            assert_eq!(grid.completed_steps(), 1);
            let report = workspace.report().unwrap();
            assert_eq!(report.sampled_at_step, 0);
            for (sample, expected) in report.samples.iter().zip(&accepted_points) {
                assert_eq!(sample.as_ref().unwrap().point(), *expected);
                assert_eq!(sample.as_ref().unwrap().transport(), &material);
            }
            assert_eq!(storage_addresses(workspace), addresses);
        };
    let subset = with_cx(|cx| {
        grid.wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[3])
            .unwrap()
    });
    let subset_drive = vec![[0.02, 0.0, 0.0]; subset.measured_links()];
    let queries = Cell::new(0);
    assert!(with_cx(|cx| grid
        .step_serial_with_channel_material_rheology_and_moving_walls(
            cx,
            &reference,
            &subset,
            &subset_drive,
            &mut workspace,
            |_, _| {
                queries.set(queries.get() + 1);
                Ok(material)
            },
        )
        .is_err()));
    assert_eq!(queries.get(), 0);
    unchanged(&grid, &workspace);
    let runner = CancelAfterStream {
        pool: TilePool::for_host(2, 0x3D_4E0),
        calls: AtomicUsize::new(0),
    };
    let changed = |_: RheologyPoint3, _: &FluidTransportState| Ok(fluid(0.03));
    let result = grid.step_pooled_with_channel_material_rheology_and_moving_walls(
        &runner,
        &CancelGate::new_clock_free(),
        &reference,
        &probe,
        &drive,
        &mut workspace,
        changed,
    );
    assert!(matches!(
        result,
        Err(RheologyError3::Grid(GridError3::Cancelled))
    ));
    unchanged(&grid, &workspace);
    let actual_load = grid
        .step_pooled_with_channel_material_rheology_and_moving_walls(
            &runner.pool,
            &CancelGate::new_clock_free(),
            &reference,
            &probe,
            &drive,
            &mut workspace,
            changed,
        )
        .unwrap()
        .1;
    let expected_load = with_cx(|cx| {
        control
            .step_serial_with_channel_material_rheology_and_moving_walls(
                cx,
                &reference,
                &probe,
                &drive,
                &mut control_workspace,
                changed,
            )
            .unwrap()
            .1
    });
    assert_eq!(actual_load, expected_load);
    assert_eq!(bits(grid.populations()), bits(control.populations()));
    assert_eq!(workspace.report().unwrap().sampled_at_step, 1);
}
