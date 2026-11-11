//! Material/source epochs and pathline history must advance together.

use super::*;

fn config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 8.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn scratch(grid: &Grid3<D3q27, 27>, channel: bool) -> LagrangianSgsWorkspace3<CollisionConfig27> {
    if channel {
        grid.channel_lagrangian_sgs_workspace(usize::MAX).unwrap()
    } else {
        grid.lagrangian_sgs_workspace(usize::MAX).unwrap()
    }
}

#[test]
fn pathline_material_transport_matches_native_composition_and_worker_counts() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let pool = fs_exec::TilePool::for_host(2, 93);
    let gate = CancelGate::new_clock_free();
    for channel in [false, true] {
        let mut serial = fluctuating(channel);
        let mut pooled = fluctuating(channel);
        let mut native = fluctuating(channel);
        let mut a = serial.rheology_workspace().unwrap();
        let mut b = pooled.rheology_workspace().unwrap();
        let mut molecular = native.rheology_workspace().unwrap();
        let mut sa = scratch(&serial, channel);
        let mut sb = scratch(&pooled, channel);
        let mut sn = scratch(&native, channel);
        let addresses = storage_addresses(&a);
        let bytes = sa.workspace_bytes();
        let pointer = sa.collisions().as_ptr();
        let mut saw_transport = false;
        with_cx(|cx| {
            for step in 0..4 {
                // Independently compose the existing source sampler and native
                // history-bearing step. The source law changes with both the
                // current physical shear and elapsed SI time.
                prepare_workspace(
                    &native,
                    &reference,
                    None,
                    &mut molecular,
                    |point, _| Ok(fluid(viscosity(point))),
                    || cx.checkpoint().map_err(|_| GridError3::Cancelled),
                )
                .unwrap();
                let expected = native
                    .step_serial_with_lagrangian_sgs(
                        cx,
                        config(),
                        Some(&molecular.trial.collisions),
                        &mut sn,
                    )
                    .unwrap();
                let first = serial
                    .step_serial_with_material_rheology_and_sgs(
                        cx,
                        &reference,
                        &mut a,
                        RheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sa,
                        },
                        |point, _| Ok(fluid(viscosity(point))),
                    )
                    .unwrap();
                let second = pooled
                    .step_pooled_with_material_rheology_and_sgs(
                        &pool,
                        &gate,
                        &reference,
                        &mut b,
                        RheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sb,
                        },
                        |point, _| Ok(fluid(viscosity(point))),
                    )
                    .unwrap();
                assert_eq!(first.sgs, RheologySgsReport3::Lagrangian(expected));
                assert_eq!(first.sgs, second.sgs);
                assert_eq!(first.rheology.collisions, second.rheology.collisions);
                assert_eq!(first.rheology.collisions, sn.collisions());
                saw_transport |= expected.fit.maximum_eddy_viscosity > 0.0;
                for (cell, sample) in first.rheology.samples.iter().enumerate() {
                    let sample = sample.as_ref().unwrap();
                    let original = molecular.trial.samples[cell].as_ref().unwrap();
                    assert_eq!(sample.point(), original.point());
                    assert_eq!(sample.point().time_s, step as f64 * reference.timestep_s());
                    assert_eq!(
                        sample.transport().dynamic_viscosity_pa_s(),
                        viscosity(sample.point())
                    );
                    assert_eq!(
                        sample.shear_relaxation_time(),
                        original.shear_relaxation_time()
                    );
                    let base =
                        D3q27::molecular_viscosity(molecular.trial.collisions[cell]).unwrap();
                    let total =
                        D3q27::molecular_viscosity(first.rheology.collisions[cell]).unwrap();
                    assert!(total >= base - 2e-16);
                }
                assert_eq!(bits(serial.populations()), bits(native.populations()));
                assert_eq!(bits(serial.populations()), bits(pooled.populations()));
                assert_eq!(sa.accepted_fit_time(), Some(step as f64));
                assert_eq!(sb.accepted_fit_time(), Some(step as f64));
            }
        });
        assert!(saw_transport);
        assert_eq!(serial.completed_steps(), 4);
        assert_eq!(addresses, storage_addresses(&a));
        assert_eq!(bytes, sa.workspace_bytes());
        assert_eq!(pointer, sa.collisions().as_ptr());
    }
}

#[test]
fn rejected_material_and_force_endpoints_preserve_pathline_retry() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    with_cx(|cx| {
        let force: Vec<_> = (0..100)
            .map(|cell| [1e-5 * (cell % 5) as f64, 0.0, 0.0])
            .collect();
        let make = || ForcedGrid3::try_new(fluctuating(false), cx, force.clone()).unwrap();
        let mut actual = make();
        let mut control = make();
        let mut a = actual.rheology_workspace().unwrap();
        let mut b = control.rheology_workspace().unwrap();
        let mut sa = actual.lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut sb = control.lagrangian_sgs_workspace(usize::MAX).unwrap();
        for (grid, source, history) in [
            (&mut actual, &mut a, &mut sa),
            (&mut control, &mut b, &mut sb),
        ] {
            grid.step_serial_with_material_rheology_and_sgs(
                cx,
                &reference,
                source,
                RheologySgs3::Lagrangian {
                    config: config(),
                    workspace: history,
                },
                |point, _| Ok(fluid(viscosity(point))),
            )
            .unwrap();
        }
        let before = bits(actual.populations());
        let accepted = a.report().unwrap().collisions.clone();
        let mut queries = 0;
        // The native endpoint is rejected AFTER source sampling and a valid
        // candidate history, exercising the real post-fit publication boundary.
        assert!(matches!(
            actual.advance_serial_with_material_rheology_and_sgs(
                cx,
                &reference,
                vec![[1.0, 0.0, 0.0]; 100],
                &mut a,
                RheologySgs3::Lagrangian {
                    config: config(),
                    workspace: &mut sa
                },
                |point, _| {
                    queries += 1;
                    Ok(fluid(viscosity(point)))
                },
            ),
            Err(RheologySgsError3::Sgs(SgsStepError::Grid(
                GridError3::Mach { .. }
            )))
        ));
        assert_eq!(queries, 100);
        assert!(
            actual
                .step_serial_with_material_rheology_and_sgs(
                    cx,
                    &reference,
                    &mut a,
                    RheologySgs3::Lagrangian {
                        config: config(),
                        workspace: &mut sa
                    },
                    |point, _| FluidTransportState::try_new(
                        1000.0,
                        4000.0,
                        0.5,
                        if point.cell == 99 {
                            -1.0
                        } else {
                            viscosity(point)
                        }
                    ),
                )
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_gate(&cancelled, |cx| {
            assert!(
                actual
                    .step_serial_with_material_rheology_and_sgs(
                        cx,
                        &reference,
                        &mut a,
                        RheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sa
                        },
                        |point, _| Ok(fluid(viscosity(point))),
                    )
                    .is_err()
            );
        });
        assert_eq!(bits(actual.populations()), before);
        assert_eq!(actual.forces(), force);
        assert_eq!(actual.completed_steps(), 1);
        assert_eq!(sa.accepted_fit_time(), Some(0.0));
        assert_eq!(a.report().unwrap().collisions, accepted);
        assert_eq!(a.report().unwrap().sampled_at_step, 0);
        let next = vec![[-2e-5, 1e-5, 0.0]; 100];
        let retry = actual
            .advance_serial_with_material_rheology_and_sgs(
                cx,
                &reference,
                next.clone(),
                &mut a,
                RheologySgs3::Lagrangian {
                    config: config(),
                    workspace: &mut sa,
                },
                |point, _| Ok(fluid(viscosity(point))),
            )
            .unwrap();
        let direct = control
            .advance_pooled_with_material_rheology_and_sgs(
                &fs_exec::TilePool::for_host(2, 93),
                &CancelGate::new_clock_free(),
                &reference,
                next.clone(),
                &mut b,
                RheologySgs3::Lagrangian {
                    config: config(),
                    workspace: &mut sb,
                },
                |point, _| Ok(fluid(viscosity(point))),
            )
            .unwrap();
        assert_eq!(retry.sgs, direct.sgs);
        assert_eq!(retry.rheology.collisions, direct.rheology.collisions);
        assert_eq!(bits(actual.populations()), bits(control.populations()));
        assert_eq!(actual.forces(), next);
        assert_eq!(sa.accepted_fit_time(), Some(1.0));
    });
}

#[test]
fn moving_material_walls_and_history_share_one_force_epoch() {
    use fs_material::fluid_transport::rheology::{
        CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
    };
    use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
    let material = fluid(0.02);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let law = ThermalCarreauYasuda::new(
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
    );
    let source = |point: RheologyPoint3, _: &FluidTransportState| {
        law.evaluate(300.0 + point.time_s, point.shear_rate_per_s)?
            .bind(material)
    };
    with_cx(|cx| {
        for normal in 0..3 {
            let tangent = (normal + 1) % 3;
            let make = || {
                let grid = grid27(
                    [5; 3],
                    3,
                    [0.0; 3],
                    core::array::from_fn(|a| a != normal),
                    |p| {
                        let cell = (p[2] * 5 + p[1]) * 5 + p[0];
                        equilibrium27(
                            1.0,
                            core::array::from_fn(|a| 0.007 * random(cell, a as u64 + 7)),
                            Equilibrium27::Factorized,
                        )
                    },
                );
                ForcedGrid3::try_new(grid, cx, vec![[1e-5, -1e-5, 2e-5]; 125]).unwrap()
            };
            let mut serial = make();
            let mut pooled = make();
            let mut native = make();
            let probe = serial
                .grid
                .wall_momentum_probe(cx, [0.0; 3], |_| true)
                .unwrap();
            let mut a = serial.rheology_workspace().unwrap();
            let mut b = pooled.rheology_workspace().unwrap();
            let mut molecular = native.rheology_workspace().unwrap();
            let mut sa = serial.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
            let mut sb = pooled.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
            let mut sn = native.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
            for step in 0..3 {
                let mut walls = [[0.0; 3]; 2];
                walls[0][tangent] = 0.003;
                walls[1][tangent] = 0.009 + 0.001 * step as f64;
                let drive: Vec<_> = probe
                    .wall_links()
                    .map(|link| walls[usize::from(link.incoming_direction[normal] < 0)])
                    .collect();
                let endpoint = vec![[1e-5 * (step + 1) as f64, 2e-5, -1e-5]; 125];
                prepare_workspace_with_walls(
                    &native.grid,
                    &reference,
                    Some(native.forces()),
                    Some((normal, walls)),
                    &mut molecular,
                    source,
                    || cx.checkpoint().map_err(|_| GridError3::Cancelled),
                )
                .unwrap();
                let expected = native
                    .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                        cx,
                        config(),
                        Some(&molecular.trial.collisions),
                        &probe,
                        &drive,
                        &mut sn,
                        endpoint.clone(),
                    )
                    .unwrap();
                let first = serial
                    .advance_serial_with_channel_material_rheology_sgs_and_moving_walls(
                        cx,
                        &reference,
                        &probe,
                        &drive,
                        endpoint.clone(),
                        &mut a,
                        ChannelRheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sa,
                        },
                        source,
                    )
                    .unwrap();
                let second = pooled
                    .advance_pooled_with_channel_material_rheology_sgs_and_moving_walls(
                        &fs_exec::TilePool::for_host(2, 93),
                        &CancelGate::new_clock_free(),
                        &reference,
                        &probe,
                        &drive,
                        endpoint.clone(),
                        &mut b,
                        ChannelRheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sb,
                        },
                        source,
                    )
                    .unwrap();
                assert_eq!(first.0.sgs, RheologySgsReport3::Lagrangian(expected.0));
                assert_eq!(first.0.sgs, second.0.sgs);
                assert_eq!(first.1, expected.1);
                assert_eq!(first.1, second.1);
                assert_eq!(first.0.rheology.collisions, sn.collisions());
                assert_eq!(bits(serial.populations()), bits(native.populations()));
                assert_eq!(bits(serial.populations()), bits(pooled.populations()));
                assert_eq!(serial.forces(), endpoint);
                assert_eq!(sa.accepted_fit_time(), Some(step as f64));
                assert_eq!(sb.accepted_fit_time(), Some(step as f64));
            }
            // Wall-normal motion is refused before it can change any source
            // or history epoch, even after several successful moving steps.
            let before = bits(serial.populations());
            let accepted = a.report().unwrap().collisions.clone();
            let bad = vec![[0.01; 3]; probe.measured_links()];
            assert!(
                serial
                    .step_serial_with_channel_material_rheology_sgs_and_moving_walls(
                        cx,
                        &reference,
                        &probe,
                        &bad,
                        &mut a,
                        ChannelRheologySgs3::Lagrangian {
                            config: config(),
                            workspace: &mut sa
                        },
                        source,
                    )
                    .is_err()
            );
            assert_eq!(bits(serial.populations()), before);
            assert_eq!(a.report().unwrap().collisions, accepted);
            assert_eq!(sa.accepted_fit_time(), Some(2.0));
        }
    });
}
