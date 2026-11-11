use super::*;
use crate::nsf3::{NsfThermalWallFace, NsfWallThermal};

fn wall_plane(axis: usize, high: bool, thermal: NsfWallThermal) -> Vec<NsfThermalWallFace> {
    let mut walls = Vec::new();
    for z in 0..SHAPE[2] {
        for y in 0..SHAPE[1] {
            for x in 0..SHAPE[0] {
                let cell = [x, y, z];
                if cell[axis] == if high { SHAPE[axis] - 1 } else { 0 } {
                    walls.push(NsfThermalWallFace {
                        cell,
                        axis,
                        high,
                        thermal,
                    });
                }
            }
        }
    }
    walls
}

fn patterned_walls() -> Vec<NsfThermalWallFace> {
    let mut walls = wall_plane(2, false, NsfWallThermal::Isothermal(293.25));
    for wall in &mut walls {
        wall.thermal = match (wall.cell[0] + 2 * wall.cell[1]) % 3 {
            0 => NsfWallThermal::Isothermal(293.25),
            1 => NsfWallThermal::HeatFlux(20.0),
            _ => NsfWallThermal::Film {
                coefficient_w_m2_k: 25.0,
                ambient_temperature_k: 293.1,
            },
        };
    }
    walls.extend(wall_plane(
        2,
        true,
        NsfWallThermal::Film {
            coefficient_w_m2_k: 25.0,
            ambient_temperature_k: 293.05,
        },
    ));
    walls
}

// The original filter, supplied independently with the physical policy for
// each plane. Its existing channel MMS checks both three-point derivatives.
fn original_fit(grid: &NsfGrid3<D3q27>, walls: [ScalarWall; 2]) -> Vec<NsfSgsTransport> {
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
    let mut flow = ChannelDynamicSgs::new(SHAPE, 1.0, 2, usize::MAX).unwrap();
    let mut heat = ChannelDynamicScalar::new(SHAPE, 2, usize::MAX).unwrap();
    flow.prepare_with_scalar(
        config().momentum,
        config().heat,
        &rho,
        &velocity,
        grid.edge_velocity[2],
        &grid.theta,
        walls,
        &mut heat,
        || false,
    )
    .unwrap();
    grid.dynamic_transport_field(
        &rho,
        &grid.theta,
        &grid.material,
        flow.viscosities(),
        heat.diffusivities(),
        config().momentum.maximum_eddy_to_molecular,
    )
    .unwrap()
}

#[test]
fn uniform_local_walls_preserve_original_dynamic_trajectory_and_distinct_plane_values() {
    let make = || {
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Isothermal(293.15)), |_| false, mixed).unwrap()
    };
    let mut walls = wall_plane(2, false, NsfWallThermal::Isothermal(293.15));
    walls.extend(wall_plane(2, true, NsfWallThermal::Isothermal(293.15)));
    let mut plain = make()
        .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap();
    let mut local = make()
        .with_thermal_walls(walls.clone())
        .unwrap()
        .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap();
    let mut reversed = make()
        .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap()
        .with_thermal_walls(walls)
        .unwrap();
    for _ in 0..3 {
        assert_eq!(local.sgs_transport, plain.sgs_transport);
        assert_eq!(reversed.sgs_transport, plain.sgs_transport);
        let ledger = plain.step().unwrap();
        assert_eq!(local.step().unwrap(), ledger);
        assert_eq!(reversed.step().unwrap(), ledger);
        assert_eq!(local.f, plain.f);
        assert_eq!(local.e, plain.e);
        assert_eq!(local.theta, plain.theta);
        assert_eq!(local.div, plain.div);
        assert_eq!(reversed.state_hash(), local.state_hash());
    }
    let mut distinct = wall_plane(2, false, NsfWallThermal::Isothermal(293.25));
    distinct.extend(wall_plane(2, true, NsfWallThermal::Isothermal(293.05)));
    local.set_thermal_walls(distinct).unwrap();
    let declared = [293.25, 293.05].map(|t| ScalarWall::Value(local.law.theta_of_kelvin(t)));
    assert_eq!(local.channel_scalar_walls(2), declared);
    assert_eq!(local.sgs_transport, original_fit(&local, declared));
    local.set_thermal_walls(Vec::new()).unwrap();
    assert_eq!(local.state_hash(), plain.state_hash());
    assert_eq!(
        local.sgs_transport,
        original_fit(
            &local,
            [ScalarWall::Value(local.law.theta_of_kelvin(293.15)); 2]
        )
    );
}

#[test]
fn heterogeneous_flux_and_film_walls_use_both_native_scalar_levels_and_one_heat_owner() {
    let make = |workers, walls_first| {
        let grid = NsfGrid3::<D3q27>::new(
            spec(workers, WallThermal::Isothermal(293.15)),
            |_| false,
            mixed,
        )
        .unwrap();
        if walls_first {
            grid.with_thermal_walls(patterned_walls())
                .unwrap()
                .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
                .unwrap()
        } else {
            grid.with_channel_dynamic_thermal_sgs(config(), usize::MAX)
                .unwrap()
                .with_thermal_walls(patterned_walls())
                .unwrap()
        }
    };
    let mut serial = make(1, true);
    let mut parallel = make(3, false);
    let expected = original_fit(&serial, [ScalarWall::Interior; 2]);
    assert_eq!(serial.sgs_transport, expected);
    assert_eq!(parallel.sgs_transport, expected);
    assert_ne!(
        expected,
        original_fit(
            &serial,
            [ScalarWall::Value(serial.law.theta_of_kelvin(293.15)); 2]
        ),
        "a uniform ghost must not silently replace the patterned wall"
    );
    let allocation = {
        let work = workspace(&serial);
        (
            work.momentum.viscosities().as_ptr(),
            work.heat.diffusivities().as_ptr(),
        )
    };
    for _ in 0..4 {
        let (rho, momentum) = serial.macro_fields(&serial.f);
        let velocity = serial.velocity_field(&rho, &momentum, &serial.force);
        let no_work = vec![[0.0; 6]; serial.n];
        let molecular = serial
            .flux_divergence(
                &serial.e,
                &rho,
                &velocity,
                &serial.theta,
                &no_work,
                &[],
                None,
            )
            .unwrap();
        let modeled = serial
            .flux_divergence(
                &serial.e,
                &rho,
                &velocity,
                &serial.theta,
                &no_work,
                &serial.sgs_transport,
                None,
            )
            .unwrap();
        assert_eq!(modeled.1, molecular.1);
        assert_eq!(modeled.2.wall_heat_out, molecular.2.wall_heat_out);
        let before = serial.total_energy_j();
        let ledger = serial.step().unwrap();
        assert_eq!(ledger, parallel.step().unwrap());
        assert_eq!(serial.state_hash(), parallel.state_hash());
        assert_eq!(serial.sgs_transport, parallel.sgs_transport);
        assert!(ledger.residual_rel < 3e-12);
        assert!(
            (serial.total_energy_j() - before - ledger.wall_heat_j - ledger.wall_work_j).abs()
                < before.abs() * 3e-12
        );
        let work = workspace(&serial);
        assert_eq!(work.momentum.viscosities().as_ptr(), allocation.0);
        assert_eq!(work.heat.diffusivities().as_ptr(), allocation.1);
    }
}

#[test]
fn local_wall_refit_is_transactional_and_uses_the_current_accepted_state() {
    let mut grid =
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Isothermal(293.15)), |_| false, mixed)
            .unwrap()
            .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
            .unwrap();
    grid.step().unwrap();
    let physical = (
        grid.f.clone(),
        grid.e.clone(),
        grid.theta.clone(),
        grid.force.clone(),
        grid.steps(),
        grid.cumulative_ledger(),
    );
    grid.set_thermal_walls(patterned_walls()).unwrap();
    assert_eq!(
        (
            grid.f.clone(),
            grid.e.clone(),
            grid.theta.clone(),
            grid.force.clone(),
            grid.steps(),
            grid.cumulative_ledger()
        ),
        physical
    );
    assert_eq!(
        grid.sgs_transport,
        original_fit(&grid, [ScalarWall::Interior; 2])
    );
    let accepted = (
        grid.state_hash(),
        grid.sgs_transport.clone(),
        grid.div.clone(),
        grid.wall_out,
        grid.open_out.wall_heat_out,
        grid.cumulative_ledger(),
    );
    let mut polls = 0;
    assert_eq!(
        grid.set_thermal_walls_cancellable(Vec::new(), || {
            polls += 1;
            polls == 3
        }),
        Err(NsfError::Input("thermal wall update cancelled"))
    );
    assert_eq!(polls, 3);
    assert_eq!(
        (
            grid.state_hash(),
            grid.sgs_transport.clone(),
            grid.div.clone(),
            grid.wall_out,
            grid.open_out.wall_heat_out,
            grid.cumulative_ledger()
        ),
        accepted
    );
    // Candidate scalar coefficients must pass the existing physical cap before
    // either they or the changed wall table can replace the accepted history.
    // Use the original nonzero initial fit as the trigger, independently of
    // whether a freely decaying endpoint retains a positive scalar coefficient.
    let mut grid =
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Isothermal(293.15)), |_| false, mixed)
            .unwrap()
            .with_thermal_walls(patterned_walls())
            .unwrap()
            .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
            .unwrap();
    assert!(
        original_fit(
            &grid,
            [ScalarWall::Value(grid.law.theta_of_kelvin(293.15)); 2]
        )
        .iter()
        .any(|value| value.turbulent_conductivity_lu > 0.0)
    );
    let accepted = (
        grid.state_hash(),
        grid.sgs_transport.clone(),
        grid.div.clone(),
        grid.wall_out,
        grid.open_out.wall_heat_out,
        grid.cumulative_ledger(),
    );
    let DynamicState::Channel { config, .. } = &mut grid.dynamic_sgs.as_mut().unwrap().state else {
        unreachable!();
    };
    config.heat.maximum_diffusivity_to_molecular = 0.0;
    assert!(matches!(
        grid.set_thermal_walls(Vec::new()),
        Err(NsfError::SgsHeatRatio { .. })
    ));
    assert_eq!(
        (
            grid.state_hash(),
            grid.sgs_transport.clone(),
            grid.div.clone(),
            grid.wall_out,
            grid.open_out.wall_heat_out,
            grid.cumulative_ledger()
        ),
        accepted
    );
}

#[test]
fn plane_policy_retains_defaults_and_exact_insulation_in_every_orientation() {
    for axis in 0..3 {
        let mut settings = spec(1, WallThermal::Adiabatic);
        settings.periodic = [true; 3];
        settings.periodic[axis] = false;
        let mut grid = NsfGrid3::<D3q27>::new(settings, |_| false, rest)
            .unwrap()
            .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
            .unwrap();
        for high in [false, true] {
            let mut partial = wall_plane(axis, high, NsfWallThermal::HeatFlux(-0.0));
            partial.truncate(1);
            grid.set_thermal_walls(partial.clone()).unwrap();
            assert_eq!(
                grid.channel_scalar_walls(axis),
                [ScalarWall::ZeroGradient; 2]
            );
            partial[0].thermal = NsfWallThermal::Film {
                coefficient_w_m2_k: 0.0,
                ambient_temperature_k: 300.0,
            };
            grid.set_thermal_walls(partial.clone()).unwrap();
            assert_eq!(
                grid.channel_scalar_walls(axis),
                [ScalarWall::ZeroGradient; 2]
            );
            partial[0].thermal = NsfWallThermal::Isothermal(300.0);
            grid.set_thermal_walls(partial).unwrap();
            let mut expected = [ScalarWall::ZeroGradient; 2];
            expected[usize::from(high)] = ScalarWall::Interior;
            assert_eq!(grid.channel_scalar_walls(axis), expected);
            grid.set_thermal_walls(wall_plane(axis, high, NsfWallThermal::Isothermal(300.0)))
                .unwrap();
            expected[usize::from(high)] = ScalarWall::Value(grid.law.theta_of_kelvin(300.0));
            assert_eq!(grid.channel_scalar_walls(axis), expected);
        }
    }
}

#[test]
fn g3_finite_wall_trials_keep_dynamic_transport_and_rollback_filter_scratch() {
    let make = || {
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Isothermal(293.15)), |_| false, mixed)
            .unwrap()
            .with_channel_dynamic_thermal_sgs(config(), 200 * SHAPE.iter().product::<usize>())
            .unwrap()
    };
    let mut grid = make();
    let mut direct = make();
    assert!(
        grid.sgs_transport
            .iter()
            .any(|t| t.turbulent_conductivity_lu > 0.0)
    );
    assert_eq!(grid.step().unwrap(), direct.step().unwrap());
    for _ in 0..2 {
        let accepted = (
            grid.state_hash(),
            grid.sgs_transport.clone(),
            grid.div.clone(),
            grid.cumulative_ledger(),
        );
        let scratch = {
            let work = workspace(&grid);
            (
                work.momentum.viscosities().to_vec(),
                work.heat.diffusivities().to_vec(),
                work.momentum.viscosities().as_ptr(),
                work.heat.diffusivities().as_ptr(),
            )
        };
        let bytes = grid.dynamic_sgs_workspace_bytes();
        // A trial rebases its starting FV history on every call (`set_thermal_walls_recorded`, record = true;
        // thermal_wall_trial.rs:73-75, "including the explicit starting-history rebase"), while the public setter returns
        // early when the table is unchanged (thermal_wall.rs:202). From the second pass on, the public setter would leave
        // this oracle on the un-rebased end-of-step history, a different numerical path that differs from the trial at
        // corrector tolerance (internal_j 1e-9 relative, measured).
        direct
            .set_thermal_walls_recorded(patterned_walls(), || false)
            .unwrap();
        let expected = direct.step().unwrap();
        let mut complete_polls = 0;
        let trial = grid
            .preview_thermal_wall_step(patterned_walls(), 100, || {
                complete_polls += 1;
                false
            })
            .unwrap();
        assert_eq!(*trial.energy_ledger(), expected);
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
        assert_eq!(
            trial.grid().channel_dynamic_thermal_sgs_config(),
            Some(config())
        );
        assert_eq!(trial.grid().dynamic_sgs_workspace_bytes(), bytes);
        let face_heat: f64 = trial
            .wall_exchange()
            .iter()
            .map(|f| f.heat_into_fluid_j)
            .sum();
        assert!((face_heat - expected.wall_heat_j).abs() < 1e-22);
        drop(trial);

        // Cancel after the complete native step, immediately before its trial
        // could be returned. Neither physical history nor original scratch moves.
        let mut polls = 0;
        assert!(
            grid.preview_thermal_wall_step(patterned_walls(), 100, || {
                polls += 1;
                polls >= complete_polls
            })
            .is_err()
        );
        assert_eq!(polls, complete_polls);
        assert_eq!(
            (
                grid.state_hash(),
                grid.sgs_transport.clone(),
                grid.div.clone(),
                grid.cumulative_ledger()
            ),
            accepted
        );
        {
            // Since ef5ab372af (2026-11-07 01:58, nine minutes after dc7e920856 wrote this test) a trial shares the grid's
            // already-admitted filter scratch instead of copying it (`shared_for_trial`, dynamic.rs: "each combined
            // preparation overwrites its deterministic filter scratch"), so a cancelled trial may leave the scratch VALUES
            // overwritten. What this block checks: the physical history is unchanged (asserted just above) and the scratch
            // is the same two allocations, with no second filter allocation (`scratch.2` and `scratch.3`).
            let work = workspace(&grid);
            assert_eq!(work.momentum.viscosities().as_ptr(), scratch.2);
            assert_eq!(work.heat.diffusivities().as_ptr(), scratch.3);
        }
        grid.preview_thermal_wall_step(patterned_walls(), 100, || false)
            .unwrap()
            .accept();
        assert_eq!(grid.state_hash(), direct.state_hash());
        assert_eq!(grid.sgs_transport, direct.sgs_transport);
    }
}
