use super::super::super::tests::{build_grid, config as frozen_config, snapshot};
use super::*;
use crate::d3q27::D3q27;
use core::cell::Cell;

fn source(grid: &NsfGrid3<D3q27>, spatial: usize, wall: bool) -> NsfEvolvingMannInflow {
    let config = frozen_config();
    let snapshot = HomogeneousMannField::at_anchor(
        snapshot(2.0).model().clone(),
        79,
        vec![2.5 * grid.law.dt_s; 16],
    )
    .unwrap();
    let field =
        FixedStepMannField::from_snapshot(snapshot, 531, grid.law.dt_s, usize::MAX).unwrap();
    let bytes = NsfEvolvingMannInflow::required_workspace_bytes(&field, 9).unwrap();
    let mut inlet = NsfEvolvingMannInflow::new(
        field,
        NsfEvolvingMannConfig {
            pressure_pa: config.pressure_pa,
            temperature_k: config.temperature_k,
            mean_velocity_m_s: config.mean_velocity_m_s,
            grid_origin_in_field_m: config.grid_origin_in_field_m,
            field_to_grid: config.field_to_grid,
            maximum_mode_samples: config.maximum_mode_samples,
        },
        grid,
        0,
        bytes + 16 * 32,
    )
    .unwrap();
    if wall {
        inlet = inlet
            .with_plane_wall(
                MannPlaneWall::new(
                    (config.grid_origin_in_field_m[2] + 3.0 * grid.law.dx_m).next_down(),
                )
                .unwrap(),
            )
            .unwrap();
    }
    match spatial {
        1 => inlet.with_face_average(16 * 32, || false).unwrap(),
        2 => inlet.with_cell_average(16 * 32, || false).unwrap(),
        _ => inlet,
    }
}

fn budget(steps: usize, catch_up: usize) -> NsfEvolvingMannAverageBudget {
    NsfEvolvingMannAverageBudget {
        maximum_cell_steps: 729 * steps,
        maximum_normal_draws: 16 * (8 * steps + 4 * catch_up),
        maximum_mode_samples: 9 * 16 * (steps + 1),
    }
}

fn drives() -> [NsfBoundaryDrive; 2] {
    [0.02, -0.01].map(NsfBoundaryDrive::IncomingCharacteristic)
}

// Explicit fixture rotation and actual ghost locations, independent of the
// adapter's generic frame/profile construction and time-loop implementation.
fn manual_profile(
    grid: &NsfGrid3<D3q27>,
    inlet: &NsfEvolvingMannInflow,
    mean: &MannStepAverage,
) -> Vec<NsfReservoirTileState> {
    let c = inlet.config();
    let gas = grid.spec.gas.at(c.temperature_k, c.pressure_pa).unwrap();
    grid.reservoir_aperture_tile_states(0)
        .unwrap()
        .map(|tile| {
            let point = [
                c.grid_origin_in_field_m[0] - 0.5 * grid.law.dx_m,
                c.grid_origin_in_field_m[1] - (tile.cell[0] as f64 + 0.5) * grid.law.dx_m,
                c.grid_origin_in_field_m[2] + (tile.cell[2] as f64 + 0.5) * grid.law.dx_m,
            ];
            let u = if let Some(filter) = inlet.wall_spatial_filter() {
                mean.sample_wall_blocked_filtered(point, filter)
            } else if let Some(wall) = inlet.plane_wall() {
                mean.sample_wall_blocked(point, wall)
            } else if let Some(filter) = inlet.spatial_filter() {
                mean.sample_filtered(point, filter)
            } else {
                mean.sample(point)
            }
            .unwrap()
            .u;
            NsfReservoirTileState {
                cell: tile.cell,
                state: NsfOpenState {
                    density_kg_m3: gas.density,
                    temperature_k: gas.temperature,
                    velocity_m_s: [-u[1], c.mean_velocity_m_s[1] + u[0], u[2]],
                },
            }
        })
        .collect()
}

#[test]
fn g1_native_time_means_match_manual_ghost_controls_and_original_mass_energy_accounts() {
    for wall in [false, true] {
        for spatial in 0..3 {
            let mut direct = build_grid(1, false);
            let mut direct_inlet = source(&direct, spatial, wall);
            let mut direct_mean =
                MannStepAverage::prepare(direct_inlet.field(), 16 * 104, || false).unwrap();
            direct.set_boundary_drives(drives()).unwrap();
            let mut exchanges = Vec::new();
            for _ in 0..3 {
                direct_inlet
                    .field
                    .advance_one_with_average(&mut direct_mean, 16 * 8, || false)
                    .unwrap();
                let profile = manual_profile(&direct, &direct_inlet, &direct_mean);
                direct.set_reservoir_aperture_profile(0, &profile).unwrap();
                assert!(direct.step().unwrap().residual_rel < 1e-10);
                exchanges.push(*direct.last_reservoir_exchange_at(0).unwrap());
            }
            direct_inlet.apply(&mut direct, 0).unwrap();

            let mut grid = build_grid(1, false);
            let mass_before = grid.total_mass_kg();
            let mut inlet = source(&grid, spatial, wall);
            let mut mean = MannStepAverage::prepare(inlet.field(), 16 * 104, || false).unwrap();
            let profile_ptr = inlet.profile.as_ptr();
            let source_bytes = inlet.workspace_bytes();
            let trial = inlet
                .preview_time_averaged_boundary_interval(
                    &mut grid,
                    &mut mean,
                    drives(),
                    3,
                    budget(3, 0),
                    || false,
                )
                .unwrap();
            assert_eq!(trial.fluid().grid().state_hash(), direct.state_hash());
            assert_eq!(trial.fluid().grid().f, direct.f);
            assert_eq!(trial.fluid().grid().e, direct.e);
            assert_eq!(trial.fluid().grid().force, direct.force);
            assert_eq!(
                trial.field().checkpoint(),
                direct_inlet.field().checkpoint()
            );
            assert_eq!(mean, direct_mean);
            let report = trial.report();
            assert_eq!(report.steps, 3);
            assert_eq!(report.profile_updates, 4);
            assert_eq!(report.normal_draws, 3 * 16 * 8);
            assert_eq!(report.mode_samples, 4 * 9 * 16);
            assert_eq!(report.atmosphere_step_range, [0, 3]);
            assert_eq!(report.endpoint.atmosphere_step, 3);
            assert_eq!(report.endpoint.ou_mode_transitions, 0);
            assert_eq!(report.cfd_time_range_s, [0.0, 3.0 * direct.law.dt_s]);
            let exchange = trial.fluid().reservoir_exchange_at(0).unwrap();
            assert_eq!(
                exchange.mass_out_kg,
                exchanges.iter().map(|x| x.mass_out_kg).sum::<f64>()
            );
            assert_eq!(
                exchange.total_energy_out_j,
                exchanges.iter().map(|x| x.total_energy_out_j).sum::<f64>()
            );
            assert_eq!(
                exchange.pressure_work_out_j,
                exchanges.iter().map(|x| x.pressure_work_out_j).sum::<f64>()
            );
            assert_eq!(
                exchange.acoustic_pressure_work_out_j,
                exchanges
                    .iter()
                    .map(|x| x.acoustic_pressure_work_out_j)
                    .sum::<f64>()
            );
            let mass_out = exchange.mass_out_kg
                + trial
                    .fluid()
                    .boundary_exchange()
                    .iter()
                    .map(|x| x.mass_out_kg)
                    .sum::<f64>();
            assert!((direct.total_mass_kg() - mass_before + mass_out).abs() < 1e-10 * mass_before);
            assert!(trial.fluid().energy_ledger().residual_rel < 1e-10);
            trial.accept();
            assert_eq!(grid.cumulative_ledger(), direct.cumulative_ledger());
            assert_eq!(inlet.profile.as_ptr(), profile_ptr);
            assert_eq!(inlet.workspace_bytes(), source_bytes);
            let hash = grid.state_hash();
            assert_eq!(inlet.apply(&mut grid, 0).unwrap().ou_mode_transitions, 0);
            assert_eq!(grid.state_hash(), hash);
        }
    }
}

#[test]
fn g4_mean_interval_rolls_back_budgets_cancellation_drop_and_replays_after_catch_up() {
    let make = || {
        let mut grid = build_grid(1, false);
        let inlet = source(&grid, 2, true);
        // Source remains at its bound epoch while the original grid advances.
        grid.step().unwrap();
        grid.step().unwrap();
        let mean = MannStepAverage::prepare(inlet.field(), 16 * 104, || false).unwrap();
        (grid, inlet, mean)
    };
    let (mut expected, mut expected_inlet, mut expected_mean) = make();
    let polls = Cell::new(0);
    let expected_trial = expected_inlet
        .preview_time_averaged_boundary_interval(
            &mut expected,
            &mut expected_mean,
            drives(),
            2,
            budget(2, 2),
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let expected_report = expected_trial.report();
    assert_eq!(expected_report.atmosphere_step_range, [2, 4]);
    assert_eq!(expected_report.normal_draws, 16 * (2 * 4 + 2 * 8));
    expected_trial.accept();
    assert!(polls.get() > 10);
    for stop in [1, 5, polls.get() / 2, polls.get()] {
        let (mut grid, mut inlet, mut mean) = make();
        let before = grid.state_hash();
        let checkpoint = inlet.field().checkpoint();
        let count = Cell::new(0);
        assert!(matches!(
            inlet.preview_time_averaged_boundary_interval(
                &mut grid,
                &mut mean,
                drives(),
                2,
                budget(2, 2),
                || {
                    count.set(count.get() + 1);
                    count.get() >= stop
                },
            ),
            Err(NsfEvolvingMannError::Cancelled)
        ));
        assert_eq!(grid.state_hash(), before);
        assert_eq!(inlet.field().checkpoint(), checkpoint);
        let trial = inlet
            .preview_time_averaged_boundary_interval(
                &mut grid,
                &mut mean,
                drives(),
                2,
                budget(2, 2),
                || false,
            )
            .unwrap();
        assert_eq!(trial.report(), expected_report);
        assert_eq!(trial.fluid().grid().state_hash(), expected.state_hash());
        assert_eq!(
            trial.field().checkpoint(),
            expected_inlet.field().checkpoint()
        );
        drop(trial);
        assert_eq!(grid.state_hash(), before);
        assert_eq!(inlet.field().checkpoint(), checkpoint);
    }
    for limit in 0..3 {
        let (mut grid, mut inlet, mut mean) = make();
        let before = grid.state_hash();
        let checkpoint = inlet.field().checkpoint();
        let mut b = budget(2, 2);
        match limit {
            0 => b.maximum_cell_steps -= 1,
            1 => b.maximum_normal_draws -= 1,
            _ => b.maximum_mode_samples -= 1,
        }
        assert!(
            inlet
                .preview_time_averaged_boundary_interval(
                    &mut grid,
                    &mut mean,
                    drives(),
                    2,
                    b,
                    || false
                )
                .is_err()
        );
        assert_eq!(grid.state_hash(), before);
        assert_eq!(inlet.field().checkpoint(), checkpoint);
        if limit == 0 {
            assert_eq!(
                inlet.trial.checkpoint(),
                checkpoint,
                "native cell-work admission must precede atmosphere catch-up"
            );
            assert_eq!(mean.end_step(), None);
        }
    }
}

#[test]
fn g1_integral_velocity_changes_native_flow_without_changing_endpoint_realization() {
    let mut ordinary = build_grid(1, false);
    let mut ordinary_inlet = source(&ordinary, 2, true);
    let ordinary_trial = ordinary_inlet
        .preview_boundary_interval(
            &mut ordinary,
            drives(),
            3,
            NsfEvolvingMannIntervalBudget {
                maximum_cell_steps: 3 * 729,
                maximum_ou_mode_transitions: 3 * 4 * 16,
                maximum_mode_samples: 4 * 9 * 16,
            },
            || false,
        )
        .unwrap();
    ordinary_trial.accept();
    let mut averaged = build_grid(1, false);
    let mut inlet = source(&averaged, 2, true);
    let mut mean = MannStepAverage::prepare(inlet.field(), 16 * 104, || false).unwrap();
    inlet
        .preview_time_averaged_boundary_interval(
            &mut averaged,
            &mut mean,
            drives(),
            3,
            budget(3, 0),
            || false,
        )
        .unwrap()
        .accept();
    assert_eq!(
        inlet.field().checkpoint(),
        ordinary_inlet.field().checkpoint()
    );
    assert_ne!(averaged.f, ordinary.f);
    assert_ne!(
        averaged.cumulative_ledger().open_energy_j.to_bits(),
        ordinary.cumulative_ledger().open_energy_j.to_bits()
    );
    assert!(averaged.cumulative_ledger().residual_rel < 1e-10);
}
