use super::super::tests::{build_grid, config as frozen_config, snapshot};
use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::NsfBoundaryDrive;
use core::cell::Cell;

fn source(grid: &NsfGrid3<D3q27>, average: usize, wall: bool) -> NsfEvolvingMannInflow {
    let c = frozen_config();
    // Deliberately resolve rapid decorrelation: a held interval must produce
    // a different physical flow, not just a different stored source clock.
    let snapshot = HomogeneousMannField::at_anchor(
        snapshot(2.0).model().clone(),
        79,
        vec![2.5 * grid.law.dt_s; 16],
    )
    .unwrap();
    let field =
        FixedStepMannField::from_snapshot(snapshot, 531, grid.law.dt_s, usize::MAX).unwrap();
    let bytes = NsfEvolvingMannInflow::required_workspace_bytes(&field, 9).unwrap();
    let config = NsfEvolvingMannConfig {
        pressure_pa: c.pressure_pa,
        temperature_k: c.temperature_k,
        mean_velocity_m_s: c.mean_velocity_m_s,
        grid_origin_in_field_m: c.grid_origin_in_field_m,
        field_to_grid: c.field_to_grid,
        maximum_mode_samples: c.maximum_mode_samples,
    };
    let mut source = NsfEvolvingMannInflow::new(field, config, grid, 0, bytes + 16 * 32).unwrap();
    if wall {
        let height = (c.grid_origin_in_field_m[2] + 3.0 * grid.law.dx_m).next_down();
        source = source
            .with_plane_wall(MannPlaneWall::new(height).unwrap())
            .unwrap();
    }
    match average {
        1 => source.with_face_average(16 * 32, || false).unwrap(),
        2 => source.with_cell_average(16 * 32, || false).unwrap(),
        _ => source,
    }
}

fn budget(steps: usize, catch_up: usize) -> NsfEvolvingMannIntervalBudget {
    NsfEvolvingMannIntervalBudget {
        maximum_cell_steps: 729 * steps,
        maximum_ou_mode_transitions: 4 * 16 * (steps + catch_up),
        maximum_mode_samples: 9 * 16 * (steps + 1),
    }
}

fn drives() -> [NsfBoundaryDrive; 2] {
    [0.02, -0.01].map(NsfBoundaryDrive::IncomingCharacteristic)
}

fn assert_same_fluid(a: &NsfGrid3<D3q27>, b: &NsfGrid3<D3q27>) {
    assert_eq!(a.state_hash(), b.state_hash());
    assert_eq!(a.steps(), b.steps());
    assert_eq!(a.f, b.f);
    assert_eq!(a.e, b.e);
    assert_eq!(a.force, b.force);
    assert_eq!(a.div, b.div);
    assert_eq!(a.cumulative_ledger(), b.cumulative_ledger());
    assert_eq!(a.last_reservoir_flows(), b.last_reservoir_flows());
    assert_eq!(
        a.last_reservoir_exchange_at(0),
        b.last_reservoir_exchange_at(0)
    );
}

#[test]
fn g1_evolving_interval_matches_native_substeps_and_mass_energy_work() {
    for wall in [false, true] {
        for average in 0..3 {
            let mut direct = build_grid(1, false);
            let mut direct_source = source(&direct, average, wall);
            let mass_before = direct.total_mass_kg();
            direct.set_boundary_drives(drives()).unwrap();
            let mut reservoir_steps = Vec::new();
            let mut face_steps = Vec::new();
            for _ in 0..3 {
                direct_source.apply(&mut direct, 4 * 16).unwrap();
                assert!(direct.step().unwrap().residual_rel < 1e-10);
                reservoir_steps.push(*direct.last_reservoir_exchange_at(0).unwrap());
                face_steps.push(*direct.last_boundary_exchange().unwrap());
            }
            direct_source.apply(&mut direct, 4 * 16).unwrap();

            let mut grid = build_grid(1, false);
            let mut inlet = source(&grid, average, wall);
            let allocation = inlet.profile.as_ptr();
            let workspace = inlet.workspace_bytes();
            let trial = inlet
                .preview_boundary_interval(&mut grid, drives(), 3, budget(3, 0), || false)
                .unwrap();
            assert_same_fluid(trial.fluid().grid(), &direct);
            assert_eq!(
                trial.field().checkpoint(),
                direct_source.field().checkpoint()
            );
            let report = trial.report();
            assert_eq!(report.steps, 3);
            assert_eq!(report.profile_updates, 4);
            assert_eq!(report.mode_samples, 4 * 9 * 16);
            assert_eq!(report.ou_mode_transitions, 3 * 4 * 16);
            assert_eq!(report.first.atmosphere_step, 0);
            assert_eq!(report.endpoint.atmosphere_step, 3);
            assert_eq!(report.endpoint.sampled_time_s, 3.0 * direct.law.dt_s);
            let exchange = trial.fluid().reservoir_exchange_at(0).unwrap();
            for (actual, expected) in [
                (
                    exchange.mass_out_kg,
                    reservoir_steps.iter().map(|x| x.mass_out_kg).sum::<f64>(),
                ),
                (
                    exchange.total_energy_out_j,
                    reservoir_steps.iter().map(|x| x.total_energy_out_j).sum(),
                ),
                (
                    exchange.pressure_work_out_j,
                    reservoir_steps.iter().map(|x| x.pressure_work_out_j).sum(),
                ),
                (
                    exchange.acoustic_pressure_work_out_j,
                    reservoir_steps
                        .iter()
                        .map(|x| x.acoustic_pressure_work_out_j)
                        .sum(),
                ),
                (
                    exchange.mass_exchange_variation_kg,
                    reservoir_steps.iter().map(|x| x.mass_out_kg.abs()).sum(),
                ),
                (
                    exchange.energy_exchange_variation_j,
                    reservoir_steps
                        .iter()
                        .map(|x| x.total_energy_out_j.abs())
                        .sum(),
                ),
                (
                    exchange.mean_acoustic_pressure_pa,
                    reservoir_steps
                        .iter()
                        .map(|x| x.mean_acoustic_pressure_pa)
                        .sum::<f64>()
                        / 3.0,
                ),
                (
                    exchange.normal_volume_flow_m3_s,
                    reservoir_steps
                        .iter()
                        .map(|x| x.normal_volume_flow_m3_s)
                        .sum::<f64>()
                        / 3.0,
                ),
            ] {
                assert_eq!(actual, expected);
            }
            assert_eq!(
                exchange.projection_covariance_j,
                exchange.acoustic_pressure_work_out_j
                    - exchange.mean_acoustic_pressure_pa
                        * exchange.normal_volume_flow_m3_s
                        * (3.0 * direct.law.dt_s)
            );
            for side in 0..2 {
                assert_eq!(
                    trial.fluid().boundary_exchange()[side].mass_out_kg,
                    face_steps
                        .iter()
                        .map(|step| step[side].mass_out_kg)
                        .sum::<f64>()
                );
                assert_eq!(
                    trial.fluid().boundary_exchange()[side].total_energy_out_j,
                    face_steps
                        .iter()
                        .map(|step| step[side].total_energy_out_j)
                        .sum::<f64>()
                );
            }
            let mass_out = trial.fluid().reservoir_mass_out_kg()
                + trial
                    .fluid()
                    .boundary_exchange()
                    .iter()
                    .map(|x| x.mass_out_kg)
                    .sum::<f64>();
            assert!((direct.total_mass_kg() - mass_before + mass_out).abs() < 1e-10 * mass_before);
            assert_eq!(
                trial.fluid().energy_ledger().open_energy_j,
                direct.cumulative_ledger().open_energy_j
            );
            assert_eq!(
                trial.fluid().energy_ledger().open_boundary_energy_j,
                direct.cumulative_ledger().open_boundary_energy_j
            );
            assert_eq!(
                trial.fluid().energy_ledger().internal_j,
                direct.cumulative_ledger().internal_j
            );
            assert_eq!(
                trial.fluid().energy_ledger().kinetic_j,
                direct.cumulative_ledger().kinetic_j
            );
            assert!(trial.fluid().energy_ledger().residual_rel < 1e-10);
            trial.accept();
            assert_same_fluid(&grid, &direct);
            assert_eq!(
                inlet.field().checkpoint(),
                direct_source.field().checkpoint()
            );
            assert_eq!(inlet.profile.as_ptr(), allocation);
            assert_eq!(inlet.workspace_bytes(), workspace);
            let accepted_hash = grid.state_hash();
            assert_eq!(inlet.apply(&mut grid, 0).unwrap().ou_mode_transitions, 0);
            assert_eq!(grid.state_hash(), accepted_hash);
        }
    }
}

#[test]
fn g1_substep_decorrelation_changes_resolved_flow_and_native_energy_exchange() {
    let mut held = build_grid(1, false);
    let mut held_source = source(&held, 2, true);
    held.set_boundary_drives(drives()).unwrap();
    held_source.apply(&mut held, 0).unwrap();
    let held_trial = held
        .preview_boundary_interval(drives(), 4, 4 * 729, || false)
        .unwrap();
    let held_energy = held_trial.energy_ledger().open_energy_j;
    held_trial.accept();

    let mut evolving = build_grid(1, false);
    let mut evolving_source = source(&evolving, 2, true);
    let trial = evolving_source
        .preview_boundary_interval(&mut evolving, drives(), 4, budget(4, 0), || false)
        .unwrap();
    assert_ne!(
        trial.fluid().energy_ledger().open_energy_j.to_bits(),
        held_energy.to_bits()
    );
    let velocity_difference: f64 = (3..=5)
        .flat_map(|x| (3..=5).map(move |z| [x, 0, z]))
        .map(|cell| {
            let a = trial.fluid().grid().cell(cell).1;
            let b = held.cell(cell).1;
            (0..3).map(|axis| (a[axis] - b[axis]).abs()).sum::<f64>()
        })
        .sum();
    assert!(
        velocity_difference > 1e-8,
        "resolved velocity difference {velocity_difference}"
    );
    assert!(trial.fluid().energy_ledger().residual_rel < 1e-10);
}

#[test]
fn g4_interval_budgets_drop_and_cancellation_preserve_both_owners_for_retry() {
    let make = || {
        let mut grid = build_grid(1, false);
        let inlet = source(&grid, 2, true);
        // The bound atmosphere intentionally lags this original fluid clock.
        grid.step().unwrap();
        grid.step().unwrap();
        (grid, inlet)
    };
    let (mut reference, mut reference_source) = make();
    let polls = Cell::new(0);
    let trial = reference_source
        .preview_boundary_interval(&mut reference, drives(), 2, budget(2, 2), || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let expected_report = trial.report();
    assert_eq!(expected_report.first.ou_mode_transitions, 2 * 4 * 16);
    assert_eq!(expected_report.ou_mode_transitions, 4 * 4 * 16);
    trial.accept();

    for refusal in 0..7 {
        let (mut grid, mut inlet) = make();
        let (before_grid, _) = make();
        let before_source = inlet.field().checkpoint();
        let allocation = inlet.profile.as_ptr();
        let mut limits = budget(2, 2);
        match refusal {
            0 => limits.maximum_cell_steps -= 1,
            1 => limits.maximum_ou_mode_transitions -= 1,
            2 => limits.maximum_mode_samples -= 1,
            _ => (),
        }
        let current = Cell::new(0);
        let stop = match refusal {
            4 => 1,
            5 => polls.get() / 2,
            6 => polls.get(),
            _ => 0,
        };
        let result = inlet.preview_boundary_interval(&mut grid, drives(), 2, limits, || {
            current.set(current.get() + 1);
            stop != 0 && current.get() == stop
        });
        if refusal == 3 {
            let candidate = result.unwrap();
            assert_eq!(candidate.report(), expected_report);
            assert_eq!(
                candidate.field().checkpoint(),
                reference_source.field().checkpoint()
            );
            drop(candidate);
        } else if refusal >= 4 {
            assert!(matches!(result, Err(NsfEvolvingMannError::Cancelled)));
            drop(result);
        } else {
            assert!(result.is_err());
            drop(result);
        }
        assert_same_fluid(&grid, &before_grid);
        assert_eq!(inlet.field().checkpoint(), before_source);
        let retry = inlet
            .preview_boundary_interval(&mut grid, drives(), 2, budget(2, 2), || false)
            .unwrap();
        assert_eq!(retry.report(), expected_report);
        retry.accept();
        assert_same_fluid(&grid, &reference);
        assert_eq!(
            inlet.field().checkpoint(),
            reference_source.field().checkpoint()
        );
        assert_eq!(inlet.profile.as_ptr(), allocation);
    }
}

#[test]
fn g4_native_fluid_refusal_after_admitted_profile_publishes_neither_owner() {
    let make = || {
        let mut grid = build_grid(1, false);
        let mut inlet = source(&grid, 0, false);
        inlet.config.temperature_k = grid.spec.temperature_k;
        // Both the initial fluid and inlet satisfy this bound. Compression
        // during an original native step must be admitted independently.
        grid.spec.maximum_temperature_ratio = 1e-12;
        (grid, inlet)
    };
    let (mut direct, mut direct_source) = make();
    direct_source.apply(&mut direct, 0).unwrap();
    assert!(matches!(
        direct.step(),
        Err(NsfError::TemperatureRatio { .. })
    ));
    let (mut grid, mut inlet) = make();
    let (before_grid, _) = make();
    let before_source = inlet.field().checkpoint();
    let held_drives = grid.boundary_drives().unwrap();
    assert!(matches!(
        inlet.preview_boundary_interval(&mut grid, held_drives, 2, budget(2, 0), || false,),
        Err(NsfEvolvingMannError::Boundary(
            NsfError::TemperatureRatio { .. }
        ))
    ));
    assert_same_fluid(&grid, &before_grid);
    assert_eq!(inlet.field().checkpoint(), before_source);
}
