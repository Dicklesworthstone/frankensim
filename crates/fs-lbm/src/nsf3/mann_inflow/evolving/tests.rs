use super::super::tests::{build_grid, config as frozen_config, snapshot};
use super::*;
use crate::d3q27::D3q27;

fn config() -> NsfEvolvingMannConfig {
    let c = frozen_config();
    NsfEvolvingMannConfig {
        pressure_pa: c.pressure_pa,
        temperature_k: c.temperature_k,
        mean_velocity_m_s: c.mean_velocity_m_s,
        grid_origin_in_field_m: c.grid_origin_in_field_m,
        field_to_grid: c.field_to_grid,
        maximum_mode_samples: c.maximum_mode_samples,
    }
}

fn field(grid: &NsfGrid3<D3q27>) -> FixedStepMannField {
    FixedStepMannField::from_snapshot(snapshot(2.0), 531, grid.law.dt_s, usize::MAX).unwrap()
}

fn inlet(grid: &NsfGrid3<D3q27>) -> NsfEvolvingMannInflow {
    let field = field(grid);
    let bytes = NsfEvolvingMannInflow::required_workspace_bytes(&field, 9).unwrap();
    NsfEvolvingMannInflow::new(field, config(), grid, 0, bytes + 16 * 8).unwrap()
}

fn check_profile(grid: &NsfGrid3<D3q27>, inlet: &NsfEvolvingMannInflow) {
    let c = config();
    let dx = grid.law.dx_m;
    for tile in grid.reservoir_aperture_tile_states(0).unwrap() {
        // Independently apply the fixture's 90-degree rotation and ghost
        // position. Atmosphere sampling supplies all phase advection.
        let x = (tile.cell[0] as f64 + 0.5) * dx;
        let z = (tile.cell[2] as f64 + 0.5) * dx;
        let at = [
            c.grid_origin_in_field_m[0] - 0.5 * dx,
            c.grid_origin_in_field_m[1] - x,
            c.grid_origin_in_field_m[2] + z,
        ];
        let u = match inlet.plane_wall() {
            Some(wall) => inlet.field().sample_wall_blocked(at, wall),
            None => match inlet.spatial_filter() {
                Some(filter) => inlet.field().sample_filtered(at, filter),
                None => inlet.field().sample(at),
            },
        }
        .unwrap()
        .u;
        let expected = [-u[1], c.mean_velocity_m_s[1] + u[0], u[2]];
        for (actual, expected) in tile.state.velocity_m_s.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12);
        }
        assert_eq!(tile.state.temperature_k, c.temperature_k);
    }
}

#[test]
fn g1_evolving_mann_uses_bound_cfd_epoch_and_single_advection_with_spatial_averages() {
    for average in 0..3 {
        let mut grid = build_grid(1, false);
        for _ in 0..2 {
            grid.step().unwrap();
        }
        let mut field = field(&grid);
        field.advance_to(4, 4 * 4 * 16, || false).unwrap();
        let bytes = NsfEvolvingMannInflow::required_workspace_bytes(&field, 9).unwrap();
        let mut inlet =
            NsfEvolvingMannInflow::new(field, config(), &grid, 0, bytes + 16 * 8).unwrap();
        inlet = match average {
            1 => inlet.with_face_average(16 * 8, || false).unwrap(),
            2 => inlet.with_cell_average(16 * 8, || false).unwrap(),
            _ => inlet,
        };
        let allocation = inlet.profile.as_ptr();
        let workspace = inlet.workspace_bytes();
        let applied = inlet.apply(&mut grid, 0).unwrap();
        assert_eq!(applied.sampled_time_s, 2.0 * grid.law.dt_s);
        assert_eq!(applied.atmosphere_step, 4);
        assert_eq!(applied.atmosphere_elapsed_time_s, 4.0 * grid.law.dt_s);
        assert_eq!(
            applied.atmosphere_time_s,
            inlet.field().origin_time_s() + 4.0 * grid.law.dt_s
        );
        assert_eq!(applied.correlation_time_range_s, [1.0; 2]);
        assert_eq!(applied.ou_mode_transitions, 0);
        assert_eq!(applied.mode_samples, 9 * 16);
        check_profile(&grid, &inlet);
        let first = grid
            .reservoir_aperture_tile_states(0)
            .unwrap()
            .collect::<Vec<_>>();
        let ledger = grid.step().unwrap();
        assert!(ledger.residual_rel < 1e-10);
        let advanced = inlet.apply(&mut grid, 4 * 16).unwrap();
        assert_eq!(advanced.sampled_time_s, 3.0 * grid.law.dt_s);
        assert_eq!(advanced.atmosphere_step, 5);
        assert_eq!(advanced.ou_mode_transitions, 4 * 16);
        assert_ne!(
            grid.reservoir_aperture_tile_states(0)
                .unwrap()
                .collect::<Vec<_>>(),
            first
        );
        check_profile(&grid, &inlet);
        let accepted = inlet.field().clone();
        let hash = grid.state_hash();
        assert_eq!(inlet.apply(&mut grid, 0).unwrap().ou_mode_transitions, 0);
        assert_eq!(inlet.field(), &accepted);
        assert_eq!(grid.state_hash(), hash);
        assert_eq!(inlet.profile.as_ptr(), allocation);
        assert_eq!(inlet.workspace_bytes(), workspace);
    }
}

#[test]
fn g0_evolving_mann_budgets_clocks_geometry_and_boundary_refuse_atomically() {
    let mut grid = build_grid(1, false);
    let source = field(&grid);
    let bytes = NsfEvolvingMannInflow::required_workspace_bytes(&source, 9).unwrap();
    assert!(matches!(
        NsfEvolvingMannInflow::new(source.clone(), config(), &grid, 0, bytes - 1),
        Err(NsfEvolvingMannError::WorkspaceBytes { .. })
    ));
    let too_small = NsfEvolvingMannConfig {
        maximum_mode_samples: 9 * 16 - 1,
        ..config()
    };
    assert!(matches!(
        NsfEvolvingMannInflow::new(source, too_small, &grid, 0, bytes),
        Err(NsfEvolvingMannError::Work { .. })
    ));
    let other_dt =
        FixedStepMannField::from_snapshot(snapshot(2.0), 531, 2.0 * grid.law.dt_s, usize::MAX)
            .unwrap();
    assert!(matches!(
        NsfEvolvingMannInflow::new(other_dt, config(), &grid, 0, bytes),
        Err(NsfEvolvingMannError::Input(_))
    ));
    let unfiltered = NsfEvolvingMannInflow::new(field(&grid), config(), &grid, 0, bytes).unwrap();
    assert!(unfiltered.with_cell_average(16 * 8, || false).is_err());

    let mut inlet = inlet(&grid);
    inlet.apply(&mut grid, 0).unwrap();
    grid.step().unwrap();
    let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
    let accepted = inlet.field().clone();
    assert!(matches!(
        inlet.apply(&mut grid, 4 * 16 - 1),
        Err(NsfEvolvingMannError::Atmosphere(error)) if error.code == "mann-fixed-step-work"
    ));
    assert_eq!(inlet.field(), &accepted);
    assert_eq!(
        (grid.state_hash(), grid.force.clone(), grid.div.clone()),
        before
    );

    inlet.config.mean_velocity_m_s[1] = 1e6;
    assert!(matches!(
        inlet.apply(&mut grid, 4 * 16),
        Err(NsfEvolvingMannError::Boundary(_))
    ));
    assert_eq!(inlet.field(), &accepted);
    assert_eq!(
        (grid.state_hash(), grid.force.clone(), grid.div.clone()),
        before
    );
    inlet.config = config();
    inlet.apply(&mut grid, 4 * 16).unwrap();

    let accepted = inlet.field().clone();
    let original_dt = grid.law.dt_s;
    grid.law.dt_s *= 2.0;
    let hash = grid.state_hash();
    assert!(matches!(
        inlet.apply(&mut grid, 0),
        Err(NsfEvolvingMannError::Input(_))
    ));
    assert_eq!(grid.state_hash(), hash);
    assert_eq!(inlet.field(), &accepted);
    grid.law.dt_s = original_dt;
    grid.steps -= 1;
    let hash = grid.state_hash();
    assert!(matches!(
        inlet.apply(&mut grid, 0),
        Err(NsfEvolvingMannError::Input(_))
    ));
    assert_eq!(grid.state_hash(), hash);
    assert_eq!(inlet.field(), &accepted);
}

#[test]
fn g4_evolving_mann_cancellation_keeps_boundary_and_noise_for_exact_retry() {
    let make = || {
        let mut grid = build_grid(1, false);
        let mut inlet = inlet(&grid).with_face_average(16 * 8, || false).unwrap();
        inlet.apply(&mut grid, 0).unwrap();
        grid.step().unwrap();
        (grid, inlet)
    };
    let (mut expected_grid, mut expected_inlet) = make();
    let mut count = 0;
    let expected_report = expected_inlet
        .apply_cancellable(&mut expected_grid, 4 * 16, || {
            count += 1;
            false
        })
        .unwrap();
    assert!(count > 12);
    for stop in [1, 4, count - 2, count] {
        let (mut grid, mut inlet) = make();
        let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
        let accepted = inlet.field().clone();
        let allocation = inlet.profile.as_ptr();
        let mut polls = 0;
        assert!(matches!(
            inlet.apply_cancellable(&mut grid, 4 * 16, || {
                polls += 1;
                polls == stop
            }),
            Err(NsfEvolvingMannError::Cancelled)
        ));
        assert_eq!(inlet.field(), &accepted);
        assert_eq!(
            (grid.state_hash(), grid.force.clone(), grid.div.clone()),
            before
        );
        assert_eq!(inlet.apply(&mut grid, 4 * 16).unwrap(), expected_report);
        assert_eq!(inlet.field(), expected_inlet.field());
        assert_eq!(grid.state_hash(), expected_grid.state_hash());
        assert_eq!(inlet.profile.as_ptr(), allocation);
    }
}

#[test]
fn g4_plane_blocked_evolving_mann_keeps_accepted_noise_and_profile_on_refusal() {
    let make = || {
        let mut grid = build_grid(1, false);
        let wall = MannPlaneWall::new(config().grid_origin_in_field_m[2] + 3.5 * grid.law.dx_m)
            .unwrap();
        let mut source = inlet(&grid).with_plane_wall(wall).unwrap();
        source.apply(&mut grid, 0).unwrap();
        check_profile(&grid, &source);
        for tile in grid.reservoir_aperture_tile_states(0).unwrap() {
            if tile.cell[2] == 3 {
                assert_eq!(tile.state.velocity_m_s[2], 0.0);
            }
        }
        grid.step().unwrap();
        (grid, source, wall)
    };
    let (mut expected_grid, mut expected_inlet, _) = make();
    let mut count = 0;
    let expected_report = expected_inlet
        .apply_cancellable(&mut expected_grid, 4 * 16, || {
            count += 1;
            false
        })
        .unwrap();
    assert_eq!(expected_report.ou_mode_transitions, 4 * 16);
    check_profile(&expected_grid, &expected_inlet);
    for outside in [false, true] {
        let (mut grid, mut source, wall) = make();
        let accepted = source.field().checkpoint();
        let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
        let profile = grid.reservoir_aperture_tile_states(0).unwrap().collect::<Vec<_>>();
        let workspace = source.workspace_bytes();
        let allocation = source.profile.as_ptr();
        if outside {
            let outside_wall = MannPlaneWall::new(wall.height_m() + 2.0 * grid.law.dx_m).unwrap();
            source = source.with_plane_wall(outside_wall).unwrap();
            assert!(matches!(
                source.apply(&mut grid, 4 * 16),
                Err(NsfEvolvingMannError::Atmosphere(_))
            ));
            source = source.with_plane_wall(wall).unwrap();
        } else {
            let mut polls = 0;
            assert!(matches!(
                source.apply_cancellable(&mut grid, 4 * 16, || {
                    polls += 1;
                    polls == count - 2
                }),
                Err(NsfEvolvingMannError::Cancelled)
            ));
        }
        assert_eq!(source.field().checkpoint(), accepted);
        assert_eq!((grid.state_hash(), grid.force.clone(), grid.div.clone()), before);
        assert_eq!(grid.reservoir_aperture_tile_states(0).unwrap().collect::<Vec<_>>(), profile);
        assert_eq!(source.apply(&mut grid, 4 * 16).unwrap(), expected_report);
        assert_eq!(source.field(), expected_inlet.field());
        assert_eq!(grid.state_hash(), expected_grid.state_hash());
        assert_eq!(source.profile.as_ptr(), allocation);
        assert_eq!(source.workspace_bytes(), workspace);
    }
}

#[test]
fn g0_plane_blocked_evolving_mann_refuses_averages_in_both_orders_and_normal_mean() {
    let grid = build_grid(1, false);
    let wall = MannPlaneWall::new(0.0).unwrap();
    for face in [false, true] {
        let averaged = if face {
            inlet(&grid).with_face_average(16 * 8, || false)
        } else {
            inlet(&grid).with_cell_average(16 * 8, || false)
        }.unwrap();
        assert!(matches!(averaged.with_plane_wall(wall), Err(NsfEvolvingMannError::Input(_))));
        let blocked = inlet(&grid).with_plane_wall(wall).unwrap();
        let result = if face {
            blocked.with_face_average(16 * 8, || false)
        } else {
            blocked.with_cell_average(16 * 8, || false)
        };
        assert!(matches!(result, Err(NsfEvolvingMannError::Input(_))));
    }
    let mut source = inlet(&grid);
    source.config.mean_velocity_m_s[2] = 0.01;
    assert!(matches!(source.with_plane_wall(wall), Err(NsfEvolvingMannError::Input(_))));
}
