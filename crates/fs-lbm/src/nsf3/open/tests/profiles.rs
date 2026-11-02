//! Spatial boundary tests use the original source and original flux owners.
use super::*;

fn profiled_states(grid: &NsfGrid3<D3q27>) -> Vec<NsfReservoirTileState> {
    grid.reservoir_aperture_tile_states(0)
        .unwrap()
        .enumerate()
        .map(|(n, mut tile)| {
            let temperature = 293.15 + 0.05 * n as f64;
            tile.state = grid
                .reservoir_aperture_state_for(0, 101325.0 + 2.0 * (n as f64 - 4.0), temperature)
                .unwrap();
            tile.state.velocity_m_s = [
                0.01 * (n as f64 - 4.0),
                0.2 + 0.015 * n as f64,
                -0.006 * n as f64,
            ];
            tile
        })
        .collect()
}

#[test]
fn g0_reservoir_profile_uniform_limit_is_bit_identical_to_original_control() {
    let build = || {
        aperture_grid(0, 1)
            .with_reservoir_aperture(aperture(1, false, 293.35))
            .unwrap()
    };
    let mut scalar = build();
    let mut profile = build();
    for _ in 0..2 {
        scalar.step().unwrap();
        profile.step().unwrap();
    }
    let mut prescribed = scalar
        .reservoir_aperture_state_at(101330.0, 293.25)
        .unwrap();
    prescribed.velocity_m_s = [0.03, 0.1, -0.02];
    let mut tiles: Vec<_> = profile
        .reservoir_aperture_tile_states(0)
        .unwrap()
        .map(|tile| NsfReservoirTileState {
            state: prescribed,
            ..tile
        })
        .collect();
    tiles.reverse();
    scalar.set_reservoir_aperture_state(prescribed).unwrap();
    profile.set_reservoir_aperture_profile(0, &tiles).unwrap();
    assert!(profile.reservoir_apertures[0].profile.is_none());
    assert_eq!(scalar.state_hash(), profile.state_hash());
    assert_eq!(scalar.force, profile.force);
    assert_eq!(scalar.div, profile.div);
    for _ in 0..3 {
        assert_eq!(scalar.step().unwrap(), profile.step().unwrap());
        assert_eq!(scalar.state_hash(), profile.state_hash());
        assert_eq!(
            scalar.last_reservoir_flows(),
            profile.last_reservoir_flows()
        );
        assert_eq!(
            scalar.last_reservoir_exchange(),
            profile.last_reservoir_exchange()
        );
    }
}

#[test]
fn g1_reservoir_profile_uses_original_energy_and_spatial_pressure_work() {
    let mut grid = aperture_grid(0, 1)
        .with_reservoir_aperture(aperture(1, false, 293.15))
        .unwrap();
    let tiles = profiled_states(&grid);
    grid.set_reservoir_aperture_profile(0, &tiles).unwrap();
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
    let (_, _, flux) = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &velocity,
            &grid.theta,
            &vec![[0.0; 6]; grid.n],
            &[],
            grid.open_faces.as_ref(),
        )
        .unwrap();
    let flows =
        grid.reservoir_volume_flows(&rho, &velocity, &grid.theta, &rho, &velocity, &grid.theta);
    let area = grid.law.dx_m.powi(2);
    let mut expected_power = 0.0;
    let mut expected_acoustic_power = 0.0;
    let mut ghost_energies = Vec::new();
    for (tile, flow) in tiles.iter().zip(flows) {
        let i = grid.index(tile.cell);
        let ghost = grid.open_at(i, 1, -1).unwrap();
        let original = grid.admit_open_state(tile.state).unwrap();
        assert_eq!(ghost.energy, original.energy);
        assert_eq!(ghost.u, original.u);
        assert_eq!(ghost.theta, original.theta);
        ghost_energies.push(ghost.energy);
        let pressure = grid.law.pressure_pa(ghost.rho, ghost.theta);
        let pressure_i = grid.law.pressure_pa(rho[i], grid.theta[i]);
        let speed_i = velocity[i][1] * grid.law.dx_m / grid.law.dt_s;
        let speed_g = tile.state.velocity_m_s[1];
        expected_power -= 0.5 * area * (pressure_i * speed_i + pressure * speed_g);
        expected_acoustic_power -=
            0.5 * area * ((pressure_i - 101325.0) * speed_i + (pressure - 101325.0) * speed_g);
        assert!((flow.outward_volume_flow_m3_s + 0.5 * area * (speed_i + speed_g)).abs() < 1e-23);
        assert!((flow.mean_pressure_pa - 0.5 * (pressure_i + pressure)).abs() < 1e-9);
    }
    assert!(ghost_energies.windows(2).any(|pair| pair[0] != pair[1]));
    assert!(
        (flux.reservoir[0].pressure_power_w - expected_power).abs() < 2e-14 * expected_power.abs()
    );
    assert!(
        (flux.reservoir[0].acoustic_pressure_power_w - expected_acoustic_power).abs()
            < 2e-11 * expected_acoustic_power.abs()
    );
    assert!(
        grid.reservoir_link_is_open([4, 0, 4], [1, -1, 1], 1, 1, 0),
        "different tile states must not create artificial walls inside one aperture"
    );
}

#[test]
fn g3_reservoir_profile_replays_conservative_steps_and_canonical_order() {
    let build = |workers| {
        aperture_grid(0, workers)
            .with_reservoir_aperture(aperture(1, false, 293.15))
            .unwrap()
    };
    let mut grid = build(1);
    let mut replay = build(3);
    let mut tiles = profiled_states(&grid);
    grid.set_reservoir_aperture_profile(0, &tiles).unwrap();
    tiles.reverse();
    replay.set_reservoir_aperture_profile(0, &tiles).unwrap();
    assert_eq!(grid.state_hash(), replay.state_hash());
    assert_eq!(
        grid.reservoir_aperture_tile_states(0)
            .unwrap()
            .collect::<Vec<_>>(),
        replay
            .reservoir_aperture_tile_states(0)
            .unwrap()
            .collect::<Vec<_>>()
    );
    for _ in 0..4 {
        let before = (grid.total_mass_kg(), grid.total_energy_j());
        let ledger = grid.step().unwrap();
        assert_eq!(ledger, replay.step().unwrap());
        assert_eq!(grid.state_hash(), replay.state_hash());
        assert_eq!(
            grid.last_reservoir_exchange(),
            replay.last_reservoir_exchange()
        );
        let mass_out = grid.last_reservoir_mass_out_kg().unwrap()
            + grid
                .last_boundary_exchange()
                .unwrap()
                .iter()
                .map(|face| face.mass_out_kg)
                .sum::<f64>();
        assert!((grid.total_mass_kg() - before.0 + mass_out).abs() < 5e-12 * before.0);
        assert!(
            (grid.total_energy_j()
                - before.1
                - ledger.wall_heat_j
                - ledger.open_boundary_energy_j
                - ledger.open_energy_j)
                .abs()
                < 5e-12 * before.1
        );
    }
}

#[test]
fn g4_reservoir_profile_refusals_preserve_warm_state_and_scalar_control_clears_it() {
    let mut grid = aperture_grid(0, 1)
        .with_reservoir_aperture(aperture(1, false, 293.15))
        .unwrap();
    grid.step().unwrap();
    let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
    let physical = (
        grid.total_energy_j(),
        grid.steps(),
        grid.cumulative_ledger(),
        grid.last_reservoir_mass_out_kg(),
    );
    let tiles = profiled_states(&grid);
    grid.set_reservoir_aperture_profile(0, &tiles).unwrap();
    for (i, before) in cells.iter().enumerate() {
        let after = grid.cell(grid.at(i));
        assert!((before.0 - after.0).abs() < 2e-14);
        assert_eq!(before.2, after.2);
        for a in 0..3 {
            assert!((before.1[a] - after.1[a]).abs() < 2e-14);
        }
    }
    assert_eq!(
        (
            grid.total_energy_j(),
            grid.steps(),
            grid.cumulative_ledger(),
            grid.last_reservoir_mass_out_kg()
        ),
        physical
    );
    let snapshot = (
        grid.state_hash(),
        grid.force.clone(),
        grid.div.clone(),
        grid.reservoir_aperture_tile_states(0)
            .unwrap()
            .collect::<Vec<_>>(),
    );
    let mut duplicate = tiles.clone();
    duplicate[1] = duplicate[0];
    let mut wrong_owner = tiles.clone();
    wrong_owner[0].cell[1] = 8;
    let mut outside = tiles.clone();
    outside[0].cell[0] = usize::MAX;
    let mut invalid = tiles.clone();
    invalid[4].state.velocity_m_s = [1000.0; 3];
    for (owner, states) in [
        (1, tiles.as_slice()),
        (0, &tiles[..8]),
        (0, &duplicate),
        (0, &wrong_owner),
        (0, &outside),
        (0, &invalid),
    ] {
        assert!(grid.set_reservoir_aperture_profile(owner, states).is_err());
        assert_eq!(
            (
                grid.state_hash(),
                grid.force.clone(),
                grid.div.clone(),
                grid.reservoir_aperture_tile_states(0)
                    .unwrap()
                    .collect::<Vec<_>>()
            ),
            snapshot
        );
    }
    let scalar = grid.reservoir_aperture_at(0).unwrap().state;
    grid.set_reservoir_aperture_state(scalar).unwrap();
    assert!(grid.reservoir_apertures[0].profile.is_none());
    assert!(
        grid.reservoir_aperture_tile_states(0)
            .unwrap()
            .all(|tile| tile.state == scalar)
    );
    assert_ne!(grid.state_hash(), snapshot.0);
}
