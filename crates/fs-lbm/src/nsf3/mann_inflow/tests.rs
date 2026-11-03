use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{
    MAX_QUALIFIED_MACH, NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfReservoirAperture,
    NsfSpec, WallThermal,
};
use fs_atmo::mann::MannParams;
use fs_atmo::mann_field::{MannQuadrature, MannSpectralModel};
use fs_material::gas::{GasSpec, GasState};

pub(super) fn build_grid(workers: usize, reverse_tiles: bool) -> NsfGrid3<D3q27> {
    build_grid_with_apertures(workers, reverse_tiles, false)
}

fn build_grid_with_apertures(
    workers: usize,
    reverse_tiles: bool,
    extra_face: bool,
) -> NsfGrid3<D3q27> {
    let spec = NsfSpec {
        dims: [9; 3],
        periodic: [false; 3],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: MAX_QUALIFIED_MACH,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, Phys. Fluids 24 066102, N2 at 293 K",
        maximum_cells: 729,
        wall: WallThermal::Adiabatic,
        workers,
        transport_points: 33,
    };
    let mut cells: Vec<_> = (3..=5)
        .flat_map(|x| (3..=5).map(move |z| [x, 0, z]))
        .collect();
    if reverse_tiles {
        cells.reverse();
    }
    let gas = GasState::try_new(&spec.gas, spec.temperature_k, spec.pressure_pa).unwrap();
    let state = NsfOpenState {
        density_kg_m3: gas.density,
        velocity_m_s: [0.0; 3],
        temperature_k: gas.temperature,
    };
    let mut apertures = vec![NsfReservoirAperture {
        axis: 1,
        high: false,
        cells,
        state,
    }];
    if extra_face {
        apertures.push(NsfReservoirAperture {
            axis: 2,
            high: false,
            cells: vec![[4, 4, 0]],
            state,
        });
    }
    NsfGrid3::new(spec, |_| false, |_| (1.0, [0.0; 3], 293.15))
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 101325.0,
                reference_temperature_k: 293.15,
                reference_velocity_m_s: [0.0; 3],
                incoming_characteristic_m_s: 0.0,
            }; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap()
        .with_reservoir_apertures(apertures)
        .unwrap()
}

pub(super) fn snapshot(max_kl: f64) -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        41,
        MannParams {
            alpha_eps23: 0.32,
            length_m: 2e-4,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: 16,
            min_kl: 0.03,
            max_kl,
        },
        [3.0, -0.2, 0.4],
    )
    .unwrap();
    // Retain a nonzero source phase; do not advance it onto the CFD clock.
    HomogeneousMannField::at_anchor(model, 79, vec![1.0; 16]).unwrap()
}

pub(super) fn config() -> NsfFrozenMannConfig {
    NsfFrozenMannConfig {
        pressure_pa: 101325.0,
        temperature_k: 293.35,
        mean_velocity_m_s: [0.0, 0.4, 0.0],
        convection_velocity_m_s: [0.0, 0.4, 0.0],
        grid_origin_in_field_m: [0.001, -0.0003, 0.0007],
        // Field streamwise x maps to grid y, keeping the authored mean aligned.
        field_to_grid: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        reference_time_s: 1e-8,
        maximum_mode_samples: 9 * 16,
    }
}

#[test]
fn g1_filtered_mann_inlet_uses_rotated_physical_face_and_volume_averages() {
    for face in [true, false] {
        let mut grid = build_grid(1, false);
        let controls = config();
        let inlet = NsfFrozenMannInflow::new(snapshot(10.0), controls).unwrap();
        let inlet = if face {
            inlet.with_face_average(&grid, 0, 16 * 8, || false).unwrap()
        } else {
            inlet.with_cell_average(&grid, 16 * 8, || false).unwrap()
        };
        let dx = grid.law.dx_m;
        let expected_edges = [
            [0.0, -0.5 * dx, 0.0],
            if face { [0.0; 3] } else { [0.5 * dx, 0.0, 0.0] },
            [0.0, 0.0, 0.5 * dx],
        ];
        let filter = inlet.spatial_filter().unwrap();
        assert_eq!(filter.half_edges_m(), expected_edges);
        let retained = inlet.snapshot().clone();
        let applied = inlet.apply(&mut grid, 0).unwrap();
        let mut differs_from_point = false;
        for tile in grid.reservoir_aperture_tile_states(0).unwrap() {
            let x = (tile.cell[0] as f64 + 0.5) * dx;
            let y = -0.5 * dx
                - controls.convection_velocity_m_s[1]
                    * (applied.sampled_time_s - controls.reference_time_s);
            let z = (tile.cell[2] as f64 + 0.5) * dx;
            let at = [
                controls.grid_origin_in_field_m[0] + y,
                controls.grid_origin_in_field_m[1] - x,
                controls.grid_origin_in_field_m[2] + z,
            ];
            let filtered = retained.sample_filtered(at, filter).unwrap().u;
            let point = retained.sample(at).unwrap().u;
            differs_from_point |= filtered
                .iter()
                .zip(point)
                .any(|(a, b)| (a - b).abs() > 1e-9);
            let expected = [
                -filtered[1],
                controls.mean_velocity_m_s[1] + filtered[0],
                filtered[2],
            ];
            for axis in 0..3 {
                assert!((tile.state.velocity_m_s[axis] - expected[axis]).abs() < 1e-12);
            }
        }
        assert!(differs_from_point);
        assert_eq!(inlet.snapshot(), &retained);
        assert_eq!(inlet.apply(&mut grid, 0).unwrap(), applied);
    }
}

#[test]
fn g4_filtered_mann_inlet_refuses_budget_cancellation_spacing_and_face_mismatch() {
    let mut grid = build_grid_with_apertures(1, false, true);
    let make = || NsfFrozenMannInflow::new(snapshot(2.0), config()).unwrap();
    assert!(matches!(
        make().with_cell_average(&grid, 16 * 8 - 1, || false),
        Err(NsfFrozenMannError::Atmosphere(error)) if error.code == "mann-filter-budget"
    ));
    let mut polls = 0;
    assert!(matches!(
        make().with_cell_average(&grid, 16 * 8, || {
            polls += 1;
            polls == 4
        }),
        Err(NsfFrozenMannError::Cancelled)
    ));
    let inlet = make()
        .with_face_average(&grid, 0, 16 * 8, || false)
        .unwrap();
    let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
    assert!(matches!(
        inlet.apply(&mut grid, 1),
        Err(NsfFrozenMannError::Input(_))
    ));
    assert_eq!(
        (grid.state_hash(), grid.force.clone(), grid.div.clone()),
        before
    );
    // Alter only the test fixture's spatial chart to exercise reuse admission;
    // the refusal must precede sampling or material/boundary publication.
    grid.law.dx_m *= 2.0;
    let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
    assert!(matches!(
        inlet.apply(&mut grid, 0),
        Err(NsfFrozenMannError::Input(_))
    ));
    assert_eq!(
        (grid.state_hash(), grid.force.clone(), grid.div.clone()),
        before
    );
}

fn sgs(maximum_eddy_to_molecular: f64) -> crate::nsf3::NsfSgsConfig {
    crate::nsf3::NsfSgsConfig {
        closure: crate::sgs3::SgsConfig {
            model: crate::sgs3::SgsModel::Smagorinsky { cs: 0.17 },
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular,
            usage: crate::sgs3::SgsUse::Exploration,
        },
        turbulent_prandtl: 0.9,
    }
}

#[test]
fn g1_mann_inlet_samples_original_ghost_centres_frame_clock_and_gas_energy() {
    let mut grid = build_grid(1, false);
    let source = snapshot(2.0);
    let retained = source.clone();
    let controls = config();
    let inlet = NsfFrozenMannInflow::new(source, controls).unwrap();
    let gas = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        controls.temperature_k,
        controls.pressure_pa,
    )
    .unwrap();
    let applied = inlet.apply(&mut grid, 0).unwrap();
    assert_eq!(applied.sampled_time_s, 0.0);
    assert_eq!(applied.frozen_tick, retained.tick());
    assert_eq!(applied.tiles, 9);
    assert_eq!(applied.mode_samples, 9 * 16);
    assert!(applied.spatial_phase_bound_rad < core::f64::consts::PI);
    let states: Vec<_> = grid.reservoir_aperture_tile_states(0).unwrap().collect();
    assert!(
        states
            .windows(2)
            .any(|pair| pair[0].state.velocity_m_s != pair[1].state.velocity_m_s)
    );
    for tile in &states {
        // Explicit 90-degree inverse rotation and the low-face ghost at y=-dx/2.
        let x = (tile.cell[0] as f64 + 0.5) * grid.law.dx_m;
        let y = -0.5 * grid.law.dx_m
            - controls.convection_velocity_m_s[1]
                * (applied.sampled_time_s - controls.reference_time_s);
        let z = (tile.cell[2] as f64 + 0.5) * grid.law.dx_m;
        let sample = retained
            .sample([
                controls.grid_origin_in_field_m[0] + y,
                controls.grid_origin_in_field_m[1] - x,
                controls.grid_origin_in_field_m[2] + z,
            ])
            .unwrap();
        let expected = [
            -sample.u[1],
            controls.mean_velocity_m_s[1] + sample.u[0],
            sample.u[2],
        ];
        for a in 0..3 {
            assert!((tile.state.velocity_m_s[a] - expected[a]).abs() < 1e-12);
        }
        assert_eq!(tile.state.density_kg_m3, gas.density);
        assert_eq!(tile.state.temperature_k, gas.temperature);
        let ghost = grid.open_at(grid.index(tile.cell), 1, -1).unwrap();
        let kinetic = 0.5
            * tile
                .state
                .velocity_m_s
                .iter()
                .map(|u| (u * grid.law.dt_s / grid.law.dx_m).powi(2))
                .sum::<f64>();
        let thermal = grid.law.theta_of_kelvin(gas.temperature) / (grid.law.gamma - 1.0);
        assert_eq!(ghost.energy, ghost.rho * (thermal + kinetic));
    }
    assert_eq!(inlet.snapshot(), &retained);
    let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
    assert_eq!(inlet.apply(&mut grid, 0).unwrap(), applied);
    assert_eq!(
        (grid.state_hash(), grid.force.clone(), grid.div.clone()),
        before
    );
}

#[test]
fn g3_mann_inlet_drives_original_conservative_steps_below_the_ou_tick() {
    let mut grid = build_grid(1, false).with_sgs(sgs(100.0)).unwrap();
    let mut replay = build_grid(3, true).with_sgs(sgs(100.0)).unwrap();
    let inlet = NsfFrozenMannInflow::new(snapshot(2.0), config()).unwrap();
    let retained = inlet.snapshot().clone();
    let initial_energy = grid.total_energy_j();
    let mut previous = None;
    for _ in 0..4 {
        let applied = inlet.apply(&mut grid, 0).unwrap();
        let replayed = inlet.apply(&mut replay, 0).unwrap();
        assert_eq!(applied, replayed);
        assert!(
            grid.sgs_transport
                .iter()
                .any(|transport| transport.eddy_viscosity_lu > 0.0
                    && transport.turbulent_conductivity_lu > 0.0)
        );
        assert!(applied.sampled_time_s < 1.0 / fs_atmo::TICK_HZ);
        let profile: Vec<_> = grid.reservoir_aperture_tile_states(0).unwrap().collect();
        if let Some(last) = &previous {
            assert_ne!(
                &profile, last,
                "convection must not freeze until the next 120 Hz OU tick"
            );
        }
        previous = Some(profile);
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
    assert_ne!(grid.total_energy_j(), initial_energy);
    assert!((0..grid.n).any(|i| grid.cell(grid.at(i)).1.iter().any(|u| u.abs() > 1e-8)));
    assert_eq!(inlet.snapshot(), &retained);
}

#[test]
fn g4_mann_inlet_refuses_aliasing_work_material_and_cancellation_atomically() {
    let mut grid = build_grid(1, false);
    let controls = config();
    let inlet = NsfFrozenMannInflow::new(snapshot(2.0), controls).unwrap();
    inlet.apply(&mut grid, 0).unwrap();
    let retained = inlet.snapshot().clone();
    let before = (
        grid.state_hash(),
        grid.total_energy_j(),
        grid.steps(),
        grid.force.clone(),
        grid.div.clone(),
        grid.reservoir_aperture_tile_states(0)
            .unwrap()
            .collect::<Vec<_>>(),
    );

    let spatial = NsfFrozenMannInflow::new(snapshot(80.0), controls).unwrap();
    assert!(matches!(spatial.apply(&mut grid, 0),
        Err(NsfFrozenMannError::SpatialNyquist { phase_bound_rad, maximum_rad })
        if phase_bound_rad > maximum_rad && maximum_rad == core::f64::consts::PI));
    let temporal = NsfFrozenMannInflow::new(
        retained.clone(),
        NsfFrozenMannConfig {
            convection_velocity_m_s: [0.0, 1e6, 0.0],
            ..controls
        },
    )
    .unwrap();
    assert!(matches!(temporal.apply(&mut grid, 0),
        Err(NsfFrozenMannError::TemporalNyquist { phase_bound_rad, maximum_rad })
        if phase_bound_rad > maximum_rad && maximum_rad == core::f64::consts::PI));
    let bounded = NsfFrozenMannInflow::new(
        retained.clone(),
        NsfFrozenMannConfig {
            maximum_mode_samples: 9 * 16 - 1,
            ..controls
        },
    )
    .unwrap();
    assert!(matches!(
        bounded.apply(&mut grid, 0),
        Err(NsfFrozenMannError::Work {
            requested: 144,
            maximum: 143
        })
    ));
    for changed in [
        NsfFrozenMannConfig {
            temperature_k: 400.0,
            ..controls
        },
        NsfFrozenMannConfig {
            mean_velocity_m_s: [1000.0; 3],
            ..controls
        },
    ] {
        let invalid = NsfFrozenMannInflow::new(retained.clone(), changed).unwrap();
        assert!(matches!(
            invalid.apply(&mut grid, 0),
            Err(NsfFrozenMannError::Boundary(_))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        inlet.apply_cancellable(&mut grid, 0, || {
            polls += 1;
            polls == 3
        }),
        Err(NsfFrozenMannError::Cancelled)
    );
    assert!(inlet.apply(&mut grid, 1).is_err());
    assert_eq!(
        (
            grid.state_hash(),
            grid.total_energy_j(),
            grid.steps(),
            grid.force.clone(),
            grid.div.clone(),
            grid.reservoir_aperture_tile_states(0)
                .unwrap()
                .collect::<Vec<_>>()
        ),
        before
    );
    assert_eq!(inlet.snapshot(), &retained);
}

#[test]
fn g0_mann_inlet_requires_an_explicit_proper_frame() {
    let controls = config();
    for rotation in [
        [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, f64::NAN, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ] {
        assert!(
            NsfFrozenMannInflow::new(
                snapshot(2.0),
                NsfFrozenMannConfig {
                    field_to_grid: rotation,
                    ..controls
                }
            )
            .is_err()
        );
    }
}

#[test]
fn g4_mann_inlet_sgs_rebind_refusal_restores_the_previous_nonuniform_owner() {
    let mut grid = build_grid(1, false);
    let previous: Vec<_> = grid
        .reservoir_aperture_tile_states(0)
        .unwrap()
        .enumerate()
        .map(|(n, mut tile)| {
            // A spatial pressure profile with zero physical velocity has no
            // modeled shear, yet exercises restoration of a real old profile.
            tile.state.density_kg_m3 *= 1.0 + 1e-5 * n as f64;
            tile
        })
        .collect();
    grid.set_reservoir_aperture_profile(0, &previous).unwrap();
    let mut grid = grid.with_sgs(sgs(1e-8)).unwrap();
    let receipt = |grid: &NsfGrid3<D3q27>| {
        (
            grid.state_hash(),
            grid.sgs_transport.clone(),
            grid.div.clone(),
            grid.force.clone(),
            grid.boundary_drives(),
            grid.cumulative_ledger(),
            grid.reservoir_aperture_tile_states(0)
                .unwrap()
                .collect::<Vec<_>>(),
        )
    };
    let before = receipt(&grid);
    let inlet = NsfFrozenMannInflow::new(snapshot(2.0), config()).unwrap();
    assert!(matches!(
        inlet.apply(&mut grid, 0),
        Err(NsfFrozenMannError::Boundary(NsfError::Sgs(
            crate::sgs3::SgsError::ViscosityRatioExceeded { .. }
        )))
    ));
    assert_eq!(receipt(&grid), before);
    let scalar = NsfOpenState {
        velocity_m_s: [0.0, 0.4, 0.0],
        ..grid.reservoir_aperture_at(0).unwrap().state
    };
    assert!(matches!(
        grid.set_reservoir_aperture_state(scalar),
        Err(NsfError::Sgs(
            crate::sgs3::SgsError::ViscosityRatioExceeded { .. }
        ))
    ));
    assert_eq!(receipt(&grid), before);
}

#[test]
fn g1_plane_blocked_frozen_mann_uses_the_original_rotated_ghost_profile() {
    let mut grid = build_grid(1, false);
    let controls = config();
    let dx = grid.law.dx_m;
    let wall = MannPlaneWall::new(controls.grid_origin_in_field_m[2] + 3.5 * dx).unwrap();
    let inlet = NsfFrozenMannInflow::new(snapshot(2.0), controls)
        .unwrap()
        .with_plane_wall(wall)
        .unwrap();
    assert_eq!(inlet.plane_wall(), Some(wall));
    let retained = inlet.snapshot().checkpoint();
    for _ in 0..2 {
        let applied = inlet.apply(&mut grid, 0).unwrap();
        for tile in grid.reservoir_aperture_tile_states(0).unwrap() {
            let at = [
                controls.grid_origin_in_field_m[0] - 0.5 * dx
                    - controls.convection_velocity_m_s[1]
                        * (applied.sampled_time_s - controls.reference_time_s),
                controls.grid_origin_in_field_m[1] - (tile.cell[0] as f64 + 0.5) * dx,
                controls.grid_origin_in_field_m[2] + (tile.cell[2] as f64 + 0.5) * dx,
            ];
            let u = inlet.snapshot().sample_wall_blocked(at, wall).unwrap().u;
            let expected = [-u[1], controls.mean_velocity_m_s[1] + u[0], u[2]];
            for i in 0..3 {
                assert!((tile.state.velocity_m_s[i] - expected[i]).abs() < 1e-12);
            }
            if tile.cell[2] == 3 {
                assert_eq!(tile.state.velocity_m_s[2], 0.0);
            }
        }
        assert!(grid.step().unwrap().residual_rel < 1e-10);
    }
    assert_eq!(inlet.snapshot().checkpoint(), retained);
}

#[test]
fn g4_plane_blocked_frozen_mann_refuses_filters_normal_mean_and_outside_samples() {
    let mut grid = build_grid(1, false);
    let wall = MannPlaneWall::new(0.0).unwrap();
    let make = || NsfFrozenMannInflow::new(snapshot(2.0), config()).unwrap();
    for face in [false, true] {
        let averaged = if face {
            make().with_face_average(&grid, 0, 16 * 8, || false)
        } else {
            make().with_cell_average(&grid, 16 * 8, || false)
        }.unwrap();
        assert!(matches!(averaged.with_plane_wall(wall), Err(NsfFrozenMannError::Input(_))));
        let blocked = make().with_plane_wall(wall).unwrap();
        let result = if face {
            blocked.with_face_average(&grid, 0, 16 * 8, || false)
        } else {
            blocked.with_cell_average(&grid, 16 * 8, || false)
        };
        assert!(matches!(result, Err(NsfFrozenMannError::Input(_))));
    }
    for convection in [false, true] {
        let mut controls = config();
        if convection {
            controls.convection_velocity_m_s[2] = 0.01;
        } else {
            controls.mean_velocity_m_s[2] = 0.01;
        }
        assert!(matches!(
            NsfFrozenMannInflow::new(snapshot(2.0), controls).unwrap().with_plane_wall(wall),
            Err(NsfFrozenMannError::Input(_))
        ));
    }
    // Invert field z with a proper rotation so the first two ghost samples
    // succeed and the third is below the plane. Partial sampling must not publish.
    let controls = NsfFrozenMannConfig {
        field_to_grid: [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        ..config()
    };
    let wall = MannPlaneWall::new(controls.grid_origin_in_field_m[2] - 4.5 * grid.law.dx_m)
        .unwrap();
    let inlet = NsfFrozenMannInflow::new(snapshot(2.0), controls)
        .unwrap()
        .with_plane_wall(wall)
        .unwrap();
    let before = (grid.state_hash(), grid.force.clone(), grid.div.clone());
    let profile = grid.reservoir_aperture_tile_states(0).unwrap().collect::<Vec<_>>();
    let source = inlet.snapshot().checkpoint();
    assert!(matches!(inlet.apply(&mut grid, 0), Err(NsfFrozenMannError::Atmosphere(_))));
    assert_eq!((grid.state_hash(), grid.force.clone(), grid.div.clone()), before);
    assert_eq!(grid.reservoir_aperture_tile_states(0).unwrap().collect::<Vec<_>>(), profile);
    assert_eq!(inlet.snapshot().checkpoint(), source);
}
