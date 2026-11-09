use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{
    NsfDynamicThermalSgsConfig, NsfSpec, THETA0, WallThermal, collide, equilibrium, raw_moments,
};
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::reservoir::{
    ReservoirDynamicScalar, ReservoirDynamicSgs, ReservoirPrimitive,
};
use crate::sgs3::dynamic::scalar::DynamicScalarConfig;
use crate::sgs3::dynamic::{DynamicAveraging, DynamicSgsConfig};
use fs_material::gas::GasSpec;

const SHAPE: [usize; 3] = [5; 3];

fn settings(normal: usize) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: core::array::from_fn(|a| a != normal),
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored complete-reservoir profile fixture",
        maximum_cells: 125,
        workers: 1,
        transport_points: 9,
    }
}

fn uniform(grid: &NsfGrid3<D3q27>) -> NsfOpenState {
    NsfOpenState {
        density_kg_m3: grid.law.rho_ref_kg_m3,
        velocity_m_s: [0.0; 3],
        temperature_k: 293.15,
    }
}

fn grid(normal: usize) -> NsfGrid3<D3q27> {
    let mut grid = NsfGrid3::<D3q27>::new(
        settings(normal),
        |_| false,
        |p| {
            let i = p[0] + 5 * (p[1] + 5 * p[2]);
            (
                1.0 + ((i * 7 + 3) % 11) as f64 / 2048.0,
                core::array::from_fn(|a| {
                    (((i * (17 + 8 * a) + i * i * (a + 1) + 13 * a) % 127) as f64 - 63.0) / 32768.0
                }),
                293.15 + ((i * i * 3 + 7 * i) % 41) as f64 / 128.0,
            )
        },
    )
    .unwrap();
    let state = uniform(&grid);
    grid.set_open_face(normal, false, Some(state)).unwrap();
    grid.set_open_face(normal, true, Some(state)).unwrap();
    grid
}

fn controls() -> NsfDynamicThermalSgsConfig {
    NsfDynamicThermalSgsConfig {
        momentum: DynamicSgsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn profile(grid: &NsfGrid3<D3q27>, normal: usize, high: bool) -> Vec<NsfReservoirTileState> {
    grid.open_face_tile_states(normal, high)
        .unwrap()
        .enumerate()
        .map(|(j, mut tile)| {
            tile.state.density_kg_m3 = grid.law.rho_ref_kg_m3 * (1.0 + (j % 7) as f64 / 2048.0);
            tile.state.velocity_m_s = core::array::from_fn(|a| {
                (((j * (3 + 5 * a) + 7 * a) % 29) as f64 - 14.0) / 8192.0 * grid.law.dx_m
                    / grid.law.dt_s
            });
            tile.state.temperature_k = 293.0 + (j % 11) as f64 / 32.0;
            tile
        })
        .collect()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-12 * b.abs().max(1e-5),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g3_uniform_face_profiles_recover_original_controls_and_canonical_state_order() {
    for normal in 0..3 {
        for high in [false, true] {
            let mut scalar = grid(normal);
            let mut tiled = grid(normal);
            let mut state = uniform(&scalar);
            state.temperature_k = 293.5;
            state.velocity_m_s[normal] = 0.001 * scalar.law.dx_m / scalar.law.dt_s;
            let mut tiles: Vec<_> = tiled
                .open_face_tile_states(normal, high)
                .unwrap()
                .map(|mut tile| {
                    tile.state = state;
                    tile
                })
                .collect();
            tiles.reverse();
            scalar.set_open_face(normal, high, Some(state)).unwrap();
            tiled
                .set_open_face_profile(normal, high, &tiles, 25, || false)
                .unwrap();
            assert_eq!(tiled.state_hash(), scalar.state_hash());
            assert_eq!(tiled.f, scalar.f);
            assert_eq!(tiled.force, scalar.force);
            assert_eq!(tiled.div, scalar.div);
            assert!(tiled.reservoir_face_profiles[normal][usize::from(high)].is_none());
            for (i, tile) in tiled
                .open_face_tile_states(normal, high)
                .unwrap()
                .enumerate()
            {
                assert_eq!(tile.cell, face_cell(SHAPE, normal, high, i));
            }
        }
    }
}

fn raw_source(q: usize, u: [f64; 3], force: [f64; 3]) -> f64 {
    let c = D3q27::E[q].map(f64::from);
    let cu = (0..3).map(|a| c[a] * u[a]).sum::<f64>();
    let cf = (0..3).map(|a| c[a] * force[a]).sum::<f64>();
    let uf = (0..3).map(|a| u[a] * force[a]).sum::<f64>();
    D3q27::W[q] * ((cf - uf) / THETA0 + cu * cf / (THETA0 * THETA0))
}

#[test]
fn g1_diagonal_incoming_populations_use_actual_profile_donors_and_native_mass_energy() {
    let mut distinguished = false;
    for normal in 0..3 {
        for high in [false, true] {
            let mut grid = grid(normal);
            let tiles = profile(&grid, normal, high);
            grid.set_open_face_profile(normal, high, &tiles, 25, || false)
                .unwrap();
            let post: Vec<_> = (0..grid.n)
                .map(|i| {
                    let (rho, _) = raw_moments::<D3q27>(&grid.f[i]);
                    let (ws, wb) = grid.rates(
                        rho,
                        grid.theta[i],
                        grid.material.get(i).and_then(Option::as_ref),
                    );
                    collide::<D3q27>(&grid.f[i], grid.force[i], ws, wb)
                })
                .collect();
            let streamed = grid.stream(&post, None);
            let [first, second] = axes(normal);
            let mut entering = 0.0;
            for i in 0..grid.n {
                let p = grid.at(i);
                if p[normal] != if high { 4 } else { 0 } {
                    continue;
                }
                let (rho, momentum) = raw_moments::<D3q27>(&post[i]);
                let u = core::array::from_fn(|a| (momentum[a] - 0.5 * grid.force[i][a]) / rho);
                let local = equilibrium::<D3q27>(rho, u, THETA0);
                for q in 0..27 {
                    let c = D3q27::E[q];
                    if c[normal] != if high { -1 } else { 1 } {
                        continue;
                    }
                    let a = (p[first] as i32 - c[first]).rem_euclid(5) as usize;
                    let b = (p[second] as i32 - c[second]).rem_euclid(5) as usize;
                    let state = tiles[a + 5 * b].state;
                    let ghost_u = state
                        .velocity_m_s
                        .map(|v| v * grid.law.dt_s / grid.law.dx_m);
                    let ghost = equilibrium::<D3q27>(
                        state.density_kg_m3 / grid.law.rho_ref_kg_m3,
                        ghost_u,
                        THETA0,
                    );
                    let expected =
                        ghost[q] + 0.5 * raw_source(q, ghost_u, grid.force[i]) + post[i][q]
                            - local[q]
                            - 0.5 * raw_source(q, u, grid.force[i]);
                    close(streamed[i][q], expected);
                    let wrong_state = tiles[p[first] + 5 * p[second]].state;
                    let wrong_u = wrong_state
                        .velocity_m_s
                        .map(|v| v * grid.law.dt_s / grid.law.dx_m);
                    let wrong = equilibrium::<D3q27>(
                        wrong_state.density_kg_m3 / grid.law.rho_ref_kg_m3,
                        wrong_u,
                        THETA0,
                    );
                    distinguished |= (wrong[q] - ghost[q]).abs() > 1e-7;
                }
            }
            for i in 0..grid.n {
                let p = grid.at(i);
                for q in 0..27 {
                    let c = D3q27::E[q][normal];
                    if (p[normal] == 0 && c == 1) || (p[normal] == 4 && c == -1) {
                        entering += streamed[i][q] - post[i][D3q27::OPP[q]];
                    }
                }
            }
            let before_mass = grid.f.iter().flatten().sum::<f64>();
            let before_energy = grid.total_energy_j();
            let ledger = grid.step().unwrap();
            assert!((grid.f.iter().flatten().sum::<f64>() - before_mass - entering).abs() < 4e-12);
            assert!(
                (grid.total_energy_j() - before_energy - ledger.open_energy_j).abs()
                    < 2e-12 * before_energy.abs()
            );
            assert!(ledger.residual_rel < 2e-12);
        }
    }
    assert!(distinguished);
}

#[test]
fn g1_dynamic_profiles_refresh_physical_ghosts_and_survive_reversible_native_steps() {
    let mut grid = grid(0)
        .with_reservoir_dynamic_thermal_sgs(controls(), 25_000)
        .unwrap();
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
    let energy = grid.e.clone();
    let tiles = profile(&grid, 0, false);
    grid.set_open_face_profile(0, false, &tiles, 25, || false)
        .unwrap();
    let (r, m) = grid.macro_fields(&grid.f);
    let u = grid.velocity_field(&r, &m, &grid.force);
    for i in 0..125 {
        close(r[i], rho[i]);
        for a in 0..3 {
            assert!((u[i][a] - velocity[i][a]).abs() < 3e-15);
        }
    }
    assert_eq!(grid.e, energy);
    let mut model = ReservoirDynamicSgs::new(SHAPE, [1.0; 3], 0, usize::MAX).unwrap();
    let mut heat = ReservoirDynamicScalar::new(SHAPE, 0, usize::MAX).unwrap();
    model
        .prepare_with_profile(
            controls().momentum,
            controls().heat,
            &r,
            &u,
            &grid.theta,
            |side, i| {
                let g = grid.reservoir_ghost(i, 0, side).unwrap();
                ReservoirPrimitive {
                    density: g.rho,
                    velocity: g.u,
                    scalar: g.theta,
                }
            },
            &mut heat,
            || false,
        )
        .unwrap();
    for i in 0..125 {
        close(
            grid.sgs_transport[i].eddy_viscosity_lu,
            model.viscosities()[i],
        );
        let cp = grid
            .material
            .get(i)
            .and_then(Option::as_ref)
            .map_or(grid.law.gamma / (grid.law.gamma - 1.0), |m| {
                m.state.specific_heat_cp / grid.law.gas_constant
            });
        close(
            grid.sgs_transport[i].turbulent_conductivity_lu,
            r[i] * cp * heat.diffusivities()[i],
        );
    }
    let saved = grid.state_hash();
    let accepted = grid.reservoir_face_profiles[0][0].clone().unwrap();
    let trial = grid
        .preview_thermal_wall_step(Vec::new(), 125, || false)
        .unwrap();
    assert!(Arc::ptr_eq(
        &accepted,
        trial.grid().reservoir_face_profiles[0][0].as_ref().unwrap()
    ));
    assert_eq!(
        trial
            .grid()
            .open_face_tile_states(0, false)
            .unwrap()
            .collect::<Vec<_>>(),
        tiles
    );
    let expected = trial.grid().state_hash();
    drop(trial);
    assert_eq!(grid.state_hash(), saved);
    grid.preview_thermal_wall_step(Vec::new(), 125, || false)
        .unwrap()
        .accept();
    assert_eq!(grid.state_hash(), expected);
}

#[test]
fn g0_g4_profile_staging_caps_errors_and_every_observed_cancellation_roll_back() {
    let mut probe = grid(2);
    let tiles = profile(&probe, 2, true);
    let mut polls = 0;
    probe
        .set_open_face_profile(2, true, &tiles, 25, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = probe.state_hash();
    for stop in 1..=polls {
        let mut tested = grid(2);
        let before = tested.state_hash();
        let force = tested.force.clone();
        let div = tested.div.clone();
        let mut seen = 0;
        assert!(
            tested
                .set_open_face_profile(2, true, &tiles, 25, || {
                    seen += 1;
                    seen == stop
                })
                .is_err()
        );
        assert_eq!(tested.state_hash(), before);
        assert_eq!(tested.force, force);
        assert_eq!(tested.div, div);
        assert!(tested.reservoir_face_profiles[2][1].is_none());
        tested
            .set_open_face_profile(2, true, &tiles, 25, || false)
            .unwrap();
        assert_eq!(tested.state_hash(), expected);
    }
    for kind in 0..5 {
        let before = probe.state_hash();
        let mut bad = tiles.clone();
        match kind {
            0 => bad[1] = bad[0],
            1 => bad[0].cell[2] = 0,
            2 => bad[0].state.temperature_k = 400.0,
            3 => {
                bad.pop();
            }
            _ => {}
        }
        assert!(
            probe
                .set_open_face_profile(2, true, &bad, if kind == 4 { 24 } else { 25 }, || false)
                .is_err()
        );
        assert_eq!(probe.state_hash(), before);
        assert_eq!(
            probe
                .open_face_tile_states(2, true)
                .unwrap()
                .collect::<Vec<_>>(),
            tiles
        );
    }
}
