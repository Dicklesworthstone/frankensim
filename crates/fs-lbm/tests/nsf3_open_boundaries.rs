//! Focused conservation/characteristic tests of the original NSF step.
//! Source packs below are authored algebra, not experimental gas validation.
#[path = "nsf3_thermochemical/source.rs"]
mod source;
#[path = "nsf3_open_boundaries/pressure.rs"]
mod pressure;
#[path = "nsf3_open_boundaries/patch.rs"]
mod patch;

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGrid3, NsfSpec, ThermochemicalNsfGas,
    WallThermal,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

fn spec<G>(gas: G, dims: [usize; 3], axis: usize, workers: usize) -> NsfSpec<G, String> {
    let mut periodic = [true; 3];
    periodic[axis] = false;
    NsfSpec {
        dims,
        periodic,
        gas,
        temperature_k: 300.0,
        pressure_pa: 100000.0,
        dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; synthetic-source cases are numerical fixtures".into(),
        maximum_cells: dims.iter().product(),
        wall: WallThermal::Adiabatic,
        workers,
        transport_points: 17,
    }
}

fn boundary(axis: usize, speed_m_s: f64) -> NsfCharacteristicAxis {
    let mut velocity = [0.0; 3];
    velocity[axis] = speed_m_s;
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: 100000.0,
        reference_temperature_k: 300.0,
        reference_velocity_m_s: velocity,
        incoming_characteristic_m_s: 0.0,
    };
    NsfCharacteristicAxis {
        axis,
        faces: [face; 2],
        maximum_relative_perturbation: 0.05,
    }
}

#[test]
fn uniform_throughflow_exchanges_actual_mass_momentum_enthalpy_and_work_on_every_axis() {
    let gas = GasSpec::dry_air_ussa1976();
    let state = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let speed = 5.0;
    let velocity_unit = (3.0 * state.specific_gas_constant * state.temperature).sqrt();
    for axis in 0..3 {
        let mut dims = [2; 3];
        dims[axis] = 8;
        let mut velocity = [0.0; 3];
        velocity[axis] = speed / velocity_unit;
        let mut grid = NsfGrid3::<D3q27>::new(
            spec(gas.clone(), dims, axis, 1),
            |_| false,
            |_| (1.0, velocity, 300.0),
        )
        .unwrap()
        .with_characteristic_axis(boundary(axis, speed))
        .unwrap();
        let energy = grid.total_energy_j();
        let mass = grid.total_mass_kg();
        let area = 4.0 * grid.law().dx_m.powi(2);
        let dt = grid.law().dt_s;
        for _ in 0..12 {
            let ledger = grid.step().unwrap();
            assert!(ledger.residual_rel < 2e-13);
            let exchange = grid.last_boundary_exchange().unwrap();
            for (side, e) in exchange.iter().enumerate() {
                let sign = if side == 0 { -1.0 } else { 1.0 };
                let mass_out = sign * state.density * speed * area * dt;
                let energy_out = mass_out * (state.specific_heat_cp * 300.0 + 0.5 * speed * speed);
                let momentum_out =
                    sign * (state.pressure + state.density * speed * speed) * area * dt;
                assert!((e.mass_out_kg - mass_out).abs() < 2e-11 * mass_out.abs());
                assert!((e.total_energy_out_j - energy_out).abs() < 2e-10 * energy_out.abs());
                assert!(
                    (e.advected_stagnation_enthalpy_out_j - energy_out).abs()
                        < 2e-10 * energy_out.abs()
                );
                assert!(e.nonadvective_energy_out_j.abs() < 2e-10 * energy_out.abs());
                assert_eq!(
                    e.total_energy_out_j,
                    e.advected_stagnation_enthalpy_out_j + e.nonadvective_energy_out_j
                );
                assert!(
                    (e.momentum_out_kg_m_s[axis] - momentum_out).abs() < 2e-11 * momentum_out.abs()
                );
                assert!(
                    (e.normal_volume_flow_m3_s - sign * speed * area).abs() < 2e-11 * speed * area
                );
                assert!(
                    (e.pressure_work_out_j - sign * state.pressure * speed * area * dt).abs()
                        < 2e-10 * state.pressure * speed * area * dt
                );
                assert!((e.mean_pressure_pa - state.pressure).abs() < 2e-10 * state.pressure);
            }
        }
        assert!((grid.total_energy_j() - energy).abs() < 2e-12 * energy);
        assert!((grid.total_mass_kg() - mass).abs() < 2e-12 * mass);
    }
}

#[test]
fn incoming_data_preserves_physical_state_then_drives_the_live_conservative_grid() {
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(GasSpec::dry_air_ussa1976(), [12, 2, 2], 0, 1),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap()
    .with_characteristic_axis(boundary(0, 0.0))
    .unwrap();
    let before = grid.cell([0, 0, 0]);
    let energy = grid.total_energy_j();
    // Low outward normal is -x. A negative incoming J- launches a +x
    // compression into the domain, without prescribing its outgoing wave.
    grid.set_incoming_characteristics([-0.2, 0.0]).unwrap();
    let after = grid.cell([0, 0, 0]);
    assert!((before.0 - after.0).abs() < 2e-15);
    assert!((before.2 - after.2).abs() < 2e-13);
    for a in 0..3 {
        assert!((before.1[a] - after.1[a]).abs() < 2e-15);
    }
    assert_eq!(grid.total_energy_j(), energy);
    let traces = grid.open_boundary_traces(false).unwrap();
    assert!(
        traces
            .iter()
            .all(|f| f.pressure_pa > 100000.0 && f.outward_normal_velocity_m_s < 0.0)
    );
    for _ in 0..8 {
        let m0 = grid.total_mass_kg();
        let p0 = grid.total_momentum_kg_m_s();
        let e0 = grid.total_energy_j();
        let ledger = grid.step().unwrap();
        let exchange = grid.last_boundary_exchange().unwrap();
        let dm = exchange.iter().map(|e| e.mass_out_kg).sum::<f64>();
        let de = exchange.iter().map(|e| e.total_energy_out_j).sum::<f64>();
        assert!((grid.total_mass_kg() - m0 + dm).abs() < 2e-12 * m0);
        assert!((grid.total_energy_j() - e0 + de).abs() < 2e-12 * e0);
        assert!((ledger.open_boundary_energy_j + de).abs() < 2e-12 * e0);
        let p1 = grid.total_momentum_kg_m_s();
        for a in 0..3 {
            let dp = exchange
                .iter()
                .map(|e| e.momentum_out_kg_m_s[a])
                .sum::<f64>();
            assert!((p1[a] - p0[a] + dp).abs() < 2e-12 * m0 * grid.law().dx_m / grid.law().dt_s);
        }
    }
    assert!(grid.cell([1, 0, 0]).1[0] > 0.0);
    assert!(grid.total_energy_j() > energy);
}

fn pulse_remainder(n: usize) -> f64 {
    let gas = GasSpec::dry_air_ussa1976();
    let g = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let c = (g.gamma / 3.0).sqrt();
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(gas, [n, 1, 1], 0, 1),
        |_| false,
        |x| {
            let z = (x[0] as f64 - 0.3 * n as f64) / (0.06 * n as f64);
            let amplitude = 1e-5 * fs_math::det::exp(-z * z);
            (
                1.0 + amplitude / g.gamma,
                [c * amplitude / g.gamma, 0.0, 0.0],
                300.0 * (1.0 + (g.gamma - 1.0) / g.gamma * amplitude),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(boundary(0, 0.0))
    .unwrap();
    let pressure_l2 = |grid: &NsfGrid3<D3q27>| -> f64 {
        (0..n)
            .map(|x| {
                let (rho, _, t) = grid.cell([x, 0, 0]);
                (rho * t / 300.0 - 1.0).powi(2)
            })
            .sum::<f64>()
    };
    let initial = pressure_l2(&grid);
    for _ in 0..(1.6 * n as f64 / c) as usize {
        grid.step().unwrap();
    }
    (pressure_l2(&grid) / initial).sqrt()
}

#[test]
fn outgoing_plane_pulse_reflection_decreases_with_boundary_resolution() {
    let coarse = pulse_remainder(64);
    let fine = pulse_remainder(128);
    // This tests the first-order boundary rung; it does not claim the plan's
    // -40 dB interface target. The initial/returned pulses are on the SAME
    // live source grid, with no external absorber or output filtering.
    assert!(
        coarse < 0.15,
        "coarse reflected/transient remainder {coarse}"
    );
    assert!(
        fine < 0.75 * coarse,
        "boundary refinement: {coarse} -> {fine}"
    );
}

#[test]
fn final_boundary_cell_domain_is_admitted_before_a_step_is_published() {
    let n = 48;
    let gas = GasSpec::dry_air_ussa1976();
    let state = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let c = (state.gamma / 3.0).sqrt();
    let mut open = boundary(0, 0.0);
    open.maximum_relative_perturbation = 1e-4;
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(gas, [n, 1, 1], 0, 1),
        |_| false,
        |x| {
            let z = (x[0] as f64 - 0.5 * n as f64) / (0.1 * n as f64);
            let pulse = 1e-3 * fs_math::det::exp(-z * z);
            (
                1.0 + pulse / state.gamma,
                [c * pulse / state.gamma, 0.0, 0.0],
                300.0 * (1.0 + (state.gamma - 1.0) / state.gamma * pulse),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(open)
    .unwrap();
    let mut refused = false;
    for _ in 0..100 {
        let before = (grid.state_hash(), grid.steps(), grid.cumulative_ledger());
        if grid.step().is_err() {
            assert_eq!(
                before,
                (grid.state_hash(), grid.steps(), grid.cumulative_ledger())
            );
            refused = true;
            break;
        }
        for x in [0, n - 1] {
            let (rho, u, t) = grid.cell([x, 0, 0]);
            let velocity = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt() / c;
            let change = (rho * t / 300.0 - 1.0)
                .abs()
                .max((t / 300.0 - 1.0).abs())
                .max(velocity);
            assert!(change <= open.maximum_relative_perturbation);
        }
    }
    assert!(
        refused,
        "the outgoing pulse must reach the declared boundary-domain limit"
    );
}

fn material_gas(offset: f64, pressure_max: f64) -> ThermochemicalNsfGas {
    ThermochemicalNsfGas::new(
        source::model(offset, pressure_max, 1.0),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap()
}

#[test]
fn source_caloric_chart_and_support_reach_open_faces_without_a_second_gas_model() {
    let build = |offset| {
        let gas = material_gas(offset, 250000.0);
        let reference = gas.model().evaluate(gas.context()).unwrap();
        let velocity_unit = (3.0 * reference.state().specific_gas_constant * 300.0).sqrt();
        NsfGrid3::<D3q27>::new(
            spec(gas, [8, 2, 2], 0, 1),
            |_| false,
            |x| {
                (
                    1.0,
                    [5.0 / velocity_unit, 0.0, 0.0],
                    300.0 + 0.01 * x[0] as f64,
                )
            },
        )
        .unwrap()
        .with_characteristic_axis(boundary(0, 5.0))
        .unwrap()
    };
    let mut a = build(-4000.0);
    let mut b = build(1_000_000.0);
    for _ in 0..5 {
        a.step().unwrap();
        b.step().unwrap();
    }
    for x in 0..8 {
        let (ra, ua, ta) = a.cell([x, 0, 0]);
        let (rb, ub, tb) = b.cell([x, 0, 0]);
        assert!((ra - rb).abs() < 1e-9);
        assert!((ta - tb).abs() < 2e-6);
        for k in 0..3 {
            assert!((ua[k] - ub[k]).abs() < 1e-9);
        }
    }
    for (ea, eb) in a
        .last_boundary_exchange()
        .unwrap()
        .iter()
        .zip(b.last_boundary_exchange().unwrap())
    {
        assert!(
            (ea.total_energy_out_j - eb.total_energy_out_j).abs()
                < 1e-8 * ea.total_energy_out_j.abs()
        );
    }
    assert!(a.thermochemical_state([0, 0, 0]).unwrap().is_some());
    let restricted = NsfGrid3::<D3q27>::new(
        spec(material_gas(-4000.0, 100500.0), [8, 2, 2], 0, 1),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap();
    let mut outside = boundary(0, 0.0);
    outside.faces[0].reference_pressure_pa = 101000.0;
    assert!(restricted.with_characteristic_axis(outside).is_err());
}

#[test]
fn boundary_data_refusals_are_atomic_and_workers_do_not_change_bits() {
    let build = |workers| {
        NsfGrid3::<D3q27>::new(
            spec(GasSpec::dry_air_ussa1976(), [8, 2, 2], 0, workers),
            |_| false,
            |_| (1.0, [0.0; 3], 300.0),
        )
        .unwrap()
        .with_characteristic_axis(boundary(0, 0.0))
        .unwrap()
    };
    let mut one = build(1);
    let snapshot = (
        one.state_hash(),
        one.total_energy_j(),
        *one.characteristic_axis().unwrap(),
        one.open_boundary_traces(false).unwrap(),
    );
    for bad in [f64::NAN, 100000.0] {
        assert!(one.set_incoming_characteristics([bad, 0.0]).is_err());
        assert_eq!(
            snapshot,
            (
                one.state_hash(),
                one.total_energy_j(),
                *one.characteristic_axis().unwrap(),
                one.open_boundary_traces(false).unwrap()
            )
        );
    }
    assert!(one.set_edge_velocity(0, false, [0.0, 0.01, 0.0]).is_err());
    let mut parallel = build(3);
    one.set_incoming_characteristics([-0.1, 0.0]).unwrap();
    parallel.set_incoming_characteristics([-0.1, 0.0]).unwrap();
    for _ in 0..4 {
        one.step().unwrap();
        parallel.step().unwrap();
    }
    assert_eq!(one.state_hash(), parallel.state_hash());
    assert_eq!(
        one.last_boundary_exchange(),
        parallel.last_boundary_exchange()
    );
    let blocked = NsfGrid3::<D3q27>::new(
        spec(GasSpec::dry_air_ussa1976(), [8, 2, 2], 0, 1),
        |x| x[0] == 1,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap();
    assert!(blocked.with_characteristic_axis(boundary(0, 0.0)).is_err());
}

// G0/G3: a finite rectangular passage uses the original wall and open owners.
// Cyclic coordinate permutations put the open axis before, between and after
// the wall axes. This exposes first-axis-wins corner streaming defects.
fn walled_passage(axis: usize, workers: usize, thin: bool) -> NsfGrid3<D3q27> {
    let canonical_dims = if thin { [8, 1, 1] } else { [10, 4, 3] };
    let dims = std::array::from_fn(|a| canonical_dims[(a + 3 - axis) % 3]);
    let gas = GasSpec::dry_air_ussa1976();
    let mut input = spec(gas, dims, axis, workers);
    input.periodic = [false; 3];
    NsfGrid3::<D3q27>::new(
        input,
        |_| false,
        |x| {
            let coordinate: [usize; 3] = std::array::from_fn(|a| x[(a + axis) % 3]);
            let u = if thin {
                [0.0; 3]
            } else {
                [
                    1e-5 * (1.0 + coordinate[1] as f64),
                    2e-6 * coordinate[2] as f64,
                    -1e-6 * coordinate[0] as f64,
                ]
            };
            let velocity = std::array::from_fn(|a| u[(a + 3 - axis) % 3]);
            let temperature = if thin {
                300.0
            } else {
                300.0 + 0.002 * (coordinate[0] + 2 * coordinate[1] + coordinate[2]) as f64
            };
            (1.0, velocity, temperature)
        },
    )
    .unwrap()
    .with_characteristic_axis(boundary(axis, 0.0))
    .unwrap()
}

#[test]
fn finite_transverse_walls_keep_actual_open_mass_and_total_energy_balanced() {
    for axis in 0..3 {
        let mut grid = walled_passage(axis, 1, false);
        grid.set_incoming_characteristics([-0.1, 0.02]).unwrap();
        for _ in 0..16 {
            let mass = grid.total_mass_kg();
            let energy = grid.total_energy_j();
            let ledger = grid.step().unwrap();
            let exchange = grid.last_boundary_exchange().unwrap();
            let mass_out = exchange.iter().map(|f| f.mass_out_kg).sum::<f64>();
            let energy_out = exchange.iter().map(|f| f.total_energy_out_j).sum::<f64>();
            assert!(
                (grid.total_mass_kg() - mass + mass_out).abs() < 3e-12 * mass,
                "actual corner mass on axis {axis}"
            );
            assert!(
                (grid.total_energy_j() - energy + energy_out - ledger.wall_heat_j).abs()
                    < 3e-12 * energy,
                "original finite-volume energy on axis {axis}"
            );
            assert!((ledger.open_boundary_energy_j + energy_out).abs() < 3e-12 * energy);
            assert_eq!(ledger.wall_heat_j, 0.0);
        }
    }
}

#[test]
fn wall_open_corner_ownership_is_independent_of_axis_order_and_worker_count() {
    let mut grids: Vec<_> = (0..3).map(|axis| walled_passage(axis, 1, false)).collect();
    let mut parallel = walled_passage(0, 3, false);
    for grid in grids.iter_mut().chain(std::iter::once(&mut parallel)) {
        grid.set_incoming_characteristics([-0.1, 0.02]).unwrap();
        for _ in 0..8 {
            grid.step().unwrap();
        }
    }
    assert_eq!(grids[0].state_hash(), parallel.state_hash());
    assert_eq!(
        grids[0].last_boundary_exchange(),
        parallel.last_boundary_exchange()
    );
    for axis in 1..3 {
        for z in 0..3 {
            for y in 0..4 {
                for x in 0..10 {
                    let coordinate = [x, y, z];
                    let rotated = std::array::from_fn(|a| coordinate[(a + 3 - axis) % 3]);
                    let (rho, u, t) = grids[0].cell(coordinate);
                    let (other_rho, other_u, other_t) = grids[axis].cell(rotated);
                    assert!((rho - other_rho).abs() < 2e-12);
                    assert!((t - other_t).abs() < 1e-9);
                    for a in 0..3 {
                        assert!((u[a] - other_u[(a + axis) % 3]).abs() < 2e-12);
                    }
                }
            }
        }
        for (base, other) in grids[0]
            .last_boundary_exchange()
            .unwrap()
            .iter()
            .zip(grids[axis].last_boundary_exchange().unwrap())
        {
            assert!(
                (base.normal_volume_flow_m3_s - other.normal_volume_flow_m3_s).abs()
                    < 1e-16 * grids[0].law().dx_m.powi(2) / grids[0].law().dt_s
            );
            assert!((base.mean_pressure_pa - other.mean_pressure_pa).abs() < 1e-6);
        }
    }
}

#[test]
fn opposing_wall_corners_do_not_open_a_one_cell_wide_passage_sideways() {
    for axis in 0..3 {
        let mut grid = walled_passage(axis, 1, true);
        for _ in 0..4 {
            grid.step().unwrap();
        }
        for coordinate in 0..8 {
            let mut x = [0; 3];
            x[axis] = coordinate;
            let (rho, u, t) = grid.cell(x);
            assert!((rho - 1.0).abs() < 2e-12);
            assert!(u.iter().all(|component| component.abs() < 2e-12));
            assert!((t - 300.0).abs() < 1e-9);
        }
        grid.set_incoming_characteristics([-0.1, 0.0]).unwrap();
        for _ in 0..6 {
            let mass = grid.total_mass_kg();
            grid.step().unwrap();
            let outgoing = grid
                .last_boundary_exchange()
                .unwrap()
                .iter()
                .map(|face| face.mass_out_kg)
                .sum::<f64>();
            assert!((grid.total_mass_kg() - mass + outgoing).abs() < 3e-12 * mass);
        }
        let mut x = [0; 3];
        x[axis] = 1;
        let (_, u, _) = grid.cell(x);
        assert!(u[axis] > 0.0);
        for a in 0..3 {
            if a != axis {
                assert!(u[a].abs() < 2e-12);
            }
        }
    }
}

#[test]
fn wall_heat_and_rejected_characteristics_remain_separate_and_atomic() {
    let mut input = spec(GasSpec::dry_air_ussa1976(), [8, 3, 3], 0, 1);
    input.periodic = [false; 3];
    input.wall = WallThermal::Isothermal(301.0);
    let mut grid = NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 300.0))
        .unwrap()
        .with_characteristic_axis(boundary(0, 0.0))
        .unwrap();
    let energy = grid.total_energy_j();
    let ledger = grid.step().unwrap();
    let outgoing = grid
        .last_boundary_exchange()
        .unwrap()
        .iter()
        .map(|face| face.total_energy_out_j)
        .sum::<f64>();
    assert!(ledger.wall_heat_j > 0.0);
    assert!(
        (grid.total_energy_j() - energy + outgoing - ledger.wall_heat_j).abs() < 3e-12 * energy
    );
    let before = (
        grid.state_hash(),
        grid.steps(),
        grid.cumulative_ledger(),
        *grid.characteristic_axis().unwrap(),
        grid.open_boundary_traces(false).unwrap(),
    );
    assert!(grid.set_incoming_characteristics([100000.0, 0.0]).is_err());
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.steps(),
            grid.cumulative_ledger(),
            *grid.characteristic_axis().unwrap(),
            grid.open_boundary_traces(false).unwrap()
        )
    );
}
