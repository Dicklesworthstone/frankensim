//! G0/G1 physical gravity, work and original compressible-source coupling.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfError, NsfGrid3, NsfSpec,
    THETA0, WallThermal,
};
use fs_material::gas::{GasSpec, GasState};
use std::cell::Cell;

fn spec() -> NsfSpec {
    NsfSpec {
        dims: [4, 3, 2], periodic: [true; 3],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0, pressure_pa: 100_000.0, dx_m: 1.0,
        temperature_window_k: [250.0, 400.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 0.3, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; gravity regression fixture",
        maximum_cells: 48, wall: WallThermal::Adiabatic,
        workers: 1, transport_points: 17,
    }
}

fn grid(rho: f64, velocity: [f64; 3]) -> NsfGrid3<D3q27> {
    NsfGrid3::new(spec(), |_| false, |_| (rho, velocity, 300.0)).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}");
}

#[test]
fn uniform_gas_accelerates_independently_of_density_without_heating() {
    let g = [2.0, -9.80665, 1.5];
    let initial = [0.003, -0.004, 0.002];
    let mut work_per_density = None;
    for density in [0.6, 1.0, 1.8] {
        let mut flow = grid(density, initial).with_gravity_si(g).unwrap();
        let unit = flow.law().dx_m / flow.law().dt_s;
        let dt = flow.law().dt_s;
        let mass = flow.total_mass_kg();
        let energy = flow.total_energy_j();
        let internal = energy - flow.kinetic_energy_j();
        let mut work = 0.0;
        for step in 1..=24 {
            let ledger = flow.step().unwrap();
            work += ledger.gravity_work_j;
            assert!(ledger.residual_rel < 2e-13);
            assert_eq!(ledger.source_energy_j, 0.0);
            let (rho, velocity, temperature) = flow.cell([1, 1, 1]);
            close(rho, density, 5e-13 * density);
            close(temperature, 300.0, 3e-9);
            for a in 0..3 {
                close(velocity[a], initial[a] + g[a] * (step as f64 * dt) / unit, 3e-13);
            }
        }
        let expected_velocity: [f64; 3] = core::array::from_fn(|a| {
            initial[a] * unit + g[a] * (24.0 * dt)
        });
        let expected_work = 0.5 * mass * (0..3).map(|a| {
            expected_velocity[a].powi(2) - (initial[a] * unit).powi(2)
        }).sum::<f64>();
        close(work, expected_work, 2e-8 * expected_work.abs().max(1.0));
        close(flow.total_mass_kg(), mass, 5e-13 * mass);
        close(flow.total_energy_j() - energy, work, 3e-12 * energy);
        close(flow.total_energy_j() - flow.kinetic_energy_j(), internal, 3e-12 * energy);
        close(flow.cumulative_ledger().viscous_heating_j, 0.0, 3e-12 * energy);
        assert_eq!(flow.cumulative_ledger().gravity_work_j, work);
        if let Some(reference) = work_per_density {
            close(work / density, reference, 2e-8 * (work / density).abs().max(1.0));
        } else {
            work_per_density = Some(work / density);
        }
    }
}

#[test]
fn gravity_uses_evolving_density_instead_of_a_frozen_force() {
    let density = 0.8;
    let mass_rate = 0.002;
    let initial = [0.002, -0.001, 0.003];
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let internal_source = mass_rate * THETA0 / (gas.gamma - 1.0);
    let mut flow = NsfGrid3::<D3q27>::new_with_sources(
        spec(), |_| false, |_| (density, initial, 300.0),
        Some(Box::new(move |_, _| (mass_rate, [0.0; 3], internal_source))),
    ).unwrap().with_gravity_si([1.0, -9.80665, 2.0]).unwrap();
    let acceleration = flow.gravity_m_s2().map(|g| {
        g * flow.law().dt_s * (flow.law().dt_s / flow.law().dx_m)
    });
    for step in 1..=12 {
        flow.step().unwrap();
        let t = step as f64;
        let (rho, velocity, _) = flow.cell([1, 1, 1]);
        close(rho, density + mass_rate * t, 5e-13);
        for a in 0..3 {
            // Exact integral of rho(t)*g, with incoming mass at rest.
            let expected = density * initial[a]
                + acceleration[a] * (density * t + 0.5 * mass_rate * t * t);
            close(rho * velocity[a], expected, 5e-13);
        }
    }
}

#[test]
fn uniform_heat_and_gravitational_work_remain_separate() {
    let heat_rate = 2e-5;
    let make = || NsfGrid3::<D3q27>::new_with_sources(
        spec(), |_| false, |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], heat_rate))),
    ).unwrap();
    let mut unforced = make();
    let mut falling = make().with_gravity_si([0.0, -9.80665, 0.0]).unwrap();
    let energy = falling.total_energy_j();
    for _ in 0..12 {
        let reference = unforced.step().unwrap();
        let actual = falling.step().unwrap();
        close(actual.source_energy_j, reference.source_energy_j, 1e-12);
        assert!(actual.gravity_work_j > 0.0);
        close(falling.cell([1, 1, 1]).2, unforced.cell([1, 1, 1]).2, 3e-9);
        close(actual.viscous_heating_j, 0.0, 3e-12 * energy);
        assert!(actual.residual_rel < 2e-13);
    }
}

#[test]
fn changing_gravity_preserves_physical_state_and_zero_is_an_exact_no_op() {
    let make = || NsfGrid3::<D3q27>::new(spec(), |_| false, |x| {
        (1.0 + 0.001 * (x[0] as f64),
            [0.005, 0.002 * (x[1] as f64), -0.003], 300.0)
    }).unwrap();
    let mut plain = make();
    let mut zero = make().with_gravity_si([0.0, -0.0, 0.0]).unwrap();
    assert_eq!(plain.state_hash(), zero.state_hash());
    for _ in 0..3 {
        assert_eq!(plain.step().unwrap(), zero.step().unwrap());
        assert_eq!(plain.state_hash(), zero.state_hash());
    }
    let snapshot = (plain.steps(), plain.cumulative_ledger(), plain.total_energy_j());
    let mut before = Vec::new();
    for z in 0..2 { for y in 0..3 { for x in 0..4 {
        before.push(([x, y, z], plain.cell([x, y, z])));
    }}}
    for gravity in [[2.0, -9.80665, 1.0], [-4.0, 3.0, -2.0], [0.0; 3]] {
        plain.set_gravity_si(gravity).unwrap();
        assert_eq!(plain.gravity_m_s2(), gravity);
        assert_eq!((plain.steps(), plain.cumulative_ledger(), plain.total_energy_j()), snapshot);
        for &(cell, (rho, velocity, temperature)) in &before {
            let actual = plain.cell(cell);
            close(actual.0, rho, 5e-14);
            for a in 0..3 { close(actual.1[a], velocity[a], 5e-14); }
            assert_eq!(actual.2, temperature);
        }
    }
}

#[test]
fn invalid_gravity_and_a_refused_acceleration_step_are_transactional() {
    let mut flow = grid(1.0, [0.0; 3]);
    flow.set_gravity_si([0.0, -9.80665, 0.0]).unwrap();
    let snapshot = (flow.state_hash(), flow.steps(), flow.gravity_m_s2(), flow.cumulative_ledger());
    for gravity in [[f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0], [f64::from_bits(1), 0.0, 0.0]] {
        assert!(flow.set_gravity_si(gravity).is_err());
        assert_eq!((flow.state_hash(), flow.steps(), flow.gravity_m_s2(), flow.cumulative_ledger()), snapshot);
    }
    let excessive = 0.5 * flow.law().dx_m / flow.law().dt_s.powi(2);
    flow.set_gravity_si([excessive, 0.0, 0.0]).unwrap();
    let snapshot = (flow.state_hash(), flow.steps(), flow.gravity_m_s2(), flow.cumulative_ledger());
    assert!(matches!(flow.step(), Err(NsfError::Mach { .. })));
    assert_eq!((flow.state_hash(), flow.steps(), flow.gravity_m_s2(), flow.cumulative_ledger()), snapshot);
    flow.set_gravity_si([0.0; 3]).unwrap();
    flow.step().unwrap();
}

#[test]
fn characteristic_previews_retain_gravity_work_and_cancel_without_publication() {
    let make = || {
        let mut input = spec();
        input.dims = [12, 2, 2];
        input.periodic = [false, true, true];
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: 100_000.0, reference_temperature_k: 300.0,
            reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
        };
        NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 300.0))
            .unwrap().with_characteristic_axis(NsfCharacteristicAxis {
                axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
            }).unwrap().with_gravity_si([0.0, -9.80665, 0.0]).unwrap()
    };
    let mut direct = make();
    let mut previewed = make();
    direct.set_incoming_characteristics([-0.1, 0.03]).unwrap();
    let mut work = 0.0;
    for _ in 0..3 { work += direct.step().unwrap().gravity_work_j; }
    let trial = previewed.preview_characteristic_interval([-0.1, 0.03], 3, 144, || false).unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.energy_ledger().gravity_work_j, work);
    assert!(work > 0.0);
    assert!(trial.energy_ledger().residual_rel < 2e-13);
    trial.accept();
    assert_eq!(previewed.cumulative_ledger(), direct.cumulative_ledger());
    assert_eq!(previewed.gravity_m_s2(), [0.0, -9.80665, 0.0]);
    let before = (previewed.state_hash(), previewed.steps(), previewed.cumulative_ledger());
    let calls = Cell::new(0);
    assert!(previewed.preview_characteristic_interval([-0.02, 0.01], 3, 144, || {
        let n = calls.get();
        calls.set(n + 1);
        n >= 2
    }).is_err());
    assert!(calls.get() >= 3);
    assert_eq!((previewed.state_hash(), previewed.steps(), previewed.cumulative_ledger()), before);
}
