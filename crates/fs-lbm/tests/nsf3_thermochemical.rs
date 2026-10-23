//! Direct behavior regressions of the existing NSF solver with manufactured
//! NASA9/transport source data. These are not experimental gas validation.
#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{NsfError, NsfGasMedium, NsfGrid3, NsfSpec, ThermochemicalNsfGas, WallThermal};
use fs_material::gas::GasSpec;
use fs_material::thermochemical_gas::energy::GasEnergySolve;

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 600.0,
        energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 64,
    }
}
fn gas(offset: f64, pressure_max: f64, k_scale: f64) -> ThermochemicalNsfGas {
    ThermochemicalNsfGas::new(
        source::model(offset, pressure_max, k_scale),
        source::point(300.0, 100000.0),
        controls(),
    )
    .unwrap()
}
fn spec<G>(gas: G) -> NsfSpec<G> {
    NsfSpec {
        dims: [8, 2, 2],
        periodic: [true; 3],
        gas,
        temperature_k: 300.0,
        pressure_pa: 100000.0,
        dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored synthetic ratio; not empirical",
        maximum_cells: 32,
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 5,
    }
}
fn varying(x: [usize; 3]) -> (f64, [f64; 3], f64) {
    let phase = 2.0 * core::f64::consts::PI * x[0] as f64 / 8.0;
    (
        1.0 + 0.01 * phase.cos(),
        [0.02, 0.0, 0.0],
        300.0 + 0.2 * phase.sin(),
    )
}

#[test]
fn nasa9_energy_heating_updates_temperature_pressure_gamma_and_transport_together() {
    let gas = gas(-4000.0, 250000.0, 1.0);
    let reference = gas.model().evaluate(gas.context()).unwrap();
    assert!(reference.specific_internal_energy_j_kg() < 0.0);
    let r = reference.state().specific_gas_constant;
    let dx: f64 = 5e-5;
    let dt = dx / (3.0 * r * 300.0).sqrt();
    let specific_unit = (dx / dt).powi(2);
    let heat_per_mass = 10000.0;
    let source_e = heat_per_mass / specific_unit;
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(gas),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
    )
    .unwrap();
    let before = grid.total_energy_j();
    let mass = 32.0 * grid.law().rho_ref_kg_m3 * dx.powi(3);
    let steps = 4;
    for _ in 0..steps {
        grid.step().unwrap();
    }
    let (rho, _, temperature) = grid.cell([0, 0, 0]);
    // Independent integral of cv/R = 2 + .001 T, not another runtime EOS.
    let change =
        r * (2.0 * (temperature - 300.0) + 0.0005 * (temperature * temperature - 300.0 * 300.0));
    assert!((change - steps as f64 * heat_per_mass).abs() < 1e-4);
    assert!(
        (grid.total_energy_j() - before - mass * steps as f64 * heat_per_mass).abs()
            < 1e-10 * before.abs()
    );
    let resolved = grid.thermochemical_state([0, 0, 0]).unwrap().unwrap();
    let state = resolved.state();
    assert!((state.density - rho * grid.law().rho_ref_kg_m3).abs() < 1e-12);
    assert!((state.pressure - rho * grid.law().rho_ref_kg_m3 * r * temperature).abs() < 1e-8);
    assert!(state.gamma < reference.state().gamma);
    assert!(state.dynamic_viscosity > reference.state().dynamic_viscosity);
    assert!(state.thermal_conductivity > reference.state().thermal_conductivity);
    assert!(
        (temperature
            - (300.0 + steps as f64 * heat_per_mass / (reference.state().specific_heat_cp - r)))
            .abs()
            > 0.1
    );
}

#[test]
fn formation_energy_offset_cannot_drive_flow_or_heat_when_density_changes() {
    let mut a =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 1.0)), |_| false, varying).unwrap();
    let mut b =
        NsfGrid3::<D3q27>::new(spec(gas(1_000_000.0, 250000.0, 1.0)), |_| false, varying).unwrap();
    let rho_before = a.cell([0, 0, 0]).0;
    let energy_a = a.total_energy_j();
    let energy_b = b.total_energy_j();
    for _ in 0..5 {
        a.step().unwrap();
        b.step().unwrap();
    }
    assert!((a.cell([0, 0, 0]).0 - rho_before).abs() > 1e-6);
    for x in 0..8 {
        let (ra, ua, ta) = a.cell([x, 0, 0]);
        let (rb, ub, tb) = b.cell([x, 0, 0]);
        assert!((ra - rb).abs() < 1e-9);
        assert!((ta - tb).abs() < 2e-6, "{ta} vs {tb}");
        for axis in 0..3 {
            assert!((ua[axis] - ub[axis]).abs() < 1e-9);
        }
    }
    assert!((a.total_energy_j() - energy_a).abs() < 1e-12 * energy_a.abs());
    assert!((b.total_energy_j() - energy_b).abs() < 1e-12 * energy_b.abs());
    assert!((energy_a - energy_b).abs() < 1e-9 * energy_a.abs());
}

#[test]
fn conductivity_data_changes_native_energy_exchange_and_worker_count_does_not() {
    let initial = |x: [usize; 3]| {
        let t = 300.0 + 2.0 * (2.0 * core::f64::consts::PI * x[0] as f64 / 8.0).cos();
        (300.0 / t, [0.0; 3], t) // uniform initial physical pressure
    };
    let mut one =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 1.0)), |_| false, initial).unwrap();
    let mut parallel_spec = spec(gas(-4000.0, 250000.0, 1.0));
    parallel_spec.workers = 3;
    let mut parallel = NsfGrid3::<D3q27>::new(parallel_spec, |_| false, initial).unwrap();
    let mut conducting =
        NsfGrid3::<D3q27>::new(spec(gas(-4000.0, 250000.0, 4.0)), |_| false, initial).unwrap();
    for _ in 0..5 {
        one.step().unwrap();
        parallel.step().unwrap();
        conducting.step().unwrap();
    }
    assert_eq!(one.state_hash(), parallel.state_hash());
    assert!((one.cell([0, 0, 0]).2 - conducting.cell([0, 0, 0]).2).abs() > 1e-4);
    assert!(conducting.cell([0, 0, 0]).2 < one.cell([0, 0, 0]).2);
}

#[test]
fn actual_local_pressure_support_refuses_initial_and_future_states_atomically() {
    assert!(matches!(
        NsfGrid3::<D3q27>::new(
            spec(gas(-4000.0, 110000.0, 1.0)),
            |_| false,
            |_| (1.2, [0.0; 3], 300.0)
        ),
        Err(NsfError::Gas(_))
    ));
    let source = gas(-4000.0, 110000.0, 1.0);
    let reference = source.model().evaluate(source.context()).unwrap();
    let r = reference.state().specific_gas_constant;
    let source_e = 100000.0 / (3.0 * r * 300.0);
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(source),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| (0.0, [0.0; 3], source_e))),
    )
    .unwrap();
    let before = (
        grid.state_hash(),
        grid.total_energy_j(),
        grid.cell([0, 0, 0]),
        grid.cumulative_ledger(),
    );
    assert!(matches!(grid.step(), Err(NsfError::Gas(_))));
    assert_eq!(grid.steps(), 0);
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.total_energy_j(),
            grid.cell([0, 0, 0]),
            grid.cumulative_ledger()
        )
    );
    assert!(grid.step().is_err());
    assert_eq!(before.0, grid.state_hash());
}

#[test]
fn local_temperature_dependent_gamma_controls_mach_admission() {
    let source = gas(-4000.0, 250000.0, 1.0);
    let gamma_reference = source
        .model()
        .evaluate(source.context())
        .unwrap()
        .state()
        .gamma;
    let speed = 0.1999 * (gamma_reference * 550.0 / 900.0).sqrt();
    assert!(matches!(
        NsfGrid3::<D3q27>::new(spec(source), |_| false, |_| (1.0, [speed, 0.0, 0.0], 550.0)),
        Err(NsfError::Mach { .. })
    ));
}

#[test]
fn perfect_gas_ingress_retains_legacy_arithmetic_and_observations() {
    let gas = GasSpec::dry_air_ussa1976();
    let mut direct = NsfGrid3::<D3q27>::new(spec(gas.clone()), |_| false, varying).unwrap();
    let mut explicit = NsfGrid3::<D3q27>::new(
        spec(NsfGasMedium::CaloricallyPerfect(gas)),
        |_| false,
        varying,
    )
    .unwrap();
    assert!(direct.thermochemical_state([0, 0, 0]).unwrap().is_none());
    for _ in 0..5 {
        assert_eq!(direct.step().unwrap(), explicit.step().unwrap());
        assert_eq!(direct.state_hash(), explicit.state_hash());
    }
}
