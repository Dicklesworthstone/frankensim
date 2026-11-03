use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::DynamicAveraging;
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[path = "../../../../tests/nsf3_thermochemical/source.rs"]
mod source;

const SHAPE: [usize; 3] = [5, 6, 7];

fn spec(workers: usize) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: [true; 3],
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored dynamic-SGS coupling fixture",
        maximum_cells: 4096,
        workers,
        transport_points: 9,
    }
}

fn config() -> NsfDynamicThermalSgsConfig {
    NsfDynamicThermalSgsConfig {
        momentum: DynamicSgsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient_squared: 16.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 16.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn mixed(p: [usize; 3]) -> (f64, [f64; 3], f64) {
    let [x, y, z] =
        core::array::from_fn(|a| core::f64::consts::TAU * (p[a] as f64 + 0.5) / SHAPE[a] as f64);
    (
        1.0 + 0.02 * x.sin() + 0.01 * y.cos(),
        [
            0.012 * x.sin() + 0.008 * (y + z).sin(),
            0.010 * (x - z).cos() + 0.006 * y.sin(),
            0.008 * (x + y).sin() + 0.004 * z.cos(),
        ],
        293.15 + 0.75 * (x + y).sin() + 0.25 * z.cos(),
    )
}

fn rest(_: [usize; 3]) -> (f64, [f64; 3], f64) {
    (1.0, [0.0; 3], 293.15)
}

fn periodic_workspace(grid: &NsfGrid3<D3q27>) -> std::sync::MutexGuard<'_, Workspace> {
    let DynamicState::Periodic { workspace, .. } = &grid.dynamic_sgs.as_ref().unwrap().state else {
        panic!("expected the periodic workspace");
    };
    workspace.lock().unwrap()
}

#[test]
fn g1_dynamic_nsf_uses_force_corrected_velocity_and_current_nasa_heat_capacity() {
    let gas = ThermochemicalNsfGas::new(
        source::model(-4000.0, 250000.0, 1.0),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    let mut settings = spec(1).map_gas(|_| gas);
    settings.temperature_k = 300.0;
    settings.pressure_pa = 100000.0;
    settings.temperature_window_k = [250.0, 600.0];
    let grid = NsfGrid3::<D3q27>::new_with_sources(
        settings,
        |_| false,
        mixed,
        Some(Box::new(|p, _| {
            let phase = core::f64::consts::TAU * (p[2] as f64 + 0.5) / SHAPE[2] as f64;
            (0.0, [0.001 * phase.sin(), 0.0, 0.0], 0.0)
        })),
    )
    .unwrap()
    .with_dynamic_thermal_sgs(config(), 200 * SHAPE.iter().product::<usize>())
    .unwrap();
    let (rho, mom) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &mom, &grid.force);
    let mut core = PeriodicDynamicSgs::new(SHAPE, 1.0, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new(SHAPE, usize::MAX).unwrap();
    let reports = core
        .prepare_with_scalar(
            config().momentum,
            config().heat,
            &rho,
            &velocity,
            &grid.theta,
            &mut scalar,
            || false,
        )
        .unwrap();
    assert!(reports.0.maximum_eddy_viscosity > 0.0);
    assert!(reports.1.maximum_diffusivity > 0.0);
    for i in 0..grid.n {
        let local = grid.material[i].as_ref().unwrap();
        let cp_over_r = local.state.specific_heat_cp / grid.law.gas_constant;
        let temperature = grid.law.kelvin_of_theta(grid.theta[i]);
        assert!((cp_over_r - (3.0 + 0.001 * temperature)).abs() < 1e-12);
        assert_eq!(
            grid.sgs_transport[i].eddy_viscosity_lu,
            core.viscosities()[i]
        );
        let expected = rho[i] * cp_over_r * scalar.diffusivities()[i];
        assert_eq!(grid.sgs_transport[i].turbulent_conductivity_lu, expected);
        let molecular = grid.rates(rho[i], grid.theta[i], Some(local));
        let fitted = grid.sgs_rates(
            rho[i],
            grid.theta[i],
            Some(local),
            Some(&grid.sgs_transport[i]),
        );
        assert_eq!(fitted.1.to_bits(), molecular.1.to_bits());
    }
    let physical = core.viscosities().to_vec();
    let raw: Vec<_> = mom
        .iter()
        .zip(&rho)
        .map(|(m, r)| m.map(|v| v / r))
        .collect();
    core.prepare(config().momentum, &rho, &raw, || false)
        .unwrap();
    assert!(
        core.viscosities()
            .iter()
            .zip(physical)
            .any(|(a, b)| (a - b).abs() > 1e-10)
    );
}

#[test]
fn g1_dynamic_nsf_heat_only_flux_is_conservative_and_uses_actual_face_limit() {
    let mut settings = spec(1);
    settings.dims = [9; 3];
    let grid = NsfGrid3::<D3q27>::new(settings, |_| false, rest)
        .unwrap()
        .with_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap();
    let rho = vec![1.0; grid.n];
    let mut velocity = Vec::new();
    let mut theta = Vec::new();
    for i in 0..grid.n {
        let p = grid
            .at(i)
            .map(|a| core::f64::consts::TAU * (a as f64 + 0.5) / 9.0);
        velocity.push([0.02 * p[1].sin(), 0.0, 0.0]);
        theta.push(THETA0 + 0.001 * (p[0] + p[1]).sin());
    }
    let transport = grid.sgs_field(&rho, &velocity, &theta, &[], None).unwrap();
    assert!(transport.iter().all(|t| t.eddy_viscosity_lu == 0.0));
    assert!(transport.iter().any(|t| t.turbulent_conductivity_lu > 0.0));
    let stress = vec![[0.0; 6]; grid.n];
    let molecular = grid
        .flux_divergence(&grid.e, &rho, &velocity, &theta, &stress, &[], None)
        .unwrap()
        .0;
    let modeled = grid
        .flux_divergence(&grid.e, &rho, &velocity, &theta, &stress, &transport, None)
        .unwrap()
        .0;
    assert!(modeled.iter().sum::<f64>().abs() < 1e-14);
    assert!(
        modeled
            .iter()
            .zip(&molecular)
            .any(|(a, b)| (a - b).abs() > 1e-10)
    );
    let constant = vec![THETA0; grid.n];
    let k = grid.law.kappa_lattice(THETA0);
    let excessive = vec![
        NsfSgsTransport {
            eddy_viscosity_lu: 0.0,
            turbulent_conductivity_lu: 1001.0 * k,
        };
        grid.n
    ];
    assert!(matches!(
        grid.flux_divergence(
            &grid.e, &rho, &velocity, &constant, &stress, &excessive, None
        ),
        Err(NsfError::SgsHeatRatio { .. })
    ));
}

#[test]
fn g0_g4_dynamic_nsf_zero_controls_budgets_and_topology_are_explicit() {
    let make = || NsfGrid3::<D3q27>::new(spec(1), |_| false, rest).unwrap();
    let bytes = make().required_dynamic_sgs_workspace_bytes().unwrap();
    assert_eq!(bytes, 200 * SHAPE.iter().product::<usize>());
    assert!(matches!(
        make().with_dynamic_thermal_sgs(config(), bytes - 1),
        Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes { .. }))
    ));
    let mut plain = make();
    let mut fitted = make().with_dynamic_thermal_sgs(config(), bytes).unwrap();
    assert_eq!(fitted.dynamic_sgs_workspace_bytes(), Some(bytes));
    assert_eq!(fitted.dynamic_thermal_sgs_config(), Some(config()));
    assert_eq!(fitted.sgs_config(), None);
    for _ in 0..4 {
        assert_eq!(plain.step().unwrap(), fitted.step().unwrap());
        assert_eq!(plain.state_hash(), fitted.state_hash());
    }
    for dims in [[2, 6, 7], SHAPE] {
        let mut settings = spec(1);
        settings.dims = dims;
        if dims == SHAPE {
            settings.periodic[0] = false;
        }
        let grid = NsfGrid3::<D3q27>::new(settings, |_| false, rest).unwrap();
        assert!(grid.with_dynamic_thermal_sgs(config(), usize::MAX).is_err());
    }
    let solid = NsfGrid3::<D3q27>::new(spec(1), |p| p == [1, 1, 1], rest).unwrap();
    assert!(
        solid
            .with_dynamic_thermal_sgs(config(), usize::MAX)
            .is_err()
    );
    let hash = fitted.state_hash();
    assert!(
        fitted
            .set_edge_velocity(0, false, [0.0, 0.01, 0.0])
            .is_err()
    );
    assert!(fitted.set_open_face(0, false, None).is_err());
    assert_eq!(fitted.state_hash(), hash);
    let local = NsfSgsConfig {
        closure: SgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.17 },
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.9,
    };
    assert!(fitted.with_sgs(local).is_err());
    assert!(
        make()
            .with_sgs(local)
            .unwrap()
            .with_dynamic_thermal_sgs(config(), bytes)
            .is_err()
    );
    let mut acceptance = config();
    acceptance.heat.usage = SgsUse::AcceptanceRender;
    assert!(make().with_dynamic_thermal_sgs(acceptance, bytes).is_err());
}

#[test]
fn g1_g5_dynamic_nsf_reuses_buffers_and_keeps_original_energy_and_worker_bits() {
    let make = |workers| {
        NsfGrid3::<D3q27>::new(spec(workers), |_| false, mixed)
            .unwrap()
            .with_dynamic_thermal_sgs(config(), usize::MAX)
            .unwrap()
    };
    let mut serial = make(1);
    let mut parallel = make(3);
    let energy = serial.total_energy_j();
    let allocation = {
        let work = periodic_workspace(&serial);
        (
            work.momentum.viscosities().as_ptr(),
            work.heat.diffusivities().as_ptr(),
        )
    };
    for _ in 0..6 {
        let step = serial.step().unwrap();
        assert_eq!(step, parallel.step().unwrap());
        assert_eq!(serial.state_hash(), parallel.state_hash());
        assert_eq!(serial.sgs_transport, parallel.sgs_transport);
        assert!(step.residual_rel < 1e-12);
        assert!((serial.total_energy_j() / energy - 1.0).abs() < 2e-12);
        assert_eq!(step.wall_heat_j, 0.0);
        let work = periodic_workspace(&serial);
        assert_eq!(work.momentum.viscosities().as_ptr(), allocation.0);
        assert_eq!(work.heat.diffusivities().as_ptr(), allocation.1);
    }
    // A newly refused fitted state must not publish scratch as accepted
    // transport or advance the original flow/energy transaction.
    let hash = serial.state_hash();
    let transport = serial.sgs_transport.clone();
    let ledger = serial.cumulative_ledger();
    let steps = serial.steps();
    let DynamicState::Periodic { config, .. } = &mut serial.dynamic_sgs.as_mut().unwrap().state
    else {
        panic!("expected periodic controls");
    };
    config.heat.maximum_diffusivity_to_molecular = 0.0;
    assert!(matches!(serial.step(), Err(NsfError::SgsHeatRatio { .. })));
    assert_eq!(serial.state_hash(), hash);
    assert_eq!(serial.sgs_transport, transport);
    assert_eq!(serial.cumulative_ledger(), ledger);
    assert_eq!(serial.steps(), steps);
}
