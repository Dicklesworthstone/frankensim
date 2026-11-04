//! Physical radiation closure on the existing prescribed-flow energy owner.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, EnergyConfig, FlowFace, FlowField, FluidProperties, RadiationConfig,
    STEFAN_BOLTZMANN, SolidMaterial, SurfaceExchange, ThermalFace, ThermalSetup,
    Voxel, VoxelDomain, solve_energy_radiating,
};

const N: usize = 4;
const DX: f64 = 0.125;
const AMBIENT: f64 = 300.0;
const INLET: f64 = 1300.0;
const POWER: f64 = 500.0;

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1e-12,
        kinematic_viscosity_m2_s: 1e-3,
    }
}

fn prescribed_plate_closes_physical_radiation(surface_exchange: bool) {
    let gate = CancelGate::new();
    let domain = VoxelDomain::from_fn(N, N, N, DX, |p| {
        if p[2] < DX { Voxel::Solid(0) } else { Voxel::Fluid }
    }).unwrap();
    let solids = [SolidMaterial::new("heated plate", 200.0)];
    let flow = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed, FlowFace::Free, FlowFace::Wall,
            FlowFace::Wall, FlowFace::Wall, FlowFace::Wall,
        ],
        |axis, _| if axis == 0 { 0.125 } else { 0.0 },
    ).unwrap();
    assert!(flow.max_flux() > 0.0);
    assert_eq!(flow.max_divergence(&domain), 0.0);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: INLET };
    faces[1] = ThermalFace::Outflow { backflow_temperature: INLET };
    let mut setup = ThermalSetup::new(faces);
    setup.add_uniform_power(&domain, POWER, |p| p[2] < DX);

    // This actual throughflow establishes a wider temperature span than the
    // plate's radiative heating. Its tiny conductivity isolates the independent
    // Stefan-Boltzmann limit while retaining native advective transport.
    let emissivity = 0.8;
    let area = (N * N) as f64 * DX * DX;
    let coefficient = emissivity * STEFAN_BOLTZMANN * area;
    let seed = AMBIENT + POWER / (4.0 * coefficient * AMBIENT.powi(3));
    let first = 0.75 * seed
        + (POWER / coefficient + AMBIENT.powi(4)) / (4.0 * seed.powi(3));
    let first_power = coefficient * (first.powi(4) - AMBIENT.powi(4));
    let loose_tolerance = 0.2;
    assert!((seed - first).abs() / (INLET - first) < loose_tolerance);
    assert!((first_power - POWER).abs() / first_power > loose_tolerance);

    let mut surroundings = [Some(AMBIENT); 6];
    surroundings[4] = None;
    let mut radiation = RadiationConfig::new(vec![emissivity], surroundings);
    radiation.surface_exchange = surface_exchange;
    radiation.patch_size = 1;
    radiation.rays_per_face = 32;
    radiation.tolerance = loose_tolerance;
    radiation.max_iterations = 1;
    match solve_energy_radiating(
        &domain, &fluid(), &solids, &flow, &setup,
        &EnergyConfig::default(), &radiation, &gate,
    ) {
        Err(ChtError::SolverNotConverged {
            system, iterations, relative_residual, tolerance,
        }) => {
            assert_eq!(system, "radiation");
            assert_eq!(iterations, 1);
            assert_eq!(tolerance, loose_tolerance);
            assert!(relative_residual > tolerance);
        }
        other => panic!("a passing global temperature gate cannot waive physical heat: {other:?}"),
    }

    radiation.tolerance = 1e-10;
    radiation.max_iterations = 50;
    let (solution, report) = solve_energy_radiating(
        &domain, &fluid(), &solids, &flow, &setup,
        &EnergyConfig::default(), &radiation, &gate,
    ).unwrap();
    let exact = (POWER / coefficient + AMBIENT.powi(4)).powf(0.25);
    let mut physical_radiation_w = 0.0;
    for cell in 0..domain.cell_count() {
        if domain.coords(cell)[2] == 0 {
            let t = solution.temperature[cell];
            assert!((t - exact).abs() < 2e-6, "{t} versus {exact}");
            physical_radiation_w += emissivity * STEFAN_BOLTZMANN * DX * DX
                * (t.powi(4) - AMBIENT.powi(4));
        }
    }
    assert!(solution.report.max_cell_peclet > 1e6);
    assert!(report.iterations > 1);
    assert!(report.temperature_change <= radiation.tolerance);
    assert!(report.relative_power_defect <= radiation.tolerance);
    assert!(report.maximum_cell_power_defect_w < 1e-7 * POWER);
    assert!((physical_radiation_w - POWER).abs() < 1e-7 * POWER);
    assert!((report.radiated_w - physical_radiation_w).abs() < 1e-10 * POWER);
    assert!(
        (report.net_outflow_correction_w
            - (physical_radiation_w - solution.report.balance.sink_outflow_w)).abs()
            < 1e-10 * POWER,
    );
    assert!(report.physical_energy_residual_w.abs() < 1e-7 * POWER);
    assert_eq!(
        report.physical_energy_residual_w,
        solution.report.balance.residual_w - report.net_outflow_correction_w,
    );
}

#[test]
fn escape_radiation_closes_physical_heat_on_prescribed_flow() {
    prescribed_plate_closes_physical_radiation(false);
}

#[test]
fn surface_exchange_closes_physical_heat_on_prescribed_flow() {
    prescribed_plate_closes_physical_radiation(true);
}

#[test]
fn fixed_temperature_radiation_correction_belongs_to_the_thermostat() {
    let gate = CancelGate::new();
    let domain = VoxelDomain::from_fn(N, N, N, DX, |p| {
        if p[2] < DX || p[2] > (N - 1) as f64 * DX {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    }).unwrap();
    let solids = [SolidMaterial::new("plates", 200.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.add_uniform_power(&domain, POWER, |p| p[2] < DX);
    let is_fixed = |cell| domain.coords(cell)[2] == N - 1;
    for cell in 0..domain.cell_count() {
        if is_fixed(cell) {
            setup.fixed_temperature.push((cell, AMBIENT));
        }
    }
    let mut radiation = RadiationConfig::new(vec![1.0], [None; 6]);
    radiation.patch_size = 1;
    radiation.rays_per_face = 64;
    // Keep a visible constitutive defect to test who owns its heat correction;
    // a fully converged tangent could conceal incorrect thermostat accounting.
    radiation.tolerance = 2.0;
    radiation.max_iterations = 1;
    let (solution, report) = solve_energy_radiating(
        &domain, &fluid(), &solids, &FlowField::quiescent(&domain), &setup,
        &EnergyConfig::default(), &radiation, &gate,
    ).unwrap();
    assert_eq!(report.iterations, 1);
    assert_eq!(report.radiated_w, 0.0);

    // Independently evaluate the heat supplied by the fixed upper plate at
    // the returned temperatures: its actual conductive links plus its original
    // radiation law. With black surfaces the irradiation is an explicit sum.
    let exchange = SurfaceExchange::build(&domain, &solids, &radiation, &gate).unwrap();
    let irradiation = exchange.irradiation(&solution.temperature);
    let physical_radiation_w = exchange.net_emission(&solution.temperature, &irradiation);
    assert!(physical_radiation_w.abs() < 1e-10 * POWER);
    let bath_radiation_w: f64 = exchange.sinks(&solution.temperature, &irradiation)
        .iter()
        .filter(|sink| is_fixed(sink.cell))
        .map(|sink| {
            sink.conductance_w_k * (solution.temperature[sink.cell] - sink.temperature_k)
        })
        .sum();
    let bath_conduction_w: f64 = solution.solid_fluid_links.iter()
        .filter(|&&(solid, _, _)| is_fixed(solid))
        .map(|&(solid, fluid, conductance)| {
            conductance * (solution.temperature[solid] - solution.temperature[fluid])
        })
        .sum();
    let physical_bath_w = bath_radiation_w + bath_conduction_w;
    let balance = solution.report.balance;
    let bath_correction_w = physical_bath_w - balance.fixed_cell_injection_w;
    assert!(bath_correction_w.abs() > 1e-4 * POWER);
    let total_radiation_correction_w = physical_radiation_w - balance.sink_outflow_w;
    assert!(
        (report.net_outflow_correction_w
            - (total_radiation_correction_w - bath_correction_w)).abs()
            < 1e-8 * POWER,
    );
    let physical_residual_w = balance.source_w + physical_bath_w
        - balance.boundary_outflow_w - physical_radiation_w;
    assert!(
        (report.physical_energy_residual_w - physical_residual_w).abs() < 1e-8 * POWER,
        "physical residual {} versus independently evaluated {physical_residual_w}",
        report.physical_energy_residual_w,
    );
}
