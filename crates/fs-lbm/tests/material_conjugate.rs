//! Shared material data through existing flow, advection and interface kernels.
//! Raw transport constants are synthetic fixtures, not measured material claims.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{EnergyConfig, FvBoundary, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, Voxel, VoxelDomain, simple_flow, solve_energy};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::plan_scaling_from_material;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn material_conductivity_reaches_actual_solid_fluid_interface_flux() {
    // Unit cross section: two metres of solid then two metres of resting
    // fluid. The exact series resistance is 2/k_s + 2/k_f, including the
    // two half-cell boundary resistances. No film coefficient is prescribed.
    let mut domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain.set(1, 0, 0, Voxel::Solid(0));
    let solids = [SolidMaterial::new("synthetic slab", 2.0)];
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(350.0);
    faces[1] = ThermalFace::Temperature(300.0);
    let thermal = ThermalSetup::new(faces);
    let before = thermal.clone();
    let gate = CancelGate::new();
    for k in [0.5, 1.0, 4.0] {
        let source = FluidTransportState::try_new(1000.0, 4000.0, k, 0.001).unwrap();
        let bound = MaterialChtFluid::try_new(&source).unwrap();
        let result = bound.solve_forced_convection(&domain, &solids, &flow, &thermal,
            &EnergyConfig::default(), &gate).unwrap();
        let heat = 50.0 / (2.0 / 2.0 + 2.0 / k);
        near(result.energy().solid_to_fluid_heat_w(&domain), heat, 1e-9);
        let resistances = [0.5 / 2.0, 1.5 / 2.0, 2.0 / 2.0 + 0.5 / k, 2.0 / 2.0 + 1.5 / k];
        for (&temperature, resistance) in result.energy().temperature.iter().zip(resistances) {
            near(temperature, 350.0 - heat * resistance, 1e-10);
        }
        assert!(result.energy().report.balance.residual_w.abs() < 1e-8);
        assert_eq!(result.transport(), &source);
        assert_eq!(result.flow().velocity_m_s, vec![[0.0; 3]; 4]);
    }
    assert_eq!(thermal, before);
}

#[test]
fn gas_pressure_changes_heat_capacity_flux_and_heated_outflow() {
    let domain = VoxelDomain::new(8, 2, 1, 0.01).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [0.02, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    let config = SimpleConfig::new(faces);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let mut thermal = ThermalSetup::new(faces);
    thermal.add_uniform_power(&domain, 0.005, |_| true);
    let snapshot = thermal.clone();
    let gate = CancelGate::new();
    let mut outlet = Vec::new();
    let mut reynolds = Vec::new();
    for pressure in [100000.0, 200000.0] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap();
        let bound = MaterialChtFluid::try_new(&gas).unwrap();
        let result = bound.solve_forced_convection(&domain, &[], &config, &thermal,
            &EnergyConfig::default(), &gate).unwrap();
        assert!(core::ptr::eq(result.material(), &gas));
        assert!(result.energy().report.balance.residual_w.abs() < 1e-8);
        assert!(result.energy().temperature.iter().all(|t| t.is_finite() && *t >= 300.0 - 1e-8));
        outlet.push(result.energy().bulk_temperature_x(&domain, &result.flow().field, 7).unwrap());
        reynolds.push(result.flow().report.max_cell_reynolds);
        // Compare to the exact existing pipeline, including its turbulence
        // transfer, not a second implementation with independently set fluid.
        let manual_flow = simple_flow(&domain, bound.properties(), &config, &gate).unwrap();
        let mut manual_setup = thermal.clone();
        manual_setup.eddy_conductivity_w_m_k = manual_flow.eddy_conductivity(bound.properties());
        let manual_heat = solve_energy(&domain, bound.properties(), &[], &manual_flow.field,
            &manual_setup, &EnergyConfig::default(), &gate).unwrap();
        assert_eq!(result.flow(), &manual_flow);
        assert_eq!(result.energy(), &manual_heat);
    }
    assert!(outlet[1] < outlet[0] && outlet[1] > 300.0);
    near(reynolds[1], 2.0 * reynolds[0], 1e-5);
    assert_eq!(thermal, snapshot);
}

#[test]
fn lattice_and_conjugate_plans_share_one_snapshot_without_resampling() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1.2, 1000.0, 0.026, 1.8e-5)
        }
    }
    let source = Counted(std::cell::Cell::new(0));
    let erased: &dyn FluidTransport = &source;
    let lattice = plan_scaling_from_material(erased, 0.001, 0.001, 0.01).unwrap();
    let bound = lattice.conjugate_fluid().unwrap();
    assert_eq!(source.0.get(), 1);
    assert!(core::ptr::eq(bound.material(), erased));
    assert_eq!(bound.transport(), lattice.transport());
    near(bound.properties().prandtl(), lattice.transport().prandtl(), 1e-14);
    near(bound.properties().kinematic_viscosity_m2_s, lattice.transport().kinematic_viscosity_m2_s(), 1e-14);
}

#[test]
fn unrepresentable_capacity_and_conflicting_eddy_input_refuse() {
    // Factored transport ratios remain finite, but CHT explicitly forms rho cp.
    let huge = FluidTransportState::try_new(1e200, 1e200, 1e300, 1e200).unwrap();
    assert!(MaterialChtFluid::try_new(&huge).is_err());
    let source = FluidTransportState::try_new(1.2, 1000.0, 0.026, 1.8e-5).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 0.01).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    setup.eddy_conductivity_w_m_k = vec![0.0; 2];
    let before = setup.clone();
    assert!(bound.solve_forced_convection(&domain, &[], &SimpleConfig::new([FvBoundary::wall(); 6]),
        &setup, &EnergyConfig::default(), &CancelGate::new()).is_err());
    assert_eq!(setup, before);
}
