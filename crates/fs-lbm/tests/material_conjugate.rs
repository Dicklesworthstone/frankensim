//! Shared material data through existing flow, advection and interface kernels.
//! Raw transport constants are synthetic fixtures, not measured material claims.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{EnergyConfig, FvBoundary, FvBuoyancyConfig, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, Voxel, VoxelDomain, fv_natural_convection, simple_flow, solve_energy};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_lbm::material::conjugate::buoyant::{MaterialBuoyancyConfig, MaterialBuoyancyError};
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_lbm::material::plan_scaling_from_material;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};

#[path = "material_conjugate/source_enthalpy.rs"]
mod source_enthalpy;

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

fn natural_material_controls() -> MaterialBuoyancyConfig {
    MaterialBuoyancyConfig {
        boussinesq: BoussinesqConfig {
            gravity_m_s2: [-9.81, 0.0, 0.0],
            max_relative_density_departure: 0.05,
            max_relative_specific_heat_departure: 0.0,
        },
        maximum_passes: 100,
        temperature_tolerance_k: 1e-10,
        transport_relative_tolerance: 1e-10,
        force_tolerance_n_m3: 1e-10,
        temperature_relaxation: 1.0,
    }
}

#[test]
fn natural_material_closes_nonlinear_conductivity_at_actual_thermal_endpoint() {
    // G1: one unit cube, six faces at 300 K. Q=12*k(T)*(T-300),
    // k=1+0.1*(T-300), hence Q=67.2 W has the independent root T=304 K.
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut thermal = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    thermal.power_w = vec![67.2];
    let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &[300.0], natural_material_controls(), &CancelGate::new(),
        |_, t, _| FluidTransportState::try_new(1.0, 2.0, 1.0 + 0.1*(t-300.0), 1.0/(1.0+0.05*(t-300.0))))
        .unwrap();
    near(result.energy.temperature[0], 304.0, 1e-11);
    assert!(result.report.passes > 2);
    assert!(result.report.temperature_defect_k <= 1e-10);
    assert!(result.report.transport_relative_defect <= 1e-10);
    assert!(result.energy.report.balance.residual_w.abs() < 1e-8);
    let final_material = result.endpoint[0].material.resolved_transport();
    near(final_material.thermal_conductivity_w_m_k(), 1.4, 1e-11);
    near(final_material.dynamic_viscosity_pa_s(), 1.0/1.2, 1e-11);
    assert_eq!(result.endpoint[0].temperature_k, result.energy.temperature[0]);
    assert_eq!(result.flow.velocity_m_s, vec![[0.0; 3]]);
}

#[test]
fn natural_material_density_and_viscosity_drive_same_existing_porous_momentum() {
    // G1/G3: both ends are pressure openings, all cells thermostat-held at
    // 310 K. In this uniform Darcy column, q=(rho-rho0)*g*K/mu exactly.
    // Changing only source data includes ordinary, zero and negative expansion.
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-12;
    flow.pressure_tolerance = 1e-12;
    flow.resistances.push(fs_lbm::conjugate::FlowResistance::Volume {
        lo: [0, 0, 0], hi: [4, 1, 1], permeability_m2: [1.0; 3], inertial_per_m: [0.0; 3],
    });
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let mut thermal = ThermalSetup::new(faces);
    thermal.fixed_temperature = (0..4).map(|c| (c, 310.0)).collect();
    for (expansion, viscosity_scale) in [(0.0, 1.0), (0.001, 1.0), (-0.001, 1.0), (0.001, 2.0)] {
        let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
            &EnergyConfig::default(), &[300.0; 4], natural_material_controls(), &CancelGate::new(),
            |_, t, _| FluidTransportState::try_new(1.0 - expansion*(t-300.0), 2.0,
                1.0 + 0.01*(t-300.0), viscosity_scale*(1.0+0.1*(t-300.0))))
            .unwrap();
        let expected = 9.81 * expansion * 10.0 / (2.0 * viscosity_scale);
        near(result.flow.field.flux_x(2, 0, 0), expected, 1e-9);
        // The frozen-coefficient public entry must admit the SAME signed
        // density slope as the source-resolved entry. Here every cell is at
        // 310 K, so the local viscosity and conductivity are exactly uniform.
        let frozen = FluidTransportState::try_new(1.0, 2.0, 1.1, 2.0 * viscosity_scale).unwrap();
        let frozen = MaterialChtFluid::try_new(&frozen).unwrap();
        let frozen_result = fv_natural_convection(&domain, frozen.properties(), &[], &thermal,
            &FvBuoyancyConfig::new([-9.81, 0.0, 0.0], expansion, 300.0, flow.clone()),
            &CancelGate::new()).unwrap();
        near(frozen_result.flow.field.flux_x(2, 0, 0), expected, 1e-9);
        assert_eq!(frozen_result.energy.temperature, result.energy.temperature);
        assert!(frozen_result.energy.report.balance.residual_w.abs() < 1e-8);
        assert!(result.energy.report.balance.residual_w.abs() < 1e-8);
        assert_eq!(result.reference, reference);
        assert_eq!(result.energy.temperature, vec![310.0; 4]);
        for sample in result.applied {
            near(sample.force_n_m3[0], 9.81*expansion*10.0, 1e-12);
            near(sample.material.resolved_transport().dynamic_viscosity_pa_s(), 2.0*viscosity_scale, 1e-12);
            near(sample.material.resolved_transport().thermal_conductivity_w_m_k(), 1.1, 1e-12);
            assert_eq!(sample.material.reference_storage_transport().density_kg_m3(), 1.0);
        }
    }
}

#[test]
fn natural_material_endpoint_refusal_cancellation_and_budget_publish_no_partial_result() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut thermal = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    thermal.power_w = vec![48.0];
    let before = thermal.clone();
    let initial = [300.0];
    let error = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &initial, natural_material_controls(), &CancelGate::new(),
        |_, t, _| if t > 301.0 { Err(MaterialError::State { what: "synthetic source ends at 301 K".into() }) }
            else { Ok(source) });
    assert!(matches!(error, Err(MaterialBuoyancyError::Material { cell: 0, temperature_k: 304.0, .. })));
    let mut controls = natural_material_controls();
    controls.maximum_passes = 2;
    controls.temperature_relaxation = 1e-8;
    controls.temperature_tolerance_k = 1e-3;
    assert!(matches!(bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &initial, controls, &CancelGate::new(), |_, _, _| Ok(source)),
        Err(MaterialBuoyancyError::NotConverged(_))));
    let gate = CancelGate::new();
    assert!(matches!(bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &initial, controls, &gate, |_, _, _| { gate.request(); Ok(source) }),
        Err(MaterialBuoyancyError::Solver(fs_lbm::conjugate::ChtError::Cancelled))));
    assert_eq!(thermal, before);
    assert_eq!(initial, [300.0]);
}

#[test]
fn natural_material_force_jump_cannot_hide_behind_small_temperature_change() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let thermal = ThermalSetup::new([ThermalFace::Temperature(310.0); 6]);
    let mut controls = natural_material_controls();
    controls.maximum_passes = 1;
    controls.temperature_tolerance_k = 1e-3;
    let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &[310.0-1e-6], controls, &CancelGate::new(),
        |_, t, _| FluidTransportState::try_new(if t < 310.0-1e-7 { 1.0 } else { 0.99 }, 2.0, 1.0, 1.0));
    match result {
        Err(MaterialBuoyancyError::NotConverged(report)) => {
            assert!(report.temperature_defect_k < controls.temperature_tolerance_k);
            near(report.force_defect_n_m3, 0.0981, 1e-12);
        }
        _ => panic!("discontinuous source force must fail its independent nonlinear gate"),
    }
}

#[test]
fn natural_material_heated_open_column_starts_from_uniform_ambient() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.1, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-12;
    flow.pressure_tolerance = 1e-12;
    flow.resistances.push(fs_lbm::conjugate::FlowResistance::Volume {
        lo: [0, 0, 0], hi: [4, 1, 1], permeability_m2: [1.0; 3], inertial_per_m: [0.0; 3],
    });
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let mut thermal = ThermalSetup::new(faces);
    thermal.power_w = vec![0.0, 1.0, 0.0, 0.0];
    let mut controls = natural_material_controls();
    controls.temperature_relaxation = 0.5;
    let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
        &EnergyConfig::default(), &[300.0; 4], controls, &CancelGate::new(),
        |_, t, _| FluidTransportState::try_new(1.0-0.001*(t-300.0), 2.0, 0.1, 1.0))
        .unwrap();
    let flux = result.flow.field.flux_x(2, 0, 0);
    assert!(flux > 0.0);
    let outlet = result.energy.temperature[3];
    assert!(outlet > 300.0);
    // All heat leaves as the SAME conservative flow's sensible energy flux.
    near(2.0 * flux * (outlet-300.0), 1.0, 1e-7);
    assert!(result.energy.report.balance.residual_w.abs() < 1e-7);
}

#[test]
fn natural_material_low_conductivity_open_column_preserves_warm_and_pressure_driven_starts() {
    use fs_lbm::conjugate::{FacePatch, FanCurve, FlowResistance, InternalFan};

    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.001, 1.0).unwrap();
    let bound = MaterialChtFluid::try_new(&source).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-12;
    flow.pressure_tolerance = 1e-12;
    flow.resistances.push(FlowResistance::Volume {
        lo: [0, 0, 0], hi: [4, 1, 1], permeability_m2: [1.0; 3], inertial_per_m: [0.0; 3],
    });
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let mut thermal = ThermalSetup::new(faces);
    thermal.power_w = vec![0.0, 1.0, 0.0, 0.0];
    let mut controls = natural_material_controls();
    controls.temperature_relaxation = 0.5;
    // A zero-flow conductive seed reaches 1237.5 K, outside the original
    // source domain. Both admitted starts already carry enough actual flow
    // to close the heat balance without that fictitious hot material state.
    for fan_driven in [false, true] {
        let mut flow = flow.clone();
        let initial = if fan_driven {
            flow.internal_fans.push(InternalFan {
                patch: FacePatch { axis: 0, index: 2, lo: [0, 0, 0], hi: [4, 1, 1] },
                blows_positive: true,
                curve: FanCurve::new(&[(0.0, 0.4), (1.0, 0.4)]).unwrap(),
            });
            [300.0; 4]
        } else {
            [300.0, 308.25, 308.25, 308.25]
        };
        let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
            &EnergyConfig::default(), &initial, controls, &CancelGate::new(), |_, t, _| {
                if !(295.0..=320.0).contains(&t) {
                    return Err(MaterialError::State { what: "synthetic source supports only 295..320 K".into() });
                }
                let density = if fan_driven { 1.0 } else { 1.0 - 0.001 * (t - 300.0) };
                FluidTransportState::try_new(density, 2.0, 0.001, 1.0)
            }).unwrap();
        let flux = result.flow.field.flux_x(2, 0, 0);
        assert!(flux > 0.0);
        let outlet = result.energy.temperature[3];
        assert!((300.0..320.0).contains(&outlet));
        near(2.0 * flux * (outlet - 300.0), 1.0, 1e-7);
        assert!(result.energy.report.balance.residual_w.abs() < 1e-7);
        if fan_driven {
            near(flux, 0.1, 1e-9);
            assert!(result.applied.iter().all(|sample| sample.force_n_m3 == [0.0; 3]));
        }
    }
}

#[test]
fn natural_material_card_swap_refreshes_same_source_claims_at_thermal_endpoint() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
    use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
    use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection,
        SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
    use fs_qty::{Density, QuantitySpec, Temperature};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("T", absolute, 300.0).unwrap()
        .with("pressure", 100000.0).unwrap();
    let card = |k: f64| {
        let mut claims = ClaimSet::new();
        for (name, dims, values) in [
            (DENSITY_PROPERTY, Density::DIMS, [1.0, 0.96]),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, [2.0, 2.0]),
            (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, [k, 5.0*k]),
            (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, [1.0, 0.5]),
        ] {
            let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
                QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
            } else { QuantitySpec::dimensional(dims) };
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (340.0, values[1])], dims },
                validity: ValidityDomain::unconstrained().with_quantity("T", absolute, 300.0, 340.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance { source: "synthetic source-card swap, not measurements".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
        }
        MaterialCard::assemble(MaterialStateId { chemistry: "synthetic liquid".into(),
            phase: "liquid".into(), process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
    };
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut thermal = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    thermal.power_w = vec![67.2];
    let mut outcomes = Vec::new();
    for k in [1.0, 2.0] {
        let card = card(k);
        let law = LiquidTemperatureLaw::try_new(&card, &point, "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
        let resolve = |_, t, _: &FluidTransportState| law.at_temperature(Temperature::new(t))
            .map_err(|e| MaterialError::State { what: e.to_string() });
        let result = bound.solve_natural_convection(&domain, &[], &flow, &thermal,
            &EnergyConfig::default(), &[300.0], natural_material_controls(), &CancelGate::new(), resolve).unwrap();
        // Independent scalar balance: 0.1*theta^2 + theta = Q/(12*k).
        let expected = 300.0 + ((1.0_f64+0.4*67.2/(12.0*k)).sqrt()-1.0)/0.2;
        near(result.energy.temperature[0], expected, 1e-11);
        let endpoint = &result.endpoint[0];
        let state = endpoint.material.material().material();
        assert_eq!(state.card_identity(), card.content_hash());
        for name in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY,
            THERMAL_CONDUCTIVITY_PROPERTY, DYNAMIC_VISCOSITY_PROPERTY] {
            let receipt = &state.property(name).unwrap().answer().receipt;
            assert_eq!(receipt.selected,
                law.reference_state().material().property(name).unwrap().answer().receipt.selected);
            assert!(receipt.query_point.contains(&("T".into(), endpoint.temperature_k)));
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
        outcomes.push((card.content_hash(), result.energy.temperature[0]));
        let mut outside = thermal.clone();
        outside.power_w = vec![1200.0];
        assert!(matches!(bound.solve_natural_convection(&domain, &[], &flow, &outside,
            &EnergyConfig::default(), &[300.0], natural_material_controls(), &CancelGate::new(), resolve),
            Err(MaterialBuoyancyError::Material { cell: 0, .. })));
    }
    assert_ne!(outcomes[0].0, outcomes[1].0);
    assert!(outcomes[1].1 < outcomes[0].1);
}

mod solid_feedback {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_lbm::material::conjugate::solid::SolidConductivityLaw;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
    use fs_qty::{QuantitySpec, Temperature};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    fn quantity(kind: QuantityKind) -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
    }
    fn card(slopes: [f64; 3], upper: f64, sharp: bool) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, slope) in ["kx", "ky", "kz"].into_iter().zip(slopes) {
            let knots = if sharp {
                vec![(300.0, 1.0), (310.0-1e-6, 1.0), (310.0, 2.0), (upper, 2.0)]
            } else { vec![(300.0, 1.0), (upper, 1.0+slope*(upper-300.0))] };
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity(QuantityKind::ThermalConductivity)),
                value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots, dims: THERMAL_CONDUCTIVITY_DIMS },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0, upper)
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0, 100000.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: Provenance { source: "synthetic solid conductivity oracle, not measurements".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
        }
        MaterialCard::assemble(MaterialStateId { chemistry: "synthetic thermal solid".into(),
            phase: "solid".into(), process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
    }
    fn law(card: &MaterialCard) -> SolidConductivityLaw<'_> {
        let point = QueryPoint::new().with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
            .unwrap().with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0).unwrap();
        let requirements = ["kx", "ky", "kz"].map(|name| ScalarPropertyRequirement::try_with_quantity(
            name, quantity(QuantityKind::ThermalConductivity), ScalarAdmissibility::StrictlyPositive).unwrap());
        let source = resolve_material_state_point(card, &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        SolidConductivityLaw::try_new(card, source, "T", ["kx", "ky", "kz"]).unwrap()
    }
    fn setup(axis: usize) -> (VoxelDomain, ThermalSetup) {
        let mut cells = [1; 3]; cells[axis] = 2;
        let mut domain = VoxelDomain::new(cells[0], cells[1], cells[2], 1.0).unwrap();
        domain.set(0, 0, 0, Voxel::Solid(0));
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[2*axis] = ThermalFace::Temperature(300.0);
        faces[2*axis+1] = ThermalFace::Temperature(400.0);
        (domain, ThermalSetup::new(faces))
    }

    #[test]
    fn g1_g3_solid_material_and_axis_changes_reach_same_conjugate_heat_solution() {
        let fluid = FluidTransportState::try_new(1.0, 2.0, 2.0, 1.0).unwrap();
        let bound = MaterialChtFluid::try_new(&fluid).unwrap();
        let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
        let mut controls = natural_material_controls();
        controls.boussinesq.gravity_m_s2 = [0.0; 3];
        for slopes in [[0.0; 3], [0.01, 0.02, 0.03]] {
            let card = card(slopes, 400.0, false);
            let source = law(&card);
            for axis in 0..3 {
                let (domain, thermal) = setup(axis);
                let result = bound.solve_natural_convection_with_material_solids(
                    &domain, std::slice::from_ref(&source), &flow, &thermal,
                    &EnergyConfig::default(), &[300.0; 2], controls, &CancelGate::new(),
                    |_, _, _| Ok(fluid)).unwrap();
                // Independent two-layer FV balance: theta*(3+a*theta)=100.
                // The solid and fluid halves carry exactly the same heat.
                let a = slopes[axis];
                let theta = if a == 0.0 { 100.0/3.0 }
                    else { 200.0 / (3.0+(9.0+400.0*a).sqrt()) };
                let heat = 2.0*(1.0+a*theta)*theta;
                near(result.energy.temperature[0], 300.0+theta, 1e-11);
                near(result.energy.temperature[1], 400.0-heat/4.0, 1e-11);
                near(result.energy.solid_to_fluid_heat_w(&domain), -heat, 1e-9);
                assert!(result.energy.report.balance.residual_w.abs() < 1e-8);
                assert!(result.flow.velocity_m_s.iter().all(|v| *v == [0.0; 3]));
                let endpoint = &result.solid_endpoint[0];
                assert_eq!(endpoint.cell, 0);
                assert_eq!(endpoint.material_index, 0);
                assert_eq!(endpoint.temperature_k, result.energy.temperature[0]);
                assert_eq!(endpoint.material.card_identity(), card.content_hash());
                for name in ["kx", "ky", "kz"] {
                    let receipt = &endpoint.material.property(name).unwrap().answer().receipt;
                    assert_eq!(receipt.selected, source.reference().property(name).unwrap().answer().receipt.selected);
                    assert!(receipt.query_point.contains(&("T".into(), endpoint.temperature_k)));
                    assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
                }
                near(endpoint.conductivity_w_m_k[axis], 1.0+a*theta, 1e-11);
                assert!(result.report.transport_relative_defect <= controls.transport_relative_tolerance);
            }
        }
    }

    #[test]
    fn g3_solid_only_substitution_changes_actual_fluid_buoyancy_and_circulation() {
        // The existing differentially heated FV cavity, with a conducting left
        // wall column. This is a coupling regression, not a mesh-convergence or
        // experimental claim. The fluid law, geometry and forcing stay fixed.
        let n = 6;
        let mut domain = VoxelDomain::new(n, n, 1, 1.0/n as f64).unwrap();
        for y in 0..n { domain.set(0, y, 0, Voxel::Solid(0)); }
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Temperature(310.0);
        faces[1] = ThermalFace::Temperature(300.0);
        let thermal = ThermalSetup::new(faces);
        let mut flow = SimpleConfig::new([FvBoundary::wall(), FvBoundary::wall(),
            FvBoundary::wall(), FvBoundary::wall(), FvBoundary::Symmetry, FvBoundary::Symmetry]);
        flow.tolerance = 1e-10;
        flow.pressure_tolerance = 1e-12;
        let fluid = FluidTransportState::try_new(1.0, 2.0, 1.0, 1.0).unwrap();
        let bound = MaterialChtFluid::try_new(&fluid).unwrap();
        let mut controls = natural_material_controls();
        controls.boussinesq.gravity_m_s2 = [0.0, -9.81, 0.0];
        controls.temperature_tolerance_k = 1e-8;
        controls.transport_relative_tolerance = 1e-9;
        controls.force_tolerance_n_m3 = 1e-9;
        let mut results = Vec::new();
        for slope in [0.0, 1.0] {
            let card = card([slope; 3], 400.0, false);
            let result = bound.solve_natural_convection_with_material_solids(
                &domain, &[law(&card)], &flow, &thermal, &EnergyConfig::default(),
                &vec![300.0; domain.cell_count()], controls, &CancelGate::new(),
                |_, t, _| FluidTransportState::try_new(1.0-0.001*(t-300.0), 2.0, 1.0, 1.0))
                .unwrap();
            assert!(result.energy.report.balance.residual_w.abs() < 1e-8);
            let peak = result.flow.velocity_m_s.iter().flat_map(|v| v.iter())
                .copied().map(f64::abs).fold(0.0, f64::max);
            assert!(peak > 1e-8);
            for sample in &result.endpoint {
                near(sample.force_n_m3[1], 9.81*0.001*(sample.temperature_k-300.0), 1e-12);
            }
            for sample in &result.solid_endpoint {
                assert_eq!(sample.temperature_k, result.energy.temperature[sample.cell]);
                near(sample.conductivity_w_m_k[0], 1.0+slope*(sample.temperature_k-300.0), 1e-12);
            }
            results.push(result);
        }
        let velocity_change = results[0].flow.velocity_m_s.iter().zip(&results[1].flow.velocity_m_s)
            .flat_map(|(a,b)| a.iter().zip(b).map(|(a,b)| (a-b).abs())).fold(0.0, f64::max);
        let force_change = results[0].endpoint.iter().zip(&results[1].endpoint)
            .map(|(a,b)| (a.force_n_m3[1]-b.force_n_m3[1]).abs()).fold(0.0, f64::max);
        assert!(force_change > 1e-4);
        assert!(velocity_change > 1e-7);
    }

    #[test]
    fn g0_solid_endpoint_refusal_and_cancel_preserve_original_inputs() {
        let card = card([0.01; 3], 320.0, false);
        let source = law(&card);
        let original = source.reference().clone();
        let (domain, thermal) = setup(0);
        let original_thermal = thermal.clone();
        let initial = [300.0; 2];
        let fluid = FluidTransportState::try_new(1.0, 2.0, 2.0, 1.0).unwrap();
        let bound = MaterialChtFluid::try_new(&fluid).unwrap();
        let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
        let gate = CancelGate::new();
        assert!(matches!(bound.solve_natural_convection_with_material_solids(
            &domain, std::slice::from_ref(&source), &flow, &thermal, &EnergyConfig::default(),
            &initial, natural_material_controls(), &gate, |_, _, _| Ok(fluid)),
            Err(MaterialBuoyancyError::SolidMaterial { cell: 0, .. })));
        assert!(matches!(bound.solve_natural_convection_with_material_solids(
            &domain, std::slice::from_ref(&source), &flow, &thermal, &EnergyConfig::default(),
            &initial, natural_material_controls(), &gate, |_, _, _| { gate.request(); Ok(fluid) }),
            Err(MaterialBuoyancyError::Solver(fs_lbm::conjugate::ChtError::Cancelled))));
        assert_eq!(source.reference(), &original);
        assert_eq!(thermal, original_thermal);
        assert_eq!(initial, [300.0; 2]);
    }

    #[test]
    fn g0_solid_conductivity_jump_has_its_own_nonlinear_acceptance_gate() {
        let card = card([0.0; 3], 400.0, true);
        let source = law(&card);
        let (domain, _) = setup(0);
        let thermal = ThermalSetup::new([ThermalFace::Temperature(310.0); 6]);
        let fluid = FluidTransportState::try_new(1.0, 2.0, 2.0, 1.0).unwrap();
        let bound = MaterialChtFluid::try_new(&fluid).unwrap();
        let mut controls = natural_material_controls();
        controls.maximum_passes = 1;
        controls.temperature_tolerance_k = 1e-3;
        let result = bound.solve_natural_convection_with_material_solids(
            &domain, &[source], &SimpleConfig::new([FvBoundary::wall(); 6]), &thermal,
            &EnergyConfig::default(), &[310.0-1e-6, 310.0], controls, &CancelGate::new(),
            |_, _, _| Ok(fluid));
        match result {
            Err(MaterialBuoyancyError::NotConverged(report)) => {
                assert!(report.temperature_defect_k < controls.temperature_tolerance_k);
                assert!(report.transport_relative_defect > 0.49);
            }
            _ => panic!("a solid material jump must fail even when temperature already agrees"),
        }
    }
}
