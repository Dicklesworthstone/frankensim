//! A generated heat/species field drives the real regional duct response.
use super::*;
use fs_material::gas_field::steady::{AxialGasCell, SteadyAxialGas, SteadyGasSolve};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[path = "../../../../render/schedule/reed/field_input/tests/support.rs"]
mod humid_source;

#[test]
fn condensation_changes_source_bound_acoustic_propagation_and_terminal_gas() {
    use fs_material::gas_field::AxialGasSampling;
    use fs_material::gas_field::steady::condensation::CondensingSteadyGasSolve;
    use fs_material::phase::condensation::WaterCondensationModel;
    use fs_material::thermochemical_gas::packs::from_packs;
    use fs_thermochem::water::SaturatedWaterState;

    let carrier = humid_source::packs("carrier", 0.028, 3.5, 200.0);
    let water = humid_source::packs("water", 0.01801528, 4.0, 200.0);
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.94, 0.06]).unwrap();
    let water_id = SpeciesId::new("water").unwrap();
    let source = from_packs(&[carrier.binding(), water.binding()], composition.clone(),
        Some(&water_id)).unwrap();
    let phase = WaterCondensationModel::new(source.clone(), water_id, 300.0).unwrap();
    let count = 8;
    let dx = LENGTH / count as f64;
    let problem = SteadyAxialGas {
        source: source.clone(), source_conditions: QueryPoint::new(),
        static_pressure_pa: 101325.0, phase: ThermochemicalPhase::LiquidWaterIapwsLimit,
        inlet: PrimitiveGasState { temperature_k: 310.0, composition: composition.clone() },
        inlet_mass_flow_kg_s: 5e-5,
        cells: (0..count).map(|index| AxialGasCell {
            end_m: LENGTH * (index + 1) as f64 / count as f64,
            volume_m3: core::f64::consts::PI * RADIUS.powi(2) * dx,
            wall_area_m2: TAU * RADIUS * dx,
            wall_conductance_w_k: 0.08 / count as f64,
            wall_temperature_k: 285.0, initial_temperature_k: 300.0,
        }).collect(), injections: vec![],
    };
    let controls = CondensingSteadyGasSolve {
        caloric: GasEnergySolve { minimum_temperature_k: 274.0,
            maximum_temperature_k: 323.0, energy_tolerance_j_kg: 1e-8, maximum_iterations: 100 },
        maximum_temperature_trials_per_cell: 100, steady_power_tolerance_w: 1e-9,
        mass_flow_tolerance_kg_s: 1e-17, maximum_mean_mach: 0.05,
    };
    let solution = problem.solve_condensing(&phase, controls, || false).unwrap();
    assert_eq!(solution.field().sampling(), AxialGasSampling::CellwiseConstant);
    assert!(solution.total_condensed_mass_flow_kg_s > 0.0);
    assert!(solution.cells.last().unwrap().gas_mass_flow_kg_s < problem.inlet_mass_flow_kg_s);
    assert!(solution.maximum_power_residual_bound_w() <= controls.steady_power_tolerance_w);
    for sample in 0..=100 {
        let resolved = solution.field().state_at(LENGTH * sample as f64 / 100.0).unwrap();
        let gas = resolved.state();
        let saturation = SaturatedWaterState::at_temperature(Temperature::new(gas.temperature))
            .unwrap().pressure().value();
        assert!(gas.water_mole_fraction * gas.pressure <= saturation * (1.0 + 1e-13));
    }

    let duct = RegionalDuct::from_gas_field(&cylinder(), solution.field(), 64).unwrap();
    let continuous = RegionalDuct::from_continuous_gas_field(
        &cylinder(), solution.field(), fs_duct::RiccatiSteps(4)).unwrap();
    assert_eq!(duct.inlet_gas(), solution.field().state_at(0.0).unwrap().state());
    assert_eq!(duct.exterior_gas(), solution.field().state_at(LENGTH).unwrap().state());
    assert_eq!(continuous.exterior_gas(), duct.exterior_gas());
    assert!(duct.exterior_gas().water_mole_fraction < duct.inlet_gas().water_mole_fraction);
    let hot = AxialGasField::new(vec![0.0, LENGTH], 101325.0,
        vec![problem.inlet.clone(); 2], source.clone(), QueryPoint::new(), problem.phase).unwrap();
    let hot_duct = RegionalDuct::from_gas_field(&cylinder(), &hot, 64).unwrap();
    assert!(resonance(&duct, 1) < resonance(&hot_duct, 1));

    // A deliberately frozen counterfactual isolates composition from cooling:
    // it has the solved temperatures but keeps all of the original vapor.
    let frozen = AxialGasField::new_piecewise_constant(solution.field().stations_m().to_vec(),
        101325.0, solution.field().primitive_states().iter().map(|state| PrimitiveGasState {
            temperature_k: state.temperature_k, composition: composition.clone(),
        }).collect(), source, QueryPoint::new(), ThermochemicalPhase::FrozenGas).unwrap();
    let frozen_duct = RegionalDuct::from_gas_field(&cylinder(), &frozen, 64).unwrap();
    assert_eq!(duct.exterior_gas().temperature, frozen_duct.exterior_gas().temperature);
    assert_ne!(duct.exterior_gas().density, frozen_duct.exterior_gas().density);
    assert_ne!(duct.exterior_gas().sound_speed, frozen_duct.exterior_gas().sound_speed);
    for frequency in [150.0, 350.0, 900.0] {
        let omega = TAU * frequency;
        let actual = duct.response(omega, LossModel::Bessel, Termination::UnflangedOpen, None).unwrap();
        let integrated = continuous.response(omega, LossModel::Bessel,
            Termination::UnflangedOpen, None).unwrap();
        let counterfactual = frozen_duct.response(omega, LossModel::Bessel,
            Termination::UnflangedOpen, None).unwrap();
        assert!(actual.impedance.re.is_finite() && actual.impedance.re >= 0.0);
        assert!(integrated.impedance.re.is_finite() && integrated.impedance.re >= 0.0);
        assert!(actual.impedance.im.is_finite() && integrated.impedance.im.is_finite());
        assert_ne!(actual.impedance, counterfactual.impedance);
        assert_ne!(actual.p_mouth_over_u_in, counterfactual.p_mouth_over_u_in);
    }
}

#[test]
fn played_wall_cooling_changes_bore_registers_and_physical_termination() {
    let count = 16;
    let dx = LENGTH / count as f64;
    let problem = SteadyAxialGas {
        source: gas_model(), source_conditions: QueryPoint::new(),
        static_pressure_pa: 101325.0, phase: ThermochemicalPhase::FrozenGas,
        inlet: PrimitiveGasState { temperature_k: 300.0,
            composition: Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap() },
        inlet_mass_flow_kg_s: 0.00005,
        cells: (0..count).map(|index| AxialGasCell {
            end_m: LENGTH * (index + 1) as f64 / count as f64,
            volume_m3: core::f64::consts::PI * RADIUS.powi(2) * dx,
            wall_area_m2: TAU * RADIUS * dx,
            wall_conductance_w_k: 0.05 / count as f64,
            wall_temperature_k: 299.0, initial_temperature_k: 299.0,
        }).collect(), injections: vec![],
    };
    let controls = SteadyGasSolve {
        caloric: GasEnergySolve { minimum_temperature_k: 290.0,
            maximum_temperature_k: 310.0, energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 },
        maximum_dt_s: 0.05, maximum_relaxations: 2000, maximum_substeps: 100000,
        steady_power_tolerance_w: 1e-9, energy_tolerance_j: 1e-8, maximum_mean_mach: 0.05,
    };
    let solution = problem.solve(controls, || false).unwrap();
    let duct = RegionalDuct::from_gas_field(&cylinder(), solution.field(), 64).unwrap();
    let uniform = RegionalDuct::from_gas_field(&cylinder(), &field(300.0), 64).unwrap();
    assert_eq!(duct.inlet_gas(), solution.field().state_at(0.0).unwrap().state());
    assert_eq!(duct.exterior_gas(), solution.field().state_at(LENGTH).unwrap().state());
    assert!(duct.exterior_gas().temperature < 299.5);
    let f1 = resonance(&duct, 1);
    let f3 = resonance(&duct, 3);
    assert!(f1 < resonance(&uniform, 1));
    let register_change_cents = 1200.0 * (f3 / (3.0 * f1)).log2();
    assert!(register_change_cents > 0.01 && register_change_cents < 1.0,
        "register shift {register_change_cents}");
    for frequency in [150.0, 350.0, 900.0] {
        let actual = duct.response(TAU * frequency, LossModel::Bessel,
            Termination::UnflangedOpen, None).unwrap();
        let original = uniform.response(TAU * frequency, LossModel::Bessel,
            Termination::UnflangedOpen, None).unwrap();
        assert!(actual.impedance.re.is_finite() && actual.impedance.re >= 0.0);
        assert!(actual.impedance.im.is_finite());
        assert_ne!(actual.impedance, original.impedance);
        assert_ne!(actual.p_mouth_over_u_in, original.p_mouth_over_u_in);
    }
}
