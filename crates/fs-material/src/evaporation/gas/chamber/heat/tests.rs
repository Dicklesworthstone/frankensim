//! Independent constant-capacity star-network solutions and original nonlinear
//! caloric/phase regressions. Coefficients here are explicit manufactured data.
use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use crate::thermochemical_gas::tests::{species_with_energy_reference, transport};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.018015;
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / M;
const CP: f64 = 4.0 * R;
const LATENT: f64 = 200000.0;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn constant_source() -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", M, 4.0, 0.0, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&SpeciesId::new("B").unwrap()),
    )
    .unwrap()
}
fn constant_temperature(h: f64) -> Result<f64, fs_thermal::evaporation::EvaporationError> {
    Ok((h + LATENT) / CP)
}
fn constant_liquid(t: f64, mass: f64, area: f64) -> LiquidSurfaceState {
    LiquidSurfaceState {
        mass_kg: mass,
        area_m2: area,
        specific_enthalpy_j_kg: CP * t - LATENT,
    }
}
fn heat_controls(constant: bool) -> GasFilmHeatSolve {
    let h = |t| {
        if constant {
            CP * t - LATENT
        } else {
            liquid(t).specific_enthalpy_j_kg
        }
    };
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [h(200.0), h(900.0)],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn phase_controls() -> GasLiquidChamberSolve {
    let mut surface = controls();
    surface.maximum_temperature_increment_k = 200.0;
    GasLiquidChamberSolve {
        surface,
        maximum_surfaces: 8,
    }
}
fn outer() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}
fn surface<'a>(
    id: u64,
    binding: &'a super::super::GasLiquidBinding,
    liquid: LiquidSurfaceState,
    coefficient: f64,
    constant: bool,
) -> GasLiquidChamberHeatSurface<'a> {
    GasLiquidChamberHeatSurface {
        surface: GasLiquidChamberSurface {
            id,
            law: &Phase,
            binding,
            liquid,
            temperature_from_enthalpy: if constant {
                &constant_temperature
            } else {
                &temperature
            },
            mass_transfer_coefficient_m_s: 0.0,
        },
        heat_transfer_coefficient_w_m2_k: coefficient,
        heat_controls: heat_controls(constant),
    }
}
fn constant_gas_capacity(gas: &ClosedGasParcel) -> f64 {
    let molar_mass = 0.99 * 0.02897 + 0.01 * M;
    gas.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass
        * (0.99 * 2.5 + 0.01 * 3.0)
}

#[test]
fn one_shared_capacity_matches_closed_form_with_initial_equality_and_heat_reversal() {
    let source = constant_source();
    let binding = binding(0.0);
    let initial = gas(&source, 300.0, 1000.0, 0.002);
    let dt = 1.0;
    let inputs = [
        surface(
            3,
            &binding,
            constant_liquid(300.0, 0.001, 0.001),
            1000.0,
            true,
        ),
        surface(
            5,
            &binding,
            constant_liquid(295.0, 0.001, 0.001),
            1000.0,
            true,
        ),
        surface(9, &binding, constant_liquid(270.0, 0.1, 0.01), 1000.0, true),
    ];
    let mut numerator = constant_gas_capacity(&initial) * 300.0;
    let mut denominator = constant_gas_capacity(&initial);
    let mut weights = Vec::new();
    for input in &inputs {
        let capacity = input.surface.liquid.mass_kg * CP;
        let k = input.heat_transfer_coefficient_w_m2_k * input.surface.liquid.area_m2 * dt;
        let weight = capacity * k / (capacity + k);
        let t = constant_temperature(input.surface.liquid.specific_enthalpy_j_kg).unwrap();
        numerator += weight * t;
        denominator += weight;
        weights.push(weight);
    }
    let expected_gas = numerator / denominator;
    assert!(expected_gas < 295.0);
    let candidate = advance_chamber_with_gas_and_heat(
        &inputs,
        &initial,
        &source,
        dt,
        phase_controls(),
        outer(),
        || false,
    )
    .unwrap();
    near(
        candidate.gas().state().state().temperature,
        expected_gas,
        2e-8,
    );
    assert!(
        candidate.heat_surfaces()[0].heat_into_liquid_j < 0.0,
        "an initially equal-temperature film must still exchange heat"
    );
    assert!(
        candidate.heat_surfaces()[1].heat_into_liquid_j < 0.0,
        "other films can reverse the initially warmer-gas heat direction"
    );
    assert!(candidate.heat_surfaces()[2].heat_into_liquid_j > 0.0);
    for ((input, report), weight) in inputs.iter().zip(candidate.heat_surfaces()).zip(weights) {
        let initial_t = constant_temperature(input.surface.liquid.specific_enthalpy_j_kg).unwrap();
        let expected_q = weight * (expected_gas - initial_t);
        near(report.heat_into_liquid_j, expected_q, 2e-8);
        near(
            report.temperature_k[1],
            initial_t + expected_q / (input.surface.liquid.mass_kg * CP),
            2e-8,
        );
        near(
            report.temperature_equation_residual_k,
            0.0,
            input.heat_controls.temperature_tolerance_k,
        );
        assert_eq!(
            report.phase_initial_liquid.mass_kg,
            input.surface.liquid.mass_kg
        );
        assert_eq!(
            report.heat_into_liquid_j,
            input.surface.liquid.mass_kg
                * (report.phase_initial_liquid.specific_enthalpy_j_kg
                    - input.surface.liquid.specific_enthalpy_j_kg)
        );
    }
    assert_eq!(candidate.phase().step().transferred_mass_kg, 0.0);
    assert_eq!(candidate.phase().step().gas_internal_energy_change_j, 0.0);
    assert_eq!(candidate.gas().mass_kg(), initial.mass_kg());
    assert_eq!(candidate.gas().model(), initial.model());
    assert!(
        candidate.step().sum_absolute_energy_residual_j
            <= phase_controls().surface.gas.energy_tolerance_j
    );
    let reversed = advance_chamber_with_gas_and_heat(
        &[inputs[2], inputs[0], inputs[1]],
        &initial,
        &source,
        dt,
        phase_controls(),
        outer(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate, reversed);
}

#[test]
fn zero_heat_is_exact_phase_path_and_one_active_film_preserves_existing_heat_owner() {
    let source = catalog(-1e7);
    let binding = binding(-1e7);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let mut a = surface(2, &binding, liquid(300.0), 0.0, false);
    let mut b = surface(7, &binding, liquid(305.0), 0.0, false);
    a.surface.mass_transfer_coefficient_m_s = 0.02;
    b.surface.mass_transfer_coefficient_m_s = 0.02;
    a.heat_controls.pure_liquid_enthalpy_range_j_kg = [f64::NAN; 2];
    b.heat_controls.maximum_iterations = 0;
    let phase = advance_chamber_with_gas(
        &[a.surface, b.surface],
        &initial,
        &source,
        0.2,
        phase_controls(),
        || false,
    )
    .unwrap();
    let zero = advance_chamber_with_gas_and_heat(
        &[b, a],
        &initial,
        &source,
        0.2,
        phase_controls(),
        GasLiquidChamberHeatSolve {
            maximum_iterations: 0,
        },
        || false,
    )
    .unwrap();
    assert_eq!(zero.phase(), &phase);
    assert_eq!(zero.surfaces(), phase.surfaces());
    assert_eq!(zero.step(), phase.step());
    assert_eq!(zero.heat().heat_into_liquid_j, 0.0);
    a.heat_transfer_coefficient_w_m2_k = 60.0;
    a.heat_controls = heat_controls(false);
    let old = crate::evaporation::gas::advance_surface_with_gas_and_heat(
        &Phase,
        &binding,
        a.surface.liquid,
        &initial,
        &source,
        temperature,
        0.02,
        60.0,
        0.2,
        phase_controls().surface,
        a.heat_controls,
        || false,
    )
    .unwrap();
    let one = advance_chamber_with_gas_and_heat(
        &[a],
        &initial,
        &source,
        0.2,
        phase_controls(),
        outer(),
        || false,
    )
    .unwrap();
    assert_eq!(one.gas(), old.gas());
    assert_eq!(one.surfaces()[0].liquid, old.liquid());
    assert_eq!(one.phase().gas_exchange(), old.phase().gas_exchange());
    assert_eq!(
        one.heat_surfaces()[0].phase_initial_liquid,
        old.phase_initial_liquid()
    );
    assert_eq!(one.heat().heat_into_liquid_j, old.heat().heat_into_liquid_j);
    assert_eq!(
        one.heat().gas_internal_energy_change_j,
        old.heat().gas_internal_energy_change_j
    );
    assert_eq!(
        one.step().liquid_enthalpy_change_j,
        old.step().liquid_enthalpy_change_j
    );
    assert_eq!(
        one.step().gas_internal_energy_change_j,
        old.step().gas_internal_energy_change_j
    );
}

#[test]
fn nonlinear_heat_precedes_one_joint_phase_selection_and_retains_exact_intermediate_states() {
    let source = catalog(-1e7);
    let binding = binding(-1e7);
    let initial = gas(&source, 330.0, 5000.0, 0.01);
    let mut inputs = [
        surface(4, &binding, liquid(300.0), 60.0, false),
        surface(8, &binding, liquid(340.0), 90.0, false),
    ];
    for input in &mut inputs {
        input.surface.mass_transfer_coefficient_m_s = 0.02;
    }
    let mut calls = Vec::new();
    let candidate = advance_chamber_with_gas_and_heat_mass_selection(
        &inputs,
        &initial,
        &source,
        0.2,
        phase_controls(),
        outer(),
        |id, mass| {
            calls.push((id, mass));
            Ok(mass)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls.len(), 2);
    assert_eq!(calls.iter().map(|row| row.0).collect::<Vec<_>>(), [4, 8]);
    for ((input, heat), phase) in inputs
        .iter()
        .zip(candidate.heat_surfaces())
        .zip(candidate.phase().surfaces())
    {
        assert_ne!(
            heat.phase_initial_liquid.specific_enthalpy_j_kg,
            input.surface.liquid.specific_enthalpy_j_kg
        );
        assert_eq!(phase.surface_temperature_k[0], heat.temperature_k[1]);
        assert_eq!(
            heat.liquid_enthalpy_change_j,
            input.surface.liquid.mass_kg
                * (heat.phase_initial_liquid.specific_enthalpy_j_kg
                    - input.surface.liquid.specific_enthalpy_j_kg)
        );
        let combined = candidate
            .surfaces()
            .iter()
            .find(|row| row.id == phase.id)
            .unwrap();
        assert_eq!(
            combined.liquid_enthalpy_change_j,
            phase.liquid_enthalpy_change_j + heat.liquid_enthalpy_change_j
        );
        assert_eq!(combined.energy_into_gas_j, phase.energy_into_gas_j);
        assert!(phase.transferred_mass_kg < 0.0);
    }
    assert_eq!(
        candidate.step().energy_into_gas_j,
        candidate.phase().step().energy_into_gas_j
    );
    near(
        candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
        0.0,
        phase_controls().surface.gas.energy_tolerance_j,
    );
    assert!(
        candidate.step().sum_absolute_energy_residual_j
            <= phase_controls().surface.gas.energy_tolerance_j
    );
}

#[test]
fn simultaneous_heat_has_first_order_physical_time_refinement() {
    let source = constant_source();
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 1000.0, 0.002);
    let start = constant_liquid(300.0, 0.01, 0.001);
    let capacity_g = constant_gas_capacity(&initial);
    let capacity_l = 2.0 * start.mass_kg * CP;
    let total_conductance = 2.0 * 60.0 * start.area_m2;
    let duration = 4.0;
    let equilibrium = (capacity_g * 320.0 + capacity_l * 300.0) / (capacity_g + capacity_l);
    let decay =
        fs_math::det::exp(-total_conductance * (1.0 / capacity_g + 1.0 / capacity_l) * duration);
    let exact_gas = equilibrium + capacity_l / (capacity_g + capacity_l) * 20.0 * decay;
    let run = |count: usize| {
        let mut current = initial.clone();
        let mut liquids = [start; 2];
        for _ in 0..count {
            let inputs = [
                surface(1, &binding, liquids[0], 60.0, true),
                surface(2, &binding, liquids[1], 60.0, true),
            ];
            let next = advance_chamber_with_gas_and_heat(
                &inputs,
                &current,
                &source,
                duration / count as f64,
                phase_controls(),
                outer(),
                || false,
            )
            .unwrap();
            let (rows, gas, _, _, _, _) = next.into_parts();
            liquids = [rows[0].liquid, rows[1].liquid];
            current = gas;
        }
        (current.state().state().temperature - exact_gas).abs()
    };
    let errors = [8, 16, 32].map(run);
    assert!(errors[0] > 1e-4, "inactive heat refinement {errors:?}");
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (0.8..1.2).contains(&order),
            "heat order {order}, errors {errors:?}"
        );
    }
}

#[test]
fn work_branch_and_late_cancellation_refuse_without_consuming_early_phase_selection() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 330.0, 1000.0, 0.01);
    let before = initial.clone();
    let inputs = [
        surface(1, &binding, liquid(280.0), 60.0, false),
        surface(2, &binding, liquid(360.0), 90.0, false),
    ];
    let mut selected = 0;
    let mut limited = inputs;
    for input in &mut limited {
        input.heat_controls.maximum_iterations = 1;
        input.heat_controls.temperature_tolerance_k = 1e-14;
    }
    assert!(
        advance_chamber_with_gas_and_heat_mass_selection(
            &limited,
            &initial,
            &source,
            0.2,
            phase_controls(),
            outer(),
            |_, mass| {
                selected += 1;
                Ok(mass)
            },
            || false,
        )
        .is_err()
    );
    assert_eq!(selected, 0);
    let mut bad_branch = inputs;
    bad_branch[0].heat_controls.pure_liquid_enthalpy_range_j_kg = [
        liquid(400.0).specific_enthalpy_j_kg,
        liquid(500.0).specific_enthalpy_j_kg,
    ];
    assert!(
        advance_chamber_with_gas_and_heat(
            &bad_branch,
            &initial,
            &source,
            0.2,
            phase_controls(),
            outer(),
            || false,
        )
        .is_err()
    );
    let mut calls = 0;
    advance_chamber_with_gas_and_heat(
        &inputs,
        &initial,
        &source,
        0.2,
        phase_controls(),
        outer(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let final_gate = calls;
    let mut calls = 0;
    assert_eq!(
        advance_chamber_with_gas_and_heat(
            &inputs,
            &initial,
            &source,
            0.2,
            phase_controls(),
            outer(),
            || {
                calls += 1;
                calls == final_gate
            },
        )
        .unwrap_err(),
        GasEvaporationError::Cancelled
    );
    assert_eq!(initial, before);
}
