use super::super::WaterCellControls;
use super::*;
use fs_blake3::ContentHash;
use fs_material::phase::water::isochoric::IsochoricWaterControls;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::thermal::conduction::ConductivityCurve;
use fs_qty::Temperature;

fn water(temperature_k: f64) -> EquilibriumWaterCell {
    EquilibriumWaterCell::at_temperature(
        0.001,
        0.02,
        Temperature::new(temperature_k),
        0.0,
        WaterCellControls {
            equilibrium: IsochoricWaterControls {
                minimum_temperature_k: temperature_k - 5.0,
                maximum_temperature_k: temperature_k + 5.0,
                minimum_pressure_pa: 100.0,
                maximum_pressure_pa: 100e6,
                volume_tolerance_m3_kg: 1e-14,
                internal_energy_tolerance_j_kg: 1e-7,
                maximum_pressure_iterations: 96,
                maximum_temperature_iterations: 96,
            },
            maximum_relative_mass_exchange: 0.1,
            maximum_inflow_parcels: 2,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-7,
            volume_tolerance_m3: 1e-13,
        },
        || false,
    )
    .unwrap()
}

fn wall_curve() -> EquilibriumEnthalpyPhaseCurve {
    // Explicit analytical fixture: cp=500 J/(kg K), fixed density, all solid.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xc7; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 280.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2700.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 110_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2700.0,
            },
            // The original chart requires its complete liquid endpoint.
            // All tested states remain below this authored melting plateau.
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 210_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 2700.0,
            },
        ],
    )
    .unwrap()
}
fn wall_at(curve: &EquilibriumEnthalpyPhaseCurve, temperature_k: f64) -> EquilibriumPhaseState {
    curve
        .state_at_specific_enthalpy(500.0 * (temperature_k - 280.0))
        .unwrap()
}
fn body<'a>(
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    initial: EquilibriumPhaseState,
    duration_s: f64,
    internal_heat_j: f64,
) -> UniformEnthalpyStepInput<'a> {
    UniformEnthalpyStepInput {
        curve,
        initial,
        mass_kg: 0.4,
        volume_m3: 0.4 / 2700.0,
        surface_area_m2: 0.01,
        internal_heat_j,
        duration_s,
    }
}
fn controls(bracket: [f64; 2]) -> WaterWallControls {
    WaterWallControls {
        exchange: HeatExchangeSolve {
            heat_bracket_j: bracket,
            source_tolerance_j: 1e-8,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 96,
        },
        conduction: LayerConductionConfig {
            max_iterations: 96,
            absolute_heat_tolerance_w: 1e-11,
            relative_heat_tolerance: 1e-12,
        },
        maximum_step_s: 2.0,
        maximum_temperature_increment_k: [5.0, 5.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance={tolerance}"
    );
}

#[test]
fn constant_capacity_wall_and_conductance_recover_an_independent_endpoint() {
    let mut water = water(350.0);
    let initial_water = water.state();
    // Specify the target temperature independently, then obtain ONLY the
    // original IF97 caloric increment. No inverse or coupled solve is used.
    let target_water = initial_water
        .source_state()
        .law()
        .state_at_temperature(Temperature::new(350.1), || false)
        .unwrap();
    let exact_heat = initial_water.mass_kg()
        * (target_water.specific_internal_energy_j_kg()
            - initial_water.source_state().specific_internal_energy_j_kg());
    let curve = wall_curve();
    let mut wall = wall_at(&curve, 370.0);
    let initial_wall = wall;
    let wall_capacity = 0.4 * 500.0;
    let exact_wall_temperature = 370.0 - exact_heat / wall_capacity;
    // G=A*k/L=10 W/K. This duration makes those independently selected
    // endpoints the exact backward-Euler solution for the finite wall.
    let duration = exact_heat / (10.0 * (exact_wall_temperature - 350.1));
    let conductivity = [ConductivityCurve::new(vec![[280.0, 2.0], [500.0, 2.0]]).unwrap()];
    let lengths = [0.001];
    let layers = LayerStackConduction::new(&conductivity, &lengths, &[], 0.005).unwrap();
    let step = advance_water_wall(
        &mut water,
        &mut wall,
        body(&curve, initial_wall, duration, 0.0),
        layers,
        0.0,
        duration,
        controls([0.0, 2.0 * exact_heat]),
        || false,
    )
    .unwrap();
    near(step.wall_to_water_heat_j, exact_heat, 2e-7);
    near(wall.temperature_k(), exact_wall_temperature, 2e-9);
    near(water.state().source_state().temperature_k(), 350.1, 2e-9);
    // The finite-capacity closed form differs from an infinite wall bath.
    let finite_wall_heat =
        duration * 10.0 * (370.0 - 350.1) / (1.0 + duration * 10.0 / wall_capacity);
    near(step.wall_to_water_heat_j, finite_wall_heat, 2e-7);
    assert!(step.wall_to_water_heat_j < duration * 10.0 * (370.0 - 350.1));
    assert_eq!(step.wall.external_heat_j, -step.wall_to_water_heat_j);
    assert_eq!(water.state().mass_kg(), initial_water.mass_kg());
    assert_eq!(water.state().volume_m3(), initial_water.volume_m3());
    near(step.combined_energy_residual_j, 0.0, 1e-7);
    near(step.source_residual_j, 0.0, 1e-8);
    near(
        step.entropy.wall_change_j_k,
        wall_capacity * (exact_wall_temperature / 370.0).ln(),
        1e-10,
    );
    assert!(step.entropy.conduction_generation_j_k > 0.0);
}

#[test]
fn nonlinear_native_layers_drive_evaporation_and_condensation_with_finite_storage() {
    let curve = wall_curve();
    let conductivity = [
        ConductivityCurve::new(vec![[280.0, 1.0], [500.0, 3.2]]).unwrap(),
        ConductivityCurve::new(vec![[280.0, 4.0], [500.0, 4.0]]).unwrap(),
    ];
    let lengths = [0.002, 0.001];
    let contacts = [0.001];
    let layers = LayerStackConduction::new(&conductivity, &lengths, &contacts, 0.01).unwrap();
    for wall_temperature in [420.0, 380.0] {
        let mut water = water(400.0);
        let initial = water.state();
        let mut wall = wall_at(&curve, wall_temperature);
        let initial_wall = wall;
        let step = advance_water_wall(
            &mut water,
            &mut wall,
            body(&curve, initial_wall, 1.0, 0.0),
            layers,
            0.0,
            1.0,
            controls([-300.0, 300.0]),
            || false,
        )
        .unwrap();
        let heating = wall_temperature > 400.0;
        assert_eq!(step.wall_to_water_heat_j > 0.0, heating);
        assert_eq!(
            water.state().phase_masses_kg()[1] > initial.phase_masses_kg()[1],
            heating
        );
        assert_eq!(
            water.state().source_state().pressure_pa() > initial.source_state().pressure_pa(),
            heating
        );
        assert_eq!(wall.temperature_k() < initial_wall.temperature_k(), heating);
        assert!(water.state().phase_masses_kg().into_iter().all(|m| m > 0.0));
        // Independent exact Kirchhoff integrals for the declared affine k(T),
        // constant second layer and intervening contact resistance.
        let response = &step.wall.report;
        let faces = response.layer_temperatures_k();
        let (ta, tb) = (faces[0][0], faces[0][1]);
        let first = 0.01 / 0.002 * (ta - tb) * (1.0 + 0.01 * (0.5 * (ta + tb) - 280.0));
        let second = 0.01 / 0.001 * 4.0 * (faces[1][0] - faces[1][1]);
        near(first, response.power_a_to_b_w, 1e-8);
        near(second, response.power_a_to_b_w, 1e-8);
        near(
            faces[0][1] - faces[1][0],
            contacts[0] * response.power_a_to_b_w / 0.01,
            1e-8,
        );
        near(step.wall_to_water_heat_j, response.power_a_to_b_w, 1e-8);
        near(step.combined_energy_residual_j, 0.0, 1e-7);
        assert!(step.entropy.conduction_generation_j_k > 0.0);
        near(step.entropy.interface_source_residual_j_k, 0.0, 1e-10);
    }
}

#[test]
fn repeated_native_history_keeps_conserved_energy_and_exposes_source_defects() {
    let curve = wall_curve();
    let conductivity = [ConductivityCurve::new(vec![[280.0, 2.0], [500.0, 2.0]]).unwrap()];
    let lengths = [0.001];
    let layers = LayerStackConduction::new(&conductivity, &lengths, &[], 0.005).unwrap();
    let run = || {
        let mut water = water(400.0);
        let mut wall = wall_at(&curve, 420.0);
        let first_water = water.state();
        let first_wall = wall;
        let mut reports = Vec::new();
        for end in [0.25, 0.5, 0.75, 1.0] {
            let initial_wall = wall;
            let input = body(&curve, initial_wall, end - water.state().time_s(), 2.0);
            let report = advance_water_wall(
                &mut water,
                &mut wall,
                input,
                layers,
                3.0,
                end,
                controls([-100.0, 100.0]),
                || false,
            )
            .unwrap();
            assert_eq!(
                report.water.state.internal_energy_j(),
                report.water.initial.internal_energy_j() + (3.0 + report.wall_to_water_heat_j)
            );
            let source_defect_change = report.water.state.source_energy_residual_j()
                - report.water.initial.source_energy_residual_j();
            near(
                report.source_combined_energy_residual_j - report.combined_energy_residual_j,
                source_defect_change,
                2e-11,
            );
            near(
                report.entropy.unaccounted_change_j_k,
                report.entropy.wall_change_j_k + report.entropy.water_source_change_j_k
                    - report.entropy.transferred_interface_entropy_j_k,
                1e-15,
            );
            reports.push(report);
        }
        near(
            (water.state().internal_energy_j() - first_water.internal_energy_j())
                + 0.4 * (wall.specific_enthalpy_j_kg() - first_wall.specific_enthalpy_j_kg()),
            4.0 * (3.0 + 2.0),
            4e-7,
        );
        (water, wall, reports)
    };
    assert_eq!(run(), run());
}

#[test]
fn late_cancellation_missing_root_and_temperature_limits_preserve_both_states() {
    let curve = wall_curve();
    let conductivity = [ConductivityCurve::new(vec![[280.0, 2.0], [500.0, 2.0]]).unwrap()];
    let lengths = [0.001];
    let layers = LayerStackConduction::new(&conductivity, &lengths, &[], 0.005).unwrap();
    let initial_water = water(400.0);
    let initial_wall = wall_at(&curve, 420.0);
    let input = body(&curve, initial_wall, 0.25, 0.0);
    let mut accepted_water = initial_water;
    let mut accepted_wall = initial_wall;
    let mut total_polls = 0usize;
    let accepted = advance_water_wall(
        &mut accepted_water,
        &mut accepted_wall,
        input,
        layers,
        0.0,
        0.25,
        controls([-100.0, 100.0]),
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    assert!(total_polls > 10);
    let mut rejected_water = initial_water;
    let mut rejected_wall = initial_wall;
    let mut polls = 0usize;
    let result = advance_water_wall(
        &mut rejected_water,
        &mut rejected_wall,
        input,
        layers,
        0.0,
        0.25,
        controls([-100.0, 100.0]),
        || {
            polls += 1;
            polls == total_polls
        },
    );
    assert_eq!(result, Err(WaterWallError::Cancelled));
    assert_eq!(rejected_water, initial_water);
    assert_eq!(rejected_wall, initial_wall);
    assert!(matches!(
        advance_water_wall(
            &mut rejected_water,
            &mut rejected_wall,
            input,
            layers,
            0.0,
            0.25,
            controls([-0.1, 0.1]),
            || false,
        ),
        Err(WaterWallError::OutsideBracket { .. })
    ));
    let mut tight = controls([-100.0, 100.0]);
    tight.maximum_temperature_increment_k[0] = 0.01;
    assert!(matches!(
        advance_water_wall(
            &mut rejected_water,
            &mut rejected_wall,
            input,
            layers,
            0.0,
            0.25,
            tight,
            || false,
        ),
        Err(WaterWallError::TemperatureIncrement { body: "wall", .. })
    ));
    assert_eq!(rejected_water, initial_water);
    assert_eq!(rejected_wall, initial_wall);
    let retry = advance_water_wall(
        &mut rejected_water,
        &mut rejected_wall,
        input,
        layers,
        0.0,
        0.25,
        controls([-100.0, 100.0]),
        || false,
    )
    .unwrap();
    assert_eq!(retry, accepted);
    assert_eq!(rejected_water, accepted_water);
    assert_eq!(rejected_wall, accepted_wall);
}
