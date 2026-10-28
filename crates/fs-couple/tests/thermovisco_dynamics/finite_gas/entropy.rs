//! Exact entropy storage and separated transfer/time errors on the original solve.
use super::*;
use fs_couple::finite_gas_thermal::FiniteGasEntropyReport;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

fn sensible_input(
    source: &Sources,
    temperature_k: f64,
    duration_s: f64,
) -> UniformEnthalpyStepInput<'_> {
    let old = Oscillator::try_new(
        &source.mechanical,
        &source.caloric,
        geometry(0.0),
        temperature_k,
        0.0,
        0.0,
    )
    .unwrap();
    UniformEnthalpyStepInput {
        curve: source.caloric.curve(),
        initial: *old.point().caloric_state(),
        mass_kg: old.point().mass_kg(),
        volume_m3: old.point().volume_m3(),
        surface_area_m2: old.geometry().thermal_surface_m2,
        internal_heat_j: 0.0,
        duration_s,
    }
}

fn environment<'a>(
    gas: &'a IsochoricGasParcel,
    body: &'a LumpedThermalEnvironment,
) -> FiniteGasThermalEnvironment<'a> {
    FiniteGasThermalEnvironment {
        gas,
        body,
        gas_heat_j: 0.0,
        gas_solve: gas_solve(),
        exchange: exchange(0.02),
        maximum_gas_temperature_increment_k: 100.0,
    }
}

fn balanced(entropy: FiniteGasEntropyReport) {
    assert!(entropy.balance_residual_j_k.abs() < 2e-18, "{entropy:?}");
    let accounted = entropy.convection_production_j_k
        + entropy.radiation_reservoir_input_j_k
        + entropy.radiation_transfer_production_j_k
        + entropy
            .unclassified_heat_at_endpoint_j_k
            .iter()
            .sum::<f64>()
        + entropy.temporal_defect_j_k.iter().sum::<f64>()
        + entropy.source_residual_j_k
        + entropy.algebraic_residual_j_k.iter().sum::<f64>();
    close(
        entropy.storage_change_j_k.iter().sum::<f64>(),
        accounted,
        2e-18,
    );
}

#[test]
fn closed_hot_and_cold_pairs_use_cv_entropy_and_separate_time_error() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let body = boundary(&source, 100.0);
        for (body_k, gas_k) in [(320.0, 300.0), (300.0, 320.0)] {
            let reservoir = gas(gas_k, 1e-6);
            let input = sensible_input(&source, body_k, 0.01);
            let step = environment(&reservoir, &body).advance(cx, input).unwrap();
            let entropy = step.report.entropy;
            let tb = step.state.temperature_k();
            let tg = step.report.next_gas.state().state().temperature;
            let cb = input.mass_kg * 2.0;
            let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
            let cg = reservoir.mass_kg() * cv;
            close(
                entropy.storage_change_j_k[0],
                cb * ((tb - body_k) / body_k).ln_1p(),
                1e-15,
            );
            close(
                entropy.storage_change_j_k[1],
                cg * ((tg - gas_k) / gas_k).ln_1p(),
                1e-15,
            );
            assert!(entropy.storage_change_j_k.iter().sum::<f64>() > 0.0);
            assert!(entropy.convection_production_j_k > 0.0);
            assert!(entropy.temporal_defect_j_k.iter().all(|value| *value > 0.0));
            assert_eq!(entropy.unclassified_heat_at_endpoint_j_k, [0.0; 2]);
            assert_eq!(entropy.radiation_reservoir_input_j_k, 0.0);
            assert_eq!(entropy.radiation_transfer_production_j_k, 0.0);
            close(
                step.report.solid_to_gas_heat_j * (tb - tg) / (tb * tg),
                entropy.convection_production_j_k,
                1e-18,
            );
            balanced(entropy);
        }
    });
}

#[test]
fn latent_enthalpy_supplies_entropy_once_without_a_temperature_jump() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            source.caloric.curve().material_card_identity(),
            [
                (-60.0, 280.0, 0.0),
                (0.0, 310.0, 0.0),
                (100.0, 310.0, 1.0),
                (160.0, 340.0, 1.0),
            ]
            .map(|(h, t, liquid)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: liquid,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
        )
        .unwrap();
        let input = UniformEnthalpyStepInput {
            curve: &curve,
            initial: curve.state_at_specific_enthalpy(25.0).unwrap(),
            mass_kg: 0.001,
            volume_m3: 1e-6,
            surface_area_m2: 0.0013,
            internal_heat_j: 0.0,
            duration_s: 0.01,
        };
        let reservoir = gas(330.0, 1e-6);
        let body = boundary(&source, 100.0);
        let step = environment(&reservoir, &body).advance(cx, input).unwrap();
        let entropy = step.report.entropy;
        assert_eq!(step.state.temperature_k(), 310.0);
        assert!(step.state.liquid_mass_fraction() > input.initial.liquid_mass_fraction());
        assert!(step.report.solid_to_gas_heat_j < 0.0);
        let stored = input.mass_kg
            * (step.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
        close(entropy.storage_change_j_k[0], stored / 310.0, 1e-18);
        close(
            entropy.storage_change_j_k[0],
            -step.report.solid_to_gas_heat_j / 310.0,
            1e-14,
        );
        close(entropy.temporal_defect_j_k[0], 0.0, 1e-18);
        assert!(entropy.convection_production_j_k > 0.0);
        balanced(entropy);
    });
}

#[test]
fn radiation_uses_actual_enclosure_and_prescribed_heat_remains_unclassified() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        for enclosure_k in [250.0, 400.0] {
            let reservoir = gas(305.0, 1e-6);
            let body = LumpedThermalEnvironment {
                temperature_k: 12345.0, // Replaced by finite gas, never an entropy bath.
                radiation_temperature_k: enclosure_k,
                transport: LumpedThermalTransport::try_declared(10.0, 0.4).unwrap(),
                ..boundary(&source, 0.0)
            };
            let mut environment = environment(&reservoir, &body);
            environment.gas_heat_j = 3e-4;
            let input = UniformEnthalpyStepInput {
                internal_heat_j: 2e-4,
                ..sensible_input(&source, 300.0, 0.001)
            };
            let step = environment.advance(cx, input).unwrap();
            let entropy = step.report.entropy;
            let tb = step.state.temperature_k();
            let tg = step.report.next_gas.state().state().temperature;
            let radiation = step.report.radiation_heat_j;
            assert_eq!(radiation.is_sign_positive(), enclosure_k > tb);
            assert_eq!(step.report.solid_to_gas_heat_j, 0.0);
            assert_eq!(entropy.convection_production_j_k, 0.0);
            assert_eq!(step.report.gas.supplied_heat_j, environment.gas_heat_j);
            close(
                step.report.gas.stored_energy_change_j,
                environment.gas_heat_j,
                2e-12,
            );
            close(
                entropy.radiation_reservoir_input_j_k,
                radiation / enclosure_k,
                1e-20,
            );
            close(
                entropy.radiation_transfer_production_j_k,
                radiation * (1.0 / tb - 1.0 / enclosure_k),
                1e-18,
            );
            assert!(entropy.radiation_transfer_production_j_k > 0.0);
            close(
                entropy.unclassified_heat_at_endpoint_j_k[0],
                input.internal_heat_j / tb,
                1e-20,
            );
            close(
                entropy.unclassified_heat_at_endpoint_j_k[1],
                environment.gas_heat_j / tg,
                1e-20,
            );
            balanced(entropy);
        }
    });
}

#[test]
fn huge_nasa_entropy_datum_does_not_erase_finite_inventory_entropy_change() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let body = boundary(&source, 100.0);
        let input = sensible_input(&source, 320.0, 0.01);
        let mut reference = None;
        for datum in [0.0, 1e20] {
            let reservoir = gas_with_entropy_reference(300.0, 1e-6, 0.025, datum);
            let step = environment(&reservoir, &body).advance(cx, input).unwrap();
            let entropy = step.report.entropy;
            assert!(entropy.storage_change_j_k[1] > 1e-6);
            balanced(entropy);
            if let Some((expected, heat, trials)) = reference {
                assert_eq!(entropy, expected);
                assert_eq!(step.report.solid_to_gas_heat_j, heat);
                assert_eq!(step.report.trials, trials);
                // Subtracting absolute entropy would lose this finite change.
                assert_eq!(
                    reservoir.state().specific_entropy_j_kg_k().unwrap(),
                    step.report
                        .next_gas
                        .state()
                        .specific_entropy_j_kg_k()
                        .unwrap()
                );
            } else {
                reference = Some((entropy, step.report.solid_to_gas_heat_j, step.report.trials));
            }
        }
    });
}

#[test]
fn fixed_duration_refinement_reduces_be_entropy_defect_and_temperature_error() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let body = boundary(&source, 100.0);
        let input = sensible_input(&source, 320.0, 0.01);
        let run = |steps: usize| {
            let mut reservoir = gas(300.0, 1e-6);
            let mut state = input.initial;
            let cb = input.mass_kg * 2.0;
            let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
            let cg = reservoir.mass_kg() * cv;
            let mut temporal = 0.0;
            let mut storage = 0.0;
            for _ in 0..steps {
                let step = environment(&reservoir, &body)
                    .advance(
                        cx,
                        UniformEnthalpyStepInput {
                            initial: state,
                            duration_s: input.duration_s / steps as f64,
                            ..input
                        },
                    )
                    .unwrap();
                balanced(step.report.entropy);
                temporal += step.report.entropy.temporal_defect_j_k.iter().sum::<f64>();
                storage += step.report.entropy.storage_change_j_k.iter().sum::<f64>();
                state = step.state;
                reservoir = step.report.next_gas;
            }
            let tb = state.temperature_k();
            let tg = reservoir.state().state().temperature;
            close(
                storage,
                cb * ((tb - 320.0) / 320.0).ln_1p() + cg * ((tg - 300.0) / 300.0).ln_1p(),
                2e-15,
            );
            let equilibrium = (cb * 320.0 + cg * 300.0) / (cb + cg);
            let decay = (-body.convection_w_per_m2_k
                * input.surface_area_m2
                * (1.0 / cb + 1.0 / cg)
                * input.duration_s)
                .exp();
            let exact_body = equilibrium + 20.0 * cg / (cb + cg) * decay;
            let exact_gas = equilibrium - 20.0 * cb / (cb + cg) * decay;
            (temporal, (tb - exact_body).abs() + (tg - exact_gas).abs())
        };
        let values = [8, 16, 32].map(run);
        for pair in values.windows(2) {
            assert!(pair[1].0 > 0.0);
            let ratio = pair[0].0 / pair[1].0;
            assert!((1.7..2.2).contains(&ratio), "{values:?}");
            assert!(pair[1].1 < pair[0].1, "{values:?}");
        }
    });
}
