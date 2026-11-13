//! Native coarse/fine cooling, latent inventory and whole-window source accounts.
use super::*;
use fs_couple::finite_gas_thermal::{ConvectionCoefficient,
    window::{FiniteGasWindowError, FiniteGasWindowSolve}};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

fn controls(duration: f64, energy_error: f64) -> FiniteGasWindowSolve {
    FiniteGasWindowSolve {
        maximum_step_s: duration, minimum_step_s: 1e-8, maximum_attempts: 8192,
        energy_absolute_tolerance_j: energy_error, temperature_absolute_tolerance_k: 1000.0,
        phase_fraction_absolute_tolerance: 1.0, relative_tolerance: 0.0,
        energy_balance_tolerance_j: 1e-7,
    }
}
fn input(source: &Sources, temperature: f64, duration_s: f64) -> UniformEnthalpyStepInput<'_> {
    let old = Oscillator::try_new(&source.mechanical, &source.caloric,
        geometry(0.0), temperature, 0.0, 0.0).unwrap();
    UniformEnthalpyStepInput {
        curve: source.caloric.curve(), initial: *old.point().caloric_state(),
        mass_kg: old.point().mass_kg(), volume_m3: old.point().volume_m3(),
        surface_area_m2: old.geometry().thermal_surface_m2,
        internal_heat_j: 0.0, duration_s,
    }
}

#[test]
fn sensible_window_refines_to_the_original_two_capacity_decay() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let reservoir = gas(300.0, 1e-6);
        let saved = reservoir.clone();
        let body = boundary(&source, 100.0);
        let input = input(&source, 320.0, 0.02);
        let environment = FiniteGasThermalEnvironment {
            gas: &reservoir, body: &body, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(1.0),
            maximum_gas_temperature_increment_k: 1.0,
        };
        assert!(environment.advance(cx, input).is_err());
        let cb = input.mass_kg * 2.0;
        let cg = 2.5 * 100000.0 * reservoir.volume_m3() / 300.0;
        let equilibrium = (cb * 320.0 + cg * 300.0) / (cb + cg);
        let decay = (-body.convection_w_per_m2_k * input.surface_area_m2
            * (1.0 / cb + 1.0 / cg) * input.duration_s).exp();
        let exact = [equilibrium + 20.0 * cg / (cb + cg) * decay,
            equilibrium - 20.0 * cb / (cb + cg) * decay];
        let mut errors = Vec::new();
        let mut defects = Vec::new();
        for tolerance in [2e-5, 5e-6, 1.25e-6] {
            let c = controls(input.duration_s, tolerance);
            let step = environment.advance_adaptive(cx, input, c).unwrap();
            let t = [step.state.temperature_k(), step.report.next_gas.state().state().temperature];
            errors.push((t[0] - exact[0]).abs() + (t[1] - exact[1]).abs());
            let entropy = step.report.entropy;
            defects.push(entropy.temporal_defect_j_k.iter().sum::<f64>());
            close(entropy.storage_change_j_k[0], cb * (t[0] / 320.0).ln(), 1e-13);
            close(entropy.storage_change_j_k[1], cg * (t[1] / 300.0).ln(), 1e-13);
            assert!(entropy.convection_production_j_k > 0.0);
            assert!(step.report.combined_energy_residual_j.abs() <= c.energy_balance_tolerance_j);
            assert!(step.report.sum_absolute_energy_defect_j <= c.energy_balance_tolerance_j);
            assert!(step.report.maximum_accepted_normalized_error <= 1.0);
            assert!(step.report.attempted_intervals <= c.maximum_attempts);
            assert!(step.report.rejected_comparisons > 0);
            assert_eq!(step.report.samples.last().unwrap().time_s, input.duration_s);
            let mut last_t = 300.0;
            for sample in &step.report.samples {
                assert!((sample.gas_temperature_k - last_t).abs() <= 1.0);
                last_t = sample.gas_temperature_k;
            }
        }
        assert!(errors.windows(2).all(|pair| pair[1] < pair[0]), "{errors:?}");
        assert!(defects.windows(2).all(|pair| pair[1] < pair[0]), "{defects:?}");
        assert_eq!(reservoir, saved);
    });
}

#[test]
fn retries_and_comparisons_spend_each_prescribed_joule_once() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let reservoir = gas(300.0, 1e-6);
        let body = boundary(&source, 0.0);
        let input = UniformEnthalpyStepInput { internal_heat_j: 0.004,
            ..input(&source, 300.0, 0.019) };
        let environment = FiniteGasThermalEnvironment {
            gas: &reservoir, body: &body, gas_heat_j: 0.006,
            gas_solve: gas_solve(), exchange: exchange(1.0),
            maximum_gas_temperature_increment_k: 1.0,
        };
        let step = environment.advance_adaptive(cx, input, controls(input.duration_s, 1e-7)).unwrap();
        close(step.state.temperature_k(), 302.0, 1e-8);
        close(step.report.next_gas.state().state().temperature, 307.2, 1e-8);
        close(step.report.internal_body_heat_j, input.internal_heat_j, 1e-15);
        close(step.report.external_gas_heat_j, environment.gas_heat_j, 1e-15);
        close(step.report.body_enthalpy_change_j, 0.004, 1e-10);
        close(step.report.gas_internal_energy_change_j, 0.006, 1e-10);
        assert_eq!(step.report.solid_to_gas_heat_j, 0.0);
        assert_eq!(step.external_heat_j, 0.0);
        assert!(step.report.rejected_comparisons > 0);
        assert!(step.report.attempted_intervals > step.report.samples.len());
        let radiating = LumpedThermalEnvironment {
            radiation_temperature_k: 400.0,
            transport: LumpedThermalTransport::try_declared(10.0, 0.4).unwrap(),
            ..body.clone()
        };
        let heated = FiniteGasThermalEnvironment { body: &radiating, ..environment }
            .advance_adaptive(cx, input, controls(input.duration_s, 1e-7)).unwrap();
        assert!(heated.report.radiation_heat_j > 0.0);
        close(heated.report.body_enthalpy_change_j,
            input.internal_heat_j + heated.report.radiation_heat_j, 1e-8);
        close(heated.report.gas_internal_energy_change_j, 0.006, 1e-10);
        close(heated.report.entropy.radiation_reservoir_input_j_k,
            heated.report.radiation_heat_j / 400.0, 1e-14);
        assert_eq!(heated.external_heat_j, heated.report.radiation_heat_j);
    });
}

#[test]
fn phase_error_controls_latent_storage_while_body_temperature_is_flat() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            source.caloric.curve().material_card_identity(),
            [(-40.0, 280.0, 0.0), (0.0, 300.0, 0.0), (100.0, 300.0, 1.0), (180.0, 340.0, 1.0)]
                .map(|(h, t, liquid)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
                    temperature_k: t, liquid_mass_fraction: liquid, bulk_density_kg_m3: 1000.0 }).to_vec(),
        ).unwrap();
        let input = UniformEnthalpyStepInput { curve: &curve,
            initial: curve.state_at_specific_enthalpy(50.0).unwrap(), mass_kg: 0.001,
            volume_m3: 1e-6, surface_area_m2: 0.0013, internal_heat_j: 0.0, duration_s: 0.01 };
        let reservoir = gas(330.0, 1e-6);
        let body = boundary(&source, 100.0);
        let environment = FiniteGasThermalEnvironment { gas: &reservoir, body: &body, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(1.0), maximum_gas_temperature_increment_k: 100.0 };
        let cg = 2.5 * 100000.0 * reservoir.volume_m3() / 330.0;
        let exact = 0.5 + cg * 30.0 * (1.0 - (-body.convection_w_per_m2_k
            * input.surface_area_m2 / cg * input.duration_s).exp()) / (input.mass_kg * 100.0);
        let run = |fraction| environment.advance_adaptive(cx, input, FiniteGasWindowSolve {
            phase_fraction_absolute_tolerance: fraction, ..controls(input.duration_s, 1000.0)
        }).unwrap();
        let loose = run(1.0);
        let tight = run(1e-5);
        assert!(tight.report.samples.len() > loose.report.samples.len());
        assert!(tight.report.samples.iter().all(|s| s.body_state.temperature_k() == 300.0));
        assert!((tight.state.liquid_mass_fraction() - exact).abs()
            < (loose.state.liquid_mass_fraction() - exact).abs());
        close(tight.report.entropy.storage_change_j_k[0], tight.report.body_enthalpy_change_j / 300.0, 1e-13);
    });
}

#[test]
fn budget_source_failure_and_cancellation_publish_no_partial_window() {
    scope(|cx, gate| {
        let source = sources(2.0, 1000.0, false, false);
        let reservoir = gas(300.0, 1e-6);
        let saved = reservoir.clone();
        let body = boundary(&source, 100.0);
        let input = input(&source, 320.0, 0.02);
        let environment = FiniteGasThermalEnvironment { gas: &reservoir, body: &body, gas_heat_j: 0.0,
            gas_solve: gas_solve(), exchange: exchange(1.0), maximum_gas_temperature_increment_k: 100.0 };
        assert!(matches!(environment.advance_adaptive(cx, input, FiniteGasWindowSolve {
            maximum_attempts: 1, ..controls(input.duration_s, 1e-6)
        }), Err(FiniteGasWindowError::WorkLimit { .. })));
        let mut calls = 0;
        let refused = environment.advance_adaptive_with_convection(cx, input, controls(input.duration_s, 1e-6), |_| {
            calls += 1;
            Err::<ConvectionCoefficient<()>, _>("unsupported convection source")
        });
        assert!(matches!(refused, Err(FiniteGasWindowError::Step(_))));
        assert_eq!(calls, 1);
        let cancelled = environment.advance_adaptive_with_convection(cx, input, controls(input.duration_s, 1e-6), |_| {
            gate.request();
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient { coefficient_w_per_m2_k: 100.0, source: () })
        });
        assert!(matches!(cancelled, Err(FiniteGasWindowError::Cancelled)));
        assert_eq!(reservoir, saved);
        assert_eq!(input.initial.temperature_k(), 320.0);
    });
}

#[test]
fn nonlinear_actual_wall_window_refines_to_independent_two_capacity_decay() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let reservoir = gas(300.0, 1e-6);
        let saved = reservoir.clone();
        let body = boundary(&source, 0.0);
        let input = input(&source, 320.0, 0.02);
        let environment = FiniteGasThermalEnvironment {
            gas: &reservoir,
            body: &body,
            gas_heat_j: 0.0,
            gas_solve: gas_solve(),
            exchange: exchange(1.0),
            maximum_gas_temperature_increment_k: 20.0,
        };
        let wall = fs_couple::finite_gas_thermal::WallConvectionSolve {
            temperature_increment_bracket_k: [-30.0, 30.0],
            temperature_tolerance_k: 1e-9,
            maximum_iterations: 128,
        };
        let cb = input.mass_kg * 2.0;
        let cg = 2.5 * 100000.0 * reservoir.volume_m3() / 300.0;
        let equilibrium = (cb * 320.0 + cg * 300.0) / (cb + cg);
        let h0 = 80.0;
        let slope = 0.02;
        // h = h0[1 + slope(Twall - Teq)] gives
        // dD/dt = -K D(1 + a D), D = Twall - Tgas.
        // This exact solution uses neither the native step nor its time controller.
        let a = slope * cg / (cb + cg);
        let k = input.surface_area_m2 * h0 * (1.0 / cb + 1.0 / cg);
        let decay = (-k * input.duration_s).exp();
        let difference = 20.0 * decay / (1.0 + a * 20.0 * (1.0 - decay));
        let exact = [
            equilibrium + cg / (cb + cg) * difference,
            equilibrium - cb / (cb + cg) * difference,
        ];
        let mut errors = Vec::new();
        let mut refinements = Vec::new();
        for tolerance in [2e-5, 5e-6, 1.25e-6] {
            let control = controls(input.duration_s, tolerance);
            let next = environment
                .advance_adaptive_with_wall_convection(
                    cx,
                    input,
                    control,
                    wall,
                    |_, temperature| {
                        Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                            coefficient_w_per_m2_k: h0
                                * (1.0 + slope * (temperature - equilibrium)),
                            source: temperature,
                        })
                    },
                )
                .unwrap();
            let r = &next.report;
            errors.push(
                (next.state.temperature_k() - exact[0]).abs()
                    + (r.next_gas.state().state().temperature - exact[1]).abs(),
            );
            refinements.push(r.samples.len());
            assert!(r.rejected_comparisons > 0);
            assert!(r.maximum_accepted_normalized_error <= 1.0);
            assert!(r.sum_absolute_energy_defect_j <= control.energy_balance_tolerance_j);
            close(
                r.body_enthalpy_change_j + r.gas_internal_energy_change_j,
                0.0,
                control.energy_balance_tolerance_j,
            );
            close(r.entropy.balance_residual_j_k, 0.0, 1e-11);
            assert_eq!(r.samples.last().unwrap().time_s, input.duration_s);
            let mut previous = 0.0;
            for sample in &r.samples {
                let dt = sample.time_s - previous;
                assert!(dt > 0.0);
                let feedback = &sample.convection.source;
                assert_eq!(feedback.source, feedback.temperature_k);
                close(
                    feedback.endpoint_source,
                    sample.body_state.temperature_k(),
                    1e-12,
                );
                assert!(feedback.temperature_residual_k.abs() <= wall.temperature_tolerance_k);
                assert!(
                    feedback.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j
                );
                close(
                    sample.convection.coefficient_w_per_m2_k,
                    h0 * (1.0 + slope * (feedback.temperature_k - equilibrium)),
                    1e-12,
                );
                close(
                    sample.solid_to_gas_heat_j,
                    dt * input.surface_area_m2
                        * sample.convection.coefficient_w_per_m2_k
                        * (sample.body_state.temperature_k() - sample.gas_temperature_k),
                    environment.exchange.source_tolerance_j,
                );
                previous = sample.time_s;
            }
            assert!(r.samples.last().unwrap().convection.source.temperature_k < 315.0);
            assert_eq!(reservoir, saved);
        }
        assert!(
            errors.windows(2).all(|pair| pair[1] < pair[0]),
            "{errors:?}"
        );
        assert!(
            refinements.windows(2).all(|pair| pair[1] > pair[0]),
            "{refinements:?}"
        );
    });
}
