//! Shared-gas windows retain individual component error and native ownership.
use super::*;
use fs_couple::finite_gas_thermal::bodies::{
    FiniteGasBodiesWindowError, FiniteGasBodiesWindowSolve,
};
use fs_couple::finite_gas_thermal::window::FiniteGasWindowSolve;
use std::cell::RefCell;

fn controls(duration: f64, energy_error: f64) -> FiniteGasBodiesWindowSolve {
    FiniteGasBodiesWindowSolve {
        time: FiniteGasWindowSolve {
            maximum_step_s: duration,
            minimum_step_s: 1e-8,
            maximum_attempts: 8192,
            energy_absolute_tolerance_j: energy_error,
            temperature_absolute_tolerance_k: 1000.0,
            phase_fraction_absolute_tolerance: 1.0,
            relative_tolerance: 0.0,
            energy_balance_tolerance_j: 1e-7,
        },
        maximum_retained_body_steps: 32_768,
    }
}

#[test]
fn g1_shared_window_refines_both_independent_decay_modes_and_is_id_order_invariant() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let wall = boundary(&source, 100.0);
        let gas = gas(300.0, 1e-6);
        let saved = gas.clone();
        let duration = 0.02;
        let mut inputs = [input(&source, 320.0, 0.0), input(&source, 290.0, 0.0)];
        for input in &mut inputs {
            input.duration_s = duration;
        }
        let bodies = [
            FiniteGasBody {
                id: 20,
                environment: &wall,
                input: inputs[0],
            },
            FiniteGasBody {
                id: 10,
                environment: &wall,
                input: inputs[1],
            },
        ];
        let environment = FiniteGasBodiesEnvironment {
            maximum_gas_temperature_increment_k: 1.0,
            exchange: exchange(1.0),
            ..shared(&gas, 0.0)
        };
        // Two identical capacities admit independent mean/gas and antisymmetric
        // eigenmodes. This reference uses no native thermal or feedback solve.
        let cb = 2.0 * inputs[0].mass_kg;
        let cg = 2.5 * 100000.0 * gas.volume_m3() / 300.0;
        let conductance = wall.convection_w_per_m2_k * inputs[0].surface_area_m2;
        let equilibrium = (2.0 * cb * 305.0 + cg * 300.0) / (2.0 * cb + cg);
        let mean_difference = 5.0 * (-conductance * (1.0 / cb + 2.0 / cg) * duration).exp();
        let half_difference = 15.0 * (-conductance / cb * duration).exp();
        let mean = equilibrium + cg / (2.0 * cb + cg) * mean_difference;
        let exact = [
            mean - half_difference,
            mean + half_difference,
            equilibrium - 2.0 * cb / (2.0 * cb + cg) * mean_difference,
        ];
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        for tolerance in [2e-5, 5e-6, 1.25e-6] {
            let c = controls(duration, tolerance);
            let result = environment.advance_adaptive(cx, &bodies, c).unwrap();
            let temperatures = [
                result.bodies()[0].state.temperature_k(),
                result.bodies()[1].state.temperature_k(),
                result.gas().state().state().temperature,
            ];
            errors.push(
                temperatures
                    .iter()
                    .zip(exact)
                    .map(|(a, b)| (a - b).abs())
                    .sum::<f64>(),
            );
            counts.push(result.report().samples.len());
            assert_eq!(
                result
                    .bodies()
                    .iter()
                    .map(|b| b.report.id)
                    .collect::<Vec<_>>(),
                vec![10, 20]
            );
            assert_eq!(result.report().samples.last().unwrap().time_s, duration);
            assert!(result.report().rejected_comparisons > 0);
            assert!(result.report().maximum_accepted_normalized_error <= 1.0);
            assert!(result.report().attempted_intervals <= c.time.maximum_attempts);
            assert!(
                result.report().sum_absolute_energy_defect_j <= c.time.energy_balance_tolerance_j
            );
            assert!(
                result.report().endpoint_absolute_energy_defect_j
                    <= c.time.energy_balance_tolerance_j
            );
            close(result.report().combined_energy_residual_j, 0.0, 1e-8);
            close(result.report().entropy.balance_residual_j_k, 0.0, 1e-13);
            close(
                result.report().entropy.storage_change_j_k[0],
                cb * ((temperatures[0] / 290.0).ln() + (temperatures[1] / 320.0).ln()),
                1e-12,
            );
            close(
                result.report().entropy.storage_change_j_k[1],
                cg * (temperatures[2] / 300.0).ln(),
                1e-12,
            );
            let mut last_gas = 300.0;
            for sample in &result.report().samples {
                assert!((sample.gas_temperature_k - last_gas).abs() <= 1.0);
                last_gas = sample.gas_temperature_k;
                assert_eq!(sample.bodies.len(), 2);
            }
        }
        assert!(
            errors.windows(2).all(|pair| pair[1] < pair[0]),
            "{errors:?}"
        );
        assert!(
            counts.windows(2).all(|pair| pair[1] > pair[0]),
            "{counts:?}"
        );
        let c = controls(duration, 2e-5);
        assert_eq!(
            environment.advance_adaptive(cx, &bodies, c).unwrap(),
            environment
                .advance_adaptive(cx, &[bodies[1], bodies[0]], c)
                .unwrap()
        );
        // Equal and opposite exchanges leave the gas exactly idle. Body error
        // must still subdivide; an error norm on only net gas heat would miss it.
        let opposing = [
            FiniteGasBody {
                input: UniformEnthalpyStepInput {
                    initial: input(&source, 310.0, 0.0).initial,
                    ..inputs[0]
                },
                ..bodies[0]
            },
            bodies[1],
        ];
        let balanced = environment.advance_adaptive(cx, &opposing, c).unwrap();
        assert_eq!(balanced.gas(), &gas);
        assert!(balanced.report().samples.len() > 2);
        assert!(balanced.report().rejected_comparisons > 0);
        assert!(balanced.bodies()[0].report.solid_to_gas_heat_j < 0.0);
        assert!(balanced.bodies()[1].report.solid_to_gas_heat_j > 0.0);
        assert_eq!(gas, saved);
    });
}

#[test]
fn g1_shared_phase_error_survives_an_insulated_second_body_and_constant_temperature() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            source.caloric.curve().material_card_identity(),
            [
                (-40.0, 280.0, 0.0),
                (0.0, 300.0, 0.0),
                (100.0, 300.0, 1.0),
                (180.0, 340.0, 1.0),
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
        let boundaries = [boundary(&source, 100.0), boundary(&source, 0.0)];
        let duration = 0.01;
        let latent = UniformEnthalpyStepInput {
            curve: &curve,
            initial: curve.state_at_specific_enthalpy(50.0).unwrap(),
            mass_kg: 0.001,
            volume_m3: 1e-6,
            surface_area_m2: 0.0013,
            internal_heat_j: 0.0,
            duration_s: duration,
        };
        let bodies = [
            FiniteGasBody {
                id: 1,
                environment: &boundaries[0],
                input: latent,
            },
            FiniteGasBody {
                id: 2,
                environment: &boundaries[1],
                input: input(&source, 310.0, 0.0),
            },
        ];
        let gas = gas(330.0, 1e-6);
        let environment = FiniteGasBodiesEnvironment {
            maximum_gas_temperature_increment_k: 100.0,
            exchange: exchange(1.0),
            ..shared(&gas, 0.0)
        };
        let cg = 2.5 * 100000.0 * gas.volume_m3() / 330.0;
        let exact = 0.5
            + cg * 30.0 * (1.0 - (-100.0 * latent.surface_area_m2 / cg * duration).exp())
                / (latent.mass_kg * 100.0);
        let loose = environment
            .advance_adaptive(cx, &bodies, controls(duration, 1000.0))
            .unwrap();
        let mut c = controls(duration, 1000.0);
        c.time.phase_fraction_absolute_tolerance = 1e-5;
        let tight = environment.advance_adaptive(cx, &bodies, c).unwrap();
        assert!(tight.report().samples.len() > loose.report().samples.len());
        assert!(
            (tight.bodies()[0].state.liquid_mass_fraction() - exact).abs()
                < (loose.bodies()[0].state.liquid_mass_fraction() - exact).abs()
        );
        assert!(
            tight
                .report()
                .samples
                .iter()
                .all(|s| s.bodies[0].state.temperature_k() == 300.0)
        );
        assert_eq!(tight.bodies()[1].state, bodies[1].input.initial);
        assert_eq!(tight.bodies()[1].report.solid_to_gas_heat_j, 0.0);
        close(
            tight.bodies()[0].report.entropy.storage_change_j_k,
            tight.bodies()[0].report.enthalpy_change_j / 300.0,
            1e-13,
        );
    });
}

#[test]
fn g3_wall_feedback_windows_preserve_original_halves_sources_and_actual_interval_admission() {
    scope(|cx, _| {
        let source = sources(2.0, 1000.0, false, false);
        let boundary = boundary(&source, 0.0);
        let duration = 0.001;
        let mut inputs = [input(&source, 320.0, 2e-4), input(&source, 290.0, -1e-4)];
        for input in &mut inputs {
            input.duration_s = duration;
        }
        let bodies = [
            FiniteGasBody {
                id: 2,
                environment: &boundary,
                input: inputs[0],
            },
            FiniteGasBody {
                id: 1,
                environment: &boundary,
                input: inputs[1],
            },
        ];
        let gas = gas(300.0, 1e-6);
        let environment = shared(&gas, 3e-4);
        let resolve = |body: &FiniteGasBody<'_>,
                       gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
                       wall_k: f64| {
            Ok::<_, Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: 80.0 + 0.5 * (wall_k - 300.0),
                source: (
                    body.id,
                    body.input.duration_s,
                    gas.state().temperature,
                    wall_k,
                ),
            })
        };
        let admissions = RefCell::new(Vec::new());
        let result = environment
            .advance_adaptive_with_wall_convection_checked(
                cx,
                &bodies,
                controls(duration, 1000.0),
                wall(),
                resolve,
                |body, state| {
                    assert!(core::ptr::eq(body.input.curve, source.caloric.curve()));
                    admissions.borrow_mut().push((
                        body.id,
                        body.input.duration_s,
                        body.input.initial,
                        state,
                    ));
                    Ok::<_, Infallible>(())
                },
            )
            .unwrap();
        assert_eq!(result.report().attempted_intervals, 3);
        assert_eq!(result.report().samples.len(), 2);
        assert_eq!(admissions.borrow().len(), 12);
        assert_eq!(
            admissions
                .borrow()
                .iter()
                .filter(|v| v.1 == duration)
                .count(),
            4
        );
        assert_eq!(
            admissions
                .borrow()
                .iter()
                .filter(|v| v.1 == duration / 2.0)
                .count(),
            8
        );

        let mut native_gas = gas.clone();
        let mut native_bodies = bodies;
        for body in &mut native_bodies {
            body.input.duration_s *= 0.5;
            body.input.internal_heat_j *= 0.5;
        }
        for sample in &result.report().samples {
            let native = FiniteGasBodiesEnvironment {
                gas: &native_gas,
                gas_heat_j: environment.gas_heat_j * 0.5,
                ..environment.clone()
            }
            .advance_with_wall_convection(cx, &native_bodies, wall(), resolve)
            .unwrap();
            assert_eq!(sample.bodies.as_slice(), native.bodies());
            assert_eq!(&sample.shared, native.report());
            assert_eq!(
                sample.gas_temperature_k,
                native.gas().state().state().temperature
            );
            for descriptor in &mut native_bodies {
                descriptor.input.initial = native
                    .bodies()
                    .iter()
                    .find(|b| b.report.id == descriptor.id)
                    .unwrap()
                    .state;
            }
            native_gas = native.into_parts().0;
            for body in &sample.bodies {
                assert_eq!(body.report.convection.source.source.1, duration / 2.0);
                assert_eq!(
                    body.report.convection.source.endpoint_source.1,
                    duration / 2.0
                );
            }
        }
        assert_eq!(result.gas(), &native_gas);
        close(
            result.report().external_gas_heat_j,
            environment.gas_heat_j,
            1e-16,
        );
        close(result.report().internal_body_heat_j, 1e-4, 1e-16);
        for body in result.bodies() {
            let original = bodies.iter().find(|b| b.id == body.report.id).unwrap();
            close(
                body.report.internal_heat_j,
                original.input.internal_heat_j,
                1e-16,
            );
            close(body.report.heat_apportionment_residual_j, 0.0, 1e-16);
            close(
                body.report.enthalpy_change_j,
                original.input.internal_heat_j - body.report.solid_to_gas_heat_j,
                1e-9,
            );
        }
        close(result.report().combined_energy_residual_j, 0.0, 1e-8);
        close(result.report().entropy.balance_residual_j_k, 0.0, 1e-13);
    });
}

#[test]
fn g4_shared_window_budgets_late_source_refusal_and_cancellation_are_atomic() {
    scope(|cx, gate| {
        let source = sources(2.0, 1000.0, false, false);
        let boundary = boundary(&source, 100.0);
        let duration = 0.001;
        let mut inputs = [input(&source, 320.0, 0.0), input(&source, 290.0, 0.0)];
        for input in &mut inputs {
            input.duration_s = duration;
        }
        let bodies = [
            FiniteGasBody {
                id: 1,
                environment: &boundary,
                input: inputs[0],
            },
            FiniteGasBody {
                id: 2,
                environment: &boundary,
                input: inputs[1],
            },
        ];
        let gas = gas(300.0, 1e-6);
        let saved = gas.clone();
        let environment = shared(&gas, 0.0);
        let mut c = controls(duration, 1000.0);
        c.time.maximum_step_s = duration / 4.0;
        let reference = environment.advance_adaptive(cx, &bodies, c).unwrap();
        let work = FiniteGasBodiesWindowSolve {
            time: FiniteGasWindowSolve {
                maximum_attempts: 4,
                ..c.time
            },
            ..c
        };
        assert!(matches!(environment.advance_adaptive(cx, &bodies, work),
            Err(FiniteGasBodiesWindowError::WorkLimit { attempted_intervals: 4, accepted_time_s })
                if accepted_time_s == duration / 4.0));
        assert!(matches!(
            environment.advance_adaptive(
                cx,
                &bodies,
                FiniteGasBodiesWindowSolve {
                    maximum_retained_body_steps: 4,
                    ..c
                }
            ),
            Err(FiniteGasBodiesWindowError::HistoryLimit {
                retained_body_steps: 4,
                ..
            })
        ));
        let mut strict = controls(duration, 1e-25);
        strict.time.minimum_step_s = duration;
        assert!(matches!(
            environment.advance_adaptive(cx, &bodies, strict),
            Err(FiniteGasBodiesWindowError::TimeAccuracy {
                accepted_time_s: 0.0,
                ..
            })
        ));
        let refused = environment.advance_adaptive_with_convection(cx, &bodies, c, |body, _| {
            if body.id == 2 && body.input.initial.temperature_k() > 290.01 {
                Err("second material domain ended")
            } else {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            }
        });
        assert!(matches!(
            refused,
            Err(FiniteGasBodiesWindowError::Step(
                FiniteGasBodiesError::Body {
                    id: 2,
                    error: ConvectionError::Convection("second material domain ended")
                }
            ))
        ));
        assert_eq!(
            environment.advance_adaptive(cx, &bodies, c).unwrap(),
            reference
        );
        let source_gate = environment.advance_adaptive_with_wall_convection_checked(
            cx,
            &bodies,
            c,
            wall(),
            |_, _, _| {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
            |body, state| {
                if body.id == 2 && state != body.input.initial {
                    Err("original endpoint conductivity unsupported")
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            source_gate,
            Err(FiniteGasBodiesWindowError::Step(
                FiniteGasBodiesError::Body {
                    id: 2,
                    error: ConvectionError::Convection(WallConvectionError::Source(
                        "original endpoint conductivity unsupported"
                    ))
                }
            ))
        ));
        let cancelled = environment.advance_adaptive_with_wall_convection_checked(
            cx,
            &bodies,
            c,
            wall(),
            |_, _, _| {
                Ok::<_, Infallible>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 100.0,
                    source: (),
                })
            },
            |body, state| {
                if body.id == 2 && state != body.input.initial {
                    gate.request();
                }
                Ok(())
            },
        );
        assert!(matches!(
            cancelled,
            Err(FiniteGasBodiesWindowError::Cancelled)
        ));
        assert_eq!(gas, saved);
        for (body, input) in bodies.iter().zip(inputs) {
            assert_eq!(body.input.initial, input.initial);
        }
    });
}
