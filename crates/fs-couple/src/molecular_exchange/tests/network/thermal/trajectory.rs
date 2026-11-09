//! Physical temporal checks for the original finite thermal/molecular owners.
use super::*;
use crate::molecular_exchange::network::{
    thermal::trajectory::{
        self as evolving, MolecularNetworkThermalAdaptiveControls,
        MolecularNetworkThermalTrajectoryCandidate, MolecularNetworkThermalTrajectoryError,
    },
    trajectory::{
        MolecularExchangeAdaptiveControls, MolecularNetworkAdaptiveControls,
        MolecularNetworkTrajectoryError,
    },
};
use std::convert::Infallible;

fn adaptive(h: f64) -> MolecularNetworkThermalAdaptiveControls {
    MolecularNetworkThermalAdaptiveControls {
        molecular: MolecularNetworkAdaptiveControls {
            states: MolecularExchangeAdaptiveControls {
                initial_step_s: h,
                minimum_step_s: 1e-8,
                maximum_step_s: h,
                maximum_attempts: 20_000,
                maximum_accepted_intervals: 10_000,
                relative_tolerance: 0.0,
                temperature_tolerance_k: 1e-3,
                pressure_tolerance_pa: 0.5,
                species_mass_tolerance_kg: 1e-9,
                velocity_tolerance_m_s: 1e-3,
            },
            stream_momentum_tolerance_n_s: 1e-7,
            stream_energy_tolerance_j: 1e-3,
            window_mass_tolerance_kg: 1e-10,
            window_momentum_tolerance_n_s: 1e-8,
            window_energy_tolerance_j: 1e-3,
        },
        body_temperature_tolerance_k: 1e-3,
        body_energy_tolerance_j: 1e-3,
        body_liquid_fraction_tolerance: 1e-5,
        heat_tolerance_j: 1e-3,
    }
}
fn fixed(h: f64) -> MolecularNetworkThermalAdaptiveControls {
    let mut a = adaptive(h);
    a.molecular.states.temperature_tolerance_k = 1e6;
    a.molecular.states.pressure_tolerance_pa = 1e12;
    a.molecular.states.species_mass_tolerance_kg = 1.0;
    a.molecular.states.velocity_tolerance_m_s = 1e6;
    a.molecular.stream_momentum_tolerance_n_s = 1e6;
    a.molecular.stream_energy_tolerance_j = 1e12;
    a.body_temperature_tolerance_k = 1e6;
    a.body_energy_tolerance_j = 1e12;
    a.body_liquid_fraction_tolerance = 1.0;
    a.heat_tolerance_j = 1e12;
    a
}
fn input(
    chart: &EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
    duration: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        duration_s: duration,
        ..wall(chart, temperature)
    }
}
fn run(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration: f64,
    groups: &[MolecularGasThermalBoundary<'_>],
    limits: MolecularNetworkThermalControls,
    adaptive: MolecularNetworkThermalAdaptiveControls,
) -> MolecularNetworkThermalTrajectoryCandidate {
    evolving::advance_candidate(
        cx,
        gases,
        links,
        duration,
        &solves(gases),
        groups,
        limits,
        adaptive,
    )
    .unwrap()
}
fn check_window(
    original: &BTreeMap<u64, MovingGasParcel>,
    groups: &[MolecularGasThermalBoundary<'_>],
    result: &MolecularNetworkThermalTrajectoryCandidate,
) {
    let gas: f64 = original
        .iter()
        .map(|(id, old)| {
            let new = &result.gases()[id];
            internal(new) - internal(old) + new.kinetic_energy_j() - old.kinetic_energy_j()
        })
        .sum();
    let walls: f64 = groups
        .iter()
        .flat_map(|group| group.bodies)
        .map(|body| {
            body.input.mass_kg
                * (result.bodies()[&body.id].specific_enthalpy_j_kg()
                    - body.input.initial.specific_enthalpy_j_kg())
        })
        .sum();
    near(gas + walls - result.step().external_heat_j, 0.0, 2e-5);
    near(walls, result.step().body_energy_change_j, 1e-10);
    assert!(result.step().maximum_accepted_error <= 1.0);
    assert_eq!(
        result.step().accepted_steps.len(),
        2 * result.step().accepted_intervals
    );
    assert_eq!(
        result.step().attempted_intervals,
        result.step().accepted_intervals + result.step().rejected_intervals
    );
}

#[test]
fn isolated_finite_gas_and_body_converge_to_the_independent_two_capacity_solution() {
    scope(|cx, _| {
        let original = gas(catalog(PURE, 0.0), T, P, [20.0, -4.0, 2.0]);
        let gases = BTreeMap::from([(2, original.clone())]);
        let chart = curve([1000.0; 2]);
        let environment = body(20.0);
        let duration = 1.0;
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: input(&chart, 700.0, duration),
        }];
        let groups = [boundary(2, &bodies)];
        let cg = original.parcel().mass_kg() * 1.5 * R / M;
        let cb = WALL_MASS * WALL_CP;
        let k = 20.0 * WALL_AREA;
        let equilibrium = (cg * T + cb * 700.0) / (cg + cb);
        let contrast = (700.0 - T) * (-k * (cg.recip() + cb.recip()) * duration).exp();
        let expected_gas = equilibrium - cb * contrast / (cg + cb);
        let expected_body = equilibrium + cg * contrast / (cg + cb);
        let errors = [0.5, 0.25, 0.125].map(|h| {
            let result = run(
                cx,
                &gases,
                &[],
                duration,
                &groups,
                thermal_controls([7.0, -2.0, 4.0]),
                fixed(h),
            );
            check_window(&gases, &groups, &result);
            (temperature(&result.gases()[&2]) - expected_gas).abs()
        });
        assert!(errors[1] < 0.6 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.6 * errors[1], "{errors:?}");
        let accepted = run(
            cx,
            &gases,
            &[],
            duration,
            &groups,
            thermal_controls([7.0, -2.0, 4.0]),
            adaptive(duration),
        );
        near(temperature(&accepted.gases()[&2]), expected_gas, 0.15);
        near(accepted.bodies()[&7].temperature_k(), expected_body, 0.015);
        assert!(accepted.step().rejected_intervals > 0);
        assert!(accepted.step().accepted_intervals > 1);
        assert_eq!(
            accepted.gases()[&2].momentum_kg_m_s(),
            original.momentum_kg_m_s()
        );
        assert_eq!(
            accepted.gases()[&2].parcel().mass_kg(),
            original.parcel().mass_kg()
        );
        check_window(&gases, &groups, &accepted);
        assert_eq!(gases[&2], original);
        assert_eq!(bodies[0].input.initial.temperature_k(), 700.0);
    });
}

#[test]
fn equal_opposite_wall_heat_is_resolved_even_when_the_gas_state_does_not_change() {
    scope(|cx, _| {
        let model = catalog(PURE, 0.0);
        let gases = BTreeMap::from([
            (2, gas(model.clone(), T, P, [0.0; 3])),
            (4, gas(model, 500.0, P, [0.0; 3])),
        ]);
        let chart = curve([1000.0; 2]);
        let environment = body(10.0);
        let duration = 1.0;
        let mut bodies = [
            FiniteGasBody {
                id: 9,
                environment: &environment,
                input: UniformEnthalpyStepInput {
                    internal_heat_j: 0.3,
                    ..input(&chart, 500.0, duration)
                },
            },
            FiniteGasBody {
                id: 3,
                environment: &environment,
                input: UniformEnthalpyStepInput {
                    internal_heat_j: -0.3,
                    ..input(&chart, 300.0, duration)
                },
            },
        ];
        let other = [FiniteGasBody {
            id: 20,
            environment: &environment,
            input: UniformEnthalpyStepInput {
                internal_heat_j: 0.11,
                ..input(&chart, 700.0, duration)
            },
        }];
        let groups = [
            boundary(2, &bodies),
            MolecularGasThermalBoundary {
                gas_heat_j: 0.7,
                ..boundary(4, &other)
            },
        ];
        let mut accuracy = fixed(duration);
        accuracy.heat_tolerance_j = 1e-3;
        let result = run(
            cx,
            &gases,
            &[],
            duration,
            &groups,
            thermal_controls([0.0; 3]),
            accuracy,
        );
        assert!(result.step().accepted_intervals > 1);
        assert!(result.step().rejected_intervals > 0);
        near(temperature(&result.gases()[&2]), T, 1e-8);
        let heat = &result.step().bodies;
        assert_eq!([heat[0].id, heat[1].id], [3, 9]);
        assert!(heat[0].solid_to_gas_heat_j < -1.0);
        assert!(heat[1].solid_to_gas_heat_j > 1.0);
        near(
            heat[0].solid_to_gas_heat_j + heat[1].solid_to_gas_heat_j,
            0.0,
            1e-8,
        );
        bodies.reverse();
        let reversed_groups = [
            MolecularGasThermalBoundary {
                gas_heat_j: 0.7,
                ..boundary(4, &other)
            },
            boundary(2, &bodies),
        ];
        let reordered = run(
            cx,
            &gases,
            &[],
            duration,
            &reversed_groups,
            thermal_controls([0.0; 3]),
            accuracy,
        );
        assert_eq!(result, reordered);
    });
}

#[test]
fn latent_fraction_controls_an_isothermal_wall_without_temperature_or_heat_estimates() {
    scope(|cx, _| {
        let chart = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([81; 32]),
            [
                (0.0, 300.0, 0.0),
                (100_000.0, 400.0, 0.0),
                (300_000.0, 400.0, 1.0),
                (400_000.0, 500.0, 1.0),
            ]
            .into_iter()
            .map(
                |(specific_enthalpy_j_kg, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k,
                    liquid_mass_fraction,
                    bulk_density_kg_m3: 1000.0,
                },
            )
            .collect(),
        )
        .unwrap();
        let original = gas(catalog(PURE, 0.0), 600.0, P, [0.0; 3]);
        let gases = BTreeMap::from([(2, original.clone())]);
        let environment = body(20.0);
        let duration = 1.0;
        let initial = chart.state_at_specific_enthalpy(100_000.0).unwrap();
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: UniformEnthalpyStepInput {
                curve: &chart,
                initial,
                mass_kg: WALL_MASS,
                volume_m3: WALL_MASS / 1000.0,
                surface_area_m2: WALL_AREA,
                internal_heat_j: 0.0,
                duration_s: duration,
            },
        }];
        let mut accuracy = fixed(duration);
        accuracy.body_liquid_fraction_tolerance = 1e-5;
        let result = run(
            cx,
            &gases,
            &[],
            duration,
            &[boundary(2, &bodies)],
            thermal_controls([0.0; 3]),
            accuracy,
        );
        let cg = original.parcel().mass_kg() * 1.5 * R / M;
        let energy = cg * 200.0 * (1.0 - (-20.0 * WALL_AREA * duration / cg).exp());
        let expected_fraction = energy / (WALL_MASS * 200_000.0);
        assert!(result.step().rejected_intervals > 0);
        assert!(result.step().accepted_intervals > 1);
        near(result.bodies()[&7].temperature_k(), 400.0, 0.0);
        near(
            result.bodies()[&7].liquid_mass_fraction(),
            expected_fraction,
            5e-4,
        );
        assert!(result.bodies()[&7].liquid_mass_fraction() > 0.01);
        for step in &result.step().accepted_steps {
            near(
                step.thermal[&2].bodies()[0].state.temperature_k(),
                400.0,
                0.0,
            );
        }
        check_window(&gases, &[boundary(2, &bodies)], &result);
    });
}

#[test]
fn composition_and_thermal_feedback_use_only_accepted_native_states_and_replay_exactly() {
    scope(|cx, _| {
        let model = catalog(MIXTURE, 0.0);
        let law = interface(&model, 2000.0);
        let only_a = model
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
            )
            .unwrap();
        let only_b = model
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap(),
            )
            .unwrap();
        let gases = BTreeMap::from([
            (2, gas(only_a, T, P, [1.0, -2.0, 0.0])),
            (4, gas(only_b, 600.0, 1.2 * P, [-2.0, 0.0, 3.0])),
        ]);
        let saved = gases.clone();
        let chart = curve([1000.0; 2]);
        let environment = body(20.0);
        let duration = 0.2;
        let original = input(&chart, 800.0, duration);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &environment,
            input: original,
        }];
        let groups = [boundary(2, &bodies)];
        let mut links = [link(5, [2, 4], &law)];
        links[0].surface.velocity_world_m_s = [0.5, 0.2, -0.3];
        let limits = thermal_controls([7.0, -2.0, 4.0]);
        let result = run(
            cx,
            &gases,
            &links,
            duration,
            &groups,
            limits,
            adaptive(duration),
        );
        assert!(result.step().accepted_intervals > 1);
        assert!(
            result.gases()[&2]
                .parcel()
                .species_mass_kg(&SpeciesId::new("B").unwrap())
                .unwrap()
                > 0.0
        );
        assert!(
            result.gases()[&4]
                .parcel()
                .species_mass_kg(&SpeciesId::new("A").unwrap())
                .unwrap()
                > 0.0
        );
        assert!(result.bodies()[&7].temperature_k() < 800.0);
        let mut replay = gases.clone();
        let mut body_state = original.initial;
        for expected in &result.step().accepted_steps {
            let bodies = [FiniteGasBody {
                id: 7,
                environment: &environment,
                input: UniformEnthalpyStepInput {
                    initial: body_state,
                    duration_s: expected.duration_s,
                    ..original
                },
            }];
            let actual = thermal::step_candidate(
                cx,
                &replay,
                &links,
                expected.duration_s,
                &solves(&replay),
                &[boundary(2, &bodies)],
                limits,
            )
            .unwrap();
            assert_eq!(actual.step(), expected);
            let (next, states, _) = actual.into_parts();
            replay = next;
            body_state = states[&7];
        }
        assert_eq!(result.gases(), &replay);
        assert_eq!(result.bodies()[&7], body_state);
        check_window(&gases, &groups, &result);
        for total in &result.step().links {
            for side in 0..2 {
                let mass: f64 = result
                    .step()
                    .accepted_steps
                    .iter()
                    .flat_map(|step| &step.molecular.links)
                    .filter(|row| row.id == total.id)
                    .flat_map(|row| &row.exchange.species)
                    .map(|row| row.mass_kg[side])
                    .sum();
                near(total.gross[side].mass_kg, mass, 1e-15);
            }
        }
        let defect = result
            .step()
            .debit
            .energy_j
            .max(result.step().debit.laboratory_energy_j);
        assert!(defect > 0.0);
        let mut tight = adaptive(duration);
        tight.molecular.window_energy_tolerance_j = 0.5 * defect;
        assert!(matches!(
            evolving::advance_candidate(
                cx,
                &gases,
                &links,
                duration,
                &solves(&gases),
                &groups,
                limits,
                tight
            ),
            Err(MolecularNetworkThermalTrajectoryError::Adaptive(
                MolecularNetworkTrajectoryError::WindowBalance { .. }
            ))
        ));
        assert_eq!(gases, saved);
        assert_eq!(bodies[0].input.initial, original.initial);
    });
}

#[test]
fn prescribed_heat_is_a_whole_window_power_and_does_not_multiply_with_substeps() {
    scope(|cx, _| {
        let gases = BTreeMap::from([(2, gas(catalog(PURE, 0.0), T, P, [0.0; 3]))]);
        let chart = curve([1000.0; 2]);
        let insulated = body(0.0);
        let duration = 2.0;
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &insulated,
            input: UniformEnthalpyStepInput {
                internal_heat_j: 10.0,
                ..input(&chart, 600.0, duration)
            },
        }];
        let groups = [MolecularGasThermalBoundary {
            gas_heat_j: 3.0,
            ..boundary(2, &bodies)
        }];
        let result = run(
            cx,
            &gases,
            &[],
            duration,
            &groups,
            thermal_controls([0.0; 3]),
            fixed(0.37),
        );
        assert!(result.step().accepted_intervals > 1);
        near(result.step().external_heat_j, 13.0, 1e-12);
        near(result.step().source_partition_residual_j, 0.0, 1e-12);
        near(result.step().bodies[0].internal_heat_j, 10.0, 1e-12);
        near(result.step().bodies[0].storage_change_j, 10.0, 1e-7);
        near(
            internal(&result.gases()[&2]) - internal(&gases[&2]),
            3.0,
            1e-7,
        );
        check_window(&gases, &groups, &result);
        let mut exhausted = fixed(0.37);
        exhausted.molecular.states.maximum_accepted_intervals = 1;
        assert!(matches!(
            evolving::advance_candidate(
                cx,
                &gases,
                &[],
                duration,
                &solves(&gases),
                &groups,
                thermal_controls([0.0; 3]),
                exhausted
            ),
            Err(MolecularNetworkThermalTrajectoryError::Adaptive(
                MolecularNetworkTrajectoryError::WorkBudget {
                    accepted_intervals: 1,
                    ..
                }
            ))
        ));
        assert_eq!(bodies[0].input.initial.temperature_k(), 600.0);
        near(temperature(&gases[&2]), T, 0.0);
    });
}

#[test]
fn genuine_gross_withdrawal_retries_but_callback_refusal_and_cancellation_are_atomic() {
    scope(|cx, gate| {
        let model = catalog(PURE, 0.0);
        let law = interface(&model, 2000.0);
        let original = gas(model, T, P, [0.0; 3]);
        let gases = BTreeMap::from([(2, original.clone()), (4, original)]);
        let saved = gases.clone();
        let chart = curve([1000.0; 2]);
        let insulated = body(0.0);
        let bodies = [FiniteGasBody {
            id: 7,
            environment: &insulated,
            input: input(&chart, 600.0, DT),
        }];
        let groups = [boundary(2, &bodies)];
        let links = [link(1, [2, 4], &law)];
        let mut limits = thermal_controls([0.0; 3]);
        limits.exchange.maximum_species_fraction = 0.0004;
        let accepted = run(cx, &gases, &links, DT, &groups, limits, fixed(DT));
        assert!(accepted.step().rejected_intervals > 0);
        assert!(accepted.step().accepted_intervals > 1);
        check_window(&gases, &groups, &accepted);

        let mut calls = 0;
        let refused = evolving::advance_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            limits,
            fixed(DT),
            |_, _, _, _| {
                calls += 1;
                Err::<ConvectionCoefficient<()>, _>("unsupported original law")
            },
        );
        assert_eq!(calls, 1);
        assert!(matches!(
            refused,
            Err(MolecularNetworkThermalTrajectoryError::Step {
                error: MolecularNetworkThermalError::Thermal {
                    error: FiniteGasBodiesError::Body {
                        error: ConvectionError::Convection("unsupported original law"),
                        ..
                    },
                    ..
                },
                ..
            })
        ));
        assert_eq!(gases, saved);
        assert_eq!(bodies[0].input.initial.temperature_k(), 600.0);

        let cancelled = evolving::advance_candidate_with_convection(
            cx,
            &gases,
            &links,
            DT,
            &solves(&gases),
            &groups,
            limits,
            fixed(DT),
            |_, _, _, _| {
                gate.request();
                Ok::<_, Infallible>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: 0.0,
                    source: (),
                })
            },
        );
        assert!(matches!(
            cancelled,
            Err(MolecularNetworkThermalTrajectoryError::Adaptive(
                MolecularNetworkTrajectoryError::Cancelled { .. }
            ))
        ));
        assert_eq!(gases, saved);
        assert_eq!(bodies[0].input.initial.temperature_k(), 600.0);
    });
}
