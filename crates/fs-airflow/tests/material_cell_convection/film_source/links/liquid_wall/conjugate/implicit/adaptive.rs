//! Physical-time accuracy for the original nonlinear finite wall/fluid pair.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::thermal::adaptive::{
    AdaptiveHeatAssemblyControls, AdaptiveHeatAssemblyError,
    advance_heat_assembly_adaptive_candidate,
};
use fs_thermal::enthalpy::conduction::implicit::adaptive::AdaptiveImplicitConductionError;
use std::cell::RefCell;

fn adaptive_controls(tolerance: f64) -> AdaptiveHeatAssemblyControls {
    let mut step = super::controls();
    step.maximum_evaluations = 30_000;
    step.newton.absolute_tolerance = 1e-10;
    AdaptiveHeatAssemblyControls {
        step,
        initial_dt_s: 0.2,
        minimum_dt_s: 1e-6,
        maximum_dt_s: 0.2,
        maximum_steps: 200,
        maximum_halvings: 24,
        maximum_trials: 500,
        enthalpy_absolute_tolerance_j_kg: CP * tolerance,
        temperature_absolute_tolerance_k: tolerance,
        phase_fraction_absolute_tolerance: Some(1e-6),
        relative_tolerance: 0.0,
    }
}

// Independent continuous energy reduction: Cwall*dTw/dt = -Q(Tw,Tfluid),
// Tfluid follows exact conservation. Integrate dt = Cwall/Q dTw by Simpson
// quadrature; no production time integrator or endpoint solve is called.
fn continuous_endpoint(time: f64) -> [f64; 2] {
    let elapsed = |wall: f64| {
        let n = 256;
        let width = (310.0 - wall) / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let t = wall + width * i as f64;
            let bulk = 290.0 + CP / FLUID_CP * (310.0 - t);
            let reciprocal = MASS * CP / surface_oracle(t, bulk, 1.0).1;
            sum += reciprocal
                * if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
        }
        sum * width / 3.0
    };
    let (mut lo, mut hi) = (294.000001, 310.0);
    for _ in 0..55 {
        let wall = 0.5 * lo + 0.5 * hi;
        if elapsed(wall) > time {
            lo = wall;
        } else {
            hi = wall;
        }
    }
    let wall = 0.5 * lo + 0.5 * hi;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn tighter_temporal_accuracy_improves_actual_cooling_and_retains_only_fine_heat() {
    let target = 0.2;
    let reference = continuous_endpoint(target);
    let mut errors = Vec::new();
    let mut steps = Vec::new();
    for tolerance in [0.05, 0.0005] {
        with_fixture(1.0, false, |cx, _, prepared, owner, liquid, thermal| {
            let saved = thermal.cells().clone();
            let assembly = SourceConjugateHeatAssembly {
                heat: prepared,
                links: owner,
                layer_config: config(),
                resolve_fluid: |_: &SourceConvectiveLink,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    liquid.at_temperature(t)
                },
            };
            let (next, report) = advance_heat_assembly_adaptive_candidate(
                &assembly,
                cx,
                thermal,
                target,
                adaptive_controls(tolerance),
                &[],
            )
            .unwrap();
            let wall = next.temperature_k(IDS[0]).unwrap();
            let bulk = next.temperature_k(IDS[1]).unwrap();
            assert_eq!(next.time_s(), target);
            assert_eq!(
                report.time.solve.energy.substeps,
                2 * report.time.accepted_comparisons
            );
            assert!(report.time.accuracy_rejections > 0);
            assert!(report.time.maximum_accepted_normalized_error <= 1.0);
            assert!(report.time.solve.residual_evaluations <= 30_000);
            near(next.energy_j(), thermal.energy_j(), 1e-7);
            near(report.time.solve.energy.external_heat_j, 0.0, 1e-12);
            near(report.time.solve.energy.residual_j, 0.0, 1e-7);
            let entropy = &report.time.solve.entropy;
            near(
                entropy.thermal_storage_change_j_k,
                MASS * CP * (wall / 310.0).ln() + MASS * FLUID_CP * (bulk / 290.0).ln(),
                1e-10,
            );
            assert!(entropy.conduction_generation_j_k >= 0.0);
            assert!(!entropy.has_unclassified_sources);
            let sample = &report.source.conjugate.samples()[0];
            assert_eq!(sample.wall_temperature_k, wall);
            assert_eq!(sample.fluid_temperature_k, bulk);
            near(
                sample.layer.surface_temperature_k,
                surface_oracle(wall, bulk, 1.0).0,
                3e-9,
            );
            errors.push((wall - reference[0]).abs().max((bulk - reference[1]).abs()));
            steps.push(report.time.solve.energy.substeps);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.time_s(), 0.0);
        });
    }
    assert!(
        errors[1] < 0.3 * errors[0],
        "tighter local accuracy must improve the actual solution"
    );
    assert!(errors[1] < 0.03);
    assert!(
        steps[1] > steps[0],
        "accuracy must cause real physical-time refinement"
    );
}

fn source_event(time_s: f64, power: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Boundary(
            IDS[0],
            Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: power,
            },
        ),
    }
}

#[test]
fn source_events_land_exactly_and_endpoint_evidence_follows_the_target_event() {
    let events = [
        source_event(0.07, 1.5),
        source_event(0.13, 0.0),
        source_event(0.2, 0.25),
    ];
    with_scheduled_fixture(
        1.0,
        false,
        &events,
        |cx, _, prepared, owner, liquid, thermal| {
            let queries = RefCell::new(Vec::new());
            let assembly = SourceConjugateHeatAssembly {
                heat: prepared,
                links: owner,
                layer_config: config(),
                resolve_fluid: |_: &SourceConvectiveLink,
                                state: &DepositionSimulation<'_>,
                                t: Temperature| {
                    queries
                        .borrow_mut()
                        .push((state.time_s(), state.cells()[&IDS[0]].boundary.source_w));
                    liquid.at_temperature(t)
                },
            };
            let (next, report) = advance_heat_assembly_adaptive_candidate(
                &assembly,
                cx,
                thermal,
                0.2,
                adaptive_controls(0.002),
                &[0.035, 0.16],
            )
            .unwrap();
            let expected = 1.5 * (0.13 - 0.07);
            near(next.energy_j() - thermal.energy_j(), expected, 1e-7);
            near(report.time.solve.energy.external_heat_j, expected, 1e-12);
            near(report.time.solve.energy.residual_j, 0.0, 1e-7);
            assert_eq!(report.time.solve.energy.events, 3);
            assert!(report.time.solve.entropy.has_unclassified_sources);
            assert_eq!(next.cells()[&IDS[0]].boundary.source_w, 0.25);
            let queries = queries.borrow();
            for time in [0.035, 0.07, 0.13, 0.16, 0.2] {
                assert!(
                    queries.iter().any(|&(seen, _)| seen == time),
                    "missed exact transition {time}"
                );
            }
            assert!(
                queries
                    .iter()
                    .any(|&(time, power)| time == 0.07 && power == 0.0)
            );
            assert!(
                queries
                    .iter()
                    .any(|&(time, power)| time == 0.13 && power == 1.5)
            );
            assert_eq!(queries.last().copied(), Some((0.2, 0.25)));
            assert_eq!(
                report.source.conjugate.samples()[0].wall_temperature_k,
                next.temperature_k(IDS[0]).unwrap()
            );
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(thermal.cells()[&IDS[0]].boundary.source_w, 0.0);
        },
    );
}

#[test]
fn work_accuracy_late_source_and_cancellation_refusals_publish_no_prefix() {
    with_fixture(1.0, false, |cx, gate, prepared, owner, liquid, thermal| {
        let saved = thermal.cells().clone();
        let assembly = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |_: &SourceConvectiveLink,
                            _: &DepositionSimulation<'_>,
                            t: Temperature| liquid.at_temperature(t),
        };
        let mut budget = adaptive_controls(0.002);
        budget.step.maximum_evaluations = 1;
        assert!(matches!(
            advance_heat_assembly_adaptive_candidate(&assembly, cx, thermal, 0.2, budget, &[],),
            Err(AdaptiveHeatAssemblyError::Solve(
                CircuitEnthalpyError::CircuitBudget
            ))
        ));
        let mut accuracy = adaptive_controls(1e-10);
        accuracy.initial_dt_s = 0.05;
        accuracy.minimum_dt_s = 0.05;
        accuracy.maximum_dt_s = 0.05;
        assert!(matches!(
            advance_heat_assembly_adaptive_candidate(&assembly, cx, thermal, 0.2, accuracy, &[],),
            Err(AdaptiveHeatAssemblyError::Time(
                AdaptiveImplicitConductionError::TimeAccuracy { .. }
            ))
        ));
        let mut retry = adaptive_controls(0.5);
        retry.initial_dt_s = 0.05;
        retry.maximum_dt_s = 0.05;
        let saw_advanced_trial = Cell::new(false);
        let limited = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |_: &SourceConvectiveLink,
                            state: &DepositionSimulation<'_>,
                            t: Temperature| {
                if state.time_s() > 0.11 {
                    saw_advanced_trial.set(true);
                    return Err("later original source refused");
                }
                liquid
                    .at_temperature(t)
                    .map_err(|_| "liquid source refused")
            },
        };
        assert!(matches!(
            advance_heat_assembly_adaptive_candidate(&limited, cx, thermal, 0.2, retry, &[],),
            Err(AdaptiveHeatAssemblyError::Solve(
                CircuitEnthalpyError::Source(SourceConjugateHeatError::Conjugate(
                    SourceConjugateError::Convection(SourceCellConvectionError::Source {
                        source: "later original source refused",
                        ..
                    })
                ))
            ))
        ));
        assert!(saw_advanced_trial.get());
        // Reuse the same prepared original owners after all rejected work.
        let (valid, _) =
            advance_heat_assembly_adaptive_candidate(&assembly, cx, thermal, 0.1, retry, &[])
                .unwrap();
        assert_eq!(valid.time_s(), 0.1);
        let interrupted = SourceConjugateHeatAssembly {
            heat: prepared,
            links: owner,
            layer_config: config(),
            resolve_fluid: |_: &SourceConvectiveLink,
                            state: &DepositionSimulation<'_>,
                            t: Temperature| {
                if state.time_s() > 0.11 {
                    gate.request();
                }
                liquid.at_temperature(t)
            },
        };
        assert!(matches!(
            advance_heat_assembly_adaptive_candidate(&interrupted, cx, thermal, 0.2, retry, &[],),
            Err(AdaptiveHeatAssemblyError::Solve(
                CircuitEnthalpyError::Cancelled
            ))
        ));
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
