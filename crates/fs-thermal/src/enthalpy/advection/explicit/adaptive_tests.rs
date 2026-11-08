use super::*;

fn adaptive(dt: f64, tolerance: f64) -> AdaptiveSsprk2Controls {
    AdaptiveSsprk2Controls {
        step: config(dt),
        initial_dt_s: dt,
        minimum_dt_s: 1e-10,
        maximum_trials: 10_000,
        enthalpy_absolute_tolerance_j_kg: tolerance,
        temperature_absolute_tolerance_k: Some(1e-10),
        relative_tolerance: 0.0,
    }
}

fn source(
    state: &DepositionSimulation<'_>,
    datum: f64,
    target_h: f64,
) -> Result<BoundaryStepCoefficients<()>, DepositionError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::from([(
            CellId(0),
            4.0 * (target_h - (state.cells()[&CellId(0)].specific_h_j_kg - datum)),
        )]),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into())
}

#[test]
fn adaptive_latent_ring_refines_enthalpy_even_when_every_temperature_is_constant() {
    let law = Latent;
    let n = 8;
    let theta = std::f64::consts::TAU / n as f64;
    let initial = ring(
        &law,
        &(0..n)
            .map(|i| 50.0 + 20.0 * (i as f64 * theta).cos())
            .collect::<Vec<_>>(),
    );
    let target = 0.2;
    for rate in [-1.0, 1.0] {
        let owner = flow(n, rate);
        let mut errors = Vec::new();
        for tolerance in [1e-2, 1e-5] {
            let (next, accepted, report) = owner
                .advance_ssprk2_adaptive_to(&initial, target, adaptive(0.1, tolerance), || Ok(()))
                .unwrap();
            let attenuation = (n as f64 * (theta.cos() - 1.0) * target).exp();
            let mut error = 0.0;
            for (i, cell) in next.cells().values().enumerate() {
                let exact = 50.0
                    + 20.0
                        * attenuation
                        * (i as f64 * theta - rate * n as f64 * theta.sin() * target).cos();
                error += (cell.specific_h_j_kg - exact).powi(2);
                assert_eq!(law.temperature_k(cell.specific_h_j_kg).unwrap(), 300.0);
            }
            errors.push((error / n as f64).sqrt());
            assert!(report.accuracy_rejections > 0);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(
                report.last_error.unwrap().maximum_temperature_difference_k,
                Some(0.0)
            );
            assert_eq!(report.work.energy.substeps, 2 * report.accepted_intervals);
            assert_eq!(accepted.audit().substeps, report.work.energy.substeps);
            near(
                accepted.audit().transported_mass_kg,
                n as f64 * target,
                2e-13,
            );
            near(report.work.energy.residual_j, 0.0, 2e-12);
        }
        assert!(errors[1] < errors[0] / 20.0, "latent refinement {errors:?}");
    }
}

#[test]
fn adaptive_nonlinear_source_control_ignores_formation_enthalpy_datum() {
    struct Shifted(f64);
    impl ThermalLaw for Shifted {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            Latent.temperature_k(h - self.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.1
        }
    }
    let owner = Advected::new((), vec![], 0).unwrap();
    let exact = 80.0 - 60.0 * (-2.0_f64).exp();
    let mut reference = Vec::new();
    for datum in [0.0, 65536.0] {
        let law = Shifted(datum);
        let initial = prepare(
            &law,
            births([cell(0, datum + 20.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0)]),
            1,
        );
        let mut errors = Vec::new();
        for (index, tolerance) in [0.02, 2e-5].into_iter().enumerate() {
            let mut controls = adaptive(0.5, tolerance);
            controls.relative_tolerance = tolerance / 100.0;
            controls.step.step.energy_tolerance_j = 1e-7;
            let (next, accepted, report) = owner
                .advance_ssprk2_adaptive_with_coefficients(
                    &initial,
                    0.5,
                    controls,
                    |state| source(state, datum, 80.0),
                    || Ok(()),
                )
                .unwrap();
            let h = next.cells()[&CellId(0)].specific_h_j_kg - datum;
            errors.push((h - exact).abs());
            assert!(report.accuracy_rejections > 0);
            assert_eq!(report.work.rejected_intervals, 0);
            assert_eq!(
                report.work.stage_evaluations,
                report.accepted_intervals + 4 * report.trials
            );
            assert_eq!(
                report.last_error.unwrap().maximum_temperature_difference_k,
                Some(0.0)
            );
            near(report.work.energy.external_heat_j, h - 20.0, 2e-9);
            near(accepted.audit().external_enthalpy_j, 0.0, 0.0);
            if datum == 0.0 {
                reference.push((h, report.trials));
            } else {
                near(h, reference[index].0, 2e-8);
                assert_eq!(report.trials, reference[index].1);
            }
        }
        assert!(errors[1] < errors[0] / 20.0, "source refinement {errors:?}");
    }
}

#[test]
fn adaptive_temperature_norm_uses_the_original_sensible_caloric_branch() {
    let law = Latent;
    let initial = prepare(
        &law,
        births([cell(0, 110.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0)]),
        1,
    );
    let owner = Advected::new((), vec![], 0).unwrap();
    let exact = 200.0 - 90.0 * (-0.8_f64).exp();
    let mut errors = Vec::new();
    for temperature in [None, Some(1e-5)] {
        let mut controls = adaptive(0.2, 1e6);
        controls.temperature_absolute_tolerance_k = temperature;
        let (next, _, report) = owner
            .advance_ssprk2_adaptive_with_coefficients(
                &initial,
                0.2,
                controls,
                |state| source(state, 0.0, 200.0),
                || Ok(()),
            )
            .unwrap();
        errors.push((next.cells()[&CellId(0)].specific_h_j_kg - exact).abs());
        if temperature.is_some() {
            assert!(report.accuracy_rejections > 0);
            assert!(
                report
                    .last_error
                    .unwrap()
                    .maximum_temperature_difference_k
                    .unwrap()
                    > 0.0
            );
        } else {
            assert_eq!(report.accepted_intervals, 1);
            assert_eq!(report.accuracy_rejections, 0);
        }
    }
    assert!(
        errors[1] < errors[0] / 100.0,
        "temperature refinement {errors:?}"
    );
}

#[test]
fn adaptive_pump_births_and_forcing_credit_only_the_accepted_fine_timeline() {
    let law = Latent;
    let first = cell(0, 20.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0);
    let mut changed = first.boundary;
    changed.source_w = 3.0;
    let mut events = births([first]);
    events.extend([
        TimedEvent {
            time_s: 0.05,
            event: Event::Boundary(CellId(0), changed),
        },
        TimedEvent {
            time_s: 0.05,
            event: Event::Birth(cell(1, 30.0, [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], 2.0)),
        },
    ]);
    let initial = prepare(&law, events, 2);
    let owner = startup();
    let mut queries = Vec::new();
    let (next, accepted, report) = owner
        .advance_ssprk2_adaptive_with_coefficients(
            &initial,
            0.1,
            adaptive(0.1, 1e-4),
            |state| {
                queries.push((
                    state.time_s(),
                    state.cells().len(),
                    state.cells()[&CellId(0)].boundary.source_w,
                ));
                let mut coefficients = unchanged(state)?;
                if state.cells().contains_key(&CellId(1)) {
                    coefficients
                        .bulk
                        .sources_w
                        .insert(CellId(1), 2.0 * state.time_s());
                }
                Ok(coefficients)
            },
            || Ok(()),
        )
        .unwrap();
    assert!(report.work.rejected_intervals > 0);
    assert!(report.accuracy_rejections > 0);
    assert!(queries.contains(&(0.05, 1, 0.0)));
    assert!(queries.contains(&(0.05, 2, 3.0)));
    assert!(queries.iter().all(|&(t, cells, _)| if t < 0.05 {
        cells == 1
    } else if t > 0.05 {
        cells == 2
    } else {
        true
    }));
    near(accepted.audit().inlet_mass_kg, 5.0, 3e-13);
    near(accepted.audit().outlet_mass_kg, 5.0, 3e-13);
    near(accepted.audit().transported_mass_kg, 10.0, 6e-13);
    near(
        report.work.energy.external_heat_j - accepted.audit().external_enthalpy_j,
        0.1575,
        1e-11,
    );
    near(next.cells()[&CellId(1)].specific_h_j_kg, 30.00375, 1e-12);
    near(report.work.energy.incoming_energy_j, 60.0, 0.0);
    near(report.work.energy.incoming_mass_kg, 2.0, 0.0);
    near(report.work.energy.residual_j, 0.0, 1e-11);
    assert_eq!(report.work.energy.events, 2);
    assert_eq!(report.work.energy.substeps, 2 * report.accepted_intervals);
    assert_eq!(accepted.audit().substeps, report.work.energy.substeps);
    assert_eq!(initial.cells().len(), 1);
    assert_eq!(owner.audit(), AdvectionAudit::default());
}

#[test]
fn adaptive_stage_trials_floor_and_cancellation_are_shared_and_transactional() {
    let law = Latent;
    let initial = prepare(
        &law,
        births([cell(0, 20.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0)]),
        1,
    );
    let owner = Advected::new((), vec![], 0).unwrap();
    let controls = adaptive(0.5, 1e-5);
    let mut stage = controls;
    stage.step.maximum_stage_evaluations = 4;
    let mut queries = 0;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_with_coefficients(
            &initial,
            0.5,
            stage,
            |state| {
                queries += 1;
                source(state, 0.0, 80.0)
            },
            || Ok(())
        ),
        Err(AdaptiveSsprk2Error::Deposition(DepositionError::Budget(
            "SSPRK2 stage evaluations"
        )))
    ));
    assert_eq!(queries, 4);
    let mut trial = controls;
    trial.maximum_trials = 1;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_with_coefficients(
            &initial,
            0.5,
            trial,
            |state| source(state, 0.0, 80.0),
            || Ok(())
        ),
        Err(AdaptiveSsprk2Error::Deposition(DepositionError::Budget(
            "adaptive SSPRK2 trials"
        )))
    ));
    let mut floor = controls;
    floor.minimum_dt_s = 0.5;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_with_coefficients(
            &initial,
            0.5,
            floor,
            |state| source(state, 0.0, 80.0),
            || Ok(())
        ),
        Err(AdaptiveSsprk2Error::TimeAccuracy { .. })
    ));
    let mut queries = 0;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_with_coefficients(
            &initial,
            0.5,
            controls,
            |state| {
                queries += 1;
                if queries == 5 {
                    Err(DepositionError::Cancelled)
                } else {
                    source(state, 0.0, 80.0)
                }
            },
            || Ok(())
        ),
        Err(AdaptiveSsprk2Error::Deposition(DepositionError::Cancelled))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 20.0);
    assert_eq!(owner.audit(), AdvectionAudit::default());
}
