//! Automatic nonlinear solves on the original event and mass timeline.
use super::*;

fn controls(maximum_dt_s: f64) -> EvaporationImplicitTimelineControls {
    EvaporationImplicitTimelineControls {
        maximum_dt_s,
        solve: EvaporationImplicitControls::new(range(0.0)),
    }
}

fn unchanged(actual: &DepositionSimulation<'_>, expected: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s, expected.time_s);
    assert_eq!(actual.cursor, expected.cursor);
    assert_eq!(actual.cells, expected.cells);
    assert_eq!(actual.links, expected.links);
    assert_eq!(actual.covered_areas, expected.covered_areas);
    assert_eq!(actual.schedule, expected.schedule);
    assert_eq!(actual.budget, expected.budget);
}

#[test]
fn event_clipped_timeline_matches_interval_composition_and_keeps_pre_event_payload() {
    let liquid = cell(1, 0.1, 300.0, 0.0, 10.0);
    let mut covered = liquid.boundary;
    covered.source_w = 3.0;
    covered.reservoir = Some(Reservoir {
        area_m2: 0.004,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let mut final_boundary = covered;
    final_boundary.source_w = 99.0;
    let old = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(liquid),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Birth(cell(2, 0.2, 310.0, 0.0, -2.0)),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.002,
                    conductance_w_k: 50.0,
                }),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), covered),
            },
            TimedEvent {
                time_s: 1.5,
                event: Event::Separation(9),
            },
            TimedEvent {
                time_s: 1.5,
                event: Event::Boundary(CellId(1), final_boundary),
            },
        ],
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 6,
            max_substeps: 3,
        },
    )
    .unwrap();
    let policy = controls(1.0);
    let (born, no_payload, birth_report) = old
        .advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            0.0,
            policy,
            empty,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(born.cells().len(), 1);
    assert!(no_payload.is_none());
    assert_eq!(birth_report.physical.events, 1);
    assert_eq!(birth_report.physical.substeps, 0);
    assert_eq!(birth_report.residual_evaluations, 0);

    let mut composed = old.advance_to(0.0, config(), || Ok(())).unwrap().0;
    let mut iterations = 0;
    let mut evaluations = 0;
    for endpoint in [0.5, 1.0, 1.5] {
        let (candidate, (), report) = composed
            .solve_implicit_evaporating_step(
                CellId(1),
                &CONSTANT,
                film(0.0),
                endpoint,
                policy.solve,
                empty,
                || Ok(()),
            )
            .unwrap();
        iterations += report.newton_iterations;
        evaluations += report.residual_evaluations;
        composed = candidate
            .advance_to(endpoint, config(), || Ok(()))
            .unwrap()
            .0;
    }
    let (next, payload, report) = old
        .advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            1.5,
            policy,
            |thermal, liquid| {
                Ok(coefficients((
                    thermal.time_s(),
                    liquid.area_m2,
                    thermal.links().len(),
                    thermal.cells()[&CellId(1)].boundary.source_w,
                )))
            },
            || Ok(()),
        )
        .unwrap();
    for id in [CellId(1), CellId(2)] {
        close(
            next.cells()[&id].mass_kg,
            composed.cells()[&id].mass_kg,
            1e-16,
        );
        close(
            next.cells()[&id].specific_h_j_kg,
            composed.cells()[&id].specific_h_j_kg,
            1e-7,
        );
        assert_eq!(next.cells()[&id].boundary, composed.cells()[&id].boundary);
    }
    let expected_transfer = [0.01, 0.008, 0.004]
        .into_iter()
        .map(|area| transfer(&CONSTANT, 300.0, film(0.0), area, 0.5))
        .sum::<f64>();
    let (endpoint, exposed, links, source) = payload.unwrap();
    assert_eq!(endpoint, 1.5);
    close(exposed, 0.004, 1e-16);
    assert_eq!(links, 1);
    assert_eq!(source, 3.0);
    assert!(next.links().is_empty());
    assert_eq!(next.cells()[&CellId(1)].boundary, final_boundary);
    close(next.exposed_area_m2(CellId(1)).unwrap(), 0.006, 1e-16);
    close(
        report.physical.transferred_mass_kg,
        expected_transfer,
        1e-16,
    );
    close(report.physical.incoming_mass_kg, 0.3, 1e-16);
    close(report.physical.incoming_energy_j, 8000.0, 1e-10);
    close(report.physical.external_heat_j, 9.5, 1e-10);
    close(
        next.energy_j() + report.physical.vapor_energy_j - 8000.0,
        9.5,
        1e-6,
    );
    assert_eq!(report.physical.events, 6);
    assert_eq!(report.physical.substeps, 3);
    assert_eq!(report.physical.phase_steps, 3);
    assert_eq!(report.newton_iterations, iterations);
    assert_eq!(report.residual_evaluations, evaluations);
    assert!(report.absolute_energy_residual_j <= policy.solve.total_energy_tolerance_j);
    assert!(report.absolute_mass_rate_residual_kg <= policy.solve.mass_tolerance_kg);
    assert_eq!(next.time_s(), 1.5);
    assert!(old.cells().is_empty());
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn nonlinear_sources_use_endpoint_time_and_proposed_mass_on_every_interval() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 2.0)], None, vec![]);
    let source = |time: f64, mass: f64, temperature: f64| {
        -1500.0 + 20_000.0 * mass + 8.0 * time + 4.0 * (temperature - 300.0).powi(2)
    };
    let raw_dm = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0);
    let mut expected_mass = 0.1;
    let mut expected_temperature = 300.0;
    let mut expected_heat = 0.0;
    for endpoint in [1.0, 2.0, 3.0] {
        let new_mass = expected_mass - raw_dm;
        let dm = expected_mass - new_mass;
        let (mut lo, mut hi) = (280.0, 330.0);
        for _ in 0..80 {
            let t = lo + 0.5 * (hi - lo);
            let residual = expected_mass * CP * (t - expected_temperature) + dm * LATENT
                - (2.0 + source(endpoint, new_mass, t));
            if residual > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        expected_temperature = lo + 0.5 * (hi - lo);
        expected_heat += 2.0 + source(endpoint, new_mass, expected_temperature);
        expected_mass = new_mass;
    }
    let (next, payload, report) = old
        .advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            3.0,
            controls(1.0),
            |thermal, liquid| {
                let time = thermal.time_s();
                assert!([1.0, 2.0, 3.0].contains(&time));
                let temperature = thermal.temperature_k(CellId(1))?;
                let original_mass = thermal.cells()[&CellId(1)].mass_kg;
                close(original_mass - liquid.mass_kg, raw_dm, 1e-16);
                let mut result = coefficients((time, temperature, original_mass, liquid.mass_kg));
                result
                    .transport
                    .heat
                    .bulk
                    .sources_w
                    .insert(CellId(1), source(time, liquid.mass_kg, temperature));
                Ok(result)
            },
            || Ok(()),
        )
        .unwrap();
    close(
        next.temperature_k(CellId(1)).unwrap(),
        expected_temperature,
        1e-8,
    );
    close(next.cells()[&CellId(1)].mass_kg, expected_mass, 1e-16);
    close(report.physical.external_heat_j, expected_heat, 1e-6);
    close(report.physical.energy_residual_j, 0.0, 1e-6);
    let (time, temperature, original_mass, proposed_mass) = payload.unwrap();
    assert_eq!(time, 3.0);
    assert_eq!(temperature, next.temperature_k(CellId(1)).unwrap());
    assert_eq!(proposed_mass, next.cells()[&CellId(1)].mass_kg);
    assert!(proposed_mass < original_mass);
    assert_eq!(
        next.cells()[&CellId(1)].boundary,
        old.cells()[&CellId(1)].boundary
    );
    assert_eq!(report.physical.substeps, 3);
    assert_eq!(report.physical.phase_steps, 3);
}

#[test]
fn numerical_and_thermal_budgets_cover_the_whole_timeline() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 2.0)], None, vec![]);
    let policy = controls(1.0);
    let (_, _, first) = old
        .advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            1.0,
            policy,
            empty,
            || Ok(()),
        )
        .unwrap();
    assert!(first.newton_iterations > 0);
    assert!(first.residual_evaluations > 1);
    for limit in 0..3 {
        let mut thermal = old.clone();
        let mut bounded = policy;
        match limit {
            0 => bounded.solve.maximum_evaluations = first.residual_evaluations,
            1 => bounded.solve.maximum_iterations = first.newton_iterations,
            _ => thermal.budget.max_substeps = 1,
        }
        let before = thermal.clone();
        let result = thermal.advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            2.0,
            bounded,
            empty,
            || Ok(()),
        );
        assert!(matches!(
            result,
            Err(EvaporationImplicitError::Physical(
                ThermalEvaporationError::Thermal(DepositionError::Budget(_))
            ))
        ));
        unchanged(&thermal, &before);
    }
}

#[test]
fn late_source_refusal_and_final_cancellation_discard_all_accepted_intervals() {
    let mut changed = cell(1, 0.1, 300.0, 0.0, 7.0).boundary;
    changed.ambient_k = 305.0;
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, 2.0)],
        None,
        vec![TimedEvent {
            time_s: 2.0,
            event: Event::Boundary(CellId(1), changed),
        }],
    );
    let before = old.clone();
    let mut seen = [false; 2];
    let refused = old.advance_implicit_evaporating_to(
        CellId(1),
        &CONSTANT,
        film(0.0),
        3.0,
        controls(1.0),
        |thermal, _| {
            match thermal.time_s() {
                1.0 => seen[0] = true,
                2.0 => seen[1] = true,
                _ => {
                    assert_eq!(thermal.time_s(), 3.0);
                    assert_eq!(thermal.cells()[&CellId(1)].boundary, changed);
                    return Err(ThermalEvaporationError::Surface(EvaporationError::Invalid(
                        "late coefficient refusal",
                    )));
                }
            }
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert_eq!(seen, [true; 2]);
    assert!(matches!(
        refused,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Surface(EvaporationError::Invalid("late coefficient refusal"))
        ))
    ));
    unchanged(&old, &before);
    let mut checkpoints = 0;
    let (_, _, completed) = old
        .advance_implicit_evaporating_to(
            CellId(1),
            &CONSTANT,
            film(0.0),
            3.0,
            controls(1.0),
            empty,
            || {
                checkpoints += 1;
                Ok(())
            },
        )
        .unwrap();
    assert_eq!(completed.physical.substeps, 3);
    let mut visited = 0;
    let cancelled = old.advance_implicit_evaporating_to(
        CellId(1),
        &CONSTANT,
        film(0.0),
        3.0,
        controls(1.0),
        empty,
        || {
            visited += 1;
            if visited == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(visited, checkpoints);
    assert!(matches!(
        cancelled,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::Cancelled)
        ))
    ));
    unchanged(&old, &before);
}

#[test]
fn complete_implicit_timeline_converges_to_the_variable_mass_heat_equation() {
    let mut liquid = cell(1, 0.1, 300.0, 0.0, 0.0);
    liquid.boundary.ambient_k = 310.0;
    liquid.boundary.htc_w_m2_k = 400.0;
    let old = simulation(&CALORIC, vec![liquid], None, vec![]);
    let duration = 200.0;
    let rate = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0);
    let conductance = 4.0;
    let equilibrium = 310.0 - rate * LATENT / conductance;
    // m(t) Cp dT/dt = G(T_ambient-T) - r L, with m(t) = m0-r t.
    let exact = equilibrium
        + (300.0 - equilibrium) * (1.0 - rate * duration / 0.1).powf(conductance / (CP * rate));
    assert!(100.0 > 0.8 * 0.1 * CP / conductance);
    let run = |dt: f64| {
        let (next, _, report) = old
            .advance_implicit_evaporating_to(
                CellId(1),
                &CONSTANT,
                film(0.0),
                duration,
                controls(dt),
                empty,
                || Ok(()),
            )
            .unwrap();
        close(next.mass_kg(), 0.1 - rate * duration, 1e-15);
        close(report.physical.energy_residual_j, 0.0, 1e-6);
        assert!(report.absolute_energy_residual_j < 1e-6);
        (next.temperature_k(CellId(1)).unwrap() - exact).abs()
    };
    let coarse = run(100.0);
    let medium = run(50.0);
    let fine = run(25.0);
    assert!(
        coarse > 1e-3,
        "the fixture must resolve physical timestep error"
    );
    assert!(medium < 0.65 * coarse, "{coarse:e}, {medium:e}");
    assert!(fine < 0.65 * medium, "{medium:e}, {fine:e}");
}
