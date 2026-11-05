use super::*;
use crate::enthalpy::conduction::implicit::adaptive::{
    AdaptiveImplicitConductionControls, AdaptiveImplicitConductionError,
};

fn source_event(time_s: f64, power_w: f64, bath: Option<f64>) -> TimedEvent {
    let mut boundary = boundary(bath);
    boundary.source_w = power_w;
    TimedEvent {
        time_s,
        event: Event::Boundary(id(0), boundary),
    }
}

fn adaptive_controls(dt: f64) -> AdaptiveImplicitConductionControls {
    AdaptiveImplicitConductionControls {
        step: controls(dt),
        initial_dt_s: dt,
        enthalpy_absolute_tolerance_j_kg: 1e-5,
        temperature_absolute_tolerance_k: 1e-5,
        relative_tolerance: 0.0,
        ..Default::default()
    }
}

#[test]
fn source_pulses_hit_exact_events_and_target_without_renewing_heat() {
    let law = Law::constant(3.0, 1.0);
    let (mut thermal, observer) = network(
        &law,
        &[300.0],
        &[2.0],
        &[],
        None,
        vec![
            source_event(0.35, -6.0, None),
            source_event(0.8, 30.0, None),
            source_event(0.8, 18.0, None),
            source_event(1.2, 50.0, None),
            source_event(1.3, -7.0, None),
        ],
    );
    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 12.0;
    let before = thermal.cells().clone();
    let direct = chain(&thermal, &observer, &law);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // Two simultaneous events retain authored order; the target-time source
    // delivers no heat before the target. These are physical interval integrals.
    let joules = 12.0 * 0.35 - 6.0 * 0.45 + 18.0 * 0.4;
    for (sparse, adaptive) in [(false, false), (true, false), (false, true), (true, true)] {
        let (next, report) = match (sparse, adaptive) {
            (false, false) => direct
                .advance_sdirk2_to(&thermal, 1.2, controls(1.0), &mut workspace, || Ok(()))
                .unwrap(),
            (true, false) => graph
                .advance_sdirk2_to(&thermal, 1.2, controls(1.0), &mut workspace, || Ok(()))
                .unwrap(),
            (false, true) => direct
                .advance_sdirk2_adaptive_to(
                    &thermal,
                    1.2,
                    adaptive_controls(1.0),
                    &mut workspace,
                    || Ok(()),
                )
                .map(|(s, r)| (s, r.solve))
                .unwrap(),
            (true, true) => graph
                .advance_sdirk2_adaptive_to(
                    &thermal,
                    1.2,
                    adaptive_controls(1.0),
                    &mut workspace,
                    || Ok(()),
                )
                .map(|(s, r)| (s, r.solve))
                .unwrap(),
        };
        near(
            next.temperature_k(id(0)).unwrap(),
            300.0 + joules / 6.0,
            1e-10,
        );
        near(report.energy.external_heat_j, joules, 1e-12);
        near(report.energy.stored_energy_change_j, joules, 1e-9);
        assert_eq!(report.energy.events, 4);
        assert_eq!(next.cells()[&id(0)].boundary.source_w, 50.0);
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(next.time_s(), 1.2);
        let (continued, continuation) = direct
            .advance_sdirk2_to(&next, 1.25, controls(1.0), &mut workspace, || Ok(()))
            .unwrap();
        near(continuation.energy.external_heat_j, 2.5, 1e-12);
        near(
            continued.temperature_k(id(0)).unwrap(),
            300.0 + (joules + 2.5) / 6.0,
            1e-10,
        );
        assert_eq!(continuation.energy.events, 0);
    }
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn driven_bath_solution_refines_across_source_discontinuities() {
    let law = Law::constant(2.0, 1.0);
    let (mut thermal, observer) = network(
        &law,
        &[400.0],
        &[1.0],
        &[],
        Some(300.0),
        vec![
            source_event(0.35, -10.0, Some(300.0)),
            source_event(0.8, 50.0, Some(300.0)),
        ],
    );
    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 20.0;
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // Independently solve C*dT/dt=P-G*(T-Tb), C=2 J/K and G=1 W/K,
    // over each constant-power physical interval.
    let mut exact = 400.0;
    for (duration, power) in [(0.35_f64, 20.0), (0.45, -10.0), (0.4, 50.0)] {
        let equilibrium = 300.0 + power;
        exact = equilibrium + (exact - equilibrium) * (-duration / 2.0).exp();
    }
    let mut errors = Vec::new();
    for dt in [0.2, 0.1, 0.05] {
        let (next, report) = solver
            .advance_sdirk2_to(&thermal, 1.2, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let temperature = next.temperature_k(id(0)).unwrap();
        errors.push((temperature - exact).abs());
        near(
            report.energy.external_heat_j,
            2.0 * (temperature - 400.0),
            1e-8,
        );
        assert_eq!(report.energy.events, 2);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.5..4.5).contains(&(pair[0] / pair[1])),
            "source-event errors: {errors:?}"
        );
    }
    let (next, report) = solver
        .advance_sdirk2_adaptive_to(
            &thermal,
            1.2,
            adaptive_controls(0.2),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), exact, 4e-4);
    assert_eq!(report.solve.energy.events, 2);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
}

#[test]
fn late_bath_or_topology_changes_refuse_before_any_source_interval() {
    let law = Law::constant(1.0, 1.0);
    for unsupported in [
        Event::Boundary(id(0), boundary(Some(350.0))),
        Event::Separation(0),
    ] {
        let (thermal, observer) = network(
            &law,
            &[400.0, 300.0],
            &[1.0, 1.0],
            &[(0, 1, 1.0)],
            Some(300.0),
            vec![
                source_event(0.25, 10.0, Some(300.0)),
                TimedEvent {
                    time_s: 1.0,
                    event: unsupported,
                },
            ],
        );
        let before = thermal.cells().clone();
        let solver = chain(&thermal, &observer, &law);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        // One evaluation cannot solve the first interval. Invalid, rather than
        // a work-budget failure, proves the entire window was admitted first.
        let step = ImplicitConductionControls {
            maximum_residual_evaluations: 1,
            ..controls(0.1)
        };
        assert!(matches!(
            solver.advance_sdirk2_to(&thermal, 1.0, step, &mut workspace, || Ok(())),
            Err(ImplicitConductionError::Invalid(_))
        ));
        assert!(matches!(
            solver.advance_sdirk2_adaptive_to(
                &thermal,
                1.0,
                AdaptiveImplicitConductionControls {
                    step,
                    ..adaptive_controls(0.1)
                },
                &mut workspace,
                || Ok(())
            ),
            Err(AdaptiveImplicitConductionError::Implicit(
                ImplicitConductionError::Invalid(_)
            ))
        ));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn source_events_keep_whole_call_step_budgets_and_cancellable_preflight() {
    let law = Law::constant(1.0, 1.0);
    let (thermal, observer) = network(
        &law,
        &[300.0],
        &[1.0],
        &[],
        None,
        vec![
            source_event(0.25, 2.0, None),
            source_event(0.5, -1.0, None),
            source_event(0.75, 4.0, None),
        ],
    );
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let step = ImplicitConductionControls {
        maximum_steps: 2,
        ..controls(1.0)
    };
    assert!(matches!(
        solver.advance_sdirk2_to(&thermal, 1.0, step, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(_)))
    ));
    assert!(matches!(
        solver.advance_sdirk2_adaptive_to(
            &thermal,
            1.0,
            AdaptiveImplicitConductionControls {
                step,
                ..adaptive_controls(1.0)
            },
            &mut workspace,
            || Ok(())
        ),
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Budget(_))
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&id(0)].boundary.source_w, 0.0);

    let (thermal, observer) = network(
        &law,
        &[300.0],
        &[1.0],
        &[],
        None,
        (1..=300)
            .map(|i| source_event(i as f64 / 512.0, 1.0, None))
            .collect(),
    );
    let solver = chain(&thermal, &observer, &law);
    let mut checks = 0;
    let result = solver.advance_sdirk2_adaptive_to(
        &thermal,
        1.0,
        adaptive_controls(1.0),
        &mut workspace,
        || {
            checks += 1;
            if checks == 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(checks, 2);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&id(0)].boundary.source_w, 0.0);
}
