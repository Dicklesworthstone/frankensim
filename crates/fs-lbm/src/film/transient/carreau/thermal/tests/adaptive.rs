//! Independent native heat/work and event references for the adaptive owner.

use super::super::adaptive::{
    ThermalFilmAdaptiveControls, ThermalFilmAdaptiveError, ThermalFilmQuery,
};
use super::*;

fn native_controls() -> ThermalFilmControls {
    ThermalFilmControls {
        maximum_endpoint_evaluations: 100_000,
        maximum_momentum_material_evaluations: 50_000_000,
        ..controls()
    }
}
fn accuracy(tolerance: f64) -> ThermalFilmAdaptiveControls {
    ThermalFilmAdaptiveControls {
        initial_dt_s: 0.08,
        minimum_dt_s: 1e-6,
        maximum_dt_s: 0.08,
        velocity_absolute_tolerance_m_s: 1e-6,
        temperature_absolute_tolerance_k: 1e-5,
        energy_absolute_tolerance_j: tolerance,
        transport_absolute_tolerance_m2: 1e-6,
        relative_tolerance: 0.0,
        window_energy_tolerance_j: 2e-8,
        ..ThermalFilmAdaptiveControls::default()
    }
}
fn stationary<'a>(thermal: DepositionSimulation<'a>) -> ThermalCarreauFilm<'a> {
    let traces = [
        thermal.temperature_k(IDS[0]).unwrap(),
        thermal.temperature_k(IDS[1]).unwrap(),
    ];
    let flow =
        TransientCarreauFilm::from_profile(source(), traces[0], SPEC, vec![0.0; 2], 0.0).unwrap();
    ThermalCarreauFilm::new(flow, thermal, &IDS, 1.0, traces, || false).unwrap()
}

fn resolution_inventory<'a>(
    law: &'a Caloric,
    temperatures: [f64; 2],
    conductance: f64,
    steps: usize,
) -> DepositionSimulation<'a> {
    // Retain the genuine original cells/links, with an explicit native work
    // allowance for these accuracy tests. A separate test exercises refusal
    // at the native ceiling; adaptive controls cannot override that ceiling.
    let original = inventory(law, temperatures, 1.0, conductance, 0.0);
    let events = original
        .cells()
        .values()
        .cloned()
        .map(Event::Birth)
        .chain(original.links().values().cloned().map(Event::Contact))
        .map(|event| TimedEvent { time_s: 0.0, event })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: steps,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn heat_and_opposing_actuator_work_control_a_state_exact_ramp() {
    // Manufactured exact state: u=.25+.5*y, T=300+10*t. The original
    // thermal Carreau stress changes, while momentum and BE storage are exact.
    // Additional heat is the independent required storage rate minus the
    // analytic uniform-shear loss. A state-only time norm cannot resolve heat.
    let law = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let make = || {
        let flow =
            TransientCarreauFilm::from_profile(source(), 300.0, SPEC, vec![0.375, 0.625], 0.0)
                .unwrap();
        ThermalCarreauFilm::new(
            flow,
            resolution_inventory(&law, [300.0; 2], 0.0, 4096),
            &IDS,
            1.0,
            [300.0; 2],
            || false,
        )
        .unwrap()
    };
    let viscosity = |t: f64| {
        source()
            .evaluate(300.0 + 10.0 * t, 0.5)
            .unwrap()
            .dynamic_viscosity_pa_s()
    };
    let forcing = |q: ThermalFilmQuery| {
        Ok::<_, DepositionError>(FilmDrive {
            wall_velocity_m_s: 0.25,
            surface_shear_pa: 0.5 * viscosity(q.endpoint_s),
            acceleration_m_s2: 0.0,
        })
    };
    let exchange = |q: ThermalFilmQuery, thermal: &DepositionSimulation<'_>, cells: &[CellId]| {
        let mut exchange = insulated(thermal, cells)?;
        let loss_w_per_cell = 0.5 * 0.25 * viscosity(q.endpoint_s);
        exchange.heat.bulk.sources_w = cells
            .iter()
            .map(|id| (*id, 0.5 * law.cp * 10.0 - loss_w_per_cell))
            .collect();
        Ok::<_, DepositionError>(exchange)
    };
    let mut loose = make();
    let mut tight = make();
    let coarse = loose
        .advance_adaptive(
            0.08,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            accuracy(5e-6),
            forcing,
            exchange,
            || false,
        )
        .unwrap();
    let fine = tight
        .advance_adaptive(
            0.08,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            accuracy(5e-8),
            forcing,
            exchange,
            || false,
        )
        .unwrap();
    // Independent Simpson quadrature of the ORIGINAL constitutive evaluations,
    // not the adaptive owner's sampled heat or a fitted apparent viscosity.
    let n = 2048;
    let dt = 0.08 / n as f64;
    let integral = dt / 3.0
        * (0..=n)
            .map(|i| {
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * viscosity(i as f64 * dt)
            })
            .sum::<f64>();
    let expected_heat = 0.25 * integral;
    assert!(fine.accuracy_rejections > 0 && fine.accepted_steps > 2);
    assert!(fine.accepted_steps > coarse.accepted_steps);
    assert!(
        (fine.viscous_heat_j - expected_heat).abs()
            < (coarse.viscous_heat_j - expected_heat).abs() / 3.0
    );
    near(fine.viscous_heat_j, expected_heat, 5e-6);
    near(fine.wall_work_j, -0.125 * integral, 2.5e-6);
    near(fine.surface_work_j, 0.375 * integral, 7.5e-6);
    assert!(fine.wall_work_j < 0.0 && fine.surface_work_j > 0.0);
    for (i, id) in IDS.iter().enumerate() {
        near(tight.flow().velocities_m_s()[i], [0.375, 0.625][i], 1e-10);
        near(tight.thermal().temperature_k(*id).unwrap(), 300.8, 1e-8);
        near(fine.deposited_cell_heat_j[i], 0.5 * expected_heat, 2.5e-6);
        near(
            fine.other_cell_heat_j[id] + fine.deposited_cell_heat_j[i],
            0.02,
            1e-9,
        );
    }
    near(fine.stored_energy_change_j, 0.04, 1e-9);
    near(fine.kinetic_energy_change_j, 0.0, 1e-10);
    assert!(fine.combined_energy_residual_j.abs() < 2e-8);
    assert!(fine.maximum_accepted_normalized_error <= 1.0);
    assert!(fine.suggested_dt_s >= accuracy(5e-8).minimum_dt_s);
}

#[test]
fn original_endpoint_link_override_matches_two_cell_exponential_and_cell_heat() {
    let law = Caloric {
        cp: 2.0,
        datum: 1000.0,
    };
    // Each cell capacity is 1 J/K; actual endpoint override is G=.5 W/K.
    // Thus mean T=300 and half-difference=10 exp(-t). The authored link
    // intentionally differs: using base-candidate rates would fail this test.
    let mut model = stationary(resolution_inventory(&law, [310.0, 290.0], 7.0, 4096));
    let exchange = |_: ThermalFilmQuery, thermal: &DepositionSimulation<'_>, cells: &[CellId]| {
        let mut result = insulated(thermal, cells)?;
        result.heat.bulk.conductances_w_k.insert(91, 0.5);
        Ok::<_, DepositionError>(result)
    };
    let report = model
        .advance_adaptive(
            0.2,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            accuracy(1e-5),
            |_| Ok(FilmDrive::default()),
            exchange,
            || false,
        )
        .unwrap();
    let half_difference = 10.0 * fs_math::det::exp(-0.2);
    near(
        model.thermal().temperature_k(IDS[0]).unwrap(),
        300.0 + half_difference,
        0.002,
    );
    near(
        model.thermal().temperature_k(IDS[1]).unwrap(),
        300.0 - half_difference,
        0.002,
    );
    for (id, initial) in [(IDS[0], 310.0), (IDS[1], 290.0)] {
        let delta = model.thermal().temperature_k(id).unwrap() - initial;
        near(report.other_cell_heat_j[&id], delta, 2e-8);
    }
    assert!(report.other_cell_heat_j[&IDS[0]] < 0.0 && report.other_cell_heat_j[&IDS[1]] > 0.0);
    near(report.other_heat_j, 0.0, 1e-15);
    near(report.deposited_viscous_heat_j, 0.0, 1e-15);
    near(report.stored_energy_change_j, 0.0, 2e-8);
    assert!(report.combined_energy_residual_j.abs() < 2e-8);
}

fn scheduled_inventory<'a>(law: &'a Caloric, birth: bool) -> DepositionSimulation<'a> {
    let base = inventory(law, [300.0; 2], 1.0, 0.0, 0.0);
    let mut events: Vec<_> = base
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if birth {
        let mut extra = base.cells()[&IDS[0]].clone();
        extra.id = CellId(99);
        extra.segment_id = 99;
        extra.body_id = 99;
        extra.geometry.start_m = [4.0, 0.0, 0.0];
        extra.geometry.end_m = [4.0, 0.5, 0.0];
        events.push(TimedEvent {
            time_s: 0.037,
            event: Event::Birth(extra),
        });
    } else {
        events.push(TimedEvent {
            time_s: 0.037,
            event: Event::Contact(ThermalLink {
                id: 91,
                a: IDS[0],
                b: IDS[1],
                area_m2: 1.0,
                conductance_w_k: 0.5,
            }),
        });
        for (time, watts) in [(0.037, 1.0), (0.074, -2.0)] {
            if time == 0.074 {
                events.push(TimedEvent {
                    time_s: time,
                    event: Event::Separation(91),
                });
            }
            for id in IDS {
                let mut boundary = base.cells()[&id].boundary;
                boundary.source_w = watts;
                events.push(TimedEvent {
                    time_s: time,
                    event: Event::Boundary(id, boundary),
                });
            }
        }
        let mut boundary = base.cells()[&IDS[0]].boundary;
        boundary.source_w = 3.0;
        events.push(TimedEvent {
            time_s: 0.111,
            event: Event::Boundary(IDS[0], boundary),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn decimal_source_and_native_events_preserve_left_limits_and_original_contacts() {
    let law = Caloric {
        cp: 2.0,
        datum: 1000.0,
    };
    let mut model = stationary(scheduled_inventory(&law, false));
    let events = [0.037, 0.074, 0.111];
    let adaptive = ThermalFilmAdaptiveControls {
        initial_dt_s: 0.037,
        maximum_dt_s: 0.037,
        ..accuracy(1e-6)
    };
    let source = |q: ThermalFilmQuery, thermal: &DepositionSimulation<'_>, cells: &[CellId]| {
        assert!(q.start_s < q.midpoint_s && q.midpoint_s < q.endpoint_s);
        let mut result = insulated(thermal, cells)?;
        let rate = [0.5, 1.0, -0.5, 0.0][q.segment];
        result.heat.bulk.sources_w = cells.iter().map(|id| (*id, rate)).collect();
        Ok::<_, DepositionError>(result)
    };
    let first = model
        .advance_adaptive(
            0.111,
            &events,
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            adaptive,
            |_| Ok(FilmDrive::default()),
            source,
            || false,
        )
        .unwrap();
    assert_eq!(model.thermal().time_s(), 0.111);
    assert_eq!(model.flow().accepted_time_s(), 0.111);
    assert_eq!(model.thermal().next_event_time_s(), Some(0.111));
    assert_eq!(first.accepted_steps, 6);
    assert_eq!(first.native_events_processed, 6);
    assert_eq!(first.suggested_dt_s, 0.037);
    assert!(model.thermal().links().is_empty());
    for id in IDS {
        near(model.thermal().temperature_k(id).unwrap(), 300.0, 1e-9);
    }
    near(first.other_heat_j, 0.0, 1e-14);
    let zero = model
        .advance_adaptive(
            0.111,
            &events,
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            adaptive,
            |_| Ok(FilmDrive::default()),
            source,
            || false,
        )
        .unwrap();
    assert_eq!(zero.accepted_steps, 0);
    assert_eq!(model.thermal().next_event_time_s(), Some(0.111));
    let second = model
        .advance_adaptive(
            0.148,
            &events,
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            adaptive,
            |_| Ok(FilmDrive::default()),
            source,
            || false,
        )
        .unwrap();
    assert_eq!(second.native_events_processed, 1);
    assert_eq!(model.thermal().next_event_time_s(), None);
    near(
        model.thermal().temperature_k(IDS[0]).unwrap(),
        300.0 + 3.0 * 0.037,
        1e-9,
    );
    near(
        model.thermal().temperature_k(IDS[1]).unwrap(),
        300.0 - 2.0 * 0.037,
        1e-9,
    );
    near(second.other_heat_j, 0.037, 1e-14);
}

#[test]
fn nonfilm_birth_source_failure_budget_and_late_cancellation_are_atomic() {
    let law = Caloric {
        cp: 2.0,
        datum: 1000.0,
    };
    let adaptive = ThermalFilmAdaptiveControls {
        initial_dt_s: 0.037,
        maximum_dt_s: 0.037,
        ..accuracy(1e-6)
    };
    let make = |birth| stationary(scheduled_inventory(&law, birth));
    let mut birth = make(true);
    let before = birth.thermal().cells().clone();
    let result = birth.advance_adaptive(
        0.074,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        native_controls(),
        adaptive,
        |_| Ok(FilmDrive::default()),
        |_, t, ids| insulated(t, ids),
        || false,
    );
    assert!(matches!(
        result,
        Err(ThermalFilmAdaptiveError::Film(ThermalFilmError::Invalid(_)))
    ));
    assert_eq!(birth.thermal().cells(), &before);
    assert_eq!(birth.thermal().next_event_time_s(), Some(0.037));
    assert_eq!(birth.flow().accepted_time_s(), 0.0);
    let mut failed = make(false);
    let before = failed.thermal().cells().clone();
    let result = failed.advance_adaptive(
        0.074,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        native_controls(),
        adaptive,
        |q| {
            if q.start_s >= 0.037 {
                Err(DepositionError::InvalidInput("authored source refusal"))
            } else {
                Ok(FilmDrive::default())
            }
        },
        |_, t, ids| insulated(t, ids),
        || false,
    );
    assert!(matches!(
        result,
        Err(ThermalFilmAdaptiveError::Film(ThermalFilmError::Source(_)))
    ));
    assert_eq!(failed.thermal().cells(), &before);
    assert_eq!(failed.thermal().next_event_time_s(), Some(0.037));
    let result = failed.advance_adaptive(
        0.074,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        ThermalFilmControls {
            maximum_endpoint_evaluations: 1,
            ..native_controls()
        },
        adaptive,
        |_| Ok(FilmDrive::default()),
        |_, t, ids| insulated(t, ids),
        || false,
    );
    assert!(matches!(
        result,
        Err(ThermalFilmAdaptiveError::Film(
            ThermalFilmError::EndpointBudget
        ))
    ));
    assert_eq!(failed.thermal().cells(), &before);
    let mut successful = make(false);
    let mut polls = 0;
    successful
        .advance_adaptive(
            0.074,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            adaptive,
            |_| Ok(FilmDrive::default()),
            |_, t, ids| insulated(t, ids),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled_polls = 0;
    let result = failed.advance_adaptive(
        0.074,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        native_controls(),
        adaptive,
        |_| Ok(FilmDrive::default()),
        |_, t, ids| insulated(t, ids),
        || {
            cancelled_polls += 1;
            cancelled_polls == polls
        },
    );
    assert!(matches!(
        result,
        Err(ThermalFilmAdaptiveError::Film(ThermalFilmError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(failed.thermal().cells(), &before);
    assert_eq!(failed.flow().velocities_m_s(), &[0.0; 2]);
    failed
        .advance_adaptive(
            0.074,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            adaptive,
            |_| Ok(FilmDrive::default()),
            |_, t, ids| insulated(t, ids),
            || false,
        )
        .unwrap();
    assert_eq!(failed.thermal().cells(), successful.thermal().cells());
    assert_eq!(
        failed.flow().velocities_m_s(),
        successful.flow().velocities_m_s()
    );
}

#[test]
fn native_substep_ceiling_and_nonuniform_zero_window_keep_original_state() {
    let law = Caloric {
        cp: 2.0,
        datum: 1000.0,
    };
    let mut model = stationary(resolution_inventory(&law, [310.0, 290.0], 0.0, 2));
    let before = model.thermal().cells().clone();
    let diffusion = model.flow().zero_shear_diffusion_time_s();
    let faces: Vec<_> = (0..=2)
        .map(|i| model.flow().material_at_face(i).unwrap().temperature_k())
        .collect();
    let zero = model
        .advance_adaptive(
            0.0,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native_controls(),
            accuracy(1e-6),
            |_| -> Result<FilmDrive, DepositionError> { panic!("zero window sampled a source") },
            |_, _, _| -> Result<ThermalFilmExchange<()>, DepositionError> {
                panic!("zero window sampled heat")
            },
            || false,
        )
        .unwrap();
    assert_eq!(zero.accepted_steps, 0);
    assert_eq!(
        model.flow().zero_shear_diffusion_time_s().to_bits(),
        diffusion.to_bits()
    );
    for (i, expected) in faces.iter().enumerate() {
        assert_eq!(
            model
                .flow()
                .material_at_face(i)
                .unwrap()
                .temperature_k()
                .to_bits(),
            expected.to_bits()
        );
    }
    let result = model.advance_adaptive(
        0.16,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        native_controls(),
        accuracy(1e-6),
        |_| Ok(FilmDrive::default()),
        |_, t, ids| insulated(t, ids),
        || false,
    );
    assert!(matches!(result, Err(ThermalFilmAdaptiveError::StepBudget)));
    assert_eq!(model.thermal().cells(), &before);
    assert_eq!(model.flow().accepted_time_s(), 0.0);
    assert_eq!(
        model.flow().zero_shear_diffusion_time_s().to_bits(),
        diffusion.to_bits()
    );
}
