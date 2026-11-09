//! Authored pump histories on the original radial/enthalpy time controller.

use super::*;

fn imposed(pressure: f64, head: f64, injection: f64) -> RegionalCarreauDrive {
    RegionalCarreauDrive {
        boundaries: pressures(pressure).to_vec(),
        injections: vec![FlowInjection {
            node: 0,
            flow: VolumetricFlowRate::new(injection),
        }],
        branch_heads: Some(vec![Pressure::new(head)]),
        branch_head_laws: None,
        branch_heat_destinations: None,
    }
}

// Reservoir and injection ports remain separate actual mass streams. In
// particular, an inward flow port must not erase an outward reservoir stream
// at the same junction just because their signed volume rates nearly cancel.
fn stream_exchange(
    thermal: &DepositionSimulation<'_>,
    flow: &TransientCarreauStep,
    drive: &RegionalCarreauDrive,
    inlet: f64,
) -> CarreauTransportExchange<f64> {
    let mut boundaries = Vec::new();
    let mut append = |id, node, q: f64| {
        if q > 0.0 {
            boundaries.push(MassBoundary::Inlet {
                id,
                cell: CellId(node as u64),
                kg_s: RHO * q,
                material_id: MaterialId(1),
                specific_h_j_kg: CP * (inlet - 300.0),
            });
        } else if q < 0.0 {
            boundaries.push(MassBoundary::Outlet {
                id,
                cell: CellId(node as u64),
                kg_s: -RHO * q,
            });
        }
    };
    for (i, (port, q)) in drive
        .boundaries
        .iter()
        .zip(&flow.reservoir_outflows)
        .enumerate()
    {
        append(100 + i as u64, port.node, q.value());
    }
    for (i, port) in drive.injections.iter().enumerate() {
        append(1000 + i as u64, port.node, port.flow.value());
    }
    CarreauTransportExchange {
        mass_boundaries: boundaries,
        heat: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: thermal.temperature_k(CellId(20)).unwrap(),
        }
        .into(),
    }
}

#[test]
fn g1_constant_drive_is_identical_to_held_entry_and_queries_once_per_native_attempt() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory_with_substeps(&law, &[300.0; 3], 2048);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let held = imposed(4000.0, 300.0, 2e-6);
        let calls = std::cell::RefCell::new(Vec::new());
        let sampled = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &[],
            controls(),
            resolution(),
            |_, time| {
                calls.borrow_mut().push(time);
                Ok::<_, ()>(held.clone())
            },
            |_, candidate, flow, sample| {
                assert_eq!(candidate.time_s(), sample.end_s);
                assert_eq!(flow.duration_s, sample.end_s - sample.start_s);
                assert_eq!(sample.sample_time_s, sample.start_s + 0.5 * flow.duration_s);
                Ok::<_, ()>(stream_exchange(candidate, flow, &sample.drive, 330.0))
            },
        )
        .unwrap();
        let reference = advance_resolved(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &held.boundaries,
            &held.injections,
            held.branch_heads.as_deref(),
            controls(),
            resolution(),
            |_, candidate, flow| Ok::<_, ()>(stream_exchange(candidate, flow, &held, 330.0)),
        )
        .unwrap();
        assert_eq!(snapshot(&sampled.momentum), snapshot(&reference.momentum));
        assert_eq!(sampled.thermal.cells(), reference.thermal.cells());
        assert_eq!(sampled.report.work, reference.report.work);
        assert_eq!(sampled.report.drive_samples, reference.report.drive_samples);
        assert_eq!(calls.borrow().len(), sampled.report.work.native_steps);
        assert_eq!(calls.borrow().len(), sampled.report.work.drive_evaluations);
        assert_eq!(
            sampled.report.drive_samples.len(),
            sampled.report.steps.len()
        );
        assert!(calls.borrow().len() > sampled.report.drive_samples.len());
        assert!(
            sampled
                .report
                .steps
                .iter()
                .map(|r| r.thermal.endpoint_evaluations)
                .sum::<usize>()
                > sampled.report.drive_samples.len()
        );
    });
}

fn discontinuous(time: f64) -> RegionalCarreauDrive {
    if time < 0.037 {
        imposed(4000.0, 200.0, 1e-6)
    } else if time < 0.074 {
        imposed(-6000.0, -300.0, -2e-6)
    } else {
        imposed(1000.0, 100.0, 3e-6)
    }
}

#[test]
fn g1_authored_nonbinary_events_split_both_comparisons_and_match_piecewise_native_drives() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let samples = std::cell::RefCell::new(Vec::new());
        let driven = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.111,
            &[-1.0, 0.0, 0.037, 0.074, 0.111, 1.0],
            controls(),
            permissive(),
            |_, time| {
                samples.borrow_mut().push(time);
                Ok::<_, ()>(discontinuous(time))
            },
            |_, candidate, flow, sample| {
                for event in [0.037, 0.074] {
                    assert!(!(sample.start_s < event && event < sample.end_s));
                    assert_ne!(sample.sample_time_s, event);
                }
                assert_eq!(sample.drive, discontinuous(sample.sample_time_s));
                Ok::<_, ()>(stream_exchange(candidate, flow, &sample.drive, 330.0))
            },
        )
        .unwrap();
        let mut manual_owner = owner.clone();
        let mut manual_state = state.clone();
        let mut manual_thermal = old.clone();
        let mut reports = Vec::new();
        for end in [0.037, 0.074, 0.111] {
            let drive =
                discontinuous(manual_thermal.time_s() + 0.5 * (end - manual_thermal.time_s()));
            let step = advance_resolved(
                cx,
                &manual_owner,
                &manual_thermal,
                &manual_state,
                &mut manual_owner.momentum().workspace().unwrap(),
                end,
                &drive.boundaries,
                &drive.injections,
                drive.branch_heads.as_deref(),
                controls(),
                permissive(),
                |_, candidate, flow| Ok::<_, ()>(stream_exchange(candidate, flow, &drive, 330.0)),
            )
            .unwrap();
            reports.extend(step.report.steps);
            manual_owner = step.next;
            manual_state = step.momentum;
            manual_thermal = step.thermal;
        }
        assert_eq!(snapshot(&driven.momentum), snapshot(&manual_state));
        assert_eq!(driven.thermal.cells(), manual_thermal.cells());
        assert_eq!(driven.report.steps.len(), reports.len());
        for event in [0.037, 0.074, 0.111] {
            assert!(
                driven
                    .report
                    .drive_samples
                    .iter()
                    .any(|sample| sample.end_s == event)
            );
        }
        for (actual, expected) in driven.report.steps.iter().zip(reports) {
            assert_eq!(actual.hydraulics.branches, expected.hydraulics.branches);
            assert_eq!(
                actual.external_stream_enthalpy_j,
                expected.external_stream_enthalpy_j
            );
        }
        assert_eq!(samples.borrow().len(), driven.report.work.native_steps);
        assert!(
            driven
                .report
                .steps
                .iter()
                .any(|r| r.hydraulics.branch_head_work_j < 0.0)
        );
        assert!(
            driven
                .report
                .steps
                .iter()
                .any(|r| r.hydraulics.branches[0].displaced_volume_m3 < 0.0)
        );
    });
}

fn smooth(time: f64) -> RegionalCarreauDrive {
    imposed(
        500.0 + 12000.0 * time,
        250.0 * (40.0 * time).sin(),
        1e-8 * time,
    )
}

fn smooth_reference(intervals: usize) -> [f64; 5] {
    let mut y = [0.0, 0.0, 300.0, 300.0, 300.0];
    let dt = 0.125 / intervals as f64;
    let rhs =
        |y: [f64; 5], t: f64| derivative(y, 500.0 + 12000.0 * t + 250.0 * (40.0 * t).sin(), 340.0);
    for index in 0..intervals {
        let time = index as f64 * dt;
        let a = rhs(y, time);
        let b = rhs(
            std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]),
            time + 0.5 * dt,
        );
        let c = rhs(
            std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]),
            time + 0.5 * dt,
        );
        let d = rhs(std::array::from_fn(|i| y[i] + dt * c[i]), time + dt);
        for i in 0..5 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn g1_smooth_midpoint_pressure_head_and_reservoir_injection_match_independent_dynamics_and_work() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let reference = smooth_reference(2048);
        let check = smooth_reference(1024);
        for i in 0..5 {
            near(reference[i], check[i], 1e-8, 0.0);
        }
        let run = |interval| {
            advance_resolved_with_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &[],
                controls(),
                RegionalCarreauResolutionControls {
                    maximum_interval_s: interval,
                    ..permissive()
                },
                |_, time| Ok::<_, ()>(smooth(time)),
                |_, candidate, flow, sample| {
                    // The very small incoming injection shares the same original
                    // inlet enthalpy as the positive reservoir stream in this oracle.
                    assert!(flow.reservoir_outflows[0].value() > 0.0);
                    Ok::<_, ()>(stream_exchange(candidate, flow, &sample.drive, 340.0))
                },
            )
            .unwrap()
        };
        let coarse = run(0.03125);
        let fine = run(0.00390625);
        let error = |result: &ResolvedRegionalCarreau<'_, f64>| {
            [CellId(0), CellId(20), CellId(1)]
                .iter()
                .enumerate()
                .map(|(i, &id)| {
                    (result.thermal.temperature_k(id).unwrap() - reference[i + 2]).abs()
                })
                .fold(0.0, f64::max)
        };
        assert!(
            error(&fine) < 0.3 * error(&coarse),
            "coarse={} fine={}",
            error(&coarse),
            error(&fine)
        );
        assert!(error(&fine) < 0.03, "thermal error {} K", error(&fine));
        let mut velocity = [0.0; 2];
        for (report, sample) in fine.report.steps.iter().zip(&fine.report.drive_samples) {
            assert_eq!(sample.drive, smooth(sample.sample_time_s));
            let dt = sample.end_s - sample.start_s;
            let pressure = sample.drive.boundaries[0].pressure.value();
            let head = sample.drive.branch_heads.as_ref().unwrap()[0].value();
            let injection = sample.drive.injections[0].flow.value();
            let (next, q, heat, kinetic) =
                momentum_reference(velocity, dt, pressure + head, report.thermal.heat);
            near(
                report.hydraulics.branches[0].midpoint_volume_flow_m3_s,
                q,
                1e-16,
                1e-10,
            );
            near(report.hydraulics.viscous_heat_j, heat, 1e-12, 1e-9);
            near(
                report.hydraulics.kinetic_energy_change_j,
                kinetic,
                1e-12,
                1e-9,
            );
            near(
                report.hydraulics.branch_head_work_j,
                head * q * dt,
                1e-12,
                1e-9,
            );
            near(
                report.hydraulics.injection_work_j,
                pressure * injection * dt,
                1e-14,
                1e-9,
            );
            near(
                report.hydraulics.reservoir_work_j,
                pressure * (q - injection) * dt,
                1e-12,
                1e-9,
            );
            assert!(report.combined_energy_residual_j.abs() <= report.combined_energy_tolerance_j);
            velocity = next;
        }
        for ((actual, analytic), continuous) in fine.momentum.branches()[0]
            .velocities_m_s()
            .iter()
            .zip(velocity)
            .zip(reference)
        {
            near(*actual, analytic, 1e-10, 1e-10);
            near(*actual, continuous, 3e-5, 0.0);
        }
    });
}

#[test]
fn g4_invalid_events_late_drive_source_and_free_junction_injection_refuse_atomically() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let initial = snapshot(&state);
        let initial_heat = old.cells().clone();
        let calls = std::cell::Cell::new(0usize);
        let mut workspace = owner.momentum().workspace().unwrap();
        for events in [&[0.0625, 0.0625][..], &[0.1, 0.05], &[f64::NAN]] {
            assert!(matches!(
                advance_resolved_with_drive(
                    cx,
                    &owner,
                    &old,
                    &state,
                    &mut workspace,
                    0.25,
                    events,
                    controls(),
                    permissive(),
                    |_, _| {
                        calls.set(calls.get() + 1);
                        Ok::<_, ()>(imposed(5000.0, 0.0, 0.0))
                    },
                    |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0))
                ),
                Err(RegionalCarreauResolutionError::Invalid(_))
            ));
        }
        assert_eq!(calls.get(), 0);
        assert!(matches!(
            advance_resolved_with_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut workspace,
                0.25,
                &[0.125],
                controls(),
                permissive(),
                |_, time| {
                    if time > 0.125 {
                        Err("pump source outside support")
                    } else {
                        Ok(imposed(5000.0, 0.0, 0.0))
                    }
                },
                |_, _, flow, _| Ok(exchange(flow, 300.0))
            ),
            Err(RegionalCarreauResolutionError::Drive {
                start_s: 0.125,
                source: "pump source outside support",
                ..
            })
        ));
        let free = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.125,
            &[0.0625],
            controls(),
            permissive(),
            |_, time| {
                Ok::<_, ()>(RegionalCarreauDrive {
                    boundaries: vec![pressures(0.0)[1]],
                    injections: vec![FlowInjection {
                        node: 0,
                        flow: VolumetricFlowRate::new(if time < 0.0625 { 0.0 } else { 1e-6 }),
                    }],
                    branch_heads: None,
                    branch_head_laws: None,
                    branch_heat_destinations: None,
                })
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        );
        assert!(matches!(
            free,
            Err(RegionalCarreauResolutionError::Step {
                start_s: 0.0625,
                source: ImplicitTransportError::Resolve(TransientCarreauTransportError::Momentum(
                    TransientCarreauError::InconsistentInitialFlow { node: 0, .. }
                )),
                ..
            })
        ));
        assert_eq!(snapshot(&state), initial);
        assert_eq!(old.cells(), &initial_heat);
        assert_eq!(old.time_s(), 0.0);
        let retry = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.125,
            &[0.0625],
            controls(),
            permissive(),
            |_, time| Ok::<_, ()>(smooth(time)),
            |_, candidate, flow, sample| {
                Ok::<_, ()>(stream_exchange(candidate, flow, &sample.drive, 330.0))
            },
        )
        .unwrap();
        let fresh = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &[0.0625],
            controls(),
            permissive(),
            |_, time| Ok::<_, ()>(smooth(time)),
            |_, candidate, flow, sample| {
                Ok::<_, ()>(stream_exchange(candidate, flow, &sample.drive, 330.0))
            },
        )
        .unwrap();
        assert_eq!(snapshot(&retry.momentum), snapshot(&fresh.momentum));
        assert_eq!(retry.thermal.cells(), fresh.thermal.cells());
        assert_eq!(retry.report.drive_samples, fresh.report.drive_samples);
    });
}

#[test]
fn g4_drive_cancellation_and_window_budget_bound_queries_after_an_accepted_event_prefix() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(&cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], &cx).unwrap();
        let initial = snapshot(&state);
        let initial_heat = old.cells().clone();
        let calls = std::cell::Cell::new(0usize);
        let result = advance_resolved_with_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.25,
            &[0.125],
            controls(),
            RegionalCarreauResolutionControls {
                maximum_native_steps: 3,
                ..permissive()
            },
            |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, ()>(imposed(5000.0, 0.0, 0.0))
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        );
        assert!(matches!(
            result,
            Err(RegionalCarreauResolutionError::WorkLimit { .. })
        ));
        assert_eq!(calls.get(), 3);
        calls.set(0);
        let cancelled = advance_resolved_with_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.25,
            &[0.125],
            controls(),
            permissive(),
            |_, time| {
                calls.set(calls.get() + 1);
                if time > 0.125 {
                    gate.request();
                }
                Ok::<_, ()>(imposed(5000.0, 0.0, 0.0))
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        );
        assert!(matches!(
            cancelled,
            Err(RegionalCarreauResolutionError::Thermal(
                DepositionError::Cancelled
            ))
        ));
        assert_eq!(calls.get(), 4);
        assert_eq!(snapshot(&state), initial);
        assert_eq!(old.cells(), &initial_heat);
        assert_eq!(old.time_s(), 0.0);
    });
}
