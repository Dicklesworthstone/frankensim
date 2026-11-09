//! Original-source thermal drive feedback and independent coupled balances.

use super::*;

const HEATER_W: f64 = 2000.0;
const INLET_K: f64 = 340.0;

fn ports(t: f64, viscosity: f64, head_scale: f64) -> RegionalCarreauDrive {
    RegionalCarreauDrive {
        boundaries: pressures(1000.0 + 10.0 * (t - 300.0).powi(2)).to_vec(),
        injections: Vec::new(),
        branch_heads: Some(vec![Pressure::new(head_scale * viscosity / 0.2)]),
        branch_head_laws: None,
        branch_heat_destinations: None,
    }
}

fn heated(
    candidate: &DepositionSimulation<'_>,
    flow: &TransientCarreauStep,
) -> CarreauTransportExchange<[f64; 3]> {
    CarreauTransportExchange {
        mass_boundaries: exchange(flow, INLET_K).mass_boundaries,
        heat: StepCoefficients {
            sources_w: BTreeMap::from([(CellId(20), HEATER_W)]),
            conductances_w_k: BTreeMap::new(),
            payload: [CellId(0), CellId(20), CellId(1)]
                .map(|id| candidate.temperature_k(id).unwrap()),
        }
        .into(),
    }
}

// Eliminate the original three backward-Euler caloric inventories, then solve
// only the scalar pipe temperature by bisection. Annular momentum uses the
// independent 2x2 physical mass/shear balance in the parent fixture. Neither
// production stepping nor its Newton residual is called by this oracle.
fn coupled_reference(
    old_velocity: [f64; 2],
    old_t: [f64; 3],
    dt: f64,
    head_scale: f64,
) -> Reference {
    let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
    let evaluate = |t: f64| {
        let pressure = 1000.0 + 10.0 * (t - 300.0).powi(2);
        let head = head_scale * mu(t) / 0.2;
        let (velocity, flow, heat, kinetic) =
            momentum_reference(old_velocity, dt, pressure + head, t);
        assert!(flow > 0.0);
        let moved = dt * RHO * flow;
        let upstream = (masses[0] * old_t[0] + moved * INLET_K) / (masses[0] + moved);
        let residual =
            masses[1] * CP * (t - old_t[1]) - moved * CP * (upstream - t) - heat - dt * HEATER_W;
        (residual, velocity, flow, heat, kinetic, upstream, moved)
    };
    let (mut lo, mut hi) = (280.0, 360.0);
    assert!(evaluate(lo).0 < 0.0 && evaluate(hi).0 > 0.0);
    for _ in 0..90 {
        let middle = 0.5 * (lo + hi);
        if evaluate(middle).0 > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    let t = 0.5 * (lo + hi);
    let (_, velocity, mean_flow, heat, kinetic, upstream, moved) = evaluate(t);
    Reference {
        velocity,
        mean_flow,
        heat,
        kinetic,
        temperatures: [
            upstream,
            t,
            (masses[2] * old_t[2] + moved * t) / (masses[2] + moved),
        ],
    }
}

#[test]
fn g1_constant_feedback_matches_prescribed_drive_without_changing_its_query_semantics() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let held = ports(300.0, 0.2, 3000.0);
        let calls = std::cell::Cell::new(0usize);
        let feedback = advance_resolved_with_thermal_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &[],
            controls(),
            permissive(),
            |_, _, candidate, source| {
                calls.set(calls.get() + 1);
                assert_eq!(
                    source.graph().branches()[0].pipe.temperature_k(),
                    candidate.temperature_k(CellId(20)).unwrap()
                );
                Ok::<_, ()>(held.clone())
            },
            |_, candidate, flow, sample| {
                assert_eq!(sample.thermal_time_s, Some(candidate.time_s()));
                Ok::<_, ()>(heated(candidate, flow))
            },
        )
        .unwrap();
        let prescribed = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &[],
            controls(),
            permissive(),
            |_, _| Ok::<_, ()>(held.clone()),
            |_, candidate, flow, sample| {
                assert_eq!(sample.thermal_time_s, None);
                Ok::<_, ()>(heated(candidate, flow))
            },
        )
        .unwrap();
        assert_eq!(snapshot(&feedback.momentum), snapshot(&prescribed.momentum));
        assert_eq!(feedback.thermal.cells(), prescribed.thermal.cells());
        assert_eq!(
            feedback.report.work.native_steps,
            prescribed.report.work.native_steps
        );
        assert_eq!(
            prescribed.report.work.drive_evaluations,
            prescribed.report.work.native_steps
        );
        assert_eq!(feedback.report.work.drive_evaluations, calls.get());
        assert!(calls.get() > feedback.report.work.native_steps);
        assert!(
            calls.get()
                <= feedback.report.work.native_steps
                    * controls().thermal.maximum_endpoint_evaluations
        );
        for (a, b) in feedback
            .report
            .drive_samples
            .iter()
            .zip(&prescribed.report.drive_samples)
        {
            assert_eq!(a.drive, b.drive);
            assert_eq!(a.sample_time_s, b.sample_time_s);
            assert_eq!(a.thermal_time_s, Some(a.end_s));
        }
    });
}

#[test]
fn g1_original_viscosity_head_feedback_matches_independent_nonlinear_heat_and_work() {
    with_cx(|cx| {
        for head_scale in [3000.0, -500.0] {
            let law = Caloric;
            let old = inventory(&law, &[300.0; 3]);
            let owner = model(cx, &old, 1, false);
            let source_law = owner.momentum().graph().branches()[0].pipe.law();
            let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
            let solved = advance_resolved_with_thermal_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.111,
                &[0.037, 0.074],
                controls(),
                permissive(),
                |_, time, candidate, source| {
                    assert!(time < candidate.time_s());
                    assert_ne!(time, 0.037);
                    assert_ne!(time, 0.074);
                    let pipe = source.graph().branches()[0].pipe;
                    let t = candidate.temperature_k(CellId(20)).unwrap();
                    assert_eq!(pipe.temperature_k(), t);
                    assert_eq!(pipe.law(), source_law);
                    assert_eq!(source.pipes()[0].density_kg_m3(), RHO);
                    // The authored controller explicitly observes original
                    // zero-rate viscosity; it is not a fitted flow resistance.
                    let viscosity = pipe
                        .law()
                        .evaluate(t, 0.0)
                        .unwrap()
                        .dynamic_viscosity_pa_s();
                    Ok::<_, ()>(ports(t, viscosity, head_scale))
                },
                |_, candidate, flow, sample| {
                    assert_eq!(sample.thermal_time_s, Some(candidate.time_s()));
                    assert_eq!(sample.end_s, candidate.time_s());
                    Ok::<_, ()>(heated(candidate, flow))
                },
            )
            .unwrap();
            let mut velocity = [0.0; 2];
            let mut temperatures = [300.0; 3];
            for (report, sample) in solved.report.steps.iter().zip(&solved.report.drive_samples) {
                for event in [0.037, 0.074] {
                    assert!(!(sample.start_s < event && event < sample.end_s));
                }
                let dt = sample.end_s - sample.start_s;
                let expected = coupled_reference(velocity, temperatures, dt, head_scale);
                let accepted = report.thermal.heat;
                for (&actual, reference) in accepted.iter().zip(expected.temperatures) {
                    near(actual, reference, 3e-8, 0.0);
                }
                let pressure = 1000.0 + 10.0 * (expected.temperatures[1] - 300.0).powi(2);
                let head = head_scale * mu(expected.temperatures[1]) / 0.2;
                near(
                    sample.drive.boundaries[0].pressure.value(),
                    pressure,
                    1e-6,
                    1e-10,
                );
                near(
                    sample.drive.branch_heads.as_ref().unwrap()[0].value(),
                    head,
                    1e-6,
                    1e-10,
                );
                near(
                    report.hydraulics.branches[0].midpoint_volume_flow_m3_s,
                    expected.mean_flow,
                    1e-14,
                    1e-8,
                );
                near(report.hydraulics.viscous_heat_j, expected.heat, 1e-10, 1e-8);
                near(
                    report.hydraulics.kinetic_energy_change_j,
                    expected.kinetic,
                    1e-10,
                    1e-8,
                );
                near(
                    report.hydraulics.reservoir_work_j,
                    pressure * expected.mean_flow * dt,
                    1e-10,
                    1e-8,
                );
                near(
                    report.hydraulics.branch_head_work_j,
                    head * expected.mean_flow * dt,
                    1e-10,
                    1e-8,
                );
                assert_eq!(
                    report.hydraulics.branch_head_work_j.is_sign_positive(),
                    head_scale.is_sign_positive()
                );
                near(report.other_heat_j, HEATER_W * dt, 1e-10, 1e-12);
                let stream =
                    dt * RHO * expected.mean_flow * CP * (INLET_K - expected.temperatures[2]);
                near(report.external_stream_enthalpy_j, stream, 3e-7, 1e-8);
                let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
                let storage = (0..3)
                    .map(|i| masses[i] * CP * (expected.temperatures[i] - temperatures[i]))
                    .sum::<f64>();
                near(
                    report.thermal.thermal.stored_energy_change_j,
                    storage,
                    3e-6,
                    1e-8,
                );
                // Head work can recover energy; ONLY positive resolved shear
                // loss is deposited in addition to the independently authored heater.
                near(storage, stream + HEATER_W * dt + expected.heat, 1e-8, 1e-10);
                assert!(
                    report.combined_energy_residual_j.abs() <= report.combined_energy_tolerance_j
                );
                velocity = expected.velocity;
                temperatures = expected.temperatures;
            }
            for (actual, expected) in solved.momentum.branches()[0]
                .velocities_m_s()
                .iter()
                .zip(velocity)
            {
                near(*actual, expected, 1e-10, 1e-8);
            }
            for event in [0.037, 0.074] {
                assert!(
                    solved
                        .report
                        .drive_samples
                        .iter()
                        .any(|sample| sample.end_s == event)
                );
            }
            assert!(
                solved.report.drive_samples.last().unwrap().drive != ports(300.0, 0.2, head_scale)
            );
            assert_eq!(
                solved
                    .thermal
                    .material_caloric_identity(MaterialId(1))
                    .unwrap(),
                old.material_caloric_identity(MaterialId(1)).unwrap()
            );
            assert_eq!(state.accepted_time_s(), 0.0);
            assert_eq!(old.time_s(), 0.0);
        }
    });
}

#[test]
fn g4_feedback_refusal_endpoint_budget_and_late_cancellation_preserve_both_histories() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 712,
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
        let before = snapshot(&state);
        let cells = old.cells().clone();
        let mut workspace = owner.momentum().workspace().unwrap();
        let calls = std::cell::Cell::new(0usize);
        let held = ports(300.0, 0.2, 3000.0);
        let refused = advance_resolved_with_thermal_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.25,
            &[0.125],
            controls(),
            permissive(),
            |_, _, candidate, _| {
                calls.set(calls.get() + 1);
                if candidate.time_s() > 0.125 {
                    Err("thermal pump source outside support")
                } else {
                    Ok(held.clone())
                }
            },
            |_, candidate, flow, _| Ok(heated(candidate, flow)),
        );
        match refused {
            Err(RegionalCarreauResolutionError::Step {
                start_s,
                work,
                source:
                    ImplicitTransportError::Resolve(TransientCarreauTransportError::Drive(message)),
                ..
            }) => {
                assert_eq!(start_s, 0.125);
                assert_eq!(work.native_steps, 4);
                assert_eq!(work.drive_evaluations, calls.get());
                assert_eq!(message, "thermal pump source outside support");
            }
            _ => panic!("late feedback refusal must remain a drive-source error"),
        }
        calls.set(0);
        let mut limited = controls();
        limited.thermal.maximum_endpoint_evaluations = 1;
        let budget = advance_resolved_with_thermal_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.125,
            &[],
            limited,
            permissive(),
            |_, _, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, ()>(held.clone())
            },
            |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
        );
        assert!(matches!(
            budget,
            Err(RegionalCarreauResolutionError::Step {
                source: ImplicitTransportError::WorkBudget {
                    performed: 1,
                    maximum: 1
                },
                ..
            })
        ));
        assert_eq!(calls.get(), 1);
        let retry = advance_resolved_with_thermal_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.125,
            &[],
            controls(),
            permissive(),
            |_, _, _, _| Ok::<_, ()>(held.clone()),
            |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
        )
        .unwrap();
        assert_eq!(retry.thermal.time_s(), 0.125);
        calls.set(0);
        let cancelled = advance_resolved_with_thermal_drive(
            &cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.25,
            &[0.125],
            controls(),
            permissive(),
            |_, _, candidate, _| {
                calls.set(calls.get() + 1);
                if candidate.time_s() > 0.125 {
                    gate.request();
                }
                Ok::<_, ()>(held.clone())
            },
            |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
        );
        match cancelled {
            Err(RegionalCarreauResolutionError::Step {
                start_s,
                work,
                source: ImplicitTransportError::Thermal(DepositionError::Cancelled),
                ..
            }) => {
                assert_eq!(start_s, 0.125);
                assert_eq!(work.drive_evaluations, calls.get());
            }
            _ => panic!("late feedback cancellation must preserve the native cancellation"),
        }
        assert_eq!(snapshot(&state), before);
        assert_eq!(old.cells(), &cells);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g1_decimal_interval_cap_reaches_authored_events_in_both_drive_modes() {
    with_cx(|cx| {
        for thermal_feedback in [false, true] {
            let law = Caloric;
            let old = inventory(&law, &[300.0; 3]);
            let owner = model(cx, &old, 1, false);
            let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
            let held = ports(300.0, 0.2, 3000.0);
            let mut policy = permissive();
            policy.maximum_interval_s = 0.037;
            let mut workspace = owner.momentum().workspace().unwrap();
            let solved = if thermal_feedback {
                advance_resolved_with_thermal_drive(
                    cx,
                    &owner,
                    &old,
                    &state,
                    &mut workspace,
                    0.111,
                    &[0.037, 0.074],
                    controls(),
                    policy,
                    |_, _, _, _| Ok::<_, ()>(held.clone()),
                    |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
                )
            } else {
                advance_resolved_with_drive(
                    cx,
                    &owner,
                    &old,
                    &state,
                    &mut workspace,
                    0.111,
                    &[0.037, 0.074],
                    controls(),
                    policy,
                    |_, _| Ok::<_, ()>(held.clone()),
                    |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
                )
            }
            .unwrap();
            assert_eq!(solved.thermal.time_s(), 0.111);
            assert_eq!(solved.momentum.accepted_time_s(), 0.111);
            assert_eq!(solved.report.work.native_steps, 9);
            assert_eq!(solved.report.estimates.len(), 3);
            assert_eq!(solved.report.steps.len(), 6);
            let mut previous_end = 0.0;
            for sample in &solved.report.drive_samples {
                assert_eq!(sample.start_s, previous_end);
                assert!(sample.start_s < sample.sample_time_s);
                assert!(sample.sample_time_s < sample.end_s);
                for event in [0.037, 0.074] {
                    assert!(!(sample.start_s < event && event < sample.end_s));
                }
                previous_end = sample.end_s;
            }
            for endpoint in [0.037, 0.074, 0.111] {
                assert!(
                    solved
                        .report
                        .drive_samples
                        .iter()
                        .any(|sample| sample.end_s == endpoint)
                );
            }
            assert_eq!(previous_end, 0.111);
            assert_eq!(state.accepted_time_s(), 0.0);
            assert_eq!(old.time_s(), 0.0);
        }
    });
}
