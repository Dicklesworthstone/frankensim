//! Native time refinement of radial history and reversing thermal transport.

mod drive;
mod equipment;
mod feedback;
mod head_laws;

use super::super::adaptive::*;
use super::*;

fn resolution() -> RegionalCarreauResolutionControls {
    RegionalCarreauResolutionControls {
        maximum_interval_s: 0.125,
        minimum_interval_s: 1e-6,
        maximum_native_steps: 4096,
        maximum_accepted_steps: 2048,
        maximum_total_branch_solves: 200_000,
        maximum_total_head_evaluations: 400_000,
        absolute_cell_error_j: 0.002,
        relative_cell_error: 1e-4,
        absolute_velocity_error_m_s: 1e-4,
        relative_velocity_error: 1e-4,
        absolute_flow_error_m3_s: 1e-8,
        relative_flow_error: 1e-4,
        absolute_transfer_error_j: 0.002,
        relative_transfer_error: 1e-4,
        maximum_reversal_step_s: 0.002,
        maximum_cancelled_volume_m3: 1e-9,
    }
}

fn permissive() -> RegionalCarreauResolutionControls {
    RegionalCarreauResolutionControls {
        absolute_cell_error_j: 1e6,
        relative_cell_error: 0.0,
        absolute_velocity_error_m_s: 1e6,
        relative_velocity_error: 0.0,
        absolute_flow_error_m3_s: 1e6,
        relative_flow_error: 0.0,
        absolute_transfer_error_j: 1e6,
        relative_transfer_error: 0.0,
        maximum_reversal_step_s: 1.0,
        maximum_cancelled_volume_m3: 1.0,
        ..resolution()
    }
}

/// Declare the original native inventory allowance for convergence fixtures;
/// adaptive half-steps must not silently renew a smaller caloric step budget.
fn inventory_with_substeps<'m>(
    law: &'m Caloric,
    temperatures: &[f64],
    maximum: usize,
) -> DepositionSimulation<'m> {
    let base = inventory(law, temperatures);
    DepositionSimulation::new(
        base.cells()
            .values()
            .cloned()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect(),
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_substeps: maximum,
            ..base.budget()
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

fn initial_profile(
    owner: &TransientRegionalCarreauTransport,
    cx: &Cx<'_>,
) -> TransientCarreauState {
    owner
        .momentum()
        .state_from_velocities(vec![vec![1.25, 0.625]], &[None; 2], 0.0, cx)
        .unwrap()
}

/// Independent two-annulus ODE and three mixed enthalpy inventories. The
/// constitutive value is evaluated analytically by the test's Arrhenius law;
/// no native momentum or thermal stepping is called by this reference.
fn derivative(y: [f64; 5], dp: f64, inlet: f64) -> [f64; 5] {
    let area = core::f64::consts::PI * R * R;
    let areas = [0.25 * area, 0.75 * area];
    let k = 2.0 * core::f64::consts::PI * mu(y[3]) * L;
    let shear = k * (y[0] - y[1]);
    let flow = areas[0] * y[0] + areas[1] * y[1];
    let mass_rate = RHO * flow.abs();
    let heat = k * (y[0] - y[1]).powi(2) + 4.0 * k * y[1] * y[1];
    let mut result = [
        (dp * areas[0] - shear) / (RHO * L * areas[0]),
        (dp * areas[1] + shear - 4.0 * k * y[1]) / (RHO * L * areas[1]),
        0.0,
        0.0,
        0.0,
    ];
    let (up, down) = if flow >= 0.0 { (2, 4) } else { (4, 2) };
    let mass = [0.0, 0.0, 0.04, RHO * area * L, 0.03];
    result[up] = mass_rate * (inlet - y[up]) / mass[up];
    result[3] = (mass_rate * CP * (y[up] - y[3]) + heat) / (mass[3] * CP);
    result[down] = mass_rate * (y[3] - y[down]) / mass[down];
    result
}

fn reference_rk4(
    mut y: [f64; 5],
    duration: f64,
    intervals: usize,
    dp: f64,
    inlet: f64,
) -> [f64; 5] {
    let dt = duration / intervals as f64;
    for _ in 0..intervals {
        let a = derivative(y, dp, inlet);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]), dp, inlet);
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]), dp, inlet);
        let d = derivative(std::array::from_fn(|i| y[i] + dt * c[i]), dp, inlet);
        for i in 0..5 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn g1_reversing_native_trajectory_converges_to_independent_radial_and_caloric_odes() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory_with_substeps(&law, &[330.0, 290.0, 300.0], 2048);
        let owner = model(cx, &old, 1, false);
        let state = initial_profile(&owner, cx);
        let initial = [1.25, 0.625, 330.0, 300.0, 290.0];
        let reference = reference_rk4(initial, 0.125, 4096, -5000.0, 340.0);
        let check = reference_rk4(initial, 0.125, 2048, -5000.0, 340.0);
        for i in 0..5 {
            near(reference[i], check[i], 2e-7, 0.0);
        }
        let run = |settings| {
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &pressures(-5000.0),
                &[],
                None,
                controls(),
                settings,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 340.0)),
            )
            .unwrap()
        };
        let coarse = run(RegionalCarreauResolutionControls {
            maximum_interval_s: 0.03125,
            ..permissive()
        });
        let fine = run(resolution());
        let error = |result: &ResolvedRegionalCarreau<'_, ()>| {
            [CellId(0), CellId(20), CellId(1)]
                .iter()
                .enumerate()
                .map(|(i, &id)| {
                    (result.thermal.temperature_k(id).unwrap() - reference[i + 2]).abs()
                })
                .fold(0.0, f64::max)
        };
        assert!(
            error(&fine) < 0.5 * error(&coarse),
            "coarse={} fine={}",
            error(&coarse),
            error(&fine)
        );
        assert!(error(&fine) < 0.05, "thermal error {} K", error(&fine));
        for (actual, expected) in fine.momentum.branches()[0]
            .velocities_m_s()
            .iter()
            .zip(reference)
        {
            near(*actual, expected, 5e-4, 0.0);
        }
        assert!(fine.report.work.rejected_intervals > 0);
        assert!(
            fine.report
                .steps
                .iter()
                .any(|report| report.hydraulics.branches[0].displaced_volume_m3 > 0.0)
        );
        assert!(
            fine.report
                .steps
                .iter()
                .any(|report| report.hydraulics.branches[0].displaced_volume_m3 < 0.0)
        );
        let mut heat = 0.0;
        let mut work = 0.0;
        let mut stream = 0.0;
        let mut other = 0.0;
        let mut kinetic = 0.0;
        let mut stored = 0.0;
        let mut tolerance = 0.0;
        let mut accepted_work = 0;
        for report in &fine.report.steps {
            let h = &report.hydraulics;
            heat += h.viscous_heat_j;
            work += h.reservoir_work_j + h.injection_work_j + h.branch_head_work_j;
            kinetic += h.kinetic_energy_change_j;
            stored += report.thermal.thermal.stored_energy_change_j;
            stream += report.external_stream_enthalpy_j;
            other += report.other_heat_j;
            tolerance += report.combined_energy_tolerance_j;
            accepted_work += report.total_branch_solves;
            assert!(h.viscous_heat_j >= 0.0);
            assert!(report.combined_energy_residual_j.abs() <= report.combined_energy_tolerance_j);
        }
        assert!((stored + kinetic - work - stream - other).abs() <= tolerance);
        near(stored, heat + stream + other, tolerance, 0.0);
        assert!(fine.report.work.completed_branch_solves > accepted_work);
        assert_eq!(fine.thermal.time_s(), 0.125);
        assert_eq!(fine.momentum.accepted_time_s(), 0.125);
        assert_eq!(
            fine.next.momentum().graph().branches()[0].pipe.law(),
            owner.momentum().graph().branches()[0].pipe.law()
        );
    });
}

#[test]
fn g3_counterflow_volume_refines_uniform_temperature_even_when_state_errors_are_disabled() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = initial_profile(&owner, cx);
        let run = |settings| {
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &pressures(-5000.0),
                &[],
                None,
                controls(),
                settings,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
            )
            .unwrap()
        };
        let unresolved = run(permissive());
        assert_eq!(unresolved.report.steps.len(), 2);
        assert!(unresolved.report.estimates[0].maximum_cancelled_volume_m3 > 1e-8);
        let settings = RegionalCarreauResolutionControls {
            maximum_cancelled_volume_m3: 1e-10,
            maximum_reversal_step_s: 0.001,
            ..permissive()
        };
        let resolved = run(settings);
        assert!(resolved.report.steps.len() > 2);
        assert!(resolved.report.work.rejected_intervals > 0);
        for estimate in &resolved.report.estimates {
            assert!(estimate.maximum_cancelled_volume_m3 <= settings.maximum_cancelled_volume_m3);
            assert!(estimate.maximum_reversal_step_s <= settings.maximum_reversal_step_s);
        }
        assert!(
            resolved
                .report
                .steps
                .iter()
                .any(|r| r.hydraulics.branches[0].displaced_volume_m3 > 0.0)
        );
        assert!(
            resolved
                .report
                .steps
                .iter()
                .any(|r| r.hydraulics.branches[0].displaced_volume_m3 < 0.0)
        );
        let limited = RegionalCarreauResolutionControls {
            minimum_interval_s: 0.125,
            ..settings
        };
        assert!(matches!(
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &pressures(-5000.0),
                &[],
                None,
                controls(),
                limited,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0))
            ),
            Err(RegionalCarreauResolutionError::Underresolved { .. })
        ));
    });
}

#[test]
fn g3_annular_velocity_resolution_cannot_be_replaced_by_loose_mean_flow_or_heat_controls() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        for velocities in [[0.8, 0.2], [0.2, 0.4]] {
            let state = owner
                .momentum()
                .state_from_velocities(vec![velocities.to_vec()], &[None; 2], 0.0, cx)
                .unwrap();
            let settings = RegionalCarreauResolutionControls {
                absolute_velocity_error_m_s: 1e-5,
                ..permissive()
            };
            let result = advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &pressures(0.0),
                &[],
                None,
                controls(),
                settings,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
            )
            .unwrap();
            assert!(result.report.steps.len() > 2);
            assert!(result.report.work.rejected_intervals > 0);
            assert!(
                result
                    .report
                    .estimates
                    .iter()
                    .all(|e| e.maximum_velocity_difference_m_s <= 1e-5 && e.error_ratio <= 1.0)
            );
            // Initial net flow is identical, but every original velocity and
            // its subsequent shear relaxation remains a resolved coordinate.
            near(
                result.report.steps[0].hydraulics.branches[0].initial_volume_flow_m3_s,
                core::f64::consts::PI * R * R * 0.35,
                1e-18,
                0.0,
            );
        }
    });
}

#[test]
fn g4_window_budget_and_late_source_refusals_preserve_inputs_and_reused_workspace() {
    with_cx(|cx| {
        let law = Caloric;
        let old = inventory(&law, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let before = snapshot(&state);
        let caloric = old.cells().clone();
        let mut workspace = owner.momentum().workspace().unwrap();
        let limited = RegionalCarreauResolutionControls {
            maximum_accepted_steps: 2,
            ..permissive()
        };
        assert!(matches!(
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut workspace,
                0.25,
                &pressures(5000.0),
                &[],
                None,
                controls(),
                limited,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0))
            ),
            Err(RegionalCarreauResolutionError::WorkLimit {
                work: RegionalCarreauResolutionWork {
                    native_steps: 3,
                    ..
                },
                ..
            })
        ));
        assert!(matches!(
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut workspace,
                0.25,
                &pressures(5000.0),
                &[],
                None,
                controls(),
                permissive(),
                |_, candidate, flow| {
                    if candidate.time_s() > 0.125 {
                        Err("source support")
                    } else {
                        Ok(exchange(flow, 300.0))
                    }
                }
            ),
            Err(RegionalCarreauResolutionError::Step { start_s: 0.125, .. })
        ));
        let limited = RegionalCarreauResolutionControls {
            maximum_total_branch_solves: 1,
            ..permissive()
        };
        assert!(
            advance_resolved(
                cx,
                &owner,
                &old,
                &state,
                &mut workspace,
                0.125,
                &pressures(5000.0),
                &[],
                None,
                controls(),
                limited,
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0))
            )
            .is_err()
        );
        assert_eq!(snapshot(&state), before);
        assert_eq!(old.cells(), &caloric);
        assert_eq!(old.time_s(), 0.0);
        let retry = advance_resolved(
            cx,
            &owner,
            &old,
            &state,
            &mut workspace,
            0.125,
            &pressures(5000.0),
            &[],
            Some(&[Pressure::new(100.0)]),
            controls(),
            permissive(),
            |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
        )
        .unwrap();
        let fresh = advance_resolved(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &pressures(5000.0),
            &[],
            Some(&[Pressure::new(100.0)]),
            controls(),
            permissive(),
            |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
        )
        .unwrap();
        assert_eq!(snapshot(&retry.momentum), snapshot(&fresh.momentum));
        assert_eq!(retry.thermal.cells(), fresh.thermal.cells());
        assert!(
            retry
                .report
                .steps
                .iter()
                .all(|r| r.hydraulics.branch_head_work_j > 0.0)
        );
    });
}

#[test]
fn g4_cancellation_after_an_accepted_local_prefix_publishes_neither_history() {
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
        let before = snapshot(&state);
        let caloric = old.cells().clone();
        let result = advance_resolved(
            &cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.25,
            &pressures(5000.0),
            &[],
            None,
            controls(),
            permissive(),
            |_, candidate, flow| {
                if candidate.time_s() > 0.125 {
                    gate.request();
                }
                Ok::<_, ()>(exchange(flow, 300.0))
            },
        );
        assert!(matches!(
            result,
            Err(RegionalCarreauResolutionError::Step { start_s: 0.125, .. })
        ));
        assert_eq!(snapshot(&state), before);
        assert_eq!(old.cells(), &caloric);
        assert_eq!(old.time_s(), 0.0);
    });
}
