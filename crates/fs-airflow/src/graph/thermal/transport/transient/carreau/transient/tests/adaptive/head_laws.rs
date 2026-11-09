//! Original nonlinear equipment, radial history and regional caloric feedback.

use super::*;
use crate::graph::rheology::transient::head::{
    FlowHeadError, FlowHeadLaw, FlowHeadResponse, FlowHeadSource,
};
use std::sync::{
    Arc,
    atomic::{AtomicU64, Ordering},
};

const HEATER: f64 = 2000.0;
const INLET: f64 = 340.0;

#[derive(Debug)]
struct Equipment {
    shutoff: f64,
    queries: AtomicU64,
}
impl Equipment {
    fn pressure(&self, q: f64, viscosity: f64) -> f64 {
        self.shutoff * viscosity / 0.2 - 2e6 * q - 4e15 * q.powi(3)
    }
}

/// A prepared view retains the exact authored equipment object and observes
/// original local viscosity; it never replaces the source with a fitted curve.
#[derive(Debug)]
struct PreparedEquipment {
    original: Arc<Equipment>,
    viscosity: f64,
}
impl FlowHeadLaw for PreparedEquipment {
    fn flow_interval_m3_s(&self) -> [f64; 2] {
        [-0.001, 0.001]
    }
    fn reference_density_kg_m3(&self) -> Option<f64> {
        Some(RHO)
    }
    fn evaluate(
        &self,
        flow: VolumetricFlowRate,
    ) -> Result<FlowHeadResponse, Box<dyn std::error::Error + Send + Sync>> {
        self.original.queries.fetch_add(1, Ordering::Relaxed);
        let q = flow.value();
        if !(-0.001..=0.001).contains(&q) {
            return Err(std::io::Error::other("original equipment support").into());
        }
        Ok(FlowHeadResponse {
            pressure: Pressure::new(self.original.pressure(q, self.viscosity)),
            pressure_slope_pa_s_m3: -2e6 - 12e15 * q * q,
            power: None,
        })
    }
}
fn drive(original: &Arc<Equipment>, source: &TransientCarreauGraph) -> RegionalCarreauDrive {
    let pipe = source.graph().branches()[0].pipe;
    let viscosity = pipe
        .law()
        .evaluate(pipe.temperature_k(), 0.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    RegionalCarreauDrive {
        boundaries: pressures(1000.0).to_vec(),
        injections: Vec::new(),
        branch_heads: None,
        branch_head_laws: Some(vec![Some(FlowHeadSource::new(Arc::new(
            PreparedEquipment {
                original: original.clone(),
                viscosity,
            },
        )))]),
        branch_heat_destinations: None,
    }
}
fn heated(
    candidate: &DepositionSimulation<'_>,
    flow: &TransientCarreauStep,
) -> CarreauTransportExchange<[f64; 3]> {
    CarreauTransportExchange {
        mass_boundaries: exchange(flow, INLET).mass_boundaries,
        heat: StepCoefficients {
            sources_w: BTreeMap::from([(CellId(20), HEATER)]),
            conductances_w_k: BTreeMap::new(),
            payload: [CellId(0), CellId(20), CellId(1)]
                .map(|id| candidate.temperature_k(id).unwrap()),
        }
        .into(),
    }
}

// Independent nested scalar balances: the inner equipment root uses the
// analytically assembled 2x2 annular matrix; the outer root eliminates the
// three original backward-Euler mass inventories. No production solve is used.
fn reference_step(
    old_velocity: [f64; 2],
    old_t: [f64; 3],
    dt: f64,
    equipment: &Equipment,
) -> Reference {
    let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
    let evaluate = |t: f64| {
        let offset = momentum_reference(old_velocity, dt, 0.0, t).1;
        let conductance = momentum_reference([0.0; 2], dt, 1.0, t).1;
        let residual = |q: f64| q - offset - conductance * (1000.0 + equipment.pressure(q, mu(t)));
        let (mut lo, mut hi) = (-0.001, 0.001);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..90 {
            let q = f64::midpoint(lo, hi);
            if residual(q) > 0.0 {
                hi = q;
            } else {
                lo = q;
            }
        }
        let pressure = 1000.0 + equipment.pressure(f64::midpoint(lo, hi), mu(t));
        let (velocity, flow, heat, kinetic) = momentum_reference(old_velocity, dt, pressure, t);
        assert!(flow > 0.0);
        let moved = dt * RHO * flow;
        let upstream = (masses[0] * old_t[0] + moved * INLET) / (masses[0] + moved);
        let residual =
            masses[1] * CP * (t - old_t[1]) - moved * CP * (upstream - t) - heat - dt * HEATER;
        (residual, velocity, flow, heat, kinetic, upstream, moved)
    };
    let (mut lo, mut hi) = (280.0, 360.0);
    assert!(evaluate(lo).0 < 0.0 && evaluate(hi).0 > 0.0);
    for _ in 0..90 {
        let t = f64::midpoint(lo, hi);
        if evaluate(t).0 > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = f64::midpoint(lo, hi);
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
fn g1_nonlinear_equipment_and_original_viscosity_close_the_same_regional_endpoint() {
    with_cx(|cx| {
        for shutoff in [2000.0, -400.0] {
            let equipment = Arc::new(Equipment {
                shutoff,
                queries: AtomicU64::new(0),
            });
            let caloric = Caloric;
            let old = inventory(&caloric, &[300.0; 3]);
            let owner = model(cx, &old, 1, false);
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
                |_, _, candidate, source| {
                    assert_eq!(
                        source.graph().branches()[0].pipe.temperature_k(),
                        candidate.temperature_k(CellId(20)).unwrap()
                    );
                    Ok::<_, ()>(drive(&equipment, source))
                },
                |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
            )
            .unwrap();
            assert_eq!(
                solved.report.work.completed_head_evaluations,
                equipment.queries.load(Ordering::Relaxed)
            );
            assert!(
                solved.report.work.completed_head_evaluations
                    > solved
                        .report
                        .steps
                        .iter()
                        .map(|step| step.total_head_evaluations)
                        .sum()
            );
            let mut velocity = [0.0; 2];
            let mut temperatures = [300.0; 3];
            for (report, sample) in solved.report.steps.iter().zip(&solved.report.drive_samples) {
                let dt = sample.end_s - sample.start_s;
                let expected = reference_step(velocity, temperatures, dt, &equipment);
                for (&actual, expected) in report.thermal.heat.iter().zip(expected.temperatures) {
                    near(actual, expected, 3e-8, 0.0);
                }
                let accepted = report.hydraulics.flow_heads[0].as_ref().unwrap();
                assert_eq!(
                    &accepted.source,
                    sample.drive.branch_head_laws.as_ref().unwrap()[0]
                        .as_ref()
                        .unwrap()
                );
                let head = equipment.pressure(expected.mean_flow, mu(expected.temperatures[1]));
                near(accepted.response.pressure.value(), head, 1e-6, 1e-9);
                near(
                    accepted.midpoint_flow.value(),
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
                    report.hydraulics.branch_head_work_j,
                    head * expected.mean_flow * dt,
                    1e-10,
                    1e-8,
                );
                assert_eq!(
                    accepted.work_j.is_sign_positive(),
                    shutoff.is_sign_positive()
                );
                near(report.other_heat_j, HEATER * dt, 1e-10, 1e-12);
                let stream =
                    dt * RHO * expected.mean_flow * CP * (INLET - expected.temperatures[2]);
                let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
                let storage = (0..3)
                    .map(|i| masses[i] * CP * (expected.temperatures[i] - temperatures[i]))
                    .sum::<f64>();
                near(report.external_stream_enthalpy_j, stream, 3e-7, 1e-8);
                near(
                    report.thermal.thermal.stored_energy_change_j,
                    storage,
                    3e-6,
                    1e-8,
                );
                near(storage, stream + HEATER * dt + expected.heat, 1e-8, 1e-10);
                assert!(
                    report.combined_energy_residual_j.abs() <= report.combined_energy_tolerance_j
                );
                velocity = expected.velocity;
                temperatures = expected.temperatures;
            }
            for (&actual, expected) in solved.momentum.branches()[0]
                .velocities_m_s()
                .iter()
                .zip(velocity)
            {
                near(actual, expected, 1e-10, 1e-8);
            }
        }
    });
}

#[test]
fn g4_head_queries_are_bounded_across_discarded_trials_and_late_failure_is_atomic() {
    with_cx(|cx| {
        let equipment = Arc::new(Equipment {
            shutoff: 2000.0,
            queries: AtomicU64::new(0),
        });
        let caloric = Caloric;
        let old = inventory(&caloric, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let initial = snapshot(&state);
        let initial_cells = old.cells().clone();
        let mut workspace = owner.momentum().workspace().unwrap();
        let mut run = |settings| {
            advance_resolved_with_thermal_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut workspace,
                0.125,
                &[0.0625],
                controls(),
                settings,
                |_, _, _, source| Ok::<_, ()>(drive(&equipment, source)),
                |_, candidate, flow, _| Ok::<_, ()>(heated(candidate, flow)),
            )
        };
        let admitted = run(permissive()).unwrap();
        let queries = admitted.report.work.completed_head_evaluations;
        assert!(queries > 1);
        equipment.queries.store(0, Ordering::Relaxed);
        let error = run(RegionalCarreauResolutionControls {
            maximum_total_head_evaluations: queries - 1,
            ..permissive()
        })
        .err()
        .unwrap();
        assert!(matches!(error, RegionalCarreauResolutionError::Step {
            start_s, source: ImplicitTransportError::Resolve(TransientCarreauTransportError::Momentum(TransientCarreauError::Head { source: FlowHeadError::SourceBudget, .. })), ..
        } if start_s >= 0.0625));
        assert_eq!(equipment.queries.load(Ordering::Relaxed), queries - 1);
        assert_eq!(snapshot(&state), initial);
        assert_eq!(old.cells(), &initial_cells);
        assert_eq!(old.time_s(), 0.0);
        equipment.queries.store(0, Ordering::Relaxed);
        let retry = run(permissive()).unwrap();
        assert_eq!(snapshot(&retry.momentum), snapshot(&admitted.momentum));
        assert_eq!(retry.thermal.cells(), admitted.thermal.cells());
    });
}

#[test]
fn g1_an_empty_equipment_roster_does_not_spend_an_equipment_budget() {
    with_cx(|cx| {
        let caloric = Caloric;
        let old = inventory(&caloric, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let solved = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.02,
            &[],
            controls(),
            RegionalCarreauResolutionControls {
                maximum_total_head_evaluations: 0,
                ..permissive()
            },
            |_, _| {
                Ok::<_, ()>(RegionalCarreauDrive {
                    boundaries: pressures(1000.0).to_vec(),
                    injections: vec![],
                    branch_heads: None,
                    branch_head_laws: Some(vec![None]),
                    branch_heat_destinations: None,
                })
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        )
        .unwrap();
        assert_eq!(solved.report.work.completed_head_evaluations, 0);
        assert!(
            solved
                .report
                .steps
                .iter()
                .all(|step| step.hydraulics.flow_heads == [None])
        );
    });
}

#[test]
fn g1_opposed_time_varying_heads_require_work_resolution_at_constant_flow_and_temperature() {
    with_cx(|cx| {
        let caloric = Caloric;
        let old = inventory(&caloric, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        // Exact discrete steady annuli for dp=1000 Pa, mu=.2 Pa s.
        let state = owner
            .momentum()
            .state_from_velocities(vec![vec![0.25, 0.125]], &[None; 2], 0.0, cx)
            .unwrap();
        let flow = core::f64::consts::PI * R * R * (0.25 * 0.25 + 0.75 * 0.125);
        let amplitude = 1e6;
        let target: f64 = 0.1;
        let run = |settings| {
            advance_resolved_with_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                target,
                &[],
                controls(),
                settings,
                |_, time| {
                    let original = Arc::new(Equipment {
                        shutoff: amplitude * time * time,
                        queries: AtomicU64::new(0),
                    });
                    let held = -original.pressure(flow, 0.2);
                    let mut sample = drive(&original, owner.momentum());
                    sample.branch_heads = Some(vec![Pressure::new(held)]);
                    Ok::<_, ()>(sample)
                },
                |_, _, response, _| {
                    // Explicit cooling cancels only physical pipe dissipation.
                    // The two independent actuator works remain separate channels.
                    let mut thermal = exchange(response, 300.0);
                    thermal
                        .heat
                        .bulk
                        .sources_w
                        .insert(CellId(20), -response.viscous_heat_j / response.duration_s);
                    Ok::<_, ()>(thermal)
                },
            )
            .unwrap()
        };
        let coarse = run(permissive());
        let fine = run(RegionalCarreauResolutionControls {
            absolute_transfer_error_j: 1e-5,
            relative_transfer_error: 0.0,
            ..permissive()
        });
        assert_eq!(coarse.report.steps.len(), 2);
        assert!(fine.report.work.rejected_intervals > 0);
        assert!(fine.report.steps.len() > coarse.report.steps.len());
        let work = |result: &ResolvedRegionalCarreau<'_, ()>| {
            result
                .report
                .steps
                .iter()
                .map(|step| step.hydraulics.flow_heads[0].as_ref().unwrap().work_j)
                .sum::<f64>()
        };
        let exact =
            flow * (amplitude * target.powi(3) / 3.0 - (2e6 * flow + 4e15 * flow.powi(3)) * target);
        assert!((work(&fine) - exact).abs() < (work(&coarse) - exact).abs() / 20.0);
        for result in [&coarse, &fine] {
            for (&actual, expected) in result.momentum.branches()[0]
                .velocities_m_s()
                .iter()
                .zip([0.25, 0.125])
            {
                near(actual, expected, 1e-11, 1e-10);
            }
            for id in [CellId(0), CellId(20), CellId(1)] {
                near(result.thermal.temperature_k(id).unwrap(), 300.0, 1e-10, 0.0);
            }
            assert!(
                result
                    .report
                    .steps
                    .iter()
                    .all(|step| step.hydraulics.branch_head_work_j.abs() < 1e-11)
            );
        }
    });
}

#[test]
fn g4_original_native_thermal_substep_limit_applies_to_the_whole_requested_window() {
    with_cx(|cx| {
        let caloric = Caloric;
        let base = inventory(&caloric, &[300.0; 3]);
        let old = DepositionSimulation::new(
            base.cells()
                .values()
                .cloned()
                .map(|cell| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                })
                .collect(),
            &[(MaterialId(1), &caloric as &dyn ThermalLaw)],
            DepositionBudget {
                max_substeps: 3,
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
        .0;
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let initial = snapshot(&state);
        let cells = old.cells().clone();
        let result = advance_resolved_with_drive(
            cx,
            &owner,
            &old,
            &state,
            &mut owner.momentum().workspace().unwrap(),
            0.125,
            &[0.0625],
            controls(),
            permissive(),
            |_, _| {
                Ok::<_, ()>(RegionalCarreauDrive {
                    boundaries: pressures(1000.0).to_vec(),
                    injections: vec![],
                    branch_heads: None,
                    branch_head_laws: None,
                    branch_heat_destinations: None,
                })
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        );
        assert!(matches!(
            result,
            Err(RegionalCarreauResolutionError::WorkLimit {
                reason: "retained native step limit",
                work: RegionalCarreauResolutionWork {
                    native_steps: 3,
                    ..
                }
            })
        ));
        assert_eq!(snapshot(&state), initial);
        assert_eq!(old.cells(), &cells);
        assert_eq!(old.time_s(), 0.0);
    });
}
