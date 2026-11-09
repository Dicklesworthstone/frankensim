//! Independent caloric/input-power balances with original retained momentum.

use super::super::super::equipment::EquipmentHeatDestinations;
use super::*;
use crate::graph::rheology::transient::head::power::{
    FlowHeadPowerReceipt, PumpLossRouting, PumpPowerState,
};
use crate::graph::rheology::transient::head::{FlowHeadLaw, FlowHeadResponse, FlowHeadSource};
use std::sync::{
    Arc,
    atomic::{AtomicU64, Ordering},
};

const INLET: f64 = 330.0;

/// Explicit synthetic equipment data; input power is independently authored,
/// never inferred from hydraulic head or an assumed efficiency.
#[derive(Debug)]
struct Calibration {
    temperature_gain: f64,
    time_quadratic: f64,
    queries: AtomicU64,
}
fn calibration(temperature_gain: f64, time_quadratic: f64) -> Arc<Calibration> {
    Arc::new(Calibration {
        temperature_gain,
        time_quadratic,
        queries: AtomicU64::new(0),
    })
}
fn head(q: f64) -> f64 {
    800.0 - 2e5 * q - 4e10 * q.powi(3)
}
fn shaft(data: &Calibration, q: f64, temperature: f64, time: f64) -> f64 {
    200.0
        + 1e5 * q
        + data.temperature_gain * (temperature - 300.0)
        + data.time_quadratic * time * time
}
#[derive(Debug)]
struct PreparedEquipment {
    original: Arc<Calibration>,
    temperature: f64,
    time: f64,
}
impl FlowHeadLaw for PreparedEquipment {
    fn flow_interval_m3_s(&self) -> [f64; 2] {
        [0.0, 0.001]
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
        if !(0.0..=0.001).contains(&q) {
            return Err(std::io::Error::other("original equipment flow span").into());
        }
        let shaft = shaft(&self.original, q, self.temperature, self.time);
        let loss = shaft - head(q) * q;
        Ok(FlowHeadResponse {
            pressure: Pressure::new(head(q)),
            pressure_slope_pa_s_m3: -2e5 - 12e10 * q * q,
            power: Some(FlowHeadPowerReceipt {
                state: PumpPowerState {
                    shaft_input_w: shaft,
                    electrical_input_w: shaft + 20.0,
                    pump_loss_w: loss,
                    drive_loss_w: 20.0,
                    reservoir_heat_w: 0.5 * loss + 5.0,
                    surroundings_heat_w: 0.5 * loss + 15.0,
                },
                routing: PumpLossRouting {
                    pump_to_reservoir: 0.5,
                    drive_to_reservoir: 0.25,
                },
            }),
        })
    }
}
fn drive(
    data: &Arc<Calibration>,
    temperature: f64,
    time: f64,
    pressure: f64,
    pump: CellId,
    motor: CellId,
) -> RegionalCarreauDrive {
    RegionalCarreauDrive {
        boundaries: pressures(pressure).to_vec(),
        injections: vec![],
        branch_heads: None,
        branch_head_laws: Some(vec![Some(FlowHeadSource::new(Arc::new(
            PreparedEquipment {
                original: data.clone(),
                temperature,
                time,
            },
        )))]),
        branch_heat_destinations: Some(vec![EquipmentHeatDestinations {
            pump: Some(pump),
            drive: Some(motor),
        }]),
    }
}
fn extra_heat(
    candidate: &DepositionSimulation<'_>,
    flow: &TransientCarreauStep,
) -> CarreauTransportExchange<[f64; 3]> {
    CarreauTransportExchange {
        mass_boundaries: exchange(flow, INLET).mass_boundaries,
        heat: StepCoefficients {
            sources_w: BTreeMap::from([(CellId(20), 3.0), (CellId(0), 1.0)]),
            conductances_w_k: BTreeMap::new(),
            payload: [CellId(0), CellId(20), CellId(1)]
                .map(|id| candidate.temperature_k(id).unwrap()),
        }
        .into(),
    }
}

#[test]
fn g1_shutoff_loss_heat_uses_distinct_or_shared_recipients_and_preserves_other_heat() {
    with_cx(|cx| {
        for motor in [CellId(0), CellId(20)] {
            let original = calibration(0.0, 0.0);
            let caloric = Caloric;
            let old = inventory(&caloric, &[300.0; 3]);
            let owner = model(cx, &old, 1, false);
            let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
            let forcing = drive(&original, 300.0, 0.125, -800.0, CellId(20), motor);
            let dt = 0.25;
            let solved = owner
                .advance_with_drive(
                    cx,
                    &old,
                    &state,
                    &mut owner.momentum().workspace().unwrap(),
                    dt,
                    &forcing,
                    controls(),
                    &[0.0; 3],
                    |_, candidate, flow| Ok::<_, ()>(extra_heat(candidate, flow)),
                )
                .unwrap();
            let report = &solved.report;
            let energy = report.equipment.branches[0].unwrap();
            assert_eq!(
                report.total_head_evaluations,
                original.queries.load(Ordering::Relaxed)
            );
            assert_eq!(
                report.hydraulics.flow_heads[0].as_ref().unwrap().source,
                *forcing.branch_head_laws.as_ref().unwrap()[0]
                    .as_ref()
                    .unwrap()
            );
            assert_eq!(report.hydraulics.viscous_heat_j, 0.0);
            assert_eq!(report.hydraulics.branch_head_work_j, 0.0);
            assert_eq!(energy.hydraulic_work_j, 0.0);
            near(energy.shaft_input_j, 200.0 * dt, 1e-12, 0.0);
            near(energy.electrical_input_j, 220.0 * dt, 1e-12, 0.0);
            near(energy.pump_heat_j, 100.0 * dt, 1e-12, 0.0);
            near(energy.drive_heat_j, 5.0 * dt, 1e-12, 0.0);
            near(energy.surroundings_heat_j, 115.0 * dt, 1e-12, 0.0);
            near(report.other_heat_j, 4.0 * dt, 1e-12, 0.0);
            let pipe_power = 103.0 + if motor == CellId(20) { 5.0 } else { 0.0 };
            let node_power = 1.0 + if motor == CellId(0) { 5.0 } else { 0.0 };
            for (cell, power) in [
                (CellId(20), pipe_power),
                (CellId(0), node_power),
                (CellId(1), 0.0),
            ] {
                let expected = 300.0 + dt * power / (old.cells()[&cell].mass_kg * CP);
                near(
                    solved.thermal.temperature_k(cell).unwrap(),
                    expected,
                    3e-9,
                    0.0,
                );
            }
            near(
                report.thermal.thermal.stored_energy_change_j,
                dt * 109.0,
                3e-7,
                0.0,
            );
            assert!(
                report.equipment.electrical_boundary_residual_j.abs()
                    <= report.equipment.electrical_boundary_tolerance_j
            );
        }
    });
}

// Eliminate the three original caloric inventories and use independent
// 2x2 annular momentum plus bisection of the original nonlinear head law.
fn reference_step(data: &Calibration, old_v: [f64; 2], old_t: [f64; 3], dt: f64) -> Reference {
    let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
    let evaluate = |t: f64| {
        let offset = momentum_reference(old_v, dt, 0.0, t).1;
        let g = momentum_reference([0.0; 2], dt, 1.0, t).1;
        let residual = |q: f64| q - offset - g * head(q);
        let (mut lo, mut hi) = (0.0, 0.001);
        assert!(residual(lo) <= 0.0 && residual(hi) > 0.0);
        for _ in 0..90 {
            let q = f64::midpoint(lo, hi);
            if residual(q) > 0.0 {
                hi = q;
            } else {
                lo = q;
            }
        }
        let (velocity, q, heat, kinetic) =
            momentum_reference(old_v, dt, head(f64::midpoint(lo, hi)), t);
        let moved = dt * RHO * q;
        let upstream = (masses[0] * CP * old_t[0] + moved * CP * INLET + dt * 6.0)
            / ((masses[0] + moved) * CP);
        let pump = 0.5 * (shaft(data, q, t, 0.0) - head(q) * q);
        let residual = masses[1] * CP * (t - old_t[1])
            - moved * CP * (upstream - t)
            - heat
            - dt * (pump + 3.0);
        (residual, velocity, q, heat, kinetic, upstream, moved)
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
fn g1_source_input_and_loss_heat_close_original_viscosity_and_transport_feedback() {
    with_cx(|cx| {
        let original = calibration(5.0, 0.0);
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
            0.125,
            &[0.0625],
            controls(),
            permissive(),
            |_, time, candidate, source| {
                let t = candidate.temperature_k(CellId(20)).unwrap();
                assert_eq!(source.graph().branches()[0].pipe.temperature_k(), t);
                Ok::<_, ()>(drive(&original, t, time, 0.0, CellId(20), CellId(0)))
            },
            |_, candidate, flow, _| Ok::<_, ()>(extra_heat(candidate, flow)),
        )
        .unwrap();
        assert_eq!(
            solved.report.work.completed_head_evaluations,
            original.queries.load(Ordering::Relaxed)
        );
        let mut velocity = [0.0; 2];
        let mut temperature = [300.0; 3];
        for (report, sample) in solved.report.steps.iter().zip(&solved.report.drive_samples) {
            let dt = sample.end_s - sample.start_s;
            let expected = reference_step(&original, velocity, temperature, dt);
            for (actual, expected) in report.thermal.heat.into_iter().zip(expected.temperatures) {
                near(actual, expected, 3e-8, 0.0);
            }
            let power = shaft(
                &original,
                expected.mean_flow,
                expected.temperatures[1],
                sample.sample_time_s,
            );
            let hydraulic = head(expected.mean_flow) * expected.mean_flow;
            let energy = report.equipment.branches[0].unwrap();
            near(energy.shaft_input_j, power * dt, 1e-9, 1e-8);
            near(energy.electrical_input_j, (power + 20.0) * dt, 1e-9, 1e-8);
            near(
                energy.pump_heat_j,
                0.5 * (power - hydraulic) * dt,
                1e-9,
                1e-8,
            );
            near(energy.drive_heat_j, 5.0 * dt, 1e-12, 0.0);
            near(report.hydraulics.viscous_heat_j, expected.heat, 1e-11, 1e-8);
            near(
                report.hydraulics.kinetic_energy_change_j,
                expected.kinetic,
                1e-11,
                1e-8,
            );
            near(report.other_heat_j, 4.0 * dt, 1e-12, 0.0);
            let stream = RHO * expected.mean_flow * CP * (INLET - expected.temperatures[2]) * dt;
            let masses = [0.04, RHO * core::f64::consts::PI * R * R * L, 0.03];
            let stored = (0..3)
                .map(|i| masses[i] * CP * (expected.temperatures[i] - temperature[i]))
                .sum::<f64>();
            near(
                stored,
                stream + expected.heat + (0.5 * (power - hydraulic) + 9.0) * dt,
                1e-8,
                1e-10,
            );
            near(
                report.thermal.thermal.stored_energy_change_j,
                stored,
                3e-7,
                1e-8,
            );
            assert!(
                report.equipment.electrical_boundary_residual_j.abs()
                    <= report.equipment.electrical_boundary_tolerance_j
            );
            velocity = expected.velocity;
            temperature = expected.temperatures;
        }
        for (&actual, expected) in solved.momentum.branches()[0]
            .velocities_m_s()
            .iter()
            .zip(velocity)
        {
            near(actual, expected, 1e-10, 1e-8);
        }
    });
}

#[test]
fn g3_input_work_is_resolved_even_when_external_cooling_keeps_all_states_constant() {
    with_cx(|cx| {
        let original = calibration(0.0, 1e5);
        let caloric = Caloric;
        let old = inventory(&caloric, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let run = |settings| {
            advance_resolved_with_drive(
                cx,
                &owner,
                &old,
                &state,
                &mut owner.momentum().workspace().unwrap(),
                0.125,
                &[],
                controls(),
                settings,
                |_, time| Ok::<_, ()>(drive(&original, 300.0, time, -800.0, CellId(20), CellId(0))),
                |_, _, flow, _| {
                    let power = flow.flow_heads[0].as_ref().unwrap().response.power.unwrap();
                    let mut heat = exchange(flow, 300.0);
                    heat.heat
                        .bulk
                        .sources_w
                        .insert(CellId(20), -0.5 * power.state.pump_loss_w);
                    heat.heat
                        .bulk
                        .sources_w
                        .insert(CellId(0), -0.25 * power.state.drive_loss_w);
                    Ok::<_, ()>(heat)
                },
            )
            .unwrap()
        };
        let coarse = run(permissive());
        let fine = run(RegionalCarreauResolutionControls {
            absolute_transfer_error_j: 0.01,
            relative_transfer_error: 0.0,
            ..permissive()
        });
        assert_eq!(coarse.report.steps.len(), 2);
        assert!(fine.report.work.rejected_intervals > 0);
        assert!(fine.report.steps.len() > coarse.report.steps.len());
        let input = |solution: &ResolvedRegionalCarreau<'_, ()>| {
            solution
                .report
                .steps
                .iter()
                .map(|r| r.equipment.electrical_input_j)
                .sum::<f64>()
        };
        let exact = 220.0 * 0.125 + 1e5 * 0.125_f64.powi(3) / 3.0;
        assert!((input(&fine) - exact).abs() < (input(&coarse) - exact).abs() / 20.0);
        for solution in [&coarse, &fine] {
            assert_eq!(solution.momentum.branches()[0].velocities_m_s(), &[0.0; 2]);
            for cell in [CellId(0), CellId(20), CellId(1)] {
                near(
                    solution.thermal.temperature_k(cell).unwrap(),
                    300.0,
                    1e-10,
                    0.0,
                );
            }
        }
    });
}

#[test]
fn g4_missing_foreign_and_late_loss_destinations_preserve_both_histories_and_retry() {
    with_cx(|cx| {
        let original = calibration(0.0, 0.0);
        let caloric = Caloric;
        let old = inventory(&caloric, &[300.0; 3]);
        let owner = model(cx, &old, 1, false);
        let state = owner.momentum().rest_state(&[None; 2], cx).unwrap();
        let initial = snapshot(&state);
        let cells = old.cells().clone();
        let good = drive(&original, 300.0, 0.01, -800.0, CellId(20), CellId(0));
        let mut workspace = owner.momentum().workspace().unwrap();
        for destinations in [
            None,
            Some(vec![]),
            Some(vec![EquipmentHeatDestinations {
                pump: Some(CellId(999)),
                drive: Some(CellId(0)),
            }]),
        ] {
            let mut bad = good.clone();
            bad.branch_heat_destinations = destinations;
            let result = owner.advance_with_drive(
                cx,
                &old,
                &state,
                &mut workspace,
                0.02,
                &bad,
                controls(),
                &[0.0; 3],
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
            );
            assert!(matches!(result, Err(ImplicitTransportError::Invalid(_))));
            assert_eq!(snapshot(&state), initial);
            assert_eq!(old.cells(), &cells);
        }
        let late = advance_resolved_with_drive(
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
                Ok::<_, ()>(drive(
                    &original,
                    300.0,
                    time,
                    -800.0,
                    if time < 0.0625 {
                        CellId(20)
                    } else {
                        CellId(999)
                    },
                    CellId(0),
                ))
            },
            |_, _, flow, _| Ok::<_, ()>(exchange(flow, 300.0)),
        );
        assert!(matches!(
            late,
            Err(RegionalCarreauResolutionError::Step {
                start_s: 0.0625,
                source: ImplicitTransportError::Invalid(_),
                ..
            })
        ));
        assert_eq!(snapshot(&state), initial);
        assert_eq!(old.cells(), &cells);
        assert_eq!(old.time_s(), 0.0);
        let retry = owner
            .advance_with_drive(
                cx,
                &old,
                &state,
                &mut workspace,
                0.02,
                &good,
                controls(),
                &[0.0; 3],
                |_, _, flow| Ok::<_, ()>(exchange(flow, 300.0)),
            )
            .unwrap();
        assert_eq!(retry.thermal.time_s(), 0.02);
        assert!(retry.thermal.temperature_k(CellId(20)).unwrap() > 300.0);
    });
}
