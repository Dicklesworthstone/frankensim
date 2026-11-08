//! G1/G3/G4: latent transport accuracy, exact forcing boundaries and refusal.

#[path = "adaptive_tests.rs"]
mod adaptive;

use super::*;
use crate::enthalpy::advection::reconstruction::{AxialCellStencil, LimitedEnthalpyReconstruction};
use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};
use crate::enthalpy::advection::{MassBoundary, MassFlow};
use crate::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, Event, MaterialId, ThermalLaw,
    TimedEvent,
};

struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(-100.0..=300.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0
            + if h < 0.0 {
                h / 10.0
            } else if h > 100.0 {
                (h - 100.0) / 10.0
            } else {
                0.0
            })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}

fn cell(id: usize, h: f64, a: [f64; 3], b: [f64; 3], mass: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id as u64),
        segment_id: id as u64,
        body_id: 0,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: a,
            end_m: b,
            volume_m3: mass,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn config(dt: f64) -> Ssprk2Controls {
    Ssprk2Controls {
        step: StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-8,
        },
        maximum_stage_evaluations: 100_000,
    }
}

fn prepare<'a>(
    law: &'a dyn ThermalLaw,
    events: Vec<TimedEvent>,
    count: usize,
) -> DepositionSimulation<'a> {
    let maximum = events.len();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: count,
            max_links: 2 * count + 4,
            max_events: maximum,
            max_substeps: 20_000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0).step, || Ok(()))
    .unwrap()
    .0
}

fn births(cells: impl IntoIterator<Item = DepositedCell>) -> Vec<TimedEvent> {
    cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect()
}

fn ring<'a>(law: &'a dyn ThermalLaw, h: &[f64]) -> DepositionSimulation<'a> {
    let n = h.len();
    let radius = 1.0 / (2.0 * n as f64 * (std::f64::consts::PI / n as f64).sin());
    let vertices: Vec<_> = (0..n)
        .map(|i| {
            let angle = std::f64::consts::TAU * i as f64 / n as f64;
            [radius * angle.cos(), radius * angle.sin(), 0.0]
        })
        .collect();
    prepare(
        law,
        births(
            h.iter()
                .enumerate()
                .map(|(i, &h)| cell(i, h, vertices[i], vertices[(i + 1) % n], 1.0 / n as f64)),
        ),
        n,
    )
}

fn flow(n: usize, rate: f64) -> Advected<()> {
    Advected::new(
        (),
        (0..n)
            .map(|i| MassFlow {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(((i + 1) % n) as u64),
                kg_s: rate,
            })
            .collect(),
        n,
    )
    .unwrap()
}

fn plan(state: &DepositionSimulation<'_>) -> LimitedEnthalpyReconstruction {
    let n = state.cells().len();
    LimitedEnthalpyReconstruction::try_new(
        state,
        &(0..n)
            .map(|i| AxialCellStencil {
                previous: CellId(((i + n - 1) % n) as u64),
                cell: CellId(i as u64),
                next: CellId(((i + 1) % n) as u64),
            })
            .collect::<Vec<_>>(),
        n,
        || Ok(()),
    )
    .unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, allowance {tolerance:e}"
    );
}

#[test]
fn latent_ring_has_second_order_time_accuracy_in_both_directions() {
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
        for count in [16, 32, 64] {
            let dt = target / count as f64;
            // Invert the existing harmonic max-dt/transport bound, so this
            // experiment refines the accepted time interval by exactly two.
            let controls = config(dt / (1.0 - dt * n as f64 / 0.8));
            let (next, accepted, report) = owner
                .advance_ssprk2_to(&initial, target, controls, || Ok(()))
                .unwrap();
            let attenuation = (n as f64 * (theta.cos() - 1.0) * target).exp();
            let mut error = 0.0;
            for (i, cell) in next.cells().values().enumerate() {
                let expected = 50.0
                    + 20.0
                        * attenuation
                        * (i as f64 * theta - rate * n as f64 * theta.sin() * target).cos();
                error += (cell.specific_h_j_kg - expected).powi(2);
                assert_eq!(law.temperature_k(cell.specific_h_j_kg).unwrap(), 300.0);
                assert!((0.0..=100.0).contains(&cell.specific_h_j_kg));
            }
            errors.push((error / n as f64).sqrt());
            near(next.energy_j(), initial.energy_j(), 2e-13);
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.residual_j, 0.0, 2e-13);
            near(
                accepted.audit().transported_mass_kg,
                n as f64 * target,
                2e-14,
            );
            assert_eq!(report.rejected_intervals, 0);
            assert_eq!(report.stage_evaluations, 2 * report.energy.substeps);
        }
        for pair in errors.windows(2) {
            assert!(
                (3.7..4.3).contains(&(pair[0] / pair[1])),
                "errors {errors:?}"
            );
        }
    }
    assert_eq!(initial.time_s(), 0.0);
}

fn startup() -> ReconstructedAdvected<()> {
    let owner = Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 0,
                    cell: CellId(0),
                    kg_s: 1.0,
                    material_id: MaterialId(1),
                    specific_h_j_kg: 80.0,
                },
                MassBoundary::Outlet {
                    id: 1,
                    cell: CellId(0),
                    kg_s: 1.0,
                },
            ],
            2,
        )
        .unwrap();
    owner
        .with_schedule(
            AdvectionSchedule::new(
                [0.0, 0.05, 0.1]
                    .into_iter()
                    .map(|time_s| AdvectionScheduleKnot {
                        time_s,
                        flow_multiplier: 1000.0 * time_s,
                        inlet_enthalpies_j_kg: BTreeMap::new(),
                    })
                    .collect(),
                3,
                0,
                || Ok(()),
            )
            .unwrap(),
        )
        .unwrap()
}

#[test]
fn pump_startup_rechecks_stage_cfl_and_integrates_actual_mass_and_heat() {
    let law = Latent;
    let initial = prepare(
        &law,
        births([cell(0, 20.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0)]),
        1,
    );
    let owner = startup();
    let exact = 80.0 - 60.0 * (-5.0_f64).exp();
    let mut errors = Vec::new();
    for dt in [0.1, 0.001] {
        let mut queries = Vec::new();
        let (next, accepted, report) = owner
            .advance_ssprk2_with_coefficients(
                &initial,
                0.1,
                config(dt),
                |state| {
                    queries.push(state.time_s());
                    unchanged(state)
                },
                || Ok(()),
            )
            .unwrap();
        let h = next.cells()[&CellId(0)].specific_h_j_kg;
        assert!((20.0..=80.0).contains(&h));
        assert_eq!(law.temperature_k(h).unwrap(), 300.0);
        assert!(
            report.rejected_intervals > 0,
            "pump acceleration must tighten the second stage"
        );
        assert!(
            queries.iter().filter(|&&t| t == 0.05).count() >= 2,
            "schedule knot must end and restart a stage"
        );
        assert!(accepted.audit().maximum_courant <= 0.8);
        near(accepted.audit().inlet_mass_kg, 5.0, 2e-13);
        near(accepted.audit().outlet_mass_kg, 5.0, 2e-13);
        near(accepted.audit().transported_mass_kg, 10.0, 4e-13);
        near(accepted.audit().external_enthalpy_j, h - 20.0, 3e-12);
        near(report.energy.external_heat_j, h - 20.0, 3e-12);
        near(report.energy.residual_j, 0.0, 3e-12);
        errors.push((h - exact).abs());
    }
    assert!(errors[1] < errors[0] / 20.0, "refinement errors {errors:?}");
    assert_eq!(owner.audit(), AdvectionAudit::default());
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 20.0);
}

#[test]
fn source_quadrature_and_births_respect_the_two_sides_of_an_event() {
    let law = Latent;
    let mut first = cell(0, 50.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0);
    first.boundary.source_w = 1.0;
    let mut changed = first.boundary;
    changed.source_w = 3.0;
    let mut events = births([first]);
    events.extend([
        TimedEvent {
            time_s: 0.3,
            event: Event::Boundary(CellId(0), changed),
        },
        TimedEvent {
            time_s: 0.3,
            event: Event::Birth(cell(1, 20.0, [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], 1.0)),
        },
    ]);
    let initial = prepare(&law, events, 2);
    let owner = Advected::new((), vec![], 0).unwrap();
    let mut queries = Vec::new();
    let (next, _, report) = owner
        .advance_ssprk2_with_coefficients(
            &initial,
            1.0,
            config(1.0),
            |state| {
                queries.push((
                    state.time_s(),
                    state.cells()[&CellId(0)].boundary.source_w,
                    state.cells().len(),
                ));
                Ok(StepCoefficients {
                    sources_w: BTreeMap::from([(CellId(0), 2.0 * state.time_s())]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                }
                .into())
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        queries,
        [(0.0, 1.0, 1), (0.3, 1.0, 1), (0.3, 3.0, 2), (1.0, 3.0, 2)]
    );
    near(next.cells()[&CellId(0)].specific_h_j_kg, 53.4, 2e-14);
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, 20.0);
    near(report.energy.external_heat_j, 3.4, 1e-14);
    near(report.energy.incoming_energy_j, 20.0, 0.0);
    near(report.energy.stored_energy_change_j, 23.4, 2e-14);
    near(report.energy.residual_j, 0.0, 2e-14);
    assert_eq!(report.energy.events, 2);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(report.stage_evaluations, 4);
    assert_eq!(initial.cells().len(), 1);
    assert_eq!(initial.cells()[&CellId(0)].boundary.source_w, 1.0);
}

#[test]
fn reconstructed_latent_pulse_stays_bounded_and_is_less_diffusive_than_donor_flow() {
    let law = Latent;
    let n = 24;
    let h: Vec<_> = (0..n)
        .map(|i| if (8..14).contains(&i) { 100.0 } else { 0.0 })
        .collect();
    let initial = ring(&law, &h);
    let reconstruction = plan(&initial);
    let target = 3.0 / n as f64;
    for rate in [-1.0, 1.0] {
        let donor = flow(n, rate);
        let limited = donor
            .clone()
            .with_reconstruction(reconstruction.clone())
            .unwrap();
        let (upwind, _, _) = donor
            .advance_ssprk2_to(&initial, target, config(0.004), || Ok(()))
            .unwrap();
        let (next, accepted, report) = limited
            .advance_ssprk2_to(&initial, target, config(0.004), || Ok(()))
            .unwrap();
        let mut errors = [0.0; 2];
        for (slot, state) in [&upwind, &next].into_iter().enumerate() {
            for (i, cell) in state.cells().values().enumerate() {
                let origin = (i as isize - (3.0 * rate) as isize).rem_euclid(n as isize) as usize;
                errors[slot] += (cell.specific_h_j_kg - h[origin]).abs() / n as f64;
                assert!((-1e-12..=100.0 + 1e-12).contains(&cell.specific_h_j_kg));
                near(
                    law.temperature_k(cell.specific_h_j_kg).unwrap(),
                    300.0,
                    2e-13,
                );
            }
        }
        assert!(errors[1] < 0.85 * errors[0], "pulse errors {errors:?}");
        near(next.energy_j(), initial.energy_j(), 2e-13);
        near(report.energy.residual_j, 0.0, 2e-13);
        assert!(accepted.audit().maximum_courant <= 0.4);
        assert_eq!(report.rejected_intervals, 0);
    }
}

#[test]
fn stage_budget_and_unsupported_entropy_refuse_without_publishing_history() {
    let law = Latent;
    let initial = ring(&law, &[20.0, 60.0, 80.0]);
    let owner = flow(3, 1.0);
    let mut controls = config(0.1);
    controls.maximum_stage_evaluations = 1;
    assert!(matches!(
        owner.advance_ssprk2_to(&initial, 0.1, controls, || Ok(())),
        Err(DepositionError::Budget("SSPRK2 stage evaluations"))
    ));
    assert!(matches!(
        owner
            .clone()
            .with_entropy_audit()
            .advance_ssprk2_to(&initial, 0.1, config(0.1), || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
    let mut calls = 0;
    assert!(matches!(
        owner.advance_ssprk2_with_coefficients(
            &initial,
            0.1,
            config(0.1),
            |state| {
                calls += 1;
                if calls == 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    unchanged(state)
                }
            },
            || Ok(())
        ),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 20.0);
    assert_eq!(owner.audit(), AdvectionAudit::default());
}

#[test]
fn repeated_out_of_chart_predictors_are_bounded_before_second_stage_queries() {
    let law = Latent;
    let initial = prepare(
        &law,
        births([cell(0, 300.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0)]),
        1,
    );
    let owner = Advected::new((), vec![], 0).unwrap();
    let mut controls = config(1.0);
    controls.maximum_stage_evaluations = 3;
    let mut queries = 0;
    assert!(matches!(
        owner.advance_ssprk2_with_coefficients(
            &initial,
            1.0,
            controls,
            |_| {
                queries += 1;
                Ok(StepCoefficients {
                    sources_w: BTreeMap::from([(CellId(0), 1.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                }
                .into())
            },
            || Ok(())
        ),
        Err(DepositionError::Budget("SSPRK2 rejected intervals"))
    ));
    assert_eq!(queries, 1);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 300.0);
}
