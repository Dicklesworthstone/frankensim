//! Prescribed transport histories on the same accepted explicit time grid.
use super::*;
use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};

fn history(inlet: Option<u64>, knots: &[(f64, f64, f64)]) -> Arc<AdvectionSchedule> {
    Arc::new(
        AdvectionSchedule::new(
            knots
                .iter()
                .map(|&(time_s, flow_multiplier, h)| AdvectionScheduleKnot {
                    time_s,
                    flow_multiplier,
                    inlet_enthalpies_j_kg: inlet.into_iter().map(|id| (id, h)).collect(),
                })
                .collect(),
            knots.len(),
            usize::from(inlet.is_some()),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn open<O>(inner: O, n: usize, rate: f64, h: f64) -> Advected<O> {
    Advected::new(
        inner,
        (0..n - 1)
            .map(|i| MassFlow {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                kg_s: rate,
            })
            .collect(),
        n - 1,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: n as u64,
                cell: CellId(0),
                kg_s: rate,
                material_id: MaterialId(1),
                specific_h_j_kg: h,
            },
            MassBoundary::Outlet {
                id: n as u64 + 1,
                cell: CellId(n as u64 - 1),
                kg_s: rate,
            },
        ],
        n + 1,
    )
    .unwrap()
}

#[test]
fn explicit_pump_ramp_refines_to_independent_non_autonomous_calorimetry() {
    for datum in [0.0, 1_048_576.0] {
        let law = Law { datum, hole: false };
        let initial = prepare(
            &law,
            vec![cell(&law, 0, 0.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)],
        );
        let base = open((), 1, 1.0, datum + 100.0);
        let schedule = history(
            Some(1),
            &[(0.0, 1.0, datum + 100.0), (1.0, 2.0, datum + 100.0)],
        );
        let flow = base.clone().with_schedule(schedule).unwrap();
        let exact = 100.0 * (1.0 - (-1.5_f64).exp());
        let mut errors = Vec::new();
        for maximum_dt in [0.1, 0.05, 0.025] {
            let (next, observed, energy) = initial
                .advance_observed(
                    &flow,
                    1.0,
                    StepConfig {
                        max_dt_s: maximum_dt,
                        energy_tolerance_j: 1e-7,
                    },
                    || Ok(()),
                )
                .unwrap();
            let h = next.cells()[&CellId(0)].specific_h_j_kg - datum;
            // m=1 and q(t)=1+t: h'=q(t)*(100-h), with h(0)=0.
            // Its integrating factor is exp(t+t^2/2), independently of dt.
            errors.push((h - exact).abs());
            near(energy.external_heat_j, h, 1e-7);
            near(observed.audit().external_enthalpy_j, h, 1e-7);
            near(
                observed.audit().inlet_mass_kg,
                observed.audit().outlet_mass_kg,
                1e-14,
            );
            assert!(observed.audit().maximum_courant <= 0.8);
        }
        assert!(
            errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
            "{errors:?}"
        );
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(flow.audit(), Default::default());
        assert_eq!(flow.boundaries(), base.boundaries());
    }
}

#[test]
fn scheduled_inlet_stop_restart_and_sources_use_the_same_native_intervals() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let geometry = ring(&law, &[100.0; 4]);
    let mut events: Vec<_> = geometry
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    let mut boundary = geometry.cells()[&CellId(2)].boundary;
    boundary.source_w = 2.0;
    events.push(TimedEvent {
        time_s: 0.05,
        event: Event::Boundary(CellId(2), boundary),
    });
    let initial = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 8,
            max_events: 5,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let schedule = history(
        Some(4),
        &[
            (0.0, 1.0, 100.0),
            (0.04, 0.0, 100.4),
            (0.06, 0.0, 100.6),
            (0.08, 2.0, 100.8),
            (0.1, 2.0, 101.0),
        ],
    );
    let flow = open(HeatingHistory::default(), 4, 0.1, 100.0)
        .with_schedule(schedule.clone())
        .unwrap();
    let (next, observed, energy) = initial
        .advance_observed(&flow, 0.1, config(1.0), || Ok(()))
        .unwrap();
    let intervals: [[f64; 2]; 6] = [
        [0.0, 0.03],
        [0.03, 0.04],
        [0.04, 0.05],
        [0.05, 0.06],
        [0.06, 0.08],
        [0.08, 0.1],
    ];
    assert_eq!(observed.inner().intervals, intervals);
    assert_eq!(observed.inner().events, 1);
    let mut expected = [100.0; 4];
    let mut external = 0.0;
    let mut mass = 0.0;
    for [start, end] in intervals {
        let dt = end - start;
        // Independent authored pump and inlet laws, on each actual START.
        let multiplier = if start < 0.04 {
            1.0 - start / 0.04
        } else if start < 0.08 {
            100.0 * (start - 0.06).max(0.0)
        } else {
            2.0
        };
        let rate = 0.1 * multiplier;
        let inlet_h = 100.0 + 10.0 * start;
        let heat = [
            HeatingHistory::watts(start) - 2.0,
            6.0,
            if start >= 0.05 { 2.0 } else { 0.0 },
            0.0,
        ];
        let old = expected;
        for i in 0..4 {
            let incoming = if i == 0 { inlet_h } else { old[i - 1] };
            expected[i] += dt * (rate * (incoming - old[i]) + heat[i]) / 0.25;
        }
        external += dt * rate * (inlet_h - old[3]);
        mass += dt * rate;
    }
    for (i, h) in expected.into_iter().enumerate() {
        near(next.cells()[&CellId(i as u64)].specific_h_j_kg, h, 1e-11);
    }
    near(observed.audit().external_enthalpy_j, external, 1e-13);
    near(observed.audit().inlet_mass_kg, mass, 1e-14);
    near(observed.audit().outlet_mass_kg, mass, 1e-14);
    near(observed.audit().transported_mass_kg, 5.0 * mass, 1e-14);
    near(
        energy.external_heat_j,
        external + observed.inner().joules + 0.1,
        1e-13,
    );
    near(energy.stored_energy_change_j, energy.external_heat_j, 1e-11);
    assert!(observed.flows().iter().all(|flow| flow.kg_s == 0.2));
    assert_eq!(observed.boundaries().len(), 2);
    // A late source refusal discards already integrated stopped/running periods.
    let failing = open(
        HeatingHistory {
            cancel_after_s: Some(0.09),
            ..HeatingHistory::default()
        },
        4,
        0.1,
        100.0,
    )
    .with_schedule(schedule)
    .unwrap();
    assert!(matches!(
        initial.advance_observed(&failing, 0.1, config(1.0), || Ok(())),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(failing.audit(), Default::default());
    assert!(failing.inner().intervals.is_empty());
    assert!(failing.flows().iter().all(|flow| flow.kg_s == 0.1));
    assert_eq!(initial.cells(), geometry.cells());
}

#[test]
fn shared_schedule_uses_original_limited_flux_and_preserves_direction() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let initial = ring(&law, &[0.0, 40.0, 100.0, 80.0, 20.0]);
    let plan = Arc::new(ring_plan(&initial));
    for rate in [-0.5, 0.5] {
        for multiplier in [0.0, 2.0] {
            let base = ring_flow(5, rate);
            let schedule = history(None, &[(0.0, multiplier, 0.0)]);
            let sampled = schedule
                .flow_at_time(&base, &initial, 0.2, || Ok(()))
                .unwrap();
            let scheduled = base
                .with_reconstruction(plan.clone())
                .unwrap()
                .with_schedule(schedule.clone())
                .unwrap();
            let fixed = ring_flow(5, rate * multiplier)
                .with_reconstruction(plan.clone())
                .unwrap();
            assert_eq!(sampled.flows(), fixed.flows());
            let (a, observer_a, energy_a) = initial
                .advance_observed(&scheduled, 0.2, config(0.05), || Ok(()))
                .unwrap();
            let (b, observer_b, energy_b) = initial
                .advance_observed(&fixed, 0.2, config(0.05), || Ok(()))
                .unwrap();
            assert_eq!(a.cells(), b.cells());
            assert_eq!(energy_a, energy_b);
            assert_eq!(observer_a.audit(), observer_b.audit());
            assert!(scheduled.with_schedule(schedule).is_err());
        }
    }
    let unknown = history(Some(99), &[(0.0, 1.0, 100.0)]);
    assert!(open((), 1, 1.0, 100.0).with_schedule(unknown).is_err());
}

#[test]
fn scheduled_upwind_entropy_uses_the_actual_inlet_and_flux_sample() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let initial = prepare(
        &law,
        vec![cell(&law, 0, 10.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)],
    );
    let scheduled = open((), 1, 1.0, 90.0)
        .with_entropy_audit()
        .with_schedule(history(Some(1), &[(0.0, 2.0, 110.0)]))
        .unwrap();
    let fixed = open((), 1, 2.0, 110.0).with_entropy_audit();
    let (a, observer_a, energy_a) = initial
        .advance_observed(&scheduled, 0.2, config(0.05), || Ok(()))
        .unwrap();
    let (b, observer_b, energy_b) = initial
        .advance_observed(&fixed, 0.2, config(0.05), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(energy_a, energy_b);
    assert_eq!(observer_a.audit(), observer_b.audit());
    assert_eq!(
        observer_a.advection().entropy_audit(),
        observer_b.entropy_audit()
    );
}

#[test]
fn a_supported_knot_does_not_authorize_an_unsupported_interpolated_inlet() {
    let law = Law {
        datum: 0.0,
        hole: true,
    };
    let initial = prepare(
        &law,
        vec![cell(&law, 0, 0.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)],
    );
    let flow = open((), 1, 1.0, 100.0)
        .with_schedule(history(Some(1), &[(0.0, 1.0, 100.0), (1.0, 1.0, 200.0)]))
        .unwrap();
    assert!(matches!(
        initial.advance_observed(&flow, 1.0, config(0.25), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    assert_eq!(flow.audit(), Default::default());
}
