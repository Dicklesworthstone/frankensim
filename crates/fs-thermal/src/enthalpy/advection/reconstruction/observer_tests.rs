//! Accepted explicit transport, shared stability and physical observer coupling.
use super::*;
use crate::enthalpy::{DepositionObserver, HeatTransfer};
use std::sync::Arc;

#[path = "schedule_tests.rs"]
mod schedule;

fn config(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-9,
    }
}

#[test]
fn explicit_limited_pulse_stays_bounded_and_reduces_diffusion_in_both_directions() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let n = 24;
    let initial_h: Vec<_> = (0..n).map(|i| if i < 8 { 100.0 } else { 0.0 }).collect();
    let initial = ring(&law, &initial_h);
    let dt = 1.0 / 256.0;
    let target = 32.0 * dt; // Exact transport moves the pulse three whole cells.
    let decay = n as f64;
    // Supply maxima which give the SAME accepted Euler dt for both spatial
    // schemes, after their respective native harmonic transport restrictions.
    let donor_config = config(dt / (1.0 - dt * decay / 0.8));
    let limited_config = config(dt / (1.0 - 2.0 * dt * decay / 0.8));
    for rate in [-1.0, 1.0] {
        let donor = ring_flow(n, rate);
        let limited = donor
            .clone()
            .with_reconstruction(ring_plan(&initial))
            .unwrap();
        let (low, low_flow, _) = initial
            .advance_observed(&donor, target, donor_config, || Ok(()))
            .unwrap();
        let (high, high_flow, audit) = initial
            .advance_observed(&limited, target, limited_config, || Ok(()))
            .unwrap();
        near(low_flow.audit().maximum_courant, dt * decay, 1e-14);
        near(high_flow.audit().maximum_courant, dt * decay, 1e-14);
        let mut low_error = 0.0;
        let mut high_error = 0.0;
        for i in 0..n {
            let origin = if rate > 0.0 {
                (i + n - 3) % n
            } else {
                (i + 3) % n
            };
            let exact = initial_h[origin];
            let a = low.cells()[&CellId(i as u64)].specific_h_j_kg;
            let b = high.cells()[&CellId(i as u64)].specific_h_j_kg;
            assert!((-1e-12..=100.0 + 1e-12).contains(&b));
            low_error += (a - exact).abs() / n as f64;
            high_error += (b - exact).abs() / n as f64;
        }
        assert!(
            high_error < 0.9 * low_error,
            "limited={high_error}, donor={low_error}"
        );
        near(high.energy_j(), initial.energy_j(), 1e-11);
        near(audit.stored_energy_change_j, 0.0, 1e-11);
        assert_eq!(audit.external_heat_j, 0.0);
        assert_eq!(high_flow.audit().external_enthalpy_j, 0.0);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(limited.audit(), Default::default());
    }
}

#[test]
fn native_conduction_and_reconstructed_convection_share_one_stability_allowance() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let geometry = ring(&law, &[100.0; 4]);
    let initial = prepare(
        &law,
        geometry
            .cells()
            .values()
            .cloned()
            .map(|mut cell| {
                // Original exposed area=1, m=1/4 and dT/dh=1: a=6 /s.
                cell.boundary.htc_w_m2_k = 1.5;
                cell
            })
            .collect(),
    );
    let flow = ring_flow(4, 1.0)
        .with_reconstruction(ring_plan(&initial))
        .unwrap();
    let (next, observed, energy) = initial
        .advance_observed(&flow, 0.15, config(1.0), || Ok(()))
        .unwrap();
    // b=4 /s; a and the two reconstructed face budgets ADD. An independent
    // min of conductive and advective limits would admit too large a step.
    let expected_dt = 0.8 / (6.0 + 2.0 * 4.0);
    near(observed.audit().maximum_courant, 4.0 * expected_dt, 1e-14);
    assert_eq!(energy.substeps, 3);
    // Uniform convection contributes no heat: the exact discrete result is
    // the original ambient Euler decay on those three admitted durations.
    let expected =
        100.0 * (1.0 - 6.0 * expected_dt).powi(2) * (1.0 - 6.0 * (0.15 - 2.0 * expected_dt));
    for cell in next.cells().values() {
        near(cell.specific_h_j_kg, expected, 1e-11);
    }
    near(energy.stored_energy_change_j, expected - 100.0, 1e-11);
    near(energy.external_heat_j, expected - 100.0, 1e-11);
}

#[test]
fn absent_and_empty_plans_preserve_original_arithmetic_and_entropy_selection() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let initial = ring(&law, &[20.0, 80.0, 100.0, 50.0]);
    let flow = ring_flow(4, -0.5).with_entropy_audit();
    let empty = LimitedEnthalpyReconstruction::try_new(&initial, &[], 0, || Ok(())).unwrap();
    let (expected, expected_flow, expected_energy) = initial
        .advance_observed(&flow, 0.2, config(0.1), || Ok(()))
        .unwrap();
    for plan in [None, Some(Arc::new(empty))] {
        let wrapped = flow.clone().with_optional_reconstruction(plan).unwrap();
        assert!(wrapped.reconstruction().is_none());
        let (actual, observed, energy) = initial
            .advance_observed(&wrapped, 0.2, config(0.1), || Ok(()))
            .unwrap();
        assert_eq!(actual.cells(), expected.cells());
        assert_eq!(energy, expected_energy);
        assert_eq!(observed.audit(), expected_flow.audit());
        assert_eq!(
            observed.advection().entropy_audit(),
            expected_flow.entropy_audit()
        );
    }
    // Existing upwind entropy histories cannot silently describe changed faces.
    assert!(flow.with_reconstruction(ring_plan(&initial)).is_err());
}

#[derive(Clone, Default)]
struct HeatingHistory {
    joules: f64,
    intervals: Vec<[f64; 2]>,
    events: usize,
    first_h: Option<[f64; 4]>,
    cancel_after_s: Option<f64>,
}

impl HeatingHistory {
    fn watts(time_s: f64) -> f64 {
        if time_s < 0.03 { 1.0 } else { 3.0 }
    }
}

impl DepositionObserver for HeatingHistory {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.joules = 0.0;
        self.intervals.clear();
        self.events = 0;
        self.first_h = None;
        Ok(())
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        Ok((thermal.time_s() < 0.03).then_some(0.03))
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok([(CellId(0), Self::watts(thermal.time_s()))].into())
    }
    fn heat_transfers(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        Ok(vec![HeatTransfer {
            a: CellId(0),
            b: CellId(1),
            power_w: 2.0,
        }])
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        _: f64,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok([(CellId(1), 4.0)].into())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events += 1;
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.joules += (Self::watts(start_s) + 4.0) * dt_s;
        self.intervals.push([start_s, thermal.time_s()]);
        self.first_h.get_or_insert_with(|| {
            std::array::from_fn(|i| thermal.cells()[&CellId(i as u64)].specific_h_j_kg)
        });
        if self
            .cancel_after_s
            .is_some_and(|time| thermal.time_s() >= time)
        {
            return Err(DepositionError::Cancelled);
        }
        Ok(())
    }
}

#[test]
fn explicit_reconstruction_preserves_inner_heat_transitions_events_and_rollback() {
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
    let flow = Advected::new(
        HeatingHistory::default(),
        ring_flow(4, 0.1).flows().to_vec(),
        4,
    )
    .unwrap()
    .with_reconstruction(ring_plan(&initial))
    .unwrap();
    assert!(initial.power_balance_observed(&flow, || Ok(())).is_err());
    let (next, observed, energy) = initial
        .advance_observed(&flow, 0.1, config(1.0), || Ok(()))
        .unwrap();
    let history = observed.inner();
    assert_eq!(
        history.intervals,
        vec![[0.0, 0.03], [0.03, 0.05], [0.05, 0.1]]
    );
    assert_eq!(history.events, 1);
    near(history.joules, 0.03 * 5.0 + 0.07 * 7.0, 1e-14);
    near(energy.external_heat_j, history.joules + 0.05 * 2.0, 1e-14);
    near(
        next.energy_j() - initial.energy_j(),
        energy.external_heat_j,
        1e-11,
    );
    // Before the first source transition all native advective powers vanish.
    // The inner owner's 2 W paired transfer must remain a separate prefix.
    let first = history.first_h.unwrap();
    near(first[0], 100.0 + 0.03 * (1.0 - 2.0) / 0.25, 1e-12);
    near(first[1], 100.0 + 0.03 * (4.0 + 2.0) / 0.25, 1e-12);
    near(first[2], 100.0, 0.0);
    near(first[3], 100.0, 0.0);
    let failing = Advected::new(
        HeatingHistory {
            cancel_after_s: Some(0.06),
            ..HeatingHistory::default()
        },
        flow.flows().to_vec(),
        4,
    )
    .unwrap()
    .with_reconstruction(ring_plan(&initial))
    .unwrap();
    assert!(matches!(
        initial.advance_observed(&failing, 0.1, config(1.0), || Ok(())),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(failing.audit(), Default::default());
    assert!(failing.inner().intervals.is_empty());
    assert_eq!(failing.inner().joules, 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), geometry.cells());
}

#[test]
fn unsupported_reconstructed_face_refuses_before_publishing_a_thermal_step() {
    let law = Law {
        datum: 0.0,
        hole: true,
    };
    let initial = ring(&law, &[0.0, 100.0, 200.0]);
    let flow = ring_flow(3, 1.0)
        .with_reconstruction(ring_plan(&initial))
        .unwrap();
    assert!(matches!(
        initial.advance_observed(&flow, 0.1, config(1.0), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 100.0);
    assert_eq!(flow.audit(), Default::default());
}
