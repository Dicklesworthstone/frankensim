//! Stage-time pump/inlet histories on the original conservative thermal solve.
use super::*;
use crate::enthalpy::advection::liquid::{LiquidCarrierKnot, LiquidEnthalpyCarrier};
use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};

fn history(inlet: u64, knots: &[(f64, f64, f64)]) -> AdvectionSchedule {
    AdvectionSchedule::new(
        knots.iter().map(|&(time_s, flow_multiplier, h)| AdvectionScheduleKnot {
            time_s,
            flow_multiplier,
            inlet_enthalpies_j_kg: BTreeMap::from([(inlet, h)]),
        }).collect(),
        knots.len(),
        1,
        || Ok(()),
    ).unwrap()
}

fn advance_one(h0: f64, start: f64, dt: f64, sdirk: bool, rate: impl Fn(f64) -> f64,
    inlet: impl Fn(f64) -> f64) -> (f64, f64, f64) {
    let gamma = 1.0 - std::f64::consts::FRAC_1_SQRT_2;
    if !sdirk {
        let r = rate(start + dt);
        let incoming = inlet(start + dt);
        let h = (h0 + dt * r * incoming) / (1.0 + dt * r);
        let boundary = dt * r * (incoming / h).ln();
        let mixing = dt * r * ((incoming - h) / h - (incoming / h).ln());
        return (h, boundary, mixing);
    }
    let t1 = start + gamma * dt;
    let r1 = rate(t1);
    let h1 = (h0 + gamma * dt * r1 * inlet(t1)) / (1.0 + gamma * dt * r1);
    let q1 = r1 * (inlet(t1) - h1);
    let r2 = rate(start + dt);
    let h2 = (h0 + (1.0 - gamma) * dt * q1 + gamma * dt * r2 * inlet(start + dt))
        / (1.0 + gamma * dt * r2);
    let mut boundary = 0.0;
    let mut mixing = 0.0;
    for (weight, r, h, incoming) in [
        (1.0 - gamma, r1, h1, inlet(t1)),
        (gamma, r2, h2, inlet(start + dt)),
    ] {
        boundary += weight * dt * r * (incoming / h).ln();
        mixing += weight * dt * r * ((incoming - h) / h - (incoming / h).ln());
    }
    (h2, boundary, mixing)
}

#[test]
fn prescribed_pump_and_inlet_use_actual_be_and_sdirk_stages_and_entropy() {
    for datum in [0.0, 1e9] {
        let law = ShiftedLaw { base: Law { nonlinear_k: false, plateau: false }, datum };
        let (thermal, conduction) = setup(&law, &[300.0 + datum], None, vec![]);
        let flow = open_flow(1, 2.0, 400.0 + datum);
        let schedule = history(1, &[(0.0, 0.5, 400.0 + datum), (1.0, 1.0, 500.0 + datum)]);
        for graph in [false, true] {
            let system = ImplicitConductionChain::prepare(
                &thermal, &conduction, &[(MaterialId(0), &law)],
                if graph { LinearTopology::Graph(ImplicitConductionLinearControls::default()) }
                else { LinearTopology::Chain }, || Ok(()),
            ).unwrap()
                .with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).unwrap()
                .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
            for sdirk in [false, true] {
                let mut policy = controls(2.0);
                policy.cell_tolerance_j = 1e-5;
                policy.total_tolerance_j = 1e-5;
                let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
                let (next, report) = if sdirk {
                    system.advance_sdirk2_to(&thermal, 1.0, policy, &mut workspace, || Ok(()))
                } else {
                    system.advance_to(&thermal, 1.0, policy, &mut workspace, || Ok(()))
                }.unwrap();
                let (expected, boundary, mixing) =
                    advance_one(300.0, 0.0, 1.0, sdirk, |t| 1.0 + t, |t| 400.0 + 100.0 * t);
                near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-6);
                near(report.energy.external_heat_j, expected - 300.0, 1e-6);
                near(report.energy.residual_j, 0.0, 1e-6);
                assert_eq!(report.energy.substeps, 1);
                let carried = report.entropy.advection.unwrap();
                near(carried.boundary_entropy_input_j_k, boundary, 1e-8);
                near(carried.upwind_mixing_j_k, mixing, 1e-8);
                near(carried.balance_residual_j_k, 0.0, 1e-8);
                assert!(!report.entropy.has_unclassified_sources);
            }
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn schedule_knots_clip_all_native_time_owners() {
    let law = ShiftedLaw { base: Law { nonlinear_k: false, plateau: false }, datum: 0.0 };
    let (thermal, conduction) = setup(&law, &[300.0], None, vec![]);
    let flow = open_flow(1, 2.0, 400.0);
    let schedule = history(1, &[(0.0, 0.5, 400.0), (0.4, 1.0, 440.0), (1.0, 0.25, 500.0)]);
    let system = ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
        .unwrap().with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).unwrap();
    assert_eq!(schedule.next_transition_time_s(0.0), Some(0.4));
    assert_eq!(schedule.next_transition_time_s(0.4), Some(1.0));
    assert_eq!(schedule.next_transition_time_s(1.0), None);
    for sdirk in [false, true] {
        let rate = |t| if t <= 0.4 { 1.0 + 2.5 * t } else { 2.0 - 2.5 * (t - 0.4) };
        let inlet = |t| 400.0 + 100.0 * t;
        let first = advance_one(300.0, 0.0, 0.4, sdirk, rate, inlet).0;
        let expected = advance_one(first, 0.4, 0.6, sdirk, rate, inlet).0;
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = if sdirk {
            system.advance_sdirk2_to(&thermal, 1.0, controls(2.0), &mut workspace, || Ok(()))
        } else {
            system.advance_to(&thermal, 1.0, controls(2.0), &mut workspace, || Ok(()))
        }.unwrap();
        assert_eq!(report.energy.substeps, 2);
        near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-9);
        near(report.energy.external_heat_j, expected - 300.0, 1e-9);
        let policy = adaptive::AdaptiveImplicitConductionControls {
            step: controls(2.0), initial_dt_s: 2.0, maximum_trials: 8,
            enthalpy_absolute_tolerance_j_kg: 1000.0,
            temperature_absolute_tolerance_k: 1000.0,
            phase_fraction_absolute_tolerance: None, relative_tolerance: 0.0,
        };
        let (_, adaptive) = if sdirk {
            system.advance_sdirk2_adaptive_to(&thermal, 1.0, policy, &mut workspace, || Ok(()))
        } else {
            system.advance_adaptive_to(&thermal, 1.0, policy, &mut workspace, || Ok(()))
        }.unwrap();
        assert_eq!(adaptive.accepted_comparisons, 2);
        assert_eq!(adaptive.solve.energy.substeps, 4);
    }
}

#[test]
fn scheduled_pump_refines_toward_the_independent_non_autonomous_solution() {
    let law = ShiftedLaw { base: Law { nonlinear_k: false, plateau: false }, datum: 0.0 };
    let (thermal, conduction) = setup(&law, &[300.0], None, vec![]);
    let flow = open_flow(1, 1.0, 400.0);
    let schedule = history(1, &[(0.0, 1.0, 400.0), (1.0, 2.0, 400.0)]);
    let system = ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
        .unwrap().with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).unwrap();
    // dh/dt = (1+t)(400-h), integrated exactly using integral(1+t) = t+t^2/2.
    let exact = 400.0 - 100.0 * (-1.5_f64).exp();
    for sdirk in [false, true] {
        let mut errors = Vec::new();
        for dt in [0.25, 0.125, 0.0625] {
            let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
            let (next, report) = if sdirk {
                system.advance_sdirk2_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            } else {
                system.advance_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            }.unwrap();
            let h = next.temperature_k(CellId(0)).unwrap();
            errors.push((h - exact).abs());
            near(report.energy.external_heat_j, h - 300.0, 1e-9);
        }
        let ratio = if sdirk { 0.3 } else { 0.6 };
        assert!(errors[1] < ratio * errors[0] && errors[2] < ratio * errors[1], "{errors:?}");
    }
}

fn liquid_carrier(law: &ShiftedLaw) -> LiquidEnthalpyCarrier {
    let knots = [(100.0, 0.0), (300.0, 0.0), (400.0, 0.5), (500.0, 1.0), (1000.0, 1.0)]
        .map(|(h, liquid_mass_fraction)| LiquidCarrierKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: law.temperature_k(h).unwrap(),
            liquid_mass_fraction,
            bulk_density_kg_m3: 1.0,
        });
    LiquidEnthalpyCarrier::from_source_chart(MaterialId(0), law, &knots, || Ok(())).unwrap()
}

#[test]
fn stopped_liquid_pump_can_heat_solid_before_admitting_actual_flow() {
    let law = ShiftedLaw { base: Law { nonlinear_k: false, plateau: true }, datum: 0.0 };
    let (mut thermal, conduction) = setup(&law, &[250.0], None, vec![]);
    let flow = open_flow(1, 2.0, 600.0).with_liquid_carrier(liquid_carrier(&law)).unwrap();
    let schedule = history(1, &[(0.0, 0.0, 600.0), (1.0, 0.0, 600.0), (2.0, 1.0, 600.0)]);
    let bind = |thermal: &DepositionSimulation<'_>| {
        ImplicitConductionChain::new(thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap().with_scheduled_implicit_advection(thermal, &flow, &schedule, || Ok(())).unwrap()
            .with_implicit_advection_entropy(thermal, || Ok(())).unwrap()
    };
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let system = bind(&thermal);
    let (stopped, report) = system.advance_to(&thermal, 1.0, controls(2.0), &mut workspace, || Ok(())).unwrap();
    near(stopped.cells()[&CellId(0)].specific_h_j_kg, 250.0, 0.0);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.entropy.advection.unwrap().upwind_mixing_j_k, 0.0, 0.0);
    // The stopped interval must not grant a later nonzero liquid flux through solid.
    assert!(system.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(())).is_err());

    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 100.0;
    let heated = bind(&thermal);
    let (next, report) = heated.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(())).unwrap();
    // First heat the stationary inventory from h=250 to 350 at zero flow.
    // Then solve  h-350 = 100 + 2*(600-h)  after crossing the liquidus.
    near(next.cells()[&CellId(0)].specific_h_j_kg, 550.0, 1e-9);
    near(report.energy.external_heat_j, 300.0, 1e-9);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
}

#[test]
fn scheduled_reconstructed_graph_uses_the_scaled_exact_reference_tangent() {
    let law = ShiftedLaw { base: Law { nonlinear_k: true, plateau: false }, datum: 0.0 };
    let (mut thermal, conduction) = setup(&law, &[250.0, 375.0, 600.0], None, vec![]);
    for (i, mass) in [1.0, 2.0, 3.0].into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    let flow = open_flow(3, 1.0, 500.0);
    let plan = axial_plan(&thermal);
    let schedule = history(3, &[(0.0, 0.5, 500.0), (1.0, 3.5, 700.0)]);
    let system = ImplicitConductionChain::prepare(
        &thermal, &conduction, &[(MaterialId(0), &law)],
        LinearTopology::Graph(ImplicitConductionLinearControls::default()), || Ok(()),
    ).unwrap().with_implicit_advection_reconstruction(&thermal, &flow, &plan, || Ok(())).unwrap()
        .with_implicit_advection_schedule(&thermal, &schedule, || Ok(())).unwrap();
    let step = BackwardEulerStep::new(&thermal, 0.5, || Ok(())).unwrap();
    let energy = [4.0, -3.0, 2.0];
    let mut report = empty_report();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let trial = system.evaluate(&step, &energy, &BTreeMap::new(), &[], false, controls(1.0),
        &mut report, &mut || Ok(())).unwrap();
    system.jacobian(&thermal, trial.payload(), 0.5, &mut workspace, &mut || Ok(())).unwrap();
    let coefficient = |row, col| if row == col { workspace.diagonal[row] } else {
        workspace.graph_entries.as_ref().unwrap().iter()
            .filter(|&&(r, c, _)| r == row && c == col).map(|entry| entry.2).sum()
    };
    for col in 0..3 {
        let mut plus = energy;
        let mut minus = energy;
        plus[col] += 1e-3;
        minus[col] -= 1e-3;
        let a = system.evaluate(&step, &plus, &BTreeMap::new(), &[], false, controls(1.0),
            &mut report, &mut || Ok(())).unwrap();
        let b = system.evaluate(&step, &minus, &BTreeMap::new(), &[], false, controls(1.0),
            &mut report, &mut || Ok(())).unwrap();
        for row in 0..3 {
            near(coefficient(row, col), (a.residual_j()[row] - b.residual_j()[row]) / 0.002, 1e-7);
        }
    }
}

#[test]
fn schedule_refusals_and_late_cancellation_leave_native_history_unchanged() {
    let knot = AdvectionScheduleKnot {
        time_s: 0.0, flow_multiplier: 1.0,
        inlet_enthalpies_j_kg: BTreeMap::from([(1, 400.0)]),
    };
    assert!(AdvectionSchedule::new(vec![knot.clone(), knot.clone()], 2, 1, || Ok(())).is_err());
    assert!(AdvectionSchedule::new(vec![knot.clone()], 0, 1, || Ok(())).is_err());
    assert!(AdvectionSchedule::new(vec![knot.clone()], 1, 0, || Ok(())).is_err());
    for multiplier in [-1.0, f64::NAN, f64::INFINITY] {
        let mut bad = knot.clone();
        bad.flow_multiplier = multiplier;
        assert!(AdvectionSchedule::new(vec![bad], 1, 1, || Ok(())).is_err());
    }
    let law = ShiftedLaw { base: Law { nonlinear_k: false, plateau: false }, datum: 0.0 };
    let (thermal, conduction) = setup(&law, &[300.0], None, vec![]);
    let flow = open_flow(1, 1.0, 400.0);
    for schedule in [history(2, &[(0.0, 1.0, 400.0)]), history(1, &[(0.0, 1.0, 1500.0)])] {
        assert!(ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap().with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).is_err());
    }
    let schedule = history(1, &[(0.0, 1.0, 400.0), (1.0, 2.0, 500.0)]);
    let system = ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
        .unwrap().with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut calls = 0;
    let (expected, report) = system.advance_sdirk2_to(&thermal, 1.0, controls(0.25), &mut workspace, || {
        calls += 1; Ok(())
    }).unwrap();
    let refuse = calls * 3 / 4;
    let mut calls = 0;
    assert!(matches!(system.advance_sdirk2_to(&thermal, 1.0, controls(0.25), &mut workspace, || {
        calls += 1;
        if calls == refuse { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))));
    let (repeated, repeated_report) = system.advance_sdirk2_to(
        &thermal, 1.0, controls(0.25), &mut workspace, || Ok(()),
    ).unwrap();
    assert_eq!(repeated.cells(), expected.cells());
    assert_eq!(repeated_report, report);
    let mut budget = controls(0.25);
    budget.maximum_residual_evaluations = 1;
    assert!(matches!(system.advance_sdirk2_to(&thermal, 1.0, budget, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(_)))));
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
    assert_eq!(flow.audit().substeps, 0);
}
