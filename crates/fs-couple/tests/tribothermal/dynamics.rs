//! Force-driven material contact, using the same source and finite heat fixtures.
use super::*;
use fs_couple::tribothermal::dynamics::*;
use fs_thermal::enthalpy::{DepositionBudget, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent};

fn motion_controls() -> MotionSolve {
    let mut sliding = controls();
    // Resolve the native compliant return map more tightly than the outer
    // source/impulse gate, so its bisection error cannot mask source feedback.
    sliding.tangential.max_iterations = 64;
    sliding.tangential.force_absolute_tolerance_n = 1e-12;
    sliding.tangential.force_relative_tolerance = 1e-12;
    MotionSolve {
        sliding, maximum_motion_trials: 80, maximum_newton_iterations: 20,
        difference_step_mps: 1e-5, impulse_tolerance_n_s: 1e-8,
        total_energy_tolerance_j: 1e-6,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-10, relative_tolerance: 1e-12,
            linear_restart: 1, ..NewtonKrylovConfig::default() },
    }
}

fn inventory<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2], mass: [f64; 2]) -> DepositionSimulation<'a> {
    let base = common::initial(curves, 0.0, None);
    let mut events: Vec<_> = CELLS.into_iter().zip(mass).map(|(id, mass_kg)| {
        let mut cell = base.cells()[&id].clone();
        cell.mass_kg = mass_kg;
        cell.geometry.volume_m3 = mass_kg; // Exact original sourced density is 1 kg/m³.
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(base.links()[&7].clone()) });
    DepositionSimulation::new(events,
        &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
          (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 20 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
}

#[test]
fn original_body_masses_and_force_work_close_momentum_heat_and_kinetic_storage() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    // This reference is Coulomb with equal static/kinetic capacity; distinct
    // stiction thresholds are exercised independently below.
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = inventory(&curves, [2.0, 3.0]);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.4, 0.6, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [1.0, 0.0]).unwrap();
    assert_eq!(pair.masses_kg(), [2.0, 3.0]);
    let next = pair.step_candidate(&thermal, 0.1,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [3.0, -1.0] },
        motion_controls(), 0.8, &[0.0; 2], other, || false).unwrap();
    let expected = [1.0 + 0.1 * (3.0 - 2.0) / 2.0, 0.1 * (-1.0 + 2.0) / 3.0];
    for (v, exact) in next.pair().velocities_mps().into_iter().zip(expected) { close(v, exact); }
    let r = next.report();
    let heat = 2.0 * 0.1 * (expected[0] - expected[1]);
    let kinetic = 0.5 * (2.0 * (expected[0] * expected[0] - 1.0) + 3.0 * expected[1] * expected[1]);
    let numerical = 0.5 * (2.0 * (expected[0] - 1.0).powi(2) + 3.0 * expected[1].powi(2));
    close(r.kinetic_energy_change_j, kinetic);
    close(r.numerical_dissipation_j, numerical);
    close(r.external_work_j, 0.1 * (3.0 * expected[0] - expected[1]));
    close(r.sliding.deposited_heat_j, heat);
    close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + 0.4 * heat / 4.0);
    close(next.thermal().temperature_k(CELLS[1]).unwrap(), 300.0 + 0.6 * heat / 9.0);
    close(r.total_energy_residual_j, 0.0);
    assert!(r.body_momentum_residual_n_s.iter().all(|r| r.abs() < 1e-8));
    assert!(r.motion_trials > 1 && r.friction_trials > r.sliding.friction_trials);
    // The scheme's kinetic defect is independently visible, and not body heat.
    assert!(r.numerical_dissipation_j > 1e-3);
    assert!((r.sliding.thermal.stored_energy_change_j - heat - numerical).abs() > 1e-3);
    close(next.pair().kinetic_energy_j().unwrap() - pair.kinetic_energy_j().unwrap(), kinetic);
    close(next.pair().contact().work().dissipated_work_j(), heat);
    assert_eq!(pair.velocities_mps(), [1.0, 0.0]);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
}

#[test]
fn heat_capacity_data_changes_actual_slowing_through_endpoint_material_friction() {
    let mut speeds = Vec::new();
    for cp_a in [0.1, 1.0] {
        let (a, ca) = common::material(1.0, cp_a, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b]; let curves = [ca, cb]; let point = context();
        // mu(Ta) = 0.2 + 0.1*(Ta-300), explicitly within its 299..310 K range.
        // Static capacity 4 N is below the 5 N reaction needed to arrest this
        // interval, and above the reached temperature-dependent kinetic force.
        let card = interface_with_static(&cards, 0.4, "Ta", Temperature::DIMS, [(299.0, 0.1), (310.0, 1.2)], 2000.0);
        let thermal = common::initial(&curves, 0.0, None);
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::new(&thermal, contact, [1.0, 0.0]).unwrap();
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [0.0; 2] },
            motion_controls(), 0.6, &[0.0; 2], other, || false).unwrap();
        // Independent elimination: mu=.2/(1-k*s), s=1-2*mu, k=.05/Ca.
        // Select the positive branch continuous from k=0 without cancellation.
        let k = 0.05 / cp_a;
        let exact = 1.2 / (1.0 + k + ((1.0 + k).powi(2) - 2.4 * k).sqrt());
        let mu = (1.0 - exact) / 2.0;
        let velocity = next.pair().velocities_mps();
        let speed = velocity[0] - velocity[1];
        close(speed, exact);
        close(velocity[0], 0.5 * (1.0 + exact));
        close(velocity[1], 0.5 * (1.0 - exact));
        close(next.report().sliding.material.speed_mps, speed);
        close(next.report().sliding.friction.traction_n[0], 10.0 * mu);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + 0.5 * mu * exact / cp_a);
        close(next.pair().contact().accumulated_distance_m(), 0.1 * exact);
        close(next.report().sliding.deposited_heat_j, mu * exact);
        close(next.report().numerical_dissipation_j, mu * mu);
        close(next.report().total_energy_residual_j, 0.0);
        assert_eq!(pair.contact().binding().selection, next.pair().contact().binding().selection);
        speeds.push(speed);
    }
    assert!(speeds[1] - speeds[0] > 0.1);
}

#[test]
fn late_failures_and_budget_refuse_before_any_motion_or_heat_is_published() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let original_cells = thermal.cells().clone();
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [1.0, 0.0]).unwrap();
    let drive = MotionDrive { normal_load_n: 10.0, body_forces_n: [0.0; 2] };
    let polls = Cell::new(0);
    let good = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), 0.8, &[0.0; 2], other,
        || { polls.set(polls.get() + 1); false }).unwrap();
    let last = polls.get(); polls.set(0);
    assert!(matches!(pair.step_candidate(&thermal, 0.1, drive, motion_controls(), 0.8, &[0.0; 2], other,
        || { polls.set(polls.get() + 1); polls.get() == last }), Err(SlidingError::Cancelled)));
    let mut limited = motion_controls(); limited.maximum_motion_trials = good.report().motion_trials - 1;
    assert!(matches!(pair.step_candidate(&thermal, 0.1, drive, limited, 0.8, &[0.0; 2], other, || false),
        Err(SlidingError::MotionBudget)));
    let refused = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), 0.8, &[0.0; 2], |state| {
        if state.temperature_k(CELLS[0]).unwrap() > 300.01 { return Err("thermal transport refused"); }
        Ok(StepCoefficients { sources_w: BTreeMap::new(), conductances_w_k: BTreeMap::new(), payload: () })
    }, || false);
    assert!(matches!(refused, Err(SlidingError::Source("thermal transport refused"))));
    let changed_inventory = inventory(&curves, [2.0, 1.0]);
    assert!(matches!(pair.step_candidate(&changed_inventory, 0.1, drive, motion_controls(), 0.8, &[0.0; 2], other, || false),
        Err(SlidingError::Invalid(_))));
    assert_eq!(pair.velocities_mps(), [1.0, 0.0]);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
    assert_eq!(pair.contact().work().dissipated_work_j(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    let retry = pair.step_candidate(&thermal, 0.1, drive, motion_controls(), 0.8, &[0.0; 2], other, || false).unwrap();
    assert_eq!(retry.thermal().cells(), good.thermal().cells());
    assert_eq!(retry.pair().velocities_mps(), good.pair().velocities_mps());
    assert_eq!(retry.pair().contact().state(), good.pair().contact().state());
    let (thermal, pair, _) = retry.into_parts();
    let continued = pair.step_candidate(&thermal, 0.2, drive, motion_controls(), 0.4,
        &[0.0; 2], other, || false).unwrap();
    close(continued.pair().velocities_mps()[0], 0.6);
    close(continued.pair().velocities_mps()[1], 0.4);
    close(continued.pair().contact().accumulated_distance_m(), 0.08);
    close(continued.pair().contact().work().dissipated_work_j(), 0.16);
    close(continued.report().total_energy_residual_j, 0.0);
}

#[test]
fn native_static_capacity_drives_rest_below_at_and_above_breakaway() {
    use fs_tribo::FrictionRegime;
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.4, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 0.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [0.0; 2]).unwrap();
    for force in [6.0, 8.0, 10.0] {
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [force, 0.0] },
            motion_controls(), 0.3, &[0.0; 2], other, || false).unwrap();
        let r = next.report();
        let v = next.pair().velocities_mps();
        if force <= 8.0 {
            // Native static traction F/2 gives the common acceleration F/(mA+mB).
            assert_eq!(r.sliding.friction.regime, FrictionRegime::Sticking);
            assert_eq!(v[0], v[1]);
            close(v[0], 0.05 * force);
            close(r.sliding.friction.traction_n[0], force / 2.0);
            assert_eq!(r.sliding.friction.slip_increment_m, [0.0; 2]);
            assert_eq!(r.sliding.material.speed_mps, 0.0);
            assert_eq!(r.sliding.deposited_heat_j, 0.0);
            assert_eq!(next.pair().contact().accumulated_distance_m(), 0.0);
            assert_eq!(next.thermal().cells(), thermal.cells());
        } else {
            assert_eq!(r.sliding.friction.regime, FrictionRegime::Sliding);
            close(v[0], 0.8); close(v[1], 0.2);
            close(r.sliding.friction.traction_n[0], 2.0);
            close(r.sliding.friction.slip_increment_m[0], 0.06);
            close(r.sliding.deposited_heat_j, 0.12);
        }
        assert!(r.body_momentum_residual_n_s.iter().all(|x| x.abs() < 1e-8));
        assert!(r.source_motion_residual_n_s.abs() < 1e-8);
        close(r.total_energy_residual_j, 0.0);
    }
}

#[test]
fn signed_source_speed_closes_against_the_native_solved_slip() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.4, "speed", SPEED, [(0.0, 0.1), (2.0, 0.5)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    for sign in [1.0, -1.0] {
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::new(&thermal, contact, [sign, 0.0]).unwrap();
        let next = pair.step_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [0.0; 2] },
            motion_controls(), sign * 0.5, &[0.0; 2], other, || false).unwrap();
        // s=1-2*(.1+.2*s), so |s|=4/7; querying released speed would differ.
        let r = next.report();
        close(next.pair().velocities_mps()[0] - next.pair().velocities_mps()[1], sign * 4.0 / 7.0);
        close(r.sliding.material.speed_mps, 4.0 / 7.0);
        close(r.sliding.friction.traction_n[0], sign * 15.0 / 7.0);
        close(r.sliding.material.accumulated_distance_m, 0.1 * 4.0 / 7.0);
        close(r.sliding.deposited_heat_j, 6.0 / 49.0);
        close(r.total_energy_residual_j, 0.0);
        assert!(r.source_motion_residual_n_s.abs() < 1e-8);
    }
}

#[test]
fn exact_native_stop_then_opposite_restart_retains_heat_and_signed_history() {
    use fs_tribo::FrictionRegime;
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [1.0, 0.0]).unwrap();
    let stopped = pair.step_candidate(&thermal, 0.25,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [0.0; 2] },
        motion_controls(), 0.2, &[0.0; 2], other, || false).unwrap();
    assert_eq!(stopped.report().sliding.friction.regime, FrictionRegime::Sticking);
    assert_eq!(stopped.pair().velocities_mps(), [0.5; 2]);
    assert_eq!(stopped.report().sliding.material.speed_mps, 0.0);
    assert_eq!(stopped.report().sliding.deposited_heat_j, 0.0);
    close(stopped.report().sliding.friction.traction_n[0], 2.0);
    close(stopped.report().kinetic_energy_change_j, -0.25);
    close(stopped.report().numerical_dissipation_j, 0.25);
    close(stopped.report().total_energy_residual_j, 0.0);
    // Endpoint stopping is BE numerical loss, not invented frictional heating.
    assert_eq!(stopped.thermal().cells(), thermal.cells());
    let (thermal, pair, _) = stopped.into_parts();
    let reversed = pair.step_candidate(&thermal, 0.35,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [-10.0, 0.0] },
        motion_controls(), -0.5, &[0.0; 2], other, || false).unwrap();
    let r = reversed.report();
    assert_eq!(r.sliding.friction.regime, FrictionRegime::Sliding);
    close(reversed.pair().velocities_mps()[0], -0.3);
    close(reversed.pair().velocities_mps()[1], 0.3);
    close(r.sliding.friction.traction_n[0], -2.0);
    close(r.sliding.friction.slip_increment_m[0], -0.06);
    close(r.sliding.material.speed_mps, 0.6);
    close(r.sliding.material.accumulated_distance_m, 0.06);
    close(reversed.pair().contact().state().slip_m()[0], -0.06);
    close(reversed.pair().contact().work().dissipated_work_j(), 0.12);
    close(r.sliding.deposited_heat_j, 0.12);
    close(r.total_energy_residual_j, 0.0);
    // A zero-speed source limit is actually queried; no nonzero speed substitute.
    let unsupported = interface_with_static(&cards, 0.2, "speed", SPEED,
        [(0.1, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &unsupported, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 0.5, InputAuthority::SyntheticFixture).unwrap();
    let bounded = SlidingPair::new(&thermal, contact, [0.5; 2]).unwrap();
    assert!(matches!(bounded.step_candidate(&thermal, 0.35,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [0.0; 2] },
        motion_controls(), 0.0, &[0.0; 2], other, || false), Err(SlidingError::Material(_))));
    assert_eq!(bounded.velocities_mps(), [0.5; 2]);
    assert_eq!(bounded.contact().work().dissipated_work_j(), 0.0);
}
