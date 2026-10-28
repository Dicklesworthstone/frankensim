//! Physical oracles through the actual source-bound contact and finite heat solve.
use super::*;
use fs_couple::tribothermal::dynamics::collinear::CollinearMotionSolve;

fn event_controls(axis: usize) -> CollinearMotionSolve {
    CollinearMotionSolve { sliding: motion_controls().sliding, tangent_axis: axis,
        maximum_segments: 3, maximum_motion_trials: 200, maximum_bisections: 80,
        impulse_tolerance_n_s: 1e-10, event_time_tolerance_s: 1e-10,
        total_energy_tolerance_j: 1e-6 }
}

#[test]
fn sourced_bodies_stop_at_the_physical_time_and_all_relative_kinetic_energy_heats_them() {
    let (a, ca, body_a) = sourced_body(2.0, 2.0);
    let (b, cb, body_b) = sourced_body(3.0, 3.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let thermal = inventory_with_geometry(&curves, [2.0, 3.0], [1.0, 1.0]);
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.4, 0.6, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::from_material_bodies(&thermal, contact, [body_a, body_b],
        [[1.0, 0.0], [0.0; 2]]).unwrap();
    let result = pair.step_collinear_candidate(&thermal, 1.0,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
        event_controls(0), &[0.0; 2], other, || false).unwrap();
    let report = result.report();
    assert_eq!(report.segments.len(), 2);
    assert!(report.segments[0].reached_stop);
    close(report.segments[0].end_time_s, 0.6); // reduced mass 1.2 kg / 2 N.
    assert_eq!(report.segments[1].sliding.friction.regime, fs_tribo::FrictionRegime::Sticking);
    close(result.pair().velocities_mps()[0][0], 0.4);
    close(result.pair().velocities_mps()[1][0], 0.4);
    close(report.kinetic_energy_change_j, -0.6);
    close(report.deposited_heat_j, 0.6);
    close(report.thermal_energy_change_j, 0.6);
    close(report.external_work_j, 0.0);
    close(report.total_energy_residual_j, 0.0);
    close(result.pair().contact().accumulated_distance_m(), 0.3);
    close(result.pair().contact().work().dissipated_work_j(), 0.6);
    close(result.thermal().temperature_k(CELLS[0]).unwrap(), 300.06);
    close(result.thermal().temperature_k(CELLS[1]).unwrap(), 300.04);
    assert!(result.pair().material_bodies().is_some());
    assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
    assert_eq!(thermal.temperature_k(CELLS[0]).unwrap(), 300.0);
}

#[test]
fn one_call_stops_then_reverses_with_exact_force_work_and_unsigned_sliding_heat() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (4.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    let next = pair.step_collinear_candidate(&thermal, 0.5,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[-6.0, 0.0], [0.0; 2]] },
        event_controls(0), &[0.0; 2], other, || false).unwrap();
    let r = next.report();
    assert_eq!(r.segments.len(), 2);
    close(r.segments[0].end_time_s, 0.1);
    assert!(r.segments[0].reached_stop);
    assert!(r.segments[1].sliding.friction.traction_n[0] < 0.0);
    close(next.pair().velocities_mps()[0][0], -1.4);
    close(next.pair().velocities_mps()[1][0], -0.6);
    close(r.body_displacements_m[0][0], -0.18);
    close(r.body_displacements_m[1][0], -0.07);
    close(r.external_work_j, 1.08);
    close(r.kinetic_energy_change_j, 0.66);
    close(r.deposited_heat_j, 0.42);
    close(r.total_energy_residual_j, 0.0);
    close(next.pair().contact().state().slip_m()[0], -0.11);
    close(next.pair().contact().accumulated_distance_m(), 0.21);
    close(next.pair().contact().work().dissipated_work_j(), 0.42);
}

#[test]
fn finite_heat_capacity_changes_midpoint_friction_and_actual_deceleration() {
    let mut final_speeds = Vec::new();
    for capacity_a in [0.1, 1.0] {
        let (a, ca) = common::material(1.0, capacity_a, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b]; let curves = [ca, cb]; let point = context();
        let thermal = inventory(&curves, [1.0, 1.0]);
        let card = interface_with_static(&cards, 0.6, "Ta", Temperature::DIMS,
            [(299.0, 0.1), (310.0, 1.2)], 2000.0);
        let contact = MaterialSlidingContact::new(&thermal,
            binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
            10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
        let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
        let next = pair.step_collinear_candidate(&thermal, 0.1,
            MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] },
            event_controls(0), &[0.0; 2], other, || false).unwrap();
        // Source law mu=.2+.1*(Ta-300), heat=mu*(1-mu), Ta-300=heat/(2*Ca).
        let k = 0.05 / capacity_a;
        let mu = 0.4 / ((1.0-k) + ((1.0-k)*(1.0-k) + 0.8*k).sqrt());
        let heat = mu * (1.0-mu);
        let relative = 1.0 - 2.0*mu;
        let velocities = next.pair().velocities_mps();
        close(velocities[0][0]-velocities[1][0], relative);
        close(next.report().deposited_heat_j, heat);
        close(next.report().kinetic_energy_change_j, -heat);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.0 + heat/(2.0*capacity_a));
        close(next.report().segments[0].sliding.material.speed_mps, 1.0-mu);
        close(next.report().segments[0].sliding.friction.traction_n[0], 10.0*mu);
        close(next.report().total_energy_residual_j, 0.0);
        assert_eq!(pair.contact().binding().selection, next.pair().contact().binding().selection);
        final_speeds.push(relative);
    }
    assert!(final_speeds[1] > final_speeds[0] + 0.1);
}

#[test]
fn common_transverse_acceleration_and_static_relative_motion_have_correct_work() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let thermal = inventory(&curves, [2.0, 3.0]);
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 0.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [1.0, 0.0]]).unwrap();
    let next = pair.step_collinear_candidate(&thermal, 0.2,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[4.0, 1.0], [6.0, 0.0]] },
        event_controls(1), &[0.0; 2], other, || false).unwrap();
    assert_eq!(next.report().segments.len(), 1);
    assert_eq!(next.report().segments[0].sliding.friction.regime, fs_tribo::FrictionRegime::Sticking);
    for velocity in next.pair().velocities_mps() { close(velocity[0], 1.4); close(velocity[1], 0.04); }
    close(next.report().deposited_heat_j, 0.0);
    close(next.report().external_work_j, 2.404);
    close(next.report().kinetic_energy_change_j, 2.404);
    close(next.report().total_energy_residual_j, 0.0);
    assert!(pair.step_collinear_candidate(&thermal, 0.2,
        MotionDrive { normal_load_n: 10.0, body_forces_n: [[4.0, 1.0], [5.0, 0.0]] },
        event_controls(1), &[0.0; 2], other, || false).is_err());
}

#[test]
fn refusal_after_a_solved_stop_preserves_all_original_histories_and_replay() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(&cards, 0.2, "speed", SPEED, [(0.0, 0.2), (2.0, 0.2)], 2000.0);
    let contact = MaterialSlidingContact::new(&thermal,
        binding(&cards, &curves, &card, &point, HeatPartition::new(0.5, 0.5, 0.0).unwrap()),
        10.0, 1.0, InputAuthority::SyntheticFixture).unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[1.0, 0.0], [0.0; 2]]).unwrap();
    let drive = MotionDrive { normal_load_n: 10.0, body_forces_n: [[0.0; 2]; 2] };
    let old_cells = thermal.cells().clone();
    let mut limited = event_controls(0); limited.maximum_segments = 1;
    assert!(matches!(pair.step_collinear_candidate(&thermal, 1.0, drive, limited,
        &[0.0; 2], other, || false), Err(SlidingError::MotionBudget)));
    let polls = Cell::new(0usize);
    let good = pair.step_collinear_candidate(&thermal, 1.0, drive, event_controls(0),
        &[0.0; 2], other, || { polls.set(polls.get()+1); false }).unwrap();
    let last = polls.get(); polls.set(0);
    assert!(matches!(pair.step_collinear_candidate(&thermal, 1.0, drive, event_controls(0),
        &[0.0; 2], other, || { polls.set(polls.get()+1); polls.get() == last }), Err(SlidingError::Cancelled)));
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(pair.contact().time_s(), 0.0);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
    assert_eq!(pair.contact().work().dissipated_work_j(), 0.0);
    assert_eq!(pair.velocities_mps(), [[1.0, 0.0], [0.0; 2]]);
    let retry = pair.step_collinear_candidate(&thermal, 1.0, drive, event_controls(0),
        &[0.0; 2], other, || false).unwrap();
    assert_eq!(retry.pair().velocities_mps(), good.pair().velocities_mps());
    assert_eq!(retry.pair().contact().state(), good.pair().contact().state());
    assert_eq!(retry.thermal().cells(), good.thermal().cells());
}
