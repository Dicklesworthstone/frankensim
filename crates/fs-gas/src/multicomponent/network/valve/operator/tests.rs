//! Actual shared-tank transfers and transactional source publication.
use super::*;
use crate::multicomponent::network::dissipation::tests::{close, fields, gas};
use crate::multicomponent::{Boundaries, MixtureEuler1d};
use crate::multicomponent::network::valve::ValveDirection;
use fs_qty::{Area, Length};

fn fixture() -> (MixtureNetwork1d, NetworkValveOperator) {
    let g = gas(2, true, 5.0);
    let branch = MixtureEuler1d::new(g.clone(), Length::new(0.5), 0.1,
        vec![g.conserved(&[0.5, 0.5], 0.0, 1.5).unwrap(); 2],
        Boundaries::Periodic, || false).unwrap();
    let plenums = [([0.8, 0.2], 3.0, 0.2), ([0.2, 0.8], 1.0, 0.1), ([0.5, 0.5], 0.8, 0.15)]
        .map(|(d, p, v)| Plenum::new(v, g.conserved(&d, 0.0, p).unwrap()).unwrap()).to_vec();
    let n = MixtureNetwork1d::new(vec![branch], plenums, vec![], || false).unwrap();
    let orifice = Orifice::new(Area::new(0.01), 0.8, 1.0, ValveDirection::Bidirectional).unwrap();
    let op = NetworkValveOperator::new(&n, vec![
        ValveLink { first: 0, second: 1, orifice },
        ValveLink { first: 1, second: 2, orifice },
    ], || false).unwrap();
    (n, op)
}

#[test]
fn shared_tank_sweep_closes_actual_species_energy_and_keeps_branch_clock() {
    let (mut n, mut op) = fixture();
    let before = fields(&n);
    let old = n.plenums().to_vec();
    let time = n.time();
    let report = op.advance_for(&mut n, Time::new(0.02), ValveSolve::default(), || false).unwrap();
    assert_eq!(report.elapsed, Time::new(0.02));
    assert_eq!(report.sweeps, 1);
    assert!(report.iterations > 0);
    assert_eq!(op.last_report(), Some(report));
    assert_eq!(fields(&n).0, before.0);
    assert_eq!(n.time(), time);
    let exchanges = op.accepted_exchanges().unwrap();
    assert!(exchanges.iter().all(|x| x.mass_to_second_kg() > 0.0));
    let mut physical = 0.0;
    for (i, p) in n.plenums().iter().enumerate() {
        let mut species = [0.0; 2];
        let mut energy = 0.0;
        for (link, exchange) in op.links().iter().zip(exchanges) {
            let sign = if link.first == i { -1.0 } else if link.second == i { 1.0 } else { 0.0 };
            for k in 0..2 { species[k] += sign * exchange.species_to_second_kg[k]; }
            energy += sign * exchange.energy_to_second_j;
        }
        for k in 0..2 {
            close(p.volume_m3() * (p.state().partial_densities[k] - old[i].state().partial_densities[k]), species[k], 3e-14);
        }
        close(p.volume_m3() * (p.state().energy - old[i].state().energy), energy, 3e-13);
        physical += p.volume_m3() * (n.gas().absolute_energy_density(p.state()).unwrap()
            - n.gas().absolute_energy_density(old[i].state()).unwrap());
    }
    close(physical, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    assert!(n.last_report().is_none());
}

#[test]
fn split_receipts_sum_only_valve_intervals_not_intervening_external_heating() {
    let (mut n, mut op) = fixture();
    let before = n.totals();
    op.begin_trial(&mut || false).unwrap();
    let first = op.advance_for_trial(&mut n, Time::new(0.01), ValveSolve::default(), false, &mut || false).unwrap();
    assert!(op.last_report().is_none());
    assert!(op.accepted_exchanges().is_none());
    // An independently supplied gas-energy increment between valve halves.
    let heat = 0.03;
    n.plenums[1].state.energy += heat / n.plenums[1].volume_m3();
    let second = op.advance_for_trial(&mut n, Time::new(0.01), ValveSolve::default(), true, &mut || false).unwrap();
    assert!(op.last_report().is_none());
    op.publish_trial_report();
    let report = op.last_report().unwrap();
    assert_eq!(report.elapsed, Time::new(0.02));
    assert_eq!(report.sweeps, 2);
    assert_eq!(report.iterations, first.iterations + second.iterations);
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    close(n.totals().energy - before.energy, heat, 3e-12);
    assert!(op.accepted_exchanges().unwrap().iter().any(|x| x.mass_to_second_kg() > 0.0));
}

#[test]
fn source_cancellation_and_root_refusal_preserve_prior_reports_and_replay() {
    let make = || {
        let (mut n, mut op) = fixture();
        op.advance_for(&mut n, Time::new(0.01), ValveSolve::default(), || false).unwrap();
        (n, op)
    };
    let (mut reference, mut reference_op) = make();
    let mut polls = 0;
    reference_op.advance_for(&mut reference, Time::new(0.02), ValveSolve::default(), || { polls += 1; false }).unwrap();
    assert!(polls > 10);
    for stop in [1, polls / 2, polls] {
        let (mut n, mut op) = make();
        let before = fields(&n);
        let time = n.time();
        let receipt = op.last_report();
        let exchanges = op.accepted_exchanges().unwrap().to_vec();
        let mut calls = 0;
        assert_eq!(op.advance_for(&mut n, Time::new(0.02), ValveSolve::default(), || {
            calls += 1; calls >= stop
        }), Err(ValveError::Cancelled));
        assert_eq!(fields(&n), before);
        assert_eq!(n.time(), time);
        assert_eq!(op.last_report(), receipt);
        assert_eq!(op.accepted_exchanges().unwrap(), exchanges);
        let controls = ValveSolve { maximum_iterations: 1, ..ValveSolve::default() };
        assert!(op.advance_for(&mut n, Time::new(0.02), controls, || false).is_err());
        assert_eq!(fields(&n), before);
        assert_eq!(op.last_report(), receipt);
        op.advance_for(&mut n, Time::new(0.02), ValveSolve::default(), || false).unwrap();
        assert_eq!(fields(&n), fields(&reference));
        assert_eq!(op.accepted_exchanges(), reference_op.accepted_exchanges());
    }
}

#[test]
fn zero_interval_is_identity_and_wrong_volume_or_source_refuses_before_publication() {
    let (mut n, mut op) = fixture();
    let original = fields(&n);
    let report = op.advance_for(&mut n, Time::new(0.0), ValveSolve::default(), || false).unwrap();
    assert_eq!(fields(&n), original);
    assert_eq!(report.iterations, 0);
    assert_eq!(report.conservation_residual.mass(), 0.0);
    let accepted = op.last_report();
    n.plenums[1].volume_m3 *= 2.0;
    assert!(op.advance_for(&mut n, Time::new(0.01), ValveSolve::default(), || false).is_err());
    assert_eq!(fields(&n), original);
    assert_eq!(op.last_report(), accepted);
    n.plenums[1].volume_m3 *= 0.5;
    op.gas = gas(2, true, 6.0);
    assert!(op.advance_for(&mut n, Time::new(0.01), ValveSolve::default(), || false).is_err());
    assert_eq!(fields(&n), original);
    assert_eq!(op.last_report(), accepted);
}
