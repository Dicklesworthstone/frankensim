//! Public-API physics and rollback checks for valves + molecular flow + walls + chemistry.
#![cfg(feature = "experimental-euler")]
use fs_gas::chemistry::IsochoricKineticsSolve;
use fs_gas::multicomponent::{Boundaries, Boundary, CaloricDomain, Cell, Conserved, MixtureEuler1d, MixtureGas, StepControl};
use fs_gas::multicomponent::dissipation::{ConstantTransport, DissipationControl};
use fs_gas::multicomponent::network::{BranchEnd, MixtureNetwork1d, Plenum, PlenumConnection};
use fs_gas::multicomponent::network::dissipation::NetworkDissipationOperator;
use fs_gas::multicomponent::network::dissipating::{DissipatingNetwork1d, DissipatingStepControl};
use fs_gas::multicomponent::network::reactive::{NetworkSources, ReactingNetwork1d, ReactiveNetworkError, ReactiveNetworkReport, ReactiveStepControl};
use fs_gas::multicomponent::network::valve::{NetworkValveOperator, Orifice, ValveDirection, ValveLink, ValveSolve};
use fs_gas::multicomponent::network::wall_exchange::NetworkWallOperator;
use fs_gas::multicomponent::wall_exchange::{WallContact, WallControl, WallDrag, WallThermal};
use fs_qty::{Area, Length, Temperature, Time};
use fs_thermochem::kinetics::IdealGasKineticsV1;

type Diffusion = NetworkDissipationOperator<ConstantTransport>;
#[path = "valved_sources/fixture.rs"]
mod fixture;
use fixture::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "actual={a:.16e}, expected={b:.16e}, tolerance={tolerance:e}");
}
fn owned(q: Cell<'_>) -> Conserved {
    Conserved { partial_densities: q.partial_densities.to_vec(), momentum: q.momentum, energy: q.energy }
}
fn fields(branches: &[MixtureEuler1d], plenums: &[Plenum]) -> Vec<Conserved> {
    branches.iter().flat_map(|b| (0..b.cell_count()).map(move |i| owned(b.cell(i).unwrap())))
        .chain(plenums.iter().map(|p| owned(p.state()))).collect()
}
fn assert_fields(actual: &[Conserved], expected: &[Conserved], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        for k in 0..2 { close(a.partial_densities[k], b.partial_densities[k], tolerance); }
        close(a.momentum, b.momentum, tolerance);
        close(a.energy, b.energy, 10.0 * tolerance);
    }
}
// Direct integration of physical fields and actual wall storage. Neither a
// network-total helper nor a reported conservation residual enters this sum.
fn inventory(g: &MixtureGas, branches: &[MixtureEuler1d], plenums: &[Plenum], wall: &NetworkWallOperator) -> ([f64; 2], f64) {
    let mut masses = [0.0; 2];
    let mut energy = 0.0;
    for b in branches {
        for (i, volume) in b.mesh().volumes().iter().enumerate() {
            let q = b.cell(i).unwrap();
            for k in 0..2 { masses[k] += volume * q.partial_densities[k]; }
            energy += volume * g.absolute_energy_density(q).unwrap();
        }
    }
    for p in plenums {
        for k in 0..2 { masses[k] += p.volume_m3() * p.state().partial_densities[k]; }
        energy += p.volume_m3() * g.absolute_energy_density(p.state()).unwrap();
    }
    for contacts in (0..branches.len()).map(|i| wall.branch_contacts(i).unwrap())
        .chain(core::iter::once(wall.plenum_contacts()))
    {
        for c in contacts { energy += wall_energy(*c); }
    }
    (masses, energy)
}
fn wall_energy(c: WallContact) -> f64 {
    match c.thermal {
        WallThermal::Lumped { temperature_k, heat_capacity_j_per_k, .. } => heat_capacity_j_per_k * temperature_k,
        _ => 0.0,
    }
}
fn wall_states(wall: &NetworkWallOperator, branches: usize) -> Vec<WallContact> {
    (0..branches).flat_map(|i| wall.branch_contacts(i).unwrap().iter().copied())
        .chain(wall.plenum_contacts().iter().copied()).collect()
}
fn advance(
    f: &mut ReactingNetwork1d, d: &mut Diffusion, w: &mut NetworkWallOperator,
    v: &mut NetworkValveOperator, dt: f64, cancelled: impl FnMut() -> bool,
) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
    f.step_with_dissipation_and_sources(d, Time::new(dt), ReactiveStepControl::default(),
        DissipationControl::default(), NetworkSources {
            wall: Some((w, WallControl::default())), valves: Some((v, ValveSolve::default())),
        }, cancelled)
}

#[test]
fn molecular_waves_valves_and_finite_walls_close_every_plenum_and_species_account() {
    let m = mechanism(0.0, 0.0);
    let n = raw(&m);
    let (mut wall, mut valves) = sources(&n);
    let law = law(n.gas(), 1.0);
    let mut f = DissipatingNetwork1d::new(n, law, vec![Length::new(0.1); 2], || false).unwrap();
    let initial = inventory(f.gas(), f.branches(), f.plenums(), &wall);
    let mut molecular_exchange = 0.0_f64;
    let mut wave_exchange = 0.0_f64;
    for _ in 0..6 {
        let old: Vec<_> = f.plenums().iter().map(|p| owned(p.state())).collect();
        let old_walls = wall.plenum_contacts().to_vec();
        let report = f.step_with_sources(Time::new(0.002), DissipatingStepControl::default(),
            NetworkSources { wall: Some((&mut wall, WallControl::default())),
                valves: Some((&mut valves, ValveSolve::default())) }, || false).unwrap();
        let v = valves.accepted_exchanges().unwrap()[0];
        assert!(v.mass_to_second_kg() > 0.0);
        let halves = report.valves.unwrap();
        assert_eq!(halves.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(halves.second.elapsed, halves.first.elapsed);
        assert_eq!(valves.last_report().unwrap().elapsed, report.elapsed);
        assert_eq!(valves.last_report().unwrap().sweeps, 2);
        assert_eq!(wall.last_report().unwrap().elapsed, report.elapsed);
        for i in 0..2 {
            let sign = if i == 0 { -1.0 } else { 1.0 };
            let p = &f.plenums()[i];
            let a = f.accepted_port_exchanges().unwrap()[i];
            let d = f.accepted_molecular_port_exchanges().unwrap()[i];
            for k in 0..2 {
                close(p.volume_m3() * (p.state().partial_densities[k] - old[i].partial_densities[k]),
                    a.species_masses_to_plenum[k] + d.species_masses_to_plenum[k] + sign * v.species_to_second_kg[k], 3e-12);
            }
            let wall_heat = wall_energy(old_walls[i]) - wall_energy(wall.plenum_contacts()[i]);
            close(p.volume_m3() * (p.state().energy - old[i].energy),
                a.energy_to_plenum + d.energy_to_plenum + sign * v.energy_to_second_j + wall_heat, 3e-10);
            molecular_exchange = molecular_exchange.max(d.energy_to_plenum.abs());
            wave_exchange = wave_exchange.max(a.energy_to_plenum.abs());
        }
        let current = inventory(f.gas(), f.branches(), f.plenums(), &wall);
        for k in 0..2 { close(current.0[k], initial.0[k], 3e-11); }
        close(current.1, initial.1, 3e-8);
        close(report.absolute_energy_residual_j, 0.0, 3e-9);
        assert_eq!(report.boundary_change.species_masses, [0.0; fs_gas::multicomponent::MAX_SPECIES]);
        close(report.boundary_change.energy, wall.last_report().unwrap().heat_to_gas_j, 3e-12);
        for b in f.branches() { assert_eq!(b.time(), f.time()); }
    }
    assert!(molecular_exchange > 0.0 && wave_exchange > 0.0);
}

#[test]
fn chemical_formation_energy_and_valve_enthalpy_are_each_counted_once() {
    let (mut f, mut d, mut w, mut v) = reacting(4.0, 10.0);
    let initial = inventory(f.gas(), f.branches(), f.plenums(), &w);
    for _ in 0..4 {
        let old: Vec<_> = f.plenums().iter().map(|p| owned(p.state())).collect();
        let old_walls = w.plenum_contacts().to_vec();
        let report = advance(&mut f, &mut d, &mut w, &mut v, 0.002, || false).unwrap();
        assert_eq!(report.chemistry.cell_solves, 20);
        assert!(report.chemistry.species_change_kg[1] > 0.0);
        let valve = v.accepted_exchanges().unwrap()[0];
        assert!(valve.mass_to_second_kg() > 0.0);
        let physical_valve = valve.energy_to_second_j + (0..2)
            .map(|k| valve.species_to_second_kg[k] * f.gas().energy_references()[k]).sum::<f64>();
        for i in 0..2 {
            let p = &f.plenums()[i];
            let a = f.accepted_port_exchanges().unwrap()[i];
            let d = f.accepted_molecular_port_exchanges().unwrap()[i];
            let sign = if i == 0 { -1.0 } else { 1.0 };
            let chemical_b = p.volume_m3() * (p.state().partial_densities[1] - old[i].partial_densities[1])
                - a.species_masses_to_plenum[1] - d.species_masses_to_plenum[1] - sign * valve.species_to_second_kg[1];
            assert!(chemical_b > 0.0, "plenum {i} did not react");
            let actual = p.volume_m3() * (f.gas().absolute_energy_density(p.state()).unwrap()
                - f.gas().absolute_energy_density(old[i].as_cell()).unwrap());
            let wall_heat = wall_energy(old_walls[i]) - wall_energy(w.plenum_contacts()[i]);
            close(actual, a.absolute_energy_to_plenum(f.gas()).unwrap()
                + d.absolute_energy_to_plenum(f.gas()).unwrap() + sign * physical_valve + wall_heat, 3e-9);
        }
        for k in 0..2 {
            close(report.after.species_masses[k] - report.before.species_masses[k], report.chemistry.species_change_kg[k], 3e-10);
        }
        let current = inventory(f.gas(), f.branches(), f.plenums(), &w);
        close(current.0.iter().sum(), initial.0.iter().sum(), 3e-11);
        close(current.1, initial.1, 3e-8);
        close(report.absolute_energy_residual_j, 0.0, 3e-9);
        close(report.boundary_change.energy, w.last_report().unwrap().heat_to_gas_j, 3e-12);
        assert!(report.valves.is_some() && report.wall.is_some() && report.dissipation.is_some());
    }
}

#[test]
fn inert_reacting_owner_matches_nonreacting_full_source_composition() {
    let (mut reactive, mut d, mut rw, mut rv) = reacting(0.0, 0.0);
    let m = mechanism(0.0, 0.0);
    let n = raw(&m);
    let (mut wall, mut valves) = sources(&n);
    let law = law(n.gas(), 1.0);
    let mut f = DissipatingNetwork1d::new(n, law, vec![Length::new(0.1); 2], || false).unwrap();
    for _ in 0..3 {
        let report = advance(&mut reactive, &mut d, &mut rw, &mut rv, 0.002, || false).unwrap();
        let reference = f.step_with_sources(report.elapsed, DissipatingStepControl::default(),
            NetworkSources { wall: Some((&mut wall, WallControl::default())),
                valves: Some((&mut valves, ValveSolve::default())) }, || false).unwrap();
        assert_eq!(report.elapsed, reference.elapsed);
        assert_fields(&fields(reactive.branches(), reactive.plenums()), &fields(f.branches(), f.plenums()), 3e-9);
        for (a, b) in wall_states(&rw, 2).into_iter().zip(wall_states(&wall, 2)) {
            close(wall_energy(a), wall_energy(b), 3e-9);
        }
    }
}

#[test]
fn heat_driven_restart_replays_every_source_at_the_actual_accepted_duration() {
    let (mut f, mut d, mut w, mut v) = reacting(40.0, 500.0);
    let report = advance(&mut f, &mut d, &mut w, &mut v, 0.2, || false).unwrap();
    assert!(report.coupling_retries > 0);
    assert_eq!(report.elapsed, report.transport.elapsed);
    assert_eq!(report.chemistry.cell_solves, 20);
    assert_eq!(v.last_report().unwrap().sweeps, 2);
    assert_eq!(v.last_report().unwrap().elapsed, report.elapsed);
    assert_eq!(w.last_report().unwrap().elapsed, report.elapsed);
    let half = 0.5 * report.elapsed.value();
    for duration in [report.valves.unwrap().first.elapsed, report.valves.unwrap().second.elapsed,
        report.wall.unwrap().first.elapsed, report.wall.unwrap().second.elapsed,
        report.dissipation.unwrap().first.elapsed, report.dissipation.unwrap().second.elapsed]
    { assert_eq!(duration.value(), half); }
    let (mut reference, mut rd, mut rw, mut rv) = reacting(40.0, 500.0);
    let accepted = advance(&mut reference, &mut rd, &mut rw, &mut rv, report.elapsed.value(), || false).unwrap();
    assert_eq!(accepted.elapsed, report.elapsed);
    assert_fields(&fields(f.branches(), f.plenums()), &fields(reference.branches(), reference.plenums()), 3e-6);
    let a = v.accepted_exchanges().unwrap()[0];
    let b = rv.accepted_exchanges().unwrap()[0];
    for k in 0..2 { close(a.species_to_second_kg[k], b.species_to_second_kg[k], 3e-8); }
    close(a.energy_to_second_j, b.energy_to_second_j, 3e-7);
    close(report.absolute_energy_residual_j, 0.0, 3e-8);
}

#[test]
fn late_full_composition_cancellation_and_root_exhaustion_publish_no_owner_prefix() {
    let make = || {
        let (mut f, mut d, mut w, mut v) = reacting(4.0, 10.0);
        advance(&mut f, &mut d, &mut w, &mut v, 0.002, || false).unwrap();
        (f, d, w, v)
    };
    let (mut reference, mut rd, mut rw, mut rv) = make();
    let mut polls = 0;
    advance(&mut reference, &mut rd, &mut rw, &mut rv, 0.002, || { polls += 1; false }).unwrap();
    assert!(polls > 100);
    for stop in [polls / 2, polls] {
        let (mut f, mut d, mut w, mut v) = make();
        let original = fields(f.branches(), f.plenums());
        let clock = f.time();
        let report = f.last_report();
        let ports = f.accepted_port_exchanges().unwrap().to_vec();
        let molecular_ports = f.accepted_molecular_port_exchanges().unwrap().to_vec();
        let branch_reports = f.accepted_branch_transport_reports().unwrap().to_vec();
        let wall = wall_states(&w, 2);
        let wall_report = w.last_report();
        let impulses = w.accepted_branch_impulses().unwrap().to_vec();
        let valve_report = v.last_report();
        let exchanges = v.accepted_exchanges().unwrap().to_vec();
        let diffusion_report = d.last_report();
        let diffusion_ports = d.accepted_port_exchanges().unwrap().to_vec();
        let mut calls = 0;
        assert_eq!(advance(&mut f, &mut d, &mut w, &mut v, 0.002, || { calls += 1; calls >= stop }), Err(ReactiveNetworkError::Cancelled));
        assert_eq!(fields(f.branches(), f.plenums()), original);
        assert_eq!(f.time(), clock);
        assert_eq!(f.last_report(), report);
        assert_eq!(f.accepted_port_exchanges().unwrap(), ports);
        assert_eq!(f.accepted_molecular_port_exchanges().unwrap(), molecular_ports);
        assert_eq!(f.accepted_branch_transport_reports().unwrap(), branch_reports);
        assert_eq!(wall_states(&w, 2), wall);
        assert_eq!(w.last_report(), wall_report);
        assert_eq!(w.accepted_branch_impulses().unwrap(), impulses);
        assert_eq!(v.last_report(), valve_report);
        assert_eq!(v.accepted_exchanges().unwrap(), exchanges);
        assert_eq!(d.last_report(), diffusion_report);
        assert_eq!(d.accepted_port_exchanges().unwrap(), diffusion_ports);
        let refused = f.step_with_dissipation_and_sources(&mut d, Time::new(0.002),
            ReactiveStepControl { maximum_coupling_retries: 0, ..ReactiveStepControl::default() },
            DissipationControl::default(), NetworkSources {
                wall: Some((&mut w, WallControl::default())),
                valves: Some((&mut v, ValveSolve { maximum_iterations: 1,
                    absolute_mass_tolerance_kg: 1e-30, relative_mass_tolerance: 1e-13 })),
            }, || false);
        assert!(refused.is_err());
        assert_eq!(fields(f.branches(), f.plenums()), original);
        assert_eq!(wall_states(&w, 2), wall);
        assert_eq!(v.accepted_exchanges().unwrap(), exchanges);
        assert_eq!(d.last_report(), diffusion_report);
        advance(&mut f, &mut d, &mut w, &mut v, 0.002, || false).unwrap();
        assert_fields(&fields(f.branches(), f.plenums()), &fields(reference.branches(), reference.plenums()), 3e-8);
        for (a, b) in wall_states(&w, 2).into_iter().zip(wall_states(&rw, 2)) {
            close(wall_energy(a), wall_energy(b), 3e-8);
        }
    }
}

#[test]
fn final_valve_cooling_respects_a_narrower_chemistry_interval() {
    let m = mechanism(0.0, 0.0);
    let make = || {
        let g = gas(&m);
        let branch = MixtureEuler1d::new(g.clone(), Length::new(1.0), 0.1,
            vec![g.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap(); 2], Boundaries::Periodic, || false).unwrap();
        let n = MixtureNetwork1d::new(vec![branch], vec![
            Plenum::new(0.2, g.conserved(&[2.0, 0.0], 0.0, 2.0).unwrap()).unwrap(),
            Plenum::new(0.1, g.conserved(&[0.2, 0.0], 0.0, 0.2).unwrap()).unwrap(),
        ], vec![], || false).unwrap();
        let v = NetworkValveOperator::new(&n, vec![ValveLink { first: 0, second: 1, orifice: orifice(0.1) }], || false).unwrap();
        (n, v)
    };
    let (mut probe, mut pv) = make();
    pv.advance_for(&mut probe, Time::new(0.01), ValveSolve::default(), || false).unwrap();
    let first = probe.gas().primitive(probe.plenums()[0].state(), || false).unwrap().temperature;
    assert_eq!(probe.step(Time::new(0.02), StepControl::default(), || false).unwrap().elapsed, Time::new(0.02));
    pv.advance_for(&mut probe, Time::new(0.01), ValveSolve::default(), || false).unwrap();
    let final_t = probe.gas().primitive(probe.plenums()[0].state(), || false).unwrap().temperature;
    assert!(first - final_t > 1e-5);
    let floor = 0.5 * (first + final_t);
    let (n, mut valves) = make();
    let mut f = ReactingNetwork1d::new(n, m.clone(),
        IsochoricKineticsSolve { minimum_temperature_k: floor, ..chemistry() }, || false).unwrap();
    let original = fields(f.branches(), f.plenums());
    let error = f.step_with_sources(Time::new(0.02),
        ReactiveStepControl { maximum_coupling_retries: 0, ..ReactiveStepControl::default() },
        NetworkSources { wall: None, valves: Some((&mut valves, ValveSolve::default())) }, || false).unwrap_err();
    match error {
        ReactiveNetworkError::RetryLimit { last, .. } => assert!(matches!(*last, ReactiveNetworkError::Cell { .. })),
        other => panic!("expected final valve caloric refusal, got {other:?}"),
    }
    assert_eq!(fields(f.branches(), f.plenums()), original);
    assert_eq!(f.time().value(), 0.0);
    assert!(valves.last_report().is_none());
    let recovered = f.step_with_sources(Time::new(0.02), ReactiveStepControl::default(),
        NetworkSources { wall: None, valves: Some((&mut valves, ValveSolve::default())) }, || false).unwrap();
    assert!(recovered.coupling_retries > 0 && recovered.elapsed.value() < 0.02);
    assert!(f.gas().primitive(f.plenums()[0].state(), || false).unwrap().temperature >= floor);
    assert_eq!(valves.last_report().unwrap().elapsed, recovered.elapsed);
}
