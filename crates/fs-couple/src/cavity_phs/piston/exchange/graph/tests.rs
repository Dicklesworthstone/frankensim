//! Actual shared gas inventory, caloric reconstruction and piston momentum owners.
use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{GasNozzle, GasNozzleStepSolve};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_ports::{ModuleWindowAudit, PortedModule};

fn solve() -> MovingGasSolve {
    MovingGasSolve { exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, total_energy_tolerance_j: 1e-7 }
}
fn heated(x: f64, heat: f64) -> ClosedGasParcel {
    let mut gas = gas(x, 0.001);
    gas.add_heat(heat, controls().gas, || false).unwrap(); gas
}
fn link(id: u64, a: u64, b: u64) -> SnapshotNozzleLink {
    SnapshotNozzleLink { id, parcels: [a,b], nozzle: GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 },
        nozzle_solve: controls().gas }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-10*b.abs().max(1.0), "{a} != {b}"); }
fn fixture() -> (BTreeMap<u64, GasPiston>, BTreeMap<u64, ClosedGasParcel>,
    Vec<SnapshotNozzleLink>, BTreeMap<u64, GasEnergySolve>, BTreeMap<u64, MovingGasDrive>) {
    let wall = GasPiston::new(heated(0.75, 30.0), 1.0, 0.01, 0.02).unwrap();
    let pistons = BTreeMap::from([(2, wall.clone()), (3, wall)]);
    let reservoirs = BTreeMap::from([(1, heated(0.25, 150.0)), (4, gas(0.75, 0.001))]);
    let links = vec![link(10,1,2), link(20,1,3), link(30,2,4), link(40,3,4)];
    let endpoints = (1..=4).map(|id| (id, controls().gas)).collect();
    let drives = [2,3].map(|id| (id, MovingGasDrive { input: drive(0.1), solve: controls() })).into();
    (pistons, reservoirs, links, endpoints, drives)
}

#[test]
fn branching_inlets_and_outlets_drive_symmetric_walls_and_one_material_history() {
    let (pistons, reservoirs, links, endpoints, drives) = fixture();
    let proposal = step_candidate(&pistons, &reservoirs, &links, drive(0.0).duration_s,
        &endpoints, &drives, solve(), || false).unwrap();
    assert_eq!(proposal.pistons()[&2], proposal.pistons()[&3]);
    assert_eq!(proposal.pistons()[&2].ports, proposal.pistons()[&3].ports);
    for id in [2,3] {
        let old = &pistons[&id]; let new = &proposal.pistons()[&id];
        assert!(new.velocity_m_s() > old.velocity_m_s());
        assert!(new.gas().volume_m3() > old.gas().volume_m3());
        assert_ne!(new.gas().state().state().specific_gas_constant, old.gas().state().state().specific_gas_constant);
        assert_ne!(new.gas().state().state().sound_speed, old.gas().state().state().sound_speed);
        assert_eq!(new.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert!(matches!(new.audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
        assert_eq!(new.last_window_transfers().unwrap().len(), 2);
        assert_eq!(proposal.step().pistons[&id].gas.supplied_heat_j, 0.1);
    }
    // Outflow carries the accepted chamber mixture, not material just received.
    let out = proposal.step().flow.links.iter().find(|e| e.id == 30).unwrap();
    assert_eq!(out.flow.upstream_enthalpy_j_kg, pistons[&2].gas().state().specific_enthalpy_j_kg());
    assert!(out.transferred_mass_kg > 0.0);
    let source_loss = reservoirs[&1].mass_kg()-proposal.reservoirs()[&1].mass_kg();
    near(source_loss, proposal.step().flow.links.iter().filter(|e| e.donor == 1)
        .map(|e| e.transferred_mass_kg).sum());
    assert!(proposal.step().total_energy_residual_j.abs() < solve().total_energy_tolerance_j);
    assert!(proposal.step().total_mass_residual_kg.abs() < solve().exchange.mass_tolerance_kg);
}

#[test]
fn link_reporting_order_does_not_select_another_physical_update() {
    let (pistons, reservoirs, links, endpoints, drives) = fixture();
    let run = |links: &[SnapshotNozzleLink]| step_candidate(&pistons, &reservoirs, links,
        drive(0.0).duration_s, &endpoints, &drives, solve(), || false).unwrap();
    let original = run(&links);
    let mut reversed = links.clone(); reversed.reverse();
    for (i, link) in reversed.iter_mut().enumerate() { link.id = 100+i as u64; link.parcels.reverse(); }
    let alternate = run(&reversed);
    assert_eq!(original.pistons(), alternate.pistons());
    assert_eq!(original.reservoirs(), alternate.reservoirs());
    assert_eq!(original.step().pistons, alternate.step().pistons);
    for id in pistons.keys() { assert_eq!(original.pistons()[id].ports, alternate.pistons()[id].ports); }
}

#[test]
fn zero_flow_retains_exact_closed_mechanics_and_fixed_only_matches_flow_owner() {
    let (pistons, reservoirs, mut links, endpoints, drives) = fixture();
    for link in &mut links { link.nozzle.area_m2 = 0.0; }
    let trial = step_candidate(&pistons, &reservoirs, &links, drive(0.0).duration_s,
        &endpoints, &drives, solve(), || false).unwrap();
    for (&id, old) in &pistons {
        let mut expected = old.clone();
        let report = expected.step(drives[&id].input, drives[&id].solve, || false).unwrap();
        assert_eq!(&trial.pistons()[&id], &expected);
        assert_eq!(trial.pistons()[&id].ports, expected.ports);
        assert_eq!(trial.step().pistons[&id], report);
    }
    assert_eq!(trial.reservoirs(), &reservoirs);
    let mut all = reservoirs.clone();
    all.extend(pistons.iter().map(|(&id, p)| (id, p.gas().clone())));
    for link in &mut links { link.nozzle.area_m2 = 1e-5; }
    let expected = snapshot::step_candidate(&all, &links, drive(0.0).duration_s,
        &endpoints, solve().exchange, || false).unwrap();
    let actual = step_candidate(&BTreeMap::new(), &all, &links, drive(0.0).duration_s,
        &endpoints, &BTreeMap::new(), solve(), || false).unwrap();
    assert!(actual.pistons().is_empty());
    assert_eq!(actual.reservoirs(), expected.parcels());
    assert_eq!(&actual.step().flow, expected.step());
}

#[test]
fn single_connection_agrees_with_existing_finite_reservoir_motion() {
    let mut wall = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    let mut reservoir = heated(0.25, 100.0);
    let input = drive(0.2); let restriction = link(1, 1, 2);
    let proposal = step_candidate(&BTreeMap::from([(2, wall.clone())]),
        &BTreeMap::from([(1, reservoir.clone())]), &[restriction], input.duration_s,
        &BTreeMap::from([(1,controls().gas),(2,controls().gas)]),
        &BTreeMap::from([(2,MovingGasDrive { input, solve: controls() })]), solve(), || false).unwrap();
    let exchange = GasNozzleStepSolve { nozzle: controls().gas, maximum_donor_fraction: 0.1,
        transfer: GasMassTransferSolve { donor: controls().gas, receiver: controls().gas,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } };
    wall.step_with_nozzle(&mut reservoir, restriction.nozzle, input, exchange, controls(), 1e-7, || false).unwrap();
    let actual = &proposal.pistons()[&2];
    near(actual.velocity_m_s(), wall.velocity_m_s());
    near(actual.gas().volume_m3(), wall.gas().volume_m3());
    near(actual.gas().state().state().pressure, wall.gas().state().state().pressure);
    near(proposal.reservoirs()[&1].state().state().temperature, reservoir.state().state().temperature);
}

#[test]
fn explicit_sources_are_not_counted_twice_and_history_survives_source_shutoff() {
    let (mut pistons, mut reservoirs, links, endpoints, mut drives) = fixture();
    let initial_pistons = pistons.clone(); let initial_reservoirs = reservoirs.clone();
    let (mut supplied, mut outside, mut applied) = (0.0, 0.0, 0.0);
    for i in 0..40 {
        for drive in drives.values_mut() {
            drive.input.heat_j = if i < 20 { 0.05 } else { 0.0 };
            drive.input.external_force_n = 2.0;
        }
        let proposal = step_candidate(&pistons, &reservoirs, &links, drive(0.0).duration_s,
            &endpoints, &drives, solve(), || false).unwrap();
        supplied += proposal.step().supplied_heat_j;
        outside += proposal.step().exterior_pressure_work_j;
        applied += proposal.step().external_force_work_j;
        let (p, r, _) = proposal.into_parts(); pistons = p; reservoirs = r;
    }
    near(supplied, 2.0);
    let du: f64 = pistons.iter().map(|(id,p)| stored_change(initial_pistons[id].gas(), p.gas())).sum::<f64>()
        + reservoirs.iter().map(|(id,r)| stored_change(&initial_reservoirs[id], r)).sum::<f64>();
    let dk: f64 = pistons.iter().map(|(id,p)| {
        let old = &initial_pistons[id];
        0.5*old.mass_kg()*(p.velocity_m_s()-old.velocity_m_s())*(p.velocity_m_s()+old.velocity_m_s())
    }).sum();
    assert!((du+dk+outside-applied-supplied).abs() < 2e-7);
    assert!(reservoirs[&1].mass_kg() < initial_reservoirs[&1].mass_kg());
    assert!(pistons[&2].gas().volume_m3() > initial_pistons[&2].gas().volume_m3());
}

#[test]
fn final_cancellation_late_wall_failure_and_bad_bindings_preserve_all_owners() {
    let (pistons, reservoirs, links, endpoints, drives) = fixture();
    let saved = (pistons.clone(), reservoirs.clone());
    let run = |drives: &BTreeMap<u64,MovingGasDrive>, cancel: &mut dyn FnMut() -> bool| {
        step_candidate(&pistons,&reservoirs,&links,drive(0.0).duration_s,&endpoints,drives,solve(),cancel)
    };
    let mut total = 0;
    run(&drives, &mut || { total += 1; false }).unwrap();
    for stop in [1, total-1, total] {
        let mut count = 0;
        assert_eq!(run(&drives, &mut || { count += 1; count == stop }), Err(Error::Cancelled));
    }
    let mut bad = drives.clone(); bad.get_mut(&3).unwrap().input.heat_j = 1e12;
    assert!(matches!(run(&bad,&mut ||false), Err(Error::Piston { id: 3, .. })));
    bad = drives.clone(); bad.get_mut(&3).unwrap().input.duration_s *= 2.0;
    assert!(matches!(run(&bad,&mut ||false), Err(Error::Invalid(_))));
    bad = drives.clone(); bad.remove(&3);
    assert!(matches!(run(&bad,&mut ||false), Err(Error::Invalid(_))));
    let mut overlapping = reservoirs.clone(); overlapping.insert(2, gas(0.75,0.001));
    assert!(step_candidate(&pistons,&overlapping,&links,drive(0.0).duration_s,
        &endpoints,&drives,solve(),||false).is_err());
    assert_eq!((&pistons,&reservoirs),(&saved.0,&saved.1));
    for id in pistons.keys() { assert_eq!(pistons[id].ports, saved.0[id].ports); }
}
