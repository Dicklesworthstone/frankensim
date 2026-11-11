//! G1/G3/G4: exact-vacuum filling with an evolving conducting wall assembly.
//! Synthetic source laws provide numerical references, not empirical validation.

use super::*;
use super::super::trajectory::{
    EffusionAdvanceControls, FiniteChannelFillAdvanceError,
    FiniteChannelFillAssemblyAdvanceCandidate, FiniteChannelFillAssemblyAdvanceStep,
    advance_candidate_with_conduction as advance,
};
use super::super::super::thermal::trajectory::FiniteChannelAdvanceError;
use crate::thermal::transport::volumetric::timeline::FixedVolumeConductionLink;
use fs_thermal::enthalpy::{CellId, ThermalLink};

const V: f64 = 1e-6;
const G: f64 = 0.1;

fn conductor<'l, 'p>(
    laws: [&'l FixedVolumeThermalLaw<'p>; 2],
) -> ChannelWallConduction<'l, 'p> {
    ChannelWallConduction::Cellwise(FixedVolumeConductionLink::try_aligned(
        &ThermalLink {
            id: 20, a: CellId(9), b: CellId(10), area_m2: 0.002,
            conductance_w_k: 0.0,
        },
        laws, [0.01; 2], 0.0,
    ).unwrap())
}
fn wall_pair<'l, 'p>(
    law: &'l FixedVolumeThermalLaw<'p>, temperatures: [f64; 2],
) -> BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
    BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(law, V, temperatures[0]).unwrap()),
        (10, FixedVolumeChannelWall::try_new_isochoric(law, V, temperatures[1]).unwrap()),
    ])
}
fn accuracy(scale: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: 0.05,
        minimum_step_s: 1e-6,
        maximum_step_s: 0.05,
        maximum_attempts: 16_384,
        maximum_accepted_intervals: 4096,
        relative_tolerance: 1e-3 * scale,
        temperature_tolerance_k: 0.05 * scale,
        pressure_tolerance_pa: 2.0 * scale,
        species_mass_tolerance_kg: 1e-10 * scale,
        stream_energy_tolerance_j: 5e-5 * scale,
    }
}
fn balance(value: &FiniteChannelFillBalance) {
    assert!(value.total_mass_residual_kg.abs() <= controls().exchange.mass_tolerance_kg);
    assert!(value.maximum_species_mass_residual_kg <= controls().exchange.mass_tolerance_kg);
    assert!(value.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(value.total_entropy_change_j_k >= -controls().entropy_tolerance_j_k);
}
fn temperatures(candidate: &FiniteChannelFillAssemblyAdvanceCandidate<'_, '_>) -> [f64; 4] {
    [
        candidate.donor().state().state().temperature,
        candidate.receiver().gas().unwrap().state().state().temperature,
        candidate.walls()[&9].temperature_k(),
        candidate.walls()[&10].temperature_k(),
    ]
}

#[test]
fn capture_pools_all_wall_heat_without_counting_conduction_as_molecular_energy() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = wall_pair(&law, [500.0, 300.0]);
    let wall_temperatures_k = [walls[&9].temperature_k(), walls[&10].temperature_k()];
    let conduction = [conductor([&law; 2])];
    let donor = gas(800.0, 100_000.0, 0.0);
    let receiver = empty(0.001);
    let path = channel(0.0);
    let baseline = step_candidate(
        &donor, &receiver, &walls[&9], link(&path), DT,
        energy(), energy(), controls(), || false,
    ).unwrap();
    let first = step_candidate_with_conduction(
        &donor, &receiver, &walls, link(&path), &conduction, DT,
        energy(), energy(), controls(), || false,
    ).unwrap();
    // Both sample the original channel wall. Side-wall heat only changes the
    // next material state; it is not an additional captured molecular beam.
    assert_eq!(first.receiver(), baseline.receiver());
    assert_eq!(first.donor(), baseline.donor());
    let FiniteChannelFillTransfer::Initial(report) = &first.step().transfer else {
        panic!("an exactly empty receiver requires first capture");
    };
    assert_eq!(report.reference_relative_energy_j, initial(&baseline).reference_relative_energy_j);
    assert_eq!(report.transport.walls.len(), 2);
    assert_eq!(report.transport.conduction.len(), 1);
    let incident = AREA * DT * 100_000.0
        / (2.0 * core::f64::consts::PI * M * R * 800.0).sqrt();
    let molecular = incident * (1.0 - path.geometry().direct_probability())
        * 2.0 * R * (800.0 - 500.0);
    let conductive = G * (500.0 - 300.0) * DT;
    near(report.transport.conduction[0].heat_a_to_b_j, conductive, 1e-14);
    let mut entropy = 0.0;
    for (id, q_molecular, q_conductive) in [(9, molecular, -conductive), (10, 0.0, conductive)] {
        let old = &walls[&id];
        let next = &first.walls()[&id];
        let row = report.transport.walls.iter().find(|row| row.id == id).unwrap();
        near(row.molecular_heat_into_wall_j, q_molecular, 1e-12);
        near(row.conduction_heat_into_wall_j, q_conductive, 1e-14);
        near(row.heat_into_wall_j, q_molecular + q_conductive, 1e-12);
        let q = q_molecular + q_conductive;
        let capacity = V * (1e6 + 700.0 * old.temperature_k());
        let delta = 2.0 * q / (capacity + (capacity * capacity + 2.0 * V * 700.0 * q).sqrt());
        near(next.temperature_k(), old.temperature_k() + delta, 2e-8);
        let ds = V * (1e6 * ((next.temperature_k() - old.temperature_k())
            / old.temperature_k()).ln_1p() + 700.0 * (next.temperature_k() - old.temperature_k()));
        near(row.entropy_change_j_k, ds, 2e-12);
        entropy += ds;
        assert_eq!(next.mass_kg(), old.mass_kg());
        assert_eq!(next.law().identity(), old.law().identity());
    }
    near(first.step().balance.wall_energy_change_j, molecular, 2e-9);
    near(first.step().balance.wall_entropy_change_j_k, entropy, 2e-12);
    balance(&first.step().balance);

    // The filled state uses the same native reciprocal transaction, with both
    // finite inventories and exactly the original prepared conductive link.
    let second = step_candidate_with_conduction(
        first.donor(), first.receiver(), first.walls(), link(&path), &conduction, DT,
        energy(), energy(), controls(), || false,
    ).unwrap();
    let native = thermal::step_candidate_with_conduction(
        &BTreeMap::from([
            (0, first.donor().clone()), (1, first.receiver().gas().unwrap().clone()),
        ]), first.walls(), GasTransportPaths::default(),
        &[FiniteChannelNetworkLink {
            id: 7, parcels: [0, 1], channel: &path, area_m2: AREA, wall: 9,
        }], &conduction, DT, &BTreeMap::from([(0, energy()), (1, energy())]),
        controls(), || false,
    ).unwrap();
    assert_eq!(second.donor(), &native.parcels()[&0]);
    assert_eq!(second.receiver().gas().unwrap(), &native.parcels()[&1]);
    assert!(matches!(&second.step().transfer,
        FiniteChannelFillTransfer::Exchange(actual) if actual == native.step()));
    for id in [9, 10] {
        assert_eq!(second.walls()[&id].specific_energy_j_kg(),
            native.walls()[&id].specific_energy_j_kg());
    }
    balance(&second.step().balance);
    assert!(matches!(receiver, EffusionReceiver::Evacuated(_)));
    // Both original walls are unchanged. try_new_isochoric stores the law's own
    // inverse of the stored energy (T to e to T, e8cf1188b7), so the requested
    // 500 K and 300 K come back within rounding, not necessarily bit for bit.
    assert_eq!(walls[&9].temperature_k(), wall_temperatures_k[0]);
    assert_eq!(walls[&10].temperature_k(), wall_temperatures_k[1]);
    assert!((wall_temperatures_k[0] - 500.0).abs() <= 2.0 * (500.0_f64.next_up() - 500.0));
    assert!((wall_temperatures_k[1] - 300.0).abs() <= 2.0 * (300.0_f64.next_up() - 300.0));
}

// Independent continuous n/U balances plus two nonlinear finite capacities.
// This RK4 reference never calls a library gas step, wall inversion or driver.
// At exact vacuum, the nonexistent receiver temperature is never evaluated.
fn continuous(path: &DiffuseMolecularChannel, duration: f64) -> [f64; 4] {
    let tau = path.geometry().transmission_probability();
    let hit = path.geometry().wall_transmitted_probability();
    let cv = 1.5 * R;
    let derivative = |y: [f64; 6]| {
        let [nl, ul, nr, ur, tw, tb] = y;
        let tl = ul / (cv * nl);
        let left = AREA * nl / 0.001
            * (R * tl / (2.0 * core::f64::consts::PI * M)).sqrt();
        let right = if nr > 0.0 {
            let tr = ur / (cv * nr);
            let rate = AREA * nr / 0.001
                * (R * tr / (2.0 * core::f64::consts::PI * M)).sqrt();
            [rate, rate * tr]
        } else {
            [0.0; 2]
        };
        let ql = 2.0 * R * ((1.0 - tau) * left * (tw - tl)
            + hit * (tw * right[0] - right[1]));
        let qr = 2.0 * R * ((1.0 - tau) * (tw * right[0] - right[1])
            + hit * left * (tw - tl));
        let amount = tau * (left - right[0]);
        let carried = 2.0 * R * tau * (left * tl - right[1]);
        let conducted = G * (tw - tb);
        [
            -amount, -carried + ql, amount, carried + qr,
            (-(ql + qr) - conducted) / (V * (1e6 + 700.0 * tw)),
            conducted / (V * (1e6 + 700.0 * tb)),
        ]
    };
    let n0 = 100_000.0 * 0.001 / (R * 400.0);
    let mut y = [n0, n0 * cv * 400.0, 0.0, 0.0, 900.0, 300.0];
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..6 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    [y[1] / (cv * y[0]), y[3] / (cv * y[2]), y[4], y[5]]
}

#[test]
fn conducting_fill_converges_to_coupled_balances_and_replays_every_original_wall() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = wall_pair(&law, [900.0, 300.0]);
    let conduction = [conductor([&law; 2])];
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = EvacuatedGasVolume::new(0.001, QueryPoint::new()).unwrap();
    let path = channel(0.0);
    let duration = 0.4;
    let run = |scale| advance(
        &donor, &receiver, &walls, link(&path), &conduction, duration,
        energy(), energy(), controls(), accuracy(scale), || false,
    ).unwrap();
    let loose = run(1.0);
    let tight = run(0.04);
    let reference = continuous(&path, duration);
    let discrepancy = |candidate: &FiniteChannelFillAssemblyAdvanceCandidate<'_, '_>| {
        temperatures(candidate).into_iter().zip(reference)
            .map(|(a, b)| (a - b).abs()).sum::<f64>()
    };
    assert!(discrepancy(&tight) < 0.65 * discrepancy(&loose));
    let FiniteChannelFillAssemblyAdvanceStep::Filling(report) = tight.step() else {
        panic!("an open channel must retain its first capture and reciprocal history");
    };
    let FiniteChannelFillAssemblyAdvanceStep::Filling(coarse) = loose.step() else {
        panic!("an open channel must fill");
    };
    assert!(report.accepted_intervals > coarse.accepted_intervals);
    assert!(report.maximum_accepted_error <= 1.0);
    assert!(report.first_reciprocal.as_ref().unwrap().channels[0]
        .transmitted.species[0].right_to_left.mass_rate_kg_s > 0.0);
    assert!(tight.walls()[&9].temperature_k() < 900.0);
    assert!(tight.walls()[&10].temperature_k() > 300.0);
    balance(tight.step().balance());
    balance(loose.step().balance());

    let first = step_candidate_with_conduction(
        &donor, &EffusionReceiver::Evacuated(receiver.clone()), &walls,
        link(&path), &conduction, report.first_fill.duration_s,
        energy(), energy(), controls(), || false,
    ).unwrap();
    assert_eq!(first.step(), &report.first_fill);
    let (mut left, mut right, mut current_walls, _) = first.into_parts();
    let mut elapsed = Dd::from_f64(report.first_fill.duration_s);
    for accepted in report.first_reciprocal.iter().chain(&report.continuation) {
        let next = step_candidate_with_conduction(
            &left, &right, &current_walls, link(&path), &conduction, accepted.duration_s,
            energy(), energy(), controls(), || false,
        ).unwrap();
        assert!(matches!(&next.step().transfer,
            FiniteChannelFillTransfer::Exchange(actual) if actual == accepted));
        (left, right, current_walls, _) = next.into_parts();
        elapsed = elapsed + Dd::from_f64(accepted.duration_s);
    }
    assert_eq!(&left, tight.donor());
    assert_eq!(&right, tight.receiver());
    for id in [9, 10] {
        assert_eq!(current_walls[&id].specific_energy_j_kg(),
            tight.walls()[&id].specific_energy_j_kg());
    }
    near(elapsed.to_f64(), duration, 1e-12);
    let actual_storage: f64 = [9, 10].into_iter().map(|id| {
        let old = &walls[&id];
        let t = tight.walls()[&id].temperature_k();
        V * (t - old.temperature_k()) * (1e6 + 350.0 * (t + old.temperature_k()))
    }).sum();
    near(report.balance.wall_energy_change_j, actual_storage, 2e-8);
}

#[test]
fn closed_passage_resolves_conduction_without_inventing_a_receiver_state() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = wall_pair(&law, [900.0, 300.0]);
    let conduction = [conductor([&law; 2])];
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = EvacuatedGasVolume::new(0.001, QueryPoint::new()).unwrap();
    let path = channel(0.0);
    let closed = FiniteChannelFillLink { area_m2: 0.0, ..link(&path) };
    let duration = 0.5;
    let run = |scale| advance(
        &donor, &receiver, &walls, closed, &conduction, duration,
        energy(), energy(), controls(),
        EffusionAdvanceControls {
            initial_step_s: duration, maximum_step_s: duration,
            minimum_step_s: 1e-7, relative_tolerance: 1e-12,
            // Only the signed conductive heat accuracy drives refinement here.
            temperature_tolerance_k: 1e6, pressure_tolerance_pa: 1e6,
            stream_energy_tolerance_j: 1e-3 * scale, ..accuracy(1.0)
        }, || false,
    ).unwrap();
    let loose = run(1.0);
    let tight = run(0.01);
    let FiniteChannelFillAssemblyAdvanceStep::ClosedPath {
        thermal: coarse, ..
    } = loose.step() else { panic!("closed passage must use the native vacuum history"); };
    let FiniteChannelFillAssemblyAdvanceStep::ClosedPath {
        thermal: report, balance: closed_balance,
    } = tight.step() else { panic!("zero captured molecules must remain exact vacuum"); };
    let reference = 600.0 + 300.0 * (-2.0 * G * duration).exp();
    let error = |result: &FiniteChannelFillAssemblyAdvanceCandidate<'_, '_>| {
        (result.walls()[&9].temperature_k() - reference).abs()
    };
    assert!(error(&tight) < 0.5 * error(&loose));
    near(tight.walls()[&9].temperature_k(), reference, 0.02);
    near(tight.walls()[&10].temperature_k(), 1200.0 - reference, 0.02);
    assert!(report.accepted_intervals > coarse.accepted_intervals);
    assert!(report.maximum_accepted_error <= 1.0);
    assert_eq!(tight.donor(), &donor);
    assert_eq!(tight.receiver(), &EffusionReceiver::Evacuated(receiver.clone()));
    assert_eq!(tight.step().duration_s(), duration);
    assert_eq!(report.balance.gas.escaped_mass_kg, 0.0);
    for step in &report.accepted_steps {
        assert_eq!(step.conduction.len(), 1);
        assert!(step.conduction[0].heat_a_to_b_j > 0.0);
        for wall in &step.walls {
            assert_eq!(wall.molecular_heat_into_wall_j, 0.0);
            assert_eq!(wall.heat_into_wall_j, wall.conduction_heat_into_wall_j);
        }
    }
    let entropy = (tight.walls()[&9].temperature_k() / 900.0).ln()
        + (tight.walls()[&10].temperature_k() / 300.0).ln();
    near(closed_balance.wall_entropy_change_j_k, entropy, 1e-10);
    assert!(closed_balance.total_entropy_change_j_k > 0.0);
    balance(closed_balance);
    assert_eq!(walls[&9].temperature_k(), 900.0);
    assert_eq!(walls[&10].temperature_k(), 300.0);
}

#[test]
fn source_binding_work_limits_and_late_cancellation_leave_every_inventory_replayable() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let other_basis = FixedVolumeThermalLaw::try_new(&source, 1200.0, 200.0)
        .unwrap().with_sensible_entropy();
    let walls = wall_pair(&law, [900.0, 300.0]);
    let conduction = [conductor([&law; 2])];
    let incompatible = [conductor([&law, &other_basis])];
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = EvacuatedGasVolume::new(0.001, QueryPoint::new()).unwrap();
    let original_donor = donor.clone();
    let original_receiver = receiver.clone();
    let path = channel(0.0);
    assert!(matches!(advance(
        &donor, &receiver, &walls, link(&path), &incompatible, 0.02,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ), Err(FiniteChannelFillAdvanceError::Fill {
        error: FiniteChannelFillError::Physical(FiniteChannelError::Invalid(_)), ..
    })));
    assert!(matches!(advance(
        &donor, &receiver, &walls, link(&path), &conduction, 0.4,
        energy(), energy(), controls(),
        EffusionAdvanceControls { maximum_attempts: 3, ..accuracy(1.0) }, || false,
    ), Err(FiniteChannelFillAdvanceError::Advance(
        FiniteChannelAdvanceError::WorkBudget { attempted_steps: 3, .. }
    ))));
    let mut total = 0;
    let successful = advance(
        &donor, &receiver, &walls, link(&path), &conduction, 0.02,
        energy(), energy(), controls(), accuracy(1.0),
        || { total += 1; false },
    ).unwrap();
    let mut seen = 0;
    assert!(matches!(advance(
        &donor, &receiver, &walls, link(&path), &conduction, 0.02,
        energy(), energy(), controls(), accuracy(1.0),
        || { seen += 1; seen == total },
    ), Err(FiniteChannelFillAdvanceError::Cancelled)));
    let replay = advance(
        &donor, &receiver, &walls, link(&path), &conduction, 0.02,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ).unwrap();
    assert_eq!(replay.step(), successful.step());
    assert_eq!(replay.donor(), successful.donor());
    assert_eq!(replay.receiver(), successful.receiver());
    for id in [9, 10] {
        assert_eq!(replay.walls()[&id].specific_energy_j_kg(),
            successful.walls()[&id].specific_energy_j_kg());
    }
    assert_eq!(donor, original_donor);
    assert_eq!(receiver, original_receiver);
    assert_eq!(walls[&9].temperature_k(), 900.0);
    assert_eq!(walls[&10].temperature_k(), 300.0);
}

#[test]
fn original_single_wall_apis_delegate_without_changing_their_results() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, V, 900.0).unwrap();
    let walls = BTreeMap::from([(9, wall.clone())]);
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = EvacuatedGasVolume::new(0.001, QueryPoint::new()).unwrap();
    let empty = EffusionReceiver::Evacuated(receiver.clone());
    let path = channel(0.0);
    let old = step_candidate(
        &donor, &empty, &wall, link(&path), DT,
        energy(), energy(), controls(), || false,
    ).unwrap();
    let new = step_candidate_with_conduction(
        &donor, &empty, &walls, link(&path), &[], DT,
        energy(), energy(), controls(), || false,
    ).unwrap();
    assert_eq!(new.step(), old.step());
    assert_eq!(new.receiver(), old.receiver());
    assert_eq!(new.donor(), old.donor());
    assert_eq!(new.walls()[&9].specific_energy_j_kg(), old.wall().specific_energy_j_kg());
    let old = super::super::trajectory::advance_candidate(
        &donor, &receiver, &wall, link(&path), 0.02,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ).unwrap();
    let new = advance(
        &donor, &receiver, &walls, link(&path), &[], 0.02,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ).unwrap();
    assert_eq!(new.step(), &FiniteChannelFillAssemblyAdvanceStep::Filling(old.step().clone()));
    assert_eq!(new.receiver(), old.receiver());
    assert_eq!(new.donor(), old.donor());
    assert_eq!(new.walls()[&9].specific_energy_j_kg(), old.wall().specific_energy_j_kg());
}
