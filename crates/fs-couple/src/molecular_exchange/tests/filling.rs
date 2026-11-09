//! G0/G1/G3/G4: retained beam drift at exact-vacuum startup, original U+K and
//! reciprocal continuation. Synthetic Maxwellian/caloric fixtures, not gas data.

use super::*;
use crate::molecular_exchange::filling::{
    self, MolecularFillCandidate, MolecularFillError, MolecularFillStage, MolecularReceiver,
};
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::parcel::evacuated::EvacuatedGasVolume;

const RECEIVER_VOLUME: f64 = 1e-6;

fn empty() -> MolecularReceiver {
    MolecularReceiver::Evacuated(
        EvacuatedGasVolume::new(
            RECEIVER_VOLUME,
            QueryPoint::new().with("receiver-lot", 19.0).unwrap(),
        )
        .unwrap(),
    )
}
fn run(
    donor: &MovingGasParcel,
    receiver: &MolecularReceiver,
    law: &EffusiveDriftInterface,
    port: MolecularExchangeSurface,
) -> MolecularFillCandidate {
    filling::step_candidate(donor, receiver, law, port, DT, controls(), || false).unwrap()
}
fn first(result: &MolecularFillCandidate) -> &filling::MolecularFirstFillStep {
    match &result.step().stage {
        MolecularFillStage::FirstFill(step) => step,
        _ => panic!("expected original first-fill receipt"),
    }
}
fn check_total(before: &MovingGasParcel, result: &MolecularFillCandidate) {
    let received = result.receiver().gas().unwrap();
    near(
        result.donor().parcel().mass_kg() + received.parcel().mass_kg(),
        before.parcel().mass_kg(),
        1e-13,
    );
    vector_near(
        core::array::from_fn(|d| {
            result.donor().momentum_kg_m_s()[d] + received.momentum_kg_m_s()[d]
        }),
        before.momentum_kg_m_s(),
        1e-12,
    );
    near(
        internal(result.donor())
            + internal(received)
            + result.donor().kinetic_energy_j()
            + received.kinetic_energy_j(),
        internal(before) + before.kinetic_energy_j(),
        1e-7,
    );
    assert_eq!(received.parcel().volume_m3(), RECEIVER_VOLUME);
    assert_eq!(
        received.parcel().point().axes().get("receiver-lot"),
        Some(&19.0)
    );
    let step = first(result);
    assert!(step.absolute_energy_defect_j <= step.energy_allowance_j);
    assert!(step.energy_residual_j.abs() <= step.energy_allowance_j);
    assert!(step.laboratory_energy_residual_j.abs() <= step.energy_allowance_j);
}

#[test]
fn stationary_maxwellian_capture_retains_bulk_motion_and_matches_independent_temperature() {
    for cp_r in [2.5, 3.5] {
        let text = format!("A 1 0.028 {cp_r} 0 1.7e-5 0.025");
        let model = catalog(&text, 0.0);
        let law = interface(&model, 2000.0);
        let donor = gas(model, T, P, [0.0; 3]);
        let receiver = empty();
        let result = run(&donor, &receiver, &law, surface());
        let filled = result.receiver().gas().unwrap();
        let mass = crossing_mass(T, P, M);
        let mean = (core::f64::consts::PI * R * T / (2.0 * M)).sqrt();
        let cv = cp_r - 1.0;
        let expected_t = T * (cv + 0.5 - core::f64::consts::PI / 4.0) / cv;
        near(filled.parcel().mass_kg(), mass, 1e-15);
        vector_near(filled.velocity_m_s(), [mean, 0.0, 0.0], 1e-9);
        near(temperature(filled), expected_t, 1e-7);
        near(
            filled.parcel().state().state().pressure,
            mass / M * R * expected_t / RECEIVER_VOLUME,
            1e-5,
        );
        let report = first(&result);
        near(
            report.receiver.kinetic_energy_j,
            mass * core::f64::consts::PI * R * T / (4.0 * M),
            1e-10,
        );
        assert_eq!(
            report.receiver.caloric.additional_energy_j,
            -report.receiver.kinetic_energy_j
        );
        assert!(temperature(filled) < T);
        check_total(&donor, &result);
        assert_eq!(receiver, empty());
        assert_eq!(temperature(&donor), T);
    }
}

#[test]
fn mixture_separation_and_tangential_drift_preserve_original_sources_under_common_boosts() {
    let shift = [1000.0, -700.0, 300.0];
    let source_velocity = [0.0, 35.0, -20.0];
    for datum in [0.0, 1e12] {
        let model = catalog(MIXTURE, datum);
        let law = interface(&model, 2000.0);
        let donor = gas(model.clone(), T, P, source_velocity);
        let result = run(&donor, &empty(), &law, surface());
        let filled = result.receiver().gas().unwrap();
        let inventory = filled.parcel().species_inventory_kg().unwrap();
        let a = SpeciesId::new("A").unwrap();
        let b = SpeciesId::new("B").unwrap();
        let ma = 0.5 * crossing_mass(T, P, 0.028);
        let mb = 0.5 * crossing_mass(T, P, 0.004);
        near(inventory[&a], ma, 1e-15);
        near(inventory[&b], mb, 1e-15);
        let normal = (ma * (core::f64::consts::PI * R * T / (2.0 * 0.028)).sqrt()
            + mb * (core::f64::consts::PI * R * T / (2.0 * 0.004)).sqrt())
            / (ma + mb);
        vector_near(filled.velocity_m_s(), [normal, 35.0, -20.0], 1e-8);
        let caloric = 1.5 * R * (ma / 0.028 + mb / 0.004);
        let crossing = 2.0 * R * T * (ma / 0.028 + mb / 0.004);
        let expected_t = (crossing - 0.5 * (ma + mb) * normal * normal) / caloric;
        near(temperature(filled), expected_t, 1e-6);
        let boosted = gas(
            model,
            T,
            P,
            core::array::from_fn(|d| source_velocity[d] + shift[d]),
        );
        let moved = run(
            &boosted,
            &empty(),
            &law,
            MolecularExchangeSurface {
                velocity_world_m_s: shift,
                ..surface()
            },
        );
        let moved_gas = moved.receiver().gas().unwrap();
        near(temperature(moved_gas), temperature(filled), 1e-7);
        near(
            temperature(moved.donor()),
            temperature(result.donor()),
            1e-7,
        );
        vector_near(
            moved_gas.velocity_m_s(),
            core::array::from_fn(|d| filled.velocity_m_s()[d] + shift[d]),
            1e-8,
        );
        assert_eq!(
            moved_gas.parcel().model().thermodynamics().components(),
            filled.parcel().model().thermodynamics().components()
        );
        if datum == 0.0 {
            check_total(&donor, &result);
            check_total(&boosted, &moved);
        }
    }
}

#[test]
fn first_capture_continues_with_the_exact_original_reciprocal_pair() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let donor = gas(model, T, P, [2.0, 3.0, 0.0]);
    let initial = run(&donor, &empty(), &law, surface());
    let receiver = initial.receiver().gas().unwrap();
    let duration = 0.1 * DT;
    let native = super::step_candidate(
        [initial.donor(), receiver],
        &law,
        surface(),
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let continued = filling::step_candidate(
        initial.donor(),
        initial.receiver(),
        &law,
        surface(),
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let MolecularFillStage::Reciprocal(report) = &continued.step().stage else {
        panic!("expected reciprocal stage")
    };
    assert_eq!(continued.donor(), &native.gases()[0]);
    assert_eq!(continued.receiver().gas().unwrap(), &native.gases()[1]);
    assert_eq!(report, native.step());
    assert!(report.exchange.species[0].mass_kg[1] > 0.0);
}

#[test]
fn closed_paths_preserve_true_vacuum_and_do_not_query_inactive_sources() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 0.0);
    let donor = gas(model, T, P, [100.0, 20.0, 0.0]);
    let receiver = empty();
    let closed = MolecularExchangeSurface {
        area_m2: 0.0,
        normal_left_to_right: [f64::NAN; 3],
        velocity_world_m_s: [f64::NAN; 3],
    };
    let result = run(&donor, &receiver, &law, closed);
    assert_eq!(result.donor(), &donor);
    assert_eq!(result.receiver(), &receiver);
    let MolecularFillStage::Closed(flux) = &result.step().stage else {
        panic!("expected exact vacuum")
    };
    assert!(flux.species.is_empty());
    let rejected =
        filling::step_candidate(&donor, &receiver, &law, surface(), DT, controls(), || false);
    assert!(matches!(
        rejected,
        Err(MolecularFillError::Exchange(
            MolecularExchangeError::Material(_)
        ))
    ));
    let resting = gas(catalog(PURE, 0.0), T, P, [0.0; 3]);
    let zero = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                resting.parcel().model().clone(),
                vec![
                    EffusiveSpeciesTransmission::try_new(
                        SpeciesId::new("A").unwrap(),
                        0.0,
                        "synthetic blocked transmission",
                    )
                    .unwrap(),
                ],
                law.momentum_interface().interface().support(),
                law.momentum_interface().interface().regime().clone(),
            )
            .unwrap(),
            "synthetic inactive angular transmission",
        )
        .unwrap(),
        0.0,
        "synthetic rest-only blocked aperture",
    )
    .unwrap();
    let blocked = run(&resting, &receiver, &zero, surface());
    assert_eq!(blocked.donor(), &resting);
    assert_eq!(blocked.receiver(), &receiver);
    assert!(matches!(
        &blocked.step().stage,
        MolecularFillStage::Closed(_)
    ));
}

#[test]
fn original_and_final_drift_support_gross_inventory_and_late_cancellation_are_atomic() {
    let model = catalog("A 1 0.028 3.5 0.001 1.7e-5 0.025", 0.0);
    let donor = gas(model.clone(), T, P, [0.0; 3]);
    let receiver = empty();
    let saved = donor.clone();
    // Initial zero drift is supported, but the actual captured mean beam is not.
    let narrow = interface(&model, 100.0);
    assert!(matches!(
        filling::step_candidate(
            &donor,
            &receiver,
            &narrow,
            surface(),
            DT,
            controls(),
            || false
        ),
        Err(MolecularFillError::Exchange(
            MolecularExchangeError::Material(_)
        ))
    ));
    let law = interface(&model, 2000.0);
    let limits = MolecularExchangeControls {
        maximum_species_fraction: 1e-6,
        ..controls()
    };
    assert!(matches!(
        filling::step_candidate(&donor, &receiver, &law, surface(), DT, limits, || false),
        Err(MolecularFillError::Exchange(
            MolecularExchangeError::Withdrawal { side: 0, .. }
        ))
    ));
    let mut count = 0;
    let accepted =
        filling::step_candidate(&donor, &receiver, &law, surface(), DT, controls(), || {
            count += 1;
            false
        })
        .unwrap();
    let defect = first(&accepted).absolute_energy_defect_j;
    assert!(defect > 0.0);
    let tight = MolecularExchangeControls {
        energy_absolute_j: 0.5 * defect,
        relative_tolerance: 0.0,
        ..controls()
    };
    assert!(matches!(
        filling::step_candidate(&donor, &receiver, &law, surface(), DT, tight, || false),
        Err(MolecularFillError::Exchange(
            MolecularExchangeError::Balance { .. }
        ))
    ));
    let mut calls = 0;
    let cancelled =
        filling::step_candidate(&donor, &receiver, &law, surface(), DT, controls(), || {
            calls += 1;
            calls == count
        });
    assert!(matches!(
        cancelled,
        Err(MolecularFillError::Exchange(
            MolecularExchangeError::Cancelled
        ))
    ));
    assert_eq!(calls, count);
    assert_eq!(donor, saved);
    assert_eq!(receiver, empty());
}
