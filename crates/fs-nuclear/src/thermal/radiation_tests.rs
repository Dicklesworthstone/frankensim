//! G0/G1/G3/G4/G5: emitted energy, nonlocal calorimetry, replay and joint refusal.
//! All response fractions and nuclear coefficients are manufactured equation data.

use super::*;
use crate::decay::thermal::DecayCellSource;
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use crate::radiation::{
    RadiationDeposit, RadiationResponse, RadiationResponseLimits, RadiationSourceResponse,
};
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;

fn row(source_id: u64, deposits: &[(u64, f64)], escape_fraction: f64) -> RadiationSourceResponse {
    RadiationSourceResponse {
        source_id,
        deposits: deposits
            .iter()
            .map(|&(receiver_id, fraction)| RadiationDeposit {
                receiver_id,
                fraction,
            })
            .collect(),
        escape_fraction,
    }
}

fn response(rows: Vec<RadiationSourceResponse>) -> Arc<RadiationResponse> {
    Arc::new(
        RadiationResponse::new(
            "manufactured-absorbed-energy-response-v1".into(),
            rows,
            RadiationResponseLimits {
                max_sources: 8,
                max_receivers: 8,
                max_transfers: 16,
            },
        )
        .unwrap(),
    )
}

fn receiver(cell_id: u64, material_id: u64) -> RadiationReceiver {
    RadiationReceiver {
        cell_id: CellId(cell_id),
        material_id: MaterialId(material_id),
    }
}

#[test]
fn g1_prompt_and_shutdown_radiation_heat_original_receiver_materials() {
    let initial = thermal(
        vec![
            cell(1, 1, 2.0, 1200.0, 0.0),
            cell(2, 2, 3.0, 600.0, 3.0),
            cell(3, 1, 1.0, 1200.0, 0.0),
        ],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = FissionHeating::new(
        &initial,
        vec![source(1, 1, 5.0, &[(2.0, 0.0)], true)],
        budget(),
    )
    .unwrap()
    .with_radiation_response(
        &initial,
        response(vec![row(1, &[(3, 0.25), (2, 0.5)], 0.25)]),
        vec![receiver(3, 1), receiver(2, 2)],
    )
    .unwrap();
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    // Ten fissions create 40 J of delayed storage. Its 1/4 radiation channel
    // decays during both the two-second exposure and the three-second shutdown.
    let stored = 40.0 * -(-1.0_f64).exp_m1() * (-1.5_f64).exp();
    let local = 80.0 + 0.5 * (40.0 - stored);
    let photons = 10.0 + 0.25 * (40.0 - stored);
    let nuclear = evolved.audit();
    close(nuclear.released.deposited_j, local, 1e-12);
    close(nuclear.released.radiation_j, photons, 1e-12);
    close(nuclear.released.escaped_j, photons, 1e-12);
    close(nuclear.radiation_deposited_j, 0.75 * photons, 1e-12);
    close(nuclear.radiation_escaped_j, 0.25 * photons, 1e-12);
    close(nuclear.thermal_credit_j, local + 0.75 * photons, 1e-12);
    close(
        audit.external_heat_j,
        nuclear.thermal_credit_j + 15.0,
        1e-11,
    );
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + local / 8.0,
        1e-12,
    );
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + (0.5 * photons + 15.0) / 6.0,
        1e-12,
    );
    close(
        heated.temperature_k(CellId(3)).unwrap(),
        300.0 + 0.25 * photons / 4.0,
        1e-12,
    );
    assert_eq!(audit.substeps, 2);
    assert!(nuclear.radiation_roundoff_j <= budget().energy_tolerance_j);
    assert_eq!(heated.cells()[&CellId(2)].boundary.source_w, 3.0);
    assert_eq!(
        evolved
            .source(CellId(1))
            .unwrap()
            .state
            .cumulative_fissions(),
        10.0
    );
    close(
        evolved.cumulative_released().unwrap().radiation_j,
        photons,
        1e-12,
    );

    let (later, cooled, _) = evolved
        .advance_to(&heated, 9.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(cooled.audit().fissions, 0.0);
    assert!(cooled.audit().radiation_deposited_j > 0.0);
    assert!(later.temperature_k(CellId(3)).unwrap() > heated.temperature_k(CellId(3)).unwrap());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .state
            .cumulative_fissions(),
        0.0
    );
}

#[test]
fn g0_complete_escape_and_self_absorption_keep_disjoint_energy_accounts() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    for absorb in [false, true] {
        let deposits = if absorb { vec![(1, 1.0)] } else { vec![] };
        let receivers = if absorb { vec![receiver(1, 1)] } else { vec![] };
        let heating = FissionHeating::new(&initial, vec![source(1, 1, 1.0, &[], false)], budget())
            .unwrap()
            .with_radiation_response(
                &initial,
                response(vec![row(1, &deposits, if absorb { 0.0 } else { 1.0 })]),
                receivers,
            )
            .unwrap();
        let (heated, evolved, audit) = heating
            .advance_to(&initial, 2.0, controls(2.0), || Ok(()))
            .unwrap();
        let absorbed = if absorb { 2.0 } else { 0.0 };
        assert_eq!(
            evolved.audit().released,
            EnergyChannels {
                deposited_j: 16.0,
                radiation_j: 2.0,
                escaped_j: 2.0,
            }
        );
        assert_eq!(evolved.audit().radiation_deposited_j, absorbed);
        assert_eq!(evolved.audit().radiation_escaped_j, 2.0 - absorbed);
        assert_eq!(evolved.audit().thermal_credit_j, 16.0 + absorbed);
        assert_eq!(audit.external_heat_j, 16.0 + absorbed);
        close(
            heated.temperature_k(CellId(1)).unwrap(),
            300.0 + (16.0 + absorbed) / 8.0,
            1e-12,
        );
    }
}

#[test]
fn g3_shared_receiver_preserves_source_transitions_and_identical_partition_restart() {
    let initial = thermal(
        vec![
            cell(1, 1, 2.0, 1200.0, 0.0),
            cell(2, 2, 3.0, 600.0, 0.0),
            cell(3, 1, 1.0, 1200.0, 2.0),
        ],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = FissionHeating::new(
        &initial,
        vec![
            source(2, 2, 0.0, &[(1.0, 2.0), (3.0, 0.0)], true),
            source(1, 1, 5.0, &[(2.0, 0.0)], true),
        ],
        budget(),
    )
    .unwrap()
    .with_radiation_response(
        &initial,
        response(vec![row(2, &[(3, 0.75)], 0.25), row(1, &[(3, 0.5)], 0.5)]),
        vec![receiver(3, 1)],
    )
    .unwrap();
    let (whole_t, whole_s, whole_audit) = heating
        .advance_to(&initial, 5.0, controls(20.0), || Ok(()))
        .unwrap();
    let (part_t, part_s, part_audit) = heating
        .advance_to(&initial, 2.0, controls(20.0), || Ok(()))
        .unwrap();
    let (split_t, split_s, split_audit) = part_s
        .clone()
        .advance_to(&part_t.clone(), 5.0, controls(20.0), || Ok(()))
        .unwrap();
    assert_eq!(whole_audit.substeps, 4);
    assert_eq!(whole_t.cells(), split_t.cells());
    for id in [CellId(1), CellId(2)] {
        assert_eq!(
            whole_s.source(id).unwrap().state,
            split_s.source(id).unwrap().state
        );
    }
    let photons_1 = whole_s
        .source(CellId(1))
        .unwrap()
        .state
        .cumulative_released()
        .radiation_j;
    let photons_2 = whole_s
        .source(CellId(2))
        .unwrap()
        .state
        .cumulative_released()
        .radiation_j;
    let absorbed = 0.5 * photons_1 + 0.75 * photons_2;
    close(whole_s.audit().radiation_deposited_j, absorbed, 1e-12);
    close(
        whole_t.temperature_k(CellId(3)).unwrap(),
        300.0 + (absorbed + 10.0) / 4.0,
        1e-12,
    );
    close(
        part_s.audit().radiation_deposited_j + split_s.audit().radiation_deposited_j,
        whole_s.audit().radiation_deposited_j,
        1e-12,
    );
    close(
        part_audit.external_heat_j + split_audit.external_heat_j,
        whole_audit.external_heat_j,
        1e-11,
    );
    let (idle_t, idle_s, _) = whole_s
        .advance_to(&whole_t, 5.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(idle_t.cells(), whole_t.cells());
    assert_eq!(idle_s.audit(), FissionHeatingAudit::default());
}

#[test]
fn g0_response_coverage_and_receiver_binding_are_checked_even_on_idle_advance() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating =
        FissionHeating::new(&initial, vec![source(1, 1, 1.0, &[], false)], budget()).unwrap();
    let transport = response(vec![row(1, &[(2, 1.0)], 0.0)]);
    for receivers in [
        vec![],
        vec![receiver(2, 2), receiver(2, 2)],
        vec![receiver(2, 1)],
        vec![receiver(3, 2)],
        vec![receiver(2, 2), receiver(1, 1)],
    ] {
        assert!(
            heating
                .clone()
                .with_radiation_response(&initial, Arc::clone(&transport), receivers)
                .is_err()
        );
    }
    for rows in [
        vec![row(9, &[(2, 1.0)], 0.0)],
        vec![row(1, &[(2, 1.0)], 0.0), row(9, &[], 1.0)],
    ] {
        assert!(
            heating
                .clone()
                .with_radiation_response(&initial, response(rows), vec![receiver(2, 2)])
                .is_err()
        );
    }
    let bound = heating
        .with_radiation_response(&initial, transport, vec![receiver(2, 2)])
        .unwrap();
    for change in 0..5 {
        let mut receiver_cell = cell(2, 2, 3.0, 600.0, 0.0);
        match change {
            0 => receiver_cell.mass_kg = 4.0,
            1 => receiver_cell.geometry.volume_m3 = 0.02,
            2 => receiver_cell.body_id += 1,
            3 => receiver_cell.segment_id += 1,
            _ => receiver_cell.material_id = MaterialId(3),
        }
        let changed = thermal(
            vec![cell(1, 1, 2.0, 1200.0, 0.0), receiver_cell],
            &[
                (MaterialId(1), &HOT),
                (MaterialId(2), &COOL),
                (MaterialId(3), &COOL),
            ],
            vec![],
        );
        assert!(
            bound
                .advance_to(&changed, 0.0, controls(1.0), || Ok(()))
                .is_err()
        );
    }
    let changed_calorics = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &LIMITED)],
        vec![],
    );
    assert!(matches!(
        bound.advance_to(&changed_calorics, 0.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::BindingChanged(CellId(2)))
    ));
    assert_eq!(bound.time_s(), 0.0);
    assert_eq!(
        bound.source(CellId(1)).unwrap().state.cumulative_fissions(),
        0.0
    );
}

#[test]
fn g4_receiver_domain_roundoff_and_late_cancellation_rollback_together() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 1.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &LIMITED)],
        vec![],
    );
    let heating = FissionHeating::new(&initial, vec![source(1, 1, 10.0, &[], false)], budget())
        .unwrap()
        .with_radiation_response(
            &initial,
            response(vec![row(1, &[(2, 1.0)], 0.0)]),
            vec![receiver(2, 2)],
        )
        .unwrap();
    let source_before = heating.source(CellId(1)).unwrap().state.clone();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(heating.source(CellId(1)).unwrap().state, source_before);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(2)).unwrap(), 300.0);

    let mut polls = 0;
    let expected = heating
        .advance_to(&initial, 0.1, controls(0.025), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for cutoff in [polls / 2, polls] {
        let mut calls = 0;
        let result = heating.advance_to(&initial, 0.1, controls(0.025), || {
            calls += 1;
            if calls == cutoff {
                Err(DepositionError::Budget("receiver checkpoint"))
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(FissionThermalError::Thermal(DepositionError::Budget(
                "receiver checkpoint"
            )))
        ));
        assert_eq!(heating.source(CellId(1)).unwrap().state, source_before);
    }
    let retried = heating
        .advance_to(&initial, 0.1, controls(0.025), || Ok(()))
        .unwrap();
    assert_eq!(retried.0.cells(), expected.0.cells());
    assert_eq!(retried.1.audit(), expected.1.audit());
    assert_eq!(
        retried.1.source(CellId(1)).unwrap().state,
        expected.1.source(CellId(1)).unwrap().state
    );
    assert_eq!(retried.2, expected.2);

    // Original source products are exact here, but the authored 0.1/0.9 row
    // has a nonzero floating-point defect. It cannot spend a zero allowance.
    let mut strict = budget();
    strict.energy_tolerance_j = 0.0;
    let strict_heating = FissionHeating::new(&initial, vec![source(1, 1, 3.0, &[], false)], strict)
        .unwrap()
        .with_radiation_response(
            &initial,
            response(vec![row(1, &[(2, 0.1)], 0.9)]),
            vec![receiver(2, 2)],
        )
        .unwrap();
    assert!(
        strict_heating
            .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
            .is_err()
    );
    assert_eq!(
        strict_heating
            .source(CellId(1))
            .unwrap()
            .state
            .cumulative_fissions(),
        0.0
    );
}

#[test]
fn g1_original_radioactive_inventory_uses_the_same_nonlocal_deposition_owner() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let model = Arc::new(
        DecayModel::new(
            "manufactured-radiating-parent-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(Time::new(1.0)).unwrap(),
                },
                DecayNuclide {
                    id: "stable-daughter".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 2.0,
                    radiation_j: 8.0,
                    escaped_j: 1.0,
                },
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let heating = NuclearHeating::new(
        &initial,
        vec![DecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: DecayInventory::new(model, vec![100.0, 0.0]).unwrap(),
        }],
        DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    )
    .unwrap()
    .with_radiation_response(
        &initial,
        response(vec![row(1, &[(2, 0.5)], 0.5)]),
        vec![receiver(2, 2)],
    )
    .unwrap();
    let (heated, decayed, audit) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    close(
        decayed.source(CellId(1)).unwrap().inventory.populations()[0],
        50.0,
        1e-12,
    );
    close(
        decayed.source(CellId(1)).unwrap().inventory.populations()[1],
        50.0,
        1e-12,
    );
    close(decayed.audit().decays, 50.0, 1e-12);
    close(decayed.audit().released.deposited_j, 100.0, 1e-12);
    close(decayed.audit().released.radiation_j, 400.0, 1e-12);
    close(decayed.audit().released.escaped_j, 50.0, 1e-12);
    close(decayed.audit().radiation_deposited_j, 200.0, 1e-12);
    close(decayed.audit().radiation_escaped_j, 200.0, 1e-12);
    close(audit.external_heat_j, 300.0, 1e-12);
    close(heated.temperature_k(CellId(1)).unwrap(), 312.5, 1e-12);
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + 200.0 / 6.0,
        1e-12,
    );
    assert_eq!(heated.cells()[&CellId(2)].mass_kg, 3.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        &[100.0, 0.0]
    );
}

#[test]
fn g1_adaptive_controller_resolves_non_source_receiver_temperature() {
    let mut receiving_cell = cell(2, 2, 3.0, 600.0, 0.0);
    receiving_cell.boundary.htc_w_m2_k = 1.5;
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), receiving_cell],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let mut cooling_source = source(1, 1, 0.0, &[], true);
    cooling_source.state = FissionState::new(model(true), 1000.0, vec![40.0]).unwrap();
    let heating = FissionHeating::new(&initial, vec![cooling_source], budget())
        .unwrap()
        .with_radiation_response(
            &initial,
            response(vec![row(1, &[(2, 1.0)], 0.0)]),
            vec![receiver(2, 2)],
        )
        .unwrap();
    // Receiver C=6 J/K, G=1.5 W/K, radiative input P=5 exp(-t/2) W.
    // Only the receiver's native convective dynamics require time refinement.
    let end = 4.0;
    let exact =
        300.0 + (5.0 / 6.0) * ((-0.5_f64 * end).exp() - (-0.25_f64 * end).exp()) / (0.25 - 0.5);
    let mut errors = Vec::new();
    for tolerance in [0.02, 0.00005] {
        let adaptive = NuclearAdaptiveControls {
            step: controls(4.0),
            initial_dt_s: 4.0,
            minimum_dt_s: 1e-8,
            maximum_trials: 20_000,
            maximum_substeps: 20_000,
            tolerances: NuclearTimeTolerances {
                specific_energy_j_kg: 4.0 * tolerance,
                temperature_k: tolerance,
                nuclei: 1.0,
                nuclear_energy_j: 1.0,
                neutron_population: 1.0,
                neutron_count: 1.0,
                relative: 0.0,
            },
        };
        let (heated, evolved, report) = heating
            .advance_adaptive_to(&initial, end, adaptive, || Ok(()))
            .unwrap();
        errors.push((heated.temperature_k(CellId(2)).unwrap() - exact).abs());
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        close(
            evolved.audit().radiation_deposited_j,
            10.0 * -(-0.5 * end).exp_m1(),
            1e-10,
        );
        close(
            evolved.audit().thermal_credit_j,
            30.0 * -(-0.5 * end).exp_m1(),
            1e-10,
        );
        // The native external account is net of receiver convection. Its
        // storage equals analytic local heat plus the receiver's C * delta T;
        // received radiation beyond that receiver storage leaves to ambient.
        let receiver_storage_j = 6.0 * (heated.temperature_k(CellId(2)).unwrap() - 300.0);
        let convected_j = evolved.audit().radiation_deposited_j - receiver_storage_j;
        assert!(convected_j > 0.0);
        close(
            report.thermal.external_heat_j,
            20.0 * -(-0.5 * end).exp_m1() + receiver_storage_j,
            1e-10,
        );
        assert!(report.thermal.residual_j.abs() < 1e-8);
        assert_eq!(evolved.audit().fissions, 0.0);
    }
    assert!(
        errors[1] < errors[0] / 5.0,
        "independent receiver ODE errors {errors:?}"
    );
    assert!(
        errors[1] < 0.01,
        "independent receiver ODE errors {errors:?}"
    );
}
