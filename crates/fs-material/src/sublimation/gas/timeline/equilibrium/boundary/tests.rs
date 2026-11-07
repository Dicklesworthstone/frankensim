//! Original clock, fusion, stream and whole-network conservation oracles.
use super::*;
use crate::sublimation::gas::{
    tests as fixture,
    timeline::{advance_deposition_with_equilibrium, equilibrium_tests as native},
};
use fs_thermal::enthalpy::{CellId, DepositionBudget, Event, MaterialId, ThermalLink, TimedEvent};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn reflection(catalog: &ThermochemicalGasModel) -> ThermalAccommodation {
    ThermalAccommodation::new(
        catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.5),
            (SpeciesId::new("B").unwrap(), 0.5),
        ],
        [200.0, 900.0],
    )
    .unwrap()
}

#[test]
fn no_event_retains_the_original_singleton_heat_and_phase_endpoint() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let curve = native::curve(&catalog, &binding, native::FUSION);
    let phase = fixture::phase();
    let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * native::FUSION;
    let thermal = native::schedule(&curve, vec![native::birth(0.0, native::cell(1, h, 0.1))]);
    let gas = fixture::gas(&catalog, 500.0, 5000.0);
    let reflection = reflection(&catalog);
    let laws = BTreeMap::from([(CellId(1), reflection.clone())]);
    for target in [0.0, 1e-5] {
        let selected = native::selection(&curve, &binding, &phase, 1);
        let ordinary = advance_deposition_with_equilibrium(
            &thermal,
            &[selected],
            Some(&laws),
            &gas,
            &catalog,
            target,
            native::heat(1.0),
            native::controls(),
            || Ok(()),
        )
        .unwrap();
        let candidate = advance_deposition_to_melting(
            &thermal,
            selected,
            Some(&reflection),
            &gas,
            &catalog,
            target,
            native::heat(1.0),
            native::controls(),
            || Ok(()),
        )
        .unwrap();
        let actual = candidate.exchange();
        assert!(!candidate.boundary().melting_completed);
        assert_eq!(candidate.boundary().elapsed_duration_s, target);
        assert_eq!(candidate.boundary().remaining_duration_s, 0.0);
        assert_eq!(actual.thermal().cells(), ordinary.thermal().cells());
        assert_eq!(actual.gas(), ordinary.gas());
        assert_eq!(actual.audit().thermal, ordinary.audit().thermal);
        assert_eq!(
            actual.audit().condensed_phases,
            ordinary.audit().condensed_phases
        );
        fixture::near(
            actual.audit().energy_residual_j,
            ordinary.audit().energy_residual_j,
            1e-12,
        );
        assert_eq!(
            actual.audit().reflected_heat_into_gas_j,
            ordinary.audit().reflected_heat_into_gas_j
        );
    }
    assert!(thermal.cells().is_empty());
}

#[test]
fn deposition_and_reflection_each_finish_on_the_same_native_neighbor_clock() {
    for reflected in [false, true] {
        let catalog = fixture::catalog(0.0);
        let mut binding = fixture::binding(0.0);
        if reflected {
            binding.accommodation_coefficient = 0.0;
        }
        let curve = native::curve(&catalog, &binding, native::FUSION);
        let phase = fixture::phase();
        let selected = native::selection(&curve, &binding, &phase, 1);
        let liquid_h = selected.caloric.liquid_at_melting_h_j_kg();
        let initial_h = liquid_h - 100.0;
        let power = 1.0;
        let thermal = native::schedule(
            &curve,
            vec![
                native::birth(0.0, native::cell(1, initial_h, power)),
                native::birth(0.0, native::cell(2, initial_h, 2.0)),
            ],
        )
        .advance_to(0.0, native::heat(1.0), || Ok(()))
        .unwrap()
        .0;
        let gas = fixture::gas(
            &catalog,
            if reflected { 500.0 } else { 400.0 },
            if reflected { 0.0 } else { 100000.0 },
        );
        let reflection = reflection(&catalog);
        let mass_rate = if reflected {
            0.0
        } else {
            0.1 * (1000.0 - 100000.0) * 1e-4
                / (2.0 * core::f64::consts::PI * (R / 0.018015) * 400.0).sqrt()
        };
        let reflection_rate = if reflected {
            0.5 * 100000.0 * 1e-4 / (2.0 * core::f64::consts::PI * (R / 0.02897) * 500.0).sqrt()
                * (3.0 * R / 0.02897)
                * (400.0 - 500.0)
        } else {
            0.0
        };
        // The source liquid endpoint includes L_f. Incoming solid must also
        // melt, so deposition supplies L_sub - R*T/2 - L_f per unit mass.
        let phase_drive = mass_rate * (native::FUSION - 200000.0 + 0.5 * (R / 0.018015) * 400.0)
            - reflection_rate;
        let expected_time = 0.01 * (liquid_h - initial_h) / (power + phase_drive);
        let target = 4.0 * expected_time;
        let candidate = advance_deposition_to_melting(
            &thermal,
            selected,
            if reflected { Some(&reflection) } else { None },
            &gas,
            &catalog,
            target,
            native::heat(target),
            native::controls(),
            || Ok(()),
        )
        .unwrap();
        let accepted = candidate.exchange();
        let audit = accepted.audit();
        assert!(candidate.boundary().melting_completed);
        fixture::near(
            accepted.thermal().time_s(),
            expected_time,
            expected_time * 1e-8,
        );
        assert_eq!(
            candidate.boundary().remaining_duration_s,
            target - accepted.thermal().time_s()
        );
        let final_cell = &accepted.thermal().cells()[&CellId(1)];
        assert_eq!(final_cell.specific_h_j_kg, liquid_h);
        assert_eq!(
            audit.condensed_phases[&CellId(1)]
                .final_inventory
                .solid_mass_kg,
            0.0
        );
        fixture::near(audit.transferred_mass_kg, mass_rate * expected_time, 1e-12);
        fixture::near(final_cell.mass_kg, 0.01 - mass_rate * expected_time, 1e-12);
        fixture::near(
            0.01 * (accepted.thermal().cells()[&CellId(2)].specific_h_j_kg - initial_h),
            2.0 * accepted.thermal().time_s(),
            1e-9,
        );
        fixture::near(
            audit.thermal.external_heat_j,
            3.0 * accepted.thermal().time_s(),
            1e-12,
        );
        fixture::near(
            audit.reflected_heat_into_gas_j,
            reflection_rate * expected_time,
            1e-8,
        );
        fixture::near(audit.mass_residual_kg, 0.0, 1e-14);
        fixture::near(audit.energy_residual_j, 0.0, 1e-7);
        fixture::near(
            audit.solid_energy_change_j + audit.gas_energy_change_j - audit.thermal.external_heat_j,
            0.0,
            1e-7,
        );
        assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial_h);
        assert!(
            advance_deposition_with_equilibrium(
                &thermal,
                &[selected],
                None,
                &gas,
                &catalog,
                target,
                native::heat(target),
                native::controls(),
                || Ok(())
            )
            .is_err()
                || reflected
        );
    }
}

#[test]
fn saturated_gross_turnover_preserves_zero_net_exchange_at_native_melting() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = native::curve(&catalog, &binding, native::FUSION);
    let phase = fixture::phase();
    let selected = native::selection(&curve, &binding, &phase, 1);
    let liquid_h = selected.caloric.liquid_at_melting_h_j_kg();
    let thermal = native::schedule(
        &curve,
        vec![native::birth(0.0, native::cell(1, liquid_h - 100.0, 4.0))],
    )
    .advance_to(0.0, native::heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let reflection = reflection(&catalog);
    let candidate = advance_deposition_to_melting(
        &thermal,
        selected,
        Some(&reflection),
        &gas,
        &catalog,
        1.0,
        native::heat(1.0),
        native::controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = candidate.exchange().audit();
    assert!(candidate.boundary().melting_completed);
    fixture::near(candidate.boundary().elapsed_duration_s, 0.25, 1e-8);
    assert_eq!(candidate.exchange().gas(), &gas);
    assert_eq!(
        candidate.exchange().thermal().cells()[&CellId(1)].mass_kg,
        0.01
    );
    assert_eq!(
        candidate.exchange().thermal().cells()[&CellId(1)].specific_h_j_kg,
        liquid_h
    );
    assert!(audit.outgoing_mass_kg > 0.0);
    assert_eq!(audit.incoming_mass_kg, audit.outgoing_mass_kg);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.energy_into_gas_j, 0.0);
    assert_eq!(audit.reflected_heat_into_gas_j, 0.0);
    fixture::near(audit.thermal.external_heat_j, 1.0, 1e-7);
    fixture::near(audit.energy_residual_j, 0.0, 1e-7);
    let wrong = GasSurfaceBinding {
        solid_enthalpy_offset_j_kg: 1000.0,
        ..binding.clone()
    };
    assert!(matches!(
        advance_deposition_to_melting(
            &thermal,
            native::selection(&curve, &wrong, &phase, 1),
            Some(&reflection),
            &gas,
            &catalog,
            1.0,
            native::heat(1.0),
            native::controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Exchange(
            GasSurfaceError::CaloricMismatch { .. }
        ))
    ));
}

#[test]
fn delayed_birth_covered_melting_and_same_time_events_retain_native_order() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = native::curve(&catalog, &binding, native::FUSION);
    let phase = fixture::phase();
    let selected = native::selection(&curve, &binding, &phase, 1);
    let liquid_h = selected.caloric.liquid_at_melting_h_j_kg();
    let initial_h = liquid_h - 100.0;
    let reversed = native::cell(1, initial_h, -4.0).boundary;
    let thermal = native::schedule(
        &curve,
        vec![
            native::birth(0.0, native::cell(2, initial_h, 2.0)),
            native::birth(0.125, native::cell(1, initial_h, 4.0)),
            TimedEvent {
                time_s: 0.25,
                event: Event::Contact(ThermalLink {
                    id: 8,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-4,
                    conductance_w_k: 0.0,
                }),
            },
            TimedEvent {
                time_s: 0.375,
                event: Event::Boundary(CellId(1), reversed),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Separation(8),
            },
        ],
    );
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let candidate = advance_deposition_to_melting(
        &thermal,
        selected,
        None,
        &gas,
        &catalog,
        1.0,
        native::heat(1.0),
        native::controls(),
        || Ok(()),
    )
    .unwrap();
    let accepted = candidate.exchange();
    assert!(candidate.boundary().melting_completed);
    fixture::near(accepted.thermal().time_s(), 0.375, 1e-9);
    assert_eq!(
        accepted.thermal().cells()[&CellId(1)].boundary.source_w,
        -4.0
    );
    assert_eq!(accepted.thermal().pending_events()[0].time_s, 0.5);
    assert_eq!(accepted.thermal().exposed_area_m2(CellId(1)).unwrap(), 0.0);
    assert_eq!(accepted.gas(), &gas);
    assert_eq!(accepted.audit().transferred_mass_kg, 0.0);
    fixture::near(accepted.audit().thermal.external_heat_j, 1.75, 1e-7);
    fixture::near(
        accepted.audit().condensed_phases[&CellId(1)].melted_mass_kg,
        1.0 / native::FUSION,
        1e-14,
    );
    fixture::near(accepted.audit().energy_residual_j, 0.0, 1e-7);
    assert!(thermal.cells().is_empty());
}

#[test]
fn rejected_native_probes_spend_the_original_shared_work_budget() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = native::curve(&catalog, &binding, native::FUSION);
    let phase = fixture::phase();
    let selected = native::selection(&curve, &binding, &phase, 1);
    let initial_h = selected.caloric.liquid_at_melting_h_j_kg() - 100.0;
    let thermal = DepositionSimulation::new(
        vec![native::birth(0.0, native::cell(1, initial_h, 4.0))],
        &[(MaterialId(1), &curve)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 16,
            max_substeps: 1,
        },
    )
    .unwrap()
    .advance_to(0.0, native::heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    assert!(matches!(
        advance_deposition_to_melting(
            &thermal,
            selected,
            None,
            &gas,
            &catalog,
            1.0,
            native::heat(1.0),
            native::controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial_h);
}

#[test]
fn original_identity_solid_seed_and_late_cancellation_remain_transactional() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = native::curve(&catalog, &binding, native::FUSION);
    let wrong_curve = native::curve(&catalog, &binding, 1.1 * native::FUSION);
    let phase = fixture::phase();
    let selected = native::selection(&curve, &binding, &phase, 1);
    let liquid_h = selected.caloric.liquid_at_melting_h_j_kg();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let scheduled = native::schedule(
        &curve,
        vec![native::birth(1.0, native::cell(1, liquid_h - 100.0, 0.0))],
    );
    assert!(matches!(
        advance_deposition_to_melting(
            &scheduled,
            native::selection(&wrong_curve, &binding, &phase, 1),
            None,
            &gas,
            &catalog,
            0.0,
            native::heat(1.0),
            native::controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Thermal(DepositionError::InvalidInput(_)))
    ));
    let molten = native::schedule(
        &curve,
        vec![native::birth(1.0, native::cell(1, liquid_h, 0.0))],
    );
    assert!(matches!(
        advance_deposition_to_melting(
            &molten,
            selected,
            None,
            &gas,
            &catalog,
            0.0,
            native::heat(1.0),
            native::controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Thermal(DepositionError::MaterialDomain))
    ));
    let thermal = native::schedule(
        &curve,
        vec![native::birth(0.0, native::cell(1, liquid_h - 100.0, 4.0))],
    );
    let mut calls = 0;
    advance_deposition_to_melting(
        &thermal,
        selected,
        None,
        &gas,
        &catalog,
        1.0,
        native::heat(1.0),
        native::controls(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut replay = 0;
    let refused = advance_deposition_to_melting(
        &thermal,
        selected,
        None,
        &gas,
        &catalog,
        1.0,
        native::heat(1.0),
        native::controls(),
        || {
            replay += 1;
            if replay == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(GasTimelineError::Thermal(DepositionError::Cancelled))
    ));
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.time_s(), 0.0);
}
