//! Original fusion inventories coupled through one finite, conserved atmosphere.
use super::*;
use crate::evaporation::gas::equilibrium::{
    advance_equilibrium_surface,
    tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
};
use crate::evaporation::gas::tests::{Phase, binding, controls, gas};
use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::evaporation::{GasFilm, surface_flux};

fn surface<'a>(
    id: u64,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    binding: &'a GasLiquidBinding,
    bulk: LiquidSurfaceState,
) -> EquilibriumChamberSurface<'a> {
    EquilibriumChamberSurface {
        id,
        law: &Phase,
        binding,
        bulk,
        inventory: EquilibriumEvaporation::new(chart).unwrap(),
        mass_transfer_coefficient_m_s: TRANSFER,
    }
}

fn solve() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 8,
    }
}

#[test]
fn singleton_matches_the_original_phase_owner_and_empty_or_zero_time_preserves_storage() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let body = bulk(&chart, 0.5);
    for (temperature, pressure) in [(MELTING, 0.0), (320.0, 5000.0)] {
        let initial = gas(&source, temperature, pressure, 0.01);
        let old = advance_equilibrium_surface(
            EquilibriumEvaporation::new(&chart).unwrap(),
            &Phase,
            &binding,
            body,
            &initial,
            &source,
            TRANSFER,
            0.5,
            controls(),
            || false,
        )
        .unwrap();
        let joint = advance_chamber_with_gas(
            &[surface(9, &chart, &binding, body)],
            &initial,
            &source,
            0.5,
            solve(),
            || false,
        )
        .unwrap();
        let local = &joint.surfaces()[0];
        assert_eq!(local.id, 9);
        assert_eq!(local.liquid, old.exchange().liquid());
        assert_eq!(joint.gas_exchange(), old.exchange().gas_exchange());
        assert_eq!(joint.phase_inventory(local.id), Some(old.phase()));
        assert_eq!(
            local.transferred_mass_kg,
            old.exchange().step().transferred_mass_kg
        );
        assert_eq!(local.latent_energy_j, old.exchange().step().latent_energy_j);
        assert_eq!(
            local.liquid_enthalpy_change_j,
            old.exchange().step().liquid_enthalpy_change_j
        );
        assert_eq!(
            joint.step().gas_internal_energy_change_j,
            old.exchange().step().gas_internal_energy_change_j
        );
        assert!(std::ptr::eq(joint.catalog(), &source));

        let empty =
            advance_chamber_with_gas(&[], &initial, &source, 0.5, solve(), || false).unwrap();
        assert_eq!(empty.gas(), &initial);
        assert!(empty.surfaces().is_empty());
        assert!(empty.step().species.is_empty());
        assert_eq!(empty.step().transferred_mass_kg, 0.0);
        assert_eq!(empty.step().sum_absolute_energy_residual_j, 0.0);
        let zero = advance_chamber_with_gas(
            &[surface(9, &chart, &binding, body)],
            &initial,
            &source,
            0.0,
            solve(),
            || false,
        )
        .unwrap();
        assert_eq!(zero.gas(), &initial);
        assert_eq!(zero.surfaces()[0].liquid, body);
        let phase = zero.phase_inventory(9).unwrap();
        assert_eq!(phase.initial, phase.final_state);
        assert_eq!(phase.freezing_mass_kg, 0.0);
    }
}

#[test]
fn unequal_mixed_bodies_keep_independent_fusion_balances_and_one_extensive_budget() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let first = bulk(&chart, 0.25);
    let second = LiquidSurfaceState {
        mass_kg: 0.02,
        area_m2: 0.003,
        ..bulk(&chart, 0.75)
    };
    for (temperature, pressure) in [(MELTING, 0.0), (320.0, 5000.0)] {
        let initial = gas(&source, temperature, pressure, 0.01);
        let joint = advance_chamber_with_gas(
            &[
                surface(7, &chart, &binding, second),
                surface(2, &chart, &binding, first),
            ],
            &initial,
            &source,
            0.5,
            solve(),
            || false,
        )
        .unwrap();
        let mut final_energy =
            joint.gas().mass_kg() * joint.gas().state().specific_internal_energy_j_kg();
        let mut final_mass = joint.gas().mass_kg();
        let mut absolute_mass_defect = joint.gas_exchange().step().mass_residual_kg.abs();
        for (local, original) in joint.surfaces().iter().zip([first, second]) {
            let phase = joint.phase_inventory(local.id).unwrap();
            let dm = local.transferred_mass_kg;
            assert_eq!(dm.is_sign_positive(), pressure == 0.0);
            // On the original plateau, signed freezing pays for vaporization
            // and, on condensation, the incoming vapor's upstream sensible heat.
            let frozen = dm
                * (VAPORIZATION
                    + if dm < 0.0 {
                        CP * (temperature - MELTING)
                    } else {
                        0.0
                    })
                / FUSION;
            near(phase.freezing_mass_kg, frozen, 4e-17);
            near(
                phase.final_state.solid_mass_kg,
                phase.initial.solid_mass_kg + frozen,
                4e-17,
            );
            near(
                phase.final_state.liquid_mass_kg,
                phase.initial.liquid_mass_kg - dm - frozen,
                4e-17,
            );
            near(phase.fusion_energy_j, frozen * FUSION, 4e-12);
            assert_eq!(local.surface_temperature_k, [MELTING, MELTING]);
            assert_eq!(phase.initial.mass_kg, original.mass_kg);
            assert_eq!(
                phase.final_state.equilibrium.phase_curve_identity(),
                chart.identity()
            );
            final_energy += local.liquid.mass_kg * local.liquid.specific_enthalpy_j_kg;
            final_mass += local.liquid.mass_kg;
            absolute_mass_defect +=
                local.mass_residual_kg.abs() + phase.liquid_partition_residual_kg.abs();
        }
        let original_energy = initial.mass_kg() * initial.state().specific_internal_energy_j_kg()
            + first.mass_kg * first.specific_enthalpy_j_kg
            + second.mass_kg * second.specific_enthalpy_j_kg;
        near(final_energy, original_energy, 1e-7);
        near(
            final_mass,
            initial.mass_kg() + first.mass_kg + second.mass_kg,
            controls().gas.mass_tolerance_kg,
        );
        near(
            joint.step().liquid_enthalpy_change_j + joint.step().gas_internal_energy_change_j,
            0.0,
            controls().gas.energy_tolerance_j,
        );
        assert!(joint.step().sum_absolute_energy_residual_j <= controls().gas.energy_tolerance_j);
        assert!(absolute_mass_defect <= controls().gas.mass_tolerance_kg);
        assert_eq!(joint.step().species.len(), 1);
    }
}

#[test]
fn opposing_streams_keep_fusion_and_sensible_heat_with_exact_canonical_mass_selection() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let initial = gas(&source, 320.0, 1600.0, 0.01);
    let hot = LiquidSurfaceState {
        specific_enthalpy_j_kg: CP * 360.0 - VAPORIZATION,
        ..bulk(&chart, 1.0)
    };
    let mut cold = bulk(&chart, 0.5);
    let duration = 0.2;
    let film = GasFilm {
        total_pressure_pa: initial.state().state().pressure,
        vapor_partial_pressure_pa: initial
            .species_partial_pressure_pa(&binding.species)
            .unwrap(),
        temperature_k: 320.0,
        mass_transfer_coefficient_m_s: TRANSFER,
    };
    let target =
        hot.area_m2 * duration * surface_flux(&Phase, 360.0, film).unwrap().mass_flux_kg_m2_s;
    cold.area_m2 = target
        / (-surface_flux(&Phase, MELTING, film)
            .unwrap()
            .mass_flux_kg_m2_s
            * duration);
    let cold_surface = surface(3, &chart, &binding, cold);
    let hot_surface = surface(8, &chart, &binding, hot);
    let mut selected_ids = Vec::new();
    let selected = advance_chamber_with_gas_mass_selection(
        &[hot_surface, cold_surface],
        &initial,
        &source,
        duration,
        solve(),
        |id, requested| {
            selected_ids.push(id);
            let represented = if id == 3 { -target } else { target };
            near(represented, requested, controls().gas.mass_tolerance_kg);
            Ok(represented)
        },
        || false,
    )
    .unwrap();
    assert_eq!(selected_ids, [3, 8]);
    let reordered = advance_chamber_with_gas_mass_selection(
        &[cold_surface, hot_surface],
        &initial,
        &source,
        duration,
        solve(),
        |id, _| Ok(if id == 3 { -target } else { target }),
        || false,
    )
    .unwrap();
    assert_eq!(selected, reordered);
    assert_eq!(selected.step().transferred_mass_kg, 0.0);
    assert_eq!(selected.step().species[0].transferred_mass_kg, 0.0);
    assert_eq!(
        selected.step().evaporated_mass_kg,
        selected.step().condensed_mass_kg
    );
    assert!(selected.step().evaporated_mass_kg > 0.0);
    near(
        selected.step().energy_into_gas_j,
        target * CP * (360.0 - 320.0),
        1e-12,
    );
    assert!(selected.step().gas_temperature_k[1] > 320.0);
    near(
        selected.gas().species_mass_kg(&binding.species).unwrap(),
        initial.species_mass_kg(&binding.species).unwrap(),
        controls().gas.mass_tolerance_kg,
    );
    let phase = selected.phase_inventory(3).unwrap();
    near(
        phase.freezing_mass_kg,
        -target * (VAPORIZATION + CP * 20.0) / FUSION,
        2e-17,
    );
    assert!(phase.freezing_mass_kg < 0.0);
    assert_eq!(selected.phase_inventory(8).unwrap().freezing_mass_kg, 0.0);
    assert!(selected.step().sum_absolute_energy_residual_j <= controls().gas.energy_tolerance_j);
}

#[test]
fn large_solid_and_other_body_inventories_cannot_finance_a_local_liquid_withdrawal() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let initial = gas(&source, MELTING, 0.0, 0.01);
    let trace_liquid = surface(1, &chart, &binding, bulk(&chart, 1e-6));
    let ample_liquid = surface(2, &chart, &binding, bulk(&chart, 0.75));
    let error = advance_chamber_with_gas(
        &[ample_liquid, trace_liquid],
        &initial,
        &source,
        1.0,
        solve(),
        || false,
    )
    .unwrap_err();
    let GasEvaporationError::StepSizeLimit {
        maximum_duration_s, ..
    } = error
    else {
        panic!("the original local liquid inventory must bound evaporation: {error:?}");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < 0.1);
    let accepted = advance_chamber_with_gas(
        &[ample_liquid, trace_liquid],
        &initial,
        &source,
        0.99 * maximum_duration_s,
        solve(),
        || false,
    )
    .unwrap();
    let local = &accepted.surfaces()[0];
    let phase = accepted.phase_inventory(local.id).unwrap();
    assert!(phase.initial.liquid_mass_kg < phase.initial.mass_kg * 1e-5);
    assert!(
        local.transferred_mass_kg
            <= controls().maximum_inventory_fraction * phase.initial.liquid_mass_kg
    );
    assert!(phase.final_state.liquid_mass_kg > 0.0);
    assert!(phase.freezing_mass_kg > 0.0);
}

#[test]
fn simultaneous_evaporation_cannot_fund_excess_gross_condensation_from_the_old_gas() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let initial = gas(&source, MELTING, 2500.0, 1e-6);
    let before = initial.clone();
    let cold = bulk(&chart, 0.5);
    let old_vapor = initial.species_mass_kg(&binding.species).unwrap();
    let film = GasFilm {
        total_pressure_pa: initial.state().state().pressure,
        vapor_partial_pressure_pa: initial
            .species_partial_pressure_pa(&binding.species)
            .unwrap(),
        temperature_k: MELTING,
        mass_transfer_coefficient_m_s: TRANSFER,
    };
    let cold_flux = surface_flux(&Phase, MELTING, film)
        .unwrap()
        .mass_flux_kg_m2_s;
    let target = 0.2 * old_vapor;
    let duration = target / (-cold_flux * cold.area_m2);
    let hot_flux = surface_flux(&Phase, 600.0, film).unwrap().mass_flux_kg_m2_s;
    assert!(cold_flux < 0.0 && hot_flux > 0.0);
    let hot = LiquidSurfaceState {
        specific_enthalpy_j_kg: CP * 600.0 - VAPORIZATION,
        area_m2: target / (hot_flux * duration),
        ..bulk(&chart, 1.0)
    };
    // Either mixed body alone withdraws 20% of old vapor and stays condensing.
    advance_chamber_with_gas(
        &[surface(1, &chart, &binding, cold)],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap();
    // Two withdraw 40%, while the hot body returns 20%. The 20% net debit
    // cannot conceal a violation of the original 25% gross donor allowance.
    let error = advance_chamber_with_gas(
        &[
            surface(3, &chart, &binding, hot),
            surface(2, &chart, &binding, cold),
            surface(1, &chart, &binding, cold),
        ],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            GasEvaporationError::StepSizeLimit {
                reason: "shared initial vapor inventory withdrawal limit",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(initial, before);
}

#[test]
fn phase_source_refusal_and_late_cancellation_never_publish_partial_bodies() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let chart = curve();
    let initial = gas(&source, MELTING, 0.0, 0.01);
    let before = initial.clone();
    let first = surface(1, &chart, &binding, bulk(&chart, 0.25));
    let second = surface(2, &chart, &binding, bulk(&chart, 0.75));
    let mut selections = 0;
    assert!(
        advance_chamber_with_gas_mass_selection(
            &[first, first],
            &initial,
            &source,
            0.2,
            solve(),
            |_, dm| {
                selections += 1;
                Ok(dm)
            },
            || false,
        )
        .is_err()
    );
    assert_eq!(selections, 0);
    let wrong_binding = GasLiquidBinding {
        liquid_enthalpy_offset_j_kg: 1000.0,
        ..binding.clone()
    };
    assert!(matches!(
        advance_chamber_with_gas(
            &[first, surface(2, &chart, &wrong_binding, second.bulk)],
            &initial,
            &source,
            0.0,
            solve(),
            || false,
        ),
        Err(GasEvaporationError::CaloricMismatch { .. })
    ));
    let mut calls = 0;
    advance_chamber_with_gas(&[first, second], &initial, &source, 0.2, solve(), || {
        calls += 1;
        false
    })
    .unwrap();
    let final_gate = calls;
    let mut calls = 0;
    assert_eq!(
        advance_chamber_with_gas(&[first, second], &initial, &source, 0.2, solve(), || {
            calls += 1;
            calls == final_gate
        })
        .unwrap_err(),
        GasEvaporationError::Cancelled
    );
    assert_eq!(initial, before);
    assert_eq!(first.bulk, bulk(&chart, 0.25));
    assert_eq!(second.bulk, bulk(&chart, 0.75));
}
