//! Original-source heat, fusion and molecular deposition in one accepted history.
use super::*;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use crate::sublimation::gas::tests as fixture;
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositionBudget, MaterialId, ThermalLink, TimedEvent,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const FUSION: f64 = 80_000.0;

fn curve(
    catalog: &ThermochemicalGasModel,
    binding: &GasSurfaceBinding,
    fusion: f64,
) -> EquilibriumEnthalpyPhaseCurve {
    let (h, cp) = fixture::chart(catalog, binding);
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([92; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h - 200.0 * cp,
                temperature_k: 200.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: 400.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h + fusion,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h + fusion + 100.0 * cp,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}

fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn controls() -> GasChamberSolve {
    GasChamberSolve {
        surface: fixture::controls(),
        maximum_surfaces: 2,
    }
}
fn cell(id: u64, h: f64, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: 0.01,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 917.0,
            surface_area_m2: 1e-4,
        },
        boundary: Boundary {
            ambient_k: 400.0,
            htc_w_m2_k: 0.0,
            source_w: power,
            reservoir: None,
        },
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn schedule(
    curve: &EquilibriumEnthalpyPhaseCurve,
    events: Vec<TimedEvent>,
) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), curve)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 16,
            max_substeps: 20000,
        },
    )
    .unwrap()
}
fn selection<'a>(
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    binding: &'a GasSurfaceBinding,
    phase: &'a dyn SublimationLaw,
    id: u64,
) -> GasTimelineEquilibriumSurface<'a> {
    GasTimelineEquilibriumSurface {
        cell_id: CellId(id),
        law: phase,
        binding,
        caloric: SharpMeltingCaloric::new(curve).unwrap(),
    }
}

#[test]
fn deposition_and_cooling_share_original_fusion_inventory_without_extra_latent_heat() {
    for datum in [0.0, -1e7] {
        let catalog = fixture::catalog(datum);
        let binding = fixture::binding(datum);
        let curve = curve(&catalog, &binding, FUSION);
        let phase = fixture::phase();
        let selected = [selection(&curve, &binding, &phase, 1)];
        let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * FUSION;
        let thermal = schedule(&curve, vec![birth(0.0, cell(1, h, -20.0))]);
        let gas = fixture::gas(&catalog, 400.0, 5000.0);
        let candidate = advance_deposition_with_equilibrium(
            &thermal,
            &selected,
            None,
            &gas,
            &catalog,
            0.01,
            heat(0.01),
            controls(),
            || Ok(()),
        )
        .unwrap();
        let audit = candidate.audit();
        let partition = &audit.condensed_phases[&CellId(1)];
        let dm = audit.cells[&CellId(1)].transferred_mass_kg;
        let effective_latent = 200_000.0 - 0.5 * (R / 0.018015) * 400.0;
        let expected_liquid_change = (-0.2 - dm * effective_latent) / FUSION;
        assert!(dm < 0.0);
        assert!(expected_liquid_change < 0.0);
        fixture::near(partition.melted_mass_kg, expected_liquid_change, 1e-14);
        fixture::near(partition.initial.liquid_mass_kg, 0.005, 1e-14);
        assert_eq!(partition.final_inventory.caloric.temperature_k(), 400.0);
        assert_eq!(
            partition.final_inventory.caloric.phase_curve_identity(),
            curve.identity()
        );
        fixture::near(audit.thermal.external_heat_j, -0.2, 1e-12);
        fixture::near(audit.mass_residual_kg, 0.0, 1e-14);
        fixture::near(audit.energy_residual_j, 0.0, 1e-7);
        fixture::near(
            audit.solid_energy_change_j + audit.gas_energy_change_j
                - audit.incoming_energy_j
                - audit.thermal.external_heat_j,
            0.0,
            1e-7,
        );
        assert!(thermal.cells().is_empty());
    }
}

#[test]
fn contact_conduction_freezes_original_plateau_and_preserves_the_cold_cell_source() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = curve(&catalog, &binding, FUSION);
    let phase = fixture::phase();
    let selected = [
        selection(&curve, &binding, &phase, 1),
        selection(&curve, &binding, &phase, 2),
    ];
    let h = curve.knots()[1].specific_enthalpy_j_kg;
    let cp = fixture::chart(&catalog, &binding).1;
    let thermal = schedule(
        &curve,
        vec![
            birth(0.0, cell(1, h + 0.5 * FUSION, 0.0)),
            birth(0.0, cell(2, h - 100.0 * cp, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 8,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-5,
                    conductance_w_k: 0.05,
                }),
            },
        ],
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let candidate = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = candidate.audit();
    fixture::near(
        audit.condensed_phases[&CellId(1)].melted_mass_kg,
        -0.05 / FUSION,
        1e-14,
    );
    assert_eq!(audit.condensed_phases[&CellId(2)].melted_mass_kg, 0.0);
    fixture::near(
        candidate.thermal().temperature_k(CellId(2)).unwrap(),
        300.0 + 0.05 / (0.01 * cp),
        1e-10,
    );
    assert_eq!(candidate.gas(), &gas);
    assert_eq!(audit.energy_into_gas_j, 0.0);
    fixture::near(audit.energy_residual_j, 0.0, 1e-7);
}

#[test]
fn heat_crosses_solid_fusion_onset_with_the_same_original_caloric_chart() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = curve(&catalog, &binding, FUSION);
    let phase = fixture::phase();
    let selected = [selection(&curve, &binding, &phase, 1)];
    let h = curve.knots()[1].specific_enthalpy_j_kg;
    let thermal = schedule(&curve, vec![birth(0.0, cell(1, h - 100.0, 200.0))]);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let candidate = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let partition = &candidate.audit().condensed_phases[&CellId(1)];
    assert_eq!(partition.initial.liquid_mass_kg, 0.0);
    fixture::near(partition.melted_mass_kg, 1.0 / FUSION, 1e-14);
    assert_eq!(partition.final_inventory.caloric.temperature_k(), 400.0);
    fixture::near(candidate.audit().thermal.external_heat_j, 2.0, 1e-12);
}

#[test]
fn reflected_gas_heat_freezes_bulk_liquid_once_without_transferring_species_mass() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = curve(&catalog, &binding, FUSION);
    let phase = fixture::phase();
    let selected = [selection(&curve, &binding, &phase, 1)];
    let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * FUSION;
    let thermal = schedule(&curve, vec![birth(0.0, cell(1, h, 2.0))]);
    let gas = fixture::gas(&catalog, 300.0, 1000.0);
    let reflection = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.5),
            (fixture::species(), 0.5),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let laws = BTreeMap::from([(CellId(1), reflection)]);
    let trial = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        Some(&laws),
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = trial.audit();
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.incoming_mass_kg, 0.0);
    assert_eq!(audit.outgoing_mass_kg, 0.0);
    assert!(audit.reflected_heat_into_gas_j > 0.002);
    assert_eq!(audit.reflected_heat_into_gas_j, audit.energy_into_gas_j);
    fixture::near(
        audit.condensed_phases[&CellId(1)].melted_mass_kg,
        (0.002 - audit.reflected_heat_into_gas_j) / FUSION,
        1e-14,
    );
    assert_eq!(
        audit.condensed_phases[&CellId(1)]
            .final_inventory
            .caloric
            .temperature_k(),
        400.0
    );
    assert!(trial.gas().state().state().temperature > 300.0);
    fixture::near(audit.energy_residual_j, 0.0, 1e-7);
}

#[test]
fn wrong_latent_identity_and_complete_melting_refuse_without_changing_history() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = curve(&catalog, &binding, FUSION);
    let other = self::curve(&catalog, &binding, 2.0 * FUSION);
    let phase = fixture::phase();
    let selected = [selection(&curve, &binding, &phase, 1)];
    let wrong = [selection(&other, &binding, &phase, 1)];
    let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * FUSION;
    let thermal = schedule(&curve, vec![birth(0.0, cell(1, h, 8000.0))]);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    assert!(matches!(
        advance_deposition_with_equilibrium(
            &thermal,
            &wrong,
            None,
            &gas,
            &catalog,
            0.0,
            heat(0.01),
            controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Thermal(DepositionError::InvalidInput(
            "partitioned timeline caloric identity mismatch"
        )))
    ));
    assert!(matches!(
        advance_deposition_with_equilibrium(
            &thermal,
            &selected,
            None,
            &gas,
            &catalog,
            0.1,
            heat(0.025),
            controls(),
            || Ok(())
        ),
        Err(GasTimelineError::Thermal(DepositionError::MaterialDomain))
    ));
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.pending_events().len(), 1);
}

#[test]
fn late_birth_and_cancellation_keep_phase_inventories_transactional() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = curve(&catalog, &binding, FUSION);
    let phase = fixture::phase();
    let selected = [selection(&curve, &binding, &phase, 1)];
    let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * FUSION;
    let thermal = schedule(&curve, vec![birth(0.05, cell(1, h, 0.0))]);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let first = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.05,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(first.audit().transferred_mass_kg, 0.0);
    assert_eq!(
        first.audit().condensed_phases[&CellId(1)].melted_mass_kg,
        0.0
    );
    let mut calls = 0;
    let whole = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.1,
        heat(0.01),
        controls(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let continued = advance_deposition_with_equilibrium(
        first.thermal(),
        &selected,
        None,
        first.gas(),
        &catalog,
        0.1,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(whole.gas(), continued.gas());
    assert_eq!(whole.thermal().cells(), continued.thermal().cells());
    assert_eq!(
        whole.audit().condensed_phases,
        continued.audit().condensed_phases
    );
    let mut count = 0;
    assert!(matches!(
        advance_deposition_with_equilibrium(
            &thermal,
            &selected,
            None,
            &gas,
            &catalog,
            0.1,
            heat(0.01),
            controls(),
            || {
                count += 1;
                if count == calls - 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(GasTimelineError::Thermal(DepositionError::Cancelled))
    ));
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.pending_events().len(), 1);
}

#[test]
fn adaptive_plateau_error_and_whole_history_phase_report_include_thermal_conversion() {
    use super::adaptive::{GasTimelineAdaptiveSolve, advance_deposition_adaptive_with_equilibrium};
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = curve(&catalog, &binding, FUSION);
    let phase = fixture::phase();
    let selected = [selection(&curve, &binding, &phase, 1)];
    let h = curve.knots()[1].specific_enthalpy_j_kg + 0.5 * FUSION;
    let thermal = schedule(&curve, vec![birth(0.0, cell(1, h, -20.0))]);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let adaptive = GasTimelineAdaptiveSolve {
        absolute_mass_tolerance_kg: 1e-12,
        absolute_temperature_tolerance_k: 1e-5,
        absolute_specific_enthalpy_tolerance_j_kg: 1e-3,
        relative_tolerance: 0.0,
        minimum_dt_s: 1e-8,
        maximum_trials: 5000,
    };
    let trial = advance_deposition_adaptive_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.2,
        heat(0.2),
        controls(),
        adaptive,
        || Ok(()),
    )
    .unwrap();
    let reference = advance_deposition_with_equilibrium(
        &thermal,
        &selected,
        None,
        &gas,
        &catalog,
        0.2,
        heat(0.0001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let report = &trial.endpoint().audit().condensed_phases[&CellId(1)];
    let reference_phase = &reference.audit().condensed_phases[&CellId(1)];
    fixture::near(report.initial.liquid_mass_kg, 0.005, 1e-14);
    assert!(report.melted_mass_kg < 0.0);
    fixture::near(report.melted_mass_kg, reference_phase.melted_mass_kg, 2e-8);
    assert!(trial.audit().accepted_intervals > 1);
    assert!(trial.audit().maximum_accepted_error_ratio <= 1.0);
    fixture::near(trial.endpoint().audit().energy_residual_j, 0.0, 1e-7);
    fixture::near(
        trial.endpoint().audit().thermal.external_heat_j,
        -4.0,
        1e-11,
    );
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        trial.endpoint().gas().species_mass_kg(&carrier).unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
}
