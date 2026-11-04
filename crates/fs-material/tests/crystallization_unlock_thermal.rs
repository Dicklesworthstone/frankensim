//! G1/G4: original elastic storage enters the original thermal/phase timeline.
//! Inputs are authored analytic fixtures, not experimental material validation.

use fs_blake3::ContentHash;
use fs_material::SmallStrainLaw;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::entropy::CrystallizationEntropy;
use fs_material::phase::crystallization::mechanics::unlocking::CrystallizationUnlock;
use fs_material::phase::crystallization::mechanics::unlocking::thermal::CrystallizationUnlockThermalError;
use fs_material::phase::crystallization::mechanics::{
    CrystallizationLock, CrystallizationMechanicsError, CrystallizationStrainLimits,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_thermal::crystallization::{CrystallizationLaw, CrystallizationSolve};
use fs_thermal::enthalpy::crystallization::{
    CrystallizationControls, CrystallizationImpulseControls, CrystallizationMaterial,
    CrystallizationSimulation, CrystallizationTimelineError,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use std::collections::BTreeMap;

const CELL: CellId = CellId(11);
const NEIGHBOR: CellId = CellId(12);
const MATERIAL: MaterialId = MaterialId(21);
const SOLID: MaterialId = MaterialId(22);
const TM: f64 = 350.0;
const LATENT: f64 = 100_000.0;
const V0_PER_MASS: f64 = 0.001;
const TOTAL: [f64; 6] = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];

fn branch(phase: SolidLiquidPhase, datum: f64, density: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([0x83; 32]),
        phase,
        datum,
        &[280.0, 380.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: 1000.0,
            bulk_density_kg_m3: density,
        }),
        1e-6,
    )
    .unwrap()
}

fn controls() -> CrystallizationControls {
    CrystallizationControls {
        maximum_dt_s: 0.01,
        caloric: CrystallizationSolve::default(),
        cell_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 1e-4,
    }
}

fn impulse() -> CrystallizationImpulseControls {
    CrystallizationImpulseControls {
        caloric: CrystallizationSolve::default(),
        cell_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 1e-7,
    }
}

fn limits() -> CrystallizationStrainLimits {
    CrystallizationStrainLimits {
        maximum_absolute_total_mandel_strain: 0.02,
        maximum_absolute_elastic_mandel_strain: 0.02,
        maximum_absolute_free_isotropic_strain: 0.01,
    }
}

fn elastic() -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(1000.0, 1e10, 4e9, ContentHash([0x84; 32]), 0.02)
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn at_melting<'a>(
    pair: &'a CrystallizationCaloricPair<'a>,
    law: &'a CrystallizationLaw<'a>,
    neighbor: bool,
) -> CrystallizationSimulation<'a> {
    let boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    };
    let cell = |id: CellId, material_id: MaterialId, h: f64, volume: f64, x: f64| DepositedCell {
        id,
        segment_id: 91,
        body_id: id.0,
        material_id,
        mass_kg: 1.0,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [x, 0.0, 0.0],
            end_m: [x + 0.1, 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 0.06,
        },
        boundary,
    };
    let mut materials: Vec<(MaterialId, &'a dyn ThermalLaw)> = vec![(MATERIAL, pair.parent())];
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(cell(
            CELL,
            MATERIAL,
            pair.parent().enthalpy_at_temperature(300.0).unwrap(),
            V0_PER_MASS,
            0.0,
        )),
    }];
    if neighbor {
        materials.push((SOLID, pair.crystal()));
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(
                NEIGHBOR,
                SOLID,
                pair.crystal().enthalpy_at_temperature(300.0).unwrap(),
                V0_PER_MASS,
                0.1,
            )),
        });
        events.push(TimedEvent {
            time_s: 3.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CELL,
                b: NEIGHBOR,
                area_m2: 0.01,
                conductance_w_k: 10.0,
            }),
        });
    }
    // The unlocking impulse must not consume a later physical event.
    events.push(TimedEvent {
        time_s: 10.0,
        event: Event::Boundary(
            CELL,
            Boundary {
                source_w: 25.0,
                ..boundary
            },
        ),
    });
    let thermal = DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 5000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let simulation = CrystallizationSimulation::new(
        thermal,
        &[CrystallizationMaterial {
            material_id: MATERIAL,
            amorphous: pair.parent(),
            law,
            maximum_temperature_slope: pair.maximum_temperature_slope(),
        }],
        || Ok(()),
    )
    .unwrap()
    .with_equilibrium_remelting(MATERIAL, law.equivalent_age_remelting().unwrap())
    .unwrap();
    // The shared solver creates actual crystal history through recalescence.
    let (aged, _) = simulation.advance_to(2.0, controls(), || Ok(())).unwrap();
    let state = aged.states()[&CELL];
    assert!(state.crystalline_mass_fraction() > 0.2);
    assert!(state.temperature_k() < TM);
    let target_h = law
        .specific_enthalpy_j_kg(TM, state.crystalline_mass_fraction())
        .unwrap();
    let heat = aged.cells()[&CELL].mass_kg * (target_h - state.specific_enthalpy_j_kg());
    let (heated, _) = aged
        .apply_heat_impulses(&BTreeMap::from([(CELL, heat)]), impulse(), || Ok(()))
        .unwrap();
    assert_eq!(heated.states()[&CELL].temperature_k(), TM);
    heated
}

fn release(
    thermal: &CrystallizationSimulation<'_>,
    pair: &CrystallizationCaloricPair<'_>,
    elastic: &ElasticAcousticState,
    reference_volume: f64,
    heat_fraction: f64,
) -> CrystallizationUnlock {
    let state = thermal.states()[&CELL];
    CrystallizationLock::new(
        pair,
        state.temperature_k(),
        state.crystalline_mass_fraction(),
        [0.0; 6],
    )
    .unwrap()
    .freeze(
        pair,
        state.temperature_k(),
        state.crystalline_mass_fraction(),
        elastic,
        limits(),
    )
    .unwrap()
    .unlock_at_fixed_strain(TOTAL, reference_volume, heat_fraction)
    .unwrap()
}

#[test]
fn g1_unlock_heats_the_original_latent_plateau_then_cools_and_relocks() {
    let parent = branch(SolidLiquidPhase::Liquid, LATENT, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1000.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([0x85; 32]), 1.0, 0.8, TM)
        .unwrap();
    let thermal = at_melting(&pair, &law, false);
    let initial = thermal.states()[&CELL];
    let original_cell = thermal.cells()[&CELL].clone();
    let elastic = elastic();
    let expected_release = 0.5 * (1e10 + (4.0 / 3.0) * 4e9) * TOTAL[0].powi(2) * V0_PER_MASS;
    let (heated, report) = release(&thermal, &pair, &elastic, V0_PER_MASS, 1.0)
        .thermalize(&pair, &thermal, CELL, impulse(), || Ok(()))
        .unwrap();
    let state = heated.states()[&CELL];
    close(report.removed_elastic_energy_j, expected_release, 1e-10);
    close(report.thermalized_energy_j, expected_release, 1e-10);
    assert_eq!(report.exported_energy_j, 0.0);
    close(
        heated.energy_j() - thermal.energy_j(),
        expected_release,
        1e-8,
    );
    close(report.balance_residual_j, 0.0, 1e-7);
    assert_eq!(state.temperature_k(), TM);
    close(
        initial.crystalline_mass_fraction() - state.crystalline_mass_fraction(),
        expected_release / LATENT,
        2e-13,
    );
    close(
        report.thermal_audit.crystallized_mass_kg,
        -expected_release / LATENT,
        2e-13,
    );
    close(
        report.thermal_audit.thermal.external_heat_j,
        expected_release,
        1e-10,
    );
    assert_eq!(report.thermal_audit.thermal.incoming_mass_kg, 0.0);
    assert_eq!(report.thermal_audit.thermal.incoming_energy_j, 0.0);
    assert_eq!(report.thermal_audit.thermal.substeps, 0);
    assert_eq!(heated.time_s(), thermal.time_s());
    assert_eq!(heated.mass_kg(), thermal.mass_kg());
    assert_eq!(heated.cells()[&CELL].geometry, original_cell.geometry);
    assert_eq!(heated.pending_events(), thermal.pending_events());
    assert_eq!(report.state.stress(&TOTAL, &()), [0.0; 6]);
    assert_eq!(report.state.tangent(&TOTAL, &()), [[0.0; 6]; 6]);
    let entropy = CrystallizationEntropy::new(&pair, TM).unwrap();
    close(
        entropy
            .state_entropy_change_j_kg_k(&law, initial, state, impulse().caloric)
            .unwrap(),
        expected_release / TM,
        1e-10,
    );

    // No elapsed-time or repeated-release heater is hidden in the handoff.
    let (held, hold_audit) = heated.advance_to(2.5, controls(), || Ok(())).unwrap();
    assert_eq!(held.states()[&CELL], state);
    assert_eq!(hold_audit.thermal.external_heat_j, 0.0);
    assert_eq!(held.energy_j(), heated.energy_j());
    let (cooled, _) = held
        .apply_heat_impulses(&BTreeMap::from([(CELL, -20_000.0)]), impulse(), || Ok(()))
        .unwrap();
    assert!(cooled.states()[&CELL].temperature_k() < TM);
    let (grown, _) = cooled.advance_to(2.6, controls(), || Ok(())).unwrap();
    let frozen = grown.states()[&CELL];
    assert!(frozen.crystalline_mass_fraction() > state.crystalline_mass_fraction());
    let relocked = report
        .state
        .relock_at(
            &pair,
            frozen.temperature_k(),
            frozen.crystalline_mass_fraction(),
            TOTAL,
        )
        .unwrap();
    let response = relocked
        .freeze(
            &pair,
            frozen.temperature_k(),
            frozen.crystalline_mass_fraction(),
            &elastic,
            limits(),
        )
        .unwrap();
    assert_eq!(response.stress(&TOTAL, &()), [0.0; 6]);
    assert!(response.tangent(&TOTAL, &())[0][0] > 0.0);
    assert_eq!(thermal.states()[&CELL], initial);
    assert_eq!(thermal.cells()[&CELL], original_cell);
}

#[test]
fn g1_export_is_separate_and_released_heat_reaches_the_original_finite_neighbor() {
    let parent = branch(SolidLiquidPhase::Liquid, LATENT, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1000.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([0x85; 32]), 1.0, 0.8, TM)
        .unwrap();
    let thermal = at_melting(&pair, &law, true);
    let (heated, report) = release(&thermal, &pair, &elastic(), V0_PER_MASS, 0.25)
        .thermalize(&pair, &thermal, CELL, impulse(), || Ok(()))
        .unwrap();
    close(
        report.thermalized_energy_j,
        0.25 * report.removed_elastic_energy_j,
        1e-12,
    );
    close(
        report.exported_energy_j,
        0.75 * report.removed_elastic_energy_j,
        1e-12,
    );
    close(
        heated.energy_j() + report.exported_energy_j,
        thermal.energy_j() + report.removed_elastic_energy_j,
        1e-8,
    );
    assert_eq!(heated.cells()[&NEIGHBOR], thermal.cells()[&NEIGHBOR]);
    let (exchanged, audit) = heated.advance_to(3.2, controls(), || Ok(())).unwrap();
    assert!(exchanged.temperature_k(NEIGHBOR).unwrap() > heated.temperature_k(NEIGHBOR).unwrap());
    assert!(exchanged.cells()[&CELL].specific_h_j_kg < heated.cells()[&CELL].specific_h_j_kg);
    assert_eq!(exchanged.mass_kg(), 2.0);
    assert_eq!(exchanged.links()[&7].a, CELL);
    assert_eq!(exchanged.links()[&7].b, NEIGHBOR);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    close(exchanged.energy_j(), heated.energy_j(), 1e-7);
    assert_eq!(exchanged.pending_events().len(), 1);
    assert_eq!(exchanged.pending_events()[0].time_s, 10.0);
}

#[test]
fn g4_binding_volume_work_and_late_cancellation_preserve_the_accepted_timeline() {
    let parent = branch(SolidLiquidPhase::Liquid, LATENT, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1000.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([0x85; 32]), 1.0, 0.8, TM)
        .unwrap();
    let thermal = at_melting(&pair, &law, false);
    let initial = thermal.states().clone();
    let cells = thermal.cells().clone();
    let events = thermal.pending_events().to_vec();
    let elastic = elastic();
    let wrong_crystal = branch(SolidLiquidPhase::Solid, 1.0, 1000.0);
    let wrong_pair =
        CrystallizationCaloricPair::new(&parent, &wrong_crystal, [280.0, 380.0]).unwrap();
    assert!(matches!(
        release(&thermal, &pair, &elastic, V0_PER_MASS, 1.0).thermalize(
            &wrong_pair,
            &thermal,
            CELL,
            impulse(),
            || Ok(())
        ),
        Err(CrystallizationUnlockThermalError::Mechanics(
            CrystallizationMechanicsError::PairMismatch
        ))
    ));
    assert!(matches!(
        release(&thermal, &pair, &elastic, 1.1 * V0_PER_MASS, 1.0).thermalize(
            &pair,
            &thermal,
            CELL,
            impulse(),
            || Ok(())
        ),
        Err(CrystallizationUnlockThermalError::Mechanics(
            CrystallizationMechanicsError::Invalid(
                "mechanical and thermal fixed reference V0/m0 differ"
            )
        ))
    ));
    // Locate the final successful checkpoint, then cancel at that same point
    // in a fresh proposal from the unchanged accepted mechanical/thermal state.
    let mut checkpoints = 0;
    let _candidate = release(&thermal, &pair, &elastic, V0_PER_MASS, 1.0)
        .thermalize(&pair, &thermal, CELL, impulse(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(checkpoints > 1);
    let mut visited = 0;
    assert!(matches!(
        release(&thermal, &pair, &elastic, V0_PER_MASS, 1.0).thermalize(
            &pair,
            &thermal,
            CELL,
            impulse(),
            || {
                visited += 1;
                if visited == checkpoints {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(CrystallizationUnlockThermalError::Thermal(
            CrystallizationTimelineError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(visited, checkpoints);
    assert_eq!(thermal.states(), &initial);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.pending_events(), events.as_slice());
    assert_eq!(thermal.time_s(), 2.0);

    // Equal source identities do not authorize unmodeled pressure/volume work.
    // This matched pair changes its free volume on latent remelting, which the
    // fixed-geometry thermalization bridge must explicitly refuse.
    let dense_crystal = branch(SolidLiquidPhase::Solid, 0.0, 1030.0);
    let unequal = CrystallizationCaloricPair::new(&parent, &dense_crystal, [280.0, 380.0]).unwrap();
    let unequal_law = unequal
        .law(&rate, ContentHash([0x85; 32]), 1.0, 0.8, TM)
        .unwrap();
    let unequal_timeline = at_melting(&unequal, &unequal_law, false);
    let previous = unequal_timeline.states()[&CELL];
    assert!(matches!(
        release(&unequal_timeline, &unequal, &elastic, V0_PER_MASS, 1.0).thermalize(
            &unequal,
            &unequal_timeline,
            CELL,
            impulse(),
            || Ok(())
        ),
        Err(CrystallizationUnlockThermalError::Mechanics(
            CrystallizationMechanicsError::Invalid(
                "density-changing unlocking requires a pressure/volume-work owner"
            )
        ))
    ));
    assert_eq!(unequal_timeline.states()[&CELL], previous);
}
