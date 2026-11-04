//! G1/G4: original material calorics, kinetic history and a constant-pressure load.
//! Authored analytic fixtures establish model behavior, not experimental validity.

use fs_blake3::ContentHash;
use fs_material::SmallStrainLaw;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::entropy::CrystallizationEntropy;
use fs_material::phase::crystallization::mechanics::unlocking::CrystallizationUnlock;
use fs_material::phase::crystallization::mechanics::unlocking::isobaric::CrystallizationUnlockIsobaricError;
use fs_material::phase::crystallization::mechanics::{
    CrystallizationLock, CrystallizationMechanicsError, CrystallizationStrainLimits,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, PhaseIntervalDensityLaw,
    SolidLiquidPhase,
};
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::pressure::{
    IsobaricCrystallizationControls, IsobaricCrystallizationStep,
};
use fs_thermal::crystallization::{
    CrystallizationAdaptivity, CrystallizationCaloric, CrystallizationError, CrystallizationSolve,
    CrystallizationState,
};
use fs_thermal::enthalpy::{ThermalLaw, ThermalStorageBasis};

const TM: f64 = 350.0;
const CP: f64 = 1000.0;
const LATENT: f64 = 100_000.0;
const PRESSURE: f64 = 2_000_000.0;
const TOTAL: [f64; 6] = [0.01, -0.01, 0.0, 0.0, 0.0, 0.0];

fn branch(
    phase: SolidLiquidPhase,
    datum: f64,
    density: [f64; 2],
    affine_density: bool,
) -> EquilibriumEnthalpyPhaseCurve {
    let values = [(280.0, density[0]), (420.0, density[1])];
    let card = ContentHash([0x94; 32]);
    if affine_density {
        EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            card,
            phase,
            datum,
            &values.map(|(temperature_k, bulk_density_kg_m3)| HeatCapacityKnot {
                temperature_k,
                specific_heat_capacity_j_kg_k: CP,
                bulk_density_kg_m3,
            }),
            1e-6,
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card,
            phase,
            values
                .map(|(temperature_k, bulk_density_kg_m3)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: datum + CP * (temperature_k - 280.0),
                    temperature_k,
                    liquid_mass_fraction: if phase == SolidLiquidPhase::Liquid {
                        1.0
                    } else {
                        0.0
                    },
                    bulk_density_kg_m3,
                })
                .to_vec(),
        )
        .unwrap()
    }
}

fn controls() -> IsobaricCrystallizationControls {
    IsobaricCrystallizationControls {
        caloric: CrystallizationSolve::default(),
        energy_tolerance_j: 1e-6,
        volume_tolerance_m3: 1e-15,
    }
}

fn adaptivity() -> CrystallizationAdaptivity {
    CrystallizationAdaptivity {
        maximum_dt_s: 0.01,
        temperature_tolerance_k: 1e-5,
        crystalline_mass_fraction_tolerance: 1e-8,
        ..CrystallizationAdaptivity::default()
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn energy_closes(step: &IsobaricCrystallizationStep) {
    close(
        step.internal_energy_change_j + step.boundary_work_j,
        step.enthalpy_change_j,
        1e-8,
    );
    close(step.enthalpy_change_j, step.supplied_heat_j, 1e-6);
    assert!(step.absolute_energy_residual_j <= 1e-6);
    assert!((step.volume_change_m3 - step.source_volume_change_m3).abs() <= 1e-15);
}

#[test]
fn g1_original_kinetics_remelting_and_entropy_close_a_pressure_load_cycle() {
    let parent = branch(SolidLiquidPhase::Liquid, LATENT, [800.0; 2], true);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, [1000.0; 2], true);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 420.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, TM).unwrap();
    let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([0x95; 32]), 1.0, 0.8, TM)
        .unwrap();
    let initial = law.state_at_temperature(300.0, 0.0).unwrap();
    let mass = 2.0;
    let parcel = IsobaricCrystallizationParcel::new(
        &law,
        initial,
        mass,
        mass / 800.0,
        PRESSURE,
        controls(),
        || false,
    )
    .unwrap()
    .with_equilibrium_remelting(law.equivalent_age_remelting().unwrap())
    .unwrap();
    let (grown, growth) = parcel
        .advance_heat(0.0, 0.2, controls(), adaptivity(), || false)
        .unwrap();
    let x = grown.state().crystalline_mass_fraction();
    assert!(x > 0.0);
    close(grown.state().temperature_k(), 300.0 + LATENT * x / CP, 1e-8);
    assert!(growth.adaptive.unwrap().accepted_substeps > 0);
    close(growth.enthalpy_change_j, 0.0, 1e-10);
    let contraction = mass * x * (1.0 / 1000.0 - 1.0 / 800.0);
    close(growth.source_volume_change_m3, contraction, 1e-18);
    close(growth.boundary_work_j, PRESSURE * contraction, 1e-9);
    assert!(growth.boundary_work_j < 0.0);
    assert!(growth.internal_energy_change_j > 0.0);
    // With no heat, contraction transfers work IN; constant H is not constant U.
    close(
        growth.internal_energy_change_j,
        -growth.boundary_work_j,
        1e-9,
    );
    assert!(
        entropy
            .state_entropy_change_j_kg_k(&law, initial, grown.state(), controls().caloric)
            .unwrap()
            > 0.0
    );
    energy_closes(&growth);

    let sensible_heat = mass * CP * (TM - grown.state().temperature_k());
    let (at_melting, sensible) = grown
        .apply_heat(sensible_heat, controls(), || false)
        .unwrap();
    close(at_melting.state().temperature_k(), TM, 1e-9);
    close(at_melting.state().crystalline_mass_fraction(), x, 1e-12);
    let heat = 0.5 * mass * LATENT * at_melting.state().crystalline_mass_fraction();
    let (partial, remelt) = at_melting.apply_heat(heat, controls(), || false).unwrap();
    close(partial.state().temperature_k(), TM, 1e-9);
    close(partial.state().crystalline_mass_fraction(), 0.5 * x, 1e-12);
    let dx = partial.state().crystalline_mass_fraction()
        - at_melting.state().crystalline_mass_fraction();
    let expansion = mass * dx * (1.0 / 1000.0 - 1.0 / 800.0);
    close(remelt.source_volume_change_m3, expansion, 1e-18);
    close(remelt.boundary_work_j, PRESSURE * expansion, 1e-9);
    assert!(remelt.boundary_work_j > 0.0);
    let ds = mass
        * entropy
            .state_entropy_change_j_kg_k(
                &law,
                at_melting.state(),
                partial.state(),
                controls().caloric,
            )
            .unwrap();
    close(ds, heat / TM, 1e-9);
    energy_closes(&remelt);

    let (liquid, complete) = partial
        .apply_heat(
            mass * LATENT * partial.state().crystalline_mass_fraction(),
            controls(),
            || false,
        )
        .unwrap();
    close(liquid.state().crystalline_mass_fraction(), 0.0, 1e-14);
    let (returned, cooling) = liquid
        .apply_heat(mass * CP * (300.0 - TM), controls(), || false)
        .unwrap();
    close(returned.state().temperature_k(), 300.0, 1e-9);
    close(returned.volume_m3(), parcel.volume_m3(), 1e-17);
    assert_eq!(returned.time_s(), grown.time_s());
    assert_eq!(returned.mass_kg(), mass);
    assert_eq!(returned.pressure_pa(), PRESSURE);
    assert_eq!(returned.reference_volume_m3(), parcel.reference_volume_m3());
    let cycle = [growth, sensible, remelt, complete, cooling];
    for step in &cycle {
        energy_closes(step);
    }
    close(cycle.iter().map(|s| s.supplied_heat_j).sum(), 0.0, 1e-6);
    close(cycle.iter().map(|s| s.boundary_work_j).sum(), 0.0, 1e-8);
    close(
        cycle.iter().map(|s| s.internal_energy_change_j).sum(),
        0.0,
        1e-6,
    );
    close(
        entropy
            .state_entropy_change_j_kg_k(&law, initial, returned.state(), controls().caloric)
            .unwrap(),
        0.0,
        1e-8,
    );
}

#[test]
fn g1_original_density_interpolation_and_common_enthalpy_datum_are_preserved() {
    let mut work_by_interpolation = Vec::new();
    for affine_density in [false, true] {
        let mut unshifted = None;
        for datum in [0.0, 1e9] {
            let parent = branch(
                SolidLiquidPhase::Liquid,
                datum + LATENT,
                [1000.0, 800.0],
                affine_density,
            );
            let crystal = branch(
                SolidLiquidPhase::Solid,
                datum,
                [1200.0, 1000.0],
                affine_density,
            );
            assert_eq!(
                parent.density_law_on_interval(0).unwrap(),
                if affine_density {
                    PhaseIntervalDensityLaw::AffineDensity
                } else {
                    PhaseIntervalDensityLaw::AffineSpecificVolume
                }
            );
            let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 420.0]).unwrap();
            let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
            let law = pair
                .law(&rate, ContentHash([0x96; 32]), 2.0, 0.8, TM)
                .unwrap();
            let initial = law.state_at_temperature(300.0, 0.3).unwrap();
            let x = initial.crystalline_mass_fraction();
            let mut budget = controls();
            budget.caloric.specific_enthalpy_tolerance_j_kg = 1e-6;
            budget.energy_tolerance_j = 1e-5;
            let mass = 1.5;
            let parcel = IsobaricCrystallizationParcel::new(
                &law,
                initial,
                mass,
                0.0015,
                PRESSURE,
                budget,
                || false,
            )
            .unwrap();
            let (warm, step) = parcel
                .apply_heat(mass * CP * 40.0, budget, || false)
                .unwrap();
            assert_eq!(warm.state().crystalline_mass_fraction(), x);
            close(warm.state().temperature_k(), 340.0, 2e-9);
            let v = |t: f64, rho: [f64; 2]| {
                let f = (t - 280.0) / 140.0;
                if affine_density {
                    (rho[0] + f * (rho[1] - rho[0])).recip()
                } else {
                    (1.0 - f) / rho[0] + f / rho[1]
                }
            };
            let mixture = |t| (1.0 - x) * v(t, [1000.0, 800.0]) + x * v(t, [1200.0, 1000.0]);
            close(parcel.volume_m3(), mass * mixture(300.0), 1e-16);
            close(warm.volume_m3(), mass * mixture(340.0), 1e-14);
            close(
                step.boundary_work_j,
                PRESSURE * mass * (mixture(340.0) - mixture(300.0)),
                2e-8,
            );
            assert_eq!(
                pair.parent().density_kg_m3(320.0).unwrap(),
                Some(pair.parent().density_at_temperature(320.0).unwrap())
            );
            assert!(pair.parent().density_kg_m3(270.0).is_err());
            assert_eq!(pair.parent().storage_basis(), ThermalStorageBasis::Enthalpy);
            if let Some(previous) = unshifted {
                close(step.boundary_work_j, previous, 2e-8);
            } else {
                unshifted = Some(step.boundary_work_j);
                work_by_interpolation.push(step.boundary_work_j);
            }
        }
    }
    // These are distinct admitted source interpolations, not interchangeable fits.
    assert!((work_by_interpolation[0] - work_by_interpolation[1]).abs() > 1.0);
}

fn release(
    pair: &CrystallizationCaloricPair<'_>,
    state: CrystallizationState,
    reference_specific_volume: f64,
) -> CrystallizationUnlock {
    release_at(pair, state, reference_specific_volume, TOTAL)
}

fn release_at(
    pair: &CrystallizationCaloricPair<'_>,
    state: CrystallizationState,
    reference_specific_volume: f64,
    total: [f64; 6],
) -> CrystallizationUnlock {
    let density = pair
        .volume_state(state.temperature_k(), state.crystalline_mass_fraction())
        .unwrap()
        .bulk_density_kg_m3;
    let elastic = ElasticAcousticState::from_bulk_and_shear(
        density,
        1e10,
        4e9,
        ContentHash([0x97; 32]),
        0.02,
    )
    .unwrap();
    let lock = CrystallizationLock::new(
        pair,
        state.temperature_k(),
        state.crystalline_mass_fraction(),
        [0.0; 6],
    )
    .unwrap();
    lock.freeze(
        pair,
        state.temperature_k(),
        state.crystalline_mass_fraction(),
        &elastic,
        CrystallizationStrainLimits {
            maximum_absolute_total_mandel_strain: 0.02,
            maximum_absolute_elastic_mandel_strain: 0.02,
            maximum_absolute_free_isotropic_strain: 0.01,
        },
    )
    .unwrap()
    .unlock_at_fixed_strain(total, reference_specific_volume, 0.75)
    .unwrap()
}

#[test]
fn g4_original_elastic_release_owns_pressure_work_and_atomic_refusals() {
    let parent = branch(SolidLiquidPhase::Liquid, LATENT, [800.0; 2], true);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, [1000.0; 2], true);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 420.0]).unwrap();
    let entropy = CrystallizationEntropy::new(&pair, TM).unwrap();
    let rate = |t: f64| Ok(0.01 * (TM - t).max(0.0));
    let law = pair
        .law(&rate, ContentHash([0x98; 32]), 2.0, 0.8, TM)
        .unwrap();
    let initial = law.state_at_temperature(TM, 0.8).unwrap();
    let parcel =
        IsobaricCrystallizationParcel::new(&law, initial, 2.0, 0.002, PRESSURE, controls(), || {
            false
        })
        .unwrap()
        .with_equilibrium_remelting(law.equivalent_age_remelting().unwrap())
        .unwrap();
    let mut checkpoints = 0;
    let (melted, report) = release(&pair, initial, 0.001)
        .thermalize_isobaric(&pair, &parcel, controls(), || {
            checkpoints += 1;
            false
        })
        .unwrap();
    close(report.removed_elastic_energy_j, 1600.0, 1e-9);
    close(report.thermalized_energy_j, 1200.0, 1e-9);
    close(report.exported_energy_j, 400.0, 1e-9);
    close(melted.state().temperature_k(), TM, 1e-9);
    close(
        melted.state().crystalline_mass_fraction(),
        initial.crystalline_mass_fraction() - 0.006,
        1e-12,
    );
    close(report.thermal.boundary_work_j, 6.0, 1e-9);
    close(report.thermal.internal_energy_change_j, 1194.0, 1e-8);
    close(
        report.thermal.internal_energy_change_j
            + report.thermal.boundary_work_j
            + report.exported_energy_j
            - report.removed_elastic_energy_j,
        0.0,
        1e-8,
    );
    assert!(report.balance_residual_j.abs() <= controls().energy_tolerance_j);
    assert_eq!(report.state.stress(&TOTAL, &()), [0.0; 6]);
    assert!(melted.volume_m3() > parcel.volume_m3());
    assert_eq!(melted.reference_volume_m3(), parcel.reference_volume_m3());
    assert_eq!(melted.time_s(), parcel.time_s());
    close(
        2.0 * entropy
            .state_entropy_change_j_kg_k(&law, initial, melted.state(), controls().caloric)
            .unwrap(),
        report.thermalized_energy_j / TM,
        1e-9,
    );
    energy_closes(&report.thermal);

    let different_parent = branch(SolidLiquidPhase::Liquid, LATENT, [801.0; 2], true);
    let different_pair =
        CrystallizationCaloricPair::new(&different_parent, &crystal, [280.0, 420.0]).unwrap();
    assert_ne!(different_pair.identity(), pair.identity());
    assert!(matches!(
        release(&pair, initial, 0.001).thermalize_isobaric(
            &different_pair,
            &parcel,
            controls(),
            || false
        ),
        Err(CrystallizationUnlockIsobaricError::Mechanics(
            CrystallizationMechanicsError::PairMismatch
        ))
    ));
    assert!(matches!(
        release(&pair, initial, 0.0011).thermalize_isobaric(&pair, &parcel, controls(), || false),
        Err(CrystallizationUnlockIsobaricError::Mechanics(
            CrystallizationMechanicsError::Invalid(_)
        ))
    ));
    // Hydrostatic relaxation needs a compressible elastic-volume owner; it
    // cannot be spent as heat while retaining this same pressure-load model.
    assert!(matches!(
        release_at(&pair, initial, 0.001, [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]).thermalize_isobaric(
            &pair,
            &parcel,
            controls(),
            || false
        ),
        Err(CrystallizationUnlockIsobaricError::Mechanics(
            CrystallizationMechanicsError::Invalid(
                "isobaric unlocking requires an owner for elastic volume relaxation"
            )
        ))
    ));
    assert!(checkpoints > 1);
    let mut cancelled_at = 0;
    let cancelled =
        release(&pair, initial, 0.001).thermalize_isobaric(&pair, &parcel, controls(), || {
            cancelled_at += 1;
            cancelled_at == checkpoints
        });
    assert!(matches!(
        cancelled,
        Err(CrystallizationUnlockIsobaricError::Thermal(
            CrystallizationError::Cancelled
        ))
    ));
    assert_eq!(parcel.state(), initial);
    assert_eq!(parcel.time_s(), 0.0);
    close(
        parcel.volume_m3(),
        2.0 * pair
            .volume_state(initial.temperature_k(), initial.crystalline_mass_fraction())
            .unwrap()
            .specific_volume_m3_kg,
        1e-18,
    );
}
