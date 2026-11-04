//! G1/G4: finite thermal contact advances original kinetic material inventories.
//! These are authored analytic material/contact models, not experimental claims.
//! Total contact conductance is an explicitly admitted input. No conductivity
//! mixture rule, contact area or geometry is inferred from kinetic crystallinity.

use fs_blake3::ContentHash;
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::entropy::CrystallizationEntropy;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, PhaseIntervalDensityLaw,
    SolidLiquidPhase,
};
use fs_thermal::crystallization::pressure::exchange::{
    IsobaricExchangeControls, IsobaricExchangeSideControls, IsobaricExchangeStep, exchange_heat,
};
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::{CrystallizationAdaptivity, CrystallizationError};

const RANGE: [f64; 2] = [280.0, 550.0];

fn branch(
    material: u8,
    phase: SolidLiquidPhase,
    cp: f64,
    datum: f64,
    densities: [f64; 2],
    affine_density: bool,
) -> EquilibriumEnthalpyPhaseCurve {
    let values = [(RANGE[0], densities[0]), (RANGE[1], densities[1])];
    if affine_density {
        EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([material; 32]),
            phase,
            datum,
            &values.map(|(temperature_k, bulk_density_kg_m3)| HeatCapacityKnot {
                temperature_k,
                specific_heat_capacity_j_kg_k: cp,
                bulk_density_kg_m3,
            }),
            1e-6,
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([material; 32]),
            phase,
            values
                .map(|(temperature_k, bulk_density_kg_m3)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: datum + cp * (temperature_k - RANGE[0]),
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

fn side_controls() -> IsobaricExchangeSideControls {
    let mut controls = IsobaricExchangeSideControls::default();
    controls.caloric.energy_tolerance_j = 1e-6;
    controls.caloric.volume_tolerance_m3 = 1e-15;
    controls.adaptivity = CrystallizationAdaptivity {
        maximum_dt_s: 0.05,
        temperature_tolerance_k: 1e-5,
        crystalline_mass_fraction_tolerance: 1e-8,
        ..CrystallizationAdaptivity::default()
    };
    controls
}

fn exchange_controls(upper_heat_j: f64) -> IsobaricExchangeControls {
    IsobaricExchangeControls {
        heat_bracket_j: [0.0, upper_heat_j],
        maximum_evaluations: 64,
        residual_tolerance_j: 1e-6,
        energy_tolerance_j: 1e-5,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn closed_contact(report: &IsobaricExchangeStep) {
    assert!(report.absolute_energy_residual_j <= 1e-5);
    assert!(report.residual_j.abs() <= 1e-5);
    close(report.combined_enthalpy_change_j, 0.0, 1e-6);
    close(
        report.combined_internal_energy_change_j + report.combined_boundary_work_j,
        0.0,
        1e-6,
    );
    close(
        report.heat_left_to_right_j,
        report.conductance_w_k
            * report.duration_s
            * (report.left.state.temperature_k() - report.right.state.temperature_k()),
        1e-5,
    );
    for step in [&report.left, &report.right] {
        close(
            step.internal_energy_change_j + step.boundary_work_j,
            step.supplied_heat_j,
            1e-6,
        );
        assert!(step.adaptive.unwrap().accepted_substeps > 0);
    }
}

#[test]
fn g1_finite_amorphous_solid_heats_a_distinct_melting_material_with_original_entropy() {
    // The first material's parent is an amorphous SOLID; its kinetic fraction
    // cannot be substituted for a liquid fraction in a conductivity source.
    let lp = branch(
        0xa1,
        SolidLiquidPhase::Solid,
        2000.0,
        100_000.0,
        [1000.0, 900.0],
        true,
    );
    let lc = branch(
        0xa1,
        SolidLiquidPhase::Solid,
        2000.0,
        0.0,
        [1200.0, 1100.0],
        true,
    );
    let rp = branch(
        0xa2,
        SolidLiquidPhase::Liquid,
        1000.0,
        80_000.0,
        [800.0, 740.0],
        false,
    );
    let rc = branch(
        0xa2,
        SolidLiquidPhase::Solid,
        1000.0,
        0.0,
        [1000.0, 940.0],
        false,
    );
    assert_eq!(
        lp.density_law_on_interval(0).unwrap(),
        PhaseIntervalDensityLaw::AffineDensity
    );
    assert_eq!(
        rp.density_law_on_interval(0).unwrap(),
        PhaseIntervalDensityLaw::AffineSpecificVolume
    );
    let left_pair = CrystallizationCaloricPair::new(&lp, &lc, RANGE).unwrap();
    let right_pair = CrystallizationCaloricPair::new(&rp, &rc, RANGE).unwrap();
    let no_growth = |_: f64| Ok(0.0);
    let growth = |t: f64| Ok(0.01 * (350.0 - t).max(0.0));
    let left_law = left_pair
        .law(&no_growth, ContentHash([0xa3; 32]), 1.0, 0.8, 500.0)
        .unwrap();
    let right_law = right_pair
        .law(&growth, ContentHash([0xa4; 32]), 1.0, 0.8, 350.0)
        .unwrap();
    let left_entropy = CrystallizationEntropy::new(&left_pair, 500.0).unwrap();
    let right_entropy = CrystallizationEntropy::new(&right_pair, 350.0).unwrap();
    let controls = side_controls();
    let left_initial = left_law.state_at_temperature(400.0, 0.0).unwrap();
    let right_initial = right_law.state_at_temperature(350.0, 1.0).unwrap();
    let left = IsobaricCrystallizationParcel::new(
        &left_law,
        left_initial,
        1.0,
        0.001,
        1e6,
        controls.caloric,
        || false,
    )
    .unwrap();
    let right = IsobaricCrystallizationParcel::new(
        &right_law,
        right_initial,
        0.5,
        0.0005,
        2e6,
        controls.caloric,
        || false,
    )
    .unwrap()
    .with_equilibrium_remelting(right_law.equivalent_age_remelting().unwrap())
    .unwrap();
    let (cooled, melted, report) = exchange_heat(
        &left,
        &right,
        100.0,
        1.0,
        controls,
        controls,
        exchange_controls(8000.0),
        || false,
    )
    .unwrap();
    // The receiver remains on its original plateau; the donor is a FINITE
    // 2000 J/K body. This is the exact backward-Euler contact solution.
    let expected_heat = 100.0 * (400.0 - 350.0) / (1.0 + 100.0 / 2000.0);
    close(report.heat_left_to_right_j, expected_heat, 1e-5);
    close(
        cooled.state().temperature_k(),
        400.0 - expected_heat / 2000.0,
        1e-8,
    );
    close(melted.state().temperature_k(), 350.0, 1e-9);
    close(
        melted.state().crystalline_mass_fraction(),
        right_initial.crystalline_mass_fraction() - expected_heat / (0.5 * 80_000.0),
        1e-9,
    );
    assert!(melted.state().crystalline_mass_fraction() > 0.0);
    assert_eq!(cooled.state().crystalline_mass_fraction(), 0.0);
    let left_rho = |t| 1000.0 - 100.0 * (t - 280.0) / 270.0;
    let left_dv = 1.0 / left_rho(cooled.state().temperature_k()) - 1.0 / left_rho(400.0);
    let original_v = |density: [f64; 2]| {
        let f = (350.0 - 280.0) / 270.0;
        (1.0 - f) / density[0] + f / density[1]
    };
    let right_dv = 0.5
        * (melted.state().crystalline_mass_fraction() - right_initial.crystalline_mass_fraction())
        * (original_v([1000.0, 940.0]) - original_v([800.0, 740.0]));
    close(report.left.boundary_work_j, 1e6 * left_dv, 1e-8);
    close(report.right.boundary_work_j, 2e6 * right_dv, 1e-8);
    assert!(report.left.boundary_work_j < 0.0);
    assert!(report.right.boundary_work_j > 0.0);
    let ds_left = left_entropy
        .state_entropy_change_j_kg_k(
            &left_law,
            left_initial,
            cooled.state(),
            controls.caloric.caloric,
        )
        .unwrap();
    let ds_right = 0.5
        * right_entropy
            .state_entropy_change_j_kg_k(
                &right_law,
                right_initial,
                melted.state(),
                controls.caloric.caloric,
            )
            .unwrap();
    let exact_ds = 2000.0 * (cooled.state().temperature_k() / 400.0).ln()
        + report.heat_left_to_right_j / 350.0;
    close(ds_left + ds_right, exact_ds, 1e-8);
    assert!(ds_left + ds_right > 0.0);
    assert_eq!(cooled.time_s(), 1.0);
    assert_eq!(melted.time_s(), 1.0);
    assert_eq!(cooled.reference_volume_m3(), left.reference_volume_m3());
    assert_eq!(melted.reference_volume_m3(), right.reference_volume_m3());
    closed_contact(&report);
}

#[test]
fn g4_contact_cooling_advances_original_crystallization_and_cancellation_is_atomic() {
    let lp = branch(
        0xb1,
        SolidLiquidPhase::Liquid,
        1000.0,
        20_000.0,
        [800.0; 2],
        false,
    );
    let lc = branch(
        0xb1,
        SolidLiquidPhase::Solid,
        1000.0,
        0.0,
        [1000.0; 2],
        false,
    );
    let rp = branch(
        0xb2,
        SolidLiquidPhase::Solid,
        2000.0,
        60_000.0,
        [1100.0, 900.0],
        true,
    );
    let rc = branch(
        0xb2,
        SolidLiquidPhase::Solid,
        2000.0,
        0.0,
        [1300.0, 1100.0],
        true,
    );
    let left_pair = CrystallizationCaloricPair::new(&lp, &lc, RANGE).unwrap();
    let right_pair = CrystallizationCaloricPair::new(&rp, &rc, RANGE).unwrap();
    // A continuous authored rate has a plateau below 340 K and vanishes at
    // Tm. Contact cooling stays on that plateau; adiabatic recalescence does not.
    let growth = |t: f64| {
        Ok(if t <= 340.0 {
            0.1
        } else {
            0.01 * (350.0 - t).max(0.0)
        })
    };
    let no_growth = |_: f64| Ok(0.0);
    let left_law = left_pair
        .law(&growth, ContentHash([0xb3; 32]), 1.0, 0.8, 350.0)
        .unwrap();
    let right_law = right_pair
        .law(&no_growth, ContentHash([0xb4; 32]), 1.0, 0.8, 500.0)
        .unwrap();
    let left_entropy = CrystallizationEntropy::new(&left_pair, 350.0).unwrap();
    let right_entropy = CrystallizationEntropy::new(&right_pair, 500.0).unwrap();
    let controls = side_controls();
    let left_initial = left_law.state_at_temperature(340.0, 0.0).unwrap();
    let right_initial = right_law.state_at_temperature(300.0, 0.0).unwrap();
    let left = IsobaricCrystallizationParcel::new(
        &left_law,
        left_initial,
        1.0,
        0.00125,
        1e6,
        controls.caloric,
        || false,
    )
    .unwrap();
    let right = IsobaricCrystallizationParcel::new(
        &right_law,
        right_initial,
        2.0,
        0.002,
        2e6,
        controls.caloric,
        || false,
    )
    .unwrap();
    let (isolated, _) = left
        .advance_heat(0.0, 2.0, controls.caloric, controls.adaptivity, || false)
        .unwrap();
    let mut checkpoints = 0;
    let (cooled, warmed, report) = exchange_heat(
        &left,
        &right,
        200.0,
        2.0,
        controls,
        controls,
        exchange_controls(18_000.0),
        || {
            checkpoints += 1;
            false
        },
    )
    .unwrap();
    let x = cooled.state().crystalline_mass_fraction();
    close(x, 0.8 * (1.0 - (-0.2_f64).exp()), 1e-10);
    close(
        report.heat_left_to_right_j,
        (16_000.0 + 8000.0 * x) / 1.5,
        1e-5,
    );
    assert!(x > isolated.state().crystalline_mass_fraction());
    assert!(cooled.state().temperature_k() < left_initial.temperature_k());
    assert!(warmed.state().temperature_k() > right_initial.temperature_k());
    close(
        cooled.state().temperature_k(),
        340.0 - report.heat_left_to_right_j / 1000.0 + 20.0 * x,
        1e-8,
    );
    close(
        warmed.state().temperature_k(),
        300.0 + report.heat_left_to_right_j / 4000.0,
        1e-8,
    );
    close(
        report.left.boundary_work_j,
        1e6 * x * (1.0 / 1000.0 - 1.0 / 800.0),
        1e-8,
    );
    let rho = |t| 1100.0 - 200.0 * (t - 280.0) / 270.0;
    close(
        report.right.boundary_work_j,
        2e6 * (2.0 / rho(warmed.state().temperature_k()) - 2.0 / rho(300.0)),
        1e-8,
    );
    let ds = left_entropy
        .state_entropy_change_j_kg_k(
            &left_law,
            left_initial,
            cooled.state(),
            controls.caloric.caloric,
        )
        .unwrap()
        + 2.0
            * right_entropy
                .state_entropy_change_j_kg_k(
                    &right_law,
                    right_initial,
                    warmed.state(),
                    controls.caloric.caloric,
                )
                .unwrap();
    let exact_ds = 1000.0 * (cooled.state().temperature_k() / 340.0).ln() - 20_000.0 * x / 350.0
        + 4000.0 * (warmed.state().temperature_k() / 300.0).ln();
    close(ds, exact_ds, 1e-8);
    assert!(ds > 0.0);
    assert_eq!(cooled.state().identity(), left_law.identity());
    assert_eq!(warmed.state().identity(), right_law.identity());
    assert_eq!(cooled.time_s(), warmed.time_s());
    closed_contact(&report);

    // Cancellation at the last checkpoint discards BOTH computed candidates.
    assert!(checkpoints > 1);
    let mut polls = 0;
    let cancelled = exchange_heat(
        &left,
        &right,
        200.0,
        2.0,
        controls,
        controls,
        exchange_controls(18_000.0),
        || {
            polls += 1;
            polls == checkpoints
        },
    );
    assert!(matches!(cancelled, Err(CrystallizationError::Cancelled)));
    assert_eq!(left.state(), left_initial);
    assert_eq!(right.state(), right_initial);
    assert_eq!(left.time_s(), 0.0);
    assert_eq!(right.time_s(), 0.0);
}
