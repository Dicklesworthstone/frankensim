use super::*;
use crate::phase::{HeatCapacityKnot, SolidLiquidPhase};

fn knot(h: f64, t: f64, f: f64, rho: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    }
}

fn matrix(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([1; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(datum, 300.0, 0.0, 8000.0),
            knot(datum + 200_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap()
}

fn pcm(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([2; 32]),
        vec![
            knot(datum, 300.0, 0.0, 1000.0),
            knot(datum + 200_000.0, 400.0, 0.0, 1000.0),
            knot(datum + 400_000.0, 400.0, 1.0, 800.0),
            knot(datum + 650_000.0, 500.0, 1.0, 800.0),
        ],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn latent_storage_and_density_use_fixed_constituent_masses() {
    let matrix = matrix(0.0);
    let pcm = pcm(0.0);
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    assert!(std::ptr::eq(law.matrix_curve(), &matrix));
    assert!(std::ptr::eq(law.pcm_curve(), &pcm));
    assert_eq!(law.enthalpy_range_j_kg(), [0.0, 312_500.0]);
    assert_eq!(
        law.knots()
            .iter()
            .map(|k| k.specific_enthalpy_j_kg)
            .collect::<Vec<_>>(),
        vec![0.0, 125_000.0, 175_000.0, 312_500.0]
    );
    close(law.maximum_temperature_slope(), 1.0 / 1250.0, 1e-18);
    close(
        law.temperature_at_specific_enthalpy(62_500.0).unwrap(),
        350.0,
        1e-12,
    );
    let half = law.values_at_specific_enthalpy(150_000.0).unwrap();
    assert_eq!(half.temperature_k, 400.0);
    assert_eq!(half.matrix.liquid_mass_fraction, 0.0);
    assert_eq!(half.pcm.liquid_mass_fraction, 0.5);
    assert_eq!(half.liquid_mass_fraction, 0.125);
    let pcm_specific_volume = 0.5 / 1000.0 + 0.5 / 800.0;
    let volume = 0.75 / 8000.0 + 0.25 * pcm_specific_volume;
    close(half.specific_volume_m3_kg, volume, 1e-19);
    close(half.bulk_density_kg_m3, 1.0 / volume, 1e-10);
    close(
        half.pcm_volume_fraction,
        0.25 * pcm_specific_volume / volume,
        1e-15,
    );
    assert!((half.pcm_volume_fraction - half.pcm_mass_fraction).abs() > 0.4);
    assert_eq!(
        law.values_at_specific_enthalpy(175_000.0)
            .unwrap()
            .liquid_mass_fraction,
        0.25
    );
    close(
        law.temperature_at_specific_enthalpy(243_750.0).unwrap(),
        450.0,
        1e-12,
    );
    let initial_volume = law.volume_m3(0.0, 2.0).unwrap();
    let melted_volume = law.volume_m3(175_000.0, 2.0).unwrap();
    close(
        melted_volume - initial_volume,
        2.0 * 0.25 * (1.0 / 800.0 - 1.0 / 1000.0),
        1e-18,
    );
    close(
        law.specific_entropy_change_j_kg_k(125_000.0, 175_000.0)
            .unwrap(),
        50_000.0 / 400.0,
        1e-12,
    );
}

#[test]
fn independent_constituent_datums_leave_the_physics_unchanged() {
    let matrix0 = matrix(0.0);
    let pcm0 = pcm(0.0);
    let shifted_matrix = matrix(1_073_741_824.0);
    let shifted_pcm = pcm(-268_435_456.0);
    let reference = PcmComposite::try_new(&matrix0, &pcm0, 0.25, -2000.0).unwrap();
    let shifted = PcmComposite::try_new(&shifted_matrix, &shifted_pcm, 0.25, -2000.0).unwrap();
    let shifted_composite =
        PcmComposite::try_new(&shifted_matrix, &shifted_pcm, 0.25, 30_000.0).unwrap();
    assert_ne!(reference.identity(), shifted.identity());
    assert_ne!(shifted.identity(), shifted_composite.identity());
    for h in [
        -2000.0, 60_500.0, 123_000.0, 148_000.0, 173_000.0, 241_750.0, 310_500.0,
    ] {
        let a = reference.values_at_specific_enthalpy(h).unwrap();
        let b = shifted.values_at_specific_enthalpy(h).unwrap();
        let c = shifted_composite
            .values_at_specific_enthalpy(h + 32_000.0)
            .unwrap();
        for state in [b, c] {
            close(state.temperature_k, a.temperature_k, 1e-12);
            close(state.bulk_density_kg_m3, a.bulk_density_kg_m3, 1e-10);
            close(state.liquid_mass_fraction, a.liquid_mass_fraction, 1e-15);
        }
    }
    let entropy = reference
        .specific_entropy_change_j_kg_k(-2000.0, 310_500.0)
        .unwrap();
    close(
        shifted
            .specific_entropy_change_j_kg_k(-2000.0, 310_500.0)
            .unwrap(),
        entropy,
        1e-10,
    );
    close(
        shifted_composite
            .specific_entropy_change_j_kg_k(30_000.0, 342_500.0)
            .unwrap(),
        entropy,
        1e-10,
    );
    close(
        entropy,
        1250.0 * (400.0_f64 / 300.0).ln() + 125.0 + 1375.0 * (500.0_f64 / 400.0).ln(),
        1e-10,
    );
}

#[test]
fn source_density_laws_and_all_caloric_breakpoints_survive_preparation() {
    let density_source = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([3; 32]),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 1000.0,
                bulk_density_kg_m3: 8000.0,
            },
            HeatCapacityKnot {
                temperature_k: 500.0,
                specific_heat_capacity_j_kg_k: 1000.0,
                bulk_density_kg_m3: 6000.0,
            },
        ],
        1e-8,
    )
    .unwrap();
    let pcm = pcm(0.0);
    let law = PcmComposite::try_new(&density_source, &pcm, 0.25, 0.0).unwrap();
    let state = law.values_at_specific_enthalpy(150_000.0).unwrap();
    // The original sourced sensible rho(T) is affine density, while the PCM's
    // authored latent coexistence uses additive specific volume.
    close(state.matrix.bulk_density_kg_m3, 7000.0, 1e-12);
    close(
        state.pcm.bulk_density_kg_m3,
        1.0 / (0.5 / 1000.0 + 0.5 / 800.0),
        1e-10,
    );
    close(
        state.bulk_density_kg_m3,
        1.0 / (0.75 / 7000.0 + 0.25 * (0.5 / 1000.0 + 0.5 / 800.0)),
        1e-10,
    );

    let kinked_matrix = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([4; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(0.0, 300.0, 0.0, 8000.0),
            knot(50_000.0, 350.0, 0.0, 8000.0),
            knot(350_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap();
    let kinked = PcmComposite::try_new(&kinked_matrix, &pcm, 0.25, 0.0).unwrap();
    assert_eq!(kinked.knots().len(), 5);
    assert_eq!(kinked.knots()[1].temperature_k, 350.0);
    assert_eq!(kinked.knots()[1].specific_enthalpy_j_kg, 62_500.0);
    close(
        kinked.temperature_at_specific_enthalpy(31_250.0).unwrap(),
        325.0,
        1e-12,
    );
    close(
        kinked.temperature_at_specific_enthalpy(112_500.0).unwrap(),
        375.0,
        1e-12,
    );
    let entropy = kinked
        .specific_entropy_change_j_kg_k(0.0, 425_000.0)
        .unwrap();
    close(
        entropy,
        1250.0 * (350.0_f64 / 300.0).ln()
            + 2000.0 * (400.0_f64 / 350.0).ln()
            + 125.0
            + 2125.0 * (500.0_f64 / 400.0).ln(),
        1e-10,
    );
    close(
        kinked
            .specific_entropy_change_j_kg_k(425_000.0, 0.0)
            .unwrap(),
        -entropy,
        1e-10,
    );
    close(
        kinked
            .specific_entropy_change_j_kg_k(0.0, 187_500.0)
            .unwrap()
            + kinked
                .specific_entropy_change_j_kg_k(187_500.0, 425_000.0)
                .unwrap(),
        entropy,
        1e-10,
    );
}

#[test]
fn incomplete_support_and_ambiguous_matrix_plateaus_refuse() {
    let matrix = matrix(0.0);
    let pcm = pcm(0.0);
    for fraction in [0.0, 1.0, -0.1, f64::NAN] {
        assert!(PcmComposite::try_new(&matrix, &pcm, fraction, 0.0).is_err());
    }
    assert!(PcmComposite::try_new(&matrix, &pcm, 0.25, f64::INFINITY).is_err());
    assert!(PcmComposite::try_new(&matrix, &matrix, 0.25, 0.0).is_err());
    let short_matrix = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([5; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(0.0, 310.0, 0.0, 8000.0),
            knot(190_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap();
    assert!(PcmComposite::try_new(&short_matrix, &pcm, 0.25, 0.0).is_err());
    let ambiguous = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([6; 32]),
        vec![
            knot(0.0, 300.0, 0.0, 8000.0),
            knot(100_000.0, 400.0, 0.0, 8000.0),
            knot(110_000.0, 400.0, 0.0, 8000.0),
            knot(210_000.0, 500.0, 0.0, 8000.0),
            knot(310_000.0, 600.0, 1.0, 7000.0),
        ],
    )
    .unwrap();
    assert!(PcmComposite::try_new(&ambiguous, &pcm, 0.25, 0.0).is_err());
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    assert!(law.temperature_at_specific_enthalpy(-1.0).is_err());
    assert!(law.temperature_at_specific_enthalpy(312_501.0).is_err());
    assert!(law.volume_m3(0.0, 0.0).is_err());
    // A pure isothermal PCM still has a unique allocation because the matrix
    // is sensible and stays at its one uniquely determined enthalpy.
    let isothermal = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([7; 32]),
        vec![
            knot(0.0, 400.0, 0.0, 1000.0),
            knot(200_000.0, 400.0, 1.0, 800.0),
        ],
    )
    .unwrap();
    let plateau = PcmComposite::try_new(&matrix, &isothermal, 0.25, 0.0).unwrap();
    assert_eq!(plateau.maximum_temperature_slope(), 0.0);
    assert_eq!(
        plateau
            .values_at_specific_enthalpy(25_000.0)
            .unwrap()
            .liquid_mass_fraction,
        0.125
    );
}
