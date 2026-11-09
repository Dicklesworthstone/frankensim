//! Analytic original-source liquid fixtures; no material qualification claimed.
use super::*;
use crate::phase::SolidLiquidPhase;
use crate::thermal::transport::isobaric::expansion::liquid::LiquidIsobaricExpansion;

fn liquid_fixture(
    strains: &[(f64, f64)],
) -> (
    MaterialCard,
    IsobaricHeatCapacityLaw,
    PreparedFreeExpansionDensity,
) {
    let keys = keys();
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim(&keys[0], &CP, true)).unwrap();
    let mut rho = claim(&keys[1], &[(300.0, 1000.0), (300.0, 1000.0)], true);
    rho.value = PropertyValue::Scalar {
        value: 1000.0,
        dims: keys[1].quantity().dims(),
    };
    rho.interpolation = InterpolationPolicy::ConstantWithinValidity;
    claims.insert_claim(rho).unwrap();
    claims.insert_claim(claim(&keys[2], strains, true)).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic Cp/free-volume fixture".into(),
            phase: "liquid".into(),
            process: "one maintained isobar".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let (cp, density) = prepare(&card, 1, 290.0, 310.0, 300.0, PRESSURE);
    (card, cp, density)
}

#[test]
fn original_variable_cp_and_cubic_volume_drive_signed_material_dilation() {
    let (card, cp, density) = liquid_fixture(&[(290.0, 0.0), (310.0, 0.02)]);
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let liquid = LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
    assert!(core::ptr::eq(liquid.caloric(), &original));
    let h = cp.specific_enthalpy_j_kg(305.0).unwrap();
    let stretch = 1.015_f64 / 1.01;
    for rate in [-6.0, 0.0, 6.0] {
        let response = liquid.dilation_at_specific_enthalpy(h, rate).unwrap();
        assert_eq!(response.caloric_identity, original.identity());
        assert_eq!(response.state.path_identity(), density.identity());
        close(response.state.temperature_k(), 305.0, 1e-12);
        close(response.state.relative_volume(), stretch.powi(3), 1e-15);
        close(
            response.state.density_kg_m3(),
            1000.0 / stretch.powi(3),
            1e-12,
        );
        close(response.specific_heat_cp_j_kg_k, 3.0, 1e-14);
        close(
            response.specific_volume_derivative_m3_j,
            3.0 * stretch.powi(2) * (0.001 / 1.01) / (1000.0 * 3.0),
            1e-20,
        );
        close(
            response.dilation_per_s,
            3.0 * (0.001 / 1.015) * (rate / 3.0),
            1e-16,
        );
    }
    // Changing only the enthalpy datum changes identity, never the physical response.
    let (other_cp, other_density) = prepare(&card, 1, 290.0, 310.0, 295.0, PRESSURE);
    let other = IsobaricFreeExpansion::try_new(&card, &other_cp, &other_density).unwrap();
    let other_liquid = LiquidIsobaricExpansion::try_new(&other, SolidLiquidPhase::Liquid).unwrap();
    assert_ne!(original.identity(), other.identity());
    close(
        other_liquid
            .dilation_at_specific_enthalpy(other_cp.specific_enthalpy_j_kg(305.0).unwrap(), 6.0)
            .unwrap()
            .dilation_per_s,
        liquid
            .dilation_at_specific_enthalpy(h, 6.0)
            .unwrap()
            .dilation_per_s,
        1e-16,
    );
}

#[test]
fn actual_source_kink_uses_left_cooling_and_right_heating_slopes() {
    let (card, cp, density) = liquid_fixture(&[(290.0, 0.0), (302.0, 0.024), (310.0, 0.016)]);
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let liquid = LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
    let h = cp.specific_enthalpy_j_kg(302.0).unwrap();
    let cooling = liquid.dilation_at_specific_enthalpy(h, -2.4).unwrap();
    let heating = liquid.dilation_at_specific_enthalpy(h, 2.4).unwrap();
    close(cooling.dilation_per_s, -3.0 * 0.002 / 1.024, 1e-16);
    close(heating.dilation_per_s, -3.0 * 0.001 / 1.024, 1e-16);
    assert!(cooling.specific_volume_derivative_m3_j > 0.0);
    assert!(heating.specific_volume_derivative_m3_j < 0.0);
    assert_eq!(cooling.state, heating.state);
    assert!(
        density
            .expansion()
            .at_temperature_with_smooth_tangent(302.0)
            .is_err()
    );

    let [low, high] = cp.specific_enthalpy_span_j_kg();
    assert!(liquid.dilation_at_specific_enthalpy(low, -1.0).is_err());
    assert!(liquid.dilation_at_specific_enthalpy(high, 1.0).is_err());
    assert!(liquid.dilation_at_specific_enthalpy(low, 1.0).is_ok());
    assert!(liquid.dilation_at_specific_enthalpy(high, -1.0).is_ok());
    for h in [low, high] {
        assert_eq!(
            liquid
                .dilation_at_specific_enthalpy(h, 0.0)
                .unwrap()
                .dilation_per_s,
            0.0
        );
    }
    assert!(liquid.dilation_at_specific_enthalpy(h, f64::NAN).is_err());
    assert!(liquid.dilation_at_specific_enthalpy(f64::NAN, 0.0).is_err());
}

#[test]
fn native_segment_tangent_survives_an_interval_without_representable_interior() {
    let (card, _, _) = liquid_fixture(&[(290.0, 0.0), (310.0, 0.02)]);
    let high = f64::from_bits(300.0_f64.to_bits() + 1);
    let (cp, density) = prepare(&card, 1, 300.0, high, 300.0, PRESSURE);
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let liquid = LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
    let [low_h, high_h] = cp.specific_enthalpy_span_j_kg();
    close(
        liquid
            .dilation_at_specific_enthalpy(low_h, 2.0)
            .unwrap()
            .dilation_per_s,
        3.0 * 0.001 / 1.01,
        1e-16,
    );
    assert!(
        liquid
            .dilation_at_specific_enthalpy(high_h, -2.0)
            .unwrap()
            .dilation_per_s
            < 0.0
    );
}

#[test]
fn liquid_declaration_and_isotropic_source_are_both_required() {
    let (card, cp, density) = liquid_fixture(&[(290.0, 0.0), (310.0, 0.02)]);
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    assert!(LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Solid).is_err());
    let (card, cp, density) = fixture(
        &CP,
        &[
            &[(290.0, 0.0), (310.0, 0.02)],
            &[(290.0, 0.0), (310.0, 0.02)],
            &[(290.0, 0.0), (310.0, 0.02)],
        ],
        true,
        300.0,
    );
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    assert!(LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).is_err());
}
