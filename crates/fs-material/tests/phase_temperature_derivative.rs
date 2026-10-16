//! G3: the derivative follows the admitted piecewise enthalpy chart exactly.
use fs_blake3::ContentHash;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase,
};

fn knot(h: f64, temperature: f64, liquid: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: temperature,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: 1000.0,
    }
}

fn melting_curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x93; 32]),
        vec![
            knot(-100.0, 200.0, 0.0),
            knot(100.0, 300.0, 0.0),
            knot(500.0, 300.0, 1.0),
            knot(1300.0, 500.0, 1.0),
        ],
    )
    .unwrap()
}

fn slope(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> f64 {
    curve
        .temperature_derivative_at_specific_enthalpy(h)
        .unwrap()
}

#[test]
fn solid_latent_and_liquid_slopes_match_both_directional_differences() {
    let curve = melting_curve();
    for (enthalpy, expected) in [(0.0, 0.5), (300.0, 0.0), (900.0, 0.25)] {
        let derivative = slope(&curve, enthalpy);
        assert_eq!(derivative, expected);
        let temperature = curve
            .state_at_specific_enthalpy(enthalpy)
            .unwrap()
            .temperature_k();
        for increment in [-8.0, 8.0] {
            let shifted = curve
                .state_at_specific_enthalpy(enthalpy + increment)
                .unwrap()
                .temperature_k();
            assert_eq!(derivative, (shifted - temperature) / increment);
        }
    }
    for enthalpy in [101.0, 299.0, 499.0] {
        assert_eq!(slope(&curve, enthalpy).to_bits(), 0.0f64.to_bits());
        assert_eq!(
            curve
                .state_at_specific_enthalpy(enthalpy)
                .unwrap()
                .temperature_k(),
            300.0
        );
    }
}

#[test]
fn internal_knots_use_the_right_slope_and_upper_endpoint_uses_the_left_slope() {
    let curve = melting_curve();
    for (enthalpy, expected) in [
        (-100.0, 0.5),
        (100.0f64.next_down(), 0.5),
        (100.0, 0.0),
        (500.0f64.next_down(), 0.0),
        (500.0, 0.25),
        (1300.0, 0.25),
    ] {
        assert_eq!(slope(&curve, enthalpy), expected);
    }
    for (phase, liquid) in [
        (SolidLiquidPhase::Solid, 0.0),
        (SolidLiquidPhase::Liquid, 1.0),
    ] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x94; 32]),
            phase,
            vec![
                knot(-20.0, 200.0, liquid),
                knot(40.0, 230.0, liquid),
                knot(100.0, 320.0, liquid),
            ],
        )
        .unwrap();
        for (enthalpy, expected) in [
            (-20.0, 0.5),
            (10.0, 0.5),
            (40.0, 1.5),
            (70.0, 1.5),
            (100.0, 1.5),
        ] {
            assert_eq!(slope(&curve, enthalpy), expected);
        }
    }
}

#[test]
fn nonfinite_and_outside_queries_return_the_existing_domain_refusals() {
    let curve = melting_curve();
    for enthalpy in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            curve.temperature_derivative_at_specific_enthalpy(enthalpy),
            Err(PhaseStateError::NonFiniteSpecificEnthalpy)
        );
    }
    for enthalpy in [-101.0, 1301.0] {
        assert_eq!(
            curve.temperature_derivative_at_specific_enthalpy(enthalpy),
            Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg: enthalpy,
                lower_j_kg: -100.0,
                upper_j_kg: 1300.0,
            })
        );
    }
}

#[test]
fn unrepresentable_positive_slopes_refuse_without_erasing_plateaus_or_subnormals() {
    let smallest = f64::from_bits(1);
    let single = |h1, t0, t1| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x95; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, t0, 0.0), knot(h1, t1, 0.0)],
        )
        .unwrap()
    };
    let overflow = single(smallest, 1.0, 2.0);
    let underflow = single(f64::MAX, smallest, f64::from_bits(2));
    for (curve, upper) in [(&overflow, smallest), (&underflow, f64::MAX)] {
        for enthalpy in [0.0, upper] {
            assert_eq!(
                curve.temperature_derivative_at_specific_enthalpy(enthalpy),
                Err(PhaseStateError::UnrepresentableTemperatureDerivative)
            );
        }
    }
    let subnormal = single(1.0, smallest, f64::from_bits(2));
    assert_eq!(slope(&subnormal, 0.5).to_bits(), smallest.to_bits());
    let plateau = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x96; 32]),
        vec![knot(0.0, 1.0, 0.0), knot(smallest, 1.0, 1.0)],
    )
    .unwrap();
    assert_eq!(slope(&plateau, 0.0).to_bits(), 0.0f64.to_bits());
}

#[test]
fn signed_zero_queries_resolve_the_same_boundary_state_and_right_slope() {
    for lower in [-0.0, 0.0] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x97; 32]),
            vec![knot(lower, 300.0, 0.0), knot(100.0, 350.0, 1.0)],
        )
        .unwrap();
        let negative = curve.state_at_specific_enthalpy(-0.0).unwrap();
        let positive = curve.state_at_specific_enthalpy(0.0).unwrap();
        assert_eq!(negative.temperature_k(), positive.temperature_k());
        assert_eq!(negative.temperature_k(), 300.0);
        assert_eq!(
            negative.solid_mass_fraction(),
            positive.solid_mass_fraction()
        );
        assert_eq!(
            negative.liquid_mass_fraction(),
            positive.liquid_mass_fraction()
        );
        assert_eq!(negative.bulk_density_kg_m3(), positive.bulk_density_kg_m3());
        assert_eq!(negative.phase(), positive.phase());
        assert_eq!(slope(&curve, -0.0), 0.5);
        assert_eq!(slope(&curve, 0.0), 0.5);
    }
}
