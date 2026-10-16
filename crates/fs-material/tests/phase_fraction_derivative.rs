//! G3: liquid-fraction slopes follow the admitted enthalpy chart without smoothing.
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

fn curve(knots: Vec<EnthalpyPhaseKnot>) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x98; 32]), knots).unwrap()
}

fn melting_curve() -> EquilibriumEnthalpyPhaseCurve {
    curve(vec![
        knot(-64.0, 200.0, 0.0),
        knot(0.0, 300.0, 0.0),
        knot(64.0, 300.0, 0.25),
        knot(128.0, 300.0, 1.0),
        knot(256.0, 500.0, 1.0),
    ])
}

fn slope(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> Result<f64, PhaseStateError> {
    curve.liquid_mass_fraction_derivative_at_specific_enthalpy(h)
}

#[test]
fn segment_slopes_match_fraction_differences_even_on_isothermal_latent_spans() {
    let curve = melting_curve();
    for (h, expected) in [
        (-32.0, 0.0),
        (32.0, 1.0 / 256.0),
        (96.0, 3.0 / 256.0),
        (192.0, 0.0),
    ] {
        let derivative = slope(&curve, h).unwrap();
        assert_eq!(derivative, expected);
        let fraction = curve
            .state_at_specific_enthalpy(h)
            .unwrap()
            .liquid_mass_fraction();
        for increment in [-8.0, 8.0] {
            let shifted = curve
                .state_at_specific_enthalpy(h + increment)
                .unwrap()
                .liquid_mass_fraction();
            assert_eq!(derivative, (shifted - fraction) / increment);
        }
        if expected > 0.0 {
            assert_eq!(
                curve.temperature_derivative_at_specific_enthalpy(h),
                Ok(0.0)
            );
            assert_eq!(
                curve.state_at_specific_enthalpy(h).unwrap().temperature_k(),
                300.0
            );
        }
    }
}

#[test]
fn interior_knots_use_the_right_segment_and_the_upper_endpoint_uses_the_left() {
    let melting = melting_curve();
    for (h, expected) in [
        (-64.0, 0.0),
        (0.0_f64.next_down(), 0.0),
        (-0.0, 1.0 / 256.0),
        (0.0, 1.0 / 256.0),
        (64.0_f64.next_down(), 1.0 / 256.0),
        (64.0, 3.0 / 256.0),
        (128.0_f64.next_down(), 3.0 / 256.0),
        (128.0, 0.0),
        (256.0, 0.0),
    ] {
        assert_eq!(slope(&melting, h), Ok(expected));
    }
    let plateau = curve(vec![knot(0.0, 300.0, 0.0), knot(128.0, 300.0, 1.0)]);
    assert_eq!(slope(&plateau, 128.0), Ok(1.0 / 128.0));
}

#[test]
fn sensible_single_phase_curves_have_exact_zero_fraction_slopes() {
    for (phase, liquid) in [
        (SolidLiquidPhase::Solid, 0.0),
        (SolidLiquidPhase::Liquid, 1.0),
    ] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x99; 32]),
            phase,
            vec![
                knot(-20.0, 200.0, liquid),
                knot(40.0, 230.0, liquid),
                knot(100.0, 320.0, liquid),
            ],
        )
        .unwrap();
        for h in [-20.0, 10.0, 40.0, 70.0, 100.0] {
            assert_eq!(slope(&curve, h).unwrap().to_bits(), 0.0_f64.to_bits());
            assert!(
                curve
                    .temperature_derivative_at_specific_enthalpy(h)
                    .unwrap()
                    > 0.0
            );
        }
    }
}

#[test]
fn nonfinite_and_outside_enthalpies_return_domain_refusals() {
    let curve = melting_curve();
    for h in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            slope(&curve, h),
            Err(PhaseStateError::NonFiniteSpecificEnthalpy)
        );
    }
    for h in [(-64.0_f64).next_down(), 256.0_f64.next_up()] {
        assert_eq!(
            slope(&curve, h),
            Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg: h,
                lower_j_kg: -64.0,
                upper_j_kg: 256.0,
            })
        );
    }
}

#[test]
fn unrepresentable_positive_slopes_refuse_and_representable_subnormals_survive() {
    let smallest = f64::from_bits(1);
    let overflow = curve(vec![knot(0.0, 300.0, 0.0), knot(smallest, 300.0, 1.0)]);
    for h in [0.0, smallest] {
        assert_eq!(
            slope(&overflow, h),
            Err(PhaseStateError::UnrepresentableLiquidMassFractionDerivative)
        );
        assert_eq!(
            overflow.temperature_derivative_at_specific_enthalpy(h),
            Ok(0.0)
        );
    }
    let huge = f64::MAX / 2.0;
    let underflow = curve(vec![
        knot(0.0, 300.0, 0.0),
        knot(huge, 300.0, smallest),
        knot(f64::MAX, 300.0, 1.0),
    ]);
    for h in [0.0, huge / 2.0] {
        assert_eq!(
            slope(&underflow, h),
            Err(PhaseStateError::UnrepresentableLiquidMassFractionDerivative)
        );
    }
    let subnormal = curve(vec![
        knot(0.0, 300.0, 0.0),
        knot(1.0, 300.0, smallest),
        knot(2.0, 300.0, 1.0),
    ]);
    assert_eq!(
        slope(&subnormal, 0.5).unwrap().to_bits(),
        smallest.to_bits()
    );
}
