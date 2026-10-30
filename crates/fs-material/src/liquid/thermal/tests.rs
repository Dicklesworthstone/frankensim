use super::*;
use crate::fluid_transport::{FluidTransport, ThermalFluidTransport};
use crate::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, VOLUMETRIC_EXPANSION_DIMS,
    VOLUMETRIC_EXPANSION_PROPERTY,
};
use crate::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use crate::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Pressure, QuantitySpec};

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), temperature)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap()
}

fn card(
    density: PropertyValue,
    interpolation: InterpolationPolicy,
    thermal_support: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let density_validity = if thermal_support {
        ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, 400.0)
    } else {
        ValidityDomain::unconstrained()
    };
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 0.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            4000.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.5,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.001),
    ] {
        let is_density = name == DENSITY_PROPERTY;
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if is_density {
                    density.clone()
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: if is_density {
                    density_validity.clone()
                } else {
                    ValidityDomain::unconstrained().with_quantity("T", absolute(), 300.0, 400.0)
                }
                .with_quantity(
                    "p",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                    100000.0,
                    200000.0,
                ),
                interpolation: if is_density {
                    interpolation
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic density/transport law; not observations".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn curve(knots: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: knots.to_vec(),
        dims: Density::DIMS,
    }
}

fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &point(300.0, 100000.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g1_density_tangent_changes_with_source_segments_and_uses_current_density() {
    let card = card(
        curve(&[(300.0, 1000.0), (350.0, 950.0), (400.0, 850.0)]),
        InterpolationPolicy::LinearInside,
        true,
    );
    let law = law(&card)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
        .unwrap();
    for (t, rho, slope) in [
        (300.0, 1000.0, -1.0),
        (325.0, 975.0, -1.0),
        (375.0, 900.0, -2.0),
        (400.0, 850.0, -2.0),
    ] {
        let state = law.at_temperature(Temperature::new(t)).unwrap();
        let expansion = state.density_expansion().unwrap();
        assert_eq!(state.density_kg_m3(), rho);
        assert_eq!(expansion.density_slope_kg_m3_k(), slope);
        assert_eq!(expansion.volumetric_expansion_per_k(), -slope / rho);
        assert_eq!(state.volumetric_expansion_per_k(), Some(-slope / rho));
        assert_eq!(expansion.source().knots.len(), 3);
        let thermal = state.thermal_transport_state().unwrap();
        assert_eq!(thermal.temperature_k(), t);
        assert_eq!(thermal.transport(), state.transport_state().unwrap());
        assert_eq!(thermal.volumetric_expansion_per_k(), Some(-slope / rho));
        let selected = state
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        assert!(
            expansion
                .source()
                .knots
                .iter()
                .all(|(_, sample)| sample.receipt.selected == selected)
        );
    }
    assert!(law.at_temperature(Temperature::new(350.0)).is_err());
    assert!(law.at_temperature(Temperature::new(350.0 + 1e-6)).is_ok());
}

#[test]
fn g1_local_volume_derivative_matches_the_same_density_law() {
    let card = card(
        curve(&[(300.0, 1000.0), (400.0, 900.0)]),
        InterpolationPolicy::LinearInside,
        true,
    );
    let law = law(&card)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
        .unwrap();
    let state = law.at_temperature(Temperature::new(340.0)).unwrap();
    let step = 0.01;
    let volume = |t| {
        2.0 / law
            .at_temperature(Temperature::new(t))
            .unwrap()
            .density_kg_m3()
    };
    let numerical = (volume(340.0 + step) - volume(340.0 - step)) / (2.0 * step * volume(340.0));
    let exact = state
        .density_expansion()
        .unwrap()
        .volumetric_expansion_per_k();
    assert!((numerical / exact - 1.0).abs() < 2e-9);
}

#[test]
fn g1_clipped_curve_interval_retains_the_native_tangent_when_density_rounds_equal() {
    let card = card(
        curve(&[(300.0, 1000.0), (400.0, 900.0)]),
        InterpolationPolicy::LinearInside,
        true,
    );
    let adjacent = f64::from_bits(300.0_f64.to_bits() + 1);
    let prepared = law(&card)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(adjacent))
        .unwrap();
    let state = prepared.at_temperature(Temperature::new(adjacent)).unwrap();
    let expansion = state.density_expansion().unwrap();
    let samples = &expansion.source().knots;
    assert_eq!(
        samples[0].1.evidence.value.value,
        samples[1].1.evidence.value.value
    );
    assert_eq!(expansion.density_slope_kg_m3_k(), -1.0);
    assert_eq!(
        expansion.volumetric_expansion_per_k(),
        1.0 / state.density_kg_m3()
    );
}

#[test]
fn g3_zero_and_negative_expansion_are_source_properties_not_gas_defaults() {
    for (a, b, expected_sign) in [(1000.0, 1000.0, 0.0), (1000.0, 1100.0, -1.0)] {
        let card = card(
            curve(&[(300.0, a), (400.0, b)]),
            InterpolationPolicy::LinearInside,
            true,
        );
        let law = law(&card)
            .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
            .unwrap();
        let value = law
            .at_temperature(Temperature::new(350.0))
            .unwrap()
            .thermal_transport_state()
            .unwrap()
            .volumetric_expansion_per_k()
            .unwrap();
        if expected_sign == 0.0 {
            assert_eq!(value, 0.0);
        } else {
            assert!(value < 0.0);
        }
    }
    let card = card(
        PropertyValue::Scalar {
            value: 1000.0,
            dims: Density::DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        true,
    );
    let liquid = law(&card)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
        .unwrap();
    assert_eq!(
        liquid
            .reference_state()
            .thermal_transport_state()
            .unwrap()
            .volumetric_expansion_per_k(),
        Some(0.0)
    );
}

#[test]
fn g0_missing_temperature_support_and_discrete_samples_cannot_invent_expansion() {
    let unqualified = card(
        PropertyValue::Scalar {
            value: 1000.0,
            dims: Density::DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        false,
    );
    let law = law(&unqualified);
    assert_eq!(
        law.reference_state()
            .thermal_transport_state()
            .unwrap()
            .volumetric_expansion_per_k(),
        None
    );
    assert!(
        law.with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
            .is_err()
    );
    let discrete = card(
        curve(&[(300.0, 1000.0), (400.0, 900.0)]),
        InterpolationPolicy::TabulatedOnly,
        true,
    );
    assert!(
        self::law(&discrete)
            .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
            .is_err()
    );
}

#[test]
fn g0_refresh_preserves_fixed_context_interval_and_last_valid_state() {
    let card = card(
        curve(&[(300.0, 1000.0), (400.0, 900.0)]),
        InterpolationPolicy::LinearInside,
        true,
    );
    let law = law(&card)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(380.0))
        .unwrap();
    let original = law.at_temperature(Temperature::new(340.0)).unwrap();
    for target in [point(390.0, 100000.0), point(340.0, 120000.0)] {
        assert!(original.material().at_query_point(&card, &target).is_ok());
        assert!(original.at_query_point(&card, &target).is_err());
    }
    let updated = original
        .at_query_point(&card, &point(360.0, 100000.0))
        .unwrap();
    assert_eq!(
        updated
            .density_expansion()
            .unwrap()
            .volumetric_expansion_per_k(),
        1.0 / 940.0
    );
    assert_eq!(
        original,
        law.at_temperature(Temperature::new(340.0)).unwrap()
    );
    for (lo, hi) in [
        (300.0, 300.0),
        (400.0, 300.0),
        (299.0, 400.0),
        (300.0, 401.0),
        (f64::NAN, 400.0),
    ] {
        assert!(
            self::law(&card)
                .with_density_expansion(Temperature::new(lo), Temperature::new(hi))
                .is_err()
        );
    }
}

#[test]
fn g1_two_coordinate_density_preserves_pressure_and_source_cell_support() {
    let density = PropertyValue::Surface2D {
        axes: ["T".into(), "p".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: vec![
            ([300.0, 100000.0], 1000.0),
            ([400.0, 100000.0], 900.0),
            ([300.0, 200000.0], 1100.0),
        ],
        triangles: vec![[0, 1, 2]],
        dims: Density::DIMS,
    };
    let card = card(density, InterpolationPolicy::LinearInside, true);
    let initial = LiquidTemperatureLaw::try_new(
        &card,
        &point(300.0, 150000.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let admitted = initial
        .clone()
        .with_density_expansion(Temperature::new(300.0), Temperature::new(350.0))
        .unwrap();
    let warm = admitted.at_temperature(Temperature::new(325.0)).unwrap();
    assert_eq!(warm.density_kg_m3(), 1025.0);
    assert_eq!(
        warm.density_expansion().unwrap().density_slope_kg_m3_k(),
        -1.0
    );
    assert_eq!(
        warm.density_expansion()
            .unwrap()
            .volumetric_expansion_per_k(),
        1.0 / 1025.0
    );
    // Clipping a native source cell must not turn its analytic gradient into
    // a subtractive difference of almost equal evaluated densities.
    let narrow = initial
        .clone()
        .with_density_expansion(Temperature::new(300.0), Temperature::new(300.000001))
        .unwrap();
    assert_eq!(
        narrow
            .reference_state()
            .density_expansion()
            .unwrap()
            .density_slope_kg_m3_k(),
        -1.0
    );
    // Without a representable interior point, matdb cannot select the native
    // cell unambiguously; a rounded endpoint is not an inferred zero tangent.
    assert!(
        initial
            .clone()
            .with_density_expansion(
                Temperature::new(300.0),
                Temperature::new(f64::from_bits(300.0_f64.to_bits() + 1)),
            )
            .is_err()
    );
    // A validity box cannot authorize the unmeasured upper-right triangle.
    assert!(
        initial
            .with_density_expansion(Temperature::new(300.0), Temperature::new(375.0))
            .is_err()
    );
}

#[test]
fn g0_declared_expansion_cannot_be_replaced_by_a_density_tangent() {
    let original = card(
        curve(&[(300.0, 1000.0), (400.0, 900.0)]),
        InterpolationPolicy::LinearInside,
        true,
    );
    for beta in [0.002, 0.0, -0.001] {
        let mut claims = original.claims().clone();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(
                    VOLUMETRIC_EXPANSION_PROPERTY,
                    QuantitySpec::dimensional(VOLUMETRIC_EXPANSION_DIMS),
                ),
                value: PropertyValue::Scalar {
                    value: beta,
                    dims: VOLUMETRIC_EXPANSION_DIMS,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute(),
                    300.0,
                    400.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic independent expansion claim; not observations".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test liquid".into(),
                phase: "liquid".into(),
                process: "independent expansion fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let law = LiquidTemperatureLaw::try_new_with_expansion(
            &card,
            &point(300.0, 100000.0),
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let state = law.reference_state();
        assert_eq!(state.volumetric_expansion_per_k(), Some(beta));
        assert_eq!(
            state
                .thermal_transport_state()
                .unwrap()
                .volumetric_expansion_per_k(),
            Some(beta),
        );
        assert!(
            state
                .clone()
                .with_density_expansion(
                    &card,
                    "T",
                    Temperature::new(300.0),
                    Temperature::new(400.0),
                )
                .is_err()
        );
        assert!(
            law.with_density_expansion(Temperature::new(300.0), Temperature::new(400.0),)
                .is_err()
        );
    }
}
