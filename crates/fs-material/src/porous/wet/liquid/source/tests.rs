//! Actual synthetic observation-backed cards; no claim of empirical validation.
use super::*;
use crate::porous::permeability::{
    DarcyGasBarrier, INTRINSIC_PERMEABILITY_DIMS, IntrinsicPermeability,
};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};

const K: &str = "intrinsic_permeability";
const KEYS: WetPoreLiquidProperties<'static> = WetPoreLiquidProperties {
    residual_liquid_saturation: "residual_liquid_pore_fraction",
    residual_gas_saturation: "residual_gas_pore_fraction",
    van_genuchten_m: "retention_m",
    connectivity_exponent: "connectivity_ell",
    pressure_scale_pa: "capillary_scale_pa",
};
#[derive(Clone, Copy)]
enum Variant {
    Constant,
    PressureHead,
    Curve,
    Moisture(MoistureBasis),
    CapillaryAxis,
    MinimumSpeed,
}
fn support() -> WetPoreLiquidSupport {
    WetPoreLiquidSupport {
        liquid_saturation: (0.0, 1.0),
        temperature_k: (270.0, 360.0),
        liquid_pressure_pa: (10000.0, 1e6),
        maximum_superficial_speed_m_s: 1.0,
    }
}
fn saturation_quantity(variant: Variant) -> QuantitySpec {
    match variant {
        Variant::Moisture(basis) => QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(basis),
            ValueForm::Static,
        )),
        _ => QuantitySpec::dimensional(Dims::NONE),
    }
}
fn fixture(ell: f64, variant: Variant) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let provenance = Provenance {
        source: "synthetic joint intrinsic/retention observation fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic porous plug with stated wetting liquid".into(),
            method: "synthetic fixed-branch VG–Mualem fit".into(),
            artifact: fs_blake3::hash_domain("fixture", b"synthetic liquid retention calibration"),
            caveats: "synthetic source-retention fixture, not material measurements".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (i, (name, value)) in [
        (K, 2e-10),
        (KEYS.residual_liquid_saturation, 0.1),
        (KEYS.residual_gas_saturation, 0.1),
        (KEYS.van_genuchten_m, 0.5),
        (KEYS.connectivity_exponent, ell),
        (KEYS.pressure_scale_pa, 800.0),
    ]
    .into_iter()
    .enumerate()
    {
        let dims = if i == 0 {
            INTRINSIC_PERMEABILITY_DIMS
        } else if i == 5 {
            if matches!(variant, Variant::PressureHead) {
                fs_qty::Length::DIMS
            } else {
                fs_qty::Pressure::DIMS
            }
        } else {
            Dims::NONE
        };
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
                290.0,
                320.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                120000.0,
                180000.0,
            );
        if i > 0 {
            validity = validity
                .with_quantity(
                    WetPoreLiquidLaw::SOURCE_SATURATION_AXIS,
                    saturation_quantity(variant),
                    0.2,
                    0.8,
                )
                .with_quantity(
                    WetPoreLiquidLaw::SOURCE_SPEED_AXIS,
                    QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
                    if matches!(variant, Variant::MinimumSpeed) {
                        0.05
                    } else {
                        0.0
                    },
                    0.2,
                );
            if matches!(variant, Variant::CapillaryAxis) {
                validity = validity.with_quantity(
                    "capillary_pressure",
                    QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                    500.0,
                    2000.0,
                );
            }
        }
        let curved = i == 5 && matches!(variant, Variant::Curve);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curved {
                    PropertyValue::Curve {
                        abscissa: "temperature".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(290.0, 800.0), (320.0, 1600.0)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity,
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.03,
                    confidence: 0.95,
                },
                interpolation: if curved {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![observation],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic porous solid".into(),
            phase: "solid".into(),
            process: "synthetic wetting-liquid characterization".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard, variant: Variant, temperature: f64) -> ResolvedMaterialStatePoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            150000.0,
        )
        .unwrap()
        .with_quantity(
            WetPoreLiquidLaw::SOURCE_SATURATION_AXIS,
            saturation_quantity(variant),
            0.58,
        )
        .unwrap()
        .with_quantity(
            WetPoreLiquidLaw::SOURCE_SPEED_AXIS,
            QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
            0.1,
        )
        .unwrap();
    if matches!(variant, Variant::CapillaryAxis) {
        point = point
            .with_quantity(
                "capillary_pressure",
                QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                1000.0,
            )
            .unwrap();
    }
    let keys = [
        K,
        KEYS.residual_liquid_saturation,
        KEYS.residual_gas_saturation,
        KEYS.van_genuchten_m,
        KEYS.connectivity_exponent,
        KEYS.pressure_scale_pa,
    ];
    let requirements: Vec<_> = keys
        .into_iter()
        .enumerate()
        .map(|(i, key)| {
            let dims = if i == 0 {
                INTRINSIC_PERMEABILITY_DIMS
            } else if i == 5 {
                if matches!(variant, Variant::PressureHead) {
                    fs_qty::Length::DIMS
                } else {
                    fs_qty::Pressure::DIMS
                }
            } else {
                Dims::NONE
            };
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::NonNegative).unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::PreferObservationBacked,
    )
    .unwrap()
}
fn barrier(source: ResolvedMaterialStatePoint) -> DarcyGasBarrier {
    DarcyGasBarrier::try_new(
        vec![PorousGasLayer::from_material_isotropic(source, K, 0.01, None).unwrap()],
        [1.0, 0.0, 0.0],
        support().flow_support(),
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-12 * a.abs().max(b.abs()), "{a} != {b}");
}

#[test]
fn selected_source_coefficients_change_original_liquid_resistance_and_keep_support() {
    for ell in [0.5, 1.5] {
        let card = fixture(ell, Variant::Constant);
        let source = resolve(&card, Variant::Constant, 300.0);
        let law = WetPoreLiquidLaw::from_material(source.clone(), KEYS, support()).unwrap();
        assert_eq!(law.source(), Some(&source));
        assert_eq!(law.source_properties(), Some(KEYS));
        assert_eq!(law.support().liquid_saturation, (0.2, 0.8));
        assert!(law.sample(0.1).is_err());
        let states = [law.sample(0.58).unwrap()];
        near(
            states[0].liquid_relative_permeability(),
            0.6_f64.powf(ell) * 0.04,
        );
        assert_eq!(states[0].law().source(), Some(&source));
        let original = barrier(source.clone());
        let geometry = original
            .continuum_geometry()
            .unwrap()
            .with_liquid_states(&states)
            .unwrap();
        near(
            geometry.geometric_resistance_per_m().unwrap(),
            0.01 / (2e-10 * 0.6_f64.powf(ell) * 0.04),
        );
        assert_eq!(geometry.barrier().layers()[0].source(), Some(&source));
        geometry
            .admit_flow(310.0, 160000.0, 140000.0, -0.15)
            .unwrap();
        for (t, p, v) in [
            (289.0, 150000.0, 0.1),
            (300.0, 180001.0, 0.1),
            (300.0, 150000.0, 0.201),
        ] {
            assert!(geometry.admit_flow(t, p, 140000.0, v).is_err());
        }
    }
}

#[test]
fn wrong_quantity_moisture_basis_curve_and_unowned_pressure_do_not_become_retention() {
    for variant in [
        Variant::PressureHead,
        Variant::Curve,
        Variant::Moisture(MoistureBasis::DryMass),
        Variant::Moisture(MoistureBasis::RelativeHumidity),
        Variant::CapillaryAxis,
        Variant::MinimumSpeed,
    ] {
        let card = fixture(0.5, variant);
        let source = resolve(&card, variant, 300.0);
        assert!(WetPoreLiquidLaw::from_material(source, KEYS, support()).is_err());
    }
    let source = resolve(&fixture(0.5, Variant::Constant), Variant::Constant, 300.0);
    let missing = WetPoreLiquidProperties {
        pressure_scale_pa: "unselected_pressure_scale",
        ..KEYS
    };
    assert!(WetPoreLiquidLaw::from_material(source, missing, support()).is_err());
}

#[test]
fn sourced_retention_cannot_attach_to_another_material_or_reference_context() {
    let card = fixture(0.5, Variant::Constant);
    let source = resolve(&card, Variant::Constant, 300.0);
    let law = WetPoreLiquidLaw::from_material(source.clone(), KEYS, support()).unwrap();
    let states = [law.sample(0.6).unwrap()];
    let different = barrier(resolve(
        &fixture(1.5, Variant::Constant),
        Variant::Constant,
        300.0,
    ));
    let other_context = barrier(resolve(&card, Variant::Constant, 310.0));
    let authored = DarcyGasBarrier::try_new(
        vec![
            PorousGasLayer::try_new(
                IntrinsicPermeability::isotropic(2e-10).unwrap(),
                0.01,
                0.0,
                "authored unrelated skeleton",
            )
            .unwrap(),
        ],
        [1.0, 0.0, 0.0],
        support().flow_support(),
    )
    .unwrap();
    for original in [&different, &other_context, &authored] {
        assert!(
            original
                .continuum_geometry()
                .unwrap()
                .with_liquid_states(&states)
                .is_err()
        );
    }
    // A changed prescribed saturation keeps the same original calibration context.
    let original = barrier(source);
    assert!(
        original
            .continuum_geometry()
            .unwrap()
            .with_liquid_states(&states)
            .is_ok()
    );
}
