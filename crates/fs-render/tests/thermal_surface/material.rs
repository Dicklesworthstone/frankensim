use super::*;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, LinearizedSurfaceRadiation, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_render::material_state::{
    MaterialOpticalBinding, MaterialOpticalError, MaterialThermalBinding, MaterialThermalError,
};

fn card(epsilon: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: epsilon,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 2500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic gray-surface regression".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "declared gray-diffuse finish".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn state(card: &MaterialCard, t: f64) -> fs_material::state_point::ResolvedMaterialStatePoint {
    state_at(card, &QueryPoint::new().with("T", t).unwrap())
}

fn state_at(
    card: &MaterialCard,
    point: &QueryPoint,
) -> fs_material::state_point::ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        point,
        &[ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn condition(axis: &str, temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            axis,
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
        .with_quantity(
            "liquid_fraction",
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap()
}

fn condition_card(axis: &str) -> MaterialCard {
    use fs_material::state_point::{
        VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
        VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY, VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
        VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES,
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, 0.2),
        (VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, Dims::NONE, 1.5),
        (
            VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
            Dims([2, 0, 0, 0, 0, 0]),
            0.0,
        ),
        (
            VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
            Dims([4, 0, 0, 0, 0, 0]),
            0.0,
        ),
        (
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            Dims([1, 0, 0, 0, 0, 0]),
            1.0,
        ),
        (
            VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[0],
            Dims::NONE,
            1.0,
        ),
        (
            VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[1],
            Dims::NONE,
            1.0,
        ),
        (
            VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[2],
            Dims::NONE,
            1.0,
        ),
    ] {
        let varying = name == SURFACE_EMISSIVITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varying {
                    PropertyValue::Curve {
                        abscissa: axis.into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(1200.0, 0.2), (1800.0, 0.6)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with(axis, 1200.0, 1800.0)
                    .with("pressure", 100_000.0, 200_000.0)
                    .with("liquid_fraction", 0.0, 1.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varying {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic complete-condition optical fixture".into(),
                    license: "internal-test-use".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "fixture".into(),
            process: "declared independent optical laws".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn condition_emissivity(card: &MaterialCard, point: &QueryPoint, axis: &str) -> SurfaceEmissivity {
    SurfaceEmissivity::from_condition(
        "body",
        ThermalMaterialPoint {
            material: card,
            context: point,
            temperature_axis: axis,
            liquid_fraction_axis: Some("liquid_fraction"),
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .0
}

#[test]
fn g3_complete_thermal_conditions_drive_radiance_independently_of_axis_name() {
    with_cx(|cx| {
        let mut reference_pixels = None;
        for axis in ["T", "temperature", "Ts"] {
            let card = condition_card(axis);
            let mut pixels = Vec::new();
            for (temperature, epsilon) in [(1200.0, 0.2), (1800.0, 0.6)] {
                let point = condition(axis, temperature);
                let accepted = state_at(&card, &point);
                let emissivity = condition_emissivity(&card, &point, axis);
                let binding = MaterialThermalBinding::try_gray_diffuse(
                    &accepted,
                    &emissivity,
                    ContentHash([19; 32]),
                )
                .unwrap();
                assert_eq!(binding.accepted_state_identity(), accepted.identity());
                assert_eq!(binding.surface().temperature_k(), temperature);
                let power = epsilon
                    * fs_conduction::radiation::STEFAN_BOLTZMANN_W_M2_K4
                    * temperature.powi(4);
                assert!((binding.surface().bolometric_exitance_w_m2() / power - 1.0).abs() < 1e-13);
                pixels.push(
                    render(
                        &scene(binding.material()),
                        cx,
                        &settings(DirectStrategy::Mis),
                    )
                    .unwrap()
                    .xyz[0],
                );
            }
            assert!(pixels[1][1] > pixels[0][1]);
            if let Some(expected) = &reference_pixels {
                assert_eq!(&pixels, expected);
            } else {
                reference_pixels = Some(pixels);
            }
        }
    });
}

#[test]
fn g0_same_temperature_cannot_hide_changed_pressure_phase_or_quantity_meaning() {
    let card = condition_card("T");
    let point = condition("T", 1800.0);
    let emissivity = condition_emissivity(&card, &point, "T");
    for other in [
        point.clone().with("pressure", 200_000.0).unwrap(),
        point.clone().with("liquid_fraction", 0.75).unwrap(),
        point
            .clone()
            .with_quantity("pressure", quantity(QuantityKind::Stress), 100_000.0)
            .unwrap(),
    ] {
        let accepted = state_at(&card, &other);
        assert!(matches!(
            MaterialThermalBinding::try_gray_diffuse(&accepted, &emissivity, ContentHash([19; 32]),),
            Err(MaterialThermalError::StateMismatch)
        ));
    }
}

#[test]
fn g0_optical_binding_retains_quantity_meaning_at_identical_numeric_coordinates() {
    use fs_material::state_point::resolve_visible_dielectric_state_point;
    let card = condition_card("Ts");
    let point = condition("Ts", 1800.0);
    let optical = resolve_visible_dielectric_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MaterialOpticalBinding::try_dielectric(
        &state_at(&card, &point),
        &optical,
        None,
        ContentHash([19; 32]),
    )
    .unwrap();
    let other = point
        .clone()
        .with_quantity("pressure", quantity(QuantityKind::Stress), 100_000.0)
        .unwrap();
    let accepted = state_at(&card, &other);
    assert_eq!(accepted.query_point(), optical.resolved().query_point());
    assert!(matches!(
        MaterialOpticalBinding::try_dielectric(&accepted, &optical, None, ContentHash([19; 32]),),
        Err(MaterialOpticalError::StateMismatch)
    ));
}

#[test]
fn g3_same_card_and_temperature_drive_cooling_and_actual_rendered_radiance() {
    with_cx(|cx| {
        let mut pixels = Vec::new();
        let mut powers = Vec::new();
        for epsilon in [0.2, 0.6] {
            let source = card(epsilon);
            let accepted = state(&source, 1800.0);
            let emission = SurfaceEmissivity::from_card(
                "body",
                &source,
                1800.0,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
            let binding = MaterialThermalBinding::try_gray_diffuse(
                &accepted,
                &emission,
                ContentHash([9; 32]),
            )
            .unwrap();
            assert_eq!(binding.emissivity(), &emission);
            let cooling =
                LinearizedSurfaceRadiation::new("body", emission.clone(), 1790.0, 1780.0, 30.0)
                    .unwrap()
                    .evaluate(1800.0)
                    .unwrap();
            let ambient =
                GrayDiffuseThermalSurface::try_new(1780.0, epsilon, ContentHash([9; 32])).unwrap();
            let net =
                binding.surface().bolometric_exitance_w_m2() - ambient.bolometric_exitance_w_m2();
            assert!((net / cooling.nonlinear_outward_flux_w_m2 - 1.0).abs() < 1.0e-13);
            powers.push(net);
            pixels.push(
                render(
                    &scene(binding.material()),
                    cx,
                    &settings(DirectStrategy::Mis),
                )
                .unwrap()
                .xyz[0],
            );
            assert!(matches!(
                MaterialThermalBinding::try_gray_diffuse(
                    &state(&source, 1801.0),
                    &emission,
                    ContentHash([9; 32])
                ),
                Err(MaterialThermalError::StateMismatch)
            ));
            assert!(matches!(
                MaterialThermalBinding::try_gray_diffuse(
                    &state(&card(0.7), 1800.0),
                    &emission,
                    ContentHash([9; 32])
                ),
                Err(MaterialThermalError::StateMismatch)
            ));
        }
        assert!((powers[1] / powers[0] - 3.0).abs() < 1.0e-12);
        for channel in 0..3 {
            assert!((pixels[1][channel] / pixels[0][channel] - 3.0).abs() < 1.0e-11);
        }
    });
}
