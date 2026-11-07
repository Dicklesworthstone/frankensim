//! Exact model selectors join original thermal measurements to native histories.
use super::*;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, SelectionPolicy, UncertaintyModel,
};
use fs_material::plastic::finite::thermal::source::tabulated::*;
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS, KIRCHHOFF_YIELD_CURVE_PROPERTY,
};

const CURVE: [(f64, f64, f64); 3] = [(0.0, 2e6, -2e3), (0.05, 3e6, -3e3), (0.2, 4e6, -5e3)];
const IDS: [&str; 4] = [
    THERMAL_FINITE_J2_TABULATED_LAW_ID,
    THERMAL_PERZYNA_FINITE_J2_TABULATED_LAW_ID,
    THERMAL_FINITE_NORTON_J2_TABULATED_LAW_ID,
    THERMAL_FINITE_PRIMARY_NORTON_J2_TABULATED_LAW_ID,
];

fn table_model(id: &str) -> ConstitutiveModelCard {
    let mut model = fixture();
    model.law = LawId(id.into());
    for name in [
        "yield_stress_pa",
        "hardening_modulus_pa",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
    ] {
        model.parameters.remove(name);
    }
    model
        .parameters
        .get_mut("maximum_accumulated_plastic_strain")
        .unwrap()
        .value = 0.2;
    let prefix = if id == IDS[1] {
        Some("overstress")
    } else if id == IDS[2] || id == IDS[3] {
        model.state_schema_version = 2;
        model.parameters.insert(
            "maximum_accumulated_creep_strain".into(),
            LawParameter {
                value: 0.2,
                dims: Dims::NONE,
            },
        );
        Some("creep")
    } else {
        None
    };
    if let Some(prefix) = prefix {
        for (suffix, value, dims) in [
            (
                "reference_viscosity_pa_s",
                VOCE_KINETICS.reference_viscosity_pa_s,
                Dims([-1, 1, -1, 0, 0, 0]),
            ),
            ("stress_scale_pa", VOCE_KINETICS.stress_scale_pa, PA),
            ("exponent", VOCE_KINETICS.exponent, Dims::NONE),
            (
                "activation_temperature_k",
                VOCE_KINETICS.activation_temperature_k,
                K,
            ),
        ] {
            model
                .parameters
                .insert(format!("{prefix}_{suffix}"), LawParameter { value, dims });
        }
    }
    if id == IDS[3] {
        for (name, value) in PRIMARY_PARAMETERS {
            model.parameters.insert(
                name.into(),
                LawParameter {
                    value,
                    dims: Dims::NONE,
                },
            );
        }
    }
    model
}

fn material(model: ConstitutiveModelCard) -> MaterialCard {
    let samples = [T_SPAN[0], T_REF, T_SPAN[1]]
        .into_iter()
        .flat_map(|t| {
            CURVE
                .into_iter()
                .map(move |(alpha, y, slope)| ([alpha, t], y + slope * (t - T_REF)))
        })
        .collect();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(KIRCHHOFF_YIELD_CURVE_PROPERTY, PA),
            value: PropertyValue::Surface2D {
                axes: [ACCUMULATED_PLASTIC_STRAIN_AXIS.into(), "temperature".into()],
                axis_dims: [Dims::NONE, K],
                samples,
                triangles: vec![],
                dims: PA,
            },
            validity: fixture().validity.with_quantity(
                ACCUMULATED_PLASTIC_STRAIN_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.2,
            ),
            interpolation: InterpolationPolicy::TabulatedOnly,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: fixture().provenance,
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal hardening source model".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap()
}

fn admitted(
    card: &MaterialCard,
) -> Result<ResolvedAffineThermalFiniteJ2, fs_material::MaterialError> {
    ResolvedAffineThermalFiniteJ2::try_from_material_card_at(
        card,
        0,
        &context(),
        "temperature",
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        9,
    )
}

fn resistance(alpha: f64, temperature: f64) -> f64 {
    let pair = if alpha <= CURVE[1].0 {
        [CURVE[0], CURVE[1]]
    } else {
        [CURVE[1], CURVE[2]]
    };
    let y = pair.map(|(_, value, slope)| value + slope * (temperature - T_REF));
    y[0] + (y[1] - y[0]) * ((alpha - pair[0].0) / (pair[1].0 - pair[0].0))
}

#[test]
fn each_explicit_table_selector_reaches_native_flow_and_retains_original_evidence() {
    for id in IDS {
        let card = material(table_model(id));
        let family = admitted(&card).unwrap();
        let source = family.tabulated_hardening_source().unwrap();
        assert_eq!(family.model_card(), &card.models()[0]);
        assert_eq!(source.source().card_identity(), card.content_hash());
        assert_eq!(source.source().observations().len(), 9);
        assert_eq!(
            source.source().claim().interpolation,
            InterpolationPolicy::TabulatedOnly
        );
        assert_eq!(family.context(), &context());
        assert!(std::ptr::eq(
            source,
            family.clone().tabulated_hardening_source().unwrap()
        ));
        for (_, observation) in source.source().observations() {
            card.claims().verify_receipt(&observation.receipt).unwrap();
        }
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
                &card.models()[0],
                &context(),
                "temperature"
            )
            .is_err()
        );
        let temperature = 380.0;
        let duration = 0.12;
        let old = family.law().initial_state();
        let response = family
            .law()
            .at_temperature_for_step(temperature, duration)
            .unwrap()
            .response(&deformation(temperature, 0.02, 0.12), &old)
            .unwrap();
        assert!(response.state.alpha > 0.05 && response.state.alpha < 0.2);
        let yield_pa = resistance(response.state.alpha, temperature);
        if id == IDS[1] {
            let eta = VOCE_KINETICS.reference_viscosity_pa_s
                * (VOCE_KINETICS.activation_temperature_k * (temperature.recip() - T_REF.recip()))
                    .exp();
            let overstress = VOCE_KINETICS.stress_scale_pa
                * (eta * response.state.alpha / (duration * VOCE_KINETICS.stress_scale_pa))
                    .powf(1.0 / VOCE_KINETICS.exponent);
            close(response.equivalent_kirchhoff_pa, yield_pa + overstress);
            assert!(family.law().at_temperature(temperature).is_err());
        } else {
            close(response.equivalent_kirchhoff_pa, yield_pa);
        }
        if id == IDS[2] || id == IDS[3] {
            assert!(response.state.creep_alpha > 0.0);
            assert_eq!(family.law().maximum_accumulated_creep_strain(), Some(0.2));
        } else {
            assert_eq!(response.state.creep_alpha, 0.0);
        }
        assert_eq!(
            family.law().norton_strain_hardening(),
            if id == IDS[3] { Some(&PRIMARY) } else { None }
        );
        assert_eq!(old.alpha, 0.0);
        assert_eq!(old.creep_alpha, 0.0);
    }
}

#[test]
fn source_table_model_requires_exact_selection_parameters_state_and_support() {
    for id in IDS {
        let original = table_model(id);
        for mutation in 0..5 {
            let mut model = original.clone();
            match mutation {
                0 => {
                    model.parameters.remove("shear_modulus_slope_pa_k");
                }
                1 => {
                    model.parameters.insert(
                        "yield_stress_pa".into(),
                        LawParameter {
                            value: 2e6,
                            dims: PA,
                        },
                    );
                }
                2 => {
                    model
                        .parameters
                        .get_mut("maximum_accumulated_plastic_strain")
                        .unwrap()
                        .value = 0.21;
                }
                3 => {
                    model.state_schema_version = if model.state_schema_version == 1 {
                        2
                    } else {
                        1
                    };
                }
                _ => {
                    model.law = LawId(THERMAL_FINITE_J2_LAW_ID.into());
                }
            }
            assert!(
                admitted(&material(model)).is_err(),
                "selector {id}, mutation {mutation}"
            );
        }
        let card = material(original);
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_material_card_at(
                &card,
                0,
                &context(),
                "temperature",
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                8,
            )
            .is_err()
        );
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_material_card_at(
                &card,
                1,
                &context(),
                "temperature",
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                9,
            )
            .is_err()
        );
        assert!(
            ResolvedAffineThermalFiniteJ2::try_from_material_card_at(
                &card,
                0,
                &context().with("temperature", 301.0).unwrap(),
                "temperature",
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                9,
            )
            .is_err()
        );
    }
}
