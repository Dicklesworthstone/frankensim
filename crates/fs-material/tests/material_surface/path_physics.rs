//! Independent temperature-path references through authored two-coordinate data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    MaterialStatePointError, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_material::thermoelectric::{
    ABSOLUTE_SEEBECK_PROPERTY as SEEBECK, ResolvedSeebeckSpan, SEEBECK_DIMS, ThermoelectricError,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

fn point(axis: &str, t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(axis, QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            pressure,
        )
        .unwrap()
}

fn card(
    axis: &str,
    property: &str,
    dims: Dims,
    last: f64,
    gap: bool,
    policy: InterpolationPolicy,
) -> MaterialCard {
    let (samples, triangles) = if gap {
        (
            vec![
                ([300.0, 1e5], 1e-6),
                ([350.0, 1e5], 2e-6),
                ([300.0, 2e5], 2e-6),
                ([350.0, 2e5], 3e-6),
                ([450.0, 1e5], 3e-6),
                ([500.0, 1e5], 4e-6),
                ([450.0, 2e5], 5e-6),
                ([500.0, 2e5], 6e-6),
            ],
            vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
        )
    } else {
        // At p=1.5e5 the line crosses a source edge at T=400:
        // values are 1.5, 2.5, (3+last)/2 micro-units at T=300,400,500.
        (
            vec![
                ([300.0, 1e5], 1e-6),
                ([500.0, 1e5], 3e-6),
                ([300.0, 2e5], 2e-6),
                ([500.0, 2e5], last * 1e-6),
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        )
    };
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            axis,
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            500.0,
        )
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            1e5,
            2e5,
        );
    let provenance = Provenance {
        source: "authored path-integration reference; not measured data".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(property, dims),
            value: PropertyValue::Surface2D {
                axes: [axis.into(), "pressure".into()],
                axis_dims: [Temperature::DIMS, Pressure::DIMS],
                samples,
                triangles,
                dims,
            },
            validity: validity.clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: policy,
            observations: vec![],
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 8000.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 200e9),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic path material".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1e-3), "{a} != {b}");
}

#[test]
fn g1_surface_seebeck_integrates_cell_breaks_and_closes_kelvin_transfer() {
    for (last, expected_v) in [(6.0, 550e-6), (12.0, 700e-6)] {
        let source = card(
            "temperature",
            SEEBECK,
            SEEBECK_DIMS,
            last,
            false,
            InterpolationPolicy::LinearInside,
        );
        let context = point("temperature", 400.0, 1.5e5);
        let span = ResolvedSeebeckSpan::resolve(
            &source,
            &context,
            "temperature",
            [300.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let left = 2e-6;
        let right = (2.5 + (3.0 + last) / 2.0) * 0.5e-6;
        near(span.edge_seebeck_v_k()[0], left);
        near(span.edge_seebeck_v_k()[1], right);
        for current in [-2.0, 0.0, 2.0] {
            let q = span.transfer(current).unwrap();
            near(q.voltage_drop_v, expected_v);
            near(q.terminal_heat_w[0], -left * 300.0 * current);
            near(q.bulk_heat_w, (left - right) * 400.0 * current);
            near(q.terminal_heat_w[1], right * 500.0 * current);
            near(q.electrical_supply_w, -expected_v * current);
            near(
                q.terminal_heat_w.iter().sum::<f64>() + q.bulk_heat_w + q.electrical_supply_w,
                0.0,
            );
            near(
                q.terminal_heat_w[0] / 300.0 + q.bulk_heat_w / 400.0 + q.terminal_heat_w[1] / 500.0,
                0.0,
            );
        }
        let selected = span.material()[0]
            .property(SEEBECK)
            .unwrap()
            .answer()
            .receipt
            .selected;
        for material in [
            &span.material()[0],
            &span.material()[1],
            span.bulk_material(),
        ] {
            assert_eq!(material.card_identity(), source.content_hash());
            assert_eq!(
                material
                    .property(SEEBECK)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                selected
            );
            assert_eq!(material.query().axes()["pressure"], 1.5e5);
        }
        let reverse = ResolvedSeebeckSpan::resolve(
            &source,
            &context,
            "temperature",
            [500.0, 300.0],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        near(reverse.transfer(2.0).unwrap().voltage_drop_v, -expected_v);
        near(
            reverse.transfer(2.0).unwrap().bulk_heat_w,
            -span.transfer(2.0).unwrap().bulk_heat_w,
        );
        let low_pressure = ResolvedSeebeckSpan::resolve(
            &source,
            &point("temperature", 400.0, 1e5),
            "temperature",
            [300.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        near(low_pressure.transfer(0.0).unwrap().voltage_drop_v, 400e-6);
    }
}

#[test]
fn g1_surface_expansion_drives_actual_restrained_section_force() {
    for (last, pressure, strain) in [
        (6.0, 1e5, 400e-6),
        (6.0, 1.5e5, 550e-6),
        (12.0, 1.5e5, 700e-6),
    ] {
        let source = card(
            "T",
            ALPHA,
            INVERSE_TEMPERATURE_DIMS,
            last,
            false,
            InterpolationPolicy::LinearInside,
        );
        let reference = point("T", 300.0, pressure);
        let current = point("T", 500.0, pressure);
        let path = integrate_isotropic_thermal_expansion(
            &source,
            &reference,
            &current,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        near(path.free_linear_strain(), strain);
        let reverse = integrate_isotropic_thermal_expansion(
            &source,
            &current,
            &reference,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        near(reverse.free_linear_strain(), -strain);
        assert_eq!(path.selected_claim(), reverse.selected_claim());
        let material = resolve_material_state_point(
            &source,
            &current,
            &[
                ScalarPropertyRequirement::try_new(
                    DENSITY_PROPERTY,
                    Density::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
                ScalarPropertyRequirement::try_new(
                    YOUNG_MODULUS_PROPERTY,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
            ],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![ConcentricMaterialLayer {
                outer_radius_m: 0.001,
                material,
            }],
        )
        .unwrap();
        let thermal = ConcentricThermoelasticSection::try_new(section, vec![path]).unwrap();
        let held = thermal.response_at_axial_strain(0.0, 0.01).unwrap();
        near(held.layers()[0].stress_pa(), -200e9 * strain);
        near(
            held.force_n(),
            -std::f64::consts::PI * 1e-6 * 200e9 * strain,
        );
        near(thermal.zero_resultant_axial_strain(), strain);
        let free = thermal.response_at_axial_strain(strain, 0.01).unwrap();
        // Cancellation is measured against the actual restrained force scale.
        assert!(free.force_n().abs() <= 64.0 * f64::EPSILON * held.force_n().abs());
    }
}

#[test]
fn g0_valid_surface_endpoints_do_not_authorize_paths_through_missing_data() {
    let source = card(
        "temperature",
        SEEBECK,
        SEEBECK_DIMS,
        6.0,
        true,
        InterpolationPolicy::LinearInside,
    );
    // All three ordinary queries succeed; only the complete path sees the hole.
    let requirement =
        ScalarPropertyRequirement::try_new(SEEBECK, SEEBECK_DIMS, ScalarAdmissibility::Finite)
            .unwrap();
    for t in [300.0, 350.0, 500.0] {
        resolve_material_state_point(
            &source,
            &point("temperature", t, 1.5e5),
            std::slice::from_ref(&requirement),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    }
    assert!(matches!(
        ResolvedSeebeckSpan::resolve(
            &source,
            &point("temperature", 350.0, 1.5e5),
            "temperature",
            [300.0, 500.0],
            MaterialPropertySelection::SingleClaimOnly
        ),
        Err(ThermoelectricError::Material(
            MaterialStatePointError::Support { .. }
        ))
    ));
    let source = card(
        "T",
        ALPHA,
        INVERSE_TEMPERATURE_DIMS,
        6.0,
        true,
        InterpolationPolicy::LinearInside,
    );
    assert!(matches!(
        integrate_isotropic_thermal_expansion(
            &source,
            &point("T", 300.0, 1.5e5),
            &point("T", 500.0, 1.5e5),
            MaterialPropertySelection::SingleClaimOnly
        ),
        Err(MaterialStatePointError::Support { .. })
    ));
    let discrete = card(
        "T",
        ALPHA,
        INVERSE_TEMPERATURE_DIMS,
        6.0,
        false,
        InterpolationPolicy::TabulatedOnly,
    );
    assert!(matches!(
        integrate_isotropic_thermal_expansion(
            &discrete,
            &point("T", 300.0, 1e5),
            &point("T", 500.0, 1e5),
            MaterialPropertySelection::SingleClaimOnly
        ),
        Err(MaterialStatePointError::Support { .. })
    ));
}

#[test]
fn g0_expansion_requires_temperature_meaning_even_when_source_and_query_agree() {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    let base = card(
        "T",
        ALPHA,
        INVERSE_TEMPERATURE_DIMS,
        6.0,
        false,
        InterpolationPolicy::LinearInside,
    );
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    for (quantity, admitted) in [
        (None, true),
        (Some(QuantitySpec::dimensional(Temperature::DIMS)), true),
        (Some(absolute), true),
        (Some(QuantitySpec::dimensional(Pressure::DIMS)), false),
        (Some(difference), false),
    ] {
        let mut claim = base.claims().claims_for(ALPHA)[0].1.clone();
        let PropertyValue::Surface2D { axis_dims, .. } = &mut claim.value else {
            unreachable!("fixture is a two-coordinate source");
        };
        axis_dims[0] = quantity.map_or(Temperature::DIMS, |q| q.dims());
        claim.validity = ValidityDomain::unconstrained()
            .with("T", 300.0, 500.0)
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                1e5,
                2e5,
            );
        if let Some(quantity) = quantity {
            claim.validity = claim.validity.with_quantity("T", quantity, 300.0, 500.0);
        }
        let mut claims = ClaimSet::new();
        claims.insert_claim(claim).unwrap();
        let source = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
        let query = |t| {
            let point = QueryPoint::new()
                .with("T", t)
                .unwrap()
                .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 1.5e5)
                .unwrap();
            match quantity {
                Some(quantity) => point.with_quantity("T", quantity, t).unwrap(),
                None => point,
            }
        };
        let reference = query(300.0);
        let current = query(500.0);
        // The source is coherent and ordinary typed queries succeed. The
        // physical consumer must additionally own what its T coordinate means.
        let requirement = ScalarPropertyRequirement::try_new(
            ALPHA,
            INVERSE_TEMPERATURE_DIMS,
            ScalarAdmissibility::Finite,
        )
        .unwrap();
        for point in [&reference, &current] {
            resolve_material_state_point(
                &source,
                point,
                std::slice::from_ref(&requirement),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
        }
        let path = integrate_isotropic_thermal_expansion(
            &source,
            &reference,
            &current,
            MaterialPropertySelection::SingleClaimOnly,
        );
        if admitted {
            near(path.unwrap().free_linear_strain(), 550e-6);
        } else {
            assert!(matches!(
                path,
                Err(MaterialStatePointError::InvalidDerived {
                    quantity: "thermal path coordinate must be an absolute temperature",
                })
            ));
        }
    }
}

#[test]
fn g0_untyped_temperature_query_cannot_override_explicit_source_axis_dimensions() {
    for (property, dims) in [(ALPHA, INVERSE_TEMPERATURE_DIMS), (SEEBECK, SEEBECK_DIMS)] {
        for surface in [false, true] {
            let base = card(
                "T",
                property,
                dims,
                6.0,
                false,
                InterpolationPolicy::LinearInside,
            );
            let mut claim = base.claims().claims_for(property)[0].1.clone();
            if surface {
                let PropertyValue::Surface2D { axis_dims, .. } = &mut claim.value else {
                    unreachable!("fixture is a two-coordinate source");
                };
                axis_dims[0] = Pressure::DIMS;
            } else {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Pressure::DIMS,
                    knots: vec![(300.0, 1e-6), (500.0, 3e-6)],
                    dims,
                };
            }
            claim.validity = ValidityDomain::unconstrained()
                .with("T", 300.0, 500.0)
                .with_quantity(
                    "pressure",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    1e5,
                    2e5,
                );
            let mut claims = ClaimSet::new();
            claims.insert_claim(claim).unwrap();
            let source = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
            let query = |t| {
                QueryPoint::new()
                    .with("T", t)
                    .unwrap()
                    .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 1.5e5)
                    .unwrap()
            };
            let requirement =
                ScalarPropertyRequirement::try_new(property, dims, ScalarAdmissibility::Finite)
                    .unwrap();
            // Legacy untyped point lookup is valid. These source coordinates
            // are nevertheless pressure, and cannot become kelvin by naming T.
            for t in [300.0, 400.0, 500.0] {
                resolve_material_state_point(
                    &source,
                    &query(t),
                    std::slice::from_ref(&requirement),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
            }
            if property == ALPHA {
                assert!(matches!(
                    integrate_isotropic_thermal_expansion(
                        &source,
                        &query(300.0),
                        &query(500.0),
                        MaterialPropertySelection::SingleClaimOnly
                    ),
                    Err(MaterialStatePointError::InvalidDerived {
                        quantity: "temperature path source coordinate must have temperature dimensions",
                    })
                ));
            } else {
                assert!(matches!(
                    ResolvedSeebeckSpan::resolve(
                        &source,
                        &query(400.0),
                        "T",
                        [300.0, 500.0],
                        MaterialPropertySelection::SingleClaimOnly
                    ),
                    Err(ThermoelectricError::Material(
                        MaterialStatePointError::InvalidDerived {
                            quantity: "temperature path source coordinate must have temperature dimensions",
                        }
                    ))
                ));
            }
        }
    }
}
