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
    MaterialOpticalBinding, MaterialOpticalError, MaterialOpticalRefreshError,
    MaterialThermalBinding, MaterialThermalError,
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
        // Retained specimen context, not an axis of these optical source laws.
        .with_quantity("load_context", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
        .with_quantity(
            "liquid_fraction",
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap()
}

fn condition_card(axis: &str) -> MaterialCard {
    condition_card_with_spectrum(axis, None)
}

fn condition_card_with_spectrum(axis: &str, absorption_per_m: Option<[f64; 9]>) -> MaterialCard {
    use fs_material::state_point::{
        VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
        VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
        VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY, VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
        VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES,
    };
    let mut claims = ClaimSet::new();
    let mut rows = vec![
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
    ];
    if let Some(coefficients) = absorption_per_m {
        rows.extend(
            VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES
                .into_iter()
                .zip(coefficients)
                .map(|(name, value)| (name, Dims([-1, 0, 0, 0, 0, 0]), value)),
        );
    } else {
        rows.push((
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            Dims([1, 0, 0, 0, 0, 0]),
            1.0,
        ));
        rows.extend(
            VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES
                .into_iter()
                .map(|name| (name, Dims::NONE, 1.0)),
        );
    }
    for (name, dims, value) in rows {
        let spectral = VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&name);
        let varying = spectral || name == SURFACE_EMISSIVITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: if spectral {
                    PropertyKey::with_quantity(name, VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY)
                } else {
                    PropertyKey::new(name, dims)
                },
                value: if varying {
                    PropertyValue::Curve {
                        abscissa: axis.into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: if spectral {
                            vec![(1200.0, value), (1800.0, 2.0 * value)]
                        } else {
                            vec![(1200.0, 0.2), (1800.0, 0.6)]
                        },
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        axis,
                        quantity(QuantityKind::AbsoluteTemperature),
                        1200.0,
                        1800.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        100_000.0,
                        200_000.0,
                    )
                    .with_quantity(
                        "liquid_fraction",
                        QuantitySpec::dimensional(Dims::NONE),
                        0.0,
                        1.0,
                    ),
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
            .with_quantity("load_context", quantity(QuantityKind::Stress), 100_000.0)
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
        .with_quantity("load_context", quantity(QuantityKind::Stress), 100_000.0)
        .unwrap();
    let accepted = state_at(&card, &other);
    assert_eq!(accepted.query_point(), optical.resolved().query_point());
    assert!(matches!(
        MaterialOpticalBinding::try_dielectric(&accepted, &optical, None, ContentHash([19; 32]),),
        Err(MaterialOpticalError::StateMismatch)
    ));
}

fn sampled_dielectric_card() -> MaterialCard {
    use fs_material::state_point::{
        VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
        VISIBLE_DIELECTRIC_IOR_PROPERTIES,
    };
    let base = card(0.2);
    let mut claims = base.claims().clone();
    for index in 0..9 {
        let eta = 1.6 - index as f64 * 0.01;
        let absorption = if index == 3 { 50.0 } else { 2.0 };
        for (name, quantity, cold, hot) in [
            (
                VISIBLE_DIELECTRIC_IOR_PROPERTIES[index],
                QuantitySpec::dimensional(Dims::NONE),
                eta,
                eta + 0.1,
            ),
            (
                VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
                absorption,
                2.0 * absorption,
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, quantity),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(1200.0, cold), (1800.0, hot)],
                        dims: quantity.dims(),
                    },
                    validity: ValidityDomain::unconstrained().with("T", 1200.0, 1800.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: "synthetic temperature-dependent dielectric spectrum".into(),
                        license: "internal-test-use".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::supersede(&base, claims, Vec::new()).unwrap()
}

#[test]
fn g1_tabulated_material_drives_chromatic_snell_and_spectral_beer_lambert() {
    use fs_material::state_point::resolve_visible_optical_state_point;
    use fs_render::dielectric::{DielectricEvent, sample_smooth_dielectric};
    let source = sampled_dielectric_card();
    let mut identities = Vec::new();
    for (temperature, delta_n, absorption_factor) in [(1200.0, 0.0, 1.0), (1800.0, 0.1, 2.0)] {
        let point = QueryPoint::new().with("T", temperature).unwrap();
        let optical = resolve_visible_optical_state_point(
            &source,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let accepted = state_at(&source, &point);
        let binding =
            MaterialOpticalBinding::try_visible(&accepted, &optical, 0.1, ContentHash([74; 32]))
                .unwrap();
        identities.push(binding.identity());
        let Material::Dielectric { glass, .. } = binding.material() else {
            panic!("sampled material must reach dielectric transport");
        };
        for (lambda, index, coefficient) in [
            (530.0, 1.57, 50.0),
            (555.0, 1.565, 26.0),
            (680.0, 1.54, 2.0),
        ] {
            let n = index + delta_n;
            assert!((glass.ior().eval(lambda).unwrap() - n).abs() < 1.0e-14);
            let refracted = sample_smooth_dielectric(
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.6, 0.0, 0.8),
                1.0,
                glass.ior().eval(lambda).unwrap(),
                0.99,
            )
            .unwrap();
            assert_eq!(refracted.event, DielectricEvent::Transmission);
            assert!((refracted.direction.x + 0.6 / n).abs() < 1.0e-14);
            let expected = (-coefficient * absorption_factor * 0.02_f64).exp();
            assert!(
                (glass.absorption().transmittance(lambda, 0.02).unwrap() - expected).abs()
                    < 1.0e-13
            );
        }
        assert_eq!(binding.accepted_state_identity(), accepted.identity());
    }
    assert_ne!(identities[0], identities[1]);
}

#[test]
fn g0_sampled_dielectric_binding_refuses_stale_state_and_missing_surface() {
    use fs_material::state_point::resolve_visible_tabulated_dielectric_state_point;
    let source = sampled_dielectric_card();
    let point = QueryPoint::new().with("T", 1500.0).unwrap();
    let optical = resolve_visible_tabulated_dielectric_state_point(
        &source,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        MaterialOpticalBinding::try_tabulated_dielectric(
            &state(&source, 1501.0),
            &optical,
            None,
            ContentHash([74; 32]),
        )
        .unwrap_err(),
        MaterialOpticalError::StateMismatch,
    );
    assert_eq!(
        MaterialOpticalBinding::try_tabulated_dielectric(
            &state_at(&source, &point),
            &optical,
            None,
            ContentHash([0; 32]),
        )
        .unwrap_err(),
        MaterialOpticalError::MissingSurfaceIdentity,
    );
    let binding = MaterialOpticalBinding::try_tabulated_dielectric(
        &state_at(&source, &point),
        &optical,
        None,
        ContentHash([74; 32]),
    )
    .unwrap();
    let Material::Dielectric { surface, .. } = binding.material() else {
        panic!("sampled dielectric binding");
    };
    assert_eq!(surface, DielectricSurface::SMOOTH);
}

#[test]
fn g3_sampled_dielectric_refresh_preserves_family_sources_and_finish() {
    use fs_material::state_point::{VisibleOpticalStatePoint, resolve_visible_optical_state_point};
    let source = sampled_dielectric_card();
    let cold = state(&source, 1200.0);
    let warm = state(&source, 1800.0);
    let optical = resolve_visible_optical_state_point(
        &source,
        &cold.query(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let VisibleOpticalStatePoint::TabulatedDielectric(sampled) = &optical else {
        panic!("sampled fixture family");
    };
    for alpha in [None, Some(0.2)] {
        let binding = MaterialOpticalBinding::try_tabulated_dielectric(
            &cold,
            sampled,
            alpha,
            ContentHash([74; 32]),
        )
        .unwrap();
        let (refreshed, refreshed_optical) =
            binding.at_material_state(&optical, &source, &warm).unwrap();
        assert!(matches!(
            &refreshed_optical,
            VisibleOpticalStatePoint::TabulatedDielectric(_)
        ));
        assert_eq!(refreshed.accepted_state_identity(), warm.identity());
        assert_ne!(refreshed.identity(), binding.identity());
        let Material::Dielectric { glass, surface } = refreshed.material() else {
            panic!("refreshed sampled dielectric material");
        };
        assert_eq!(surface.roughness_alpha(), alpha);
        assert!((glass.ior().eval(530.0).unwrap() - 1.67).abs() < 1.0e-14);
        assert!(
            (glass.absorption().transmittance(530.0, 0.02).unwrap() - (-2.0_f64).exp()).abs()
                < 1.0e-13
        );
        for (before, after) in optical
            .resolved()
            .properties()
            .iter()
            .zip(refreshed_optical.resolved().properties())
        {
            assert_eq!(
                before.answer().receipt.selected,
                after.answer().receipt.selected
            );
        }
        // The accepted mechanical state can lie beyond optical source support.
        // Refresh must refuse instead of keeping stale values or changing law.
        assert!(
            binding
                .at_material_state(&optical, &source, &state(&source, 1801.0))
                .is_err()
        );
        assert_eq!(binding.accepted_state_identity(), cold.identity());
    }
}

#[test]
fn g1_spectral_material_absorption_reaches_the_same_accepted_dielectric() {
    use fs_material::state_point::{
        VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM, VisibleOpticalStatePoint,
        resolve_visible_dielectric_state_point,
    };
    use fs_render::dielectric::BeerLambertParameters;

    let coefficients = [0.0, 2.0, 4.0, 8.0, 16.0, 12.0, 8.0, 4.0, 2.0];
    let card = condition_card_with_spectrum("Ts", Some(coefficients));
    let initial_point = condition("Ts", 1200.0);
    let initial_state = state_at(&card, &initial_point);
    let initial_optical = resolve_visible_dielectric_state_point(
        &card,
        &initial_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut binding = MaterialOpticalBinding::try_dielectric(
        &initial_state,
        &initial_optical,
        None,
        ContentHash([19; 32]),
    )
    .unwrap();
    let mut optical = VisibleOpticalStatePoint::Dielectric(initial_optical);
    let mut identities = Vec::new();
    for (temperature_k, factor) in [(1200.0, 1.0), (1800.0, 2.0)] {
        let point = condition("Ts", temperature_k);
        let accepted = initial_state.at_query_point(&card, &point).unwrap();
        (binding, optical) = binding
            .at_material_state(&optical, &card, &accepted)
            .expect("original spectral claims follow the accepted thermal endpoint");
        let Material::Dielectric { glass, surface } = binding.material() else {
            panic!("spectral glass must reach the native dielectric material");
        };
        assert!(
            surface.is_delta(),
            "refresh must preserve the smooth finish"
        );
        let BeerLambertParameters::Spectral { spectrum } = glass.absorption().parameters() else {
            panic!("material spectrum was reduced to RGB");
        };
        assert_eq!(spectrum.knots().len(), 9);
        for ((&wavelength_nm, &coefficient), &knot) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM
            .iter()
            .zip(&coefficients)
            .zip(spectrum.knots())
        {
            assert_eq!(knot, (wavelength_nm, factor * coefficient));
            let expected = (-factor * coefficient * 0.1).exp();
            let actual = glass
                .absorption()
                .transmittance(wavelength_nm, 0.1)
                .unwrap();
            assert!((actual - expected).abs() < 3.0e-15);
        }
        assert_eq!(binding.accepted_state_identity(), accepted.identity());
        identities.push(binding.identity());
    }
    assert_ne!(
        identities[0], identities[1],
        "changed accepted absorption must change identity"
    );
}

#[test]
fn g3_optical_continuation_changes_native_conductor_response_and_rejects_substitution() {
    use fs_material::state_point::{
        VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, VISIBLE_COMPLEX_IOR_K_PROPERTIES,
        VisibleOpticalStatePoint, resolve_visible_conductor_state_point,
    };
    let base = condition_card("Ts");
    let mut claims = base.claims().clone();
    let template = base.claims().claims_ordered().next().unwrap().1.clone();
    for (names, values) in [
        (VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, [1.5, 2.0]),
        (VISIBLE_COMPLEX_IOR_K_PROPERTIES, [2.0, 3.0]),
    ] {
        for name in names {
            let mut source = template.clone();
            source.key = PropertyKey::new(name, Dims::NONE);
            source.value = PropertyValue::Curve {
                abscissa: "Ts".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(1200.0, values[0]), (1800.0, values[1])],
                dims: Dims::NONE,
            };
            source.validity = ValidityDomain::unconstrained().with_quantity(
                "Ts",
                quantity(QuantityKind::AbsoluteTemperature),
                1200.0,
                1600.0,
            );
            source.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(source).unwrap();
        }
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, Vec::new()).unwrap();
    let cold = state_at(&card, &condition("Ts", 1200.0));
    let original = VisibleOpticalStatePoint::Conductor(
        resolve_visible_conductor_state_point(
            &card,
            &cold.query(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    );
    let binding =
        MaterialOpticalBinding::try_visible(&cold, &original, 0.15, ContentHash([29; 32])).unwrap();
    // Refresh receipts honestly record pinned selection, so returning to the
    // initial endpoint matches its pinned evaluation, not its discovery policy.
    let pinned_at_cold = binding.at_material_state(&original, &card, &cold).unwrap().0;
    assert_eq!(pinned_at_cold.accepted_state_identity(), cold.identity());
    let hot = cold
        .at_query_point(&card, &condition("Ts", 1600.0))
        .unwrap();
    let (heated, optical) = binding.at_material_state(&original, &card, &hot).unwrap();
    assert_eq!(heated.accepted_state_identity(), hot.identity());
    let mut reflectances = Vec::new();
    for (bound, eta, k) in [(binding, 1.5_f64, 2.0_f64), (heated, 11.0 / 6.0, 8.0 / 3.0)] {
        let Material::Conductor { optics, surface } = bound.material() else {
            panic!("an available dielectric family cannot replace the original conductor");
        };
        assert_eq!(surface.roughness_alpha(), 0.15);
        let expected = ((eta - 1.0).powi(2) + k * k) / ((eta + 1.0).powi(2) + k * k);
        let actual = optics.fresnel(580.0, 1.0, 1.0).unwrap();
        assert!((actual - expected).abs() < 2e-15);
        reflectances.push(actual);
    }
    assert_ne!(reflectances[0], reflectances[1]);
    with_cx(|cx| {
        let pixel = |bound: MaterialOpticalBinding| {
            let mut s = scene(bound.material());
            s.environment =
                Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap());
            render(&s, cx, &settings(DirectStrategy::BsdfOnly))
                .unwrap()
                .xyz[0][1]
        };
        assert!(
            pixel(heated) > pixel(binding),
            "the same native scene must observe changed source optics"
        );
    });
    assert!(matches!(
        heated.at_material_state(&original, &card, &hot),
        Err(MaterialOpticalRefreshError::SourceMismatch)
    ));
    assert!(matches!(
        binding.at_material_state(&optical, &card, &hot),
        Err(MaterialOpticalRefreshError::SourceMismatch)
    ));
    let foreign = MaterialCard::supersede(&card, card.claims().clone(), Vec::new()).unwrap();
    assert!(matches!(
        binding.at_material_state(&original, &foreign, &hot),
        Err(MaterialOpticalRefreshError::Source(_))
    ));
    for point in [
        condition("Ts", 1800.0),
        condition("Ts", 1600.0)
            .with_quantity("load_context", quantity(QuantityKind::Stress), 100_000.0)
            .unwrap(),
    ] {
        let accepted = state_at(&card, &point);
        assert!(matches!(
            binding.at_material_state(&original, &card, &accepted),
            Err(MaterialOpticalRefreshError::Source(_))
        ));
    }
    assert_eq!(
        heated.at_material_state(&optical, &card, &cold).unwrap().0,
        pinned_at_cold
    );
}

#[test]
fn g0_sellmeier_continuation_preserves_spectral_state_and_readmits_the_transparent_band() {
    use fs_material::state_point::{
        VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
        VisibleOpticalStatePoint, resolve_visible_sellmeier_state_point,
    };
    use fs_render::dielectric::{BeerLambertParameters, DielectricError, DielectricIor};

    let base = condition_card_with_spectrum("Ts", Some([2.0; 9]));
    let mut claims = base.claims().clone();
    let template = base.claims().claims_ordered().next().unwrap().1.clone();
    for (names, dims, cold, hot) in [
        (
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES,
            Dims::NONE,
            [1.0, 0.2, 1.0],
            [1.1, 0.2, 1.0],
        ),
        (
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
            Dims([2, 0, 0, 0, 0, 0]),
            [0.0, 2.0e-14, 1.0e-10],
            [2.0e-13, 2.0e-14, 1.0e-10],
        ),
    ] {
        for (index, name) in names.into_iter().enumerate() {
            let mut source = template.clone();
            source.key = PropertyKey::new(name, dims);
            source.value = PropertyValue::Curve {
                abscissa: "Ts".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(1200.0, cold[index]), (1800.0, hot[index])],
                dims,
            };
            source.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(source).unwrap();
        }
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, Vec::new()).unwrap();
    let initial = state_at(&card, &condition("Ts", 1200.0));
    let optical = VisibleOpticalStatePoint::DielectricSellmeier(
        resolve_visible_sellmeier_state_point(
            &card,
            &initial.query(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    );
    let binding =
        MaterialOpticalBinding::try_visible(&initial, &optical, 0.12, ContentHash([31; 32]))
            .unwrap();
    let hot = initial
        .at_query_point(&card, &condition("Ts", 1500.0))
        .unwrap();
    let (heated, successor) = binding.at_material_state(&optical, &card, &hot).unwrap();
    let Material::Dielectric { glass, surface } = heated.material() else {
        panic!("Sellmeier continuation must retain dielectric transport");
    };
    let DielectricIor::Sellmeier(law) = glass.ior() else {
        panic!("available Cauchy properties cannot replace the selected Sellmeier law");
    };
    assert_eq!(surface.roughness_alpha(), Some(0.12));
    assert_eq!(heated.accepted_state_identity(), hot.identity());
    let square = (580.0e-9_f64).powi(2);
    let expected_squared = 1.0
        + [1.05, 0.2, 1.0]
            .into_iter()
            .zip([1.0e-13, 2.0e-14, 1.0e-10])
            .map(|(strength, resonance)| strength * square / (square - resonance))
            .sum::<f64>();
    assert!((law.eval(580.0).unwrap().powi(2) - expected_squared).abs() < 3.0e-14);
    assert!(matches!(
        glass.absorption().parameters(),
        BeerLambertParameters::Spectral { .. }
    ));
    assert!(
        (glass.absorption().transmittance(580.0, 0.1).unwrap() - (-0.3_f64).exp()).abs() < 3e-15
    );
    for (before, after) in optical
        .resolved()
        .properties()
        .iter()
        .zip(successor.resolved().properties())
    {
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
    }
    let unsupported = initial
        .at_query_point(&card, &condition("Ts", 1800.0))
        .unwrap();
    assert!(
        successor
            .at_query_point(&card, &unsupported.query())
            .is_ok(),
        "material coefficient support alone cannot admit the optical band"
    );
    assert!(matches!(
        heated.at_material_state(&successor, &card, &unsupported),
        Err(MaterialOpticalRefreshError::Binding(
            MaterialOpticalError::Dielectric(DielectricError::InvalidIor)
        ))
    ));
    assert_eq!(
        heated
            .at_material_state(&successor, &card, &initial)
            .unwrap()
            .0,
        binding.at_material_state(&optical, &card, &initial).unwrap().0
    );
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
