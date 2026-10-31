//! Original sources -> accepted local heating/deformation -> native scene pixels.

use super::*;
use fs_material::state_point::{
    VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, VISIBLE_COMPLEX_IOR_K_PROPERTIES, VisibleOpticalStatePoint,
    resolve_visible_conductor_state_point,
};
use fs_qty::QuantitySpec;
use fs_render::lighting::EnvironmentMap;
use fs_render::material_geometry::regional::{
    MaterialRegionContinuation, MaterialRegionContinuationError, MaterialRegionSource,
    continue_thermal_body_primitives,
};
use fs_render::material_state::{MaterialOpticalBinding, MaterialOpticalRefreshError};

fn surface_card(name: &str) -> MaterialCard {
    let original = card(name, 2000.0, 1.0e-5);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        if claim.key.name() != SURFACE_EMISSIVITY_PROPERTY {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    let curve = |name: &str, knots: Vec<(f64, f64)>| {
        let low = knots.first().unwrap().0;
        let high = knots.last().unwrap().0;
        PropertyClaim {
            key: PropertyKey::new(name, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", low, high)
                .with("pressure", 100_000.0, 100_000.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic temperature-dependent surface; not measured data".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        }
    };
    for (index, name) in VISIBLE_COMPLEX_IOR_ETA_PROPERTIES.into_iter().enumerate() {
        let limit = if index == 0 { 1000.0 } else { 2000.0 };
        claims
            .insert_claim(curve(name, vec![(300.0, 1.5), (900.0, 2.0), (limit, 2.0)]))
            .unwrap();
    }
    for name in VISIBLE_COMPLEX_IOR_K_PROPERTIES {
        claims
            .insert_claim(curve(name, vec![(300.0, 2.0), (900.0, 5.0), (2000.0, 5.0)]))
            .unwrap();
    }
    // Fresh discovery is ambiguous at the hot endpoint and would substitute
    // these alternate claims above 1000 K. Continuation must retain its pins.
    claims
        .insert_claim(curve(
            VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[0],
            vec![(700.0, 7.0), (2000.0, 7.0)],
        ))
        .unwrap();
    claims
        .insert_claim(curve(
            SURFACE_EMISSIVITY_PROPERTY,
            vec![(300.0, 0.2), (1000.0, 0.6)],
        ))
        .unwrap();
    claims
        .insert_claim(curve(
            SURFACE_EMISSIVITY_PROPERTY,
            vec![(700.0, 0.9), (2000.0, 0.9)],
        ))
        .unwrap();
    MaterialCard::assemble(original.id().clone(), claims, Vec::new()).unwrap()
}

fn typed_reference(cards: &[MaterialCard]) -> ResolvedBody {
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            100_000.0,
        )
        .unwrap();
    let bindings: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, material)| RegionBinding {
            key: KEYS[i],
            material,
            point: &point,
            geometry: VolumeMoments::tetrahedron(TETS[i].map(|node| {
                let [x, y, z] = NODES[node];
                MassVec3::new(x, y, z)
            }))
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    resolve_body(&bindings, 2, || false).unwrap()
}

fn initial_sources<'a>(
    cards: &'a [MaterialCard],
    body: &ThermalBody<'_>,
) -> Vec<MaterialRegionContinuation<'a>> {
    // Reverse canonical order intentionally; neither caller nor solver element
    // order is an appearance assignment.
    KEYS.iter()
        .enumerate()
        .map(|(input, &key)| {
            let index = body
                .body()
                .reference()
                .regions()
                .iter()
                .position(|region| region.key() == key)
                .unwrap();
            let accepted = body.sources()[index].stiffness();
            let point = accepted.query();
            let source = if key == "upper" {
                let optics = resolve_visible_conductor_state_point(
                    &cards[input],
                    &point,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
                MaterialRegionSource::Optical {
                    binding: MaterialOpticalBinding::try_conductor(
                        accepted,
                        &optics,
                        0.05,
                        ContentHash([83; 32]),
                    )
                    .unwrap(),
                    source: VisibleOpticalStatePoint::Conductor(optics),
                }
            } else {
                let (emissivity, _) = SurfaceEmissivity::from_condition(
                    key,
                    ThermalMaterialPoint {
                        material: &cards[input],
                        context: &point,
                        temperature_axis: "T",
                        liquid_fraction_axis: None,
                    },
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
                MaterialRegionSource::Thermal {
                    binding: MaterialThermalBinding::try_gray_diffuse(
                        accepted,
                        &emissivity,
                        ContentHash([97; 32]),
                    )
                    .unwrap(),
                    temperature_axis: "T".into(),
                }
            };
            MaterialRegionContinuation {
                key,
                card: &cards[input],
                source,
            }
        })
        .collect()
}

fn observed_pixels(scene: &mut Scene, cx: &Cx<'_>) -> [f64; 2] {
    let settings = Settings {
        width: 1,
        height: 1,
        spp: 64,
        max_depth: 2,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::BsdfOnly,
        seed: 41,
    };
    scene.environment =
        Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap());
    scene.camera.eye = Point3::new(0.1, 0.1, 2.0);
    scene.camera.forward = Vec3::new(0.0, 0.0, -1.0);
    let reflected = render(scene, cx, &settings).unwrap().xyz[0][1];
    scene.environment = None;
    scene.camera.eye = Point3::new(0.1, 0.1, -2.0);
    scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
    let emitted = render(scene, cx, &settings).unwrap().xyz[0][1];
    [reflected, emitted]
}

#[test]
fn accepted_local_heating_refreshes_original_sources_and_native_pixels_together() {
    let cards = [surface_card("upper"), surface_card("lower")];
    let reference = typed_reference(&cards);
    let mass = TetrahedralMass::bind(&reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    with_cx(|cx| {
        let cold = mass
            .solve_thermal(
                &cells(&cards, &[300.0, 300.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let hot = mass
            .solve_thermal(
                &cells(&cards, &[900.0, 800.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let original = initial_sources(&cards, &cold);
        let (cold_primitives, cold_sources) =
            continue_thermal_body_primitives(&cold, &original, BUDGET, 8, cx).unwrap();
        let cold_pixels = observed_pixels(&mut scene(cold_primitives), cx);
        let (hot_primitives, hot_sources) =
            continue_thermal_body_primitives(&hot, &cold_sources, BUDGET, 8, cx).unwrap();
        assert_ne!(hot.body().nodes_m(), cold.body().nodes_m());
        assert_eq!(
            hot_sources
                .iter()
                .map(|source| source.key)
                .collect::<Vec<_>>(),
            ["lower", "upper"]
        );
        for (index, (primitive, successor)) in hot_primitives.iter().zip(&hot_sources).enumerate() {
            let accepted = hot.sources()[index].stiffness();
            assert_eq!(accepted.state_coordinate("pressure"), Some(100_000.0));
            assert_eq!(accepted.axis_quantities().len(), 2);
            let Shape::Mesh(mesh) = &primitive.shape else {
                panic!("accepted exterior mesh required")
            };
            let nodes = TETS[hot.body().region_elements()[index]];
            assert_eq!(mesh.triangles.len(), 3);
            assert_eq!(mesh.vertices, nodes.map(|node| hot.body().nodes_m()[node]));
            assert!(
                primitive.emission.is_none(),
                "thermal source must not be counted twice"
            );
            match (
                &cold_sources[index].source,
                &successor.source,
                primitive.material,
            ) {
                (
                    MaterialRegionSource::Optical {
                        binding: old,
                        source: before,
                    },
                    MaterialRegionSource::Optical { binding, source },
                    Material::Conductor { optics, surface },
                ) => {
                    assert_eq!(binding.accepted_state_identity(), accepted.identity());
                    assert_eq!(source.resolved().query_point(), accepted.query_point());
                    assert_eq!(
                        source.resolved().axis_quantities(),
                        accepted.axis_quantities()
                    );
                    let pins = |state: &VisibleOpticalStatePoint| {
                        state
                            .resolved()
                            .properties()
                            .iter()
                            .map(|property| property.answer().receipt.selected)
                            .collect::<Vec<_>>()
                    };
                    assert_eq!(pins(source), pins(before));
                    let Material::Conductor {
                        surface: original_finish,
                        ..
                    } = old.material()
                    else {
                        unreachable!()
                    };
                    assert_eq!(surface, original_finish);
                    // Independent normal-incidence complex-index Fresnel oracle.
                    near(
                        optics.fresnel(550.0, 1.0, 1.0).unwrap(),
                        (1.0 + 25.0) / (9.0 + 25.0),
                    );
                    assert!(
                        resolve_visible_conductor_state_point(
                            successor.card,
                            &accepted.query(),
                            MaterialPropertySelection::SingleClaimOnly
                        )
                        .is_err()
                    );
                }
                (
                    MaterialRegionSource::Thermal { binding: old, .. },
                    MaterialRegionSource::Thermal {
                        binding,
                        temperature_axis,
                    },
                    Material::GrayDiffuseThermal { surface },
                ) => {
                    assert_eq!(binding.accepted_state_identity(), accepted.identity());
                    assert_eq!(
                        binding.emissivity().receipt().query_point.as_slice(),
                        accepted.query_point()
                    );
                    assert_eq!(
                        &binding.emissivity().receipt().axis_quantities,
                        accepted.axis_quantities()
                    );
                    assert_eq!(
                        binding.emissivity().receipt().selected,
                        old.emissivity().receipt().selected
                    );
                    assert_eq!(
                        binding.surface_state_identity(),
                        old.surface_state_identity()
                    );
                    assert_eq!(temperature_axis, "T");
                    assert_eq!(surface.temperature_k(), 800.0);
                    near(surface.emissivity(), 0.2 + 0.4 * 500.0 / 700.0);
                    assert!(
                        SurfaceEmissivity::from_condition(
                            successor.key,
                            ThermalMaterialPoint {
                                material: successor.card,
                                context: &accepted.query(),
                                temperature_axis,
                                liquid_fraction_axis: None,
                            },
                            MaterialPropertySelection::SingleClaimOnly
                        )
                        .is_err()
                    );
                }
                _ => panic!("a region's original appearance model changed"),
            }
        }
        let hot_pixels = observed_pixels(&mut scene(hot_primitives), cx);
        assert!(hot_pixels[0].is_finite() && hot_pixels[0] > cold_pixels[0] * 1.1);
        assert!(hot_pixels[1].is_finite() && hot_pixels[1] > 0.0);
        assert!(cold_pixels[1] >= 0.0 && cold_pixels[1] < hot_pixels[1] * 1.0e-6);
        // Continue from returned sources back to the same accepted cold state.
        // Pinning policy is already canonical after the first continuation.
        let (returned, returned_sources) =
            continue_thermal_body_primitives(&cold, &hot_sources, BUDGET, 8, cx).unwrap();
        assert_eq!(observed_pixels(&mut scene(returned), cx), cold_pixels);
        assert_eq!(
            returned_sources
                .iter()
                .map(|s| &s.source)
                .collect::<Vec<_>>(),
            cold_sources.iter().map(|s| &s.source).collect::<Vec<_>>()
        );
    });
}

#[test]
fn continuation_refuses_late_source_substitution_without_publishing_partial_regions() {
    let cards = [surface_card("upper"), surface_card("lower")];
    let reference = typed_reference(&cards);
    let mass = TetrahedralMass::bind(&reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    with_cx(|cx| {
        let cold = mass
            .solve_thermal(
                &cells(&cards, &[300.0, 300.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let original = initial_sources(&cards, &cold);
        let (primitives, sources) =
            continue_thermal_body_primitives(&cold, &original, BUDGET, 8, cx).unwrap();
        let retained: Vec<_> = sources.iter().map(|s| s.source.clone()).collect();
        let mut published = scene(primitives);
        let before = observed_pixels(&mut published, cx);

        // Lower's thermal source succeeds first. Upper's original eta claim
        // then expires even though fresh discovery could select a replacement.
        let optical_outside = mass
            .solve_thermal(
                &cells(&cards, &[1200.0, 900.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        assert!(
            resolve_visible_conductor_state_point(
                &cards[0],
                &optical_outside.sources()[1].stiffness().query(),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_ok()
        );
        assert!(
            matches!(continue_thermal_body_primitives(&optical_outside, &sources, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Optical { region, .. }) if region == "upper")
        );

        let thermal_outside = mass
            .solve_thermal(
                &cells(&cards, &[900.0, 1200.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let point = thermal_outside.sources()[0].stiffness().query();
        assert!(
            SurfaceEmissivity::from_condition(
                "lower",
                ThermalMaterialPoint {
                    material: &cards[1],
                    context: &point,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_ok()
        );
        assert!(
            matches!(continue_thermal_body_primitives(&thermal_outside, &sources, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::ThermalSource { region, .. }) if region == "lower")
        );

        let mut foreign = sources.clone();
        let MaterialRegionSource::Optical { source, .. } = &mut foreign[1].source else {
            unreachable!()
        };
        *source = VisibleOpticalStatePoint::Conductor(
            resolve_visible_conductor_state_point(
                &cards[1],
                &cold.sources()[0].stiffness().query(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        );
        assert!(matches!(
            continue_thermal_body_primitives(&cold, &foreign, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Optical {
                error: MaterialOpticalRefreshError::SourceMismatch,
                ..
            })
        ));
        let mut wrong_axis = sources.clone();
        let MaterialRegionSource::Thermal {
            temperature_axis, ..
        } = &mut wrong_axis[0].source
        else {
            unreachable!()
        };
        *temperature_axis = "pressure".into();
        assert!(matches!(
            continue_thermal_body_primitives(&cold, &wrong_axis, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Binding(
                MaterialRegionBindingError::StateMismatch { .. }
            ))
        ));
        let mut foreign_card = sources.clone();
        foreign_card[0].card = &cards[0];
        assert!(matches!(
            continue_thermal_body_primitives(&cold, &foreign_card, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Binding(
                MaterialRegionBindingError::StateMismatch { .. }
            ))
        ));
        // Equal numeric coordinates do not authorize dropping their quantities.
        let untyped_reference = super::super::reference(&cards);
        let untyped_mass =
            TetrahedralMass::bind(&untyped_reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
        let untyped = untyped_mass
            .solve_thermal(
                &cells(&cards, &[300.0, 300.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            continue_thermal_body_primitives(&untyped, &sources, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Binding(
                MaterialRegionBindingError::StateMismatch { .. }
            ))
        ));
        let duplicate = [sources[0].clone(), sources[0].clone()];
        assert!(matches!(
            continue_thermal_body_primitives(&cold, &duplicate, BUDGET, 8, cx),
            Err(MaterialRegionContinuationError::Binding(
                MaterialRegionBindingError::Bindings(_)
            ))
        ));
        assert!(continue_thermal_body_primitives(&cold, &sources[..1], BUDGET, 8, cx).is_err());
        let gate = CancelGate::new_clock_free();
        gate.request();
        assert!(
            with_gate(&gate, |cancelled| continue_thermal_body_primitives(
                &cold, &sources, BUDGET, 8, cancelled
            ))
            .is_err()
        );
        assert_eq!(
            sources.iter().map(|s| &s.source).collect::<Vec<_>>(),
            retained.iter().collect::<Vec<_>>()
        );
        assert_eq!(observed_pixels(&mut published, cx), before);
    });
}
