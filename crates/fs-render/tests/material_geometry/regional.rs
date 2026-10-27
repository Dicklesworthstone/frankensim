use super::*;
use fs_mbd_material::deformation::thermoelastic::ThermalBody;
use fs_render::material_geometry::regional::{
    MaterialRegionAppearance, MaterialRegionBinding, MaterialRegionBindingError,
    deformed_body_regions, thermal_body_primitives,
};
use fs_render::tracer::Material;
use std::collections::BTreeSet;

fn appearance(
    body: &ThermalBody<'_>,
    cards: &[MaterialCard],
    temperatures: &[f64],
) -> Vec<(String, MaterialThermalBinding)> {
    body.body()
        .reference()
        .regions()
        .iter()
        .zip(body.sources())
        .map(|(region, source)| {
            let index = KEYS.iter().position(|key| *key == region.key()).unwrap();
            let point = source.stiffness().query();
            let (emissivity, _) = SurfaceEmissivity::from_condition(
                region.key(),
                ThermalMaterialPoint {
                    material: &cards[index],
                    context: &point,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            assert_eq!(emissivity.temperature_k(), temperatures[index]);
            let binding = MaterialThermalBinding::try_gray_diffuse(
                source.stiffness(),
                &emissivity,
                ContentHash([73; 32]),
            )
            .unwrap();
            assert_eq!(
                binding.accepted_state_identity(),
                source.stiffness().identity()
            );
            (region.key().to_owned(), binding)
        })
        .collect()
}

// Deliberately supply reverse canonical order: input position cannot decide
// which cell's material or temperature appears on a physical surface.
fn bindings(appearance: &[(String, MaterialThermalBinding)]) -> Vec<MaterialRegionBinding<'_>> {
    appearance
        .iter()
        .rev()
        .map(|(key, binding)| MaterialRegionBinding {
            key,
            appearance: MaterialRegionAppearance::Thermal(binding),
        })
        .collect()
}

fn pixels(scene: &mut Scene, cx: &Cx<'_>) -> [f64; 2] {
    let settings = Settings {
        width: 1,
        height: 1,
        spp: 64,
        max_depth: 1,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::BsdfOnly,
        seed: 41,
    };
    [2.0_f64, -2.0].map(|z| {
        scene.camera.eye = Point3::new(0.1, 0.1, z);
        scene.camera.forward = Vec3::new(0.0, 0.0, -z.signum());
        render(scene, cx, &settings).unwrap().xyz[0][1]
    })
}

fn scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.1, 0.1, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

#[test]
fn regional_exterior_keeps_material_assignment_when_element_order_changes() {
    let cards = [card("upper", 1000.0, 0.0), card("lower", 2000.0, 0.0)];
    let reference = reference(&cards);
    let deformed = NODES.map(|[x, y, z]| [x + 0.1 * z, y + 0.05 * z, 1.1 * z]);
    for (tets, keys) in [(TETS, KEYS), ([TETS[1], TETS[0]], [KEYS[1], KEYS[0]])] {
        let mass = TetrahedralMass::bind(&reference, &NODES, &tets, &keys, 5, || false).unwrap();
        let body = mass.deform(&deformed, || false).unwrap();
        assert_eq!(body.region_elements(), mass.region_elements());
        with_cx(|cx| {
            let regions = deformed_body_regions(&body, BUDGET, 8, cx).unwrap();
            let whole = deformed_body_mesh(&body, BUDGET, 8, cx).unwrap();
            assert_eq!(regions.len(), 2);
            let mut actual_faces = BTreeSet::new();
            for (index, region) in regions.iter().enumerate() {
                assert_eq!(region.region_index, index);
                let element = body.region_elements()[index];
                assert_eq!(reference.regions()[index].key(), keys[element]);
                assert_eq!(region.source_nodes, tets[element]);
                assert_eq!(region.mesh.vertices.len(), 4);
                assert_eq!(region.mesh.triangles.len(), 3);
                for (&node, point) in region.source_nodes.iter().zip(&region.mesh.vertices) {
                    assert_eq!(point.map(f64::to_bits), deformed[node].map(f64::to_bits));
                }
                for triangle in &region.mesh.triangles {
                    let ids = triangle.map(|i| region.source_nodes[i as usize]);
                    let opposite = *region
                        .source_nodes
                        .iter()
                        .find(|id| !ids.contains(id))
                        .unwrap();
                    let [a, b, c] = ids.map(|id| {
                        let [x, y, z] = deformed[id];
                        MassVec3::new(x, y, z)
                    });
                    let [x, y, z] = deformed[opposite];
                    assert!((b - a).cross(c - a).dot(MassVec3::new(x, y, z) - a) < 0.0);
                    let mut key = ids;
                    key.sort_unstable();
                    assert!(actual_faces.insert(key));
                }
            }
            let expected: BTreeSet<_> = whole
                .triangles
                .iter()
                .map(|face| {
                    let mut key = face.map(|i| i as usize);
                    key.sort_unstable();
                    key
                })
                .collect();
            assert_eq!(actual_faces, expected);
            assert!(!actual_faces.contains(&[0, 1, 2]));
        });
    }
}

#[test]
fn local_accepted_temperatures_choose_the_actual_visible_thermal_region() {
    let same_card = card("one-shared-material", 2000.0, 0.0);
    let cards = [same_card.clone(), same_card];
    let reference = reference(&cards);
    let mass = TetrahedralMass::bind(&reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    with_cx(|cx| {
        let mut observations = Vec::new();
        for temperatures in [[1800.0, 300.0], [300.0, 1800.0]] {
            let body = mass
                .solve_thermal(&cells(&cards, &temperatures), &SUPPORT, solve_config(), cx)
                .unwrap();
            let appearance = appearance(&body, &cards, &temperatures);
            let primitives =
                thermal_body_primitives(&body, &bindings(&appearance), BUDGET, 8, cx).unwrap();
            assert_eq!(primitives.len(), 2);
            for (primitive, region) in primitives.iter().zip(reference.regions()) {
                let input = KEYS.iter().position(|key| *key == region.key()).unwrap();
                let Material::GrayDiffuseThermal { surface } = primitive.material else {
                    panic!("the original material thermal binding was not retained");
                };
                assert_eq!(surface.temperature_k(), temperatures[input]);
                let Shape::Mesh(mesh) = &primitive.shape else {
                    panic!("regional appearance must use the accepted material mesh");
                };
                assert_eq!(mesh.triangles.len(), 3);
            }
            // Alpha is zero, so geometry and emissivity stay fixed. Only the
            // two accepted physical temperatures exchange places.
            assert_eq!(body.body().nodes_m(), NODES);
            observations.push(pixels(&mut scene(primitives), cx));
        }
        let [top_hot, bottom_cold] = observations[0];
        let [top_cold, bottom_hot] = observations[1];
        assert!(top_hot.is_finite() && top_hot > 0.0);
        assert!(bottom_hot.is_finite() && bottom_hot > 0.0);
        assert!(top_cold >= 0.0 && top_cold < top_hot * 1.0e-6);
        assert!(bottom_cold >= 0.0 && bottom_cold < bottom_hot * 1.0e-6);
        assert!((top_hot / bottom_hot - 1.0).abs() < 1.0e-12);
    });
}

#[test]
fn regional_binding_refuses_stale_foreign_and_incomplete_state_atomically() {
    let cards = [card("upper", 2000.0, 0.0), card("lower", 2000.0, 0.0)];
    let reference = reference(&cards);
    let mass = TetrahedralMass::bind(&reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    with_cx(|cx| {
        let old = mass
            .solve_thermal(
                &cells(&cards, &[900.0, 900.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let body = mass
            .solve_thermal(
                &cells(&cards, &[1800.0, 900.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let stale = appearance(&old, &cards, &[900.0, 900.0]);
        let current = appearance(&body, &cards, &[1800.0, 900.0]);
        let good = bindings(&current);
        let mut published = scene(thermal_body_primitives(&body, &good, BUDGET, 8, cx).unwrap());
        let before = pixels(&mut published, cx);
        assert!(matches!(
            thermal_body_primitives(&body, &bindings(&stale), BUDGET, 8, cx),
            Err(MaterialRegionBindingError::StateMismatch { .. })
        ));
        assert!(matches!(
            thermal_body_primitives(&body, &good[..1], BUDGET, 8, cx),
            Err(MaterialRegionBindingError::Bindings(_))
        ));
        let duplicate = [
            MaterialRegionBinding {
                key: "lower",
                appearance: MaterialRegionAppearance::Thermal(&current[0].1),
            },
            MaterialRegionBinding {
                key: "lower",
                appearance: MaterialRegionAppearance::Thermal(&current[0].1),
            },
        ];
        assert!(matches!(
            thermal_body_primitives(&body, &duplicate, BUDGET, 8, cx),
            Err(MaterialRegionBindingError::Bindings(_))
        ));
        let foreign_cards = [card("foreign-upper", 2000.0, 0.0), cards[1].clone()];
        let foreign_reference = super::reference(&foreign_cards);
        let foreign_mass =
            TetrahedralMass::bind(&foreign_reference, &NODES, &TETS, &KEYS, 5, || false).unwrap();
        let foreign = foreign_mass
            .solve_thermal(
                &cells(&foreign_cards, &[1800.0, 900.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
            .unwrap();
        let foreign_appearance = appearance(&foreign, &foreign_cards, &[1800.0, 900.0]);
        assert!(matches!(
            thermal_body_primitives(&body, &bindings(&foreign_appearance), BUDGET, 8, cx),
            Err(MaterialRegionBindingError::StateMismatch { .. })
        ));
        assert_eq!(pixels(&mut published, cx), before);
    });
}
