//! Material data -> accepted solid shape -> real ray hits and thermal pixels.
#[path = "material_geometry/regional.rs"]
mod regional;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::radiation::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Vec3 as MassVec3;
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_mbd_material::deformation::TetrahedralMass;
use fs_mbd_material::deformation::thermoelastic::{ThermalCell, ThermalSolveConfig};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_render::charts::{Ray, TriMesh};
use fs_render::material_geometry::{MaterialGeometryError, deformed_body_mesh, thermal_body_mesh};
use fs_render::material_state::MaterialThermalBinding;
use fs_render::tracer::{
    Camera, DirectStrategy, Primitive, Sampler, Scene, Settings, Shape, render,
};
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_solid::loads3::LoadError3;

const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 4]];
const KEYS: [&str; 2] = ["upper", "lower"];
const SUPPORT: [usize; 6] = [0, 1, 2, 4, 5, 8];
const BUDGET: TetAssemblyBudget = TetAssemblyBudget {
    maximum_nodes: 5,
    maximum_tetrahedra: 2,
    maximum_free_dofs: 15,
    minimum_scaled_jacobian: 1.0e-10,
};

fn near(got: f64, expected: f64) {
    assert!((got - expected).abs() < 2.0e-9, "{got} != {expected}");
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new_clock_free(), f)
}

fn card(name: &str, young: f64, alpha: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 6.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.25),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            alpha,
        ),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, 0.6),
    ] {
        let density = key == DENSITY_PROPERTY;
        let emissivity = key == SURFACE_EMISSIVITY_PROPERTY;
        let validity =
            ValidityDomain::unconstrained().with("T", 300.0, if density { 300.0 } else { 2000.0 });
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: if density || emissivity {
                    PropertyValue::Scalar { value, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (2000.0, value)],
                        dims,
                    }
                },
                validity: if emissivity {
                    validity
                } else {
                    validity.with("pressure", 100_000.0, 100_000.0)
                },
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if density || emissivity {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic {name}/{key}; not measured material data"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic material-geometry fixture".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn reference(cards: &[MaterialCard]) -> ResolvedBody {
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap();
    let bindings: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, material)| RegionBinding {
            key: KEYS[i],
            geometry: VolumeMoments::tetrahedron(TETS[i].map(|node| {
                let [x, y, z] = NODES[node];
                MassVec3::new(x, y, z)
            }))
            .unwrap(),
            material,
            point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    resolve_body(&bindings, 2, || false).unwrap()
}

fn cells<'a>(cards: &'a [MaterialCard], temperatures: &[f64]) -> Vec<ThermalCell<'a>> {
    cards
        .iter()
        .zip(temperatures)
        .enumerate()
        .map(|(i, (card, &temperature_k))| ThermalCell {
            stiffness_model: fs_mbd_material::deformation::thermoelastic::ThermalStiffnessModel::Isotropic,
            expansion_model:fs_mbd_material::deformation::thermoelastic::ThermalExpansionModel::Isotropic,
            key: KEYS[i],
            card,
            temperature_k,
            stiffness_selection: MaterialPropertySelection::SingleClaimOnly,
            expansion_selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect()
}

fn solve_config() -> ThermalSolveConfig {
    ThermalSolveConfig {
        assembly: BUDGET,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1.0e-11,
            maximum_iterations: 100,
        },
        mesh_frame: ContentHash([41; 32]),
        maximum_absolute_mandel_strain: 0.05,
        maximum_nodal_displacement_m: 0.1,
    }
}

fn above(x: f64, y: f64) -> Ray {
    Ray {
        origin: Point3::new(x, y, 2.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    }
}

#[test]
fn shared_boundary_omits_the_internal_face_and_winds_outward() {
    let cards = [card("upper", 1000.0, 1.0e-4), card("lower", 2000.0, 1.0e-4)];
    let source = reference(&cards);
    let mass = TetrahedralMass::bind(&source, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    let body = mass.deform(&NODES, || false).unwrap();
    let mesh = with_cx(|cx| deformed_body_mesh(&body, BUDGET, 8, cx)).unwrap();
    assert_eq!(mesh.vertices, body.nodes_m());
    assert_eq!(mesh.triangles.len(), 6);
    for triangle in &mesh.triangles {
        let ids = triangle.map(|v| v as usize);
        let owners: Vec<_> = TETS
            .iter()
            .filter(|tet| ids.iter().all(|id| tet.contains(id)))
            .collect();
        assert_eq!(owners.len(), 1, "an internal face leaked into the surface");
        let opposite = *owners[0].iter().find(|id| !ids.contains(id)).unwrap();
        let [a, b, c] = ids.map(|id| {
            let [x, y, z] = mesh.vertices[id];
            MassVec3::new(x, y, z)
        });
        let [x, y, z] = mesh.vertices[opposite];
        assert!((b - a).cross(c - a).dot(MassVec3::new(x, y, z) - a) < 0.0);
    }
    // The first ray crossing from inside the upper cell must be the LOWER
    // exterior, not the shared z=0 interface.
    let through_shared_face = Ray {
        origin: Point3::new(0.1, 0.1, 0.1),
        dir: Vec3::new(0.0, 0.0, -1.0),
    };
    let hit = mesh.intersect(&through_shared_face).unwrap();
    near(hit.point.z, -0.8);
    near(hit.t, 0.9);
    assert!(hit.normal.unwrap().z < 0.0);
}

#[test]
fn material_stiffness_changes_the_actual_thermoelastic_ray_hit() {
    let mut observed_z = Vec::new();
    for top_e in [1000.0, 2000.0] {
        let cards = [card("upper", top_e, 1.0e-4), card("lower", 2000.0, 1.0e-4)];
        let source = reference(&cards);
        let mass = TetrahedralMass::bind(&source, &NODES, &TETS, &KEYS, 5, || false).unwrap();
        let thermal = with_cx(|cx| {
            mass.solve_thermal(
                &cells(&cards, &[400.0, 300.0]),
                &SUPPORT,
                solve_config(),
                cx,
            )
        })
        .unwrap();
        let mesh = with_cx(|cx| thermal_body_mesh(&thermal, BUDGET, 8, cx)).unwrap();
        assert_eq!(mesh.vertices, thermal.body().nodes_m());
        assert_eq!(mesh.vertices, thermal.geometry().nodes_m());
        // Independent compatible two-cell P1 solution with nu=1/4. Shared
        // nodes stretch by s; the heated apex moves by 0.01+2/3*(0.01-s).
        let shared = top_e * 0.01 / (top_e + 2000.0);
        let apex = 1.0 + 0.01 + (2.0 / 3.0) * (0.01 - shared);
        let expected_z = apex * (1.0 - 0.2 / (1.0 + shared));
        let ray = above(0.1, 0.1);
        let hit = mesh.intersect(&ray).unwrap();
        near(hit.point.z, expected_z);
        near(hit.t, 2.0 - expected_z);
        assert!(hit.normal.unwrap().z > 0.0);
        near(mesh.intersect_bruteforce(&ray).unwrap().t, hit.t);
        assert!(hit.point.z > 0.8); // Undeformed top face has z=1-x-y.
        observed_z.push(hit.point.z);
    }
    assert!((observed_z[1] - observed_z[0]).abs() > 1.0e-5);
}

#[test]
fn accepted_hot_state_drives_both_real_silhouette_and_thermal_radiance() {
    let mut brightness = Vec::new();
    for alpha in [0.0, 1.0e-5] {
        let cards = [card("expansion-swap", 2000.0, alpha)];
        let source = reference(&cards);
        let mass = TetrahedralMass::bind(&source, &NODES[..4], &TETS[..1], &KEYS[..1], 4, || false)
            .unwrap();
        with_cx(|cx| {
            let thermal = mass
                .solve_thermal(&cells(&cards, &[1800.0]), &SUPPORT, solve_config(), cx)
                .unwrap();
            let mesh = thermal_body_mesh(&thermal, BUDGET, 4, cx).unwrap();
            for (original, actual) in NODES[..4].iter().zip(&mesh.vertices) {
                for axis in 0..3 {
                    near(actual[axis], original[axis] * (1.0 + alpha * 1500.0));
                }
            }
            let point = thermal.sources()[0].stiffness().query();
            let (emissivity, _) = SurfaceEmissivity::from_condition(
                "upper",
                ThermalMaterialPoint {
                    material: &cards[0],
                    context: &point,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let binding = MaterialThermalBinding::try_gray_diffuse(
                thermal.sources()[0].stiffness(),
                &emissivity,
                ContentHash([42; 32]),
            )
            .unwrap();
            // Both cards have the SAME hot radiance. Only alpha changes. The
            // physical expansion moves x+y=1 to x+y=1.015 across this pixel.
            let ray = above(0.505, 0.5);
            assert_eq!(mesh.intersect(&ray).is_some(), alpha > 0.0);
            let scene = Scene {
                primitives: vec![Primitive {
                    shape: Shape::Mesh(mesh),
                    material: binding.material(),
                    emission: None,
                }],
                lights: Vec::new(),
                environment: None,
                camera: Camera {
                    eye: ray.origin,
                    forward: ray.dir,
                    up: Vec3::new(0.0, 1.0, 0.0),
                    half_tan: 0.0,
                },
            };
            let image = render(
                &scene,
                cx,
                &Settings {
                    width: 1,
                    height: 1,
                    spp: 64,
                    max_depth: 1,
                    sampler: Sampler::OwenSobol,
                    strategy: DirectStrategy::BsdfOnly,
                    seed: 41,
                },
            )
            .unwrap();
            brightness.push(image.xyz[0][1]);
        });
    }
    assert_eq!(brightness[0], 0.0);
    assert!(brightness[1].is_finite() && brightness[1] > 0.0);
}

#[test]
fn budgets_scale_and_cancellation_refuse_before_a_scene_replacement() {
    let cards = [card("upper", 1000.0, 1.0e-4), card("lower", 2000.0, 1.0e-4)];
    let source = reference(&cards);
    let mass = TetrahedralMass::bind(&source, &NODES, &TETS, &KEYS, 5, || false).unwrap();
    let body = mass.deform(&NODES, || false).unwrap();
    let published: TriMesh = with_cx(|cx| deformed_body_mesh(&body, BUDGET, 8, cx)).unwrap();
    let before = published.bvh_fingerprint();
    assert!(matches!(
        with_cx(|cx| deformed_body_mesh(&body, BUDGET, 7, cx)),
        Err(MaterialGeometryError::Boundary(LoadError3::Budget))
    ));
    let tiny_nodes = NODES.map(|node| node.map(|x| x * 1.0e-7));
    let tiny = mass.deform(&tiny_nodes, || false).unwrap();
    assert!(matches!(
        with_cx(|cx| deformed_body_mesh(&tiny, BUDGET, 8, cx)),
        Err(MaterialGeometryError::RendererNormal { .. })
    ));
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        with_gate(&gate, |cx| deformed_body_mesh(&body, BUDGET, 8, cx)),
        Err(MaterialGeometryError::Cancelled)
    ));
    assert_eq!(published.bvh_fingerprint(), before);
    assert_eq!(published.vertices, NODES);
    near(published.intersect(&above(0.2, 0.2)).unwrap().point.z, 0.6);
    assert_eq!(body.nodes_m(), NODES);
}
