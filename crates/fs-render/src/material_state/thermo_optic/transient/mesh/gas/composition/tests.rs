use super::transport::GasCompositionMeshError;
use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, SecondOrderIndexField, propagate_bundle};
use crate::gradient_index::event::{BoundaryTermination, GradientEventOptions};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, GradientIndexPropagation};
use crate::material_state::thermo_optic::transient::mesh::bundle::P1MeshBundleOptions;
use crate::material_state::thermo_optic::transient::mesh::{
    P1ConductionMesh, P1MeshError, P1MeshOptions,
};
use fs_conduction::ConductionMesh;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::composition::MolecularCauchyDispersion;
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_rep_mesh::TetComplex;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 191,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.dot(v).sqrt())
}
fn model(equal: bool) -> MolecularRefractionModel {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    if !equal {
        specs[1].molar_mass *= 2.0;
    }
    // Manufactured optical references, with distinct original row identities.
    let optical: Vec<_> = [0.00015, if equal { 0.00015 } else { 0.0006 }]
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            MolecularRayleigh::try_new(r, 2.5e25, 550.0, 0.03, ContentHash([i as u8 + 191; 32]))
                .unwrap()
        })
        .collect();
    let source = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optical,
        None,
    )
    .unwrap();
    let original = MolecularRefractionModel::from_mixture(
        &source,
        source.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    let dispersion: Vec<_> = original
        .components()
        .iter()
        .enumerate()
        .map(|(i, row)| {
            MolecularCauchyDispersion::try_new(
                row.optical(),
                if equal || i == 0 { 2e-6 } else { 8e-6 },
                if equal || i == 0 { 1e-7 } else { 3e-7 },
                [380.0, 780.0],
                ContentHash([i as u8 + 193; 32]),
            )
            .unwrap()
        })
        .collect();
    original.with_dispersion(&dispersion).unwrap()
}
fn bipyramid(contact: bool, reverse: bool) -> ConductionMesh {
    let mut positions = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
    ];
    let mut tets = if contact {
        positions.extend([positions[0], positions[1], positions[2]]);
        vec![[0, 1, 2, 3], [5, 6, 7, 4]]
    } else {
        vec![[0, 1, 2, 3], [0, 1, 2, 4]]
    };
    if reverse {
        for tet in &mut tets {
            tet.swap(0, 1);
        }
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn options() -> P1MeshBundleOptions {
    P1MeshBundleOptions {
        mesh: P1MeshOptions {
            boundary_id: 191,
            max_internal_faces: 32,
            event: GradientEventOptions {
                integration: GradientIndexOptions {
                    initial_step_m: 0.01,
                    max_step_m: 0.05,
                    relative_tolerance: 2e-12,
                    position_tolerance_m: 2e-12,
                    direction_tolerance: 2e-12,
                    optical_path_tolerance_m: 2e-12,
                    ..GradientIndexOptions::default()
                },
                surface_tolerance_m: 1e-8,
                ..GradientEventOptions::default()
            },
        },
        position_derivative_tolerance_m: 2e-12,
        direction_derivative_tolerance: 2e-12,
    }
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(-0.06, -0.04, -0.6),
        dir: unit(Vec3::new(0.18, 0.07, 1.0)),
    }
}
fn same_path(a: &GradientIndexPropagation, b: &GradientIndexPropagation, tolerance: f64) {
    vector_near(
        a.ray.origin.delta_from(b.ray.origin),
        Vec3::new(0.0, 0.0, 0.0),
        tolerance,
    );
    vector_near(a.ray.dir, b.ray.dir, tolerance);
    near(a.optical_path_length_m, b.optical_path_length_m, tolerance);
    near(
        a.group_optical_path_length_m.unwrap(),
        b.group_optical_path_length_m.unwrap(),
        tolerance,
    );
    assert_eq!(a.power_transmittance(), 1.0);
}
const EMPTY: GasPartialDensitySample = GasPartialDensitySample {
    density: 0.0,
    gradient_per_m: Vec3::new(0.0, 0.0, 0.0),
};
const EMPTY_JET: GasPartialDensityJet = GasPartialDensityJet {
    sample: EMPTY,
    hessian_diagonal_per_m2: [0.0; 3],
    hessian_off_diagonal_per_m2: [0.0; 3],
};

#[test]
fn g1_original_species_rows_supply_native_gradients_and_nonlinear_optical_hessians() {
    with_cx(|_, cx| {
        let model = model(false);
        let position = Point3::new(0.07, -0.08, -0.27);
        let gradient = [0.2, -0.1, 0.05];
        let fraction = |p: Point3| 1.0 + gradient[0] * p.x + gradient[1] * p.y + gradient[2] * p.z;
        for reverse in [false, true] {
            let geometry = bipyramid(false, reverse);
            let temperature = [300.0; 5];
            let thermal =
                P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &temperature).unwrap();
            let mut first_index = None;
            for basis in [
                PartialDensityBasis::MassKgPerM3,
                PartialDensityBasis::MoleculesPerM3,
            ] {
                let scales = core::array::from_fn::<_, 2, _>(|i| match basis {
                    PartialDensityBasis::MassKgPerM3 => 1.0,
                    PartialDensityBasis::MoleculesPerM3 => {
                        AVOGADRO_PER_MOL / model.components()[i].molar_mass_kg_mol()
                    }
                });
                let values: [Vec<f64>; 2] = core::array::from_fn(|i| {
                    geometry
                        .positions()
                        .iter()
                        .map(|p| {
                            let a = fraction(Point3::new(p[0], p[1], p[2]));
                            (if i == 0 { a } else { 2.0 - a }) * scales[i]
                        })
                        .collect()
                });
                let rows = [values[0].as_slice(), values[1].as_slice()];
                let snapshot = P1GasCompositionSnapshot::try_new(
                    cx,
                    thermal,
                    &model,
                    &rows,
                    basis,
                    ContentHash([195; 32]),
                )
                .unwrap();
                assert!(core::ptr::eq(
                    snapshot.nodal_partial_densities()[0],
                    rows[0]
                ));
                assert_eq!(
                    snapshot.model().content_identity(),
                    model.content_identity()
                );
                let element = snapshot.element(0).unwrap();
                let mut samples = [EMPTY_JET; 2];
                element.sample_jets(position, &mut samples).unwrap();
                for (i, sample) in samples.into_iter().enumerate() {
                    let a = fraction(position);
                    near(
                        sample.sample.density / scales[i],
                        if i == 0 { a } else { 2.0 - a },
                        5e-15,
                    );
                    let sign = if i == 0 { 1.0 } else { -1.0 };
                    vector_near(
                        sample.sample.gradient_per_m.scale(1.0 / scales[i]),
                        Vec3::new(gradient[0], gradient[1], gradient[2]).scale(sign),
                        5e-15,
                    );
                    assert_eq!(sample.hessian_diagonal_per_m2, [0.0; 3]);
                    assert_eq!(sample.hessian_off_diagonal_per_m2, [0.0; 3]);
                }
                let field = element.bind_refraction().unwrap();
                let jet = field.sample_second_order(position, 430.0).unwrap();
                assert_eq!(jet.optical(), field.sample(position, 430.0).unwrap());
                let scalar = |p: [f64; 3]| {
                    let a = fraction(Point3::new(p[0], p[1], p[2]));
                    model
                        .refractive_index(&[a * scales[0], (2.0 - a) * scales[1]], basis, 430.0)
                        .unwrap()
                        .phase_index
                };
                let p = [position.x, position.y, position.z];
                let h = 0.02;
                let g = jet.optical().phase_gradient_per_m();
                let diagonal = jet.hessian().diagonal();
                for axis in 0..3 {
                    let mut lower = p;
                    let mut upper = p;
                    lower[axis] -= h;
                    upper[axis] += h;
                    near(
                        [g.x, g.y, g.z][axis],
                        (scalar(upper) - scalar(lower)) / (2.0 * h),
                        2e-11,
                    );
                    near(
                        diagonal[axis],
                        (scalar(upper) - 2.0 * scalar(p) + scalar(lower)) / (h * h),
                        2e-11,
                    );
                    assert!(diagonal[axis] > 0.0);
                }
                for (entry, (a, b)) in
                    jet.hessian()
                        .off_diagonal()
                        .into_iter()
                        .zip([(0, 1), (0, 2), (1, 2)])
                {
                    let mut difference = 0.0;
                    for (sa, sb) in [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)] {
                        let mut point = p;
                        point[a] += sa * h;
                        point[b] += sb * h;
                        difference += sa * sb * scalar(point);
                    }
                    near(entry, difference / (4.0 * h * h), 2e-11);
                }
                if let Some((index, identity)) = first_index {
                    near(jet.optical().phase_index(), index, 3e-15);
                    assert_ne!(snapshot.content_identity(), identity);
                } else {
                    first_index = Some((jet.optical().phase_index(), snapshot.content_identity()));
                }
            }
        }
    });
}

// Independent global affine reconstruction, without mesh coordinates or caches.
struct AffineMixture<'a> {
    model: &'a MolecularRefractionModel,
}
impl GasCompositionGradientField for AffineMixture<'_> {
    fn model_identity(&self) -> ContentHash {
        self.model.content_identity()
    }
    fn component_count(&self) -> usize {
        2
    }
    fn density_basis(&self) -> PartialDensityBasis {
        PartialDensityBasis::MassKgPerM3
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        for (i, row) in output.iter_mut().enumerate() {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            *row = GasPartialDensitySample {
                density: 1.0 + sign * 0.3 * p.x,
                gradient_per_m: Vec3::new(sign * 0.3, 0.0, 0.0),
            };
        }
        Ok(())
    }
}
impl GasCompositionSecondOrderField for AffineMixture<'_> {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        let mut samples = [EMPTY; 2];
        self.sample_partial_densities(p, &mut samples)?;
        for (out, sample) in output.iter_mut().zip(samples) {
            *out = GasPartialDensityJet {
                sample,
                ..EMPTY_JET
            };
        }
        Ok(())
    }
}

#[test]
fn g1_constant_total_density_still_gives_chromatic_mesh_bending_and_group_delay() {
    with_cx(|_, cx| {
        let model = model(false);
        let independent = AffineMixture { model: &model };
        let reference = MaterialGasCompositionRefractiveField::try_new(
            &model,
            &independent,
            ContentHash([196; 32]),
        )
        .unwrap();
        for reverse in [false, true] {
            let geometry = bipyramid(false, reverse);
            let temperature = [300.0; 5];
            let thermal =
                P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &temperature).unwrap();
            let values: [Vec<f64>; 2] = core::array::from_fn(|i| {
                geometry
                    .positions()
                    .iter()
                    .map(|p| 1.0 + if i == 0 { 0.3 * p[0] } else { -0.3 * p[0] })
                    .collect()
            });
            for vertex in 0..geometry.vertex_count() {
                assert_eq!(values[0][vertex] + values[1][vertex], 2.0);
            }
            let rows = [values[0].as_slice(), values[1].as_slice()];
            let source = P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &rows,
                PartialDensityBasis::MassKgPerM3,
                ContentHash([197; 32]),
            )
            .unwrap();
            let element = source.element(0).unwrap();
            let mesh = P1ConductionMesh::try_new(cx, thermal).unwrap();
            let mut chromatic_direction = Vec::new();
            for wavelength in [400.0, 550.0, 750.0] {
                let result = mesh
                    .propagate_gas_composition_bundle(
                        cx,
                        &element,
                        ray(),
                        Vec3::new(1.0, 0.0, 0.0),
                        wavelength,
                        1.0,
                        options(),
                    )
                    .unwrap();
                let expected = propagate_bundle(
                    &reference,
                    ray(),
                    Vec3::new(1.0, 0.0, 0.0),
                    wavelength,
                    1.0,
                    BundleOptions {
                        integration: options().mesh.event.integration,
                        position_derivative_tolerance_m: options().position_derivative_tolerance_m,
                        direction_derivative_tolerance: options().direction_derivative_tolerance,
                    },
                )
                .unwrap();
                let actual = result.bundle.propagation();
                same_path(actual, expected.propagation(), 5e-10);
                for (a, b) in result
                    .bundle
                    .angular_derivatives()
                    .into_iter()
                    .zip(expected.angular_derivatives())
                {
                    vector_near(a.position_m_per_rad, b.position_m_per_rad, 5e-9);
                    vector_near(a.direction_per_rad, b.direction_per_rad, 5e-9);
                }
                assert_eq!(result.internal_faces, 1);
                assert_eq!(result.bundle.interface_events(), 0);
                assert!(actual.ray.dir.x < ray().dir.x - 1e-5);
                assert!(actual.group_optical_path_length_m.unwrap() > actual.optical_path_length_m);
                // Translational symmetry in y,z preserves those optical momenta.
                let start = reference
                    .sample(ray().origin, wavelength)
                    .unwrap()
                    .phase_index();
                let end = reference
                    .sample(actual.ray.origin, wavelength)
                    .unwrap()
                    .phase_index();
                near(end * actual.ray.dir.z, start * ray().dir.z, 2e-11);
                near(end * actual.ray.dir.y, start * ray().dir.y, 2e-11);
                let primal = mesh
                    .propagate_gas_composition(
                        cx,
                        &element,
                        ray(),
                        Vec3::new(1.0, 0.0, 0.0),
                        wavelength,
                        1.0,
                        options().mesh,
                    )
                    .unwrap();
                same_path(&primal.propagation, actual, 5e-10);
                chromatic_direction.push(actual.ray.dir.x);
            }
            assert!(chromatic_direction[2] > chromatic_direction[0] + 1e-7);
        }
    });
}

fn jumping_density(p: [f64; 3]) -> f64 {
    1.0 + 0.24 * p[0] + 0.04 * p[1] + (if p[2] < 0.0 { -0.13 } else { 0.34 }) * p[2]
}

#[test]
fn g1_species_gradient_jump_preserves_the_original_launch_angle_family() {
    with_cx(|_, cx| {
        let model = model(false);
        let geometry = bipyramid(false, false);
        let temperature = [300.0; 5];
        let thermal =
            P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &temperature).unwrap();
        let values: [Vec<f64>; 2] = core::array::from_fn(|i| {
            geometry
                .positions()
                .iter()
                .map(|&p| {
                    let a = jumping_density(p);
                    if i == 0 { a } else { 2.0 - a }
                })
                .collect()
        });
        let rows = [values[0].as_slice(), values[1].as_slice()];
        let source = P1GasCompositionSnapshot::try_new(
            cx,
            thermal,
            &model,
            &rows,
            PartialDensityBasis::MassKgPerM3,
            ContentHash([198; 32]),
        )
        .unwrap();
        let element = source.element(0).unwrap();
        let mesh = P1ConductionMesh::try_new(cx, thermal).unwrap();
        let result = mesh
            .propagate_gas_composition_bundle(
                cx,
                &element,
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                430.0,
                1.0,
                options(),
            )
            .unwrap();
        assert_eq!(result.internal_faces, 1);
        assert_eq!(result.element_index, 1);
        assert_eq!(result.bundle.interface_events(), 0);
        let epsilon = 2e-5;
        for (basis, actual) in result
            .bundle
            .launch_basis()
            .into_iter()
            .zip(result.bundle.angular_derivatives())
        {
            let neighbour = |sign: f64| {
                let d = ray().dir;
                let direction = unit(Vec3::new(
                    d.x + sign * epsilon * basis.x,
                    d.y + sign * epsilon * basis.y,
                    d.z + sign * epsilon * basis.z,
                ));
                mesh.propagate_gas_composition(
                    cx,
                    &element,
                    Ray {
                        dir: direction,
                        ..ray()
                    },
                    Vec3::new(1.0, 0.0, 0.0),
                    430.0,
                    1.0,
                    options().mesh,
                )
                .unwrap()
                .propagation
                .ray
            };
            let (upper, lower) = (neighbour(1.0), neighbour(-1.0));
            vector_near(
                actual.position_m_per_rad,
                upper.origin.delta_from(lower.origin).scale(0.5 / epsilon),
                5e-7,
            );
            vector_near(
                actual.direction_per_rad,
                Vec3::new(
                    upper.dir.x - lower.dir.x,
                    upper.dir.y - lower.dir.y,
                    upper.dir.z - lower.dir.z,
                )
                .scale(0.5 / epsilon),
                5e-7,
            );
        }
    });
}

#[test]
fn g1_compensating_species_gradients_preserve_constant_index_and_vacuum() {
    with_cx(|_, cx| {
        let model = model(true);
        let geometry = bipyramid(false, false);
        let temperature = [300.0; 5];
        let thermal =
            P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &temperature).unwrap();
        let mesh = P1ConductionMesh::try_new(cx, thermal).unwrap();
        for vacuum in [false, true] {
            let values: [Vec<f64>; 2] = core::array::from_fn(|i| {
                geometry
                    .positions()
                    .iter()
                    .map(|&p| {
                        if vacuum {
                            0.0
                        } else {
                            let a = jumping_density(p);
                            if i == 0 { a } else { 2.0 - a }
                        }
                    })
                    .collect()
            });
            let rows = [values[0].as_slice(), values[1].as_slice()];
            let source = P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &rows,
                PartialDensityBasis::MassKgPerM3,
                ContentHash([199; 32]),
            )
            .unwrap();
            let element = source.element(0).unwrap();
            let result = mesh
                .propagate_gas_composition_bundle(
                    cx,
                    &element,
                    ray(),
                    Vec3::new(1.0, 0.0, 0.0),
                    550.0,
                    1.0,
                    options(),
                )
                .unwrap();
            let propagation = result.bundle.propagation();
            vector_near(
                propagation.ray.origin.delta_from(ray().origin),
                ray().dir,
                5e-10,
            );
            vector_near(propagation.ray.dir, ray().dir, 5e-10);
            let density = if vacuum { 0.0 } else { 1.0 };
            let index = model
                .refractive_index(&[density; 2], PartialDensityBasis::MassKgPerM3, 550.0)
                .unwrap();
            near(propagation.optical_path_length_m, index.phase_index, 5e-10);
            near(
                propagation.group_optical_path_length_m.unwrap(),
                index.group_index,
                5e-10,
            );
            for (basis, derivative) in result
                .bundle
                .launch_basis()
                .into_iter()
                .zip(result.bundle.angular_derivatives())
            {
                vector_near(derivative.position_m_per_rad, basis, 5e-9);
                vector_near(derivative.direction_per_rad, basis, 5e-9);
            }
            assert_eq!(result.internal_faces, 1);
            assert_eq!(result.bundle.interface_events(), 0);
        }
    });
}

#[test]
fn g0_g3_original_rows_geometry_support_work_and_cancellation_remain_explicit() {
    with_cx(|gate, cx| {
        let model = model(false);
        let geometry = bipyramid(false, false);
        let temperature = [300.0; 5];
        let thermal =
            P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &temperature).unwrap();
        let a = [0.0, 0.7, 0.2, 1.0, 0.3];
        let b = [1.0; 5];
        let rows = [a.as_slice(), b.as_slice()];
        let identity = ContentHash([200; 32]);
        let source = P1GasCompositionSnapshot::try_new(
            cx,
            thermal,
            &model,
            &rows,
            PartialDensityBasis::MassKgPerM3,
            identity,
        )
        .unwrap();
        let element = source.element(0).unwrap();
        let mesh = P1ConductionMesh::try_new(cx, thermal).unwrap();
        let axis = Vec3::new(1.0, 0.0, 0.0);
        assert!(matches!(
            P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &rows[..1],
                PartialDensityBasis::MassKgPerM3,
                identity
            ),
            Err(GasCompositionSnapshotError::ComponentCount { .. })
        ));
        let short = [a[..4].as_ref(), b.as_slice()];
        assert!(matches!(
            P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &short,
                PartialDensityBasis::MassKgPerM3,
                identity
            ),
            Err(GasCompositionSnapshotError::FieldLength { component: 0, .. })
        ));
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let mut invalid = a;
            invalid[2] = bad;
            let rows = [invalid.as_slice(), b.as_slice()];
            assert!(matches!(
                P1GasCompositionSnapshot::try_new(
                    cx,
                    thermal,
                    &model,
                    &rows,
                    PartialDensityBasis::MassKgPerM3,
                    identity
                ),
                Err(GasCompositionSnapshotError::InvalidDensity {
                    component: 0,
                    vertex: 2
                })
            ));
        }
        assert!(matches!(
            P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &rows,
                PartialDensityBasis::MassKgPerM3,
                ContentHash([0; 32])
            ),
            Err(GasCompositionSnapshotError::MissingSnapshotIdentity)
        ));
        let mut samples = [EMPTY; 2];
        for (vertex, p) in element
            .geometry_element()
            .vertex_positions_m()
            .iter()
            .enumerate()
        {
            element.sample(*p, &mut samples).unwrap();
            let original = geometry.complex().tets[0][vertex] as usize;
            assert_eq!(samples[0].density, a[original]);
            assert_eq!(samples[1].density, 1.0);
            assert_eq!(samples[1].gradient_per_m, Vec3::new(0.0, 0.0, 0.0));
        }
        let before = samples;
        assert!(element.sample(ray().origin, &mut samples[..1]).is_err());
        assert_eq!(samples, before);
        assert!(
            element
                .sample(Point3::new(0.0, 0.0, 0.1), &mut samples)
                .is_err()
        );
        assert_eq!(samples, before);
        let foreign_model = self::model(true);
        assert!(
            MaterialGasCompositionRefractiveField::try_new(&foreign_model, &element, identity)
                .is_err()
        );
        let copied = temperature.to_vec();
        let foreign_geometry =
            P1ConductionTemperatureSnapshot::try_new(cx, &geometry, &copied).unwrap();
        let foreign = P1GasCompositionSnapshot::try_new(
            cx,
            foreign_geometry,
            &model,
            &rows,
            PartialDensityBasis::MassKgPerM3,
            identity,
        )
        .unwrap()
        .element(0)
        .unwrap();
        assert!(matches!(
            mesh.propagate_gas_composition_bundle(cx, &foreign, ray(), axis, 550.0, 1.0, options()),
            Err(GasCompositionMeshError::Mesh(P1MeshError::ForeignSnapshot))
        ));
        let mut limited = options();
        limited.mesh.max_internal_faces = 0;
        assert!(matches!(
            mesh.propagate_gas_composition_bundle(cx, &element, ray(), axis, 550.0, 1.0, limited),
            Err(GasCompositionMeshError::Mesh(
                P1MeshError::FaceBudgetExceeded
            ))
        ));
        assert!(
            mesh.propagate_gas_composition(cx, &element, ray(), axis, 900.0, 1.0, options().mesh)
                .is_err()
        );
        let contact = bipyramid(true, false);
        let contact_temperature = vec![300.0; contact.vertex_count()];
        let contact_thermal =
            P1ConductionTemperatureSnapshot::try_new(cx, &contact, &contact_temperature).unwrap();
        let left = vec![1.0; contact.vertex_count()];
        let mut right = vec![0.0; contact.vertex_count()];
        for vertex in [4, 5, 6, 7] {
            right[vertex] = 2.0;
        }
        let contact_rows = [left.as_slice(), right.as_slice()];
        let contact_source = P1GasCompositionSnapshot::try_new(
            cx,
            contact_thermal,
            &model,
            &contact_rows,
            PartialDensityBasis::MassKgPerM3,
            identity,
        )
        .unwrap();
        let contact_element = contact_source.element(0).unwrap();
        let contact_mesh = P1ConductionMesh::try_new(cx, contact_thermal).unwrap();
        let stop = contact_mesh
            .propagate_gas_composition_bundle(
                cx,
                &contact_element,
                ray(),
                axis,
                550.0,
                1.0,
                options(),
            )
            .unwrap();
        assert!(matches!(stop.termination, BoundaryTermination::Surface(_)));
        assert_eq!(stop.internal_faces, 0);
        assert_eq!(stop.bundle.interface_events(), 0);
        gate.request();
        assert!(matches!(
            P1GasCompositionSnapshot::try_new(
                cx,
                thermal,
                &model,
                &rows,
                PartialDensityBasis::MassKgPerM3,
                identity
            ),
            Err(GasCompositionSnapshotError::Cancelled)
        ));
        assert!(
            mesh.propagate_gas_composition_bundle(cx, &element, ray(), axis, 550.0, 1.0, options())
                .is_err()
        );
    });
}
