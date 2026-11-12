//! Spatial pressure/context maps on the original nonlinear common refinement.
use super::*;
use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};

fn source_value(a: f64, b: f64) -> f64 {
    0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0
}

fn point(pressure: f64) -> QueryPoint {
    temperature::source_point()
        .with_quantity(
            "normal_pressure",
            fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            pressure,
        )
        .unwrap()
}

fn card(reverse: bool) -> InterfaceSystemCard {
    let dims = AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    // Two independently declared pressure charts on one ordered interface
    // card. Their selected claims and temperature gradients are different.
    for (pressure, scale) in [(1e6, 1.0), (2e6, 2.0)] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
                value: PropertyValue::Surface2D {
                    axes: ["TA".into(), "TB".into()],
                    axis_dims: [fs_qty::Temperature::DIMS; 2],
                    samples: [
                        [256.0, 256.0],
                        [384.0, 256.0],
                        [256.0, 384.0],
                        [384.0, 384.0],
                    ]
                    .into_iter()
                    .map(|p| (p, scale * source_value(p[0], p[1])))
                    .collect(),
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims,
                },
                // Typed as point() queries them, as source_card() declares:
                // fs-matdb compares an axis's quantity with the query's exactly.
                validity: ValidityDomain::unconstrained()
                    .with_quantity("TA", temperature::kelvin(), 256.0, 384.0)
                    .with_quantity("TB", temperature::kelvin(), 256.0, 384.0)
                    .with_quantity(
                        "normal_pressure",
                        fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
                        pressure,
                        pressure,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured two-pressure contact charts, not measured data".into(),
                    license: "internal-test-use".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "fixed ordered contact normal".into(),
    };
    let names = if reverse { ["B", "A"] } else { ["A", "B"] };
    InterfaceSystemCard::assemble(
        side(names[0]),
        side(names[1]),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "test".into(),
            history: "unaged".into(),
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn bind<'m, 'c>(
    cx: &fs_exec::Cx<'_>,
    mesh: &'m ConductionMesh,
    bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard,
    a: Vec<usize>,
    b: Vec<usize>,
    points: Vec<(usize, QueryPoint)>,
    reverse: bool,
) -> TemperatureContacts<'m, 'c> {
    let axes = if reverse { ["TB", "TA"] } else { ["TA", "TB"] };
    let law = TemperatureContact::new_nonmatching_mapped(
        cx,
        "bond",
        a,
        b,
        card,
        points,
        TemperatureContactCoordinates::TwoSided(axes.map(String::from)),
        SelectionPolicy::SingleClaimOnly,
        options(),
    )
    .unwrap();
    assert!(
        law.reference_face_states()
            .windows(2)
            .all(|pair| pair[0].side_a < pair[1].side_a)
    );
    assert_eq!(
        law.reference_point(),
        &law.reference_face_states()[0].reference_point
    );
    TemperatureContacts::with_nonmatching(
        cx,
        mesh,
        bc,
        Vec::new(),
        Vec::new(),
        vec![law],
        10_000,
        10_000,
    )
    .unwrap()
}

fn patch_points(mesh: &ConductionMesh, faces: &[usize]) -> Vec<(usize, QueryPoint)> {
    faces
        .iter()
        .rev()
        .map(|&face| {
            (
                face,
                point(if mesh.boundary()[face].centroid[1] < 1.5 {
                    1e6
                } else {
                    2e6
                }),
            )
        })
        .collect()
}

fn affine_field(mesh: &ConductionMesh, on_a: &[bool]) -> Vec<f64> {
    mesh.positions()
        .iter()
        .zip(on_a)
        .map(|(p, &a)| {
            let y = if p[1] > 1.5 { p[1] - 2.0 } else { p[1] };
            if a {
                344.0 + 19.0 * y
            } else {
                312.0 + 15.0 * y
            }
        })
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn g3_uniform_context_maps_preserve_original_pointwise_operators_and_receipts() {
    with_cx(|cx| {
        for right in [[1, 1, 1], [1, 2, 3]] {
            let (mesh, split) = mesh(right, 0.0);
            let bc = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let (a, b) = sides(&mesh, split, 0.0);
            let card = temperature::source_card(1.0, false);
            let uniform = temperature::bind(cx, &mesh, &bc, &card, a.clone(), b.clone(), false);
            let contexts = a.iter().rev().map(|&face| (face, point(1e6))).collect();
            let mapped = bind(cx, &mesh, &bc, &card, a, b, contexts, false);
            let on_a: Vec<_> = (0..mesh.vertex_count()).map(|i| i < split).collect();
            let t = affine_field(&mesh, &on_a);
            let old = uniform.linearize_temperature(cx, &bc, &t, 10_000).unwrap();
            let new = mapped.linearize_temperature(cx, &bc, &t, 10_000).unwrap();
            assert_eq!(old.nodal_outflow_w(), new.nodal_outflow_w());
            assert_eq!(
                old.physical_contacts().quadrature_states(),
                new.physical_contacts().quadrature_states()
            );
            let direction: Vec<_> = (0..t.len()).map(|i| (i % 7) as f64 - 3.0).collect();
            assert_eq!(
                old.apply(cx, &direction).unwrap(),
                new.apply(cx, &direction).unwrap()
            );
            assert_eq!(
                old.apply_transpose(cx, &direction).unwrap(),
                new.apply_transpose(cx, &direction).unwrap()
            );
        }
    });
}

#[test]
fn g1_each_a_face_retains_its_pressure_chart_across_independent_b_refinement() {
    with_cx(|cx| {
        for right in [[1, 1, 1], [1, 2, 3], [1, 3, 4]] {
            let (mesh, split) = mesh(right, 0.0);
            let bc = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let (a, b) = sides(&mesh, split, 0.0);
            let card = card(false);
            let scale = |face: usize| {
                let p = mesh.boundary()[face].centroid;
                if p[1] > p[2] { 2.0 } else { 1.0 }
            };
            let points: Vec<_> = a
                .iter()
                .rev()
                .map(|&face| (face, point(1e6 * scale(face))))
                .collect();
            let set = bind(cx, &mesh, &bc, &card, a.clone(), b, points, false);
            let on_a: Vec<_> = (0..mesh.vertex_count()).map(|i| i < split).collect();
            let t = affine_field(&mesh, &on_a);
            let linear = set.linearize_temperature(cx, &bc, &t, 10_000).unwrap();
            // TA-TB=32+4y and R_base=.5+y/16. The physical flux is
            // independently 64/scale W/m² on each original A triangle.
            let mut expected_a = vec![0.0; split];
            let mut heat = 0.0;
            for &face in &a {
                let triangle = &mesh.boundary()[face];
                let rate = 64.0 * triangle.area / scale(face);
                heat += rate;
                for vertex in triangle.vertices {
                    expected_a[vertex as usize] += rate / 3.0;
                }
            }
            near(heat, 48.0, 2e-12);
            for (&actual, expected) in linear.nodal_outflow_w()[..split].iter().zip(expected_a) {
                near(actual, expected, 3e-10);
            }
            near(linear.nodal_outflow_w().iter().sum(), 0.0, 3e-10);
            let physical = linear.physical_contacts();
            near(
                physical.frozen_interfaces().fluxes(&t).unwrap()[0].heat_rate_a_to_b_w,
                heat,
                3e-10,
            );
            let mut claims = std::collections::BTreeSet::new();
            for state in physical.quadrature_states() {
                let factor = scale(state.pair.side_a);
                near(
                    state.resistance.value_m2_k_per_w(),
                    factor * source_value(state.side_a_temperature_k, state.side_b_temperature_k),
                    2e-14,
                );
                assert!(
                    state
                        .resistance
                        .receipt()
                        .query_point
                        .contains(&("normal_pressure".into(), 1e6 * factor))
                );
                assert_eq!(state.resistance.card_identity(), card.content_hash());
                claims.insert(state.resistance.receipt().selected);
            }
            assert_eq!(claims.len(), 2);
        }
    });
}

#[test]
fn g3_mapped_nonlinear_tangents_and_ordered_surface_reversal_preserve_physics() {
    with_cx(|cx| {
        let (mesh, on_a) = parallel_contact_mesh([[1, 2, 3], [1, 3, 2]]);
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let (a, b) = parallel_sides(&mesh);
        let card = card(false);
        let set = bind(
            cx,
            &mesh,
            &bc,
            &card,
            a.clone(),
            b.clone(),
            patch_points(&mesh, &a),
            false,
        );
        let t = affine_field(&mesh, &on_a);
        let direction: Vec<_> = (0..t.len()).map(|i| (i % 7) as f64 / 5.0 - 0.4).collect();
        let weights: Vec<_> = (0..t.len()).map(|i| (i % 11) as f64 / 7.0 - 0.7).collect();
        let linear = set.linearize_temperature(cx, &bc, &t, 10_000).unwrap();
        let tangent = linear.apply(cx, &direction).unwrap();
        let dual = linear.apply_transpose(cx, &weights).unwrap();
        near(dot(&weights, &tangent), dot(&direction, &dual), 2e-12);
        near(tangent.iter().sum(), 0.0, 2e-12);
        let h = 1e-3;
        let plus: Vec<_> = t.iter().zip(&direction).map(|(t, d)| t + h * d).collect();
        let minus: Vec<_> = t.iter().zip(&direction).map(|(t, d)| t - h * d).collect();
        let plus = set.linearize_temperature(cx, &bc, &plus, 10_000).unwrap();
        let minus = set.linearize_temperature(cx, &bc, &minus, 10_000).unwrap();
        for ((&actual, &plus), &minus) in tangent
            .iter()
            .zip(plus.nodal_outflow_w())
            .zip(minus.nodal_outflow_w())
        {
            near(actual, (plus - minus) / (2.0 * h), 2e-8);
        }
        let reverse_card = self::card(true);
        let reverse = bind(
            cx,
            &mesh,
            &bc,
            &reverse_card,
            b.clone(),
            a,
            patch_points(&mesh, &b),
            true,
        );
        let reverse = reverse.linearize_temperature(cx, &bc, &t, 10_000).unwrap();
        for (&forward, &backward) in linear
            .nodal_outflow_w()
            .iter()
            .zip(reverse.nodal_outflow_w())
        {
            near(forward, backward, 3e-10);
        }
        for (&forward, backward) in tangent.iter().zip(reverse.apply(cx, &direction).unwrap()) {
            near(forward, backward, 3e-12);
        }
        near(
            linear
                .physical_contacts()
                .frozen_interfaces()
                .fluxes(&t)
                .unwrap()[0]
                .heat_rate_a_to_b_w,
            96.0,
            3e-10,
        );
        near(
            reverse
                .physical_contacts()
                .frozen_interfaces()
                .fluxes(&t)
                .unwrap()[0]
                .heat_rate_a_to_b_w,
            -96.0,
            3e-10,
        );
    });
}

#[test]
fn g1_mapped_contact_sources_advance_the_original_conservative_thermal_owner() {
    with_cx(|cx| {
        let (mesh, on_a) = parallel_contact_mesh([[1, 2, 3], [1, 3, 2]]);
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let (a, b) = parallel_sides(&mesh);
        let card = card(false);
        let set = bind(
            cx,
            &mesh,
            &bc,
            &card,
            a.clone(),
            b.clone(),
            patch_points(&mesh, &a),
            false,
        );
        let target: Vec<_> = on_a
            .iter()
            .map(|&a| if a { 344.0 } else { 312.0 })
            .collect();
        let mut load = vec![0.0; target.len()];
        for (faces, sign) in [(&a, 1.0), (&b, -1.0)] {
            for &slot in faces {
                let face = &mesh.boundary()[slot];
                let flux = if face.centroid[1] < 1.5 { 64.0 } else { 32.0 };
                for vertex in face.vertices {
                    load[vertex as usize] += sign * flux * face.area / 3.0;
                }
            }
        }
        let mut mass = vec![0.0; target.len()];
        for (e, tet) in mesh.complex().tets.iter().enumerate() {
            for &vertex in tet {
                mass[vertex as usize] += 10.0 * mesh.element_volume(e) / 4.0;
            }
        }
        let dt = 0.01;
        // Independent endpoint: bulk gradients vanish and each original
        // pressure chart supplies its own constant contact flux.
        let old: Vec<_> = target
            .iter()
            .zip(&load)
            .zip(&mass)
            .map(|((t, q), m)| t + dt * q / m)
            .collect();
        let before = old.clone();
        let engine =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(10.0).unwrap())
                .unwrap();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let mut linear = SolveConfig::default().linear;
        linear.tolerance = 1e-12;
        linear.max_iterations = 10_000;
        let step = StepConfig {
            linear,
            energy_tolerance_j: 1e-8,
        };
        let coupling = ContactIterationConfig {
            max_iterations: 64,
            relaxation: 1.0,
            residual_tolerance_j: 1e-9,
        };
        let next = set
            .advance_capacity(cx, &engine, problem, &old, dt, step, None, coupling)
            .unwrap();
        for (&actual, expected) in next.conduction.temperature.iter().zip(target) {
            near(actual, expected, 3e-8);
        }
        near(next.physical.stored_energy_change_j, 0.0, 1e-8);
        near(next.physical.energy_residual_j, 0.0, 1e-8);
        assert!(next.physical.residual_norm_j <= coupling.residual_tolerance_j);
        let end = engine
            .linearize_temperature_contact_endpoint(
                cx,
                problem,
                &set,
                &old,
                &next.conduction.temperature,
                dt,
                step,
                coupling.residual_tolerance_j,
                &[],
                10_000,
            )
            .unwrap();
        assert_eq!(
            end.contact().physical_contacts().quadrature_states(),
            next.physical_contacts.quadrature_states()
        );
        assert_eq!(old, before);
    });
}

#[test]
fn g0_mapped_source_coverage_selection_support_budget_and_cancellation_refuse() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let (a, b) = sides(&mesh, split, 0.0);
        let card = card(false);
        let points: Vec<_> = a.iter().map(|&face| (face, point(1e6))).collect();
        let create = |points, limits| {
            TemperatureContact::new_nonmatching_mapped(
                cx,
                "bond",
                a.clone(),
                b.clone(),
                &card,
                points,
                TemperatureContactCoordinates::TwoSided(["TA".into(), "TB".into()]),
                SelectionPolicy::SingleClaimOnly,
                limits,
            )
        };
        let mut duplicate = points.clone();
        duplicate[1].0 = duplicate[0].0;
        let mut foreign = points.clone();
        foreign[1].0 = b[0];
        let mut unsupported = points.clone();
        unsupported[1].1 = point(3e6);
        for bad in [
            Vec::new(),
            points[..1].to_vec(),
            duplicate,
            foreign,
            unsupported,
        ] {
            assert!(create(bad, options()).is_err());
        }
        assert!(
            create(
                points.clone(),
                NonmatchingOptions {
                    max_pair_tests: 1,
                    ..options()
                }
            )
            .is_err()
        );
        assert!(
            support::with_cancelled_cx(|cancelled| TemperatureContact::new_nonmatching_mapped(
                cancelled,
                "bond",
                a.clone(),
                b.clone(),
                &card,
                points.clone(),
                TemperatureContactCoordinates::TwoSided(["TA".into(), "TB".into()]),
                SelectionPolicy::SingleClaimOnly,
                options()
            ))
            .is_err()
        );
        let set = bind(cx, &mesh, &bc, &card, a.clone(), b, points, false);
        let good = vec![320.0; mesh.vertex_count()];
        let before = set.evaluate(cx, &bc, &good).unwrap();
        let mut bad = good.clone();
        for face in a {
            for vertex in mesh.boundary()[face].vertices {
                bad[vertex as usize] = 400.0;
            }
        }
        assert!(set.evaluate(cx, &bc, &bad).is_err());
        assert_eq!(
            before.quadrature_states(),
            set.evaluate(cx, &bc, &good).unwrap().quadrature_states()
        );
    });
}
