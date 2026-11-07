//! G1/G3/G4: original frame force envelopes, material layer sharing and refusals.
//! Synthetic constitutive sources verify implementation, not measured failure.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::strength::FrameTensileError;
use fs_solid::frame2d::{FrameMassFormulation, FrameViscousSection, SectionModel};
use fs_solid::{Frame2d, Geometry, Rayleigh};

const RADIUS: f64 = 0.02;
const LIMIT: &str = "specified_tensile_limit";

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 761,
                kernel_id: 761,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn section(tag: &str, data: &[(f64, f64)]) -> ConcentricElasticSection {
    rate_section(tag, data, &[], None)
}

fn rate_section(
    tag: &str,
    data: &[(f64, f64)],
    viscosities: &[(f64, f64)],
    shear: Option<(f64, f64)>,
) -> ConcentricElasticSection {
    assert!(viscosities.is_empty() || viscosities.len() == data.len());
    let mut layers = Vec::new();
    for (index, &(young, strength)) in data.iter().enumerate() {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        let mut properties = vec![
            ("density", 1000.0, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            (LIMIT, strength, Pressure::DIMS),
        ];
        if let Some(&(axial, bending)) = viscosities.get(index) {
            properties.extend([
                ("axial_eta", axial, DynViscosity::DIMS),
                ("bending_eta", bending, DynViscosity::DIMS),
            ]);
        }
        if let Some((modulus, viscosity)) = shear {
            properties.extend([
                ("shear_modulus", modulus, Pressure::DIMS),
                ("shear_eta", viscosity, DynViscosity::DIMS),
            ]);
        }
        for (name, value, dims) in properties {
            let validity = if name.ends_with("_eta") {
                ValidityDomain::unconstrained().with("omega", 0.0, 1e7)
            } else {
                ValidityDomain::unconstrained()
            };
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity,
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: format!("synthetic frame tensile {tag}:{index}"),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
            requirements.push(
                ScalarPropertyRequirement::try_new(
                    name,
                    dims,
                    if name.ends_with("_eta") {
                        ScalarAdmissibility::Finite
                    } else {
                        ScalarAdmissibility::StrictlyPositive
                    },
                )
                .unwrap(),
            );
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic {tag}:{index}"),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let point = if viscosities.is_empty() {
            QueryPoint::new()
        } else {
            QueryPoint::new().with("omega", 0.0).unwrap()
        };
        let material = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        layers.push(ConcentricMaterialLayer {
            outer_radius_m: RADIUS * ((index + 1) as f64 / data.len() as f64).sqrt(),
            material,
        });
    }
    ConcentricElasticSection::try_new(0.0, layers).unwrap()
}

fn cantilever(sections: &[ConcentricElasticSection], angle: f64, geometry: Geometry) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(geometry, FrameMassFormulation::ConsistentReference);
    let mut previous = frame.add_node(0.0, 0.0);
    frame.fix(previous, [true; 3]);
    for (index, section) in sections.iter().enumerate() {
        let length = (index + 1) as f64 / sections.len() as f64;
        let next = frame.add_node(length * angle.cos(), length * angle.sin());
        frame
            .add_material_element(previous, next, 5, section.clone())
            .unwrap();
        previous = next;
    }
    frame
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1e-8),
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn g1_compression_and_reversed_bending_preserve_tensile_outer_fibers() {
    with_cx(false, |cx| {
        let source = section("layers", &[(2e9, 4e5), (8e9, 5e4)]);
        let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
        for sign in [-1.0, 1.0] {
            let mut frame = cantilever(&[source.clone()], 0.0, Geometry::Linear);
            let mut load = vec![0.0; frame.ndof()];
            load[3] = -100.0;
            load[4] = sign;
            frame.static_load(&load, 1).unwrap();
            let result = frame.assess_tensile_strength(&capacities, 4, cx).unwrap();
            assert_eq!(result.samples().len(), 2);
            let root = &result.samples()[0];
            near(root.axial_force_n, -100.0);
            near(root.bending_moment_n_m, sign);
            for (index, layer) in source.layer_integrals().iter().enumerate() {
                let radius = source.layers()[index].outer_radius_m;
                let center = layer.young_modulus_pa() * (-100.0 / source.axial_stiffness_n());
                let bending = layer.young_modulus_pa() * radius / source.bending_stiffness_n_m2();
                near(
                    root.assessment.layers()[index].maximum_axial_stress_pa,
                    center + bending,
                );
                near(
                    root.assessment.layers()[index].minimum_axial_stress_pa,
                    center - bending,
                );
            }
            assert_eq!(root.assessment.governing_layer(), 1);
            assert!(result.reaches_limit());
            assert_eq!(result.governing_sample().unwrap().end, 0);
            near(result.samples()[1].assessment.utilization(), 0.0);
            assert_eq!(root.capacity, &capacities[0]);
        }
    });
}

#[test]
fn g3_heterogeneous_members_use_local_strength_and_ignore_rigid_rotation() {
    with_cx(false, |cx| {
        let sources = [
            section("strong", &[(2e9, 1e6)]),
            section("weak", &[(2e9, 5e4)]),
        ];
        let capacities: Vec<_> = sources
            .iter()
            .map(|s| s.axial_tensile_capacity(LIMIT).unwrap())
            .collect();
        let area = std::f64::consts::PI * RADIUS * RADIUS;
        for angle in [0.0, 0.731] {
            let mut frame = cantilever(&sources, angle, Geometry::Corotational);
            let mut load = vec![0.0; frame.ndof()];
            load[6] = 100.0 * angle.cos();
            load[7] = 100.0 * angle.sin();
            frame.static_load(&load, 1).unwrap();
            let result = frame.assess_tensile_strength(&capacities, 4, cx).unwrap();
            let governing = result.governing_sample().unwrap();
            assert_eq!(governing.element, 1);
            near(governing.assessment.utilization(), 100.0 / (area * 5e4));
            for sample in result.samples() {
                near(
                    sample.assessment.layers()[0].maximum_axial_stress_pa,
                    100.0 / area,
                );
            }
            assert!(result.reaches_limit());
        }
    });
}

#[test]
fn g1_native_newmark_double_curvature_is_not_averaged_to_zero() {
    with_cx(false, |cx| {
        let source = section("dynamic", &[(2e9, 100.0)]);
        let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
        let mut frame = cantilever(&[source], 0.0, Geometry::Linear);
        frame.fix(1, [true, false, true]);
        let mut force = vec![0.0; frame.ndof()];
        force[4] = 0.1;
        frame
            .newmark_with_load(
                &[0.0; 30],
                0.001,
                (0.0, 0.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
                &force,
            )
            .unwrap();
        let displacement = frame.displacements()[4];
        assert!(displacement.abs() > 1e-7);
        let result = frame.assess_tensile_strength(&capacities, 2, cx).unwrap();
        let a = &result.samples()[0];
        let b = &result.samples()[1];
        near(a.bending_moment_n_m, -b.bending_moment_n_m);
        // Fixed-guided cubic displacement: endpoint curvature is +/- 6*y/L^2.
        let peak = 6.0 * 2e9 * RADIUS * displacement.abs();
        near(a.assessment.layers()[0].maximum_axial_stress_pa, peak);
        near(b.assessment.layers()[0].maximum_axial_stress_pa, peak);
        assert!(result.reaches_limit());
    });
}

#[test]
fn g1_accepted_axial_rates_preserve_distinct_layer_viscosities_and_orientation() {
    with_cx(false, |cx| {
        let data = [(2e9, 1e3), (8e9, 1e3)];
        let viscosities = [(2e6, 0.0), (1e5, 0.0)];
        let source = rate_section("axial-rate", &data, &viscosities, None);
        let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
        let area = std::f64::consts::PI * RADIUS.powi(2);
        let mass = 1000.0 * area / 3.0;
        let stiffness = 0.5 * (data[0].0 + data[1].0) * area;
        let damping = 0.5 * (viscosities[0].0 + viscosities[1].0) * area;
        let dt = 0.02 / (stiffness / mass).sqrt();
        let force = 100.0;
        // Independent constant-force, single-coordinate Newmark solution.
        let denominator = mass + 0.5 * dt * damping + 0.25 * dt * dt * stiffness;
        let expected_velocity = force * dt / denominator;
        let expected_displacement = 0.5 * dt * expected_velocity;
        for angle in [0.0_f64, 0.731] {
            let (sine, cosine) = angle.sin_cos();
            let mut frame = Frame2d::with_mass_formulation(
                Geometry::Linear,
                FrameMassFormulation::ConsistentReference,
            );
            frame.add_node(0.0, 0.0);
            frame.add_node(cosine, sine);
            frame.fix(0, [true; 3]);
            frame.fix(1, [false, false, true]);
            frame
                .add_material_viscous_element(
                    0,
                    1,
                    5,
                    FrameViscousSection::try_new(source.clone(), "axial_eta", "bending_eta")
                        .unwrap(),
                )
                .unwrap();
            let load = [0.0, 0.0, 0.0, force * cosine, force * sine, 0.0];
            frame
                .newmark_with_load(
                    &[0.0; 2],
                    dt,
                    (0.0, 0.0),
                    Rayleigh { a0: 0.0, a1: 0.0 },
                    &load,
                )
                .unwrap();
            let displacement = cosine * frame.displacements()[3] + sine * frame.displacements()[4];
            let velocity = cosine * frame.velocities()[3] + sine * frame.velocities()[4];
            near(displacement, expected_displacement);
            near(velocity, expected_velocity);
            let result = frame.assess_tensile_strength(&capacities, 4, cx).unwrap();
            for sample in result.samples() {
                for layer in 0..data.len() {
                    let expected = data[layer].0 * expected_displacement
                        + viscosities[layer].0 * expected_velocity;
                    near(
                        sample.assessment.layers()[layer].maximum_axial_stress_pa,
                        expected,
                    );
                    near(
                        sample.assessment.layers()[layer].minimum_axial_stress_pa,
                        expected,
                    );
                }
                near(
                    sample.axial_force_n,
                    stiffness * expected_displacement + damping * expected_velocity,
                );
                // The softer inner layer governs because its own eta is larger.
                assert_eq!(sample.assessment.governing_layer(), 0);
            }
            assert!(result.reaches_limit());
            let accepted_u = frame.displacements().to_vec();
            let accepted_v = frame.velocities().to_vec();
            assert_eq!(
                frame.assess_tensile_strength(&capacities, 3, cx),
                Err(FrameTensileError::Budget)
            );
            with_cx(true, |cancelled| {
                assert_eq!(
                    frame.assess_tensile_strength(&capacities, 4, cancelled),
                    Err(FrameTensileError::Cancelled)
                );
            });
            assert_eq!(frame.displacements(), accepted_u);
            assert_eq!(frame.velocities(), accepted_v);
            assert_eq!(
                frame.assess_tensile_strength(&capacities, 4, cx).unwrap(),
                result
            );
        }
    });
}

#[test]
fn g1_unloading_combines_signed_elastic_and_viscous_bending_before_extrema() {
    with_cx(false, |cx| {
        let data = [(2e9, 1e3), (8e9, 5e3)];
        let viscosities = [(0.0, 1e8), (0.0, 1e5)];
        let source = rate_section("bending-unload", &data, &viscosities, None);
        let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
        let area = std::f64::consts::PI * RADIUS.powi(2);
        let inertia = area * RADIUS.powi(2) / 4.0;
        let ei = inertia * (0.25 * data[0].0 + 0.75 * data[1].0);
        let eta_i = inertia * (0.25 * viscosities[0].1 + 0.75 * viscosities[1].1);
        let mass = 1000.0 * (area * 13.0 / 35.0 + 6.0 * inertia / 5.0);
        let (stiffness, damping) = (12.0 * ei, 12.0 * eta_i);
        let dt = 0.001;
        let denominator = mass + 0.5 * dt * damping + 0.25 * dt * dt * stiffness;
        for sign in [-1.0, 1.0] {
            let initial = sign * 1e-5;
            let mut frame = Frame2d::with_mass_formulation(
                Geometry::Linear,
                FrameMassFormulation::ConsistentReference,
            );
            frame.add_node(0.0, 0.0);
            frame.add_node(1.0, 0.0);
            frame.fix(0, [true; 3]);
            frame.fix(1, [true, false, true]);
            frame
                .add_material_viscous_element(
                    0,
                    1,
                    5,
                    FrameViscousSection::try_new(source.clone(), "axial_eta", "bending_eta")
                        .unwrap(),
                )
                .unwrap();
            frame
                .static_load(&[0.0, 0.0, 0.0, 0.0, stiffness * initial, 0.0], 1)
                .unwrap();
            frame
                .newmark_with_load(
                    &[0.0; 2],
                    dt,
                    (0.0, 0.0),
                    Rayleigh { a0: 0.0, a1: 0.0 },
                    &[0.0; 6],
                )
                .unwrap();
            let velocity = -initial * stiffness * dt / denominator;
            let displacement = initial + 0.5 * dt * velocity;
            near(frame.velocities()[4], velocity);
            near(frame.displacements()[4], displacement);
            assert!(displacement * velocity < 0.0);
            let inner_elastic = data[0].0 * displacement;
            let inner_viscous = viscosities[0].1 * velocity;
            assert!(
                (inner_elastic + inner_viscous).abs()
                    < 0.05 * (inner_elastic.abs() + inner_viscous.abs())
            );
            let result = frame.assess_tensile_strength(&capacities, 4, cx).unwrap();
            for sample in result.samples() {
                for layer in 0..data.len() {
                    let radius = RADIUS * ((layer + 1) as f64 / data.len() as f64).sqrt();
                    let peak = 6.0
                        * radius
                        * (data[layer].0 * displacement + viscosities[layer].1 * velocity).abs();
                    near(
                        sample.assessment.layers()[layer].maximum_axial_stress_pa,
                        peak,
                    );
                    near(
                        sample.assessment.layers()[layer].minimum_axial_stress_pa,
                        -peak,
                    );
                }
                let direction = if sample.end == 0 { 1.0 } else { -1.0 };
                near(
                    sample.bending_moment_n_m,
                    direction * 6.0 * (ei * displacement + eta_i * velocity),
                );
                assert_eq!(sample.assessment.governing_layer(), 1);
            }
            assert!(result.reaches_limit());
        }
    });
}

#[test]
fn g1_viscous_timoshenko_strength_uses_original_curvature_not_shear_rate() {
    with_cx(false, |cx| {
        let (young, shear, correction, eta_bending) = (2e9, 8e8, 0.85, 5e4);
        let source = rate_section(
            "shear-rate",
            &[(young, 1e3)],
            &[(0.0, eta_bending)],
            Some((shear, 2e4)),
        );
        let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
        let transverse = source
            .with_shear_modulus("shear_modulus", correction)
            .unwrap()
            .with_transverse_viscosity("shear_eta")
            .unwrap();
        let mut frame = Frame2d::with_mass_formulation(
            Geometry::Linear,
            FrameMassFormulation::ConsistentReference,
        );
        let length = 0.1;
        frame.add_node(0.0, 0.0);
        frame.add_node(length, 0.0);
        frame.fix(0, [true; 3]);
        frame.fix(1, [true, false, true]);
        frame
            .add_material_viscous_element(
                0,
                1,
                5,
                FrameViscousSection::try_new_with_shear(transverse, "axial_eta", "bending_eta")
                    .unwrap(),
            )
            .unwrap();
        frame
            .newmark_with_load(
                &[0.0; 3],
                2e-5,
                (0.0, 0.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
                &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            )
            .unwrap();
        let area = std::f64::consts::PI * RADIUS.powi(2);
        let inertia = area * RADIUS.powi(2) / 4.0;
        let beta =
            1.0 / (1.0 + 12.0 * young * inertia / (correction * shear * area * length * length));
        assert!(beta < 0.8);
        let stress_gradient = 6.0 * beta / (length * length)
            * (young * frame.displacements()[4] + eta_bending * frame.velocities()[4]);
        assert!(stress_gradient > 0.0);
        let result = frame.assess_tensile_strength(&capacities, 2, cx).unwrap();
        for sample in result.samples() {
            near(
                sample.assessment.layers()[0].maximum_axial_stress_pa,
                RADIUS * stress_gradient,
            );
            near(
                sample.assessment.layers()[0].minimum_axial_stress_pa,
                -RADIUS * stress_gradient,
            );
            let direction = if sample.end == 0 { 1.0 } else { -1.0 };
            near(
                sample.bending_moment_n_m,
                direction * inertia * stress_gradient,
            );
        }
    });
}

#[test]
fn g4_exact_budget_sources_and_cancellation_are_atomic() {
    let source = section("original", &[(2e9, 1e5)]);
    let capacities = [source.axial_tensile_capacity(LIMIT).unwrap()];
    let foreign = [section("foreign", &[(2e9, 1e5)])
        .axial_tensile_capacity(LIMIT)
        .unwrap()];
    let mut frame = cantilever(&[source], 0.0, Geometry::Linear);
    frame
        .static_load(&[0.0, 0.0, 0.0, -10.0, 0.0, 0.0], 1)
        .unwrap();
    let before = frame.displacements().to_vec();
    with_cx(false, |cx| {
        let original = frame.assess_tensile_strength(&capacities, 2, cx).unwrap();
        assert!(!original.reaches_limit());
        assert_eq!(original.governing_sample().unwrap().end, 0);
        assert_eq!(
            original
                .governing_sample()
                .unwrap()
                .assessment
                .utilization(),
            0.0
        );
        assert_eq!(
            frame.assess_tensile_strength(&capacities, 1, cx),
            Err(FrameTensileError::Budget)
        );
        assert_eq!(
            frame.assess_tensile_strength(&foreign, 2, cx),
            Err(FrameTensileError::SourceMismatch { element: 0 })
        );
        assert!(matches!(
            frame.assess_tensile_strength(&[], 2, cx),
            Err(FrameTensileError::Invalid(_))
        ));
        assert_eq!(
            frame.assess_tensile_strength(&capacities, 2, cx).unwrap(),
            original
        );
        let mut raw = Frame2d::new(Geometry::Linear);
        raw.add_node(0.0, 0.0);
        raw.add_node(1.0, 0.0);
        raw.add_element(0, 1, 5, 1.0, &|| SectionModel::Elastic { ea: 1.0, ei: 1.0 })
            .unwrap();
        assert_eq!(
            raw.assess_tensile_strength(&capacities, 2, cx),
            Err(FrameTensileError::SourceMismatch { element: 0 })
        );
    });
    with_cx(true, |cx| {
        assert_eq!(
            frame.assess_tensile_strength(&capacities, 2, cx),
            Err(FrameTensileError::Cancelled)
        )
    });
    assert_eq!(frame.displacements(), before);
}
