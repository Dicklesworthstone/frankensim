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
use fs_qty::{Density, Pressure};
use fs_solid::frame2d::strength::FrameTensileError;
use fs_solid::frame2d::{FrameMassFormulation, SectionModel};
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
    let mut layers = Vec::new();
    for (index, &(young, strength)) in data.iter().enumerate() {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [
            ("density", 1000.0, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            (LIMIT, strength, Pressure::DIMS),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
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
                    ScalarAdmissibility::StrictlyPositive,
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
        let material = resolve_material_state_point(
            &card,
            &QueryPoint::new(),
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
