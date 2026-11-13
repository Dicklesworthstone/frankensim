//! G1/G3/G4: native free vibration, original inertia and atomic fatigue.
//! All material and fatigue numbers are synthetic implementation fixtures.
#![cfg(feature = "rod-dynamics")]

use std::sync::Arc;

use fs_couple::rod_motion::fatigue::{RodMotionFatigueError, RodMotionFatiguePath};
use fs_couple::rod_motion::{
    PreparedRodLoads, RodMotion, RodMotionControls, RodNodalWrench, RodPointForce, SpatialRodMotion,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fatigue::{
    FatigueScope, MeanStressCorrection, NormalStressDirection, RainflowLimits, RainflowMiner,
    SnCurve, SnPoint,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, LayerIsotropicParameters,
};
use fs_material::state_point::{
    IsotropicElasticParameters, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::fatigue::{
    RodFatigueConfig, RodFatigueError, RodFatigueGauge, RodFatigueObservation, RodFatiguePoint,
};
use fs_solid::rod::material_field::motion::RodRootBoundary;
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};

fn rod(supported: bool) -> SpatialMaterialRod {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let domain = ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 300.0);
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e4, Pressure::DIMS),
        ("shear_modulus", 4e3, Pressure::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic motion-fatigue test".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic rod".into(),
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
        &QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.05,
            material,
        }],
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.8)
    .unwrap();
    let rod = SpatialMaterialRod::try_straight(0.8, vec![section]).unwrap();
    if supported {
        rod.with_supports(
            &[RodSupport {
                node: 1,
                position_m: [None, Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &CancelGate::new(),
        )
        .unwrap()
    } else {
        rod
    }
}

fn gauge(correction: MeanStressCorrection) -> RodFatigueGauge {
    let curve = Arc::new(
        SnCurve::new(
            FatigueScope {
                source_reference: "synthetic implementation calibration".into(),
                // An S-N scope names the material "chemistry; phase".
                material: "synthetic rod; solid".into(),
                process: "test".into(),
                surface_and_geometry: "synthetic shaft".into(),
                environment: "synthetic".into(),
                load_basis: "fixed-plane normal stress amplitude".into(),
                temperature_k: 300.0,
            },
            vec![
                SnPoint {
                    amplitude_pa: 1e-8,
                    cycles_to_failure: 1e16,
                },
                SnPoint {
                    amplitude_pa: 1e4,
                    cycles_to_failure: 1.0,
                },
            ],
            Some(1e-9),
        )
        .unwrap(),
    );
    RodFatigueGauge {
        point: RodFatiguePoint {
            segment: 0,
            layer: 0,
            coordinates_yz_m: [0.05, 0.0],
        },
        normal_in_material_frame: NormalStressDirection::new([1.0, 0.0, 0.0]).unwrap(),
        temperature_axis: "T".into(),
        curve,
        mean_stress: correction,
    }
}
fn config(samples: u64, gauges: usize) -> RodFatigueConfig {
    RodFatigueConfig {
        history: RainflowLimits {
            maximum_samples: samples,
            maximum_residue: 64,
        },
        maximum_gauges: gauges,
        maximum_total_residue: gauges * 64,
        maximum_total_findley_samples: 0,
        maximum_absolute_axial_strain: 0.05,
        maximum_engineering_shear: 0.05,
        maximum_unassessed_transverse_force_n: 0.0,
    }
}
fn controls(dt: f64) -> RodMotionControls {
    RodMotionControls {
        maximum_step_s: dt,
        maximum_steps: 1000,
        energy_defect_limit_j: 1e-4,
        maximum_rotation_rad: 0.25,
        balance_origin_m: [0.0; 3],
    }
}
fn initial(outgoing_velocity: f64) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        rod(true)
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.008;
    RodMotion::new(
        model,
        pose,
        &[[0.0; 3], [outgoing_velocity, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap()
}
fn same_motion(a: &RodMotion, b: &RodMotion) {
    assert_eq!(a.pose(), b.pose());
    assert_eq!(a.linear_momenta_kg_m_s(), b.linear_momenta_kg_m_s());
    assert_eq!(a.body_spins_kg_m2_s(), b.body_spins_kg_m2_s());
    assert_eq!(a.time_s(), b.time_s());
    assert_eq!(a.accepted_steps(), b.accepted_steps());
    assert_eq!(a.observe([0.0; 3]).unwrap(), b.observe([0.0; 3]).unwrap());
}

#[test]
fn transverse_shear_recovery_follows_actual_motion_without_changing_native_work() {
    let gate = CancelGate::new();
    let original = rod(false);
    let layers = original.sections()[0].elastic().layers().to_vec();
    let requirement = |name| {
        layers[0]
            .material
            .property(name)
            .unwrap()
            .requirement()
            .clone()
    };
    let binding = LayerIsotropicParameters {
        density: requirement("density"),
        parameters: IsotropicElasticParameters::YoungShear {
            young_modulus: requirement("young_modulus"),
            shear_modulus: requirement("shear_modulus"),
        },
    };
    let section =
        ConcentricElasticSection::try_new_with_isotropic_parameters(0.0, layers, &[binding])
            .unwrap()
            .with_isotropic_shear(0.8)
            .unwrap();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section])
            .unwrap()
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][1] = 0.0008;
    let mut native =
        RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, &gate).unwrap();
    let mut point = gauge(MeanStressCorrection::Uncorrected);
    point.point.coordinates_yz_m = [0.0; 2];
    point.normal_in_material_frame = NormalStressDirection::new([1.0, 1.0, 0.0]).unwrap();
    let limits = config(17, 1);
    assert!(
        RodMotionFatiguePath::new(native.clone(), vec![point.clone()], vec![], limits, &gate)
            .is_err()
    );
    let mut path = RodMotionFatiguePath::new_with_transverse_shear(
        native.clone(),
        vec![point],
        vec![],
        limits,
        &gate,
    )
    .unwrap();
    let first = path.observation().normal_samples()[0];
    // Native gamma_y=0.001, kappa=0.8, G=4000 Pa. The exact circular
    // centerline multiplier at nu=1/4 is 7/5, independently of that kappa.
    assert!((first.stress_pa[3] - (7.0 / 5.0) * 0.8 * 4000.0 * 0.001).abs() < 1e-12);
    assert_eq!(first.unassessed_transverse_force_n, [0.0; 2]);
    for _ in 0..16 {
        let (next, expected) = native
            .step_candidate(
                0.001,
                &[RodNodalWrench::default(); 2],
                controls(0.001),
                &gate,
            )
            .unwrap();
        let report = path
            .advance(
                0.001,
                &[RodNodalWrench::default(); 2],
                controls(0.001),
                &gate,
            )
            .unwrap();
        native = next;
        same_motion(path.motion(), &native);
        assert_eq!(report.after, expected.after);
        assert_eq!(report.integration_defect_j, expected.integration_defect_j);
        assert_eq!(report.support_impulses, expected.support_impulses);
        assert_eq!(
            path.observation().normal_samples()[0].unassessed_transverse_force_n,
            [0.0; 2]
        );
    }
    assert_ne!(
        path.observation().normal_samples()[0].stress_pa[3],
        first.stress_pa[3]
    );
    let before = path.clone();
    assert!(
        path.advance(
            0.001,
            &[RodNodalWrench::default(); 2],
            controls(0.001),
            &gate
        )
        .is_err()
    );
    same_motion(path.motion(), before.motion());
    assert_eq!(
        path.observation().normal_samples(),
        before.observation().normal_samples()
    );
    assert_eq!(path.summaries(), before.summaries());
}

#[test]
fn unforced_vibration_counts_real_stress_and_preserves_native_motion() {
    let gate = CancelGate::new();
    let mut native = initial(0.0);
    let spec = gauge(MeanStressCorrection::Uncorrected);
    let limits = config(361, 1);
    let mut independent =
        RainflowMiner::new(spec.curve.clone(), spec.mean_stress, limits.history).unwrap();
    let mut path =
        RodMotionFatiguePath::new(native.clone(), vec![spec], vec![], limits, &gate).unwrap();
    let stress = |state: &RodMotion| 1e4 * (state.pose().positions_m[1][0] / 0.8 - 1.0);
    independent.push(stress(&native)).unwrap();
    assert_eq!(path.sample_count(), 1);
    assert!((path.observation().normal_samples()[0].stress_pa[0] - 100.0).abs() < 1e-10);
    let mut minimum = 100.0_f64;
    for step in 0..360 {
        let (next, expected) = native
            .step_candidate(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .unwrap();
        let report = path
            .advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .unwrap();
        native = next;
        same_motion(path.motion(), &native);
        assert_eq!(report.after, expected.after);
        assert_eq!(report.integration_defect_j, expected.integration_defect_j);
        assert_eq!(report.support_impulses, expected.support_impulses);
        assert_eq!(
            path.observation().normal_samples()[0].stress_pa[0],
            stress(&native)
        );
        minimum = minimum.min(stress(&native));
        independent.push(stress(&native)).unwrap();
        if step == 130 {
            let before = path.summaries();
            assert_eq!(
                path.preview_fatigue(&gate).unwrap()[0],
                independent.preview_finish().unwrap()
            );
            assert_eq!(path.summaries(), before);
        }
    }
    assert!(
        minimum < -99.0,
        "unforced motion must actually reverse stress"
    );
    assert_eq!(path.sample_count(), 361);
    let expected = independent.finish().unwrap();
    assert!(expected.miner_damage > 0.0);
    assert_eq!(path.finish(&gate).unwrap(), vec![expected]);
    assert!(
        path.advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .is_err()
    );
    same_motion(path.motion(), &native);
}

#[test]
fn late_second_gauge_refusal_rolls_back_motion_and_first_counter() {
    let gate = CancelGate::new();
    let mut path = RodMotionFatiguePath::new(
        initial(0.1),
        vec![
            gauge(MeanStressCorrection::Uncorrected),
            gauge(MeanStressCorrection::GoodmanTensionOnly {
                ultimate_strength_pa: 105.0,
            }),
        ],
        vec![],
        config(10, 2),
        &gate,
    )
    .unwrap();
    let old = path.clone();
    let error = path
        .advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
        .unwrap_err();
    assert!(matches!(
        error,
        RodMotionFatigueError::Fatigue(RodFatigueError::Gauge { gauge: 1, .. })
    ));
    same_motion(path.motion(), old.motion());
    assert_eq!(path.summaries(), old.summaries());
    assert_eq!(
        path.observation().normal_samples(),
        old.observation().normal_samples()
    );
    assert_eq!(path.sample_count(), 1);
    // A smaller physical step is still admissible from that exact prefix.
    path.advance(
        0.001,
        &[RodNodalWrench::default(); 2],
        controls(0.001),
        &gate,
    )
    .unwrap();
    assert_eq!(path.sample_count(), 2);
    assert_eq!(path.motion().accepted_steps(), 1);
}

#[test]
fn budget_and_cancellation_preserve_cloned_continuation() {
    let gate = CancelGate::new();
    let mut path = RodMotionFatiguePath::new(
        initial(0.0),
        vec![gauge(MeanStressCorrection::Uncorrected)],
        vec![],
        config(2, 1),
        &gate,
    )
    .unwrap();
    let mut other = path.clone();
    path.advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
        .unwrap();
    other
        .advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
        .unwrap();
    same_motion(path.motion(), other.motion());
    assert_eq!(path.summaries(), other.summaries());
    assert!(
        path.advance(0.01, &[RodNodalWrench::default(); 2], controls(0.01), &gate)
            .is_err()
    );
    same_motion(path.motion(), other.motion());
    assert_eq!(
        path.observation().normal_samples(),
        other.observation().normal_samples()
    );
    assert_eq!(path.sample_count(), 2);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        path.advance(
            0.01,
            &[RodNodalWrench::default(); 2],
            controls(0.01),
            &cancelled
        )
        .is_err()
    );
    same_motion(path.motion(), other.motion());
    assert_eq!(
        path.preview_fatigue(&gate).unwrap(),
        other.preview_fatigue(&gate).unwrap()
    );
}

#[test]
fn material_frame_follows_rigid_pose_and_retains_original_source_owner() {
    let gate = CancelGate::new();
    let model = rod(false)
        .into_motion_model(RodRootBoundary::Free, &gate)
        .unwrap();
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.008;
    let mut limits = config(8, 1);
    // Explicit numerical allowance on native transverse resultant, not Q/A.
    limits.maximum_unassessed_transverse_force_n = 1e-10;
    let mut bank = RodFatigueObservation::new_motion(
        &model,
        &pose,
        vec![gauge(MeanStressCorrection::Uncorrected)],
        vec![],
        limits,
        &gate,
    )
    .unwrap();
    for angle in [0.3, -0.7, 1.2] {
        let rotation =
            So3::exp(So3Tangent::new(Vec3::new(angle, 0.4 * angle, -0.3 * angle))).unwrap();
        let mut transformed = pose.clone();
        for (position, frame) in transformed
            .positions_m
            .iter_mut()
            .zip(&mut transformed.frames)
        {
            let p = rotation
                .rotate(Vec3::new(position[0], position[1], position[2]))
                .unwrap();
            *position = [p.x + 0.1, p.y - 0.2, p.z + 0.3];
            *frame = rotation;
        }
        bank.observe_motion(&model.clone(), &transformed, &gate)
            .unwrap();
        assert!((bank.normal_samples()[0].stress_pa[0] - 100.0).abs() < 1e-9);
    }
    let old = bank.normal_samples().to_vec();
    let other = rod(false)
        .into_motion_model(RodRootBoundary::Free, &gate)
        .unwrap();
    assert!(bank.observe_motion(&other, &pose, &gate).is_err());
    assert_eq!(bank.normal_samples(), old);
    assert_eq!(bank.sample_count(), 4);
    assert_eq!(bank.finish(&gate).unwrap()[0].miner_damage, 0.0);
}

#[test]
fn eccentric_reference_inertia_and_spatial_force_use_the_same_native_step() {
    let gate = CancelGate::new();
    let source = rod(true);
    let link = MaterialTensionLink::try_new(
        source.sections()[0].elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: [0.0, 0.01, 0.0],
            },
            TensionLinkAttachment {
                node: 1,
                offset_body_m: [0.0, 0.01, 0.0],
            },
        ],
        0.81,
        1,
        0.05,
    )
    .unwrap();
    let model = Arc::new(
        source
            .with_tension_links(vec![link], &gate)
            .unwrap()
            .into_spatial_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    assert!(model.inertia()[1].center_of_mass_body_m()[1] > 0.0);
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.004;
    let native = SpatialRodMotion::new_spatial(
        model,
        pose,
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let loads = PreparedRodLoads::new(
        &[RodNodalWrench::default(); 2],
        &[[0.0; 3]; 2],
        &[RodPointForce {
            node: 1,
            offset_body_m: [0.0, 0.02, 0.0],
            force_world_n: [0.1, 0.0, 0.0],
        }],
        &gate,
    )
    .unwrap();
    let (expected, report) = native
        .step_spatial_candidate(0.01, &loads, controls(0.01), &gate)
        .unwrap();
    let mut path = RodMotionFatiguePath::new(
        native,
        vec![gauge(MeanStressCorrection::Uncorrected)],
        vec![],
        config(4, 1),
        &gate,
    )
    .unwrap();
    let actual = path
        .advance_spatial(0.01, &loads, controls(0.01), &gate)
        .unwrap();
    assert_eq!(path.motion().pose(), expected.pose());
    assert_eq!(
        path.motion().linear_momenta_kg_m_s(),
        expected.linear_momenta_kg_m_s()
    );
    assert_eq!(
        path.motion().body_spins_kg_m2_s(),
        expected.body_spins_kg_m2_s()
    );
    assert_eq!(actual.after, report.after);
    assert_eq!(actual.support_impulses, report.support_impulses);
    assert_eq!(actual.integration_defect_j, report.integration_defect_j);
    assert_eq!(path.sample_count(), 2);
}
