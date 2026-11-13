//! Native frame motion and accepted metal state drive bounded empirical fatigue.
//! Every calibration here is synthetic; these are implementation checks, not service-life validation.
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::elastic::IsotropicElastic;
use fs_material::fatigue::strain_life::{
    CmbParameters, PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainStressSample,
};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::fiber::MenegottoPintoSteel;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::reduction::{UniaxialJ2, UniaxialVoce};
use fs_material::plastic::voce::VoceJ2;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::fatigue::{
    FrameFatigueConfig, FrameFatigueError, FrameFatiguePath, FrameFiberFatigueGauge,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

#[path = "frame_fatigue/shear.rs"]
mod shear;

// Material-defined annuli retain separate elastic stresses and original S-N sources.
mod elastic_layers {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_material::fatigue::critical_plane::findley::FindleyCalibration;
    use fs_material::fatigue::critical_plane::findley::source::{
        FindleyProperty, FindleySourceRequest, resolve_findley,
    };
    use fs_material::fatigue::source::{
        SnLifeMeasure, SnSourceRequest, SnStressMeasure, resolve_log_log_sn,
    };
    use fs_material::fatigue::{MeanStressCorrection, SnCurve};
    use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
    use fs_material::state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    };
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Density, Dims, Pressure, QuantitySpec};
    use fs_solid::frame2d::fatigue::FrameElasticFatigueGauge;

    fn source(young: f64, revision: u32) -> (ResolvedMaterialStatePoint, Arc<SnCurve>) {
        let (material, curve, _) = source_bundle(young, revision);
        (material, curve)
    }
    pub(super) fn source_bundle(
        young: f64,
        revision: u32,
    ) -> (
        ResolvedMaterialStatePoint,
        Arc<SnCurve>,
        Arc<FindleyCalibration>,
    ) {
        let request = SnSourceRequest {
            stress_measure: SnStressMeasure::Amplitude,
            life_measure: SnLifeMeasure::Cycles,
            temperature_axis: "absolute-temperature".into(),
            stress_ratio_axis: "calibration-R".into(),
            surface_and_geometry: "synthetic layer coupon".into(),
            environment: "synthetic fixed condition".into(),
            load_basis: "fully reversed normal stress amplitude".into(),
            maximum_points: 5,
        };
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let dimensionless = QuantitySpec::dimensional(Dims::NONE);
        let domain = ValidityDomain::unconstrained().with_quantity(
            &request.temperature_axis,
            temperature,
            300.0,
            300.0,
        );
        let point = QueryPoint::new()
            .with_quantity(&request.temperature_axis, temperature, 300.0)
            .unwrap();
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, 7800.0),
            ("young_modulus", Pressure::DIMS, young),
            ("shear_modulus", Pressure::DIMS, young / 2.6),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: domain.clone(),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic frame layer".into(),
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
        let life = request.life_measure.axis_name();
        claims
            .insert_claim(PropertyClaim {
                key: request.stress_measure.property_key(),
                value: PropertyValue::Curve {
                    abscissa: life.into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![
                        (1.0, 160e6),
                        (10.0, 80e6),
                        (100.0, 40e6),
                        (1000.0, 20e6),
                        (10000.0, 10e6),
                    ],
                    dims: Pressure::DIMS,
                },
                validity: domain
                    .clone()
                    .with_quantity(&request.stress_ratio_axis, dimensionless, -1.0, -1.0)
                    .with_quantity(life, dimensionless, 1.0, 10000.0),
                interpolation: InterpolationPolicy::TabulatedOnly,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic S-N layer coupon".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        for (property, dims, value) in [
            (FindleyProperty::NormalStressSensitivity, Dims::NONE, 0.2),
            (FindleyProperty::LimitFactor, Pressure::DIMS, 10e6),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: property.property_key(),
                    value: PropertyValue::Scalar { value, dims },
                    validity: domain.clone(),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic common Findley fit".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        // A lineage starts at revision 0; each later revision supersedes the one before.
        let mut card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic layered alloy".into(),
                phase: "solid".into(),
                process: "synthetic coupon".into(),
                revision: 0,
            },
            claims.clone(),
            vec![],
        )
        .unwrap();
        for _ in 0..revision {
            card = MaterialCard::supersede(&card, claims.clone(), vec![]).unwrap();
        }
        let fit = point
            .clone()
            .with_quantity(&request.stress_ratio_axis, dimensionless, -1.0)
            .unwrap()
            .with_quantity(life, dimensionless, 100.0)
            .unwrap();
        let curve = Arc::new(
            resolve_log_log_sn(
                &card,
                &fit,
                &request,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap(),
        );
        let material = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let findley = resolve_findley(
            &card,
            &point,
            &FindleySourceRequest {
                temperature_axis: request.temperature_axis,
                surface_and_geometry: "synthetic circular coupon".into(),
                environment: "synthetic fixed condition".into(),
                load_basis: "synthetic common normal/transverse-shear calibration".into(),
            },
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        )
        .unwrap();
        (material, curve, Arc::new(findley))
    }
    fn native(youngs: [f64; 2], bending: bool) -> (Frame2d, Vec<FrameElasticFatigueGauge>, f64) {
        let mut layers = Vec::new();
        let mut gauges = Vec::new();
        for (index, young) in youngs.into_iter().enumerate() {
            let (material, curve) = source(young, index as u32);
            let radius = [0.01, 0.02][index];
            layers.push(ConcentricMaterialLayer {
                outer_radius_m: radius,
                material,
            });
            gauges.push(FrameElasticFatigueGauge {
                element: 0,
                section: 0,
                layer: index,
                cross_section_position_m: [radius, 0.0],
                temperature_axis: "absolute-temperature".into(),
                curve,
                mean_stress: MeanStressCorrection::Uncorrected,
            });
        }
        let section = ConcentricElasticSection::try_new(0.0, layers).unwrap();
        let ea = std::f64::consts::PI
            * (youngs[0] * 0.01_f64.powi(2) + youngs[1] * (0.02_f64.powi(2) - 0.01_f64.powi(2)));
        let mut frame = Frame2d::new(Geometry::Linear);
        let left = frame.add_node(0.0, 0.0);
        let right = frame.add_node(1.0, 0.0);
        frame.fix(left, [true; 3]);
        frame.fix(
            right,
            if bending {
                [true, false, false]
            } else {
                [false, true, true]
            },
        );
        frame.add_material_element(left, right, 3, section).unwrap();
        (frame, gauges, ea)
    }

    #[test]
    fn original_layer_moduli_drive_distinct_stress_and_hand_counted_damage() {
        with_cx(&CancelGate::new(), |cx| {
            let (frame, gauges, ea) = native([E, E / 2.0], false);
            let identities: Vec<_> = gauges
                .iter()
                .map(|g| g.curve.source().unwrap().card_identity())
                .collect();
            let mut run =
                FrameFatiguePath::new_with_elastic(frame, vec![], gauges, config(), cx).unwrap();
            let steps = run
                .displacement_control(&PATTERN, 3, &[0.0004, -0.0004, 0.0004], cx)
                .unwrap();
            close(steps.last().unwrap().lambda, ea * 0.0004, 1e-9);
            for (i, expected) in [80e6, 40e6].into_iter().enumerate() {
                close(run.elastic_samples()[i].normal_stress_pa, expected, 1e-10);
                close(run.elastic_samples()[i].mechanical_strain, 0.0004, 1e-10);
                assert_eq!(
                    run.elastic_material_source(i).unwrap().card_identity(),
                    identities[i]
                );
            }
            let before = run.summaries();
            let preview = run.preview_finish(cx).unwrap();
            assert_eq!(run.summaries(), before);
            assert!(!run.is_finished());
            assert!(preview.fibers.is_empty());
            for (summary, damage) in preview.elastic.iter().zip([0.105, 0.0105]) {
                close(summary.miner_damage, damage, 2e-9);
                assert_eq!(summary.half_cycles, 3);
            }
            assert_eq!(preview, run.finish(cx).unwrap());
        });
    }

    #[test]
    fn opposite_layer_surfaces_retain_signed_bending_stress() {
        with_cx(&CancelGate::new(), |cx| {
            let (frame, mut gauges, _) = native([E, E / 2.0], true);
            gauges.remove(0);
            let mut opposite = gauges[0].clone();
            opposite.cross_section_position_m[0] *= -1.0;
            gauges.push(opposite);
            let mut run =
                FrameFatiguePath::new_with_elastic(frame, vec![], gauges, config(), cx).unwrap();
            run.displacement_control(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 5, &[0.02, -0.02, 0.02], cx)
                .unwrap();
            let points = run.elastic_samples();
            close(points[0].mechanical_strain, -0.0004, 2e-9);
            close(points[0].normal_stress_pa, -40e6, 2e-9);
            close(points[1].normal_stress_pa, 40e6, 2e-9);
            let done = run.finish(cx).unwrap();
            close(done.elastic[0].miner_damage, 0.0105, 2e-9);
            assert_eq!(done.elastic[0], done.elastic[1]);
        });
    }

    #[test]
    fn layer_geometry_temperature_and_original_source_are_not_interchangeable() {
        with_cx(&CancelGate::new(), |cx| {
            let (frame, gauges, _) = native([E, E / 2.0], false);
            for position in [[0.015, 0.0], [f64::NAN, 0.0], [0.0, f64::INFINITY]] {
                let mut wrong = gauges[0].clone();
                wrong.cross_section_position_m = position;
                assert!(
                    FrameFatiguePath::new_with_elastic(
                        frame.clone(),
                        vec![],
                        vec![wrong],
                        config(),
                        cx
                    )
                    .is_err()
                );
            }
            let mut wrong = gauges[1].clone();
            wrong.cross_section_position_m = [0.005, 0.0];
            assert!(
                FrameFatiguePath::new_with_elastic(
                    frame.clone(),
                    vec![],
                    vec![wrong],
                    config(),
                    cx
                )
                .is_err()
            );
            let mut wrong = gauges[0].clone();
            wrong.temperature_axis = "T".into();
            assert!(matches!(
                FrameFatiguePath::new_with_elastic(
                    frame.clone(),
                    vec![],
                    vec![wrong],
                    config(),
                    cx
                ),
                Err(FrameFatigueError::SourceMismatch { .. })
            ));
            let mut wrong = gauges[0].clone();
            wrong.curve = source(E, 8).1;
            assert!(matches!(
                FrameFatiguePath::new_with_elastic(frame, vec![], vec![wrong], config(), cx),
                Err(FrameFatigueError::SourceMismatch { .. })
            ));
        });
    }

    #[test]
    fn late_layer_life_refusal_rolls_back_all_accepted_native_steps() {
        with_cx(&CancelGate::new(), |cx| {
            let (frame, gauges, _) = native([E / 2.0, E], false);
            let mut run =
                FrameFatiguePath::new_with_elastic(frame, vec![], gauges, config(), cx).unwrap();
            let before = run.clone();
            assert!(matches!(
                run.displacement_control(&PATTERN, 3, &[0.0012, -0.0012, 0.0012, -0.0012], cx),
                Err(FrameFatigueError::Gauge { gauge: 1, .. })
            ));
            assert_eq!(run.frame().displacements(), before.frame().displacements());
            assert_eq!(run.sample_count(), before.sample_count());
            assert_eq!(run.summaries(), before.summaries());
            assert_eq!(run.elastic_samples(), before.elastic_samples());
            assert_eq!(
                run.preview_finish(cx).unwrap(),
                before.preview_finish(cx).unwrap()
            );
            let mut uninterrupted = before;
            let good = [0.0004, -0.0004, 0.0004];
            run.displacement_control(&PATTERN, 3, &good, cx).unwrap();
            uninterrupted
                .displacement_control(&PATTERN, 3, &good, cx)
                .unwrap();
            assert_eq!(run.finish(cx).unwrap(), uninterrupted.finish(cx).unwrap());
        });
    }
}

const E: f64 = 200e9;
const SY: f64 = 200e6;
const HI: f64 = 2e9;
const HK: f64 = 18e9;
const AREA: f64 = 1e-3;
const ARM: f64 = 0.1;
const PATTERN: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn backbone() -> J2Plasticity {
    J2Plasticity::new_mixed(IsotropicElastic::new(E, 0.3, 0.04).unwrap(), SY, HI, HK).unwrap()
}
fn law() -> FiberLaw {
    FiberLaw::J2(UniaxialJ2::new(backbone()).unwrap())
}
fn branches<const N: usize>() -> FiberLaw {
    let voce = VoceJ2::new(backbone(), 80e6, 500.0).unwrap();
    let native = ChabocheJ2::<N>::new(
        voce,
        std::array::from_fn(|i| {
            BackstressBranch::new(6e9 + 2e9 * i as f64, 50.0 * (i + 1).pow(2) as f64).unwrap()
        }),
    )
    .unwrap();
    FiberLaw::chaboche_branches(UniaxialChaboche::new(native).unwrap()).unwrap()
}
fn frame(law: &FiberLaw, bending: bool) -> Frame2d {
    let mut f = Frame2d::new(Geometry::Linear);
    let left = f.add_node(0.0, 0.0);
    let right = f.add_node(1.0, 0.0);
    f.fix(left, [true; 3]);
    f.fix(
        right,
        if bending {
            [true, false, false]
        } else {
            [false, true, true]
        },
    );
    f.add_mass(right, 100.0);
    f.add_element(left, right, 3, 0.0, &|| {
        SectionModel::Fiber(Section::new(
            [-ARM, ARM]
                .map(|y| Fiber {
                    y,
                    area: AREA,
                    law: law.clone(),
                })
                .to_vec(),
        ))
    })
    .unwrap();
    f
}
fn curve(minimum_cycles: f64) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic axial CMB".into(),
                material: "synthetic metal".into(),
                process: "synthetic fixed condition".into(),
                surface_and_geometry: "axial coupon".into(),
                environment: "synthetic dry atmosphere".into(),
                load_basis: "total mechanical strain amplitude".into(),
                temperature_k: 300.0,
            },
            CmbParameters {
                youngs_modulus_pa: E,
                fatigue_strength_coefficient_pa: 1e9,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.005,
                fatigue_ductility_exponent: -0.5,
            },
            minimum_cycles,
            1e12,
        )
        .unwrap(),
    )
}
fn config() -> FrameFatigueConfig {
    FrameFatigueConfig {
        history: RainflowLimits {
            maximum_samples: 1000,
            maximum_residue: 32,
        },
        maximum_gauges: 8,
        maximum_total_residue: 256,
        maximum_absolute_axial_strain: 0.04,
    }
}
fn gauge(fiber: usize, correction: StrainLifeMeanStress) -> FrameFiberFatigueGauge {
    FrameFiberFatigueGauge {
        element: 0,
        section: 0,
        fiber,
        material_source: None,
        curve: curve(1.0),
        mean_stress: correction,
    }
}
fn path(f: Frame2d, correction: StrainLifeMeanStress, cx: &Cx<'_>) -> FrameFatiguePath {
    FrameFatiguePath::new(f, vec![gauge(0, correction)], config(), cx).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual:.16e} != {expected:.16e}"
    );
}

// Independent scalar associative return: no invocation of the material owner.
fn prager(strain: f64, plastic: &mut f64, accumulated: &mut f64) -> f64 {
    let trial = E * (strain - *plastic) - HK * *plastic;
    let excess = trial.abs() - (SY + HI * *accumulated);
    if excess > 0.0 {
        let increment = excess / (E + HI + HK);
        *plastic += trial.signum() * increment;
        *accumulated += increment;
    }
    E * (strain - *plastic)
}

#[test]
fn native_axial_cycling_matches_independent_prager_and_paired_fatigue() {
    with_cx(&CancelGate::new(), |cx| {
        let mut run = path(frame(&law(), false), StrainLifeMeanStress::Morrow, cx);
        let mut reference =
            PairedStrainLifeMiner::new(curve(1.0), StrainLifeMeanStress::Morrow, config().history)
                .unwrap();
        reference
            .push(StrainStressSample {
                mechanical_strain: 0.0,
                normal_stress_pa: 0.0,
            })
            .unwrap();
        let (mut plastic, mut accumulated) = (0.0, 0.0);
        for strain in [0.003, -0.003, 0.004, -0.002, 0.003] {
            run.displacement_control(&PATTERN, 3, &[strain], cx)
                .unwrap();
            let stress = prager(strain, &mut plastic, &mut accumulated);
            let actual = run.samples()[0];
            close(actual.mechanical_strain, strain, 1e-12);
            close(actual.normal_stress_pa, stress, 2e-9);
            close(run.frame().element_forces(0)[0], 2.0 * AREA * stress, 2e-9);
            assert!((actual.normal_stress_pa - E * strain).abs() > 10e6);
            reference
                .push(StrainStressSample {
                    mechanical_strain: strain,
                    normal_stress_pa: stress,
                })
                .unwrap();
        }
        let live = run.summaries();
        let preview = run.preview_finish(cx).unwrap();
        assert_eq!(run.summaries(), live);
        assert!(!run.is_finished());
        close(
            preview.fibers[0].miner_damage,
            reference.finish().unwrap().miner_damage,
            2e-9,
        );
        assert_eq!(run.finish(cx).unwrap(), preview);
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.0], cx),
            Err(FrameFatigueError::Finished)
        ));
    });
}

#[test]
fn hand_counted_cycle_retains_signed_mean_stress_and_bending_fibers() {
    with_cx(&CancelGate::new(), |cx| {
        for sign in [-1.0, 1.0] {
            let amplitude = 0.0005 * sign;
            let mut run = FrameFatiguePath::new(
                frame(&law(), false),
                vec![
                    gauge(0, StrainLifeMeanStress::Uncorrected),
                    gauge(0, StrainLifeMeanStress::Morrow),
                ],
                config(),
                cx,
            )
            .unwrap();
            run.displacement_control(&PATTERN, 3, &[amplitude, 0.0], cx)
                .unwrap();
            let result = run.finish(cx).unwrap();
            // Two residual half cycles form one cycle with epsilon_a=|amplitude|/2.
            let epsilon_a = amplitude.abs() * 0.5;
            let mean_stress = E * amplitude * 0.5;
            let ordinary_life = 0.5 * (0.01 / epsilon_a).powi(2);
            let morrow_life = 0.5 * ((0.01 - mean_stress / E) / epsilon_a).powi(2);
            close(result.fibers[0].miner_damage, 1.0 / ordinary_life, 1e-12);
            close(result.fibers[1].miner_damage, 1.0 / morrow_life, 1e-12);
            if sign > 0.0 {
                assert!(result.fibers[1].miner_damage > result.fibers[0].miner_damage);
            } else {
                assert!(result.fibers[1].miner_damage < result.fibers[0].miner_damage);
            }
        }
        let mut run = FrameFatiguePath::new(
            frame(&law(), true),
            vec![
                gauge(0, StrainLifeMeanStress::Morrow),
                gauge(1, StrainLifeMeanStress::Morrow),
            ],
            config(),
            cx,
        )
        .unwrap();
        run.displacement_control(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 5, &[0.02], cx)
            .unwrap();
        let samples = run.samples();
        close(samples[0].mechanical_strain, 0.002, 1e-12);
        close(samples[1].mechanical_strain, -0.002, 1e-12);
        close(
            samples[0].normal_stress_pa,
            -samples[1].normal_stress_pa,
            1e-10,
        );
        assert!(samples[0].normal_stress_pa > 0.0);
        assert!(samples[1].normal_stress_pa < 0.0);
    });
}

#[test]
fn every_native_metal_variant_exposes_accepted_stress_without_an_elastic_replay() {
    with_cx(&CancelGate::new(), |cx| {
        let laws = [
            law(),
            FiberLaw::Voce(
                UniaxialVoce::new(VoceJ2::new(backbone(), 80e6, 500.0).unwrap()).unwrap(),
            ),
            branches::<1>(),
            branches::<2>(),
            branches::<3>(),
            branches::<4>(),
            branches::<5>(),
            branches::<6>(),
            branches::<7>(),
            branches::<8>(),
            FiberLaw::Steel(MenegottoPintoSteel::new(E, SY, 0.05).unwrap()),
        ];
        for law in laws {
            let mut run = path(frame(&law, false), StrainLifeMeanStress::Uncorrected, cx);
            for strain in [0.003, -0.003, 0.002, -0.002] {
                run.displacement_control(&PATTERN, 3, &[strain], cx)
                    .unwrap();
                let actual = run.samples()[0];
                close(actual.mechanical_strain, strain, 1e-11);
                close(
                    actual.normal_stress_pa,
                    run.frame().element_forces(0)[0] / (2.0 * AREA),
                    2e-8,
                );
                assert!((actual.normal_stress_pa - E * strain).abs() > 10e6);
            }
            assert!(run.finish(cx).unwrap().fibers[0].miner_damage > 0.0);
        }
    });
}

#[test]
fn native_load_ramps_static_chunks_clones_and_plateaus_preserve_history() {
    with_cx(&CancelGate::new(), |cx| {
        let mut run = path(frame(&law(), false), StrainLifeMeanStress::Uncorrected, cx);
        let load = [0.0, 0.0, 0.0, 100_000.0, 0.0, 0.0];
        let mut original = run.frame().clone();
        let expected = original.static_load(&load, 4).unwrap();
        assert_eq!(run.static_load(&load, 4, cx).unwrap(), expected);
        assert_eq!(run.sample_count(), 5);
        let plateau = run.samples()[0].mechanical_strain;
        let preview = run.preview_finish(cx).unwrap().fibers[0].miner_damage;
        run.displacement_control(&PATTERN, 3, &[plateau; 4], cx)
            .unwrap();
        close(
            run.preview_finish(cx).unwrap().fibers[0].miner_damage,
            preview,
            1e-12,
        );
        let mut split = run.clone();
        let targets = [0.003, -0.003, 0.002, -0.002, 0.003];
        run.displacement_control(&PATTERN, 3, &targets, cx).unwrap();
        split
            .displacement_control(&PATTERN, 3, &targets[..2], cx)
            .unwrap();
        let mut resumed = split.clone();
        for target in &targets[2..] {
            resumed
                .displacement_control(&PATTERN, 3, &[*target], cx)
                .unwrap();
        }
        assert_eq!(run.sample_count(), resumed.sample_count());
        close(
            run.finish(cx).unwrap().fibers[0].miner_damage,
            resumed.finish(cx).unwrap().fibers[0].miner_damage,
            1e-11,
        );
        close(
            run.frame().displacements()[3],
            resumed.frame().displacements()[3],
            1e-12,
        );
    });
}

#[test]
fn late_gauge_refusal_rolls_back_motion_plasticity_and_all_counters() {
    with_cx(&CancelGate::new(), |cx| {
        let mut narrow = gauge(1, StrainLifeMeanStress::Uncorrected);
        narrow.curve = curve(100.0);
        let mut run = FrameFatiguePath::new(
            frame(&law(), false),
            vec![gauge(0, StrainLifeMeanStress::Uncorrected), narrow],
            config(),
            cx,
        )
        .unwrap();
        let mut unchanged = run.clone();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.003, -0.003, 0.003, -0.003], cx),
            Err(FrameFatigueError::Gauge { gauge: 1, .. })
        ));
        assert_eq!(run.samples(), unchanged.samples());
        assert_eq!(run.summaries(), unchanged.summaries());
        assert_eq!(
            run.frame().displacements(),
            unchanged.frame().displacements()
        );
        assert_eq!(run.sample_count(), 1);
        // Hidden plastic history must also survive the refusal unchanged.
        for path in [&mut run, &mut unchanged] {
            path.displacement_control(&PATTERN, 3, &[0.0005, 0.0], cx)
                .unwrap();
        }
        assert_eq!(run.samples(), unchanged.samples());
        assert_eq!(run.finish(cx).unwrap(), unchanged.finish(cx).unwrap());
    });
}

#[test]
fn native_newmark_motion_and_energy_are_unchanged_and_pending_cycles_span_blocks() {
    with_cx(&CancelGate::new(), |cx| {
        let mut initial = frame(&law(), false);
        initial
            .displacement_control(&PATTERN, 3, &[0.0005])
            .unwrap();
        let mut unobserved = initial.clone();
        let mut run = path(initial, StrainLifeMeanStress::Uncorrected, cx);
        let mut split = run.clone();
        let ground = vec![0.0; 301];
        let damping = Rayleigh { a0: 0.0, a1: 0.0 };
        let original = unobserved
            .newmark_with_load(&ground, 1e-5, (1.0, 0.0), damping, &[0.0; 6])
            .unwrap();
        let observed = run
            .newmark_with_load(&ground, 1e-5, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(observed.u, original.u);
        assert_eq!(observed.kinetic_energy, original.kinetic_energy);
        assert_eq!(observed.internal_work, original.internal_work);
        assert_eq!(observed.damping_energy, original.damping_energy);
        assert_eq!(run.sample_count(), 301);
        let mut direct = StrainLifeMiner::new(
            curve(1.0),
            StrainLifeMeanStress::Uncorrected,
            config().history,
        )
        .unwrap();
        for u in &original.u {
            direct.push(u[3]).unwrap();
        }
        close(
            run.preview_finish(cx).unwrap().fibers[0].miner_damage,
            direct.finish().unwrap().miner_damage,
            1e-11,
        );
        split
            .newmark_with_load(&ground[..151], 1e-5, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        split
            .newmark_with_load(&ground[150..], 1e-5, (1.0, 0.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(split.sample_count(), run.sample_count());
        close(
            split.frame().displacements()[3],
            run.frame().displacements()[3],
            1e-12,
        );
        close(
            split.finish(cx).unwrap().fibers[0].miner_damage,
            run.finish(cx).unwrap().fibers[0].miner_damage,
            1e-11,
        );
    });
}

#[test]
fn native_budget_geometry_modulus_and_cancellation_refuse_without_publication() {
    with_cx(&CancelGate::new(), |cx| {
        let mut out_of_bounds = gauge(0, StrainLifeMeanStress::Uncorrected);
        out_of_bounds.section = 99;
        assert!(
            FrameFatiguePath::new(frame(&law(), false), vec![out_of_bounds], config(), cx).is_err()
        );
        let mut wrong = gauge(0, StrainLifeMeanStress::Uncorrected);
        let mut p = wrong.curve.parameters();
        p.youngs_modulus_pa *= 2.0;
        wrong.curve =
            Arc::new(StrainLifeCurve::new(wrong.curve.scope().clone(), p, 1.0, 1e12).unwrap());
        assert!(matches!(
            FrameFatiguePath::new(frame(&law(), false), vec![wrong], config(), cx),
            Err(FrameFatigueError::SourceMismatch { .. })
        ));
        let mut limited = config();
        limited.history.maximum_samples = 2;
        let mut run = FrameFatiguePath::new(
            frame(&law(), false),
            vec![gauge(0, StrainLifeMeanStress::Uncorrected)],
            limited,
            cx,
        )
        .unwrap();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.0005, 0.0], cx),
            Err(FrameFatigueError::Budget)
        ));
        assert_eq!(run.sample_count(), 1);
        assert_eq!(run.frame().displacements(), &[0.0; 6]);
    });
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut run = path(frame(&law(), false), StrainLifeMeanStress::Uncorrected, cx);
        gate.request();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.003], cx),
            Err(FrameFatigueError::Cancelled)
        ));
        assert!(matches!(
            run.preview_finish(cx),
            Err(FrameFatigueError::Cancelled)
        ));
        assert_eq!(run.sample_count(), 1);
        assert!(!run.is_finished());
    });
}

mod source_binding {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_material::fatigue::strain_life::source::{
        CmbPowerBasis, CmbProperty, CmbSourceRequest, resolve_cmb,
    };
    use fs_material::plastic::source::HardeningPropertyBasis;
    use fs_material::state_point::{MaterialPropertySelection, MaterialStatePointError};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
    use fs_solid::fiber::material::{FiberMaterial, RcSectionError};

    fn request() -> CmbSourceRequest {
        CmbSourceRequest {
            power_basis: CmbPowerBasis::Reversals,
            temperature_axis: "T".into(),
            stress_ratio_axis: "calibration-stress-ratio".into(),
            surface_and_geometry: "synthetic axial coupon".into(),
            environment: "synthetic dry atmosphere".into(),
            load_basis: "total mechanical strain amplitude".into(),
        }
    }
    fn point() -> QueryPoint {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap()
    }
    fn card(revision: u32) -> MaterialCard {
        card_with_temperature_quantity(revision, QuantitySpec::dimensional(Temperature::DIMS))
    }
    fn card_with_temperature_quantity(revision: u32, quantity: QuantitySpec) -> MaterialCard {
        let request = request();
        let temperature =
            ValidityDomain::unconstrained().with_quantity("T", quantity, 300.0, 300.0);
        let mut claims = ClaimSet::new();
        for property in CmbProperty::ALL {
            let key = property.property_key();
            let value = [E, 1e9, -0.5, 0.005, -0.5][property.index()];
            let validity = if property == CmbProperty::YoungsModulus {
                temperature.clone()
            } else {
                temperature
                    .clone()
                    .with_quantity(
                        &request.stress_ratio_axis,
                        QuantitySpec::dimensional(Dims::NONE),
                        -1.0,
                        -1.0,
                    )
                    .with_quantity(
                        request.power_basis.life_axis(),
                        QuantitySpec::dimensional(Dims::NONE),
                        2.0,
                        2e12,
                    )
            };
            claims
                .insert_claim(PropertyClaim {
                    value: PropertyValue::Scalar {
                        value,
                        dims: key.dims(),
                    },
                    key,
                    validity,
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: if property == CmbProperty::YoungsModulus {
                            "synthetic modulus"
                        } else {
                            "one synthetic CMB fit"
                        }
                        .into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        for (name, value, dims) in [
            ("density", 7850.0, Density::DIMS),
            ("poisson_ratio", 0.3, Dims::NONE),
            ("yield_stress", SY, Pressure::DIMS),
            ("isotropic_hardening_modulus", HI, Pressure::DIMS),
            ("kinematic_hardening_modulus", HK, Pressure::DIMS),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: temperature.clone(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic cyclic material".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        // A lineage starts at revision 0; each later revision supersedes the one before.
        let mut card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic frame alloy".into(),
                phase: "solid".into(),
                process: "synthetic cyclic coupon".into(),
                revision: 0,
            },
            claims.clone(),
            vec![],
        )
        .unwrap();
        for _ in 0..revision {
            card = MaterialCard::supersede(&card, claims.clone(), vec![]).unwrap();
        }
        card
    }
    fn material(card: &MaterialCard) -> Arc<FiberMaterial> {
        Arc::new(
            FiberMaterial::resolve_j2(
                card,
                &point(),
                "density",
                HardeningPropertyBasis::MixedPlasticStrainModuli,
                0.04,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        )
    }

    #[test]
    fn authored_curve_requires_a_native_absolute_temperature_quantity() {
        with_cx(&CancelGate::new(), |cx| {
            // An untyped pressure axis is a sampled stress, so the J2 source
            // refuses to freeze the density on it before any gauge exists
            // (fiber/material.rs, since e628b0746c). A static typed pressure is
            // a fixed external condition: the source resolves at it, and the
            // gauge must refuse it as the curve's temperature.
            let untyped = QuantitySpec::dimensional(Pressure::DIMS);
            assert!(matches!(
                FiberMaterial::resolve_j2(
                    &card_with_temperature_quantity(0, untyped),
                    &QueryPoint::new()
                        .with_quantity("T", untyped, 300.0)
                        .unwrap(),
                    "density",
                    HardeningPropertyBasis::MixedPlasticStrainModuli,
                    0.04,
                    MaterialPropertySelection::SingleClaimOnly,
                ),
                Err(RcSectionError::Source(
                    MaterialStatePointError::InvalidRequirement { .. }
                ))
            ));
            let typed_pressure = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Static,
            ));
            for (quantity, accepted) in [
                (QuantitySpec::dimensional(Temperature::DIMS), true),
                (typed_pressure, false),
            ] {
                let card = card_with_temperature_quantity(0, quantity);
                let point = QueryPoint::new()
                    .with_quantity("T", quantity, 300.0)
                    .unwrap();
                let material = Arc::new(
                    FiberMaterial::resolve_j2(
                        &card,
                        &point,
                        "density",
                        HardeningPropertyBasis::MixedPlasticStrainModuli,
                        0.04,
                        MaterialPropertySelection::SingleClaimOnly,
                    )
                    .unwrap(),
                );
                let authored = curve(0.5);
                let mut scope = authored.scope().clone();
                scope.material = "synthetic frame alloy".into();
                scope.process = "synthetic cyclic coupon".into();
                let gauge = FrameFiberFatigueGauge {
                    element: 0,
                    section: 0,
                    fiber: 0,
                    material_source: Some(material.clone()),
                    curve: Arc::new(
                        StrainLifeCurve::new(scope, authored.parameters(), 0.5, 1e12).unwrap(),
                    ),
                    mean_stress: StrainLifeMeanStress::Uncorrected,
                };
                let result =
                    FrameFatiguePath::new(frame(material.law(), false), vec![gauge], config(), cx);
                if accepted {
                    assert!(result.is_ok());
                } else {
                    assert!(matches!(
                        result,
                        Err(FrameFatigueError::SourceMismatch { .. })
                    ));
                }
            }
        });
    }

    #[test]
    fn original_cmb_card_state_receipts_and_native_law_must_all_agree() {
        with_cx(&CancelGate::new(), |cx| {
            let card = card(0);
            let material = material(&card);
            let request = request();
            let query = point()
                .with_quantity(
                    &request.stress_ratio_axis,
                    QuantitySpec::dimensional(Dims::NONE),
                    -1.0,
                )
                .unwrap()
                .with_quantity(
                    request.power_basis.life_axis(),
                    QuantitySpec::dimensional(Dims::NONE),
                    100.0,
                )
                .unwrap();
            let curve = Arc::new(
                resolve_cmb(
                    &card,
                    &query,
                    &request,
                    [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 5],
                )
                .unwrap(),
            );
            let correct = FrameFiberFatigueGauge {
                element: 0,
                section: 0,
                fiber: 0,
                material_source: Some(material.clone()),
                curve: curve.clone(),
                mean_stress: StrainLifeMeanStress::Morrow,
            };
            let mut run = FrameFatiguePath::new(
                frame(material.law(), false),
                vec![correct.clone()],
                config(),
                cx,
            )
            .unwrap();
            run.displacement_control(&PATTERN, 3, &[0.003, -0.003, 0.003], cx)
                .unwrap();
            assert!(run.finish(cx).unwrap().fibers[0].miner_damage > 0.0);
            assert_eq!(
                run.material_source(0)
                    .unwrap()
                    .source()
                    .unwrap()
                    .card_identity(),
                card.content_hash()
            );
            assert_eq!(
                run.curve(0).unwrap().source().unwrap().card_identity(),
                card.content_hash()
            );

            let mut missing = correct.clone();
            missing.material_source = None;
            assert!(matches!(
                FrameFatiguePath::new(frame(material.law(), false), vec![missing], config(), cx),
                Err(FrameFatigueError::SourceMismatch { .. })
            ));
            let mut wrong = correct.clone();
            wrong.material_source = Some(self::material(&self::card(1)));
            assert!(matches!(
                FrameFatiguePath::new(frame(material.law(), false), vec![wrong], config(), cx),
                Err(FrameFatigueError::SourceMismatch { .. })
            ));
            let different = FiberLaw::J2(
                UniaxialJ2::new(
                    J2Plasticity::new_mixed(
                        IsotropicElastic::new(E, 0.3, 0.04).unwrap(),
                        SY * 2.0,
                        HI,
                        HK,
                    )
                    .unwrap(),
                )
                .unwrap(),
            );
            assert!(matches!(
                FrameFatiguePath::new(frame(&different, false), vec![correct], config(), cx),
                Err(FrameFatigueError::SourceMismatch { .. })
            ));
        });
    }
}
