//! Original transverse resultants reach fixed material-point Findley histories.
use super::*;
use fs_material::fatigue::MeanStressCorrection;
use fs_material::fatigue::critical_plane::findley::{
    FindleyError, FindleyLimits, StressClosureTolerance,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, LayerIsotropicParameters,
};
use fs_material::state_point::{
    IsotropicElasticParameters, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::{Density, Pressure};
use fs_solid::frame2d::fatigue::{FrameElasticFatigueGauge, FrameFindleyGauge};

const R: f64 = 0.05;
const L: f64 = 0.5;
const POINT: [f64; 2] = [0.04, 0.02];
const LOAD: [f64; 6] = [0.0, 0.0, 0.0, 2000.0, 10000.0, 0.0];

fn native(
    correction: f64,
    isotropic: bool,
) -> (Frame2d, FrameElasticFatigueGauge, FrameFindleyGauge) {
    let (material, curve, calibration) = elastic_layers::source_bundle(E, 0);
    let layers = vec![ConcentricMaterialLayer {
        outer_radius_m: R,
        material,
    }];
    let positive = ScalarAdmissibility::StrictlyPositive;
    let section = if isotropic {
        ConcentricElasticSection::try_new_with_isotropic_parameters(
            0.0,
            layers,
            &[LayerIsotropicParameters {
                density: ScalarPropertyRequirement::try_new("density", Density::DIMS, positive)
                    .unwrap(),
                parameters: IsotropicElasticParameters::YoungShear {
                    young_modulus: ScalarPropertyRequirement::try_new(
                        "young_modulus",
                        Pressure::DIMS,
                        positive,
                    )
                    .unwrap(),
                    shear_modulus: ScalarPropertyRequirement::try_new(
                        "shear_modulus",
                        Pressure::DIMS,
                        positive,
                    )
                    .unwrap(),
                },
            }],
        )
    } else {
        ConcentricElasticSection::try_new(0.0, layers)
    }
    .unwrap();
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_node(L, 0.0);
    frame.fix(0, [true; 3]);
    frame
        .add_material_shear_element(
            0,
            1,
            3,
            section
                .with_shear_modulus("shear_modulus", correction)
                .unwrap(),
        )
        .unwrap();
    let elastic = FrameElasticFatigueGauge {
        element: 0,
        section: 1,
        layer: 0,
        cross_section_position_m: POINT,
        temperature_axis: "absolute-temperature".into(),
        curve,
        mean_stress: MeanStressCorrection::Uncorrected,
    };
    let findley = FrameFindleyGauge {
        element: 0,
        section: 1,
        layer: 0,
        cross_section_position_m: POINT,
        temperature_axis: "absolute-temperature".into(),
        calibration,
        normals_in_material_frame: vec![[1.0, 0.0, 0.0]],
        limits: FindleyLimits {
            maximum_samples: 128,
            maximum_planes: 4,
            maximum_projection_work: 512,
            maximum_circle_work: 100_000,
        },
        closure: StressClosureTolerance {
            absolute_stress_pa: 1e-4,
            relative: 1e-8,
        },
        maximum_engineering_shear: 0.01,
    };
    (frame, elastic, findley)
}

fn observed(correction: f64, cx: &Cx<'_>) -> FrameFatiguePath {
    let (frame, elastic, findley) = native(correction, true);
    FrameFatiguePath::new_with_findley(
        frame,
        vec![],
        vec![elastic],
        vec![findley],
        128,
        config(),
        cx,
    )
    .unwrap()
}

// Independent Cartesian solid-circle solution and elementary cantilever statics.
// No call to the production stress recovery or fatigue evaluator.
fn stress(force_scale: f64) -> [f64; 6] {
    let [y, z] = POINT;
    let area = std::f64::consts::PI * R * R;
    let inertia = area * R * R / 4.0;
    let nu = 0.3;
    let q = force_scale * LOAD[4];
    let normal = force_scale * LOAD[3] / area - q * (L / 2.0) * y / inertia;
    let denominator = 8.0 * inertia * (1.0 + nu);
    let isotropic = (3.0 + 2.0 * nu) * R * R - (1.0 - 2.0 * nu) * (y * y + z * z);
    [
        normal,
        0.0,
        0.0,
        q * (isotropic - 2.0 * (1.0 + 2.0 * nu) * y * y) / denominator,
        0.0,
        -q * 2.0 * (1.0 + 2.0 * nu) * y * z / denominator,
    ]
}

#[test]
fn fixed_off_axis_stress_and_findley_use_actual_q_independently_of_authored_correction() {
    with_cx(&CancelGate::new(), |cx| {
        let mut deflections = Vec::new();
        for correction in [0.4, 0.9] {
            let mut run = observed(correction, cx);
            for amplitude in [1.0, -1.0, 0.0] {
                run.static_load(&LOAD.map(|v| amplitude * v), 1, cx)
                    .unwrap();
                let actual = run.findley_samples()[0];
                let expected = stress(amplitude);
                for (component, (a, b)) in actual.stress_pa.into_iter().zip(expected).enumerate() {
                    assert!((a - b).abs() <= 3e-9 * stress(1.0)[component].abs().max(1.0));
                }
                close(actual.transverse_force_n, amplitude * LOAD[4], 3e-9);
                close(actual.axial_mechanical_strain, expected[0] / E, 1e-12);
                assert_eq!(
                    actual.stress_pa[0],
                    run.elastic_samples()[0].normal_stress_pa
                );
                if amplitude == 1.0 {
                    deflections.push(run.frame().displacements()[4]);
                }
            }
            assert_eq!(run.sample_count(), 4);
            let source = run.findley_calibration(0).unwrap().source().unwrap();
            assert_eq!(
                source.card_identity(),
                run.findley_material_source(0).unwrap().card_identity()
            );
            let result = run.close_findley_cycle(cx).unwrap()[0];
            let tensor = stress(1.0);
            close(result.shear_amplitude_pa, tensor[3].hypot(tensor[5]), 3e-9);
            close(result.maximum_normal_stress_pa, tensor[0].abs(), 3e-9);
            close(
                result.usage_factor,
                (tensor[3].hypot(tensor[5]) + 0.2 * tensor[0].abs()) / 10e6,
                3e-9,
            );
            assert!(matches!(
                run.static_load(&LOAD, 1, cx),
                Err(FrameFatigueError::Findley {
                    error: FindleyError::Finished,
                    ..
                })
            ));
            run.restart_findley_cycle(cx).unwrap();
            for amplitude in [1.0, -1.0, 0.0] {
                run.static_load(&LOAD.map(|v| amplitude * v), 1, cx)
                    .unwrap();
            }
            close(
                run.close_findley_cycle(cx).unwrap()[0].usage_factor,
                result.usage_factor,
                3e-9,
            );
            assert_eq!(run.sample_count(), 7);
            run.finish(cx).unwrap();
        }
        assert!(deflections[0] > deflections[1]);
    });
}

#[test]
fn newmark_keeps_original_motion_energy_and_continuous_tensor_history() {
    with_cx(&CancelGate::new(), |cx| {
        let (mut frame, _, findley) = native(0.85, true);
        frame.static_load(&LOAD, 1).unwrap();
        let mut direct = frame.clone();
        let mut run = FrameFatiguePath::new_with_findley(
            frame,
            vec![],
            vec![],
            vec![findley],
            128,
            config(),
            cx,
        )
        .unwrap();
        let mut split = run.clone();
        let ground = [0.0; 41];
        let damping = Rayleigh { a0: 0.0, a1: 0.0 };
        let expected = direct
            .newmark_with_load(&ground, 1e-6, (0.0, 1.0), damping, &[0.0; 6])
            .unwrap();
        let actual = run
            .newmark_with_load(&ground, 1e-6, (0.0, 1.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(actual.u, expected.u);
        assert_eq!(actual.kinetic_energy, expected.kinetic_energy);
        assert_eq!(actual.internal_work, expected.internal_work);
        assert_eq!(actual.damping_energy, expected.damping_energy);
        split
            .newmark_with_load(&ground[..21], 1e-6, (0.0, 1.0), damping, &[0.0; 6], cx)
            .unwrap();
        split
            .newmark_with_load(&ground[20..], 1e-6, (0.0, 1.0), damping, &[0.0; 6], cx)
            .unwrap();
        assert_eq!(run.sample_count(), 41);
        assert_eq!(split.sample_count(), 41);
        for (a, b) in run.findley_samples()[0]
            .stress_pa
            .into_iter()
            .zip(split.findley_samples()[0].stress_pa)
        {
            close(a, b, 1e-8);
        }
        close(
            run.findley_samples()[0].transverse_force_n,
            direct.element_shear_response(0).unwrap().force_n,
            1e-12,
        );
    });
}

#[test]
fn late_shear_refusal_rolls_back_mechanics_scalar_history_and_every_tensor_suffix() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, mut elastic, mut findley) = native(0.85, true);
        elastic.section = 0; // Root stress keeps the synthetic S-N half cycles in its fitted band.
        let peak = (3.0 + 0.6) / (2.0 * 1.3) * LOAD[4] / (std::f64::consts::PI * R * R * (E / 2.6));
        findley.maximum_engineering_shear = 0.75 * peak;
        let mut early = findley.clone();
        early.maximum_engineering_shear = 0.01;
        findley.cross_section_position_m = [R, 0.0]; // Zero local transverse traction cannot hide the annular peak.
        let mut run = FrameFatiguePath::new_with_findley(
            frame,
            vec![],
            vec![elastic],
            vec![early, findley],
            256,
            config(),
            cx,
        )
        .unwrap();
        let mut unchanged = run.clone();
        assert!(run.static_load(&LOAD, 2, cx).is_err()); // First half-load is accepted privately.
        assert_eq!(
            run.frame().displacements(),
            unchanged.frame().displacements()
        );
        assert_eq!(run.sample_count(), unchanged.sample_count());
        assert_eq!(run.elastic_samples(), unchanged.elastic_samples());
        assert_eq!(run.findley_samples(), unchanged.findley_samples());
        assert_eq!(run.summaries(), unchanged.summaries());
        for target in [&mut run, &mut unchanged] {
            for amplitude in [0.5, -0.5, 0.0] {
                target
                    .static_load(&LOAD.map(|v| amplitude * v), 1, cx)
                    .unwrap();
            }
        }
        assert_eq!(
            run.close_findley_cycle(cx).unwrap(),
            unchanged.close_findley_cycle(cx).unwrap()
        );
    });
}

#[test]
fn source_geometry_temperature_and_declared_window_budgets_are_enforced() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, _, gauge) = native(0.85, true);
        for kind in 0..5 {
            let mut bad = gauge.clone();
            match kind {
                0 => bad.calibration = elastic_layers::source_bundle(E, 1).2,
                1 => bad.temperature_axis = "ambient-temperature".into(),
                2 => bad.cross_section_position_m = [2.0 * R, 0.0],
                3 => bad.section = usize::MAX,
                _ => bad.maximum_engineering_shear = f64::NAN,
            }
            assert!(
                FrameFatiguePath::new_with_findley(
                    frame.clone(),
                    vec![],
                    vec![],
                    vec![bad],
                    128,
                    config(),
                    cx
                )
                .is_err()
            );
        }
        assert!(matches!(
            FrameFatiguePath::new_with_findley(
                frame.clone(),
                vec![],
                vec![],
                vec![gauge.clone()],
                127,
                config(),
                cx
            ),
            Err(FrameFatigueError::Budget)
        ));
        let (nonisotropic, _, binding) = native(0.85, false);
        assert!(
            FrameFatiguePath::new_with_findley(
                nonisotropic,
                vec![],
                vec![],
                vec![binding],
                128,
                config(),
                cx
            )
            .is_err()
        );
        let mut limited = gauge;
        limited.limits.maximum_samples = 3;
        let mut run = FrameFatiguePath::new_with_findley(
            frame,
            vec![],
            vec![],
            vec![limited],
            3,
            config(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            run.static_load(&LOAD, 3, cx),
            Err(FrameFatigueError::Findley {
                error: FindleyError::Budget,
                ..
            })
        ));
        assert_eq!(run.sample_count(), 1);
    });
}

#[test]
fn failed_finalization_and_cancellation_leave_all_windows_live() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut run = observed(0.85, cx);
        run.static_load(&LOAD, 1, cx).unwrap();
        run.static_load(&LOAD.map(|v| -v), 1, cx).unwrap();
        let before = run.summaries();
        assert!(run.finish(cx).is_err());
        assert!(!run.is_finished());
        assert_eq!(run.findley_results(), vec![None]);
        assert_eq!(run.summaries(), before);
        run.static_load(&[0.0; 6], 1, cx).unwrap();
        run.finish(cx).unwrap();
        assert!(run.findley_results()[0].is_some());
    });
    let mut run = with_cx(&CancelGate::new(), |cx| observed(0.85, cx));
    gate.request();
    with_cx(&gate, |cx| {
        assert_eq!(
            run.static_load(&LOAD, 1, cx).unwrap_err(),
            FrameFatigueError::Cancelled
        );
        assert_eq!(run.sample_count(), 1);
        assert_eq!(run.findley_results(), vec![None]);
    });
}
