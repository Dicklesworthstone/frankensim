//! Original elastic/Newmark stress tensors feeding the Findley fatigue limit.
//! Every material and fatigue parameter is a synthetic implementation fixture.
#[path = "tet10/support.rs"]
mod support;

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::findley::source::{
    FindleyProperty, FindleySourceRequest, resolve_findley,
};
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyError, FindleyLimits, StressClosureTolerance,
};
use fs_material::fatigue::{FatigueScope, MeanStressCorrection, RainflowLimits, SnCurve, SnPoint};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::fatigue::{
    TetFatigueConfig, TetFatigueError, TetFatigueGauge, TetFatigueHistoryConfig, TetFatiguePath,
    TetFindleyGauge,
};
use fs_solid::linear3::source::SolidStiffnessModel;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use support::{block, with_cx};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([101; 32]);

fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(6.0, 60.0, 0.0, [101; 32]).unwrap()
}
fn thermal() -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        [101; 32],
        ContentHash([102; 32]),
        0.01,
    )
    .unwrap()
}
fn scope() -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic Findley/S-N arithmetic fixture; not measured".into(),
        material: "synthetic elastic solid".into(),
        process: "fixture".into(),
        surface_and_geometry: "smooth synthetic coupon".into(),
        environment: "synthetic fixed environment".into(),
        load_basis: "explicit native stress tensor cycle".into(),
        temperature_k: 300.0,
    }
}
fn sn_gauge(normal: [f64; 3]) -> TetFatigueGauge {
    TetFatigueGauge {
        element: 0,
        normal_in_mesh_frame: normal,
        material_state_identity: [101; 32],
        curve: Arc::new(
            SnCurve::new(
                scope(),
                vec![
                    SnPoint {
                        amplitude_pa: 1e-6,
                        cycles_to_failure: 1e7,
                    },
                    SnPoint {
                        amplitude_pa: 1.0,
                        cycles_to_failure: 10.0,
                    },
                ],
                Some(1e-8),
            )
            .unwrap(),
        ),
        mean_stress: MeanStressCorrection::Uncorrected,
    }
}
fn findley_gauge(normals: Vec<[f64; 3]>) -> TetFindleyGauge {
    TetFindleyGauge {
        element: 0,
        material_state_identity: [101; 32],
        material_source: None,
        calibration: Arc::new(FindleyCalibration::new(scope(), 0.2, 0.1).unwrap()),
        normals_in_mesh_frame: normals,
        limits: FindleyLimits {
            maximum_samples: 256,
            maximum_planes: 64,
            maximum_projection_work: 16384,
            maximum_circle_work: 1_000_000,
        },
        closure: StressClosureTolerance {
            absolute_stress_pa: 1e-12,
            relative: 1e-10,
        },
    }
}
fn config() -> TetFatigueHistoryConfig {
    TetFatigueHistoryConfig {
        maximum_absolute_mandel_strain: 0.01,
        history: RainflowLimits {
            maximum_samples: 1000,
            maximum_residue: 32,
        },
        maximum_gauges: 8,
        maximum_total_residue: 256,
    }
}
fn static_config() -> TetFatigueConfig {
    let config = config();
    TetFatigueConfig {
        solver: Default::default(),
        maximum_absolute_mandel_strain: config.maximum_absolute_mandel_strain,
        history: config.history,
        maximum_gauges: config.maximum_gauges,
        maximum_total_residue: config.maximum_total_residue,
    }
}
fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::default(),
    }
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn axial_initial() -> [[f64; 3]; 4] {
    let mut u = ZERO;
    u[1][0] = 0.001;
    u
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn quasistatic_tensors_select_shear_plane_and_keep_original_solve_and_sn_damage() {
    let mesh = block([1; 3], [1.0; 3]);
    let material = material();
    let thermal = thermal();
    let fixed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .flat_map(|(i, p)| (0..3).filter_map(move |axis| (p[axis] == 0.0).then_some(3 * i + axis)))
        .collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &mesh.nodes,
        tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    with_cx(|cx| {
        let geometry =
            TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 100, cx).unwrap();
        let right: Vec<_> = geometry
            .faces()
            .iter()
            .enumerate()
            .filter_map(|(i, face)| (face.outward_normal[0] > 0.99).then_some(i))
            .collect();
        let mut path = TetFatiguePath::new_with_findley(
            &problem,
            TetThermalStrainField::Uniform(&thermal),
            FRAME,
            vec![sn_gauge([1.0, 0.0, 0.0])],
            vec![findley_gauge(vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])],
            static_config(),
            cx,
        )
        .unwrap();
        assert!(path.restart_findley_cycle(cx).is_err());
        for (i, stress) in [0.06, -0.06, 0.06].into_iter().enumerate() {
            let load = geometry
                .assemble(
                    &[SurfaceLoad3::Traction {
                        faces: &right,
                        pascals: [stress, 0.0, 0.0],
                    }],
                    VolumeLoad3::None,
                    &[],
                    cx,
                )
                .unwrap();
            let observed = if i == 2 {
                path.advance_and_close_findley_cycle(i as f64, &load, cx)
                    .unwrap()
            } else {
                path.advance(i as f64, &load, cx).unwrap()
            };
            let original = load
                .solve(
                    &problem,
                    TetThermalStrainField::Uniform(&thermal),
                    Default::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap();
            assert_eq!(observed.mechanics.stress, original.stress);
            assert_eq!(observed.mechanics.displacement, original.displacement);
            close(observed.signed_stress_pa[0], stress, 1e-12);
        }
        let result = path.findley_results()[0].unwrap();
        assert_eq!(result.samples, 3);
        assert_eq!(result.selected_plane, 1);
        close(result.shear_amplitude_pa, 0.03, 1e-12);
        close(result.maximum_normal_stress_pa, 0.03, 1e-12);
        close(result.findley_parameter_pa, 0.036, 1e-12);
        close(result.usage_factor, 0.36, 1e-11);
        assert_eq!(path.close_findley_cycle(cx).unwrap()[0], result);
        close(path.finish(cx).unwrap()[0].miner_damage, 0.006, 1e-12);
        assert_eq!(path.findley_results()[0], Some(result));
        assert!(matches!(
            path.restart_findley_cycle(cx),
            Err(TetFatigueError::Finished)
        ));
    });
}

#[test]
fn rotating_shear_keeps_both_tangential_components_and_restarts_at_real_checkpoint() {
    with_cx(|cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut initial = ZERO;
        initial[1] = [0.0002, 0.002, 0.0];
        let mut native = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial,
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut run = native
            .clone()
            .with_fatigue_and_findley(
                FRAME,
                vec![sn_gauge([1.0, 1.0, 0.0])],
                vec![findley_gauge(vec![[1.0, 0.0, 0.0]])],
                config(),
                cx,
            )
            .unwrap();
        // n=x: sigma_n=.012; (tau_xy,tau_xz) follows a circle of radius .06.
        // Its scalar shear magnitude is CONSTANT, so range(norm(tau))/2=0 is wrong.
        let ring = [(0.0, 0.002), (-0.002, 0.0), (0.0, -0.002), (0.002, 0.0)];
        for cycle in 0..2 {
            for (i, (uy, uz)) in ring.into_iter().enumerate() {
                let supports = [(3, 0.0002), (4, uy), (5, uz)];
                let original = native
                    .advance_with_supports(&load, &supports, 0.05, cx)
                    .unwrap();
                let observed = if i == 3 {
                    run.advance_with_supports_and_close_findley_cycle(&load, &supports, 0.05, cx)
                        .unwrap()
                } else {
                    run.advance_with_supports(&load, &supports, 0.05, cx)
                        .unwrap()
                };
                assert_eq!(observed.mechanics, original);
                assert_eq!(
                    observed.response,
                    native.recover_response(FRAME, 0.01, cx).unwrap()
                );
                assert_eq!(run.mechanics().state(), native.state());
            }
            let result = run.findley_results()[0].unwrap();
            assert_eq!(result.samples, 5);
            close(result.shear_amplitude_pa, 0.06, 1e-14);
            close(result.maximum_normal_stress_pa, 0.012, 1e-14);
            close(result.findley_parameter_pa, 0.0624, 1e-14);
            close(result.usage_factor, 0.624, 1e-13);
            if cycle == 0 {
                let history = run.summaries();
                let samples = run.observed_samples();
                let state = run.mechanics().state().clone();
                assert!(
                    run.advance_with_supports(&load, &[(3, 0.0002), (5, 0.002)], 0.05, cx)
                        .is_err()
                );
                run.restart_findley_cycle(cx).unwrap();
                assert_eq!(run.mechanics().state(), &state);
                assert_eq!(run.observed_samples(), samples);
                assert_eq!(run.summaries(), history);
                assert_eq!(run.findley_results(), vec![None]);
            }
        }
        assert_eq!(run.observed_samples(), 9);
        close(run.finish(cx).unwrap()[0].miner_damage, 0.012, 1e-13);
    });
}

#[test]
fn free_vibration_uses_actual_newmark_stress_including_zero_external_load() {
    with_cx(|cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        for per_cycle in [16, 32] {
            // V=1/6, K=10, consistent M=.1, omega=10. Newmark phase per step
            // is 2*atan(omega*h/2); this is an exactly closed discrete period.
            let h = (std::f64::consts::PI / f64::from(per_cycle)).tan() / 5.0;
            let mut native = load
                .begin_transient(
                    &problem,
                    Some(TetThermalStrainField::Uniform(&thermal)),
                    0.0,
                    &axial_initial(),
                    &ZERO,
                    TransientConfig3::default(),
                    cx,
                )
                .unwrap();
            let mut run = native
                .clone()
                .with_fatigue_and_findley(
                    FRAME,
                    Vec::new(),
                    vec![findley_gauge(vec![[1.0, 1.0, 0.0]])],
                    config(),
                    cx,
                )
                .unwrap();
            for step in 1..=per_cycle {
                let original = native.advance(&load, h, cx).unwrap();
                let observed = if step == per_cycle {
                    run.advance_and_close_findley_cycle(&load, h, cx).unwrap()
                } else {
                    run.advance(&load, h, cx).unwrap()
                };
                assert_eq!(observed.mechanics, original);
                assert_eq!(
                    observed.response,
                    native.recover_response(FRAME, 0.01, cx).unwrap()
                );
                assert!(
                    observed
                        .response
                        .external_force_n
                        .iter()
                        .flatten()
                        .all(|f| *f == 0.0)
                );
                assert!(observed.fatigue.is_empty());
            }
            let result = run.findley_results()[0].unwrap();
            assert_eq!(result.samples, u64::try_from(per_cycle + 1).unwrap());
            close(result.shear_amplitude_pa, 0.03, 1e-13);
            close(result.findley_parameter_pa, 0.036, 1e-13);
            close(result.usage_factor, 0.36, 1e-12);
            assert!(run.finish(cx).unwrap().is_empty());
        }
    });
}

#[test]
fn late_findley_closure_refusal_rolls_back_native_motion_and_every_gauge_then_retries() {
    with_cx(|cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut native = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &axial_initial(),
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let strict = findley_gauge(vec![[1.0, 1.0, 0.0]]);
        let mut loose = strict.clone();
        loose.closure.absolute_stress_pa = 0.01;
        let mut run = native
            .clone()
            .with_fatigue_and_findley(
                FRAME,
                vec![sn_gauge([1.0, 0.0, 0.0])],
                vec![loose, strict],
                config(),
                cx,
            )
            .unwrap();
        native
            .advance_with_supports(&load, &[(3, -0.001)], 0.05, cx)
            .unwrap();
        run.advance_with_supports(&load, &[(3, -0.001)], 0.05, cx)
            .unwrap();
        let before = run.mechanics().state().clone();
        let history = run.summaries();
        // The first gauge admits a .006-Pa closure defect; the SECOND refuses.
        assert!(matches!(
            run.advance_with_supports_and_close_findley_cycle(&load, &[(3, 0.0009)], 0.05, cx),
            Err(TetFatigueError::Findley {
                gauge: 1,
                error: FindleyError::NotClosed { .. }
            })
        ));
        assert_eq!(run.mechanics().state(), &before);
        assert_eq!(run.summaries(), history);
        assert_eq!(run.observed_samples(), 2);
        assert_eq!(run.findley_results(), vec![None, None]);
        // Failure did not consume a native time step or leave gauge 0 closed.
        let original = native
            .advance_with_supports(&load, &[(3, 0.001)], 0.05, cx)
            .unwrap();
        let accepted = run
            .advance_with_supports_and_close_findley_cycle(&load, &[(3, 0.001)], 0.05, cx)
            .unwrap();
        assert_eq!(accepted.mechanics, original);
        assert_eq!(run.mechanics().state(), native.state());
        for result in run.findley_results().into_iter().map(Option::unwrap) {
            assert_eq!(result.samples, 3);
            close(result.usage_factor, 0.36, 1e-12);
        }
        close(run.finish(cx).unwrap()[0].miner_damage, 0.006, 1e-13);
    });
}

#[test]
fn findley_only_warmup_restart_preserves_the_global_native_sample_budget() {
    with_cx(|cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let native = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &axial_initial(),
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut bounded = config();
        bounded.history.maximum_samples = 3;
        let mut run = native
            .with_fatigue_and_findley(
                FRAME,
                Vec::new(),
                vec![findley_gauge(vec![[1.0, 1.0, 0.0]])],
                bounded,
                cx,
            )
            .unwrap();
        run.advance_with_supports(&load, &[(3, -0.001)], 0.05, cx)
            .unwrap();
        run.restart_findley_cycle(cx).unwrap();
        assert_eq!(run.observed_samples(), 2);
        run.advance_with_supports(&load, &[(3, 0.001)], 0.05, cx)
            .unwrap();
        run.restart_findley_cycle(cx).unwrap();
        let before = run.mechanics().state().clone();
        assert!(matches!(
            run.advance_with_supports(&load, &[(3, -0.001)], 0.05, cx),
            Err(TetFatigueError::Budget)
        ));
        assert_eq!(run.observed_samples(), 3);
        assert_eq!(run.mechanics().state(), &before);
    });
}

struct SourceFixture {
    card: MaterialCard,
    point: QueryPoint,
    bundle: Arc<ResolvedMaterialStatePoint>,
    material: TetElasticMaterial,
    thermal: TetThermalStrainState,
    calibration: Arc<FindleyCalibration>,
}
fn source_request() -> FindleySourceRequest {
    FindleySourceRequest {
        temperature_axis: "absolute_temperature".into(),
        surface_and_geometry: "synthetic smooth coupon".into(),
        environment: "synthetic humidity context".into(),
        load_basis: "synthetic calibrated multiaxial limit".into(),
    }
}
fn source_fixture(revision: u32) -> SourceFixture {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let request = source_request();
    let domain = ValidityDomain::unconstrained()
        .with_quantity(&request.temperature_axis, temperature, 290.0, 330.0)
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    let point = QueryPoint::new()
        .with_quantity(&request.temperature_axis, temperature, 300.0)
        .unwrap()
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.2)
        .unwrap();
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (PropertyKey::new("density", Density::DIMS), 6.0),
        (PropertyKey::new("young_modulus", Pressure::DIMS), 60.0),
        (PropertyKey::new("poisson_ratio", Dims::NONE), 0.0),
        (FindleyProperty::NormalStressSensitivity.property_key(), 0.2),
        (FindleyProperty::LimitFactor.property_key(), 0.1),
    ] {
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity: domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native Findley source fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "same synthetic named chemistry".into(),
            phase: "solid".into(),
            process: "same synthetic process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    let bundle = Arc::new(
        SolidStiffnessModel::Isotropic
            .resolve_reference(&card, &point, MaterialPropertySelection::SingleClaimOnly)
            .unwrap(),
    );
    let material = SolidStiffnessModel::Isotropic
        .material(
            &bundle,
            bundle.property("density").unwrap().value_si(),
            bundle.identity(),
            0.01,
        )
        .unwrap();
    let thermal = TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        bundle.identity().0,
        ContentHash([103; 32]),
        0.01,
    )
    .unwrap();
    let calibration = Arc::new(
        resolve_findley(
            &card,
            &point,
            &request,
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        )
        .unwrap(),
    );
    SourceFixture {
        card,
        point,
        bundle,
        material,
        thermal,
        calibration,
    }
}

#[test]
fn sourced_findley_requires_original_native_bundle_card_and_typed_operating_context() {
    with_cx(|cx| {
        let source = source_fixture(0);
        let foreign = source_fixture(1);
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&source.material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let native = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&source.thermal)),
                0.0,
                &axial_initial(),
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut bound = findley_gauge(vec![[1.0, 1.0, 0.0]]);
        bound.material_state_identity = source.material.material_state_identity;
        bound.material_source = Some(source.bundle.clone());
        bound.calibration = source.calibration.clone();
        let accepted = native
            .clone()
            .with_fatigue_and_findley(FRAME, Vec::new(), vec![bound.clone()], config(), cx)
            .unwrap();
        assert_eq!(
            accepted.findley_material_source(0).unwrap().identity(),
            source.bundle.identity()
        );
        assert_eq!(
            accepted
                .findley_calibration(0)
                .unwrap()
                .source()
                .unwrap()
                .card_identity(),
            source.card.content_hash()
        );
        let mut missing_bundle = bound.clone();
        missing_bundle.material_source = None;
        let mut foreign_bundle = bound.clone();
        foreign_bundle.material_source = Some(foreign.bundle);
        let mut foreign_fit = bound.clone();
        foreign_fit.calibration = foreign.calibration;
        let other_context = source
            .point
            .clone()
            .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.8)
            .unwrap();
        let mut foreign_context = bound;
        foreign_context.calibration = Arc::new(
            resolve_findley(
                &source.card,
                &other_context,
                &source_request(),
                [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
            )
            .unwrap(),
        );
        // Material labels and temperature match in all four cases. Immutable
        // material/source lineage and the full typed operating context do not.
        for bad in [missing_bundle, foreign_bundle, foreign_fit, foreign_context] {
            assert!(matches!(
                native
                    .clone()
                    .with_fatigue_and_findley(FRAME, Vec::new(), vec![bad], config(), cx,),
                Err(TetFatigueError::SourceMismatch { .. })
            ));
        }
    });
}
