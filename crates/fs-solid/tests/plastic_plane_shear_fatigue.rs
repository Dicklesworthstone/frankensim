//! Native P1/Q1 shear-mode life; coefficients are synthetic implementation data.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::CycleClosureTolerance;
use fs_material::fatigue::critical_plane::fatemi_socie::source::{
    FatemiSociePowerBasis, FatemiSocieProperty, FatemiSocieSourceRequest,
    resolve_isotropic_fatemi_socie,
};
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieError, FatemiSocieLimits, FatemiSocieNormalStress,
};
use fs_material::fatigue::critical_plane::miner::CriticalPlaneMinerLimits;
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};
use fs_solid::plastic2::fatigue::{
    PlaneFatemiSocieGauge, PlaneFatigueConfig, PlaneFatigueError, PlaneFatigueGauge,
    PlaneFatiguePath, PlaneShearFatigueBudget,
};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const G: f64 = E / (2.0 * (1.0 + NU));
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 37,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn scope() -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic shear/normal fatigue fixture".into(),
        material: "plane shear test alloy".into(),
        process: "numerical coupon".into(),
        surface_and_geometry: "smooth synthetic coupon".into(),
        environment: "declared synthetic humidity".into(),
        load_basis: "original signed mechanical strain history".into(),
        temperature_k: 300.0,
    }
}

fn fixture(
    yield_stress: f64,
    hardening: f64,
    sensitivity: f64,
    humidity: f64,
) -> (ResolvedJ2Plasticity, Arc<FatemiSocieCalibration>) {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    // Humidity is a fixed external condition only as a typed static relative humidity; J2 and
    // every other source law refuse an untyped extra axis as a possibly evolving coordinate.
    let humidity_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::RelativeHumidity),
        ValueForm::Static,
    ));
    let domain = ValidityDomain::unconstrained()
        .with_quantity("T", temperature, 300.0, 300.0)
        .with_quantity("humidity", humidity_quantity, 0.0, 1.0);
    let request = FatemiSocieSourceRequest {
        power_basis: FatemiSociePowerBasis::Reversals,
        normal_stress: FatemiSocieNormalStress::Signed,
        temperature_axis: "T".into(),
        shear_stress_ratio_axis: "torsional-calibration-R".into(),
        surface_and_geometry: "smooth synthetic coupon".into(),
        environment: "declared synthetic humidity".into(),
        load_basis: "original reversed engineering shear calibration".into(),
    };
    let values = [E, NU, yield_stress, 0.02 * G, -0.5, 0.03, -0.5, sensitivity];
    let mut claims = ClaimSet::new();
    for property in FatemiSocieProperty::ALL {
        let key = property.property_key(request.normal_stress);
        let mut validity = domain.clone();
        if property.index() >= 3 {
            validity = validity
                .with_quantity(
                    &request.shear_stress_ratio_axis,
                    QuantitySpec::dimensional(Dims::NONE),
                    -1.0,
                    -1.0,
                )
                .with_quantity(
                    request.power_basis.life_axis(),
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                    2e12,
                );
        }
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value: values[property.index()],
                    dims: key.dims(),
                },
                key,
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic original elastic/yield/torsion calibration".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    for (name, value) in [
        ("isotropic_hardening_modulus", hardening),
        ("kinematic_hardening_modulus", 80.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Pressure::DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: Pressure::DIMS,
                },
                validity: domain.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native Prager/isotropic calibration".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: scope().material,
            phase: "solid".into(),
            process: scope().process,
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", temperature, 300.0)
        .unwrap()
        .with_quantity("humidity", humidity_quantity, humidity)
        .unwrap();
    let law = resolve_j2_plasticity(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
    )
    .unwrap();
    let fit = point
        .with_quantity(
            &request.shear_stress_ratio_axis,
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
    let calibration = resolve_isotropic_fatemi_socie(&card, &fit, &request, [SINGLE; 8]).unwrap();
    (law, Arc::new(calibration))
}

fn shear_gauge(
    law: &ResolvedJ2Plasticity,
    calibration: &Arc<FatemiSocieCalibration>,
) -> PlaneFatemiSocieGauge {
    PlaneFatemiSocieGauge {
        element: 0,
        quadrature_point: 0,
        material_state_identity: law.material().identity(),
        calibration: calibration.clone(),
        normals_in_mesh_frame: vec![[1.0, 0.0, 0.0]],
        limits: FatemiSocieLimits {
            maximum_samples: 8,
            maximum_planes: 8,
            maximum_projection_work: 64,
            maximum_circle_work: 4096,
        },
        closure: CycleClosureTolerance {
            absolute_strain: 1e-10,
            absolute_stress_pa: 1e-9,
            relative: 1e-10,
        },
    }
}
fn normal_gauge(law: &ResolvedJ2Plasticity) -> PlaneFatigueGauge {
    PlaneFatigueGauge {
        element: 0,
        quadrature_point: 0,
        material_state_identity: law.material().identity(),
        curve: Arc::new(
            StrainLifeCurve::new(
                scope(),
                CmbParameters {
                    youngs_modulus_pa: E,
                    fatigue_strength_coefficient_pa: 10.0,
                    fatigue_strength_exponent: -0.5,
                    fatigue_ductility_coefficient: 0.04,
                    fatigue_ductility_exponent: -0.5,
                },
                0.5,
                1e12,
            )
            .unwrap(),
        ),
        mean_stress: StrainLifeMeanStress::Uncorrected,
        normals_in_mesh_frame: vec![[1.0, 1.0, 0.0]],
        limits: CriticalPlaneMinerLimits {
            rainflow: RainflowLimits {
                maximum_samples: 32,
                maximum_residue: 16,
            },
            maximum_planes: 8,
            maximum_projection_work: 256,
        },
    }
}
fn config() -> PlaneFatigueConfig {
    PlaneFatigueConfig {
        maximum_gauges: 8,
        maximum_total_residue: 1024,
        maximum_total_projection_work: 2048,
    }
}
fn shear_budget() -> PlaneShearFatigueBudget {
    PlaneShearFatigueBudget {
        maximum_retained_samples: 64,
        maximum_planes: 16,
        maximum_projection_work: 512,
        maximum_circle_work: 32768,
    }
}
fn problem<'a>(
    mesh: &'a Mesh2,
    law: &'a ResolvedJ2Plasticity,
    fixed: &'a [(usize, f64)],
    loads: &'a [f64],
) -> SmallStrainPlaneProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainPlaneProblem {
        mesh,
        laws: SmallStrainLawField::Uniform(law),
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        thickness_m: 0.2,
        reduction: ReductionConfig {
            max_iterations: 40,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.02,
        },
        rotation_limit_rad: 0.02,
        prescribed_m: fixed,
        nodal_forces_n: &[],
        load_path: loads,
        budget: PlaneAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}
fn supports(mesh: &Mesh2, normal_to_shear: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(n, p)| [(2 * n, p[1] + normal_to_shear * p[0]), (2 * n + 1, 0.0)])
        .collect()
}
fn known_life(parameter: f64) -> f64 {
    0.5 * (0.05 / parameter).powi(2)
}

#[test]
fn g1_shear_miner_accumulates_explicit_cycles_across_native_restart_once() {
    let (law, calibration) = fixture(10.0, 0.0, 0.5, 0.2);
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed = supports(&mesh, 0.0);
    let loads = [0.004, -0.004, 0.0, 0.002, -0.002, 0.0];
    let p = problem(&mesh, &law, &fixed, &loads);
    with_cx(&CancelGate::new(), |cx| {
        let mut gauge = shear_gauge(&law, &calibration);
        gauge.normals_in_mesh_frame = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let mut path = PlaneFatiguePath::new_with_shear(
            &p,
            vec![normal_gauge(&law)],
            vec![gauge],
            config(),
            shear_budget(),
            cx,
        )
        .unwrap();
        assert_eq!(path.fatemi_socie_miner_summaries(), vec![None]);
        assert!(path.enable_fatemi_socie_miner(0, cx).is_err());
        assert_eq!(path.fatemi_socie_miner_summaries(), vec![None]);
        path.enable_fatemi_socie_miner(2, cx).unwrap();
        path.advance(3, cx).unwrap();
        let result = path.finish_fatemi_socie_cycle(cx).unwrap();
        let once = path.fatemi_socie_miner_summaries();
        let summary = once[0].unwrap();
        assert_eq!(summary.completed_cycles, 1);
        assert_eq!(summary.selected_plane, 0);
        close(summary.miner_damage, 1.0 / known_life(0.004), 1e-12);
        assert_eq!(path.finish_fatemi_socie_cycle(cx).unwrap(), result);
        assert_eq!(path.fatemi_socie_miner_summaries(), once);

        let mut resumed = path.clone();
        for run in [&mut path, &mut resumed] {
            run.restart_fatemi_socie_cycle(cx).unwrap();
            assert_eq!(run.fatemi_socie_miner_summaries(), once);
            assert!(run.enable_fatemi_socie_miner(100, cx).is_err());
            assert_eq!(run.fatemi_socie_miner_summaries(), once);
        }
        path.advance(3, cx).unwrap();
        resumed.advance(1, cx).unwrap();
        resumed.advance(2, cx).unwrap();
        assert_eq!(path.finish(cx).unwrap(), resumed.finish(cx).unwrap());
        assert_eq!(path.solution(), resumed.solution());
        assert_eq!(
            path.fatemi_socie_miner_summaries(),
            resumed.fatemi_socie_miner_summaries()
        );
        let summary = path.fatemi_socie_miner_summaries()[0].unwrap();
        assert_eq!(summary.completed_cycles, 2);
        // Pure shear has zero normal opening stress on the x plane. Each
        // explicitly closed window contributes one independently known life.
        close(
            summary.miner_damage,
            1.0 / known_life(0.004) + 1.0 / known_life(0.002),
            1e-12,
        );
        let planes = path.fatemi_socie_plane_damage(0).unwrap();
        close(planes[0].miner_damage, summary.miner_damage, 1e-14);
        assert_eq!(planes[1].miner_damage, 0.0);
        let closed = path.fatemi_socie_miner_summaries();
        path.finish(cx).unwrap();
        assert_eq!(path.fatemi_socie_miner_summaries(), closed);
    });
}

#[test]
fn native_p1_q1_shear_cycles_keep_plane_reduction_normal_stress_and_live_miner_history() {
    let (law, calibration) = fixture(10.0, 0.0, 0.5, 0.2);
    let amplitude = 0.004;
    let loads = [amplitude, -amplitude, amplitude];
    for mesh in [
        Mesh2::triangles(1.0, 1.0, 1, 1),
        Mesh2::quads(1.0, 1.0, 1, 1),
    ] {
        let fixed = supports(&mesh, 0.5);
        for plane in [PlaneKind::Stress, PlaneKind::Strain] {
            let mut p = problem(&mesh, &law, &fixed, &loads);
            p.plane = plane;
            with_cx(&CancelGate::new(), |cx| {
                let mut path = PlaneFatiguePath::new_with_shear(
                    &p,
                    vec![normal_gauge(&law)],
                    vec![shear_gauge(&law, &calibration)],
                    config(),
                    shear_budget(),
                    cx,
                )
                .unwrap();
                let mut normal =
                    PlaneFatiguePath::new(&p, vec![normal_gauge(&law)], config(), cx).unwrap();
                path.advance(1, cx).unwrap();
                let conditioned = path.solution().clone();
                path.restart_fatemi_socie_cycle(cx).unwrap();
                assert_eq!(path.solution(), &conditioned);
                let mut replay = path.clone();
                path.advance(2, cx).unwrap();
                replay.advance(1, cx).unwrap();
                replay.advance(1, cx).unwrap();
                normal.advance(3, cx).unwrap();
                assert_eq!(path.solution(), normal.solution());
                assert_eq!(path.solution(), replay.solution());
                let sample = path.fatemi_socie_samples()[0];
                close(sample.total_mechanical_strain[3], amplitude / 2.0, 1e-12);
                let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
                let normal_maximum = match plane {
                    PlaneKind::Stress => {
                        close(sample.stress_pa[2], 0.0, 1e-10);
                        close(
                            sample.total_mechanical_strain[2],
                            -NU / (1.0 - NU) * amplitude / 2.0,
                            1e-11,
                        );
                        E / (1.0 - NU * NU) * amplitude / 2.0
                    }
                    PlaneKind::Strain => {
                        assert_eq!(sample.total_mechanical_strain[2], 0.0);
                        close(sample.stress_pa[2], lambda * amplitude / 2.0, 1e-11);
                        (lambda + 2.0 * G) * amplitude / 2.0
                    }
                };
                let result = path.finish_fatemi_socie_cycle(cx).unwrap()[0];
                close(result.engineering_shear_strain_amplitude, amplitude, 1e-12);
                close(result.maximum_normal_stress_pa, normal_maximum, 1e-9);
                let correction = 1.0 + 0.5 * normal_maximum / 10.0;
                close(result.normal_stress_correction, correction, 1e-10);
                close(
                    result.cycles_to_failure.unwrap(),
                    known_life(amplitude * correction),
                    1e-6,
                );
                assert!(result.cycles_to_failure.unwrap() < known_life(amplitude));
                assert_eq!(path.fatemi_socie_sample_count(0), Some(3));
                assert_eq!(path.sample_count(), 4);
                assert_eq!(path.finish(cx).unwrap(), normal.finish(cx).unwrap());
                replay.finish(cx).unwrap();
                assert_eq!(path.fatemi_socie_results(), replay.fatemi_socie_results());
                assert_eq!(path.finish_fatemi_socie_cycle(cx).unwrap(), vec![result]);
            });
        }
    }
}

#[test]
fn native_plastic_shear_keeps_prager_history_and_refuses_unclosed_hardening_stress() {
    let amplitude = 0.01;
    let loads = [amplitude, -amplitude, amplitude, -amplitude, amplitude];
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed = supports(&mesh, 0.0);
    for hardening in [0.0, 20.0] {
        let (law, calibration) = fixture(1.0, hardening, 0.5, 0.2);
        let p = problem(&mesh, &law, &fixed, &loads);
        with_cx(&CancelGate::new(), |cx| {
            let mut path = PlaneFatiguePath::new_with_shear(
                &p,
                vec![],
                vec![shear_gauge(&law, &calibration)],
                config(),
                shear_budget(),
                cx,
            )
            .unwrap();
            path.advance(1, cx).unwrap();
            path.restart_fatemi_socie_cycle(cx).unwrap();
            path.advance(2, cx).unwrap();
            assert!(path.solution().quadrature_states[0].history.alpha > 0.01);
            if hardening > 0.0 {
                let native = path.solution().clone();
                assert!(matches!(
                    path.finish_fatemi_socie_cycle(cx),
                    Err(PlaneFatigueError::FatemiSocie {
                        error: FatemiSocieError::NotClosed {
                            quantity: "stress",
                            ..
                        },
                        ..
                    })
                ));
                assert_eq!(path.fatemi_socie_results(), vec![None]);
                assert_eq!(path.solution(), &native);
                return;
            }
            // Prager pure shear: tau - (2 Hk/3) epsilon_p_xy = sigma_y/sqrt(3).
            let plastic = (G * amplitude - 1.0 / 3.0_f64.sqrt()) / (2.0 * G + 160.0 / 3.0);
            close(
                path.solution().quadrature_states[0].history.plastic_strain[3],
                plastic,
                1e-10,
            );
            close(
                path.fatemi_socie_samples()[0].stress_pa[3],
                G * (amplitude - 2.0 * plastic),
                1e-9,
            );
            let result = path.finish_fatemi_socie_cycle(cx).unwrap()[0];
            close(result.engineering_shear_strain_amplitude, amplitude, 1e-12);
            close(
                result.cycles_to_failure.unwrap(),
                known_life(amplitude),
                1e-7,
            );
            assert!(matches!(
                path.advance(1, cx),
                Err(PlaneFatigueError::Invalid(_))
            ));
            path.restart_fatemi_socie_cycle(cx).unwrap();
            path.advance(2, cx).unwrap();
            assert!(path.finish(cx).unwrap().is_empty());
            close(
                path.fatemi_socie_results()[0]
                    .unwrap()
                    .cycles_to_failure
                    .unwrap(),
                result.cycles_to_failure.unwrap(),
                1e-7,
            );
            assert_eq!(
                path.fatemi_socie_calibration(0)
                    .unwrap()
                    .parameters()
                    .initial_yield_stress_pa,
                1.0
            );
        });
    }
}

#[test]
fn source_context_and_each_aggregate_shear_budget_are_checked_before_observation() {
    let (law, calibration) = fixture(10.0, 0.0, 0.5, 0.2);
    let (_, other_context) = fixture(10.0, 0.0, 0.5, 0.7);
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed = supports(&mesh, 0.0);
    let loads = [0.001, -0.001, 0.001];
    let p = problem(&mesh, &law, &fixed, &loads);
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(
            PlaneFatiguePath::new_with_shear(
                &p,
                vec![],
                vec![shear_gauge(&law, &other_context)],
                config(),
                shear_budget(),
                cx
            ),
            Err(PlaneFatigueError::SourceMismatch { .. })
        ));
        let exact = PlaneShearFatigueBudget {
            maximum_retained_samples: 8,
            maximum_planes: 1,
            maximum_projection_work: 64,
            maximum_circle_work: 4096,
        };
        let admitted = PlaneFatiguePath::new_with_shear(
            &p,
            vec![],
            vec![shear_gauge(&law, &calibration)],
            config(),
            exact,
            cx,
        )
        .unwrap();
        assert_eq!(
            admitted.fatemi_socie_material_source(0).unwrap().identity(),
            law.material().identity()
        );
        assert_eq!(admitted.fatemi_socie_sample_count(0), Some(1));
        for index in 0..4 {
            let mut short = exact;
            match index {
                0 => short.maximum_retained_samples -= 1,
                1 => short.maximum_planes -= 1,
                2 => short.maximum_projection_work -= 1,
                _ => short.maximum_circle_work -= 1,
            }
            assert!(matches!(
                PlaneFatiguePath::new_with_shear(
                    &p,
                    vec![],
                    vec![shear_gauge(&law, &calibration)],
                    config(),
                    short,
                    cx
                ),
                Err(PlaneFatigueError::Budget)
            ));
        }
        for wrong_yield in [false, true] {
            let mut parameters = calibration.parameters();
            if wrong_yield {
                parameters.initial_yield_stress_pa *= 2.0;
            } else {
                parameters.shear_modulus_pa *= 2.0;
            }
            let wrong =
                Arc::new(FatemiSocieCalibration::new(scope(), parameters, 0.5, 1e12).unwrap());
            assert!(matches!(
                PlaneFatiguePath::new_with_shear(
                    &p,
                    vec![],
                    vec![shear_gauge(&law, &wrong)],
                    config(),
                    shear_budget(),
                    cx
                ),
                Err(PlaneFatigueError::SourceMismatch { .. })
            ));
        }
    });
}

#[test]
fn late_shear_append_failure_restores_mechanics_miner_and_every_tensor_prefix() {
    let (law, calibration) = fixture(10.0, 0.0, 0.5, 0.2);
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed = supports(&mesh, 0.0);
    let loads = [0.004, -0.004, 0.004];
    let p = problem(&mesh, &law, &fixed, &loads);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let good = shear_gauge(&law, &calibration);
        let mut short = good.clone();
        short.limits.maximum_samples = 3;
        let mut path = PlaneFatiguePath::new_with_shear(
            &p,
            vec![normal_gauge(&law)],
            vec![good, short],
            config(),
            shear_budget(),
            cx,
        )
        .unwrap();
        path.enable_fatemi_socie_miner(2, cx).unwrap();
        let before = path.solution().clone();
        let normal = path.plane_summaries(0);
        let shear_damage = path.fatemi_socie_miner_summaries();
        assert!(matches!(
            path.advance(3, cx),
            Err(PlaneFatigueError::FatemiSocie {
                gauge: 1,
                error: FatemiSocieError::Budget
            })
        ));
        assert_eq!(path.solution(), &before);
        assert_eq!(path.plane_summaries(0), normal);
        assert_eq!(path.fatemi_socie_miner_summaries(), shear_damage);
        assert_eq!(path.sample_count(), 1);
        for index in 0..2 {
            assert_eq!(path.fatemi_socie_sample_count(index), Some(1));
        }
        // Continue from the accepted state and deliberately start a 3-sample cycle.
        path.advance(1, cx).unwrap();
        path.restart_fatemi_socie_cycle(cx).unwrap();
        path.advance(2, cx).unwrap();
        let observed = path.solution().clone();
        let normal = path.plane_summaries(0);
        gate.request();
        assert_eq!(path.finish(cx), Err(PlaneFatigueError::Cancelled));
        assert_eq!(
            path.restart_fatemi_socie_cycle(cx),
            Err(PlaneFatigueError::Cancelled)
        );
        assert_eq!(path.solution(), &observed);
        assert_eq!(path.plane_summaries(0), normal);
        assert_eq!(path.fatemi_socie_results(), vec![None, None]);
        assert_eq!(path.fatemi_socie_miner_summaries(), shear_damage);
    });
}

#[test]
fn late_shear_life_refusal_does_not_seal_either_fatigue_bank() {
    let (law, calibration) = fixture(10.0, 0.0, 0.5, 0.2);
    let narrow = Arc::new(
        FatemiSocieCalibration::new(scope(), calibration.parameters(), 1000.0, 1e12).unwrap(),
    );
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed = supports(&mesh, 0.0);
    let loads = [0.004, -0.004, 0.004];
    let p = problem(&mesh, &law, &fixed, &loads);
    with_cx(&CancelGate::new(), |cx| {
        let mut path = PlaneFatiguePath::new_with_shear(
            &p,
            vec![normal_gauge(&law)],
            vec![shear_gauge(&law, &calibration), shear_gauge(&law, &narrow)],
            config(),
            shear_budget(),
            cx,
        )
        .unwrap();
        path.enable_fatemi_socie_miner(2, cx).unwrap();
        path.advance(1, cx).unwrap();
        path.restart_fatemi_socie_cycle(cx).unwrap();
        path.advance(2, cx).unwrap();
        let native = path.solution().clone();
        let normal = path.plane_summaries(0);
        let shear_damage = path.fatemi_socie_miner_summaries();
        assert!(matches!(
            path.finish(cx),
            Err(PlaneFatigueError::FatemiSocie { gauge: 1, .. })
        ));
        assert_eq!(path.solution(), &native);
        assert_eq!(path.plane_summaries(0), normal);
        assert!(!path.is_finished());
        assert_eq!(path.results(), vec![None]);
        assert_eq!(path.fatemi_socie_results(), vec![None, None]);
        assert_eq!(path.fatemi_socie_miner_summaries(), shear_damage);
        assert_eq!(path.fatemi_socie_sample_count(0), Some(3));
        assert_eq!(path.fatemi_socie_sample_count(1), Some(3));
    });
}
