//! Calibrated shear-strain life from actual native J2 tensor trajectories.
//! All coefficients are synthetic implementation references, not measured life.
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::fatemi_socie::source::{
    FatemiSociePowerBasis, FatemiSocieProperty, FatemiSocieSourceRequest,
    resolve_isotropic_fatemi_socie,
};
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieError, FatemiSocieLimits, FatemiSocieNormalStress,
};
use fs_material::fatigue::critical_plane::{CriticalPlaneLimits, CycleClosureTolerance};
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::source::{HardeningPropertyBasis, ResolvedJ2Plasticity, j2_requirements};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::plastic3::dynamics::fatigue::{
    J2DynamicFatigueConfig, J2DynamicFatigueRun, J2DynamicFatigueTarget,
};
use fs_solid::plastic3::fatigue::{
    J2CriticalPlaneGauge, J2FatemiSocieGauge, J2FatigueConfig, J2FatigueError, J2FatigueGauge,
    J2FatiguePath,
};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FRAME: ContentHash = ContentHash([111; 32]);
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];

struct Fixture {
    law: ResolvedJ2Plasticity,
    calibration: Arc<FatemiSocieCalibration>,
    card: MaterialCard,
    point: QueryPoint,
    request: FatemiSocieSourceRequest,
}
fn fixture(plastic: bool, revision: u32) -> Fixture {
    let request = FatemiSocieSourceRequest {
        power_basis: FatemiSociePowerBasis::Reversals,
        normal_stress: FatemiSocieNormalStress::TensionOnly,
        temperature_axis: "T".into(),
        shear_stress_ratio_axis: "torsion-R".into(),
        surface_and_geometry: "synthetic smooth torsion specimen".into(),
        environment: "synthetic fixed humidity".into(),
        load_basis: "engineering shear-strain amplitude; original torsion fit".into(),
    };
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let scalar = QuantitySpec::dimensional(Dims::NONE);
    let domain = ValidityDomain::unconstrained()
        .with_quantity("T", temperature, 290.0, 330.0)
        .with_quantity("humidity", scalar, 0.0, 1.0);
    let point = QueryPoint::new()
        .with_quantity("T", temperature, 300.0)
        .unwrap()
        .with_quantity("humidity", scalar, 0.2)
        .unwrap();
    let fit_domain = domain
        .clone()
        .with_quantity(&request.shear_stress_ratio_axis, scalar, -1.0, -1.0)
        .with_quantity(request.power_basis.life_axis(), scalar, 2.0, 2e8);
    let youngs = 1000.0;
    let poisson = 0.3;
    let shear = youngs / (2.0 * (1.0 + poisson));
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value: f64, fit: bool| {
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity: if fit {
                    fit_domain.clone()
                } else {
                    domain.clone()
                },
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic source-bound torsion fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, youngs),
        ("poisson_ratio", Dims::NONE, poisson),
        (
            "yield_stress",
            Pressure::DIMS,
            if plastic { 1.0 } else { 100.0 },
        ),
        ("isotropic_hardening_modulus", Pressure::DIMS, 0.0),
        (
            "kinematic_hardening_modulus",
            Pressure::DIMS,
            if plastic { 80.0 } else { 0.0 },
        ),
        ("density", Density::DIMS, 6.0),
    ] {
        insert(PropertyKey::new(name, dims), value, false);
    }
    for (property, value) in [
        (
            FatemiSocieProperty::ShearFatigueStrengthCoefficient,
            0.1 * shear,
        ),
        (FatemiSocieProperty::ShearFatigueStrengthExponent, -0.5),
        (FatemiSocieProperty::ShearFatigueDuctilityCoefficient, 0.1),
        (FatemiSocieProperty::ShearFatigueDuctilityExponent, -0.5),
        (FatemiSocieProperty::NormalStressSensitivity, 0.2),
    ] {
        insert(property.property_key(request.normal_stress), value, true);
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shear metal".into(),
            phase: "solid".into(),
            process: "synthetic torsion coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    let basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
    let mut requirements = j2_requirements(basis).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law = ResolvedJ2Plasticity::try_from_material(bundle, basis, 0.05).unwrap();
    let fit_point = point
        .clone()
        .with_quantity(&request.shear_stress_ratio_axis, scalar, -1.0)
        .unwrap()
        .with_quantity(request.power_basis.life_axis(), scalar, 100.0)
        .unwrap();
    let calibration = Arc::new(
        resolve_isotropic_fatemi_socie(
            &card,
            &fit_point,
            &request,
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 8],
        )
        .unwrap(),
    );
    Fixture {
        law,
        calibration,
        card,
        point,
        request,
    }
}
fn scope(law: &ResolvedJ2Plasticity) -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic axial reference; not torsional measurements".into(),
        material: law.material().material().chemistry.clone(),
        process: law.material().material().process.clone(),
        surface_and_geometry: "synthetic smooth coupon".into(),
        environment: "synthetic fixed humidity".into(),
        load_basis: "axial total mechanical strain amplitude".into(),
        temperature_k: 300.0,
    }
}
fn axial_curve(law: &ResolvedJ2Plasticity) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            scope(law),
            CmbParameters {
                youngs_modulus_pa: law.law().elastic.youngs,
                fatigue_strength_coefficient_pa: 0.1 * law.law().elastic.youngs,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.1,
                fatigue_ductility_exponent: -0.5,
            },
            1.0,
            1e8,
        )
        .unwrap(),
    )
}
fn closure() -> CycleClosureTolerance {
    CycleClosureTolerance {
        absolute_strain: 1e-12,
        absolute_stress_pa: 1e-10,
        relative: 1e-10,
    }
}
fn fs_gauge(f: &Fixture) -> J2FatemiSocieGauge {
    J2FatemiSocieGauge {
        element: 0,
        material_state_identity: f.law.material().identity(),
        calibration: f.calibration.clone(),
        normals_in_mesh_frame: vec![[1.0, 0.0, 0.0]],
        limits: FatemiSocieLimits {
            maximum_samples: 128,
            maximum_planes: 32,
            maximum_projection_work: 4096,
            maximum_circle_work: 1_000_000,
        },
        closure: closure(),
    }
}
fn fixed_gauge(law: &ResolvedJ2Plasticity) -> J2FatigueGauge {
    J2FatigueGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 1.0, 0.0],
        material_state_identity: law.material().identity(),
        curve: axial_curve(law),
        mean_stress: StrainLifeMeanStress::Uncorrected,
    }
}
fn swt_gauge(law: &ResolvedJ2Plasticity) -> J2CriticalPlaneGauge {
    J2CriticalPlaneGauge {
        element: 0,
        material_state_identity: law.material().identity(),
        curve: axial_curve(law),
        normals_in_mesh_frame: vec![[1.0, 1.0, 0.0]],
        limits: CriticalPlaneLimits {
            maximum_samples: 128,
            maximum_planes: 32,
            maximum_projection_work: 4096,
        },
        closure: closure(),
    }
}
fn config(steps: usize) -> J2DynamicFatigueConfig {
    J2DynamicFatigueConfig {
        fatigue: J2FatigueConfig {
            history: RainflowLimits {
                maximum_samples: (steps + 1) as u64,
                maximum_residue: 64,
            },
            maximum_gauges: 8,
            maximum_total_residue: 512,
        },
        maximum_steps: steps,
    }
}
fn supports() -> Vec<(usize, f64)> {
    (0..12)
        .map(|d| (d, if d == 4 { 0.02 } else { 0.0 }))
        .collect()
}
fn problem<'a>(
    law: &'a ResolvedJ2Plasticity,
    supports: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: supports,
        nodal_forces_n: &[],
        load_path: path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 111,
                kernel_id: 111,
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
fn same(a: &J2DynamicFatigueRun<'_, '_>, b: &J2DynamicFatigueRun<'_, '_>) {
    assert_eq!(a.solution(), b.solution());
    assert_eq!(a.velocity_mps(), b.velocity_mps());
    assert_eq!(a.time_s(), b.time_s());
    assert_eq!(a.accepted_steps(), b.accepted_steps());
    assert_eq!(a.last_step(), b.last_step());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(a.summaries(), b.summaries());
    assert_eq!(a.critical_plane_results(), b.critical_plane_results());
    assert_eq!(a.fatemi_socie_results(), b.fatemi_socie_results());
}

#[test]
fn plastic_shear_uses_engineering_amplitude_and_keeps_cmb_swt_and_restart_continuous() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(true, 0);
        let supports = supports();
        let targets = [1.0, -1.0, 1.0, -1.0, 1.0];
        let p = problem(&f.law, &supports, &targets);
        let mut run = J2FatiguePath::new_with_fatemi_socie(
            &p,
            FRAME,
            vec![fixed_gauge(&f.law)],
            vec![swt_gauge(&f.law)],
            vec![fs_gauge(&f)],
            config(targets.len()).fatigue,
            cx,
        )
        .unwrap();
        run.advance(1, cx).unwrap();
        assert!(run.solution().element_states[0].alpha > 0.0);
        assert!(
            f.law
                .law()
                .backstress_pa(&run.solution().element_states[0])
                .iter()
                .any(|v| *v != 0.0)
        );
        close(run.solution().element_strain[0][3], 0.01, 1e-14);
        for cycle in 0..2 {
            let history = run.summaries();
            let increments = run.accepted_increments();
            run.restart_critical_plane_cycle(cx).unwrap();
            assert_eq!(run.summaries(), history);
            assert_eq!(run.accepted_increments(), increments);
            run.advance_cycle(2, cx).unwrap();
            let shear = run.fatemi_socie_results()[0].unwrap();
            assert_eq!(shear.samples, 3);
            // gamma_a=.02; tau_f/G=.1, gamma_f=.1, b=c=-.5:
            // .02=.2/sqrt(2N), so N=50. Tensor shear .01 would give N=200.
            close(shear.engineering_shear_strain_amplitude, 0.02, 1e-13);
            close(shear.fatemi_socie_parameter, 0.02, 1e-13);
            close(shear.cycles_to_failure.unwrap(), 50.0, 1e-8);
            assert!(
                run.critical_plane_results()[0]
                    .unwrap()
                    .cycles_to_failure
                    .unwrap()
                    > 0.0
            );
            assert_eq!(run.summaries()[0].samples, 4 + 2 * cycle);
        }
        let expected = p.solve(cx).unwrap();
        assert_eq!(run.solution().displacement_m, expected.displacement_m);
        assert_eq!(run.solution().element_states, expected.element_states);
        let result = run.fatemi_socie_results();
        assert!(run.finish(cx).unwrap()[0].finished);
        assert_eq!(run.fatemi_socie_results(), result);
        assert!(matches!(
            run.restart_critical_plane_cycle(cx),
            Err(J2FatigueError::Finished)
        ));
    });
}

#[test]
fn inertial_source_bound_shear_matches_the_original_backward_euler_solution_and_energy() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(true, 0);
        let supports = supports();
        let p = problem(&f.law, &supports, &[1.0]);
        let mut run = J2DynamicFatigueRun::new_with_fatemi_socie(
            &p,
            FRAME,
            Vec::new(),
            Vec::new(),
            vec![fs_gauge(&f)],
            config(3),
            0.0,
            ZERO.to_vec(),
            cx,
        )
        .unwrap();
        let mut solution = problem(&f.law, &supports, &[0.0]).solve(cx).unwrap();
        let mut velocity = ZERO.to_vec();
        let mut previous_time_s = 0.0;
        for (index, factor) in [1.0, -1.0, 1.0].into_iter().enumerate() {
            let prescribed: Vec<_> = supports
                .iter()
                .map(|&(d, value)| (d, factor * value))
                .collect();
            let endpoint = problem(&f.law, &prescribed, &[1.0]);
            let target = J2DynamicFatigueTarget {
                time_s: 0.05 * (index + 1) as f64,
                load_factor: factor,
            };
            let native = endpoint
                .advance_backward_euler(
                    solution.displacement_m.as_flattened(),
                    &velocity,
                    &solution.element_states,
                    &[6.0],
                    target.time_s - previous_time_s,
                    cx,
                )
                .unwrap();
            let report = if index == 2 {
                run.advance_cycle(&[target], cx).unwrap()
            } else {
                run.advance(&[target], cx).unwrap()
            };
            assert_eq!(run.solution(), &native.solution);
            assert_eq!(run.velocity_mps(), native.velocity_mps);
            assert_eq!(report.energy, native.energy);
            if index == 0 {
                run.restart_critical_plane_cycle(cx).unwrap();
            }
            solution = native.solution;
            velocity = native.velocity_mps;
            previous_time_s = target.time_s;
        }
        let result = run.fatemi_socie_results()[0].unwrap();
        close(result.cycles_to_failure.unwrap(), 50.0, 1e-8);
        assert_eq!(
            run.fatemi_socie_calibration(0)
                .unwrap()
                .source()
                .unwrap()
                .card_identity(),
            f.law.material().card_identity()
        );
        assert!(run.finish(cx).unwrap().is_empty());
    });
}

#[test]
fn late_second_shear_life_refusal_restores_every_bank_and_native_checkpoint_before_retry() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let f = fixture(true, 0);
        let supports = supports();
        let p = problem(&f.law, &supports, &[1.0]);
        let mut narrow = fs_gauge(&f);
        narrow.calibration = Arc::new(
            FatemiSocieCalibration::new(
                f.calibration.scope().clone(),
                f.calibration.parameters(),
                100.0,
                1e8,
            )
            .unwrap(),
        );
        let mut run = J2DynamicFatigueRun::new_with_fatemi_socie(
            &p,
            FRAME,
            vec![fixed_gauge(&f.law)],
            vec![swt_gauge(&f.law)],
            vec![fs_gauge(&f), narrow],
            config(3),
            0.0,
            ZERO.to_vec(),
            cx,
        )
        .unwrap();
        run.advance_to(0.05, 1.0, cx).unwrap();
        run.restart_critical_plane_cycle(cx).unwrap();
        let before = run.clone();
        assert!(matches!(
            run.advance_cycle(
                &[
                    J2DynamicFatigueTarget {
                        time_s: 0.1,
                        load_factor: -1.0
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.15,
                        load_factor: 1.0
                    },
                ],
                cx
            ),
            Err(J2FatigueError::FatemiSocie { gauge: 1, .. })
        ));
        same(&run, &before);
        // Same initial checkpoint and times, but half the engineering-shear
        // range: gamma_a=.01 => N=200, inside both declared life domains.
        run.advance_cycle(
            &[
                J2DynamicFatigueTarget {
                    time_s: 0.1,
                    load_factor: 0.0,
                },
                J2DynamicFatigueTarget {
                    time_s: 0.15,
                    load_factor: 1.0,
                },
            ],
            cx,
        )
        .unwrap();
        for result in run.fatemi_socie_results().into_iter().map(Option::unwrap) {
            assert_eq!(result.samples, 3);
            close(result.cycles_to_failure.unwrap(), 200.0, 1e-7);
        }
        let mut exhausted = run.clone();
        exhausted.restart_critical_plane_cycle(cx).unwrap();
        let checkpoint = exhausted.clone();
        assert!(matches!(
            exhausted.advance_to(0.2, -1.0, cx),
            Err(J2FatigueError::Budget)
        ));
        same(&exhausted, &checkpoint);
        let before_cancel = run.clone();
        gate.request();
        assert!(matches!(run.finish(cx), Err(J2FatigueError::Cancelled)));
        same(&run, &before_cancel);
    });
}

#[test]
fn source_card_context_initial_yield_and_native_shear_modulus_cannot_be_rebound() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(true, 0);
        let foreign = fixture(true, 1);
        let supports = supports();
        let p = problem(&f.law, &supports, &[1.0]);
        let mut wrong_yield = fs_gauge(&f);
        let mut parameters = f.calibration.parameters();
        parameters.initial_yield_stress_pa *= 2.0;
        wrong_yield.calibration = Arc::new(
            FatemiSocieCalibration::new(f.calibration.scope().clone(), parameters, 1.0, 1e8)
                .unwrap(),
        );
        let mut wrong_modulus = fs_gauge(&f);
        parameters = f.calibration.parameters();
        parameters.shear_modulus_pa *= 2.0;
        wrong_modulus.calibration = Arc::new(
            FatemiSocieCalibration::new(f.calibration.scope().clone(), parameters, 1.0, 1e8)
                .unwrap(),
        );
        let mut wrong_card = fs_gauge(&f);
        wrong_card.calibration = foreign.calibration;
        let scalar = QuantitySpec::dimensional(Dims::NONE);
        let changed_context = f
            .point
            .clone()
            .with_quantity("humidity", scalar, 0.8)
            .unwrap()
            .with_quantity(&f.request.shear_stress_ratio_axis, scalar, -1.0)
            .unwrap()
            .with_quantity(f.request.power_basis.life_axis(), scalar, 100.0)
            .unwrap();
        let mut wrong_context = fs_gauge(&f);
        wrong_context.calibration = Arc::new(
            resolve_isotropic_fatemi_socie(
                &f.card,
                &changed_context,
                &f.request,
                [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 8],
            )
            .unwrap(),
        );
        for bad in [wrong_yield, wrong_modulus, wrong_card, wrong_context] {
            assert!(matches!(
                J2DynamicFatigueRun::new_with_fatemi_socie(
                    &p,
                    FRAME,
                    Vec::new(),
                    Vec::new(),
                    vec![bad],
                    config(3),
                    0.0,
                    ZERO.to_vec(),
                    cx,
                ),
                Err(J2FatigueError::SourceMismatch { .. })
            ));
        }
    });
}

#[test]
fn decaying_unforced_native_shear_is_observed_but_cannot_be_fabricated_as_a_closed_cycle() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(false, 0);
        let supports: Vec<_> = (0..12).filter(|d| *d != 4).map(|d| (d, 0.0)).collect();
        let p = problem(&f.law, &supports, &[1.0]);
        let mut velocity = ZERO.to_vec();
        velocity[1][1] = 0.01;
        let mut run = J2DynamicFatigueRun::new_with_fatemi_socie(
            &p,
            FRAME,
            Vec::new(),
            Vec::new(),
            vec![fs_gauge(&f)],
            config(8),
            0.0,
            velocity,
            cx,
        )
        .unwrap();
        run.advance_to(0.005, 0.0, cx).unwrap();
        assert!(run.solution().element_strain[0][3] != 0.0);
        assert!(run.solution().element_stress[0][3] != 0.0);
        let before = run.clone();
        assert!(matches!(
            run.advance_cycle(
                &[
                    J2DynamicFatigueTarget {
                        time_s: 0.01,
                        load_factor: 0.0
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.015,
                        load_factor: 0.0
                    },
                ],
                cx
            ),
            Err(J2FatigueError::FatemiSocie {
                error: FatemiSocieError::NotClosed { .. },
                ..
            })
        ));
        same(&run, &before);
        run.advance_to(0.01, 0.0, cx).unwrap();
        assert_eq!(run.accepted_steps(), 2);
        assert_eq!(run.fatemi_socie_results(), vec![None]);
    });
}
