//! Native nonlinear cyclic hardening reaches fatigue without erasing history.
//! All coefficients are synthetic implementation references, not measured fits.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::miner::CriticalPlaneMinerLimits;
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::fiber::Uniaxial;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, ResolvedChabocheBranches, chaboche_branch_requirements,
};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::voce::source::ResolvedVoceJ2;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};
use fs_solid::plastic3::dynamics::fatigue::{
    J2DynamicFatigueConfig, J2DynamicFatigueRun, J2DynamicFatigueTarget,
};
use fs_solid::plastic3::fatigue::{
    J2FatigueConfig, J2FatigueGauge, J2FatiguePath, J2VariableAmplitudeGauge,
};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const E: f64 = 1200.0;
const SY: f64 = 1.0;
const HI: f64 = 30.0;
const HK: f64 = 20.0;
const Q: f64 = 0.6;
const RATE: f64 = 80.0;
const RHO: f64 = 6.0;
const A: f64 = 0.006;
const LOADS: [f64; 5] = [A, -A, A, -A, A];
const FRAME: ContentHash = ContentHash([41; 32]);
const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn branch_properties() -> [BackstressProperties; 3] {
    [
        BackstressProperties::new("slow_modulus", "slow_recovery"),
        BackstressProperties::new("medium_modulus", "medium_recovery"),
        BackstressProperties::new("fast_modulus", "fast_recovery"),
    ]
}

fn source(saturation: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, E),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, SY),
        ("isotropic_hardening_modulus", Pressure::DIMS, HI),
        ("kinematic_hardening_modulus", Pressure::DIMS, HK),
        ("voce_saturation_increment", Pressure::DIMS, saturation),
        ("voce_saturation_rate", Dims::NONE, RATE),
        ("slow_modulus", Pressure::DIMS, 80.0),
        ("slow_recovery", Dims::NONE, 30.0),
        ("medium_modulus", Pressure::DIMS, 300.0),
        ("medium_recovery", Dims::NONE, 250.0),
        ("fast_modulus", Pressure::DIMS, 140.0),
        ("fast_recovery", Dims::NONE, 1000.0),
        ("density", Density::DIMS, RHO),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic cyclic hardening fatigue fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear metal".into(),
            phase: "solid".into(),
            process: "synthetic cyclic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut requirements = chaboche_branch_requirements(
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        &branch_properties(),
    )
    .unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn voce(saturation: f64) -> ResolvedVoceJ2 {
    ResolvedVoceJ2::try_from_material(
        source(saturation),
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
    )
    .unwrap()
}

fn chaboche() -> ResolvedChabocheBranches<3> {
    ResolvedChabocheBranches::try_from_material(
        source(Q),
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
        branch_properties(),
    )
    .unwrap()
}

fn curve() -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic equal-exponent CMB; not measured".into(),
                material: "synthetic nonlinear metal".into(),
                process: "synthetic cyclic coupon".into(),
                surface_and_geometry: "smooth axial coupon".into(),
                environment: "synthetic fixed-condition fixture".into(),
                load_basis: "total mechanical strain amplitude; named mean-stress policy".into(),
                temperature_k: 300.0,
            },
            CmbParameters {
                youngs_modulus_pa: E,
                fatigue_strength_coefficient_pa: 0.01 * E,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.04,
                fatigue_ductility_exponent: -0.5,
            },
            0.5,
            1e12,
        )
        .unwrap(),
    )
}

fn gauge(law: &impl SourceBoundFatigueLaw, policy: StrainLifeMeanStress) -> J2FatigueGauge {
    J2FatigueGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: law.fatigue_source().identity(),
        curve: curve(),
        mean_stress: policy,
    }
}

fn config() -> J2FatigueConfig {
    J2FatigueConfig {
        history: RainflowLimits {
            maximum_samples: 32,
            maximum_residue: 16,
        },
        maximum_gauges: 4,
        maximum_total_residue: 128,
    }
}

fn prescribed() -> Vec<(usize, f64)> {
    // Exact free-sided axial response: yy and zz can contract; shear is zero.
    (0..12)
        .filter(|dof| *dof != 7 && *dof != 11)
        .map(|dof| (dof, if dof == 3 { 1.0 } else { 0.0 }))
        .collect()
}

fn problem<'a, L: SourceBoundFatigueLaw>(
    law: &'a L,
    fixed: &'a [(usize, f64)],
    loads: &'a [f64],
) -> SmallStrainTetProblem<'a, L> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: fixed,
        nodal_forces_n: &[],
        load_path: loads,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 41,
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
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn voce_cyclic_stress_matches_scalar_consistency_and_changes_swt_damage() {
    let law = voce(Q);
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    with_cx(|cx| {
        let mut path = J2FatiguePath::new(
            &p,
            FRAME,
            vec![
                gauge(&law, StrainLifeMeanStress::Uncorrected),
                gauge(&law, StrainLifeMeanStress::SmithWatsonTopper),
            ],
            config(),
            cx,
        )
        .unwrap();
        let (mut plastic, mut alpha) = (0.0_f64, 0.0_f64);
        for epsilon in LOADS {
            // Independent scalar consistency for free-sided isotropic/Prager
            // plasticity with Q(1-exp(-b alpha)); no native return map is used.
            let trial = E * (epsilon - plastic) - HK * plastic;
            let residual = |delta: f64| {
                trial.abs()
                    - (E + HK) * delta
                    - SY
                    - HI * (alpha + delta)
                    - Q * (1.0 - (-RATE * (alpha + delta)).exp())
            };
            if residual(0.0) > 0.0 {
                let (mut low, mut high) = (0.0, trial.abs() / E);
                for _ in 0..80 {
                    let mid = 0.5 * (low + high);
                    if residual(mid) > 0.0 {
                        low = mid;
                    } else {
                        high = mid;
                    }
                }
                let delta = 0.5 * (low + high);
                plastic += trial.signum() * delta;
                alpha += delta;
            }
            path.advance(1, cx).unwrap();
            close(
                path.samples()[0].normal_stress_pa,
                E * (epsilon - plastic),
                2e-8,
            );
            close(
                path.solution().element_states[0].plastic_strain[0],
                plastic,
                2e-10,
            );
            close(path.solution().element_states[0].alpha, alpha, 2e-10);
            assert_eq!(path.samples()[0], path.samples()[1]);
        }
        let native = p.solve(cx).unwrap();
        assert_eq!(path.solution().element_states, native.element_states);
        assert_eq!(path.solution().element_stress, native.element_stress);
        assert_eq!(path.material(0), Some(&law));
        let result = path.finish(cx).unwrap();
        close(
            result[0].miner_damage,
            4.25 * A.powi(2) / 0.05_f64.powi(2),
            1e-12,
        );

        let unsaturated = voce(0.0);
        let linear_problem = problem(&unsaturated, &fixed, &LOADS);
        let mut linear = J2FatiguePath::new(
            &linear_problem,
            FRAME,
            vec![gauge(&unsaturated, StrainLifeMeanStress::SmithWatsonTopper)],
            config(),
            cx,
        )
        .unwrap();
        linear.advance(LOADS.len(), cx).unwrap();
        let linear_damage = linear.finish(cx).unwrap()[0].miner_damage;
        assert!(result[1].miner_damage > linear_damage);
    });
}

#[test]
fn three_recovery_branches_survive_reversal_and_variable_amplitude_observation() {
    let law = chaboche();
    let axial = UniaxialChaboche::new(law.law().clone()).unwrap();
    let fixed = prescribed();
    let loads = [A, -0.004, 0.007, -0.005, 0.003];
    let p = problem(&law, &fixed, &loads);
    with_cx(|cx| {
        let mut path = J2FatiguePath::new_with_variable_amplitude(
            &p,
            FRAME,
            vec![gauge(&law, StrainLifeMeanStress::Uncorrected)],
            vec![],
            vec![],
            vec![J2VariableAmplitudeGauge {
                element: 0,
                material_state_identity: law.material().identity(),
                curve: curve(),
                mean_stress: StrainLifeMeanStress::Uncorrected,
                normals_in_mesh_frame: vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
                limits: CriticalPlaneMinerLimits {
                    rainflow: config().history,
                    maximum_planes: 3,
                    maximum_projection_work: 96,
                },
            }],
            config(),
            cx,
        )
        .unwrap();
        let mut expected = axial.initial_state();
        for epsilon in loads {
            let response = axial.response(epsilon, &expected).unwrap();
            path.advance(1, cx).unwrap();
            close(path.samples()[0].normal_stress_pa, response.stress_pa, 2e-8);
            let actual = &path.solution().element_states[0];
            close(
                actual.plastic.alpha,
                response.state.equivalent_plastic_strain(),
                2e-10,
            );
            for branch in 0..3 {
                for component in 0..6 {
                    close(
                        actual.backstress_pa[branch][component],
                        response.state.tensor_history().backstress_pa[branch][component],
                        2e-8,
                    );
                }
            }
            expected = response.state;
        }
        assert!(expected.equivalent_plastic_strain() > 0.01);
        assert!(
            expected
                .backstress_pa()
                .iter()
                .all(|stress| stress.abs() > 1e-4)
        );
        let preview = path.preview_variable_amplitude_fatigue(cx).unwrap();
        let final_fixed = path.finish(cx).unwrap();
        let result = path.variable_amplitude_results()[0].clone().unwrap();
        assert_eq!(preview[0], result);
        assert_eq!(result.selected_plane, 0);
        assert_eq!(result.summary, final_fixed[0]);
        assert_eq!(path.material(0), Some(&law));
    });
}

#[test]
fn nonlinear_dynamic_fatigue_retains_original_inertia_and_all_native_branches() {
    let law = chaboche();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &[1.0]);
    with_cx(|cx| {
        let mut run = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            vec![gauge(&law, StrainLifeMeanStress::Uncorrected)],
            J2DynamicFatigueConfig {
                fatigue: config(),
                maximum_steps: LOADS.len(),
            },
            0.0,
            vec![[0.0; 3]; NODES.len()],
            cx,
        )
        .unwrap();
        let mut native = problem(&law, &fixed, &[0.0]).solve(cx).unwrap();
        let mut velocity = vec![[0.0; 3]; NODES.len()];
        let mut time = 0.0;
        for (index, factor) in LOADS.into_iter().enumerate() {
            let next_time = 0.1 * (index + 1) as f64;
            let endpoint: Vec<_> = fixed
                .iter()
                .map(|&(dof, value)| (dof, value * factor))
                .collect();
            let motion = problem(&law, &endpoint, &[1.0])
                .advance_backward_euler(
                    native.displacement_m.as_flattened(),
                    &velocity,
                    &native.element_states,
                    &[RHO],
                    next_time - time,
                    cx,
                )
                .unwrap();
            let report = run.advance_to(next_time, factor, cx).unwrap();
            assert_eq!(run.solution(), &motion.solution);
            assert_eq!(run.velocity_mps(), &motion.velocity_mps);
            assert_eq!(report.energy, motion.energy);
            assert_eq!(
                run.samples()[0].normal_stress_pa,
                motion.solution.element_stress[0][0]
            );
            native = motion.solution;
            velocity = motion.velocity_mps;
            time = next_time;
        }
        assert_eq!(run.reference_density_kg_m3(), &[RHO]);
        assert_eq!(run.material(0), Some(&law));
        assert!(run.solution().element_states[0].plastic.alpha > 0.01);
        for branch in run.solution().element_states[0].backstress_pa {
            assert!(branch.iter().any(|stress| stress.abs() > 1e-4));
        }
        assert!(run.finish(cx).unwrap()[0].miner_damage > 0.0);
    });
}

#[test]
fn rejected_nonlinear_static_and_dynamic_blocks_restore_every_branch_and_counter() {
    let law = chaboche();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &[A, -A, 0.08]);
    let dynamic = problem(&law, &fixed, &[1.0]);
    with_cx(|cx| {
        let mut path = J2FatiguePath::new(
            &p,
            FRAME,
            vec![gauge(&law, StrainLifeMeanStress::Uncorrected)],
            config(),
            cx,
        )
        .unwrap();
        path.advance(1, cx).unwrap();
        let accepted = path.clone();
        assert!(path.advance(2, cx).is_err());
        assert_eq!(path.solution(), accepted.solution());
        assert_eq!(path.summaries(), accepted.summaries());
        assert_eq!(path.samples(), accepted.samples());
        assert_eq!(path.accepted_increments(), 1);

        let mut run = J2DynamicFatigueRun::new(
            &dynamic,
            FRAME,
            vec![gauge(&law, StrainLifeMeanStress::Uncorrected)],
            J2DynamicFatigueConfig {
                fatigue: config(),
                maximum_steps: 4,
            },
            0.0,
            vec![[0.0; 3]; NODES.len()],
            cx,
        )
        .unwrap();
        run.advance_to(0.1, A, cx).unwrap();
        let accepted = run.clone();
        assert!(
            run.advance(
                &[
                    J2DynamicFatigueTarget {
                        time_s: 0.2,
                        load_factor: -A
                    },
                    J2DynamicFatigueTarget {
                        time_s: 0.3,
                        load_factor: 0.08
                    },
                ],
                cx
            )
            .is_err()
        );
        assert_eq!(run.solution(), accepted.solution());
        assert_eq!(run.velocity_mps(), accepted.velocity_mps());
        assert_eq!(run.time_s(), accepted.time_s());
        assert_eq!(run.last_step(), accepted.last_step());
        assert_eq!(run.summaries(), accepted.summaries());
        assert_eq!(run.samples(), accepted.samples());
        assert_eq!(
            run.preview_fatigue(cx).unwrap(),
            accepted.preview_fatigue(cx).unwrap()
        );
    });
}
