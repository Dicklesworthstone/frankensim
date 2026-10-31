//! G1/G3/G4: native cyclic J2 state and bounded strain-life observation.
//! All property and fatigue data are synthetic implementation references.
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Pressure};
use fs_solid::plastic3::fatigue::{J2FatigueConfig, J2FatigueError, J2FatigueGauge, J2FatiguePath};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [2.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [2.0, 1.0, 1.0],
];
const TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 7, 5],
    [0, 2, 7, 3],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 7, 6],
];
const AMPLITUDE: f64 = 0.006;
const LOADS: [f64; 5] = [AMPLITUDE, -AMPLITUDE, AMPLITUDE, -AMPLITUDE, AMPLITUDE];

fn material() -> ResolvedJ2Plasticity {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 1000.0),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 1.0),
        ("isotropic_hardening_modulus", Pressure::DIMS, 20.0),
        ("kinematic_hardening_modulus", Pressure::DIMS, 80.0),
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
                    source: "synthetic cyclic J2 fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic metal".into(),
            phase: "solid".into(),
            process: "synthetic cyclic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_j2_plasticity(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
    )
    .unwrap()
}

fn curve(temperature: f64, maximum_cycles: f64) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic equal-exponent CMB analytic reference; not measured"
                    .into(),
                material: "synthetic metal".into(),
                process: "synthetic cyclic coupon".into(),
                surface_and_geometry: "smooth axial coupon".into(),
                environment: "synthetic fixture".into(),
                load_basis:
                    "total mechanical axial strain amplitude; full cycles; mean stress uncorrected"
                        .into(),
                temperature_k: temperature,
            },
            CmbParameters {
                youngs_modulus_pa: 1000.0,
                fatigue_strength_coefficient_pa: 10.0,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.04,
                fatigue_ductility_exponent: -0.5,
            },
            0.5,
            maximum_cycles,
        )
        .unwrap(),
    )
}
fn gauge(
    law: &ResolvedJ2Plasticity,
    curve: Arc<StrainLifeCurve>,
    element: usize,
) -> J2FatigueGauge {
    J2FatigueGauge {
        element,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: law.material().identity(),
        curve,
        mean_stress: StrainLifeMeanStress::Uncorrected,
    }
}
fn config() -> J2FatigueConfig {
    J2FatigueConfig {
        history: RainflowLimits {
            maximum_samples: 32,
            maximum_residue: 16,
        },
        maximum_gauges: 8,
        maximum_total_residue: 128,
    }
}
fn prescribed() -> Vec<(usize, f64)> {
    let mut result = Vec::new();
    for (node, x) in NODES.iter().enumerate() {
        result.push((3 * node, x[0]));
        if x[1] == 0.0 {
            result.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            result.push((3 * node + 2, 0.0));
        }
    }
    result
}
fn problem<'a>(
    law: &'a ResolvedJ2Plasticity,
    fixed: &'a [(usize, f64)],
    loads: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
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
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 5,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn g1_real_cyclic_plastic_state_and_hand_counted_strain_life() {
    let law = material();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    with_cx(&CancelGate::new(), |cx| {
        let mut path = J2FatiguePath::new(
            &p,
            ContentHash([9; 32]),
            vec![gauge(&law, curve(300.0, 10_000.0), 0)],
            config(),
            cx,
        )
        .unwrap();
        assert_eq!(path.summaries()[0].samples, 1);
        assert_eq!(path.samples()[0].mechanical_strain, 0.0);
        let (mut plastic, mut accumulated) = (0.0_f64, 0.0_f64);
        for epsilon in LOADS {
            // Independent 1D Prager / Hooke return mapping, compared against
            // the real 3D solve with free lateral displacements.
            let trial = 1000.0 * (epsilon - plastic) - 80.0 * plastic;
            let excess = trial.abs() - (1.0 + 20.0 * accumulated);
            if excess > 0.0 {
                let delta = excess / 1100.0;
                plastic += trial.signum() * delta;
                accumulated += delta;
            }
            path.advance(1, cx).unwrap();
            close(path.samples()[0].mechanical_strain, epsilon, 1e-12);
            close(
                path.samples()[0].normal_stress_pa,
                1000.0 * (epsilon - plastic),
                2e-8,
            );
            for state in &path.solution().element_states {
                close(state.plastic_strain[0], plastic, 1e-9);
                close(state.alpha, accumulated, 1e-9);
            }
        }
        assert!(accumulated > 0.01);
        let native = p.solve(cx).unwrap();
        assert_eq!(path.solution().element_states, native.element_states);
        assert_eq!(path.solution().element_strain, native.element_strain);
        let result = path.finish(cx).unwrap();
        assert_eq!(result[0].samples, 6);
        assert_eq!(2 * result[0].full_cycles + result[0].half_cycles, 5);
        // CMB b=c=-1/2 gives N_f = (0.01+0.04)^2/(2 epsilon_a^2).
        // Record 0,+a,-a,+a,-a,+a has one a/2 half cycle and two a cycles.
        close(
            result[0].miner_damage,
            4.25 * AMPLITUDE.powi(2) / 0.05_f64.powi(2),
            1e-12,
        );
        assert_eq!(path.finish(cx).unwrap(), result);
        assert_eq!(path.advance(1, cx), Err(J2FatigueError::Finished));
        assert_eq!(path.material(0).unwrap(), &law);
    });
}

#[test]
fn g3_chunking_and_cloning_preserve_plastic_and_open_cycle_histories() {
    let law = material();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    with_cx(&CancelGate::new(), |cx| {
        let mut whole = J2FatiguePath::new(
            &p,
            ContentHash([9; 32]),
            vec![gauge(&law, curve(300.0, 10_000.0), 0)],
            config(),
            cx,
        )
        .unwrap();
        let mut split = whole.clone();
        whole.advance(5, cx).unwrap();
        split.advance(2, cx).unwrap();
        assert_eq!(split.finish(cx), Err(J2FatigueError::PathIncomplete));
        let mut resumed = split.clone();
        split.advance(3, cx).unwrap();
        resumed.advance(1, cx).unwrap();
        resumed.advance(2, cx).unwrap();
        assert_eq!(whole.solution(), split.solution());
        assert_eq!(whole.solution(), resumed.solution());
        assert_eq!(whole.finish(cx).unwrap(), split.finish(cx).unwrap());
        assert_eq!(split.summaries(), resumed.finish(cx).unwrap());
    });
}

#[test]
fn g4_late_gauge_refusal_preserves_all_mechanics_and_all_counters() {
    let law = material();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    with_cx(&CancelGate::new(), |cx| {
        let mut path = J2FatiguePath::new(
            &p,
            ContentHash([9; 32]),
            vec![
                gauge(&law, curve(300.0, 10_000.0), 0),
                gauge(&law, curve(300.0, 10.0), 5),
            ],
            config(),
            cx,
        )
        .unwrap();
        let before = path.clone();
        assert!(matches!(
            path.advance(5, cx),
            Err(J2FatigueError::Gauge { gauge: 1, .. })
        ));
        assert_eq!(path.solution(), before.solution());
        assert_eq!(path.summaries(), before.summaries());
        assert_eq!(path.samples(), before.samples());
        assert_eq!(path.accepted_increments(), 0);
    });
}

#[test]
fn g4_cancellation_and_constitutive_strain_bound_roll_back() {
    let law = material();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut path = J2FatiguePath::new(
            &p,
            ContentHash([9; 32]),
            vec![gauge(&law, curve(300.0, 10_000.0), 0)],
            config(),
            cx,
        )
        .unwrap();
        path.advance(1, cx).unwrap();
        let before = path.clone();
        gate.request();
        assert_eq!(path.advance(1, cx), Err(J2FatigueError::Cancelled));
        assert_eq!(path.solution(), before.solution());
        assert_eq!(path.summaries(), before.summaries());
    });
    let outside = problem(&law, &fixed, &[0.06]);
    with_cx(&CancelGate::new(), |cx| {
        let mut path = J2FatiguePath::new(
            &outside,
            ContentHash([9; 32]),
            vec![gauge(&law, curve(300.0, 10_000.0), 0)],
            config(),
            cx,
        )
        .unwrap();
        let before = path.clone();
        assert!(matches!(
            path.advance(1, cx),
            Err(J2FatigueError::Invalid(_))
        ));
        assert_eq!(path.solution(), before.solution());
        assert_eq!(path.summaries(), before.summaries());
    });
}

#[test]
fn g0_source_temperature_and_history_budgets_refuse_before_start() {
    let law = material();
    let fixed = prescribed();
    let p = problem(&law, &fixed, &LOADS);
    with_cx(&CancelGate::new(), |cx| {
        let mut wrong = gauge(&law, curve(300.0, 10_000.0), 0);
        wrong.material_state_identity = ContentHash([0; 32]);
        assert!(matches!(
            J2FatiguePath::new(&p, ContentHash([9; 32]), vec![wrong], config(), cx),
            Err(J2FatigueError::SourceMismatch { .. })
        ));
        assert!(matches!(
            J2FatiguePath::new(
                &p,
                ContentHash([9; 32]),
                vec![gauge(&law, curve(301.0, 10_000.0), 0)],
                config(),
                cx
            ),
            Err(J2FatigueError::SourceMismatch { .. })
        ));
        let mut limited = config();
        limited.history.maximum_samples = 5;
        assert!(matches!(
            J2FatiguePath::new(
                &p,
                ContentHash([9; 32]),
                vec![gauge(&law, curve(300.0, 10_000.0), 0)],
                limited,
                cx
            ),
            Err(J2FatigueError::Budget)
        ));
        let mut limited = config();
        limited.maximum_total_residue = 15;
        assert!(matches!(
            J2FatiguePath::new(
                &p,
                ContentHash([9; 32]),
                vec![gauge(&law, curve(300.0, 10_000.0), 0)],
                limited,
                cx
            ),
            Err(J2FatigueError::Budget)
        ));
    });
}
