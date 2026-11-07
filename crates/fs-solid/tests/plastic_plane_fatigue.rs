//! G1/G3/G4: native plane plasticity to fatigue, using synthetic material data.
//! Numerical agreement below does not establish experimental component life.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::miner::CriticalPlaneMinerLimits;
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::plastic2::fatigue::{
    PlaneFatigueConfig, PlaneFatigueError, PlaneFatigueGauge, PlaneFatiguePath,
};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};
use std::sync::Arc;

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const A: f64 = 0.006;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 31,
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
fn material() -> ResolvedJ2Plasticity {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, E),
        ("poisson_ratio", Dims::NONE, NU),
        ("yield_stress", Pressure::DIMS, 1.0),
        ("isotropic_hardening_modulus", Pressure::DIMS, 20.0),
        ("kinematic_hardening_modulus", Pressure::DIMS, 80.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    300.0,
                    300.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic plane cyclic material".into(),
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
            process: "synthetic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_j2_plasticity(
        &card,
        &QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
    )
    .unwrap()
}
fn curve(minimum_cycles: f64) -> Arc<StrainLifeCurve> {
    Arc::new(
        StrainLifeCurve::new(
            FatigueScope {
                source_reference: "synthetic equal-exponent CMB formula, not measured".into(),
                material: "synthetic metal".into(),
                process: "synthetic coupon".into(),
                surface_and_geometry: "smooth numerical plane coupon".into(),
                environment: "synthetic fixed environment".into(),
                load_basis: "signed total mechanical normal strain".into(),
                temperature_k: 300.0,
            },
            CmbParameters {
                youngs_modulus_pa: E,
                fatigue_strength_coefficient_pa: 10.0,
                fatigue_strength_exponent: -0.5,
                fatigue_ductility_coefficient: 0.04,
                fatigue_ductility_exponent: -0.5,
            },
            minimum_cycles,
            1e12,
        )
        .unwrap(),
    )
}
fn gauge(law: &ResolvedJ2Plasticity) -> PlaneFatigueGauge {
    PlaneFatigueGauge {
        element: 0,
        quadrature_point: 0,
        material_state_identity: law.material().identity(),
        curve: curve(0.5),
        mean_stress: StrainLifeMeanStress::Uncorrected,
        normals_in_mesh_frame: vec![[1.0, 0.0, 0.0]],
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
fn axial_supports(mesh: &Mesh2) -> Vec<(usize, f64)> {
    let mut fixed = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        if p[0] == 0.0 || p[0] == 2.0 {
            fixed.push((2 * node, p[0]));
        }
        if p[1] == 0.0 {
            fixed.push((2 * node + 1, 0.0));
        }
    }
    fixed
}

#[test]
fn p1_and_q1_plane_stress_keep_native_plastic_reversal_and_hand_counted_damage() {
    let law = material();
    let loads = [A, -A, A, -A, A];
    for mesh in [
        Mesh2::triangles(2.0, 1.0, 2, 1),
        Mesh2::quads(2.0, 1.0, 2, 1),
    ] {
        let fixed = axial_supports(&mesh);
        let p = problem(&mesh, &law, &fixed, &loads);
        with_cx(&CancelGate::new(), |cx| {
            let mut path = PlaneFatiguePath::new(&p, vec![gauge(&law)], config(), cx).unwrap();
            let (mut plastic, mut accumulated) = (0.0_f64, 0.0_f64);
            for strain in loads {
                let trial = E * (strain - plastic) - 80.0 * plastic;
                let increment = (trial.abs() - 1.0 - 20.0 * accumulated).max(0.0) / 1100.0;
                plastic += trial.signum() * increment;
                accumulated += increment;
                path.advance(1, cx).unwrap();
                let sample = path.samples()[0];
                close(sample.total_mechanical_strain[0], strain, 1e-12);
                close(sample.stress_pa[0], E * (strain - plastic), 2e-8);
                close(sample.stress_pa[2], 0.0, 2e-10);
                close(
                    sample.total_mechanical_strain[2],
                    -NU * (strain - plastic) - 0.5 * plastic,
                    2e-9,
                );
                for point in &path.solution().quadrature_states {
                    close(point.history.plastic_strain[0], plastic, 2e-9);
                    close(point.history.alpha, accumulated, 2e-9);
                }
            }
            assert!(accumulated > 0.01);
            let native = p.solve(cx).unwrap();
            assert_eq!(native.quadrature_states, path.solution().quadrature_states);
            assert_eq!(
                native.quadrature_stress_pa,
                path.solution().quadrature_stress_pa
            );
            let preview = path.preview_finish(cx).unwrap();
            assert!(!path.is_finished());
            let result = path.finish(cx).unwrap();
            assert_eq!(result, preview);
            assert_eq!(result, path.finish(cx).unwrap());
            assert_eq!(result[0].samples, 6);
            // N_f = 0.05^2/(2 epsilon_a^2). One half cycle at a/2, two full at a.
            close(
                result[0].summary.miner_damage,
                4.25 * A.powi(2) / 0.05_f64.powi(2),
                1e-11,
            );
        });
    }
}

#[test]
fn plane_strain_swt_retains_constrained_transverse_stress_on_tilted_plane() {
    let law = material();
    let mesh = Mesh2::quads(2.0, 1.0, 1, 1);
    let fixed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .flat_map(|(n, p)| [(2 * n, p[0]), (2 * n + 1, p[1])])
        .collect();
    let a = 0.0002;
    let loads = [a, -a, a];
    let mut p = problem(&mesh, &law, &fixed, &loads);
    p.plane = PlaneKind::Strain;
    let mut g = gauge(&law);
    g.normals_in_mesh_frame = vec![[1.0, 0.0, 1.0]];
    g.mean_stress = StrainLifeMeanStress::SmithWatsonTopper;
    with_cx(&CancelGate::new(), |cx| {
        let mut path = PlaneFatiguePath::new(&p, vec![g], config(), cx).unwrap();
        path.advance(3, cx).unwrap();
        let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
        let mu = E / (2.0 * (1.0 + NU));
        close(path.samples()[0].stress_pa[2], 2.0 * lambda * a, 1e-12);
        assert_eq!(path.samples()[0].total_mechanical_strain[2], 0.0);
        let result = path.finish(cx).unwrap();
        // Plane epsilon_a=a/2, sigma_max=(mu+2 lambda)a. SWT b=c=-1/2:
        // 1/N = 4 sigma_max epsilon_a, plus the initial half cycle.
        close(
            result[0].summary.miner_damage,
            2.5 * (mu + 2.0 * lambda) * a * a,
            1e-11,
        );
    });
}

#[test]
fn q1_bbar_fatigue_uses_each_original_assumed_material_strain() {
    let law = material();
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let fixed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .flat_map(|(n, p)| [(2 * n, p[0] * p[1]), (2 * n + 1, 0.0)])
        .collect();
    let a = 0.0002;
    let loads = [a, -a, a];
    let mut p = problem(&mesh, &law, &fixed, &loads);
    p.plane = PlaneKind::Strain;
    p.formulation = Formulation::BBar;
    let gauges: Vec<_> = (0..4)
        .map(|point| {
            let mut g = gauge(&law);
            g.quadrature_point = point;
            g.normals_in_mesh_frame = vec![[0.0, 0.0, 1.0]];
            g
        })
        .collect();
    with_cx(&CancelGate::new(), |cx| {
        let mut path = PlaneFatiguePath::new(&p, gauges, config(), cx).unwrap();
        path.advance(3, cx).unwrap();
        let amplitude = a / (6.0 * 3.0_f64.sqrt());
        for sample in path.samples() {
            close(sample.total_mechanical_strain[2].abs(), amplitude, 1e-14);
        }
        for result in path.finish(cx).unwrap() {
            close(
                result.summary.miner_damage,
                2.25 * amplitude.powi(2) / 0.05_f64.powi(2),
                1e-12,
            );
        }
    });
}

#[test]
fn chunking_checkpoint_and_late_native_refusal_preserve_all_histories() {
    let law = material();
    let mesh = Mesh2::quads(2.0, 1.0, 1, 1);
    let fixed = axial_supports(&mesh);
    let loads = [A, -A, A, -A, A];
    let p = problem(&mesh, &law, &fixed, &loads);
    with_cx(&CancelGate::new(), |cx| {
        let mut one = PlaneFatiguePath::new(&p, vec![gauge(&law)], config(), cx).unwrap();
        let mut blocks = one.clone();
        one.advance(5, cx).unwrap();
        blocks.advance(2, cx).unwrap();
        let mut checkpoint = blocks.clone();
        let preview = blocks.preview_finish(cx).unwrap();
        assert_eq!(blocks.samples(), checkpoint.samples());
        blocks.advance(3, cx).unwrap();
        checkpoint.advance(1, cx).unwrap();
        checkpoint.advance(2, cx).unwrap();
        assert!(
            preview[0].summary.miner_damage
                < one.preview_finish(cx).unwrap()[0].summary.miner_damage
        );
        assert_eq!(one.solution(), blocks.solution());
        assert_eq!(one.finish(cx).unwrap(), blocks.finish(cx).unwrap());
        assert_eq!(one.finish(cx).unwrap(), checkpoint.finish(cx).unwrap());
    });
    let bad_loads = [A, -A, 0.2];
    let p = problem(&mesh, &law, &fixed, &bad_loads);
    with_cx(&CancelGate::new(), |cx| {
        let mut path = PlaneFatiguePath::new(&p, vec![gauge(&law)], config(), cx).unwrap();
        let before = path.clone();
        assert!(path.advance(3, cx).is_err());
        assert_eq!(path.solution(), before.solution());
        assert_eq!(path.plane_summaries(0), before.plane_summaries(0));
        assert_eq!(path.sample_count(), 1);
    });
}

#[test]
fn invalid_bindings_budgets_and_late_finalization_refuse_atomically() {
    let law = material();
    let mesh = Mesh2::quads(2.0, 1.0, 1, 1);
    let fixed = axial_supports(&mesh);
    let loads = [0.008];
    let p = problem(&mesh, &law, &fixed, &loads);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut bad = gauge(&law);
        bad.material_state_identity = ContentHash([7; 32]);
        assert!(matches!(
            PlaneFatiguePath::new(&p, vec![bad], config(), cx),
            Err(PlaneFatigueError::SourceMismatch { .. })
        ));
        let mut bad = gauge(&law);
        bad.quadrature_point = 4;
        assert!(matches!(
            PlaneFatiguePath::new(&p, vec![bad], config(), cx),
            Err(PlaneFatigueError::Invalid(_))
        ));
        let mut bounds = config();
        bounds.maximum_total_projection_work = 1;
        assert!(matches!(
            PlaneFatiguePath::new(&p, vec![gauge(&law)], bounds, cx),
            Err(PlaneFatigueError::Budget)
        ));
        let mut narrow = gauge(&law);
        narrow.curve = curve(100.0);
        let mut path = PlaneFatiguePath::new(&p, vec![gauge(&law), narrow], config(), cx).unwrap();
        assert_eq!(path.finish(cx), Err(PlaneFatigueError::PathIncomplete));
        path.advance(1, cx).unwrap();
        let before = path.clone();
        assert!(path.finish(cx).is_err());
        assert!(!path.is_finished());
        assert_eq!(path.results(), vec![None, None]);
        assert_eq!(path.solution(), before.solution());
        assert_eq!(path.plane_summaries(0), before.plane_summaries(0));
        gate.request();
        assert_eq!(path.preview_finish(cx), Err(PlaneFatigueError::Cancelled));
        assert_eq!(path.solution(), before.solution());
    });
}
