//! The ordinary product projection must obtain a real bound beyond 256 DOFs.
use super::*;
use crate::json_read::JsonValue;
use fs_airflow::conjugate::{AirSegment, goal::maximum::prepare_linear_maximum};
use fs_conduction::{
    ConductionMesh, ConductivityModel, ScalarField, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackBoundStatus};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 7306, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx);
    });
}

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
    paths: Vec<AirPath>,
}
impl Fixture {
    fn new(cells: usize, shear: f64) -> Self {
        let (complex, mut points) = fs_conduction::fixtures::unit_cube(cells);
        for point in &mut points { point[0] += shear * point[1]; }
        let mesh = ConductionMesh::new(complex, points).unwrap();
        let area: f64 = mesh.boundary().iter().map(|face| face.area).sum();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap().finish().unwrap();
        let paths = vec![AirPath::new(330.0, 1.0, 20.0 * area,
            vec![AirSegment::new("air", area, 2.0).unwrap()]).unwrap()];
        Self { mesh, boundary, paths,
            material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(0.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, source: &self.source, element_materials: None }
    }
    fn linear(&self) -> LinearConfig {
        LinearConfig { tolerance: 1e-12, max_iterations: 2000, restart: 32 }
    }
    fn solid(&self) -> LinearGoalAnalysisConfig {
        LinearGoalAnalysisConfig {
            residual_limits: GoalResidualLimits { max_rows: 1000, max_nonzeros: 100_000 },
            // Deliberately disable numerical scaling/dense inverse proposals.
            max_stability_iterations: 0,
        }
    }
    fn assess(&self, cx: &Cx<'_>, field: &[f64], memory: u64) -> MaximumEvidence {
        maximum_evidence(cx, self.problem(), None, &self.paths, self.linear(), field,
            &(0..field.len()).collect::<Vec<_>>(), memory, self.solid(), 1e-4).unwrap()
    }
}
fn bound(evidence: &MaximumEvidence) -> f64 {
    match evidence.term.as_ref().unwrap() {
        PropagatedTerm::Measured { half_width_k, method, vertices, .. } => {
            assert_eq!(*method, METHOD);
            assert!(vertices.is_empty());
            *half_width_k
        }
        other => panic!("no measured coupled bound: {other:?}"),
    }
}

#[test]
fn product_uses_spectral_proof_on_an_obtuse_mesh_beyond_dense_cap() {
    let fixture = Fixture::new(6, 1.5);
    assert_eq!(fixture.mesh.vertex_count(), 343);
    with_cx(|_, cx| {
        let field = vec![300.0; fixture.mesh.vertex_count()];
        let vertices: Vec<_> = (0..field.len()).collect();
        let ordinary = prepare_linear_maximum(cx, fixture.problem(), None, &fixture.paths,
            fixture.linear(), &field, fixture.solid(),
            policy(64 * 1024 * 1024, fixture.solid(), fixture.linear())).unwrap();
        let previous = ordinary.analyze_maximum(cx, &field, &vertices).unwrap();
        assert_eq!(previous.coupled().status(), FeedbackBoundStatus::SolidInverseUnavailable,
            "fixture must need the newly connected proof");
        let evidence = fixture.assess(cx, &field, 64 * 1024 * 1024);
        assert!(bound(&evidence) >= 30.0, "the coupled equilibrium is 330 K, not 300 K");
        assert!(field.iter().all(|v| *v == 300.0));
        assert_eq!(evidence.primal_iterations, 0);
        let json = JsonValue::parse(evidence.control_json.as_ref().unwrap()).unwrap();
        let proof = json.get("solid_spectral").expect("actual sparse witness retained");
        assert_eq!(proof.str_field("stop"), Some("Certified"));
        assert!(proof.f64_field("coercivity_lower").unwrap() > 0.0);
        assert!(proof.f64_field("work_entries").unwrap() > 0.0);
        assert!(proof.f64_field("work_entries").unwrap()
            <= proof.f64_field("max_shared_verification_entries").unwrap());
        assert!(proof.f64_field("peak_storage_entries").unwrap()
            <= proof.f64_field("max_storage_entries").unwrap());
        assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(false)));
        assert_eq!(json.get("candidate_accepted"), Some(&JsonValue::Bool(false)));
        let near = fixture.assess(cx, &vec![330.0; field.len()], 64 * 1024 * 1024);
        assert!(bound(&near) < 1e-4, "a near-equilibrium field must have a useful bound");
        let repeat = fixture.assess(cx, &field, 64 * 1024 * 1024);
        assert_eq!(evidence.control_json, repeat.control_json);
        let terms = super::super::super::propagation_term_receipts(
            None, None, evidence.term.as_ref(), fs_blake3::hash_bytes(b"spectral-product"),
        ).unwrap();
        assert_eq!(terms.iter().filter(|term|
            term.kind() == fs_evidence::uncertainty::EngineeringUncertaintyKind::SolverAlgebraic
        ).count(), 1);
    });
}

#[test]
fn product_sparse_budget_stop_is_no_data_with_paid_work_not_zero_error() {
    let fixture = Fixture::new(6, 1.5);
    with_cx(|_, cx| {
        let evidence = fixture.assess(cx, &vec![300.0; fixture.mesh.vertex_count()], 1024 * 1024);
        assert!(matches!(evidence.term, Some(PropagatedTerm::Unmeasured { .. })));
        let json = JsonValue::parse(evidence.control_json.as_ref().unwrap()).unwrap();
        let proof = json.get("solid_spectral").expect("failed proof diagnostics retained");
        assert!(matches!(proof.str_field("stop"), Some("WorkLimit" | "StorageLimit")));
        assert!(proof.f64_field("work_entries").unwrap() > 0.0);
        assert_eq!(proof.get("coercivity_lower"), Some(&JsonValue::Null));
        assert_eq!(json.get("final_bound_k"), Some(&JsonValue::Null));
        assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(false)));
    });
}

#[test]
fn spectral_policy_does_not_expand_solid_or_shared_work_limits() {
    let fixture = Fixture::new(1, 0.0);
    for memory in [0, 127, 128, 1024 * 1024, u64::MAX] {
        let config = policy(memory, fixture.solid(), fixture.linear());
        let spectral = spectral_policy(memory, config);
        assert_eq!(spectral.system, config.residual.solid);
        assert_eq!(spectral.max_work_entries, config.residual.max_verification_entries);
        assert_eq!(spectral.max_storage_entries,
            usize::try_from(memory / 128).unwrap_or(usize::MAX));
        assert_eq!(spectral.max_shift_attempts, 16);
    }
}

#[test]
fn successful_existing_bounds_do_not_acquire_a_spectral_claim() {
    let fixture = Fixture::new(1, 0.0);
    with_cx(|_, cx| {
        let field = vec![300.0; fixture.mesh.vertex_count()];
        let analyzer = prepare_linear_maximum(cx, fixture.problem(), None, &fixture.paths,
            fixture.linear(), &field, fixture.solid(),
            policy(64 * 1024 * 1024, fixture.solid(), fixture.linear())).unwrap();
        let ordinary = analyzer.analyze_maximum(cx, &field, &(0..field.len()).collect::<Vec<_>>()).unwrap();
        assert!(ordinary.algebraic_half_width_k().is_some());
        let evidence = fixture.assess(cx, &field, 64 * 1024 * 1024);
        assert_eq!(bound(&evidence).to_bits(), ordinary.algebraic_half_width_k().unwrap().to_bits());
        let json = JsonValue::parse(evidence.control_json.as_ref().unwrap()).unwrap();
        assert!(json.get("solid_spectral").is_none());
    });
}
