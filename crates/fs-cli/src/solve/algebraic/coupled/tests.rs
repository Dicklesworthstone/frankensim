use super::*;
use crate::json_read::JsonValue;
use fs_airflow::conjugate::AirSegment;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solver::goal::GoalResidualLimits;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 7305, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
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
    fn new() -> Self {
        let (complex, positions) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let area = mesh.boundary().iter().map(|face| face.area).sum();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap().finish().unwrap();
        let paths = vec![AirPath::new(330.0, 1.0, 24.0,
            vec![AirSegment::new("air", area, 2.0).unwrap()]).unwrap()];
        Self { mesh, boundary, paths,
            material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(0.0) }
    }
    fn assess(&self, cx: &Cx<'_>, field: &[f64], memory: u64, requested: f64)
        -> Result<MaximumEvidence, SolveRefusal>
    {
        maximum_evidence(cx, ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            element_materials: None, source: &self.source,
        }, None, &self.paths,
        LinearConfig { tolerance: 1e-11, max_iterations: 100, restart: 20 },
        field, &(0..field.len()).collect::<Vec<_>>(), memory,
        LinearGoalAnalysisConfig {
            residual_limits: GoalResidualLimits { max_rows: 100, max_nonzeros: 10_000 },
            max_stability_iterations: 100,
        }, requested)
    }
}
fn measured(evidence: &MaximumEvidence) -> f64 {
    match evidence.term.as_ref().unwrap() {
        PropagatedTerm::Measured { half_width_k, method, vertices, .. } => {
            assert_eq!(*method, METHOD);
            assert!(vertices.is_empty(), "no invented re-solve samples");
            *half_width_k
        }
        other => panic!("expected the real coupled bound, got {other:?}"),
    }
}

#[test]
fn coupled_product_term_sees_the_error_hidden_by_frozen_robin_and_reaches_qoi() {
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        let field = vec![300.0; fixture.mesh.vertex_count()];
        let before = field.clone();
        let evidence = fixture.assess(cx, &field, 64 * 1024 * 1024, 1.0).unwrap();
        // The no-source, uniform-inlet coupled equilibrium is 330 K, not
        // the 300 K equilibrium of the frozen solid boundary.
        assert!(measured(&evidence) >= 30.0);
        assert_eq!(field, before);
        assert_eq!(evidence.primal_iterations, 0);
        let json = JsonValue::parse(evidence.control_json.as_ref().unwrap()).unwrap();
        assert_eq!(json.str_field("schema"), Some(SCHEMA));
        assert_eq!(json.str_field("mode"), Some("assessment-only"));
        assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(false)));
        assert_eq!(json.f64_field("final_bound_k"), Some(measured(&evidence)));
        assert!(json.f64_field("response_iterations").unwrap() <= 100.0);
        // Exercise the existing actual QoI bridge: exactly one solver term,
        // bound to the conduction receipt rather than a made-up baseline.
        let terms = super::super::super::propagation_term_receipts(None, None,
            evidence.term.as_ref(), fs_blake3::hash_bytes(b"coupled-conduction-receipt")).unwrap();
        assert_eq!(terms.len(), 1);
        let again = fixture.assess(cx, &field, 64 * 1024 * 1024, 1.0).unwrap();
        assert_eq!(evidence.control_json, again.control_json);
        assert_eq!(measured(&evidence).to_bits(), measured(&again).to_bits());
        let near = fixture.assess(cx, &vec![330.0; field.len()], 64 * 1024 * 1024, 1e-6).unwrap();
        assert!(measured(&near) < 1e-6);
        assert!(measured(&near) < measured(&evidence));
        let json = JsonValue::parse(near.control_json.as_ref().unwrap()).unwrap();
        assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(true)));
    });
}

#[test]
fn coupled_product_gaps_are_not_replaced_by_frozen_or_tolerance_estimates() {
    let mut fixture = Fixture::new();
    with_cx(|_, cx| {
        let field = vec![300.0; fixture.mesh.vertex_count()];
        let refused = fixture.assess(cx, &field, 0, 1.0).unwrap();
        assert!(matches!(refused.term, Some(PropagatedTerm::Unmeasured { .. })));
        let json = JsonValue::parse(refused.control_json.as_ref().unwrap()).unwrap();
        assert_eq!(json.str_field("status"), Some("coupled-analysis-unavailable"));
        assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(false)));
    });
    fixture.paths = vec![AirPath::new(330.0, 1.0, 24.0,
        vec![AirSegment::new("air", 7.0, 2.0).unwrap()]).unwrap()];
    with_cx(|_, cx| {
        let refused = fixture.assess(cx, &vec![300.0; fixture.mesh.vertex_count()],
            64 * 1024 * 1024, 1.0).unwrap();
        assert!(matches!(refused.term, Some(PropagatedTerm::Unmeasured { .. })));
        assert!(refused.control_json.unwrap().contains("matching solid/air"));
    });
}

#[test]
fn coupled_product_does_not_invent_accuracy_or_publish_after_cancellation() {
    let fixture = Fixture::new();
    with_cx(|gate, cx| {
        let field = vec![330.0; fixture.mesh.vertex_count()];
        for request in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let evidence = fixture.assess(cx, &field, 64 * 1024 * 1024, request).unwrap();
            assert!(measured(&evidence).is_finite());
            let json = JsonValue::parse(evidence.control_json.as_ref().unwrap()).unwrap();
            assert_eq!(json.get("requested_tolerance_k"), Some(&JsonValue::Null));
            assert_eq!(json.get("goal_met"), Some(&JsonValue::Bool(false)));
        }
        gate.request();
        assert!(fixture.assess(cx, &field, 64 * 1024 * 1024, 1.0).is_err());
    });
}
