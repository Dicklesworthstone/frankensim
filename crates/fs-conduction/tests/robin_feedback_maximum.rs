//! Real FEM Robin feedback, with independent dense and constant-field oracles.
mod support;
use fs_conduction::adjoint::{LinearGoalAnalyzer, LinearGoalAnalysisConfig, RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBoundary, ThermalBoundaryBuilder, ThermalBc};
use fs_conduction::fixtures::{unit_cube, on_box_face};
use fs_solver::goal::{GoalResidualLimits};
use fs_solver::goal::feedback::{FeedbackBoundStatus, FeedbackResidualLimits};
use fs_sparse::Csr;
use support::{with_cx, with_cancelled_cx};

struct Fixture { mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel, source: ScalarField }
impl Fixture {
    fn new(fixed: bool) -> Self {
        let (complex, positions) = unit_cube(1);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let mut builder = ThermalBoundaryBuilder::new(&mesh);
        if fixed {
            builder = builder.region("cold", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        }
        let boundary = builder.region("air", |f| !fixed || !on_box_face(f.centroid[0], 0.0),
            ThermalBc::robin(2.0, 300.0).unwrap()).unwrap().finish().unwrap();
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(10.0).unwrap(), source: ScalarField::Uniform(0.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            element_materials: None, source: &self.source }
    }
    fn base(&self, cx: &fs_exec::Cx<'_>) -> LinearGoalAnalyzer<'_> {
        LinearGoalAnalyzer::new_for_maximum(cx, self.problem(), None,
            LinearConfig { tolerance: 1e-11, max_iterations: 100, restart: 20 },
            &vec![300.0; self.mesh.vertex_count()],
            LinearGoalAnalysisConfig { residual_limits: limits().residual.solid, max_stability_iterations: 100 }).unwrap()
    }
}
fn limits() -> RobinFeedbackAnalysisConfig {
    RobinFeedbackAnalysisConfig {
        residual: FeedbackResidualLimits {
            solid: GoalResidualLimits { max_rows: 100, max_nonzeros: 10_000 },
            max_ports: 4, max_transfer_nonzeros: 10_000,
            max_response_entries: 400, max_verification_entries: 100_000,
        },
        max_response_iterations: 100, max_lowering_entries: 100_000,
    }
}
fn solve_dense(a: &Csr, rhs: &[f64], b: &Csr, c: &Csr, d: &[f64]) -> Vec<f64> {
    let n = a.nrows(); let p = d.len();
    let mut rows = vec![vec![0.0; n+1]; n];
    for i in 0..n {
        rows[i][n] = rhs[i];
        for k in 0..p { rows[i][n] += b.get(i,k)*d[k]; }
        for j in 0..n {
            rows[i][j] = a.get(i,j);
            for k in 0..p { rows[i][j] -= b.get(i,k)*c.get(k,j); }
        }
    }
    for k in 0..n {
        let pivot = (k..n).max_by(|&i,&j| rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k,pivot); let diagonal=rows[k][k]; assert!(diagonal.abs()>1e-12);
        for j in k..=n { rows[k][j] /= diagonal; }
        for i in 0..n { if i!=k {
            let factor=rows[i][k];
            for j in k..=n { let source=rows[k][j]; rows[i][j] -= factor*source; }
        } }
    }
    rows.iter().map(|row| row[n]).collect()
}

#[test]
fn a_frozen_solid_solution_does_not_hide_coupled_maximum_error() {
    let fixture=Fixture::new(false);
    with_cx(|cx| {
        let temperature=vec![300.0; fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..temperature.len()).collect();
        let base=fixture.base(cx);
        let frozen=base.analyze_maximum(cx,&temperature,&vertices).unwrap();
        assert!(frozen.algebraic_half_width_k().unwrap()<1e-6);
        let coupled=base.with_robin_feedback(cx,&["air"],&[198.0],&[0.4],limits()).unwrap();
        let report=coupled.analyze_maximum(cx,&temperature,&vertices).unwrap();
        assert_eq!(report.coupled().status(),FeedbackBoundStatus::Enclosed);
        let (a,rhs,b,c,d)=coupled.stored_system();
        let exact=solve_dense(a,rhs,b,c,d);
        for &value in &exact { assert!((value-330.0).abs()<1e-8); }
        let maximum=exact.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        let band=report.interval_k().unwrap();
        assert!(band[0]<=maximum && maximum<=band[1]);
        assert!(report.algebraic_half_width_k().unwrap()>=29.999999);
        assert!(report.response_iterations()>0 && report.response_iterations()<=100);
        assert!(!report.meets_absolute_tolerance(1.0));
    });
}

#[test]
fn cached_feedback_assessment_replays_and_never_repeats_response_work() {
    let fixture=Fixture::new(false);
    with_cx(|cx| {
        let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],limits()).unwrap();
        let vertices:Vec<_>=(0..fixture.mesh.vertex_count()).collect();
        let temperature=vec![320.0; vertices.len()];
        let first=analyzer.analyze_maximum(cx,&temperature,&vertices).unwrap();
        let mut reversed=vertices.clone();reversed.reverse();
        assert_eq!(first,analyzer.analyze_maximum(cx,&temperature,&reversed).unwrap());
        let closer=analyzer.analyze_maximum(cx,&vec![329.0;vertices.len()],&vertices).unwrap();
        assert_eq!(first.response_iterations(),closer.response_iterations());
        assert!(closer.algebraic_half_width_k().unwrap()<first.algebraic_half_width_k().unwrap());
    });
}

#[test]
fn prescribed_wall_contributions_survive_reduction_and_fixed_goals_are_exact() {
    let fixture=Fixture::new(true);
    with_cx(|cx| {
        let base=fixture.base(cx);
        let fixed=base.dofs().fixed().to_vec();let dofs=base.dofs().free().to_vec();
        let analyzer=base.with_robin_feedback(cx,&["air"],&[198.0],&[0.4],limits()).unwrap();
        let (a,rhs,b,c,d)=analyzer.stored_system();
        let exact=solve_dense(a,rhs,b,c,d);
        let temperature=vec![300.0;fixture.mesh.vertex_count()];
        let report=analyzer.analyze_maximum(cx,&temperature,&dofs).unwrap();
        let bound=report.interval_k().unwrap();
        let maximum=exact.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        assert!(bound[0]<=maximum && maximum<=bound[1]);
        assert!(maximum>300.0);
        let fixed_report=analyzer.analyze_maximum(cx,&temperature,&fixed).unwrap();
        assert_eq!(fixed_report.interval_k(),Some([300.0,300.0]));
        assert_eq!(fixed_report.algebraic_half_width_k(),Some(0.0));
        let mut changed=temperature;changed[fixed[0]]=301.0;
        assert!(analyzer.analyze_maximum(cx,&changed,&fixed).is_err());
    });
}

#[test]
fn incomplete_response_work_cannot_be_mistaken_for_a_certificate() {
    let fixture=Fixture::new(false);
    with_cx(|cx| {
        for cap in [0,1,3] {
            let mut config=limits();config.max_response_iterations=cap;
            let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],config).unwrap();
            let temperature=vec![300.0;fixture.mesh.vertex_count()];
            let vertices:Vec<_>=(0..temperature.len()).collect();
            let r=analyzer.analyze_maximum(cx,&temperature,&vertices).unwrap();
            assert!(r.response_iterations()<=cap);
            if let Some(band)=r.interval_k() { assert!(band[0]<=330.0 && 330.0<=band[1]); }
            else { assert!(!r.meets_absolute_tolerance(1e9)); }
        }
        let a=fixture.base(cx).with_robin_feedback(cx,&["air"],&[0.0],&[2.0],limits()).unwrap();
        let temperature=vec![300.0;fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..temperature.len()).collect();
        let r=a.analyze_maximum(cx,&temperature,&vertices).unwrap();
        assert_eq!(r.coupled().status(),FeedbackBoundStatus::ContractionNotEstablished);
        assert!(r.interval_k().is_none());
    });
}

#[test]
fn malformed_feedback_limits_and_cancellation_refuse_without_publication() {
    let fixture=Fixture::new(false);
    with_cx(|cx| {
        assert!(fixture.base(cx).with_robin_feedback(cx,&["missing"],&[0.0],&[0.5],limits()).is_err());
        assert!(fixture.base(cx).with_robin_feedback(cx,&["air","air"],&[0.0;2],&[0.0;4],limits()).is_err());
        assert!(fixture.base(cx).with_robin_feedback(cx,&["air"],&[f64::NAN],&[0.5],limits()).is_err());
        let mut config=limits();config.residual.max_transfer_nonzeros=0;
        assert!(fixture.base(cx).with_robin_feedback(cx,&["air"],&[0.0],&[0.5],config).is_err());
        let mut config=limits();config.max_lowering_entries=0;
        assert!(fixture.base(cx).with_robin_feedback(cx,&["air"],&[0.0],&[0.5],config).is_err());
        let prepared=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],limits()).unwrap();
        let temperature=vec![300.0;fixture.mesh.vertex_count()];
        with_cancelled_cx(|cancelled| {
            assert!(matches!(prepared.analyze_maximum(cancelled,&temperature,&[0]),Err(fs_conduction::ConductionError::Cancelled{..})));
        });
        assert!(prepared.analyze_maximum(cx,&temperature,&[0]).unwrap().interval_k().is_some());
    });
}
