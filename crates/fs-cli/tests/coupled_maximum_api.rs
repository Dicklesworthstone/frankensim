//! Public solid/air consumer path, with no mocked solver or precomputed output.
use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::{affine_reference_law, prepare_linear_maximum};
use fs_conduction::adjoint::{LinearGoalAnalyzer, LinearGoalAnalysisConfig,
    LinearRobinFeedbackAnalyzer, RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBoundary, ThermalBoundaryBuilder, ThermalBc};
use fs_solver::goal::{GoalResidualLimits};
use fs_solver::goal::feedback::{FeedbackBoundStatus, FeedbackResidualLimits};
use fs_alloc::{ArenaPool, ArenaConfig};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 7304, kernel_id: 73, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn path(inlet:f64, capacity:f64, rows:&[(&str,f64,f64)]) -> AirPath {
    AirPath::new(inlet,1.0,capacity,rows.iter().map(|(name,area,h)|
        AirSegment::new(name,*area,*h).unwrap()).collect()).unwrap()
}

#[test]
fn affine_air_law_matches_real_marches_and_retains_all_upstream_dependencies() {
    with_cx(|_,cx| {
        let paths=[path(290.0,8.0,&[("a",1.0,2.0),("b",2.0,3.0),("c",1.0,1.0)]),
            path(310.0,5.0,&[("d",1.0,4.0),("e",2.0,1.0)])];
        let walls=[350.0,330.0,320.0,310.0,370.0];
        let law=affine_reference_law(cx,&paths,5,25).unwrap();
        let got=law.evaluate(cx,&walls).unwrap();
        let mut start=0;
        for p in &paths {
            let end=start+p.segments().len();
            let marched=p.march(&walls[start..end]).unwrap();
            for (actual,step) in got[start..end].iter().zip(&marched.segments) {
                assert!((*actual-step.reference_temperature_k).abs()<1e-10);
            }
            start=end;
        }
        assert_eq!(law.regions(),&["a","b","c","d","e"]);
        let d=law.wall_matrix();
        assert!(d[2*5]>0.0 && d[2*5+1]>0.0); // Every upstream wall matters.
        for i in 0..5 { for j in 0..5 {
            if (i<3)!=(j<3) || j>i { assert_eq!(d[i*5+j],0.0); }
        } }
        let changed=[370.0,330.0,320.0,310.0,370.0];
        let next=law.evaluate(cx,&changed).unwrap();
        assert!(next[2]>got[2]);
        assert_eq!(&next[3..],&got[3..]); // Independent inlets stay independent.
    });
}

#[test]
fn air_reference_limits_cancellation_and_extreme_ntu_are_explicit() {
    with_cx(|gate,cx| {
        for ntu in [1e-12,1e-6,0.5,1000.0] {
            let p=path(300.0,1.0,&[("air",1.0,ntu)]);
            let law=affine_reference_law(cx,std::slice::from_ref(&p),1,1).unwrap();
            let got=law.evaluate(cx,&[360.0]).unwrap()[0];
            let want=p.march(&[360.0]).unwrap().segments[0].reference_temperature_k;
            assert!((got-want).abs()<1e-10);
            assert!(law.evaluate(cx,&[]).is_err());
            assert!(law.evaluate(cx,&[f64::NAN]).is_err());
        }
        let p=path(300.0,1.0,&[("air",1.0,1.0)]);
        assert!(affine_reference_law(cx,&[],1,1).is_err());
        assert!(affine_reference_law(cx,std::slice::from_ref(&p),0,1).is_err());
        assert!(affine_reference_law(cx,std::slice::from_ref(&p),1,0).is_err());
        assert!(affine_reference_law(cx,&[p.clone(),p.clone()],2,4).is_err());
        gate.request();
        assert!(affine_reference_law(cx,&[p],1,1).is_err());
    });
}

struct Fixture { mesh:ConductionMesh,boundary:ThermalBoundary,material:ConductivityModel,source:ScalarField }
impl Fixture {
    fn new() -> Self {
        let (complex,positions)=fs_conduction::fixtures::unit_cube(1);
        let mesh=ConductionMesh::new(complex,positions).unwrap();
        let boundary=ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_|true,ThermalBc::robin(2.0,300.0).unwrap()).unwrap().finish().unwrap();
        Self { mesh,boundary,material:ConductivityModel::isotropic_declared(10.0).unwrap(),source:ScalarField::Uniform(0.0) }
    }
    fn problem(&self)->ConductionProblem<'_> {
        ConductionProblem {mesh:&self.mesh,boundary:&self.boundary,material:&self.material,
            element_materials:None,source:&self.source}
    }
    fn path(&self,inlet:f64)->AirPath {
        let area=self.mesh.boundary().iter().map(|face|face.area).sum();
        path(inlet,24.0,&[("air",area,2.0)])
    }
}
fn configs()->(LinearConfig,LinearGoalAnalysisConfig,RobinFeedbackAnalysisConfig) {
    let solid=GoalResidualLimits {max_rows:100,max_nonzeros:10_000};
    (LinearConfig {tolerance:1e-11,max_iterations:100,restart:20},
        LinearGoalAnalysisConfig {residual_limits:solid,max_stability_iterations:100},
        RobinFeedbackAnalysisConfig {residual:FeedbackResidualLimits {solid,max_ports:4,
            max_transfer_nonzeros:10_000,max_response_entries:400,max_verification_entries:100_000},
            max_response_iterations:100,max_lowering_entries:100_000})
}
fn dense_maximum(analyzer:&LinearRobinFeedbackAnalyzer<'_>)->f64 {
    let (a,rhs,b,c,d)=analyzer.stored_system();
    let n=a.nrows();let mut rows=vec![vec![0.0;n+1];n];
    for i in 0..n {
        rows[i][n]=rhs[i];
        for k in 0..d.len() {rows[i][n]+=b.get(i,k)*d[k];}
        for j in 0..n {rows[i][j]=a.get(i,j);
            for k in 0..d.len() {rows[i][j]-=b.get(i,k)*c.get(k,j);}
        }
    }
    for k in 0..n {
        let pivot=(k..n).max_by(|&i,&j|rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k,pivot);let diagonal=rows[k][k];assert!(diagonal.abs()>1e-12);
        for j in k..=n {rows[k][j]/=diagonal;}
        for i in 0..n {if i!=k {let factor=rows[i][k];
            for j in k..=n {let value=rows[k][j];rows[i][j]-=factor*value;}
        }}
    }
    rows.iter().map(|row|row[n]).fold(f64::NEG_INFINITY,f64::max)
}

#[test]
fn real_air_feedback_covers_the_error_a_frozen_solid_would_miss_and_replays() {
    let fixture=Fixture::new();
    with_cx(|_,cx| {
        let (linear,solid,feedback)=configs();
        let temperature=vec![300.0;fixture.mesh.vertex_count()];
        let vertices:Vec<_>=(0..temperature.len()).collect();
        let frozen=LinearGoalAnalyzer::new_for_maximum(cx,fixture.problem(),None,linear,&temperature,solid).unwrap();
        assert!(frozen.analyze_maximum(cx,&temperature,&vertices).unwrap().algebraic_half_width_k().unwrap()<1e-6);
        let paths=[fixture.path(330.0)];
        let analyzer=prepare_linear_maximum(cx,fixture.problem(),None,&paths,linear,&temperature,solid,feedback).unwrap();
        let actual=analyzer.analyze_maximum(cx,&temperature,&vertices).unwrap();
        assert_eq!(actual.coupled().status(),FeedbackBoundStatus::Enclosed);
        let want=dense_maximum(&analyzer);assert!((want-330.0).abs()<1e-8);
        let band=actual.interval_k().unwrap();assert!(band[0]<=want && want<=band[1]);
        assert!(!actual.meets_absolute_tolerance(1.0));
        assert_eq!(actual,analyzer.analyze_maximum(cx,&temperature,&vertices).unwrap());
        let nearer=analyzer.analyze_maximum(cx,&vec![329.0;temperature.len()],&vertices).unwrap();
        assert!(nearer.algebraic_half_width_k().unwrap()<actual.algebraic_half_width_k().unwrap());
        assert_eq!(nearer.response_iterations(),actual.response_iterations());
        assert!(actual.response_iterations()<=feedback.max_response_iterations);
    });
}

#[test]
fn mismatched_physics_and_nonlinear_material_cannot_silently_freeze_into_a_bound() {
    let mut fixture=Fixture::new();
    with_cx(|gate,cx| {
        let (linear,solid,feedback)=configs();
        let temperature=vec![300.0;fixture.mesh.vertex_count()];
        for p in [path(330.0,24.0,&[("air",6.0,3.0)]),path(330.0,24.0,&[("air",7.0,2.0)]),
            path(330.0,24.0,&[("missing",6.0,2.0)])] {
            assert!(prepare_linear_maximum(cx,fixture.problem(),None,&[p],linear,&temperature,solid,feedback).is_err());
        }
        let paths=[fixture.path(330.0)];
        let prepared=prepare_linear_maximum(cx,fixture.problem(),None,&paths,linear,&temperature,solid,feedback).unwrap();
        gate.request();
        assert!(prepared.analyze_maximum(cx,&temperature,&[0]).is_err());
    });
    fixture.material=ConductivityModel::isotropic(fs_conduction::material::ConductivityTable::declared_curve(
        vec![(290.0,9.0),(350.0,11.0)]).unwrap());
    with_cx(|_,cx| {
        let (linear,solid,feedback)=configs();
        assert!(prepare_linear_maximum(cx,fixture.problem(),None,&[fixture.path(330.0)],linear,
            &vec![300.0;fixture.mesh.vertex_count()],solid,feedback).is_err());
    });
}
