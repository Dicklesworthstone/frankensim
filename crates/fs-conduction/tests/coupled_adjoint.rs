//! G0/G3/G4: real FEM adjoints of A-B*C, checked against dense solves and
//! parameter perturbations. No frozen-solid or solver-iteration derivative.
mod support;

use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer,
    LinearRobinFeedbackAnalyzer, RobinFeedbackAnalysisConfig};
use fs_conduction::assemble::DofMap;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};
use support::{with_cx, with_gate};

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-10, max_iterations: 300, restart: 16 }
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary,
    material: ConductivityModel, source: ScalarField,
}
impl Fixture {
    fn new(fixed: bool, h: f64) -> Self {
        let (complex, positions) = fs_conduction::fixtures::unit_cube(2);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let mut b = ThermalBoundaryBuilder::new(&mesh);
        if fixed {
            b = b.region("fixed", |f| f.centroid[0].abs() < 1e-10,
                ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        }
        let boundary = b.region("first", |f| f.centroid[1] < 0.5
                && (!fixed || f.centroid[0] > 1e-10), ThermalBc::robin(h, 300.0).unwrap()).unwrap()
            .region("second", |f| f.centroid[1] >= 0.5
                && (!fixed || f.centroid[0] > 1e-10), ThermalBc::robin(4.0, 300.0).unwrap()).unwrap()
            .finish().unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 3.0+p[0]+2.0*p[1]).collect());
        Self { mesh, boundary, source, material: ConductivityModel::isotropic_declared(10.0).unwrap() }
    }
    fn analyzer<'a>(&'a self, cx: &fs_exec::Cx<'_>, offset: f64, mixing: [f64;4])
        -> LinearRobinFeedbackAnalyzer<'a>
    {
        let limits = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
        LinearGoalAnalyzer::new_for_maximum(cx, ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            source: &self.source, element_materials: None,
        }, None, linear(), &vec![300.0;self.mesh.vertex_count()],
            LinearGoalAnalysisConfig { residual_limits: limits, max_stability_iterations: 300 })
            .unwrap().with_robin_feedback(cx, &["first", "second"], &[250.0+offset, 180.0],
                &mixing, RobinFeedbackAnalysisConfig { residual: FeedbackResidualLimits {
                    solid: limits, max_ports: 4, max_transfer_nonzeros: 10_000,
                    max_response_entries: 1000, max_verification_entries: 10_000_000,
                }, max_response_iterations: 300, max_lowering_entries: 100_000 }).unwrap()
    }
    fn dofs(&self) -> DofMap { DofMap::new(&self.boundary, self.mesh.vertex_count()).unwrap() }
}
fn dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    for j in 0..b.len() {
        let p = (j..b.len()).max_by(|&i,&k| a[i][j].abs().total_cmp(&a[k][j].abs())).unwrap();
        a.swap(p,j); b.swap(p,j);
        let diagonal=a[j][j]; assert!(diagonal.abs()>1e-12);
        for k in j..b.len() { a[j][k]/=diagonal; } b[j]/=diagonal;
        for i in 0..b.len() { if i!=j {
            let factor=a[i][j];
            for k in j..b.len() { let v=a[j][k]; a[i][k]-=factor*v; }
            let v=b[j]; b[i]-=factor*v;
        } }
    }
    b
}
fn operator(a: &LinearRobinFeedbackAnalyzer<'_>) -> Vec<Vec<f64>> {
    let (m,_,b,c,d)=a.stored_system();
    (0..m.nrows()).map(|i| (0..m.nrows()).map(|j|
        m.get(i,j)-(0..d.len()).map(|k| b.get(i,k)*c.get(k,j)).sum::<f64>()).collect()).collect()
}
fn primal(f: &Fixture, a: &LinearRobinFeedbackAnalyzer<'_>) -> Vec<f64> {
    let (_,rhs,b,_,d)=a.stored_system();
    let rhs=(0..rhs.len()).map(|i| rhs[i]+(0..d.len()).map(|j| b.get(i,j)*d[j]).sum::<f64>()).collect();
    let free=dense(operator(a),rhs); let dofs=f.dofs(); let mut full=dofs.prescribed().to_vec();
    for (&v,&t) in dofs.free().iter().zip(&free) { full[v]=t; } full
}
fn dot(a: &[f64],b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
const MIXING: [f64;4] = [0.2, 0.0, 0.3, 0.1];

#[test]
fn full_transpose_and_control_gradients_match_independent_dense_system() {
    for fixed in [false,true] {
        let f=Fixture::new(fixed,2.0);
        with_cx(|cx| {
            let a=f.analyzer(cx,0.0,MIXING); let temperature=primal(&f,&a);
            let weights: Vec<_>=(0..temperature.len()).map(|i| (i as f64*0.73).sin()).collect();
            let gradient=a.pullback_affine_controls(cx,&temperature,&weights,linear()).unwrap();
            let m=operator(&a); let dofs=f.dofs();
            let transposed=(0..m.len()).map(|i| (0..m.len()).map(|j| m[j][i]).collect()).collect();
            let expected=dense(transposed,dofs.gather(&weights));
            for (&v,&value) in dofs.free().iter().zip(&expected) {
                assert!((gradient.nodal_load[v]-value).abs()<1e-8);
            }
            for &v in dofs.fixed() { assert_eq!(gradient.nodal_load[v],0.0); }
            let (_,_,b,_,_)=a.stored_system();
            for j in 0..2 {
                let expected=(0..expected.len()).map(|i| b.get(i,j)*expected[i]).sum::<f64>();
                assert!((gradient.references[j]-expected).abs()<1e-8);
            }
            assert!(gradient.relative_residual<linear().tolerance);
            assert!(gradient.iterations>0 && gradient.iterations<=linear().max_iterations);
            let frozen=dense((0..m.len()).map(|j| (0..m.len())
                .map(|i| a.stored_system().0.get(i,j)).collect()).collect(),dofs.gather(&weights));
            assert!(expected.iter().zip(frozen).any(|(a,b)| (a-b).abs()>1e-5));
        });
    }
}

#[test]
fn offsets_and_log_coefficient_match_complete_perturbed_solves() {
    let f=Fixture::new(true,2.0);
    with_cx(|cx| {
        let a=f.analyzer(cx,0.0,MIXING); let temperature=primal(&f,&a);
        let weights: Vec<_>=f.mesh.positions().iter().map(|p| 1.0+p[2]).collect();
        let g=a.pullback_affine_controls(cx,&temperature,&weights,linear()).unwrap();
        let step=1e-4;
        let plus=primal(&f,&f.analyzer(cx,step,MIXING));
        let minus=primal(&f,&f.analyzer(cx,-step,MIXING));
        assert!(((dot(&weights,&plus)-dot(&weights,&minus))/(2.0*step)-g.references[0]).abs()<1e-6);
        let plus_f=Fixture::new(true,2.0*step.exp());
        let minus_f=Fixture::new(true,2.0*(-step).exp());
        let plus=primal(&plus_f,&plus_f.analyzer(cx,0.0,MIXING));
        let minus=primal(&minus_f,&minus_f.analyzer(cx,0.0,MIXING));
        assert!(((dot(&weights,&plus)-dot(&weights,&minus))/(2.0*step)-g.log_htc[0]).abs()<1e-5);
        // A pullback neither mutates the immutable preparation nor spends a new response solve.
        let prepared=a.response_iterations();
        let repeated=a.pullback_affine_controls(cx,&temperature,&weights,linear()).unwrap();
        assert_eq!(g,repeated); assert_eq!(a.response_iterations(),prepared);
    });
}

#[test]
fn signed_scaled_and_prescribed_only_functionals_keep_their_semantics() {
    let f=Fixture::new(true,2.0);
    with_cx(|cx| {
        let a=f.analyzer(cx,0.0,MIXING); let temperature=primal(&f,&a);
        let weights=vec![1.0;temperature.len()];
        let g=a.pullback_affine_controls(cx,&temperature,&weights,linear()).unwrap();
        for scale in [-1e-120,2.0,1e120] {
            let weights=vec![scale;temperature.len()];
            let scaled=a.pullback_affine_controls(cx,&temperature,&weights,linear()).unwrap();
            for (x,y) in g.nodal_load.iter().zip(&scaled.nodal_load) {
                assert!((x-y/scale).abs()<1e-9);
            }
        }
        let mut fixed=vec![0.0;temperature.len()];
        for &v in f.dofs().fixed() { fixed[v]=7.0; }
        let zero=a.pullback_affine_controls(cx,&temperature,&fixed,linear()).unwrap();
        assert_eq!(zero.iterations,0); assert_eq!(zero.relative_residual,0.0);
        assert!(zero.nodal_load.iter().chain(&zero.references).chain(&zero.log_htc).all(|&v| v==0.0));
    });
}

#[test]
fn invalid_inputs_exhaustion_and_cancellation_never_publish_a_gradient() {
    let f=Fixture::new(true,2.0);
    with_cx(|cx| {
        let a=f.analyzer(cx,0.0,MIXING); let temperature=primal(&f,&a);
        let weights: Vec<_>=(0..temperature.len()).map(|i| (i as f64).cos()).collect();
        assert!(matches!(a.pullback_affine_controls(cx,&temperature,&weights,
            LinearConfig { max_iterations:1, restart:1, tolerance:1e-14 }),
            Err(fs_conduction::ConductionError::LinearSolveFailed { krylov_iterations:1,.. })));
        assert!(a.pullback_affine_controls(cx,&temperature,&weights[..2],linear()).is_err());
        let mut malformed=weights.clone();malformed[0]=f64::NAN;
        assert!(a.pullback_affine_controls(cx,&temperature,&malformed,linear()).is_err());
        assert!(a.pullback_affine_controls(cx,&temperature,&weights,
            LinearConfig { restart:0,..linear() }).is_err());
        let mut changed=temperature.clone();changed[f.dofs().fixed()[0]]=301.0;
        assert!(a.pullback_affine_controls(cx,&changed,&weights,linear()).is_err());
        with_gate(|gate,cancelled| {
            gate.request();
            assert!(matches!(a.pullback_affine_controls(cancelled,&temperature,&weights,linear()),
                Err(fs_conduction::ConductionError::Cancelled {..})));
        });
        assert!(a.pullback_affine_controls(cx,&temperature,&weights,linear()).is_ok());
    });
}

#[test]
fn compatible_singular_feedback_does_not_authorize_an_adjoint() {
    let f=Fixture::new(false,2.0);
    with_cx(|cx| {
        let a=f.analyzer(cx,0.0,[1.0,0.0,0.0,1.0]);
        let mut weights=vec![0.0;f.mesh.vertex_count()];weights[0]=1.0;weights[1]=-1.0;
        assert!(a.pullback_affine_controls(cx,&vec![300.0;weights.len()],&weights,linear()).is_err());
    });
}
