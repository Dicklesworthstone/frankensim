//! Independent directional and augmented-Lagrangian identities on the live IR.
use fs_opt::{AugmentedLagrangian, ConstraintKind, Manifold, OptError, Problem,
    ProblemBuilder, ReverseProblem, ReverseProblemError, Sense};
use fs_opt::reverse::{ReverseError, ReverseLimits, ReverseProgram};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian,
    PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_qty::Dims;
use std::cell::Cell;

fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 10_000, max_scalar_slots: 20_000 } }
fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        assert!((a-b).abs() <= tolerance*(1.0+b.abs()), "{a} != {b}");
    }
}
fn fixture() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let f = b.norm_sq(r).unwrap(); b.objective(f, Sense::Minimize, 1.0).unwrap();
    let xy = b.mul(x,y).unwrap(); let one = b.konst(1.0,Dims::NONE).unwrap();
    let c = b.sub(xy,one).unwrap(); b.constraint(c,ConstraintKind::EqZero,"product").unwrap();
    let four = b.konst(4.0,Dims::NONE).unwrap(); let c = b.sub(f,four).unwrap();
    b.constraint(c,ConstraintKind::LeZero,"ball").unwrap(); b.finish()
}

#[test]
fn nonlinear_active_and_inactive_constraints_have_complete_hessians() {
    let problem = fixture(); let oracle = ReverseProblem::new(&problem,limits()).unwrap();
    let lambda = [-0.3,0.4]; let mu = 2.0;
    let model = AugmentedLagrangian::new(&problem,&lambda,mu).unwrap();
    for point in [[0.5,0.25],[1.4,1.2],[-2.0,0.5]] {
        let [x,y] = point; let v = [0.7,-0.2];
        let c = [x*y-1.0,x*x+y*y-4.0];
        let w = [lambda[0]+mu*c[0],(lambda[1]+mu*c[1]).max(0.0)];
        let on = w[1] > 0.0;
        let jv = [y*v[0]+x*v[1],2.0*x*v[0]+2.0*y*v[1]];
        let expected = [(2.0+2.0*w[1])*v[0]+w[0]*v[1]+mu*jv[0]*y+if on {mu*jv[1]*2.0*x} else {0.0},
            (2.0+2.0*w[1])*v[1]+w[0]*v[0]+mu*jv[0]*x+if on {mu*jv[1]*2.0*y} else {0.0}];
        let tape = oracle.evaluate(&[point.to_vec()]).unwrap();
        close(&tape.constraint_jacobian_vector_product(&[v.to_vec()],None).unwrap(),&jv,1e-14);
        close(&model.hessian_vector_product(&tape,&[v.to_vec()],None).unwrap()[0],&expected,1e-13);
        let (value,g) = model.value_gradient(&tape,None).unwrap();
        let expected_value = x*x+y*y+lambda[0]*c[0]+0.5*mu*c[0]*c[0]
            +(w[1]*w[1]-lambda[1]*lambda[1])/(2.0*mu);
        close(&[value],&[expected_value],1e-14);
        close(&g[0],&[2.0*x+w[0]*y+w[1]*2.0*x,2.0*y+w[0]*x+w[1]*2.0*y],1e-14);
        assert_eq!(model.updated_multipliers(&tape,None).unwrap(),w);
        let h = 1e-5;
        let plus = oracle.evaluate(&[vec![x+h*v[0],y+h*v[1]]]).unwrap();
        let minus = oracle.evaluate(&[vec![x-h*v[0],y-h*v[1]]]).unwrap();
        let gp = model.value_gradient(&plus,None).unwrap().1;
        let gm = model.value_gradient(&minus,None).unwrap().1;
        close(&gp[0].iter().zip(&gm[0]).map(|(p,m)|(p-m)/(2.0*h)).collect::<Vec<_>>(),&expected,1e-8);
    }
}

#[test]
fn directional_sweep_reuses_shared_vector_primitives_and_duplicate_roots() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x",Manifold::Rn {dim:3},Dims::NONE).unwrap();
    let s = b.var("s",Manifold::Rn {dim:1},Dims::NONE).unwrap();
    let xr = b.var_ref(x).unwrap(); let sr = b.var_ref(s).unwrap();
    let s = b.component(sr,0).unwrap(); let sx = b.mul(s,xr).unwrap();
    let dot = b.dot(sx,xr).unwrap(); let norm = b.norm_sq(sx).unwrap();
    let root = b.add(dot,norm).unwrap(); let p = b.finish();
    let program = ReverseProgram::new(&p,&[root,root],limits()).unwrap();
    let tape = program.evaluate(&[vec![1.0,2.0,3.0],vec![2.0]],None).unwrap();
    let direction = [vec![3.0,-1.0,2.0],vec![0.5]];
    // d[(s+s^2)|x|^2] = (1+2s)*ds*14 + 2*(s+s^2)*(x dot dx).
    assert_eq!(tape.jacobian_vector_product(&direction,None).unwrap(),vec![119.0,119.0]);
    let jt = tape.pullback(&[2.0,-0.5],None).unwrap();
    let paired: f64 = jt.iter().zip(&direction).flat_map(|(a,b)|a.iter().zip(b)).map(|(a,b)|a*b).sum();
    assert_eq!(paired,1.5*119.0);
}

#[test]
fn inactive_singular_constraint_does_not_poison_augmented_products() {
    let mut b = ProblemBuilder::new();
    let v=b.var("x",Manifold::Rn {dim:1},Dims::NONE).unwrap();let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();
    let square=b.powi(x,2).unwrap();b.objective(square,Sense::Minimize,1.0).unwrap();
    let sqrt=b.sqrt(x).unwrap();let one=b.konst(1.0,Dims::NONE).unwrap();let c=b.sub(sqrt,one).unwrap();
    b.constraint(c,ConstraintKind::LeZero,"inactive sqrt").unwrap();let p=b.finish();
    let oracle=ReverseProblem::new(&p,limits()).unwrap();let tape=oracle.evaluate(&[vec![0.0]]).unwrap();
    let model=AugmentedLagrangian::new(&p,&[0.0],1.0).unwrap();
    assert_eq!(model.value_gradient(&tape,None).unwrap(),(0.0,vec![vec![0.0]]));
    assert_eq!(model.hessian_vector_product(&tape,&[vec![1.0]],None).unwrap(),vec![vec![2.0]]);
    assert!(tape.constraint_jacobian_vector_product(&[vec![1.0]],None).is_err());
}

#[derive(Debug,Default)]
struct Physical { primals:Cell<usize>, products:Cell<usize> }
impl PhysicsModel for Physical {
    fn signature(&self)->PhysicsSignature<'_> {PhysicsSignature {study:"q",parameter_count:2,parameter_dims:Dims::NONE,residual_dims:Dims::NONE}}
    fn value_gradient(&self,x:&[f64],_:Option<&fs_exec::Cx<'_>>)->Result<PhysicsSample,PhysicsError> {
        self.primals.set(self.primals.get()+1);
        Ok(PhysicsSample {value:x[0]*x[0]+3.0*x[1],gradient:vec![2.0*x[0],3.0]})
    }
    fn second_order(&self)->Option<&dyn PhysicsHessian> {Some(self)}
}
impl PhysicsHessian for Physical {
    fn hessian_vector_product(&self,_:&[f64],v:&[f64],_:Option<&fs_exec::Cx<'_>>)->Result<Vec<f64>,PhysicsError> {
        self.products.set(self.products.get()+1);Ok(vec![2.0*v[0],0.0])
    }
}
#[test]
fn physical_constraints_use_one_hessian_action_and_no_extra_primal() {
    let mut b=ProblemBuilder::new();let x=b.var("x",Manifold::Rn {dim:2},Dims::NONE).unwrap();
    let r=b.var_ref(x).unwrap();let f=b.norm_sq(r).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let node=b.pde_residual("q",x,true,Dims::NONE).unwrap();b.constraint(node,ConstraintKind::EqZero,"physical").unwrap();
    let p=b.finish();let physics=Physical::default();let oracle=ReverseProblem::new_with_physics(&p,limits(),&[PhysicsBinding {node,model:&physics}],None).unwrap();
    let tape=oracle.evaluate(&[vec![2.0,-1.0]]).unwrap();let v=[vec![1.0,2.0]];
    assert_eq!(tape.constraint_jacobian_vector_product(&v,None).unwrap(),vec![10.0]);
    assert_eq!(physics.products.get(),0);
    let model=AugmentedLagrangian::new(&p,&[1.0],2.0).unwrap();
    assert_eq!(model.hessian_vector_product(&tape,&v,None).unwrap(),vec![vec![88.0,64.0]]);
    assert_eq!(physics.products.get(),1);assert_eq!(physics.primals.get(),1);
}

#[test]
fn shifted_inequality_zero_uses_the_documented_inactive_generalized_action() {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn {dim:1},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let f=b.powi(x,2).unwrap();
    b.objective(f,Sense::Minimize,1.0).unwrap();b.constraint(x,ConstraintKind::LeZero,"cap").unwrap();let p=b.finish();
    let oracle=ReverseProblem::new(&p,limits()).unwrap();let model=AugmentedLagrangian::new(&p,&[2.0],4.0).unwrap();
    let tape=oracle.evaluate(&[vec![-0.5]]).unwrap();
    assert_eq!(model.value_gradient(&tape,None).unwrap(),(-0.25,vec![vec![-1.0]]));
    assert_eq!(model.hessian_vector_product(&tape,&[vec![1.0]],None).unwrap(),vec![vec![2.0]]);
}

#[test]
fn admission_identity_direction_and_cancellation_fail_without_mutating_tapes() {
    let p=fixture();let oracle=ReverseProblem::new(&p,limits()).unwrap();let tape=oracle.evaluate(&[vec![1.0,1.0]]).unwrap();
    for mu in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(AugmentedLagrangian::new(&p,&[0.0,0.0],mu).is_err());}
    for w in [vec![],vec![0.0,-1.0],vec![f64::NAN,0.0]] {assert!(AugmentedLagrangian::new(&p,&w,1.0).is_err());}
    let copy=fixture();let foreign=AugmentedLagrangian::new(&copy,&[0.0,0.0],1.0).unwrap();
    assert!(matches!(foreign.value_gradient(&tape,None),Err(ReverseProblemError::Augmentation {..})));
    let model=AugmentedLagrangian::new(&p,&[0.0,0.0],1.0).unwrap();
    for direction in [vec![],vec![vec![1.0]],vec![vec![f64::INFINITY,1.0]]] {assert!(model.hessian_vector_product(&tape,&direction,None).is_err());}
    let gate=fs_exec::CancelGate::new();gate.request();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=fs_exec::Cx::new(&gate,arena,fs_exec::StreamKey {seed:0,kernel_id:1,tile:0,iteration:0},fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        assert!(matches!(tape.constraint_jacobian_vector_product(&[vec![1.0,1.0]],Some(&cx)),Err(ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::Cancelled)))));
        assert!(model.value_gradient(&tape,Some(&cx)).is_err());
    });
    assert_eq!(tape.objective_value(),2.0);
    // Dimensional residuals need explicit normalization, even with a zero multiplier.
    let mut b=ProblemBuilder::new();let f=b.konst(1.0,Dims::NONE).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let c=b.konst(1.0,Dims([1,0,0,0,0,0])).unwrap();b.constraint(c,ConstraintKind::EqZero,"meters").unwrap();
    assert!(AugmentedLagrangian::new(&b.finish(),&[0.0],1.0).is_err());
}

#[test]
fn wide_constraint_fanout_needs_only_linear_tape_storage() {
    let mut b=ProblemBuilder::new();let x=b.var("x",Manifold::Rn {dim:4096},Dims::NONE).unwrap();let r=b.var_ref(x).unwrap();
    let q=b.norm_sq(r).unwrap();b.objective(q,Sense::Minimize,1.0).unwrap();
    for i in 0..128 {b.constraint(q,ConstraintKind::EqZero,&format!("copy {i}")).unwrap();}
    let p=b.finish();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let tape=oracle.evaluate(&[vec![0.25;4096]]).unwrap();
    let model=AugmentedLagrangian::new(&p,&vec![0.0;128],1.0).unwrap();
    let hv=model.hessian_vector_product(&tape,&[vec![1.0;4096]],None).unwrap();
    assert_eq!(hv,vec![vec![196610.0;4096]]);
}
