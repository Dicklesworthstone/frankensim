//! Exact outer/inner chain rules across explicitly bound solved physics.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel,
    PhysicsSample, PhysicsSignature};
use fs_opt::reverse::{HessianError, ReverseError, ReverseLimits, ReverseProgram};
use fs_opt::{ConstraintKind, Manifold, NodeId, OptError, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::cell::{Cell, RefCell};

const LIMITS: ReverseLimits = ReverseLimits { max_nodes: 128, max_scalar_slots: 1024 };

#[derive(Debug, Default)]
struct Model<'a> {
    primals: Cell<usize>,
    products: Cell<usize>,
    error: RefCell<Option<PhysicsError>>,
    malformed: Cell<u8>,
    gate: Option<&'a CancelGate>,
    points: RefCell<Vec<Vec<f64>>>,
}
impl PhysicsModel for Model<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "q", parameter_count: 2,
            parameter_dims: Dims::NONE, residual_dims: Dims::NONE }
    }
    fn value_gradient(&self, x: &[f64], _: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.primals.set(self.primals.get()+1);
        Ok(PhysicsSample { value: x[0]*x[0]+3.0*x[1], gradient: vec![2.0*x[0],3.0] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Model<'_> {
    fn hessian_vector_product(&self, x: &[f64], v: &[f64], _: Option<&Cx<'_>>)
        -> Result<Vec<f64>, PhysicsError>
    {
        self.products.set(self.products.get()+1);
        self.points.borrow_mut().push(x.to_vec());
        if let Some(error) = self.error.borrow().as_ref() { return Err(error.clone()); }
        let mut product = vec![2.0*v[0],0.0];
        match self.malformed.get() {
            1 => { product.pop(); }
            2 => product.push(0.0),
            3 => product[1] = f64::NAN,
            4 => product[0] = f64::INFINITY,
            _ => {}
        }
        if let Some(gate) = self.gate { gate.request(); }
        Ok(product)
    }
}
fn fixture() -> (Problem, NodeId) {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let xr = b.var_ref(x).unwrap();
    let norm = b.norm_sq(xr).unwrap();
    let pde = b.pde_residual("q", x, true, Dims::NONE).unwrap();
    let square = b.mul(pde,pde).unwrap();
    let half = b.konst(0.5,Dims::NONE).unwrap();
    let reg = b.mul(half,norm).unwrap();
    let total = b.add(square,reg).unwrap();
    let two = b.konst(2.0,Dims::NONE).unwrap();
    let residual = b.sub(pde,two).unwrap();
    b.objective(total,Sense::Minimize,1.0).unwrap();
    b.constraint(residual,ConstraintKind::LeZero,"q cap").unwrap();
    (b.finish(),pde)
}

#[test]
fn composed_objective_and_lagrangian_use_both_chain_terms_and_retained_point() {
    let (problem,node) = fixture();
    let model = Model::default();
    let oracle = ReverseProblem::new_with_physics(&problem,LIMITS,
        &[PhysicsBinding { node, model: &model }],None).unwrap();
    let mut input = [vec![2.0,-1.0]];
    let tape = oracle.evaluate(&input).unwrap();
    input[0][0] = 99.0;
    // H(q^2 + |x|^2/2) = 2 grad(q)grad(q)^T + 2q Hq + I.
    // At (2,-1), H = [[37,24],[24,19]]. Constraint adds multiplier*Hq.
    let v = [vec![-0.5,2.0]];
    assert_eq!(tape.objective_hessian_vector_product(&v,None).unwrap(),vec![vec![29.5,26.0]]);
    assert_eq!(tape.lagrangian_hessian_vector_product(&[-3.0],&v,None).unwrap(),vec![vec![32.5,26.0]]);
    assert_eq!(tape.objective_gradient().unwrap(),vec![vec![10.0,5.0]]);
    assert_eq!(model.primals.get(),1);
    assert_eq!(model.products.get(),2,"one action per shared physical node, not per edge");
    assert_eq!(*model.points.borrow(),vec![vec![2.0,-1.0];2]);
}

#[test]
fn separate_physics_variables_preserve_cross_hessian_blocks() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x",Manifold::Rn { dim: 2 },Dims::NONE).unwrap();
    let y = b.var("y",Manifold::Rn { dim: 2 },Dims::NONE).unwrap();
    let px = b.pde_residual("q",x,true,Dims::NONE).unwrap();
    let py = b.pde_residual("q",y,true,Dims::NONE).unwrap();
    let root = b.mul(px,py).unwrap();
    let problem = b.finish(); let model = Model::default();
    let program = ReverseProgram::new_with_physics(&problem,&[root],LIMITS,
        &[PhysicsBinding { node: px, model: &model },PhysicsBinding { node: py, model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0],vec![-1.0,4.0]],None).unwrap();
    // qx=1, qy=13, gx=(4,3), gy=(-2,3), dx=(1,2), dy=(3,-1).
    // Hx=13*(2,0) + gx*(-9); Hy=1*(6,0) + gy*10.
    assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0],vec![3.0,-1.0]],None).unwrap(),
        vec![vec![-10.0,-27.0],vec![-14.0,30.0]]);
    assert_eq!(model.primals.get(),2); assert_eq!(model.products.get(),2);
}

#[test]
fn duplicate_and_zero_seeds_accumulate_without_reissuing_physics_work() {
    let (problem,node) = fixture(); let model = Model::default();
    let program = ReverseProgram::new_with_physics(&problem,&[node,node],LIMITS,
        &[PhysicsBinding { node, model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    assert_eq!(tape.hessian_vector_product(&[1.0,2.0],&[vec![3.0,4.0]],None).unwrap(),vec![vec![18.0,0.0]]);
    assert_eq!(model.products.get(),1);
    assert_eq!(tape.hessian_vector_product(&[1.0,-1.0],&[vec![3.0,4.0]],None).unwrap(),vec![vec![0.0,0.0]]);
    assert_eq!(model.products.get(),1);
}

#[test]
fn zero_outer_gradient_keeps_gauss_term_without_a_local_hessian_call() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x",Manifold::Rn { dim: 2 },Dims::NONE).unwrap();
    let node = b.pde_residual("q",x,true,Dims::NONE).unwrap();
    let one = b.konst(1.0,Dims::NONE).unwrap();
    let r = b.sub(node,one).unwrap(); let root = b.powi(r,2).unwrap();
    let problem = b.finish(); let model = Model::default();
    *model.error.borrow_mut() = Some(PhysicsError::Budget { resource: "linear solves",used: 2,limit: 2 });
    let program = ReverseProgram::new_with_physics(&problem,&[root],LIMITS,
        &[PhysicsBinding { node,model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap(),vec![vec![80.0,60.0]]);
    assert_eq!(model.products.get(),0);
}

#[test]
fn provider_faults_keep_node_attribution_and_the_tape_can_retry() {
    let (problem,node) = fixture(); let model = Model::default();
    let program = ReverseProgram::new_with_physics(&problem,&[node],LIMITS,
        &[PhysicsBinding { node,model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    for source in [
        PhysicsError::Budget { resource: "tangent solves",used: 4,limit: 4 },
        PhysicsError::NotConverged { phase: "incremental adjoint",iterations: 3,
            residual_bits: 0.1f64.to_bits(),tolerance_bits: 1e-9f64.to_bits() },
        PhysicsError::Domain("invalid coefficient"),
    ] {
        *model.error.borrow_mut() = Some(source.clone());
        assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap_err(),
            HessianError::Reverse(ReverseError::Physics { node,source }));
    }
    *model.error.borrow_mut() = None;
    assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap(),vec![vec![2.0,0.0]]);
    assert_eq!(model.products.get(),4); assert_eq!(model.primals.get(),1);
}

#[test]
fn malformed_and_nonfinite_products_are_never_published() {
    let (problem,node) = fixture(); let model = Model::default();
    let program = ReverseProgram::new_with_physics(&problem,&[node],LIMITS,
        &[PhysicsBinding { node,model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    for malformed in 1..=4 {
        model.malformed.set(malformed);
        let error = tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap_err();
        if malformed<=2 {
            assert!(matches!(error,HessianError::PhysicsProductLength { node: n,expected: 2,actual }
                if n==node && actual==if malformed==1 {1} else {3}));
        } else { assert!(matches!(error,HessianError::NonFiniteDerivative { node: n,.. } if n==node)); }
    }
    model.malformed.set(0);
    assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap(),vec![vec![2.0,0.0]]);
}

#[test]
fn invalid_direction_is_refused_before_any_physics_hessian_action() {
    let (problem,node) = fixture(); let model = Model::default();
    let program = ReverseProgram::new_with_physics(&problem,&[node],LIMITS,
        &[PhysicsBinding { node,model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    for direction in [vec![],vec![vec![1.0]],vec![vec![f64::NAN,2.0]]] {
        assert!(tape.hessian_vector_product(&[1.0],&direction,None).is_err());
    }
    assert_eq!(model.products.get(),0);
}

#[test]
fn cancellation_inside_provider_and_before_retry_returns_no_partial_product() {
    let (problem,node) = fixture(); let gate = CancelGate::new();
    let model = Model { gate: Some(&gate),..Model::default() };
    let program = ReverseProgram::new_with_physics(&problem,&[node],LIMITS,
        &[PhysicsBinding { node,model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0,-1.0]],None).unwrap();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate,arena,StreamKey { seed: 0,kernel_id: 1,tile: 0,iteration: 0 },
            Budget::INFINITE,ExecMode::Deterministic);
        for _ in 0..2 {
            assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],Some(&cx)).unwrap_err(),
                HessianError::Reverse(ReverseError::Evaluation(OptError::Cancelled)));
        }
    });
    assert_eq!(model.products.get(),1);
    *model.error.borrow_mut() = Some(PhysicsError::Cancelled);
    assert_eq!(tape.hessian_vector_product(&[1.0],&[vec![1.0,2.0]],None).unwrap_err(),
        HessianError::Reverse(ReverseError::Evaluation(OptError::Cancelled)));
}
