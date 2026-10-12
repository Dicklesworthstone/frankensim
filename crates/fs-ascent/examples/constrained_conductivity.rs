//! cargo run --release -p fs-ascent --features conduction-design --example constrained_conductivity
//! Synthetic native P1 design with an explicit average conductivity allowance.
use fs_ascent::conduction_design::ConductivityNewtonModel;
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStudy};
use fs_ascent::trust::KrylovOptions;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundaryBuilder, TEMPERATURE_DIMS};
use fs_conduction::fixtures::{unit_cube, on_box_face};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (complex, positions)=unit_cube(1);let mesh=ConductionMesh::new(complex,positions)?;
    let boundary=ThermalBoundaryBuilder::new(&mesh)
        .region("base",|f|on_box_face(f.centroid[2],0.0),ThermalBc::dirichlet(300.0)?)?
        .adiabatic_remainder().finish()?;
    let material=ConductivityModel::constant_tensor([[3.0,0.5,0.25],[0.5,2.0,0.3],[0.25,0.3,1.5]])?;
    let source=ScalarField::Uniform(4000.0);
    let design=ConductivityDesign::new(ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
        source:&source,element_materials:None},LinearConfig {tolerance:1e-12,max_iterations:2000,restart:30})?;
    let n=design.parameter_count();let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
    let model=ConductivityNewtonModel::new(&design,"heated cube",weights,20000)?;
    let mut b=ProblemBuilder::new();let v=b.var("log element multipliers",Manifold::Rn {dim:n as u32},Dims::NONE)?;
    let x=b.var_ref(v)?;
    let temperature=b.pde_residual("heated cube",v,true,TEMPERATURE_DIMS)?;
    let reference=b.konst(300.0,TEMPERATURE_DIMS)?;let scale=b.konst(1000.0,TEMPERATURE_DIMS)?;
    let rise=b.sub(temperature,reference)?;let rise=b.div(rise,scale)?;
    // Mild regularization discourages extreme log coefficients; the allowance
    // below is still an actual constraint, not substituted by this penalty.
    let norm=b.norm_sq(x)?;let regularization=b.konst(0.02/n as f64,Dims::NONE)?;
    let regularization=b.mul(regularization,norm)?;let objective=b.add(rise,regularization)?;
    b.objective(objective,Sense::Minimize,1.0)?;
    let mut sum=b.konst(0.0,Dims::NONE)?;
    for i in 0..n { let xi=b.component(x,i as u32)?;let rho=b.exp(xi)?;sum=b.add(sum,rho)?; }
    let reciprocal=b.konst(1.0/n as f64,Dims::NONE)?;let average=b.mul(reciprocal,sum)?;
    let allowance=b.konst(1.5,Dims::NONE)?;let c=b.sub(average,allowance)?;
    b.constraint(c,ConstraintKind::LeZero,"mean conductivity multiplier <= 1.5")?;
    let problem=b.finish();
    let oracle=ReverseProblem::new_with_physics(&problem,ReverseLimits {max_nodes:512,max_scalar_slots:4096},
        &[PhysicsBinding {node:temperature,model:&model}],None)?;
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| -> Result<(),Box<dyn std::error::Error>> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:1,kernel_id:82,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let mut study=AugmentedNewtonStudy::new(&oracle,&vec![0.0;n],AugmentedNewtonOptions::default(),Some(&cx))?;
        study.set_krylov_options(KrylovOptions::new(0.01,12)?);
        let initial=study.inner().objective_value();
        let report=study.run(1e-5,150,AugmentedNewtonBudget {
            objective_evals:500,hessian_products:2000,multiplier_updates:50,
        },Some(&cx))?;
        let rho:Vec<_>=report.x.iter().map(|&x|fs_math::det::exp(x)).collect();
        let sample=model.value_gradient(&report.x,Some(&cx))?;
        println!("original objective: {initial:.9e} -> {:.9e}; stop: {:?}",report.objective,report.stop);
        println!("temperature: {:.6} K; mean multiplier: {:.9} (allowance 1.5)",sample.value,rho.iter().sum::<f64>()/n as f64);
        println!("original constraint residuals: {:?}; candidate multipliers: {:?}",report.constraints,report.multipliers);
        println!("KKT numerical residuals: {:?}",report.kkt);
        println!("Newton steps: {}; multiplier updates: {}; IR samples: {}; augmented Hessians: {}; native reservations: {}",
            report.newton_iterations,report.multiplier_updates,report.objective_evals,report.hessian_products,model.reserved_solves());
        println!("element multipliers: {rho:?}");
        // Do not relabel a bounded or stalled solve as a converged design.
        Ok(())
    })
}
