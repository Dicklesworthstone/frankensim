//! cargo run --release -p fs-ascent --features conduction-design --example conductivity_newton
//! Synthetic anisotropic heated cube; no measured-material or continuum claim.
use fs_ascent::conduction_design::ConductivityNewtonModel;
use fs_ascent::{ReverseNewtonStudy,StopRule};
use fs_ascent::trust::KrylovOptions;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::{ConductionMesh,ConductionProblem,ConductivityModel,LinearConfig,ScalarField,
    ThermalBc,ThermalBoundaryBuilder,TEMPERATURE_DIMS};
use fs_conduction::fixtures::{unit_cube,on_box_face};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding,PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{Manifold,ProblemBuilder,ReverseProblem,Sense};
use fs_qty::Dims;

fn main()->Result<(),Box<dyn std::error::Error>> {
    let (complex,positions)=unit_cube(2);let mesh=ConductionMesh::new(complex,positions)?;
    let boundary=ThermalBoundaryBuilder::new(&mesh)
        .region("cooled base",|f| on_box_face(f.centroid[2],0.0),ThermalBc::dirichlet(300.0)?)?
        .region("ambient side",|f| on_box_face(f.centroid[0],1.0),ThermalBc::robin(8.0,290.0)?)?
        .adiabatic_remainder().finish()?;
    let material=ConductivityModel::constant_tensor([[3.0,0.5,0.25],[0.5,2.0,0.3],[0.25,0.3,1.5]])?;
    let source=ScalarField::Uniform(4000.0);
    let design=ConductivityDesign::new(ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
        source:&source,element_materials:None},LinearConfig {tolerance:1e-12,max_iterations:2000,restart:30})?;
    let n=design.parameter_count();let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
    let model=ConductivityNewtonModel::new(&design,"heated cube",weights,4000)?;
    let mut b=ProblemBuilder::new();
    let x=b.var("log conductivity multipliers",Manifold::Rn {dim:n as u32},Dims::NONE)?;
    let xr=b.var_ref(x)?;let penalty=b.norm_sq(xr)?;
    let penalty_weight=b.konst(0.05/n as f64,Dims::NONE)?;let penalty=b.mul(penalty_weight,penalty)?;
    let temperature=b.pde_residual("heated cube",x,true,TEMPERATURE_DIMS)?;
    let reference=b.konst(300.0,TEMPERATURE_DIMS)?;let scale=b.konst(1000.0,TEMPERATURE_DIMS)?;
    let rise=b.sub(temperature,reference)?;let rise=b.div(rise,scale)?;let fit=b.powi(rise,2)?;
    let objective=b.add(fit,penalty)?;b.objective(objective,Sense::Minimize,1.0)?;let problem=b.finish();
    let oracle=ReverseProblem::new_with_physics(&problem,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
        &[PhysicsBinding {node:temperature,model:&model}],None)?;
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| -> Result<(),Box<dyn std::error::Error>> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:1,kernel_id:81,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let mut study=ReverseNewtonStudy::new(&oracle,&vec![0.0;n],Some(&cx))?;
        study.set_krylov_options(KrylovOptions::new(0.05,4)?);
        let initial=study.snapshot().f;
        // Distinct ceilings: objective samples, IR Hessian actions, and the
        // provider's conservative native linear-solve reservations.
        let report=study.run(&StopRule::Any(vec![StopRule::GradNorm(1e-8),StopRule::Budget(80)]),20,200,Some(&cx))?;
        let sample=model.value_gradient(&report.solution.x,Some(&cx))?;
        let rho:Vec<_>=report.solution.x.iter().map(|&x|fs_math::det::exp(x)).collect();
        println!("objective: {initial:.9e} -> {:.9e}; stop: {:?}",report.solution.f,report.stop);
        println!("weighted free-node temperature: {:.6} K; gradient norm: {:.6e}",sample.value,report.solution.grad_norm);
        println!("IR samples: {}; Hessian actions: {}; native solve reservations: {} / {}",
            report.solution.evals,report.solution.hv_evals,model.reserved_solves(),model.solve_limit());
        println!("positive element conductivity multipliers: {rho:?}");
        Ok(())
    })
}
