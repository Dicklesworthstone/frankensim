//! Material, density and damping gradients through the native P1 dynamics.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::transient::{TransientConfig3, TetTransientRun3};
use fs_solid::loads3::transient::adjoint::{MotionCotangent3, TransientAdjointConfig3, TransientAdjointError3};
use fs_solid::loads3::transient::adjoint::material::*;

const NODES: [[f64; 3]; 4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize; 4]; 1] = [[0,1,2,3]];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey { seed: 91,kernel_id: 91,tile: 0,iteration: 0 },Budget::INFINITE,ExecMode::Deterministic)))
}
fn config(eta: f64, zeta: f64) -> TransientConfig3 {
    TransientConfig3 { mass_damping_per_s: eta, stiffness_damping_s: zeta,
        maximum_displacement_m: 1.0, solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-12, maximum_iterations: 10_000 }, ..Default::default() }
}
fn thermal(material: &TetElasticMaterial, eps: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(301.0,300.0,[eps,0.0,0.0,0.0,0.0,0.0],
        material.material_state_identity,fs_blake3::ContentHash([92;32]),0.01).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*(1.0+a.abs().max(b.abs())), "{a} != {b}");
}
fn seed(n: usize) -> MotionCotangent3 {
    MotionCotangent3 { q: (0..n).map(|i| 0.3+i as f64*0.07).collect(),
        v: (0..n).map(|i| -0.2+i as f64*0.03).collect(),
        a: (0..n).map(|i| 0.04-i as f64*0.002).collect() }
}
fn objective(run: &TetTransientRun3<'_, '_, '_, '_>, s: &MotionCotangent3) -> f64 {
    s.q.iter().zip(&run.state().displacement).chain(s.v.iter().zip(&run.state().velocity))
        .chain(s.a.iter().zip(&run.state().acceleration)).map(|(a,b)| a*b).sum()
}
fn one_step(p: [f64;4], cx: &Cx<'_>) -> (f64, [f64;4]) {
    let m = TetElasticMaterial::try_new(6.0*p[1],60.0*p[0],0.0,[91;32]).unwrap();
    let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
    let problem = TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,
        materials:TetMaterialField::Uniform(&m),fixed_dofs:&fixed,budget:Default::default() };
    let geometry = TetLoadGeometry3::new(&NODES,&TETS,problem.budget,4,cx).unwrap();
    let f0 = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.017)],cx).unwrap();
    let f1 = geometry.assemble(&[],VolumeLoad3::None,&[(3,-0.023)],cx).unwrap();
    let thermal = thermal(&m,0.0007);
    let mut q = [[0.0;3];4]; q[1][0]=0.001;
    let mut v = [[0.0;3];4]; v[1][0]=0.002;
    let mut run = f0.begin_transient(&problem,Some(TetThermalStrainField::Uniform(&thermal)),
        0.0,&q,&v,config(p[2],p[3]),cx).unwrap();
    let initial = run.clone();
    let adjoint = run.adjoint(TransientAdjointConfig3::default(),cx).unwrap();
    run.advance(&f1,0.037,cx).unwrap();
    let s = seed(1);
    let step = adjoint.pullback_step(0.037,&s,cx).unwrap();
    let a = run.material_scale_pullback(&step.force,Default::default(),cx).unwrap();
    let start = adjoint.pullback_acceleration(&step.previous,cx).unwrap();
    let b = initial.material_scale_pullback(&start.force,Default::default(),cx).unwrap();
    (objective(&run,&s),[(a.element_scales[0][0]+b.element_scales[0][0])/p[0],
        (a.element_scales[0][1]+b.element_scales[0][1])/p[1],
        a.rayleigh[0]+b.rayleigh[0],a.rayleigh[1]+b.rayleigh[1]])
}

#[test]
fn g1_initial_acceleration_includes_thermal_and_both_damping_derivatives() {
    with_cx(false, |cx| {
        let m = TetElasticMaterial::try_new(6.0,60.0,0.0,[91;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,
            materials:TetMaterialField::Uniform(&m),fixed_dofs:&fixed,budget:Default::default() };
        let geometry = TetLoadGeometry3::new(&NODES,&TETS,problem.budget,4,cx).unwrap();
        let load = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.017)],cx).unwrap();
        let thermal = thermal(&m,0.0007);
        let mut q = [[0.0;3];4]; q[1][0]=0.001;
        let mut v = [[0.0;3];4]; v[1][0]=0.002;
        let run = load.begin_transient(&problem,Some(TetThermalStrainField::Uniform(&thermal)),
            0.0,&q,&v,config(0.3,0.004),cx).unwrap();
        let s = seed(1);
        let adjoint = run.adjoint(Default::default(),cx).unwrap();
        let pb = adjoint.pullback_acceleration(&s,cx).unwrap();
        let g = run.material_scale_pullback(&pb.force,Default::default(),cx).unwrap();
        let lambda = s.a[0]/0.1;
        let a = (0.017+10.0*0.0007-10.0*0.001-(0.3*0.1+0.004*10.0)*0.002)/0.1;
        close(run.state().acceleration[0],a,1e-13);
        close(g.element_scales[0][0],-lambda*10.0*(0.001+0.004*0.002-0.0007),1e-13);
        close(g.element_scales[0][1],-lambda*0.1*(a+0.3*0.002),1e-13);
        close(g.rayleigh[0],-lambda*0.1*0.002,1e-13);
        close(g.rayleigh[1],-lambda*10.0*0.002,1e-13);
        assert!((g.element_scales[0][0]+lambda*10.0*(0.001+0.004*0.002)).abs()>1e-4);
    });
}

#[test]
fn g3_full_native_step_and_initialization_match_material_finite_differences() {
    with_cx(false, |cx| {
        for p in [[1.0,1.0,0.3,0.004],[1.7,0.6,0.8,0.009]] {
            let (_,g) = one_step(p,cx);
            for axis in 0..4 {
                let h = 1e-5*p[axis];
                let mut plus=p; let mut minus=p; plus[axis]+=h; minus[axis]-=h;
                let fd=(one_step(plus,cx).0-one_step(minus,cx).0)/(2.0*h);
                close(g[axis],fd,2e-8);
            }
        }
    });
}

#[test]
fn g3_anisotropic_heterogeneous_residual_contraction_matches_original_assembly() {
    with_cx(false, |cx| {
        let nodes = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,1.,1.]];
        let tets = [[0,1,2,3],[1,2,3,4]]; let fixed = [0,1,2,4,5,8];
        let make = |scales: [[f64;2];2]| -> Vec<TetElasticMaterial> {
            (0..2).map(|e| {
                let c = std::array::from_fn(|i| std::array::from_fn(|j|
                    scales[e][0]*((if i==j {40.0+10.0*i as f64} else {0.0})
                        +(i+1) as f64*(j+1) as f64*(e+1) as f64)));
                TetElasticMaterial::try_new_mandel((6.0+3.0*e as f64)*scales[e][1],c,[91+e as u8;32]).unwrap()
            }).collect()
        };
        let materials=make([[1.0;2];2]);
        let thermal: Vec<_> = materials.iter().enumerate().map(|(i,m)| thermal(m,0.0002*(i+1) as f64)).collect();
        let problem = TetLinearElasticProblem { nodes_m:&nodes,tetrahedra:&tets,
            materials:TetMaterialField::PerElement(&materials),fixed_dofs:&fixed,budget:Default::default() };
        let geometry = TetLoadGeometry3::new(&nodes,&tets,problem.budget,8,cx).unwrap();
        let load = geometry.assemble(&[],VolumeLoad3::Uniform([0.02,-0.03,0.01]),&[],cx).unwrap();
        let mut q=[[0.0;3];5]; let mut v=q;
        for d in 0..15 { if !fixed.contains(&d) { q[d/3][d%3]=1e-4*(d+1) as f64; v[d/3][d%3]=-2e-4*(d+1) as f64; } }
        let mut run=load.begin_transient(&problem,Some(TetThermalStrainField::PerElement(&thermal)),
            0.0,&q,&v,config(0.3,0.004),cx).unwrap();
        run.advance(&load,0.02,cx).unwrap();
        let lambda: Vec<_> = (0..run.free_dofs().len()).map(|i| 0.2-0.03*i as f64).collect();
        let gradient=run.material_scale_pullback(&lambda,Default::default(),cx).unwrap();
        let residual_objective=|scales,eta:f64,zeta:f64| {
            let ms=make(scales);
            let problem=TetLinearElasticProblem { materials:TetMaterialField::PerElement(&ms),..problem };
            let operators=problem.assemble(cx).unwrap();
            let ft=problem.assemble_thermal_load(TetThermalStrainField::PerElement(&thermal),cx).unwrap();
            let n=lambda.len(); let mut ma=vec![0.0;n]; let mut mv=ma.clone(); let mut ku=ma.clone(); let mut kv=ma.clone();
            operators.mass.spmv(&run.state().acceleration,&mut ma);
            operators.mass.spmv(&run.state().velocity,&mut mv);
            operators.stiffness.spmv(&run.state().displacement,&mut ku);
            operators.stiffness.spmv(&run.state().velocity,&mut kv);
            (0..n).map(|i| -lambda[i]*(ma[i]+eta*mv[i]+ku[i]+zeta*kv[i]-ft.reduced_force_n[i])).sum::<f64>()
        };
        for e in 0..2 { for axis in 0..2 {
            let mut plus=[[1.0;2];2]; let mut minus=plus; plus[e][axis]+=1e-5; minus[e][axis]-=1e-5;
            let fd=(residual_objective(plus,0.3,0.004)-residual_objective(minus,0.3,0.004))/2e-5;
            close(gradient.element_scales[e][axis],fd,1e-9);
        } }
        close(gradient.rayleigh[0],(residual_objective([[1.0;2];2],0.30001,0.004)
            -residual_objective([[1.0;2];2],0.29999,0.004))/2e-5,1e-9);
        close(gradient.rayleigh[1],(residual_objective([[1.0;2];2],0.3,0.00401)
            -residual_objective([[1.0;2];2],0.3,0.00399))/2e-5,1e-9);
    });
}

#[test]
fn g3_density_dependent_acceleration_load_needs_its_explicit_force_chain() {
    with_cx(false, |cx| {
        let m=TetElasticMaterial::try_new(6.0,60.0,0.25,[91;32]).unwrap();
        let problem=TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,
            materials:TetMaterialField::Uniform(&m),fixed_dofs:&[],budget:Default::default() };
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,problem.budget,4,cx).unwrap();
        let load=geometry.assemble(&[],VolumeLoad3::Uniform([0.0,0.0,3.0]),&[],cx).unwrap();
        let run=load.begin_transient(&problem,None,0.0,&[[0.0;3];4],&[[0.0;3];4],config(0.0,0.0),cx).unwrap();
        let mut s=MotionCotangent3 {q:vec![0.0;12],v:vec![0.0;12],a:vec![0.0;12]}; s.a[2]=1.0;
        let adjoint=run.adjoint(Default::default(),cx).unwrap();
        let pb=adjoint.pullback_acceleration(&s,cx).unwrap();
        let g=run.material_scale_pullback(&pb.force,Default::default(),cx).unwrap();
        close(g.element_scales[0][1],-0.5,1e-11);
        let force_chain=adjoint.contract_load_direction(&pb.force,&load,cx).unwrap();
        close(g.element_scales[0][1]+force_chain,0.0,1e-11);
    });
}

#[test]
fn g0_g4_material_observation_bounds_and_cancellation_preserve_the_primal() {
    with_cx(false, |cx| {
        let m=TetElasticMaterial::try_new(6.0,60.0,0.0,[91;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let problem=TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,
            materials:TetMaterialField::Uniform(&m),fixed_dofs:&fixed,budget:Default::default() };
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,problem.budget,4,cx).unwrap();
        let load=geometry.assemble(&[],VolumeLoad3::None,&[(3,0.1)],cx).unwrap();
        let run=load.begin_transient(&problem,None,0.0,&[[0.0;3];4],&[[0.0;3];4],config(0.0,0.0),cx).unwrap();
        let before=run.state().clone();
        assert!(run.material_scale_pullback(&[],Default::default(),cx).is_err());
        assert!(run.material_scale_pullback(&[f64::NAN],Default::default(),cx).is_err());
        for budget in [MaterialGradientBudget3 {maximum_elements:0,..Default::default()},
            MaterialGradientBudget3 {maximum_full_dofs:11,..Default::default()}] {
            assert_eq!(run.material_scale_pullback(&[1.0],budget,cx),Err(TransientAdjointError3::Budget));
        }
        with_cx(true,|cancelled| assert_eq!(run.material_scale_pullback(&[1.0],Default::default(),cancelled),
            Err(TransientAdjointError3::Cancelled)));
        let gradient=run.material_scale_pullback(&[1.0],Default::default(),cx).unwrap();
        assert!(gradient.contract(&[],[0.0;2],cx).is_err());
        assert!(gradient.contract(&[[f64::INFINITY,0.0]],[0.0;2],cx).is_err());
        close(gradient.contract(&[[2.0,3.0]],[0.0;2],cx).unwrap(),
            2.0*gradient.element_scales[0][0]+3.0*gradient.element_scales[0][1],1e-14);
        assert_eq!(run.state(),&before);
    });
}
