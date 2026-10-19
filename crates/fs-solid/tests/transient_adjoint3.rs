//! G0/G1/G3/G4 pullbacks through the native P1 mass solve and Newmark owner.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3, SurfaceLoad3};
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};
use fs_solid::loads3::transient::adjoint::*;

const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const ZERO: [[f64;3];4] = [[0.;3];4];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 91, kernel_id: 91, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES,&TETS,TetAssemblyBudget::default(),4,cx).unwrap()
}
fn problem<'a>(mat: &'a TetElasticMaterial, fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(mat),
        fixed_dofs:fixed,budget:TetAssemblyBudget::default() }
}
fn config() -> TransientConfig3 {
    TransientConfig3 { maximum_displacement_m:1., mass_damping_per_s:0.4, stiffness_damping_s:0.006,
        solve:TetStaticSolveConfig { relative_residual_tolerance:1e-12,..Default::default() },..Default::default() }
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x,y)|x*y).sum() }
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-7*(1.+a.abs()+b.abs()), "{a:e} != {b:e}"); }

#[test]
fn g1_initial_motion_and_both_force_endpoints_match_native_finite_differences() {
    with_cx(false, |cx| {
        let mat = TetElasticMaterial::try_new(6.,60.,0.25,[91;32]).unwrap();
        let fixed = [0,1,2,4,5,8]; let problem = problem(&mat,&fixed); let geometry = geometry(cx);
        let free = [3,6,7,9,10,11];
        for case in 0..12 {
            let h = 0.005 + case as f64*0.003;
            let q0: Vec<_> = (0..6).map(|i| 0.0002*(i as f64-2.)).collect();
            let v0: Vec<_> = (0..6).map(|i| 0.0003*(3.-i as f64)).collect();
            let f0: Vec<_> = (0..6).map(|i| 0.0001*(i as f64+1.)).collect();
            let f1: Vec<_> = (0..6).map(|i| -0.0002*(i as f64+1.)).collect();
            let directions: Vec<Vec<f64>> = (0..4).map(|j| (0..6)
                .map(|i| 0.01*((i*7+j*11+case*3) as f64).sin()).collect()).collect();
            let seed = MotionCotangent3 { q:vec![0.3,-0.2,0.4,0.1,-0.3,0.5],
                v:vec![0.1,0.2,-0.3,0.4,0.2,-0.1], a:vec![0.05,-0.02,0.01,0.03,-0.04,0.02] };
            let run_at = |epsilon: f64| {
                let mut q = ZERO; let mut v = ZERO;
                for i in 0..6 { q[free[i]/3][free[i]%3]=q0[i]+epsilon*directions[0][i];
                    v[free[i]/3][free[i]%3]=v0[i]+epsilon*directions[1][i]; }
                let first: Vec<_> = free.iter().enumerate().map(|(i,&d)|(d,f0[i]+epsilon*directions[2][i])).collect();
                let last: Vec<_> = free.iter().enumerate().map(|(i,&d)|(d,f1[i]+epsilon*directions[3][i])).collect();
                let initial = geometry.assemble(&[],VolumeLoad3::None,&first,cx).unwrap();
                let next = geometry.assemble(&[],VolumeLoad3::None,&last,cx).unwrap();
                let mut run = initial.begin_transient(&problem,None,0.,&q,&v,config(),cx).unwrap();
                run.advance(&next,h,cx).unwrap(); run
            };
            let baseline = run_at(0.); let adjoint = baseline.adjoint(Default::default(),cx).unwrap();
            let next = adjoint.pullback_step(h,&seed,cx).unwrap();
            let initial = adjoint.pullback_acceleration(&next.previous,cx).unwrap();
            let predicted = dot(&initial.previous.q,&directions[0])+dot(&initial.previous.v,&directions[1])
                +dot(&initial.force,&directions[2])+dot(&next.force,&directions[3]);
            let evaluate = |epsilon| { let r=run_at(epsilon); dot(&seed.q,&r.state().displacement)
                +dot(&seed.v,&r.state().velocity)+dot(&seed.a,&r.state().acceleration) };
            close(predicted,(evaluate(1e-5)-evaluate(-1e-5))/2e-5);
            assert!(initial.previous.a.iter().all(|x| *x==0.));
            assert!(next.true_relative_residual <= 1e-11);
            assert!(initial.true_relative_residual <= 1e-11);
        }
    });
}

#[test]
fn g1_initial_force_gradient_is_not_silently_discarded() {
    with_cx(false, |cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[91;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d| *d!=3).collect(); let p=problem(&mat,&fixed); let g=geometry(cx);
        let load=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run=load.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let a=run.adjoint(Default::default(),cx).unwrap();
        let h=0.1; let b=h*h/4.; let d=h/2.;
        // One free coordinate: M=.1, K=10, C=.1. Initial a0=f0/M.
        let effective=0.1+d*0.1+b*10.;
        let pb=a.pullback_step(h,&MotionCotangent3 { q:vec![1.],v:vec![0.],a:vec![0.] },cx).unwrap();
        let init=a.pullback_acceleration(&pb.previous,cx).unwrap();
        assert!((pb.force[0]-b/effective).abs()<1e-12);
        let expected=(b-(b*10.+d*0.1)*b/effective)/0.1;
        assert!((init.force[0]-expected).abs()<1e-12 && init.force[0].abs()>0.01);
    });
}

#[test]
fn g3_force_cotangents_chain_through_consistent_pressure_integration() {
    with_cx(false, |cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.25,[91;32]).unwrap();
        let fixed=[0,1,2,4,5,8]; let p=problem(&mat,&fixed); let g=geometry(cx);
        let faces:Vec<_>=(0..g.faces().len()).collect();
        let direction=g.assemble(&[SurfaceLoad3::Pressure { faces:&faces,pascals:1. }],VolumeLoad3::None,&[],cx).unwrap();
        let empty=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let a=run.adjoint(Default::default(),cx).unwrap();
        let seed=MotionCotangent3 { q:vec![1.;6],v:vec![0.;6],a:vec![0.;6] };
        let pb=a.pullback_step(0.01,&seed,cx).unwrap();
        let derivative=a.contract_load_direction(&pb.force,&direction,cx).unwrap();
        let mut driven=run.clone(); driven.advance(&direction,0.01,cx).unwrap();
        close(derivative,driven.state().displacement.iter().sum());
        assert!(derivative<0.);
        let constrained=g.assemble(&[],VolumeLoad3::None,&[(0,100.)],cx).unwrap();
        assert_eq!(a.contract_load_direction(&pb.force,&constrained,cx).unwrap(),0.);
    });
}

#[test]
fn g3_cotangent_scaling_and_zero_seed_preserve_the_linear_map() {
    with_cx(false, |cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.25,[91;32]).unwrap();
        let fixed=[0,1,2,4,5,8]; let p=problem(&mat,&fixed); let g=geometry(cx);
        let empty=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let a=run.adjoint(Default::default(),cx).unwrap();
        let seed=MotionCotangent3 { q:vec![0.3;6],v:vec![-0.4;6],a:vec![0.2;6] };
        let base=a.pullback_step(0.013,&seed,cx).unwrap();
        for scale in [1e-120,1e120] {
            let scaled=MotionCotangent3 { q:seed.q.iter().map(|x|x*scale).collect(),
                v:seed.v.iter().map(|x|x*scale).collect(),a:seed.a.iter().map(|x|x*scale).collect() };
            let pb=a.pullback_step(0.013,&scaled,cx).unwrap();
            for (got,expected) in pb.force.iter().zip(&base.force) {
                assert!((got/scale-expected).abs()<1e-9*expected.abs().max(1.));
            }
        }
        let zero=MotionCotangent3 { q:vec![0.;6],v:vec![0.;6],a:vec![0.;6] };
        let pb=a.pullback_step(0.01,&zero,cx).unwrap();
        assert_eq!(pb.previous,zero); assert_eq!(pb.force,vec![0.;6]); assert_eq!(pb.iterations,0);
    });
}

#[test]
fn g1_jump_pullback_overwrites_initial_acceleration_instead_of_adding_impulse() {
    with_cx(false, |cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[91;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d| *d!=3).collect(); let p=problem(&mat,&fixed); let g=geometry(cx);
        let old=g.assemble(&[],VolumeLoad3::None,&[(3,0.03)],cx).unwrap();
        let new=g.assemble(&[],VolumeLoad3::None,&[(3,-0.02)],cx).unwrap();
        let mut q=ZERO; q[1][0]=0.001; let mut v=ZERO; v[1][0]=0.002;
        let mut run=old.begin_transient(&p,None,0.,&q,&v,config(),cx).unwrap();
        let a=run.adjoint(Default::default(),cx).unwrap();
        run.apply_load_jump(&new,cx).unwrap();
        let seed=MotionCotangent3 { q:vec![0.3],v:vec![-0.4],a:vec![0.2] };
        let pb=a.pullback_acceleration(&seed,cx).unwrap();
        assert!((pb.force[0]-2.).abs()<1e-12);
        assert!((pb.previous.q[0]-(0.3-20.)).abs()<1e-12);
        assert!((pb.previous.v[0]-(-0.4-0.2)).abs()<1e-12);
        assert_eq!(pb.previous.a,vec![0.]);
        let init=a.pullback_acceleration(&pb.previous,cx).unwrap();
        assert_eq!(init.force,vec![0.],"overwritten acceleration must not depend on the old force");
    });
}

#[test]
fn g0_g4_bad_dimensions_limits_foreign_geometry_and_cancellation_refuse() {
    with_cx(false, |cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.25,[91;32]).unwrap();
        let p=problem(&mat,&[]); let g=geometry(cx); let empty=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        assert!(matches!(run.adjoint(TransientAdjointConfig3 { maximum_transpose_entries:1,..Default::default() },cx),Err(TransientAdjointError3::Budget)));
        assert!(matches!(run.adjoint(TransientAdjointConfig3 { maximum_free_dofs:1,..Default::default() },cx),Err(TransientAdjointError3::Budget)));
        let a=run.adjoint(Default::default(),cx).unwrap();
        let seed=MotionCotangent3 { q:vec![1.;12],v:vec![0.;12],a:vec![0.;12] };
        for h in [0.,-1.,f64::NAN,f64::INFINITY,1e-300] { assert!(a.pullback_step(h,&seed,cx).is_err()); }
        let mut bad=seed.clone(); bad.v.pop(); assert!(a.pullback_step(0.01,&bad,cx).is_err());
        let mut bad=seed.clone(); bad.a[3]=f64::NAN; assert!(a.pullback_acceleration(&bad,cx).is_err());
        let before=run.state().clone();
        with_cx(true,|cancelled| {
            assert!(matches!(run.adjoint(Default::default(),cancelled),Err(TransientAdjointError3::Cancelled)));
            assert!(matches!(a.pullback_step(0.01,&seed,cancelled),Err(TransientAdjointError3::Cancelled)));
        });
        assert_eq!(run.state(),&before);
        let moved=NODES.map(|p|p.map(|v|v+1.));
        let foreign=TetLoadGeometry3::new(&moved,&TETS,TetAssemblyBudget::default(),4,cx).unwrap();
        let direction=foreign.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(a.contract_load_direction(&vec![1.;12],&direction,cx).is_err());
        assert_eq!(TransientAdjointError3::from(TransientError3::Cancelled),TransientAdjointError3::Cancelled);
    });
}
