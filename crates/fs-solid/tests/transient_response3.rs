//! G1/G3/G4 actual transient stress, full-row support forces and thermal energy.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};

const NODES: [[f64; 3]; 4] = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
const TETS: [[usize; 4]; 1] = [[0,1,2,3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([81;32]);
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed:81,kernel_id:81,tile:0,iteration:0 }, Budget::INFINITE,ExecMode::Deterministic)))
}
fn material(rho: f64, e: f64) -> TetElasticMaterial {
    TetElasticMaterial::try_new(rho,e,0.0,[81;32]).unwrap()
}
fn problem<'a>(m: &'a TetElasticMaterial, fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(m),
        fixed_dofs:fixed,budget:TetAssemblyBudget::default() }
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES,&TETS,TetAssemblyBudget::default(),4,cx).unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs()<1e-10,"{a:e} != {b:e}"); }
fn thermal(m: &TetElasticMaterial, alpha: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(301.0,300.0,[alpha,0.0,0.0,0.0,0.0,0.0],
        m.material_state_identity,ContentHash([82;32]),0.1).unwrap()
}

#[test]
fn g1_reactions_include_consistent_mass_cross_rows_and_both_rayleigh_terms() {
    with_cx(false,|cx| {
        let m=material(6.0,60.0); let fixed:Vec<_>=(0..12).filter(|d| *d!=3).collect();
        let p=problem(&m,&fixed); let g=geometry(cx);
        let load=g.assemble(&[],VolumeLoad3::None,&[(0,0.3),(3,0.01)],cx).unwrap();
        let mut u=ZERO; u[1][0]=0.001; let mut v=ZERO; v[1][0]=0.002;
        let cfg=TransientConfig3 { mass_damping_per_s:0.4,stiffness_damping_s:0.006,..Default::default() };
        let mut run=load.begin_transient(&p,None,0.0,&u,&v,cfg,cx).unwrap();
        for _ in 0..20 {
            let r=run.recover_response(FRAME,0.01,cx).unwrap(); let s=run.state();
            let (q,v,a)=(s.displacement[0],s.velocity[0],s.acceleration[0]);
            // Independent one-tet equations: V=1/6, M_11=.1, M_i1=.05, K_11=10.
            close(r.stress.elements()[0].stress_mandel_pa()[0],60.0*q);
            close(r.stress.elastic_energy_j(),5.0*q*q);
            close(r.kinetic_plus_elastic_energy_j,s.total_energy_j);
            close(r.inertia_force_n[0][0],0.05*a); close(r.inertia_force_n[1][0],0.1*a);
            close(r.damping_force_n[0][0],-0.04*v); close(r.damping_force_n[1][0],0.1*v);
            close(r.support_reactions_n[0][0],0.05*a-0.04*v-10.0*q-0.3);
            for n in [2,3] { close(r.support_reactions_n[n][0],0.05*a+0.02*v); }
            assert_eq!(r.support_reactions_n[1][0],0.0);
            assert!(r.free_force_residual_n<1e-11);
            let balance:f64=(0..4).map(|n| r.external_force_n[n][0]+r.support_reactions_n[n][0]
                -r.inertia_force_n[n][0]-r.damping_force_n[n][0]).sum(); close(balance,0.0);
            let rotated=r.stress.stress_in_frame(0,StressTensorBasis {
                notation:StressTensorNotation::Mandel,order:ElasticTensorOrder::XxYyZzXyYzZx,
                frame:ContentHash([83;32]) },[[0.0,-1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,1.0]]).unwrap();
            close(rotated.stress_pa()[1],60.0*q);
            run.advance(&load,0.01,cx).unwrap();
        }
    });
}

#[test]
fn g1_free_acceleration_has_inertia_but_no_artificial_stress_or_supports() {
    with_cx(false,|cx| {
        let m=material(6.0,60.0); let p=problem(&m,&[]); let g=geometry(cx);
        let a=[0.2,-0.1,0.3];
        let load=g.assemble(&[],VolumeLoad3::Uniform(a.map(|v|6.0*v)),&[],cx).unwrap();
        let mut run=load.begin_transient(&p,None,0.0,&ZERO,&ZERO,TransientConfig3::default(),cx).unwrap();
        for _ in 0..10 { run.advance(&load,0.01,cx).unwrap(); }
        let r=run.recover_response(FRAME,0.01,cx).unwrap();
        assert!(r.stress.elements()[0].stress_mandel_pa().iter().all(|v|v.abs()<1e-10));
        assert_eq!(r.support_reactions_n,ZERO.to_vec());
        for n in 0..4 { for c in 0..3 {
            close(r.inertia_force_n[n][c],0.25*a[c]);
            close(r.inertia_force_n[n][c],r.external_force_n[n][c]);
        } }
        assert!(r.free_force_residual_n<1e-10);
    });
}

#[test]
fn g3_frozen_thermal_stress_and_energy_survive_source_mutation_and_fork() {
    with_cx(false,|cx| {
        let m=material(6.0,60.0); let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed); let g=geometry(cx);
        let load=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let mut source=[thermal(&m,0.001)];
        let mut run=load.begin_transient(&p,Some(TetThermalStrainField::PerElement(&source)),
            0.0,&ZERO,&ZERO,TransientConfig3::default(),cx).unwrap();
        source[0]=thermal(&m,0.02);
        let r=run.recover_response(FRAME,0.01,cx).unwrap();
        close(r.stress.elements()[0].stress_mandel_pa()[0],-0.06);
        close(r.kinetic_plus_elastic_energy_j,5e-6);
        assert_eq!(run.state().total_energy_j,0.0);
        assert!(r.kinetic_plus_elastic_energy_j>run.state().total_energy_j);
        let mut fork=run.clone();
        for _ in 0..8 {
            run.advance(&load,0.01,cx).unwrap(); fork.advance(&load,0.01,cx).unwrap();
            let r=run.recover_response(FRAME,0.01,cx).unwrap();
            assert_eq!(r,fork.recover_response(FRAME,0.01,cx).unwrap());
            let q=run.state().displacement[0];
            close(r.stress.elements()[0].stress_mandel_pa()[0],60.0*(q-0.001));
            close(r.stress.elastic_energy_j(),5.0*(q-0.001).powi(2));
        }
    });
}

#[test]
fn g3_heterogeneous_fields_retain_each_elements_density_stiffness_and_eigenstrain() {
    with_cx(false,|cx| {
        let mut nodes=NODES.to_vec(); nodes.extend(NODES.map(|mut p|{p[0]+=2.0;p}));
        let tets=[[0,1,2,3],[4,5,6,7]];
        let m=[material(6.0,60.0),material(12.0,180.0)];
        let strains=[thermal(&m[0],0.001),thermal(&m[1],0.002)];
        let fixed:Vec<_>=(0..24).filter(|d|![3,15].contains(d)).collect();
        let p=TetLinearElasticProblem { nodes_m:&nodes,tetrahedra:&tets,
            materials:TetMaterialField::PerElement(&m),fixed_dofs:&fixed,budget:TetAssemblyBudget::default() };
        let g=TetLoadGeometry3::new(&nodes,&tets,p.budget,8,cx).unwrap();
        let load=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let mut u=vec![[0.0;3];8]; u[1][0]=0.003;u[5][0]=0.0015;
        let run=load.begin_transient(&p,Some(TetThermalStrainField::PerElement(&strains)),0.0,
            &u,&vec![[0.0;3];8],TransientConfig3::default(),cx).unwrap();
        let r=run.recover_response(FRAME,0.01,cx).unwrap();
        for e in 0..2 {
            let delta=u[4*e+1][0]-strains[e].free_strain_mandel()[0]; let young=[60.0,180.0][e];
            let expected_a=-(young/6.0)*delta/([6.0,12.0][e]/60.0);
            close(r.stress.elements()[e].stress_mandel_pa()[0],young*delta);
            close(r.inertia_force_n[4*e][0],([6.0,12.0][e]/120.0)*expected_a);
        }
        assert!(r.free_force_residual_n<1e-11);
    });
}

#[test]
fn g3_fixed_node_loads_change_reactions_not_motion_and_stay_bound_to_each_fork() {
    with_cx(false,|cx| {
        let m=material(6.0,60.0);let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed);let g=geometry(cx);
        let a=g.assemble(&[],VolumeLoad3::None,&[(3,0.01)],cx).unwrap();
        let b=g.assemble(&[],VolumeLoad3::None,&[(0,0.7),(3,0.01)],cx).unwrap();
        let mut left=a.begin_transient(&p,None,0.0,&ZERO,&ZERO,TransientConfig3::default(),cx).unwrap();
        let mut right=left.clone(); left.advance(&a,0.01,cx).unwrap();right.advance(&b,0.01,cx).unwrap();
        assert_eq!(left.state(),right.state());
        let l=left.recover_response(FRAME,0.01,cx).unwrap();let r=right.recover_response(FRAME,0.01,cx).unwrap();
        close(r.support_reactions_n[0][0]-l.support_reactions_n[0][0],-0.7);
        assert_eq!(l.external_force_n[0][0],0.0);assert_eq!(r.external_force_n[0][0],0.7);
        assert_eq!(l.stress,r.stress);
    });
}

#[test]
fn g4_observation_refusal_and_failed_steps_leave_motion_and_retained_load_unchanged() {
    with_cx(false,|cx| {
        let m=material(6.0,60.0);let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed);let g=geometry(cx);
        let load=g.assemble(&[],VolumeLoad3::None,&[(0,0.2),(3,0.01)],cx).unwrap();
        let mut u=ZERO;u[1][0]=0.001;
        let mut run=load.begin_transient(&p,None,0.0,&u,&ZERO,TransientConfig3::default(),cx).unwrap();
        let before=run.recover_response(FRAME,0.01,cx).unwrap();let state=run.state().clone();
        for ceiling in [0.0,f64::NAN,1e-5] {assert!(run.recover_response(FRAME,ceiling,cx).is_err());}
        assert!(run.recover_response(ContentHash([0;32]),0.01,cx).is_err());
        with_cx(true,|cancelled|assert_eq!(run.recover_response(FRAME,0.01,cancelled),Err(TransientError3::Cancelled)));
        let huge=g.assemble(&[],VolumeLoad3::None,&[(0,1.0),(3,1e6)],cx).unwrap();
        assert!(run.advance(&huge,1.0,cx).is_err());
        assert_eq!(run.state(),&state);assert_eq!(run.recover_response(FRAME,0.01,cx).unwrap(),before);
    });
}
