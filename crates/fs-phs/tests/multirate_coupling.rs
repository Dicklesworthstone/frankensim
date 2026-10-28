use fs_phs::{PortHamiltonian, QuadraticStorage, Storage, duffing_oscillator, mass_spring_damper};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::multirate::{Clock, IntegrationMethod, MultirateError, MultiratePair, PortPair, SolveConfig};

fn exact(substeps: usize) -> Clock {
    Clock { substeps, method: IntegrationMethod::ExactQuadratic, maximum_contact_step_s: None }
}
fn nonlinear(substeps: usize) -> Clock {
    Clock { substeps, method: IntegrationMethod::DiscreteGradient {
        maximum_rate_rad_s: 5.0, maximum_rate_dt: 0.5 }, maximum_contact_step_s: None }
}
fn pair() -> [PortPair;1] { [PortPair { left:0, right:0, input_scale:1.0 }] }
fn close(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a:e} != {b:e}"); }
fn distance(a:&[f64],b:&[f64])->f64 {
    a.iter().zip(b).map(|(x,y)|(x-y)*(x-y)).sum::<f64>().sqrt()
}

#[test]
fn lossless_energy_remains_closed_for_ten_thousand_windows_on_different_clocks() {
    let left=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let right=mass_spring_damper(1.0,1.0,0.0).unwrap();
    for substeps in [1,2,4,8] {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),0.005,
            exact(1),exact(substeps),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.7,0.2],vec![-0.3,0.4]);
        let initial=left.hamiltonian(&a)+right.hamiltonian(&b);
        for _ in 0..10_000 {
            let r=coupling.step(&a,&b,&[0.0],&[0.0]).unwrap();
            close(r.exchange_work_j.iter().sum(),0.0,2e-13);
            close(r.balance_residual_j,0.0,2e-13);
            assert_eq!(r.external_work_j,[0.0,0.0]);
            assert_eq!(r.left.substeps,1);
            assert_eq!(r.right.substeps,substeps);
            a=r.left.state; b=r.right.state;
        }
        close(left.hamiltonian(&a)+right.hamiltonian(&b),initial,2e-10);
    }
}

#[test]
fn source_work_and_physical_damping_are_separate_from_interface_exchange() {
    let left=mass_spring_damper(1.0,2.0,0.3).unwrap();
    let right=mass_spring_damper(2.0,3.0,0.4).unwrap();
    let mut coupling=MultiratePair::new(&left,&right,&pair(),0.01,
        exact(2),exact(8),SolveConfig::default()).unwrap();
    let (mut a,mut b)=(vec![0.5,0.3],vec![0.1,-0.2]);
    let initial=left.hamiltonian(&a)+right.hamiltonian(&b);
    let (mut work,mut heat)=(0.0,0.0);
    for _ in 0..300 {
        let r=coupling.step(&a,&b,&[0.7],&[-0.2]).unwrap();
        assert!(r.left.dissipated_j>=0.0 && r.right.dissipated_j>=0.0);
        close(r.exchange_work_j.iter().sum(),0.0,1e-13);
        work+=r.external_work_j.iter().sum::<f64>();
        heat+=r.left.dissipated_j+r.right.dissipated_j;
        a=r.left.state; b=r.right.state;
    }
    assert!(heat>0.01);
    close(left.hamiltonian(&a)+right.hamiltonian(&b)-initial+heat-work,0.0,2e-11);
}

fn monolithic() -> PortHamiltonian {
    // Original states [qL,pL,qR,pR], including the full reciprocal reaction.
    PortHamiltonian::new(4,0,
        vec![0.0,1.0,0.0,0.0, -1.0,0.0,0.0,-1.0,
            0.0,0.0,0.0,1.0, 0.0,1.0,-1.0,0.0],vec![0.0;16],vec![],
        Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,
            0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0],4).unwrap())).unwrap()
}

#[test]
fn refining_the_common_window_converges_to_the_full_coupled_exact_flow() {
    let left=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let right=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let whole=monolithic();
    let truth=PreparedExactZoh::new(&whole,2.0,LinearSplitBudget {
        max_modes:4,max_omega_dt:0.04 }).unwrap().step(&[0.7,0.2,-0.3,0.4],&[]).unwrap().x;
    let mut errors=Vec::new();
    for windows in [20,40,80] {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),2.0/windows as f64,
            exact(1),exact(4),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.7,0.2],vec![-0.3,0.4]);
        for _ in 0..windows {
            let r=coupling.step(&a,&b,&[0.0],&[0.0]).unwrap();
            a=r.left.state; b=r.right.state;
        }
        a.extend(b); errors.push(distance(&a,&truth));
    }
    assert!(errors[0]/errors[1]>3.8,"{errors:?}");
    assert!(errors[1]/errors[2]>3.8,"{errors:?}");
}

#[test]
fn nonlinear_native_subcycles_resolve_the_original_duffing_law() {
    let left=duffing_oscillator(1.0,1.0,0.8,0.05).unwrap();
    let right=duffing_oscillator(1.2,1.3,0.6,0.02).unwrap();
    let evolve=|substeps| {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),0.04,
            nonlinear(substeps),nonlinear(substeps),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.6,0.2],vec![-0.4,0.1]);
        for _ in 0..25 {
            let r=coupling.step(&a,&b,&[0.1],&[-0.05]).unwrap();
            close(r.balance_residual_j,0.0,2e-12);
            a=r.left.state; b=r.right.state;
        }
        a.extend(b); a
    };
    let reference=evolve(16);
    let errors:Vec<f64>=[1,2,4,8].into_iter().map(|n|distance(&evolve(n),&reference)).collect();
    assert!(errors.windows(2).all(|e|e[0]>3.0*e[1]),"{errors:?}");
}

fn two_modes() -> PortHamiltonian {
    PortHamiltonian::new(4,3,
        vec![0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0],
        vec![0.0;16],vec![0.0,0.0,0.3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,1.0,0.0],
        Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
            0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.5],4).unwrap())).unwrap()
}

#[test]
fn multiple_ports_keep_independent_sources_and_invariant_numerical_scaling() {
    let (left,right)=(two_modes(),two_modes());
    let normal=[PortPair{left:0,right:1,input_scale:1.0},PortPair{left:1,right:0,input_scale:1.0}];
    let scaled=[PortPair{left:1,right:0,input_scale:100.0},PortPair{left:0,right:1,input_scale:0.01}];
    let solve=|pairs:&[PortPair]| {
        MultiratePair::new(&left,&right,pairs,0.03,exact(2),exact(8),SolveConfig::default())
            .unwrap().step(&[0.3,0.4,-0.2,0.1],&[-0.1,0.2,0.7,-0.3],
                &[0.4,-0.1,0.2],&[-0.2,0.3,-0.1]).unwrap()
    };
    let (a,b)=(solve(&normal),solve(&scaled));
    assert!(distance(&a.left.state,&b.left.state)<1e-13);
    assert!(distance(&a.right.state,&b.right.state)<1e-13);
    close(a.external_work_j.iter().sum(),b.external_work_j.iter().sum(),1e-13);
    assert!(a.external_work_j.iter().any(|x|x.abs()>1e-6));
    close(a.exchange_work_j.iter().sum(),0.0,1e-13);
}

#[test]
fn clock_contact_and_iteration_limits_are_explicit() {
    let (left,right)=(mass_spring_damper(1.0,1.0,0.0).unwrap(),mass_spring_damper(1.0,1.0,0.0).unwrap());
    for n in [0,17] {
        assert!(MultiratePair::new(&left,&right,&pair(),0.01,exact(n),exact(1),SolveConfig::default()).is_err());
    }
    let unresolved=Clock{maximum_contact_step_s:Some(0.001),..exact(1)};
    assert!(MultiratePair::new(&left,&right,&pair(),0.01,unresolved,exact(1),SolveConfig::default()).is_err());
    assert!(MultiratePair::new(&left,&right,&pair(),1.0,nonlinear(1),exact(1),SolveConfig::default()).is_err());
    let mut zero=MultiratePair::new(&left,&right,&pair(),0.01,exact(1),exact(1),SolveConfig {
        max_iterations:0,..SolveConfig::default() }).unwrap();
    assert!(matches!(zero.step(&[0.1,0.3],&[0.2,0.4],&[0.0],&[0.0]),Err(MultirateError::NotConverged{..})));
}

#[test]
fn cancellation_publishes_nothing_and_retry_matches_a_fresh_coupler() {
    let (left,right)=(mass_spring_damper(1.0,1.0,0.0).unwrap(),mass_spring_damper(1.0,1.0,0.0).unwrap());
    let make=||MultiratePair::new(&left,&right,&pair(),0.02,exact(2),exact(8),SolveConfig::default()).unwrap();
    let mut candidate=make();
    let mut calls=0;
    let failed=candidate.step_controlled(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1],||{
        calls+=1; calls>=7
    });
    assert!(matches!(failed,Err(MultirateError::Cancelled)));
    let retry=candidate.step(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1]).unwrap();
    let fresh=make().step(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1]).unwrap();
    assert_eq!(retry,fresh);
}

#[test]
fn prepared_stiff_loss_integral_matches_the_analytic_held_input_solution() {
    let system=PortHamiltonian::new(1,1,vec![0.0],vec![100.0],vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0],1).unwrap())).unwrap();
    let r=PreparedExactZoh::new(&system,0.1,LinearSplitBudget{max_modes:1,max_omega_dt:0.04})
        .unwrap().step(&[1.0],&[1.0]).unwrap();
    let decay=(-10.0_f64).exp();
    let endpoint=0.01+0.99*decay;
    let work=0.001+0.99*(1.0-decay)/100.0;
    let heat=100.0*(0.0001*0.1+2.0*0.01*0.99*(1.0-decay)/100.0
        +0.99*0.99*(1.0-decay*decay)/200.0);
    close(r.x[0],endpoint,2e-14); close(r.supplied,work,2e-14);
    close(r.dissipated,heat,2e-13); close(r.delta_h,0.5*(endpoint*endpoint-1.0),2e-13);
}

struct DatumStorage { datum:f64, defect:f64 }
impl Storage for DatumStorage {
    fn hamiltonian(&self,x:&[f64])->f64 { self.datum+0.5*x[0]*x[0] }
    fn gradient(&self,x:&[f64],out:&mut[f64]) { out[0]=x[0]; }
    fn midpoint_energy_remainder(&self,a:&[f64],b:&[f64],mid:&[f64])->Option<f64> {
        let dx=b[0]-a[0];
        Some(0.5*dx*((a[0]-mid[0])+(b[0]-mid[0]))+self.defect*dx*dx*dx)
    }
}

#[test]
fn energy_admission_cannot_be_diluted_by_an_arbitrary_hamiltonian_datum() {
    for datum in [0.0,1e30] {
        for defect in [0.0,100.0] {
            let law=||PortHamiltonian::new(1,1,vec![0.0],vec![0.0],vec![1.0],
                Box::new(DatumStorage{datum,defect})).unwrap();
            let (left,right)=(law(),law());
            let mut coupling=MultiratePair::new(&left,&right,&pair(),0.04,
                nonlinear(2),nonlinear(4),SolveConfig::default()).unwrap();
            let result=coupling.step(&[0.7],&[0.8],&[0.2],&[0.3]);
            if defect==0.0 { assert!(result.is_ok(),"{result:?}"); }
            else { assert!(matches!(result,Err(MultirateError::Energy{..})),"{result:?}"); }
        }
    }
}
