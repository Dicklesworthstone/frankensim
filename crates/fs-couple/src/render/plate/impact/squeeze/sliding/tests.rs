use super::*;
use crate::render::plate::impact::{ImpactBody, ImpactConfig};
use crate::render::plate::impact::squeeze::{FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm};
use fs_exec::CancelGate;
use fs_tribo::resistive_film::isothermal::AmbientGas;

const M1:f64=0.2;
const M2:f64=0.4;
const GAP:f64=0.001;
const WIDTH:f64=0.02;
const LENGTH:f64=0.03;
const MU:f64=0.2;
const DT:f64=0.01;
const LAUNCH:f64=0.01;

fn bare() -> ImpactSystem {
    let (a,_)=ImpactBody::free_mass(M1,0.0,LAUNCH).unwrap();
    let (b,_)=ImpactBody::free_mass(M2,0.0,0.0).unwrap();
    ImpactSystem::new(vec![a,b],vec![],vec![],vec![],ImpactConfig {
        dt_s:DT,max_steps:128,maximum_energy_j:1.0,
        energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8,
        maximum_generalized_force:1.0,
    }).unwrap()
}
fn graph(channel_gap:f64) -> ResistiveFilm {
    let cell_gap=GapPort {reference_m:GAP,closure:vec![0.0;2]};
    let edge_gap=GapPort {reference_m:channel_gap,closure:vec![0.0;2]};
    ResistiveFilm::new(vec![FilmCell {area_m2:WIDTH*LENGTH,gap:cell_gap}],
        (0..2).map(|_|FilmChannel {from:0,to:None,width_m:WIDTH,length_m:LENGTH,
            gap:edge_gap.clone()}).collect(),2,MU,
        FilmLimits {minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,
            maximum_pressure_pa:1e6}).unwrap()
}
fn film(channel_gap:f64) -> ShearedFilm {
    let rows=[1.0,-1.0].into_iter().map(|sign|WallVelocityPorts {
        lower:vec![sign/M1.sqrt(),0.0],upper:vec![0.0,sign/M2.sqrt()],
    }).collect();
    ShearedFilm::new(graph(channel_gap),rows).unwrap()
}
fn pair() -> ImpactSystem {bare().with_sliding_film(film(GAP)).unwrap()}
fn momentum(x:&[f64])->f64 {M1.sqrt()*x[1]+M2.sqrt()*x[3]}
fn relative(x:&[f64])->f64 {x[1]/M1.sqrt()-x[3]/M2.sqrt()}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} != {b}; tolerance {tol}");}

#[test]
fn moving_wall_drag_advances_both_bodies_on_the_original_clock() {
    // G1: two opposed open passages give p=0 and exact Couette drag.
    // Quadratic free-mass storage makes the native DG update the trapezoidal
    // solution of d(U1-U2)/dt=-S*(1/m1+1/m2)*(U1-U2).
    let drag=2.0*WIDTH*LENGTH*MU/GAP;
    let rate=drag*(1.0/M1+1.0/M2);
    let ratio=(1.0-0.5*DT*rate)/(1.0+0.5*DT*rate);
    let mut system=pair();
    assert_eq!(system.sliding_film_count(),1);
    assert_eq!(system.squeeze_film_count(),0);
    assert_eq!(system.sliding_film(0).unwrap().wall_ports().len(),2);
    let initial=system.stored_energy_j();let p0=momentum(system.state());
    let mut expected=LAUNCH;let mut loss=0.0;
    for tick in 1..=32 {
        let before=relative(system.state());
        let frame=system.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
        let after=relative(system.state());expected*=ratio;
        near(after,expected,2e-11);
        near(frame.dissipated_energy_j,DT*drag*(0.5*(before+after)).powi(2),2e-12);
        near(momentum(system.state()),p0,2e-11);
        assert_eq!(frame.sample,tick);near(frame.time_s,tick as f64*DT,1e-15);
        loss+=frame.dissipated_energy_j;
    }
    assert!(system.state()[3]>0.0);
    near(system.stored_energy_j()+loss,initial,2e-11);
}

#[test]
fn prepared_sliding_force_and_tangent_keep_joint_motion_and_atomic_retry() {
    // G3: the two prepared Jacobian images must advance the same physical law.
    let mut reference=pair();
    for _ in 0..32 {reference.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();}
    for analytic in [false,true] {
        let mut system=if analytic {pair().prepare_analytic().unwrap()} else {pair().prepare().unwrap()};
        let initial=system.stored_energy_j();let mut heat=0.0;
        for tick in 0..32 {
            if tick==8 {
                let state=system.state().to_vec();let samples=system.samples();
                let stopped=CancelGate::new_clock_free();stopped.request();
                assert!(system.step(&[0.0;2],&stopped).is_err());
                assert!(system.step(&[2.0,0.0],&CancelGate::new_clock_free()).is_err());
                assert_eq!(system.state(),state);assert_eq!(system.samples(),samples);
            }
            heat+=system.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap().dissipated_energy_j;
        }
        for (&actual,&expected) in system.state().iter().zip(reference.state()) {
            near(actual,expected,2e-10);
        }
        near(system.stored_energy_j()+heat,initial,2e-11);
    }
}

#[test]
fn sliding_and_squeeze_share_bounded_islands_and_gas_domain_exclusion() {
    let mut system=pair();
    for _ in 0..3 {system=system.with_squeeze_film(graph(GAP)).unwrap();}
    assert!(system.with_sliding_film(film(GAP)).is_err());
    assert!(bare().with_sliding_film(film(0.0)).is_err());
    let gas=AmbientGas {pressure_pa:101325.0,temperature_k:300.0,specific_gas_constant_j_kg_k:287.05};
    assert!(pair().with_compressible_squeeze_film(graph(GAP),gas).is_err());
    let gas_system=bare().with_compressible_squeeze_film(graph(GAP),gas).unwrap();
    assert!(gas_system.with_sliding_film(film(GAP)).is_err());
}
