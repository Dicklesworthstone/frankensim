//! Pressure/flow duality against independent physical RC equations.
use crate::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};
use crate::waveguide::network::admittance::{AdmittanceTerm,RelaxationAdmittance,RelaxationAdmittanceSpec};
use crate::impedance::{SeriesImpedance,SeriesImpedanceSpec};
fn near(a:f64,b:f64,scale:f64) {assert!((a-b).abs()<=2e-11*scale.max(1e-30),"{a:e} != {b:e}");}

#[test]
fn rc_memory_matches_independent_midpoint_pressure_flow_and_work() {
    let(g,c,z,dt)=(1e-8,1e-10,1e6,1.0/96000.0);
    let arms=[AdmittanceTerm{conductance_m3_pa_s:2e-7,rate_per_s:500.0},
        AdmittanceTerm{conductance_m3_pa_s:5e-8,rate_per_s:2000.0}];
    let mut load=RelaxationAdmittance::new(RelaxationAdmittanceSpec::new(g,c,&arms).unwrap(),z,dt).unwrap();
    let(mut old,mut history)=(0.0,[0.0;2]);
    for n in 0..2000 {
        let a=if n<1500 {((n%43) as f64-21.0)*0.3}else{0.0};
        let mut denominator=1.0+z*(g+2.0*c/dt);
        let mut rhs=2.0*a+z*2.0*c/dt*old;
        for (j,t) in arms.iter().enumerate() {
            let k=t.conductance_m3_pa_s/(1.0+t.rate_per_s*dt/2.0);
            denominator+=z*k;rhs+=z*k*history[j];
        }
        let p=rhs/denominator;let new_pressure=2.0*p-old;
        let before=load.stored_energy_j();
        let mut energy=0.5*c*new_pressure*new_pressure;
        let mut loss=g*p*p*dt;
        for(j,t)in arms.iter().enumerate(){
            let next=history[j]+dt*t.rate_per_s/(1.0+t.rate_per_s*dt/2.0)*(p-history[j]);
            let diff=p-0.5*(history[j]+next);
            loss+=t.conductance_m3_pa_s*diff*diff*dt;
            energy+=0.5*t.conductance_m3_pa_s/t.rate_per_s*next*next;
            history[j]=next;
        }
        let f=load.step(a).unwrap();
        near(f.pressure_pa,p,p.abs()+a.abs());near(f.flow_m3_s,(2.0*a-p)/z,1e-4);
        near(f.stored_energy_j,energy,energy);near(f.dissipated_energy_j,loss,loss);
        near(f.balance_residual_j(),0.0,before+energy+loss+f.supplied_work_j.abs());
        for (actual,expected) in load.branch_pressures_pa().iter().zip(history) {near(*actual,expected,expected.abs().max(1.0));}
        near(load.base_pressure_pa(),new_pressure,new_pressure.abs().max(1.0));old=new_pressure;
    }
}

#[test]
fn a_single_parallel_arm_is_the_same_physical_series_rc_load() {
    let(g,p,z,dt)=(1e-7,1000.0,2e6,1.0/48000.0);
    let mut a=RelaxationAdmittance::new(RelaxationAdmittanceSpec::new(0.0,0.0,
        &[AdmittanceTerm{conductance_m3_pa_s:g,rate_per_s:p}]).unwrap(),z,dt).unwrap();
    let mut b=SeriesImpedance::new(SeriesImpedanceSpec{resistance_pa_s_m3:1.0/g,
        inertance_pa_s2_m3:0.0,compliance_m3_pa:Some(g/p)},z,dt).unwrap();
    for n in 0..500 {
        let incident=if n<333 {(n%29) as f64-14.0}else{0.0};
        let x=a.step(incident).unwrap();let y=b.step(incident).unwrap();
        near(x.pressure_pa,y.pressure_pa,30.0);near(x.flow_m3_s,y.flow_m3_s,30.0/z);
        near(x.stored_energy_j,y.stored_energy_j,1e-8);
        near(x.dissipated_energy_j,y.dissipated_energy_j,1e-8);
    }
    let pressures=a.branch_pressures_pa().to_vec();let energy=a.stored_energy_j();
    for bad in [f64::NAN,f64::INFINITY,f64::MAX] {assert!(a.step(bad).is_err());}
    assert_eq!(a.branch_pressures_pa(),pressures);assert_eq!(a.stored_energy_j(),energy);
    assert!(RelaxationAdmittanceSpec::new(-1.0,0.0,&[]).is_err());
    assert!(RelaxationAdmittanceSpec::new(0.0,f64::NAN,&[]).is_err());
}

#[test]
fn thermal_shunt_shares_network_pressure_without_lag_and_retains_retry_state() {
    let y=RelaxationAdmittanceSpec::new(0.0,1e-10,
        &[AdmittanceTerm{conductance_m3_pa_s:1e-7,rate_per_s:800.0}]).unwrap();
    let make=||WaveguideNetwork::new(&[NetworkNode::Inlet,NetworkNode::ShuntAdmittance{load:y},
        NetworkNode::Termination{reflection:-0.8}],&[
        NetworkSegment{nodes:[0,1],one_way_samples:5,impedance_pa_s_m3:1e6},
        NetworkSegment{nodes:[1,2],one_way_samples:7,impedance_pa_s_m3:2e6}],1.0/48000.0,1<<20).unwrap();
    let(mut a,mut b)=(make(),make());let mut seen_loss=false;
    for n in 0..700 {
        let input=if n<500 {(n%23) as f64-11.0}else{0.0};let before=a.stored_energy_j();
        if n==300 {
            let p=a.shunt_relaxation_pressures(1).unwrap().to_vec();let node=*a.node_frame(1).unwrap();
            assert!(a.step(f64::MAX).is_err());assert_eq!(a.node_frame(1),Some(&node));
            assert_eq!(a.shunt_relaxation_pressures(1).unwrap(),p);
        }
        let f=a.step(input).unwrap();assert_eq!(f,b.step(input).unwrap());
        near(f.balance_residual_j(),0.0,before+f.stored_energy_j+f.inlet_work_j.abs()+f.interior_loss_j+f.terminal_loss_j);
        let node=a.node_frame(1).unwrap();near(node.net_flow_into_node_m3_s,node.load_flow_m3_s,1e-4);
        seen_loss|=node.absorbed_energy_j>0.0;
    }
    assert!(seen_loss);assert!(a.terminal_state(1).is_none());
    assert!(a.terminal_relaxation_flows(1).is_none());
    assert!(a.set_terminal_reflection(1,0.0).is_err());
}
