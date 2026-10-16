use fs_vfit::waveguide::{PassiveWaveguide, WaveguideSpec};
use fs_vfit::waveguide::fractional::FractionalTransit;
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};
use fs_vfit::impedance::SeriesImpedanceSpec;
const DT: f64=1.0/48000.0;

fn network(zero: bool) -> WaveguideNetwork {
    let nodes=[NetworkNode::Inlet,NetworkNode::Junction,
        NetworkNode::Termination {reflection:-0.7},NetworkNode::Impedance {
            load:SeriesImpedanceSpec {resistance_pa_s_m3:2e5,inertance_pa_s2_m3:20.0,compliance_m3_pa:Some(1e-9)}}];
    let segments=[NetworkSegment {nodes:[0,1],one_way_samples:3,impedance_pa_s_m3:1e6},
        NetworkSegment {nodes:[1,2],one_way_samples:5,impedance_pa_s_m3:2e6},
        NetworkSegment {nodes:[1,3],one_way_samples:7,impedance_pa_s_m3:3e6}];
    let d=if zero {[0.0;3]} else {[0.2,0.4,0.8]};
    WaveguideNetwork::with_fractional_transits(&nodes,&segments,&d.map(|d|FractionalTransit::new(d).unwrap()),DT,1<<20).unwrap()
}
#[test]
fn fractional_network_uses_the_original_line_and_actual_storage() {
    let transit=FractionalTransit::new(0.37).unwrap();
    let mut line=PassiveWaveguide::with_fractional_transit(WaveguideSpec {
        one_way_samples:7,impedance_pa_s_m3:1e6,time_step_s:DT,reflection:-0.8,max_memory_bytes:1<<20},transit).unwrap();
    let mut net=WaveguideNetwork::with_fractional_transits(&[NetworkNode::Inlet,NetworkNode::Termination {reflection:-0.8}],
        &[NetworkSegment {nodes:[0,1],one_way_samples:7,impedance_pa_s_m3:1e6}],&[transit],DT,1<<20).unwrap();
    for n in 0..1024 {
        let x=if n<100 {(n as f64*0.3).sin()} else {0.0};
        let a=line.step(x).unwrap();let b=net.step(x).unwrap();
        assert_eq!(a.incoming_pressure_pa,b.incoming_pressure_pa);
        assert_eq!(a.terminal_pressure_pa,net.node_frame(1).unwrap().pressure_pa);
        assert_eq!(a.stored_energy_j,b.stored_energy_j);
        assert_eq!(line.fractional_state_pa(),net.section_fractional_state_pa(0).unwrap());
    }
}
#[test]
fn branched_reactive_network_balances_allpass_load_and_wave_storage() {
    let mut net=network(false);
    for n in 0..4096 {
        let old=net.stored_energy_j();
        let f=net.step(if n<128 { f64::from(n%7)-3.0 } else {0.0}).unwrap();
        let scale=old+f.stored_energy_j+f.inlet_work_j.abs()+f.terminal_loss_j+f.interior_loss_j;
        assert!(f.balance_residual_j().abs()<=3e-12*scale.max(1e-300));
        assert_eq!(net.stored_energy_j(),f.stored_energy_j);
        assert!(f.terminal_loss_j>=0.0);
    }
    assert!(net.section_fractional_transit(1).unwrap().samples()>0.0);
    assert!(net.section_fractional_transit(3).is_none());
}
#[test]
fn paired_network_refusal_cannot_advance_fractional_memory_on_either_side() {
    let (mut a,mut b,mut ra,mut rb)=(network(false),network(false),network(false),network(false));
    for n in 0..20 {
        let x=[f64::from(n%3),-f64::from(n%5)];
        assert_eq!(a.step_pair(&mut b,x).unwrap(),ra.step_pair(&mut rb,x).unwrap());
    }
    assert!(a.step_pair(&mut b,[1.0,f64::MAX]).is_err());
    for i in 0..3 {
        assert_eq!(a.section_fractional_state_pa(i),ra.section_fractional_state_pa(i));
        assert_eq!(b.section_fractional_state_pa(i),rb.section_fractional_state_pa(i));
    }
    for _ in 0..100 {assert_eq!(a.step_pair(&mut b,[0.0;2]).unwrap(),ra.step_pair(&mut rb,[0.0;2]).unwrap());}
}
