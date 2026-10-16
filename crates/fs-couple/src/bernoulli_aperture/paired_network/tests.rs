use super::*;
use super::super::paired::{ApertureBlade, BladeState, PairedApertureLimits, PairedApertureSpec};
use super::super::network::TubeSection;
use fs_dcontact::Obstacle;
use fs_material::gas::GasSpec;

fn gas(t:f64,p:f64)->GasState {GasState::try_new(&GasSpec::dry_air_ussa1976(),t,p).unwrap()}
fn make(second_mass:f64,temperature:f64)->PairedApertureNetwork {
    let gases=vec![gas(293.15,101325.0),gas(temperature,101325.0)];
    let graph=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,
        NetworkNode::Termination {reflection:-0.8}],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.12,radius_m:0.006,max_length_error_m:0.005},
            TubeSection {nodes:[1,2],length_m:0.18,radius_m:0.007,max_length_error_m:0.005}],
        sound_speed_m_s:gases[0].sound_speed,max_wave_memory_bytes:1<<20};
    let z=graph.inlet_impedance_with_gases(&gases).unwrap();
    let spec=PairedApertureSpec {rest_opening_m:4e-4,width_m:0.013,
        blades:[ApertureBlade {mass_kg:1e-5,stiffness_n_m:350.0,damping_n_s_m:0.03,pressure_area_m2:1e-4},
            ApertureBlade {mass_kg:second_mass,stiffness_n_m:700.0,damping_n_s_m:0.04,pressure_area_m2:1.2e-4}],
        density_kg_m3:gases[0].density,impedance_pa_s_m3:z,time_step_s:1.0/48000.0,max_steps:4096,
        limits:PairedApertureLimits {maximum_abs_displacement_m:0.05,maximum_abs_velocity_m_s:100.0,
            maximum_energy_j:10.0,max_evaluations:256,velocity_tolerance_m_s:1e-15,relative_tolerance:1e-13}};
    let law=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],2e9,2.0,"synthetic shared contact".into()).unwrap();
    let a=PairedAperture::new(spec,[BladeState::default();2],law).unwrap();
    PairedApertureNetwork::new(a,graph,gases,&CancelGate::new_clock_free()).unwrap()
}

#[test]
fn independent_blades_and_regional_gas_change_actual_outlet_motion() {
    let mut base=make(2e-5,293.15);let mut mass=make(5e-5,293.15);let mut warm=make(2e-5,323.15);
    assert_eq!(base.represented_sections()[0],warm.represented_sections()[0]);
    assert_ne!(base.represented_sections()[1].one_way_samples,warm.represented_sections()[1].one_way_samples);
    let mut mass_delta=0.0;let mut gas_delta=0.0;
    for i in 0..512 {
        let drive=TubeDrive {upstream_pressure_pa:if i<180 {500.0} else {0.0},body_flow_m3_s:0.0};
        let f=base.step(drive).unwrap();mass.step(drive).unwrap();warm.step(drive).unwrap();
        assert!(f.coupled.balance_residual_j().abs()<1e-12);
        assert_ne!(f.mechanics.blades[0],f.mechanics.blades[1]);
        let q=base.node_frame(2).unwrap().net_flow_into_node_m3_s;
        mass_delta+=(q-mass.node_frame(2).unwrap().net_flow_into_node_m3_s).abs();
        gas_delta+=(q-warm.node_frame(2).unwrap().net_flow_into_node_m3_s).abs();
    }
    assert!(mass_delta>1e-8 && gas_delta>1e-8);
    assert!(base.stored_energy_j()>0.0,"pressure release must retain ringdown");
}

#[test]
fn block_partition_and_invalid_inputs_preserve_the_complete_network() {
    let mut a=make(2e-5,293.15);let mut b=make(2e-5,293.15);
    let drive=TubeDrive {upstream_pressure_pa:500.0,body_flow_m3_s:0.0};
    let inputs=vec![drive;127];let mut block=vec![PairedNetworkFrame::default();127];
    let gate=CancelGate::new_clock_free();
    a.advance_block(&inputs,&mut block,&gate).unwrap();
    for expected in &block { assert_eq!(&b.step(drive).unwrap(),expected); }
    let before=*a.aperture().blade_states();let energy=a.stored_energy_j();let node=*a.node_frame(2).unwrap();
    assert!(a.step_with_forces(drive,[f64::NAN,0.0]).is_err());
    assert_eq!(a.aperture().blade_states(),&before);assert_eq!(a.stored_energy_j(),energy);
    assert_eq!(*a.node_frame(2).unwrap(),node);
    assert_eq!(a.step(drive).unwrap(),b.step(drive).unwrap());
    // Rejected active terminal controls do not replace the existing passive load.
    assert!(a.set_terminal_reflection(2,1.1).is_err());
    assert_eq!(a.step(drive).unwrap(),b.step(drive).unwrap());
}

#[test]
fn an_exterior_receiver_uses_the_outlet_gas_and_retained_physical_flow() {
    let mut model=make(2e-5,323.15);
    let position=CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,
        angular_points:16,maximum_frequency_hz:2000.0};
    let mut receiver=model.baffled_receiver(2,48000,position).unwrap();
    assert!(model.baffled_receiver(1,48000,position).is_err());
    assert!(model.baffled_receiver(2,24000,position).is_err());
    let g=model.section_gases()[1];
    let mut reference=BaffledPressure::circular_outlet(0.007,48000,position,
        RayleighMedium {density:g.density,sound_speed:g.sound_speed}).unwrap();
    let mut nonzero=false;
    for i in 0..256 {
        model.step(TubeDrive {upstream_pressure_pa:if i<128 {500.0} else {0.0},body_flow_m3_s:0.0}).unwrap();
        let q=model.node_frame(2).unwrap().net_flow_into_node_m3_s;
        let p=receiver.step(&[q]).unwrap();
        assert_eq!(p.to_bits(),reference.step(&[q]).unwrap().to_bits());
        nonzero|=p!=0.0;
    }
    assert!(nonzero);assert!(model.baffled_receiver(2,48000,position).is_err());
}
