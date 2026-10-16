//! The selected loss image must keep solid/gas source data through real pressure output.
use super::*;
use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;

fn fractional_source(remote: GasState) -> FlueNetworkSpec {
    let mut spec=source();
    spec.graph.sections[0].length_m=0.1501;
    spec.graph.sections[1].length_m=0.1502;
    for section in &mut spec.graph.sections {section.max_length_error_m=1e-12;}
    spec.gases[1]=remote;
    spec.viscothermal=(0..2).map(|section|ViscothermalSectionSpec {section,cells:4,
        minimum_frequency_hz:300.0,maximum_frequency_hz:500.0}).collect();
    spec
}
fn fractional_voice(spec: FlueNetworkSpec, walls: MaterialWallInput) -> FlueVoice {
    let card=mint_tonal_interim_card().unwrap();
    FlueVoice::from_fractional_network(Some(&card),90.0,1e-6,geometry(),spec,Some(walls),
        FractionalTubePolicy {maximum_frequency_hz:1000.0,maximum_summed_phase_error_rad:0.01},
        RATE,3,0,&CancelGate::new_clock_free()).unwrap()
}

#[test]
fn physical_wall_material_and_regional_losses_share_fractional_feedback_and_microphone() {
    let remote=GasState::try_new(&GasSpec::dry_air_ussa1976(),320.0,gas().pressure).unwrap();
    let spec=fractional_source(remote);
    let mut signals=Vec::new();
    for young in [1e6,4e6] {
        let walls=input(young);
        let mut direct=fractional_voice(spec.clone(),walls.clone());
        let n=direct.network().unwrap();
        assert_eq!(n.source().graph,spec.graph);
        assert_eq!(n.source().gases,spec.gases);
        assert_eq!(n.wall_input(),Some(&walls));
        assert_eq!(n.material_walls()[0].source(),&walls.assignments[0]);
        close(n.material_walls()[0].length_m(),0.5*(0.1501+0.1502));
        assert_eq!(n.fractional_sections().unwrap().1.len(),32);
        for (original,report) in spec.graph.sections.iter().zip(n.viscothermal_sections()) {
            assert_eq!(original.length_m,report.requested_length_m);
            assert!((report.represented_length_m-original.length_m).abs()<1e-12);
        }
        direct.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
        let mut receiver=BaffledPressure::circular_outlet(0.006,RATE,config().receiver,
            RayleighMedium {density:remote.density,sound_speed:remote.sound_speed}).unwrap();
        let mut expected=Vec::new();
        for _ in 0..1024 {
            direct.step_block(&mut [0.0]).unwrap();
            let n=direct.network().unwrap();
            let f=n.last_frame().unwrap();
            let scale=(f.stored_energy_j+f.storage_change_j.abs()+f.inlet_work_j.abs()
                +f.interior_loss_j+f.terminal_loss_j).max(1e-30);
            assert!(f.balance_residual_j().abs()<1e-10*scale);
            expected.push(receiver.step(&[n.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap());
        }
        assert!(direct.network().unwrap().wall_state(1).unwrap().stored_energy_j>0.0);
        for partition in [1,37,1024] {
            let mut performance=FlueOutletPerformance::new(
                fractional_voice(spec.clone(),walls.clone()),schedule(),config()).unwrap();
            let mut output=vec![0.0;1024];
            for block in output.chunks_mut(partition) {performance.block(block).unwrap();}
            assert_eq!(output,expected);
            let n=performance.voice().network().unwrap();
            assert_eq!(n.wall_input(),Some(&walls));
            assert_eq!(n.wall_state(1).unwrap(),direct.network().unwrap().wall_state(1).unwrap());
        }
        assert!(expected.iter().any(|p|p.abs()>1e-12));
        signals.push(expected);
    }
    assert_ne!(signals[0],signals[1]);
}

#[test]
fn remote_atmosphere_changes_sound_without_changing_the_wall_or_bore_geometry() {
    let walls=input(1e6);
    let mut signals=Vec::new();
    let mut geometries=Vec::new();
    for t in [293.15,320.0] {
        let remote=GasState::try_new(&GasSpec::dry_air_ussa1976(),t,gas().pressure).unwrap();
        let spec=fractional_source(remote);
        let v=fractional_voice(spec,walls.clone());
        let n=v.network().unwrap();
        assert_eq!(n.wall_input(),Some(&walls));
        assert_eq!(n.section_gases()[16],remote);
        assert_eq!(n.viscothermal_sections()[1].loss.dynamic_viscosity_pa_s(),remote.dynamic_viscosity);
        geometries.push(n.graph().sections.clone());
        let mut performance=FlueOutletPerformance::new(v,schedule(),config()).unwrap();
        let mut out=vec![0.0;1024];
        performance.block(&mut out).unwrap();signals.push(out);
    }
    assert_eq!(geometries[0],geometries[1]);
    assert_ne!(signals[0],signals[1]);
}
