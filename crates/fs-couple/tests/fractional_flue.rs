//! The original jet island and wave network consume the optional phase image.
use fs_couple::flue_loop::{FlueControl, FlueError, FlueGeometry, FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
const DT: f64 = 1.0 / 48_000.0;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap()
}
fn policy() -> FractionalTubePolicy {
    FractionalTubePolicy { maximum_frequency_hz:3000.0, maximum_summed_phase_error_rad:0.01 }
}
fn geometry() -> FlueGeometry {
    FlueGeometry { flue_width_m:0.010,flue_height_m:0.001,cut_up_m:0.004,
        labium_offset_m:0.0,jet_angle_rad:0.0,
        bore_area_m2:core::f64::consts::PI*0.006*0.006,mouth_loss_rel:0.02 }
}
fn source(delay: f64) -> FlueNetworkSpec {
    let g=gas();
    FlueNetworkSpec {
        graph:TubeNetworkSpec { nodes:vec![NetworkNode::Inlet,NetworkNode::Termination {reflection:-0.9}],
            sections:vec![TubeSection {nodes:[0,1],length_m:delay*(g.sound_speed*DT),
                radius_m:0.006,max_length_error_m:1e-12}],
            sound_speed_m_s:g.sound_speed,max_wave_memory_bytes:1_000_000 },
        gases:vec![g],viscothermal:vec![],minimum_jet_speed_m_s:1.0,
        maximum_jet_mach:0.3,max_jet_history_samples:4096,
    }
}
fn voice(delay: f64, fractional: bool) -> FlueVoice {
    let card=mint_tonal_interim_card().unwrap();
    let gate=CancelGate::new_clock_free();
    if fractional {
        FlueVoice::from_fractional_network(Some(&card),90.0,1e-6,geometry(),source(delay),
            None,policy(),RATE,3,0,&gate).unwrap()
    } else {
        FlueVoice::from_network(Some(&card),90.0,1e-6,geometry(),source(delay),RATE,3,0,&gate).unwrap()
    }
}

#[test]
fn integer_selection_retains_original_jet_equations_and_block_partition() {
    // Power-of-two scaling makes this an exactly integral dimensionless transit.
    let mut old=voice(32.0,false);
    let mut new=voice(32.0,true);
    assert_eq!(new.network().unwrap().fractional_sections().unwrap().1[0].transit.samples(),0.0);
    for v in [&mut old,&mut new] {v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();}
    let mut a=[0.0;513];let mut b=[0.0;513];
    old.step_block(&mut a).unwrap();
    for block in b.chunks_mut(37) {new.step_block(block).unwrap();}
    assert_eq!(a,b);
    assert_eq!(old.network().unwrap().last_frame(),new.network().unwrap().last_frame());
    assert_eq!(old.island(),new.island());
    assert!(old.network().unwrap().fractional_sections().is_none());
}

#[test]
fn subcell_bore_changes_actual_jet_feedback_without_artificial_wave_loss() {
    let mut a=voice(42.2,true);let mut b=voice(42.3,true);
    for v in [&mut a,&mut b] {v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();}
    let mut different=false;
    for n in 0..768 {
        if n==256 {
            a.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
            b.apply(FlueControl::SetBlowingPressure(0.0)).unwrap();
        }
        let mut x=[0.0];let mut y=[0.0];
        a.step_block(&mut x).unwrap();b.step_block(&mut y).unwrap();
        if n<84 {assert_eq!(x,y,"no source response before round-trip return");}
        else {different|=(x[0]-y[0]).abs()>1e-8;}
        let f=a.network().unwrap().last_frame().unwrap();
        let d=a.diagnostics().last().unwrap();
        assert!((f.inlet_work_j-d.pipe_work_j).abs()<1e-12);
        assert!(f.balance_residual_j().abs()<1e-12);
        assert!(d.ledger_defect_j.abs()<1e-12);
    }
    assert!(different,"formerly identical rounded bore transits must affect the jet loop");
    assert!(a.network().unwrap().stored_energy_j()>0.0,"source release must not reset waves");
    assert_eq!(a.network().unwrap().source().graph.sections[0].length_m,
        source(42.2).graph.sections[0].length_m);
}

#[test]
fn fractional_selection_does_not_bypass_card_geometry_or_phase_admission() {
    let card=mint_tonal_interim_card().unwrap();let gate=CancelGate::new_clock_free();
    assert!(matches!(FlueVoice::from_fractional_network(None,90.0,1e-6,geometry(),
        source(42.2),None,policy(),RATE,3,0,&gate),Err(FlueError::NoCard)));
    let mut tight=policy();tight.maximum_summed_phase_error_rad=0.0;
    assert!(FlueVoice::from_fractional_network(Some(&card),90.0,1e-6,geometry(),
        source(42.2),None,tight,RATE,3,0,&gate).is_err());
    let mut bad=geometry();bad.bore_area_m2*=2.0;
    assert!(FlueVoice::from_fractional_network(Some(&card),90.0,1e-6,bad,
        source(42.2),None,policy(),RATE,3,0,&gate).is_err());
    gate.request();
    assert!(FlueVoice::from_fractional_network(Some(&card),90.0,1e-6,geometry(),
        source(42.2),None,policy(),RATE,3,0,&gate).is_err());
}
