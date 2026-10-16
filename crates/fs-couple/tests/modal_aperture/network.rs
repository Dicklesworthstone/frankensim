use super::*;
use fs_couple::bernoulli_aperture::two_sided::{ApertureDuctSpec, LoadedAperture,
    TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, ValveSide, combined_inlet_state};
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::baffled::openings::duct::{DuctOpening, DuctRadiationSpec, DuctReconstruction};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_material::gas::{GasSpec, GasState};
use fs_exec::CancelGate;
use fs_vfit::waveguide::{PassiveWaveguide, WaveguideSpec};

pub(super) fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.).unwrap()
}
pub(super) fn sides() -> [ApertureDuctSpec;2] {
    let g=gas(); let dt=spec().time_step_s;
    let side=|delay:usize,radius:f64,reflection:f64| ApertureDuctSpec {
        graph:TubeNetworkSpec {
            nodes:vec![NetworkNode::Inlet,NetworkNode::Termination{reflection}],
            sections:vec![TubeSection{nodes:[0,1],length_m:delay as f64*(g.sound_speed*dt),
                radius_m:radius,max_length_error_m:1e-12}],
            sound_speed_m_s:g.sound_speed,max_wave_memory_bytes:1<<20,
        },gases:vec![g],viscothermal:vec![],
    };
    [side(12,0.007,0.7),side(32,0.006,-0.8)]
}
fn loaded_spec(count:usize) -> ModalApertureSpec {
    let (g,z)=combined_inlet_state(&sides()).unwrap();
    let mut s=spec();s.modes.truncate(count);s.density_kg_m3=g.density;s.impedance_pa_s_m3=z;s
}
fn network(count:usize) -> TwoSidedAperture {
    TwoSidedAperture::new(model(loaded_spec(count)).into(),sides(),TwoSidedLimits{
        maximum_energy_j:1.,energy_absolute_tolerance_j:1e-16,energy_relative_tolerance:1e-7,
        flow_absolute_tolerance_m3_s:1e-16,flow_relative_tolerance:1e-8,
    },&CancelGate::new_clock_free()).unwrap()
}
fn supply(n:usize) -> TwoSidedDrive {
    TwoSidedDrive {source_pressure_pa:if n<123 {50.}else{0.},
        blade_forces_n:if n<99 {[0.001,-0.0003]}else{[0.;2]} }
}

#[test]
fn complete_modal_source_matches_independent_original_wave_owners() {
    let mut coupled=network(3);let mut direct=model(loaded_spec(3));
    let dt=direct.spec().time_step_s;
    let [up,down]=sides();
    let mut a=PassiveWaveguide::new(WaveguideSpec{one_way_samples:12,
        impedance_pa_s_m3:up.graph.inlet_impedance_with_gases(&up.gases).unwrap(),
        time_step_s:dt,reflection:0.7,max_memory_bytes:1<<20}).unwrap();
    let mut b=PassiveWaveguide::new(WaveguideSpec{one_way_samples:32,
        impedance_pa_s_m3:down.graph.inlet_impedance_with_gases(&down.gases).unwrap(),
        time_step_s:dt,reflection:-0.8,max_memory_bytes:1<<20}).unwrap();
    let gate=CancelGate::new_clock_free();
    for n in 0..512 {
        let supply=supply(n);let ai=a.incoming_pressure_pa();let bi=b.incoming_pressure_pa();
        let expected=direct.step(ApertureDrive{upstream_pressure_pa:supply.source_pressure_pa+2.0*ai,
            incoming_pressure_pa:bi,body_flow_m3_s:0.},supply.blade_forces_n).unwrap();
        let flow=expected.aperture.bore_flow_m3_s;
        a.step(ai-a.spec().impedance_pa_s_m3*flow).unwrap();
        b.step(bi+b.spec().impedance_pa_s_m3*flow).unwrap();
        let actual=coupled.step_under_gate(supply,&gate).unwrap();
        let LoadedAperture::Modal(valve)=coupled.valve() else {panic!("modal source retained")};
        assert_eq!(valve.modal_states(),direct.modal_states());
        assert!((actual.stored_energy_j-direct.stored_energy_j()-a.stored_energy_j()
            -b.stored_energy_j()).abs()<1e-12);
        assert!(actual.balance_residual_j().abs()<1e-12);
    }
    assert!(coupled.stored_energy_j()>0.);
}

#[test]
fn higher_modes_reach_actual_multirate_microphone_without_callback_dependence() {
    let render=|count| {
        TwoSidedPressure::with_radiation(network(count),DuctRadiationSpec{
            openings:vec![DuctOpening{network:1,node:1,center_m:[0.;2]}],exterior:Some(gas()),
            receiver:CircularOutletReceiver{position_m:[0.,0.,0.25],radial_rings:4,
                angular_points:16,maximum_frequency_hz:1000.},
            reconstruction:DuctReconstruction::Differentiating{
                mechanics_rate_hz:96000,output_rate_hz:48000},
        },64).unwrap()
    };
    let mut reference=Vec::new();
    for size in [1,17,64] {
        let mut p=render(3);p.set_drive(supply(0)).unwrap();
        let mut out=vec![0.;512];
        for piece in out.chunks_mut(size) {p.block(piece).unwrap();}
        if reference.is_empty(){reference=out;}else{assert_eq!(out,reference);}
        assert_eq!(p.mechanical_samples_rendered(),1024);
        assert!(p.receiver().filter_delay_samples()>0.);
    }
    assert!(reference.iter().any(|x|x.abs()>1e-12));
    let mut lower=render(2);lower.set_drive(supply(0)).unwrap();let mut out=vec![0.;512];
    for piece in out.chunks_mut(64){lower.block(piece).unwrap();}
    assert!(out.iter().zip(&reference).any(|(a,b)|(a-b).abs()>1e-9));
}

#[test]
fn coupled_refusal_preserves_modal_coordinates_and_future_return_waves() {
    let mut a=network(3);let mut b=network(3);let gate=CancelGate::new_clock_free();
    a.step_under_gate(supply(0),&gate).unwrap();b.step_under_gate(supply(0),&gate).unwrap();
    for value in [f64::NAN,f64::INFINITY,f64::MAX] {
        let mut bad=supply(0);bad.source_pressure_pa=value;
        assert!(a.step_under_gate(bad,&gate).is_err());
        assert_eq!(a.valve().accepted_steps(),1);
        assert_eq!(a.stored_energy_j(),b.stored_energy_j());
    }
    for n in 1..128 {assert_eq!(a.step_under_gate(supply(n),&gate).unwrap(),
        b.step_under_gate(supply(n),&gate).unwrap());}
    assert!(a.node_frame(ValveSide::Downstream,1).is_some());
}
