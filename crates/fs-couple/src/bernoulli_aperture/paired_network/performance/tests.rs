use crate::render::{ControlDelta, schedule::{PressureGestureBinding, compile_pressure_gestures}};
use super::*;
use crate::bernoulli_aperture::paired::{ApertureBlade, BladeState, PairedAperture, PairedApertureLimits, PairedApertureSpec};
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::pcm_wav::baffled::{CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::DecimatedRenderer;
use crate::pcm_wav::decimate::Decimator;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
const SAMPLES: u64 = 513;
fn config() -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz: RATE, samples: SAMPLES, max_block: 128,
        max_compile_work: 100_000, max_controls: 128 }
}
fn schedule() -> GestureSchedule {
    GestureSchedule::from_canonical_bytes(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t3\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t0e0\tpressure-pa\t0e0\nevent\t7e-3\t0e0\tpressure-pa\t2e0\n").unwrap()
}
fn forces() -> Vec<BladeForceControl> {
    vec![BladeForceControl { sample: 3, blade: 0, force_n: 1e-5 },
        BladeForceControl { sample: 3, blade: 1, force_n: -2e-5 },
        BladeForceControl { sample: 150, blade: 1, force_n: 0.0 },
        BladeForceControl { sample: 151, blade: 0, force_n: 0.0 }]
}
fn network() -> PairedApertureNetwork {
    let gas = |t| GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101325.0).unwrap();
    let gases = vec![gas(293.15), gas(323.15)];
    let graph = TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.8 }],
        sections: vec![TubeSection { nodes: [0,1], length_m: 0.12, radius_m: 0.006, max_length_error_m: 0.005 },
            TubeSection { nodes: [1,2], length_m: 0.18, radius_m: 0.007, max_length_error_m: 0.005 }],
        sound_speed_m_s: gases[0].sound_speed, max_wave_memory_bytes: 1<<20 };
    let spec = PairedApertureSpec { rest_opening_m: 4e-4, width_m: 0.013,
        blades: [ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 350.0, damping_n_s_m: 0.03, pressure_area_m2: 1e-4 },
            ApertureBlade { mass_kg: 2e-5, stiffness_n_m: 700.0, damping_n_s_m: 0.04, pressure_area_m2: 1.2e-4 }],
        density_kg_m3: gases[0].density, impedance_pa_s_m3: graph.inlet_impedance_with_gases(&gases).unwrap(),
        time_step_s: 1.0/f64::from(RATE), max_steps: SAMPLES,
        limits: PairedApertureLimits { maximum_abs_displacement_m: 0.05, maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 10.0, max_evaluations: 256, velocity_tolerance_m_s: 1e-15, relative_tolerance: 1e-13 } };
    let contact = Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],2e9,2.0,"authored paired performance contact".into()).unwrap();
    let aperture = PairedAperture::new(spec,[BladeState::default();2],contact).unwrap();
    PairedApertureNetwork::new(aperture,graph,gases,&CancelGate::new_clock_free()).unwrap()
}
fn renderer(observation: ApertureObservation) -> PairedAperturePerformance {
    PairedAperturePerformance::new(network(),observation,schedule(),config()).unwrap()
        .with_blade_forces(forces()).unwrap()
}
fn render(mut p: PairedAperturePerformance, size: usize) -> (PairedAperturePerformance, Vec<f64>) {
    let mut out=vec![0.0;SAMPLES as usize];
    for block in out.chunks_mut(size) { p.block(block).unwrap(); }
    (p,out)
}

#[test]
fn pressure_and_independent_force_events_match_direct_physics_at_exact_samples() {
    let mut p=renderer(ApertureObservation::Inlet);
    let mut direct=network();
    let controls=compile_pressure_gestures(&schedule(),&[PressureGestureBinding {
        track: "blow".into(), voice: 0,
    }],RATE,SAMPLES,100_000).unwrap();
    let f=forces(); let mut pressure=0.0; let mut held=[0.0;2];
    for sample in 0..SAMPLES {
        for e in controls.iter().filter(|e|e.sample==sample) {
            if let ControlDelta::SetBlowingPressure {pressure_pa,..}=e.delta {pressure=pressure_pa;}
        }
        for e in f.iter().filter(|e|e.sample==sample) {held[e.blade]=e.force_n;}
        let expected=direct.step_with_forces(TubeDrive {upstream_pressure_pa:pressure,body_flow_m3_s:0.0},held).unwrap();
        let mut out=[f64::NAN];p.block(&mut out).unwrap();
        assert_eq!(out[0],expected.mechanics.aperture.bore_pressure_pa);
        assert_eq!(*p.last_frame().unwrap(),expected);
        assert_eq!(p.system().aperture().blade_states(),direct.aperture().blade_states());
    }
    assert!(p.pending_controls().is_empty() && p.pending_force_controls().is_empty());
    assert_eq!(p.applied_force_controls(),f);
}

#[test]
fn arbitrary_callbacks_preserve_the_entire_pair_and_release_history() {
    let (reference,expected)=render(renderer(ApertureObservation::NetworkNode(2)),1);
    for size in [7,64,128] {
        let (p,out)=render(renderer(ApertureObservation::NetworkNode(2)),size);
        assert_eq!(out,expected);
        assert_eq!(p.last_frame(),reference.last_frame());
        assert_eq!(p.system().aperture().blade_states(),reference.system().aperture().blade_states());
        assert_eq!(p.remaining_samples(),0);
        assert!(p.validate_sample_count(1).is_err());
    }
    assert!(expected.iter().any(|p|*p!=0.0));
    assert!(reference.system().stored_energy_j()>0.0);
}

#[test]
fn exterior_observation_uses_real_outlet_flow_gas_and_the_original_receiver() {
    let location=CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,
        angular_points:16,maximum_frequency_hz:2000.0};
    let mut external=renderer(ApertureObservation::NetworkBaffled {node:2,receiver:location});
    let mut internal=renderer(ApertureObservation::Inlet);
    let gas=internal.system().section_gases()[1];
    let mut direct=BaffledPressure::circular_outlet(0.007,RATE,location,
        RayleighMedium {density:gas.density,sound_speed:gas.sound_speed}).unwrap();
    let mut heard=false;
    for _ in 0..SAMPLES {
        internal.block(&mut [0.0]).unwrap();
        let flow=internal.system().node_frame(2).unwrap().net_flow_into_node_m3_s;
        let expected=direct.step(&[flow]).unwrap();
        let mut out=[0.0];external.block(&mut out).unwrap();
        assert_eq!(out[0].to_bits(),expected.to_bits());heard|=out[0]!=0.0;
        assert_eq!(external.last_frame(),internal.last_frame());
    }
    assert!(heard);
    assert_eq!(external.baffled_receiver().unwrap().delay_samples,direct.delay_samples);
}

#[test]
fn shared_decimator_changes_only_observation_not_force_event_timing() {
    for ratio in [1,2,3,4] {
        let mut direct=renderer(ApertureObservation::Inlet);
        let source=renderer(ApertureObservation::Inlet);
        let mut decimated=DecimatedRenderer::new(source,RATE,RATE/ratio,64).unwrap();
        let mut filter=Decimator::new(ratio as usize,1).unwrap();
        assert_eq!(decimated.info().delay_output_samples,filter.delay_output_frames());
        for _ in 0..(SAMPLES/u64::from(ratio)) {
            let mut input=vec![0.0;ratio as usize];direct.block(&mut input).unwrap();
            let expected=filter.preview(&input).unwrap()[0];filter.commit();
            let mut actual=[0.0];decimated.block(&mut actual).unwrap();
            assert_eq!(actual[0].to_bits(),expected.to_bits());
            assert_eq!(decimated.source().last_frame(),direct.last_frame());
        }
    }
}

#[test]
fn bad_requests_leave_state_and_outputs_untouched_and_failed_physics_cannot_resume() {
    let make=||PairedAperturePerformance::new(network(),ApertureObservation::Inlet,schedule(),config()).unwrap();
    let mut p=make();let mut out=[123.0;129];
    assert!(p.block(&mut out).is_err());assert_eq!(out,[123.0;129]);
    assert_eq!(p.samples_rendered(),0);assert!(p.block(&mut []).is_err());
    assert!(p.validate_sample_rate(RATE/2).is_err());
    for bad in [
        vec![BladeForceControl {sample:1,blade:2,force_n:0.0}],
        vec![BladeForceControl {sample:SAMPLES,blade:0,force_n:0.0}],
        vec![BladeForceControl {sample:1,blade:0,force_n:f64::NAN}],
        vec![forces()[0],forces()[0]],vec![forces()[1],forces()[0]],
    ] {assert!(make().with_blade_forces(bad).is_err());}
    let mut p=make().with_blade_forces(vec![BladeForceControl {sample:2,blade:1,force_n:f64::MAX}]).unwrap();
    let mut prefix=[123.0;4];assert!(p.block(&mut prefix).is_err());
    assert_eq!(p.samples_rendered(),2);assert_eq!(&prefix[2..],&[123.0,123.0]);
    assert!(p.applied_force_controls().is_empty());assert_eq!(p.pending_force_controls().len(),1);
    assert!(matches!(p.validate_sample_count(0),Err(RenderError::Poisoned)));
    let mut sentinel=[42.0];assert!(p.block(&mut sentinel).is_err());assert_eq!(sentinel,[42.0]);
    assert!(PairedAperturePerformance::new(network(),ApertureObservation::TubeTerminal,schedule(),config()).is_err());
    assert!(PairedAperturePerformance::new(network(),ApertureObservation::NetworkNode(99),schedule(),config()).is_err());
    let mut bad=config();bad.sample_rate_hz/=2;
    assert!(PairedAperturePerformance::new(network(),ApertureObservation::Inlet,schedule(),bad).is_err());
    bad=config();bad.max_controls=1;
    assert!(PairedAperturePerformance::new(network(),ApertureObservation::Inlet,schedule(),bad).is_err());
}
