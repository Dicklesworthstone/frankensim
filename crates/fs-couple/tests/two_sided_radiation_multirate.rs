//! Real single/independent-blade owners, two ducts and the original Rayleigh kernel.
//! Synthetic mechanics isolate data propagation; these are not measured reeds.
use fs_couple::bernoulli_aperture::{BernoulliAperture, dynamic::{
    ApertureState, DynamicAperture, DynamicApertureSpec,
}, network::{NetworkNode, TubeNetworkSpec, TubeSection}, paired::{
    ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits,
    PairedApertureSpec,
}, performance::AperturePerformanceConfig, two_sided::{
    ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedLimits,
    ValveSide, combined_inlet_state, observer::{TwoSidedPressure, schedule::ScheduledTwoSidedPressure},
}};
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_couple::pcm_wav::baffled::openings::{CircularOpening, duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction,
}};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::{ControlDelta, RenderError};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::gesture::GestureSchedule;

const AUDIO: u32 = 48_000;
const MECHANICS: u32 = 192_000;
const OUTPUT: usize = 600;
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),t,101_325.0).unwrap()
}
fn duct(delay: usize, reflection: f64) -> ApertureDuctSpec {
    let g = gas(293.15);
    ApertureDuctSpec { graph: TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
        sections: vec![TubeSection { nodes: [0,1], length_m: delay as f64*g.sound_speed/f64::from(AUDIO),
            radius_m: 0.007, max_length_error_m: 1e-12 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
    }, gases: vec![g], viscothermal: vec![] }
}
fn system(paired: bool, rate: u32, mass_ratio: f64, steps: u64) -> TwoSidedAperture {
    let sides = [duct(12,0.7),duct(32,-0.8)];
    let (g,z) = combined_inlet_state(&sides).unwrap();
    let contact = Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic source for multi-rate radiation integration".into()).unwrap();
    let valve: LoadedAperture = if paired {
        let a = ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 500.0,
            damping_n_s_m: 0.03, pressure_area_m2: 1.25e-4 };
        // Outside the midpoint one-cent band (DEF-04): exact blade blocks.
        PairedAperture::new_with_propagator(PairedApertureSpec {
            rest_opening_m: 4e-4, width_m: 0.013,
            blades: [a,ApertureBlade { mass_kg: a.mass_kg*mass_ratio,..a }],
            density_kg_m3: g.density, impedance_pa_s_m3: z,
            time_step_s: 1.0/f64::from(rate), max_steps: steps,
            limits: PairedApertureLimits { maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
                max_evaluations: 256, velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
        },[BladeState::default();2],contact,BladePropagator::ExactZoh).unwrap().into()
    } else {
        DynamicAperture::new_with_propagator(DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 1000.0 },
            mass_kg: 1e-5*mass_ratio, stiffness_n_m: 500.0, damping_ratio: 0.05,
            density_kg_m3: g.density, impedance_pa_s_m3: z,
            time_step_s: 1.0/f64::from(rate), max_steps: steps,
        },ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 },contact,
            BladePropagator::ExactZoh).unwrap().into()
    };
    TwoSidedAperture::new(valve,sides,TwoSidedLimits { maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-16, energy_relative_tolerance: 1e-7,
        flow_absolute_tolerance_m3_s: 1e-16, flow_relative_tolerance: 1e-8 },
        &CancelGate::new_clock_free()).unwrap()
}
fn spec() -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: vec![DuctOpening { network: 0,node: 1,center_m: [-0.025,0.0] },
            DuctOpening { network: 1,node: 1,center_m: [0.025,0.0] }],
        exterior: Some(gas(310.0)),
        receiver: CircularOutletReceiver { position_m: [0.08,0.0,0.3], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        reconstruction: DuctReconstruction::Differentiating {
            mechanics_rate_hz: MECHANICS,output_rate_hz: AUDIO },
    }
}
fn independent_receiver(s: &DuctRadiationSpec) -> BaffledPressure {
    let geometry: Vec<_> = s.openings.iter().map(|p| CircularOpening {
        center_m: p.center_m,radius_m: 0.007,
    }).collect();
    let g = s.exterior.unwrap();
    BaffledPressure::circular_openings_multirate_differentiating(&geometry,MECHANICS,AUDIO,
        s.receiver,RayleighMedium { density: g.density,sound_speed: g.sound_speed }).unwrap()
}
fn renderer(paired: bool, mass: f64, steps: u64) -> TwoSidedPressure {
    TwoSidedPressure::with_radiation(system(paired,MECHANICS,mass,steps),spec(),128).unwrap()
}
fn drive(paired: bool, pressure: f64) -> TwoSidedDrive {
    TwoSidedDrive { source_pressure_pa: pressure, blade_forces_n: [1e-5,if paired {-2e-5} else {0.0}] }
}

#[test]
fn single_and_independent_blade_physics_reach_both_outlets_without_skipping_substeps() {
    for paired in [false,true] {
        let mut source = system(paired,MECHANICS,1.3,2400);
        let mut output = renderer(paired,1.3,2400);
        output.set_drive(drive(paired,50.0)).unwrap();
        assert_eq!(output.radiation(),&spec());
        let mut receiver = independent_receiver(&spec());
        let gate = CancelGate::new_clock_free();
        let mut heard = false;
        for i in 0..OUTPUT {
            let mut trace = [0.0;8];
            let mut last = None;
            for step in 0..4 {
                last = Some(source.step_under_gate(drive(paired,50.0),&gate).unwrap());
                for (side,offset) in [(ValveSide::Upstream,0),(ValveSide::Downstream,1)] {
                    trace[2*step+offset] = source.node_frame(side,1).unwrap().net_flow_into_node_m3_s;
                }
            }
            let expected = receiver.step_trace(&trace).unwrap();
            let mut actual = [0.0];output.block(&mut actual).unwrap();
            assert_eq!(actual[0],expected);
            assert_eq!(output.last_frame(),last.as_ref());
            assert_eq!(output.mechanical_samples_rendered(),4*(i as u64+1));
            heard |= actual[0].abs()>1e-9;
        }
        assert!(heard);
        assert_eq!(output.samples_rendered(),OUTPUT as u64);
        assert!(output.validate_sample_count(1).is_err());
        assert_eq!(output.system().stored_energy_j(),source.stored_energy_j());
    }
}

fn phrase() -> GestureSchedule {
    GestureSchedule::from_canonical_bytes(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t44100\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e1\nevent\t3e-3\t0e0\tpressure-pa\t0e0\n").unwrap()
}
fn scheduled() -> ScheduledTwoSidedPressure {
    let mut p = renderer(true,1.3,2400);
    p.set_drive(drive(true,0.0)).unwrap();
    ScheduledTwoSidedPressure::new(p,phrase(),AperturePerformanceConfig {
        sample_rate_hz: AUDIO,samples: OUTPUT as u64,max_block: 128,
        max_compile_work: 100_000,max_controls: 1000,
    }).unwrap()
}

#[test]
fn controls_inside_audio_packets_match_direct_mechanics_and_survive_callbacks() {
    let mut scheduled = scheduled();
    let controls = scheduled.pending_controls().to_vec();
    assert!(controls.iter().any(|c| c.sample%4!=0));
    let mut direct = system(true,MECHANICS,1.3,2400);
    let mut receiver = independent_receiver(&spec());
    let gate = CancelGate::new_clock_free();
    let mut held = 0.0;
    let mut expected = Vec::new();
    let mut last = None;
    for frame in 0..OUTPUT {
        let mut trace = [0.0;8];
        for sub in 0..4 {
            let sample = (4*frame+sub) as u64;
            for c in controls.iter().filter(|c| c.sample==sample) {
                if let ControlDelta::SetBlowingPressure { pressure_pa,.. } = c.delta { held=pressure_pa; }
            }
            last = Some(direct.step_under_gate(drive(true,held),&gate).unwrap());
            trace[2*sub] = direct.node_frame(ValveSide::Upstream,1).unwrap().net_flow_into_node_m3_s;
            trace[2*sub+1] = direct.node_frame(ValveSide::Downstream,1).unwrap().net_flow_into_node_m3_s;
        }
        expected.push(receiver.step_trace(&trace).unwrap());
    }
    let mut actual = vec![0.0;OUTPUT];
    for piece in actual.chunks_mut(37) { scheduled.block(piece).unwrap(); }
    assert_eq!(actual,expected);
    assert_eq!(scheduled.source().last_frame(),last.as_ref());
    assert_eq!(scheduled.applied_controls(),&controls);
    assert_eq!(scheduled.remaining_samples(),0);
    assert_eq!(scheduled.source().last_frame().unwrap().source_work_j,0.0);
    for chunk in [1,128] {
        let mut p = self::scheduled();let mut actual=vec![0.0;OUTPUT];
        for piece in actual.chunks_mut(chunk) { p.block(piece).unwrap(); }
        assert_eq!(actual,expected);
        assert_eq!(p.source().last_frame(),scheduled.source().last_frame());
    }
}

#[test]
fn lifetime_clock_and_failure_gates_apply_to_mechanics_not_output_count() {
    let mut p = renderer(true,1.3,9); // Two full audio frames; one unusable substep.
    let mut too_long = [987.0;3];assert!(p.block(&mut too_long).is_err());
    assert_eq!(too_long,[987.0;3]);assert_eq!(p.mechanical_samples_rendered(),0);
    p.block(&mut [0.0;2]).unwrap();
    assert_eq!(p.mechanical_samples_rendered(),8);assert!(p.validate_sample_count(1).is_err());
    assert!(p.validate_sample_count(u64::MAX).is_err());
    assert!(p.validate_sample_rate(MECHANICS).is_err());
    assert!(TwoSidedPressure::with_radiation(system(true,AUDIO,1.3,2400),spec(),128).is_err());
    let mut p = renderer(false,1.3,2400);
    assert!(p.set_drive(drive(true,0.0)).is_err()); // No fictitious second force port.
    let mut p = self::scheduled();let pending = p.pending_controls().len();
    p.set_blade_forces([f64::MAX,0.0]).unwrap();
    let mut untouched = [987.0;2];assert!(p.block(&mut untouched).is_err());
    assert_eq!(untouched,[987.0;2]);assert_eq!(p.pending_controls().len(),pending);
    assert!(matches!(p.validate_sample_count(0),Err(RenderError::Poisoned)));
}

#[test]
fn independent_blade_mass_data_changes_real_pressure_not_an_observer_gain() {
    let mut signals = Vec::new();
    for mass in [1.0,2.0] {
        let mut p = renderer(true,mass,2400);p.set_drive(drive(true,50.0)).unwrap();
        let mut signal=vec![0.0;OUTPUT];for block in signal.chunks_mut(128) {p.block(block).unwrap();}
        if let LoadedAperture::Paired(a) = p.system().valve() {
            assert_eq!(a.spec().blades[1].mass_kg,1e-5*mass);
        } else { panic!("independent blades were lost"); }
        signals.push(signal);
    }
    assert_ne!(signals[0],signals[1]);
}
