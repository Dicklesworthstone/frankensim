//! Original reduced plates -> source-owned force ports -> both ducts -> pressure.
//! Numerical fixtures exercise coupling, not measured reed or lip/tongue fidelity.
use super::*;
use fs_couple::bernoulli_aperture::{dynamic::{ApertureState, DynamicAperture},
    dynamic::force::{PlateForceFootprint, PlateForcePort},
    network::{NetworkNode, TubeNetworkSpec, TubeSection},
    paired::{BladeState, PairedAperture, PairedApertureLimits},
    performance::{AperturePerformanceConfig, force::{ApertureForceEvent, ApertureForceProgram}},
    two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive,
        TwoSidedLimits, ValveSide, combined_inlet_state,
        observer::{TwoSidedPressure, schedule::ScheduledTwoSidedPressure}}};
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_couple::pcm_wav::baffled::openings::{CircularOpening, duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction}};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::{ControlDelta, RenderError};
use fs_dcontact::Obstacle;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const RATE: u32 = 48_000;
const MECHANICS: u32 = 192_000;
const SAMPLES: usize = 600;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap()
}
fn model(plates: &[PlateApertureReduction; 2], paired: bool) -> TwoSidedAperture {
    let g = gas();
    let sides = [12,32].map(|delay| ApertureDuctSpec { graph: TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }],
        sections: vec![TubeSection { nodes: [0,1], length_m: delay as f64*g.sound_speed/f64::from(RATE),
            radius_m: 0.007, max_length_error_m: 1e-12 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
    }, gases: vec![g], viscothermal: vec![] });
    let (g,z) = combined_inlet_state(&sides).unwrap();
    let lay = Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic source plate force integration".into()).unwrap();
    let dt = 1.0/f64::from(MECHANICS);
    let valve: LoadedAperture = if paired {
        PairedAperture::from_plates((*plates).clone(),g,z,dt,2400,
            PairedApertureLimits { maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
                max_evaluations: 256, velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
            [BladeState::default();2],lay).unwrap().into()
    } else {
        DynamicAperture::from_plate(plates[0].clone(),g.density,z,dt,2400,
            ApertureState { opening_m: plates[0].options().rest_opening_m,
                opening_velocity_m_s: 0.0 },lay).unwrap().into()
    };
    TwoSidedAperture::new(valve,sides,TwoSidedLimits { maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-16,energy_relative_tolerance: 1e-7,
        flow_absolute_tolerance_m3_s: 1e-16,flow_relative_tolerance: 1e-8 },
        &CancelGate::new_clock_free()).unwrap()
}
fn radiation() -> DuctRadiationSpec {
    DuctRadiationSpec { openings: vec![
        DuctOpening { network: 0,node: 1,center_m: [-0.025,0.0] },
        DuctOpening { network: 1,node: 1,center_m: [0.025,0.0] }],
        exterior: Some(gas()),
        receiver: CircularOutletReceiver { position_m: [0.08,0.0,0.3],radial_rings: 4,
            angular_points: 16,maximum_frequency_hz: 1000.0 },
        reconstruction: DuctReconstruction::Differentiating {
            mechanics_rate_hz: MECHANICS,output_rate_hz: RATE },
    }
}
fn config() -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz: RATE,samples: SAMPLES as u64,
        max_block: 128,max_compile_work: 100_000,max_controls: 1024 }
}
fn phrase() -> GestureSchedule {
    GestureSchedule::try_new(700,vec![GestureTrack { id: "blow".into(),
        target: GestureTarget::BlowingPressure,initial: GestureValue::PressurePa(3.0),
        events: vec![GestureEvent { time_s: 3.0/700.0,transition_s: 0.0,
            value: GestureValue::PressurePa(0.0) }] }]).unwrap()
}
fn player(plates: &[PlateApertureReduction;2], paired: bool, config: AperturePerformanceConfig)
    -> ScheduledTwoSidedPressure {
    let source = TwoSidedPressure::with_radiation(model(plates,paired),radiation(),128).unwrap();
    ScheduledTwoSidedPressure::new(source,phrase(),config).unwrap()
}
fn program(plate: &PlateApertureReduction, scale: f64) -> ApertureForceProgram {
    let event = |sample,port,force| ApertureForceEvent { sample,port,force_n: scale*force };
    ApertureForceProgram { footprints: vec![
        PlateForceFootprint::Node(plate.options().slit_edges[0][0]),
        PlateForceFootprint::Patch((0..plate.chart().mesh.tris.len()).collect())],
        // Unsorted, overlapping loads, duplicate assignments at the same sample,
        // and a valid mechanical event LATER than the number of output frames.
        events: vec![event(1001,0,0.0),event(1,0,1e-5),event(3,1,-2e-5),
            event(1,0,2e-5),event(601,0,-1e-5),event(801,1,0.0)] }
}
fn ports(plate: &PlateApertureReduction, program: &ApertureForceProgram) -> Vec<PlateForcePort> {
    program.footprints.iter().map(|f| plate.force_port(f.clone()).unwrap()).collect()
}
fn plates() -> [PlateApertureReduction;2] {
    [reduction(4e9,900.0),reduction(4e9,1200.0)]
}

#[test]
fn physical_plate_loads_on_single_and_independent_blades_match_direct_substep_physics() {
    let plates = plates();
    let programs = [program(&plates[0],1.0),program(&plates[1],-0.7)];
    let projected = [ports(&plates[0],&programs[0]),ports(&plates[1],&programs[1])];
    for paired in [false,true] {
        let count = if paired {2} else {1};
        let mut source = model(&plates,paired);
        let mut p = player(&plates,paired,config());
        let initial_energy = p.source().system().stored_energy_j();
        let pressure_controls = p.pending_controls().to_vec();
        for (i, program) in programs[..count].iter().enumerate() {
            p = p.with_plate_forces(i,program.clone()).unwrap();
            for (a,b) in p.plate_force_ports(i).iter().zip(&projected[i]) {
                assert_eq!(a.coefficient(),b.coefficient());
                assert_eq!(a.footprint(),b.footprint());
            }
        }
        assert_eq!(p.source().mechanical_samples_rendered(),0);
        assert_eq!(p.source().system().stored_energy_j(),initial_energy);
        assert!(p.set_blade_forces([0.0;2]).is_err()); // No silent override.
        let s = radiation();let g = gas();
        let geometry: Vec<_> = s.openings.iter().map(|o| CircularOpening {
            center_m: o.center_m,radius_m: 0.007 }).collect();
        let mut mic = BaffledPressure::circular_openings_multirate_differentiating(
            &geometry,MECHANICS,RATE,s.receiver,
            RayleighMedium { density: g.density,sound_speed: g.sound_speed }).unwrap();
        let gate = CancelGate::new_clock_free();
        let mut held = [[0.0;2];2];let mut pressure = 0.0;
        let mut expected = Vec::new();let mut last = None;let mut work = 0.0;
        for frame in 0..SAMPLES {
            let mut trace = [0.0;8];
            for sub in 0..4 {
                let sample = (4*frame+sub) as u64;
                for c in pressure_controls.iter().filter(|c| c.sample == sample) {
                    if let ControlDelta::SetBlowingPressure { pressure_pa,.. } = c.delta { pressure = pressure_pa; }
                }
                let mut forces = [0.0;2];
                for i in 0..count {
                    for e in programs[i].events.iter().filter(|e| e.sample == sample) { held[i][e.port] = e.force_n; }
                    // Separate per-sample physical projection, not the schedule compiler.
                    for (port,force) in projected[i].iter().zip(held[i]) {
                        forces[i] += port.generalized_force_n(force).unwrap();
                    }
                }
                let f = source.step_under_gate(TwoSidedDrive {
                    source_pressure_pa: pressure,blade_forces_n: forces },&gate).unwrap();
                assert!(f.balance_residual_j().abs() < 1e-11);
                work += f.aperture.mechanical_work_j.abs();
                last = Some(f);
                for (i,side) in [ValveSide::Upstream,ValveSide::Downstream].into_iter().enumerate() {
                    trace[2*sub+i] = source.node_frame(side,1).unwrap().net_flow_into_node_m3_s;
                }
            }
            expected.push(mic.step_trace(&trace).unwrap());
        }
        assert!(work > 0.0);
        assert!(expected.iter().any(|v| v.abs() > 1e-12));
        let mut actual = vec![0.0;SAMPLES];
        for piece in actual.chunks_mut(37) { p.block(piece).unwrap(); }
        assert_eq!(actual,expected);
        assert_eq!(p.source().last_frame(),last.as_ref());
        assert_eq!(p.held_generalized_forces_n(),[0.0;2]);
        for (i,program) in programs[..count].iter().enumerate() {
            let mut sorted = program.events.clone();sorted.sort_by_key(|e|e.sample);
            assert_eq!(p.applied_plate_force_controls(i),sorted.as_slice());
            assert!(p.pending_plate_force_controls(i).is_empty());
        }
        for chunk in [1,128] {
            let mut p = player(&plates,paired,config());
            for (i,program) in programs[..count].iter().enumerate() {
                p = p.with_plate_forces(i,program.clone()).unwrap();
            }
            let mut actual = vec![0.0;SAMPLES];
            for piece in actual.chunks_mut(chunk) { p.block(piece).unwrap(); }
            assert_eq!(actual,expected);
            assert_eq!(p.source().last_frame(),last.as_ref());
        }
    }
}

#[test]
fn bad_attachments_missing_geometry_and_combined_budgets_refuse_cold() {
    let plates = plates();let a = program(&plates[0],1.0);let b = program(&plates[1],1.0);
    assert!(player(&plates,false,config()).with_plate_forces(1,b.clone()).is_err());
    assert!(player(&plates,true,config()).with_plate_forces(2,b.clone()).is_err());
    assert!(player(&plates,true,config()).with_plate_forces(0,a.clone()).unwrap()
        .with_plate_forces(0,a.clone()).is_err());
    let mut p = player(&plates,true,config());p.set_blade_forces([1e-5,0.0]).unwrap();
    assert!(p.with_plate_forces(0,a.clone()).is_err());
    let mut p = player(&plates,true,config());p.block(&mut [0.0]).unwrap();
    assert!(p.with_plate_forces(0,a.clone()).is_err());
    let mut wrong = a.clone();wrong.events[0].sample = 2400;
    assert!(player(&plates,true,config()).with_plate_forces(0,wrong).is_err());
    let mut wrong = a.clone();wrong.events[0].force_n = f64::NAN;
    assert!(player(&plates,true,config()).with_plate_forces(0,wrong).is_err());
    let mut wrong = a.clone();wrong.footprints[0] = PlateForceFootprint::Node(usize::MAX);
    assert!(player(&plates,true,config()).with_plate_forces(0,wrong).is_err());
    let mut small = config();small.max_controls = 10;
    let p = player(&plates,true,small).with_plate_forces(0,a.clone()).unwrap();
    assert!(p.with_plate_forces(1,b.clone()).is_err()); // Combined, not per blade.
    let mut small = config();small.max_compile_work = 80;
    let p = player(&plates,true,small).with_plate_forces(0,a.clone()).unwrap();
    assert!(p.with_plate_forces(1,b).is_err()); // Each costs 24 + 1 + 6*3 = 43.
    // A mechanically valid scalar valve without a retained mesh must not get fake ports.
    let physical = model(&plates,false);
    let LoadedAperture::Single(valve) = physical.valve() else { unreachable!() };
    let bare = DynamicAperture::new(*valve.spec(),valve.state(),valve.contact_law().clone()).unwrap();
    let bare = TwoSidedAperture::new(bare.into(),
        [physical.source(ValveSide::Upstream).clone(),physical.source(ValveSide::Downstream).clone()],
        TwoSidedLimits { maximum_energy_j: 1.0,energy_absolute_tolerance_j: 1e-16,
            energy_relative_tolerance: 1e-7,flow_absolute_tolerance_m3_s: 1e-16,flow_relative_tolerance: 1e-8 },
        &CancelGate::new_clock_free()).unwrap();
    let bare = TwoSidedPressure::with_radiation(bare,radiation(),128).unwrap();
    assert!(ScheduledTwoSidedPressure::new(bare,phrase(),config()).unwrap().with_plate_forces(0,a).is_err());
}

#[test]
fn failed_substep_does_not_publish_future_plate_loads_or_an_incomplete_audio_frame() {
    let plates = plates();
    let mut program = program(&plates[0],1.0);
    program.events = vec![ApertureForceEvent { sample: 0,port: 0,force_n: 1e-5 },
        ApertureForceEvent { sample: 1,port: 0,force_n: 1e100 },
        ApertureForceEvent { sample: 2,port: 0,force_n: 0.0 }];
    let mut p = player(&plates,true,config()).with_plate_forces(0,program.clone()).unwrap();
    let first = p.plate_force_ports(0)[0].generalized_force_n(1e-5).unwrap();
    let mut untouched = [987.0;2];assert!(p.block(&mut untouched).is_err());
    assert_eq!(untouched,[987.0;2]);assert_eq!(p.samples_rendered(),0);
    assert_eq!(p.source().mechanical_samples_rendered(),1);
    assert_eq!(p.applied_plate_force_controls(0),&program.events[..1]);
    assert_eq!(p.pending_plate_force_controls(0),&program.events[1..]);
    assert_eq!(p.held_generalized_forces_n(),[first,0.0]);
    assert!(matches!(p.validate_sample_count(0),Err(RenderError::Poisoned)));
}

#[test]
fn a_stationary_support_cannot_inject_sound_and_symmetric_ports_keep_their_original_projection() {
    let original = plates();
    for symmetric in [false,true] {
        let mut plates = original.clone();
        if symmetric { plates[0] = plates[0].clone().into_symmetric_opposed_pair().unwrap(); }
        let root = plates[0].chart().mesh.nodes.iter().position(|p| p.0 == 0.0).unwrap();
        let tip = plates[0].options().slit_edges[0][0];
        for (node, audible) in [(root,false),(tip,true)] {
            let quiet = GestureSchedule::try_new(700,vec![GestureTrack { id: "blow".into(),
                target: GestureTarget::BlowingPressure,initial: GestureValue::PressurePa(0.0),
                events: vec![] }]).unwrap();
            let source = TwoSidedPressure::with_radiation(model(&plates,false),radiation(),128).unwrap();
            let mut p = ScheduledTwoSidedPressure::new(source,quiet,config()).unwrap()
                .with_plate_forces(0,ApertureForceProgram { footprints: vec![PlateForceFootprint::Node(node)],
                    events: vec![ApertureForceEvent { sample: 1,port: 0,force_n: 1e-5 },
                        ApertureForceEvent { sample: 17,port: 0,force_n: 0.0 }] }).unwrap();
            let port = &p.plate_force_ports(0)[0];
            let expected = plates[0].force_port(PlateForceFootprint::Node(node)).unwrap();
            assert_eq!(port.coefficient(),expected.coefficient());
            if symmetric {
                let single = original[0].force_port(PlateForceFootprint::Node(node)).unwrap();
                assert_eq!(port.coefficient(),0.5*single.coefficient());
            }
            let mut output = vec![0.0;SAMPLES];
            for block in output.chunks_mut(37) { p.block(block).unwrap(); }
            assert_eq!(output.iter().any(|x| x.abs()>1e-15),audible);
            if !audible {
                assert!(output.iter().all(|x| *x==0.0));
                assert_eq!(p.source().system().stored_energy_j(),0.0);
            }
        }
    }
}
