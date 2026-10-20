//! Key programs reach the existing physical receivers and sole PCM encoder.
use super::*;
use super::keyed::{geometry, key, key_plate};
use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use fs_couple::bernoulli_aperture::network::boundaries::BoundaryDrive;
use fs_couple::bernoulli_aperture::paired_network::performance::{
    BladeForceControl, PairedAperturePerformance,
};
use fs_couple::bernoulli_aperture::performance::{
    ApertureObservation, AperturePerformance, AperturePerformanceConfig, CoupledAperture,
    boundaries::{BoundaryControl, BoundaryProgram},
    force::{ApertureForceEvent, ApertureForceProgram},
};
use fs_couple::bernoulli_aperture::plate::PlateApertureReduction;
use fs_couple::pcm_wav::{encode_pcm16_wav,
    baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium},
    observation::PressureRenderer};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const RATE: u32 = 100_000;
const SAMPLES: u64 = 2048;
fn phrase() -> GestureSchedule {
    GestureSchedule::try_new(1, vec![GestureTrack { id: "mouth".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(100.0),
        events: vec![] }]).unwrap()
}
fn config(block: usize, controls: usize) -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz: RATE, samples: SAMPLES, max_block: block,
        max_compile_work: 10_000, max_controls: controls }
}
fn location() -> CircularOutletReceiver {
    CircularOutletReceiver { position_m: [0.0, 0.0, 0.2], radial_rings: 4,
        angular_points: 16, maximum_frequency_hz: 2000.0 }
}
fn observation() -> ApertureObservation {
    ApertureObservation::NetworkBaffled { node: 2, receiver: location() }
}
fn microphone() -> BaffledPressure {
    let (graph, gases) = geometry();
    BaffledPressure::circular_outlet(graph.sections[1].radius_m, RATE, location(),
        RayleighMedium { density: gases[1].density, sound_speed: gases[1].sound_speed }).unwrap()
}
fn boundary_at(sample: u64, scale: f64) -> BoundaryDrive {
    BoundaryDrive { node: 3, exterior_pressure_pa: if (257..1001).contains(&sample) { 1.0 } else { 0.0 },
        opening_force_n: scale * match sample { 0..257 => -0.001, 257..1001 => -0.03,
            1001..1536 => -0.015, _ => 0.0 } }
}
fn program(scale: f64) -> BoundaryProgram {
    BoundaryProgram { initial: vec![BoundaryDrive { node: 3, exterior_pressure_pa: 0.0, opening_force_n: 0.0 }],
        // Cold sorting must preserve the physical sample address, including zero.
        controls: [1536, 0, 1001, 257].into_iter().map(|sample|
            BoundaryControl { sample, drive: boundary_at(sample, scale) }).collect() }
}
fn single(plate: &PlateApertureReduction, block: usize, controls: usize) -> AperturePerformance {
    AperturePerformance::new(CoupledAperture::Network(keyed::make(plate, SAMPLES)),
        observation(), phrase(), config(block, controls)).unwrap()
}
fn paired(plate: &PlateApertureReduction, block: usize, controls: usize) -> PairedAperturePerformance {
    PairedAperturePerformance::new(keyed_paired::make(plate, SAMPLES),
        observation(), phrase(), config(block, controls)).unwrap()
}
fn blade_forces() -> Vec<BladeForceControl> {
    vec![BladeForceControl { sample: 0, blade: 0, force_n: 0.001 },
        BladeForceControl { sample: 0, blade: 1, force_n: -0.0007 },
        BladeForceControl { sample: 701, blade: 0, force_n: 0.0 },
        BladeForceControl { sample: 701, blade: 1, force_n: 0.0 }]
}
fn pressure(mut source: impl PressureRenderer, block: usize) -> Vec<f64> {
    let mut result = vec![0.0; SAMPLES as usize];
    for chunk in result.chunks_mut(block) { source.block(chunk).unwrap(); }
    result
}

#[test]
fn single_key_program_matches_direct_receiver_and_material_force_changes_reach_pcm() {
    let plate = key_plate();
    let mut reference = Vec::new();
    for block in [1, 37, 128] {
        let mut renderer = single(&plate, block, 64).with_boundary_controls(program(1.0)).unwrap();
        assert!(renderer.applied_boundary_controls().is_empty());
        assert_eq!(renderer.held_boundary_drives(), program(1.0).initial);
        let mut direct = keyed::make(&plate, SAMPLES);
        let mut mic = microphone();
        let mut actual = vec![0.0; SAMPLES as usize];
        let mut sample = 0;
        let mut key_work = 0.0;
        for chunk in actual.chunks_mut(block) {
            renderer.block(chunk).unwrap();
            let mut last = ApertureNetworkFrame::default();
            for value in chunk {
                last = direct.step_with_boundaries(TubeDrive { upstream_pressure_pa: 100.0,
                    body_flow_m3_s: 0.0 }, 0.0, &[boundary_at(sample, 1.0)]).unwrap();
                let expected = mic.step(&[direct.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap();
                assert_eq!(value.to_bits(), expected.to_bits());
                key_work += last.boundaries.unwrap().mechanical_work_j.abs();
                sample += 1;
            }
            let CoupledAperture::Network(network) = renderer.system() else { unreachable!() };
            assert_eq!(network.aperture().state(), direct.aperture().state());
            assert_eq!(network.boundary_apertures().unwrap().last_frame(),
                direct.boundary_apertures().unwrap().last_frame());
            assert_eq!(renderer.last_mechanical_work_j(), last.aperture.mechanical_work_j);
            assert_eq!(renderer.last_boundary_work_j(), last.boundaries.unwrap().mechanical_work_j);
            assert_eq!(renderer.held_boundary_drives(), &[boundary_at(sample - 1, 1.0)]);
        }
        assert!(key_work > 0.0);
        assert_eq!(renderer.applied_boundary_controls().len(), 4);
        assert!(renderer.pending_boundary_controls().is_empty());
        assert!(actual.iter().any(|p| p.abs() > 1e-10));
        if reference.is_empty() { reference = actual; } else { assert_eq!(actual, reference); }
    }
    let mut dense_chart = plate.chart().clone();
    dense_chart.section.density *= 2.0;
    let dense = PlateApertureReduction::from_chart(dense_chart, plate.options().clone(),
        &CancelGate::new_clock_free()).unwrap();
    let base_pcm = encode_pcm16_wav(&reference, RATE, 1.0).unwrap().0;
    for changed in [pressure(single(&dense, 37, 64).with_boundary_controls(program(1.0)).unwrap(), 37),
        pressure(single(&plate, 37, 64).with_boundary_controls(program(0.0)).unwrap(), 37)] {
        assert_ne!(encode_pcm16_wav(&changed, RATE, 1.0).unwrap().0, base_pcm,
            "key material and force changes must survive the original PCM quantizer");
    }
}

#[test]
fn double_reed_and_key_sample_zero_controls_match_direct_frames_and_receiver_across_blocks() {
    let plate = key_plate();
    let mut reference_pcm = Vec::new();
    for block in [1, 37, 128] {
        let mut renderer = paired(&plate, block, 64).with_boundary_controls(program(1.0)).unwrap()
            .with_blade_forces(blade_forces()).unwrap();
        let mut direct = keyed_paired::make(&plate, SAMPLES);
        let mut mic = microphone();
        let mut actual = vec![0.0; SAMPLES as usize];
        let mut sample = 0;
        for chunk in actual.chunks_mut(block) {
            renderer.block(chunk).unwrap();
            let mut last = None;
            for value in chunk {
                let forces = if sample < 701 { [0.001, -0.0007] } else { [0.0; 2] };
                last = Some(direct.step_with_boundaries(TubeDrive { upstream_pressure_pa: 100.0,
                    body_flow_m3_s: 0.0 }, forces, &[boundary_at(sample, 1.0)]).unwrap());
                let expected = mic.step(&[direct.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap();
                assert_eq!(value.to_bits(), expected.to_bits());
                sample += 1;
            }
            assert_eq!(renderer.last_frame(), last.as_ref());
            assert_eq!(renderer.system().aperture().blade_states(), direct.aperture().blade_states());
            assert_eq!(renderer.held_boundary_drives(), &[boundary_at(sample - 1, 1.0)]);
            assert_eq!(renderer.last_boundary_work_j(), last.unwrap().coupled.boundaries.unwrap().mechanical_work_j);
        }
        assert_eq!(renderer.applied_force_controls().len(), 4);
        assert_eq!(renderer.applied_boundary_controls().len(), 4);
        assert!(renderer.pending_boundary_controls().is_empty());
        let pcm = encode_pcm16_wav(&actual, RATE, 1.0).unwrap().0;
        if reference_pcm.is_empty() { reference_pcm = pcm; } else { assert_eq!(pcm, reference_pcm); }
    }
}

#[test]
fn missing_key_program_preserves_output_and_combined_budget_is_attachment_order_independent() {
    let plate = key_plate();
    let mut one = single(&plate, 8, 64);
    let mut two = paired(&plate, 8, 64);
    for source in [&mut one as &mut dyn PressureRenderer, &mut two as &mut dyn PressureRenderer] {
        let mut output = [123.0; 8];
        assert!(source.validate_sample_count(8).is_err());
        assert!(source.block(&mut output).is_err());
        assert_eq!(output, [123.0; 8]);
        assert_eq!(source.samples_rendered(), 0);
    }
    for invalid in [BoundaryProgram { controls: vec![BoundaryControl { sample: SAMPLES,
        drive: boundary_at(0, 1.0) }], ..program(1.0) },
        BoundaryProgram { controls: vec![BoundaryControl { sample: 0, drive: boundary_at(0, 1.0) }; 2],
            ..program(1.0) }] {
        assert!(single(&plate, 8, 64).with_boundary_controls(invalid).is_err());
    }
    // One pressure assignment + five key assignments + four blade assignments.
    assert!(paired(&plate, 8, 9).with_boundary_controls(program(1.0)).unwrap()
        .with_blade_forces(blade_forces()).is_err());
    assert!(paired(&plate, 8, 9).with_blade_forces(blade_forces()).unwrap()
        .with_boundary_controls(program(1.0)).is_err());
    // A retained source plate also shares the same budget with physical force ports.
    let plate_source = || {
        let (graph, gases) = geometry();
        let reed = key(&plate, gases[0], graph.inlet_impedance_with_gases(&gases).unwrap(), SAMPLES).aperture;
        let network = ApertureNetwork::with_section_gases(reed, graph, gases.clone()).unwrap();
        let boundary = key(&plate, gases[2], network.represented_sections()[2].impedance_pa_s_m3, SAMPLES);
        let network = network.with_boundary_apertures(&gases, vec![boundary], 1).unwrap();
        AperturePerformance::new(CoupledAperture::Network(network), observation(), phrase(), config(8, 6)).unwrap()
    };
    let force = || ApertureForceProgram { footprints: vec![PlateForceFootprint::Node(6)],
        events: vec![ApertureForceEvent { sample: 0, port: 0, force_n: 0.001 }] };
    assert!(plate_source().with_boundary_controls(program(1.0)).unwrap().with_plate_forces(force()).is_err());
    assert!(plate_source().with_plate_forces(force()).unwrap().with_boundary_controls(program(1.0)).is_err());
}
