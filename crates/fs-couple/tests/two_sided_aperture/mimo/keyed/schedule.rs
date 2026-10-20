//! Original material keys, mechanical sample controls and the actual microphone.
use super::*;
use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use fs_couple::bernoulli_aperture::performance::{
    AperturePerformanceConfig,
    boundaries::{BoundaryControl, BoundaryProgram},
    force::{ApertureForceEvent, ApertureForceProgram},
};
use fs_couple::bernoulli_aperture::two_sided::observer::schedule::ScheduledTwoSidedPressure;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

fn radiation() -> DuctRadiationSpec {
    let mut radiation = observation();
    radiation.reconstruction = DuctReconstruction::Differentiating {
        mechanics_rate_hz: 48_000, output_rate_hz: 3_000,
    };
    radiation
}
fn phrase() -> GestureSchedule {
    GestureSchedule::try_new(1, vec![GestureTrack { id: "mouth".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(20.0),
        events: vec![] }]).unwrap()
}
fn config() -> AperturePerformanceConfig {
    AperturePerformanceConfig { sample_rate_hz: 3_000, samples: 64, max_block: 37,
        max_compile_work: 10_000, max_controls: 64 }
}
fn boundary(sample: u64, scale: f64) -> BoundaryDrive {
    drives(scale * match sample { 0..7 => 0.0, 7..193 => -0.01,
        193..511 => -0.02, _ => 0.0 }, if (193..511).contains(&sample) { 0.3 } else { 0.0 })[0]
}
fn program(scale: f64) -> BoundaryProgram {
    BoundaryProgram { initial: vec![boundary(0, scale)],
        controls: [511, 7, 193].into_iter().map(|sample| BoundaryControl {
            sample, drive: boundary(sample, scale),
        }).collect() }
}
fn performer(plate: &PlateApertureReduction, config: AperturePerformanceConfig)
    -> ScheduledTwoSidedPressure
{
    let source = TwoSidedPressure::with_radiation(make(plate, Kind::Modal, 1024, limits()),
        radiation(), 37).unwrap();
    ScheduledTwoSidedPressure::new(source, phrase(), config).unwrap()
}

#[test]
fn subframe_key_events_match_direct_physics_receiver_and_pcm_for_all_callback_sizes() {
    let plate = plate(900.0);
    let mut direct = make(&plate, Kind::Modal, 1024, limits());
    let mut microphone = DuctPressure::new(&SIDES.map(|side| DuctView {
        graph: direct.graph(side), gases: direct.section_gases(side),
    }), radiation()).unwrap();
    let mut expected = Vec::new();
    let mut last = None;
    for sample in 0..1024 {
        let up = [boundary(sample, 0.7)];
        let down = [boundary(sample, 1.0)];
        last = Some(direct.step_with_boundaries(drive(20.0), [&up, &down]).unwrap());
        if let Some(pressure) = microphone.push_frame(|network, node|
            direct.node_frame(SIDES[network], node).map(|f| f.net_flow_into_node_m3_s)).unwrap() {
            expected.push(pressure);
        }
    }
    assert_eq!(expected.len(), 64);
    assert!(expected.iter().any(|p| p.abs() > 1e-12));
    for block in [1, 7, 37] {
        let mut performer = performer(&plate, config())
            .with_boundary_controls(ValveSide::Upstream, program(0.7)).unwrap()
            .with_boundary_controls(ValveSide::Downstream, program(1.0)).unwrap();
        let mut actual = vec![0.0; 64];
        for chunk in actual.chunks_mut(block) { performer.block(chunk).unwrap(); }
        assert_eq!(actual, expected);
        assert_eq!(encode_pcm16_wav(&actual, 3_000, 1.0).unwrap().0,
            encode_pcm16_wav(&expected, 3_000, 1.0).unwrap().0);
        assert_eq!(performer.source().last_frame(), last.as_ref());
        same_state(performer.source().system(), &direct);
        for (side, scale) in SIDES.into_iter().zip([0.7, 1.0]) {
            assert_eq!(performer.held_boundary_drives(side), &[boundary(1023, scale)]);
            assert_eq!(performer.applied_boundary_controls(side).len(), 3);
            assert!(performer.pending_boundary_controls(side).is_empty());
            assert_eq!(performer.last_boundary_frame(side),
                direct.boundary_apertures(side).unwrap().last_frame());
        }
    }
}

#[test]
fn incomplete_program_and_failed_second_side_do_not_publish_key_controls() {
    let plate = plate(900.0);
    let mut partial = performer(&plate, config())
        .with_boundary_controls(ValveSide::Upstream, program(0.7)).unwrap();
    let mut output = [123.0];
    assert!(partial.block(&mut output).is_err());
    assert_eq!(output, [123.0]);
    assert_eq!(partial.source().mechanical_samples_rendered(), 0);
    // Preflight refusal is recoverable: the absent side can still attach cold.
    let mut complete = partial.with_boundary_controls(ValveSide::Downstream, program(1.0)).unwrap();
    complete.block(&mut output).unwrap();
    assert_eq!(complete.source().mechanical_samples_rendered(), 16);
    assert_eq!(complete.applied_boundary_controls(ValveSide::Upstream)[0].sample, 7);

    let bad = BoundaryProgram { initial: vec![boundary(0, 1.0)],
        controls: vec![BoundaryControl { sample: 7,
            drive: BoundaryDrive { opening_force_n: f64::MAX, ..boundary(7, 1.0) } }] };
    let mut failing = performer(&plate, config())
        .with_boundary_controls(ValveSide::Upstream, program(0.7)).unwrap()
        .with_boundary_controls(ValveSide::Downstream, bad).unwrap();
    output = [123.0];
    assert!(failing.block(&mut output).is_err());
    assert_eq!(output, [123.0]);
    assert_eq!(failing.samples_rendered(), 0);
    assert_eq!(failing.source().mechanical_samples_rendered(), 7);
    for side in SIDES {
        assert!(failing.applied_boundary_controls(side).is_empty());
        assert_eq!(failing.held_boundary_drives(side), &[boundary(0, 1.0)]);
        assert_eq!(failing.source().system().boundary_apertures(side).unwrap()
            .aperture(4).unwrap().accepted_steps(), 7);
    }
}

fn force_program() -> ApertureForceProgram {
    ApertureForceProgram { footprints: vec![PlateForceFootprint::Node(6)],
        events: vec![ApertureForceEvent { sample: 0, port: 0, force_n: 0.001 },
            ApertureForceEvent { sample: 23, port: 0, force_n: 0.0 }] }
}
fn plate_performer(plate: &PlateApertureReduction, config: AperturePerformanceConfig)
    -> ScheduledTwoSidedPressure
{
    let ducts = ducts();
    let (gas, z) = combined_inlet_state(&ducts).unwrap();
    let valve = key(plate, gas, z, 1024).aperture;
    let system = TwoSidedAperture::new(valve.into(), ducts, limits(), &CancelGate::new_clock_free()).unwrap();
    let system = add_keys(attach(system, 0.04), plate);
    let source = TwoSidedPressure::with_radiation(system, radiation(), 37).unwrap();
    ScheduledTwoSidedPressure::new(source, phrase(), config).unwrap()
}

#[test]
fn both_key_sides_and_source_force_share_storage_and_work_in_either_attachment_order() {
    let plate = plate(900.0);
    for config in [AperturePerformanceConfig { max_controls: 10, ..config() },
        AperturePerformanceConfig { max_compile_work: 84, ..config() }] {
        assert!(plate_performer(&plate, config)
            .with_boundary_controls(ValveSide::Upstream, program(0.7)).unwrap()
            .with_boundary_controls(ValveSide::Downstream, program(1.0)).unwrap()
            .with_plate_forces(0, force_program()).is_err());
        assert!(plate_performer(&plate, config).with_plate_forces(0, force_program()).unwrap()
            .with_boundary_controls(ValveSide::Upstream, program(0.7)).unwrap()
            .with_boundary_controls(ValveSide::Downstream, program(1.0)).is_err());
    }
}
