//! Original multimode force ports reach the same physics, work and microphone.
use super::*;
use crate::bernoulli_aperture::network::boundaries::BoundaryDrive;
use crate::bernoulli_aperture::performance::force::ApertureForceEvent;
use crate::bernoulli_aperture::two_sided::TwoSidedDrive;
use crate::pcm_wav::baffled::openings::duct::{DuctPressure, DuctView};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const FORCED: &str = include_str!("../../../../../../examples/material_modal_valve/paired-forced.performance");
const SIDES: [ValveSide; 2] = [ValveSide::Upstream, ValveSide::Downstream];

fn source() -> String {
    FORCED.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
}
fn without_forces(text: &str) -> String {
    let (before, after) = text.split_once("source_forces ").unwrap();
    let (_, after) = after.split_once("compile_limits ").unwrap();
    format!("{before}compile_limits {after}")
}
fn event(sample: u64, force_n: f64) -> ApertureForceEvent {
    ApertureForceEvent { sample, port: 0, force_n }
}

#[test]
fn g1_subframe_source_forces_match_direct_original_modes_work_and_receiver() {
    let phrase = GestureSchedule::try_new(1, vec![GestureTrack { id: "mouth".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(5.0),
        events: vec![] }]).unwrap();
    let text = format!("{}\nschedule\n{}", source().split_once("\nschedule\n").unwrap().0,
        String::from_utf8(phrase.to_canonical_bytes()).unwrap());
    let native = build(&text, 37);
    let original = valve(&native);
    let modes = original.spec().modes.clone();
    assert_eq!(original.material_plates().len(), 2);
    assert_eq!(modes.len(), 6);
    for (blade, plate) in original.material_plates().iter().enumerate() {
        for (local, mode) in modes[plate.mode_range()].iter().enumerate() {
            let projected: f64 = plate.options().reduction.unit_force_weights.iter()
                .zip(plate.displacement_shapes()).map(|(weight, shape)| weight * shape[local]).sum();
            assert!((projected - mode.force_projection[blade]).abs() < 1e-13);
        }
    }
    let mut direct = super::keyed::direct_system(&native);
    let mut receiver = DuctPressure::new(&SIDES.map(|side| DuctView {
        graph: direct.graph(side), gases: direct.section_gases(side),
    }), native.renderer().source().radiation().clone()).unwrap();
    let mut expected = Vec::new();
    let mut mechanical_work = 0.0;
    for sample in 0..4096 {
        let forces = [if (97..769).contains(&sample) { -0.0001 }
            else if (1025..1537).contains(&sample) { 0.0001 } else { 0.0 },
            if (225..2049).contains(&sample) { -0.00015 } else { 0.0 }];
        let up = [BoundaryDrive { node: 3,
            exterior_pressure_pa: if (193..769).contains(&sample) { 0.1 } else { 0.0 },
            opening_force_n: if (193..769).contains(&sample) { -0.002 } else { -0.0 } }];
        let down = [BoundaryDrive { node: 3,
            exterior_pressure_pa: if (321..1505).contains(&sample) { 0.3 } else { 0.0 },
            opening_force_n: if (321..1505).contains(&sample) { -0.003 } else { -0.0 } }];
        let LoadedAperture::Modal(reed) = direct.valve() else { panic!("original modal source"); };
        let before = reed.modal_states().to_vec();
        let frame = direct.step_with_boundaries(TwoSidedDrive { source_pressure_pa: 5.0,
            blade_forces_n: forces }, [&up, &down]).unwrap();
        let LoadedAperture::Modal(reed) = direct.valve() else { panic!("original modal source"); };
        let virtual_work: f64 = modes.iter().zip(&before).zip(reed.modal_states())
            .map(|((mode, before), after)| {
                let force = mode.force_projection[0] * forces[0] + mode.force_projection[1] * forces[1];
                force * (after.displacement_m - before.displacement_m)
            }).sum();
        assert!((virtual_work - frame.aperture.mechanical_work_j).abs()
            <= 1e-20 + 1e-10 * virtual_work.abs());
        assert!(frame.balance_residual_j().abs() < 1e-12);
        mechanical_work += frame.aperture.mechanical_work_j.abs();
        if let Some(pressure) = receiver.push_frame(|network, node|
            direct.node_frame(SIDES[network], node).map(|frame| frame.net_flow_into_node_m3_s)).unwrap() {
            expected.push(pressure);
        }
    }
    assert!(mechanical_work > 0.0);
    let (actual, rendered) = render(native);
    assert_eq!(actual, expected);
    for block in [1, 128] { assert_eq!(actual, render(build(&text, block)).0); }
    assert_eq!(rendered.held_generalized_forces_n(), [0.0; 2]);
    for (blade, count) in [(0, 4), (1, 2)] {
        assert_eq!(rendered.applied_plate_force_controls(blade).len(), count);
        assert!(rendered.pending_plate_force_controls(blade).is_empty());
    }
    let LoadedAperture::Modal(expected_reed) = direct.valve() else { panic!("original modal source"); };
    let LoadedAperture::Modal(actual_reed) = rendered.source().system().valve() else { panic!("original modal source"); };
    assert_eq!(actual_reed.modal_states(), expected_reed.modal_states());
}

#[test]
fn g3_material_force_site_and_force_release_change_the_actual_native_pressure() {
    let text = source();
    let (pressure, _) = render(build(&text, 37));
    let pcm = encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0;
    for changed in [
        without_forces(&text),
        text.replace("source_force 769 0 0", "source_force 769 0 -0.0001"),
        text.replace("section blade.matpack", "section altered.matpack"),
        text.replace("node 0.025000000000000005 0.005 0 1", "node 0.025000000000000005 0.005 0 0")
            .replace("node 0.020833333333333332 0.005 0 0", "node 0.020833333333333332 0.005 0 1"),
    ] {
        let (other, _) = render(build(&changed, 37));
        assert_ne!(pressure, other);
        assert_ne!(pcm, encode_pcm16_wav(&other, 48_000, 0.005).unwrap().0);
    }
    // Same-time updates do not invent a physical interval or reorder a port.
    let equivalent = text.replace("source_forces 6", "source_forces 7")
        .replace("source_force 97 0 -0.0001", "source_force 97 0 0.0003\nsource_force 97 0 -0.0001");
    assert_eq!(pressure, render(build(&equivalent, 37)).0);
    let ordinary = without_forces(&text);
    assert_eq!(render(build(&ordinary, 37)).0,
        render(build(&ordinary.replace("compile_limits", "source_forces none\ncompile_limits"), 37)).0);
}

#[test]
fn g0_source_force_program_preserves_material_walls_supply_and_mutual_radiation() {
    let text = include_str!("../../../../../../examples/material_modal_valve/paired-mutual.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
        .replace("network 2 1 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_section 1 0 0.12 0.007 0.002",
            "network 3 2 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_node junction\nduct_section 1 2 0.06 0.007 0.002\nduct_section 2 0 0.06 0.007 0.002")
        .replace("valve 0.0004", "walls 1 20 1000\nwall upstream 2 isotropic blade.matpack 293.15 100 0.0002 stress-free elastic\nvalve 0.0004")
        .replace("supply none", "supply slit 0.00002 0.01 0.001 16")
        .replace("compile_limits", "source_forces 2\nsource_force 97 0 -0.0001\nsource_force 1025 0 0\ncompile_limits");
    let load = |text: &str, block| MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(),
        block, &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)),
        |_| Ok(include_bytes!("../../../../../../examples/material_modal_valve/mutual.radiation").to_vec())).unwrap();
    let loaded = load(&text, 37);
    assert_eq!(loaded.wall_records().len(), 1);
    assert_eq!(loaded.mimo_radiation_records().len(), 1);
    let (pressure, rendered) = render(loaded);
    assert_eq!(pressure, render(load(&text, 1)).0);
    assert_ne!(pressure, render(load(&without_forces(&text), 37)).0);
    let system = rendered.source().system();
    assert!(system.supply_channel().is_some());
    assert_eq!(system.material_walls(ValveSide::Upstream).len(), 1);
    assert!(system.mimo_pole_states(ValveSide::Downstream).unwrap().iter()
        .any(|state| state.re != 0.0 || state.im != 0.0));
    assert!(rendered.source().last_frame().unwrap().balance_residual_j().abs() < 1e-12);
}

#[test]
fn g0_force_addresses_conflicts_and_combined_control_budgets_refuse() {
    let text = source();
    let load = |text: &str| MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)));
    for invalid in [
        text.replace("source_force 97 0", "source_force 97 2"),
        text.replace("source_force 97 0", "source_force 4096 0"),
        text.replace("source_force 97 0 -0.0001", "source_force 97 0 NaN"),
        text.replace("source_forces 6", "source_forces 262145"),
        text.replace("source_forces 6", "source_forces 7"),
    ] { assert!(load(&invalid).is_err()); }
    let ordinary = without_forces(&text);
    let baseline = load(&ordinary).unwrap();
    let r = baseline.renderer();
    let used = r.pending_controls().len() + SIDES.iter().map(|side|
        r.pending_boundary_controls(*side).len() + r.held_boundary_drives(*side).len()).sum::<usize>();
    let limit = format!("compile_limits 100000 {used}");
    assert!(load(&ordinary.replace("compile_limits 100000 8192", &limit)).is_ok());
    assert!(load(&text.replace("compile_limits 100000 8192", &limit)).is_err());
    assert!(load(&text.replace("compile_limits 100000 8192",
        &format!("compile_limits 100000 {}", used + 6))).is_ok());
    assert!(build(&ordinary, 37).into_renderer().with_material_plate_force_controls(2, vec![event(0, 0.0)]).is_err());
    let wrong_port = ApertureForceEvent { port: 1, ..event(0, 0.0001) };
    assert!(build(&ordinary, 37).into_renderer().with_material_plate_force_controls(0, vec![wrong_port]).is_err());
    let bound = build(&ordinary, 37).into_renderer()
        .with_material_plate_force_controls(0, vec![event(0, -0.0001)]).unwrap();
    assert!(bound.with_material_plate_force_controls(0, vec![event(0, 0.0)]).is_err());
    let mut held = build(&ordinary, 37).into_renderer();
    held.set_blade_forces([0.0001, 0.0]).unwrap();
    assert!(held.with_material_plate_force_controls(0, vec![event(0, 0.0)]).is_err());
    let mut advanced = build(&ordinary, 37).into_renderer();
    advanced.block(&mut [0.0]).unwrap();
    assert!(advanced.with_material_plate_force_controls(0, vec![event(0, 0.0)]).is_err());
}
