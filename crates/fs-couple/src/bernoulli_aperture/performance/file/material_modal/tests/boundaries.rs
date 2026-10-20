//! Original material key -> both coupled media -> actual multirate pressure.
use super::*;
use crate::bernoulli_aperture::material_boundary::boundary_closing_force_path;
use crate::bernoulli_aperture::two_sided::observer::schedule::boundaries::source_pressure_path;
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};

fn keyed(force: bool) -> String {
    let original = input();
    let header = original.split_once("\nschedule\n").unwrap().0
        .replace("network 3 2 1048576\nduct_node inlet\nduct_node junction\nduct_node reflection 0\n",
            "network 4 3 1048576\nduct_node inlet\nduct_node junction\nduct_node reflection 0\nduct_node driven\n")
        .replace("\npropagation ", "\nduct_section 1 3 0.017 0.003 0.002\npropagation ")
        .replace("\ncompile_limits ", concat!(
            "\ndownstream_boundaries\nboundaries 1\nboundary 3 instrument key-1 0 0\n",
            "boundary_plate key.matpack 293.15 1000 0.025 0.01 0.0003 6 2 0.0002 0.0002 0\n",
            "boundary_mode 0 600 0 0.25 0.1\nboundary_contact 10000000000 2 0.5\n",
            "boundary_slit 0.0003\nboundary_loss elastic\ncontrol_source instrument reed\ncompile_limits "));
    let schedule = GestureSchedule::try_new_v2(
        96000,
        vec![
            GestureTrack {
                id: "breath".into(),
                target: GestureTarget::Parameter(source_pressure_path("instrument", "reed")),
                initial: GestureValue::PressurePa(5.0),
                events: vec![GestureEvent {
                    time_s: 0.004,
                    transition_s: 0.001,
                    value: GestureValue::PressurePa(0.0),
                }],
            },
            GestureTrack {
                id: "key".into(),
                target: GestureTarget::Parameter(boundary_closing_force_path(
                    "instrument",
                    "key-1",
                )),
                initial: GestureValue::ForceN(0.0),
                events: if force {
                    vec![
                        GestureEvent {
                            time_s: 1.0 / 96000.0,
                            transition_s: 0.0,
                            value: GestureValue::ForceN(0.0005),
                        },
                        GestureEvent {
                            time_s: 0.003,
                            transition_s: 0.001,
                            value: GestureValue::ForceN(0.0),
                        },
                    ]
                } else {
                    vec![]
                },
            },
        ],
    )
    .unwrap();
    format!(
        "{header}\nschedule\n{}",
        String::from_utf8(schedule.to_canonical_bytes()).unwrap()
    )
}
fn load(text: &str, block: usize, key_density: f64) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            Ok(pack(
                if path == "key.matpack" {
                    key_density
                } else {
                    900.0
                },
                false,
            ))
        },
    )
    .unwrap()
}
fn key(source: &ScheduledTwoSidedPressure) -> &crate::bernoulli_aperture::dynamic::DynamicAperture {
    source
        .source()
        .system()
        .boundary_apertures(ValveSide::Downstream)
        .unwrap()
        .aperture(3)
        .unwrap()
}

#[test]
fn physical_material_keys_and_subframe_force_events_reach_the_original_modal_microphone() {
    let text = keyed(true);
    let loaded = load(&text, 37, 900.0);
    assert_eq!(loaded.info().modes, 6);
    assert_eq!(loaded.boundary_records().len(), 1);
    assert_eq!(loaded.boundary_records()[0].0, ValveSide::Downstream);
    assert_eq!(loaded.boundary_records()[0].1.node(), 3);
    assert_eq!(
        loaded.boundary_records()[0].1.closing_force(),
        &boundary_closing_force_path("instrument", "key-1")
    );
    let owner = key(loaded.renderer());
    assert_eq!(owner.spec().max_steps, 512 * 16);
    assert_eq!(
        owner.viscous_slit().unwrap().gas(),
        loaded
            .renderer()
            .source()
            .system()
            .boundary_apertures(ValveSide::Downstream)
            .unwrap()
            .gas(3)
    );
    assert!(
        loaded
            .renderer()
            .pending_boundary_controls(ValveSide::Downstream)
            .iter()
            .any(|event| event.sample % 16 == 8),
        "force changes must retain sub-audio-frame timing"
    );
    let (pressure, rendered) = render(loaded);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(key(&rendered).accepted_steps(), 512 * 16);
    assert!(
        rendered
            .pending_boundary_controls(ValveSide::Downstream)
            .is_empty()
    );
    assert_eq!(
        rendered.held_boundary_drives(ValveSide::Downstream)[0].opening_force_n,
        0.0
    );
    assert!(
        rendered
            .last_boundary_frame(ValveSide::Downstream)
            .is_some()
    );
    assert_eq!(pressure, render(load(&text, 1, 900.0)).0);
    assert_eq!(pressure, render(load(&text, 128, 900.0)).0);
    let (unchanged_force, _) = render(load(&keyed(false), 37, 900.0));
    let (material, changed) = render(load(&text, 37, 1800.0));
    assert_ne!(pressure, unchanged_force);
    assert_ne!(pressure, material);
    assert_ne!(key(&rendered).state(), key(&changed).state());
    assert_ne!(
        encode_pcm16_wav(&pressure, 48000, 0.005).unwrap().0,
        encode_pcm16_wav(&material, 48000, 0.005).unwrap().0
    );
}

#[test]
fn keys_share_original_material_cache_and_single_blade_owner() {
    let text = keyed(true).replace("boundary_plate key.matpack", "boundary_plate blade.matpack");
    let mut loads = 0;
    let source = MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        },
    )
    .unwrap();
    assert_eq!(loads, 1);
    assert_eq!(
        source.boundary_records()[0].1.pack(),
        source.material_records()[0].pack()
    );
    let first = text.find("\nblade ").unwrap();
    let second = first + 1 + text[first + 1..].find("\nblade ").unwrap();
    let tail = second + text[second..].find("\ndownstream_boundaries\n").unwrap();
    let single = format!("{}{}", &text[..second], &text[tail..]).replace("blades 2", "blades 1");
    let mut source = load(&single, 37, 900.0).into_renderer();
    let LoadedAperture::Modal(valve) = source.source().system().valve() else {
        panic!("original modal valve")
    };
    assert_eq!(valve.material_plates().len(), 1);
    assert_eq!(valve.modal_states().len(), 3);
    source.block(&mut [0.0; 37]).unwrap();
    assert_eq!(key(&source).accepted_steps(), 37 * 16);
}

#[test]
fn incomplete_keys_and_unbound_typed_consumers_refuse_before_playback() {
    let text = keyed(true);
    let gate = CancelGate::new_clock_free();
    for altered in [
        text.replace("boundary 3 instrument", "boundary 2 instrument"),
        text.replace(
            "control_source instrument reed",
            "control_source instrument unbound",
        ),
        text.replace(
            "boundary 3 instrument key-1",
            "boundary 3 instrument unknown-key",
        ),
        text.replace(
            "boundary_plate key.matpack",
            "boundary_plate ../key.matpack",
        ),
        text.replace("compile_limits 100000 8192", "compile_limits 1 8192"),
    ] {
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                altered.as_bytes(),
                37,
                &gate,
                None,
                |_| Ok(pack(900.0, false))
            )
            .is_err()
        );
    }
    let declared = "downstream_boundaries\nboundaries 1";
    let start = text.find(declared).unwrap();
    let end = start + text[start..].find("control_source").unwrap();
    let missing = format!("{}{}", &text[..start], &text[end..]);
    assert!(
        MaterialModalValvePerformance::from_bytes_with_sources(
            missing.as_bytes(),
            37,
            &gate,
            None,
            |_| Ok(pack(900.0, false))
        )
        .is_err()
    );
}
