//! G3: admitted axial waves reach physical pressure on the existing score clock.
use super::*;
use std::{
    io::Write,
    sync::atomic::{AtomicU64, Ordering},
};

fn input(text: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "frankensim-longitudinal-render-{}-{stamp}-{}.fspx",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    file.write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}

#[test]
fn finite_axial_waves_change_pressure_and_preserve_midi_csv_and_block_timing() {
    let course = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let waves = input(&format!(
        "{}\nwaves,69,150000,0.2,60000,16,0.02,0.05\n",
        linear::string_stretching::HEADER
    ));
    let averaged = input(&format!(
        "{}\nstretch,69,150000,0.2\n",
        linear::string_stretching::HEADER
    ));
    let faces = input("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n");
    let mut options = Options::parse(
        &[
            "--render",
            "unused.wav",
            "--modes",
            "12",
            "--substeps",
            "8",
            "--dampers",
            "estimated",
            "--hammer-footprints",
            &faces,
            "--string-stretching",
            &waves,
        ]
        .map(str::to_owned),
    )
    .unwrap();
    // The waves row's axial coupling comes from the board; these demonstration modes carry no rotations, so the
    // test supplies its fixture row (a fixed fraction of the vertical bridge row) explicitly.
    let material = load_string_stretching(&[course], &options)
        .unwrap()
        .unwrap()
        .with_fixture_axial_rows(&modes);
    let first = prepare_instrument_with_string_material(vec![course], &modes, &options, Some(&material)).unwrap();
    let second = prepare_instrument_with_string_material(vec![course], &modes, &options, Some(&material)).unwrap();
    options.string_stretching = Some(averaged);
    let reference = prepare_instrument(vec![course], &modes, &options).unwrap();
    assert!(first.bank.longitudinal_mode_count() > 0);
    assert_eq!(reference.bank.longitudinal_mode_count(), 0);
    assert!(first.damper_resolution().is_some());
    assert_eq!(first.hammer_contact_count(), 2 * course.unison);
    assert_eq!(
        first.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>(),
        reference
            .bank
            .modes
            .iter()
            .map(|m| m.omega)
            .collect::<Vec<_>>()
    );

    let bytes = b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk\0\0\0\x10\
        \0\xb0\x40\x7f\0\x90\x45\x7f\x04\x90\x45\0\x08\xff\x2f\0";
    let imported = midi::read(
        bytes,
        &[69],
        48_000,
        1500,
        midi::Mapping {
            maximum_velocity_m_s: 1.0,
            ..midi::Mapping::default()
        },
    )
    .unwrap();
    let imported = performance::Performance::from_events(imported.events, &[69], 1500).unwrap();
    let csv = || {
        performance::Performance::read(
            "sample,event,key,value\n0,sustain,0,1\n\
        0,note_on,69,1\n200,note_off,69,0\n600,sustain,0,0\n",
            &[69],
            1500,
        )
        .unwrap()
    };

    // The existing authored Rayleigh patch: actual Pa, without a claim that
    // these demonstration board coordinates are measured FE eigenpairs.
    let surface = [board_geometry::SurfaceSample {
        position_m: [0.; 3],
        area_m2: 0.2,
        mode_shape: modes.iter().map(|m| m.volume / 0.2).collect(),
    }];
    let make_stream = |piano, score| {
        audio::AudioStream::new_stereo(
            piano,
            score,
            &surface,
            [[0., 0., 1.], [0.3, 0., 1.2]],
            fs_bem::helmholtz::Medium::air(),
        )
        .unwrap()
    };
    let mut first = make_stream(first, imported);
    let mut second = make_stream(second, csv());
    let mut reference = make_stream(reference, csv());
    let mut whole = vec![0.; 3000];
    let mut split = whole.clone();
    let mut averaged_pressure = whole.clone();
    first.render_interleaved_block(&mut whole).unwrap();
    for block in split.chunks_mut(274) {
        second.render_interleaved_block(block).unwrap();
    }
    reference
        .render_interleaved_block(&mut averaged_pressure)
        .unwrap();
    assert_eq!(whole, split);
    assert!(whole.iter().all(|x| x.is_finite()));
    assert!(whole.iter().any(|x| x.abs() > 1e-10));
    assert!(whole.chunks_exact(2).any(|f| f[0] != f[1]));
    assert_ne!(whole, averaged_pressure);
    assert_eq!(first.sample_position(), 1500);
    assert_eq!(second.sample_position(), 1500);
    let first = first.instrument();
    let second = second.instrument();
    assert_eq!(first.bank.q, second.bank.q);
    assert_eq!(first.bank.v, second.bank.v);
    let mut moving_axial_coordinate = false;
    for string in 0..first.bank.strings.len() {
        let a = first.bank.string_longitudinal_state(string);
        assert_eq!(a, second.bank.string_longitudinal_state(string));
        if let Some((q, v)) = a {
            moving_axial_coordinate |= q.iter().chain(v).any(|x| *x != 0.);
        }
    }
    assert!(moving_axial_coordinate);
    assert!(first.bank.longitudinal_energy_j() > 0.);
    assert!(first.accounting.felt_loss_j > 0. && first.accounting.damper_loss_j > 0.);
}
