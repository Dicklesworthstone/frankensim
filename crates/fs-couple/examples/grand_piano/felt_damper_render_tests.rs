//! G3: supplied finite-gap pads reach Rayleigh pressure on the score clock.
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
        "frankensim-felt-damper-render-{}-{stamp}-{}.fspd",
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
fn finite_gap_file_changes_release_pressure_and_matches_direct_split_playback() {
    let course = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let card = "# Supplied numerical fixture, not a measured damper.\n\
        frankensim-piano-felt-dampers-v1 # header selection accepts comments\n\
        pad,69,0.12,0.14,0.00001,0.004,0,0.001,0.05,1,0,2\n\
        felt,69,60000,0.2,2.5,3.2,0.1,0.8,300000\n\
        branch,69,200000,0.002\n";
    let path = input(card);
    let mut options = Options::parse(
        &[
            "--render",
            "unused.wav",
            "--modes",
            "12",
            "--substeps",
            "8",
            "--dampers",
            &path,
        ]
        .map(str::to_owned),
    )
    .unwrap();
    let stretching = linear::string_stretching::Specification::read(
        "frankensim-piano-string-stretching-v1\nstretch,69,150000,0.2\n",
        &[course],
    )
    .unwrap();
    let build = |options: &Options| {
        prepare_instrument_with_string_material(vec![course], &modes, options, Some(&stretching))
            .unwrap()
    };
    let configured = build(&options);
    options.dampers = None;
    let mut direct = build(&options);
    direct
        .configure_damper_spec(&engine::DamperSpecification::read(card, &[course]).unwrap())
        .unwrap();
    let mut free = build(&options);
    free.configure_damper_spec(
        &engine::DamperSpecification::read(
            "frankensim-piano-felt-dampers-v1\nfree,69\n",
            &[course],
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        configured.felt_damper_resolution(),
        Some((course.unison, 2 * course.unison))
    );
    assert!(configured.damper_resolution().is_none());
    assert!(configured.bank.has_string_stretching());
    assert_eq!(configured.damper_energy_j(), 0.0);
    let score = || {
        performance::Performance::read(
            "sample,event,key,value\n0,note_on,69,1\n\
            480,note_off,69,0\n480,sustain,0,0.25\n\
            960,sustain,0,0\n1800,sustain,0,1\n",
            &[69],
            2400,
        )
        .unwrap()
    };
    let surface = [board_geometry::SurfaceSample {
        position_m: [0.; 3],
        area_m2: 0.2,
        mode_shape: modes.iter().map(|m| m.volume / 0.2).collect(),
    }];
    let stream = |piano| {
        audio::AudioStream::new_stereo(
            piano,
            score(),
            &surface,
            [[0., 0., 1.], [0.3, 0., 1.2]],
            fs_bem::helmholtz::Medium::air(),
        )
        .unwrap()
    };
    let mut configured = stream(configured);
    let mut direct = stream(direct);
    let mut free = stream(free);
    let mut whole = vec![0.; 4800];
    let mut split = whole.clone();
    let mut free_pressure = whole.clone();
    configured.render_interleaved_block(&mut whole).unwrap();
    for block in split.chunks_mut(274) {
        direct.render_interleaved_block(block).unwrap();
    }
    free.render_interleaved_block(&mut free_pressure).unwrap();
    assert_eq!(whole, split);
    assert!(whole.iter().all(|x| x.is_finite()));
    assert!(whole.iter().any(|x| x.abs() > 1e-10));
    assert!(whole.chunks_exact(2).any(|frame| frame[0] != frame[1]));
    assert_ne!(whole, free_pressure);
    assert_eq!(configured.sample_position(), 2400);
    assert_eq!(direct.sample_position(), 2400);
    let configured = configured.instrument();
    let direct = direct.instrument();
    assert_eq!(configured.bank.q, direct.bank.q);
    assert_eq!(configured.bank.v, direct.bank.v);
    assert_eq!(configured.damper_energy_j(), direct.damper_energy_j());
    assert_eq!(
        configured.accounting.damper_actuator_work_j,
        direct.accounting.damper_actuator_work_j
    );
    assert_eq!(
        configured.accounting.damper_loss_j,
        direct.accounting.damper_loss_j
    );
    assert_eq!(
        configured.accounting.damper_relaxation_loss_j,
        direct.accounting.damper_relaxation_loss_j
    );
    assert!(configured.accounting.damper_loss_j > 0.0);
    assert!(configured.accounting.damper_relaxation_loss_j > 0.0);
    assert_ne!(configured.accounting.damper_actuator_work_j, 0.0);
    assert_eq!(free.instrument().accounting.damper_loss_j, 0.0);
    assert_eq!(free.instrument().damper_energy_j(), 0.0);
    for piano in [configured, direct, free.instrument()] {
        assert!(
            (piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j())
                .abs()
                < 1e-7
        );
    }
}
