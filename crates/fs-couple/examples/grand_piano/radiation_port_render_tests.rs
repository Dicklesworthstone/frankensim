//! G3: a complete 40-coordinate piano reaches the real BEM/feedback/stereo path.
use super::*;
use fs_math::c64::C64;

fn many_mode_scene(radiation_ports: Option<usize>) -> Scene {
    use board_geometry::motion::MotionSurface;
    // Forty independently sprung, mass-normalized coordinates intentionally
    // share two interface fields. These are authored kinematics, not FE
    // eigenpairs or a high-frequency soundboard convergence claim.
    let points = vec![[0., 0., 0.], [0.1, 0., 0.], [0.1, 0.1, 0.], [0., 0.1, 0.]];
    let mesh = fs_plate::ShellMesh::new(points.clone(), vec![[0, 1, 2], [0, 2, 3]]).unwrap();
    let mut modes = Vec::new();
    let mut shapes = Vec::new();
    for j in 0..40 {
        let a = 0.05 * (1 + j % 3) as f64;
        let b = if j % 2 == 0 { 0.15 } else { -0.15 };
        shapes.push(
            points
                .iter()
                .map(|p| [0., 0., a + b * (p[0] / 0.05 - 1.), 0., 0., 0.])
                .collect(),
        );
        let mut bridge = [0.; 88];
        bridge[48] = a; // the same field at the center bridge point
        modes.push(linear::BoardMode {
            frequency_hz: 90. + 4. * j as f64,
            damping_ratio: 0.01,
            bridge,
            volume: 0.01 * a,
        });
    }
    let motion = MotionSurface::new(mesh, shapes).unwrap();
    let (_, courses, obj, spec) = tests::small_source_inputs();
    let options = playback::Options {
        radiation_ports,
        ..playback::Options::default()
    };
    let controls = playback::Controls::from_texts(&courses, None, None, None).unwrap();
    let piano = controls
        .instrument_with_motion(courses, &modes, Some(&motion), &options)
        .unwrap();
    let boundary = Boundary::from_obj(&obj, &spec, &motion)
        .unwrap()
        .loaded(&piano.bank)
        .unwrap();
    let board = board_geometry::PreparedBoard {
        surface: vec![],
        modes,
        motion: Some(motion),
        provenance:
            "estimated,authored 40-coordinate rank-two interface fixture; no FE eigenvalue claim"
                .into(),
        area_m2: 0.01,
        mass_kg: 40.,
        frequency_intervals_hz: vec![],
        physical_damping: None,
        reduction: None,
        free_dofs: 40,
    };
    Scene {
        piano,
        board,
        boundary,
        spec,
        radiation_ports,
    }
}

#[test]
fn radiation_rank_control_is_explicit_and_refused_admission_keeps_all_mechanics() {
    let options = playback::Options::parse(&["--radiation-ports".into(), "2".into()]).unwrap();
    assert_eq!(options.radiation_ports, Some(2));
    assert!(
        playback::Options::harmonic(&["--radiation-ports".into(), "2".into()])
            .unwrap_err()
            .contains("render-loaded")
    );
    for count in ["0", "33", "NaN", "-1"] {
        assert!(playback::Options::parse(&["--radiation-ports".into(), count.into()]).is_err());
    }
    assert!(playback::Options::parse(&["--radiation-ports".into()]).is_err());
    assert!(
        playback::Options::parse(&[
            "--radiation-ports".into(),
            "2".into(),
            "--radiation-ports".into(),
            "3".into()
        ])
        .is_err()
    );
    let args = [
        "render",
        "missing.fsb",
        "missing.csv",
        "missing.obj",
        "missing.fspe",
        "unused.wav",
        "0.05",
        "--radiation-ports",
        "2",
    ]
    .map(str::to_owned);
    assert!(run(&args).unwrap_err().contains("render-loaded"));

    // The parsed control survives the real geometry/physical-control front door.
    let (board, courses, obj, spec) = tests::small_source_inputs();
    let controls = playback::Controls::from_texts(&courses, None, None, None).unwrap();
    let prepared = prepare_controlled(&board, courses, &obj, spec, &options, controls).unwrap();
    assert_eq!(prepared.radiation_ports, Some(2));

    let mut scene = many_mode_scene(None);
    let q = scene.piano.bank.q.clone();
    let v = scene.piano.bank.v.clone();
    let energy = scene.piano.energy_j();
    let modes = scene.piano.bank.board_count;
    assert_eq!(modes, 40);
    assert!(
        bake(&mut scene, true)
            .err()
            .unwrap()
            .contains("<=32 complete board modes")
    );
    assert!(!scene.piano.has_radiation());
    scene.radiation_ports = Some(2);
    assert!(
        bake(&mut scene, false)
            .err()
            .unwrap()
            .contains("render-loaded")
    );
    scene.spec.frequencies = 17;
    assert!(bake(&mut scene, true).is_err());
    assert!(!scene.piano.has_radiation());
    assert_eq!(scene.piano.bank.board_count, modes);
    assert_eq!(scene.piano.bank.q, q);
    assert_eq!(scene.piano.bank.v, v);
    assert_eq!(scene.piano.energy_j(), energy);
}

#[test]
fn forty_board_coordinates_reach_projected_bem_feedback_and_complete_stereo_receivers() {
    let mut scene = many_mode_scene(Some(2));
    let mut manual = many_mode_scene(Some(2));
    let mut bare = many_mode_scene(None);
    let (baked, samples, report) = bake(&mut scene, true).unwrap();
    assert!(report.contains("40 complete board coordinates retained, 2 acoustic ports"));
    assert!(report.contains("All 40 board inputs remain in every receiver transfer"));
    assert!(scene.piano.has_radiation());
    assert_eq!(scene.piano.bank.board_count, 40);
    assert_eq!(
        scene.piano.board_trace_len(),
        40 * playback::Options::default().substeps
    );

    // Prepare the same physical owner directly, independently of bake's
    // selection/attachment branch, and retain the complete receiver matrix.
    let (projected, direct_samples) =
        radiation_fit::prepare_projected(&manual.boundary, &manual.spec, manual.piano.bank.rate, 2)
            .unwrap();
    assert_eq!(projected.basis.board_ports(), 40);
    assert_eq!(projected.basis.ports(), 2);
    assert!(projected.surface_rms_error < 1e-10);
    assert!(projected.fit.peak_error <= 0.05);
    assert!(projected.fit.resistance_peak_error <= 0.10);
    manual
        .piano
        .configure_projected_radiation(&projected.fit.model, &projected.basis)
        .unwrap();
    assert_eq!(samples.values, direct_samples.values);
    assert_eq!(samples.values.len(), 2);
    assert!(samples.values.iter().all(|rows| rows.len() == 40));
    let f = samples.omega.len() / 2;
    let w = samples.omega[f];
    let direct = exterior_loading::sample(&scene.boundary, &scene.spec, w).unwrap();
    for (channel, rows) in samples.values.iter().enumerate() {
        for (input, row) in rows.iter().enumerate() {
            assert_eq!(row.len(), samples.omega.len());
            assert!(row.iter().any(|h| h.abs() > 0.));
            let recovered = row[f] * C64::new(0., -w);
            let expected = direct.receiver_transfer[channel][input];
            assert!((recovered - expected).abs() < 1e-10 * (1. + expected.abs()));
        }
    }

    let score = || {
        performance::Performance::read(
            "sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",
            &[69],
            2400,
        )
        .unwrap()
    };
    let audio = exterior_audio::render(&mut scene.piano, score(), 2400, &baked, 4.).unwrap();
    let mut events = score();
    let mut unreacted = score();
    for n in 0..2400 {
        events.dispatch(n, &mut manual.piano).unwrap();
        manual.piano.step().unwrap();
        unreacted.dispatch(n, &mut bare.piano).unwrap();
        bare.piano.step().unwrap();
    }
    assert_eq!(scene.piano.bank.q, manual.piano.bank.q);
    assert_eq!(scene.piano.bank.v, manual.piano.bank.v);
    assert_eq!(
        scene.piano.radiation_energy_j(),
        manual.piano.radiation_energy_j()
    );
    assert_eq!(
        scene.piano.accounting.radiation_loss_j,
        manual.piano.accounting.radiation_loss_j
    );
    assert_ne!(scene.piano.bank.q, bare.piano.bank.q);
    assert!(
        scene.piano.bank.q[scene.piano.bank.modes.len() + 32..]
            .iter()
            .any(|q| *q != 0.)
    );
    assert!(scene.piano.accounting.felt_loss_j > 0.);
    assert!(scene.piano.accounting.radiation_loss_j > 0.);
    let defect = scene.piano.accounting.input_work_j
        - scene.piano.energy_j()
        - scene.piano.accounting.dissipated_j();
    assert!(defect.abs() < 1e-7, "complete piano work defect {defect:e}");
    assert!(audio.peak_pa > 1e-14);
    assert_eq!(u16::from_le_bytes([audio.wav[22], audio.wav[23]]), 2);
    let data = audio.wav.windows(4).position(|w| w == b"data").unwrap() + 8;
    assert_eq!(audio.wav[data..].len(), 2400 * 4);
    for frame in audio.wav[data..].chunks_exact(4) {
        assert_eq!(&frame[..2], &frame[2..]);
    }
}
