//! G3: finite longitudinal dynamics reach the existing BEM stereo consumer.
use super::*;
use crate::{exterior_audio::ExteriorStream, performance::Performance};

fn material(waves: bool) -> String {
    format!(
        "{}\n{}\n",
        linear::string_stretching::HEADER,
        if waves {
            "waves,69,150000,0.2,60000,16,0.02,0.05"
        } else {
            "stretch,69,150000,0.2"
        }
    )
}

#[test]
fn axial_state_reaches_one_way_and_loaded_stereo_with_one_contact_clock() {
    for loaded in [false, true] {
        let options = Options {
            modes: 12,
            substeps: 8,
            ..Options::default()
        };
        let (mut board, courses, obj, _) = crate::tests::small_source_inputs();
        // The string-side bridge runs diagonally from node 0 to the
        // attachment at node 4, with positive eccentricity 2 mm.
        board.push_str("stiffener,1000000,400000,0.00001,1e-11,2e-11,0.002,450,0,4\n");
        let make = |selection: Option<&str>| {
            let spec = crate::tests::small_source_inputs().3;
            let mut controls = Controls::from_texts(
                &courses,
                None,
                Some("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n"),
                Some("estimated"),
            )
            .unwrap();
            if let Some(text) = selection {
                controls = controls.with_string_stretching(text, &courses).unwrap();
            }
            crate::prepare_controlled(&board, courses.clone(), &obj, spec, &options, controls).unwrap()
        };
        let waves = material(true);
        let mut scene = make(Some(&waves));
        let mut manual = make(None);
        let mut reference = make(Some(&material(false)));
        // Attach through the direct physical API independently of the inline
        // Controls path used by render preparation.
        let mut direct = linear::string_stretching::Specification::read(&waves, &courses).unwrap();
        let motion = manual.board.motion.as_ref().unwrap();
        let unit = 1.0 / 2.0_f64.sqrt();
        let arm = 0.002 + 0.5 * (12.0_f64 * 1e-11 / 1e-5).sqrt();
        let analytic = motion.project_at(0, [0., 0., 1.], [0., 0., arm], [-unit, unit, 0.]).unwrap().0;
        let attachment = crate::board_geometry::BoardGeometry::read(&board).unwrap().bridge_attachments(&[69]).unwrap()[0];
        let axial = attachment.axial_row(motion).unwrap();
        assert!(analytic.iter().zip(&axial).all(|(a,b)| (a-b).abs() < 1e-12 * (1.0 + b.abs())));
        assert!(axial.iter().zip(&manual.board.modes)
            .any(|(&a, mode)| (a - 0.06 * mode.bridge[48]).abs() > 1e-8));
        direct.set_bridge_axial_rows(std::collections::BTreeMap::from([(69, axial)])).unwrap();
        manual.piano.configure_string_stretching(&direct).unwrap();
        let count = scene.piano.bank.longitudinal_mode_count();
        assert!(count > 0);
        assert_eq!(reference.piano.bank.longitudinal_mode_count(), 0);
        assert_eq!(scene.piano.hammer_contact_count(), 6);
        assert!(scene.piano.damper_resolution().is_some());
        assert!(
            options
                .report(&scene.piano)
                .contains(&format!("retained longitudinal coordinates: {count}"))
        );
        let (baked, _, _) = crate::bake(&mut scene, loaded).unwrap();
        crate::bake(&mut manual, loaded).unwrap();
        crate::bake(&mut reference, loaded).unwrap();
        let score = || {
            Performance::read(
                "sample,event,key,value\n0,note_on,69,1\n\
            1200,sustain,0,0.4\n1200,note_off,69,0\n1800,sustain,0,0\n",
                &[69],
                2400,
            )
            .unwrap()
        };
        let mut whole = vec![0.; 4800];
        let mut split = whole.clone();
        let mut averaged_pressure = whole.clone();
        {
            let mut stream = ExteriorStream::new(&mut scene.piano, score(), &baked).unwrap();
            stream.render_interleaved_block(&mut whole).unwrap();
            assert_eq!(stream.sample_position(), 2400);
        }
        {
            let mut stream = ExteriorStream::new(&mut manual.piano, score(), &baked).unwrap();
            for block in split.chunks_mut(274) {
                stream.render_interleaved_block(block).unwrap();
            }
            assert_eq!(stream.sample_position(), 2400);
        }
        {
            let mut stream = ExteriorStream::new(&mut reference.piano, score(), &baked).unwrap();
            stream
                .render_interleaved_block(&mut averaged_pressure)
                .unwrap();
        }
        assert_eq!(whole, split);
        assert!(whole.iter().all(|x| x.is_finite()));
        assert!(whole.iter().any(|x| x.abs() > 1e-14));
        assert!(whole.chunks_exact(2).all(|f| f[0] == f[1]));
        assert_ne!(whole, averaged_pressure);
        assert_eq!(scene.piano.bank.q, manual.piano.bank.q);
        assert_eq!(scene.piano.bank.v, manual.piano.bank.v);
        assert_eq!(
            scene.piano.radiation_energy_j(),
            manual.piano.radiation_energy_j()
        );
        let mut moving_axial_coordinate = false;
        for string in 0..scene.piano.bank.strings.len() {
            let state = scene.piano.bank.string_longitudinal_state(string);
            assert_eq!(state, manual.piano.bank.string_longitudinal_state(string));
            if let Some((q, v)) = state {
                moving_axial_coordinate |= q.iter().chain(v).any(|x| *x != 0.);
            }
        }
        assert!(moving_axial_coordinate);
        assert!(scene.piano.accounting.felt_loss_j > 0.);
        assert!(scene.piano.accounting.damper_loss_j > 0.);
        if loaded {
            assert!(scene.piano.accounting.radiation_loss_j > 0.);
        }
    }
}

#[test]
fn missing_attachment_requires_geometry_or_an_explicit_legacy_comparison() {
    let make = |options: &Options, waves| {
        let (board, courses, obj, spec) = crate::tests::small_source_inputs();
        let controls = Controls::from_texts(&courses, None, None, None).unwrap()
            .with_string_stretching(&material(waves), &courses).unwrap();
        crate::prepare_controlled(&board, courses, &obj, spec, options, controls)
    };
    let parsed = |args: &[&str]| Options::parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    let options = parsed(&["--string-stretching", "waves.fspx", "--modes", "12", "--substeps", "8"]).unwrap();
    assert!(make(&options, true).err().unwrap().contains("no string-side stiffener"));
    let legacy = parsed(&["--string-stretching", "waves.fspx", "--legacy-axial-ratio", "--modes", "12", "--substeps", "8"]).unwrap();
    let old = make(&legacy, true).unwrap();
    assert!(legacy.report(&old.piano).contains("historical authored axial ratio 0.06"));
    assert!(make(&legacy, false).err().unwrap().contains("at least one waves key"));
    assert!(parsed(&["--legacy-axial-ratio"]).is_err());
    assert!(parsed(&["--legacy-axial-ratio", "--string-polarization", "frames", "--string-stretching", "waves"]).is_err());
}
