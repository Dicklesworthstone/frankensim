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
        let make = |selection: Option<&str>| {
            let (board, courses, obj, spec) = crate::tests::small_source_inputs();
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
            crate::prepare_controlled(&board, courses, &obj, spec, &options, controls).unwrap()
        };
        let waves = material(true);
        let mut scene = make(Some(&waves));
        let mut manual = make(None);
        let mut reference = make(Some(&material(false)));
        // Attach through the direct physical API independently of the inline
        // Controls path used by render preparation.
        let courses = crate::tests::small_source_inputs().1;
        manual
            .piano
            .configure_string_stretching(
                &linear::string_stretching::Specification::read(&waves, &courses).unwrap(),
            )
            .unwrap();
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
