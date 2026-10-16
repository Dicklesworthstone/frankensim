//! G3: finite-gap felt contact drives the radiation-loaded BEM stereo path.
use super::*;
use crate::{exterior_audio::ExteriorStream, performance::Performance};

#[test]
fn supplied_felt_pads_reach_loaded_bem_and_preserve_split_block_state() {
    let options = Options {
        modes: 12,
        substeps: 8,
        ..Options::default()
    };
    let card = "# Explicit numerical fixture.\n\
        frankensim-piano-felt-dampers-v1\n\
        pad,69,0.12,0.14,0.00001,0.004,0,0.001,0.05,1,0,2\n\
        felt,69,60000,0.2,2.5,3.2,0.1,0.8,300000\n\
        branch,69,200000,0.002\n";
    let build = |dampers| {
        let (board, courses, obj, spec) = crate::tests::small_source_inputs();
        let controls = Controls::from_texts(
            &courses,
            None,
            Some("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n"),
            dampers,
        )
        .unwrap();
        crate::prepare_controlled(&board, courses, &obj, spec, &options, controls).unwrap()
    };
    let mut scene = build(Some(card));
    let mut direct = build(None);
    let mut free = build(Some("frankensim-piano-felt-dampers-v1\nfree,69\n"));
    let courses = crate::tests::small_source_inputs().1;
    direct
        .piano
        .configure_damper_spec(&engine::DamperSpecification::read(card, &courses).unwrap())
        .unwrap();
    assert_eq!(scene.piano.felt_damper_resolution(), Some((3, 6)));
    assert!(scene.piano.damper_resolution().is_none());
    assert_eq!(scene.piano.hammer_contact_count(), 6);
    let (baked, _, _) = crate::bake(&mut scene, true).unwrap();
    crate::bake(&mut direct, true).unwrap();
    crate::bake(&mut free, true).unwrap();
    let score = || {
        Performance::read(
            "sample,event,key,value\n0,note_on,69,1\n\
            480,note_off,69,0\n480,sustain,0,0.25\n\
            960,sustain,0,0\n1800,sustain,0,1\n",
            &[69],
            2400,
        )
        .unwrap()
    };
    let mut whole = vec![0.; 4800];
    let mut split = whole.clone();
    let mut free_pressure = whole.clone();
    {
        let mut stream = ExteriorStream::new(&mut scene.piano, score(), &baked).unwrap();
        stream.render_interleaved_block(&mut whole).unwrap();
        assert_eq!(stream.sample_position(), 2400);
    }
    {
        let mut stream = ExteriorStream::new(&mut direct.piano, score(), &baked).unwrap();
        for block in split.chunks_mut(274) {
            stream.render_interleaved_block(block).unwrap();
        }
        assert_eq!(stream.sample_position(), 2400);
    }
    {
        let mut stream = ExteriorStream::new(&mut free.piano, score(), &baked).unwrap();
        stream.render_interleaved_block(&mut free_pressure).unwrap();
    }
    assert_eq!(whole, split);
    assert!(whole.iter().all(|x| x.is_finite()));
    assert!(whole.iter().any(|x| x.abs() > 1e-14));
    assert_ne!(whole, free_pressure);
    assert_eq!(scene.piano.bank.q, direct.piano.bank.q);
    assert_eq!(scene.piano.bank.v, direct.piano.bank.v);
    assert_eq!(
        scene.piano.damper_energy_j(),
        direct.piano.damper_energy_j()
    );
    assert_eq!(
        scene.piano.radiation_energy_j(),
        direct.piano.radiation_energy_j()
    );
    assert_eq!(
        scene.piano.accounting.damper_loss_j,
        direct.piano.accounting.damper_loss_j
    );
    assert_eq!(
        scene.piano.accounting.damper_actuator_work_j,
        direct.piano.accounting.damper_actuator_work_j
    );
    assert_eq!(
        scene.piano.accounting.damper_relaxation_loss_j,
        direct.piano.accounting.damper_relaxation_loss_j
    );
    assert!(scene.piano.accounting.damper_loss_j > 0.0);
    assert!(scene.piano.accounting.damper_relaxation_loss_j > 0.0);
    assert_ne!(scene.piano.accounting.damper_actuator_work_j, 0.0);
    assert!(scene.piano.accounting.radiation_loss_j > 0.0);
    assert_eq!(free.piano.accounting.damper_loss_j, 0.0);
    for piano in [&scene.piano, &direct.piano, &free.piano] {
        assert!(
            (piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j())
                .abs()
                < 1e-7
        );
    }
    let report = options.report(&scene.piano);
    assert!(report.contains("finite-gap felt, 3 speaking-string pads / 6 contact sites"));
    assert!(report.contains("Prony relaxation"));
    assert!(report.contains("signed actuator work"));
}
