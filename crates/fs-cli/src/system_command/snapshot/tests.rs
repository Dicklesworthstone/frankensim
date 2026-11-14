use super::*;

const FIXED: &str = include_str!("../../../../../examples/system/variable_transmission.json");
const ADAPTIVE: &str = include_str!("../../../../../examples/system/adaptive_oscillator.json");

fn with_cx(body: impl FnOnce(&fs_exec::Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0,
                kernel_id: 0x5953,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        body(&cx);
    });
}

fn scene(text: &str) -> SystemScene {
    SystemScene::from_text(text, "retained/scene.json").unwrap()
}

fn encode_image(scene: &SystemScene, image: &SceneImage) -> (Vec<u8>, [u8; 32]) {
    let sealed = image
        .seal_v2(
            &context::<SceneImage>(scene.checkpoint_problem()).unwrap(),
            limits(),
            || false,
        )
        .unwrap();
    let root = *sealed.authority_subject_receipt().id().as_bytes();
    (sealed.into_parts().0, root)
}

#[test]
fn native_scene_checkpoints_preserve_complete_fixed_and_adaptive_reports_and_forks() {
    for text in [FIXED, ADAPTIVE] {
        with_cx(|cx| {
            let mut uninterrupted = scene(text);
            for _ in 0..19 {
                assert!(uninterrupted.advance(f64::MAX, cx).unwrap());
            }
            let prefix = uninterrupted.render_json();
            let (bytes, root) = uninterrupted.checkpoint().unwrap();
            assert_eq!(prefix, uninterrupted.render_json());
            assert_eq!(uninterrupted.checkpoint().unwrap(), (bytes.clone(), root));

            let mut resumed = scene(text);
            resumed.advance(0.001, cx).unwrap();
            resumed.resume_checkpoint(&bytes, root).unwrap();
            let mut fork = scene(text);
            fork.resume_checkpoint(&bytes, root).unwrap();
            assert_eq!(resumed.render_json(), prefix);
            assert_eq!(fork.render_json(), prefix);
            uninterrupted.run(cx).unwrap();
            resumed.run(cx).unwrap();
            fork.run(cx).unwrap();
            assert_eq!(resumed.render_json(), uninterrupted.render_json());
            assert_eq!(fork.render_json(), uninterrupted.render_json());
            assert_eq!(resumed.render_text(), uninterrupted.render_text());

            // A completed checkpoint remains completed without a duplicate
            // final sample, and restoring it overwrites an accepted prefix.
            let (complete, complete_root) = resumed.checkpoint().unwrap();
            let mut finished = scene(text);
            finished
                .resume_checkpoint(&complete, complete_root)
                .unwrap();
            assert!(!finished.advance(1.0, cx).unwrap());
            assert_eq!(finished.render_json(), resumed.render_json());
        });
    }
}

#[test]
fn initial_unknown_samples_and_unsampled_scenes_round_trip() {
    let unsampled = FIXED.replace("  \"sample_every_s\": 0.25,\n", "");
    for text in [FIXED, ADAPTIVE, unsampled.as_str()] {
        let original = scene(text);
        let (bytes, root) = original.checkpoint().unwrap();
        let mut resumed = scene(text);
        resumed.resume_checkpoint(&bytes, root).unwrap();
        assert_eq!(resumed.render_json(), original.render_json());
        let original_bits: Vec<_> = original
            .device
            .samples
            .iter()
            .flat_map(|(_, row)| row.iter().map(|(_, value)| value.to_bits()))
            .collect();
        let resumed_bits: Vec<_> = resumed
            .device
            .samples
            .iter()
            .flat_map(|(_, row)| row.iter().map(|(_, value)| value.to_bits()))
            .collect();
        assert_eq!(original_bits, resumed_bits);
        with_cx(|cx| resumed.run(cx).unwrap());
        let mut reference = scene(text);
        with_cx(|cx| reference.run(cx).unwrap());
        assert_eq!(resumed.render_json(), reference.render_json());
    }
}

#[test]
fn foreign_boundary_responses_and_candidate_roots_refuse_without_mutation() {
    with_cx(|cx| {
        let mut original = scene(FIXED);
        original.advance(0.01, cx).unwrap();
        let (bytes, root) = original.checkpoint().unwrap();
        let mut target = scene(FIXED);
        target.advance(0.002, cx).unwrap();
        let before = target.render_json();
        let mut wrong_root = root;
        wrong_root[0] ^= 1;
        assert!(target.resume_checkpoint(&bytes, wrong_root).is_err());
        assert_eq!(before, target.render_json());
        let mut corrupt = bytes.clone();
        let last = corrupt.len() - 1;
        corrupt[last] ^= 1;
        assert!(target.resume_checkpoint(&corrupt, root).is_err());
        assert_eq!(before, target.render_json());
        assert!(
            target
                .resume_checkpoint(&bytes[..bytes.len() - 1], root)
                .is_err()
        );
        assert_eq!(before, target.render_json());

        let changed = FIXED.replace("[[0, 0.02]", "[[0, 0.03]");
        assert_ne!(changed, FIXED);
        let mut foreign = scene(&changed);
        // Kernel graphs deliberately exclude the opaque boundary functions;
        // the scene owner must bind their original authored input instead.
        assert_eq!(
            foreign.device.system.problem_id(),
            original.device.system.problem_id()
        );
        assert_ne!(foreign.checkpoint_problem(), original.checkpoint_problem());
        foreign.advance(0.002, cx).unwrap();
        let before = foreign.render_json();
        assert!(foreign.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(before, foreign.render_json());

        // The source's label is not its content: relocation is permitted.
        let mut relocated = SystemScene::from_text(FIXED, "elsewhere/scene.json").unwrap();
        relocated.resume_checkpoint(&bytes, root).unwrap();
        assert_eq!(relocated.time_s().to_bits(), original.time_s().to_bits());
        assert_eq!(relocated.steps(), original.steps());
    });
}

#[test]
fn admitted_but_malformed_driver_and_sample_states_refuse_atomically() {
    with_cx(|cx| {
        for text in [FIXED, ADAPTIVE] {
            let mut target = scene(text);
            for _ in 0..83 {
                target.advance(f64::MAX, cx).unwrap();
            }
            let valid = SceneImage::capture(&target).unwrap();
            let before = target.render_json();
            for mutation in 0..9 {
                let mut image = valid.clone();
                match mutation {
                    0 => image.iterations = u64::MAX,
                    1 => image.shortenings = u64::MAX,
                    2 => image.sample_width += 1,
                    3 => image.sample_times[0] = f64::NAN,
                    4 => image.next_sample_s = f64::INFINITY,
                    5 => image.sample_values[0] += 1.0,
                    6 => image.sample_times.push(image.system.time_s() + 1.0),
                    7 => image.sample_values.push(0.0),
                    8 => image.adaptive.push(0.0),
                    _ => unreachable!(),
                }
                let (bytes, root) = encode_image(&target, &image);
                let error = target.resume_checkpoint(&bytes, root).unwrap_err();
                assert_eq!(error.code, "system-checkpoint", "mutation {mutation}");
                assert_eq!(before, target.render_json(), "mutation {mutation}");
            }
            if !valid.adaptive.is_empty() {
                for (index, value) in [
                    (0, f64::from_bits(u64::MAX)),
                    (1, f64::from_bits(1)),
                    (2, f64::NAN),
                    (3, 2.0),
                    (4, 0.0),
                ] {
                    let mut image = valid.clone();
                    image.adaptive[index] = value;
                    let (bytes, root) = encode_image(&target, &image);
                    assert!(
                        target.resume_checkpoint(&bytes, root).is_err(),
                        "adaptive {index}"
                    );
                    assert_eq!(before, target.render_json());
                }
            }
        }
    });
}

#[test]
fn unsupported_owners_and_excessive_sample_storage_are_named_refusals() {
    let phase = scene(include_str!(
        "../../../../../examples/system/phase_change_store.json"
    ));
    let before = phase.render_json();
    let error = phase.checkpoint().unwrap_err();
    assert!(error.message.contains("store"), "{error}");
    assert!(error.message.contains("codec"), "{error}");
    assert_eq!(phase.render_json(), before);

    assert!(sample_storage(&["coordinate".repeat(1024)], MAX_SAMPLES).is_err());
    assert!(sample_storage(&["x".into()], usize::MAX).is_err());
}
