//! A cold voice-set change that fails after it has changed the set leaves the
//! engine exactly as it was (frankensim-ahffu). The failure used to leave the
//! Woodbury core sized for the old set, and the next sample indexed past it (a
//! panic, which aborts the browser's WASM engine).
use fs_live_piano::{Config, Description, LivePiano};

fn engine() -> LivePiano {
    let desc = Description::embedded().unwrap();
    let mut cfg = Config::live(48_000);
    cfg.max_dynamic_modes = cfg.max_dynamic_modes.min(1200);
    LivePiano::new(&desc, cfg).unwrap()
}

fn play(piano: &mut LivePiano, blocks: usize) -> Vec<u32> {
    let mut out = Vec::new();
    let mut block = [0.0_f32; 128];
    for _ in 0..blocks {
        piano.render(&mut block);
        out.extend(block.iter().map(|x| x.to_bits()));
    }
    out
}

/// The strike whose rebuild fails is refused, the old set steps bit for bit as
/// if the strike never came, and a later strike of the same key is admitted as
/// on an engine that never saw the failure.
#[test]
fn a_failed_rebuild_leaves_the_engine_as_it_was() {
    let (mut failed, mut reference) = (engine(), engine());
    for piano in [&mut failed, &mut reference] {
        piano.note_on(48, 2.5).unwrap();
    }
    assert!(play(&mut failed, 20) == play(&mut reference, 20));
    assert_eq!((failed.take_error(), reference.take_error()), (None, None));
    assert_eq!(failed.sounding_keys(), 1);
    failed.note_on(64, 2.5).unwrap();
    failed.fail_next_rebuild();
    // The strike is admitted (and its rebuild fails) at the next block's first sample.
    let mut after_failure = play(&mut failed, 1);
    assert_eq!(failed.take_error().as_deref(), Some("injected rebuild failure"));
    after_failure.extend(play(&mut failed, 39));
    assert!(after_failure == play(&mut reference, 40), "the old set must step as if the strike never came");
    assert_eq!(failed.take_error(), reference.take_error());
    assert_eq!(failed.sounding_keys(), 1);
    assert_eq!(format!("{:?}", failed.ledger), format!("{:?}", reference.ledger));
    for piano in [&mut failed, &mut reference] {
        piano.note_on(64, 2.5).unwrap();
    }
    let admitted = play(&mut failed, 60);
    assert!(admitted == play(&mut reference, 60), "a later strike must sound as on an untouched engine");
    assert_eq!(failed.take_error(), reference.take_error());
    assert_eq!(failed.sounding_keys(), 2);
    assert!(admitted.iter().any(|&bits| f32::from_bits(bits) != 0.0));
}
