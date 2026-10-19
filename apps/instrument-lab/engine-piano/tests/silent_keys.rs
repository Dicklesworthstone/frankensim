//! A key held down without a strike (engine.rs `silent_key_down`). Against the offline reference
//! engine, with C4 held silently and C3 struck, both bridges agree to roundoff: the live kernel
//! attaches C4's strings from rest and the board drives them as the reference's always-dynamic
//! strings. In the live budget the silent key stays attached while held, however quiet, and
//! returns to the static image after it is released.
use fs_live_piano::refmods::{engine, geometry, steinway_scale};
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

#[test]
fn a_silently_held_key_rings_along_as_in_the_reference() {
    let scale = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).unwrap();
    let courses = vec![scale[48 - 21], scale[60 - 21]];
    let desc = Description::with_courses(courses.clone(), 800.0).unwrap();
    let materials = courses.iter().map(|c| steinway_scale::hammer_material(c).unwrap()).collect();
    // 256 partials: the table reaches the partial band, so no cap binds (see tests/pedals.rs).
    let mut reference = engine::Instrument::new_with_string_damping(courses, &desc.board, RATE, 1, 256, true,
        materials, Some(engine::ShankGeometry::published()), None, None, false).unwrap();
    let mut cfg = Config::live(RATE);
    cfg.board_band_hz = 800.0;
    cfg.max_partials = 256;
    cfg.course_mode_budget = usize::MAX / 4;
    cfg.static_attachment = false;
    cfg.damper_points = 1;
    let mut live = LivePiano::new(&desc, cfg).unwrap();
    reference.silent_key_down(60).unwrap();
    live.silent_key_down(60).unwrap();
    reference.note_on(48, 4.0).unwrap();
    live.note_on(48, 4.0).unwrap();
    live.apply_pending_notes();
    assert_eq!(live.sounding_keys(), 2);
    let n = reference.bank.modes.len();
    let row = |course: usize| {
        let s = reference.bank.strings.iter().position(|s| s.course == course).unwrap();
        reference.bank.strings[s].bridge.clone()
    };
    let rows = [row(0), row(1)];
    let (mut peak, mut worst) = ([0.0_f64; 2], [0.0_f64; 2]);
    for _ in 0..RATE / 4 {
        reference.step().unwrap();
        live.step();
        for (k, midi) in [48u8, 60].into_iter().enumerate() {
            let vr: f64 = rows[k].iter().zip(&reference.bank.v[n..]).map(|(g, v)| g * v).sum();
            let vl = live.bridge_velocity(midi).unwrap();
            peak[k] = peak[k].max(vr.abs());
            worst[k] = worst[k].max((vr - vl).abs());
        }
    }
    assert!(live.take_error().is_none());
    for k in 0..2 {
        assert!(worst[k] < 1.0e-9 * peak[k], "bridge {k}: {:e} of {:e}", worst[k], peak[k]);
    }
}

#[test]
fn a_silent_key_stays_attached_while_held_and_returns_when_released() {
    let mut p = LivePiano::new(&Description::embedded().unwrap(), Config::live(RATE)).unwrap();
    p.silent_key_down(60).unwrap();
    p.note_on(48, 3.0).unwrap();
    let mut block = vec![0.0_f32; RATE as usize];
    p.render(&mut block);
    // Far past the 2048-sample quiet window: still attached, though it was never struck.
    assert_eq!(p.sounding_keys(), 2);
    p.note_off(60);
    for _ in 0..3 {
        p.render(&mut block);
    }
    assert_eq!(p.sounding_keys(), 1, "released, C4 returns to the static image once its dampers have quieted it");
    // A strike queued in the same sample wins over a silent press.
    p.note_on(64, 2.0).unwrap();
    p.silent_key_down(64).unwrap();
    p.render(&mut block[..64]);
    assert!(p.hammer_force(64).is_some());
    assert_eq!(p.counters.note_ons, 2);
}

#[test]
fn a_silent_key_stays_attached_while_anything_holds_its_dampers_up() {
    let mut block = vec![0.0_f32; RATE as usize / 4];
    // Pressed silently, caught by the sostenuto, released: nothing else plays, so it has no
    // energy at all, yet its dampers are up and it stays attached until the rod rises.
    let mut p = LivePiano::new(&Description::embedded().unwrap(), Config::live(RATE)).unwrap();
    p.silent_key_down(62).unwrap();
    p.render(&mut block[..64]);
    p.set_sostenuto(true);
    p.note_off(62);
    for _ in 0..4 {
        p.render(&mut block);
    }
    assert_eq!(p.sounding_keys(), 1, "caught by the sostenuto");
    p.set_sostenuto(false);
    p.render(&mut block);
    assert_eq!(p.sounding_keys(), 0, "the rod rose: damped and quiet, it returns");
    // The same under the full sustain pedal.
    p.set_sustain(1.0).unwrap();
    p.silent_key_down(62).unwrap();
    p.render(&mut block[..64]);
    p.note_off(62);
    p.render(&mut block);
    assert_eq!(p.sounding_keys(), 1, "held up by the sustain pedal");
    p.set_sustain(0.0).unwrap();
    p.render(&mut block);
    assert_eq!(p.sounding_keys(), 0);
    // Above E6 there is no damper: once released, a quiet silent key returns as any key does.
    p.silent_key_down(100).unwrap();
    p.render(&mut block[..64]);
    p.note_off(100);
    p.render(&mut block);
    assert_eq!(p.sounding_keys(), 0, "a released key without dampers is not kept forever");
}
