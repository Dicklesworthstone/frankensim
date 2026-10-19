//! The live una corda and sostenuto pedals. Una corda against the offline reference engine
//! (same source files, one course on the bridge, as examples/validate.rs): the strikes agree
//! to roundoff. Sostenuto against its definition: a key the rod caught sounds exactly as a key
//! held down, and a key pressed after the rod rose is damped exactly as without it.
use fs_live_piano::refmods::{engine, geometry, steinway_scale};
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;
const PARTIALS: usize = 256;

/// Relative bridge-velocity difference (max over the run) and the bridge velocity's peak.
fn una_corda_run(midi: u8, una_corda: bool) -> (f64, f64, Vec<f64>) {
    let scale = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).unwrap();
    let course = scale[usize::from(midi - 21)];
    let desc = Description::with_courses(vec![course], 800.0).unwrap();
    let material = steinway_scale::hammer_material(&course).unwrap();
    // 256 partials: the table reaches the partial band, so no cap binds. A binding cap ends the
    // live duplex segments where the speaking strings end, which the reference does not do.
    let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, RATE, 1, PARTIALS, true,
        vec![material], Some(engine::ShankGeometry::published()), None, None, false).unwrap();
    let mut cfg = Config::live(RATE);
    cfg.board_band_hz = 800.0;
    cfg.max_partials = PARTIALS;
    cfg.course_mode_budget = usize::MAX / 4;
    cfg.static_attachment = false;
    cfg.damper_points = 1;
    let mut live = LivePiano::new(&desc, cfg).unwrap();
    reference.set_una_corda(una_corda);
    live.set_una_corda(una_corda);
    reference.note_on(midi, 2.5).unwrap();
    live.note_on(midi, 2.5).unwrap();
    live.apply_pending_notes();
    let row: Vec<f64> = reference.bank.strings[0].bridge.clone();
    let n = reference.bank.modes.len();
    let (mut peak, mut worst, mut trace) = (0.0_f64, 0.0_f64, Vec::new());
    for _ in 0..RATE / 5 {
        reference.step().unwrap();
        live.step();
        let vr: f64 = row.iter().zip(&reference.bank.v[n..]).map(|(g, v)| g * v).sum();
        let vl = live.bridge_velocity(midi).unwrap();
        peak = peak.max(vr.abs());
        worst = worst.max((vr - vl).abs());
        trace.push(vl);
    }
    assert!(live.take_error().is_none());
    (worst / peak, peak, trace)
}

#[test]
fn una_corda_strikes_match_the_reference_engine() {
    // A trichord keeps two strings, a bichord one, a single string its own.
    for (midi, changes) in [(57u8, true), (35, true), (25, false)] {
        let (rel, _, shifted) = una_corda_run(midi, true);
        assert!(rel < 1.0e-9, "key {midi}: live vs reference with una corda, {rel:e}");
        let (_, peak, normal) = una_corda_run(midi, false);
        let moved = shifted.iter().zip(&normal).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / peak;
        assert_eq!(moved > 1.0e-3, changes, "key {midi}: una corda moved the bridge by {moved:e}");
    }
}

fn piano() -> LivePiano {
    LivePiano::new(&Description::embedded().unwrap(), Config::live(RATE)).unwrap()
}
fn play(p: &mut LivePiano, seconds: f64) -> Vec<u32> {
    let mut out = vec![0.0_f32; (seconds * f64::from(RATE)) as usize];
    p.render(&mut out);
    out.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn a_key_caught_by_the_sostenuto_sounds_as_a_held_key() {
    let mut held = piano();
    held.note_on(57, 2.0).unwrap();
    let mut a = play(&mut held, 0.3);
    a.extend(play(&mut held, 1.0));

    let mut caught = piano();
    caught.note_on(57, 2.0).unwrap();
    let mut b = play(&mut caught, 0.3);
    caught.set_sostenuto(true);
    caught.note_off(57);
    b.extend(play(&mut caught, 1.0));
    assert!(a == b, "a latched key must ring exactly as a held one");

    let mut released = piano();
    released.note_on(57, 2.0).unwrap();
    let mut c = play(&mut released, 0.3);
    released.note_off(57);
    c.extend(play(&mut released, 1.0));
    let tail = |x: &[u32]| x[x.len() - 4800..].iter().map(|v| f64::from(f32::from_bits(*v)).powi(2)).sum::<f64>();
    assert!(tail(&a) > 100.0 * tail(&c), "the dampers must fall on a released key without the sostenuto");

    // The rod rising lets the dampers fall on the caught key.
    caught.set_sostenuto(false);
    assert!(play(&mut caught, 0.5) != play(&mut held, 0.5));
}

#[test]
fn the_sostenuto_catches_only_what_is_up_when_it_is_pressed() {
    let released = |prepare: &dyn Fn(&mut LivePiano)| {
        let mut p = piano();
        prepare(&mut p);
        p.note_on(57, 2.0).unwrap();
        let mut x = play(&mut p, 0.3);
        p.note_off(57);
        x.extend(play(&mut p, 1.0));
        x
    };
    let plain = released(&|_| {});
    // Pressed before the key went down: the key is not caught.
    assert!(released(&|p| p.set_sostenuto(true)) == plain);
    // Pressed while the sustain pedal is fully down: every damper is caught.
    let caught_by_sustain = released(&|p| {
        p.set_sustain(1.0).unwrap();
        p.set_sostenuto(true);
        p.set_sustain(0.0).unwrap();
    });
    let mut held = piano();
    held.note_on(57, 2.0).unwrap();
    let held = play(&mut held, 1.3);
    assert!(caught_by_sustain == held);
    // Half a sustain pedal raises no damper fully: nothing is caught.
    assert!(released(&|p| {
        p.set_sustain(0.5).unwrap();
        p.set_sostenuto(true);
        p.set_sustain(0.0).unwrap();
    }) == plain);
    // Every key up lifts the sostenuto too.
    let mut p = piano();
    p.set_sostenuto(true);
    p.set_una_corda(true);
    p.all_notes_off();
    assert!(!p.sostenuto() && !p.una_corda());
}
