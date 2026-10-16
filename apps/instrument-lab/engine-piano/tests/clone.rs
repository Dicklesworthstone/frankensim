//! A copy of the live piano renders bit for bit what the original would: a
//! copy of a fresh engine (what the browser plays after timing the original)
//! and a copy taken mid-ring, with keys sounding, dampers and the pedal moving.
use fs_live_piano::{Config, Description, LivePiano};

fn play(piano: &mut LivePiano, from: usize, blocks: usize) -> Vec<u32> {
    let mut out = Vec::new();
    let mut block = [0.0_f32; 128];
    for b in from..from + blocks {
        // a bass note, a four-key chord and a run, released and re-struck
        match b % 40 {
            0 => {
                for (k, v) in [(36_u8, 3.0), (55, 2.6), (59, 2.6), (62, 2.6), (66, 2.4)] {
                    piano.note_on(k, v).unwrap();
                }
            }
            10 => piano.set_sustain(0.8).unwrap(),
            20 => {
                for k in [55_u8, 59, 62, 66] {
                    piano.note_off(k);
                }
                piano.note_on(74 + (b / 40 % 5) as u8, 2.2).unwrap();
            }
            30 => {
                piano.set_sustain(0.0).unwrap();
                piano.note_off(36);
            }
            _ => {}
        }
        piano.render(&mut block);
        out.extend(block.iter().map(|x| x.to_bits()));
    }
    out
}

fn books(piano: &LivePiano) -> String {
    format!(
        "{:?} {:?} {} {}",
        piano.ledger,
        piano.counters,
        piano.sounding_keys(),
        piano.dynamic_modes()
    )
}

#[test]
fn a_copy_renders_what_the_original_would() {
    let desc = Description::embedded().unwrap();
    let mut cfg = Config::live(48_000);
    cfg.max_dynamic_modes = cfg.max_dynamic_modes.min(1200);
    // fresh: a copy of a new engine against a second new engine
    let mut fresh = LivePiano::new(&desc, cfg.clone()).unwrap();
    let mut copy = fresh.clone();
    let mut reference = LivePiano::new(&desc, cfg.clone()).unwrap();
    // the browser lowers the cap after building, then copies
    for piano in [&mut fresh, &mut reference] {
        piano.set_quality(40, cfg.partial_band_hz, 900).unwrap();
    }
    copy.set_quality(40, cfg.partial_band_hz, 900).unwrap();
    let a = play(&mut copy, 0, 150);
    assert_eq!(a, play(&mut reference, 0, 150), "copy of a fresh engine");
    assert_eq!(books(&copy), books(&reference));
    assert!(a.iter().any(|&x| f32::from_bits(x).abs() > 1e-4), "silent");
    // mid-ring: keys sounding and damping, the pedal down
    let _ = play(&mut fresh, 0, 150);
    let mut ringing = fresh.clone();
    assert!(ringing.sounding_keys() > 0);
    assert_eq!(
        play(&mut ringing, 150, 200),
        play(&mut fresh, 150, 200),
        "copy taken mid-ring"
    );
    assert_eq!(books(&ringing), books(&fresh));
}
