//! Step-0 measurements that decide the live design (run natively):
//! (b) `PianoVertical` fixture speed, cold board preparation, and the cost of
//! the small dense eigen-solves a dynamic voice set needs at each note-on.
use fs_couple::piano_vertical::{HammerLaw, PedalState, PianoStringSpec, PianoVertical, RATE};
use std::time::Instant;

fn main() {
    // (b) PianoVertical: one note, 24 kHz fixture clock.
    for n_modes in [8usize, 16, 32] {
        let spec = PianoStringSpec { f0_hz: 261.6, b_inharmonicity: 4.0e-4, detune_cents: 0.0,
            n_modes, damping_ratio: 2.0e-4 };
        let mut p = PianoVertical::new(spec, 0.8, HammerLaw::Felt,
            PedalState { sustain: true, una_corda: false }, 0.02).expect("fixture");
        p.strike(2.5);
        let steps = 24_000usize;
        let t = Instant::now();
        let mut peak = 0.0_f64;
        for _ in 0..steps { peak = peak.max(p.step().expect("step").abs()); }
        let s = t.elapsed().as_secs_f64();
        println!("piano_vertical n_modes={n_modes}: {steps} steps in {s:.4} s = {:.0} steps/s; real-time factor at its {RATE} Hz clock {:.2}x; peak {peak:.3e}",
            steps as f64 / s, steps as f64 / s / f64::from(RATE));
    }
    // Cold board preparation at several bands.
    for band in [400.0, 800.0, 1000.0, 1200.0] {
        let t = Instant::now();
        let d = fs_live_piano::Description::model_d(band).expect("board");
        println!("board band {band} Hz: {} modes, {} surface samples, f0 {:.3} Hz; {:.4} s",
            d.board.len(), d.surface.len(), d.board[0].frequency_hz, t.elapsed().as_secs_f64());
    }
    // Dense generalized eigen-solve cost (voice-set change).
    for r in [8usize, 16, 25, 32, 43, 62] {
        let mut k = vec![0.0; r * r];
        let mut m = vec![0.0; r * r];
        for i in 0..r {
            for j in 0..r {
                let g = ((i * 7 + j * 3) % 11) as f64 * 0.01;
                k[i * r + j] = g * 1e4; k[j * r + i] = g * 1e4;
                m[i * r + j] = g * 0.01; m[j * r + i] = g * 0.01;
            }
            k[i * r + i] = 1e6 * (i + 1) as f64; m[i * r + i] = 1.0;
        }
        let t = Instant::now();
        let reps = 20;
        for _ in 0..reps { let e = fs_modal::eigh_gen_dense(&k, &m, r).expect("eig"); std::hint::black_box(e); }
        println!("eigh_gen_dense r={r}: {:.1} us", t.elapsed().as_secs_f64() / reps as f64 * 1e6);
    }
}
