//! One hi-hat foot chick on the shared live hat, to a 48 kHz float WAV, the way QuarterTone drives a GM44
//! foot note: the hat settled at 80 N, the pedal released over 5 ms, then 40 ms later pressed back to 80 N
//! over RAMP_S (QuarterTone's 5 ms divided by the note's strength). The WAV starts at the closing command.
//! With STICK_M_S the pedal stays down and the stick strikes the same settled hat at 0.85 R instead, at the
//! same frame: the closed stroke that places the foot chick's level against the hat's stick strokes.
//! usage: footwav SCALE RAMP_S SECONDS OUT.wav [STICK_M_S]
use fs_vkplate::live::Live;

const RATE: usize = 48_000;
const BLOCK: usize = 120;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 && a.len() != 6 {
        eprintln!("usage: footwav SCALE RAMP_S SECONDS OUT.wav [STICK_M_S]");
        std::process::exit(2);
    }
    let scale: f64 = a[1].parse().expect("scale");
    let ramp: f64 = a[2].parse().expect("ramp seconds");
    let secs: f64 = a[3].parse().expect("seconds");
    let stick: Option<f64> = a.get(5).map(|s| s.parse().expect("stick speed"));
    let t0 = std::time::Instant::now();
    let mut p = Live::with_hihat_pedal(scale, 80.0).unwrap_or_else(|e| {
        eprintln!("refused: {e}");
        std::process::exit(1)
    });
    let release = RATE / 200; // 5 ms of rest before the release command
    let close = release + RATE * 40 / 1000; // the closing command, 40 ms after the release
    let total = close + (secs * RATE as f64) as usize;
    let mut out = Vec::with_capacity(total);
    let mut blk = [0.0f32; BLOCK];
    while out.len() < total {
        if let Some(speed) = stick {
            if out.len() == close {
                p.strike(0.85 * p.radius, 0.0, speed);
            }
        } else if out.len() == release {
            p.set_hihat_pedal(0.0, 0.005).expect("release");
        } else if out.len() == close {
            p.set_hihat_pedal(80.0, ramp).expect("close");
        }
        p.render(&mut blk, 1.0);
        out.extend(blk.iter().map(|&v| f64::from(v)));
    }
    let chick: Vec<f64> = out[close..total].to_vec();
    let peak = chick.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(peak.is_finite(), "non-finite output");
    let what = stick.map_or(format!("foot ramp {ramp:.4} s"), |v| {
        format!("closed stick {v} m/s")
    });
    eprintln!(
        "hihat {what} x{scale}: peak {peak:.3} Pa, {:.1} s wall",
        t0.elapsed().as_secs_f64()
    );
    std::fs::write(&a[4], fs_vkplate::wav_f32(&[chick], RATE as u32)).expect("write");
}
