//! Hi-hat with rim friction: one stick stroke on a grid `scale` times the offline one.
//! usage: hatfric SCALE FRICTION RIM_MU SECONDS OUT.wav [SPEED_M_S] [open]
//! (RIM_MU: the rims' Hunt-Crossley loss, 0.08 offline, 3.0 in the live module)
use fs_vkplate::*;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 6 {
        eprintln!("usage: hatfric SCALE FRICTION RIM_MU SECONDS OUT.wav [SPEED_M_S] [open]");
        std::process::exit(2);
    }
    let scale: f64 = a[1].parse().unwrap();
    let friction: f64 = a[2].parse().unwrap();
    let rim_mu: f64 = a[3].parse().unwrap();
    let secs: f64 = a[4].parse().unwrap();
    let speed: f64 = a.get(6).and_then(|v| v.parse().ok()).unwrap_or(2.5);
    let rim = Law::Power { k: 2.0e9, alpha: 1.5, mu: rim_mu };
    let t0 = std::time::Instant::now();
    let open = a.get(7).is_some_and(|v| v == "open");
    let (gap, pedal) = if open { (0.002, 0.0) } else { (0.0, 80.0) };
    let mut p = kit::hihat_grid_mu(gap, pedal, &[(0.005, speed)], rim, 0.8, scale, friction).unwrap_or_else(|e| panic!("{e}"));
    for n in &p.notes {
        eprintln!("{n}");
    }
    let out = kit::render(&mut p, secs, |_| {});
    let (mut stuck, mut fmax) = (0usize, 0.0f64);
    for f in &p.instrument.frictions {
        stuck += f.stuck as usize;
        fmax = fmax.max(f.force);
    }
    eprintln!("{} frictions ({stuck} stuck at the end, largest force {fmax:.2} N); {:.1} s wall for {secs} s", p.instrument.frictions.len(), t0.elapsed().as_secs_f64());
    std::fs::write(&a[5], wav_f32(&out, 48000)).unwrap();
}
