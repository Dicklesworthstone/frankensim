//! One live stroke (or a train of strokes) to a WAV, exactly as the browser lane renders it.
//! usage: livewav PIECE SCALE SPEED SECONDS OUT.wav [RADIUS_FRACTION] [PERIOD_S]
use fs_vkplate::live::Live;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let name = a[1].as_str();
    let scale: f64 = a[2].parse().unwrap();
    let speed: f64 = a[3].parse().unwrap();
    let secs: f64 = a[4].parse().unwrap();
    let path = &a[5];
    let rf: f64 = a.get(6).and_then(|v| v.parse().ok()).unwrap_or(match name {
        "ride" => 0.6,
        _ => 0.85,
    });
    let period: f64 = a.get(7).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let t0 = std::time::Instant::now();
    let mut p = Live::new(name, scale).unwrap_or_else(|e| panic!("{e}"));
    let n = (secs * 48000.0) as usize;
    let mut out = Vec::with_capacity(n);
    let mut blk = vec![0.0f32; 128];
    let mut next = 0.005;
    let mut t = 0.0;
    let mut peak_f = 0.0f64;
    while out.len() < n {
        if t >= next {
            p.strike(rf * p.radius, 0.0, speed);
            next = if period > 0.0 { next + period } else { f64::INFINITY };
        }
        p.render(&mut blk, 1.0);
        peak_f = peak_f.max(p.force);
        out.extend(blk.iter().map(|&v| v as f64));
        t += 128.0 / 48000.0;
    }
    out.truncate(n);
    let wall = t0.elapsed().as_secs_f64();
    let pk = out.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let rms50 = (out[240..240 + 2400].iter().map(|v| v * v).sum::<f64>() / 2400.0).sqrt();
    eprintln!("{name} x{scale}: nodes {} h {:.1} mm, peak {pk:.3} Pa, peak force {peak_f:.0} N, {:.1}x real time", p.surface_xy().len() / 2, p.surface_h() * 1e3, secs / wall);
    let _ = rms50;
    std::fs::write(path, fs_vkplate::wav_f32(&[out], 48000)).unwrap();
}
