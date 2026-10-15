//! One kit piece on a coarser grid for a profiler: `liveprof PIECE GRID_SCALE SECONDS`.
use fs_vkplate::kit::{self, CymbalSpec, DrumSpec};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let s: f64 = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(3.5);
    let secs: f64 = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(4.0);
    let coarse = |mut d: DrumSpec| { d.grid *= s; d.oversample = 1; d };
    let cym = |mut c: CymbalSpec| { c.grid *= s; c };
    let hits: Vec<(f64, f64)> = (0..(secs * 2.0) as usize).map(|k| (0.005 + 0.5 * k as f64, 3.0)).collect();
    let mut p = match a.get(1).map(String::as_str).unwrap_or("kick") {
        "kick" => kit::drum(&coarse(kit::kick()), &hits),
        "snare" => kit::drum(&coarse(kit::snare()), &hits),
        "tom-12" => kit::drum(&coarse(kit::tom(12.0)), &hits),
        "crash" => kit::cymbal(&cym(kit::crash()), &hits.iter().map(|&(t, v)| (t, v, 0.95, 0.008)).collect::<Vec<_>>()),
        other => panic!("piece {other}"),
    }
    .unwrap();
    let nodes: usize = p.instrument.bodies.iter().map(|b| b.nodes()).sum();
    let t = Instant::now();
    let out = kit::render(&mut p, secs, |_| {});
    let wall = t.elapsed().as_secs_f64();
    println!("{nodes} nodes, {:.2}x real time, peak {:.2}", secs / wall, out[0].iter().fold(0.0f64, |m, x| m.max(x.abs())));
}
