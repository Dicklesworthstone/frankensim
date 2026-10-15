//! Can each kit piece run in real time on a coarser grid? Native timing at 48 kHz output.
//! usage: livebench [grid_scale]
use fs_vkplate::kit::{self, CymbalSpec, DrumSpec};
use std::time::Instant;

fn time(name: &str, mut p: kit::Piece, secs: f64) {
    let rate = p.instrument.rate;
    let nodes: usize = p.instrument.bodies.iter().map(|b| b.nodes()).sum();
    let t = Instant::now();
    let out = kit::render(&mut p, secs, |_| {});
    let wall = t.elapsed().as_secs_f64();
    let peak = out[0].iter().fold(0.0f64, |m, x| m.max(x.abs()));
    println!("{name:14} nodes {nodes:6} internal {:4.0} kHz  {:6.2}x real time  peak {peak:.2} Pa", rate / 1000.0, secs / wall);
}

fn main() {
    let s: f64 = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(2.5);
    let coarse = |mut d: DrumSpec| { d.grid *= s; d.oversample = 1; d };
    let cym = |mut c: CymbalSpec| { c.grid *= s; c };
    time("kick", kit::drum(&coarse(kit::kick()), &[(0.005, 3.0)]).unwrap(), 0.5);
    time("snare", kit::drum(&coarse(kit::snare()), &[(0.005, 4.5)]).unwrap(), 0.5);
    time("snare-off", kit::drum(&coarse(kit::snare_wires_off()), &[(0.005, 4.5)]).unwrap(), 0.5);
    time("tom-12", kit::drum(&coarse(kit::tom(12.0)), &[(0.005, 3.5)]).unwrap(), 0.5);
    time("tom-14", kit::drum(&coarse(kit::tom(14.0)), &[(0.005, 3.5)]).unwrap(), 0.5);
    time("crash", kit::cymbal(&cym(kit::crash()), &[(0.005, 6.0, 0.95, 0.008)]).unwrap(), 0.5);
    time("ride", kit::cymbal(&cym(kit::ride()), &[(0.005, 2.5, 0.6, 0.004)]).unwrap(), 0.5);
}
