//! Velocity fields of the crash and the snare batter shortly after the stick lands (share-image art).
//! usage: ogfield OUT.json
use fs_vkplate::*;
use std::fmt::Write;

fn snapshot(p: &mut kit::Piece, body: usize, t: f64) -> (Vec<[f64; 2]>, Vec<f64>, f64) {
    let ins = &mut p.instrument;
    let mut out = vec![0.0; ins.mics.len()];
    let steps = (t * ins.rate) as usize;
    for _ in 0..steps {
        ins.step(&mut out);
    }
    let b = &ins.bodies[body];
    let k = 1.0 / ins.rate;
    let v: Vec<f64> = (0..b.n).map(|i| (b.w[i] - b.wp[i]) / k).collect();
    let r = b.xy.iter().map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt()).fold(0.0, f64::max);
    (b.xy.clone(), v, r)
}

fn main() {
    let path = std::env::args().nth(1).expect("OUT.json");
    let mut crash = kit::cymbal(&kit::crash(), &[(0.0005, 6.0, 0.9, 0.008)]).unwrap();
    let tc: f64 = std::env::var("OG_TC").ok().and_then(|s| s.parse().ok()).unwrap_or(0.002);
    let ts: f64 = std::env::var("OG_TS").ok().and_then(|s| s.parse().ok()).unwrap_or(0.0035);
    let (cxy, cv, cr) = snapshot(&mut crash, 0, tc);
    let mut snare = kit::drum(&kit::snare_wires_off(), &[(0.0005, 4.5)]).unwrap();
    let (sxy, sv, sr) = snapshot(&mut snare, 0, ts);
    let mut s = String::from("{");
    for (name, xy, v, r) in [("crash", &cxy, &cv, cr), ("snare", &sxy, &sv, sr)] {
        write!(s, "\"{name}\":{{\"r\":{r},\"x\":[").unwrap();
        s += &xy.iter().map(|p| format!("{:.4}", p[0])).collect::<Vec<_>>().join(",");
        s += "],\"y\":[";
        s += &xy.iter().map(|p| format!("{:.4}", p[1])).collect::<Vec<_>>().join(",");
        s += "],\"v\":[";
        s += &v.iter().map(|x| format!("{:.4e}", x)).collect::<Vec<_>>().join(",");
        s += "]},";
    }
    s.pop();
    s += "}";
    std::fs::write(path, s).unwrap();
}
