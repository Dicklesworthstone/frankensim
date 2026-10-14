//! Snare-wire experiments: which wire physics makes the buzz?
use fs_vkplate::*;
use fs_vkplate::kit::Wires;
fn main() {
    let v = std::env::args().nth(1).unwrap_or_default();
    let path = std::env::args().nth(2).unwrap();
    let mut s = kit::snare();
    let mut w = Wires::estimated();
    let hard = Law::Power { k: 1.0e9, alpha: 1.5, mu: 0.1 };
    match v.as_str() {
        "hard" => w.contact = hard,
        "hard-light-preload" => { w.contact = hard; w.preload_per_m = 0.08; }
        "hard-tense" => { w.contact = hard; w.tension = 8.0; }
        "hard-light" => { w.contact = hard; w.mass_per_m = 0.0015; }
        "taut-light" => { w.tension = 8.0; w.mass_per_m = 0.0015; }
        "preload" => { w.preload_per_m = 0.8; }
        "slack" => { w.tension = 0.5; w.preload_per_m = 0.15; }
        "slack-heavy" => { w.tension = 1.0; w.mass_per_m = 0.008; }
        "mid" => { w.contact = Law::Power { k: 1.0e8, alpha: 1.5, mu: 0.5 }; w.tension = 8.0; w.mass_per_m = 0.0015; }
        _ => {}
    }
    s.wires = Some(w);
    let mut p = kit::drum(&s, &[(0.005, 4.5)]).unwrap();
    let out = kit::render(&mut p, 0.45, |_| {});
    std::fs::write(path, wav_f32(&out, 48000)).unwrap();
}
