//! Closed hi-hat experiments: rim contact stiffness, loss and contact area.
use fs_vkplate::*;
fn main() {
    let v = std::env::args().nth(1).unwrap_or_default();
    let path = std::env::args().nth(2).unwrap();
    let (law, from) = match v.as_str() {
        "stiff-lossy" => (Law::Power { k: 2.0e9, alpha: 1.5, mu: 3.0 }, 0.8),
        "wide" => (Law::Power { k: 2.0e9, alpha: 1.5, mu: 3.0 }, 0.6),
        "stiff-elastic" => (Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.2 }, 0.8),
        _ => (Law::Power { k: 3.0e8, alpha: 1.5, mu: 0.8 }, 0.8),
    };
    let mut p = kit::hihat_with(0.0, 80.0, &[(0.005, 2.5)], law, from).unwrap();
    let out = kit::render(&mut p, 0.4, |_| {});
    std::fs::write(path, wav_f32(&out, 48000)).unwrap();
}
