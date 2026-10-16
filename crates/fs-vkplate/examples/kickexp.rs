//! Which mechanism puts the kick's energy at ~170 Hz? One variant per run, mic spectra printed by the caller.
use fs_vkplate::*;
fn main() {
    let v = std::env::args().nth(1).unwrap_or_default();
    let out_path = std::env::args().nth(2).unwrap();
    let mut s = kit::kick();
    match v.as_str() {
        "noport" => s.port = None,
        "soft" => s.striker.1 = Law::Power { k: 3.0e7, alpha: 2.5, mu: 0.3 },
        "nomuffle" => s.muffle = None,
        "batterside" => s.mic = [0.0, 0.1, 1.0],
        "offcentre" => s.strike_xy = (0.07, 0.0),
        _ => {}
    }
    let mut p = kit::drum(&s, &[(0.005, 3.0)]).unwrap();
    let out = kit::render(&mut p, 0.6, |_| {});
    std::fs::write(out_path, wav_f32(&out, 48000)).unwrap();
}
