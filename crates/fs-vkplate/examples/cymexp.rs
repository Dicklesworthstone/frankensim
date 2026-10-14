//! Cymbal grid-convergence experiment: the same ride at two grid spacings.
use fs_vkplate::*;
fn main() {
    let grid: f64 = std::env::args().nth(1).unwrap().parse().unwrap();
    let path = std::env::args().nth(2).unwrap();
    let mut c = kit::ride();
    c.grid = grid;
    let mut p = kit::cymbal(&c, &[(0.005, 2.5, 0.6, 0.004)]).unwrap();
    eprintln!("{:?} os {}", p.notes, p.oversample);
    let out = kit::render(&mut p, 0.6, |_| {});
    std::fs::write(path, wav_f32(&out, 48000)).unwrap();
}
