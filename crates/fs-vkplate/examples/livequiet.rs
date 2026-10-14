//! How a live piece dies away: kinetic energy and output level per 0.25 s after one stroke.
//! usage: livequiet PIECE [GRID_SCALE] [SECONDS]
use fs_vkplate::live::Live;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let name = a.get(1).map(String::as_str).unwrap_or("hihat-closed");
    let scale: f64 = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(5.0);
    let secs: f64 = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(6.0);
    let mut p = Live::new(name, scale).unwrap();
    // before any stroke: is the settled piece at rest?
    let mut blk = vec![0.0f32; 12000];
    p.render(&mut blk, 1.0);
    println!("unstruck: kinetic {:.3e} J, peak {:.3e} Pa", p.kinetic, blk.iter().fold(0.0f32, |m, v| m.max(v.abs())));
    p.strike(0.3 * p.radius, 0.0, 3.0);
    for q in 0..(secs * 4.0) as usize {
        p.render(&mut blk, 1.0);
        let pk = blk.iter().fold(0.0f32, |m, v| m.max(v.abs()));
        println!("{:5.2} s  kinetic {:.3e} J  peak {:.3e} Pa  asleep {}", 0.25 * (q + 1) as f64, p.kinetic, pk, p.asleep());
    }
}
