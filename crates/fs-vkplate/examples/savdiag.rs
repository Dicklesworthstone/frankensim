//! SAV drift diagnostic on a live piece: psi against sqrt(2 V(w)) and kinetic energy over time.
//! usage: savdiag PIECE SCALE SPEED SECONDS
use fs_vkplate::live::Live;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let name = a[1].as_str();
    let scale: f64 = a[2].parse().unwrap();
    let speed: f64 = a[3].parse().unwrap();
    let secs: f64 = a[4].parse().unwrap();
    let mut p = Live::new(name, scale).unwrap();
    let rf = if name == "ride" { 0.6 } else { 0.85 };
    p.strike(rf * p.radius, 0.0, speed);
    let mut blk = vec![0.0f32; 1200];
    for q in 0..(secs / 0.025) as usize {
        p.render(&mut blk, 1.0);
        let pk = blk.iter().fold(0.0f32, |m, v| m.max(v.abs()));
        if q % 4 == 3 {
            let k = p.ins.k();
            let mut line = format!("{:5.2}s pk {pk:.3e} ", 0.025 * (q + 1) as f64);
            for b in p.ins.bodies.iter_mut().filter(|b| b.n > 1) {
                let ke: f64 = b.mass.iter().zip(&b.w).zip(&b.wp).map(|((m, w), wp)| 0.5 * m * ((w - wp) / k).powi(2)).sum();
                let (psi, r) = b.sav_state();
                line += &format!("| KE {ke:.3e} psi {psi:.3e} sqrt2V {r:.3e} ratio {:.3} ", psi / r.max(1e-300));
            }
            println!("{line}");
        }
    }
}
