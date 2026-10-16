//! Leapfrog CFL use of each cymbal shell with and without its linearised shell (curvature) stiffness.
//! usage: cfl PIECE SCALE...
use fs_vkplate::kit;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    for s in &a[2..] {
        let scale: f64 = s.parse().unwrap();
        let mut p = match a[1].as_str() {
            "crash" => { let mut c = kit::crash(); c.grid *= scale; kit::cymbal(&c, &[]).unwrap() }
            "ride" => { let mut c = kit::ride(); c.grid *= scale; kit::cymbal(&c, &[]).unwrap() }
            _ => kit::hihat_grid(0.002, 0.0, &[], fs_vkplate::Law::Power { k: 2e9, alpha: 1.5, mu: 3.0 }, 0.8, scale).unwrap(),
        };
        let k = 1.0 / p.instrument.rate;
        for b in p.instrument.bodies.iter_mut().filter(|b| b.n > 1) {
            let l0 = b.max_rate2();
            let l1 = b.max_rate2_with_shell();
            println!("{} x{scale} rate {} kHz: K alone {:.3}  K+B {:.3}   (f_max {:.0} / {:.0} Hz)", b.name, 1e-3 / k, k * k * l0 / 4.0, k * k * l1 / 4.0, l0.sqrt() / (2.0 * std::f64::consts::PI), l1.sqrt() / (2.0 * std::f64::consts::PI));
        }
    }
}
