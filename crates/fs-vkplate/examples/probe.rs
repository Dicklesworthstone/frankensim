//! Diagnostic: per-body net volume and port flow of a drum, written as CSV (one row per output sample).
use fs_vkplate::*;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let spec = match a.get(1).map(|s| s.as_str()).unwrap_or("kick") {
        "kick" => kit::kick(),
        "tom-12" => kit::tom(12.0),
        _ => kit::snare_wires_off(),
    };
    let secs: f64 = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.5);
    let mut p = kit::drum(&spec, &[(0.005, 3.0)]).unwrap();
    let ins = &mut p.instrument;
    let mut out = vec![0.0; ins.mics.len()];
    println!("t,batter_volume,reso_volume,port_volume,mic,strike_force");
    let n = (secs * ins.rate) as usize;
    for _ in 0..n {
        ins.step(&mut out);
        let vol = |b: &Body| b.w.iter().sum::<f64>() * b.h * b.h;
        let px = ins.cavities[0].port.as_ref().map_or(0.0, |q| q.x);
        let fs = ins.contacts.last().map_or(0.0, |c| c.force);
        println!("{:.6},{:.6e},{:.6e},{:.6e},{:.6e},{:.4e}", ins.time(), vol(&ins.bodies[0]), vol(&ins.bodies[1]), px, out[0], fs);
    }
}
