//! Diagnostic: uniform versus physical two-lobe rim gap on the same coarse grid.
//! Does not change kit defaults. Output directories must be new, and are retained.
use fs_vkplate::{kit, wav_f32};
use std::{fs, io::Write, time::Instant};

fn main() {
    let destination = std::env::args().nth(1).expect("output directory");
    let geometry = std::env::args().nth(2).unwrap_or_else(|| "waviness".into());
    assert!(matches!(geometry.as_str(), "waviness" | "bow"));
    fs::create_dir(&destination).expect("new output directory");
    for (name, gap, pedal) in [("closed", 0.0, 80.0), ("open", 0.002, 0.0)] {
        for scale in [4.5, 3.0] {
            let original = kit::hihat_scaled(gap, pedal, &[(0.005, 2.5)], scale).unwrap();
            assert_eq!(original.oversample, 1);
            let effective_gap = if pedal == 0.0 { gap.min(1.4e-4) } else { gap };
            for rep in 0..2 {
                for wavy in if rep == 0 { [false, true] } else { [true, false] } {
                    let mut p = kit::Piece {
                        instrument: original.instrument.clone(),
                        oversample: original.oversample,
                        settle: original.settle,
                        notes: original.notes.clone(),
                    };
                    let mut changed = 0;
                    for c in &mut p.instrument.contacts {
                        if c.a.len() != 1 || c.b.len() != 1 || c.a[0].0 != 1 || c.b[0].0 != 0 {
                            continue;
                        }
                        assert_eq!(c.gap, effective_gap);
                        if wavy {
                            let [x, y] = p.instrument.bodies[0].xy[c.b[0].1];
                            let r2 = x*x + y*y;
                            let cos2 = if r2 > 0.0 { (x*x-y*y)/r2 } else { 1.0 };
                            if geometry == "waviness" {
                                c.gap += 7.5e-5 * (1.0-cos2);
                            } else {
                                // The opposed bows separate toward the centre; use the
                                // declared rest shape, without changing its stiffness.
                                let (top, bottom) = kit::hihat_spec();
                                for shell in [top, bottom] {
                                    let u = (1.0-r2.sqrt()/shell.radius).max(0.0);
                                    c.gap += shell.bow_rise * (0.7*u+0.3*u*u);
                                }
                            }
                        }
                        changed += 1;
                    }
                    assert!(changed > 0);
                    let before = Instant::now();
                    let out = kit::render(&mut p, 0.6, |_| {});
                    let seconds = before.elapsed().as_secs_f64();
                    assert_eq!(out[0].len(), 28_800);
                    assert!(out[0].iter().all(|x| x.is_finite()));
                    let peak = out[0].iter().map(|x| x.abs()).fold(0.0, f64::max);
                    assert!(peak > 1e-5);
                    let path = format!("{destination}/{name}-{scale}-{wavy}-{rep}.wav");
                    let mut file = fs::OpenOptions::new().create_new(true).write(true).open(path).unwrap();
                    file.write_all(&wav_f32(&out, 48_000)).unwrap();
                    println!("{{\"piece\":\"{name}\",\"scale\":{scale},\"wavy\":{wavy},\"rep\":{rep},\"contacts\":{changed},\"seconds\":{seconds},\"peak_pa\":{peak}}}");
                    std::io::stdout().flush().unwrap();
                }
            }
        }
    }
}
