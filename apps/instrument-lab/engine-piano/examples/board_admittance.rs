//! A bridge-informed board reduction against the FE slice it reduces (run NATIVELY).
//!
//! For each reduction spec and each listed key's bridge point, in octave bands
//! from 125 Hz to 8 kHz, in dB against the complete source slice:
//! - G: driving-point conductance Re Y (Y = bridge velocity per bridge force,
//!   averaged over linear frequency), the power a string can feed the board;
//! - P: Rayleigh half-space pressure per bridge force, mean-square over the
//!   listed microphones (default: the live engine's), what the board radiates.
//!
//! and the mass line sum(b_j^2): above every board mode a unit bridge force
//! accelerates the bridge by it, so it sets G and P above the slice's band.
//! Both models carry the source damping: the slice its modal ratios, the
//! reduction its full projected matrix.
//!
//! usage: board_admittance [--mesh 16,0.03|preset] [--band HZ] [--keys 33,45,...]
//!        [--mics X,Y,Z;X,Y,Z;...] SPEC...
use fs_la::eigen_complex::lu_complex;
use fs_live_piano::refmods::{board_geometry, geometry, linear::BoardMode, steinway_d};
use fs_live_piano::{Config, Description};
use fs_math::c64::C64;
use std::f64::consts::TAU;

/// Air at 20 C. The pressure ratios depend on it only through travel phases.
const SOUND_SPEED_M_S: f64 = 343.2;
const OCTAVES_HZ: [f64; 7] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0];
const POINTS_PER_OCTAVE: usize = 48;

/// One board in mass-normalized coordinates: stiffness, damping, bridge
/// coefficients per key, and each microphone's radiation weight per mode at
/// each analysis frequency.
struct Board {
    lambda: Vec<f64>,
    /// Row-major viscous matrix [1/s].
    damping: Vec<f64>,
    bridge: Vec<[f64; 88]>,
    /// [microphone][frequency][mode]
    radiation: Vec<Vec<Vec<C64>>>,
}

impl Board {
    fn new(modes: &[BoardMode], damping: Option<&[f64]>, surface: &[board_geometry::SurfaceSample],
        mics: &[[f64; 3]], hz: &[f64]) -> Self {
        let r = modes.len();
        let lambda: Vec<f64> = modes.iter().map(|m| (TAU * m.frequency_hz).powi(2)).collect();
        let damping = damping.map_or_else(|| {
            let mut c = vec![0.0; r * r];
            for (j, m) in modes.iter().enumerate() { c[j * r + j] = 2.0 * m.damping_ratio * lambda[j].sqrt(); }
            c
        }, <[f64]>::to_vec);
        // Rayleigh: p ~ sum_s A_s a_s exp(i k R_s) / R_s (exp(-i w t)); constants cancel in the ratios.
        let radiation = mics.iter().map(|mic| {
            let distance: Vec<f64> = surface.iter().map(|s| {
                let d = [s.position_m[0] - mic[0], s.position_m[1] - mic[1], s.position_m[2] - mic[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            }).collect();
            hz.iter().map(|&f| {
                let k = TAU * f / SOUND_SPEED_M_S;
                let mut w = vec![C64::ZERO; r];
                for (s, &d) in surface.iter().zip(&distance) {
                    let (sin, cos) = (k * d).sin_cos();
                    let g = C64::new(cos, sin).scale(s.area_m2 / d);
                    for (wj, &shape) in w.iter_mut().zip(&s.mode_shape) { *wj = *wj + g.scale(shape); }
                }
                w
            }).collect()
        }).collect();
        Self { lambda, damping, bridge: modes.iter().map(|m| m.bridge).collect(), radiation }
    }

    /// Bridge displacement and mean-square microphone pressure (up to one constant)
    /// per unit force at each key.
    fn responses(&self, keys: &[usize], at: usize, hz: f64) -> Vec<(C64, f64)> {
        let r = self.lambda.len();
        let w = TAU * hz;
        let mut a = vec![C64::ZERO; r * r];
        for i in 0..r {
            for j in 0..r { a[i * r + j] = C64::new(0.0, -w * self.damping[i * r + j]); }
            a[i * r + i] = a[i * r + i] + C64::from_re(self.lambda[i] - w * w);
        }
        let lu = lu_complex(&a, r).expect("board dynamic stiffness is singular");
        keys.iter().map(|&key| {
            let mut x: Vec<C64> = self.bridge.iter().map(|b| C64::from_re(b[key])).collect();
            lu.solve(&mut x);
            let h = (0..r).fold(C64::ZERO, |h, j| h + x[j].scale(self.bridge[j][key]));
            let pressure = self.radiation.iter().map(|mic| {
                (0..r).fold(C64::ZERO, |p, j| p + x[j] * mic[at][j]).scale(-w * w).norm_sq()
            }).sum::<f64>() / self.radiation.len() as f64;
            (h, pressure)
        }).collect()
    }

    fn mass_line(&self, key: usize) -> f64 { self.bridge.iter().map(|b| b[key] * b[key]).sum() }
}

/// Per key, octave-band means of conductance and of mean-square pressure, both over linear frequency.
fn bands(board: &Board, keys: &[usize], hz: &[f64]) -> Vec<(Vec<f64>, Vec<f64>)> {
    let mut out = vec![(vec![0.0; OCTAVES_HZ.len()], vec![0.0; OCTAVES_HZ.len()]); keys.len()];
    let mut weight = vec![0.0; OCTAVES_HZ.len()];
    for (at, &f) in hz.iter().enumerate() {
        let band = at / POINTS_PER_OCTAVE;
        // Log-spaced samples: weight by f for a mean over linear frequency.
        weight[band] += f;
        for ((g, p), (h, pressure)) in out.iter_mut().zip(board.responses(keys, at, f)) {
            g[band] += f * TAU * f * h.im; // Re(-i w H)
            p[band] += f * pressure;
        }
    }
    for (g, p) in &mut out {
        for b in 0..weight.len() { g[b] /= weight[b]; p[b] /= weight[b]; }
    }
    out
}

fn db(ratio: f64) -> f64 { 10.0 * ratio.log10() }

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (mut mesh, mut band, mut keys, mut specs) = ("16,0.03".to_owned(), 2400.0, vec![33u8, 45, 57, 60, 72, 84, 96], Vec::new());
    let mut mics = vec![Config::live(48_000).microphone_m];
    let mut i = 0;
    while i < args.len() {
        let value = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--mesh" => { mesh = value(i)?; i += 1; }
            "--band" => { band = value(i)?.parse().map_err(|_| "invalid --band")?; i += 1; }
            "--keys" => {
                keys = value(i)?.split(',').map(|k| k.parse().ok().filter(|m| (21..=108).contains(m))
                    .ok_or("a key is a MIDI note 21-108")).collect::<Result<_, _>>()?;
                i += 1;
            }
            "--mics" => {
                mics = value(i)?.split(';').map(|m| {
                    let v: Vec<f64> = m.split(',').map(|x| x.parse().map_err(|_| "invalid microphone")).collect::<Result<_, _>>()?;
                    <[f64; 3]>::try_from(v).map_err(|_| "a microphone is X,Y,Z")
                }).collect::<Result<_, _>>()?;
                i += 1;
            }
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => specs.push(s.to_owned()),
        }
        i += 1;
    }
    let preset = if mesh == "preset" { steinway_d::build(8)? } else {
        let (d, s) = mesh.split_once(',').ok_or("--mesh is DIVISIONS,STEP or preset")?;
        steinway_d::build_resolved(d.parse().map_err(|_| "divisions")?, s.parse().map_err(|_| "rib-normal step")?, true)?
    };
    let courses = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV)?;
    let all: Vec<u8> = courses.iter().map(|c| c.midi).collect();
    let hz: Vec<f64> = (0..OCTAVES_HZ.len() * POINTS_PER_OCTAVE).map(|n| {
        OCTAVES_HZ[0] * 2f64.powf((n as f64 + 0.5) / POINTS_PER_OCTAVE as f64 - 0.5)
    }).collect();
    let start = std::time::Instant::now();
    let source = board_geometry::BoardGeometry::read(&preset.geometry)?.prepare_source_slice(&all, band)?;
    let reference = Board::new(&source.modes, None, &source.surface, &mics, &hz);
    println!("source slice: {} modes to {:.0} Hz, {} surface points, {} microphone(s) ({:.1} s)", source.modes.len(),
        source.modes.last().map_or(0.0, |m| m.frequency_hz), source.surface.len(), mics.len(), start.elapsed().as_secs_f64());
    // Rows are labelled #spec-number/coordinates.
    let reduced: Vec<(String, Board)> = specs.iter().enumerate().map(|(s, spec)| {
        let desc = Description::with_reduction_on(&preset, courses.clone(), band, spec)?;
        let top = desc.board.iter().map(|m| m.frequency_hz).fold(0.0, f64::max);
        println!("#{} {spec}: {} coordinates ({} exact), top {top:.0} Hz", s + 1, desc.board.len(), desc.exact_low_modes);
        Ok((format!("#{}/{}", s + 1, desc.board.len()), Board::new(&desc.board, desc.damping.as_deref(), &desc.surface, &mics, &hz)))
    }).collect::<Result<_, String>>()?;
    let header: Vec<String> = OCTAVES_HZ.iter().map(|f| format!("{:>6}", if *f < 1000.0 { format!("{f:.0}") } else { format!("{:.0}k", f / 1000.0) })).collect();
    println!("dB against the source slice; G = driving-point conductance, P = microphone pressure; mass = sum b^2");
    // Tilt: pressure in the 2k-8k octaves over the 250-1k octaves (about 1.4-11 kHz over 0.18-1.4 kHz).
    println!("key  #/N     | G {} | P {} | mass  P 2k-8k/250-1k", header.join(""), header.join(""));
    let index: Vec<usize> = keys.iter().map(|&m| usize::from(m - 21)).collect();
    let source_bands = bands(&reference, &index, &hz);
    let reduced_bands: Vec<_> = reduced.iter().map(|(_, board)| bands(board, &index, &hz)).collect();
    let tilt = |p: &[f64]| db((p[4] + p[5] + p[6]) / (p[1] + p[2] + p[3]));
    let cells = |x: &[f64], x0: &[f64]| x.iter().zip(x0).map(|(a, b)| format!("{:>6.1}", db(a / b))).collect::<String>();
    for (k, &midi) in keys.iter().enumerate() {
        let (g0, p0) = &source_bands[k];
        for ((label, board), spectra) in reduced.iter().zip(&reduced_bands) {
            let (g, p) = &spectra[k];
            println!("{midi:<4} {label:<7} | G {} | P {} | {:>5.1} {:>6.1}", cells(g, g0), cells(p, p0),
                db(board.mass_line(index[k]) / reference.mass_line(index[k])), tilt(p) - tilt(p0));
        }
    }
    Ok(())
}
