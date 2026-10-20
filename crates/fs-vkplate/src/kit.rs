//! A standard drum kit built from one solver: every head is a tensioned Mylar film (membrane
//! tension plus the film's own bending stiffness) and every cymbal a free-edged von Karman shallow
//! shell. All numbers below are declared ESTIMATES of ordinary instruments (dimensions from common
//! sizes, film and bronze from handbook values), not measurements of a particular specimen.

use crate::*;

pub const OUT_RATE: f64 = 48_000.0;

/// Effective tip mass of a hickory stick in a short impact: FrankenSim's own Z5A stick estimate
/// (crates/fs-couple/examples/percussion: contact_effective_mass_kg = 0.00724314).
pub const STICK_TIP_KG: f64 = 0.00724314;

/// Polyester (Mylar) film.
const MYLAR: (f64, f64, f64) = (4.0e9, 0.38, 1390.0);
/// B20 cymbal bronze.
pub(crate) const BRONZE: (f64, f64, f64) = (1.05e11, 0.33, 8600.0);

pub struct Piece {
    pub instrument: Instrument,
    pub oversample: usize,
    /// Seconds of simulation to discard at the start (snare-wire preload settling, hi-hat closing).
    pub settle: f64,
    pub notes: Vec<String>,
}

#[derive(Clone, Copy)]
pub struct Head {
    pub thickness: f64,
    pub tension: f64,
    pub t60: [(f64, f64); 2],
}

/// Snare wires: coiled steel strands under the resonant head, pressed up into it by the strainer.
#[derive(Clone, Copy, Debug)]
pub struct Wires {
    pub strands: i32,
    pub spacing: f64,
    pub length: f64,
    pub mass_per_m: f64,
    pub tension: f64,
    pub preload_per_m: f64,
    pub t60: f64,
    pub contact: Law,
}

impl Wires {
    pub fn estimated() -> Wires {
        // a loosely set strainer: slack strands slap the head instead of tracking it (measured: the most buzz)
        Wires { strands: 18, spacing: 0.004, length: 0.36, mass_per_m: 0.004, tension: 0.8, preload_per_m: 0.22, t60: 0.35,
                contact: Law::Power { k: 3.0e7, alpha: 1.5, mu: 0.35 } }
    }
}

pub struct DrumSpec {
    pub name: &'static str,
    pub radius: f64,
    pub depth: f64,
    pub batter: Head,
    pub reso: Head,
    pub grid: f64,
    pub oversample: usize,
    pub strike_xy: (f64, f64),
    /// (mass kg, law, footprint radius m)
    pub striker: (f64, Law, f64),
    pub muffle: Option<(f64, f64, f64, f64)>,
    pub wires: Option<Wires>,
    /// Front-head port radius (kick drums), m.
    pub port: Option<f64>,
    pub mic: [f64; 3],
    /// Floor height relative to the batter head (None: free field).
    pub floor: Option<f64>,
}

pub(crate) fn hertz_stick(tip_radius: f64, e_target: f64, nu_target: f64) -> Law {
    // hickory tip on the target: K = 4/3 E* sqrt(R)
    let (e_w, nu_w) = (12.0e9, 0.3);
    let es = 1.0 / ((1.0 - nu_w * nu_w) / e_w + (1.0 - nu_target * nu_target) / e_target);
    Law::Power { k: 4.0 / 3.0 * es * tip_radius.sqrt(), alpha: 1.5, mu: 0.05 }
}

fn head_body(name: &str, r: f64, grid: f64, hd: Head, muffle: Option<(f64, f64, f64, f64)>) -> Result<Body, String> {
    let t = hd.thickness;
    let th = move |_r: f64| t;
    let z = |_r: f64| 0.0;
    Body::plate(name, &PlateSpec {
        radius: r, hole: 0.0, grid, edge: Edge::Clamped, thickness: &th, rest_height: &z,
        young: MYLAR.0, poisson: MYLAR.1, density: MYLAR.2, tension: hd.tension, bending: true,
        nonlinear: Nonlinear::Berger, t60: hd.t60, support: None, muffle,
    })
}

/// Lowest (0,1) membrane frequency of a head, ignoring air (for the radiation resistance).
fn f01(r: f64, hd: Head) -> f64 {
    2.405 * (hd.tension / (MYLAR.2 * hd.thickness)).sqrt() / (2.0 * PI * r)
}

pub fn drum(s: &DrumSpec, hits: &[(f64, f64)]) -> Result<Piece, String> {
    let rate = OUT_RATE * s.oversample as f64;
    let mut ins = Instrument::new(rate);
    let mut notes = vec![];
    let mut batter = head_body(&format!("{} batter", s.name), s.radius, s.grid, s.batter, s.muffle)?;
    batter.origin = [0.0, 0.0, 0.0];
    batter.facing = 1.0;
    let mut reso = head_body(&format!("{} resonant", s.name), s.radius, s.grid, s.reso, s.muffle)?;
    reso.origin = [0.0, 0.0, -s.depth];
    reso.facing = -1.0;
    notes.push(format!("heads: {} nodes each at {:.1} mm; (0,1) in vacuo {:.0} / {:.0} Hz", batter.n, s.grid * 1e3, f01(s.radius, s.batter), f01(s.radius, s.reso)));
    // outside air on each head's volume motion: a baffled piston's radiation mass, 8 rho a^3 / 3
    batter.air_mass = 8.0 * RHO_AIR * s.radius.powi(3) / 3.0;
    reso.air_mass = batter.air_mass;
    let bi = ins.add(batter);
    let ri = ins.add(reso);
    let port = s.port.map(|pr| {
        let mut q = Port::new(pr, s.reso.thickness, [0.0, -0.6 * s.radius, -s.depth], -1.0);
        // small-ka radiation resistance at the Helmholtz frequency of cavity plus port
        let v = PI * s.radius * s.radius * s.depth;
        let fh = C_AIR / (2.0 * PI) * (q.area / (v * q.length)).sqrt();
        q.resistance = RHO_AIR * (2.0 * PI * fh).powi(2) / (2.0 * PI * C_AIR);
        notes.push(format!("port {:.0} mm: Helmholtz {fh:.0} Hz with the cavity alone", 2e3 * pr));
        q
    });
    ins.cavities.push(Cavity { heads: vec![bi, ri], volume: PI * s.radius * s.radius * s.depth, port });
    ins.rads.push(RadDamp::monopole(vec![bi, ri], f01(s.radius, s.batter)));
    let mut settle = 0.0;
    if let Some(wr) = s.wires {
        let (mass_per_m, tension, preload_per_m, len) = (wr.mass_per_m, wr.tension, wr.preload_per_m, wr.length);
        let h = s.grid;
        let n = (len / h).round() as usize;
        let x0 = -(n as f64 - 1.0) * h / 2.0;
        let mut count = 0;
        let mut rows: Vec<i64> = vec![];
        for k in 0..wr.strands {
            let y = (k as f64 - (wr.strands as f64 - 1.0) / 2.0) * wr.spacing;
            let row = (y / h).round() as i64;
            if rows.contains(&row) {
                continue; // two strands on one grid row would share contact nodes
            }
            rows.push(row);
            let y = row as f64 * h;
            let mut wb = Body::string(&format!("wire {row}"), n, h, mass_per_m, tension, wr.t60, x0, y);
            wb.origin = [0.0, 0.0, -s.depth];
            for k in 0..n {
                wb.bias[k] = -preload_per_m * h; // toward the head (outward is +, the head is above)
            }
            let wi = ins.add(wb);
            for k in 0..n {
                let x = x0 + k as f64 * h;
                if let Some(hn) = ins.bodies[ri].node_at(x, y) {
                    let mut c = Contact::new(vec![(ri, hn, 1.0)], vec![(wi, k, 1.0)], 0.0, wr.contact);
                    c.fixed_point_exit = true;
                    ins.contacts.push(c);
                }
            }
            count += 1;
        }
        notes.push(format!("snare: {count} wires of {n} nodes, {:.0} contact points", ins.contacts.len()));
        settle = 0.25;
    }
    let (m, law, foot) = s.striker;
    let st = ins.add(Body::point("striker", m));
    let fp = ins.bodies[bi].footprint(s.strike_xy.0, s.strike_xy.1, foot);
    notes.push(format!("striker {m} kg over {} nodes", fp.len()));
    let mut c = Contact::new(fp.into_iter().map(|(k, w)| (bi, k, w)).collect(), vec![(st, 0, 1.0)], 0.0, law);
    c.active = false;
    ins.contacts.push(c);
    let ci = ins.contacts.len() - 1;
    for &(t, v) in hits {
        ins.strikes.push(Strike { time: t + settle, striker: st, contact: ci, speed: v, clear: 2e-4, done: false });
    }
    let mut mic = Mic::new(s.mic);
    mic.floor = s.floor;
    ins.mics.push(mic);
    Ok(Piece { instrument: ins, oversample: s.oversample, settle, notes })
}

// Sizes and heads follow the reference kit we grade against (Salamander Drumkit, A. Holm, CC BY-SA
// 3.0): 18x18 kick with a 2-ply batter and a coated single-ply front head, 14x6 snare, 12x7 rack tom,
// 14x14 floor tom. Tensions are the drummer's tuning: chosen so each fundamental lands near the
// recording's. Loss points come from the recordings' own band decay times (see README).

const PLY1: f64 = 2.54e-4; // single-ply 10 mil
const PLY2: f64 = 3.56e-4; // two plies of 7 mil
const SNARE_SIDE: f64 = 7.62e-5; // 3 mil
const RESO: f64 = 1.9e-4; // 7.5 mil clear resonant head

pub fn snare() -> DrumSpec {
    DrumSpec {
        name: "snare 14x6",
        radius: 0.1703,
        depth: 0.1524,
        // recording: ~0.55 s T60 across the mid bands, with off-centre strike exciting (1,1) at ~235 Hz and wire rattle
        batter: Head { thickness: PLY1, tension: 1450.0, t60: [(150.0, 0.50), (3000.0, 0.45)] },
        reso: Head { thickness: SNARE_SIDE, tension: 2250.0, t60: [(150.0, 0.50), (3000.0, 0.45)] },
        grid: 3.5e-3,
        oversample: 2,
        strike_xy: (0.072, 0.025),
        striker: (0.0060, hertz_stick(0.005, MYLAR.0, MYLAR.1), 0.0),
        muffle: None,
        wires: Some(Wires::estimated()),
        port: None,
        mic: [0.35, 0.25, 0.9], // an overhead about 1 m away, 25 degrees off the drum's axis
        floor: Some(-0.65),
    }
}

pub fn snare_wires_off() -> DrumSpec {
    DrumSpec { wires: None, name: "snare 14x6, wires off", ..snare() }
}

pub fn kick() -> DrumSpec {
    DrumSpec {
        name: "kick 18x18",
        radius: 0.221,
        depth: 0.4572,
        // a 2-ply batter and a pillow: fundamental on the recording's 44.5 Hz, the beater crown's
        // click decaying within ~60-100 ms
        batter: Head { thickness: PLY2, tension: 358.0, t60: [(45.0, 0.48), (4000.0, 0.18)] },
        reso: Head { thickness: PLY1, tension: 253.0, t60: [(45.0, 0.48), (4000.0, 0.18)] },
        grid: 4.5e-3,
        oversample: 1,
        strike_xy: (0.055, 0.015),
        // convex felt beater crown: point contact at initial impact launches the 0.5-4 kHz beater click
        striker: (0.032, Law::Power { k: 5.5e9, alpha: 2.25, mu: 0.04 }, 0.0),
        // a pillow resting against the lower part of the heads
        muffle: Some((0.0, -0.16, 0.07, 22.0)),
        wires: None,
        port: None,
        mic: [0.35, 0.25, 0.9],
        floor: Some(-0.65),
    }
}

pub fn tom(size_in: f64) -> DrumSpec {
    let (radius, depth, batter, reso, strike_xy, grid, name) = if size_in <= 12.5 {
        (0.146, 0.1778,
         Head { thickness: PLY1, tension: 1300.0, t60: [(139.0, 1.55), (2500.0, 0.25)] },
         Head { thickness: RESO, tension: 975.0, t60: [(139.0, 1.55), (2500.0, 0.25)] },
         (0.275 * 0.146, 0.0),
         4.0e-3,
         "rack tom 12x7")
    } else {
        (0.1703, 0.3556,
         Head { thickness: PLY2, tension: 410.0, t60: [(69.0, 1.5), (3500.0, 0.28)] },
         Head { thickness: RESO, tension: 215.0, t60: [(69.0, 1.5), (3500.0, 0.28)] },
         (0.43 * 0.1703, 0.0),
         4.2e-3,
         "floor tom 14x14")
    };
    DrumSpec {
        name,
        radius,
        depth,
        batter,
        reso,
        grid,
        oversample: 1,
        strike_xy,
        striker: (STICK_TIP_KG, hertz_stick(0.005, MYLAR.0, MYLAR.1), 0.0),
        muffle: None,
        wires: None,
        port: None,
        mic: [0.35, 0.25, 0.9],
        floor: Some(if size_in <= 12.5 { -0.75 } else { -0.55 }),
    }
}

/// Scale both heads' tension (a drummer's tuning turn), keeping everything else.
pub fn retune(mut s: DrumSpec, factor: f64) -> DrumSpec {
    s.batter.tension *= factor;
    s.reso.tension *= factor;
    s
}

pub struct CymbalSpec {
    pub name: &'static str,
    pub radius: f64,
    pub bow_thickness: f64,
    pub bell_thickness: f64,
    pub bell_radius: f64,
    pub bow_rise: f64,
    pub bell_rise: f64,
    pub grid: f64,
    pub t60: [(f64, f64); 2],
}

/// Grid points (spacing h) with r in [r0, r1]. A grid too coarse to put four of them in the
/// annulus (a real-time grid) widens r1 to the first ring that does; the caller then spreads the
/// annulus's total stiffness and damping over those nodes. Offline grids never widen.
fn ring_nodes(h: f64, r0: f64, r1: f64) -> (f64, usize) {
    let count = |r1: f64| {
        let m = (r1 / h).ceil() as i64 + 1;
        let mut n = 0;
        for j in -m..=m {
            for i in -m..=m {
                let r = h * ((i * i + j * j) as f64).sqrt();
                if r >= r0 && r <= r1 + 1e-12 {
                    n += 1;
                }
            }
        }
        n
    };
    let mut r = r1;
    loop {
        let n = count(r);
        if n >= 4 {
            return (r, n);
        }
        r += 0.25 * h;
    }
}

fn cymbal_body(name: &str, c: &CymbalSpec) -> Result<Body, String> {
    cymbal_shell(name, c, true)
}

/// A cymbal's thickness at radius `x` outside the bell (the bow law of `cymbal_shell`).
fn bow_thickness(c: &CymbalSpec, x: f64) -> f64 {
    c.bow_thickness * (1.15 - 0.15 * x / c.radius)
}

/// A cymbal shell; `felt` puts it on a stand's felt washers (a hi-hat's top shell hangs from the
/// clutch on the pull rod instead, so all the pedal force reaches the rims).
fn cymbal_shell(name: &str, c: &CymbalSpec, felt: bool) -> Result<Body, String> {
    let (rb, tb, tw, r) = (c.bell_radius, c.bell_thickness, c.bow_thickness, c.radius);
    // thickness: bell, a blend over 2 cm, then the bow thinning slightly toward the rim
    let th = move |x: f64| {
        let bow = tw * (1.15 - 0.15 * x / r);
        if x <= rb { tb } else if x < rb + 0.02 { let u = (x - rb) / 0.02; tb + (bow - tb) * (3.0 * u * u - 2.0 * u * u * u) } else { bow }
    };
    let (ab, al) = (c.bow_rise, c.bell_rise);
    // rest height above the rim: a slightly convex bow plus a raised bell (cosine cap)
    let z = move |x: f64| {
        let u = (1.0 - x / r).max(0.0);
        let bow = ab * (0.7 * u + 0.3 * u * u);
        let bell = if x < rb { al * 0.5 * (1.0 + (PI * x / rb).cos()) } else { 0.0 };
        bow + bell
    };
    // felt washers clamp the annulus round the 13 mm centre hole
    let (r1, n) = ring_nodes(c.grid, 0.0065, 0.02);
    let mount_scale = if r < 0.15 { 0.25 } else { 1.0 };
    let (ks, cs) = if r1 > 0.02 {
        let per = PI * (0.02f64.powi(2) - 0.0065f64.powi(2)) / (n as f64 * c.grid * c.grid);
        (9.0e7 * per * mount_scale, 2.7e3 * per * mount_scale)
    } else {
        (9.0e7 * mount_scale, 2.7e3 * mount_scale)
    };
    Body::plate(name, &PlateSpec {
        radius: r, hole: 0.0065, grid: c.grid, edge: Edge::Free, thickness: &th, rest_height: &z,
        young: BRONZE.0, poisson: BRONZE.1, density: BRONZE.2, tension: 0.0, bending: true,
        nonlinear: Nonlinear::VonKarman, t60: c.t60,
        // felt washers: ~1e5 N/m and ~3 N s/m in all, spread over the clamped annulus
        support: if felt { Some((0.0065, r1, ks, cs)) } else { None }, muffle: None,
    })
}

pub fn splash() -> CymbalSpec {
    CymbalSpec { name: "splash 8", radius: 0.1016, bow_thickness: 0.68e-3, bell_thickness: 1.5e-3, bell_radius: 0.039, bow_rise: 0.008, bell_rise: 0.014, grid: 4.5e-3, t60: [(500.0, 6.8), (10000.0, 3.1)] }
}

pub fn crash() -> CymbalSpec {
    CymbalSpec { name: "crash 18", radius: 0.2286, bow_thickness: 0.95e-3, bell_thickness: 1.6e-3, bell_radius: 0.075, bow_rise: 0.025, bell_rise: 0.015, grid: 6.0e-3, t60: [(500.0, 10.0), (10000.0, 2.0)] }
}

pub fn ride() -> CymbalSpec {
    CymbalSpec { name: "ride 20", radius: 0.254, bow_thickness: 1.65e-3, bell_thickness: 2.25e-3, bell_radius: 0.085, bow_rise: 0.022, bell_rise: 0.020, grid: 6.2e-3, t60: [(500.0, 14.0), (10000.0, 3.2)] }
}

fn pick_oversample(ins: &mut Instrument, max_ratio: f64) -> Result<usize, String> {
    let (lam, who) = ins.admission_rate2();
    for os in [1usize, 2, 3, 4, 6, 8] {
        let k = 1.0 / (OUT_RATE * os as f64);
        if k * k * lam / 4.0 <= max_ratio {
            ins.rate = OUT_RATE * os as f64;
            return Ok(os);
        }
    }
    Err(format!("{who}: leapfrog unstable even at 8x oversampling"))
}

/// A suspended cymbal struck once per hit: (time s, speed m/s, radius fraction, stick tip radius m).
pub fn cymbal(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)]) -> Result<Piece, String> {
    let mut ins = Instrument::new(OUT_RATE);
    let mut body = cymbal_body(c.name, c)?;
    body.origin = [0.0, 0.0, 0.0];
    let n = body.n;
    let ci = ins.add(body);
    let st = ins.add(Body::point("stick", STICK_TIP_KG));
    let mut notes = vec![format!("{}: {n} nodes at {:.1} mm", c.name, c.grid * 1e3)];
    for (k, &(t, v, rf, tip)) in hits.iter().enumerate() {
        let fp = ins.bodies[ci].footprint(rf * c.radius, 0.0, 0.0);
        // the stick comes down onto the upper surface: plate on side a
        let mut con = Contact::new(fp.into_iter().map(|(n_, w)| (ci, n_, w)).collect(), vec![(st, 0, 1.0)], 0.0, hertz_stick(tip, BRONZE.0, BRONZE.1));
        con.active = false;
        ins.contacts.push(con);
        ins.strikes.push(Strike { time: t, striker: st, contact: ins.contacts.len() - 1, speed: v, clear: 2e-4, done: false });
        let _ = k;
    }
    let mut mic = Mic::new([0.3, 0.2, 1.0]);
    mic.unbaffled = Some(c.radius);
    ins.mics.push(mic);
    let os = pick_oversample(&mut ins, 0.6)?;
    let (s0, s1) = ins.bodies[ci].loss();
    notes.push(format!("internal rate {} kHz; loss s0 {s0:.3}/s s1 {s1:.3e} m^2/s", OUT_RATE * os as f64 / 1000.0));
    Ok(Piece { instrument: ins, oversample: os, settle: 0.0, notes })
}

pub fn hihat_spec() -> (CymbalSpec, CymbalSpec) {
    let top = CymbalSpec { name: "hi-hat top 14", radius: 0.1778, bow_thickness: 0.85e-3, bell_thickness: 1.4e-3, bell_radius: 0.055, bow_rise: 0.015, bell_rise: 0.01, grid: 5.5e-3, t60: [(500.0, 8.0), (10000.0, 1.5)] };
    let bottom = CymbalSpec { name: "hi-hat bottom 14", bow_thickness: 1.1e-3, ..top };
    (top, bottom)
}

/// Hi-hat: opposing shells on one rod. `gap` is the undeformed face-to-face rim opening [m].
/// Bow clearance follows the supplied profiles; pedal loading deforms them from this geometry.
pub fn hihat(gap: f64, pedal_n: f64, hits: &[(f64, f64)]) -> Result<Piece, String> {
    hihat_scaled(gap, pedal_n, hits, 1.0)
}

/// The same physical hi-hat and rim loss on an explicitly coarsened grid.
pub fn hihat_scaled(gap: f64, pedal_n: f64, hits: &[(f64, f64)], scale: f64) -> Result<Piece, String> {
    // bronze rims pressed together: stiff and lightly damped so rim chatter radiates 4-16 kHz sizzle
    hihat_grid(gap, pedal_n, hits, Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.08 }, 0.8, scale)
}

/// The same with an explicit rim contact law and the inner edge of the contact annulus (fraction of the radius).
pub fn hihat_with(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64) -> Result<Piece, String> {
    hihat_grid(gap, pedal_n, hits, rim, from, 1.0)
}

/// Dry friction coefficient between the hi-hat's bronze rims. Unlubricated copper alloys on themselves
/// span about 0.3-1; with 0.3 and 0.5 the closed hat's octave-band decay grades inside the real-vs-real
/// ruler (0.40 and 0.32 octaves from the Salamander stroke, two real strokes 0.47; without friction 0.95),
/// and 0.5 is the closer of the two. The open hat, whose rims seldom touch, barely changes.
pub const RIM_FRICTION: f64 = 0.5;

/// The same with grid spacing multiplied by `scale`.
/// The reference-cell rim force law scales with area to preserve traction under refinement.
pub fn hihat_grid(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64) -> Result<Piece, String> {
    hihat_grid_mu(gap, pedal_n, hits, rim, from, scale, RIM_FRICTION)
}

/// The same with an explicit rim friction coefficient (0: the rims press but do not rub).
pub fn hihat_grid_mu(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64, friction: f64) -> Result<Piece, String> {
    if !gap.is_finite() || gap < 0.0 {
        return Err(format!("hi-hat: rim opening must be finite and nonnegative, got {gap} m"));
    }
    let (mut ts, mut bs) = hihat_spec();
    let Law::Power { k, alpha, mu } = rim;
    let rim = Law::Power { k: k * scale * scale, alpha, mu };
    ts.grid *= scale;
    bs.grid *= scale;
    // Displacements remain upward; the bottom's bell and bow face away from the top.
    bs.bow_rise = -bs.bow_rise;
    bs.bell_rise = -bs.bell_rise;
    let mut ins = Instrument::new(OUT_RATE);
    let mut top = cymbal_shell(ts.name, &ts, false)?;
    top.origin = [0.0, 0.0, 0.0];
    let mut bot = cymbal_body(bs.name, &bs)?;
    // Profiles describe midplanes, while the requested opening is between their faces.
    let rim_half_thickness = 0.5 * (bow_thickness(&ts, ts.radius) + bow_thickness(&bs, bs.radius));
    bot.origin = [0.0, 0.0, -rim_half_thickness - gap];
    let ti = ins.add(top);
    let bi = ins.add(bot);
    // the top hangs from the clutch on the pull rod; the bottom sits on the stand's felt cup
    let rod = ins.add(Body::point("rod", 0.3));
    ins.links.push(Link { a: (rod, 0), b: None, k: 2000.0, c: 20.0 });
    // Preserve the physical totals on both coarse and refined discretizations.
    let (r1, _) = ring_nodes(ins.bodies[ti].h, 0.0065, 0.02);
    let clutch: Vec<usize> = (0..ins.bodies[ti].n).filter(|&k| {
        let p = ins.bodies[ti].xy[k];
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        r >= 0.0065 && r <= r1 + 1e-12
    }).collect();
    if clutch.is_empty() { return Err("hi-hat: grid leaves no clutch nodes".into()); }
    let (kc, cc) = (8.0e6 / clutch.len() as f64, 80.0 / clutch.len() as f64);
    for &k in &clutch {
        ins.links.push(Link { a: (ti, k), b: Some((rod, 0)), k: kc, c: cc });
    }
    // rim contact between the shells: matching nodes in the outer 20 % of the radius, with slight
    // 2-lobe rim waviness (~0.15 mm) so the pedal rests on high spots while adjacent rim arcs sizzle
    let mut pairs = 0;
    for k in 0..ins.bodies[ti].n {
        let p = ins.bodies[ti].xy[k];
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        if r < from * ts.radius {
            continue;
        }
        if let Some(kb) = ins.bodies[bi].node_at(p[0], p[1]) {
            let cos2 = if r > 0.0 { (p[0] * p[0] - p[1] * p[1]) / (r * r) } else { 1.0 };
            let (tb, tt) = (bow_thickness(&bs, r), bow_thickness(&ts, r));
            let top_z = ins.bodies[ti].origin[2] + ins.bodies[ti].rest_z[k];
            let bottom_z = ins.bodies[bi].origin[2] + ins.bodies[bi].rest_z[kb];
            let local_gap = top_z - bottom_z - 0.5 * (tt + tb) + 7.5e-5 * (1.0 - cos2);
            ins.contacts.push(Contact::new(vec![(bi, kb, 1.0)], vec![(ti, k, 1.0)], local_gap, rim));
            if friction > 0.0 {
                // the faces the contact presses together rub as the shells bend (see Friction)
                let (tb, tt) = (bow_thickness(&bs, r), bow_thickness(&ts, r));
                let fr = Friction::between(&ins.bodies, ins.contacts.len() - 1, bi, kb, tb, ti, k, tt, friction);
                ins.frictions.push(fr);
            }
            pairs += 1;
        }
    }
    // the pedal pulls the hanging top shell ~1 cm down onto the bottom; let the suspension settle
    let settle = if pedal_n != 0.0 { 0.5 } else { 0.0 };
    if pedal_n != 0.0 {
        ins.drives.push(Drive { body: rod, node: 0, knots: vec![(0.0, 0.0), (0.12, -pedal_n), (1e9, -pedal_n)] });
    }
    let st = ins.add(Body::point("stick", STICK_TIP_KG));
    for &(t, v) in hits {
        let fp = ins.bodies[ti].footprint(0.85 * ts.radius, 0.0, 0.0);
        let tip = 0.004;
        let mut con = Contact::new(fp.into_iter().map(|(n_, w)| (ti, n_, w)).collect(), vec![(st, 0, 1.0)], 0.0, hertz_stick(tip, BRONZE.0, BRONZE.1));
        con.active = false;
        ins.contacts.push(con);
        ins.strikes.push(Strike { time: t + settle, striker: st, contact: ins.contacts.len() - 1, speed: v, clear: 2e-4, done: false });
    }
    let mut mic = Mic::new([0.3, 0.2, 1.0]);
    mic.unbaffled = Some(ts.radius);
    ins.mics.push(mic);
    let os = pick_oversample(&mut ins, 0.6)?;
    let notes = vec![format!("hi-hat: {} + {} nodes, {pairs} rim contact pairs (friction {friction}), rim opening {:.2} mm, pedal {pedal_n} N, internal {} kHz", ins.bodies[ti].n, ins.bodies[bi].n, gap * 1e3, OUT_RATE * os as f64 / 1000.0)];
    Ok(Piece { instrument: ins, oversample: os, settle, notes })
}

/// Run a piece for `seconds` of output after its settling time; returns 48 kHz pressure [Pa].
pub fn render(p: &mut Piece, seconds: f64, mut progress: impl FnMut(f64)) -> Vec<Vec<f64>> {
    let ins = &mut p.instrument;
    let nm = ins.mics.len();
    let total = ((seconds + p.settle) * ins.rate).round() as usize;
    let skip = (p.settle * ins.rate).round() as usize;
    let mut raw = vec![Vec::with_capacity(total - skip); nm];
    let mut buf = vec![0.0; nm];
    for s in 0..total {
        ins.step(&mut buf);
        if s >= skip {
            for (c, v) in raw.iter_mut().zip(&buf) {
                c.push(*v);
            }
        }
        if s % 4800 == 0 {
            progress(s as f64 / total as f64);
        }
    }
    raw.iter().map(|c| decimate(c, p.oversample)).collect()
}
