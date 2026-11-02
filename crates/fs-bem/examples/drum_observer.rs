//! Bake `crates/fs-vkplate/data/closed-body/<piece>.tsv` (pexfp.378): what a
//! kit drum's overhead microphone hears from the drum as ONE rigid closed body
//! (batter, resonant head and a rigid shell) standing on a rigid floor, for
//! each retained head shape, below the closed body's CBIE band top f_x.
//!
//! Geometry: a rigid cylinder of the card's radius and depth, the batter cap
//! at z = 0 (outward +z) and the resonant cap at z = -depth (outward -z), as
//! fs-vkplate's `kit::build_drum` places the heads. The floor is the exact
//! image plane: the drum and its mirror image are two disjoint closed bodies
//! whose heads move as mirror images, observed at the microphone alone.
//!
//! Sources: every clamped-membrane shape `J_m(j_mn r/a) cos|sin(m theta)`
//! with `j_mn <= ZERO_LIMIT` on each head, normalized to unit mean square as
//! fs-vkplate's `radiation::ModeShape` (duplicated here, since fs-vkplate is
//! a standalone workspace; the table's `probe` rows let fs-vkplate refuse any
//! difference). A source's trace is the shape's value at each head panel's
//! centroid: outward normal velocity per unit modal rate [m/s].
//!
//! Producer: `fs_bem::observation_transfers::BemObservation::new_routed`
//! with the triangle rule, used ONLY below the geometric CBIE bound (one
//! discrete operator across the band; Burton-Miller is refused, since its
//! low-ka dipole resistance bias sits on the drums' together lines). The taper
//! stops at f_x = `GeometryPolicy::plain_cbie_limit` (0.8 of the box
//! Dirichlet bound), so no bake bin reaches the other formulation; a bin that
//! would is refused, never rerouted.
//!
//! The bake: `bake_observation_transfers` with one extra analytic source
//! (`COMPLEMENT_SOURCE`, 2 i tan(omega / 2 fs)): its kernel is the bake's own
//! C4 taper and latency on the ratio of the centred second difference to the
//! centred half difference, which fs-vkplate needs for the exact complement
//! (delta_D - taper) on the baffled Rayleigh path of the retained shapes above
//! f_x.
//!
//! The table stores the producer's own pressure rows (exp(-i omega t)) at
//! every finest-rung bin below the taper stop, the declaration (rate, band,
//! latency, FFT length, integer delays) and, per kernel, the energy, energy
//! centroid and one in-band transfer of fs-bem's own FIR, against which
//! fs-vkplate checks its reconstruction at load. Self-check: every transfer is
//! computed both ways (each shape radiating to the microphone; a point source
//! at the microphone scattered by the bodies, integrated against the shape)
//! at four frequencies up to f_x; the bake refuses above `RECIPROCITY_LIMIT`.
//!
//! Run:  cargo run --release -p fs-bem --example drum_observer -- PIECE RATE_HZ OUT.tsv [--workers N]
//!       cargo run --release -p fs-bem --example drum_observer -- --verify TABLE.tsv [--workers N]
//! PIECE: kick | snare | snare-off | tom-12 | tom-14. RATE_HZ is the piece's
//! internal clock (the kit example prints it). `--verify` re-bakes from the
//! table's own `bake-args` and refuses unless the result is byte-identical.
//! `--workers` declares fs-bem's dense workers (bit-identical for every
//! count). A bake refused after its solves leaves the rows in
//! OUT.tsv.refused-rows (exact bits).
use fs_bem::helmholtz::bc::Admittance;
use fs_bem::helmholtz::bc::incident::{IncidentField, Monopole, solve_scattering};
use fs_bem::helmholtz::{Formulation, HelmholtzError, Medium, PanelRule};
use fs_bem::observation_transfers::{
    BemObservation, ObservationTransferBake, ProducerUncertainty, ReceiverSpec, SourceId,
    TransferDeclaration, TransferError, bake_observation_transfers,
};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_blake3::DomainHasher;
use fs_math::{c64::C64, det};
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::fmt::Write as _;

const TABLE_FORMAT: &str = "frankensim.vkplate-closed-body-observer.v1";
const ROUTE: &str = "closed-body-triangle-cbie";
const MODE_SET: &str = "clamped-membrane-bessel";
/// Every head shape with a Bessel zero at most this is a source. It keeps
/// (0,1)-(0,3), the lumped model's whole set behind the stamped predictions,
/// and every m <= 5 shape below j_03.
const ZERO_LIMIT: f64 = 9.0;
/// Rim vertices (caps and shell) and radial rings of each cap: about 0.11
/// head radius per panel, six panels per radial wavelength of J_m(9 r).
const SECTORS: usize = 56;
const CAP_RINGS: usize = 9;
/// Shell rings are at most this many head-panel sizes tall.
const SHELL_ASPECT: f64 = 2.0;
/// Receiver clearance admitted by the observation owner [m].
const CLEARANCE_M: f64 = 0.02;
/// The bake refuses a both-ways transfer difference above this (relative to
/// the largest transfer at that frequency): .376's reciprocal-part number.
const RECIPROCITY_LIMIT: f64 = 0.02;
/// FIR gates: the fs-bem maxima.
const TAIL_TOLERANCE: f64 = 1e-10;
const DTFT_TOLERANCE: f64 = 1e-4;
/// The highest both-ways check and the crossover step, as a fraction of f_x.
const TOP_CHECK: f64 = 0.995;
/// Highest power of i omega among the FIR sizing probes. The first bakes,
/// sized on the complement alone, refused all five drums with tails and DTFT
/// changes that (i omega)^5..7 reproduces at their FFT and latency: the
/// high-order shapes rise that steeply to f_x. One power of margin.
const STEEP_PROBE_ORDER: usize = 8;
/// The analytic source whose kernel carries the complement: 2 i tan(omega / 2 fs),
/// the centred second difference over the centred half difference, behind the
/// bake's own taper and latency. A unit transfer cannot serve: the bake zeroes
/// the DC bin, and a truncated kernel with a DC hole leaves a floor that its
/// own wrap gate refuses.
const COMPLEMENT_SOURCE: &str = "complement: 2 i tan(omega / 2 fs), zero delay (the bake's own taper and latency)";

/// The complement source's producer row (exp(-i omega t)) at `omega` [rad/s].
/// fs-vkplate's `closed_body::complement_row` evaluates the same expression.
fn complement_row(omega: f64, fs: f64) -> C64 {
    C64::new(0.0, -2.0 * det::tan(0.5 * omega / fs))
}
/// The kit's gas (fs-vkplate's `standard_gas`): dry air at 293.15 K, 101325 Pa.
const GAS_STATE: (f64, f64, f64) = (293.15, 101_325.0, 0.0);

// ------------------------------------------------------------------ cards

/// One kit piece's drum as fs-vkplate's card builds it.
#[derive(Clone, Copy, Debug)]
struct Card {
    piece: &'static str,
    radius: f64,
    depth: f64,
    mic: [f64; 3],
    /// (point on the floor, unit normal into the air).
    floor: ([f64; 3], [f64; 3]),
    /// Pass-band top as a fraction of f_x: every stamped prediction line of
    /// the piece lies in the pass band, where the closed body is heard alone.
    pass_fraction: f64,
}

fn card(piece: &str) -> Option<Card> {
    let overhead = [0.35, 0.25, 0.9];
    let z_floor = |h: f64| ([0.0, 0.0, h], [0.0, 0.0, 1.0]);
    Some(match piece {
        // kit.rs kick(): lying on its side (y up), the floor tangent to the
        // hoops (KICK_HOOP_RISE_M = 0.02 beyond the shell radius).
        "kick" => Card {
            piece: "kick",
            radius: 0.221,
            depth: 0.4572,
            mic: [0.0, 0.221 + 1.0, 0.1],
            floor: ([0.0, -(0.221 + 0.02), 0.0], [0.0, 1.0, 0.0]),
            pass_fraction: 0.4,
        },
        "snare" | "snare-off" => Card {
            piece: if piece == "snare" { "snare" } else { "snare-off" },
            radius: 0.1703,
            depth: 0.1524,
            mic: overhead,
            floor: z_floor(-0.65),
            pass_fraction: 0.6,
        },
        "tom-12" => Card {
            piece: "tom-12",
            radius: 0.146,
            depth: 0.1778,
            mic: overhead,
            floor: z_floor(-0.75),
            pass_fraction: 0.6,
        },
        "tom-14" => Card {
            piece: "tom-14",
            radius: 0.1703,
            depth: 0.3556,
            mic: overhead,
            floor: z_floor(-0.55),
            pass_fraction: 0.4,
        },
        _ => return None,
    })
}

// ------------------------------------------------------------------ shapes

/// `J_m(x)`, `x >= 0`, by Miller's backward recurrence normalized with
/// `1 = J_0 + 2 sum J_2k` (fs-vkplate's `radiation::bessel_j`).
fn bessel_j(order: usize, x: f64) -> f64 {
    if x == 0.0 {
        return if order == 0 { 1.0 } else { 0.0 };
    }
    let span = order.max(x as usize);
    let top = 2 * ((span + 24 + (40.0 * span as f64).sqrt() as usize) / 2);
    let (mut above, mut at) = (0.0, 1.0);
    let mut even = if top % 2 == 0 { at } else { 0.0 };
    let mut value = if order == top { at } else { 0.0 };
    for k in (1..=top).rev() {
        let below = 2.0 * k as f64 / x * at - above;
        above = at;
        at = below;
        let index = k - 1;
        if index == order {
            value = at;
        }
        if index >= 2 && index % 2 == 0 {
            even += at;
        }
        if at.abs() > 1e250 {
            at *= 1e-250;
            above *= 1e-250;
            even *= 1e-250;
            value *= 1e-250;
        }
    }
    value / (at + 2.0 * even)
}

fn bessel_zero(order: usize, index: usize) -> f64 {
    let mut found = 0;
    let (mut x, step) = (order as f64 + 0.5, 0.05);
    let mut previous = bessel_j(order, x);
    loop {
        let next = x + step;
        let value = bessel_j(order, next);
        if previous == 0.0 || previous.signum() != value.signum() {
            found += 1;
            if found == index {
                let (mut lo, mut hi) = (x, next);
                let lo_sign = previous.signum();
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if mid <= lo || mid >= hi {
                        break;
                    }
                    if bessel_j(order, mid).signum() == lo_sign {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return 0.5 * (lo + hi);
            }
        }
        x = next;
        previous = value;
    }
}

/// Unit mean-square shape on the unit disk (fs-vkplate's `ModeShape`).
#[derive(Clone, Debug)]
struct Shape {
    m: usize,
    n: usize,
    zero: f64,
    scale: f64,
}

impl Shape {
    fn new(m: usize, n: usize) -> Shape {
        let zero = bessel_zero(m, n);
        let angular = if m == 0 { 1.0 } else { 2.0_f64.sqrt() };
        Shape { m, n, zero, scale: angular / bessel_j(m + 1, zero).abs() }
    }
    fn value(&self, x: f64, y: f64, sine: bool) -> f64 {
        let r2 = x * x + y * y;
        if r2 >= 1.0 || (self.m == 0 && sine) {
            return 0.0;
        }
        if r2 == 0.0 {
            return if self.m == 0 { self.scale } else { 0.0 };
        }
        let r = r2.sqrt();
        let (c, s) = (x / r, y / r);
        let (mut re, mut im) = (1.0, 0.0);
        for _ in 0..self.m {
            (re, im) = (re * c - im * s, re * s + im * c);
        }
        self.scale * bessel_j(self.m, self.zero * r) * if sine { im } else { re }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Part {
    Batter,
    Resonant,
    Shell,
}

impl Part {
    fn label(self) -> &'static str {
        match self {
            Part::Batter => "batter",
            Part::Resonant => "resonant",
            Part::Shell => "shell",
        }
    }
}

/// One retained source: a shape on one head, in one parity.
#[derive(Clone, Debug)]
struct Source {
    head: Part,
    shape: Shape,
    sine: bool,
}

impl Source {
    fn id(&self) -> SourceId {
        SourceId::StructuralMode(format!(
            "{} ({},{}) {}: outward normal velocity per unit modal rate [m/s]",
            self.head.label(),
            self.shape.m,
            self.shape.n,
            if self.sine { "sin" } else { "cos" }
        ))
    }
}

fn retained_sources() -> Vec<Source> {
    let mut shapes = Vec::new();
    for m in 0.. {
        if Shape::new(m, 1).zero > ZERO_LIMIT {
            break;
        }
        for n in 1.. {
            let shape = Shape::new(m, n);
            if shape.zero > ZERO_LIMIT {
                break;
            }
            shapes.push(shape);
        }
    }
    shapes.sort_by(|a, b| a.zero.total_cmp(&b.zero));
    let mut out = Vec::new();
    for head in [Part::Batter, Part::Resonant] {
        for shape in &shapes {
            for sine in [false, true] {
                if sine && shape.m == 0 {
                    continue;
                }
                out.push(Source { head, shape: shape.clone(), sine });
            }
        }
    }
    out
}

// ------------------------------------------------------------------ meshes

/// A closed triangle surface with each panel's part and its body (0: the
/// drum, 1: its floor image).
#[derive(Clone, Debug, Default)]
struct Mesh {
    tris: Vec<[[f64; 3]; 3]>,
    part: Vec<Part>,
    body: Vec<usize>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn centroid(t: &[[f64; 3]; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (t[0][i] + t[1][i] + t[2][i]) / 3.0)
}

impl Mesh {
    /// Push a triangle oriented so its normal has a positive component on `outward`.
    fn push(&mut self, mut t: [[f64; 3]; 3], outward: [f64; 3], part: Part, body: usize) {
        if dot(cross(sub(t[1], t[0]), sub(t[2], t[0])), outward) < 0.0 {
            t.swap(1, 2);
        }
        self.tris.push(t);
        self.part.push(part);
        self.body.push(body);
    }
    fn append(&mut self, other: &Mesh) {
        self.tris.extend_from_slice(&other.tris);
        self.part.extend_from_slice(&other.part);
        self.body.extend_from_slice(&other.body);
    }
}

/// `count` points of the unit circle, the first at angle 0.
fn circle(count: usize) -> Vec<(f64, f64)> {
    (0..count)
        .map(|j| {
            let (s, c) = det::sin_cos(TAU * j as f64 / count as f64);
            (c, s)
        })
        .collect()
}

/// Triangulate the annulus between two concentric rings (inner may be one point).
fn zipper(mesh: &mut Mesh, inner: &[[f64; 3]], outer: &[[f64; 3]], outward: [f64; 3], part: Part) {
    let (p, q) = (inner.len(), outer.len());
    if p == 1 {
        for o in 0..q {
            mesh.push([inner[0], outer[o], outer[(o + 1) % q]], outward, part, 0);
        }
        return;
    }
    let (mut i, mut o) = (0, 0);
    while i < p || o < q {
        // advance the ring whose next vertex comes first in angle: (o+1)/q <= (i+1)/p
        if o < q && (i >= p || (o + 1) * p <= (i + 1) * q) {
            mesh.push([inner[i % p], outer[o % q], outer[(o + 1) % q]], outward, part, 0);
            o += 1;
        } else {
            mesh.push([inner[i % p], outer[o % q], inner[(i + 1) % p]], outward, part, 0);
            i += 1;
        }
    }
}

/// The rigid closed drum: two polar-ring caps welded to a ringed shell at
/// exactly shared rim vertices.
fn drum_body(radius: f64, depth: f64, sectors: usize, cap_rings: usize) -> Mesh {
    let rim = circle(sectors);
    let mut mesh = Mesh::default();
    for (z, sign, part) in [(0.0, 1.0, Part::Batter), (-depth, -1.0, Part::Resonant)] {
        let outward = [0.0, 0.0, sign];
        let mut inner = vec![[0.0, 0.0, z]];
        for i in 1..=cap_rings {
            let ring: Vec<[f64; 3]> = if i == cap_rings {
                rim.iter().map(|&(c, s)| [radius * c, radius * s, z]).collect()
            } else {
                let r = radius * i as f64 / cap_rings as f64;
                let count = ((sectors * i) as f64 / cap_rings as f64).round().max(3.0) as usize;
                circle(count).into_iter().map(|(c, s)| [r * c, r * s, z]).collect()
            };
            zipper(&mut mesh, &inner, &ring, outward, part);
            inner = ring;
        }
    }
    let panel = TAU * radius / sectors as f64;
    let rings = (depth / (SHELL_ASPECT * panel)).ceil().max(1.0) as usize;
    let height = |k: usize| {
        if k == 0 {
            0.0
        } else if k == rings {
            -depth
        } else {
            -depth * k as f64 / rings as f64
        }
    };
    for k in 0..rings {
        let (z0, z1) = (height(k), height(k + 1));
        for j in 0..sectors {
            let (c0, s0) = rim[j];
            let (c1, s1) = rim[(j + 1) % sectors];
            let a = [radius * c0, radius * s0, z0];
            let b = [radius * c1, radius * s1, z0];
            let c = [radius * c1, radius * s1, z1];
            let d = [radius * c0, radius * s0, z1];
            let (sm, cm) = det::sin_cos(TAU * (j as f64 + 0.5) / sectors as f64);
            let outward = [cm, sm, 0.0];
            mesh.push([a, b, c], outward, Part::Shell, 0);
            mesh.push([a, c, d], outward, Part::Shell, 0);
        }
    }
    mesh
}

/// Mirror image across the floor plane, winding reversed (still outward).
fn mirror_point(p: [f64; 3], point: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let d = dot(sub(p, point), normal);
    std::array::from_fn(|i| p[i] - 2.0 * d * normal[i])
}

fn mirrored(mesh: &Mesh, point: [f64; 3], normal: [f64; 3]) -> Mesh {
    let mut out = mesh.clone();
    for t in out.tris.iter_mut() {
        let m = t.map(|p| mirror_point(p, point, normal));
        *t = [m[0], m[2], m[1]];
    }
    out.body.iter_mut().for_each(|b| *b = 1);
    out
}

/// The drum and (with a floor) its image body.
fn scene(card: &Card, with_image: bool) -> Mesh {
    let one = drum_body(card.radius, card.depth, SECTORS, CAP_RINGS);
    let mut all = one.clone();
    if with_image {
        all.append(&mirrored(&one, card.floor.0, card.floor.1));
    }
    all
}

/// A source's trace on every panel (the image body repeats the drum's).
fn trace(mesh: &Mesh, radius: f64, f: &dyn Fn(Part, f64, f64) -> f64) -> Vec<f64> {
    mesh.tris
        .iter()
        .zip(&mesh.part)
        .zip(&mesh.body)
        .enumerate()
        .map(|(i, ((t, &part), &body))| {
            if part == Part::Shell {
                return 0.0;
            }
            // the image body's panel i - half is the mirror of the drum's
            let source = if body == 0 { *t } else { mesh.tris[i - mesh.tris.len() / 2] };
            let c = centroid(&source);
            f(part, c[0] / radius, c[1] / radius)
        })
        .collect()
}

fn source_trace(mesh: &Mesh, radius: f64, s: &Source) -> Vec<f64> {
    trace(mesh, radius, &|part, x, y| if part == s.head { s.shape.value(x, y, s.sine) } else { 0.0 })
}

fn medium(gas: (f64, f64)) -> Medium {
    Medium { density: gas.0, sound_speed: gas.1 }
}

fn receiver(mic: [f64; 3]) -> ReceiverSpec {
    ReceiverSpec { id: "overhead pressure+Cartesian-gradient [Pa per unit modal rate]".into(), position_m: mic }
}

/// Pressure [source][receiver] for each field (p rows only), routed CBIE,
/// refused above the CBIE band.
fn forward(
    surface: &SpherePanels,
    gas: (f64, f64),
    ids: &[SourceId],
    traces: &[Vec<f64>],
    receivers: &[[f64; 3]],
    omega: f64,
) -> Result<Vec<Vec<C64>>, String> {
    let n = surface.areas().len();
    let specs: Vec<ReceiverSpec> = receivers
        .iter()
        .enumerate()
        .map(|(i, &p)| ReceiverSpec { id: format!("receiver {i} pressure+Cartesian-gradient"), position_m: p })
        .collect();
    let obs = BemObservation::new_routed(
        surface,
        medium(gas),
        ids,
        traces,
        &specs,
        CLEARANCE_M,
        PanelRule::Triangle { max_pairs: n * n },
    )
    .map_err(|e| format!("observation admission: {e}"))?;
    let policy = GeometryPolicy::new(surface).map_err(|e| e.to_string())?;
    if policy.formulation(omega / gas.1).map_err(|e| e.to_string())? != Formulation::PlainCbie {
        return Err(format!("{:.1} Hz is above the CBIE band; refused", omega / TAU));
    }
    let row = obs.frequency_sample(omega).map_err(|e| e.to_string())?;
    Ok(row.chunks(4 * receivers.len()).map(|s| s.chunks(4).map(|c| c[0]).collect()).collect())
}

/// Baffled Rayleigh hearing (today's fs-vkplate `Mic`, on panels): each head
/// panel facing a target radiates 2G, the others are skipped.
fn baffled(mesh: &Mesh, surface: &SpherePanels, gas: (f64, f64), field: &[f64], targets: &[[f64; 3]], omega: f64) -> C64 {
    let k = omega / gas.1;
    let mut total = C64::ZERO;
    for i in 0..mesh.tris.len() {
        if field[i] == 0.0 || mesh.body[i] != 0 {
            continue;
        }
        let (c, n, area) = (surface.centroids()[i], surface.normals()[i], surface.areas()[i]);
        for &x in targets {
            let d = sub(x, c);
            if dot(d, n) <= 0.0 {
                continue;
            }
            let r = norm(d);
            let (s, co) = det::sin_cos(k * r);
            total = total + C64::new(co, s).scale(field[i] * area * 2.0 / (4.0 * PI * r));
        }
    }
    C64::new(0.0, -omega * gas.0) * total
}

/// The same transfers by reciprocity: a unit point source at the microphone
/// scattered by every body (PlainCbie, the finite-source owner), its total
/// surface pressure integrated against each trace.
fn reverse(surface: &SpherePanels, gas: (f64, f64), traces: &[Vec<f64>], mic: [f64; 3], omega: f64) -> Result<Vec<C64>, String> {
    let monopole = [Monopole { position_m: mic, volume_velocity_m3_s: C64::new(1.0, 0.0) }];
    let scattered = solve_scattering(
        surface,
        omega / gas.1,
        medium(gas),
        Admittance::Rigid,
        IncidentField::Monopoles(&monopole),
        CLEARANCE_M,
        Formulation::PlainCbie,
    )
    .map_err(|e| format!("reciprocal scattering: {e}"))?;
    let p = scattered.solution().total_pressure_pa();
    Ok(traces
        .iter()
        .map(|t| {
            t.iter()
                .zip(&p)
                .zip(surface.areas())
                .fold(C64::ZERO, |s, ((&v, &pp), &a)| s + pp.scale(v * a))
        })
        .collect())
}

fn db(z: C64) -> f64 {
    20.0 * z.abs().log10()
}

// ------------------------------------------------------------------ the bake

/// Sizing probes for the bake's own gates, before any boundary solve, each
/// alone and behind a pure delay of `delay_samples` (exp(-i omega t)): the
/// complement source, then (i omega / omega_stop)^p for p = 1..=`order`. The
/// gates measure every source against its own energy, and a high-order head
/// shape's row rises steeply to f_x, so its energy sits in the taper's
/// transition: the complement alone under-sizes the FIR for it. Like every
/// physical row the probes vanish at DC.
fn sizing_spectrum(fs: f64, stop: f64, delay_samples: f64, order: usize) -> impl Fn(f64) -> Result<Vec<C64>, HelmholtzError> {
    move |omega: f64| {
        let (s, c) = det::sin_cos(omega * delay_samples / fs);
        let late = C64::new(c, s);
        let step = C64::new(0.0, omega / (TAU * stop));
        let mut power = C64::new(1.0, 0.0);
        let mut row = Vec::with_capacity(8 * (order + 1));
        for p in 0..=order {
            let probe = if p == 0 {
                complement_row(omega, fs)
            } else {
                power = power * step;
                power
            };
            row.extend([probe, C64::ZERO, C64::ZERO, C64::ZERO, probe * late, C64::ZERO, C64::ZERO, C64::ZERO]);
        }
        Ok(row)
    }
}

/// The FIR declaration of `sources` analytic probe sources at one microphone.
fn probe_declaration(fs: f64, pass: f64, stop: f64, n: usize, latency: usize, sources: usize) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: fs_blake3::hash_domain("frankensim.drum-observer.sizing", b"taper"),
        sources: (0..sources).map(|i| SourceId::StructuralMode(format!("sizing probe {i}"))).collect(),
        receivers: vec![receiver([0.0; 3])],
        fs_hz: fs,
        geometric_delays_s: vec![0.0; sources],
        pass_band_upper_hz: pass,
        taper_stop_hz: stop,
        filter_latency_samples: latency,
        fft_lengths: vec![n / 8, n / 4, n / 2, n],
        tail_energy_tolerance: TAIL_TOLERANCE,
        dtft_relative_tolerance: DTFT_TOLERANCE,
        producer_uncertainty: ProducerUncertainty::ExactAnalytic,
        maximum_frequency_entries: usize::MAX,
    }
}

/// Smallest power-of-two FFT and latency for which the bake's own gates pass
/// on every probe of [`sizing_spectrum`] up to `order`, alone and behind
/// `spread` samples of response.
fn size_fir(fs: f64, pass: f64, stop: f64, spread: f64, order: usize) -> Result<(usize, usize), String> {
    for log2 in 12..=15 {
        let n = 1usize << log2;
        for step in 1..16 {
            let latency = n * step / 64;
            let decl = probe_declaration(fs, pass, stop, n, latency, 2 * (order + 1));
            if bake_observation_transfers(sizing_spectrum(fs, stop, spread, order), &decl, || false).is_ok() {
                return Ok((n, latency));
            }
        }
    }
    Err(format!("no FFT up to 32768 holds the taper, (i omega)^{order} and {spread:.0} samples of response"))
}

/// Write every solved row (exact bits: the omega key, then re and im of each
/// channel) to `path`, so a refused bake keeps its solves for diagnosis.
fn dump_rows(path: &str, rows: &BTreeMap<u64, Vec<C64>>) -> String {
    let mut t = String::new();
    for (key, row) in rows {
        let _ = write!(t, "{key:016x}");
        for z in row {
            let _ = write!(t, "\t{:016x}\t{:016x}", z.re.to_bits(), z.im.to_bits());
        }
        t.push('\n');
    }
    match std::fs::write(path, t) {
        Ok(()) => format!("; the {} solved rows are in {path}", rows.len()),
        Err(e) => format!("; the {} solved rows could not be written to {path}: {e}", rows.len()),
    }
}

/// Print a receipt line as soon as it is known, and keep it for the table.
fn note(receipt: &mut Vec<String>, line: String) {
    eprintln!("{line}");
    receipt.push(line);
}

struct Bake {
    text: String,
}

/// Bake one piece's table at `fs` [Hz].
/// A refusal after the solves writes the solved rows to `dump`, if given.
fn bake(piece: &str, fs: f64, dump: Option<&str>) -> Result<Bake, String> {
    let card = card(piece).ok_or_else(|| format!("unknown piece {piece}"))?;
    let gas_state = fs_material::gas::GasState::try_new_moist_air(GAS_STATE.0, GAS_STATE.1, GAS_STATE.2)
        .map_err(|e| format!("gas: {e:?}"))?;
    let gas = (gas_state.density, gas_state.sound_speed);
    let mesh = scene(&card, true);
    let surface = SpherePanels::from_triangles(mesh.tris.clone()).map_err(|e| format!("surface: {e:?}"))?;
    let policy = GeometryPolicy::new(&surface).map_err(|e| e.to_string())?;
    let k_x = policy.plain_cbie_limit();
    let f_x = k_x * gas.1 / TAU;
    let stop = f_x;
    let pass = card.pass_fraction * f_x;
    if !(stop < 0.5 * fs) {
        return Err(format!("f_x {f_x:.1} Hz is not below Nyquist at {fs} Hz"));
    }
    let sources = retained_sources();
    let ids: Vec<SourceId> = sources.iter().map(Source::id).collect();
    let traces: Vec<Vec<f64>> = sources.iter().map(|s| source_trace(&mesh, card.radius, s)).collect();
    let mut receipt = Vec::new();
    note(
        &mut receipt,
        format!(
            "{}: radius {} m, depth {} m; {} panels on the drum, {} with the image body; f_x {:.2} Hz (0.8 of the box Dirichlet bound), pass {:.2} Hz",
            card.piece,
            card.radius,
            card.depth,
            mesh.tris.len() / 2,
            mesh.tris.len(),
            f_x,
            pass
        ),
    );

    // Both ways, at four frequencies up to f_x.
    let mut eps = 0.0_f64;
    let mut eps_lines = Vec::new();
    let mut at_fx = Vec::new();
    // The top check sits just inside the band: omega / c at exactly f_x may
    // round above the CBIE limit, and that bin is refused, not rerouted.
    for fraction in [0.25, 0.5, 0.75, TOP_CHECK] {
        let omega = TAU * fraction * f_x;
        let fwd: Vec<C64> = forward(&surface, gas, &ids, &traces, &[card.mic], omega)?.into_iter().map(|r| r[0]).collect();
        let rev = reverse(&surface, gas, &traces, card.mic, omega)?;
        let scale = fwd.iter().fold(0.0_f64, |m, z| m.max(z.abs()));
        let worst = fwd
            .iter()
            .zip(&rev)
            .enumerate()
            .map(|(i, (a, b))| ((*a - *b).abs() / scale, i))
            .fold((0.0_f64, 0), |m, x| if x.0 > m.0 { x } else { m });
        eps = eps.max(worst.0);
        eps_lines.push(format!("{:.2} Hz: {:.3e} ({})", fraction * f_x, worst.0, ids[worst.1].name_text()));
        if fraction == TOP_CHECK {
            at_fx = fwd;
        }
    }
    note(&mut receipt, format!("reciprocity (both ways, re the largest transfer): {}", eps_lines.join("; ")));
    if !(eps <= RECIPROCITY_LIMIT) {
        return Err(format!("reciprocity {eps:.3e} exceeds {RECIPROCITY_LIMIT}: {}", eps_lines.join("; ")));
    }

    // Crossover step: closed body over baffled Rayleigh for the batter's
    // (0,1) and (0,2) at TOP_CHECK f_x, heard at the microphone and its image.
    let image_mic = mirror_point(card.mic, card.floor.0, card.floor.1);
    let mut steps = Vec::new();
    for n in [1, 2] {
        let i = sources
            .iter()
            .position(|s| s.head == Part::Batter && s.shape.m == 0 && s.shape.n == n)
            .ok_or("missing batter (0,n)")?;
        let r = baffled(&mesh, &surface, gas, &traces[i], &[card.mic, image_mic], TAU * TOP_CHECK * f_x);
        steps.push(db(at_fx[i]) - db(r));
    }
    note(
        &mut receipt,
        format!(
            "crossover step at {TOP_CHECK} f_x, closed body re baffled Rayleigh, batter (0,1) {:+.2} dB, (0,2) {:+.2} dB",
            steps[0], steps[1]
        ),
    );

    // Size the FIR on the taper, the complement and steep rows, each alone
    // and behind the response spread over both bodies.
    let distances: Vec<f64> = surface.centroids().iter().map(|c| norm(sub(card.mic, *c))).collect();
    let (d_min, d_max) = distances.iter().fold((f64::INFINITY, 0.0_f64), |(a, b), &d| (a.min(d), b.max(d)));
    let spread = ((d_max - d_min) + 2.0 * (2.0 * card.radius + card.depth)) / gas.1 * fs;
    let (n_fft, probe_latency) = size_fir(fs, pass, stop, spread, STEEP_PROBE_ORDER)?;
    note(
        &mut receipt,
        format!(
            "FIR sizing: spread {spread:.0} samples; the complement and (i omega / omega_stop)^p, p <= {STEEP_PROBE_ORDER}, alone and behind the spread, hold at FFT {n_fft}, latency {probe_latency}"
        ),
    );

    let n_panels = surface.areas().len();
    let obs = BemObservation::new_routed(
        &surface,
        medium(gas),
        &ids,
        &traces,
        &[receiver(card.mic)],
        CLEARANCE_M,
        PanelRule::Triangle { max_pairs: n_panels * n_panels },
    )
    .map_err(|e| format!("observation admission: {e}"))?;
    let mut all_ids = ids.clone();
    all_ids.push(SourceId::StructuralMode(COMPLEMENT_SOURCE.into()));
    let mut delays = obs.geometric_delays_s();
    delays.push(0.0);
    let mut provenance = DomainHasher::new("frankensim.vkplate.closed-body-observer.producer");
    provenance.update(obs.operator_hash().as_bytes());
    provenance.update(COMPLEMENT_SOURCE.as_bytes());
    let producer_hash = provenance.finalize();
    let evidence = fs_blake3::hash_domain("frankensim.drum-observer.reciprocity", eps_lines.join("\n").as_bytes());
    let declare = |latency: usize| TransferDeclaration {
        producer_hash,
        sources: all_ids.clone(),
        receivers: vec![receiver(card.mic)],
        fs_hz: fs,
        geometric_delays_s: delays.clone(),
        pass_band_upper_hz: pass,
        taper_stop_hz: stop,
        filter_latency_samples: latency,
        fft_lengths: vec![n_fft / 8, n_fft / 4, n_fft / 2, n_fft],
        tail_energy_tolerance: TAIL_TOLERANCE,
        dtft_relative_tolerance: DTFT_TOLERANCE,
        producer_uncertainty: ProducerUncertainty::RelativeEstimate { relative: eps, source_hash: evidence },
        maximum_frequency_entries: usize::MAX,
    };
    // Every bin is solved once, at the probes' latency. If the rows refuse
    // there, the latency steps up by n/64 on the same rows (no new solve)
    // until fs-bem's own gates pass; the FFT, and so every bin, stays.
    let mut rows: BTreeMap<u64, Vec<C64>> = BTreeMap::new();
    let started = std::time::Instant::now();
    let mut refused_at = Vec::new();
    let mut found = None;
    for step in probe_latency * 64 / n_fft..16 {
        let latency = n_fft * step / 64;
        let attempt = bake_observation_transfers(
            |omega| {
                if let Some(row) = rows.get(&omega.to_bits()) {
                    return Ok(row.clone());
                }
                if policy.formulation(omega / gas.1)? != Formulation::PlainCbie {
                    return Err(HelmholtzError::BadParameter { what: "a bake bin above the CBIE band" });
                }
                let mut row = obs.frequency_sample(omega)?;
                row.extend([complement_row(omega, fs), C64::ZERO, C64::ZERO, C64::ZERO]);
                rows.insert(omega.to_bits(), row.clone());
                if rows.len() % 10 == 0 {
                    eprintln!("  {} bins, {:.1} Hz, {:.0} s", rows.len(), omega / TAU, started.elapsed().as_secs_f64());
                }
                Ok(row)
            },
            &declare(latency),
            || false,
        );
        match attempt {
            Ok(baked) => {
                found = Some((baked, latency));
                break;
            }
            Err(TransferError::NotConverged { tail, wrap, dtft }) => {
                let line = format!("latency {latency}: tail {tail:.2e} wrap {wrap:.2e} dtft {dtft:.2e}");
                eprintln!("  refused at {line}");
                refused_at.push(line);
            }
            Err(e) => {
                let kept = dump.map_or(String::new(), |path| dump_rows(path, &rows));
                return Err(format!("bake refused: {e}{kept}"));
            }
        }
    }
    let Some((baked, latency)) = found else {
        let kept = dump.map_or(String::new(), |path| dump_rows(path, &rows));
        return Err(format!(
            "bake refused at every latency up to {} on FFT {n_fft}: {}{kept}",
            n_fft * 15 / 64,
            refused_at.join("; ")
        ));
    };
    if !refused_at.is_empty() {
        note(&mut receipt, format!("the rows refused at the probes' latency and above: {}", refused_at.join("; ")));
    }
    note(&mut receipt, format!(
        "FIR: rate {fs} Hz, FFT {n_fft}, latency {latency} samples, {} bins below the stop, ladder {}",
        rows.len(),
        baked
            .ladder
            .iter()
            .map(|r| format!(
                "{}: tail {:.1e} wrap {:.1e} dtft {}",
                r.n_fft,
                r.tail_energy_fraction,
                r.wrap_energy_fraction,
                r.dtft_relative_change.map_or("-".into(), |c| format!("{c:.1e}"))
            ))
            .collect::<Vec<_>>()
            .join(", ")
    ));
    let text = write_table(&card, gas, fs, &mesh, &sources, &all_ids, &baked, &rows, n_fft, (f_x, pass), eps, &steps, &receipt)?;
    Ok(Bake { text })
}

trait NameText {
    fn name_text(&self) -> String;
}
impl NameText for SourceId {
    fn name_text(&self) -> String {
        match self {
            SourceId::StructuralMode(s)
            | SourceId::ApertureVolumeFlow(s)
            | SourceId::PointVolumeFlow(s)
            | SourceId::PointForce(s) => s.split(':').next().unwrap_or("").to_string(),
        }
    }
}

/// Energy, energy centroid and the transfer at `omega` [rad/sample] of a kernel.
fn kernel_stats(taps: &[f64], integer: usize, omega: f64) -> (f64, f64, C64) {
    let energy: f64 = taps.iter().map(|t| t * t).sum();
    let moment: f64 = taps.iter().enumerate().map(|(i, t)| i as f64 * t * t).sum();
    let mut h = C64::ZERO;
    for (i, &t) in taps.iter().enumerate() {
        let (s, c) = det::sin_cos(-omega * (integer + i) as f64);
        h = h + C64::new(c, s).scale(t);
    }
    (energy, if energy > 0.0 { moment / energy } else { 0.0 }, h)
}

#[allow(clippy::too_many_arguments)]
fn write_table(
    card: &Card,
    gas: (f64, f64),
    fs: f64,
    mesh: &Mesh,
    sources: &[Source],
    ids: &[SourceId],
    baked: &ObservationTransferBake,
    rows: &BTreeMap<u64, Vec<C64>>,
    n_fft: usize,
    (f_x, pass): (f64, f64),
    eps: f64,
    steps: &[f64],
    receipt: &[String],
) -> Result<String, String> {
    let mut t = String::new();
    let _ = writeln!(t, "{TABLE_FORMAT}");
    let _ = writeln!(
        t,
        "route\t{ROUTE}\tfs_bem::observation_transfers routed owner, triangle rule, PlainCbie only (crates/fs-bem/examples/drum_observer.rs)"
    );
    let _ = writeln!(t, "bake-args\t{}\t{:?}", card.piece, fs);
    let _ = writeln!(t, "piece\t{}", card.piece);
    let _ = writeln!(t, "drum\t{:?}\t{:?}", card.radius, card.depth);
    let _ = writeln!(t, "heads\tbatter\t0.0\t0.0\t0.0\t1.0\tresonant\t0.0\t0.0\t{:?}\t-1.0", -card.depth);
    let _ = writeln!(t, "mic\t{:?}\t{:?}\t{:?}", card.mic[0], card.mic[1], card.mic[2]);
    let (q, nf) = card.floor;
    let _ = writeln!(t, "floor\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}", q[0], q[1], q[2], nf[0], nf[1], nf[2]);
    let _ = writeln!(t, "gas\t{:?}\t{:?}", gas.0, gas.1);
    let _ = writeln!(t, "rate\t{fs:?}");
    let _ = writeln!(
        t,
        "mesh\t{SECTORS} rim sectors\t{CAP_RINGS} cap rings\t{} panels on the drum\t{} with the image body",
        mesh.tris.len() / 2,
        mesh.tris.len()
    );
    let _ = writeln!(t, "band\t{pass:?}\t{f_x:?}");
    let _ = writeln!(t, "latency\t{}", baked.filter_latency_samples);
    let _ = writeln!(t, "fft\t{n_fft}");
    let _ = writeln!(t, "declaration\t{}", baked.declaration_hash.to_hex());
    let _ = writeln!(t, "reciprocity\t{eps:e}");
    let _ = writeln!(t, "crossover-step-db\t{:?}\t{:?}", steps[0], steps[1]);
    let _ = writeln!(t, "mode-set\t{MODE_SET}\tj_mn <= {ZERO_LIMIT}, both heads, cos and sin parities, unit mean square");
    for line in receipt {
        let _ = writeln!(t, "note\t{line}");
    }
    for (i, s) in sources.iter().enumerate() {
        let _ = writeln!(
            t,
            "source\t{i}\t{}\t{}\t{}\t{}\t{:?}",
            s.head.label(),
            s.shape.m,
            s.shape.n,
            if s.sine { "sin" } else { "cos" },
            s.shape.zero
        );
        for (x, y) in [(0.31, 0.12), (-0.55, 0.47), (0.08, -0.83)] {
            let _ = writeln!(t, "probe\t{i}\t{x:?}\t{y:?}\t{:?}", s.shape.value(x, y, s.sine));
        }
    }
    let probe_omega = TAU * 0.5 * pass / fs;
    for i in 0..ids.len() {
        let kernel = baked
            .kernels
            .iter()
            .find(|k| k.input == i && k.output == 0)
            .ok_or("missing pressure kernel")?;
        let (energy, centre, h) = kernel_stats(&kernel.taps, kernel.integer_delay_samples, probe_omega);
        let _ = writeln!(
            t,
            "kernel\t{i}\t{}\t{}\t{energy:?}\t{centre:?}\t{:?}\t{:?}\t{}",
            kernel.integer_delay_samples,
            kernel.taps.len(),
            h.re,
            h.im,
            if i == sources.len() { "complement" } else { "shape" }
        );
    }
    let mut bin = 1usize;
    loop {
        let f = fs * bin as f64 / n_fft as f64;
        if f >= baked.taper_stop_hz {
            break;
        }
        let row = rows
            .get(&(TAU * f).to_bits())
            .ok_or_else(|| format!("bin {bin} ({f:.2} Hz) was not solved"))?;
        for i in 0..sources.len() {
            let z = row[4 * i];
            let _ = writeln!(t, "row\t{i}\t{bin}\t{:?}\t{:?}", z.re, z.im);
        }
        bin += 1;
    }
    let _ = writeln!(t, "end\t{}", bin - 1);
    Ok(t)
}

fn main() {
    let mut a: Vec<String> = std::env::args().collect();
    // A scheduling budget for the dense assembly and LU only: every entry and
    // factor is bit-identical for every count, so the table never records it.
    let workers = match a.iter().position(|x| x == "--workers") {
        Some(i) if i + 1 < a.len() => {
            let w: usize = a[i + 1].parse().unwrap_or_else(|_| {
                eprintln!("refused: --workers takes a count");
                std::process::exit(2)
            });
            a.drain(i..i + 2);
            w.max(1)
        }
        Some(_) => {
            eprintln!("refused: --workers takes a count");
            std::process::exit(2)
        }
        None => 1,
    };
    eprintln!("dense workers: {workers}");
    fs_bem::helmholtz::with_dense_workers(workers, || run(&a));
}

fn run(a: &[String]) {
    if a.len() == 3 && a[1] == "--verify" {
        let committed = std::fs::read_to_string(&a[2]).unwrap_or_else(|e| {
            eprintln!("refused: cannot read {}: {e}", a[2]);
            std::process::exit(1)
        });
        let args: Vec<&str> = committed
            .lines()
            .find_map(|l| l.strip_prefix("bake-args\t"))
            .unwrap_or_else(|| {
                eprintln!("refused: no bake-args line");
                std::process::exit(1)
            })
            .split('\t')
            .collect();
        let fs: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or_else(|| {
            eprintln!("refused: bad bake-args");
            std::process::exit(1)
        });
        match bake(args[0], fs, None) {
            Ok(b) if b.text == committed => {
                eprintln!("verified: re-baked {} bytes are identical to {}", b.text.len(), a[2]);
            }
            Ok(b) => {
                let line = b.text.lines().zip(committed.lines()).position(|(x, y)| x != y);
                eprintln!(
                    "refused: the re-bake differs from {} ({} vs {} bytes; first differing line {:?})",
                    a[2],
                    b.text.len(),
                    committed.len(),
                    line.map(|l| l + 1)
                );
                std::process::exit(1);
            }
            Err(e) => {
                eprintln!("refused: {e}");
                std::process::exit(1);
            }
        }
        return;
    }
    if a.len() != 4 {
        eprintln!("usage: drum_observer PIECE RATE_HZ OUT.tsv [--workers N] | drum_observer --verify TABLE.tsv [--workers N]");
        std::process::exit(2);
    }
    let fs: f64 = a[2].parse().expect("RATE_HZ");
    let dump = format!("{}.refused-rows", a[3]);
    match bake(&a[1], fs, Some(dump.as_str())) {
        Ok(b) => {
            std::fs::write(&a[3], &b.text).expect("write the table");
            eprintln!("{} bytes written to {}", b.text.len(), a[3]);
        }
        Err(e) => {
            eprintln!("refused: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
#[path = "drum_observer/tests.rs"]
mod tests;
