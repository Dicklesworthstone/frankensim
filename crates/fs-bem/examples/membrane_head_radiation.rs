//! Bake `crates/fs-vkplate/data/head-radiation-planar-baffled.tsv`: the
//! dimensionless exterior loading `zhat_rs(ka) = Z_rs / (rho c a^2)` of
//! clamped-membrane shapes `J_m(j_mn r) cos(m theta)` (Bessel zero
//! `j_mn <= 20`) on a unit disk in an infinite rigid baffle, through
//! `fs_bem::modal_loading`'s power-consistent Rayleigh route. fs-vkplate
//! loads drum heads' `m >= 1` motion from it; its `m = 0` family only states
//! the error of fs-vkplate's compact volume model.
//!
//! With `rho = c = a = 1` the angular frequency IS `ka`. Polar cells carry
//! each shape at their centroid. The shapes here duplicate fs-vkplate's
//! `radiation::ModeShape` (fs-vkplate is a standalone workspace); the table's
//! `probe` rows let fs-vkplate refuse any difference.
//!
//! Measured here and written into the table:
//! - `zmid` rows: diagonal samples at every geometric midpoint of the ka
//!   grid, from which fs-vkplate measures its own interpolation error;
//! - discretization: the largest change of a diagonal between the two finest
//!   rungs of a 24/33/45-ring ladder (observed order printed);
//! - resolution: per family, the first ka from which every diagonal
//!   resistance is 1000 times above its own rounding floor (a high `m`
//!   family's small-ka resistance is a multipole cancellation);
//! - symmetry: cos/sin parity and cross-family reactance, printed.
//!
//! Run: `cargo run --release -p fs-bem --example membrane_head_radiation -- OUT.tsv`

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::rayleigh_modal::AreaSample;
use fs_bem::modal_loading::{
    LoadingRoute, ModalLoadingBake, SurfaceModeTraces, modal_radiation_loading,
};
use std::f64::consts::PI;
use std::fmt::Write as _;

const TABLE_FORMAT: &str = "frankensim.vkplate-head-radiation.v1";
const ROUTE: &str = "planar-baffled-rayleigh";
const GEOMETRY: &str = "unit-disk-in-infinite-rigid-baffle";
const MODE_SET: &str = "clamped-membrane-bessel";
const ZERO_LIMIT: f64 = 20.0;
/// The ka grid: log-spaced (PER_DECADE per decade) from KA_MIN below ka = 1,
/// where the loads are power laws, then uniform steps of KA_STEP up to
/// KA_MAX, where the resistance oscillates through coincidence on a scale of
/// order one in ka. Kit heads in air reach ka of about 8.5; the 64 cm timpano
/// oracle about 8.3.
const KA_MIN: f64 = 0.02;
const PER_DECADE: usize = 20;
const KA_STEP: f64 = 0.2;
const KA_MAX: f64 = 12.0;
/// Ring counts of the ladder. The coarsest keeps the Rayleigh owner's six
/// equivalent-disc cells per wavelength at the largest ka (20 rings do not).
const RUNGS: [usize; 3] = [24, 33, 45];
/// Largest relative change of a resistance between the two finest rungs for
/// it to count as resolved.
const RESOLVED_LADDER_CHANGE: f64 = 0.05;

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
        Shape {
            m,
            n,
            zero,
            scale: angular / bessel_j(m + 1, zero).abs(),
        }
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

/// Polar cells of the unit disk with `rings` rings: (centroid, area).
/// Every ring's sector count is a multiple of this. The cell set is then
/// invariant under rotation by `2 pi / SYMMETRY` and reflection, so families
/// with `m + m' < SYMMETRY` (all of the table's) cannot couple through the
/// quadrature, and the cos and sin parities sample identically.
const SYMMETRY: usize = 32;

fn cells(rings: usize) -> Vec<([f64; 2], f64)> {
    let mut out = Vec::new();
    for i in 0..rings {
        let (r1, r2) = (i as f64 / rings as f64, (i + 1) as f64 / rings as f64);
        // Rounded up: no cell is larger than the square-cell nominal, so the
        // Rayleigh owner's cells-per-wavelength floor holds on every ring.
        let nominal = 2.0 * PI * (i as f64 + 0.5) / SYMMETRY as f64;
        let sectors = SYMMETRY * (nominal.ceil() as usize).max(1);
        let dt = 2.0 * PI / sectors as f64;
        let area = 0.5 * (r2 * r2 - r1 * r1) * dt;
        let rc = 2.0 / 3.0 * (r2.powi(3) - r1.powi(3)) / (r2 * r2 - r1 * r1) * (dt / 2.0).sin()
            / (dt / 2.0);
        for j in 0..sectors {
            let t = (j as f64 + 0.5) * dt;
            out.push(([rc * t.cos(), rc * t.sin()], area));
        }
    }
    out
}

/// Fields: every mode in cos parity, then sin-parity checks of (1,1), (2,1).
fn fields(modes: &[Shape]) -> Vec<(usize, bool)> {
    let mut f: Vec<(usize, bool)> = (0..modes.len()).map(|i| (i, false)).collect();
    for (i, s) in modes.iter().enumerate() {
        if s.n == 1 && (s.m == 1 || s.m == 2) {
            f.push((i, true));
        }
    }
    f
}

fn bake(modes: &[Shape], rings: usize, omegas: &[f64]) -> (ModalLoadingBake, Vec<f64>) {
    let fields = fields(modes);
    let samples: Vec<AreaSample> = cells(rings)
        .into_iter()
        .map(|(p, area)| AreaSample {
            position_m: p,
            area_m2: area,
            mode_shape: fields
                .iter()
                .map(|&(i, sine)| modes[i].value(p[0], p[1], sine))
                .collect(),
        })
        .collect();
    // Per field, sum |area * shape|: the resistance rounding floor's scale.
    let abs_sums: Vec<f64> = (0..fields.len())
        .map(|f| {
            samples
                .iter()
                .map(|s| (s.area_m2 * s.mode_shape[f]).abs())
                .sum()
        })
        .collect();
    let mut provenance =
        format!("{TABLE_FORMAT}\n{ROUTE}\n{GEOMETRY}\n{MODE_SET}\nrings {rings}\n");
    for s in modes {
        let _ = writeln!(provenance, "{} {} {:?}", s.m, s.n, s.zero);
    }
    let provenance = fs_blake3::hash_domain(
        "frankensim.vkplate.head-radiation-bake",
        provenance.as_bytes(),
    );
    let started = std::time::Instant::now();
    let mut polls = 0usize;
    let bake = modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: &samples },
        LoadingRoute::PlanarBaffledRayleigh,
        omegas,
        Medium {
            density: 1.0,
            sound_speed: 1.0,
        },
        provenance,
        1 << 26,
        || {
            polls += 1;
            if polls % 20 == 0 {
                eprintln!(
                    "  rings {rings}: poll {polls}, {:.1} s",
                    started.elapsed().as_secs_f64()
                );
            }
            false
        },
    )
    .unwrap_or_else(|e| panic!("bake refused at {rings} rings: {e}"));
    eprintln!(
        "rings {rings}: {} cells, {} fields, {} frequencies, {:.1} s",
        samples.len(),
        fields.len(),
        omegas.len(),
        started.elapsed().as_secs_f64()
    );
    (bake, abs_sums)
}

fn entry(bake: &ModalLoadingBake, row: usize, r: usize, s: usize) -> (f64, f64) {
    let z = bake.rows()[row].impedance()[r * bake.modes() + s];
    (z.re, z.im)
}

fn main() {
    let out = std::env::args()
        .nth(1)
        .expect("usage: membrane_head_radiation OUT.tsv");
    let mut modes = Vec::new();
    for m in 0.. {
        if Shape::new(m, 1).zero > ZERO_LIMIT {
            break;
        }
        for n in 1.. {
            let shape = Shape::new(m, n);
            if shape.zero > ZERO_LIMIT {
                break;
            }
            assert!(
                bessel_j(m, shape.zero).abs() < 1e-12,
                "J_{m}({}) is not a zero",
                shape.zero
            );
            modes.push(shape);
        }
    }
    let step = 10f64.powf(1.0 / PER_DECADE as f64);
    let mut ka: Vec<f64> = (0..)
        .map(|i| KA_MIN * step.powi(i))
        .take_while(|&k| k < 1.0)
        .collect();
    let uniform = ((KA_MAX - 1.0) / KA_STEP).round() as usize;
    ka.extend((0..=uniform).map(|i| 1.0 + i as f64 * KA_STEP));
    let mids: Vec<f64> = ka.windows(2).map(|w| (w[0] * w[1]).sqrt()).collect();
    let mut both: Vec<f64> = ka.iter().chain(&mids).copied().collect();
    both.sort_by(f64::total_cmp);
    let fine = RUNGS[2];
    // The cheap coarse rungs first: a refusal there costs seconds.
    let coarse: Vec<ModalLoadingBake> = RUNGS[..2]
        .iter()
        .map(|&rings| bake(&modes, rings, &ka).0)
        .collect();
    let (main, abs_sums) = bake(&modes, fine, &both);
    let grid_rows: Vec<usize> = ka
        .iter()
        .map(|k| both.iter().position(|b| b == k).unwrap())
        .collect();
    let mid_rows: Vec<usize> = mids
        .iter()
        .map(|k| both.iter().position(|b| b == k).unwrap())
        .collect();

    let mut families: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, s) in modes.iter().enumerate() {
        match families.iter_mut().find(|(m, _)| *m == s.m) {
            Some((_, f)) => f.push(i),
            None => families.push((s.m, vec![i])),
        }
    }
    let resolved = |row: usize, r: usize| {
        let k = main.rows()[row].omega();
        let floor = 64.0 * f64::EPSILON * k * k / (2.0 * PI) * abs_sums[r] * abs_sums[r];
        entry(&main, row, r, r).0 > 1e3 * floor
    };
    // A small-ka resistance of a high family is a cancellation of lower
    // multipoles: below roundoff, and on a coarse quadrature below its
    // aliasing error. It is trusted where it clears the floor on the fine
    // rung AND the two finest rungs agree within RESOLVED_LADDER_CHANGE.
    let converged = |i: usize, r: usize| {
        let (middle, high) = (
            entry(&coarse[1], i, r, r).0,
            entry(&main, grid_rows[i], r, r).0,
        );
        resolved(grid_rows[i], r)
            && middle > 0.0
            && (middle - high).abs() <= RESOLVED_LADDER_CHANGE * high
    };
    let resolved_from: Vec<usize> = families
        .iter()
        .map(|(_, members)| {
            (0..ka.len())
                .find(|&first| (first..ka.len()).all(|i| members.iter().all(|&r| converged(i, r))))
                .unwrap_or(ka.len() - 1)
        })
        .collect();

    // Symmetry at the finest rung over the grid.
    let (mut parity, mut cross) = (0.0_f64, 0.0_f64);
    let checks: Vec<(usize, usize)> = fields(&modes)
        .iter()
        .enumerate()
        .filter(|(_, f)| f.1)
        .map(|(col, f)| (f.0, col))
        .collect();
    for &row in &grid_rows {
        for &(cos, sin) in &checks {
            let (a, b) = (entry(&main, row, cos, cos), entry(&main, row, sin, sin));
            parity = parity.max((a.1 - b.1).abs() / a.1.abs());
            if resolved(row, cos) {
                parity = parity.max((a.0 - b.0).abs() / a.0);
            }
        }
        for r in 0..modes.len() {
            for s in 0..modes.len() {
                if modes[r].m != modes[s].m {
                    let scale = (entry(&main, row, r, r).1 * entry(&main, row, s, s).1)
                        .abs()
                        .sqrt();
                    cross = cross.max(entry(&main, row, r, s).1.abs() / scale);
                }
            }
        }
    }

    // Discretization ladder on the ka grid (diagonals): reactance everywhere,
    // resistance where resolved. The worst entry of each part is printed.
    let (mut d_mid, mut d_fine) = ([0.0_f64; 2], [0.0_f64; 2]);
    let mut worst = [(0.0_f64, String::new()), (0.0_f64, String::new())];
    for i in 0..ka.len() {
        for (f, (_, members)) in families.iter().enumerate() {
            for &r in members {
                let (low, middle, high) = (
                    entry(&coarse[0], i, r, r),
                    entry(&coarse[1], i, r, r),
                    entry(&main, grid_rows[i], r, r),
                );
                // Reactance relative to the load magnitude |Z| (it crosses
                // zero where a shape turns supersonic); resolved resistance
                // relative to itself.
                let magnitude = high.0.hypot(high.1);
                let mut parts = vec![(1, low.1, middle.1, magnitude)];
                if i >= resolved_from[f] {
                    parts.push((0, low.0, middle.0, high.0));
                }
                for (part, l, m, scale) in parts {
                    let h = if part == 0 { high.0 } else { high.1 };
                    d_mid[part] = d_mid[part].max((l - m).abs() / scale);
                    let change = (m - h).abs() / scale;
                    d_fine[part] = d_fine[part].max(change);
                    if change >= worst[part].0 {
                        worst[part] = (
                            change,
                            format!("({},{}) at ka {:.4}", modes[r].m, modes[r].n, ka[i]),
                        );
                    }
                }
            }
        }
    }
    // The two ring ratios are within 2 % of each other; the finer one is used.
    let ratio = (RUNGS[2] as f64 / RUNGS[1] as f64).ln();
    let order = [0, 1].map(|p| (d_mid[p] / d_fine[p]).ln() / ratio);
    let bound = d_fine[0].max(d_fine[1]);
    let [a, b, c] = RUNGS;
    let notes = [
        format!(
            "ladder {a}/{b}/{c} rings (sector counts multiples of {SYMMETRY}), largest change of a diagonal: reactance (relative to |Z|) {a}->{b} {:.3e}, {b}->{c} {:.3e} (observed order {:.2}, worst {}); resolved resistance {a}->{b} {:.3e}, {b}->{c} {:.3e} (observed order {:.2}, worst {}); the bound is the larger {b}->{c} change",
            d_mid[1], d_fine[1], order[1], worst[1].1, d_mid[0], d_fine[0], order[0], worst[0].1
        ),
        format!(
            "symmetry at {fine} rings: cos/sin parity <= {parity:.2e}, cross-family reactance <= {cross:.2e} of the diagonal scale"
        ),
        format!(
            "resistance resolution: a family's resistance is used from the first ka at which every diagonal exceeds 1000 rounding floors, 64 eps (ka^2 / 2 pi) (sum over cells of area |phi|)^2, and the {b}->{c} rungs agree within {RESOLVED_LADDER_CHANGE}"
        ),
    ];

    let mut t = String::new();
    let _ = writeln!(t, "{TABLE_FORMAT}");
    let _ = writeln!(t, "route\t{ROUTE}");
    let _ = writeln!(
        t,
        "geometry\t{GEOMETRY}\tpolar centroid cells: {fine} rings, {} cells\tdeclaration {}",
        cells(fine).len(),
        main.declaration_hash()
    );
    let _ = writeln!(
        t,
        "mode-set\t{MODE_SET}\tj_mn <= {ZERO_LIMIT}, m >= 0 (m = 0 is a diagnostic family), cos parity tabulated (sin identical by rotation), shapes normalized to unit mean square"
    );
    let _ = writeln!(
        t,
        "declaration\tfs_bem::modal_loading power-consistent Rayleigh (crates/fs-bem/examples/membrane_head_radiation.rs), rho = c = a = 1 so omega = ka; zhat = Z / (rho c a^2); exp(-i omega t): mass-like loads have negative imaginary part"
    );
    for note in &notes {
        let _ = writeln!(t, "note\t{note}");
    }
    let _ = writeln!(t, "discretization-bound\t{bound:e}");
    for s in &modes {
        let _ = writeln!(t, "mode\t{}\t{}\t{:?}", s.m, s.n, s.zero);
    }
    for (i, s) in modes.iter().enumerate() {
        for (x, y) in [(0.31, 0.12), (-0.55, 0.47), (0.08, -0.83)] {
            let _ = writeln!(t, "probe\t{i}\t{x:?}\t{y:?}\t{:?}", s.value(x, y, false));
        }
    }
    for k in &ka {
        let _ = writeln!(t, "ka\t{k:?}");
    }
    for ((m, _), first) in families.iter().zip(&resolved_from) {
        let _ = writeln!(t, "resistance-resolved-from\t{m}\t{first}");
    }
    for (i, &row) in grid_rows.iter().enumerate() {
        for (_, members) in &families {
            for (a, &r) in members.iter().enumerate() {
                for &s in &members[a..] {
                    let (re, im) = entry(&main, row, r, s);
                    let _ = writeln!(t, "z\t{i}\t{r}\t{s}\t{re:?}\t{im:?}");
                }
            }
        }
    }
    for (j, &row) in mid_rows.iter().enumerate() {
        for r in 0..modes.len() {
            let (re, im) = entry(&main, row, r, r);
            let _ = writeln!(t, "zmid\t{:?}\t{r}\t{re:?}\t{im:?}", mids[j]);
        }
    }
    std::fs::write(&out, &t).expect("write the table");
    eprintln!(
        "{} modes, {} ka samples, {} midpoints written to {out}",
        modes.len(),
        ka.len(),
        mids.len()
    );
    for note in &notes {
        eprintln!("{note}");
    }
    for ((m, members), first) in families.iter().zip(&resolved_from) {
        eprintln!(
            "family m = {m}: {} modes, resistance resolved from ka {:.4}",
            members.len(),
            ka[*first]
        );
    }
}
