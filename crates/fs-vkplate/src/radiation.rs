//! Exterior air loading of a drum head's non-volume (m >= 1) motion from a
//! baked, dimensionless modal radiation table (bead .370, ruling 6310).
//!
//! The table holds `zhat_rs(ka) = Z_rs / (rho c a^2)` for clamped-membrane
//! shapes on a unit disk in an infinite rigid baffle, computed with
//! `fs_bem::modal_loading` by `examples/bake_head_radiation.rs`. A head of
//! radius `a` in the current gas sees `Z_rs(omega) = rho c a^2 zhat_rs(omega a / c)`,
//! so density and sound speed always come from the instrument's gas. Only
//! `m >= 1` families load the head: the volume (`m = 0`) motion keeps the
//! kit builder's compact monopole, piston air mass and cavity, and its baked
//! family only states that compact model's error in the receipt. A head whose
//! own losses are modelled (film, viscous layer) extends the load to every
//! shape its grid resolves (`ModalAirLoad::new_with_losses`), so that no
//! resolved mode is left without loss.
//!
//! The time-domain load is narrowband. Each retained mode gets the diagonal
//! resistance `Re Z_rr` and air mass `-Im Z_rr / omega` at its own loaded
//! frequency, found by fixed-point iteration on the head's Rayleigh-quotient
//! stiffness and mass. The neglected coupling between modes of one family is
//! bounded to first order and reported, never silently dropped.
//!
//! Only arithmetic, `sqrt` and `fs_math::det`'s deterministic `ln`/`exp` are
//! used here (no platform libm), so native and WASM builds evaluate the same
//! shapes and loads.

use crate::GasState;
use std::fmt::Write as _;

/// The shipped table, baked by `examples/bake_head_radiation.rs`.
pub const SHIPPED_TABLE: &str = include_str!("../data/head-radiation-planar-baffled.tsv");

pub const TABLE_FORMAT: &str = "frankensim.vkplate-head-radiation.v1";
pub const ROUTE: &str = "planar-baffled-rayleigh";
pub const GEOMETRY: &str = "unit-disk-in-infinite-rigid-baffle";
/// Families `m >= 1` load the head; the `m = 0` family is baked only to state
/// the compact volume model's error (`R_compact / R_baffled`), never applied.
pub const MODE_SET: &str = "clamped-membrane-bessel";

/// `J_m(x)` for `x >= 0` by Miller's backward recurrence, normalized with
/// `1 = J_0 + 2 sum J_2k`. Accurate to about 1e-14 of the largest value for
/// the orders and arguments of the table (`x` below about 60).
pub fn bessel_j(order: usize, x: f64) -> f64 {
    if !(x.is_finite() && x >= 0.0) {
        return f64::NAN;
    }
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

/// The `index`-th positive zero of `J_order` (`index >= 1`), by bracketing
/// and bisection to adjacent floating-point values.
pub fn bessel_zero(order: usize, index: usize) -> f64 {
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

/// One clamped-membrane shape `J_m(j r) cos(m theta)` (or `sin`) on the unit
/// disk, scaled so its mean square over the disk is one: the integral of its
/// square over a disk of radius `a` is `pi a^2`.
#[derive(Clone, Debug, PartialEq)]
pub struct ModeShape {
    pub m: usize,
    pub n: usize,
    pub zero: f64,
    scale: f64,
}

impl ModeShape {
    pub fn new(m: usize, n: usize) -> ModeShape {
        let zero = bessel_zero(m, n);
        // integral_0^1 J_m(j r)^2 r dr = J_{m+1}(j)^2 / 2, and the angular
        // mean square of cos(m theta) is 1/2 for m >= 1 (1 for m = 0).
        let angular = if m == 0 { 1.0 } else { 2.0_f64.sqrt() };
        let scale = angular / bessel_j(m + 1, zero).abs();
        ModeShape { m, n, zero, scale }
    }

    /// Value at `(x, y)` in disk radii; zero on and outside the rim. An
    /// `m = 0` shape has no sine parity.
    pub fn value(&self, x: f64, y: f64, sine: bool) -> f64 {
        let r2 = x * x + y * y;
        if r2 >= 1.0 || (self.m == 0 && sine) {
            return 0.0;
        }
        if r2 == 0.0 {
            return if self.m == 0 { self.scale } else { 0.0 };
        }
        let r = r2.sqrt();
        // cos(m theta) + i sin(m theta) = ((x + i y) / r)^m, by multiplication.
        let (c, s) = (x / r, y / r);
        let (mut re, mut im) = (1.0, 0.0);
        for _ in 0..self.m {
            (re, im) = (re * c - im * s, re * s + im * c);
        }
        self.scale * bessel_j(self.m, self.zero * r) * if sine { im } else { re }
    }
}

/// A parsed table: family blocks of `zhat` at increasing `ka` samples.
#[derive(Clone, Debug)]
pub struct RadiationTable {
    modes: Vec<ModeShape>,
    ka: Vec<f64>,
    /// Per family `m`: the indices into `modes`, in table order.
    families: Vec<(usize, Vec<usize>)>,
    /// Per family: the first sample from which every diagonal resistance is
    /// above the bake's rounding floor (the multipole cancellation of a high
    /// `m` family at small `ka` leaves only roundoff below it).
    resolved_from: Vec<usize>,
    /// `[ka sample][family][r * size + s]` as baked: `(Re zhat, Im zhat)`.
    blocks: Vec<Vec<Vec<(f64, f64)>>>,
    pub identity: String,
    /// Largest measured interpolation error over all intervals and modes.
    pub interpolation_bound: f64,
    /// `[interval][mode]`, measured at parse from the bake's midpoints.
    interval_error: Vec<Vec<f64>>,
    pub discretization_bound: f64,
    pub content_hash: u64,
}

/// One interpolated entry. `resistance` is `None` below the family's
/// resolved range, where only its bound `resistance_bound` (the value at the
/// first resolved sample; low-ka resistance increases with ka) is known.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Zhat {
    pub resistance: Option<f64>,
    pub resistance_bound: f64,
    pub reactance: f64,
}

fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(*b)).wrapping_mul(0x100_0000_01b3)
    })
}

impl RadiationTable {
    pub fn shipped() -> Result<RadiationTable, String> {
        RadiationTable::parse(SHIPPED_TABLE)
    }

    /// Parse and verify a table. It refuses another format, route, geometry
    /// or mode family; a mode whose zero differs from this crate's own shape;
    /// and an incomplete, nonfinite or unordered sample grid.
    pub fn parse(text: &str) -> Result<RadiationTable, String> {
        let refuse = |what: &str| format!("head radiation table: {what}");
        let mut lines = text.lines();
        if lines.next() != Some(TABLE_FORMAT) {
            return Err(refuse("not a frankensim.vkplate-head-radiation.v1 table"));
        }
        let (mut route, mut geometry, mut mode_set, mut identity) = (None, None, None, Vec::new());
        let mut discretization = None;
        let mut modes: Vec<ModeShape> = Vec::new();
        let mut ka: Vec<f64> = Vec::new();
        let mut entries: Vec<(usize, usize, usize, f64, f64)> = Vec::new();
        let mut resolved: Vec<(usize, usize)> = Vec::new();
        let mut probes: Vec<(usize, f64, f64, f64)> = Vec::new();
        let mut midpoints: Vec<(f64, usize, f64, f64)> = Vec::new();
        let number = |s: &str| s.parse::<f64>().map_err(|_| refuse("malformed number"));
        let index = |s: &str| s.parse::<usize>().map_err(|_| refuse("malformed index"));
        for line in lines {
            let cols: Vec<&str> = line.split('\t').collect();
            match cols[0] {
                "route" if cols.len() == 2 => route = Some(cols[1]),
                "geometry" if cols.len() >= 2 => {
                    geometry = Some(cols[1]);
                    identity.push(line.to_string());
                }
                "mode-set" if cols.len() >= 2 => {
                    mode_set = Some(cols[1]);
                    identity.push(line.to_string());
                }
                "declaration" | "note" => identity.push(line.to_string()),
                "discretization-bound" if cols.len() == 2 => {
                    discretization = Some(number(cols[1])?)
                }
                "probe" if cols.len() == 5 => {
                    probes.push((
                        index(cols[1])?,
                        number(cols[2])?,
                        number(cols[3])?,
                        number(cols[4])?,
                    ));
                }
                "zmid" if cols.len() == 5 => {
                    midpoints.push((
                        number(cols[1])?,
                        index(cols[2])?,
                        number(cols[3])?,
                        number(cols[4])?,
                    ));
                }
                "mode" if cols.len() == 4 => {
                    let (m, n, zero) = (index(cols[1])?, index(cols[2])?, number(cols[3])?);
                    if n == 0 || modes.iter().any(|s| s.m == m && s.n == n) {
                        return Err(refuse("modes start at n = 1 and are unique"));
                    }
                    let shape = ModeShape::new(m, n);
                    if (shape.zero - zero).abs() > 1e-9 * zero {
                        return Err(refuse(
                            "a mode's Bessel zero differs from this crate's shape",
                        ));
                    }
                    modes.push(shape);
                }
                "ka" if cols.len() == 2 => {
                    let value = number(cols[1])?;
                    if !(value.is_finite() && value > 0.0)
                        || ka.last().is_some_and(|&last| value <= last)
                    {
                        return Err(refuse("ka samples must be positive and increasing"));
                    }
                    ka.push(value);
                }
                "resistance-resolved-from" if cols.len() == 3 => {
                    resolved.push((index(cols[1])?, index(cols[2])?))
                }
                "z" if cols.len() == 6 => entries.push((
                    index(cols[1])?,
                    index(cols[2])?,
                    index(cols[3])?,
                    number(cols[4])?,
                    number(cols[5])?,
                )),
                _ => return Err(refuse("unknown or malformed row")),
            }
        }
        if route != Some(ROUTE) || geometry != Some(GEOMETRY) || mode_set != Some(MODE_SET) {
            return Err(refuse(
                "route, geometry or mode set is not the one this crate applies",
            ));
        }
        let Some(discretization_bound) = discretization else {
            return Err(refuse("missing discretization bound"));
        };
        if modes.is_empty() || ka.len() < 2 {
            return Err(refuse("no modes or fewer than two ka samples"));
        }
        // The bake's own shape evaluations must be this crate's shapes.
        for (i, shape) in modes.iter().enumerate() {
            let checked: Vec<_> = probes.iter().filter(|p| p.0 == i).collect();
            if checked.is_empty() {
                return Err(refuse("a mode has no probe values"));
            }
            for &&(_, x, y, value) in &checked {
                if (shape.value(x, y, false) - value).abs() > 1e-12 * value.abs().max(1.0) {
                    return Err(refuse("a probe value differs from this crate's shape"));
                }
            }
        }
        if probes.iter().any(|p| p.0 >= modes.len()) || midpoints.is_empty() {
            return Err(refuse(
                "a probe names an undeclared mode, or no midpoint samples",
            ));
        }
        let mut families: Vec<(usize, Vec<usize>)> = Vec::new();
        for (i, shape) in modes.iter().enumerate() {
            match families.iter_mut().find(|(m, _)| *m == shape.m) {
                Some((_, members)) => members.push(i),
                None => families.push((shape.m, vec![i])),
            }
        }
        let mut blocks: Vec<Vec<Vec<Option<(f64, f64)>>>> = ka
            .iter()
            .map(|_| {
                families
                    .iter()
                    .map(|(_, f)| vec![None; f.len() * f.len()])
                    .collect()
            })
            .collect();
        for (sample, r, s, re, im) in entries {
            let (Some(shape_r), Some(shape_s)) = (modes.get(r), modes.get(s)) else {
                return Err(refuse("entry names an undeclared mode"));
            };
            if sample >= ka.len() || shape_r.m != shape_s.m || !re.is_finite() || !im.is_finite() {
                return Err(refuse(
                    "entry outside the sample grid, across families, or nonfinite",
                ));
            }
            let f = families
                .iter()
                .position(|(m, _)| *m == shape_r.m)
                .unwrap_or(0);
            let members = &families[f].1;
            let size = members.len();
            let (a, b) = (
                members.iter().position(|&i| i == r).unwrap_or(0),
                members.iter().position(|&i| i == s).unwrap_or(0),
            );
            for (x, y) in [(a, b), (b, a)] {
                if blocks[sample][f][x * size + y].is_some_and(|v| v != (re, im)) {
                    return Err(refuse("a family block is not symmetric"));
                }
                blocks[sample][f][x * size + y] = Some((re, im));
            }
        }
        let mut resolved_from = vec![usize::MAX; families.len()];
        for (m, first) in resolved {
            let Some(f) = families.iter().position(|(fm, _)| *fm == m) else {
                return Err(refuse("a resolution row names an undeclared family"));
            };
            if first >= ka.len() || resolved_from[f] != usize::MAX {
                return Err(refuse("a resolution row is outside the grid or repeated"));
            }
            resolved_from[f] = first;
        }
        if resolved_from.contains(&usize::MAX) {
            return Err(refuse(
                "every family needs its resistance-resolved-from row",
            ));
        }
        let blocks = blocks
            .into_iter()
            .map(|families| {
                families
                    .into_iter()
                    .map(|block| block.into_iter().collect::<Option<Vec<_>>>())
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| refuse("a family block is incomplete"))?;
        for (f, (_, members)) in families.iter().enumerate() {
            let size = members.len();
            for sample in &blocks[resolved_from[f]..] {
                if (0..size).any(|a| sample[f][a * size + a].0 <= 0.0) {
                    return Err(refuse("a resolved diagonal resistance is not positive"));
                }
            }
        }
        let mut table = RadiationTable {
            modes,
            ka,
            families,
            resolved_from,
            blocks,
            identity: identity.join(" | "),
            interpolation_bound: 0.0,
            interval_error: Vec::new(),
            discretization_bound,
            content_hash: fnv1a(text.as_bytes()),
        };
        // This crate's interpolation rule against the bake at midpoints, per
        // interval and mode: resolved resistance relative to itself, reactance
        // relative to |Z| (it crosses zero where a shape turns supersonic).
        // Spans `zhat` refuses are not measured; they are never used.
        let intervals = table.ka.len() - 1;
        table.interval_error = vec![vec![f64::INFINITY; table.modes.len()]; intervals];
        for (k, r, re, im) in midpoints {
            if r >= table.modes.len() || !(re.is_finite() && im.is_finite()) {
                return Err(refuse("a midpoint sample is undeclared or nonfinite"));
            }
            let interval = table
                .interval(k)
                .ok_or_else(|| refuse("a midpoint lies outside the ka grid"))?;
            let Ok(z) = table.zhat(r, r, k) else { continue };
            let mut error = (z.reactance - im).abs() / re.hypot(im);
            if let Some(resistance) = z.resistance {
                error = error.max((resistance - re).abs() / re.abs());
            }
            table.interval_error[interval][r] = error;
        }
        table.interpolation_bound = table
            .interval_error
            .iter()
            .flatten()
            .filter(|e| e.is_finite())
            .fold(0.0, |a, &e| a.max(e));
        Ok(table)
    }

    /// The sample interval containing `ka`, if inside the grid.
    fn interval(&self, ka: f64) -> Option<usize> {
        let (lo, hi) = self.ka_range();
        (ka >= lo && ka <= hi).then(|| {
            self.ka
                .partition_point(|&x| x < ka)
                .clamp(1, self.ka.len() - 1)
                - 1
        })
    }

    /// The measured interpolation error of `mode`'s diagonal in the interval
    /// containing `ka` (infinite where unmeasured or outside the grid).
    pub fn interpolation_error(&self, mode: usize, ka: f64) -> f64 {
        self.interval(ka)
            .map_or(f64::INFINITY, |i| self.interval_error[i][mode])
    }

    pub fn modes(&self) -> &[ModeShape] {
        &self.modes
    }
    pub fn ka_range(&self) -> (f64, f64) {
        (self.ka[0], self.ka[self.ka.len() - 1])
    }

    /// `zhat_rs(ka)` for two modes of one family; zero across families.
    /// Between samples it is linear in `ln ka`: of `ln Re` and `ln(-Im)` on
    /// the diagonal (exact for the low-ka power laws), and of the entry
    /// normalized by the geometric mean of its two diagonals off it. Outside
    /// the sampled range it refuses; it never clamps or extrapolates.
    pub fn zhat(&self, r: usize, s: usize, ka: f64) -> Result<Zhat, String> {
        use fs_math::det::{exp, ln};
        let (lo, hi) = self.ka_range();
        if !(ka >= lo && ka <= hi) {
            return Err(format!(
                "head radiation: ka {ka} is outside the table's [{lo}, {hi}]"
            ));
        }
        let m = self.modes[r].m;
        if self.modes[s].m != m {
            return Ok(Zhat {
                resistance: Some(0.0),
                resistance_bound: 0.0,
                reactance: 0.0,
            });
        }
        let f = self
            .families
            .iter()
            .position(|(fm, _)| *fm == m)
            .unwrap_or(0);
        let members = &self.families[f].1;
        let size = members.len();
        let a = members.iter().position(|&i| i == r).unwrap_or(0);
        let b = members.iter().position(|&i| i == s).unwrap_or(0);
        let upper = self
            .ka
            .partition_point(|&x| x < ka)
            .clamp(1, self.ka.len() - 1);
        let t = (ln(ka) - ln(self.ka[upper - 1])) / (ln(self.ka[upper]) - ln(self.ka[upper - 1]));
        let at = |sample: usize, x: usize, y: usize| self.blocks[sample][f][x * size + y];
        let log_lerp = |v0: f64, v1: f64| exp(ln(v0) + t * (ln(v1) - ln(v0)));
        let diagonal = |x: usize| -> Result<(Option<f64>, f64, f64), String> {
            let (r0, r1) = (at(upper - 1, x, x), at(upper, x, x));
            if r0.1 >= 0.0 || r1.1 >= 0.0 {
                let shape = &self.modes[members[x]];
                return Err(format!(
                    "head radiation: the ({},{}) reactance is not mass-like near ka {ka} (a supersonic shape); refused",
                    shape.m, shape.n
                ));
            }
            let resolved = upper - 1 >= self.resolved_from[f];
            let resistance = resolved.then(|| log_lerp(r0.0, r1.0));
            let bound = at(self.resolved_from[f], x, x).0;
            Ok((resistance, bound, -log_lerp(-r0.1, -r1.1)))
        };
        let (ra, bound_a, xa) = diagonal(a)?;
        if a == b {
            return Ok(Zhat {
                resistance: ra,
                resistance_bound: bound_a,
                reactance: xa,
            });
        }
        let (rb, bound_b, xb) = diagonal(b)?;
        let normalized = |sample: usize, part: fn((f64, f64)) -> f64| {
            part(at(sample, a, b))
                / (part(at(sample, a, a)) * part(at(sample, b, b)))
                    .abs()
                    .sqrt()
        };
        let lerp = |part: fn((f64, f64)) -> f64| {
            let (c0, c1) = (normalized(upper - 1, part), normalized(upper, part));
            c0 + t * (c1 - c0)
        };
        let resistance = match (ra, rb) {
            (Some(ra), Some(rb)) => Some(lerp(|v| v.0) * (ra * rb).sqrt()),
            _ => None,
        };
        Ok(Zhat {
            resistance,
            resistance_bound: (bound_a * bound_b).sqrt(),
            reactance: lerp(|v| v.1) * (xa * xb).sqrt(),
        })
    }
}

/// One head's retained narrowband modal load, fixed for one gas and grid.
///
/// With area-weighted shapes `U = W Phi`, Gram `G = Phi' W Phi`, retained
/// resistances `R` and air masses `Ma` (diagonal), the generalized velocity
/// is `G^-1 U' v`, so the nodal damping is `U B U'` with `B = G^-1 R G^-1`
/// and the nodal added mass `U A U'` with `A = G^-1 Ma G^-1`.
#[derive(Clone, Debug)]
pub struct ModalAirLoad {
    /// Area-weighted shapes, one per retained mode and parity: `h^2 phi(node)`.
    pub(crate) cols: Vec<Vec<f64>>,
    gram: Vec<f64>,
    air: Vec<f64>,
    /// `B` and `A` (row-major, `cols.len()` square).
    pub(crate) damping: Vec<f64>,
    pub(crate) mass: Vec<f64>,
    /// `trace(B U' M^-1 U)`, an upper bound on the damping's largest
    /// `M^-1 C` eigenvalue for the explicit-step admission.
    pub(crate) rate_bound: f64,
    /// With a film loss, `damping` is radiation plus film and this is the
    /// film part alone, so `radiated` can say how much of a step was film.
    film_damping: Option<Vec<f64>>,
    film_step_j: f64,
    film_dissipated_j: f64,
    /// Woodbury data for the implicit air mass at the step it was prepared for.
    pub(crate) prepared: Option<Prepared>,
    /// Per-step buffers: nodal right-hand side and correction, projections.
    pub(crate) rhs: Vec<f64>,
    pub(crate) corr: Vec<f64>,
    pub(crate) p_wp: Vec<f64>,
    pub(crate) y: Vec<f64>,
    pub(crate) b: Vec<f64>,
    pub(crate) w2: Vec<f64>,
    za: Vec<f64>,
    zb: Vec<f64>,
    /// One line per retained mode, then the summary and stated scope.
    pub receipt: Vec<String>,
    /// Largest first-order share of the dropped in-family coupling, relative
    /// to the retained resistance and to the retained air mass.
    pub dropped_loss_share: f64,
    pub dropped_mass_share: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct Prepared {
    pub(crate) k: f64,
    /// Cholesky factor (lower, row-major) of `C^-1 + V' D^-1 V` over the
    /// columns `[u?, modal...]`, with `C^-1 = k^2 A^-1 = k^2 G Ma^-1 G` for
    /// the modal block and `1 / mu` for the volume column.
    pub(crate) chol: Vec<f64>,
    pub(crate) with_volume: bool,
}

/// Mode `(m, n, parity)` retained on a head, with its narrowband load.
struct Retained {
    shape: usize,
    sine: bool,
    stiffness: f64,
    mass: f64,
    omega: f64,
    ka: f64,
    resistance: f64,
    /// Below the table's resolved range: no resistance is applied, and the
    /// true one is under this bound.
    below_resolution: Option<f64>,
    air_mass: f64,
    iterations: usize,
    /// MAC against the family's Rayleigh-Ritz estimate of the FD mode.
    mac: f64,
    /// `phi' K_bend phi / phi' K phi`: the share of the mode's stiffness
    /// energy held in bending, from the body's own operators.
    bending_share: f64,
    /// The film's loss tangent at the loaded frequency; `None` without a
    /// film loss, or outside the film source's band.
    film_tan_delta: Option<f64>,
    /// Film loss as a modal damping coefficient,
    /// `tan delta * phi' K_bend phi / omega` at the loaded frequency.
    film_resistance: f64,
    /// Stokes face-layer resistance at the loaded frequency.
    viscous_resistance: f64,
    /// `shape` indexes the table's modes; false: a shape beyond the table's
    /// set (evanescent air mass, no radiation resistance).
    in_table: bool,
    /// A volume shape's share already on the compact net-volume model,
    /// `(R_c V^2, piston air mass V^2 / A^2)`, not in `resistance`/`air_mass`.
    compact: Option<(f64, f64)>,
    /// A table shape's evanescent air mass over its table air mass, minus 1.
    evanescent_check: Option<f64>,
}

fn cholesky(a: &mut [f64], n: usize) -> Result<(), String> {
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j {
            d -= a[j * n + k] * a[j * n + k];
        }
        if !(d > 0.0 && d.is_finite()) {
            return Err("head radiation: the air-mass system is not positive definite".into());
        }
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut v = a[i * n + j];
            for k in 0..j {
                v -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = v / d;
        }
        for i in 0..j {
            a[i * n + j] = 0.0;
        }
    }
    Ok(())
}

pub(crate) fn cholesky_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for i in 0..n {
        let mut v = b[i];
        for k in 0..i {
            v -= l[i * n + k] * b[k];
        }
        b[i] = v / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut v = b[i];
        for k in i + 1..n {
            v -= l[k * n + i] * b[k];
        }
        b[i] = v / l[i * n + i];
    }
}

fn invert_spd(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut l = a.to_vec();
    cholesky(&mut l, n)?;
    let mut inverse = vec![0.0; n * n];
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        cholesky_solve(&l, n, &mut e);
        for i in 0..n {
            inverse[i * n + j] = e[i];
        }
    }
    Ok(inverse)
}

/// Shapes are retained when the grid resolves them: at least six nodes per
/// radial wavelength and per azimuthal wavelength at the rim. With spacing
/// `h` on a head of radius `a` that is `j_mn <= 2 pi a / (6 h)` and
/// `m <= 2 pi a / (6 h)`.
pub const NODES_PER_WAVELENGTH: f64 = 6.0;

/// The loaded-frequency fixed point refuses beyond this many iterations.
pub const MAX_LOADED_ITERATIONS: usize = 100;

/// A shape whose MAC against its family's Rayleigh-Ritz FD mode is below
/// this carries a caveat in the receipt.
pub const LOW_MAC: f64 = 0.9;

/// Eigenpairs of a small symmetric matrix by cyclic Jacobi rotations
/// (columns of the returned row-major matrix are the eigenvectors).
fn jacobi_eigen(a: &mut [f64], n: usize) -> Vec<f64> {
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    let scale = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum::<f64>()
            .sqrt();
        if off <= 1e-15 * scale {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    v
}

/// MAC of each trial shape against the FD head's own modes as estimated by
/// Rayleigh-Ritz on the family's trials and one Krylov step `M^-1 K phi`
/// (mass-orthonormalized; nearly dependent directions are dropped).
fn ritz_macs(
    phis: &[&Vec<f64>],
    kphis: &[&Vec<f64>],
    mass: &[f64],
    stiffness: &dyn Fn(&[f64]) -> Vec<f64>,
) -> Vec<f64> {
    let m_dot = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .zip(mass)
            .map(|((x, y), m)| x * y * m)
            .sum::<f64>()
    };
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let candidates = phis.iter().map(|p| (*p).clone()).chain(
        kphis
            .iter()
            .map(|kp| kp.iter().zip(mass).map(|(k, m)| k / m).collect()),
    );
    for mut v in candidates {
        let norm0 = m_dot(&v, &v).sqrt();
        for _ in 0..2 {
            for q in &basis {
                let c = m_dot(q, &v);
                for (x, y) in v.iter_mut().zip(q) {
                    *x -= c * y;
                }
            }
        }
        let norm = m_dot(&v, &v).sqrt();
        if norm > 1e-8 * norm0 {
            basis.push(v.iter().map(|x| x / norm).collect());
        }
    }
    let b = basis.len();
    let images: Vec<Vec<f64>> = basis.iter().map(|q| stiffness(q)).collect();
    let mut a = vec![0.0; b * b];
    for i in 0..b {
        for j in 0..b {
            a[i * b + j] =
                0.5 * (crate::dot(&basis[i], &images[j]) + crate::dot(&basis[j], &images[i]));
        }
    }
    let vectors = jacobi_eigen(&mut a, b);
    phis.iter()
        .map(|phi| {
            let coefficients: Vec<f64> = basis.iter().map(|q| m_dot(q, phi)).collect();
            let norm = m_dot(phi, phi);
            (0..b)
                .map(|k| {
                    let c: f64 = (0..b).map(|i| coefficients[i] * vectors[i * b + k]).sum();
                    c * c / norm
                })
                .fold(0.0, f64::max)
        })
        .collect()
}

/// A head film's loss tangent at a mode's loaded frequency [Hz]: `Ok(None)`
/// outside the source's band (that mode then carries no film loss, and the
/// receipt names it), an error when the source cannot answer.
pub type FilmTanDelta<'a> = &'a dyn Fn(f64) -> Result<Option<f64>, String>;

impl ModalAirLoad {

    /// Build the load of a clamped full circular head of `radius` whose nodes
    /// are `xy` (centred) with spacing `h`, nodal `mass` and linear
    /// `stiffness(v) -> K v`, in `gas`. Every `m >= 1` table mode the grid
    /// resolves is retained in both parities; the rest are named in the
    /// receipt as upper bounds. `compact_volume_resistance` is the volume
    /// model's resistance on net volume velocity [Pa s/m^3], compared with
    /// the baked `m = 0` self-resistance for the receipt only.
    pub fn new(
        table: &RadiationTable,
        radius: f64,
        h: f64,
        xy: &[[f64; 2]],
        mass: &[f64],
        stiffness: &dyn Fn(&[f64]) -> Vec<f64>,
        gas: &GasState,
        compact_volume_resistance: Option<f64>,
    ) -> Result<ModalAirLoad, String> {
        Self::new_with_film(
            table,
            radius,
            h,
            xy,
            mass,
            stiffness,
            gas,
            compact_volume_resistance,
            None,
        )
    }

    /// As `new`, with an optional film loss `(tan delta, K_bend)` whose loss
    /// tangent is held at one value for every mode.
    ///
    /// Each retained mode gets the hysteretic loss of its bending energy,
    /// `c = tan delta * phi' K_bend phi / omega` at its own loaded frequency
    /// (`eta = tan delta * phi' K_bend phi / phi' K phi`), applied through the
    /// same modal damping as the radiation. No fit: each retained mode's loss
    /// is exact for its shape. The volume (m = 0) shapes and every shape the
    /// table lacks or the grid cannot resolve carry no film loss; the receipt
    /// names the covered band. `K_bend` is the body's bending operator alone,
    /// without the tension: tension energy is pre-stress work, not strain.
    /// `new_with_losses` takes the loss tangent at each mode's frequency.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_film(
        table: &RadiationTable,
        radius: f64,
        h: f64,
        xy: &[[f64; 2]],
        mass: &[f64],
        stiffness: &dyn Fn(&[f64]) -> Vec<f64>,
        gas: &GasState,
        compact_volume_resistance: Option<f64>,
        film: Option<(f64, &dyn Fn(&[f64]) -> Vec<f64>)>,
    ) -> Result<ModalAirLoad, String> {
        let held = |_: f64| Ok(film.map(|(tan_delta, _)| tan_delta));
        let film = film.map(|(_, bending)| (&held as FilmTanDelta, bending));
        Self::new_with_losses(table, radius, h, xy, mass, stiffness, gas, compact_volume_resistance, film, None, None)
    }

    /// The general constructor: an optional film loss whose loss tangent is
    /// read at each retained mode's own loaded frequency, and an optional
    /// Stokes face layer `(areal density, coincidence frequency)`.
    ///
    /// `extend` (the piston air mass over the head's area squared [kg/m^4],
    /// 0 without one) retains every shape the grid resolves, not only the
    /// table's `m >= 1` set, so that the film and viscous losses reach every
    /// resolved mode. A volume shape `(0,n)`, `n >= 2`, takes the table's
    /// baffled load less what the compact net-volume model already gives it
    /// (`compact_volume_resistance` and the piston, times its volume
    /// squared); `(0,1)` stays on the compact model alone. A shape beyond the
    /// table's set takes the air mass of its untruncated evanescent field,
    /// `rho A / sqrt(k_s^2 - k^2)` (A its area integral of phi^2, k_s its
    /// Bessel zero over the radius), and no radiation resistance: its decay
    /// is an upper bound.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new_with_losses(
        table: &RadiationTable,
        radius: f64,
        h: f64,
        xy: &[[f64; 2]],
        mass: &[f64],
        stiffness: &dyn Fn(&[f64]) -> Vec<f64>,
        gas: &GasState,
        compact_volume_resistance: Option<f64>,
        film: Option<(FilmTanDelta, &dyn Fn(&[f64]) -> Vec<f64>)>,
        viscous: Option<(f64, f64)>,
        extend: Option<f64>,
    ) -> Result<ModalAirLoad, String> {
        if let Some((areal_density, coincidence_hz)) = viscous
            && !(areal_density.is_finite() && areal_density > 0.0 && coincidence_hz.is_finite() && coincidence_hz > 0.0)
        {
            return Err("head viscous layer needs positive areal density and coincidence frequency".into());
        }
        if !(radius > 0.0 && h > 0.0 && radius.is_finite() && h.is_finite())
            || xy.len() != mass.len()
        {
            return Err("head radiation: invalid head geometry".into());
        }
        let (rho, c) = (gas.density, gas.sound_speed);
        let scale = rho * c * radius * radius;
        let limit = 2.0 * std::f64::consts::PI * radius / (NODES_PER_WAVELENGTH * h);
        let h2 = h * h;
        let sample = |shape: &ModeShape, sine: bool| -> Vec<f64> {
            xy.iter()
                .map(|p| shape.value(p[0] / radius, p[1] / radius, sine))
                .collect()
        };
        let quotient = |phi: &[f64], k_phi: &[f64]| -> (f64, f64) {
            (
                phi.iter().zip(k_phi).map(|(a, b)| a * b).sum(),
                phi.iter().zip(mass).map(|(a, m)| a * a * m).sum(),
            )
        };
        let mut retained: Vec<Retained> = Vec::new();
        let (mut phis, mut kphis): (Vec<Vec<f64>>, Vec<Vec<f64>>) = (Vec::new(), Vec::new());
        let mut dropped = Vec::new();
        // In-vacuo frequency [Hz] and label of every shape the table refused.
        let mut refused_shapes: Vec<(f64, String)> = Vec::new();
        let mut volume_lines = Vec::new();
        // The table's shapes, then (extended) every shape the grid resolves
        // that the table's set lacks.
        let tabulated = table.modes().len();
        let mut shapes: Vec<ModeShape> = table.modes().to_vec();
        if extend.is_some() {
            for m in 0..=(limit as usize) {
                for n in 1.. {
                    if table.modes().iter().any(|s| s.m == m && s.n == n) {
                        continue;
                    }
                    let shape = ModeShape::new(m, n);
                    if shape.zero > limit {
                        break;
                    }
                    shapes.push(shape);
                }
            }
        }
        let (compact_r, piston) = (compact_volume_resistance.unwrap_or(0.0), extend.unwrap_or(0.0));
        for (index, shape) in shapes.iter().enumerate() {
            let in_table = index < tabulated;
            let label = format!("({},{})", shape.m, shape.n);
            let unresolved = shape.zero > limit || shape.m as f64 > limit;
            if unresolved && shape.m == 0 {
                volume_lines.push(format!(
                    "volume shape {label}: not resolved by this grid, no comparison"
                ));
                continue;
            }
            // The receipt's main claim is the (1,1) shape's load: without it
            // the drum refuses whole.
            if unresolved && shape.m == 1 && shape.n == 1 {
                return Err(format!(
                    "head radiation: this grid cannot resolve the (1,1) shape (j_mn {:.3} > {limit:.3}); the drum refuses",
                    shape.zero
                ));
            }
            if unresolved {
                dropped.push(format!(
                    "{label} dropped: j_mn = {:.3} exceeds the grid's {limit:.3}; NO exterior load (no radiation loss, no air mass), so its T60 is an upper bound, not a computed decay",
                    shape.zero
                ));
                continue;
            }
            if shape.m == 0 && in_table {
                // The compact route's error at this ka: a diagnostic, and
                // (extended) the share subtracted from the shape's own load.
                let applied = extend.is_some() && shape.n >= 2;
                let phi = sample(shape, false);
                let (stiff, modal_mass) = quotient(&phi, &stiffness(&phi));
                let omega = (stiff / modal_mass).sqrt();
                let ka = omega * radius / c;
                let volume: f64 = phi.iter().sum::<f64>() * h2;
                let baffled = table
                    .zhat(index, index, ka)
                    .map(|z| (z.resistance.map(|r| scale * r), scale * z.resistance_bound));
                let line = match (compact_volume_resistance, baffled) {
                    (Some(compact), Ok((Some(r), _))) => format!(
                        "{label}: ka {ka:.4} in vacuo; R_compact / R_baffled = {:.3} (compact: the kit's net-volume resistance times this shape's volume squared; baffled: Re Z at this frequency); diagnostic, not applied",
                        compact * volume * volume / r
                    ),
                    (Some(compact), Ok((None, bound))) => format!(
                        "{label}: ka {ka:.4} in vacuo; R_baffled below the table's resolution (< {bound:.3e} kg/s), R_compact {:.3e} kg/s; diagnostic, not applied",
                        compact * volume * volume
                    ),
                    (_, Err(e)) => format!("{label}: no baffled comparison ({e})"),
                    (None, _) => format!(
                        "{label}: ka {ka:.4} in vacuo; this head has no compact volume resistance to compare"
                    ),
                };
                let line = if applied {
                    line.replace("diagnostic, not applied", "applied: the shape takes the baffled load less this compact share")
                } else {
                    line
                };
                volume_lines.push(format!("volume shape {line}"));
                if !applied {
                    continue;
                }
            }
            // The area integral of phi^2 and wavenumber of a shape beyond the
            // table: its untruncated field decays as exp(-sqrt(k_s^2 - k^2) z).
            let k_shape = shape.zero / radius;
            let evanescent = |phi: &[f64], omega: f64| -> Result<f64, String> {
                let kz2 = k_shape * k_shape - (omega / c) * (omega / c);
                if !(kz2 > 0.0) {
                    return Err(format!("the shape is supersonic at {omega:.1} rad/s, no evanescent air mass"));
                }
                Ok(rho * phi.iter().map(|v| v * v).sum::<f64>() * h2 / kz2.sqrt())
            };
            let parities: &[bool] = if shape.m == 0 { &[false] } else { &[false, true] };
            for &sine in parities {
                let phi = sample(shape, sine);
                let k_phi = stiffness(&phi);
                let (stiff, modal_mass) = quotient(&phi, &k_phi);
                if !(stiff > 0.0 && modal_mass > 0.0) {
                    return Err(format!(
                        "head radiation: mode {label} has no positive stiffness and mass on this grid"
                    ));
                }
                // The shape's exterior load at omega: its total air mass and,
                // in the table, the table's entry.
                let air_at = |omega: f64| -> Result<(f64, Option<Zhat>), String> {
                    if in_table {
                        let z = table.zhat(index, index, omega * radius / c)?;
                        Ok((-scale * z.reactance / omega, Some(z)))
                    } else {
                        Ok((evanescent(&phi, omega)?, None))
                    }
                };
                // Loaded frequency: omega^2 (M + Ma(omega)) = K, by fixed point.
                // A table refusal along the way (outside its ka range, or a
                // supersonic non-mass-like span) refuses THIS shape: it gets
                // no exterior load and is named as an upper bound. Nothing is
                // clamped or extrapolated.
                let mut omega = (stiff / modal_mass).sqrt();
                let mut iterations = 0;
                let refused = loop {
                    iterations += 1;
                    let air_mass = match air_at(omega) {
                        Ok((air_mass, _)) => air_mass,
                        Err(e) => break Some(e),
                    };
                    let next = (stiff / (modal_mass + air_mass)).sqrt();
                    let done = (next - omega).abs() <= 1e-13 * omega;
                    omega = next;
                    if done {
                        break None;
                    }
                    if iterations == MAX_LOADED_ITERATIONS {
                        return Err(format!(
                            "head radiation: mode {label} loaded frequency did not converge in {MAX_LOADED_ITERATIONS} iterations"
                        ));
                    }
                };
                let ka = omega * radius / c;
                let (air_total, z) = match refused.map_or_else(|| air_at(omega), Err) {
                    Ok(load) => load,
                    Err(e) => {
                        let parity = if sine { "sin" } else { "cos" };
                        if shape.m == 1 && shape.n == 1 {
                            return Err(format!(
                                "head radiation: the (1,1) {parity} shape is refused ({e}); the drum refuses"
                            ));
                        }
                        let hz = (stiff / modal_mass).sqrt() / (2.0 * std::f64::consts::PI);
                        refused_shapes.push((hz, format!("{label} {parity}")));
                        dropped.push(format!(
                            "{label} {parity} refused for this shape at {hz:.1} Hz in vacuo ({e}); NO exterior load (no radiation loss, no air mass), so its T60 is an upper bound, not a computed decay"
                        ));
                        continue;
                    }
                };
                let mut air_mass = air_total;
                let (mut resistance, below_resolution) = match z.as_ref().map(|z| z.resistance) {
                    Some(Some(r)) if r >= 0.0 => (scale * r, None),
                    Some(Some(_)) => {
                        return Err(format!(
                            "head radiation: mode {label} has an active load at ka {ka}"
                        ));
                    }
                    Some(None) => (0.0, z.as_ref().map(|z| scale * z.resistance_bound)),
                    None => (0.0, None),
                };
                // A table shape's own air mass against its untruncated
                // evanescent field (the law applied beyond the table).
                let evanescent_check = (extend.is_some() && in_table)
                    .then(|| evanescent(&phi, omega).ok().map(|e| e / air_total - 1.0))
                    .flatten();
                // A volume shape's net volume already carries the compact
                // monopole's resistance and the piston's air mass: apply the
                // rest, so that its own load is the baffled table's (beyond
                // the table, the evanescent air mass; the monopole stays).
                let compact = if shape.m == 0 {
                    let volume: f64 = phi.iter().sum::<f64>() * h2;
                    let (r_c, m_c) = (compact_r * volume * volume, piston * volume * volume);
                    let r_off = if z.is_some() && below_resolution.is_none() { r_c } else { 0.0 };
                    if !(resistance >= r_off && air_mass > m_c) {
                        dropped.push(format!(
                            "{label} kept on the compact model alone: its compact share (R {r_c:.3e} kg/s, air mass {m_c:.3e} kg) exceeds its own load (R {resistance:.3e} kg/s, air mass {air_mass:.3e} kg)"
                        ));
                        continue;
                    }
                    resistance -= r_off;
                    air_mass -= m_c;
                    Some((r_c, m_c))
                } else {
                    None
                };
                let (bending_share, film_tan_delta, film_resistance) = match film {
                    Some((tan_delta, bending)) => {
                        let bend: f64 = phi.iter().zip(&bending(&phi)).map(|(a, b)| a * b).sum();
                        if !(bend >= 0.0 && bend <= stiff * (1.0 + 1e-9)) {
                            return Err(format!(
                                "head film loss: mode {label} has bending energy {bend:e} outside [0, {stiff:e}]"
                            ));
                        }
                        let hz = omega / (2.0 * std::f64::consts::PI);
                        match tan_delta(hz).map_err(|e| format!("head film loss: mode {label} at {hz:.1} Hz: {e}"))? {
                            Some(t) if t.is_finite() && t > 0.0 && t < 1.0 => (bend / stiff, Some(t), t * bend / omega),
                            Some(t) => {
                                return Err(format!(
                                    "head film loss: mode {label} tan delta {t} must be finite and in (0, 1)"
                                ));
                            }
                            None => (bend / stiff, None, 0.0),
                        }
                    }
                    None => (0.0, None, 0.0),
                };
                let viscous_resistance = match viscous {
                    Some((areal_density, coincidence_hz)) => {
                        let layer = crate::air_loss::plate_viscous_layer(
                            areal_density,
                            omega,
                            omega / (2.0 * std::f64::consts::PI * coincidence_hz),
                            2,
                            gas,
                        )?;
                        layer.loss_factor * modal_mass * omega
                    }
                    None => 0.0,
                };
                retained.push(Retained {
                    shape: index,
                    sine,
                    stiffness: stiff,
                    mass: modal_mass,
                    omega,
                    ka,
                    resistance,
                    below_resolution,
                    air_mass,
                    iterations,
                    mac: 1.0,
                    bending_share,
                    film_tan_delta,
                    film_resistance,
                    viscous_resistance,
                    in_table,
                    compact,
                    evanescent_check,
                });
                phis.push(phi);
                kphis.push(k_phi);
            }
        }
        let n = retained.len();
        if n == 0 {
            return Err(
                "head radiation: the grid resolves none of the table's m >= 1 modes".into(),
            );
        }
        // MAC per (family, parity) against Rayleigh-Ritz FD modes.
        let mut groups: Vec<(usize, bool)> = retained
            .iter()
            .map(|r| (shapes[r.shape].m, r.sine))
            .collect();
        groups.sort_unstable();
        groups.dedup();
        for (m, sine) in groups {
            let members: Vec<usize> = (0..n)
                .filter(|&i| shapes[retained[i].shape].m == m && retained[i].sine == sine)
                .collect();
            let macs = ritz_macs(
                &members.iter().map(|&i| &phis[i]).collect::<Vec<_>>(),
                &members.iter().map(|&i| &kphis[i]).collect::<Vec<_>>(),
                mass,
                stiffness,
            );
            for (&i, mac) in members.iter().zip(macs) {
                retained[i].mac = mac;
            }
        }
        let cols: Vec<Vec<f64>> = phis
            .iter()
            .map(|phi| phi.iter().map(|v| v * h2).collect())
            .collect();
        drop(kphis);
        // Gram of the sampled shapes in the area metric: G = Phi' W Phi.
        let mut gram = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let g: f64 = cols[i]
                    .iter()
                    .zip(&cols[j])
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    / h2;
                gram[i * n + j] = g;
                gram[j * n + i] = g;
            }
        }
        let inverse = invert_spd(&gram, n)?;
        let sandwich = |diag: &dyn Fn(usize) -> f64| {
            let mut out = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    out[i * n + j] = (0..n)
                        .map(|r| inverse[i * n + r] * diag(r) * inverse[r * n + j])
                        .sum();
                }
            }
            out
        };
        let film_damping = film.map(|_| sandwich(&|r| retained[r].film_resistance));
        let damping = sandwich(&|r| retained[r].resistance + retained[r].film_resistance + retained[r].viscous_resistance);
        let mass_matrix = sandwich(&|r| retained[r].air_mass);
        let mut inverse_mass_gram = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v: f64 = cols[i]
                    .iter()
                    .zip(&cols[j])
                    .zip(mass)
                    .map(|((a, b), m)| a * b / m)
                    .sum();
                inverse_mass_gram[i * n + j] = v;
                inverse_mass_gram[j * n + i] = v;
            }
        }
        let rate_bound: f64 = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| damping[i * n + j] * inverse_mass_gram[j * n + i])
                    .sum::<f64>()
            })
            .sum();
        // First-order back-reaction of the dropped in-family coupling at each
        // retained mode's loaded frequency (same parity only: the families
        // of a circular baffled disk do not couple cos to sin).
        let (mut loss_share, mut mass_share) = (0.0_f64, 0.0_f64);
        let mut shares = Vec::with_capacity(n);
        for a in &retained {
            // Beyond the table there is no coupling entry to bound.
            if !a.in_table {
                shares.push((f64::NAN, f64::NAN));
                continue;
            }
            let (mut d_re, mut d_im) = (0.0, 0.0);
            for b in &retained {
                if std::ptr::eq(a, b)
                    || !b.in_table
                    || b.sine != a.sine
                    || shapes[b.shape].m != shapes[a.shape].m
                {
                    continue;
                }
                let w = a.omega;
                // A partner the table cannot evaluate here leaves the share unbounded.
                let (Ok(z), Ok(zs)) = (
                    table.zhat(a.shape, b.shape, a.ka),
                    table.zhat(b.shape, b.shape, a.ka),
                ) else {
                    (d_re, d_im) = (f64::INFINITY, f64::INFINITY);
                    break;
                };
                let (zr, zi) = (scale * z.resistance.unwrap_or(0.0), scale * z.reactance);
                let (bz_re, bz_im) = (zs.resistance.unwrap_or(0.0), zs.reactance);
                // D_s = K_s - w^2 M_s - i w Z_ss(w); delta = i w Z_rs Z_sr / D_s.
                let (dr, di) = (
                    b.stiffness - w * w * b.mass + w * scale * bz_im,
                    -w * scale * bz_re,
                );
                let (pr, pi) = (zr * zr - zi * zi, 2.0 * zr * zi);
                let (nr, ni) = (-w * pi, w * pr);
                let den = dr * dr + di * di;
                d_re += (nr * dr + ni * di) / den;
                d_im += (ni * dr - nr * di) / den;
            }
            let share = (
                if a.resistance > 0.0 {
                    d_re.abs() / a.resistance
                } else {
                    0.0
                },
                d_im.abs() / (a.omega * a.air_mass),
            );
            loss_share = loss_share.max(share.0);
            mass_share = mass_share.max(share.1);
            shares.push(share);
        }
        let mut receipt = Vec::new();
        let mut low_mac = 0;
        for (r, share) in retained.iter().zip(&shares) {
            let shape = &shapes[r.shape];
            // A volume shape's compact share is part of its exterior load.
            let (r_c, m_c) = r.compact.unwrap_or((0.0, 0.0));
            let t60 = |resistance: f64| 6.907755 * 2.0 * (r.mass + r.air_mass + m_c) / resistance;
            let loss = match (r.in_table, r.below_resolution) {
                (true, None) => format!("radiation T60 {:.3} s", t60(r.resistance + r_c)),
                (true, Some(bound)) => format!(
                    "radiation below the table's resolution, none applied: its T60 is an upper bound (true radiation T60 > {:.3e} s)",
                    t60(bound)
                ),
                (false, _) => format!(
                    "radiation NOT TABULATED (j_mn {:.3} lies beyond the table's set): no radiation resistance{}, so its T60 is an upper bound; air mass of the untruncated shape's evanescent field",
                    shape.zero,
                    if r_c > 0.0 { format!(" but the compact monopole's (T60 {:.3e} s)", t60(r_c)) } else { String::new() }
                ),
            };
            let tag = match (shape.m == 0, r.in_table) {
                (false, true) => "",
                (true, true) => " [volume]",
                (false, false) => " [untabulated]",
                (true, false) => " [volume untabulated]",
            };
            let interpolation = if r.in_table {
                format!("ka interpolation <= {:.1e}; dropped in-family coupling <= {:.1e} of its resistance, {:.1e} of its air mass",
                    table.interpolation_error(r.shape, r.ka), share.0, share.1)
            } else {
                "ka interpolation n/a; dropped in-family coupling not tabulated".to_string()
            };
            let extended_note = if extend.is_some() {
                let mut note = String::new();
                if r.compact.is_some() {
                    let _ = write!(note, "; compact net-volume share: R_c V^2 {r_c:.3e} kg/s, piston air mass {m_c:.3e} kg (exterior air mass in all {:.3e} kg)", r.air_mass + m_c);
                }
                if let Some(check) = r.evanescent_check {
                    let _ = write!(note, "; evanescent air mass {:+.1} % against the table's", 100.0 * check);
                }
                if viscous.is_some() {
                    let _ = write!(note, "; viscous T60 {:.3} s", t60(r.viscous_resistance));
                }
                let _ = write!(note, "; total T60 {:.3} s", t60(r.resistance + r_c + r.film_resistance + r.viscous_resistance));
                note
            } else {
                String::new()
            };
            let caveat = if r.mac < LOW_MAC {
                low_mac += 1;
                "; LOW MAC: the FD mode differs from this shape, so its load is approximate"
            } else {
                ""
            };
            let film_note = match (film, r.film_tan_delta) {
                (Some(_), Some(tan_delta)) => format!(
                    "; film: tan delta {tan_delta:.4}, bending share {:.3e}, eta {:.3e}, film T60 {:.4} s",
                    r.bending_share,
                    tan_delta * r.bending_share,
                    t60(r.film_resistance)
                ),
                (Some(_), None) => format!(
                    "; film: NONE, the loaded frequency is outside the film source's band (bending share {:.3e}), so its T60 is an upper bound",
                    r.bending_share
                ),
                (None, _) => String::new(),
            };
            receipt.push(format!(
                "({},{}) {}{tag}: in vacuo {:.1} Hz, loaded {:.1} Hz (ka {:.4}, {} iterations), air mass {:.3e} kg on {:.3e} kg, {loss}; MAC {:.4}{caveat}; {interpolation}{film_note}{extended_note}",
                shape.m, shape.n, if r.sine { "sin" } else { "cos" },
                (r.stiffness / r.mass).sqrt() / (2.0 * std::f64::consts::PI),
                r.omega / (2.0 * std::f64::consts::PI), r.ka, r.iterations, r.air_mass, r.mass, r.mac,
            ));
        }
        // At a glance: does a refusal reach into the band that dominates the decay?
        let refused_note = refused_shapes
            .iter()
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .map_or_else(
                || "no shape refused by the table".to_string(),
                |(hz, label)| {
                    format!(
                        "{} shapes refused by the table, the lowest {label} at {hz:.1} Hz in vacuo",
                        refused_shapes.len()
                    )
                },
            );
        let mut summary = String::new();
        let _ = write!(
            summary,
            "head radiation: {n} shapes from table {:016x} ({ROUTE}, {GEOMETRY}, {MODE_SET}); largest dropped in-family coupling, first order: {:.2e} of a retained resistance, {:.2e} of a retained air mass; narrowband: Re Z_rr and -Im Z_rr/omega at each mode's loaded frequency (fixed point, refused beyond {MAX_LOADED_ITERATIONS} iterations); ka interpolation <= {:.1e} over the whole table (each shape states its own), discretization <= {:.1e} (relative); {} m >= 1 entries dropped as upper bounds (modes the grid cannot resolve, shapes the table refused); {refused_note}; {low_mac} shapes below MAC {LOW_MAC}",
            table.content_hash,
            loss_share,
            mass_share,
            table.interpolation_bound,
            table.discretization_bound,
            dropped.len()
        );
        if extend.is_some() {
            let volume = retained.iter().filter(|r| r.compact.is_some()).count();
            let beyond = retained.iter().filter(|r| !r.in_table).count();
            let checks = retained.iter().filter_map(|r| r.evanescent_check.filter(|_| shapes[r.shape].zero >= 15.0));
            let (lo, hi) = checks.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)));
            let check = if lo <= hi {
                format!("{:+.1} % to {:+.1} %", 100.0 * lo, 100.0 * hi)
            } else {
                "no such shape retained on this grid".to_string()
            };
            let _ = write!(
                summary,
                "; extended to every shape the grid resolves (the head's own losses are modelled): {volume} volume shapes (m = 0, n >= 2) take the table's baffled load (or beyond it the evanescent air mass) less the compact net-volume model's own share, (0,1) staying on the compact model alone; {beyond} shapes beyond the table's set take the air mass of their untruncated evanescent field, rho A / sqrt(k_s^2 - k^2), and no radiation resistance, so their radiation is an upper bound; that air mass against the table's own shapes with j_mn >= 15: {check}"
            );
        }
        receipt.push(summary);
        if film.is_some() {
            let hz = |r: &Retained| r.omega / (2.0 * std::f64::consts::PI);
            let span = |values: &mut dyn Iterator<Item = f64>| {
                values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)))
            };
            let covered: Vec<&Retained> = retained.iter().filter(|r| r.film_tan_delta.is_some()).collect();
            let (low, high) = span(&mut covered.iter().map(|r| hz(r)));
            let (t_low, t_high) = span(&mut covered.iter().filter_map(|r| r.film_tan_delta));
            let tangent = if t_low == t_high {
                format!("tan delta {t_low}")
            } else {
                format!("tan delta {t_low:.4} to {t_high:.4}")
            };
            let outside: Vec<String> = retained
                .iter()
                .filter(|r| r.film_tan_delta.is_none())
                .map(|r| {
                    let shape = &shapes[r.shape];
                    format!("({},{}) {} at {:.1} Hz", shape.m, shape.n, if r.sine { "sin" } else { "cos" }, hz(r))
                })
                .collect();
            let outside = if outside.is_empty() {
                String::new()
            } else {
                format!("; {} retained shapes lie outside the film source's band and carry NO film loss: {}", outside.len(), outside.join(", "))
            };
            receipt.push(if covered.is_empty() {
                format!("head film loss: NONE applied, no retained shape lies inside the film source's band{outside}")
            } else if extend.is_some() {
                format!(
                    "head film loss: {tangent} at each mode's own loaded frequency on {} of the {n} retained shapes (every shape the grid resolves but (0,1)), loaded {low:.1} to {high:.1} Hz (eta = tan delta x bending share){outside}; the (0,1) volume shape and every shape the grid cannot resolve carry NO film loss, so their decays stay upper bounds",
                    covered.len()
                )
            } else {
                format!(
                    "head film loss: {tangent} at each mode's own loaded frequency on {} of the {n} retained m >= 1 shapes, loaded {low:.1} to {high:.1} Hz (eta = tan delta x bending share){outside}; the volume (m = 0) shapes and every dropped or refused shape carry NO film loss, so their decays stay upper bounds",
                    covered.len()
                )
            });
        }
        if viscous.is_some() {
            let (low, high) = retained.iter().map(|r| r.viscous_resistance / (r.mass * r.omega))
                .fold((f64::INFINITY, 0.0_f64), |(lo, hi), eta| (lo.min(eta), hi.max(eta)));
            let uncovered = if extend.is_some() {
                "(0,1) and the shapes the grid cannot resolve carry none"
            } else {
                "m = 0 and dropped shapes carry none"
            };
            receipt.push(format!(
                "head viscous layer: exact linear Stokes face loss on both faces, eta {low:.3e} to {high:.3e} at each retained mode's loaded frequency; {uncovered}, so their decays stay upper bounds; every retained mode is below coincidence"
            ));
        }
        receipt.push(format!(
            "head radiation retention: a shape needs {NODES_PER_WAVELENGTH} grid spacings per radial and per azimuthal wavelength, j_mn <= 2 pi a / ({NODES_PER_WAVELENGTH} h) = {limit:.3} and m <= {limit:.3}"
        ));
        receipt.extend(dropped.iter().map(|line| format!("head radiation {line}")));
        receipt.extend(volume_lines);
        receipt.push(if extend.is_some() {
            "head radiation scope: each head in its own infinite rigid baffle; no shell diffraction; no exterior coupling between the batter and resonant heads' m >= 1 motion; no cavity loading of m >= 1 modes; shapes beyond the table's set carry their evanescent air mass and no radiation resistance; shapes the grid cannot resolve carry no exterior load".into()
        } else {
            "head radiation scope: each head in its own infinite rigid baffle; no shell diffraction; no exterior coupling between the batter and resonant heads' m >= 1 motion; no cavity loading of m >= 1 modes; shapes outside the table's set carry no exterior load".into()
        });
        let nodes = xy.len();
        Ok(ModalAirLoad {
            cols,
            gram,
            air: retained.iter().map(|r| r.air_mass).collect(),
            damping,
            mass: mass_matrix,
            rate_bound,
            film_damping,
            film_step_j: 0.0,
            film_dissipated_j: 0.0,
            prepared: None,
            rhs: vec![0.0; nodes],
            corr: vec![0.0; nodes],
            p_wp: vec![0.0; n],
            y: vec![0.0; n + 1],
            b: vec![0.0; n + 1],
            w2: vec![0.0; n + 1],
            za: vec![0.0; n],
            zb: vec![0.0; n],
            receipt,
            dropped_loss_share: loss_share,
            dropped_mass_share: mass_share,
        })
    }

    pub fn shapes(&self) -> usize {
        self.cols.len()
    }

    /// `out = M v` for a modal-space matrix stored row-major.
    pub(crate) fn apply(matrix: &[f64], v: &[f64], out: &mut [f64]) {
        let n = v.len();
        for i in 0..n {
            out[i] = (0..n).map(|j| matrix[i * n + j] * v[j]).sum();
        }
    }

    /// The modal part of a step's right-hand side,
    /// `rhs = U (A (2 P_w - P_wp) / k^2 - B (P_w - P_wp) / k)` with `P = U' w`:
    /// the air mass's share of the leapfrog inertia and the backward-velocity
    /// radiation force. Leaves `U' wp` in `p_wp` for the energy account.
    pub(crate) fn load_rhs(&mut self, w: &[f64], wp: &[f64], k: f64) {
        let n = self.cols.len();
        for r in 0..n {
            let (pw, pwp) = (crate::dot(&self.cols[r], w), crate::dot(&self.cols[r], wp));
            self.p_wp[r] = pwp;
            self.y[r] = (2.0 * pw - pwp) / (k * k);
            self.b[r] = (pw - pwp) / k;
        }
        ModalAirLoad::apply(&self.mass, &self.y[..n], &mut self.za);
        ModalAirLoad::apply(&self.damping, &self.b[..n], &mut self.zb);
        self.rhs.fill(0.0);
        for r in 0..n {
            let z = self.za[r] - self.zb[r];
            for (out, u) in self.rhs.iter_mut().zip(&self.cols[r]) {
                *out += z * u;
            }
        }
    }

    /// Energy the modal damping removes over one step at centred velocity
    /// `(wn - wp) / 2k`: `(U'(wn - wp))' B (U'(wn - wp)) / 4k`. With a film
    /// loss this is radiation plus film; `film_step_j` is the film part.
    pub(crate) fn radiated(&mut self, wn: &[f64], k: f64) -> f64 {
        let n = self.cols.len();
        for r in 0..n {
            self.y[r] = crate::dot(&self.cols[r], wn) - self.p_wp[r];
        }
        if let Some(film) = &self.film_damping {
            ModalAirLoad::apply(film, &self.y[..n], &mut self.zb);
            self.film_step_j = crate::dot(&self.y[..n], &self.zb) / (4.0 * k);
            self.film_dissipated_j += self.film_step_j;
        }
        ModalAirLoad::apply(&self.damping, &self.y[..n], &mut self.zb);
        crate::dot(&self.y[..n], &self.zb) / (4.0 * k)
    }

    /// The film share of the last `radiated` step [J]; 0 without a film loss.
    pub fn film_step_j(&self) -> f64 {
        self.film_step_j
    }

    /// Energy the film loss has removed so far [J].
    pub fn film_dissipated_j(&self) -> f64 {
        self.film_dissipated_j
    }

    /// Prepare the Woodbury factor for step `k`, nodal compliances `ia`
    /// (`D^-1`) and an optional volume column `u = h^2` with weight `mu`.
    pub(crate) fn prepare(
        &mut self,
        k: f64,
        ia: &[f64],
        volume: Option<(f64, f64)>,
    ) -> Result<(), String> {
        let n = self.cols.len();
        let size = n + usize::from(volume.is_some());
        let mut q = vec![0.0; size * size];
        let off = usize::from(volume.is_some());
        for i in 0..n {
            for j in 0..n {
                let v: f64 = (0..n)
                    .map(|r| self.gram[i * n + r] * self.gram[r * n + j] / self.air[r])
                    .sum();
                q[(i + off) * size + j + off] = k * k * v;
            }
        }
        let column = |c: usize| -> Vec<f64> {
            if let (Some((h2, _)), 0) = (volume, c) {
                vec![h2; ia.len()]
            } else {
                self.cols[c - off].clone()
            }
        };
        if let Some((_, mu)) = volume {
            q[0] += 1.0 / mu;
        }
        let columns: Vec<Vec<f64>> = (0..size).map(column).collect();
        for i in 0..size {
            for j in 0..=i {
                let v: f64 = columns[i]
                    .iter()
                    .zip(&columns[j])
                    .zip(ia)
                    .map(|((a, b), d)| a * b * d)
                    .sum();
                q[i * size + j] += v;
                if i != j {
                    q[j * size + i] += v;
                }
            }
        }
        cholesky(&mut q, size)?;
        self.prepared = Some(Prepared {
            k,
            chol: q,
            with_volume: volume.is_some(),
        });
        Ok(())
    }
}
