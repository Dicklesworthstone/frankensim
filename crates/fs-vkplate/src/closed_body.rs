//! What an overhead microphone hears from a drum as ONE rigid closed body
//! (pexfp.378): batter, resonant head and a rigid shell on a rigid floor.
//!
//! The table is baked by fs-bem's `examples/drum_observer.rs` for one kit
//! piece, its gas and its internal clock: the closed body's pressure at the
//! microphone per unit modal rate of every clamped-membrane shape with
//! `j_mn <= 9` on either head, from fs-bem's routed triangle CBIE, below the
//! body's CBIE band top f_x; the floor is an exact image body. The rows are
//! the producer's exp(-i omega t) samples on the bake's finest FFT grid. Here
//! they become the bake's own FIR kernels (conjugate, integer-delay and
//! latency phase, C4 taper, inverse real DFT); each kernel's energy, energy
//! centroid and one in-band transfer must match fs-bem's within 1e-9, or the
//! table is refused.
//!
//! The microphone output is then
//!   y = R(t - D) - (g * U)(t) + sum_s (h_s * qdot_s)(t),
//! where R is today's baffled Rayleigh sum, qdot_s the least-squares modal
//! rates of the retained shapes, h_s the closed-body kernels, U the baffled
//! Rayleigh sum of the retained shapes driven by k qdot_s (a centred half
//! difference) instead of their centred second difference, and g the bake's
//! kernel of 2 i tan(omega / 2), the ratio of the two: so g * U is the taper
//! on R_ret, the baffled sum of the retained shapes' projection, delayed by D.
//! (A unit-transfer source cannot carry the taper: the bake zeroes the DC bin
//! and refuses the floor that leaves.) Below f_x the retained
//! shapes are heard as the closed body; above it through the exact complement
//! of the same taper; head motion outside the retained shapes stays on the
//! baffled path at every frequency, and its share is reported. Where the two
//! hearings agree the bands sum to a pure delay of D samples, which the
//! render trims. The observer is outside the loop: the head motion is
//! unchanged. No table: the microphone is untouched.
//!
//! Only arithmetic, `sqrt` and `fs_math::det`'s sin, cos and tan are used, so
//! every build reconstructs the same kernels.

use crate::Body;
use crate::radiation::{MODE_SET, ModeShape};
use fs_math::det;
use std::f64::consts::TAU;

pub const TABLE_FORMAT: &str = "frankensim.vkplate-closed-body-observer.v1";
pub const ROUTE: &str = "closed-body-triangle-cbie";
/// The bake's both-ways transfer limit (.376's reciprocal-part number).
pub const RECIPROCITY_LIMIT: f64 = 0.02;
/// A reconstructed kernel's statistics must match fs-bem's within this.
pub const KERNEL_MATCH: f64 = 1e-9;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeadSide {
    Batter,
    Resonant,
}

#[derive(Clone, Debug)]
pub struct TableSource {
    pub head: HeadSide,
    pub shape: ModeShape,
    pub sine: bool,
}

/// One kernel: integer delay and taps (fs-vfit `FirKernel` order).
#[derive(Clone, Debug)]
pub(crate) struct Kernel {
    pub(crate) integer: usize,
    pub(crate) taps: Vec<f64>,
}

/// A parsed, reconstructed and checked table.
#[derive(Clone, Debug)]
pub struct ClosedBodyTable {
    pub piece: String,
    pub radius: f64,
    pub depth: f64,
    pub mic: [f64; 3],
    pub floor: ([f64; 3], [f64; 3]),
    pub gas: (f64, f64),
    pub rate: f64,
    pub pass_hz: f64,
    pub stop_hz: f64,
    pub latency: usize,
    pub n_fft: usize,
    pub declaration: String,
    pub reciprocity: f64,
    pub crossover_db: (f64, f64),
    pub notes: Vec<String>,
    pub sources: Vec<TableSource>,
    /// FNV-1a of the table text: the render receipt's identity.
    pub digest: u64,
    pub(crate) kernels: Vec<Kernel>,
    pub(crate) complement: Kernel,
}

fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325_u64, |h, b| (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3))
}

/// fs-bem's C4 observation taper: 1 through `pass`, 0 from `stop`.
pub(crate) fn taper(f: f64, pass: f64, stop: f64) -> f64 {
    if f <= pass {
        return 1.0;
    }
    if f >= stop {
        return 0.0;
    }
    let low = |t: f64| t * t * t * t * t * (126.0 + t * (-420.0 + t * (540.0 + t * (-315.0 + 70.0 * t))));
    let x = (f - pass) / (stop - pass);
    if x <= 0.5 { 1.0 - low(x) } else { low(1.0 - x) }
}

/// cos and sin of 2 pi j / n, j < n.
pub(crate) fn twiddles(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|j| {
            let (s, c) = det::sin_cos(TAU * j as f64 / n as f64);
            (c, s)
        })
        .collect()
}

/// The bake's complement source at finest bin `b`: 2 i tan(omega / 2 rate) in
/// exp(+i omega t), its producer row the conjugate, evaluated exactly as fs-bem's
/// `examples/drum_observer.rs` does from the bin's angular frequency.
pub(crate) fn complement_row(rate: f64, n: usize, b: usize) -> (f64, f64) {
    let f = rate * b as f64 / n as f64;
    let omega = TAU * f;
    (0.0, -2.0 * det::tan(0.5 * omega / rate))
}

/// The bake's kernel from its producer rows `rows[b - 1]` (bins 1..=B):
/// X_b = conj(row) e^{i 2 pi b (integer - latency) / n} W(f_b), taps the first
/// n/2 samples of the inverse real DFT (1/n normalized).
pub(crate) fn reconstruct(
    rows: &[(f64, f64)],
    integer: usize,
    latency: usize,
    n: usize,
    rate: f64,
    band: (f64, f64),
    tw: &[(f64, f64)],
) -> Vec<f64> {
    let shift = integer as i64 - latency as i64;
    let mut spectrum = Vec::with_capacity(rows.len());
    for (i, &(re, im)) in rows.iter().enumerate() {
        let b = i + 1;
        let f = rate * b as f64 / n as f64;
        let w = taper(f, band.0, band.1);
        let (c, s) = tw[((b as i64 * shift).rem_euclid(n as i64)) as usize];
        // conj(row) * (c + i s)
        spectrum.push(((re * c + im * s) * w, (re * s - im * c) * w));
    }
    let scale = 2.0 / n as f64;
    (0..n / 2)
        .map(|t| {
            let mut sum = 0.0;
            for (i, &(xr, xi)) in spectrum.iter().enumerate() {
                let (c, s) = tw[((i + 1) * t) % n];
                sum += xr * c - xi * s;
            }
            sum * scale
        })
        .collect()
}

/// Energy, energy centroid and transfer at `omega` [rad/sample] (exp(-i omega t)).
pub(crate) fn kernel_stats(k: &Kernel, omega: f64) -> (f64, f64, (f64, f64)) {
    let energy: f64 = k.taps.iter().map(|t| t * t).sum();
    let moment: f64 = k.taps.iter().enumerate().map(|(i, t)| i as f64 * t * t).sum();
    let (mut re, mut im) = (0.0, 0.0);
    for (i, &t) in k.taps.iter().enumerate() {
        let (s, c) = det::sin_cos(-omega * (k.integer + i) as f64);
        re += t * c;
        im += t * s;
    }
    (energy, if energy > 0.0 { moment / energy } else { 0.0 }, (re, im))
}

impl ClosedBodyTable {
    pub fn parse(text: &str) -> Result<ClosedBodyTable, String> {
        let refuse = |what: &str| format!("closed-body table: {what}");
        let mut lines = text.lines();
        if lines.next() != Some(TABLE_FORMAT) {
            return Err(refuse("not a frankensim.vkplate-closed-body-observer.v1 table"));
        }
        let number = |s: &str| s.parse::<f64>().map_err(|_| refuse("malformed number"));
        let index = |s: &str| s.parse::<usize>().map_err(|_| refuse("malformed index"));
        let mut once: std::collections::BTreeMap<&str, Vec<&str>> = std::collections::BTreeMap::new();
        let mut notes = Vec::new();
        let mut sources: Vec<(usize, TableSource, f64)> = Vec::new();
        let mut probes: Vec<(usize, f64, f64, f64)> = Vec::new();
        let mut stats: Vec<(usize, usize, usize, f64, f64, f64, f64, String)> = Vec::new();
        let mut rows: Vec<(usize, usize, f64, f64)> = Vec::new();
        for line in lines {
            let cols: Vec<&str> = line.split('\t').collect();
            match cols[0] {
                "note" | "mesh" | "bake-args" => notes.push(line.to_string()),
                "source" if cols.len() == 7 => {
                    let head = match cols[2] {
                        "batter" => HeadSide::Batter,
                        "resonant" => HeadSide::Resonant,
                        _ => return Err(refuse("unknown head")),
                    };
                    let sine = match cols[5] {
                        "cos" => false,
                        "sin" => true,
                        _ => return Err(refuse("unknown parity")),
                    };
                    let shape = ModeShape::new(index(cols[3])?, index(cols[4])?);
                    sources.push((index(cols[1])?, TableSource { head, shape, sine }, number(cols[6])?));
                }
                "probe" if cols.len() == 5 => {
                    probes.push((index(cols[1])?, number(cols[2])?, number(cols[3])?, number(cols[4])?))
                }
                "kernel" if cols.len() == 9 => stats.push((
                    index(cols[1])?,
                    index(cols[2])?,
                    index(cols[3])?,
                    number(cols[4])?,
                    number(cols[5])?,
                    number(cols[6])?,
                    number(cols[7])?,
                    cols[8].to_string(),
                )),
                "row" if cols.len() == 5 => {
                    rows.push((index(cols[1])?, index(cols[2])?, number(cols[3])?, number(cols[4])?))
                }
                key @ ("route" | "piece" | "drum" | "heads" | "mic" | "floor" | "gas" | "rate" | "band"
                | "latency" | "fft" | "declaration" | "reciprocity" | "crossover-step-db" | "mode-set"
                | "end") => {
                    if once.insert(key, cols[1..].to_vec()).is_some() {
                        return Err(refuse("a header line repeats"));
                    }
                }
                _ => return Err(refuse("unknown or malformed line")),
            }
        }
        let field = |key: &str, count: usize| -> Result<Vec<String>, String> {
            match once.get(key) {
                Some(v) if v.len() >= count => Ok(v.iter().map(|c| c.to_string()).collect()),
                _ => Err(refuse(&format!("missing or short `{key}` line"))),
            }
        };
        if field("route", 1)?[0] != ROUTE {
            return Err(refuse("not the closed-body triangle-CBIE route"));
        }
        if field("mode-set", 1)?[0] != MODE_SET {
            return Err(refuse("not this crate's clamped-membrane shape set"));
        }
        let drum = field("drum", 2)?;
        let (radius, depth) = (number(&drum[0])?, number(&drum[1])?);
        let heads = field("heads", 10)?;
        let expected = ["batter", "0.0", "0.0", "0.0", "1.0", "resonant", "0.0", "0.0"];
        if heads[..8] != expected || number(&heads[8])? != -depth || heads[9] != "-1.0" {
            return Err(refuse("the heads are not the kit builder's (batter at the origin facing +z, resonant at -depth facing -z)"));
        }
        let mic = field("mic", 3)?;
        let mic = [number(&mic[0])?, number(&mic[1])?, number(&mic[2])?];
        let floor = field("floor", 6)?;
        let f: Vec<f64> = floor[..6].iter().map(|s| number(s.as_str())).collect::<Result<_, _>>()?;
        let floor = ([f[0], f[1], f[2]], [f[3], f[4], f[5]]);
        let gas = field("gas", 2)?;
        let gas = (number(&gas[0])?, number(&gas[1])?);
        let rate = number(&field("rate", 1)?[0])?;
        let band = field("band", 2)?;
        let (pass_hz, stop_hz) = (number(&band[0])?, number(&band[1])?);
        let latency = index(&field("latency", 1)?[0])?;
        let n_fft = index(&field("fft", 1)?[0])?;
        let reciprocity = number(&field("reciprocity", 1)?[0])?;
        let step = field("crossover-step-db", 2)?;
        let crossover_db = (number(&step[0])?, number(&step[1])?);
        let last_bin = index(&field("end", 1)?[0])?;
        if !(radius > 0.0 && depth > 0.0 && rate > 0.0 && pass_hz > 0.0 && pass_hz < stop_hz && stop_hz < 0.5 * rate)
            || latency == 0
            || !n_fft.is_power_of_two()
            || n_fft < 64
            || latency > n_fft / 4
        {
            return Err(refuse("invalid band, rate, latency or FFT length"));
        }
        if !(reciprocity <= RECIPROCITY_LIMIT) {
            return Err(refuse("its both-ways transfer check exceeds the reciprocity limit"));
        }
        // Bins 1..=last_bin are exactly those below the stop.
        let bin_hz = |b: usize| rate * b as f64 / n_fft as f64;
        if last_bin == 0 || bin_hz(last_bin) >= stop_hz || bin_hz(last_bin + 1) < stop_hz {
            return Err(refuse("its rows do not end at the taper stop"));
        }
        sources.sort_by_key(|s| s.0);
        if sources.is_empty() || sources.iter().enumerate().any(|(i, s)| s.0 != i) {
            return Err(refuse("source indices are not 0..S"));
        }
        for (i, s, zero) in &sources {
            if (s.shape.zero - zero).abs() > 1e-12 * zero {
                return Err(refuse("a Bessel zero differs from this crate's"));
            }
            let checked: Vec<_> = probes.iter().filter(|p| p.0 == *i).collect();
            if checked.is_empty() {
                return Err(refuse("a source has no probe values"));
            }
            for &&(_, x, y, value) in &checked {
                if (s.shape.value(x, y, s.sine) - value).abs() > 1e-12 * value.abs().max(1.0) {
                    return Err(refuse("a probe value differs from this crate's shape"));
                }
            }
        }
        if probes.iter().any(|p| p.0 >= sources.len()) {
            return Err(refuse("a probe names an undeclared source"));
        }
        let count = sources.len();
        let mut table_rows = vec![vec![None; last_bin]; count];
        for &(i, b, re, im) in &rows {
            if i >= count || b == 0 || b > last_bin || table_rows[i][b - 1].is_some() || !re.is_finite() || !im.is_finite() {
                return Err(refuse("a row is out of range, repeated or nonfinite"));
            }
            table_rows[i][b - 1] = Some((re, im));
        }
        let tw = twiddles(n_fft);
        let probe_omega = TAU * 0.5 * pass_hz / rate;
        stats.sort_by_key(|s| s.0);
        if stats.len() != count + 1 || stats.iter().enumerate().any(|(i, s)| s.0 != i) {
            return Err(refuse("kernel lines are not 0..=S"));
        }
        let mut kernels = Vec::with_capacity(count);
        let mut complement = None;
        for (i, &(_, integer, length, energy, centre, re, im, ref kind)) in stats.iter().enumerate() {
            if length != n_fft / 2 || (kind == "complement") != (i == count) {
                return Err(refuse("a kernel's length or kind is not the bake's"));
            }
            let source_rows: Vec<(f64, f64)> = if i == count {
                if integer != 0 {
                    return Err(refuse("the complement kernel has a delay"));
                }
                (1..=last_bin).map(|b| complement_row(rate, n_fft, b)).collect()
            } else {
                table_rows[i].iter().map(|r| r.ok_or_else(|| refuse("a source misses a bin"))).collect::<Result<_, _>>()?
            };
            let kernel = Kernel {
                integer,
                taps: reconstruct(&source_rows, integer, latency, n_fft, rate, (pass_hz, stop_hz), &tw),
            };
            let (e, c, (hr, hi)) = kernel_stats(&kernel, probe_omega);
            let abs_sum: f64 = kernel.taps.iter().map(|t| t.abs()).sum();
            if (e - energy).abs() > KERNEL_MATCH * energy
                || (c - centre).abs() > KERNEL_MATCH * centre.max(1.0)
                || (hr - re).hypot(hi - im) > KERNEL_MATCH * abs_sum
            {
                return Err(refuse(&format!(
                    "kernel {i} reconstructs differently from fs-bem's (energy {e:e} vs {energy:e}, centroid {c} vs {centre}, transfer ({hr:e}, {hi:e}) vs ({re:e}, {im:e}))"
                )));
            }
            if i == count {
                complement = Some(kernel);
            } else {
                kernels.push(kernel);
            }
        }
        Ok(ClosedBodyTable {
            piece: field("piece", 1)?[0].clone(),
            radius,
            depth,
            mic,
            floor,
            gas,
            rate,
            pass_hz,
            stop_hz,
            latency,
            n_fft,
            declaration: field("declaration", 1)?[0].clone(),
            reciprocity,
            crossover_db,
            notes,
            sources: sources.into_iter().map(|s| s.1).collect(),
            digest: fnv(text),
            kernels,
            complement: complement.ok_or_else(|| refuse("no complement kernel"))?,
        })
    }

    /// The route receipt for the render's notes.
    pub fn receipt(&self) -> String {
        format!(
            "overhead hears the drum as one closed body (pexfp.378): table {} fnv {:016x}, declaration {}; {} shapes (j_mn <= 9, both heads) from fs-bem's triangle CBIE below f_x {:.0} Hz (pass {:.0} Hz), the floor an exact image body, reciprocity {:.2e}; above f_x the baffled Rayleigh path through the taper's exact complement, crossover step (0,1) {:+.1} dB, (0,2) {:+.1} dB; head motion outside the shapes stays on the baffled path (share reported after the render); latency {} samples trimmed",
            self.piece,
            self.digest,
            self.declaration,
            self.sources.len(),
            self.stop_hz,
            self.pass_hz,
            self.reciprocity,
            self.crossover_db.0,
            self.crossover_db.1,
            self.latency
        )
    }
}

/// One head's retained shapes on its grid.
#[derive(Clone, Debug)]
struct HeadState {
    body: usize,
    /// Global source indices of this head's shapes.
    sources: Vec<usize>,
    /// Least-squares projection rows `G^-1 h^2 phi` (one per source).
    proj: Vec<Vec<f64>>,
    /// The shapes at the nodes.
    phi: Vec<Vec<f64>>,
}

/// A history with a contiguous window: each sample is written twice.
#[derive(Clone, Debug)]
struct History {
    buf: Vec<f64>,
    len: usize,
}

impl History {
    fn new(span: usize) -> History {
        let len = span.next_power_of_two().max(2);
        History { buf: vec![0.0; 2 * len], len }
    }
    fn write(&mut self, t: usize, x: f64) {
        let p = t & (self.len - 1);
        self.buf[p] = x;
        self.buf[p + self.len] = x;
    }
    /// sum_k taps[k] x[t - delay - k], `reversed` = taps reversed.
    fn dot(&self, t: usize, delay: usize, reversed: &[f64]) -> f64 {
        let p = (t & (self.len - 1)) + self.len;
        let start = p - delay - (reversed.len() - 1);
        reversed.iter().zip(&self.buf[start..start + reversed.len()]).map(|(a, b)| a * b).sum()
    }
    fn at(&self, t: usize, lag: usize) -> f64 {
        self.buf[(t & (self.len - 1)) + self.len - lag]
    }
}

/// The observer on one microphone.
#[derive(Clone, Debug)]
pub struct ClosedBodyObserver {
    digest: u64,
    latency: usize,
    heads: Vec<HeadState>,
    /// Per source: integer delay and reversed taps; the complement's reversed taps.
    kernels: Vec<(usize, Vec<f64>)>,
    complement: Vec<f64>,
    /// Per source: first delay and the aggregated baffled response of the
    /// shape (built from the microphone's own taps by `prepare`).
    rho: Vec<(usize, Vec<f64>)>,
    /// The baffled sums R_ret (second differences) and U (half differences).
    ring_ret: Vec<f64>,
    ring_u: Vec<f64>,
    qdot: Vec<History>,
    u: History,
    full: History,
    step: usize,
    alpha: Vec<f64>,
    half: Vec<f64>,
    rate_now: Vec<f64>,
    /// Sums of R^2, R_ret^2 and (R - R_ret)^2 over the steps.
    energies: [f64; 3],
}

fn cholesky_inverse(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut l = a.to_vec();
    for j in 0..n {
        let mut d = l[j * n + j];
        for k in 0..j {
            d -= l[j * n + k] * l[j * n + k];
        }
        if !(d > 0.0 && d.is_finite()) {
            return Err("closed-body observer: the head grid does not resolve the retained shapes (singular Gram matrix)".into());
        }
        let d = d.sqrt();
        l[j * n + j] = d;
        for i in j + 1..n {
            let mut v = l[i * n + j];
            for k in 0..j {
                v -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = v / d;
        }
    }
    let mut inverse = vec![0.0; n * n];
    for c in 0..n {
        let mut x = vec![0.0; n];
        x[c] = 1.0;
        for i in 0..n {
            let mut v = x[i];
            for k in 0..i {
                v -= l[i * n + k] * x[k];
            }
            x[i] = v / l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut v = x[i];
            for k in i + 1..n {
                v -= l[k * n + i] * x[k];
            }
            x[i] = v / l[i * n + i];
        }
        for i in 0..n {
            inverse[i * n + c] = x[i];
        }
    }
    Ok(inverse)
}

impl ClosedBodyObserver {
    /// Bind a checked table to the head bodies `[batter, resonant]` (their
    /// indices and bodies). The caller has checked clock, gas, microphone,
    /// floor and head placement against the table.
    pub(crate) fn new(table: &ClosedBodyTable, heads: [(usize, &Body); 2]) -> Result<ClosedBodyObserver, String> {
        let mut states = Vec::new();
        for (side, (index, body)) in [HeadSide::Batter, HeadSide::Resonant].into_iter().zip(heads) {
            let sources: Vec<usize> = (0..table.sources.len()).filter(|&s| table.sources[s].head == side).collect();
            let r_max = body.xy.iter().map(|p| p[0].hypot(p[1])).fold(0.0_f64, f64::max);
            if !(r_max < table.radius && r_max > table.radius - 1.5 * body.h) {
                return Err(format!(
                    "closed-body observer: {}'s grid reaches {r_max:.4} m, not the table's {} m head",
                    body.name, table.radius
                ));
            }
            let phi: Vec<Vec<f64>> = sources
                .iter()
                .map(|&s| {
                    let t = &table.sources[s];
                    body.xy.iter().map(|p| t.shape.value(p[0] / table.radius, p[1] / table.radius, t.sine)).collect()
                })
                .collect();
            let n = sources.len();
            let h2 = body.h * body.h;
            let mut gram = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    gram[i * n + j] = h2 * phi[i].iter().zip(&phi[j]).map(|(a, b)| a * b).sum::<f64>();
                }
            }
            let inverse = cholesky_inverse(&gram, n)?;
            let proj: Vec<Vec<f64>> = (0..n)
                .map(|i| {
                    (0..body.n)
                        .map(|node| h2 * (0..n).map(|j| inverse[i * n + j] * phi[j][node]).sum::<f64>())
                        .collect()
                })
                .collect();
            states.push(HeadState { body: index, sources, proj, phi });
        }
        let span = table.kernels.iter().map(|k| k.integer + k.taps.len()).max().unwrap_or(1);
        let reversed = |taps: &[f64]| taps.iter().rev().copied().collect::<Vec<f64>>();
        Ok(ClosedBodyObserver {
            digest: table.digest,
            latency: table.latency,
            heads: states,
            kernels: table.kernels.iter().map(|k| (k.integer, reversed(&k.taps))).collect(),
            complement: reversed(&table.complement.taps),
            rho: vec![],
            ring_ret: vec![],
            ring_u: vec![],
            qdot: vec![History::new(span + 1); table.sources.len()],
            u: History::new(table.complement.taps.len() + 1),
            full: History::new(table.latency + 1),
            step: 0,
            alpha: vec![0.0; table.sources.len()],
            half: vec![0.0; table.sources.len()],
            rate_now: vec![0.0; table.sources.len()],
            energies: [0.0; 3],
        })
    }

    /// The output delay [internal samples] the render trims.
    pub fn latency(&self) -> usize {
        self.latency
    }

    /// Aggregate the microphone's baffled taps per retained shape.
    pub(crate) fn prepare(&mut self, taps: &[(usize, usize, usize, f64, f64)], ring_len: usize) {
        let total: usize = self.heads.iter().map(|h| h.sources.len()).sum();
        let mut rho = vec![(0usize, Vec::new()); total];
        for head in &self.heads {
            let mine: Vec<&(usize, usize, usize, f64, f64)> = taps.iter().filter(|t| t.0 == head.body).collect();
            let (d0, d1) = mine.iter().fold((usize::MAX, 0), |(a, b), t| (a.min(t.2), b.max(t.2 + 1)));
            for (j, &s) in head.sources.iter().enumerate() {
                if mine.is_empty() {
                    rho[s] = (0, vec![]);
                    continue;
                }
                let mut r = vec![0.0; d1 - d0 + 1];
                for &&(_, node, di, lo, hi) in &mine {
                    let phi = head.phi[j][node];
                    r[di - d0] += lo * phi;
                    r[di - d0 + 1] += hi * phi;
                }
                rho[s] = (d0, r);
            }
        }
        self.rho = rho;
        self.ring_ret = vec![0.0; ring_len];
        self.ring_u = vec![0.0; ring_len];
    }

    /// One step: `full` is today's baffled sample at this step.
    pub(crate) fn step(&mut self, full: f64, bodies: &[Body], k: f64, base: usize) -> f64 {
        for head in &self.heads {
            let b = &bodies[head.body];
            for (j, &s) in head.sources.iter().enumerate() {
                let row = &head.proj[j];
                let (mut a, mut v) = (0.0, 0.0);
                for i in 0..row.len() {
                    a += row[i] * (b.wn[i] - 2.0 * b.w[i] + b.wp[i]);
                    v += row[i] * (b.wn[i] - b.wp[i]);
                }
                self.alpha[s] = a;
                self.half[s] = 0.5 * v;
                self.rate_now[s] = v / (2.0 * k);
            }
        }
        let mask = self.ring_ret.len() - 1;
        for (s, (d0, r)) in self.rho.iter().enumerate() {
            let (a, u) = (self.alpha[s], self.half[s]);
            for (j, &w) in r.iter().enumerate() {
                let slot = (base + d0 + j) & mask;
                self.ring_ret[slot] += w * a;
                self.ring_u[slot] += w * u;
            }
        }
        let slot = base & mask;
        let (retained, u) = (self.ring_ret[slot], self.ring_u[slot]);
        self.ring_ret[slot] = 0.0;
        self.ring_u[slot] = 0.0;
        let t = self.step;
        self.full.write(t, full);
        self.u.write(t, u);
        for (h, &v) in self.qdot.iter_mut().zip(&self.rate_now) {
            h.write(t, v);
        }
        let mut y = if t >= self.latency { self.full.at(t, self.latency) } else { 0.0 };
        y -= self.u.dot(t, 0, &self.complement);
        for ((delay, taps), h) in self.kernels.iter().zip(&self.qdot) {
            y += h.dot(t, *delay, taps);
        }
        self.energies[0] += full * full;
        self.energies[1] += retained * retained;
        self.energies[2] += (full - retained) * (full - retained);
        self.step += 1;
        y
    }

    /// After a render: the baffled signal's split between the retained shapes
    /// and the head motion outside them (heard as today at every frequency).
    pub fn report(&self) -> String {
        let [full, retained, residual] = self.energies;
        let db = |x: f64| if full > 0.0 && x > 0.0 { 10.0 * (x / full).log10() } else { f64::NEG_INFINITY };
        format!(
            "closed-body observer (table fnv {:016x}): of the baffled overhead signal, the retained shapes' projection carries {:+.1} dB and the motion outside them {:+.1} dB (energy re the whole); the latter stays baffled",
            self.digest,
            db(retained),
            db(residual)
        )
    }
}

#[cfg(test)]
#[path = "closed_body_tests.rs"]
mod tests;
