//! T2 (replay and complement) and the planted negatives N2 (a tampered
//! table) and N3 (any header mismatch refuses) of the closed-body observer.
use super::*;
use crate::kit::{self, DrumSpec};
use crate::{Instrument, floor_z};
use std::f64::consts::PI;

/// A synthetic table in the bake's format, its kernel statistics computed by
/// this module's reconstruction (so a correct parse is the only gate).
#[derive(Clone)]
struct Synthetic {
    radius: f64,
    depth: f64,
    mic: [f64; 3],
    floor: ([f64; 3], [f64; 3]),
    gas: (f64, f64),
    rate: f64,
    pass: f64,
    stop: f64,
    latency: usize,
    n: usize,
    sources: Vec<(HeadSide, usize, usize, bool)>,
    integer: Vec<usize>,
}

impl Synthetic {
    fn for_drum(ins: &Instrument, spec: &DrumSpec, sources: Vec<(HeadSide, usize, usize, bool)>) -> Synthetic {
        let floor = spec.floor.unwrap();
        Synthetic {
            radius: spec.radius,
            depth: spec.depth,
            mic: spec.mic,
            floor: (floor.point, floor.normal),
            gas: (ins.gas().density, ins.gas().sound_speed),
            rate: ins.rate,
            pass: 3000.0,
            stop: 6000.0,
            latency: 128,
            n: 2048,
            integer: vec![0; sources.len()],
            sources,
        }
    }

    fn last_bin(&self) -> usize {
        (1..).take_while(|&b| self.rate * b as f64 / (self.n as f64) < self.stop).last().unwrap()
    }

    fn write(&self, rows: &dyn Fn(usize, usize) -> (f64, f64)) -> String {
        let mut t = String::new();
        let mut line = |s: String| {
            t.push_str(&s);
            t.push('\n');
        };
        line(TABLE_FORMAT.into());
        line(format!("route\t{ROUTE}\tsynthetic test table"));
        line(format!("bake-args\ttest\t{:?}", self.rate));
        line("piece\ttest".into());
        line(format!("drum\t{:?}\t{:?}", self.radius, self.depth));
        line(format!("heads\tbatter\t0.0\t0.0\t0.0\t1.0\tresonant\t0.0\t0.0\t{:?}\t-1.0", -self.depth));
        line(format!("mic\t{:?}\t{:?}\t{:?}", self.mic[0], self.mic[1], self.mic[2]));
        let (q, nf) = self.floor;
        line(format!("floor\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}", q[0], q[1], q[2], nf[0], nf[1], nf[2]));
        line(format!("gas\t{:?}\t{:?}", self.gas.0, self.gas.1));
        line(format!("rate\t{:?}", self.rate));
        line(format!("band\t{:?}\t{:?}", self.pass, self.stop));
        line(format!("latency\t{}", self.latency));
        line(format!("fft\t{}", self.n));
        line("declaration\tsynthetic".into());
        line("reciprocity\t0e0".into());
        line("crossover-step-db\t0.0\t0.0".into());
        line(format!("mode-set\t{MODE_SET}\tsynthetic"));
        for (i, &(head, m, n, sine)) in self.sources.iter().enumerate() {
            let shape = ModeShape::new(m, n);
            line(format!(
                "source\t{i}\t{}\t{m}\t{n}\t{}\t{:?}",
                if head == HeadSide::Batter { "batter" } else { "resonant" },
                if sine { "sin" } else { "cos" },
                shape.zero
            ));
            line(format!("probe\t{i}\t0.31\t0.12\t{:?}", shape.value(0.31, 0.12, sine)));
        }
        let last = self.last_bin();
        let tw = twiddles(self.n);
        let probe = TAU * 0.5 * self.pass / self.rate;
        for i in 0..=self.sources.len() {
            let complement = i == self.sources.len();
            let integer = if complement { 0 } else { self.integer[i] };
            let r: Vec<(f64, f64)> =
                (1..=last).map(|b| if complement { complement_row(self.rate, self.n, b) } else { rows(i, b) }).collect();
            let kernel = Kernel { integer, taps: reconstruct(&r, integer, self.latency, self.n, self.rate, (self.pass, self.stop), &tw) };
            let (e, c, (hr, hi)) = kernel_stats(&kernel, probe);
            line(format!(
                "kernel\t{i}\t{integer}\t{}\t{e:?}\t{c:?}\t{hr:?}\t{hi:?}\t{}",
                self.n / 2,
                if complement { "complement" } else { "shape" }
            ));
        }
        for i in 0..self.sources.len() {
            for b in 1..=last {
                let (re, im) = rows(i, b);
                line(format!("row\t{i}\t{b}\t{re:?}\t{im:?}"));
            }
        }
        line(format!("end\t{last}"));
        t
    }
}

/// A small drum on the 12-inch tom's card (48 kHz, 6 mm grid), so the tests run in seconds.
fn test_drum() -> DrumSpec {
    DrumSpec {
        name: "test drum",
        radius: 0.06,
        depth: 0.05,
        grid: 0.006,
        strike_xy: (0.013, 0.004),
        mic: [0.2, 0.1, 0.5],
        floor: Some(floor_z(-0.3)),
        provenance: "closed-body observer test drum",
        ..kit::tom(12.0)
    }
}

fn shapes() -> Vec<(HeadSide, usize, usize, bool)> {
    let mut out = Vec::new();
    for head in [HeadSide::Batter, HeadSide::Resonant] {
        out.extend([(head, 0, 1, false), (head, 0, 2, false), (head, 1, 1, false), (head, 1, 1, true)]);
    }
    out
}

fn heads(ins: &Instrument) -> (usize, usize) {
    let find = |suffix: &str| ins.bodies.iter().position(|b| b.name.ends_with(suffix)).unwrap();
    (find(" batter"), find(" resonant"))
}

/// T2: driving a kernel's input with a sine through the observer's own
/// history replay gives the table's transfer (a delayed difference here) with the
/// latency, within 0.05 dB and 0.5 degrees, at in-band frequencies off the
/// FFT grid; and the complement kernel is 2 i tan(omega / 2) behind the latency.
#[test]
fn t2_replay_matches_the_table_transfer() {
    let rate = 48_000.0;
    let delay = 37.25; // samples, 30 of them integer
    let s = Synthetic {
        radius: 0.1,
        depth: 0.1,
        mic: [0.0, 0.0, 1.0],
        floor: ([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]),
        gas: (1.2, 343.0),
        rate,
        pass: 3000.0,
        stop: 6000.0,
        latency: 128,
        n: 2048,
        sources: vec![(HeadSide::Batter, 0, 1, false)],
        integer: vec![30],
    };
    // A delay of `delay` samples times 2 i sin(omega / 2), a half-sample difference: it vanishes
    // at DC, as every physical row does (velocity to pressure), because the bake zeroes bin 0. A
    // unit-gain delay is not a table the bake can make: its missing DC bin leaves a 1/n step
    // over the n/2 taps, 0.035 of the transfer at 211.7 Hz, which this test first measured.
    // Producer rows are conjugates (exp(-i omega t)): the delay is e^{+i omega delay}, the
    // difference -2 i sin(omega / 2).
    let text = s.write(&|_, b| {
        let w = TAU * b as f64 / s.n as f64;
        let (sn, cs) = det::sin_cos(w * delay);
        let g = 2.0 * det::sin(0.5 * w);
        (g * sn, -g * cs)
    });
    let table = ClosedBodyTable::parse(&text).unwrap();
    // (kernel, expected gain and phase at omega [rad/sample])
    let delayed = |w: f64| (2.0 * det::sin(0.5 * w), -w * (delay + 128.0) + 0.5 * PI);
    let complement = |w: f64| (2.0 * det::tan(0.5 * w), -w * 128.0 + 0.5 * PI);
    let cases: [(&Kernel, &dyn Fn(f64) -> (f64, f64)); 2] = [(&table.kernels[0], &delayed), (&table.complement, &complement)];
    for (kernel, expected) in cases {
        let reversed: Vec<f64> = kernel.taps.iter().rev().copied().collect();
        for f in [211.7, 1234.5, 2777.0] {
            let w = TAU * f / rate;
            let mut h = History::new(kernel.integer + kernel.taps.len() + 1);
            // least squares y = a sin(w t) + b cos(w t) = A sin(w t + phi) over the steady state
            let (mut ss, mut sc, mut cc, mut ys, mut yc) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let warm = kernel.integer + kernel.taps.len();
            for t in 0..warm + 8192 {
                let (sn, cs) = det::sin_cos(w * t as f64);
                h.write(t, sn);
                let y = h.dot(t, kernel.integer, &reversed);
                if t >= warm {
                    (ss, sc, cc, ys, yc) = (ss + sn * sn, sc + sn * cs, cc + cs * cs, ys + y * sn, yc + y * cs);
                }
            }
            let det2 = ss * cc - sc * sc;
            let (a, b) = ((ys * cc - yc * sc) / det2, (yc * ss - ys * sc) / det2);
            let (gain, phase) = expected(w);
            let level = 20.0 * (a.hypot(b) / gain).log10();
            let dphi = (b.atan2(a) - phase + PI).rem_euclid(TAU) - PI;
            eprintln!("T2 replay {f} Hz: level {level:+.5} dB re the table, phase {:+.5} deg", dphi.to_degrees());
            assert!(level.abs() <= 0.05 && dphi.to_degrees().abs() <= 0.5, "T2 replay at {f} Hz");
        }
    }
}

/// T2: with a table equal to the baffled path of its retained shapes, the
/// observer's output is today's output delayed by the latency: the closed-body
/// band and the complement sum to a pure delay.
#[test]
fn t2_complement_with_a_rayleigh_table_is_a_pure_delay() {
    let spec = test_drum();
    let hits = [(0.0005, 2.0)];
    let piece = kit::drum(&spec, &hits).unwrap();
    let control = piece.instrument.clone();
    let mut ins = piece.instrument;
    let (batter, resonant) = heads(&ins);
    let syn = Synthetic::for_drum(&ins, &spec, shapes());
    let zero = ClosedBodyTable::parse(&syn.write(&|_, _| (0.0, 0.0))).unwrap();
    ins.attach_closed_body(0, &zero, batter, resonant).unwrap();
    let _ = ins.image(); // prepares the taps and the observer's aggregated responses
    let rho = ins.mics[0].closed.as_ref().unwrap().rho.clone();
    let k = 1.0 / ins.rate;
    let rows = |s: usize, b: usize| {
        let w = TAU * b as f64 / syn.n as f64;
        let (d0, r) = &rho[s];
        let (mut re, mut im) = (0.0, 0.0);
        for (j, &v) in r.iter().enumerate() {
            let (sn, cs) = det::sin_cos(-w * (d0 + j) as f64);
            re += v * cs;
            im += v * sn;
        }
        // times 2 i k tan(w/2): the centred acceleration over the centred rate
        let g = 2.0 * k * det::tan(0.5 * w);
        let (tr, ti) = (-im * g, re * g);
        (tr, -ti) // the producer's row is the conjugate
    };
    let table = ClosedBodyTable::parse(&syn.write(&rows)).unwrap();
    ins.attach_closed_body(0, &table, batter, resonant).unwrap();
    let mut today = control;
    let steps = 3000;
    let (mut a, mut b) = (vec![0.0; 1], vec![0.0; 1]);
    let (mut y, mut y0) = (Vec::new(), Vec::new());
    for _ in 0..steps {
        ins.step(&mut a);
        today.step(&mut b);
        y.push(a[0]);
        y0.push(b[0]);
    }
    let d = syn.latency;
    let peak = y0.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let worst = (d..steps).map(|t| (y[t] - y0[t - d]).abs()).fold(0.0_f64, f64::max);
    eprintln!("T2 complement: peak {peak:.3e} Pa, largest difference from the delayed baffled output {worst:.3e} ({:.2e} of peak)", worst / peak);
    assert!(peak > 0.0 && worst <= 1e-6 * peak, "T2 complement: {worst:e} of {peak:e}");
    assert!(y[..d].iter().all(|v| v.abs() <= 1e-6 * peak), "nothing before the latency");
}

/// N2: one row's sign flipped (a resonant shape heard inverted) no longer
/// reconstructs fs-bem's kernel statistics, and the table is refused.
#[test]
fn n2_a_tampered_row_is_refused() {
    let spec = test_drum();
    let piece = kit::drum(&spec, &[(0.0005, 2.0)]).unwrap();
    let syn = Synthetic::for_drum(&piece.instrument, &spec, shapes());
    let text = syn.write(&|s, b| (1e-3 * (s + 1) as f64 / b as f64, 2e-4));
    assert!(ClosedBodyTable::parse(&text).is_ok());
    let row = text.lines().find(|l| l.starts_with("row\t5\t3\t")).unwrap().to_string();
    let cols: Vec<&str> = row.split('\t').collect();
    let flipped = format!("row\t5\t3\t{:?}\t{:?}", -cols[3].parse::<f64>().unwrap(), -cols[4].parse::<f64>().unwrap());
    let err = ClosedBodyTable::parse(&text.replace(&row, &flipped)).unwrap_err();
    assert!(err.contains("reconstructs differently"), "{err}");
}

/// N3: every header the table was baked for must match the instrument.
#[test]
fn n3_any_header_mismatch_refuses() {
    let spec = test_drum();
    let piece = kit::drum(&spec, &[(0.0005, 2.0)]).unwrap();
    let base = Synthetic::for_drum(&piece.instrument, &spec, shapes());
    let rows = |s: usize, b: usize| (1e-3 * (s + 1) as f64 / b as f64, 2e-4);
    let (batter, resonant) = heads(&piece.instrument);
    let attach = |syn: &Synthetic| -> Result<String, String> {
        let table = ClosedBodyTable::parse(&syn.write(&rows))?;
        let mut ins = piece.instrument.clone();
        ins.attach_closed_body(0, &table, batter, resonant)
    };
    attach(&base).unwrap();
    let mut cases: Vec<(&str, Synthetic)> = Vec::new();
    cases.push(("rate", Synthetic { rate: 2.0 * base.rate, ..base.clone() }));
    cases.push(("gas", Synthetic { gas: (base.gas.0 * 1.01, base.gas.1), ..base.clone() }));
    cases.push(("microphone", Synthetic { mic: [base.mic[0] + 0.01, base.mic[1], base.mic[2]], ..base.clone() }));
    cases.push(("floor", Synthetic { floor: ([0.0, 0.0, -0.31], base.floor.1), ..base.clone() }));
    cases.push(("depth", Synthetic { depth: base.depth * 1.1, ..base.clone() }));
    cases.push(("radius", Synthetic { radius: base.radius * 1.2, ..base.clone() }));
    for (what, syn) in &cases {
        let err = attach(syn).unwrap_err();
        eprintln!("N3 {what}: {err}");
    }
    // The kit's own check also binds the card's radius and depth.
    let table = base.write(&rows);
    let mut p = kit::drum(&spec, &[(0.0005, 2.0)]).unwrap();
    assert!(kit::attach_closed_body(&mut p, &DrumSpec { radius: 0.061, ..spec }, &table).is_err());
    assert!(kit::attach_closed_body(&mut p, &spec, &table).is_ok());
    // A table of another format or route refuses at parse.
    assert!(ClosedBodyTable::parse(&table.replacen(TABLE_FORMAT, "frankensim.vkplate-head-radiation.v1", 1)).is_err());
    assert!(ClosedBodyTable::parse(&table.replacen(ROUTE, "planar-baffled-rayleigh", 1)).is_err());
    // Attaching after the first step refuses: the rings are already running.
    let mut stepped = piece.instrument.clone();
    stepped.step(&mut [0.0]);
    let parsed = ClosedBodyTable::parse(&table).unwrap();
    assert!(stepped.attach_closed_body(0, &parsed, batter, resonant).is_err());
}

/// The dipole's length from the committed tables, on z-floors and the kick's y-floor, against an independent axisymmetric
/// boundary-integral solve of the closed cylinder in free field (its net-flow power reaches 0.9987 of the monopole law).
/// The script's read of the table's own rows is the first number; the solve's power-equivalent length is the second.
#[test]
fn the_tables_give_the_independent_solves_dipole_length() {
    let tables = [
        ("snare-off", include_str!("../data/closed-body/snare-off.tsv"), 0.3694, 0.3705, 0.01),
        ("tom-12", include_str!("../data/closed-body/tom-12.tsv"), 0.3675, 0.3666, 0.01),
        ("tom-14", include_str!("../data/closed-body/tom-14.tsv"), 0.6017, 0.5764, 0.05),
        ("kick", include_str!("../data/closed-body/kick.tsv"), 0.7882, 0.7436, 0.07),
    ];
    for (name, text, script, solve, limit) in tables {
        let d = ClosedBodyTable::parse(text).unwrap().dipole_length().unwrap();
        eprintln!("{name}: D {d:.4} m (script's read of the rows {script}, solve {solve})");
        assert!((d - script).abs() < 1e-3, "{name}: {d} m, the rows read {script} m");
        assert!((d / solve - 1.0).abs() < limit, "{name}: {d} m against the solve's {solve} m");
    }
    // Planted negative: the heads' sources swapped make a negative length, which is refused, not used.
    let swapped: String = tables[0]
        .1
        .lines()
        .map(|line| match line.starts_with("source\t") {
            true if line.contains("\tbatter\t") => line.replace("\tbatter\t", "\tresonant\t"),
            true => line.replace("\tresonant\t", "\tbatter\t"),
            false => line.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n");
    let error = ClosedBodyTable::parse(&swapped).unwrap().dipole_length().unwrap_err();
    assert!(error.contains("not a real positive length"), "{error}");
}
