//! Finite-footprint viscous dampers on the existing moving-bridge strings.
//!
//! This is a spatial refinement of the point-drag image, NOT a wool-felt impact
//! or falling damper action. For a pad [a,b], c [N s/m] is its TOTAL drag on
//! EACH speaking string: power = c/(b-a) * integral(v(x)^2 dx). Averaging v
//! first would incorrectly cancel opposite-moving portions of the string.
//! Positive midpoint weights retain the full off-diagonal modal damping and
//! the moving endpoint. Symmetric composition of exact rank-one flows is
//! passive and second-order in time, not the exact full damping exponential.
use super::{Bank, Course};
use fs_math::det;
use std::{collections::BTreeMap, io::Read, ops::Range};

const HEADER: &str = "frankensim-piano-dampers-v1";
const MAX_BYTES: u64 = 1024 * 1024;
const MAX_CELLS: usize = 128;
const MAX_TERMS: usize = 4_000_000;

#[derive(Clone, Copy, Debug)]
struct Pad {
    start_m: f64,
    end_m: f64,
    drag_ns_m: f64,
}

/// Every supplied course has either a finite pad or an explicit `free` row.
/// Missing, duplicate and extra keys never acquire an estimated damper.
#[derive(Debug)]
pub struct Specification {
    pads: BTreeMap<u8, Option<Pad>>,
}
impl Specification {
    pub fn read(text: &str, courses: &[Course]) -> Result<Self, String> {
        if text.len() as u64 > MAX_BYTES { return Err("damper file exceeds 1 MiB".into()); }
        let mut header = false;
        let mut pads = BTreeMap::new();
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            if !header {
                if row != HEADER { return Err(format!("line {}: expected {HEADER}", line + 1)); }
                header = true;
                continue;
            }
            let parse = (|| -> Result<(u8, Option<Pad>), String> {
                let fields: Vec<_> = row.split(',').map(str::trim).collect();
                if fields.len() != 2 && fields.len() != 5 {
                    return Err("expected free,key or pad,key,start_m,end_m,drag_ns_m".into());
                }
                let key: u8 = fields[1].parse().map_err(|_| "invalid damper key")?;
                if !courses.iter().any(|c| c.midi == key) { return Err("damper key absent from scale".into()); }
                let pad = match (fields[0], fields.len()) {
                    ("free", 2) => None,
                    ("pad", 5) => Some(Pad {
                        start_m: fields[2].parse().map_err(|_| "invalid pad start")?,
                        end_m: fields[3].parse().map_err(|_| "invalid pad end")?,
                        drag_ns_m: fields[4].parse().map_err(|_| "invalid pad drag")?,
                    }),
                    _ => return Err("unknown damper row or wrong field count".into()),
                };
                Ok((key, pad))
            })();
            let (key, pad) = parse.map_err(|e| format!("line {}: {e}", line + 1))?;
            if pads.insert(key, pad).is_some() { return Err(format!("duplicate damper key {key}")); }
        }
        if !header { return Err("missing damper header".into()); }
        let spec = Self { pads };
        spec.validate(courses)?;
        Ok(spec)
    }

    pub fn load(path: &str, courses: &[Course]) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?
            .take(MAX_BYTES + 1).read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
        Self::read(&text, courses)
    }

    /// Explicit estimates, not Steinway dimensions or identified felt losses.
    /// Same 35%-length centre and 0.4 N s/m per-string drag as the point image;
    /// width = min(40 mm, 8% of L). The existing estimated upper break is 88.
    pub fn estimated(courses: &[Course]) -> Result<Self, String> {
        let pads = courses.iter().map(|c| {
            let width = (0.08 * c.length_m).min(0.04);
            (c.midi, (c.midi <= 88).then_some(Pad {
                start_m: 0.35 * c.length_m - 0.5 * width,
                end_m: 0.35 * c.length_m + 0.5 * width,
                drag_ns_m: 0.4,
            }))
        }).collect();
        let spec = Self { pads };
        spec.validate(courses)?;
        Ok(spec)
    }

    fn validate(&self, courses: &[Course]) -> Result<(), String> {
        if courses.is_empty() || courses.len() > 88 || self.pads.len() != courses.len()
            || courses.iter().enumerate().any(|(i,c)| courses[..i].iter().any(|p| p.midi == c.midi)) {
            return Err("damper rows must cover a nonempty unique scale exactly".into());
        }
        for c in courses {
            c.validate()?;
            let pad = self.pads.get(&c.midi).ok_or_else(|| format!("missing damper key {}", c.midi))?;
            if let Some(p) = pad {
                if [p.start_m, p.end_m, p.drag_ns_m].iter().any(|x| !x.is_finite())
                    || p.start_m <= 0.0 || p.end_m <= p.start_m || p.end_m >= c.length_m
                    || p.drag_ns_m <= 0.0 || p.drag_ns_m > 1e6 {
                    return Err(format!("key {}: pad must lie strictly inside its speaking string with finite drag in (0,1e6] N s/m", c.midi));
                }
            }
        }
        Ok(())
    }
}

struct Point { shape: Vec<f64>, lift: f64 }
struct StringPad {
    course: usize,
    modes: Range<usize>,
    bridge: Vec<f64>,
    drag_ns_m: f64,
    points: Vec<Point>,
}

/// `apply_reduced`'s constants: per pad |bridge|^2, per station its norm |shape|^2 + lift^2 |bridge|^2,
/// and per station expm1(-drag norm dt) / norm for the engagement and interval they were last built for
/// (they change only with the pedal).
struct Reduced {
    bridge_norm2: Vec<f64>,
    norm: Vec<Vec<f64>>,
    factor: Vec<Vec<f64>>,
    key: Option<(u64, u64)>,
    /// `apply_superposed`'s per-direction unit responses (built with `factor`) and the pads' bridge Gram
    /// matrix b_p . b_q (constant).
    unit: [Vec<Unit>; 2],
    gram: Vec<f64>,
}

/// One pad's chain driven by beta = 1 with its string modes at rest, in one direction: the modes it
/// leaves behind, the board amount it moves (sum of lift change), its loss, and per station (in pass
/// order) the speed and change, for the cross term of the loss.
#[derive(Clone, Default)]
struct Unit {
    w: Vec<f64>,
    alpha: f64,
    loss: f64,
    speed: Vec<f64>,
    change: Vec<f64>,
}

/// `Prepared::gpu_image`: the pads (string mode range, bridge row and |row|^2; per station its shape,
/// lift, norm and factor in station order) with their unit responses per direction (0 forward, 1
/// reverse: board amount moved, loss, and per station the speed and change in pass order), and the
/// bridge Gram matrix.
#[derive(Clone, Debug)]
#[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
pub struct DamperImage {
    pub board_start: usize,
    pub gram: Vec<f64>,
    pub pads: Vec<PadImage>,
}
#[derive(Clone, Debug)]
#[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
pub struct PadImage {
    pub course: usize,
    pub modes: Range<usize>,
    pub bridge: Vec<f64>,
    pub bridge_norm2: f64,
    pub shape: Vec<Vec<f64>>,
    pub lift: Vec<f64>,
    pub norm: Vec<f64>,
    pub factor: Vec<f64>,
    pub unit: [UnitImage; 2],
}
#[derive(Clone, Debug)]
#[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
pub struct UnitImage {
    pub alpha: f64,
    pub loss: f64,
    pub speed: Vec<f64>,
    pub change: Vec<f64>,
}

/// Frozen displacement/force projections in the Bank's loaded coordinates.
/// Does not contain, duplicate, reset or truncate any resonator state.
pub struct Prepared {
    pads: Vec<StringPad>,
    dimension: usize,
    board_start: usize,
    cells: usize,
    reduced: Option<Reduced>,
}
impl Prepared {
    pub fn new(spec: &Specification, courses: &[Course], bank: &Bank) -> Result<Self, String> {
        Self::new_with_transverse_drag(spec,courses,bank,None)
    }

    /// Keep the original vertical pad law and explicitly scale its lateral
    /// drag for the second transverse direction. Neither direction of a duplex
    /// is touched. Missing lateral material response is refused, not guessed.
    pub fn new_with_transverse_drag(spec: &Specification, courses: &[Course], bank: &Bank,
        lateral_ratios: Option<&[f64]>) -> Result<Self, String> {
        spec.validate(courses)?;
        if bank.has_secondary_polarization()!=lateral_ratios.is_some()
            || lateral_ratios.is_some_and(|r|r.len()!=courses.len()
                || r.iter().any(|v|!v.is_finite() || !(0.0..=10.0).contains(v))) {
            return Err("damper preparation needs the complete explicit transverse drag selection".into());
        }
        let board_start = bank.modes.len();
        let dimension = board_start + bank.board_count;
        if bank.v.len() != dimension || bank.q.len() != dimension {
            return Err("damper preparation needs the original bank coordinates".into());
        }
        let mut pads = Vec::new();
        let mut terms = 0usize;
        let mut cells = 0usize;
        for (si, s) in bank.strings.iter().enumerate() {
            // Dampers touch speaking strings only, never duplex segments.
            if s.duplex { continue; }
            let c = courses.get(s.course).ok_or("invalid damper course address")?;
            let Some(p) = spec.pads[&c.midi] else { continue; };
            let ratio=if s.polarization==0 {1.}else {lateral_ratios.ok_or("missing lateral pad response")?[s.course]};
            let drag_ns_m=p.drag_ns_m*ratio;
            if !drag_ns_m.is_finite() || drag_ns_m>1e6 {return Err("directional pad drag exceeds 1e6 N s/m".into());}
            if drag_ns_m==0. {continue;} // explicit zero drag, not an absent vibrating string
            if s.modes.is_empty() || s.modes.end > board_start || s.bridge.len() != bank.board_count
                || s.bridge.iter().any(|x| !x.is_finite())
                || bank.modes[s.modes.clone()].iter().any(|m| m.string != si || !m.beta.is_finite()) {
                return Err("damper projection does not match the string/bridge basis".into());
            }
            let span = (p.end_m - p.start_m) / c.length_m;
            // v^2 contains spatial frequencies up to twice the highest partial.
            // At least four points per shortest product wavelength, and eight
            // per footprint. This is a resolution guard, not an error bound.
            let count = (4.0 * s.modes.len() as f64 * span).ceil().max(8.0) as usize;
            if count > MAX_CELLS { return Err("damper spatial budget exceeded; refine the declared footprint/basis study".into()); }
            terms = terms.checked_add(count * (s.modes.len() + bank.board_count))
                .ok_or("damper setup size overflow")?;
            if terms > MAX_TERMS { return Err("damper projection exceeds four million terms".into()); }
            let root_mass = det::sqrt(c.modal_mass_kg());
            let mut points = Vec::with_capacity(count);
            let mut previous = p.start_m / c.length_m;
            for i in 0..count {
                let fraction = p.start_m / c.length_m + span * (i as f64 + 0.5) / count as f64;
                if !fraction.is_finite() || fraction <= previous || fraction >= p.end_m / c.length_m {
                    return Err("damper footprint has no representable spatial resolution".into());
                }
                previous = fraction;
                let shape: Vec<f64> = (1..=s.modes.len()).map(|n|
                    det::sin(n as f64 * std::f64::consts::PI * fraction) / root_mass).collect();
                // z=q_fixed+beta*b: y=phi*z+(x/L-sum(phi*beta))*b.
                // Omitting this lift silently applies damping in the wrong
                // coordinate system and loses the reciprocal bridge reaction.
                let lift = fraction - shape.iter().zip(&bank.modes[s.modes.clone()])
                    .map(|(g,m)| g*m.beta).sum::<f64>();
                let norm = shape.iter().map(|g| g*g).sum::<f64>()
                    + s.bridge.iter().map(|b| (lift*b).powi(2)).sum::<f64>();
                if !norm.is_finite() || norm <= 0.0 || !lift.is_finite()
                    || !(drag_ns_m * norm).is_finite() {
                    return Err("damper projection or damping rate is unrepresentable".into());
                }
                points.push(Point { shape, lift });
            }
            cells += count;
            pads.push(StringPad { course: s.course, modes: s.modes.clone(), bridge: s.bridge.clone(),
                drag_ns_m, points });
        }
        Ok(Self { pads, dimension, board_start, cells, reduced: None })
    }

    /// Select `apply_reduced` for this damper set (`--fast-dampers`). The default `apply` keeps its
    /// order of arithmetic; this one carries the same flows on fewer coordinates.
    pub fn use_reduced(&mut self) {
        let bridge_norm2: Vec<f64> = self.pads.iter().map(|p| p.bridge.iter().map(|b| b * b).sum()).collect();
        let norm: Vec<Vec<f64>> = self.pads.iter().zip(&bridge_norm2).map(|(p, &nb2)| {
            p.points.iter().map(|pt| pt.shape.iter().map(|g| g * g).sum::<f64>() + pt.lift * pt.lift * nb2).collect()
        }).collect();
        let factor = norm.iter().map(|n| vec![0.0; n.len()]).collect();
        let np = self.pads.len();
        let mut gram = vec![0.0; np * np];
        for p in 0..np {
            for q in 0..np {
                gram[p * np + q] = self.pads[p].bridge.iter().zip(&self.pads[q].bridge).map(|(a, b)| a * b).sum();
            }
        }
        self.reduced = Some(Reduced { bridge_norm2, norm, factor, key: None, unit: [vec![], vec![]], gram });
    }

    /// Rebuild the station factors (and the unit responses) when the pedal or the interval changed.
    fn refresh(&mut self, engagement: f64, dt: f64) {
        let Some(red) = self.reduced.as_mut() else { return };
        let key = (engagement.to_bits(), dt.to_bits());
        if red.key == Some(key) { return; }
        for (pi, pad) in self.pads.iter().enumerate() {
            let drag = pad.drag_ns_m * engagement / pad.points.len() as f64;
            for (f, &n) in red.factor[pi].iter_mut().zip(&red.norm[pi]) {
                *f = if drag == 0.0 { 0.0 } else { det::expm1(-drag*n*(0.5*dt)) / n };
            }
        }
        for (d, units) in red.unit.iter_mut().enumerate() {
            *units = self.pads.iter().enumerate().map(|(pi, pad)| {
                let n = pad.points.len();
                let mut u = Unit { w: vec![0.0; pad.modes.len()], ..Unit::default() };
                let mut beta = 1.0;
                for t in 0..n {
                    let i = if d == 0 { t } else { n - 1 - t };
                    let pt = &pad.points[i];
                    let speed = pt.shape.iter().zip(&u.w).map(|(g, x)| g * x).sum::<f64>() + pt.lift * beta;
                    let change = red.factor[pi][i] * speed;
                    for (x, g) in u.w.iter_mut().zip(&pt.shape) { *x += g * change; }
                    beta += pt.lift * red.bridge_norm2[pi] * change;
                    u.alpha += pt.lift * change;
                    u.loss += -change*speed-0.5*change*change*red.norm[pi][i];
                    u.speed.push(speed);
                    u.change.push(change);
                }
                u
            }).collect();
        }
        red.key = Some(key);
    }

    /// `apply_reduced` reorganised so that every pad's work is independent (for parallel hardware), with
    /// the same flows and loss in exact arithmetic. A pad's chain is linear in its string modes v and the
    /// scalar beta = b . v_board it starts from, so it is the chain run from (v, 0) plus beta times its
    /// unit response (the chain run from (0, 1), fixed until the pedal moves). Pads touch disjoint string
    /// modes and couple only through the board: the beta pad p starts from is b_p . v_board plus
    /// sum over the pads before it of (b_p . b_q) times the board amount q moved, a recurrence over the
    /// bridge Gram matrix. The loss, quadratic in (v, beta), splits the same way.
    pub fn apply_superposed(&mut self, velocity: &mut [f64], dt: f64, sustain: f64,
        lifted: impl Fn(usize) -> bool) -> Result<f64, &'static str> {
        if velocity.len() != self.dimension || velocity.iter().any(|x| !x.is_finite())
            || !dt.is_finite() || dt <= 0.0 || dt > 1.0
            || !sustain.is_finite() || !(0.0..=1.0).contains(&sustain) {
            return Err("invalid damper velocity, duration or pedal control");
        }
        if sustain == 1.0 { return Ok(0.0); }
        if self.reduced.is_none() { return Err("reduced dampers were not selected"); }
        self.refresh((1.0-sustain).powi(2), dt);
        let red = self.reduced.as_ref().expect("selected above");
        let np = self.pads.len();
        let active: Vec<usize> = (0..np).filter(|&p| !lifted(self.pads[p].course)).collect();
        let (strings, board) = velocity.split_at_mut(self.board_start);
        let mut loss = 0.0;
        for d in 0..2 {
            let order: Vec<usize> = if d == 0 { active.clone() } else { active.iter().rev().copied().collect() };
            // every pad from beta = 0 (independent), with the cross term of its loss
            let mut a = vec![0.0; np];
            let mut l01 = vec![0.0; np];
            for &p in &order {
                let pad = &self.pads[p];
                let unit = &red.unit[d][p];
                let v = &mut strings[pad.modes.clone()];
                let n = pad.points.len();
                // the run from beta = 0 still moves its own beta (by lift |b|^2 change per station)
                let mut own = 0.0;
                for t in 0..n {
                    let i = if d == 0 { t } else { n - 1 - t };
                    let pt = &pad.points[i];
                    let speed = pt.shape.iter().zip(v.iter()).map(|(g, x)| g * x).sum::<f64>() + pt.lift * own;
                    let change = red.factor[p][i] * speed;
                    for (x, g) in v.iter_mut().zip(&pt.shape) { *x += g * change; }
                    own += pt.lift * red.bridge_norm2[p] * change;
                    a[p] += pt.lift * change;
                    loss += -change*speed-0.5*change*change*red.norm[p][i];
                    l01[p] += -(change*unit.speed[t] + unit.change[t]*speed) - change*unit.change[t]*red.norm[p][i];
                }
            }
            // the betas, pad by pad in pass order: the board as it was, plus the pads already through
            let mut beta = vec![0.0; np];
            for (k, &p) in order.iter().enumerate() {
                let mut x: f64 = self.pads[p].bridge.iter().zip(board.iter()).map(|(b, v)| b * v).sum();
                for &q in &order[..k] {
                    x += red.gram[p * np + q] * (a[q] + red.unit[d][q].alpha * beta[q]);
                }
                beta[p] = x;
            }
            for &p in &order {
                let unit = &red.unit[d][p];
                let b = beta[p];
                for (x, w) in strings[self.pads[p].modes.clone()].iter_mut().zip(&unit.w) { *x += w * b; }
                loss += b * l01[p] + b * b * unit.loss;
                let moved = a[p] + unit.alpha * b;
                for (x, g) in board.iter_mut().zip(&self.pads[p].bridge) { *x += g * moved; }
            }
        }
        Ok(loss)
    }

    /// Whether `use_reduced` was selected.
    pub fn is_reduced(&self) -> bool { self.reduced.is_some() }

    /// The superposed pass as plain data for the pedal value and interval of one damper call (selects
    /// the reduced coordinates if they were not): for an external stepper (render-gpu-piano).
    #[allow(dead_code)] // used by apps/instrument-lab/render-gpu-piano
    pub fn gpu_image(&mut self, sustain: f64, dt: f64) -> DamperImage {
        if self.reduced.is_none() { self.use_reduced(); }
        let engagement = (1.0-sustain).powi(2);
        self.refresh(engagement, dt);
        let red = self.reduced.as_ref().expect("selected above");
        DamperImage {
            board_start: self.board_start,
            gram: red.gram.clone(),
            pads: self.pads.iter().enumerate().map(|(p, pad)| PadImage {
                course: pad.course, modes: pad.modes.clone(), bridge: pad.bridge.clone(),
                bridge_norm2: red.bridge_norm2[p],
                shape: pad.points.iter().map(|pt| pt.shape.clone()).collect(),
                lift: pad.points.iter().map(|pt| pt.lift).collect(),
                norm: red.norm[p].clone(), factor: red.factor[p].clone(),
                unit: [0, 1].map(|d| {
                    let u = &red.unit[d][p];
                    UnitImage { alpha: u.alpha, loss: u.loss, speed: u.speed.clone(), change: u.change.clone() }
                }),
            }).collect(),
        }
    }

    /// `apply`'s palindrome of exact rank-one flows, carried on fewer coordinates. Every station of a pad
    /// moves the board along that pad's one bridge direction b, scaled by the station's lift, so a pad's
    /// run of stations needs only its string modes and the scalar beta = b . v_board: a station's speed is
    /// shape . v_string + lift beta, its update moves beta by lift |b|^2 change, and the board takes one
    /// update b sum(lift change) at the end of the run. The flows, their order and the reported loss are
    /// those of `apply` in exact arithmetic; the rounding differs. Per station the work falls from three
    /// passes over (string modes + board coordinates) to two over the string modes, and the norms and the
    /// expm1 factors are computed once per pedal value instead of at every station.
    pub fn apply_reduced(&mut self, velocity: &mut [f64], dt: f64, sustain: f64,
        lifted: impl Fn(usize) -> bool) -> Result<f64, &'static str> {
        if velocity.len() != self.dimension || velocity.iter().any(|x| !x.is_finite())
            || !dt.is_finite() || dt <= 0.0 || dt > 1.0
            || !sustain.is_finite() || !(0.0..=1.0).contains(&sustain) {
            return Err("invalid damper velocity, duration or pedal control");
        }
        if sustain == 1.0 { return Ok(0.0); }
        let engagement = (1.0-sustain).powi(2);
        if self.reduced.is_none() { return Err("reduced dampers were not selected"); }
        self.refresh(engagement, dt);
        let red = self.reduced.as_mut().expect("selected above");
        let (strings, board) = velocity.split_at_mut(self.board_start);
        let mut loss = 0.0;
        let mut run = |pi: usize, forward: bool| {
            let pad = &self.pads[pi];
            let v = &mut strings[pad.modes.clone()];
            let nb2 = red.bridge_norm2[pi];
            let mut beta: f64 = pad.bridge.iter().zip(board.iter()).map(|(b, x)| b * x).sum();
            let mut moved = 0.0;
            let n = pad.points.len();
            for t in 0..n {
                let i = if forward { t } else { n - 1 - t };
                let pt = &pad.points[i];
                let speed = pt.shape.iter().zip(v.iter()).map(|(g, x)| g * x).sum::<f64>() + pt.lift * beta;
                let change = red.factor[pi][i] * speed;
                for (x, g) in v.iter_mut().zip(&pt.shape) { *x += g * change; }
                beta += pt.lift * nb2 * change;
                moved += pt.lift * change;
                loss += -change*speed-0.5*change*change*red.norm[pi][i];
            }
            for (x, b) in board.iter_mut().zip(&pad.bridge) { *x += b * moved; }
        };
        for (pi, pad) in self.pads.iter().enumerate() {
            if !lifted(pad.course) { run(pi, true); }
        }
        for (pi, pad) in self.pads.iter().enumerate().rev() {
            if !lifted(pad.course) { run(pi, false); }
        }
        Ok(loss)
    }

    /// Number of damped scalar direction fields; a two-plane string may count twice.
    pub fn string_count(&self) -> usize { self.pads.len() }
    pub fn cell_count(&self) -> usize { self.cells }

    /// One symmetric damping interval. `lifted(course)` includes key hold and
    /// sostenuto latch; sustain retains the existing squared drag-travel map.
    /// No hot allocation. q and all hammer/material histories remain untouched.
    pub fn apply(&self, velocity: &mut [f64], dt: f64, sustain: f64,
        lifted: impl Fn(usize) -> bool) -> Result<f64, &'static str> {
        if velocity.len() != self.dimension || velocity.iter().any(|x| !x.is_finite())
            || !dt.is_finite() || dt <= 0.0 || dt > 1.0
            || !sustain.is_finite() || !(0.0..=1.0).contains(&sustain) {
            return Err("invalid damper velocity, duration or pedal control");
        }
        if sustain == 1.0 { return Ok(0.0); }
        let engagement = (1.0-sustain).powi(2);
        let mut loss = 0.0;
        let apply = |pad: &StringPad, point: &Point, velocity: &mut [f64]| {
            rank_one(velocity, pad.modes.clone(), self.board_start,
                |k| point.shape[k-pad.modes.start], &pad.bridge, point.lift,
                pad.drag_ns_m * engagement / pad.points.len() as f64, 0.5*dt)
        };
        // Reverse the ENTIRE point sequence, not merely points within each pad:
        // different strings share board coordinates, so their flows can fail to
        // commute. This palindrome retains second-order splitting for the set.
        for pad in &self.pads {
            if !lifted(pad.course) { for point in &pad.points { loss += apply(pad, point, velocity); } }
        }
        for pad in self.pads.iter().rev() {
            if !lifted(pad.course) { for point in pad.points.iter().rev() { loss += apply(pad, point, velocity); } }
        }
        Ok(loss)
    }
}

/// One owner for the point image AND each finite-footprint quadrature port.
/// Exact flow dv/dt=-drag*g*g^T*v for mass-normalized v; loss is its actual
/// kinetic-energy decrement. Callers provide validated disjoint string/board
/// ranges and nonnegative finite drag/time. Preserve the point image's order
/// of arithmetic, including expm1 for small losses.
#[allow(clippy::too_many_arguments)]
pub(super) fn rank_one(velocity: &mut [f64], modes: Range<usize>, board_start: usize,
    shape: impl Fn(usize) -> f64, bridge: &[f64], lift: f64, drag: f64, dt: f64) -> f64 {
    let mut norm = 0.0;
    let mut speed = 0.0;
    for k in modes.clone() { let g = shape(k); norm += g*g; speed += g*velocity[k]; }
    for (j,b) in bridge.iter().enumerate() { let g = lift*b; norm += g*g; speed += g*velocity[board_start+j]; }
    if norm == 0.0 || drag == 0.0 { return 0.0; }
    let change = det::expm1(-drag*norm*dt)*speed/norm;
    for k in modes { velocity[k] += shape(k)*change; }
    for (j,b) in bridge.iter().enumerate() { velocity[board_start+j] += lift*b*change; }
    -change*speed-0.5*change*change*norm
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::BoardMode;
    fn course() -> Course {
        Course { midi: 69, unison: 1, length_m: 1.0, linear_density_kg_m: 0.02,
            tension_n: 100.0, flexural_rigidity_nm2: 0.0, strike_fraction: 0.12,
            bonded_section: None,
            hammer_mass_kg: 0.008, felt_thickness_m: 0.008, felt_area_m2: 0.0001,
            duplex_length_m: 0.0, detune_cents: 0.0 }
    }
    fn bank(bridge: f64) -> Bank {
        Bank::new(&[course()], &[BoardMode { frequency_hz: 100.0, damping_ratio: 0.0,
            bridge: [bridge;88], volume: 0.1 }], 192_000, 21_600.0, 32, false).unwrap()
    }
    fn specification() -> Specification {
        Specification::read(&format!("{HEADER}\npad,69,0.33,0.37,0.4\n"), &[course()]).unwrap()
    }
    #[test]
    fn distributed_loss_reaches_a_partial_at_the_point_dampers_node() {
        let mut point = bank(0.0);
        let mut span = bank(0.0);
        point.v[19] = 0.2; span.v[19] = 0.2;
        let before = span.energy(); let q = span.q.clone();
        let legacy = point.damp_string(0, 0.4, 0.001);
        let prepared = Prepared::new(&specification(), &[course()], &span).unwrap();
        let loss = prepared.apply(&mut span.v, 0.001, 0.0, |_| false).unwrap();
        assert!(legacy.abs() < 1e-24, "sin(20*pi*0.35) is a point node");
        assert!(loss > 1e-5, "a pad must integrate speed squared, not square mean speed");
        assert_eq!(span.q, q);
        assert!((span.energy()+loss-before).abs() < 1e-12);
    }
    #[test]
    fn spatial_projection_includes_the_loaded_moving_endpoint() {
        let b = bank(0.1);
        let d = Prepared::new(&specification(), &[course()], &b).unwrap();
        let pad = &d.pads[0]; let bridge_motion = pad.bridge[0]*0.03;
        for (i,p) in pad.points.iter().enumerate() {
            let fraction = 0.33 + 0.04*(i as f64+0.5)/pad.points.len() as f64;
            let observed = p.shape.iter().zip(&b.modes[pad.modes.clone()])
                .map(|(g,m)| g*m.beta*bridge_motion).sum::<f64>()+p.lift*bridge_motion;
            assert!((observed-fraction*bridge_motion).abs() < 1e-14);
        }
    }
    #[test]
    fn pedal_key_and_explicit_free_rows_do_not_erase_vibration() {
        let mut b = bank(0.1); b.v.fill(0.03);
        let d = Prepared::new(&specification(), &[course()], &b).unwrap();
        let original = b.v.clone();
        assert_eq!(d.apply(&mut b.v, 0.001, 1.0, |_| false).unwrap(), 0.0);
        assert_eq!(d.apply(&mut b.v, 0.001, 0.0, |_| true).unwrap(), 0.0);
        assert_eq!(b.v, original);
        let free = Specification::read(&format!("{HEADER}\nfree,69\n"), &[course()]).unwrap();
        let free = Prepared::new(&free, &[course()], &b).unwrap();
        assert_eq!(free.string_count(), 0);
        assert_eq!(free.apply(&mut b.v, 0.001, 0.0, |_| false).unwrap(), 0.0);
        assert_eq!(b.v, original);
        let before = b.energy();
        let loss = d.apply(&mut b.v, 0.001, 0.5, |_| false).unwrap();
        assert!(loss > 0.0 && (b.energy()+loss-before).abs() < 1e-12);
    }
    #[test]
    fn every_unison_member_is_damped_but_duplex_modes_are_not_directly_touched() {
        let c = Course { unison: 3, duplex_length_m: 0.3, ..course() };
        let mut b = Bank::new(&[c], &[BoardMode { frequency_hz: 100.0, damping_ratio: 0.0,
            bridge: [0.0;88], volume: 0.1 }], 192_000, 21_600.0, 32, false).unwrap();
        let d = Prepared::new(&Specification::estimated(&[c]).unwrap(), &[c], &b).unwrap();
        assert_eq!(d.string_count(), 3); assert!(d.cell_count() >= 24);
        b.v.fill(0.03); let original = b.v.clone();
        d.apply(&mut b.v, 0.001, 0.0, |_| false).unwrap();
        for s in &b.strings {
            if s.contact.is_none() { assert_eq!(&b.v[s.modes.clone()], &original[s.modes.clone()]); }
            else { assert_ne!(&b.v[s.modes.clone()], &original[s.modes.clone()]); }
        }
    }
    #[test]
    fn invalid_specs_and_controls_refuse_without_defaults_or_state_changes() {
        for rows in ["", "free,68", "free,69\nfree,69", "pad,69,0,0.1,1",
            "pad,69,0.1,1,1", "pad,69,0.3,0.2,1", "pad,69,0.2,0.3,-1",
            "pad,69,NaN,0.3,1", "pad,69,0.2,0.3,inf", "free,69,0,0,0"] {
            assert!(Specification::read(&format!("{HEADER}\n{rows}"), &[course()]).is_err(), "{rows}");
        }
        let mut b = bank(0.1); b.v.fill(0.03); let saved = b.v.clone();
        let d = Prepared::new(&specification(), &[course()], &b).unwrap();
        for (dt,pedal) in [(0.0,0.0),(0.001,f64::NAN),(0.001,-0.1),(f64::INFINITY,0.0)] {
            assert!(d.apply(&mut b.v, dt, pedal, |_| false).is_err()); assert_eq!(b.v, saved);
        }
    }
    #[test]
    fn reduced_coordinates_carry_the_same_flows_and_losses_as_the_full_chain() {
        // two pads sharing a two-mode board through different bridge rows: the board couples them
        let low = Course { midi: 57, length_m: 1.4, ..course() };
        let courses = [low, course()];
        let rows: [f64; 88] = std::array::from_fn(|i| 0.02*i as f64 - 0.3);
        let board = [BoardMode { frequency_hz: 100.0, damping_ratio: 0.0, bridge: [0.1;88], volume: 0.1 },
            BoardMode { frequency_hz: 173.0, damping_ratio: 0.0, bridge: rows, volume: 0.05 }];
        let b = Bank::new(&courses, &board, 192_000, 21_600.0, 32, false).unwrap();
        let spec = Specification::read(&format!("{HEADER}\npad,57,0.40,0.50,0.4\npad,69,0.33,0.37,0.4\n"), &courses).unwrap();
        let full = Prepared::new(&spec, &courses, &b).unwrap();
        let mut fast = Prepared::new(&spec, &courses, &b).unwrap();
        fast.use_reduced();
        let mut superposed = Prepared::new(&spec, &courses, &b).unwrap();
        superposed.use_reduced();
        assert!(fast.is_reduced() && !full.is_reduced());
        let v0: Vec<f64> = (0..b.v.len()).map(|i| ((i*37 % 11) as f64 - 5.0)*0.01).collect();
        for (pedal, lifted_course) in [(0.0, None), (0.4, None), (0.0, Some(0usize)), (0.0, Some(1)), (0.4, None), (1.0, None)] {
            let lifted = |ci: usize| Some(ci) == lifted_course;
            let (mut a, mut r, mut u) = (v0.clone(), v0.clone(), v0.clone());
            let la = full.apply(&mut a, 0.002, pedal, lifted).unwrap();
            let lr = fast.apply_reduced(&mut r, 0.002, pedal, lifted).unwrap();
            let lu = superposed.apply_superposed(&mut u, 0.002, pedal, lifted).unwrap();
            for (k, ((x, y), z)) in a.iter().zip(&r).zip(&u).enumerate() {
                assert!((x-y).abs() <= 1e-14, "pedal {pedal}, lifted {lifted_course:?}: coordinate {k} {x} vs {y}");
                assert!((x-z).abs() <= 1e-14, "superposed, pedal {pedal}, lifted {lifted_course:?}: coordinate {k} {x} vs {z}");
            }
            assert!((la-lr).abs() <= 1e-12*la.abs().max(1e-30), "loss {la} vs {lr}");
            assert!((la-lu).abs() <= 1e-12*la.abs().max(1e-30), "superposed loss {la} vs {lu}");
            if pedal < 1.0 { assert!(la > 0.0 && a != v0, "the dampers did nothing"); }
        }
        // the full chain is still required for an unprepared set
        let mut v = v0.clone();
        let mut plain = Prepared::new(&spec, &courses, &b).unwrap();
        assert!(plain.apply_reduced(&mut v, 0.002, 0.0, |_| false).is_err());
        assert_eq!(v, v0);
    }
    #[test]
    fn symmetric_shared_port_flow_converges_to_the_full_damping_operator() {
        // Deliberately noncommuting rows, with eigenanalysis ONLY as test oracle.
        let d = Prepared { pads: vec![StringPad { course: 0, modes: 0..2, bridge: vec![],
            drag_ns_m: 1.0, points: vec![Point { shape: vec![1.0,0.0], lift: 0.0 },
                Point { shape: vec![1.0,1.0], lift: 0.0 }] }], dimension: 2, board_start: 2, cells: 2, reduced: None };
        let eig = fs_modal::eigh_gen_dense(&[1.0,0.5,0.5,0.5], &[1.0,0.0,0.0,1.0], 2).unwrap();
        let initial = [0.3,-0.2]; let time = 0.8; let mut exact = [0.0;2];
        for e in eig {
            let a = e.phi.iter().zip(initial).map(|(g,v)| g*v).sum::<f64>()*det::exp(-e.lambda*time);
            for (v,g) in exact.iter_mut().zip(e.phi) { *v += a*g; }
        }
        let mut previous = f64::INFINITY;
        for steps in [4,8,16] {
            let mut v = initial; let mut loss = 0.0;
            for _ in 0..steps { loss += d.apply(&mut v, time/f64::from(steps), 0.0, |_| false).unwrap(); }
            let error = v.iter().zip(exact).map(|(a,b)| (a-b).powi(2)).sum::<f64>().sqrt();
            assert!(error < previous/3.5); previous = error;
            assert!((loss+0.5*v.iter().map(|x|x*x).sum::<f64>()-0.065).abs() < 1e-13);
        }
    }
}
