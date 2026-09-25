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

/// Frozen displacement/force projections in the Bank's loaded coordinates.
/// Does not contain, duplicate, reset or truncate any resonator state.
pub struct Prepared {
    pads: Vec<StringPad>,
    dimension: usize,
    board_start: usize,
    cells: usize,
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
        Ok(Self { pads, dimension, board_start, cells })
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
    fn symmetric_shared_port_flow_converges_to_the_full_damping_operator() {
        // Deliberately noncommuting rows, with eigenanalysis ONLY as test oracle.
        let d = Prepared { pads: vec![StringPad { course: 0, modes: 0..2, bridge: vec![],
            drag_ns_m: 1.0, points: vec![Point { shape: vec![1.0,0.0], lift: 0.0 },
                Point { shape: vec![1.0,1.0], lift: 0.0 }] }], dimension: 2, board_start: 2, cells: 2 };
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
