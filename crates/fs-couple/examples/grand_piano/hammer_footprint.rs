//! Finite longitudinal hammer faces in the unchanged moving-bridge basis.
//!
//! Each quadrature site is a separate work-conjugate contact, not an averaged
//! displacement passed to one nonlinear felt law. This module owns only cold
//! geometry/projection. The engine owns each site's existing felt/Prony history.
use std::{collections::BTreeMap, io::Read, ops::Range};
use fs_math::det;
use super::Bank;
use super::super::geometry::Course;

const MAX_BYTES: u64 = 64 * 1024;
/// Four sites on each of three strings on every one of the 88 keys.
const MAX_CONTACTS: usize = 88 * 3 * 4;
const MAX_SHAPE_TERMS: usize = MAX_CONTACTS * super::MAX_STRING_MODES;
pub const HEADER: &str = "frankensim-hammer-footprints-v1";

/// Uniform longitudinal face centred at the scale's existing strike station.
/// Total area and thickness still belong to that course's physical felt card.
#[derive(Clone, Copy, Debug)]
pub struct Footprint { pub length_m: f64, pub sites: usize }

#[derive(Clone, Debug)]
pub struct Specification { pub footprints: BTreeMap<u8, Option<Footprint>> }
impl Specification {
    /// Complete per-key selection: `point,key` or `span,key,length_m,sites`.
    /// Span sites are 2 or 4 Gauss points. No default width/material is inferred.
    pub fn read(text: &str, courses: &[Course]) -> Result<Self, String> {
        if text.len() as u64 > MAX_BYTES { return Err("hammer footprint file exceeds 64 KiB".into()); }
        let mut footprints = BTreeMap::new(); let mut header = false;
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            let error = || format!("hammer footprint line {}: expected point,key or span,key,length_m,sites", line + 1);
            if !header {
                if row != HEADER { return Err(error()); }
                header = true; continue;
            }
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            if fields.len() != 2 && fields.len() != 4 { return Err(error()); }
            let key = fields[1].parse::<u8>().map_err(|_| error())?;
            let footprint = match (fields[0], fields.len()) {
                ("point", 2) => None,
                ("span", 4) => Some(Footprint {
                    length_m: fields[2].parse().map_err(|_| error())?,
                    sites: fields[3].parse().map_err(|_| error())?,
                }),
                _ => return Err(error()),
            };
            if footprints.insert(key, footprint).is_some() { return Err(format!("duplicate hammer footprint key {key}")); }
        }
        if !header { return Err("missing hammer footprint header".into()); }
        let result = Self { footprints }; result.validate(courses)?; Ok(result)
    }

    pub fn load(path: &str, courses: &[Course]) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?
            .take(MAX_BYTES + 1).read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
        Self::read(&text, courses)
    }

    fn validate(&self, courses: &[Course]) -> Result<(), String> {
        if courses.is_empty() || courses.len() > 88 || self.footprints.len() != courses.len()
            || courses.iter().enumerate().any(|(i,c)| courses[..i].iter().any(|p| p.midi == c.midi)) {
            return Err("hammer footprints must cover a nonempty unique scale exactly".into());
        }
        for c in courses {
            c.validate()?;
            let entry = self.footprints.get(&c.midi).ok_or_else(|| format!("missing hammer footprint key {}", c.midi))?;
            if let Some(p) = entry {
                let centre = c.strike_fraction * c.length_m;
                let left = centre - 0.5 * p.length_m;
                let right = centre + 0.5 * p.length_m;
                if !p.length_m.is_finite() || p.length_m <= 0.0 || !matches!(p.sites, 2 | 4)
                    || !left.is_finite() || !right.is_finite() || left <= 0.0
                    || left >= centre || right <= centre || right >= c.length_m {
                    return Err(format!("key {}: complete hammer span must lie inside the speaking string; sites must be 2 or 4", c.midi));
                }
            }
        }
        Ok(())
    }
}

/// Abscissas on [-1,1] and positive weights normalized to total face area.
fn quadrature(sites: usize) -> Vec<(f64, f64)> {
    if sites == 2 {
        let x = 1.0 / det::sqrt(3.0); vec![(-x, 0.5), (x, 0.5)]
    } else {
        let inner = det::sqrt((3.0 - 2.0 * det::sqrt(1.2)) / 7.0);
        let outer = det::sqrt((3.0 + 2.0 * det::sqrt(1.2)) / 7.0);
        let small = (18.0 - det::sqrt(30.0)) / 72.0;
        let large = (18.0 + det::sqrt(30.0)) / 72.0;
        vec![(-outer, small), (-inner, large), (inner, large), (outer, small)]
    }
}

pub(super) struct Point { pub shape: Vec<f64>, pub lift: f64, pub fraction: f64 }
pub(super) struct Prepared {
    pub points: Vec<Point>,
    ranges: Vec<Range<usize>>,
    pub modal_force: Vec<f64>,
    pub bridge_force: Vec<f64>,
}
impl Prepared {
    /// Project all site forces BEFORE the unchanged modal propagation. A point
    /// force is in newtons, already area-scaled by its constitutive owner.
    pub fn gather(&mut self, strings: &[super::StringPort], force: &[f64]) {
        self.modal_force.fill(0.0); self.bridge_force.fill(0.0);
        for (si, s) in strings.iter().enumerate() {
            for c in self.ranges[si].clone() {
                if force[c] == 0.0 { continue; }
                let p = &self.points[c];
                self.bridge_force[si] += p.lift * force[c];
                for (k,g) in s.modes.clone().zip(&p.shape) { self.modal_force[k] += g * force[c]; }
            }
        }
    }
}

impl Bank {
    /// Prepare a new bank with finite contact geometry. The same courses feed
    /// BOTH mass/stiffness assembly and contact projection, so no mismatched
    /// string scale can be attached to a running bank.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_hammer_footprints(courses: &[Course], board: &[super::BoardMode],
        rate: u32, band_hz: f64, max_modes: usize, damping: bool,
        spec: &Specification) -> Result<Self, String> {
        Self::new_with_hammer_footprints_and_transverse_bridge(courses,board,rate,band_hz,
            max_modes,damping,spec,None)
    }

    /// Project the same primary-plane face after completing both directional
    /// string mass forms. The extra polarization has no duplicate felt sites.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_hammer_footprints_and_transverse_bridge(courses: &[Course], board: &[super::BoardMode],
        rate: u32, band_hz: f64, max_modes: usize, damping: bool,
        spec: &Specification, secondary: Option<&[Vec<f64>]>) -> Result<Self, String> {
        spec.validate(courses)?;
        let mut bank = Self::new_with_transverse_bridge(courses,board,rate,band_hz,max_modes,damping,secondary)?;
        bank.configure_hammer_footprints(spec, courses)?;
        Ok(bank)
    }

    /// Cold contact preparation on a fresh bank. Oscillators, mass loading,
    /// stiffness, damping, board basis and physical motion are NOT replaced.
    /// Refuse a moving or previously configured bank rather than discarding
    /// history. All-point input retains the original image bit for bit.
    fn configure_hammer_footprints(&mut self, spec: &Specification, courses: &[Course]) -> Result<(), String> {
        spec.validate(courses)?;
        if self.footprints.is_some() || self.q.iter().chain(&self.v).any(|x| *x != 0.0)
            || self.strings.iter().any(|s| s.course >= courses.len()) {
            return Err("hammer footprints require a fresh unconfigured bank and its original course order".into());
        }
        if spec.footprints.values().all(Option::is_none) { return Ok(()); }
        let mut points = Vec::new(); let mut contacts = Vec::new();
        let mut ranges = Vec::with_capacity(self.strings.len()); let mut terms = 0usize;
        for (si, s) in self.strings.iter().enumerate() {
            let first = points.len();
            if s.contact.is_some() {
                let c = &courses[s.course];
                match spec.footprints[&c.midi] {
                    None => {
                        points.push(Point { shape: s.modes.clone().map(|k| self.modes[k].hammer_shape).collect(),
                            lift: s.hammer_lift, fraction: 1.0 });
                        contacts.push(si);
                    }
                    Some(p) => {
                        let half = 0.5 * p.length_m / c.length_m;
                        let mut previous = 0.0;
                        for (x, fraction) in quadrature(p.sites) {
                            let station = c.strike_fraction + half * x;
                            if !station.is_finite() || station <= previous || station >= 1.0 {
                                return Err("hammer footprint quadrature has no representable separation".into());
                            }
                            previous = station;
                            let shape: Vec<f64> = (1..=s.modes.len()).map(|n|
                                det::sin(n as f64 * std::f64::consts::PI * station) / det::sqrt(c.modal_mass_kg())).collect();
                            let lift = station - shape.iter().zip(&self.modes[s.modes.clone()]).map(|(g,m)|g*m.beta).sum::<f64>();
                            if !lift.is_finite() || shape.iter().any(|x| !x.is_finite()) {
                                return Err("hammer footprint projection overflow".into());
                            }
                            points.push(Point { shape, lift, fraction }); contacts.push(si);
                        }
                    }
                }
            }
            let end = points.len();
            terms = terms.checked_add((end-first) * s.modes.len()).ok_or("hammer footprint extent overflow")?;
            if end > MAX_CONTACTS || terms > MAX_SHAPE_TERMS { return Err("hammer footprint contact/shape budget exceeded".into()); }
            ranges.push(first..end);
        }
        let nc = points.len(); let r = self.board_count;
        let mut w = vec![0.0; nc*r]; let mut compliance = vec![0.0; nc*nc];
        for (si,s) in self.strings.iter().enumerate() {
            for i in ranges[si].clone() {
                let p = &points[i];
                let lift = p.lift + s.modes.clone().zip(&p.shape).map(|(k,g)|
                    0.5*self.modes[k].a*self.transition[k].bq*g).sum::<f64>();
                for j in 0..r { w[i*r+j] = lift*s.bridge[j]; }
                // Different stations on the SAME string share its direct
                // displacement compliance as well as its reciprocal board.
                for j in ranges[si].clone().filter(|j| *j >= i) {
                    let value: f64 = s.modes.clone().zip(p.shape.iter().zip(&points[j].shape))
                        .map(|(k,(a,b))| self.transition[k].bq*a*b).sum();
                    compliance[i*nc+j] = value; compliance[j*nc+i] = value;
                }
            }
        }
        let mut response = vec![0.0; nc*r];
        for c in 0..nc { for a in 0..r {
            response[c*r+a] = (0..r).map(|b| self.schur_inverse[a*r+b]*w[c*r+b]).sum();
        } }
        for i in 0..nc { for j in i..nc {
            compliance[i*nc+j] += (0..r).map(|a| w[i*r+a]*response[j*r+a]).sum::<f64>();
            compliance[j*nc+i] = compliance[i*nc+j];
        } }
        if compliance.iter().chain(&response).any(|x| !x.is_finite())
            || (0..nc).any(|i| compliance[i*nc+i] <= 0.0) {
            return Err("hammer footprint compliance is not finite positive on its diagonal".into());
        }
        // Publication follows complete admission; no live state was touched.
        for (si,s) in self.strings.iter_mut().enumerate() {
            s.contact = (!ranges[si].is_empty()).then_some(ranges[si].start);
        }
        self.contact_strings = contacts; self.contact_compliance = compliance;
        self.contact_response = response; self.free_contact = vec![0.0; nc];
        #[cfg(test)] { self.contact_board = w; }
        self.footprints = Some(Prepared { points, ranges,
            modal_force: vec![0.0; self.modes.len()], bridge_force: vec![0.0; self.strings.len()] });
        Ok(())
    }

    /// Fraction of ONE string's allocated felt area at a contact site.
    /// The course's original area is divided among unison strings separately.
    pub fn contact_area_fraction(&self, c: usize) -> f64 {
        self.footprints.as_ref().map_or(1.0, |p| p.points[c].fraction)
    }
}

#[cfg(test)]
#[path = "hammer_footprint/tests.rs"]
mod tests;
