//! Supplied, moving finite-gap felt pads on the original string/board bank.
//!
//! A station observes y = g.q along its supplied transverse normal. Its
//! overlap is delta = y - height, and its force on the bank is -g*F. The
//! prescribed surface therefore does actuator work -F*delta(height). This
//! owner reuses the existing WoolFelt quadrature, Prony realization and block
//! contact solver. It owns no mechanical integrator, damper mass or gravity.
//! All trial solves keep surface travel and material history frozen.

use super::super::hammer_materials;
use super::{Bank, Course, felt, relaxation};
use felt::block::{MAX_SITES, Site as FeltSite, Workspace};
use fs_material::{Uniaxial, WoolFelt};
use fs_math::det;
use std::{collections::BTreeMap, io::Read, ops::Range};

pub const HEADER: &str = "frankensim-piano-felt-dampers-v1";
const MAX_BYTES: u64 = 1024 * 1024;
const MAX_CONTACTS: usize = 88 * MAX_SITES;
const MAX_PROJECTION_TERMS: usize = 4_000_000;
const FORCE_ABS: f64 = 1e-5;
const FORCE_REL: f64 = 1e-8;

#[derive(Clone, Copy, Debug)]
struct Pad {
    start_m: f64,
    end_m: f64,
    area_m2: f64,
    thickness_m: f64,
    closed_gap_m: f64,
    lift_m: f64,
    speed_m_s: f64,
    normal: [f64; 2],
    stations: usize,
}

/// Complete scale coverage; felt/branch records belong only to padded keys.
/// Area is the total on EACH physical speaking string, divided among stations.
#[derive(Debug)]
pub struct Specification {
    pads: BTreeMap<u8, Option<Pad>>,
    materials: BTreeMap<u8, hammer_materials::Material>,
}

impl Specification {
    pub fn read(text: &str, courses: &[Course]) -> Result<Self, String> {
        if text.len() as u64 > MAX_BYTES {
            return Err("felt damper file exceeds 1 MiB".into());
        }
        let mut header = false;
        let mut pads = BTreeMap::new();
        let mut materials = format!("{}\n", hammer_materials::HEADER);
        let mut material_rows = 0;
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() {
                continue;
            }
            if !header {
                if row != HEADER {
                    return Err(format!("line {}: expected {HEADER}", line + 1));
                }
                header = true;
                continue;
            }
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            if matches!(fields[0], "felt" | "branch") {
                materials.push_str(row);
                materials.push('\n');
                material_rows += 1;
                continue;
            }
            let parse = (|| -> Result<(u8, Option<Pad>), String> {
                if !matches!((fields[0], fields.len()), ("free", 2) | ("pad", 12)) {
                    return Err("expected free,key or pad,key,start_m,end_m,area_m2,thickness_m,closed_gap_m,lift_m,speed_m_s,normal_u,normal_v,stations".into());
                }
                let key: u8 = fields[1].parse().map_err(|_| "invalid felt damper key")?;
                if !courses.iter().any(|c| c.midi == key) {
                    return Err("felt damper key absent from scale".into());
                }
                if fields[0] == "free" {
                    return Ok((key, None));
                }
                let scalar = |i: usize| {
                    fields[i]
                        .parse::<f64>()
                        .map_err(|_| "invalid felt damper scalar".to_owned())
                };
                Ok((
                    key,
                    Some(Pad {
                        start_m: scalar(2)?,
                        end_m: scalar(3)?,
                        area_m2: scalar(4)?,
                        thickness_m: scalar(5)?,
                        closed_gap_m: scalar(6)?,
                        lift_m: scalar(7)?,
                        speed_m_s: scalar(8)?,
                        normal: [scalar(9)?, scalar(10)?],
                        stations: fields[11]
                            .parse()
                            .map_err(|_| "invalid felt damper station count")?,
                    }),
                ))
            })();
            let (key, pad) = parse.map_err(|e| format!("line {}: {e}", line + 1))?;
            if pads.insert(key, pad).is_some() {
                return Err(format!("duplicate felt damper key {key}"));
            }
        }
        if !header {
            return Err("missing felt damper header".into());
        }
        let keys: Vec<_> = courses
            .iter()
            .filter(|c| pads.get(&c.midi).is_some_and(Option::is_some))
            .map(|c| c.midi)
            .collect();
        let cards = if keys.is_empty() {
            if material_rows != 0 {
                return Err("free courses cannot carry felt damper materials".into());
            }
            Vec::new()
        } else {
            hammer_materials::read(&materials, &keys)
                .map_err(|e| format!("damper material: {e}"))?
        };
        let spec = Self {
            pads,
            materials: keys.into_iter().zip(cards).collect(),
        };
        spec.validate(courses)?;
        Ok(spec)
    }

    pub fn load(path: &str, courses: &[Course]) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path)
            .map_err(|e| format!("{path}: {e}"))?
            .take(MAX_BYTES + 1)
            .read_to_string(&mut text)
            .map_err(|e| format!("{path}: {e}"))?;
        Self::read(&text, courses)
    }

    fn validate(&self, courses: &[Course]) -> Result<(), String> {
        if courses.is_empty()
            || courses.len() > 88
            || self.pads.len() != courses.len()
            || courses
                .iter()
                .enumerate()
                .any(|(i, c)| courses[..i].iter().any(|p| p.midi == c.midi))
        {
            return Err("felt damper rows must cover a nonempty unique scale exactly".into());
        }
        for c in courses {
            c.validate()?;
            let pad = self
                .pads
                .get(&c.midi)
                .ok_or_else(|| format!("missing felt damper key {}", c.midi))?;
            let Some(p) = pad else {
                continue;
            };
            if [
                p.start_m,
                p.end_m,
                p.area_m2,
                p.thickness_m,
                p.closed_gap_m,
                p.lift_m,
                p.speed_m_s,
                p.normal[0],
                p.normal[1],
            ]
            .iter()
            .any(|x| !x.is_finite())
                || p.start_m <= 0.
                || p.end_m <= p.start_m
                || p.end_m >= c.length_m
                || p.area_m2 <= 0.
                || p.thickness_m <= 0.
                || p.closed_gap_m < 0.
                || p.lift_m <= 0.
                || p.speed_m_s <= 0.
                || !(1..=4).contains(&p.stations)
                || (p.normal[0].hypot(p.normal[1]) - 1.).abs() > 1e-10
                || !(p.closed_gap_m + p.lift_m).is_finite()
                || !(p.area_m2 * p.thickness_m).is_finite()
                || p.area_m2 * p.thickness_m <= 0.
            {
                return Err(format!(
                    "key {}: invalid finite pad geometry, positive area/thickness/lift/speed, nonnegative closed gap, unit transverse normal or 1..4 stations",
                    c.midi
                ));
            }
            if !self.materials.contains_key(&c.midi) {
                return Err(format!("missing felt damper material for key {}", c.midi));
            }
        }
        Ok(())
    }
}

struct Projection {
    indices: Vec<usize>,
    values: Vec<f64>,
}
impl Projection {
    fn position(&self, q: &[f64]) -> f64 {
        self.indices
            .iter()
            .zip(&self.values)
            .map(|(&i, &g)| g * q[i])
            .sum()
    }
    fn gather(&self, force: f64, out: &mut [f64]) {
        for (&i, &g) in self.indices.iter().zip(&self.values) {
            out[i] += g * force;
        }
    }
}

struct Site {
    projection: Projection,
    course: usize,
    member: usize,
    station: usize,
}

#[derive(Clone)]
struct History {
    felt: felt::State,
    memory: relaxation::Memory,
    overlap: f64,
    force: f64,
}

struct Group {
    course: usize,
    pad: Pad,
    law: WoolFelt,
    material: relaxation::Prepared,
    sites: Range<usize>,
    compliance: Vec<f64>,
    height: f64,
    next_height: f64,
    saved_height: f64,
}

/// An accepted station snapshot; numerical force iterates are not published.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactObservation {
    pub course: usize,
    pub member: usize,
    pub station: usize,
    pub overlap_m: f64,
    pub force_n: f64,
    pub maximum_felt_strain: f64,
    pub maximum_felt_stress_pa: f64,
    pub creep_coordinates_m: [f64; relaxation::MAX_TERMS],
}

/// Cold projections and reciprocal per-course compliance, with independent
/// station history. All hot storage is prepared here, including rollback.
pub struct Prepared {
    groups: Vec<Group>,
    sites: Vec<Site>,
    history: Vec<History>,
    next_history: Vec<History>,
    saved_history: Vec<History>,
    force: Vec<f64>,
    corrected_force: Vec<f64>,
    generalized: Vec<f64>,
    workspace: Workspace,
    dt: f64,
    physical_pads: usize,
    in_step: bool,
}

impl Prepared {
    pub fn new(spec: &Specification, courses: &[Course], bank: &Bank) -> Result<Self, String> {
        spec.validate(courses)?;
        let board_start = bank.modes.len();
        let dimension = board_start + bank.board_count;
        if bank.rate == 0
            || bank.q.len() != dimension
            || bank.v.len() != dimension
            || bank.q.iter().chain(&bank.v).any(|x| *x != 0.)
        {
            return Err("felt dampers require the original, unexcited bank coordinates".into());
        }
        let dt = 1. / f64::from(bank.rate);
        let mut groups = Vec::new();
        let mut sites = Vec::new();
        let mut history = Vec::new();
        let mut terms = 0usize;
        let mut physical_pads = 0;
        for (ci, c) in courses.iter().enumerate() {
            let Some(pad) = spec.pads[&c.midi] else {
                continue;
            };
            if pad.normal[1] != 0. && !bank.has_secondary_polarization() {
                return Err("a lateral felt damper normal needs the supplied second transverse polarization".into());
            }
            if !(pad.speed_m_s * dt).is_finite()
                || pad.speed_m_s * dt == 0.
                || pad.closed_gap_m + pad.lift_m <= pad.closed_gap_m
            {
                return Err(
                    "felt damper travel is not representable at this mechanics rate".into(),
                );
            }
            let (law, prony) = &spec.materials[&c.midi];
            let area = pad.area_m2 / pad.stations as f64;
            let material =
                relaxation::Spectrum::from_prony(prony)?.prepare(area, pad.thickness_m, dt)?;
            let first = sites.len();
            for member in 0..c.unison {
                let primary: Vec<_> = bank
                    .strings
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| {
                        s.course == ci && s.member == member && !s.duplex && s.polarization == 0
                    })
                    .collect();
                if primary.len() != 1 {
                    return Err("felt damper speaking-string identity is ambiguous".into());
                }
                let secondary: Vec<_> = bank
                    .strings
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| {
                        s.course == ci && s.member == member && !s.duplex && s.polarization == 1
                    })
                    .collect();
                if pad.normal[1] != 0. && secondary.len() != 1 {
                    return Err(
                        "felt damper needs both directions of one physical speaking segment".into(),
                    );
                }
                let root_mass = det::sqrt(c.modal_mass_kg());
                for station in 0..pad.stations {
                    let fraction = (pad.start_m
                        + (pad.end_m - pad.start_m) * (station as f64 + 0.5) / pad.stations as f64)
                        / c.length_m;
                    if !fraction.is_finite()
                        || fraction <= pad.start_m / c.length_m
                        || fraction >= pad.end_m / c.length_m
                    {
                        return Err(
                            "felt damper station has no representable span resolution".into()
                        );
                    }
                    let mut indices = Vec::new();
                    let mut values = Vec::new();
                    let mut bridge = vec![0.; bank.board_count];
                    for plane in 0..2 {
                        let normal = pad.normal[plane];
                        if normal == 0. {
                            continue;
                        }
                        let (si, s) = if plane == 0 { primary[0] } else { secondary[0] };
                        if s.modes.is_empty()
                            || s.modes.end > board_start
                            || s.bridge.len() != bank.board_count
                            || bank.modes[s.modes.clone()].iter().any(|m| m.string != si)
                        {
                            return Err(
                                "felt damper projection does not match the string/bridge basis"
                                    .into(),
                            );
                        }
                        let mut lift = fraction;
                        for (mode, k) in s.modes.clone().enumerate() {
                            let shape =
                                det::sin((mode + 1) as f64 * std::f64::consts::PI * fraction)
                                    / root_mass;
                            lift -= shape * bank.modes[k].beta;
                            let value = normal * shape;
                            if value != 0. {
                                indices.push(k);
                                values.push(value);
                            }
                        }
                        for (out, b) in bridge.iter_mut().zip(&s.bridge) {
                            *out += normal * lift * b;
                        }
                    }
                    for (j, value) in bridge.into_iter().enumerate() {
                        if value != 0. {
                            indices.push(board_start + j);
                            values.push(value);
                        }
                    }
                    let norm: f64 = values.iter().map(|x| x * x).sum();
                    if values.iter().any(|x| !x.is_finite()) || !norm.is_finite() || norm <= 0. {
                        return Err("felt damper projection is nonfinite or unobservable".into());
                    }
                    terms = terms
                        .checked_add(values.len())
                        .ok_or("felt damper projection size overflow")?;
                    if terms > MAX_PROJECTION_TERMS || sites.len() >= MAX_CONTACTS {
                        return Err(
                            "felt dampers exceed 1056 stations or four million projection terms"
                                .into(),
                        );
                    }
                    sites.push(Site {
                        projection: Projection { indices, values },
                        course: ci,
                        member,
                        station,
                    });
                    history.push(History {
                        felt: law.initial_state(),
                        memory: relaxation::Memory::default(),
                        overlap: -pad.closed_gap_m,
                        force: 0.,
                    });
                }
                physical_pads += 1;
            }
            let n = sites.len() - first;
            if n == 0 || n > MAX_SITES {
                return Err("felt damper course exceeds twelve contact stations".into());
            }
            groups.push(Group {
                course: ci,
                pad,
                law: law.clone(),
                material,
                sites: first..sites.len(),
                compliance: vec![0.; n * n],
                height: pad.closed_gap_m,
                next_height: pad.closed_gap_m,
                saved_height: pad.closed_gap_m,
            });
        }
        // Only two dense cold vectors; each station retains its sparse row.
        // The response is the bank's actual coupled discrete displacement map,
        // not an isolated-string or static compliance approximation.
        let mut applied = vec![0.; dimension];
        let mut response = vec![0.; dimension];
        for group in &mut groups {
            let n = group.sites.len();
            for col in 0..n {
                applied.fill(0.);
                sites[group.sites.start + col]
                    .projection
                    .gather(1., &mut applied);
                bank.force_response(&applied, &mut response)?;
                for row in 0..n {
                    group.compliance[row * n + col] = sites[group.sites.start + row]
                        .projection
                        .position(&response);
                }
            }
            if group.compliance.iter().any(|x| !x.is_finite())
                || (0..n).any(|i| group.compliance[i * n + i] <= 0.)
            {
                return Err(
                    "felt damper response has no positive finite contact compliance".into(),
                );
            }
        }
        let count = sites.len();
        Ok(Self {
            groups,
            sites,
            next_history: history.clone(),
            saved_history: history.clone(),
            history,
            force: vec![0.; count],
            corrected_force: vec![0.; count],
            generalized: vec![0.; dimension],
            workspace: Workspace::new()?,
            dt,
            physical_pads,
            in_step: false,
        })
    }

    pub fn resolution(&self) -> (usize, usize) {
        (self.physical_pads, self.sites.len())
    }
    pub fn pad_height_m(&self, course: usize) -> Option<f64> {
        self.groups
            .iter()
            .find(|g| g.course == course)
            .map(|g| g.height)
    }
    pub fn contact_state(&self, index: usize) -> Option<ContactObservation> {
        let site = self.sites.get(index)?;
        let h = &self.history[index];
        Some(ContactObservation {
            course: site.course,
            member: site.member,
            station: site.station,
            overlap_m: h.overlap,
            force_n: h.force,
            maximum_felt_strain: h.felt.eps_max,
            maximum_felt_stress_pa: h.felt.sig_max,
            creep_coordinates_m: h.memory.0,
        })
    }

    /// Predict prescribed travel once, before numerical force iteration. A held
    /// key or sostenuto latch raises the pad fully; sustain sets fractional lift.
    pub fn begin_step(
        &mut self,
        sustain: f64,
        lifted: impl Fn(usize) -> bool,
    ) -> Result<(), &'static str> {
        if self.in_step || !sustain.is_finite() || !(0. ..=1.).contains(&sustain) {
            return Err("invalid felt damper step or pedal control");
        }
        for group in &mut self.groups {
            let travel = if lifted(group.course) { 1. } else { sustain };
            let target = group.pad.closed_gap_m + travel * group.pad.lift_m;
            let change = target - group.height;
            let limit = group.pad.speed_m_s * self.dt;
            group.next_height = if change.abs() <= limit {
                target
            } else {
                group.height + change.signum() * limit
            };
            if !group.next_height.is_finite() {
                return Err("felt damper travel overflow");
            }
        }
        for (force, h) in self.force.iter_mut().zip(&self.history) {
            *force = h.force;
        }
        self.regather()?;
        self.in_step = true;
        Ok(())
    }

    pub fn generalized_force(&self) -> &[f64] {
        &self.generalized
    }

    fn regather(&mut self) -> Result<(), &'static str> {
        self.generalized.fill(0.);
        for (site, &force) in self.sites.iter().zip(&self.force) {
            site.projection.gather(-force, &mut self.generalized);
        }
        if self.generalized.iter().any(|x| !x.is_finite()) {
            return Err("felt damper generalized force overflow");
        }
        Ok(())
    }

    fn residual_passes(&self, candidate_q: &[f64]) -> Result<bool, &'static str> {
        if !self.in_step
            || candidate_q.len() != self.generalized.len()
            || candidate_q.iter().any(|x| !x.is_finite())
        {
            return Err("invalid felt damper candidate state");
        }
        let mut converged = true;
        for group in &self.groups {
            for i in group.sites.clone() {
                let h = &self.history[i];
                let total = self.sites[i].projection.position(candidate_q) - group.next_height;
                let end = total
                    - group.material.free_deformation(&h.memory)
                    - group.material.compliance() * self.force[i];
                let start = h.overlap - group.material.deformation(&h.memory);
                let expected = felt::average(
                    &group.law,
                    &h.felt,
                    start,
                    end,
                    group.pad.thickness_m,
                    group.pad.area_m2 / group.pad.stations as f64,
                )
                .0;
                if !total.is_finite()
                    || !end.is_finite()
                    || !expected.is_finite()
                    || expected < 0.
                    || !self.force[i].is_finite()
                    || self.force[i] < 0.
                {
                    return Err("nonfinite or tensile felt damper contact force");
                }
                if (expected == 0. && self.force[i] != 0.)
                    || (self.force[i] - expected).abs() > FORCE_ABS + FORCE_REL * expected.abs()
                {
                    converged = false;
                }
            }
        }
        Ok(converged)
    }

    /// True only when the SAME held forces that produced candidate_q satisfy
    /// their constitutive equations. A corrected force always requires another
    /// bank prediction; material state and travel never advance in this method.
    pub fn correct(&mut self, candidate_q: &[f64]) -> Result<bool, &'static str> {
        if self.residual_passes(candidate_q)? {
            return Ok(true);
        }
        self.corrected_force.copy_from_slice(&self.force);
        for group in &self.groups {
            let n = group.sites.len();
            let first = group.sites.start;
            let area = group.pad.area_m2 / group.pad.stations as f64;
            let mut sites = [FeltSite {
                law: &group.law,
                history: &self.history[first].felt,
                start_m: 0.,
                thickness_m: group.pad.thickness_m,
                area_m2: area,
            }; MAX_SITES];
            let mut a = [0.; MAX_SITES * MAX_SITES];
            let mut free = [0.; MAX_SITES];
            let mut warm = [0.; MAX_SITES];
            let mut output = [0.; MAX_SITES];
            for row in 0..n {
                let i = first + row;
                let h = &self.history[i];
                sites[row].history = &h.felt;
                sites[row].start_m = h.overlap - group.material.deformation(&h.memory);
                free[row] = self.sites[i].projection.position(candidate_q)
                    - group.next_height
                    - group.material.free_deformation(&h.memory);
                warm[row] = self.force[i];
                for col in 0..n {
                    let compliance = group.compliance[row * n + col];
                    a[row * n + col] = compliance;
                    free[row] += compliance * self.force[first + col];
                }
                a[row * n + row] += group.material.compliance();
            }
            self.workspace.solve(
                &sites[..n],
                &a[..n * n],
                &free[..n],
                &warm[..n],
                &mut output[..n],
            )?;
            for row in 0..n {
                // An exactly separated constitutive branch carries no force.
                // Do not leave a sub-tolerance warm force doing spurious work
                // indefinitely after lift; re-predict with exact zero instead.
                self.corrected_force[first + row] = if output[row] == 0. {
                    0.
                } else {
                    0.5 * self.force[first + row] + 0.5 * output[row]
                };
            }
        }
        self.force.copy_from_slice(&self.corrected_force);
        self.regather()?;
        Ok(false)
    }

    /// Commit one accepted mechanics tick. Returns total damper loss, its Prony
    /// subset, and signed actuator work. The subset is not counted twice.
    pub fn accept(&mut self, q: &[f64]) -> Result<(f64, f64, f64), &'static str> {
        if !self.residual_passes(q)? {
            return Err("felt damper force was not converged before acceptance");
        }
        let mut loss = 0.;
        let mut relaxation_loss = 0.;
        let mut actuator_work = 0.;
        for group in &self.groups {
            let thickness = group.pad.thickness_m;
            let volume = group.pad.area_m2 / group.pad.stations as f64 * thickness;
            for i in group.sites.clone() {
                let old = &self.history[i];
                let force = self.force[i];
                let overlap = self.sites[i].projection.position(q) - group.next_height;
                let (memory, viscous) = group.material.advance(&old.memory, force);
                let start = old.overlap - group.material.deformation(&old.memory);
                let end = overlap - group.material.deformation(&memory);
                if !overlap.is_finite()
                    || !end.is_finite()
                    || overlap / thickness > group.law.eps_densify + 1e-10
                {
                    return Err("felt damper exceeded total densification strain");
                }
                let state = group.law.update_state(end / thickness, &old.felt);
                let stored_start = volume * felt::stored(&group.law, start / thickness, &old.felt);
                let stored_end = volume * felt::stored(&group.law, end / thickness, &state);
                let crush = force * (end - start) - (stored_end - stored_start);
                let tolerance = 1e-12
                    + 1e-8 * stored_start.abs().max(stored_end.abs())
                    + (FORCE_ABS + FORCE_REL * force.abs()) * (end - start).abs();
                if !crush.is_finite()
                    || crush < -tolerance
                    || !viscous.is_finite()
                    || viscous < 0.
                    || !stored_end.is_finite()
                    || memory.0.iter().any(|x| !x.is_finite())
                {
                    return Err("felt damper failed its material work balance");
                }
                loss += crush + viscous;
                relaxation_loss += viscous;
                actuator_work -= force * (group.next_height - group.height);
                self.next_history[i] = History {
                    felt: state,
                    memory,
                    overlap,
                    force,
                };
            }
        }
        if !loss.is_finite() || !relaxation_loss.is_finite() || !actuator_work.is_finite() {
            return Err("felt damper work accounting overflow");
        }
        self.history.clone_from_slice(&self.next_history);
        for group in &mut self.groups {
            group.height = group.next_height;
        }
        self.in_step = false;
        Ok((loss, relaxation_loss, actuator_work))
    }

    pub fn energy(&self) -> f64 {
        let mut energy = 0.;
        for group in &self.groups {
            let thickness = group.pad.thickness_m;
            let volume = group.pad.area_m2 / group.pad.stations as f64 * thickness;
            for i in group.sites.clone() {
                let h = &self.history[i];
                let elastic = h.overlap - group.material.deformation(&h.memory);
                energy += group.material.stored(&h.memory)
                    + volume * felt::stored(&group.law, elastic / thickness, &h.felt);
            }
        }
        energy
    }

    pub fn checkpoint(&mut self) {
        self.saved_history.clone_from_slice(&self.history);
        for group in &mut self.groups {
            group.saved_height = group.height;
        }
    }
    pub fn restore(&mut self) {
        self.history.clone_from_slice(&self.saved_history);
        for group in &mut self.groups {
            group.height = group.saved_height;
            group.next_height = group.height;
        }
        for (force, h) in self.force.iter_mut().zip(&self.history) {
            *force = h.force;
        }
        // Saved accepted forces/projections were already finite. Restoring
        // them never depends on an abandoned Newton iterate or target travel.
        self.generalized.fill(0.);
        for (site, &force) in self.sites.iter().zip(&self.force) {
            site.projection.gather(-force, &mut self.generalized);
        }
        self.in_step = false;
    }
}

#[cfg(test)]
#[path = "felt_damper_tests.rs"]
mod tests;
