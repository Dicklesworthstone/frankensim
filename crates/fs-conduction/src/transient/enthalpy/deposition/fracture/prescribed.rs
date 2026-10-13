//! Event-resolved cooling under explicitly prescribed, fixed-frame separation.
//!
//! A straight jump segment has affine normal gap and convex effective opening
//! radius. Starting inside the failure surface, it can cross that surface
//! outward at most once; each gap threshold can also be crossed at most once.
//! Locate these three possible crossings independently, rather than detecting
//! a change only from the final thermal-contact state (which can miss an open
//! and reclose in the same segment). Queries use the actual cohesive law.
//!
//! Brackets describe numerical switching locations on the DECLARED linear jump
//! path, not a proof about an unknown structural trajectory. Split prescribed
//! paths at reversals. No equilibrium, contact discovery, displacement feedback,
//! fracture heating or simultaneous healing is supplied by this operation.

use super::*;

/// One endpoint of a straight local-jump segment, in metres. The start is the
/// accepted point with this ID; caller arrival order has no physical meaning.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Endpoint {
    pub interface: u64,
    pub jump_m: [f64; 3],
}

/// Whole-segment work bounds and independent thermal/event numerical controls.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub thermal: StepConfig,
    /// Width allowed for each numerical event-time bracket, seconds.
    pub event_time_tolerance_s: f64,
    pub local_energy_tolerance_j_m2: f64,
    /// Total predicate queries across all failure and gap searches.
    pub max_event_trials: usize,
    /// Maximum retained failure/gap crossings, including inactive gap crossings.
    pub max_crossings: usize,
    /// Shared across ALL event-split thermal intervals. Each interval also
    /// obeys the thermal simulation's own (possibly tighter) substep ceiling.
    pub max_thermal_substeps: usize,
}

/// A threshold crossing is not necessarily a change in thermal contact: an
/// intact cohesive interface remains closed even after crossing the open gap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CrossingKind {
    Failure,
    OpeningGap,
    ReclosingGap,
}

/// False/true predicate bracket on the declared path. Its upper endpoint is a
/// thermal cut. A simultaneous update at another cut may detect the switch
/// earlier inside this bracket. Switching delay is at most the bracket width
/// relative to the numerical predicate transition, not an interval-verified
/// physical accuracy bound. Time rounding can give equal bracket endpoints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    pub interface: u64,
    pub kind: CrossingKind,
    pub before_s: f64,
    pub after_s: f64,
}

/// Whole-path accounting. Only a complete successful result can be published.
#[derive(Debug, Clone)]
pub struct Report {
    pub thermal: EnergyAudit,
    pub contacts: ContactReport,
    pub crossings: Vec<Crossing>,
    pub event_trials: usize,
    pub thermal_intervals: usize,
}

struct Located {
    id: u64,
    kind: CrossingKind,
    low: f64,
    high: f64,
}

/// Canonical interpolation used by both event predicates and material updates.
/// Exact endpoints stay bitwise unchanged. A convex sum avoids end-start overflow.
fn jump_at(start: [f64; 3], end: [f64; 3], fraction: f64) -> [f64; 3] {
    if fraction == 0.0 { return start; }
    if fraction == 1.0 { return end; }
    std::array::from_fn(|i| (1.0-fraction)*start[i]+fraction*end[i])
}
fn stamp(start: f64, end: f64, fraction: f64) -> f64 {
    if fraction == 0.0 { start }
    else if fraction == 1.0 { end }
    else { (start+fraction*(end-start)).clamp(start,end) }
}

fn locate(
    start_s: f64, end_s: f64, config: Config, trials: &mut usize,
    mut predicate: impl FnMut(f64) -> Result<bool, Error>,
) -> Result<Option<(f64, f64)>, Error> {
    let mut query = |fraction| {
        if *trials >= config.max_event_trials {
            return Err(Error::Thermal(DepositionError::Budget("prescribed contact event queries")));
        }
        *trials += 1;
        predicate(fraction)
    };
    // Every caller supplies a predicate whose false set is an initial interval.
    // Already active predicates have no new crossing on this segment.
    if query(0.0)? || !query(1.0)? { return Ok(None); }
    let mut low = 0.0;
    let mut high = 1.0;
    while stamp(start_s,end_s,high)-stamp(start_s,end_s,low) > config.event_time_tolerance_s {
        let middle = low+0.5*(high-low);
        if middle == low || middle == high {
            return Err(Error::Invalid("event-time tolerance is below representable path resolution"));
        }
        if query(middle)? { high = middle; } else { low = middle; }
    }
    Ok(Some((low,high)))
}

impl<'t, 'm> FractureCooling<'t, 'm> {
    /// Evolve a complete straight prescribed separation segment and its heat.
    ///
    /// All interfaces move simultaneously from their retained jumps to supplied
    /// endpoints. Locate complete failure and directed normal opening/reclosure
    /// threshold crossings, then march the existing thermal kernel between them.
    /// Scheduled unrelated births/boundary changes keep their exact event times;
    /// scheduled changes to mechanically owned interfaces still refuse.
    ///
    /// At equal times the thermal schedule is applied first, then every bound
    /// mechanical point is updated in stable ID order. No energy is added by a
    /// closure change. Mechanical and thermal histories publish together only
    /// after the whole path and its global thermal-energy gate succeed.
    ///
    /// Unlike `synchronize`, endpoints are prescribed JUMPS, not claims that an
    /// independently solved endpoint encodes the same intervening path. Do not
    /// pass sparse equilibrium samples as an exact motion history. Material-law
    /// dependence on temperature and trajectory reconstruction remain outside
    /// this fixed-law, fixed-frame, prescribed-motion model.
    pub fn advance_prescribed_to(
        &self, target_s: f64, endpoints: &[Endpoint], config: Config,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Report), Error> {
        checkpoint()?;
        let start_s = self.thermal.time_s();
        if !target_s.is_finite() || !positive(target_s-start_s)
            || !positive(config.event_time_tolerance_s)
            || !positive(config.local_energy_tolerance_j_m2)
            || !positive(config.thermal.max_dt_s) || !nonnegative(config.thermal.energy_tolerance_j)
            || endpoints.len() != self.contacts.0.len()
        { return Err(Error::Invalid("require positive duration, valid numerical controls and complete prescribed endpoints")); }
        let mut ends = BTreeMap::new();
        for endpoint in endpoints {
            checkpoint()?;
            let old = self.contacts.0.get(&endpoint.interface)
                .ok_or(Error::Invalid("unknown prescribed interface"))?;
            if ends.insert(endpoint.interface,endpoint.jump_m).is_some() {
                return Err(Error::Invalid("duplicate prescribed interface"));
            }
            // Endpoints bound the norm and slip of a straight segment by convexity.
            // Validate the complete material path before any thermal integration.
            let end = old.point.trial(endpoint.jump_m,config.local_energy_tolerance_j_m2)
                .map_err(|error| Error::Material { interface: endpoint.interface,error })?;
            old.next_closed(end.point)?;
        }
        let mut located = Vec::new();
        let mut trials = 0usize;
        for (&id, old) in &self.contacts.0 {
            let end = ends[&id];
            for kind in [CrossingKind::Failure,CrossingKind::OpeningGap,CrossingKind::ReclosingGap] {
                checkpoint()?;
                if kind == CrossingKind::Failure && old.point.failed() { continue; }
                let bracket = locate(start_s,target_s,config,&mut trials,|fraction| {
                    checkpoint()?;
                    let jump = jump_at(old.point.jump_m(),end,fraction);
                    Ok(match kind {
                        CrossingKind::Failure => old.point.trial(jump,config.local_energy_tolerance_j_m2)
                            .map_err(|error| Error::Material { interface: id,error })?.point.failed(),
                        CrossingKind::OpeningGap => jump[0] >= old.policy.open_gap_m,
                        CrossingKind::ReclosingGap => jump[0] <= old.policy.reclose_gap_m,
                    })
                })?;
                if let Some((low,high)) = bracket {
                    if located.len() >= config.max_crossings {
                        return Err(Error::Thermal(DepositionError::Budget("prescribed contact crossings")));
                    }
                    located.push(Located { id,kind,low,high });
                }
            }
        }
        located.sort_by(|a,b| a.high.total_cmp(&b.high).then(a.id.cmp(&b.id)).then(a.kind.cmp(&b.kind)));
        let crossings = located.iter().map(|c| Crossing { interface: c.id,kind: c.kind,
            before_s: stamp(start_s,target_s,c.low),after_s: stamp(start_s,target_s,c.high) }).collect();
        let mut cuts: Vec<_> = located.iter().map(|c| c.high).collect();
        cuts.push(1.0);
        cuts.dedup(); // Already sorted; the last crossing may be the final endpoint.
        let initial_energy = self.thermal.energy_j();
        let mut next = self.clone();
        let mut report = Report { thermal: EnergyAudit { incoming_mass_kg: 0.0,incoming_energy_j: 0.0,
            external_heat_j: 0.0,stored_energy_change_j: 0.0,residual_j: 0.0,substeps: 0,events: 0 },
            contacts: ContactReport::default(),crossings,event_trials: trials,thermal_intervals: 0 };
        for fraction in cuts {
            checkpoint()?;
            let time_s = stamp(start_s,target_s,fraction);
            let remaining = config.max_thermal_substeps-report.thermal.substeps;
            let (thermal,contacts,audit) = next.thermal.advance_observed_limited(
                &next.contacts,time_s,config.thermal,remaining,&mut checkpoint)?;
            if time_s > next.thermal.time_s() { report.thermal_intervals += 1; }
            add_thermal(&mut report.thermal,audit)?;
            next = Self { thermal,contacts };
            let mut proposed = Vec::with_capacity(ends.len());
            for (&id,&end) in &ends {
                checkpoint()?;
                let jump = jump_at(self.contacts.0[&id].point.jump_m(),end,fraction);
                let point = next.contacts.0[&id].point.trial(jump,config.local_energy_tolerance_j_m2)
                    .map_err(|error| Error::Material { interface: id,error })?.point;
                proposed.push((id,point));
            }
            let (updated,contact) = next.synchronize(time_s,&proposed,
                config.local_energy_tolerance_j_m2,&mut checkpoint)?;
            add_contacts(&mut report.contacts,contact)?;
            next = updated;
        }
        report.thermal.stored_energy_change_j = finite_value(next.thermal.energy_j()-initial_energy)?;
        report.thermal.residual_j = finite_value(report.thermal.stored_energy_change_j
            -report.thermal.incoming_energy_j-report.thermal.external_heat_j)?;
        if report.thermal.residual_j.abs() > config.thermal.energy_tolerance_j {
            return Err(Error::Thermal(DepositionError::EnergyBalance(report.thermal.residual_j)));
        }
        checkpoint()?;
        Ok((next,report))
    }
}

fn finite_value(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Thermal(DepositionError::Arithmetic)) }
}
fn add_count(a: usize,b: usize) -> Result<usize, Error> {
    a.checked_add(b).ok_or(Error::Thermal(DepositionError::Budget("prescribed path count overflow")))
}
fn add_thermal(sum: &mut EnergyAudit, part: EnergyAudit) -> Result<(), Error> {
    sum.incoming_mass_kg = finite_value(sum.incoming_mass_kg+part.incoming_mass_kg)?;
    sum.incoming_energy_j = finite_value(sum.incoming_energy_j+part.incoming_energy_j)?;
    sum.external_heat_j = finite_value(sum.external_heat_j+part.external_heat_j)?;
    sum.substeps = add_count(sum.substeps,part.substeps)?;
    sum.events = add_count(sum.events,part.events)?;
    Ok(())
}
fn add_contacts(sum: &mut ContactReport,part: ContactReport) -> Result<(), Error> {
    sum.opened = add_count(sum.opened,part.opened)?;
    sum.reclosed = add_count(sum.reclosed,part.reclosed)?;
    sum.changed_conductances = add_count(sum.changed_conductances,part.changed_conductances)?;
    sum.fracture_dissipation_j = finite_value(sum.fracture_dissipation_j+part.fracture_dissipation_j)?;
    Ok(())
}

#[cfg(test)]
mod tests;
