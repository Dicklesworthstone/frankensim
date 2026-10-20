//! Sample-addressed total forces on the retained original plate mesh.
//! Compilation is cold. Callbacks only select an admitted held generalized force;
//! pressure and force controls publish together AFTER their physical sample.
use super::{AperturePerformance, GestureCompileError, invalid};
use crate::bernoulli_aperture::dynamic::force::{PlateForceFootprint, PlateForcePort};
use crate::render::RenderError;
use crate::bernoulli_aperture::plate::PlateApertureReduction;

/// One signed physical force assignment, applied BEFORE its named source sample.
/// It remains held until the same port receives another assignment. Zero releases
/// the input only; it does not remove contact or erase material/wave history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ApertureForceEvent {
    /// Integer mechanical sample in the half-open finite performance window.
    pub sample: u64,
    /// Zero-based index in the supplied footprint list.
    pub port: usize,
    /// Total transverse force [N], positive along positive plate displacement.
    pub force_n: f64,
}

/// Fixed footprints and finite force histories. Distinct ports may overlap;
/// their physical loads add, in port order. All ports initially hold zero force.
#[derive(Clone, Debug)]
pub struct ApertureForceProgram {
    /// Exact original node or triangle addresses, never precompiled for another specimen.
    pub footprints: Vec<PlateForceFootprint>,
    /// May be unsorted. Equal-time assignments retain input order (last write wins).
    pub events: Vec<ApertureForceEvent>,
}
#[derive(Clone, Copy)]
struct Group { sample: u64, force_n: f64, end: usize }
pub(in crate::bernoulli_aperture) struct ForceSchedule {
    pub(in crate::bernoulli_aperture) ports: Vec<PlateForcePort>,
    pub(in crate::bernoulli_aperture) events: Vec<ApertureForceEvent>,
    groups: Vec<Group>,
    next: usize,
    pub(in crate::bernoulli_aperture) applied: usize,
    pub(in crate::bernoulli_aperture) held: f64,
    pub(in crate::bernoulli_aperture) work: u64,
}
impl ForceSchedule {
    // One projection/grouping owner. Callers pass REMAINING shared control/work
    // capacity, so independent blades do not each consume the entire allowance.
    pub(in crate::bernoulli_aperture) fn compile(plate: &PlateApertureReduction, program: ApertureForceProgram,
        samples: u64, max_events: usize, max_work: u64) -> Result<Self, GestureCompileError>
    {
        let count = program.footprints.len();
        if count == 0 || count > 32 || program.events.len() > max_events {
            return Err(invalid("mechanical program requires 1..=32 ports and combined controls within budget"));
        }
        let sites: u128 = program.footprints.iter().map(|p| match p {
            PlateForceFootprint::Node(_) => 1, PlateForceFootprint::Patch(t) => t.len() as u128,
        }).sum();
        let work = sites + program.events.len() as u128 * (count as u128 + 1);
        if work > u128::from(max_work) {
            return Err(invalid("mechanical force projection exceeds the compilation work budget"));
        }
        let physics = |e| GestureCompileError::Render(RenderError::Voice(e));
        let ports: Vec<_> = program.footprints.into_iter().map(|p| plate.force_port(p))
            .collect::<Result<_, _>>().map_err(physics)?;
        let (events, groups) = Self::groups(program.events, count, samples,
            |port, force| ports[port].generalized_force_n(force).map_err(physics))?;
        Ok(Self { ports, events, groups, next: 0, applied: 0, held: 0.0, work: work as u64 })
    }

    /// Keep force in the existing source port's physical newtons. Its retained
    /// multimode basis owns the force/velocity projection during the same solve;
    /// projecting to an opening coordinate here would discard the other modes.
    pub(in crate::bernoulli_aperture) fn source_port(events: Vec<ApertureForceEvent>,
        samples: u64, max_events: usize, max_work: u64) -> Result<Self, GestureCompileError>
    {
        let work = 2 * events.len() as u128;
        if events.len() > max_events || work > u128::from(max_work) {
            return Err(invalid("source force assignments exceed the shared compilation budget"));
        }
        let (events, groups) = Self::groups(events, 1, samples, |_, force| {
            if force.is_finite() { Ok(force) }
            else { Err(invalid("source force requires finite physical newtons")) }
        })?;
        Ok(Self { ports: Vec::new(), events, groups, next: 0, applied: 0, held: 0.0,
            work: work as u64 })
    }

    // The original stable grouping and preview/accept rules serve both explicit
    // one-coordinate footprints and already-projected material source ports.
    fn groups(mut events: Vec<ApertureForceEvent>, count: usize, samples: u64,
        mut project: impl FnMut(usize, f64) -> Result<f64, GestureCompileError>)
        -> Result<(Vec<ApertureForceEvent>, Vec<Group>), GestureCompileError>
    {
        for event in &events {
            if event.sample >= samples || event.port >= count {
                return Err(invalid("mechanical control names an absent port or a sample outside the performance"));
            }
            project(event.port, event.force_n)?;
        }
        events.sort_by_key(|e| e.sample);
        let mut groups = Vec::new();
        groups.try_reserve_exact(events.len())
            .map_err(|_| invalid("cannot reserve mechanical force controls"))?;
        let mut held = vec![0.0; count];
        let mut i = 0;
        while i < events.len() {
            let sample = events[i].sample;
            while i < events.len() && events[i].sample == sample {
                let e = events[i]; held[e.port] = e.force_n; i += 1;
            }
            // No intermediate same-sample assignment is a physical time step.
            // Sum only the final values, retaining signed cancellation and order.
            let mut force_n = 0.0;
            for (port, force) in held.iter().enumerate() {
                force_n += project(port, *force)?;
                if !force_n.is_finite() { return Err(invalid("sum of simultaneous mechanical forces overflowed")); }
            }
            groups.push(Group { sample, force_n, end: i });
        }
        Ok((events, groups))
    }

    /// Immutable sample preview. Clock admission guarantees no skipped group.
    pub(in crate::bernoulli_aperture) fn candidate(&self, sample: u64) -> f64 {
        self.groups.get(self.next).filter(|g| g.sample == sample).map_or(self.held, |g| g.force_n)
    }
    pub(in crate::bernoulli_aperture) fn accept(&mut self, sample: u64) {
        if let Some(g) = self.groups.get(self.next).filter(|g| g.sample == sample) {
            self.held = g.force_n; self.applied = g.end; self.next += 1;
        }
    }
}
impl AperturePerformance {
    /// Bind independent mechanical actuation to THIS unadvanced physical specimen.
    /// Pressure remains owned by the existing gesture schedule. Up to 32 ports
    /// are projected from its retained plate, not an externally substituted basis.
    /// Complete input admission and stable equal-time grouping occur before any
    /// state changes. The generalized sum enters the existing nonlinear solve.
    ///
    /// The existing max_controls bounds pressure, raw source-force and boundary
    /// assignments. max_compile_work bounds combined force projection/summation
    /// and boundary admission; pressure compilation retains its independent bound.
    /// No physical force magnitude is silently clipped to satisfy a solver.
    ///
    /// # Errors
    /// Missing plate, advanced/poisoned/already-bound performance, invalid site,
    /// event/clock/force, exhausted budget, or failed storage allocation.
    pub fn with_plate_forces(mut self, program: ApertureForceProgram)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.system.aperture().accepted_steps() != 0 || self.forces.is_some() {
            return Err(invalid("mechanical forces must bind once before performance playback"));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        let plate = self.system.aperture().plate_reduction()
            .ok_or_else(|| invalid("physical force footprints require a retained plate specimen"))?;
        self.forces = Some(ForceSchedule::compile(plate, program, self.config.samples,
            assignments, work)?);
        Ok(self)
    }
    /// Specimen-derived work-conjugate ports; empty for pressure-only playback.
    #[must_use]
    pub fn force_ports(&self) -> &[PlateForcePort] {
        self.forces.as_ref().map_or(&[], |f| f.ports.as_slice())
    }
    /// Raw physical assignments committed with successful source samples.
    #[must_use]
    pub fn applied_force_controls(&self) -> &[ApertureForceEvent] {
        self.forces.as_ref().map_or(&[], |f| &f.events[..f.applied])
    }
    /// Raw physical assignments not yet committed, in stable execution order.
    #[must_use]
    pub fn pending_force_controls(&self) -> &[ApertureForceEvent] {
        self.forces.as_ref().map_or(&[], |f| &f.events[f.applied..])
    }
    /// Current accepted opening-coordinate force [N], not a pressure or PCM gain.
    #[must_use]
    pub fn held_generalized_force_n(&self) -> f64 {
        self.forces.as_ref().map_or(0.0, |f| f.held)
    }
    /// Mechanical supply in the last accepted sample [J]; zero before playback.
    /// It is already included in the coupled owner's complete work balance.
    #[must_use]
    pub const fn last_mechanical_work_j(&self) -> f64 { self.mechanical_work_j }
}
