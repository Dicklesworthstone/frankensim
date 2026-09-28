//! Sample-addressed total forces on the retained original plate mesh.
//! Compilation is cold. Callbacks only select an admitted held generalized force;
//! pressure and force controls publish together AFTER their physical sample.
use super::{AperturePerformance, GestureCompileError, invalid};
use crate::bernoulli_aperture::dynamic::force::{PlateForceFootprint, PlateForcePort};
use crate::render::RenderError;

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
pub(super) struct ForceSchedule {
    ports: Vec<PlateForcePort>,
    events: Vec<ApertureForceEvent>,
    groups: Vec<Group>,
    next: usize,
    applied: usize,
    held: f64,
}
impl ForceSchedule {
    /// Immutable sample preview. Clock admission guarantees no skipped group.
    pub(super) fn candidate(&self, sample: u64) -> f64 {
        self.groups.get(self.next).filter(|g| g.sample == sample).map_or(self.held, |g| g.force_n)
    }
    pub(super) fn accept(&mut self, sample: u64) {
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
    /// The existing max_controls bounds pressure assignments plus raw force
    /// assignments. max_compile_work separately bounds force projection/summation
    /// visits; it also remains the pressure compiler's independent work bound.
    /// No physical force magnitude is silently clipped to satisfy a solver.
    ///
    /// # Errors
    /// Missing plate, advanced/poisoned/already-bound performance, invalid site,
    /// event/clock/force, exhausted budget, or failed storage allocation.
    pub fn with_plate_forces(mut self, mut program: ApertureForceProgram)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.system.aperture().accepted_steps() != 0 || self.forces.is_some() {
            return Err(invalid("mechanical forces must bind once before performance playback"));
        }
        let count = program.footprints.len();
        if count == 0 || count > 32 || program.events.len() > self.config.max_controls.saturating_sub(self.controls.len()) {
            return Err(invalid("mechanical program requires 1..=32 ports and combined controls within budget"));
        }
        let plate = self.system.aperture().plate_reduction()
            .ok_or_else(|| invalid("physical force footprints require a retained plate specimen"))?;
        let sites: u128 = program.footprints.iter().map(|p| match p {
            PlateForceFootprint::Node(_) => 1, PlateForceFootprint::Patch(t) => t.len() as u128,
        }).sum();
        let work = sites + program.events.len() as u128 * (count as u128 + 1);
        if work > u128::from(self.config.max_compile_work) {
            return Err(invalid("mechanical force projection exceeds the compilation work budget"));
        }
        let physics = |e| GestureCompileError::Render(RenderError::Voice(e));
        let ports: Vec<_> = program.footprints.into_iter().map(|p| plate.force_port(p))
            .collect::<Result<_, _>>().map_err(physics)?;
        for event in &program.events {
            if event.sample >= self.config.samples || event.port >= count {
                return Err(invalid("mechanical control names an absent port or a sample outside the performance"));
            }
            ports[event.port].generalized_force_n(event.force_n).map_err(physics)?;
        }
        program.events.sort_by_key(|e| e.sample);
        let mut groups = Vec::new();
        groups.try_reserve_exact(program.events.len())
            .map_err(|_| invalid("cannot reserve mechanical force controls"))?;
        let mut held = vec![0.0; count];
        let mut i = 0;
        while i < program.events.len() {
            let sample = program.events[i].sample;
            while i < program.events.len() && program.events[i].sample == sample {
                let e = program.events[i]; held[e.port] = e.force_n; i += 1;
            }
            // No intermediate same-sample assignment is a physical time step.
            // Sum only the final values, retaining signed cancellation and order.
            let mut force_n = 0.0;
            for (port, force) in ports.iter().zip(&held) {
                force_n += port.generalized_force_n(*force).map_err(physics)?;
                if !force_n.is_finite() { return Err(invalid("sum of simultaneous mechanical forces overflowed")); }
            }
            groups.push(Group { sample, force_n, end: i });
        }
        self.forces = Some(ForceSchedule { ports, events: program.events, groups, next: 0, applied: 0, held: 0.0 });
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
