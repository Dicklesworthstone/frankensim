//! Physical string loads and load-release plucks on the SAME bowed instrument.
//!
//! Source geometry supplies force/displacement maps. The original force compiler
//! supplies stable grouped controls; the original coupled solver supplies static
//! preload, motion, friction, work and pressure. No additional physics stepper.

use super::{MultiBowPerformance, MultiBowStringBank, render_invalid};
use crate::bowed_string::runtime::reciprocal::{ReciprocalStringError, bank::StringBowLocation};
use crate::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirection;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent};
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::program::{PhysicalForceInput, ProjectedForceProgram};
use fs_exec::CancelGate;

/// Which physical motion is work-conjugate to the applied transverse force.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StringForceCoupling {
    /// Transverse centreline force along the source's bow direction in both
    /// selected planes, with zero applied torque.
    Centreline,
    /// Centreline force along its own unit direction in the source's transverse
    /// frame, independent of bow/bridge directions. It applies no torque.
    /// A nonzero second component requires the original two-plane source.
    CentrelineAlong(TransverseDirection),
    /// Force at the same signed outer surface used by bow contact. Includes
    /// radius times angular velocity/torque when the source has torsional modes.
    Surface,
}
/// A fixed physical load location, never an authored modal participation vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringForcePort {
    /// Original material-string index and interior station on its COMPLETE span.
    pub location: StringBowLocation,
    /// Centreline translation or the original bow-facing surface motion.
    pub coupling: StringForceCoupling,
}
/// Original physical inputs retained alongside material, gas and interface data.
#[derive(Clone, Debug)]
pub struct StringForceProgram {
    /// Independent load ports. Distinct ports at the same location SUM forces.
    pub ports: Vec<StringForcePort>,
    /// Initial held newtons, exactly one per port. Not a pitch or strike velocity.
    pub initial_force_n: Vec<f64>,
    /// Existing physical-force events. `voice` must be ZERO: this is one coupled
    /// instrument, even with many strings. Samples are local to the performance;
    /// `port` addresses `ports`. Equal-time assignments to one port use the last.
    pub events: Vec<ModalForceEvent>,
    /// Retain motion (a real load-on transient) or settle the entire initially
    /// unexcited instrument BEFORE sample zero. A sample-zero event may release
    /// the latter. Preload includes attached stop laws; it does not infer tissue.
    pub initialization: ForceInitialization,
    /// Existing force compiler's rate, callback, authored-event, emitted-control
    /// and projection-work limits. Rate/callback must match this performance.
    /// Raw input events and compiled controls are both retained within these caps.
    pub limits: ForceRenderConfig,
}

pub(super) struct BoundStringForces {
    pub(super) source: StringForceProgram,
    pub(super) projected: ProjectedForceProgram,
    pub(super) candidate: Vec<f64>,
    pub(super) next: usize,
    pub(super) preload_energy_j: Option<f64>,
}
fn poll(gate: &CancelGate) -> Result<(), RenderError> {
    if gate.is_requested() { Err(RenderError::Coupled(ModalCouplingError::Cancelled)) } else { Ok(()) }
}
fn source_error(error: ReciprocalStringError) -> RenderError {
    match error {
        ReciprocalStringError::Material(e) => RenderError::Voice(e),
        ReciprocalStringError::Network(e) => RenderError::Coupled(e),
        ReciprocalStringError::Friction(e) => RenderError::Friction(e),
    }
}
impl MultiBowStringBank {
    /// Obtain physical load/displacement columns in this instrument's actual
    /// flattened basis, including torsion where explicitly requested. Body and
    /// cavity entries are zero DIRECT loads; their reactions remain solved by
    /// the shared network. `F dot (B^T dq)` and `(B F) dot dq` are the same work.
    /// No source/material queries or modal reductions are repeated at runtime.
    ///
    /// This is also usable with direct `step` calls and the shared force compiler.
    /// Setup size is checked BEFORE allocating the full columns. Refuses unknown
    /// strings, endpoint/nonfinite stations, dimension/work limits or cancellation.
    pub fn string_force_columns(&self, ports: &[StringForcePort], max_terms: usize,
        gate: &CancelGate) -> Result<Vec<Vec<f64>>, RenderError>
    {
        poll(gate)?;
        let n = self.system.mode_count();
        let terms = n.checked_mul(ports.len()).ok_or_else(|| render_invalid("string force map size overflow"))?;
        if ports.is_empty() || terms > max_terms {
            return Err(RenderError::Sizing { what: "string force maps exceed the declared nonempty projection budget" });
        }
        let sources = self.source.bank().strings();
        let mut columns = Vec::new();
        columns.try_reserve(ports.len()).map_err(|_| render_invalid("cannot reserve string force maps"))?;
        for port in ports {
            poll(gate)?;
            let index = port.location.string;
            let station = port.location.station_fraction;
            let source = sources.get(index).ok_or_else(|| render_invalid("force names an unknown material string"))?;
            if !station.is_finite() || station <= 0.0 || station >= 1.0 {
                return Err(render_invalid("string force needs an actual interior point, not a pinned endpoint"));
            }
            let count = self.system.components()[index].modes().len();
            let card = &self.cards[index];
            let local = match port.coupling {
                StringForceCoupling::Surface => source.source.surface_shape(card, station).map_err(source_error)?,
                StringForceCoupling::Centreline => source.source.transverse_shape(card, station).map_err(source_error)?,
                StringForceCoupling::CentrelineAlong(direction) => source.source
                    .transverse_shape_along(card, station, direction).map_err(source_error)?,
            };
            if local.len() != count || local.iter().any(|x| !x.is_finite()) {
                return Err(render_invalid("string force map differs from its actual finite basis"));
            }
            let offset: usize = self.system.components()[..index].iter().map(|m| m.modes().len()).sum();
            let mut column = vec![0.0; n];
            column[offset..offset+count].copy_from_slice(&local);
            columns.push(column);
        }
        poll(gate)?;
        Ok(columns)
    }
}
impl MultiBowPerformance {
    /// Add material-derived plucks, taps or other prescribed held loads to the
    /// existing scheduled bow performance. All due bow AND force controls act
    /// before one joint step. No mode, body, cavity or interface is reconstructed.
    ///
    /// A static pluck uses the original full-network equilibrium solver, including
    /// bridge/body/cavity and attached stop reactions, with all bows released before the window.
    /// It never initializes each component independently. `RetainState` instead
    /// preserves any existing vibration and applies a load-on transient.
    ///
    /// All authored events, including future events, are admitted before preload.
    /// Bind only once before local sample zero. Preload additionally requires an
    /// entirely unexcited instrument at absolute sample zero. The normal contact
    /// model, nonlinear tension and physical finger disengagement are not inferred
    /// from a prescribed load program. Pressure inherits the existing observer.
    pub fn with_string_forces(mut self, program: StringForceProgram, gate: &CancelGate)
        -> Result<Self, RenderError>
    {
        poll(gate)?;
        self.validate_sample_count(0)?;
        if self.samples_rendered() != 0 || self.string_forces.is_some()
            || program.limits.sample_rate_hz != self.bank.source.bank().options().sample_rate_hz
            || program.limits.max_block != self.config.max_block
            || program.ports.len() != program.initial_force_n.len()
            || program.initial_force_n.iter().any(|f| !f.is_finite()) {
            return Err(render_invalid("string forces need one initial load per port and a single unstarted binding on the original clock"));
        }
        if program.events.len() > program.limits.max_events {
            return Err(RenderError::Sizing { what: "string force source event budget exceeded" });
        }
        if program.events.iter().any(|e| e.voice != 0 || e.sample >= self.config.samples
            || e.port >= program.ports.len() || !e.force_n.is_finite()) {
            return Err(render_invalid("string force event must name this instrument, an actual load port and a sample inside its finite window"));
        }
        let columns = self.bank.string_force_columns(&program.ports, program.limits.max_projection_terms, gate)?;
        let input = PhysicalForceInput::new(columns, program.initial_force_n.clone())?;
        let projected = ProjectedForceProgram::compile(vec![input], program.events.clone(), program.limits)?;
        let initial = projected.initial_forces()[0].clone();
        let candidate = initial.clone();
        poll(gate)?;
        // No fallible work after the original transactional preload publishes.
        let preload_energy_j = if program.initialization == ForceInitialization::StaticPreload {
            Some(self.bank.system.initialize_released_equilibrium(&initial, gate).map_err(RenderError::Friction)?)
        } else { None };
        self.forces.copy_from_slice(&initial);
        self.string_forces = Some(BoundStringForces { source: program, projected, candidate,
            next: 0, preload_energy_j });
        Ok(self)
    }
    /// Original prescribed load locations, forces, initialization and limits.
    #[must_use]
    pub fn string_force_program(&self) -> Option<&StringForceProgram> {
        self.string_forces.as_ref().map(|p| &p.source)
    }
    /// Exact controls prepared by the common force compiler; not mutable gains.
    #[must_use]
    pub fn projected_string_forces(&self) -> Option<&ProjectedForceProgram> {
        self.string_forces.as_ref().map(|p| &p.projected)
    }
    /// Full initial stored energy only when explicit static preload was selected.
    /// This is not claimed as work supplied inside the subsequent audio window.
    #[must_use]
    pub fn preload_energy_j(&self) -> Option<f64> {
        self.string_forces.as_ref().and_then(|p| p.preload_energy_j)
    }
}