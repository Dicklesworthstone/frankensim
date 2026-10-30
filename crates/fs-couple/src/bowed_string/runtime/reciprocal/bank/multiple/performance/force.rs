//! Physical string and shell-bridge loads on the SAME bowed instrument.
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
/// A source-bound physical load on the original coupled instrument.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialForcePort {
    /// A fixed centreline or surface point on an actual material string.
    String(StringForcePort),
    /// Original shell or shell-assembly bridge force map, in source order.
    /// Input is newtons along that declared load pattern, including any original
    /// moment-per-newton lever arms; no direction or torque is inferred.
    ShellBridge { bridge: usize },
}
/// Original physical inputs retained alongside material, gas and interface data.
#[derive(Clone, Debug)]
pub struct PhysicalForceProgram<P> {
    /// Independent load ports. Distinct ports at the same location SUM forces.
    pub ports: Vec<P>,
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

/// Existing string-only authoring API, retaining its original source type.
pub type StringForceProgram = PhysicalForceProgram<StringForcePort>;
/// Mixed source-bound string and shell-bridge loads on one physical clock.
pub type MaterialForceProgram = PhysicalForceProgram<MaterialForcePort>;

enum BoundForceSource {
    Strings(StringForceProgram),
    Material(MaterialForceProgram),
}
struct ForceControls<'a> {
    ports: usize,
    initial_force_n: &'a [f64],
    events: &'a [ModalForceEvent],
    initialization: ForceInitialization,
    limits: ForceRenderConfig,
}
impl<P> PhysicalForceProgram<P> {
    fn controls(&self) -> ForceControls<'_> {
        ForceControls { ports: self.ports.len(), initial_force_n: &self.initial_force_n,
            events: &self.events, initialization: self.initialization, limits: self.limits }
    }
}
impl BoundForceSource {
    fn controls(&self) -> ForceControls<'_> {
        match self { Self::Strings(p) => p.controls(), Self::Material(p) => p.controls() }
    }
}

pub(super) struct BoundStringForces {
    source: BoundForceSource,
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
    /// the shared network. Force and displacement maps are work-conjugate.
    ///
    /// Setup size is checked BEFORE allocating full columns. Source material
    /// maps are prepared once; no queries or reductions enter the audio callback.
    pub fn string_force_columns(&self, ports: &[StringForcePort], max_terms: usize,
        gate: &CancelGate) -> Result<Vec<Vec<f64>>, RenderError>
    {
        self.force_columns(ports.iter().copied().map(MaterialForcePort::String), max_terms, gate)
    }

    /// Physical columns for mixed string and original shell bridge actuators.
    /// A bridge row loads only the retained structural prefix of the actual body
    /// component. Strings, appended cavity coordinates and all other components
    /// receive exactly zero DIRECT load; the original reciprocal network solves
    /// their responses. No modal participation values may be authored here.
    ///
    /// Refuses a missing source shell, absent bridge, incompatible extent,
    /// cancellation or a setup budget overrun. Residual-flexibility condensation
    /// requires its own complete direct-force columns and is explicitly refused.
    pub fn material_force_columns(&self, ports: &[MaterialForcePort], max_terms: usize,
        gate: &CancelGate) -> Result<Vec<Vec<f64>>, RenderError>
    {
        self.force_columns(ports.iter().copied(), max_terms, gate)
    }

    fn force_columns(&self, ports: impl ExactSizeIterator<Item = MaterialForcePort>,
        max_terms: usize, gate: &CancelGate) -> Result<Vec<Vec<f64>>, RenderError>
    {
        poll(gate)?;
        let n = self.system.mode_count();
        let terms = n.checked_mul(ports.len()).ok_or_else(|| render_invalid("physical force map size overflow"))?;
        if ports.len() == 0 || terms > max_terms {
            return Err(RenderError::Sizing { what: "physical force maps exceed the declared nonempty projection budget" });
        }
        let mut columns = Vec::new();
        columns.try_reserve(ports.len()).map_err(|_| render_invalid("cannot reserve physical force maps"))?;
        for port in ports {
            poll(gate)?;
            let (index, local) = match port {
                MaterialForcePort::String(port) => self.string_force_map(port)?,
                MaterialForcePort::ShellBridge { bridge } => self.shell_bridge_force_map(bridge)?,
            };
            let components = self.system.components();
            let count = components.get(index).ok_or_else(|| render_invalid("physical force lost its original component"))?
                .modes().len();
            if local.is_empty() || local.len() > count || local.iter().any(|x| !x.is_finite()) {
                return Err(render_invalid("physical force map differs from its actual finite basis"));
            }
            let offset = components[..index].iter().try_fold(0usize,
                |sum, m| sum.checked_add(m.modes().len()))
                .ok_or_else(|| render_invalid("physical force component offset overflow"))?;
            let end = offset.checked_add(local.len())
                .filter(|&end| end <= n).ok_or_else(|| render_invalid("physical force extent exceeds the original system"))?;
            let mut column = vec![0.0; n];
            column[offset..end].copy_from_slice(&local);
            columns.push(column);
        }
        poll(gate)?;
        Ok(columns)
    }

    fn string_force_map(&self, port: StringForcePort) -> Result<(usize, Vec<f64>), RenderError> {
        let index = port.location.string;
        let station = port.location.station_fraction;
        let source = self.source.bank().strings().get(index)
            .ok_or_else(|| render_invalid("force names an unknown material string"))?;
        if !station.is_finite() || station <= 0.0 || station >= 1.0 {
            return Err(render_invalid("string force needs an actual interior point, not a pinned endpoint"));
        }
        let count = self.system.components().get(index)
            .ok_or_else(|| render_invalid("string force lost its original component"))?.modes().len();
        let card = self.cards.get(index).ok_or_else(|| render_invalid("string force lost its original material card"))?;
        let local = match port.coupling {
            StringForceCoupling::Surface => source.source.surface_shape(card, station).map_err(source_error)?,
            StringForceCoupling::Centreline => source.source.transverse_shape(card, station).map_err(source_error)?,
            StringForceCoupling::CentrelineAlong(direction) => source.source
                .transverse_shape_along(card, station, direction).map_err(source_error)?,
        };
        if local.len() != count {
            return Err(render_invalid("string force map differs from its actual finite basis"));
        }
        Ok((index, local))
    }

    fn shell_bridge_force_map(&self, bridge: usize) -> Result<(usize, Vec<f64>), RenderError> {
        let source = self.source.bank();
        if source.residual_flexibility() {
            return Err(render_invalid("bridge loads require explicit residual-flexibility direct-force columns"));
        }
        let (ports, exterior, cavity) = match (source.material_shell_image(), source.material_shell_assembly_image()) {
            (Some(image), None) => (image.force_ports(), image.volume_shapes_m2_per_sqrt_kg(),
                image.cavity_volume_shapes_m2_per_sqrt_kg()),
            (None, Some(image)) => (image.force_ports(), image.volume_shapes_m2_per_sqrt_kg(),
                image.cavity_volume_shapes_m2_per_sqrt_kg()),
            _ => return Err(render_invalid("bridge load requires one original material shell image")),
        };
        let row = ports.get(bridge).ok_or_else(|| render_invalid("bridge force names an absent original shell port"))?;
        if row.is_empty() || row.len() != exterior.len() || row.len() != cavity.len()
            || ports.iter().any(|port| port.len() != row.len()) {
            return Err(render_invalid("bridge load lost its original shell force and volume maps"));
        }
        Ok((source.strings().len(), row.clone()))
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
    pub fn with_string_forces(self, program: StringForceProgram, gate: &CancelGate)
        -> Result<Self, RenderError>
    {
        self.bind_force_program(BoundForceSource::Strings(program), gate)
    }

    /// Add prescribed forces at original shell bridges, string points or both.
    /// Held loads and finite pulse/ramp controls use the SAME force compiler,
    /// preload and joint mechanics as string plucks. All due forces and bow
    /// controls act before one ordinary or thermal mechanical step, with the
    /// existing shell observer and coating admission guard remaining in force.
    ///
    /// A finite pulse is expressed by distinct load-on/load-off samples in
    /// newtons. A ramp uses authored sampled force assignments; neither silently
    /// substitutes an ideal impulse or rescales force by the output block size.
    /// Source maps are original mass-normalized displacement/load columns.
    pub fn with_material_forces(self, program: MaterialForceProgram, gate: &CancelGate)
        -> Result<Self, RenderError>
    {
        self.bind_force_program(BoundForceSource::Material(program), gate)
    }

    fn bind_force_program(mut self, source: BoundForceSource, gate: &CancelGate)
        -> Result<Self, RenderError>
    {
        let program = source.controls();
        poll(gate)?;
        self.validate_sample_count(0)?;
        if self.samples_rendered() != 0 || self.string_forces.is_some()
            || program.limits.sample_rate_hz != self.bank.source.bank().options().sample_rate_hz
            || program.limits.max_block != self.config.max_block
            || program.ports != program.initial_force_n.len()
            || program.initial_force_n.iter().any(|f| !f.is_finite()) {
            return Err(render_invalid("physical forces need one initial load per port and a single unstarted binding on the original clock"));
        }
        if program.events.len() > program.limits.max_events {
            return Err(RenderError::Sizing { what: "physical force source event budget exceeded" });
        }
        if program.events.iter().any(|e| e.voice != 0 || e.sample >= self.config.samples
            || e.port >= program.ports || !e.force_n.is_finite()) {
            return Err(render_invalid("physical force event must name this instrument, an actual load port and a sample inside its finite window"));
        }
        let columns = match &source {
            BoundForceSource::Strings(p) => self.bank.string_force_columns(&p.ports, p.limits.max_projection_terms, gate)?,
            BoundForceSource::Material(p) => self.bank.material_force_columns(&p.ports, p.limits.max_projection_terms, gate)?,
        };
        let input = PhysicalForceInput::new(columns, program.initial_force_n.to_vec())?;
        let projected = ProjectedForceProgram::compile(vec![input], program.events.to_vec(), program.limits)?;
        let initial = projected.initial_forces()[0].clone();
        let candidate = initial.clone();
        poll(gate)?;
        // No fallible work after the original transactional preload publishes.
        let preload_energy_j = if program.initialization == ForceInitialization::StaticPreload {
            Some(self.bank.system.initialize_released_equilibrium(&initial, gate).map_err(RenderError::Friction)?)
        } else { None };
        self.forces.copy_from_slice(&initial);
        self.string_forces = Some(BoundStringForces { source, projected, candidate,
            next: 0, preload_energy_j });
        Ok(self)
    }
    /// Original prescribed load locations, forces, initialization and limits.
    #[must_use]
    pub fn string_force_program(&self) -> Option<&StringForceProgram> {
        match &self.string_forces.as_ref()?.source {
            BoundForceSource::Strings(p) => Some(p),
            BoundForceSource::Material(_) => None,
        }
    }
    /// Original mixed physical inputs, if bound through the material-force API.
    #[must_use]
    pub fn material_force_program(&self) -> Option<&MaterialForceProgram> {
        match &self.string_forces.as_ref()?.source {
            BoundForceSource::Material(p) => Some(p),
            BoundForceSource::Strings(_) => None,
        }
    }
    /// Exact controls prepared by the common force compiler; not mutable gains.
    #[must_use]
    pub fn projected_string_forces(&self) -> Option<&ProjectedForceProgram> {
        self.projected_material_forces()
    }
    /// Common compiled controls for either physical authoring API.
    #[must_use]
    pub fn projected_material_forces(&self) -> Option<&ProjectedForceProgram> {
        self.string_forces.as_ref().map(|p| &p.projected)
    }
    /// Full initial stored energy only when explicit static preload was selected.
    /// This is not claimed as work supplied inside the subsequent audio window.
    #[must_use]
    pub fn preload_energy_j(&self) -> Option<f64> {
        self.string_forces.as_ref().and_then(|p| p.preload_energy_j)
    }
}