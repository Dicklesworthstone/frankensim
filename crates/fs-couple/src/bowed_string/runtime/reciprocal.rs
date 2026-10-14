//! Material-bound strings in the existing two-way modal/friction network.
//!
//! The specimen spans TWO fixed end supports. An explicitly located interior
//! connection can represent a seated bridge with an afterlength, a compliant
//! fret, or another bilateral support. It cannot represent a moving endpoint:
//! pinned shapes vanish there, and pretending otherwise would add a false port.
//! The same signed displacement shape distributes the connection reaction back
//! to the string and body. The body is not a one-way sound filter.
//!
//! Dispersion, material bending loss and cylinder-air drag come from the existing
//! material bowed-card lowering. Bilateral reactions and friction come from the
//! existing coupled-modal owners. No new string, friction, or integration law.
//! This is fixed-state small-amplitude bending with a fixed bow station. It does
//! not add torsion, polarization, nonlinear tension, detached bridge contact,
//! radiation back-loading, or experimentally validated instrument timbre.

use crate::acoustic_realize::AcousticRealizeError;
use crate::bowed_string::{BowedStringCard, unit_shape};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::{GestureCompileError, ScheduledRenderer};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionConfig, ModalFrictionError, ModalFrictionPort,
};
use crate::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use crate::string_specimen::ResolvedStringSpecimen;
use fs_exec::{CancelGate, Cx};
use fs_scenario::{AmbientGas, gesture::GestureSchedule};
use fs_tribo::FrictionLaw;

/// A physical bilateral connection at an INTERIOR point on the full string.
/// The supplied body shape must use its actual mass-normalized displacement
/// basis. Its provenance/geometry remains the body producer's responsibility.
#[derive(Clone, Debug, PartialEq)]
pub struct InteriorStringConnection {
    /// Fraction of the complete supported span, strictly inside (0,1).
    pub station_fraction: f64,
    /// Index into the supplied body list, NOT the final network component list.
    /// The assembled network orders the string first, then these bodies.
    pub body: ModalAttachment,
    /// Nonnegative transverse stiffness [N/m], independent of bulk string E.
    pub stiffness_n_m: f64,
    /// Nonnegative transverse viscous resistance [N s/m].
    pub damping_n_s_m: f64,
    /// Stress-free string-minus-body displacement [m].
    pub rest_extension_m: f64,
}

/// Explicit resolution, attachment geometry and original numerical-owner limits.
#[derive(Clone, Debug)]
pub struct ReciprocalStringOptions {
    /// Shared mechanical/audio rate; every body must already use this rate.
    pub sample_rate_hz: u32,
    /// Fixed bow station on the complete string, strictly inside (0,1).
    pub bow_station_fraction: f64,
    /// At least one nonzero physical connection to a supplied body.
    pub connections: Vec<InteriorStringConnection>,
    /// Limits for the material-derived string component, without fake pressure.
    pub string_budget: ModalAcousticTimeBudget,
    /// Existing network's setup, coupled-bandwidth, force and energy admission.
    pub network: ModalCouplingConfig,
    /// Existing implicit friction's physical and nonlinear-solve admission.
    pub friction: ModalFrictionConfig,
}

/// Exact material query/geometry and authored environment/connection inputs.
/// This is not new source authority for the body modes or friction coefficients.
#[derive(Clone, Debug)]
pub struct ReciprocalStringSource {
    specimen: ResolvedStringSpecimen,
    ambient: AmbientGas,
    options: ReciprocalStringOptions,
}
impl ReciprocalStringSource {
    /// Complete original specimen, including material identities and receipts.
    #[must_use]
    pub const fn specimen(&self) -> &ResolvedStringSpecimen { &self.specimen }
    /// Original fixed environment supplied to the common drag/loss lowering.
    #[must_use]
    pub const fn ambient(&self) -> &AmbientGas { &self.ambient }
    /// Original geometry, sample clock and declared numerical/physical limits.
    #[must_use]
    pub const fn options(&self) -> &ReciprocalStringOptions { &self.options }
}

/// Admission failures retain their original physical owner.
#[derive(Debug)]
pub enum ReciprocalStringError {
    /// Material/geometry/loss lowering refused the specimen.
    Material(AcousticRealizeError),
    /// Bilateral network or cancellation refused admission.
    Network(ModalCouplingError),
    /// The unchanged implicit friction owner refused admission.
    Friction(ModalFrictionError),
}
impl core::fmt::Display for ReciprocalStringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "material string: {error:?}"),
            Self::Network(error) => write!(f, "string/body network: {error}"),
            Self::Friction(error) => write!(f, "string/body friction: {error}"),
        }
    }
}
impl std::error::Error for ReciprocalStringError {}

fn invalid(what: &'static str) -> ReciprocalStringError {
    ReciprocalStringError::Network(ModalCouplingError::Invalid { what })
}
fn poll(gate: &CancelGate) -> Result<(), ReciprocalStringError> {
    if gate.is_requested() {
        Err(ReciprocalStringError::Network(ModalCouplingError::Cancelled))
    } else { Ok(()) }
}
fn interior(station: f64) -> bool { station.is_finite() && station > 0.0 && station < 1.0 }

/// Retained two-way network, not an alternative bow or modal time stepper.
/// String material changes rebuild its modes AND force/displacement ports, then
/// the original network owner rebuilds the reaction solve over those modes.
pub struct ReciprocalMaterialString {
    source: ReciprocalStringSource,
    system: FrictionModalSystem,
}
impl ReciprocalMaterialString {
    /// Assemble the material string, supplied body modes and actual connections.
    ///
    /// The string initially rests; any existing body vibration is preserved.
    /// Pressure comes only from the supplied body observations: transverse string
    /// velocity is never relabeled as pascals. The complete network reacts back
    /// during both bilateral and friction solves. Radiation loading is whatever
    /// the body producer modeled, not something inferred from its pressure port.
    ///
    /// For a bridge, use the complete fixed-end span INCLUDING afterlength and
    /// an interior bridge station. A speaking-length-only pinned model refuses
    /// an endpoint connection instead of pretending its bridge can move.
    ///
    /// # Errors
    /// Missing bodies/connections, endpoint/zero ports, incompatible bases/rates,
    /// material-law applicability, aliasing, numerical budgets or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_specimen(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, ambient: AmbientGas,
        bodies: Vec<ModalAcousticTimeModel>, law: FrictionLaw,
        options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if bodies.is_empty() || options.connections.is_empty()
            || !interior(options.bow_station_fraction) {
            return Err(invalid("reciprocal string requires bodies, connections and an interior bow station"));
        }
        let n = specimen.string().n_modes;
        let total = bodies.iter().try_fold(n, |sum, body| sum.checked_add(body.modes().len()))
            .ok_or_else(|| invalid("reciprocal string mode count overflow"))?;
        if total > options.network.max_modes || options.connections.len() > options.network.max_connections {
            return Err(invalid("reciprocal string exceeds its declared setup dimensions"));
        }
        let card = BowedStringCard::from_material_specimen(cx, &specimen, ambient, options.sample_rate_hz)
            .map_err(ReciprocalStringError::Material)?;
        let shape = |station| -> Vec<f64> {
            (0..card.mode_count).map(|k| unit_shape(k, station, card.length_m, card.linear_density_kg_m)).collect()
        };
        let mut links = Vec::with_capacity(options.connections.len());
        for connection in &options.connections {
            poll(gate)?;
            if !interior(connection.station_fraction) {
                return Err(invalid("a pinned string endpoint cannot be a moving body attachment; supply an interior station and afterlength"));
            }
            let body = bodies.get(connection.body.component)
                .ok_or_else(|| invalid("string connection names an unknown supplied body"))?;
            if connection.body.shapes.len() != body.modes().len()
                || !connection.body.shapes.iter().all(|x| x.is_finite())
                || !connection.body.shapes.iter().any(|x| *x != 0.0) {
                return Err(invalid("reciprocal body attachment must be a nonzero finite shape in its original basis"));
            }
            if !(connection.stiffness_n_m > 0.0 || connection.damping_n_s_m > 0.0) {
                return Err(invalid("reciprocal string connection cannot have both stiffness and damping zero"));
            }
            links.push(ModalConnection {
                left: ModalAttachment { component: 0, shapes: shape(connection.station_fraction) },
                right: ModalAttachment {
                    component: connection.body.component + 1,
                    shapes: connection.body.shapes.clone(),
                },
                stiffness_n_m: connection.stiffness_n_m,
                damping_n_s_m: connection.damping_n_s_m,
                rest_extension_m: connection.rest_extension_m,
            });
        }
        let modes = (0..card.mode_count).map(|k| ModalAcousticMode {
            angular_frequency_rad_s: card.mode_omega_rad_s(k + 1),
            damping_ratio: card.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        }).collect();
        let string = ModalAcousticTimeModel::try_new(options.sample_rate_hz, modes, options.string_budget)
            .map_err(|source| ReciprocalStringError::Network(ModalCouplingError::Component { component: 0, source }))?;
        let mut components = Vec::with_capacity(bodies.len() + 1);
        components.push(string);
        components.extend(bodies);
        let network = CoupledModalSystem::new(components, links, options.network, gate)
            .map_err(ReciprocalStringError::Network)?;
        let port = ModalFrictionPort {
            left: ModalAttachment { component: 0, shapes: shape(options.bow_station_fraction) },
            right: None, law,
        };
        let system = FrictionModalSystem::new(network, port, options.friction, gate)
            .map_err(ReciprocalStringError::Friction)?;
        poll(gate)?;
        Ok(Self { source: ReciprocalStringSource { specimen, ambient, options }, system })
    }

    /// Material/source data remain attached to the accepted mechanics.
    #[must_use]
    pub const fn source(&self) -> &ReciprocalStringSource { &self.source }
    /// Existing owner exposes accepted modes, pressure, work/loss and slip checks.
    #[must_use]
    pub const fn system(&self) -> &FrictionModalSystem { &self.system }
    /// Advance through the original owner's transactional step or gated step.
    pub fn system_mut(&mut self) -> &mut FrictionModalSystem { &mut self.system }
    /// Transfer to another existing renderer while explicitly retaining sources.
    #[must_use]
    pub fn into_parts(self) -> (ReciprocalStringSource, FrictionModalSystem) { (self.source, self.system) }

    /// Bind an existing physical bow track without dropping the material source.
    /// The already-assembled bow station and clock must match the requested track
    /// compilation. All state, pressure observations and control semantics remain
    /// in ScheduledRenderer; this adapter adds no sample or interpolation loop.
    pub fn into_scheduled(
        self, schedule: &GestureSchedule, track: &str, config: FrictionGestureConfig,
    ) -> Result<ReciprocalStringPerformance, GestureCompileError> {
        if config.sample_rate_hz != self.source.options.sample_rate_hz
            || config.station_fraction.to_bits() != self.source.options.bow_station_fraction.to_bits() {
            return Err(GestureCompileError::Render(crate::render::RenderError::Control {
                what: "material string schedule must match its actual fixed bow port and sample clock",
            }));
        }
        let external = vec![0.0; self.system.mode_count()];
        let renderer = ScheduledRenderer::from_friction_gesture(self.system, external, schedule, track, config)?;
        Ok(ReciprocalStringPerformance { source: self.source, renderer })
    }
}

/// Original material/geometry paired with the common sample-accurate renderer.
pub struct ReciprocalStringPerformance {
    source: ReciprocalStringSource,
    renderer: ScheduledRenderer,
}
impl ReciprocalStringPerformance {
    /// Exact material and mechanical input binding used to build this performance.
    #[must_use]
    pub const fn source(&self) -> &ReciprocalStringSource { &self.source }
    /// Read-only access to the unchanged scheduled renderer.
    #[must_use]
    pub const fn renderer(&self) -> &ScheduledRenderer { &self.renderer }
    /// Render with the existing callback, control and cancellation semantics.
    pub fn renderer_mut(&mut self) -> &mut ScheduledRenderer { &mut self.renderer }
}

#[cfg(test)]
mod tests;
