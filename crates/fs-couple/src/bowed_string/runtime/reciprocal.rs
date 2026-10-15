//! Material-bound strings in the existing two-way modal/friction network.
//!
//! The specimen spans TWO fixed end supports. An interior connection can
//! represent a seated bridge with an afterlength or another bilateral support.
//! Pinned shapes vanish at endpoints, so an endpoint cannot be a moving port.
//! The same signed displacement shape distributes reactions back to both
//! string and body. The body is not a one-way sound filter.
//!
//! Dispersion, material bending loss and cylinder-gas drag use the existing
//! material bowed-card lowering. Single strings and sympathetic banks share
//! one assembly and the existing coupled-modal/friction owners. No new string,
//! friction or integration law. This is fixed-state small-amplitude bending
//! at a fixed bow station, not torsion, nonlinear tension, moving contact,
//! broadband exterior back-loading or experimentally validated instrument timbre.

use crate::acoustic_realize::AcousticRealizeError;
use crate::bowed_string::BowedStringCard;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::{GestureCompileError, ScheduledRenderer};
use crate::render::schedule::force::coupled::{
    ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionConfig, ModalFrictionError,
};
use crate::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use crate::string_specimen::ResolvedStringSpecimen;
use fs_exec::{CancelGate, Cx};
use fs_material::gas::GasState;
use fs_scenario::{AmbientGas, gesture::GestureSchedule};
use fs_tribo::FrictionLaw;

/// A bilateral connection at an INTERIOR point on the full string.
/// Body shapes must use the actual mass-normalized displacement basis;
/// their source and geometry remain the body producer's responsibility.
#[derive(Clone, Debug, PartialEq)]
pub struct InteriorStringConnection {
    /// Fraction of the complete supported span, strictly inside (0,1).
    pub station_fraction: f64,
    /// Index into the supplied body list, before strings are prepended.
    pub body: ModalAttachment,
    /// Nonnegative transverse stiffness [N/m], independent of bulk string E.
    pub stiffness_n_m: f64,
    /// Nonnegative transverse viscous resistance [N s/m].
    pub damping_n_s_m: f64,
    /// Stress-free string-minus-body displacement [m].
    pub rest_extension_m: f64,
}

/// Explicit resolution, geometry and original physics-owner limits.
#[derive(Clone, Debug)]
pub struct ReciprocalStringOptions {
    /// Shared mechanical/audio rate; every body must already use this rate.
    pub sample_rate_hz: u32,
    /// Fixed bow station on the complete string, strictly inside (0,1).
    pub bow_station_fraction: f64,
    /// At least one nonzero physical connection to a supplied body.
    pub connections: Vec<InteriorStringConnection>,
    /// Limits for the material-derived string, without fake pressure.
    pub string_budget: ModalAcousticTimeBudget,
    /// Complete network's setup, coupled bandwidth, force and energy limits.
    pub network: ModalCouplingConfig,
    /// Original implicit friction's physical and nonlinear-solve limits.
    pub friction: ModalFrictionConfig,
}

/// Original material, environment and connection inputs. Numerical receiver
/// coefficients and supplied friction laws gain no new source authority.
#[derive(Clone, Debug)]
pub struct ReciprocalStringSource {
    specimen: ResolvedStringSpecimen,
    ambient: Option<AmbientGas>,
    gas: GasState,
    options: ReciprocalStringOptions,
    body_connections: Vec<ModalConnection>,
    material_plate: Option<plate::MaterialPlateBody>,
    material_cavity: Option<plate::MaterialBodyCavity>,
}
impl ReciprocalStringSource {
    /// Complete specimen, including material identities and receipts.
    #[must_use]
    pub const fn specimen(&self) -> &ResolvedStringSpecimen { &self.specimen }
    /// Original air input, only for legacy air constructors. An explicit gas
    /// has no fabricated terrestrial-air or humidity description.
    #[must_use]
    pub const fn ambient(&self) -> Option<&AmbientGas> { self.ambient.as_ref() }
    /// Actual admitted medium. Gas source receipts remain with its producer;
    /// numerical coherence is not measured-composition or gas-phase validity.
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.gas }
    /// Original geometry, clock and physical/numerical limits.
    #[must_use]
    pub const fn options(&self) -> &ReciprocalStringOptions { &self.options }
    /// Receiver links in the ORIGINAL supplied-body index basis.
    #[must_use]
    pub fn body_connections(&self) -> &[ModalConnection] { &self.body_connections }
    /// Original material body and observer inputs, when constructed from a plate.
    #[must_use]
    pub const fn material_plate(&self) -> Option<&plate::MaterialPlateBody> { self.material_plate.as_ref() }
    /// Original lumped cavity and neck observation, when present.
    #[must_use]
    pub const fn material_cavity(&self) -> Option<&plate::MaterialBodyCavity> { self.material_cavity.as_ref() }
}

/// Admission failures retain their original physical owner.
#[derive(Debug)]
pub enum ReciprocalStringError {
    /// Material, geometry or loss lowering refused the specimen.
    Material(AcousticRealizeError),
    /// Bilateral network or cancellation refused admission.
    Network(ModalCouplingError),
    /// Original implicit friction owner refused admission.
    Friction(ModalFrictionError),
}
impl core::fmt::Display for ReciprocalStringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "material string: {e:?}"),
            Self::Network(e) => write!(f, "string/body network: {e}"),
            Self::Friction(e) => write!(f, "string/body friction: {e}"),
        }
    }
}
impl std::error::Error for ReciprocalStringError {}
fn invalid(what: &'static str) -> ReciprocalStringError {
    ReciprocalStringError::Network(ModalCouplingError::Invalid { what })
}
fn poll(gate: &CancelGate) -> Result<(), ReciprocalStringError> {
    if gate.is_requested() { Err(ReciprocalStringError::Network(ModalCouplingError::Cancelled)) }
    else { Ok(()) }
}
fn interior(station: f64) -> bool { station.is_finite() && station > 0.0 && station < 1.0 }

/// Material string in a complete reciprocal network, not another bow stepper.
/// Material changes rebuild both modes and physical force/displacement ports.
pub struct ReciprocalMaterialString {
    source: ReciprocalStringSource,
    system: FrictionModalSystem,
}
impl ReciprocalMaterialString {
    /// Legacy air convenience form of `from_specimen_with_gas`.
    /// Body vibration and observations survive. A bridge must be interior to
    /// the complete span INCLUDING afterlength, never a fake moving pinned end.
    #[allow(clippy::too_many_arguments)]
    pub fn from_specimen(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, ambient: AmbientGas,
        bodies: Vec<ModalAcousticTimeModel>, law: FrictionLaw,
        options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let gas = crate::acoustic_realize::gas_state(ambient).map_err(ReciprocalStringError::Material)?;
        let mut model = Self::from_specimen_with_gas(cx, specimen, gas, bodies, law, options, gate)?;
        model.source.ambient = Some(ambient);
        Ok(model)
    }

    /// Same mechanics in an explicit coherent gas/mixture. Caller-supplied body
    /// modes must already represent that gas; source-built plate/cavity paths
    /// derive them together. No inferred air absorption or live medium change.
    #[allow(clippy::too_many_arguments)]
    pub fn from_specimen_with_gas(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, gas: GasState,
        bodies: Vec<ModalAcousticTimeModel>, law: FrictionLaw,
        options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::from_specimen_and_body_network(cx, specimen, gas, bodies, Vec::new(), law, options, gate)
    }

    /// Include receiver-to-receiver physics in the SAME bow and bridge solve.
    /// Extra links use original supplied-body indices, remapped once by the
    /// common string-bank assembler. Their storage, reactions and losses are
    /// not applied as delayed or one-way body drives.
    ///
    /// # Errors
    /// Missing bodies/bridges, invalid ports or links, incompatible clocks/bases,
    /// material/gas/loss limits, complete-network budgets or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_specimen_and_body_network(
        cx: &Cx<'_>, specimen: ResolvedStringSpecimen, gas: GasState,
        bodies: Vec<ModalAcousticTimeModel>, body_connections: Vec<ModalConnection>,
        law: FrictionLaw, options: ReciprocalStringOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if bodies.is_empty() || options.connections.is_empty() || !interior(options.bow_station_fraction) {
            return Err(invalid("reciprocal string requires bodies, connections and an interior bow station"));
        }
        let total = bodies.iter().try_fold(specimen.string().n_modes,
            |sum, body| sum.checked_add(body.modes().len()))
            .ok_or_else(|| invalid("reciprocal string mode count overflow"))?;
        let link_count = options.connections.len().checked_add(body_connections.len())
            .ok_or_else(|| invalid("reciprocal string connection count overflow"))?;
        if total > options.network.max_modes || link_count > options.network.max_connections {
            return Err(invalid("reciprocal string exceeds its declared setup dimensions"));
        }
        let card = BowedStringCard::from_material_specimen_with_gas(cx, &specimen, &gas, options.sample_rate_hz)
            .map_err(ReciprocalStringError::Material)?;
        let shared = bank::StringBankOptions {
            sample_rate_hz: options.sample_rate_hz,
            bow: bank::StringBowLocation { string: 0, station_fraction: options.bow_station_fraction },
            string_budget: options.string_budget, network: options.network, friction: options.friction,
        };
        let system = bank::assemble(&[card], bodies, &[&options.connections], &body_connections,
            law, &shared, gate)?;
        poll(gate)?;
        Ok(Self { source: ReciprocalStringSource { specimen, ambient: None, gas, options,
            body_connections, material_plate: None, material_cavity: None }, system })
    }

    /// Original source binding.
    #[must_use]
    pub const fn source(&self) -> &ReciprocalStringSource { &self.source }
    /// Accepted state and original work/loss diagnostics.
    #[must_use]
    pub const fn system(&self) -> &FrictionModalSystem { &self.system }
    /// Advance using the original owner's transactional step or gated step.
    pub fn system_mut(&mut self) -> &mut FrictionModalSystem { &mut self.system }
    /// Transfer both accepted mechanics and its original sources.
    #[must_use]
    pub fn into_parts(self) -> (ReciprocalStringSource, FrictionModalSystem) { (self.source, self.system) }

    /// Bind a physical bow track without dropping material or receiver sources.
    /// Only scheduling changes; the existing renderer owns every physical step.
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

/// Material/geometry paired with the common sample-accurate renderer.
pub struct ReciprocalStringPerformance {
    source: ReciprocalStringSource,
    renderer: ScheduledRenderer,
}
impl ReciprocalStringPerformance {
    /// Exact source binding used to build this performance.
    #[must_use]
    pub const fn source(&self) -> &ReciprocalStringSource { &self.source }
    /// Read-only common renderer.
    #[must_use]
    pub const fn renderer(&self) -> &ScheduledRenderer { &self.renderer }
    /// Render with existing control, callback and cancellation semantics.
    pub fn renderer_mut(&mut self) -> &mut ScheduledRenderer { &mut self.renderer }
}

#[cfg(test)]
mod tests;
/// Source-derived plate and optional cavity in the same reciprocal network.
pub mod plate;
/// Shared-body sympathetic strings with heterogeneous materials and explicit gas.
pub mod bank;
