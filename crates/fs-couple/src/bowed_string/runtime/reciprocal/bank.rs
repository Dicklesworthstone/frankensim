//! Several material strings on shared, reciprocally loaded bodies.
//!
//! Unplayed strings are mechanical components, not muted voices. Their motion
//! loads the same receiver and returns through its original bilateral solve to
//! every other string. One prescribed bow drives one selected material string;
//! the existing friction owner accounts for its work and the whole network's
//! loss. No independent voice mixing, per-note normalization or new stepper.
//!
//! All spans include both fixed endpoints and any afterlength. Body attachment
//! stations are interior points, not fictitious moving ports on pinned ends.
//! Sources remain frozen-state linear bending models; concentric layers mean
//! perfect bonding, not helical winding/slip. Optional sourced torsion shares
//! the surface friction solve; large-amplitude tension and broadband radiation
//! back-loading are absent. Multiple fixed bow contacts
//! use the existing simultaneous-friction owner through `multiple`.

use super::{InteriorStringConnection, ReciprocalStringError, invalid, interior, poll};
use crate::bowed_string::{BowedStringCard, unit_shape};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::{GestureCompileError, ScheduledRenderer};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalFrame, FrictionModalSystem, ModalFrictionConfig,
    ModalFrictionError, ModalFrictionPort,
};
use crate::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use crate::string_specimen::{ResolvedStringSpecimen, StringPrestress};
use fs_exec::{CancelGate, Cx};
use fs_material::gas::GasState;
use fs_material::section::ConcentricKelvinVoigtSection;
use fs_scenario::{PrestressedString, gesture::{GestureSchedule, GestureTarget}};
use fs_tribo::FrictionLaw;

/// Original source data, not an averaged material or an authored modal preset.
#[derive(Clone, Debug)]
pub enum MaterialStringSource {
    /// Original bending source extended by material-derived small-angle twist.
    WithTorsion(Box<torsion::TorsionalStringSource>),
    /// Homogeneous circular specimen, with its complete material/loss receipts.
    Homogeneous(ResolvedStringSpecimen),
    /// Perfectly bonded annuli with per-layer elasticity, density and viscosity.
    Concentric {
        /// Span and modal resolution. Inherited internal damping must be absent.
        template: PrestressedString,
        /// Geometry, every source state and the intersected loss applicability.
        section: ConcentricKelvinVoigtSection,
        /// Independent support prescription; no hidden pitch normalization.
        prestress: StringPrestress,
    },
}
impl MaterialStringSource {
    fn modes(&self) -> usize {
        match self { Self::WithTorsion(t) => t.total_modes(),
            Self::Homogeneous(s) => s.string().n_modes,
            Self::Concentric { template, .. } => template.n_modes }
    }
    fn lower(&self, cx: &Cx<'_>, gas: &GasState, rate: u32)
        -> Result<BowedStringCard, ReciprocalStringError>
    {
        match self {
            Self::WithTorsion(t) => return t.bending().lower(cx, gas, rate),
            Self::Homogeneous(s) => BowedStringCard::from_material_specimen_with_gas(cx, s, gas, rate),
            Self::Concentric { template, section, prestress } =>
                BowedStringCard::from_concentric_material_section(cx, template.clone(), section, *prestress, gas, rate),
        }.map_err(ReciprocalStringError::Material)
    }
}

/// One full supported span and its independently specified body connections.
#[derive(Clone, Debug)]
pub struct MaterialStringMember {
    /// Source-bound solid geometry/materials and support prescription.
    pub source: MaterialStringSource,
    /// Body indices are relative to the supplied body list, not network order.
    /// Every member needs at least one nonzero connection to that shared system.
    pub connections: Vec<InteriorStringConnection>,
}

/// Actual bow contact, not an instrument or note identifier.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringBowLocation {
    /// Zero-based index in the original material string list.
    pub string: usize,
    /// Interior station along that member's full supported span.
    pub station_fraction: f64,
}

/// Shared clock and original physics-owner budgets.
#[derive(Clone, Debug)]
pub struct StringBankOptions {
    /// Every body and string must use this physical rate [Hz].
    pub sample_rate_hz: u32,
    /// One selected material string and its fixed bow station.
    pub bow: StringBowLocation,
    /// Original single-component modal admission.
    pub string_budget: ModalAcousticTimeBudget,
    /// Original complete-network dimension, bandwidth, force and energy limits.
    pub network: ModalCouplingConfig,
    /// Original implicit friction's constitutive/solve limits.
    pub friction: ModalFrictionConfig,
}

/// Original frozen material inputs; no new authority for supplied body modes.
#[derive(Clone, Debug)]
pub struct MaterialStringBankSource {
    strings: Vec<MaterialStringMember>,
    gas: GasState,
    options: StringBankOptions,
    body_connections: Vec<ModalConnection>,
    material_plate: Option<plate::SharedMaterialPlateBody>,
    material_chart: Option<chart::SharedMaterialPlateChartBody>,
    material_chart_viscosity: Option<crate::plate_viscosity::ViscousPlateImage>,
    material_chart_bridges: Option<chart::PlateBridgeFootprints>,
    material_cavity: Option<super::plate::MaterialBodyCavity>,
}
impl MaterialStringBankSource {
    /// Original string order, material/loss receipts and physical connections.
    #[must_use]
    pub fn strings(&self) -> &[MaterialStringMember] { &self.strings }
    /// Same coherent gas used to lower all string loss models.
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.gas }
    /// Original selected bow, shared clock and explicit model budgets.
    #[must_use]
    pub const fn options(&self) -> &StringBankOptions { &self.options }
    /// Original receiver-to-receiver links before string-prefix remapping.
    #[must_use]
    pub fn body_connections(&self) -> &[ModalConnection] { &self.body_connections }
    /// Original source-bound plate and complete physical bridge/observer inputs.
    #[must_use]
    pub const fn material_plate(&self) -> Option<&plate::SharedMaterialPlateBody> { self.material_plate.as_ref() }
    /// Original regional geometry, material receipts and physical force footprint.
    #[must_use]
    pub const fn material_chart(&self) -> Option<&chart::SharedMaterialPlateChartBody> { self.material_chart.as_ref() }
    /// Full tensor material loss and its SAME-basis projection, when selected.
    /// The original viscous sources survive single/multiple-bow scheduling.
    #[must_use]
    pub const fn material_chart_viscosity(&self) -> Option<&crate::plate_viscosity::ViscousPlateImage> {
        self.material_chart_viscosity.as_ref()
    }
    /// Original body-air law, geometry, gas and observer, when supplied.
    #[must_use]
    pub const fn material_cavity(&self) -> Option<&super::plate::MaterialBodyCavity> { self.material_cavity.as_ref() }
}

/// All string/body vibration persists in ONE mechanical/friction state.
/// Components are ordered as all strings, followed by all supplied bodies.
pub struct ReciprocalStringBank {
    source: MaterialStringBankSource,
    cards: Vec<BowedStringCard>,
    system: FrictionModalSystem,
}
impl ReciprocalStringBank {
    /// Compile homogeneous and layered strings into the same reciprocal system.
    ///
    /// Strings begin at rest; existing body modal states are preserved. The
    /// selected bow's traction is solved against the complete coupled response,
    /// including unplayed strings. External forcing can excite any component.
    /// Supplied body pressure observations and material validity remain their
    /// producer's responsibility; no velocity is relabeled as pressure here.
    ///
    /// # Errors
    /// Missing or detached members, bad indices/stations, incompatible clocks,
    /// material/gas/loss-domain refusal, network limits or cancellation.
    pub fn new(
        cx: &Cx<'_>, strings: Vec<MaterialStringMember>, gas: GasState,
        bodies: Vec<ModalAcousticTimeModel>, law: FrictionLaw,
        options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        Self::with_body_connections(cx, strings, gas, bodies, Vec::new(), law, options, gate)
    }

    /// Retain connected receiver physics such as a plate and shared cavity neck.
    /// Extra links address the ORIGINAL body list, not the string-prefixed bank.
    /// Their storage, force reactions and loss participate in the same bow solve
    /// as every played and sympathetic string. Body source metadata remains the
    /// caller's responsibility; a raw connection cannot certify its material law.
    ///
    /// # Errors
    /// Same source admission as `new`, plus invalid body indices/links and the
    /// total budget including every receiver connection. No partial bank escapes.
    #[allow(clippy::too_many_arguments)]
    pub fn with_body_connections(
        cx: &Cx<'_>, strings: Vec<MaterialStringMember>, gas: GasState,
        bodies: Vec<ModalAcousticTimeModel>, body_connections: Vec<ModalConnection>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        preflight(&strings, &bodies, &options, body_connections.len())?;
        let mut cards = Vec::with_capacity(strings.len());
        for member in &strings {
            poll(gate)?;
            cards.push(member.source.lower(cx, &gas, options.sample_rate_hz)?);
        }
        let connections: Vec<&[InteriorStringConnection]> = strings.iter()
            .map(|member| member.connections.as_slice()).collect();
        let system = assemble_sources(&cards, Some(&strings), bodies, &connections, &body_connections, law, &options, gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        poll(gate)?;
        Ok(Self { source: MaterialStringBankSource { strings, gas, options, body_connections, material_plate: None, material_chart: None, material_chart_viscosity: None, material_chart_bridges: None, material_cavity: None }, cards, system })
    }

    /// Complete source assignments survive rendering and sympathetic vibration.
    #[must_use]
    pub const fn source(&self) -> &MaterialStringBankSource { &self.source }
    /// Derived TRANSVERSE coefficients; optional torsion and receipts remain in `source()`.
    #[must_use]
    pub fn cards(&self) -> &[BowedStringCard] { &self.cards }
    /// Accepted mechanics and work/loss diagnostics, without mutable source drift.
    #[must_use]
    pub const fn system(&self) -> &FrictionModalSystem { &self.system }
    /// Advance all members using the original transactional physical step.
    pub fn step(&mut self, external: &[f64], drive: FrictionDrive)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.system.step(external, drive) }
    /// Same step with the original cancellation semantics.
    pub fn step_under_gate(&mut self, external: &[f64], drive: FrictionDrive, gate: &CancelGate)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.system.step_under_gate(external, drive, gate) }
    /// Explicitly transfer the complete state and its original source bindings.
    #[must_use]
    pub fn into_parts(self) -> (MaterialStringBankSource, FrictionModalSystem) { (self.source, self.system) }

    /// Bind an existing bow track without dropping unplayed strings or body state.
    /// Its string index must name the ACTUAL selected member; a matching station
    /// on a different string is not the same physical contact.
    pub fn into_scheduled(self, schedule: &GestureSchedule, track: &str,
        config: FrictionGestureConfig) -> Result<MaterialStringBankPerformance, GestureCompileError>
    {
        let bow = self.source.options.bow;
        let target = schedule.tracks().iter().find(|t| t.id == track).map(|t| &t.target);
        if !matches!(target, Some(GestureTarget::BowStroke { string }) if usize::from(*string) == bow.string)
            || config.station_fraction.to_bits() != bow.station_fraction.to_bits()
            || config.sample_rate_hz != self.source.options.sample_rate_hz {
            return Err(GestureCompileError::Invalid {
                what: "string-bank track must match its actual string, fixed bow station and clock",
            });
        }
        let external = vec![0.0; self.system.mode_count()];
        let renderer = ScheduledRenderer::from_friction_gesture(self.system, external, schedule, track, config)?;
        Ok(MaterialStringBankPerformance { source: self.source, renderer })
    }
}

/// Source-retaining common renderer; no independent per-string render loops.
pub struct MaterialStringBankPerformance {
    source: MaterialStringBankSource,
    renderer: ScheduledRenderer,
}
impl MaterialStringBankPerformance {
    /// Original materials, geometry, connections and common gas.
    #[must_use]
    pub const fn source(&self) -> &MaterialStringBankSource { &self.source }
    /// Original scheduled renderer and accepted clock.
    #[must_use]
    pub const fn renderer(&self) -> &ScheduledRenderer { &self.renderer }
    /// Render blocks using existing scheduling, cancellation and refusal rules.
    pub fn renderer_mut(&mut self) -> &mut ScheduledRenderer { &mut self.renderer }
}

fn preflight(strings: &[MaterialStringMember], bodies: &[ModalAcousticTimeModel],
    options: &StringBankOptions, body_links: usize) -> Result<(), ReciprocalStringError>
{
    if strings.is_empty() || bodies.is_empty() || strings.len() > options.network.max_modes
        || bodies.len() > options.network.max_modes || options.bow.string >= strings.len()
        || !interior(options.bow.station_fraction) || options.sample_rate_hz == 0 {
        return Err(invalid("string bank requires bounded members, bodies, a clock and an actual interior bow"));
    }
    let mut modes = 0usize;
    let mut links = body_links;
    for member in strings {
        if member.source.modes() == 0 || member.connections.is_empty() {
            return Err(invalid("every material string requires modes and a shared-body connection"));
        }
        modes = modes.checked_add(member.source.modes()).ok_or_else(|| invalid("string bank mode count overflow"))?;
        links = links.checked_add(member.connections.len()).ok_or_else(|| invalid("string bank connection count overflow"))?;
    }
    for body in bodies {
        modes = modes.checked_add(body.modes().len()).ok_or_else(|| invalid("string bank mode count overflow"))?;
    }
    if modes > options.network.max_modes || links > options.network.max_connections {
        return Err(invalid("string bank exceeds declared complete-network dimensions"));
    }
    Ok(())
}

// One assembly for single/multiple strings, including receiver-to-receiver links.
// All force condensation, stepping, friction and energy accounting remain in the
// original CoupledModalSystem/FrictionModalSystem. Body links are remapped once.
pub(super) fn assemble(
    cards: &[BowedStringCard], bodies: Vec<ModalAcousticTimeModel>,
    connections: &[&[InteriorStringConnection]], body_connections: &[ModalConnection],
    law: FrictionLaw, options: &StringBankOptions, gate: &CancelGate,
) -> Result<FrictionModalSystem, ReciprocalStringError> {
    assemble_sources(cards, None, bodies, connections, body_connections, law, options, gate)
}

// Keep the original one-string API in the same assembler. Its absent source
// extensions retain the original basis and arithmetic exactly.
#[allow(clippy::too_many_arguments)]
fn assemble_sources(
    cards: &[BowedStringCard], sources: Option<&[MaterialStringMember]>,
    bodies: Vec<ModalAcousticTimeModel>, connections: &[&[InteriorStringConnection]],
    body_connections: &[ModalConnection], law: FrictionLaw,
    options: &StringBankOptions, gate: &CancelGate,
) -> Result<FrictionModalSystem, ReciprocalStringError> {
    poll(gate)?;
    if cards.is_empty() || cards.len() != connections.len() || bodies.is_empty()
        || sources.is_some_and(|s| s.len() != cards.len())
        || options.bow.string >= cards.len() || !interior(options.bow.station_fraction) {
        return Err(invalid("reciprocal assembly requires aligned strings, bodies and an interior bow"));
    }
    let mut components = Vec::with_capacity(cards.len() + bodies.len());
    let mut links = Vec::new();
    for (index, (card, connections)) in cards.iter().zip(connections).enumerate() {
        poll(gate)?;
        let source = sources.map(|members| &members[index].source);
        let total_modes = source.map_or(card.mode_count, MaterialStringSource::modes);
        for connection in *connections {
            let body = bodies.get(connection.body.component)
                .ok_or_else(|| invalid("string connection names an unknown supplied body"))?;
            if !interior(connection.station_fraction) || connection.body.shapes.len() != body.modes().len()
                || !connection.body.shapes.iter().all(|x| x.is_finite())
                || !connection.body.shapes.iter().any(|x| *x != 0.0)
                || !(connection.stiffness_n_m > 0.0 || connection.damping_n_s_m > 0.0) {
                return Err(invalid("string/body connection needs an interior station and nonzero physical ports"));
            }
            // This support transmits transverse motion, not an invented torque.
            let mut bridge_shape = shape(card, connection.station_fraction);
            bridge_shape.resize(total_modes, 0.0);
            links.push(ModalConnection {
                left: ModalAttachment { component: index, shapes: bridge_shape },
                right: ModalAttachment { component: cards.len() + connection.body.component,
                    shapes: connection.body.shapes.clone() },
                stiffness_n_m: connection.stiffness_n_m, damping_n_s_m: connection.damping_n_s_m,
                rest_extension_m: connection.rest_extension_m,
            });
        }
        let mut modes: Vec<_> = (0..card.mode_count).map(|k| ModalAcousticMode {
            angular_frequency_rad_s: card.mode_omega_rad_s(k+1), damping_ratio: card.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        }).collect();
        if let Some(MaterialStringSource::WithTorsion(twist)) = source {
            twist.append_modes(card, &mut modes, gate)?;
        }
        components.push(ModalAcousticTimeModel::try_new(options.sample_rate_hz, modes, options.string_budget)
            .map_err(|source| ReciprocalStringError::Network(ModalCouplingError::Component { component: index, source }))?);
    }
    for connection in body_connections {
        poll(gate)?;
        if connection.left.component >= bodies.len() || connection.right.component >= bodies.len() {
            return Err(invalid("receiver connection names an unknown original body component"));
        }
        let mut link = connection.clone();
        link.left.component += cards.len();
        link.right.component += cards.len();
        links.push(link);
    }
    components.extend(bodies);
    let network = CoupledModalSystem::new(components, links, options.network, gate)
        .map_err(ReciprocalStringError::Network)?;
    let bow_card = &cards[options.bow.string];
    let bow_shape = match sources {
        Some(members) => members[options.bow.string].source.surface_shape(bow_card, options.bow.station_fraction)?,
        None => shape(bow_card, options.bow.station_fraction),
    };
    let port = ModalFrictionPort { left: ModalAttachment { component: options.bow.string,
        shapes: bow_shape }, right: None, law };
    FrictionModalSystem::new(network, port, options.friction, gate).map_err(ReciprocalStringError::Friction)
}

pub(super) fn shape(card: &BowedStringCard, station: f64) -> Vec<f64> {
    (0..card.mode_count).map(|k| unit_shape(k, station, card.length_m, card.linear_density_kg_m)).collect()
}

/// One material-derived body and physical pressure shared by every string.
pub mod plate;

/// Regional material geometry in the same reciprocal string/body network.
pub mod chart;

/// Several source-bound bow contacts without resetting the shared mechanics.
pub mod multiple;

/// Material-derived surface twisting in the same reciprocal contact solve.
pub mod torsion;

/// The same material sources in the existing nonlinear impact/contact owner.
pub mod nonlinear;
