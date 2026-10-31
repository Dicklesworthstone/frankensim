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
//! the surface friction solve. Source-derived geometric stretching or retained
//! longitudinal waves can join that same network through the nonlinear module.
//! Explicit passive exterior radiation reacts through that same physical body.
//! Interior rotational connections can also
//! transfer material-string torque to the shared body. Multiple fixed bow contacts
//! use the existing simultaneous-friction owner through `multiple`.

use super::{InteriorStringConnection, ReciprocalStringError, invalid, interior, poll};
use crate::bowed_string::{BowedStringCard, unit_shape};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::{GestureCompileError, ScheduledRenderer};
use crate::render::schedule::force::coupled::{
    ConstantLossConnection, CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
    ModalCouplingError, ModalPreparation, PreparedModalBudget,
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
use fs_scenario::{PrestressedString, gesture::GestureSchedule};
use fs_tribo::FrictionLaw;

/// Original source data, not an averaged material or an authored modal preset.
#[derive(Clone, Debug)]
pub enum MaterialStringSource {
    /// Both physical transverse planes of the same source, with explicit directions.
    WithPolarizations(Box<polarization::PolarizedStringSource>),
    /// Original bending source extended by material-derived small-angle twist.
    WithTorsion(Box<torsion::TorsionalStringSource>),
    /// Homogeneous circular specimen, with its complete material/loss receipts.
    Homogeneous(ResolvedStringSpecimen),
    /// Direct measured whole-string mass, EI, diameter, length and tension.
    /// No intrinsic loss or hidden core, axial, shear or surface law is inferred.
    Measured(crate::material_assembly::measured_string::ResolvedMeasuredStringSpecimen),
    /// Sliding close-coiled winding: actual wire mass/spring stiffness and core flexure.
    Wound(crate::material_assembly::wound::ResolvedWoundStringSpecimen),
    /// Perfectly bonded annuli with per-layer elasticity, density and viscosity.
    Concentric {
        /// Span and modal resolution. Inherited internal damping must be absent.
        template: PrestressedString,
        /// Geometry, every source state and the intersected loss applicability.
        section: ConcentricKelvinVoigtSection,
        /// Independent support prescription; no hidden pitch normalization.
        prestress: StringPrestress,
    },
    /// Perfectly bonded annuli with exclusive per-layer bending memory.
    /// Every source arm survives in the shared exact mechanical owner.
    ConcentricProny {
        /// Span and modal resolution, without inherited internal damping.
        template: PrestressedString,
        /// Original layer states, geometry, selectors and independent spectra.
        section: fs_material::section::relaxation::ConcentricPronySection,
        /// Applied force or fixed extension, with no implicit pitch tuning.
        prestress: StringPrestress,
    },
}
impl MaterialStringSource {
    // Reconstruct the same transverse source when another medium supplies
    // inertia. Keep every intrinsic constitutive law instead of recovering it
    // by subtracting a sampled air loss from a rounded modal damping ratio.
    fn resolved_transverse_string(&self) -> Result<PrestressedString, ReciprocalStringError> {
        match self {
            Self::WithPolarizations(p) => p.source().resolved_transverse_string(),
            Self::WithTorsion(t) => t.bending().resolved_transverse_string(),
            Self::Homogeneous(s) => Ok(s.string()),
            Self::Measured(s) => Ok(s.string()),
            Self::Wound(s) => Ok(s.string()),
            Self::Concentric { template, section, prestress } =>
                crate::material_assembly::concentric::with_concentric_kelvin_voigt_section(
                    template.clone(), section, *prestress,
                ).map_err(ReciprocalStringError::Material),
            Self::ConcentricProny { template, section, prestress } =>
                crate::material_assembly::concentric::with_concentric_prony_section(
                    template.clone(), section, *prestress,
                ).map_err(ReciprocalStringError::Material),
        }
    }
    fn modes(&self) -> usize {
        match self { Self::WithPolarizations(p) => p.total_modes(),
            Self::WithTorsion(t) => t.total_modes(),
            Self::Homogeneous(s) => s.string().n_modes,
            Self::Measured(s) => s.string().n_modes,
            Self::Wound(s) => s.string().n_modes,
            Self::Concentric { template, .. } | Self::ConcentricProny { template, .. } => template.n_modes }
    }
    fn lower(&self, cx: &Cx<'_>, gas: &GasState, rate: u32)
        -> Result<BowedStringCard, ReciprocalStringError>
    {
        if self.torsion_source().is_some_and(|twist| twist.relaxation().is_some()) {
            return Err(invalid("this scalar string-card consumer cannot retain torsional material history"));
        }
        match self {
            Self::WithPolarizations(p) => return p.source().lower(cx, gas, rate),
            Self::WithTorsion(t) => return t.bending().lower(cx, gas, rate),
            Self::Homogeneous(s) => BowedStringCard::from_material_specimen_with_gas(cx, s, gas, rate),
            Self::Measured(s) => BowedStringCard::from_measured_string_specimen(cx, s, gas, rate),
            Self::Wound(s) => BowedStringCard::from_wound_material_specimen(cx, s, gas, rate),
            Self::Concentric { template, section, prestress } =>
                BowedStringCard::from_concentric_material_section(cx, template.clone(), section, *prestress, gas, rate),
            Self::ConcentricProny { .. } =>
                return Err(invalid("this scalar string-card consumer cannot retain layered bending history")),
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
    material_bridge_interfaces: Vec<interfaces::MaterialBridgeInterface>,
    translational_body_connections: Vec<polarization::StringBodyTranslationConnection>,
    rotational_body_connections: Vec<torsion::StringBodyRotationConnection>,
    fluid_loading: Option<fluid::StringBankFluidLoading>,
    geometric_stretching: Option<nonlinear::bank::StringBankStretching>,
    longitudinal_waves: Option<nonlinear::longitudinal::StringBankLongitudinal>,
    prepared_linear: ModalPreparation,
    material_plate: Option<plate::SharedMaterialPlateBody>,
    material_chart: Option<chart::SharedMaterialPlateChartBody>,
    material_chart_basis: Option<crate::thin_plate::RetainedPlateBasis>,
    material_chart_viscosity: Option<crate::plate_viscosity::ViscousPlateImage>,
    material_chart_relaxation: Option<crate::plate_viscosity::relaxation::ModalPlateRelaxation>,
    material_chart_tensor_relaxation: Option<crate::plate_viscosity::relaxation::tensor::TensorPlateRelaxation>,
    string_bending_relaxation: Option<relaxation::StringBankBendingRelaxation>,
    material_chart_bridges: Option<chart::PlateBridgeFootprints>,
    material_shell: Option<shell::SharedMaterialShellBody>,
    material_shell_image: Option<crate::shear_shell::MaterialShellImage>,
    material_shell_assembly: Option<shell::SharedMaterialShellAssemblyBody>,
    material_shell_assembly_image: Option<crate::shear_shell::assembly::MaterialShellAssemblyImage>,
    material_shell_coating_relaxation: Option<crate::shear_shell::coating::relaxation::ShellCoatingRelaxation>,
    material_cavity: Option<super::plate::MaterialBodyCavity>,
    body_radiation: Option<radiation::MaterialBodyRadiation>,
    spatial_cavity: Option<cavity::MaterialShellCavity>,
    residual_flexibility: bool,
}
impl MaterialStringBankSource {
    /// Whether bridges and stops carry the omitted bending modes' static
    /// compliance (`ReciprocalStringBank::with_residual_flexibility`).
    #[must_use]
    pub const fn residual_flexibility(&self) -> bool { self.residual_flexibility }
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
    /// Independently directed centreline/body laws in the original supplied-body basis.
    #[must_use]
    pub fn translational_body_connections(&self) -> &[polarization::StringBodyTranslationConnection] {
        &self.translational_body_connections
    }
    /// Original torque/rotation laws and angular body maps, in supplied order.
    /// Their string maps use the same source-bound polar inertia as bow torsion.
    #[must_use]
    pub fn rotational_body_connections(&self) -> &[torsion::StringBodyRotationConnection] {
        &self.rotational_body_connections
    }
    /// Explicit complete-operator preparation budget, retained through fluid,
    /// torque, multiple-bow and scheduled-performance ownership.
    #[must_use]
    pub const fn prepared_linear_budget(&self) -> Option<PreparedModalBudget> { self.prepared_linear.prepared_budget() }
    /// The explicit construction purpose survives material, gas and connection rebuilds.
    #[must_use]
    pub const fn preparation(&self) -> ModalPreparation { self.prepared_linear }
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
    /// Complete material-derived curved shear shell and physical bridge geometry.
    #[must_use]
    pub const fn material_shell(&self) -> Option<&shell::SharedMaterialShellBody> {
        self.material_shell.as_ref()
    }
    /// One common shell basis for every force, full-loss and volume projection.
    #[must_use]
    pub const fn material_shell_image(&self) -> Option<&crate::shear_shell::MaterialShellImage> {
        self.material_shell_image.as_ref()
    }
    /// Complete material parts, ties and physical ports before body reduction.
    #[must_use]
    pub const fn material_shell_assembly(&self) -> Option<&shell::SharedMaterialShellAssemblyBody> {
        self.material_shell_assembly.as_ref()
    }
    /// Shared assembled basis for all strings, material losses and fluid rows.
    #[must_use]
    pub const fn material_shell_assembly_image(&self)
        -> Option<&crate::shear_shell::assembly::MaterialShellAssemblyImage>
    {
        self.material_shell_assembly_image.as_ref()
    }
    /// Original body-air law, geometry, gas and observer, when supplied.
    #[must_use]
    pub const fn material_cavity(&self) -> Option<&super::plate::MaterialBodyCavity> { self.material_cavity.as_ref() }
}

/// All string/body vibration persists in ONE mechanical/friction state.
/// Components are strings, supplied physical bodies, optional fluid memory,
/// then explicitly selected longitudinal coordinates.
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
    /// Original homogeneous Prony specimens keep one material history per
    /// positive arm and transverse mode in this same exact owner. Their
    /// selected bending loss is never collapsed to a modal decay ratio.
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
        Self::with_body_connections_prepared(cx, strings, gas, bodies, body_connections, None, law, options, gate)
    }

    /// `with_body_connections` with an explicit construction choice. An
    /// `Option<PreparedModalBudget>` preserves ordinary/prepared playback;
    /// `ModalPreparation::HarmonicOnly` retains the complete physical bank
    /// without a time propagator and refuses playback. The choice survives
    /// material, gas and connection rebuilding.
    /// A connected network above the ordinary dense block (42 coordinates) is refused by its own constructor, so the
    /// record has to be present at construction; `with_prepared_linear_propagation` afterwards cannot rescue it.
    /// The budget is kept in the bank's source for every later fluid/torsion rebuild, and must not be bound again.
    ///
    /// # Errors
    /// Those of `with_body_connections`, plus a prepared budget the complete network does not fit.
    #[allow(clippy::too_many_arguments)]
    pub fn with_body_connections_prepared(
        cx: &Cx<'_>, strings: Vec<MaterialStringMember>, gas: GasState,
        bodies: Vec<ModalAcousticTimeModel>, body_connections: Vec<ModalConnection>,
        prepared: impl Into<ModalPreparation>,
        law: FrictionLaw, options: StringBankOptions, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        preflight(&strings, &bodies, &options, body_connections.len())?;
        let prepared = prepared.into();
        let mut cards = Vec::with_capacity(strings.len());
        for member in &strings {
            poll(gate)?;
            cards.push(member.source.lower_for_bank(cx, &gas, options.sample_rate_hz)?);
        }
        let connections: Vec<&[InteriorStringConnection]> = strings.iter()
            .map(|member| member.connections.as_slice()).collect();
        let system = assemble_sources(&cards, Some(&strings), bodies, &connections, &body_connections, &[], &[], None, prepared, &[], false, law, &options, gate)?;
        cx.checkpoint().map_err(|_| ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Cancelled))?;
        poll(gate)?;
        let source = MaterialStringBankSource { strings, gas, options, body_connections, material_bridge_interfaces: Vec::new(), translational_body_connections: Vec::new(), rotational_body_connections: Vec::new(), fluid_loading: None, geometric_stretching: None, longitudinal_waves: None, prepared_linear: prepared, material_plate: None, material_chart: None, material_chart_basis: None, material_chart_viscosity: None, material_chart_relaxation: None, material_chart_tensor_relaxation: None, string_bending_relaxation: None, material_chart_bridges: None, material_shell: None, material_shell_image: None, material_shell_assembly: None, material_shell_assembly_image: None, material_shell_coating_relaxation: None, material_cavity: None, body_radiation: None, spatial_cavity: None, residual_flexibility: false };
        let system = source.restore_material_memory(system, &cards, None, gate)?;
        Ok(Self { source, cards, system })
    }

    /// Select bounded cold exact propagation before playback or before adding
    /// a larger material-gas memory network. The same complete coordinates,
    /// source receipts, physical limits, initial state, pressure and friction
    /// interface remain present; each later fluid/torque rebuild reuses this
    /// explicit budget. No ordinary constructor silently selects this path.
    ///
    /// The initial bank must first fit its ordinary constructor. For a bank
    /// whose added cylinder memory exceeds that capacity, call this method
    /// BEFORE with_cylinder_gas_loading. Larger byte/work budgets can refuse
    /// without substituting a reduced model.
    pub fn with_prepared_linear_propagation(
        mut self, budget: PreparedModalBudget, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        self.system = self.system.with_prepared_linear_propagation(budget, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.prepared_linear = ModalPreparation::Prepared(budget);
        poll(gate)?;
        Ok(self)
    }

    /// Put the omitted bending modes' static compliance in series with every
    /// string's bridge, and later with its stops (`residual`; pexfp.299).
    ///
    /// Each string's single spring/dashpot bridge connection is rebuilt as the
    /// retained modes see it: k/(1 + k G_bb), with the dashpot to first order.
    /// G_bb is the exact static Green's function of that string's own tension
    /// and EI minus its retained modes' static sum. The source keeps k and c;
    /// later fluid, torsion and cavity rebuilds keep the condensation. Stops
    /// bound later carry the cross terms. Directed bridges (polarization.rs)
    /// condense along their own axes, admitted only orthogonal to every other
    /// connection on their string. Fixed bows opt into the same correction with
    /// `MultiBowStringBank::with_bow_residual_flexibility` after their assembly.
    ///
    /// # Errors
    /// After motion, playback or bound nonlinear storage; repeated; beside
    /// rotational bridges, material bridge interfaces, a directed bridge with a
    /// nonzero cross term or a string with other than one bridge connection; or
    /// the original assembly's own refusals.
    pub fn with_residual_flexibility(mut self, gate: &CancelGate) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if self.system.has_nonlinear_potentials() || self.source.residual_flexibility
            || !self.source.material_bridge_interfaces.is_empty() {
            return Err(invalid("residual flexibility binds once, before nonlinear storage and without material bridge interfaces"));
        }
        if self.system.samples_rendered() != 0
            || self.system.components()[..self.cards.len()].iter().flat_map(|c| c.states())
                .any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("residual flexibility binds before string motion or playback"));
        }
        // The rotational rebuild's own recipe: same cards, bodies, links and memory.
        let body_end = self.source.fluid_loading.as_ref().map_or(
            self.system.components().len(), |fluid| self.cards.len() + fluid.body_components,
        );
        let bodies = self.system.components()[self.cards.len()..body_end].to_vec();
        let bridges: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(&self.cards, Some(&self.source.strings), bodies,
            &bridges, &self.source.body_connections, &self.source.translational_body_connections,
            &self.source.rotational_body_connections,
            self.source.fluid_loading.as_ref(), self.source.prepared_linear,
            &self.source.constant_loss_bridges(), true, self.system.port().law.clone(),
            &self.source.options, gate)?;
        if let Some(weights) = self.system.acceleration_pressure_weights() {
            system = system.with_acceleration_pressure(weights.to_vec())
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(binding) = self.system.interface() {
            system = system.with_interface(binding.clone()).map_err(ReciprocalStringError::Friction)?;
        }
        system = self.source.restore_material_memory(system, &self.cards,
            Some(self.system.material_memory().to_vec()), gate)?;
        self.source.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.system = system;
        self.source.residual_flexibility = true;
        Ok(self)
    }

    /// Complete source assignments survive rendering and sympathetic vibration.
    #[must_use]
    pub const fn source(&self) -> &MaterialStringBankSource { &self.source }
    /// Derived TRANSVERSE coefficients, including instantaneous added gas inertia
    /// when selected. Original solid cards and fluid coefficients remain in
    /// `source().cylinder_gas_loading()`; torsion retains its solid polar inertia.
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

    /// Bind an exact declared path to the actual source-selected material member.
    pub fn into_scheduled_path(
        self,
        source: &fs_scenario::gesture::ResolvedSchedule,
        path: &fs_scenario::assembly::ParameterPath,
        config: FrictionGestureConfig,
    ) -> Result<
        (
            MaterialStringBankPerformance,
            crate::gesture_path::BoundPaths,
        ),
        GestureCompileError,
    > {
        use crate::gesture_path::{PathBinding, PathInput, bind_paths};
        let bound = bind_paths(
            source,
            &[PathBinding {
                path: path.clone(),
                input: PathInput::Bow {
                    member: self.source.options.bow.string,
                },
            }],
        )
        .map_err(GestureCompileError::Gesture)?;
        let track = source
            .track(path)
            .map_err(GestureCompileError::Gesture)?
            .id
            .clone();
        let performance = self.into_scheduled_admitted(source.schedule(), &track, config,Some(path))?;
        Ok((performance, bound))
    }

    /// Bind an existing bow track without dropping unplayed strings or body state.
    /// Its string index must name the ACTUAL selected member; a matching station
    /// on a different string is not the same physical contact.
    pub fn into_scheduled(
        self,
        schedule: &GestureSchedule,
        track: &str,
        config: FrictionGestureConfig,
    ) -> Result<MaterialStringBankPerformance, GestureCompileError> {
        self.into_scheduled_admitted(schedule,track,config,None)
    }

    fn into_scheduled_admitted(self,schedule:&GestureSchedule,track:&str,
        config:FrictionGestureConfig,path:Option<&fs_scenario::assembly::ParameterPath>)
        ->Result<MaterialStringBankPerformance,GestureCompileError> {
        let bow = self.source.options.bow;
        let source_track = schedule.tracks().iter().find(|t| t.id == track);
        let accepts=match path {
            Some(path)=>source_track.is_some_and(|t|crate::gesture_path::matches_input(t,crate::gesture_path::PathInput::Bow {member:bow.string},Some(path))),
            None=>source_track.and_then(|t|fs_scenario::gesture::path::legacy_member(&t.target,"bow"))==Some(bow.string),
        };
        if !accepts
            || config.station_fraction.to_bits() != bow.station_fraction.to_bits()
            || config.sample_rate_hz != self.source.options.sample_rate_hz {
            return Err(GestureCompileError::Invalid {
                what: "string-bank track must match its actual string, fixed bow station and clock",
            });
        }
        let external = vec![0.0; self.system.mode_count()];
        let renderer = ScheduledRenderer::from_friction_gesture_admitted(self.system, external, schedule, track, config,path)?;
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
    assemble_sources(cards, None, bodies, connections, body_connections, &[], &[], None, None, &[], false, law, options, gate)
}

// Keep the original one-string API in the same assembler. Its absent source
// extensions retain the original basis and arithmetic exactly. Constant-loss
// bridges name their connections in the assembled order: every string's own
// connections come first, in string order.
#[allow(clippy::too_many_arguments)]
fn assemble_sources(
    cards: &[BowedStringCard], sources: Option<&[MaterialStringMember]>,
    bodies: Vec<ModalAcousticTimeModel>, connections: &[&[InteriorStringConnection]],
    body_connections: &[ModalConnection],
    translational_body_connections: &[polarization::StringBodyTranslationConnection],
    rotational_body_connections: &[torsion::StringBodyRotationConnection],
    fluid_loading: Option<&fluid::StringBankFluidLoading>, prepared: impl Into<ModalPreparation>,
    constant_loss: &[ConstantLossConnection], residual: bool, law: FrictionLaw,
    options: &StringBankOptions, gate: &CancelGate,
) -> Result<FrictionModalSystem, ReciprocalStringError> {
    poll(gate)?;
    if cards.is_empty() || cards.len() != connections.len() || bodies.is_empty()
        || sources.is_some_and(|s| s.len() != cards.len())
        || options.bow.string >= cards.len() || !interior(options.bow.station_fraction) {
        return Err(invalid("reciprocal assembly requires aligned strings, bodies and an interior bow"));
    }
    if residual && (sources.is_none()
        || !rotational_body_connections.is_empty() || !constant_loss.is_empty()
        || connections.iter().any(|c| c.len() != 1)) {
        // Only a single spring/dashpot bridge per string condenses here.
        return Err(invalid("residual flexibility condenses one plain bridge connection per material string"));
    }
    if let (true, Some(members)) = (residual, sources) {
        // Directed springs condense one by one (polarization.rs), which is exact only when the
        // omitted modes couple none of a string's connections: d_i.d_j G(s_i, s_j) = 0, so every
        // axis must be orthogonal to its string's bridge axis and to each other directed axis.
        for (i, a) in translational_body_connections.iter().enumerate() {
            let member = members.get(a.string)
                .ok_or_else(|| invalid("directed bridge names an unknown material string"))?;
            let others = translational_body_connections.iter().enumerate()
                .filter(|(j, b)| *j != i && b.string == a.string).map(|(_, b)| b.direction.components());
            let [ax, ay] = a.direction.components();
            if core::iter::once(self::residual::directions(&member.source).1).chain(others)
                .any(|[bx, by]| (ax*bx + ay*by).abs() > 1e-12) {
                return Err(invalid("with residual flexibility a directed bridge must be orthogonal to its string's other connections"));
            }
        }
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
            // This connection is transverse. Separately declared rotation ports
            // below supply torque; a centreline connection never invents one.
            let bridge_shape = match source {
                Some(source) => source.bridge_shape(card, connection.station_fraction)?,
                None => shape(card, connection.station_fraction),
            };
            // The omitted modes are in series with the bridge: the retained
            // modes see it softened (residual.rs). The source keeps k and c.
            let (stiffness_n_m, damping_n_s_m) = match source {
                Some(source) if residual => {
                    let condensed = self::residual::condense_bridge(source, card, connection)?;
                    (condensed.stiffness_n_m, condensed.damping_n_s_m)
                },
                _ => (connection.stiffness_n_m, connection.damping_n_s_m),
            };
            links.push(ModalConnection {
                left: ModalAttachment { component: index, shapes: bridge_shape },
                right: ModalAttachment { component: cards.len() + connection.body.component,
                    shapes: connection.body.shapes.clone() },
                stiffness_n_m, damping_n_s_m,
                rest_extension_m: connection.rest_extension_m,
            });
        }
        let mut modes: Vec<_> = (0..card.mode_count).map(|k| ModalAcousticMode {
            angular_frequency_rad_s: card.mode_omega_rad_s(k+1), damping_ratio: card.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        }).collect();
        if let Some(source) = source {
            source.append_extra_modes(card, &mut modes, gate)?;
            if let Some(loading) = fluid_loading {
                loading.add_rotary_damping(index, source, card, &mut modes)?;
            }
        }
        if modes.len() != total_modes {
            return Err(invalid("material string source differs from its complete modal basis"));
        }
        components.push(ModalAcousticTimeModel::try_new(options.sample_rate_hz, modes, options.string_budget)
            .map_err(|source| ReciprocalStringError::Network(ModalCouplingError::Component { component: index, source }))?);
    }
    for connection in translational_body_connections {
        poll(gate)?;
        let members = sources.ok_or_else(|| invalid("a directed bridge requires original material string sources"))?;
        links.push(connection.lower(cards, members, &bodies, residual)?);
    }
    for connection in rotational_body_connections {
        poll(gate)?;
        let members = sources.ok_or_else(|| invalid("a rotational bridge requires original material string sources"))?;
        links.push(connection.lower(cards, members, &bodies)?);
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
    if let Some(loading) = fluid_loading {
        let sources = sources.ok_or_else(|| invalid("cylinder loading requires original material string sources"))?;
        loading.append_components(cards, sources, &mut components, &mut links, options, gate)?;
    }
    let network = CoupledModalSystem::new_with_preparation(components, links, constant_loss,
        options.network, prepared.into(), gate).map_err(ReciprocalStringError::Network)?;
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
/// Curved shear-deformable material body in the shared reciprocal mechanics.
pub mod shell;

/// Released bridge response of the complete original string/body/air performer.
pub mod harmonic;
/// The shared material body's exterior impedance in the original joint solve.
pub mod radiation;
/// Distributed volume-air modes on the exact retained material-shell basis.
pub mod cavity;

/// Original material-interface laws at actual string/body bridge footprints.
pub mod interfaces;
pub use interfaces::MaterialBridgeInterface;

/// Several source-bound bow contacts without resetting the shared mechanics.
pub mod multiple;

/// Physical transverse directions in the same material/body/friction solve.
pub mod polarization;

/// Material-derived surface twisting in the same reciprocal contact solve.
pub mod torsion;

/// Material-gas inertia and causal viscous memory in the same reciprocal solve.
pub mod fluid;

/// The same material sources in the existing nonlinear impact/contact owner.
pub mod nonlinear;

/// Source-bound bending memory in the same string/body/contact operator.
pub mod relaxation;

/// The omitted bending modes' static compliance at bridges and stops.
pub mod residual;
