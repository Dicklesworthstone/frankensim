//! Several material-derived bow contacts on the SAME retained string/body state.
//!
//! Contact shapes come from the selected material strings, not supplied modal
//! vectors. The original simultaneous-friction owner computes all reactions,
//! work, sliding loss and pressure together. This module adds no force law or
//! time integration. Adding contacts preserves existing vibration and clocks.
//!
//! Loads are prescribed by default. Explicit co-located normal-contact bindings
//! instead resolve each load from contact mechanics. Optional receiving-body
//! attachments also return tangential traction and resolve surface motion;
//! a shared bow-hair structure is not inferred by this adapter.
//! Explicit station paths can translate prescribed-load bow contacts along the
//! original strings; independently admitted compliant stops can also migrate.
//! The underlying bounded nonlinear solve may
//! refuse a step; global uniqueness and measured instrument timbre are not claimed.

use super::*;
use crate::render::schedule::force::coupled::contact::friction::ModalFrictionPort;
use crate::render::schedule::force::coupled::contact::friction::multiple::{
    MultiFrictionConfig, MultiFrictionFrame, MultiFrictionModalSystem,
};

/// One physical contact on a selected material string and its original law.
/// The law retains caller authority; string material data do not certify rosin.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBowContact {
    /// Original string index and fixed interior station along its full span.
    pub location: StringBowLocation,
    /// Original dry-friction construction law. An explicitly bound material
    /// shear layer replaces its active traction mechanism on that contact.
    pub law: FrictionLaw,
    /// Original per-contact physical and nonlinear-solve limits.
    pub friction: ModalFrictionConfig,
}

/// Retained material, body, cavity and gas bindings plus all contact inputs.
#[derive(Clone, Debug)]
pub struct MultiBowStringBankSource {
    bank: MaterialStringBankSource,
    bows: Vec<MaterialBowContact>,
    joint: MultiFrictionConfig,
    stops: Vec<stops::MaterialStringStop>,
    bow_contact_loads: Vec<contact_load::BowContactLoad>,
    ribbon: Option<ribbon::UniformBowRibbon>,
    compliant_ribbon: Option<compliant_ribbon::CompliantBowRibbon>,
    fiber_ribbon: Option<ribbon::MaterialFiberRibbon>,
    body_surfaces: Vec<surfaces::BowBodySurface>,
}
impl MultiBowStringBankSource {
    /// Original strings, layers, regional/uniform body, cavity, gas and supports.
    #[must_use]
    pub const fn bank(&self) -> &MaterialStringBankSource { &self.bank }
    /// Actual drive order: the original bow FIRST, then the added contacts.
    #[must_use]
    pub fn bows(&self) -> &[MaterialBowContact] { &self.bows }
    /// Original simultaneous-friction owner's declared setup and sweep limits.
    #[must_use]
    pub const fn joint(&self) -> &MultiFrictionConfig { &self.joint }
}

/// One mechanical state and pressure observer with jointly solved bow contacts.
/// A released bow removes only traction; no string or body vibration is muted.
pub struct MultiBowStringBank {
    source: MultiBowStringBankSource,
    cards: Vec<BowedStringCard>,
    system: MultiFrictionModalSystem,
    stop_stations: Option<stops::stations::StationState>,
    bow_stations: Option<stations::BowStationState>,
}
impl ReciprocalStringBank {
    /// Add material-derived contacts without reducing or resetting any component.
    ///
    /// The original bow remains port zero. Additional locations use the original
    /// string declaration order and its actual material-derived displacement
    /// basis. They may contact different strings or distinct stations on the
    /// same string. All string, bridge, body and cavity reactions participate in
    /// the original owner's simultaneous solve, not sequential voice steps.
    ///
    /// Existing modal state, accepted clock and physical pressure weights survive
    /// the conversion. Empty additions use the original single-port step exactly.
    /// No source claims, state histories or observer normalizations are replaced.
    ///
    /// # Errors
    /// Unknown strings, endpoint/nonfinite stations, port/setup limits, invalid
    /// friction laws or per-contact budgets, numerical response or cancellation.
    pub fn with_additional_bows(
        self, additional: Vec<MaterialBowContact>, joint: MultiFrictionConfig,
        gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        let count = additional.len().checked_add(1)
            .ok_or_else(|| invalid("material bow count overflow"))?;
        // Enforce the existing owner's port cap before allocating modal shapes.
        if !(1..=64).contains(&joint.max_ports) || count > joint.max_ports {
            return Err(invalid("material bow contacts exceed the joint port budget"));
        }
        for bow in &additional {
            if bow.location.string >= self.cards.len() || !interior(bow.location.station_fraction) {
                return Err(invalid("each added bow must name an actual material string and interior station"));
            }
        }
        let mut ports = Vec::with_capacity(additional.len());
        for bow in &additional {
            poll(gate)?;
            ports.push((ModalFrictionPort {
                left: ModalAttachment { component: bow.location.string,
                    shapes: self.source.strings[bow.location.string].source
                        .surface_shape(&self.cards[bow.location.string], bow.location.station_fraction)? },
                right: None, law: bow.law.clone(),
            }, bow.friction));
        }
        let mut bows = Vec::with_capacity(count);
        bows.push(MaterialBowContact { location: self.source.options.bow,
            law: self.system.port().law.clone(), friction: self.source.options.friction });
        bows.extend(additional);
        let system = self.system.with_additional_ports(ports, joint, gate)
            .map_err(ReciprocalStringError::Friction)?;
        poll(gate)?;
        Ok(MultiBowStringBank { source: MultiBowStringBankSource {
            bank: self.source, bows, joint, stops: Vec::new(), bow_contact_loads: Vec::new(),
            ribbon: None, compliant_ribbon: None, fiber_ribbon: None, body_surfaces: Vec::new(),
        }, cards: self.cards, system, stop_stations: None, bow_stations: None })
    }
}
impl MultiBowStringBank {
    /// Bind complete material strain histories after optional static preload.
    /// The relaxed initial history uses that actual solved displacement; all
    /// bow, stop, surface and original mechanical state remain in their owner.
    pub fn with_material_finite_strain_longitudinal_relaxation(mut self,
        selections:Vec<super::nonlinear::longitudinal::relaxation::LongitudinalRelaxationSelection>,
        initial_memory:Option<Vec<f64>>,max_memory_states:usize,gate:&CancelGate)
        -> Result<Self,crate::render::plate::impact::ImpactError> {
        use super::nonlinear::longitudinal::relaxation::finite_strain;
        use crate::render::plate::impact::ImpactError;
        let (image,bindings)=finite_strain::bindings(&self.source.bank,&self.cards,
            self.system.original_core(),selections,initial_memory,max_memory_states,gate)?;
        if !bindings.is_empty() {
            self.system=self.system.with_mapped_material_histories(bindings,max_memory_states,gate)
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
        }
        finite_strain::retain_source(&mut self.source.bank,image);
        finite_strain::validate_sources(&self.source.bank,self.system.original_core(),gate)?;
        // Fibers were assembled before the force/preload program. Their
        // original source spectrum must also admit the expanded string band.
        if let Some(memory)=self.source.fiber_ribbon.as_ref()
            .and_then(|ribbon|ribbon.model().material_relaxation()) {
            memory.validate_band(self.system.linear_frequency_upper_bound_rad_s(),
                self.system.sample_period_s()).map_err(|e|ImpactError::Owner(e.to_string()))?;
        }
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok(self)
    }

    /// Original material/geometry bindings and every contact in drive order.
    #[must_use]
    pub const fn source(&self) -> &MultiBowStringBankSource { &self.source }
    /// Original transverse coefficients; torsional sources remain attached separately.
    #[must_use]
    pub fn cards(&self) -> &[BowedStringCard] { &self.cards }
    /// Read-only state, physical pressure weights and per-contact work reports.
    #[must_use]
    pub const fn system(&self) -> &MultiFrictionModalSystem { &self.system }
    /// Override the active-set rescue pivot budget (frankensim-o7qbp) so a test
    /// can plant a budget-exhaustion negative; production uses the default.
    pub fn force_active_set_budget(&mut self, budget: usize) { self.system.force_active_set_budget(budget); }
    /// One physical sample for ALL strings, bodies and contacts. The drive slice
    /// must match `source().bows()` exactly. A prescribed zero load releases its
    /// bow unless a retained material shear layer is bound (opening then refuses).
    /// For contact-loaded bows, zero means no additional prescribed load
    /// and separation of the bound normal contact supplies the actual release.
    pub fn step(&mut self, external: &[f64], drives: &[FrictionDrive])
        -> Result<&MultiFrictionFrame, ModalFrictionError> {
        if self.bow_stations.is_some() {
            self.step_with_bow_geometry(external, drives, None, None, None, &CancelGate::new_clock_free())
        } else { self.system.step(external, drives) }
    }
    /// Same atomic physical sample with the original owner's cancellation polls.
    pub fn step_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive], gate: &CancelGate)
        -> Result<&MultiFrictionFrame, ModalFrictionError> {
        if self.bow_stations.is_some() {
            self.step_with_bow_geometry(external, drives, None, None, None, gate)
        } else { self.system.step_under_gate(external, drives, gate) }
    }
    /// Transfer the original sources together with the complete accepted state.
    #[must_use]
    pub fn into_parts(self) -> (MultiBowStringBankSource, MultiFrictionModalSystem) { (self.source, self.system) }
}

/// Common gesture-clock playback of all contacts and the shared pressure observer.
pub mod performance;

/// Material-derived unilateral stops in the same bow/string/body mechanics.
pub mod stops;

/// Source-derived bow contact migration on the unchanged material bases.
pub mod stations;


/// Source/interface applicability retained through the existing bow performers.
mod interface;


/// Co-located normal and tangential material-string contacts in one solve.
pub mod contact_load;

/// Finite-width prescribed ribbons resolved by the original joint contact solver.
pub mod ribbon;

/// Compliant bow-hair ribbon: per-cell spring-dashpot hair compliance derived
/// from hair data (Pitteroff & Woodhouse 1998); differential slip emerges.
pub mod compliant_ribbon;

/// Resolved receiving-body motion and opposite tangential contact reactions.
pub mod surfaces;

/// Source-bound isothermal material shear in the original pressure renderer.
mod shear_layer;

/// Explicit finite thermal strips following the source string/fiber contact.
pub use shear_layer::dry::{BowDryContactHeat, DryReceiverHeatPosition};
