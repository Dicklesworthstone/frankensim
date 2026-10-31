//! Heat-connected articulated links on one source-bound enthalpy timeline.
use super::*;
use crate::articulated::ResolvedArticulated;
use crate::thermal::{ExpandedArticulated, JointMaterialAttachment};

/// Continuing joint motion with source-resolved inertia on the thermal clock.
pub mod motion;

#[cfg(feature = "electrothermal-coupling")]
mod electrical;

/// Complete thermal description for one named mechanical link.
/// Region keys are local to the link; cell and material IDs are globally unique.
#[derive(Clone, Debug)]
pub struct LinkHeating<'a> {
    /// Exact reference link name, including fixed links and the root.
    pub link: &'a str,
    /// One thermal cell and original material description per reference region.
    pub regions: Vec<ThermalRegion<'a>>,
}

/// One thermal history and one coherent configuration of a complete mechanism.
///
/// Inter-link heat transport uses ordinary `Event::Contact` rows between the
/// globally named cells. Mechanical joints do NOT imply thermal contact, nor
/// does a thermal link create a mechanical constraint. Both graphs are explicit.
/// Every link uses its own material, temperature and fixed mass inventory.
/// This inherits the single-body adapter's solid free-expansion approximation.
/// Joint placements may follow explicit material attachments; other placements
/// and thermal geometry stay prescribed. Forward dynamics consumes the model;
/// q/qdot and their
/// time integration, expansion work and mechanical back-reaction are not owned
/// by this adapter. The energy audit covers the thermal timeline only.
#[derive(Clone)]
pub struct HeatedArticulated<'a> {
    thermal: DepositionSimulation<'a>,
    mechanics: MechanismExpansion<'a>,
}

impl<'a> HeatedArticulated<'a> {
    /// Compile all links into ONE thermal timeline using their reference masses.
    ///
    /// Supply exactly one binding for every link; input order is immaterial.
    /// Cell, law, event, contact and substep limits apply to the whole mechanism,
    /// not independently per body. Births are generated before caller events at
    /// t=0, so contacts may cross links immediately. Later material births refuse.
    /// Each caloric source must match its own link's original card and state;
    /// different links need not share reference temperature or composition.
    pub fn new(
        reference: &'a ResolvedArticulated,
        bindings: Vec<LinkHeating<'a>>,
        events: Vec<TimedEvent>,
        budget: DepositionBudget,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, HeatingError> {
        checkpoint()?;
        if bindings.len() != reference.model().link_count() {
            return Err(HeatingError::Binding("every mechanical link needs exactly one thermal assignment"));
        }
        // Every admitted reference link has at least one mass region. Check the
        // global lower bound before creating the binding lookup or new cells.
        if bindings.len() > budget.max_cells {
            return Err(DepositionError::Budget("heated mechanism initialization").into());
        }
        let mut by_name = BTreeMap::new();
        for binding in bindings {
            checkpoint()?;
            if binding.link.trim().is_empty() || binding.link.len() > 128 {
                return Err(HeatingError::Binding("thermal link names must be nonblank and bounded"));
            }
            if by_name.insert(binding.link, binding.regions).is_some() {
                return Err(HeatingError::Binding("duplicate thermal link assignment"));
            }
        }
        let mut regions = Vec::with_capacity(reference.model().link_count());
        for (link, body) in reference.model().links().iter().zip(reference.bodies()) {
            checkpoint()?;
            let mut assigned = by_name.remove(link.name())
                .ok_or(HeatingError::Binding("missing or foreign thermal link assignment"))?;
            if assigned.len() != body.regions().len() {
                return Err(HeatingError::Binding("every mechanical region needs exactly one thermal cell"));
            }
            assigned.sort_by_key(|r| r.key);
            regions.push(assigned);
        }
        let inventories: Vec<_> = reference.bodies().iter()
            .zip(regions.iter().map(Vec::as_slice)).collect();
        let thermal = initialize(&inventories, events, budget, config, &mut checkpoint)?;
        let configuration =
            configure_mechanism(reference, &regions, &thermal, &[], &mut checkpoint)?;
        checkpoint()?;
        Ok(Self { thermal, mechanics: MechanismExpansion { reference, regions, configuration } })
    }

    /// Original mechanism with all immutable material and reference mass records.
    pub const fn reference(&self) -> &'a ResolvedArticulated { self.mechanics.reference }
    /// Shared time, control volumes, links, enthalpy and thermal observations.
    pub const fn thermal(&self) -> &DepositionSimulation<'a> { &self.thermal }
    /// Current full mechanism, directly usable by the existing dynamics solver.
    pub const fn configuration(&self) -> &ExpandedArticulated<'a> { &self.mechanics.configuration }
    /// Original thermal descriptions in link order, then canonical region order.
    pub fn regions(&self) -> &[Vec<ThermalRegion<'a>>] { &self.mechanics.regions }

    /// Attach selected joint origins to their actual source-expanded material.
    /// Each subsequent thermal substep rebuilds those placements from the same
    /// original reference geometry as inertia. The thermal history is unchanged.
    /// Configure attachments before constructing a continuing motion owner; this
    /// does not remap a live mechanical state or supply its parameter work.
    pub fn with_joint_attachments(
        &self,
        attachments: &[JointMaterialAttachment],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, HeatingError> {
        checkpoint()?;
        let mut error = None;
        let configuration = self
            .configuration()
            .with_joint_attachments(attachments, || match checkpoint() {
                Ok(()) => false,
                Err(e) => {
                    error = Some(e);
                    true
                }
            });
        if let Some(e) = error {
            return Err(e.into());
        }
        let configuration = configuration?;
        let mut next = self.clone();
        next.mechanics.configuration = configuration;
        checkpoint()?;
        Ok(next)
    }

    /// Evolve shared contact/boundary heat and reconfigure EVERY actual substep.
    /// A late material refusal rolls back all links and the thermal timeline,
    /// not merely the failing link. Clone retains an immutable restart point.
    pub fn advance_to(
        &self, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), HeatingError> {
        let (thermal, mechanics, audit) = self.thermal.advance_observed(
            &self.mechanics, target_s, config, checkpoint,
        )?;
        Ok((Self { thermal, mechanics }, audit))
    }

    /// Add a rollback-safe heat/transport observer on the SAME thermal substeps.
    /// Sources on any link are integrated once by the original thermal owner.
    /// Original source errors, observer state and total thermal work are retained.
    /// This does not solve an implicit mechanical/thermal feedback problem.
    pub fn advance_with<O: DepositionObserver>(
        &self, observer: &O, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, EnergyAudit), ObserverCouplingError<O::Error, HeatingError>> {
        let pair = CoupledObservers::new(observer.clone(), self.mechanics.clone());
        let (thermal, pair, audit) = self.thermal.advance_observed(
            &pair, target_s, config, checkpoint,
        )?;
        let (observer, mechanics) = pair.into_parts();
        Ok((Self { thermal, mechanics }, observer, audit))
    }
}

#[derive(Clone)]
struct MechanismExpansion<'a> {
    reference: &'a ResolvedArticulated,
    regions: Vec<Vec<ThermalRegion<'a>>>,
    configuration: ExpandedArticulated<'a>,
}
impl DepositionObserver for MechanismExpansion<'_> {
    type Error = HeatingError;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if matches!(event, Event::Birth(_)) {
            return Err(HeatingError::Binding("new mass is outside this fixed-mechanism model"));
        }
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.configuration = configure_mechanism(
            self.reference,
            &self.regions,
            thermal,
            self.configuration.joint_attachments(),
            checkpoint,
        )?;
        Ok(())
    }
}

fn configure_mechanism<'a>(
    reference: &'a ResolvedArticulated,
    regions: &[Vec<ThermalRegion<'_>>],
    thermal: &DepositionSimulation<'_>,
    attachments: &[JointMaterialAttachment],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<ExpandedArticulated<'a>, HeatingError> {
    checkpoint()?;
    let mut bodies = Vec::with_capacity(reference.bodies().len());
    for (original, regions) in reference.bodies().iter().zip(regions) {
        // Reuse the exact single-body source/state/geometry path. No law is
        // resampled through a second model, and reference mass is never updated.
        bodies.push(configure(original, regions, thermal, checkpoint)?);
    }
    let mut error = None;
    let configuration = ExpandedArticulated::from_bodies_with_joint_attachments(
        reference,
        bodies,
        attachments,
        || {
        match checkpoint() { Ok(()) => false, Err(e) => { error = Some(e); true } }
    },
    );
    if let Some(e) = error { return Err(e.into()); }
    Ok(configuration?)
}
