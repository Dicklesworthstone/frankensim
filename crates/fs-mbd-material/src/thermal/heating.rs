//! Shared enthalpy transport driving source-bound geometry and mechanical inertia.
//!
//! Every cell stores the SAME reference mass as its mechanical region. The
//! original material Cp chart supplies thermal storage, and its expansion law
//! supplies geometry at each actual transport endpoint. Internal heat, boundary
//! events and optional generated heat use the existing fs-thermal timeline.
//!
//! This is one-way, uniform-region, solid free-expansion coupling. Thermal
//! areas, links and volumes remain reference-configuration approximations;
//! no constrained stress, moving-body energy transfer, melt flow or thermal
//! geometry feedback is inferred. A configured inertia is not a live-state
//! momentum update. Uniform-cell/Biot and time-refinement obligations remain.

use super::{
    ExpandedBody,
    temperature::{RegionTemperature, expand_body_at_temperatures},
};
use crate::{MaterialInertiaError, ResolvedBody};
use fs_ga::Vec3;
use fs_matdb::MaterialCard;
use fs_material::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::state_point::MaterialPropertySelection;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use std::collections::{BTreeMap, BTreeSet};

/// Solved heat drives compatible constrained FEM shape, stress and inertia.
#[cfg(feature = "thermoelastic-coupling")]
pub mod constrained;
/// Continuing free-body motion on the same source-bound thermal timeline.
pub mod motion;

mod mechanism;
mod source;
pub use mechanism::{HeatedArticulated, LinkHeating};

/// Thermal description for one complete material region; no hand-entered mass.
#[derive(Clone, Debug)]
pub struct ThermalRegion<'a> {
    /// Exact reference region key. All regions are required exactly once.
    pub key: &'a str,
    /// Original immutable mechanical material card, also owning Cp and alpha.
    pub card: &'a MaterialCard,
    /// Source-resolved solid Cp/density chart with the same state context.
    pub caloric: &'a ResolvedHeatCapacityCurve,
    /// Stable control-volume identity in the thermal timeline.
    pub cell_id: CellId,
    /// Stable thermal-law identity. This compiler requires one per region.
    pub material_id: MaterialId,
    /// Original source segment identity, retained by thermal observations.
    pub segment_id: u64,
    /// Original part identity, retained by thermal observations.
    pub body_id: u64,
    /// Reference thermal geometry. Volume must match the mechanical region;
    /// line location and surface/contact coverage remain caller declarations.
    pub geometry: CellGeometry,
    /// Initial specific enthalpy [J/kg], on the caloric chart's explicit datum.
    /// Its temperature must agree with the mechanical reference temperature.
    pub initial_specific_enthalpy_j_kg: f64,
    /// Initial environment and applied watts; later changes are timed events.
    pub boundary: Boundary,
    /// Expansion claim selection at the unchanged reference temperature.
    pub expansion_selection: MaterialPropertySelection,
    /// Explicit free-expansion fixed point in the reference body frame [m].
    pub fixed_point_m: Vec3,
}

/// The original failing owner is retained; no partial configuration is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum HeatingError {
    /// Thermal input, constitutive domain, budget, energy or cancellation error.
    Thermal(DepositionError),
    /// Material-source, expansion, geometry or mechanical-inertia failure.
    Mechanical(MaterialInertiaError),
    /// Original caloric-state failure.
    Caloric(PhaseStateError),
    /// Incompatible original material, reference inventory or source context.
    Binding(&'static str),
}
impl From<DepositionError> for HeatingError {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl From<MaterialInertiaError> for HeatingError {
    fn from(e: MaterialInertiaError) -> Self {
        Self::Mechanical(e)
    }
}
impl From<PhaseStateError> for HeatingError {
    fn from(e: PhaseStateError) -> Self {
        Self::Caloric(e)
    }
}
impl core::fmt::Display for HeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(e) => write!(f, "heated body thermal state: {e}"),
            Self::Mechanical(e) => write!(f, "heated body configuration: {e}"),
            Self::Caloric(e) => write!(f, "heated body caloric state: {e}"),
            Self::Binding(e) => write!(f, "heated body binding: {e}"),
        }
    }
}
impl std::error::Error for HeatingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(e) => Some(e),
            Self::Mechanical(e) => Some(e),
            Self::Caloric(e) => Some(e),
            Self::Binding(_) => None,
        }
    }
}

/// One accepted thermal history and its current material-defined body geometry.
/// All laws and source cards are immutable; Clone is a restartable checkpoint.
/// Neither the private thermal state nor its material map can be replaced alone.
#[derive(Clone)]
pub struct HeatedBody<'a> {
    thermal: DepositionSimulation<'a>,
    mechanics: ExpansionObserver<'a>,
}
impl<'a> HeatedBody<'a> {
    /// Compile source-bound storage, reference mass and actual thermal transport.
    ///
    /// Births at t=0 are derived from the reference inventory, before the caller's
    /// ordered contact/boundary events. Later births refuse: adding matter is a
    /// different model. Caloric data may use a wider temperature window than the
    /// expansion law, but every accepted thermal endpoint must support BOTH.
    /// All original thermal budgets apply, including generated initialization
    /// events. Source tables retain their own existing size/error limits.
    pub fn new(
        reference: &'a ResolvedBody,
        mut regions: Vec<ThermalRegion<'a>>,
        events: Vec<TimedEvent>,
        budget: DepositionBudget,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, HeatingError> {
        checkpoint()?;
        if regions.len() != reference.regions().len() {
            return Err(HeatingError::Binding(
                "every mechanical region needs exactly one thermal cell",
            ));
        }
        regions.sort_by_key(|r| r.key);
        let thermal = initialize(
            &[(reference, regions.as_slice())],
            events,
            budget,
            config,
            &mut checkpoint,
        )?;
        let configuration = configure(reference, &regions, &thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(Self {
            thermal,
            mechanics: ExpansionObserver {
                reference,
                regions,
                configuration,
            },
        })
    }

    /// Original source-resolved mechanical mass inventory, unchanged by heating.
    pub const fn reference(&self) -> &'a ResolvedBody {
        self.mechanics.reference
    }
    /// Actual current thermal cells, links, enthalpies and time, read-only.
    pub const fn thermal(&self) -> &DepositionSimulation<'a> {
        &self.thermal
    }
    /// Current geometry, full inertia and integrated alpha source receipts.
    pub const fn configuration(&self) -> &ExpandedBody<'a> {
        &self.mechanics.configuration
    }
    /// Immutable original caloric and expansion source descriptions.
    pub fn regions(&self) -> &[ThermalRegion<'a>] {
        &self.mechanics.regions
    }

    /// Advance authored heating/cooling through the shared thermal solver.
    /// A rejected intermediate material state, even followed by cooling back
    /// into range, refuses the WHOLE advance. No resampling hides that state.
    pub fn advance_to(
        &self,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), HeatingError> {
        let (thermal, mechanics, audit) =
            self.thermal
                .advance_observed(&self.mechanics, target_s, config, checkpoint)?;
        Ok((Self { thermal, mechanics }, audit))
    }

    /// Compose existing electrical, frictional or other heat/transport observers
    /// on the SAME substeps, not a second thermal solve. Original typed errors
    /// and accepted observer state are preserved. Generated heat is counted once
    /// by the original thermal audit. Observer Clone must be rollback-safe.
    /// This is explicit thermal staggering, not strongly coupled mechanics.
    pub fn advance_with<O: DepositionObserver>(
        &self,
        observer: &O,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, EnergyAudit), ObserverCouplingError<O::Error, HeatingError>> {
        let pair = CoupledObservers::new(observer.clone(), self.mechanics.clone());
        let (thermal, pair, audit) = self
            .thermal
            .advance_observed(&pair, target_s, config, checkpoint)?;
        let (observer, mechanics) = pair.into_parts();
        Ok((Self { thermal, mechanics }, observer, audit))
    }
}

#[derive(Clone)]
struct ExpansionObserver<'a> {
    reference: &'a ResolvedBody,
    regions: Vec<ThermalRegion<'a>>,
    configuration: ExpandedBody<'a>,
}
impl DepositionObserver for ExpansionObserver<'_> {
    type Error = HeatingError;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if matches!(event, Event::Birth(_)) {
            return Err(HeatingError::Binding(
                "new mass is outside this fixed-body model",
            ));
        }
        // The admitted remaining events change contacts/boundaries, not h or T.
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
        self.configuration = configure(self.reference, &self.regions, thermal, checkpoint)?;
        Ok(())
    }
}

fn configure<'a>(
    reference: &'a ResolvedBody,
    regions: &[ThermalRegion<'_>],
    thermal: &DepositionSimulation<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<ExpandedBody<'a>, HeatingError> {
    checkpoint()?;
    let mut temperatures = Vec::with_capacity(regions.len());
    for region in regions {
        checkpoint()?;
        // The private timeline was constructed with this exact immutable chart.
        let cell = thermal
            .cells()
            .get(&region.cell_id)
            .ok_or(HeatingError::Binding("mapped thermal cell is absent"))?;
        let state = region
            .caloric
            .curve()
            .state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if !state.is_fully_solid() {
            return Err(HeatingError::Binding(
                "non-solid enthalpy cannot configure this rigid body",
            ));
        }
        temperatures.push(RegionTemperature {
            key: region.key,
            card: region.card,
            temperature_k: state.temperature_k(),
            selection: region.expansion_selection.clone(),
            fixed_point_m: region.fixed_point_m,
        });
    }
    // Preserve arbitrary original checkpoint errors across the bool adapter.
    let mut error = None;
    let configuration =
        expand_body_at_temperatures(reference, &temperatures, || match checkpoint() {
            Ok(()) => false,
            Err(e) => {
                error = Some(e);
                true
            }
        });
    if let Some(e) = error {
        return Err(e.into());
    }
    Ok(configuration?)
}

// One physical timeline, whether the inventory belongs to one body or many
// articulated links. Region keys are local; cell/law IDs and budgets are global.
fn initialize<'a>(
    inventories: &[(&ResolvedBody, &[ThermalRegion<'a>])],
    events: Vec<TimedEvent>,
    budget: DepositionBudget,
    config: StepConfig,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<DepositionSimulation<'a>, HeatingError> {
    let mut count = 0usize;
    for &(reference, regions) in inventories {
        checkpoint()?;
        if regions.len() != reference.regions().len() {
            return Err(HeatingError::Binding(
                "every mechanical region needs exactly one thermal cell",
            ));
        }
        count = count
            .checked_add(regions.len())
            .ok_or(DepositionError::Budget("heated body schedule"))?;
    }
    let total = count
        .checked_add(events.len())
        .ok_or(DepositionError::Budget("heated body schedule"))?;
    if total > budget.max_events || count > budget.max_cells {
        return Err(DepositionError::Budget("heated body initialization").into());
    }
    let mut cells = BTreeSet::new();
    let mut materials = BTreeSet::new();
    let mut schedule = Vec::with_capacity(total);
    let mut laws: Vec<(MaterialId, &'a dyn ThermalLaw)> = Vec::with_capacity(count);
    for &(reference, regions) in inventories {
        for (original, region) in reference.regions().iter().zip(regions) {
            checkpoint()?;
            if original.key() != region.key
                || !cells.insert(region.cell_id)
                || !materials.insert(region.material_id)
            {
                return Err(HeatingError::Binding(
                    "region keys, cell IDs and material IDs must be unique and complete",
                ));
            }
            source::admit(original, region, checkpoint)?;
            laws.push((region.material_id, region.caloric.curve()));
            schedule.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: region.cell_id,
                    segment_id: region.segment_id,
                    body_id: region.body_id,
                    material_id: region.material_id,
                    mass_kg: original.inertia().mass(),
                    specific_h_j_kg: region.initial_specific_enthalpy_j_kg,
                    geometry: region.geometry,
                    boundary: region.boundary,
                }),
            });
        }
    }
    for event in events {
        checkpoint()?;
        if matches!(&event.event, Event::Birth(_)) {
            return Err(HeatingError::Binding(
                "a fixed reference body cannot accept new mass",
            ));
        }
        schedule.push(event);
    }
    let initial = DepositionSimulation::new(schedule, &laws, budget)?;
    let (thermal, _) = initial.advance_to(0.0, config, &mut *checkpoint)?;
    Ok(thermal)
}
