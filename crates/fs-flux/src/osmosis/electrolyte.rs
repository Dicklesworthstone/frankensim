//! Finite ionic crystals on the original osmotic water timeline.
//!
//! Every membrane retains all analytical solutes. Candidate solvent volumes
//! re-speciate the original fast aqueous chemistry with the actual mineral
//! phases frozen. The original mineral kernel supplies the surrounding phase
//! half-stages. No charged permeation, reaction water, heat or concentrated
//! electrolyte activity law is inferred.

use super::*;
use fs_material::osmosis::solubility::electrolyte::{
    ElectrolyteSolubilityError, kinetics::KineticElectrolyteSolubilitySolution,
};
use fs_qty::{Amount, Dimensionless};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::thermodynamics::acid_base::{
        AcidBaseError, precipitation::kinetics::MineralKineticControlsV1,
    },
};

/// An analytical inventory, distinct from its changing aqueous species forms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElectrolyteInventory {
    /// Selected mineral cation, including its original solid coefficient.
    Cation(SpeciesId),
    /// Original acid/base family, including any mineral contribution.
    Family(String),
    /// Nonreacting original spectator; its charge is never inferred from a name.
    Spectator { species: SpeciesId, charge: i32 },
}

/// One original analytical ledger and the independently recomputed phase balance.
#[derive(Clone, Debug)]
pub struct ElectrolyteInventoryBalance {
    pub inventory: ElectrolyteInventory,
    /// Complete original amount [mol]; the source ledger is retained exactly.
    pub total_mol: f64,
    /// Sum of actual final-minus-initial dissolved and solid contributions [mol].
    /// This is not the difference of two rounded complete totals.
    pub residual_mol: f64,
}

/// Per-cell conservation for a fully retaining reactive ionic compartment.
#[derive(Clone, Debug)]
pub struct ElectrolyteCompartmentBalance {
    pub compartment: usize,
    pub inventories: Vec<ElectrolyteInventoryBalance>,
    /// Actual final-minus-initial formula stock [mol], never an osmotic particle.
    pub solid_formula_change_mol: f64,
    /// Original final aqueous charge residual [mol charge equivalents].
    pub charge_balance_residual_mol: f64,
}

#[derive(Clone, Debug)]
pub(super) struct KineticElectrolyte {
    pub(super) material: KineticElectrolyteSolubilitySolution,
    pub(super) controls: MineralKineticControlsV1,
}

impl From<ElectrolyteSolubilityError> for OsmoticNetworkError {
    fn from(error: ElectrolyteSolubilityError) -> Self {
        match error {
            ElectrolyteSolubilityError::Speciation(AcidBaseError::Cancelled) => Self::Cancelled,
            error => Self::ElectrolyteSolubility(error),
        }
    }
}

impl OsmoticCompartment {
    /// Retain an original finite-rate ionic mineral and its actual aqueous view.
    ///
    /// All cells of the resulting network must share the original aqueous
    /// model and activity basis; clone that source when preparing other cells.
    /// Bare membranes must fully reflect every dissolved species and have zero
    /// solute permeance. Water alone moves; acid/base speciation and finite-rate
    /// mineral exchange change the dissolved particle counts in each cell.
    ///
    /// Native kinetic controls bound phase work. Their amount/relative accuracy
    /// is tightened when necessary to sit below the outer coupled-step scales.
    /// Root accuracy uses the actual before/endpoint phase scales; at most
    /// three native trials of the same interval qualify that allowance. Every
    /// attempted call is counted, and only the qualified candidate is retained.
    /// Speciation and phase publication allocate original chemical payloads;
    /// this constructor does not select an allocation-free transport path.
    /// The network byte budget counts retained numerical state; immutable
    /// pre-existing shared sources and transient chemical/report allocations
    /// are outside that estimate.
    pub fn with_kinetic_electrolyte_solubility(
        material: KineticElectrolyteSolubilitySolution,
        kinetics: MineralKineticControlsV1,
        pressure: OsmoticPressureVolume,
        minimum_volume_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        // The original owner validates every phase/aqueous work control. A
        // zero-time validation neither moves crystals nor projects equilibrium.
        material.advance(Time::new(0.0), kinetics, || false)?;
        let mut cell = Self::new(material.solution().clone(), pressure, minimum_volume_m3)?;
        cell.kinetic_electrolyte = Some(KineticElectrolyte {
            material,
            controls: kinetics,
        });
        Ok(cell)
    }

    /// Independent checkpoint with the exact original source, analytical
    /// ledgers, actual phases and fast aqueous speciation retained.
    pub fn kinetic_electrolyte_solubility_snapshot(
        &self,
    ) -> Option<KineticElectrolyteSolubilitySolution> {
        self.kinetic_electrolyte
            .as_ref()
            .map(|phase| phase.material.clone())
    }
}

pub(super) fn admit(
    species: &[String],
    cells: &[OsmoticCompartment],
    membranes: &[OsmoticMembrane],
) -> Result<(), OsmoticNetworkError> {
    let Some(reference) = cells
        .iter()
        .find_map(|cell| cell.kinetic_electrolyte.as_ref())
    else {
        return Ok(());
    };
    let original = reference.material.aqueous().model();
    for cell in cells {
        let phase = cell
            .kinetic_electrolyte
            .as_ref()
            .ok_or(OsmoticNetworkError::Invalid(
                "reactive ionic water networks require original ionic owners in every compartment",
            ))?;
        let material = &phase.material;
        let model = material.aqueous().model();
        // Both borrowed fields live inside the original Arc<Layout>. Pointer
        // identity admits its clones, not a new water/pKa source with matching
        // labels or a numerically similar equilibrium constant.
        if !core::ptr::eq(model.water(), original.water())
            || model.families().as_ptr() != original.families().as_ptr()
            || model.families().len() != original.families().len()
            || material.basis() != reference.material.basis()
            || material.species().len() != species.len()
            || material
                .species()
                .iter()
                .zip(species)
                .any(|(id, name)| id.as_str() != name)
            || cell.solubility.is_some()
            || cell.kinetic_solubility.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "reactive ionic cells need the same original aqueous model, activity basis and actual species ordering",
            ));
        }
    }
    if membranes.iter().any(|membrane| {
        !matches!(membrane.material, MembraneBoundary::WellMixed(_))
            || membrane
                .membrane_law()
                .solutes()
                .iter()
                .any(|(reflection, permeance)| {
                    reflection.value() != 1.0 || permeance.value() != 0.0
                })
    }) {
        return Err(OsmoticNetworkError::Invalid(
            "reactive ionic minerals require bare fully reflecting membranes with zero dissolved-species permeance",
        ));
    }
    Ok(())
}

pub(super) fn install(
    cell: &mut OsmoticCompartment,
    material: KineticElectrolyteSolubilitySolution,
) -> Result<(), OsmoticNetworkError> {
    let owner = cell
        .kinetic_electrolyte
        .as_mut()
        .ok_or(OsmoticNetworkError::Invalid(
            "reactive ionic workspace layout",
        ))?;
    cell.solution.replace_inventories(
        material.solution().solvent_volume_m3(),
        material.solution().solute_moles(),
    )?;
    owner.material = material;
    Ok(())
}

pub(super) fn advance_phase(
    cell: &mut OsmoticCompartment,
    dt: f64,
    outer: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    let Some(phase) = &cell.kinetic_electrolyte else {
        return Ok(());
    };
    poll(cancelled)?;
    let mut controls = phase.controls;
    let absolute = 0.01 * outer.absolute_amount_tolerance_mol;
    let relative = 0.01 * outer.relative_tolerance;
    if !nonnegative(absolute) || !positive(relative) {
        return Err(OsmoticNetworkError::Invalid(
            "ionic phase accuracy loses representable range",
        ));
    }
    controls.amount_absolute_tolerance =
        Amount::new(controls.amount_absolute_tolerance.value().min(absolute));
    controls.relative_tolerance =
        Dimensionless::new(controls.relative_tolerance.value().min(relative));
    // Match the outer norm's max(before, endpoint) stock scale, including a
    // phase that first appears in this stage. An initial zero stock must not
    // collapse useful relative accuracy to an arbitrary 1e-20 mol root limit.
    // Start with the caller's extent ceiling; ordinary settings qualify after
    // one native call. Only an insufficient endpoint allowance repeats the
    // SAME duration from the unchanged initial state, at tighter root accuracy.
    // This is bounded accuracy refinement, not a second time integrator.
    for _ in 0..3 {
        poll(cancelled)?;
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(OsmoticNetworkError::WorkBudget)?;
        let (material, _) = phase
            .material
            .advance(Time::new(dt), controls, &mut *cancelled)?
            .into_parts();
        let extent = endpoint_extent_allowance(&phase.material, &material, absolute, relative);
        if !positive(extent) {
            return Err(OsmoticNetworkError::Invalid(
                "ionic formula-root accuracy loses representable range",
            ));
        }
        if controls.extent_absolute_tolerance.value() <= extent {
            return install(cell, material);
        }
        // Leave margin for the smaller endpoint shift on the refined call.
        // This can only tighten the original caller-selected ceiling.
        controls.extent_absolute_tolerance = Amount::new(0.5 * extent);
        if !positive(controls.extent_absolute_tolerance.value()) {
            return Err(OsmoticNetworkError::Invalid(
                "ionic formula-root refinement loses range",
            ));
        }
    }
    Err(OsmoticNetworkError::NonlinearConvergence)
}

fn endpoint_extent_allowance(
    before: &KineticElectrolyteSolubilitySolution,
    after: &KineticElectrolyteSolubilitySolution,
    absolute: f64,
    relative: f64,
) -> f64 {
    // Map physical cation/family and independent crystal scales to formula
    // moles through the original reaction coefficients. This is root-location
    // accuracy, not a pH error certificate: original aqueous uncertainty and
    // the outer actual-particle comparison remain separate acceptance gates.
    let source = before.model().source();
    let p = f64::from(source.cation_stoichiometry());
    let solid = before
        .solid_amount()
        .value()
        .max(after.solid_amount().value());
    let metal = before
        .dissolved_cation()
        .value()
        .max(after.dissolved_cation().value());
    let mut extent = (absolute + relative * solid).min(absolute / p + relative * (metal / p));
    if let (Some(before), Some(after)) = (before.dissolved_family(), after.dissolved_family()) {
        let q = f64::from(source.anion_stoichiometry());
        extent = extent.min(absolute / q + relative * (before.value().max(after.value()) / q));
    }
    extent
}

fn at_volume(
    cell: &OsmoticCompartment,
    volume: f64,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<KineticElectrolyteSolubilitySolution, OsmoticNetworkError> {
    let phase = cell
        .kinetic_electrolyte
        .as_ref()
        .ok_or(OsmoticNetworkError::Invalid(
            "reactive ionic water exchange requires two original ionic owners",
        ))?;
    poll(cancelled)?;
    if volume == cell.solution.solvent_volume_m3() {
        return Ok(phase.material.clone());
    }
    *evaluations = evaluations
        .checked_add(1)
        .ok_or(OsmoticNetworkError::WorkBudget)?;
    Ok(phase
        .material
        .with_solvent_volume_m3_cancelled(volume, &mut *cancelled)?)
}

pub(super) fn pair_endpoint(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    area_time: f64,
    water: f64,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    let volumes = [
        left.solution.solvent_volume_m3() - water,
        right.solution.solvent_volume_m3() + water,
    ];
    if !positive(volumes[0])
        || !positive(volumes[1])
        || volumes[0] < left.minimum_volume_m3
        || volumes[1] < right.minimum_volume_m3
    {
        return Err(OsmoticNetworkError::Unbracketed);
    }
    let phases = [
        at_volume(left, volumes[0], evaluations, cancelled)?,
        at_volume(right, volumes[1], evaluations, cancelled)?,
    ];
    let count = left.solution.solute_moles().len();
    let (mut cl, mut cr) = ([0.0; MAX_OSMOTIC_SPECIES], [0.0; MAX_OSMOTIC_SPECIES]);
    phases[0].solution().concentrations_into(&mut cl[..count])?;
    phases[1].solution().concentrations_into(&mut cr[..count])?;
    let flux = membrane
        .membrane_law()
        .water_flux(
            left.solution.temperature_k(),
            left.pressure
                .pressure_difference_pa(volumes[0], right.pressure, volumes[1])?,
            0.0,
            &cl[..count],
            &cr[..count],
        )?
        .volume_flux_m_s;
    let residual_m3 = water - area_time * flux;
    if !residual_m3.is_finite() {
        return Err(OsmoticNetworkError::Invalid(
            "ionic implicit water residual is not representable",
        ));
    }
    let (mut nl, mut nr) = ([0.0; MAX_OSMOTIC_SPECIES], [0.0; MAX_OSMOTIC_SPECIES]);
    nl[..count].copy_from_slice(phases[0].solution().solute_moles());
    nr[..count].copy_from_slice(phases[1].solution().solute_moles());
    Ok(PairEndpoint {
        volumes,
        left: nl,
        right: nr,
        residual_m3,
        electrolytes: Some(phases),
    })
}

pub(super) fn balances(
    before: &[OsmoticCompartment],
    after: &[OsmoticCompartment],
    controls: OsmoticStepControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<ElectrolyteCompartmentBalance>, OsmoticNetworkError> {
    let mut result = Vec::with_capacity(before.len());
    for (compartment, (old, new)) in before.iter().zip(after).enumerate() {
        poll(cancelled)?;
        let old = &old
            .kinetic_electrolyte
            .as_ref()
            .ok_or(OsmoticNetworkError::Invalid("reactive ionic ledger layout"))?
            .material;
        let new = &new
            .kinetic_electrolyte
            .as_ref()
            .ok_or(OsmoticNetworkError::Invalid("reactive ionic ledger layout"))?
            .material;
        if old.component_total() != new.component_total()
            || old.family_totals() != new.family_totals()
            || old.spectators() != new.spectators()
        {
            return Err(OsmoticNetworkError::Invalid(
                "original reactive ionic analytical ledger changed across a retaining membrane",
            ));
        }
        let solid = new.solid_amount().value() - old.solid_amount().value();
        let source = old.model().source();
        let mut inventories =
            Vec::with_capacity(1 + old.family_totals().len() + old.spectators().len());
        let mut record = |inventory,
                          total: f64,
                          aqueous_delta: f64,
                          solid_delta: f64|
         -> Result<(), OsmoticNetworkError> {
            let (mut residual, mut correction) = (0.0, 0.0);
            compensated_add(&mut residual, &mut correction, aqueous_delta);
            compensated_add(&mut residual, &mut correction, solid_delta);
            let scale =
                controls.absolute_amount_tolerance_mol + controls.relative_tolerance * total;
            if !residual.is_finite() || residual.abs() > scale {
                return Err(OsmoticNetworkError::Invalid(
                    "reactive ionic phase balance exceeds the declared amount scale",
                ));
            }
            inventories.push(ElectrolyteInventoryBalance {
                inventory,
                total_mol: total,
                residual_mol: residual,
            });
            Ok(())
        };
        record(
            ElectrolyteInventory::Cation(source.cation().clone()),
            old.component_total().value(),
            new.dissolved_cation().value() - old.dissolved_cation().value(),
            f64::from(source.cation_stoichiometry()) * solid,
        )?;
        for ((id, total), (a, b)) in old
            .family_totals()
            .iter()
            .zip(old.aqueous().families.iter().zip(&new.aqueous().families))
        {
            if id != &a.id || id != &b.id {
                return Err(OsmoticNetworkError::Invalid(
                    "reactive ionic family axis changed",
                ));
            }
            let solid_delta = if source.family_id() == Some(id.as_str()) {
                f64::from(source.anion_stoichiometry()) * solid
            } else {
                0.0
            };
            record(
                ElectrolyteInventory::Family(id.clone()),
                total.value(),
                b.total.value() - a.total.value(),
                solid_delta,
            )?;
        }
        for spectator in old.spectators() {
            // The immutable ledger equality above is source retention, not an
            // independent check of actual aqueous/publication inventories.
            let a = published_spectator(old, &spectator.species, spectator.charge)?;
            let b = published_spectator(new, &spectator.species, spectator.charge)?;
            record(
                ElectrolyteInventory::Spectator {
                    species: spectator.species.clone(),
                    charge: spectator.charge,
                },
                spectator.amount.value(),
                b - a,
                0.0,
            )?;
        }
        result.push(ElectrolyteCompartmentBalance {
            compartment,
            inventories,
            solid_formula_change_mol: solid,
            charge_balance_residual_mol: new.aqueous().charge_balance_residual.value(),
        });
    }
    Ok(result)
}

fn published_spectator(
    material: &KineticElectrolyteSolubilitySolution,
    species: &SpeciesId,
    charge: i32,
) -> Result<f64, OsmoticNetworkError> {
    let actual = material
        .aqueous()
        .spectators
        .iter()
        .find(|ion| &ion.species == species && ion.charge == charge)
        .ok_or(OsmoticNetworkError::Invalid(
            "retained spectator missing from original aqueous state",
        ))?
        .amount
        .value();
    let index = material
        .species()
        .iter()
        .position(|id| id == species)
        .ok_or(OsmoticNetworkError::Invalid(
            "retained spectator missing from published species axis",
        ))?;
    if material.solution().solute_moles()[index] != actual {
        return Err(OsmoticNetworkError::Invalid(
            "published spectator differs from original aqueous amount",
        ));
    }
    Ok(actual)
}

#[cfg(test)]
mod tests;
