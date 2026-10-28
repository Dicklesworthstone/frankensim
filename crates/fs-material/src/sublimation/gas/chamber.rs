//! Joint nonisothermal molecular exchange with one original gas parcel.
//!
//! Every surface sees the same initial gas temperature, composition and
//! species partial pressure. Gross incoming streams of each species share its
//! OLD inventory budget; outgoing streams are bounded by their own solids.
//! Solid endpoints use the same preparation as the single-surface operation.
//! All signed species/energy streams are then aggregated canonically and the
//! original gas owner reconstructs ONE endpoint. Sequential surface updates
//! would introduce order dependence and spend the same vapor more than once.
//!
//! Surfaces may use different explicit species, phase laws and caloric charts.
//! Each keeps its fixed datum offset and must satisfy latent compatibility at
//! both endpoints. The local molecular assumptions in the parent module apply.
//! The opt-in accommodation path adds reflected heat from the same old gas to
//! the same solid endpoints and performs one combined gas reconstruction.

/// Bounded adaptive histories using the same conservative joint interval.
pub mod advance;
/// Explicit solid/liquid inventories on original sharp-melting caloric charts.
pub mod mixed;

use super::{
    ClosedGasParcel, GasSpeciesBoundaryCandidate, GasSpeciesSnapshot, GasSurfaceBinding,
    GasSurfaceError, GasSurfaceSolve, PreparedSurface, SolidSurfaceState, SpeciesBoundaryRelativeTransfer,
    SpeciesId, SublimationError, SublimationLaw, ThermochemicalGasModel, invalid, inventory_limit,
    poll, prepare_surface_with_heat, snapshot_surface, species_snapshot, validate_controls,
};
use std::collections::BTreeMap;
use super::accommodation::{ReflectedHeatStep, ThermalAccommodation};

/// One caller-owned physical surface and its original constitutive charts.
/// Identity is explicit and stable across intervals; input order is immaterial.
/// On the opt-in mixed path, `solid` stores total condensed mass and bulk h;
/// the separate admitted phase chart owns its actual solid/liquid partition.
#[derive(Clone, Copy)]
pub struct GasChamberSurface<'a> {
    pub id: u64,
    pub law: &'a dyn SublimationLaw,
    pub binding: &'a GasSurfaceBinding,
    pub solid: SolidSurfaceState,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
}

/// The original single-surface controls and a caller-owned work/allocation cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasChamberSolve {
    pub surface: GasSurfaceSolve,
    pub maximum_surfaces: usize,
}

/// One accepted solid and its allocated stream, in canonical surface-id order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasChamberSurfaceStep {
    pub id: u64,
    pub solid: SolidSurfaceState,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub surface_temperature_k: [f64; 2],
    pub outgoing_specific_energy_j_kg: f64,
    pub incoming_specific_energy_j_kg: Option<f64>,
    pub energy_into_gas_j: f64,
    /// Original-datum report: intended transferred-mass reference energy plus
    /// actual represented caloric storage. Inventory rounding is independently
    /// bounded. On the mixed path this includes the complete condensed body.
    pub solid_energy_change_j: f64,
    /// Solid storage change plus its allocated stream energy, with conserved
    /// formation references canceled before evaluation. The actual gas storage
    /// residual is accounted for only in the joint chamber report.
    pub energy_residual_j: f64,
    /// Phase-transfer retreat only; mixed-path fusion is reported separately.
    pub recession_m: f64,
}
impl GasChamberSurfaceStep {
    fn new(id: u64, surface: PreparedSurface) -> Self {
        Self {
            id,
            solid: surface.solid,
            outgoing_mass_kg: surface.outgoing_mass_kg,
            incoming_mass_kg: surface.incoming_mass_kg,
            transferred_mass_kg: surface.transferred_mass_kg,
            surface_temperature_k: surface.surface_temperature_k,
            outgoing_specific_energy_j_kg: surface.outgoing_specific_energy_j_kg,
            incoming_specific_energy_j_kg: surface.incoming_specific_energy_j_kg,
            energy_into_gas_j: surface.energy_into_gas_j,
            solid_energy_change_j: surface.solid_energy_change_j,
            energy_residual_j: surface.energy_residual_j,
            recession_m: surface.recession_m,
        }
    }
}

/// One exchanging gas species, canonicalized by the original SpeciesId.
/// A carrier may exchange reflected heat with exactly zero phase mass streams.
#[derive(Clone, Debug, PartialEq)]
pub struct GasChamberSpeciesStep {
    pub species: SpeciesId,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub energy_into_gas_j: f64,
    pub vapor_mass_kg: [f64; 2],
    /// Actual gas species change plus actual same-species solid mass changes.
    pub mass_residual_kg: f64,
}

/// Global accounts use one actual accepted gas endpoint and every accepted solid.
#[derive(Clone, Debug, PartialEq)]
pub struct GasChamberStep {
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub energy_into_gas_j: f64,
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    /// Sum of actual local thermal balance defects after conserved formation
    /// terms cancel; subtracting the displayed absolute reports may lose heat.
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub gas_temperature_k: [f64; 2],
    pub species: Vec<GasChamberSpeciesStep>,
}

/// Publish all solids and the one gas together; no trial mutates accepted input.
#[derive(Clone, Debug, PartialEq)]
pub struct GasChamberCandidate<'catalog> {
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    surfaces: Vec<GasChamberSurfaceStep>,
    step: GasChamberStep,
}
impl<'catalog> GasChamberCandidate<'catalog> {
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    #[must_use]
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    #[must_use]
    pub fn surfaces(&self) -> &[GasChamberSurfaceStep] {
        &self.surfaces
    }
    #[must_use]
    pub const fn step(&self) -> &GasChamberStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasChamberSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasChamberStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.surfaces, gas, catalog, self.step)
    }
}

/// Deterministic compensated accumulation, using the shared arithmetic owner.
#[derive(Default)]
struct Total {
    sum: f64,
    correction: f64,
}
impl Total {
    fn add(&mut self, term: f64) -> Result<(), GasSurfaceError> {
        let (sum, error) = fs_math::eft::two_sum(self.sum, term);
        let correction = self.correction + error;
        if ![term, sum, error, correction]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(invalid("joint molecular accounting is unrepresentable"));
        }
        self.sum = sum;
        self.correction = correction;
        Ok(())
    }
    fn value(&self) -> Result<f64, GasSurfaceError> {
        let value = self.sum + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("joint molecular total is unrepresentable"))
        }
    }
}

struct SpeciesAccount {
    gas: GasSpeciesSnapshot,
    outgoing: Total,
    incoming: Total,
    transferred: Total,
    solid_mass_change: Total,
    reference_energy_j_kg: f64,
    relative_energy: Total,
    phase_relative_energy: Total,
    reflected_heat: Total,
    active_streams: usize,
    single_stream_energy_j: f64,
}
impl SpeciesAccount {
    fn new(gas: GasSpeciesSnapshot, outgoing_energy: f64) -> Self {
        Self {
            gas,
            outgoing: Total::default(),
            incoming: Total::default(),
            transferred: Total::default(),
            solid_mass_change: Total::default(),
            reference_energy_j_kg: gas.incoming_specific_energy_j_kg.unwrap_or(outgoing_energy),
            relative_energy: Total::default(),
            phase_relative_energy: Total::default(),
            reflected_heat: Total::default(),
            active_streams: 0,
            single_stream_energy_j: 0.0,
        }
    }
    fn energy(&self) -> Result<f64, GasSurfaceError> {
        // Preserve the exact single-surface operation tree. For real multi-
        // surface exchange, transport the common formation reference only with
        // NET species mass; accumulate opposing sensible contributions separately.
        let value = if self.active_streams <= 1 {
            self.single_stream_energy_j
        } else {
            self.transferred.value()? * self.reference_energy_j_kg + self.relative_energy.value()?
        };
        let reflected = self.reflected_heat.value()?;
        let value = if reflected == 0.0 { value } else { value + reflected };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("joint species stream energy is unrepresentable"))
        }
    }
}

fn merge_limit(limit: &mut Option<f64>, candidate: Option<f64>) {
    if let Some(candidate) = candidate {
        *limit = Some(limit.map_or(candidate, |previous| previous.min(candidate)));
    }
}

/// Advance heterogeneous surfaces jointly against one evolving gas snapshot.
///
/// Rates freeze at the initial temperatures and partial pressures. Both gross
/// outgoing/solid and AGGREGATE incoming/species inventory bounds must hold.
/// The returned InventoryLimit supplies a conservative duration for retrying
/// the whole interval; no component is clipped or partially published.
///
/// Inputs are canonicalized by stable surface id. Duplicate ids, excessive work,
/// conflicting sources, failed endpoint charts/balances, and cancellation refuse
/// every proposed endpoint. An empty active surface set preserves the gas.
pub fn advance_chamber_with_gas<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasChamberCandidate<'catalog>, GasSurfaceError> {
    advance_chamber(surfaces, None, None, gas, catalog, duration_s, controls, cancelled)
        .map(|(candidate, _)| candidate)
}

/// Joint phase exchange and per-surface reflected molecular heat.
///
/// The accommodation map must cover EXACTLY the active surface identities.
/// Supply explicit zero species coefficients for a specular surface; absent
/// entries never become guessed thermal coefficients. All rates sample the
/// same original gas, not successive updated gases. Each surface's reflected
/// subaccount is included in its total energy transfer and must not be added
/// again. Carrier heat does not consume or replenish species mass allowances.
///
/// Returns the original joint candidate and canonical surface reflection
/// accounts. A failure in any surface or final gas rejects the entire interval.
/// This explicit step still requires physical timestep refinement; the existing
/// adaptive and heat-timeline entry points retain their specular default.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_chamber_with_accommodation<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    accommodation: &BTreeMap<u64, ThermalAccommodation>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<(GasChamberCandidate<'catalog>, BTreeMap<u64, ReflectedHeatStep>), GasSurfaceError> {
    advance_chamber(surfaces, Some(accommodation), None, gas, catalog, duration_s, controls, cancelled)
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_chamber<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    calorics: Option<&BTreeMap<u64, mixed::SharpMeltingCaloric<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(GasChamberCandidate<'catalog>, BTreeMap<u64, ReflectedHeatStep>), GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls.surface, duration_s)?;
    if controls.maximum_surfaces == 0 || surfaces.len() > controls.maximum_surfaces {
        return Err(invalid(
            "joint molecular surface count exceeds its declared work budget",
        ));
    }
    let mut ordered = BTreeMap::new();
    for surface in surfaces {
        poll(&mut cancelled)?;
        if ordered.insert(surface.id, surface).is_some() {
            return Err(invalid(
                "joint molecular surfaces require unique stable identities",
            ));
        }
    }
    if let Some(laws) = accommodation {
        if laws.len() != ordered.len() || ordered.keys().any(|id| !laws.contains_key(id)) {
            return Err(invalid("thermal accommodation must cover exactly the active chamber surfaces"));
        }
    }
    if let Some(laws) = calorics { mixed::validate_roster(surfaces, laws)?; }
    let mut reflections = BTreeMap::new();
    let mut species: BTreeMap<SpeciesId, SpeciesAccount> = BTreeMap::new();
    let mut snapshots = Vec::with_capacity(ordered.len());
    let mut limit = None;
    for surface in ordered.values() {
        poll(&mut cancelled)?;
        let id = &surface.binding.species;
        let context = match species.get(id) {
            Some(account) => account.gas,
            None => species_snapshot(gas, catalog, id)?,
        };
        let snapshot = match calorics {
            Some(laws) => mixed::snapshot(surface, laws[&surface.id], context, catalog, duration_s, controls.surface)?,
            None => snapshot_surface(
                surface.law,
                surface.binding,
                surface.solid,
                context,
                catalog,
                surface.temperature_from_enthalpy,
                duration_s,
                controls.surface,
            )?,
        };
        let account = species.entry(id.clone()).or_insert_with(|| {
            SpeciesAccount::new(context, snapshot.outgoing_specific_energy_j_kg)
        });
        account.outgoing.add(snapshot.outgoing_mass_kg)?;
        account.incoming.add(snapshot.incoming_mass_kg)?;
        account.transferred.add(snapshot.transferred_mass_kg)?;
        merge_limit(
            &mut limit,
            inventory_limit(
                duration_s,
                snapshot.outgoing_mass_kg,
                // In a mixture this is only its OLD solid inventory.
                snapshot.solid.mass_kg,
                controls.surface.maximum_gross_inventory_fraction,
            ),
        );
        snapshots.push((*surface, snapshot));
    }
    for account in species.values() {
        poll(&mut cancelled)?;
        merge_limit(
            &mut limit,
            inventory_limit(
                duration_s,
                account.incoming.value()?,
                account.gas.vapor_mass_kg,
                controls.surface.maximum_gross_inventory_fraction,
            ),
        );
    }
    if let Some(maximum_duration_s) = limit {
        return Err(GasSurfaceError::InventoryLimit { maximum_duration_s });
    }
    let mut outputs = Vec::with_capacity(snapshots.len());
    let mut solid_energy = Total::default();
    let mut solid_residual = Total::default();
    let mut solid_mass_change = Total::default();
    for (surface, snapshot) in snapshots {
        poll(&mut cancelled)?;
        let mut reflected_heat = 0.0;
        if let Some(laws) = accommodation {
            let reflection = laws[&surface.id].evaluate(
                surface.binding, gas, snapshot.temperature_k, surface.solid.area_m2,
                duration_s, &mut cancelled,
            )?;
            reflected_heat = reflection.heat_into_gas_j;
            for row in &reflection.species {
                if row.heat_into_gas_j == 0.0 { continue; }
                // Reflection transfers energy, not mass. Create only the needed
                // carrier account, retaining the original sources and inventory.
                if !species.contains_key(&row.species) {
                    let context = species_snapshot(gas, catalog, &row.species)?;
                    species.insert(row.species.clone(), SpeciesAccount::new(context, 0.0));
                }
                species.get_mut(&row.species).expect("original reflected species account")
                    .reflected_heat.add(row.heat_into_gas_j)?;
            }
            reflections.insert(surface.id, reflection);
        }
        let prepared = match calorics {
            Some(laws) => mixed::prepare(surface, laws[&surface.id], snapshot, catalog, controls.surface, reflected_heat)?,
            None => prepare_surface_with_heat(
                surface.law,
                surface.binding,
                snapshot,
                catalog,
                surface.temperature_from_enthalpy,
                controls.surface,
                reflected_heat,
            )?,
        };
        let account = species
            .get_mut(&surface.binding.species)
            .expect("snapshot retains its original species account");
        let dm = prepared.solid.mass_kg - surface.solid.mass_kg;
        account.solid_mass_change.add(dm)?;
        solid_mass_change.add(dm)?;
        solid_energy.add(prepared.solid_energy_change_j)?;
        solid_residual.add(prepared.energy_residual_j)?;
        account.phase_relative_energy.add(prepared.phase_relative_energy_into_gas_j)?;
        let difference = snapshot.partial_enthalpy_change_j_kg
            - 0.5 * snapshot.gas.specific_gas_constant
                * (snapshot.temperature_k - snapshot.gas.temperature_k);
        account.relative_energy.add(
            prepared.transferred_mass_kg
                * (prepared.outgoing_specific_energy_j_kg - account.reference_energy_j_kg)
                + prepared.incoming_mass_kg * difference,
        )?;
        // Reconstruct phase-only energy directly. Subtracting reflected heat
        // from a large total formation energy can erase its sensible increment.
        let phase_energy = prepared.transferred_mass_kg * prepared.outgoing_specific_energy_j_kg
            + prepared.incoming_mass_kg * difference;
        if prepared.transferred_mass_kg != 0.0 || phase_energy != 0.0 {
            account.active_streams += 1;
            account.single_stream_energy_j = phase_energy;
        }
        outputs.push(GasChamberSurfaceStep::new(surface.id, prepared));
    }
    let mut rows = Vec::with_capacity(species.len());
    let mut total_outgoing = Total::default();
    let mut total_incoming = Total::default();
    let mut total_transferred = Total::default();
    let mut total_energy = Total::default();
    for (id, account) in &species {
        poll(&mut cancelled)?;
        let mass = account.transferred.value()?;
        let energy = account.energy()?;
        total_outgoing.add(account.outgoing.value()?)?;
        total_incoming.add(account.incoming.value()?)?;
        total_transferred.add(mass)?;
        total_energy.add(energy)?;
        rows.push(SpeciesBoundaryRelativeTransfer {
            species: id.clone(),
            mass_into_gas_kg: mass,
            relative_energy_into_gas_j: account.phase_relative_energy.value()?
                + account.reflected_heat.value()?,
        });
    }
    // The ONLY gas reconstruction in this joint interval.
    let endpoint =
        gas.exchange_species_boundary_relative(catalog, &rows, controls.surface.gas, &mut cancelled)?;
    let solid_energy_change_j = solid_energy.value()?;
    let gas_energy_change_j = endpoint.step().stored_energy_change_j;
    let energy_residual_j = solid_residual.value()? + endpoint.step().energy_residual_j;
    if !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.surface.gas.energy_tolerance_j
    {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual_j,
            tolerance_j: controls.surface.gas.energy_tolerance_j,
        });
    }
    let mass_residual_kg = solid_mass_change.value()? + (endpoint.gas().mass_kg() - gas.mass_kg());
    if !mass_residual_kg.is_finite()
        || mass_residual_kg.abs() > controls.surface.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "joint molecular endpoints failed total mass balance",
        ));
    }
    let mut species_reports = Vec::with_capacity(species.len());
    for (id, account) in species {
        poll(&mut cancelled)?;
        let final_mass = endpoint.gas().species_mass_kg(&id)?;
        let mass_residual_kg =
            account.solid_mass_change.value()? + (final_mass - account.gas.vapor_mass_kg);
        if !mass_residual_kg.is_finite()
            || mass_residual_kg.abs() > controls.surface.gas.mass_tolerance_kg
        {
            return Err(invalid(
                "joint molecular endpoints failed individual species mass balance",
            ));
        }
        species_reports.push(GasChamberSpeciesStep {
            species: id,
            outgoing_mass_kg: account.outgoing.value()?,
            incoming_mass_kg: account.incoming.value()?,
            transferred_mass_kg: account.transferred.value()?,
            energy_into_gas_j: account.energy()?,
            vapor_mass_kg: [account.gas.vapor_mass_kg, final_mass],
            mass_residual_kg,
        });
    }
    let final_temperature = endpoint.gas().state().state().temperature;
    if let Some(laws) = accommodation {
        for surface in &outputs {
            poll(&mut cancelled)?;
            laws[&surface.id].admit_temperatures(surface.surface_temperature_k[1], final_temperature)?;
        }
    }
    poll(&mut cancelled)?;
    Ok((GasChamberCandidate {
        gas: endpoint,
        surfaces: outputs,
        step: GasChamberStep {
            outgoing_mass_kg: total_outgoing.value()?,
            incoming_mass_kg: total_incoming.value()?,
            transferred_mass_kg: total_transferred.value()?,
            energy_into_gas_j: total_energy.value()?,
            solid_energy_change_j,
            gas_energy_change_j,
            energy_residual_j,
            mass_residual_kg,
            gas_temperature_k: [gas.state().state().temperature, final_temperature],
            species: species_reports,
        },
    }, reflections))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "chamber/accommodation_tests.rs"]
mod accommodation_tests;
