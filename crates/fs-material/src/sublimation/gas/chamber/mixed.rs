//! Solid-vapor exchange from an explicitly resolved solid/liquid inventory.
//!
//! The original equilibrium enthalpy chart owns bulk storage and fusion. The
//! original sublimation law owns the SOLID-to-vapor latent difference. These
//! are different quantities: h_v-h_bulk is not the sublimation latent heat.
//! Gross emission spends only the OLD solid inventory, never the liquid or
//! material frozen later in the interval. Incoming deposition and reflected
//! heat enter the same bulk caloric endpoint; fusion energy is not added again.
//!
//! This reduction admits one pure substance with a sharp isothermal melting
//! plateau and local solid/liquid equilibrium. The caller supplies actual
//! exposed SOLID area; mass fraction alone does not determine surface coverage.
//! Liquid evaporation, nucleation on an all-liquid body, alloy mushy ranges,
//! separated phase temperatures, moving geometry, pressure-dependent melting
//! and gas-volume work require other owners. The original molecular and fixed
//! gas-volume assumptions still apply. Caloric/source admission is not phase-
//! diagram validation: the caller must admit the solid/vapor coexistence regime.

use std::collections::BTreeMap;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use super::{GasChamberCandidate, GasChamberSolve, GasChamberSurface};
use super::super::{
    ClosedGasParcel, GasSpeciesSnapshot, GasSurfaceError, GasSurfaceSolve,
    PreparedSurface, SolidSurfaceState, SurfaceSnapshot, ThermochemicalGasModel,
    compatibility, invalid, phase_at, poll, positive, prepare_surface_with_endpoint,
    snapshot_surface,
};
use super::super::accommodation::{ReflectedHeatStep, ThermalAccommodation};

/// Borrowed original caloric source, admitted once without copying its knots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharpMeltingCaloric<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    solid_at_melting_h_j_kg: f64,
    liquid_at_melting_h_j_kg: f64,
    solid_at_melting_density_kg_m3: f64,
}
impl<'a> SharpMeltingCaloric<'a> {
    /// Admit a strictly increasing solid branch followed by ONE sharp fusion
    /// plateau with linear enthalpy/phase fraction. Liquid sensible knots may
    /// remain on the original source, but this solid-interface path cannot use
    /// enthalpy above the liquid end of the plateau. No mushy data are fitted.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, GasSurfaceError> {
        let knots = curve.knots();
        if knots.len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(invalid("sharp fusion exceeds its original caloric knot budget"));
        }
        let first_liquid = knots.iter().position(|k| k.liquid_mass_fraction > 0.0)
            .ok_or(invalid("sharp fusion requires an original solid/liquid caloric chart"))?;
        let end = knots.iter().position(|k| k.liquid_mass_fraction == 1.0)
            .ok_or(invalid("sharp fusion requires the original fully liquid endpoint"))?;
        if first_liquid == 0 {
            return Err(invalid("sharp fusion has no solid endpoint"));
        }
        let start = first_liquid - 1;
        let solid = knots[start];
        let liquid = knots[end];
        let fusion = liquid.specific_enthalpy_j_kg - solid.specific_enthalpy_j_kg;
        if !positive(fusion) || solid.temperature_k != liquid.temperature_k {
            return Err(invalid("nonisothermal mushy intervals need independent constituent calorics"));
        }
        for pair in knots[..=start].windows(2) {
            if pair[0].liquid_mass_fraction != 0.0
                || pair[1].liquid_mass_fraction != 0.0
                || pair[0].temperature_k >= pair[1].temperature_k {
                return Err(invalid("solid sensible branch must be single phase and strictly increasing"));
            }
        }
        for knot in &knots[start..=end] {
            let fraction = (knot.specific_enthalpy_j_kg - solid.specific_enthalpy_j_kg) / fusion;
            if knot.temperature_k != solid.temperature_k
                || (fraction - knot.liquid_mass_fraction).abs() > 64.0 * f64::EPSILON {
                return Err(invalid("fusion knots must retain one temperature and a single latent lever rule"));
            }
        }
        Ok(Self {
            curve,
            solid_at_melting_h_j_kg: solid.specific_enthalpy_j_kg,
            liquid_at_melting_h_j_kg: liquid.specific_enthalpy_j_kg,
            solid_at_melting_density_kg_m3: solid.bulk_density_kg_m3,
        })
    }

    #[must_use]
    pub const fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.curve }

    /// Original plateau energy difference, not a separately supplied heat term.
    #[must_use]
    pub fn fusion_latent_heat_j_kg(&self) -> f64 {
        self.liquid_at_melting_h_j_kg - self.solid_at_melting_h_j_kg
    }

    /// Original-source state, including its identity. Exact all-liquid completion
    /// is a valid endpoint; a subsequent solid-interface interval must refuse.
    pub fn resolve(&self, h_j_kg: f64) -> Result<EquilibriumPhaseState, GasSurfaceError> {
        if !h_j_kg.is_finite() || h_j_kg > self.liquid_at_melting_h_j_kg {
            return Err(invalid("solid-interface caloric state exceeded its sharp fusion endpoint"));
        }
        self.curve.state_at_specific_enthalpy(h_j_kg)
            .map_err(|_| invalid("original solid/liquid enthalpy chart refused state"))
    }

    fn solid_enthalpy(&self, state: EquilibriumPhaseState) -> f64 {
        if state.is_fully_solid() { state.specific_enthalpy_j_kg() }
        else { self.solid_at_melting_h_j_kg }
    }
    fn solid_density(&self, state: EquilibriumPhaseState) -> f64 {
        if state.is_fully_solid() { state.bulk_density_kg_m3() }
        else { self.solid_at_melting_density_kg_m3 }
    }
    fn admit_density(&self, state: EquilibriumPhaseState, density: f64) -> Result<(), GasSurfaceError> {
        if !positive(density) || (density / self.solid_density(state) - 1.0).abs() > 64.0 * f64::EPSILON {
            Err(invalid("sublimation density disagrees with the original SOLID caloric branch"))
        } else { Ok(()) }
    }
}

/// Source-resolved phase inventories in a single well-mixed condensed body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensedInventory {
    pub caloric: EquilibriumPhaseState,
    pub solid_mass_kg: f64,
    pub liquid_mass_kg: f64,
}
impl CondensedInventory {
    pub fn new(source: SharpMeltingCaloric<'_>, state: SolidSurfaceState) -> Result<Self, GasSurfaceError> {
        if !positive(state.mass_kg) {
            return Err(invalid("positive total condensed mass required"));
        }
        let caloric = source.resolve(state.specific_enthalpy_j_kg)?;
        let solid_mass_kg = state.mass_kg * caloric.solid_mass_fraction();
        let liquid_mass_kg = state.mass_kg * caloric.liquid_mass_fraction();
        if !solid_mass_kg.is_finite() || !liquid_mass_kg.is_finite()
            || (caloric.solid_mass_fraction() > 0.0 && solid_mass_kg == 0.0)
            || (caloric.liquid_mass_fraction() > 0.0 && liquid_mass_kg == 0.0) {
            return Err(invalid("condensed constituent inventory is unrepresentable"));
        }
        Ok(Self { caloric, solid_mass_kg, liquid_mass_kg })
    }
}

/// Fusion is the signed LIQUID inventory change, not another vapor mass stream.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensedPhaseStep {
    pub initial: CondensedInventory,
    pub final_inventory: CondensedInventory,
    /// Positive melting; negative freezing, after solid-vapor exchange.
    pub melted_mass_kg: f64,
    /// Delta m_s + mass_into_vapor + melted_mass, in kg.
    pub solid_mass_residual_kg: f64,
}

/// Existing joint endpoint plus explicit original-source constituent accounts.
/// In exchange().surfaces(), `solid` carries TOTAL condensed mass and bulk h.
#[derive(Clone, Debug, PartialEq)]
pub struct PartitionedChamberCandidate<'catalog> {
    exchange: GasChamberCandidate<'catalog>,
    phases: BTreeMap<u64, CondensedPhaseStep>,
    reflection: BTreeMap<u64, ReflectedHeatStep>,
}
impl<'catalog> PartitionedChamberCandidate<'catalog> {
    #[must_use]
    pub const fn exchange(&self) -> &GasChamberCandidate<'catalog> { &self.exchange }
    #[must_use]
    pub fn phases(&self) -> &BTreeMap<u64, CondensedPhaseStep> { &self.phases }
    #[must_use]
    pub fn reflection(&self) -> &BTreeMap<u64, ReflectedHeatStep> { &self.reflection }
    #[must_use]
    pub fn into_parts(self) -> (GasChamberCandidate<'catalog>, BTreeMap<u64, CondensedPhaseStep>, BTreeMap<u64, ReflectedHeatStep>) {
        (self.exchange, self.phases, self.reflection)
    }
}

pub(super) fn validate_roster(
    surfaces: &[GasChamberSurface<'_>], calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
) -> Result<(), GasSurfaceError> {
    if calorics.len() != surfaces.len() || surfaces.iter().any(|s| !calorics.contains_key(&s.id)) {
        return Err(invalid("sharp fusion sources must cover exactly the selected surfaces"));
    }
    Ok(())
}

// Shared by the scalar facade and the joint chamber. Only the ACTUAL old solid
// inventory and constituent h enter this original molecular snapshot.
pub(in crate::sublimation::gas) fn snapshot(
    surface: &GasChamberSurface<'_>, caloric: SharpMeltingCaloric<'_>, gas: GasSpeciesSnapshot,
    catalog: &ThermochemicalGasModel, dt: f64, controls: GasSurfaceSolve,
) -> Result<SurfaceSnapshot, GasSurfaceError> {
    let inventory = CondensedInventory::new(caloric, surface.solid)?;
    if inventory.solid_mass_kg == 0.0 {
        return Err(invalid("no solid interface remains; evaporation or nucleation needs its own owner"));
    }
    let temperature = inventory.caloric.temperature_k();
    if (surface.temperature_from_enthalpy)(surface.solid.specific_enthalpy_j_kg)? != temperature {
        return Err(invalid("surface temperature callback disagrees with its original phase chart"));
    }
    let phase = phase_at(surface.law, temperature, gas.molar_mass)?;
    caloric.admit_density(inventory.caloric, phase.solid_density_kg_m3)?;
    snapshot_surface(surface.law, surface.binding, SolidSurfaceState {
        mass_kg: inventory.solid_mass_kg,
        specific_enthalpy_j_kg: caloric.solid_enthalpy(inventory.caloric),
        area_m2: surface.solid.area_m2,
    }, gas, catalog, &|_| Ok(temperature), dt, controls)
}

// One extensive update shared with the pure-solid owner. Only the caloric
// endpoint admission differs; no duplicate storage or fusion-energy arithmetic.
pub(in crate::sublimation::gas) fn prepare(
    surface: &GasChamberSurface<'_>, caloric: SharpMeltingCaloric<'_>, mut snapshot: SurfaceSnapshot,
    catalog: &ThermochemicalGasModel, controls: GasSurfaceSolve, reflected_heat: f64,
) -> Result<PreparedSurface, GasSurfaceError> {
    snapshot.solid = surface.solid;
    prepare_surface_with_endpoint(surface.binding, snapshot, controls, reflected_heat, |h| {
        let final_state = caloric.resolve(h)?;
        let temperature = final_state.temperature_k();
        if (surface.temperature_from_enthalpy)(h)? != temperature {
            return Err(invalid("endpoint callback disagrees with the original phase chart"));
        }
        let phase = phase_at(surface.law, temperature, snapshot.gas.molar_mass)?;
        caloric.admit_density(final_state, phase.solid_density_kg_m3)?;
        let vapor_h = catalog.species_specific_enthalpy_j_kg(&surface.binding.species, temperature)?;
        compatibility(vapor_h, caloric.solid_enthalpy(final_state),
            surface.binding.solid_enthalpy_offset_j_kg, phase.latent_heat_j_kg,
            temperature, controls.latent_compatibility_tolerance_j_kg)?;
        Ok((temperature, phase))
    })
}

pub(super) fn phase_reports(
    surfaces: &[GasChamberSurface<'_>], calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    endpoint: &GasChamberCandidate<'_>, mass_tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<u64, CondensedPhaseStep>, GasSurfaceError> {
    let originals: BTreeMap<_, _> = surfaces.iter().map(|s| (s.id, s)).collect();
    let mut reports = BTreeMap::new();
    for row in endpoint.surfaces() {
        poll(cancelled)?;
        let original = originals[&row.id];
        let source = calorics[&row.id];
        let initial = CondensedInventory::new(source, original.solid)?;
        let final_inventory = CondensedInventory::new(source, row.solid)?;
        let melted_mass_kg = final_inventory.liquid_mass_kg - initial.liquid_mass_kg;
        // Use actual bulk inventory loss, not only the last interval's flux:
        // this routine is also valid for whole adaptive-history endpoints.
        let loss = original.solid.mass_kg - row.solid.mass_kg;
        let residual = (final_inventory.solid_mass_kg - initial.solid_mass_kg) + loss + melted_mass_kg;
        if !residual.is_finite() || residual.abs() > mass_tolerance {
            return Err(invalid("solid/liquid/vapor constituent inventory balance failed"));
        }
        reports.insert(row.id, CondensedPhaseStep { initial, final_inventory, melted_mass_kg, solid_mass_residual_kg: residual });
    }
    Ok(reports)
}

/// One jointly reconstructed gas endpoint for sharply melting condensed bodies.
///
/// `surface.solid.mass_kg` is TOTAL condensed mass and its h is BULK enthalpy.
/// Its callback must evaluate the same retained chart. Area is exposed solid
/// area, not total geometric area or solid mass fraction times total area.
/// The original chamber owns canonical species aggregation, old gas inventory
/// limits, reflected heat, one gas reconstruction and global conservation.
/// An exact all-liquid endpoint is returned explicitly; a later call refuses
/// until an evaporation/interface-topology owner takes over. No partial state
/// is published on failure. Refinement is required for this explicit interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_partitioned_chamber<'catalog>(
    surfaces: &[GasChamberSurface<'_>], calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>, gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: GasChamberSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PartitionedChamberCandidate<'catalog>, GasSurfaceError> {
    let (exchange, reflection) = super::advance_chamber(surfaces, accommodation, Some(calorics),
        gas, catalog, duration_s, controls, &mut cancelled)?;
    let phases = phase_reports(surfaces, calorics, &exchange, controls.surface.gas.mass_tolerance_kg, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(PartitionedChamberCandidate { exchange, phases, reflection })
}

#[cfg(test)]
mod tests;
