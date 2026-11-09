//! Source-derived molecular films refreshed inside the original chamber history.
//!
//! One selected binary diffusivity law is retained for each bound planar strip.
//! Each physical trial re-evaluates that SAME law at the declared representative
//! temperature and the current gas's total pressure. The original parcel model
//! supplies transport at its current composition; no carrier preset or Lewis
//! number substitutes for these sources. The existing material chamber owns
//! mass transfer, heat, time adaptation and rollback.
//!
//! The prescribed strip remains measured from its physical leading edge.
//! Surface area belongs to the chamber/growing-geometry owner. This adapter
//! does not infer a new streamwise length from area, resolve velocity or add
//! another gas-side heat source. The original dilute, weak-transfer, binary,
//! attached laminar continuum assumptions remain caller declarations.
//! Both growing wrappers additionally require each film-bound layer to have
//! the original planar morphology; whole spherical/cylindrical layers refuse.

use super::laminar_plate::LaminarPlateGeometry;
use super::{
    BinaryDiffusionMaterial, GasFilmError, LaminarPlateFilmConditions,
    ResolvedMaterialStatePoint, SpeciesId, Temperature, Pressure, Length, typed,
};
use fs_material::sublimation::gas::accommodation::ThermalAccommodation;
use fs_material::sublimation::gas::chamber::coating::resistance::{
    GasResistiveSurface,
    advance::{GasResistiveAdvanceCandidate, GasResistiveAdvanceSolve},
    driven::GasResistivePower,
    film::{GasResistiveFilmError, advance_resistive_chamber_adaptive_resolved},
    growing::{
        GasGrowingCoating, GasGrowingCoatingCandidate, GasGrowingCoatingSolve,
        GasIntegratedGrowingCoating, advance_growing_integrated_chamber_resolved,
        advance_growing_resistive_chamber_resolved,
    },
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::ClosedGasParcel,
};
use fs_qty::semantic::QuantityKind;
use fs_thermal::sublimation::coating::resistance::geometry::ConductiveCoatingGeometry;
use fs_thermal::sublimation::SublimationError;
use std::collections::BTreeMap;

/// Conservative exclusion of high vapor loading from the dilute plate model.
/// This five-percent ceiling is an implementation admission policy, not a
/// published universal error bound or proof of weak Stefan/blowing effects.
pub const MAX_DILUTE_VAPOR_MOLE_FRACTION: f64 = 0.05;

/// Explicit property-evaluation convention, frozen only within one physical
/// trial. Both choices retain the actual bulk total pressure and composition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChamberFilmTemperature {
    Bulk,
    /// Arithmetic mean of current deposit and current bulk gas temperatures.
    MeanSurfaceAndBulk,
}

/// Original binary D(T,p) source and physical planar strip for one surface.
#[derive(Clone, Debug)]
pub struct LaminarPlateFilmBinding<'m> {
    pub source: BinaryDiffusionMaterial<'m>,
    pub geometry: LaminarPlateGeometry,
    pub property_temperature: ChamberFilmTemperature,
    /// Authored ceiling in (0, MAX_DILUTE_VAPOR_MOLE_FRACTION]. Checked against
    /// actual bulk vapor, saturation/total pressure and the original
    /// two-temperature kinetic-equilibrium pressure/total pressure at every
    /// trial, including zero-time admission. This does not establish continuum,
    /// thin-boundary-layer or weak normal-flow validity by itself.
    pub maximum_vapor_mole_fraction: f64,
}

#[derive(Debug)]
struct SelectedFilm<'m> {
    id: u64,
    vapor: SpeciesId,
    binding: LaminarPlateFilmBinding<'m>,
    selected: ResolvedMaterialStatePoint,
}

/// Prepared source selection for a bounded, possibly sparse surface roster.
///
/// Preparation pins the actual diffusivity claim once. Rejected trials never
/// change that selection. The original current gas model continues to own its
/// thermodynamic and pure-species transport sources. A different active carrier
/// or a third gas species refuses instead of relabeling a binary film.
///
/// The binary source context may contain only canonical temperature/pressure
/// coordinates. Additional composition, humidity or other state coordinates
/// cannot remain silently frozen as the gas changes; they need their own
/// explicitly state-dependent source binding.
#[derive(Debug)]
pub struct SourceLaminarPlateChamberFilms<'m> {
    films: Vec<SelectedFilm<'m>>,
    maximum_surfaces: usize,
}

impl<'m> SourceLaminarPlateChamberFilms<'m> {
    /// Prepare each original source at the actual initial gas/deposit state.
    ///
    /// The maximum_surfaces bound applies to preparation and later resolution.
    /// The sparse binding map references existing, uniquely identified rows.
    /// A carrier-only initial gas is admitted: zero vapor is not replaced by
    /// an artificial trace. Pure vapor and foreign active species refuse.
    pub fn new(
        bindings: &BTreeMap<u64, LaminarPlateFilmBinding<'m>>,
        surfaces: &[GasResistiveSurface<'_>],
        gas: &ClosedGasParcel,
        maximum_surfaces: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasFilmError> {
        poll(&mut cancelled)?;
        if surfaces.len() > maximum_surfaces || bindings.len() > surfaces.len() {
            return Err(GasFilmError::Invalid(
                "live gas-film roster exceeds its declared surface bound",
            ));
        }
        let mut rows = BTreeMap::new();
        for surface in surfaces {
            poll(&mut cancelled)?;
            if rows.insert(surface.id, surface).is_some() {
                return Err(GasFilmError::Invalid("duplicate gas-film surface identity"));
            }
        }
        let mut films = Vec::with_capacity(bindings.len());
        for (&id, binding) in bindings {
            poll(&mut cancelled)?;
            let surface = rows.get(&id).ok_or(GasFilmError::Invalid(
                "gas-film binding references an unknown surface",
            ))?;
            if binding.source.context.axes().keys().any(|axis| {
                axis != "temperature" && axis != "pressure"
            }) {
                return Err(GasFilmError::Invalid(
                    "live binary diffusivity requires only canonical temperature/pressure context axes",
                ));
            }
            let (_, selected) = evaluate(
                binding, surface, gas, None, &mut cancelled,
            )?;
            films.push(SelectedFilm {
                id,
                vapor: surface.binding.species.clone(),
                binding: binding.clone(),
                selected,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self { films, maximum_surfaces })
    }

    /// Original physical strip and source binding; no material-name inference.
    pub fn binding(&self, id: u64) -> Option<&LaminarPlateFilmBinding<'m>> {
        self.films.binary_search_by_key(&id, |film| film.id)
            .ok().map(|index| &self.films[index].binding)
    }

    /// Actual selected diffusivity law at preparation, retained across trials.
    pub fn selected_diffusivity(&self, id: u64) -> Option<&ResolvedMaterialStatePoint> {
        self.films.binary_search_by_key(&id, |film| film.id)
            .ok().map(|index| &self.films[index].selected)
    }

    /// Fill reusable coefficients from current canonical rows and current gas.
    ///
    /// Rows must be in strictly increasing ID order, as supplied by the original
    /// adaptive owner. Every bound entry is overwritten, including zero-duration
    /// admission. The history wrappers leave unbound rows free-molecular.
    /// No physical state or prepared source is mutated. The coefficient map is
    /// trial scratch; an error may leave it partially filled and it must then
    /// be discarded. Material wrappers already enforce that transaction.
    pub fn resolve(
        &self,
        surfaces: &[GasResistiveSurface<'_>],
        gas: &ClosedGasParcel,
        coefficients: &mut BTreeMap<u64, f64>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), GasFilmError> {
        poll(cancelled)?;
        if surfaces.len() > self.maximum_surfaces {
            return Err(GasFilmError::Invalid(
                "live gas-film roster exceeds its declared surface bound",
            ));
        }
        let mut previous = None;
        for surface in surfaces {
            poll(cancelled)?;
            if previous.is_some_and(|id| id >= surface.id) {
                return Err(GasFilmError::Invalid(
                    "live gas-film rows must retain canonical unique identities",
                ));
            }
            previous = Some(surface.id);
        }
        // One monotone walk through the already ordered surface roster. The
        // coefficient map reuses existing nodes across original physical trials.
        let mut rows = surfaces.iter().peekable();
        for film in &self.films {
            poll(cancelled)?;
            while rows.peek().is_some_and(|surface| surface.id < film.id) {
                poll(cancelled)?;
                rows.next();
            }
            let surface = rows.next().filter(|surface| surface.id == film.id)
                .ok_or(GasFilmError::Invalid(
                    "prepared gas-film surface is missing from the current roster",
                ))?;
            if surface.binding.species != film.vapor {
                return Err(GasFilmError::Invalid(
                    "prepared gas film cannot change its volatile species",
                ));
            }
            let (coefficient, _) = evaluate(
                &film.binding, surface, gas, Some(&film.selected), cancelled,
            )?;
            coefficients.insert(film.id, coefficient);
        }
        poll(cancelled)
    }

    fn admit_planar_layers(
        &self,
        mut geometry: impl FnMut(u64) -> Option<ConductiveCoatingGeometry>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), GasFilmError> {
        poll(cancelled)?;
        for film in &self.films {
            poll(cancelled)?;
            let geometry = geometry(film.id).ok_or(GasFilmError::Invalid(
                "prepared gas-film surface needs its growing layer model",
            ))?;
            // Structural identity against the original validated planar shape.
            // Equal exposed areas or rounded curvature tangents are insufficient.
            let planar = ConductiveCoatingGeometry::planar(
                geometry.geometry().density_kg_m3(), geometry.contact_area_m2(),
            ).map_err(GasFilmError::Surface)?;
            if geometry != planar {
                return Err(GasFilmError::Invalid(
                    "laminar flat-plate gas films require planar growing layers",
                ));
            }
        }
        Ok(())
    }

    /// Evolve fixed-geometry resistive coatings with live original gas films.
    /// The material owner retains its existing controller and energy accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_resistive<'catalog>(
        &self,
        surfaces: &[GasResistiveSurface<'_>],
        accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
        powers: &BTreeMap<u64, GasResistivePower<'_>>,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        duration_s: f64,
        controls: GasResistiveAdvanceSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveFilmError<GasFilmError>> {
        advance_resistive_chamber_adaptive_resolved(
            surfaces,
            accommodation,
            |rows, gas, coefficients, cancelled| self.resolve(rows, gas, coefficients, cancelled),
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Evolve growing geometry with an original enthalpy-chart conductivity law.
    /// Conductivity is sampled per interval by the existing frozen-law owner;
    /// the molecular film is independently refreshed from the current gas.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_growing_resistive<'catalog>(
        &self,
        surfaces: &[GasResistiveSurface<'_>],
        layers: &BTreeMap<u64, GasGrowingCoating<'_>>,
        accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
        powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        duration_s: f64,
        controls: GasGrowingCoatingSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasGrowingCoatingCandidate<'catalog>, GasResistiveFilmError<GasFilmError>> {
        self.admit_planar_layers(
            |id| layers.get(&id).map(|layer| layer.geometry), &mut cancelled,
        ).map_err(GasResistiveFilmError::Source)?;
        advance_growing_resistive_chamber_resolved(
            surfaces,
            layers,
            accommodation,
            |rows, gas, coefficients, cancelled| self.resolve(rows, gas, coefficients, cancelled),
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Evolve growing coatings, original nonlinear conductivity and live films.
    /// The material owner refreshes mass-derived geometry BEFORE this resolver.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_growing_integrated<'catalog>(
        &self,
        surfaces: &[GasResistiveSurface<'_>],
        layers: &BTreeMap<u64, GasIntegratedGrowingCoating<'_>>,
        accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
        powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        duration_s: f64,
        controls: GasGrowingCoatingSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasGrowingCoatingCandidate<'catalog>, GasResistiveFilmError<GasFilmError>> {
        self.admit_planar_layers(
            |id| layers.get(&id).map(|layer| layer.geometry), &mut cancelled,
        ).map_err(GasResistiveFilmError::Source)?;
        advance_growing_integrated_chamber_resolved(
            surfaces,
            layers,
            accommodation,
            |rows, gas, coefficients, cancelled| self.resolve(rows, gas, coefficients, cancelled),
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    }
}

fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), GasFilmError> {
    if cancelled() { Err(GasFilmError::Cancelled) } else { Ok(()) }
}

fn evaluate(
    binding: &LaminarPlateFilmBinding<'_>,
    surface: &GasResistiveSurface<'_>,
    gas: &ClosedGasParcel,
    previous: Option<&ResolvedMaterialStatePoint>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<(f64, ResolvedMaterialStatePoint), GasFilmError> {
    poll(cancelled)?;
    let ceiling = binding.maximum_vapor_mole_fraction;
    if !ceiling.is_finite() || ceiling <= 0.0 || ceiling > MAX_DILUTE_VAPOR_MOLE_FRACTION {
        return Err(GasFilmError::Invalid(
            "live dilute film needs an authored vapor-fraction ceiling in (0,0.05]",
        ));
    }
    let vapor = &surface.binding.species;
    let vapor_fraction = binary_vapor_fraction(&binding.source, vapor, gas, cancelled)?;
    let bulk = gas.state().state();
    let surface_temperature = (surface.deposit_temperature)(
        surface.state.deposit.specific_enthalpy_j_kg,
    ).map_err(GasFilmError::Surface)?;
    let [minimum_k, maximum_k] = surface.law.temperature_range_k();
    if !minimum_k.is_finite() || minimum_k <= 0.0
        || !maximum_k.is_finite() || maximum_k < minimum_k
    {
        return Err(GasFilmError::Invalid("invalid original coating-temperature domain"));
    }
    if !surface_temperature.is_finite()
        || surface_temperature < minimum_k || surface_temperature > maximum_k
    {
        return Err(GasFilmError::Surface(SublimationError::TemperatureOutOfRange {
            temperature_k: surface_temperature, minimum_k, maximum_k,
        }));
    }
    poll(cancelled)?;
    let phase = surface.law.properties(surface_temperature).map_err(GasFilmError::Surface)?;
    if !phase.saturation_pressure_pa.is_finite() || phase.saturation_pressure_pa < 0.0 {
        return Err(GasFilmError::Invalid("invalid original coating saturation pressure"));
    }
    let saturation_fraction = phase.saturation_pressure_pa / bulk.pressure;
    // The original two-temperature molecular owner balances at
    // p_v = p_sat(T_d) sqrt(T_g/T_d), not simply p_sat(T_d).
    let kinetic_fraction = saturation_fraction
        * fs_math::det::sqrt(bulk.temperature / surface_temperature);
    if [vapor_fraction, saturation_fraction, kinetic_fraction].into_iter()
        .any(|fraction| !fraction.is_finite() || fraction < 0.0 || fraction > ceiling)
    {
        return Err(GasFilmError::OutsideDilute {
            bulk_vapor_mole_fraction: vapor_fraction,
            saturation_pressure_fraction: saturation_fraction,
            kinetic_equilibrium_pressure_fraction: kinetic_fraction,
            maximum_vapor_mole_fraction: ceiling,
        });
    }
    evaluate_original_film(binding, vapor, surface_temperature, gas, previous, cancelled)
        .map(|(_, coefficient, selected)| (coefficient, selected))
}

// The liquid and solid adapters share chemical admission and original
// transport evaluation. Their different equilibrium/kinetic laws remain with
// their respective phase owners.
pub(super) fn binary_vapor_fraction(
    source: &BinaryDiffusionMaterial<'_>,
    vapor: &SpeciesId,
    gas: &ClosedGasParcel,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<f64, GasFilmError> {
    poll(cancelled)?;
    let pair = &source.species;
    if pair[0] == pair[1] || !pair.contains(vapor) {
        return Err(GasFilmError::Invalid(
            "surface vapor must belong to two distinct binary-film species",
        ));
    }
    let carrier = if &pair[0] == vapor { &pair[1] } else { &pair[0] };
    let species = gas.model().thermodynamics().components();
    let fractions = gas.model().thermodynamics().mole_composition().fractions();
    if fractions.len() != species.len() {
        return Err(GasFilmError::Invalid("gas composition and species roster disagree"));
    }
    let mut has_carrier = false;
    let mut vapor_fraction = 0.0;
    for (component, &fraction) in species.iter().zip(fractions) {
        poll(cancelled)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(GasFilmError::Invalid("invalid current gas mole fraction"));
        }
        // A catalog entry with zero inventory is not an active constituent.
        if fraction == 0.0 { continue; }
        if !pair.contains(component.species()) {
            return Err(GasFilmError::Invalid("gas contains a foreign active binary-film species"));
        }
        has_carrier |= component.species() == carrier;
        if component.species() == vapor { vapor_fraction = fraction; }
    }
    if !has_carrier {
        return Err(GasFilmError::Invalid(
            "binary gas film requires a positive carrier inventory",
        ));
    }
    Ok(vapor_fraction)
}

pub(super) fn evaluate_original_film(
    binding: &LaminarPlateFilmBinding<'_>,
    vapor: &SpeciesId,
    surface_temperature: f64,
    gas: &ClosedGasParcel,
    previous: Option<&ResolvedMaterialStatePoint>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<(f64, f64, ResolvedMaterialStatePoint), GasFilmError> {
    let bulk = gas.state().state();
    let temperature = match binding.property_temperature {
        ChamberFilmTemperature::Bulk => bulk.temperature,
        ChamberFilmTemperature::MeanSurfaceAndBulk =>
            0.5 * surface_temperature + 0.5 * bulk.temperature,
    };
    let conditions = LaminarPlateFilmConditions {
        start: Length::new(binding.geometry.start_m),
        end: Length::new(binding.geometry.end_m),
        speed_m_s: binding.geometry.speed_m_s,
        temperature: Temperature::new(temperature),
        total_pressure: Pressure::new(bulk.pressure),
    };
    // Re-evaluate the actual CURRENT mixture at film temperature and BULK total
    // pressure. Heating a fixed-mass gas parcel would incorrectly change that
    // pressure; no fictitious parcel or energy exchange is constructed here.
    let film_state = if temperature == bulk.temperature {
        // The current parcel already owns this exact source evaluation. Reuse
        // it for the bulk convention and the equal-temperature mean limit.
        None
    } else {
        let point = gas.point().clone()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), temperature)
            .and_then(|point| point.with_quantity(
                "pressure", typed(QuantityKind::Pressure), bulk.pressure,
            ))
            .map_err(GasFilmError::Query)?;
        poll(cancelled)?;
        Some(gas.model().evaluate(&point).map_err(GasFilmError::Material)?)
    };
    poll(cancelled)?;
    let fluid = film_state.as_ref().unwrap_or_else(|| gas.state());
    let mut checkpoint = || cancelled();
    let result = super::evaluate(
        fluid, &binding.source, vapor, conditions, previous, &mut checkpoint,
    )?;
    Ok((temperature, result.mass_transfer_coefficient_m_s(), result.diffusion.material))
}

#[cfg(test)]
mod tests;
