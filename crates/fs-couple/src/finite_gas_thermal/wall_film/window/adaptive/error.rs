//! Original-coordinate error observations; dry film has no caloric coordinate.
use super::*;

#[derive(Clone, Copy)]
pub(super) struct Endpoint<'a> {
    pub wall: EquilibriumPhaseState,
    pub gas: &'a ClosedGasParcel,
    pub film: GasWallFilmInventory,
}

/// Coarse/fine differences before the original method-order correction.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WallFilmTimeError {
    pub wall_enthalpy_difference_j_kg: f64,
    pub wall_phase_fraction_difference: f64,
    pub wall_temperature_difference_k: f64,
    /// Present only when BOTH compared endpoints contain actual liquid.
    pub film_enthalpy_difference_j_kg: Option<f64>,
    pub film_mass_difference_kg: f64,
    pub gas_temperature_difference_k: f64,
    pub maximum_species_mass_difference_kg: f64,
    pub normalized: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn compare(
    environment: &WallFilmThermalEnvironment<'_>,
    cx: &Cx<'_>,
    initial: Endpoint<'_>,
    coarse: Endpoint<'_>,
    fine: Endpoint<'_>,
    controls: WallFilmAdaptiveControls,
    correction: f64,
) -> Result<WallFilmTimeError, WallFilmAdaptiveError> {
    poll(cx)?;
    let mut result = WallFilmTimeError::default();
    let mut component = |old: f64, a: f64, b: f64, absolute: f64| {
        let difference = finite((a - b).abs())?;
        let change = finite((a - old).abs().max((b - old).abs()))?;
        let allowance = finite(absolute + controls.relative_tolerance * change)?;
        let normalized = (difference / allowance) * correction;
        if normalized.is_nan() {
            return Err(WallFilmAdaptiveError::Invalid(
                "unrepresentable wall-film time error",
            ));
        }
        result.normalized = result.normalized.max(normalized);
        Ok::<_, WallFilmAdaptiveError>(difference)
    };
    result.wall_enthalpy_difference_j_kg = component(
        initial.wall.specific_enthalpy_j_kg(),
        coarse.wall.specific_enthalpy_j_kg(),
        fine.wall.specific_enthalpy_j_kg(),
        controls.enthalpy_absolute_tolerance_j_kg,
    )?;
    result.wall_phase_fraction_difference = component(
        initial.wall.liquid_mass_fraction(),
        coarse.wall.liquid_mass_fraction(),
        fine.wall.liquid_mass_fraction(),
        controls.phase_fraction_absolute_tolerance,
    )?;
    result.wall_temperature_difference_k = component(
        initial.wall.temperature_k(),
        coarse.wall.temperature_k(),
        fine.wall.temperature_k(),
        controls.temperature_absolute_tolerance_k,
    )?;
    result.gas_temperature_difference_k = component(
        initial.gas.state().state().temperature,
        coarse.gas.state().state().temperature,
        fine.gas.state().state().temperature,
        controls.temperature_absolute_tolerance_k,
    )?;
    result.film_mass_difference_kg = component(
        initial.film.mass_kg(),
        coarse.film.mass_kg(),
        fine.film.mass_kg(),
        controls.mass_absolute_tolerance_kg,
    )?;
    if let (Some(a), Some(b)) = (
        coarse.film.specific_enthalpy_j_kg(),
        fine.film.specific_enthalpy_j_kg(),
    ) {
        let old = match initial.film.specific_enthalpy_j_kg() {
            Some(value) => value,
            None => {
                (environment.boundary.liquid_enthalpy_at_temperature)(initial.wall.temperature_k())
                    .map_err(|error| {
                        WallFilmWindowError::Step(WallFilmThermalError::Phase(error.into()))
                    })?
            }
        };
        result.film_enthalpy_difference_j_kg = Some(component(
            old,
            a,
            b,
            controls.enthalpy_absolute_tolerance_j_kg,
        )?);
    }
    let original = initial
        .gas
        .species_inventory_kg()
        .map_err(WallFilmWindowError::from)?;
    let coarse_mass = coarse
        .gas
        .species_inventory_kg()
        .map_err(WallFilmWindowError::from)?;
    let fine_mass = fine
        .gas
        .species_inventory_kg()
        .map_err(WallFilmWindowError::from)?;
    for source in environment.catalog.thermodynamics().components() {
        poll(cx)?;
        let id = source.species();
        let difference = component(
            original.get(id).copied().unwrap_or(0.0),
            coarse_mass.get(id).copied().unwrap_or(0.0),
            fine_mass.get(id).copied().unwrap_or(0.0),
            controls.mass_absolute_tolerance_kg,
        )?;
        result.maximum_species_mass_difference_kg =
            result.maximum_species_mass_difference_kg.max(difference);
    }
    poll(cx)?;
    Ok(result)
}
