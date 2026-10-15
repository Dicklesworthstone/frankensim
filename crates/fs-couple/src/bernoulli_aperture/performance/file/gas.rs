//! One source adapter for both ambient and regional gas declarations.
//! It selects material data and a composition basis, never gas coefficients.
use super::{PlateValveInputError, Row, bad};
use fs_material::gas::GasState;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_material::MaterialError;
use fs_qty::semantic::{Composition, CompositionBasis};

fn material(line: usize, source: MaterialError) -> PlateValveInputError {
    PlateValveInputError::Gas { line, source }
}

// None: ambient T P RH; Some(P): regional T RH, preserving the old dry/humid
// route. Either position also admits mixture T P (mole|mass) N SPECIES FRACTION...
// The caller owns the remainder (end of ambient, or the section's loss suffix).
pub(super) fn read(row: &mut Row<'_>, inherited_pressure: Option<f64>)
    -> Result<GasState, PlateValveInputError>
{
    if row.fields.clone().next() != Some("mixture") {
        let temperature = row.scalar()?;
        let pressure = match inherited_pressure { Some(p) => p, None => row.scalar()? };
        let humidity = row.scalar()?;
        return GasState::try_new_moist_air(temperature, pressure, humidity)
            .map_err(|_| bad(row.line, if inherited_pressure.is_some() {
                "section gas is outside the shared moist-air domain"
            } else { "ambient outside the shared moist-air domain" }));
    }
    row.word()?;
    let temperature = row.scalar()?;
    let pressure = row.scalar()?;
    let basis = match row.word()? {
        "mole" => CompositionBasis::MoleFraction,
        "mass" => CompositionBasis::MassFraction,
        _ => return Err(bad(row.line, "gas mixture basis must be mole or mass, not percentages or volume fractions")),
    };
    let count = row.count(GasSpecies::ALL.len())?;
    if count == 0 { return Err(bad(row.line, "gas mixture requires at least one explicit species")); }
    let mut species = Vec::with_capacity(count);
    let mut fractions = Vec::with_capacity(count);
    for _ in 0..count {
        species.push(row.word()?.parse::<GasSpecies>().map_err(|e| material(row.line, e))?);
        fractions.push(row.scalar()?);
    }
    let composition = Composition::new(basis, fractions).map_err(|e|
        material(row.line, MaterialError::Parameters { what: format!("gas composition: {e}") }))?;
    let gas = resolve_gas_mixture(&species, &composition, temperature, pressure)
        .map_err(|e| material(row.line, e))?;
    // The existing reciprocal network is linearized around a stationary
    // common-pressure background. A static jump would require a mean-flow model.
    // Check against ambient even when EVERY physical section supplies an override.
    if inherited_pressure.is_some_and(|p| p.to_bits() != gas.pressure.to_bits()) {
        return Err(bad(row.line, "section and ambient must share one static pressure; pressure-driven mean flow is not modeled"));
    }
    Ok(gas)
}

#[cfg(test)]
mod tests;
