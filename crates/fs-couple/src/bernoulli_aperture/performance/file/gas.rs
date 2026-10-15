//! Gas records use the material-owned parser; connectivity remains network-owned.
use super::{PlateValveInputError, Row, bad};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasInputError, read_gas_state};
#[cfg(test)]
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
#[cfg(test)]
use fs_qty::semantic::{Composition, CompositionBasis};

pub(super) fn read(row: &mut Row<'_>, inherited_pressure: Option<f64>)
    -> Result<GasState, PlateValveInputError>
{
    let gas = read_gas_state(&mut row.fields, inherited_pressure).map_err(|error| match error {
        GasInputError::Input(what) => bad(row.line, what),
        GasInputError::Material(source) => PlateValveInputError::Gas { line: row.line, source },
    })?;
    // The reciprocal graph is linearized around a stationary common-pressure
    // background. This consumer-specific check also covers all-section overrides.
    if inherited_pressure.is_some_and(|p| p.to_bits() != gas.pressure.to_bits()) {
        return Err(bad(row.line, "section and ambient must share one static pressure; pressure-driven mean flow is not modeled"));
    }
    Ok(gas)
}

#[cfg(test)]
mod tests;
