//! Bounded textual gas admission shared by file-driven physics consumers.
//! This is a source adapter, not a constitutive law or an instrument schema.

use super::{GasSpecies, resolve_gas_mixture};
use crate::MaterialError;
use crate::gas::GasState;
use fs_qty::semantic::{Composition, CompositionBasis};
use std::str::SplitAsciiWhitespace;

/// Syntax errors are separate from the original material owner's refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum GasInputError {
    /// Missing, malformed, unsupported or over-budget source fields.
    Input(&'static str),
    /// Species, composition or physical state refused by its material owner.
    Material(MaterialError),
}
impl core::fmt::Display for GasInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input(what) => write!(f, "gas input: {what}"),
            Self::Material(source) => write!(f, "gas material: {source}"),
        }
    }
}
impl std::error::Error for GasInputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Material(source) => Some(source), Self::Input(_) => None }
    }
}

/// Consume exactly one gas declaration, leaving any caller-owned suffix intact.
///
/// Supported records (SI kelvin, pascals, dimensionless fractions):
/// - `T P RH`, or `T RH` when `inherited_pressure_pa` is supplied;
/// - `mixture T P (mole|mass) COUNT SPECIES FRACTION ...` in either position.
///
/// The legacy form still uses `GasState::try_new_moist_air` without conversion.
/// The mixture form uses the shared composition and named-material resolvers.
/// Counts are checked before allocation. On error the iterator is unchanged.
/// Callers must reject unconsumed fields or parse their own explicit suffix.
///
/// Pressure inheritance affects only the legacy form. The consumer owns any
/// additional common-pressure requirement across connected regions; material
/// parsing neither assumes connectivity nor silently repairs pressure jumps.
/// This adapter does not infer eligibility for air-only absorption laws.
///
/// # Errors
/// Malformed or over-budget source, unsupported species/composition, or any
/// domain/phase refusal of the existing dry-air, moist-air or mixture owner.
pub fn read_gas_state(
    fields: &mut SplitAsciiWhitespace<'_>, inherited_pressure_pa: Option<f64>,
) -> Result<GasState, GasInputError> {
    let mut remaining = fields.clone();
    let state = read(&mut remaining, inherited_pressure_pa)?;
    *fields = remaining;
    Ok(state)
}

fn word<'a>(fields: &mut SplitAsciiWhitespace<'a>) -> Result<&'a str, GasInputError> {
    fields.next().ok_or(GasInputError::Input("missing gas field"))
}
fn scalar(fields: &mut SplitAsciiWhitespace<'_>) -> Result<f64, GasInputError> {
    let value: f64 = word(fields)?.parse()
        .map_err(|_| GasInputError::Input("invalid numeric gas field"))?;
    if !value.is_finite() { return Err(GasInputError::Input("gas scalars must be finite")); }
    Ok(value)
}
fn read(fields: &mut SplitAsciiWhitespace<'_>, inherited: Option<f64>) -> Result<GasState, GasInputError> {
    if fields.clone().next() != Some("mixture") {
        let temperature = scalar(fields)?;
        let pressure = match inherited { Some(p) => p, None => scalar(fields)? };
        let humidity = scalar(fields)?;
        return GasState::try_new_moist_air(temperature, pressure, humidity)
            .map_err(GasInputError::Material);
    }
    word(fields)?;
    let temperature = scalar(fields)?;
    let pressure = scalar(fields)?;
    let basis = match word(fields)? {
        "mole" => CompositionBasis::MoleFraction,
        "mass" => CompositionBasis::MassFraction,
        _ => return Err(GasInputError::Input("gas mixture basis must be mole or mass")),
    };
    let count: usize = word(fields)?.parse()
        .map_err(|_| GasInputError::Input("invalid gas species count"))?;
    if count == 0 || count > GasSpecies::ALL.len() {
        return Err(GasInputError::Input("gas mixture count exceeds the supported species budget"));
    }
    let mut species = Vec::with_capacity(count);
    let mut fractions = Vec::with_capacity(count);
    for _ in 0..count {
        species.push(word(fields)?.parse::<GasSpecies>().map_err(GasInputError::Material)?);
        fractions.push(scalar(fields)?);
    }
    let composition = Composition::new(basis, fractions).map_err(|error|
        GasInputError::Material(MaterialError::Parameters { what: format!("gas composition: {error}") }))?;
    resolve_gas_mixture(&species, &composition, temperature, pressure).map_err(GasInputError::Material)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn complete(text: &str, inherited: Option<f64>) -> Result<GasState, GasInputError> {
        let mut fields = text.split_ascii_whitespace();
        let state = read_gas_state(&mut fields, inherited)?;
        if fields.next().is_some() { return Err(GasInputError::Input("extra gas field")); }
        Ok(state)
    }

    #[test]
    fn legacy_dry_and_humid_records_retain_the_original_owner_bits() {
        for (t, p, rh) in [(293.15, 101325.0, 0.0), (303.15, 90000.0, 0.6), (283.15, 101325.0, 1.0)] {
            let expected = GasState::try_new_moist_air(t, p, rh).unwrap();
            assert_eq!(complete(&format!("{t} {p} {rh}"), None).unwrap(), expected);
            assert_eq!(complete(&format!("{t} {rh}"), Some(p)).unwrap(), expected);
        }
    }

    #[test]
    fn both_fraction_bases_and_aliases_use_the_same_material_resolver() {
        for (name, basis) in [("mole", CompositionBasis::MoleFraction), ("mass", CompositionBasis::MassFraction)] {
            let species = [GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::Argon];
            let composition = Composition::new(basis, vec![0.7, 0.2, 0.1]).unwrap();
            let expected = resolve_gas_mixture(&species, &composition, 300.0, 100000.0).unwrap();
            assert_eq!(complete(&format!("mixture 300 100000 {name} 3 nitrogen 0.7 O2 0.2 Ar 0.1"), None).unwrap(), expected);
        }
    }

    #[test]
    fn prefix_consumption_is_exact_and_refusal_does_not_consume_source() {
        let mut fields = "mixture 300 100000 mole 1 N2 1 viscothermal 100 1000 6 8".split_ascii_whitespace();
        read_gas_state(&mut fields, None).unwrap();
        assert_eq!(fields.collect::<Vec<_>>(), ["viscothermal", "100", "1000", "6", "8"]);
        let text = "mixture 300 100000 mole 2 N2 0.5 unknown 0.5";
        let mut fields = text.split_ascii_whitespace();
        assert!(read_gas_state(&mut fields, None).is_err());
        assert_eq!(fields.collect::<Vec<_>>(), text.split_ascii_whitespace().collect::<Vec<_>>());
    }

    #[test]
    fn explicit_pressure_is_not_replaced_by_an_inherited_default() {
        let state = complete("mixture 300 50000 mole 1 N2 1", Some(100000.0)).unwrap();
        assert_eq!(state.pressure, 50000.0);
        // Connected acoustic consumers separately refuse this pressure jump.
        assert!(complete("300 0", Some(f64::NAN)).is_err());
    }

    #[test]
    fn malformed_and_unsupported_physics_never_become_air() {
        for text in ["", "300", "NaN 100000 0", "300 100000 2",
            "mixture 300 100000 volume 1 N2 1", "mixture 300 100000 mole 0",
            "mixture 300 100000 mole 18446744073709551615",
            "mixture 300 100000 mole 1 He 1", "mixture 300 100000 mole 2 N2 0.5 nitrogen 0.5",
            "mixture 300 100000 mole 2 N2 1 O2 0", "mixture 300 100000 mole 2 N2 -0.1 O2 1.1",
            "mixture 300 100000 mole 2 N2 70 O2 30", "mixture 300 100000 mole 1 N2 NaN",
            "mixture 300 100000 mole 1 N2 1 extra", "mixture 500 100000 mole 1 N2 1",
            "mixture 300 0 mole 1 N2 1", "mixture 300 100000 mole 1 H2O 1"] {
            assert!(complete(text, None).is_err(), "accepted {text}");
        }
    }
}
