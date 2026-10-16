//! Bounded textual gas admission shared by file-driven physics consumers.
//! This is a source adapter, not a constitutive law or an instrument schema.

use super::{GasSpecies, resolve_gas_mixture};
use crate::MaterialError;
use crate::gas::GasState;
use fs_qty::semantic::{Composition, CompositionBasis};
use crate::thermochemical_gas::{
    ResolvedThermochemicalGasState,
    packs::{SpeciesPackRef, from_packs},
};
use fs_thermochem::{SpeciesId, mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};
use std::str::SplitAsciiWhitespace;

/// Phase admission is caller policy, separate from frozen ideal-gas chemistry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermochemicalPhase {
    /// Frozen gas; no statement about condensation or phase equilibrium.
    FrozenGas,
    /// Enforce the existing liquid-water saturation law in its stated window.
    LiquidWaterLimit,
}

/// Already-loaded, explicitly associated artifacts. This L3 adapter opens no
/// files or stores and never interprets a species name as transport data.
pub struct ThermochemicalSources<'a> {
    /// Available canonical species and their explicitly paired source packs.
    pub species: &'a [SpeciesPackRef<'a>],
    /// Caller-declared water identity; absent means no water-fraction authority.
    pub water_species: Option<&'a SpeciesId>,
    /// Explicit phase-admission condition for thermochemical records.
    pub phase: ThermochemicalPhase,
}

/// A numerical view retains its thermochemical source authority when present.
#[derive(Clone, Debug, PartialEq)]
pub enum GasRecord {
    /// The unchanged legacy (T,p,RH) owner.
    Legacy(GasState),
    /// The unchanged bounded, calorically perfect named-species owner.
    Mixture(GasState),
    /// NASA-9 evaluation, composition, transport and original pack identities.
    Thermochemical(ResolvedThermochemicalGasState),
}
impl GasRecord {
    /// Borrow the exact coefficients; no conversion or gas-law reevaluation.
    #[must_use]
    pub fn state(&self) -> &GasState {
        match self {
            Self::Legacy(gas) | Self::Mixture(gas) => gas,
            Self::Thermochemical(resolved) => resolved.state(),
        }
    }
}

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
    fields: &mut SplitAsciiWhitespace<'_>,
    inherited_pressure_pa: Option<f64>,
) -> Result<GasState, GasInputError> {
    Ok(*read_gas_record(fields, inherited_pressure_pa, None)?.state())
}

/// Read the legacy records or `thermochemical T P (mole|mass) COUNT ID FRACTION…`.
/// Thermochemical input requires an already-loaded registry; missing sources,
/// ambiguous ids and out-of-domain coefficients refuse without consuming input.
/// The returned record, rather than its GasState view, owns source receipts.
pub fn read_gas_record(
    fields: &mut SplitAsciiWhitespace<'_>,
    inherited_pressure_pa: Option<f64>,
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<GasRecord, GasInputError> {
    let mut remaining = fields.clone();
    let kind = remaining.clone().next();
    let record = if kind == Some("thermochemical") {
        let sources = sources.ok_or(GasInputError::Input(
            "thermochemical gas requires an explicit source registry",
        ))?;
        word(&mut remaining)?;
        let temperature = scalar(&mut remaining)?;
        let pressure = scalar(&mut remaining)?;
        let basis = match word(&mut remaining)? {
            "mole" => CompositionBasis::MoleFraction,
            "mass" => CompositionBasis::MassFraction,
            _ => {
                return Err(GasInputError::Input(
                    "thermochemical basis must be mole or mass",
                ));
            }
        };
        let count: usize = word(&mut remaining)?
            .parse()
            .map_err(|_| GasInputError::Input("invalid thermochemical component count"))?;
        if count == 0 || count > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
            return Err(GasInputError::Input(
                "thermochemical count exceeds MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1",
            ));
        }
        let mut selected = Vec::with_capacity(count);
        let mut fractions = Vec::with_capacity(count);
        for _ in 0..count {
            let name = word(&mut remaining)?;
            let mut candidates = sources
                .species
                .iter()
                .filter(|source| source.species.association().species().as_str() == name);
            let source = candidates.next().ok_or(GasInputError::Input(
                "thermochemical species has no source registry entry",
            ))?;
            if candidates.next().is_some() {
                return Err(GasInputError::Input(
                    "ambiguous thermochemical species registry entry",
                ));
            }
            selected.push(source.clone());
            fractions.push(scalar(&mut remaining)?);
        }
        let water = sources.water_species.filter(|id| {
            selected
                .iter()
                .any(|source| source.species.association().species() == *id)
        });
        let composition = Composition::new(basis, fractions).map_err(|error| {
            GasInputError::Material(MaterialError::Parameters {
                what: format!("thermochemical composition: {error}"),
            })
        })?;
        let model = from_packs(&selected, composition, water).map_err(GasInputError::Material)?;
        use fs_qty::{
            QuantitySpec,
            semantic::{QuantityKind, SemanticType, ValueForm},
        };
        let mut point = fs_matdb::QueryPoint::new();
        for (axis, kind, value) in [
            (
                "temperature",
                QuantityKind::AbsoluteTemperature,
                temperature,
            ),
            ("pressure", QuantityKind::Pressure, pressure),
        ] {
            point = point
                .with_quantity(
                    axis,
                    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)),
                    value,
                )
                .map_err(|_| GasInputError::Input("invalid thermochemical physical coordinate"))?;
        }
        let resolved = match sources.phase {
            ThermochemicalPhase::FrozenGas => model.evaluate(&point),
            ThermochemicalPhase::LiquidWaterLimit => model.evaluate_with_liquid_water_limit(&point),
        }
        .map_err(GasInputError::Material)?;
        GasRecord::Thermochemical(resolved)
    } else {
        let gas = read(&mut remaining, inherited_pressure_pa)?;
        if kind == Some("mixture") {
            GasRecord::Mixture(gas)
        } else {
            GasRecord::Legacy(gas)
        }
    };
    *fields = remaining;
    Ok(record)
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

    #[test]
    fn thermochemical_without_registry_refuses_without_consuming() {
        let text = "thermochemical 300 101325 mole 1 N2 1";
        let mut fields = text.split_ascii_whitespace();
        assert!(
            matches!(read_gas_state(&mut fields, None), Err(GasInputError::Input(message))
            if message.contains("explicit source registry"))
        );
        assert_eq!(
            fields.collect::<Vec<_>>(),
            text.split_ascii_whitespace().collect::<Vec<_>>()
        );
    }

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
