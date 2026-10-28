//! Bounded textual gas admission shared by file-driven physics consumers.
//! This is a source adapter, not a constitutive law or an instrument schema.

/// Canonical source bytes for browser and other file-independent callers.
pub mod packs;

use super::{GasSpecies, resolve_gas_mixture};
use crate::MaterialError;
use crate::gas::GasState;
use fs_qty::semantic::{Composition, CompositionBasis};
use crate::thermochemical_gas::{
    ResolvedThermochemicalGasState, ThermochemicalGasModel,
    energy::{GasEnergyError, IsochoricGasParcel},
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
    /// A record excluding the explicitly identified water species is known dry;
    /// its own thermodynamic and transport source domains still apply.
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
    /// The unchanged bounded, calorically perfect named-species owner, with
    /// the record's species and their mole fractions in record order (a mass
    /// record is converted by the species' own molar masses). A consumer that
    /// needs the composition, such as an observer law, reads it here; it is
    /// never inferred from the state's gamma or molar mass.
    Mixture(GasState, Vec<(GasSpecies, f64)>),
    /// NASA-9 evaluation, composition, transport and original pack identities.
    Thermochemical(ResolvedThermochemicalGasState),
}
impl GasRecord {
    /// Borrow the exact coefficients; no conversion or gas-law reevaluation.
    #[must_use]
    pub fn state(&self) -> &GasState {
        match self {
            Self::Legacy(gas) | Self::Mixture(gas, _) => gas,
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

/// Original source parsing or conservative parcel initialization refusal.
/// This keeps the snapshot reader's error type unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum ParcelInputError {
    /// The shared grammar, registry, source law or phase policy refused.
    Input(GasInputError),
    /// The original finite-volume gas owner refused.
    Energy(GasEnergyError),
}
impl core::fmt::Display for ParcelInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input(error) => write!(f, "finite gas input: {error}"),
            Self::Energy(error) => write!(f, "finite gas input: {error}"),
        }
    }
}
impl std::error::Error for ParcelInputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input(error) => Some(error),
            Self::Energy(error) => Some(error),
        }
    }
}
impl From<GasInputError> for ParcelInputError {
    fn from(error: GasInputError) -> Self { Self::Input(error) }
}
impl From<GasEnergyError> for ParcelInputError {
    fn from(error: GasEnergyError) -> Self { Self::Energy(error) }
}

/// Admit one finite, closed gas reservoir from the same thermochemical record.
///
/// The original model, composition, typed conditions and source packs are
/// retained by `IsochoricGasParcel`; later heat steps reevaluate those source
/// laws instead of treating this initial snapshot as a constant-capacity gas.
/// Volume is explicit [m3], and the original owner derives mass as rho*V.
///
/// Only `FrozenGas` is admitted: an initial liquid-water saturation check
/// cannot guarantee the phase of later heated/cooled states. This path does
/// not silently weaken a requested liquid-water policy or claim condensation
/// tracking. The caller must retain a valid frozen-gas domain during evolution.
///
/// The iterator is unchanged on every refusal. On success any caller-owned
/// suffix remains; complete-record consumers must reject unexpected fields.
pub fn read_isochoric_parcel(
    fields: &mut SplitAsciiWhitespace<'_>,
    sources: &ThermochemicalSources<'_>,
    volume_m3: f64,
) -> Result<IsochoricGasParcel, ParcelInputError> {
    let mut remaining = fields.clone();
    let (model, point) = read_frozen_thermochemical_model(&mut remaining, sources)?;
    let parcel = IsochoricGasParcel::new(model, point, volume_m3)?;
    *fields = remaining;
    Ok(parcel)
}

/// Retain the original source model for spatial or finite-volume evolution.
///
/// This is the same registry, composition and typed-condition parser used by
/// snapshots and closed parcels. The caller owns subsequent state evolution;
/// a single resolved snapshot must not replace its temperature-dependent laws.
/// Only an explicit frozen-gas phase policy is admitted. A one-time saturation
/// check cannot establish the phase of future cell states. On refusal the input
/// iterator is unchanged; successful callers must consume any remaining suffix.
pub fn read_frozen_thermochemical_model(
    fields: &mut SplitAsciiWhitespace<'_>,
    sources: &ThermochemicalSources<'_>,
) -> Result<(ThermochemicalGasModel, fs_matdb::QueryPoint), GasInputError> {
    if sources.phase != ThermochemicalPhase::FrozenGas {
        return Err(GasInputError::Input(
            "finite gas evolution requires an explicit FrozenGas phase policy",
        ));
    }
    let mut remaining = fields.clone();
    let result = read_thermochemical_model(&mut remaining, sources)?;
    *fields = remaining;
    Ok(result)
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
/// With LiquidWaterLimit, the registry must explicitly identify water. A dry
/// selection need not fabricate a zero-fraction component or evaluate a water
/// saturation curve; records containing that identity retain the original phase
/// check. This allows a dry exterior and humid bore to share one source policy.
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
        let (model, point) = read_thermochemical_model(&mut remaining, sources)?;
        let resolved = match sources.phase {
            ThermochemicalPhase::FrozenGas => model.evaluate(&point),
            // The source explicitly identifies water, and no selected species
            // has that canonical identity: the record is dry, not unknown.
            // Do not add a zero-fraction species or query an irrelevant water
            // saturation window. All selected source domains remain enforced.
            ThermochemicalPhase::LiquidWaterLimit
                if sources.water_species.is_some() && model.water_species().is_none() => model.evaluate(&point),
            ThermochemicalPhase::LiquidWaterLimit => model.evaluate_with_liquid_water_limit(&point),
        }
        .map_err(GasInputError::Material)?;
        GasRecord::Thermochemical(resolved)
    } else {
        match read(&mut remaining, inherited_pressure_pa)? {
            (gas, Some(composition)) => GasRecord::Mixture(gas, composition),
            (gas, None) => GasRecord::Legacy(gas),
        }
    };
    *fields = remaining;
    Ok(record)
}

/// Parse the original source model and typed conditions once for both
/// snapshot evaluation and conservative reservoir initialization.
fn read_thermochemical_model(
    fields: &mut SplitAsciiWhitespace<'_>,
    sources: &ThermochemicalSources<'_>,
) -> Result<(ThermochemicalGasModel, fs_matdb::QueryPoint), GasInputError> {
    if word(fields)? != "thermochemical" {
        return Err(GasInputError::Input(
            "finite source gas requires the thermochemical grammar",
        ));
    }
    let temperature = scalar(fields)?;
    let pressure = scalar(fields)?;
    let basis = match word(fields)? {
        "mole" => CompositionBasis::MoleFraction,
        "mass" => CompositionBasis::MassFraction,
        _ => {
            return Err(GasInputError::Input(
                "thermochemical basis must be mole or mass",
            ));
        }
    };
    let count: usize = word(fields)?
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
        let name = word(fields)?;
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
        fractions.push(scalar(fields)?);
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
    Ok((model, point))
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
/// A legacy record's state, or a mixture record's state with its species and
/// mole fractions.
type ReadGas = (GasState, Option<Vec<(GasSpecies, f64)>>);

const BREATH_WOODWIND_MANIFEST_TSV: &str = include_str!(
    "../../../../data/matdb/seed-v1/breath-woodwind-exhaled-coltman1966-nederveen1969/manifest.tsv"
);
const BREATH_WOODWIND_PROPERTIES_TSV: &str = include_str!(
    "../../../../data/matdb/seed-v1/breath-woodwind-exhaled-coltman1966-nederveen1969/properties.tsv"
);

fn resolve_seed_gas_pack(
    pack_id: &str,
    inherited_pressure_pa: Option<f64>,
) -> Result<ReadGas, GasInputError> {
    let manifest_id = BREATH_WOODWIND_MANIFEST_TSV
        .lines()
        .find_map(|line| line.strip_prefix("pack_id\t"))
        .ok_or(GasInputError::Input("missing pack_id in seed gas manifest"))?;
    if pack_id != manifest_id {
        return Err(GasInputError::Input("unknown gas pack id"));
    }
    let mut temperature = None;
    let mut pressure = None;
    let mut n2 = None;
    let mut o2 = None;
    let mut ar = None;
    let mut co2 = None;
    let mut h2o = None;
    for line in BREATH_WOODWIND_PROPERTIES_TSV.lines() {
        let mut cols = line.split('\t');
        if cols.next() != Some("scalar") {
            continue;
        }
        let _claim = cols.next();
        let _obs = cols.next();
        let Some(prop) = cols.next() else { continue };
        let Some(val_str) = cols.next() else { continue };
        let val: f64 = val_str
            .parse()
            .map_err(|_| GasInputError::Input("invalid scalar in gas pack properties"))?;
        match prop {
            "temperature" => temperature = Some(val),
            "pressure" => pressure = Some(val),
            "mole_fraction_n2" => n2 = Some(val),
            "mole_fraction_o2" => o2 = Some(val),
            "mole_fraction_ar" => ar = Some(val),
            "mole_fraction_co2" => co2 = Some(val),
            "mole_fraction_h2o" => h2o = Some(val),
            _ => {}
        }
    }
    let temperature = temperature.ok_or(GasInputError::Input("gas pack missing temperature"))?;
    let pack_pressure = pressure.ok_or(GasInputError::Input("gas pack missing pressure"))?;
    let pressure = inherited_pressure_pa.unwrap_or(pack_pressure);
    let species = [
        GasSpecies::Nitrogen,
        GasSpecies::Oxygen,
        GasSpecies::Argon,
        GasSpecies::CarbonDioxide,
        GasSpecies::WaterVapor,
    ];
    let fractions = vec![
        n2.ok_or(GasInputError::Input("gas pack missing mole_fraction_n2"))?,
        o2.ok_or(GasInputError::Input("gas pack missing mole_fraction_o2"))?,
        ar.ok_or(GasInputError::Input("gas pack missing mole_fraction_ar"))?,
        co2.ok_or(GasInputError::Input("gas pack missing mole_fraction_co2"))?,
        h2o.ok_or(GasInputError::Input("gas pack missing mole_fraction_h2o"))?,
    ];
    let composition = Composition::new(CompositionBasis::MoleFraction, fractions.clone()).map_err(
        |error| {
            GasInputError::Material(MaterialError::Parameters {
                what: format!("gas pack composition: {error}"),
            })
        },
    )?;
    let gas = resolve_gas_mixture(&species, &composition, temperature, pressure)
        .map_err(GasInputError::Material)?;
    Ok((gas, Some(species.into_iter().zip(fractions).collect())))
}

fn read(fields: &mut SplitAsciiWhitespace<'_>, inherited: Option<f64>) -> Result<ReadGas, GasInputError> {
    let first = fields.clone().next();
    if first == Some("pack") {
        word(fields)?;
        let pack_id = word(fields)?;
        return resolve_seed_gas_pack(pack_id, inherited);
    }
    if first != Some("mixture") {
        let temperature = scalar(fields)?;
        let pressure = match inherited { Some(p) => p, None => scalar(fields)? };
        let humidity = scalar(fields)?;
        return GasState::try_new_moist_air(temperature, pressure, humidity)
            .map(|gas| (gas, None))
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
    let gas = resolve_gas_mixture(&species, &composition, temperature, pressure).map_err(GasInputError::Material)?;
    let masses: Vec<_> = species.iter().map(|s| fs_qty::MolarMass::new(s.spec().molar_mass)).collect();
    let moles = composition.to_mole_fractions(&masses).map_err(|error|
        GasInputError::Material(MaterialError::Parameters { what: format!("gas composition: {error}") }))?;
    Ok((gas, Some(species.into_iter().zip(moles.fractions().iter().copied()).collect())))
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
    fn seed_breath_woodwind_pack_resolves_identically_to_explicit_mixture() {
        let mut pack_fields =
            "pack breath-woodwind-exhaled-coltman1966-nederveen1969".split_ascii_whitespace();
        let mut mix_fields = "mixture 300.15 101325 mole 5 N2 0.7534 O2 0.1591 Ar 0.0090 CO2 0.0434 H2O 0.0351"
            .split_ascii_whitespace();
        let pack_rec = read_gas_record(&mut pack_fields, None, None).unwrap();
        let mix_rec = read_gas_record(&mut mix_fields, None, None).unwrap();
        assert_eq!(pack_rec, mix_rec);
        let bad = "pack unknown-breath-pack";
        let mut bad_fields = bad.split_ascii_whitespace();
        assert!(read_gas_state(&mut bad_fields, None).is_err());
        assert_eq!(
            bad_fields.collect::<Vec<_>>(),
            bad.split_ascii_whitespace().collect::<Vec<_>>()
        );
    }

    #[test]
    fn malformed_and_unsupported_physics_never_become_air() {
        for text in ["", "300", "NaN 100000 0", "300 100000 2",
            "mixture 300 100000 volume 1 N2 1", "mixture 300 100000 mole 0",
            "mixture 300 100000 mole 18446744073709551615",
            "mixture 300 100000 mole 1 Xe 1", "mixture 300 100000 mole 2 N2 0.5 nitrogen 0.5",
            "mixture 300 100000 mole 2 N2 1 O2 0", "mixture 300 100000 mole 2 N2 -0.1 O2 1.1",
            "mixture 300 100000 mole 2 N2 70 O2 30", "mixture 300 100000 mole 1 N2 NaN",
            "mixture 300 100000 mole 1 N2 1 extra", "mixture 500 100000 mole 1 N2 1",
            "mixture 300 0 mole 1 N2 1", "mixture 300 100000 mole 1 H2O 1"] {
            assert!(complete(text, None).is_err(), "accepted {text}");
        }
    }
}
