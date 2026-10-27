//! Shared explicit gas-source loading for native render entry points.
//! Physics libraries still receive already-loaded sources and never open stores.

use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack,
    NormalizedSpeciesPack,
};
use fs_material::{
    gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources, read_gas_record},
    state_point::MaterialPropertySelection,
    thermochemical_gas::packs::SpeciesPackRef,
};
use fs_qty::chemistry::SpeciesId;
use std::{
    ffi::OsString,
    io::Read,
    path::{Path, PathBuf},
};

const MAX_COMPONENTS: usize = 128;
const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
const MAX_TOTAL_BYTES: usize = 64 * 1024 * 1024;

#[derive(Debug, Default)]
pub(super) struct Options {
    pub record: Option<PathBuf>,
    components: Vec<[PathBuf; 3]>,
    water: Option<Option<SpeciesId>>,
    phase: Option<ThermochemicalPhase>,
}

impl Options {
    /// Remove only explicitly named gas options, preserving legacy arguments.
    /// Each component is SPECIES.pack NASA9.pack TRANSPORT.pack. Transport may
    /// be the original material-card pack or a legacy pure-transport claim pack;
    /// the declaration explicitly associates it with this species in both cases.
    pub fn extract(args: &[String]) -> Result<(Vec<String>, Self), String> {
        let words = args.iter().map(OsString::from).collect::<Vec<_>>();
        let (remaining, options) = Self::extract_os(&words)?;
        Ok((remaining.into_iter().map(|word| word.into_string()
            .expect("remaining arguments came from UTF-8 strings")).collect(), options))
    }

    /// The same options with native paths, retaining unrelated non-UTF-8 paths.
    pub fn extract_os(args: &[OsString]) -> Result<(Vec<OsString>, Self), String> {
        let mut result = Self::default();
        let mut remaining = Vec::new();
        let mut words = args.iter();
        while let Some(word) = words.next() {
            match word.to_str().unwrap_or("") {
                "--gas-component" => {
                    if result.components.len() == MAX_COMPONENTS {
                        return Err("--gas-component exceeds 128 components".into());
                    }
                    let mut path = || {
                        words.next().filter(|s| !s.to_string_lossy().starts_with("--"))
                        .map(PathBuf::from).ok_or_else(|| "--gas-component needs species, NASA-9 and pure-transport pack paths".to_string())
                    };
                    result.components.push([path()?, path()?, path()?]);
                }
                "--gas-record" => {
                    if result.record.is_some() {
                        return Err("--gas-record may only be supplied once".into());
                    }
                    result.record = Some(PathBuf::from(
                        words
                            .next()
                            .filter(|s| !s.to_string_lossy().starts_with("--"))
                            .ok_or("--gas-record needs a gas declaration file")?,
                    ));
                }
                "--gas-water" => {
                    if result.water.is_some() {
                        return Err("--gas-water may only be supplied once".into());
                    }
                    result.water = Some(match words.next().and_then(|word| word.to_str()) {
                        Some("none") => None,
                        Some(id) => {
                            Some(SpeciesId::new(id).map_err(|e| format!("--gas-water: {e}"))?)
                        }
                        None => {
                            return Err("--gas-water needs none or a canonical species id".into());
                        }
                    });
                }
                "--gas-phase" => {
                    if result.phase.is_some() {
                        return Err("--gas-phase may only be supplied once".into());
                    }
                    result.phase = Some(match words.next().and_then(|word| word.to_str()) {
                        Some("frozen") => ThermochemicalPhase::FrozenGas,
                        Some("liquid-water") => ThermochemicalPhase::LiquidWaterLimit,
                        _ => return Err("--gas-phase needs frozen or liquid-water".into()),
                    });
                }
                _ => remaining.push(word.clone()),
            }
        }
        if !result.components.is_empty() && (result.water.is_none() || result.phase.is_none()) {
            return Err("--gas-component requires explicit --gas-water and --gas-phase; no water identity or phase condition is defaulted".into());
        }
        if result.components.is_empty() && (result.water.is_some() || result.phase.is_some()) {
            return Err(
                "gas water/phase declarations require explicit component source packs".into(),
            );
        }
        if result.phase == Some(ThermochemicalPhase::LiquidWaterLimit)
            && result.water.as_ref().is_some_and(Option::is_none)
        {
            return Err("liquid-water phase admission needs an explicit water species".into());
        }
        Ok((remaining, result))
    }

    /// Whether these explicit options would resolve source data outside the input.
    pub fn uses_external_files(&self) -> bool {
        self.record.is_some() || !self.components.is_empty()
    }

    pub fn load(&self) -> Result<LoadedSources, String> {
        self.load_with(|path| bounded_read(path, MAX_PACK_BYTES))
    }

    /// Decode the same native sources through an explicitly supplied byte reader.
    /// Restricted callers can use an immutable attachment map without permitting
    /// fallback filesystem access. Every successful read still obeys the native
    /// per-file and total-byte limits; the callback cannot weaken those bounds.
    pub fn load_with(
        &self,
        mut read_source: impl FnMut(&Path) -> Result<Vec<u8>, String>,
    ) -> Result<LoadedSources, String> {
        let mut loaded = Vec::with_capacity(self.components.len());
        let mut total_bytes = 0usize;
        for paths in &self.components {
            let mut read = |path: &Path| -> Result<Vec<u8>, String> {
                let bytes = read_source(path)?;
                if bytes.is_empty() || bytes.len() > MAX_PACK_BYTES {
                    return Err(format!("gas source {} must contain 1 byte..16 MiB", path.display()));
                }
                total_bytes = total_bytes.checked_add(bytes.len())
                    .ok_or_else(|| "gas source byte count overflow".to_string())?;
                if total_bytes > MAX_TOTAL_BYTES {
                    return Err("gas source packs exceed the aggregate 64 MiB budget".into());
                }
                Ok(bytes)
            };
            let species = NormalizedSpeciesPack::from_bytes(&read(&paths[0])?)
                .map_err(|e| format!("species pack: {e}"))?;
            let nasa9 = NormalizedModelPack::from_bytes(&read(&paths[1])?)
                .map_err(|e| format!("NASA-9 pack: {e}"))?;
            let transport_bytes = read(&paths[2])?;
            let transport = decode_transport(&transport_bytes, &species)?;
            loaded.push(Component { species, nasa9, transport });
        }
        Ok(LoadedSources {
            components: loaded,
            water: self.water.clone().flatten(),
            phase: self.phase,
        })
    }
}

fn decode_transport(
    bytes: &[u8],
    species: &NormalizedSpeciesPack,
) -> Result<NormalizedMaterialCardPack, String> {
    // Select only by the canonical wire header, never filename or chemistry.
    // An invalid card must not be retried as legacy data or replaced by a preset.
    if bytes.starts_with(b"FSMCDPK\0") {
        return NormalizedMaterialCardPack::from_bytes(bytes)
            .map_err(|e| format!("pure transport material-card pack: {e}"));
    }
    let claims = NormalizedPack::from_bytes(bytes)
        .map_err(|e| format!("pure transport claim pack: {e}"))?;
    // Preserve the historical CLI association and its identity for raw claim
    // packs. Canonical material cards above retain their own exact identity,
    // process, model association and provenance, matching the live gas owner.
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: species.association().species().as_str().into(),
            phase: species.association().standard_state_phase().into(),
            process: "caller-declared gas-component source association".into(),
            revision: 0,
        },
        claims,
    ).map_err(|e| format!("transport association: {e}"))
}

struct Component {
    species: NormalizedSpeciesPack,
    nasa9: NormalizedModelPack,
    transport: NormalizedMaterialCardPack,
}
pub(super) struct LoadedSources {
    components: Vec<Component>,
    water: Option<SpeciesId>,
    phase: Option<ThermochemicalPhase>,
}
impl LoadedSources {
    pub fn bindings(&self) -> Vec<SpeciesPackRef<'_>> {
        self.components
            .iter()
            .map(|source| SpeciesPackRef {
                species: &source.species,
                nasa9: Some(&source.nasa9),
                transport: Some((
                    &source.transport,
                    MaterialPropertySelection::SingleClaimOnly,
                )),
            })
            .collect()
    }
    pub fn registry<'a>(
        &'a self,
        bindings: &'a [SpeciesPackRef<'a>],
    ) -> Option<ThermochemicalSources<'a>> {
        self.phase.map(|phase| ThermochemicalSources {
            species: bindings,
            water_species: self.water.as_ref(),
            phase,
        })
    }
}

fn bounded_read(path: &Path, limit: usize) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| format!("gas source {}: {e}", path.display()))?
        .take((limit + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("gas source {}: {e}", path.display()))?;
    if bytes.len() > limit {
        return Err(format!(
            "gas source {} exceeds {limit} bytes",
            path.display()
        ));
    }
    Ok(bytes)
}

pub(super) fn read_record(
    path: &Path,
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<GasRecord, String> {
    let bytes = bounded_read(path, 1 << 20)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| "gas record must be UTF-8".to_string())?;
    let mut words = text.split_ascii_whitespace();
    let record = read_gas_record(&mut words, None, sources).map_err(|e| e.to_string())?;
    if words.next().is_some() {
        return Err("gas record contains trailing fields".into());
    }
    Ok(record)
}

/// Empty for legacy records, preserving the original deterministic sidecars.
pub(super) fn provenance(records: &[GasRecord]) -> String {
    let mut rows = Vec::new();
    for record in records {
        let GasRecord::Thermochemical(resolved) = record else {
            continue;
        };
        let gas = resolved.state();
        let receipt = resolved.thermodynamics().receipt();
        let components = receipt
            .species()
            .iter()
            .zip(receipt.mole_fraction_bits())
            .map(|(id, bits)| {
                format!(
                    "{{\"species\":{},\"mole_fraction\":{:e}}}",
                    super::json_string(id.as_str()),
                    f64::from_bits(*bits)
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let packs = resolved.source_packs().iter().map(|pack| format!(
            "{{\"species\":{},\"species_pack\":\"{}\",\"nasa9_pack\":\"{}\",\"transport_pack\":\"{}\"}}",
            super::json_string(pack.species.as_str()),pack.species_pack,pack.nasa9_pack,pack.transport_pack)).collect::<Vec<_>>().join(",");
        rows.push(format!("{{\"law\":\"NASA-9 frozen ideal gas; Wilke/WMS transport\",\"temperature_K\":{:e},\"pressure_Pa\":{:e},\"composition\":[{components}],\"packs\":[{packs}]}}",gas.temperature,gas.pressure));
    }
    if rows.is_empty() {
        String::new()
    } else {
        format!(",\"thermochemical_gases\":[{}]", rows.join(","))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn legacy_options_are_byte_unchanged_and_source_policy_is_explicit() {
        let args = vec![
            "reed".into(),
            "out.wav".into(),
            "--block".into(),
            "37".into(),
        ];
        let (remaining, options) = Options::extract(&args).unwrap();
        assert_eq!(remaining, args);
        let loaded = options.load().unwrap();
        assert!(loaded.registry(&[]).is_none());
        assert!(provenance(&[]).is_empty());
        for args in [
            vec!["--gas-component", "species", "nasa9", "transport"],
            vec!["--gas-water", "none"],
            vec!["--gas-phase", "frozen"],
            vec![
                "--gas-component",
                "species",
                "nasa9",
                "transport",
                "--gas-water",
                "none",
                "--gas-phase",
                "liquid-water",
            ],
        ] {
            assert!(
                Options::extract(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err()
            );
        }
    }
    #[test]
    fn component_cap_is_checked_before_source_io() {
        let args = |count| {
            let mut words: Vec<String> = ["--gas-water", "none", "--gas-phase", "frozen"]
                .into_iter()
                .map(String::from)
                .collect();
            for _ in 0..count {
                words.extend(
                    ["--gas-component", "species", "nasa9", "transport"]
                        .into_iter()
                        .map(String::from),
                );
            }
            words
        };
        assert!(Options::extract(&args(MAX_COMPONENTS)).is_ok());
        assert!(
            Options::extract(&args(MAX_COMPONENTS + 1))
                .unwrap_err()
                .contains("128")
        );
    }
}
