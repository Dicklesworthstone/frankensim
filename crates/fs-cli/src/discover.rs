//! Bounded CLI adapter for compound material discovery. Selection, support,
//! and executable admission remain owned by fs-matdb-store/fs-material.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use fs_blake3::ContentHash;
use fs_matdb::{
    LawId, MaterialStateId, NormalizedInterfacePack, NormalizedMaterialCardPack,
    NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack, PropertyKey, QueryPoint,
    SelectionPolicy,
};
use fs_matdb_store::{
    CatalogPack, DiscoveryDomain, DiscoveryReport, DiscoveryRequest, DiscoveryStatus,
    DiscoveryTarget, MaterialStore, ModelRequirement,
};
use fs_material::graph::LawRegistry;
use fs_qty::parse::parse_qty;
use fs_qty::semantic::{QuantityKind, QuantitySpec, SemanticType, ValueForm};

use crate::json_read::JsonValue;
use crate::{CommandOutput, Diagnostic, OutputMode, escape_text, exit, push_json_string, refusal};

const REQUEST_BYTES: u64 = 64 * 1024;
const PACK_BYTES: u64 = 16 * 1024 * 1024;
const TOTAL_PACK_BYTES: u64 = 64 * 1024 * 1024;
const MAX_PACKS: usize = 32;
const MAX_ITEMS: usize = 64;
const NO_CLAIM: &str = "declared scalar data/model coverage only; no specimen inference, solver initialization, trajectory prediction or physical validation";

/// Discover through canonical pack files and an explicit executable registry.
/// The command-line binary supplies an empty registry; embedders may supply
/// their actual factories. The request and packs are read with fixed byte caps,
/// then ingested transactionally into a temporary in-memory store. No input
/// file or persistent catalog is modified. A successfully produced gap report
/// exits zero even when no candidate supplies the complete requested bundle.
#[must_use]
pub fn discover_paths(
    request: &Path,
    packs: &[PathBuf],
    json: bool,
    registry: &LawRegistry,
) -> CommandOutput {
    let mode = if json {
        OutputMode::Json
    } else {
        OutputMode::Text
    };
    match execute(request, packs, registry) {
        Ok((report, identities)) => render(&report, &identities, mode),
        Err((class, code, detail)) => refusal(
            mode,
            class,
            &Diagnostic::new(
                "discover",
                code,
                detail,
                "use the documented discovery request and admitted canonical pack files; correct the named input or requirement",
            ),
            None,
        ),
    }
}

type Failure = (u8, &'static str, String);

fn execute(
    request: &Path,
    paths: &[PathBuf],
    registry: &LawRegistry,
) -> Result<(DiscoveryReport, BTreeMap<String, String>), Failure> {
    if paths.is_empty() || paths.len() > MAX_PACKS {
        return Err((
            exit::INPUT,
            "cli-discover-pack-count",
            format!("expected 1..={MAX_PACKS} packs"),
        ));
    }
    let bytes = read(request, REQUEST_BYTES)?;
    let source = std::str::from_utf8(&bytes)
        .map_err(|e| (exit::INPUT, "cli-discover-encoding", e.to_string()))?;
    let request = parse_request(source).map_err(|e| (exit::REFUSED, "cli-discover-request", e))?;
    let mut packs = Vec::new();
    let mut total = 0_u64;
    for path in paths {
        let bytes = read(path, PACK_BYTES.min(TOTAL_PACK_BYTES - total))?;
        total += bytes.len() as u64;
        let pack = match bytes.get(..8) {
            Some(b"FSMATPK\0") => NormalizedPack::from_bytes(&bytes).map(CatalogPack::Properties),
            Some(b"FSMCDPK\0") => {
                NormalizedMaterialCardPack::from_bytes(&bytes).map(CatalogPack::MaterialCard)
            }
            Some(b"FSINTPK\0") => {
                NormalizedInterfacePack::from_bytes(&bytes).map(CatalogPack::Interface)
            }
            Some(b"FSMODPK\0") => NormalizedModelPack::from_bytes(&bytes).map(CatalogPack::Model),
            Some(b"FSSPCPK\0") => {
                NormalizedSpeciesPack::from_bytes(&bytes).map(CatalogPack::Species)
            }
            _ => {
                return Err((
                    exit::REFUSED,
                    "cli-discover-pack",
                    format!("{}: unknown canonical pack family", path.display()),
                ));
            }
        }
        .map_err(|e| {
            (
                exit::REFUSED,
                "cli-discover-pack",
                format!("{}: {e}", path.display()),
            )
        })?;
        // Exact duplicates are the same input; conflicting pack names remain
        // the store's transactional refusal, never last-one-wins selection.
        if !packs
            .iter()
            .any(|p: &CatalogPack| p.content_hash() == pack.content_hash())
        {
            packs.push(pack);
        }
    }
    let identities = packs
        .iter()
        .map(|pack| {
            let identity = match pack {
                CatalogPack::MaterialCard(p) => p.card().id().to_string(),
                CatalogPack::Interface(p) => format!(
                    "{} -> {}",
                    p.card().surface_a().material,
                    p.card().surface_b().material
                ),
                _ => "unbound to a named material state".to_owned(),
            };
            (pack.pack_id().to_owned(), identity)
        })
        .collect();
    let store = MaterialStore::open(":memory:")
        .map_err(|e| (exit::REFUSED, "cli-discover-store", e.to_string()))?;
    store
        .ingest_bundle(&packs)
        .map_err(|e| (exit::REFUSED, "cli-discover-store", e.to_string()))?;
    store
        .seal_corpus()
        .map_err(|e| (exit::REFUSED, "cli-discover-store", e.to_string()))?;
    store
        .discover(&request, registry)
        .map(|report| (report, identities))
        .map_err(|e| (exit::REFUSED, "cli-discover-request", e.to_string()))
}

fn read(path: &Path, cap: u64) -> Result<Vec<u8>, Failure> {
    let file = std::fs::File::open(path).map_err(|e| {
        (
            exit::INPUT,
            "cli-discover-read",
            format!("{}: {e}", path.display()),
        )
    })?;
    let metadata = file
        .metadata()
        .map_err(|e| (exit::INPUT, "cli-discover-read", e.to_string()))?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err((
            exit::INPUT,
            "cli-discover-size",
            format!(
                "{}: expected a regular file of at most {cap} bytes",
                path.display()
            ),
        ));
    }
    let mut bytes = Vec::new();
    file.take(cap + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| (exit::INPUT, "cli-discover-read", e.to_string()))?;
    if bytes.len() as u64 > cap {
        return Err((
            exit::INPUT,
            "cli-discover-size",
            format!("{} grew beyond {cap} bytes", path.display()),
        ));
    }
    Ok(bytes)
}

fn fields(value: &JsonValue, allowed: &[&str]) -> Result<(), String> {
    let object = value.as_object().ok_or("expected an object")?;
    if let Some((name, _)) = object
        .iter()
        .find(|(name, _)| !allowed.contains(&name.as_str()))
    {
        return Err(format!(
            "unsupported field {name:?}; supported fields: {allowed:?}"
        ));
    }
    Ok(())
}

fn string<'a>(value: &'a JsonValue, name: &str) -> Result<&'a str, String> {
    value
        .str_field(name)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| format!("{name}: expected a nonempty string"))
}

fn items<'a>(value: &'a JsonValue, name: &str) -> Result<&'a [JsonValue], String> {
    let values = value
        .get(name)
        .and_then(JsonValue::as_array)
        .ok_or_else(|| format!("{name}: expected an array"))?;
    if values.len() > MAX_ITEMS {
        return Err(format!("{name}: exceeds {MAX_ITEMS} entries"));
    }
    Ok(values)
}

fn integer(value: &JsonValue, name: &str) -> Result<u32, String> {
    value
        .get(name)
        .and_then(JsonValue::number_raw)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| format!("{name}: expected an unsigned 32-bit integer"))
}

fn state(value: &JsonValue) -> Result<MaterialStateId, String> {
    fields(value, &["chemistry", "phase", "process", "revision"])?;
    Ok(MaterialStateId {
        chemistry: string(value, "chemistry")?.into(),
        phase: string(value, "phase")?.into(),
        process: string(value, "process")?.into(),
        revision: integer(value, "revision")?,
    })
}

fn quantity(kind: &str, dims: fs_qty::Dims) -> Result<QuantitySpec, String> {
    if kind == "dimensional" {
        return Ok(QuantitySpec::dimensional(dims));
    }
    let kind = match kind {
        "absolute-temperature" => QuantityKind::AbsoluteTemperature,
        "temperature-difference" => QuantityKind::TemperatureDifference,
        "energy" => QuantityKind::Energy,
        "torque" => QuantityKind::Torque,
        "pressure" => QuantityKind::Pressure,
        "stress" => QuantityKind::Stress,
        "heat-capacity" => QuantityKind::HeatCapacity,
        "entropy" => QuantityKind::Entropy,
        name => QuantityKind::from_material_convention_name(name)
            .ok_or_else(|| format!("unsupported scalar quantity kind {name:?}"))?,
    };
    if kind.expected_dims() != dims {
        return Err(format!(
            "quantity kind {kind:?} disagrees with unit dimensions {dims:?}"
        ));
    }
    Ok(QuantitySpec::semantic(SemanticType::new(
        kind,
        ValueForm::Static,
    )))
}

fn parse_request(source: &str) -> Result<DiscoveryRequest, String> {
    let root = JsonValue::parse(source).map_err(|e| e.to_string())?;
    fields(
        &root,
        &[
            "schema",
            "target",
            "properties",
            "models",
            "domain",
            "selection",
        ],
    )?;
    if string(&root, "schema")? != "frankensim.discovery.v1" {
        return Err("unsupported discovery schema".into());
    }
    let target = root.get("target").ok_or("target is required")?;
    let target = match target.as_str() {
        Some("materials") => DiscoveryTarget::Materials,
        Some("properties") => DiscoveryTarget::Properties,
        Some(_) => {
            return Err("target must be materials, properties or an ordered surface pair".into());
        }
        None => {
            fields(target, &["surface_a", "surface_b"])?;
            DiscoveryTarget::Interfaces {
                surface_a: state(target.get("surface_a").ok_or("surface_a required")?)?,
                surface_b: state(target.get("surface_b").ok_or("surface_b required")?)?,
            }
        }
    };
    let mut properties = Vec::new();
    for item in items(&root, "properties")? {
        fields(item, &["name", "unit", "kind"])?;
        let unit = string(item, "unit")?;
        let literal = if unit == "1" {
            "1".to_owned()
        } else {
            format!("1 {unit}")
        };
        let unit = parse_qty(&literal).map_err(|e| e.to_string())?;
        properties.push(PropertyKey::with_quantity(
            string(item, "name")?,
            quantity(string(item, "kind")?, unit.dims)?,
        ));
    }
    let mut models = Vec::new();
    for item in items(&root, "models")? {
        fields(item, &["law", "version", "pin"])?;
        let pin = match item.get("pin") {
            None => None,
            Some(_) => Some(
                ContentHash::from_hex(string(item, "pin")?)
                    .ok_or("pin: expected a 64-digit content hash")?,
            ),
        };
        models.push(ModelRequirement {
            law: LawId(string(item, "law")?.into()),
            law_version: integer(item, "version")?,
            pin,
        });
    }
    let domain = root.get("domain").ok_or("domain is required")?;
    fields(domain, &["mode", "axes"])?;
    let envelope = match string(domain, "mode")? {
        "local-state" => false,
        "envelope" => true,
        _ => return Err("domain.mode must be local-state or envelope".into()),
    };
    let mut lower = QueryPoint::new();
    let mut upper = QueryPoint::new();
    for axis in items(domain, "axes")? {
        let keys: &[&str] = if envelope {
            &["name", "kind", "lower", "upper"]
        } else {
            &["name", "kind", "value"]
        };
        fields(axis, keys)?;
        let name = string(axis, "name")?;
        if lower.axes().contains_key(name) {
            return Err(format!("duplicate domain axis {name:?}"));
        }
        let lo = parse_qty(string(axis, if envelope { "lower" } else { "value" })?)
            .map_err(|e| e.to_string())?;
        let hi = if envelope {
            parse_qty(string(axis, "upper")?).map_err(|e| e.to_string())?
        } else {
            lo
        };
        if lo.dims != hi.dims || lo.value > hi.value {
            return Err(format!(
                "{name}: require equal dimensions and lower <= upper"
            ));
        }
        let kind = string(axis, "kind")?;
        if kind == "legacy" {
            lower = lower.with(name, lo.value).map_err(|e| e.to_string())?;
            upper = upper.with(name, hi.value).map_err(|e| e.to_string())?;
        } else {
            let spec = quantity(kind, lo.dims)?;
            lower = lower
                .with_quantity(name, spec, lo.value)
                .map_err(|e| e.to_string())?;
            upper = upper
                .with_quantity(name, spec, hi.value)
                .map_err(|e| e.to_string())?;
        }
    }
    let selection = match string(&root, "selection")? {
        "single-claim-only" => SelectionPolicy::SingleClaimOnly,
        "prefer-observation-backed" => SelectionPolicy::PreferObservationBacked,
        _ => return Err("unsupported selection policy".into()),
    };
    Ok(DiscoveryRequest {
        target,
        properties,
        models,
        domain: if envelope {
            DiscoveryDomain::Envelope { lower, upper }
        } else {
            DiscoveryDomain::LocalState(lower)
        },
        selection,
    })
}

fn status(value: DiscoveryStatus) -> &'static str {
    match value {
        DiscoveryStatus::Complete => "complete",
        DiscoveryStatus::Partial => "partial",
        DiscoveryStatus::Unavailable => "unavailable",
    }
}

fn strings(out: &mut String, values: impl IntoIterator<Item = String>) {
    out.push('[');
    for (i, value) in values.into_iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        push_json_string(out, &value);
    }
    out.push(']');
}

fn render(
    report: &DiscoveryReport,
    identities: &BTreeMap<String, String>,
    mode: OutputMode,
) -> CommandOutput {
    let domain = match report.request.domain {
        DiscoveryDomain::LocalState(_) => "local-state",
        DiscoveryDomain::Envelope { .. } => "envelope",
    };
    let mut out = format!(
        "{{\"schema\":\"frankensim.discovery-result.v1\",\"command\":\"discover\",\"status\":\"ok\",\"domain\":\"{domain}\",\"unknown_properties\":"
    );
    strings(
        &mut out,
        report
            .unknown_properties
            .iter()
            .map(|p| p.name().to_owned()),
    );
    out.push_str(",\"unknown_models\":");
    strings(
        &mut out,
        report
            .unknown_models
            .iter()
            .map(|m| format!("{}@{}", m.law.0, m.law_version)),
    );
    out.push_str(",\"requested_state\":");
    push_json_string(&mut out, &format!("{:?}", report.request.domain));
    out.push_str(",\"candidates\":[");
    let mut text = format!(
        "command=discover domain={domain} candidates={}\n",
        report.candidates.len()
    );
    for (i, candidate) in report.candidates.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        out.push_str("{\"pack\":");
        push_json_string(&mut out, &candidate.pack.pack_id);
        let identity = &identities[&candidate.pack.pack_id];
        out.push_str(",\"identity\":");
        push_json_string(&mut out, identity);
        let _ = write!(
            out,
            ",\"hash\":\"{}\",\"kind\":\"{}\",\"status\":\"{}\",\"properties\":[",
            candidate.pack.content_hash,
            candidate.pack.kind.as_str(),
            status(candidate.status)
        );
        let _ = writeln!(
            text,
            "pack={} kind={} status={}",
            escape_text(&candidate.pack.pack_id),
            candidate.pack.kind.as_str(),
            status(candidate.status)
        );
        let _ = writeln!(text, "  identity={}", escape_text(identity));
        for (j, property) in candidate.properties.iter().enumerate() {
            if j != 0 {
                out.push(',');
            }
            out.push_str("{\"name\":");
            push_json_string(&mut out, property.property.name());
            out.push_str(",\"quantity\":");
            push_json_string(&mut out, &format!("{:?}", property.property.quantity()));
            let detail = match &property.support {
                Ok(support) => {
                    let _ = write!(
                        out,
                        ",\"status\":\"supported\",\"lower_si\":{},\"upper_si\":{},\"claim\":\"{}\"",
                        support.lower.evidence.value.value,
                        support.upper.as_ref().map_or("null".into(), |s| s
                            .evidence
                            .value
                            .value
                            .to_string()),
                        support.lower.receipt.selected.0.to_hex()
                    );
                    format!("supported lower_si={}", support.lower.evidence.value.value)
                }
                Err(gap) => {
                    out.push_str(",\"status\":\"gap\",\"detail\":");
                    let detail = format!("{gap:?}");
                    push_json_string(&mut out, &detail);
                    detail
                }
            };
            out.push('}');
            let _ = writeln!(
                text,
                "  property={} {}",
                escape_text(property.property.name()),
                escape_text(&detail)
            );
        }
        out.push_str("],\"models\":[");
        for (j, model) in candidate.models.iter().enumerate() {
            if j != 0 {
                out.push(',');
            }
            out.push_str("{\"law\":");
            push_json_string(&mut out, &model.requirement.law.0);
            let _ = write!(
                out,
                ",\"version\":{},\"status\":\"{}\",\"detail\":",
                model.requirement.law_version,
                if model.support.is_ok() {
                    "supported"
                } else {
                    "gap"
                }
            );
            let detail = match &model.support {
                Ok(card) => format!(
                    "model={} initial_state={:?}",
                    card.content_hash(),
                    card.initial_state
                ),
                Err(gap) => format!("{gap:?}"),
            };
            push_json_string(&mut out, &detail);
            out.push('}');
            let _ = writeln!(
                text,
                "  model={}@{} {}",
                escape_text(&model.requirement.law.0),
                model.requirement.law_version,
                escape_text(&detail)
            );
        }
        out.push_str("]}");
    }
    out.push_str("],\"no_claim\":");
    push_json_string(&mut out, NO_CLAIM);
    out.push_str("}\n");
    for property in &report.unknown_properties {
        let _ = writeln!(text, "unknown_property={}", escape_text(property.name()));
    }
    for model in &report.unknown_models {
        let _ = writeln!(
            text,
            "unknown_model={}@{}",
            escape_text(&model.law.0),
            model.law_version
        );
    }
    let _ = writeln!(text, "no_claim={NO_CLAIM}");
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout: if mode == OutputMode::Json { out } else { text },
        stderr: String::new(),
    }
}
