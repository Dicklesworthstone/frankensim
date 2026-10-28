//! Search the shipped source catalog and inspect the actual compiled claims.
//! Property presence is a discovery aid, not admission of a compatible model.

use super::{
    CompiledPack, MAX_MANIFEST_BYTES, compile_manifest, parse_manifest, read_bounded_regular,
};
use fs_matdb::{ClaimSet, PropertyValue};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::path::Path;

const USAGE: &str = "usage: cargo run -p xtask -- matdb-catalog [--search TEXT]... [--property NAME]... [--pack ID]\n\
Search is case-insensitive across source/pack IDs, citations, specimens and properties.\n\
All terms and required properties must match. --pack selects an exact source\n\
directory name (internal pack IDs may repeat across sources) and prints source\n\
conditions, SI values, quantity types and observation caveats for one exact ID.\n\
Species/model inspection is metadata-only. Refused sources remain visible in\n\
unfiltered listings with their error; property filters select compiled claims.\n\
Property presence alone does not establish a common valid state or model.";

#[derive(Default)]
struct Options {
    search: Vec<String>,
    properties: BTreeSet<String>,
    pack: Option<String>,
}

fn options(args: &[String]) -> Result<Options, String> {
    let mut result = Options::default();
    let mut args = args.iter();
    while let Some(flag) = args.next() {
        if !matches!(flag.as_str(), "--search" | "--property" | "--pack") {
            return Err(format!("unknown catalog option {flag:?}\n{USAGE}"));
        }
        let value = args
            .next()
            .filter(|value| !value.trim().is_empty() && !value.starts_with("--"))
            .ok_or_else(|| format!("{flag} requires a nonempty value\n{USAGE}"))?;
        match flag.as_str() {
            "--search" => result.search.push(value.to_lowercase()),
            "--property" => {
                result.properties.insert(value.clone());
            }
            "--pack" => {
                if result.pack.replace(value.clone()).is_some() {
                    return Err("--pack may be supplied only once".into());
                }
            }
            _ => unreachable!(),
        }
    }
    Ok(result)
}

fn claims(pack: &CompiledPack) -> Option<&ClaimSet> {
    match pack {
        CompiledPack::Material(pack) => Some(pack.claims()),
        CompiledPack::MaterialCard(pack) => Some(pack.claims_pack().claims()),
        CompiledPack::Interface(pack) => Some(pack.claims_pack().claims()),
        CompiledPack::Model(_) | CompiledPack::Species(_) => None,
    }
}

fn family(pack: &CompiledPack) -> &'static str {
    match pack {
        CompiledPack::Material(_) => "material",
        CompiledPack::MaterialCard(_) => "material-card",
        CompiledPack::Interface(_) => "interface",
        CompiledPack::Model(_) => "model",
        CompiledPack::Species(_) => "species",
    }
}

/// Escape text without allowing source annotations to create extra TSV records.
fn cell(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

fn detail(output: &mut String, pack: &CompiledPack) {
    let Some(claims) = claims(pack) else {
        return;
    };
    for (id, claim) in claims.claims_ordered() {
        writeln!(
            output,
            "claim\t{}\t{}\t{:?}\t{:?}\t{:?}",
            id.0,
            cell(claim.key.name()),
            claim.key.quantity(),
            claim.interpolation,
            claim.uncertainty
        )
        .unwrap();
        match &claim.value {
            PropertyValue::Scalar { value, .. } => {
                writeln!(output, "value\t{}\t{value}", id.0).unwrap();
            }
            PropertyValue::Curve {
                abscissa, knots, ..
            } => {
                for (x, y) in knots {
                    writeln!(output, "knot\t{}\t{}\t{x}\t{y}", id.0, cell(abscissa)).unwrap();
                }
            }
            PropertyValue::Surface2D {
                axes,
                samples,
                triangles,
                ..
            } => {
                for (index, sample) in samples.iter().enumerate() {
                    writeln!(
                        output,
                        "sample\t{}\t{index}\t{}\t{}\t{}\t{}\t{}",
                        id.0,
                        cell(&axes[0]),
                        sample.0[0],
                        cell(&axes[1]),
                        sample.0[1],
                        sample.1
                    )
                    .unwrap();
                }
                for triangle in triangles {
                    writeln!(
                        output,
                        "triangle\t{}\t{}\t{}\t{}",
                        id.0, triangle[0], triangle[1], triangle[2]
                    )
                    .unwrap();
                }
            }
        }
        for (axis, (lo, hi)) in claim.validity.bounds() {
            writeln!(
                output,
                "axis\t{}\t{}\t{lo}\t{hi}\t{:?}",
                id.0,
                cell(axis),
                claim.validity.axis_quantities().get(axis)
            )
            .unwrap();
        }
        for observation in &claim.observations {
            writeln!(output, "claim-observation\t{}\t{}", id.0, observation.0).unwrap();
        }
    }
    for id in claims.observation_ids() {
        if let Some(observation) = claims.observation(id) {
            writeln!(
                output,
                "observation\t{}\t{}\t{}\t{}",
                id.0,
                cell(&observation.specimen),
                cell(&observation.method),
                cell(&observation.caveats)
            )
            .unwrap();
        }
    }
}

fn catalog(root: &Path, options: &Options) -> Result<(String, Vec<String>), String> {
    let directory = root.join("data/matdb/seed-v1");
    let mut paths = Vec::new();
    for entry in
        std::fs::read_dir(&directory).map_err(|e| format!("{}: {e}", directory.display()))?
    {
        let entry = entry.map_err(|e| e.to_string())?;
        if entry.file_type().map_err(|e| e.to_string())?.is_dir() {
            let path = entry.path().join("manifest.tsv");
            if path.is_file() {
                paths.push(path);
            }
        }
    }
    paths.sort();
    let mut output = String::new();
    let mut warnings = Vec::new();
    let mut found = false;
    if options.pack.is_none() {
        output.push_str(
            "source-id\tpack-id\tfamily\tclaims\tproperties\tmanifest\tstatus\tdiagnostic\n",
        );
    }
    for path in paths {
        let source_id = path
            .parent()
            .and_then(Path::file_name)
            .and_then(|s| s.to_str())
            .ok_or_else(|| format!("non-UTF-8 source directory: {}", path.display()))?;
        if options.pack.as_ref().is_some_and(|id| id != source_id) {
            continue;
        }
        let bytes = read_bounded_regular(&path, MAX_MANIFEST_BYTES, "catalog manifest")
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let text = std::str::from_utf8(&bytes).map_err(|e| format!("{}: {e}", path.display()))?;
        let manifest = parse_manifest(text).map_err(|e| format!("{}: {e}", path.display()))?;
        let relative = path.strip_prefix(root).unwrap_or(&path).to_string_lossy();
        let compiled = match compile_manifest(&path) {
            Ok(compiled) => compiled,
            Err(error) => {
                let diagnostic = format!("{}: {error}", relative);
                if options.pack.is_some() {
                    return Err(diagnostic);
                }
                // One existing unsupported source must not hide unrelated usable materials.
                // Metadata can identify a refused source, but cannot satisfy a property filter.
                let metadata = format!("{source_id} {text}").to_lowercase();
                if options.properties.is_empty()
                    && options.search.iter().all(|term| metadata.contains(term))
                {
                    writeln!(
                        output,
                        "{}\t{}\tunavailable\t-\t\t{}\trefused\t{}",
                        cell(source_id),
                        cell(&manifest.pack_id),
                        cell(&relative),
                        cell(&diagnostic)
                    )
                    .unwrap();
                }
                warnings.push(diagnostic);
                continue;
            }
        };
        let claim_set = claims(&compiled.pack);
        let properties: BTreeSet<&str> = claim_set
            .into_iter()
            .flat_map(|set| set.claims_ordered().map(|(_, claim)| claim.key.name()))
            .collect();
        if !options
            .properties
            .iter()
            .all(|name| properties.contains(name.as_str()))
        {
            continue;
        }
        let property_text = properties.into_iter().collect::<Vec<_>>().join(",");
        let mut searchable = format!(
            "{} {} {} {}",
            source_id, manifest.pack_id, manifest.citation, property_text
        );
        if let Some(set) = claim_set {
            for id in set.observation_ids() {
                if let Some(observation) = set.observation(id) {
                    write!(
                        searchable,
                        " {} {}",
                        observation.specimen, observation.caveats
                    )
                    .unwrap();
                }
            }
        }
        let searchable = searchable.to_lowercase();
        if !options.search.iter().all(|term| searchable.contains(term)) {
            continue;
        }
        found = true;
        if options.pack.is_some() {
            writeln!(
                output,
                "source\t{}\npack\t{}\nfamily\t{}\nmanifest\t{}\ncontent-hash\t{}\ncitation\t{}\nlicense\t{}",
                cell(source_id),
                cell(&manifest.pack_id),
                family(&compiled.pack),
                cell(&relative),
                compiled.pack.content_hash(),
                cell(&manifest.citation),
                cell(&manifest.license)
            )
            .unwrap();
            if claim_set.is_some() {
                output.push_str("inspection-scope\tcompiled material claims and observations\n");
                output.push_str(
                    "value-basis\tcoherent SI; quantity descriptors retain semantic kinds\n",
                );
            } else {
                output.push_str("inspection-scope\tmetadata-only; species/model parameters require their dedicated pack reader\n");
            }
            detail(&mut output, &compiled.pack);
        } else {
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\tcompiled\t",
                cell(source_id),
                cell(&manifest.pack_id),
                family(&compiled.pack),
                claim_set.map_or(0, ClaimSet::claim_count),
                cell(&property_text),
                cell(&relative)
            )
            .unwrap();
        }
    }
    if !found && options.pack.is_some() {
        return Err("no source pack matches the exact ID and supplied filters".into());
    }
    Ok((output, warnings))
}

pub(crate) fn cmd_matdb_catalog(root: &Path, args: &[String]) -> Result<(), String> {
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    // Buffer results; an exact-pack refusal cannot leave a successful-looking prefix.
    let (output, warnings) = catalog(root, &options(args)?)?;
    for warning in warnings {
        eprintln!("catalog source refused: {warning}");
    }
    print!("{output}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_refused_source_is_visible_without_hiding_usable_materials() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let source = root.join("data/matdb/seed-v1/ensinger-tecamid-66-natural");
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let scratch =
            std::env::temp_dir().join(format!("matdb-catalog-{}-{nonce}", std::process::id()));
        for directory in ["usable", "refused"] {
            let target = scratch.join("data/matdb/seed-v1").join(directory);
            std::fs::create_dir_all(&target).unwrap();
            std::fs::copy(source.join("manifest.tsv"), target.join("manifest.tsv")).unwrap();
            if directory == "usable" {
                std::fs::copy(source.join("properties.tsv"), target.join("properties.tsv"))
                    .unwrap();
            } else {
                std::fs::write(target.join("properties.tsv"), "unsupported source schema\n")
                    .unwrap();
            }
        }
        let (listing, warnings) = catalog(&scratch, &Options::default()).unwrap();
        assert_eq!(listing.lines().count(), 3);
        assert!(
            listing
                .lines()
                .any(|row| row.starts_with("refused\t") && row.contains("\trefused\t"))
        );
        assert!(
            listing
                .lines()
                .any(|row| row.starts_with("usable\t") && row.ends_with("\tcompiled\t"))
        );
        assert_eq!(warnings.len(), 1);
        let filtered = Options {
            properties: BTreeSet::from(["thermal-conductivity".into()]),
            ..Options::default()
        };
        let (listing, warnings) = catalog(&scratch, &filtered).unwrap();
        assert_eq!(listing.lines().count(), 2);
        assert_eq!(warnings.len(), 1);
        let exact = Options {
            pack: Some("refused".into()),
            ..Options::default()
        };
        assert!(catalog(&scratch, &exact).is_err());
    }

    #[test]
    fn filters_are_explicit_and_tsv_cannot_inject_rows() {
        let parsed = options(&[
            "--search".into(),
            "MDF".into(),
            "--property".into(),
            "density".into(),
        ])
        .unwrap();
        assert_eq!(parsed.search, ["mdf"]);
        assert!(parsed.properties.contains("density"));
        for args in [
            vec!["--search"],
            vec!["--search", "--pack"],
            vec!["--pack", "a", "--pack", "b"],
            vec!["--unknown", "x"],
        ] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
        assert_eq!(cell("a\tb\nc\\d\r"), "a\\tb\\nc\\\\d\\r");
    }
}
