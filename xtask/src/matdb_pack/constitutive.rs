//! Bounded generic constitutive source -> existing FSMODPK v1.
//!
//! Header: `frankensim.constitutive-source.v1`. Tab-separated rows:
//! - `model id law law_version state_version zero|declared bounded|unconstrained`
//! - `parameter model_id name number unit`
//! - `validity model_id axis lower upper unit`
//!
//! Every selector, count and coefficient is dimensioned explicitly. Parameters
//! require linear units; validity endpoints may use affine units. This parser
//! carries arbitrary law identities without claiming to implement their laws;
//! the runtime registry owns parameter physics. No defaults or inferred tissue
//! properties. Duplicate/dangling rows, empty bounded domains and implicit
//! unbounded domains refuse. Parent source/path/licence/byte gates are reused.

use super::*;

pub(super) const PROFILE: &str = "constitutive-v1";
pub(super) const COMPILER_ID: &str = "frankensim-matdb-constitutive-model-pack-compiler-v1";
const HEADER: &str = "frankensim.constitutive-source.v1";
const MODEL_CAP: usize = 128;
const PARAMETER_CAP: usize = 128;
const AXIS_CAP: usize = 32;

struct Model {
    law: String,
    version: u32,
    state: u32,
    initial: InitialStatePolicy,
    bounded: bool,
    source: ContentHash,
}

fn error(subject: &str, detail: impl Into<String>) -> CompileError {
    CompileError::new("constitutive_source_refused", subject, detail)
}

pub(super) fn compile(
    manifest: Manifest,
    sources: &[LoadedSource],
    artifact: ContentHash,
) -> Result<CompileOutput, CompileError> {
    let mut models = BTreeMap::new();
    let mut parameters: BTreeMap<(String, String), ParsedQuantity> = BTreeMap::new();
    let mut axes: BTreeMap<(String, String), (ParsedQuantity, ParsedQuantity)> = BTreeMap::new();
    let mut records = 0usize;
    for source in sources {
        if source.text.as_bytes().contains(&0) {
            return Err(error(&source.spec.id, "NUL byte in source"));
        }
        let mut lines = source.text.lines();
        if lines.next() != Some(HEADER) {
            return Err(error(&source.spec.id, format!("expected {HEADER}")));
        }
        for (offset, line) in lines.enumerate() {
            records += 1;
            let subject = format!("source:{}:line:{}", source.spec.id, offset + 2);
            if line.len() > MAX_LINE_BYTES || records > MAX_RECORDS {
                return Err(error(&subject, "source work cap exceeded"));
            }
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "model" => {
                    require_field_count(&f, 7, &source.spec.id, offset + 2)?;
                    let id = require_identifier(f[1], "model id", &subject)?;
                    let law = require_identifier(f[2], "law id", &subject)?;
                    let version = parse_canonical_u32(f[3], "law version", &subject)?;
                    let state = parse_canonical_u32(f[4], "state version", &subject)?;
                    let initial = match f[5] {
                        "zero" => InitialStatePolicy::ZeroInternalState,
                        "declared" => InitialStatePolicy::RequiresDeclaredState,
                        _ => {
                            return Err(error(
                                &subject,
                                "explicit zero|declared initial-state policy required",
                            ));
                        }
                    };
                    let bounded = match f[6] {
                        "bounded" => true,
                        "unconstrained" => false,
                        _ => {
                            return Err(error(
                                &subject,
                                "explicit bounded|unconstrained validity required",
                            ));
                        }
                    };
                    if models.len() == MODEL_CAP || models.contains_key(&id) {
                        return Err(error(&subject, "duplicate model or model cap exceeded"));
                    }
                    models.insert(
                        id,
                        Model {
                            law,
                            version,
                            state,
                            initial,
                            bounded,
                            source: source.hash,
                        },
                    );
                }
                "parameter" => {
                    require_field_count(&f, 5, &source.spec.id, offset + 2)?;
                    let id = require_identifier(f[1], "model id", &subject)?;
                    let name = require_identifier(f[2], "parameter name", &subject)?;
                    if parameters.len() >= MODEL_CAP * PARAMETER_CAP
                        || parameters.keys().filter(|(model, _)| model == &id).count()
                            == PARAMETER_CAP
                    {
                        return Err(error(&subject, "per-model parameter cap exceeded"));
                    }
                    let (_, _, offset) = unit_transform(f[4], &subject)?;
                    if offset != 0. {
                        return Err(error(
                            &subject,
                            "parameters require linear units (no affine offsets)",
                        ));
                    }
                    let value = parse_linear_quantity(f[3], f[4], &subject)?;
                    if parameters.insert((id, name), value).is_some() {
                        return Err(error(&subject, "duplicate parameter"));
                    }
                }
                "validity" => {
                    require_field_count(&f, 6, &source.spec.id, offset + 2)?;
                    let id = require_identifier(f[1], "model id", &subject)?;
                    let axis = require_identifier(f[2], "axis name", &subject)?;
                    if axes.len() >= MODEL_CAP * AXIS_CAP
                        || axes.keys().filter(|(model, _)| model == &id).count() == AXIS_CAP
                    {
                        return Err(error(&subject, "per-model axis cap exceeded"));
                    }
                    let lower = parse_quantity(f[3], f[5], &subject)?;
                    let upper = parse_quantity(f[4], f[5], &subject)?;
                    if lower.dims != upper.dims || lower.value > upper.value {
                        return Err(error(&subject, "invalid validity endpoints"));
                    }
                    if axes.insert((id, axis), (lower, upper)).is_some() {
                        return Err(error(&subject, "duplicate validity axis"));
                    }
                }
                _ => return Err(error(&subject, "unknown or blank constitutive record")),
            }
        }
    }
    if models.is_empty() {
        return Err(error("sources", "at least one model is required"));
    }
    let repeated_bytes =
        models.len() * (manifest.citation.len() + manifest.license.len() + MAX_IDENTIFIER_BYTES);
    if repeated_bytes > MAX_REPEATED_PROVENANCE_BYTES {
        return Err(error(
            "provenance",
            "repeated citation/licence text exceeds the existing aggregate cap",
        ));
    }
    for (id, _) in parameters.keys().chain(axes.keys()) {
        if !models.contains_key(id) {
            return Err(error(id, "parameter/axis references an undeclared model"));
        }
    }
    let mut cards = Vec::new();
    let mut receipts = Vec::new();
    let mut decisions = vec![Decision::admit(
        "manifest:redistribution",
        "redistribution_permitted",
        "explicit permitted policy and retained source licence admitted by the existing manifest parser",
        Some(artifact),
    )];
    for (id, model) in models {
        let mut validity = ValidityDomain::unconstrained();
        let mut bounds = Vec::new();
        for ((model_id, axis), (lower, upper)) in &axes {
            if model_id == &id {
                validity = validity.with(axis, lower.value, upper.value);
                bounds.push((axis, lower, upper));
            }
        }
        if model.bounded == bounds.is_empty() {
            return Err(error(
                &id,
                "bounded model requires axes; unconstrained model forbids axes",
            ));
        }
        let fields: BTreeMap<_, _> = parameters
            .iter()
            .filter(|((m, _), _)| m == &id)
            .map(|((_, name), value)| {
                (
                    name.clone(),
                    LawParameter {
                        value: value.value,
                        dims: value.dims,
                    },
                )
            })
            .collect();
        if fields.is_empty() {
            return Err(error(&id, "model has no parameters"));
        }
        let card = ConstitutiveModelCard {
            law: LawId(model.law),
            law_version: model.version,
            parameters: fields,
            state_schema_version: model.state,
            initial_state: model.initial,
            validity,
            sources: BTreeSet::from([artifact, model.source])
                .into_iter()
                .collect(),
            provenance: Provenance {
                source: format!("{} [model:{id}]", manifest.citation),
                license: manifest.license.clone(),
                artifact: Some(artifact),
            },
        };
        card.validate().map_err(|e| error(&id, e.to_string()))?;
        let hash = card.content_hash();
        for ((m, name), value) in &parameters {
            if m == &id {
                receipts.push(model_normalization(
                    ModelNormalizationTarget::Parameter {
                        model: hash,
                        parameter: name.clone(),
                    },
                    value,
                ));
            }
        }
        for (axis, lower, upper) in bounds {
            for (side, value) in [
                (ValidityBoundSide::Lower, lower),
                (ValidityBoundSide::Upper, upper),
            ] {
                receipts.push(model_normalization(
                    ModelNormalizationTarget::ValidityBound {
                        model: hash,
                        axis: axis.clone(),
                        side,
                    },
                    value,
                ));
            }
        }
        decisions.push(Decision::admit(format!("model:{id}"),"constitutive_card_admitted","explicit law, state, units and validity transported; runtime registry owns physical admissibility",Some(artifact)));
        cards.push(card);
    }
    let pack = NormalizedModelPack::new(
        manifest.pack_id,
        COMPILER_ID,
        artifact,
        manifest.redistribution_terms,
        cards,
        receipts,
    )
    .map_err(|e| error("pack", e.to_string()))?;
    let bytes = pack.to_bytes();
    let hash = pack.content_hash();
    let decoded = NormalizedModelPack::from_bytes_verified(hash, &bytes)
        .map_err(|e| error("pack", e.to_string()))?;
    if decoded.to_bytes() != bytes {
        return Err(error("pack", "canonical decode changed bytes"));
    }
    for decision in &mut decisions {
        decision.pack_hash = Some(hash);
    }
    sort_decisions(&mut decisions);
    Ok(CompileOutput {
        pack: CompiledPack::Model(pack),
        bytes,
        decisions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const OGDEN: &str = concat!(
        "frankensim.constitutive-source.v1\n",
        "model\tsoft\togden-decoupled\t1\t1\tzero\tbounded\n",
        "parameter\tsoft\tstretch_limit\t3\t1\n",
        "parameter\tsoft\tvolumetric_kind\t0\t1\n",
        "parameter\tsoft\tkappa\t5000\tkPa\n",
        "parameter\tsoft\tterm_count\t1\t1\n",
        "parameter\tsoft\tmu_0\t200\tkPa\n",
        "parameter\tsoft\talpha_0\t2\t1\n",
        "validity\tsoft\tstretch\t0.3333333333333333\t3\t1\n",
    );

    fn run(text: &str) -> Result<CompileOutput, CompileError> {
        let manifest = parse_manifest(&super::super::tests::manifest(PROFILE, true))?;
        let source = LoadedSource {
            spec: manifest.sources[0].clone(),
            bytes: text.as_bytes().to_vec(),
            text: text.to_string(),
            hash: hash_domain(SOURCE_FILE_DOMAIN, text.as_bytes()),
        };
        let artifact = source_envelope_hash(&manifest, std::slice::from_ref(&source));
        compile(manifest, &[source], artifact)
    }

    #[test]
    fn g1_constitutive_source_compiles_to_real_ogden_and_hgo_nodes() {
        use fs_material::graph::LawRegistry;
        use fs_material::hyper::cards::register_hyperelastic_laws;
        use fs_material::hyper::{Hyperelastic, HyperelasticModel};
        let hgo = concat!(
            "frankensim.constitutive-source.v1\n",
            "model\tsoft\thgo-dispersed\t1\t1\tzero\tunconstrained\n",
            "parameter\tsoft\tstretch_limit\t3\t1\n",
            "parameter\tsoft\tvolumetric_kind\t0\t1\n",
            "parameter\tsoft\tkappa\t5000\tkPa\n",
            "parameter\tsoft\tc10\t100\tkPa\n",
            "parameter\tsoft\tk1\t20\tkPa\n",
            "parameter\tsoft\tk2\t2\t1\n",
            "parameter\tsoft\tdispersion\t0\t1\n",
            "parameter\tsoft\tfibre_invariants\t1\t1\n",
            "parameter\tsoft\tfibre_count\t1\t1\n",
            "parameter\tsoft\tfibre_0_x\t1\t1\n",
            "parameter\tsoft\tfibre_0_y\t0\t1\n",
            "parameter\tsoft\tfibre_0_z\t0\t1\n",
        );
        for source in [OGDEN, hgo] {
            let output = run(source).unwrap();
            let CompiledPack::Model(pack) = &output.pack else {
                panic!("expected model pack")
            };
            let decoded =
                NormalizedModelPack::from_bytes_verified(pack.content_hash(), &output.bytes)
                    .unwrap();
            let card = &decoded.models()[0];
            let law = Hyperelastic::from_model_card(card).unwrap();
            if let HyperelasticModel::Ogden { terms, .. } = &law.model {
                assert_eq!(terms[0].mu_pa, 200000.);
            }
            let mut registry = LawRegistry::new();
            register_hyperelastic_laws(&mut registry);
            let node = registry.instantiate("soft", card).unwrap();
            let h = [0.2, 0., 0., 0., 0., 0., 0., 0., 0.];
            let p = node.evaluate(&[], &h).unwrap();
            assert!(p.outputs[0] > 0.);
            assert_eq!(
                p.outputs,
                law.piola(&[1.2, 0., 0., 0., 1., 0., 0., 0., 1.])
                    .unwrap()
                    .to_vec()
            );
            assert_eq!(run(source).unwrap().bytes, output.bytes);
        }
        let (path, _directory) = super::super::tests::write_fixture(
            &super::super::tests::manifest(PROFILE, true),
            OGDEN,
        );
        assert_eq!(
            compile_manifest(&path).unwrap().bytes,
            run(OGDEN).unwrap().bytes
        );
    }

    #[test]
    fn g0_constitutive_source_refuses_ambiguous_or_incomplete_rows() {
        for source in [
            HEADER.to_string() + "\n",
            HEADER.to_string() + "\nmodel\n",
            OGDEN.replace("zero\tbounded", "virgin\tbounded"),
            OGDEN.replace("zero\tbounded", "zero\tunconstrained"),
            OGDEN.replace("mu_0\t200\tkPa", "mu_0\t200\tdegC"),
            OGDEN.replace("parameter\tsoft\tmu_0", "parameter\tunknown\tmu_0"),
            OGDEN.to_string() + "parameter\tsoft\tmu_0\t200\tkPa\n",
            OGDEN.replace("validity\tsoft\tstretch\t0.3333333333333333\t3\t1\n", ""),
            OGDEN.replace("mu_0\t200", "mu_0\tNaN"),
            OGDEN.to_string() + "\n",
        ] {
            assert!(run(&source).is_err(), "{source}");
        }
    }

    #[test]
    fn g0_constitutive_model_parameter_and_axis_caps() {
        let header = "frankensim.constitutive-source.v1\n";
        for count in [MODEL_CAP, MODEL_CAP + 1] {
            let mut source = header.to_string();
            for i in 0..count {
                source += &format!(
                    "model\tm{i}\tsoftware-fixture\t1\t1\tzero\tunconstrained\nparameter\tm{i}\tx\t{i}\t1\n"
                );
            }
            assert_eq!(run(&source).is_ok(), count == MODEL_CAP);
        }
        for count in [PARAMETER_CAP, PARAMETER_CAP + 1] {
            let mut source =
                format!("{header}model\tm\tsoftware-fixture\t1\t1\tzero\tunconstrained\n");
            for i in 0..count {
                source += &format!("parameter\tm\tp{i}\t{i}\t1\n");
            }
            assert_eq!(run(&source).is_ok(), count == PARAMETER_CAP);
        }
        for count in [AXIS_CAP, AXIS_CAP + 1] {
            let mut source = format!(
                "{header}model\tm\tsoftware-fixture\t1\t1\tzero\tbounded\nparameter\tm\tx\t1\t1\n"
            );
            for i in 0..count {
                source += &format!("validity\tm\ta{i}\t0\t1\t1\n");
            }
            assert_eq!(run(&source).is_ok(), count == AXIS_CAP);
        }
    }
}
