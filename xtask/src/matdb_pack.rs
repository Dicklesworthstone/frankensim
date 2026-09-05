//! Strict offline compiler for canonical fs-matdb material, interface, and model packs.
//!
//! This tool boundary owns source-format and redistribution decisions. The
//! runtime fs-matdb crate only admits the normalized artifact. The material
//! profile emits [`NormalizedPack`], while interface-tsv-v1 adds an ordered
//! pair of named surfaces and one complete system context before emitting
//! [`NormalizedInterfacePack`]. NASA-9 and kinetics emit a separate
//! [`NormalizedModelPack`], while species-v1 emits [`NormalizedSpeciesPack`]
//! with the explicit source-declared thermochemical association that generic
//! model cards cannot carry. None of these artifacts authenticates chemistry,
//! supplies an evaluator, or proves a reaction mechanism or conservation law.
//!
//! The bounded manifest is tab-separated: one pack_id, redistribution,
//! citation, and license record, followed by source records containing logical
//! id, safe relative path, and profile. A material-tsv-v1 source starts with
//! SOURCE_HEADER and then uses observation, scalar, curve, uncertainty,
//! validity, frame, and joint records. Interface-tsv-v1 uses that same claim
//! grammar plus exactly one `surface_a`, `surface_b`, and `context` record.
//! Every numeric token has a separate explicit unit field; confidence and
//! correlation use the exact basis 1.
//! Manifest v2 additionally admits explicit `property <source> <target>
//! <kind|dimensional>` and `axis <source> <target>` records. These map names
//! without changing the source grammar. Quantity kinds apply to values;
//! axis aliases do not claim typed validity axes. V1 output stays frozen.
//! Manifest v3 adds an explicit kind to each declared axis. Curve coordinates
//! and validity endpoints share its conversion and retain it in pack v3.
//! A kinetics-v1 source declares exactly one explicitly first-order reaction
//! plus `pre_exponential` and `activation_temperature`, representing only the
//! immutable coefficient schema `k(T) = A exp(-T_a / T)` with `A` in `s^-1`
//! and `T_a` in kelvin. No evaluator, stoichiometry, or mechanism is inferred.
//! A species-v1 source contains exactly one `species` record with, in order,
//! canonical `SpeciesId`, molar-mass value/unit, standard-state phase,
//! reference EOS, reference-pressure value/unit, and an opaque elemental-
//! reference convention id. Version 1 admits only gas/ideal-gas associations.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{Read as _, Write as _};
use std::path::{Component, Path, PathBuf};

#[cfg(unix)]
use std::os::unix::fs::MetadataExt as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    JointStatistics, LawId, LawParameter, MATDB_PACK_TARGET_BASIS, MODEL_PACK_TARGET_BASIS,
    MaterialStateId, ModelNormalizationReceipt, ModelNormalizationTarget, NormalizationReceipt,
    NormalizationTarget, NormalizedInterfacePack, NormalizedModelPack, NormalizedPack,
    NormalizedSpeciesPack, ObservationDataset, ObservationId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SPECIES_MOLAR_MASS_DIMS, SPECIES_PACK_TARGET_BASIS,
    SPECIES_REFERENCE_PRESSURE_DIMS, SpeciesAssociation, SpeciesNormalizationReceipt,
    SpeciesNormalizationTarget, StatisticMember, SurfaceSpec, SystemContext, UncertaintyModel,
    ValidityBoundSide,
};
use fs_qty::chemistry::{ReactionId, SpeciesId};
use fs_qty::parse::{ParseBudget, parse_qty_with_budget};
use fs_qty::semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v1";
const INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v1";
const MAPPED_COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v2";
const MAPPED_INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v2";
const TYPED_AXIS_COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v3";
const TYPED_AXIS_INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v3";
const HARDNESS_COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v4";
const HARDNESS_INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v4";
const NASA9_COMPILER_ID: &str = "frankensim-matdb-nasa9-model-pack-compiler-v1";
const KINETICS_COMPILER_ID: &str = "frankensim-matdb-kinetics-model-pack-compiler-v1";
const SPECIES_COMPILER_ID: &str = "frankensim-matdb-species-pack-compiler-v1";
/// Semantic contract version for normalized material-pack compilation.
///
/// Bump this whenever parsing, admission, normalization, or provenance
/// semantics can change the canonical compiler fixture.
#[allow(dead_code)] // consumed textually by `xtask check-goldens`
pub const MATDB_PACK_COMPILER_SEMANTICS_VERSION: u32 = 7;
const MANIFEST_HEADER: &str = "frankensim.matdb-manifest.v1";
const MAPPED_MANIFEST_HEADER: &str = "frankensim.matdb-manifest.v2";
const TYPED_AXIS_MANIFEST_HEADER: &str = "frankensim.matdb-manifest.v3";
const HARDNESS_MANIFEST_HEADER: &str = "frankensim.matdb-manifest.v4";
const SOURCE_HEADER: &str = "frankensim.matdb-source.v1";
const NASA9_SOURCE_HEADER: &str = "frankensim.nasa9-source.v1";
const KINETICS_SOURCE_HEADER: &str = "frankensim.kinetics-source.v1";
const SPECIES_SOURCE_HEADER: &str = "frankensim.species-source.v1";
const MATERIAL_PROFILE: &str = "material-tsv-v1";
const INTERFACE_PROFILE: &str = "interface-tsv-v1";
const NASA9_PROFILE: &str = "nasa9-v1";
const KINETICS_PROFILE: &str = "kinetics-v1";
const SPECIES_PROFILE: &str = "species-v1";
const NASA9_LAW_ID: &str = "nasa9-standard-state";
const NASA9_LAW_VERSION: u32 = 1;
const NASA9_STATE_SCHEMA_VERSION: u32 = 0;
const NASA9_COEFFICIENT_NAMES: [&str; 9] = ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8"];
const KINETICS_LAW_ID: &str = "arrhenius-first-order-rate";
const KINETICS_LAW_VERSION: u32 = 1;
const KINETICS_STATE_SCHEMA_VERSION: u32 = 0;
const KINETICS_PARAMETER_NAMES: [&str; 2] = ["activation_temperature", "pre_exponential"];
const SOURCE_ENVELOPE_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-envelope.v1";
const SOURCE_FILE_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-file.v1";
const SOURCE_RECORD_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-record.v1";
const SOURCE_LITERAL_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-literal.v1";
const DECISION_ID_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.decision.v1";
const MAX_MANIFEST_BYTES: usize = 1_048_576;
const MAX_SOURCE_BYTES: usize = 64 * 1_048_576;
const MAX_TOTAL_SOURCE_BYTES: usize = 256 * 1_048_576;
const MAX_SOURCES: usize = 256;
const MAX_RECORDS: usize = 200_000;
const MAX_LINE_BYTES: usize = 1_048_576;
const MAX_IDENTIFIER_BYTES: usize = 128;
const MAX_TEXT_BYTES: usize = 1_048_576;
const MAX_CURVE_KNOTS: usize = 1_000_000;
const MAX_JOINT_MEMBERS: usize = 256;
const MAX_NORMALIZATION_RECEIPTS: usize = 100_000;
const MAX_REPEATED_PROVENANCE_BYTES: usize = 64 * 1_048_576;
const MAX_NASA9_REGIONS: usize = 16;
const MAX_NASA9_COEFFICIENTS: usize = MAX_NASA9_REGIONS * NASA9_COEFFICIENT_NAMES.len();
const MAX_KINETICS_REACTIONS: usize = 1;
const MAX_KINETICS_PARAMETERS: usize = KINETICS_PARAMETER_NAMES.len();
const QTY_BUDGET: ParseBudget = ParseBudget::new(4_096, 256, 64, 256);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CompileError {
    compiler_id: String,
    code: &'static str,
    subject: String,
    detail: String,
    input_hash: Option<ContentHash>,
    prior_decisions: Vec<Decision>,
}

impl CompileError {
    fn new(code: &'static str, subject: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            compiler_id: COMPILER_ID.to_string(),
            code,
            subject: subject.into(),
            detail: detail.into(),
            input_hash: None,
            prior_decisions: Vec::new(),
        }
    }

    fn with_compiler_id(mut self, compiler_id: impl Into<String>) -> Self {
        self.compiler_id = compiler_id.into();
        self
    }

    fn with_input_hash(mut self, input_hash: ContentHash) -> Self {
        if self.input_hash.is_none() {
            self.input_hash = Some(input_hash);
        }
        self
    }

    fn with_prior_decisions(mut self, mut decisions: Vec<Decision>) -> Self {
        if !self.prior_decisions.is_empty() {
            decisions.append(&mut self.prior_decisions);
        }
        sort_decisions(&mut decisions);
        self.prior_decisions = decisions;
        self
    }
}

impl std::fmt::Display for CompileError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "matdb pack refused [{}] {}: {}",
            self.code, self.subject, self.detail
        )
    }
}

impl std::error::Error for CompileError {}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Decision {
    subject: String,
    verdict: &'static str,
    reason_code: &'static str,
    detail: String,
    source_hash: Option<ContentHash>,
    pack_hash: Option<ContentHash>,
}

impl Decision {
    fn admit(
        subject: impl Into<String>,
        reason_code: &'static str,
        detail: impl Into<String>,
        source_hash: Option<ContentHash>,
    ) -> Self {
        Self {
            subject: subject.into(),
            verdict: "admit",
            reason_code,
            detail: detail.into(),
            source_hash,
            pack_hash: None,
        }
    }

    fn refusal(error: &CompileError, input_hash: Option<ContentHash>) -> Self {
        Self {
            subject: error.subject.clone(),
            verdict: "refuse",
            reason_code: error.code,
            detail: error.detail.clone(),
            source_hash: input_hash,
            pack_hash: None,
        }
    }

    fn refusal_for_pack(error: &CompileError, pack_hash: ContentHash) -> Self {
        let mut decision = Self::refusal(error, None);
        decision.pack_hash = Some(pack_hash);
        decision
    }

    fn canonical_preimage(&self, compiler_id: &str) -> Vec<u8> {
        let mut bytes = Vec::new();
        for part in [
            compiler_id.as_bytes(),
            self.subject.as_bytes(),
            self.verdict.as_bytes(),
            self.reason_code.as_bytes(),
            self.detail.as_bytes(),
        ] {
            push_part(&mut bytes, part);
        }
        push_part(&mut bytes, b"source_hash");
        push_part(&mut bytes, &[u8::from(self.source_hash.is_some())]);
        if let Some(hash) = self.source_hash {
            push_part(&mut bytes, hash.as_bytes());
        }
        push_part(&mut bytes, b"pack_hash");
        push_part(&mut bytes, &[u8::from(self.pack_hash.is_some())]);
        if let Some(hash) = self.pack_hash {
            push_part(&mut bytes, hash.as_bytes());
        }
        bytes
    }
}

fn refusal_decisions(error: &CompileError, input_hash: Option<ContentHash>) -> Vec<Decision> {
    let mut decisions = error.prior_decisions.clone();
    decisions.push(Decision::refusal(error, input_hash));
    sort_decisions(&mut decisions);
    decisions
}

fn print_refusal_decisions(error: &CompileError, input_hash: Option<ContentHash>) {
    for decision in refusal_decisions(error, input_hash) {
        println!("{}", render_decision(&error.compiler_id, &decision));
    }
}

#[derive(Debug)]
struct CompileOutput {
    pack: CompiledPack,
    bytes: Vec<u8>,
    decisions: Vec<Decision>,
}

#[derive(Debug)]
enum CompiledPack {
    Material(NormalizedPack),
    Interface(NormalizedInterfacePack),
    Model(NormalizedModelPack),
    Species(NormalizedSpeciesPack),
}

impl CompiledPack {
    fn compiler_id(&self) -> &str {
        match self {
            Self::Material(pack) => pack.compiler(),
            Self::Interface(pack) => pack.compiler(),
            Self::Model(pack) => pack.compiler(),
            Self::Species(pack) => pack.compiler(),
        }
    }

    fn source_artifact(&self) -> ContentHash {
        match self {
            Self::Material(pack) => pack.source_artifact(),
            Self::Interface(pack) => pack.source_artifact(),
            Self::Model(pack) => pack.source_artifact(),
            Self::Species(pack) => pack.source_artifact(),
        }
    }

    fn content_hash(&self) -> ContentHash {
        match self {
            Self::Material(pack) => pack.content_hash(),
            Self::Interface(pack) => pack.content_hash(),
            Self::Model(pack) => pack.content_hash(),
            Self::Species(pack) => pack.content_hash(),
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        match self {
            Self::Material(pack) => pack.to_bytes(),
            Self::Interface(pack) => pack.to_bytes(),
            Self::Model(pack) => pack.to_bytes(),
            Self::Species(pack) => pack.to_bytes(),
        }
    }

    fn verify_bytes(&self, expected: ContentHash, bytes: &[u8]) -> Result<(), String> {
        match self {
            Self::Material(_) => NormalizedPack::from_bytes_verified(expected, bytes)
                .map(|_| ())
                .map_err(|error| error.to_string()),
            Self::Interface(_) => NormalizedInterfacePack::from_bytes_verified(expected, bytes)
                .map(|_| ())
                .map_err(|error| error.to_string()),
            Self::Model(_) => NormalizedModelPack::from_bytes_verified(expected, bytes)
                .map(|_| ())
                .map_err(|error| error.to_string()),
            Self::Species(_) => NormalizedSpeciesPack::from_bytes_verified(expected, bytes)
                .map(|_| ())
                .map_err(|error| error.to_string()),
        }
    }

    #[cfg(test)]
    fn material(&self) -> &NormalizedPack {
        let Self::Material(pack) = self else {
            panic!("expected a material pack");
        };
        pack
    }

    #[cfg(test)]
    fn model(&self) -> &NormalizedModelPack {
        let Self::Model(pack) = self else {
            panic!("expected a model pack");
        };
        pack
    }

    #[cfg(test)]
    fn species(&self) -> &NormalizedSpeciesPack {
        let Self::Species(pack) = self else {
            panic!("expected a species pack");
        };
        pack
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceProfile {
    Material,
    Interface,
    Nasa9,
    Kinetics,
    Species,
}

impl SourceProfile {
    const fn compiler_id(self) -> &'static str {
        match self {
            Self::Material => COMPILER_ID,
            Self::Interface => INTERFACE_COMPILER_ID,
            Self::Nasa9 => NASA9_COMPILER_ID,
            Self::Kinetics => KINETICS_COMPILER_ID,
            Self::Species => SPECIES_COMPILER_ID,
        }
    }
}

#[derive(Debug)]
struct Manifest {
    version: u32,
    pack_id: String,
    redistribution_terms: String,
    citation: String,
    license: String,
    sources: Vec<SourceSpec>,
    properties: BTreeMap<String, (String, Option<QuantityKind>)>,
    axes: BTreeMap<String, String>,
    axis_kinds: BTreeMap<String, Option<QuantityKind>>,
    hardness_tests: BTreeMap<String, HardnessSourceTest>,
}

#[derive(Debug)]
struct HardnessSourceTest {
    observation: String,
    indenter: String,
    protocol: String,
    // Source-spelled force, unit, dwell and unit; order is physical history.
    holds: Vec<[String; 4]>,
}

impl Manifest {
    fn axis<'a>(&'a self, source: &'a str) -> &'a str {
        self.axes.get(source).map_or(source, String::as_str)
    }
}

#[derive(Debug, Clone)]
struct SourceSpec {
    id: String,
    relative_path: PathBuf,
    profile: String,
}

#[derive(Debug)]
struct LoadedSource {
    spec: SourceSpec,
    bytes: Vec<u8>,
    text: String,
    hash: ContentHash,
}

#[derive(Debug)]
struct RawObservation {
    id: String,
    specimen: String,
    method: String,
    caveats: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
enum RawClaimValue {
    Scalar {
        number: String,
        unit: String,
    },
    Curve {
        abscissa: String,
        abscissa_unit: String,
        ordinate_unit: String,
        points: Vec<(String, String)>,
    },
}

#[derive(Debug)]
struct RawClaim {
    id: String,
    observations: Vec<String>,
    property: String,
    value: RawClaimValue,
    interpolation: InterpolationPolicy,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
enum RawUncertaintyKind {
    Unstated,
    Absolute {
        number: String,
        unit: String,
        confidence: String,
        confidence_unit: String,
    },
    Relative {
        number: String,
        unit: String,
        confidence: String,
        confidence_unit: String,
    },
}

#[derive(Debug)]
struct RawUncertainty {
    kind: RawUncertaintyKind,
    source_hash: ContentHash,
}

#[derive(Debug)]
struct RawValidity {
    axis: String,
    lower: String,
    upper: String,
    unit: String,
    source_hash: ContentHash,
}

#[derive(Debug, Clone)]
struct RawFrame {
    source: String,
    target: String,
    source_hash: ContentHash,
}

#[derive(Debug)]
struct RawJoint {
    observation: String,
    block_id: String,
    members: Vec<String>,
    covariance: Vec<String>,
    correlation: Option<Vec<String>>,
    source_hash: ContentHash,
}

#[derive(Debug)]
struct RawSurfaceSpec {
    chemistry: String,
    phase: String,
    process: String,
    revision: u32,
    texture_frame: String,
    record_hash: ContentHash,
}

#[derive(Debug)]
struct RawSystemContext {
    medium: String,
    third_body: Option<String>,
    environment: String,
    history: String,
    record_hash: ContentHash,
}

#[derive(Debug, Default)]
struct RawDatabase {
    observations: BTreeMap<String, RawObservation>,
    claims: BTreeMap<String, RawClaim>,
    uncertainties: BTreeMap<String, RawUncertainty>,
    validities: BTreeMap<String, BTreeMap<String, RawValidity>>,
    frames: BTreeMap<String, RawFrame>,
    joints: BTreeMap<(String, String), RawJoint>,
    surface_a: Option<RawSurfaceSpec>,
    surface_b: Option<RawSurfaceSpec>,
    context: Option<RawSystemContext>,
    decisions: Vec<Decision>,
    records: usize,
    curve_knots: usize,
}

#[derive(Debug)]
struct RawNasa9Region {
    species: SpeciesId,
    id: String,
    temperature_lower: String,
    temperature_upper: String,
    temperature_unit: String,
    reference_pressure: String,
    reference_pressure_unit: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
struct RawNasa9Coefficient {
    value: String,
    unit: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug, Default)]
struct RawNasa9Database {
    regions: BTreeMap<(SpeciesId, String), RawNasa9Region>,
    coefficients: BTreeMap<(SpeciesId, String, String), RawNasa9Coefficient>,
    decisions: Vec<Decision>,
    records: usize,
}

#[derive(Debug)]
struct RawKineticsReaction {
    reaction: ReactionId,
    temperature_lower: String,
    temperature_upper: String,
    temperature_unit: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
struct RawKineticsParameter {
    value: String,
    unit: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug, Default)]
struct RawKineticsDatabase {
    reactions: BTreeMap<ReactionId, RawKineticsReaction>,
    parameters: BTreeMap<(ReactionId, String), RawKineticsParameter>,
    decisions: Vec<Decision>,
    records: usize,
}

#[derive(Debug)]
struct RawSpeciesAssociation {
    species: SpeciesId,
    molar_mass: String,
    molar_mass_unit: String,
    standard_state_phase: String,
    reference_eos: String,
    reference_pressure: String,
    reference_pressure_unit: String,
    elemental_reference: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug, Default)]
struct RawSpeciesDatabase {
    association: Option<RawSpeciesAssociation>,
    decisions: Vec<Decision>,
    records: usize,
}

#[derive(Debug, Clone)]
struct ParsedQuantity {
    value: f64,
    dims: Dims,
    scale: f64,
    offset: f64,
    unit: String,
    literal_hash: ContentHash,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum LocalComponent {
    Scalar,
    CurveAbscissa(u32),
    CurveOrdinate(u32),
}

#[derive(Debug, Clone)]
struct ComponentMeta {
    member: StatisticMember,
    dims: Dims,
    scale: f64,
    unit: String,
    observations: BTreeSet<ObservationId>,
}

/// Compile a strict source manifest twice and publish one verified new artifact.
pub(super) fn cmd_matdb_pack(root: &Path, args: &[String]) -> Result<(), String> {
    let (manifest_path, output_path) = match parse_cli(root, args) {
        Ok(paths) => paths,
        Err(error) => {
            print_refusal_decisions(&error, Some(cli_input_hash(args)));
            return Err(error.to_string());
        }
    };
    let mut first = match compile_manifest(&manifest_path) {
        Ok(output) => output,
        Err(error) => {
            print_refusal_decisions(&error, error.input_hash);
            return Err(error.to_string());
        }
    };
    let second = match compile_manifest(&manifest_path) {
        Ok(output) => output,
        Err(error) => {
            print_refusal_decisions(&error, error.input_hash);
            return Err(format!("second deterministic pass failed: {error}"));
        }
    };
    let first_source = first.pack.source_artifact();
    let second_source = second.pack.source_artifact();
    if first_source != second_source {
        let mut pair = Vec::new();
        push_part(&mut pair, b"first_source_artifact");
        push_part(&mut pair, first_source.as_bytes());
        push_part(&mut pair, b"second_source_artifact");
        push_part(&mut pair, second_source.as_bytes());
        let error = CompileError::new(
            "input_changed_between_passes",
            "sources",
            "source envelope changed between the two complete compilation passes",
        )
        .with_compiler_id(first.pack.compiler_id())
        .with_input_hash(hash_domain(SOURCE_ENVELOPE_DOMAIN, &pair));
        println!(
            "{}",
            render_decision(
                &error.compiler_id,
                &Decision::refusal(&error, error.input_hash)
            )
        );
        return Err(error.to_string());
    }
    if first.pack.compiler_id() != second.pack.compiler_id()
        || first.bytes != second.bytes
        || first.decisions != second.decisions
    {
        let error = CompileError::new(
            "nondeterministic_compilation",
            "pack",
            "two complete source-to-pack passes produced different artifacts or decisions",
        )
        .with_compiler_id(first.pack.compiler_id())
        .with_input_hash(first_source);
        println!(
            "{}",
            render_decision(
                &error.compiler_id,
                &Decision::refusal(&error, error.input_hash,)
            )
        );
        return Err(error.to_string());
    }

    let compiler_id = first.pack.compiler_id().to_string();
    let pack_hash = first.pack.content_hash();
    write_new_verified(&output_path, &first.bytes, &first.pack).map_err(|error| {
        let error = error.with_compiler_id(compiler_id.clone());
        println!(
            "{}",
            render_decision(
                &error.compiler_id,
                &Decision::refusal_for_pack(&error, pack_hash)
            )
        );
        error.to_string()
    })?;
    first.decisions.push(Decision {
        subject: "output".to_string(),
        verdict: "admit",
        reason_code: "published_new_verified_artifact",
        detail: "created one previously absent output and re-admitted its exact bytes".to_string(),
        source_hash: None,
        pack_hash: Some(pack_hash),
    });
    sort_decisions(&mut first.decisions);
    for decision in &first.decisions {
        println!("{}", render_decision(&compiler_id, decision));
    }
    eprintln!(
        "matdb pack OK: {} bytes, hash {}, {} decisions",
        first.bytes.len(),
        pack_hash,
        first.decisions.len()
    );
    Ok(())
}

fn parse_cli(root: &Path, args: &[String]) -> Result<(PathBuf, PathBuf), CompileError> {
    let mut manifest = None;
    let mut output = None;
    let mut index = 0;
    while index < args.len() {
        let option = &args[index];
        let value = args.get(index + 1).ok_or_else(|| {
            CompileError::new(
                "invalid_cli",
                "arguments",
                format!("option {option:?} requires one path"),
            )
        })?;
        let slot = match option.as_str() {
            "--manifest" => &mut manifest,
            "--out" => &mut output,
            _ => {
                return Err(CompileError::new(
                    "invalid_cli",
                    "arguments",
                    format!("unknown option {option:?}; expected --manifest and --out"),
                ));
            }
        };
        if slot.replace(resolve_cli_path(root, value)).is_some() {
            return Err(CompileError::new(
                "invalid_cli",
                "arguments",
                format!("option {option:?} was supplied more than once"),
            ));
        }
        index += 2;
    }
    let manifest = manifest.ok_or_else(|| {
        CompileError::new(
            "invalid_cli",
            "arguments",
            "missing required --manifest path",
        )
    })?;
    let output = output.ok_or_else(|| {
        CompileError::new("invalid_cli", "arguments", "missing required --out path")
    })?;
    if manifest == output {
        return Err(CompileError::new(
            "invalid_cli",
            "arguments",
            "manifest and output paths must differ",
        ));
    }
    Ok((manifest, output))
}

fn resolve_cli_path(root: &Path, raw: &str) -> PathBuf {
    let path = PathBuf::from(raw);
    if path.is_absolute() {
        path
    } else {
        root.join(path)
    }
}

fn cli_input_hash(args: &[String]) -> ContentHash {
    let mut payload = Vec::new();
    for argument in args {
        push_part(&mut payload, argument.as_bytes());
    }
    hash_domain(SOURCE_ENVELOPE_DOMAIN, &payload)
}

fn parse_manifest(text: &str) -> Result<Manifest, CompileError> {
    if text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_manifest_encoding",
            "manifest",
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = text.lines();
    let version = match lines.next() {
        Some(MANIFEST_HEADER) => 1,
        Some(MAPPED_MANIFEST_HEADER) => 2,
        Some(TYPED_AXIS_MANIFEST_HEADER) => 3,
        Some(HARDNESS_MANIFEST_HEADER) => 4,
        _ => {
            return Err(CompileError::new(
                "unsupported_manifest_schema",
                "manifest",
                format!(
                    "first line must be {MANIFEST_HEADER:?}, {MAPPED_MANIFEST_HEADER:?}, {TYPED_AXIS_MANIFEST_HEADER:?}, or {HARDNESS_MANIFEST_HEADER:?}"
                ),
            ));
        }
    };
    let mut pack_id = None;
    let mut redistribution_terms = None;
    let mut citation = None;
    let mut license = None;
    let mut sources = Vec::new();
    let mut properties = BTreeMap::new();
    let mut axes = BTreeMap::new();
    let mut axis_kinds = BTreeMap::new();
    let mut hardness_tests: BTreeMap<String, HardnessSourceTest> = BTreeMap::new();
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                format!("manifest:line:{line_number}"),
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_manifest_record",
                format!("manifest:line:{line_number}"),
                "blank records are forbidden",
            ));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("hardness-test") if version >= 4 => {
                require_field_count(&fields, 5, "manifest", line_number)?;
                let claim = require_identifier(fields[1], "hardness claim", "manifest")?;
                let test = HardnessSourceTest {
                    observation: require_identifier(fields[2], "hardness observation", "manifest")?,
                    indenter: require_text(fields[3], "hardness indenter", "manifest")?,
                    protocol: require_text(fields[4], "hardness protocol", "manifest")?,
                    holds: Vec::new(),
                };
                if hardness_tests.insert(claim, test).is_some() {
                    return Err(CompileError::new(
                        "duplicate_hardness_test",
                        "manifest",
                        "one claim may declare only one hardness test",
                    ));
                }
            }
            Some("hardness-hold") if version >= 4 => {
                require_field_count(&fields, 6, "manifest", line_number)?;
                let test = hardness_tests.get_mut(fields[1]).ok_or_else(|| {
                    CompileError::new(
                        "unknown_hardness_test",
                        "manifest",
                        "declare hardness-test before its ordered holds",
                    )
                })?;
                if test.holds.len() == 64 {
                    return Err(CompileError::new(
                        "resource_limit",
                        "hardness-hold",
                        "at most 64 force/dwell holds are admitted",
                    ));
                }
                test.holds.push([
                    fields[2].into(),
                    fields[3].into(),
                    fields[4].into(),
                    fields[5].into(),
                ]);
            }
            Some("property") if version >= 2 => {
                require_field_count(&fields, 4, "manifest", line_number)?;
                let source = require_identifier(fields[1], "source property", "manifest")?;
                let target = require_identifier(fields[2], "target property", "manifest")?;
                let kind = manifest_quantity_kind(fields[3])?;
                if properties.insert(source, (target, kind)).is_some() {
                    return Err(CompileError::new(
                        "duplicate_property_mapping",
                        "manifest",
                        "one source property may be mapped only once",
                    ));
                }
            }
            Some("axis") if version >= 2 => {
                require_field_count(
                    &fields,
                    if version >= 3 { 4 } else { 3 },
                    "manifest",
                    line_number,
                )?;
                let source = require_identifier(fields[1], "source axis", "manifest")?;
                let target = require_identifier(fields[2], "target axis", "manifest")?;
                if version >= 3 {
                    axis_kinds.insert(source.clone(), manifest_quantity_kind(fields[3])?);
                }
                if axes.insert(source, target).is_some() {
                    return Err(CompileError::new(
                        "duplicate_axis_mapping",
                        "manifest",
                        "one source axis may be mapped only once",
                    ));
                }
            }
            Some("pack_id") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut pack_id,
                    require_text(fields[1], "pack_id", "manifest")?,
                    "duplicate_pack_id",
                    "manifest",
                )?;
            }
            Some("redistribution") => {
                require_field_count(&fields, 3, "manifest", line_number)?;
                if fields[1] != "permitted" {
                    return Err(CompileError::new(
                        "redistribution_refused",
                        "manifest",
                        "redistribution policy must be explicitly permitted",
                    ));
                }
                set_once(
                    &mut redistribution_terms,
                    require_text(fields[2], "redistribution terms", "manifest")?,
                    "duplicate_redistribution",
                    "manifest",
                )?;
            }
            Some("citation") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut citation,
                    require_text(fields[1], "citation", "manifest")?,
                    "duplicate_citation",
                    "manifest",
                )?;
            }
            Some("license") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut license,
                    require_text(fields[1], "license", "manifest")?,
                    "duplicate_license",
                    "manifest",
                )?;
            }
            Some("source") => {
                require_field_count(&fields, 4, "manifest", line_number)?;
                if sources.len() == MAX_SOURCES {
                    return Err(CompileError::new(
                        "resource_limit",
                        "manifest:sources",
                        format!("at most {MAX_SOURCES} source files are admitted"),
                    ));
                }
                let id = require_identifier(fields[1], "source id", "manifest")?;
                let relative_path = safe_relative_path(fields[2])?;
                let profile = require_text(fields[3], "source profile", "manifest")?;
                sources.push(SourceSpec {
                    id,
                    relative_path,
                    profile,
                });
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_manifest_record",
                    format!("manifest:line:{line_number}"),
                    format!("unknown record type {other:?}"),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    if sources.is_empty() {
        return Err(CompileError::new(
            "missing_sources",
            "manifest",
            "at least one source record is required",
        ));
    }
    sources.sort_by(|left, right| left.id.cmp(&right.id));
    let mut paths = BTreeSet::new();
    for pair in sources.windows(2) {
        if pair[0].id == pair[1].id {
            return Err(CompileError::new(
                "duplicate_source_id",
                format!("source:{}", pair[0].id),
                "source ids must be unique",
            ));
        }
    }
    for source in &sources {
        if !paths.insert(source.relative_path.clone()) {
            return Err(CompileError::new(
                "duplicate_source_path",
                format!("source:{}", source.id),
                "one physical source path may not be admitted under two ids",
            ));
        }
    }
    Ok(Manifest {
        version,
        properties,
        axes,
        pack_id: pack_id.ok_or_else(|| {
            CompileError::new("missing_pack_id", "manifest", "pack_id record is required")
        })?,
        redistribution_terms: redistribution_terms.ok_or_else(|| {
            CompileError::new(
                "missing_redistribution_terms",
                "manifest",
                "an explicit permitted redistribution record is required",
            )
        })?,
        citation: citation.ok_or_else(|| {
            CompileError::new(
                "missing_citation",
                "manifest",
                "citation record is required",
            )
        })?,
        license: license.ok_or_else(|| {
            CompileError::new("missing_license", "manifest", "license record is required")
        })?,
        sources,
        axis_kinds,
        hardness_tests,
    })
}

fn manifest_quantity_kind(name: &str) -> Result<Option<QuantityKind>, CompileError> {
    Ok(match name {
        "dimensional" => None,
        "absolute-temperature" => Some(QuantityKind::AbsoluteTemperature),
        "temperature-difference" => Some(QuantityKind::TemperatureDifference),
        "energy" => Some(QuantityKind::Energy),
        "torque" => Some(QuantityKind::Torque),
        name => Some(
            QuantityKind::from_material_convention_name(name).ok_or_else(|| {
                CompileError::new(
                    "unknown_quantity_kind",
                    "manifest",
                    format!("unsupported quantity kind {name:?}"),
                )
            })?,
        ),
    })
}

fn source_profile(manifest: &Manifest) -> Result<SourceProfile, CompileError> {
    let Some((first, rest)) = manifest.sources.split_first() else {
        return Err(CompileError::new(
            "missing_sources",
            "manifest",
            "at least one source record is required",
        ));
    };
    for source in rest {
        if source.profile != first.profile {
            return Err(CompileError::new(
                "mixed_source_profiles",
                format!("source:{}", source.id),
                format!(
                    "one output may contain exactly one artifact profile; source {:?} uses {:?} while source {:?} uses {:?}",
                    first.id, first.profile, source.id, source.profile
                ),
            ));
        }
    }
    if manifest.version >= 2
        && !matches!(first.profile.as_str(), MATERIAL_PROFILE | INTERFACE_PROFILE)
    {
        return Err(CompileError::new(
            "unsupported_manifest_profile",
            "manifest",
            "v2/v3/v4 declarations apply only to material and interface sources",
        ));
    }
    match first.profile.as_str() {
        MATERIAL_PROFILE => Ok(SourceProfile::Material),
        INTERFACE_PROFILE => Ok(SourceProfile::Interface),
        NASA9_PROFILE => Ok(SourceProfile::Nasa9),
        KINETICS_PROFILE => Ok(SourceProfile::Kinetics),
        SPECIES_PROFILE => Ok(SourceProfile::Species),
        _ => Err(CompileError::new(
            "unsupported_source_profile",
            format!("source:{}", first.id),
            format!(
                "profile {:?} has no runtime-loadable compiler path; supported profiles are {MATERIAL_PROFILE:?}, {INTERFACE_PROFILE:?}, {NASA9_PROFILE:?}, {KINETICS_PROFILE:?}, and {SPECIES_PROFILE:?}",
                first.profile
            ),
        )),
    }
}

fn load_sources(
    manifest_path: &Path,
    manifest: &Manifest,
) -> Result<Vec<LoadedSource>, CompileError> {
    let parent = manifest_path.parent().unwrap_or_else(|| Path::new("."));
    let canonical_parent = parent.canonicalize().map_err(|error| {
        CompileError::new(
            "source_root_unavailable",
            "manifest",
            format!("cannot resolve manifest directory: {error}"),
        )
    })?;
    let mut loaded = Vec::with_capacity(manifest.sources.len());
    let mut total = 0usize;
    for spec in &manifest.sources {
        let path = resolve_safe_source(&canonical_parent, &spec.relative_path, &spec.id)?;
        let bytes = read_bounded_regular(&path, MAX_SOURCE_BYTES, &format!("source:{}", spec.id))?;
        let hash = hash_domain(SOURCE_FILE_DOMAIN, &bytes);
        total = total.checked_add(bytes.len()).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                "sources",
                "aggregate source byte count overflowed",
            )
            .with_input_hash(hash)
        })?;
        if total > MAX_TOTAL_SOURCE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                "sources",
                format!("aggregate source bytes exceed {MAX_TOTAL_SOURCE_BYTES}"),
            )
            .with_input_hash(hash));
        }
        let text = String::from_utf8(bytes.clone()).map_err(|_| {
            CompileError::new(
                "invalid_source_encoding",
                format!("source:{}", spec.id),
                "source must be valid UTF-8",
            )
            .with_input_hash(hash)
        })?;
        loaded.push(LoadedSource {
            spec: spec.clone(),
            bytes,
            text,
            hash,
        });
    }
    Ok(loaded)
}

fn resolve_safe_source(root: &Path, relative: &Path, id: &str) -> Result<PathBuf, CompileError> {
    let mut cursor = root.to_path_buf();
    for component in relative.components() {
        let Component::Normal(part) = component else {
            return Err(CompileError::new(
                "unsafe_source_path",
                format!("source:{id}"),
                "source path must contain only ordinary relative components",
            ));
        };
        cursor.push(part);
        let metadata = fs::symlink_metadata(&cursor).map_err(|error| {
            CompileError::new(
                "source_unavailable",
                format!("source:{id}"),
                format!("cannot inspect source path component: {error}"),
            )
        })?;
        if metadata.file_type().is_symlink() {
            return Err(CompileError::new(
                "unsafe_source_path",
                format!("source:{id}"),
                "symlinks are forbidden in source paths",
            ));
        }
    }
    let canonical = cursor.canonicalize().map_err(|error| {
        CompileError::new(
            "source_unavailable",
            format!("source:{id}"),
            format!("cannot resolve source file: {error}"),
        )
    })?;
    if !canonical.starts_with(root) {
        return Err(CompileError::new(
            "unsafe_source_path",
            format!("source:{id}"),
            "resolved source escapes the manifest directory",
        ));
    }
    Ok(canonical)
}

fn read_bounded_regular(path: &Path, limit: usize, subject: &str) -> Result<Vec<u8>, CompileError> {
    let file = File::open(path).map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot open input: {error}"),
        )
    })?;
    let opened = file.metadata().map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot inspect opened input: {error}"),
        )
    })?;
    let path_metadata = fs::symlink_metadata(path).map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot inspect input: {error}"),
        )
    })?;
    if !opened.file_type().is_file()
        || !path_metadata.file_type().is_file()
        || path_metadata.file_type().is_symlink()
    {
        return Err(CompileError::new(
            "invalid_input_type",
            subject,
            "input must be a regular non-symlink file",
        ));
    }
    #[cfg(unix)]
    if opened.dev() != path_metadata.dev() || opened.ino() != path_metadata.ino() {
        return Err(CompileError::new(
            "input_changed_during_admission",
            subject,
            "opened descriptor no longer names the inspected source inode",
        ));
    }
    let declared = usize::try_from(opened.len()).unwrap_or(usize::MAX);
    if declared > limit {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("input exceeds {limit} bytes"),
        ));
    }
    let take_limit = u64::try_from(limit).unwrap_or(u64::MAX).saturating_add(1);
    let mut bytes = Vec::with_capacity(declared);
    file.take(take_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| {
            CompileError::new(
                "input_unavailable",
                subject,
                format!("cannot read opened input: {error}"),
            )
        })?;
    if bytes.len() > limit {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("input grew beyond {limit} bytes while being read"),
        ));
    }
    Ok(bytes)
}

fn safe_relative_path(raw: &str) -> Result<PathBuf, CompileError> {
    if raw.is_empty() || raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "unsafe_source_path",
            "manifest",
            "source path must be nonempty and bounded",
        ));
    }
    let path = PathBuf::from(raw);
    if path.is_absolute()
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(CompileError::new(
            "unsafe_source_path",
            "manifest",
            "source paths must be relative and contain no dot or parent components",
        ));
    }
    Ok(path)
}

fn set_once<T>(
    slot: &mut Option<T>,
    value: T,
    code: &'static str,
    subject: &str,
) -> Result<(), CompileError> {
    if slot.replace(value).is_some() {
        Err(CompileError::new(
            code,
            subject,
            "record may appear exactly once",
        ))
    } else {
        Ok(())
    }
}

fn require_field_count(
    fields: &[&str],
    expected: usize,
    source_id: &str,
    line: usize,
) -> Result<(), CompileError> {
    if fields.len() == expected {
        Ok(())
    } else {
        Err(CompileError::new(
            "invalid_record_shape",
            format!("source:{source_id}:line:{line}"),
            format!(
                "expected {expected} tab-separated fields, found {}",
                fields.len()
            ),
        ))
    }
}

fn require_text(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.trim().is_empty() {
        return Err(CompileError::new(
            "missing_required_field",
            subject,
            format!("{field} must be nonblank"),
        ));
    }
    if raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("{field} exceeds {MAX_TEXT_BYTES} bytes"),
        ));
    }
    if raw.chars().any(char::is_control) {
        return Err(CompileError::new(
            "invalid_text_field",
            subject,
            format!("{field} contains a control character"),
        ));
    }
    Ok(raw.to_string())
}

fn require_identifier(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() || raw.len() > MAX_IDENTIFIER_BYTES {
        return Err(CompileError::new(
            "invalid_identifier",
            subject,
            format!("{field} must contain 1..={MAX_IDENTIFIER_BYTES} bytes"),
        ));
    }
    if !raw
        .bytes()
        .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.'))
    {
        return Err(CompileError::new(
            "invalid_identifier",
            subject,
            format!("{field} may contain only ASCII letters, digits, dot, dash, and underscore"),
        ));
    }
    Ok(raw.to_string())
}

fn push_part(output: &mut Vec<u8>, part: &[u8]) {
    output.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
    output.extend_from_slice(part);
}

fn parse_canonical_u32(raw: &str, field: &str, subject: &str) -> Result<u32, CompileError> {
    let value = raw.parse::<u32>().map_err(|_| {
        CompileError::new(
            "invalid_unsigned_integer",
            subject,
            format!("{field} must be a canonical base-10 u32"),
        )
    })?;
    if value.to_string() != raw {
        return Err(CompileError::new(
            "noncanonical_unsigned_integer",
            subject,
            format!("{field} must not contain signs, whitespace, or leading zeroes"),
        ));
    }
    Ok(value)
}

fn parse_interface_surface(
    fields: &[&str],
    source_id: &str,
    line_number: usize,
    subject: &str,
    record_hash: ContentHash,
) -> Result<RawSurfaceSpec, CompileError> {
    require_field_count(fields, 6, source_id, line_number)?;
    Ok(RawSurfaceSpec {
        chemistry: require_text(fields[1], "surface chemistry", subject)?,
        phase: require_text(fields[2], "surface phase", subject)?,
        process: require_text(fields[3], "surface process", subject)?,
        revision: parse_canonical_u32(fields[4], "material revision", subject)?,
        texture_frame: require_identifier(fields[5], "texture-frame id", subject)?,
        record_hash,
    })
}

#[allow(clippy::too_many_lines)] // Exact record grammar and per-record refusals stay co-located.
fn parse_source(
    source: &LoadedSource,
    raw: &mut RawDatabase,
    profile: SourceProfile,
) -> Result<(), CompileError> {
    if source.text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_source_encoding",
            format!("source:{}", source.spec.id),
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = source.text.lines();
    if lines.next() != Some(SOURCE_HEADER) {
        return Err(CompileError::new(
            "unsupported_source_schema",
            format!("source:{}", source.spec.id),
            format!("first line must be {SOURCE_HEADER:?}"),
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let subject = format!("source:{}:line:{line_number}", source.spec.id);
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                subject,
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_source_record",
                subject,
                "blank records are forbidden",
            ));
        }
        raw.records = raw.records.checked_add(1).ok_or_else(|| {
            CompileError::new("resource_limit", "records", "record count overflowed")
        })?;
        if raw.records > MAX_RECORDS {
            return Err(CompileError::new(
                "resource_limit",
                "records",
                format!("at most {MAX_RECORDS} source records are admitted"),
            ));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let record_hash = hash_domain(SOURCE_RECORD_DOMAIN, line.as_bytes());
        match fields.first().copied() {
            Some("surface_a") => {
                if profile != SourceProfile::Interface {
                    return Err(CompileError::new(
                        "unsupported_source_record",
                        subject,
                        "surface identity records require interface-tsv-v1",
                    ));
                }
                let surface = parse_interface_surface(
                    &fields,
                    &source.spec.id,
                    line_number,
                    &subject,
                    record_hash,
                )?;
                set_once(
                    &mut raw.surface_a,
                    surface,
                    "duplicate_interface_surface",
                    "interface:surface_a",
                )?;
            }
            Some("surface_b") => {
                if profile != SourceProfile::Interface {
                    return Err(CompileError::new(
                        "unsupported_source_record",
                        subject,
                        "surface identity records require interface-tsv-v1",
                    ));
                }
                let surface = parse_interface_surface(
                    &fields,
                    &source.spec.id,
                    line_number,
                    &subject,
                    record_hash,
                )?;
                set_once(
                    &mut raw.surface_b,
                    surface,
                    "duplicate_interface_surface",
                    "interface:surface_b",
                )?;
            }
            Some("context") => {
                if profile != SourceProfile::Interface {
                    return Err(CompileError::new(
                        "unsupported_source_record",
                        subject,
                        "system context records require interface-tsv-v1",
                    ));
                }
                require_field_count(&fields, 5, &source.spec.id, line_number)?;
                let third_body = if fields[2] == "-" {
                    None
                } else {
                    Some(require_text(fields[2], "third body", &subject)?)
                };
                let context = RawSystemContext {
                    medium: require_text(fields[1], "medium", &subject)?,
                    third_body,
                    environment: require_text(fields[3], "environment", &subject)?,
                    history: require_text(fields[4], "history", &subject)?,
                    record_hash,
                };
                set_once(
                    &mut raw.context,
                    context,
                    "duplicate_interface_context",
                    "interface:context",
                )?;
            }
            Some("observation") => {
                require_field_count(&fields, 5, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "observation id", &subject)?;
                let observation = RawObservation {
                    id: id.clone(),
                    specimen: require_text(fields[2], "specimen", &subject)?,
                    method: require_text(fields[3], "method", &subject)?,
                    caveats: bounded_optional_text(fields[4], "caveats", &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.observations.insert(id.clone(), observation).is_some() {
                    return Err(CompileError::new(
                        "duplicate_observation",
                        format!("observation:{id}"),
                        "observation ids are globally unique",
                    ));
                }
            }
            Some("scalar") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "claim id", &subject)?;
                let claim = RawClaim {
                    id: id.clone(),
                    observations: parse_identifier_list(
                        fields[2],
                        "observation references",
                        &subject,
                    )?,
                    property: require_identifier(fields[3], "property", &subject)?,
                    value: RawClaimValue::Scalar {
                        number: require_number_token(fields[4], "scalar value", &subject)?,
                        unit: require_unit(fields[5], "scalar unit", &subject)?,
                    },
                    interpolation: parse_interpolation(fields[6], &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                insert_claim(raw, id, claim)?;
            }
            Some("curve") => {
                require_field_count(&fields, 9, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "claim id", &subject)?;
                let points = parse_points(fields[7], &subject)?;
                raw.curve_knots = raw.curve_knots.checked_add(points.len()).ok_or_else(|| {
                    CompileError::new(
                        "resource_limit",
                        "curve_knots",
                        "aggregate knot count overflowed",
                    )
                })?;
                if raw.curve_knots > MAX_CURVE_KNOTS {
                    return Err(CompileError::new(
                        "resource_limit",
                        "curve_knots",
                        format!("at most {MAX_CURVE_KNOTS} curve knots are admitted"),
                    ));
                }
                let claim = RawClaim {
                    id: id.clone(),
                    observations: parse_identifier_list(
                        fields[2],
                        "observation references",
                        &subject,
                    )?,
                    property: require_identifier(fields[3], "property", &subject)?,
                    value: RawClaimValue::Curve {
                        abscissa: require_identifier(fields[4], "abscissa axis", &subject)?,
                        abscissa_unit: require_unit(fields[5], "abscissa unit", &subject)?,
                        ordinate_unit: require_unit(fields[6], "ordinate unit", &subject)?,
                        points,
                    },
                    interpolation: parse_interpolation(fields[8], &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                insert_claim(raw, id, claim)?;
            }
            Some("uncertainty") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let kind = match fields[2] {
                    "unstated" => {
                        if fields[3..] != ["-", "-", "-", "-"] {
                            return Err(CompileError::new(
                                "invalid_uncertainty",
                                format!("claim:{claim}:uncertainty"),
                                "unstated uncertainty requires four '-' placeholders",
                            ));
                        }
                        RawUncertaintyKind::Unstated
                    }
                    "absolute" => RawUncertaintyKind::Absolute {
                        number: require_number_token(fields[3], "absolute uncertainty", &subject)?,
                        unit: require_unit(fields[4], "uncertainty unit", &subject)?,
                        confidence: require_number_token(fields[5], "confidence", &subject)?,
                        confidence_unit: require_unit(fields[6], "confidence unit", &subject)?,
                    },
                    "relative" => RawUncertaintyKind::Relative {
                        number: require_number_token(fields[3], "relative uncertainty", &subject)?,
                        unit: require_unit(fields[4], "uncertainty unit", &subject)?,
                        confidence: require_number_token(fields[5], "confidence", &subject)?,
                        confidence_unit: require_unit(fields[6], "confidence unit", &subject)?,
                    },
                    _ => {
                        return Err(CompileError::new(
                            "invalid_uncertainty",
                            format!("claim:{claim}:uncertainty"),
                            "kind must be unstated, absolute, or relative",
                        ));
                    }
                };
                if raw
                    .uncertainties
                    .insert(
                        claim.clone(),
                        RawUncertainty {
                            kind,
                            source_hash: record_hash,
                        },
                    )
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_uncertainty",
                        format!("claim:{claim}:uncertainty"),
                        "each claim requires exactly one uncertainty record",
                    ));
                }
            }
            Some("validity") => {
                require_field_count(&fields, 6, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let axis = require_identifier(fields[2], "validity axis", &subject)?;
                let validity = RawValidity {
                    axis: axis.clone(),
                    lower: require_number_token(fields[3], "validity lower bound", &subject)?,
                    upper: require_number_token(fields[4], "validity upper bound", &subject)?,
                    unit: require_unit(fields[5], "validity unit", &subject)?,
                    source_hash: record_hash,
                };
                if raw
                    .validities
                    .entry(claim.clone())
                    .or_default()
                    .insert(axis.clone(), validity)
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_validity_axis",
                        format!("claim:{claim}:validity:{axis}"),
                        "one claim may name each validity axis once",
                    ));
                }
            }
            Some("frame") => {
                require_field_count(&fields, 4, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let frame = RawFrame {
                    source: require_identifier(fields[2], "source frame", &subject)?,
                    target: require_identifier(fields[3], "target frame", &subject)?,
                    source_hash: record_hash,
                };
                if raw.frames.insert(claim.clone(), frame).is_some() {
                    return Err(CompileError::new(
                        "duplicate_frame",
                        format!("claim:{claim}:frame"),
                        "one claim may declare at most one paired frame mapping",
                    ));
                }
            }
            Some("joint") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let observation = require_identifier(fields[1], "observation id", &subject)?;
                let block_id = require_identifier(fields[2], "joint block id", &subject)?;
                let members = parse_member_list(fields[3], &subject)?;
                if members.len() > MAX_JOINT_MEMBERS {
                    return Err(CompileError::new(
                        "resource_limit",
                        format!("joint:{observation}:{block_id}"),
                        format!("at most {MAX_JOINT_MEMBERS} members are admitted per block"),
                    ));
                }
                let covariance = parse_numeric_list(fields[4], "covariance", &subject)?;
                let correlation = if fields[5] == "-" {
                    if fields[6] != "-" {
                        return Err(CompileError::new(
                            "invalid_correlation_unit",
                            format!("joint:{observation}:{block_id}"),
                            "absent correlation requires a '-' unit placeholder",
                        ));
                    }
                    None
                } else {
                    if fields[6] != "1" {
                        return Err(CompileError::new(
                            "invalid_correlation_unit",
                            format!("joint:{observation}:{block_id}"),
                            "correlation must explicitly declare the dimensionless unit '1'",
                        ));
                    }
                    Some(parse_numeric_list(fields[5], "correlation", &subject)?)
                };
                let joint = RawJoint {
                    observation: observation.clone(),
                    block_id: block_id.clone(),
                    members,
                    covariance,
                    correlation,
                    source_hash: record_hash,
                };
                if raw
                    .joints
                    .insert((observation.clone(), block_id.clone()), joint)
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_joint_block",
                        format!("joint:{observation}:{block_id}"),
                        "joint block ids are unique within an observation",
                    ));
                }
            }
            Some("species" | "nasa9" | "kinetics") => {
                return Err(CompileError::new(
                    "unsupported_source_record",
                    subject,
                    "material-tsv-v1 admits only material records; species metadata must use species-v1, NASA-9 data must use nasa9-v1, and kinetics data must use kinetics-v1",
                ));
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_source_record",
                    subject,
                    format!("unknown record type {other:?}"),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    raw.decisions.push(Decision::admit(
        format!("source:{}", source.spec.id),
        "source_schema_admitted",
        match profile {
            SourceProfile::Material => "bounded material-tsv-v1 source parsed without inference",
            SourceProfile::Interface => {
                "bounded interface-tsv-v1 source parsed without inferred system identity"
            }
            SourceProfile::Nasa9 | SourceProfile::Kinetics | SourceProfile::Species => {
                unreachable!("specialized profiles use dedicated parsers")
            }
        },
        Some(source.hash),
    ));
    Ok(())
}

fn parse_species_id(raw: &str, subject: &str) -> Result<SpeciesId, CompileError> {
    SpeciesId::new(raw.to_string()).map_err(|error| {
        CompileError::new(
            "invalid_species_id",
            subject,
            format!("fs-qty refused the canonical species id: {error}"),
        )
    })
}

fn parse_reaction_id(raw: &str, subject: &str) -> Result<ReactionId, CompileError> {
    ReactionId::new(raw.to_string()).map_err(|error| {
        CompileError::new(
            "invalid_reaction_id",
            subject,
            format!("fs-qty refused the canonical reaction id: {error}"),
        )
    })
}

fn require_elemental_reference(raw: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() || raw.len() > MAX_IDENTIFIER_BYTES {
        return Err(CompileError::new(
            "invalid_elemental_reference",
            subject,
            format!("elemental-reference id must contain 1..={MAX_IDENTIFIER_BYTES} bytes"),
        ));
    }
    let bytes = raw.as_bytes();
    if !bytes[0].is_ascii_alphanumeric()
        || !bytes.iter().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(CompileError::new(
            "invalid_elemental_reference",
            subject,
            "elemental-reference id must use compact ASCII without whitespace",
        ));
    }
    Ok(raw.to_string())
}

fn parse_species_source(
    source: &LoadedSource,
    raw: &mut RawSpeciesDatabase,
) -> Result<(), CompileError> {
    if source.text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_source_encoding",
            format!("source:{}", source.spec.id),
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = source.text.lines();
    if lines.next() != Some(SPECIES_SOURCE_HEADER) {
        return Err(CompileError::new(
            "unsupported_source_schema",
            format!("source:{}", source.spec.id),
            format!("first line must be {SPECIES_SOURCE_HEADER:?}"),
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let subject = format!("source:{}:line:{line_number}", source.spec.id);
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                subject,
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_source_record",
                subject,
                "blank records are forbidden",
            ));
        }
        raw.records = raw.records.checked_add(1).ok_or_else(|| {
            CompileError::new("resource_limit", "records", "record count overflowed")
        })?;
        if raw.records > 1 {
            return Err(CompileError::new(
                "duplicate_species_association",
                subject,
                "species-v1 admits exactly one complete species association record",
            ));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.first().copied() != Some("species") {
            return Err(CompileError::new(
                "unknown_source_record",
                subject,
                "species-v1 expects exactly one 'species' record",
            ));
        }
        require_field_count(&fields, 9, &source.spec.id, line_number)?;
        let association = RawSpeciesAssociation {
            species: parse_species_id(fields[1], &subject)?,
            molar_mass: require_number_token(fields[2], "molar mass", &subject)?,
            molar_mass_unit: require_unit(fields[3], "molar-mass unit", &subject)?,
            standard_state_phase: require_identifier(fields[4], "standard-state phase", &subject)?,
            reference_eos: require_identifier(fields[5], "reference EOS", &subject)?,
            reference_pressure: require_number_token(fields[6], "reference pressure", &subject)?,
            reference_pressure_unit: require_unit(fields[7], "reference-pressure unit", &subject)?,
            elemental_reference: require_elemental_reference(fields[8], &subject)?,
            source_id: source.spec.id.clone(),
            source_hash: source.hash,
            record_hash: hash_domain(SOURCE_RECORD_DOMAIN, line.as_bytes()),
        };
        if raw.association.replace(association).is_some() {
            return Err(CompileError::new(
                "duplicate_species_association",
                subject,
                "species-v1 admits exactly one complete species association record",
            ));
        }
    }
    raw.decisions.push(Decision::admit(
        format!("source:{}", source.spec.id),
        "source_schema_admitted",
        "bounded species-v1 source parsed with explicit units and thermochemical convention fields",
        Some(source.hash),
    ));
    Ok(())
}

fn parse_nasa9_source(
    source: &LoadedSource,
    raw: &mut RawNasa9Database,
) -> Result<(), CompileError> {
    if source.text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_source_encoding",
            format!("source:{}", source.spec.id),
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = source.text.lines();
    if lines.next() != Some(NASA9_SOURCE_HEADER) {
        return Err(CompileError::new(
            "unsupported_source_schema",
            format!("source:{}", source.spec.id),
            format!("first line must be {NASA9_SOURCE_HEADER:?}"),
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let subject = format!("source:{}:line:{line_number}", source.spec.id);
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                subject,
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_source_record",
                subject,
                "blank records are forbidden",
            ));
        }
        raw.records = raw.records.checked_add(1).ok_or_else(|| {
            CompileError::new("resource_limit", "records", "record count overflowed")
        })?;
        if raw.records > MAX_RECORDS {
            return Err(CompileError::new(
                "resource_limit",
                "records",
                format!("at most {MAX_RECORDS} source records are admitted"),
            ));
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let record_hash = hash_domain(SOURCE_RECORD_DOMAIN, line.as_bytes());
        match fields.first().copied() {
            Some("region") => {
                require_field_count(&fields, 8, &source.spec.id, line_number)?;
                let species = parse_species_id(fields[1], &subject)?;
                let id = require_identifier(fields[2], "region id", &subject)?;
                let key = (species.clone(), id.clone());
                if !raw.regions.contains_key(&key) && raw.regions.len() == MAX_NASA9_REGIONS {
                    return Err(CompileError::new(
                        "resource_limit",
                        subject,
                        format!(
                            "one species artifact admits at most {MAX_NASA9_REGIONS} NASA-9 regions"
                        ),
                    ));
                }
                let region = RawNasa9Region {
                    species,
                    id,
                    temperature_lower: require_number_token(
                        fields[3],
                        "temperature lower bound",
                        &subject,
                    )?,
                    temperature_upper: require_number_token(
                        fields[4],
                        "temperature upper bound",
                        &subject,
                    )?,
                    temperature_unit: require_unit(fields[5], "temperature unit", &subject)?,
                    reference_pressure: require_number_token(
                        fields[6],
                        "reference pressure",
                        &subject,
                    )?,
                    reference_pressure_unit: require_unit(
                        fields[7],
                        "reference pressure unit",
                        &subject,
                    )?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.regions.insert(key, region).is_some() {
                    return Err(CompileError::new(
                        "duplicate_nasa9_region",
                        subject,
                        "a species may name each NASA-9 region id once",
                    ));
                }
            }
            Some("coefficient") => {
                require_field_count(&fields, 6, &source.spec.id, line_number)?;
                let species = parse_species_id(fields[1], &subject)?;
                let region = require_identifier(fields[2], "region id", &subject)?;
                let name = require_identifier(fields[3], "coefficient name", &subject)?;
                if !NASA9_COEFFICIENT_NAMES.contains(&name.as_str()) {
                    return Err(CompileError::new(
                        "unexpected_nasa9_coefficient",
                        subject,
                        format!(
                            "coefficient {name:?} is not one of the exact a0..a8 NASA-9 fields"
                        ),
                    ));
                }
                let key = (species, region, name);
                if !raw.coefficients.contains_key(&key)
                    && raw.coefficients.len() == MAX_NASA9_COEFFICIENTS
                {
                    return Err(CompileError::new(
                        "resource_limit",
                        subject,
                        format!(
                            "one species artifact admits at most {MAX_NASA9_COEFFICIENTS} NASA-9 coefficient records"
                        ),
                    ));
                }
                let coefficient = RawNasa9Coefficient {
                    value: require_number_token(fields[4], "coefficient value", &subject)?,
                    unit: require_unit(fields[5], "coefficient unit", &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.coefficients.insert(key, coefficient).is_some() {
                    return Err(CompileError::new(
                        "duplicate_nasa9_coefficient",
                        subject,
                        "a NASA-9 region may declare each coefficient exactly once",
                    ));
                }
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_source_record",
                    subject,
                    format!("unknown NASA-9 record type {other:?}; expected region or coefficient"),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    raw.decisions.push(Decision::admit(
        format!("source:{}", source.spec.id),
        "source_schema_admitted",
        "bounded nasa9-v1 source parsed without unit or species inference",
        Some(source.hash),
    ));
    Ok(())
}

fn parse_kinetics_source(
    source: &LoadedSource,
    raw: &mut RawKineticsDatabase,
) -> Result<(), CompileError> {
    if source.text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_source_encoding",
            format!("source:{}", source.spec.id),
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = source.text.lines();
    if lines.next() != Some(KINETICS_SOURCE_HEADER) {
        return Err(CompileError::new(
            "unsupported_source_schema",
            format!("source:{}", source.spec.id),
            format!("first line must be {KINETICS_SOURCE_HEADER:?}"),
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let subject = format!("source:{}:line:{line_number}", source.spec.id);
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                subject,
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_source_record",
                subject,
                "blank records are forbidden",
            ));
        }
        raw.records = raw.records.checked_add(1).ok_or_else(|| {
            CompileError::new("resource_limit", "records", "record count overflowed")
        })?;
        if raw.records > MAX_RECORDS {
            return Err(CompileError::new(
                "resource_limit",
                "records",
                format!("at most {MAX_RECORDS} source records are admitted"),
            ));
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let record_hash = hash_domain(SOURCE_RECORD_DOMAIN, line.as_bytes());
        match fields.first().copied() {
            Some("reaction") => {
                require_field_count(&fields, 6, &source.spec.id, line_number)?;
                let reaction = parse_reaction_id(fields[1], &subject)?;
                if fields[2] != "first-order" {
                    return Err(CompileError::new(
                        "unsupported_kinetics_rate_basis",
                        subject,
                        "kinetics-v1 admits only the explicit first-order rate basis",
                    ));
                }
                if !raw.reactions.contains_key(&reaction)
                    && raw.reactions.len() == MAX_KINETICS_REACTIONS
                {
                    return Err(CompileError::new(
                        "resource_limit",
                        subject,
                        "one kinetics-v1 artifact admits exactly one reaction",
                    ));
                }
                let entry = RawKineticsReaction {
                    reaction: reaction.clone(),
                    temperature_lower: require_number_token(
                        fields[3],
                        "temperature lower bound",
                        &subject,
                    )?,
                    temperature_upper: require_number_token(
                        fields[4],
                        "temperature upper bound",
                        &subject,
                    )?,
                    temperature_unit: require_unit(fields[5], "temperature unit", &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.reactions.insert(reaction, entry).is_some() {
                    return Err(CompileError::new(
                        "duplicate_kinetics_reaction",
                        subject,
                        "a kinetics-v1 artifact may declare its reaction exactly once",
                    ));
                }
            }
            Some("parameter") => {
                require_field_count(&fields, 5, &source.spec.id, line_number)?;
                let reaction = parse_reaction_id(fields[1], &subject)?;
                let name = require_identifier(fields[2], "kinetics parameter name", &subject)?;
                if !KINETICS_PARAMETER_NAMES.contains(&name.as_str()) {
                    return Err(CompileError::new(
                        "unexpected_kinetics_parameter",
                        subject,
                        format!(
                            "parameter {name:?} is not one of the exact activation_temperature and pre_exponential fields"
                        ),
                    ));
                }
                let key = (reaction, name);
                if !raw.parameters.contains_key(&key)
                    && raw.parameters.len() == MAX_KINETICS_PARAMETERS
                {
                    return Err(CompileError::new(
                        "resource_limit",
                        subject,
                        format!(
                            "one kinetics-v1 artifact admits exactly {MAX_KINETICS_PARAMETERS} parameter records"
                        ),
                    ));
                }
                let parameter = RawKineticsParameter {
                    value: require_number_token(fields[3], "kinetics parameter value", &subject)?,
                    unit: require_unit(fields[4], "kinetics parameter unit", &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.parameters.insert(key, parameter).is_some() {
                    return Err(CompileError::new(
                        "duplicate_kinetics_parameter",
                        subject,
                        "a reaction may declare each kinetics parameter exactly once",
                    ));
                }
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_source_record",
                    subject,
                    format!(
                        "unknown kinetics record type {other:?}; expected reaction or parameter"
                    ),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    raw.decisions.push(Decision::admit(
        format!("source:{}", source.spec.id),
        "source_schema_admitted",
        "bounded kinetics-v1 source parsed with an explicit first-order basis and no inferred units",
        Some(source.hash),
    ));
    Ok(())
}

fn insert_claim(raw: &mut RawDatabase, id: String, claim: RawClaim) -> Result<(), CompileError> {
    if raw.claims.insert(id.clone(), claim).is_some() {
        Err(CompileError::new(
            "duplicate_claim",
            format!("claim:{id}"),
            "claim ids are globally unique",
        ))
    } else {
        Ok(())
    }
}

fn bounded_optional_text(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("{field} exceeds {MAX_TEXT_BYTES} bytes"),
        ));
    }
    if raw.chars().any(char::is_control) {
        return Err(CompileError::new(
            "invalid_text_field",
            subject,
            format!("{field} contains a control character"),
        ));
    }
    Ok(raw.to_string())
}

fn require_number_token(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() || raw.trim() != raw || raw.len() > 128 {
        return Err(CompileError::new(
            "invalid_number",
            subject,
            format!("{field} must be a bounded whitespace-free numeric token"),
        ));
    }
    parse_source_f64(raw, field, subject)?;
    Ok(raw.to_string())
}

fn require_unit(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "ambiguous_unit",
            subject,
            format!("{field} is required; dimensionless values must explicitly use '1'"),
        ));
    }
    if raw.trim() != raw || raw.len() > 256 || raw.chars().any(char::is_whitespace) {
        return Err(CompileError::new(
            "invalid_unit",
            subject,
            format!("{field} must be a bounded whitespace-free fs-qty expression"),
        ));
    }
    Ok(raw.to_string())
}

fn parse_identifier_list(
    raw: &str,
    field: &str,
    subject: &str,
) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "missing_observation_reference",
            subject,
            format!("{field} must name at least one observation"),
        ));
    }
    let mut values = Vec::new();
    for item in raw.split(',') {
        values.push(require_identifier(item, field, subject)?);
    }
    values.sort();
    if values.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(CompileError::new(
            "duplicate_reference",
            subject,
            format!("{field} contains a duplicate id"),
        ));
    }
    Ok(values)
}

fn parse_points(raw: &str, subject: &str) -> Result<Vec<(String, String)>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "malformed_curve",
            subject,
            "curve point list is empty",
        ));
    }
    let mut points = Vec::new();
    for point in raw.split(',') {
        let (x, y) = point.split_once(':').ok_or_else(|| {
            CompileError::new(
                "malformed_curve",
                subject,
                "each curve point must be encoded as x:y",
            )
        })?;
        if y.contains(':') {
            return Err(CompileError::new(
                "malformed_curve",
                subject,
                "curve points contain too many ':' separators",
            ));
        }
        points.push((
            require_number_token(x, "curve abscissa", subject)?,
            require_number_token(y, "curve ordinate", subject)?,
        ));
    }
    if points.len() < 2 {
        return Err(CompileError::new(
            "malformed_curve",
            subject,
            "a curve requires at least two points",
        ));
    }
    Ok(points)
}

fn parse_member_list(raw: &str, subject: &str) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "invalid_joint_members",
            subject,
            "joint block must name at least one member",
        ));
    }
    let members: Vec<String> = raw.split(',').map(str::to_string).collect();
    if members.iter().any(String::is_empty) {
        return Err(CompileError::new(
            "invalid_joint_members",
            subject,
            "joint member references may not be blank",
        ));
    }
    Ok(members)
}

fn parse_numeric_list(raw: &str, field: &str, subject: &str) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "invalid_matrix",
            subject,
            format!("{field} matrix is empty"),
        ));
    }
    raw.split(',')
        .map(|item| require_number_token(item, field, subject))
        .collect()
}

fn parse_interpolation(raw: &str, subject: &str) -> Result<InterpolationPolicy, CompileError> {
    match raw {
        "linear" => Ok(InterpolationPolicy::LinearInside),
        "constant" => Ok(InterpolationPolicy::ConstantWithinValidity),
        "tabulated" => Ok(InterpolationPolicy::TabulatedOnly),
        _ => Err(CompileError::new(
            "invalid_interpolation",
            subject,
            "interpolation must be linear, constant, or tabulated",
        )),
    }
}

fn parse_source_f64(raw: &str, field: &str, subject: &str) -> Result<f64, CompileError> {
    let value = raw.parse::<f64>().map_err(|_| {
        CompileError::new(
            "invalid_number",
            subject,
            format!("{field} is not an exact finite f64 token"),
        )
    })?;
    if !value.is_finite() {
        return Err(CompileError::new(
            "non_finite_number",
            subject,
            format!("{field} must be finite"),
        ));
    }
    let exponent = raw.find(|character| matches!(character, 'e' | 'E'));
    let significand = exponent.map_or(raw, |index| &raw[..index]);
    if value == 0.0 && significand.bytes().any(|byte| matches!(byte, b'1'..=b'9')) {
        return Err(CompileError::new(
            "numeric_underflow",
            subject,
            format!("{field} has a nonzero significand that underflows f64"),
        ));
    }
    if value.to_bits() == (-0.0f64).to_bits() {
        return Err(CompileError::new(
            "negative_zero",
            subject,
            format!("{field} uses forbidden signed zero"),
        ));
    }
    Ok(value)
}

fn parse_confidence(raw: &str, unit: &str, subject: &str) -> Result<f64, CompileError> {
    if unit != "1" {
        return Err(CompileError::new(
            "invalid_confidence_unit",
            subject,
            "confidence must use the exact dimensionless unit '1' because schema v1 has no \
             separate confidence normalization receipt",
        ));
    }
    let quantity = parse_linear_quantity(raw, unit, subject)?;
    if quantity.dims != Dims::NONE {
        return Err(CompileError::new(
            "invalid_confidence_unit",
            subject,
            "confidence must explicitly use a dimensionless unit",
        ));
    }
    let value = quantity.value;
    if value > 0.0 && value < 1.0 {
        Ok(value)
    } else {
        Err(CompileError::new(
            "invalid_confidence",
            subject,
            "confidence must be strictly between zero and one",
        ))
    }
}

fn unit_transform(unit: &str, subject: &str) -> Result<(Dims, f64, f64), CompileError> {
    if unit == "1" {
        return Ok((Dims::NONE, 1.0, 0.0));
    }
    if !unit
        .chars()
        .next()
        .is_some_and(|first| first.is_alphabetic() || first == '%' || first == 'µ')
    {
        return Err(CompileError::new(
            "invalid_unit",
            subject,
            "unit must begin with a unit symbol; numeric continuations are forbidden",
        ));
    }
    let zero = parse_qty_with_budget(&format!("0 {unit}"), QTY_BUDGET).map_err(|_| {
        CompileError::new(
            "invalid_unit",
            subject,
            "fs-qty refused the explicit unit expression",
        )
    })?;
    let one = parse_qty_with_budget(&format!("1 {unit}"), QTY_BUDGET).map_err(|_| {
        CompileError::new(
            "invalid_unit",
            subject,
            "fs-qty refused the explicit unit expression",
        )
    })?;
    if zero.dims != one.dims {
        return Err(CompileError::new(
            "unstable_unit_dimensions",
            subject,
            "unit dimensions differ between zero and one probes",
        ));
    }
    let scale = one.value - zero.value;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(CompileError::new(
            "invalid_unit_scale",
            subject,
            "unit transform must have a finite positive scale",
        ));
    }
    let offset = positive_zero(zero.value);
    Ok((zero.dims, scale, offset))
}

fn parse_quantity(number: &str, unit: &str, subject: &str) -> Result<ParsedQuantity, CompileError> {
    let source_value = parse_source_f64(number, "quantity", subject)?;
    let (dims, scale, offset) = unit_transform(unit, subject)?;
    let parsed_value = if unit == "1" {
        source_value
    } else {
        let parsed =
            parse_qty_with_budget(&format!("{number} {unit}"), QTY_BUDGET).map_err(|_| {
                CompileError::new(
                    "invalid_quantity",
                    subject,
                    "fs-qty refused the explicit number and unit",
                )
            })?;
        if parsed.dims != dims {
            return Err(CompileError::new(
                "unstable_unit_dimensions",
                subject,
                "literal dimensions differ from the unit probes",
            ));
        }
        parsed.value
    };
    if !parsed_value.is_finite() {
        return Err(CompileError::new(
            "non_finite_normalization",
            subject,
            "normalized value is not finite",
        ));
    }
    let receipt_value = source_value * scale + offset;
    let tolerance = f64::EPSILON * 8.0 * parsed_value.abs().max(receipt_value.abs()).max(1.0);
    if !receipt_value.is_finite() || (receipt_value - parsed_value).abs() > tolerance {
        return Err(CompileError::new(
            "non_reproducible_unit_transform",
            subject,
            "derived affine receipt disagrees materially with the independently parsed quantity",
        ));
    }
    if offset == 0.0 && source_value != 0.0 && receipt_value == 0.0 {
        return Err(CompileError::new(
            "normalization_underflow",
            subject,
            "nonzero source quantity underflows the normalized representation",
        ));
    }
    Ok(ParsedQuantity {
        value: positive_zero(receipt_value),
        dims,
        scale,
        offset,
        unit: unit.to_string(),
        literal_hash: source_literal_hash(number, unit),
    })
}

fn parse_linear_quantity(
    number: &str,
    unit: &str,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    let source_value = parse_source_f64(number, "quantity", subject)?;
    let (dims, scale, _) = unit_transform(unit, subject)?;
    let value = source_value * scale;
    if !value.is_finite() {
        return Err(CompileError::new(
            "non_finite_normalization",
            subject,
            "normalized linear magnitude is not finite",
        ));
    }
    if source_value != 0.0 && value == 0.0 {
        return Err(CompileError::new(
            "normalization_underflow",
            subject,
            "nonzero source magnitude underflows the normalized representation",
        ));
    }
    Ok(ParsedQuantity {
        value: positive_zero(value),
        dims,
        scale,
        offset: 0.0,
        unit: unit.to_string(),
        literal_hash: source_literal_hash(number, unit),
    })
}

fn parse_declared_quantity(
    number: &str,
    unit: &str,
    kind: Option<QuantityKind>,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    if kind == Some(QuantityKind::TemperatureDifference) {
        return parse_linear_quantity(number, unit, subject);
    }
    let mut quantity = parse_quantity(number, unit, subject)?;
    if let Some(QuantityKind::Frequency(target)) = kind {
        // Frequency meaning comes from this closed source-unit convention,
        // never from inverse-time dimensions alone.
        let source = match unit {
            "Hz" | "kHz" | "MHz" | "GHz" | "1/s" | "s^-1" => FrequencyConvention::Cyclic,
            "rad/s" => FrequencyConvention::Angular,
            _ => {
                return Err(CompileError::new(
                    "frequency_unit_convention",
                    subject,
                    "frequency values require Hz, kHz, MHz, GHz, 1/s, s^-1, or rad/s",
                ));
            }
        };
        let scale = fs_qty::semantic::SemanticQty::new(
            fs_qty::QtyAny::new(quantity.scale, quantity.dims),
            SemanticType::new(QuantityKind::Frequency(source), ValueForm::Static),
        )
        .map_err(|error| {
            CompileError::new("frequency_conversion_refused", subject, error.to_string())
        })?;
        quantity.scale = fs_qty::semantic::convert_frequency(scale, target)
            .map_err(|error| {
                CompileError::new("frequency_conversion_refused", subject, error.to_string())
            })?
            .quantity()
            .value;
        let source_value = parse_source_f64(number, "frequency", subject)?;
        quantity.value = positive_zero(source_value * quantity.scale);
        if !quantity.value.is_finite() || (source_value != 0.0 && quantity.value == 0.0) {
            return Err(CompileError::new(
                "frequency_conversion_refused",
                subject,
                "frequency normalization overflowed or underflowed",
            ));
        }
    }
    Ok(quantity)
}

fn source_literal_hash(number: &str, unit: &str) -> ContentHash {
    let mut bytes = Vec::new();
    push_part(&mut bytes, number.as_bytes());
    push_part(&mut bytes, unit.as_bytes());
    hash_domain(SOURCE_LITERAL_DOMAIN, &bytes)
}

fn positive_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

const NASA9_TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const NASA9_PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn nasa9_coefficient_dims(name: &str) -> Dims {
    match name {
        "a0" => Dims([0, 0, 0, 2, 0, 0]),
        "a1" => Dims([0, 0, 0, 1, 0, 0]),
        "a2" | "a8" => Dims::NONE,
        "a3" => Dims([0, 0, 0, -1, 0, 0]),
        "a4" => Dims([0, 0, 0, -2, 0, 0]),
        "a5" => Dims([0, 0, 0, -3, 0, 0]),
        "a6" => Dims([0, 0, 0, -4, 0, 0]),
        "a7" => Dims([0, 0, 0, 1, 0, 0]),
        _ => unreachable!("parser admits only a0..a8"),
    }
}

fn require_nasa9_dims(
    field: &str,
    expected: Dims,
    found: Dims,
    subject: &str,
) -> Result<(), CompileError> {
    if found == expected {
        Ok(())
    } else {
        Err(CompileError::new(
            "nasa9_dims_mismatch",
            subject,
            format!("{field} has dimensions {found:?}, expected {expected:?}"),
        ))
    }
}

fn parse_nasa9_linear_quantity(
    number: &str,
    unit: &str,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    let (_, _, offset) = unit_transform(unit, subject)?;
    if offset != 0.0 {
        return Err(CompileError::new(
            "affine_nasa9_parameter_unit",
            subject,
            "NASA-9 coefficients and reference pressure require linear units; affine units are valid only for temperature-region endpoints",
        ));
    }
    parse_linear_quantity(number, unit, subject)
}

fn model_normalization(
    target: ModelNormalizationTarget,
    quantity: &ParsedQuantity,
) -> ModelNormalizationReceipt {
    ModelNormalizationReceipt::new(
        target,
        quantity.literal_hash,
        quantity.dims,
        quantity.scale,
        quantity.offset,
        quantity.unit.clone(),
        MODEL_PACK_TARGET_BASIS,
        None,
        None,
    )
}

fn parse_species_linear_quantity(
    number: &str,
    unit: &str,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    let (_, _, offset) = unit_transform(unit, subject)?;
    if offset != 0.0 {
        return Err(CompileError::new(
            "affine_species_unit",
            subject,
            "molar mass and reference pressure require linear units",
        ));
    }
    parse_linear_quantity(number, unit, subject)
}

fn species_normalization(
    target: SpeciesNormalizationTarget,
    quantity: &ParsedQuantity,
) -> SpeciesNormalizationReceipt {
    SpeciesNormalizationReceipt::new(
        target,
        quantity.literal_hash,
        quantity.dims,
        quantity.scale,
        quantity.offset,
        quantity.unit.clone(),
        SPECIES_PACK_TARGET_BASIS,
    )
}

#[allow(clippy::too_many_lines)] // Ordered source-to-runtime admission is intentionally auditable.
fn compile_species_manifest(
    manifest: Manifest,
    sources: &[LoadedSource],
    source_artifact: ContentHash,
) -> Result<CompileOutput, CompileError> {
    if sources.len() != 1 {
        return Err(CompileError::new(
            "species_source_count",
            "sources",
            "species-v1 requires exactly one source file containing one complete association",
        )
        .with_input_hash(source_artifact));
    }
    let mut raw = RawSpeciesDatabase::default();
    for source in sources {
        if let Err(error) = parse_species_source(source, &mut raw) {
            return Err(error
                .with_input_hash(source_artifact)
                .with_prior_decisions(std::mem::take(&mut raw.decisions)));
        }
    }
    let mut decisions = std::mem::take(&mut raw.decisions);
    let result = (|| -> Result<CompileOutput, CompileError> {
        let association = raw.association.ok_or_else(|| {
            CompileError::new(
                "missing_species_association",
                "sources",
                "species-v1 requires exactly one complete species record",
            )
        })?;
        if manifest.pack_id != association.species.as_str() {
            return Err(CompileError::new(
                "species_pack_id_mismatch",
                "manifest:pack_id",
                format!(
                    "species-v1 requires pack_id {:?} to equal canonical SpeciesId {:?}",
                    manifest.pack_id,
                    association.species.as_str()
                ),
            ));
        }
        if association.standard_state_phase != "gas" {
            return Err(CompileError::new(
                "unsupported_species_phase",
                format!("species:{}", association.species),
                "species-v1 admits exactly the explicit gas standard-state phase",
            ));
        }
        if association.reference_eos != "ideal-gas" {
            return Err(CompileError::new(
                "unsupported_species_reference_eos",
                format!("species:{}", association.species),
                "species-v1 admits exactly the explicit ideal-gas reference EOS",
            ));
        }

        let subject = format!("species:{}", association.species);
        let molar_mass = parse_species_linear_quantity(
            &association.molar_mass,
            &association.molar_mass_unit,
            &format!("{subject}:molar_mass"),
        )?;
        if molar_mass.dims != SPECIES_MOLAR_MASS_DIMS {
            return Err(CompileError::new(
                "species_molar_mass_dims_mismatch",
                &subject,
                format!(
                    "molar mass has dimensions {:?}, expected {:?}",
                    molar_mass.dims, SPECIES_MOLAR_MASS_DIMS
                ),
            ));
        }
        if molar_mass.value <= 0.0 {
            return Err(CompileError::new(
                "invalid_species_molar_mass",
                &subject,
                "normalized molar mass must be positive",
            ));
        }
        let reference_pressure = parse_species_linear_quantity(
            &association.reference_pressure,
            &association.reference_pressure_unit,
            &format!("{subject}:reference_pressure"),
        )?;
        if reference_pressure.dims != SPECIES_REFERENCE_PRESSURE_DIMS {
            return Err(CompileError::new(
                "species_reference_pressure_dims_mismatch",
                &subject,
                format!(
                    "reference pressure has dimensions {:?}, expected {:?}",
                    reference_pressure.dims, SPECIES_REFERENCE_PRESSURE_DIMS
                ),
            ));
        }
        if reference_pressure.value <= 0.0 {
            return Err(CompileError::new(
                "invalid_species_reference_pressure",
                &subject,
                "normalized reference pressure must be positive",
            ));
        }

        decisions.push(Decision::admit(
            "manifest:redistribution",
            "redistribution_permitted",
            "explicit permitted policy and nonblank retained terms admitted before compilation",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            "manifest:provenance",
            "licensed_citation_admitted",
            "nonblank source citation and license attached to the emitted species association",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            &subject,
            "species_pack_id_bound",
            "manifest pack_id exactly matches the fs-qty-validated SpeciesId; this is a source association, not chemical authentication",
            Some(association.record_hash),
        ));

        let source_id = association.source_id.clone();
        let source_hash = association.source_hash;
        let record_hash = association.record_hash;
        let runtime_association = SpeciesAssociation::new(
            association.species,
            molar_mass.value,
            association.standard_state_phase,
            association.reference_eos,
            reference_pressure.value,
            association.elemental_reference,
            vec![source_hash, record_hash],
            Provenance {
                source: format!(
                    "{} [source:{}] [species-association]",
                    manifest.citation, source_id
                ),
                license: manifest.license.clone(),
                artifact: Some(source_hash),
            },
        )
        .map_err(|error| {
            CompileError::new(
                "species_association_refused",
                &subject,
                format!("fs-matdb refused the source association: {error}"),
            )
        })?;
        let normalizations = vec![
            species_normalization(SpeciesNormalizationTarget::MolarMass, &molar_mass),
            species_normalization(
                SpeciesNormalizationTarget::ReferencePressure,
                &reference_pressure,
            ),
        ];
        let pack = NormalizedSpeciesPack::new(
            manifest.pack_id,
            SPECIES_COMPILER_ID,
            source_artifact,
            manifest.redistribution_terms,
            runtime_association,
            normalizations,
        )
        .map_err(|error| {
            CompileError::new(
                "normalized_species_pack_refused",
                "pack",
                format!("runtime species-pack admission refused the compiled artifact: {error}"),
            )
        })?;
        let bytes = pack.to_bytes();
        let pack_hash = pack.content_hash();
        let decoded =
            NormalizedSpeciesPack::from_bytes_verified(pack_hash, &bytes).map_err(|error| {
                CompileError::new(
                    "self_verification_failed",
                    "pack",
                    format!("fresh species artifact failed verified runtime decode: {error}"),
                )
            })?;
        if decoded.to_bytes() != bytes {
            return Err(CompileError::new(
                "self_verification_failed",
                "pack",
                "verified runtime decode did not reproduce canonical species-pack bytes",
            ));
        }
        decisions.push(Decision::admit(
            &subject,
            "species_association_normalized",
            "molar mass and the explicit gas/ideal-gas standard-state convention compiled with complete numeric normalization receipts",
            Some(record_hash),
        ));
        decisions.push(Decision::admit(
            "pack",
            "runtime_species_pack_self_verified",
            "canonical FSSPCPK bytes decoded under their externally pinned content hash",
            Some(source_artifact),
        ));
        for decision in &mut decisions {
            decision.pack_hash = Some(pack_hash);
        }
        sort_decisions(&mut decisions);
        Ok(CompileOutput {
            pack: CompiledPack::Species(pack),
            bytes,
            decisions: std::mem::take(&mut decisions),
        })
    })();
    result.map_err(|error| {
        error
            .with_compiler_id(SPECIES_COMPILER_ID)
            .with_input_hash(source_artifact)
            .with_prior_decisions(decisions)
    })
}

#[allow(clippy::too_many_lines)] // The NASA-9 admission pipeline is intentionally auditable in order.
fn compile_nasa9_manifest(
    manifest: Manifest,
    sources: &[LoadedSource],
    source_artifact: ContentHash,
) -> Result<CompileOutput, CompileError> {
    let mut raw = RawNasa9Database::default();
    for source in sources {
        if let Err(error) = parse_nasa9_source(source, &mut raw) {
            return Err(error
                .with_input_hash(source_artifact)
                .with_prior_decisions(std::mem::take(&mut raw.decisions)));
        }
    }
    let mut decisions = std::mem::take(&mut raw.decisions);
    let result = (|| -> Result<CompileOutput, CompileError> {
        if raw.regions.is_empty() {
            return Err(CompileError::new(
                "missing_nasa9_regions",
                "sources",
                "NASA-9 model pack requires at least one region record",
            ));
        }
        for (species, region, coefficient) in raw.coefficients.keys() {
            if !raw.regions.contains_key(&(species.clone(), region.clone())) {
                return Err(CompileError::new(
                    "unknown_nasa9_region",
                    format!("species:{species}:region:{region}:coefficient:{coefficient}"),
                    "coefficient references a region that was not declared",
                ));
            }
        }
        let species_set: BTreeSet<SpeciesId> = raw
            .regions
            .keys()
            .map(|(species, _)| species.clone())
            .collect();
        if species_set.len() != 1 {
            return Err(CompileError::new(
                "multiple_species_in_model_pack",
                "sources",
                "nasa9-v1 emits one explicitly associated species artifact per manifest",
            ));
        }
        let species = species_set.into_iter().next().ok_or_else(|| {
            CompileError::new(
                "missing_nasa9_regions",
                "sources",
                "NASA-9 model pack requires one species",
            )
        })?;
        if manifest.pack_id != species.as_str() {
            return Err(CompileError::new(
                "species_pack_id_mismatch",
                "manifest:pack_id",
                format!(
                    "nasa9-v1 requires pack_id {:?} to equal the sole canonical SpeciesId {:?}",
                    manifest.pack_id,
                    species.as_str()
                ),
            ));
        }

        decisions.push(Decision::admit(
            "manifest:redistribution",
            "redistribution_permitted",
            "explicit permitted policy and nonblank retained terms admitted before compilation",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            "manifest:provenance",
            "licensed_citation_admitted",
            "nonblank source citation and license attached to every emitted model card",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            format!("species:{species}"),
            "species_pack_id_bound",
            "manifest pack_id exactly matches the fs-qty-validated source SpeciesId; this is a source association, not physical authentication",
            Some(source_artifact),
        ));

        let mut models = Vec::with_capacity(raw.regions.len());
        let mut normalizations = Vec::with_capacity(raw.regions.len().saturating_mul(12));
        let mut intervals: Vec<(f64, f64, String)> = Vec::with_capacity(raw.regions.len());
        let mut reference_pressure_bits = None;

        for ((region_species, region_id), region) in &raw.regions {
            debug_assert_eq!(region_species, &region.species);
            debug_assert_eq!(region_id, &region.id);
            let region_subject = format!("species:{region_species}:region:{region_id}");
            let temperature_lower = parse_quantity(
                &region.temperature_lower,
                &region.temperature_unit,
                &format!("{region_subject}:temperature:lower"),
            )?;
            let temperature_upper = parse_quantity(
                &region.temperature_upper,
                &region.temperature_unit,
                &format!("{region_subject}:temperature:upper"),
            )?;
            require_nasa9_dims(
                "temperature lower bound",
                NASA9_TEMPERATURE_DIMS,
                temperature_lower.dims,
                &region_subject,
            )?;
            require_nasa9_dims(
                "temperature upper bound",
                NASA9_TEMPERATURE_DIMS,
                temperature_upper.dims,
                &region_subject,
            )?;
            if temperature_lower.value <= 0.0 || temperature_lower.value >= temperature_upper.value
            {
                return Err(CompileError::new(
                    "invalid_nasa9_temperature_region",
                    &region_subject,
                    "normalized temperature bounds must be positive and strictly increasing",
                ));
            }
            let reference_pressure = parse_nasa9_linear_quantity(
                &region.reference_pressure,
                &region.reference_pressure_unit,
                &format!("{region_subject}:reference_pressure"),
            )?;
            require_nasa9_dims(
                "reference pressure",
                NASA9_PRESSURE_DIMS,
                reference_pressure.dims,
                &region_subject,
            )?;
            if reference_pressure.value <= 0.0 {
                return Err(CompileError::new(
                    "invalid_nasa9_reference_pressure",
                    &region_subject,
                    "normalized reference pressure must be positive",
                ));
            }
            match reference_pressure_bits {
                Some(expected) if expected != reference_pressure.value.to_bits() => {
                    return Err(CompileError::new(
                        "nasa9_reference_pressure_mismatch",
                        &region_subject,
                        "all regions in one species artifact require bit-identical normalized reference pressure",
                    ));
                }
                None => reference_pressure_bits = Some(reference_pressure.value.to_bits()),
                Some(_) => {}
            }

            let mut parameters = BTreeMap::new();
            let mut parsed_coefficients = Vec::with_capacity(NASA9_COEFFICIENT_NAMES.len());
            let mut card_sources = BTreeSet::from([region.source_hash, region.record_hash]);
            for name in NASA9_COEFFICIENT_NAMES {
                let coefficient = raw
                    .coefficients
                    .get(&(region_species.clone(), region_id.clone(), name.to_string()))
                    .ok_or_else(|| {
                        CompileError::new(
                            "missing_nasa9_coefficient",
                            format!("{region_subject}:coefficient:{name}"),
                            "every NASA-9 region requires exactly a0 through a8",
                        )
                    })?;
                if coefficient.source_id != region.source_id
                    || coefficient.source_hash != region.source_hash
                {
                    return Err(CompileError::new(
                        "cross_source_nasa9_region",
                        format!("{region_subject}:coefficient:{name}"),
                        "a region and all of its coefficient records must share one admitted source file",
                    ));
                }
                let parsed = parse_nasa9_linear_quantity(
                    &coefficient.value,
                    &coefficient.unit,
                    &format!("{region_subject}:coefficient:{name}"),
                )?;
                let expected = nasa9_coefficient_dims(name);
                require_nasa9_dims(
                    &format!("coefficient {name}"),
                    expected,
                    parsed.dims,
                    &region_subject,
                )?;
                parameters.insert(
                    name.to_string(),
                    LawParameter {
                        value: parsed.value,
                        dims: parsed.dims,
                    },
                );
                card_sources.insert(coefficient.source_hash);
                card_sources.insert(coefficient.record_hash);
                parsed_coefficients.push((name.to_string(), parsed));
            }
            parameters.insert(
                "reference_pressure".to_string(),
                LawParameter {
                    value: reference_pressure.value,
                    dims: reference_pressure.dims,
                },
            );
            let card = ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID.to_string()),
                law_version: NASA9_LAW_VERSION,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with(
                    "T",
                    temperature_lower.value,
                    temperature_upper.value,
                ),
                sources: card_sources.into_iter().collect(),
                provenance: Provenance {
                    source: format!(
                        "{} [source:{}] [species:{}] [region:{}]",
                        manifest.citation, region.source_id, region_species, region_id
                    ),
                    license: manifest.license.clone(),
                    artifact: Some(region.source_hash),
                },
            };
            card.validate().map_err(|error| {
                CompileError::new(
                    "nasa9_card_refused",
                    &region_subject,
                    format!("fs-matdb refused the compiled model card: {error}"),
                )
            })?;
            let model = card.content_hash();
            for (parameter, parsed) in parsed_coefficients {
                normalizations.push(model_normalization(
                    ModelNormalizationTarget::Parameter { model, parameter },
                    &parsed,
                ));
            }
            normalizations.push(model_normalization(
                ModelNormalizationTarget::Parameter {
                    model,
                    parameter: "reference_pressure".to_string(),
                },
                &reference_pressure,
            ));
            for (side, parsed) in [
                (ValidityBoundSide::Lower, &temperature_lower),
                (ValidityBoundSide::Upper, &temperature_upper),
            ] {
                normalizations.push(model_normalization(
                    ModelNormalizationTarget::ValidityBound {
                        model,
                        axis: "T".to_string(),
                        side,
                    },
                    parsed,
                ));
            }
            intervals.push((
                temperature_lower.value,
                temperature_upper.value,
                region_subject.clone(),
            ));
            decisions.push(Decision::admit(
                region_subject,
                "nasa9_region_normalized",
                "complete a0..a8 block, temperature interval, and reference pressure compiled into one immutable model card",
                Some(region.record_hash),
            ));
            models.push(card);
        }

        intervals.sort_by(|left, right| {
            left.0
                .total_cmp(&right.0)
                .then_with(|| left.1.total_cmp(&right.1))
                .then_with(|| left.2.cmp(&right.2))
        });
        for pair in intervals.windows(2) {
            if pair[0].1 > pair[1].0 {
                return Err(CompileError::new(
                    "overlapping_nasa9_regions",
                    &pair[1].2,
                    format!(
                        "temperature interval overlaps preceding region {}",
                        pair[0].2
                    ),
                ));
            }
        }

        let pack = NormalizedModelPack::new(
            manifest.pack_id,
            NASA9_COMPILER_ID,
            source_artifact,
            manifest.redistribution_terms,
            models,
            normalizations,
        )
        .map_err(|error| {
            CompileError::new(
                "normalized_model_pack_refused",
                "pack",
                format!("runtime model-pack admission refused the compiled artifact: {error}"),
            )
        })?;
        let bytes = pack.to_bytes();
        let pack_hash = pack.content_hash();
        let decoded =
            NormalizedModelPack::from_bytes_verified(pack_hash, &bytes).map_err(|error| {
                CompileError::new(
                    "self_verification_failed",
                    "pack",
                    format!("fresh artifact failed verified runtime decode: {error}"),
                )
            })?;
        if decoded.to_bytes() != bytes {
            return Err(CompileError::new(
                "self_verification_failed",
                "pack",
                "verified runtime decode did not reproduce canonical model-pack bytes",
            ));
        }
        decisions.push(Decision::admit(
            "pack",
            "runtime_model_pack_self_verified",
            "canonical FSMODPK bytes decoded under their externally pinned content hash",
            Some(source_artifact),
        ));
        for decision in &mut decisions {
            decision.pack_hash = Some(pack_hash);
        }
        sort_decisions(&mut decisions);
        Ok(CompileOutput {
            pack: CompiledPack::Model(pack),
            bytes,
            decisions: std::mem::take(&mut decisions),
        })
    })();
    result.map_err(|error| {
        error
            .with_compiler_id(NASA9_COMPILER_ID)
            .with_input_hash(source_artifact)
            .with_prior_decisions(decisions)
    })
}

const KINETICS_TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const KINETICS_FIRST_ORDER_RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn kinetics_parameter_dims(name: &str) -> Dims {
    match name {
        "activation_temperature" => KINETICS_TEMPERATURE_DIMS,
        "pre_exponential" => KINETICS_FIRST_ORDER_RATE_DIMS,
        _ => unreachable!("parser admits only the exact kinetics-v1 parameter names"),
    }
}

fn require_kinetics_dims(
    field: &str,
    expected: Dims,
    found: Dims,
    subject: &str,
) -> Result<(), CompileError> {
    if found == expected {
        Ok(())
    } else {
        Err(CompileError::new(
            "kinetics_dims_mismatch",
            subject,
            format!("{field} has dimensions {found:?}, expected {expected:?}"),
        ))
    }
}

fn parse_kinetics_linear_quantity(
    number: &str,
    unit: &str,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    let (_, _, offset) = unit_transform(unit, subject)?;
    if offset != 0.0 {
        return Err(CompileError::new(
            "affine_kinetics_parameter_unit",
            subject,
            "Arrhenius parameters require linear units; affine units are valid only for temperature-validity endpoints",
        ));
    }
    parse_linear_quantity(number, unit, subject)
}

#[allow(clippy::too_many_lines)] // The kinetics admission pipeline is intentionally auditable in order.
fn compile_kinetics_manifest(
    manifest: Manifest,
    sources: &[LoadedSource],
    source_artifact: ContentHash,
) -> Result<CompileOutput, CompileError> {
    let mut raw = RawKineticsDatabase::default();
    for source in sources {
        if let Err(error) = parse_kinetics_source(source, &mut raw) {
            return Err(error
                .with_input_hash(source_artifact)
                .with_prior_decisions(std::mem::take(&mut raw.decisions)));
        }
    }
    let mut decisions = std::mem::take(&mut raw.decisions);
    let result = (|| -> Result<CompileOutput, CompileError> {
        if raw.reactions.is_empty() {
            return Err(CompileError::new(
                "missing_kinetics_reaction",
                "sources",
                "kinetics-v1 requires exactly one first-order reaction record",
            ));
        }
        for (reaction, parameter) in raw.parameters.keys() {
            if !raw.reactions.contains_key(reaction) {
                return Err(CompileError::new(
                    "unknown_kinetics_reaction",
                    format!("reaction:{reaction}:parameter:{parameter}"),
                    "parameter references a reaction that was not declared",
                ));
            }
        }
        let (reaction_id, reaction) = raw.reactions.iter().next().ok_or_else(|| {
            CompileError::new(
                "missing_kinetics_reaction",
                "sources",
                "kinetics-v1 requires exactly one reaction",
            )
        })?;
        debug_assert_eq!(reaction_id, &reaction.reaction);
        if manifest.pack_id != reaction_id.as_str() {
            return Err(CompileError::new(
                "reaction_pack_id_mismatch",
                "manifest:pack_id",
                format!(
                    "kinetics-v1 requires pack_id {:?} to equal the sole canonical ReactionId {:?}",
                    manifest.pack_id,
                    reaction_id.as_str()
                ),
            ));
        }

        decisions.push(Decision::admit(
            "manifest:redistribution",
            "redistribution_permitted",
            "explicit permitted policy and nonblank retained terms admitted before compilation",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            "manifest:provenance",
            "licensed_citation_admitted",
            "nonblank source citation and license attached to the emitted kinetics card",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            format!("reaction:{reaction_id}"),
            "reaction_pack_id_bound",
            "manifest pack_id exactly matches the fs-qty-validated source ReactionId; this is a source association, not a mechanism or conservation proof",
            Some(source_artifact),
        ));

        let reaction_subject = format!("reaction:{reaction_id}");
        let temperature_lower = parse_quantity(
            &reaction.temperature_lower,
            &reaction.temperature_unit,
            &format!("{reaction_subject}:temperature:lower"),
        )?;
        let temperature_upper = parse_quantity(
            &reaction.temperature_upper,
            &reaction.temperature_unit,
            &format!("{reaction_subject}:temperature:upper"),
        )?;
        require_kinetics_dims(
            "temperature lower bound",
            KINETICS_TEMPERATURE_DIMS,
            temperature_lower.dims,
            &reaction_subject,
        )?;
        require_kinetics_dims(
            "temperature upper bound",
            KINETICS_TEMPERATURE_DIMS,
            temperature_upper.dims,
            &reaction_subject,
        )?;
        if temperature_lower.value <= 0.0 || temperature_lower.value >= temperature_upper.value {
            return Err(CompileError::new(
                "invalid_kinetics_temperature_range",
                &reaction_subject,
                "normalized temperature bounds must be positive and strictly increasing",
            ));
        }

        let mut parameters = BTreeMap::new();
        let mut parsed_parameters = Vec::with_capacity(KINETICS_PARAMETER_NAMES.len());
        let mut card_sources = BTreeSet::from([reaction.source_hash, reaction.record_hash]);
        for name in KINETICS_PARAMETER_NAMES {
            let parameter = raw
                .parameters
                .get(&(reaction_id.clone(), name.to_string()))
                .ok_or_else(|| {
                    CompileError::new(
                        "missing_kinetics_parameter",
                        format!("{reaction_subject}:parameter:{name}"),
                        "kinetics-v1 requires exactly activation_temperature and pre_exponential",
                    )
                })?;
            if parameter.source_id != reaction.source_id
                || parameter.source_hash != reaction.source_hash
            {
                return Err(CompileError::new(
                    "cross_source_kinetics_reaction",
                    format!("{reaction_subject}:parameter:{name}"),
                    "a reaction and both of its parameter records must share one admitted source file",
                ));
            }
            let parsed = parse_kinetics_linear_quantity(
                &parameter.value,
                &parameter.unit,
                &format!("{reaction_subject}:parameter:{name}"),
            )?;
            let expected = kinetics_parameter_dims(name);
            require_kinetics_dims(name, expected, parsed.dims, &reaction_subject)?;
            match name {
                "activation_temperature" if parsed.value < 0.0 => {
                    return Err(CompileError::new(
                        "invalid_activation_temperature",
                        &reaction_subject,
                        "activation temperature must be nonnegative",
                    ));
                }
                "pre_exponential" if parsed.value <= 0.0 => {
                    return Err(CompileError::new(
                        "invalid_pre_exponential",
                        &reaction_subject,
                        "first-order pre-exponential factor must be positive",
                    ));
                }
                _ => {}
            }
            parameters.insert(
                name.to_string(),
                LawParameter {
                    value: parsed.value,
                    dims: parsed.dims,
                },
            );
            card_sources.insert(parameter.source_hash);
            card_sources.insert(parameter.record_hash);
            parsed_parameters.push((name.to_string(), parsed));
        }

        let card = ConstitutiveModelCard {
            law: LawId(KINETICS_LAW_ID.to_string()),
            law_version: KINETICS_LAW_VERSION,
            parameters,
            state_schema_version: KINETICS_STATE_SCHEMA_VERSION,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with(
                "T",
                temperature_lower.value,
                temperature_upper.value,
            ),
            sources: card_sources.into_iter().collect(),
            provenance: Provenance {
                source: format!(
                    "{} [source:{}] [reaction:{}] [rate-basis:first-order]",
                    manifest.citation, reaction.source_id, reaction_id
                ),
                license: manifest.license.clone(),
                artifact: Some(reaction.source_hash),
            },
        };
        card.validate().map_err(|error| {
            CompileError::new(
                "kinetics_card_refused",
                &reaction_subject,
                format!("fs-matdb refused the compiled model card: {error}"),
            )
        })?;
        let model = card.content_hash();
        let mut normalizations = Vec::with_capacity(KINETICS_PARAMETER_NAMES.len() + 2);
        for (parameter, parsed) in parsed_parameters {
            normalizations.push(model_normalization(
                ModelNormalizationTarget::Parameter { model, parameter },
                &parsed,
            ));
        }
        for (side, parsed) in [
            (ValidityBoundSide::Lower, &temperature_lower),
            (ValidityBoundSide::Upper, &temperature_upper),
        ] {
            normalizations.push(model_normalization(
                ModelNormalizationTarget::ValidityBound {
                    model,
                    axis: "T".to_string(),
                    side,
                },
                parsed,
            ));
        }
        decisions.push(Decision::admit(
            reaction_subject,
            "kinetics_reaction_normalized",
            "explicit first-order basis, complete Arrhenius parameter block, and temperature interval compiled into one immutable model card; no executor or conservation claim is implied",
            Some(reaction.record_hash),
        ));

        let pack = NormalizedModelPack::new(
            manifest.pack_id,
            KINETICS_COMPILER_ID,
            source_artifact,
            manifest.redistribution_terms,
            vec![card],
            normalizations,
        )
        .map_err(|error| {
            CompileError::new(
                "normalized_model_pack_refused",
                "pack",
                format!("runtime model-pack admission refused the compiled artifact: {error}"),
            )
        })?;
        let bytes = pack.to_bytes();
        let pack_hash = pack.content_hash();
        let decoded =
            NormalizedModelPack::from_bytes_verified(pack_hash, &bytes).map_err(|error| {
                CompileError::new(
                    "self_verification_failed",
                    "pack",
                    format!("fresh artifact failed verified runtime decode: {error}"),
                )
            })?;
        if decoded.to_bytes() != bytes {
            return Err(CompileError::new(
                "self_verification_failed",
                "pack",
                "verified runtime decode did not reproduce canonical model-pack bytes",
            ));
        }
        decisions.push(Decision::admit(
            "pack",
            "runtime_model_pack_self_verified",
            "canonical FSMODPK bytes decoded under their externally pinned content hash",
            Some(source_artifact),
        ));
        for decision in &mut decisions {
            decision.pack_hash = Some(pack_hash);
        }
        sort_decisions(&mut decisions);
        Ok(CompileOutput {
            pack: CompiledPack::Model(pack),
            bytes,
            decisions: std::mem::take(&mut decisions),
        })
    })();
    result.map_err(|error| {
        error
            .with_compiler_id(KINETICS_COMPILER_ID)
            .with_input_hash(source_artifact)
            .with_prior_decisions(decisions)
    })
}

fn compile_interface_identity(
    raw: &RawDatabase,
) -> Result<(SurfaceSpec, SurfaceSpec, SystemContext, [ContentHash; 3]), CompileError> {
    let surface_a = raw.surface_a.as_ref().ok_or_else(|| {
        CompileError::new(
            "missing_interface_surface",
            "interface:surface_a",
            "interface-tsv-v1 requires exactly one surface_a record",
        )
    })?;
    let surface_b = raw.surface_b.as_ref().ok_or_else(|| {
        CompileError::new(
            "missing_interface_surface",
            "interface:surface_b",
            "interface-tsv-v1 requires exactly one surface_b record",
        )
    })?;
    let context = raw.context.as_ref().ok_or_else(|| {
        CompileError::new(
            "missing_interface_context",
            "interface:context",
            "interface-tsv-v1 requires exactly one context record",
        )
    })?;
    let to_surface = |surface: &RawSurfaceSpec| SurfaceSpec {
        material: MaterialStateId {
            chemistry: surface.chemistry.clone(),
            phase: surface.phase.clone(),
            process: surface.process.clone(),
            revision: surface.revision,
        },
        texture_frame: surface.texture_frame.clone(),
    };
    Ok((
        to_surface(surface_a),
        to_surface(surface_b),
        SystemContext {
            medium: context.medium.clone(),
            third_body: context.third_body.clone(),
            environment: context.environment.clone(),
            history: context.history.clone(),
        },
        [
            surface_a.record_hash,
            surface_b.record_hash,
            context.record_hash,
        ],
    ))
}

#[allow(clippy::too_many_lines)] // The ordered admission pipeline is easier to audit in one pass.
fn compile_manifest(manifest_path: &Path) -> Result<CompileOutput, CompileError> {
    let manifest_bytes = read_bounded_regular(manifest_path, MAX_MANIFEST_BYTES, "manifest")?;
    let manifest_snapshot = hash_domain(SOURCE_ENVELOPE_DOMAIN, &manifest_bytes);
    let manifest_text = std::str::from_utf8(&manifest_bytes).map_err(|_| {
        CompileError::new(
            "invalid_manifest_encoding",
            "manifest",
            "manifest must be valid UTF-8",
        )
        .with_input_hash(manifest_snapshot)
    })?;
    let manifest =
        parse_manifest(manifest_text).map_err(|error| error.with_input_hash(manifest_snapshot))?;
    let profile =
        source_profile(&manifest).map_err(|error| error.with_input_hash(manifest_snapshot))?;
    let compiler_id = match (manifest.version, profile) {
        (4, SourceProfile::Material) => HARDNESS_COMPILER_ID,
        (4, SourceProfile::Interface) => HARDNESS_INTERFACE_COMPILER_ID,
        (3, SourceProfile::Material) => TYPED_AXIS_COMPILER_ID,
        (3, SourceProfile::Interface) => TYPED_AXIS_INTERFACE_COMPILER_ID,
        (2, SourceProfile::Material) => MAPPED_COMPILER_ID,
        (2, SourceProfile::Interface) => MAPPED_INTERFACE_COMPILER_ID,
        _ => profile.compiler_id(),
    };
    let sources = load_sources(manifest_path, &manifest).map_err(|error| {
        error
            .with_compiler_id(compiler_id)
            .with_input_hash(manifest_snapshot)
    })?;
    let source_artifact = source_envelope_hash(&manifest, &sources);
    match profile {
        SourceProfile::Nasa9 => {
            return compile_nasa9_manifest(manifest, &sources, source_artifact)
                .map_err(|error| error.with_compiler_id(compiler_id));
        }
        SourceProfile::Kinetics => {
            return compile_kinetics_manifest(manifest, &sources, source_artifact)
                .map_err(|error| error.with_compiler_id(compiler_id));
        }
        SourceProfile::Species => {
            return compile_species_manifest(manifest, &sources, source_artifact)
                .map_err(|error| error.with_compiler_id(compiler_id));
        }
        SourceProfile::Material | SourceProfile::Interface => {}
    }
    let mut raw = RawDatabase::default();
    for source in &sources {
        if let Err(error) = parse_source(source, &mut raw, profile) {
            return Err(error
                .with_compiler_id(compiler_id)
                .with_input_hash(source_artifact)
                .with_prior_decisions(std::mem::take(&mut raw.decisions)));
        }
    }
    let mut decisions = std::mem::take(&mut raw.decisions);
    let result = (|| -> Result<CompileOutput, CompileError> {
        validate_raw_references(&raw)?;
        validate_name_mappings(&manifest, &raw)?;
        preflight_provenance_budget(&manifest, &raw)?;
        preflight_normalization_budget(&raw)?;

        let mut claims = ClaimSet::new();
        let mut observation_ids = BTreeMap::new();
        let mut components = BTreeMap::new();
        let mut normalizations = Vec::new();
        decisions.push(Decision::admit(
            "manifest:redistribution",
            "redistribution_permitted",
            "explicit permitted policy and nonblank retained terms admitted before compilation",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            "manifest:provenance",
            "licensed_citation_admitted",
            "nonblank source citation and license attached to every emitted datum",
            Some(source_artifact),
        ));
        let mut validity_dims: BTreeMap<String, Dims> = BTreeMap::new();

        for (local_id, observation) in &raw.observations {
            debug_assert_eq!(local_id, &observation.id);
            let id = claims
                .register_observation(ObservationDataset {
                    specimen: observation.specimen.clone(),
                    method: observation.method.clone(),
                    artifact: observation.record_hash,
                    caveats: observation.caveats.clone(),
                    provenance: provenance_for(
                        &manifest,
                        &observation.source_id,
                        observation.source_hash,
                    ),
                })
                .map_err(|error| {
                    CompileError::new(
                        "observation_admission_refused",
                        format!("observation:{local_id}"),
                        error.to_string(),
                    )
                })?;
            observation_ids.insert(local_id.clone(), id);
            decisions.push(Decision::admit(
                format!("observation:{local_id}"),
                "observation_normalized",
                "registered provenance-complete observation dataset",
                Some(observation.record_hash),
            ));
        }

        for (local_id, raw_claim) in &raw.claims {
            debug_assert_eq!(local_id, &raw_claim.id);
            let frame = raw.frames.get(local_id);
            let mapping = manifest.properties.get(&raw_claim.property);
            let kind = mapping.and_then(|(_, kind)| *kind);
            let (mut value, parsed_components) = compile_value(raw_claim, kind, &manifest)?;
            let value_dims = value.dims();
            let (validity, validity_receipts) = compile_validity(
                local_id,
                raw.validities.get(local_id),
                &mut validity_dims,
                &manifest,
            )?;
            if let PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                ..
            } = &value
            {
                let declared = raw
                    .validities
                    .get(local_id)
                    .and_then(|axes| axes.get(abscissa))
                    .ok_or_else(|| {
                        CompileError::new(
                            "missing_curve_validity",
                            format!("claim:{local_id}"),
                            "curve abscissa must have an explicit validity interval",
                        )
                    })?;
                let parsed = parse_declared_quantity(
                    &declared.lower,
                    &declared.unit,
                    manifest.axis_kinds.get(abscissa).copied().flatten(),
                    &format!("claim:{local_id}:validity:{abscissa}:lower"),
                )?;
                if parsed.dims != *abscissa_dims {
                    return Err(CompileError::new(
                        "curve_validity_dims_mismatch",
                        format!("claim:{local_id}"),
                        "curve abscissa dimensions disagree with its validity axis",
                    ));
                }
            }
            if let PropertyValue::Curve { abscissa, .. } = &mut value {
                *abscissa = manifest.axis(abscissa).to_owned();
            }
            let raw_uncertainty = raw.uncertainties.get(local_id).ok_or_else(|| {
                CompileError::new(
                    "missing_uncertainty",
                    format!("claim:{local_id}:uncertainty"),
                    "each source claim must explicitly declare its uncertainty, including unstated",
                )
            })?;
            let (uncertainty, uncertainty_quantity) =
                compile_uncertainty(local_id, raw_uncertainty, value_dims, kind)?;
            let mut observations = raw_claim
                .observations
                .iter()
                .map(|reference| {
                    observation_ids.get(reference).copied().ok_or_else(|| {
                        CompileError::new(
                            "unknown_observation",
                            format!("claim:{local_id}"),
                            format!("observation reference {reference:?} is not registered"),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            observations.sort_unstable();
            let component_observations: BTreeSet<ObservationId> =
                observations.iter().copied().collect();
            let mut key = PropertyKey::with_quantity(
                mapping.map_or(raw_claim.property.as_str(), |(name, _)| name.as_str()),
                kind.map_or(QuantitySpec::dimensional(value_dims), |kind| {
                    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
                }),
            );
            if let Some(test) = manifest.hardness_tests.get(local_id) {
                let subject = format!("claim:{local_id}:hardness-test");
                let observation =
                    observation_ids
                        .get(&test.observation)
                        .copied()
                        .ok_or_else(|| {
                            CompileError::new(
                                "unknown_observation",
                                &subject,
                                "hardness test observation is not registered",
                            )
                        })?;
                let mut holds = Vec::with_capacity(test.holds.len());
                for [force, force_unit, dwell, dwell_unit] in &test.holds {
                    let force = parse_linear_quantity(force, force_unit, &subject)?;
                    let dwell = parse_linear_quantity(dwell, dwell_unit, &subject)?;
                    holds.push(
                        fs_matdb::HardnessLoadStep::new(
                            fs_qty::QtyAny::new(force.value, force.dims),
                            fs_qty::QtyAny::new(dwell.value, dwell.dims),
                        )
                        .map_err(|error| {
                            CompileError::new("invalid_hardness_test", &subject, error.to_string())
                        })?,
                    );
                }
                let context = fs_matdb::HardnessTestContext::new(
                    &test.indenter,
                    holds,
                    &test.protocol,
                    observation,
                )
                .map_err(|error| {
                    CompileError::new("invalid_hardness_test", &subject, error.to_string())
                })?;
                key = key.with_hardness_test(context).map_err(|error| {
                    CompileError::new("invalid_hardness_test", &subject, error.to_string())
                })?;
            }
            let id = claims
                .insert_claim(PropertyClaim {
                    key,
                    value,
                    validity,
                    uncertainty,
                    interpolation: raw_claim.interpolation,
                    observations,
                    provenance: provenance_for(
                        &manifest,
                        &raw_claim.source_id,
                        raw_claim.source_hash,
                    ),
                })
                .map_err(|error| {
                    CompileError::new(
                        "claim_admission_refused",
                        format!("claim:{local_id}"),
                        error.to_string(),
                    )
                })?;
            for (component, quantity) in parsed_components {
                let member = statistic_member(id, component);
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ClaimValue(member),
                    &quantity,
                    frame,
                ));
                components.insert(
                    (local_id.clone(), component),
                    ComponentMeta {
                        member,
                        dims: quantity.dims,
                        scale: quantity.scale,
                        unit: quantity.unit.clone(),
                        observations: component_observations.clone(),
                    },
                );
                decisions.push(Decision::admit(
                    normalization_subject(local_id, component),
                    "quantity_normalized",
                    "explicit source quantity converted to the six-base SI target",
                    Some(quantity.literal_hash),
                ));
            }
            if let Some(quantity) = uncertainty_quantity {
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ClaimUncertainty { claim: id },
                    &quantity,
                    frame,
                ));
                decisions.push(Decision::admit(
                    format!("normalization:{local_id}:uncertainty"),
                    "uncertainty_normalized",
                    "uncertainty magnitude converted without an affine translation",
                    Some(quantity.literal_hash),
                ));
            }
            decisions.push(Decision::admit(
                format!("claim:{local_id}:uncertainty"),
                "uncertainty_policy_admitted",
                "source uncertainty kind and confidence were admitted without inference",
                Some(raw_uncertainty.source_hash),
            ));
            for (axis, lower, upper) in validity_receipts {
                let validity_source_hash = raw
                    .validities
                    .get(local_id)
                    .and_then(|values| values.get(&axis))
                    .map(|value| value.source_hash)
                    .ok_or_else(|| {
                        CompileError::new(
                            "internal_linkage_error",
                            format!("claim:{local_id}:validity:{axis}"),
                            "compiled validity receipt lost its admitted source record",
                        )
                    })?;
                let axis = manifest.axis(&axis).to_owned();
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ValidityBound {
                        claim: id,
                        axis: axis.clone(),
                        side: ValidityBoundSide::Lower,
                    },
                    &lower,
                    None,
                ));
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ValidityBound {
                        claim: id,
                        axis: axis.clone(),
                        side: ValidityBoundSide::Upper,
                    },
                    &upper,
                    None,
                ));
                decisions.push(Decision::admit(
                    format!("normalization:{local_id}:validity:{axis}"),
                    "validity_normalized",
                    "both validity endpoints converted under one explicit unit transform",
                    Some(validity_source_hash),
                ));
            }
            decisions.push(Decision::admit(
                format!("claim:{local_id}"),
                "claim_normalized",
                "inserted one immutable dims-checked scalar/curve claim",
                Some(raw_claim.record_hash),
            ));
            if let Some(frame) = frame {
                decisions.push(Decision::admit(
                    format!("claim:{local_id}:frame"),
                    "frame_pair_retained",
                    "paired source and target frames retained as normalization metadata",
                    Some(frame.source_hash),
                ));
            }
        }

        if manifest.version >= 2 {
            for (source, (target, kind)) in &manifest.properties {
                decisions.push(Decision::admit(
                    format!("mapping:property:{source}"),
                    "property_schema_normalized",
                    format!(
                        "source property {source:?} maps to {target:?} with declared kind {kind:?}"
                    ),
                    Some(source_artifact),
                ));
            }
            for (source, target) in &manifest.axes {
                decisions.push(Decision::admit(
                    format!("mapping:axis:{source}"), "axis_name_normalized",
                    format!("source axis {source:?} maps to {target:?}; declared quantity {:?}; values retain explicit normalization", manifest.axis_kinds.get(source)),
                    Some(source_artifact),
                ));
            }
        }

        let mut joint_statistics = Vec::new();
        for ((observation, block_id), raw_joint) in &raw.joints {
            let observation_id = observation_ids.get(observation).copied().ok_or_else(|| {
                CompileError::new(
                    "unknown_observation",
                    format!("joint:{observation}:{block_id}"),
                    "joint block references an unknown observation",
                )
            })?;
            let (joint, receipts, joint_decisions) =
                compile_joint(raw_joint, observation_id, &components)?;
            joint_statistics.push(joint);
            normalizations.extend(receipts);
            decisions.extend(joint_decisions);
            decisions.push(Decision::admit(
                format!("joint:{observation}:{block_id}"),
                "joint_statistics_normalized",
                "member order and packed matrices were permuted together before SI scaling",
                Some(raw_joint.source_hash),
            ));
        }

        let material_pack = NormalizedPack::new(
            manifest.pack_id,
            compiler_id,
            source_artifact,
            manifest.redistribution_terms,
            claims,
            joint_statistics,
            normalizations,
        )
        .map_err(|error| {
            CompileError::new(
                "normalized_pack_refused",
                "pack",
                format!("runtime pack admission refused the compiled artifact: {error}"),
            )
        })?;
        let pack = match profile {
            SourceProfile::Material => CompiledPack::Material(material_pack),
            SourceProfile::Interface => {
                let (surface_a, surface_b, context, identity_hashes) =
                    compile_interface_identity(&raw)?;
                decisions.push(Decision::admit(
                    "interface:surface_a",
                    "ordered_surface_identity_admitted",
                    "first material state and texture frame retained as the ordered A surface",
                    Some(identity_hashes[0]),
                ));
                decisions.push(Decision::admit(
                    "interface:surface_b",
                    "ordered_surface_identity_admitted",
                    "second material state and texture frame retained as the ordered B surface",
                    Some(identity_hashes[1]),
                ));
                decisions.push(Decision::admit(
                    "interface:context",
                    "interface_context_admitted",
                    "medium, optional third body, environment, and named history retained without inference",
                    Some(identity_hashes[2]),
                ));
                let interface_pack = NormalizedInterfacePack::new(
                    surface_a,
                    surface_b,
                    context,
                    material_pack,
                )
                .map_err(|error| {
                    CompileError::new(
                        "normalized_interface_pack_refused",
                        "pack",
                        format!(
                            "runtime interface-pack admission refused the compiled artifact: {error}"
                        ),
                    )
                })?;
                CompiledPack::Interface(interface_pack)
            }
            SourceProfile::Nasa9 | SourceProfile::Kinetics | SourceProfile::Species => {
                unreachable!("specialized profiles return before material claim compilation")
            }
        };
        let bytes = pack.to_bytes();
        let pack_hash = pack.content_hash();
        pack.verify_bytes(pack_hash, &bytes).map_err(|error| {
            CompileError::new(
                "self_verification_failed",
                "pack",
                format!("fresh artifact failed verified runtime decode: {error}"),
            )
        })?;
        if pack.to_bytes() != bytes {
            return Err(CompileError::new(
                "self_verification_failed",
                "pack",
                "verified runtime decode did not reproduce canonical bytes",
            ));
        }
        decisions.push(Decision::admit(
            "pack",
            match profile {
                SourceProfile::Material => "runtime_pack_self_verified",
                SourceProfile::Interface => "runtime_interface_pack_self_verified",
                SourceProfile::Nasa9 | SourceProfile::Kinetics | SourceProfile::Species => {
                    unreachable!("specialized profiles return before material claim compilation")
                }
            },
            match profile {
                SourceProfile::Material => {
                    "canonical bytes decoded under their externally pinned content hash"
                }
                SourceProfile::Interface => {
                    "canonical FSINTPK bytes decoded under their externally pinned content hash"
                }
                SourceProfile::Nasa9 | SourceProfile::Kinetics | SourceProfile::Species => {
                    unreachable!("specialized profiles return before material claim compilation")
                }
            },
            Some(source_artifact),
        ));
        for decision in &mut decisions {
            decision.pack_hash = Some(pack_hash);
        }
        sort_decisions(&mut decisions);
        Ok(CompileOutput {
            pack,
            bytes,
            decisions: std::mem::take(&mut decisions),
        })
    })();
    result.map_err(|error| {
        error
            .with_compiler_id(compiler_id)
            .with_input_hash(source_artifact)
            .with_prior_decisions(decisions)
    })
}

fn source_envelope_hash(manifest: &Manifest, sources: &[LoadedSource]) -> ContentHash {
    let mut payload = Vec::new();
    if manifest.version >= 2 {
        push_part(&mut payload, &(sources.len() as u64).to_le_bytes());
        push_part(
            &mut payload,
            &(manifest.properties.len() as u64).to_le_bytes(),
        );
        push_part(&mut payload, &(manifest.axes.len() as u64).to_le_bytes());
    }
    for part in [
        manifest.pack_id.as_bytes(),
        manifest.redistribution_terms.as_bytes(),
        manifest.citation.as_bytes(),
        manifest.license.as_bytes(),
    ] {
        push_part(&mut payload, part);
    }
    for source in sources {
        push_part(&mut payload, source.spec.id.as_bytes());
        push_part(&mut payload, source.spec.profile.as_bytes());
        push_part(&mut payload, &source.bytes);
    }
    if manifest.version >= 2 {
        for (source, (target, kind)) in &manifest.properties {
            push_part(&mut payload, source.as_bytes());
            push_part(&mut payload, target.as_bytes());
            match kind {
                Some(kind) => push_part(
                    &mut payload,
                    &QuantitySpec::semantic(SemanticType::new(*kind, ValueForm::Static))
                        .canonical_bytes(),
                ),
                None => push_part(&mut payload, b"dimensional"),
            }
        }
        for (source, target) in &manifest.axes {
            push_part(&mut payload, source.as_bytes());
            push_part(&mut payload, target.as_bytes());
            if manifest.version >= 3 {
                match manifest.axis_kinds[source] {
                    Some(kind) => push_part(
                        &mut payload,
                        &QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
                            .canonical_bytes(),
                    ),
                    None => push_part(&mut payload, b"dimensional"),
                }
            }
        }
        if manifest.version >= 4 {
            push_part(
                &mut payload,
                &(manifest.hardness_tests.len() as u64).to_le_bytes(),
            );
            for (claim, test) in &manifest.hardness_tests {
                for text in [claim, &test.observation, &test.indenter, &test.protocol] {
                    push_part(&mut payload, text.as_bytes());
                }
                push_part(&mut payload, &(test.holds.len() as u64).to_le_bytes());
                for hold in &test.holds {
                    for text in hold {
                        push_part(&mut payload, text.as_bytes());
                    }
                }
            }
        }
        hash_domain(
            if manifest.version >= 4 {
                "org.frankensim.xtask.matdb-pack.source-envelope.v4"
            } else if manifest.version == 3 {
                "org.frankensim.xtask.matdb-pack.source-envelope.v3"
            } else {
                "org.frankensim.xtask.matdb-pack.source-envelope.v2"
            },
            &payload,
        )
    } else {
        hash_domain(SOURCE_ENVELOPE_DOMAIN, &payload)
    }
}

fn validate_name_mappings(manifest: &Manifest, raw: &RawDatabase) -> Result<(), CompileError> {
    for claim in manifest.hardness_tests.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unused_hardness_test",
                format!("claim:{claim}"),
                "hardness test must refer to an imported claim",
            ));
        }
    }
    let properties: BTreeSet<&str> = raw
        .claims
        .values()
        .map(|claim| claim.property.as_str())
        .collect();
    let axes: BTreeSet<&str> = raw
        .validities
        .values()
        .flat_map(|axes| axes.keys().map(String::as_str))
        .collect();
    for (label, names, mappings) in [
        (
            "property",
            &properties,
            manifest
                .properties
                .iter()
                .map(|(source, (target, _))| (source.as_str(), target.as_str()))
                .collect::<BTreeMap<_, _>>(),
        ),
        (
            "axis",
            &axes,
            manifest
                .axes
                .iter()
                .map(|(source, target)| (source.as_str(), target.as_str()))
                .collect(),
        ),
    ] {
        if let Some(source) = mappings.keys().find(|source| !names.contains(**source)) {
            return Err(CompileError::new(
                "unused_name_mapping",
                format!("mapping:{label}:{source}"),
                "mapping must reference a declared source name",
            ));
        }
        let mut targets = BTreeSet::new();
        for source in names {
            let target = mappings.get(source).copied().unwrap_or(source);
            if !targets.insert(target) {
                return Err(CompileError::new(
                    "normalization_name_collision",
                    format!("mapping:{label}:{target}"),
                    "distinct source names may not collapse into one runtime name",
                ));
            }
        }
    }
    Ok(())
}

fn provenance_for(manifest: &Manifest, source_id: &str, source_hash: ContentHash) -> Provenance {
    Provenance {
        source: format!("{} [source:{source_id}]", manifest.citation),
        license: manifest.license.clone(),
        artifact: Some(source_hash),
    }
}

fn validate_raw_references(raw: &RawDatabase) -> Result<(), CompileError> {
    if raw.observations.is_empty() {
        return Err(CompileError::new(
            "missing_observations",
            "sources",
            "material pack must contain at least one observation",
        ));
    }
    if raw.claims.is_empty() {
        return Err(CompileError::new(
            "missing_claims",
            "sources",
            "material pack must contain at least one scalar or curve claim",
        ));
    }
    for claim in raw.uncertainties.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:uncertainty"),
                "uncertainty record references an unknown claim",
            ));
        }
    }
    for claim in raw.validities.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:validity"),
                "validity record references an unknown claim",
            ));
        }
    }
    for claim in raw.frames.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:frame"),
                "frame record references an unknown claim",
            ));
        }
    }
    for raw_claim in raw.claims.values() {
        for observation in &raw_claim.observations {
            if !raw.observations.contains_key(observation) {
                return Err(CompileError::new(
                    "unknown_observation",
                    format!("claim:{}", raw_claim.id),
                    format!("observation reference {observation:?} is not declared"),
                ));
            }
        }
    }
    for raw_joint in raw.joints.values() {
        if !raw.observations.contains_key(&raw_joint.observation) {
            return Err(CompileError::new(
                "unknown_observation",
                format!("joint:{}:{}", raw_joint.observation, raw_joint.block_id),
                "joint block references an unknown observation",
            ));
        }
    }
    Ok(())
}

fn preflight_provenance_budget(manifest: &Manifest, raw: &RawDatabase) -> Result<(), CompileError> {
    let source_ids = raw
        .observations
        .values()
        .map(|value| value.source_id.as_str())
        .chain(raw.claims.values().map(|value| value.source_id.as_str()));
    let mut repeated = 0usize;
    for source_id in source_ids {
        // Two wire string lengths, the artifact tag/hash, and the literal
        // " [source:...]" expansion are included. This is a conservative
        // lower-level budget solely for repeated global provenance; ordinary
        // source payloads remain bounded independently.
        let record = manifest
            .citation
            .len()
            .checked_add(manifest.license.len())
            .and_then(|value| value.checked_add(source_id.len()))
            .and_then(|value| value.checked_add(64))
            .ok_or_else(|| {
                CompileError::new(
                    "resource_limit",
                    "provenance",
                    "repeated provenance byte estimate overflowed",
                )
            })?;
        repeated = repeated.checked_add(record).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                "provenance",
                "aggregate repeated provenance bytes overflowed",
            )
        })?;
        if repeated > MAX_REPEATED_PROVENANCE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                "provenance",
                format!(
                    "repeated citation/license expansion exceeds \
                     {MAX_REPEATED_PROVENANCE_BYTES} bytes"
                ),
            ));
        }
    }
    Ok(())
}

fn preflight_normalization_budget(raw: &RawDatabase) -> Result<(), CompileError> {
    let mut receipts = 0usize;
    let mut add = |count: usize, subject: &str| -> Result<(), CompileError> {
        receipts = receipts.checked_add(count).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                subject,
                "normalization receipt count overflowed",
            )
        })?;
        if receipts > MAX_NORMALIZATION_RECEIPTS {
            return Err(CompileError::new(
                "resource_limit",
                "normalizations",
                format!(
                    "compiled source requires {receipts} receipts; schema v1 admits at most \
                     {MAX_NORMALIZATION_RECEIPTS}"
                ),
            ));
        }
        Ok(())
    };
    for (claim_id, claim) in &raw.claims {
        let value_receipts = match &claim.value {
            RawClaimValue::Scalar { .. } => 1,
            RawClaimValue::Curve { points, .. } => {
                points.len().checked_mul(2).ok_or_else(|| {
                    CompileError::new(
                        "resource_limit",
                        format!("claim:{claim_id}"),
                        "curve normalization count overflowed",
                    )
                })?
            }
        };
        add(value_receipts, &format!("claim:{claim_id}"))?;
        if raw
            .uncertainties
            .get(claim_id)
            .is_some_and(|value| !matches!(&value.kind, RawUncertaintyKind::Unstated))
        {
            add(1, &format!("claim:{claim_id}:uncertainty"))?;
        }
        let validity_receipts = raw
            .validities
            .get(claim_id)
            .map_or(0, |axes| axes.len().saturating_mul(2));
        add(validity_receipts, &format!("claim:{claim_id}:validity"))?;
    }
    for ((observation, block_id), joint) in &raw.joints {
        let packed = joint
            .members
            .len()
            .checked_mul(joint.members.len() + 1)
            .and_then(|value| value.checked_div(2))
            .ok_or_else(|| {
                CompileError::new(
                    "resource_limit",
                    format!("joint:{observation}:{block_id}"),
                    "joint normalization count overflowed",
                )
            })?;
        add(packed, &format!("joint:{observation}:{block_id}"))?;
    }
    Ok(())
}

fn compile_value(
    claim: &RawClaim,
    kind: Option<QuantityKind>,
    manifest: &Manifest,
) -> Result<(PropertyValue, Vec<(LocalComponent, ParsedQuantity)>), CompileError> {
    let subject = format!("claim:{}", claim.id);
    match &claim.value {
        RawClaimValue::Scalar { number, unit } => {
            if claim.interpolation == InterpolationPolicy::LinearInside {
                return Err(CompileError::new(
                    "incompatible_interpolation",
                    &subject,
                    "scalar claims admit only constant or tabulated interpolation",
                ));
            }
            let quantity =
                parse_declared_quantity(number, unit, kind, &format!("{subject}:value"))?;
            Ok((
                PropertyValue::Scalar {
                    value: quantity.value,
                    dims: quantity.dims,
                },
                vec![(LocalComponent::Scalar, quantity)],
            ))
        }
        RawClaimValue::Curve {
            abscissa,
            abscissa_unit,
            ordinate_unit,
            points,
        } => {
            if claim.interpolation == InterpolationPolicy::ConstantWithinValidity {
                return Err(CompileError::new(
                    "incompatible_interpolation",
                    &subject,
                    "curve claims admit only linear or tabulated interpolation",
                ));
            }
            let mut knots = Vec::with_capacity(points.len());
            let mut components = Vec::with_capacity(points.len() * 2);
            let mut abscissa_dims = None;
            let mut ordinate_dims = None;
            for (index, (x, y)) in points.iter().enumerate() {
                let knot = u32::try_from(index).map_err(|_| {
                    CompileError::new(
                        "resource_limit",
                        &subject,
                        "curve knot index exceeds the runtime wire format",
                    )
                })?;
                let axis_kind = manifest.axis_kinds.get(abscissa).copied().flatten();
                let parsed_x = parse_declared_quantity(
                    x,
                    abscissa_unit,
                    axis_kind,
                    &format!("{subject}:curve:x:{index}"),
                )?;
                let parsed_y = parse_declared_quantity(
                    y,
                    ordinate_unit,
                    kind,
                    &format!("{subject}:curve:y:{index}"),
                )?;
                if let Some(expected) = abscissa_dims {
                    if expected != parsed_x.dims {
                        return Err(CompileError::new(
                            "curve_dims_mismatch",
                            &subject,
                            "curve abscissa dimensions changed between knots",
                        ));
                    }
                } else {
                    abscissa_dims = Some(parsed_x.dims);
                }
                if let Some(expected) = ordinate_dims {
                    if expected != parsed_y.dims {
                        return Err(CompileError::new(
                            "curve_dims_mismatch",
                            &subject,
                            "curve ordinate dimensions changed between knots",
                        ));
                    }
                } else {
                    ordinate_dims = Some(parsed_y.dims);
                }
                knots.push((parsed_x.value, parsed_y.value));
                components.push((LocalComponent::CurveAbscissa(knot), parsed_x));
                components.push((LocalComponent::CurveOrdinate(knot), parsed_y));
            }
            if !knots.windows(2).all(|pair| pair[0].0 < pair[1].0) {
                return Err(CompileError::new(
                    "malformed_curve",
                    &subject,
                    "normalized curve abscissae must be strictly increasing",
                ));
            }
            let abscissa_dims = abscissa_dims.ok_or_else(|| {
                CompileError::new(
                    "malformed_curve",
                    &subject,
                    "curve has no parsed abscissa dimensions",
                )
            })?;
            let ordinate_dims = ordinate_dims.ok_or_else(|| {
                CompileError::new(
                    "malformed_curve",
                    &subject,
                    "curve has no parsed ordinate dimensions",
                )
            })?;
            Ok((
                PropertyValue::Curve {
                    abscissa: abscissa.clone(),
                    abscissa_dims,
                    knots,
                    dims: ordinate_dims,
                },
                components,
            ))
        }
    }
}

fn compile_validity(
    claim_id: &str,
    axes: Option<&BTreeMap<String, RawValidity>>,
    global_dims: &mut BTreeMap<String, Dims>,
    manifest: &Manifest,
) -> Result<
    (
        ValidityDomain,
        Vec<(String, ParsedQuantity, ParsedQuantity)>,
    ),
    CompileError,
> {
    let mut domain = ValidityDomain::unconstrained();
    let mut receipts = Vec::new();
    for (axis, raw) in axes.into_iter().flatten() {
        debug_assert_eq!(axis, &raw.axis);
        let kind = manifest.axis_kinds.get(axis).copied().flatten();
        let lower = parse_declared_quantity(
            &raw.lower,
            &raw.unit,
            kind,
            &format!("claim:{claim_id}:validity:{axis}:lower"),
        )?;
        let upper = parse_declared_quantity(
            &raw.upper,
            &raw.unit,
            kind,
            &format!("claim:{claim_id}:validity:{axis}:upper"),
        )?;
        if lower.dims != upper.dims {
            return Err(CompileError::new(
                "validity_dims_mismatch",
                format!("claim:{claim_id}:validity:{axis}"),
                "validity endpoints have different dimensions",
            ));
        }
        if lower.value > upper.value {
            return Err(CompileError::new(
                "reversed_validity",
                format!("claim:{claim_id}:validity:{axis}"),
                "normalized lower endpoint exceeds upper endpoint",
            ));
        }
        if let Some(previous) = global_dims.insert(manifest.axis(axis).to_owned(), lower.dims)
            && previous != lower.dims
        {
            return Err(CompileError::new(
                "validity_axis_dims_mismatch",
                format!("validity:{axis}"),
                "one named validity axis has contradictory dimensions across claims",
            ));
        }
        domain = if manifest.axis_kinds.contains_key(axis) {
            let quantity = kind.map_or_else(
                || QuantitySpec::dimensional(lower.dims),
                |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)),
            );
            domain.with_quantity(manifest.axis(axis), quantity, lower.value, upper.value)
        } else {
            domain.with(manifest.axis(axis), lower.value, upper.value)
        };
        receipts.push((axis.clone(), lower, upper));
    }
    Ok((domain, receipts))
}

fn compile_uncertainty(
    claim_id: &str,
    raw: &RawUncertainty,
    value_dims: Dims,
    kind: Option<QuantityKind>,
) -> Result<(UncertaintyModel, Option<ParsedQuantity>), CompileError> {
    let subject = format!("claim:{claim_id}:uncertainty");
    match &raw.kind {
        RawUncertaintyKind::Unstated => Ok((UncertaintyModel::Unstated, None)),
        RawUncertaintyKind::Absolute {
            number,
            unit,
            confidence,
            confidence_unit,
        } => {
            let quantity = if matches!(kind, Some(QuantityKind::Frequency(_))) {
                parse_declared_quantity(number, unit, kind, &subject)?
            } else {
                parse_linear_quantity(number, unit, &subject)?
            };
            if quantity.dims != value_dims {
                return Err(CompileError::new(
                    "uncertainty_dims_mismatch",
                    &subject,
                    "absolute uncertainty dimensions disagree with the claim value",
                ));
            }
            Ok((
                UncertaintyModel::HalfWidth {
                    half_width: quantity.value,
                    confidence: parse_confidence(confidence, confidence_unit, &subject)?,
                },
                Some(quantity),
            ))
        }
        RawUncertaintyKind::Relative {
            number,
            unit,
            confidence,
            confidence_unit,
        } => {
            let quantity = parse_linear_quantity(number, unit, &subject)?;
            if quantity.dims != Dims::NONE {
                return Err(CompileError::new(
                    "uncertainty_dims_mismatch",
                    &subject,
                    "relative uncertainty must use an explicit dimensionless unit",
                ));
            }
            Ok((
                UncertaintyModel::RelativeHalfWidth {
                    fraction: quantity.value,
                    confidence: parse_confidence(confidence, confidence_unit, &subject)?,
                },
                Some(quantity),
            ))
        }
    }
}

fn statistic_member(claim: ClaimId, component: LocalComponent) -> StatisticMember {
    match component {
        LocalComponent::Scalar => StatisticMember::scalar(claim),
        LocalComponent::CurveAbscissa(knot) => StatisticMember::curve_abscissa(claim, knot),
        LocalComponent::CurveOrdinate(knot) => StatisticMember::curve_ordinate(claim, knot),
    }
}

fn quantity_receipt(
    target: NormalizationTarget,
    quantity: &ParsedQuantity,
    frame: Option<&RawFrame>,
) -> NormalizationReceipt {
    NormalizationReceipt::new(
        target,
        quantity.literal_hash,
        quantity.dims,
        quantity.scale,
        quantity.offset,
        quantity.unit.clone(),
        MATDB_PACK_TARGET_BASIS,
        frame.map(|value| value.source.clone()),
        frame.map(|value| value.target.clone()),
    )
}

fn normalization_subject(claim: &str, component: LocalComponent) -> String {
    match component {
        LocalComponent::Scalar => format!("normalization:{claim}:scalar"),
        LocalComponent::CurveAbscissa(knot) => {
            format!("normalization:{claim}:curve:x:{knot}")
        }
        LocalComponent::CurveOrdinate(knot) => {
            format!("normalization:{claim}:curve:y:{knot}")
        }
    }
}

fn compile_joint(
    raw: &RawJoint,
    observation: ObservationId,
    components: &BTreeMap<(String, LocalComponent), ComponentMeta>,
) -> Result<(JointStatistics, Vec<NormalizationReceipt>, Vec<Decision>), CompileError> {
    let subject = format!("joint:{}:{}", raw.observation, raw.block_id);
    let size = raw.members.len();
    let packed_len = size
        .checked_mul(size + 1)
        .and_then(|value| value.checked_div(2))
        .ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                &subject,
                "packed matrix length overflowed",
            )
        })?;
    if raw.covariance.len() != packed_len {
        return Err(CompileError::new(
            "invalid_covariance_shape",
            &subject,
            format!(
                "{} members require {packed_len} packed covariance entries, found {}",
                size,
                raw.covariance.len()
            ),
        ));
    }
    if let Some(correlation) = &raw.correlation
        && correlation.len() != packed_len
    {
        return Err(CompileError::new(
            "invalid_correlation_shape",
            &subject,
            format!(
                "{} members require {packed_len} packed correlation entries, found {}",
                size,
                correlation.len()
            ),
        ));
    }

    let mut resolved = Vec::with_capacity(size);
    for (old_index, reference) in raw.members.iter().enumerate() {
        let (claim, component) = parse_member_reference(reference, &subject)?;
        let meta = components
            .get(&(claim.clone(), component))
            .cloned()
            .ok_or_else(|| {
                CompileError::new(
                    "unknown_joint_member",
                    &subject,
                    format!("member {reference:?} does not name a compiled claim component"),
                )
            })?;
        if !meta.observations.contains(&observation) {
            return Err(CompileError::new(
                "joint_observation_mismatch",
                &subject,
                format!("member {reference:?} is not backed by the joint block observation"),
            ));
        }
        resolved.push((old_index, meta));
    }
    resolved.sort_by_key(|(_, meta)| meta.member);
    if resolved
        .windows(2)
        .any(|pair| pair[0].1.member == pair[1].1.member)
    {
        return Err(CompileError::new(
            "duplicate_joint_member",
            &subject,
            "two source member references resolve to the same claim component",
        ));
    }

    let raw_covariance = raw
        .covariance
        .iter()
        .map(|token| parse_source_f64(token, "covariance entry", &subject))
        .collect::<Result<Vec<_>, _>>()?;
    let raw_correlation = raw
        .correlation
        .as_ref()
        .map(|values| {
            values
                .iter()
                .map(|token| parse_source_f64(token, "correlation entry", &subject))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let mut members = Vec::with_capacity(size);
    let mut covariance = Vec::with_capacity(packed_len);
    let mut correlation = raw_correlation
        .as_ref()
        .map(|_| Vec::with_capacity(packed_len));
    let mut receipts = Vec::with_capacity(packed_len);
    let mut decisions = Vec::with_capacity(packed_len);
    for (_, meta) in &resolved {
        members.push(meta.member);
    }
    for row in 0..size {
        for column in 0..=row {
            let (old_row, row_meta) = &resolved[row];
            let (old_column, column_meta) = &resolved[column];
            let old_index = packed_index(*old_row, *old_column);
            let raw_value = raw_covariance[old_index];
            let scale = row_meta.scale * column_meta.scale;
            let value = raw_value * scale;
            if !scale.is_finite() || scale <= 0.0 || !value.is_finite() {
                return Err(CompileError::new(
                    "non_finite_covariance_normalization",
                    &subject,
                    "covariance scaling produced a non-finite or non-positive transform",
                ));
            }
            if raw_value != 0.0 && value == 0.0 {
                return Err(CompileError::new(
                    "normalization_underflow",
                    &subject,
                    "nonzero covariance entry underflows the normalized representation",
                ));
            }
            let dims = row_meta
                .dims
                .checked_plus(column_meta.dims)
                .ok_or_else(|| {
                    CompileError::new(
                        "covariance_dimension_overflow",
                        &subject,
                        "member dimension product exceeds the six-base exponent range",
                    )
                })?;
            let basis = format!("({})*({})", row_meta.unit, column_meta.unit);
            let literal_hash = source_literal_hash(&raw.covariance[old_index], &basis);
            let row = u32::try_from(row).map_err(|_| {
                CompileError::new(
                    "resource_limit",
                    &subject,
                    "joint covariance row exceeds the runtime wire index",
                )
            })?;
            let column = u32::try_from(column).map_err(|_| {
                CompileError::new(
                    "resource_limit",
                    &subject,
                    "joint covariance column exceeds the runtime wire index",
                )
            })?;
            covariance.push(positive_zero(value));
            receipts.push(NormalizationReceipt::new(
                NormalizationTarget::JointCovariance {
                    observation,
                    block_id: raw.block_id.clone(),
                    row,
                    column,
                },
                literal_hash,
                dims,
                scale,
                0.0,
                basis,
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ));
            decisions.push(Decision::admit(
                format!("normalization:{subject}:covariance:{row}:{column}"),
                "covariance_normalized",
                "packed covariance entry scaled by both canonically permuted member transforms",
                Some(literal_hash),
            ));
            if let (Some(source), Some(output)) = (&raw_correlation, &mut correlation) {
                output.push(source[old_index]);
            }
        }
    }
    Ok((
        JointStatistics::new(
            observation,
            raw.block_id.clone(),
            members,
            covariance,
            correlation,
        ),
        receipts,
        decisions,
    ))
}

fn parse_member_reference(
    raw: &str,
    subject: &str,
) -> Result<(String, LocalComponent), CompileError> {
    let fields: Vec<&str> = raw.split(':').collect();
    match fields.as_slice() {
        [claim, "scalar"] => Ok((
            require_identifier(claim, "joint claim id", subject)?,
            LocalComponent::Scalar,
        )),
        [claim, "x", knot] | [claim, "y", knot] => {
            let claim = require_identifier(claim, "joint claim id", subject)?;
            let knot = knot.parse::<u32>().map_err(|_| {
                CompileError::new(
                    "invalid_joint_member",
                    subject,
                    "curve member knot must be a canonical u32",
                )
            })?;
            if knot.to_string() != fields[2] {
                return Err(CompileError::new(
                    "invalid_joint_member",
                    subject,
                    "curve member knot must use canonical unsigned decimal spelling",
                ));
            }
            let component = if fields[1] == "x" {
                LocalComponent::CurveAbscissa(knot)
            } else {
                LocalComponent::CurveOrdinate(knot)
            };
            Ok((claim, component))
        }
        _ => Err(CompileError::new(
            "invalid_joint_member",
            subject,
            "member must be claim:scalar, claim:x:knot, or claim:y:knot",
        )),
    }
}

fn packed_index(row: usize, column: usize) -> usize {
    let (row, column) = if row >= column {
        (row, column)
    } else {
        (column, row)
    };
    row * (row + 1) / 2 + column
}

fn sort_decisions(decisions: &mut [Decision]) {
    decisions.sort_by(|left, right| {
        (
            left.subject.as_str(),
            left.reason_code,
            left.verdict,
            left.source_hash,
            left.detail.as_str(),
        )
            .cmp(&(
                right.subject.as_str(),
                right.reason_code,
                right.verdict,
                right.source_hash,
                right.detail.as_str(),
            ))
    });
}

/// Render deterministic admission evidence for the selected compiler semantics.
///
/// The compiler and case identifiers bind provenance; they are not signatures,
/// source authentication, legal determinations, or proof that a named binary ran.
fn render_decision(compiler_id: &str, decision: &Decision) -> String {
    let case_id = hash_domain(
        DECISION_ID_DOMAIN,
        &decision.canonical_preimage(compiler_id),
    )
    .to_hex();
    let source_hash = decision
        .source_hash
        .map_or_else(String::new, |hash| hash.to_hex());
    let pack_hash = decision
        .pack_hash
        .map_or_else(String::new, |hash| hash.to_hex());
    format!(
        "{{\"check\":\"matdb-pack\",\"compiler\":\"{}\",\"case_id\":\"{}\",\"subject\":\"{}\",\"verdict\":\"{}\",\
         \"reason_code\":\"{}\",\"source_hash\":\"{}\",\"pack_hash\":\"{}\",\"detail\":\"{}\"}}",
        compiler_id,
        case_id,
        super::json_escape(&decision.subject),
        decision.verdict,
        decision.reason_code,
        source_hash,
        pack_hash,
        super::json_escape(&decision.detail),
    )
}

#[allow(clippy::too_many_lines)] // Publication barriers are intentionally explicit and sequential.
fn write_new_verified(
    output: &Path,
    bytes: &[u8],
    pack: &CompiledPack,
) -> Result<(), CompileError> {
    let expected = pack.content_hash();
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let file_name = output.file_name().ok_or_else(|| {
        CompileError::new(
            "invalid_output_path",
            "output",
            "output must name a file within an existing directory",
        )
    })?;
    let metadata = fs::symlink_metadata(parent).map_err(|error| {
        CompileError::new(
            "output_unavailable",
            "output",
            format!("cannot inspect output directory: {error}"),
        )
    })?;
    if !metadata.file_type().is_dir() || metadata.file_type().is_symlink() {
        return Err(CompileError::new(
            "invalid_output_directory",
            "output",
            "output parent must be an existing non-symlink directory",
        ));
    }
    let canonical_parent = parent.canonicalize().map_err(|error| {
        CompileError::new(
            "output_unavailable",
            "output",
            format!("cannot resolve output directory: {error}"),
        )
    })?;
    let canonical_output = canonical_parent.join(file_name);
    match fs::symlink_metadata(&canonical_output) {
        Ok(_) => {
            return Err(CompileError::new(
                "output_exists",
                "output",
                "refusing to replace an existing file, directory, or symlink",
            ));
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(CompileError::new(
                "output_unavailable",
                "output",
                format!("cannot inspect output path: {error}"),
            ));
        }
    }

    // The verified candidate deliberately remains as a same-inode retained
    // evidence link. This obeys the repository's no-deletion rule while
    // allowing hard_link to publish the final name atomically and no-clobber.
    let mut staged = None;
    for sequence in 0..1_024u16 {
        let candidate = canonical_parent.join(format!(
            ".frankensim-matdb-{}-{}-{sequence}.verified",
            expected,
            std::process::id()
        ));
        if candidate == canonical_output {
            continue;
        }
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                staged = Some((candidate, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => {
                return Err(CompileError::new(
                    "output_stage_failed",
                    "output",
                    format!("cannot create same-directory candidate: {error}"),
                ));
            }
        }
    }
    let (stage_path, mut file) = staged.ok_or_else(|| {
        CompileError::new(
            "output_stage_failed",
            "output",
            "exhausted the bounded same-directory candidate namespace",
        )
    })?;
    file.write_all(bytes).map_err(|error| {
        CompileError::new(
            "output_stage_failed",
            "output",
            format!("cannot write complete candidate: {error}"),
        )
    })?;
    file.sync_all().map_err(|error| {
        CompileError::new(
            "output_stage_failed",
            "output",
            format!("cannot sync complete candidate: {error}"),
        )
    })?;
    drop(file);
    let staged_bytes = read_bounded_regular(&stage_path, bytes.len(), "output-stage")?;
    if staged_bytes != bytes {
        return Err(CompileError::new(
            "output_validation_failed",
            "output",
            "candidate bytes differ from the verified compiler artifact",
        ));
    }
    pack.verify_bytes(expected, &staged_bytes)
        .map_err(|error| {
            CompileError::new(
                "output_validation_failed",
                "output",
                format!("candidate failed runtime admission: {error}"),
            )
        })?;
    fs::hard_link(&stage_path, &canonical_output).map_err(|error| {
        let (code, detail) = if error.kind() == std::io::ErrorKind::AlreadyExists {
            (
                "output_exists",
                "refusing to replace an output created during compilation".to_string(),
            )
        } else {
            (
                "output_publish_failed",
                format!("atomic no-clobber publication failed: {error}"),
            )
        };
        CompileError::new(code, "output", detail)
    })?;
    let published = read_bounded_regular(&canonical_output, bytes.len(), "output")?;
    if published != bytes {
        return Err(CompileError::new(
            "output_validation_failed",
            "output",
            "published hard link does not retain the verified candidate bytes",
        ));
    }
    pack.verify_bytes(expected, &published).map_err(|error| {
        CompileError::new(
            "output_validation_failed",
            "output",
            format!("published output failed runtime admission: {error}"),
        )
    })?;
    File::open(&canonical_parent)
        .and_then(|directory| directory.sync_all())
        .map_err(|error| {
            CompileError::new(
                "output_sync_failed",
                "output",
                format!("cannot sync output directory: {error}"),
            )
        })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

    const MATDB_PACK_FIXTURE_V1_BYTES_GOLDEN: usize = 3_177;
    const MATDB_PACK_FIXTURE_V1_CONTENT_HASH_GOLDEN: &str =
        "c1fb2f443708d297423179f4ac6024ee26b1d0c940a229d1d9084726ccbd2bc5";

    const SOURCE: &str = concat!(
        "frankensim.matdb-source.v1\n",
        "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
        "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant\n",
        "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1\n",
        "validity\tdensity\ttemperature\t0\t100\tdegC\n",
        "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear\n",
        "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
        "validity\tmodulus\ttemperature\t0\t100\tdegC\n",
        "frame\tmodulus\tspecimen\tlab\n",
        "joint\tcoupon\tdensity-modulus\tdensity:scalar,modulus:y:0\t0.000025,0,0.000009\t1,0,1\t1\n",
    );

    const NASA9_SOURCE: &str = concat!(
        "frankensim.nasa9-source.v1\n",
        "region\tN2\tlow\t200\t1000\tK\t100\tkPa\n",
        "coefficient\tN2\tlow\ta0\t0\tK^2\n",
        "coefficient\tN2\tlow\ta1\t0\tK\n",
        "coefficient\tN2\tlow\ta2\t3.5\t1\n",
        "coefficient\tN2\tlow\ta3\t0.001\tK^-1\n",
        "coefficient\tN2\tlow\ta4\t0\tK^-2\n",
        "coefficient\tN2\tlow\ta5\t0\tK^-3\n",
        "coefficient\tN2\tlow\ta6\t0\tK^-4\n",
        "coefficient\tN2\tlow\ta7\t100\tK\n",
        "coefficient\tN2\tlow\ta8\t1\t1\n",
        "region\tN2\thigh\t1000\t6000\tK\t100000\tPa\n",
        "coefficient\tN2\thigh\ta0\t0\tK^2\n",
        "coefficient\tN2\thigh\ta1\t0\tK\n",
        "coefficient\tN2\thigh\ta2\t4\t1\n",
        "coefficient\tN2\thigh\ta3\t0.0001\tK^-1\n",
        "coefficient\tN2\thigh\ta4\t0\tK^-2\n",
        "coefficient\tN2\thigh\ta5\t0\tK^-3\n",
        "coefficient\tN2\thigh\ta6\t0\tK^-4\n",
        "coefficient\tN2\thigh\ta7\t200\tK\n",
        "coefficient\tN2\thigh\ta8\t2\t1\n",
    );

    const KINETICS_SOURCE: &str = concat!(
        "frankensim.kinetics-source.v1\n",
        "reaction\twater-formation\tfirst-order\t300\t2500\tK\n",
        "parameter\twater-formation\tactivation_temperature\t12000\tK\n",
        "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1\n",
    );

    const SPECIES_SOURCE: &str = concat!(
        "frankensim.species-source.v1\n",
        "species\tN2\t28.0134\tg/mol\tgas\tideal-gas\t100\tkPa\tNASA-TP-2002-211556\n",
    );

    fn fixture_dir() -> PathBuf {
        loop {
            let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "frankensim-matdb-pack-test-{}-{sequence}",
                std::process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => return path,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(error) => panic!("unique fixture directory: {error}"),
            }
        }
    }

    fn manifest(profile: &str, redistribution: bool) -> String {
        let redistribution = if redistribution {
            "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n"
        } else {
            ""
        };
        format!(
            "{MANIFEST_HEADER}\n\
             pack_id\tfixture-alloy-x\n\
             {redistribution}\
             citation\tfixture handbook table 7\n\
             license\tCC-BY-4.0\n\
             source\tprimary\tsource.tsv\t{profile}\n"
        )
    }

    fn nasa9_manifest(pack_id: &str) -> String {
        format!(
            "{MANIFEST_HEADER}\n\
             pack_id\t{pack_id}\n\
             redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n\
             citation\tfixture NASA-9 species table\n\
             license\tCC-BY-4.0\n\
             source\tprimary\tsource.tsv\t{NASA9_PROFILE}\n"
        )
    }

    fn kinetics_manifest(pack_id: &str) -> String {
        format!(
            "{MANIFEST_HEADER}\n\
             pack_id\t{pack_id}\n\
             redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n\
             citation\tfixture first-order kinetics table\n\
             license\tCC-BY-4.0\n\
             source\tprimary\tsource.tsv\t{KINETICS_PROFILE}\n"
        )
    }

    fn species_manifest(pack_id: &str) -> String {
        format!(
            "{MANIFEST_HEADER}\n\
             pack_id\t{pack_id}\n\
             redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n\
             citation\tfixture licensed species metadata\n\
             license\tCC-BY-4.0\n\
             source\tprimary\tsource.tsv\t{SPECIES_PROFILE}\n"
        )
    }

    fn write_fixture(manifest_text: &str, source_text: &str) -> (PathBuf, PathBuf) {
        let directory = fixture_dir();
        let manifest_path = directory.join("manifest.tsv");
        let source_path = directory.join("source.tsv");
        fs::write(&manifest_path, manifest_text).expect("manifest fixture");
        fs::write(&source_path, source_text).expect("source fixture");
        (manifest_path, directory)
    }

    fn mapped_manifest(records: &str) -> String {
        manifest(MATERIAL_PROFILE, true).replacen(MANIFEST_HEADER, MAPPED_MANIFEST_HEADER, 1)
            + records
    }

    #[test]
    fn g3_manifest_v4_imports_hardness_test_history_without_conflating_loads() {
        let source = concat!(
            "frankensim.matdb-source.v1\n",
            "observation\ttest\tsynthetic specimen\tsynthetic hardness apparatus\tno empirical claim\n",
            "scalar\th1\ttest\thardness\t210\t1\tconstant\n",
            "uncertainty\th1\tunstated\t-\t-\t-\t-\n",
            "scalar\th2\ttest\thardness\t220\t1\tconstant\n",
            "uncertainty\th2\tunstated\t-\t-\t-\t-\n",
        );
        let manifest = mapped_manifest(concat!(
            "property\thardness\thardness\thardness-vickers\n",
            "hardness-test\th1\ttest\tsynthetic diamond pyramid geometry A\tsynthetic protocol revision 1\n",
            "hardness-hold\th1\t2\tN\t5\ts\n",
            "hardness-hold\th1\t20\tN\t10\ts\n",
            "hardness-test\th2\ttest\tsynthetic diamond pyramid geometry A\tsynthetic protocol revision 1\n",
            "hardness-hold\th2\t40\tN\t10\ts\n",
        )).replacen(MAPPED_MANIFEST_HEADER, HARDNESS_MANIFEST_HEADER, 1);
        let (path, _) = write_fixture(&manifest, source);
        let output = compile_manifest(&path).unwrap();
        let pack = NormalizedPack::from_bytes(&output.bytes).unwrap();
        assert_eq!(pack.schema_version(), 4);
        assert_eq!(pack.compiler(), HARDNESS_COMPILER_ID);
        let claim = pack
            .claims()
            .claims_for("hardness")
            .into_iter()
            .find(|(_, claim)| claim.key.hardness_test().unwrap().loading().len() == 2)
            .unwrap()
            .1;
        let key = claim.key.clone();
        let context = key.hardness_test().unwrap();
        assert_eq!(context.indenter(), "synthetic diamond pyramid geometry A");
        assert_eq!(context.protocol(), "synthetic protocol revision 1");
        assert_eq!(context.loading()[0].force_newtons(), 2.0);
        assert_eq!(context.loading()[0].dwell_seconds(), 5.0);
        assert_eq!(context.loading()[1].force_newtons(), 20.0);
        assert_eq!(context.loading()[1].dwell_seconds(), 10.0);
        assert_eq!(
            pack.claims()
                .observation(context.observation())
                .unwrap()
                .specimen,
            "synthetic specimen"
        );
        let point = fs_matdb::QueryPoint::new();
        let answer = pack
            .claims()
            .query_typed(&key, &point, fs_matdb::SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert_eq!(answer.evidence.value.value, 210.0);
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        let equivalent = manifest.replace("20\tN\t10\ts", "0.02\tkN\t10000\tms");
        let (path, _) = write_fixture(&equivalent, source);
        let converted = compile_manifest(&path).unwrap();
        let converted = NormalizedPack::from_bytes(&converted.bytes).unwrap();
        let equivalent_answer = converted
            .claims()
            .query_typed(&key, &point, fs_matdb::SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert_eq!(equivalent_answer.receipt.selected, answer.receipt.selected);
        assert_eq!(equivalent_answer.evidence.value.value, 210.0);
        assert_ne!(
            converted.content_hash(),
            pack.content_hash(),
            "source-spelled units remain bound in the source envelope"
        );
        for (changed, code) in [
            (
                manifest.replace(HARDNESS_MANIFEST_HEADER, TYPED_AXIS_MANIFEST_HEADER),
                "unknown_manifest_record",
            ),
            (
                manifest
                    .replace("hardness-test\th1", "hardness-test\tabsent")
                    .replace("hardness-hold\th1", "hardness-hold\tabsent"),
                "unused_hardness_test",
            ),
            (
                manifest.replace("20\tN\t10\ts", "20\tkg\t10\ts"),
                "invalid_hardness_test",
            ),
            (
                manifest.replace("hardness-hold\th2\t40\tN\t10\ts\n", ""),
                "invalid_hardness_test",
            ),
            (
                manifest.replace("hardness-test\th1\ttest", "hardness-test\th1\tunknown"),
                "unknown_observation",
            ),
            (
                manifest.clone() + "hardness-test\th1\ttest\tduplicate\tduplicate\n",
                "duplicate_hardness_test",
            ),
        ] {
            let (path, _) = write_fixture(&changed, source);
            assert_eq!(compile_manifest(&path).unwrap_err().code, code);
        }
    }

    #[test]
    fn g3_manifest_v3_normalizes_typed_temperature_and_frequency_axes() {
        let source = concat!(
            "frankensim.matdb-source.v1\n",
            "observation\ttest\tsynthetic\ttyped axis adapter fixture\tno empirical claim\n",
            "curve\tk\ttest\tresponse\tf\tHz\t1\t2:10,4:20\tlinear\n",
            "uncertainty\tk\tunstated\t-\t-\t-\t-\n",
            "validity\tk\tf\t2\t4\tHz\n",
            "curve\tt\ttest\tthermal_response\ttemperature\tdegC\t1\t10:1,20:3\tlinear\n",
            "uncertainty\tt\tunstated\t-\t-\t-\t-\n",
            "validity\tt\ttemperature\t10\t20\tdegC\n",
        );
        let manifest = mapped_manifest(
            "axis\tf\tfrequency\tangular-frequency\naxis\ttemperature\tT\tabsolute-temperature\n",
        )
        .replacen(MAPPED_MANIFEST_HEADER, TYPED_AXIS_MANIFEST_HEADER, 1);
        let (path, _) = write_fixture(&manifest, source);
        let output = compile_manifest(&path).unwrap();
        let pack = NormalizedPack::from_bytes(&output.bytes).unwrap();
        assert_eq!(pack.schema_version(), 3);
        assert_eq!(pack.compiler(), TYPED_AXIS_COMPILER_ID);
        let spec = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let angular = spec(QuantityKind::Frequency(FrequencyConvention::Angular));
        let temperature = spec(QuantityKind::AbsoluteTemperature);
        let tau = core::f64::consts::TAU;
        let point = fs_matdb::QueryPoint::new()
            .with_quantity("frequency", angular, 3.0 * tau)
            .unwrap()
            .with_quantity("T", temperature, 288.15)
            .unwrap();
        for (name, expected) in [("response", 15.0), ("thermal_response", 2.0)] {
            let answer = pack
                .claims()
                .query(name, &point, fs_matdb::SelectionPolicy::SingleClaimOnly)
                .unwrap();
            assert!((answer.evidence.value.value - expected).abs() < 1e-12);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let response = pack.claims().claims_for("response")[0];
        assert_eq!(
            response.1.validity.bound("frequency"),
            Some((2.0 * tau, 4.0 * tau))
        );
        for receipt in pack.normalizations().iter().filter(|r| matches!(r.target(), NormalizationTarget::ValidityBound { claim, .. } if *claim == response.0)) {
            assert_eq!(receipt.scale(), tau);
            assert_eq!(receipt.offset(), 0.0);
        }
        let alternate = source
            .replace("f\tHz\t1\t2:10,4:20", "f\tkHz\t1\t0.002:10,0.004:20")
            .replace("f\t2\t4\tHz", "f\t0.002\t0.004\tkHz");
        let (alternate_path, _) = write_fixture(&manifest, &alternate);
        let alternate_output = compile_manifest(&alternate_path).unwrap();
        let alternate_pack = alternate_output.pack.material();
        assert_ne!(alternate_pack.source_artifact(), pack.source_artifact());
        let alternate_answer = alternate_pack
            .claims()
            .query(
                "response",
                &point,
                fs_matdb::SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!((alternate_answer.evidence.value.value - 15.0).abs() < 1e-12);
        alternate_pack
            .claims()
            .verify_receipt(&alternate_answer.receipt)
            .unwrap();
        let wrong = fs_matdb::QueryPoint::new()
            .with_quantity(
                "frequency",
                spec(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
                3.0 * tau,
            )
            .unwrap();
        assert!(matches!(
            pack.claims().query(
                "response",
                &wrong,
                fs_matdb::SelectionPolicy::SingleClaimOnly
            ),
            Err(fs_matdb::MatDbError::AxisQuantityMismatch { .. })
        ));
        let interval_manifest = manifest.replace("absolute-temperature", "temperature-difference");
        let (path, _) = write_fixture(&interval_manifest, source);
        let interval_output = compile_manifest(&path).unwrap();
        assert_ne!(
            interval_output.pack.material().source_artifact(),
            pack.source_artifact()
        );
        let interval_pack = interval_output.pack.material();
        assert_eq!(
            interval_pack.claims().claims_for("thermal_response")[0]
                .1
                .validity
                .bound("T"),
            Some((10.0, 20.0))
        );
        let point = fs_matdb::QueryPoint::new()
            .with_quantity("T", spec(QuantityKind::TemperatureDifference), 15.0)
            .unwrap();
        assert_eq!(
            interval_pack
                .claims()
                .query(
                    "thermal_response",
                    &point,
                    fs_matdb::SelectionPolicy::SingleClaimOnly
                )
                .unwrap()
                .evidence
                .value
                .value,
            2.0
        );
        // The v2 parser refuses the extra declaration field; it cannot silently
        // reinterpret typed axes as its legacy name-only mapping.
        let (path, _) = write_fixture(
            &manifest.replace(TYPED_AXIS_MANIFEST_HEADER, MAPPED_MANIFEST_HEADER),
            source,
        );
        assert!(compile_manifest(&path).is_err());
    }

    #[test]
    fn g3_manifest_v2_renames_curves_and_receipts_without_reinterpreting_v1() {
        let (old_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let old = compile_manifest(&old_path).unwrap();
        let (path, _) = write_fixture(
            &mapped_manifest(
                "axis\ttemperature\tT\nproperty\tdensity\tbulk_density\tdimensional\n",
            ),
            SOURCE,
        );
        let output = compile_manifest(&path).unwrap();
        let pack = output.pack.material();
        let point = fs_matdb::QueryPoint::new().with("T", 323.15).unwrap();
        let answer = pack
            .claims()
            .query(
                "young_modulus",
                &point,
                fs_matdb::SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        let old_answer = old
            .pack
            .material()
            .claims()
            .query(
                "young_modulus",
                &fs_matdb::QueryPoint::new()
                    .with("temperature", 323.15)
                    .unwrap(),
                fs_matdb::SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(answer.evidence.value, old_answer.evidence.value);
        assert_ne!(pack.source_artifact(), old.pack.source_artifact());
        assert_ne!(answer.receipt.selected, old_answer.receipt.selected);
        assert_eq!(pack.compiler(), MAPPED_COMPILER_ID);
        assert!(pack.claims().claims_for("density").is_empty());
        assert_eq!(pack.claims().claims_for("bulk_density").len(), 1);
        let receipts: Vec<_> = pack.normalizations().iter().filter(|receipt| matches!(receipt.target(), NormalizationTarget::ValidityBound { axis, .. } if axis == "T")).collect();
        assert_eq!(receipts.len(), 4);
        assert!(
            receipts
                .iter()
                .all(|receipt| receipt.source_basis() == "degC" && receipt.offset() == 273.15)
        );
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(old.bytes.len(), MATDB_PACK_FIXTURE_V1_BYTES_GOLDEN);
        assert_eq!(
            old.pack.content_hash().to_string(),
            MATDB_PACK_FIXTURE_V1_CONTENT_HASH_GOLDEN
        );
        for records in [
            "axis\tmissing\tT\n",
            "property\tyoung_modulus\tdensity\tdimensional\n",
            "axis\ttemperature\tT\naxis\ttemperature\tother\n",
        ] {
            let (path, _) = write_fixture(&mapped_manifest(records), SOURCE);
            assert!(
                compile_manifest(&path).is_err(),
                "mapping must refuse: {records}"
            );
        }
        let v1_with_mapping = manifest(MATERIAL_PROFILE, true) + "axis\ttemperature\tT\n";
        assert_eq!(
            parse_manifest(&v1_with_mapping).unwrap_err().code,
            "unknown_manifest_record"
        );
    }

    #[test]
    fn g0_manifest_v2_distinguishes_absolute_temperature_and_interval_width() {
        let source = concat!(
            "frankensim.matdb-source.v1\n",
            "observation\ttest\tsynthetic\ttemperature convention test\tno empirical claim\n",
            "scalar\tabs\ttest\tabs\t10\tdegC\tconstant\n",
            "uncertainty\tabs\tabsolute\t2\tdegC\t0.95\t1\n",
            "scalar\tdelta\ttest\tdelta\t10\tdegC\tconstant\n",
            "uncertainty\tdelta\tunstated\t-\t-\t-\t-\n",
        );
        let (path, _) = write_fixture(
            &mapped_manifest(
                "property\tabs\tabs\tabsolute-temperature\nproperty\tdelta\tdelta\ttemperature-difference\n",
            ),
            source,
        );
        let output = compile_manifest(&path).unwrap();
        let pack = output.pack.material();
        for (name, kind, expected) in [
            ("abs", QuantityKind::AbsoluteTemperature, 283.15),
            ("delta", QuantityKind::TemperatureDifference, 10.0),
        ] {
            let key = PropertyKey::with_quantity(
                name,
                QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)),
            );
            let answer = pack
                .claims()
                .query_typed(
                    &key,
                    &fs_matdb::QueryPoint::new(),
                    fs_matdb::SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert_eq!(answer.evidence.value.value, expected);
            if name == "abs" {
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::HalfWidth {
                        half_width: 2.0,
                        confidence: 0.95
                    }
                );
            }
        }
    }

    #[test]
    fn g3_manifest_v2_frequency_conversion_scales_values_uncertainty_and_covariance() {
        let source = concat!(
            "frankensim.matdb-source.v1\n",
            "observation\ttest\tsynthetic\tfrequency convention test\tno empirical claim\n",
            "scalar\tf\ttest\tf\t2\tHz\tconstant\n",
            "uncertainty\tf\tabsolute\t0.1\tHz\t0.95\t1\n",
            "scalar\tg\ttest\tg\t3\trad/s\tconstant\n",
            "uncertainty\tg\tabsolute\t0.2\trad/s\t0.95\t1\n",
            "joint\ttest\tfg\tf:scalar,g:scalar\t0.01,0.002,0.04\t-\t-\n",
        );
        let (path, _) = write_fixture(
            &mapped_manifest(
                "property\tf\tf\tangular-frequency\nproperty\tg\tg\tangular-frequency\n",
            ),
            source,
        );
        let output = compile_manifest(&path).unwrap();
        let pack = NormalizedPack::from_bytes(&output.bytes).unwrap();
        let spec = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Angular),
            ValueForm::Static,
        ));
        let answer = pack
            .query_joint_typed(
                &[
                    PropertyKey::with_quantity("f", spec),
                    PropertyKey::with_quantity("g", spec),
                ],
                &fs_matdb::QueryPoint::new(),
                fs_matdb::SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        let tau = core::f64::consts::TAU;
        assert!((answer.members[0].evidence.value.value - 2.0 * tau).abs() < 1e-14);
        assert_eq!(answer.members[1].evidence.value.value, 3.0);
        assert!(
            matches!(answer.members[0].evidence.value.uncertainty, UncertaintyModel::HalfWidth { half_width, .. } if (half_width - 0.1 * tau).abs() < 1e-14)
        );
        let fs_matdb::JointCorrelation::Covariance { covariance, .. } = &answer.receipt.correlation
        else {
            panic!("declared joint covariance must survive conversion");
        };
        for (found, expected) in covariance.iter().zip([0.01 * tau * tau, 0.002 * tau, 0.04]) {
            assert!((found - expected).abs() < 1e-14);
        }
        pack.verify_joint_receipt(&answer.receipt).unwrap();
        let uncertainty = pack.normalizations().iter().find(|receipt| matches!(receipt.target(), NormalizationTarget::ClaimUncertainty { claim } if *claim == answer.receipt.selected[0])).unwrap();
        assert_eq!(uncertainty.source_basis(), "Hz");
        assert_eq!(uncertainty.scale(), tau);
        // An equivalent source-unit expression yields the same physical
        // answer but retains a different source-literal provenance.
        let alternate = source
            .replace("\tf\t2\tHz\t", "\tf\t0.002\tkHz\t")
            .replace("\tf\tabsolute\t0.1\tHz\t", "\tf\tabsolute\t0.0001\tkHz\t")
            .replace("\t0.01,0.002,0.04\t", "\t0.00000001,0.000002,0.04\t");
        let (other_path, _) = write_fixture(
            &mapped_manifest(
                "property\tf\tf\tangular-frequency\nproperty\tg\tg\tangular-frequency\n",
            ),
            &alternate,
        );
        let other = compile_manifest(&other_path).unwrap();
        let other_answer = other
            .pack
            .material()
            .claims()
            .query_typed(
                &PropertyKey::with_quantity("f", spec),
                &fs_matdb::QueryPoint::new(),
                fs_matdb::SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!(
            (other_answer.evidence.value.value - answer.members[0].evidence.value.value).abs()
                < 1e-14
        );
        assert_ne!(other.pack.source_artifact(), pack.source_artifact());
        let (malformed_path, _) = write_fixture(
            &mapped_manifest("property\tf\tf\tangular-frequency\n"),
            &source.replace("\t-\t-\n", "\t-\t1\n"),
        );
        let refusal = compile_manifest(&malformed_path).unwrap_err();
        assert_eq!(refusal.code, "invalid_correlation_unit");
        assert_eq!(refusal.compiler_id, MAPPED_COMPILER_ID);
    }

    #[test]
    fn g1_nasa_seed_compiles_into_the_actual_thermoelastic_state_resolver() {
        use fs_material::state_point::{
            MaterialPropertySelection, resolve_isotropic_thermoelastic_state_point,
        };
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../data/matdb/seed-v1/aluminum-2024-t3-nasa-tn-d6448/manifest.tsv");
        let compiled = compile_manifest(&path).unwrap();
        let pack = NormalizedPack::from_bytes(&compiled.bytes).unwrap();
        let card = fs_matdb::MaterialCard::assemble(
            MaterialStateId {
                chemistry: "aluminum-2024".to_owned(),
                phase: "solid".to_owned(),
                process: "T3; NASA TN D-6448".to_owned(),
                revision: 0,
            },
            pack.claims().clone(),
            Vec::new(),
        )
        .unwrap();
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = fs_matdb::QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap();
        let state = resolve_isotropic_thermoelastic_state_point(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let law = state.law();
        assert_eq!(
            (law.rho, law.e, law.cp, law.conductivity, law.t0),
            (2705.0, 72_400_000_000.0, 840.0, 126.0, 300.0)
        );
        assert!((law.alpha_t - 23e-6).abs() < 1e-18);
        assert_eq!(state.poisson_ratio(), 0.33);
        assert_eq!(state.resolved().properties().len(), 6);
        for property in state.resolved().properties() {
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
        assert!(
            resolve_isotropic_thermoelastic_state_point(
                &card,
                &fs_matdb::QueryPoint::new()
                    .with_quantity("T", temperature, 301.0)
                    .unwrap(),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        let omega = 1000.0;
        let thickness: f64 = 0.002;
        let relaxation_time =
            thickness.powi(2) * 2705.0 * 840.0 / (core::f64::consts::PI.powi(2) * 126.0);
        let x = omega * relaxation_time;
        let expected = 72.4e9 * (23e-6_f64).powi(2) * 300.0 / (2705.0 * 840.0) * x / (1.0 + x * x);
        assert!((law.loss_factor(omega, thickness) / expected - 1.0).abs() < 1e-12);
    }

    #[test]
    fn compiles_scalar_curve_covariance_and_is_two_run_deterministic() {
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let first = compile_manifest(&manifest_path).expect("fixture compiles");
        let second = compile_manifest(&manifest_path).expect("independent repeat compiles");

        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        assert_eq!(first.pack.content_hash(), second.pack.content_hash());
        assert_eq!(
            first.bytes.len(),
            MATDB_PACK_FIXTURE_V1_BYTES_GOLDEN,
            "normalized compiler fixture byte length moved; actual content hash {}",
            first.pack.content_hash()
        );
        assert_eq!(
            first.pack.content_hash().to_string(),
            MATDB_PACK_FIXTURE_V1_CONTENT_HASH_GOLDEN,
            "normalized compiler fixture content hash moved"
        );
        let decoded = NormalizedPack::from_bytes_verified(first.pack.content_hash(), &first.bytes)
            .expect("compiler output is runtime-loadable");
        assert_eq!(decoded.to_bytes(), first.bytes);
        assert_eq!(decoded.claims().claim_count(), 2);
        assert_eq!(decoded.joint_statistics().len(), 1);
        let (modulus_id, modulus) = decoded.claims().claims_for("young_modulus")[0];
        let PropertyValue::Curve { knots, .. } = &modulus.value else {
            panic!("young_modulus remains a normalized curve");
        };
        assert_eq!(knots.as_slice(), &[(273.15, 210.0e9), (373.15, 202.0e9)]);
        assert_eq!(
            modulus.validity.bound("temperature"),
            Some((273.15, 373.15))
        );

        let temperature_value_receipts: Vec<_> = decoded
            .normalizations()
            .iter()
            .filter(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ClaimValue(member)
                        if member.claim() == modulus_id && receipt.source_basis() == "degC"
                )
            })
            .collect();
        assert_eq!(temperature_value_receipts.len(), 2);
        for receipt in temperature_value_receipts {
            assert_eq!(receipt.scale().to_bits(), 1.0f64.to_bits());
            assert_eq!(receipt.offset().to_bits(), 273.15f64.to_bits());
            assert_eq!(receipt.source_frame(), Some("specimen"));
            assert_eq!(receipt.target_frame(), Some("lab"));
        }

        let temperature_validity_receipts: Vec<_> = decoded
            .normalizations()
            .iter()
            .filter(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ValidityBound { claim, axis, .. }
                        if *claim == modulus_id && axis == "temperature"
                )
            })
            .collect();
        assert_eq!(temperature_validity_receipts.len(), 2);
        for receipt in temperature_validity_receipts {
            assert_eq!(receipt.source_basis(), "degC");
            assert_eq!(receipt.scale().to_bits(), 1.0f64.to_bits());
            assert_eq!(receipt.offset().to_bits(), 273.15f64.to_bits());
            assert_eq!(receipt.source_frame(), None);
            assert_eq!(receipt.target_frame(), None);
        }

        let modulus_uncertainty_receipt = decoded
            .normalizations()
            .iter()
            .find(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ClaimUncertainty { claim } if *claim == modulus_id
                )
            })
            .expect("modulus uncertainty receipt");
        assert_eq!(
            modulus_uncertainty_receipt.offset().to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(modulus_uncertainty_receipt.source_frame(), Some("specimen"));
        assert_eq!(modulus_uncertainty_receipt.target_frame(), Some("lab"));
        assert!(decoded.normalizations().iter().all(|receipt| {
            !matches!(
                receipt.target(),
                NormalizationTarget::JointCovariance { .. }
            ) || (receipt.offset() == 0.0
                && receipt.source_frame().is_none()
                && receipt.target_frame().is_none())
        }));

        let joint = &decoded.joint_statistics()[0];
        assert_eq!(joint.members().len(), 2);
        assert_eq!(
            joint.correlation().expect("source correlation")[1].to_bits(),
            0.0f64.to_bits()
        );
        for (row, member) in joint.members().iter().copied().enumerate() {
            let claim = decoded
                .claims()
                .claim(member.claim())
                .expect("joint member claim");
            let diagonal = joint.covariance()[packed_index(row, row)];
            let receipt = decoded
                .normalizations()
                .iter()
                .find(|receipt| {
                    matches!(
                        receipt.target(),
                        NormalizationTarget::JointCovariance {
                            row: target_row,
                            column,
                            ..
                        } if usize::try_from(*target_row).ok() == Some(row)
                            && *column == *target_row
                    )
                })
                .expect("diagonal covariance receipt");
            match claim.key.name() {
                "density" => {
                    assert!((diagonal - 25.0).abs() <= 25.0 * f64::EPSILON);
                    assert_eq!(receipt.dims(), Dims([-6, 2, 0, 0, 0, 0]));
                    // The covariance transform composes two independently rounded
                    // fs-qty unit scales, so allow both conversion roundoff steps.
                    assert!(
                        (receipt.scale() - 1.0e6).abs() <= 4.0e6 * f64::EPSILON,
                        "density covariance receipt scale mismatch at row {row}: {receipt:?}"
                    );
                }
                "young_modulus" => {
                    assert!((diagonal - 9.0e12).abs() <= 9.0e12 * f64::EPSILON);
                    assert_eq!(receipt.dims(), Dims([-2, 2, -4, 0, 0, 0]));
                    assert!(
                        (receipt.scale() - 1.0e18).abs() <= 4.0e18 * f64::EPSILON,
                        "modulus covariance receipt scale mismatch at row {row}: {receipt:?}"
                    );
                }
                other => panic!("unexpected covariance member property {other}"),
            }
        }
        assert_eq!(decoded.normalizations().len(), 14);
        assert!(first.decisions.iter().all(|decision| {
            decision.verdict == "admit" && decision.pack_hash == Some(first.pack.content_hash())
        }));
    }

    #[test]
    fn compiles_complete_nasa9_regions_into_a_verified_model_pack() {
        let (manifest_path, _) = write_fixture(&nasa9_manifest("N2"), NASA9_SOURCE);
        let first = compile_manifest(&manifest_path).expect("NASA-9 fixture compiles");
        let second = compile_manifest(&manifest_path).expect("independent repeat compiles");

        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        assert_eq!(first.pack.content_hash(), second.pack.content_hash());
        let pack = first.pack.model();
        assert_eq!(pack.pack_id(), "N2");
        assert_eq!(pack.compiler(), NASA9_COMPILER_ID);
        assert_eq!(pack.models().len(), 2);
        assert_eq!(pack.normalizations().len(), 24);
        assert!(pack.models().iter().all(|card| {
            card.law.0 == NASA9_LAW_ID
                && card.law_version == NASA9_LAW_VERSION
                && card.parameters.len() == 10
                && card.validity.bound("T").is_some()
        }));
        let decoded = NormalizedModelPack::from_bytes_verified(pack.content_hash(), &first.bytes)
            .expect("compiler output is runtime-loadable");
        assert_eq!(decoded.to_bytes(), first.bytes);
        assert!(first.decisions.iter().all(|decision| {
            decision.verdict == "admit" && decision.pack_hash == Some(pack.content_hash())
        }));
    }

    #[test]
    fn nasa9_profile_refuses_incomplete_dims_mismatched_and_overlapping_regions() {
        let cases = [
            (
                NASA9_SOURCE.replacen("coefficient\tN2\tlow\ta8\t1\t1\n", "", 1),
                "missing_nasa9_coefficient",
            ),
            (
                NASA9_SOURCE.replacen(
                    "coefficient\tN2\tlow\ta0\t0\tK^2",
                    "coefficient\tN2\tlow\ta0\t0\t1",
                    1,
                ),
                "nasa9_dims_mismatch",
            ),
            (
                NASA9_SOURCE.replacen(
                    "coefficient\tN2\tlow\ta1\t0\tK",
                    "coefficient\tN2\tlow\ta1\t0\tdegC",
                    1,
                ),
                "affine_nasa9_parameter_unit",
            ),
            (
                NASA9_SOURCE.replacen(
                    "region\tN2\thigh\t1000\t6000\tK",
                    "region\tN2\thigh\t900\t6000\tK",
                    1,
                ),
                "overlapping_nasa9_regions",
            ),
        ];
        for (source, code) in cases {
            let (manifest_path, _) = write_fixture(&nasa9_manifest("N2"), &source);
            let error = compile_manifest(&manifest_path).expect_err("malformed NASA-9 refuses");
            assert_eq!(error.code, code);
        }

        let (manifest_path, _) = write_fixture(&nasa9_manifest("O2"), NASA9_SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("pack id binds species");
        assert_eq!(error.code, "species_pack_id_mismatch");
    }

    #[test]
    fn compiles_first_order_kinetics_into_a_verified_model_pack() {
        let (manifest_path, _) =
            write_fixture(&kinetics_manifest("water-formation"), KINETICS_SOURCE);
        let first = compile_manifest(&manifest_path).expect("kinetics fixture compiles");
        let second = compile_manifest(&manifest_path).expect("independent repeat compiles");

        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        assert_eq!(first.pack.content_hash(), second.pack.content_hash());
        let pack = first.pack.model();
        assert_eq!(pack.pack_id(), "water-formation");
        assert_eq!(pack.compiler(), KINETICS_COMPILER_ID);
        assert_eq!(pack.models().len(), 1);
        assert_eq!(pack.normalizations().len(), 4);
        let card = &pack.models()[0];
        assert_eq!(card.law.0, KINETICS_LAW_ID);
        assert_eq!(card.law_version, KINETICS_LAW_VERSION);
        assert_eq!(card.parameters.len(), KINETICS_PARAMETER_NAMES.len());
        assert_eq!(
            card.parameters["activation_temperature"].dims,
            KINETICS_TEMPERATURE_DIMS
        );
        assert_eq!(
            card.parameters["pre_exponential"].dims,
            KINETICS_FIRST_ORDER_RATE_DIMS
        );
        assert!(card.validity.bound("T").is_some());
        assert!(
            card.provenance
                .source
                .contains("[reaction:water-formation]")
        );
        assert!(card.provenance.source.contains("[rate-basis:first-order]"));
        let decoded = NormalizedModelPack::from_bytes_verified(pack.content_hash(), &first.bytes)
            .expect("compiler output is runtime-loadable");
        assert_eq!(decoded.to_bytes(), first.bytes);
        assert!(first.decisions.iter().all(|decision| {
            decision.verdict == "admit" && decision.pack_hash == Some(pack.content_hash())
        }));
    }

    #[test]
    fn kinetics_profile_refuses_incomplete_ambiguous_or_mismatched_records() {
        let cases = [
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1\n",
                    "",
                    1,
                ),
                "missing_kinetics_parameter",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1",
                    "parameter\twater-formation\tpre_exponential\t2.5e7\tK",
                    1,
                ),
                "kinetics_dims_mismatch",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tactivation_temperature\t12000\tK",
                    "parameter\twater-formation\tactivation_temperature\t12000\tdegC",
                    1,
                ),
                "affine_kinetics_parameter_unit",
            ),
            (
                KINETICS_SOURCE.replacen("\tfirst-order\t", "\tsecond-order\t", 1),
                "unsupported_kinetics_rate_basis",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tactivation_temperature\t12000\tK",
                    "parameter\twater-formation\tactivation_temperature\t-1\tK",
                    1,
                ),
                "invalid_activation_temperature",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1",
                    "parameter\twater-formation\tpre_exponential\t0\ts^-1",
                    1,
                ),
                "invalid_pre_exponential",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "reaction\twater-formation\tfirst-order\t300\t2500\tK",
                    "reaction\twater-formation\tfirst-order\t2500\t300\tK",
                    1,
                ),
                "invalid_kinetics_temperature_range",
            ),
            (
                KINETICS_SOURCE.replacen(
                    "parameter\twater-formation\tpre_exponential",
                    "parameter\tother-reaction\tpre_exponential",
                    1,
                ),
                "unknown_kinetics_reaction",
            ),
        ];
        for (source, code) in cases {
            let (manifest_path, _) = write_fixture(&kinetics_manifest("water-formation"), &source);
            let error = compile_manifest(&manifest_path).expect_err("malformed kinetics refuses");
            assert_eq!(error.code, code);
        }

        let (manifest_path, _) = write_fixture(&kinetics_manifest("other"), KINETICS_SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("pack id binds reaction");
        assert_eq!(error.code, "reaction_pack_id_mismatch");
    }

    #[test]
    fn multi_observation_claims_are_canonicalized_by_content_id() {
        let source = SOURCE
            .replacen(
                "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
                concat!(
                    "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
                    "observation\tsecond\talloy-X-retest\tASTM-fixture\trepeat coupon series\n",
                ),
                1,
            )
            .replacen(
                "scalar\tdensity\tcoupon\tdensity",
                "scalar\tdensity\tcoupon,second\tdensity",
                1,
            );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let compiled = compile_manifest(&manifest_path).expect("two observations compile");
        let density = compiled.pack.material().claims().claims_for("density")[0].1;

        assert_eq!(density.observations.len(), 2);
        assert!(
            density
                .observations
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
    }

    #[test]
    fn absent_redistribution_terms_refuse_before_source_loading() {
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, false), SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("terms are load-bearing");

        assert_eq!(error.code, "missing_redistribution_terms");
        assert_eq!(error.subject, "manifest");
    }

    #[test]
    fn missing_license_refuses_before_source_loading() {
        let manifest = manifest(MATERIAL_PROFILE, true).replacen("license\tCC-BY-4.0\n", "", 1);
        let (manifest_path, _) = write_fixture(&manifest, SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("license is load-bearing");

        assert_eq!(error.code, "missing_license");
        assert_eq!(error.subject, "manifest");
    }

    #[test]
    fn blank_unit_refuses_instead_of_inference() {
        let source = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\t\tconstant",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("unit inference is forbidden");

        assert_eq!(error.code, "ambiguous_unit");
        assert!(error.subject.contains("source:primary:line:"));
    }

    #[test]
    fn numeric_continuation_cannot_masquerade_as_a_unit() {
        let source = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\te3\tconstant",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("e3 is not a unit token");

        assert_eq!(error.code, "invalid_unit");
    }

    #[test]
    fn confidence_requires_the_unconverted_dimensionless_basis() {
        let source = SOURCE.replacen(
            "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1",
            "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t%",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("confidence has no receipt target");

        assert_eq!(error.code, "invalid_confidence_unit");
    }

    #[test]
    fn unevaluable_value_interpolation_pairs_refuse() {
        let scalar = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tlinear",
            1,
        );
        let curve = SOURCE.replacen(
            "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear",
            "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tconstant",
            1,
        );
        for source in [scalar, curve] {
            let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
            let error =
                compile_manifest(&manifest_path).expect_err("claim could never be evaluated");
            assert_eq!(error.code, "incompatible_interpolation");
        }
    }

    #[test]
    fn decimal_underflow_refuses_before_signed_zero_can_be_laundered() {
        for number in ["1e-999", "-1e-999"] {
            let source = SOURCE.replacen(
                "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
                &format!("scalar\tdensity\tcoupon\tdensity\t{number}\t1\tconstant"),
                1,
            );
            let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
            let error =
                compile_manifest(&manifest_path).expect_err("nonzero significand underflows");
            assert_eq!(error.code, "numeric_underflow");
        }
    }

    #[test]
    fn compiles_species_association_into_a_verified_species_pack() {
        let (manifest_path, _) = write_fixture(&species_manifest("N2"), SPECIES_SOURCE);
        let first = compile_manifest(&manifest_path).expect("species fixture compiles");
        let second = compile_manifest(&manifest_path).expect("independent repeat compiles");

        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        assert_eq!(first.pack.content_hash(), second.pack.content_hash());
        let pack = first.pack.species();
        assert_eq!(pack.pack_id(), "N2");
        assert_eq!(pack.compiler(), SPECIES_COMPILER_ID);
        assert_eq!(pack.association().species().as_str(), "N2");
        assert_eq!(
            pack.association().molar_mass().to_bits(),
            0.028_013_4f64.to_bits()
        );
        assert_eq!(pack.association().standard_state_phase(), "gas");
        assert_eq!(pack.association().reference_eos(), "ideal-gas");
        assert_eq!(
            pack.association().reference_pressure().to_bits(),
            100_000.0f64.to_bits()
        );
        assert_eq!(
            pack.association().elemental_reference(),
            "NASA-TP-2002-211556"
        );
        assert_eq!(pack.association().sources().len(), 2);
        assert_eq!(pack.normalizations().len(), 2);
        assert_eq!(
            pack.normalizations()[0].target(),
            SpeciesNormalizationTarget::MolarMass
        );
        assert_eq!(pack.normalizations()[0].dims(), SPECIES_MOLAR_MASS_DIMS);
        assert_eq!(pack.normalizations()[0].source_basis(), "g/mol");
        assert_eq!(
            pack.normalizations()[1].target(),
            SpeciesNormalizationTarget::ReferencePressure
        );
        assert_eq!(
            pack.normalizations()[1].dims(),
            SPECIES_REFERENCE_PRESSURE_DIMS
        );
        assert_eq!(pack.normalizations()[1].source_basis(), "kPa");

        let decoded = NormalizedSpeciesPack::from_bytes_verified(pack.content_hash(), &first.bytes)
            .expect("verified species pack decodes");
        assert_eq!(&decoded, pack);
        assert_eq!(decoded.to_bytes(), first.bytes);
        assert!(first.decisions.iter().all(|decision| {
            decision.verdict == "admit" && decision.pack_hash == Some(pack.content_hash())
        }));
    }

    #[test]
    fn unknown_source_profile_remains_a_typed_refusal() {
        let (manifest_path, _) = write_fixture(&manifest("unknown-v1", true), SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("unknown profile refuses");

        assert_eq!(error.code, "unsupported_source_profile");
        assert_eq!(error.subject, "source:primary");
    }

    #[test]
    fn malformed_species_associations_refuse_without_inference() {
        let cases = [
            (
                SPECIES_SOURCE.replacen("28.0134\tg/mol", "28.0134\tkg", 1),
                "species_molar_mass_dims_mismatch",
            ),
            (
                SPECIES_SOURCE.replacen("\tgas\t", "\tliquid\t", 1),
                "unsupported_species_phase",
            ),
            (
                SPECIES_SOURCE.replacen("\tideal-gas\t", "\tcubic\t", 1),
                "unsupported_species_reference_eos",
            ),
            (
                SPECIES_SOURCE.replacen("28.0134\tg/mol", "0\tg/mol", 1),
                "invalid_species_molar_mass",
            ),
            (
                SPECIES_SOURCE.replacen("100\tkPa", "100\tdegC", 1),
                "affine_species_unit",
            ),
            (
                SPECIES_SOURCE.replacen("NASA-TP-2002-211556", "bad reference", 1),
                "invalid_elemental_reference",
            ),
            (
                SPECIES_SOURCE.replacen(
                    "species\tN2\t28.0134\tg/mol\tgas\tideal-gas\t100\tkPa\tNASA-TP-2002-211556\n",
                    "",
                    1,
                ),
                "missing_species_association",
            ),
            (
                format!(
                    "{SPECIES_SOURCE}species\tN2\t28.0134\tg/mol\tgas\tideal-gas\t100\tkPa\tNASA-TP-2002-211556\n"
                ),
                "duplicate_species_association",
            ),
        ];
        for (source, code) in cases {
            let (manifest_path, _) = write_fixture(&species_manifest("N2"), &source);
            let error = compile_manifest(&manifest_path).expect_err("malformed species refuses");
            assert_eq!(error.code, code);
        }

        let (manifest_path, _) = write_fixture(&species_manifest("O2"), SPECIES_SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("pack id binds species");
        assert_eq!(error.code, "species_pack_id_mismatch");
    }

    #[test]
    fn covariance_member_permutation_preserves_the_matrix_semantics() {
        let density = ClaimId(ContentHash([1; 32]));
        let modulus = ClaimId(ContentHash([2; 32]));
        let observation = ObservationId(ContentHash([3; 32]));
        let mut components = BTreeMap::new();
        components.insert(
            ("density".to_string(), LocalComponent::Scalar),
            ComponentMeta {
                member: StatisticMember::scalar(density),
                dims: Dims([-3, 1, 0, 0, 0, 0]),
                scale: 1_000.0,
                unit: "g/cm3".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        components.insert(
            ("modulus".to_string(), LocalComponent::CurveAbscissa(0)),
            ComponentMeta {
                member: StatisticMember::curve_abscissa(modulus, 0),
                dims: Dims([0, 0, 0, 1, 0, 0]),
                scale: 1.0,
                unit: "degC".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        components.insert(
            ("modulus".to_string(), LocalComponent::CurveOrdinate(0)),
            ComponentMeta {
                member: StatisticMember::curve_ordinate(modulus, 0),
                dims: Dims([-1, 1, -2, 0, 0, 0]),
                scale: 1.0e9,
                unit: "GPa".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        let first = RawJoint {
            observation: "coupon".to_string(),
            block_id: "block".to_string(),
            members: vec![
                "density:scalar".to_string(),
                "modulus:x:0".to_string(),
                "modulus:y:0".to_string(),
            ],
            covariance: vec![
                "4".to_string(),
                "1".to_string(),
                "9".to_string(),
                "2".to_string(),
                "3".to_string(),
                "16".to_string(),
            ],
            correlation: Some(vec![
                "1".to_string(),
                "0.1".to_string(),
                "1".to_string(),
                "0.2".to_string(),
                "0.3".to_string(),
                "1".to_string(),
            ]),
            source_hash: ContentHash([4; 32]),
        };
        let second = RawJoint {
            observation: "coupon".to_string(),
            block_id: "block".to_string(),
            members: vec![
                "modulus:y:0".to_string(),
                "density:scalar".to_string(),
                "modulus:x:0".to_string(),
            ],
            covariance: vec![
                "16".to_string(),
                "2".to_string(),
                "4".to_string(),
                "3".to_string(),
                "1".to_string(),
                "9".to_string(),
            ],
            correlation: Some(vec![
                "1".to_string(),
                "0.2".to_string(),
                "1".to_string(),
                "0.3".to_string(),
                "0.1".to_string(),
                "1".to_string(),
            ]),
            source_hash: ContentHash([4; 32]),
        };
        let (first, _, _) =
            compile_joint(&first, observation, &components).expect("first member order");
        let (second, _, _) =
            compile_joint(&second, observation, &components).expect("permuted member order");

        assert_eq!(first.members(), second.members());
        assert_eq!(first.covariance(), second.covariance());
        assert_eq!(first.correlation(), second.correlation());
    }

    #[test]
    fn publication_never_replaces_an_existing_path() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join("pack.fsmatpk");
        fs::write(&output, b"sentinel").expect("preexisting output");

        let error = write_new_verified(&output, &compiled.bytes, &compiled.pack)
            .expect_err("replacement is forbidden");
        assert_eq!(error.code, "output_exists");
        assert_eq!(fs::read(output).expect("sentinel retained"), b"sentinel");
    }

    #[test]
    fn publication_exposes_only_a_complete_runtime_verified_hard_link() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join("pack.fsmatpk");

        write_new_verified(&output, &compiled.bytes, &compiled.pack)
            .expect("verified no-clobber publication");
        let retained = read_bounded_regular(&output, compiled.bytes.len(), "test-output")
            .expect("published output");
        assert_eq!(retained, compiled.bytes);
        let decoded = NormalizedPack::from_bytes_verified(compiled.pack.content_hash(), &retained)
            .expect("published runtime artifact");
        assert_eq!(&decoded, compiled.pack.material());
    }

    #[test]
    fn caller_cannot_alias_the_hidden_candidate_namespace() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join(format!(
            ".frankensim-matdb-{}-{}-0.verified",
            compiled.pack.content_hash(),
            std::process::id()
        ));

        write_new_verified(&output, &compiled.bytes, &compiled.pack)
            .expect("candidate alias is skipped");
        assert_eq!(
            read_bounded_regular(&output, compiled.bytes.len(), "aliased-output")
                .expect("published aliased output"),
            compiled.bytes
        );
    }

    #[test]
    fn receipt_budget_refuses_before_normalized_allocations() {
        let mut raw = RawDatabase::default();
        raw.claims.insert(
            "oversized".to_string(),
            RawClaim {
                id: "oversized".to_string(),
                observations: vec!["observation".to_string()],
                property: "oversized".to_string(),
                value: RawClaimValue::Curve {
                    abscissa: "x".to_string(),
                    abscissa_unit: "1".to_string(),
                    ordinate_unit: "1".to_string(),
                    points: vec![("0".to_string(), "0".to_string()); 50_001],
                },
                interpolation: InterpolationPolicy::TabulatedOnly,
                source_id: "source".to_string(),
                source_hash: ContentHash([1; 32]),
                record_hash: ContentHash([2; 32]),
            },
        );

        let error =
            preflight_normalization_budget(&raw).expect_err("two receipts per knot exceed cap");
        assert_eq!(error.code, "resource_limit");
        assert_eq!(error.subject, "normalizations");
    }

    #[test]
    fn repeated_provenance_budget_refuses_before_claim_set_allocations() {
        let manifest = Manifest {
            version: 1,
            properties: BTreeMap::new(),
            axes: BTreeMap::new(),
            axis_kinds: BTreeMap::new(),
            hardness_tests: BTreeMap::new(),
            pack_id: "fixture".to_string(),
            redistribution_terms: "permitted".to_string(),
            citation: "x".repeat(1_048_000),
            license: "MIT".to_string(),
            sources: Vec::new(),
        };
        let mut raw = RawDatabase::default();
        for index in 0..64 {
            let id = format!("observation-{index}");
            raw.observations.insert(
                id.clone(),
                RawObservation {
                    id,
                    specimen: String::new(),
                    method: String::new(),
                    caveats: String::new(),
                    source_id: "source".to_string(),
                    source_hash: ContentHash([1; 32]),
                    record_hash: ContentHash([2; 32]),
                },
            );
        }
        preflight_provenance_budget(&manifest, &raw).expect("64 records remain under the cap");
        raw.observations.insert(
            "observation-64".to_string(),
            RawObservation {
                id: "observation-64".to_string(),
                specimen: String::new(),
                method: String::new(),
                caveats: String::new(),
                source_id: "source".to_string(),
                source_hash: ContentHash([1; 32]),
                record_hash: ContentHash([2; 32]),
            },
        );

        let error = preflight_provenance_budget(&manifest, &raw)
            .expect_err("65 repeated megabyte-scale citations exceed the cap");
        assert_eq!(error.code, "resource_limit");
        assert_eq!(error.subject, "provenance");
    }

    #[test]
    fn invalid_utf8_refusals_are_bound_to_the_offending_source_bytes() {
        let mut hashes = Vec::new();
        for invalid in [0xff_u8, 0xfe] {
            let directory = fixture_dir();
            let manifest_path = directory.join("manifest.tsv");
            let source_path = directory.join("source.tsv");
            fs::write(&manifest_path, manifest(MATERIAL_PROFILE, true)).expect("manifest fixture");
            let mut source = format!("{SOURCE_HEADER}\n").into_bytes();
            source.push(invalid);
            fs::write(&source_path, &source).expect("invalid UTF-8 source fixture");

            let error = compile_manifest(&manifest_path).expect_err("invalid UTF-8 refuses");
            assert_eq!(error.code, "invalid_source_encoding");
            let source_hash = error.input_hash.expect("offending byte hash");
            assert_eq!(source_hash, hash_domain(SOURCE_FILE_DOMAIN, &source));
            hashes.push(source_hash);
        }

        assert_ne!(hashes[0], hashes[1]);
    }

    #[test]
    fn refusal_evidence_is_bound_to_the_rejected_input_envelope() {
        let first_source = SOURCE.replacen("\t7.85\tg/cm3\t", "\t7.85\t\t", 1);
        let second_source = SOURCE.replacen("\t7.85\tg/cm3\t", "\t8.00\t\t", 1);
        let (first_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &first_source);
        let (second_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &second_source);
        let first_error = compile_manifest(&first_path).expect_err("blank unit refuses");
        let second_error = compile_manifest(&second_path).expect_err("blank unit refuses");
        let first_hash = first_error.input_hash.expect("first envelope hash");
        let second_hash = second_error.input_hash.expect("second envelope hash");

        assert_ne!(first_hash, second_hash);
        let first = Decision::refusal(&first_error, Some(first_hash));
        let second = Decision::refusal(&second_error, Some(second_hash));
        assert_ne!(
            first.canonical_preimage(COMPILER_ID),
            second.canonical_preimage(COMPILER_ID)
        );
        assert!(render_decision(COMPILER_ID, &first).contains(&first_hash.to_hex()));
    }

    #[test]
    fn later_refusal_retains_prior_admission_decisions_deterministically() {
        let rejected_source = SOURCE.replacen(
            "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
            "uncertainty\tmodulus\tabsolute\t2\tkg\t0.95\t1\n",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &rejected_source);
        let first = compile_manifest(&manifest_path).expect_err("later claim must refuse");
        let second = compile_manifest(&manifest_path).expect_err("repeat must refuse identically");

        assert_eq!(first, second);
        assert_eq!(first.code, "uncertainty_dims_mismatch");
        assert_eq!(first.subject, "claim:modulus:uncertainty");
        assert!(
            first
                .prior_decisions
                .iter()
                .all(|row| row.verdict == "admit")
        );
        assert!(
            first
                .prior_decisions
                .iter()
                .all(|row| row.pack_hash.is_none())
        );
        for (subject, reason_code) in [
            ("source:primary", "source_schema_admitted"),
            ("manifest:redistribution", "redistribution_permitted"),
            ("observation:coupon", "observation_normalized"),
            ("claim:density", "claim_normalized"),
        ] {
            assert!(
                first
                    .prior_decisions
                    .iter()
                    .any(|row| row.subject == subject && row.reason_code == reason_code),
                "missing prior admission {reason_code} for {subject}"
            );
        }

        let rows = refusal_decisions(&first, first.input_hash);
        let terminal: Vec<_> = rows.iter().filter(|row| row.verdict == "refuse").collect();
        assert_eq!(terminal.len(), 1);
        assert_eq!(terminal[0].reason_code, "uncertainty_dims_mismatch");
        assert_eq!(terminal[0].source_hash, first.input_hash);
    }

    #[test]
    fn decision_identity_distinguishes_source_and_pack_hash_roles() {
        let error = CompileError::new("fixture_refusal", "fixture", "typed refusal");
        let hash = ContentHash([9; 32]);
        let source = Decision::refusal(&error, Some(hash));
        let pack = Decision::refusal_for_pack(&error, hash);

        assert_ne!(
            source.canonical_preimage(COMPILER_ID),
            pack.canonical_preimage(COMPILER_ID)
        );
        assert_ne!(
            hash_domain(DECISION_ID_DOMAIN, &source.canonical_preimage(COMPILER_ID)),
            hash_domain(DECISION_ID_DOMAIN, &pack.canonical_preimage(COMPILER_ID))
        );
    }

    #[test]
    fn decision_identity_is_bound_to_the_selected_profile_compiler() {
        let decision = Decision::admit("fixture", "fixture_admitted", "typed admission", None);

        let material_preimage = decision.canonical_preimage(COMPILER_ID);
        let species_preimage = decision.canonical_preimage(SPECIES_COMPILER_ID);
        assert_ne!(material_preimage, species_preimage);
        assert_ne!(
            hash_domain(DECISION_ID_DOMAIN, &material_preimage),
            hash_domain(DECISION_ID_DOMAIN, &species_preimage)
        );
        assert_ne!(
            render_decision(COMPILER_ID, &decision),
            render_decision(SPECIES_COMPILER_ID, &decision)
        );
    }

    #[test]
    fn source_profiles_select_their_exact_compiler_identity() {
        for (profile, expected) in [
            (SourceProfile::Material, COMPILER_ID),
            (SourceProfile::Interface, INTERFACE_COMPILER_ID),
            (SourceProfile::Nasa9, NASA9_COMPILER_ID),
            (SourceProfile::Kinetics, KINETICS_COMPILER_ID),
            (SourceProfile::Species, SPECIES_COMPILER_ID),
        ] {
            assert_eq!(profile.compiler_id(), expected);
        }
    }

    #[test]
    fn decision_rows_are_stable_json_lines_without_host_paths() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let rows: Vec<String> = compiled
            .decisions
            .iter()
            .map(|decision| render_decision(compiled.pack.compiler_id(), decision))
            .collect();

        assert_eq!(
            rows,
            compiled
                .decisions
                .iter()
                .map(|decision| render_decision(compiled.pack.compiler_id(), decision))
                .collect::<Vec<_>>()
        );
        assert!(
            rows.iter()
                .all(|row| row.starts_with("{\"check\":\"matdb-pack\""))
        );
        assert!(
            rows.iter()
                .all(|row| !row.contains(&directory.display().to_string()))
        );
        assert!(rows.iter().all(|row| row.contains("\"pack_hash\":\"")));
    }
}
