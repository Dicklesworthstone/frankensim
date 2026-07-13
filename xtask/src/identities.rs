//! Generated semantic-identity registry and owner-declaration policy gate.

use super::depgraph::{JsonParser, JsonValue};
use super::{Violation, json_escape};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt::Write as _;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, ExitCode};

const DECLARATION_MARKER: &str = concat!("frankensim-identity-", "schema-v1");
const REGISTRY_FILE: &str = "identity-schemas.json";
const AUTHORITY_FILE: &str = "identity-authorities.json";
const AUTHORITY_SCHEMA: &str = "frankensim-identity-authorities-v1";
const GOLDEN_COUPLING_SCHEMA: &str = "frankensim-golden-couplings-v1";
const SCHEMA_FINGERPRINT_DOMAIN: &str =
    "org.frankensim.xtask.semantic-identity-schema-fingerprint.v1";
const SCHEMA_FINGERPRINT_ALGORITHM: &str = "blake3-256-domain-separated-v1";
const GATE_IMPLEMENTATION_PATH: &str = "xtask/src/identities.rs";
const MAX_SOURCE_BYTES: u64 = 4 * 1024 * 1024;
const MAX_DECLARATIONS: usize = 4_096;
const REQUIRED_IDENTITY_IDS: &[&str] = &[
    "ci:content-snapshot",
    "ci:quality-proof-record",
    "ci:x86-cross-verdict",
    "ci:x86-runtime-verdict",
    "fs-adjoint:explanation-node",
    "fs-adjoint:explanation-receipt",
    "fs-alloc:hugepage-decision",
    "fs-bisect:compound-family",
    "fs-checker:decision-report",
    "fs-checker:semantic-plugin",
    "fs-checker:semantic-registry",
    "fs-checker:semantic-report",
    "fs-exec:gemm-tune-key",
    "fs-exec:tilepool-placement",
    "fs-exec:tune-row",
    "fs-exec:tuning-decision",
    "fs-ir:planner-cache-key",
    "fs-la:depgraph-receipt",
    "fs-la:gemm-build-fingerprint",
    "fs-ledger:artifact-content",
    "fs-ledger:color-admission-policy",
    "fs-ledger:color-node",
    "fs-ledger:derived-color-waiver-subject",
    "fs-ledger:physical-instance",
    "fs-ledger:session-flush-batch",
    "fs-ledger:session-mutation-claim",
    "fs-ledger:session-terminal-events",
    "fs-ledger:source-color-waiver-subject",
    "fs-ledger:source-origin-request",
    "fs-ledger:vcs-commit-envelope",
    "fs-ledger:vcs-commit-leaf",
    "fs-ledger:vcs-commit-root",
    "fs-ledger:vcs-ledger-lineage",
    "fs-obs:event-content",
    "fs-obs:replay-identity-frame",
    "fs-package:claim-declaration",
    "fs-package:claim-verification-subject",
    "fs-package:coverage-decision",
    "fs-package:package-root",
    "fs-package:presence-decision",
    "fs-package:release-admission-context",
    "fs-package:semantic-witness",
    "fs-package:signature-subject",
    "fs-package:source-certificate-subject",
    "fs-package:verification-receipt",
    "fs-package:waiver-authorization-subject",
    "fs-plan:voi-audit-context",
    "fs-plan:voi-ranked-menu",
    "fs-plan:voi-ranked-source",
    "fs-rand:stream-checkpoint",
    "fs-rand:stream-position",
    "fs-recompute:artifact-content",
    "fs-recompute:node-record",
    "fs-roofline:baseline-record",
    "fs-roofline:executable-build",
    "fs-roofline:execution-binding",
    "fs-roofline:finalized-run",
    "fs-roofline:production-axes-receipt",
    "fs-session:durable-governor-id",
    "fs-session:gate-binding-id",
    "fs-session:gemm-execution-receipt",
    "fs-session:gemm-tune-row-receipt",
    "fs-session:meter-receipt",
    "fs-session:meter-report-id",
    "fs-session:pause-acknowledgement-id",
    "fs-session:pause-acknowledgement-receipt",
    "fs-session:pressure-action-id",
    "fs-session:pressure-receipt",
    "fs-session:resume-activation-id",
    "fs-session:resume-activation-receipt",
    "fs-session:retained-evidence",
    "fs-session:session-open-id",
    "fs-session:session-open-receipt",
    "fs-session:session-token-identity",
    "fs-session:submission-agent-key",
    "fs-session:submission-program",
    "fs-session:submission-receipt",
    "fs-session:submission-request-id",
    "fs-verify:fem1d-mms-class",
    "fs-verify:fem1d-mms-problem",
    "fs-vskeleton:artifact-content",
    "xtask:constellation-bootstrap-provenance",
    "xtask:constellation-lock",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FieldClass {
    Semantic,
    Derived,
    Nonsemantic,
}

impl FieldClass {
    fn parse(value: &str) -> Option<Self> {
        match value {
            "semantic" => Some(Self::Semantic),
            "derived" => Some(Self::Derived),
            "nonsemantic" => Some(Self::Nonsemantic),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Semantic => "semantic",
            Self::Derived => "derived",
            Self::Nonsemantic => "nonsemantic",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SourceField {
    qualified: String,
    class: FieldClass,
    reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Mutation {
    field: String,
    test_path: String,
    test_symbol: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SourceBinding {
    source_field: String,
    semantic_fields: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SchemaConstant {
    path: Option<String>,
    symbol: String,
}

impl SchemaConstant {
    fn canonical(&self) -> String {
        self.path.as_ref().map_or_else(
            || self.symbol.clone(),
            |path| format!("{path}#{}", self.symbol),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SchemaFunction {
    path: Option<String>,
    symbol: String,
}

impl SchemaFunction {
    fn canonical(&self) -> String {
        self.path.as_ref().map_or_else(
            || self.symbol.clone(),
            |path| format!("{path}#{}", self.symbol),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IdentityDecl {
    id: String,
    owner: String,
    version_const: String,
    version: u32,
    domain: String,
    domain_const: Option<String>,
    encoder: String,
    encoder_helpers: Vec<String>,
    schema_functions: Vec<SchemaFunction>,
    schema_constants: Vec<SchemaConstant>,
    schema_dependencies: Vec<String>,
    digest: String,
    encoding: String,
    sources: Vec<String>,
    source_fields: Vec<SourceField>,
    source_bindings: Vec<SourceBinding>,
    external_semantic_fields: Vec<String>,
    semantic_fields: Vec<String>,
    excluded_fields: Vec<(String, String)>,
    consumers: Vec<String>,
    mutations: Vec<Mutation>,
    nonsemantic_mutations: Vec<Mutation>,
    field_guard: String,
    transport_guard: String,
    version_guard: String,
    coupling_surface: String,
    schema_base_hash: Option<[u8; 32]>,
    schema_fingerprint: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CouplingSurface {
    Rust {
        file: String,
        version_const: String,
        version: u32,
        domain_const: Option<String>,
        domain: Option<String>,
        schema_fingerprint: Option<String>,
    },
    External {
        file: String,
        symbol: String,
        version: u32,
        domain: String,
        schema_fingerprint: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExternalOwner {
    id: String,
    path: String,
    symbol: String,
    version: u32,
    domain: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IdentityExemption {
    path: String,
    symbol: String,
    reason: String,
    covered_by: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AuthorityManifest {
    required_ids: BTreeSet<String>,
    external_owners: Vec<ExternalOwner>,
    exemptions: Vec<IdentityExemption>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct IdentityCandidate {
    path: String,
    symbol: String,
    identity_signal: String,
    sink_signal: String,
}

fn identity_violation(owner: &str, detail: impl Into<String>) -> Violation {
    Violation {
        check: "semantic-identities",
        crate_name: owner.to_string(),
        detail: detail.into(),
    }
}

fn source_files(root: &Path) -> (Vec<PathBuf>, Vec<Violation>) {
    fn visit(
        root: &Path,
        path: &Path,
        allow_missing: bool,
        files: &mut Vec<PathBuf>,
        violations: &mut Vec<Violation>,
    ) {
        let context = path
            .strip_prefix(root)
            .ok()
            .filter(|relative| !relative.as_os_str().is_empty())
            .map_or_else(
                || "<repo>".to_string(),
                |relative| relative.display().to_string(),
            );
        let metadata = match std::fs::symlink_metadata(path) {
            Ok(metadata) => metadata,
            Err(error) if allow_missing && error.kind() == std::io::ErrorKind::NotFound => return,
            Err(error) => {
                violations.push(identity_violation(
                    &context,
                    format!("identity source discovery cannot inspect path: {error}"),
                ));
                return;
            }
        };
        if metadata.file_type().is_symlink() {
            violations.push(identity_violation(
                &context,
                "identity source discovery refuses symlinked paths",
            ));
            return;
        }
        if !metadata.is_dir() {
            violations.push(identity_violation(
                &context,
                "identity source discovery root is not a directory",
            ));
            return;
        }
        let read_dir = match std::fs::read_dir(path) {
            Ok(read_dir) => read_dir,
            Err(error) => {
                violations.push(identity_violation(
                    &context,
                    format!("identity source discovery cannot read directory: {error}"),
                ));
                return;
            }
        };
        let mut entries = Vec::new();
        for entry in read_dir {
            match entry {
                Ok(entry) => entries.push(entry),
                Err(error) => violations.push(identity_violation(
                    &context,
                    format!("identity source discovery cannot read directory entry: {error}"),
                )),
            }
        }
        entries.sort_by_key(std::fs::DirEntry::file_name);
        for entry in entries {
            let path = entry.path();
            let relative = path
                .strip_prefix(root)
                .unwrap_or(&path)
                .display()
                .to_string();
            let file_type = match entry.file_type() {
                Ok(file_type) => file_type,
                Err(error) => {
                    violations.push(identity_violation(
                        &relative,
                        format!("identity source discovery cannot inspect file type: {error}"),
                    ));
                    continue;
                }
            };
            if file_type.is_symlink() {
                violations.push(identity_violation(
                    &relative,
                    "identity source discovery refuses symlinked paths",
                ));
            } else if file_type.is_dir() {
                if entry.file_name() != "target" {
                    visit(root, &path, false, files, violations);
                }
            } else if file_type.is_file() && path.extension().is_some_and(|ext| ext == "rs") {
                files.push(path);
            } else if !file_type.is_file() {
                violations.push(identity_violation(
                    &relative,
                    "identity source discovery refuses non-file paths",
                ));
            }
        }
    }

    let mut files = Vec::new();
    let mut violations = Vec::new();
    visit(
        root,
        &root.join("crates"),
        true,
        &mut files,
        &mut violations,
    );
    visit(root, &root.join("xtask"), true, &mut files, &mut violations);
    files.sort();
    (files, violations)
}

fn quoted_literal(line: &str) -> Option<&str> {
    let start = line.find('"')? + 1;
    let end = line.rfind('"')?;
    (end >= start).then_some(&line[start..end])
}

fn identity_marker_lines(text: &str) -> BTreeSet<usize> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted {
            escaped: bool,
            start: usize,
            line: usize,
        },
        Raw {
            hashes: usize,
        },
        LineComment,
        BlockComment {
            depth: usize,
        },
    }

    let bytes = text.as_bytes();
    let mut lines = BTreeSet::new();
    let mut state = State::Normal;
    let mut line = 0usize;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted {
                        escaped: false,
                        start: index + 1,
                        line,
                    };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if byte == b'\n' {
                    line += 1;
                }
                index += 1;
            }
            State::Quoted {
                escaped,
                start,
                line: start_line,
            } => {
                if escaped {
                    state = State::Quoted {
                        escaped: false,
                        start,
                        line: start_line,
                    };
                } else if byte == b'\\' {
                    state = State::Quoted {
                        escaped: true,
                        start,
                        line: start_line,
                    };
                } else if byte == b'"' {
                    if text.get(start..index) == Some(DECLARATION_MARKER) {
                        lines.insert(start_line);
                    }
                    state = State::Normal;
                } else {
                    if byte == b'\n' {
                        line += 1;
                    }
                    state = State::Quoted {
                        escaped: false,
                        start,
                        line: start_line,
                    };
                }
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'\n' {
                    line += 1;
                }
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    line += 1;
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth - 1;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    if byte == b'\n' {
                        line += 1;
                    }
                    index += 1;
                }
            }
        }
    }
    lines
}

fn list(value: &str) -> Vec<String> {
    if value == "none" {
        Vec::new()
    } else {
        value.split(',').map(str::to_string).collect()
    }
}

fn safe_relative(path: &str) -> bool {
    !path.is_empty()
        && !Path::new(path).is_absolute()
        && Path::new(path)
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
}

fn checked_repo_file(root: &Path, relative: &str, purpose: &str) -> Result<PathBuf, String> {
    if !safe_relative(relative) {
        return Err(format!(
            "{purpose} path {relative:?} must be a canonical repo-relative path"
        ));
    }
    let canonical_root = root
        .canonicalize()
        .map_err(|error| format!("repository root is not canonicalizable: {error}"))?;
    let mut target = root.to_path_buf();
    for component in Path::new(relative).components() {
        let Component::Normal(component) = component else {
            return Err(format!(
                "{purpose} path {relative:?} contains a forbidden component"
            ));
        };
        target.push(component);
        let metadata = std::fs::symlink_metadata(&target).map_err(|error| {
            format!("{purpose} path {relative:?} is missing or unreadable: {error}")
        })?;
        if metadata.file_type().is_symlink() {
            return Err(format!(
                "{purpose} path {relative:?} traverses symlink component {:?}",
                target.strip_prefix(root).unwrap_or(&target)
            ));
        }
    }
    let canonical_target = target
        .canonicalize()
        .map_err(|error| format!("{purpose} path {relative:?} is not canonicalizable: {error}"))?;
    if !canonical_target.starts_with(&canonical_root) {
        return Err(format!(
            "{purpose} path {relative:?} escapes the repository root"
        ));
    }
    let metadata = std::fs::metadata(&canonical_target)
        .map_err(|error| format!("{purpose} path {relative:?} is unreadable: {error}"))?;
    if !metadata.is_file() {
        return Err(format!("{purpose} path {relative:?} is not a regular file"));
    }
    if metadata.len() > MAX_SOURCE_BYTES {
        return Err(format!(
            "{purpose} path {relative:?} exceeds the {MAX_SOURCE_BYTES}-byte scan cap"
        ));
    }
    Ok(canonical_target)
}

fn read_repo_utf8(root: &Path, relative: &str, purpose: &str) -> Result<String, String> {
    let path = checked_repo_file(root, relative, purpose)?;
    std::fs::read_to_string(path)
        .map_err(|error| format!("{purpose} path {relative:?} is not readable UTF-8: {error}"))
}

fn canonical_symbol(symbol: &str) -> bool {
    let mut chars = symbol.chars();
    chars
        .next()
        .is_some_and(|first| first == '_' || first.is_ascii_alphabetic())
        && chars.all(|ch| ch == '_' || ch.is_ascii_alphanumeric())
}

fn canonical_function_reference(reference: &str) -> bool {
    !reference.is_empty() && reference.split("::").all(canonical_symbol)
}

fn parse_schema_functions(value: &str) -> Result<Vec<SchemaFunction>, String> {
    let mut functions = Vec::new();
    for item in list(value) {
        let function = if let Some((path, symbol)) = item.split_once('#') {
            if path.contains('#') || symbol.contains('#') || !safe_relative(path) {
                return Err(format!(
                    "schema function {item:?} must be repo-relative path#function-or-Type::method"
                ));
            }
            SchemaFunction {
                path: Some(path.to_string()),
                symbol: symbol.to_string(),
            }
        } else {
            SchemaFunction {
                path: None,
                symbol: item.clone(),
            }
        };
        if !canonical_function_reference(&function.symbol) {
            return Err(format!(
                "schema function {item:?} must name a canonical function or Type::method"
            ));
        }
        functions.push(function);
    }
    functions.sort();
    if functions.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("schema_functions contains duplicates".to_string());
    }
    Ok(functions)
}

fn parse_schema_constants(value: &str) -> Result<Vec<SchemaConstant>, String> {
    let mut constants = Vec::new();
    for item in list(value) {
        let constant = if let Some((path, symbol)) = item.split_once('#') {
            if path.contains('#') || symbol.contains('#') || !safe_relative(path) {
                return Err(format!(
                    "schema constant {item:?} must be repo-relative path#SYMBOL"
                ));
            }
            SchemaConstant {
                path: Some(path.to_string()),
                symbol: symbol.to_string(),
            }
        } else {
            SchemaConstant {
                path: None,
                symbol: item.clone(),
            }
        };
        if !canonical_symbol(&constant.symbol) {
            return Err(format!(
                "schema constant {item:?} must name a bare Rust identifier"
            ));
        }
        constants.push(constant);
    }
    constants.sort();
    if constants.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("schema_constants contains duplicates".to_string());
    }
    Ok(constants)
}

fn parse_schema_constant_reference(value: &str) -> Result<SchemaConstant, String> {
    let constants = parse_schema_constants(value)?;
    let [constant] = constants.as_slice() else {
        return Err(format!(
            "schema constant reference {value:?} must name exactly one constant"
        ));
    };
    Ok(constant.clone())
}

fn parse_schema_dependencies(value: &str, id: &str) -> Result<Vec<String>, String> {
    let mut dependencies = list(value);
    for dependency in &dependencies {
        if !dependency.contains(':') || dependency.chars().any(char::is_whitespace) {
            return Err(format!(
                "identity {id}: schema dependency {dependency:?} is not a canonical identity id"
            ));
        }
        if dependency == id {
            return Err(format!(
                "identity {id}: schema_dependencies may not contain a self-dependency"
            ));
        }
    }
    dependencies.sort();
    if dependencies.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(format!(
            "identity {id}: schema_dependencies contains duplicates"
        ));
    }
    Ok(dependencies)
}

fn domain_carries_version(domain: &str, version: u32) -> bool {
    let version = format!("v{version}");
    domain == "fsid"
        || domain
            .split(['.', ':', '-'])
            .any(|segment| segment == version.as_str())
}

fn parse_source_fields(value: &str) -> Result<Vec<SourceField>, String> {
    let mut fields = Vec::new();
    for item in list(value) {
        let mut parts = item.splitn(3, ':');
        let qualified = parts.next().unwrap_or_default();
        let class_name = parts.next().unwrap_or_default();
        let reason = parts.next().map(str::to_string);
        let Some(class) = FieldClass::parse(class_name) else {
            return Err(format!(
                "source field {qualified:?} has unknown class {class_name:?}"
            ));
        };
        if qualified.split_once('.').is_none() {
            return Err(format!(
                "source field {qualified:?} must be qualified as Type.field"
            ));
        }
        if class == FieldClass::Semantic && reason.is_some() {
            return Err(format!(
                "semantic source field {qualified:?} must not carry an exclusion reason"
            ));
        }
        if class != FieldClass::Semantic && reason.as_deref().is_none_or(str::is_empty) {
            return Err(format!(
                "{} source field {qualified:?} needs a non-empty reason",
                class.name()
            ));
        }
        fields.push(SourceField {
            qualified: qualified.to_string(),
            class,
            reason,
        });
    }
    Ok(fields)
}

fn parse_reasoned(value: &str) -> Result<Vec<(String, String)>, String> {
    let mut fields = Vec::new();
    for item in list(value) {
        let Some((field, reason)) = item.split_once(':') else {
            return Err(format!("excluded field {item:?} needs field:reason"));
        };
        if field.is_empty() || reason.is_empty() {
            return Err(format!(
                "excluded field {item:?} contains an empty component"
            ));
        }
        fields.push((field.to_string(), reason.to_string()));
    }
    Ok(fields)
}

fn parse_mutations(value: &str) -> Result<Vec<Mutation>, String> {
    let mut mutations = Vec::new();
    for item in list(value) {
        let Some((field, target)) = item.split_once(':') else {
            return Err(format!("mutation {item:?} needs field:path#test_symbol"));
        };
        let Some((test_path, test_symbol)) = target.split_once('#') else {
            return Err(format!("mutation {item:?} needs path#test_symbol"));
        };
        if field.is_empty()
            || test_symbol.is_empty()
            || !safe_relative(test_path)
            || !test_symbol
                .chars()
                .all(|ch| ch == '_' || ch.is_ascii_alphanumeric())
        {
            return Err(format!("mutation target {item:?} is not canonical"));
        }
        mutations.push(Mutation {
            field: field.to_string(),
            test_path: test_path.to_string(),
            test_symbol: test_symbol.to_string(),
        });
    }
    Ok(mutations)
}

fn parse_source_bindings(value: &str) -> Result<Vec<SourceBinding>, String> {
    let mut bindings = Vec::new();
    for item in list(value) {
        let Some((source_field, semantic_fields)) = item.split_once('>') else {
            return Err(format!(
                "source binding {item:?} needs Type.field>semantic-field[+semantic-field]"
            ));
        };
        let semantic_fields = semantic_fields
            .split('+')
            .map(str::to_string)
            .collect::<Vec<_>>();
        if source_field.split_once('.').is_none()
            || semantic_fields.is_empty()
            || semantic_fields.iter().any(String::is_empty)
        {
            return Err(format!("source binding {item:?} is not canonical"));
        }
        bindings.push(SourceBinding {
            source_field: source_field.to_string(),
            semantic_fields,
        });
    }
    Ok(bindings)
}

fn take_required(fields: &mut BTreeMap<String, String>, key: &str) -> Result<String, String> {
    fields
        .remove(key)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| format!("missing or empty declaration key {key:?}"))
}

fn parse_declaration(owner: &str, raw: BTreeMap<String, String>) -> Result<IdentityDecl, String> {
    let mut fields = raw;
    let id = take_required(&mut fields, "id")?;
    let version_const = take_required(&mut fields, "version_const")?;
    let version = take_required(&mut fields, "version")?
        .parse::<u32>()
        .map_err(|_| "version must be an unsigned 32-bit integer".to_string())?;
    let domain = take_required(&mut fields, "domain")?;
    let domain_const = match take_required(&mut fields, "domain_const")?.as_str() {
        "none" => None,
        value => Some(parse_schema_constant_reference(value)?.canonical()),
    };
    if domain_const.is_none() {
        return Err(format!(
            "identity {id}: domain_const must name the encoded domain source of truth"
        ));
    }
    let encoder = take_required(&mut fields, "encoder")?;
    let encoder_helpers = list(&take_required(&mut fields, "encoder_helpers")?);
    let schema_functions =
        parse_schema_functions(&take_required(&mut fields, "schema_functions")?)?;
    let schema_constants =
        parse_schema_constants(&take_required(&mut fields, "schema_constants")?)?;
    let schema_dependencies =
        parse_schema_dependencies(&take_required(&mut fields, "schema_dependencies")?, &id)?;
    let digest = take_required(&mut fields, "digest")?;
    let encoding = take_required(&mut fields, "encoding")?;
    if !matches!(
        encoding.as_str(),
        "typed-binary" | "fixed-width-key" | "canonical-transport-exact-bits"
    ) {
        return Err(format!("identity {id}: unsupported encoding {encoding:?}"));
    }
    let sources = list(&take_required(&mut fields, "sources")?);
    if sources.is_empty() {
        return Err(format!(
            "identity {id}: at least one source struct is required"
        ));
    }
    let source_fields = parse_source_fields(&take_required(&mut fields, "source_fields")?)?;
    let source_bindings = parse_source_bindings(&take_required(&mut fields, "source_bindings")?)?;
    let external_semantic_fields = list(&take_required(&mut fields, "external_semantic_fields")?);
    let semantic_fields = list(&take_required(&mut fields, "semantic_fields")?);
    let excluded_fields = parse_reasoned(&take_required(&mut fields, "excluded_fields")?)?;
    let consumers = list(&take_required(&mut fields, "consumers")?);
    let mutations = parse_mutations(&take_required(&mut fields, "mutations")?)?;
    let nonsemantic_mutations =
        parse_mutations(&take_required(&mut fields, "nonsemantic_mutations")?)?;
    let field_guard = take_required(&mut fields, "field_guard")?;
    let transport_guard = take_required(&mut fields, "transport_guard")?;
    let version_guard = take_required(&mut fields, "version_guard")?;
    let coupling_surface = take_required(&mut fields, "coupling_surface")?;
    if let Some((unknown, _)) = fields.pop_first() {
        return Err(format!(
            "identity {id}: unknown declaration key {unknown:?}"
        ));
    }
    if !id.contains(':') || id.chars().any(char::is_whitespace) {
        return Err(format!("identity id {id:?} is not canonical"));
    }
    if !domain_carries_version(&domain, version) {
        return Err(format!(
            "identity {id}: domain {domain:?} must carry an exact v{version} dot/colon segment; domain rotations require a version/coupling bump"
        ));
    }
    if semantic_fields.is_empty() || consumers.is_empty() || mutations.is_empty() {
        return Err(format!(
            "identity {id}: semantic fields, consumers, and mutations must be non-empty"
        ));
    }
    Ok(IdentityDecl {
        id,
        owner: owner.to_string(),
        version_const,
        version,
        domain,
        domain_const,
        encoder,
        encoder_helpers,
        schema_functions,
        schema_constants,
        schema_dependencies,
        digest,
        encoding,
        sources,
        source_fields,
        source_bindings,
        external_semantic_fields,
        semantic_fields,
        excluded_fields,
        consumers,
        mutations,
        nonsemantic_mutations,
        field_guard,
        transport_guard,
        version_guard,
        coupling_surface,
        schema_base_hash: None,
        schema_fingerprint: String::new(),
    })
}

fn declaration_blocks_with_scopes(
    owner: &str,
    text: &str,
    module_scopes: &[RustOwnerScope],
) -> (Vec<IdentityDecl>, Vec<Violation>) {
    let lines = text.lines().collect::<Vec<_>>();
    let mut line_offsets = Vec::with_capacity(lines.len());
    let mut offset = 0_usize;
    for line in text.split_inclusive('\n') {
        line_offsets.push(offset);
        offset += line.len();
    }
    let marker_lines = identity_marker_lines(text);
    let mut declarations = Vec::new();
    let mut violations = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        if !marker_lines.contains(&index) {
            index += 1;
            continue;
        }
        let marker_offset = line_offsets.get(index).copied().unwrap_or(text.len());
        let item_start = text[..marker_offset]
            .rfind("const ")
            .unwrap_or(marker_offset);
        let runtime_active =
            rust_item_is_runtime_active_with_scopes(text, item_start, module_scopes);
        let start = index + 1;
        index = start;
        let mut raw = BTreeMap::new();
        let mut terminated = false;
        while index < lines.len() {
            if lines[index].trim_start().starts_with("];") {
                terminated = true;
                break;
            }
            let Some(literal) = quoted_literal(lines[index]) else {
                violations.push(identity_violation(
                    owner,
                    format!(
                        "{owner}:{}: declaration entries must be one quoted key=value literal per line",
                        index + 1
                    ),
                ));
                index += 1;
                continue;
            };
            let Some((key, value)) = literal.split_once('=') else {
                violations.push(identity_violation(
                    owner,
                    format!(
                        "{owner}:{}: declaration entry {literal:?} lacks =",
                        index + 1
                    ),
                ));
                index += 1;
                continue;
            };
            if raw.insert(key.to_string(), value.to_string()).is_some() {
                violations.push(identity_violation(
                    owner,
                    format!("{owner}:{}: duplicate declaration key {key:?}", index + 1),
                ));
            }
            index += 1;
        }
        if !terminated {
            violations.push(identity_violation(
                owner,
                format!("{owner}:{}: unterminated identity declaration", start),
            ));
            break;
        }
        match parse_declaration(owner, raw) {
            Ok(declaration) if runtime_active => declarations.push(declaration),
            Ok(declaration) => violations.push(identity_violation(
                owner,
                format!(
                    "identity {} declaration is behind cfg/cfg_attr and cannot be a runtime authority",
                    declaration.id
                ),
            )),
            Err(detail) => violations.push(identity_violation(owner, detail)),
        }
        index += 1;
    }
    (declarations, violations)
}

#[cfg(test)]
fn declaration_blocks(owner: &str, text: &str) -> (Vec<IdentityDecl>, Vec<Violation>) {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    declaration_blocks_with_scopes(owner, text, &module_scopes)
}

#[cfg(test)]
fn source_const_u32(text: &str, name: &str) -> Option<u32> {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    source_const_u32_with_scopes(text, name, &module_scopes)
}

fn source_const_u32_with_scopes(
    text: &str,
    name: &str,
    module_scopes: &[RustOwnerScope],
) -> Option<u32> {
    let declarations = runtime_const_declarations_with_scopes(text, name, module_scopes);
    let [declaration] = declarations.as_slice() else {
        return None;
    };
    let (left, right) = declaration.split_once('=')?;
    left.trim_end().ends_with(": u32").then_some(())?;
    right
        .trim()
        .strip_suffix(';')?
        .trim()
        .replace('_', "")
        .parse()
        .ok()
}

#[cfg(test)]
fn source_const_str<'a>(text: &'a str, name: &str) -> Option<&'a str> {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    source_const_str_with_scopes(text, name, &module_scopes)
}

#[cfg(test)]
fn source_const_str_with_scopes<'a>(
    text: &'a str,
    name: &str,
    module_scopes: &[RustOwnerScope],
) -> Option<&'a str> {
    let declarations = runtime_const_declarations_with_scopes(text, name, module_scopes);
    let [declaration] = declarations.as_slice() else {
        return None;
    };
    let (left, right) = declaration.split_once('=')?;
    left.trim_end().ends_with(": &str").then_some(())?;
    let literal = right.trim().strip_suffix(';')?.trim();
    literal.strip_prefix('"')?.strip_suffix('"')
}

fn runtime_const_declarations_with_scopes<'a>(
    text: &'a str,
    symbol: &str,
    module_scopes: &[RustOwnerScope],
) -> Vec<&'a str> {
    const_declarations(text, symbol)
        .into_iter()
        .filter(|declaration| {
            let start = declaration.as_ptr() as usize - text.as_ptr() as usize;
            rust_item_is_runtime_active_with_scopes(text, start, module_scopes)
        })
        .collect()
}

fn char_literal_end(bytes: &[u8], start: usize) -> Option<usize> {
    let first = *bytes.get(start + 1)?;
    let payload_end = if first == b'\\' {
        match *bytes.get(start + 2)? {
            b'x' => start.checked_add(5)?,
            b'u' if bytes.get(start + 3) == Some(&b'{') => bytes[start + 4..]
                .iter()
                .position(|byte| *byte == b'}')?
                .checked_add(start + 5)?,
            _ => start.checked_add(3)?,
        }
    } else {
        let width = match first {
            0x00..=0x7f => 1,
            0xc0..=0xdf => 2,
            0xe0..=0xef => 3,
            0xf0..=0xf7 => 4,
            _ => return None,
        };
        start.checked_add(1 + width)?
    };
    (bytes.get(payload_end) == Some(&b'\'')).then_some(payload_end)
}

fn rust_function_starts(text: &str) -> Vec<(String, usize)> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = text.as_bytes();
    let mut functions = Vec::new();
    let mut state = State::Normal;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if bytes[index..].starts_with(b"fn") {
                    let before_is_ident = index > 0
                        && matches!(bytes[index - 1], b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z');
                    let after_is_ident = bytes.get(index + 2).is_some_and(
                        |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                    );
                    if !before_is_ident && !after_is_ident {
                        let mut cursor = index + 2;
                        while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                            cursor += 1;
                        }
                        let symbol_start = cursor;
                        while bytes.get(cursor).is_some_and(
                            |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                        ) {
                            cursor += 1;
                        }
                        if let Some(symbol) = text.get(symbol_start..cursor)
                            && canonical_symbol(symbol)
                        {
                            functions.push((symbol.to_string(), index));
                            index = cursor;
                            continue;
                        }
                    }
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth - 1;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    functions
}

/// Return every matching declaration body while ignoring braces in Rust
/// literals and comments. `direct_only` limits matches to the current lexical
/// scope, which lets qualified function lookup distinguish modules from impls.
#[allow(clippy::too_many_lines)]
fn braced_bodies_limited<'a>(
    text: &'a str,
    declaration: &str,
    direct_only: bool,
    limit: usize,
) -> Vec<&'a str> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = text.as_bytes();
    let declaration_bytes = declaration.as_bytes();
    let mut state = State::Normal;
    let mut index = 0usize;
    let mut awaiting_brace = false;
    let mut open = None;
    let mut brace_depth = 0usize;
    let mut scope_depth = 0usize;
    let mut signature_paren_depth = 0usize;
    let mut signature_bracket_depth = 0usize;
    let mut signature_brace_depth = 0usize;
    let mut signature_angle_depth = 0usize;
    let mut bodies = Vec::new();
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if open.is_none()
                    && !awaiting_brace
                    && (!direct_only || scope_depth == 0)
                    && bytes[index..].starts_with(declaration_bytes)
                {
                    let before_is_ident = index > 0
                        && matches!(bytes[index - 1], b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z');
                    let after = index + declaration_bytes.len();
                    let after_is_ident = bytes.get(after).is_some_and(
                        |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                    );
                    if !before_is_ident && !after_is_ident {
                        awaiting_brace = true;
                        signature_paren_depth = 0;
                        signature_bracket_depth = 0;
                        signature_brace_depth = 0;
                        signature_angle_depth = 0;
                        index = after;
                        continue;
                    }
                }
                if awaiting_brace {
                    match byte {
                        b'(' => signature_paren_depth += 1,
                        b')' => {
                            signature_paren_depth = signature_paren_depth.saturating_sub(1);
                        }
                        b'[' => signature_bracket_depth += 1,
                        b']' => {
                            signature_bracket_depth = signature_bracket_depth.saturating_sub(1);
                        }
                        b'<' => signature_angle_depth += 1,
                        b'>' if signature_angle_depth > 0 => signature_angle_depth -= 1,
                        b'{' if signature_paren_depth > 0
                            || signature_bracket_depth > 0
                            || signature_angle_depth > 0
                            || signature_brace_depth > 0 =>
                        {
                            signature_brace_depth += 1;
                        }
                        b'}' if signature_brace_depth > 0 => signature_brace_depth -= 1,
                        b';' if signature_paren_depth == 0
                            && signature_bracket_depth == 0
                            && signature_angle_depth == 0
                            && signature_brace_depth == 0 =>
                        {
                            awaiting_brace = false;
                        }
                        b'{' if signature_paren_depth == 0
                            && signature_bracket_depth == 0
                            && signature_angle_depth == 0
                            && signature_brace_depth == 0 =>
                        {
                            open = Some(index);
                            brace_depth = 1;
                            awaiting_brace = false;
                        }
                        _ => {}
                    }
                } else if byte == b'{' {
                    if open.is_none() {
                        if direct_only {
                            scope_depth += 1;
                        }
                    } else {
                        brace_depth += 1;
                    }
                } else if byte == b'}' {
                    if let Some(body_start) = open {
                        let Some(next_depth) = brace_depth.checked_sub(1) else {
                            return bodies;
                        };
                        brace_depth = next_depth;
                        if brace_depth == 0 {
                            bodies.push(&text[body_start + 1..index]);
                            if bodies.len() == limit {
                                return bodies;
                            }
                            open = None;
                        }
                    } else if direct_only {
                        let Some(next_depth) = scope_depth.checked_sub(1) else {
                            return bodies;
                        };
                        scope_depth = next_depth;
                    }
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let Some(depth) = depth.checked_sub(1) else {
                        return bodies;
                    };
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    bodies
}

fn braced_bodies<'a>(text: &'a str, declaration: &str, direct_only: bool) -> Vec<&'a str> {
    braced_bodies_limited(text, declaration, direct_only, usize::MAX)
}

fn braced_body<'a>(text: &'a str, declaration: &str) -> Option<&'a str> {
    braced_bodies(text, declaration, false).into_iter().next()
}

#[cfg(test)]
fn direct_braced_bodies<'a>(text: &'a str, declaration: &str) -> Vec<&'a str> {
    braced_bodies(text, declaration, true)
}

fn first_direct_braced_body<'a>(text: &'a str, declaration: &str) -> Option<&'a str> {
    braced_bodies_limited(text, declaration, true, 1)
        .into_iter()
        .next()
}

#[allow(clippy::too_many_lines)] // The lexer prevents nested Rust syntax from splitting fields.
fn split_rust_top_level(fragment: &str, delimiter: u8) -> Vec<&str> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = fragment.as_bytes();
    let mut state = State::Normal;
    let mut parens = 0usize;
    let mut brackets = 0usize;
    let mut braces = 0usize;
    let mut angles = 0usize;
    let mut start = 0usize;
    let mut index = 0usize;
    let mut items = Vec::new();
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if byte == delimiter && parens == 0 && brackets == 0 && braces == 0 && angles == 0 {
                    items.push(&fragment[start..index]);
                    start = index + 1;
                    index += 1;
                    continue;
                }
                match byte {
                    b'(' => parens += 1,
                    b')' => parens = parens.saturating_sub(1),
                    b'[' => brackets += 1,
                    b']' => brackets = brackets.saturating_sub(1),
                    b'{' => braces += 1,
                    b'}' => braces = braces.saturating_sub(1),
                    b'<' => angles += 1,
                    b'>' if angles > 0 => angles -= 1,
                    _ => {}
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    state = if depth == 1 {
                        State::Normal
                    } else {
                        State::BlockComment { depth: depth - 1 }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    items.push(&fragment[start..]);
    items
}

fn last_identifier(fragment: &str) -> Option<&str> {
    fragment
        .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
        .filter(|token| !token.is_empty())
        .next_back()
}

fn named_fields(body: &str) -> BTreeSet<String> {
    let mut fields = BTreeSet::new();
    for item in split_rust_top_level(body, b',') {
        let parts = split_rust_top_level(item, b':');
        let Some(prefix) = parts.first() else {
            continue;
        };
        let Some(field) = last_identifier(prefix) else {
            continue;
        };
        if field
            .chars()
            .all(|ch| ch == '_' || ch.is_ascii_alphanumeric())
        {
            fields.insert(field.to_string());
        }
    }
    fields
}

#[allow(clippy::too_many_lines)] // Attributes/comments must not masquerade as variants.
fn first_top_level_identifier(fragment: &str) -> Option<&str> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = fragment.as_bytes();
    let mut state = State::Normal;
    let mut brackets = 0usize;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                match byte {
                    b'[' => brackets += 1,
                    b']' => brackets = brackets.saturating_sub(1),
                    b'_' | b'a'..=b'z' | b'A'..=b'Z' if brackets == 0 => {
                        let start = index;
                        index += 1;
                        while bytes.get(index).is_some_and(
                            |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                        ) {
                            index += 1;
                        }
                        return fragment.get(start..index);
                    }
                    _ => {}
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    state = if depth == 1 {
                        State::Normal
                    } else {
                        State::BlockComment { depth: depth - 1 }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    None
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SourceShape {
    Struct(BTreeSet<String>),
    Enum {
        fields: BTreeSet<String>,
        variants: BTreeSet<String>,
        tuple_variants: BTreeSet<String>,
    },
}

impl SourceShape {
    fn fields(&self) -> &BTreeSet<String> {
        match self {
            Self::Struct(fields) | Self::Enum { fields, .. } => fields,
        }
    }
}

#[cfg(test)]
fn source_shape(text: &str, name: &str) -> Option<SourceShape> {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    source_shape_with_scopes(text, name, &module_scopes)
}

fn source_shape_with_scopes(
    text: &str,
    name: &str,
    module_scopes: &[RustOwnerScope],
) -> Option<SourceShape> {
    let structs = runtime_braced_bodies_with_scopes(text, &format!("struct {name}"), module_scopes);
    let enums = runtime_braced_bodies_with_scopes(text, &format!("enum {name}"), module_scopes);
    if let ([body], []) = (structs.as_slice(), enums.as_slice()) {
        return Some(SourceShape::Struct(named_fields(body)));
    }
    let ([], [body]) = (structs.as_slice(), enums.as_slice()) else {
        return None;
    };
    let mut fields = BTreeSet::from(["variant".to_string()]);
    let mut variants = BTreeSet::new();
    let mut tuple_variants = BTreeSet::new();
    for item in split_rust_top_level(body, b',') {
        let Some(variant) = first_top_level_identifier(item) else {
            continue;
        };
        variants.insert(variant.to_string());
        if let Some(payload) = braced_body(item, variant) {
            fields.extend(named_fields(payload));
        } else if item
            .trim_start()
            .strip_prefix(variant)
            .is_some_and(|rest| rest.trim_start().starts_with('('))
        {
            tuple_variants.insert(variant.to_string());
        }
    }
    (!variants.is_empty()).then_some(SourceShape::Enum {
        fields,
        variants,
        tuple_variants,
    })
}

fn symbol_body_has_no_rest_pattern(text: &str, symbol: &str) -> bool {
    let Some(body) = braced_body(text, &format!("fn {symbol}")) else {
        return false;
    };
    !body.contains("..")
}

#[cfg(test)]
fn function_segments(reference: &str) -> Option<Vec<&str>> {
    let segments = reference.split("::").collect::<Vec<_>>();
    canonical_function_reference(reference).then_some(segments)
}

#[cfg(test)]
fn module_scopes<'a>(text: &'a str, modules: &[&str]) -> Vec<&'a str> {
    let mut scopes = vec![text];
    for module in modules {
        let declaration = format!("mod {module}");
        scopes = scopes
            .into_iter()
            .flat_map(|scope| braced_bodies(scope, &declaration, false))
            .collect();
    }
    scopes
}

fn leading_implementation_type(rest: &str) -> Option<String> {
    let rest = rest.trim_start_matches(|character: char| {
        character.is_ascii_whitespace() || matches!(character, '&' | '\'' | '!')
    });
    let self_type = rest
        .split(|character: char| {
            character.is_ascii_whitespace() || matches!(character, '<' | '{' | '(')
        })
        .next()?;
    let owner = self_type.rsplit("::").next()?;
    canonical_symbol(owner).then(|| owner.to_string())
}

fn implementation_owner(header: &str) -> Option<String> {
    let tokens = header
        .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
        .filter(|token| !token.is_empty())
        .collect::<Vec<_>>();
    let impl_index = tokens.iter().position(|token| *token == "impl")?;
    if let Some(for_index) = tokens.iter().rposition(|token| *token == "for") {
        if for_index <= impl_index {
            return None;
        }
        let (at, _) = header
            .match_indices("for")
            .filter(|(at, _)| {
                let before = at
                    .checked_sub(1)
                    .and_then(|index| header.as_bytes().get(index));
                let after = header.as_bytes().get(at + "for".len());
                before.is_none_or(
                    |byte| !matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                ) && after.is_none_or(
                    |byte| !matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                )
            })
            .last()?;
        return leading_implementation_type(&header[at + "for".len()..]);
    }

    let mut rest = header
        .get(header.find("impl")? + "impl".len()..)?
        .trim_start();
    if rest.starts_with('<') {
        let mut depth = 0_usize;
        let mut end = None;
        for (index, character) in rest.char_indices() {
            match character {
                '<' => depth += 1,
                '>' => {
                    depth = depth.checked_sub(1)?;
                    if depth == 0 {
                        end = Some(index + character.len_utf8());
                        break;
                    }
                }
                _ => {}
            }
        }
        rest = rest.get(end?..)?.trim_start();
    }
    leading_implementation_type(rest)
}

#[cfg(test)]
fn implementation_bodies<'a>(text: &'a str, owner: &str) -> Vec<&'a str> {
    let mut bodies = Vec::new();
    for fragment in braced_declaration_fragments(text, "impl") {
        let Some(body) = braced_bodies(fragment, "impl", false).into_iter().next() else {
            continue;
        };
        let body_start = body.as_ptr() as usize - fragment.as_ptr() as usize;
        let Some(open) = body_start.checked_sub(1) else {
            continue;
        };
        let Some(header) = fragment.get(..open) else {
            continue;
        };
        if implementation_owner(header).as_deref() == Some(owner) {
            bodies.push(body);
        }
    }
    bodies.sort_by_key(|body| body.as_ptr() as usize);
    bodies.dedup_by_key(|body| body.as_ptr() as usize);
    bodies
}

#[cfg(test)]
fn function_bodies<'a>(text: &'a str, reference: &str) -> Vec<&'a str> {
    let Some(segments) = function_segments(reference) else {
        return Vec::new();
    };
    let symbol = segments[segments.len() - 1];
    let declaration = format!("fn {symbol}");
    if segments.len() == 1 {
        return direct_braced_bodies(text, &declaration);
    }

    let mut bodies = module_scopes(text, &segments[..segments.len() - 1])
        .into_iter()
        .flat_map(|scope| direct_braced_bodies(scope, &declaration))
        .collect::<Vec<_>>();

    let owner = segments[segments.len() - 2];
    for scope in module_scopes(text, &segments[..segments.len() - 2]) {
        for implementation in implementation_bodies(scope, owner) {
            bodies.extend(direct_braced_bodies(implementation, &declaration));
        }
    }
    bodies.sort_by_key(|body| body.as_ptr() as usize);
    bodies.dedup_by_key(|body| body.as_ptr() as usize);
    bodies
}

#[cfg(test)]
fn function_body<'a>(text: &'a str, reference: &str) -> Option<&'a str> {
    let bodies = function_bodies(text, reference);
    let [body] = bodies.as_slice() else {
        return None;
    };
    Some(*body)
}

fn braced_declaration_fragments<'a>(text: &'a str, declaration: &str) -> Vec<&'a str> {
    let mut fragments = Vec::new();
    for body in braced_bodies(text, declaration, false) {
        let body_start = body.as_ptr() as usize - text.as_ptr() as usize;
        let close = body_start + body.len();
        let start = text[..body_start]
            .match_indices(declaration)
            .filter_map(|(start, _)| {
                let candidate = text.get(start..=close)?;
                braced_bodies(candidate, declaration, false)
                    .first()
                    .is_some_and(|candidate_body| candidate_body.as_ptr() == body.as_ptr())
                    .then_some(start)
            })
            .last();
        if let Some(start) = start
            && let Some(fragment) = text.get(start..=close)
        {
            fragments.push(fragment);
        }
    }
    fragments
}

fn runtime_braced_bodies_with_scopes<'a>(
    text: &'a str,
    declaration: &str,
    module_scopes: &[RustOwnerScope],
) -> Vec<&'a str> {
    braced_declaration_fragments(text, declaration)
        .into_iter()
        .filter_map(|fragment| {
            let start = fragment.as_ptr() as usize - text.as_ptr() as usize;
            rust_item_is_runtime_active_with_scopes(text, start, module_scopes)
                .then(|| {
                    braced_bodies(fragment, declaration, false)
                        .into_iter()
                        .next()
                })
                .flatten()
        })
        .collect()
}

#[derive(Clone)]
struct RustOwnerScope {
    start: usize,
    end: usize,
    name: String,
    cfg_attributes: Vec<String>,
}

fn rust_item_locator(
    path: &str,
    module_scopes: &[RustOwnerScope],
    item_start: usize,
    kind: &str,
    symbol: &str,
) -> String {
    let mut scopes = module_scopes
        .iter()
        .filter(|scope| scope.start <= item_start && item_start < scope.end)
        .collect::<Vec<_>>();
    scopes.sort_by(|left, right| {
        (left.start, std::cmp::Reverse(left.end)).cmp(&(right.start, std::cmp::Reverse(right.end)))
    });
    let module = if scopes.is_empty() {
        "<root>".to_string()
    } else {
        scopes
            .into_iter()
            .map(|scope| scope.name.as_str())
            .collect::<Vec<_>>()
            .join("::")
    };
    format!("{path}#{module}::{kind}:{symbol}")
}

fn attached_rust_attributes(text: &str, item_start: usize) -> Vec<String> {
    let prefix = &text[..item_start.min(text.len())];
    let lines = prefix.lines().collect::<Vec<_>>();
    let mut cursor = lines.len();
    let mut attributes = Vec::new();
    while cursor > 0 {
        cursor -= 1;
        let line = lines[cursor].trim();
        if line.is_empty()
            || line.starts_with("//")
            || line.starts_with("/*")
            || line.starts_with('*')
            || line.ends_with("*/")
        {
            continue;
        }
        let attribute = if line.starts_with("#[") {
            Some(line.to_string())
        } else if line.ends_with(']') {
            let end = cursor;
            while cursor > 0 && !lines[cursor].trim().starts_with("#[") {
                cursor -= 1;
            }
            lines[cursor]
                .trim()
                .starts_with("#[")
                .then(|| lines[cursor..=end].join(""))
        } else {
            None
        };
        let Some(attribute) = attribute else {
            break;
        };
        attributes.push(
            attribute
                .chars()
                .filter(|character| !character.is_ascii_whitespace())
                .collect(),
        );
    }
    attributes
}

fn cfg_attributes(text: &str, item_start: usize) -> Vec<String> {
    attached_rust_attributes(text, item_start)
        .into_iter()
        .filter(|attribute| attribute.starts_with("#[cfg(") || attribute.starts_with("#[cfg_attr("))
        .collect()
}

fn rust_owner_scopes(
    text: &str,
    declaration: &str,
    owner: impl Fn(&str) -> Option<String>,
) -> Vec<RustOwnerScope> {
    let mut scopes = Vec::new();
    for fragment in braced_declaration_fragments(text, declaration) {
        let Some(body) = braced_bodies(fragment, declaration, false)
            .into_iter()
            .next()
        else {
            continue;
        };
        let body_start_in_fragment = body.as_ptr() as usize - fragment.as_ptr() as usize;
        let Some(open) = body_start_in_fragment.checked_sub(1) else {
            continue;
        };
        let Some(name) = fragment.get(..open).and_then(&owner) else {
            continue;
        };
        let fragment_start = fragment.as_ptr() as usize - text.as_ptr() as usize;
        let start = body.as_ptr() as usize - text.as_ptr() as usize;
        scopes.push(RustOwnerScope {
            start,
            end: start + body.len(),
            name,
            cfg_attributes: cfg_attributes(text, fragment_start),
        });
    }
    scopes.sort_by(|left, right| (left.start, left.end).cmp(&(right.start, right.end)));
    scopes.dedup_by(|left, right| {
        left.start == right.start && left.end == right.end && left.name == right.name
    });
    scopes
}

fn rust_module_owner(header: &str) -> Option<String> {
    let tokens = header
        .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
        .filter(|token| !token.is_empty())
        .collect::<Vec<_>>();
    let module = tokens.iter().position(|token| *token == "mod")?;
    let name = *tokens.get(module + 1)?;
    canonical_symbol(name).then(|| name.to_string())
}

fn rust_item_is_runtime_active_with_scopes(
    text: &str,
    item_start: usize,
    module_scopes: &[RustOwnerScope],
) -> bool {
    cfg_attributes(text, item_start).is_empty()
        && module_scopes
            .iter()
            .filter(|scope| scope.start <= item_start && item_start < scope.end)
            .all(|scope| scope.cfg_attributes.is_empty())
}

fn rust_item_is_test_active_with_scopes(
    item_start: usize,
    module_scopes: &[RustOwnerScope],
) -> bool {
    module_scopes
        .iter()
        .filter(|scope| scope.start <= item_start && item_start < scope.end)
        .flat_map(|scope| scope.cfg_attributes.iter())
        .all(|attribute| attribute == "#[cfg(test)]")
}

fn rust_function_fragment_index_with_scopes<'a>(
    text: &'a str,
    modules: &[RustOwnerScope],
    implementations: &[RustOwnerScope],
) -> BTreeMap<String, Vec<&'a str>> {
    struct FunctionSpan<'a> {
        symbol: String,
        start: usize,
        body_start: usize,
        body_end: usize,
        fragment: &'a str,
    }

    let starts = rust_function_starts(text);
    let mut spans = Vec::new();
    for (index, (symbol, start)) in starts.iter().enumerate() {
        let declaration = format!("fn {symbol}");
        let Some(body) = first_direct_braced_body(&text[*start..], &declaration) else {
            continue;
        };
        let body_start = body.as_ptr() as usize - text.as_ptr() as usize;
        if starts
            .get(index + 1)
            .is_some_and(|(_, next_start)| *next_start < body_start)
        {
            continue;
        }
        let body_end = body_start + body.len();
        let Some(fragment) = text.get(*start..=body_end) else {
            continue;
        };
        spans.push(FunctionSpan {
            symbol: symbol.clone(),
            start: *start,
            body_start,
            body_end,
            fragment,
        });
    }

    let mut index = BTreeMap::<String, Vec<&str>>::new();
    let mut active_function_ends = Vec::<usize>::new();
    for span in &spans {
        active_function_ends.retain(|end| span.start < *end);
        let nested_in_function = !active_function_ends.is_empty();
        let mut module_scopes = modules
            .iter()
            .filter(|scope| scope.start <= span.start && span.start < scope.end)
            .collect::<Vec<_>>();
        module_scopes.sort_by(|left, right| {
            (left.start, std::cmp::Reverse(left.end))
                .cmp(&(right.start, std::cmp::Reverse(right.end)))
        });
        let modules = module_scopes
            .into_iter()
            .map(|scope| scope.name.as_str())
            .collect::<Vec<_>>();
        let implementation = (!nested_in_function)
            .then(|| {
                implementations
                    .iter()
                    .filter(|scope| scope.start <= span.start && span.start < scope.end)
                    .min_by_key(|scope| scope.end - scope.start)
            })
            .flatten();

        let mut full = modules.clone();
        if let Some(implementation) = implementation {
            full.push(implementation.name.as_str());
        }
        full.push(span.symbol.as_str());
        let full = full.join("::");
        index.entry(full.clone()).or_default().push(span.fragment);
        if let Some(implementation) = implementation {
            let alias = format!("{}::{}", implementation.name, span.symbol);
            if alias != full {
                index.entry(alias).or_default().push(span.fragment);
            }
        }
        active_function_ends.push(span.body_end);
        active_function_ends.sort_unstable();
        let _ = span.body_start;
    }
    for fragments in index.values_mut() {
        fragments.sort_by_key(|fragment| fragment.as_ptr() as usize);
        fragments.dedup_by_key(|fragment| fragment.as_ptr() as usize);
    }
    index
}

struct RustSourceIndex<'a> {
    functions: BTreeMap<String, Vec<&'a str>>,
    module_scopes: Vec<RustOwnerScope>,
}

impl<'a> RustSourceIndex<'a> {
    fn new(text: &'a str) -> Self {
        let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
        let implementations = rust_owner_scopes(text, "impl", implementation_owner);
        Self {
            functions: rust_function_fragment_index_with_scopes(
                text,
                &module_scopes,
                &implementations,
            ),
            module_scopes,
        }
    }
}

fn has_function_with_index(text: &str, index: &RustSourceIndex<'_>, reference: &str) -> bool {
    let Some(fragments) = index.functions.get(reference) else {
        return false;
    };
    let [fragment] = fragments.as_slice() else {
        return false;
    };
    let start = fragment.as_ptr() as usize - text.as_ptr() as usize;
    rust_item_is_runtime_active_with_scopes(text, start, &index.module_scopes)
}

fn has_function(text: &str, reference: &str) -> bool {
    let index = RustSourceIndex::new(text);
    has_function_with_index(text, &index, reference)
}

fn has_test_function_with_scopes(
    text: &str,
    symbol: &str,
    module_scopes: &[RustOwnerScope],
) -> bool {
    let declaration = format!("fn {symbol}");
    let starts = rust_function_starts(text)
        .into_iter()
        .filter_map(|(candidate, start)| (candidate == symbol).then_some(start))
        .collect::<Vec<_>>();
    let [at] = starts.as_slice() else {
        return false;
    };
    if !rust_item_is_test_active_with_scopes(*at, module_scopes) {
        return false;
    }
    let bodies = braced_bodies(text, &declaration, false);
    let [body] = bodies.as_slice() else {
        return false;
    };
    if normalized_rust_fragment(body).is_empty() {
        return false;
    }
    let prefix = &text[..*at];
    let line_start = prefix.rfind('\n').map_or(0, |index| index + 1);
    let mut saw_test = false;
    let lines = prefix[..line_start].lines().collect::<Vec<_>>();
    let mut cursor = lines.len();
    while cursor > 0 {
        cursor -= 1;
        let line = lines[cursor].trim();
        if line.is_empty()
            || line.starts_with("//")
            || line.starts_with("/*")
            || line.starts_with('*')
            || line.ends_with("*/")
        {
            continue;
        }
        let attribute = if line.starts_with("#[") {
            Some(line.to_string())
        } else if line.ends_with(']') {
            let end = cursor;
            while cursor > 0 && !lines[cursor].trim().starts_with("#[") {
                cursor -= 1;
            }
            lines[cursor]
                .trim()
                .starts_with("#[")
                .then(|| lines[cursor..=end].join(""))
        } else {
            None
        };
        if let Some(attribute) = attribute {
            if attribute.starts_with("#[ignore")
                || attribute.starts_with("#[should_panic")
                || attribute.starts_with("#[cfg(")
                || attribute.starts_with("#[cfg_attr(")
            {
                return false;
            }
            saw_test |= attribute == "#[test]";
            continue;
        }
        break;
    }
    saw_test
}

fn has_test_function(text: &str, symbol: &str) -> bool {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    has_test_function_with_scopes(text, symbol, &module_scopes)
}

fn guard_destructures_sources_with_scopes(
    text: &str,
    symbol: &str,
    sources: &[String],
    module_scopes: &[RustOwnerScope],
) -> bool {
    let Some(body) = braced_body(text, &format!("fn {symbol}")) else {
        return false;
    };
    sources.iter().all(
        |source| match source_shape_with_scopes(text, source, module_scopes) {
            Some(SourceShape::Struct(_)) => body.contains(&format!("let {source} {{")),
            Some(SourceShape::Enum { variants, .. }) => {
                body.contains("match ")
                    && variants
                        .iter()
                        .all(|variant| body.contains(&format!("{source}::{variant}")))
            }
            None => false,
        },
    )
}

#[cfg(test)]
fn guard_destructures_sources(text: &str, symbol: &str, sources: &[String]) -> bool {
    let module_scopes = rust_owner_scopes(text, "mod", rust_module_owner);
    guard_destructures_sources_with_scopes(text, symbol, sources, &module_scopes)
}

#[allow(clippy::too_many_lines)] // One small lexer keeps schema hashes formatting-insensitive.
fn normalized_rust_fragment(fragment: &str) -> Vec<u8> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = fragment.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut state = State::Normal;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte.is_ascii_whitespace() {
                    index += 1;
                    continue;
                }
                if bytes.get(index + 1) == Some(&b'/') && byte == b'/' {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if bytes.get(index + 1) == Some(&b'*') && byte == b'/' {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut cursor = index + 1;
                    while bytes.get(cursor) == Some(&b'#') {
                        cursor += 1;
                    }
                    if bytes.get(cursor) == Some(&b'"') {
                        out.extend_from_slice(&bytes[index..=cursor]);
                        state = State::Raw {
                            hashes: cursor - index - 1,
                        };
                        index = cursor + 1;
                        continue;
                    }
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    out.extend_from_slice(&bytes[index..=end]);
                    index = end + 1;
                    continue;
                }
                out.push(byte);
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                out.push(byte);
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                out.push(byte);
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    out.extend_from_slice(&bytes[index..index + hashes]);
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth - 1;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    out
}

fn identifier_tokens(fragment: &str) -> BTreeSet<String> {
    fragment
        .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
        .filter(|token| canonical_symbol(token))
        .map(str::to_string)
        .collect()
}

fn append_schema_frame(out: &mut Vec<u8>, label: &str, bytes: &[u8]) {
    let label_len = u64::try_from(label.len()).expect("schema label length fits u64");
    let bytes_len = u64::try_from(bytes.len()).expect("schema fragment length fits u64");
    out.extend_from_slice(&label_len.to_le_bytes());
    out.extend_from_slice(label.as_bytes());
    out.extend_from_slice(&bytes_len.to_le_bytes());
    out.extend_from_slice(bytes);
}

fn normalized_rust_function_closure_with_symbols_and_index(
    text: &str,
    index: &RustSourceIndex<'_>,
    roots: impl IntoIterator<Item = String>,
) -> Result<(Vec<u8>, BTreeSet<String>), String> {
    let mut pending = roots.into_iter().collect::<BTreeSet<_>>();
    let mut seen = BTreeSet::new();
    let mut constants = BTreeSet::new();
    let mut out = Vec::new();
    while let Some(reference) = pending.pop_first() {
        if !seen.insert(reference.clone()) {
            continue;
        }
        let fragments = index
            .functions
            .get(&reference)
            .map_or(&[][..], Vec::as_slice);
        let [fragment] = fragments else {
            return Err(format!(
                "function {reference:?} must resolve to exactly one item; found {}",
                fragments.len()
            ));
        };
        let fragment_start = fragment.as_ptr() as usize - text.as_ptr() as usize;
        if !rust_item_is_runtime_active_with_scopes(text, fragment_start, &index.module_scopes) {
            return Err(format!(
                "function {reference:?} is behind cfg/cfg_attr and cannot be a runtime identity schema authority"
            ));
        }
        append_schema_frame(
            &mut out,
            &format!("fn:{reference}"),
            &normalized_rust_fragment(fragment),
        );
        let scope = reference.rsplit_once("::").map(|(scope, _)| scope);
        for token in identifier_tokens(fragment) {
            if token
                .bytes()
                .all(|byte| byte == b'_' || byte.is_ascii_digit() || byte.is_ascii_uppercase())
                && runtime_const_declarations_with_scopes(text, &token, &index.module_scopes).len()
                    == 1
            {
                constants.insert(token.clone());
            }
            let mut candidates = BTreeSet::from([token.clone()]);
            if let Some(scope) = scope {
                candidates.insert(format!("{scope}::{token}"));
            }
            for candidate in candidates {
                if !seen.contains(&candidate)
                    && index
                        .functions
                        .get(&candidate)
                        .is_some_and(|fragments| fragments.len() == 1)
                {
                    pending.insert(candidate);
                }
            }
        }
    }
    for constant in constants {
        let declarations =
            runtime_const_declarations_with_scopes(text, &constant, &index.module_scopes);
        let [declaration] = declarations.as_slice() else {
            continue;
        };
        append_schema_frame(
            &mut out,
            &format!("const:{constant}"),
            &normalized_rust_fragment(declaration),
        );
    }
    Ok((out, seen))
}

fn normalized_rust_function_closure_with_symbols(
    text: &str,
    roots: impl IntoIterator<Item = String>,
) -> Result<(Vec<u8>, BTreeSet<String>), String> {
    let index = RustSourceIndex::new(text);
    normalized_rust_function_closure_with_symbols_and_index(text, &index, roots)
}

fn normalized_rust_function_closure_with_index(
    text: &str,
    index: &RustSourceIndex<'_>,
    roots: impl IntoIterator<Item = String>,
) -> Result<Vec<u8>, String> {
    normalized_rust_function_closure_with_symbols_and_index(text, index, roots)
        .map(|(bytes, _)| bytes)
}

fn normalized_rust_function_closure(
    text: &str,
    roots: impl IntoIterator<Item = String>,
) -> Result<Vec<u8>, String> {
    normalized_rust_function_closure_with_symbols(text, roots).map(|(bytes, _)| bytes)
}

#[allow(clippy::too_many_lines)] // Exact const spans need the same literal/comment discipline.
fn const_item_end(text: &str, start: usize) -> Option<usize> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = text.as_bytes();
    let mut state = State::Normal;
    let mut parens = 0usize;
    let mut brackets = 0usize;
    let mut braces = 0usize;
    let mut index = start;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                match byte {
                    b'(' => parens += 1,
                    b')' => parens = parens.checked_sub(1)?,
                    b'[' => brackets += 1,
                    b']' => brackets = brackets.checked_sub(1)?,
                    b'{' => braces += 1,
                    b'}' => braces = braces.checked_sub(1)?,
                    b';' if parens == 0 && brackets == 0 && braces == 0 => {
                        return Some(index + 1);
                    }
                    _ => {}
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth.checked_sub(1)?;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    None
}

#[allow(clippy::too_many_lines)] // Rejecting comment/string decoys is part of exact resolution.
fn const_declarations<'a>(text: &'a str, symbol: &str) -> Vec<&'a str> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = text.as_bytes();
    let marker = format!("const {symbol}");
    let marker = marker.as_bytes();
    let mut declarations = Vec::new();
    let mut state = State::Normal;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if bytes[index..].starts_with(marker) {
                    let before_is_ident = index > 0
                        && matches!(bytes[index - 1], b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z');
                    let after = index + marker.len();
                    let after_is_ident = bytes.get(after).is_some_and(
                        |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                    );
                    let has_type = text[after..].trim_start().starts_with(':');
                    if !before_is_ident
                        && !after_is_ident
                        && has_type
                        && let Some(end) = const_item_end(text, index)
                    {
                        declarations.push(&text[index..end]);
                        index = end;
                        continue;
                    }
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth - 1;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    declarations
}

fn identity_test_bytes_with_index(
    text: &str,
    index: &RustSourceIndex<'_>,
    path: &str,
    symbol: &str,
) -> Result<Vec<u8>, String> {
    if !has_test_function_with_scopes(text, symbol, &index.module_scopes) {
        return Err(format!(
            "identity evidence {path}#{symbol} must resolve to one active nonempty #[test] function"
        ));
    }
    let mut fragments = index
        .functions
        .get(symbol)
        .into_iter()
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    if fragments.is_empty() {
        fragments.extend(
            index
                .functions
                .iter()
                .filter(|(reference, _)| {
                    reference.rsplit_once("::").map(|(_, tail)| tail) == Some(symbol)
                })
                .flat_map(|(_, fragments)| fragments.iter().copied()),
        );
    }
    fragments.sort_by_key(|fragment| fragment.as_ptr() as usize);
    fragments.dedup_by_key(|fragment| fragment.as_ptr() as usize);
    let [fragment] = fragments.as_slice() else {
        return Err(format!(
            "identity evidence {path}#{symbol} must resolve to exactly one function item; found {}",
            fragments.len()
        ));
    };
    Ok(normalized_rust_fragment(fragment))
}

#[derive(Default)]
struct IdentityReferenceRequests {
    tests: BTreeSet<String>,
    functions: BTreeSet<String>,
    constants: BTreeSet<String>,
}

#[derive(Default)]
struct IdentityReferenceCache {
    tests: BTreeMap<(String, String), Result<Vec<u8>, String>>,
    functions: BTreeMap<(String, String), Result<Vec<u8>, String>>,
    constants: BTreeMap<(String, String), Result<IdentityConstantReference, String>>,
}

struct IdentityConstantReference {
    locator: String,
    declaration: Vec<u8>,
}

impl IdentityReferenceCache {
    fn build<'a>(
        root: &Path,
        declarations: impl IntoIterator<Item = &'a IdentityDecl>,
        inline_sources: &BTreeMap<String, &str>,
    ) -> Self {
        let mut requests = BTreeMap::<String, IdentityReferenceRequests>::new();
        for declaration in declarations {
            for mutation in declaration
                .mutations
                .iter()
                .chain(declaration.nonsemantic_mutations.iter())
            {
                requests
                    .entry(mutation.test_path.clone())
                    .or_default()
                    .tests
                    .insert(mutation.test_symbol.clone());
            }
            if let Some((path, symbol)) = declaration.version_guard.split_once('#') {
                requests
                    .entry(path.to_string())
                    .or_default()
                    .tests
                    .insert(symbol.to_string());
            }
            for function in declaration
                .schema_functions
                .iter()
                .filter(|function| function.path.is_some())
            {
                requests
                    .entry(function.path.clone().expect("filtered path is present"))
                    .or_default()
                    .functions
                    .insert(function.symbol.clone());
            }
            for constant in &declaration.schema_constants {
                requests
                    .entry(
                        constant
                            .path
                            .clone()
                            .unwrap_or_else(|| declaration.owner.clone()),
                    )
                    .or_default()
                    .constants
                    .insert(constant.symbol.clone());
            }
            if let Some(domain_const) = declaration
                .domain_const
                .as_deref()
                .and_then(|value| parse_schema_constant_reference(value).ok())
            {
                requests
                    .entry(
                        domain_const
                            .path
                            .clone()
                            .unwrap_or_else(|| declaration.owner.clone()),
                    )
                    .or_default()
                    .constants
                    .insert(domain_const.symbol);
            }
        }

        let mut cache = Self::default();
        for (path, requested) in requests {
            let loaded;
            let text = if let Some(text) = inline_sources.get(&path) {
                *text
            } else {
                loaded = match read_repo_utf8(root, &path, "identity schema reference") {
                    Ok(text) => text,
                    Err(detail) => {
                        for symbol in requested.tests {
                            cache
                                .tests
                                .insert((path.clone(), symbol), Err(detail.clone()));
                        }
                        for symbol in requested.functions {
                            cache
                                .functions
                                .insert((path.clone(), symbol), Err(detail.clone()));
                        }
                        for symbol in requested.constants {
                            cache
                                .constants
                                .insert((path.clone(), symbol), Err(detail.clone()));
                        }
                        continue;
                    }
                };
                loaded.as_str()
            };
            let index = (!requested.tests.is_empty() || !requested.functions.is_empty())
                .then(|| RustSourceIndex::new(text));
            for symbol in requested.tests {
                let result = identity_test_bytes_with_index(
                    text,
                    index.as_ref().expect("test requests build a Rust index"),
                    &path,
                    &symbol,
                );
                cache.tests.insert((path.clone(), symbol), result);
            }
            for symbol in requested.functions {
                let index = index
                    .as_ref()
                    .expect("schema function requests build a Rust index");
                let result =
                    normalized_rust_function_closure_with_index(text, index, [symbol.clone()])
                        .map_err(|detail| {
                            format!("schema function {path}#{symbol} closure is invalid: {detail}")
                        });
                cache.functions.insert((path.clone(), symbol), result);
            }
            let module_scopes = index
                .as_ref()
                .map(|index| index.module_scopes.clone())
                .unwrap_or_else(|| rust_owner_scopes(text, "mod", rust_module_owner));
            for symbol in requested.constants {
                let declarations =
                    runtime_const_declarations_with_scopes(text, &symbol, &module_scopes);
                let result = match declarations.as_slice() {
                    [declaration] => {
                        let normalized = normalized_rust_fragment(declaration);
                        if !normalized.contains(&b'=') {
                            Err(format!(
                                "schema constant {path}#{symbol} declares no exact value"
                            ))
                        } else {
                            let start = declaration.as_ptr() as usize - text.as_ptr() as usize;
                            Ok(IdentityConstantReference {
                                locator: rust_item_locator(
                                    &path,
                                    &module_scopes,
                                    start,
                                    "const",
                                    &symbol,
                                ),
                                declaration: normalized,
                            })
                        }
                    }
                    declarations => Err(format!(
                        "schema constant {path}#{symbol} must resolve to exactly one const declaration; found {}",
                        declarations.len()
                    )),
                };
                cache.constants.insert((path.clone(), symbol), result);
            }
        }
        cache
    }

    fn test(&self, path: &str, symbol: &str) -> Result<&[u8], String> {
        cached_reference(&self.tests, path, symbol, "identity evidence")
    }

    fn function(&self, path: &str, symbol: &str) -> Result<&[u8], String> {
        cached_reference(&self.functions, path, symbol, "schema function")
    }

    fn constant(&self, path: &str, symbol: &str) -> Result<&IdentityConstantReference, String> {
        match self.constants.get(&(path.to_string(), symbol.to_string())) {
            Some(Ok(reference)) => Ok(reference),
            Some(Err(detail)) => Err(detail.clone()),
            None => Err(format!("uncached schema constant {path}#{symbol}")),
        }
    }
}

fn cached_reference<'a>(
    cache: &'a BTreeMap<(String, String), Result<Vec<u8>, String>>,
    path: &str,
    symbol: &str,
    purpose: &str,
) -> Result<&'a [u8], String> {
    match cache.get(&(path.to_string(), symbol.to_string())) {
        Some(Ok(bytes)) => Ok(bytes),
        Some(Err(detail)) => Err(detail.clone()),
        None => Err(format!("uncached {purpose} {path}#{symbol}")),
    }
}

fn normalized_string_constant_value(bytes: &[u8]) -> Option<&str> {
    let declaration = std::str::from_utf8(bytes).ok()?;
    let (left, value) = declaration.split_once('=')?;
    left.ends_with(":&str").then_some(())?;
    value
        .strip_suffix(';')?
        .strip_prefix('"')?
        .strip_suffix('"')
}

fn fingerprint_part(payload: &mut Vec<u8>, bytes: &[u8]) {
    let len = u64::try_from(bytes.len()).expect("Rust slice length fits the u64 registry frame");
    payload.extend_from_slice(&len.to_le_bytes());
    payload.extend_from_slice(bytes);
}

fn schema_fingerprint_digest(payload: &[u8]) -> [u8; 32] {
    *fs_blake3::hash_domain(SCHEMA_FINGERPRINT_DOMAIN, payload).as_bytes()
}

fn schema_fingerprint_hex(digest: &[u8; 32]) -> String {
    use std::fmt::Write as _;

    let mut hex = String::with_capacity(64);
    for byte in digest {
        let _ = write!(hex, "{byte:02x}");
    }
    hex
}

fn canonical_schema_fingerprint(value: &str, version: u32, digest_bytes: usize) -> bool {
    let Some(digest) = value.strip_prefix(&format!("v{version}-")) else {
        return false;
    };
    digest.len() == digest_bytes * 2
        && digest
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
}

fn identity_schema_base_hash_with_index(
    root: &Path,
    decl: &IdentityDecl,
    text: &str,
    index: &RustSourceIndex<'_>,
    references: &IdentityReferenceCache,
    exemptions: &[IdentityExemption],
) -> Result<[u8; 32], String> {
    let mut payload = Vec::new();
    fingerprint_part(&mut payload, b"owner-path");
    fingerprint_part(&mut payload, decl.owner.as_bytes());
    for part in [
        decl.domain.as_bytes(),
        decl.encoder.as_bytes(),
        decl.digest.as_bytes(),
        decl.encoding.as_bytes(),
    ] {
        fingerprint_part(&mut payload, part);
    }
    let domain_const = decl.domain_const.as_deref().ok_or_else(|| {
        format!(
            "identity {}: domain_const is required for schema fingerprinting",
            decl.id
        )
    })?;
    let domain_constant = parse_schema_constant_reference(domain_const).map_err(|detail| {
        format!(
            "identity {}: domain_const {domain_const:?} is invalid: {detail}",
            decl.id
        )
    })?;
    fingerprint_part(&mut payload, b"domain-const");
    fingerprint_part(&mut payload, domain_constant.canonical().as_bytes());
    let domain_path = domain_constant.path.as_deref().unwrap_or(&decl.owner);
    let domain_reference = references.constant(domain_path, &domain_constant.symbol)?;
    fingerprint_part(&mut payload, domain_reference.locator.as_bytes());
    fingerprint_part(&mut payload, &domain_reference.declaration);
    for value in decl.sources.iter().chain(decl.semantic_fields.iter()) {
        fingerprint_part(&mut payload, value.as_bytes());
    }
    for helper in &decl.encoder_helpers {
        fingerprint_part(&mut payload, helper.as_bytes());
    }
    for field in &decl.source_fields {
        fingerprint_part(&mut payload, field.qualified.as_bytes());
        fingerprint_part(&mut payload, field.class.name().as_bytes());
        fingerprint_part(
            &mut payload,
            field.reason.as_deref().unwrap_or("").as_bytes(),
        );
    }
    for binding in &decl.source_bindings {
        fingerprint_part(&mut payload, binding.source_field.as_bytes());
        for semantic in &binding.semantic_fields {
            fingerprint_part(&mut payload, semantic.as_bytes());
        }
    }
    for semantic in &decl.external_semantic_fields {
        fingerprint_part(&mut payload, semantic.as_bytes());
    }
    for (field, reason) in &decl.excluded_fields {
        fingerprint_part(&mut payload, field.as_bytes());
        fingerprint_part(&mut payload, reason.as_bytes());
    }
    for source in &decl.sources {
        let structs = braced_declaration_fragments(text, &format!("struct {source}"))
            .into_iter()
            .filter(|fragment| {
                let start = fragment.as_ptr() as usize - text.as_ptr() as usize;
                rust_item_is_runtime_active_with_scopes(text, start, &index.module_scopes)
            })
            .collect::<Vec<_>>();
        let enums = braced_declaration_fragments(text, &format!("enum {source}"))
            .into_iter()
            .filter(|fragment| {
                let start = fragment.as_ptr() as usize - text.as_ptr() as usize;
                rust_item_is_runtime_active_with_scopes(text, start, &index.module_scopes)
            })
            .collect::<Vec<_>>();
        let (kind, declaration) = match (structs.as_slice(), enums.as_slice()) {
            ([declaration], []) => ("struct", *declaration),
            ([], [declaration]) => ("enum", *declaration),
            _ => {
                return Err(format!(
                    "identity {}: source {source:?} must resolve to exactly one struct or enum declaration",
                    decl.id
                ));
            }
        };
        let declaration_start = declaration.as_ptr() as usize - text.as_ptr() as usize;
        fingerprint_part(&mut payload, source.as_bytes());
        fingerprint_part(
            &mut payload,
            rust_item_locator(
                &decl.owner,
                &index.module_scopes,
                declaration_start,
                kind,
                source,
            )
            .as_bytes(),
        );
        fingerprint_part(&mut payload, &normalized_rust_fragment(declaration));
    }
    let mut evidence_tests = BTreeSet::new();
    for mutation in decl
        .mutations
        .iter()
        .chain(decl.nonsemantic_mutations.iter())
    {
        fingerprint_part(&mut payload, mutation.field.as_bytes());
        fingerprint_part(&mut payload, mutation.test_path.as_bytes());
        fingerprint_part(&mut payload, mutation.test_symbol.as_bytes());
        evidence_tests.insert((mutation.test_path.as_str(), mutation.test_symbol.as_str()));
    }
    for (path, symbol) in evidence_tests {
        fingerprint_part(&mut payload, references.test(path, symbol)?);
    }
    let (version_guard_path, version_guard_symbol) =
        decl.version_guard.split_once('#').ok_or_else(|| {
            format!(
                "identity {}: version_guard must be path#test_symbol",
                decl.id
            )
        })?;
    fingerprint_part(&mut payload, decl.version_guard.as_bytes());
    fingerprint_part(
        &mut payload,
        references.test(version_guard_path, version_guard_symbol)?,
    );
    let owner_roots = std::iter::once(decl.encoder.clone())
        .chain(std::iter::once(decl.transport_guard.clone()))
        .chain(decl.encoder_helpers.iter().cloned())
        .chain(
            decl.schema_functions
                .iter()
                .filter(|function| function.path.is_none())
                .map(|function| function.symbol.clone()),
        )
        .collect::<BTreeSet<_>>();
    let owner_closure = normalized_rust_function_closure_with_index(text, index, owner_roots)
        .map_err(|detail| {
            format!(
                "identity {}: owner-local function closure is invalid: {detail}",
                decl.id
            )
        })?;
    fingerprint_part(&mut payload, &owner_closure);
    for function in decl
        .schema_functions
        .iter()
        .filter(|function| function.path.is_some())
    {
        fingerprint_part(&mut payload, function.canonical().as_bytes());
        let path = function.path.as_deref().expect("filtered path is present");
        fingerprint_part(&mut payload, references.function(path, &function.symbol)?);
    }
    for constant in &decl.schema_constants {
        fingerprint_part(&mut payload, constant.canonical().as_bytes());
        let path = constant.path.as_deref().unwrap_or(&decl.owner);
        let reference = references.constant(path, &constant.symbol)?;
        fingerprint_part(&mut payload, reference.locator.as_bytes());
        fingerprint_part(&mut payload, &reference.declaration);
    }
    for exemption in exemptions
        .iter()
        .filter(|exemption| exemption.covered_by == decl.id)
    {
        fingerprint_part(&mut payload, b"covered-exemption");
        fingerprint_part(&mut payload, exemption.path.as_bytes());
        fingerprint_part(&mut payload, exemption.symbol.as_bytes());
        fingerprint_part(&mut payload, exemption.reason.as_bytes());
        fingerprint_part(&mut payload, exemption.covered_by.as_bytes());
        fingerprint_part(
            &mut payload,
            &external_exemption_schema_bytes(root, exemption)?,
        );
    }
    Ok(schema_fingerprint_digest(&payload))
}

#[cfg(test)]
fn identity_schema_base_hash(
    root: &Path,
    decl: &IdentityDecl,
    text: &str,
    exemptions: &[IdentityExemption],
) -> Result<[u8; 32], String> {
    let index = RustSourceIndex::new(text);
    let inline_sources = BTreeMap::from([(decl.owner.clone(), text)]);
    let references = IdentityReferenceCache::build(root, std::iter::once(decl), &inline_sources);
    identity_schema_base_hash_with_index(root, decl, text, &index, &references, exemptions)
}

fn resolve_schema_fingerprint(
    index: usize,
    declarations: &[IdentityDecl],
    indices: &BTreeMap<String, usize>,
    states: &mut [u8],
    resolved: &mut [Option<[u8; 32]>],
    stack: &mut Vec<usize>,
) -> Result<[u8; 32], String> {
    if states[index] == 2 {
        return resolved[index].ok_or_else(|| {
            format!(
                "identity {}: resolved dependency fingerprint is unavailable",
                declarations[index].id
            )
        });
    }
    if states[index] == 1 {
        let cycle_start = stack
            .iter()
            .position(|candidate| *candidate == index)
            .unwrap_or(0);
        let mut cycle = stack[cycle_start..]
            .iter()
            .map(|candidate| declarations[*candidate].id.clone())
            .collect::<Vec<_>>();
        cycle.push(declarations[index].id.clone());
        return Err(format!(
            "schema dependency cycle is forbidden: {}",
            cycle.join(" -> ")
        ));
    }

    states[index] = 1;
    stack.push(index);
    let declaration = &declarations[index];
    let base_hash = declaration.schema_base_hash.ok_or_else(|| {
        format!(
            "identity {}: schema fingerprint base is unavailable",
            declaration.id
        )
    })?;
    let mut payload = Vec::new();
    fingerprint_part(&mut payload, b"owner-schema-base");
    fingerprint_part(&mut payload, &base_hash);
    let mut dependencies = declaration.schema_dependencies.clone();
    dependencies.sort();
    for dependency in dependencies {
        let dependency_index = *indices.get(&dependency).ok_or_else(|| {
            format!(
                "identity {}: schema dependency {dependency:?} does not exist",
                declaration.id
            )
        })?;
        if dependency_index == index {
            return Err(format!(
                "identity {}: schema self-dependency is forbidden",
                declaration.id
            ));
        }
        let dependency_hash = resolve_schema_fingerprint(
            dependency_index,
            declarations,
            indices,
            states,
            resolved,
            stack,
        )?;
        fingerprint_part(&mut payload, dependency.as_bytes());
        fingerprint_part(
            &mut payload,
            &declarations[dependency_index].version.to_le_bytes(),
        );
        fingerprint_part(&mut payload, &dependency_hash);
    }
    let hash = schema_fingerprint_digest(&payload);
    let popped = stack.pop();
    debug_assert_eq!(popped, Some(index));
    states[index] = 2;
    resolved[index] = Some(hash);
    Ok(hash)
}

fn resolve_schema_fingerprints(declarations: &mut [IdentityDecl]) -> Vec<Violation> {
    for declaration in declarations.iter_mut() {
        declaration.schema_fingerprint.clear();
    }
    let mut indices = BTreeMap::new();
    for (index, declaration) in declarations.iter().enumerate() {
        if indices.insert(declaration.id.clone(), index).is_some() {
            return vec![identity_violation(
                "<repo>",
                format!(
                    "duplicate semantic identity id {:?} prevents dependency resolution",
                    declaration.id
                ),
            )];
        }
    }

    let mut preflight = Vec::new();
    for declaration in declarations.iter() {
        if declaration.schema_base_hash.is_none() {
            preflight.push(identity_violation(
                &declaration.owner,
                format!(
                    "identity {}: schema fingerprint base is unavailable",
                    declaration.id
                ),
            ));
        }
        for dependency in &declaration.schema_dependencies {
            if dependency == &declaration.id {
                preflight.push(identity_violation(
                    &declaration.owner,
                    format!(
                        "identity {}: schema self-dependency is forbidden",
                        declaration.id
                    ),
                ));
            } else if !indices.contains_key(dependency) {
                preflight.push(identity_violation(
                    &declaration.owner,
                    format!(
                        "identity {}: schema dependency {dependency:?} does not exist",
                        declaration.id
                    ),
                ));
            }
        }
    }
    if !preflight.is_empty() {
        return preflight;
    }

    let mut states = vec![0; declarations.len()];
    let mut resolved = vec![None; declarations.len()];
    let mut stack = Vec::new();
    for index in indices.values().copied() {
        if let Err(detail) = resolve_schema_fingerprint(
            index,
            declarations,
            &indices,
            &mut states,
            &mut resolved,
            &mut stack,
        ) {
            return vec![identity_violation("<repo>", detail)];
        }
    }
    for (declaration, hash) in declarations.iter_mut().zip(resolved) {
        let hash = hash.expect("each declaration was resolved by the complete DFS pass");
        declaration.schema_fingerprint =
            format!("v{}-{}", declaration.version, schema_fingerprint_hex(&hash));
    }
    Vec::new()
}

fn strict_json_object<'a>(
    value: &'a JsonValue,
    context: &str,
) -> Result<&'a BTreeMap<String, JsonValue>, String> {
    if let JsonValue::Object(object) = value {
        Ok(object)
    } else {
        Err(format!("{context} must be a JSON object"))
    }
}

fn strict_json_array<'a>(value: &'a JsonValue, context: &str) -> Result<&'a [JsonValue], String> {
    if let JsonValue::Array(values) = value {
        Ok(values)
    } else {
        Err(format!("{context} must be a JSON array"))
    }
}

fn strict_json_string<'a>(value: &'a JsonValue, context: &str) -> Result<&'a str, String> {
    if let JsonValue::String(value) = value {
        Ok(value)
    } else {
        Err(format!("{context} must be a JSON string"))
    }
}

fn strict_json_u32(value: &JsonValue, context: &str) -> Result<u32, String> {
    if let JsonValue::Number(value) = value
        && value.bytes().all(|byte| byte.is_ascii_digit())
    {
        return value
            .parse::<u32>()
            .map_err(|_| format!("{context} must be an unsigned 32-bit integer"));
    }
    Err(format!("{context} must be an unsigned 32-bit JSON integer"))
}

fn strict_json_keys(
    object: &BTreeMap<String, JsonValue>,
    expected: &[&str],
    context: &str,
) -> Result<(), String> {
    let actual = object.keys().map(String::as_str).collect::<BTreeSet<_>>();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    if actual == expected {
        Ok(())
    } else {
        Err(format!(
            "{context} has noncanonical keys: expected {expected:?}, found {actual:?}"
        ))
    }
}

fn strict_json_field<'a>(
    object: &'a BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<&'a JsonValue, String> {
    object
        .get(key)
        .ok_or_else(|| format!("{context} is missing field {key:?}"))
}

fn authority_symbol_is_canonical(symbol: &str) -> bool {
    symbol == "<script>" || canonical_function_reference(symbol)
}

#[allow(clippy::too_many_lines)] // Strict row validation stays beside the tiny manifest parser.
fn load_authority_manifest(root: &Path) -> Result<AuthorityManifest, Vec<Violation>> {
    let text = read_repo_utf8(root, AUTHORITY_FILE, "authority manifest").map_err(|detail| {
        vec![identity_violation(
            "<repo>",
            format!("{AUTHORITY_FILE} is invalid: {detail}"),
        )]
    })?;

    let mut violations = Vec::new();
    let mut required_ids = BTreeSet::new();
    let mut external_owners = Vec::new();
    let mut exemptions = Vec::new();
    let parsed = JsonParser::new(&text).finish().map_err(|detail| {
        vec![identity_violation(
            "<repo>",
            format!("{AUTHORITY_FILE} is not strict JSON: {detail}"),
        )]
    })?;
    let object = strict_json_object(&parsed, AUTHORITY_FILE)
        .and_then(|object| {
            strict_json_keys(
                object,
                &["schema", "required_ids", "external_owners", "exemptions"],
                AUTHORITY_FILE,
            )?;
            Ok(object)
        })
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    let schema = strict_json_field(object, "schema", AUTHORITY_FILE)
        .and_then(|value| strict_json_string(value, "authority manifest schema"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    if schema != AUTHORITY_SCHEMA {
        violations.push(identity_violation(
            "<repo>",
            format!("{AUTHORITY_FILE} must declare schema {AUTHORITY_SCHEMA:?}"),
        ));
    }

    let required_rows = strict_json_field(object, "required_ids", AUTHORITY_FILE)
        .and_then(|value| strict_json_array(value, "authority required_ids"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    for (index, value) in required_rows.iter().enumerate() {
        let context = format!("{AUTHORITY_FILE} required_ids row {}", index + 1);
        let row = match strict_json_object(value, &context).and_then(|row| {
            strict_json_keys(row, &["id"], &context)?;
            Ok(row)
        }) {
            Ok(row) => row,
            Err(detail) => {
                violations.push(identity_violation("<repo>", detail));
                continue;
            }
        };
        let id = match strict_json_field(row, "id", &context)
            .and_then(|value| strict_json_string(value, &format!("{context} id")))
        {
            Ok(id) => id,
            Err(detail) => {
                violations.push(identity_violation("<repo>", detail));
                continue;
            }
        };
        if !id.contains(':') || id.chars().any(char::is_whitespace) {
            violations.push(identity_violation(
                "<repo>",
                format!("{context} id {id:?} is not canonical"),
            ));
        } else if !required_ids.insert(id.to_string()) {
            violations.push(identity_violation(
                "<repo>",
                format!("{context} duplicates required id {id:?}"),
            ));
        }
    }

    let external_rows = strict_json_field(object, "external_owners", AUTHORITY_FILE)
        .and_then(|value| strict_json_array(value, "authority external_owners"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    for (index, value) in external_rows.iter().enumerate() {
        let context = format!("{AUTHORITY_FILE} external_owners row {}", index + 1);
        let parsed = (|| -> Result<ExternalOwner, String> {
            let row = strict_json_object(value, &context)?;
            strict_json_keys(
                row,
                &["id", "path", "symbol", "version", "domain"],
                &context,
            )?;
            let string = |key: &str| {
                strict_json_field(row, key, &context)
                    .and_then(|value| strict_json_string(value, &format!("{context} {key}")))
            };
            let id = string("id")?;
            let path = string("path")?;
            let symbol = string("symbol")?;
            let domain = string("domain")?;
            let version = strict_json_field(row, "version", &context)
                .and_then(|value| strict_json_u32(value, &format!("{context} version")))?;
            if !id.contains(':') || id.chars().any(char::is_whitespace) {
                return Err(format!("{context} id {id:?} is not canonical"));
            }
            if !safe_relative(path) || !authority_symbol_is_canonical(symbol) {
                return Err(format!(
                    "{context} has unsafe/noncanonical target {path}#{symbol}"
                ));
            }
            if !domain_carries_version(domain, version) {
                return Err(format!(
                    "{context} domain {domain:?} must carry exact version {version}"
                ));
            }
            Ok(ExternalOwner {
                id: id.to_string(),
                path: path.to_string(),
                symbol: symbol.to_string(),
                version,
                domain: domain.to_string(),
            })
        })();
        match parsed {
            Ok(owner) => external_owners.push(owner),
            Err(detail) => violations.push(identity_violation("<repo>", detail)),
        }
    }

    let exemption_rows = strict_json_field(object, "exemptions", AUTHORITY_FILE)
        .and_then(|value| strict_json_array(value, "authority exemptions"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    for (index, value) in exemption_rows.iter().enumerate() {
        let context = format!("{AUTHORITY_FILE} exemptions row {}", index + 1);
        let parsed = (|| -> Result<IdentityExemption, String> {
            let row = strict_json_object(value, &context)?;
            strict_json_keys(row, &["path", "symbol", "reason", "covered_by"], &context)?;
            let string = |key: &str| {
                strict_json_field(row, key, &context)
                    .and_then(|value| strict_json_string(value, &format!("{context} {key}")))
            };
            let path = string("path")?;
            let symbol = string("symbol")?;
            let reason = string("reason")?;
            let covered_by = string("covered_by")?;
            if !safe_relative(path) || !authority_symbol_is_canonical(symbol) {
                return Err(format!(
                    "{context} target {path}#{symbol} is unsafe/noncanonical"
                ));
            }
            if reason.is_empty()
                || (covered_by != "none"
                    && (!covered_by.contains(':') || covered_by.chars().any(char::is_whitespace)))
            {
                return Err(format!(
                    "{context} {path}#{symbol} has a noncanonical reason/covered_by pair"
                ));
            }
            Ok(IdentityExemption {
                path: path.to_string(),
                symbol: symbol.to_string(),
                reason: reason.to_string(),
                covered_by: covered_by.to_string(),
            })
        })();
        match parsed {
            Ok(exemption) => exemptions.push(exemption),
            Err(detail) => violations.push(identity_violation("<repo>", detail)),
        }
    }
    external_owners.sort_by(|left, right| {
        (&left.id, &left.path, &left.symbol).cmp(&(&right.id, &right.path, &right.symbol))
    });
    exemptions.sort_by(|left, right| (&left.path, &left.symbol).cmp(&(&right.path, &right.symbol)));
    for pair in external_owners.windows(2) {
        if pair[0] == pair[1] || pair[0].id == pair[1].id {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "{AUTHORITY_FILE} contains duplicate external owner {:?}",
                    pair[1].id
                ),
            ));
        }
    }
    for pair in exemptions.windows(2) {
        if pair[0].path == pair[1].path && pair[0].symbol == pair[1].symbol {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "{AUTHORITY_FILE} contains duplicate exemption {}#{}",
                    pair[1].path, pair[1].symbol
                ),
            ));
        }
    }
    if root.join(".git").exists() {
        let pinned = REQUIRED_IDENTITY_IDS
            .iter()
            .map(|id| (*id).to_string())
            .collect::<BTreeSet<_>>();
        for id in pinned.difference(&required_ids) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "{AUTHORITY_FILE} removed pinned identity authority {id:?}; retained authorities require an explicit migration and gate update"
                ),
            ));
        }
        for id in required_ids.difference(&pinned) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "{AUTHORITY_FILE} adds identity authority {id:?} without updating the independent REQUIRED_IDENTITY_IDS gate baseline"
                ),
            ));
        }
    }
    if violations.is_empty() {
        Ok(AuthorityManifest {
            required_ids,
            external_owners,
            exemptions,
        })
    } else {
        Err(violations)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScriptDialect {
    Shell,
    Python,
}

impl ScriptDialect {
    const fn name(self) -> &'static str {
        match self {
            Self::Shell => "shell",
            Self::Python => "Python",
        }
    }
}

fn interpreter_dialect(word: &str) -> Option<ScriptDialect> {
    let executable = word.rsplit('/').next().unwrap_or(word);
    if matches!(executable, "sh" | "bash" | "dash" | "ksh" | "zsh") {
        Some(ScriptDialect::Shell)
    } else if executable == "python"
        || executable == "python3"
        || executable.strip_prefix("python3.").is_some_and(|version| {
            !version.is_empty()
                && version
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || byte == b'.')
        })
    {
        Some(ScriptDialect::Python)
    } else {
        None
    }
}

fn script_dialect(path: &str, text: &str) -> Option<ScriptDialect> {
    match Path::new(path).extension().and_then(|value| value.to_str()) {
        Some("sh") => return Some(ScriptDialect::Shell),
        Some("py") => return Some(ScriptDialect::Python),
        _ => {}
    }
    let shebang = text.lines().next()?.strip_prefix("#!")?.trim();
    let mut words = shebang.split_ascii_whitespace();
    let first = words.next()?;
    if first.rsplit('/').next() != Some("env") {
        return interpreter_dialect(first);
    }
    for word in words {
        if word == "-S" || word.starts_with('-') || shell_assignment_word(word) {
            continue;
        }
        return interpreter_dialect(word);
    }
    None
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ScriptStructureKind {
    UnterminatedQuote {
        dialect: ScriptDialect,
        delimiter: char,
    },
    UnterminatedTripleString {
        delimiter: &'static str,
    },
    UnterminatedHeredoc {
        delimiter: String,
    },
    AmbiguousHeredocExecution {
        delimiter: String,
    },
    IncompleteFunctionScope {
        dialect: ScriptDialect,
        symbol: Option<String>,
        expected: &'static str,
    },
    IncompleteFStringReplacement,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ScriptStructureError {
    opening_offset: usize,
    kind: ScriptStructureKind,
}

impl ScriptStructureError {
    fn shifted(mut self, offset: usize) -> Self {
        self.opening_offset = self.opening_offset.saturating_add(offset);
        self
    }

    fn describe(&self) -> String {
        match &self.kind {
            ScriptStructureKind::UnterminatedQuote { dialect, delimiter } => format!(
                "unterminated {} quote opened with {delimiter:?}",
                dialect.name()
            ),
            ScriptStructureKind::UnterminatedTripleString { delimiter } => {
                format!("unterminated Python triple string opened with {delimiter:?}")
            }
            ScriptStructureKind::UnterminatedHeredoc { delimiter } => {
                format!("unterminated heredoc {delimiter:?}")
            }
            ScriptStructureKind::AmbiguousHeredocExecution { delimiter } => format!(
                "identity-bearing heredoc {delimiter:?} has no unambiguous shell or Python execution context"
            ),
            ScriptStructureKind::IncompleteFunctionScope {
                dialect,
                symbol,
                expected,
            } => format!(
                "incomplete {} function scope{}; expected {expected}",
                dialect.name(),
                symbol
                    .as_deref()
                    .map_or_else(String::new, |symbol| format!(" {symbol:?}"))
            ),
            ScriptStructureKind::IncompleteFStringReplacement => {
                "incomplete Python f-string replacement expression; expected a closing '}'"
                    .to_string()
            }
        }
    }
}

type ScriptResult<T> = Result<T, ScriptStructureError>;

#[derive(Debug, Clone, PartialEq, Eq)]
struct ScriptBlock<'a> {
    symbol: String,
    text: &'a str,
    dialect: ScriptDialect,
    start: usize,
}

fn script_structure_detail(path: &str, text: &str, error: &ScriptStructureError) -> String {
    let offset = error.opening_offset.min(text.len());
    let prefix = &text[..offset];
    let line = prefix.bytes().filter(|byte| *byte == b'\n').count() + 1;
    let column = prefix
        .rsplit('\n')
        .next()
        .map_or(1, |line| line.chars().count() + 1);
    format!(
        "identity discovery refused malformed script {path} at line {line}, column {column}: {}",
        error.describe()
    )
}

fn python_signature_colon(code: &str, opening: usize) -> Option<usize> {
    let bytes = code.as_bytes();
    let mut stack = vec![b')'];
    let mut index = opening + 1;
    while let Some(byte) = bytes.get(index).copied() {
        match byte {
            b'(' => stack.push(b')'),
            b'[' => stack.push(b']'),
            b'{' => stack.push(b'}'),
            b')' | b']' | b'}' => {
                if stack.pop() != Some(byte) {
                    return None;
                }
                if stack.is_empty() {
                    index += 1;
                    break;
                }
            }
            _ => {}
        }
        index += 1;
    }
    if !stack.is_empty() {
        return None;
    }
    while bytes.get(index).is_some_and(u8::is_ascii_whitespace) {
        index += 1;
    }
    if bytes.get(index..index + 2) == Some(b"->") {
        index += 2;
        let mut annotation = Vec::new();
        while let Some(byte) = bytes.get(index).copied() {
            match byte {
                b'(' => annotation.push(b')'),
                b'[' => annotation.push(b']'),
                b'{' => annotation.push(b'}'),
                b')' | b']' | b'}' => {
                    if annotation.pop() != Some(byte) {
                        return None;
                    }
                }
                b':' if annotation.is_empty() => return Some(index),
                _ => {}
            }
            index += 1;
        }
        None
    } else {
        (bytes.get(index) == Some(&b':')).then_some(index)
    }
}

fn python_decorator_start(
    projection_lines: &[&str],
    offsets: &[usize],
    definition_line: usize,
    indentation: usize,
) -> usize {
    let mut start_line = definition_line;
    let mut line_index = definition_line;
    while line_index > 0 {
        line_index -= 1;
        let line = projection_lines[line_index];
        let trimmed = line.trim_start();
        if trimmed.trim().is_empty() {
            break;
        }
        let candidate_indentation = line.len() - trimmed.len();
        if candidate_indentation < indentation {
            break;
        }
        if candidate_indentation == indentation && trimmed.starts_with('@') {
            start_line = line_index;
            continue;
        }
        if candidate_indentation == indentation
            && !trimmed
                .trim()
                .bytes()
                .all(|byte| matches!(byte, b')' | b']' | b'}' | b','))
        {
            break;
        }
    }
    offsets[start_line]
}

fn python_function_blocks<'a>(path: &str, text: &'a str) -> ScriptResult<Vec<ScriptBlock<'a>>> {
    let projection = executable_script_views(path, text)?.0;
    let lines = projection.split_inclusive('\n').collect::<Vec<_>>();
    let source_lines = text.split_inclusive('\n').collect::<Vec<_>>();
    let mut offsets = Vec::with_capacity(lines.len() + 1);
    let mut offset = 0_usize;
    for line in &lines {
        offsets.push(offset);
        offset += line.len();
    }
    offsets.push(offset);

    let mut blocks = Vec::new();
    for (line_index, line) in lines.iter().enumerate() {
        let trimmed = line.trim_start();
        let Some(rest) = trimmed
            .strip_prefix("def ")
            .or_else(|| trimmed.strip_prefix("async def "))
        else {
            continue;
        };
        let Some((symbol, _)) = rest.split_once('(') else {
            return Err(ScriptStructureError {
                opening_offset: offsets[line_index] + line.len() - trimmed.len(),
                kind: ScriptStructureKind::IncompleteFunctionScope {
                    dialect: ScriptDialect::Python,
                    symbol: None,
                    expected: "a canonical function name followed by '('",
                },
            });
        };
        let symbol = symbol.trim();
        if !canonical_symbol(symbol) {
            continue;
        }
        let definition_start = offsets[line_index] + line.len() - trimmed.len();
        let opening = definition_start
            + trimmed
                .find('(')
                .expect("a split_once match retains its opening parenthesis");
        let Some(colon) = python_signature_colon(&projection, opening) else {
            return Err(ScriptStructureError {
                opening_offset: definition_start,
                kind: ScriptStructureKind::IncompleteFunctionScope {
                    dialect: ScriptDialect::Python,
                    symbol: Some(symbol.to_string()),
                    expected: "a balanced closing ')' and ':' after any return annotation",
                },
            });
        };
        let indentation = line.len() - trimmed.len();
        let header_line = offsets[..lines.len()]
            .partition_point(|start| *start <= colon)
            .saturating_sub(1);
        let inline_suite = !projection[colon + 1..offsets[header_line + 1]]
            .trim()
            .is_empty();
        let mut end_line = header_line + 1;
        let mut block_end = offsets[header_line + 1];
        let mut has_indented_suite = false;
        while !inline_suite && end_line < lines.len() {
            let candidate = source_lines[end_line];
            let projected_candidate = lines[end_line];
            let candidate_trimmed = candidate.trim_start();
            let projected_trimmed = projected_candidate.trim_start();
            let candidate_indentation = candidate.len() - candidate_trimmed.len();
            if candidate_trimmed.starts_with('#') {
                if candidate_indentation <= indentation {
                    break;
                }
                block_end = offsets[end_line + 1];
                end_line += 1;
                continue;
            }
            if projected_trimmed.trim().is_empty() {
                end_line += 1;
                continue;
            }
            if candidate_indentation <= indentation {
                break;
            }
            has_indented_suite = true;
            block_end = offsets[end_line + 1];
            end_line += 1;
        }
        if !inline_suite && !has_indented_suite {
            return Err(ScriptStructureError {
                opening_offset: offsets[line_index] + line.len() - trimmed.len(),
                kind: ScriptStructureKind::IncompleteFunctionScope {
                    dialect: ScriptDialect::Python,
                    symbol: Some(symbol.to_string()),
                    expected: "an inline or indented function suite",
                },
            });
        }
        let block_start = python_decorator_start(&lines, &offsets, line_index, indentation);
        blocks.push(ScriptBlock {
            symbol: symbol.to_string(),
            text: &text[block_start..block_end],
            dialect: ScriptDialect::Python,
            start: block_start,
        });
    }
    Ok(blocks)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ShellHeredoc {
    delimiter: String,
    strip_tabs: bool,
    operator_start: usize,
    operator_end: usize,
}

fn shell_heredoc_word(command: &str, mut index: usize) -> Option<(String, usize)> {
    let bytes = command.as_bytes();
    let mut delimiter = Vec::new();
    let mut quote = None;
    let mut saw_word = false;
    while index < bytes.len() {
        let byte = bytes[index];
        if let Some(quote_delimiter) = quote {
            if byte == quote_delimiter {
                quote = None;
                saw_word = true;
                index += 1;
            } else if byte == b'\\' && quote_delimiter == b'"' {
                let escaped = *bytes.get(index + 1)?;
                if matches!(escaped, b'$' | b'`' | b'"' | b'\\' | b'\n') {
                    if escaped != b'\n' {
                        delimiter.push(escaped);
                    }
                } else {
                    delimiter.extend([byte, escaped]);
                }
                saw_word = true;
                index += 2;
            } else {
                delimiter.push(byte);
                saw_word = true;
                index += 1;
            }
            continue;
        }
        match byte {
            b'\'' | b'"' => {
                quote = Some(byte);
                saw_word = true;
                index += 1;
            }
            b'\\' => {
                let escaped = *bytes.get(index + 1)?;
                if escaped != b'\n' {
                    delimiter.push(escaped);
                }
                saw_word = true;
                index += 2;
            }
            byte if byte.is_ascii_whitespace()
                || matches!(byte, b';' | b'|' | b'&' | b'(' | b')' | b'<' | b'>') =>
            {
                break;
            }
            _ => {
                delimiter.push(byte);
                saw_word = true;
                index += 1;
            }
        }
    }
    if quote.is_some() || !saw_word {
        return None;
    }
    String::from_utf8(delimiter)
        .ok()
        .map(|delimiter| (delimiter, index))
}

fn shell_heredocs(command: &str) -> Vec<ShellHeredoc> {
    let bytes = command.as_bytes();
    let mut heredocs = Vec::new();
    let mut index = 0_usize;
    let mut quote = None;
    let mut arithmetic_depth = 0_usize;
    while index < bytes.len() {
        let byte = bytes[index];
        if let Some(delimiter) = quote {
            if byte == b'\\' && delimiter == b'"' {
                index = (index + 2).min(bytes.len());
            } else {
                index += 1;
                if byte == delimiter {
                    quote = None;
                }
            }
            continue;
        }
        if arithmetic_depth > 0 {
            match byte {
                b'\\' => index = (index + 2).min(bytes.len()),
                b'\'' | b'"' => {
                    quote = Some(byte);
                    index += 1;
                }
                b'(' => {
                    arithmetic_depth += 1;
                    index += 1;
                }
                b')' => {
                    arithmetic_depth -= 1;
                    index += 1;
                }
                _ => index += 1,
            }
            continue;
        }
        match byte {
            b'\\' => index = (index + 2).min(bytes.len()),
            b'\'' | b'"' => {
                quote = Some(byte);
                index += 1;
            }
            b'$' if bytes
                .get(index + 1..index + 3)
                .is_some_and(|candidate| candidate == b"((") =>
            {
                arithmetic_depth = 2;
                index += 3;
            }
            b'(' if bytes.get(index + 1) == Some(&b'(') => {
                arithmetic_depth = 2;
                index += 2;
            }
            b'#' if index == 0
                || bytes[index - 1].is_ascii_whitespace()
                || matches!(bytes[index - 1], b';' | b'|' | b'&' | b'(') =>
            {
                while index < bytes.len() && bytes[index] != b'\n' {
                    index += 1;
                }
            }
            b'<' if (index == 0 || bytes[index - 1] != b'<')
                && bytes.get(index + 1) == Some(&b'<')
                && bytes.get(index + 2) != Some(&b'<') =>
            {
                let operator_start = index;
                index += 2;
                let strip_tabs = bytes.get(index) == Some(&b'-');
                if strip_tabs {
                    index += 1;
                }
                while bytes
                    .get(index)
                    .is_some_and(|byte| matches!(byte, b' ' | b'\t'))
                {
                    index += 1;
                }
                let Some((delimiter, operator_end)) = shell_heredoc_word(command, index) else {
                    continue;
                };
                index = operator_end;
                heredocs.push(ShellHeredoc {
                    delimiter,
                    strip_tabs,
                    operator_start,
                    operator_end,
                });
            }
            _ => index += 1,
        }
    }
    heredocs
}

fn shell_function_signature(line: &str) -> Option<(String, Option<usize>)> {
    let trimmed = line.trim_start();
    let (symbol, rest) = if let Some(after_keyword) = trimmed.strip_prefix("function ") {
        let after_keyword = after_keyword.trim_start();
        let symbol_end = after_keyword
            .find(|character: char| {
                character.is_ascii_whitespace() || matches!(character, '(' | '{')
            })
            .unwrap_or(after_keyword.len());
        let symbol = &after_keyword[..symbol_end];
        let mut rest = after_keyword[symbol_end..].trim_start();
        if let Some(after_parens) = rest.strip_prefix("()") {
            rest = after_parens.trim_start();
        }
        (symbol, rest)
    } else {
        let (symbol, rest) = trimmed.split_once("()")?;
        (symbol.trim(), rest)
    };
    if !canonical_symbol(symbol) {
        return None;
    }
    if rest.trim().is_empty() {
        return Some((symbol.to_string(), None));
    }
    let brace = rest.find('{')?;
    if !rest[..brace].trim().is_empty() {
        return None;
    }
    let opening = rest.as_ptr() as usize - line.as_ptr() as usize + brace;
    Some((symbol.to_string(), Some(opening)))
}

fn shell_braced_block_end(text: &str, opening_brace: usize, symbol: &str) -> ScriptResult<usize> {
    let bytes = text.as_bytes();
    if bytes.get(opening_brace) != Some(&b'{') {
        return Err(ScriptStructureError {
            opening_offset: opening_brace.min(text.len()),
            kind: ScriptStructureKind::IncompleteFunctionScope {
                dialect: ScriptDialect::Shell,
                symbol: Some(symbol.to_string()),
                expected: "an opening '{'",
            },
        });
    }
    let mut index = opening_brace;
    let mut depth = 0_usize;
    let mut quote = None::<(u8, usize)>;
    let mut in_comment = false;
    let mut at_line_start = false;
    let mut heredocs = VecDeque::<(ShellHeredoc, usize)>::new();
    let mut logical_command_start = opening_brace;
    let mut scope_closed = false;
    while index < bytes.len() {
        if at_line_start && let Some((heredoc, _)) = heredocs.front() {
            let line_end = text[index..]
                .find('\n')
                .map_or(text.len(), |relative| index + relative);
            let mut line = &text[index..line_end];
            if heredoc.strip_tabs {
                line = line.trim_start_matches('\t');
            }
            index = (line_end + usize::from(line_end < text.len())).min(text.len());
            if line.trim_end_matches('\r') == heredoc.delimiter {
                heredocs.pop_front();
            }
            at_line_start = true;
            if heredocs.is_empty() {
                logical_command_start = index;
            }
            if heredocs.is_empty() && scope_closed {
                return Ok(index);
            }
            continue;
        }

        let byte = bytes[index];
        if in_comment {
            index += 1;
            if byte == b'\n' {
                in_comment = false;
                let logical_command = &text[logical_command_start..index - 1];
                if !shell_line_continues(logical_command) {
                    heredocs.extend(shell_heredocs(logical_command).into_iter().map(|heredoc| {
                        let opening = logical_command_start + heredoc.operator_start;
                        (heredoc, opening)
                    }));
                    logical_command_start = index;
                    if heredocs.is_empty() && scope_closed {
                        return Ok(index - 1);
                    }
                }
                at_line_start = true;
            }
            continue;
        }
        if let Some((delimiter, _)) = quote {
            if byte == b'\\' && delimiter == b'"' {
                index = (index + 2).min(bytes.len());
            } else {
                index += 1;
                if byte == delimiter {
                    quote = None;
                }
            }
            continue;
        }
        match byte {
            b'\\' => index = (index + 2).min(bytes.len()),
            b'\'' | b'"' => {
                quote = Some((byte, index));
                index += 1;
            }
            b'#' if index == 0
                || bytes[index - 1].is_ascii_whitespace()
                || matches!(bytes[index - 1], b';' | b'|' | b'&' | b'(') =>
            {
                in_comment = true;
                index += 1;
            }
            b'{' if !scope_closed => {
                depth += 1;
                index += 1;
            }
            b'}' if !scope_closed => {
                if depth == 0 {
                    return Err(ScriptStructureError {
                        opening_offset: index,
                        kind: ScriptStructureKind::IncompleteFunctionScope {
                            dialect: ScriptDialect::Shell,
                            symbol: Some(symbol.to_string()),
                            expected: "a balanced function body",
                        },
                    });
                }
                depth -= 1;
                index += 1;
                if depth == 0 {
                    let line_end = text[index..]
                        .find('\n')
                        .map_or(text.len(), |relative| index + relative);
                    let logical_command = &text[logical_command_start..line_end];
                    if shell_heredocs(logical_command).is_empty() {
                        return Ok(index);
                    }
                    scope_closed = true;
                }
            }
            b'\n' => {
                let logical_command = &text[logical_command_start..index];
                if !shell_line_continues(logical_command) {
                    heredocs.extend(shell_heredocs(logical_command).into_iter().map(|heredoc| {
                        let opening = logical_command_start + heredoc.operator_start;
                        (heredoc, opening)
                    }));
                    logical_command_start = index + 1;
                    if heredocs.is_empty() && scope_closed {
                        return Ok(index);
                    }
                }
                at_line_start = true;
                index += 1;
            }
            _ => {
                at_line_start = false;
                index += 1;
            }
        }
    }
    if heredocs.is_empty() {
        let logical_command = &text[logical_command_start..];
        heredocs.extend(shell_heredocs(logical_command).into_iter().map(|heredoc| {
            let opening = logical_command_start + heredoc.operator_start;
            (heredoc, opening)
        }));
    }
    if let Some((heredoc, opening_offset)) = heredocs.front() {
        return Err(ScriptStructureError {
            opening_offset: *opening_offset,
            kind: ScriptStructureKind::UnterminatedHeredoc {
                delimiter: heredoc.delimiter.clone(),
            },
        });
    }
    if let Some((delimiter, opening_offset)) = quote {
        return Err(ScriptStructureError {
            opening_offset,
            kind: ScriptStructureKind::UnterminatedQuote {
                dialect: ScriptDialect::Shell,
                delimiter: char::from(delimiter),
            },
        });
    }
    Err(ScriptStructureError {
        opening_offset: opening_brace,
        kind: ScriptStructureKind::IncompleteFunctionScope {
            dialect: ScriptDialect::Shell,
            symbol: Some(symbol.to_string()),
            expected: "a closing '}'",
        },
    })
}

fn shell_function_blocks(text: &str) -> ScriptResult<Vec<(String, &str)>> {
    let heredoc_projection = mask_non_executable_heredocs(text)?;
    let projection = script_views(ScriptDialect::Shell, &heredoc_projection.masked)?.0;
    let mut blocks = Vec::new();
    for embedded in heredoc_projection
        .embedded
        .iter()
        .filter(|embedded| embedded.dialect == ScriptDialect::Shell)
    {
        let end = embedded.start + embedded.code.len();
        blocks.extend(shell_function_blocks(&text[embedded.start..end])?);
    }
    let mut offset = 0_usize;
    while offset < text.len() {
        let line_end = text[offset..]
            .find('\n')
            .map_or(text.len(), |relative| offset + relative);
        let projected_line = &projection[offset..line_end];
        let next_line = (line_end + usize::from(line_end < text.len())).min(text.len());
        if let Some((symbol, opening_in_line)) = shell_function_signature(projected_line) {
            let opening = opening_in_line.map(|opening| offset + opening).or_else(|| {
                let mut candidate_start = next_line;
                while candidate_start < text.len() {
                    let candidate_end = text[candidate_start..]
                        .find('\n')
                        .map_or(text.len(), |relative| candidate_start + relative);
                    let candidate = &projection[candidate_start..candidate_end];
                    let trimmed = candidate.trim_start();
                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        candidate_start = (candidate_end + usize::from(candidate_end < text.len()))
                            .min(text.len());
                        continue;
                    }
                    let after = trimmed.strip_prefix('{')?.trim_start();
                    if !after.is_empty() && !after.starts_with('#') {
                        return None;
                    }
                    return Some(candidate_start + candidate.len() - trimmed.len());
                }
                None
            });
            let Some(opening) = opening else {
                return Err(ScriptStructureError {
                    opening_offset: offset + projected_line.len()
                        - projected_line.trim_start().len(),
                    kind: ScriptStructureKind::IncompleteFunctionScope {
                        dialect: ScriptDialect::Shell,
                        symbol: Some(symbol),
                        expected: "an opening '{'",
                    },
                });
            };
            let end = shell_braced_block_end(text, opening, &symbol)?;
            blocks.push((symbol, &text[offset..end]));
            // Continue through the validated block so executable shell heredocs
            // and nested shell function scopes retain their own exact owners.
            offset = next_line;
            continue;
        }
        offset = next_line;
    }
    blocks.sort_by(|left, right| {
        (left.1.as_ptr() as usize, &left.0).cmp(&(right.1.as_ptr() as usize, &right.0))
    });
    blocks.dedup_by(|left, right| left.1.as_ptr() == right.1.as_ptr() && left.0 == right.0);
    Ok(blocks)
}

fn script_has_function_signature(path: &str, text: &str) -> ScriptResult<bool> {
    let dialect = script_dialect(path, text);
    if dialect == Some(ScriptDialect::Python) {
        return Ok(!python_function_blocks(path, text)?.is_empty());
    }
    let projection = executable_script_views(path, text)?.0;
    Ok(projection
        .lines()
        .any(|line| shell_function_signature(line).is_some())
        || dialect == Some(ScriptDialect::Shell) && !python_function_blocks(path, text)?.is_empty())
}

fn primary_script_blocks<'a>(
    path: &str,
    text: &'a str,
    symbol: &str,
) -> ScriptResult<Vec<ScriptBlock<'a>>> {
    if script_dialect(path, text) == Some(ScriptDialect::Python) {
        return Ok(python_function_blocks(path, text)?
            .into_iter()
            .filter(|block| block.symbol == symbol)
            .collect());
    }
    Ok(shell_function_blocks(text)?
        .into_iter()
        .filter_map(|(candidate, block)| {
            (candidate == symbol).then(|| ScriptBlock {
                symbol: candidate,
                text: block,
                dialect: ScriptDialect::Shell,
                start: block.as_ptr() as usize - text.as_ptr() as usize,
            })
        })
        .collect())
}

fn script_symbol_blocks<'a>(
    path: &str,
    text: &'a str,
    symbol: &str,
) -> ScriptResult<Vec<ScriptBlock<'a>>> {
    let mut blocks = Vec::new();
    let dialect = script_dialect(path, text);
    if dialect == Some(ScriptDialect::Shell) {
        blocks.extend(
            shell_function_blocks(text)?
                .into_iter()
                .filter_map(|(candidate, block)| {
                    (candidate == symbol).then(|| ScriptBlock {
                        symbol: candidate,
                        text: block,
                        dialect: ScriptDialect::Shell,
                        start: block.as_ptr() as usize - text.as_ptr() as usize,
                    })
                }),
        );
    }
    if dialect.is_some() {
        blocks.extend(
            python_function_blocks(path, text)?
                .into_iter()
                .filter(|block| block.symbol == symbol),
        );
    }
    Ok(blocks)
}

fn script_block_views(_path: &str, block: &ScriptBlock<'_>) -> ScriptResult<(String, String)> {
    let views = match block.dialect {
        ScriptDialect::Shell => shell_executable_script_views(block.text),
        ScriptDialect::Python => script_views(ScriptDialect::Python, block.text),
    };
    views.map_err(|error| error.shifted(block.start))
}

fn external_exemption_schema_bytes(
    root: &Path,
    exemption: &IdentityExemption,
) -> Result<Vec<u8>, String> {
    let path = checked_repo_file(root, &exemption.path, "external exemption source")?;
    let text = read_repo_utf8(root, &exemption.path, "external exemption source")?;
    if exemption.symbol == "<script>" {
        return Ok(text.into_bytes());
    }
    if path.extension().is_some_and(|extension| extension == "rs") {
        return normalized_rust_function_closure(&text, [exemption.symbol.clone()]).map_err(
            |detail| {
                format!(
                    "external exemption {}#{} has an invalid function closure: {detail}",
                    exemption.path, exemption.symbol
                )
            },
        );
    }
    let blocks = script_symbol_blocks(&exemption.path, &text, &exemption.symbol)
        .map_err(|error| script_structure_detail(&exemption.path, &text, &error))?;
    let [block] = blocks.as_slice() else {
        return Err(format!(
            "external exemption {}#{} resolves to {} script blocks",
            exemption.path,
            exemption.symbol,
            blocks.len()
        ));
    };
    Ok(block.text.as_bytes().to_vec())
}

fn exemption_schema_fingerprint(
    root: &Path,
    exemption: &IdentityExemption,
) -> Result<String, String> {
    let mut payload = Vec::new();
    let schema_bytes = external_exemption_schema_bytes(root, exemption)?;
    for part in [
        exemption.path.as_bytes(),
        exemption.symbol.as_bytes(),
        exemption.reason.as_bytes(),
        exemption.covered_by.as_bytes(),
        schema_bytes.as_slice(),
    ] {
        fingerprint_part(&mut payload, part);
    }
    Ok(format!(
        "v1-{}",
        schema_fingerprint_hex(&schema_fingerprint_digest(&payload))
    ))
}

fn external_owner_schema_fingerprint(
    root: &Path,
    owner: &ExternalOwner,
    exemptions: &[IdentityExemption],
) -> Result<String, String> {
    let path = checked_repo_file(root, &owner.path, "external owner source")?;
    let text = read_repo_utf8(root, &owner.path, "external owner source")?;
    let is_rust = path.extension().is_some_and(|extension| extension == "rs");
    let (producer, producer_block) = if owner.symbol == "<script>" && !is_rust {
        (text.as_bytes().to_vec(), None)
    } else if is_rust {
        (
            normalized_rust_function_closure(&text, [owner.symbol.clone()]).map_err(|detail| {
                format!(
                    "external owner {} target {}#{} has an invalid function closure: {detail}",
                    owner.id, owner.path, owner.symbol
                )
            })?,
            None,
        )
    } else {
        // CI shell functions commonly embed Python heredocs whose helpers are
        // part of the exact producer closure. Keep the entire outer function
        // block so nested-language helpers cannot move invisibly.
        let blocks = primary_script_blocks(&owner.path, &text, &owner.symbol)
            .map_err(|error| script_structure_detail(&owner.path, &text, &error))?;
        let [block] = blocks.as_slice() else {
            return Err(format!(
                "external owner {} target {}#{} resolves to {} primary script blocks",
                owner.id,
                owner.path,
                owner.symbol,
                blocks.len()
            ));
        };
        (
            block.text.as_bytes().to_vec(),
            Some((block.dialect, block.start)),
        )
    };
    let domain_is_bound = if is_rust {
        producer
            .windows(owner.domain.len())
            .any(|window| window == owner.domain.as_bytes())
    } else {
        let producer_text = std::str::from_utf8(&producer).map_err(|error| {
            format!(
                "external owner {} target {}#{} producer is not UTF-8: {error}",
                owner.id, owner.path, owner.symbol
            )
        })?;
        let views = match producer_block {
            Some((ScriptDialect::Python, start)) => {
                script_views(ScriptDialect::Python, producer_text)
                    .map_err(|error| error.shifted(start))
            }
            _ => executable_script_views(&owner.path, producer_text)
                .map_err(|error| error.shifted(producer_block.map_or(0, |(_, start)| start))),
        };
        let (code, source) =
            views.map_err(|error| script_structure_detail(&owner.path, &text, &error))?;
        script_declared_domains(&code, &source)
            .iter()
            .any(|domain| domain == &owner.domain)
    };
    if !domain_is_bound {
        return Err(format!(
            "external owner {} target {}#{} does not bind declared domain {:?} in its exact producer closure",
            owner.id, owner.path, owner.symbol, owner.domain
        ));
    }
    let mut payload = Vec::new();
    let version_bytes = owner.version.to_le_bytes();
    for part in [
        owner.id.as_bytes(),
        owner.path.as_bytes(),
        owner.symbol.as_bytes(),
        owner.domain.as_bytes(),
        version_bytes.as_slice(),
        producer.as_slice(),
    ] {
        fingerprint_part(&mut payload, part);
    }
    for exemption in exemptions
        .iter()
        .filter(|exemption| exemption.covered_by == owner.id && exemption.path == owner.path)
    {
        fingerprint_part(&mut payload, exemption.symbol.as_bytes());
        fingerprint_part(&mut payload, exemption.reason.as_bytes());
        fingerprint_part(
            &mut payload,
            &external_exemption_schema_bytes(root, exemption)?,
        );
    }
    Ok(format!(
        "v{}-{}",
        owner.version,
        schema_fingerprint_hex(&schema_fingerprint_digest(&payload))
    ))
}

fn external_owner_covers_exemption(
    root: &Path,
    owner: &ExternalOwner,
    exemption: &IdentityExemption,
) -> bool {
    if owner.path != exemption.path {
        return false;
    }
    external_exemption_schema_bytes(root, exemption).is_ok()
}

fn first_pattern<'a>(text: &str, patterns: &'a [&str]) -> Option<&'a str> {
    patterns
        .iter()
        .copied()
        .find(|pattern| text.contains(pattern))
}

#[allow(clippy::too_many_lines)] // Executable-token discovery must ignore Rust literals/comments exactly.
fn rust_code_view(text: &str) -> String {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    fn blank(byte: &mut u8) {
        if *byte != b'\n' && *byte != b'\r' {
            *byte = b' ';
        }
    }

    let source = text.as_bytes();
    let mut code = source.to_vec();
    let mut state = State::Normal;
    let mut index = 0usize;
    while index < source.len() {
        let byte = source[index];
        match state {
            State::Normal => {
                if byte == b'/' && source.get(index + 1) == Some(&b'/') {
                    blank(&mut code[index]);
                    blank(&mut code[index + 1]);
                    state = State::LineComment;
                    index += 2;
                } else if byte == b'/' && source.get(index + 1) == Some(&b'*') {
                    blank(&mut code[index]);
                    blank(&mut code[index + 1]);
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                } else if byte == b'r' {
                    let mut quote = index + 1;
                    while source.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if source.get(quote) == Some(&b'"') {
                        for byte in &mut code[index..=quote] {
                            blank(byte);
                        }
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                    } else {
                        index += 1;
                    }
                } else if byte == b'"' {
                    blank(&mut code[index]);
                    state = State::Quoted { escaped: false };
                    index += 1;
                } else if byte == b'\''
                    && let Some(end) = char_literal_end(source, index)
                {
                    for byte in &mut code[index..=end] {
                        blank(byte);
                    }
                    index = end + 1;
                } else {
                    index += 1;
                }
            }
            State::Quoted { escaped } => {
                blank(&mut code[index]);
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                blank(&mut code[index]);
                index += 1;
                if byte == b'"'
                    && source
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    for byte in &mut code[index..index + hashes] {
                        blank(byte);
                    }
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                blank(&mut code[index]);
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                blank(&mut code[index]);
                if byte == b'/' && source.get(index + 1) == Some(&b'*') {
                    blank(&mut code[index + 1]);
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && source.get(index + 1) == Some(&b'/') {
                    blank(&mut code[index + 1]);
                    state = if depth == 1 {
                        State::Normal
                    } else {
                        State::BlockComment { depth: depth - 1 }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    String::from_utf8(code).expect("blanking Rust UTF-8 bytes with ASCII spaces preserves UTF-8")
}

fn python_f_string_prefix(source: &[u8], quote: usize) -> Option<usize> {
    for width in [2_usize, 1] {
        let Some(start) = quote.checked_sub(width) else {
            continue;
        };
        let prefix = source.get(start..quote)?;
        let canonical = prefix
            .iter()
            .map(u8::to_ascii_lowercase)
            .collect::<Vec<_>>();
        if !matches!(canonical.as_slice(), b"f" | b"fr" | b"rf") {
            continue;
        }
        if start > 0
            && matches!(
                source[start - 1],
                b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'
            )
        {
            continue;
        }
        return Some(start);
    }
    None
}

fn python_quoted_end(source: &[u8], opening: usize) -> Option<usize> {
    let delimiter = *source.get(opening)?;
    let width = if source.get(opening..opening + 3) == Some(&[delimiter, delimiter, delimiter]) {
        3
    } else {
        1
    };
    let mut index = opening + width;
    while index < source.len() {
        if source[index] == b'\\' {
            index += if source.get(index + 1..index + 3) == Some(b"\r\n") {
                3
            } else {
                2
            };
            continue;
        }
        let closes = if width == 3 {
            source.get(index..index + 3) == Some(&[delimiter, delimiter, delimiter])
        } else {
            source.get(index) == Some(&delimiter)
        };
        if closes {
            return Some(index + width);
        }
        if width == 1 && source[index] == b'\n' {
            return None;
        }
        index += 1;
    }
    None
}

fn python_f_string_replacement_end(
    source: &[u8],
    start: usize,
    opening: usize,
) -> ScriptResult<usize> {
    let mut depth = 1_usize;
    let mut index = start;
    while index < source.len() {
        match source[index] {
            b'\'' | b'"' => {
                let Some(end) = python_quoted_end(source, index) else {
                    return Err(ScriptStructureError {
                        opening_offset: opening,
                        kind: ScriptStructureKind::IncompleteFStringReplacement,
                    });
                };
                index = end;
            }
            b'#' => {
                index += 1;
                while source.get(index).is_some_and(|byte| *byte != b'\n') {
                    index += 1;
                }
            }
            b'{' => {
                depth += 1;
                index += 1;
            }
            b'}' => {
                depth -= 1;
                if depth == 0 {
                    return Ok(index);
                }
                index += 1;
            }
            _ => index += 1,
        }
    }
    Err(ScriptStructureError {
        opening_offset: opening,
        kind: ScriptStructureKind::IncompleteFStringReplacement,
    })
}

fn python_f_string_span(text: &str, opening: usize) -> ScriptResult<(usize, Vec<(usize, usize)>)> {
    let source = text.as_bytes();
    let delimiter = source[opening];
    let width = if source.get(opening..opening + 3) == Some(&[delimiter, delimiter, delimiter]) {
        3
    } else {
        1
    };
    let mut replacements = Vec::new();
    let mut index = opening + width;
    while index < source.len() {
        if source[index] == b'\\' {
            index += if source.get(index + 1..index + 3) == Some(b"\r\n") {
                3
            } else {
                2
            };
            continue;
        }
        let closes = if width == 3 {
            source.get(index..index + 3) == Some(&[delimiter, delimiter, delimiter])
        } else {
            source.get(index) == Some(&delimiter)
        };
        if closes {
            return Ok((index + width, replacements));
        }
        match source[index] {
            b'{' if source.get(index + 1) == Some(&b'{') => index += 2,
            b'{' => {
                let expression_start = index + 1;
                let expression_end =
                    python_f_string_replacement_end(source, expression_start, index)?;
                replacements.push((expression_start, expression_end));
                index = expression_end + 1;
            }
            b'}' if source.get(index + 1) == Some(&b'}') => index += 2,
            b'}' => {
                return Err(ScriptStructureError {
                    opening_offset: index,
                    kind: ScriptStructureKind::IncompleteFStringReplacement,
                });
            }
            b'\n' if width == 1 => {
                return Err(ScriptStructureError {
                    opening_offset: opening,
                    kind: ScriptStructureKind::UnterminatedQuote {
                        dialect: ScriptDialect::Python,
                        delimiter: char::from(delimiter),
                    },
                });
            }
            _ => index += 1,
        }
    }
    let kind = if width == 3 {
        let delimiter = if delimiter == b'\'' { "'''" } else { "\"\"\"" };
        ScriptStructureKind::UnterminatedTripleString { delimiter }
    } else {
        ScriptStructureKind::UnterminatedQuote {
            dialect: ScriptDialect::Python,
            delimiter: char::from(delimiter),
        }
    };
    Err(ScriptStructureError {
        opening_offset: opening,
        kind,
    })
}

fn script_views(dialect: ScriptDialect, text: &str) -> ScriptResult<(String, String)> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted {
            delimiter: u8,
            escaped: bool,
            opening: usize,
        },
        Triple {
            delimiter: u8,
            escaped: bool,
            opening: usize,
        },
        Comment,
    }

    fn blank(byte: &mut u8) {
        if *byte != b'\n' && *byte != b'\r' {
            *byte = b' ';
        }
    }

    let source = text.as_bytes();
    let mut code = source.to_vec();
    let mut uncommented = source.to_vec();
    let mut state = State::Normal;
    let mut index = 0usize;
    while index < source.len() {
        let byte = source[index];
        match state {
            State::Normal if byte == b'#' => {
                blank(&mut code[index]);
                blank(&mut uncommented[index]);
                state = State::Comment;
                index += 1;
            }
            State::Normal
                if dialect == ScriptDialect::Python
                    && matches!(byte, b'\'' | b'"')
                    && python_f_string_prefix(source, index).is_some() =>
            {
                let prefix = python_f_string_prefix(source, index)
                    .expect("the guarded Python f-string prefix remains present");
                let (end, replacements) = python_f_string_span(text, index)?;
                for target in &mut code[prefix..end] {
                    blank(target);
                }
                for (start, expression_end) in replacements {
                    let (expression_code, expression_source) =
                        script_views(ScriptDialect::Python, &text[start..expression_end])
                            .map_err(|error| error.shifted(start))?;
                    code[start..expression_end].copy_from_slice(expression_code.as_bytes());
                    uncommented[start..expression_end]
                        .copy_from_slice(expression_source.as_bytes());
                }
                index = end;
            }
            State::Normal
                if dialect == ScriptDialect::Python
                    && matches!(byte, b'\'' | b'"')
                    && source.get(index..index + 3) == Some(&[byte, byte, byte]) =>
            {
                for target in &mut code[index..index + 3] {
                    blank(target);
                }
                state = State::Triple {
                    delimiter: byte,
                    escaped: false,
                    opening: index,
                };
                index += 3;
            }
            State::Normal if matches!(byte, b'\'' | b'"') => {
                blank(&mut code[index]);
                state = State::Quoted {
                    delimiter: byte,
                    escaped: false,
                    opening: index,
                };
                index += 1;
            }
            State::Normal => index += 1,
            State::Quoted {
                delimiter,
                escaped,
                opening,
            } => {
                blank(&mut code[index]);
                if dialect == ScriptDialect::Python && byte == b'\n' && !escaped {
                    return Err(ScriptStructureError {
                        opening_offset: opening,
                        kind: ScriptStructureKind::UnterminatedQuote {
                            dialect,
                            delimiter: char::from(delimiter),
                        },
                    });
                }
                let escaping_is_active = dialect == ScriptDialect::Python || delimiter == b'"';
                state = if escaped
                    && dialect == ScriptDialect::Python
                    && byte == b'\r'
                    && source.get(index + 1) == Some(&b'\n')
                {
                    State::Quoted {
                        delimiter,
                        escaped: true,
                        opening,
                    }
                } else if escaped {
                    State::Quoted {
                        delimiter,
                        escaped: false,
                        opening,
                    }
                } else if byte == b'\\' && escaping_is_active {
                    State::Quoted {
                        delimiter,
                        escaped: true,
                        opening,
                    }
                } else if byte == delimiter {
                    State::Normal
                } else {
                    State::Quoted {
                        delimiter,
                        escaped: false,
                        opening,
                    }
                };
                index += 1;
            }
            State::Triple {
                delimiter,
                escaped,
                opening,
            } => {
                blank(&mut code[index]);
                if escaped && byte == b'\r' && source.get(index + 1) == Some(&b'\n') {
                    state = State::Triple {
                        delimiter,
                        escaped: true,
                        opening,
                    };
                    index += 1;
                } else if escaped {
                    state = State::Triple {
                        delimiter,
                        escaped: false,
                        opening,
                    };
                    index += 1;
                } else if byte == b'\\' {
                    state = State::Triple {
                        delimiter,
                        escaped: true,
                        opening,
                    };
                    index += 1;
                } else if source.get(index..index + 3) == Some(&[delimiter, delimiter, delimiter]) {
                    for target in &mut code[index..index + 3] {
                        blank(target);
                    }
                    state = State::Normal;
                    index += 3;
                } else {
                    state = State::Triple {
                        delimiter,
                        escaped: false,
                        opening,
                    };
                    index += 1;
                }
            }
            State::Comment => {
                blank(&mut code[index]);
                blank(&mut uncommented[index]);
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
        }
    }
    match state {
        State::Quoted {
            delimiter, opening, ..
        } => {
            return Err(ScriptStructureError {
                opening_offset: opening,
                kind: ScriptStructureKind::UnterminatedQuote {
                    dialect,
                    delimiter: char::from(delimiter),
                },
            });
        }
        State::Triple {
            delimiter, opening, ..
        } => {
            let delimiter = if delimiter == b'\'' { "'''" } else { "\"\"\"" };
            return Err(ScriptStructureError {
                opening_offset: opening,
                kind: ScriptStructureKind::UnterminatedTripleString { delimiter },
            });
        }
        State::Normal | State::Comment => {}
    }
    Ok((
        String::from_utf8(code)
            .expect("blanking script UTF-8 bytes with ASCII spaces preserves UTF-8"),
        String::from_utf8(uncommented)
            .expect("blanking script comments with ASCII spaces preserves UTF-8"),
    ))
}

fn shell_line_continues(command: &str) -> bool {
    let code = match script_views(ScriptDialect::Shell, command) {
        Ok((code, _)) => code,
        Err(ScriptStructureError {
            kind:
                ScriptStructureKind::UnterminatedQuote {
                    dialect: ScriptDialect::Shell,
                    ..
                },
            ..
        }) => return true,
        Err(_) => return false,
    };
    let trimmed = code.trim_end_matches(char::is_whitespace);
    trimmed
        .as_bytes()
        .iter()
        .rev()
        .take_while(|byte| **byte == b'\\')
        .count()
        % 2
        == 1
}

fn shell_assignment_word(word: &str) -> bool {
    let Some((name, _)) = word.split_once('=') else {
        return false;
    };
    let mut bytes = name.bytes();
    bytes
        .next()
        .is_some_and(|byte| byte == b'_' || byte.is_ascii_alphabetic())
        && bytes.all(|byte| byte == b'_' || byte.is_ascii_alphanumeric())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HeredocBody {
    Opaque,
    Shell,
    Python,
}

fn interpreter_body(executable: &str) -> HeredocBody {
    match executable.rsplit('/').next().unwrap_or(executable) {
        "python" | "python3" => HeredocBody::Python,
        "bash" | "sh" | "zsh" => HeredocBody::Shell,
        _ => HeredocBody::Opaque,
    }
}

fn simple_command_body(command: &str) -> HeredocBody {
    let mut wrapper = None::<&str>;
    let mut skip_wrapper_argument = false;
    for raw_word in command.split_ascii_whitespace() {
        let word = raw_word.trim_matches('\\');
        if word.is_empty() || word == "!" || shell_assignment_word(word) {
            continue;
        }
        if skip_wrapper_argument {
            skip_wrapper_argument = false;
            continue;
        }
        if let Some(active_wrapper) = wrapper {
            if word == "--" {
                wrapper = None;
                continue;
            }
            if word.starts_with('-') {
                skip_wrapper_argument = matches!(
                    (active_wrapper, word),
                    ("env", "-u" | "--unset" | "-C" | "--chdir") | ("exec", "-a")
                );
                continue;
            }
        }
        let executable = word.rsplit('/').next().unwrap_or(word);
        if matches!(executable, "command" | "env" | "exec") {
            wrapper = Some(executable);
            continue;
        }
        return interpreter_body(executable);
    }
    HeredocBody::Opaque
}

fn heredoc_command_body(command: &str, heredoc: &ShellHeredoc) -> HeredocBody {
    let Ok((code, _)) = script_views(ScriptDialect::Shell, command) else {
        return HeredocBody::Opaque;
    };
    let Some(prefix) = code.get(..heredoc.operator_start) else {
        return HeredocBody::Opaque;
    };
    let bytes = prefix.as_bytes();
    let mut segment_start = 0usize;
    for (index, byte) in bytes.iter().copied().enumerate() {
        if matches!(byte, b';' | b'|' | b'&' | b'(' | b'{') {
            segment_start = index + 1;
        } else if byte == b'\n' {
            let backslashes = bytes[..index]
                .iter()
                .rev()
                .take_while(|candidate| **candidate == b'\\')
                .count();
            if backslashes % 2 == 0 {
                segment_start = index + 1;
            }
        }
    }
    let direct = simple_command_body(&prefix[segment_start..]);
    if direct != HeredocBody::Opaque {
        return direct;
    }

    let Some(suffix) = code.get(heredoc.operator_end..) else {
        return HeredocBody::Opaque;
    };
    let suffix_bytes = suffix.as_bytes();
    let mut index = 0usize;
    while index < suffix_bytes.len() {
        match suffix_bytes[index] {
            b'|' if suffix_bytes.get(index + 1) != Some(&b'|') => {
                let pipeline_stage = &suffix[index + 1..];
                let stage_end = pipeline_stage
                    .find(|character| matches!(character, ';' | '|' | '&' | '\n' | '{' | '}'))
                    .unwrap_or(pipeline_stage.len());
                return simple_command_body(&pipeline_stage[..stage_end]);
            }
            b';' | b'\n' | b'{' | b'}' => break,
            b'&' | b'|' => break,
            _ => index += 1,
        }
    }
    HeredocBody::Opaque
}

struct EmbeddedScriptViews {
    start: usize,
    dialect: ScriptDialect,
    code: String,
    source: String,
}

struct HeredocProjection {
    masked: String,
    embedded: Vec<EmbeddedScriptViews>,
}

struct PendingHeredoc {
    heredoc: ShellHeredoc,
    body: HeredocBody,
    opening_offset: usize,
    body_start: usize,
}

fn mask_non_executable_heredocs(text: &str) -> ScriptResult<HeredocProjection> {
    fn blank(bytes: &mut [u8]) {
        for byte in bytes {
            if *byte != b'\n' && *byte != b'\r' {
                *byte = b' ';
            }
        }
    }

    let mut bytes = text.as_bytes().to_vec();
    let mut pending = VecDeque::<PendingHeredoc>::new();
    let mut embedded = Vec::new();
    let mut offset = 0usize;
    let mut command_start = 0usize;
    while offset < text.len() {
        let line_end = text[offset..]
            .find('\n')
            .map_or(text.len(), |relative| offset + relative);
        let line = &text[offset..line_end];
        let next_line = (line_end + usize::from(line_end < text.len())).min(text.len());
        if let Some(pending_heredoc) = pending.front() {
            let candidate = if pending_heredoc.heredoc.strip_tabs {
                line.trim_start_matches('\t')
            } else {
                line
            };
            if candidate.trim_end_matches('\r') == pending_heredoc.heredoc.delimiter {
                let completed = pending
                    .pop_front()
                    .expect("the pending heredoc inspected above is present");
                let body_text = &text[completed.body_start..offset];
                let views = match completed.body {
                    HeredocBody::Opaque => {
                        if possibly_contains_identity_signal(body_text)
                            && script_sink_signal(body_text).is_some()
                        {
                            return Err(ScriptStructureError {
                                opening_offset: completed.opening_offset,
                                kind: ScriptStructureKind::AmbiguousHeredocExecution {
                                    delimiter: completed.heredoc.delimiter,
                                },
                            });
                        }
                        None
                    }
                    HeredocBody::Shell => {
                        let views = executable_script_views("<heredoc>.sh", body_text)
                            .map_err(|error| error.shifted(completed.body_start))?;
                        shell_function_blocks(body_text)
                            .map_err(|error| error.shifted(completed.body_start))?;
                        Some(views)
                    }
                    HeredocBody::Python => Some(
                        executable_script_views("<heredoc>.py", body_text)
                            .map_err(|error| error.shifted(completed.body_start))?,
                    ),
                };
                if let Some((code, source)) = views {
                    embedded.push(EmbeddedScriptViews {
                        start: completed.body_start,
                        dialect: match completed.body {
                            HeredocBody::Shell => ScriptDialect::Shell,
                            HeredocBody::Python => ScriptDialect::Python,
                            HeredocBody::Opaque => {
                                unreachable!("opaque heredocs never expose executable script views")
                            }
                        },
                        code,
                        source,
                    });
                }
                blank(&mut bytes[offset..next_line]);
                if let Some(next) = pending.front_mut() {
                    next.body_start = next_line;
                }
            } else {
                blank(&mut bytes[offset..next_line]);
            }
            command_start = next_line;
            offset = next_line;
            continue;
        }
        let logical_command = &text[command_start..line_end];
        if shell_line_continues(logical_command) {
            offset = next_line;
            continue;
        }
        pending.extend(shell_heredocs(logical_command).into_iter().map(|heredoc| {
            let body = heredoc_command_body(logical_command, &heredoc);
            PendingHeredoc {
                opening_offset: command_start + heredoc.operator_start,
                heredoc,
                body,
                body_start: next_line,
            }
        }));
        command_start = next_line;
        offset = next_line;
    }
    if let Some(pending) = pending.front() {
        return Err(ScriptStructureError {
            opening_offset: pending.opening_offset,
            kind: ScriptStructureKind::UnterminatedHeredoc {
                delimiter: pending.heredoc.delimiter.clone(),
            },
        });
    }
    Ok(HeredocProjection {
        masked: String::from_utf8(bytes)
            .expect("blanking heredoc bytes with ASCII spaces preserves UTF-8"),
        embedded,
    })
}

fn shell_executable_script_views(text: &str) -> ScriptResult<(String, String)> {
    let projection = mask_non_executable_heredocs(text)?;
    let (code, source) = script_views(ScriptDialect::Shell, &projection.masked)?;
    let mut code = code.into_bytes();
    let mut source = source.into_bytes();
    for embedded in projection.embedded {
        let end = embedded.start + embedded.code.len();
        debug_assert_eq!(embedded.code.len(), embedded.source.len());
        debug_assert!(end <= code.len());
        code[embedded.start..end].copy_from_slice(embedded.code.as_bytes());
        source[embedded.start..end].copy_from_slice(embedded.source.as_bytes());
    }
    Ok((
        String::from_utf8(code).expect("embedded script code preserves UTF-8"),
        String::from_utf8(source).expect("embedded script source preserves UTF-8"),
    ))
}

fn executable_script_views(path: &str, text: &str) -> ScriptResult<(String, String)> {
    match script_dialect(path, text).unwrap_or(ScriptDialect::Python) {
        ScriptDialect::Shell => shell_executable_script_views(text),
        ScriptDialect::Python => script_views(ScriptDialect::Python, text),
    }
}

fn quoted_literal_before(source: &str, at: usize) -> Option<&str> {
    let bytes = source.as_bytes();
    let mut end = at;
    while end > 0 && bytes[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    let delimiter = *bytes.get(end.checked_sub(1)?)?;
    if !matches!(delimiter, b'\'' | b'"') {
        return None;
    }
    let mut start = end - 1;
    while start > 0 {
        start -= 1;
        if bytes[start] == delimiter {
            let escapes = bytes[..start]
                .iter()
                .rev()
                .take_while(|byte| **byte == b'\\')
                .count();
            if escapes % 2 == 0 {
                return source.get(start + 1..end - 1);
            }
        }
    }
    None
}

fn quoted_literal_after(source: &str, at: usize) -> Option<&str> {
    let bytes = source.as_bytes();
    let mut start = at;
    while bytes.get(start).is_some_and(u8::is_ascii_whitespace) {
        start += 1;
    }
    let delimiter = *bytes.get(start)?;
    if !matches!(delimiter, b'\'' | b'"') {
        return None;
    }
    let mut end = start + 1;
    let mut escaped = false;
    while let Some(byte) = bytes.get(end).copied() {
        if escaped {
            escaped = false;
        } else if byte == b'\\' {
            escaped = true;
        } else if byte == delimiter {
            return source.get(start + 1..end);
        }
        end += 1;
    }
    None
}

fn structured_json_identity_domains(code: &str, source: &str) -> Vec<String> {
    let bytes = code.as_bytes();
    let mut domains = Vec::new();
    let mut search = 0usize;
    while let Some(relative) = code[search..].find("json.dumps") {
        let call = search + relative;
        let Some(open_relative) = code[call..].find('{') else {
            break;
        };
        let open = call + open_relative;
        let mut depth = 0usize;
        let mut domain = None;
        let mut has_version = false;
        let mut index = open;
        while index < bytes.len() {
            match bytes[index] {
                b'{' => depth += 1,
                b'}' => {
                    let Some(next_depth) = depth.checked_sub(1) else {
                        break;
                    };
                    depth = next_depth;
                    if depth == 0 {
                        if has_version && let Some(domain) = domain {
                            domains.push(domain);
                        }
                        index += 1;
                        break;
                    }
                }
                b':' if depth == 1 => {
                    let key = quoted_literal_before(source, index);
                    if key == Some("identity_domain") {
                        domain = quoted_literal_after(source, index + 1).map(str::to_string);
                    } else if key == Some("identity_version") {
                        let mut value = index + 1;
                        while bytes.get(value).is_some_and(u8::is_ascii_whitespace) {
                            value += 1;
                        }
                        has_version = bytes.get(value).is_some_and(u8::is_ascii_digit);
                    }
                }
                _ => {}
            }
            index += 1;
        }
        search = index.max(call + "json.dumps".len());
    }
    domains.sort();
    domains.dedup();
    domains
}

fn script_declared_domains(code: &str, source: &str) -> Vec<String> {
    fn token_at(bytes: &[u8], start: usize, token: &[u8]) -> bool {
        bytes.get(start..start + token.len()) == Some(token)
            && (start == 0
                || !matches!(
                    bytes[start - 1],
                    b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'
                ))
            && bytes
                .get(start + token.len())
                .is_none_or(|byte| !matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'))
    }

    fn bare_value_after(source: &str, at: usize) -> Option<&str> {
        let bytes = source.as_bytes();
        let mut start = at;
        while bytes.get(start).is_some_and(u8::is_ascii_whitespace) {
            start += 1;
        }
        let mut end = start;
        while bytes.get(end).is_some_and(|byte| {
            !byte.is_ascii_whitespace() && !matches!(byte, b',' | b';' | b')' | b'}')
        }) {
            end += 1;
        }
        (end > start).then(|| &source[start..end])
    }

    let bytes = code.as_bytes();
    let mut domains = structured_json_identity_domains(code, source);
    for token in [b"IDENTITY_DOMAIN".as_slice(), b"DOMAIN".as_slice()] {
        let mut start = 0usize;
        while start + token.len() <= bytes.len() {
            if !token_at(bytes, start, token) {
                start += 1;
                continue;
            }
            let mut equals = start + token.len();
            while bytes.get(equals).is_some_and(u8::is_ascii_whitespace) {
                equals += 1;
            }
            if bytes.get(equals) == Some(&b'=') && bytes.get(equals + 1) != Some(&b'=') {
                let value = quoted_literal_after(source, equals + 1)
                    .or_else(|| bare_value_after(source, equals + 1));
                if let Some(value) = value.filter(|value| !value.is_empty()) {
                    domains.push(value.to_string());
                }
            }
            start += token.len();
        }
    }

    let mut search = 0usize;
    while let Some(relative) = code[search..].find("add") {
        let add = search + relative;
        if !token_at(bytes, add, b"add") {
            search = add + 1;
            continue;
        }
        let mut open = add + 3;
        while bytes.get(open).is_some_and(u8::is_ascii_whitespace) {
            open += 1;
        }
        if bytes.get(open) != Some(&b'(') {
            search = add + 3;
            continue;
        }
        let Some(comma_relative) = code[open + 1..].find(',') else {
            break;
        };
        let comma = open + 1 + comma_relative;
        if quoted_literal_after(source, open + 1) == Some("identity-domain")
            && let Some(domain) = quoted_literal_after(source, comma + 1)
        {
            domains.push(domain.to_string());
        }
        search = comma + 1;
    }
    domains.sort();
    domains.dedup();
    domains
}

fn script_identity_signal(code: &str, source: &str) -> Option<String> {
    let domains = script_declared_domains(code, source);
    if domains.is_empty() {
        return None;
    }
    candidate_identity_signal(code).or_else(|| {
        (!structured_json_identity_domains(code, source).is_empty())
            .then(|| "structured-json-identity".to_string())
    })
}

fn candidate_identity_signal(block: &str) -> Option<String> {
    let block = block.to_ascii_lowercase();
    let digest = first_pattern(
        &block,
        &[
            "hash_domain",
            "hashlib.",
            "sha256",
            "sha-256",
            "shasum -a",
            "blake3",
            "fnv1a",
            ".digest(",
            ".hexdigest(",
        ],
    );
    let vocabulary = first_pattern(
        &block,
        &[
            "identity_domain",
            "identity-domain",
            "identity_version",
            "identity-version",
            "receipt_identity",
            "snapshot_identity",
            "schema_fingerprint",
            "canonical_bytes",
            "to_canonical_",
            "cache_key",
            "cache-key",
            "content_address",
            "content-address",
            "idempotency",
            "merkle_root",
            "claim_subject",
            "semantic_key",
            "content_hash",
        ],
    );
    if digest.is_none()
        && vocabulary.is_some()
        && first_pattern(
            &block,
            &[
                "additive_identity",
                "multiplicative_identity",
                "matrix_identity",
                "identity_element",
            ],
        )
        .is_some()
    {
        return None;
    }
    if digest.is_none()
        && first_pattern(
            &block,
            &["default_hasher", "process_local", "process-local hash"],
        )
        .is_some()
    {
        return None;
    }
    if digest.is_none()
        && block.contains("diagnostic")
        && first_pattern(&block, &["process::id", "pid", "log_path", "log-path"]).is_some()
    {
        return None;
    }
    if digest.is_none()
        && block.lines().next().is_some_and(|symbol| {
            symbol.starts_with("classify_") && symbol.ends_with("_identity_fields")
        })
    {
        return None;
    }
    Some(digest.or(vocabulary)?.to_string())
}

fn possibly_contains_identity_signal(text: &str) -> bool {
    let text = text.to_ascii_lowercase();
    first_pattern(
        &text,
        &[
            "hash_domain",
            "hashlib.",
            "sha256",
            "sha-256",
            "shasum -a",
            "blake3",
            "fnv1a",
            ".digest(",
            ".hexdigest(",
            "identity_domain",
            "identity-domain",
            "identity_version",
            "identity-version",
            "receipt_identity",
            "snapshot_identity",
            "schema_fingerprint",
            "canonical_bytes",
            "to_canonical_",
            "cache_key",
            "cache-key",
            "content_address",
            "content-address",
            "idempotency",
            "merkle_root",
            "claim_subject",
            "semantic_key",
            "content_hash",
        ],
    )
    .is_some()
}

fn rust_sink_signal(block: &str) -> Option<String> {
    let block = block.to_ascii_lowercase();
    first_pattern(
        &block,
        &[
            "std::fs::write",
            "fs::write",
            ".write_all(",
            ".append(",
            ".put(",
            ".persist(",
            ".save(",
            ".commit(",
            "insert into",
            "verificationdecision::reject",
        ],
    )
    .map(str::to_string)
    .or_else(|| {
        (block.contains(".insert(")
            && first_pattern(
                &block,
                &[
                    "cache",
                    "store",
                    "ledger",
                    "idempotency",
                    "session",
                    "receipt",
                ],
            )
            .is_some())
        .then(|| ".insert(identity-store-context)".to_string())
    })
}

fn script_sink_signal(block: &str) -> Option<String> {
    let block = block.to_ascii_lowercase();
    first_pattern(
        &block,
        &[
            ".write_text(",
            ".write_bytes(",
            "json.dump(",
            "tee -a",
            ">>",
            "print(json.dumps({",
        ],
    )
    .map(str::to_string)
    .or_else(|| {
        (block.contains("snapshot") && block.contains("hexdigest") && block.contains("print("))
            .then(|| "snapshot-hexdigest-output".to_string())
    })
}

#[allow(dead_code)] // Retained for a future resolved-call discovery tier.
fn rust_called_symbols(block: &str) -> BTreeSet<String> {
    #[derive(Clone, Copy)]
    enum State {
        Normal,
        Quoted { escaped: bool },
        Raw { hashes: usize },
        LineComment,
        BlockComment { depth: usize },
    }

    let bytes = block.as_bytes();
    let mut state = State::Normal;
    let mut index = 0_usize;
    let mut symbols = BTreeSet::new();
    while index < bytes.len() {
        let byte = bytes[index];
        match state {
            State::Normal => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'/') {
                    state = State::LineComment;
                    index += 2;
                    continue;
                }
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: 1 };
                    index += 2;
                    continue;
                }
                if byte == b'r' {
                    let mut quote = index + 1;
                    while bytes.get(quote) == Some(&b'#') {
                        quote += 1;
                    }
                    if bytes.get(quote) == Some(&b'"') {
                        state = State::Raw {
                            hashes: quote - index - 1,
                        };
                        index = quote + 1;
                        continue;
                    }
                }
                if byte == b'"' {
                    state = State::Quoted { escaped: false };
                    index += 1;
                    continue;
                }
                if byte == b'\''
                    && let Some(end) = char_literal_end(bytes, index)
                {
                    index = end + 1;
                    continue;
                }
                if matches!(byte, b'_' | b'a'..=b'z' | b'A'..=b'Z') {
                    let start = index;
                    index += 1;
                    while bytes.get(index).is_some_and(
                        |byte| matches!(byte, b'_' | b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z'),
                    ) {
                        index += 1;
                    }
                    let symbol = block.get(start..index).unwrap_or_default();
                    while bytes.get(index).is_some_and(u8::is_ascii_whitespace) {
                        index += 1;
                    }
                    if bytes.get(index) == Some(&b'!') {
                        index += 1;
                        while bytes.get(index).is_some_and(u8::is_ascii_whitespace) {
                            index += 1;
                        }
                    }
                    if bytes.get(index) == Some(&b'(') && canonical_symbol(symbol) {
                        symbols.insert(symbol.to_string());
                    }
                    continue;
                }
                index += 1;
            }
            State::Quoted { escaped } => {
                state = if escaped {
                    State::Quoted { escaped: false }
                } else if byte == b'\\' {
                    State::Quoted { escaped: true }
                } else if byte == b'"' {
                    State::Normal
                } else {
                    State::Quoted { escaped: false }
                };
                index += 1;
            }
            State::Raw { hashes } => {
                index += 1;
                if byte == b'"'
                    && bytes
                        .get(index..index.saturating_add(hashes))
                        .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                {
                    index += hashes;
                    state = State::Normal;
                }
            }
            State::LineComment => {
                index += 1;
                if byte == b'\n' {
                    state = State::Normal;
                }
            }
            State::BlockComment { depth } => {
                if byte == b'/' && bytes.get(index + 1) == Some(&b'*') {
                    state = State::BlockComment { depth: depth + 1 };
                    index += 2;
                } else if byte == b'*' && bytes.get(index + 1) == Some(&b'/') {
                    let depth = depth - 1;
                    state = if depth == 0 {
                        State::Normal
                    } else {
                        State::BlockComment { depth }
                    };
                    index += 2;
                } else {
                    index += 1;
                }
            }
        }
    }
    symbols
}

fn rust_without_test_modules(text: &str) -> String {
    let mut module_names = BTreeSet::from(["tests".to_string()]);
    let mut cfg_test = false;
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed == "#[cfg(test)]" {
            cfg_test = true;
            continue;
        }
        if cfg_test
            && (trimmed.is_empty()
                || trimmed.starts_with("//")
                || trimmed.starts_with("/*")
                || trimmed.starts_with('*')
                || trimmed.ends_with("*/")
                || trimmed.starts_with("#["))
        {
            continue;
        }
        if cfg_test
            && let Some(rest) = trimmed
                .strip_prefix("mod ")
                .or_else(|| trimmed.strip_prefix("pub mod "))
            && let Some(name) = rest
                .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
                .find(|token| !token.is_empty())
            && canonical_symbol(name)
        {
            module_names.insert(name.to_string());
        }
        cfg_test = false;
    }

    // Test-only helpers and their writes must contribute neither candidates nor
    // module-wide sink evidence. Preserve byte offsets/newlines for the lexer.
    let mut bytes = text.as_bytes().to_vec();
    for module in module_names {
        for body in braced_bodies(text, &format!("mod {module}"), false) {
            let start = body.as_ptr() as usize - text.as_ptr() as usize;
            for byte in &mut bytes[start..start + body.len()] {
                if *byte != b'\n' {
                    *byte = b' ';
                }
            }
        }
    }
    String::from_utf8(bytes).expect("blanking Rust source bytes preserves UTF-8")
}

fn discover_rust_candidates(path: &str, text: &str) -> Vec<IdentityCandidate> {
    #[derive(Clone)]
    struct FunctionItem<'a> {
        reference: String,
        symbol: String,
        start: usize,
        body_start: usize,
        body_end: usize,
        body: &'a str,
    }

    #[derive(Clone)]
    struct OwnerScope {
        start: usize,
        end: usize,
        name: String,
    }

    fn owner_scopes(
        text: &str,
        declaration: &str,
        owner: impl Fn(&str) -> Option<String>,
    ) -> Vec<OwnerScope> {
        let mut scopes = Vec::new();
        for fragment in braced_declaration_fragments(text, declaration) {
            let Some(body) = braced_bodies(fragment, declaration, false)
                .into_iter()
                .next()
            else {
                continue;
            };
            let body_start_in_fragment = body.as_ptr() as usize - fragment.as_ptr() as usize;
            let Some(open) = body_start_in_fragment.checked_sub(1) else {
                continue;
            };
            let Some(name) = fragment.get(..open).and_then(&owner) else {
                continue;
            };
            let start = body.as_ptr() as usize - text.as_ptr() as usize;
            scopes.push(OwnerScope {
                start,
                end: start + body.len(),
                name,
            });
        }
        scopes.sort_by(|left, right| (left.start, left.end).cmp(&(right.start, right.end)));
        scopes.dedup_by(|left, right| {
            left.start == right.start && left.end == right.end && left.name == right.name
        });
        scopes
    }

    fn module_owner(header: &str) -> Option<String> {
        let tokens = header
            .split(|character: char| character != '_' && !character.is_ascii_alphanumeric())
            .filter(|token| !token.is_empty())
            .collect::<Vec<_>>();
        let module = tokens.iter().position(|token| *token == "mod")?;
        let name = *tokens.get(module + 1)?;
        canonical_symbol(name).then(|| name.to_string())
    }

    let production = rust_without_test_modules(text);
    let text = production.as_str();
    let starts = rust_function_starts(text);
    let modules = owner_scopes(text, "mod", module_owner);
    let implementations = owner_scopes(text, "impl", implementation_owner);
    let mut items = Vec::<FunctionItem<'_>>::new();
    for (index, (symbol, start)) in starts.iter().enumerate() {
        let suffix = &text[*start..];
        let declaration = format!("fn {symbol}");
        let Some(body) = first_direct_braced_body(suffix, &declaration) else {
            continue;
        };
        let body_start = body.as_ptr() as usize - text.as_ptr() as usize;
        if starts
            .get(index + 1)
            .is_some_and(|(_, next_start)| *next_start < body_start)
        {
            // The current item ended in a semicolon. A later same-name body
            // must not be attributed to this trait declaration.
            continue;
        }
        let body_end = body_start + body.len();
        items.push(FunctionItem {
            reference: String::new(),
            symbol: symbol.clone(),
            start: *start,
            body_start,
            body_end,
            body,
        });
    }

    for index in 0..items.len() {
        let nested_in_function = items.iter().enumerate().any(|(other_index, other)| {
            other_index != index
                && other.body_start <= items[index].start
                && items[index].start < other.body_end
        });
        let mut segments = modules
            .iter()
            .filter(|scope| scope.start <= items[index].start && items[index].start < scope.end)
            .collect::<Vec<_>>();
        segments.sort_by(|left, right| {
            (left.start, std::cmp::Reverse(left.end))
                .cmp(&(right.start, std::cmp::Reverse(right.end)))
        });
        let mut reference = segments
            .into_iter()
            .map(|scope| scope.name.as_str())
            .collect::<Vec<_>>();
        if !nested_in_function
            && let Some(implementation) = implementations
                .iter()
                .filter(|scope| scope.start <= items[index].start && items[index].start < scope.end)
                .min_by_key(|scope| scope.end - scope.start)
        {
            reference.push(implementation.name.as_str());
        }
        reference.push(items[index].symbol.as_str());
        items[index].reference = reference.join("::");
    }
    items.sort_by(|left, right| {
        (left.body.as_ptr() as usize, &left.reference)
            .cmp(&(right.body.as_ptr() as usize, &right.reference))
    });
    items.dedup_by(|left, right| {
        left.body.as_ptr() == right.body.as_ptr() && left.reference == right.reference
    });

    let mut sinks = BTreeMap::<usize, (String, String)>::new();
    for (index, item) in items.iter().enumerate() {
        let code = rust_code_view(item.body);
        if let Some(sink) = rust_sink_signal(&code) {
            sinks.insert(index, (sink, code));
        }
    }

    let mut candidates = BTreeMap::<String, IdentityCandidate>::new();
    for (index, (sink_signal, code)) in sinks {
        let item = &items[index];
        let Some(identity_signal) = candidate_identity_signal(&code) else {
            continue;
        };
        if has_test_function(text, &item.symbol) {
            continue;
        }
        candidates
            .entry(item.reference.clone())
            .or_insert_with(|| IdentityCandidate {
                path: path.to_string(),
                symbol: item.reference.clone(),
                identity_signal,
                sink_signal,
            });
    }
    candidates.into_values().collect()
}

fn discover_script_candidates(path: &str, text: &str) -> ScriptResult<Vec<IdentityCandidate>> {
    let mut blocks = Vec::new();
    let dialect = script_dialect(path, text).unwrap_or(ScriptDialect::Python);
    let shell = if dialect == ScriptDialect::Shell {
        shell_function_blocks(text)?
    } else {
        Vec::new()
    };
    let python = python_function_blocks(path, text)?;
    if dialect == ScriptDialect::Python {
        blocks.extend(python.iter().cloned());
    } else {
        blocks.extend(shell.iter().map(|(symbol, block)| ScriptBlock {
            symbol: symbol.clone(),
            text: block,
            dialect: ScriptDialect::Shell,
            start: block.as_ptr() as usize - text.as_ptr() as usize,
        }));
    }
    if dialect == ScriptDialect::Shell {
        blocks.extend(python.iter().cloned());
    }
    let excluded = blocks.iter().map(|block| block.text).collect::<Vec<_>>();
    let mut candidates = Vec::new();
    for block in &blocks {
        let block_end = block.start.saturating_add(block.text.len());
        let nested = blocks
            .iter()
            .filter(|other| {
                other.start > block.start
                    && other.start.saturating_add(other.text.len()) <= block_end
            })
            .map(|other| other.text)
            .collect::<Vec<_>>();
        let scoped_text = text_outside_blocks(block.text, &nested);
        let scoped_block = ScriptBlock {
            symbol: block.symbol.clone(),
            text: &scoped_text,
            dialect: block.dialect,
            start: block.start,
        };
        let (code, source) = script_block_views(path, &scoped_block)?;
        if let Some(identity_signal) = script_identity_signal(&code, &source)
            && let Some(sink_signal) = script_sink_signal(&code)
        {
            candidates.push(IdentityCandidate {
                path: path.to_string(),
                symbol: block.symbol.clone(),
                identity_signal,
                sink_signal,
            });
        }
    }
    let top_level = text_outside_blocks(text, &excluded);
    let (top_level_code, top_level_source) = executable_script_views(path, &top_level)?;
    if let Some(identity_signal) = script_identity_signal(&top_level_code, &top_level_source)
        && let Some(sink_signal) = script_sink_signal(&top_level_code)
    {
        candidates.push(IdentityCandidate {
            path: path.to_string(),
            symbol: "<script>".to_string(),
            identity_signal,
            sink_signal,
        });
    }
    Ok(candidates)
}

fn text_outside_blocks(text: &str, blocks: &[&str]) -> String {
    let mut bytes = text.as_bytes().to_vec();
    for block in blocks {
        let start = block.as_ptr() as usize - text.as_ptr() as usize;
        let end = start.saturating_add(block.len()).min(bytes.len());
        for byte in &mut bytes[start..end] {
            if *byte != b'\n' && *byte != b'\r' {
                *byte = b' ';
            }
        }
    }
    String::from_utf8(bytes).expect("replacing UTF-8 bytes with ASCII spaces preserves UTF-8")
}

fn discover_identity_candidates(
    root: &Path,
    manifest: &AuthorityManifest,
) -> Result<Vec<IdentityCandidate>, Vec<Violation>> {
    fn visit_scripts(
        root: &Path,
        path: &Path,
        scripts: &mut Vec<PathBuf>,
        violations: &mut Vec<Violation>,
    ) {
        if path != root && path.join(".git").exists() {
            return;
        }
        let context = path
            .strip_prefix(root)
            .ok()
            .filter(|relative| !relative.as_os_str().is_empty())
            .map_or_else(
                || "<repo>".to_string(),
                |relative| relative.display().to_string(),
            );
        let read_dir = match std::fs::read_dir(path) {
            Ok(entries) => entries,
            Err(error) => {
                violations.push(identity_violation(
                    &context,
                    format!("identity script discovery cannot read directory: {error}"),
                ));
                return;
            }
        };
        let mut entries = Vec::new();
        for entry in read_dir {
            match entry {
                Ok(entry) => entries.push(entry),
                Err(error) => violations.push(identity_violation(
                    &context,
                    format!("identity script discovery cannot read directory entry: {error}"),
                )),
            }
        }
        entries.sort_by_key(std::fs::DirEntry::file_name);
        for entry in entries {
            let path = entry.path();
            let relative = path
                .strip_prefix(root)
                .unwrap_or(&path)
                .display()
                .to_string();
            let file_type = match entry.file_type() {
                Ok(file_type) => file_type,
                Err(error) => {
                    violations.push(identity_violation(
                        &relative,
                        format!("identity script discovery cannot inspect file type: {error}"),
                    ));
                    continue;
                }
            };
            if file_type.is_symlink() {
                violations.push(identity_violation(
                    &relative,
                    "identity script discovery refuses symlinked paths",
                ));
            } else if file_type.is_dir() {
                if matches!(
                    entry.file_name().to_str(),
                    Some(".git" | ".beads" | ".bv" | ".doctor" | ".wrangler" | "target")
                ) {
                    continue;
                }
                visit_scripts(root, &path, scripts, violations);
            } else if file_type.is_file() {
                let extension = path.extension().and_then(|value| value.to_str());
                if matches!(extension, Some("sh" | "py")) || extension.is_none() {
                    scripts.push(path);
                }
            } else {
                violations.push(identity_violation(
                    &relative,
                    "identity script discovery refuses non-file paths",
                ));
            }
        }
    }

    let mut candidates = Vec::new();
    let (sources, mut violations) = source_files(root);
    for path in sources {
        let relative = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .display()
            .to_string();
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => {
                violations.push(identity_violation(
                    &relative,
                    format!("identity source discovery cannot read Rust source: {error}"),
                ));
                continue;
            }
        };
        if relative == GATE_IMPLEMENTATION_PATH
            || relative.split('/').any(|segment| segment == "tests")
        {
            continue;
        }
        if !possibly_contains_identity_signal(&text) || rust_sink_signal(&text).is_none() {
            continue;
        }
        candidates.extend(discover_rust_candidates(&relative, &text));
    }
    let mut scripts = Vec::new();
    visit_scripts(root, root, &mut scripts, &mut violations);
    scripts.extend(
        manifest
            .external_owners
            .iter()
            .map(|owner| owner.path.as_str())
            .chain(
                manifest
                    .exemptions
                    .iter()
                    .map(|exemption| exemption.path.as_str()),
            )
            .filter(|path| !path.ends_with(".rs"))
            .map(|path| root.join(path)),
    );
    scripts.sort();
    scripts.dedup();
    for path in scripts {
        let relative = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .display()
            .to_string();
        let registered_target = manifest
            .external_owners
            .iter()
            .any(|owner| owner.path == relative)
            || manifest
                .exemptions
                .iter()
                .any(|exemption| exemption.path == relative);
        let known_extension = path
            .extension()
            .and_then(|value| value.to_str())
            .is_some_and(|extension| matches!(extension, "sh" | "py"));
        let metadata = match std::fs::metadata(&path) {
            Ok(metadata) => metadata,
            Err(error) => {
                if registered_target || known_extension {
                    violations.push(identity_violation(
                        &relative,
                        format!("identity script source metadata is unavailable: {error}"),
                    ));
                }
                continue;
            }
        };
        if metadata.len() > MAX_SOURCE_BYTES {
            violations.push(identity_violation(
                &relative,
                format!("candidate script exceeds identity scan cap of {MAX_SOURCE_BYTES} bytes"),
            ));
            continue;
        }
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) => {
                violations.push(identity_violation(
                    &relative,
                    format!("candidate script is unreadable: {error}"),
                ));
                continue;
            }
        };
        let text = match String::from_utf8(bytes) {
            Ok(text) => text,
            Err(error) => {
                let bytes = error.as_bytes();
                let prefix = String::from_utf8_lossy(&bytes[..bytes.len().min(256)]);
                if registered_target
                    || known_extension
                    || script_dialect(&relative, &prefix).is_some()
                {
                    violations.push(identity_violation(
                        &relative,
                        format!("candidate script is not valid UTF-8: {error}"),
                    ));
                }
                continue;
            }
        };
        let Some(_) = script_dialect(&relative, &text) else {
            if registered_target {
                violations.push(identity_violation(
                    &relative,
                    "registered script target needs a .sh/.py suffix or supported shell/Python shebang",
                ));
            }
            continue;
        };
        let raw_producer_hint =
            possibly_contains_identity_signal(&text) && script_sink_signal(&text).is_some();
        if !registered_target && !raw_producer_hint {
            continue;
        }
        match discover_script_candidates(&relative, &text) {
            Ok(discovered) => candidates.extend(discovered),
            Err(error) => violations.push(identity_violation(
                &relative,
                script_structure_detail(&relative, &text, &error),
            )),
        }
    }
    candidates.sort();
    candidates.dedup_by(|left, right| left.path == right.path && left.symbol == right.symbol);
    if violations.is_empty() {
        Ok(candidates)
    } else {
        Err(violations)
    }
}

fn authority_target_matches(candidate: &IdentityCandidate, path: &str, symbol: &str) -> bool {
    candidate.path == path && candidate.symbol == symbol
}

fn authority_target_error(root: &Path, path: &str, symbol: &str) -> Option<String> {
    let target = match checked_repo_file(root, path, "authority target") {
        Ok(target) => target,
        Err(error) => return Some(error),
    };
    let text = match std::fs::read_to_string(&target) {
        Ok(text) => text,
        Err(error) => return Some(format!("target is unreadable UTF-8: {error}")),
    };
    let is_rust = target
        .extension()
        .is_some_and(|extension| extension == "rs");
    if symbol == "<script>" {
        return is_rust.then(|| "<script> is valid only for script sources".to_string());
    }
    let present = if is_rust {
        has_function(&text, symbol)
    } else {
        let blocks = match script_symbol_blocks(path, &text, symbol) {
            Ok(blocks) => blocks,
            Err(error) => return Some(script_structure_detail(path, &text, &error)),
        };
        blocks.len() == 1
    };
    (!present).then(|| "target symbol does not exist exactly once".to_string())
}

fn top_level_external_owner_error(root: &Path, path: &str) -> Option<String> {
    let target = match checked_repo_file(root, path, "top-level external authority target") {
        Ok(target) => target,
        Err(error) => return Some(error),
    };
    let text = match std::fs::read_to_string(&target) {
        Ok(text) => text,
        Err(error) => return Some(format!("target is unreadable UTF-8: {error}")),
    };
    match script_has_function_signature(path, &text) {
        Ok(true) => Some(
            "<script> cannot own a source containing function scopes; register the exact producer function and cover any true top-level child explicitly"
                .to_string(),
        ),
        Ok(false) => None,
        Err(error) => Some(script_structure_detail(path, &text, &error)),
    }
}

fn declaration_authority_targets(declaration: &IdentityDecl) -> Vec<(String, String)> {
    let mut targets = Vec::new();
    for symbol in std::iter::once(&declaration.encoder)
        .chain(declaration.encoder_helpers.iter())
        .chain(std::iter::once(&declaration.transport_guard))
    {
        targets.push((declaration.owner.clone(), symbol.clone()));
    }
    targets.sort();
    targets.dedup();
    targets
}

fn declaration_schema_targets(declaration: &IdentityDecl) -> Vec<(String, String)> {
    let mut targets = declaration_authority_targets(declaration);
    for function in &declaration.schema_functions {
        targets.push((
            function
                .path
                .clone()
                .unwrap_or_else(|| declaration.owner.clone()),
            function.symbol.clone(),
        ));
    }
    targets.sort();
    targets.dedup();
    targets
}

#[allow(clippy::too_many_lines)] // One pass keeps inventory, target, and coverage diagnostics coherent.
fn authority_violations_against(
    root: &Path,
    declarations: &[IdentityDecl],
    manifest: &AuthorityManifest,
    candidates: &[IdentityCandidate],
) -> Vec<Violation> {
    const EXEMPTION_REASONS: &[&str] = &[
        "child-digest-helper",
        "ephemeral-staging-key",
        "identity-consumer-not-producer",
        "transient-scheduling-child",
        "diagnostic-only",
        "generated-registry-output",
        "negative-refusal-sentinel",
        "mathematical-identity",
    ];

    let mut violations = Vec::new();
    let declaration_ids = declarations
        .iter()
        .map(|declaration| declaration.id.clone())
        .collect::<BTreeSet<_>>();
    let external_ids = manifest
        .external_owners
        .iter()
        .map(|owner| owner.id.clone())
        .collect::<BTreeSet<_>>();
    let registered_ids = declaration_ids
        .union(&external_ids)
        .cloned()
        .collect::<BTreeSet<_>>();

    for id in manifest.required_ids.difference(&registered_ids) {
        violations.push(identity_violation(
            "<repo>",
            format!(
                "{AUTHORITY_FILE} requires identity authority {id:?}, but no owner declaration or external owner row registers it"
            ),
        ));
    }
    for id in registered_ids.difference(&manifest.required_ids) {
        violations.push(identity_violation(
            "<repo>",
            format!(
                "registered identity authority {id:?} is absent from {AUTHORITY_FILE} required_ids"
            ),
        ));
    }
    for id in declaration_ids.intersection(&external_ids) {
        violations.push(identity_violation(
            "<repo>",
            format!(
                "identity authority {id:?} is registered by both an owner declaration and an external owner row"
            ),
        ));
    }

    for owner in &manifest.external_owners {
        if let Some(detail) = authority_target_error(root, &owner.path, &owner.symbol) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "external authority {} target {}#{} is invalid: {detail}",
                    owner.id, owner.path, owner.symbol
                ),
            ));
        }
        if owner.symbol == "<script>"
            && let Some(detail) = top_level_external_owner_error(root, &owner.path)
        {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "external authority {} target {}#{} is invalid: {detail}",
                    owner.id, owner.path, owner.symbol
                ),
            ));
        }
        let discovered = candidates
            .iter()
            .filter(|candidate| authority_target_matches(candidate, &owner.path, &owner.symbol))
            .count();
        if discovered != 1 {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "external authority {} target {}#{} matches {discovered} independently discovered producers; exactly one is required",
                    owner.id, owner.path, owner.symbol
                ),
            ));
        }
    }

    for exemption in &manifest.exemptions {
        if !EXEMPTION_REASONS.contains(&exemption.reason.as_str()) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "exemption {}#{} has unsupported reason {:?}",
                    exemption.path, exemption.symbol, exemption.reason
                ),
            ));
        }
        if exemption.reason == "child-digest-helper" && exemption.covered_by == "none" {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "child digest exemption {}#{} must name a registered parent in covered_by",
                    exemption.path, exemption.symbol
                ),
            ));
        }
        if exemption.reason != "child-digest-helper" && exemption.covered_by != "none" {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "exemption {}#{} reason {:?} describes an independent exclusion and may not name parent {:?}",
                    exemption.path, exemption.symbol, exemption.reason, exemption.covered_by
                ),
            ));
        }
        if exemption.covered_by != "none" {
            if !manifest.required_ids.contains(&exemption.covered_by) {
                violations.push(identity_violation(
                    "<repo>",
                    format!(
                        "exemption {}#{} names covered_by {:?}, which is absent from required_ids",
                        exemption.path, exemption.symbol, exemption.covered_by
                    ),
                ));
            } else if !registered_ids.contains(&exemption.covered_by) {
                violations.push(identity_violation(
                    "<repo>",
                    format!(
                        "exemption {}#{} names unregistered parent {:?}",
                        exemption.path, exemption.symbol, exemption.covered_by
                    ),
                ));
            } else {
                let declaration_covers = declarations
                    .iter()
                    .filter(|declaration| declaration.id == exemption.covered_by)
                    .flat_map(declaration_schema_targets)
                    .any(|(path, symbol)| path == exemption.path && symbol == exemption.symbol);
                let external_covers = manifest.external_owners.iter().any(|owner| {
                    owner.id == exemption.covered_by
                        && external_owner_covers_exemption(root, owner, exemption)
                });
                if !declaration_covers && !external_covers {
                    violations.push(identity_violation(
                        "<repo>",
                        format!(
                            "exemption {}#{} names parent {:?}, but that parent does not include the child target in its exact schema closure",
                            exemption.path, exemption.symbol, exemption.covered_by
                        ),
                    ));
                }
            }
        }
        if let Some(detail) = authority_target_error(root, &exemption.path, &exemption.symbol) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "exemption target {}#{} is invalid: {detail}",
                    exemption.path, exemption.symbol
                ),
            ));
        }
        if !matches!(
            exemption.reason.as_str(),
            "child-digest-helper" | "ephemeral-staging-key"
        ) && !candidates.iter().any(|candidate| {
            authority_target_matches(candidate, &exemption.path, &exemption.symbol)
        }) {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "exemption {}#{} is stale: the target is not a discovered identity producer",
                    exemption.path, exemption.symbol
                ),
            ));
        }
    }

    let declaration_authority_targets = declarations
        .iter()
        .flat_map(|declaration| {
            declaration_authority_targets(declaration)
                .into_iter()
                .map(move |(path, symbol)| (declaration.id.as_str(), path, symbol))
        })
        .collect::<Vec<_>>();
    let declaration_schema_targets = declarations
        .iter()
        .flat_map(|declaration| {
            declaration_schema_targets(declaration)
                .into_iter()
                .map(move |(path, symbol)| (declaration.id.as_str(), path, symbol))
        })
        .collect::<Vec<_>>();
    for candidate in candidates {
        let authority_owners = declaration_authority_targets
            .iter()
            .filter_map(|(id, path, symbol)| {
                (manifest.required_ids.contains(*id)
                    && authority_target_matches(candidate, path, symbol))
                .then_some(*id)
            })
            .chain(manifest.external_owners.iter().filter_map(|owner| {
                (manifest.required_ids.contains(&owner.id)
                    && authority_target_matches(candidate, &owner.path, &owner.symbol))
                .then_some(owner.id.as_str())
            }))
            .collect::<BTreeSet<_>>();
        let schema_owners = declaration_schema_targets
            .iter()
            .filter_map(|(id, path, symbol)| {
                (manifest.required_ids.contains(*id)
                    && authority_target_matches(candidate, path, symbol))
                .then_some(*id)
            })
            .collect::<BTreeSet<_>>();
        let owners = if authority_owners.is_empty() {
            &schema_owners
        } else {
            &authority_owners
        };
        let exemptions = manifest
            .exemptions
            .iter()
            .filter(|exemption| {
                authority_target_matches(candidate, &exemption.path, &exemption.symbol)
            })
            .collect::<Vec<_>>();
        if owners.is_empty() && exemptions.is_empty() {
            violations.push(identity_violation(
                &candidate.path,
                format!(
                    "discovered identity producer {}#{} is unregistered (identity signal {:?}, durable sink {:?}); add a required owner declaration/external owner or a structured exemption",
                    candidate.path,
                    candidate.symbol,
                    candidate.identity_signal,
                    candidate.sink_signal
                ),
            ));
        } else if authority_owners.len() > 1 {
            violations.push(identity_violation(
                &candidate.path,
                format!(
                    "discovered identity producer {}#{} is claimed by multiple owners {:?}",
                    candidate.path, candidate.symbol, authority_owners
                ),
            ));
        } else if !owners.is_empty() && !exemptions.is_empty() {
            violations.push(identity_violation(
                &candidate.path,
                format!(
                    "discovered identity producer {}#{} is both owned by {:?} and exempted; remove the redundant exemption",
                    candidate.path, candidate.symbol, owners
                ),
            ));
        }
    }
    violations
}

fn identity_authority_violations(root: &Path, declarations: &[IdentityDecl]) -> Vec<Violation> {
    let manifest = match load_authority_manifest(root) {
        Ok(manifest) => manifest,
        Err(violations) => return violations,
    };
    let candidates = match discover_identity_candidates(root, &manifest) {
        Ok(candidates) => candidates,
        Err(violations) => return violations,
    };
    authority_violations_against(root, declarations, &manifest, &candidates)
}

fn coupling_surfaces(root: &Path) -> Result<BTreeMap<String, CouplingSurface>, Vec<Violation>> {
    let text = read_repo_utf8(root, "golden-couplings.json", "golden coupling registry")
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    let parsed = JsonParser::new(&text).finish().map_err(|detail| {
        vec![identity_violation(
            "<repo>",
            format!("golden-couplings.json is not strict JSON: {detail}"),
        )]
    })?;
    let object = strict_json_object(&parsed, "golden-couplings.json")
        .and_then(|object| {
            strict_json_keys(
                object,
                &["schema", "note", "surfaces", "goldens"],
                "golden-couplings.json",
            )?;
            Ok(object)
        })
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    let schema = strict_json_field(object, "schema", "golden-couplings.json")
        .and_then(|value| strict_json_string(value, "golden coupling schema"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    if schema != GOLDEN_COUPLING_SCHEMA {
        return Err(vec![identity_violation(
            "<repo>",
            format!("golden-couplings.json must declare schema {GOLDEN_COUPLING_SCHEMA:?}"),
        )]);
    }
    strict_json_field(object, "note", "golden-couplings.json")
        .and_then(|value| strict_json_string(value, "golden coupling note"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    strict_json_field(object, "goldens", "golden-couplings.json")
        .and_then(|value| strict_json_array(value, "golden coupling goldens"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    let rows = strict_json_field(object, "surfaces", "golden-couplings.json")
        .and_then(|value| strict_json_array(value, "golden coupling surfaces"))
        .map_err(|detail| vec![identity_violation("<repo>", detail)])?;
    let mut surfaces = BTreeMap::new();
    let mut violations = Vec::new();
    for (index, value) in rows.iter().enumerate() {
        let context = format!("golden-couplings.json surface row {}", index + 1);
        let parsed = (|| -> Result<(String, CouplingSurface), String> {
            let row = strict_json_object(value, &context)?;
            let keys = row.keys().map(String::as_str).collect::<BTreeSet<_>>();
            let base = ["id", "file", "const", "version"]
                .into_iter()
                .collect::<BTreeSet<_>>();
            let with_domain = ["id", "file", "const", "version", "domain_const", "domain"]
                .into_iter()
                .collect::<BTreeSet<_>>();
            let with_fingerprint = [
                "id",
                "file",
                "const",
                "version",
                "domain_const",
                "domain",
                "schema_fingerprint",
            ]
            .into_iter()
            .collect::<BTreeSet<_>>();
            let external = [
                "id",
                "file",
                "symbol",
                "version",
                "domain",
                "schema_fingerprint",
            ]
            .into_iter()
            .collect::<BTreeSet<_>>();
            if keys != base && keys != with_domain && keys != with_fingerprint && keys != external {
                return Err(format!(
                    "{context} has noncanonical keys {keys:?}; use one complete Rust const surface or external symbol surface"
                ));
            }
            let string = |key: &str| {
                strict_json_field(row, key, &context)
                    .and_then(|value| strict_json_string(value, &format!("{context} {key}")))
            };
            let id = string("id")?;
            let file = string("file")?;
            let version = strict_json_field(row, "version", &context)
                .and_then(|value| strict_json_u32(value, &format!("{context} version")))?;
            if !id.contains(':') || id.chars().any(char::is_whitespace) {
                return Err(format!("{context} id {id:?} is not canonical"));
            }
            if !safe_relative(file) {
                return Err(format!(
                    "{context} has unsafe/noncanonical source path {file:?}"
                ));
            }
            if keys == external {
                let symbol = string("symbol")?;
                let domain = string("domain")?;
                let schema_fingerprint = string("schema_fingerprint")?;
                if !authority_symbol_is_canonical(symbol) {
                    return Err(format!(
                        "{context} has noncanonical external symbol {file}#{symbol}"
                    ));
                }
                if !domain_carries_version(domain, version) {
                    return Err(format!(
                        "{context} external domain {domain:?} must carry exact version {version}"
                    ));
                }
                if !canonical_schema_fingerprint(schema_fingerprint, version, 32) {
                    return Err(format!(
                        "{context} schema_fingerprint must be a canonical version-{version} BLAKE3-256 digest"
                    ));
                }
                return Ok((
                    id.to_string(),
                    CouplingSurface::External {
                        file: file.to_string(),
                        symbol: symbol.to_string(),
                        version,
                        domain: domain.to_string(),
                        schema_fingerprint: schema_fingerprint.to_string(),
                    },
                ));
            }
            let version_const = string("const")?;
            if !canonical_symbol(version_const) {
                return Err(format!(
                    "{context} has noncanonical Rust source {file}#{version_const}"
                ));
            }
            let domain_const = row
                .get("domain_const")
                .map(|value| strict_json_string(value, &format!("{context} domain_const")))
                .transpose()?
                .map(str::to_string);
            let domain = row
                .get("domain")
                .map(|value| strict_json_string(value, &format!("{context} domain")))
                .transpose()?
                .map(str::to_string);
            let schema_fingerprint = row
                .get("schema_fingerprint")
                .map(|value| strict_json_string(value, &format!("{context} schema_fingerprint")))
                .transpose()?
                .map(str::to_string);
            if domain_const.as_deref().is_some_and(|value| {
                parse_schema_constant_reference(value)
                    .ok()
                    .is_none_or(|constant| constant.canonical() != value)
            }) {
                return Err(format!("{context} domain_const is not canonical"));
            }
            if schema_fingerprint
                .as_deref()
                .is_some_and(|fingerprint| !canonical_schema_fingerprint(fingerprint, version, 32))
            {
                return Err(format!(
                    "{context} schema_fingerprint must be a canonical version-{version} BLAKE3-256 digest"
                ));
            }
            Ok((
                id.to_string(),
                CouplingSurface::Rust {
                    file: file.to_string(),
                    version_const: version_const.to_string(),
                    version,
                    domain_const,
                    domain,
                    schema_fingerprint,
                },
            ))
        })();
        match parsed {
            Ok((id, surface)) => {
                if surfaces.insert(id.clone(), surface).is_some() {
                    violations.push(identity_violation(
                        "<repo>",
                        format!("golden-couplings.json contains duplicate surface id {id:?}"),
                    ));
                }
            }
            Err(detail) => violations.push(identity_violation("<repo>", detail)),
        }
    }
    if violations.is_empty() {
        Ok(surfaces)
    } else {
        Err(violations)
    }
}

fn validate_owner_items(
    decl: &IdentityDecl,
    text: &str,
    index: &RustSourceIndex<'_>,
    references: &IdentityReferenceCache,
) -> Vec<Violation> {
    let mut violations = Vec::new();
    let fail = |detail| identity_violation(&decl.owner, format!("identity {}: {detail}", decl.id));

    if source_const_u32_with_scopes(text, &decl.version_const, &index.module_scopes)
        != Some(decl.version)
    {
        violations.push(fail(format!(
            "{} must declare u32 version {}; schema changes require a coupling bump",
            decl.version_const, decl.version
        )));
    }
    if let Some(domain_const) = &decl.domain_const {
        let constant = parse_schema_constant_reference(domain_const).ok();
        let value = constant.as_ref().and_then(|constant| {
            let path = constant.path.as_deref().unwrap_or(&decl.owner);
            references
                .constant(path, &constant.symbol)
                .ok()
                .and_then(|reference| normalized_string_constant_value(&reference.declaration))
        });
        if value != Some(decl.domain.as_str()) {
            violations.push(fail(format!(
                "{domain_const} must declare exact domain {:?}",
                decl.domain
            )));
        }
    }
    for (role, symbol) in [
        ("encoder", &decl.encoder),
        ("transport guard", &decl.transport_guard),
    ] {
        if !has_function_with_index(text, index, symbol) {
            violations.push(fail(format!(
                "owner source does not define exact {role} function {symbol:?}"
            )));
        }
    }
    for helper in &decl.encoder_helpers {
        if !has_function_with_index(text, index, helper) {
            violations.push(fail(format!(
                "owner source does not define exact encoder helper function {helper:?}"
            )));
        }
    }
    for function in decl
        .schema_functions
        .iter()
        .filter(|function| function.path.is_none())
    {
        if !has_function_with_index(text, index, &function.symbol) {
            violations.push(fail(format!(
                "owner source does not define one runtime-active schema function {:?}",
                function.symbol
            )));
        }
    }
    if !symbol_body_has_no_rest_pattern(text, &decl.field_guard) {
        violations.push(fail(format!(
            "field guard {:?} is missing or contains `..`; it must exhaustively destructure the owner type",
            decl.field_guard
        )));
    }
    if !guard_destructures_sources_with_scopes(
        text,
        &decl.field_guard,
        &decl.sources,
        &index.module_scopes,
    ) {
        violations.push(fail(format!(
            "field guard {:?} must explicitly destructure every declared source type {:?}",
            decl.field_guard, decl.sources
        )));
    }
    violations
}

fn validate_source_classification(
    decl: &IdentityDecl,
    text: &str,
    index: &RustSourceIndex<'_>,
) -> Vec<Violation> {
    let mut violations = Vec::new();
    let fail = |detail| identity_violation(&decl.owner, format!("identity {}: {detail}", decl.id));

    let declared_sources = decl.sources.iter().cloned().collect::<BTreeSet<_>>();
    let mut classified = BTreeSet::new();
    for field in &decl.source_fields {
        if !classified.insert(field.qualified.clone()) {
            violations.push(fail(format!(
                "source field {:?} is classified twice",
                field.qualified
            )));
        }
        let source = field.qualified.split_once('.').map(|(source, _)| source);
        if source.is_none_or(|source| !declared_sources.contains(source)) {
            violations.push(fail(format!(
                "source field {:?} names an undeclared source struct",
                field.qualified
            )));
        }
    }
    for source in &decl.sources {
        let Some(shape) = source_shape_with_scopes(text, source, &index.module_scopes) else {
            violations.push(fail(format!("source struct/enum {source:?} was not found")));
            continue;
        };
        if let SourceShape::Enum { tuple_variants, .. } = &shape
            && !tuple_variants.is_empty()
        {
            violations.push(fail(format!(
                "source enum {source:?} has tuple-payload variants {tuple_variants:?}; tuple positions need an explicit variant-qualified identity schema before this source can be registered"
            )));
            continue;
        }
        let actual = shape.fields();
        let expected = decl
            .source_fields
            .iter()
            .filter_map(|field| {
                field
                    .qualified
                    .strip_prefix(&format!("{source}."))
                    .map(str::to_string)
            })
            .collect::<BTreeSet<_>>();
        for field in actual.difference(&expected) {
            violations.push(fail(format!(
                "unclassified source field {source}.{field}; classify it as semantic, derived, or nonsemantic"
            )));
        }
        for field in expected.difference(actual) {
            violations.push(fail(format!(
                "declared source field {source}.{field} does not exist"
            )));
        }
    }

    let classes = decl
        .source_fields
        .iter()
        .map(|field| (field.qualified.as_str(), field.class))
        .collect::<BTreeMap<_, _>>();
    let mut bound_sources = BTreeSet::new();
    let mut bound_semantics = BTreeSet::new();
    for binding in &decl.source_bindings {
        if !bound_sources.insert(binding.source_field.as_str()) {
            violations.push(fail(format!(
                "semantic source field {:?} has multiple bindings",
                binding.source_field
            )));
        }
        if classes.get(binding.source_field.as_str()) != Some(&FieldClass::Semantic) {
            violations.push(fail(format!(
                "source binding {:?} must name a source field classified semantic",
                binding.source_field
            )));
        }
        for semantic in &binding.semantic_fields {
            if !bound_semantics.insert(semantic.as_str()) {
                violations.push(fail(format!(
                    "logical semantic field {semantic:?} is bound more than once"
                )));
            }
        }
    }
    for field in decl
        .source_fields
        .iter()
        .filter(|field| field.class == FieldClass::Semantic)
    {
        if !bound_sources.contains(field.qualified.as_str()) {
            violations.push(fail(format!(
                "semantic source field {:?} has no source_bindings entry",
                field.qualified
            )));
        }
    }
    let mut all_semantics = bound_semantics;
    for semantic in &decl.external_semantic_fields {
        if !all_semantics.insert(semantic.as_str()) {
            violations.push(fail(format!(
                "external semantic field {semantic:?} duplicates a source binding"
            )));
        }
    }
    let declared_semantics = decl
        .semantic_fields
        .iter()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    for semantic in declared_semantics.difference(&all_semantics) {
        violations.push(fail(format!(
            "semantic field {semantic:?} has no source binding or external input classification"
        )));
    }
    for semantic in all_semantics.difference(&declared_semantics) {
        violations.push(fail(format!(
            "source/external binding names undeclared semantic field {semantic:?}"
        )));
    }
    violations
}

fn validate_mutations(decl: &IdentityDecl, references: &IdentityReferenceCache) -> Vec<Violation> {
    let mut violations = Vec::new();
    let fail = |detail| identity_violation(&decl.owner, format!("identity {}: {detail}", decl.id));

    let semantic = decl
        .semantic_fields
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    if semantic.len() != decl.semantic_fields.len() {
        violations.push(fail("semantic_fields contains duplicates".to_string()));
    }
    let mutation_fields = decl
        .mutations
        .iter()
        .map(|mutation| mutation.field.clone())
        .collect::<BTreeSet<_>>();
    if mutation_fields.len() != decl.mutations.len() {
        violations.push(fail(
            "mutations contains duplicate semantic fields".to_string(),
        ));
    }
    for field in semantic.difference(&mutation_fields) {
        violations.push(fail(format!(
            "semantic field {field:?} has no mutation test"
        )));
    }
    for field in mutation_fields.difference(&semantic) {
        violations.push(fail(format!(
            "mutation names unknown semantic field {field:?}"
        )));
    }
    let nonsemantic_sources = decl
        .source_fields
        .iter()
        .filter(|field| field.class == FieldClass::Nonsemantic)
        .map(|field| field.qualified.clone())
        .collect::<BTreeSet<_>>();
    let reasoned_exclusions = decl
        .excluded_fields
        .iter()
        .map(|(field, _)| field.clone())
        .collect::<BTreeSet<_>>();
    if reasoned_exclusions.len() != decl.excluded_fields.len() {
        violations.push(fail(
            "excluded_fields contains duplicate field names".to_string(),
        ));
    }
    for field in nonsemantic_sources.intersection(&reasoned_exclusions) {
        violations.push(fail(format!(
            "nonsemantic field {field:?} is declared both as a source classification and a reasoned external exclusion"
        )));
    }
    let expected_nonmovement = nonsemantic_sources
        .union(&reasoned_exclusions)
        .cloned()
        .collect::<BTreeSet<_>>();
    let mutation_counts = decl.nonsemantic_mutations.iter().fold(
        BTreeMap::<String, usize>::new(),
        |mut counts, mutation| {
            *counts.entry(mutation.field.clone()).or_default() += 1;
            counts
        },
    );
    for field in &expected_nonmovement {
        match mutation_counts.get(field).copied().unwrap_or(0) {
            1 => {}
            0 => violations.push(fail(format!(
                "nonsemantic/excluded field {field:?} has no non-movement test"
            ))),
            count => violations.push(fail(format!(
                "nonsemantic/excluded field {field:?} has {count} non-movement tests; exactly one is required"
            ))),
        }
    }
    for field in mutation_counts.keys() {
        if !expected_nonmovement.contains(field) {
            violations.push(fail(format!(
                "nonsemantic mutation names unknown source/excluded field {field:?}"
            )));
        }
    }
    for mutation in decl
        .mutations
        .iter()
        .chain(decl.nonsemantic_mutations.iter())
    {
        if let Err(detail) = references.test(&mutation.test_path, &mutation.test_symbol) {
            violations.push(fail(format!(
                "mutation target {}#{} must resolve to an exact #[test] function: {detail}",
                mutation.test_path, mutation.test_symbol,
            )));
        }
    }
    let Some((guard_path, guard_symbol)) = decl.version_guard.split_once('#') else {
        violations.push(fail("version_guard must be path#test_symbol".to_string()));
        return violations;
    };
    if let Err(detail) = references.test(guard_path, guard_symbol) {
        violations.push(fail(format!(
            "version guard {:?} must resolve to an exact #[test] function: {detail}",
            decl.version_guard,
        )));
    }
    if decl.encoding == "canonical-transport-exact-bits" && decl.transport_guard == "none" {
        violations.push(fail(
            "canonical transport encoding needs a closed exact-bit transport guard".to_string(),
        ));
    }
    violations
}

fn validate_coupling(
    decl: &IdentityDecl,
    surfaces: &BTreeMap<String, CouplingSurface>,
) -> Vec<Violation> {
    let fail = |detail| identity_violation(&decl.owner, format!("identity {}: {detail}", decl.id));
    let expected = CouplingSurface::Rust {
        file: decl.owner.clone(),
        version_const: decl.version_const.clone(),
        version: decl.version,
        domain_const: decl.domain_const.clone(),
        domain: Some(decl.domain.clone()),
        schema_fingerprint: Some(decl.schema_fingerprint.clone()),
    };
    let expected_row = format!(
        "{{\"id\": \"{}\", \"file\": \"{}\", \"const\": \"{}\", \"version\": {}, \"domain_const\": \"{}\", \"domain\": \"{}\", \"schema_fingerprint\": \"{}\"}}",
        json_escape(&decl.coupling_surface),
        json_escape(&decl.owner),
        json_escape(&decl.version_const),
        decl.version,
        json_escape(
            decl.domain_const
                .as_deref()
                .expect("identity declarations require domain_const")
        ),
        json_escape(&decl.domain),
        json_escape(&decl.schema_fingerprint),
    );
    let Some(surface) = surfaces.get(&decl.coupling_surface) else {
        return vec![fail(format!(
            "golden surface {:?} is missing from golden-couplings.json; add exact row {expected_row}",
            decl.coupling_surface
        ))];
    };
    if *surface == expected {
        Vec::new()
    } else {
        vec![fail(format!(
            "golden surface {:?} must exactly match owner/version/domain declaration: replace it with exact row {expected_row}; found {surface:?}",
            decl.coupling_surface
        ))]
    }
}

fn validate_external_couplings(
    root: &Path,
    manifest: &AuthorityManifest,
    surfaces: &BTreeMap<String, CouplingSurface>,
) -> (BTreeMap<String, u32>, Vec<Violation>) {
    let mut versions = BTreeMap::new();
    let mut violations = Vec::new();
    let external_ids = manifest
        .external_owners
        .iter()
        .map(|owner| owner.id.as_str())
        .collect::<BTreeSet<_>>();

    for owner in &manifest.external_owners {
        let fingerprint = match external_owner_schema_fingerprint(root, owner, &manifest.exemptions)
        {
            Ok(fingerprint) => fingerprint,
            Err(detail) => {
                violations.push(identity_violation(
                    &owner.path,
                    format!(
                        "identity {}: external golden coupling fingerprint is unavailable: {detail}",
                        owner.id
                    ),
                ));
                continue;
            }
        };
        let expected = CouplingSurface::External {
            file: owner.path.clone(),
            symbol: owner.symbol.clone(),
            version: owner.version,
            domain: owner.domain.clone(),
            schema_fingerprint: fingerprint.clone(),
        };
        let expected_row = format!(
            "{{\"id\": \"{}\", \"file\": \"{}\", \"symbol\": \"{}\", \"version\": {}, \"domain\": \"{}\", \"schema_fingerprint\": \"{}\"}}",
            json_escape(&owner.id),
            json_escape(&owner.path),
            json_escape(&owner.symbol),
            owner.version,
            json_escape(&owner.domain),
            json_escape(&fingerprint),
        );
        match surfaces.get(&owner.id) {
            None => violations.push(identity_violation(
                &owner.path,
                format!(
                    "identity {}: external golden surface is missing from golden-couplings.json; add exact row {expected_row}",
                    owner.id
                ),
            )),
            Some(surface) if *surface != expected => violations.push(identity_violation(
                &owner.path,
                format!(
                    "identity {}: external golden surface must exactly match authority path/symbol/version/domain declaration and producer fingerprint: replace it with exact row {expected_row}; found {surface:?}",
                    owner.id
                ),
            )),
            Some(_) => {
                versions.insert(owner.id.clone(), owner.version);
            }
        }
    }

    for (id, surface) in surfaces {
        if matches!(surface, CouplingSurface::External { .. })
            && !external_ids.contains(id.as_str())
        {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "golden-couplings.json contains stale external surface {id:?} with no matching identity-authorities.json external owner"
                ),
            ));
        }
    }

    (versions, violations)
}

pub(super) fn external_coupling_versions(
    root: &Path,
) -> Result<BTreeMap<String, u32>, Vec<Violation>> {
    let manifest = load_authority_manifest(root)?;
    let surfaces = coupling_surfaces(root)?;
    let (versions, violations) = validate_external_couplings(root, &manifest, &surfaces);
    if violations.is_empty() {
        Ok(versions)
    } else {
        Err(violations)
    }
}

fn validate_declaration(
    decl: &IdentityDecl,
    text: &str,
    index: &RustSourceIndex<'_>,
    references: &IdentityReferenceCache,
) -> Vec<Violation> {
    let mut violations = validate_owner_items(decl, text, index, references);
    violations.extend(validate_source_classification(decl, text, index));
    violations.extend(validate_mutations(decl, references));
    violations
}

fn load_declarations(root: &Path) -> (Vec<IdentityDecl>, Vec<Violation>) {
    struct PendingSource {
        owner: String,
        text: String,
        declarations: Vec<IdentityDecl>,
        violations: Vec<Violation>,
    }

    let mut declarations = Vec::new();
    let mut violations = Vec::new();
    let mut pending = Vec::new();
    let exemptions = load_authority_manifest(root)
        .map(|manifest| manifest.exemptions)
        .unwrap_or_default();
    let (sources, source_violations) = source_files(root);
    violations.extend(source_violations);
    for path in sources {
        let owner = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .display()
            .to_string();
        let metadata = match std::fs::metadata(&path) {
            Ok(metadata) => metadata,
            Err(error) => {
                violations.push(identity_violation(
                    &owner,
                    format!("identity declaration scan cannot inspect Rust source: {error}"),
                ));
                continue;
            }
        };
        if metadata.len() > MAX_SOURCE_BYTES {
            violations.push(identity_violation(
                &owner,
                format!("source exceeds identity scan cap of {MAX_SOURCE_BYTES} bytes"),
            ));
            continue;
        }
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => {
                violations.push(identity_violation(
                    &owner,
                    format!("identity declaration scan cannot read Rust source: {error}"),
                ));
                continue;
            }
        };
        if identity_marker_lines(&text).is_empty() {
            continue;
        }
        let module_scopes = rust_owner_scopes(&text, "mod", rust_module_owner);
        let (found, parse_violations) =
            declaration_blocks_with_scopes(&owner, &text, &module_scopes);
        pending.push(PendingSource {
            owner,
            text,
            declarations: found,
            violations: parse_violations,
        });
    }
    let inline_sources = pending
        .iter()
        .map(|source| (source.owner.clone(), source.text.as_str()))
        .collect::<BTreeMap<_, _>>();
    let references = IdentityReferenceCache::build(
        root,
        pending.iter().flat_map(|source| source.declarations.iter()),
        &inline_sources,
    );
    drop(inline_sources);
    for mut source in pending {
        let source_index = RustSourceIndex::new(&source.text);
        for declaration in &mut source.declarations {
            match identity_schema_base_hash_with_index(
                root,
                declaration,
                &source.text,
                &source_index,
                &references,
                &exemptions,
            ) {
                Ok(hash) => declaration.schema_base_hash = Some(hash),
                Err(detail) => {
                    source
                        .violations
                        .push(identity_violation(&source.owner, detail));
                }
            }
            source.violations.extend(validate_declaration(
                declaration,
                &source.text,
                &source_index,
                &references,
            ));
        }
        declarations.append(&mut source.declarations);
        violations.append(&mut source.violations);
    }
    declarations.sort_by(|left, right| left.id.cmp(&right.id));
    if declarations.len() > MAX_DECLARATIONS {
        violations.push(identity_violation(
            "<repo>",
            format!(
                "identity declaration count {} exceeds cap {MAX_DECLARATIONS}",
                declarations.len()
            ),
        ));
    }
    let mut duplicate_ids = false;
    for pair in declarations.windows(2) {
        if pair[0].id == pair[1].id {
            duplicate_ids = true;
            violations.push(identity_violation(
                "<repo>",
                format!("duplicate semantic identity id {:?}", pair[0].id),
            ));
        }
    }
    let fingerprints_resolved = if duplicate_ids {
        false
    } else {
        let dependency_violations = resolve_schema_fingerprints(&mut declarations);
        let resolved = dependency_violations.is_empty();
        violations.extend(dependency_violations);
        resolved
    };
    if fingerprints_resolved {
        match coupling_surfaces(root) {
            Ok(surfaces) => {
                for declaration in &declarations {
                    violations.extend(validate_coupling(declaration, &surfaces));
                }
                if let Ok(manifest) = load_authority_manifest(root) {
                    let (_, external_violations) =
                        validate_external_couplings(root, &manifest, &surfaces);
                    violations.extend(external_violations);
                }
            }
            Err(coupling_violations) => violations.extend(coupling_violations),
        }
    }
    (declarations, violations)
}

fn render_string_array(out: &mut String, values: impl IntoIterator<Item = String>) {
    out.push('[');
    for (index, value) in values.into_iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\"{}\"", json_escape(&value));
    }
    out.push(']');
}

fn render_registry(
    root: &Path,
    declarations: &[IdentityDecl],
    external_owners: &[ExternalOwner],
    exemptions: &[IdentityExemption],
) -> Result<String, String> {
    let mut out = format!(
        "{{\n\"schema\": \"frankensim-identity-schemas-v1\",\n\"schema_fingerprint_algorithm\": \"{SCHEMA_FINGERPRINT_ALGORITHM}\",\n\"identities\": [\n"
    );
    for (index, declaration) in declarations.iter().enumerate() {
        if index > 0 {
            out.push_str(",\n");
        }
        let _ = write!(
            out,
            "{{\"id\":\"{}\",\"owner\":\"{}\",\"version_const\":\"{}\",\"version\":{},\"domain\":\"{}\",\"encoder\":\"{}\",\"digest\":\"{}\",\"encoding\":\"{}\",\"schema_fingerprint\":\"{}\",\"coupling_surface\":\"{}\",\"sources\":",
            json_escape(&declaration.id),
            json_escape(&declaration.owner),
            json_escape(&declaration.version_const),
            declaration.version,
            json_escape(&declaration.domain),
            json_escape(&declaration.encoder),
            json_escape(&declaration.digest),
            json_escape(&declaration.encoding),
            json_escape(&declaration.schema_fingerprint),
            json_escape(&declaration.coupling_surface),
        );
        render_string_array(&mut out, declaration.sources.clone());
        out.push_str(",\"encoder_helpers\":");
        render_string_array(&mut out, declaration.encoder_helpers.clone());
        out.push_str(",\"schema_functions\":");
        render_string_array(
            &mut out,
            declaration
                .schema_functions
                .iter()
                .map(SchemaFunction::canonical),
        );
        out.push_str(",\"schema_constants\":");
        render_string_array(
            &mut out,
            declaration
                .schema_constants
                .iter()
                .map(SchemaConstant::canonical),
        );
        out.push_str(",\"schema_dependencies\":");
        render_string_array(&mut out, declaration.schema_dependencies.clone());
        out.push_str(",\"source_fields\":");
        render_string_array(
            &mut out,
            declaration.source_fields.iter().map(|field| {
                field.reason.as_ref().map_or_else(
                    || format!("{}:{}", field.qualified, field.class.name()),
                    |reason| format!("{}:{}:{reason}", field.qualified, field.class.name()),
                )
            }),
        );
        out.push_str(",\"source_bindings\":");
        render_string_array(
            &mut out,
            declaration.source_bindings.iter().map(|binding| {
                format!(
                    "{}>{}",
                    binding.source_field,
                    binding.semantic_fields.join("+")
                )
            }),
        );
        out.push_str(",\"external_semantic_fields\":");
        render_string_array(&mut out, declaration.external_semantic_fields.clone());
        out.push_str(",\"semantic_fields\":");
        render_string_array(&mut out, declaration.semantic_fields.clone());
        out.push_str(",\"excluded_fields\":");
        render_string_array(
            &mut out,
            declaration
                .excluded_fields
                .iter()
                .map(|(field, reason)| format!("{field}:{reason}")),
        );
        out.push_str(",\"consumers\":");
        let mut consumers = declaration.consumers.clone();
        consumers.sort();
        render_string_array(&mut out, consumers);
        out.push_str(",\"mutations\":");
        render_string_array(
            &mut out,
            declaration.mutations.iter().map(|mutation| {
                format!(
                    "{}:{}#{}",
                    mutation.field, mutation.test_path, mutation.test_symbol
                )
            }),
        );
        out.push_str(",\"nonsemantic_mutations\":");
        render_string_array(
            &mut out,
            declaration.nonsemantic_mutations.iter().map(|mutation| {
                format!(
                    "{}:{}#{}",
                    mutation.field, mutation.test_path, mutation.test_symbol
                )
            }),
        );
        let _ = write!(
            out,
            ",\"field_guard\":\"{}\",\"transport_guard\":\"{}\",\"version_guard\":\"{}\"}}",
            json_escape(&declaration.field_guard),
            json_escape(&declaration.transport_guard),
            json_escape(&declaration.version_guard),
        );
    }
    out.push_str("\n],\n\"external_authorities\": [\n");
    for (index, owner) in external_owners.iter().enumerate() {
        if index > 0 {
            out.push_str(",\n");
        }
        let fingerprint = external_owner_schema_fingerprint(root, owner, exemptions)?;
        let _ = write!(
            out,
            "{{\"id\":\"{}\",\"owner\":\"{}\",\"symbol\":\"{}\",\"version\":{},\"domain\":\"{}\",\"schema_fingerprint\":\"{}\",\"schema_children\":",
            json_escape(&owner.id),
            json_escape(&owner.path),
            json_escape(&owner.symbol),
            owner.version,
            json_escape(&owner.domain),
            json_escape(&fingerprint),
        );
        render_string_array(
            &mut out,
            exemptions
                .iter()
                .filter(|exemption| {
                    exemption.covered_by == owner.id && exemption.path == owner.path
                })
                .map(|exemption| {
                    format!(
                        "{}#{}:{}",
                        exemption.path, exemption.symbol, exemption.reason
                    )
                }),
        );
        out.push('}');
    }
    out.push_str("\n],\n\"exemptions\": [\n");
    let mut rendered_exemptions = exemptions.to_vec();
    rendered_exemptions.sort_by(|left, right| {
        (&left.path, &left.symbol, &left.reason, &left.covered_by).cmp(&(
            &right.path,
            &right.symbol,
            &right.reason,
            &right.covered_by,
        ))
    });
    for (index, exemption) in rendered_exemptions.iter().enumerate() {
        if index > 0 {
            out.push_str(",\n");
        }
        let fingerprint = exemption_schema_fingerprint(root, exemption)?;
        let _ = write!(
            out,
            "{{\"path\":\"{}\",\"symbol\":\"{}\",\"reason\":\"{}\",\"covered_by\":\"{}\",\"schema_fingerprint\":\"{}\"}}",
            json_escape(&exemption.path),
            json_escape(&exemption.symbol),
            json_escape(&exemption.reason),
            json_escape(&exemption.covered_by),
            json_escape(&fingerprint),
        );
    }
    out.push_str("\n]\n}\n");
    Ok(out)
}

fn registry_fingerprints(text: &str) -> Result<BTreeMap<(String, u32), Option<String>>, String> {
    let parsed = JsonParser::new(text)
        .finish()
        .map_err(|detail| format!("identity registry is not strict JSON: {detail}"))?;
    let object = strict_json_object(&parsed, "identity registry")?;
    let keys = object.keys().map(String::as_str).collect::<BTreeSet<_>>();
    let minimal_legacy = ["schema", "identities"]
        .into_iter()
        .collect::<BTreeSet<_>>();
    let legacy = ["schema", "identities", "external_authorities"]
        .into_iter()
        .collect::<BTreeSet<_>>();
    let current = [
        "schema",
        "schema_fingerprint_algorithm",
        "identities",
        "external_authorities",
        "exemptions",
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    if keys != minimal_legacy && keys != legacy && keys != current {
        return Err(format!(
            "identity registry has noncanonical top-level keys {keys:?}"
        ));
    }
    let schema = strict_json_field(object, "schema", "identity registry")
        .and_then(|value| strict_json_string(value, "identity registry schema"))?;
    if schema != "frankensim-identity-schemas-v1" {
        return Err(format!(
            "identity registry schema {schema:?} is not supported"
        ));
    }
    let current_algorithm = if let Some(value) = object.get("schema_fingerprint_algorithm") {
        let algorithm = strict_json_string(value, "identity registry fingerprint algorithm")?;
        if algorithm != SCHEMA_FINGERPRINT_ALGORITHM {
            return Err(format!(
                "identity registry fingerprint algorithm {algorithm:?} is not supported"
            ));
        }
        true
    } else {
        false
    };
    let identities = strict_json_field(object, "identities", "identity registry")
        .and_then(|value| strict_json_array(value, "identity registry identities"))?;
    let empty = Vec::new();
    let external = object
        .get("external_authorities")
        .map(|value| strict_json_array(value, "identity registry external_authorities"))
        .transpose()?
        .unwrap_or(&empty);
    let mut fingerprints = BTreeMap::new();
    let mut ids = BTreeSet::new();
    for (section, rows) in [
        ("identities", identities),
        ("external_authorities", external),
    ] {
        for (index, value) in rows.iter().enumerate() {
            let context = format!("identity registry {section} row {}", index + 1);
            let row = strict_json_object(value, &context)?;
            let id = strict_json_field(row, "id", &context)
                .and_then(|value| strict_json_string(value, &format!("{context} id")))?;
            let version = strict_json_field(row, "version", &context)
                .and_then(|value| strict_json_u32(value, &format!("{context} version")))?;
            let fingerprint = row
                .get("schema_fingerprint")
                .map(|value| strict_json_string(value, &format!("{context} schema_fingerprint")))
                .transpose()?;
            if current_algorithm && fingerprint.is_none() {
                return Err(format!(
                    "{context} is missing schema_fingerprint under {SCHEMA_FINGERPRINT_ALGORITHM}"
                ));
            }
            if !id.contains(':') || id.chars().any(char::is_whitespace) {
                return Err(format!("{context} id {id:?} is not canonical"));
            }
            if let Some(fingerprint) = fingerprint {
                let canonical = if current_algorithm {
                    canonical_schema_fingerprint(fingerprint, version, 32)
                } else {
                    canonical_schema_fingerprint(fingerprint, version, 8)
                        || canonical_schema_fingerprint(fingerprint, version, 32)
                };
                if !canonical {
                    return Err(format!(
                        "{context} schema_fingerprint is not canonical for version {version}"
                    ));
                }
            }
            if !ids.insert(id.to_string()) {
                return Err(format!("identity registry duplicates id {id:?}"));
            }
            if fingerprints
                .insert((id.to_string(), version), fingerprint.map(str::to_string))
                .is_some()
            {
                return Err(format!(
                    "identity registry duplicates id/version {id:?} v{version}"
                ));
            }
        }
    }
    if let Some(value) = object.get("exemptions") {
        let rows = strict_json_array(value, "identity registry exemptions")?;
        let mut targets = BTreeSet::new();
        for (index, value) in rows.iter().enumerate() {
            let context = format!("identity registry exemptions row {}", index + 1);
            let row = strict_json_object(value, &context)?;
            strict_json_keys(
                row,
                &[
                    "path",
                    "symbol",
                    "reason",
                    "covered_by",
                    "schema_fingerprint",
                ],
                &context,
            )?;
            let string = |key: &str| {
                strict_json_field(row, key, &context)
                    .and_then(|value| strict_json_string(value, &format!("{context} {key}")))
            };
            let path = string("path")?;
            let symbol = string("symbol")?;
            string("reason")?;
            string("covered_by")?;
            let fingerprint = string("schema_fingerprint")?;
            if !safe_relative(path) || !authority_symbol_is_canonical(symbol) {
                return Err(format!(
                    "{context} target {path}#{symbol} is unsafe/noncanonical"
                ));
            }
            if !canonical_schema_fingerprint(fingerprint, 1, 32) {
                return Err(format!(
                    "{context} schema_fingerprint is not a canonical BLAKE3-256 digest"
                ));
            }
            if !targets.insert((path, symbol)) {
                return Err(format!(
                    "identity registry duplicates exemption {path}#{symbol}"
                ));
            }
        }
    }
    Ok(fingerprints)
}

fn git_registry(root: &Path, revision: &str) -> Result<String, String> {
    let output = Command::new("git")
        .args(["show", &format!("{revision}:{REGISTRY_FILE}")])
        .current_dir(root)
        .output()
        .map_err(|error| format!("cannot spawn git show for retained registry: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git show {revision}:{REGISTRY_FILE} failed with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    String::from_utf8(output.stdout)
        .map_err(|error| format!("retained registry at {revision} is not UTF-8: {error}"))
}

fn git_registry_history(root: &Path) -> Result<Vec<String>, String> {
    #[cfg(test)]
    if !root.join(".git").exists() {
        return Ok(Vec::new());
    }

    let output = Command::new("git")
        .args([
            "log",
            "--format=%H",
            "--diff-filter=AM",
            "--",
            REGISTRY_FILE,
        ])
        .current_dir(root)
        .output()
        .map_err(|error| format!("cannot spawn git log for schema history: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git log for {REGISTRY_FILE} failed with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    let revisions = String::from_utf8(output.stdout)
        .map_err(|error| format!("git log emitted non-UTF-8 revision data: {error}"))?;
    revisions
        .lines()
        .map(|revision| git_registry(root, revision))
        .collect()
}

fn schema_history_violations(
    root: &Path,
    expected_registry: &str,
) -> Result<Vec<Violation>, String> {
    let mut violations = Vec::new();
    let mut baselines = BTreeSet::new();
    for baseline in git_registry_history(root)? {
        if baselines.insert(baseline.clone()) {
            violations.extend(schema_history_against(expected_registry, &baseline));
        }
    }
    violations.sort_by(|left, right| {
        (&left.crate_name, &left.detail).cmp(&(&right.crate_name, &right.detail))
    });
    violations
        .dedup_by(|left, right| left.crate_name == right.crate_name && left.detail == right.detail);
    Ok(violations)
}

fn schema_history_against(expected_registry: &str, baseline: &str) -> Vec<Violation> {
    let previous = match registry_fingerprints(baseline) {
        Ok(previous) => previous,
        Err(detail) => {
            return vec![identity_violation(
                "<repo>",
                format!("retained identity registry is malformed: {detail}"),
            )];
        }
    };
    let current = match registry_fingerprints(expected_registry) {
        Ok(current) => current,
        Err(detail) => {
            return vec![identity_violation(
                "<repo>",
                format!("generated identity registry is malformed: {detail}"),
            )];
        }
    };
    let previous_by_id = previous
        .iter()
        .map(|((id, version), fingerprint)| (id, (*version, fingerprint)))
        .collect::<BTreeMap<_, _>>();
    let current_by_id = current
        .iter()
        .map(|((id, version), fingerprint)| (id, (*version, fingerprint)))
        .collect::<BTreeMap<_, _>>();
    let mut violations = Vec::new();
    for (id, (previous_version, previous_fingerprint)) in previous_by_id {
        let Some((version, fingerprint)) = current_by_id.get(id) else {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "retained identity {id} v{previous_version} was removed from the generated registry; keep a tombstone declaration or provide an explicit migration before removing replay authority"
                ),
            ));
            continue;
        };
        if *version < previous_version {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "identity {id} version regressed from retained v{previous_version} to v{version}"
                ),
            ));
        } else if *version == previous_version
            && let Some(previous_fingerprint) = previous_fingerprint
            && fingerprint.as_ref() != Some(previous_fingerprint)
        {
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "identity {id} changed schema fingerprint at retained version {version} ({previous_fingerprint} -> {}); bump the schema/domain/coupling version before regeneration",
                    fingerprint.as_deref().unwrap_or("<missing>")
                ),
            ));
        }
    }
    violations
}

pub(super) fn check_identities(root: &Path) -> Vec<Violation> {
    let (declarations, mut violations) = load_declarations(root);
    violations.extend(identity_authority_violations(root, &declarations));
    if !violations.is_empty() {
        return violations;
    }
    let manifest = match load_authority_manifest(root) {
        Ok(manifest) => manifest,
        Err(manifest_violations) => return manifest_violations,
    };
    let expected = match render_registry(
        root,
        &declarations,
        &manifest.external_owners,
        &manifest.exemptions,
    ) {
        Ok(expected) => expected,
        Err(detail) => return vec![identity_violation("<repo>", detail)],
    };
    match schema_history_violations(root, &expected) {
        Ok(history_violations) => violations.extend(history_violations),
        Err(detail) => violations.push(identity_violation(
            "<repo>",
            format!("schema history is unavailable: {detail}"),
        )),
    }
    match read_repo_utf8(root, REGISTRY_FILE, "generated identity registry") {
        Ok(actual) if actual == expected => {}
        Ok(actual) => {
            let mismatch = actual
                .bytes()
                .zip(expected.bytes())
                .position(|(left, right)| left != right)
                .unwrap_or_else(|| actual.len().min(expected.len()));
            violations.push(identity_violation(
                "<repo>",
                format!(
                    "{REGISTRY_FILE} is stale at byte {mismatch}; run `cargo run -p xtask -- generate-identities` after reviewing schema/domain coupling bumps"
                ),
            ));
        }
        Err(error) => violations.push(identity_violation(
            "<repo>",
            format!("{REGISTRY_FILE} is missing or unreadable: {error}"),
        )),
    }
    violations
}

pub(super) fn generate_identities(root: &Path) -> ExitCode {
    let (declarations, mut violations) = load_declarations(root);
    violations.extend(identity_authority_violations(root, &declarations));
    if !violations.is_empty() {
        for violation in violations {
            eprintln!("VIOLATION [{}] {}", violation.check, violation.detail);
        }
        return ExitCode::FAILURE;
    }
    let manifest = match load_authority_manifest(root) {
        Ok(manifest) => manifest,
        Err(violations) => {
            for violation in violations {
                eprintln!("VIOLATION [{}] {}", violation.check, violation.detail);
            }
            return ExitCode::FAILURE;
        }
    };
    let registry = match render_registry(
        root,
        &declarations,
        &manifest.external_owners,
        &manifest.exemptions,
    ) {
        Ok(registry) => registry,
        Err(detail) => {
            eprintln!("VIOLATION [semantic-identities] {detail}");
            return ExitCode::FAILURE;
        }
    };
    let history_violations = match schema_history_violations(root, &registry) {
        Ok(violations) => violations,
        Err(detail) => {
            eprintln!("VIOLATION [semantic-identities] schema history is unavailable: {detail}");
            return ExitCode::FAILURE;
        }
    };
    if !history_violations.is_empty() {
        for violation in history_violations {
            eprintln!("VIOLATION [{}] {}", violation.check, violation.detail);
        }
        return ExitCode::FAILURE;
    }
    if let Err(error) = std::fs::write(root.join(REGISTRY_FILE), registry) {
        eprintln!("error writing {REGISTRY_FILE}: {error}");
        return ExitCode::FAILURE;
    }
    eprintln!(
        "semantic identity registry updated: {} owner declarations -> {REGISTRY_FILE}",
        declarations.len()
    );
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_root(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "frankensim-identity-registry-{tag}-{}",
            std::process::id()
        ))
    }

    fn empty_authority_manifest() -> AuthorityManifest {
        AuthorityManifest {
            required_ids: BTreeSet::new(),
            external_owners: Vec::new(),
            exemptions: Vec::new(),
        }
    }

    fn discover_identity_candidates(root: &Path) -> Vec<IdentityCandidate> {
        super::discover_identity_candidates(root, &empty_authority_manifest())
            .expect("valid identity-discovery fixture")
    }

    fn discover_script_candidates(path: &str, text: &str) -> Vec<IdentityCandidate> {
        super::discover_script_candidates(path, text).expect("valid script fixture")
    }

    fn shell_function_blocks(text: &str) -> Vec<(String, &str)> {
        super::shell_function_blocks(text).expect("valid shell function fixture")
    }

    fn primary_script_blocks<'a>(path: &str, text: &'a str, symbol: &str) -> Vec<&'a str> {
        super::primary_script_blocks(path, text, symbol)
            .expect("valid primary-script fixture")
            .into_iter()
            .map(|block| block.text)
            .collect()
    }

    fn executable_script_views(path: &str, text: &str) -> (String, String) {
        super::executable_script_views(path, text).expect("valid executable-script fixture")
    }

    fn owner_source(extra_field: &str, semantic_fields: &str, mutations: &str) -> String {
        format!(
            r#"pub const MINI_VERSION: u32 = 1;
pub const MINI_DOMAIN: &str = "org.frankensim.mini.v1";
pub const MINI_TAG: u8 = 1;
pub const MINI_DECL: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=mini:identity",
    "version_const=MINI_VERSION",
    "version=1",
    "domain=org.frankensim.mini.v1",
    "domain_const=MINI_DOMAIN",
    "encoder=encode_mini",
    "encoder_helpers=none",
    "schema_functions=none",
    "schema_constants=MINI_TAG",
    "schema_dependencies=none",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=Mini",
    "source_fields=Mini.a:semantic",
    "source_bindings=Mini.a>a",
    "external_semantic_fields=none",
    "semantic_fields={semantic_fields}",
    "excluded_fields=none",
    "consumers=fixture-consumer",
    "mutations={mutations}",
    "nonsemantic_mutations=none",
    "field_guard=classify_mini",
    "transport_guard=encode_mini",
    "version_guard=crates/mini/src/lib.rs#version_refuses",
    "coupling_surface=mini:identity",
];
pub struct Mini {{
    pub a: u64,
    {extra_field}
}}
fn classify_mini(value: &Mini) {{
    let Mini {{ a: _ }} = value;
}}
fn encode_mini() {{}}
#[test]
fn mutation_a() {{ assert_eq!(1_u8, 1_u8); }}
#[test]
fn version_refuses() {{ assert_eq!(MINI_VERSION, 1); }}
"#
        )
    }

    fn seed_fixture(tag: &str, source: &str, coupling_version: u32) -> PathBuf {
        let root = fixture_root(tag);
        let source_dir = root.join("crates/mini/src");
        std::fs::create_dir_all(&source_dir).expect("fixture dirs");
        std::fs::write(source_dir.join("lib.rs"), source).expect("fixture source");
        std::fs::write(
            root.join(AUTHORITY_FILE),
            concat!(
                "{\n",
                "\"schema\": \"frankensim-identity-authorities-v1\",\n",
                "\"required_ids\": [\n",
                "{\"id\":\"mini:identity\"}\n",
                "],\n",
                "\"external_owners\": [\n",
                "],\n",
                "\"exemptions\": [\n",
                "]\n",
                "}\n",
            ),
        )
        .expect("fixture authority inventory");
        let (mut declarations, violations) = declaration_blocks("crates/mini/src/lib.rs", source);
        assert!(
            violations.is_empty(),
            "fixture declaration invalid: {violations:?}"
        );
        assert_eq!(declarations.len(), 1, "fixture must declare one identity");
        declarations[0].schema_base_hash = Some(
            identity_schema_base_hash(&root, &declarations[0], source, &[])
                .expect("fixture encoder, transport, and constants are fingerprintable"),
        );
        let resolution_violations = resolve_schema_fingerprints(&mut declarations);
        assert!(
            resolution_violations.is_empty(),
            "fixture dependencies resolve: {resolution_violations:?}"
        );
        let fingerprint = &declarations[0].schema_fingerprint;
        let coupling_fingerprint = if coupling_version == declarations[0].version {
            fingerprint.clone()
        } else {
            format!(
                "v{coupling_version}-{}",
                fingerprint
                    .split_once('-')
                    .map_or(fingerprint.as_str(), |(_, digest)| digest)
            )
        };
        std::fs::write(
            root.join("golden-couplings.json"),
            format!(
                "{{\n\"schema\": \"frankensim-golden-couplings-v1\",\n\"note\": \"fixture\",\n\"surfaces\": [\n{{\"id\": \"mini:identity\", \"file\": \"crates/mini/src/lib.rs\", \"const\": \"MINI_VERSION\", \"version\": {coupling_version}, \"domain_const\": \"MINI_DOMAIN\", \"domain\": \"org.frankensim.mini.v1\", \"schema_fingerprint\": \"{coupling_fingerprint}\"}}\n],\n\"goldens\": [\n]\n}}\n"
            ),
        )
        .expect("fixture coupling");
        root
    }

    fn identity_source(id: &str, dependencies: &str, constant_value: u8) -> String {
        owner_source("", "a", "a:crates/mini/src/lib.rs#mutation_a")
            .replace("id=mini:identity", &format!("id={id}"))
            .replace(
                "schema_dependencies=none",
                &format!("schema_dependencies={dependencies}"),
            )
            .replace(
                "pub const MINI_TAG: u8 = 1;",
                &format!("pub const MINI_TAG: u8 = {constant_value};"),
            )
    }

    fn identity_source_with_schema_functions(id: &str, functions: &str) -> String {
        identity_source(id, "none", 1).replace(
            "schema_functions=none",
            &format!("schema_functions={functions}"),
        )
    }

    fn fixture_declaration(root: &Path, owner: &str, source: &str) -> IdentityDecl {
        let owner_path = root.join(owner);
        std::fs::create_dir_all(owner_path.parent().expect("fixture owner parent"))
            .expect("fixture owner directory");
        std::fs::write(&owner_path, source).expect("fixture owner source");
        let shared_evidence = root.join("crates/mini/src/lib.rs");
        if !shared_evidence.exists() {
            std::fs::create_dir_all(shared_evidence.parent().expect("fixture evidence parent"))
                .expect("fixture evidence directory");
            std::fs::write(&shared_evidence, source).expect("fixture evidence source");
        }
        let (mut declarations, violations) = declaration_blocks(owner, source);
        assert!(
            violations.is_empty(),
            "fixture declaration invalid: {violations:?}"
        );
        assert_eq!(declarations.len(), 1, "fixture must declare one identity");
        let mut declaration = declarations.pop().expect("one fixture declaration");
        declaration.schema_base_hash = Some(
            identity_schema_base_hash(root, &declaration, source, &[])
                .expect("fixture schema base is fingerprintable"),
        );
        declaration
    }

    fn resolved_fixture(
        root: &Path,
        declarations: impl IntoIterator<Item = (&'static str, String)>,
    ) -> Vec<IdentityDecl> {
        let mut declarations = declarations
            .into_iter()
            .map(|(owner, source)| fixture_declaration(root, owner, &source))
            .collect::<Vec<_>>();
        let violations = resolve_schema_fingerprints(&mut declarations);
        assert!(
            violations.is_empty(),
            "fixture dependencies resolve: {violations:?}"
        );
        declarations.sort_by(|left, right| left.id.cmp(&right.id));
        declarations
    }

    fn write_current_registry(root: &Path) {
        let (declarations, violations) = load_declarations(root);
        assert!(violations.is_empty(), "fixture invalid: {violations:?}");
        let manifest = load_authority_manifest(root).expect("fixture authority manifest");
        let registry = render_registry(
            root,
            &declarations,
            &manifest.external_owners,
            &manifest.exemptions,
        )
        .expect("fixture registry renders");
        std::fs::write(root.join(REGISTRY_FILE), registry).expect("fixture registry");
    }

    #[test]
    fn clean_generated_registry_is_accepted() {
        let source = owner_source("", "a", "a:crates/mini/src/lib.rs#mutation_a");
        let root = seed_fixture("clean", &source, 1);
        write_current_registry(&root);
        assert!(check_identities(&root).is_empty());
    }

    #[test]
    fn rust_body_lookup_ignores_literal_and_comment_braces() {
        let source = r##"
impl Demo {
    fn earlier() {
        let _json = "}";
        let _format = "{{";
        let _raw = r#"}"#;
        let _brace = '}';
        /* } nested /* { */ } */
    }

    fn retained_encoder() {
        let _ = 7_u64;
    }
}
"##;
        let body = function_body(source, "Demo::retained_encoder")
            .expect("later method remains inside the real impl body");
        assert!(body.contains("7_u64"));
    }

    #[test]
    fn source_inventory_covers_one_line_structs_and_enum_variant_fields() {
        let source = r#"
struct OneLine { pub a: u64, pub b: Option<(u64, u64)> }
enum Observation {
    Empty,
    Count { label: String, count: u64 },
    Samples { label: String, values: Vec<u64> },
}
enum TupleObservation { Pair(u64, u64) }
fn classify(value: &OneLine, observation: &Observation) {
    let OneLine { a, b } = value;
    match observation {
        Observation::Empty => {}
        Observation::Count { label, count } => { let _ = (label, count); }
        Observation::Samples { label, values } => { let _ = (label, values); }
    }
    let _ = (a, b);
}
"#;
        assert_eq!(
            source_shape(source, "OneLine"),
            Some(SourceShape::Struct(BTreeSet::from([
                "a".to_string(),
                "b".to_string(),
            ])))
        );
        assert_eq!(
            source_shape(source, "TupleObservation"),
            Some(SourceShape::Enum {
                fields: BTreeSet::from(["variant".to_string()]),
                variants: BTreeSet::from(["Pair".to_string()]),
                tuple_variants: BTreeSet::from(["Pair".to_string()]),
            })
        );
        assert_eq!(
            source_shape(source, "Observation"),
            Some(SourceShape::Enum {
                fields: BTreeSet::from([
                    "variant".to_string(),
                    "label".to_string(),
                    "count".to_string(),
                    "values".to_string(),
                ]),
                variants: BTreeSet::from([
                    "Empty".to_string(),
                    "Count".to_string(),
                    "Samples".to_string(),
                ]),
                tuple_variants: BTreeSet::new(),
            })
        );
        assert!(guard_destructures_sources(
            source,
            "classify",
            &["OneLine".to_string(), "Observation".to_string()]
        ));
        assert!(symbol_body_has_no_rest_pattern(source, "classify"));
    }

    #[test]
    fn qualified_free_functions_and_ambiguous_impl_methods_are_exact() {
        let module_source = r#"
mod codec {
    fn encode() {
        let _ = 11_u64;
    }
}
"#;
        let module_body = function_body(module_source, "codec::encode")
            .expect("module-qualified free function resolves as a free function");
        assert!(module_body.contains("11_u64"));

        let ambiguous_impls = r#"
struct Codec;
#[cfg(feature = "left")]
impl Codec {
    fn encode() { let _ = 1_u64; }
}
#[cfg(not(feature = "left"))]
impl Codec {
    fn encode() { let _ = 2_u64; }
}
"#;
        assert_eq!(function_bodies(ambiguous_impls, "Codec::encode").len(), 2);
        assert!(function_body(ambiguous_impls, "Codec::encode").is_none());

        let trait_impl = r#"
trait Verifier { fn verify(&self); }
struct DenyAll;
impl Verifier for DenyAll {
    fn verify(&self) { let _ = "deny"; }
}
"#;
        assert!(function_body(trait_impl, "DenyAll::verify").is_some());
    }

    #[test]
    fn compact_registry_fields_drive_history_checks() {
        let baseline = concat!(
            "{\n\"schema\":\"frankensim-identity-schemas-v1\",\n\"identities\":[\n",
            "{\"id\":\"mini:identity\",\"version\":1,\"schema_fingerprint\":\"v1-deadbeefdeadbeef\"}\n",
            "]\n}\n",
        );
        assert_eq!(
            registry_fingerprints(baseline)
                .expect("strict baseline registry")
                .get(&("mini:identity".to_string(), 1))
                .and_then(Option::as_deref),
            Some("v1-deadbeefdeadbeef")
        );
        let fingerprintless_legacy = concat!(
            "{\n\"schema\":\"frankensim-identity-schemas-v1\",\n\"identities\":[\n",
            "{\"id\":\"mini:identity\",\"version\":1}\n",
            "]\n}\n",
        );
        assert_eq!(
            registry_fingerprints(fingerprintless_legacy)
                .expect("pre-fingerprint legacy registry")
                .get(&("mini:identity".to_string(), 1)),
            Some(&None)
        );
        assert!(
            schema_history_against(baseline, fingerprintless_legacy).is_empty(),
            "a legacy row without an invented fingerprint bootstraps the first fingerprinted epoch"
        );
        let current = concat!(
            "{\n\"schema\":\"frankensim-identity-schemas-v1\",\n\"identities\":[\n",
            "],\n\"external_authorities\":[\n]\n}\n",
        );
        let violations = schema_history_against(current, baseline);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("was removed")),
            "{violations:?}"
        );

        let advanced = concat!(
            "{\n\"schema\":\"frankensim-identity-schemas-v1\",\n\"identities\":[\n",
            "{\"id\":\"mini:identity\",\"version\":2,\"schema_fingerprint\":\"v2-feedfacefeedface\"}\n",
            "]\n}\n",
        );
        let violations = schema_history_against(advanced, baseline);
        assert!(
            violations.is_empty(),
            "a higher version is a deliberate schema migration: {violations:?}"
        );

        let drifted = baseline.replace("deadbeefdeadbeef", "feedfacefeedface");
        let violations = schema_history_against(&drifted, baseline);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("changed schema fingerprint at retained version 1")),
            "{violations:?}"
        );

        let regressed = baseline
            .replace("\"version\":1", "\"version\":0")
            .replace("v1-deadbeefdeadbeef", "v0-deadbeefdeadbeef");
        let violations = schema_history_against(&regressed, baseline);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("version regressed")),
            "{violations:?}"
        );

        let empty_baseline = concat!(
            "{\n\"schema\":\"frankensim-identity-schemas-v1\",\n\"identities\":[\n",
            "],\n\"external_authorities\":[\n]\n}\n",
        );
        assert!(
            schema_history_against(baseline, empty_baseline).is_empty(),
            "the first retained registry epoch must bootstrap from empty history"
        );

        let duplicate = baseline.replace(
            "\"id\":\"mini:identity\"",
            "\"id\":\"mini:identity\",\"id\":\"hidden:identity\"",
        );
        let violations = schema_history_against(current, &duplicate);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("retained identity registry is malformed")),
            "{violations:?}"
        );
    }

    #[test]
    fn authority_and_coupling_manifests_reject_noncanonical_json_shapes() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "frankensim-identity-strict-json-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir_all(&root).expect("strict JSON fixture root");
        std::fs::write(
            root.join(AUTHORITY_FILE),
            concat!(
                "{\n",
                "\"schema\":\"frankensim-identity-authorities-v1\",\n",
                "\"schema\":\"shadow\",\n",
                "\"required_ids\":[],\n\"external_owners\":[],\n\"exemptions\":[]\n",
                "}\n",
            ),
        )
        .expect("duplicate authority manifest");
        let authority = load_authority_manifest(&root).expect_err("duplicate key must fail");
        assert!(
            authority
                .iter()
                .any(|violation| violation.detail.contains("not strict JSON")),
            "{authority:?}"
        );

        std::fs::write(
            root.join("golden-couplings.json"),
            concat!(
                "{\n",
                "\"schema\":\"frankensim-golden-couplings-v1\",\n",
                "\"note\":\"fixture\",\n",
                "\"surfaces\":[],\n\"goldens\":[],\n",
                "\"shadow\":[]\n",
                "}\n",
            ),
        )
        .expect("unknown coupling key");
        let coupling = coupling_surfaces(&root).expect_err("unknown key must fail");
        assert!(
            coupling
                .iter()
                .any(|violation| violation.detail.contains("noncanonical keys")),
            "{coupling:?}"
        );
    }

    #[test]
    fn top_level_script_external_authority_is_rendered_and_history_guarded() {
        let root = fixture_root("top-level-script-authority");
        let path = root.join("scripts/proof.py");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let baseline_source = concat!(
            "import hashlib\n",
            "import pathlib\n",
            "DOMAIN = \"org.frankensim.ci.top-level-proof.v1\"\n",
            "digest = hashlib.sha256(f\"{DOMAIN}:payload-v1\".encode()).hexdigest()\n",
            "pathlib.Path(\"proof.receipt\").write_text(digest)\n",
        );
        std::fs::write(&path, baseline_source).expect("baseline top-level script");
        let owner = ExternalOwner {
            id: "ci:top-level-proof".to_string(),
            path: "scripts/proof.py".to_string(),
            symbol: "<script>".to_string(),
            version: 1,
            domain: "org.frankensim.ci.top-level-proof.v1".to_string(),
        };
        let candidates = discover_identity_candidates(&root);
        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from([owner.id.clone()]),
            external_owners: vec![owner.clone()],
            exemptions: Vec::new(),
        };
        let violations = authority_violations_against(&root, &[], &manifest, &candidates);
        assert!(violations.is_empty(), "{violations:?}");
        let baseline = render_registry(&root, &[], std::slice::from_ref(&owner), &[])
            .expect("top-level external authority renders");
        assert!(baseline.contains("\"symbol\":\"<script>\""));

        std::fs::write(&path, baseline_source.replace("payload-v1", "payload-v2"))
            .expect("moved top-level script");
        let moved =
            render_registry(&root, &[], &[owner], &[]).expect("moved top-level authority renders");
        let violations = schema_history_against(&moved, &baseline);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("changed schema fingerprint at retained version 1")),
            "{violations:?}"
        );
    }

    #[test]
    fn shell_python_blocks_are_not_rediscovered_as_top_level_producers() {
        let script = concat!(
            "python3 - <<'PY'\n",
            "def emit_receipt():\n",
            "    import hashlib\n",
            "    import pathlib\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.heredoc-proof.v1'\n",
            "    digest = hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest()\n",
            "    pathlib.Path('proof.receipt').write_text(digest)\n",
            "emit_receipt()\n",
            "PY\n",
            "printf '%s\\n' unrelated >> proof.log\n",
        );
        let candidates = discover_script_candidates("scripts/ci/heredoc.sh", script);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["emit_receipt"],
            "a Python heredoc producer must not combine with an unrelated shell sink into a synthetic top-level producer: {candidates:?}"
        );
    }

    #[test]
    fn alternate_script_function_syntax_keeps_exact_ownership() {
        let shell = concat!(
            "producer ()\n",
            "{\n",
            "    DOMAIN=org.frankensim.ci.shell-proof.v1\n",
            "    printf '%s\\n' \"$DOMAIN\" | shasum -a 256 >> proof.receipt\n",
            "}\n",
        );
        let shell_candidates = discover_script_candidates("scripts/ci/proof.sh", shell);
        assert_eq!(
            shell_candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"],
            "spaced and multiline shell function headers must not fall into <script>: {shell_candidates:?}"
        );

        let python = concat!(
            "import hashlib\n",
            "import pathlib\n",
            "async def producer():\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.async-proof.v1'\n",
            "    digest = hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest()\n",
            "    pathlib.Path('proof.receipt').write_text(digest)\n",
        );
        let python_candidates = discover_script_candidates("scripts/ci/proof.py", python);
        assert_eq!(
            python_candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"],
            "async Python producers must retain their exact function authority: {python_candidates:?}"
        );
    }

    #[test]
    fn python_signatures_balance_nested_defaults_multiline_forms_and_return_annotations() {
        let source = concat!(
            "import hashlib\n",
            "import pathlib\n",
            "def annotated() -> None:\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.annotated.v1'\n",
            "    pathlib.Path('annotated.receipt').write_text(hashlib.sha256(b'a').hexdigest())\n",
            "def nested(value=(1, (2, 3))):\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.nested.v1'\n",
            "    pathlib.Path('nested.receipt').write_text(hashlib.sha256(b'b').hexdigest())\n",
            "def multiline(\n",
            "    value: tuple[int, int] = (1, 2),\n",
            ") -> None:\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.multiline.v1'\n",
            "    pathlib.Path('multiline.receipt').write_text(hashlib.sha256(b'c').hexdigest())\n",
        );
        let symbols = discover_script_candidates("scripts/ci/signatures.py", source)
            .into_iter()
            .map(|candidate| candidate.symbol)
            .collect::<BTreeSet<_>>();
        assert_eq!(
            symbols,
            BTreeSet::from([
                "annotated".to_string(),
                "multiline".to_string(),
                "nested".to_string(),
            ])
        );
    }

    #[test]
    fn python_f_string_replacements_are_executable_and_malformed_fields_fail_closed() {
        let executable = concat!(
            "def producer():\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.f-string.v1'\n",
            "    payload = f\"{pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())}\"\n",
        );
        let candidates = discover_script_candidates("scripts/ci/f_string.py", executable);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"],
            "an identity sink inside a replacement expression must remain executable: {candidates:?}"
        );

        let malformed = concat!(
            "IDENTITY_DOMAIN = 'org.frankensim.ci.f-string.v1'\n",
            "pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
            "payload = f\"{missing\"\n",
        );
        let error = super::discover_script_candidates("scripts/ci/f_string.py", malformed)
            .expect_err("an unclosed f-string replacement must fail closed");
        assert!(matches!(
            error.kind,
            ScriptStructureKind::IncompleteFStringReplacement
        ));
        assert_eq!(error.opening_offset, malformed.find("{missing").unwrap());
    }

    #[test]
    fn embedded_python_blocks_retain_dialect_and_outer_source_offsets() {
        let source = concat!(
            "outer() {\n",
            "  python3 <<'PY'\n",
            "@audit\n",
            "def embedded():\n",
            "    payload = 'can\\'t'\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.embedded.v1'\n",
            "    pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
            "PY\n",
            "}\n",
        );
        let blocks = super::python_function_blocks("scripts/ci/embedded.sh", source)
            .expect("the executable heredoc is valid Python");
        let [block] = blocks.as_slice() else {
            panic!("embedded Python function must resolve exactly once: {blocks:?}");
        };
        assert_eq!(block.dialect, ScriptDialect::Python);
        assert_eq!(block.start, source.find("@audit").unwrap());
        assert!(block.text.starts_with("@audit\ndef embedded():"));
        assert_eq!(
            discover_script_candidates("scripts/ci/embedded.sh", source)
                .into_iter()
                .map(|candidate| candidate.symbol)
                .collect::<Vec<_>>(),
            vec!["embedded".to_string()],
            "a Python apostrophe escape must not be interpreted with shell quote rules"
        );

        let malformed = source.replace("'can\\'t'", "'unterminated");
        let error = super::discover_script_candidates("scripts/ci/embedded.sh", &malformed)
            .expect_err("an embedded Python error must fail at its outer-source offset");
        assert_eq!(
            error.opening_offset,
            malformed.find("'unterminated").unwrap(),
            "the heredoc body offset must be shifted into the governing shell source"
        );
        assert!(matches!(
            error.kind,
            ScriptStructureKind::UnterminatedQuote {
                dialect: ScriptDialect::Python,
                delimiter: '\'',
            }
        ));
    }

    #[test]
    fn python_fingerprint_slices_include_decorators_but_exclude_dedented_comments() {
        let root = fixture_root("python-decorator-fingerprint-slice");
        let path = root.join("scripts/proof.py");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let baseline_source = concat!(
            "@audit('v1')\n",
            "@deterministic\n",
            "def producer() -> None:\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.decorated.v1'\n",
            "    pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
            "# unrelated trailing comment v1\n",
            "unrelated = 1\n",
        );
        std::fs::write(&path, baseline_source).expect("baseline decorated producer");
        let owner = ExternalOwner {
            id: "ci:decorated".to_string(),
            path: "scripts/proof.py".to_string(),
            symbol: "producer".to_string(),
            version: 1,
            domain: "org.frankensim.ci.decorated.v1".to_string(),
        };
        let blocks = super::primary_script_blocks(&owner.path, baseline_source, &owner.symbol)
            .expect("decorated producer resolves");
        let [block] = blocks.as_slice() else {
            panic!("decorated producer must resolve exactly once: {blocks:?}");
        };
        assert!(block.text.starts_with("@audit('v1')\n@deterministic\n"));
        assert!(!block.text.contains("unrelated trailing comment"));

        let baseline =
            external_owner_schema_fingerprint(&root, &owner, &[]).expect("baseline fingerprint");
        std::fs::write(
            &path,
            baseline_source.replace("@audit('v1')", "@audit('v2')"),
        )
        .expect("moved decorator");
        let moved = external_owner_schema_fingerprint(&root, &owner, &[])
            .expect("moved decorator fingerprint");
        assert_ne!(
            baseline, moved,
            "decorators are part of the producer schema"
        );

        std::fs::write(
            &path,
            baseline_source.replace("trailing comment v1", "trailing comment v2"),
        )
        .expect("moved dedented comment");
        let trailing = external_owner_schema_fingerprint(&root, &owner, &[])
            .expect("trailing-comment fingerprint");
        assert_eq!(
            baseline, trailing,
            "the first dedented comment is outside the exact producer slice"
        );
    }

    #[test]
    fn python_crlf_escape_continuations_and_utf8_columns_are_exact() {
        let continued = "payload = \"left\\\r\nright\"\r\n";
        let (code, source) = super::script_views(ScriptDialect::Python, continued)
            .expect("a backslash-CRLF pair continues a Python string literal");
        assert_eq!(code.len(), continued.len());
        assert_eq!(source, continued);

        let malformed = "prefix éé payload = 'unterminated";
        let opening = malformed.find("'unterminated").unwrap();
        let error = super::script_views(ScriptDialect::Python, malformed)
            .expect_err("the malformed UTF-8 fixture must fail closed");
        let expected_column = malformed[..opening].chars().count() + 1;
        let detail = script_structure_detail("scripts/ci/utf8.py", malformed, &error);
        assert!(
            detail.contains(&format!("line 1, column {expected_column}:")),
            "diagnostic columns must count Unicode scalar values: {detail}"
        );
    }

    #[test]
    fn comments_and_inert_literals_are_not_external_identity_producers() {
        let root = fixture_root("external-authority-decoys");
        let script_path = root.join("scripts/decoy.py");
        let rust_path = root.join("crates/decoy/src/lib.rs");
        std::fs::create_dir_all(script_path.parent().expect("script parent"))
            .expect("script directory");
        std::fs::create_dir_all(rust_path.parent().expect("Rust parent")).expect("Rust directory");
        std::fs::write(
            &script_path,
            "# org.frankensim.ci.decoy.v1 sha256 >> receipt\n",
        )
        .expect("script comment decoy");
        std::fs::write(
            &rust_path,
            concat!(
                "fn decoy() {\n",
                "    let _ = \"org.frankensim.rust.decoy.v1 hash_domain std::fs::write\";\n",
                "}\n",
            ),
        )
        .expect("Rust literal decoy");
        let candidates = discover_identity_candidates(&root);
        assert!(
            candidates.is_empty(),
            "inert tokens are not producers: {candidates:?}"
        );

        for owner in [
            ExternalOwner {
                id: "ci:decoy".to_string(),
                path: "scripts/decoy.py".to_string(),
                symbol: "<script>".to_string(),
                version: 1,
                domain: "org.frankensim.ci.decoy.v1".to_string(),
            },
            ExternalOwner {
                id: "rust:decoy".to_string(),
                path: "crates/decoy/src/lib.rs".to_string(),
                symbol: "decoy".to_string(),
                version: 1,
                domain: "org.frankensim.rust.decoy.v1".to_string(),
            },
        ] {
            let manifest = AuthorityManifest {
                required_ids: BTreeSet::from([owner.id.clone()]),
                external_owners: vec![owner.clone()],
                exemptions: Vec::new(),
            };
            let violations = authority_violations_against(&root, &[], &manifest, &candidates);
            assert!(
                violations.iter().any(|violation| violation
                    .detail
                    .contains("matches 0 independently discovered producers")),
                "an inert external target must fail independent producer coverage: {violations:?}"
            );
        }
    }

    #[test]
    fn quoted_function_decoys_and_data_heredocs_are_not_producers() {
        let python = concat!(
            "payload = \"\"\"\n",
            "def forged():\n",
            "    digest = hashlib.sha256(b'identity-domain').hexdigest()\n",
            "    pathlib.Path('receipt').write_text(digest)\n",
            "\"\"\"\n",
        );
        assert!(
            discover_script_candidates("scripts/decoy.py", python).is_empty(),
            "a def-shaped string is not a Python producer"
        );

        let shell_quote = concat!(
            "payload='\n",
            "forged() { shasum -a 256 >> receipt; }\n",
            "'\n",
        );
        assert!(
            discover_script_candidates("scripts/decoy.sh", shell_quote).is_empty(),
            "a function-shaped multiline shell string is not a producer"
        );

        let data_heredoc = concat!(
            "cat <<'EOF'\n",
            "def forged():\n",
            "    digest = hashlib.sha256(b'identity-domain').hexdigest()\n",
            "    print(digest)\n",
            "EOF\n",
            "printf '%s\\n' unrelated >> proof.log\n",
        );
        assert!(
            discover_script_candidates("scripts/decoy.sh", data_heredoc).is_empty(),
            "data heredoc tokens must not combine with an unrelated executable sink"
        );
    }

    #[test]
    fn identity_bearing_malformed_scripts_fail_closed() {
        let shell_quote = super::discover_script_candidates(
            "scripts/proof.sh",
            concat!(
                "producer() {\n",
                "  DOMAIN=org.frankensim.ci.proof.v1\n",
                "  printf '%s\\n' \"$DOMAIN\" >> proof.receipt\n",
                "  payload='unterminated\n",
            ),
        )
        .expect_err("unterminated shell quote must fail closed");
        assert!(matches!(
            shell_quote.kind,
            ScriptStructureKind::UnterminatedQuote {
                dialect: ScriptDialect::Shell,
                delimiter: '\'',
            }
        ));

        let python_quote = super::discover_script_candidates(
            "scripts/proof.py",
            concat!(
                "IDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\n",
                "pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
                "payload = \"unterminated",
            ),
        )
        .expect_err("unterminated Python quote must fail closed");
        assert!(matches!(
            python_quote.kind,
            ScriptStructureKind::UnterminatedQuote {
                dialect: ScriptDialect::Python,
                delimiter: '"',
            }
        ));

        let python_triple = super::discover_script_candidates(
            "scripts/proof.py",
            concat!(
                "IDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\n",
                "pathlib.Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
                "payload = \"\"\"unterminated",
            ),
        )
        .expect_err("unterminated Python triple string must fail closed");
        assert!(matches!(
            python_triple.kind,
            ScriptStructureKind::UnterminatedTripleString {
                delimiter: "\"\"\""
            }
        ));

        let heredoc = super::discover_script_candidates(
            "scripts/proof.sh",
            concat!(
                "producer() { python3 <<'PY'\n",
                "import hashlib\n",
                "from pathlib import Path\n",
                "IDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\n",
                "Path('proof.receipt').write_text(hashlib.sha256(b'x').hexdigest())\n",
            ),
        )
        .expect_err("unterminated executable heredoc must fail closed");
        assert!(matches!(
            heredoc.kind,
            ScriptStructureKind::UnterminatedHeredoc { ref delimiter }
                if delimiter == "PY"
        ));

        let shell_scope = super::discover_script_candidates(
            "scripts/proof.sh",
            concat!(
                "producer() {\n",
                "  DOMAIN=org.frankensim.ci.proof.v1\n",
                "  printf '%s\\n' \"$DOMAIN\" >> proof.receipt\n",
            ),
        )
        .expect_err("unclosed shell function must fail closed");
        assert!(matches!(
            shell_scope.kind,
            ScriptStructureKind::IncompleteFunctionScope {
                dialect: ScriptDialect::Shell,
                ref symbol,
                expected: "a closing '}'",
            } if symbol.as_deref() == Some("producer")
        ));

        for source in [
            "def producer(\n",
            "def producer()\n",
            "def producer():\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\n",
        ] {
            let error = super::discover_script_candidates("scripts/proof.py", source)
                .expect_err("incomplete Python function must fail closed");
            assert!(matches!(
                error.kind,
                ScriptStructureKind::IncompleteFunctionScope {
                    dialect: ScriptDialect::Python,
                    ..
                }
            ));
        }
    }

    #[test]
    fn malformed_manifest_target_is_relevant_but_unrelated_script_is_not_linted() {
        let root = fixture_root("malformed-script-relevance");
        let path = root.join("scripts/hidden.py");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        std::fs::write(&path, "safe = 1\npayload = 'unterminated")
            .expect("malformed script fixture");

        let unrelated = super::discover_identity_candidates(&root, &empty_authority_manifest())
            .expect("unrelated malformed scripts are outside this policy gate");
        assert!(unrelated.is_empty(), "{unrelated:?}");

        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from(["ci:hidden".to_string()]),
            external_owners: vec![ExternalOwner {
                id: "ci:hidden".to_string(),
                path: "scripts/hidden.py".to_string(),
                symbol: "<script>".to_string(),
                version: 1,
                domain: "org.frankensim.ci.hidden.v1".to_string(),
            }],
            exemptions: Vec::new(),
        };
        let violations = super::discover_identity_candidates(&root, &manifest)
            .expect_err("a registered malformed target must fail closed");
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert!(
            violations[0]
                .detail
                .contains("scripts/hidden.py at line 2, column 11: unterminated Python quote"),
            "{violations:?}"
        );
    }

    #[test]
    fn repository_wide_script_discovery_includes_supported_shebangs_and_excludes_generated_trees() {
        let root = fixture_root("repository-wide-script-discovery");
        let write = |relative: &str, text: &str| {
            let path = root.join(relative);
            std::fs::create_dir_all(path.parent().expect("script parent"))
                .expect("script directory");
            std::fs::write(path, text).expect("script fixture");
        };
        let shell = concat!(
            "#!/usr/bin/env bash\n",
            "emit_shell() {\n",
            "  IDENTITY_DOMAIN=org.frankensim.ci.shell-proof.v1\n",
            "  printf '%s\\n' \"$IDENTITY_DOMAIN\" | shasum -a 256 >> proof.receipt\n",
            "}\n",
        );
        let python = concat!(
            "import hashlib\n",
            "import pathlib\n",
            "def emit_python():\n",
            "    IDENTITY_DOMAIN = 'org.frankensim.ci.python-proof.v1'\n",
            "    pathlib.Path('proof.receipt').write_text(hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest())\n",
        );
        write("tools/proof", shell);
        write("ci/proof.py", python);
        write(
            "root_proof",
            &format!("#!/usr/bin/env -S python3 -u\n{python}"),
        );
        write("target/hidden.sh", shell);
        write("vendor/nested/.git/HEAD", "ref: refs/heads/main\n");
        write("vendor/nested/hidden.py", python);

        let candidates = super::discover_identity_candidates(&root, &empty_authority_manifest())
            .expect("controlled-tree scripts are readable and structurally valid");
        let actual = candidates
            .iter()
            .map(|candidate| (candidate.path.as_str(), candidate.symbol.as_str()))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            actual,
            BTreeSet::from([
                ("ci/proof.py", "emit_python"),
                ("root_proof", "emit_python"),
                ("tools/proof", "emit_shell"),
            ]),
            "scripts outside scripts/, including supported extensionless shebangs, must be independently owned: {candidates:?}"
        );
    }

    #[test]
    fn relevant_script_io_and_utf8_failures_are_violations() {
        let root = fixture_root("script-io-failures");
        let path = root.join("ci/bad.sh");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let mut bytes = b"#!/bin/sh\nIDENTITY_DOMAIN=org.frankensim.ci.bad.v1\n".to_vec();
        bytes.push(0xff);
        std::fs::write(&path, bytes).expect("invalid UTF-8 fixture");
        let violations = super::discover_identity_candidates(&root, &empty_authority_manifest())
            .expect_err("a known script source with invalid UTF-8 must fail closed");
        assert!(
            violations
                .iter()
                .any(|violation| violation.crate_name == "ci/bad.sh"
                    && violation.detail.contains("not valid UTF-8")),
            "{violations:?}"
        );

        let missing = fixture_root("missing-script-tree").join("absent");
        let violations = super::discover_identity_candidates(&missing, &empty_authority_manifest())
            .expect_err("an unreadable controlled root must not become an empty inventory");
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("cannot read directory")),
            "{violations:?}"
        );
    }

    #[cfg(unix)]
    #[test]
    fn rust_source_discovery_rejects_symlinked_sources() {
        use std::os::unix::fs::symlink;

        let root = fixture_root("symlinked-rust-source");
        let source_dir = root.join("crates/mini/src");
        std::fs::create_dir_all(&source_dir).expect("source fixture directory");
        std::fs::write(source_dir.join("real.rs"), "pub const VALUE: u8 = 1;\n")
            .expect("source fixture target");
        symlink("real.rs", source_dir.join("lib.rs")).expect("source fixture symlink");

        let (_, violations) = source_files(&root);
        assert!(
            violations.iter().any(|violation| {
                violation.crate_name == "crates/mini/src/lib.rs"
                    && violation.detail.contains("refuses symlinked paths")
            }),
            "{violations:?}"
        );
    }

    #[test]
    fn rust_source_read_failures_are_violations() {
        let root = fixture_root("invalid-rust-source-utf8");
        let path = root.join("crates/mini/src/lib.rs");
        std::fs::create_dir_all(path.parent().expect("source fixture parent"))
            .expect("source fixture directory");
        let mut bytes = b"pub const IDENTITY_DOMAIN: &str = \"org.frankensim.bad.v1\";\n".to_vec();
        bytes.push(0xff);
        std::fs::write(&path, bytes).expect("invalid UTF-8 source fixture");

        let violations = super::discover_identity_candidates(&root, &empty_authority_manifest())
            .expect_err("an unreadable Rust source must fail closed");
        assert!(
            violations.iter().any(|violation| {
                violation.crate_name == "crates/mini/src/lib.rs"
                    && violation.detail.contains("cannot read Rust source")
            }),
            "{violations:?}"
        );
    }

    #[test]
    fn executable_python_heredoc_uses_python_string_boundaries() {
        let script = concat!(
            "row() {\n",
            "  python3 <<'PY'\n",
            "payload = \"\"\"unmatched shell ' and } stay Python data\"\"\"\n",
            "import json\n",
            "print(json.dumps({\n",
            "    \"identity_domain\": \"org.frankensim.ci.row.v1\",\n",
            "    \"identity_version\": 1,\n",
            "}, separators=(\",\", \":\")))\n",
            "PY\n",
            "}\n",
        );
        let candidates = discover_script_candidates("scripts/row.sh", script);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["row"]
        );

        let malformed = script.replacen("stay Python data\"\"\"", "stay Python data", 1);
        let error = super::discover_script_candidates("scripts/row.sh", &malformed)
            .expect_err("an unterminated embedded Python triple string must fail closed");
        assert!(matches!(
            error.kind,
            ScriptStructureKind::UnterminatedTripleString {
                delimiter: "\"\"\""
            }
        ));
    }

    #[test]
    fn pipeline_and_interpreter_wrappers_classify_executable_heredocs() {
        for opener in [
            "cat <<'PY' | python3",
            "env -i python3 - <<'PY'",
            "command -p python3 - <<'PY'",
            "exec python3 - <<'PY'",
        ] {
            let script = format!(
                "row() {{\n  {opener}\nimport hashlib\nimport pathlib\nIDENTITY_DOMAIN = 'org.frankensim.ci.pipeline.v1'\npathlib.Path('proof.receipt').write_text(hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest())\nPY\n}}\n"
            );
            let candidates = discover_script_candidates("scripts/pipeline.sh", &script);
            assert_eq!(
                candidates
                    .iter()
                    .map(|candidate| candidate.symbol.as_str())
                    .collect::<Vec<_>>(),
                vec!["row"],
                "the governing interpreter must classify {opener:?}: {candidates:?}"
            );
        }
    }

    #[test]
    fn identity_bearing_opaque_heredoc_has_typed_ambiguity() {
        let script = concat!(
            "cat <<'DATA'\n",
            "import hashlib\n",
            "import pathlib\n",
            "IDENTITY_DOMAIN = 'org.frankensim.ci.ambiguous.v1'\n",
            "pathlib.Path('proof.receipt').write_text(hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest())\n",
            "DATA\n",
        );
        let error = super::discover_script_candidates("scripts/ambiguous.sh", script)
            .expect_err("identity-bearing opaque heredocs must fail closed");
        assert!(matches!(
            &error.kind,
            ScriptStructureKind::AmbiguousHeredocExecution { delimiter }
                if delimiter == "DATA"
        ));
        assert!(
            script_structure_detail("scripts/ambiguous.sh", script, &error)
                .contains("has no unambiguous shell or Python execution context")
        );
    }

    #[test]
    fn heredoc_activation_waits_for_command_completion_and_same_line_scope_close() {
        let continued = concat!(
            "producer() {\n",
            "  python3 - <<'PY' \\\n",
            "    --sentinel\n",
            "import hashlib\n",
            "import pathlib\n",
            "IDENTITY_DOMAIN = 'org.frankensim.ci.continued.v1'\n",
            "pathlib.Path('proof.receipt').write_text(hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest())\n",
            "PY\n",
            "}\n",
        );
        assert_eq!(
            discover_script_candidates("scripts/continued.sh", continued)
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"]
        );

        let same_line_close = concat!(
            "producer() { python3 - <<'PY'; }\n",
            "import hashlib\n",
            "import pathlib\n",
            "IDENTITY_DOMAIN = 'org.frankensim.ci.same-line.v1'\n",
            "pathlib.Path('proof.receipt').write_text(hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest())\n",
            "PY\n",
        );
        assert_eq!(
            discover_script_candidates("scripts/same-line.sh", same_line_close)
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"]
        );
        let blocks = primary_script_blocks("scripts/same-line.sh", same_line_close, "producer");
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].ends_with("PY\n"));
    }

    #[test]
    fn shell_shift_expressions_and_delimiter_words_are_not_confused() {
        assert!(shell_heredocs("value=$((1 << 2))").is_empty());
        assert!(shell_heredocs("done <<<\"$entries\"").is_empty());
        assert_eq!(
            shell_heredocs("cat <<\\EOF")
                .into_iter()
                .map(|heredoc| heredoc.delimiter)
                .collect::<Vec<_>>(),
            vec!["EOF"]
        );
        assert_eq!(
            shell_heredocs("cat <<'END.JSON'")
                .into_iter()
                .map(|heredoc| heredoc.delimiter)
                .collect::<Vec<_>>(),
            vec!["END.JSON"]
        );
        for (opener, closer) in [("<<\\EOF", "EOF"), ("<<'END.JSON'", "END.JSON")] {
            let script = format!("cat {opener}\nproducer() {{\n{closer}\n");
            let (code, _) = executable_script_views("scripts/data.sh", &script);
            assert!(
                !code.contains("producer"),
                "the body for delimiter {closer:?} must remain opaque"
            );
        }
    }

    #[test]
    fn executable_shell_heredoc_scopes_are_discovered_and_validated_recursively() {
        let valid = concat!(
            "bash <<'SH'\n",
            "producer() {\n",
            "  DOMAIN=org.frankensim.ci.nested-shell.v1\n",
            "  printf '%s\\n' \"$DOMAIN\" | shasum -a 256 >> proof.receipt\n",
            "}\n",
            "SH\n",
        );
        assert_eq!(
            discover_script_candidates("scripts/nested.sh", valid)
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["producer"]
        );

        let malformed = concat!("bash <<'SH'\n", "producer() {\n", "SH\n");
        let error = super::discover_script_candidates("scripts/nested.sh", malformed)
            .expect_err("an executable shell heredoc must validate its own scopes");
        assert!(matches!(
            error.kind,
            ScriptStructureKind::IncompleteFunctionScope {
                dialect: ScriptDialect::Shell,
                ref symbol,
                expected: "a closing '}'",
            } if symbol.as_deref() == Some("producer")
        ));
    }

    #[test]
    fn structured_json_domain_fields_are_executable_identity_signals() {
        let script = concat!(
            "row() {\n",
            "  python3 - <<'PY'\n",
            "import json\n",
            "print(json.dumps({\n",
            "    \"identity_domain\": \"org.frankensim.ci.row.v1\",\n",
            "    \"identity_version\": 1,\n",
            "    \"status\": \"pass\",\n",
            "}, separators=(\",\", \":\")))\n",
            "PY\n",
            "}\n",
        );
        let candidates = discover_script_candidates("scripts/row.sh", script);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| (&*candidate.symbol, &*candidate.identity_signal))
                .collect::<Vec<_>>(),
            vec![("row", "structured-json-identity")],
            "an executable canonical JSON identity row remains discoverable: {candidates:?}"
        );
    }

    #[test]
    fn continued_python_heredocs_keep_their_governing_interpreter() {
        let quality = concat!(
            "emit_proof_seal() {\n",
            "  if ! seal_row=$(python3 - \"$HEAD\" \\\n",
            "      \"$STATUS\" <<'PY'\n",
            "import json\n",
            "print(json.dumps({\n",
            "    \"identity_domain\": \"org.frankensim.ci.quality-proof-record.v1\",\n",
            "    \"identity_version\": 1,\n",
            "}, separators=(\",\", \":\")))\n",
            "PY\n",
            "  ); then return 1; fi\n",
            "  printf '%s\\n' \"$seal_row\" | tee -a \"$VERDICTS\"\n",
            "}\n",
        );
        let cross = concat!(
            "row() {\n",
            "  python3 - \"$STATUS\" \\\n",
            "    \"$DETAIL\" <<'PY' |\n",
            "import json\n",
            "print(json.dumps({\n",
            "    \"identity_domain\": \"org.frankensim.ci.x86-cross-verdict.v1\",\n",
            "    \"identity_version\": 1,\n",
            "}, separators=(\",\", \":\")))\n",
            "PY\n",
            "    tee -a \"$VERDICTS\"\n",
            "}\n",
        );
        for (path, source, symbol) in [
            ("scripts/ci/quality_lanes.sh", quality, "emit_proof_seal"),
            ("scripts/ci/x86_cross_check.sh", cross, "row"),
        ] {
            let candidates = discover_script_candidates(path, source);
            assert_eq!(
                candidates
                    .iter()
                    .map(|candidate| (&*candidate.symbol, &*candidate.identity_signal))
                    .collect::<Vec<_>>(),
                vec![(symbol, "structured-json-identity")],
                "continued interpreter command must retain its executable heredoc: {candidates:?}"
            );
        }
    }

    #[test]
    fn interpreter_decoys_do_not_silently_activate_data_heredocs() {
        for opener in [
            "python3 --version; cat <<'PY'",
            "python3 --version && cat <<'PY'",
            "cat python3 <<'PY'",
            "PYTHON=python3 cat <<'PY'",
            "note='python3'; cat <<'PY'",
            "cat <<'PY' # python3",
            "stamp=$(python3 --version); cat <<'PY'",
        ] {
            let script = format!(
                "row() {{\n  {opener}\nimport json\nprint(json.dumps({{\n    \"identity_domain\": \"org.frankensim.ci.decoy.v1\",\n    \"identity_version\": 1,\n}}, separators=(\",\", \":\")))\nPY\n  printf '%s\\n' unrelated >> proof.log\n}}\n"
            );
            let error = super::discover_script_candidates("scripts/ci/decoy.sh", &script)
                .expect_err("identity-bearing data heredocs require an execution context");
            assert!(
                matches!(
                    &error.kind,
                    ScriptStructureKind::AmbiguousHeredocExecution { delimiter }
                        if delimiter == "PY"
                ),
                "an interpreter outside the governing simple command must not activate a data heredoc ({opener:?}): {error:?}"
            );
        }
    }

    #[test]
    fn live_external_script_and_bootstrap_producers_remain_discoverable() {
        for (path, source, symbol, domain) in [
            (
                "scripts/ci/quality_lanes.sh",
                include_str!("../../scripts/ci/quality_lanes.sh"),
                "emit_proof_seal",
                "org.frankensim.ci.quality-proof-record.v1",
            ),
            (
                "scripts/ci/x86_cross_check.sh",
                include_str!("../../scripts/ci/x86_cross_check.sh"),
                "row",
                "org.frankensim.ci.x86-cross-verdict.v1",
            ),
        ] {
            let candidates = discover_script_candidates(path, source)
                .into_iter()
                .filter(|candidate| candidate.symbol == symbol)
                .collect::<Vec<_>>();
            assert_eq!(candidates.len(), 1, "{path}#{symbol}: {candidates:?}");
            let block = primary_script_blocks(path, source, symbol);
            let [block] = block.as_slice() else {
                panic!("{path}#{symbol} must resolve exactly once");
            };
            let (code, uncommented) = executable_script_views(path, block);
            assert_eq!(
                script_declared_domains(&code, &uncommented),
                vec![domain.to_string()],
                "{path}#{symbol} must bind its exact domain"
            );
        }

        let bootstrap = discover_rust_candidates("xtask/src/main.rs", include_str!("main.rs"))
            .into_iter()
            .filter(|candidate| candidate.symbol == "cmd_bootstrap")
            .collect::<Vec<_>>();
        assert_eq!(bootstrap.len(), 1, "{bootstrap:?}");
        assert_eq!(bootstrap[0].identity_signal, "identity_domain");
        assert_eq!(bootstrap[0].sink_signal, "std::fs::write");
    }

    #[test]
    fn top_level_external_owner_refuses_sources_with_function_scopes() {
        let root = fixture_root("top-level-owner-function-scope");
        let path = root.join("scripts/proof.py");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        std::fs::write(
            &path,
            concat!(
                "async def producer():\n",
                "    import hashlib\n",
                "    import pathlib\n",
                "    digest = hashlib.sha256(b'org.frankensim.ci.proof.v1').hexdigest()\n",
                "    pathlib.Path('proof.receipt').write_text(digest)\n",
            ),
        )
        .expect("function-scoped producer");
        let owner = ExternalOwner {
            id: "ci:proof".to_string(),
            path: "scripts/proof.py".to_string(),
            symbol: "<script>".to_string(),
            version: 1,
            domain: "org.frankensim.ci.proof.v1".to_string(),
        };
        let candidates = discover_identity_candidates(&root);
        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from([owner.id.clone()]),
            external_owners: vec![owner],
            exemptions: Vec::new(),
        };
        let violations = authority_violations_against(&root, &[], &manifest, &candidates);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("<script> cannot own a source containing function scopes")),
            "a top-level row may not absorb a function-scoped producer: {violations:?}"
        );
    }

    #[test]
    fn external_authority_schema_is_generated_and_history_guarded() {
        let root = fixture_root("external-registry-history");
        let path = root.join("scripts/ci/proof.sh");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let owner = ExternalOwner {
            id: "ci:proof".to_string(),
            path: "scripts/ci/proof.sh".to_string(),
            symbol: "row".to_string(),
            version: 1,
            domain: "org.frankensim.ci.proof.v1".to_string(),
        };
        std::fs::write(
            &path,
            "row() { python3 - <<'PY'\nimport hashlib\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\nprint(hashlib.sha256(f'{IDENTITY_DOMAIN}:identity-v1'.encode()).hexdigest())\nPY\n}\n",
        )
        .expect("baseline script");
        let baseline = render_registry(&root, &[], std::slice::from_ref(&owner), &[])
            .expect("baseline external registry");
        assert!(baseline.contains("\"external_authorities\""));
        assert!(baseline.contains("\"id\":\"ci:proof\""));

        std::fs::write(
            &path,
            "row() { python3 - <<'PY'\nimport hashlib\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\nprint(hashlib.sha256(f'{IDENTITY_DOMAIN}:identity-v2'.encode()).hexdigest())\nPY\n}\n",
        )
        .expect("moved script");
        let moved = render_registry(&root, &[], &[owner], &[]).expect("moved external registry");
        let violations = schema_history_against(&moved, &baseline);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("changed schema fingerprint at retained version 1")),
            "{violations:?}"
        );
    }

    #[test]
    fn external_authority_coupling_requires_exact_producer_epoch() {
        let root = fixture_root("external-coupling");
        let path = root.join("scripts/ci/proof.sh");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let mut owner = ExternalOwner {
            id: "ci:proof".to_string(),
            path: "scripts/ci/proof.sh".to_string(),
            symbol: "row".to_string(),
            version: 1,
            domain: "org.frankensim.ci.proof.v1".to_string(),
        };
        let mut manifest = AuthorityManifest {
            required_ids: BTreeSet::from([owner.id.clone()]),
            external_owners: vec![owner.clone()],
            exemptions: Vec::new(),
        };
        std::fs::write(
            &path,
            "row() { python3 - <<'PY'\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\nprint(f'{IDENTITY_DOMAIN}:identity-v1')\nPY\n}\n",
        )
        .expect("baseline script");

        let (_, missing) = validate_external_couplings(&root, &manifest, &BTreeMap::new());
        assert!(
            missing.iter().any(|violation| violation
                .detail
                .contains("external golden surface is missing")),
            "an external authority needs an explicit coupling row: {missing:?}"
        );

        let baseline_fingerprint =
            external_owner_schema_fingerprint(&root, &owner, &[]).expect("baseline fingerprint");
        let mut surfaces = BTreeMap::from([(
            owner.id.clone(),
            CouplingSurface::External {
                file: owner.path.clone(),
                symbol: owner.symbol.clone(),
                version: owner.version,
                domain: owner.domain.clone(),
                schema_fingerprint: baseline_fingerprint,
            },
        )]);
        let (versions, baseline_violations) =
            validate_external_couplings(&root, &manifest, &surfaces);
        assert!(baseline_violations.is_empty(), "{baseline_violations:?}");
        assert_eq!(versions.get("ci:proof"), Some(&1));

        std::fs::write(
            &path,
            "row() { python3 - <<'PY'\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v1'\nprint(f'{IDENTITY_DOMAIN}:identity-v1-revised')\nPY\n}\n",
        )
        .expect("moved script");
        let (_, moved_violations) = validate_external_couplings(&root, &manifest, &surfaces);
        assert!(
            moved_violations
                .iter()
                .any(|violation| violation.detail.contains("producer fingerprint")),
            "moving the producer at a retained epoch must invalidate its coupling: {moved_violations:?}"
        );

        owner.version = 2;
        owner.domain = "org.frankensim.ci.proof.v2".to_string();
        manifest.external_owners = vec![owner.clone()];
        std::fs::write(
            &path,
            "row() { python3 - <<'PY'\nIDENTITY_DOMAIN = 'org.frankensim.ci.proof.v2'\nprint(f'{IDENTITY_DOMAIN}:identity-v2')\nPY\n}\n",
        )
        .expect("version-two script");
        let (_, bumped_violations) = validate_external_couplings(&root, &manifest, &surfaces);
        assert!(
            bumped_violations
                .iter()
                .any(|violation| violation.detail.contains("replace it with exact row")),
            "an authority bump must not silently reuse its old coupling: {bumped_violations:?}"
        );

        let bumped_fingerprint =
            external_owner_schema_fingerprint(&root, &owner, &[]).expect("bumped fingerprint");
        surfaces.insert(
            owner.id.clone(),
            CouplingSurface::External {
                file: owner.path.clone(),
                symbol: owner.symbol.clone(),
                version: owner.version,
                domain: owner.domain.clone(),
                schema_fingerprint: bumped_fingerprint,
            },
        );
        let (versions, updated_violations) =
            validate_external_couplings(&root, &manifest, &surfaces);
        assert!(updated_violations.is_empty(), "{updated_violations:?}");
        assert_eq!(versions.get("ci:proof"), Some(&2));
    }

    #[test]
    fn stale_external_coupling_without_authority_is_rejected() {
        let surfaces = BTreeMap::from([(
            "ci:retired".to_string(),
            CouplingSurface::External {
                file: "scripts/retired.sh".to_string(),
                symbol: "row".to_string(),
                version: 1,
                domain: "org.frankensim.ci.retired.v1".to_string(),
                schema_fingerprint:
                    "v1-0000000000000000000000000000000000000000000000000000000000000000"
                        .to_string(),
            },
        )]);
        let (_, violations) =
            validate_external_couplings(Path::new("."), &empty_authority_manifest(), &surfaces);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("stale external surface")),
            "retired external rows cannot retain replay authority: {violations:?}"
        );
    }

    #[test]
    fn shell_function_blocks_stop_at_the_exact_closing_brace() {
        let script = concat!(
            "outer() { python3 - <<'PY'\n",
            "embedded() {\n",
            "    print({\"digest\": \"inside-heredoc\"})\n",
            "}\n",
            "PY\n",
            "}\n",
            "seal_on_exit() { printf '%s\\n' done; }\n",
            "python3 - <<'PY'\n",
            "import hashlib\n",
            "from pathlib import Path\n",
            "IDENTITY_DOMAIN = 'org.frankensim.ci.top-level-proof.v1'\n",
            "digest = hashlib.sha256(IDENTITY_DOMAIN.encode()).hexdigest()\n",
            "Path('receipt').write_text(digest)\n",
            "PY\n",
        );
        let blocks = shell_function_blocks(script);
        assert_eq!(
            blocks
                .iter()
                .map(|(symbol, _)| symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["outer", "seal_on_exit"]
        );
        assert!(blocks[0].1.contains("inside-heredoc"));
        assert!(!blocks[1].1.contains("write_text"));

        let candidates = discover_script_candidates("scripts/ci/proof.sh", script);
        assert!(
            candidates
                .iter()
                .all(|candidate| candidate.symbol != "seal_on_exit"),
            "top-level script tail must not be attributed to the final function: {candidates:?}"
        );
        assert!(
            candidates
                .iter()
                .any(|candidate| candidate.symbol == "<script>"),
            "top-level identity producers must remain discoverable beside functions: {candidates:?}"
        );
    }

    #[test]
    fn external_child_schema_movement_is_history_guarded() {
        let root = fixture_root("external-child-movement");
        let path = root.join("scripts/ci/proof.sh");
        std::fs::create_dir_all(path.parent().expect("script parent")).expect("script directory");
        let owner = ExternalOwner {
            id: "ci:proof".to_string(),
            path: "scripts/ci/proof.sh".to_string(),
            symbol: "parent".to_string(),
            version: 1,
            domain: "org.frankensim.ci.proof.v1".to_string(),
        };
        let exemption = IdentityExemption {
            path: owner.path.clone(),
            symbol: "child".to_string(),
            reason: "child-digest-helper".to_string(),
            covered_by: owner.id.clone(),
        };
        std::fs::write(
            &path,
            "parent() { DOMAIN=org.frankensim.ci.proof.v1; printf '%s\\n' \"$DOMAIN\" | shasum -a 256 >> parent.receipt; }\nchild() { printf '%s\\n' child-v1 | shasum -a 256; }\n",
        )
        .expect("baseline child script");
        let baseline = render_registry(
            &root,
            &[],
            std::slice::from_ref(&owner),
            std::slice::from_ref(&exemption),
        )
        .expect("baseline child registry");

        std::fs::write(
            &path,
            "parent() { DOMAIN=org.frankensim.ci.proof.v1; printf '%s\\n' \"$DOMAIN\" | shasum -a 256 >> parent.receipt; }\nchild() { printf '%s\\n' child-v2 | shasum -a 256; }\n",
        )
        .expect("moved child script");
        let moved = render_registry(
            &root,
            &[],
            std::slice::from_ref(&owner),
            std::slice::from_ref(&exemption),
        )
        .expect("moved child registry");
        let violations = schema_history_against(&moved, &baseline);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("changed schema fingerprint at retained version 1")),
            "{violations:?}"
        );

        std::fs::write(
            &path,
            "parent() { DOMAIN=org.frankensim.ci.proof.v2; printf '%s\\n' \"$DOMAIN\" | shasum -a 256 >> parent.receipt; }\nchild() { printf '%s\\n' child-v2 | shasum -a 256; }\n",
        )
        .expect("versioned child script");
        let bumped_owner = ExternalOwner {
            version: 2,
            domain: "org.frankensim.ci.proof.v2".to_string(),
            ..owner
        };
        let bumped = render_registry(&root, &[], &[bumped_owner], &[exemption])
            .expect("versioned child registry");
        assert!(
            schema_history_against(&bumped, &baseline).is_empty(),
            "a parent version bump must authorize classified child movement"
        );
    }

    #[test]
    fn internal_child_exemption_moves_its_parent_schema_fingerprint() {
        let root = fixture_root("internal-child-fingerprint");
        let path = root.join("crates/mini/src/lib.rs");
        std::fs::create_dir_all(path.parent().expect("owner parent")).expect("owner directory");
        let baseline_source = format!(
            "{}\nfn child_digest() -> u64 {{ 1 }}\n",
            identity_source("mini:identity", "none", 1)
        );
        std::fs::write(&path, &baseline_source).expect("baseline owner");
        let (baseline_declarations, violations) =
            declaration_blocks("crates/mini/src/lib.rs", &baseline_source);
        assert!(violations.is_empty(), "{violations:?}");
        let exemption = IdentityExemption {
            path: "crates/mini/src/lib.rs".to_string(),
            symbol: "child_digest".to_string(),
            reason: "child-digest-helper".to_string(),
            covered_by: "mini:identity".to_string(),
        };
        let baseline = identity_schema_base_hash(
            &root,
            &baseline_declarations[0],
            &baseline_source,
            std::slice::from_ref(&exemption),
        )
        .expect("baseline parent schema");

        let moved_source = baseline_source.replace(
            "fn child_digest() -> u64 { 1 }",
            "fn child_digest() -> u64 { 2 }",
        );
        std::fs::write(&path, &moved_source).expect("moved owner");
        let (moved_declarations, violations) =
            declaration_blocks("crates/mini/src/lib.rs", &moved_source);
        assert!(violations.is_empty(), "{violations:?}");
        let moved =
            identity_schema_base_hash(&root, &moved_declarations[0], &moved_source, &[exemption])
                .expect("moved parent schema");
        assert_ne!(baseline, moved);
    }

    #[test]
    fn mutation_targets_require_the_exact_test_symbol() {
        let prefixed_only = "#[test]\nfn mutation_a_extended() { assert!(true); }\nfn mutation_a() { assert!(true); }\n";
        assert!(!has_test_function(prefixed_only, "mutation_a"));
        let exact = "#[test]\nfn mutation_a() { assert!(true); }\n";
        assert!(has_test_function(exact, "mutation_a"));
    }

    #[test]
    fn mutation_test_attribute_survives_adjacent_docs_comments_and_attributes() {
        let decorated = concat!(
            "#[test]\n",
            "/// Exact mutation fixture.\n",
            "#[allow(clippy::let_unit_value)]\n",
            "/* retained association\n",
            " * across a block comment\n",
            " */\n",
            "fn mutation_a() { assert!(true); }\n",
        );
        assert!(has_test_function(decorated, "mutation_a"));

        let detached = concat!(
            "#[test]\n",
            "fn earlier() {}\n",
            "/// This comment cannot carry the prior attribute forward.\n",
            "fn mutation_a() { assert!(true); }\n",
        );
        assert!(!has_test_function(detached, "mutation_a"));
    }

    #[test]
    fn ignored_cfg_gated_and_empty_tests_are_not_identity_evidence() {
        assert!(!has_test_function(
            "#[test]\n#[ignore]\nfn mutation_a() { assert!(true); }\n",
            "mutation_a"
        ));
        assert!(!has_test_function(
            "#[test]\n#[should_panic]\nfn mutation_a() { panic!(\"vacuous\"); }\n",
            "mutation_a"
        ));
        assert!(!has_test_function(
            "#[test]\n#[cfg(any())]\nfn mutation_a() { assert!(true); }\n",
            "mutation_a"
        ));
        assert!(!has_test_function(
            "#[test]\n#[cfg(\n    any()\n)]\nfn mutation_a() { assert!(true); }\n",
            "mutation_a"
        ));
        assert!(!has_test_function(
            "#[test]\nfn mutation_a() {}\n",
            "mutation_a"
        ));
        assert!(!has_test_function(
            "#[cfg(any())]\nmod disabled {\n#[test]\nfn mutation_a() { assert!(true); }\n}\n",
            "mutation_a"
        ));
        assert!(has_test_function(
            "#[cfg(test)]\nmod tests {\n#[test]\nfn mutation_a() { assert!(true); }\n}\n",
            "mutation_a"
        ));
    }

    #[test]
    fn cfg_disabled_parent_cannot_hide_an_identity_authority() {
        let source = format!(
            "#[cfg(any())]\nmod disabled {{\n{}\n}}\n",
            owner_source("", "a", "a:crates/mini/src/lib.rs#mutation_a")
        );
        let (declarations, violations) = declaration_blocks("crates/mini/src/lib.rs", &source);
        assert!(declarations.is_empty(), "{declarations:?}");
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("behind cfg/cfg_attr")),
            "{violations:?}"
        );
    }

    #[test]
    fn raw_fixture_mirrors_do_not_become_owner_declarations() {
        let mirrored = format!(
            "const MIRROR: &str = r###\"{}\"###;\n",
            owner_source("", "a", "none")
        );
        let (declarations, violations) = declaration_blocks("crates/mirror/src/lib.rs", &mirrored);
        assert!(declarations.is_empty());
        assert!(violations.is_empty());
        assert!(identity_marker_lines(&owner_source("", "a", "none")).contains(&4));
    }

    #[test]
    fn gate_source_neither_declares_nor_discovers_itself() {
        let gate_source = include_str!("identities.rs");
        let (declarations, violations) = declaration_blocks(GATE_IMPLEMENTATION_PATH, gate_source);
        assert!(declarations.is_empty(), "{declarations:?}");
        assert!(violations.is_empty(), "{violations:?}");

        let root = fixture_root("gate-self-discovery");
        let gate_path = root.join(GATE_IMPLEMENTATION_PATH);
        let producer_path = root.join("crates/demo/src/lib.rs");
        std::fs::create_dir_all(gate_path.parent().expect("gate parent")).expect("gate directory");
        std::fs::create_dir_all(producer_path.parent().expect("producer parent"))
            .expect("producer directory");
        let producer = r#"
fn persist_receipt(ledger: &mut Ledger, payload: &[u8]) {
    let root = hash_domain("demo", payload);
    ledger.insert(root, payload);
}
"#;
        std::fs::write(&gate_path, producer).expect("fixture gate source");
        std::fs::write(&producer_path, producer).expect("fixture producer source");
        let candidates = discover_identity_candidates(&root);
        assert_eq!(candidates.len(), 1, "{candidates:?}");
        assert_eq!(candidates[0].path, "crates/demo/src/lib.rs");
    }

    #[test]
    fn schema_constant_grammar_accepts_local_and_repo_relative_symbols() {
        let constants = parse_schema_constants("LOCAL_TAG,crates/shared/src/schema.rs#SHARED_TAG")
            .expect("local and repo-relative constants are canonical");
        assert_eq!(
            constants
                .iter()
                .map(SchemaConstant::canonical)
                .collect::<Vec<_>>(),
            vec![
                "LOCAL_TAG".to_string(),
                "crates/shared/src/schema.rs#SHARED_TAG".to_string(),
            ]
        );
        assert!(
            parse_schema_constants("none")
                .expect("none is allowed")
                .is_empty()
        );
        assert!(parse_schema_constants("../schema.rs#TAG").is_err());
    }

    #[test]
    fn qualified_domain_constant_resolves_and_mismatch_fails_closed() {
        let root = fixture_root("qualified-domain-constant");
        let shared_path = root.join("crates/shared/src/domain.rs");
        std::fs::create_dir_all(shared_path.parent().expect("shared parent"))
            .expect("shared directory");
        std::fs::write(
            &shared_path,
            "pub const SHARED_DOMAIN: &str = \"org.frankensim.mini.v1\";\n",
        )
        .expect("shared domain");
        let owner = identity_source("mini:qualified-domain", "none", 1)
            .replace(
                "pub const MINI_DOMAIN: &str = \"org.frankensim.mini.v1\";\n",
                "",
            )
            .replace(
                "domain_const=MINI_DOMAIN",
                "domain_const=crates/shared/src/domain.rs#SHARED_DOMAIN",
            );
        let (declarations, parse_violations) = declaration_blocks("crates/mini/src/lib.rs", &owner);
        assert!(parse_violations.is_empty(), "{parse_violations:?}");
        let index = RustSourceIndex::new(&owner);
        let inline_sources =
            BTreeMap::from([("crates/mini/src/lib.rs".to_string(), owner.as_str())]);
        let references = IdentityReferenceCache::build(&root, declarations.iter(), &inline_sources);
        let baseline = validate_owner_items(&declarations[0], &owner, &index, &references);
        assert!(
            baseline.is_empty(),
            "qualified domain source of truth must resolve: {baseline:?}"
        );

        std::fs::write(
            &shared_path,
            "pub const SHARED_DOMAIN: &str = \"org.frankensim.other.v1\";\n",
        )
        .expect("moved shared domain");
        let moved_references =
            IdentityReferenceCache::build(&root, declarations.iter(), &inline_sources);
        let moved = validate_owner_items(&declarations[0], &owner, &index, &moved_references);
        assert!(
            moved.iter().any(|violation| violation
                .detail
                .contains("must declare exact domain \"org.frankensim.mini.v1\"")),
            "same-version external domain mismatch must fail closed: {moved:?}"
        );

        std::fs::write(
            &shared_path,
            "pub const SHARED_DOMAIN: &[u8] = \"org.frankensim.mini.v1\";\n",
        )
        .expect("wrongly typed shared domain");
        let typed_references =
            IdentityReferenceCache::build(&root, declarations.iter(), &inline_sources);
        let typed = validate_owner_items(&declarations[0], &owner, &index, &typed_references);
        assert!(
            typed.iter().any(|violation| violation
                .detail
                .contains("must declare exact domain \"org.frankensim.mini.v1\"")),
            "qualified domain authority must be declared with exact &str type: {typed:?}"
        );
    }

    #[test]
    fn same_value_qualified_domain_target_swap_changes_history_fingerprint() {
        let root = fixture_root("qualified-domain-target-swap");
        let shared_path = root.join("crates/shared/src/domain.rs");
        std::fs::create_dir_all(shared_path.parent().expect("shared parent"))
            .expect("shared directory");
        std::fs::write(
            &shared_path,
            concat!(
                "pub const PRIMARY_DOMAIN: &str = \"org.frankensim.mini.v1\";\n",
                "pub const ALIAS_DOMAIN: &str = \"org.frankensim.mini.v1\";\n",
            ),
        )
        .expect("same-valued domain authorities");
        let baseline_source = identity_source("mini:qualified-domain-swap", "none", 1)
            .replace(
                "pub const MINI_DOMAIN: &str = \"org.frankensim.mini.v1\";\n",
                "",
            )
            .replace(
                "domain_const=MINI_DOMAIN",
                "domain_const=crates/shared/src/domain.rs#PRIMARY_DOMAIN",
            );
        let moved_source = baseline_source.replace("#PRIMARY_DOMAIN", "#ALIAS_DOMAIN");
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_source)]);
        let baseline_registry =
            render_registry(&root, &baseline, &[], &[]).expect("baseline registry renders");
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_source)]);
        let moved_registry =
            render_registry(&root, &moved, &[], &[]).expect("moved registry renders");
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "the canonical qualified domain target is part of the schema fingerprint"
        );
        let history = schema_history_against(&moved_registry, &baseline_registry);
        assert!(
            history.iter().any(|violation| violation
                .detail
                .contains("changed schema fingerprint at retained version 1")),
            "same-valued qualified domain authority swaps require a version bump: {history:?}"
        );
    }

    #[test]
    fn owner_and_lexical_item_moves_change_history_fingerprint() {
        let root = fixture_root("identity-target-locators");
        let source = identity_source("mini:target-locator", "none", 1);
        let baseline_owner =
            resolved_fixture(&root, [("crates/owner_a/src/lib.rs", source.clone())]);
        let moved_owner = resolved_fixture(&root, [("crates/owner_b/src/lib.rs", source.clone())]);
        assert_ne!(
            baseline_owner[0].schema_fingerprint, moved_owner[0].schema_fingerprint,
            "moving an owner file changes the retained schema target"
        );

        let baseline_const = source.replace(
            "pub const MINI_DOMAIN: &str = \"org.frankensim.mini.v1\";",
            "mod first { pub const MINI_DOMAIN: &str = \"org.frankensim.mini.v1\"; }",
        );
        let moved_const = baseline_const.replace("mod first", "mod second");
        let baseline_const = resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_const)]);
        let moved_const = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_const)]);
        assert_ne!(
            baseline_const[0].schema_fingerprint, moved_const[0].schema_fingerprint,
            "moving an identical domain const between modules changes its locator"
        );

        let baseline_source = source
            .replace("pub struct Mini {", "mod first { pub struct Mini {")
            .replace("}\nfn classify_mini", "} }\nfn classify_mini");
        let moved_source = baseline_source.replace("mod first", "mod second");
        let baseline_source =
            resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_source)]);
        let moved_source = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_source)]);
        assert_ne!(
            baseline_source[0].schema_fingerprint, moved_source[0].schema_fingerprint,
            "moving an identical source type between modules changes its locator"
        );
    }

    #[cfg(unix)]
    #[test]
    fn explicit_identity_inputs_reject_symlink_components() {
        use std::os::unix::fs::symlink;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "frankensim-identity-symlink-root-{}-{nonce}",
            std::process::id()
        ));
        let outside = std::env::temp_dir().join(format!(
            "frankensim-identity-symlink-outside-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir_all(root.join("schema")).expect("fixture root");
        std::fs::create_dir_all(&outside).expect("outside fixture");
        std::fs::write(outside.join("owner.rs"), "const TAG: u8 = 1;\n").expect("outside source");
        symlink(outside.join("owner.rs"), root.join("schema/owner.rs")).expect("file symlink");
        symlink(&outside, root.join("linked-directory")).expect("directory symlink");

        for relative in ["schema/owner.rs", "linked-directory/owner.rs"] {
            let error = checked_repo_file(&root, relative, "identity input")
                .expect_err("symlinked identity input must fail closed");
            assert!(error.contains("symlink component"), "{error}");
        }
    }

    #[test]
    fn identity_domain_accepts_only_an_exact_version_segment() {
        assert!(domain_carries_version("org.frankensim.demo.v8", 8));
        assert!(domain_carries_version("fs-package:v8:semantic-witness", 8));
        assert!(domain_carries_version("fs-recompute-node-v2", 2));
        assert!(!domain_carries_version(
            "fs-package:v80:semantic-witness",
            8
        ));
        assert!(!domain_carries_version(
            "fs-package:rev8:semantic-witness",
            8
        ));
    }

    #[test]
    fn exact_constant_lookup_rejects_comment_prefix_and_cfg_decoys() {
        let exact = r#"
// const ID_VERSION: u32 = 99;
const ID_VERSION_EXTRA: u32 = 88;
const ID_VERSION: u32 = 7;
const ID_DOMAIN: &str = "org.frankensim.demo.v7";
"#;
        assert_eq!(source_const_u32(exact, "ID_VERSION"), Some(7));
        assert_eq!(
            source_const_str(exact, "ID_DOMAIN"),
            Some("org.frankensim.demo.v7")
        );
        assert_eq!(const_declarations(exact, "ID_VERSION").len(), 1);

        let ambiguous = r#"
#[cfg(feature = "left")]
const ID_VERSION: u32 = 7;
#[cfg(not(feature = "left"))]
const ID_VERSION: u32 = 8;
"#;
        assert_eq!(const_declarations(ambiguous, "ID_VERSION").len(), 2);
        assert_eq!(source_const_u32(ambiguous, "ID_VERSION"), None);
    }

    #[test]
    fn schema_function_grammar_accepts_local_method_and_repo_relative_symbols() {
        let functions = parse_schema_functions(
            "local_helper,Codec::method,crates/shared/src/schema.rs#codec::helper",
        )
        .expect("local, method, and repo-relative function references are canonical");
        assert_eq!(functions.len(), 3);
        assert!(
            functions
                .iter()
                .any(|function| function.canonical() == "local_helper")
        );
        assert!(functions.iter().any(|function| {
            function.canonical() == "crates/shared/src/schema.rs#codec::helper"
        }));
        assert!(
            parse_schema_functions("none")
                .expect("none is allowed")
                .is_empty()
        );
        assert!(parse_schema_functions("../schema.rs#helper").is_err());
        assert!(parse_schema_functions("helper,helper").is_err());
    }

    #[test]
    fn cross_file_schema_function_body_movement_changes_fingerprint() {
        let root = fixture_root("cross-file-function-movement");
        let helper_path = root.join("crates/shared/src/schema.rs");
        std::fs::create_dir_all(helper_path.parent().expect("helper parent"))
            .expect("helper directory");
        let owner = identity_source_with_schema_functions(
            "mini:function",
            "crates/shared/src/schema.rs#helpers::semantic_tag",
        );
        std::fs::write(
            &helper_path,
            "pub mod helpers { pub fn semantic_tag() -> u64 { 1 } }\n",
        )
        .expect("baseline helper");
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner.clone())]);
        std::fs::write(
            &helper_path,
            "pub mod helpers { pub fn semantic_tag() -> u64 { 2 } }\n",
        )
        .expect("moved helper");
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner)]);
        assert_eq!(baseline[0].version, moved[0].version);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "same-version cross-file helper movement must change the fingerprint"
        );
    }

    #[test]
    fn cross_file_schema_function_signature_movement_changes_fingerprint() {
        let root = fixture_root("cross-file-function-signature-movement");
        let helper_path = root.join("crates/shared/src/schema.rs");
        std::fs::create_dir_all(helper_path.parent().expect("helper parent"))
            .expect("helper directory");
        let owner = identity_source_with_schema_functions(
            "mini:function-signature",
            "crates/shared/src/schema.rs#helpers::semantic_tag",
        );
        std::fs::write(
            &helper_path,
            "pub mod helpers { pub fn semantic_tag(value: u32) -> u64 { let _ = value; 1 } }\n",
        )
        .expect("baseline helper signature");
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner.clone())]);
        std::fs::write(
            &helper_path,
            "pub mod helpers { pub fn semantic_tag(value: u64) -> u64 { let _ = value; 1 } }\n",
        )
        .expect("moved helper signature");
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner)]);
        assert_eq!(baseline[0].version, moved[0].version);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "same-version function signature movement must change the fingerprint"
        );
    }

    #[test]
    fn cross_file_transitive_schema_callee_movement_changes_fingerprint() {
        let root = fixture_root("cross-file-transitive-function-movement");
        let helper_path = root.join("crates/shared/src/schema.rs");
        std::fs::create_dir_all(helper_path.parent().expect("helper parent"))
            .expect("helper directory");
        let owner = identity_source_with_schema_functions(
            "mini:transitive-function",
            "crates/shared/src/schema.rs#helpers::semantic_tag",
        );
        std::fs::write(
            &helper_path,
            "pub mod helpers { fn nested() -> u64 { 1 } pub fn semantic_tag() -> u64 { nested() } }\n",
        )
        .expect("baseline transitive helper");
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner.clone())]);
        std::fs::write(
            &helper_path,
            "pub mod helpers { fn nested() -> u64 { 2 } pub fn semantic_tag() -> u64 { nested() } }\n",
        )
        .expect("moved transitive helper");
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", owner)]);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "an unlisted but reachable local callee must remain in the schema fingerprint"
        );
    }

    #[test]
    fn missing_and_ambiguous_schema_functions_are_refused() {
        let root = fixture_root("invalid-schema-functions");
        let helper_path = root.join("crates/shared/src/schema.rs");
        std::fs::create_dir_all(helper_path.parent().expect("helper parent"))
            .expect("helper directory");
        std::fs::write(&helper_path, "pub fn present() {}\n").expect("present helper");

        let missing_source = identity_source_with_schema_functions(
            "mini:missing-function",
            "crates/shared/src/schema.rs#missing",
        );
        let (missing, missing_parse_violations) =
            declaration_blocks("crates/mini/src/lib.rs", &missing_source);
        assert!(missing_parse_violations.is_empty());
        let missing_error = identity_schema_base_hash(&root, &missing[0], &missing_source, &[])
            .expect_err("missing schema function must fail closed");
        assert!(missing_error.contains("found 0"), "{missing_error}");

        std::fs::write(
            &helper_path,
            concat!(
                "#[cfg(feature = \"left\")]\n",
                "pub fn helper() { let _ = 1_u64; }\n",
                "#[cfg(not(feature = \"left\"))]\n",
                "pub fn helper() { let _ = 2_u64; }\n",
            ),
        )
        .expect("ambiguous helpers");
        let ambiguous_source = identity_source_with_schema_functions(
            "mini:ambiguous-function",
            "crates/shared/src/schema.rs#helper",
        );
        let (ambiguous, ambiguous_parse_violations) =
            declaration_blocks("crates/mini/src/lib.rs", &ambiguous_source);
        assert!(ambiguous_parse_violations.is_empty());
        let ambiguous_error =
            identity_schema_base_hash(&root, &ambiguous[0], &ambiguous_source, &[])
                .expect_err("cfg-alternative schema functions must be ambiguous");
        assert!(ambiguous_error.contains("found 2"), "{ambiguous_error}");
    }

    #[test]
    fn same_version_schema_constant_movement_changes_fingerprint() {
        let root = fixture_root("constant-movement");
        let baseline = resolved_fixture(
            &root,
            [(
                "crates/mini/src/lib.rs",
                identity_source("mini:constant", "none", 1),
            )],
        );
        let moved = resolved_fixture(
            &root,
            [(
                "crates/mini/src/lib.rs",
                identity_source("mini:constant", "none", 2),
            )],
        );
        assert_eq!(baseline[0].version, moved[0].version);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "same-version schema constant movement must change the final fingerprint"
        );
    }

    #[test]
    fn mutation_evidence_body_movement_changes_fingerprint() {
        let root = fixture_root("mutation-evidence-movement");
        let baseline_source = identity_source("mini:mutation-evidence", "none", 1);
        let moved_source = baseline_source.replacen(
            "fn mutation_a() { assert_eq!(1_u8, 1_u8); }",
            "fn mutation_a() { assert_ne!(1_u8, 2_u8); }",
            1,
        );
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_source)]);
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_source)]);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "mutation evidence is part of the retained schema proof"
        );
    }

    #[test]
    fn same_field_name_with_moved_source_type_changes_fingerprint() {
        let root = fixture_root("source-field-type-movement");
        let baseline_source = identity_source("mini:source-type", "none", 1);
        let moved_source = baseline_source.replacen("pub a: u64", "pub a: u32", 1);
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_source)]);
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_source)]);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "source field type movement must not retain an identity fingerprint"
        );
    }

    #[test]
    fn moved_source_generic_header_changes_fingerprint() {
        let root = fixture_root("source-generic-header-movement");
        let baseline_source = identity_source("mini:source-header", "none", 1);
        let moved_source =
            baseline_source.replacen("pub struct Mini {", "pub struct Mini<T = u64> {", 1);
        let baseline = resolved_fixture(&root, [("crates/mini/src/lib.rs", baseline_source)]);
        let moved = resolved_fixture(&root, [("crates/mini/src/lib.rs", moved_source)]);
        assert_ne!(
            baseline[0].schema_fingerprint, moved[0].schema_fingerprint,
            "source generic/default movement must not retain an identity fingerprint"
        );
    }

    #[test]
    fn same_version_dependency_movement_changes_transitive_fingerprint() {
        let root = fixture_root("dependency-movement");
        let baseline = resolved_fixture(
            &root,
            [
                (
                    "crates/base/src/lib.rs",
                    identity_source("mini:base", "none", 1),
                ),
                (
                    "crates/dependent/src/lib.rs",
                    identity_source("mini:dependent", "mini:base", 1),
                ),
            ],
        );
        let moved = resolved_fixture(
            &root,
            [
                (
                    "crates/base/src/lib.rs",
                    identity_source("mini:base", "none", 2),
                ),
                (
                    "crates/dependent/src/lib.rs",
                    identity_source("mini:dependent", "mini:base", 1),
                ),
            ],
        );
        let baseline_dependent = baseline
            .iter()
            .find(|declaration| declaration.id == "mini:dependent")
            .expect("baseline dependent");
        let moved_dependent = moved
            .iter()
            .find(|declaration| declaration.id == "mini:dependent")
            .expect("moved dependent");
        assert_eq!(baseline_dependent.version, moved_dependent.version);
        assert_ne!(
            baseline_dependent.schema_fingerprint, moved_dependent.schema_fingerprint,
            "dependency movement must propagate into the dependent final fingerprint"
        );
    }

    #[test]
    fn missing_cycle_and_self_dependencies_are_refused() {
        let root = fixture_root("invalid-dependencies");
        let mut missing = vec![fixture_declaration(
            &root,
            "crates/dependent/src/lib.rs",
            &identity_source("mini:dependent", "mini:missing", 1),
        )];
        let missing_violations = resolve_schema_fingerprints(&mut missing);
        assert!(
            missing_violations
                .iter()
                .any(|violation| violation.detail.contains("does not exist")),
            "{missing_violations:?}"
        );

        let mut cycle = vec![
            fixture_declaration(
                &root,
                "crates/a/src/lib.rs",
                &identity_source("mini:a", "mini:b", 1),
            ),
            fixture_declaration(
                &root,
                "crates/b/src/lib.rs",
                &identity_source("mini:b", "mini:a", 1),
            ),
        ];
        let cycle_violations = resolve_schema_fingerprints(&mut cycle);
        assert!(
            cycle_violations
                .iter()
                .any(|violation| violation.detail.contains("cycle is forbidden")),
            "{cycle_violations:?}"
        );

        let (_, self_violations) = declaration_blocks(
            "crates/self/src/lib.rs",
            &identity_source("mini:self", "mini:self", 1),
        );
        assert!(
            self_violations
                .iter()
                .any(|violation| violation.detail.contains("self-dependency")),
            "{self_violations:?}"
        );
    }

    #[test]
    fn two_signal_discovery_rejects_unowned_rust_and_embedded_python_producers() {
        let hashed = discover_rust_candidates(
            "crates/demo/src/lib.rs",
            r#"
fn persist_receipt(ledger: &mut Ledger, payload: &[u8]) {
    let root = hash_domain("demo", payload);
    ledger.insert(root, payload);
}
"#,
        );
        assert_eq!(hashed.len(), 1, "hash-domain ledger producer is discovered");

        let canonical = discover_rust_candidates(
            "crates/demo/src/cache.rs",
            r#"
fn remember(cache: &mut Cache, value: Value) {
    let canonical_key = value.to_canonical_bytes();
    cache.insert(canonical_key, value);
}
"#,
        );
        assert_eq!(
            canonical.len(),
            1,
            "canonical cache-key producer is discovered"
        );

        let python = discover_script_candidates(
            "scripts/ci/snapshot.sh",
            r#"
snapshot_identity() {
  python3 <<'PY'
import hashlib
import pathlib
IDENTITY_DOMAIN = "org.frankensim.ci.snapshot.v1"
digest = hashlib.sha256(IDENTITY_DOMAIN.encode() + b"payload").hexdigest()
pathlib.Path("snapshot.json").write_text(digest)
PY
}
"#,
        );
        assert_eq!(python.len(), 1, "embedded Python snapshot is discovered");
        assert_eq!(python[0].symbol, "snapshot_identity");

        let manifest = AuthorityManifest {
            required_ids: BTreeSet::new(),
            external_owners: Vec::new(),
            exemptions: Vec::new(),
        };
        for candidates in [&hashed, &canonical, &python] {
            let violations =
                authority_violations_against(Path::new("."), &[], &manifest, candidates);
            assert!(
                violations
                    .iter()
                    .any(|violation| violation.detail.contains("is unregistered")),
                "{violations:?}"
            );
        }
    }

    #[test]
    fn mathematical_process_local_and_pid_log_names_are_not_authorities() {
        let candidates = discover_rust_candidates(
            "crates/demo/src/lib.rs",
            r#"
fn additive_identity() -> i64 { 0 }
fn process_local_identity() {
    let _state = std::collections::hash_map::DefaultHasher::new();
}
fn diagnostic_identity_log_path() -> String {
    format!("diagnostic-{}.log", std::process::id())
}
fn classify_demo_identity_fields(value: &Demo) {
    let Demo { identity: _ } = value;
}
fn unrelated_storage(ledger: &mut Ledger) {
    ledger.insert("ordinary", 1);
}
#[cfg(test)]
mod tests {
    fn gemm_identity(cache: &mut Cache) {
        cache.save("test-only-canonical-identity");
    }
}
"#,
        );
        assert!(candidates.is_empty(), "{candidates:?}");
    }

    #[test]
    fn duplicate_method_candidates_keep_their_qualified_owner() {
        let candidates = discover_rust_candidates(
            "crates/demo/src/policy.rs",
            r#"
trait Verify { fn verify(&self) -> Decision; }
struct DenyA;
struct DenyB;
impl Verify for DenyA {
    fn verify(&self) -> Decision {
        VerificationDecision::reject(hash_domain("deny-a", b"sentinel"))
    }
}
impl Verify for DenyB {
    fn verify(&self) -> Decision {
        VerificationDecision::reject(hash_domain("deny-b", b"sentinel"))
    }
}
"#,
        );
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["DenyA::verify", "DenyB::verify"])
        );
        assert!(!authority_target_matches(
            &candidates[0],
            "crates/demo/src/policy.rs",
            "verify"
        ));
    }

    #[test]
    fn array_signatures_and_one_line_generic_impls_keep_their_owner() {
        let source = r#"
struct ArrayOwner<const N: usize>;
impl<const N: usize> ArrayOwner<{ N }> { fn persist_identity<const M: usize>(cache: &mut Cache, payload: &[u8; M]) -> [u8; M] { let digest = hash_domain("array", payload); cache.insert(digest, payload); *payload } }
"#;
        let body = function_body(source, "ArrayOwner::persist_identity")
            .expect("array semicolons and const generics do not terminate the signature");
        assert!(body.contains("cache.insert"));
        let candidates = discover_rust_candidates("crates/demo/src/array.rs", source);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["ArrayOwner::persist_identity"]
        );
    }

    #[test]
    fn unrelated_module_sink_does_not_promote_a_pure_identity_helper() {
        let candidates = discover_rust_candidates(
            "crates/demo/src/modules.rs",
            r#"
mod pure {
    fn canonical_identity(payload: &[u8]) -> Digest {
        hash_domain("pure", payload)
    }
}
mod storage {
    fn remember_unrelated(cache: &mut Cache, value: Value) {
        cache.insert("ordinary", value);
    }
}
"#,
        );
        assert!(candidates.is_empty(), "{candidates:?}");
    }

    #[test]
    fn python_function_sinks_are_scoped_to_their_exact_definition() {
        let candidates = discover_script_candidates(
            "scripts/proof.py",
            r#"
def pure_snapshot_identity(payload):
    return hashlib.sha256(payload).hexdigest()

def write_unrelated(payload):
    pathlib.Path("ordinary.txt").write_text(payload)
"#,
        );
        assert!(candidates.is_empty(), "{candidates:?}");
    }

    #[test]
    fn external_owner_must_match_an_independently_discovered_producer() {
        let root = fixture_root("external-owner-independent-discovery");
        let path = root.join("crates/demo/src/lib.rs");
        std::fs::create_dir_all(path.parent().expect("demo parent")).expect("demo directory");
        std::fs::write(&path, "fn producer() { let _ = 1_u8; }\n").expect("producer source");
        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from(["demo:producer".to_string()]),
            external_owners: vec![ExternalOwner {
                id: "demo:producer".to_string(),
                path: "crates/demo/src/lib.rs".to_string(),
                symbol: "producer".to_string(),
                version: 1,
                domain: "org.frankensim.demo.producer.v1".to_string(),
            }],
            exemptions: Vec::new(),
        };
        let violations = authority_violations_against(&root, &[], &manifest, &[]);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("matches 0 independently discovered producers")),
            "{violations:?}"
        );
    }

    #[test]
    fn registered_child_helper_exemption_requires_and_accepts_its_parent() {
        let root = fixture_root("registered-child-authority");
        let path = root.join("crates/demo/src/lib.rs");
        std::fs::create_dir_all(path.parent().expect("demo parent")).expect("demo directory");
        let source = r#"
fn parent(cache: &mut Cache, value: Value) {
    let canonical_key = value.to_canonical_bytes();
    let digest = child(&canonical_key);
    cache.insert(digest, value);
}
fn child(payload: &[u8]) -> Digest {
    hash_domain("demo-child", payload)
}
"#;
        std::fs::write(&path, source).expect("authority source");
        let candidates = discover_rust_candidates("crates/demo/src/lib.rs", source);
        assert_eq!(candidates.len(), 1, "{candidates:?}");
        assert_eq!(candidates[0].symbol, "parent");
        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from(["demo:parent".to_string()]),
            external_owners: vec![ExternalOwner {
                id: "demo:parent".to_string(),
                path: "crates/demo/src/lib.rs".to_string(),
                symbol: "parent".to_string(),
                version: 1,
                domain: "org.frankensim.demo.parent.v1".to_string(),
            }],
            exemptions: vec![IdentityExemption {
                path: "crates/demo/src/lib.rs".to_string(),
                symbol: "child".to_string(),
                reason: "child-digest-helper".to_string(),
                covered_by: "demo:parent".to_string(),
            }],
        };
        let violations = authority_violations_against(&root, &[], &manifest, &candidates);
        assert!(violations.is_empty(), "{violations:?}");

        let mut unparented = manifest.clone();
        unparented.exemptions[0].covered_by = "none".to_string();
        let violations = authority_violations_against(&root, &[], &unparented, &candidates);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("must name a registered parent")),
            "{violations:?}"
        );

        let mut misclassified = manifest;
        misclassified.exemptions[0].reason = "identity-consumer-not-producer".to_string();
        let violations = authority_violations_against(&root, &[], &misclassified, &candidates);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("may not name parent")),
            "{violations:?}"
        );
    }

    #[test]
    fn stale_and_missing_exemption_targets_fail_closed() {
        let root = fixture_root("stale-authority-exemptions");
        let path = root.join("crates/demo/src/lib.rs");
        std::fs::create_dir_all(path.parent().expect("demo parent")).expect("demo directory");
        std::fs::write(&path, "fn harmless() {}\n").expect("stale source");
        let manifest = AuthorityManifest {
            required_ids: BTreeSet::from(["demo:missing".to_string()]),
            external_owners: Vec::new(),
            exemptions: vec![
                IdentityExemption {
                    path: "crates/demo/src/lib.rs".to_string(),
                    symbol: "harmless".to_string(),
                    reason: "diagnostic-only".to_string(),
                    covered_by: "none".to_string(),
                },
                IdentityExemption {
                    path: "crates/demo/src/lib.rs".to_string(),
                    symbol: "missing".to_string(),
                    reason: "diagnostic-only".to_string(),
                    covered_by: "none".to_string(),
                },
            ],
        };
        let violations = authority_violations_against(&root, &[], &manifest, &[]);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("demo:missing")
                    && violation.detail.contains("no owner declaration")),
            "{violations:?}"
        );
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("harmless is stale")),
            "{violations:?}"
        );
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("missing is invalid")),
            "{violations:?}"
        );
    }

    #[test]
    fn reasoned_external_exclusion_requires_one_exact_nonmovement_test() {
        let source = owner_source("", "a", "a:crates/mini/src/lib.rs#mutation_a")
            .replace(
                "\"excluded_fields=none\"",
                "\"excluded_fields=runtime.pid:diagnostic-envelope-only\"",
            )
            .replace(
                "\"nonsemantic_mutations=none\"",
                "\"nonsemantic_mutations=runtime.pid:crates/mini/src/lib.rs#external_pid_does_not_move\"",
            )
            .replace(
                "#[test]\nfn version_refuses() { assert_eq!(MINI_VERSION, 1); }",
                "#[test]\nfn external_pid_does_not_move() { assert!(true); }\n#[test]\nfn version_refuses() { assert_eq!(MINI_VERSION, 1); }",
            );
        let root = seed_fixture("reasoned-external-exclusion", &source, 1);
        let (_, violations) = load_declarations(&root);
        assert!(violations.is_empty(), "{violations:?}");

        let missing = source.replace(
            "\"nonsemantic_mutations=runtime.pid:crates/mini/src/lib.rs#external_pid_does_not_move\"",
            "\"nonsemantic_mutations=none\"",
        );
        let root = seed_fixture("reasoned-external-exclusion-missing", &missing, 1);
        let (_, violations) = load_declarations(&root);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("runtime.pid\" has no non-movement test")),
            "{violations:?}"
        );
    }

    #[test]
    fn seeded_unclassified_field_is_named() {
        let source = owner_source("pub b: u64,", "a", "a:crates/mini/src/lib.rs#mutation_a");
        let root = seed_fixture("unclassified", &source, 1);
        let (_, violations) = load_declarations(&root);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("unclassified source field Mini.b")),
            "{violations:?}"
        );
    }

    #[test]
    fn seeded_missing_mutation_is_named() {
        let source = owner_source("", "a,b", "a:crates/mini/src/lib.rs#mutation_a");
        let root = seed_fixture("missing-mutation", &source, 1);
        let (_, violations) = load_declarations(&root);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("semantic field \"b\" has no mutation")),
            "{violations:?}"
        );
    }

    #[test]
    fn stale_generated_bytes_and_coupling_version_are_refused() {
        let source = owner_source("", "a", "a:crates/mini/src/lib.rs#mutation_a");
        let stale_root = seed_fixture("stale", &source, 1);
        std::fs::write(stale_root.join(REGISTRY_FILE), "stale\n").expect("stale registry");
        assert!(
            check_identities(&stale_root)
                .iter()
                .any(|violation| violation.detail.contains("stale at byte"))
        );

        let coupling_root = seed_fixture("coupling", &source, 2);
        let (_, violations) = load_declarations(&coupling_root);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("version: 2")),
            "{violations:?}"
        );

        let missing_root = seed_fixture("missing-coupling", &source, 1);
        std::fs::write(
            missing_root.join("golden-couplings.json"),
            concat!(
                "{\n",
                "\"schema\": \"frankensim-golden-couplings-v1\",\n",
                "\"note\": \"fixture\",\n",
                "\"surfaces\": [\n],\n",
                "\"goldens\": [\n]\n",
                "}\n",
            ),
        )
        .expect("missing coupling fixture");
        let (_, violations) = load_declarations(&missing_root);
        assert!(
            violations.iter().any(|violation| {
                violation.detail.contains("add exact row")
                    && violation.detail.contains("\"schema_fingerprint\"")
            }),
            "missing coupling diagnostics must be directly actionable: {violations:?}"
        );
    }
}
