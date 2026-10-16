//! Shared canonical format for fs-la dependency-graph receipts.
//!
//! This module is independent of higher FrankenSim crates: xtask mints the
//! receipt and fs-la's build script parses and hashes the exact same grammar.
//! Its only workspace dependency is the shared `fs-blake3` implementation.

use std::collections::BTreeSet;
use std::fmt::Write as _;

pub const SCHEMA: &str = "fs-la-depgraph-receipt-v1";
pub const SCOPE: &str = "single-root-normal-build-fs-la-closure-v2";
/// Semantic version of the canonical dependency-graph receipt identity.
pub const DEPGRAPH_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Domain separating dependency-graph receipts from other BLAKE3 inputs.
pub const DEPGRAPH_RECEIPT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-la.depgraph-receipt.v1";
pub const MAX_RECEIPT_BYTES: usize = 1_048_576;
pub const MAX_PACKAGES: usize = 8_192;
pub const MAX_FEATURES: usize = 1_024;
pub const MAX_STRING_BYTES: usize = 8_192;

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)] // consumed as source text by xtask in build-script mode
pub const DEPGRAPH_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-la:depgraph-receipt",
    "version_const=DEPGRAPH_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-la.depgraph-receipt.v1",
    "domain_const=DEPGRAPH_RECEIPT_IDENTITY_DOMAIN",
    "encoder=content_identity",
    "encoder_helpers=content_identity_with_domain,emit",
    "schema_functions=validate,parse,is_lower_hex_64,validate_string,validate_feature_set,validate_identity,push_json_string,push_optional_string,push_string_array,push_identity,Parser::new,Parser::expect,Parser::byte,Parser::hex4,Parser::string,Parser::optional_string,Parser::boolean,Parser::strings,Parser::identity,depgraph_receipt_identity_version_is_supported,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=DEPGRAPH_RECEIPT_IDENTITY_VERSION,DEPGRAPH_RECEIPT_IDENTITY_DOMAIN,SCHEMA,SCOPE,MAX_RECEIPT_BYTES,MAX_PACKAGES,MAX_FEATURES,MAX_STRING_BYTES",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=Receipt,CargoIdentity,RootRow,SelectionRow,PackageRow,PackageIdentity",
    "source_fields=Receipt.cargo:semantic,Receipt.root:semantic,Receipt.selection:semantic,Receipt.packages:semantic,CargoIdentity.executable_digest:derived:nested-under-receipt-cargo,CargoIdentity.version:derived:nested-under-receipt-cargo,RootRow.identity:derived:nested-under-receipt-root,RootRow.features:derived:nested-under-receipt-root,SelectionRow.target:derived:nested-under-receipt-selection,SelectionRow.features:derived:nested-under-receipt-selection,SelectionRow.all_features:derived:nested-under-receipt-selection,SelectionRow.default_features:derived:nested-under-receipt-selection,PackageRow.identity:derived:nested-under-receipt-packages,PackageRow.features:derived:nested-under-receipt-packages,PackageIdentity.name:derived:nested-under-root-or-package-row,PackageIdentity.version:derived:nested-under-root-or-package-row,PackageIdentity.package_id:derived:nested-under-root-or-package-row,PackageIdentity.source_id:derived:nested-under-root-or-package-row,PackageIdentity.path_digest:derived:nested-under-root-or-package-row",
    "source_bindings=Receipt.cargo>cargo-identity,Receipt.root>root-identity,Receipt.selection>root-selection,Receipt.packages>ordered-package-closure",
    "external_semantic_fields=artifact-domain,receipt-schema,receipt-scope",
    "semantic_fields=artifact-domain,receipt-schema,receipt-scope,cargo-identity,root-identity,root-selection,ordered-package-closure",
    "excluded_fields=none",
    "consumers=crates/fs-la/build.rs#add_depgraph_evidence,FRANKENSIM_DEPGRAPH_RECEIPT,FS_LA_GEMM_DEPGRAPH_RECEIPT_DIGEST",
    "mutations=artifact-domain:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_domain_moves_identity,receipt-schema:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_schema_moves_identity,receipt-scope:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_scope_moves_identity,cargo-identity:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_cargo_identity_moves_identity,root-identity:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_root_identity_moves_identity,root-selection:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_selection_moves_identity,ordered-package-closure:crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_package_closure_moves_identity",
    "nonsemantic_mutations=none",
    "field_guard=classify_depgraph_receipt_identity_fields",
    "transport_guard=parse",
    "version_guard=crates/fs-la/depgraph_receipt_format.rs#depgraph_receipt_identity_version_fails_closed",
    "coupling_surface=fs-la:depgraph-receipt",
];

/// Whether a retained dependency-graph receipt uses the one identity version
/// accepted by this build and its shared xtask emitter.
#[must_use]
#[allow(dead_code)] // exercised by the module identity-version guard
pub const fn depgraph_receipt_identity_version_is_supported(declared: u32) -> bool {
    declared == DEPGRAPH_RECEIPT_IDENTITY_VERSION
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CargoIdentity {
    pub executable_digest: String,
    pub version: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct PackageIdentity {
    pub name: String,
    pub version: String,
    pub package_id: String,
    pub source_id: Option<String>,
    pub path_digest: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootRow {
    pub identity: PackageIdentity,
    pub features: BTreeSet<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionRow {
    pub target: Option<String>,
    pub features: BTreeSet<String>,
    pub all_features: bool,
    pub default_features: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct PackageRow {
    pub identity: PackageIdentity,
    pub features: BTreeSet<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Receipt {
    pub cargo: CargoIdentity,
    pub root: RootRow,
    pub selection: SelectionRow,
    pub packages: Vec<PackageRow>,
}

#[allow(dead_code)] // exhaustive source-shape guard consumed by xtask
fn classify_depgraph_receipt_identity_fields(
    receipt: &Receipt,
    cargo: &CargoIdentity,
    root: &RootRow,
    selection: &SelectionRow,
    package: &PackageRow,
    identity: &PackageIdentity,
) {
    let Receipt {
        cargo: _,
        root: _,
        selection: _,
        packages: _,
    } = receipt;
    let CargoIdentity {
        executable_digest: _,
        version: _,
    } = cargo;
    let RootRow {
        identity: _,
        features: _,
    } = root;
    let SelectionRow {
        target: _,
        features: _,
        all_features: _,
        default_features: _,
    } = selection;
    let PackageRow {
        identity: _,
        features: _,
    } = package;
    let PackageIdentity {
        name: _,
        version: _,
        package_id: _,
        source_id: _,
        path_digest: _,
    } = identity;
}

fn is_lower_hex_64(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn validate_string(label: &str, value: &str, allow_controls: bool) -> Result<(), String> {
    if value.is_empty() || value.len() > MAX_STRING_BYTES {
        return Err(format!(
            "{label} must contain 1..={MAX_STRING_BYTES} UTF-8 bytes"
        ));
    }
    if !allow_controls && value.chars().any(char::is_control) {
        return Err(format!("{label} contains a control character"));
    }
    Ok(())
}

fn validate_feature_set(label: &str, features: &BTreeSet<String>) -> Result<(), String> {
    if features.len() > MAX_FEATURES {
        return Err(format!(
            "{label} exceeds the {MAX_FEATURES}-feature receipt bound"
        ));
    }
    for feature in features {
        if feature.is_empty()
            || feature.len() > 256
            || !feature.bytes().all(|byte| byte.is_ascii_graphic())
        {
            return Err(format!("{label} contains invalid feature {feature:?}"));
        }
    }
    Ok(())
}

fn validate_identity(identity: &PackageIdentity) -> Result<(), String> {
    validate_string("package name", &identity.name, false)?;
    validate_string("package version", &identity.version, false)?;
    if !identity
        .name
        .bytes()
        .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
        || !identity.version.bytes().all(|byte| byte.is_ascii_graphic())
    {
        return Err("package name/version is outside Cargo's receipt grammar".to_string());
    }
    validate_string("Cargo metadata package id", &identity.package_id, false)?;
    let expected_suffix = format!("#{}@{}", identity.name, identity.version);
    if !identity.package_id.ends_with(&expected_suffix) {
        return Err("package id does not bind its declared name/version".to_string());
    }
    if let Some(source) = &identity.source_id {
        validate_string("Cargo metadata source id", source, false)?;
        if identity.path_digest.is_some() {
            return Err("non-path package carries a path digest".to_string());
        }
        if !identity.package_id.starts_with(source) {
            return Err("package id does not bind its Cargo metadata source id".to_string());
        }
    } else {
        let digest = identity
            .path_digest
            .as_deref()
            .ok_or_else(|| "path package lacks its source/build digest".to_string())?;
        if !is_lower_hex_64(digest) {
            return Err("path package digest is not canonical lowercase BLAKE3".to_string());
        }
        let expected = format!(
            "path+blake3:{digest}#{}@{}",
            identity.name, identity.version
        );
        if identity.package_id != expected {
            return Err("normalized path package id is not bound to its digest".to_string());
        }
    }
    Ok(())
}

pub fn validate(receipt: &Receipt) -> Result<(), String> {
    if !is_lower_hex_64(&receipt.cargo.executable_digest) {
        return Err("Cargo executable digest is not canonical lowercase BLAKE3".to_string());
    }
    validate_string("Cargo version identity", &receipt.cargo.version, true)?;
    validate_identity(&receipt.root.identity)?;
    validate_feature_set("root resolved features", &receipt.root.features)?;
    if receipt.root.identity.name.is_empty() {
        return Err("receipt root has no package name".to_string());
    }
    if let Some(target) = &receipt.selection.target {
        validate_string("target triple", target, false)?;
        if !target
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
        {
            return Err("target triple is outside the canonical machine-token grammar".to_string());
        }
    }
    validate_feature_set("selected features", &receipt.selection.features)?;
    if receipt.selection.features.iter().any(|feature| {
        !feature.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'+' | b'.' | b'/' | b'?')
        })
    }) {
        return Err("selected feature is outside the canonical Cargo selector grammar".to_string());
    }
    if receipt.selection.all_features && !receipt.selection.features.is_empty() {
        return Err("all_features and explicit features cannot coexist".to_string());
    }
    if receipt.packages.is_empty() || receipt.packages.len() > MAX_PACKAGES {
        return Err(format!(
            "receipt must contain 1..={MAX_PACKAGES} package-unit rows"
        ));
    }
    let mut previous = None;
    let mut saw_fs_la = false;
    let mut saw_root_unit = false;
    for row in &receipt.packages {
        validate_identity(&row.identity)?;
        validate_feature_set("resolved package features", &row.features)?;
        if previous.is_some_and(|prior| prior >= row) {
            return Err("package-unit rows are not strictly sorted/unique".to_string());
        }
        if row.identity.name == "fs-la" {
            saw_fs_la = true;
        }
        if row.identity == receipt.root.identity && row.features == receipt.root.features {
            saw_root_unit = true;
        }
        previous = Some(row);
    }
    if !saw_fs_la {
        return Err("receipt contains no metadata-resolved fs-la unit".to_string());
    }
    if receipt.root.identity.name == "fs-la" && !saw_root_unit {
        return Err("fs-la root row is absent from its package-unit rows".to_string());
    }
    Ok(())
}

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{20}'..='\u{7e}' => out.push(character),
            _ => {
                let codepoint = u32::from(character);
                if codepoint <= 0xffff {
                    let _ = write!(out, "\\u{codepoint:04x}");
                } else {
                    let scalar = codepoint - 0x1_0000;
                    let high = 0xd800 + (scalar >> 10);
                    let low = 0xdc00 + (scalar & 0x3ff);
                    let _ = write!(out, "\\u{high:04x}\\u{low:04x}");
                }
            }
        }
    }
    out.push('"');
}

fn push_optional_string(out: &mut String, value: Option<&str>) {
    if let Some(value) = value {
        push_json_string(out, value);
    } else {
        out.push_str("null");
    }
}

fn push_string_array(out: &mut String, values: &BTreeSet<String>) {
    out.push('[');
    let mut separator = "";
    for value in values {
        out.push_str(separator);
        push_json_string(out, value);
        separator = ",";
    }
    out.push(']');
}

fn push_identity(out: &mut String, identity: &PackageIdentity) {
    out.push_str("\"name\":");
    push_json_string(out, &identity.name);
    out.push_str(",\"version\":");
    push_json_string(out, &identity.version);
    out.push_str(",\"package_id\":");
    push_json_string(out, &identity.package_id);
    out.push_str(",\"source_id\":");
    push_optional_string(out, identity.source_id.as_deref());
    out.push_str(",\"path_digest\":");
    push_optional_string(out, identity.path_digest.as_deref());
}

pub fn emit(receipt: &Receipt) -> Result<String, String> {
    validate(receipt)?;
    let mut out = String::new();
    out.push_str("{\"schema\":");
    push_json_string(&mut out, SCHEMA);
    out.push_str(",\"scope\":");
    push_json_string(&mut out, SCOPE);
    out.push_str(",\"cargo\":{\"executable_digest\":");
    push_json_string(&mut out, &receipt.cargo.executable_digest);
    out.push_str(",\"version\":");
    push_json_string(&mut out, &receipt.cargo.version);
    out.push_str("},\"root\":{");
    push_identity(&mut out, &receipt.root.identity);
    out.push_str(",\"features\":");
    push_string_array(&mut out, &receipt.root.features);
    out.push_str("},\"selection\":{\"target\":");
    push_optional_string(&mut out, receipt.selection.target.as_deref());
    out.push_str(",\"features\":");
    push_string_array(&mut out, &receipt.selection.features);
    let _ = write!(
        out,
        ",\"all_features\":{},\"default_features\":{}",
        receipt.selection.all_features, receipt.selection.default_features
    );
    out.push_str("},\"packages\":[");
    let mut separator = "";
    for row in &receipt.packages {
        out.push_str(separator);
        out.push('{');
        push_identity(&mut out, &row.identity);
        out.push_str(",\"features\":");
        push_string_array(&mut out, &row.features);
        out.push('}');
        separator = ",";
    }
    out.push_str("]}");
    if out.len() > MAX_RECEIPT_BYTES {
        return Err(format!(
            "canonical receipt is {} bytes, exceeding the {MAX_RECEIPT_BYTES}-byte bound",
            out.len()
        ));
    }
    debug_assert!(out.is_ascii());
    Ok(out)
}

/// Content identity of one validated receipt under the canonical emitter.
pub fn content_identity(receipt: &Receipt) -> Result<fs_blake3::ContentHash, String> {
    content_identity_with_domain(DEPGRAPH_RECEIPT_IDENTITY_DOMAIN, receipt)
}

fn content_identity_with_domain(
    domain: &str,
    receipt: &Receipt,
) -> Result<fs_blake3::ContentHash, String> {
    let canonical = emit(receipt)?;
    Ok(fs_blake3::hash_domain(domain, canonical.as_bytes()))
}

struct Parser<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Result<Self, String> {
        if input.is_empty() || input.len() > MAX_RECEIPT_BYTES || !input.is_ascii() {
            return Err(format!(
                "receipt must be non-empty canonical ASCII within {MAX_RECEIPT_BYTES} bytes"
            ));
        }
        Ok(Self {
            bytes: input.as_bytes(),
            cursor: 0,
        })
    }

    fn expect(&mut self, literal: &[u8]) -> Result<(), String> {
        if self.bytes.get(self.cursor..self.cursor + literal.len()) == Some(literal) {
            self.cursor += literal.len();
            Ok(())
        } else {
            Err(format!(
                "expected canonical token {:?} at byte {}",
                String::from_utf8_lossy(literal),
                self.cursor
            ))
        }
    }

    fn byte(&mut self) -> Result<u8, String> {
        let byte = *self
            .bytes
            .get(self.cursor)
            .ok_or_else(|| "unexpected end of receipt".to_string())?;
        self.cursor += 1;
        Ok(byte)
    }

    fn hex4(&mut self) -> Result<u16, String> {
        let mut value = 0u16;
        for _ in 0..4 {
            let byte = self.byte()?;
            let digit = match byte {
                b'0'..=b'9' => u16::from(byte - b'0'),
                b'a'..=b'f' => u16::from(byte - b'a' + 10),
                _ => return Err("Unicode escape requires lowercase canonical hex".to_string()),
            };
            value = (value << 4) | digit;
        }
        Ok(value)
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect(b"\"")?;
        let mut out = String::new();
        loop {
            let byte = self.byte()?;
            match byte {
                b'"' => break,
                b'\\' => match self.byte()? {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'b' => out.push('\u{08}'),
                    b'f' => out.push('\u{0c}'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        let high = self.hex4()?;
                        let codepoint = if (0xd800..=0xdbff).contains(&high) {
                            self.expect(b"\\u")?;
                            let low = self.hex4()?;
                            if !(0xdc00..=0xdfff).contains(&low) {
                                return Err("invalid low surrogate".to_string());
                            }
                            0x1_0000
                                + ((u32::from(high) - 0xd800) << 10)
                                + (u32::from(low) - 0xdc00)
                        } else {
                            if (0xdc00..=0xdfff).contains(&high) {
                                return Err("unpaired low surrogate".to_string());
                            }
                            u32::from(high)
                        };
                        if (0x20..=0x7e).contains(&codepoint)
                            || matches!(codepoint, 8 | 9 | 10 | 12 | 13)
                        {
                            return Err("non-canonical JSON Unicode escape".to_string());
                        }
                        out.push(
                            char::from_u32(codepoint)
                                .ok_or_else(|| "invalid Unicode scalar".to_string())?,
                        );
                    }
                    _ => return Err("unsupported JSON escape".to_string()),
                },
                0x20..=0x7e => out.push(char::from(byte)),
                _ => return Err("receipt JSON must be canonical ASCII".to_string()),
            }
            if out.len() > MAX_STRING_BYTES {
                return Err(format!("receipt string exceeds {MAX_STRING_BYTES} bytes"));
            }
        }
        Ok(out)
    }

    fn optional_string(&mut self) -> Result<Option<String>, String> {
        if self.bytes.get(self.cursor..self.cursor + 4) == Some(b"null") {
            self.cursor += 4;
            Ok(None)
        } else {
            self.string().map(Some)
        }
    }

    fn boolean(&mut self) -> Result<bool, String> {
        if self.bytes.get(self.cursor..self.cursor + 4) == Some(b"true") {
            self.cursor += 4;
            Ok(true)
        } else if self.bytes.get(self.cursor..self.cursor + 5) == Some(b"false") {
            self.cursor += 5;
            Ok(false)
        } else {
            Err("expected canonical boolean".to_string())
        }
    }

    fn strings(&mut self) -> Result<BTreeSet<String>, String> {
        self.expect(b"[")?;
        let mut values = BTreeSet::new();
        if self.bytes.get(self.cursor) == Some(&b']') {
            self.cursor += 1;
            return Ok(values);
        }
        loop {
            if values.len() >= MAX_FEATURES {
                return Err(format!("string array exceeds {MAX_FEATURES} entries"));
            }
            let value = self.string()?;
            if values.last().is_some_and(|previous| previous >= &value) || !values.insert(value) {
                return Err("string array is not strictly sorted/unique".to_string());
            }
            match self.byte()? {
                b',' => {}
                b']' => break,
                _ => return Err("expected comma or array terminator".to_string()),
            }
        }
        Ok(values)
    }

    fn identity(&mut self) -> Result<PackageIdentity, String> {
        self.expect(b"\"name\":")?;
        let name = self.string()?;
        self.expect(b",\"version\":")?;
        let version = self.string()?;
        self.expect(b",\"package_id\":")?;
        let package_id = self.string()?;
        self.expect(b",\"source_id\":")?;
        let source_id = self.optional_string()?;
        self.expect(b",\"path_digest\":")?;
        let path_digest = self.optional_string()?;
        Ok(PackageIdentity {
            name,
            version,
            package_id,
            source_id,
            path_digest,
        })
    }
}

pub fn parse(input: &str) -> Result<Receipt, String> {
    let mut parser = Parser::new(input)?;
    parser.expect(b"{\"schema\":")?;
    if parser.string()? != SCHEMA {
        return Err("unexpected receipt schema".to_string());
    }
    parser.expect(b",\"scope\":")?;
    if parser.string()? != SCOPE {
        return Err("unexpected receipt scope".to_string());
    }
    parser.expect(b",\"cargo\":{\"executable_digest\":")?;
    let executable_digest = parser.string()?;
    parser.expect(b",\"version\":")?;
    let cargo_version = parser.string()?;
    parser.expect(b"},\"root\":{")?;
    let root_identity = parser.identity()?;
    parser.expect(b",\"features\":")?;
    let root_features = parser.strings()?;
    parser.expect(b"},\"selection\":{\"target\":")?;
    let target = parser.optional_string()?;
    parser.expect(b",\"features\":")?;
    let selected_features = parser.strings()?;
    parser.expect(b",\"all_features\":")?;
    let all_features = parser.boolean()?;
    parser.expect(b",\"default_features\":")?;
    let default_features = parser.boolean()?;
    parser.expect(b"},\"packages\":[")?;
    let mut packages = Vec::new();
    if parser.bytes.get(parser.cursor) == Some(&b']') {
        parser.cursor += 1;
    } else {
        loop {
            if packages.len() >= MAX_PACKAGES {
                return Err(format!("package array exceeds {MAX_PACKAGES} rows"));
            }
            parser.expect(b"{")?;
            let identity = parser.identity()?;
            parser.expect(b",\"features\":")?;
            let features = parser.strings()?;
            parser.expect(b"}")?;
            packages.push(PackageRow { identity, features });
            match parser.byte()? {
                b',' => {}
                b']' => break,
                _ => return Err("expected comma or package-array terminator".to_string()),
            }
        }
    }
    parser.expect(b"}")?;
    if parser.cursor != parser.bytes.len() {
        return Err("trailing bytes after receipt".to_string());
    }
    let receipt = Receipt {
        cargo: CargoIdentity {
            executable_digest,
            version: cargo_version,
        },
        root: RootRow {
            identity: root_identity,
            features: root_features,
        },
        selection: SelectionRow {
            target,
            features: selected_features,
            all_features,
            default_features,
        },
        packages,
    };
    validate(&receipt)?;
    if emit(&receipt)? != input {
        return Err("receipt is valid JSON but not the canonical byte form".to_string());
    }
    Ok(receipt)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path_identity(name: &str) -> PackageIdentity {
        let digest = "11".repeat(32);
        PackageIdentity {
            name: name.to_string(),
            version: "0.0.1".to_string(),
            package_id: format!("path+blake3:{digest}#{name}@0.0.1"),
            source_id: None,
            path_digest: Some(digest),
        }
    }

    fn sample() -> Receipt {
        let identity = path_identity("fs-la");
        Receipt {
            cargo: CargoIdentity {
                executable_digest: "22".repeat(32),
                version: "cargo 1.2.3\ncontrol=\u{0001} unicode=\u{2603}".to_string(),
            },
            root: RootRow {
                identity: identity.clone(),
                features: BTreeSet::new(),
            },
            selection: SelectionRow {
                target: None,
                features: BTreeSet::new(),
                all_features: false,
                default_features: true,
            },
            packages: vec![PackageRow {
                identity,
                features: BTreeSet::new(),
            }],
        }
    }

    #[test]
    fn depgraph_receipt_domain_moves_identity() {
        let receipt = sample();
        assert_ne!(
            content_identity(&receipt).expect("current identity"),
            content_identity_with_domain(
                "org.frankensim.fs-la.depgraph-receipt.v1.alternate",
                &receipt,
            )
            .expect("alternate identity")
        );
    }

    #[test]
    fn depgraph_receipt_schema_moves_identity() {
        let canonical = emit(&sample()).expect("canonical receipt");
        let changed = canonical.replacen(SCHEMA, "fs-la-depgraph-receipt-v1-alternate", 1);
        assert_ne!(
            fs_blake3::hash_domain(DEPGRAPH_RECEIPT_IDENTITY_DOMAIN, canonical.as_bytes()),
            fs_blake3::hash_domain(DEPGRAPH_RECEIPT_IDENTITY_DOMAIN, changed.as_bytes())
        );
    }

    #[test]
    fn depgraph_receipt_scope_moves_identity() {
        let canonical = emit(&sample()).expect("canonical receipt");
        let changed = canonical.replacen(SCOPE, "single-root-normal-build-fs-la-alternate", 1);
        assert_ne!(
            fs_blake3::hash_domain(DEPGRAPH_RECEIPT_IDENTITY_DOMAIN, canonical.as_bytes()),
            fs_blake3::hash_domain(DEPGRAPH_RECEIPT_IDENTITY_DOMAIN, changed.as_bytes())
        );
    }

    #[test]
    fn depgraph_receipt_cargo_identity_moves_identity() {
        let baseline = sample();
        let mut changed = baseline.clone();
        changed.cargo.executable_digest = "44".repeat(32);
        assert_ne!(
            content_identity(&baseline).expect("baseline"),
            content_identity(&changed).expect("changed Cargo identity")
        );
    }

    #[test]
    fn depgraph_receipt_root_identity_moves_identity() {
        let baseline = sample();
        let mut changed = baseline.clone();
        changed.root.identity = path_identity("workspace-root");
        assert_ne!(
            content_identity(&baseline).expect("baseline"),
            content_identity(&changed).expect("changed root identity")
        );
    }

    #[test]
    fn depgraph_receipt_selection_moves_identity() {
        let baseline = sample();
        let mut changed = baseline.clone();
        changed.selection.target = Some("aarch64-apple-darwin".to_string());
        assert_ne!(
            content_identity(&baseline).expect("baseline"),
            content_identity(&changed).expect("changed selection")
        );
    }

    #[test]
    fn depgraph_receipt_package_closure_moves_identity() {
        let baseline = sample();
        let mut changed = baseline.clone();
        changed.packages.push(PackageRow {
            identity: path_identity("fs-simd"),
            features: BTreeSet::new(),
        });
        changed.packages.sort();
        assert_ne!(
            content_identity(&baseline).expect("baseline"),
            content_identity(&changed).expect("changed package closure")
        );
    }

    #[test]
    fn depgraph_receipt_identity_version_fails_closed() {
        assert_eq!(DEPGRAPH_RECEIPT_IDENTITY_VERSION, 1);
        assert_eq!(
            DEPGRAPH_RECEIPT_IDENTITY_DOMAIN,
            "org.frankensim.fs-la.depgraph-receipt.v1"
        );
        assert!(depgraph_receipt_identity_version_is_supported(1));
        assert!(!depgraph_receipt_identity_version_is_supported(0));
        assert!(!depgraph_receipt_identity_version_is_supported(2));
    }

    #[test]
    fn emitter_and_build_parser_interoperate_with_control_escapes() {
        let receipt = sample();
        let encoded = emit(&receipt).expect("emit");
        assert!(encoded.is_ascii());
        assert!(encoded.contains("\\ncontrol=\\u0001 unicode=\\u2603"));
        assert_eq!(parse(&encoded).expect("parse"), receipt);
    }

    #[test]
    fn parser_rejects_noncanonical_escapes_and_bounds() {
        let encoded = emit(&sample()).expect("emit");
        let noncanonical = encoded.replacen("cargo", "\\u0063argo", 1);
        assert!(parse(&noncanonical).is_err());
        let oversized = "x".repeat(MAX_RECEIPT_BYTES + 1);
        assert!(parse(&oversized).is_err());
        let mut oversized_field = sample();
        oversized_field.cargo.version = "x".repeat(MAX_STRING_BYTES + 1);
        assert!(emit(&oversized_field).is_err());
    }

    #[test]
    fn parser_rejects_missing_path_digest_and_duplicate_units() {
        let mut missing = sample();
        missing.root.identity.path_digest = None;
        assert!(emit(&missing).is_err());
        let mut duplicate = sample();
        duplicate.packages.push(duplicate.packages[0].clone());
        assert!(emit(&duplicate).is_err());
    }

    #[test]
    fn feature_source_and_cargo_drift_change_stable_replay_bytes() {
        let baseline = sample();
        let first = emit(&baseline).expect("first");
        assert_eq!(first, emit(&baseline).expect("stable replay"));

        let mut feature_drift = baseline.clone();
        feature_drift.root.features.insert("simd".to_string());
        feature_drift.packages[0]
            .features
            .insert("simd".to_string());
        assert_ne!(first, emit(&feature_drift).expect("feature drift"));

        let mut source_drift = baseline.clone();
        let digest = "33".repeat(32);
        source_drift.root.identity.path_digest = Some(digest.clone());
        source_drift.root.identity.package_id = format!("path+blake3:{digest}#fs-la@0.0.1");
        source_drift.packages[0].identity = source_drift.root.identity.clone();
        assert_ne!(first, emit(&source_drift).expect("source drift"));

        let mut cargo_drift = baseline;
        cargo_drift.cargo.executable_digest = "44".repeat(32);
        assert_ne!(first, emit(&cargo_drift).expect("Cargo drift"));
    }
}
