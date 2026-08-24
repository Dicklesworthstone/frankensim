//! Workspace solver-state charter collision check (bead
//! frankensim-sj31i.52.5.2.1, close-item "workspace-level xtask
//! collision check").
//!
//! Every `SolverStateV2::charter()` implementation declares a
//! `StateIdentityCharterV2` literal somewhere in workspace sources.
//! The charter registry refuses duplicate (owner, family) pairs with
//! drifted grammar at RUNTIME; this lane enforces the same contract at
//! REPO level so two lanes cannot even land colliding declarations.
//!
//! Scope honesty: this is a deterministic SOURCE scan, not a Rust
//! type-system proof. It resolves `..Base` struct-update syntax only
//! against same-file `const`/`static` charter bindings; anything it
//! cannot resolve is a violation (fail-closed), never a silent pass.
//!
//! Checks:
//! 1. duplicate (owner, state_family) with differing
//!    (schema_grammar, codec_grammar, codec_version) — refuse;
//! 2. struct-update bases that do not resolve to a same-file charter
//!    const/static binding — refuse;
//! 3. malformed or incomplete charter literals — refuse.

use std::collections::BTreeMap;
use std::path::Path;

/// One parsed charter declaration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CharterDeclaration {
    /// Workspace-relative source path.
    pub path: String,
    /// 1-based line of the charter value's opening brace.
    pub line: usize,
    /// Owning crate/module path.
    pub owner: String,
    /// State-family name.
    pub state_family: String,
    /// Schema grammar string.
    pub schema_grammar: String,
    /// Codec grammar string.
    pub codec_grammar: String,
    /// Codec version.
    pub codec_version: u32,
}

/// Extract a string-typed field value from one literal's text.
fn string_field(literal: &str, field: &str) -> Option<String> {
    let marker = format!("{field}:");
    let field_offset = literal.find(&marker)?;
    let rest = literal[field_offset + marker.len()..].trim_start();
    // A leading `&` (string-literal reference) is permitted but never
    // required; `?` here would make it mandatory and reject every plain
    // `owner: "..."` field as missing.
    let rest = rest.strip_prefix('&').unwrap_or(rest).trim_start();
    let quoted = rest.strip_prefix('"')?;
    let end = quoted.find('"')?;
    Some(quoted[..end].to_string())
}

/// Extract a u32-typed field value from one literal's text.
fn uint_field(literal: &str, field: &str) -> Option<u32> {
    let marker = format!("{field}:");
    let field_offset = literal.find(&marker)?;
    let rest = literal[field_offset + marker.len()..].trim_start();
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        None
    } else {
        digits.parse().ok()
    }
}

/// Parse every charter CONSTRUCTOR literal in one source file's text.
/// Constructor = the `StateIdentityCharterV2 {` token form (the type-
/// annotation occurrence is followed by `=`, never `{`, so requiring
/// the brace disambiguates without a full tokenizer).
///
/// `..Base` updates resolve ONLY against same-file `const`/`static`
/// bindings captured earlier in the file (deterministic order).
pub fn scan_source(path: &str, text: &str) -> Result<Vec<CharterDeclaration>, String> {
    let mut declarations = Vec::new();
    let mut base_bindings: BTreeMap<String, CharterDeclaration> = BTreeMap::new();
    let bytes = text.as_bytes();
    let mut cursor = 0usize;
    while let Some(offset) = text[cursor..].find("StateIdentityCharterV2") {
        let token_start = cursor + offset;
        let after_token = token_start + "StateIdentityCharterV2".len();
        let rest_trim = text[after_token..].trim_start();
        // Declarations (`pub struct X {`, `impl X {`) also place `{`
        // directly after the type name; they are not constructor
        // literals, and treating them as such misreads their member
        // lists as charter fields (fs-exec defines the type beside its
        // constructors). Skip when the same-line prefix ends in a
        // declaration keyword.
        let line_start = text[..token_start].rfind('\n').map_or(0, |pos| pos + 1);
        let same_line_prefix = text[line_start..token_start].trim_end();
        if same_line_prefix.ends_with("struct") || same_line_prefix.ends_with("impl")
        {
            cursor = after_token;
            continue;
        }
        // Constructor only: next non-space character opens the literal.
        if !rest_trim.starts_with('{') {
            cursor = after_token;
            continue;
        }
        let open = after_token + (text[after_token..].len() - rest_trim.len());
        let line = text[..open].bytes().filter(|b| *b == b'\n').count() + 1;
        // Brace-match the literal body.
        let mut depth = 0usize;
        let mut close = None;
        for (index, byte) in bytes[open..].iter().enumerate() {
            match byte {
                b'{' => depth += 1,
                b'}' => {
                    depth -= 1;
                    if depth == 0 {
                        close = Some(open + index);
                        break;
                    }
                }
                _ => {}
            }
        }
        let Some(close) = close else {
            return Err(format!("{path}:{line}: unterminated charter literal"));
        };
        let literal = &text[open..=close];
        // A `{ ..Base }` update places the `..` after the opening
        // brace (and possibly no commas at all), so strip the brace
        // before scanning fields and take the base name as its leading
        // identifier rather than trusting surrounding punctuation.
        let body = literal.trim_start().strip_prefix('{').unwrap_or(literal);
        let update_base = body.split(',').find_map(|field| {
            let field = field.trim().strip_prefix("..")?;
            let name: String = field
                .chars()
                .take_while(|c| c.is_alphanumeric() || *c == '_')
                .collect();
            (!name.is_empty()).then_some(name)
        });

        let declaration = if let Some(base) = update_base {
            let Some(resolved) = base_bindings.get(base.as_str()) else {
                return Err(format!(
                    "{path}:{line}: charter struct-update base {base:?} does not \
                     resolve to a same-file const/static charter binding"
                ));
            };
            let mut inherited = resolved.clone();
            inherited.path = path.to_string();
            inherited.line = line;
            if let Some(value) = string_field(literal, "owner") {
                inherited.owner = value;
            }
            if let Some(value) = string_field(literal, "state_family") {
                inherited.state_family = value;
            }
            if let Some(value) = string_field(literal, "schema_grammar") {
                inherited.schema_grammar = value;
            }
            if let Some(value) = string_field(literal, "codec_grammar") {
                inherited.codec_grammar = value;
            }
            if let Some(value) = uint_field(literal, "codec_version") {
                inherited.codec_version = value;
            }
            inherited
        } else {
            let missing = [
                ("owner", string_field(literal, "owner")),
                ("state_family", string_field(literal, "state_family")),
                ("schema_grammar", string_field(literal, "schema_grammar")),
                ("codec_grammar", string_field(literal, "codec_grammar")),
            ]
            .into_iter()
            .find_map(|(name, value)| value.is_none().then_some(name));
            if let Some(name) = missing {
                return Err(format!("{path}:{line}: charter literal lacks {name}"));
            }
            let Some(codec_version) = uint_field(literal, "codec_version") else {
                return Err(format!(
                    "{path}:{line}: charter literal lacks codec_version"
                ));
            };
            CharterDeclaration {
                path: path.to_string(),
                line,
                owner: string_field(literal, "owner")
                    .expect("checked above"),
                state_family: string_field(literal, "state_family")
                    .expect("checked above"),
                schema_grammar: string_field(literal, "schema_grammar")
                    .expect("checked above"),
                codec_grammar: string_field(literal, "codec_grammar")
                    .expect("checked above"),
                codec_version,
            }
        };

        // Capture `const NAME` / `static NAME` bindings (the binding
        // keyword precedes the name, which precedes the type and `=`).
        let prefix = text[..token_start].trim_end();
        if prefix.ends_with('}') || prefix.ends_with(';') || prefix.is_empty() {
            // Not a binding position? Still probe: bindings read
            // `const NAME: <type> = <value>`, so search backwards for
            // the nearest `const `/`static ` keyword.
        }
        if let Some(keyword_position) = prefix.rfind("const ").or_else(|| prefix.rfind("static ")) {
            let between = &prefix[keyword_position..];
            // Only a DIRECT binding (keyword then name then colon)
            // counts; a distant keyword belongs to another item.
            if let Some(colon) = between.find(':') {
                let name_part = between[..colon].trim();
                let mut name_tokens = name_part.split_whitespace();
                let _kind = name_tokens.next();
                if let Some(name) = name_tokens.next() {
                    if name.chars().all(|c| c.is_alphanumeric() || c == '_') {
                        base_bindings.insert(name.to_string(), declaration.clone());
                    }
                }
            }
        }
        declarations.push(declaration);
        cursor = close + 1;
    }
    Ok(declarations)
}

/// Deterministic collision check over the parsed declarations.
#[must_use]
pub fn check_collisions(declarations: &[CharterDeclaration]) -> Vec<String> {
    let mut violations = Vec::new();
    let mut by_key: BTreeMap<(&str, &str), &CharterDeclaration> = BTreeMap::new();
    for declaration in declarations {
        let key = (declaration.owner.as_str(), declaration.state_family.as_str());
        match by_key.get(&key) {
            Some(existing) => {
                let same_grammar = existing.schema_grammar == declaration.schema_grammar
                    && existing.codec_grammar == declaration.codec_grammar
                    && existing.codec_version == declaration.codec_version;
                if !same_grammar {
                    violations.push(format!(
                        "charter collision: {}:{} and {}:{} declare owner {:?} family {:?} \
                         with drifted grammar (schema {:?} vs {:?}, codec {:?} vs {:?}, \
                         version {} vs {}); rotate the intended ids or pick a new family",
                        existing.path,
                        existing.line,
                        declaration.path,
                        declaration.line,
                        declaration.owner,
                        declaration.state_family,
                        existing.schema_grammar,
                        declaration.schema_grammar,
                        existing.codec_grammar,
                        declaration.codec_grammar,
                        existing.codec_version,
                        declaration.codec_version,
                    ));
                }
            }
            None => {
                by_key.insert(key, declaration);
            }
        }
    }
    violations
}

/// Remove `#[cfg(test)] mod ... { ... }` blocks (brace-matched) from one
/// source text. The workspace gate audits the production identity
/// surface; test modules deliberately construct colliding and drifted
/// charters (impostor owners, rotated codec versions) to exercise the
/// collision registry, so scanning them reports fixture data as
/// violations. Stripping is deterministic: occurrences are processed
/// left to right, each removing its whole brace-matched module body.
fn strip_test_modules(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut cursor = 0usize;
    while let Some(attribute) = text[cursor..].find("#[cfg(test)]") {
        let attribute_start = cursor + attribute;
        let after_attribute = attribute_start + "#[cfg(test)]".len();
        let rest = text[after_attribute..].trim_start();
        // Only a following `mod` opens a removable block; other uses of
        // the attribute (e.g. on individual test helpers) stay.
        if !rest.starts_with("mod") {
            out.push_str(&text[cursor..after_attribute]);
            cursor = after_attribute;
            continue;
        }
        out.push_str(&text[cursor..attribute_start]);
        // rest starts after trimmed whitespace; convert its internal
        // offsets with the same length-difference form scan_source uses.
        let body_start = after_attribute + (text[after_attribute..].len() - rest.len());
        let brace = body_start + rest.find('{').unwrap_or(0);
        if brace >= text.len() || text.as_bytes()[brace] != b'{' {
            cursor = after_attribute;
            continue;
        }
        let mut depth = 0usize;
        let mut close = None;
        for (index, byte) in bytes[brace..].iter().enumerate() {
            match byte {
                b'{' => depth += 1,
                b'}' => {
                    depth -= 1;
                    if depth == 0 {
                        close = Some(brace + index);
                        break;
                    }
                }
                _ => {}
            }
        }
        match close {
            Some(close) => cursor = close + 1,
            None => {
                // Unterminated block: keep everything from here verbatim
                // rather than guessing; scan_source still brace-checks it.
                out.push_str(&text[attribute_start..]);
                return out;
            }
        }
    }
    out.push_str(&text[cursor..]);
    out
}

/// Walk the workspace `crates/` tree and run the full lane over the
/// production identity surface (test modules stripped). Returns the
/// sorted violation list (empty = green).
pub fn run(root: &Path) -> Result<Vec<String>, String> {
    let mut declarations = Vec::new();
    let crates_dir = root.join("crates");
    let mut stack = vec![crates_dir];
    let mut sources = Vec::new();
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir)
            .map_err(|error| format!("cannot read {}: {error}", dir.display()))?;
        for entry in entries.flatten() {
            let path = entry.path();
            if entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
                if matches!(
                    entry.file_name().to_str(),
                    Some("target" | "node_modules" | ".git" | ".beads")
                ) {
                    continue;
                }
                stack.push(path);
            } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
                sources.push(path);
            }
        }
    }
    sources.sort();
    for source in sources {
        let relative = source
            .strip_prefix(root)
            .unwrap_or(&source)
            .display()
            .to_string();
        let text = std::fs::read_to_string(&source)
            .map_err(|error| format!("cannot read {relative}: {error}"))?;
        if text.contains("StateIdentityCharterV2") {
            let production = strip_test_modules(&text);
            declarations.extend(scan_source(&relative, &production)?);
        }
    }
    let mut violations = check_collisions(&declarations);
    violations.sort();
    Ok(violations)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FULL_LITERAL: &str = concat!(
        "const ALPHA_CHARTER: snapshot_v2::StateIdentityCharterV2 = snapshot_v2::StateIdentityCharterV2 {\n",
        "            owner: \"fs-exec::solver\",\n",
        "            state_family: \"jacobi-iteration\",\n",
        "            schema_grammar: \"iter: u64 le; x: f64 le[]\",\n",
        "            codec_grammar: \"v2 framed: u64 count, f64 le slice\",\n",
        "            codec_version: 1,\n",
        "        };"
    );

    #[test]
    fn full_literals_parse_and_bind_for_updates() {
        let source = format!(
            "{FULL_LITERAL}\n let derived = StateIdentityCharterV2 {{ ..ALPHA_CHARTER }};\n"
        );
        let declarations = scan_source("crates/x/src/lib.rs", &source)
            .expect("both literals parse");
        assert_eq!(declarations.len(), 2);
        assert_eq!(declarations[0].owner, "fs-exec::solver");
        assert_eq!(declarations[0].state_family, "jacobi-iteration");
        // The update inherits every identity field but records its own
        // site, so compare the identity fields rather than whole struct
        // (path/line are per-occurrence by design).
        let (base, derived) = (&declarations[0], &declarations[1]);
        assert_eq!(derived.owner, base.owner);
        assert_eq!(derived.state_family, base.state_family);
        assert_eq!(derived.schema_grammar, base.schema_grammar);
        assert_eq!(derived.codec_grammar, base.codec_grammar);
        assert_eq!(derived.codec_version, base.codec_version);
        assert_ne!(derived.line, base.line, "the update records its own site");
        assert!(check_collisions(&declarations).is_empty());
    }

    #[test]
    fn duplicate_pair_with_drifted_grammar_refuses_deterministically() {
        let drifted = FULL_LITERAL
            .replace("codec_version: 1,", "codec_version: 2,")
            .replace(
                "const ALPHA_CHARTER:",
                "const BETA_CHARTER:",
            )
            .replace(
                "snapshot_v2::StateIdentityCharterV2 = snapshot_v2::StateIdentityCharterV2 {",
                "snapshot_v2::StateIdentityCharterV2 = snapshot_v2::StateIdentityCharterV2 {",
            )
            .replacen("StateIdentityCharterV2 {", "StateIdentityCharterV2 {", 1);
        let source =
            format!("{FULL_LITERAL}\n{drifted}\n");
        let declarations = scan_source("crates/x/src/lib.rs", &source).expect("parses");
        assert_eq!(declarations.len(), 2);
        let violations = check_collisions(&declarations);
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert!(violations[0].contains("charter collision"), "{violations:?}");
        assert!(violations[0].contains("jacobi-iteration"), "{violations:?}");
    }

    #[test]
    fn identical_duplicates_are_registry_idempotent() {
        let source = format!("{FULL_LITERAL}\n{FULL_LITERAL}");
        let declarations = scan_source("crates/x/src/lib.rs", &source).expect("parses");
        assert!(check_collisions(&declarations).is_empty());
    }

    #[test]
    fn unresolvable_update_base_fails_closed() {
        let source =
            "let d = StateIdentityCharterV2 { ..NO_SUCH_BASE };".to_string();
        let error = scan_source("crates/x/src/lib.rs", &source)
            .expect_err("unresolvable base refuses");
        assert!(error.contains("NO_SUCH_BASE"), "{error}");
    }

    #[test]
    fn incomplete_literal_names_the_missing_field() {
        let source = "const C: StateIdentityCharterV2 = StateIdentityCharterV2 { owner: \"a\" };"
            .to_string();
        let error = scan_source("crates/x/src/lib.rs", &source)
            .expect_err("missing fields refuse");
        assert!(error.contains("lacks"), "{error}");
    }

    #[test]
    fn type_declarations_are_not_constructor_literals() {
        // `pub struct StateIdentityCharterV2 {` and
        // `impl StateIdentityCharterV2 {` place `{` after the type name
        // exactly like constructors; both must be skipped so the
        // definition beside its constructors in fs-exec/src/solver.rs
        // does not read as a field-less charter literal.
        let source = concat!(
            "pub struct StateIdentityCharterV2 {\n",
            "    pub owner: &'static str,\n",
            "}\n",
            "impl StateIdentityCharterV2 {\n",
            "    pub fn owner(&self) -> &str { self.owner }\n",
            "}\n",
            "let real = StateIdentityCharterV2 {\n",
            "    owner: \"fs-exec::solver\",\n",
            "    state_family: \"jacobi-iteration\",\n",
            "    schema_grammar: \"g\",\n",
            "    codec_grammar: \"c\",\n",
            "    codec_version: 1,\n",
            "};\n",
        );
        let declarations =
            scan_source("crates/x/src/lib.rs", source).expect("parses");
        assert_eq!(declarations.len(), 1, "{declarations:?}");
        assert_eq!(declarations[0].owner, "fs-exec::solver");
    }

    #[test]
    fn type_annotation_occurrence_is_not_a_constructor() {
        // The annotation `: StateIdentityCharterV2 =` must not produce a
        // phantom declaration or consume the real one.
        let declarations =
            scan_source("crates/x/src/lib.rs", FULL_LITERAL).expect("parses");
        assert_eq!(declarations.len(), 1);
        assert_eq!(declarations[0].line, 1);
    }

    #[test]
    fn test_modules_are_stripped_before_scanning() {
        // Deliberately colliding fixture charters live inside
        // `#[cfg(test)] mod` blocks; the workspace gate audits the
        // production surface, so the stripper removes those blocks whole.
        let source = concat!(
            "const REAL: StateIdentityCharterV2 = StateIdentityCharterV2 {\n",
            "    owner: \"c::real\",\n",
            "    state_family: \"fam\",\n",
            "    schema_grammar: \"g\",\n",
            "    codec_grammar: \"k\",\n",
            "    codec_version: 1,\n",
            "};\n",
            "#[cfg(test)]\n",
            "mod fixtures {\n",
            "    const IMPOSTOR_A: StateIdentityCharterV2 = StateIdentityCharterV2 {\n",
            "        owner: \"t::impostor\",\n",
            "        state_family: \"fam\",\n",
            "        schema_grammar: \"drifted\",\n",
            "        codec_grammar: \"k\",\n",
            "        codec_version: 2,\n",
            "    };\n",
            "}\n",
        );
        let stripped = strip_test_modules(source);
        assert!(!stripped.contains("IMPOSTOR_A"), "{stripped}");
        let declarations =
            scan_source("crates/x/src/lib.rs", &stripped).expect("parses");
        assert_eq!(declarations.len(), 1);
        assert_eq!(declarations[0].owner, "c::real");
    }
}
