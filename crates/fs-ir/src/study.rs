//! Study-form recognition (plan §11.1 / Appendix C): extract the Five
//! Explicits and structure from a `(study "name" clauses...)` form.
//! Extraction only — presence/validity POLICY is the admission bead's
//! (gp3.5); this module gives it spans to point at.

use crate::ast::{Node, NodeKind};
use crate::{IrError, IrErrorKind};

/// A recognized study: borrowed views into the AST.
#[derive(Debug)]
pub struct Study<'a> {
    /// Study name.
    pub name: &'a str,
    /// The `(seed 0x…)` value, if present.
    pub seed: Option<u64>,
    /// The `(versions …)` clause, if present.
    pub versions: Option<&'a Node>,
    /// The `(budget …)` clause, if present.
    pub budget: Option<&'a Node>,
    /// The `(capability …)` clause, if present.
    pub capability: Option<&'a Node>,
    /// `(let name expr)` bindings in order.
    pub lets: Vec<(&'a str, &'a Node)>,
    /// Every remaining body clause in order.
    pub body: Vec<&'a Node>,
}

impl<'a> Study<'a> {
    /// Recognize a study form.
    ///
    /// # Errors
    /// Structured [`IrError`] pointing at the malformed clause.
    pub fn from_node(node: &'a Node) -> Result<Study<'a>, IrError> {
        let items = match node.head() {
            Some("study") => node.items().expect("head implies list"),
            _ => {
                return Err(IrError {
                    span: node.span,
                    kind: IrErrorKind::NotAStudy,
                    detail: "expected a (study \"name\" ...) form".to_string(),
                    hint: "wrap the program in (study \"name\" ...)".to_string(),
                });
            }
        };
        let Some(name_node) = items.get(1) else {
            return Err(IrError {
                span: node.span,
                kind: IrErrorKind::NotAStudy,
                detail: "study has no name".to_string(),
                hint: "the first argument is the study name string".to_string(),
            });
        };
        let NodeKind::Str(name) = &name_node.kind else {
            return Err(IrError {
                span: name_node.span,
                kind: IrErrorKind::NotAStudy,
                detail: "study name must be a string".to_string(),
                hint: "e.g. (study \"spout-laminar-v3\" ...)".to_string(),
            });
        };
        let mut study = Study {
            name,
            seed: None,
            versions: None,
            budget: None,
            capability: None,
            lets: Vec::new(),
            body: Vec::new(),
        };
        let mut seen_seed = false;
        let mut seen_versions = false;
        let mut seen_budget = false;
        let mut seen_capability = false;
        for clause in &items[2..] {
            match clause.head() {
                Some("seed") => {
                    reject_duplicate_pillar("seed", clause, &mut seen_seed)?;
                    let values = clause.items().expect("a clause head implies a list");
                    if values.len() != 2 {
                        return Err(malformed_clause(
                            clause,
                            "seed pillar takes exactly one value",
                            "write (seed 0x...) with no extra operands",
                        ));
                    }
                    study.seed = Some(seed_value(clause).ok_or_else(|| {
                        malformed_clause(
                            clause,
                            "seed pillar does not contain a non-negative u64 seed",
                            "supply one hexadecimal seed such as (seed 0x5EED0001)",
                        )
                    })?);
                }
                Some("versions") => {
                    reject_duplicate_pillar("versions", clause, &mut seen_versions)?;
                    validate_versions_clause(clause)?;
                    study.versions = Some(clause);
                }
                Some("budget") => {
                    reject_duplicate_pillar("budget", clause, &mut seen_budget)?;
                    study.budget = Some(clause);
                }
                Some("capability") => {
                    reject_duplicate_pillar("capability", clause, &mut seen_capability)?;
                    study.capability = Some(clause);
                }
                Some("let") => {
                    if let Some(list) = clause.items()
                        && list.len() == 3
                        && let (Some(sym), Some(expr)) = (list.get(1), list.get(2))
                        && let NodeKind::Symbol(s) = &sym.kind
                    {
                        if study.lets.iter().any(|(name, _)| *name == s.as_str()) {
                            return Err(IrError {
                                span: clause.span,
                                kind: IrErrorKind::MalformedClause,
                                detail: format!("duplicate let binding {s:?} is ambiguous"),
                                hint: "give every let binding a unique name".to_string(),
                            });
                        }
                        study.lets.push((s, expr));
                    } else {
                        return Err(IrError {
                            span: clause.span,
                            kind: IrErrorKind::MalformedClause,
                            detail: "malformed let".to_string(),
                            hint: "(let name expr)".to_string(),
                        });
                    }
                }
                _ => study.body.push(clause),
            }
        }
        Ok(study)
    }

    /// The pinned constellation lock string from
    /// `(versions (constellation :lock "…"))`, if present (the Five
    /// Explicits' versions pillar; pinning must round-trip).
    #[must_use]
    pub fn constellation_lock(&self) -> Option<&'a str> {
        let versions = self.versions?;
        for clause in versions.items()? {
            if clause.head() == Some("constellation")
                && let Some(items) = clause.items()
            {
                for pair in items.windows(2) {
                    if let (NodeKind::Keyword(k), NodeKind::Str(v)) = (&pair[0].kind, &pair[1].kind)
                        && k == "lock"
                    {
                        return Some(v);
                    }
                }
            }
        }
        None
    }
}

fn malformed_clause(clause: &Node, detail: &str, hint: &str) -> IrError {
    IrError {
        span: clause.span,
        kind: IrErrorKind::MalformedClause,
        detail: detail.to_string(),
        hint: hint.to_string(),
    }
}

fn validate_versions_clause(clause: &Node) -> Result<(), IrError> {
    let items = clause.items().expect("a clause head implies a list");
    let mut seen_constellation = false;
    for entry in &items[1..] {
        let Some(entry_items) = entry.items() else {
            return Err(malformed_clause(
                entry,
                "version entries must be parenthesized clauses",
                "write entries such as (constellation :lock \"2026-07\")",
            ));
        };
        if entry.head().is_none() {
            return Err(malformed_clause(
                entry,
                "version entries must have a symbolic name",
                "remove the empty entry or name the versioned component",
            ));
        }
        if entry_items.len() < 2 {
            return Err(malformed_clause(
                entry,
                "version entry has no identity value",
                "supply the version or content identity after the component name",
            ));
        }
        if entry.head() != Some("constellation") {
            continue;
        }
        if seen_constellation {
            return Err(malformed_clause(
                entry,
                "duplicate constellation version entry is ambiguous",
                "retain exactly one (constellation :lock ...) entry",
            ));
        }
        seen_constellation = true;
        let fields = &entry_items[1..];
        if !fields.len().is_multiple_of(2) {
            return Err(malformed_clause(
                entry,
                "constellation version fields must be exact keyword/value pairs",
                "supply one value after every version keyword",
            ));
        }
        let mut seen_lock = false;
        for pair in fields.chunks_exact(2) {
            let NodeKind::Keyword(field) = &pair[0].kind else {
                return Err(malformed_clause(
                    &pair[0],
                    "constellation version field names must be keywords",
                    "write the pin as :lock \"...\"",
                ));
            };
            if field != "lock" {
                continue;
            }
            if seen_lock {
                return Err(malformed_clause(
                    &pair[0],
                    "duplicate constellation :lock field is ambiguous",
                    "retain exactly one :lock field",
                ));
            }
            seen_lock = true;
            if !matches!(&pair[1].kind, NodeKind::Str(value) if !value.trim().is_empty()) {
                return Err(malformed_clause(
                    &pair[1],
                    "constellation :lock must be a non-blank string",
                    "supply the exact constellation lock identity",
                ));
            }
        }
    }
    Ok(())
}

fn reject_duplicate_pillar(name: &str, clause: &Node, seen: &mut bool) -> Result<(), IrError> {
    if *seen {
        return Err(IrError {
            span: clause.span,
            kind: IrErrorKind::MalformedClause,
            detail: format!("duplicate {name} pillar is ambiguous"),
            hint: format!("retain exactly one ({name} ...) clause"),
        });
    }
    *seen = true;
    Ok(())
}

fn seed_value(clause: &Node) -> Option<u64> {
    let items = clause.items()?;
    match items.get(1).map(|n| &n.kind) {
        Some(NodeKind::Seed(v)) => Some(*v),
        Some(NodeKind::Int(i)) if *i >= 0 => Some(u64::try_from(*i).ok()?),
        _ => None,
    }
}
