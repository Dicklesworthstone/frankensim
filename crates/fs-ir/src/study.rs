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
        for clause in &items[2..] {
            match clause.head() {
                Some("seed") => study.seed = seed_value(clause),
                Some("versions") => study.versions = Some(clause),
                Some("budget") => study.budget = Some(clause),
                Some("capability") => study.capability = Some(clause),
                Some("let") => {
                    if let Some(list) = clause.items()
                        && let (Some(sym), Some(expr)) = (list.get(1), list.get(2))
                        && let NodeKind::Symbol(s) = &sym.kind
                    {
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
                    if let (NodeKind::Keyword(k), NodeKind::Str(v)) =
                        (&pair[0].kind, &pair[1].kind)
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

fn seed_value(clause: &Node) -> Option<u64> {
    let items = clause.items()?;
    match items.get(1).map(|n| &n.kind) {
        Some(NodeKind::Seed(v)) => Some(*v),
        Some(NodeKind::Int(i)) if *i >= 0 => Some(u64::try_from(*i).ok()?),
        _ => None,
    }
}
