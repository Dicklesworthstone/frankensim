//! Model cards and the registration lint (patch Rev B mechanism 1): EVERY
//! solver/closure/constitutive law carries a card stating its assumptions,
//! validity domain, known failures, calibration provenance, and promotion
//! level — and a solver WITHOUT a card cannot register. The named
//! model-form risks these exist for: LES closure adequacy, contact-line
//! behavior (vessel flagship), constitutive calibration vs actual
//! construction (frame flagship), surrogate out-of-distribution failure,
//! ground-motion model assumptions.

use crate::{ProvenanceHash, ValidityDomain};
use core::fmt;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Promotion level (the plan's ambition tags, machine-checkable here).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ambition {
    /// `[S]` established mathematics and engineering.
    Solid,
    /// `[F]` research-backed frontier work.
    Frontier,
    /// `[M]` moonshot; must stay feature-flagged until certified.
    Moonshot,
}

impl Ambition {
    /// The plan's tag string.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            Ambition::Solid => "S",
            Ambition::Frontier => "F",
            Ambition::Moonshot => "M",
        }
    }
}

/// One model card: the ledger `model_cards` row.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelCard {
    /// Stable model name (e.g. "les-smagorinsky").
    pub name: String,
    /// Semantic version of the model implementation.
    pub version: String,
    /// Promotion level.
    pub ambition: Ambition,
    /// Stated assumptions (sorted at construction for determinism).
    pub assumptions: Vec<String>,
    /// Parameter box the model is calibrated/valid for.
    pub validity: ValidityDomain,
    /// Known failure modes (teaching text for reports).
    pub known_failures: Vec<String>,
    /// Content-address of the calibration artifact, when one exists.
    pub calibration: Option<ProvenanceHash>,
    /// Headline relative model-form discrepancy band inside the validity
    /// domain (e.g. 0.10 for an LES closure).
    pub discrepancy_rel: f64,
}

impl ModelCard {
    /// Build a card (assumptions/failures are sorted + deduplicated so
    /// every downstream rendering is deterministic).
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        version: impl Into<String>,
        ambition: Ambition,
        mut assumptions: Vec<String>,
        validity: ValidityDomain,
        mut known_failures: Vec<String>,
        discrepancy_rel: f64,
    ) -> Self {
        assumptions.sort_unstable();
        assumptions.dedup();
        known_failures.sort_unstable();
        known_failures.dedup();
        ModelCard {
            name: name.into(),
            version: version.into(),
            ambition,
            assumptions,
            validity,
            known_failures,
            calibration: None,
            discrepancy_rel,
        }
    }

    /// Attach the calibration artifact's provenance.
    #[must_use]
    pub fn with_calibration(mut self, calibration: ProvenanceHash) -> Self {
        self.calibration = Some(calibration);
        self
    }

    /// The ledger `model_cards` row (canonical, deterministic).
    #[must_use]
    pub fn to_ledger_row_json(&self) -> String {
        let list = |xs: &[String]| {
            xs.iter()
                .map(|x| format!("\"{x}\""))
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut s = String::with_capacity(256);
        let _ = write!(
            s,
            "{{\"name\":\"{}\",\"version\":\"{}\",\"ambition\":\"{}\",\"assumptions\":[{}],\
             \"known_failures\":[{}],\"discrepancy_rel\":{},\"calibration\":{},\"validity\":",
            self.name,
            self.version,
            self.ambition.tag(),
            list(&self.assumptions),
            list(&self.known_failures),
            self.discrepancy_rel,
            self.calibration
                .map_or_else(|| "null".to_string(), |c| format!("\"{:016x}\"", c.0)),
        );
        s.push_str(&self.validity_json());
        s.push('}');
        s
    }

    fn validity_json(&self) -> String {
        // Reuse the domain's canonical rendering through a tiny shim (the
        // field is private to keep the invariants in one place).
        let mut s = String::from("{");
        let mut first = true;
        for name in self.validity_params() {
            if let Some((lo, hi)) = self.validity.bound(&name) {
                if !first {
                    s.push(',');
                }
                first = false;
                let _ = write!(s, "\"{name}\":[{lo},{hi}]");
            }
        }
        s.push('}');
        s
    }

    /// The constrained parameter names, sorted (deterministic rendering
    /// and audit order).
    #[must_use]
    pub fn validity_params(&self) -> Vec<String> {
        self.validity.param_names()
    }
}

/// Structured registration failure — the LINT of the acceptance criteria:
/// a solver without a card cannot register.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegistryError {
    /// The named card was never registered.
    NoSuchCard {
        /// The solver attempting to register.
        solver: String,
        /// The card it referenced.
        card: String,
    },
}

impl fmt::Display for RegistryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RegistryError::NoSuchCard { solver, card } => write!(
                f,
                "solver `{solver}` cannot register: model card `{card}` does not exist — every \
                 solver/closure/constitutive law needs a card (assumptions, validity domain, \
                 discrepancy band) BEFORE it can produce Evidence; register the card first"
            ),
        }
    }
}

impl core::error::Error for RegistryError {}

/// The model registry: cards + solver→card bindings. Registration is the
/// enforcement point ("no card, no solver").
#[derive(Debug, Default)]
pub struct ModelRegistry {
    cards: BTreeMap<String, ModelCard>,
    solvers: BTreeMap<String, String>,
}

impl ModelRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        ModelRegistry::default()
    }

    /// Register (or version-bump) a card.
    pub fn register_card(&mut self, card: ModelCard) {
        self.cards.insert(card.name.clone(), card);
    }

    /// Bind a solver to its card. THE lint: refuses when the card is
    /// missing.
    ///
    /// # Errors
    /// [`RegistryError::NoSuchCard`] with teaching text.
    pub fn register_solver(
        &mut self,
        solver: impl Into<String>,
        card: impl Into<String>,
    ) -> Result<(), RegistryError> {
        let (solver, card) = (solver.into(), card.into());
        if !self.cards.contains_key(&card) {
            return Err(RegistryError::NoSuchCard { solver, card });
        }
        self.solvers.insert(solver, card);
        Ok(())
    }

    /// The card a solver is bound to, if registered.
    #[must_use]
    pub fn card_for_solver(&self, solver: &str) -> Option<&ModelCard> {
        self.solvers.get(solver).and_then(|c| self.cards.get(c))
    }

    /// Look a card up by name.
    #[must_use]
    pub fn card(&self, name: &str) -> Option<&ModelCard> {
        self.cards.get(name)
    }

    /// All rows, sorted by card name (the `model_cards` table dump).
    #[must_use]
    pub fn to_ledger_rows_json(&self) -> Vec<String> {
        self.cards
            .values()
            .map(ModelCard::to_ledger_row_json)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn les_card() -> ModelCard {
        ModelCard::new(
            "les-smagorinsky",
            "1.0.0",
            Ambition::Frontier,
            vec![
                "resolved-eddy regime".to_string(),
                "wall model valid for y+ < 100".to_string(),
            ],
            ValidityDomain::unconstrained().with("Re", 1e4, 1e6),
            vec!["under-predicts separation on smooth adverse gradients".to_string()],
            0.10,
        )
    }

    #[test]
    fn solver_without_a_card_cannot_register() {
        let mut reg = ModelRegistry::new();
        let err = reg
            .register_solver("flux.lbm-les", "les-smagorinsky")
            .expect_err("no card yet");
        assert!(err.to_string().contains("register the card first"), "{err}");
        reg.register_card(les_card());
        reg.register_solver("flux.lbm-les", "les-smagorinsky")
            .expect("card exists now");
        let card = reg.card_for_solver("flux.lbm-les").expect("bound");
        assert_eq!(card.ambition.tag(), "F");
        assert!((card.discrepancy_rel - 0.10).abs() < 1e-12);
    }

    #[test]
    fn card_rows_are_canonical_and_sorted() {
        let card = les_card().with_calibration(ProvenanceHash::of_bytes(b"cal"));
        let row = card.to_ledger_row_json();
        assert_eq!(row, card.to_ledger_row_json(), "deterministic");
        for key in [
            "\"name\":\"les-smagorinsky\"",
            "\"ambition\":\"F\"",
            "\"discrepancy_rel\":0.1",
            "\"validity\":{\"Re\":[10000,1000000]}",
        ] {
            assert!(row.contains(key), "missing {key} in {row}");
        }
        // Assumptions are sorted regardless of construction order.
        let a = row.find("resolved-eddy").expect("assumption present");
        let b = row.find("wall model").expect("assumption present");
        assert!(a < b);
    }
}
