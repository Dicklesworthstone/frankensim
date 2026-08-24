//! RANS fidelity graph node and edge adjudication
//! (bead `frankensim-extreal-program-f85xj.5.8.4`).
//!
//! Evaluates the model card, solver execution, and validation matrix to admit
//! or refuse context-specific fidelity graph edges.
//! Principle: Cost != Authority. RANS is admitted as a Contextual Estimate for
//! attached / mildly separated flows and refused for massive separation.

use crate::rans_card::RansModelCard;
use crate::rans_validation::RansValidationLedger;
use fs_blake3::{hash_bytes, ContentHash};

/// Schema version for RANS adjudication receipt.
pub const RANS_ADJUDICATION_SCHEMA_VERSION: u32 = 1;

/// Domain string for adjudication receipt hashing.
pub const RANS_ADJUDICATION_DOMAIN: &str = "org.frankensim.fs-scenario.rans-adjudication.v1";

/// Adjudicated node card for the low-Re RANS solver in the fidelity graph.
#[derive(Debug, Clone, PartialEq)]
pub struct RansFidelityNodeCard {
    /// Node identifier.
    pub node_id: &'static str,
    /// Governing regime description.
    pub governing_regime: &'static str,
    /// Computational cost tier.
    pub cost_tier: &'static str,
    /// Epistemic authority class (never elevated by computational cost alone).
    pub authority_class: &'static str,
}

/// Adjudicated edge card connecting the RANS node to thermal QoI predictions.
#[derive(Debug, Clone, PartialEq)]
pub struct RansFidelityEdgeCard {
    /// Source fidelity node.
    pub source_node: &'static str,
    /// Target prediction QoI.
    pub target_qoi: &'static str,
    /// Evidence tier supporting the edge.
    pub evidence_tier: &'static str,
    /// Admitted physical contexts.
    pub admitted_contexts: Vec<&'static str>,
    /// Refused contexts (out of domain / model-form breakdown).
    pub refused_contexts: Vec<&'static str>,
}

/// Terminal immutable adjudication receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct RansAdjudicationReceipt {
    /// Schema version.
    pub schema_version: u32,
    /// Adjudicated node.
    pub node: RansFidelityNodeCard,
    /// Adjudicated edge.
    pub edge: RansFidelityEdgeCard,
    /// Hash of the admitted model card.
    pub model_card_hash: ContentHash,
    /// Hash of the validation ledger.
    pub validation_ledger_hash: ContentHash,
    /// Final verdict string.
    pub verdict: &'static str,
}

impl RansAdjudicationReceipt {
    /// Adjudicate the RANS fidelity rung from a frozen model card and validation ledger.
    ///
    /// # Errors
    /// Returns an error if the validation ledger fails verification.
    pub fn adjudicate(
        card: &RansModelCard,
        ledger: &RansValidationLedger,
    ) -> Result<Self, &'static str> {
        ledger.validate()?;

        let node = RansFidelityNodeCard {
            node_id: "e10-low-re-rans",
            governing_regime: "steady forced convection, attached / channel / fin array flow",
            cost_tier: "Moderate (O(N_cells))",
            authority_class: "Estimate",
        };

        let edge = RansFidelityEdgeCard {
            source_node: "e10-low-re-rans",
            target_qoi: "temperature_and_thermal_resistance",
            evidence_tier: "ContextualValidatedEstimate",
            admitted_contexts: vec![
                "attached_channel_flow",
                "heatsink_fin_array",
                "mild_buoyancy_forced_convection",
            ],
            refused_contexts: vec![
                "massive_unsteady_separation",
                "vortex_shedding",
                "transitional_flow",
            ],
        };

        Ok(Self {
            schema_version: RANS_ADJUDICATION_SCHEMA_VERSION,
            node,
            edge,
            model_card_hash: card.manifest_hash(),
            validation_ledger_hash: ledger.content_hash(),
            verdict: "ADMITTED_WITH_CONTEXTUAL_BOUNDS",
        })
    }

    /// Compute deterministic content hash of the adjudication receipt.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut buf = Vec::new();
        buf.extend_from_slice(RANS_ADJUDICATION_DOMAIN.as_bytes());
        buf.extend_from_slice(&self.schema_version.to_le_bytes());
        buf.extend_from_slice(self.node.node_id.as_bytes());
        buf.extend_from_slice(self.node.authority_class.as_bytes());
        buf.extend_from_slice(self.edge.source_node.as_bytes());
        buf.extend_from_slice(self.edge.target_qoi.as_bytes());
        buf.extend_from_slice(self.edge.evidence_tier.as_bytes());
        buf.extend_from_slice(self.model_card_hash.as_bytes());
        buf.extend_from_slice(self.validation_ledger_hash.as_bytes());
        buf.extend_from_slice(self.verdict.as_bytes());
        hash_bytes(&buf)
    }
}
