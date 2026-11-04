//! Resolve one reported Bond ball-mill test without manufacturing ore defaults.
//!
//! Four scalar claims must identify the same original observation: work index,
//! closing screen, feed F80 and product P80. The selected claims, evidence,
//! uncertainty, receipts, observation and typed query remain in the original
//! `BondGrindingCard`, so blends, inversion and circuit use retain their source.
//! Original metric-tonne/short-ton property names prevent a denominator swap.
//!
//! A reported test point is not a fitted operational validity box. The caller
//! supplies the model's `BondSizeDomain` explicitly. Resolving the observation
//! does not establish transfer to another ore, grind, preparation or circuit,
//! nor does it turn a calculated energy or throughput into a measurement.

use fs_blake3::ContentHash;
use fs_matdb::{
    MatDbError, MaterialAnswer, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy,
};
use fs_qty::Length;

use super::{
    BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex, ComminutionError,
    JOULES_PER_KWH, KG_PER_SHORT_TON, SpecificGrindingEnergy, WorkIndexMassBasis,
};

/// Exact property names preserve the work-index test and original mass convention.
#[must_use]
pub fn work_index_key(basis: WorkIndexMassBasis) -> PropertyKey {
    PropertyKey::new(
        match basis {
            WorkIndexMassBasis::MetricTonne => "bond-ball-mill-work-index-metric-tonne",
            WorkIndexMassBasis::ShortTon => "bond-ball-mill-work-index-short-ton",
        },
        SpecificGrindingEnergy::DIMS,
    )
}

/// Immutable material source retained by a resolved grinding card.
#[derive(Debug, Clone, PartialEq)]
pub struct BondCalibrationSource {
    pack_id: String,
    pack_hash: ContentHash,
    query: QueryPoint,
    observation: ObservationDataset,
    claims: [PropertyClaim; 4],
    answers: [MaterialAnswer; 4],
}

impl BondCalibrationSource {
    /// Identifier of the original normalized material pack.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }

    /// Content identity of the complete original source pack.
    #[must_use]
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }

    /// Source-selection conditions, not the model's operational size window.
    #[must_use]
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }

    /// Original specimen, test method, caveats, artifact and source terms.
    #[must_use]
    pub const fn observation(&self) -> &ObservationDataset {
        &self.observation
    }

    /// Original claims in work-index, closing-screen, F80, P80 order.
    #[must_use]
    pub const fn claims(&self) -> &[PropertyClaim; 4] {
        &self.claims
    }

    /// Original evaluated evidence and verified usage receipts, in claim order.
    #[must_use]
    pub const fn answers(&self) -> &[MaterialAnswer; 4] {
        &self.answers
    }

    /// Reported laboratory F80 and P80 in metres; these are not a validity range.
    #[must_use]
    pub fn test_sizes(&self) -> [Length; 2] {
        [
            Length::new(self.answers[2].evidence.value.value),
            Length::new(self.answers[3].evidence.value.value),
        ]
    }
}

/// Original material/model refusal or an incompatible test row.
#[derive(Debug, Clone, PartialEq)]
pub enum BondSourceError {
    /// Original typed query or receipt-verification failure.
    Material(MatDbError),
    /// Original physical Bond-card admission failure.
    Grinding(ComminutionError),
    /// Missing or incompatible original test information.
    Invalid(&'static str),
}

impl core::fmt::Display for BondSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Bond test source: {self:?}")
    }
}
impl std::error::Error for BondSourceError {}
impl From<MatDbError> for BondSourceError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}
impl From<ComminutionError> for BondSourceError {
    fn from(value: ComminutionError) -> Self {
        Self::Grinding(value)
    }
}

fn require(condition: bool, message: &'static str) -> Result<(), BondSourceError> {
    if condition {
        Ok(())
    } else {
        Err(BondSourceError::Invalid(message))
    }
}

/// Resolve a complete scalar test row and build the existing empirical Bond card.
///
/// Every source field must have exactly one common original observation and the
/// same source validity. Citation-only values, interpolation, mixed specimens,
/// unsupported query coordinates and unreported values refuse. The work-index
/// property is selected by its original mass convention, while its stored value
/// remains coherent SI J/kg. Source labels/receipts establish provenance, not
/// independent authentication of laboratory results or plant qualification.
///
/// `model_domain` is an explicit caller declaration. It is not widened to include
/// the source's test sizes, and the test sizes are not widened into an ore model.
pub fn resolve_bond_ball_mill(
    pack: &NormalizedPack,
    point: &QueryPoint,
    basis: WorkIndexMassBasis,
    policy: SelectionPolicy,
    model_domain: BondSizeDomain,
) -> Result<BondGrindingCard, BondSourceError> {
    let keys = [
        work_index_key(basis),
        PropertyKey::new("bond-test-closing-screen", Length::DIMS),
        PropertyKey::new("bond-test-feed-f80", Length::DIMS),
        PropertyKey::new("bond-test-product-p80", Length::DIMS),
    ];
    let resolve = |key: &PropertyKey| -> Result<_, BondSourceError> {
        let answer = pack.claims().query_typed(key, point, policy)?;
        pack.claims().verify_receipt(&answer.receipt)?;
        let claim = pack
            .claims()
            .claim(answer.receipt.selected)
            .ok_or(BondSourceError::Invalid("selected source claim is absent"))?;
        require(
            matches!(claim.value, PropertyValue::Scalar { .. }) && claim.observations.len() == 1,
            "each test field requires a scalar and one original specimen observation",
        )?;
        require(
            point
                .axes()
                .keys()
                .all(|axis| claim.validity.bounds().contains_key(axis)),
            "query includes conditions unsupported by the original test",
        )?;
        let value = answer.evidence.value.value;
        require(
            value.is_finite() && value > 0.0,
            "reported work index and sizes must be finite and positive",
        )?;
        Ok((claim.clone(), answer))
    };
    let [(work, wi), (screen, closing), (feed, f80), (product, p80)] = [
        resolve(&keys[0])?,
        resolve(&keys[1])?,
        resolve(&keys[2])?,
        resolve(&keys[3])?,
    ];
    for claim in [&screen, &feed, &product] {
        require(
            claim.observations == work.observations && claim.validity == work.validity,
            "work index and sizes must retain the same specimen observation and source support",
        )?;
    }
    let observation = pack
        .claims()
        .observation(work.observations[0])
        .ok_or(BondSourceError::Invalid(
            "original test observation is absent",
        ))?
        .clone();
    require(
        f80.evidence.value.value > p80.evidence.value.value
            && p80.evidence.value.value <= closing.evidence.value.value,
        "reported test must reduce F80 to P80 no larger than the closing screen",
    )?;

    // The original normalized energy is retained exactly. Reconstruct only the
    // displayed kWh/source-mass value; do not round-trip it back into the kernel.
    let kg = match basis {
        WorkIndexMassBasis::MetricTonne => 1000.0,
        WorkIndexMassBasis::ShortTon => KG_PER_SHORT_TON,
    };
    let source_value = wi.evidence.value.value / (JOULES_PER_KWH / kg);
    require(
        source_value.is_finite() && source_value > 0.0,
        "work-index source-unit reconstruction is unrepresentable",
    )?;
    let mut card = BondGrindingCard::new(
        BondWorkIndex {
            source_value,
            source_basis: basis,
            energy: SpecificGrindingEnergy::new(wi.evidence.value.value),
        },
        BondTestContext {
            specimen: observation.specimen.clone(),
            method: observation.method.clone(),
            closing_screen: Length::new(closing.evidence.value.value),
            provenance: observation.provenance.clone(),
        },
        model_domain,
    )?;
    card.source = Some(Box::new(BondCalibrationSource {
        pack_id: pack.pack_id().to_owned(),
        pack_hash: pack.content_hash(),
        query: point.clone(),
        observation,
        claims: [work, screen, feed, product],
        answers: [wi, closing, f80, p80],
    }));
    Ok(card)
}

#[cfg(test)]
#[path = "source/tests.rs"]
mod tests;
