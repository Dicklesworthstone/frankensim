//! Preserve the original material artifact across every structural consumer.
//! Canonical cards retain chemistry, process, claims and model associations;
//! legacy claim packs keep their historical explicit synthetic association.
use super::*;
use fs_matdb::{ClaimSet, NormalizedMaterialCardPack};
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

#[derive(Clone)]
pub enum MaterialPack {
    Claims(NormalizedPack),
    Card(NormalizedMaterialCardPack),
}

impl From<NormalizedPack> for MaterialPack {
    fn from(pack: NormalizedPack) -> Self {
        Self::Claims(pack)
    }
}
impl From<NormalizedMaterialCardPack> for MaterialPack {
    fn from(pack: NormalizedMaterialCardPack) -> Self {
        Self::Card(pack)
    }
}

impl MaterialPack {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        // Dispatch by the wire identity, never the filename. A damaged card
        // cannot silently become an anonymous claim pack.
        if bytes.starts_with(b"FSMCDPK\0") {
            Ok(Self::Card(NormalizedMaterialCardPack::from_bytes(bytes)?))
        } else {
            Ok(Self::Claims(NormalizedPack::from_bytes(bytes)?))
        }
    }
    pub fn card(&self) -> Result<MaterialCard> {
        match self {
            Self::Card(pack) => Ok(pack.card().clone()),
            Self::Claims(pack) => Ok(MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("source-pack:{}", pack.pack_id()),
                    phase: "declared solid model".into(),
                    process: "frozen source-pack realization".into(),
                    revision: 0,
                },
                pack.claims().clone(),
                vec![],
            )?),
        }
    }
    pub fn claims(&self) -> &ClaimSet {
        match self {
            Self::Claims(pack) => pack.claims(),
            Self::Card(pack) => pack.card().claims(),
        }
    }
    pub fn pack_id(&self) -> &str {
        match self {
            Self::Claims(pack) => pack.pack_id(),
            Self::Card(pack) => pack.pack_id(),
        }
    }
    pub fn content_hash(&self) -> fs_blake3::ContentHash {
        match self {
            Self::Claims(pack) => pack.content_hash(),
            Self::Card(pack) => pack.content_hash(),
        }
    }
    pub fn source_artifact(&self) -> fs_blake3::ContentHash {
        match self {
            Self::Claims(pack) => pack.source_artifact(),
            Self::Card(pack) => pack.source_artifact(),
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        match self {
            Self::Claims(pack) => pack.to_bytes(),
            Self::Card(pack) => pack.to_bytes(),
        }
    }
}

/// The performance declares kelvin and radians/second. Preserve compatible
/// source quantity metadata without reinterpreting these numbers as delta-T,
/// hertz, or a different physical axis. Bare historical axes remain bare.
pub fn query(card: &MaterialCard, temperature: f64, omega: f64) -> Result<QueryPoint> {
    for (axis, kind) in [
        ("T", QuantityKind::AbsoluteTemperature),
        (
            "omega",
            QuantityKind::Frequency(FrequencyConvention::Angular),
        ),
    ] {
        let expected = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        for (_, claim) in card.claims().claims_ordered() {
            if let Some(quantity) = claim.validity.axis_quantities().get(axis)
                && *quantity != expected
                && *quantity != QuantitySpec::dimensional(kind.expected_dims())
            {
                return Err(format!("material axis {axis} must retain the performance's kelvin/radian-per-second convention").into());
            }
        }
    }
    // Use the same typed-axis construction as material reeds, walls and plates.
    Ok(fs_couple::material_source::MaterialRecord {
        pack: String::new(),
        state: card.id().process.clone(),
        axes: vec![("T".into(), temperature), ("omega".into(), omega)],
    }
    .query_point(card)?)
}
