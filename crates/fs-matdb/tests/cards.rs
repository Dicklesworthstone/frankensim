//! fs-matdb PR-2 conformance: material/constitutive cards, revision
//! lineage, and content identity. Cards are immutable — supersession
//! creates a successor that links its predecessor's hash; nothing is
//! ever edited in place.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MATDB_SCHEMA_VERSION, MatDbError, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::Dims;

const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const STRESS_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn provenance() -> Provenance {
    Provenance {
        source: "calibration report LAB-42".to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn genesis_id() -> MaterialStateId {
    MaterialStateId {
        chemistry: "AA6061".to_string(),
        phase: "wrought".to_string(),
        process: "T6".to_string(),
        revision: 0,
    }
}

fn density_claims() -> ClaimSet {
    let mut set = ClaimSet::new();
    set.insert_claim(PropertyClaim {
        key: PropertyKey::new("density", DENSITY_DIMS),
        value: PropertyValue::Scalar {
            value: 2700.0,
            dims: DENSITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.01,
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(),
    })
    .expect("density claim inserts");
    set
}

fn j2_card() -> ConstitutiveModelCard {
    let mut parameters = BTreeMap::new();
    parameters.insert(
        "yield_stress".to_string(),
        LawParameter {
            value: 276.0e6,
            dims: STRESS_DIMS,
        },
    );
    parameters.insert(
        "hardening_modulus".to_string(),
        LawParameter {
            value: 1.2e9,
            dims: STRESS_DIMS,
        },
    );
    ConstitutiveModelCard {
        law: LawId("j2-plasticity-voce".to_string()),
        law_version: 1,
        parameters,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
        sources: vec![hash_bytes(b"tensile calibration v1")],
        provenance: provenance(),
    }
}

fn parameter_hash(card: &ConstitutiveModelCard) -> ContentHash {
    card.canonical_parameters_hash()
        .expect("valid card mints canonical parameter identity")
}

#[test]
fn genesis_card_assembles_and_indexes_claims_and_models() {
    let card = MaterialCard::assemble(genesis_id(), density_claims(), vec![j2_card()])
        .expect("genesis assembles");
    assert_eq!(card.schema_version(), MATDB_SCHEMA_VERSION);
    assert_eq!(card.supersedes(), None);
    assert_eq!(card.claims_for("density").len(), 1);
    assert_eq!(
        card.models_for(&LawId("j2-plasticity-voce".to_string()))
            .len(),
        1
    );
    assert!(
        card.models_for(&LawId("neo-hookean".to_string()))
            .is_empty(),
        "unknown law yields an empty slice, not an invented card"
    );
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"card-assemble\",\"verdict\":\"pass\",\
         \"detail\":\"genesis card assembles; by-key and by-law indexes answer\"}}"
    );
}

#[test]
fn genesis_refuses_nonzero_revision_and_bad_models() {
    let mut nonzero = genesis_id();
    nonzero.revision = 3;
    assert!(matches!(
        MaterialCard::assemble(nonzero, ClaimSet::new(), Vec::new()),
        Err(MatDbError::RevisionNotZero { offered: 3 })
    ));

    let mut empty_block = j2_card();
    empty_block.parameters.clear();
    assert!(matches!(
        MaterialCard::assemble(genesis_id(), ClaimSet::new(), vec![empty_block]),
        Err(MatDbError::EmptyParameterBlock { .. })
    ));

    let mut nan_parameter = j2_card();
    nan_parameter.parameters.insert(
        "yield_stress".to_string(),
        LawParameter {
            value: f64::NAN,
            dims: STRESS_DIMS,
        },
    );
    assert!(matches!(
        MaterialCard::assemble(genesis_id(), ClaimSet::new(), vec![nan_parameter]),
        Err(MatDbError::NonFiniteParameter { .. })
    ));

    let mut unlicensed = j2_card();
    unlicensed.provenance.license = String::new();
    assert!(matches!(
        MaterialCard::assemble(genesis_id(), ClaimSet::new(), vec![unlicensed]),
        Err(MatDbError::MissingLicense { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"card-gates\",\"verdict\":\"pass\",\
         \"detail\":\"nonzero genesis revision and model-card pathologies refuse typed\"}}"
    );
}

#[test]
fn supersession_links_predecessor_hash_and_advances_revision() {
    let genesis = MaterialCard::assemble(genesis_id(), density_claims(), vec![j2_card()])
        .expect("genesis assembles");
    let genesis_hash = genesis.content_hash();

    let successor = MaterialCard::supersede(&genesis, density_claims(), vec![j2_card()])
        .expect("successor builds");
    assert_eq!(successor.id().revision, 1);
    assert_eq!(successor.id().chemistry, genesis.id().chemistry);
    assert_eq!(successor.supersedes(), Some(genesis_hash));
    assert_ne!(
        successor.content_hash(),
        genesis_hash,
        "revision + lineage are identity-bearing"
    );

    // The predecessor is untouched and both remain valid.
    assert_eq!(genesis.id().revision, 0);
    assert_eq!(genesis.content_hash(), genesis_hash);

    let third = MaterialCard::supersede(&successor, density_claims(), Vec::new())
        .expect("second supersession");
    assert_eq!(third.id().revision, 2);
    assert_eq!(third.supersedes(), Some(successor.content_hash()));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"supersession\",\"verdict\":\"pass\",\
         \"detail\":\"revision 0->1->2 with predecessor hashes bound; predecessors immutable\"}}"
    );
}

#[test]
fn model_card_content_identity_is_field_sensitive() {
    let base = j2_card();
    assert_eq!(base.content_hash(), j2_card().content_hash());

    let mut moved_value = j2_card();
    moved_value
        .parameters
        .get_mut("yield_stress")
        .expect("parameter exists")
        .value = 276.0e6 + 1.0;
    assert_ne!(base.content_hash(), moved_value.content_hash());

    let mut moved_version = j2_card();
    moved_version.law_version = 2;
    assert_ne!(base.content_hash(), moved_version.content_hash());

    let mut moved_state = j2_card();
    moved_state.initial_state = InitialStatePolicy::RequiresDeclaredState;
    assert_ne!(base.content_hash(), moved_state.content_hash());

    let mut moved_validity = j2_card();
    moved_validity.validity = ValidityDomain::unconstrained().with("T", 200.0, 500.0);
    assert_ne!(base.content_hash(), moved_validity.content_hash());

    let mut moved_sources = j2_card();
    moved_sources.sources = vec![hash_bytes(b"tensile calibration v2")];
    assert_ne!(base.content_hash(), moved_sources.content_hash());
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"model-identity\",\"verdict\":\"pass\",\
         \"detail\":\"model-card hash stable on equal content, moves on every semantic field\"}}"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one mutation per identity field, in field order
fn canonical_parameter_block_hash_is_ordered_and_narrowly_scoped() {
    let base = j2_card();
    assert_eq!(
        parameter_hash(&base).to_hex(),
        "bea281f14bd420545be1cdacee1b0fdc4b63b29197d97639a4a734b3b05432ef"
    );
    assert_eq!(parameter_hash(&base), parameter_hash(&j2_card()));
    let mut reordered = base.clone();
    reordered.parameters = base
        .parameters
        .iter()
        .rev()
        .map(|(name, parameter)| (name.clone(), parameter.clone()))
        .collect();
    assert_eq!(
        parameter_hash(&base),
        parameter_hash(&reordered),
        "BTreeMap insertion history cannot move the canonical order"
    );

    let mut moved_value = j2_card();
    moved_value
        .parameters
        .get_mut("yield_stress")
        .expect("parameter exists")
        .value = 276.0e6 + 1.0;
    assert_ne!(parameter_hash(&base), parameter_hash(&moved_value));

    let mut positive_zero = j2_card();
    positive_zero
        .parameters
        .get_mut("yield_stress")
        .expect("parameter exists")
        .value = 0.0;
    let mut negative_zero = j2_card();
    negative_zero
        .parameters
        .get_mut("yield_stress")
        .expect("parameter exists")
        .value = -0.0;
    assert_ne!(
        parameter_hash(&positive_zero),
        parameter_hash(&negative_zero),
        "parameter identity preserves exact floating-point bits"
    );

    let mut moved_dims = j2_card();
    moved_dims
        .parameters
        .get_mut("yield_stress")
        .expect("parameter exists")
        .dims = Dims([0, 0, 0, 0, 0, 0]);
    assert_ne!(parameter_hash(&base), parameter_hash(&moved_dims));

    let mut renamed = j2_card();
    let yield_stress = renamed
        .parameters
        .remove("yield_stress")
        .expect("parameter exists");
    renamed
        .parameters
        .insert("yield_strength".to_string(), yield_stress);
    assert_ne!(parameter_hash(&base), parameter_hash(&renamed));

    let mut extra_parameter = j2_card();
    extra_parameter.parameters.insert(
        "viscosity".to_string(),
        LawParameter {
            value: 2.0,
            dims: Dims([0, 0, 0, 0, 0, 0]),
        },
    );
    assert_ne!(parameter_hash(&base), parameter_hash(&extra_parameter));

    let base_hash = parameter_hash(&base);
    let mut moved_law = j2_card();
    moved_law.law = LawId("different-law".to_string());
    assert_eq!(base_hash, parameter_hash(&moved_law));
    let mut moved_law_version = j2_card();
    moved_law_version.law_version = 9;
    assert_eq!(base_hash, parameter_hash(&moved_law_version));
    let mut moved_state_schema = j2_card();
    moved_state_schema.state_schema_version = 7;
    assert_eq!(base_hash, parameter_hash(&moved_state_schema));
    let mut moved_initial_state = j2_card();
    moved_initial_state.initial_state = InitialStatePolicy::RequiresDeclaredState;
    assert_eq!(base_hash, parameter_hash(&moved_initial_state));
    let mut moved_validity = j2_card();
    moved_validity.validity = ValidityDomain::unconstrained().with("T", 1.0, 2.0);
    assert_eq!(base_hash, parameter_hash(&moved_validity));
    let mut moved_sources = j2_card();
    moved_sources.sources = vec![hash_bytes(b"different source")];
    assert_eq!(base_hash, parameter_hash(&moved_sources));
    let mut moved_provenance = j2_card();
    moved_provenance.provenance = Provenance {
        source: "independent parameter audit".to_string(),
        license: "test-only".to_string(),
        artifact: Some(hash_bytes(b"parameter audit")),
    };
    assert_eq!(
        base_hash,
        parameter_hash(&moved_provenance),
        "non-parameter semantics and provenance bind outside this narrow identity"
    );

    let mut invalid = j2_card();
    invalid.parameters.clear();
    assert!(matches!(
        invalid.canonical_parameters_hash(),
        Err(MatDbError::EmptyParameterBlock { .. })
    ));
}

#[test]
fn material_card_hash_binds_claims_models_and_lineage() {
    let base = MaterialCard::assemble(genesis_id(), density_claims(), vec![j2_card()])
        .expect("genesis assembles");

    let fewer_claims = MaterialCard::assemble(genesis_id(), ClaimSet::new(), vec![j2_card()])
        .expect("no-claims card");
    assert_ne!(base.content_hash(), fewer_claims.content_hash());

    let fewer_models =
        MaterialCard::assemble(genesis_id(), density_claims(), Vec::new()).expect("no-models card");
    assert_ne!(base.content_hash(), fewer_models.content_hash());

    let mut other_id = genesis_id();
    other_id.process = "T4".to_string();
    let other_process =
        MaterialCard::assemble(other_id, density_claims(), vec![j2_card()]).expect("other process");
    assert_ne!(base.content_hash(), other_process.content_hash());
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"material-identity\",\"verdict\":\"pass\",\
         \"detail\":\"card hash binds claims, models, and the named-state id\"}}"
    );
}
