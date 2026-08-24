//! Freeze-integrity gate battery for the RANS model card
//! (bead frankensim-extreal-program-f85xj.5.8.1).
//!
//! The controlling idea: a card that omits a governing term, flips a
//! coefficient sign, silently enables transition, or erases the no-claim
//! surface MUST be refused at freeze. Every mutant below is constructed
//! from the canonical draft by exactly one deliberate defect.

use fs_scenario::rans_card::{
    AdmissionError, BoussinesqOption, MAX_MANIFEST_BYTES, PorousFinSink, REQUIRED_TERMS,
    RansCardDraft, RansModelCard, WallTreatment,
};

fn canonical() -> RansCardDraft {
    RansCardDraft::launder_sharma_channel("electronics-cooling/e10-rans")
}

fn mutate<F: FnOnce(&mut RansCardDraft)>(f: F) -> Result<RansModelCard, AdmissionError> {
    let mut d = canonical();
    f(&mut d);
    d.freeze()
}

fn assert_ambiguous_list<F: FnOnce(&mut RansCardDraft)>(field: &'static str, f: F) {
    let error = mutate(f)
        .err()
        .expect("ambiguous flattened-list input must be refused");
    assert_eq!(error, AdmissionError::AmbiguousManifestList { field });
}

#[test]
fn canonical_draft_freezes_with_stable_binding() {
    let card = canonical().freeze().expect("canonical card admits");
    assert_eq!(
        card.schema(),
        fs_scenario::rans_card::RANS_MODEL_CARD_SCHEMA
    );
    let again = canonical().freeze().expect("re-freeze");
    assert_eq!(
        card.manifest_hash().to_hex(),
        again.manifest_hash().to_hex(),
        "identical drafts must bind identical manifest bytes"
    );
    // The frozen card is opaque: no field of the struct is public, so the
    // only mutation path is a new draft + new freeze.
}

#[test]
fn embedded_nul_cannot_alias_manifest_pair_boundaries() {
    let mut smuggled = canonical();
    smuggled.boundary_conditions.clear();
    smuggled
        .boundary_conditions
        .insert("a".to_string(), "v\0bc/b\0w".to_string());

    let mut split = canonical();
    split.boundary_conditions.clear();
    split
        .boundary_conditions
        .insert("a".to_string(), "v".to_string());
    split
        .boundary_conditions
        .insert("b".to_string(), "w".to_string());
    let split = split.freeze().expect("ordinary separate rows admit");

    let smuggled = smuggled.freeze();
    if let Ok(card) = &smuggled {
        assert_ne!(
            card.statement_manifest(),
            split.statement_manifest(),
            "the cards are semantically distinct"
        );
        assert_ne!(
            card.manifest_hash(),
            split.manifest_hash(),
            "embedded delimiters must not alias separate manifest rows"
        );
    }
    assert!(
        matches!(smuggled, Err(AdmissionError::EmbeddedManifestDelimiter)),
        "the frozen v1 delimiter format must refuse embedded NULs"
    );
}

#[test]
fn embedded_list_delimiters_cannot_alias_distinct_items() {
    let mut smuggled = canonical();
    smuggled.governing_terms.push("extra-a,extra-b".to_string());

    let mut split = canonical();
    split.governing_terms.push("extra-a".to_string());
    split.governing_terms.push("extra-b".to_string());
    let split = split.freeze().expect("ordinary separate list items admit");

    let smuggled = smuggled.freeze();
    if let Ok(card) = &smuggled {
        assert_ne!(
            card.manifest_hash(),
            split.manifest_hash(),
            "one item containing the list separator must not alias two items"
        );
    }
    assert_eq!(
        smuggled.err(),
        Some(AdmissionError::AmbiguousManifestList {
            field: "governing_terms",
        }),
        "the v1 flattened-list format must refuse ambiguous list items"
    );
}

#[test]
fn every_caller_controlled_flattened_list_has_injective_admission() {
    assert_ambiguous_list("governing_terms", |draft| {
        draft.governing_terms.push(String::new());
    });
    assert_ambiguous_list("validation_case_families", |draft| {
        draft.validation_case_families = vec!["case-a,case-b".to_string()];
    });
    assert_ambiguous_list("falsifiers", |draft| {
        draft.falsifiers = vec!["claim-a;claim-b".to_string()];
    });
    assert_ambiguous_list("exclusions", |draft| {
        draft.exclusions.push("claim-a;claim-b".to_string());
    });
}

#[test]
fn manifest_size_cap_counts_hash_framing_bytes() {
    const CANONICAL_FAMILY: &str = "electronics-cooling/e10-rans";
    let canonical_card = RansCardDraft::launder_sharma_channel(CANONICAL_FAMILY)
        .freeze()
        .expect("canonical card admits");
    let base_len = canonical_card
        .statement_manifest()
        .iter()
        .map(|(key, value)| key.len() + 1 + value.len() + 1)
        .sum::<usize>()
        - CANONICAL_FAMILY.len();
    assert!(base_len < MAX_MANIFEST_BYTES);

    let oversized_family = "x".repeat(MAX_MANIFEST_BYTES - base_len + 1);
    let error = RansCardDraft::launder_sharma_channel(oversized_family)
        .freeze()
        .err()
        .expect("one byte beyond the framed manifest cap must be refused");
    assert_eq!(
        error,
        AdmissionError::ManifestTooLarge {
            actual: MAX_MANIFEST_BYTES + 1,
        }
    );
}

#[test]
fn mutant_omitting_a_governing_term_is_refused() {
    for term in REQUIRED_TERMS {
        let err = mutate(|d| {
            d.governing_terms.retain(|g| g != term);
        })
        .err()
        .expect("omitting a required term must refuse");
        assert_eq!(err, AdmissionError::MissingTerm { term }, "term {term}");
    }
}

#[test]
fn mutants_changing_coefficient_value_or_sign_are_refused() {
    // Out-of-bound high.
    let err = mutate(|d| d.coefficients.c_eps_2 = 2.5).err().unwrap();
    assert!(matches!(err, AdmissionError::CoefficientOutOfBounds { .. }));
    // Sign flip on the eddy-viscosity coefficient.
    let err = mutate(|d| d.coefficients.c_mu = -0.09).err().unwrap();
    assert!(matches!(
        err,
        AdmissionError::CoefficientOutOfBounds { .. } | AdmissionError::NonFinite { .. }
    ));
    // Non-finite injection.
    let err = mutate(|d| d.coefficients.sigma_k = f64::NAN).err().unwrap();
    assert_eq!(err, AdmissionError::NonFinite { field: "sigma_k" });
    // Turbulent Prandtl outside the constant-model band.
    let err = mutate(|d| d.turbulent_prandtl = 2.0).err().unwrap();
    assert!(matches!(
        err,
        AdmissionError::CoefficientOutOfBounds {
            name: "turbulent_prandtl"
        }
    ));
}

#[test]
fn mutants_enabling_transition_or_erasing_no_claim_are_refused() {
    // Silently enable transition: replace the refusal with its opposite.
    let err = mutate(|d| {
        d.exclusions = d
            .exclusions
            .iter()
            .map(|e| {
                if e.to_ascii_lowercase().contains("transition") {
                    "transitional flows are claimed valid".to_string()
                } else {
                    e.clone()
                }
            })
            .collect();
    })
    .err()
    .expect("claiming transition must refuse");
    assert!(matches!(err, AdmissionError::NoClaimViolation { .. }));

    // Erase the no-claim surface entirely.
    let err = mutate(|d| d.exclusions.clear()).err().unwrap();
    assert!(matches!(err, AdmissionError::NoClaimViolation { .. }));

    // A card without falsifiers cannot be frozen.
    let err = mutate(|d| d.falsifiers.clear()).err().unwrap();
    assert!(matches!(err, AdmissionError::NoClaimViolation { .. }));
}

#[test]
fn mutants_enabling_unavailable_capabilities_are_refused() {
    // Porous sink declared enabled but solver capability does not exist yet.
    let err = mutate(|d| {
        d.porous_fin = PorousFinSink {
            enabled: true,
            permeability_m2: Some(1.0e-9),
            forchheimer_c_f: Some(0.2),
        };
    })
    .err()
    .unwrap();
    assert_eq!(
        err,
        AdmissionError::CapabilityUnavailable {
            capability: "porous-sink"
        }
    );
    // Boussinesq same.
    let err = mutate(|d| {
        d.boussinesq = BoussinesqOption {
            enabled: true,
            beta_per_k: Some(3.4e-3),
            reference_temperature_k: 300.0,
        };
    })
    .err()
    .unwrap();
    assert_eq!(
        err,
        AdmissionError::CapabilityUnavailable {
            capability: "buoyancy-source"
        }
    );
    // Enabled but out-of-bound coefficients refuse BEFORE the capability
    // check matters (bounds are checked first).
    let err = mutate(|d| {
        d.porous_fin = PorousFinSink {
            enabled: true,
            permeability_m2: Some(1.0e-3),
            forchheimer_c_f: Some(0.2),
        };
    })
    .err()
    .unwrap();
    assert!(matches!(
        err,
        AdmissionError::CoefficientOutOfBounds { name: "porous_fin" }
    ));
}

#[test]
fn boundary_regimes_and_disabled_options_hold_on_the_frozen_card() {
    let card = canonical().freeze().expect("canonical admits");
    let (lo, hi) = card.reynolds_band();
    // Laminar-equivalent limit stays inside the applicability band.
    assert!(lo <= 500.0 && 500.0 <= hi);
    assert!(hi >= 10_000.0, "mid-band turbulence must be covered");
    // Wall treatment is resolve-to-sublayer and the y+ target travels.
    assert_eq!(card_wall(&card), WallTreatment::ResolveToViscousSublayer);
    // Disabled options stay honestly disabled.
    assert!(!card.boussinesq_enabled());
    assert!(!card.porous_fin_enabled());
    // Solver budgets are bounded and finite.
    assert!(card.max_iterations() > 0 && card.max_iterations() <= 100_000);
    assert!(card.residual_tolerance_rel() > 0.0 && card.residual_tolerance_rel() < 1.0);

    fn card_wall(card: &RansModelCard) -> WallTreatment {
        // Access through the manifest so the wall clause is part of the
        // adjudicator's binding bytes.
        let m = card.statement_manifest();
        let wall = m
            .iter()
            .find(|(k, _)| k == "wall_treatment")
            .map(|(_, v)| v.clone())
            .expect("wall treatment statement present");
        if wall.contains("ResolveToViscousSublayer") {
            WallTreatment::ResolveToViscousSublayer
        } else {
            panic!("unexpected wall treatment {wall}");
        }
    }
}

#[test]
fn manifest_is_sorted_complete_and_symmetry_consistent() {
    let card = canonical().freeze().expect("canonical admits");
    let m = card.statement_manifest();
    // Sorted and duplicate-free keys.
    let keys: Vec<&String> = m.iter().map(|(k, _)| k).collect();
    let mut sorted = keys.clone();
    sorted.sort();
    assert_eq!(keys, sorted, "manifest must be deterministically ordered");
    let mut uniq = sorted.clone();
    uniq.dedup();
    assert_eq!(sorted, uniq, "manifest keys must be unique");
    // Symmetry/metamorphic: reversing the Reynolds band order changes the
    // binding (the band is directional), while re-freezing an identical
    // draft does not (covered elsewhere).
    let flipped = mutate(|d| {
        d.reynolds_band = (d.reynolds_band.1 * 10.0, d.reynolds_band.0);
    });
    assert!(flipped.is_err(), "inverted band must refuse");
    // Completeness: every required term appears in the manifest terms row.
    let terms_row = m
        .iter()
        .find(|(k, _)| k == "governing_terms")
        .map(|(_, v)| v.clone())
        .unwrap_or_default();
    for term in REQUIRED_TERMS {
        assert!(terms_row.contains(term), "missing {term} in manifest");
    }
}

#[test]
fn oversized_manifest_blocks_freeze() {
    let err = mutate(|d| {
        d.system_family = "x".repeat(64 * 1024);
    })
    .err()
    .expect("oversized manifest must refuse");
    assert!(matches!(
        err,
        AdmissionError::ManifestTooLarge { actual: _ }
    ));
}
