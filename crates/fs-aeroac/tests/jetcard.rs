//! Jet-card battery (bead frankensim-music-v8-root-3ez8g.10.2):
//! schema round trip, content hashing, the authority law, minting
//! refusals, and the tonal interim card's recorded facts.

use fs_aeroac::jetcard::{
    CardAuthority, CardResidual, JET_CARD_SCHEMA, JetCard, JetCardClaim, MeanJetProfile,
    jet_labium_fingerprint, mint_broadband_card, mint_refusal_boundary_card,
    mint_tonal_interim_card, momentum_thickness_smoothed_tophat, staging_rig_config,
};
use fs_aeroac::noisetable::N_BANDS;
use fs_aeroac::slot_jet_3d::SlotJet3dRung;
use fs_aeroac::{AeroacError, SCOPE_STATEMENT};

/// A synthetic admitted rung (in-regime, amplitude-qualified).
fn rung(reynolds: f64, tonal: bool) -> SlotJet3dRung {
    SlotJet3dRung {
        reynolds,
        second_order_rate: 1.9,
        higher_order_rate: 2.0,
        flatness: if tonal { 1.0e-17 } else { 3.0e-3 },
        tonal,
        strouhal: 0.2,
        peak_bin: 40,
        prominence: 1.0e3,
        force_rms: 1.0e-3,
        amplitude_qualified: true,
        strouhal_bin_width: 1.0e-3,
        mach_max_lattice: 0.07,
        flux_imbalance: 0.01,
    }
}

fn profile() -> MeanJetProfile {
    MeanJetProfile {
        u_centerline: 0.04,
        slot_half: 2.5,
        momentum_thickness: 1.1,
    }
}

#[test]
fn jc_001_tonal_interim_card_mints_the_recorded_stage_one_point() {
    let card = mint_tonal_interim_card().expect("the pinned staging record mints");
    card.validate().expect("minted card validates");
    // The narrow claim: edge-tone class, stage I, the recorded lock.
    let JetCardClaim::EdgeToneTonal { feedback } = &card.claim else {
        panic!("tonal interim card must carry the edge-tone claim");
    };
    assert_eq!(feedback.stage, 1);
    assert!((feedback.locked_strouhal - 0.036_62).abs() < 1e-12);
    assert!(feedback.hysteresis_recorded);
    assert_eq!(feedback.multi_stable_strouhal.len(), 2);
    // Authority law: the residual against Brown is present, so the
    // card is X-Struct, and the residual reproduces the recorded
    // +3.0% inside the ±6% bin quantization.
    assert_eq!(card.authority, CardAuthority::XStruct);
    let residual = card.residual.as_ref().expect("X-Struct carries a residual");
    let dev = residual.relative_deviation();
    assert!((dev - 0.030_4).abs() < 5e-4, "deviation {dev:.4}");
    assert!(dev.abs() < residual.bin_halfwidth_rel);
    // NOT flute broadband: no band noise content on a tonal claim.
    assert!(card.band_db.is_none());
    // Provenance binds the exact staging rig.
    assert_eq!(
        card.provenance.rig_fingerprint,
        jet_labium_fingerprint(&staging_rig_config())
    );
    // Scope law embedded verbatim (marketing-mutation guard).
    assert_eq!(card.scope, SCOPE_STATEMENT);
    // Validity: the demonstrated ramp range only.
    card.admit_query(144.0).expect("inside");
    card.admit_query(264.0).expect("upper edge");
    assert!(card.admit_query(500.0).is_err(), "no extrapolation");
    assert!(card.admit_query(f64::NAN).is_err());
    // The receipt: card fields against the lab record, as the
    // canonical serialized document plus its content hash.
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"jc-001-tonal-interim-card\",\
         \"content_hash\":\"{:016x}\",\"verdict\":\"pass\"}}",
        card.content_hash()
    );
    println!("{}", card.to_json());
}

#[test]
fn jc_002_round_trip_is_exact_and_content_hash_is_stable() {
    let card = mint_tonal_interim_card().expect("mints");
    let json = card.to_json();
    let back = JetCard::from_json(&json).expect("round trip parses");
    assert_eq!(card, back, "round trip must be exact");
    assert_eq!(card.content_hash(), back.content_hash());
    // Deterministic across independent mints.
    let again = mint_tonal_interim_card().expect("mints again");
    assert_eq!(card.content_hash(), again.content_hash());
    // A changed field changes the hash (the hash is over content,
    // not identity).
    let mut mutated = card.clone();
    mutated.validity.reynolds_hi = 265.0;
    assert_ne!(card.content_hash(), mutated.content_hash());
    // The optional saturation amplitude survives the round trip in
    // both states.
    let mut with_sat = card.clone();
    if let JetCardClaim::EdgeToneTonal { feedback } = &mut with_sat.claim {
        feedback.saturated_force_rms = Some(2.5e-4);
    }
    let back = JetCard::from_json(&with_sat.to_json()).expect("saturated variant parses");
    assert_eq!(with_sat, back);
}

#[test]
fn jc_003_round_trip_covers_band_content_and_broadband_claim() {
    let bands = {
        let mut b = [0.0f64; N_BANDS];
        for (i, v) in b.iter_mut().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            {
                *v = -(i as f64) * 1.5;
            }
        }
        b
    };
    let card = mint_broadband_card(
        &[rung(240.0, false), rung(120.0, true)],
        profile(),
        bands,
        Some(CardResidual {
            quantity: "band shape".to_owned(),
            measured: 1.0,
            reference: 1.05,
            reference_source: "synthetic lab row \"quoted\" and back\\slashed".to_owned(),
            bin_halfwidth_rel: 0.1,
        }),
        0xDEAD_BEEF,
        vec!["synthetic receipt".to_owned()],
    )
    .expect("broadband mints from an admitted broadband rung");
    assert_eq!(card.authority, CardAuthority::XStruct);
    let json = card.to_json();
    let back = JetCard::from_json(&json).expect("parses");
    assert_eq!(card, back);
}

#[test]
fn jc_004_authority_law_caps_no_residual_at_x_est() {
    let bands = [-3.0f64; N_BANDS];
    let card = mint_broadband_card(&[rung(240.0, false)], profile(), bands, None, 1, vec![])
        .expect("mints without a residual");
    assert_eq!(
        card.authority,
        CardAuthority::XEst,
        "no residual against lab data must cap authority at X-Est"
    );
    // And a hand-built X-Struct card without a residual refuses.
    let mut forged = card;
    forged.authority = CardAuthority::XStruct;
    let err = forged.validate().expect_err("forged authority refuses");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
}

#[test]
fn jc_005_broadband_mint_refuses_without_a_broadband_rung() {
    let err = mint_broadband_card(
        &[rung(120.0, true), rung(240.0, true)],
        profile(),
        [-3.0f64; N_BANDS],
        None,
        1,
        vec![],
    )
    .expect_err("all-tonal rungs cannot back a broadband card");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
    // An out-of-regime broadband rung is equally refused.
    let mut hot = rung(240.0, false);
    hot.mach_max_lattice = 0.4;
    let err = mint_broadband_card(&[hot], profile(), [-3.0f64; N_BANDS], None, 1, vec![])
        .expect_err("out-of-regime rung cannot back a card");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
}

#[test]
fn jc_006_refusal_boundary_card_is_the_honest_all_tonal_artifact() {
    let card = mint_refusal_boundary_card(
        &[rung(80.0, true), rung(160.0, true), rung(240.0, true)],
        profile(),
        7,
        vec!["sweep receipts".to_owned()],
    )
    .expect("all-tonal admitted sweep mints the boundary");
    let JetCardClaim::BroadbandRefusalBoundary {
        max_reynolds_probed,
        rungs_probed,
    } = card.claim
    else {
        panic!("wrong claim kind");
    };
    assert!((max_reynolds_probed - 240.0).abs() < 1e-12);
    assert_eq!(rungs_probed, 3);
    assert_eq!(card.authority, CardAuthority::XEst);
    // Round trip this kind too.
    let back = JetCard::from_json(&card.to_json()).expect("parses");
    assert_eq!(card, back);
    // A broadband rung in the set means this mint is a lie: refuse.
    let err = mint_refusal_boundary_card(
        &[rung(80.0, true), rung(240.0, false)],
        profile(),
        7,
        vec![],
    )
    .expect_err("broadband rung present: the boundary card refuses");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
    // One rung is not a boundary.
    let err = mint_refusal_boundary_card(&[rung(80.0, true)], profile(), 7, vec![])
        .expect_err("a single rung is not a boundary");
    assert!(matches!(err, AeroacError::InvalidParameter { .. }));
}

#[test]
fn jc_007_parser_refuses_foreign_and_tampered_bytes() {
    let card = mint_tonal_interim_card().expect("mints");
    let json = card.to_json();
    // Foreign schema refuses by name.
    let foreign = json.replacen("jet-card/v1", "jet-card/v2", 1);
    assert!(JetCard::from_json(&foreign).is_err());
    // A mutated scope statement refuses (marketing-mutation guard
    // holds through deserialization).
    let laundered = json.replacen("NOT absolute SPL predictions", "absolute SPL", 1);
    assert!(JetCard::from_json(&laundered).is_err());
    // Trailing bytes refuse.
    let mut trailing = json.clone();
    trailing.push('x');
    assert!(JetCard::from_json(&trailing).is_err());
    // Truncation refuses.
    assert!(JetCard::from_json(&json[..json.len() - 2]).is_err());
    // Arbitrary JSON refuses.
    assert!(JetCard::from_json("{\"schema\":\"nope\"}").is_err());
    assert!(JetCard::from_json("").is_err());
}

#[test]
fn jc_008_momentum_thickness_matches_the_tanh_edge_law() {
    // For b >> w the smoothed top-hat's momentum thickness
    // approaches w (two tanh edges contributing w/2 each).
    let theta = momentum_thickness_smoothed_tophat(12.0, 1.0, 256).expect("computes");
    assert!(
        (theta - 1.0).abs() < 0.05,
        "theta {theta:.4} should approach the smoothing width"
    );
    // The staging geometry's value is positive and below the slot
    // width (a sanity envelope, not a claim).
    let cfg = staging_rig_config();
    let staging = momentum_thickness_smoothed_tophat(cfg.slot_half, cfg.slot_smoothing, cfg.ny)
        .expect("staging geometry computes");
    assert!(staging > 0.0 && staging < 2.0 * cfg.slot_half);
    // Refusals: non-positive parameters and a too-small domain.
    assert!(momentum_thickness_smoothed_tophat(0.0, 1.0, 64).is_err());
    assert!(momentum_thickness_smoothed_tophat(3.0, 0.0, 64).is_err());
    assert!(momentum_thickness_smoothed_tophat(30.0, 1.0, 64).is_err());
}

#[test]
fn jc_009_schema_constant_is_pinned() {
    assert_eq!(JET_CARD_SCHEMA, "fs-aeroac.jet-card/v1");
    let card = mint_tonal_interim_card().expect("mints");
    assert!(
        card.to_json()
            .starts_with("{\"schema\":\"fs-aeroac.jet-card/v1\"")
    );
    assert!(
        card.to_json()
            .contains(&format!("\"claim_kind\":\"{}\"", card.claim.kind()))
    );
}
