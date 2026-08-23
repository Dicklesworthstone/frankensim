//! Selector gates (music bead `frankensim-music-v8-root-3ez8g.14`).
//!
//! The DONE-WHEN demos, driven through the REAL tracked registry
//! (`instrument-claims.json`) and the REAL measured hop-policy artifact
//! (`data/claims/wind-hop-policy.tsv`) — never fixture stand-ins for those
//! two inputs:
//!
//! 1. a two-image wind voice (char-line <-> vfit-hold) switches under a
//!    gesture schedule with receipts at every boundary;
//! 2. a forced budget squeeze drops to the cheapest admitted image, and a
//!    POISONED registry (ungated fallback / live-default-without-budget)
//!    is refused, proving the selector can never select off-menu or
//!    ungated as a live default;
//! 3. an FD-claim request routes to the oracle without rendering audio;
//! 4. selection receipts round-trip byte-stably and carry a content digest
//!    callers cite from render provenance.

use fs_couple::selector::{
    EntryClass, GesturePhase, HopPolicy, HopVerdict, MenuEntry, RegistrySnapshot, SelectionRequest,
    Selector, SelectorError, Serving, SessionPosture, VoiceMenu,
};

fn repo_root() -> std::path::PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    std::path::Path::new(manifest)
        .join("../../")
        .canonicalize()
        .expect("repo root")
}

fn tracked_registry() -> String {
    std::fs::read_to_string(repo_root().join("instrument-claims.json"))
        .expect("tracked instrument-claims.json")
}

fn tracked_hop_policy() -> String {
    std::fs::read_to_string(repo_root().join("data/claims/wind-hop-policy.tsv"))
        .expect("committed wind-hop-policy.tsv")
}

/// The real wind-reed menu: phrase rides char-line, settled notes hop to
/// vfit-hold, spectral questions go to the TMM oracle. Cost tiers follow
/// the doctrine: hold images are cheaper than full characteristic lines.
fn wind_menu() -> VoiceMenu {
    VoiceMenu {
        filling: "wind-reed",
        entries: vec![
            MenuEntry {
                image: "tmm",
                qoi: "fingering-peak-cents",
                cost_tier: 0,
                class: EntryClass::SpectralOracle,
                lift_supported: false,
            },
            MenuEntry {
                image: "vfit-hold",
                qoi: "held-fingering-cents",
                cost_tier: 1,
                class: EntryClass::SettledHold,
                lift_supported: true,
            },
            MenuEntry {
                image: "char-line",
                qoi: "junction-match-tmm-r",
                cost_tier: 3,
                class: EntryClass::SpatialOrIsland,
                lift_supported: true,
            },
        ],
    }
}

#[test]
fn registry_and_policy_parse_from_tracked_artifacts() {
    let reg = RegistrySnapshot::parse(&tracked_registry()).expect("registry parses");
    assert!(reg.len() >= 40, "tracked registry rows present");
    // The exact claims the menu cites exist and are green in the tracked file.
    for (image, qoi) in [
        ("tmm", "fingering-peak-cents"),
        ("vfit-hold", "held-fingering-cents"),
        ("char-line", "junction-match-tmm-r"),
    ] {
        let row = reg
            .row("wind-reed", image, qoi)
            .unwrap_or_else(|| panic!("row wind-reed/{image}/{qoi} present"));
        assert_eq!(
            row.gate,
            fs_couple::selector::Gate::Green,
            "{image} must be gated green for this demo"
        );
    }
    let policy = HopPolicy::parse(&tracked_hop_policy()).expect("policy parses");
    assert!((policy.drift_threshold - 0.05).abs() < f64::EPSILON);
    assert_eq!(policy.consecutive_blocks, 4);
    assert!((policy.window_ms - 25.0).abs() < f64::EPSILON);
    assert_eq!(policy.first_settled_block, 6);
}

#[test]
fn two_image_wind_voice_switches_under_measured_schedule() {
    let reg = RegistrySnapshot::parse(&tracked_registry()).expect("registry");
    let policy = HopPolicy::parse(&tracked_hop_policy()).expect("policy");
    let menu = wind_menu();
    let mut receipts = Vec::new();
    let mut current: Option<&'static str> = None;

    // Measured-schedule drive: attack blocks (drift above threshold), then a
    // settling tail. The policy decides Held, not this test's intuition.
    let drifts = [0.90, 0.72, 0.55, 0.31, 0.12, 0.04, 0.03, 0.02, 0.02, 0.01];
    for (block, drift) in drifts.iter().enumerate() {
        let block_idx = block as u32;
        let streak = drifts[..=block]
            .iter()
            .rev()
            .take_while(|d| **d < policy.drift_threshold)
            .count() as u32;
        let phase = policy.classify(*drift, streak, block_idx);
        let boundary = current.is_none();
        let receipt = Selector::select(
            &reg,
            &menu,
            &SelectionRequest {
                serving: Serving::Play,
                phase,
                budget_headroom_tier: 9,
                current_image: current,
                at_phrase_boundary: boundary || block_idx == 0,
                posture: SessionPosture::LiveDefault,
            },
        )
        .unwrap_or_else(|e| panic!("block {block}: {e}"));
        assert!(receipt.render_audio);
        if let Some((image, _)) = &receipt.chosen {
            current = Some(
                menu.entries
                    .iter()
                    .find(|e| e.image == image.as_str())
                    .expect("chosen image on menu")
                    .image,
            );
        }
        receipts.push((receipt.to_jsonl(), boundary));
    }

    // Phrase blocks ride char-line; held blocks hop to the cheaper hold.
    let first_char = receipts
        .iter()
        .position(|(l, _)| l.contains("\"chosen_image\":\"char-line\""))
        .expect("char chosen");
    let first_hold = receipts
        .iter()
        .position(|(l, _)| l.contains("\"chosen_image\":\"vfit-hold\""))
        .expect("hold chosen");
    assert!(first_char < first_hold, "phrase rides char before hold");

    // Every decision logged; the held-phase receipt names its reason.
    let held_line = &receipts[first_hold].0;
    assert!(held_line.contains("held-prefers-gated-hold"));

    // Determinism: same schedule twice => byte-identical receipt stream.
    let mut again = Vec::new();
    let mut current2: Option<&'static str> = None;
    for (block, drift) in drifts.iter().enumerate() {
        let block_idx = block as u32;
        let streak = drifts[..=block]
            .iter()
            .rev()
            .take_while(|d| **d < policy.drift_threshold)
            .count() as u32;
        let phase = policy.classify(*drift, streak, block_idx);
        let boundary = current2.is_none();
        let receipt = Selector::select(
            &reg,
            &menu,
            &SelectionRequest {
                serving: Serving::Play,
                phase,
                budget_headroom_tier: 9,
                current_image: current2,
                at_phrase_boundary: boundary,
                posture: SessionPosture::LiveDefault,
            },
        )
        .expect("second pass");
        if let Some((image, _)) = &receipt.chosen {
            current2 = Some(
                menu.entries
                    .iter()
                    .find(|e| e.image == image.as_str())
                    .expect("on menu")
                    .image,
            );
        }
        again.push(receipt.to_jsonl());
    }
    let first_pass: Vec<String> = receipts.into_iter().map(|(l, _)| l).collect();
    assert_eq!(first_pass, again, "byte-identical decisions across runs");
}

#[test]
fn d25_live_default_requires_green_and_budget_row_together() {
    // Focused unit fixture for D25's conjunction on a green row CLAIMING the
    // live-default slot: no budget row -> refused; budget row added ->
    // admitted. Non-default posture admits regardless (it never claims the
    // slot).
    let without_budget = r#"{
  "schema": "frankensim-instrument-claims-v1",
  "rows": [
    {
      "filling": "test-filling",
      "image": "probe-image",
      "qoi": "probe-qoi",
      "owner_crates": ["fs-couple"],
      "exactness": ["X-Struct"],
      "gate": "green",
      "live_default": "yes",
      "determinism": "one-host",
      "evidence": [],
      "budget_row": null,
      "corpus_refs": [],
      "notes": ""
    }
  ]
}"#;
    let reg = RegistrySnapshot::parse(without_budget).expect("fixture parses");
    let menu = VoiceMenu {
        filling: "test-filling",
        entries: vec![MenuEntry {
            image: "probe-image",
            qoi: "probe-qoi",
            cost_tier: 1,
            class: EntryClass::SettledHold,
            lift_supported: false,
        }],
    };
    let err = Selector::select(
        &reg,
        &menu,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Held,
            budget_headroom_tier: 9,
            current_image: None,
            at_phrase_boundary: true,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect_err("live default without budget row refuses");
    match &err {
        SelectorError::NoAdmittedImage(d) => {
            assert!(d.contains("live-default-without-budget-row"), "{d}");
        }
        other => panic!("{other}"),
    }

    let with_budget = without_budget.replace(
        "\"budget_row\": null",
        "\"budget_row\": \"data/budget-rows/probe.budgetrow\"",
    );
    let reg_ok = RegistrySnapshot::parse(&with_budget).expect("fixture parses");
    let receipt = Selector::select(
        &reg_ok,
        &menu,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Held,
            budget_headroom_tier: 9,
            current_image: None,
            at_phrase_boundary: true,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect("budget row completes the conjunction");
    assert_eq!(
        receipt.chosen.as_ref().map(|c| c.0.as_str()),
        Some("probe-image")
    );
}

#[test]
fn budget_squeeze_drops_to_cheapest_admitted_and_never_off_menu() {
    let reg = RegistrySnapshot::parse(&tracked_registry()).expect("registry");
    let menu = wind_menu();
    // Headroom only covers the oracle tier: play candidates over budget.
    let receipt = Selector::select(
        &reg,
        &menu,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Phrase,
            budget_headroom_tier: 2,
            current_image: None,
            at_phrase_boundary: false,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect("squeeze still admits something");
    assert!(receipt.fallback_used, "fallback recorded");
    let (img, _) = receipt.chosen.expect("chosen");
    assert_eq!(img, "vfit-hold", "cheapest admitted non-oracle entry");
    assert_eq!(
        receipt.hop,
        HopVerdict::MidNoteAllowed,
        "lift supported mid-note"
    );

    // Poisoned variant: make EVERY candidate refuse (ungated fallback under
    // live-default posture) and watch a typed refusal come back instead of
    // an off-menu pick.
    let poisoned = tracked_registry().replace("\"gate\": \"green\"", "\"gate\": \"ungated\"");
    let reg_bad = RegistrySnapshot::parse(&poisoned).expect("structurally fine");
    let err = Selector::select(
        &reg_bad,
        &menu,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Phrase,
            budget_headroom_tier: 9,
            current_image: None,
            at_phrase_boundary: true,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect_err("all-ungated refuses under live default");
    assert!(matches!(err, SelectorError::NoAdmittedImage(_)));
    let detail = match &err {
        SelectorError::NoAdmittedImage(d) => d.clone(),
        other => panic!("{other}"),
    };
    // The oracle is refused by serving discipline before gates are even
    // consulted; the audio images are refused by the live-default gate on
    // ungated rows.
    assert!(
        detail.contains("tmm=oracle-not-for-audio"),
        "oracle never renders audio: {detail}"
    );
    for entry in ["vfit-hold", "char-line"] {
        assert!(
            detail.contains(&format!("{entry}=ungated-not-live-default")),
            "entry {entry} refused as ungated: {detail}"
        );
    }

    // Unknown-filling refusals name the filling.
    let bogus = VoiceMenu {
        filling: "kazoo",
        entries: menu.entries.clone(),
    };
    let err = Selector::select(
        &reg,
        &bogus,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Phrase,
            budget_headroom_tier: 9,
            current_image: None,
            at_phrase_boundary: true,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect_err("unknown filling refuses");
    assert_eq!(
        err,
        SelectorError::UnknownFilling("kazoo".to_string()),
        "typed unknown-filling refusal"
    );
}

#[test]
fn fd_claim_routes_to_oracle_without_audio() {
    let reg = RegistrySnapshot::parse(&tracked_registry()).expect("registry");
    let menu = wind_menu();
    let receipt = Selector::select(
        &reg,
        &menu,
        &SelectionRequest {
            serving: Serving::PeakLocations,
            phase: GesturePhase::Phrase,
            budget_headroom_tier: 9,
            current_image: Some("char-line"),
            at_phrase_boundary: false,
            posture: SessionPosture::LiveDefault,
        },
    )
    .expect("oracle route");
    let (img, _) = receipt.chosen.expect("chosen");
    assert_eq!(img, "tmm", "spectral question routes to the FD oracle");
    assert!(
        !receipt.render_audio,
        "the oracle answers WITHOUT rendering audio"
    );
    assert!(
        receipt.chosen_because.contains("no-audio"),
        "reason names the no-audio routing"
    );
}

#[test]
fn receipts_round_trip_and_digest_stably() {
    let reg = RegistrySnapshot::parse(&tracked_registry()).expect("registry");
    let menu = wind_menu();
    let receipt = Selector::select(
        &reg,
        &menu,
        &SelectionRequest {
            serving: Serving::Play,
            phase: GesturePhase::Held,
            budget_headroom_tier: 9,
            current_image: Some("char-line"),
            at_phrase_boundary: false,
            posture: SessionPosture::NonDefaultDeclared,
        },
    )
    .expect("receipt");
    let line = receipt.to_jsonl();
    let decoded = fs_couple::selector::SelectionReceipt::from_jsonl(&line).expect("round trip");
    assert_eq!(decoded, receipt, "lossless round trip");
    assert_eq!(decoded.digest().to_string(), receipt.digest().to_string());
    assert_eq!(decoded.to_jsonl(), line, "canonical bytes stable");

    // A mutated receipt changes its digest (tamper-evident provenance).
    let mut tampered = decoded.clone();
    tampered.chosen_because = "forged".to_string();
    assert_ne!(tampered.digest().to_string(), receipt.digest().to_string());

    // Schema mismatch refuses.
    let bad = line.replace(fs_couple::selector::SELECTION_RECEIPT_SCHEMA, "other-v9");
    assert!(fs_couple::selector::SelectionReceipt::from_jsonl(&bad).is_err());
}
