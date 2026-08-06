//! G0/G3 tests for the complete cinematic configuration transport.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    cinematic_budget::{
        CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION, CinematicQualityProfile, CinematicQualityTier,
    },
    cinematic_config::{CinematicAssetBinding, CinematicAssetInterpretation},
    cinematic_config_codec::{
        CINEMATIC_CONFIG_DOCUMENT_SCHEMA, CinematicAssetAccessError, CinematicConfigDocument,
        CinematicConfigDocumentError, CinematicDocumentAssetClass,
        MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES, MAX_CINEMATIC_CONFIG_DOCUMENT_LINES,
    },
};

const MATERIAL_A: &[u8] = b"measured steel spectral reflectance v1";
const MATERIAL_B: &[u8] = b"measured tungsten spectral reflectance v1";
const LIGHT: &[u8] = b"area light spectral emission v1";
const ENVIRONMENT: &[u8] = b"studio environment spectral emission v1";

fn id(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.tests.cinematic-config-document.v1",
        label.as_bytes(),
    )
}

fn asset_id(bytes: &[u8], interpretation: CinematicAssetInterpretation) -> ContentHash {
    CinematicAssetBinding::from_bytes(bytes, interpretation, 1, "relocatable".to_owned())
        .expect("fixture asset")
        .content_identity()
}

fn component(label: &str) -> String {
    format!("1:{}", id(label).to_hex())
}

fn document_text(
    profile: CinematicQualityTier,
    material_a_locator: &str,
    material_b_locator: &str,
    mux: &str,
) -> String {
    let profile = CinematicQualityProfile::canonical(profile).expect("canonical profile");
    let profile_name = match profile.input().tier {
        CinematicQualityTier::StoryboardSmoke => "storyboard-smoke",
        CinematicQualityTier::Daily1080p => "daily-1080p",
        CinematicQualityTier::Qualification4kFrame => "qualification-4k-frame",
        CinematicQualityTier::Final4k => "final-4k",
    };
    let capabilities = if mux == "none" {
        "render,audio"
    } else {
        "render,audio,quarantined-mux"
    };
    let profile_ref = format!(
        "{}:{}",
        CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION,
        profile.identity().to_hex()
    );
    let material_a = asset_id(
        MATERIAL_A,
        CinematicAssetInterpretation::SpectralReflectance,
    );
    let material_b = asset_id(
        MATERIAL_B,
        CinematicAssetInterpretation::SpectralReflectance,
    );
    let light = asset_id(LIGHT, CinematicAssetInterpretation::SpectralEmission);
    let environment = asset_id(ENVIRONMENT, CinematicAssetInterpretation::SpectralEmission);
    format!(
        "schema={CINEMATIC_CONFIG_DOCUMENT_SCHEMA}\n\
         quality_profile={profile_name}\n\
         units=si-m-kg-s-rad\n\
         seed=7\n\
         capabilities={capabilities}\n\
         render_budget_profile={profile_ref}\n\
         audio_budget_profile={profile_ref}\n\
         trajectory={}\n\
         timeline={}\n\
         camera={}\n\
         scene_geometry={}\n\
         instance_mapping={}\n\
         renderer={}\n\
         image_pipeline={}\n\
         audio_excitation={}\n\
         sound_model={}\n\
         microphone={}\n\
         room={}\n\
         material_asset=spectral-reflectance:1:{material_a}:{material_a_locator}\n\
         material_asset=spectral-reflectance:1:{material_b}:{material_b_locator}\n\
         light_asset=spectral-emission:1:{light}:assets/light.spectrum\n\
         environment_asset=spectral-emission:1:{environment}:assets/environment.spectrum\n\
         artifact_namespace=euler/reference-v1\n\
         artifact_root=artifacts/euler/reference-v1\n\
         mux={mux}\n",
        component("trajectory"),
        component("timeline"),
        component("camera"),
        component("scene-geometry"),
        component("instance-mapping"),
        component("renderer"),
        component("image-pipeline"),
        component("audio-excitation"),
        component("sound-model"),
        component("microphone"),
        component("room"),
    )
}

fn resolve(
    _class: CinematicDocumentAssetClass,
    _index: usize,
    declaration: &fs_evidence::cinematic_config_codec::CinematicAssetDeclaration,
) -> Result<Vec<u8>, CinematicAssetAccessError> {
    let bytes = match declaration.locator_hint() {
        "assets/steel.spectrum" | "relocated/steel.spectrum" => MATERIAL_A,
        "assets/tungsten.spectrum" | "relocated/tungsten.spectrum" => MATERIAL_B,
        "assets/light.spectrum" => LIGHT,
        "assets/environment.spectrum" => ENVIRONMENT,
        _ => return Err(CinematicAssetAccessError::Unavailable),
    };
    Ok(bytes.to_vec())
}

#[test]
fn complete_document_round_trips_and_reenters_authoritative_config() {
    let source = document_text(
        CinematicQualityTier::StoryboardSmoke,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        "none",
    );
    let document = CinematicConfigDocument::from_bytes(source.as_bytes()).expect("decode");
    let canonical = document.to_canonical_text();
    let decoded = CinematicConfigDocument::from_str(&canonical).expect("canonical decode");
    assert_eq!(decoded, document);

    let config = decoded
        .admit_with_asset_resolver(resolve)
        .expect("asset-backed config admission");
    assert_eq!(config.input().seed, Some(7));
    assert_eq!(
        config.input().material_assets.len(),
        2,
        "both material declarations survive reconstruction"
    );
    assert_ne!(config.composition_identity(), ContentHash([0; 32]));
}

#[test]
fn comments_declaration_order_and_locator_relocation_do_not_change_composition() {
    let first = CinematicConfigDocument::from_str(&document_text(
        CinematicQualityTier::StoryboardSmoke,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        "none",
    ))
    .expect("first document");
    let reordered = document_text(
        CinematicQualityTier::StoryboardSmoke,
        "relocated/steel.spectrum",
        "relocated/tungsten.spectrum",
        "none",
    )
    .replace(
        "schema=frankensim.cinematic-config-document.v1\n",
        "# relocation must be operational only\nseed=7\nschema=frankensim.cinematic-config-document.v1\n",
    )
    .replacen("seed=7\n", "", 1);
    let second = CinematicConfigDocument::from_str(&reordered).expect("reordered document");

    let first = first
        .admit_with_asset_resolver(resolve)
        .expect("first admit");
    let second = second
        .admit_with_asset_resolver(resolve)
        .expect("second admit");
    assert_eq!(first.trajectory_identity(), second.trajectory_identity());
    assert_eq!(first.image_identity(), second.image_identity());
    assert_eq!(first.audio_identity(), second.audio_identity());
    assert_eq!(first.mux_identity(), second.mux_identity());
    assert_eq!(first.composition_identity(), second.composition_identity());
}

#[test]
fn closed_grammar_reports_stable_non_disclosing_field_paths() {
    let valid = document_text(
        CinematicQualityTier::Daily1080p,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        "none",
    );
    let unknown = format!("secret-token-value=do-not-echo\n{valid}");
    let error = CinematicConfigDocument::from_str(&unknown).expect_err("unknown field");
    assert_eq!(error.code(), "cinematic-document-unknown-field");
    assert_eq!(error.field_path(), "config.<unknown>");
    assert!(!error.to_string().contains("secret-token-value"));

    let duplicate = valid.replace("seed=7\n", "seed=7\nseed=8\n");
    assert!(matches!(
        CinematicConfigDocument::from_str(&duplicate),
        Err(CinematicConfigDocumentError::DuplicateField { field }) if field == "seed"
    ));

    let missing = valid.replace(&format!("room={}\n", component("room")), "");
    assert!(matches!(
        CinematicConfigDocument::from_str(&missing),
        Err(CinematicConfigDocumentError::MissingField { field }) if field == "room"
    ));
}

#[test]
fn named_profile_must_match_both_versioned_budget_references() {
    let source = document_text(
        CinematicQualityTier::Final4k,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        "none",
    );
    let storyboard = CinematicQualityProfile::canonical(CinematicQualityTier::StoryboardSmoke)
        .expect("storyboard profile");
    let wrong = source.replacen(
        &format!(
            "render_budget_profile={}:{}",
            CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION,
            CinematicQualityProfile::canonical(CinematicQualityTier::Final4k)
                .expect("final profile")
                .identity()
                .to_hex()
        ),
        &format!(
            "render_budget_profile={}:{}",
            CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION,
            storyboard.identity().to_hex()
        ),
        1,
    );
    let error = CinematicConfigDocument::from_str(&wrong).expect_err("profile mismatch");
    assert_eq!(
        error,
        CinematicConfigDocumentError::BudgetProfileMismatch {
            field: "render_budget_profile"
        }
    );
}

#[test]
fn stale_or_unavailable_asset_bytes_refuse_at_the_canonical_asset_index() {
    let document = CinematicConfigDocument::from_str(&document_text(
        CinematicQualityTier::StoryboardSmoke,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        "none",
    ))
    .expect("document");

    let stale = document
        .admit_with_asset_resolver(|class, index, declaration| {
            if class == CinematicDocumentAssetClass::Material && index == 0 {
                Ok(b"substituted bytes".to_vec())
            } else {
                resolve(class, index, declaration)
            }
        })
        .expect_err("stale bytes");
    assert!(matches!(
        stale,
        CinematicConfigDocumentError::AssetIdentityMismatch {
            class: CinematicDocumentAssetClass::Material,
            index: 0
        }
    ));

    let unavailable = document
        .admit_with_asset_resolver(|class, index, declaration| {
            if class == CinematicDocumentAssetClass::Environment {
                Err(CinematicAssetAccessError::Unavailable)
            } else {
                resolve(class, index, declaration)
            }
        })
        .expect_err("missing environment");
    assert!(matches!(
        unavailable,
        CinematicConfigDocumentError::AssetAccess {
            class: CinematicDocumentAssetClass::Environment,
            index: 0,
            kind: CinematicAssetAccessError::Unavailable,
        }
    ));
}

#[test]
fn byte_and_line_ceilings_fail_before_unbounded_parsing() {
    let oversized = vec![b'x'; MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES + 1];
    assert!(matches!(
        CinematicConfigDocument::from_bytes(&oversized),
        Err(CinematicConfigDocumentError::DocumentTooLarge { .. })
    ));

    let too_many_lines = "#\n".repeat(MAX_CINEMATIC_CONFIG_DOCUMENT_LINES + 1);
    assert_eq!(
        CinematicConfigDocument::from_str(&too_many_lines),
        Err(CinematicConfigDocumentError::TooManyLines)
    );
}

#[test]
fn mux_request_requires_an_explicit_mux_capability() {
    let mux = format!("av1-opus-matroska:1:{}", id("mux-adapter").to_hex());
    let source = document_text(
        CinematicQualityTier::StoryboardSmoke,
        "assets/steel.spectrum",
        "assets/tungsten.spectrum",
        &mux,
    );
    let without_capability = source.replace(
        "capabilities=render,audio,quarantined-mux",
        "capabilities=render,audio",
    );
    assert_eq!(
        CinematicConfigDocument::from_str(&without_capability),
        Err(CinematicConfigDocumentError::InvalidField {
            field: "capabilities"
        })
    );
}
