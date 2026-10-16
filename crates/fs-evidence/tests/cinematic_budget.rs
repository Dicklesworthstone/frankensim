//! G0/G3 boundary and dry-admission tests for cinematic resource envelopes.

use fs_evidence::cinematic_budget::{
    AovPreset, CinematicBudgetError, CinematicBudgetRepair, CinematicQualityProfile,
    CinematicQualityProfileInput, CinematicQualityTier, CinematicResourceAvailability,
    CinematicResourceDeficit, CinematicResourceKind, DenoisePolicy, ResourceLimitSource,
    admit_cinematic_budget,
};

const GIB: u64 = 1024 * 1024 * 1024;
type Mutation = fn(&mut CinematicQualityProfileInput);

fn abundant() -> CinematicResourceAvailability {
    CinematicResourceAvailability {
        memory_bytes: 128 * GIB,
        free_storage_bytes: 2 * 1024 * GIB,
        wall_time_available_s: 365 * 86_400,
        worker_capacity: 256,
        measured_camera_paths_per_second: 10_000_000,
    }
}

fn shortage_parts(
    error: CinematicBudgetError,
) -> (Vec<CinematicResourceDeficit>, Vec<CinematicBudgetRepair>) {
    match error {
        CinematicBudgetError::InsufficientResources {
            deficits, repairs, ..
        } => (deficits, repairs),
        other => {
            assert_eq!(other.code(), "cinematic-budget-insufficient-resources");
            (Vec::new(), Vec::new())
        }
    }
}

#[test]
fn all_four_frozen_profiles_are_explicit_and_admit_with_abundant_resources() {
    let mut identities = std::collections::BTreeSet::new();
    for tier in [
        CinematicQualityTier::StoryboardSmoke,
        CinematicQualityTier::Daily1080p,
        CinematicQualityTier::Qualification4kFrame,
        CinematicQualityTier::Final4k,
    ] {
        let profile = CinematicQualityProfile::canonical(tier).expect("canonical profile");
        let admitted = admit_cinematic_budget(&profile, abundant()).expect("abundant admission");
        assert!(admitted.estimate().camera_paths > 0);
        assert!(admitted.estimate().total_storage_bytes > 0);
        assert!(admitted.summary_json().contains("\"verdict\":\"admitted\""));
        assert!(identities.insert(profile.identity()));
        assert_eq!(profile.canonical_bytes().len(), 91);
    }
}

#[test]
fn profile_identity_binds_every_budget_and_quality_field() {
    let baseline =
        CinematicQualityProfile::canonical(CinematicQualityTier::Final4k).expect("canonical final");
    let baseline_identity = baseline.identity();
    let mut encoded = String::new();
    for byte in baseline.canonical_bytes() {
        use std::fmt::Write as _;
        let _ = write!(encoded, "{byte:02x}");
    }
    assert_eq!(
        encoded,
        "010003000f0000700800001800000000000000f000000000010000000400001000e8030000020210002000200040004000000000000000020000008051010000000000001a4f00000000000000000040000000000000800c000000",
        "v1 profile layout and explicit enum tags are a compatibility boundary",
    );
    assert_eq!(
        baseline_identity.to_hex(),
        "fad299a8bd59ce85bbedced35212e408ad1be28c4a34a6e21d8085458ac4af04",
        "domain-separated v1 profile identity is a known-answer vector",
    );

    let mutations: Vec<Mutation> = vec![
        |p| p.first_frame = 1,
        |p| p.frame_count -= 1,
        |p| p.spp_floor += 1,
        |p| p.spp_ceiling += 1,
        |p| p.max_path_depth += 1,
        |p| p.adaptive_error_ppm += 1,
        |p| p.shutter_samples += 1,
        |p| p.tile_width += 1,
        |p| p.tile_height += 1,
        |p| p.worker_limit += 1,
        |p| p.checkpoint_cadence_spp += 1,
        |p| p.memory_ceiling_bytes += 1,
        |p| p.per_frame_wall_time_ceiling_s += 1,
        |p| p.sequence_wall_time_ceiling_s += 1,
        |p| p.output_ceiling_bytes += 1,
        |p| p.minimum_free_space_reserve_bytes += 1,
    ];
    for mutate in mutations {
        let mut changed = baseline.input().clone();
        mutate(&mut changed);
        let changed = CinematicQualityProfile::try_new(changed).expect("valid one-field change");
        assert_ne!(baseline_identity, changed.identity());
    }

    let storyboard = CinematicQualityProfile::canonical(CinematicQualityTier::StoryboardSmoke)
        .expect("canonical storyboard");
    for mutate in [
        (|p: &mut CinematicQualityProfileInput| p.width_pixels -= 1) as Mutation,
        |p: &mut CinematicQualityProfileInput| p.height_pixels -= 1,
    ] {
        let mut changed = storyboard.input().clone();
        mutate(&mut changed);
        let changed = CinematicQualityProfile::try_new(changed).expect("valid dimension change");
        assert_ne!(storyboard.identity(), changed.identity());
    }

    assert_eq!(
        baseline.identity(),
        baseline.identity(),
        "identity is stable"
    );
}

#[test]
fn final_profile_cannot_inherit_preview_quality_or_overwrite_raw_master() {
    let final_profile =
        CinematicQualityProfile::canonical(CinematicQualityTier::Final4k).expect("canonical final");
    let mut candidate = final_profile.input().clone();
    candidate.spp_floor = 16;
    assert_eq!(
        CinematicQualityProfile::try_new(candidate),
        Err(CinematicBudgetError::UnsupportedTierCombination(
            CinematicQualityTier::Final4k
        ))
    );

    let mut candidate = final_profile.input().clone();
    candidate.denoise_policy = DenoisePolicy::PreviewBiased;
    assert_eq!(
        CinematicQualityProfile::try_new(candidate),
        Err(CinematicBudgetError::UnsupportedTierCombination(
            CinematicQualityTier::Final4k
        ))
    );

    let mut candidate = final_profile.input().clone();
    candidate.aov_preset = AovPreset::BeautyXyz;
    assert!(matches!(
        CinematicQualityProfile::try_new(candidate),
        Err(CinematicBudgetError::UnsupportedTierCombination(_))
    ));
}

#[test]
fn zero_and_max_plus_one_scalar_boundaries_refuse_before_estimation() {
    let base =
        CinematicQualityProfile::canonical(CinematicQualityTier::Final4k).expect("canonical final");
    let mutations: Vec<Mutation> = vec![
        |p| p.width_pixels = 0,
        |p| p.height_pixels = 7_681,
        |p| p.frame_count = 0,
        |p| p.frame_count = 289,
        |p| {
            p.first_frame = u32::MAX;
            p.frame_count = 2;
        },
        |p| p.spp_floor = 0,
        |p| p.spp_ceiling = 4_097,
        |p| p.max_path_depth = 65,
        |p| p.tile_width = 1_025,
        |p| p.worker_limit = 1_025,
        |p| p.shutter_samples = 0,
        |p| p.shutter_samples = 1_025,
        |p| p.checkpoint_cadence_spp = 0,
    ];
    for mutate in mutations {
        let mut candidate = base.input().clone();
        mutate(&mut candidate);
        assert!(CinematicQualityProfile::try_new(candidate).is_err());
    }
}

#[test]
fn exact_structural_maxima_remain_valid_before_resource_admission() {
    let mut candidate = CinematicQualityProfile::canonical(CinematicQualityTier::Final4k)
        .expect("canonical final")
        .input()
        .clone();
    candidate.frame_count = 288;
    candidate.spp_ceiling = 4_096;
    candidate.max_path_depth = 64;
    candidate.adaptive_error_ppm = 1_000_000;
    candidate.shutter_samples = 1_024;
    candidate.tile_width = 1_024;
    candidate.tile_height = 1_024;
    candidate.worker_limit = 1_024;
    candidate.checkpoint_cadence_spp = 4_096;
    assert!(CinematicQualityProfile::try_new(candidate).is_ok());
}

#[test]
fn exact_resource_boundary_admits_and_max_plus_one_refuses() {
    let profile = CinematicQualityProfile::canonical(CinematicQualityTier::Qualification4kFrame)
        .expect("qualification profile");
    let initial = admit_cinematic_budget(&profile, abundant()).expect("estimate");
    let estimate = initial.estimate();
    let reserve = profile.input().minimum_free_space_reserve_bytes;
    let exact = CinematicResourceAvailability {
        memory_bytes: estimate.live_memory_bytes,
        free_storage_bytes: estimate.total_storage_bytes + reserve,
        wall_time_available_s: estimate.sequence_wall_time_s,
        worker_capacity: profile.input().worker_limit,
        measured_camera_paths_per_second: abundant().measured_camera_paths_per_second,
    };
    assert!(admit_cinematic_budget(&profile, exact).is_ok());

    let one_short = CinematicResourceAvailability {
        memory_bytes: estimate.live_memory_bytes - 1,
        ..exact
    };
    let error = admit_cinematic_budget(&profile, one_short).expect_err("one byte must refuse");
    let (deficits, _) = shortage_parts(error);
    assert!(deficits.iter().any(|deficit| {
        deficit.kind == CinematicResourceKind::LiveMemoryBytes
            && deficit.source == ResourceLimitSource::HostAvailability
            && deficit.required == estimate.live_memory_bytes
            && deficit.available + 1 == deficit.required
    }));
}

#[test]
fn checked_estimate_accounts_for_film_aovs_checkpoints_audio_and_mux() {
    let profile =
        CinematicQualityProfile::canonical(CinematicQualityTier::Final4k).expect("canonical final");
    let estimate = admit_cinematic_budget(&profile, abundant())
        .expect("admit")
        .estimate();
    assert_eq!(estimate.pixels_per_frame, 3_840 * 2_160);
    assert_eq!(estimate.film_bytes, estimate.pixels_per_frame * 3 * 8);
    assert_eq!(estimate.staging_bytes, estimate.film_bytes);
    assert!(estimate.aov_bytes > estimate.film_bytes);
    assert_eq!(estimate.wav_bytes, 480_000 * 2 * 4);
    assert!(estimate.checkpoint_bytes > estimate.film_bytes);
    assert_eq!(
        estimate.image_sequence_bytes,
        estimate.exr_sequence_bytes + estimate.png_sequence_bytes
    );
    assert!(estimate.total_storage_bytes > estimate.image_sequence_bytes);
}

#[test]
fn missing_measurement_and_each_host_shortage_fail_closed_with_ranked_repairs() {
    let profile =
        CinematicQualityProfile::canonical(CinematicQualityTier::Final4k).expect("canonical final");
    let mut unavailable = abundant();
    unavailable.measured_camera_paths_per_second = 0;
    assert_eq!(
        admit_cinematic_budget(&profile, unavailable),
        Err(CinematicBudgetError::MissingThroughputMeasurement)
    );

    let unavailable = CinematicResourceAvailability {
        memory_bytes: 1,
        free_storage_bytes: 1,
        wall_time_available_s: 1,
        worker_capacity: 1,
        measured_camera_paths_per_second: 1,
    };
    let error = admit_cinematic_budget(&profile, unavailable).expect_err("must refuse");
    let (deficits, repairs) = shortage_parts(error);
    assert!(deficits.len() >= 4);
    assert!(repairs.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(repairs.contains(&CinematicBudgetRepair::IncreaseHostMemory));
    assert!(repairs.contains(&CinematicBudgetRepair::IncreaseFreeStorage));
    assert!(repairs.contains(&CinematicBudgetRepair::ExtendWallTime));
    assert!(repairs.contains(&CinematicBudgetRepair::IncreaseWorkerCapacity));
    assert!(!repairs.contains(&CinematicBudgetRepair::LowerPreviewSppWithNewConfiguration));
}

#[test]
fn preview_shortage_may_suggest_explicit_quality_changes_without_applying_them() {
    let profile = CinematicQualityProfile::canonical(CinematicQualityTier::StoryboardSmoke)
        .expect("storyboard");
    let error = admit_cinematic_budget(
        &profile,
        CinematicResourceAvailability {
            memory_bytes: 1,
            free_storage_bytes: 1,
            wall_time_available_s: 1,
            worker_capacity: 1,
            measured_camera_paths_per_second: 1,
        },
    )
    .expect_err("must refuse rather than degrade");
    let (_, repairs) = shortage_parts(error);
    assert!(repairs.contains(&CinematicBudgetRepair::LowerPreviewSppWithNewConfiguration));
    assert!(repairs.contains(&CinematicBudgetRepair::ReducePreviewAovsWithNewConfiguration));
    assert!(repairs.contains(&CinematicBudgetRepair::ShortenRangeWithNewConfiguration));
}
