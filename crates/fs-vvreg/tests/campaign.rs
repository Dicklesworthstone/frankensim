//! Level-E campaign manifest and operating matrix test suite
//! (bead `frankensim-extreal-program-f85xj.4.5.4`).
//!
//! Validates:
//! - Complete operating matrix definition with adversarial variations (vent, fan, TIM)
//! - Mandatory pre-flight safety checklist sign-off
//! - Enforced candidate blind-holdout partition
//! - Explicit `NoData` disposition prior to physical hardware availability
//! - Deterministic content hashing and mutant sensitivity

use fs_vvreg::campaign::{
    CampaignDisposition, FanState, LevelECampaignManifest, OperatingMatrixPoint, SafetyChecklist,
    TimCondition, VentConfig,
};
use fs_vvreg::corpus::DatasetPartition;
use fs_vvreg::rig::RigError;

fn valid_safety_checklist() -> SafetyChecklist {
    SafetyChecklist {
        operator_id: "lead-thermo-engineer".to_string(),
        overtemp_cutoff_k: 363.15, // 90 °C
        emergency_stop_verified: true,
        wiring_verified: true,
        calibrations_verified: true,
        signoff_date: "2026-06-01".to_string(),
    }
}

fn valid_operating_matrix() -> Vec<OperatingMatrixPoint> {
    vec![
        OperatingMatrixPoint {
            point_id: "matrix-01-nominal".to_string(),
            heater_power_w: 120.0,
            flow_rate_m3_s: 0.01,
            vent: VentConfig::NominalOpen,
            fan: FanState::Nominal100,
            tim: TimCondition::StandardGrease,
            partition: DatasetPartition::Training,
        },
        OperatingMatrixPoint {
            point_id: "matrix-02-restricted-vent".to_string(),
            heater_power_w: 120.0,
            flow_rate_m3_s: 0.005,
            vent: VentConfig::Restricted50,
            fan: FanState::Nominal100,
            tim: TimCondition::StandardGrease,
            partition: DatasetPartition::Validation,
        },
        OperatingMatrixPoint {
            point_id: "matrix-03-adversarial-blocked".to_string(),
            heater_power_w: 80.0,
            flow_rate_m3_s: 0.001,
            vent: VentConfig::Blocked100,
            fan: FanState::Reduced50,
            tim: TimCondition::DegradedPad,
            partition: DatasetPartition::Validation,
        },
        OperatingMatrixPoint {
            point_id: "matrix-04-blind-holdout".to_string(),
            heater_power_w: 150.0,
            flow_rate_m3_s: 0.012,
            vent: VentConfig::NominalOpen,
            fan: FanState::Nominal100,
            tim: TimCondition::StandardGrease,
            partition: DatasetPartition::BlindHoldout,
        },
    ]
}

#[test]
fn valid_campaign_manifest_with_no_data_validates_and_hashes() {
    let manifest = LevelECampaignManifest {
        campaign_id: "level-e-cooling-campaign-v1".to_string(),
        rig_id: "inhouse-thermal-rig-01".to_string(),
        matrix: valid_operating_matrix(),
        safety: valid_safety_checklist(),
        disposition: CampaignDisposition::NoData {
            reason: "Awaiting physical rig hardware assembly and metrology sign-off".to_string(),
        },
    };

    assert_eq!(manifest.validate(), Ok(()));
    let hash1 = manifest.content_hash();
    let hash2 = manifest.content_hash();
    assert_eq!(hash1, hash2);
    assert!(!hash1.to_hex().is_empty());
}

#[test]
fn matrix_lacking_blind_holdout_refuses_validation() {
    let training_only_matrix = vec![OperatingMatrixPoint {
        point_id: "matrix-01-nominal".to_string(),
        heater_power_w: 120.0,
        flow_rate_m3_s: 0.01,
        vent: VentConfig::NominalOpen,
        fan: FanState::Nominal100,
        tim: TimCondition::StandardGrease,
        partition: DatasetPartition::Training,
    }];

    let manifest = LevelECampaignManifest {
        campaign_id: "level-e-cooling-campaign-v1".to_string(),
        rig_id: "inhouse-thermal-rig-01".to_string(),
        matrix: training_only_matrix,
        safety: valid_safety_checklist(),
        disposition: CampaignDisposition::NoData {
            reason: "Planning phase".to_string(),
        },
    };

    assert!(matches!(
        manifest.validate(),
        Err(RigError::InvalidScalar {
            field: "matrix.partition",
            ..
        })
    ));
}

#[test]
fn unverified_safety_checklist_refuses_validation() {
    let mut bad_safety = valid_safety_checklist();
    bad_safety.emergency_stop_verified = false;

    let manifest = LevelECampaignManifest {
        campaign_id: "level-e-cooling-campaign-v1".to_string(),
        rig_id: "inhouse-thermal-rig-01".to_string(),
        matrix: valid_operating_matrix(),
        safety: bad_safety,
        disposition: CampaignDisposition::NoData {
            reason: "Planning phase".to_string(),
        },
    };

    assert!(matches!(
        manifest.validate(),
        Err(RigError::InvalidScalar {
            field: "safety",
            ..
        })
    ));
}

#[test]
fn negative_heater_power_refuses_validation() {
    let mut bad_matrix = valid_operating_matrix();
    bad_matrix[0].heater_power_w = -10.0;

    let manifest = LevelECampaignManifest {
        campaign_id: "level-e-cooling-campaign-v1".to_string(),
        rig_id: "inhouse-thermal-rig-01".to_string(),
        matrix: bad_matrix,
        safety: valid_safety_checklist(),
        disposition: CampaignDisposition::NoData {
            reason: "Planning phase".to_string(),
        },
    };

    assert!(matches!(
        manifest.validate(),
        Err(RigError::InvalidScalar {
            field: "matrix.heater_power_w",
            ..
        })
    ));
}

#[test]
fn content_hash_is_sensitive_to_matrix_and_safety_mutants() {
    let base = LevelECampaignManifest {
        campaign_id: "level-e-cooling-campaign-v1".to_string(),
        rig_id: "inhouse-thermal-rig-01".to_string(),
        matrix: valid_operating_matrix(),
        safety: valid_safety_checklist(),
        disposition: CampaignDisposition::NoData {
            reason: "Planning phase".to_string(),
        },
    };
    let base_hash = base.content_hash();

    // Mutate fan condition on one point
    let mut mutant1 = base.clone();
    mutant1.matrix[0].fan = FanState::Stalled;
    assert_ne!(base_hash, mutant1.content_hash());

    // Mutate overtemp cutoff
    let mut mutant2 = base.clone();
    mutant2.safety.overtemp_cutoff_k = 373.15;
    assert_ne!(base_hash, mutant2.content_hash());

    // Mutate vent configuration
    let mut mutant3 = base;
    mutant3.matrix[0].vent = VentConfig::Blocked100;
    assert_ne!(base_hash, mutant3.content_hash());
}
