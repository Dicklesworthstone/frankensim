//! Hunt–Crossley vs hysteresis-island felt bake-off (music bead
//! `frankensim-music-v8-root-3ez8g.5.2`): same strike fixture (the
//! 87zbd hammer–string island), two hammer force laws — (a) the
//! `fs_material::WoolFelt` hysteretic island and (b) fs-contact's
//! EXISTING `NormalPatchLaw::HuntCrossleySphere` (D24: no third HC
//! implementation anywhere; every HC force here comes from
//! `NormalPatchRequest::evaluate`). QoIs across three hammer
//! velocities: contact duration, mode-3/mode-1 spectral tilt, residual
//! strain after the strike, and energy dissipated per strike. The
//! default test validates the COMMITTED receipt at
//! `tests/receipts/piano-felt-vs-hc.bakeoff` and asserts the
//! structural discrimination (HC holds NO residual strain; its
//! tilt-vs-velocity trend is flatter than the felt's). The `--ignored`
//! mint re-measures everything.

use std::collections::BTreeMap;

use fs_blake3::hash_domain;
use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchGeometry,
    NormalPatchIdentity, NormalPatchLaw, NormalPatchRequest,
};
use fs_couple::bakeoff::{BakeoffOutcome, BakeoffReceipt, ContenderResult};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_material::{Uniaxial, WoolFelt};
use fs_tribo::{InputAuthority::SyntheticFixture, InterfaceMedium::Dry, InterfaceSystemRef};

const RATE: u32 = 48_000;
const SUBSTEPS: usize = 64;
const T_FELT: f64 = 8.0e-3;
const A_CONTACT: f64 = 1.0e-4;
const M_HAMMER: f64 = 8.0e-3;
/// HC effective crown radius [m] (the felt crown at rest).
const R_CROWN: f64 = 8.0e-3;
/// HC reduced modulus [Pa] chosen so the elastic force matches the felt
/// secant at 30% strain (the same matched-fixture rule the 87zbd island
/// used for its linear-spring control).
const HC_MODULUS: f64 = 2.2e6;
/// HC dissipation coefficient [s/m], authored at the standard
/// restitution scale `c ~ 3(1-e)/(2 v_ref)` for e ~ 0.75 at the 2 m/s
/// reference strike. (The first mint used c = 8 s/m — a 33x loading
/// multiplier at forte that inverted the contact-duration trend; the
/// executed lesson that HC's linear-in-rate factor needs restitution-
/// scale coefficients to stay physical.)
const HC_DISSIPATION: f64 = 0.19;

fn receipt_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/receipts/piano-felt-vs-hc.bakeoff")
}

fn string_model() -> ModalAcousticTimeModel {
    let wave_speed = (60.0f64 / 6.0e-4).sqrt();
    let modes = (1..=3)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: k as f64 * core::f64::consts::PI * wave_speed / 0.65,
            damping_ratio: 5.0e-4 * k as f64,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect();
    ModalAcousticTimeModel::try_new(RATE, modes, ModalAcousticTimeBudget::audible_reference())
        .expect("string modes admit")
}

fn felt() -> WoolFelt {
    WoolFelt::new(4.0e5, 0.2, 2.5, 3.2, 0.25, 0.8).expect("felt admits")
}

/// The Hunt–Crossley force through fs-contact's EXISTING law (D24: the
/// only HC evaluation in this file).
fn hc_force(overlap_m: f64, rate_m_per_s: f64) -> f64 {
    if overlap_m <= 0.0 {
        return 0.0;
    }
    let request = NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: "bakeoff-felt-hc".into(),
            source_id: "authored-fixture".into(),
            state_id: "strike".into(),
        },
        interface: InterfaceSystemRef::new("felt-hc", "hammer", "string", SyntheticFixture, Dry)
            .expect("interface"),
        law: NormalPatchLaw::HuntCrossleySphere {
            effective_radius_m: R_CROWN,
            reduced_modulus_pa: HC_MODULUS,
            dissipation_s_per_m: HC_DISSIPATION,
        },
        geometry: NormalPatchGeometry::SpherePlane,
        indentation_m: overlap_m,
        indentation_rate_m_per_s: rate_m_per_s,
        step_s: 1.0 / (f64::from(RATE) * SUBSTEPS as f64),
        line_load_n_per_m: 0.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 10.0,
            layer_thickness_m: 10.0,
            yield_strength_pa: 1.0e9,
            characteristic_rate_m_per_s: 100.0,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 10.0,
            max_strain: 10.0,
            max_patch_to_depth: 10.0,
            max_patch_to_layer: 10.0,
            max_pressure_to_yield: 10.0,
            max_rate_ratio: 10.0,
            min_temperature_k: 200.0,
            max_temperature_k: 400.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.0,
            modulus_relative: 0.0,
            load_relative: 0.0,
        },
    };
    match request.evaluate() {
        Ok(fs_contact::normal_patch::NormalPatchReceipt::Point(receipt)) => {
            receipt.normal_force_n.max(0.0)
        }
        _ => 0.0,
    }
}

enum Hammer {
    Felt {
        law: WoolFelt,
        state: <WoolFelt as Uniaxial>::State,
    },
    HuntCrossley,
}

struct StrikeOutcome {
    contact_samples: usize,
    tilt_3_over_1: f64,
    residual_strain: f64,
    dissipated_j: f64,
}

/// The 87zbd strike island with a pluggable hammer force law.
fn strike(mut hammer: Hammer, v0: f64) -> StrikeOutcome {
    let mut string = string_model();
    let phi: Vec<f64> = (1..=3)
        .map(|k| fs_math::det::sin(k as f64 * core::f64::consts::PI * 0.12))
        .collect();
    let dt = 1.0 / f64::from(RATE);
    let h = dt / SUBSTEPS as f64;
    let mut y_hammer = 0.0f64;
    let mut v_hammer = v0;
    let ke_in = 0.5 * M_HAMMER * v0 * v0;
    let mut work_into_string = 0.0f64;
    let mut contact_samples = 0usize;
    let mut frames = 0usize;
    loop {
        let string_disp: f64 = string
            .states()
            .iter()
            .zip(&phi)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum();
        let mut force_avg = 0.0f64;
        let mut peak_overlap = 0.0f64;
        for _ in 0..SUBSTEPS {
            let overlap = y_hammer - string_disp;
            let f = match &hammer {
                Hammer::Felt { law, state } => {
                    if overlap <= 0.0 {
                        0.0
                    } else {
                        let strain = (overlap / T_FELT).min(0.65);
                        A_CONTACT * law.stress(strain, state)
                    }
                }
                Hammer::HuntCrossley => hc_force(overlap, v_hammer),
            };
            force_avg += f;
            peak_overlap = peak_overlap.max(overlap);
            v_hammer -= h * f / M_HAMMER;
            y_hammer += h * v_hammer;
        }
        force_avg /= SUBSTEPS as f64;
        if let Hammer::Felt { law, state } = &mut hammer {
            let strain = (peak_overlap.max(y_hammer - string_disp).max(0.0) / T_FELT).min(0.65);
            *state = law.update_state(strain, state);
        }
        if force_avg > 0.0 {
            contact_samples += 1;
        }
        let generalized: Vec<f64> = phi.iter().map(|p| force_avg * p).collect();
        let frame = string.step(&generalized).expect("string step");
        work_into_string += frame.input_work_j;
        frames += 1;
        let separated = v_hammer < 0.0 && (y_hammer - string_disp) < -1.0e-4;
        if separated || frames > 4 * RATE as usize / 100 {
            let modal = frame.modal_energy_j.clone();
            let ke_out = 0.5 * M_HAMMER * v_hammer * v_hammer;
            let dissipated_j = ke_in - ke_out - work_into_string;
            let residual_strain = match &hammer {
                Hammer::Felt { law, state } => law.eps_residual(state),
                Hammer::HuntCrossley => 0.0,
            };
            return StrikeOutcome {
                contact_samples,
                tilt_3_over_1: modal[2] / modal[0].max(f64::MIN_POSITIVE),
                residual_strain,
                dissipated_j,
            };
        }
    }
}

fn shared_cards_digest() -> fs_blake3::ContentHash {
    let mut bytes = Vec::new();
    for v in [
        4.0e5f64,
        0.2,
        2.5,
        3.2,
        0.25,
        0.8, // felt
        R_CROWN,
        HC_MODULUS,
        HC_DISSIPATION, // HC
        T_FELT,
        A_CONTACT,
        M_HAMMER, // island
    ] {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    hash_domain("org.frankensim.fs-couple.felt-bakeoff-cards.v1", &bytes)
}

fn measure(kind: &str) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for &v0 in &[0.8f64, 2.0, 4.0] {
        let hammer = match kind {
            "felt" => Hammer::Felt {
                law: felt(),
                state: felt().initial_state(),
            },
            _ => Hammer::HuntCrossley,
        };
        let s = strike(hammer, v0);
        let tag = format!("v{v0:.1}");
        out.insert(format!("contact-samples-{tag}"), s.contact_samples as f64);
        out.insert(format!("tilt-3-1-{tag}"), s.tilt_3_over_1);
        out.insert(format!("residual-strain-{tag}"), s.residual_strain);
        out.insert(format!("dissipated-j-{tag}"), s.dissipated_j);
    }
    out
}

#[test]
#[ignore = "minting run: re-measures both hammers and writes fresh receipt bytes"]
fn mint_felt_bakeoff_receipt() {
    let felt_m = measure("felt");
    let hc_m = measure("hc");
    // Reference = the hysteretic island (the physics-stack contender);
    // residuals then read as each image's distance from the felt stack.
    let receipt = BakeoffReceipt {
        filling: "hammer-felt".to_string(),
        fixture: "crates/fs-couple/tests/bakeoff_felt.rs 87zbd strike island, v0 in \
                  {0.8, 2.0, 4.0} m/s, 3-mode exact-ZOH string"
            .to_string(),
        shared_cards: shared_cards_digest(),
        reference: felt_m.clone(),
        contenders: [
            ContenderResult {
                image: "hysteresis-island".to_string(),
                owner_crates: vec!["fs-material".to_string(), "fs-couple".to_string()],
                measured: felt_m,
                states: 2,
                steps: 3 * SUBSTEPS * 1920,
                solver_iterations: 0,
                failure_modes: vec![
                    "rate-independent envelope (rate half lives in the Prony chain, a \
                     separate island)"
                        .to_string(),
                ],
            },
            ContenderResult {
                image: "hunt-crossley".to_string(),
                owner_crates: vec!["fs-contact".to_string()],
                measured: hc_m,
                states: 0,
                steps: 3 * SUBSTEPS * 1920,
                solver_iterations: 0,
                failure_modes: vec![
                    "structurally CANNOT hold residual crush (memoryless law: residual \
                     strain identically zero at every velocity)"
                        .to_string(),
                    "1 + c*rate factor can go negative at fast rebound (clamped to zero \
                     force here)"
                        .to_string(),
                ],
            },
        ],
        outcome: BakeoffOutcome::KeepForSubset {
            narrowed: "hunt-crossley".to_string(),
            subset: "debug/second-hammer image and generic transient point contacts; never \
                     the piano felt claim (no residual crush, flatter tilt trend)"
                .to_string(),
        },
        rationale: "the hysteretic island is the felt authority (residual crush + \
                    velocity-steepened tilt are the piano mechanisms; HC's tilt trend is \
                    flatter and its residual strain is identically zero); HC stays on the \
                    menu for the claims it still owns (D21)"
            .to_string(),
        listening_receipts: Vec::new(),
    };
    receipt.validate().expect("receipt validates");
    let bytes = receipt.to_canonical_bytes().expect("encode");
    std::fs::write(receipt_path(), &bytes).expect("write receipt");
    println!(
        "minted {} ({} bytes), hash {}",
        receipt_path().display(),
        bytes.len(),
        receipt.content_hash().expect("hash").to_hex()
    );
}

#[test]
fn committed_felt_bakeoff_receipt_shows_the_structural_gap() {
    let bytes = std::fs::read(receipt_path())
        .expect("tests/receipts/piano-felt-vs-hc.bakeoff must be committed (mint test)");
    let receipt = BakeoffReceipt::from_canonical_bytes(&bytes).expect("receipt decodes");
    receipt.validate().expect("receipt validates");
    assert_eq!(receipt.filling, "hammer-felt");
    let felt_m = &receipt.contenders[0].measured;
    let hc_m = &receipt.contenders[1].measured;
    // Structural gap 1: HC holds NO residual strain at any velocity; the
    // felt's residual GROWS with velocity.
    for tag in ["v0.8", "v2.0", "v4.0"] {
        assert_eq!(
            hc_m[&format!("residual-strain-{tag}")],
            0.0,
            "HC must be memoryless"
        );
    }
    let felt_res_soft = felt_m["residual-strain-v0.8"];
    let felt_res_forte = felt_m["residual-strain-v4.0"];
    assert!(
        felt_res_forte > felt_res_soft && felt_res_forte > 0.05,
        "felt residual must grow with velocity ({felt_res_soft} -> {felt_res_forte})"
    );
    // Structural gap 2: the felt's tilt-vs-velocity trend is steeper.
    let trend = |m: &BTreeMap<String, f64>| m["tilt-3-1-v4.0"] / m["tilt-3-1-v0.8"].max(1e-300);
    let felt_trend = trend(felt_m);
    let hc_trend = trend(hc_m);
    assert!(
        felt_trend > 1.2 * hc_trend,
        "felt tilt trend {felt_trend:.3} must exceed HC's {hc_trend:.3} by a real margin"
    );
    // Both contenders behave like hammers: contact shortens with velocity.
    for m in [felt_m, hc_m] {
        assert!(
            m["contact-samples-v4.0"] < m["contact-samples-v0.8"],
            "contact duration must fall with velocity for both laws"
        );
    }
    // The verdict narrows HC, never deletes it.
    let narrowed_ok = matches!(
        &receipt.outcome,
        BakeoffOutcome::KeepForSubset { narrowed, .. } if narrowed == "hunt-crossley"
    );
    assert!(
        narrowed_ok,
        "outcome must narrow hunt-crossley: {:?}",
        receipt.outcome
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-bakeoff-receipt\",\"verdict\":\"pass\",\
         \"felt_residual_forte\":{felt_res_forte:.4},\"felt_tilt_trend\":{felt_trend:.3},\
         \"hc_tilt_trend\":{hc_trend:.3},\"hash\":\"{}\"}}",
        receipt.content_hash().expect("hash").to_hex()
    );
}
