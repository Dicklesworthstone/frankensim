//! File/command tests for the same two-body contact, not a prescribed receiver sound.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::{ControlDelta, GatedRenderOutcome, RenderContext, RenderError, RenderVoice};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig,
};
use fs_couple::render::schedule::force::coupled::contact::{
    ContactModalSystem, ModalContact, ModalContactConfig,
};
use fs_couple::render::schedule::force::coupled::render::contact::ContactModalVoice;
use fs_couple::render::schedule::force::file::{
    ModalPerformance, MODAL_CONTACT_PERFORMANCE_SCHEMA, MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN,
    MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN,
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/modal-contact.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn direct() -> ContactModalSystem {
    let modes = [800.0, 1100.0].into_iter().enumerate().map(|(i, omega)| {
        let mut model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
            angular_frequency_rad_s: omega, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(i as f64, 0.0),
        }], ModalAcousticTimeBudget {
            nyquist_guard_fraction: 0.9, maximum_abs_displacement_m_sqrt_kg: 1.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1000.0,
            maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
        }).unwrap();
        model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: if i == 0 { 1.0 } else { 0.0 } }]).unwrap();
        model
    }).collect();
    let gate = CancelGate::new();
    let network = CoupledModalSystem::new(modes, vec![], ModalCouplingConfig {
        max_modes: 4096, max_connections: 4, max_setup_terms: 1024,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-10,
        energy_relative_tolerance: 1e-8,
    }, &gate).unwrap();
    ContactModalSystem::new(network, ModalContact {
        left: ModalAttachment { component: 0, shapes: vec![1.0] },
        right: ModalAttachment { component: 1, shapes: vec![1.0] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0005], vec![1.0], 3e8, 1.5,
            "authored:two-flexible-bodies".into()).unwrap().with_internal_loss(0.2).unwrap(),
    }, ModalContactConfig { max_iterations: 80, maximum_force_n: 1e5,
        maximum_penetration_m: 0.002, force_absolute_tolerance_n: 1e-10,
        force_relative_tolerance: 1e-10 }, &gate).unwrap()
}
fn reference(samples: usize) -> Vec<f64> {
    let mut system = direct();
    (0..samples).map(|i| system.step(&[if (71..120).contains(&i) { 3.0 } else { 0.0 }, 0.0])
        .unwrap().observer_pressure_pa).collect()
}
fn renderer(input: &str, block: usize) -> ScheduledRenderer {
    ModalPerformance::from_bytes(input.as_bytes(), block).unwrap().into_renderer()
}
fn pressure(input: &str, block: usize) -> Vec<f64> {
    let admitted = ModalPerformance::from_bytes(input.as_bytes(), block).unwrap();
    let mut out = vec![0.0; admitted.info().samples as usize];
    let mut render = admitted.into_renderer();
    for chunk in out.chunks_mut(block) { render.block(chunk).unwrap(); }
    out
}

#[test]
fn file_render_equals_direct_contact_physics_for_every_partition_and_tail() {
    let expected = reference(1201);
    assert!(expected[..20].iter().all(|p| *p == 0.0), "no response before contact");
    assert!(expected.iter().any(|p| p.abs() > 0.1), "only the unforced receiver is observed");
    for block in [1, 7, 37, 64, 512] {
        let actual = pressure(INPUT, block);
        assert_eq!(actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    }
    let admitted = ModalPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap();
    let info = admitted.info();
    assert_eq!(info.schema, MODAL_CONTACT_PERFORMANCE_SCHEMA);
    assert_eq!((info.voices, info.modes, info.connections, info.force_events), (2, 2, 0, 2));
    assert_eq!(info.input_hash, hash_domain(MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()));
    assert_ne!(info.input_hash, hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()));
    let render = admitted.into_renderer();
    assert_eq!(render.pending_controls().iter().map(|c| c.sample).collect::<Vec<_>>(), vec![71, 120]);
    assert!(render.validate_sample_rate(44_100).is_err());
}

#[test]
fn cancellation_preserves_contact_state_and_does_not_repeat_pending_actuator_events() {
    let expected = reference(192);
    let mut render = renderer(INPUT, 64);
    let mut actual = [0.0; 192];
    render.block(&mut actual[..64]).unwrap();
    let gate = CancelGate::new(); gate.request();
    let mut untouched = [-7.0; 64];
    assert_eq!(render.render_under_gate(&gate, &mut untouched, 64, 1).unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(untouched, [-7.0; 64]);
    assert_eq!(render.samples_rendered(), 64);
    assert_eq!(render.pending_controls()[0].sample, 71);
    render.render_under_gate(&CancelGate::new(), &mut actual[64..], 64, 2).unwrap();
    assert_eq!(actual.map(f64::to_bits).as_slice(),
        expected.iter().map(|p| p.to_bits()).collect::<Vec<_>>().as_slice());
    assert_eq!(render.applied_controls().len(), 2);
}

#[test]
fn disabled_or_unreachable_contact_cannot_excite_the_unforced_receiver() {
    for inactive in [
        INPUT.replace("contact 300000000", "contact 0"),
        INPUT.replace("0.0005 1 authored:", "1 1 authored:"),
    ] {
        assert!(pressure(&inactive, 37).iter().all(|p| *p == 0.0));
    }
    let base = pressure(INPUT, 37);
    let lossless = pressure(&INPUT.replace("1.5 0.2 0.0005", "1.5 0 0.0005"), 37);
    assert!(base.iter().zip(&lossless).any(|(a, b)| (a-b).abs() > 1e-3));
}

#[test]
fn contact_records_reject_malformed_shapes_budgets_and_unimplemented_preloads() {
    for (end, _) in INPUT.match_indices('\n') {
        if end+1 < INPUT.len() {
            assert!(ModalPerformance::from_bytes(&INPUT.as_bytes()[..end+1], 37).is_err());
        }
    }
    for invalid in [
        INPUT.replace("contact_left 0 1", "contact_left 2 1"),
        INPUT.replace("contact_right 1 1", "contact_right 1"),
        INPUT.replace("contact_right 1 1", "contact_right 1 1 2"),
        INPUT.replace("contact_limits 80", "contact_limits 129"),
        INPUT.replace("contact_limits 80", "contact_limits 0"),
        INPUT.replace("0.002 1e-10", "-0.002 1e-10"),
        INPUT.replace("1.5 0.2 0.0005", "1.5 -0.2 0.0005"),
        INPUT.replace("contact 300000000", "contact NaN"),
        INPUT.replace("authored:two-flexible-bodies", &"s".repeat(1025)),
        INPUT.replace("voice retain-state", "voice static-preload")
            .replace("mode 800 0 0 0 0 1", "mode 800 0 0 0 0 0"),
        format!("{INPUT}contact 1 1 0 0 1 ignored\n"),
    ] {
        assert!(ModalPerformance::from_bytes(invalid.as_bytes(), 37).is_err(), "{invalid}");
    }
}

#[test]
fn control_batch_refusals_are_atomic_and_physics_failure_poisons_the_callback_host() {
    let host = || RenderContext::new(vec![RenderVoice::ContactModal(Box::new(
        ContactModalVoice::new(direct(), vec![0.0; 2]).unwrap()))], 32);
    let valid = ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 3.0 };
    for invalid in [
        ControlDelta::SetModalForce { voice: 0, mode: 2, force_n_per_sqrt_kg: 1.0 },
        ControlDelta::SetBlowingPressure { voice: 0, pressure_pa: 3.0 },
        ControlDelta::SetPlateForce { voice: 0, force_n: 3.0 },
    ] {
        let mut context = host();
        assert!(context.apply_controls(&[valid, invalid]).is_err());
        assert!(context.control_log().is_empty());
        let mut actual = [0.0; 32]; context.block(&mut actual).unwrap();
        assert_eq!(actual.map(f64::to_bits).as_slice(),
            reference(32).iter().map(|v| v.to_bits()).collect::<Vec<_>>().as_slice());
    }
    let mut context = host();
    context.apply_controls(&[ControlDelta::SetModalForce {
        voice: 0, mode: 0, force_n_per_sqrt_kg: 1e100,
    }]).unwrap();
    assert!(matches!(context.block(&mut [0.0; 32]), Err(RenderError::Coupled(_))));
    let mut untouched = [7.0; 32];
    assert!(matches!(context.block(&mut untouched), Err(RenderError::Poisoned)));
    assert_eq!(untouched, [7.0; 32]);
}

fn directory() -> PathBuf {
    let sequence = NEXT.fetch_add(1, Ordering::Relaxed);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-contact-render-{}-{stamp}-{sequence}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    dir
}
fn run(input: &Path, output: &Path, block: usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(output)
        .arg("--block").arg(block.to_string()).output().unwrap()
}

#[test]
fn command_exports_direct_contact_audio_and_relocatable_identity_without_overwriting() {
    let dir = directory();
    let input = dir.join("contact.performance"); std::fs::write(&input, INPUT).unwrap();
    let (expected, clips) = encode_pcm16_wav(&reference(1201), 48_000, 1.0).unwrap();
    let wav_hash = hash_domain(WAV_DOMAIN, &expected).to_hex();
    let input_hash = hash_domain(MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()).to_hex();
    for block in [37, 512] {
        let output = dir.join(format!("contact-{block}.wav"));
        let result = run(&input, &output, block);
        assert!(result.status.success(), "stdout={} stderr={}",
            String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read(output.with_extension("provenance.json")).unwrap();
        let sidecar = std::str::from_utf8(&metadata).unwrap();
        for value in [format!("\"schema\":\"{MODAL_CONTACT_PERFORMANCE_SCHEMA}\""),
            "\"contacts\":1".into(), "\"contact\":\"implicit-nonadhesive-power-law\"".into(),
            format!("\"wav_blake3\":\"{wav_hash}\""), format!("\"blake3\":\"{input_hash}\""),
            format!("\"clipped_samples\":{clips}")]
        { assert!(sidecar.contains(&value), "missing {value}: {sidecar}"); }
        assert!(!run(&input, &output, block).status.success());
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        assert_eq!(std::fs::read(output.with_extension("provenance.json")).unwrap(), metadata);
    }
    let relocated = dir.join("relocated.performance"); std::fs::write(&relocated, INPUT).unwrap();
    let replay = dir.join("replay.wav");
    assert!(run(&relocated, &replay, 37).status.success());
    assert_eq!(std::fs::read(&replay).unwrap(), expected);
    assert_eq!(std::fs::read(replay.with_extension("provenance.json")).unwrap(),
        std::fs::read(dir.join("contact-37.provenance.json")).unwrap());
    assert_eq!(expected.len(), 44+2*1201);
}

#[test]
fn command_refuses_invalid_contact_before_output_creation_and_renders_the_null_control() {
    let dir = directory(); let input = dir.join("bad.performance"); let output = dir.join("new.wav");
    std::fs::write(&input, INPUT.replace("contact_left 0 1", "contact_left 99 1")).unwrap();
    assert!(!run(&input, &output, 37).status.success());
    assert!(!output.exists()); assert!(!output.with_extension("provenance.json").exists());
    std::fs::write(&input, INPUT.replace("contact 300000000", "contact 0")).unwrap();
    assert!(run(&input, &output, 37).status.success());
    assert!(std::fs::read(&output).unwrap()[44..].iter().all(|b| *b == 0));
}
