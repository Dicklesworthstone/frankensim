//! Friction input -> real modal network -> scheduler -> existing command/encoder.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::{ControlDelta, GatedRenderOutcome, RenderError};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_couple::render::schedule::force::coupled::contact::multiple::friction::ModalFriction;
use fs_couple::render::schedule::force::file::{ModalPerformance, MODAL_FRICTION_PERFORMANCE_SCHEMA,
    MODAL_FRICTION_PERFORMANCE_HASH_DOMAIN, MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA,
    MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/modal-friction.performance");
const SAMPLES: usize = 601;
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn direct() -> MultiContactModalSystem {
    let models = [1.0, 0.0, -0.3].into_iter().enumerate().map(|(i, velocity)| {
        let modes = [800.0, 600.0].into_iter().enumerate().map(|(mode, omega)| ModalAcousticMode {
            angular_frequency_rad_s: omega, damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(if i == 1 && mode == 1 { 1.0 } else { 0.0 }, 0.0),
        }).collect();
        let mut model = ModalAcousticTimeModel::try_new(48000, modes, ModalAcousticTimeBudget {
            nyquist_guard_fraction: 0.9, maximum_abs_displacement_m_sqrt_kg: 1.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1000.0, maximum_total_energy_j: 1000.0,
            maximum_abs_pressure_pa: 1000.0,
        }).unwrap();
        model.restore_states(&[ModalAcousticState::default(), ModalAcousticState {
            displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: velocity,
        }]).unwrap();
        model
    }).collect();
    let network = CoupledModalSystem::new(models, vec![], ModalCouplingConfig {
        max_modes: 4096, max_connections: 4, max_setup_terms: 4096, nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1000.0, maximum_abs_pressure_pa: 1000.0,
        maximum_abs_connection_force_n: 1e6, solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8,
    }, &CancelGate::new()).unwrap();
    let contacts = [(0, 1), (1, 2)].into_iter().map(|(a, b)| (ModalContact {
        left: ModalAttachment { component: a, shapes: vec![1.0, 0.0] },
        right: ModalAttachment { component: b, shapes: vec![1.0, 0.0] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![-0.0006], vec![1.0], 1e5, 1.0,
            "synthetic direct comparison".into()).unwrap(),
    }, ModalContactConfig { max_iterations: 80, maximum_force_n: 1e5,
        maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-10, force_relative_tolerance: 1e-10 })).collect();
    let friction = [(0.3, 0.4), (0.5, -0.2)].into_iter().map(|(coefficient, mixed)| Some(ModalFriction {
        left_shapes: vec![mixed, 1.0], right_shapes: vec![mixed, 1.0], coefficient,
        regularization_speed_m_s: 0.01, maximum_force_n: 1e4, source: "synthetic direct comparison".into(),
    })).collect();
    MultiContactModalSystem::new(network, contacts, MultiContactConfig {
        max_contacts: 4, max_sweeps: 128, max_setup_terms: 4096,
    }, &CancelGate::new()).unwrap().with_friction(friction, &CancelGate::new()).unwrap()
}
fn renderer(text: &str, block: usize) -> ScheduledRenderer {
    ModalPerformance::from_bytes(text.as_bytes(), block).unwrap().into_renderer()
}
fn pressure(text: &str, block: usize) -> Vec<f64> {
    let mut runtime = renderer(text, block);
    let mut out = vec![0.0; SAMPLES];
    for chunk in out.chunks_mut(block) { runtime.block(chunk).unwrap(); }
    assert!(runtime.pending_controls().is_empty());
    out
}
fn replace_friction(records: &str, schema: &str) -> String {
    let (before, rest) = INPUT.split_once("frictions 2\n").unwrap();
    let (_, events) = rest.split_once("events 4\n").unwrap();
    format!("{}{records}events 4\n{events}", before.replace(MODAL_FRICTION_PERFORMANCE_SCHEMA, schema))
}
fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-friction-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap(); path
}
fn run(input: &Path, out: &Path, block: usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(out)
        .arg("--block").arg(block.to_string()).output().unwrap()
}
fn succeeded(result: &Output) {
    assert!(result.status.success(), "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
}

#[test]
fn v5_matches_direct_joint_friction_with_short_tails_and_all_callback_partitions() {
    let mut system = direct();
    let mut expected = Vec::new();
    let mut simultaneous = false;
    let mut loss = 0.0;
    for sample in 0..SAMPLES {
        let external = [0.0, if (70..110).contains(&sample) { 0.5 } else { 0.0 }, 0.0, 0.0,
            0.0, if (200..260).contains(&sample) { -0.4 } else { 0.0 }];
        let frame = system.step(&external).unwrap();
        expected.push(frame.observer_pressure_pa.to_bits());
        simultaneous |= frame.contacts.iter().all(|c| c.normal_force_n > 1.0)
            && frame.friction.iter().all(|c| c.unwrap().reaction_n.abs() > 0.1);
        loss += frame.friction_dissipation_j;
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert!(simultaneous && loss > 0.01);
    for block in [1, 7, 37, 64, 512, SAMPLES] {
        assert_eq!(pressure(INPUT, block).iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
    }
    let info = ModalPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap().info();
    assert_eq!(info.schema, MODAL_FRICTION_PERFORMANCE_SCHEMA);
    assert_eq!((info.voices, info.modes, info.contacts, info.friction_contacts, info.connections, info.force_events),
        (3, 6, 2, 2, 0, 4));
    assert_eq!(info.input_hash, hash_domain(MODAL_FRICTION_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()));
}

#[test]
fn explicit_none_and_zero_friction_preserve_normal_only_motion_and_schema_identity() {
    let none = replace_friction("frictions 2\nfriction none\nfriction none\n", MODAL_FRICTION_PERFORMANCE_SCHEMA);
    let v4 = replace_friction("", MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA);
    let zero = INPUT.replace("regularized-coulomb 0.3", "regularized-coulomb 0")
        .replace("regularized-coulomb 0.5", "regularized-coulomb 0");
    let expected = pressure(&v4, 37);
    assert!(expected.iter().all(|p| *p == 0.0));
    assert_eq!(pressure(&none, 37), expected);
    assert_eq!(pressure(&zero, 37), expected);
    assert!(pressure(INPUT, 37).iter().any(|p| p.abs() > 0.01));
    let info = |text: &str| ModalPerformance::from_bytes(text.as_bytes(), 37).unwrap().info();
    assert_eq!(info(&none).friction_contacts, 0);
    assert_eq!(info(&zero).friction_contacts, 2);
    assert_eq!(info(&v4).input_hash, hash_domain(MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN, v4.as_bytes()));
    for changed in [INPUT.replace("synthetic-left-friction", "different-source"),
        INPUT.replace("regularized-coulomb 0.3", "regularized-coulomb 0.31")] {
        assert_ne!(info(&changed).input_hash, info(INPUT).input_hash);
    }
    let left_none = INPUT.replace(
        "friction regularized-coulomb 0.3 0.01 10000 synthetic-left-friction\nfriction_left 0 0.4 1\nfriction_right 1 0.4 1\n",
        "friction none\n");
    assert_eq!(info(&left_none).friction_contacts, 1);
}

#[test]
fn malformed_or_incomplete_friction_sets_never_receive_defaults() {
    for text in [
        INPUT.replace("frictions 2", "frictions 0"), INPUT.replace("frictions 2", "frictions 3"),
        INPUT.replace("frictions 2", "frictions 18446744073709551615"),
        INPUT.replace("regularized-coulomb 0.3", "coulomb 0.3"),
        INPUT.replace("regularized-coulomb 0.3", "regularized-coulomb -0.3"),
        INPUT.replace("regularized-coulomb 0.3", "regularized-coulomb NaN"),
        INPUT.replace("0.3 0.01 10000", "0.3 0 10000"),
        INPUT.replace("0.3 0.01 10000", "0.3 0.01 0"),
        INPUT.replace("friction_left 0 0.4 1", "friction_left 2 0.4 1"),
        INPUT.replace("friction_right 2 -0.2 1", "friction_right 0 -0.2 1"),
        INPUT.replace("friction_left 0 0.4 1", "friction_left 999 0.4 1"),
        INPUT.replace("friction_left 0 0.4 1", "friction_left 0 0.4"),
        INPUT.replace("friction_left 0 0.4 1", "friction_left 0 0.4 1 0"),
        INPUT.replace("friction_left 0 0.4 1", "friction_left 0 NaN 1"),
        INPUT.replace("synthetic-left-friction", &"a".repeat(1025)),
        INPUT.replace("multi_contact_limits 4 128 4096", "multi_contact_limits 4 128 100"),
        INPUT.replace(MODAL_FRICTION_PERFORMANCE_SCHEMA, MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA),
        INPUT.replacen("voice retain-state 2 1", "voice static-preload 2 1", 1),
        replace_friction("frictions 2\nfriction none ignored\nfriction none\n", MODAL_FRICTION_PERFORMANCE_SCHEMA),
        format!("{INPUT}ignored\n"),
    ] { assert!(ModalPerformance::from_bytes(text.as_bytes(), 37).is_err(), "accepted {text}"); }
    for (end, _) in INPUT.match_indices('\n') {
        if end+1 < INPUT.len() { assert!(ModalPerformance::from_bytes(&INPUT.as_bytes()[..end+1], 37).is_err()); }
    }
}

#[test]
fn cancelled_friction_callbacks_resume_without_replaying_or_dropping_forces() {
    let expected = pressure(INPUT, 37);
    let mut runtime = renderer(INPUT, 37);
    let mut output = vec![0.0; SAMPLES];
    runtime.block(&mut output[..37]).unwrap();
    let pending = runtime.pending_controls().to_vec();
    let gate = CancelGate::new(); gate.request();
    let mut untouched = [7.0; 37];
    assert_eq!(runtime.render_under_gate(&gate, &mut untouched, 37, 1).unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 });
    assert_eq!(untouched, [7.0; 37]); assert_eq!(runtime.samples_rendered(), 37);
    assert_eq!(runtime.pending_controls(), pending.as_slice());
    for chunk in output[37..].chunks_mut(37) { runtime.block(chunk).unwrap(); }
    assert_eq!(output, expected); assert!(runtime.pending_controls().is_empty());
}

#[test]
fn friction_hosts_keep_atomic_controls_and_poison_after_a_refused_callback() {
    let make = || renderer(INPUT, 32).into_context();
    let mut baseline = make(); let mut expected = [0.0; 32]; baseline.block(&mut expected).unwrap();
    let valid = ControlDelta::SetModalForce { voice: 0, mode: 1, force_n_per_sqrt_kg: 100.0 };
    let invalid = ControlDelta::SetModalForce { voice: 0, mode: 6, force_n_per_sqrt_kg: 1.0 };
    let mut context = make(); assert!(context.apply_controls(&[valid, invalid]).is_err());
    assert!(context.control_log().is_empty());
    let mut got = [0.0; 32]; context.block(&mut got).unwrap(); assert_eq!(got, expected);
    let capped = INPUT.replace("0.3 0.01 10000", "0.3 0.01 0.000000000001");
    let mut context = renderer(&capped, 32).into_context();
    assert!(context.block(&mut [0.0; 32]).is_err());
    let mut untouched = [7.0; 32];
    assert!(matches!(context.block(&mut untouched), Err(RenderError::Poisoned)));
    assert_eq!(untouched, [7.0; 32]);
}

#[test]
fn command_retains_v5_identity_exact_wav_short_tail_and_nonoverwrite() {
    let dir = directory(); let input = dir.join("source.performance"); std::fs::write(&input, INPUT).unwrap();
    let wave = pressure(INPUT, 37);
    let (expected, clipped) = encode_pcm16_wav(&wave, 48000, 2.0).unwrap(); assert_eq!(clipped, 0);
    let wav_hash = hash_domain(WAV_DOMAIN, &expected).to_hex();
    let source_hash = hash_domain(MODAL_FRICTION_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()).to_hex();
    for block in [37, 512] {
        let out = dir.join(format!("block-{block}.wav")); succeeded(&run(&input, &out, block));
        let bytes = std::fs::read(&out).unwrap(); assert_eq!(bytes, expected);
        assert_eq!(bytes.len(), 44+2*SAMPLES);
        assert_eq!(u32::from_le_bytes(bytes[40..44].try_into().unwrap()), 1202);
        let metadata = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        for token in [format!("\"schema\":\"{MODAL_FRICTION_PERFORMANCE_SCHEMA}\""),
            "\"contacts\":2".to_string(), "\"friction_contacts\":2".to_string(),
            "\"friction\":\"simultaneous-1d-regularized-coulomb\"".to_string(),
            format!("\"wav_blake3\":\"{wav_hash}\""), format!("\"blake3\":\"{source_hash}\"")] {
            assert!(metadata.contains(&token), "missing {token}: {metadata}");
        }
    }
    let relocated = dir.join("relocated.performance"); std::fs::write(&relocated, INPUT).unwrap();
    let out = dir.join("replay.wav"); succeeded(&run(&relocated, &out, 37));
    assert_eq!(std::fs::read(&out).unwrap(), expected);
    let metadata = std::fs::read(out.with_extension("provenance.json")).unwrap();
    assert_eq!(metadata, std::fs::read(dir.join("block-37.provenance.json")).unwrap());
    assert!(!run(&input, &out, 37).status.success());
    assert_eq!(std::fs::read(&out).unwrap(), expected);
    assert_eq!(std::fs::read(out.with_extension("provenance.json")).unwrap(), metadata);
    std::fs::write(&input, INPUT.replace("frictions 2", "frictions 1")).unwrap();
    let refused = dir.join("refused.wav"); assert!(!run(&input, &refused, 37).status.success());
    assert!(!refused.exists() && !refused.with_extension("provenance.json").exists());
}
