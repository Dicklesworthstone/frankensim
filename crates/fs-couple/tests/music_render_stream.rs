//! The actual command must publish the same waveform/hash as direct physical
//! stepping, including a tail shorter than the requested output block.
use std::process::Command;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::encode_pcm16_wav;

#[test]
fn real_command_streams_the_pinned_string_without_changing_samples_or_clobbering_output() {
    let rate = 48_000;
    let samples = 257;
    let wave_speed = (60.0_f64 / 6.0e-4).sqrt();
    let modes = (1..=3).map(|k| ModalAcousticMode {
        angular_frequency_rad_s: f64::from(k) * core::f64::consts::PI * wave_speed / 0.65,
        damping_ratio: 1.0e-3 * f64::from(k),
        pressure_per_modal_velocity: fs_math::c64::C64::new(2.0, 0.0),
    }).collect();
    let mut model = ModalAcousticTimeModel::try_new(rate, modes, ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&(1..=3).map(|k| ModalAcousticState {
        displacement_m_sqrt_kg: 1.0e-3 / f64::from(k), velocity_m_sqrt_kg_per_s: 0.0,
    }).collect::<Vec<_>>()).unwrap();
    let pressure: Vec<_> = (0..samples).map(|_| model.step(&[0.0; 3]).unwrap().observer_pressure_pa).collect();
    let (expected, clips) = encode_pcm16_wav(&pressure, rate, 200.0).unwrap();
    let hash = fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1", &expected).to_hex();
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let directory = std::env::temp_dir().join(format!("frankensim-music-stream-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    let seconds = (f64::from(samples) / f64::from(rate)).to_string();
    for block in [37, 128] {
        let wav = directory.join(format!("string-{block}.wav"));
        let run = || Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("string").arg(&wav).args(["--seconds", &seconds, "--block", &block.to_string()])
            .output().unwrap();
        let output = run();
        assert!(output.status.success(), "stdout={} stderr={}",
            String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
        let sidecar = wav.with_extension("provenance.json");
        let provenance = std::fs::read_to_string(&sidecar).unwrap();
        assert_eq!(std::fs::read(&wav).unwrap(), expected);
        assert!(provenance.contains(&format!("\"wav_blake3\":\"{hash}\"")));
        assert!(provenance.contains("\"samples\":257,"));
        assert!(provenance.contains(&format!("\"clipped_samples\":{clips},")));
        assert!(!run().status.success(), "existing artifacts must refuse");
        assert_eq!(std::fs::read(&wav).unwrap(), expected);
        assert_eq!(std::fs::read_to_string(&sidecar).unwrap(), provenance);
    }
}
