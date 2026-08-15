//! Music-lane render CLI (bead `frankensim-music-t-out-render-ib15w`,
//! program root `frankensim-music-v8-root-3ez8g`): the first binary that
//! puts music-lane audio ON DISK.
//!
//! ```text
//! music_render <fixture> <out.wav> [--seconds S] [--block N] [--full-scale-pa P]
//! ```
//!
//! Fixtures are PINNED compositions of gated machinery (`reed`: the
//! massless-reed 2.2 mm characteristic-line voice; `string`: a plucked
//! three-mode exact-ZOH modal string). Rendering goes through the block
//! render API (`fs_couple::render`) — the same path the budget lane
//! measures — and encoding goes through the ONE pascals→PCM owner,
//! `fs_couple::pcm_wav::encode_pcm16_wav`: mono PCM16, physically scaled
//! by the declared full-scale, NEVER peak-normalized (normalization would
//! hide a material or temperature change), clips COUNTED and reported,
//! never rewritten.
//!
//! Seam decision (recorded here and on beads ib15w + h7xu5.7.8): the
//! music lane's encoder owner is `fs_couple::pcm_wav`; the cinematic
//! stack's receipt-hashed stereo encoder stays cinematic; the 7.8 adapter
//! consumes music-lane pressure through its own admission contract when
//! it lands. No third RIFF writer exists in the music lane.
//!
//! Output contract (campaign-script discipline): the output path and its
//! `.provenance.json` sidecar are REFUSED if they already exist; every
//! emitted JSON-lines field is deterministic — same arguments produce
//! bit-identical WAV bytes and provenance on any single host (the WAV
//! content hash is the replay check). Sample rate is pinned at 48 kHz to
//! keep the ecosystem coherent (fs-psycho refuses other rates).

use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::{ModalStringVoice, ReedBoreVoice, RenderContext, RenderVoice};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;

const RATE: u32 = 48_000;
const WAV_HASH_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";

fn fail(what: &str) -> ! {
    // Structured refusal on stdout so agents parse one stream.
    println!("{{\"suite\":\"music-render\",\"verdict\":\"refused\",\"what\":\"{what}\"}}");
    std::process::exit(1)
}

fn reed_context(samples: usize, block: usize) -> RenderContext {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)
        .unwrap_or_else(|_| fail("gas state refused"));
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.50,
        }],
    };
    let reed = BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
    };
    let voice = ReedBoreVoice::new(
        &duct,
        &air,
        reed,
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        RATE,
        samples,
        None,
    )
    .unwrap_or_else(|e| fail(&format!("reed voice refused: {e:?}")));
    RenderContext::new(vec![RenderVoice::ReedBore(voice)], block)
}

fn string_context(block: usize) -> RenderContext {
    // A plucked steel-ish string: three modes of the 0.65 m / 60 N /
    // 0.6 g/m card (the bake-off fixture's card, reused deliberately so
    // artifacts cross-reference).
    let wave_speed = (60.0f64 / 6.0e-4).sqrt();
    let modes = (1..=3)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: k as f64 * core::f64::consts::PI * wave_speed / 0.65,
            damping_ratio: 1.0e-3 * k as f64,
            pressure_per_modal_velocity: fs_math::c64::C64::new(2.0, 0.0),
        })
        .collect::<Vec<_>>();
    let mut model =
        ModalAcousticTimeModel::try_new(RATE, modes, ModalAcousticTimeBudget::audible_reference())
            .unwrap_or_else(|e| fail(&format!("modal model refused: {e:?}")));
    let pluck: Vec<ModalAcousticState> = (1..=3)
        .map(|k| ModalAcousticState {
            displacement_m_sqrt_kg: 1.0e-3 / f64::from(k),
            velocity_m_sqrt_kg_per_s: 0.0,
        })
        .collect();
    model
        .restore_states(&pluck)
        .unwrap_or_else(|e| fail(&format!("pluck refused: {e:?}")));
    let voice = ModalStringVoice::new(model, vec![0.0; 3])
        .unwrap_or_else(|e| fail(&format!("string voice refused: {e:?}")));
    RenderContext::new(vec![RenderVoice::ModalString(voice)], block)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut positional = Vec::new();
    let mut seconds = 1.0f64;
    let mut block = 512usize;
    let mut full_scale_pa = 200.0f64;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--seconds" => {
                seconds = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--seconds needs a positive number"));
            }
            "--block" => {
                block = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--block needs a positive integer"));
            }
            "--full-scale-pa" => {
                full_scale_pa = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--full-scale-pa needs a positive number"));
            }
            other => positional.push(other.to_string()),
        }
    }
    let [fixture, out_path] = positional.as_slice() else {
        fail(
            "usage: music_render <reed|string> <out.wav> [--seconds S] [--block N] [--full-scale-pa P]",
        );
    };
    if !(seconds > 0.0 && seconds <= 600.0) {
        fail("--seconds must be in (0, 600]");
    }
    if block == 0 || block > 1 << 16 {
        fail("--block must be in 1..=65536");
    }
    let out = std::path::Path::new(out_path);
    let sidecar = out.with_extension("provenance.json");
    if out.exists() || sidecar.exists() {
        fail("output or sidecar already exists; this lane refuses to overwrite evidence");
    }

    let samples = (seconds * f64::from(RATE)).round() as usize;
    let mut context = match fixture.as_str() {
        "reed" => reed_context(samples, block),
        "string" => string_context(block),
        _ => fail("fixture must be `reed` or `string`"),
    };

    // Render block by block through the same API the budget lane measures.
    let mut pressure = vec![0.0f64; samples];
    let mut cursor = 0;
    while cursor < samples {
        let len = block.min(samples - cursor);
        context
            .block(&mut pressure[cursor..cursor + len])
            .unwrap_or_else(|e| fail(&format!("render refused: {e}")));
        cursor += len;
    }

    let (wav, clipped) = encode_pcm16_wav(&pressure, RATE, full_scale_pa)
        .unwrap_or_else(|e| fail(&format!("encode refused: {e}")));
    let wav_hash = hash_domain(WAV_HASH_DOMAIN, &wav);
    let hash_hex: String = wav_hash.0.iter().map(|b| format!("{b:02x}")).collect();
    let peak = pressure.iter().fold(0.0f64, |m, p| m.max(p.abs()));
    let rms = (pressure.iter().map(|p| p * p).sum::<f64>() / pressure.len() as f64).sqrt();

    // Deterministic provenance sidecar: everything a replayer needs. No
    // wall-clock, no commit stamp (the git history of committed artifacts
    // carries those); the WAV content hash is the replay check.
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"{fixture}\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{full_scale_pa:e},\"clipped_samples\":{clipped},\
         \"peak_pa\":{peak:e},\"rms_pa\":{rms:e},\"wav_blake3\":\"{hash_hex}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\"}}"
    );
    std::fs::write(out, &wav).unwrap_or_else(|e| fail(&format!("wav write failed: {e}")));
    std::fs::write(&sidecar, format!("{provenance}\n"))
        .unwrap_or_else(|e| fail(&format!("sidecar write failed: {e}")));

    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"{fixture}\",\
         \"wav\":\"{}\",\"samples\":{samples},\"clipped\":{clipped},\"peak_pa\":{peak:.3},\
         \"rms_pa\":{rms:.3},\"wav_blake3\":\"{hash_hex}\"}}",
        out.display()
    );
}
