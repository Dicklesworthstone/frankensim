//! Jet-island loop phase: the first register at flutists' Strouhal numbers (IRP bead pexfp.320).
//!
//! Gauntlet G1/G3. The amplification band of the jet island is a magnitude; its phase is the
//! convective delay alone. A lip-jet card whose band is placed by the published Strouhal number of
//! maximum jet amplification must therefore speak the pipe's first register while the jet length
//! W and speed U keep Str_W = f1 W / U between about 0.1 and 0.25 (de la Cuadra, Fabre,
//! Montgermont and Chafe 2008, Acta Acustica 94(5), Secs. 5.1-5.2: first register at 0.1 to 0.25,
//! octave jumps near 0.1) and must not do so by seed filtering.
//!
//! NOT ASSERTED, open: the published octave jump. A Python mirror of the island (not native) holds
//! the first register down to Str_W 0.05 with only isolated higher-register points, and the point
//! moves with the sample rate, so no threshold is claimed. The strong-jet lock is logged below.
//!
//! The pipe is the plain 0.3 m, 6 mm bore network of examples/flue.performance (end reflection
//! -0.9, no end correction). Nothing here names a frequency: the reference f1 = c / 2L is only used
//! to label which register the loop locked, never fed to the loop. This is a numerical fixture,
//! not a flute: the card below is a published-parameter estimate (X-Est), not a measured lip-jet
//! card, and the authored receptivity gain is the demonstrator value of the existing gates.

use fs_aeroac::jetcard::JetCard;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::{FlueControl, FlueGeometry, FlueVoice};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
const BLOCK: usize = 4096;
const PIPE_M: f64 = 0.30;
const BORE_R: f64 = 0.006;
const FLUE_H: f64 = 0.001;
/// Jet length (cut-up) [m]. Long enough that the band's low-frequency lag, 16 theta / U, stays
/// under the convective transit W / (c_v U) at every pressure used below (ratio about 0.8).
const CUT_UP: f64 = 0.007;
const GAIN: f64 = 90.0;
const SEED_REL: f64 = 1e-6;
/// Convection velocity over jet centreline speed (Auvray et al. 2012; Terrien et al. 2014 cite
/// 0.4; de la Cuadra 2008 cites 0.3 to 0.5).
const CONVECTION: f64 = 0.4;
/// theta / b of the island's band for a lip jet. The island's band peaks at 0.55 omega_c with
/// omega_c = U / (4 theta). Placing that peak at the Strouhal number of maximum jet amplification,
/// Str_h = f h / U = 0.025 (de la Cuadra 2008, Sec. 7), gives theta = 0.55 h / (8 pi 0.025) =
/// 0.875 h, so theta / b = 1.75 for b = h / 2. This replaces the interim card's 0.4, which is
/// the thickness of a 2-D lattice slot jet at Re 144 to 264.
const THETA_OVER_B: f64 = 1.75;

fn lip_card() -> JetCard {
    // c_p / U = St (h / delta) / (stage + 1/4) with stage 1 and h / delta 10, so St = c_p / U / 8.
    let st = CONVECTION / 8.0;
    let json = format!(
        "{{\"schema\":\"fs-aeroac.jet-card/v1\",\"claim_kind\":\"edge-tone-tonal\",\"claim\":{{\"stage\":1,\
         \"locked_strouhal\":{st},\"ladder_band_lo\":0.7,\"ladder_band_hi\":1.4,\"multi_stable_strouhal\":[{st}],\
         \"hysteresis_recorded\":true,\"saturated_force_rms\":null}},\"authority\":\"X-Est\",\
         \"profile\":{{\"u_centerline\":0.08,\"slot_half\":1,\"momentum_thickness\":{THETA_OVER_B}}},\
         \"band_db\":null,\"residual\":null,\"validity\":{{\"reynolds_lo\":1000,\"reynolds_hi\":2500,\
         \"h_over_delta\":10,\"requires_receptivity_edge\":true,\"seeded_amplitude_claims_only\":true}},\
         \"provenance\":{{\"rig_fingerprint\":0,\"receipts\":[\"de la Cuadra et al. 2008 Acta Acustica \
         94(5) Sec. 7: maximum amplification at Str_h 0.025; Re 1000 to 2500 (Eq. 5)\"],\
         \"classification\":\"published-parameter estimate for a lip jet, not a measured jet card\"}},\
         \"scope\":\"{scope}\"}}",
        scope = fs_aeroac::SCOPE_STATEMENT
    );
    JetCard::from_json(&json).expect("lip-jet card parses and validates")
}

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn voice(card: &JetCard, seed_rel: f64, seed: u64) -> FlueVoice {
    let gas = air();
    let geometry = FlueGeometry {
        flue_width_m: 0.010,
        flue_height_m: FLUE_H,
        cut_up_m: CUT_UP,
        labium_offset_m: 0.0,
        jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * BORE_R * BORE_R,
        mouth_loss_rel: 0.02,
    };
    let source = FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![
                NetworkNode::Inlet,
                NetworkNode::Termination { reflection: -0.9 },
            ],
            sections: vec![TubeSection {
                nodes: [0, 1],
                length_m: PIPE_M,
                radius_m: BORE_R,
                max_length_error_m: 0.005,
            }],
            sound_speed_m_s: gas.sound_speed,
            max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![gas],
        viscothermal: vec![],
        minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3,
        max_jet_history_samples: 4096,
    };
    FlueVoice::from_network(
        Some(card),
        GAIN,
        seed_rel,
        geometry,
        source,
        RATE,
        seed,
        0,
        &CancelGate::new_clock_free(),
    )
    .expect("voice")
}

/// Goertzel magnitude at `hz` over `x` (Hann windowed).
fn tone(x: &[f64], hz: f64) -> f64 {
    let n = x.len() as f64;
    let w = core::f64::consts::TAU * hz / f64::from(RATE);
    let (mut s1, mut s2) = (0.0, 0.0);
    for (i, v) in x.iter().enumerate() {
        let hann = 0.5 - 0.5 * (core::f64::consts::TAU * i as f64 / (n - 1.0)).cos();
        let s0 = v * hann + 2.0 * w.cos() * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    (s1 * s1 + s2 * s2 - 2.0 * w.cos() * s1 * s2).sqrt()
}

/// Settle at one blowing pressure; return (dominant tone [Hz], RMS [Pa]) of the last two blocks.
fn settle(v: &mut FlueVoice, blow_pa: f64) -> (f64, f64) {
    v.apply(FlueControl::SetBlowingPressure(blow_pa))
        .expect("pressure admits");
    let mut tail = Vec::new();
    for block in 0..6 {
        let mut out = vec![0.0; BLOCK];
        v.step_block(&mut out).expect("block steps");
        if block >= 4 {
            tail.extend_from_slice(&out);
        }
    }
    let rms = (tail.iter().map(|p| p * p).sum::<f64>() / tail.len() as f64).sqrt();
    let (mut best, mut best_hz) = (0.0, 0.0);
    let mut hz = 150.0;
    while hz <= 6000.0 {
        let m = tone(&tail, hz);
        if m > best {
            best = m;
            best_hz = hz;
        }
        hz += 4.0;
    }
    (best_hz, rms)
}

fn f1_ref_hz() -> f64 {
    air().sound_speed / (2.0 * PIPE_M)
}

/// The published first register at the flutist's Str_W, from a limit cycle that does not depend on
/// the onset seed. The overblowing threshold is logged, not asserted (see the module note).
#[test]
fn fj_001_lip_jet_speaks_the_first_register_and_ignores_the_seed() {
    let card = lip_card();
    let f1 = f1_ref_hz();
    let u_of = |p: f64| (2.0 * p / air().density).sqrt();
    // Str_W = f1 W / U for the pressures asserted: 0.19, 0.15, 0.12.
    let first_register = [250.0, 400.0, 600.0];
    let mut seen = Vec::new();
    for &p in &first_register {
        let strw = f1 * CUT_UP / u_of(p);
        assert!(
            (0.1..0.25).contains(&strw),
            "fixture Str_W {strw} must sit in the published first-register range"
        );
        let mut v = voice(&card, SEED_REL, 3);
        let (hz, rms) = settle(&mut v, p);
        println!(
            "{{\"fixture\":\"fj-001\",\"blow_pa\":{p},\"str_w\":{strw:.3},\"lock_hz\":{hz},\"ratio\":{:.3},\"rms_pa\":{rms}}}",
            hz / f1
        );
        seen.push((p, hz, rms));
        assert!(rms > 1.0, "the loop must self-sustain at {p} Pa: rms {rms}");
        assert!(
            (0.85..1.15).contains(&(hz / f1)),
            "{p} Pa (Str_W {strw:.3}) must speak the first register: {hz} Hz vs f1 {f1}"
        );
    }
    // A strong jet (Str_W under 0.09, U above about 40 m/s): logged, not asserted, that it moves up.
    let (hz_hi, rms_hi) = settle(&mut voice(&card, SEED_REL, 3), 1500.0);
    println!(
        "{{\"fixture\":\"fj-001\",\"blow_pa\":1500,\"str_w\":{:.3},\"lock_hz\":{hz_hi},\"ratio\":{:.3},\"rms_pa\":{rms_hi}}}",
        f1 * CUT_UP / u_of(1500.0),
        hz_hi / f1
    );
    assert!(
        rms_hi > 1.0,
        "the strong jet must still sound: rms {rms_hi}"
    );

    // SELF-SUSTAINED, not seed-filtered: a limit cycle does not scale with the onset seed and the
    // lock does not depend on the seed stream. Ten times smaller seed, another stream.
    let (p, hz, rms) = seen[1];
    let (hz_q, rms_q) = settle(&mut voice(&card, SEED_REL / 10.0, 11), p);
    assert!(
        (rms_q / rms - 1.0).abs() < 0.2,
        "a limit cycle does not scale with the seed: {rms} vs {rms_q}"
    );
    assert!(
        (hz_q / hz - 1.0).abs() < 0.03,
        "the lock does not depend on the seed: {hz} vs {hz_q}"
    );
}
