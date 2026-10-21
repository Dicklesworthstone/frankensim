//! G2 loop-closure regression of the opt-in Nolle/Price jet transfer on the real island.
//!
//! The fixture tube (0.3 m, r 6 mm, far-end reflection -0.7, no viscothermal loss) has, from the closed form of
//! Price's J with the module's constants (phase-zero crossing of `Lambda = i w M_d Phi / (i w L + R + Z_in)`, scanned
//! to 6 kHz): a linear margin of about 14 near 1.03 kHz for flue 10 mm x 1 mm with cut-up 8 mm at 800 Pa with the
//! Bernoulli jet, and about 20 near 1.68 kHz for the same flue at 1200 Pa behind a 24.55 mm viscous slit windway
//! (Poiseuille centre speed 1.5 x bulk). A loop that is wired must sustain both. NEGATIVE CONTROL: for flue 10 mm x 2 mm
//! with cut-up 4 mm at 20 Pa the same closed form gives max |Lambda| = 0.055 over the WHOLE band 40 Hz to 24 kHz (Nyquist
//! of the 48 kHz fixture; 8 Hz step), so the loop is stable by small gain and must stay silent; a host that always
//! oscillates fails it. (The same flue as the sounding case, 10 x 1 mm with 8 mm cut-up, has no pressure below its
//! threshold with max |Lambda| < 1, so the control keeps the width and changes height, cut-up and pressure.) These are
//! fixture geometries chosen for a large margin, not instrument claims. (An earlier "no crossing" control was wrong: its scan
//! stopped at 3 kHz and that flue is unstable above it.)
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::nolle_jet::{ChannelClass, NolleSpec};
use fs_couple::flue_loop::{FlueControl, FlueGeometry, FlueVoice};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

fn last_block_rms(width: f64, height: f64, cut_up: f64, pressure: f64, windway_m: Option<f64>, hosted: bool) -> f64 {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let geometry = FlueGeometry {
        flue_width_m: width,
        flue_height_m: height,
        cut_up_m: cut_up,
        labium_offset_m: 0.0,
        jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * 0.006 * 0.006,
        mouth_loss_rel: 0.02,
    };
    let source = FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
            sections: vec![TubeSection { nodes: [0, 1], length_m: 0.3, radius_m: 0.006, max_length_error_m: 0.005 }],
            sound_speed_m_s: gas.sound_speed,
            max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![gas],
        viscothermal: vec![],
        minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3,
        max_jet_history_samples: 4096,
    };
    let card = mint_tonal_interim_card().unwrap();
    let mut voice = FlueVoice::from_network(Some(&card), 90.0, 1e-6, geometry, source, 48_000, 3, 0, &CancelGate::new_clock_free())
        .unwrap();
    if let Some(length) = windway_m {
        voice = voice.with_viscous_windway(length).unwrap();
    }
    if hosted {
        voice = voice.with_nolle_jet(NolleSpec::price_2015(ChannelClass::Long, 0.0).unwrap()).unwrap();
    }
    voice.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
    let mut block = vec![0.0; 2400];
    let mut last = 0.0;
    for _ in 0..16 {
        voice.step_block(&mut block).unwrap();
        last = (block.iter().map(|p| p * p).sum::<f64>() / block.len() as f64).sqrt();
    }
    last
}

#[test]
fn g2_hosted_jet_transfer_sustains_the_loop_where_the_linear_margin_exceeds_one() {
    let rms = last_block_rms(0.010, 0.001, 0.008, 800.0, None, true);
    println!("Bernoulli jet, linear margin about 14: rms {rms:.3e} Pa");
    assert!(rms > 1.0, "predicted unstable (margin about 14) but rms {rms} Pa");
}

#[test]
fn g2_hosted_jet_transfer_sustains_the_loop_behind_a_viscous_windway() {
    let rms = last_block_rms(0.010, 0.001, 0.008, 1200.0, Some(0.02455), true);
    println!("viscous windway (Poiseuille U0 = 1.5 u), linear margin about 20: rms {rms:.3e} Pa");
    assert!(rms > 1.0, "predicted unstable (margin about 20) but rms {rms} Pa");
}

#[test]
fn g2_hosted_jet_transfer_is_silent_where_the_loop_gain_stays_below_one_to_nyquist() {
    let sounding = last_block_rms(0.010, 0.001, 0.008, 800.0, None, true);
    let silent = last_block_rms(0.010, 0.002, 0.004, 20.0, None, true);
    println!("negative control (max |Lambda| 0.055 to 24 kHz): rms {silent:.3e} Pa against sounding {sounding:.3e} Pa");
    assert!(sounding > 1.0, "positive case did not sound: {sounding} Pa");
    assert!(silent < 1e-3 * sounding, "predicted stable everywhere but rms {silent} Pa against {sounding} Pa");
}

#[test]
fn g2_default_island_is_not_the_hosted_transfer() {
    // Same fixture, transfer not named: the island's own chain runs (opt-in, default untouched).
    let island = last_block_rms(0.010, 0.001, 0.008, 800.0, None, false);
    let hosted = last_block_rms(0.010, 0.001, 0.008, 800.0, None, true);
    println!("island chain rms {island:.3e} Pa, hosted transfer rms {hosted:.3e} Pa");
    assert!((island - hosted).abs() > 1e-6 * hosted.max(island), "{island} vs {hosted}");
}
