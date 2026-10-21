//! G2 loop-closure regression of the opt-in Nolle/Price jet transfer on the real island.
//!
//! The fixture tube (0.3 m, r 6 mm, far-end reflection -0.7, Bernoulli jet, no viscothermal loss) has, from the
//! closed form of Price's J with the module's constants (phase-zero crossing of
//! `Lambda = i w M_d Phi / (i w L + R + Z_in)`), a linear margin of about 14 near 1.03 kHz for flue 10 mm x 1 mm
//! with cut-up 8 mm at 800 Pa, and no crossing at all for flue 4 mm x 0.5 mm with cut-up 5 mm at 100 Pa. A loop that
//! is wired must sustain the first and ring down the second. These are fixture geometries chosen for a large
//! margin, not instrument claims.
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::nolle_jet::{ChannelClass, NolleSpec};
use fs_couple::flue_loop::{FlueControl, FlueGeometry, FlueVoice};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

fn last_block_rms(width: f64, height: f64, cut_up: f64, pressure: f64, hosted: bool) -> f64 {
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
    let sounding = last_block_rms(0.010, 0.001, 0.008, 800.0, true);
    let silent = last_block_rms(0.004, 0.0005, 0.005, 100.0, true);
    println!("hosted transfer: margin-14 fixture rms {sounding:.3e} Pa, no-crossing fixture rms {silent:.3e} Pa");
    assert!(sounding > 1.0, "predicted unstable (margin about 14) but rms {sounding} Pa");
    assert!(silent < 1e-3 * sounding, "predicted stable but rms {silent} Pa against {sounding} Pa");
}

#[test]
fn g2_default_island_is_not_the_hosted_transfer() {
    // Same fixture, transfer not named: the island's own chain runs (opt-in, default untouched). Only
    // records the two levels; the assertion is that naming the transfer changes the jet response.
    let island = last_block_rms(0.010, 0.001, 0.008, 800.0, false);
    let hosted = last_block_rms(0.010, 0.001, 0.008, 800.0, true);
    println!("island chain rms {island:.3e} Pa, hosted transfer rms {hosted:.3e} Pa");
    assert!((island - hosted).abs() > 1e-6 * hosted.max(island), "{island} vs {hosted}");
}
