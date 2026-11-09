use super::*;
use crate::brass_loop::brass_loop_tests::{air, authored_lip, brass_combo, config};
use crate::brass_loop::{BrassControl, BrassVoice};
use crate::mm_line::MmLoad;
use crate::turbulence_noise::{ConstrictionNoiseCard, JetNoiseCard};
use crate::upstream_airway::AirwaySection;
use fs_duct::{LossModel, Termination};
use fs_math::c64::C64;

/// Kaburagi, Yamada, Fukui and Minamiya 2011 (JASA 130(1):536-545), Fig. 3, the LOW panel (a trumpeter at Eb4, 311 Hz), as
/// digitized (+-0.032 cm2): the sectional area [cm2] of 40 equal sections from the lips to the glottis, each of
/// `LOW_LENGTH_CM / 40`. The data file this repeats is a scratch table; the repeat pins the numbers a lip card carries.
const LOW_AREAS_CM2: [f64; 40] = [
    0.142, 0.679, 4.216, 5.289, 4.784, 3.205, 3.016, 2.763, 2.511, 2.574, 1.658, 1.089, 0.584,
    0.458, 0.363, 0.679, 0.584, 0.489, 0.426, 1.279, 0.963, 0.995, 1.153, 1.626, 2.068, 3.458,
    3.047, 2.858, 2.732, 2.984, 4.626, 4.342, 3.995, 2.037, 1.247, 0.995, 0.774, 0.932, 0.995,
    1.216,
];
/// The printed vocal-tract length of the LOW panel [cm].
const LOW_LENGTH_CM: f64 = 17.4;
/// The printed glottis: one uniform tube of 0.14 cm2 and 0.12 cm.
const GLOTTIS_AREA_CM2: f64 = 0.14;
const GLOTTIS_LENGTH_CM: f64 = 0.12;

fn low_chain() -> AirwayChain {
    let section_m = LOW_LENGTH_CM / 40.0 * 1.0e-2;
    let mut tubes: Vec<AirwaySection> = LOW_AREAS_CM2
        .iter()
        .map(|&area_cm2| {
            AirwaySection::from_area(section_m, area_cm2 * 1.0e-4, 1, None).expect("a section")
        })
        .collect();
    tubes.push(
        AirwaySection::from_area(
            GLOTTIS_LENGTH_CM * 1.0e-2,
            GLOTTIS_AREA_CM2 * 1.0e-4,
            1,
            None,
        )
        .expect("the glottis"),
    );
    AirwayChain::new(tubes, LossModel::Bessel).expect("a chain")
}

/// The probe voice of the brass loop's own tests (authored lip, 1e-4 m2 face, tension 1.6), with the tract or without it.
fn voice(tract: Option<(&AirwayChain, TractPressure)>) -> Result<BrassVoice, String> {
    let gas = air(293.15);
    let load = MmLoad::Analytic(Termination::FlangedOpen);
    let mut lip = authored_lip();
    lip.face_area_m2 = 1.0e-4;
    let voice = BrassVoice::new(
        &[brass_combo(0.0)],
        &["open"],
        &gas,
        &load,
        &config(),
        lip,
        1.5e-6,
    )
    .map_err(|e| e.to_string())?;
    let mut voice = match tract {
        Some((chain, pressure)) => {
            voice.with_player_tract(chain, pressure, &gas, &LineLadder::default())?
        }
        None => voice,
    };
    voice
        .apply(BrassControl::SetLipTension(1.6))
        .map_err(|e| e.to_string())?;
    Ok(voice)
}

/// `samples` samples of mouthpiece pressure at a steady `blow` [Pa], one sample per block.
fn render(voice: &mut BrassVoice, blow: f64, samples: usize) -> Vec<f64> {
    voice
        .apply(BrassControl::SetBlowingPressure(blow))
        .expect("blow");
    let mut out = Vec::with_capacity(samples);
    let mut one = [0.0f64];
    for _ in 0..samples {
        voice.step_block(&mut one).expect("block");
        out.push(one[0]);
    }
    out
}

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// G0: a voice without the tract renders the same bits twice, and the opt-in is refused by name where it cannot apply: a
/// voice that has advanced, one that has a tract, one with the jet noise, a bad `tau`.
#[test]
fn no_tract_is_the_original_voice_and_the_opt_in_is_refused_by_name() {
    let plain = render(&mut voice(None).expect("voice"), 3000.0, 2400);
    assert!(
        plain.iter().any(|p| p.abs() > 1.0),
        "the probe cup must carry pressure for the comparison to mean anything"
    );
    assert!(
        same_bits(
            &plain,
            &render(&mut voice(None).expect("voice"), 3000.0, 2400)
        ),
        "the same voice renders the same bits"
    );
    let chain = low_chain();
    let gas = air(293.15);
    let ladder = LineLadder::default();
    let mut advanced = voice(None).expect("voice");
    render(&mut advanced, 3000.0, 1);
    let error = advanced
        .with_player_tract(&chain, TractPressure::Lung, &gas, &ladder)
        .err()
        .expect("an advanced voice refuses");
    assert!(error.contains("has not advanced"), "{error}");
    let held = voice(Some((&chain, TractPressure::Lung))).expect("a lung tract");
    assert!(held.player_tract().is_some());
    let error = held
        .with_player_tract(&chain, TractPressure::Lung, &gas, &ladder)
        .err()
        .expect("a second tract refuses");
    assert!(error.contains("has no tract yet"), "{error}");
    let card = JetNoiseCard {
        seed: 7,
        gain: 1.3e-3,
        strouhal_lo: 0.045,
        strouhal_hi: 0.115,
        reynolds_onset: 1700.0,
    };
    let noisy = voice(None)
        .expect("voice")
        .with_jet_noise(card)
        .expect("a usable card");
    let error = noisy
        .with_player_tract(&chain, TractPressure::Lung, &gas, &ladder)
        .err()
        .expect("the jet noise refuses");
    assert!(error.contains("not combined"), "{error}");
    for tau_s in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = voice(None)
            .expect("voice")
            .with_player_tract(&chain, TractPressure::Mouth { tau_s }, &gas, &ladder)
            .err()
            .expect("a bad tau refuses");
        assert!(
            error.contains("positive finite time constant"),
            "tau {tau_s}: {error}"
        );
    }
}

/// G3: the mouth-pressure tract of the digitized LOW trumpeter tract, `tau` 20 ms. The load is admitted, `R_dc` is the
/// airway's Poiseuille sum, the mean mouth pressure over the last 0.1 s of a 0.3 s hold is the card's pressure (the DC drop
/// comes out through the one-pole), and the sound is not the ideal source's.
#[test]
fn a_mouth_pressure_tract_is_admitted_holds_the_card_pressure_and_changes_the_sound() {
    let chain = low_chain();
    let gas = air(293.15);
    let plain = render(&mut voice(None).expect("voice"), 3000.0, 7200);
    let mut with = voice(Some((&chain, TractPressure::Mouth { tau_s: 0.02 })))
        .expect("the 20 ms one-pole admits the LOW tract");
    let (margin, hz) = with
        .player_tract()
        .expect("tract")
        .passivity_margin()
        .expect("a mouth tract reads its margin");
    assert!(margin >= 0.0 && hz > 0.0, "margin {margin} at {hz} Hz");
    let r_dc = with.player_tract().expect("tract").dc_resistance();
    assert!(
        r_dc > 0.0 && (r_dc - chain.dc_resistance(&gas)).abs() <= 1.0e-9 * r_dc,
        "R_dc {r_dc}"
    );
    println!(
        "LOW tract, tau 20 ms: R_dc {r_dc:.6e} Pa s/m^3; smallest Re Z_in - R_dc Re H {margin:.6e} Pa s/m^3 at {hz:.3} Hz"
    );
    with.apply(BrassControl::SetBlowingPressure(3000.0))
        .expect("blow");
    let (mut out, mut one, mut sum, mut count) = (Vec::new(), [0.0f64], 0.0, 0u32);
    for i in 0..7200 {
        with.step_block(&mut one).expect("block");
        out.push(one[0]);
        if i >= 4800 {
            sum += with.player_tract().expect("tract").mouth_pressure_pa();
            count += 1;
        }
    }
    assert!(out.iter().all(|p| p.is_finite()));
    let mean = sum / f64::from(count);
    println!("mean mouth pressure over the last 0.1 s: {mean:.3} Pa against the card's 3000");
    assert!(
        (mean - 3000.0).abs() <= 30.0,
        "the mean mouth pressure {mean} is not the card's 3000 Pa within 1 percent"
    );
    assert!(
        !same_bits(&plain, &out),
        "the tract must reach the mouthpiece pressure"
    );
}

/// Planted (MB 18707): the flow into the tract is the aperture flow, the flow the cup junction took, `u = (p+ - p-)/Zc0 +
/// u_cup`. The check reads both from the voice's own outputs; a flow 1 percent off fails the same comparison.
#[test]
fn the_flow_into_the_tract_is_the_aperture_flow() {
    let chain = low_chain();
    let mut with =
        voice(Some((&chain, TractPressure::Mouth { tau_s: 0.02 }))).expect("a mouth tract");
    with.apply(BrassControl::SetBlowingPressure(3000.0))
        .expect("blow");
    let zc0 = with.port_impedance();
    let (mut out, mut wave) = ([0.0f64], [0.0f64]);
    let (mut worst, mut planted, mut counted) = (0.0f64, 0.0f64, 0u32);
    for _ in 0..2400 {
        with.step_block_outgoing(&mut out, &mut wave)
            .expect("block");
        let (p_plus, p_now) = (wave[0], out[0]);
        let p_minus = p_now - p_plus;
        let aperture = (p_plus - p_minus) / zc0 + with.cup.u_prev;
        let tract = with.player_tract().expect("tract").flow_m3_s();
        if aperture.abs() > 1.0e-9 {
            counted += 1;
            worst = worst.max((tract - aperture).abs() / aperture.abs());
            planted = planted.max((1.01 * tract - aperture).abs() / aperture.abs());
        }
    }
    assert!(
        counted > 1000,
        "the lip must pass flow for the comparison to mean anything ({counted} samples)"
    );
    assert!(
        worst < 1.0e-9,
        "the tract's flow is not the aperture flow: worst relative difference {worst}"
    );
    assert!(
        planted > 5.0e-3,
        "the planted 1 percent error must be caught by the same comparison ({planted})"
    );
}

/// A negative passivity margin is refused by name with its value and frequency; the same read admits a load that takes
/// nothing out (`H = 0`) and the real one-pole of the LOW tract.
#[test]
fn a_negative_passivity_margin_is_refused_by_name() {
    let chain = low_chain();
    let gas = air(293.15);
    let error = admitted_margin(&chain, &gas, 24_000, 4096, |_| C64::new(1.0e6, 0.0))
        .expect_err("a response of 1e6 takes out far more than the tract holds")
        .to_string();
    assert!(
        error.contains("player tract refused") && error.contains("would be active"),
        "{error}"
    );
    let (margin, _) = admitted_margin(&chain, &gas, 24_000, 4096, |_| C64::ZERO)
        .expect("nothing taken out is passive");
    assert!(margin >= 0.0, "{margin}");
    let pole = 1.0 - fs_math::det::exp(-(1.0 / 24_000.0) / 0.02);
    let dt = 1.0 / 24_000.0;
    admitted_margin(&chain, &gas, 24_000, 4096, |omega| {
        one_pole_response(pole, omega, dt)
    })
    .expect("the 20 ms one-pole is admitted");
    let h0 = one_pole_response(pole, 0.0, dt);
    assert!(
        (h0.re - 1.0).abs() < 1.0e-12 && h0.im.abs() < 1.0e-12,
        "H(0) = 1"
    );
}

/// The constriction noise rides on the tract as a series source (it reads only the held pressure): the voice with both repeats bit
/// for bit under one seed, differs under another, and is not the tract alone's.
#[test]
fn the_constriction_noise_rides_on_the_tract_seeded() {
    let chain = low_chain();
    let gas = air(293.15);
    let build = |seed: u64| -> BrassVoice {
        let mut lip = authored_lip();
        lip.face_area_m2 = 1.0e-4;
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut voice = BrassVoice::new(
            &[brass_combo(0.0)],
            &["open"],
            &gas,
            &load,
            &config(),
            lip,
            1.5e-6,
        )
        .expect("voice")
        .with_constriction_noise(ConstrictionNoiseCard {
            seed,
            duct_radius_m: 8.4e-3,
        })
        .expect("a usable card")
        .with_player_tract(
            &chain,
            TractPressure::Mouth { tau_s: 0.02 },
            &gas,
            &LineLadder::default(),
        )
        .expect("a tract beside the constriction noise");
        voice
            .apply(BrassControl::SetLipTension(1.6))
            .expect("tension");
        voice
    };
    let seven = render(&mut build(7), 3000.0, 2400);
    assert!(seven.iter().all(|p| p.is_finite()));
    assert!(
        same_bits(&seven, &render(&mut build(7), 3000.0, 2400)),
        "a declared seed repeats bit for bit"
    );
    assert!(
        !same_bits(&seven, &render(&mut build(8), 3000.0, 2400)),
        "another seed gives another noise"
    );
    let alone = render(
        &mut voice(Some((&chain, TractPressure::Mouth { tau_s: 0.02 }))).expect("voice"),
        3000.0,
        2400,
    );
    assert!(
        !same_bits(&seven, &alone),
        "the noise must reach the pressure through the tract"
    );
}
