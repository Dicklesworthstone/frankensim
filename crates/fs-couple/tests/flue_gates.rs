//! Flue playing-loop fixtures (music bead `frankensim-music-v8-root-3ez8g.10.3`):
//! the jet-drive island composed with the wind line. Every fixture
//! input is a chart gesture, a medium state, or a card-backed number;
//! no frequency is assigned anywhere (grep this file: no `hz` input
//! reaches the loop). The regime maps printed here are the bead's
//! logged evidence; the assertions are the falsifiers.

use fs_aeroac::jetcard::{JetCard, mint_tonal_interim_card};
use fs_couple::flue_loop::forced_jet::PRODUCTION_SPEC;
use fs_couple::flue_loop::{FlueControl, FlueError, FlueGeometry, FlueRank, FlueVoice};
use fs_duct::{Duct, Fingering, FingeringTable, HoleState, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};

/// Line sample rate. The duct's radiation law admits `ka < 1` up to the
/// line's Nyquist; a flute-class 6 mm bore therefore needs a Nyquist
/// under 9 kHz. 12 kHz keeps the first five modes of the fixture pipe
/// (567 Hz fundamental) inside the line, which is what the ladder
/// needs; audio bandwidth is not this fixture's subject.
const RATE: u32 = 12_000;
const BLOCK: usize = 2048;
/// Bore radius [m] (`ka = 0.66` at the 6 kHz Nyquist).
const BORE_R: f64 = 0.006;
/// Default cut-up [m].
const CUT_UP: f64 = 0.004;
/// Authored spatial gain of the jet deflection over the cut-up (no
/// receptivity gain on the tonal card; labeled authored in the island's
/// provenance). Estimate: sinuous-mode spatial growth of a Bickley jet
/// is about `mu b = 0.45` near its most amplified Strouhal, so over a
/// cut-up of `W / b = 8` the deflection grows by `exp(3.6) = 37`; the
/// two-section band limit costs about half of that at the fundamental,
/// so 90 is the estimate restored to the band's passband.
const GAIN: f64 = 90.0;
/// Onset seed relative to jet speed (lab practice; not a noise claim).
/// MEASURED: at 1e-3 the seed alone, integrated through the receptivity,
/// flaps the jet by its own half-width and the loop never locks (the
/// labium's slope averages to nothing under the noise); 1e-6 leaves the
/// seed three decades under the jet width and the loop free to grow.
const SEED_REL: f64 = 1e-6;

/// Lip-jet cut-up [m] of the overblowing and voicing fixtures. The island's
/// band lags the loop by `16 theta / U` on top of the convective transit
/// `W / (c_v U)`; the early tap that removes that lag needs the transit to
/// be the longer of the two, i.e. `W > 16 theta c_v` (5.6 mm for the lip
/// card below at a 1 mm flue). The interim fixture's 4 mm cut-up is under
/// that limit for a lip jet.
const LIP_CUT_UP: f64 = 0.007;
/// Flue height [m] of `geometry` (kept in step with the literal there).
const FLUE_H: f64 = 0.001;
/// Convection velocity over jet centreline speed (Auvray et al. 2012;
/// Terrien et al. 2014 cite 0.4; de la Cuadra 2008 cites 0.3 to 0.5).
const LIP_CONVECTION: f64 = 0.4;
/// theta / b of the island's band for a lip jet: the band peaks at
/// `0.55 omega_c` with `omega_c = U / (4 theta)`; placing that peak at the
/// Strouhal number of maximum jet amplification, `Str_h = f h / U = 0.025`
/// (de la Cuadra 2008, Sec. 7), gives `theta = 0.875 h`, so `theta / b =
/// 1.75` for `b = h / 2`. The interim card's 0.4 is the thickness of a 2-D
/// lattice slot jet at Re 144 to 264.
const LIP_THETA_OVER_B: f64 = 1.75;
/// Published first-register window `Str_W = f1 W / U` (de la Cuadra 2008:
/// the first register holds for 0.1 to 0.25; the octave jump sits near 0.1).
const FIRST_REGISTER_WINDOW: (f64, f64) = (0.10, 0.25);

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

/// The lip-jet card: a published-parameter estimate (authority X-Est), the
/// same one `flue_jet_phase.rs` pins its first register with. Its Reynolds
/// band is de la Cuadra's 1000 to 2500, so the ladder rows between roughly
/// 135 and 845 Pa are inside it; the card gate (fg_007) is unchanged and
/// stays on the interim card.
fn lip_card() -> JetCard {
    // c_p / U = St (h / delta) / (stage + 1/4) with stage 1 and h / delta 10,
    // so St = c_p / U / 8.
    let st = LIP_CONVECTION / 8.0;
    let json = format!(
        "{{\"schema\":\"fs-aeroac.jet-card/v1\",\"claim_kind\":\"edge-tone-tonal\",\"claim\":{{\"stage\":1,\
         \"locked_strouhal\":{st},\"ladder_band_lo\":0.7,\"ladder_band_hi\":1.4,\"multi_stable_strouhal\":[{st}],\
         \"hysteresis_recorded\":true,\"saturated_force_rms\":null}},\"authority\":\"X-Est\",\
         \"profile\":{{\"u_centerline\":0.08,\"slot_half\":1,\"momentum_thickness\":{LIP_THETA_OVER_B}}},\
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

/// The pair of whole-millimetre cut-ups that both lie inside the published
/// first-register window at jet speed `u` (`Str_W = f1 W / U` in
/// `FIRST_REGISTER_WINDOW`, `f1` the DESIGNED fundamental) and above the
/// early-tap limit `16 theta c_v`. The ends of the usable window rounded
/// inward: derived from the window, never from which pair happens to pass.
fn window_cut_ups(f1_hz: f64, u: f64) -> (f64, f64) {
    let theta = LIP_THETA_OVER_B * 0.5 * FLUE_H;
    let tap_limit = 16.0 * theta * LIP_CONVECTION;
    let lo = (FIRST_REGISTER_WINDOW.0 * u / f1_hz).max(tap_limit);
    let hi = FIRST_REGISTER_WINDOW.1 * u / f1_hz;
    let mm = 1e-3;
    ((lo / mm).ceil() * mm, (hi / mm).floor() * mm)
}

/// A recorder-class open cylinder with one tone hole near the foot,
/// closed in the only fingering (the chart law needs at least one
/// hole; a closed hole leaves the bore an open-open pipe of the full
/// length up to its small chimney compliance).
fn open_pipe(length_m: f64, radius_m: f64) -> FingeringTable {
    FingeringTable::try_new(
        Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: radius_m,
                    length: 0.9 * length_m,
                },
                Segment::ToneHole {
                    hole_radius: 0.0015,
                    chimney_height: 0.0015,
                    bore_radius: radius_m,
                    state: HoleState::Closed,
                },
                Segment::Cylinder {
                    radius: radius_m,
                    length: 0.1 * length_m,
                },
            ],
        },
        vec![Fingering {
            label: "closed".to_string(),
            holes: vec![HoleState::Closed],
        }],
    )
    .expect("table admits")
}

/// Flute-class mouth for the 6 mm bore: a 10 mm wide, 1 mm high flue
/// with a 4 mm cut-up (mouth area 0.35 of the bore area).
fn geometry(cut_up_m: f64, labium_offset_m: f64, jet_angle_rad: f64) -> FlueGeometry {
    FlueGeometry {
        flue_width_m: 0.010,
        flue_height_m: FLUE_H,
        cut_up_m,
        labium_offset_m,
        jet_angle_rad,
        bore_area_m2: core::f64::consts::PI * BORE_R * BORE_R,
        mouth_loss_rel: 0.02,
    }
}

fn voice(
    card: Option<&JetCard>,
    geom: FlueGeometry,
    seed: u64,
    index: u32,
) -> Result<FlueVoice, FlueError> {
    FlueVoice::new(
        card,
        GAIN,
        SEED_REL,
        geom,
        &open_pipe(0.30, BORE_R),
        &air20(),
        Termination::UnflangedOpen,
        RATE,
        seed,
        index,
    )
}

/// Run `blocks` blocks at one blowing pressure and return the last
/// block's lock estimate and RMS (after the transient).
fn settle(v: &mut FlueVoice, blow_pa: f64, blocks: usize) -> (f64, f64) {
    v.apply(FlueControl::SetBlowingPressure(blow_pa))
        .expect("pressure admits");
    let mut out = vec![0.0; BLOCK];
    for _ in 0..blocks {
        v.step_block(&mut out).expect("block steps");
    }
    let d = v.diagnostics().last().expect("diag");
    (d.lock_hz, d.p_rms_pa)
}

/// Open-open pipe fundamental estimate with end corrections, for the
/// regime map's reference column only (never fed to the loop).
fn pipe_mode_hz(gas: &GasState, length_m: f64, radius_m: f64, n: f64) -> f64 {
    n * gas.sound_speed / (2.0 * (length_m + 2.0 * 0.61 * radius_m))
}

#[test]
fn fg_001_no_card_refuses_and_a_refusal_boundary_card_cannot_drive() {
    let geom = geometry(CUT_UP, 0.0, 0.0);
    match voice(None, geom, 1, 0) {
        Err(FlueError::NoCard) => {}
        other => panic!("no card must refuse: {:?}", other.err()),
    }
    // The real refusal-boundary card the 3-D campaign minted (retained
    // by fs-aeroac's jc_014): it says where broadband does NOT exist and
    // must not drive a tone.
    let retained = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../fs-aeroac/tests/receipts/jet-card-refusal-boundary.json");
    let boundary = JetCard::from_json(&std::fs::read_to_string(&retained).expect("retained card"))
        .expect("retained card parses");
    match voice(Some(&boundary), geom, 1, 0) {
        Err(FlueError::CardClassCannotDrive { kind }) => {
            assert_eq!(kind, "broadband-refusal-boundary");
        }
        other => panic!("refusal-boundary card must not drive: {:?}", other.err()),
    }
    let tonal = mint_tonal_interim_card().expect("tonal interim card mints");
    let v = voice(Some(&tonal), geom, 1, 0).expect("tonal card drives");
    let island = v.island();
    assert_eq!(island.claim_kind, "edge-tone-tonal");
    assert!(
        island.convection_ratio > 0.2 && island.convection_ratio < 0.4,
        "{island:?}"
    );
    assert!(island.provenance.contains("card-backed") && island.provenance.contains("authored"));
    println!(
        "{{\"fixture\":\"fg-001\",\"island\":{:?}}}",
        island.provenance
    );
}

#[test]
fn fg_002_energy_ledger_closes_to_roundoff_every_block() {
    let card = mint_tonal_interim_card().expect("card");
    let mut v = voice(Some(&card), geometry(CUT_UP, 0.0, 0.0), 2, 0).expect("voice");
    let _ = settle(&mut v, 300.0, 40);
    let mut worst = 0.0f64;
    for d in v.diagnostics() {
        let scale = d.source_work_j.abs()
            + d.pipe_work_j.abs()
            + d.radiation_loss_j
            + d.stored_delta_j.abs()
            + 1e-30;
        worst = worst.max(d.ledger_defect_j.abs() / scale);
        assert!(d.numerical_dissipation_j >= 0.0);
        assert!(d.radiation_loss_j >= 0.0);
    }
    println!("{{\"fixture\":\"fg-002\",\"worst_relative_ledger_defect\":{worst:e}}}");
    assert!(worst < 1e-9, "ledger must close to roundoff: {worst:e}");
}

/// THE OVERBLOWING LADDER: sweeping blowing pressure (jet speed) walks
/// the lock from the pipe's first mode to its second. Nothing in the
/// sweep names a frequency; the reference modes are computed for the
/// log only.
///
/// On the LIP-JET card (see `lip_card`): the interim slot-jet card is
/// outside its own Reynolds band (144 to 264) at every row of this sweep.
/// The sweep keeps its 1.5-ratio ladder and is continued past the old top
/// (360 Pa) to the island's own incompressibility bound (jet Mach 0.2, about
/// 2.8 kPa), where the old top-of-sweep bar (`>= 1.7 f1`) is asked of the
/// physics instead of being lowered: a lip jet's first register holds to a
/// much faster jet than the interim card's did.
#[test]
fn fg_003_overblowing_ladder_walks_fundamental_to_octave() {
    let card = lip_card();
    let gas = air20();
    let f1 = pipe_mode_hz(&gas, 0.30, BORE_R, 1.0);
    let f2 = pipe_mode_hz(&gas, 0.30, BORE_R, 2.0);
    let mut v = voice(Some(&card), geometry(LIP_CUT_UP, 0.0, 0.0), 3, 0).expect("voice");
    let pressures = [
        10.0, 15.0, 22.0, 33.0, 50.0, 75.0, 110.0, 160.0, 240.0, 360.0, 540.0, 810.0, 1200.0,
        1800.0, 2700.0,
    ];
    let mach_top = (2.0 * pressures[pressures.len() - 1] / gas.density).sqrt() / gas.sound_speed;
    assert!(
        mach_top < 0.2,
        "the sweep stays incompressible: Mach {mach_top}"
    );
    let mut map = Vec::new();
    println!(
        "REGIME MAP (blow_pa, U, Re, inside_card, transit_ms, lock_hz, lock/f1, p_rms, eta_rms/b, duty):"
    );
    for &p in &pressures {
        let (lock, rms) = settle(&mut v, p, 60);
        let d = *v.diagnostics().last().expect("diag");
        println!(
            "{p:.0}\t{:.2}\t{:.0}\t{}\t{:.3}\t{lock:.1}\t{:.3}\t{rms:.3e}\t{:.3}\t{:.3}",
            d.jet_speed_m_s,
            d.reynolds,
            d.inside_card_validity,
            d.transit_s * 1e3,
            lock / f1,
            d.eta_rms_over_b,
            d.duty_pipe_side
        );
        map.push((p, lock, rms));
    }
    println!("reference open-open modes f1={f1:.1} f2={f2:.1} (log only)");
    let oscillating: Vec<&(f64, f64, f64)> = map.iter().filter(|(_, _, rms)| *rms > 1.0).collect();
    assert!(
        oscillating.len() >= 3,
        "the loop must self-sustain at several pressures: {map:?}"
    );
    let low = oscillating.first().expect("first");
    let high = oscillating.last().expect("last");
    let ratio_low = low.1 / f1;
    let ratio_high = high.1 / f1;
    assert!(
        (0.8..1.25).contains(&ratio_low),
        "lowest oscillating pressure locks near the fundamental: {ratio_low}"
    );
    // The ladder: some higher pressure locks the octave, and the top of
    // the sweep sits at or above it (a strong jet overblows further).
    let octave = oscillating
        .iter()
        .find(|(p, lock, _)| *p > low.0 && (1.7..2.4).contains(&(lock / f1)))
        .expect("a higher pressure locks near the octave");
    println!(
        "{{\"fixture\":\"fg-003\",\"fundamental_pa\":{},\"octave_pa\":{},\"top_pa\":{},\"top_ratio\":{ratio_high}}}",
        low.0, octave.0, high.0
    );
    assert!(
        ratio_high >= 1.7,
        "the top of the sweep stays at or above the octave: {ratio_high}"
    );

    // SELF-SUSTAINED, not seed-filtered: a limit cycle's amplitude does
    // not depend on the onset seed level, while a resonant pipe merely
    // filtering the seed scales with it. Re-run the lowest oscillating
    // pressure with the seed ten times smaller; the RMS must agree
    // within 20 percent and the lock must not move.
    let mut quiet = FlueVoice::new(
        Some(&card),
        GAIN,
        SEED_REL / 10.0,
        geometry(LIP_CUT_UP, 0.0, 0.0),
        &open_pipe(0.30, BORE_R),
        &gas,
        Termination::UnflangedOpen,
        RATE,
        3,
        0,
    )
    .expect("voice");
    let (lock_quiet, rms_quiet) = settle(&mut quiet, low.0, 60);
    println!(
        "{{\"fixture\":\"fg-003\",\"seed_check_pa\":{},\"rms_seeded\":{},\"rms_quiet_seed\":{rms_quiet},\"lock_seeded\":{},\"lock_quiet_seed\":{lock_quiet}}}",
        low.0, low.2, low.1
    );
    assert!(
        (rms_quiet / low.2 - 1.0).abs() < 0.2,
        "a limit cycle does not scale with the seed: {} vs {rms_quiet}",
        low.2
    );
    assert!(
        (lock_quiet / low.1 - 1.0).abs() < 0.05,
        "the lock does not depend on the seed: {} vs {lock_quiet}",
        low.1
    );
}

/// The gate's voice with the opt-in forced jet on the production row (`PRODUCTION_SPEC`),
/// integrated at the pressure it is held at.
fn forced_voice(
    card: &JetCard,
    geom: FlueGeometry,
    seed_rel: f64,
    seed: u64,
    index: u32,
    held_pa: f64,
) -> Result<FlueVoice, FlueError> {
    FlueVoice::new(
        Some(card),
        GAIN,
        seed_rel,
        geom,
        &open_pipe(0.30, BORE_R),
        &air20(),
        Termination::UnflangedOpen,
        RATE,
        seed,
        index,
    )?
    .with_forced_jet(PRODUCTION_SPEC, held_pa)
}

/// fg_003 once more under the opt-in forced jet (ruling 7518), with the same assertions. Each row
/// is a voice built cold at its own held pressure, since the response is integrated at the
/// Reynolds number of the held pressure and claims no sweep between rows. It fails today and
/// stays ignored until it passes unchanged; the default switches only then.
#[test]
#[ignore = "frankensim-irp-v5-instrument-realism-pexfp.135.1: forced jet lands high: the onset row, 75 Pa, locks 2.00 f1; 110-360 Pa 0.97-1.00, 540-1200 Pa 1.95-1.99, 1800-2700 Pa 2.97-2.99"]
fn fg_003_forced_jet_overblowing_ladder_walks_fundamental_to_octave() {
    let card = lip_card();
    let gas = air20();
    let f1 = pipe_mode_hz(&gas, 0.30, BORE_R, 1.0);
    let f2 = pipe_mode_hz(&gas, 0.30, BORE_R, 2.0);
    let pressures = [
        10.0, 15.0, 22.0, 33.0, 50.0, 75.0, 110.0, 160.0, 240.0, 360.0, 540.0, 810.0, 1200.0,
        1800.0, 2700.0,
    ];
    let mach_top = (2.0 * pressures[pressures.len() - 1] / gas.density).sqrt() / gas.sound_speed;
    assert!(
        mach_top < 0.2,
        "the sweep stays incompressible: Mach {mach_top}"
    );
    let mut map = Vec::new();
    println!(
        "FORCED-JET REGIME MAP (blow_pa, U, Re, inside_card, transit_ms, lock_hz, lock/f1, p_rms, eta_rms/b, duty):"
    );
    // Each row integrates the production response at its own Reynolds number (minutes each), so
    // the rows run concurrently; the map is printed in pressure order.
    let rows: Vec<_> = std::thread::scope(|s| {
        let rows: Vec<_> = pressures
            .iter()
            .map(|&p| {
                s.spawn(move || {
                    let geom = geometry(LIP_CUT_UP, 0.0, 0.0);
                    let mut v = forced_voice(&lip_card(), geom, SEED_REL, 3, 0, p).expect("voice");
                    let (lock, rms) = settle(&mut v, p, 60);
                    (p, lock, rms, *v.diagnostics().last().expect("diag"))
                })
            })
            .collect();
        rows.into_iter()
            .map(|row| row.join().expect("row"))
            .collect()
    });
    for (p, lock, rms, d) in rows {
        println!(
            "{p:.0}\t{:.2}\t{:.0}\t{}\t{:.3}\t{lock:.1}\t{:.3}\t{rms:.3e}\t{:.3}\t{:.3}",
            d.jet_speed_m_s,
            d.reynolds,
            d.inside_card_validity,
            d.transit_s * 1e3,
            lock / f1,
            d.eta_rms_over_b,
            d.duty_pipe_side
        );
        map.push((p, lock, rms));
    }
    println!("reference open-open modes f1={f1:.1} f2={f2:.1} (log only)");
    let oscillating: Vec<&(f64, f64, f64)> = map.iter().filter(|(_, _, rms)| *rms > 1.0).collect();
    assert!(
        oscillating.len() >= 3,
        "the loop must self-sustain at several pressures: {map:?}"
    );
    let low = oscillating.first().expect("first");
    let high = oscillating.last().expect("last");
    let ratio_low = low.1 / f1;
    let ratio_high = high.1 / f1;
    assert!(
        (0.8..1.25).contains(&ratio_low),
        "lowest oscillating pressure locks near the fundamental: {ratio_low}"
    );
    // The ladder: some higher pressure locks the octave, and the top of
    // the sweep sits at or above it (a strong jet overblows further).
    let octave = oscillating
        .iter()
        .find(|(p, lock, _)| *p > low.0 && (1.7..2.4).contains(&(lock / f1)))
        .expect("a higher pressure locks near the octave");
    println!(
        "{{\"fixture\":\"fg-003-forced\",\"fundamental_pa\":{},\"octave_pa\":{},\"top_pa\":{},\"top_ratio\":{ratio_high}}}",
        low.0, octave.0, high.0
    );
    assert!(
        ratio_high >= 1.7,
        "the top of the sweep stays at or above the octave: {ratio_high}"
    );

    // SELF-SUSTAINED, not seed-filtered: a limit cycle's amplitude does
    // not depend on the onset seed level, while a resonant pipe merely
    // filtering the seed scales with it. Re-run the lowest oscillating
    // pressure with the seed ten times smaller; the RMS must agree
    // within 20 percent and the lock must not move.
    let mut quiet = forced_voice(
        &card,
        geometry(LIP_CUT_UP, 0.0, 0.0),
        SEED_REL / 10.0,
        3,
        0,
        low.0,
    )
    .expect("voice");
    let (lock_quiet, rms_quiet) = settle(&mut quiet, low.0, 60);
    println!(
        "{{\"fixture\":\"fg-003-forced\",\"seed_check_pa\":{},\"rms_seeded\":{},\"rms_quiet_seed\":{rms_quiet},\"lock_seeded\":{},\"lock_quiet_seed\":{lock_quiet}}}",
        low.0, low.2, low.1
    );
    assert!(
        (rms_quiet / low.2 - 1.0).abs() < 0.2,
        "a limit cycle does not scale with the seed: {} vs {rms_quiet}",
        low.2
    );
    assert!(
        (lock_quiet / low.1 - 1.0).abs() < 0.05,
        "the lock does not depend on the seed: {} vs {lock_quiet}",
        low.1
    );
}

/// VOICING: a larger cut-up lengthens the jet transit, so onset needs
/// more pressure; a labium offset towards the pipe side shifts the duty
/// of the cut jet. Both logged and asserted in direction only.
///
/// On the LIP-JET card (see `lip_card`). The cut-up pair is derived, not
/// typed: both cut-ups must lie inside the published first-register window
/// (`Str_W = f1 W / U` in 0.10 to 0.25, `f1` the designed fundamental)
/// because "a longer cut-up locks the same or a lower mode" is a law of one
/// edge-tone stage only (the stage ladder is a sawtooth: frequency falls with
/// W inside a stage and jumps up at a stage boundary), and above the island's
/// early-tap limit `16 theta c_v`. At `LAW_BLOW_PA` that window is 5.6 to 8.9
/// mm, so the pair is 6 mm and 8 mm; 240 Pa is the first row of the fg_003
/// ladder whose window holds two whole-millimetre cut-ups at least 2 mm apart
/// (at 160 Pa it holds 6 and 7 mm only).
#[test]
fn fg_004_voicing_gestures_move_onset_and_duty_in_the_expected_direction() {
    const LAW_BLOW_PA: f64 = 240.0;
    let card = lip_card();
    let gas = air20();
    let f1 = pipe_mode_hz(&gas, 0.30, BORE_R, 1.0);
    let u_law = (2.0 * LAW_BLOW_PA / gas.density).sqrt();
    let (cut_short, cut_long) = window_cut_ups(f1, u_law);
    assert!(
        cut_long - cut_short >= 2e-3 - 1e-12,
        "the window holds two distinct cut-ups: {cut_short} {cut_long}"
    );
    for w in [cut_short, cut_long] {
        let st_w = f1 * w / u_law;
        assert!(
            (FIRST_REGISTER_WINDOW.0..=FIRST_REGISTER_WINDOW.1).contains(&st_w),
            "cut-up {w} m sits inside the first-register window: Str_W {st_w}"
        );
    }
    let onset = |cut_up: f64| -> Option<f64> {
        let mut v = voice(Some(&card), geometry(cut_up, 0.0, 0.0), 4, 0).expect("voice");
        for &p in &[
            20.0, 35.0, 50.0, 70.0, 100.0, 140.0, 200.0, 280.0, 400.0, 560.0, 800.0,
        ] {
            let (_, rms) = settle(&mut v, p, 40);
            if rms > 1.0 {
                return Some(p);
            }
        }
        None
    };
    let short = onset(cut_short);
    let long = onset(cut_long);
    println!(
        "{{\"fixture\":\"fg-004\",\"cut_up_short_m\":{cut_short},\"cut_up_long_m\":{cut_long},\"onset_pa_short\":{short:?},\"onset_pa_long\":{long:?}}}"
    );
    assert!(
        short.is_some() && long.is_some(),
        "both cut-ups must reach onset in the sweep"
    );
    // The effect this island models for the cut-up is the transit delay:
    // at equal blowing pressure a longer cut-up lengthens the delay and
    // the loop locks the same or a lower pipe mode. (Onset pressure vs
    // cut-up also depends on jet spreading beyond the card's profile,
    // which is not modeled; it is logged above, not asserted.)
    let lock_at = |cut_up: f64, blow_pa: f64| -> f64 {
        let mut v = voice(Some(&card), geometry(cut_up, 0.0, 0.0), 4, 2).expect("voice");
        let (lock, rms) = settle(&mut v, blow_pa, 60);
        assert!(rms > 1.0, "cut-up {cut_up} must oscillate at {blow_pa} Pa");
        lock
    };
    let lock_short = lock_at(cut_short, LAW_BLOW_PA);
    let lock_long = lock_at(cut_long, LAW_BLOW_PA);
    println!(
        "{{\"fixture\":\"fg-004\",\"blow_pa\":{LAW_BLOW_PA},\"lock_hz_short\":{lock_short},\"lock_hz_long\":{lock_long},\"lock_over_f1_short\":{},\"lock_over_f1_long\":{}}}",
        lock_short / f1,
        lock_long / f1
    );
    assert!(
        lock_long <= lock_short * 1.05,
        "a longer cut-up (longer transit) locks the same or a lower mode: {lock_short} vs {lock_long}"
    );

    let duty = |offset: f64| -> f64 {
        let mut v = voice(Some(&card), geometry(LIP_CUT_UP, offset, 0.0), 4, 1).expect("voice");
        let _ = settle(&mut v, 110.0, 40);
        v.diagnostics().last().expect("diag").duty_pipe_side
    };
    let centered = duty(0.0);
    let inward = duty(0.0003);
    println!(
        "{{\"fixture\":\"fg-004\",\"duty_centered\":{centered},\"duty_labium_inward\":{inward}}}"
    );
    assert!(
        inward < centered,
        "moving the labium into the jet cuts less of it onto the pipe side"
    );
}

/// fg_004 once more under the opt-in forced jet (ruling 7518; mail 7859 makes it the fourth gate of
/// the switch), with the same assertions, seeds and voice indices. The onset sweep is one voice held
/// at the sweep's top pressure, as a card's `jet_response` row holds the gesture's highest blowing
/// pressure; the lock and duty voices are held at the pressure they run at. The six voices are
/// independent, so they run concurrently.
#[test]
fn fg_004_forced_jet_voicing_gestures_move_onset_and_duty_in_the_expected_direction() {
    const LAW_BLOW_PA: f64 = 240.0;
    const ONSET_SWEEP_PA: [f64; 11] = [
        20.0, 35.0, 50.0, 70.0, 100.0, 140.0, 200.0, 280.0, 400.0, 560.0, 800.0,
    ];
    let gas = air20();
    let f1 = pipe_mode_hz(&gas, 0.30, BORE_R, 1.0);
    let u_law = (2.0 * LAW_BLOW_PA / gas.density).sqrt();
    let (cut_short, cut_long) = window_cut_ups(f1, u_law);
    assert!(
        cut_long - cut_short >= 2e-3 - 1e-12,
        "the window holds two distinct cut-ups: {cut_short} {cut_long}"
    );
    for w in [cut_short, cut_long] {
        let st_w = f1 * w / u_law;
        assert!(
            (FIRST_REGISTER_WINDOW.0..=FIRST_REGISTER_WINDOW.1).contains(&st_w),
            "cut-up {w} m sits inside the first-register window: Str_W {st_w}"
        );
    }
    let onset = |cut_up: f64| -> Option<f64> {
        let top = ONSET_SWEEP_PA[ONSET_SWEEP_PA.len() - 1];
        let mut v = forced_voice(&lip_card(), geometry(cut_up, 0.0, 0.0), SEED_REL, 4, 0, top)
            .expect("voice");
        for &p in &ONSET_SWEEP_PA {
            let (_, rms) = settle(&mut v, p, 40);
            if rms > 1.0 {
                return Some(p);
            }
        }
        None
    };
    let lock_at = |cut_up: f64, blow_pa: f64| -> (f64, f64) {
        let mut v = forced_voice(
            &lip_card(),
            geometry(cut_up, 0.0, 0.0),
            SEED_REL,
            4,
            2,
            blow_pa,
        )
        .expect("voice");
        settle(&mut v, blow_pa, 60)
    };
    let duty = |offset: f64| -> f64 {
        let mut v = forced_voice(
            &lip_card(),
            geometry(LIP_CUT_UP, offset, 0.0),
            SEED_REL,
            4,
            1,
            110.0,
        )
        .expect("voice");
        let _ = settle(&mut v, 110.0, 40);
        v.diagnostics().last().expect("diag").duty_pipe_side
    };
    let (short, long, (lock_short, rms_short), (lock_long, rms_long), centered, inward) =
        std::thread::scope(|s| {
            let short = s.spawn(move || onset(cut_short));
            let long = s.spawn(move || onset(cut_long));
            let at_short = s.spawn(move || lock_at(cut_short, LAW_BLOW_PA));
            let at_long = s.spawn(move || lock_at(cut_long, LAW_BLOW_PA));
            let centered = s.spawn(move || duty(0.0));
            let inward = s.spawn(move || duty(0.0003));
            (
                short.join().expect("onset short"),
                long.join().expect("onset long"),
                at_short.join().expect("lock short"),
                at_long.join().expect("lock long"),
                centered.join().expect("duty centered"),
                inward.join().expect("duty inward"),
            )
        });
    println!(
        "{{\"fixture\":\"fg-004-forced\",\"cut_up_short_m\":{cut_short},\"cut_up_long_m\":{cut_long},\"onset_pa_short\":{short:?},\"onset_pa_long\":{long:?}}}"
    );
    println!(
        "{{\"fixture\":\"fg-004-forced\",\"blow_pa\":{LAW_BLOW_PA},\"lock_hz_short\":{lock_short},\"lock_hz_long\":{lock_long},\"lock_over_f1_short\":{},\"lock_over_f1_long\":{},\"rms_short\":{rms_short},\"rms_long\":{rms_long}}}",
        lock_short / f1,
        lock_long / f1
    );
    println!(
        "{{\"fixture\":\"fg-004-forced\",\"duty_centered\":{centered},\"duty_labium_inward\":{inward}}}"
    );
    assert!(
        short.is_some() && long.is_some(),
        "both cut-ups must reach onset in the sweep"
    );
    for (cut_up, rms) in [(cut_short, rms_short), (cut_long, rms_long)] {
        assert!(
            rms > 1.0,
            "cut-up {cut_up} must oscillate at {LAW_BLOW_PA} Pa"
        );
    }
    assert!(
        lock_long <= lock_short * 1.05,
        "a longer cut-up (longer transit) locks the same or a lower mode: {lock_short} vs {lock_long}"
    );
    assert!(
        inward < centered,
        "moving the labium into the jet cuts less of it onto the pipe side"
    );
}

/// ORGAN RANK: N voices from one chart with per-voice logging; every
/// voice's output is bitwise the solo run of the same voice, and the
/// stepping order does not matter.
#[test]
fn fg_005_rank_voices_are_independent_and_order_free() {
    let card = mint_tonal_interim_card().expect("card");
    let gas = air20();
    let table = open_pipe(0.30, BORE_R);
    let geoms = [
        geometry(CUT_UP, 0.0, 0.0),
        geometry(0.0045, 0.0001, 0.0),
        geometry(0.0035, -0.0001, 0.02),
        geometry(0.005, 0.0, -0.02),
    ];
    let run_rank = |order: &[usize]| -> Vec<Vec<f64>> {
        let mut rank = FlueRank::new(
            Some(&card),
            GAIN,
            SEED_REL,
            &geoms,
            &table,
            &gas,
            Termination::UnflangedOpen,
            RATE,
            5,
        )
        .expect("rank");
        for (i, v) in rank.voices_mut().iter_mut().enumerate() {
            v.apply(FlueControl::SetBlowingPressure(200.0 + 50.0 * i as f64))
                .expect("pressure");
        }
        let mut outs = vec![vec![0.0; BLOCK]; geoms.len()];
        let mut all = vec![Vec::new(); geoms.len()];
        for _ in 0..12 {
            rank.step_block(order, &mut outs).expect("rank steps");
            for (k, o) in outs.iter().enumerate() {
                all[k].extend_from_slice(o);
            }
        }
        for (i, v) in rank.voices().iter().enumerate() {
            let d = v.diagnostics().last().expect("diag");
            println!(
                "{{\"fixture\":\"fg-005\",\"voice\":{i},\"seed_tile\":{},\"lock_hz\":{:.1},\"p_rms\":{:.3e}}}",
                v.seed_key().tile,
                d.lock_hz,
                d.p_rms_pa
            );
        }
        all
    };
    let forward = run_rank(&[0, 1, 2, 3]);
    let reversed = run_rank(&[3, 2, 1, 0]);
    assert_eq!(
        forward, reversed,
        "stepping order is a scheduling choice, not a physics input"
    );
    for (i, g) in geoms.iter().enumerate() {
        let mut solo = voice(Some(&card), *g, 5, i as u32).expect("solo");
        solo.apply(FlueControl::SetBlowingPressure(200.0 + 50.0 * i as f64))
            .expect("pressure");
        let mut out = vec![0.0; BLOCK];
        let mut all = Vec::new();
        for _ in 0..12 {
            solo.step_block(&mut out).expect("solo steps");
            all.extend_from_slice(&out);
        }
        assert_eq!(
            all, forward[i],
            "voice {i} in the rank is bitwise its solo self"
        );
    }
}

/// DETERMINISM: the same logical voice identity replays bitwise; a
/// different voice identity under the same seed does not collide.
#[test]
fn fg_006_replay_is_bitwise_and_keyed_by_voice_identity() {
    let card = mint_tonal_interim_card().expect("card");
    let render = |index: u32| -> Vec<f64> {
        let mut v = voice(Some(&card), geometry(CUT_UP, 0.0, 0.0), 9, index).expect("voice");
        v.apply(FlueControl::SetBlowingPressure(250.0))
            .expect("pressure");
        let mut out = vec![0.0; BLOCK];
        let mut all = Vec::new();
        for _ in 0..8 {
            v.step_block(&mut out).expect("steps");
            all.extend_from_slice(&out);
        }
        all
    };
    let a = render(0);
    let b = render(0);
    let c = render(1);
    assert_eq!(a, b, "same voice identity replays bitwise");
    assert_ne!(
        a, c,
        "a different voice identity draws a different onset seed"
    );
}

/// CLAIM BOUNDARY: blocks outside the card's Reynolds band are refused
/// by the claim gate even though the loop ran them. MEASURED: the tonal
/// interim card admits jet Re 144..264 (the 2-D lattice rig at
/// h/delta = 10); a 1.2 mm flue in air at 300 Pa runs at Re ~ 1700, so
/// the whole playing range of this fixture is OUTSIDE the card — the
/// loop still runs (it is physics), but no claim may be minted from it
/// until a card covering that band exists. That is the boundary this
/// bead's dependents inherit, not a defect to paper over.
#[test]
fn fg_007_claims_outside_the_card_validity_region_are_refused() {
    let card = mint_tonal_interim_card().expect("card");
    let mut v = voice(Some(&card), geometry(CUT_UP, 0.0, 0.0), 11, 0).expect("voice");
    let _ = settle(&mut v, 300.0, 2);
    match v.claim_check() {
        Err(FlueError::OutsideCardValidity {
            blocks,
            reynolds_seen,
            reynolds_card,
        }) => {
            assert_eq!(blocks, 2);
            assert!(
                reynolds_seen.0 > reynolds_card.1,
                "{reynolds_seen:?} vs {reynolds_card:?}"
            );
            println!(
                "{{\"fixture\":\"fg-007\",\"reynolds_seen\":[{},{}],\"reynolds_card\":[{},{}]}}",
                reynolds_seen.0, reynolds_seen.1, reynolds_card.0, reynolds_card.1
            );
        }
        other => panic!("claims outside the card must refuse: {other:?}"),
    }
    // Inside the band the gate admits: the blowing pressure that puts the
    // jet at the middle of the card's Reynolds band, derived from the
    // card, the medium, and the flue height (no literal).
    let gas = air20();
    let geom = geometry(CUT_UP, 0.0, 0.0);
    let (re_lo, re_hi) = v.island().reynolds_band;
    let nu = gas.dynamic_viscosity / gas.density;
    let u_mid = 0.5 * (re_lo + re_hi) * nu / geom.flue_height_m;
    let p_mid = 0.5 * gas.density * u_mid * u_mid;
    let mut slow = voice(Some(&card), geom, 11, 1).expect("voice");
    let _ = settle(&mut slow, p_mid, 1);
    let d = slow.diagnostics()[0];
    println!(
        "{{\"fixture\":\"fg-007\",\"slow_reynolds\":{},\"inside\":{}}}",
        d.reynolds, d.inside_card_validity
    );
    assert!(d.inside_card_validity, "{d:?}");
    slow.claim_check()
        .expect("in-band blocks pass the claim gate");
}

/// GATE SUMMARY (bead 3ez8g.10.4 logging clause): the committed flue gate
/// summary names every flue and jet row of the registry, the refused
/// broadband row reads refused in BOTH files, and the summary's digest line
/// is the registry's actual blake2b-128 (a stale summary is a silently
/// wrong gate table).
#[test]
fn fg_008_flue_gate_summary_enumerates_every_registry_row() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/flue-gate-summary.tsv"))
        .expect("committed flue gate summary");
    assert!(summary.starts_with("# frankensim-flue-gate-summary-v1"));
    // The header's blake2b-128 is the maintainer's regeneration receipt
    // (computed outside the tree); here it must at least be present and
    // well-formed. The row enumeration below is the gate.
    let header = summary.lines().nth(1).expect("digest line");
    let stated = header
        .split("blake2b-128 ")
        .nth(1)
        .and_then(|rest| rest.split(')').next())
        .expect("stated digest");
    assert_eq!(stated.len(), 32, "{header}");
    let mut keys = Vec::new();
    for filling in ["\"filling\": \"flue\"", "\"filling\": \"jet\""] {
        let mut cursor = 0usize;
        while let Some(hit) = registry[cursor..].find(filling) {
            let at = cursor + hit;
            let tail = &registry[at..];
            let image = field(tail, "\"image\": \"");
            let qoi = field(tail, "\"qoi\": \"");
            let gate = field(tail, "\"gate\": \"");
            keys.push((image, qoi, gate));
            cursor = at + filling.len();
        }
    }
    assert!(
        keys.len() >= 10,
        "the flue menu plus the jet labs: {keys:?}"
    );
    for (image, qoi, gate) in &keys {
        let row = format!("\t{image}\t{qoi}\t{gate}\t");
        assert!(summary.contains(&row), "summary must carry {row:?}");
    }
    assert!(
        summary.contains("\tjet-island-char-line\tbroadband-flute-noise\trefused\t"),
        "the broadband refusal reads refused in the summary"
    );
    assert!(
        registry.contains("\"qoi\": \"broadband-flute-noise\"")
            && keys
                .iter()
                .any(|(_, q, g)| q == "broadband-flute-noise" && g == "refused"),
        "the broadband refusal reads refused in the registry"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"flue-gate-summary\",\"verdict\":\"pass\",\"rows\":{}}}",
        keys.len()
    );
}

fn field(tail: &str, tag: &str) -> String {
    let at = tail.find(tag).expect("field") + tag.len();
    let end = tail[at..].find('"').expect("end");
    tail[at..at + end].to_string()
}
