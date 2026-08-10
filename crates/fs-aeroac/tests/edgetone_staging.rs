//! Edge-tone Strouhal STAGING against published data (bead 9ok02's
//! acceptance item). Literature provenance, fetched 2026-08-09 (two
//! sources, saved to the session scratchpad as
//! edgetone/vaik-paal-part{1,2}.pdf):
//!
//! - Vaik, Varga & Paal, "Frequency and Phase Characteristics of the
//!   Edge Tone Part I", Period. Polytech. Mech. Eng. 58(1) 2014,
//!   Table 1: Brown's (1937) stage-I coefficients for
//!   `St = (c1 - c2/Re)((delta/h)^k - c3)`, c1 = 0.4659,
//!   c2 = 12.06, c3 = 0.007, k = 1 (validity h/delta in 3.1..60,
//!   Re in 75..1300, claimed max 6% deviation).
//! - Part II, Table 1 (same journal issue): the authors' own top-hat
//!   measurements and CFD at h/delta ~ 10,
//!   `St = c/Re + St_inf`: experiment c = -1.150,
//!   St_inf = 0.04522; CFD c = -0.7387, St_inf = 0.04010 (pure
//!   stage I).
//!
//! The fixture runs the jet-labium rig at the canonical staging
//! geometry h/delta = 10, Re = 144 (inside the stage-I band: onset
//! ~75, stage II ~220 per Part I) and compares the measured
//! `St = f delta / u` against all three published predictions. The
//! DEFAULT jetlab fixture (h/delta = 3) is deliberately NOT used
//! here: it sits at the very edge of Brown's validity band and its
//! oscillation is not on the stage ladder (measured St 0.64 —
//! recorded as a negative result in the bead).
//!
//! `#[ignore]`d for runtime (~4 min debug): executed on demand and on
//! the record — the JSON line below was produced by a real run.
//!
//! SEED PROVENANCE (executed, all three states): unseeded, the
//! mirror-symmetric rig amplifies ROUNDOFF at the physically
//! selected frequency (same St 0.03662, but the amplitude is
//! machine noise — a frequency-selection result only); seeded at
//! 0.005 the SATURATED limit cycle (prominence ~1e20, real
//! amplitude) locks at the same St 0.03662, +3.0% vs Brown within
//! the +-6% bin quantization — the claim below; seeded at 0.02 the
//! rig settles on a NEIGHBORING locked state (St 0.0458, +29%) —
//! multi-stability/hysteresis of edge tones is documented in the
//! literature (Brown 1937; Part I's stage-overlap discussion), and
//! the strong-seed state is recorded, not asserted against.
//!
//! OPEN SCOPE (probed, recorded on the bead, not yet passing —
//! refusal beats fabrication): stage II at Re 243 measured
//! St 0.0549, BETWEEN rungs (the literature's documented
//! transition/coexistence band near the stage-II onset); grid
//! convergence at delta = 9 lu is blocked by a slit-lip mode
//! (St_delta ~ 2-3) that outgrows the edge tone — refining by
//! slowing u weakens the edge-tone force ~u^4 (executed), and
//! viscosity-scaled refinement still leaves the fast mode dominant
//! (executed); the suspected mechanism is the profile/wall mismatch
//! at the nozzle-slit shoulder, a recorded follow-up.

use fs_aeroac::jetlab::{
    JetLabiumConfig, RampConfig, RampDirection, run_adiabatic_ramp, run_jet_labium,
};
use fs_fft::{C64, Fft};

/// Brown's (1937) stage-I Strouhal prediction (Vaik/Paal Part I,
/// Table 1 coefficients) at h/delta = 10.
fn st_brown_stage_one(re: f64) -> f64 {
    (0.4659 - 12.06 / re) * (0.1 - 0.007)
}

#[test]
#[ignore = "heavy staging run (~4 min); execute explicitly"]
fn edge_tone_stage_one_strouhal_matches_published() {
    let cfg = JetLabiumConfig {
        nx: 192,
        ny: 64,
        slot_half: 3.0,
        slot_smoothing: 1.2,
        u_jet: 0.08,
        tau: 0.51,         // nu = 1/300 -> Re = u * 2*slot_half / nu = 144
        edge_distance: 60, // h/delta = 60/6 = 10
        plate_length: 50,
        fringe_width: 32,
        fringe_sigma: 0.3,
        steps_settle: 4000,
        steps_record: 16_384,
        seed_amplitude: 0.005,
        nozzle_thickness: 2,
    };
    let run = run_jet_labium(&cfg).expect("run");
    let d = &run.diagnostics;
    assert!(d.mach_max_lattice < 0.25, "Mach {}", d.mach_max_lattice);
    let imbalance = (d.flux_plate_plane - d.flux_fringe_plane).abs() / d.flux_plate_plane.abs();
    assert!(imbalance < 0.02, "flux imbalance {imbalance:.4}");
    // Transverse-force spectrum (Hann).
    let n = run.force_series.len();
    let mean = run.force_series.iter().map(|f| f[1]).sum::<f64>() / n as f64;
    let fft = Fft::new(n);
    let mut buf: Vec<C64> = run
        .force_series
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let w = 0.5 - 0.5 * ((2.0 * core::f64::consts::PI * i as f64) / (n as f64 - 1.0)).cos();
            C64::new((f[1] - mean) * w, 0.0)
        })
        .collect();
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    let power: Vec<f64> = buf[..n / 2].iter().map(|c| c.norm_sq()).collect();
    let (peak_bin, peak_pow) = power
        .iter()
        .enumerate()
        .skip(8)
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("spectrum");
    let mut sorted = power[8..].to_vec();
    sorted.sort_by(f64::total_cmp);
    let prominence = peak_pow / sorted[sorted.len() / 2].max(1e-300);
    assert!(
        prominence > 50.0,
        "no oscillation: prominence {prominence:.1}"
    );
    let freq = peak_bin as f64 / n as f64;
    let delta = 2.0 * cfg.slot_half;
    let st = freq * delta / cfg.u_jet;
    // Published stage-I predictions at Re = 144, h/delta = 10:
    let re = d.reynolds;
    let st_brown = (0.4659 - 12.06 / re) * (delta / cfg.edge_distance as f64 - 0.007);
    let st_vaik_exp = -1.150 / re + 0.045_22;
    let st_vaik_cfd = -0.738_7 / re + 0.040_10;
    // Authored envelope: within 15% of the Brown prediction
    // (measured +3.0% on the recorded run, INSIDE the published
    // sources' own ~8% spread — but note the spectral peak sits at
    // bin ~8 of the 16384-step record, so the frequency quantization
    // alone is +-6%: the honest claim is agreement within
    // quantization, not 3% precision. The two-cell plate vs wedge,
    // the fringe closure, and delta = 6 lu discretization are the
    // rig's other error sources — re-measure if re-dimensioned).
    let dev = (st - st_brown) / st_brown;
    assert!(
        dev.abs() < 0.15,
        "stage-I Strouhal {st:.4} vs Brown {st_brown:.4} (dev {dev:.2}); Vaik exp {st_vaik_exp:.4}, CFD {st_vaik_cfd:.4}"
    );
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"edge-tone-staging\",\"re\":{re:.0},\"h_over_delta\":10,\"st_measured\":{st:.5},\"st_brown\":{st_brown:.5},\"st_vaik_exp\":{st_vaik_exp:.5},\"st_vaik_cfd\":{st_vaik_cfd:.5},\"deviation_vs_brown\":{dev:.3},\"prominence\":{prominence:.0},\"flux_imbalance\":{imbalance:.5},\"verdict\":\"pass\"}}"
    );
}

/// Shared measurement: run a config, return
/// (St, prominence, imbalance, mach, force_rms).
fn measure_st(cfg: &JetLabiumConfig) -> (f64, f64, f64, f64, f64) {
    let run = run_jet_labium(cfg).expect("run");
    let d = &run.diagnostics;
    let imbalance = (d.flux_plate_plane - d.flux_fringe_plane).abs() / d.flux_plate_plane.abs();
    let n = run.force_series.len();
    let mean = run.force_series.iter().map(|f| f[1]).sum::<f64>() / n as f64;
    let rms = (run
        .force_series
        .iter()
        .map(|f| (f[1] - mean) * (f[1] - mean))
        .sum::<f64>()
        / n as f64)
        .sqrt();
    let fft = Fft::new(n);
    let mut buf: Vec<C64> = run
        .force_series
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let w = 0.5 - 0.5 * ((2.0 * core::f64::consts::PI * i as f64) / (n as f64 - 1.0)).cos();
            C64::new((f[1] - mean) * w, 0.0)
        })
        .collect();
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    let power: Vec<f64> = buf[..n / 2].iter().map(|c| c.norm_sq()).collect();
    let (peak_bin, peak_pow) = power
        .iter()
        .enumerate()
        .skip(8)
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("spectrum");
    let mut sorted = power[8..].to_vec();
    sorted.sort_by(f64::total_cmp);
    let prominence = peak_pow / sorted[sorted.len() / 2].max(1e-300);
    let st = (peak_bin as f64 / n as f64) * 2.0 * cfg.slot_half / cfg.u_jet;
    (st, prominence, imbalance, d.mach_max_lattice, rms)
}

/// SLIT-FIX REGRESSION at the finer lattice (delta = 7.5 lu,
/// geometric similarity at 1.25x, viscosity-scaled to Re 144): the
/// fine rig must lock a LADDER-SCALE mode — before the binary-slit
/// fringe fix it locked a slit-lip mode at St_delta ~ 2-3 (and a
/// non-similar config locked St 0.82). Recorded execution:
/// St 0.04578 (1.29x Brown stage-I), prominence ~2e20, imbalance
/// 0.02%, real amplitude. HONEST LIMIT (recorded, not hidden):
/// strict cross-resolution convergence of the SELECTED attractor is
/// NOT demonstrated — the rig is multi-stable (coarse at seed 0.005
/// locks St 0.0366; coarse at seed 0.02 AND fine at seed 0.005 lock
/// the neighboring St ~ 0.0458 state), so the gate is the ladder
/// BAND [0.7, 1.4] x Brown, and attractor-selection convergence
/// remains open scope on the bead.
#[test]
#[ignore = "heavy staging run (~5 min); execute explicitly"]
fn edge_tone_fine_lattice_slitfix_regression() {
    let fine = JetLabiumConfig {
        // GEOMETRIC SIMILARITY at 1.25x (delta = 7.5 lu): every
        // length scales with delta (executed lesson: a fine config
        // with a shorter dimensionless domain and plate locked onto
        // a different attractor, St 0.82 — convergence studies must
        // scale ALL lengths, and the slow stage-I mode needs settle
        // measured in ITS periods, ~2.7 here).
        nx: 240,
        ny: 80,
        slot_half: 3.75,
        slot_smoothing: 1.5,
        u_jet: 0.08,
        tau: 0.512_5,      // nu = 0.08 * 7.5 / 144 -> Re = 144
        edge_distance: 75, // h/delta = 10
        plate_length: 62,
        fringe_width: 40,
        fringe_sigma: 0.3,
        steps_settle: 7000,
        steps_record: 16_384,
        seed_amplitude: 0.005,
        nozzle_thickness: 2,
    };
    let (st_f, prom_f, imb_f, mach_f, rms_f) = measure_st(&fine);
    let st_brown = (0.4659 - 12.06 / 144.0) * (0.1 - 0.007);
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"edge-tone-slitfix-regression\",\"st_fine\":{st_f:.5},\"st_brown\":{st_brown:.5},\"prominence\":{prom_f:.0},\"imbalance\":{imb_f:.5},\"mach\":{mach_f:.3},\"force_rms\":{rms_f:.3e}}}"
    );
    assert!(prom_f > 50.0 && imb_f < 0.02 && mach_f < 0.25);
    assert!(rms_f > 1.0e-6, "force at machine-noise scale: {rms_f:.3e}");
    // Ladder-band gate (recorded execution: St 0.04578 = 1.29x
    // Brown): the slit-lip (St_delta 2-3, i.e. >5x Brown) and
    // free-jet (St 0.82, >20x) failure modes sit far outside.
    let ratio = st_f / st_brown;
    assert!(
        (0.7..1.4).contains(&ratio),
        "fine lattice off the stage-I ladder band: St {st_f:.4} = {ratio:.2}x Brown"
    );
}

/// ADIABATIC HYSTERESIS-FOLLOWING RAMP (the bead's recorded
/// follow-up to the multi-stability finding): lock the stage-I
/// attractor at Re 144, then quasi-statically sweep Reynolds
/// (viscosity at fixed jet speed — Mach, fringe, slit, and seed all
/// unchanged) up to Re 264 and back, measuring a spectral rung after
/// every transition. Fixed-parameter runs cannot demonstrate stage
/// selection on a multi-stable rig; branch-following can.
/// Transition budget: ramp 30000 steps (~15 stage-I periods) +
/// settle 8000 (~4 periods) per rung; record 16384 (~8 periods,
/// +-6% quantization at the stage-I bin).
///
/// EXECUTED FINDINGS (2026-08-10, two ramp rates — the committed
/// gentle protocol AND a 5x-faster variant with rungs at every
/// Re 20, ramp 6000/settle 4000 — agree qualitatively, so the
/// result is robust to ramp rate over that range):
///
/// 1. NO STAGE-II LOCK by Re 264: every rung on both legs stays in
///    the stage-I NEIGHBORHOOD (this run: 0.72-1.29x Brown stage I;
///    fast variant: 0.72-1.55x). The literature's stage-II onset
///    (~Re 220-250 for a wedge) does not reproduce on this
///    two-cell-plate, periodic-fringe, 2D rig by Re 264 — recorded
///    as a NO-CLAIM boundary, not asserted against. Candidate
///    mechanisms (unproven): plate-vs-wedge receptivity, fringe
///    recirculation phase, 2D confinement (ny = 64).
/// 2. The rig WANDERS between neighboring locked states (bins 6-8 =
///    the 0.0366 family vs bins 9-11 = the 0.0458 family) along the
///    sweep — mode competition, consistent with the fixed-run
///    multi-stability finding.
/// 3. HYSTERESIS IS REAL AND REPRODUCIBLE: the down-leg returns to
///    Re 144 on a DIFFERENT attractor than the up-leg started from
///    (this run: bin 9 / St 0.0412 vs start bin 7 / St 0.0320;
///    fast variant: bin 12 vs 7). Both runs bitwise deterministic,
///    so the path dependence is a pinnable same-ISA measurement.
///
/// Pinned gates: per-rung health, the rung-0 stage-I anchor, the
/// ladder-NEIGHBORHOOD band [0.6, 1.7]x Brown on every rung (under
/// branch-following, the slit-lip St_delta 2-3 (>5x) and free-jet
/// St 0.46 (>12x) pathologies never appear at any Reynolds in
/// 144..264, either direction — a regression surface fixed-Re runs
/// cannot provide), and the witnessed start-vs-end attractor
/// difference at Re 144.
#[test]
#[ignore = "heavy ramp protocol (~340k lattice steps); execute explicitly, prefer release"]
fn edge_tone_adiabatic_ramp_hysteresis() {
    let cfg = RampConfig {
        base: JetLabiumConfig {
            nx: 192,
            ny: 64,
            slot_half: 3.0,
            slot_smoothing: 1.2,
            u_jet: 0.08,
            tau: 0.51, // Re0 = 144
            edge_distance: 60,
            plate_length: 50,
            fringe_width: 32,
            fringe_sigma: 0.3,
            steps_settle: 4000,
            steps_record: 16_384, // unused by the ramp
            seed_amplitude: 0.005,
            nozzle_thickness: 2,
        },
        reynolds_end: 264.0,
        rungs: 4, // Re = 144, 184, 224, 264
        steps_ramp: 30_000,
        steps_rung_settle: 8000,
        steps_rung_record: 16_384,
        skip_bins: 6,
    };
    let report = run_adiabatic_ramp(&cfg).expect("ramp");
    assert_eq!(report.rungs.len(), 7);
    for r in &report.rungs {
        let dir = match r.direction {
            RampDirection::Up => "up",
            RampDirection::Down => "down",
        };
        let ratio = r.peak.strouhal / st_brown_stage_one(r.reynolds);
        println!(
            "{{\"suite\":\"fs-aeroac\",\"case\":\"edge-tone-ramp\",\"direction\":\"{dir}\",\"re\":{:.0},\"tau\":{:.6},\"st\":{:.5},\"bin\":{},\"ratio_vs_brown_stage1\":{:.3},\"prominence\":{:.0},\"force_rms\":{:.3e},\"imbalance\":{:.5},\"mach\":{:.3}}}",
            r.reynolds,
            r.tau,
            r.peak.strouhal,
            r.peak.bin,
            ratio,
            r.peak.prominence,
            r.peak.force_rms,
            r.flux_imbalance,
            r.mach_max_lattice
        );
    }
    // Health gates on every rung: in-regime flow, a real (seeded,
    // saturated) oscillation, and a prominent peak.
    for r in &report.rungs {
        assert!(r.mach_max_lattice < 0.25, "Mach {}", r.mach_max_lattice);
        assert!(r.flux_imbalance < 0.02, "imbalance {}", r.flux_imbalance);
        assert!(
            r.peak.prominence > 50.0,
            "no oscillation at Re {}: prominence {:.1}",
            r.reynolds,
            r.peak.prominence
        );
        assert!(
            r.peak.force_rms > 1.0e-6,
            "machine-noise amplitude at Re {}: {:.3e}",
            r.reynolds,
            r.peak.force_rms
        );
    }
    // Rung 0 anchors the protocol on the validated stage-I state
    // (executed: St 0.03204 = 0.90x Brown at Re 144 — one bin below
    // the staging fixture's 0.03662, inside quantization).
    let first = &report.rungs[0];
    let ratio0 = first.peak.strouhal / st_brown_stage_one(first.reynolds);
    assert!(
        (0.7..1.4).contains(&ratio0),
        "ramp did not start on the stage-I ladder: {ratio0:.2}x Brown"
    );
    // Ladder-neighborhood band on EVERY rung, both legs (executed
    // extremes across both ramp rates: 0.717 and 1.546): under
    // adiabatic branch-following the rig never leaves the stage-I
    // neighborhood — no slit-lip (>5x), no free-jet (>12x), and no
    // stage-II (~2.3x) excursion anywhere in Re 144..264.
    for r in &report.rungs {
        let ratio = r.peak.strouhal / st_brown_stage_one(r.reynolds);
        assert!(
            (0.6..1.7).contains(&ratio),
            "off the stage-I neighborhood at Re {} ({:?}): {ratio:.2}x Brown",
            r.reynolds,
            r.direction
        );
    }
    // Hysteresis witness (executed, bitwise-deterministic same-ISA):
    // the down-leg ends at Re 144 on a different attractor (bin 9,
    // St 0.0412) than the up-leg started from (bin 7, St 0.0320).
    let last = &report.rungs[report.rungs.len() - 1];
    assert!(
        (last.reynolds - first.reynolds).abs() < 1e-9,
        "protocol must close the loop at the base Reynolds"
    );
    assert_ne!(
        first.peak.bin, last.peak.bin,
        "hysteresis witness vanished: start and end select the same bin \
         ({}) — if a rig change made the loop reversible, re-measure and \
         re-pin the multi-stability findings",
        first.peak.bin
    );
}
