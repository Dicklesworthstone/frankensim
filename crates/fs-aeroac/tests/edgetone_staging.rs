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

use fs_aeroac::jetlab::{JetLabiumConfig, run_jet_labium};
use fs_fft::{C64, Fft};

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
