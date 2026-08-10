//! Noise-table conformance: a real velocity sweep produces a fitted
//! table (real-amplitude entries, measured power exponent reported
//! as data, in-regime diagnostics per entry), the JSON export
//! carries the scope statement (the marketing-claim mutation guard),
//! the demo synth consumes the table and reproduces its spectral
//! shape (round-trip), refusals are typed, and everything is
//! deterministic.

use fs_aeroac::jetlab::JetLabiumConfig;
use fs_aeroac::noisetable::{N_BANDS, fit_noise_table};
use fs_aeroac::{AeroacError, SCOPE_STATEMENT};
use fs_lbm::core2::CollisionModel2;

fn geometry() -> JetLabiumConfig {
    JetLabiumConfig {
        nx: 128,
        ny: 64,
        slot_half: 5.0,
        slot_smoothing: 1.5,
        u_jet: 0.08, // overwritten per sweep point
        tau: 0.51,
        edge_distance: 30,
        plate_length: 40,
        fringe_width: 24,
        fringe_sigma: 0.3,
        steps_settle: 2500,
        steps_record: 4096,
        // NO nozzle: with a nozzle at h/delta = 3 the flow is STEADY
        // (executed: transverse force at machine noise 1e-15 — the
        // first exponent 'fit' was drift). The no-nozzle rig's
        // free-jet-mode oscillation is the broadband source this
        // table catalogs; it is NOT an edge-tone stage (see the
        // staging suite) and the table claims shape/scaling only.
        seed_amplitude: 0.02,
        nozzle_thickness: 0,
        collision: CollisionModel2::Bgk,
    }
}

#[test]
fn noise_table_sweep_export_and_synth_round_trip() {
    let sweep = [0.06, 0.08, 0.105];
    let table = fit_noise_table(&geometry(), &sweep).expect("table");
    assert_eq!(table.entries.len(), 3);
    // REAL amplitudes: the band-limited force RMS must be a
    // physical scale, not amplified roundoff (executed trap: the
    // unseeded symmetric rig read 1e-15 'oscillations' whose
    // spectral ratios passed every shape gate).
    for e in &table.entries {
        assert!(
            e.force_rms > 1.0e-6,
            "force RMS {:.3e} at machine-noise scale — vacuous run",
            e.force_rms
        );
    }
    // The dipole strength vs velocity is REPORTED, not prescribed:
    // the saturated limit-cycle amplitude of this low-Re tonal rig
    // is NOT monotone in u near mode switches (executed:
    // [4.9e-5, 3.5e-5, 8.0e-5] across the sweep — multi-stability,
    // consistent with the staging suite's hysteresis findings), and
    // a broadband U^6-style power law does not apply to a coherent
    // limit cycle. The table carries the measured exponent as DATA;
    // the assertion is only that the fit is finite.
    assert!(
        table.power_exponent.is_finite(),
        "exponent fit failed: {}",
        table.power_exponent
    );
    // Every entry's diagnostics are in regime.
    for e in &table.entries {
        assert!(e.mach_max_lattice < 0.25 && e.flux_imbalance < 0.02);
        // Shape sanity: normalized (max band = 0 dB) and non-flat.
        let max = e.band_db.iter().copied().fold(f64::MIN, f64::max);
        let min = e.band_db.iter().copied().fold(f64::MAX, f64::min);
        assert!(max.abs() < 1e-12 && min < -10.0, "degenerate shape");
    }
    // Export: the marketing-claim mutation guard — the JSON must
    // embed the scope statement verbatim.
    let json = table.to_json();
    assert!(
        json.contains(SCOPE_STATEMENT),
        "scope statement missing from export"
    );
    assert!(json.contains("\"power_exponent_measured\""));
    // Demo-synth round trip: synthesize at the middle sweep point
    // and verify the output's band spectrum matches the table shape.
    let pcm = table.synthesize(0.08, 16_384, 4242).expect("synth");
    assert_eq!(pcm.len(), 16_384);
    // Band analysis of the synthesized signal (same folding).
    let n = pcm.len();
    let fft = fs_fft::Fft::new(n);
    let mut buf: Vec<fs_fft::C64> = pcm.iter().map(|&v| fs_fft::C64::new(v, 0.0)).collect();
    let mut scratch = vec![fs_fft::C64::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    let delta = 2.0 * table.geometry.slot_half;
    // Fold with the SAME edge-based banding the synth applies
    // (executed: nearest-center folding disagreed by 10 dB in the
    // few-bin low bands).
    let mut band_pow = [0.0f64; N_BANDS];
    let mut band_bins = [0usize; N_BANDS];
    for (k, c) in buf[..n / 2].iter().enumerate().skip(1) {
        let st = (k as f64 / n as f64) * delta / 0.08;
        if let Some(b) = fs_aeroac::noisetable::band_of(st) {
            band_pow[b] += c.norm_sq();
            band_bins[b] += 1;
        }
    }
    // Density normalization, matching the table's convention.
    for (p, &m) in band_pow.iter_mut().zip(&band_bins) {
        if m > 0 {
            *p /= m as f64;
        }
    }
    let entry = &table.entries[1];
    let peak = band_pow.iter().copied().fold(f64::MIN, f64::max);
    let mut worst = 0.0f64;
    let mut checked = 0usize;
    for (b, (&p, &want_db)) in band_pow.iter().zip(&entry.band_db).enumerate() {
        // Skip deep-floor bands and bands with too few bins for a
        // stable chi-squared power estimate.
        if want_db < -60.0 || p <= 0.0 || band_bins[b] < 16 {
            continue;
        }
        checked += 1;
        let got_db = 10.0 * (p / peak).log10();
        let err = (got_db - want_db).abs();
        worst = worst.max(err);
        assert!(
            err < 6.0,
            "band {b}: synth {got_db:.1} dB vs table {want_db:.1} dB"
        );
    }
    assert!(
        checked >= 6,
        "round-trip must check enough bands: {checked}"
    );
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"noise-table\",\"power_exponent\":{:.2},\"rms\":{:?},\"synth_worst_band_db_err\":{worst:.2},\"verdict\":\"pass\"}}",
        table.power_exponent,
        table
            .entries
            .iter()
            .map(|e| e.force_rms)
            .collect::<Vec<_>>()
    );
}

/// Typed refusals of the table and the synth.
#[test]
fn noise_table_refusals_are_typed() {
    let g = geometry();
    assert!(matches!(
        fit_noise_table(&g, &[0.08]),
        Err(AeroacError::InvalidParameter { .. })
    ));
    assert!(matches!(
        fit_noise_table(&g, &[0.08, 0.06]),
        Err(AeroacError::InvalidParameter { .. })
    ));
    // Synth refusals against a hand-built minimal table.
    let table = fs_aeroac::noisetable::NoiseTable {
        entries: vec![],
        power_exponent: 0.0,
        geometry: g.clone(),
        scope: SCOPE_STATEMENT,
    };
    assert!(matches!(
        table.synthesize(0.08, 4096, 1),
        Err(AeroacError::InvalidParameter { .. })
    ));
}
