//! Jet-labium fixture: the edge-tone base flow produces a REAL
//! self-sustained oscillation (prominent transverse-force spectral
//! peak), diagnostics stay in regime, the Curle radiation of the
//! measured force carries the scope statement, refusals are typed,
//! and the run is bitwise deterministic.

use fs_aeroac::jetlab::{JetLabiumConfig, dipole_spectrum_line, run_jet_labium};
use fs_aeroac::{AeroacError, SCOPE_STATEMENT};
use fs_math::c64::C64;

fn base_config() -> JetLabiumConfig {
    JetLabiumConfig {
        nx: 192,
        ny: 80,
        slot_half: 5.0,
        slot_smoothing: 1.5,
        u_jet: 0.08,
        tau: 0.51, // nu = 1/300 -> Re = 240
        edge_distance: 30,
        plate_length: 60,
        fringe_width: 32,
        fringe_sigma: 0.3,
        steps_settle: 3000,
        steps_record: 4096,
    }
}

/// Hann-windowed power spectrum of the transverse force.
fn fy_power_spectrum(force: &[[f64; 2]]) -> Vec<f64> {
    let n = force.len();
    let mean = force.iter().map(|f| f[1]).sum::<f64>() / n as f64;
    let fft = fs_fft_shim::Fft::new(n);
    let mut buf: Vec<fs_fft_shim::C64f> = force
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let w = 0.5 - 0.5 * ((2.0 * core::f64::consts::PI * i as f64) / (n as f64 - 1.0)).cos();
            fs_fft_shim::C64f::new((f[1] - mean) * w, 0.0)
        })
        .collect();
    let mut scratch = vec![fs_fft_shim::C64f::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    buf[..n / 2].iter().map(|c| c.norm_sq()).collect()
}

mod fs_fft_shim {
    pub use fs_fft::{C64 as C64f, Fft};
}

#[test]
fn jet_labium_edge_tone_oscillates_and_radiates() {
    let cfg = base_config();
    let run = run_jet_labium(&cfg).expect("run");
    assert_eq!(run.force_series.len(), 4096);
    // Regime diagnostics.
    let d = &run.diagnostics;
    assert!(
        d.mach_max_lattice < 0.25,
        "low-Mach diagnostic: {}",
        d.mach_max_lattice
    );
    assert!((d.reynolds - 240.0).abs() < 1.0, "Re: {}", d.reynolds);
    // Both measurement planes carry positive mean through-flux (the
    // fringe recycles the jet, it does not stall it).
    assert!(
        d.flux_plate_plane > 0.0 && d.flux_fringe_plane > 0.0,
        "through-flux: {} / {}",
        d.flux_plate_plane,
        d.flux_fringe_plane
    );
    // Flux imbalance between the plate and fringe planes: measured
    // 0.13% (the outlet-reflection pilot pathology read 6%) — gated
    // at 1%.
    let imbalance = (d.flux_plate_plane - d.flux_fringe_plane).abs() / d.flux_plate_plane.abs();
    assert!(imbalance < 0.01, "flux imbalance {imbalance:.4}");
    // The transverse force oscillates: a prominent spectral peak
    // (exclude the near-DC drift bins), measured against the median
    // power. Prominence and location printed for the record.
    let power = fy_power_spectrum(&run.force_series);
    let (peak_bin, peak_pow) = power
        .iter()
        .enumerate()
        .skip(4)
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("spectrum");
    let mut sorted = power[4..].to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = sorted[sorted.len() / 2];
    let prominence = peak_pow / median.max(1e-300);
    assert!(
        prominence > 100.0,
        "no self-sustained oscillation: peak/median = {prominence:.1}"
    );
    let freq = peak_bin as f64 / 4096.0; // cycles per step
    let strouhal = freq * cfg.edge_distance as f64 / cfg.u_jet;
    assert!(
        freq > 0.0 && strouhal.is_finite() && strouhal > 0.0,
        "St: {strouhal}"
    );
    // Radiate the peak line through the Curle dipole at a far
    // observer: finite pressure, scope statement embedded end-to-end.
    let k = 2.0 * core::f64::consts::PI * freq / (1.0 / 3.0f64.sqrt());
    let p = dipole_spectrum_line(
        [C64::new(0.0, 0.0), C64::new(peak_pow.sqrt(), 0.0)],
        k,
        [0.0, 60.0],
        [cfg.edge_distance as f64, 0.0],
    )
    .expect("radiation");
    assert!(p.re.is_finite() && p.im.is_finite() && p.abs() > 0.0);
    assert_eq!(run.scope, SCOPE_STATEMENT);
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"jet-labium\",\"peak_bin\":{peak_bin},\"freq_per_step\":{freq:.6},\"strouhal\":{strouhal:.4},\"prominence\":{prominence:.1},\"mach_max\":{:.4},\"flux_plate\":{:.4},\"flux_fringe\":{:.4},\"verdict\":\"pass\"}}",
        d.mach_max_lattice, d.flux_plate_plane, d.flux_fringe_plane
    );
}

/// Config refusals, typed by name.
#[test]
fn jet_labium_refusals_are_typed() {
    let mut c = base_config();
    c.steps_record = 1000; // not a power of two
    assert!(matches!(
        run_jet_labium(&c),
        Err(AeroacError::InvalidParameter { .. })
    ));
    let mut c = base_config();
    c.plate_length = 200; // reaches into the fringe
    assert!(matches!(
        run_jet_labium(&c),
        Err(AeroacError::InvalidParameter { .. })
    ));
    let mut c = base_config();
    c.u_jet = 0.5; // beyond low-Mach
    assert!(matches!(
        run_jet_labium(&c),
        Err(AeroacError::InvalidParameter { .. })
    ));
    let mut c = base_config();
    c.slot_half = 30.0; // slot taller than half the domain
    assert!(matches!(
        run_jet_labium(&c),
        Err(AeroacError::InvalidParameter { .. })
    ));
    let mut c = base_config();
    c.tau = f64::NAN;
    assert!(matches!(
        run_jet_labium(&c),
        Err(AeroacError::NonFinite { .. })
    ));
}

/// Bitwise determinism of a short run.
#[test]
fn jet_labium_determinism_bitwise() {
    let mut c = base_config();
    c.steps_settle = 100;
    c.steps_record = 64;
    let a = run_jet_labium(&c).expect("a");
    let b = run_jet_labium(&c).expect("b");
    for (x, y) in a.force_series.iter().zip(&b.force_series) {
        assert_eq!(x[0].to_bits(), y[0].to_bits());
        assert_eq!(x[1].to_bits(), y[1].to_bits());
    }
}
