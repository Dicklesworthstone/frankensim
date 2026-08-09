//! E2E casebook: a TMM clarinet-class input impedance (fs-duct, first
//! principles: USSA-1976 air + Zwikker-Kosten losses + unflanged
//! radiation) lowered to a <=24-pole PASSIVE rational filter and
//! bilinear-discretized to an audio-rate biquad bank, with every stage
//! gated in musical units (cents at the impedance peaks) and logged as
//! JSON lines.
//!
//! This is the bridge the musical program needs: offline
//! frequency-domain physics -> runtime-realizable filter, with
//! passivity certified (an active radiation load in a waveguide
//! feedback loop CREATES energy).

use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::discretize::bilinear;
use fs_vfit::passivity::{check_passivity, repair_passivity};
use fs_vfit::vf::fit_auto_order;
use fs_vfit::{FitOptions, RationalModel};

use fs_duct::{Duct, LossModel, Segment, Termination, impedance_peaks, impedance_sweep};

const TWO_PI: f64 = 2.0 * core::f64::consts::PI;

/// Cylindrical clarinet-class bore: 14.6 mm bore, 570 mm, closed reed
/// end (impedance measured AT the closed end looking into the bore
/// with an open unflanged far end).
fn clarinet_duct() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 7.3e-3,
            length: 0.57,
        }],
    }
}

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

/// Refine the peak of `|f|` near a bracketing triple by golden-section
/// on a fine grid (0.05-cent steps): measurement infrastructure, no
/// model knowledge.
fn refine_peak(mut lo: f64, mut hi: f64, f: &dyn Fn(f64) -> f64) -> f64 {
    for _ in 0..64 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        if f(m1) > f(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
        if hi / lo < 1.0 + 1.0e-7 {
            break;
        }
    }
    0.5 * (lo + hi)
}

fn cents(a: f64, b: f64) -> f64 {
    1200.0 * fs_math::det::ln(a / b) / fs_math::det::ln(2.0)
}

#[test]
fn clarinet_impedance_to_passive_audio_filter() {
    let duct = clarinet_duct();
    let state = air();
    let (f_lo, f_hi) = (60.0, 3000.0);
    let sweep = impedance_sweep(
        &duct,
        &state,
        TWO_PI * f_lo,
        TWO_PI * f_hi,
        2400,
        LossModel::WideTube,
        Termination::UnflangedOpen,
    )
    .expect("sweep");
    let omega: Vec<f64> = sweep.iter().map(|r| r.omega).collect();
    let h: Vec<C64> = sweep.iter().map(|r| r.impedance).collect();
    // TMM reference peaks, refined to sub-cent.
    let peak_idx = impedance_peaks(&sweep);
    assert!(
        peak_idx.len() >= 8,
        "clarinet band should carry >= 8 odd-harmonic peaks, got {}",
        peak_idx.len()
    );
    let z_tmm = |w: f64| {
        fs_duct::input_impedance(
            &duct,
            &state,
            w,
            LossModel::WideTube,
            Termination::UnflangedOpen,
        )
        .expect("z")
        .impedance
        .abs()
    };
    let tmm_peaks: Vec<f64> = peak_idx
        .iter()
        .map(|&i| refine_peak(omega[i - 1], omega[i + 1], &z_tmm))
        .collect();
    // STAGE 1: rational lowering, ascending order, <= 24 poles.
    let base = FitOptions {
        fit_e: false,
        ..FitOptions::new(8)
    };
    let (fit, curve) =
        fit_auto_order(&omega, &h, &[8, 12, 16, 20, 24], &base, 0.2, 0.0).expect("fit");
    assert!(fit.model.order() <= 24, "order cap violated");
    assert!(fit.model.is_stable());
    println!(
        "{{\"suite\":\"fs-vfit-clarinet\",\"stage\":\"order-curve\",\"curve\":{:?},\"selected\":{}}}",
        curve,
        fit.model.order()
    );
    // Continuous-fit peaks within the authored cents gate.
    let z_fit = |w: f64| fit.model.eval_iw(w).abs();
    let mut worst_fit_cents = 0.0f64;
    for &wp in &tmm_peaks {
        let refined = refine_peak(wp * 0.97, wp * 1.03, &z_fit);
        worst_fit_cents = worst_fit_cents.max(cents(refined, wp).abs());
    }
    // Authored: measured 0.1-cent class at order 24 on clean TMM data;
    // 2 cents is the musical-transparency bar with wide headroom.
    assert!(
        worst_fit_cents < 2.0,
        "continuous fit peaks off by {worst_fit_cents:.3} cents"
    );
    // STAGE 2: passivity certification (+ repair if the raw fit is
    // active — either way the CERTIFIED model proceeds).
    let band = (TWO_PI * f_lo, TWO_PI * f_hi);
    let pre_report = check_passivity(&fit.model, band).expect("check");
    let (certified, margin_before, rounds): (RationalModel, f64, usize) = if pre_report.passive {
        (fit.model.clone(), pre_report.worst.0, 0)
    } else {
        let worst = pre_report.worst.0;
        let (repaired, rep) = repair_passivity(&fit.model, band).expect("repair");
        assert!(rep.certificate.passive);
        (repaired, worst, rep.rounds)
    };
    let post_report = check_passivity(&certified, band).expect("post check");
    assert!(post_report.passive, "certified model must be passive");
    println!(
        "{{\"suite\":\"fs-vfit-clarinet\",\"stage\":\"passivity\",\"class\":\"{:?}\",\"worst_re_before\":{margin_before:.6e},\"worst_re_after\":{:.6e},\"repair_rounds\":{rounds}}}",
        post_report.class, post_report.worst.0
    );
    // The certificate must not have cost musical fidelity: re-measure
    // the certified model's peaks against TMM.
    let z_cert = |w: f64| certified.eval_iw(w).abs();
    let mut worst_cert_cents = 0.0f64;
    for &wp in &tmm_peaks {
        let refined = refine_peak(wp * 0.97, wp * 1.03, &z_cert);
        worst_cert_cents = worst_cert_cents.max(cents(refined, wp).abs());
    }
    assert!(
        worst_cert_cents < 2.0,
        "certified model peaks off by {worst_cert_cents:.3} cents"
    );
    // STAGE 3: bilinear discretization at 192 kHz (the internal
    // oversampled rate physical-model synthesis runs at; bilinear warp
    // is (omega*T)^2/12 ~ 8e-4 relative at 3 kHz, i.e. ~1.4 cents
    // worst-case, split by prewarping mid-band).
    let fs_hz = 192_000.0;
    let t_s = 1.0 / fs_hz;
    let prewarp = fs_math::det::sqrt(TWO_PI * f_lo * TWO_PI * f_hi);
    let filt = bilinear(&certified, t_s, prewarp).expect("bilinear");
    assert!(filt.is_stable());
    // Peaks re-measured FROM THE DIGITAL FILTER.
    let z_dig = |w: f64| filt.eval(w).abs();
    let mut worst_dig_cents = 0.0f64;
    let mut peak_rows = String::new();
    for &wp in &tmm_peaks {
        let refined = refine_peak(wp * 0.97, wp * 1.03, &z_dig);
        let c = cents(refined, wp);
        worst_dig_cents = worst_dig_cents.max(c.abs());
        peak_rows.push_str(&format!(
            "{{\"f_tmm_hz\":{:.3},\"cents_dig\":{c:.4}}},",
            wp / TWO_PI
        ));
    }
    // Authored: measured sub-cent at 192 kHz with mid-band prewarp;
    // 2 cents keeps the same musical-transparency bar as the fit.
    assert!(
        worst_dig_cents < 2.0,
        "digital filter peaks off by {worst_dig_cents:.3} cents"
    );
    // f32 coefficient quantization: report the peak-band response
    // sensitivity (biquad banks ship to f32 DSP paths).
    let mut worst_f32 = 0.0f64;
    for &wp in &tmm_peaks {
        let dq = (filt.eval_f32_quantized(wp) - filt.eval(wp)).abs() / filt.eval(wp).abs();
        worst_f32 = worst_f32.max(dq);
    }
    println!(
        "{{\"suite\":\"fs-vfit-clarinet\",\"stage\":\"digital\",\"fs_hz\":{fs_hz},\"prewarp_rad\":{prewarp:.1},\"worst_fit_cents\":{worst_fit_cents:.4},\"worst_cert_cents\":{worst_cert_cents:.4},\"worst_dig_cents\":{worst_dig_cents:.4},\"worst_f32_rel\":{worst_f32:.3e},\"peaks\":[{peak_rows}],\"verdict\":\"pass\"}}"
    );
}
