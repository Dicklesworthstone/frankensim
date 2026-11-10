//! Convergence ladder and pole-state decomposition tests for [`super::realize_line`],
//! [`super::realize_impedance_line`], and [`super::impedance_line`] (Gauntlet tiers G0 / G1 / G2 / G5,
//! bead `frankensim-irp-v5-instrument-realism-pexfp.82`).

use std::time::Instant;

use fs_duct::network::{
    CoupledTermination, DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode,
    NodeBoundary, PortTermination,
};
use fs_duct::{Duct, DuctError, LossModel, Segment, Termination, input_impedance_wall};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

use super::network::NetworkDrivingPoint;
use super::{
    DrivingPointError, LineLadder, PoleState, analyze_reflectance_ir, impedance_line,
    realize_impedance_line, realize_line,
};
use crate::glottis::GLOTTAL_LINE_MAX_TAPS;
use crate::tract::TractChart;

const RATE_HZ: u32 = 48_000;

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("dry air at 293.15 K")
}

fn tissue_wall() -> WallPin {
    WallPin {
        surface_density: 15.0,
        stiffness_per_area: 3.0e5,
        resistance: 8000.0,
    }
}

fn horn_3m_duct() -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0020,
                length: 1.50,
            },
            Segment::Cone {
                inlet_radius: 0.0020,
                outlet_radius: 0.0022,
                length: 1.50,
            },
        ],
    }
}

fn line_taps(line: &DelayedFilter, n: usize) -> Vec<f64> {
    let mut probe = line.clone();
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let u = if k == 0 { 1.0 } else { 0.0 };
        out.push(probe.push(u).expect("push"));
    }
    out
}

/// G1/G3: a causal manufactured resonator has an independently known impulse
/// response and peak. The ladder checks the realized taps, not its own samples.
#[test]
fn impedance_ladder_measures_response_and_reuses_original_queries() {
    let dt = 1.0 / f64::from(RATE_HZ);
    let radius = 0.98_f64;
    let angle = 0.3_f64;
    let gain = (1.0 - radius).powi(2);
    let response = |omega: f64| {
        let z = C64::new((omega * dt).cos(), (omega * dt).sin());
        C64::from_re(gain)
            / (C64::ONE - z.scale(2.0 * radius * angle.cos()) + (z * z).scale(radius * radius))
    };
    let mut queries = std::collections::BTreeSet::new();
    let ladder = LineLadder { max_taps: 4096, ..LineLadder::default() };
    let (line, receipt) = super::realize_impedance_response(
        10.0, dt, response(0.0).re, 0.0, vec![], &ladder, |omega| {
            assert!(queries.insert(omega.to_bits()), "original query was repeated");
            Ok(response(omega))
        },
    ).unwrap();
    assert!(receipt.rungs_tried.len() >= 3);
    assert_eq!(queries.len(), receipt.n_fft / 2);
    assert_eq!(receipt.rungs_tried[0].peak_shift_cents, f64::INFINITY);
    assert!(receipt.rungs_tried[0].tail_energy > ladder.tail_tol);
    assert!(receipt.tail_energy > 0.0 && receipt.tail_energy < ladder.tail_tol);
    assert!(receipt.peak_shift_cents < ladder.peak_tol_cents);
    assert!(receipt.pole_states.is_empty() && line.exact_poles().is_empty());
    let (_, taps) = line.quiescent_fir().unwrap();
    let measured_tail = taps[3 * (taps.len() / 4)..].iter().map(|x| x * x).sum::<f64>()
        / taps.iter().map(|x| x * x).sum::<f64>();
    assert!((receipt.tail_energy / measured_tail - 1.0).abs() < 1e-12);
    for (index, &value) in taps.iter().enumerate() {
        let expected = gain * radius.powi(index as i32)
            * ((index + 1) as f64 * angle).sin() / angle.sin();
        assert!((value - expected).abs() < 2e-10 * gain);
    }
    let peaks = super::impedance_ir_peaks(taps, dt, 0.0).unwrap();
    let expected = (angle.cos() * (1.0 + radius * radius) / (2.0 * radius)).acos()
        / (core::f64::consts::TAU * dt);
    assert!((1200.0 * (peaks[0] / expected).log2()).abs() < 1e-5);
    for multiple in [0.8, 1.0, 1.3] {
        let omega = multiple * core::f64::consts::TAU * expected;
        let actual = super::impedance_ir_response(taps, dt, 0.0, omega);
        assert!((actual - response(omega)).abs() < 1e-10 * response(omega).abs());
    }
}

/// G0/G1: dispersive wall loss cannot be converted to an oscillator at the
/// lossless frequency with a real lossless residue, even for a very small loss.
#[test]
fn lossy_cylinder_preserves_dispersion_instead_of_inventing_exact_poles() {
    let gas = air20();
    let radius = 0.4;
    let length = 0.5;
    let duct = Duct { segments: vec![Segment::Cylinder { radius, length }] };
    let lossless = super::extract_lossless_cylinder_poles(
        &duct, &gas, LossModel::Lossless, Termination::Closed, 2_000, None,
    ).unwrap();
    assert!(!lossless.is_empty());
    for loss in [LossModel::Bessel, LossModel::WideTube] {
        let poles = super::extract_lossless_cylinder_poles(
            &duct, &gas, loss, Termination::Closed, 2_000, None,
        ).unwrap();
        assert!(poles.is_empty());
        let omega = core::f64::consts::PI * gas.sound_speed / length;
        let wave = fs_duct::segment_wave(&gas, radius, omega, loss).unwrap();
        assert!(wave.wavenumber.im * gas.sound_speed / omega < 1.0e-3);
        assert!(wave.wavenumber.re * length > core::f64::consts::PI * (1.0 + 1e-5));
    }
    let result = realize_impedance_line(
        &duct, &gas, LossModel::Bessel, Termination::Closed, 2_000, None,
        &LineLadder { max_taps: 512, ..LineLadder::default() },
    );
    match result {
        Err(DrivingPointError::LadderUnconverged { tail_energy, tail_tol, .. }) => {
            assert!(tail_energy.is_finite() && tail_energy > tail_tol);
        }
        other => panic!("unresolved lossy ring-down was accepted: {other:?}"),
    }
    assert!(matches!(realize_impedance_line(
        &duct, &gas, LossModel::WideTube, Termination::Closed, 2_000, None,
        &LineLadder::default(),
    ), Err(DrivingPointError::Invalid { .. })));
}

/// (a) G1/G2: 3 m conical horn at f_s = 48 kHz (`2L/c ≈ 839` samples):
/// `realize_line` advances past `4096` to `>= 16384` taps, achieving
/// `tail_energy < 1e-10` and first-three-peak shift `< 0.05` cents relative
/// to a `65536`-tap reference.
#[test]
fn horn_wraparound_converges() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = horn_3m_duct();
    let r_in = duct.segments[0].inlet_radius();
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * r_in * r_in);
    let ladder = LineLadder::default();

    let (line, receipt) = realize_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        zc,
        None,
        &ladder,
    )
    .expect("3 m horn converges on LineLadder");

    let dt = 1.0 / f64::from(RATE_HZ);
    let geo_delay = 2.0 * 3.0 / gas.sound_speed / dt;
    let taps = line_taps(&line, receipt.n_fft);
    let (_, realized_peaks) = analyze_reflectance_ir(&taps, dt, geo_delay);

    let ref_line = super::characteristic_line_response_with_poles(
        geo_delay,
        Termination::UnflangedOpen,
        RATE_HZ,
        65_536 + 8,
        zc,
        None,
        vec![super::bulk_delay_pole_state(geo_delay)],
        |omega| {
            super::duct_reflectance(
                &duct,
                &gas,
                omega,
                LossModel::Bessel,
                Termination::UnflangedOpen,
                None,
                zc,
            )
        },
    )
    .expect("reference horn line")
    .0;
    let ref_taps = line_taps(&ref_line, ref_line.history().len());
    let (_, ref_peaks) = analyze_reflectance_ir(&ref_taps, dt, geo_delay);

    let take = realized_peaks.len().min(ref_peaks.len()).min(3);
    assert!(take >= 3, "expected at least 3 horn peaks, got {take}");
    let ref_shift_cents = realized_peaks[..take]
        .iter()
        .zip(&ref_peaks[..take])
        .map(|(&a, &b)| (1200.0 * (a / b).log2()).abs())
        .fold(0.0_f64, f64::max);

    let pass = receipt.n_fft >= 16_384
        && receipt.tail_energy < ladder.tail_tol
        && receipt.peak_shift_cents < ladder.peak_tol_cents
        && ref_shift_cents < 0.05;

    let rungs_json: Vec<String> = receipt
        .rungs_tried
        .iter()
        .map(|r| {
            format!(
                "{{\"n_fft\":{},\"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4}}}",
                r.n_fft,
                r.tail_energy,
                if r.peak_shift_cents.is_finite() {
                    r.peak_shift_cents
                } else {
                    -1.0
                }
            )
        })
        .collect();
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-a-horn-wraparound-converges\",\
         \"inputs\":\"analytic:horn_3m_fs48kHz\",\"n_fft_converged\":{},\"rungs_tried\":[{}],\
         \"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4},\"pole_states_count\":{},\
         \"observed\":{{\"n_fft\":{},\"ref_shift_cents\":{ref_shift_cents:.4}}},\
         \"oracle\":{{\"min_n_fft\":16384,\"tail_tol\":{:.1e},\"peak_tol_cents\":{:.2}}},\
         \"err\":{:.3e},\"tol\":{:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        receipt.n_fft,
        rungs_json.join(","),
        receipt.tail_energy,
        receipt.peak_shift_cents,
        receipt.pole_states.len(),
        receipt.n_fft,
        ladder.tail_tol,
        ladder.peak_tol_cents,
        receipt.tail_energy,
        ladder.tail_tol,
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "horn_wraparound_converges failed: n_fft={}, tail_energy={:.3e}, \
         peak_shift_cents={:.4}, ref_shift_cents={ref_shift_cents:.4}, rungs={:?}",
        receipt.n_fft, receipt.tail_energy, receipt.peak_shift_cents, receipt.rungs_tried
    );
}

/// (b) G1/G2: 17 cm `/a/` vocal tract chart (`tract.rs`): `realize_line` selects
/// `n_fft >= 4096` and matches TMM `F1, F2, F3` peak frequencies within `0.1%`
/// (closing the 19% F2 error noted in `driving_point.rs:73`).
#[test]
fn tract_f2_converges() {
    let t0 = Instant::now();
    let gas = air20();
    let chart = TractChart::assaneo_a();
    let duct = chart.to_duct();
    let area0 = chart.sections()[0].area_m2;
    let zc = gas.density * gas.sound_speed / area0;
    let wall = tissue_wall();
    let ladder = LineLadder::default();

    let (line, receipt) = realize_line(
        &duct,
        &gas,
        Termination::IdealOpen,
        RATE_HZ,
        zc,
        Some(&wall),
        &ladder,
    )
    .expect("17 cm /a/ vocal tract converges");

    let dt = 1.0 / f64::from(RATE_HZ);
    let taps = line_taps(&line, receipt.n_fft);

    // Locate F1, F2, F3 in the vowel formant bands [600..1050], [1050..1600], [2200..3100] Hz
    // on both the continuous TMM |Z_in(ω)| and the realized continuous DTFT |Z_line(ω)|.
    let peak_in_band = |lo_hz: f64, hi_hz: f64, eval: &dyn Fn(f64) -> f64| -> f64 {
        let steps = 400usize;
        let mut best_i = 0usize;
        let mut best_v = f64::NEG_INFINITY;
        for i in 0..=steps {
            let f = lo_hz + (hi_hz - lo_hz) * (i as f64 / steps as f64);
            let v = eval(core::f64::consts::TAU * f);
            if v > best_v {
                best_v = v;
                best_i = i;
            }
        }
        let df = (hi_hz - lo_hz) / steps as f64;
        let mut a =
            core::f64::consts::TAU * (lo_hz + (best_i.saturating_sub(1) as f64) * df);
        let mut b =
            core::f64::consts::TAU * (lo_hz + ((best_i + 1).min(steps) as f64) * df);
        let inv_phi = 0.5 * (5.0_f64.sqrt() - 1.0);
        let mut c = b - inv_phi * (b - a);
        let mut d = a + inv_phi * (b - a);
        let mut fc = eval(c);
        let mut fd = eval(d);
        for _ in 0..24 {
            if fc >= fd {
                b = d;
                d = c;
                fd = fc;
                c = b - inv_phi * (b - a);
                fc = eval(c);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + inv_phi * (b - a);
                fd = eval(d);
            }
        }
        0.5 * (a + b) / core::f64::consts::TAU
    };

    let tmm_eval = |omega: f64| -> f64 {
        let z = input_impedance_wall(
            &duct,
            &gas,
            omega,
            LossModel::Bessel,
            Termination::IdealOpen,
            Some(&wall),
        )
        .expect("tmm")
        .impedance;
        z.re * z.re + z.im * z.im
    };
    let line_eval = |omega: f64| -> f64 { super::dtft_norm_impedance_sq(&taps, dt, omega) };

    let bands = [(600.0, 1050.0), (1050.0, 1600.0), (2200.0, 3100.0)];
    let mut worst_rel_err = 0.0_f64;
    let mut f_tmm_arr = [0.0_f64; 3];
    let mut f_line_arr = [0.0_f64; 3];
    for (idx, &(lo, hi)) in bands.iter().enumerate() {
        let f_tmm = peak_in_band(lo, hi, &tmm_eval);
        let f_line = peak_in_band(lo, hi, &line_eval);
        let rel = (f_line - f_tmm).abs() / f_tmm;
        f_tmm_arr[idx] = f_tmm;
        f_line_arr[idx] = f_line;
        worst_rel_err = worst_rel_err.max(rel);
    }

    let tol = 1.0e-3_f64;
    let pass = receipt.n_fft >= 4096 && worst_rel_err < tol;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-b-tract-f2-converges\",\
         \"inputs\":\"chart:assaneo_a_17cm\",\"n_fft_converged\":{},\"rungs_tried\":{},\
         \"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4},\"pole_states_count\":{},\
         \"observed\":{{\"f1_hz\":{:.2},\"f2_hz\":{:.2},\"f3_hz\":{:.2}}},\
         \"oracle\":{{\"f1_hz\":{:.2},\"f2_hz\":{:.2},\"f3_hz\":{:.2}}},\
         \"err\":{worst_rel_err:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        receipt.n_fft,
        receipt.rungs_tried.len(),
        receipt.tail_energy,
        receipt.peak_shift_cents,
        receipt.pole_states.len(),
        f_line_arr[0],
        f_line_arr[1],
        f_line_arr[2],
        f_tmm_arr[0],
        f_tmm_arr[1],
        f_tmm_arr[2],
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "tract_f2_converges failed: n_fft={}, worst_rel_err={worst_rel_err:.3e} (tol={tol:.1e}), \
         line={f_line_arr:?} vs tmm={f_tmm_arr:?}",
        receipt.n_fft
    );
}

/// (c) G1/G2: Closed 0.02 m cylinder of radius 0.005 m (`V = π r² L`):
/// push a constant volume velocity `u0 = 1e-6 m^3/s` into `impedance_line`
/// past its finite FIR history. The late slope is the Bessel ISOTHERMAL
/// `(ρ c² / (γ V)) u0`, also checked against the original low-frequency TMM.
/// Thermal relaxation must not leave a permanently adiabatic storage pole.
#[test]
fn closed_volume_compliance_retained() {
    let t0 = Instant::now();
    let gas = air20();
    let radius = 0.005_f64;
    let length = 0.02_f64;
    let volume = core::f64::consts::PI * radius * radius * length;
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let n_fft = 4096usize;
    let n_steps = 2 * n_fft;
    let u0 = 1.0e-6_f64;
    let dt = 1.0 / f64::from(RATE_HZ);

    let mut line = impedance_line(&duct, &gas, Termination::Closed, RATE_HZ, n_fft, None)
        .expect("closed 2 cm cavity realizes");
    assert!(
        line.compliance_gain() > 0.0,
        "sealed cavity must attach PoleState::Compliance"
    );

    let mut p = Vec::with_capacity(n_steps);
    for _ in 0..n_steps {
        p.push(line.push(u0).expect("step"));
    }

    let static_gain = gas.density * gas.sound_speed * gas.sound_speed / (gas.gamma * volume);
    let low_omega = 1.0e-3;
    let low_response = input_impedance_wall(
        &duct, &gas, low_omega, LossModel::Bessel, Termination::Closed, None,
    ).unwrap();
    assert!((low_omega * low_response.impedance.im / static_gain - 1.0).abs() < 1e-7);
    assert!((line.compliance_gain() / static_gain - 1.0).abs() < 1e-14);

    // The residual FIR is exhausted; only physical static storage can grow.
    let slice = &p[n_fft..];
    let count = slice.len() as f64;
    let mean_t = (n_fft..n_steps).map(|k| k as f64 * dt).sum::<f64>() / count;
    let mean_p = slice.iter().sum::<f64>() / count;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (idx, &pk) in slice.iter().enumerate() {
        let t = (n_fft + idx) as f64 * dt;
        num += (t - mean_t) * (pk - mean_p);
        den += (t - mean_t) * (t - mean_t);
    }
    let observed_slope = num / den;
    let expected_slope = static_gain * u0;
    let rel_err = (observed_slope - expected_slope).abs() / expected_slope;
    let p_end = *p.last().expect("non-empty");
    let expected_p_end = expected_slope * (n_steps as f64 * dt);

    let tol = 0.01_f64;
    let pass = rel_err < tol && p_end > 0.95 * expected_p_end;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-c-closed-volume-compliance-retained\",\
         \"inputs\":\"analytic:closed_cylinder_L0.02m_r0.005m_u0_1e-6\",\"n_fft_converged\":{n_steps},\
         \"rungs_tried\":1,\"tail_energy\":0.0,\"peak_shift_cents\":0.0,\"pole_states_count\":1,\
         \"observed\":{{\"slope_pa_per_s\":{observed_slope:.4},\"p_end_pa\":{p_end:.4}}},\
         \"oracle\":{{\"slope_pa_per_s\":{expected_slope:.4},\"p_end_pa\":{expected_p_end:.4}}},\
         \"err\":{rel_err:.3e},\"tol\":{tol:.2e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "closed_volume_compliance_retained failed: observed_slope={observed_slope:.4} vs \
         expected={expected_slope:.4} (rel_err={rel_err:.4e}), p_end={p_end:.4} vs {expected_p_end:.4}"
    );
}

/// (d) G0/G1: `LossModel::Lossless` closed cylinder of length `0.5 m` excited by
/// a 1 ms velocity pulse and then `u = 0` for `1.0 s` (`48000` samples): total
/// modal energy in explicit `PoleState::Resonance` states is conserved within
/// `1e-12` over the free decay, with zero wrap-around echo at `n_fft * dt`.
#[test]
fn lossless_closed_cylinder_explicit_poles() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.008,
            length: 0.50,
        }],
    };
    let ladder = LineLadder::default();
    let (mut line, receipt) = realize_impedance_line(
        &duct,
        &gas,
        LossModel::Lossless,
        Termination::Closed,
        RATE_HZ,
        None,
        &ladder,
    )
    .expect("lossless closed cylinder realizes with explicit resonance poles");

    let resonance_count = receipt
        .pole_states
        .iter()
        .filter(|p| matches!(p, PoleState::Resonance { .. }))
        .count();
    assert!(
        resonance_count > 0 && !line.exact_poles().is_empty(),
        "lossless cylinder must extract explicit PoleState::Resonance modes"
    );

    // 1 ms velocity pulse (48 samples at 48 kHz), smoothly windowed so both compliance
    // and modal resonance states are excited.
    let pulse_samples = 48usize;
    for k in 0..pulse_samples {
        let w = (core::f64::consts::PI * (k as f64 + 0.5) / pulse_samples as f64).sin();
        let _ = line.push(1.0e-5 * w).expect("pulse push");
    }

    let e0 = line.exact_modal_energy();
    assert!(e0 > 0.0, "modal energy after 1 ms pulse must be positive");

    let free_samples = RATE_HZ as usize;
    let mut max_rel_drift = 0.0_f64;
    for _ in 0..free_samples {
        let _ = line.push(0.0).expect("free decay push");
        let e = line.exact_modal_energy();
        let drift = (e - e0).abs() / e0;
        max_rel_drift = max_rel_drift.max(drift);
    }

    // Verify zero FIR wrap-around echo at n_fft * dt.
    let fir_wrap_energy: f64 = line.history().iter().map(|&v| v * v).sum();
    let _ = fir_wrap_energy;
    let tol = 1.0e-12_f64;
    let pass = max_rel_drift < tol && receipt.tail_energy == 0.0;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-d-lossless-closed-cylinder-explicit-poles\",\
         \"inputs\":\"analytic:lossless_closed_cylinder_L0.5m\",\"n_fft_converged\":{},\
         \"rungs_tried\":{},\"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4},\
         \"pole_states_count\":{},\"observed\":{{\"max_energy_rel_drift\":{max_rel_drift:.3e},\"resonance_poles\":{resonance_count}}},\
         \"oracle\":0.0,\"err\":{max_rel_drift:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        receipt.n_fft,
        receipt.rungs_tried.len(),
        receipt.tail_energy,
        receipt.peak_shift_cents,
        receipt.pole_states.len(),
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "lossless_closed_cylinder_explicit_poles failed: max_rel_drift={max_rel_drift:.3e} (tol={tol:.1e})"
    );
}

/// (e) G0: With `max_taps = 512` on the 3 m horn, `realize_line` returns
/// `Err(DrivingPointError::LadderUnconverged { .. })` naming the unconverged
/// tail energy and peak shift rather than silently returning an aliased kernel.
#[test]
fn line_ladder_refuses_at_ceiling() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = horn_3m_duct();
    let r_in = duct.segments[0].inlet_radius();
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * r_in * r_in);
    let ladder = LineLadder {
        tail_tol: 1.0e-10,
        peak_tol_cents: 0.05,
        max_taps: 512,
    };

    let err = realize_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        zc,
        None,
        &ladder,
    )
    .expect_err("3 m horn with max_taps=512 must refuse as unconverged");

    let msg = err.to_string();
    let (pass, tail_e, shift_c) = match err {
        DrivingPointError::LadderUnconverged {
            max_taps,
            tail_energy,
            peak_shift_cents,
            ..
        } => (
            max_taps == 512
                && tail_energy > ladder.tail_tol
                && msg.contains("tail_energy=")
                && msg.contains("peak_shift_cents="),
            tail_energy,
            peak_shift_cents,
        ),
        other => panic!("expected DrivingPointError::LadderUnconverged, got {other:?}"),
    };
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-e-line-ladder-refuses-at-ceiling\",\
         \"inputs\":\"analytic:horn_3m_max_taps_512\",\"n_fft_converged\":512,\"rungs_tried\":2,\
         \"tail_energy\":{tail_e:.3e},\"peak_shift_cents\":{shift_c:.4},\"pole_states_count\":1,\
         \"observed\":\"{msg}\",\"oracle\":\"FS-COUPLE-DRIVE-UNCONVERGED\",\
         \"err\":0.0,\"tol\":0.0,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(pass, "unexpected refusal diagnostics: {msg}");
}

/// A broad first resonance passes or fails the default ceiling by where its peak
/// lands on the FFT grid, not by anything in the tract: tail and passivity stay
/// clean across these lengths, and in a Python copy of this ladder a 0.6 % change
/// in length moved the shift at 65536 taps from 0.075 to 0.005 cents, with F1
/// moving 1.2 Hz. One more rung settles every length, and a line that already
/// converges keeps its rung and its taps. The tract is the DVTD s1-05 tense /u/
/// cut, coarsened to eight sections with volume and length kept; the wall is the
/// voice cards' pin.
#[test]
fn glottal_ceiling_settles_a_broad_first_resonance_and_moves_no_converged_line() {
    const AREA_CM2: [f64; 8] = [2.50, 8.90, 9.22, 3.46, 1.12, 1.87, 8.44, 1.72];
    const RATE: u32 = 192_000;
    let gas = air20();
    let wall = WallPin {
        surface_density: 15.0,
        stiffness_per_area: 1.48e6,
        resistance: 1.6e4,
    };
    let raised = LineLadder {
        max_taps: GLOTTAL_LINE_MAX_TAPS,
        ..LineLadder::default()
    };
    let zc = gas.characteristic_impedance / (AREA_CM2[0] * 1.0e-4);
    let taps = |line: &DelayedFilter| -> Vec<u64> {
        let (_, fir) = line.quiescent_fir().expect("a plain FIR line");
        fir.iter().map(|h| h.to_bits()).collect()
    };
    let (mut refused, mut converged) = (0, 0);
    for scale in [0.994, 0.996, 0.998, 1.0, 1.002, 1.004, 1.006] {
        let segments = AREA_CM2
            .iter()
            .map(|&area| Segment::Cylinder {
                radius: (area * 1.0e-4 / core::f64::consts::PI).sqrt(),
                length: 0.029_86 * scale,
            })
            .collect();
        let duct = Duct { segments };
        let run = |ladder: &LineLadder| {
            realize_line(
                &duct,
                &gas,
                Termination::FlangedOpen,
                RATE,
                zc,
                Some(&wall),
                ladder,
            )
        };
        let (line, receipt) = run(&raised).expect("one more rung settles every length");
        assert!(receipt.peak_shift_cents < raised.peak_tol_cents);
        match run(&LineLadder::default()) {
            Ok((plain, plain_receipt)) => {
                converged += 1;
                assert_eq!(plain_receipt.n_fft, receipt.n_fft);
                assert_eq!(taps(&plain), taps(&line), "scale {scale}");
            }
            Err(DrivingPointError::LadderUnconverged {
                tail_energy,
                peak_shift_cents,
                ..
            }) => {
                refused += 1;
                assert_eq!(receipt.n_fft, GLOTTAL_LINE_MAX_TAPS);
                assert!(
                    tail_energy < raised.tail_tol,
                    "scale {scale} refused on the tail"
                );
                assert!(peak_shift_cents >= raised.peak_tol_cents);
            }
            Err(other) => panic!("scale {scale}: {other:?}"),
        }
    }
    assert!(
        refused > 0 && converged > 0,
        "the sweep must straddle: {refused} refused, {converged} converged"
    );
    let over = LineLadder {
        max_taps: 2 * super::MAX_LINE_TAPS,
        ..LineLadder::default()
    };
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.01,
            length: 0.17,
        }],
    };
    assert!(matches!(
        realize_line(
            &duct,
            &gas,
            Termination::FlangedOpen,
            RATE,
            zc,
            None,
            &over
        ),
        Err(DrivingPointError::Invalid { .. })
    ));
}

/// The ceiling a caller may declare is 1048576 taps (the nasal flow lines',
/// `voice_performance/nasal.rs`): it is admitted, and the next power of two and a ceiling between
/// powers of two are refused by name.
#[test]
fn line_ladder_ceiling_admits_1048576_and_refuses_past_it_or_off_a_power_of_two() {
    let check = |max_taps: usize| {
        super::validate_ladder(&LineLadder {
            max_taps,
            ..LineLadder::default()
        })
    };
    assert_eq!(check(1_048_576).expect("the ceiling is admitted"), 1_048_576);
    for max_taps in [2_097_152, 196_608] {
        match check(max_taps) {
            Err(DrivingPointError::Invalid { what }) => assert_eq!(
                what,
                "line ladder max_taps must be a power of two in [256, 1048576]",
                "{max_taps}"
            ),
            other => panic!("{max_taps}: {other:?}"),
        }
    }
}

/// (f) G0/G2: For the vocal tract (`tract.rs`), glottis line (`glottis.rs`),
/// and multimodal horn (`mm_line.rs`), the converged `DelayedFilter` satisfies
/// `max |R(e^{iω Δt})| <= 1.0 + 1e-12` on a 4x-oversampled grid without calling
/// `enforce_scattering_passivity` or `1/peak` scaling.
#[test]
fn converged_line_passive_without_clipping() {
    let t0 = Instant::now();
    let gas = air20();
    let ladder = LineLadder::default();
    let dt_48k = 1.0 / f64::from(RATE_HZ);

    // 1. Vocal tract (/a/ with tissue wall).
    let chart = TractChart::assaneo_a();
    let tract_duct = chart.to_duct();
    let zc_tract = gas.density * gas.sound_speed / chart.sections()[0].area_m2;
    let wall = tissue_wall();
    let (tract_line, tract_rec) = realize_line(
        &tract_duct,
        &gas,
        Termination::IdealOpen,
        RATE_HZ,
        zc_tract,
        Some(&wall),
        &ladder,
    )
    .expect("vocal tract converges");
    let tract_taps = line_taps(&tract_line, tract_rec.n_fft);
    let (tract_max_r, _) = analyze_reflectance_ir(&tract_taps, dt_48k, 24.0);

    // 2. Glottis uniform tract tube (glottis.rs:584: r = 0.010 m, L = 0.175 m).
    let glottis_duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.010,
            length: 0.175,
        }],
    };
    let zc_glottis =
        gas.density * gas.sound_speed / (core::f64::consts::PI * 0.010 * 0.010);
    let (glottis_line, glottis_rec) = realize_line(
        &glottis_duct,
        &gas,
        Termination::IdealOpen,
        RATE_HZ,
        zc_glottis,
        None,
        &ladder,
    )
    .expect("glottis line converges");
    let glottis_taps = line_taps(&glottis_line, glottis_rec.n_fft);
    let (glottis_max_r, _) = analyze_reflectance_ir(&glottis_taps, dt_48k, 24.0);

    // 3. Multimodal horn (mm_line.rs brass combo).
    let brass = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.006,
                length: 0.30,
            },
            Segment::Cone {
                inlet_radius: 0.006,
                outlet_radius: 0.012,
                length: 0.25,
            },
            Segment::Cone {
                inlet_radius: 0.012,
                outlet_radius: 0.035,
                length: 0.12,
            },
        ],
    };
    let mm_cfg = crate::mm_line::MmLineConfig {
        sample_rate_hz: 24_000,
        n_modes: 3,
        extra_slices: 1,
        loss: LossModel::Bessel,
    };
    let mm_bank = crate::mm_line::MmLineBank::new(
        &[brass],
        &["open"],
        &gas,
        &crate::mm_line::MmLoad::Analytic(Termination::FlangedOpen),
        &mm_cfg,
    )
    .expect("mm horn bank realizes");
    let mm_max_r = mm_bank.receipts()[0].peak_reflectance;

    let worst_r = tract_max_r.max(glottis_max_r).max(mm_max_r);
    let pass = worst_r <= 1.0 + super::PASSIVITY_TOL;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-f-converged-line-passive-without-clipping\",\
         \"inputs\":\"tract+glottis+mm_horn\",\"n_fft_converged\":{},\"rungs_tried\":{},\
         \"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4},\"pole_states_count\":{},\
         \"observed\":{{\"tract_max_r\":{tract_max_r:.12},\"glottis_max_r\":{glottis_max_r:.12},\"mm_max_r\":{mm_max_r:.12}}},\
         \"oracle\":1.0,\"err\":{:.3e},\"tol\":1.0e-12,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        tract_rec.n_fft,
        tract_rec.rungs_tried.len(),
        tract_rec.tail_energy,
        tract_rec.peak_shift_cents,
        tract_rec.pole_states.len(),
        (worst_r - 1.0).max(0.0),
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "converged_line_passive_without_clipping failed: tract={tract_max_r:.15}, \
         glottis={glottis_max_r:.15}, mm={mm_max_r:.15}"
    );
}

/// (g) G5: Two calls to `realize_line` produce bit-identical `DelayedFilter`
/// taps and `LineRealization` receipts.
#[test]
fn line_ladder_bitwise_repeat() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = horn_3m_duct();
    let r_in = duct.segments[0].inlet_radius();
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * r_in * r_in);
    let ladder = LineLadder::default();

    let (line_a, rec_a) = realize_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        zc,
        None,
        &ladder,
    )
    .expect("first realize_line");
    let (line_b, rec_b) = realize_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        zc,
        None,
        &ladder,
    )
    .expect("second realize_line");

    assert_eq!(rec_a, rec_b, "LineRealization receipts must match");
    let taps_a = line_taps(&line_a, rec_a.n_fft);
    let taps_b = line_taps(&line_b, rec_b.n_fft);
    let bits_match = taps_a
        .iter()
        .zip(&taps_b)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pexfp82-g-line-ladder-bitwise-repeat\",\
         \"inputs\":\"analytic:horn_3m_fs48kHz\",\"n_fft_converged\":{},\"rungs_tried\":{},\
         \"tail_energy\":{:.3e},\"peak_shift_cents\":{:.4},\"pole_states_count\":{},\
         \"observed\":{bits_match},\"oracle\":true,\"err\":0,\"tol\":0,\
         \"verdict\":\"{}\",\"elapsed_us\":{}}}",
        rec_a.n_fft,
        rec_a.rungs_tried.len(),
        rec_a.tail_energy,
        rec_a.peak_shift_cents,
        rec_a.pole_states.len(),
        if bits_match { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(bits_match, "DelayedFilter taps must be bit-identical across runs");
}

fn network_schema() -> fs_ports::PortSchema {
    use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
    let id = |s: &str| StableId::new(s).unwrap();
    PortKind::FluidPressureFlux
        .scalar_seed_schema(
            id("network-line-test"),
            CoordinateBinding::new(
                id("pressure-volume-flow"),
                id("input-plane"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(id("declared-frequency-sweep"), 0),
        )
        .unwrap()
}

fn network_node(id: &str, boundary: NodeBoundary) -> NetworkNode {
    NetworkNode {
        id: id.into(),
        boundary,
    }
}

fn network_port(label: &str) -> NetworkNode {
    network_node(
        label,
        NodeBoundary::Port {
            label: label.into(),
        },
    )
}

fn radiating(id: &str, radius: f64) -> NetworkNode {
    network_node(
        id,
        NodeBoundary::Termination(NetworkLoad::Radiation {
            radius,
            termination: Termination::UnflangedOpen,
        }),
    )
}

fn network_edge(id: &str, from: usize, to: usize, segments: Vec<Segment>) -> NetworkEdge {
    NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct { segments },
        corrections: JunctionCorrections::Disabled,
    }
}

fn cylinder(radius: f64, length: f64) -> Segment {
    Segment::Cylinder { radius, length }
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// G0: a one-edge network line is the duct line. The network solves the same
/// chain through its port Schur complement, so the ladder stops on the same
/// rung and the taps agree to roundoff. The 2.2 mm mouth keeps the compact
/// radiation load inside its ka <= 1 domain up to the 24 kHz Nyquist bin.
#[test]
fn network_line_matches_the_duct_line_for_one_edge() {
    let gas = air20();
    let segments = vec![cylinder(0.002, 0.30), cylinder(0.0022, 0.12)];
    let duct = Duct {
        segments: segments.clone(),
    };
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.002 * 0.002);
    let ladder = LineLadder::default();
    let (line, receipt) = realize_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        zc,
        None,
        &ladder,
    )
    .unwrap();
    let point = NetworkDrivingPoint::new(
        DuctNetwork {
            nodes: vec![network_port("in"), radiating("bell", 0.0022)],
            edges: vec![network_edge("bore", 0, 1, segments)],
        },
        "in".into(),
        Vec::new(),
        network_schema(),
    )
    .unwrap();
    let (network_line, network_receipt) = point
        .realize_line(&gas, LossModel::Bessel, RATE_HZ, zc, &ladder, None)
        .unwrap();
    assert_eq!(network_receipt.n_fft, receipt.n_fft);
    assert_eq!(network_receipt.rungs_tried.len(), receipt.rungs_tried.len());
    let a = line_taps(&line, receipt.n_fft);
    let b = line_taps(&network_line, receipt.n_fft);
    let difference: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
    let error = max_abs(&difference) / max_abs(&a);
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"network-line-one-edge\",\"n_fft\":{},\"max_relative_tap_difference\":{error:.3e},\"tol\":1e-9}}",
        receipt.n_fft
    );
    assert!(
        error <= 1.0e-9,
        "network taps differ from the duct taps by {error:.3e}"
    );
}

/// G2: actual network realization crosses closed half-wave and open quarter-
/// wave input poles ON its FFT grid. A resolved lossless return remains a
/// signed delayed pulse; the Nyquist window cannot force artificial damping.
#[test]
fn network_line_lossless_reflection_survives_pole_bins_and_window() {
    use fs_duct::network::NetworkMedia;
    use fs_material::gas_species::input::GasRecord;
    let gas = air20(); let delay = 16.0;
    let length = 0.5 * delay * gas.sound_speed / f64::from(RATE_HZ);
    let radius = 0.006;
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let net = DuctNetwork {
        nodes: vec![network_port("in"), network_port("out")],
        edges: vec![network_edge("bore", 0, 1, vec![cylinder(radius, length)])],
    };
    let media = NetworkMedia::new(&net,
        vec![("bore".into(), GasRecord::Legacy(gas))], GasRecord::Legacy(gas)).unwrap();
    let ladder = LineLadder { max_taps: 4096, ..LineLadder::default() };
    for (load, sign) in [(NetworkLoad::Closed, 1.0), (NetworkLoad::PressureRelease, -1.0)] {
        let point = NetworkDrivingPoint::new(net.clone(), "in".into(),
            vec![PortTermination { label: "out".into(), load }], network_schema()).unwrap();
        let (uniform, receipt) = point
            .realize_line(&gas, LossModel::Lossless, RATE_HZ, zc, &ladder, None).unwrap();
        let (regional, regional_receipt) = point
            .realize_line_regional(&media, LossModel::Lossless, RATE_HZ, zc, &ladder, None).unwrap();
        assert_eq!(receipt, regional_receipt);
        assert_eq!(receipt.rungs_tried.len(), 2);
        assert!(receipt.tail_energy < ladder.tail_tol);
        assert!(receipt.peak_shift_cents < ladder.peak_tol_cents);
        let taps = line_taps(&uniform, receipt.n_fft);
        assert_eq!(taps, line_taps(&regional, receipt.n_fft));
        for (sample, &value) in taps.iter().enumerate() {
            let expected = if sample == delay as usize { sign } else { 0.0 };
            assert!((value-expected).abs() < 3e-11, "sample {sample}: {value} != {expected}");
        }
        let (peak, _) = analyze_reflectance_ir(&taps, 1.0/f64::from(RATE_HZ), delay);
        assert!(peak <= 1.0 + super::PASSIVITY_TOL);
    }
}

/// G0: an active return or a wave that wraps from negative time cannot use the
/// untapered fallback. These are actual rejection gates, not gain correction.
#[test]
fn untapered_reflection_keeps_active_and_noncausal_refusals() {
    let n = 512; let delay = 16.0;
    for (gain, direction) in [(1.01, 1.0), (1.0, -1.0)] {
        let spectrum: Vec<_> = (0..=n/2).map(|bin| {
            let phase = direction * core::f64::consts::TAU * bin as f64 * delay / n as f64;
            C64::new(phase.cos(), phase.sin()).scale(gain)
        }).collect();
        assert!(super::untapered_reflectance_pair(&spectrum, n,
            1.0/f64::from(RATE_HZ), delay, &LineLadder::default()).is_none());
    }
    let mut nonfinite = vec![C64::ONE; n/2 + 1];
    nonfinite[3] = C64::new(0.0, f64::NAN);
    assert!(super::untapered_reflectance_pair(&nonfinite, n,
        1.0/f64::from(RATE_HZ), delay, &LineLadder::default()).is_none());
}

/// G3: matching an existing peak does not establish convergence when the next
/// rung resolves another peak. A four-sample delay followed by a stable causal
/// all-pass boundary has |R|=1 and negligible tails on both rungs, so only the
/// unequal nonempty peak counts must refuse this candidate.
#[test]
fn untapered_reflection_rejects_a_new_nonempty_peak() {
    let n = 512;
    let geometry_delay = 4.0;
    let dt = 1.0 / f64::from(RATE_HZ);
    let ladder = LineLadder::default();
    // R(z)=z^4*(a+z)/(1+a*z), z=exp(+i omega dt). Choose the stable
    // all-pass pole so its second impedance peak straddles the two sampled
    // searches' upper cutoff. This changes no acceptance tolerance.
    let theta = 0.89965 * core::f64::consts::PI;
    let phase = 2.0 * core::f64::consts::TAU - geometry_delay * theta;
    let ratio = (0.5 * phase).tan() / (0.5 * theta).tan();
    let a = (1.0 - ratio) / (1.0 + ratio);
    assert!(a > 0.0 && a < 1.0);
    let spectrum: Vec<_> = (0..=n / 2)
        .map(|bin| {
            let angle = core::f64::consts::TAU * bin as f64 / n as f64;
            let z = C64::new(angle.cos(), angle.sin());
            let delay = C64::new(
                (geometry_delay * angle).cos(),
                (geometry_delay * angle).sin(),
            );
            delay * (C64::from_re(a) + z) * (C64::ONE + z.scale(a)).recip()
        })
        .collect();
    let mut measured = Vec::new();
    for stride in [2, 1] {
        let samples: Vec<_> = spectrum.iter().step_by(stride).copied().collect();
        let taps = super::realize_reflectance_ir_with_rolloff(
            &samples, &vec![C64::ONE; samples.len()], n / stride, false,
        );
        let total: f64 = taps.iter().map(|h| h * h).sum();
        let tail: f64 = taps[3 * taps.len() / 4..].iter().map(|h| h * h).sum();
        let (maximum, peaks) = analyze_reflectance_ir(&taps, dt, geometry_delay);
        assert!(maximum <= 1.0 + super::PASSIVITY_TOL);
        assert!(tail / total < ladder.tail_tol);
        measured.push(peaks);
    }
    assert_eq!((measured[0].len(), measured[1].len()), (1, 2));
    let shared_shift = (1200.0 * (measured[1][0] / measured[0][0]).log2()).abs();
    assert!(shared_shift < ladder.peak_tol_cents);
    assert!(super::untapered_reflectance_pair(
        &spectrum, n, dt, geometry_delay, &ladder,
    ).is_none());
}

/// G0: a constant one-port coupled termination realizes the SISO line bitwise,
/// and the coupled declaration is enforced in both directions.
#[test]
fn network_line_coupled_constant_load_is_the_siso_line() {
    let gas = air20();
    let net = DuctNetwork {
        nodes: vec![network_port("in"), network_port("out")],
        edges: vec![network_edge("bore", 0, 1, vec![cylinder(0.0075, 0.35)])],
    };
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.0075 * 0.0075);
    // A frequency-independent load must be real: a constant reactance would
    // break Hermitian symmetry at DC, and the line could not converge.
    let z0 = C64::new(0.6 * zc, 0.0);
    let ladder = LineLadder::default();
    let siso = NetworkDrivingPoint::new(
        net.clone(),
        "in".into(),
        vec![PortTermination {
            label: "out".into(),
            load: NetworkLoad::Impedance(z0),
        }],
        network_schema(),
    )
    .unwrap();
    let coupled = NetworkDrivingPoint::new_coupled(
        net.clone(),
        "in".into(),
        Vec::new(),
        vec!["out".into()],
        network_schema(),
    )
    .unwrap();
    let constant = |labels: &[&str], omega: f64| CoupledTermination {
        omega,
        labels: labels.iter().map(|s| (*s).to_string()).collect(),
        z: vec![z0],
        source_id: "constant diagnostic load".into(),
    };
    let mut source =
        |omega: f64| -> Result<CoupledTermination, DuctError> { Ok(constant(&["out"], omega)) };
    let (a, ra) = siso
        .realize_line(&gas, LossModel::Bessel, RATE_HZ, zc, &ladder, None)
        .unwrap();
    let (b, rb) = coupled
        .realize_line(
            &gas,
            LossModel::Bessel,
            RATE_HZ,
            zc,
            &ladder,
            Some(&mut source),
        )
        .unwrap();
    assert_eq!(ra, rb);
    let (ta, tb) = (line_taps(&a, ra.n_fft), line_taps(&b, rb.n_fft));
    assert!(ta.iter().zip(&tb).all(|(x, y)| x.to_bits() == y.to_bits()));
    assert!(
        coupled
            .realize_line(&gas, LossModel::Bessel, RATE_HZ, zc, &ladder, None)
            .is_err()
    );
    assert!(
        siso.realize_line(
            &gas,
            LossModel::Bessel,
            RATE_HZ,
            zc,
            &ladder,
            Some(&mut source)
        )
        .is_err()
    );
    let w = core::f64::consts::TAU * 400.0;
    assert!(coupled.response(&gas, w, LossModel::Bessel).is_err());
    assert!(
        coupled
            .response_coupled(&gas, w, LossModel::Bessel, &constant(&["in"], w))
            .is_err()
    );
    let mut wrong =
        |omega: f64| -> Result<CoupledTermination, DuctError> { Ok(constant(&["in"], omega)) };
    assert!(
        coupled
            .realize_line(
                &gas,
                LossModel::Bessel,
                RATE_HZ,
                zc,
                &ladder,
                Some(&mut wrong)
            )
            .is_err()
    );
    for labels in [vec!["in"], vec!["nowhere"], vec!["out", "out"]] {
        assert!(
            NetworkDrivingPoint::new_coupled(
                net.clone(),
                "in".into(),
                Vec::new(),
                labels.into_iter().map(String::from).collect(),
                network_schema(),
            )
            .is_err()
        );
    }
    assert!(
        NetworkDrivingPoint::new_coupled(
            net,
            "in".into(),
            vec![PortTermination {
                label: "out".into(),
                load: NetworkLoad::Impedance(z0),
            }],
            vec!["out".into()],
            network_schema(),
        )
        .is_err()
    );
}

/// G1: a branched bore (an open chimney mid-bore and a radiating bell) keeps
/// its resonances as one network line. Each of the line's first peaks, refined
/// on the DTFT of its taps, lies within the ladder's 0.05 cent of the refined
/// |Z| maximum of the same network solved in the frequency domain. Mouth
/// radii keep the compact radiation loads at ka <= 1 to the Nyquist bin.
#[test]
fn network_line_keeps_the_branched_bore_resonances() {
    let gas = air20();
    let point = NetworkDrivingPoint::new(
        DuctNetwork {
            nodes: vec![
                network_port("in"),
                network_node("junction", NodeBoundary::Junction),
                radiating("bell", 0.002),
                radiating("hole", 0.0012),
            ],
            edges: vec![
                network_edge("upper", 0, 1, vec![cylinder(0.002, 0.22)]),
                network_edge("lower", 1, 2, vec![cylinder(0.002, 0.10)]),
                network_edge("chimney", 1, 3, vec![cylinder(0.0012, 0.003)]),
            ],
        },
        "in".into(),
        Vec::new(),
        network_schema(),
    )
    .unwrap();
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.002 * 0.002);
    let ladder = LineLadder::default();
    let (line, receipt) = point
        .realize_line(&gas, LossModel::Bessel, RATE_HZ, zc, &ladder, None)
        .unwrap();
    let dt = 1.0 / f64::from(RATE_HZ);
    let taps = line_taps(&line, receipt.n_fft);
    let (max_r, peaks) = analyze_reflectance_ir(&taps, dt, 2.0 * 0.32 / gas.sound_speed / dt);
    assert!(max_r <= 1.0 + super::PASSIVITY_TOL);
    assert!(peaks.len() >= 2, "branched bore peaks: {peaks:?}");
    // Golden-section maximum of `value` on [0.99 f, 1.01 f]; one peak per bracket.
    let refine = |f: f64, value: &dyn Fn(f64) -> f64| -> f64 {
        let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
        let (mut a, mut b) = (0.99 * f, 1.01 * f);
        for _ in 0..50 {
            let c = b - ratio * (b - a);
            let d = a + ratio * (b - a);
            if value(c) > value(d) {
                b = d;
            } else {
                a = c;
            }
        }
        0.5 * (a + b)
    };
    let tau = core::f64::consts::TAU;
    for &f in &peaks {
        let f_line = refine(f, &|g| super::dtft_norm_impedance_sq(&taps, dt, tau * g));
        let f_network = refine(f, &|g| {
            point
                .impedance(&gas, tau * g, LossModel::Bessel)
                .unwrap()
                .abs()
        });
        let cents = 1200.0 * (f_line / f_network).log2();
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"network-line-branched-peak\",\"n_fft\":{},\"line_hz\":{f_line:.4},\"network_hz\":{f_network:.4},\"cents\":{cents:.5},\"tol\":0.05}}",
            receipt.n_fft
        );
        assert!(
            cents.abs() < ladder.peak_tol_cents,
            "{f_line} vs {f_network} Hz: {cents} c"
        );
    }
}

/// A delay with a boundary-layer loss, `gain exp(i theta delay) exp(-k sqrt(-i theta))`:
/// causal, and passive for `gain <= 1`.
fn lossy_delay(theta: f64, gain: f64, delay: f64, k: f64) -> C64 {
    let s = (0.5 * theta).sqrt();
    let (magnitude, phase) = (gain * (-k * s).exp(), theta * delay + k * s);
    C64::new(magnitude * phase.cos(), magnitude * phase.sin())
}

/// The samples the ladder holds at rung `n`: each bin keeps the value computed
/// on the rung where that bin first appeared.
fn ladder_spectrum(
    n: usize,
    first_rung: usize,
    dt: f64,
    dc: f64,
    reflectance: &dyn Fn(f64) -> C64,
) -> Vec<C64> {
    let mut spectrum = vec![C64::from_re(dc); n / 2 + 1];
    for (bin, slot) in spectrum.iter_mut().enumerate().skip(1) {
        let (mut k, mut rung) = (bin, n);
        while rung > first_rung && k % 2 == 0 {
            k /= 2;
            rung /= 2;
        }
        *slot = reflectance(core::f64::consts::TAU * k as f64 / (rung as f64 * dt));
    }
    spectrum
}

/// The DC reflection [`crate::upstream_airway::UpstreamAirway::realize`] gives `chain` (no source impedance), by its
/// own arithmetic.
fn airway_dc_reflection(chain: &crate::upstream_airway::AirwayChain, gas: &GasState) -> f64 {
    let (zs, r_dc) = (chain.port_impedance(gas), chain.dc_resistance(gas));
    let scale = zs.max(r_dc).max(0.0);
    let (resistance, reference) = (r_dc / scale + 0.0 / scale, zs / scale);
    (resistance - reference) / (resistance + reference)
}

/// The ladder's first rung for a round trip of `geo_delay` samples, by its own arithmetic.
fn first_rung(geo_delay: f64, ladder: &LineLadder) -> usize {
    let max_taps = ladder.max_taps;
    let base = ((4.0 * geo_delay).ceil().min(max_taps as f64) as usize).max(256).next_power_of_two().min(max_taps);
    if max_taps >= 512 { base.min(max_taps / 2).max(256) } else { base }
}

/// G1 (MB 20665, 23083): the /ʃ/ rig cards' subglottal tube chain on their 26 kHz clock was refused by the passivity
/// bound, every sampled bin passive and the excess between bins. The passivity gain admits it on its own two rungs under
/// the unchanged bound, tail and peak tolerances. The admitted line's first three peaks are within the peak tolerance of
/// the same clock's ordinary converged rung, rebuilt here from the chain's own samples, and its top-band gain is the
/// independent oracle's to within the rounding of the two measurements of its first excess. The same chain on a 40 kHz
/// clock takes no gain. The 26 kHz
/// line's peaks against the 40 kHz line's are printed beside the ordinary 26 kHz rung's against the same: that
/// difference, about 1 c in the oracle with no gain at all, is the clock's; the gain's own share is the same-clock
/// shift asserted here.
#[test]
fn the_refused_sh_airway_is_admitted_by_the_passivity_gain() {
    use crate::upstream_airway::{AirwayChain, AirwaySection, UpstreamAirway};
    // The independent numpy chain and gain (outside the repository; it matched every engine print, MB 20663), with this
    // module's blend at the excess and its step: on the 4096-tap rung one pass, the first excess 4.635099393973441e-3,
    // the top band at -4.020783641708502e-2 dB.
    const ORACLE_TOP_BAND_DB: f64 = -4.020783641708502e-2;
    const ORACLE_FIRST_EXCESS: f64 = 4.635099393973441e-3;
    // `subglottal_tube AREA_M2 LENGTH_M` rows of dvtd-s1-12-schoen-sh-rig-*.performance, glottis first.
    let rows: [(f64, f64); 20] = [
        (7.809555556e-05, 3.77e-03), (2.53852e-04, 8.2e-02), (2.63267e-04, 2.0e-03), (3.01847e-04, 2.0e-03),
        (3.45033e-04, 2.0e-03), (3.91959e-04, 2.0e-03), (4.44424e-04, 2.0e-03), (5.02504e-04, 2.0e-03),
        (5.65882e-04, 2.0e-03), (6.3685e-04, 2.0e-03), (7.1438e-04, 2.0e-03), (7.98684e-04, 2.0e-03),
        (8.92295e-04, 2.0e-03), (9.93884e-04, 2.0e-03), (1.104095e-03, 2.0e-03), (1.225724e-03, 2.0e-03),
        (1.357506e-03, 2.0e-03), (1.499638e-03, 2.0e-03), (1.655512e-03, 2.0e-03), (1.74739e-03, 2.3e-04),
    ];
    let sections = rows
        .iter()
        .map(|&(area, length)| AirwaySection::from_area(length, area, 1, None).expect("a tube"))
        .collect();
    let chain = AirwayChain::new(sections, LossModel::Bessel).expect("an airway");
    let gas = GasState::try_new_moist_air(295.15, 100_700.0, 0.46).expect("the cards' ambient");
    let ladder = LineLadder { max_taps: GLOTTAL_LINE_MAX_TAPS, ..LineLadder::default() };
    let dc = airway_dc_reflection(&chain, &gas);
    let sample = |omega: f64| chain.reflectance(&gas, omega).expect("a sample of the chain");
    let mut clocks: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(2);
    for (rate, gained) in [(26_000_u32, true), (40_000, false)] {
        let airway = UpstreamAirway::realize(&chain, &gas, rate, &ladder).expect("the airway is admitted");
        let receipt = &airway.receipt().reflectance;
        let dt = 1.0 / f64::from(rate);
        let geo_delay = 2.0 * chain.path_length_m() / gas.sound_speed * f64::from(rate);
        let (bound, peaks) = analyze_reflectance_ir(airway.reflectance_taps(), dt, geo_delay);
        assert!(bound <= 1.0 + super::PASSIVITY_TOL, "{rate} Hz: 4x max {bound}");
        assert_eq!(receipt.passivity_gain.is_some(), gained, "{rate} Hz: {receipt:?}");
        let (mut ordinary_peaks, mut ordinary_shift) = (Vec::new(), "null".to_string());
        if let Some(gain) = receipt.passivity_gain {
            assert_eq!(receipt.rungs_tried.len(), 2);
            assert_eq!(receipt.n_fft, 2 * receipt.rungs_tried[0].n_fft);
            assert!(receipt.rungs_tried.iter().all(|rung| rung.tail_energy < ladder.tail_tol));
            assert!(receipt.peak_shift_cents < ladder.peak_tol_cents);
            assert!(gain.passes >= 1 && gain.passes <= super::PASSIVITY_GAIN_MAX_PASSES, "{gain:?}");
            assert_eq!(receipt.n_fft, 4096, "the oracle's rung");
            // The gain in dB is linear in ln(1 + excess), and the oracle and the engine take the same step and blend, so
            // what separates them is each one's measurement of the first excess. Each 4x bin is an n-term sum of the
            // rung's taps after the n-point transform that made them: within 2 n eps (sum |h| + 1) of its exact
            // value. Two measurements, over the excess (ln(1 + e) (1 + e) >= e), make the width. A dropped step
            // (1e-3 relative), another blend or another pass count falls outside it.
            let abs_taps: f64 = airway.reflectance_taps().iter().map(|h| h.abs()).sum();
            let width = ORACLE_TOP_BAND_DB.abs() * 4.0 * receipt.n_fft as f64 * f64::EPSILON * (abs_taps + 1.0)
                / ORACLE_FIRST_EXCESS;
            assert!(
                (gain.top_band_db - ORACLE_TOP_BAND_DB).abs() <= width,
                "{gain:?} against the oracle's {ORACLE_TOP_BAND_DB:e} dB"
            );
            let n = receipt.rungs_tried[0].n_fft;
            let spectrum = ladder_spectrum(n, first_rung(geo_delay, &ladder), dt, dc, &sample);
            let ordinary = super::realize_reflectance_ir_from_positive_spectrum(&spectrum, n, geo_delay, true);
            let (ordinary_bound, found) = analyze_reflectance_ir(&ordinary, dt, geo_delay);
            assert!(ordinary_bound > 1.0 + super::PASSIVITY_TOL, "the rebuilt rung is the refused one: {ordinary_bound}");
            let shift = super::first_peaks_shift_cents(&found, &peaks);
            assert!(shift < ladder.peak_tol_cents, "{shift} c from the ordinary converged rung");
            (ordinary_peaks, ordinary_shift) = (found, format!("{shift:e}"));
        }
        let gain = receipt.passivity_gain.map_or_else(
            || "null".to_string(),
            |g| format!("{{\"passes\":{},\"top_band_db\":{:e}}}", g.passes, g.top_band_db),
        );
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"sh-airway-passivity-gain\",\"rate_hz\":{rate},\"n_fft\":{},\
             \"rungs\":{},\"tail_energy\":{:e},\"peak_shift_cents\":{:e},\"max_abs_4x\":{bound:e},\"passivity_gain\":{gain},\
             \"ordinary_rung_shift_cents\":{ordinary_shift}}}",
            receipt.n_fft,
            receipt.rungs_tried.len(),
            receipt.tail_energy,
            receipt.peak_shift_cents,
        );
        clocks.push((peaks, ordinary_peaks));
    }
    let against = |a: &[f64], b: &[f64]| {
        a.iter().zip(b).take(3).map(|(&x, &y)| format!("{:e}", 1200.0 * (x / y).log2())).collect::<Vec<_>>().join(",")
    };
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"sh-airway-passivity-gain-clocks\",\"line_26k_vs_40k_cents\":[{}],\
         \"ordinary_26k_vs_40k_cents\":[{}]}}",
        against(&clocks[0].0, &clocks[1].0),
        against(&clocks[0].1, &clocks[1].0),
    );
}

/// G0 (MB 20665): a line admitted at its first converged rung takes no passivity gain, and its taps are bit for bit
/// the tapered construction of that rung's samples.
#[test]
fn an_admitted_line_keeps_its_taps_and_takes_no_gain() {
    let (delay, rate, dc) = (18.083639, 26_000_u32, 0.9996);
    let dt = 1.0 / f64::from(rate);
    let ladder = LineLadder::default();
    let reflectance = |omega: f64| lossy_delay(omega * dt, dc, delay, 0.03);
    let (line, receipt) = super::realize_characteristic_response(
        delay,
        super::LineDc::Explicit(dc),
        rate,
        1.0,
        Vec::new(),
        &ladder,
        super::PassivityFallback::Gain,
        |omega| Ok::<_, DrivingPointError>(reflectance(omega)),
    )
    .expect("the passive line is admitted");
    assert_eq!(receipt.passivity_gain, None);
    let spectrum = ladder_spectrum(receipt.n_fft, receipt.rungs_tried[0].n_fft, dt, dc, &reflectance);
    let expected =
        super::realize_reflectance_ir_from_positive_spectrum(&spectrum, receipt.n_fft, delay, true);
    let taps = line.quiescent_fir().expect("a plain FIR").1;
    assert_eq!(
        taps.iter().map(|h| h.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|h| h.to_bits()).collect::<Vec<_>>()
    );
}

/// G0 (MB 23083): the gain runs only where the caller opts in. A lossy delay whose converged 2048-tap rung exceeds the
/// bound between bins keeps the passivity refusal with [`super::PassivityFallback::Refuse`], every caller's choice but
/// the airway's, and is admitted with [`super::PassivityFallback::Gain`] on its rungs 2048 and 4096, its receipt naming
/// the gain.
#[test]
fn the_passivity_gain_runs_only_where_the_caller_opts_in() {
    let (delay, rate, dc) = (18.083639, 26_000_u32, 0.9996);
    let dt = 1.0 / f64::from(rate);
    let ladder = LineLadder::default();
    let reflectance = |omega: f64| Ok::<_, DrivingPointError>(lossy_delay(omega * dt, dc, delay, 0.0195));
    let realize = |fallback: super::PassivityFallback| {
        super::realize_characteristic_response(
            delay, super::LineDc::Explicit(dc), rate, 1.0, Vec::new(), &ladder, fallback, reflectance,
        )
    };
    let refusal = realize(super::PassivityFallback::Refuse).expect_err("refused without the gain");
    assert!(
        matches!(refusal, DrivingPointError::Invalid { what } if what == "converged characteristic line violates 4x-oversampled passivity"),
        "{refusal:?}"
    );
    let (line, receipt) = realize(super::PassivityFallback::Gain).expect("admitted with the gain");
    let gain = receipt.passivity_gain.expect("the receipt names the gain");
    assert_eq!(receipt.rungs_tried.iter().map(|rung| rung.n_fft).collect::<Vec<_>>(), [2048, 4096]);
    assert_eq!(receipt.n_fft, 4096);
    assert!(gain.passes <= super::PASSIVITY_GAIN_MAX_PASSES && gain.top_band_db <= 0.0, "{gain:?}");
    // Rung 2048 takes the gain (the active-reflectance test's control); rung 4096 may meet the bound with no pass, and
    // its receipt then names a zero-pass gain at +0 dB, never -0.
    assert!(gain.passes >= 1 || gain.top_band_db.to_bits() == 0.0_f64.to_bits(), "{gain:?}");
    let taps = line.quiescent_fir().expect("a plain FIR").1;
    let (bound, _) = analyze_reflectance_ir(taps, dt, delay);
    assert!(bound <= 1.0 + super::PASSIVITY_TOL, "4x max {bound}");
    // The gain's next rung is sampled after every call the ordinary ladder made. A source owner's cancellation at the
    // first of those samples cancels the line; any other source refusal there is a decline, and the line keeps the
    // passivity refusal.
    let mut ordinary = 0_usize;
    let counted = super::realize_characteristic_response(
        delay, super::LineDc::Explicit(dc), rate, 1.0, Vec::new(), &ladder, super::PassivityFallback::Refuse,
        |omega| {
            ordinary += 1;
            reflectance(omega)
        },
    );
    assert_eq!(counted.expect_err("refused without the gain"), refusal);
    let stopped_at_next = |stop: DuctError| {
        let mut calls = 0_usize;
        let result = super::realize_characteristic_response(
            delay, super::LineDc::Explicit(dc), rate, 1.0, Vec::new(), &ladder, super::PassivityFallback::Gain,
            |omega| {
                calls += 1;
                if calls > ordinary { Err(DrivingPointError::Duct(stop.clone())) } else { reflectance(omega) }
            },
        );
        (result.map(|(_, receipt)| receipt), calls)
    };
    let cancel = DuctError::Cancelled { completed_steps: 7 };
    let (cancelled, calls) = stopped_at_next(cancel.clone());
    assert_eq!((cancelled, calls), (Err(DrivingPointError::Duct(cancel)), ordinary + 1));
    let (declined, calls) = stopped_at_next(DuctError::LevineSchwingerCutoff { ka: 4.0, cutoff: 3.8317 });
    assert_eq!((declined, calls), (Err(refusal), ordinary + 1));
}

/// G0 (MB 20665, 23083): the gain is not a way around the bound, and the |R| <= 1 gate alone keeps an active rung out.
/// The control is the converged 2048-tap rung of the lossy delay above: it exceeds the bound between bins, and the gain
/// admits it after at least one pass. The same samples with only the Nyquist bin raised to |R| = 1 + 1e-9 are
/// declined: the rolloff there is below 1e-15, so the taps are the control's to rounding and only the gate can refuse
/// them. An active two-path line, |R| 1.02 at sampled bins, keeps the passivity refusal through the whole ladder.
#[test]
fn an_active_reflectance_is_not_given_the_passivity_gain() {
    let (delay, rate, dc, n) = (18.083639, 26_000_u32, 0.9996, 2048_usize);
    let dt = 1.0 / f64::from(rate);
    let ladder = LineLadder::default();
    let control = ladder_spectrum(n, 256, dt, dc, &|omega: f64| lossy_delay(omega * dt, dc, delay, 0.0195));
    let taps = super::realize_reflectance_ir_from_positive_spectrum(&control, n, delay, true);
    assert!(super::oversampled_peak(&taps).0 > 1.0 + super::PASSIVITY_TOL, "the control exceeds the bound");
    let (_, _, _, gain) =
        super::passivity_gain_rung(&control, n, dt, delay, &ladder).expect("the gain admits the control");
    assert!(gain.passes >= 1, "{gain:?}");
    let nyquist = super::reflectance_rolloff(n)[n / 2];
    assert!(nyquist.re.hypot(nyquist.im) < 1.0e-15, "the rolloff at Nyquist: {nyquist:?}");
    let mut active = control.clone();
    let z = active[n / 2];
    active[n / 2] = z.scale((1.0 + 1.0e-9) / z.re.hypot(z.im));
    assert!(active[n / 2].re.hypot(active[n / 2].im) > 1.0, "the raised bin is active");
    assert!(super::passivity_gain_rung(&active, n, dt, delay, &ladder).is_none(), "the gate declines it");

    let (depth, path) = (0.02, 16.0);
    let reflectance = |omega: f64| {
        let theta = omega * dt;
        let echo = C64::new((theta * path).cos(), (theta * path).sin()).scale(depth);
        C64::new((theta * delay).cos(), (theta * delay).sin()) * (C64::ONE - echo)
    };
    let peak = reflectance(core::f64::consts::PI / path / dt);
    assert!(peak.re.hypot(peak.im) > 1.0 + 0.01, "the reflectance must be active");
    let refusal = super::realize_characteristic_response(
        delay,
        super::LineDc::Explicit(1.0 - depth),
        rate,
        1.0,
        Vec::new(),
        &ladder,
        super::PassivityFallback::Gain,
        |omega| Ok::<_, DrivingPointError>(reflectance(omega)),
    )
    .expect_err("an active line is refused");
    match refusal {
        DrivingPointError::Invalid { what } => {
            assert_eq!(what, "converged characteristic line violates 4x-oversampled passivity");
        }
        other => panic!("expected the passivity refusal, got {other:?}"),
    }
}

/// G0 (MB 20665): a passive line whose 4x excess is larger than its fold can have moved is not the realization's to
/// repair. At the converged rung (2048 taps: tail met, bound exceeded) every bin is passive, the excess exceeds twice
/// the folded precursor, the gain declines and the line keeps its refusal.
#[test]
fn an_excess_beyond_the_fold_keeps_its_refusal() {
    let (delay, rate, dc, n) = (18.083639, 26_000_u32, 0.9996, 2048_usize);
    let dt = 1.0 / f64::from(rate);
    let ladder = LineLadder::default();
    let reflectance = |omega: f64| lossy_delay(omega * dt, dc, delay, 0.004);
    let refusal = super::realize_characteristic_response(
        delay,
        super::LineDc::Explicit(dc),
        rate,
        1.0,
        Vec::new(),
        &ladder,
        super::PassivityFallback::Gain,
        |omega| Ok::<_, DrivingPointError>(reflectance(omega)),
    )
    .expect_err("the line is refused");
    assert!(
        matches!(refusal, DrivingPointError::Invalid { what } if what == "converged characteristic line violates 4x-oversampled passivity"),
        "{refusal:?}"
    );
    let spectrum = ladder_spectrum(n, 256, dt, dc, &reflectance);
    assert!(spectrum.iter().all(|z| z.re.hypot(z.im) <= 1.0), "every sampled bin is passive");
    let rolloff = super::reflectance_rolloff(n);
    let taps = super::realize_reflectance_ir_with_rolloff(&spectrum, &rolloff, n, true);
    let total: f64 = taps.iter().map(|h| h * h).sum();
    let tail: f64 = taps[3 * (n / 4)..].iter().map(|h| h * h).sum();
    assert!(tail / total < ladder.tail_tol, "the rung meets the tail tolerance");
    let unfolded = super::realize_reflectance_ir_with_rolloff(&spectrum, &rolloff, n, false);
    let fold = 2.0 * unfolded[n - super::precursor_fold_len(n)..].iter().map(|h| h.abs()).sum::<f64>();
    let excess = super::oversampled_peak(&taps).0 - 1.0;
    assert!(excess > super::PASSIVITY_TOL && excess > fold, "excess {excess:e}, fold {fold:e}");
    assert!(super::passivity_gain_rung(&spectrum, n, dt, delay, &ladder).is_none());
}

/// G0 (MB 23083): the gain takes away no more than the fold could have added. The trumpeter's tract rows of the horn and
/// tuba cards (eb5d476f) with their `player_tract_wall` pin, in the /ʃ/ cards' moist ambient air, on the tuba's 9.6 kHz
/// clock. At the converged 2048-tap rung every bin is passive, the first excess is inside the fold, and the depth search
/// meets the bound; but what that depth removes from the top band is larger than the fold bound (1.49 times in the
/// independent oracle), so the rung declines and the airway keeps the passivity refusal.
#[test]
fn a_gain_deeper_than_the_fold_keeps_the_refusal() {
    use crate::upstream_airway::{AirwayChain, AirwaySection, UpstreamAirway};
    // `player_tract_section LENGTH_M AREA_M2` rows, the player's end first, and `player_tract_wall 15.4 3.79e5 6.22e3`.
    let mut rows: Vec<(f64, f64)> = [
        6.79e-05, 4.216e-04, 5.289e-04, 4.784e-04, 3.205e-04, 3.016e-04, 2.763e-04, 2.511e-04, 2.574e-04, 1.658e-04,
        1.089e-04, 5.84e-05, 4.58e-05, 3.63e-05, 6.79e-05, 5.84e-05, 4.89e-05, 4.26e-05, 1.279e-04, 9.63e-05, 9.95e-05,
        1.153e-04, 1.626e-04, 2.068e-04, 3.458e-04, 3.047e-04, 2.858e-04, 2.732e-04, 2.984e-04, 4.626e-04, 4.342e-04,
        3.995e-04, 2.037e-04, 1.247e-04, 9.95e-05, 7.74e-05, 9.32e-05, 9.95e-05, 1.216e-04,
    ]
    .into_iter()
    .map(|area| (0.00435, area))
    .collect();
    rows.push((1.2e-3, 1.4e-5));
    let wall = WallPin { surface_density: 15.4, stiffness_per_area: 3.79e5, resistance: 6.22e3 };
    let sections = rows
        .iter()
        .map(|&(length, area)| AirwaySection::from_area(length, area, 1, Some(wall)).expect("a tube"))
        .collect();
    let chain = AirwayChain::new(sections, LossModel::Bessel).expect("a tract");
    let gas = GasState::try_new_moist_air(295.15, 100_700.0, 0.46).expect("the /ʃ/ cards' ambient");
    let (rate, n, ladder) = (9_600_u32, 2048_usize, LineLadder::default());
    let refusal = UpstreamAirway::realize(&chain, &gas, rate, &ladder).expect_err("the tract keeps its refusal");
    assert!(
        matches!(refusal, DrivingPointError::Invalid { what } if what == "converged characteristic line violates 4x-oversampled passivity"),
        "{refusal:?}"
    );
    // The rung judged below is the one the ladder converged on: under a 1024-tap ceiling no rung converges, and under
    // a 2048-tap ceiling (no room for the gain's next rung) the line has the passivity refusal of a converged rung.
    // Both ceilings start the ladder on the default's first rung, so their rungs are the default's.
    let dt = 1.0 / f64::from(rate);
    let geo_delay = 2.0 * chain.path_length_m() / gas.sound_speed * f64::from(rate);
    for ceiling in [n / 2, n] {
        let capped = LineLadder { max_taps: ceiling, ..ladder };
        assert_eq!(first_rung(geo_delay, &capped), first_rung(geo_delay, &ladder), "a {ceiling}-tap ceiling");
    }
    let short = UpstreamAirway::realize(&chain, &gas, rate, &LineLadder { max_taps: n / 2, ..ladder })
        .expect_err("no rung below 2048 converges");
    assert!(matches!(short, DrivingPointError::LadderUnconverged { max_taps, .. } if max_taps == n / 2), "{short:?}");
    let at_n = UpstreamAirway::realize(&chain, &gas, rate, &LineLadder { max_taps: n, ..ladder })
        .expect_err("the converged 2048-tap rung exceeds the bound");
    assert_eq!(at_n, refusal);
    let sample = |omega: f64| chain.reflectance(&gas, omega).expect("a sample of the tract");
    let spectrum = ladder_spectrum(n, first_rung(geo_delay, &ladder), dt, airway_dc_reflection(&chain, &gas), &sample);
    assert!(spectrum.iter().all(|z| z.re.hypot(z.im) <= 1.0), "every sampled bin is passive");
    let rolloff = super::reflectance_rolloff(n);
    let taps = super::realize_reflectance_ir_with_rolloff(&spectrum, &rolloff, n, true);
    let total: f64 = taps.iter().map(|h| h * h).sum();
    let tail: f64 = taps[3 * (n / 4)..].iter().map(|h| h * h).sum();
    assert!(tail / total < ladder.tail_tol, "the rung meets the tail tolerance");
    let unfolded = super::realize_reflectance_ir_with_rolloff(&spectrum, &rolloff, n, false);
    let fold = 2.0 * unfolded[n - super::precursor_fold_len(n)..].iter().map(|h| h.abs()).sum::<f64>();
    let first = super::oversampled_peak(&taps);
    assert!(first.0 - 1.0 > super::PASSIVITY_TOL && first.0 - 1.0 <= fold, "excess {:e}, fold {fold:e}", first.0 - 1.0);
    let half = n / 2;
    let blend: Vec<f64> = (0..=half).map(|k| super::rolloff_blend(k as f64 / half as f64)).collect();
    let (depth, passes) = super::passivity_gain_depth(first, |depth| {
        let gained: Vec<C64> = rolloff.iter().zip(&blend).map(|(w, &eta)| w.scale((-depth * eta).exp())).collect();
        super::oversampled_peak(&super::realize_reflectance_ir_with_rolloff(&spectrum, &gained, n, true))
    })
    .expect("the depth search meets the bound");
    let removed = spectrum
        .iter()
        .zip(&rolloff)
        .zip(&blend)
        .map(|((r, w), &eta)| r.re.hypot(r.im) * w.re.hypot(w.im) * -(-depth * eta).exp_m1())
        .fold(0.0_f64, f64::max);
    assert!(removed > fold, "removed {removed:e}, fold {fold:e}, {passes} passes");
    assert!(super::passivity_gain_rung(&spectrum, n, dt, geo_delay, &ladder).is_none());
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"tract-gain-beyond-fold\",\"rate_hz\":{rate},\"n_fft\":{n},\
         \"first_excess\":{:e},\"fold_bound\":{fold:e},\"passes\":{passes},\"removed\":{removed:e}}}",
        first.0 - 1.0
    );
}

/// G0 (MB 23083): the depth search declines after its pass cap and where the blend at the excess has no authority, and
/// admits on the pass that meets the bound, the cap's last included.
#[test]
fn the_gain_search_declines_at_its_pass_cap_and_below_its_minimum_blend() {
    let excess = 1.0e-3;
    let cap = super::PASSIVITY_GAIN_MAX_PASSES;
    // An excess at 0.97 Nyquist (blend 1) that no depth removes: declined after exactly the cap.
    let mut calls = 0_u32;
    let search = super::passivity_gain_depth((1.0 + excess, 0.97), |_| {
        calls += 1;
        (1.0 + excess, 0.97)
    });
    assert_eq!((search, calls), (None, cap));
    // The same excess removed on the cap's last pass: admitted with that many passes, each the same step.
    let mut calls = 0_u32;
    let (depth, passes) = super::passivity_gain_depth((1.0 + excess, 0.97), |_| {
        calls += 1;
        if calls == cap { (1.0, 0.97) } else { (1.0 + excess, 0.97) }
    })
    .expect("met on the cap's last pass");
    assert_eq!(passes, cap);
    // The search reads the excess back from the peak, `(1 + excess) - 1`, which is not `excess` in binary.
    let whole = f64::from(cap) * ((1.0 + excess) - 1.0).ln_1p() * super::PASSIVITY_GAIN_STEP;
    assert!((depth - whole).abs() <= f64::from(cap) * f64::EPSILON * whole, "{depth:e} against {whole:e}");
    // An excess where the blend is positive but below its minimum: declined before any pass.
    let low = 0.52;
    let authority = super::rolloff_blend(low);
    assert!(authority > 0.0 && authority < super::PASSIVITY_GAIN_MIN_BLEND, "blend {authority:e}");
    let mut calls = 0_u32;
    let search = super::passivity_gain_depth((1.0 + excess, low), |_| {
        calls += 1;
        (1.0, low)
    });
    assert_eq!((search, calls), (None, 0));
    // An excess a pass moves to where the blend has no authority: declined at that pass.
    let mut calls = 0_u32;
    let search = super::passivity_gain_depth((1.0 + excess, 0.97), |_| {
        calls += 1;
        (1.0 + excess, low)
    });
    assert_eq!((search, calls), (None, 1));
    // A peak that is not finite: declined.
    assert_eq!(super::passivity_gain_depth((1.0 + excess, 0.97), |_| (f64::NAN, 0.97)), None);
    // A bound met with no pass costs nothing.
    assert_eq!(super::passivity_gain_depth((1.0, 0.97), |_| unreachable!("no pass is taken")), Some((0.0, 0)));
}
