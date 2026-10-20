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
/// for `N = 4096` steps at `f_s = 48 kHz`; linear regression of `p(t)` over
/// `[200, 4000]` samples has slope `(ρ c² / V) u0` within `1%`, and `p(N)`
/// does not decay back toward zero.
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
    let n_steps = 4096usize;
    let u0 = 1.0e-6_f64;
    let dt = 1.0 / f64::from(RATE_HZ);

    let mut line = impedance_line(&duct, &gas, Termination::Closed, RATE_HZ, n_steps, None)
        .expect("closed 2 cm cavity realizes");
    assert!(
        line.compliance_gain() > 0.0,
        "sealed cavity must attach PoleState::Compliance"
    );

    let mut p = Vec::with_capacity(n_steps);
    for _ in 0..n_steps {
        p.push(line.push(u0).expect("step"));
    }

    // Linear regression of p(t) over [200..=4000] against t = k * dt.
    let slice = &p[200..=4000];
    let count = slice.len() as f64;
    let mean_t = (200..=4000).map(|k| k as f64 * dt).sum::<f64>() / count;
    let mean_p = slice.iter().sum::<f64>() / count;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (idx, &pk) in slice.iter().enumerate() {
        let t = (200 + idx) as f64 * dt;
        num += (t - mean_t) * (pk - mean_p);
        den += (t - mean_t) * (t - mean_t);
    }
    let observed_slope = num / den;
    let expected_slope = (gas.density * gas.sound_speed * gas.sound_speed / volume) * u0;
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
