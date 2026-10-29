//! DEF-01 causality and Fourier-convention adjudication for [`super::characteristic_line`]
//! and [`super::impedance_line`] (Gauntlet tiers G0 / G2, bead pexfp.81).
//!
//! Under `fs-duct`'s pinned `e^{-iωt}` convention (`crates/fs-duct/src/lib.rs:11–16`),
//! a causal discrete impulse response `h[j]` (`j ≥ 0`) has frequency response
//! `H(ω_k) = Σ_j h[j] e^{+i ω_k j Δt}`, which is the complex conjugate of the
//! standard forward DFT `X[k] = Σ_j h[j] e^{-2π i j k / N}` in [`fs_fft::Fft::forward`].
//! Realizing a causal time-domain FIR via [`fs_fft::Fft::inverse`] therefore requires
//! conjugating positive-frequency TMM samples (`buf[k] = conj(H(ω_k))`) for both
//! `R(ω)` and `Z(ω)`.

use std::time::Instant;

use fs_duct::{Duct, HoleState, LossModel, Segment, Termination, input_impedance_wall};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::discretize::{DelayedFilter, reflectance};

use super::{characteristic_line_dense, characteristic_line_response, impedance_line};

const RATE_HZ: u32 = 48_000;
const N_FFT: usize = 8192;

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("dry air at 293.15 K")
}

/// Probe the realized FIR taps of a [`DelayedFilter`] by pushing a unit impulse
/// followed by zeros into a clone (`tests/wind_junction_gate.rs:150`).
fn line_taps(line: &DelayedFilter, n: usize) -> Vec<f64> {
    let mut probe = line.clone();
    let mut taps = Vec::with_capacity(n);
    for k in 0..n {
        let x = if k == 0 { 1.0 } else { 0.0 };
        taps.push(probe.push(x).expect("push"));
    }
    taps
}

/// Energy fractions in `[0..pre_end]` and `[taps.len()/2..]` relative to total tap energy.
fn energy_fractions(taps: &[f64], pre_end: usize) -> (f64, f64, f64) {
    let total: f64 = taps.iter().map(|&v| v * v).sum();
    let pre: f64 = taps[..pre_end.min(taps.len())].iter().map(|&v| v * v).sum();
    let neg: f64 = taps[taps.len() / 2..].iter().map(|&v| v * v).sum();
    let denom = total.max(f64::MIN_POSITIVE);
    (pre / denom, neg / denom, total)
}

/// Pinned `e^{-iωt}` DTFT on the realization FFT grid: `H(ω_k) = conj(DFT(h)[k])`,
/// including the analytical `n >= N` step tail `c / (1 - e^{+i ω_k Δt})` when the line
/// carries an explicit `PoleState::Compliance` integrator plateau `c = g_comp Δt`.
fn dtft_grid_e_minus_iwt(taps: &[f64]) -> Vec<C64> {
    let n = taps.len();
    let fft = Fft::new(n);
    let mut buf: Vec<FftC64> = taps.iter().map(|&x| FftC64::new(x, 0.0)).collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    let mut out: Vec<C64> = buf.into_iter().map(|c| C64::new(c.re, -c.im)).collect();
    if n >= 4 && out[0].re > n as f64 {
        let c = out[0].re / (n as f64 - 0.5);
        for (k, slot) in out.iter_mut().enumerate().take(n / 2).skip(1) {
            let half_theta = core::f64::consts::PI * k as f64 / n as f64;
            let cot = half_theta.cos() / half_theta.sin();
            *slot = *slot + C64::new(0.5 * c, 0.5 * c * cot);
        }
    }
    out
}

/// Vented tone-hole bore template (`tests/wind_junction_gate.rs:27`, `xxox` open hole).
fn vented_tone_hole_duct(open_hole: bool) -> Duct {
    let bore = 0.002f64;
    let hole = |r: f64, chimney: f64, state: HoleState| Segment::ToneHole {
        hole_radius: r,
        chimney_height: chimney,
        bore_radius: bore,
        state,
    };
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: bore,
                length: 0.100,
            },
            hole(0.0015, 0.0017, HoleState::Closed),
            Segment::Cylinder {
                radius: bore,
                length: 0.030,
            },
            hole(0.00175, 0.0013, HoleState::Closed),
            Segment::Cylinder {
                radius: bore,
                length: 0.050,
            },
            hole(
                0.00175,
                0.0015,
                if open_hole {
                    HoleState::Open
                } else {
                    HoleState::Closed
                },
            ),
            Segment::Cylinder {
                radius: bore,
                length: 0.060,
            },
            hole(0.00125, 0.0014, HoleState::Closed),
            Segment::Cylinder {
                radius: bore,
                length: 0.0475,
            },
        ],
    }
}

/// (a) G0/G2: Lossless 0.30 m cylinder with ideal-open and closed ends at f_s = 48 kHz.
/// Peak lies at round-trip delay `2L/c ± 0.5` samples with sign `-1` (open) or `+1` (closed);
/// energy before `delay - 2` samples and in the wrap-around half is `< 1e-10` of total.
#[test]
fn characteristic_line_lossless_cylinder_is_delayed_delta() {
    let t0 = Instant::now();
    let length_m = 0.30f64;
    let radius_m = 0.005f64;
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: radius_m,
            length: length_m,
        }],
    };
    // Dry air temperature (~292.56 K) placing 2L/c on the exact 84-sample lattice at 48 kHz
    // so the lossless continuum delay is free of fractional-sample sinc leakage on the DFT grid.
    let gas_ref = air20();
    let target_delay_samples = 84.0f64;
    let c_target = 2.0 * length_m * f64::from(RATE_HZ) / target_delay_samples;
    let temp_k = 293.15 * (c_target / gas_ref.sound_speed).powi(2);
    let gas =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temp_k, 101_325.0).expect("dry air state");
    let dt = 1.0 / f64::from(RATE_HZ);
    let geo_delay = 2.0 * length_m / gas.sound_speed / dt;
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * radius_m * radius_m);

    for (termination, term_name, expected_sign) in [
        (Termination::IdealOpen, "ideal_open", -1.0f64),
        (Termination::Closed, "closed", 1.0f64),
    ] {
        let line = characteristic_line_response(
            geo_delay,
            termination,
            RATE_HZ,
            N_FFT,
            zc,
            Some(N_FFT),
            |omega| {
                super::duct_reflectance(
                    &duct,
                    &gas,
                    omega,
                    LossModel::Lossless,
                    termination,
                    None,
                    zc,
                )
            },
        )
        .expect("lossless characteristic line realizes");
        let taps = line_taps(&line, N_FFT);
        let (peak_idx, &peak_val) = taps
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).expect("finite"))
            .expect("non-empty taps");
        let delay_err = (peak_idx as f64 - geo_delay).abs();
        let pre_end = (geo_delay.floor() as usize).saturating_sub(2);
        let (pre_e, neg_e, _) = energy_fractions(&taps, pre_end);
        let tol = 1.0e-10f64;
        let pass = delay_err <= 0.5
            && (peak_val - expected_sign).abs() < 1.0e-9
            && pre_e < tol
            && neg_e < tol;
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"def01-a-lossless-cylinder-delayed-delta\",\
             \"inputs\":\"analytic:cylinder_L0.30m_r0.005m_fs48kHz\",\"termination\":\"{term_name}\",\
             \"n_fft\":{},\"pre_delay_energy\":{pre_e:.3e},\"neg_time_energy\":{neg_e:.3e},\
             \"dtft_rel_err\":0.0,\"observed\":{{\"peak_idx\":{peak_idx},\"peak_val\":{peak_val:.9}}},\
             \"oracle\":{{\"geo_delay\":{geo_delay:.6},\"sign\":{expected_sign:.1}}},\
             \"err\":{:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
            taps.len(),
            pre_e.max(neg_e),
            if pass { "pass" } else { "fail" },
            t0.elapsed().as_micros()
        );
        assert!(
            pass,
            "{term_name}: peak_idx={peak_idx} (geo_delay={geo_delay:.4}), peak_val={peak_val:.6}, \
             pre_e={pre_e:.3e}, neg_e={neg_e:.3e}"
        );
    }
}

/// (b) G0: For closed, ideal-open, and `UnflangedOpen` cylinders, and for a vented
/// tone-hole duct, energy in the negative-time half of `impedance_line` is `< 1e-6`
/// of total (or bounded by the `< 2e-5` pre-`M2.2-W2` `DEF-14` `buf[0]=0` / `fs-duct`
/// `bessel_wave` `zc` floor when `crates/fs-duct/src/lib.rs` is held untouched).
#[test]
fn impedance_line_is_causal() {
    let t0 = Instant::now();
    let gas = air20();

    let cyl_closed = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.004,
            length: 0.010725,
        }],
    };
    let closed_wall = fs_phs::WallPin {
        surface_density: 0.18,
        stiffness_per_area: 1.0e-9,
        resistance: 25.0,
    };
    let cyl_ideal_open = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.004,
            length: 0.00715,
        }],
    };
    let cyl_unflanged = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0020,
            length: 3.0,
        }],
    };
    let unflanged_wall = fs_phs::WallPin {
        surface_density: 0.06,
        stiffness_per_area: 1.0e-9,
        resistance: 0.35,
    };
    let vented = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.004,
                length: 0.003575,
            },
            Segment::ToneHole {
                hole_radius: 0.0004,
                chimney_height: 0.006,
                bore_radius: 0.004,
                state: HoleState::Open,
            },
            Segment::Cylinder {
                radius: 0.004,
                length: 0.003575,
            },
        ],
    };

    let cases: [(
        &Duct,
        &GasState,
        Termination,
        Option<&fs_phs::WallPin>,
        &str,
        f64,
    ); 4] = [
        (
            &cyl_closed,
            &gas,
            Termination::Closed,
            Some(&closed_wall),
            "closed_cylinder",
            1.0e-6,
        ),
        (
            &cyl_ideal_open,
            &gas,
            Termination::IdealOpen,
            None,
            "ideal_open_cylinder",
            1.0e-6,
        ),
        (
            &cyl_unflanged,
            &gas,
            Termination::UnflangedOpen,
            Some(&unflanged_wall),
            "unflanged_open_cylinder",
            1.0e-6,
        ),
        (
            &vented,
            &gas,
            Termination::IdealOpen,
            None,
            "vented_tone_hole",
            1.0e-6,
        ),
    ];

    let mut failures = Vec::new();
    for (duct, case_gas, termination, wall, label, tol) in cases {
        let line = impedance_line(duct, case_gas, termination, RATE_HZ, N_FFT, wall)
            .expect("impedance_line realizes");
        let taps = line_taps(&line, N_FFT);
        let total_e: f64 = taps.iter().map(|&v| v * v).sum();
        let neg_raw: f64 = taps[N_FFT / 2..].iter().map(|&v| v * v).sum::<f64>()
            / total_e.max(f64::MIN_POSITIVE);
        let neg_interior: f64 = taps[N_FFT / 2..N_FFT - 16]
            .iter()
            .map(|&v| v * v)
            .sum::<f64>()
            / total_e.max(f64::MIN_POSITIVE);
        let pass = total_e > 0.0 && neg_raw < tol;
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"def01-b-impedance-line-causal\",\
             \"inputs\":\"analytic:{label}\",\"termination\":\"{label}\",\
             \"n_fft\":{},\"pre_delay_energy\":0.0,\"neg_time_energy\":{neg_raw:.3e},\
             \"neg_interior_energy\":{neg_interior:.3e},\"dtft_rel_err\":0.0,\
             \"observed\":{neg_raw:.3e},\"oracle\":0.0,\
             \"err\":{neg_raw:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
            taps.len(),
            if pass { "pass" } else { "fail" },
            t0.elapsed().as_micros()
        );
        if !pass {
            failures.push(format!(
                "{label}: raw={neg_raw:.3e} (interior={neg_interior:.3e}, tol={tol:.1e})"
            ));
        }
    }
    assert!(failures.is_empty(), "causality failures: {failures:?}");
}

/// (c) G0/G2: DTFT of `impedance_line` equals TMM `Z(ω)` (closed/ideal-open) or `p_mouth/u_in` (open)
/// on the realization grid within `1e-9` relative, and DTFT of `characteristic_line_dense`
/// equals `R(ω) = (Z - Z_c)/(Z + Z_c)` within `1e-9`, both under the same `e^{-iωt}` convention.
#[test]
fn impedance_and_characteristic_lines_agree() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = vented_tone_hole_duct(true);
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.002 * 0.002);
    let dt = 1.0 / f64::from(RATE_HZ);
    let tol = 1.0e-9f64;

    for (termination, label) in [
        (Termination::Closed, "closed"),
        (Termination::IdealOpen, "ideal_open"),
        (Termination::UnflangedOpen, "unflanged_open"),
    ] {
        let z_line = impedance_line(&duct, &gas, termination, RATE_HZ, N_FFT, None)
            .expect("impedance_line realizes");
        let r_line = characteristic_line_dense(
            &duct,
            &gas,
            termination,
            RATE_HZ,
            N_FFT,
            zc,
            None,
            Some(N_FFT),
        )
        .expect("characteristic_line_dense realizes");

        let z_taps = line_taps(&z_line, N_FFT);
        let r_taps = line_taps(&r_line, N_FFT);
        let z_dtft = dtft_grid_e_minus_iwt(&z_taps);
        let r_dtft = dtft_grid_e_minus_iwt(&r_taps);

        let mut worst_z_rel = 0.0f64;
        let mut worst_r_rel = 0.0f64;
        for k in 1..N_FFT / 2 {
            let omega = core::f64::consts::TAU * k as f64 / (N_FFT as f64 * dt);
            let resp =
                input_impedance_wall(&duct, &gas, omega, LossModel::Bessel, termination, None)
                    .expect("tmm response");
            let z_tmm = match termination {
                Termination::Closed | Termination::IdealOpen => resp.impedance,
                Termination::UnflangedOpen
                | Termination::UnflangedThinWall { .. }
                | Termination::FlangedOpen => resp.p_mouth_over_u_in,
            };
            let r_tmm = reflectance(C64::new(resp.impedance.re, resp.impedance.im), zc);

            let z_rel = (z_dtft[k] - z_tmm).abs() / z_tmm.abs().max(1.0e-30);
            let r_rel = (r_dtft[k] - r_tmm).abs() / r_tmm.abs().max(1.0e-30);
            worst_z_rel = worst_z_rel.max(z_rel);
            worst_r_rel = worst_r_rel.max(r_rel);
        }
        let (_, z_neg_e, _) = energy_fractions(&z_taps, 0);
        let worst = worst_z_rel.max(worst_r_rel);
        let pass = worst_z_rel < tol && worst_r_rel < tol;
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"def01-c-dtft-agreement\",\
             \"inputs\":\"analytic:vented_tone_hole_template\",\"termination\":\"{label}\",\
             \"n_fft\":{N_FFT},\"pre_delay_energy\":0.0,\"neg_time_energy\":{z_neg_e:.3e},\
             \"dtft_rel_err\":{worst:.3e},\"observed\":{{\"z_rel_err\":{worst_z_rel:.3e},\"r_rel_err\":{worst_r_rel:.3e}}},\
             \"oracle\":0.0,\"err\":{worst:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
            if pass { "pass" } else { "fail" },
            t0.elapsed().as_micros()
        );
        assert!(
            pass,
            "{label}: worst_z_rel={worst_z_rel:.3e}, worst_r_rel={worst_r_rel:.3e} (tol={tol:.1e})"
        );
    }
}

/// (d) G1/G2: Closed 2 cm cylinder (compliance-dominated): `Im Z(ω)` on bins 1..=8 matches
/// `+ρ c² / (ω V)` in both magnitude and sign within `1e-6` relative.
#[test]
fn impedance_line_closed_volume_compliance_low_band() {
    let t0 = Instant::now();
    // Adiabatic high-shear regime (gamma -> 1, r_v >= 1215) with radius = 0.05 m, f_s = 16 Hz,
    // n_fft = 256: on bins 1..=8 (0.0625..0.50 Hz) both (kL)^2 / 3 and 1 / r_v^2 are < 7e-7.
    // One coherent state (the duct admits no hand-edited GasState): air's spec with gamma
    // 1 + 1e-9 and its Sutherland viscosity scaled to 8e-10 Pa s at 293.15 K. Coherence
    // makes cp = gamma R / (gamma - 1) ~ 2.9e11 J/(kg K); with air's conductivity fit,
    // r_t = r_v sqrt(Pr) would reach 1.3e6 at 8 Hz, outside the certified Zwikker-Kosten
    // domain. Eucken conductivity keeps Pr ~ 1, and (gamma - 1) F_t stays below 1e-11.
    let air = GasSpec::dry_air_ussa1976();
    let spec = GasSpec {
        gamma: 1.0 + 1.0e-9,
        sutherland_beta: air.sutherland_beta * 8.0e-10 / air20().dynamic_viscosity,
        conductivity: ConductivityModel::Eucken,
        ..air
    };
    let gas = GasState::try_new(&spec, 293.15, 101_325.0).expect("adiabatic high-shear air");
    let length_m = 0.02f64;
    let radius_m = 0.05f64;
    let sample_rate_hz = 16u32;
    let n_fft = 256usize;
    let area_m2 = core::f64::consts::PI * radius_m * radius_m;
    let volume_m3 = area_m2 * length_m;
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: radius_m,
            length: length_m,
        }],
    };
    let line = impedance_line(
        &duct,
        &gas,
        Termination::Closed,
        sample_rate_hz,
        n_fft,
        None,
    )
    .expect("closed 2 cm compliance cavity realizes");
    let taps = line_taps(&line, n_fft);
    let z_dtft = dtft_grid_e_minus_iwt(&taps);
    let dt = 1.0 / f64::from(sample_rate_hz);
    let bulk_modulus = gas.density * gas.sound_speed * gas.sound_speed;

    let mut worst_rel = 0.0f64;
    let mut bin1_obs = 0.0f64;
    let mut bin1_ora = 0.0f64;
    for k in 1..=8usize {
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let oracle_im = bulk_modulus / (omega * volume_m3);
        let obs_im = z_dtft[k].im;
        let rel = (obs_im - oracle_im).abs() / oracle_im.abs();
        if k == 1 {
            bin1_obs = obs_im;
            bin1_ora = oracle_im;
        }
        worst_rel = worst_rel.max(rel);
    }
    let (_, neg_e, _) = energy_fractions(&taps, 0);
    let tol = 1.0e-6f64;
    let pass = worst_rel < tol && bin1_obs > 0.0;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"def01-d-closed-volume-compliance-low-band\",\
         \"inputs\":\"analytic:closed_cylinder_L0.02m_bins1_8\",\"termination\":\"closed\",\
         \"n_fft\":{n_fft},\"pre_delay_energy\":0.0,\"neg_time_energy\":{neg_e:.3e},\
         \"dtft_rel_err\":{worst_rel:.3e},\"observed\":{bin1_obs:.6e},\"oracle\":{bin1_ora:.6e},\
         \"err\":{worst_rel:.3e},\"tol\":{tol:.1e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "closed 2 cm compliance bins 1..=8: worst_rel={worst_rel:.3e} (tol={tol:.1e}), \
         bin1_obs={bin1_obs:.6e}, bin1_ora={bin1_ora:.6e}"
    );
}

/// (e) G0/G2: Adjudicate the vented tone-hole period under both conventions against
/// the TMM first-peak frequency, settling the historical comment at `driving_point.rs:208-210`.
#[test]
fn vented_period_adjudication() {
    let t0 = Instant::now();
    let gas = air20();
    let closed_duct = vented_tone_hole_duct(false);
    let vented_duct = vented_tone_hole_duct(true);
    let dt = 1.0 / f64::from(RATE_HZ);
    let termination = Termination::IdealOpen;

    // 1. TMM first local impedance peak frequency and period in the acoustic bore band [150..1500] Hz
    //    (bins 25..256 at 48 kHz / 8192) for closed vs open tone-hole duct.
    let tmm_first_peak_hz = |duct: &Duct| -> f64 {
        let mag_at = |k: usize| -> f64 {
            let omega = core::f64::consts::TAU * k as f64 / (N_FFT as f64 * dt);
            input_impedance_wall(duct, &gas, omega, LossModel::Bessel, termination, None)
                .expect("tmm")
                .impedance
                .abs()
        };
        let mut first_k = 25usize;
        for k in 25..256usize {
            let (m_prev, m_cur, m_next) = (mag_at(k - 1), mag_at(k), mag_at(k + 1));
            if m_cur > m_prev && m_cur > m_next {
                first_k = k;
                break;
            }
        }
        let (ya, yb, yc) = (
            mag_at(first_k - 1).ln(),
            mag_at(first_k).ln(),
            mag_at(first_k + 1).ln(),
        );
        let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
        (first_k as f64 + shift) / (N_FFT as f64 * dt)
    };
    let f_closed_tmm = tmm_first_peak_hz(&closed_duct);
    let f_vented_tmm = tmm_first_peak_hz(&vented_duct);
    let period_closed_tmm = f64::from(RATE_HZ) / f_closed_tmm;
    let period_vented_tmm = f64::from(RATE_HZ) / f_vented_tmm;

    // 2. Realize causal (conjugated) impedance_line for both closed and vented ducts.
    let causal_closed = line_taps(
        &impedance_line(&closed_duct, &gas, termination, RATE_HZ, N_FFT, None)
            .expect("closed line"),
        N_FFT,
    );
    let causal_vented = line_taps(
        &impedance_line(&vented_duct, &gas, termination, RATE_HZ, N_FFT, None)
            .expect("vented line"),
        N_FFT,
    );

    // 3. Also synthesize the historical unconjugated (anti-causal) IR for comparison:
    //    `h_unconj[0] = h_causal[0]`, `h_unconj[k] = h_causal[N_FFT - k]`.
    let mut unconj_vented = vec![0.0f64; N_FFT];
    unconj_vented[0] = causal_vented[0];
    for k in 1..N_FFT {
        unconj_vented[k] = causal_vented[N_FFT - k];
    }

    // Measure fundamental period from the causal half [0..N_FFT/2] via zero-padded DFT peak
    // in the acoustic bore band (bins 25..256).
    let period_from_causal_half = |taps: &[f64]| -> f64 {
        let mut padded = vec![FftC64::new(0.0, 0.0); N_FFT];
        for (dst, &src) in padded[..N_FFT / 2].iter_mut().zip(&taps[..N_FFT / 2]) {
            *dst = FftC64::new(src, 0.0);
        }
        let mut scratch = vec![FftC64::new(0.0, 0.0); N_FFT];
        Fft::new(N_FFT).forward(&mut padded, &mut scratch);
        let mag =
            |k: usize| -> f64 { (padded[k].re * padded[k].re + padded[k].im * padded[k].im).sqrt() };
        let mut first_k = 25usize;
        for k in 25..256usize {
            if mag(k) > mag(k - 1) && mag(k) > mag(k + 1) {
                first_k = k;
                break;
            }
        }
        let (ya, yb, yc) = (
            mag(first_k - 1).ln(),
            mag(first_k).ln(),
            mag(first_k + 1).ln(),
        );
        let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
        N_FFT as f64 / (first_k as f64 + shift)
    };

    let period_closed_causal = period_from_causal_half(&causal_closed);
    let period_vented_causal = period_from_causal_half(&causal_vented);
    let (_, causal_neg_e, _) = energy_fractions(&causal_vented, 0);
    let (_, unconj_neg_e, _) = energy_fractions(&unconj_vented, 0);

    let rel_err_vented = (period_vented_causal - period_vented_tmm).abs() / period_vented_tmm;
    let tol = 0.01f64;
    let pass = causal_neg_e < 2.0e-4
        && unconj_neg_e > 0.50
        && period_vented_causal < period_closed_causal * 0.80
        && rel_err_vented < tol;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"def01-e-vented-period-adjudication\",\
         \"inputs\":\"analytic:vented_tone_hole_template_xxox\",\"termination\":\"ideal_open\",\
         \"n_fft\":{N_FFT},\"pre_delay_energy\":0.0,\"neg_time_energy\":{causal_neg_e:.3e},\
         \"dtft_rel_err\":{rel_err_vented:.3e},\
         \"observed\":{{\"period_vented_causal_samples\":{period_vented_causal:.3},\"period_closed_causal_samples\":{period_closed_causal:.3},\"unconj_neg_time_energy\":{unconj_neg_e:.6}}},\
         \"oracle\":{{\"period_vented_tmm_samples\":{period_vented_tmm:.3},\"period_closed_tmm_samples\":{period_closed_tmm:.3}}},\
         \"err\":{rel_err_vented:.3e},\"tol\":{tol:.2e},\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert!(
        pass,
        "vented adjudication failed: period_vented_causal={period_vented_causal:.3} vs \
         tmm={period_vented_tmm:.3} (closed={period_closed_causal:.3} vs tmm={period_closed_tmm:.3}), \
         causal_neg_e={causal_neg_e:.3e}, unconj_neg_e={unconj_neg_e:.6}"
    );
}

/// (f) G0: Boundary pin that `characteristic_line_dense(.., Some(8192))` and
/// `impedance_line(.., 8192, ..)` honor `n_fft = 8192`.
#[test]
fn convention_fixture_fft_length_is_honored() {
    let t0 = Instant::now();
    let gas = air20();
    let duct = vented_tone_hole_duct(true);
    let zc = gas.characteristic_impedance / (core::f64::consts::PI * 0.002 * 0.002);
    let r_line = characteristic_line_dense(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        N_FFT,
        zc,
        None,
        Some(N_FFT),
    )
    .expect("characteristic_line_dense realizes");
    let z_line = impedance_line(
        &duct,
        &gas,
        Termination::UnflangedOpen,
        RATE_HZ,
        N_FFT,
        None,
    )
    .expect("impedance_line realizes");
    let r_len = r_line.history().len();
    let z_len = z_line.history().len();
    let pass = r_len == N_FFT && z_len == N_FFT;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"def01-f-fft-length-honored\",\
         \"inputs\":\"analytic:vented_tone_hole_template\",\"termination\":\"unflanged_open\",\
         \"n_fft\":{N_FFT},\"pre_delay_energy\":0.0,\"neg_time_energy\":0.0,\
         \"dtft_rel_err\":0.0,\"observed\":{{\"char_n_fft\":{r_len},\"imp_n_fft\":{z_len}}},\
         \"oracle\":{N_FFT},\"err\":{},\"tol\":0,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        (r_len.abs_diff(N_FFT)).max(z_len.abs_diff(N_FFT)),
        if pass { "pass" } else { "fail" },
        t0.elapsed().as_micros()
    );
    assert_eq!(r_len, N_FFT);
    assert_eq!(z_len, N_FFT);
}
