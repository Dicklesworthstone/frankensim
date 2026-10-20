//! IRP M2.5-W3b G3: the passive, delay-extracted multiport radiation model in
//! a closed reed + bore + radiation loop, against the same loop with the exact
//! radiation matrix.
//!
//! Declared synthetic fixture, no instrument claim: a cylindrical bore with one
//! open side chimney, whose bell and chimney mouths end flush in one infinite
//! planar baffle 45 mm apart. There the coplanar baffled-piston matrix (fs-bem
//! `baffled_piston_array_impedance`: closed-form self terms, converged Rayleigh
//! mutual terms) is the exact Z_rad at every frequency. The model is fitted
//! from 200 Hz to 8 kHz (relaxed vector fitting collapsed at every order on
//! a 10 Hz start, where the mass-like low band dominates the samples). The
//! loop runs at 16 kHz and samples the load from its lowest line bin
//! (0.24 Hz at 2^16 taps), so grid passivity is repaired from 1 Hz to 8 kHz
//! (grid 0.1 Hz to 80 kHz). A fit-band-only repair left an active Hermitian
//! part at a low bin, which the coupled-load admission refused. The repair's
//! change is gated over the fit band. The realizer's near-DC point
//! (1e-3 rad/s) falls back to R = -1 for both sources on a refusal.
//! The bake passes fs-vfit's dispersion-relation causality test, with its
//! gap delays removed, before it is fitted.
//! "Loaded modes" are the radiation-loaded resonances of the
//! realized line below the realizer's 0.5-Nyquist taper; the valve enters
//! through the oscillation threshold and the rendered sound.

use std::collections::HashMap;

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::aperture::{BaffledPiston, baffled_piston_array_impedance};
use fs_couple::driving_point::LineLadder;
use fs_couple::driving_point::network::NetworkDrivingPoint;
use fs_couple::render::ReedBoreVoice;
use fs_couple::thin_plate::PlateBank;
use fs_duct::network::{
    CoupledTermination, DuctNetwork, JunctionCorrections, NetworkEdge, NetworkNode, NodeBoundary,
};
use fs_duct::{Duct, DuctError, LossModel, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_math::det;
use fs_scenario::BeatingReed;
use fs_vfit::FitOptions;
use fs_vfit::causality::{CausalityOptions, check_mimo_causality};
use fs_vfit::discretize::DelayedFilter;
use fs_vfit::mimo::{DelayedMimoModel, delayed_mimo_vector_fit, repair_delayed_mimo_passivity};

const RATE_HZ: u32 = 16_000;
const BORE_RADIUS: f64 = 0.0075;
const HOLE_RADIUS: f64 = 0.004;
const SPACING: f64 = 0.045;
const BAND_HZ: (f64, f64) = (200.0, 8000.0);
const TAU: f64 = core::f64::consts::TAU;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn area(radius: f64) -> f64 {
    core::f64::consts::PI * radius * radius
}

/// Exact acoustic matrix [Pa s/m^3] under `exp(-i w t)`, bell then hole:
/// the mechanical piston matrix over `A_a A_b`.
fn exact(omega: f64, medium: Medium) -> Result<Vec<C64>, DuctError> {
    let pistons = [
        BaffledPiston {
            center_m: [0.0, 0.0],
            radius_m: BORE_RADIUS,
        },
        BaffledPiston {
            center_m: [SPACING, 0.0],
            radius_m: HOLE_RADIUS,
        },
    ];
    let z = baffled_piston_array_impedance(&pistons, omega / medium.sound_speed, medium).map_err(
        |_| DuctError::BadParameter {
            what: "exact baffled-piston matrix refused",
        },
    )?;
    Ok((0..4)
        .map(|i| z[i].scale(1.0 / (area(pistons[i / 2].radius_m) * area(pistons[i % 2].radius_m))))
        .collect())
}

fn coupled(omega: f64, z: Vec<C64>, source_id: &str) -> CoupledTermination {
    CoupledTermination {
        omega,
        labels: vec!["bell".into(), "hole".into()],
        z,
        source_id: source_id.into(),
    }
}

fn model_matrix(model: &DelayedMimoModel, omega: f64) -> Vec<C64> {
    // fs-vfit's s = i w axis back onto exp(-i w t).
    model.eval_iw(omega).into_iter().map(|z| z.conj()).collect()
}

fn driving_point() -> NetworkDrivingPoint {
    use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
    let id = |s: &str| StableId::new(s).unwrap();
    let schema = PortKind::FluidPressureFlux
        .scalar_seed_schema(
            id("closed-loop-radiation-audit"),
            CoordinateBinding::new(
                id("pressure-volume-flow"),
                id("mouthpiece-plane"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(id("declared-frequency-sweep"), 0),
        )
        .unwrap();
    let node = |name: &str, boundary| NetworkNode {
        id: name.into(),
        boundary,
    };
    let port = |name: &str| node(name, NodeBoundary::Port { label: name.into() });
    let edge = |name: &str, from, to, radius, length| NetworkEdge {
        id: name.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    };
    NetworkDrivingPoint::new_coupled(
        DuctNetwork {
            nodes: vec![
                port("in"),
                node("junction", NodeBoundary::Junction),
                port("bell"),
                port("hole"),
            ],
            edges: vec![
                edge("upper", 0, 1, BORE_RADIUS, 0.25),
                edge("lower", 1, 2, BORE_RADIUS, 0.06),
                edge("chimney", 1, 3, HOLE_RADIUS, 0.006),
            ],
        },
        "in".into(),
        Vec::new(),
        vec!["bell".into(), "hole".into()],
        schema,
    )
    .unwrap()
}

/// The first shared-pole delayed model of the exact matrix that meets the
/// 1e-3 matrix-relative operator tolerance, residue-repaired to grid passivity.
fn fitted_model(medium: Medium) -> DelayedMimoModel {
    let omegas: Vec<f64> = (0..160)
        .map(|i| TAU * BAND_HZ.0 * det::exp(det::ln(BAND_HZ.1 / BAND_HZ.0) * f64::from(i) / 159.0))
        .collect();
    let data: Vec<Vec<C64>> = omegas
        .iter()
        .map(|&w| {
            exact(w, medium)
                .unwrap()
                .into_iter()
                .map(|z| z.conj())
                .collect()
        })
        .collect();
    // The mutual response starts with the closest-point arrival across the
    // gap; extracting the centre distance would leave an advanced, non-causal
    // remainder of (a + b) / c.
    let tau = (SPACING - BORE_RADIUS - HOLE_RADIUS) / medium.sound_speed;
    let delays = [0.0, tau, tau, 0.0];
    // The bake must pass the dispersion-relation causality test, delays
    // removed, before it is fitted; a failure refuses here. Beyond the band a
    // piston entry stays below 3x its in-band maximum (self terms tend to
    // rho c / A with a 1.17 overshoot); the mutual ladder is good to 1e-9.
    let step = TAU * BAND_HZ.1 / 512.0;
    let uniform: Vec<f64> = (0..512).map(|j| (f64::from(j) + 0.5) * step).collect();
    let samples: Vec<Vec<C64>> = uniform
        .iter()
        .map(|&w| {
            exact(w, medium)
                .unwrap()
                .into_iter()
                .map(|z| z.conj())
                .collect()
        })
        .collect();
    let causality = check_mimo_causality(
        &uniform,
        &samples,
        2,
        &delays,
        &CausalityOptions {
            subtractions: 4,
            tail_bound_ratio: 3.0,
            relative_data_tolerance: 1.0e-9,
        },
    )
    .unwrap();
    for (entry, report) in causality.iter().enumerate() {
        println!(
            "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"causality\",\"entry\":{entry},\"tested\":{},\"worst_ratio\":{:.3e}}}",
            report.tested, report.worst_ratio
        );
    }
    let band = (omegas[0], omegas[omegas.len() - 1]);
    for order in [6, 10, 14, 18, 24] {
        let opts = FitOptions {
            fit_d: false,
            ..FitOptions::new(order)
        };
        // A refused rung (for example a collapsed relaxation sigma) is
        // reported, and the ladder moves to the next order.
        let (model, report) = match delayed_mimo_vector_fit(&omegas, &data, 2, &delays, &opts, true)
        {
            Ok(fit) => fit,
            Err(error) => {
                println!(
                    "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"fit\",\"order\":{order},\"refused\":\"{error}\"}}"
                );
                continue;
            }
        };
        println!(
            "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"fit\",\"order\":{order},\"max_rel\":{:.3e},\"stable\":{}}}",
            report.max_relative_error,
            model.model.is_stable()
        );
        if report.max_relative_error <= 1.0e-3 && model.model.is_stable() {
            // The loop samples the load from its lowest line bin up, so
            // passivity is repaired from 1 Hz (grid 0.1 Hz to 80 kHz); the
            // repair's change is gated over the fit band itself.
            let (repaired, repair) =
                repair_delayed_mimo_passivity(&model, (TAU * 1.0, band.1)).unwrap();
            let frobenius = |z: &[C64]| det::sqrt(z.iter().map(|v| v.norm_sq()).sum());
            let change = omegas
                .iter()
                .map(|&w| {
                    let (a, b) = (model.eval_iw(w), repaired.eval_iw(w));
                    let gap: Vec<C64> = a.iter().zip(&b).map(|(&x, &y)| x - y).collect();
                    frobenius(&gap) / frobenius(&a)
                })
                .fold(0.0, f64::max);
            println!(
                "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"repair\",\"rounds\":{},\"fit_band_change\":{change:.3e},\"passive\":{}}}",
                repair.rounds, repair.certificate.passive
            );
            assert!(repair.certificate.passive);
            assert!(
                change <= 1.0e-4,
                "repair change {change:.3e} in the fit band"
            );
            return repaired;
        }
    }
    panic!("no model order up to 24 meets the 1e-3 operator tolerance");
}

/// Golden-section maximum of `value` on `[lo, hi]`.
fn maximize(lo: f64, hi: f64, value: &mut dyn FnMut(f64) -> f64) -> f64 {
    let ratio = 0.5 * (det::sqrt(5.0) - 1.0);
    let (mut a, mut b) = (lo, hi);
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
}

/// Half-power bandwidth [Hz] of the peak at `f` of `power`, when both
/// half-power points lie within 20 % of `f`.
fn bandwidth(f: f64, power: &mut dyn FnMut(f64) -> f64) -> Option<f64> {
    let half = 0.5 * power(f);
    let mut edge = |outer: f64| -> Option<f64> {
        if power(outer) >= half {
            return None;
        }
        let (mut inside, mut outside) = (f, outer);
        for _ in 0..50 {
            let mid = 0.5 * (inside + outside);
            if power(mid) >= half {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        Some(0.5 * (inside + outside))
    };
    Some(edge(1.2 * f)? - edge(0.8 * f)?)
}

fn cents(a: f64, b: f64) -> f64 {
    1200.0 * det::ln(a / b) / core::f64::consts::LN_2
}

fn taps(line: &DelayedFilter, n: usize) -> Vec<f64> {
    let mut probe = line.clone();
    (0..n)
        .map(|k| probe.push(if k == 0 { 1.0 } else { 0.0 }).unwrap())
        .collect()
}

/// `|Z / zc|^2` of a realized line from the DTFT of its taps (the realizer's
/// `exp(+i w n dt)` convention; the magnitude is conjugation invariant).
fn line_power(ir: &[f64], f: f64) -> f64 {
    let theta = TAU * f / f64::from(RATE_HZ);
    let r = ir.iter().enumerate().fold(C64::ZERO, |acc, (n, &h)| {
        let (s, c) = det::sin_cos(theta * n as f64);
        acc + C64::new(h * c, h * s)
    });
    let num = C64::ONE + r;
    let den = C64::ONE - r;
    (num.re * num.re + num.im * num.im) / (den.re * den.re + den.im * den.im)
}

fn reed(blowing_pressure_pa: f64) -> BeatingReed {
    BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    }
}

fn ac_rms(x: &[f64]) -> f64 {
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    det::sqrt(x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / x.len() as f64)
}

/// Fundamental [Hz] and harmonic levels [dB] below the band edge, from a
/// Hann-windowed spectrum of the last 8192 samples.
fn harmonics(signal: &[f64]) -> (f64, Vec<f64>) {
    let n = 8192;
    let tail = &signal[signal.len() - n..];
    let windowed: Vec<f64> = tail
        .iter()
        .enumerate()
        .map(|(i, &v)| v * 0.5 * (1.0 - det::sin_cos(TAU * i as f64 / n as f64).1))
        .collect();
    let spectrum: Vec<f64> = fs_fft::RealFft::new(n)
        .forward(&windowed)
        .iter()
        .map(|c| det::sqrt(c.re * c.re + c.im * c.im))
        .collect();
    let df = f64::from(RATE_HZ) / n as f64;
    // Fundamental: the strongest bin from 100 Hz to 2 kHz, parabolic in log.
    let k = (51..1025)
        .max_by(|&a, &b| spectrum[a].total_cmp(&spectrum[b]))
        .unwrap();
    let (l, c, r) = (
        det::ln(spectrum[k - 1]),
        det::ln(spectrum[k]),
        det::ln(spectrum[k + 1]),
    );
    let f0 = (k as f64 + 0.5 * (l - r) / (l - 2.0 * c + r)) * df;
    let levels = (1..)
        .map(|h| f64::from(h) * f0)
        .take_while(|&f| f < BAND_HZ.1 - 4.0 * df)
        .map(|f| {
            let centre = (f / df).round() as usize;
            let peak = spectrum[centre - 2..=centre + 2]
                .iter()
                .fold(0.0_f64, |m, &v| m.max(v));
            20.0 * det::ln(peak) / core::f64::consts::LN_10
        })
        .collect();
    (f0, levels)
}

#[test]
fn closed_loop_audit_reed_bore_radiation() {
    let gas = gas();
    let medium = Medium::from(&gas);
    let point = driving_point();
    let model = fitted_model(medium);
    let loss = LossModel::Bessel;
    let mut cache: HashMap<u64, Vec<C64>> = HashMap::new();
    let mut exact_at = |omega: f64| -> Result<CoupledTermination, DuctError> {
        let z = match cache.get(&omega.to_bits()) {
            Some(z) => z.clone(),
            None => {
                let z = exact(omega, medium)?;
                cache.insert(omega.to_bits(), z.clone());
                z
            }
        };
        Ok(coupled(omega, z, "exact coplanar baffled-piston matrix"))
    };
    let mut model_at = |omega: f64| -> Result<CoupledTermination, DuctError> {
        Ok(coupled(
            omega,
            model_matrix(&model, omega),
            "repaired delayed MIMO model",
        ))
    };

    // Complex response: input impedance peaks (1 cent) and levels (0.5 dB).
    let mut z_exact = |f: f64| {
        let c = exact_at(TAU * f).unwrap();
        point
            .response_coupled(&gas, TAU * f, loss, &c)
            .unwrap()
            .impedance
            .abs()
    };
    let mut z_model = |f: f64| {
        let c = model_at(TAU * f).unwrap();
        point
            .response_coupled(&gas, TAU * f, loss, &c)
            .unwrap()
            .impedance
            .abs()
    };
    let scan: Vec<f64> = (0..600)
        .map(|i| 100.0 * det::exp(det::ln(80.0) * f64::from(i) / 599.0))
        .collect();
    let values: Vec<f64> = scan.iter().map(|&f| z_exact(f)).collect();
    let mut peaks = Vec::new();
    for i in 1..scan.len() - 1 {
        if values[i] > values[i - 1] && values[i] > values[i + 1] {
            let f_exact = maximize(scan[i - 1], scan[i + 1], &mut z_exact);
            let f_model = maximize(scan[i - 1], scan[i + 1], &mut z_model);
            let shift = cents(f_model, f_exact);
            let level =
                20.0 * det::ln(z_model(f_model) / z_exact(f_exact)) / core::f64::consts::LN_10;
            println!(
                "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"input-impedance-peak\",\"exact_hz\":{f_exact:.3},\"cents\":{shift:.5},\"db\":{level:.5}}}"
            );
            assert!(
                shift.abs() <= 1.0 && level.abs() <= 0.5,
                "{f_exact} Hz: {shift} c, {level} dB"
            );
            peaks.push(f_exact);
        }
    }
    assert!(peaks.len() >= 4, "input impedance peaks in band: {peaks:?}");

    // Loaded modes and decay rates of the realized model line (1 cent, 2 %).
    let zc = gas.characteristic_impedance / area(BORE_RADIUS);
    let ladder = LineLadder::default();
    let (line, receipt) = point
        .realize_line(&gas, loss, RATE_HZ, zc, &ladder, Some(&mut model_at))
        .unwrap();
    let ir = taps(&line, receipt.n_fft);
    let mut audited = 0;
    for &f_exact in peaks.iter().filter(|&&f| f < 0.25 * f64::from(RATE_HZ)) {
        let mut power_line = |f: f64| line_power(&ir, f);
        let f_line = maximize(0.99 * f_exact, 1.01 * f_exact, &mut power_line);
        let mut power_exact = |f: f64| {
            let z = z_exact(f);
            z * z
        };
        let (Some(width_line), Some(width_exact)) = (
            bandwidth(f_line, &mut power_line),
            bandwidth(f_exact, &mut power_exact),
        ) else {
            continue;
        };
        let shift = cents(f_line, f_exact);
        let decay = width_line / width_exact - 1.0;
        println!(
            "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"loaded-mode\",\"n_fft\":{},\"exact_hz\":{f_exact:.3},\"cents\":{shift:.5},\"decay_rel\":{decay:.3e}}}",
            receipt.n_fft
        );
        assert!(
            shift.abs() <= 1.0 && decay.abs() <= 0.02,
            "{f_exact} Hz: {shift} c, {decay}"
        );
        audited += 1;
    }
    assert!(audited >= 3, "loaded modes audited: {audited}");

    // Valve threshold (1 %) and rendered sound in the closed loop.
    let render = |source: &mut dyn FnMut(f64) -> Result<CoupledTermination, DuctError>,
                  pressure: f64,
                  samples: usize| {
        let mut voice = ReedBoreVoice::new_with_network(
            &point,
            &gas,
            loss,
            reed(pressure),
            PlateBank::default(),
            1.0,
            RATE_HZ,
            &ladder,
            Some(source),
        )
        .unwrap();
        let mut out = vec![0.0; samples];
        voice.step_block(&mut out).unwrap();
        out
    };
    let threshold = |source: &mut dyn FnMut(f64) -> Result<CoupledTermination, DuctError>| {
        let mut sounding = |pressure: f64| {
            let p = render(&mut *source, pressure, 9600);
            let (late, later) = (ac_rms(&p[4800..7200]), ac_rms(&p[7200..]));
            later > 0.01 * pressure || (later > 1.0e-6 * pressure && later > late)
        };
        let (mut lo, mut hi) = (1000.0, 5000.0);
        assert!(!sounding(lo) && sounding(hi), "threshold bracket");
        for _ in 0..12 {
            let mid = 0.5 * (lo + hi);
            if sounding(mid) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let threshold_exact = threshold(&mut exact_at);
    let threshold_model = threshold(&mut model_at);
    let threshold_rel = threshold_model / threshold_exact - 1.0;
    let pressure = 1.3 * threshold_exact;
    let (f0_exact, levels_exact) = harmonics(&render(&mut exact_at, pressure, 16_000));
    let (f0_model, levels_model) = harmonics(&render(&mut model_at, pressure, 16_000));
    // Harmonics within 60 dB of the exact loop's strongest, on both spectra.
    let strongest = levels_exact
        .iter()
        .fold(f64::NEG_INFINITY, |m, &v| m.max(v));
    let compared: Vec<(f64, f64)> = levels_exact
        .iter()
        .zip(&levels_model)
        .filter(|(a, _)| **a > strongest - 60.0)
        .map(|(&a, &b)| (a, b))
        .collect();
    let count = compared.len();
    let distance =
        det::sqrt(compared.iter().map(|(a, b)| (a - b) * (a - b)).sum::<f64>() / count as f64);
    let pitch = cents(f0_model, f0_exact);
    println!(
        "{{\"suite\":\"fs-couple-w3b-g3\",\"case\":\"closed-loop\",\"threshold_exact_pa\":{threshold_exact:.2},\"threshold_rel\":{threshold_rel:.3e},\"f0_hz\":{f0_exact:.3},\"pitch_cents\":{pitch:.4},\"harmonics\":{count},\"harmonic_rms_db\":{distance:.4}}}"
    );
    assert!(
        threshold_rel.abs() <= 0.01,
        "threshold differs by {threshold_rel}"
    );
    assert!(count >= 5, "harmonics compared: {count}");
    assert!(pitch.abs() <= 1.0, "playing frequency differs by {pitch} c");
    assert!(
        distance <= 0.5,
        "harmonic spectrum differs by {distance} dB rms"
    );
}
