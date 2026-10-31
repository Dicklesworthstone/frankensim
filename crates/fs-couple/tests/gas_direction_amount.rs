//! The declared gas reaches each self-oscillating family's resonator and its
//! sound, with a direction and an amount, through the render's own route
//! (.362 item 12, ruling 6349). Each family's shipped description is read by
//! the same entry point its `music_render` command uses; only the `ambient`
//! row, and where stated the clock or the open-end load, change.
//!
//! The amount is on the passive resonator the render builds from the
//! description and the gas record. With rigid walls and an ideally open end,
//! a lossless piecewise-cylinder or conic resonator's resonances scale exactly
//! with c, in the same mode, in every gas. To first order the losses shift
//! each by -(S_k + S_z):
//! - S_k, the mode-weighted mean of Kirchhoff's phase excess, lies between its
//!   values at the widest and narrowest sections;
//! - S_z, the Zc correction acting as an area perturbation, is at most half its
//!   largest value in size.
//! Both scale by each gas's own factor (the radius dependence is shared), so
//! the formant ratio over the c ratio lies between the corners of
//! (1 - kS_k - kzS_z)/(1 - S_k - S_z). Also allowed: second order, s^2 per gas
//! (twice a uniform tube's s^2/2), plus the stated resolution of the peak
//! finder. Every bound is printed before its assert.
//!
//! The rendered direction is asserted. The amount is printed against the c
//! ratio with no pass band: source-resonator pulling sets it, and its bound is
//! research-grade.
use fs_couple::brass_loop::performance::LipDescription;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::voice_performance::VoiceDescription;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use std::f64::consts::{PI, TAU};

/// Kirchhoff's wide-tube corrections at radius r and frequency f, from the
/// record's own rho, mu, gamma and Pr: (phase excess of k, real part of the
/// relative Zc correction), both (1 +- (gamma-1)/sqrt(Pr)) / (sqrt(2) r_v).
fn corrections(gas: &GasState, r: f64, f: f64) -> (f64, f64) {
    let thermal = (gas.gamma - 1.0) / gas.prandtl.sqrt();
    let scale = 1.0 / (2.0_f64.sqrt() * r * (gas.density * TAU * f / gas.dynamic_viscosity).sqrt());
    (scale * (1.0 + thermal), scale * (1.0 - thermal))
}

/// Autocorrelation f0 over a steady window, 60 Hz to `highest` Hz. A periodic
/// signal also correlates at two and three periods, so the period is the
/// shortest lag whose normalized score is within 5 % of the best (stated rule).
/// The period is not that whole number of samples: it is the vertex of the
/// parabola through the scores at that lag and its two neighbours.
///
/// What this replaced: the period was the integer lag itself, which is off by
/// up to half a sample, f0 / (2 lag) in f0 (0.13 Hz at 113 Hz and 48 kHz).
/// `harmonic_centroid` samples an unwindowed window at exact multiples of the
/// f0 it is given, so at harmonic k the comb sat k times that far from the true
/// line. Past about the 20th harmonic that is a whole bin of a 0.35 s window
/// (2.86 Hz), and the centroid, which weighs those lines, moved with the
/// rounding of the lag alone: a gas verdict could pass or fail on where the
/// integer lag fell, not on the sound.
fn f0(signal: &[f64], rate: f64, highest: f64) -> f64 {
    let lags = (rate / highest) as usize..=(rate / 60.0) as usize;
    let score = |lag: usize| {
        let (a, b) = (&signal[..signal.len() - lag], &signal[lag..]);
        a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>()
            / (a.iter().map(|x| x * x).sum::<f64>() * b.iter().map(|y| y * y).sum::<f64>()).sqrt()
    };
    let scores: Vec<(usize, f64)> = lags.map(|lag| (lag, score(lag))).collect();
    let best = scores.iter().map(|(_, s)| *s).fold(f64::NEG_INFINITY, f64::max);
    // The first local maximum within 5 % of the best, not the best itself.
    let at = scores.windows(3).position(|w| w[1].1 >= w[0].1 && w[1].1 >= w[2].1 && w[1].1 >= 0.95 * best)
        .map_or_else(|| scores.iter().enumerate().max_by(|p, q| p.1.1.total_cmp(&q.1.1)).unwrap().0, |i| i + 1);
    // Vertex of the parabola through (-1, l), (0, c), (1, r), in samples from the peak lag; none at an
    // edge of the searched range, or where the three scores do not bend downward.
    let offset = match (at.checked_sub(1).map(|i| scores[i].1), scores.get(at + 1).map(|s| s.1)) {
        (Some(l), Some(r)) => {
            let bend = l - 2.0 * scores[at].1 + r;
            if bend < 0.0 { 0.5 * (l - r) / bend } else { 0.0 }
        }
        _ => 0.0,
    };
    rate / (scores[at].0 as f64 + offset)
}

// ---- Voice: the Story (2008) [ɑ] description, VoiceDescription::from_bytes -> build.

const FIXTURE: &str = include_str!("../../../data/instruments/descriptions/voice/story2008-a-C3.performance");
const VOICE_AIR: &str = "ambient moist_air 310.15 101325 1.0";
const VOICE_HELIOX: &str = "ambient mixture 310.15 101325 mole 2 He 0.8 O2 0.2";
const VOICE_CO2: &str = "ambient mixture 310.15 101325 mole 1 CO2 1";
/// TractChart::tmm_formants_terminated samples |Z_in| every 2 Hz and refines
/// the largest sample's parabola, whose vertex stays within 1 Hz of it; the
/// true peak lies within one step of that sample, so the error is under 3 Hz.
const FORMANT_SCAN_HZ: f64 = 3.0;
/// Stated, not derived: the description tunes f0 in air to 3 cents, and a gas
/// reaches the folds only through the tract load, which the description says
/// moves f0 by a vowel-dependent amount. A semitone is that move's tolerance.
const F0_TOLERANCE: f64 = 0.059_463_094_359_295_26; // 2^(1/12) - 1

fn voice(ambient: &str, termination: &str) -> VoiceDescription {
    assert_eq!(FIXTURE.matches(VOICE_AIR).count(), 1);
    assert_eq!(FIXTURE.matches("termination flanged\n").count(), 1);
    let text = FIXTURE.replace(VOICE_AIR, ambient).replace("termination flanged\n", &format!("termination {termination}\n"));
    VoiceDescription::from_bytes(text.as_bytes()).unwrap()
}

/// Section radii of the description, glottis to lips.
fn tract_radii() -> Vec<f64> {
    FIXTURE.lines().filter_map(|line| line.strip_prefix("section ")).map(|row| {
        let area: f64 = row.split_whitespace().next().unwrap().parse().unwrap();
        (area / PI).sqrt()
    }).collect()
}


#[test]
fn the_tract_resonances_follow_the_declared_gas_by_its_sound_speed() {
    // A rigid tract (the description declares no wall) with an ideal open
    // mouth: its lossless resonances scale exactly with c, the same mode shape
    // in every gas. To first order, the losses shift each by -(S_k + S_z):
    // - S_k, a mode-weighted mean of the phase excess, lies between its values
    //   at the widest and narrowest sections;
    // - S_z, the Zc correction acting as an area perturbation, is at most half
    //   its largest value in size.
    // In the second gas both scale by that gas's own factor (the radius
    // dependence is shared), so the ratio of the two formants over the c ratio
    // lies between the corners of (1 - kS_k - kzS_z)/(1 - S_k - S_z). Allowed
    // on top: second order, s^2 per gas (twice a uniform tube's s^2/2), and the
    // finder's 3 Hz.
    let radii = tract_radii();
    let (r_min, r_max) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), r| (lo.min(*r), hi.max(*r)));
    let air = voice(VOICE_AIR, "ideal_open");
    let air_formants = air.build().unwrap().info().passive_formants_hz.clone().unwrap();
    assert!(air_formants.len() >= 2, "{air_formants:?}");
    for (name, record, up) in [("heliox", VOICE_HELIOX, true), ("co2", VOICE_CO2, false)] {
        let gas = voice(record, "ideal_open");
        let formants = gas.build().unwrap().info().passive_formants_hz.clone().unwrap();
        let (g, a) = (gas.ambient_gas(), air.ambient_gas());
        let c_ratio = g.sound_speed / a.sound_speed;
        assert_eq!(c_ratio > 1.0, up, "{name}: c {} against {}", g.sound_speed, a.sound_speed);
        let pairs = formants.len().min(air_formants.len());
        assert!(pairs >= 2, "{name}: {formants:?}");
        for n in 0..pairs {
            let (fg, fa) = (formants[n], air_formants[n]);
            // Direction: every resonance moves with c.
            assert_eq!(fg > fa, up, "{name} F{}: {fg} Hz against {fa} Hz", n + 1);
            let (k_ref, z_ref) = corrections(a, 1.0, fa);
            let (k_gas, z_gas) = corrections(g, 1.0, fg);
            let (kappa, kappa_z) = (k_gas / k_ref, z_gas / z_ref);
            let (k_lo, k_hi) = (corrections(a, r_max, fa).0, corrections(a, r_min, fa).0);
            let z_half = 0.5 * corrections(a, r_min, fa).1;
            let corners: Vec<f64> = [k_lo, k_hi].iter().flat_map(|&sk| [-z_half, z_half].map(move |sz|
                (1.0 - kappa * sk - kappa_z * sz) / (1.0 - sk - sz))).collect();
            let second = corrections(g, r_min, fg).0.powi(2) + corrections(a, r_min, fa).0.powi(2);
            let scan = FORMANT_SCAN_HZ / fg + FORMANT_SCAN_HZ / fa;
            let lo = corners.iter().copied().fold(f64::INFINITY, f64::min) - second - scan;
            let hi = corners.iter().copied().fold(0.0_f64, f64::max) + second + scan;
            let measured = (fg / fa) / c_ratio;
            println!("{{\"suite\":\"fs-couple\",\"case\":\"voice-resonance-amount\",\"gas\":\"{name}\",\"formant\":{},\"f_air_hz\":{fa:.2},\"f_gas_hz\":{fg:.2},\"c_ratio\":{c_ratio:.6},\"measured_over_c_ratio\":{measured:.6},\"band\":[{lo:.6},{hi:.6}],\"entry\":\"VoiceDescription::from_bytes -> build -> info().passive_formants_hz\"}}", n + 1);
            assert!(lo <= measured && measured <= hi, "{name} F{}: {measured} outside [{lo}, {hi}]", n + 1);
        }
    }
}


/// Amplitude-weighted mean frequency of the harmonics of f0 up to 5 kHz: the
/// formant region's position in the radiated sound.
fn harmonic_centroid(signal: &[f64], rate: f64, f0: f64) -> f64 {
    let (mut weighted, mut total) = (0.0, 0.0);
    let mut k = 1.0;
    while k * f0 < 5000.0 {
        let w = TAU * k * f0 / rate;
        let (re, im) = signal.iter().enumerate().fold((0.0, 0.0), |(re, im), (n, x)| {
            (re + x * (w * n as f64).cos(), im - x * (w * n as f64).sin())
        });
        let amplitude = re.hypot(im);
        weighted += k * f0 * amplitude;
        total += amplitude;
        k += 1.0;
    }
    weighted / total
}

fn render_voice(record: &str) -> Result<Vec<f64>, String> {
    // The shipped description, 0.75 s instead of 3 s; the same ramp to the same
    // drive, held, then released.
    let text = FIXTURE.replace(VOICE_AIR, record).replace("audio 48000 144000 2", "audio 48000 36000 2")
        .replace("subglottal 2.85 411.859\nsubglottal 3 0", "subglottal 0.7 411.859\nsubglottal 0.75 0");
    assert!(text.contains("subglottal 0.7 411.859\n"));
    let mut voice = VoiceDescription::from_bytes(text.as_bytes())?.build()?;
    let mut pressure = vec![0.0; 36_000];
    for block in pressure.chunks_mut(480) { voice.block(block).map_err(|e| format!("{e:?}"))?; }
    Ok(pressure)
}

#[test]
fn the_rendered_voice_moves_its_formant_region_with_the_gas_while_the_folds_keep_f0() {
    // Direction asserted; amount printed against the c ratio, with no pass
    // band: how far the sound's formant region moves is set by source-tract
    // pulling, whose bound is research-grade (ruling 6349).
    let rate = 48_000.0;
    let steady = 14_400..31_200; // 0.30-0.65 s
    let air = render_voice(VOICE_AIR).unwrap();
    let f0_air = f0(&air[steady.clone()], rate, 400.0);
    let centroid_air = harmonic_centroid(&air[steady.clone()], rate, f0_air);
    let c_air = voice(VOICE_AIR, "flanged").ambient_gas().sound_speed;
    for (name, record, up) in [("heliox", VOICE_HELIOX, true), ("co2", VOICE_CO2, false)] {
        let pressure = render_voice(record).unwrap_or_else(|e| panic!("{name}: the voice render refused the gas: {e}"));
        let f0_gas = f0(&pressure[steady.clone()], rate, 400.0);
        let centroid = harmonic_centroid(&pressure[steady.clone()], rate, f0_gas);
        let c_ratio = voice(record, "flanged").ambient_gas().sound_speed / c_air;
        println!("{{\"suite\":\"fs-couple\",\"case\":\"voice-rendered-direction\",\"gas\":\"{name}\",\"f0_air_hz\":{f0_air:.2},\"f0_gas_hz\":{f0_gas:.2},\"centroid_air_hz\":{centroid_air:.1},\"centroid_gas_hz\":{centroid:.1},\"centroid_ratio\":{:.4},\"c_ratio\":{c_ratio:.4},\"amount\":\"printed, no band\"}}",
            centroid / centroid_air);
        assert_eq!(centroid > centroid_air, up, "{name}: formant region {centroid} Hz against {centroid_air} Hz");
        assert!((f0_gas / f0_air - 1.0).abs() <= F0_TOLERANCE, "{name}: f0 {f0_gas} Hz against {f0_air} Hz");
    }
}

// ---- Lip: the shipped trumpet, composed as its render.sh does, LipDescription::from_bytes.

const HEAD: &str = include_str!("../../../data/instruments/descriptions/trumpet/trumpet-open.head");
const LIP_AIR: &str = "ambient moist_air 293.15 101325 0";
const LIP_HELIOX: &str = "ambient mixture 293.15 101325 mole 2 He 0.8 O2 0.2";
const LIP_CO2: &str = "ambient mixture 293.15 101325 mole 1 CO2 1";

/// render.sh's hold at lip pre-stress 1.42: the blowing pressure rises to
/// 1.4e4 Pa in 50 ms and is held.
fn compose(ambient: &str, rate: u32, seconds: f64, termination: &str) -> String {
    assert_eq!(HEAD.matches(LIP_AIR).count(), 1);
    assert_eq!(HEAD.matches("termination flanged\n").count(), 1);
    let head = HEAD.replace(LIP_AIR, ambient).replace("termination flanged\n", &format!("termination {termination}\n"));
    format!("frankensim-lip-performance-v1\naudio {rate} {}\n{head}schedule\nfrankensim-gesture-schedule-v1\n\
        control_rate_hz\t700\ntracks\t2\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t1\n\
        event\t0e0\t5e-2\tpressure-pa\t1.4e4\ntrack\tlip\tpre-stress\ninitial\tmultiplier\t1.42e0\nevents\t0\n",
        (f64::from(rate) * seconds).round() as u64)
}

fn lip(ambient: &str, rate: u32, seconds: f64, termination: &str) -> Result<LipDescription, String> {
    LipDescription::from_bytes(compose(ambient, rate, seconds, termination).as_bytes())
}

/// Every bore radius the head declares (cylinder, cone and bessel end radii).
fn bore_radii() -> Vec<f64> {
    HEAD.lines().flat_map(|line| {
        let words: Vec<&str> = line.split_whitespace().collect();
        let count = match words.first() { Some(&"cylinder") => 1, Some(&"cone" | &"bessel") => 2, _ => 0 };
        words.iter().skip(1).take(count).map(|w| w.parse::<f64>().unwrap()).collect::<Vec<_>>()
    }).collect()
}


/// The first `count` plane-wave |Z_in| peaks of the render's own bore
/// (`bore_scan`), located on a 2 Hz grid and refined by golden section.
fn bore_peaks(description: &LipDescription, count: usize) -> Vec<f64> {
    let magnitude = |f: f64| {
        let row = &description.bore_scan(&[f]).unwrap().rows[0];
        row.impedance_scalar[0].hypot(row.impedance_scalar[1])
    };
    let grid: Vec<f64> = (25..2000).map(|i| 2.0 * f64::from(i)).collect();
    let values: Vec<f64> = description.bore_scan(&grid).unwrap().rows.iter()
        .map(|row| row.impedance_scalar[0].hypot(row.impedance_scalar[1])).collect();
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    (1..grid.len() - 1).filter(|&i| values[i] > values[i - 1] && values[i] >= values[i + 1]).take(count).map(|i| {
        let (mut a, mut b) = (grid[i - 1], grid[i + 1]);
        while b - a > 1e-9 * grid[i] {
            let (c, d) = (b - ratio * (b - a), a + ratio * (b - a));
            if magnitude(c) > magnitude(d) { b = d } else { a = c }
        }
        0.5 * (a + b)
    }).collect()
}

#[test]
fn the_bore_resonances_follow_the_declared_gas_by_its_sound_speed() {
    // An ideally open bell: the lossless piecewise-conic bore's resonances
    // scale exactly with c, the same mode in every gas. The loss interval is
    // the voice test's (first order: the mode-weighted phase excess between
    // its widest- and narrowest-section values, the Zc correction at most
    // half its largest value; second order s^2 per gas). The peaks are
    // refined to 1e-9, so no scan term.
    let radii = bore_radii();
    let (r_min, r_max) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), r| (lo.min(*r), hi.max(*r)));
    let air = lip(LIP_AIR, 32_000, 0.1, "ideal_open").unwrap();
    let air_peaks = bore_peaks(&air, 6);
    assert_eq!(air_peaks.len(), 6, "{air_peaks:?}");
    for (name, record, up) in [("heliox", LIP_HELIOX, true), ("co2", LIP_CO2, false)] {
        let gas = lip(record, 32_000, 0.1, "ideal_open").unwrap();
        let peaks = bore_peaks(&gas, 6);
        let (g, a) = (gas.ambient_gas(), air.ambient_gas());
        let c_ratio = g.sound_speed / a.sound_speed;
        assert_eq!(c_ratio > 1.0, up);
        for (n, (&fg, &fa)) in peaks.iter().zip(&air_peaks).enumerate() {
            assert_eq!(fg > fa, up, "{name} peak {}: {fg} Hz against {fa} Hz", n + 1);
            let (kappa, kappa_z) = {
                let (kg, zg) = corrections(g, 1.0, fg);
                let (ka, za) = corrections(a, 1.0, fa);
                (kg / ka, zg / za)
            };
            let (k_lo, k_hi) = (corrections(a, r_max, fa).0, corrections(a, r_min, fa).0);
            let z_half = 0.5 * corrections(a, r_min, fa).1;
            let corners: Vec<f64> = [k_lo, k_hi].iter().flat_map(|&sk| [-z_half, z_half].map(move |sz|
                (1.0 - kappa * sk - kappa_z * sz) / (1.0 - sk - sz))).collect();
            let second = corrections(g, r_min, fg).0.powi(2) + corrections(a, r_min, fa).0.powi(2);
            let lo = corners.iter().copied().fold(f64::INFINITY, f64::min) - second;
            let hi = corners.iter().copied().fold(0.0_f64, f64::max) + second;
            let measured = (fg / fa) / c_ratio;
            println!("{{\"suite\":\"fs-couple\",\"case\":\"lip-bore-resonance-amount\",\"gas\":\"{name}\",\"peak\":{},\"f_air_hz\":{fa:.4},\"f_gas_hz\":{fg:.4},\"c_ratio\":{c_ratio:.6},\"measured_over_c_ratio\":{measured:.6},\"band\":[{lo:.6},{hi:.6}],\"entry\":\"LipDescription::from_bytes -> bore_scan\"}}", n + 1);
            assert!(lo <= measured && measured <= hi, "{name} peak {}: {measured} outside [{lo}, {hi}]", n + 1);
        }
    }
}


fn render_lip(record: &str) -> Vec<f64> {
    let mut performance = lip(record, 24_000, 1.0, "flanged").unwrap().build().unwrap();
    let mut pressure = vec![0.0; 24_000];
    for block in pressure.chunks_mut(240) { performance.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn the_played_pitch_moves_with_the_gas_in_the_rendered_sound() {
    // The render refuses CO2 by name, and that is a finding of this test, not
    // a skip. The shipped microphone integrates the 61 mm bell on 128 rings of
    // 512 points up to 6900 Hz. Its quadrature phase 2 pi f cell / c must stay
    // under pi/4: it is 0.776 in air (a 1 % margin) and 0.997 at CO2's 267 m/s.
    // The refusal is the receiver's admission, independent of the clock.
    let refused = lip(LIP_CO2, 24_000, 1.0, "flanged").and_then(|d| d.build().map(|_| ())).unwrap_err();
    assert!(refused.contains("quadrature underresolves the declared acoustic band"), "{refused}");
    // Direction asserted in heliox on a 24 kHz clock. The amount is printed
    // against the c ratio with no pass band: lip-bore pulling sets it, and its
    // bound is research-grade (ruling 6349).
    let rate = 24_000.0;
    let steady = 12_000..22_800; // 0.50-0.95 s
    let air = render_lip(LIP_AIR);
    let f0_air = f0(&air[steady.clone()], rate, 1500.0);
    let c_air = lip(LIP_AIR, 24_000, 1.0, "flanged").unwrap().ambient_gas().sound_speed;
    for (name, record, up) in [("heliox", LIP_HELIOX, true)] {
        let pressure = render_lip(record);
        let f0_gas = f0(&pressure[steady.clone()], rate, 1500.0);
        let c_ratio = lip(record, 24_000, 1.0, "flanged").unwrap().ambient_gas().sound_speed / c_air;
        println!("{{\"suite\":\"fs-couple\",\"case\":\"lip-rendered-direction\",\"gas\":\"{name}\",\"f0_air_hz\":{f0_air:.2},\"f0_gas_hz\":{f0_gas:.2},\"f0_ratio\":{:.4},\"c_ratio\":{c_ratio:.4},\"amount\":\"printed, no band\"}}",
            f0_gas / f0_air);
        assert_eq!(f0_gas > f0_air, up, "{name}: played {f0_gas} Hz against {f0_air} Hz");
    }
}

// ---- Flue: the shipped examples/flue.performance, FlueFilePerformance::from_bytes.

const FLUE: &str = include_str!("../examples/flue.performance");
const FLUE_AIR: &str = "gas 293.15 101325 0\n";
const FLUE_HELIOX: &str = "gas mixture 293.15 101325 mole 2 He 0.8 O2 0.2\n";
const FLUE_CO2: &str = "gas mixture 293.15 101325 mole 1 CO2 1\n";

fn flue(gas: &str) -> FlueFilePerformance {
    assert_eq!(FLUE.matches(FLUE_AIR).count(), 1);
    FlueFilePerformance::from_bytes(FLUE.replace(FLUE_AIR, gas).as_bytes(), 512, &CancelGate::new()).unwrap()
}

/// Per realized section of the render's own tube network: the one-way transit
/// it realizes (integer buffer plus the first-order allpass's DC delay), the
/// prediction L f_mech / c from the section's length and the record's c, and
/// the binary64 bound on their difference. The bound covers the allpass
/// coefficient a = (1-d)/(1+d) (three roundings) and its DC inverse
/// (1-a)/(1+a), which moves d by 2|da|/(1+a)^2, plus two roundings of the
/// request.
fn flue_transits(performance: &FlueFilePerformance) -> (f64, Vec<(f64, f64, f64)>) {
    let rate = f64::from(performance.renderer().mechanics_rate_hz());
    let network = performance.renderer().voice().network().expect("the shipped flue is a tube network");
    let c = network.section_gases()[0].sound_speed;
    assert!(network.section_gases().iter().all(|g| g.sound_speed == c), "one declared gas");
    let (_, sections) = network.fractional_sections().expect("the shipped flue selects the fractional image");
    let eps = f64::EPSILON;
    (c, sections.iter().map(|s| {
        let realized = s.integer_samples as f64 + s.transit.realized_dc_samples();
        let predicted = s.requested_length_m * rate / c;
        let d = s.transit.samples();
        let a = if d == 0.0 { 1.0 } else { (1.0 - d) / (1.0 + d) };
        let bound = 2.0 * 3.0 * eps * a / (1.0 + a).powi(2) + 4.0 * eps + 2.0 * eps * predicted;
        (realized, predicted, bound)
    }).collect())
}

#[test]
fn the_flue_tube_transits_follow_the_declared_gas_by_its_sound_speed() {
    let (c_air, air) = flue_transits(&flue(FLUE_AIR));
    assert!(!air.is_empty());
    for (name, record, up) in [("heliox", FLUE_HELIOX, true), ("co2", FLUE_CO2, false)] {
        let (c_gas, gas) = flue_transits(&flue(record));
        assert_eq!(c_gas > c_air, up);
        assert_eq!(gas.len(), air.len(), "{name}: the same lowered sections");
        for (i, (&(realized, predicted, bound), &(realized_air, _, _))) in gas.iter().zip(&air).enumerate() {
            println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-transit-amount\",\"gas\":\"{name}\",\"section\":{i},\"realized_samples\":{realized:.12},\"predicted_samples\":{predicted:.12},\"bound\":{bound:.3e},\"transit_ratio\":{:.9},\"inverse_c_ratio\":{:.9},\"entry\":\"FlueFilePerformance::from_bytes -> renderer().voice().network().fractional_sections()\"}}",
                realized / realized_air, c_air / c_gas);
            assert!((realized - predicted).abs() <= bound, "{name} section {i}: realized {realized} against {predicted} (bound {bound:e})");
            // Direction: a faster gas crosses the same section in fewer samples.
            assert_eq!(realized < realized_air, up, "{name} section {i}");
        }
    }
}

fn render_flue(record: &str) -> Vec<f64> {
    let mut renderer = flue(record).into_renderer();
    let mut pressure = vec![0.0; 24_000];
    for block in pressure.chunks_mut(480) { renderer.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn the_flue_sounds_higher_in_heliox_and_lower_in_carbon_dioxide() {
    // 0.5 s at 48 kHz, as shipped; steady window 0.20-0.34 s, before the
    // shipped release at 0.35 s. Direction asserted; amount printed against
    // the c ratio with no pass band (jet-tube pulling, research-grade).
    let rate = 48_000.0;
    let steady = 9_600..16_320;
    let rms = |p: &[f64]| (p.iter().map(|x| x * x).sum::<f64>() / p.len() as f64).sqrt();
    let air = render_flue(FLUE_AIR);
    assert!(rms(&air[steady.clone()]) > 0.0);
    let f0_air = f0(&air[steady.clone()], rate, 2500.0);
    let c_air = flue(FLUE_AIR).renderer().voice().network().unwrap().section_gases()[0].sound_speed;
    for (name, record, up) in [("heliox", FLUE_HELIOX, true), ("co2", FLUE_CO2, false)] {
        let pressure = render_flue(record);
        let level = rms(&pressure[steady.clone()]) / rms(&air[steady.clone()]);
        let f0_gas = f0(&pressure[steady.clone()], rate, 2500.0);
        let c_gas = flue(record).renderer().voice().network().unwrap().section_gases()[0].sound_speed;
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-rendered-direction\",\"gas\":\"{name}\",\"f0_air_hz\":{f0_air:.2},\"f0_gas_hz\":{f0_gas:.2},\"f0_ratio\":{:.4},\"c_ratio\":{:.4},\"rms_over_air\":{level:.3e},\"amount\":\"printed, no band\"}}",
            f0_gas / f0_air, c_gas / c_air);
        assert!(level > 1e-3, "{name}: the flue fell silent (rms {level:e} of air's)");
        assert_eq!(f0_gas > f0_air, up, "{name}: sounded {f0_gas} Hz against {f0_air} Hz");
    }
}
