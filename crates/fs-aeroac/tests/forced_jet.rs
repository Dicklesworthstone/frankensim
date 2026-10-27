//! The forced jet's causal response (fs_aeroac::forced_jet): the published
//! spatial Bickley mode downstream, the step kernel against independent
//! harmonic forcing, a three-rung grid ladder, box and band invariance, the
//! audio-rate hat weights, and the causality admission against a pure delay.
//! Run with --nocapture: every number the flue cards will use is printed.

use fs_aeroac::forced_jet::{
    ForcedJet, ForcedJetResolution, JetStepResponse, bickley_step_response,
};
use fs_fft::{C64 as Bin, RealFft};
use fs_math::c64::C64;
use fs_math::det;

const STATION: f64 = 11.0;
const WORKERS: usize = 8;

fn resolution(
    spacing: f64,
    cells: usize,
    half_height: f64,
    window: f64,
    reynolds: f64,
) -> ForcedJetResolution {
    ForcedJetResolution {
        spacing,
        streamwise_cells: cells,
        half_height,
        window,
        courant: 0.8,
        reynolds,
        workers: WORKERS,
    }
}

fn polar(g: C64) -> (f64, f64) {
    (g.abs(), det::atan2(g.im, g.re).to_degrees())
}

/// Harmonic forcing `v_ac = env(t) cos(w t)` from rest; returns the complex
/// displacement over velocity at `stations`, projected on `e^(-i w t)` over
/// whole periods ending at the window.
fn harmonic(
    res: ForcedJetResolution,
    omega: f64,
    ramp_periods: f64,
    periods: f64,
    stations: &[f64],
) -> Vec<C64> {
    let mut jet = ForcedJet::new(res).expect("solver");
    let period = 2.0 * core::f64::consts::PI / omega;
    let ramp = ramp_periods * period;
    let a = core::f64::consts::PI / ramp;
    let env = |t: f64| {
        if t >= ramp {
            1.0
        } else {
            0.5 * (1.0 - det::cos(a * t))
        }
    };
    let vac = |t: f64| env(t) * det::cos(omega * t);
    // Exact antiderivative of the raised-cosine ramp times cos(w t) (a != w).
    let ramp_part = |t: f64| {
        0.5 * det::sin(omega * t) / omega
            - 0.25
                * (det::sin((omega + a) * t) / (omega + a)
                    + det::sin((omega - a) * t) / (omega - a))
    };
    let xi = |t: f64| {
        if t < ramp {
            ramp_part(t)
        } else {
            ramp_part(ramp) + (det::sin(omega * t) - det::sin(omega * ramp)) / omega
        }
    };
    let dt = res.time_step();
    let steps = (res.window / dt).ceil() as usize;
    let t_end = steps as f64 * dt;
    let t_start = t_end - periods * period;
    let mut acc = vec![C64::new(0.0, 0.0); stations.len()];
    let mut vacc = C64::new(0.0, 0.0);
    let mut prev: Option<(Vec<C64>, C64)> = None;
    let forcing = |t: f64| (vac(t), xi(t));
    for _ in 0..steps {
        jet.step(&forcing).expect("step");
        let t = jet.time();
        if t > t_start {
            let (s, c) = det::sin_cos(omega * t);
            let w = C64::new(c, -s);
            let eta: Vec<C64> = stations
                .iter()
                .map(|&x| w.scale(jet.displacement_at(x).expect("station")))
                .collect();
            let v = w.scale(vac(t));
            if let Some((pe, pv)) = &prev {
                for (a, (e, p)) in acc.iter_mut().zip(eta.iter().zip(pe)) {
                    *a = *a + (*e + *p).scale(0.5 * dt);
                }
                vacc = vacc + (v + *pv).scale(0.5 * dt);
            }
            prev = Some((eta, v));
        }
    }
    acc.into_iter().map(|a| a / vacc).collect()
}

const STROUHALS: [f64; 10] = [0.03, 0.04, 0.06, 0.1, 0.16, 0.3, 0.5, 1.0, 2.0, 3.0];

fn print_transfer(label: &str, s: &JetStepResponse) -> Vec<C64> {
    let gs: Vec<C64> = STROUHALS
        .iter()
        .map(|&w| s.transfer(w).expect("transfer"))
        .collect();
    let row: Vec<String> = STROUHALS
        .iter()
        .zip(&gs)
        .map(|(w, g)| {
            let (m, p) = polar(*g);
            format!("{w}:{m:.3}/{p:.1}")
        })
        .collect();
    println!("{label} held {:.4} | {}", s.held_gain(), row.join(" "));
    gs
}

#[test]
fn exit_receptivity_grows_into_the_published_spatial_bickley_mode() {
    // Sabatini & Bailly (2015) Table 1, sinuous Bickley spatial mode at w = 0.8.
    let pin = C64::new(1.449709, -0.134109);
    let res = resolution(0.1, 2048, 8.0, 125.0, f64::INFINITY);
    let stations: Vec<f64> = (0..=40).map(|k| 20.0 + 0.5 * k as f64).collect();
    let g = harmonic(res, 0.8, 2.0, 4.0, &stations);
    // Least squares of log G = c - i alpha x over the stations.
    let n = stations.len() as f64;
    let mut phase = Vec::with_capacity(stations.len());
    let mut last = 0.0;
    for (k, gk) in g.iter().enumerate() {
        let mut p = det::atan2(gk.im, gk.re);
        if k > 0 {
            while p - last > core::f64::consts::PI {
                p -= 2.0 * core::f64::consts::PI;
            }
            while p - last < -core::f64::consts::PI {
                p += 2.0 * core::f64::consts::PI;
            }
        }
        phase.push(p);
        last = p;
    }
    let mx = stations.iter().sum::<f64>() / n;
    let slope = |ys: &[f64]| {
        let my = ys.iter().sum::<f64>() / n;
        let num: f64 = stations
            .iter()
            .zip(ys)
            .map(|(x, y)| (x - mx) * (y - my))
            .sum();
        let den: f64 = stations.iter().map(|x| (x - mx) * (x - mx)).sum();
        num / den
    };
    let logmag: Vec<f64> = g.iter().map(|gk| det::ln(gk.abs())).collect();
    // The projection on e^(-i w t) conjugates the published e^(i(alpha x - w t))
    // convention: G ~ e^(-i conj(alpha) x), so alpha = -d(phase)/dx - i d(log|G|)/dx.
    let alpha = C64::new(-slope(&phase), -slope(&logmag));
    println!(
        "spatial fit x in [20, 40] at w = 0.8: alpha = {:.6} {:+.6}i (pin {:.6} {:+.6}i)",
        alpha.re, alpha.im, pin.re, pin.im
    );
    // The wavenumber is exact to 0.1 %. The growth rate converges from below
    // as the critical layer (width c_i / |U'| ~ 0.07 b) is resolved: 3.4 %
    // low at dx 0.1 and 1.2 % at dx 0.05 (prototype, fit [20, 40], window
    // 200), 4.4 % here with this shorter window. The test runs the production
    // spacing and holds that growth to 6 %.
    assert!(
        (alpha.re - pin.re).abs() < 2e-3 * pin.re,
        "downstream wave does not carry the published wavenumber: {alpha:?}"
    );
    assert!(
        (alpha.im - pin.im).abs() < 0.06 * pin.im.abs(),
        "downstream wave does not grow at the published rate: {alpha:?}"
    );
}

#[test]
fn step_kernel_is_the_harmonic_steady_state_up_to_its_held_tail() {
    let res = resolution(0.1, 1024, 8.0, 90.0, 600.0);
    let step = bickley_step_response(STATION, &res).expect("step response");
    // The harmonic side is converged: 16 and 32 periods, and windows 90 and 150, agree to 1e-3.
    // Holding the step response after window 90 cuts a slowly decaying tail. Against harmonic
    // forcing that costs 2.5-3.9 % at w 1-2 and 1.2 % at w 3; at w 2 window 150 halves it, while
    // spacing 0.2/0.1/0.05 and Courant 0.8/0.4 leave it at 3.2-3.8 %.
    for (omega, ramp, periods, bound) in [(0.5, 2.0, 4.0, 0.01), (2.0, 4.0, 16.0, 0.04)] {
        let from_step = step.transfer(omega).expect("transfer");
        let forced = harmonic(res, omega, ramp, periods, &[STATION])[0];
        let rel = (from_step - forced).abs() / forced.abs();
        let (ms, ps) = polar(from_step);
        let (mh, ph) = polar(forced);
        println!(
            "w {omega}: step kernel {ms:.4}/{ps:.2} deg, harmonic forcing {mh:.4}/{ph:.2} deg, rel {rel:.2e}"
        );
        assert!(
            rel < bound,
            "the step kernel is not the harmonic response at w = {omega}"
        );
    }
}

#[test]
fn kernel_converges_on_a_grid_ladder_and_ignores_the_worker_count() {
    let window = 90.0;
    let rungs = [(0.2, 512), (0.1, 1024), (0.05, 2048)];
    let mut transfers = Vec::new();
    for (spacing, cells) in rungs {
        let res = resolution(spacing, cells, 8.0, window, 600.0);
        let t0 = std::time::Instant::now();
        let s = bickley_step_response(STATION, &res).expect("step response");
        let label = format!("dx {spacing} ({:.1} s)", t0.elapsed().as_secs_f64());
        transfers.push(print_transfer(&label, &s));
    }
    for (k, w) in STROUHALS.iter().enumerate() {
        let coarse = (transfers[0][k] - transfers[1][k]).abs() / transfers[1][k].abs();
        let fine = (transfers[1][k] - transfers[2][k]).abs() / transfers[2][k].abs();
        let order = (coarse / fine).log2();
        println!(
            "w {w}: rel change 0.2->0.1 {coarse:.4}, 0.1->0.05 {fine:.4}, observed order {order:.2}"
        );
        if (0.1..=0.5).contains(w) {
            assert!(
                fine < 0.03 && fine < coarse,
                "kernel not converging at w = {w}"
            );
        }
    }
    let one = bickley_step_response(
        STATION,
        &ForcedJetResolution {
            workers: 1,
            ..resolution(0.2, 512, 8.0, window, 600.0)
        },
    )
    .expect("serial");
    let many = bickley_step_response(
        STATION,
        &ForcedJetResolution {
            workers: 5,
            ..resolution(0.2, 512, 8.0, window, 600.0)
        },
    )
    .expect("parallel");
    assert!(
        one.values()
            .iter()
            .zip(many.values())
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "worker count changed the bits"
    );
}

#[test]
fn box_length_and_band_height_leave_the_kernel_unchanged() {
    let base =
        bickley_step_response(STATION, &resolution(0.1, 1024, 8.0, 90.0, 600.0)).expect("base");
    let long =
        bickley_step_response(STATION, &resolution(0.1, 2048, 8.0, 90.0, 600.0)).expect("long box");
    let narrow = bickley_step_response(STATION, &resolution(0.1, 1024, 6.0, 90.0, 600.0))
        .expect("narrow band");
    let gb = print_transfer("box 102.4, band 8", &base);
    let gl = print_transfer("box 204.8, band 8", &long);
    let gn = print_transfer("box 102.4, band 6", &narrow);
    for (k, w) in STROUHALS.iter().enumerate() {
        let box_rel = (gb[k] - gl[k]).abs() / gl[k].abs();
        let band_rel = (gb[k] - gn[k]).abs() / gb[k].abs();
        println!("w {w}: box rel {box_rel:.2e}, band rel {band_rel:.2e}");
        if *w >= 0.06 {
            assert!(
                box_rel < 0.01 && band_rel < 0.01,
                "kernel depends on the box at w = {w}"
            );
        }
    }
}

#[test]
fn hat_weights_are_causal_and_close_on_the_held_gain() {
    let s = bickley_step_response(STATION, &resolution(0.1, 1024, 8.0, 90.0, 600.0))
        .expect("step response");
    // A recorder-like jet: centre speed 22.5 m/s, half-width 0.4 mm, 44.1 kHz.
    let (speed, width, period) = (22.5, 0.4e-3, 1.0 / 44_100.0);
    let g = s.hat_weights(speed, width, period).expect("weights");
    let c = width / speed;
    let sum: f64 = g.iter().sum();
    println!(
        "{} taps over {:.3} ms; g0 = {:.4e} s (period/2 = {:.4e} s); sum = {:.6e} s = {:.6} b/U0; peak |g| = {:.4e} s",
        g.len(),
        s.window() * c * 1e3,
        g[0],
        0.5 * period,
        sum,
        sum / c,
        g.iter().fold(0.0f64, |m, v| m.max(v.abs()))
    );
    assert!(
        (sum / c - s.held_gain()).abs() <= 1e-9 * s.held_gain().abs(),
        "weights do not sum to the held gain"
    );
    assert!(
        (g[0] / (0.5 * period) - 1.0).abs() < 0.02,
        "current-sample weight is not the early s ~ t response"
    );
    assert!(
        g.last().expect("taps").abs() < 1e-12 * sum.abs(),
        "weights must end after the window"
    );
    // Nothing precedes the current sample by construction. SnowySummit's
    // admission check (ruling 7376) asks the same of the transfer: sampled on
    // [0, 10] b/U0 and inverted, band truncation rings before time zero for
    // any response, so the kernel's share there must not exceed that of the
    // probe's pure-delay controls (5, 10, 15 and 20 b) on the same grid. A delay that
    // lands on a sample (11 b is 35.01 samples here) barely rings and is not
    // a control.
    let (count, omega_max) = (2048usize, 10.0);
    let fft = RealFft::new(count);
    let before_zero = |sigma: f64, transfer: &dyn Fn(usize, f64) -> C64| {
        let margin = (5.0 * sigma * omega_max / core::f64::consts::PI).ceil() as usize;
        let spectrum: Vec<Bin> = (0..=count / 2)
            .map(|k| {
                let w = omega_max * k as f64 / (count / 2) as f64;
                let g = transfer(k, w);
                let taper = det::exp(-0.5 * sigma * sigma * w * w);
                Bin::new(g.re * taper, g.im * taper)
            })
            .collect();
        let impulse = fft.inverse(&spectrum);
        let before = &impulse[count / 2..count - margin];
        let energy =
            before.iter().map(|v| v * v).sum::<f64>() / impulse.iter().map(|v| v * v).sum::<f64>();
        let l1 = before.iter().map(|v| v.abs()).sum::<f64>()
            / impulse.iter().map(|v| v.abs()).sum::<f64>();
        (energy, l1)
    };
    for sigma in [0.0, 0.25] {
        let (ek, lk) = before_zero(sigma, &|k, w| {
            if k == 0 {
                C64::new(s.held_gain(), 0.0)
            } else {
                s.transfer(w).expect("transfer")
            }
        });
        let (mut ed, mut ld) = (f64::INFINITY, f64::INFINITY);
        for length in [5.0, 10.0, 15.0, 20.0] {
            let (e, l) = before_zero(sigma, &|_, w| {
                let (sn, cs) = det::sin_cos(w * length);
                C64::new(cs, -sn)
            });
            ed = ed.min(e);
            ld = ld.min(l);
        }
        println!(
            "taper {sigma}: before time zero, kernel energy {ek:.2e} L1 {lk:.2e}; least pure delay energy {ed:.2e} L1 {ld:.2e}"
        );
        assert!(
            ek <= ed && lk <= ld,
            "the kernel rings before time zero more than a pure delay (taper {sigma})"
        );
    }
}

#[test]
fn interval_weights_recompose_uniform_hats_and_refuse_invalid_clocks() {
    // Small numerical response; no measured jet or register claim.
    let response = bickley_step_response(3.0, &resolution(0.5, 32, 8.0, 8.0, 600.0)).unwrap();
    for rate in [48_000, 96_000, 192_000] {
        let dt = 1.0 / f64::from(rate);
        let advance = 12.0 * dt / 0.00075;
        let hats = response.hat_weights(12.0, 0.00075, dt).unwrap();
        for (k, expected) in hats.iter().enumerate() {
            let newer = response.interval_hat_weights(k as f64 * advance, advance, dt).unwrap().0;
            let older = if k == 0 { 0.0 } else {
                response.interval_hat_weights((k - 1) as f64 * advance, advance, dt).unwrap().1
            };
            assert!((newer + older - expected).abs() <= 1e-10 * expected.abs().max(1e-6),
                "intervals do not recompose hat {k} at {rate} Hz");
        }
    }
    assert_eq!(response.interval_hat_weights(response.window(), 1.0, 1e-5).unwrap(), (0.0, 0.0));
    for (age, advance, duration) in [
        (-1.0, 1.0, 1e-5), (f64::NAN, 1.0, 1e-5),
        (0.0, 0.0, 1e-5), (0.0, f64::INFINITY, 1e-5),
        (0.0, 1.0, 0.0), (0.0, 1.0, f64::NAN),
    ] {
        assert!(response.interval_hat_weights(age, advance, duration).is_err());
    }
}
