//! G1 disclosure through public PianoAudio and board-coordinate observations.
//! The authored two-coordinate fixture isolates integration error. It is not a
//! piano fidelity render, a measured specimen, or evidence for missing modes.
use fs_live_piano::refmods::{
    board_geometry::SurfaceSample, geometry::Course, linear::BoardMode, steinway_scale,
};
use fs_live_piano::{Config, Description, PianoAudio};
use fs_math::det;

const OUTPUT_RATE: u32 = 48_000;
const BOARD_HZ: f64 = 5_000.0;
const BRIDGE: f64 = 0.4;

fn fixture() -> (Description, Config, [f64; 2]) {
    let source = steinway_scale::courses()
        .unwrap()
        .into_iter()
        .find(|c| c.midi == 69)
        .unwrap();
    let course = Course {
        unison: 1,
        duplex_length_m: 0.0,
        detune_cents: 0.0,
        ..source
    };
    let mut bridge = [0.0; 88];
    bridge[48] = BRIDGE;
    let board = BoardMode {
        frequency_hz: BOARD_HZ,
        damping_ratio: 0.0,
        bridge,
        volume: 0.02,
    };
    // Independent continuous pencil in [string modal q, physical bare board b].
    // H = w_s²(q-beta*g*b)²/2 + (w_b²+T*g²/L)b²/2;
    // M = diag(1, 1+(mu*L/3-beta²)*g²). No fitted frequency enters it.
    let omega = core::f64::consts::TAU * course.partial_hz(1, course.tension_n);
    let beta =
        det::sqrt(2.0 * course.linear_density_kg_m * course.length_m) / core::f64::consts::PI;
    let mass =
        1.0 + (course.linear_density_kg_m * course.length_m / 3.0 - beta * beta) * BRIDGE * BRIDGE;
    let a = omega * omega;
    let d = ((core::f64::consts::TAU * BOARD_HZ).powi(2)
        + (course.tension_n / course.length_m + a * beta * beta) * BRIDGE * BRIDGE)
        / mass;
    let off = -a * beta * BRIDGE / det::sqrt(mass);
    let gap = det::sqrt((a - d) * (a - d) + 4.0 * off * off);
    let exact = [
        det::sqrt(0.5 * (a + d - gap)),
        det::sqrt(0.5 * (a + d + gap)),
    ];
    let desc = Description {
        courses: vec![course],
        board: vec![board],
        damping: None,
        surface: vec![SurfaceSample {
            position_m: [0.0; 3],
            area_m2: 0.02,
            mode_shape: vec![1.0],
        }],
        board_provenance: "authored two-coordinate integrator audit; source-backed string scale"
            .into(),
        board_mass_kg: 1.0,
        board_area_m2: 0.02,
        exact_low_modes: 1,
    };
    let mut cfg = Config::live(OUTPUT_RATE);
    cfg.max_partials = 1;
    cfg.course_mode_budget = 1;
    cfg.max_dynamic_modes = 1;
    cfg.partial_band_hz = 12_000.0;
    cfg.board_band_hz = 10_000.0;
    cfg.duplex = false;
    cfg.retire_ringing_j = 1e-30;
    cfg.retire_damped_j = 1e-30;
    (desc, cfg, exact)
}

fn phase_from_observed_board(samples: &[f64], low_omega: f64) -> f64 {
    // Twice apply the analytic low-mode annihilator, solely as an observation
    // filter. No sample or force fed back to the physical engine is modified.
    // Two applications suppress leakage from its small damped/warped offset.
    let cosine = det::cos(low_omega / f64::from(OUTPUT_RATE));
    let filter = |input: &[f64]| {
        input
            .windows(3)
            .map(|x| x[2] - 2.0 * cosine * x[1] + x[0])
            .collect::<Vec<_>>()
    };
    let isolated = filter(&filter(samples));
    let mut gram = [0.0; 3];
    let mut rhs = [0.0; 2];
    for x in isolated.windows(3) {
        gram[0] += x[1] * x[1];
        gram[1] += x[1] * x[0];
        gram[2] += x[0] * x[0];
        rhs[0] += x[1] * x[2];
        rhs[1] += x[0] * x[2];
    }
    let determinant = gram[0] * gram[2] - gram[1] * gram[1];
    assert!(
        determinant.is_finite() && determinant > 0.0,
        "board mode must be observable"
    );
    let a = (rhs[0] * gram[2] - rhs[1] * gram[1]) / determinant;
    let b = (rhs[1] * gram[0] - rhs[0] * gram[1]) / determinant;
    let real = 0.5 * a;
    assert!(
        b < 0.0 && -b > real * real,
        "observed mode must have a complex pole pair"
    );
    det::atan2(det::sqrt(-b - real * real), real).abs()
}

#[test]
fn coupled_board_modes_per_factor() {
    let git = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let host = std::process::Command::new("uname")
        .arg("-n")
        .output()
        .unwrap();
    assert!(git.status.success() && host.status.success());
    let git = String::from_utf8(git.stdout).unwrap();
    let host = String::from_utf8(host.stdout).unwrap();
    for value in [git.trim(), host.trim()] {
        assert!(
            !value.is_empty()
                && value
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b".-_".contains(&b))
        );
    }
    for factor in [1, 2, 4] {
        let (desc, cfg, exact) = fixture();
        let mut audio = PianoAudio::new(&desc, cfg, factor).unwrap();
        audio.note_on(69, 1.0).unwrap();
        let mut samples = Vec::with_capacity(12_288);
        let mut coordinate = [0.0];
        for n in 0..14_336 {
            audio.step_pressure().unwrap();
            audio.piano().board_coordinates(&mut coordinate);
            if n >= 2_048 {
                samples.push(coordinate[0]);
            }
        }
        assert_eq!(audio.mechanical_rate(), OUTPUT_RATE * factor as u32);
        assert_eq!(
            audio.piano().dynamic_modes(),
            1,
            "audit retained exactly its declared string coordinate"
        );
        assert!(
            samples.iter().any(|x| *x != 0.0),
            "actual coupled board motion required"
        );
        let a = phase_from_observed_board(&samples[..6_144], exact[0]) * f64::from(OUTPUT_RATE);
        let b = phase_from_observed_board(&samples[6_144..], exact[0]) * f64::from(OUTPUT_RATE);
        let measured = 0.5 * (a + b);
        let cents = 1200.0 * det::ln(measured / exact[1]) / core::f64::consts::LN_2;
        let window_difference = (1200.0 * det::ln(a / b) / core::f64::consts::LN_2).abs();
        assert!(measured.is_finite() && cents.is_finite());
        assert!(
            window_difference < 1.0,
            "high-mode estimate is not stationary: {window_difference} cents"
        );
        // This is the measured two-coordinate domain, not a universal band
        // for the full live piano. Include half the inter-window discrepancy
        // in its one-cent allowance; never relabel finite-step coupling exact.
        assert!(
            cents.abs() + 0.5 * window_difference <= 1.0,
            "factor {factor}: {cents} cents, window discrepancy {window_difference}; migrate the coupled solve"
        );
        println!(
            "{{\"path\":\"engine-piano::PianoAudio\",\"integrator\":\"exact local ZOH plus finite-step coupled solve\",\"factor\":{},\"omega_hz\":{},\"sample_rate_hz\":{},\"omega_dt\":{},\"window_a_omega_rad_s\":{},\"window_b_omega_rad_s\":{},\"measured_cents\":{},\"window_difference_cents\":{},\"predicted_cents\":null,\"tolerance_cents\":1,\"verdict\":\"measured-coupled-fixture\",\"git_sha\":\"{}\",\"host\":\"{}\"}}",
            factor,
            exact[1] / core::f64::consts::TAU,
            OUTPUT_RATE * factor as u32,
            exact[1] / f64::from(OUTPUT_RATE * factor as u32),
            a,
            b,
            cents,
            window_difference,
            git.trim(),
            host.trim()
        );
    }
}
