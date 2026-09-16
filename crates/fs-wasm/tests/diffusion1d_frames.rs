//! Tests for `diffusion1d_frames`. Bead: am-fs-export-diffusion1d-cew.
//!
//! The TypeScript owner (`src/physics/reference/diffusion/ftcs.ts`) implements
//! the same scheme. Field bits are required to match the documented unfused
//! update. If they disagree, the test names which engine is wrong.

use fs_wasm::{
    DIFFUSION1D_MAX_OUTPUT_LEN, DIFFUSION1D_MAX_TOTAL_STEPS, admit_diffusion1d_frames,
    assemble_zero_flux_laplacian, diffusion1d_frames, stability_ratio,
};

/// Prespecified before `convergence_against_gaussian` runs. Expected order is 2.
const MIN_OBSERVED_ORDER: f64 = 1.8;

/// Hand-verifiable golden from the TypeScript owner (`ftcs.test.mjs`):
/// n=5, frames=3, steps_per_frame=1, D=1, dx=1, dt=0.25, profile=0.
const TS_GOLDEN: [f64; 15] = [
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.25, 0.0, 0.0625, 0.25, 0.375, 0.25, 0.0625,
];

fn verdict(case: &str, outcome: &str) {
    eprintln!(
        "{{\"suite\":\"diffusion1d_frames\",\"beadId\":\"am-fs-export-diffusion1d-cew\",\"caseId\":\"{case}\",\"outcome\":\"{outcome}\"}}"
    );
}

fn fail_repro(case: &str, args: &str, message: &str) -> ! {
    eprintln!(
        "{{\"suite\":\"diffusion1d_frames\",\"beadId\":\"am-fs-export-diffusion1d-cew\",\"caseId\":\"{case}\",\"outcome\":\"fail\",\"message\":\"{message}\",\"repro\":\"cargo test -p fs-wasm --test diffusion1d_frames {case} -- --nocapture\",\"args\":{args}}}"
    );
    panic!("{case}: {message}");
}

#[test]
fn layout_and_profiles() {
    let spike = diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.25, 0).expect("valid spike");
    assert_eq!(spike.len(), 15);
    assert_eq!(&spike[..5], &[0.0, 0.0, 1.0, 0.0, 0.0]);
    let mass0: f64 = spike[..5].iter().sum::<f64>() * 1.0;
    assert!((mass0 - 1.0).abs() < 1e-15);

    let step = diffusion1d_frames(8, 1, 1, 1.0, 0.5, 0.01, 1).expect("valid step");
    assert_eq!(step.len(), 8);
    assert_eq!(&step[..4], &[1.0, 1.0, 1.0, 1.0]);
    assert_eq!(&step[4..], &[0.0, 0.0, 0.0, 0.0]);
    let mass_step: f64 = step.iter().sum::<f64>() * 0.5;
    assert!((mass_step - 4.0 * 0.5).abs() < 1e-15);

    let two = diffusion1d_frames(8, 1, 1, 1.0, 2.0, 0.01, 2).expect("valid two-spike");
    assert_eq!(two[8 / 4], 0.5 / 2.0);
    assert_eq!(two[(3 * 8) / 4], 0.5 / 2.0);
    let mass_two: f64 = two.iter().sum::<f64>() * 2.0;
    assert!((mass_two - 1.0).abs() < 1e-15);
    verdict("layout_and_profiles", "pass");
}

#[test]
fn operator_entries_and_column_sums() {
    let n = 7usize;
    let l = assemble_zero_flux_laplacian(n);
    assert_eq!(l.get(0, 0), -1.0);
    assert_eq!(l.get(0, 1), 1.0);
    for i in 1..n - 1 {
        assert_eq!(l.get(i, i - 1), 1.0);
        assert_eq!(l.get(i, i), -2.0);
        assert_eq!(l.get(i, i + 1), 1.0);
    }
    assert_eq!(l.get(n - 1, n - 2), 1.0);
    assert_eq!(l.get(n - 1, n - 1), -1.0);
    for c in 0..n {
        let mut sum = 0.0;
        for r in 0..n {
            sum += l.get(r, c);
        }
        assert!(
            sum.abs() <= 0.0,
            "column {c} sum {sum} must be exactly zero"
        );
    }
    verdict("operator_entries_and_column_sums", "pass");
}

#[test]
fn mass_conservation() {
    for profile in [0u32, 1, 2] {
        let n = 41usize;
        let dx = 0.2;
        let buf = diffusion1d_frames(n, 101, 5, 1.0, dx, 0.019, profile).expect("stable");
        let initial = &buf[..n];
        let mass: f64 = initial.iter().sum::<f64>() * dx;
        for frame in 1..101 {
            let sl = &buf[frame * n..(frame + 1) * n];
            let m: f64 = sl.iter().sum::<f64>() * dx;
            let rel = (m - mass).abs() / mass.max(1e-30);
            if rel > 1e-12 {
                fail_repro(
                    "mass_conservation",
                    &format!("{{\"profile\":{profile},\"frame\":{frame}}}"),
                    &format!("relative mass drift {rel}"),
                );
            }
        }
    }
    verdict("mass_conservation", "pass");
}

#[test]
fn maximum_principle() {
    for profile in [0u32, 1, 2] {
        let n = 41usize;
        let buf = diffusion1d_frames(n, 101, 5, 1.0, 0.2, 0.019, profile).expect("stable");
        let initial = &buf[..n];
        let max0 = initial.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min0 = initial.iter().copied().fold(f64::INFINITY, f64::min);
        for v in &buf {
            if *v < min0 - 1e-15 || *v > max0 + 1e-15 {
                fail_repro(
                    "maximum_principle",
                    &format!("{{\"profile\":{profile}}}"),
                    "new extremum",
                );
            }
        }
    }
    verdict("maximum_principle", "pass");
}

#[test]
fn second_moment_growth() {
    let n = 101usize;
    let dx = 0.1;
    let d = 1.0;
    let dt = 0.0025;
    let steps = 8usize;
    let buf = diffusion1d_frames(n, steps + 1, 1, d, dx, dt, 0).expect("stable");
    let ic = n / 2;
    let xc = (ic as f64 + 0.5) * dx;
    let m2 = |frame: usize| -> f64 {
        let sl = &buf[frame * n..(frame + 1) * n];
        assert!(
            sl[0].abs() < 1e-18 && sl[n - 1].abs() < 1e-18,
            "boundary contact"
        );
        sl.iter()
            .enumerate()
            .map(|(i, u)| {
                let x = (i as f64 + 0.5) * dx;
                u * (x - xc) * (x - xc) * dx
            })
            .sum()
    };
    for s in 0..steps {
        let delta = m2(s + 1) - m2(s);
        let expected = 2.0 * d * dt;
        if (delta - expected).abs() > 1e-12 {
            fail_repro(
                "second_moment_growth",
                &format!("{{\"step\":{s},\"delta\":{delta},\"expected\":{expected}}}"),
                "second moment did not grow by 2 D dt",
            );
        }
    }
    verdict("second_moment_growth", "pass");
}

#[test]
fn convergence_against_gaussian() {
    // Fixed physical time T = steps*dt. n doubles ⇒ dx halves ⇒ dt /4 ⇒ steps *4.
    let d = 1.0;
    let r = 0.25;
    let mut errors: Vec<f64> = Vec::new();
    let mut dxs: Vec<f64> = Vec::new();
    let mut t_used = 0.0_f64;
    for refine in 0..3 {
        let n = 200 * (1 << refine);
        let dx = 0.2 / (1 << refine) as f64;
        let dt = r * dx * dx / d;
        let steps = 20 * (1 << (2 * refine));
        t_used = steps as f64 * dt;
        let buf = diffusion1d_frames(n, 2, steps, d, dx, dt, 0).expect("stable");
        let field = &buf[n..];
        let ic = n / 2;
        let xc = (ic as f64 + 0.5) * dx;
        let denom = (4.0 * std::f64::consts::PI * d * t_used).sqrt();
        let mut max_err: f64 = 0.0;
        for (i, u) in field.iter().enumerate() {
            let x = (i as f64 + 0.5) * dx;
            let analytic = (-(x - xc) * (x - xc) / (4.0 * d * t_used)).exp() / denom;
            max_err = max_err.max((u - analytic).abs());
        }
        errors.push(max_err);
        dxs.push(dx);
    }
    let e0 = errors[0];
    let e1 = errors[1];
    let e2 = errors[2];
    let order01 = (e0 / e1).log2();
    let order12 = (e1 / e2).log2();
    eprintln!(
        "{{\"suite\":\"diffusion1d_frames\",\"caseId\":\"convergence_against_gaussian\",\"T\":{t_used},\"dx\":[{},{},{}],\"errors\":[{e0},{e1},{e2}],\"observedOrder01\":{order01},\"observedOrder12\":{order12},\"threshold\":{MIN_OBSERVED_ORDER}}}",
        dxs[0], dxs[1], dxs[2]
    );
    if order01 < MIN_OBSERVED_ORDER || order12 < MIN_OBSERVED_ORDER {
        fail_repro(
            "convergence_against_gaussian",
            &format!(
                "{{\"observedOrder01\":{order01},\"observedOrder12\":{order12},\"threshold\":{MIN_OBSERVED_ORDER}}}"
            ),
            "observed order below prespecified 1.8",
        );
    }
    verdict("convergence_against_gaussian", "pass");
}

#[test]
fn stability_boundary() {
    diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.5, 0).expect("r = 0.5 is admissible");
    let err = diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.5000001, 0).expect_err("must refuse");
    assert_eq!(err.code, "ftcs-unstable");
    assert!(err.details.contains("\"limit\":0.5"));
    assert!(err.details.contains("\"dtMax\":"));
    assert!(err.details.contains("0.5000001") || err.details.contains("ratio"));
    let r = stability_ratio(1.0, 0.5000001, 1.0);
    assert!(r > 0.5);
    verdict("stability_boundary", "pass");
}

#[test]
fn ratio_expression_order() {
    let d = 0.1;
    let dx = 0.1;
    let dt = 0.05;
    let documented = (d * dt) / (dx * dx);
    let chained = d * dt / dx / dx;
    let grouped = d * (dt / (dx * dx));
    assert_eq!(documented, 0.5);
    assert!(
        chained > 0.5,
        "planted chained order must refuse: {chained}"
    );
    assert!(grouped <= 0.5, "grouped order {grouped}");
    assert_eq!(stability_ratio(d, dt, dx), documented);
    diffusion1d_frames(5, 2, 1, d, dx, dt, 0).expect("documented order accepts exactly 0.5");
    verdict("ratio_expression_order", "pass");
}

#[test]
fn unfused_update_reconstruction() {
    let n = 5usize;
    let r = 0.25;
    let mut u = vec![0.0, 0.0, 1.0, 0.0, 0.0];
    let mut y = vec![0.0; n];
    y[0] = -u[0] + u[1];
    for i in 1..n - 1 {
        let mut sum = u[i - 1];
        sum += -2.0 * u[i];
        sum += u[i + 1];
        y[i] = sum;
    }
    y[n - 1] = u[n - 2] - u[n - 1];
    for i in 0..n {
        let increment = r * y[i];
        u[i] = u[i] + increment;
    }
    let exported = diffusion1d_frames(n, 2, 1, 1.0, 1.0, 0.25, 0).expect("valid");
    assert_eq!(&exported[n..], u.as_slice());
    for (i, (a, b)) in exported[n..].iter().zip(u.iter()).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "cell {i} not bitwise equal");
    }
    verdict("unfused_update_reconstruction", "pass");
}

#[test]
fn ts_golden_frames_bitwise() {
    let buf = diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.25, 0).expect("valid");
    assert_eq!(buf.len(), TS_GOLDEN.len());
    for (i, (a, b)) in buf.iter().zip(TS_GOLDEN.iter()).enumerate() {
        if a.to_bits() != b.to_bits() {
            fail_repro(
                "ts_golden_frames_bitwise",
                &format!("{{\"index\":{i},\"rust\":{a:?},\"typescript\":{b:?}}}"),
                "Rust field disagrees with the TypeScript golden; do not pick a winner silently",
            );
        }
    }
    verdict("ts_golden_frames_bitwise", "pass");
}

#[test]
fn refusal_payloads() {
    let cases: &[(&str, Result<Vec<f64>, _>)] = &[
        ("n=2", diffusion1d_frames(2, 2, 1, 1.0, 1.0, 0.1, 0)),
        ("frames=0", diffusion1d_frames(5, 0, 1, 1.0, 1.0, 0.1, 0)),
        ("steps=0", diffusion1d_frames(5, 2, 0, 1.0, 1.0, 0.1, 0)),
        ("D<0", diffusion1d_frames(5, 2, 1, -1.0, 1.0, 0.1, 0)),
        ("dt=0", diffusion1d_frames(5, 2, 1, 1.0, 1.0, 0.0, 0)),
        ("dx=0", diffusion1d_frames(5, 2, 1, 1.0, 0.0, 0.1, 0)),
        ("profile=3", diffusion1d_frames(5, 2, 1, 1.0, 1.0, 0.1, 3)),
        ("NaN D", diffusion1d_frames(5, 2, 1, f64::NAN, 1.0, 0.1, 0)),
    ];
    for (name, result) in cases {
        let err = result.as_ref().expect_err(name);
        assert!(
            err.code == "invalid-parameter"
                || err.code == "nonfinite-input"
                || err.code == "unsupported-kernel",
            "{name} code {}",
            err.code
        );
        assert!(!err.message.is_empty());
    }
    let huge = diffusion1d_frames(2000, 2000, 1, 1.0, 1.0, 0.1, 0).expect_err("budget");
    assert_eq!(huge.code, "budget-exhausted");
    assert!(DIFFUSION1D_MAX_OUTPUT_LEN < 2000 * 2000);
    let _ = DIFFUSION1D_MAX_TOTAL_STEPS;
    verdict("refusal_payloads", "pass");
}

#[test]
fn zero_diffusion() {
    let buf = diffusion1d_frames(5, 4, 7, 0.0, 1.0, 100.0, 0).expect("D=0 is valid");
    for frame in 1..4 {
        assert_eq!(&buf[frame * 5..(frame + 1) * 5], &buf[..5]);
    }
    verdict("zero_diffusion", "pass");
}

#[test]
fn golden_hash() {
    let buf = diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.25, 0).expect("valid");
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for v in &buf {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    }
    const GOLDEN_HASH: u64 = 0x7d40_609c_183c_41d0;
    if h != GOLDEN_HASH {
        eprintln!(
            "{{\"suite\":\"diffusion1d_frames\",\"caseId\":\"golden_hash\",\"actual\":\"{h:#x}\",\"expected\":\"{GOLDEN_HASH:#x}\",\"outcome\":\"fail\"}}"
        );
        panic!(
            "golden hash {h:#x} != {GOLDEN_HASH:#x}; re-pin only with a documented scheme change"
        );
    }
    verdict("golden_hash", "pass");
}

#[test]
fn admit_matches_run() {
    let spec = admit_diffusion1d_frames(5, 3, 1, 1.0, 1.0, 0.25, 0).expect("admit");
    assert_eq!(spec.n(), 5);
    assert_eq!(spec.frames(), 3);
    assert_eq!(spec.r(), 0.25);
    verdict("admit_matches_run", "pass");
}
