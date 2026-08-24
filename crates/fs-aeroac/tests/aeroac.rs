//! fs-aeroac conformance battery: Bessel/Hankel certification (exact
//! Wronskian, independent-path derivative cross-checks, fsci-special
//! cephes-heritage oracle), 2D Curle dipole physics (cylindrical
//! spreading exponent, exact directivity null, the 3D-Green's-swap
//! mutation), the self-verified Bickley analytic pins, the Rayleigh
//! shooting solver against them, and typed refusals.

use fs_aeroac::bessel::{hankel0_outgoing, j0, j1, y0, y1};
use fs_aeroac::bickley::{JetSymmetry, bickley_rayleigh_mode, rayleigh_residual_closed_form};
use fs_aeroac::curle2d::{
    dipole_pressure, modulate_observer_by_tone, strouhal_at_reynolds, tonal_dipole_observer,
    tonal_lock_frequency,
};
use fs_aeroac::{AeroacError, SCOPE_STATEMENT};
use fs_math::c64::C64;
use fs_math::det;

/// Log grid spanning series regime, crossover, and asymptotic regime.
fn log_grid() -> Vec<f64> {
    let mut xs = Vec::new();
    let mut x = 1.0e-3;
    while x < 3.0e3 {
        xs.push(x);
        x *= 1.17;
    }
    // The crossover itself, both sides.
    xs.extend_from_slice(&[17.9, 17.999, 18.0, 18.001, 18.1]);
    xs
}

/// The EXACT Wronskian identity `J1 Y0 - J0 Y1 = 2/(pi x)`: one
/// equation certifying all four functions jointly (an error in any
/// regime or at the crossover breaks it).
#[test]
fn bessel_wronskian_identity() {
    let mut worst = 0.0f64;
    for x in log_grid() {
        let w = j1(x) * y0(x) - j0(x) * y1(x);
        let err = (w * core::f64::consts::PI * x / 2.0 - 1.0).abs();
        worst = worst.max(err);
        assert!(err < 5.0e-13, "Wronskian at x={x}: err={err:.3e}");
    }
    println!("{{\"suite\":\"fs-aeroac\",\"case\":\"wronskian\",\"worst\":{worst:.3e}}}");
}

/// Derivative identities via central finite differences —
/// INDEPENDENT code paths (J0 vs J1, Y0 vs Y1): `J0' = -J1`,
/// `Y0' = -Y1`. FD truncation limits the gate, not the functions.
#[test]
fn bessel_derivative_cross_checks() {
    let h = 1.0e-5;
    for x in [0.3, 1.0, 4.0, 9.0, 17.5, 18.5, 60.0, 400.0] {
        let dj0 = (j0(x + h) - j0(x - h)) / (2.0 * h);
        assert!(
            (dj0 + j1(x)).abs() < 1.0e-8,
            "J0' = -J1 at x={x}: {dj0} vs {}",
            -j1(x)
        );
        let dy0 = (y0(x + h) - y0(x - h)) / (2.0 * h);
        assert!(
            (dy0 + y1(x)).abs() < 1.0e-8,
            "Y0' = -Y1 at x={x}: {dy0} vs {}",
            -y1(x)
        );
    }
}

/// Cross-implementation oracle: fsci-special's cephes-heritage
/// j0/j1/y0/y1 (dev-dependency, BSD lineage) over the log grid.
/// Gate at 5e-12 relative-to-envelope: fsci's Y-series accumulate in
/// PLAIN f64 and carry ~1e-12 cancellation near their own x=14
/// crossover (measured: y1(12.34) differs 9.5e-13), while this
/// crate's double-double series hold ~1e-15 — so the oracle is a
/// gross-error cross-check and the EXACT Wronskian identity test is
/// the precision arbiter (5e-13 there).
#[test]
fn bessel_matches_fsci_special_oracle() {
    use fsci_runtime::RuntimeMode;
    use fsci_special::SpecialTensor;
    let eval = |f: fn(&SpecialTensor, RuntimeMode) -> fsci_special::SpecialResult, x: f64| -> f64 {
        match f(&SpecialTensor::RealScalar(x), RuntimeMode::Strict).expect("fsci") {
            SpecialTensor::RealScalar(v) => v,
            other => panic!("expected scalar, got {other:?}"),
        }
    };
    let mut worst = 0.0f64;
    for x in log_grid() {
        let envelope = det::sqrt(2.0 / (core::f64::consts::PI * x)).max(1.0);
        for (name, mine, theirs) in [
            ("j0", j0(x), eval(fsci_special::j0, x)),
            ("j1", j1(x), eval(fsci_special::j1, x)),
            ("y0", y0(x), eval(fsci_special::y0, x)),
            ("y1", y1(x), eval(fsci_special::y1, x)),
        ] {
            let err = (mine - theirs).abs() / envelope.max(theirs.abs());
            worst = worst.max(err);
            assert!(err < 5.0e-12, "{name}({x}): {mine:.17} vs {theirs:.17}");
        }
    }
    println!("{{\"suite\":\"fs-aeroac\",\"case\":\"fsci-oracle\",\"worst\":{worst:.3e}}}");
}

/// Small-argument limits with exact leading behavior.
#[test]
fn bessel_small_argument_limits() {
    assert!((j0(0.0) - 1.0).abs() < 1.0e-15);
    assert!((j1(1.0e-8) - 0.5e-8).abs() < 1.0e-20);
    // Y0 log singularity: Y0(x) - (2/pi) ln(x) bounded as x -> 0.
    let a = y0(1.0e-6) - core::f64::consts::FRAC_2_PI * det::ln(1.0e-6);
    let b = y0(1.0e-8) - core::f64::consts::FRAC_2_PI * det::ln(1.0e-8);
    assert!((a - b).abs() < 1.0e-9, "Y0 log constant drifts: {a} vs {b}");
    // Domain: Y at non-positive arguments is NaN, not a fabrication.
    assert!(y0(0.0).is_nan() && y0(-1.0).is_nan() && y1(-2.0).is_nan());
    // Hankel domain policy: BOTH components NaN outside the domain
    // (review-caught: a half-valid C64(1.0, NaN) once leaked).
    let h = hankel0_outgoing(-1.0);
    assert!(h.re.is_nan() && h.im.is_nan());
    let h = fs_aeroac::bessel::hankel1_outgoing(0.0);
    assert!(h.re.is_nan() && h.im.is_nan());
    // Odd symmetry including signed zero.
    assert!(j1(-0.0).is_sign_negative());
}

/// Far-field amplitude of the outgoing Hankel function:
/// `|H0(x)| -> sqrt(2/(pi x))`.
#[test]
fn hankel_far_field_amplitude() {
    for x in [50.0, 300.0, 2000.0] {
        let h = hankel0_outgoing(x);
        let expect = det::sqrt(2.0 / (core::f64::consts::PI * x));
        // The envelope correction is O(1/(16 x^2)) (from P^2 + Q^2
        // = 1 + O(x^-2)); gate at ~1.6x that.
        assert!(
            (h.abs() / expect - 1.0).abs() < 0.1 / (x * x) + 1.0e-12,
            "|H0({x})| = {} vs {expect}",
            h.abs()
        );
    }
}

/// 2D Curle dipole physics: cylindrical spreading (log-slope -1/2 —
/// the 3D-Green's-swap mutation would read -1), the EXACT directivity
/// null perpendicular to the force, cos-theta directivity shape, and
/// the embedded honest-scope statement.
#[test]
// float_cmp: the perpendicular dipole null is EXACT by construction
// (rhat . F = 0 in exact arithmetic), not a computed near-zero.
#[allow(clippy::float_cmp)]
fn curle_dipole_spreading_and_directivity() {
    let k = 30.0;
    let fx = [C64::new(1.0, 0.0), C64::new(0.0, 0.0)];
    // Far-field decay exponent between r = 10 and r = 20 (kr = 300+).
    let p1 = dipole_pressure(fx, k, [10.0, 0.0], [0.0, 0.0]).expect("p1");
    let p2 = dipole_pressure(fx, k, [20.0, 0.0], [0.0, 0.0]).expect("p2");
    let slope = det::ln(p2.pressure.abs() / p1.pressure.abs()) / det::ln(2.0f64);
    assert!(
        (slope + 0.5).abs() < 0.01,
        "cylindrical spreading exponent: {slope} (a 3D Green's function reads -1)"
    );
    // Exact null perpendicular to the force (rhat . F = 0 exactly).
    let perp = dipole_pressure(fx, k, [0.0, 7.0], [0.0, 0.0]).expect("perp");
    assert_eq!(perp.pressure.abs(), 0.0, "dipole null must be exact");
    // cos(theta) directivity at fixed radius (45 degrees).
    let diag = dipole_pressure(
        fx,
        k,
        [10.0 / det::sqrt(2.0), 10.0 / det::sqrt(2.0)],
        [0.0, 0.0],
    )
    .expect("diag");
    let ratio = diag.pressure.abs() / p1.pressure.abs();
    assert!(
        (ratio - 1.0 / det::sqrt(2.0)).abs() < 1.0e-3,
        "cos-theta directivity: {ratio}"
    );
    // The honest-scope statement rides on every output and is the
    // crate constant (the marketing-mutation guard).
    assert_eq!(p1.scope, SCOPE_STATEMENT);
    assert!(SCOPE_STATEMENT.contains("NOT absolute SPL"));
    assert!(SCOPE_STATEMENT.contains("2D-to-3D span correction"));
}

/// Curle refusals, typed by name.
#[test]
fn curle_refusals_are_typed() {
    let f = [C64::new(1.0, 0.0), C64::new(0.0, 0.0)];
    assert!(matches!(
        dipole_pressure(f, 0.0, [1.0, 0.0], [0.0, 0.0]),
        Err(AeroacError::InvalidParameter { .. })
    ));
    assert!(matches!(
        dipole_pressure(f, -3.0, [1.0, 0.0], [0.0, 0.0]),
        Err(AeroacError::InvalidParameter { .. })
    ));
    assert!(matches!(
        dipole_pressure(f, 5.0, [2.0, 3.0], [2.0, 3.0]),
        Err(AeroacError::InvalidParameter { .. })
    ));
    assert!(matches!(
        dipole_pressure(f, f64::NAN, [1.0, 0.0], [0.0, 0.0]),
        Err(AeroacError::NonFinite { .. })
    ));
    let bad = [C64::new(f64::INFINITY, 0.0), C64::new(0.0, 0.0)];
    assert!(matches!(
        dipole_pressure(bad, 5.0, [1.0, 0.0], [0.0, 0.0]),
        Err(AeroacError::NonFinite { .. })
    ));
}

#[test]
fn pinned_strouhal_ladder_is_interpolated_not_invented() {
    let st_lo = strouhal_at_reynolds(432.0).expect("432");
    let st_hi = strouhal_at_reynolds(2_304.0).expect("2304");
    assert!((st_lo - 0.092).abs() < 1.0e-12);
    assert!((st_hi - 0.467).abs() < 1.0e-12);
    let mid = strouhal_at_reynolds(500.0).expect("mid");
    assert!(mid > st_lo && mid < 0.101);
    assert!(strouhal_at_reynolds(0.0).is_none());
    let f = tonal_lock_frequency(432.0, 20.0, 0.006).expect("f");
    assert!((f - 0.092 * 20.0 / 0.006).abs() < 1.0e-12);
}

#[test]
fn tonal_dipole_is_observer_side_and_nulls_off_axis() {
    let on = tonal_dipole_observer(1.2, 20.0, 0.006, 432.0, 0.4, 343.0, [2.0, 0.0], [0.0, 0.0])
        .expect("on");
    let off = tonal_dipole_observer(1.2, 20.0, 0.006, 432.0, 0.4, 343.0, [0.0, 2.0], [0.0, 0.0])
        .expect("off");
    assert!(on.pressure.abs() > 0.0);
    assert_eq!(off.pressure.abs(), 0.0);
    assert_eq!(on.scope, SCOPE_STATEMENT);
}

#[test]
fn observer_tone_modulation_does_not_invent_broadband() {
    let n = 512;
    let dt = 1.0 / 8_192.0;
    let mut p = vec![1.0; n];
    modulate_observer_by_tone(&mut p, dt, 256.0, 0.5);
    // Mean stays 1 (zero-mean sine). Peak-to-peak is the depth.
    let max = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let min = p.iter().copied().fold(f64::INFINITY, f64::min);
    assert!((max - 1.5).abs() < 1.0e-9);
    assert!((min - 0.5).abs() < 1.0e-9);
}

/// The analytic pins, SELF-VERIFIED: `phi = sech^2` at
/// `(alpha, c) = (2, 2/3)` and `phi = sech tanh` at `(1, 2/3)` make
/// the Rayleigh residual machine-zero at every probe point — and the
/// falsifier half proves the residual is LIVE (perturbing alpha or c
/// moves it far from zero).
#[test]
fn bickley_analytic_pins_self_verify() {
    let probes = [0.15, 0.5, 0.9, 1.4, 2.2, 3.5];
    for &y in &probes {
        let r_sin = rayleigh_residual_closed_form(JetSymmetry::Sinuous, 2.0, 2.0 / 3.0, y);
        assert!(
            r_sin.abs() < 1.0e-14,
            "sinuous pin residual at y={y}: {r_sin:e}"
        );
        let r_var = rayleigh_residual_closed_form(JetSymmetry::Varicose, 1.0, 2.0 / 3.0, y);
        assert!(
            r_var.abs() < 1.0e-14,
            "varicose pin residual at y={y}: {r_var:e}"
        );
    }
    // Falsifiers: the residual detects wrong eigenpairs.
    let off_alpha = rayleigh_residual_closed_form(JetSymmetry::Sinuous, 2.05, 2.0 / 3.0, 0.9);
    let off_c = rayleigh_residual_closed_form(JetSymmetry::Sinuous, 2.0, 0.68, 0.9);
    assert!(
        off_alpha.abs() > 1.0e-3,
        "alpha falsifier inert: {off_alpha:e}"
    );
    assert!(off_c.abs() > 1.0e-3, "c falsifier inert: {off_c:e}");
}

/// The shooting solver against the self-verified pins and physical
/// structure: near the sinuous neutral point `c -> 2/3` with
/// vanishing growth; the jet is UNSTABLE (positive growth) across the
/// sinuous band; the eigenvalue is grid-converged.
#[test]
fn bickley_shooting_solver_matches_pins() {
    // Near-neutral sinuous: alpha = 1.95 must sit close to the exact
    // neutral point (alpha = 2, c = 2/3), approached from the
    // unstable side (Im c > 0 keeps the critical layer off the real
    // axis; AT the neutral point the real-axis singularity makes
    // shooting ill-posed, which is why the pin is checked NEARBY).
    let near = bickley_rayleigh_mode(1.95, JetSymmetry::Sinuous, C64::new(0.66, 0.02), 14.0, 2048)
        .expect("near-neutral sinuous");
    assert!(
        (near.c.re - 2.0 / 3.0).abs() < 0.02 && near.c.im.abs() < 0.02 && near.c.im > 0.0,
        "sinuous near-neutral c = {:?}",
        near.c
    );
    // Unstable sinuous band: growth strictly positive and larger at
    // mid-band than near the neutral point.
    let mid = bickley_rayleigh_mode(1.0, JetSymmetry::Sinuous, C64::new(0.6, 0.15), 14.0, 2048)
        .expect("mid-band sinuous");
    assert!(
        mid.growth_rate > 0.05,
        "mid-band growth: {}",
        mid.growth_rate
    );
    assert!(
        mid.growth_rate > near.growth_rate,
        "growth must fall toward the neutral point: {} vs {}",
        mid.growth_rate,
        near.growth_rate
    );
    // Varicose near ITS exact neutral point (alpha = 1, c = 2/3):
    // approach from below in alpha where the mode is unstable.
    let vnear = bickley_rayleigh_mode(
        0.95,
        JetSymmetry::Varicose,
        C64::new(0.66, 0.02),
        14.0,
        2048,
    )
    .expect("near-neutral varicose");
    assert!(
        (vnear.c.re - 2.0 / 3.0).abs() < 0.03 && vnear.c.im.abs() < 0.03 && vnear.c.im > 0.0,
        "varicose near-neutral c = {:?}",
        vnear.c
    );
    // Grid convergence of the eigenvalue (2048 vs 4096 steps).
    let fine = bickley_rayleigh_mode(1.0, JetSymmetry::Sinuous, C64::new(0.6, 0.15), 14.0, 4096)
        .expect("fine");
    let dc = (fine.c - mid.c).abs();
    assert!(dc < 1.0e-8, "eigenvalue grid drift: {dc:e}");
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"bickley\",\"sinuous_mid_c\":[{},{}],\"near_neutral_c\":[{},{}],\"verdict\":\"pass\"}}",
        mid.c.re, mid.c.im, near.c.re, near.c.im
    );
}

/// Rayleigh solver refusals, typed by name.
#[test]
fn bickley_refusals_are_typed() {
    assert!(matches!(
        bickley_rayleigh_mode(0.0, JetSymmetry::Sinuous, C64::new(0.6, 0.1), 14.0, 2048),
        Err(AeroacError::InvalidParameter { .. })
    ));
    assert!(matches!(
        bickley_rayleigh_mode(
            f64::NAN,
            JetSymmetry::Sinuous,
            C64::new(0.6, 0.1),
            14.0,
            2048
        ),
        Err(AeroacError::NonFinite { .. })
    ));
    // Finite inputs can still overflow the RK4 state. A non-finite
    // boundary mismatch must refuse rather than slipping through the
    // `mismatch > tolerance` check (which is false for NaN).
    assert!(matches!(
        bickley_rayleigh_mode(f64::MAX, JetSymmetry::Sinuous, C64::new(0.6, 0.1), 14.0, 64),
        Err(AeroacError::NonFinite { .. })
    ));
    assert!(matches!(
        bickley_rayleigh_mode(1.0, JetSymmetry::Sinuous, C64::new(0.6, 0.1), 1.0, 2048),
        Err(AeroacError::InvalidParameter { .. })
    ));
    // A hopeless guess must refuse, not return a fabricated mode.
    assert!(matches!(
        bickley_rayleigh_mode(1.0, JetSymmetry::Sinuous, C64::new(50.0, 40.0), 14.0, 256),
        Err(AeroacError::NotConverged { .. } | AeroacError::NonFinite { .. })
    ));
}

/// Determinism: bitwise-identical reruns (table stakes for the
/// workspace).
#[test]
fn aeroac_determinism_bitwise() {
    let a = bickley_rayleigh_mode(1.2, JetSymmetry::Sinuous, C64::new(0.6, 0.12), 14.0, 1024)
        .expect("a");
    let b = bickley_rayleigh_mode(1.2, JetSymmetry::Sinuous, C64::new(0.6, 0.12), 14.0, 1024)
        .expect("b");
    assert_eq!(a.c.re.to_bits(), b.c.re.to_bits());
    assert_eq!(a.c.im.to_bits(), b.c.im.to_bits());
    let p = dipole_pressure(
        [C64::new(1.0, 0.5), C64::new(-0.25, 0.0)],
        12.0,
        [3.0, 4.0],
        [0.0, 0.0],
    )
    .expect("p");
    let q = dipole_pressure(
        [C64::new(1.0, 0.5), C64::new(-0.25, 0.0)],
        12.0,
        [3.0, 4.0],
        [0.0, 0.0],
    )
    .expect("q");
    assert_eq!(p.pressure.re.to_bits(), q.pressure.re.to_bits());
    assert_eq!(p.pressure.im.to_bits(), q.pressure.im.to_bits());
}
