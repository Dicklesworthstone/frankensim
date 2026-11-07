use super::super::{
    DgBoundary, DgDuct, DgMesh, DgPropagator, PortId, PropagatorKernel, scale_matrix, series_phis,
    series_shift,
};
use super::{REFINEMENTS, SERIES_DD_MAX_DOFS, engages, phis, refine};
use crate::{Duct, DuctError, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;

const DT: f64 = 1.0 / 96_000.0;

fn refinements() -> usize {
    REFINEMENTS.with(std::cell::Cell::get)
}

fn reset_refinements() {
    REFINEMENTS.with(|count| count.set(0));
}

/// One chimney as the music_render DG route builds it: a cylinder, one element of order 1,
/// ports on both ends, the parent's wall-memory band.
fn chimney(radius_m: f64, length_m: f64, poles: usize, band_hz: [f64; 2], error: f64) -> DgDuct {
    let source = Duct {
        segments: vec![Segment::Cylinder { radius: radius_m, length: length_m }],
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 305.0, 101_325.0).unwrap();
    DgDuct::new(
        &source,
        &gas,
        DgMesh { elements: 1, order: 1 },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap()
    .with_viscothermal(poles, band_hz, error)
    .unwrap()
}

/// The oboe card's chimney: 60 poles, 50 Hz to 20 kHz, error 0.01 (244 DOFs).
fn card_chimney(radius_m: f64) -> DgDuct {
    chimney(radius_m, 0.006, 60, [50.0, 20_000.0], 0.01)
}

/// A small chimney (36 DOFs: 8 poles over 50 Hz to 2 kHz, the pole count the port-response tests
/// use) with a step `dt` doubled from the 96 kHz clock until it takes at least 13 squarings, as
/// the oboe's refusing chimney does. The fallback's logic in milliseconds, in debug builds too.
fn small_case() -> (DgDuct, f64) {
    let duct = chimney(0.0015, 0.006, 8, [50.0, 2_000.0], 0.01);
    let n = duct.total_dofs();
    let mut dt = DT;
    for _ in 0..12 {
        if series_shift(&duct.l_tilde, n, dt) >= 13 {
            break;
        }
        dt *= 2.0;
    }
    (duct, dt)
}

struct F64Series {
    shift: u32,
    carried: f64,
    halving: f64,
    phi0: Vec<f64>,
}

impl F64Series {
    fn certificate(&self) -> f64 {
        self.carried.max(self.halving)
    }
}

/// The f64 path's own series and certificate terms at a given number of squarings.
fn f64_series(duct: &DgDuct, dt: f64, shift: u32) -> F64Series {
    let n = duct.total_dofs();
    let (phi0, _, _, truncation) = series_phis(&duct.l_tilde, n, dt, shift).unwrap();
    let (finer, _, _, _) = series_phis(&duct.l_tilde, n, dt, shift + 1).unwrap();
    let halving = det::sqrt(
        phi0.iter().zip(&finer).map(|(a, b)| (a - b) * (a - b)).sum::<f64>(),
    );
    let carried = f64::from(1u32 << shift) * det::sqrt(n as f64) * truncation;
    F64Series { shift, carried, halving, phi0 }
}

fn dense(propagator: &DgPropagator) -> (&[f64], &[f64]) {
    match &propagator.kernel {
        PropagatorKernel::DenseSeries { phi0, dt_phi1, .. } => (phi0, dt_phi1),
        _ => panic!("a dense series kernel is expected"),
    }
}

fn frobenius_difference(a: &[f64], b: &[f64]) -> f64 {
    det::sqrt(a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>())
}

#[test]
fn the_fallback_engages_only_on_a_finite_refusal_of_a_small_duct() {
    assert!(engages(244, 1.136e-11, 1e-11), "today's refusal engages");
    assert!(engages(SERIES_DD_MAX_DOFS, 1.136e-11, 1e-11));
    assert!(!engages(SERIES_DD_MAX_DOFS + 1, 1.136e-11, 1e-11), "above the DOF ceiling");
    assert!(!engages(244, 1e-11, 1e-11), "a passing certificate is not a refusal");
    assert!(!engages(244, 0.5e-11, 1e-11));
    assert!(!engages(244, f64::NAN, 1e-11), "a NaN certificate is a bad operator");
    assert!(!engages(244, f64::INFINITY, 1e-11));
}

#[test]
fn scalar_decay_in_double_double_matches_the_closed_forms_and_itself_at_more_squarings() {
    for z in [-0.25_f64, -0.5, -4.0, -32.0] {
        let exp = det::exp(z);
        let expected = [exp, (exp - 1.0) / z, (exp - 1.0 - z) / (z * z)];
        let shift = series_shift(&[z], 1, 1.0);
        let series = refine(&[z], 1, 1.0, shift).expect("a scalar decay is finite");
        for (order, (actual, reference)) in [series.phi0[0], series.phi1[0], series.phi2[0]]
            .into_iter()
            .zip(expected)
            .enumerate()
        {
            assert!(
                (actual - reference).abs() <= 2e-13 * reference.abs() + 1e-30,
                "phi_{order}({z}): {actual:e} != {reference:e}"
            );
        }
        assert!(series.certificate() < 1e-14, "z = {z}: certificate {:e}", series.certificate());
        // The same value with three more squarings, in double-double: equal to far below f64.
        let base = phis(&[z], 1, 1.0, shift, true);
        let more = phis(&[z], 1, 1.0, shift + 3, true);
        for (a, b) in [
            (base.phi0[0], more.phi0[0]),
            (base.phi1[0], more.phi1[0]),
            (base.phi2[0], more.phi2[0]),
        ] {
            let d = (a - b).abs().to_f64();
            assert!(d <= 1e-26 * a.abs().to_f64().max(1e-300), "z = {z}: {d:e}");
        }
    }
}

#[test]
fn a_small_chimney_refused_by_f64_is_certified_in_double_double() {
    let (duct, dt) = small_case();
    let n = duct.total_dofs();
    let shift = series_shift(&duct.l_tilde, n, dt);
    let f = f64_series(&duct, dt, shift);
    let f64_certificate = f.certificate();
    assert!(
        f64_certificate > 1e-13,
        "the f64 certificate {f64_certificate:e} (n {n}, shift {shift}) is not well above the double-double rounding floor; the refusal cannot be forced"
    );
    // A tolerance the f64 certificate misses and the double-double one meets.
    let tol = 0.5 * f64_certificate;
    reset_refinements();
    let propagator = duct.propagator(dt, tol).expect("the fallback certifies the chimney");
    assert_eq!(refinements(), 1, "the f64 certificate refused, so exactly one refinement ran");
    let certificate = propagator.certificate_residual();
    assert!(certificate <= tol, "certificate {certificate:e} > tol {tol:e}");
    assert!(
        certificate < 0.1 * f64_certificate,
        "certificate {certificate:e} against the f64 {f64_certificate:e}"
    );
    // The double-double result is stable under more squarings far below the tolerance, and the
    // f64 matrices differ from it by the order of the f64 certificate.
    let (phi0, _) = dense(&propagator);
    let more = phis(&duct.l_tilde, n, dt, shift + 2, false);
    let more_f64: Vec<f64> = more.phi0.iter().map(|v| v.to_f64()).collect();
    let stability = frobenius_difference(phi0, &more_f64);
    assert!(stability <= 4.0 * f64::EPSILON * (n as f64).sqrt(), "stability {stability:e}");
    let f64_error = frobenius_difference(&f.phi0, phi0);
    assert!(
        f64_error <= 10.0 * f64_certificate,
        "the f64 matrices are {f64_error:e} from the double-double ones; the f64 certificate says {f64_certificate:e}"
    );
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"series-dd-small\",\"n\":{n},\"shift\":{shift},\"f64_certificate\":{f64_certificate:e},\"dd_certificate\":{certificate:e},\"f64_phi0_vs_dd\":{f64_error:e},\"dd_stability\":{stability:e}}}"
    );
}

#[test]
fn a_duct_that_passes_in_f64_keeps_the_f64_bits_and_never_reaches_the_fallback() {
    let (duct, dt) = small_case();
    let n = duct.total_dofs();
    let shift = series_shift(&duct.l_tilde, n, dt);
    let f = f64_series(&duct, dt, shift);
    let tol = 2.0 * f.certificate();
    reset_refinements();
    let propagator = duct.propagator(dt, tol).expect("the f64 certificate passes");
    assert_eq!(refinements(), 0, "no refinement on a duct that passes today");
    let (phi0, dt_phi1) = dense(&propagator);
    let (expected_phi0, expected_phi1, _, _) = series_phis(&duct.l_tilde, n, dt, shift).unwrap();
    assert!(
        phi0.iter().zip(&expected_phi0).all(|(a, b)| a.to_bits() == b.to_bits()),
        "phi0 is the f64 series' bits"
    );
    let expected_dt_phi1 = scale_matrix(&expected_phi1, dt);
    assert!(
        dt_phi1.iter().zip(&expected_dt_phi1).all(|(a, b)| a.to_bits() == b.to_bits()),
        "dt_phi1 is the f64 series' bits"
    );
    assert_eq!(
        propagator.certificate_residual().to_bits(),
        f.certificate().to_bits(),
        "the certificate is the f64 one"
    );
}

#[test]
fn a_tolerance_the_double_double_series_cannot_meet_returns_the_f64_refusal() {
    let (duct, dt) = small_case();
    let n = duct.total_dofs();
    let shift = series_shift(&duct.l_tilde, n, dt);
    let f64_certificate = f64_series(&duct, dt, shift).certificate();
    // Below one ulp per entry: the rounding term alone exceeds it.
    let tol = 1e-30;
    reset_refinements();
    let error = duct.propagator(dt, tol).expect_err("nothing certifies 1e-30");
    assert_eq!(refinements(), 1, "the fallback ran and could not certify");
    match error {
        DuctError::PropagatorCertificate { residual, tol: refused_tol } => {
            assert_eq!(
                residual.to_bits(),
                f64_certificate.to_bits(),
                "the refusal carries the f64 certificate, not the double-double one"
            );
            assert_eq!(refused_tol.to_bits(), tol.to_bits());
        }
        other => panic!("expected the certificate refusal, got {other:?}"),
    }
}

/// 244 DOFs in double-double is minutes unoptimized: release builds only.
#[cfg(not(debug_assertions))]
#[test]
fn the_244_dof_oboe_chimney_certifies_in_double_double_where_f64_refuses() {
    use std::time::Instant;
    for radius in [1.45e-3, 1.5e-3, 1.55e-3, 1.6e-3] {
        let duct = card_chimney(radius);
        let n = duct.total_dofs();
        assert_eq!(n, 244, "the card's chimney has 244 DOFs");
        let shift = series_shift(&duct.l_tilde, n, DT);
        let f = f64_series(&duct, DT, shift);
        let f64_certificate = f.certificate();
        reset_refinements();
        let started = Instant::now();
        let result = duct.propagator(DT, 1e-11);
        let seconds = started.elapsed().as_secs_f64();
        let ran = refinements();
        let refused_in_f64 = f64_certificate > 1e-11;
        match result {
            Ok(propagator) => {
                assert_eq!(ran, usize::from(refused_in_f64), "radius {radius}");
                let certificate = propagator.certificate_residual();
                assert!(certificate <= 1e-11, "radius {radius}: {certificate:e}");
                let (phi0, _) = dense(&propagator);
                let f64_error = frobenius_difference(&f.phi0, phi0);
                if refused_in_f64 {
                    assert!(certificate < 1e-13, "radius {radius}: dd certificate {certificate:e}");
                } else {
                    assert!(
                        phi0.iter().zip(&f.phi0).all(|(a, b)| a.to_bits() == b.to_bits()),
                        "radius {radius}: a passing chimney keeps the f64 bits"
                    );
                }
                println!(
                    "{{\"suite\":\"fs-duct\",\"case\":\"series-dd-oboe-sweep\",\"radius_m\":{radius:e},\"n\":{n},\"shift\":{shift},\"carried\":{:e},\"halving\":{:e},\"f64_refused\":{refused_in_f64},\"refinements\":{ran},\"certificate\":{certificate:e},\"f64_phi0_vs_final\":{f64_error:e},\"seconds\":{seconds:.1}}}",
                    f.carried, f.halving
                );
            }
            Err(error) => panic!(
                "radius {radius}: f64 certificate {f64_certificate:e}, refinements {ran}: {error:?}"
            ),
        }
    }
    // Whatever the sweep did, one real 244-DOF fallback: a tolerance between the double-double
    // and the f64 certificate at 1.5 mm.
    let duct = card_chimney(1.5e-3);
    let n = duct.total_dofs();
    let shift = series_shift(&duct.l_tilde, n, DT);
    let f64_certificate = f64_series(&duct, DT, shift).certificate();
    assert!(f64_certificate > 1e-14, "f64 certificate {f64_certificate:e}");
    reset_refinements();
    let started = Instant::now();
    let tol = 0.5 * f64_certificate;
    let propagator = duct.propagator(DT, tol).expect("the fallback certifies");
    let seconds = started.elapsed().as_secs_f64();
    assert_eq!(refinements(), 1);
    assert!(propagator.certificate_residual() <= tol);
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"series-dd-oboe-forced\",\"radius_m\":1.5e-3,\"n\":{n},\"shift\":{shift},\"f64_certificate\":{f64_certificate:e},\"tol\":{tol:e},\"certificate\":{:e},\"seconds\":{seconds:.1}}}",
        propagator.certificate_residual()
    );
}

/// F2, measured beside F1: the f64 path with a larger scaled-norm threshold `theta` (fewer
/// squarings, more Taylor terms). One certificate per threshold on the oboe chimney; no
/// assertion beyond finiteness. Release only (244 DOFs).
#[cfg(not(debug_assertions))]
#[test]
fn f2_the_f64_certificate_against_the_scaled_norm_threshold() {
    for radius in [1.45e-3, 1.5e-3, 1.55e-3, 1.6e-3] {
        let duct = card_chimney(radius);
        let n = duct.total_dofs();
        let mut norm1 = 0.0f64;
        for j in 0..n {
            let col: f64 = (0..n).map(|i| duct.l_tilde[i * n + j].abs() * DT).sum();
            norm1 = norm1.max(col);
        }
        let mut line = format!(
            "{{\"suite\":\"fs-duct\",\"case\":\"series-f2\",\"radius_m\":{radius:e},\"n\":{n},\"norm1\":{norm1:e},\"theta\":["
        );
        for (index, theta) in [0.5_f64, 1.0, 2.0, 4.0].into_iter().enumerate() {
            let mut shift = 0u32;
            let mut scaled = norm1;
            while scaled > theta && shift < 29 {
                scaled *= 0.5;
                shift += 1;
            }
            if index == 0 {
                assert_eq!(shift, series_shift(&duct.l_tilde, n, DT), "theta 0.5 is today's shift");
            }
            let f = f64_series(&duct, DT, shift);
            assert!(f.certificate().is_finite());
            if index > 0 {
                line.push(',');
            }
            line.push_str(&format!(
                "{{\"theta\":{theta},\"shift\":{},\"carried\":{:e},\"halving\":{:e},\"certificate\":{:e}}}",
                f.shift, f.carried, f.halving, f.certificate()
            ));
        }
        line.push_str("]}");
        println!("{line}");
    }
}

/// The 9.6 kHz Db4 page's chimneys (48 poles, 50 Hz to 9.6 kHz, error 0.01; 196 DOFs) at the 96 and
/// 48 kHz clocks, over the card's radius range: every one the f64 certificate passes keeps the
/// f64 series' bits and never reaches the fallback. Release only.
#[cfg(not(debug_assertions))]
#[test]
fn the_9p6_khz_band_chimneys_that_pass_today_keep_their_bits() {
    let mut passed = 0;
    for (radius, length) in [(1.5e-3, 6.0e-3), (2.1e-3, 6.0e-3), (3.0e-3, 6.0e-3), (5.4e-3, 4.5e-3)] {
        let duct = chimney(radius, length, 48, [50.0, 9_600.0], 0.01);
        let n = duct.total_dofs();
        assert_eq!(n, 196);
        for dt in [DT, 2.0 * DT] {
            let shift = series_shift(&duct.l_tilde, n, dt);
            let f = f64_series(&duct, dt, shift);
            let certificate = f.certificate();
            reset_refinements();
            let result = duct.propagator(dt, 1e-11);
            let ran = refinements();
            println!(
                "{{\"suite\":\"fs-duct\",\"case\":\"series-dd-9p6-identity\",\"radius_m\":{radius:e},\"dt\":{dt:e},\"shift\":{shift},\"f64_certificate\":{certificate:e},\"f64_passes\":{},\"refinements\":{ran}}}",
                certificate <= 1e-11
            );
            if certificate <= 1e-11 {
                let propagator = result.expect("the f64 certificate passes");
                assert_eq!(ran, 0, "radius {radius}, dt {dt:e}: no refinement on a passing chimney");
                let (phi0, dt_phi1) = dense(&propagator);
                let (expected_phi0, expected_phi1, _, _) =
                    series_phis(&duct.l_tilde, n, dt, shift).unwrap();
                assert!(phi0.iter().zip(&expected_phi0).all(|(a, b)| a.to_bits() == b.to_bits()));
                let expected_dt_phi1 = scale_matrix(&expected_phi1, dt);
                assert!(dt_phi1.iter().zip(&expected_dt_phi1).all(|(a, b)| a.to_bits() == b.to_bits()));
                assert_eq!(propagator.certificate_residual().to_bits(), certificate.to_bits());
                passed += 1;
            }
        }
    }
    assert!(passed > 0, "no 9.6 kHz chimney passed in f64; the identity was not exercised");
}
