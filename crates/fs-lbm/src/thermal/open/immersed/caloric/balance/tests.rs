use super::*;

type Error = CaloricBalanceError<&'static str>;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn sensible_constant_capacity_matches_the_independent_backward_euler_balance() {
    for (old_t, trace, source) in [(0.8, 0.2, 0.03), (0.2, 1.4, -0.06), (1.0, 1.0, 0.0)] {
        let cp = 4.0;
        let mass = 2.0;
        let conductance = 1.3;
        let expected =
            (mass * cp * old_t + source + conductance * trace) / (mass * cp + conductance);
        let result = solve_balance(
            cp * old_t,
            [-8.0, 16.0],
            mass,
            source,
            conductance,
            trace,
            2e-14,
            32,
            &mut |h| Ok::<_, Error>((h / cp, cp.recip())),
            &mut || Ok(()),
            &|e| e,
        )
        .unwrap();
        near(result.1, expected, 5e-15);
        near(
            mass * (result.0 - cp * old_t) + conductance * (result.1 - trace),
            source,
            2e-14,
        );
    }
}

#[test]
fn nonlinear_capacity_uses_the_enthalpy_primitive_and_preserves_datum_shifts() {
    // Cp(T)=1+T, hence h=T+T^2/2. The independent endpoint equation is
    // T^2+6T-18=0, not a fixed-C update evaluated at either endpoint.
    let expected = (108.0_f64.sqrt() - 6.0) / 2.0;
    for datum in [0.0, 1000.0] {
        let (h, t) = solve_balance(
            datum + 1.5,
            [datum, datum + 17.5],
            2.0,
            3.0,
            4.0,
            3.0,
            1e-12,
            64,
            &mut |h| {
                let t = (1.0 + 2.0 * (h - datum)).sqrt() - 1.0;
                Ok::<_, Error>((t, 1.0 / (1.0 + t)))
            },
            &mut || Ok(()),
            &|e| e,
        )
        .unwrap();
        near(t, expected, 2e-13);
        near(h - datum, expected + 0.5 * expected * expected, 2e-13);
    }
}

fn latent(h: f64) -> Result<(f64, f64), Error> {
    Ok(if h < 2.0 {
        (h / 2.0, 0.5)
    } else if h <= 10.0 {
        (1.0, 0.0)
    } else {
        (1.0 + (h - 10.0) / 3.0, 1.0 / 3.0)
    })
}

#[test]
fn latent_storage_holds_exact_temperature_and_crosses_original_phase_knots() {
    for (source, trace) in [(3.0, 0.2), (-3.0, 0.2), (0.0, 1.0)] {
        let (h, t) = solve_balance(
            6.0,
            [0.0, 16.0],
            4.0,
            source,
            2.0,
            trace,
            1e-13,
            32,
            &mut latent,
            &mut || Ok(()),
            &|e| e,
        )
        .unwrap();
        assert_eq!(t, 1.0);
        near(h, 6.0 + (source - 2.0 * (1.0 - trace)) / 4.0, 1e-14);
    }
    for source in [-24.0, 32.0] {
        let (h, t) = solve_balance(
            6.0,
            [0.0, 16.0],
            4.0,
            source,
            2.0,
            0.2,
            1e-12,
            32,
            &mut latent,
            &mut || Ok(()),
            &|e| e,
        )
        .unwrap();
        assert_eq!(h < 2.0, source < 0.0);
        assert_eq!(h > 10.0, source > 0.0);
        near(4.0 * (h - 6.0) + 2.0 * (t - 0.2), source, 1e-12);
    }
}

#[test]
fn source_refusal_cancellation_work_and_representability_remain_terminal() {
    let mut calls = 0;
    let failed = solve_balance(
        6.0,
        [0.0, 16.0],
        4.0,
        100.0,
        2.0,
        0.2,
        1e-12,
        32,
        &mut |h| {
            calls += 1;
            latent(h)
        },
        &mut || Ok(()),
        &|e| e,
    );
    assert!(matches!(
        failed,
        Err(Error::OutsideEnthalpyRange {
            bound_j_kg: 16.0,
            ..
        })
    ));
    assert_eq!(calls, 2);

    let source = solve_balance(
        6.0,
        [0.0, 16.0],
        4.0,
        1.0,
        2.0,
        0.2,
        1e-12,
        32,
        &mut |_| Err(Error::Material("source refused")),
        &mut || Ok(()),
        &|e| e,
    );
    assert_eq!(source, Err(Error::Material("source refused")));
    calls = 0;
    let cancelled = solve_balance(
        6.0,
        [0.0, 16.0],
        4.0,
        1.0,
        2.0,
        0.2,
        1e-12,
        32,
        &mut |h| {
            calls += 1;
            latent(h)
        },
        &mut || Err(Error::Material("cancelled")),
        &|e| e,
    );
    assert_eq!(cancelled, Err(Error::Material("cancelled")));
    assert_eq!(calls, 0);

    let exhausted = solve_balance(
        6.0,
        [0.0, 16.0],
        4.0,
        1.0,
        2.0,
        0.2,
        1e-12,
        1,
        &mut latent,
        &mut || Ok(()),
        &|e| e,
    );
    assert!(matches!(
        exhausted,
        Err(Error::NonConvergence { iterations: 1, .. })
    ));

    // Required energy is smaller than one representable enthalpy increment.
    // An unrounded hidden increment must not masquerade as accepted storage.
    let unrepresented = solve_balance(
        1e30,
        [0.0, 2e30],
        1.0,
        1.0,
        1.0,
        0.0,
        1e-12,
        256,
        &mut |_| Ok::<_, Error>((0.0, 0.0)),
        &mut || Ok(()),
        &|e| e,
    );
    assert!(matches!(unrepresented, Err(Error::NonConvergence { .. })));
}
