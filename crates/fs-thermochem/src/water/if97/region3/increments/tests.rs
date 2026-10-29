use super::*;

fn state(rho: f64, t: f64) -> If97DenseState {
    If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)).unwrap()
}
fn prepared(rho: f64, lo: f64, hi: f64) -> DenseIsochoricInterval {
    DenseIsochoricInterval::new(
        Density::new(rho),
        Temperature::new(lo),
        Temperature::new(hi),
        DenseIsochoricIntervalSolve {
            maximum_intervals: MAX_INTERVAL_EVALUATIONS,
        },
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, abs: f64, rel: f64) {
    assert!(
        (a - b).abs() <= abs + rel * b.abs(),
        "{a:.17e} versus {b:.17e}"
    );
}

#[test]
fn whole_isochore_bounds_contain_independent_source_queries() {
    for (rho, lo, hi) in [
        (500.0, 650.0, 750.0),
        (200.0, 650.0, 660.0),
        (700.0, 625.0, 635.0),
    ] {
        let path = prepared(rho, lo, hi);
        assert!(path.minimum_heat_capacity_j_kg_k() > 0.0);
        assert!(path.interval_evaluations() > 1);
        for i in 0..=128 {
            let t = lo + (hi - lo) * i as f64 / 128.0;
            let cv = state(rho, t).cv().value();
            assert!(cv >= path.minimum_heat_capacity_j_kg_k());
            assert!(cv <= path.maximum_heat_capacity_j_kg_k());
            assert!(1.0 / cv <= path.maximum_temperature_slope());
        }
    }
}

#[test]
fn stable_increments_match_independent_cv_and_cv_over_t_integrals() {
    for (rho, lo, hi) in [
        (500.0, 650.0, 750.0),
        (200.0, 650.0, 660.0),
        (700.0, 625.0, 635.0),
    ] {
        let path = prepared(rho, lo, hi);
        let difference = path
            .difference(Temperature::new(lo), Temperature::new(hi))
            .unwrap();
        const N: usize = 1024;
        let dt = (hi - lo) / N as f64;
        let (mut u, mut s, mut p) = (0.0, 0.0, 0.0);
        for i in 0..=N {
            let t = if i == N { hi } else { lo + dt * i as f64 };
            let point = state(rho, t);
            let weight = if i == 0 || i == N {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            u += weight * point.cv().value();
            s += weight * point.cv().value() / t;
            p += weight * point.dp_dt_at_density().value();
        }
        close(difference.internal_energy_j_kg, u * dt / 3.0, 1e-6, 2e-10);
        close(difference.entropy_j_kg_k, s * dt / 3.0, 1e-9, 2e-10);
        close(difference.pressure_pa, p * dt / 3.0, 1e-3, 2e-10);
        let (a, b) = (state(rho, lo), state(rho, hi));
        close(
            difference.internal_energy_j_kg,
            b.specific_internal_energy().value() - a.specific_internal_energy().value(),
            2e-6,
            1e-10,
        );
        close(
            difference.enthalpy_j_kg,
            b.specific_enthalpy().value() - a.specific_enthalpy().value(),
            2e-6,
            1e-10,
        );
        close(
            difference.pressure_pa,
            b.pressure().value() - a.pressure().value(),
            1e-3,
            1e-10,
        );
        close(
            difference.enthalpy_j_kg,
            difference.internal_energy_j_kg + difference.pressure_pa / rho,
            0.0,
            2e-15,
        );
    }
}

#[test]
fn adjacent_temperatures_retain_internal_energy_pressure_and_entropy_changes() {
    let mut naive_subtraction_fails = false;
    for (rho, t) in [(500.0, 650.0_f64), (200.0, 650.0), (700.0, 630.0)] {
        let b = t.next_up();
        let path = prepared(rho, t, b);
        let d = path
            .difference(Temperature::new(t), Temperature::new(b))
            .unwrap();
        let a = state(rho, t);
        let expected = a.cv().value() * (b - t);
        assert!(d.internal_energy_j_kg > 0.0);
        close(d.internal_energy_j_kg, expected, 0.0, 1e-10);
        close(d.entropy_j_kg_k, expected / t, 0.0, 1e-10);
        close(
            d.pressure_pa,
            a.dp_dt_at_density().value() * (b - t),
            0.0,
            1e-10,
        );
        let raw =
            state(rho, b).specific_internal_energy().value() - a.specific_internal_energy().value();
        naive_subtraction_fails |= (raw - expected).abs() > 0.1 * expected;
    }
    assert!(
        naive_subtraction_fails,
        "fixture must expose rounded endpoint subtraction loss"
    );
}

#[test]
fn increments_are_additive_reversible_and_exact_zero() {
    let path = prepared(500.0, 650.0, 750.0);
    let diff = |a, b| {
        path.difference(Temperature::new(a), Temperature::new(b))
            .unwrap()
    };
    let (ab, bc, ac, ca) = (
        diff(650.0, 700.0),
        diff(700.0, 750.0),
        diff(650.0, 750.0),
        diff(750.0, 650.0),
    );
    for (x, y, sum, reverse) in [
        (
            ab.internal_energy_j_kg,
            bc.internal_energy_j_kg,
            ac.internal_energy_j_kg,
            ca.internal_energy_j_kg,
        ),
        (
            ab.enthalpy_j_kg,
            bc.enthalpy_j_kg,
            ac.enthalpy_j_kg,
            ca.enthalpy_j_kg,
        ),
        (
            ab.entropy_j_kg_k,
            bc.entropy_j_kg_k,
            ac.entropy_j_kg_k,
            ca.entropy_j_kg_k,
        ),
        (
            ab.pressure_pa,
            bc.pressure_pa,
            ac.pressure_pa,
            ca.pressure_pa,
        ),
    ] {
        close(x + y, sum, 1e-8, 2e-11);
        assert_eq!((-sum).to_bits(), reverse.to_bits());
    }
    let zero = diff(700.0, 700.0);
    assert_eq!(
        [
            zero.internal_energy_j_kg,
            zero.enthalpy_j_kg,
            zero.entropy_j_kg_k,
            zero.pressure_pa
        ],
        [0.0; 4]
    );
    assert_eq!(
        isochoric_difference(
            Density::new(500.0),
            Temperature::new(650.0),
            Temperature::new(750.0)
        )
        .unwrap(),
        ac
    );
}

#[test]
fn source_domain_work_and_cancellation_refuse_without_partial_path() {
    let controls = DenseIsochoricIntervalSolve {
        maximum_intervals: MAX_INTERVAL_EVALUATIONS,
    };
    for (rho, lo, hi) in [
        (500.0, 600.0, 750.0),
        (500.0, 650.0, 900.0),
        (1000.0, 650.0, 750.0),
        (322.0, 647.096, 650.0),
    ] {
        assert!(
            DenseIsochoricInterval::new(
                Density::new(rho),
                Temperature::new(lo),
                Temperature::new(hi),
                controls,
                || false
            )
            .is_err()
        );
    }
    let low_work = DenseIsochoricInterval::new(
        Density::new(500.0),
        Temperature::new(650.0),
        Temperature::new(750.0),
        DenseIsochoricIntervalSolve {
            maximum_intervals: 1,
        },
        || false,
    );
    assert!(matches!(
        low_work,
        Err(DenseIsochoricIntervalError::NotProven { .. })
    ));
    let mut polls = 0;
    let cancelled = DenseIsochoricInterval::new(
        Density::new(500.0),
        Temperature::new(650.0),
        Temperature::new(750.0),
        controls,
        || {
            polls += 1;
            polls == 4
        },
    );
    assert_eq!(cancelled, Err(DenseIsochoricIntervalError::Cancelled));
    assert_eq!(polls, 4);
    let path = prepared(500.0, 650.0, 750.0);
    assert!(
        path.difference(Temperature::new(649.0), Temperature::new(700.0))
            .is_err()
    );
}
