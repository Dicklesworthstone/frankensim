use super::*;
use fs_qty::{Density, Temperature};

fn state(rho: f64, t: f64) -> If97DenseState {
    If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)).unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let limit = absolute + relative * expected.abs();
    assert!(
        actual.is_finite() && (actual - expected).abs() <= limit,
        "actual={actual:.16e}, expected={expected:.16e}, limit={limit:.5e}"
    );
}
fn values(d: DenseStateDifference) -> [f64; 6] {
    [
        d.specific_internal_energy_change().value(),
        d.specific_enthalpy_change().value(),
        d.specific_entropy_change().value(),
        d.specific_volume_change().value(),
        d.pressure_change().value(),
        d.flow_work_change().value(),
    ]
}

#[test]
fn adjacent_density_and_temperature_increments_retain_analytic_response() {
    // G1: adjacent binary64 inputs are much closer than the rounding error of
    // separately evaluated absolute energies. Derivatives supply an independent
    // local limit, including simultaneously changing density and temperature.
    for (rho, t) in [(500.0_f64, 650.0_f64), (200.0, 650.0), (700.0, 630.0)] {
        let a = state(rho, t);
        for (r, temp) in [
            (rho.next_up(), t),
            (rho, t.next_up()),
            (rho.next_up(), t.next_up()),
        ] {
            let b = state(r, temp);
            let (dr, dt) = (r - rho, temp - t);
            let dp = a.dp_drho_at_temperature().value() * dr + a.dp_dt_at_density().value() * dt;
            let du = a.du_drho_at_temperature().value() * dr + a.cv().value() * dt;
            let dv = -dr / (rho * rho);
            let dpv = dp / rho + a.pressure().value() * dv;
            let ds = a.cv().value() / t * dt - a.dp_dt_at_density().value() / (rho * rho) * dr;
            let expected = [du, du + dpv, ds, dv, dp, dpv];
            for (got, want) in values(difference(a, b)).into_iter().zip(expected) {
                near(got, want, 1e-29, 5e-10);
            }
            for (forward, reverse) in values(difference(a, b))
                .into_iter()
                .zip(values(difference(b, a)))
            {
                assert_eq!(forward, -reverse);
            }
        }
        assert!(values(difference(a, a)).into_iter().all(|x| x == 0.0));
    }
}

#[test]
fn original_endpoint_values_and_prepared_isochore_are_preserved() {
    // G0/G3: changes keep the existing released caloric references, even for
    // endpoints whose connecting straight line is not an admitted source path.
    let a = state(500.0, 650.0);
    for b in [
        state(500.0, 750.0),
        state(200.0, 650.0),
        state(700.0, 630.0),
    ] {
        let raw = [
            b.specific_internal_energy().value() - a.specific_internal_energy().value(),
            b.specific_enthalpy().value() - a.specific_enthalpy().value(),
            b.specific_entropy().value() - a.specific_entropy().value(),
            b.specific_volume().value() - a.specific_volume().value(),
            b.pressure().value() - a.pressure().value(),
            b.pressure().value() / b.density().value() - a.pressure().value() / a.density().value(),
        ];
        for ((got, want), abs) in values(difference(a, b))
            .into_iter()
            .zip(raw)
            .zip([2e-6, 2e-6, 5e-9, 1e-18, 0.001, 1e-6])
        {
            near(got, want, abs, 3e-12);
        }
    }
    let existing = super::super::increments::isochoric_difference(
        Density::new(500.0),
        Temperature::new(650.0),
        Temperature::new(750.0),
    )
    .unwrap();
    let both = difference(a, state(500.0, 750.0));
    near(
        both.specific_internal_energy_change().value(),
        existing.specific_internal_energy_change().value(),
        5e-8,
        1e-13,
    );
    near(
        both.specific_entropy_change().value(),
        existing.specific_entropy_change().value(),
        1e-10,
        1e-13,
    );
}

#[test]
fn simultaneous_changes_match_first_law_and_entropy_path_integrals() {
    // G1: independently integrate the EOS derivatives along a stable path.
    // The implementation differences monomials; it never performs this integral.
    let (ra, ta, rb, tb) = (500.0, 650.0, 450.0, 700.0);
    let (dr, dt) = (rb - ra, tb - ta);
    let n = 512;
    let mut integrated = [0.0; 6];
    for i in 0..=n {
        let fraction = i as f64 / n as f64;
        let (rho, t) = (ra + fraction * dr, ta + fraction * dt);
        let s = state(rho, t);
        let dp = s.dp_drho_at_temperature().value() * dr + s.dp_dt_at_density().value() * dt;
        let du = s.du_drho_at_temperature().value() * dr + s.cv().value() * dt;
        let dv = -dr / (rho * rho);
        let dpv = dp / rho + s.pressure().value() * dv;
        let ds = s.cv().value() / t * dt - s.dp_dt_at_density().value() / (rho * rho) * dr;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in integrated.iter_mut().zip([du, du + dpv, ds, dv, dp, dpv]) {
            *sum += weight * value / (3.0 * n as f64);
        }
    }
    let a = state(ra, ta);
    let b = state(rb, tb);
    for (got, want) in values(difference(a, b)).into_iter().zip(integrated) {
        near(got, want, 1e-16, 2e-10);
    }
    let mid = state((ra + rb) * 0.5, (ta + tb) * 0.5);
    for ((whole, first), second) in values(difference(a, b))
        .into_iter()
        .zip(values(difference(a, mid)))
        .zip(values(difference(mid, b)))
    {
        near(whole, first + second, 1e-16, 3e-11);
    }
}
