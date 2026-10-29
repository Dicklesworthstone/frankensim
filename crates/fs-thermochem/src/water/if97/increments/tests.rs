use super::*;

fn difference(region: If97Region, p: f64, from: f64, to: f64) -> If97IsobaricDifference {
    isobaric_difference(
        Pressure::new(p),
        region,
        Temperature::new(from),
        Temperature::new(to),
    )
    .unwrap()
}
fn state(region: If97Region, p: f64, t: f64) -> If97State {
    If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region).unwrap()
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:.17e} versus {expected:.17e}"
    );
}

#[test]
fn regional_increments_match_independent_cp_and_cp_over_t_integrals() {
    for (region, p, lo, hi) in [
        (If97Region::Region1, 3e6, 300.0, 500.0),
        (If97Region::Region1, 80e6, 290.0, 610.0),
        (If97Region::Region2, 1e6, 500.0, 900.0),
        (If97Region::Region2, 30e6, 700.0, 1000.0),
    ] {
        let result = difference(region, p, lo, hi);
        // Composite Simpson integration uses independently evaluated original
        // source Cp. It does not use the increment's factored Gibbs algebra.
        const N: usize = 2048;
        let dt = (hi - lo) / N as f64;
        let (mut h_integral, mut s_integral) = (0.0, 0.0);
        for i in 0..=N {
            let t = if i == N { hi } else { lo + i as f64 * dt };
            let cp = state(region, p, t).cp().value();
            let weight = if i == 0 || i == N {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            h_integral += weight * cp;
            s_integral += weight * cp / t;
        }
        close(
            result.specific_enthalpy_change().value(),
            h_integral * dt / 3.0,
            1e-6,
            2e-10,
        );
        close(
            result.specific_entropy_change().value(),
            s_integral * dt / 3.0,
            1e-9,
            2e-10,
        );
        let (a, b) = (state(region, p, lo), state(region, p, hi));
        close(
            result.specific_enthalpy_change().value(),
            b.specific_enthalpy().value() - a.specific_enthalpy().value(),
            2e-7,
            1e-11,
        );
        close(
            result.specific_entropy_change().value(),
            b.specific_entropy().value() - a.specific_entropy().value(),
            1e-9,
            1e-11,
        );
        close(
            result.specific_volume_change().value(),
            b.specific_volume().value() - a.specific_volume().value(),
            1e-15,
            1e-10,
        );
        close(
            result.specific_internal_energy_change().value(),
            b.specific_internal_energy().value() - a.specific_internal_energy().value(),
            2e-7,
            1e-11,
        );
    }
}

#[test]
fn adjacent_temperatures_keep_sensible_heat_and_entropy_beside_absolute_references() {
    let mut absolute_subtraction_is_inaccurate = false;
    for (region, p, t) in [
        (If97Region::Region1, 3e6, 300.0_f64),
        (If97Region::Region1, 80e6, 600.0),
        (If97Region::Region2, 100_000.0, 600.0),
        (If97Region::Region2, 30e6, 700.0),
    ] {
        let next = t.next_up();
        let initial = state(region, p, t);
        let expected_h = initial.cp().value() * (next - t);
        let result = difference(region, p, t, next);
        assert!(result.specific_enthalpy_change().value() > 0.0);
        close(
            result.specific_enthalpy_change().value(),
            expected_h,
            0.0,
            1e-10,
        );
        close(
            result.specific_entropy_change().value(),
            expected_h / t,
            0.0,
            1e-10,
        );
        close(
            result.specific_volume_change().value(),
            initial.dv_dt_at_pressure().value() * (next - t),
            0.0,
            1e-9,
        );
        close(
            result.specific_internal_energy_change().value()
                + p * result.specific_volume_change().value(),
            result.specific_enthalpy_change().value(),
            0.0,
            2e-15,
        );
        let raw = state(region, p, next).specific_enthalpy().value()
            - initial.specific_enthalpy().value();
        absolute_subtraction_is_inaccurate |= (raw - expected_h).abs() > 0.1 * expected_h;
    }
    assert!(
        absolute_subtraction_is_inaccurate,
        "fixture must expose the absolute-subtraction failure"
    );
}

#[test]
fn isobaric_increments_are_additive_reversible_and_zero_at_equal_endpoints() {
    for (region, p, a, b, c) in [
        (If97Region::Region1, 80e6, 300.0, 450.0, 600.0),
        (If97Region::Region2, 100_000.0, 400.0, 600.0, 1000.0),
    ] {
        let (ab, bc, ac, ca) = (
            difference(region, p, a, b),
            difference(region, p, b, c),
            difference(region, p, a, c),
            difference(region, p, c, a),
        );
        for (x, y, whole, reverse) in [
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
                ab.volume_m3_kg,
                bc.volume_m3_kg,
                ac.volume_m3_kg,
                ca.volume_m3_kg,
            ),
            (
                ab.internal_energy_j_kg,
                bc.internal_energy_j_kg,
                ac.internal_energy_j_kg,
                ca.internal_energy_j_kg,
            ),
        ] {
            close(x + y, whole, 1e-15, 1e-11);
            assert_eq!((-whole).to_bits(), reverse.to_bits());
        }
        let zero = difference(region, p, b, b);
        assert_eq!(
            [
                zero.enthalpy_j_kg,
                zero.entropy_j_kg_k,
                zero.volume_m3_kg,
                zero.internal_energy_j_kg
            ],
            [0.0; 4]
        );
    }
}

#[test]
fn regional_changes_refuse_unimplemented_dense_or_metastable_paths() {
    for (region, p, from, to) in [
        (If97Region::Region1, 80e6, 300.0, 650.0),
        (If97Region::Region2, 3e6, 300.0, 600.0),
        (If97Region::Region2, 100_000.0, 400.0, 1500.0),
        (If97Region::Region2, 0.0, 400.0, 500.0),
    ] {
        assert!(
            isobaric_difference(
                Pressure::new(p),
                region,
                Temperature::new(from),
                Temperature::new(to)
            )
            .is_err()
        );
    }
}
