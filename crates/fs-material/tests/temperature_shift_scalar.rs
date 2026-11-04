//! G0/G1/G3: one calibrated shift reaches the original scalar physics owners.
//! Numerical fixture coefficients make no claim of measured material accuracy.

use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco::{
    FractionalZener, GeneralizedMaxwell, ModifiedFractionalZener, ViscoError, lower_to_prony,
};

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 4.0,
            c2_k: 80.0,
        },
        300.0,
        [270.0, 380.0],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_scalar_spectra_use_the_same_temperature_frequency_map() {
    let maxwell = GeneralizedMaxwell::new(2e8, vec![(3e7, 0.01), (7e7, 0.2)]).unwrap();
    let fractional = FractionalZener::new(2e8, 3e8, 0.4, 0.02).unwrap();
    for temperature in [280.0, 300.0, 330.0] {
        let a = shift().evaluate(temperature).unwrap().factor();
        let gm = shift().shift_maxwell(&maxwell, temperature).unwrap();
        let fz = shift()
            .shift_fractional_zener(&fractional, temperature)
            .unwrap();
        assert_eq!(gm.law().e_inf.to_bits(), maxwell.e_inf.to_bits());
        assert_eq!(fz.law().e0.to_bits(), fractional.e0.to_bits());
        assert_eq!(fz.law().e_inf.to_bits(), fractional.e_inf.to_bits());
        assert_eq!(fz.law().alpha.to_bits(), fractional.alpha.to_bits());
        for (original, shifted) in maxwell.terms.iter().zip(&gm.law().terms) {
            assert_eq!(original.0.to_bits(), shifted.0.to_bits());
            assert_eq!(shifted.1, a * original.1);
        }
        for omega in [-10000.0, -10.0, 0.0, 0.1, 10.0, 1000.0, 10000.0] {
            for (observed, expected) in [
                (gm.law().modulus(omega), maxwell.modulus(a * omega)),
                (fz.law().modulus(omega), fractional.modulus(a * omega)),
            ] {
                close(observed.0, expected.0);
                close(observed.1, expected.1);
            }
        }
    }
    assert_eq!(
        shift().shift_maxwell(&maxwell, 300.0).unwrap().law(),
        &maxwell
    );
}

#[test]
fn g1_original_scalar_history_has_reduced_time_and_nonnegative_heat() {
    let reference = GeneralizedMaxwell::new(2e8, vec![(3e7, 0.01), (7e7, 0.2)]).unwrap();
    for temperature in [280.0, 330.0] {
        let shifted = shift().shift_maxwell(&reference, temperature).unwrap();
        let a = shifted.state().factor();
        let mut actual = shifted.law().state();
        let mut expected = reference.state();
        for (strain, dt) in [(0.002, 0.0), (0.002, 0.03), (-0.001, 0.02), (0.0, 0.1)] {
            shifted.law().step(&mut actual, strain, dt);
            reference.step(&mut expected, strain, dt / a);
            close(actual.stress, expected.stress);
            close(actual.work, expected.work);
            close(actual.dissipated, expected.dissipated);
            for (&observed, &control) in actual.q.iter().zip(&expected.q) {
                close(observed, control);
            }
            let stored = 0.5 * reference.e_inf * strain * strain
                + actual
                    .q
                    .iter()
                    .zip(&reference.terms)
                    .map(|(q, (modulus, _))| 0.5 * modulus * q * q)
                    .sum::<f64>();
            close(actual.work, stored + actual.dissipated);
            assert!(actual.dissipated >= 0.0);
        }
    }
    let cold = shift().shift_maxwell(&reference, 280.0).unwrap();
    let hot = shift().shift_maxwell(&reference, 330.0).unwrap();
    let mut cold_state = cold.law().state();
    let mut hot_state = hot.law().state();
    for dt in [0.0, 0.03] {
        cold.law().step(&mut cold_state, 0.002, dt);
        hot.law().step(&mut hot_state, 0.002, dt);
    }
    assert!(hot_state.stress < cold_state.stress);
    assert!(hot_state.dissipated > cold_state.dissipated);
}

#[test]
fn g3_lowered_and_modified_laws_keep_their_source_bands_and_fit_measurements() {
    // Same admitted lowering fixture used by the existing owner's own battery.
    let fractional = FractionalZener::new(8e9, 16e9, 0.35, 1e-3).unwrap();
    let reference = lower_to_prony(&fractional, 20.0, 2e5, 8, 0.02).unwrap();
    let modified = ModifiedFractionalZener::new(2e8, 3e8, 0.7, 0.4, 0.01, (20.0, 2e4)).unwrap();
    for temperature in [280.0, 330.0] {
        let a = shift().evaluate(temperature).unwrap().factor();
        let lowered = shift().shift_lowered(&reference, temperature).unwrap();
        let original_band = reference.band;
        assert_eq!(
            lowered.law().band,
            (original_band.0 / a, original_band.1 / a)
        );
        assert_eq!(
            lowered.law().sup_rel_err.to_bits(),
            reference.sup_rel_err.to_bits()
        );
        assert_eq!(
            lowered.law().sup_rel_storage_err.to_bits(),
            reference.sup_rel_storage_err.to_bits()
        );
        assert_eq!(
            lowered.law().sup_rel_loss_err.to_bits(),
            reference.sup_rel_loss_err.to_bits()
        );
        assert_eq!(
            lowered.law().verification_points,
            reference.verification_points
        );
        for fraction in [0.01, 0.2, 0.5, 0.9, 0.99] {
            let omega_ref = original_band.0 * (original_band.1 / original_band.0).powf(fraction);
            let actual = lowered.law().modulus_checked(omega_ref / a).unwrap();
            let expected = reference.modulus_checked(omega_ref).unwrap();
            close(actual.0, expected.0);
            close(actual.1, expected.1);
        }
        assert!(matches!(
            lowered
                .law()
                .modulus_checked(lowered.law().band.0.next_down()),
            Err(ViscoError::OutOfBand { .. })
        ));
        assert!(matches!(
            lowered
                .law()
                .modulus_checked(lowered.law().band.1.next_up()),
            Err(ViscoError::OutOfBand { .. })
        ));

        let shifted = shift()
            .shift_modified_zener(&modified, temperature)
            .unwrap();
        let band = modified.band_hz();
        assert_eq!(shifted.law().band_hz(), (band.0 / a, band.1 / a));
        for omega_ref in [200.0, 1000.0, 10000.0, 100000.0] {
            let actual = shifted.law().modulus(omega_ref / a).unwrap();
            let expected = modified.modulus(omega_ref).unwrap();
            close(actual.0, expected.0);
            close(actual.1, expected.1);
        }
        assert!(matches!(
            shifted
                .law()
                .modulus(shifted.law().band_rad_s().0.next_down()),
            Err(ViscoError::OutOfBand { .. })
        ));
    }
}

#[test]
fn g0_forged_public_coefficients_and_unrepresentable_shifts_refuse() {
    let valid = GeneralizedMaxwell::new(2e8, vec![(3e7, 0.01)]).unwrap();
    let saved = valid.clone();
    assert!(shift().shift_maxwell(&valid, 269.0).is_err());
    let invalid = GeneralizedMaxwell {
        e_inf: -1.0,
        terms: vec![(1.0, 0.1)],
    };
    assert!(shift().shift_maxwell(&invalid, 300.0).is_err());
    for (tau, temperature) in [(f64::MAX, 280.0), (f64::from_bits(1), 330.0)] {
        let original = GeneralizedMaxwell::new(1.0, vec![(1.0, tau)]).unwrap();
        assert!(shift().shift_maxwell(&original, temperature).is_err());
    }
    let fractional = FractionalZener {
        e0: 1.0,
        e_inf: 2.0,
        alpha: 1.5,
        tau: 0.1,
    };
    assert!(shift().shift_fractional_zener(&fractional, 300.0).is_err());
    assert_eq!(valid, saved);
}
