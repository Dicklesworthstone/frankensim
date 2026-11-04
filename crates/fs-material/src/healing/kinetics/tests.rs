use super::*;
use crate::healing::ContactHistory;
use crate::visco::temperature::TemperatureShiftModel;

fn controls() -> MobilityIntegration {
    MobilityIntegration {
        absolute_dose_tolerance: 1e-12,
        relative_dose_tolerance: 1e-10,
        maximum_evaluations: 4097,
        maximum_depth: 24,
    }
}

fn arrhenius() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 60_000.0,
        },
        400.0,
        [330.0, 600.0],
    )
    .unwrap()
}

fn wlf() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 100.0,
        },
        400.0,
        [330.0, 600.0],
    )
    .unwrap()
}

fn law(shift: TemperatureShift, cutoff: Option<f64>) -> ShiftedMobility {
    ShiftedMobility::new(
        "synthetic-weld-kinetics-v1".into(),
        shift,
        0.25,
        cutoff,
        controls(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn original_shift_accelerates_the_welding_clock_and_preserves_exact_reference() {
    for shift in [arrhenius(), wlf()] {
        let kinetics = law(shift, None);
        assert_eq!(kinetics.shift(), &shift);
        assert_eq!(kinetics.identity(), "synthetic-weld-kinetics-v1");
        assert_eq!(kinetics.reference_rate_per_s(), 0.25);
        let reference = kinetics.integrate(400.0, 400.0, 8.0).unwrap();
        assert_eq!(
            reference,
            MobilityDose {
                dose: 2.0,
                estimated_absolute_error: 0.0,
                evaluations: 1
            }
        );
        let cold = kinetics.rate_per_s(350.0).unwrap();
        let hot = kinetics.rate_per_s(500.0).unwrap();
        assert!(cold < 0.25 && hot > 0.25);
        for temperature in [350.0, 400.0, 500.0] {
            let canonical = 0.25 / shift.evaluate(temperature).unwrap().factor();
            assert_eq!(kinetics.rate_per_s(temperature).unwrap(), canonical);
            assert_eq!(
                kinetics
                    .dose_increment(temperature, temperature, 0.125)
                    .unwrap(),
                canonical * 0.125
            );
        }
    }
}

#[test]
fn nonlinear_ramps_agree_with_independent_high_precision_integrals() {
    let steep = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 8.0,
            c2_k: 60.0,
        },
        400.0,
        [360.0, 550.0],
    )
    .unwrap();
    // Independent 65-digit Decimal integration of the defining scalar laws,
    // using 64- and 96-point Gauss-Legendre rules. Their discrepancies were
    // < 1e-38 dose. The last ramp spans more than 21 decades of molecular rate.
    // Arrhenius uses the exact binary64 gas constant 8.31446261815324 J/(mol K).
    for (shift, start, end, duration, expected) in [
        (arrhenius(), 330.0, 550.0, 2.0, 11.417_675_046_586_117),
        (wlf(), 350.0, 600.0, 0.6, 1.254_286_367_050_357_2),
        (steep, 360.0, 550.0, 0.002, 40.701_253_700_441_85),
    ] {
        let kinetics = law(shift, None);
        let result = kinetics.integrate(start, end, duration).unwrap();
        close(result.dose, expected, 2e-10 * expected);
        assert!(
            result.estimated_absolute_error
                <= controls().absolute_dose_tolerance
                    + controls().relative_dose_tolerance * result.dose
        );
        assert!(result.evaluations > 5 && result.evaluations <= controls().maximum_evaluations);
        assert_eq!(result, kinetics.integrate(end, start, duration).unwrap());
        assert_eq!(result, kinetics.integrate(start, end, duration).unwrap());
    }
}

#[test]
fn explicit_cutoff_splits_the_path_without_smearing_the_active_rate() {
    let kinetics = law(wlf(), Some(350.0));
    assert_eq!(kinetics.cutoff_temperature_k(), Some(350.0));
    assert_eq!(kinetics.rate_per_s(350.0).unwrap(), 0.0);
    assert!(kinetics.rate_per_s(350.0001).unwrap() > 0.0);
    assert_eq!(
        kinetics.integrate(330.0, 350.0, 1000.0).unwrap(),
        MobilityDose {
            dose: 0.0,
            estimated_absolute_error: 0.0,
            evaluations: 0
        }
    );
    assert_eq!(kinetics.dose_increment(350.0, 350.0, 1000.0).unwrap(), 0.0);
    let ramp = kinetics.integrate(330.0, 450.0, 0.6).unwrap();
    close(ramp.dose, 0.183_067_380_238_716_9, 2e-11);
    assert_eq!(ramp, kinetics.integrate(450.0, 330.0, 0.6).unwrap());
    close(
        ramp.dose,
        law(wlf(), None).dose_increment(350.0, 450.0, 0.5).unwrap(),
        2e-14,
    );
    // The cutoff does not license WLF extrapolation or bypass its calibration.
    assert_eq!(
        kinetics.rate_per_s(329.0),
        Err(HealingError::Temperature(329.0))
    );
    assert!(kinetics.integrate(329.0, 340.0, 0.0).is_err());
}

#[test]
fn dose_survives_thermal_path_subdivision_and_cold_restarts() {
    for shift in [arrhenius(), wlf()] {
        let kinetics = law(shift, Some(350.0));
        let whole = kinetics.dose_increment(330.0, 550.0, 0.8).unwrap();
        for count in [2, 7, 31] {
            let mut state = ContactHistory::new(&kinetics).unwrap();
            for i in 0..count {
                let start = 330.0 + 220.0 * f64::from(i) / f64::from(count);
                let end = 330.0 + 220.0 * f64::from(i + 1) / f64::from(count);
                state = state.advance(start, end, 0.8 / f64::from(count)).unwrap();
            }
            close(state.dose(), whole, 2e-9 * whole);
            let cold = state.advance(340.0, 340.0, 20.0).unwrap();
            assert_eq!(cold.dose().to_bits(), state.dose().to_bits());
            let resumed = cold.advance(400.0, 400.0, 0.5).unwrap();
            close(resumed.dose(), state.dose() + 0.125, 1e-14);
            close(resumed.contact_time_s(), 21.3, 1e-13);
        }
    }
}

#[test]
fn exhausted_work_or_invalid_paths_do_not_age_accepted_contact() {
    let limited = ShiftedMobility::new(
        "limited-weld-v1".into(),
        wlf(),
        0.25,
        None,
        MobilityIntegration {
            maximum_evaluations: 5,
            ..controls()
        },
    )
    .unwrap();
    let accepted = ContactHistory::new(&limited)
        .unwrap()
        .advance(400.0, 400.0, 0.5)
        .unwrap();
    assert_eq!(
        accepted.advance(350.0, 600.0, 1.0).unwrap_err(),
        HealingError::IntegrationBudget
    );
    assert_eq!(accepted.dose(), 0.125);
    assert_eq!(accepted.contact_time_s(), 0.5);
    let shallow = ShiftedMobility::new(
        "shallow-weld-v1".into(),
        wlf(),
        0.25,
        None,
        MobilityIntegration {
            maximum_depth: 1,
            ..controls()
        },
    )
    .unwrap();
    assert_eq!(
        shallow.integrate(350.0, 600.0, 1.0),
        Err(HealingError::IntegrationBudget)
    );
    for (start, end, duration) in [
        (329.0, 400.0, 1.0),
        (400.0, 601.0, 1.0),
        (f64::NAN, 400.0, 1.0),
        (400.0, f64::INFINITY, 1.0),
        (400.0, 400.0, -1.0),
        (400.0, 400.0, f64::INFINITY),
    ] {
        assert!(accepted.advance(start, end, duration).is_err());
        assert_eq!(accepted.dose(), 0.125);
    }
}

#[test]
fn invalid_calibration_cutoff_and_controls_are_rejected() {
    assert!(
        TemperatureShift::new(
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 2.0,
                c2_k: 100.0
            },
            400.0,
            [300.0, 500.0],
        )
        .is_err()
    );
    for cutoff in [Some(329.0), Some(400.0), Some(f64::NAN)] {
        assert!(ShiftedMobility::new("bad".into(), wlf(), 0.25, cutoff, controls()).is_err());
    }
    for integration in [
        MobilityIntegration {
            absolute_dose_tolerance: 0.0,
            ..controls()
        },
        MobilityIntegration {
            relative_dose_tolerance: -1.0,
            ..controls()
        },
        MobilityIntegration {
            relative_dose_tolerance: f64::NAN,
            ..controls()
        },
        MobilityIntegration {
            maximum_evaluations: 4,
            ..controls()
        },
        MobilityIntegration {
            maximum_evaluations: 16_386,
            ..controls()
        },
        MobilityIntegration {
            maximum_depth: 33,
            ..controls()
        },
    ] {
        assert!(ShiftedMobility::new("bad".into(), wlf(), 0.25, None, integration).is_err());
    }
    for rate in [0.0, -1.0, f64::INFINITY] {
        assert!(ShiftedMobility::new("bad".into(), wlf(), rate, None, controls()).is_err());
    }
    assert!(ShiftedMobility::new("".into(), wlf(), 0.25, None, controls()).is_err());
}
