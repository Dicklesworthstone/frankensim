//! Full-circle impedance admission of the actual stored digital coefficients.
//! The oracles are elementary real-part formulas and a discrete storage law.

use fs_vfit::discretize::{Biquad, DigitalFilter, check_digital_passivity};
use fs_vfit::passivity::{
    GlobalPassivityBudget, GlobalPassivityUnresolved, GlobalPassivityVerdict,
};

fn filter(direct: f64, sections: Vec<Biquad>) -> DigitalFilter {
    DigitalFilter {
        sections,
        direct,
        t_s: 1.0 / 48_000.0,
        prewarp: 0.0,
    }
}

fn assert_certified(filter: &DigitalFilter) {
    let report = check_digital_passivity(filter, GlobalPassivityBudget::default()).unwrap();
    assert_eq!(
        report.verdict,
        GlobalPassivityVerdict::Certified,
        "{report:?}"
    );
}

fn differentiator(gain: f64) -> Biquad {
    Biquad {
        b: [gain, -gain, 0.0],
        a: [1.0, 0.0],
    }
}

fn integrator(gain: f64) -> Biquad {
    Biquad {
        b: [gain, gain, 0.0],
        a: [-1.0, 0.0],
    }
}

#[test]
fn g1_global_admission_keeps_compensating_terms_and_impedance_units() {
    for gain in [1.0e-250, 1.0, 1.0e250] {
        // Total impedance is gain * (2 - .25*(1+z^-1)/(1-.5*z^-1)).
        // Its real part ranges from gain at DC to 2*gain at Nyquist.
        // Both the direct slot and one complete section are negative.
        let model = filter(
            -gain,
            vec![
                Biquad {
                    b: [3.0 * gain, 0.0, 0.0],
                    a: [0.0, 0.0],
                },
                Biquad {
                    b: [-0.25 * gain, -0.25 * gain, 0.0],
                    a: [-0.5, 0.0],
                },
            ],
        );
        assert_certified(&model);
    }
    assert_certified(&filter(0.0, vec![]));
    assert_certified(&filter(
        -1.0,
        vec![Biquad {
            b: [1.0, 0.0, 0.0],
            a: [0.0, 0.0],
        }],
    ));
    assert_eq!(
        check_digital_passivity(&filter(-1.0, vec![]), GlobalPassivityBudget::default())
            .unwrap()
            .verdict,
        GlobalPassivityVerdict::Violated,
    );
}

#[test]
fn g1_passive_relaxation_and_bandpass_keep_exact_lossless_endpoints() {
    // Z(s)=1/(s+1), with K=96000: positive real on the complete circle,
    // including its exactly zero Nyquist value. No endpoint floor is needed.
    let k = 96_000.0;
    assert_certified(&filter(
        0.0,
        vec![Biquad {
            b: [1.0 / (k + 1.0), 1.0 / (k + 1.0), 0.0],
            a: [(1.0 - k) / (1.0 + k), 0.0],
        }],
    ));

    // Re Z = 2*g*(1-a2)*sin(theta)^2 / |D|^2 for b=[g,0,-g].
    // Strict Jury stability makes the denominator positive; both endpoints
    // are exactly lossless, while the interior has strictly positive loss.
    let a2 = 0.9801;
    assert_certified(&filter(
        0.0,
        vec![Biquad {
            b: [0.125, 0.0, -0.125],
            a: [-1.8, a2],
        }],
    ));
}

#[test]
fn g1_lossless_inertance_uses_summed_boundary_storage_not_a_frequency_grid() {
    assert_certified(&filter(
        0.25,
        vec![differentiator(-3.0), differentiator(4.0)],
    ));
    assert_certified(&filter(
        0.0,
        vec![differentiator(-1.0), differentiator(1.0)],
    ));
    // A tiny NEGATIVE remainder is still active. Ordinary summation of
    // these gains loses -1 between the two large opposite numbers.
    let active = filter(
        2.0,
        vec![
            differentiator(1.0e16),
            differentiator(-1.0),
            differentiator(-1.0e16),
        ],
    );
    let report = check_digital_passivity(&active, GlobalPassivityBudget::default()).unwrap();
    assert_eq!(report.verdict, GlobalPassivityVerdict::Violated);
    assert!(report.violation_frequency_rad_s.is_none());

    // Independent storage identity for a positive inertance. For zero
    // prehistory, its boundary coordinate is the alternating sum of flows.
    let model = filter(0.0, vec![differentiator(1.0)]);
    assert_certified(&model);
    let mut state = model.zero_state();
    let flows = [1.0, -0.5, 0.25, 0.0, -0.125, 0.75, 0.0];
    let mut work = 0.0;
    for (n, &flow) in flows.iter().enumerate() {
        work += model.t_s * model.step(&mut state, flow).unwrap() * flow;
        let coordinate: f64 = flows[..=n]
            .iter()
            .enumerate()
            .map(|(j, u)| if (n - j) % 2 == 0 { *u } else { -*u })
            .sum();
        let storage = model.t_s * coordinate * coordinate;
        assert!((work - storage).abs() <= 2.0e-14 * model.t_s);
    }
}

#[test]
fn g1_lossless_compliance_has_independent_positive_dc_storage() {
    for gain in [1.0e-250, 1.0, 1.0e250] {
        let model = filter(0.0, vec![integrator(gain)]);
        assert!(model.is_stable());
        assert_certified(&model);
    }
    assert_certified(&filter(
        0.25,
        vec![integrator(-3.0), integrator(4.0), differentiator(2.0)],
    ));
    assert_certified(&filter(0.0, vec![integrator(-1.0), integrator(1.0)]));
    for sections in [
        // Positive inertance cannot hide negative compliance, or vice versa.
        vec![integrator(-1.0), differentiator(100.0)],
        vec![integrator(100.0), differentiator(-1.0)],
        // Retain the small negative DC residue lost by ordinary summation.
        vec![integrator(1.0e16), integrator(-1.0), integrator(-1.0e16)],
    ] {
        let report =
            check_digital_passivity(&filter(2.0, sections), GlobalPassivityBudget::default())
                .unwrap();
        assert_eq!(report.verdict, GlobalPassivityVerdict::Violated);
        assert!(report.violation_frequency_rad_s.is_none());
    }

    // Independent compliance work identity: the charge/displaced-volume
    // coordinate is the cumulative flow, without an alternating sign.
    let gain = 0.75;
    let model = filter(0.0, vec![integrator(gain)]);
    let mut state = model.zero_state();
    let (mut coordinate, mut work) = (0.0, 0.0);
    for flow in [1.0, -0.5, 0.25, 0.0, -0.125, 0.75, 0.0] {
        work += model.t_s * model.step(&mut state, flow).unwrap() * flow;
        coordinate += flow;
        let storage = model.t_s * gain * coordinate * coordinate;
        assert!((work - storage).abs() <= 2.0e-14 * model.t_s);
    }
}

#[test]
fn g0_noncanonical_dc_poles_are_not_silently_repaired() {
    // This change admits only the exact Tustin storage identity. Neither an
    // arbitrary numerator nor a pole outside the circle gets that exception.
    for section in [
        Biquad {
            b: [1.0, 0.0, 0.0],
            a: [-1.0, 0.0],
        },
        Biquad {
            b: [1.0, 1.0_f64.next_down(), 0.0],
            a: [-1.0, 0.0],
        },
        Biquad {
            b: [1.0, 1.0, 0.0],
            a: [-1.0_f64.next_up(), 0.0],
        },
    ] {
        let report = check_digital_passivity(
            &filter(1.0, vec![section]),
            GlobalPassivityBudget::default(),
        )
        .unwrap();
        assert_ne!(report.verdict, GlobalPassivityVerdict::Certified);
    }
}

#[test]
fn g1_narrow_active_band_is_found_between_passive_endpoints() {
    let center = -0.75_f64;
    let width = 1.0e-4_f64;
    // Exact reference: Re Z = (cos(theta)-center)^2-width^2.
    // Both endpoints and the instantaneous coefficient are positive.
    let model = filter(
        0.0,
        vec![Biquad {
            b: [center * center - width * width + 0.5, -2.0 * center, 0.5],
            a: [0.0, 0.0],
        }],
    );
    assert!((1.0 - center).powi(2) > width * width);
    assert!((-1.0 - center).powi(2) > width * width);
    assert!(model.sections[0].b[0] > 0.0);
    let report = check_digital_passivity(&model, GlobalPassivityBudget::default()).unwrap();
    assert_eq!(
        report.verdict,
        GlobalPassivityVerdict::Violated,
        "{report:?}"
    );
    if let Some(omega) = report.violation_frequency_rad_s {
        let theta = omega * model.t_s;
        assert!((theta.cos() - center).powi(2) < width * width, "{theta:e}");
    }

    let mut budget = GlobalPassivityBudget::default();
    budget.max_subdivisions = 0;
    let bounded = check_digital_passivity(&model, budget).unwrap();
    assert_ne!(bounded.verdict, GlobalPassivityVerdict::Certified);
    assert_eq!(bounded.subdivisions, 0);
    if bounded.verdict == GlobalPassivityVerdict::Unresolved {
        assert!(bounded.unresolved.is_some());
    }
}

#[test]
fn g0_continuous_passivity_does_not_certify_rounded_digital_coefficients() {
    // These are the emitted binary64 coefficients for K=96000 and
    // Z(s)=(s+2)/((s+2)^2+1000^2), a sum of two conjugate stable poles
    // with positive real residues. The continuous real part is positive.
    //
    // EXACT dyadic arithmetic on the stored coefficients gives
    // b0-b1+b2 = -20646 * 2^-84 < 0, and 1-a1+a2 > 0.
    // The actual filter therefore has negative Nyquist resistance.
    let model = filter(
        0.0,
        vec![Biquad {
            b: [
                0.000010415319569884331,
                4.3396260785751677e-10,
                -0.000010414885607276475,
            ],
            a: [-1.9994827148355834, 0.9999166791792914],
        }],
    );
    let report = check_digital_passivity(&model, GlobalPassivityBudget::default()).unwrap();
    assert_eq!(
        report.verdict,
        GlobalPassivityVerdict::Violated,
        "{report:?}"
    );
    let theta = report.violation_frequency_rad_s.unwrap() * model.t_s;
    assert!((theta - core::f64::consts::PI).abs() < 1.0e-12);
}

#[test]
fn g0_invalid_coefficients_unestablished_stability_and_work_limits_cannot_pass() {
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut model = filter(1.0, vec![]);
        model.t_s = duration;
        assert!(check_digital_passivity(&model, GlobalPassivityBudget::default()).is_err());
    }
    let malformed = filter(
        0.0,
        vec![Biquad {
            b: [f64::NAN, 0.0, 0.0],
            a: [0.0, 0.0],
        }],
    );
    assert!(check_digital_passivity(&malformed, GlobalPassivityBudget::default()).is_err());
    let unstable = filter(
        1.0,
        vec![Biquad {
            b: [1.0, 0.0, 0.0],
            a: [-1.25, 0.0],
        }],
    );
    assert_ne!(
        check_digital_passivity(&unstable, GlobalPassivityBudget::default())
            .unwrap()
            .verdict,
        GlobalPassivityVerdict::Certified,
    );
    let mut budget = GlobalPassivityBudget::default();
    budget.max_order = 1;
    let second_order = filter(
        0.0,
        vec![Biquad {
            b: [0.125, 0.0, -0.125],
            a: [-1.8, 0.9801],
        }],
    );
    let bounded = check_digital_passivity(&second_order, budget).unwrap();
    assert_eq!(bounded.verdict, GlobalPassivityVerdict::Unresolved);
    assert_eq!(
        bounded.unresolved,
        Some(GlobalPassivityUnresolved::OrderLimit)
    );
}
