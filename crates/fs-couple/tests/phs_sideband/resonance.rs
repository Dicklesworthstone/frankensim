use super::{controls, derivative_controls, near, orbit, scales};
use core::cell::Cell;
use fs_couple::phs_orbit::{
    PeriodicPortInput,
    sideband::{
        PhsPortResonanceRequest, PhsPortSidebandControls, PhsPortSidebandError,
        PhsPortSidebandSource,
    },
};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{ResonanceControls, ResonanceError, TransferExtremumKind};

const MASS: f64 = 2.0;
const STIFFNESS: f64 = 8.0;
const DAMPING: f64 = 3.0;
const OMEGA: f64 = 10.0;

fn physical_controls() -> PhsPortSidebandControls {
    let mut c = controls(2, 1, 1, 1);
    c.max_hessian_actions = 4096;
    c
}

fn request<'a>(
    forcing: &'a [C64],
    derivative_tolerances: &'a [f64],
) -> PhsPortResonanceRequest<'a> {
    PhsPortResonanceRequest {
        forcing,
        output_port: 0,
        output_sideband: 1,
        peak_bracket_rad_s: (1.6, 2.4),
        half_power_brackets_rad_s: Some(((0.5, 1.8), (2.2, 4.5))),
        controls: ResonanceControls {
            absolute_angular_tolerance_rad_s: 1e-8,
            relative_angular_tolerance: 1e-10,
            maximum_evaluations: 256,
        },
        derivative_controls: derivative_controls(),
        derivative_envelope_tolerances: derivative_tolerances,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-8 * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn original_velocity_peak_and_bandwidth_use_output_frequency_with_a_negative_carrier() {
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING).unwrap();
    let base = orbit(OMEGA, vec![vec![C64::ZERO]; 2]);
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO]]).unwrap();
    let c = physical_controls();
    let force = C64::new(0.8, 0.6);
    // The probe occupies +1, so an observed resonance at w=2 uses carrier=-8.
    let forcing = [C64::ZERO, C64::ZERO, force];
    let tolerances = [2e-9];
    let request = request(&forcing, &tolerances);
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let result = prepared.refine_port_resonance(&request, || false).unwrap();
    let peak_frequency = (STIFFNESS / MASS).sqrt();
    assert_eq!(result.peak.kind, TransferExtremumKind::Maximum);
    close(result.peak.angular_frequency_rad_s, peak_frequency);
    assert!(result.peak.angular_frequency_rad_s - OMEGA < 0.0);

    let w = result.peak.angular_frequency_rad_s;
    let denominator = C64::new(STIFFNESS - MASS * w * w, DAMPING * w);
    // Original velocity under exp(+i*w*t), followed by the exact convention
    // conversion required by the existing exp(-i*w*t) resonance owner.
    let velocity = C64::new(0.0, w) * force / denominator;
    let derivative = C64::new(0.0, STIFFNESS + MASS * w * w) * force / (denominator * denominator);
    near(result.peak.sample.amplitude, velocity.conj(), 2e-9);
    near(
        result.peak.sample.angular_derivative,
        derivative.conj(),
        2e-9,
    );

    // For the actual velocity response, |K-M*w^2|=C*w gives BOTH half-power
    // crossings exactly. This checks physical output units, not a fitted pole.
    let discriminant = (DAMPING * DAMPING + 4.0 * MASS * STIFFNESS).sqrt();
    let expected_lower = (discriminant - DAMPING) / (2.0 * MASS);
    let expected_upper = (discriminant + DAMPING) / (2.0 * MASS);
    let bandwidth = result.bandwidth.as_ref().unwrap();
    close(bandwidth.lower_angular_frequency_rad_s, expected_lower);
    close(bandwidth.upper_angular_frequency_rad_s, expected_upper);
    close(bandwidth.angular_bandwidth_rad_s, 1.5);
    close(bandwidth.response_quality_factor, 4.0 / 3.0);
    assert_eq!(
        result.evaluations,
        result.peak.evaluations + bandwidth.evaluations
    );
    assert!(result.evaluations <= request.controls.maximum_evaluations);
}

#[test]
fn commensurate_spans_refuse_before_physics_and_peak_and_flanks_share_one_cap() {
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING).unwrap();
    let base = orbit(OMEGA, vec![vec![C64::ZERO]; 2]);
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO]]).unwrap();
    let calls = Cell::new(0_usize);
    let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out[0] = STIFFNESS * direction[0];
        out[1] = direction[1] / MASS;
        true
    };
    let c = physical_controls();
    let forcing = [C64::ZERO, C64::ZERO, C64::new(0.8, 0.6)];
    let tolerances = [2e-9];
    let original = request(&forcing, &tolerances);
    let mut source =
        PhsPortSidebandSource::forced_analytic(&system, &input, &base, &c, &hessian).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let initial_calls = calls.get();

    let mut crosses = original.clone();
    crosses.peak_bracket_rad_s = (4.9, 5.1);
    crosses.half_power_brackets_rad_s = None;
    let mut at_zero = original.clone();
    at_zero.peak_bracket_rad_s = (0.0, 2.4);
    let mut flank_touches = original.clone();
    flank_touches.half_power_brackets_rad_s = Some(((0.5, 1.8), (2.2, 5.0)));
    for invalid in [crosses, at_zero, flank_touches] {
        assert!(matches!(
            prepared.refine_port_resonance(&invalid, || false),
            Err(ResonanceError::Invalid(_))
        ));
        assert_eq!(calls.get(), initial_calls);
    }

    let full = prepared.refine_port_resonance(&original, || false).unwrap();
    let mut limited = original.clone();
    limited.controls.maximum_evaluations = full.evaluations - 1;
    // Each phase individually fits this cap. Only cumulative accounting can
    // refuse: restarting the full allowance for bandwidth would incorrectly pass.
    assert!(limited.controls.maximum_evaluations > full.peak.evaluations);
    assert!(limited.controls.maximum_evaluations > full.bandwidth.as_ref().unwrap().evaluations);
    let before_limited = calls.get();
    assert!(matches!(
        prepared.refine_port_resonance(&limited, || false),
        Err(ResonanceError::EvaluationBudget { maximum_evaluations })
            if maximum_evaluations == limited.controls.maximum_evaluations
    ));
    assert!(calls.get() > before_limited);
}

#[test]
fn physical_refusal_pair_budgets_and_cancellation_never_publish_a_partial_peak() {
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING).unwrap();
    let base = orbit(OMEGA, vec![vec![C64::ZERO]; 2]);
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO]]).unwrap();
    let refused = Cell::new(false);
    let calls = Cell::new(0_usize);
    let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        if refused.get() {
            return false;
        }
        out[0] = STIFFNESS * direction[0];
        out[1] = direction[1] / MASS;
        true
    };
    let c = physical_controls();
    let forcing = [C64::ZERO, C64::ZERO, C64::new(0.8, 0.6)];
    let tolerances = [2e-9];
    let request = request(&forcing, &tolerances);
    {
        let mut source =
            PhsPortSidebandSource::forced_analytic(&system, &input, &base, &c, &hessian).unwrap();
        let mut prepared = source.prepare(&scales(2), || false).unwrap();
        let before_cancel = calls.get();
        assert!(matches!(
            prepared.refine_port_resonance(&request, || true),
            Err(ResonanceError::Cancelled)
        ));
        assert_eq!(calls.get(), before_cancel);
        let mut polls = 0_usize;
        assert!(matches!(
            prepared.refine_port_resonance(&request, || {
                polls += 1;
                polls >= 40
            }),
            Err(ResonanceError::Cancelled)
        ));
        assert!(polls >= 40);
        refused.set(true);
        assert!(matches!(
            prepared.refine_port_resonance(&request, || false),
            Err(ResonanceError::Evaluation(
                PhsPortSidebandError::SourceRefused
            ))
        ));
    }

    // One ordinary physical observation fits 1024 actions, while its paired
    // derivative requires 1792 on these grids. The campaign must retain this
    // original physical budget refusal, not report a valid peak from the primal.
    let mut limited = c;
    limited.max_hessian_actions = 1024;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &limited).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    assert!(matches!(
        prepared.refine_port_resonance(&request, || false),
        Err(ResonanceError::Evaluation(PhsPortSidebandError::Budget(_)))
    ));
}
