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

fn joint_refinement(
    states: usize,
    maximum_harmonics: usize,
) -> fs_orbit::sideband::sensitivity::refinement::SidebandCarrierRefinementControls {
    use fs_orbit::sideband::sensitivity::refinement::{
        SidebandCarrierRefinementControls, SidebandCarrierRefinementTolerances,
    };
    SidebandCarrierRefinementControls {
        maximum_harmonics,
        tolerances: SidebandCarrierRefinementTolerances {
            equation_tolerances: vec![2e-9; states],
            state_tolerances: vec![2e-9; states],
            derivative_equation_tolerances: vec![2e-9; states],
            derivative_state_tolerances: vec![2e-9; states],
        },
    }
}

#[test]
fn adaptive_velocity_peak_and_width_require_complete_rounds_and_share_all_query_budgets() {
    use fs_couple::phs_orbit::sideband::PhsPortAdaptiveResonanceRequest;
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING).unwrap();
    let base = orbit(OMEGA, vec![vec![C64::ZERO]; 2]);
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO]]).unwrap();
    let mut c = physical_controls();
    c.sideband.maximum_unknowns = 30;
    c.max_hessian_actions = 2_000_000;
    c.max_observation_work = 10_000_000_000;
    let force = C64::new(0.8, 0.6);
    let forcing = [C64::ZERO, C64::ZERO, force];
    let tolerances = [2e-9];
    let refinement = joint_refinement(2, 7);
    let mut fixed = request(&forcing, &tolerances);
    fixed.controls.maximum_evaluations = 512;
    let request = PhsPortAdaptiveResonanceRequest {
        resonance: fixed,
        refinement: &refinement,
        port_change_tolerances: &tolerances,
        derivative_port_change_tolerances: &tolerances,
        absolute_angular_change_tolerance_rad_s: 1e-7,
        relative_angular_change_tolerance: 1e-9,
    };
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    assert!(matches!(
        source.refine_port_resonance_adaptive(&request, &scales(2), || true),
        Err(ResonanceError::Cancelled)
    ));
    let resolved = source
        .refine_port_resonance_adaptive(&request, &scales(2), || false)
        .unwrap();
    assert_eq!(resolved.state_harmonics, 7);
    assert_eq!(
        resolved
            .refinement_trace
            .iter()
            .map(|c| (c.lower_harmonics, c.harmonics, c.completed))
            .collect::<Vec<_>>(),
        [(1, 3, true), (3, 7, true)]
    );
    assert!(resolved.refinement_trace[0].peak_change_ratio.is_none());
    assert!(resolved.refinement_trace[1].peak_change_ratio.unwrap() <= 1.0);
    assert!(resolved.refinement_trace[1].bandwidth_change_ratio.unwrap() <= 1.0);
    assert_eq!(
        resolved.evaluations,
        resolved
            .refinement_trace
            .iter()
            .map(|c| c.evaluations)
            .sum::<usize>()
    );
    assert_eq!(resolved.pair_solves, 2 * resolved.evaluations);
    assert_eq!(resolved.total_hessian_actions, 3584 * resolved.evaluations);
    close(resolved.response.peak.angular_frequency_rad_s, 2.0);
    let bandwidth = resolved.response.bandwidth.as_ref().unwrap();
    close(bandwidth.angular_bandwidth_rad_s, 1.5);
    close(bandwidth.response_quality_factor, 4.0 / 3.0);
    let w = resolved.response.peak.angular_frequency_rad_s;
    let denominator = C64::new(STIFFNESS - MASS * w * w, DAMPING * w);
    near(
        resolved.response.peak.sample.amplitude,
        (C64::new(0.0, w) * force / denominator).conj(),
        2e-9,
    );
    near(
        resolved.response.peak.sample.angular_derivative,
        (C64::new(0.0, STIFFNESS + MASS * w * w) * force / (denominator * denominator)).conj(),
        2e-9,
    );

    let short_refinement = joint_refinement(2, 3);
    let mut short = request.clone();
    short.refinement = &short_refinement;
    assert!(
        matches!(source.refine_port_resonance_adaptive(&short, &scales(2), || false),
        Err(ResonanceError::Evaluation(PhsPortSidebandError::ResonanceUnresolved { last }))
            if last.completed && last.peak_change_ratio.is_none())
    );
    // Two remaining evaluations cannot start another peak search. This must
    // stay a TOTAL campaign budget error, not become invalid local controls.
    let mut limited = request.clone();
    limited.resonance.controls.maximum_evaluations = resolved.refinement_trace[0].evaluations + 2;
    assert!(
        matches!(source.refine_port_resonance_adaptive(&limited, &scales(2), || false),
        Err(ResonanceError::EvaluationBudget { maximum_evaluations })
            if maximum_evaluations == limited.resonance.controls.maximum_evaluations)
    );
    c.max_hessian_actions = resolved.total_hessian_actions - 1;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    assert!(matches!(
        source.refine_port_resonance_adaptive(&request, &scales(2), || false),
        Err(ResonanceError::Evaluation(PhsPortSidebandError::Budget(_)))
    ));

    // The same smooth curve and same midpoint estimates at two orders do NOT
    // establish fine peak/linewidth location if the root brackets remain wide.
    c.max_hessian_actions = 2_000_000;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let mut coarse = request;
    coarse.resonance.controls.absolute_angular_tolerance_rad_s = 0.02;
    coarse.resonance.controls.relative_angular_tolerance = 0.0;
    coarse.absolute_angular_change_tolerance_rad_s = 1e-8;
    coarse.relative_angular_change_tolerance = 0.0;
    assert!(
        matches!(source.refine_port_resonance_adaptive(&coarse, &scales(2), || false),
        Err(ResonanceError::Evaluation(PhsPortSidebandError::ResonanceUnresolved { last }))
            if last.completed && (last.peak_change_ratio.unwrap() > 1.0
                || last.bandwidth_change_ratio.unwrap() > 1.0))
    );
}

#[test]
fn adaptive_periodic_stiffness_restarts_common_order_and_matches_continued_fraction_peak() {
    use fs_couple::phs_orbit::sideband::PhsPortAdaptiveResonanceRequest;
    use fs_modal::harmonic::resonance::{
        HarmonicTransferSample, refine_transfer_extremum_with_cancel,
    };
    use fs_nlmodal::{SosModalStorage, StressChannel};
    use fs_phs::PortHamiltonian;
    let (stiffness, damping, cubic, amplitude) = (4.0, 0.4, 0.1, 0.5);
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 1.0],
        Box::new(SosModalStorage {
            omegas: vec![2.0],
            channels: vec![StressChannel {
                coefficient: cubic,
                coupling: vec![1.0],
            }],
        }),
    )
    .unwrap();
    let q1 = C64::from_re(0.5 * amplitude);
    let base = orbit(
        OMEGA,
        vec![
            vec![C64::ZERO, q1],
            vec![C64::ZERO, C64::new(0.0, OMEGA) * q1],
        ],
    );
    // Exact Duffing sustaining force for q=A*cos(Omega*t), including the third
    // harmonic. The small-probe stiffness then mixes sidebands separated by 2.
    let drive1 = C64::new(
        stiffness - OMEGA * OMEGA + 0.75 * cubic * amplitude * amplitude,
        damping * OMEGA,
    ) * q1;
    let drive3 = (q1 * q1 * q1).scale(cubic);
    let input =
        PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, drive1, C64::ZERO, drive3]]).unwrap();
    let mut c = controls(2, 1, 1, 7);
    c.sideband.maximum_unknowns = 62;
    c.max_hessian_actions = 2_000_000;
    c.max_observation_work = 10_000_000_000;
    let force = C64::new(0.8, 0.6);
    let forcing = [C64::ZERO, C64::ZERO, force];
    let tolerances = [2e-10];
    let changes = [2e-7];
    let mut refinement = joint_refinement(2, 15);
    refinement.tolerances.state_tolerances.fill(2e-7);
    refinement.tolerances.derivative_state_tolerances.fill(2e-7);
    let mut fixed = request(&forcing, &tolerances);
    fixed.peak_bracket_rad_s = (1.5, 2.5);
    fixed.half_power_brackets_rad_s = None;
    fixed.controls.maximum_evaluations = 256;
    let request = PhsPortAdaptiveResonanceRequest {
        resonance: fixed,
        refinement: &refinement,
        port_change_tolerances: &changes,
        derivative_port_change_tolerances: &changes,
        absolute_angular_change_tolerance_rad_s: 1e-7,
        relative_angular_change_tolerance: 1e-9,
    };
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let resolved = source
        .refine_port_resonance_adaptive(&request, &scales(2), || false)
        .unwrap();
    assert_eq!(
        resolved
            .refinement_trace
            .iter()
            .map(|c| (c.harmonics, c.completed))
            .collect::<Vec<_>>(),
        [(3, false), (7, true), (15, true)]
    );
    assert_eq!(resolved.refinement_trace[0].evaluations, 1);
    assert!(
        resolved.refinement_trace[0].maximum_equation_ratio > 1.0
            || resolved.refinement_trace[0].maximum_derivative_equation_ratio > 1.0
            || resolved.refinement_trace[0].maximum_response_change_ratio > 1.0
            || resolved.refinement_trace[0].maximum_derivative_change_ratio > 1.0
    );
    assert_eq!(resolved.pair_solves, 2 * resolved.evaluations);

    // Independent scalar Schur complement of the exact tridiagonal Mathieu
    // frequency equations, with ANALYTIC derivatives of each continued fraction.
    // Twelve eliminated pairs are far beyond the physically relevant bands here;
    // this is a high-accuracy reference, not an asserted infinite-tail enclosure.
    let mean_stiffness = stiffness + 1.5 * cubic * amplitude * amplitude;
    let coupling_squared = (0.75 * cubic * amplitude * amplitude).powi(2);
    let reference = |w: f64| {
        let mut denominator = C64::new(mean_stiffness - w * w, damping * w);
        let mut slope = C64::new(-2.0 * w, damping);
        for sign in [-1.0, 1.0] {
            let (mut tail, mut tail_slope) = (C64::ZERO, C64::ZERO);
            for j in (1..=12).rev() {
                let frequency = w + sign * 2.0 * OMEGA * j as f64;
                let d =
                    C64::new(mean_stiffness - frequency * frequency, damping * frequency) - tail;
                let ds = C64::new(-2.0 * frequency, damping) - tail_slope;
                tail = C64::from_re(coupling_squared) / d;
                tail_slope = (ds / (d * d)).scale(-coupling_squared);
            }
            denominator = denominator - tail;
            slope = slope - tail_slope;
        }
        let q = force / denominator;
        HarmonicTransferSample {
            amplitude: (C64::new(0.0, w) * q).conj(),
            angular_derivative: (C64::new(0.0, 1.0) * q
                - C64::new(0.0, w) * q * slope / denominator)
                .conj(),
        }
    };
    let exact = refine_transfer_extremum_with_cancel(
        request.resonance.peak_bracket_rad_s,
        request.resonance.controls,
        || false,
        |w| Ok::<_, ()>(reference(w)),
    )
    .unwrap();
    close(
        resolved.response.peak.angular_frequency_rad_s,
        exact.angular_frequency_rad_s,
    );
    let expected = reference(resolved.response.peak.angular_frequency_rad_s);
    near(
        resolved.response.peak.sample.amplitude,
        expected.amplitude,
        2e-10,
    );
    near(
        resolved.response.peak.sample.angular_derivative,
        expected.angular_derivative,
        2e-10,
    );
}
