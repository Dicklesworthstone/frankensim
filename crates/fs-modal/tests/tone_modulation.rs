//! Prescribed source sweeps: phase, force reconstruction and physical response.
use fs_math::c64::C64;
use fs_modal::time::{
    HarmonicForceTransition, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace, ModalToneBudget, PreparedModalTones,
    harmonic::ToneModulationStep,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn phase_integral_and_same_basis_endpoints_cover_signed_gain_and_dc() {
    for (h, carrier, frequency, gain) in [
        (0.13, 11.0, [7.0, 19.0], [-0.4, 1.2]),
        (0.17, 0.0, [0.0, 17.0], [1.0, -1.0]),
        (0.19, 13.0, [19.0, 7.0], [0.0, 0.0]),
        (0.23, 0.0, [0.0, 0.0], [-0.3, 0.9]),
    ] {
        let step = ToneModulationStep::try_new(h, carrier, frequency, gain).unwrap();
        let (c, s) = (1.3, -0.7);
        let [c0, s0, c1, s1] = step.endpoints(c, s).unwrap();
        near(c0, c * gain[0], 1e-15);
        near(s0, s * gain[0], 1e-15);
        let phase = h * (frequency[0] + frequency[1]) / 2.0;
        let (ps, pc) = phase.sin_cos();
        let (nc, ns) = step.advance_force_phase(c, s).unwrap();
        near(nc, c * pc + s * ps, 3e-15);
        near(ns, s * pc - c * ps, 3e-15);
        let (cs, cc) = (h * carrier).sin_cos();
        near(c1 * cc + s1 * cs, gain[1] * nc, 3e-15);

        // Check the complete complex envelope, a stronger check than one
        // scalar force projection at a potentially favourable carrier phase.
        for index in 0..=64 {
            let t = index as f64 / 64.0;
            let exact_gain = gain[0] * (1.0 - t) + gain[1] * t;
            let relative =
                h * ((frequency[0] - carrier) * t + 0.5 * (frequency[1] - frequency[0]) * t * t);
            let (rs, rc) = relative.sin_cos();
            let exact_c = exact_gain * (c * rc + s * rs);
            let exact_s = exact_gain * (s * rc - c * rs);
            let approximate_c = c0 * (1.0 - t) + c1 * t;
            let approximate_s = s0 * (1.0 - t) + s1 * t;
            let error = (approximate_c - exact_c).hypot(approximate_s - exact_s);
            let bound = c.hypot(s) * step.coefficient_error_bound();
            assert!(error <= bound + 3e-15, "{error} > {bound}");
        }
    }
}

#[test]
fn zero_detuning_retains_exact_affine_coefficients_and_silent_phase() {
    let h = 1.0 / 48_000.0;
    let w = 817.3;
    let step = ToneModulationStep::try_new(h, w, [w, w], [-0.2, 0.7]).unwrap();
    assert_eq!(step.coefficient_error_bound().to_bits(), 0.0_f64.to_bits());
    let (c, s) = (0.3, -0.8);
    assert_eq!(
        step.endpoints(c, s).unwrap(),
        [c * -0.2, s * -0.2, c * 0.7, s * 0.7]
    );
    let original = HarmonicForceTransition::try_new(1000.0, 0.01, w, h).unwrap();
    assert_eq!(
        step.advance_force_phase(c, s).unwrap(),
        original.advance_force_phase(c, s).unwrap()
    );

    let mut phase = (c, s);
    let mut integrated_phase = 0.0;
    for index in 0..1000 {
        let omega0 = w + 400.0 * index as f64 / 1000.0;
        let omega1 = w + 400.0 * (index + 1) as f64 / 1000.0;
        let silent = ToneModulationStep::try_new(h, w, [omega0, omega1], [0.0, 0.0]).unwrap();
        assert!(
            silent
                .endpoints(phase.0, phase.1)
                .unwrap()
                .iter()
                .all(|x| *x == 0.0)
        );
        assert_eq!(silent.coefficient_error_bound(), 0.0);
        phase = silent.advance_force_phase(phase.0, phase.1).unwrap();
        integrated_phase += h * (omega0 + omega1) / 2.0;
    }
    let (ps, pc) = integrated_phase.sin_cos();
    near(phase.0, c * pc + s * ps, 3e-13);
    near(phase.1, s * pc - c * ps, 3e-13);
    let attack = ToneModulationStep::try_new(h, w, [w + 400.0; 2], [0.0, 1.0]).unwrap();
    let endpoints = attack.endpoints(phase.0, phase.1).unwrap();
    assert!(endpoints[2].hypot(endpoints[3]) > 0.8);
}

fn model() -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 11.0,
            damping_ratio: 0.06,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 10.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 100.0,
            maximum_total_energy_j: 1e6,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap()
}

const DURATION: f64 = 0.7;
const INITIAL: [f64; 4] = [0.013, -0.02, 0.0, 0.0];

// Independent continuous physical ODE. Its force phase is evaluated directly
// from the analytic chirp at absolute time, never using modulation endpoints.
fn reference() -> [f64; 4] {
    let derivative = |t: f64, x: [f64; 4]| {
        let gain = -0.2 + 1.3 * t / DURATION;
        let phase = 5.0 * t + 6.0 * t * t / DURATION;
        let force = gain * (1.3 * phase.cos() - 0.4 * phase.sin());
        [
            x[1],
            force - 121.0 * x[0] - 1.32 * x[1],
            force * x[1],
            1.32 * x[1] * x[1],
        ]
    };
    let mut x = INITIAL;
    let h = DURATION / 16_384.0;
    for index in 0..16_384 {
        let t = index as f64 * h;
        let a = derivative(t, x);
        let b = derivative(t + h * 0.5, core::array::from_fn(|j| x[j] + h * a[j] * 0.5));
        let c = derivative(t + h * 0.5, core::array::from_fn(|j| x[j] + h * b[j] * 0.5));
        let d = derivative(t + h, core::array::from_fn(|j| x[j] + h * c[j]));
        for j in 0..4 {
            x[j] += h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}

fn render(steps: usize) -> [f64; 4] {
    let mut source = model();
    source
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: INITIAL[0],
            velocity_m_sqrt_kg_per_s: INITIAL[1],
        }])
        .unwrap();
    let h = DURATION / steps as f64;
    let prepared = PreparedModalTones::try_new_affine(
        &source,
        &[11.0],
        h,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&source);
    let mut phase = (1.3, -0.4);
    let (mut work, mut heat) = (0.0, 0.0);
    for index in 0..steps {
        let [t0, t1] = [
            index as f64 / steps as f64,
            (index + 1) as f64 / steps as f64,
        ];
        let proposal = ToneModulationStep::try_new(
            h,
            11.0,
            [5.0 + 12.0 * t0, 5.0 + 12.0 * t1],
            [-0.2 + 1.3 * t0, -0.2 + 1.3 * t1],
        )
        .unwrap();
        let [c0, s0, c1, s1] = proposal.endpoints(phase.0, phase.1).unwrap();
        let frame = source
            .step_tones_affine_into(&prepared, &[c0], &[s0], &[c1], &[s1], &mut workspace)
            .unwrap();
        work += frame.input_work_j;
        heat += frame.viscous_dissipation_j;
        assert!(frame.viscous_dissipation_j >= 0.0);
        phase = proposal.advance_force_phase(phase.0, phase.1).unwrap();
    }
    let q = source.states()[0].displacement_m_sqrt_kg;
    let v = source.states()[0].velocity_m_sqrt_kg_per_s;
    let energy_change =
        0.5 * (v * v + 121.0 * q * q - INITIAL[1] * INITIAL[1] - 121.0 * INITIAL[0] * INITIAL[0]);
    near(energy_change, work - heat, 2e-13);
    [q, v, work, heat]
}

#[test]
fn swept_force_crossing_resonance_converges_to_independent_motion_work_and_heat() {
    let expected = reference();
    let results = [render(32), render(64), render(128)];
    for component in 0..4 {
        let errors = results.map(|result| (result[component] - expected[component]).abs());
        eprintln!(
            "chirp component={component} steps=32,64,128 errors={errors:?} ratios={:.8},{:.8}",
            errors[0] / errors[1],
            errors[1] / errors[2]
        );
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!(
                (3.6..4.4).contains(&ratio),
                "component {component}: errors {errors:?}, ratio {ratio}"
            );
        }
        assert!(errors[2] < 3e-5, "component {component}: {errors:?}");
    }
}

#[test]
fn invalid_or_unrepresentable_proposals_refuse_without_changing_the_input() {
    for (h, carrier, frequency, gain) in [
        (0.0, 1.0, [1.0; 2], [1.0; 2]),
        (-1.0, 1.0, [1.0; 2], [1.0; 2]),
        (f64::NAN, 1.0, [1.0; 2], [1.0; 2]),
        (1.0, -1.0, [1.0; 2], [1.0; 2]),
        (1.0, 1.0, [-1.0, 1.0], [1.0; 2]),
        (1.0, 1.0, [1.0, f64::INFINITY], [1.0; 2]),
        (1.0, 1.0, [1.0; 2], [f64::NAN, 1.0]),
        (f64::MAX, 1.0, [f64::MAX; 2], [1.0; 2]),
        (1.0, 0.0, [1e100; 2], [1e308; 2]),
    ] {
        assert!(ToneModulationStep::try_new(h, carrier, frequency, gain).is_err());
    }
    let step = ToneModulationStep::try_new(0.1, 11.0, [7.0, 19.0], [0.2, 1.2]).unwrap();
    let saved = step;
    assert!(step.endpoints(f64::NAN, 0.0).is_err());
    assert!(step.advance_force_phase(0.0, f64::INFINITY).is_err());
    assert!(step.endpoints(f64::MAX, f64::MAX).is_err());
    assert_eq!(step, saved);

    // Neither a large gain nor small dimensionless detuning alone justifies
    // overflowing or underflowing an otherwise representable bound.
    let scaled = ToneModulationStep::try_new(1.0, 0.0, [1e-200; 2], [1e308; 2]).unwrap();
    near(scaled.coefficient_error_bound(), 1.25e-93, 1e-108);
    let cancelling =
        ToneModulationStep::try_new(1.0, 0.0, [1e-308; 2], [-f64::MAX, f64::MAX]).unwrap();
    assert!(cancelling.coefficient_error_bound().is_finite());
    assert!(cancelling.coefficient_error_bound() > 0.8);
}
