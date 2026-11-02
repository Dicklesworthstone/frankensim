use super::super::{BroadbandModalControls, BroadbandModalReport, Delay, ExactFilter};
use super::*;
use crate::broadband_radiation::BroadbandRadiationAuthority;
use fs_math::{c64::C64, det};
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_vfit::{PoleTerm, RationalModel};

const RATE: u32 = 8000;
const FILTER_POLE: f64 = -700.0;
const FILTER_RESIDUE: f64 = 350.0;
const DIRECT: f64 = 0.25;

// An exact supplied real-rational receiver isolates trajectory composition.
// The existing broadband integration suite independently tests vector fitting
// and held-out admission; no fit or natural-frequency response is the oracle.
fn observer(
    modes: Vec<ModalAcousticMode>,
    budget: ModalAcousticTimeBudget,
    delay: usize,
) -> BroadbandModalObserver {
    let model = ModalAcousticTimeModel::try_new_with_free_coordinates(RATE, modes, budget).unwrap();
    let filters = model
        .modes()
        .iter()
        .map(|&mode| {
            ExactFilter::try_new(
                mode,
                RationalModel {
                    terms: vec![PoleTerm::Real {
                        pole: FILTER_POLE,
                        residue: FILTER_RESIDUE,
                    }],
                    d: DIRECT,
                    e: 0.0,
                },
                model.sample_period_s(),
            )
            .unwrap()
        })
        .collect();
    let workspace = ModalAcousticWorkspace::new(&model);
    let candidate_workspace = ModalAcousticWorkspace::new(&model);
    BroadbandModalObserver {
        model,
        workspace,
        candidate_workspace,
        filters,
        constant: None,
        delay: Delay::new(delay as f64, 1).unwrap(),
        report: BroadbandModalReport {
            source_id: "analytic continuous one-pole receiver".into(),
            authority: BroadbandRadiationAuthority::EstimateOnly,
            producer_relative_error_estimate: 0.0,
            controls: BroadbandModalControls {
                fit_order: 1,
                fit_iterations: 1,
                delay_kernel_half_len: 1,
                max_held_out_rel_error: 1e-10,
            },
            sample_rate_hz: RATE,
            held_out_max_rel_error: 0.0,
            delay_samples: delay as f64,
            constant_narrow_band_compatibility: false,
            input_batches: 1,
        },
        forcing_owner: Arc::new(()),
        linear_force_gains: None,
    }
}

fn mode(omega: f64, zeta: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: zeta,
        pressure_per_modal_velocity: C64::ZERO,
    }
}

fn exponential(z: C64, t: f64) -> C64 {
    C64::new(det::cos(z.im * t), det::sin(z.im * t)).scale(det::exp(z.re * t))
}

// Independent partial fractions of the physical damped oscillator followed
// by H(s)=DIRECT+FILTER_RESIDUE/(s-FILTER_POLE). Includes startup and original
// mechanical displacement/velocity with zero initial receiver history.
fn analytical_pressure(
    t: f64,
    mode: ModalAcousticMode,
    omega: f64,
    force: C64,
    q0: f64,
    v0: f64,
) -> f64 {
    let w = mode.angular_frequency_rad_s;
    let gamma = mode.damping_ratio * w;
    let wd = w * det::sqrt(1.0 - mode.damping_ratio * mode.damping_ratio);
    let p1 = C64::new(-gamma, wd);
    let p2 = p1.conj();
    let s = C64::new(0.0, -omega);
    let particular = force * (s * s + s.scale(2.0 * gamma) + C64::from_re(w * w)).recip();
    let first = (C64::from_re(v0) - p2.scale(q0) + (p2 - s) * particular) * (p1 - p2).recip();
    let second = C64::from_re(q0) - particular - first;
    [(s, s * particular), (p1, p1 * first), (p2, p2 * second)]
        .iter()
        .fold(C64::ZERO, |sum, &(pole, velocity)| {
            let e = exponential(pole, t);
            let memory = (e - C64::from_re(det::exp(FILTER_POLE * t)))
                * (pole - C64::from_re(FILTER_POLE)).recip();
            sum + velocity * (e.scale(DIRECT) + memory.scale(FILTER_RESIDUE))
        })
        .re
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn fitted_public_receiver_follows_the_continuous_force_at_its_actual_frequency() {
    let physical_mode = mode(1200.0, 0.02);
    let source = ModalAcousticTimeModel::try_new(
        RATE,
        vec![physical_mode],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let transfer =
        |w: f64| C64::from_re(DIRECT) + C64::new(-FILTER_POLE, -w).recip().scale(FILTER_RESIDUE);
    let omega: Vec<_> = (1..=12).map(|i| 200.0 * f64::from(i)).collect();
    let samples = super::super::ModalTransferSamples {
        source_id: "analytic one-pole public receiver".into(),
        per_mode: vec![omega.iter().map(|&w| transfer(w)).collect()],
        omega_rad_s: omega,
        held_out: (1..=11)
            .map(|i| {
                let w = 200.0 * f64::from(i) + 50.0;
                (w, vec![transfer(w)])
            })
            .collect(),
        propagation_delay_s: 0.0,
        producer_relative_error_estimate: 0.0,
    };
    let (mut observed, report) = BroadbandModalObserver::build(
        &source,
        &samples,
        BroadbandModalControls {
            fit_order: 1,
            fit_iterations: 8,
            delay_kernel_half_len: 1,
            max_held_out_rel_error: 1e-8,
        },
        RATE,
    )
    .unwrap();
    assert!(!report.constant_narrow_band_compatibility);
    assert!(report.held_out_max_rel_error < 1e-8);
    let prepared = observed.prepare_harmonic(740.3).unwrap();
    let (mut cosine, mut sine) = (1.7, -0.8);
    for sample in 1..=128 {
        close(
            observed
                .step_harmonic(&prepared, &[cosine], &[sine])
                .unwrap(),
            analytical_pressure(
                sample as f64 / f64::from(RATE),
                physical_mode,
                740.3,
                C64::new(1.7, -0.8),
                0.0,
                0.0,
            ),
        );
        (cosine, sine) = prepared.advance_force_phase(cosine, sine).unwrap();
    }
}

#[test]
fn continuous_harmonic_receiver_matches_transient_convolution_and_mechanical_bits() {
    let physical_mode = mode(1200.0, 0.02);
    let mut observed = observer(
        vec![physical_mode],
        ModalAcousticTimeBudget::audible_reference(),
        7,
    );
    observed
        .model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0002,
            velocity_m_sqrt_kg_per_s: -0.03,
        }])
        .unwrap();
    let prepared = observed.prepare_harmonic(740.3).unwrap();
    let mut original = observed.model.clone();
    let mut workspace = ModalAcousticWorkspace::new(&original);
    let (mut cosine, mut sine) = (1.7, -0.8);
    for sample in 1..=512 {
        let pressure = observed
            .step_harmonic(&prepared, &[cosine], &[sine])
            .unwrap();
        let expected = if sample <= 7 {
            0.0
        } else {
            analytical_pressure(
                (sample - 7) as f64 / f64::from(RATE),
                physical_mode,
                740.3,
                C64::new(1.7, -0.8),
                0.0002,
                -0.03,
            )
        };
        close(pressure, expected);
        let frame = original
            .step_harmonic_into(&prepared.mechanics, &[cosine], &[sine], &mut workspace)
            .unwrap();
        assert_eq!(observed.model.states(), original.states());
        assert_eq!(observed.last_modal_frame(), frame);
        (cosine, sine) = prepared.advance_force_phase(cosine, sine).unwrap();
    }
}

#[test]
fn simultaneous_tones_preserve_interference_dc_and_duplicate_cancellation() {
    let physical_modes = vec![mode(1200.0, 0.02), mode(1700.0, 0.04)];
    let mut observed = observer(
        physical_modes.clone(),
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    let prepared = observed
        .prepare_tones(
            &[0.0, 600.0, 1300.0, 600.0, 600.0],
            ModalToneBudget::realtime_reference(),
        )
        .unwrap();
    let mut cosine = vec![0.4, 1e16, 0.25, 1.0, -1e16, -0.2, -3e15, -0.17, 0.6, 3e15];
    let mut sine = vec![0.0, 2e16, -0.3, -0.3, -2e16, 0.0, -1e16, 0.11, -0.15, 1e16];
    let forces = [
        [
            C64::new(0.4, 0.0),
            C64::new(1.0, -0.3),
            C64::new(0.25, -0.3),
        ],
        [
            C64::new(-0.2, 0.0),
            C64::new(0.6, -0.15),
            C64::new(-0.17, 0.11),
        ],
    ];
    let mut original = observed.model.clone();
    let mut workspace = ModalAcousticWorkspace::new(&original);
    for sample in 1..=384 {
        let pressure = observed.step_tones(&prepared, &cosine, &sine).unwrap();
        let t = sample as f64 / f64::from(RATE);
        let expected = physical_modes
            .iter()
            .zip(forces)
            .map(|(&mode, row)| {
                [0.0, 600.0, 1300.0]
                    .iter()
                    .zip(row)
                    .map(|(&w, force)| analytical_pressure(t, mode, w, force, 0.0, 0.0))
                    .sum::<f64>()
            })
            .sum::<f64>();
        close(pressure, expected);
        let frame = original
            .step_tones_into(&prepared.mechanics, &cosine, &sine, &mut workspace)
            .unwrap();
        assert_eq!(observed.model.states(), original.states());
        assert_eq!(observed.last_modal_frame(), frame);
        for index in 0..cosine.len() {
            (cosine[index], sine[index]) = prepared
                .advance_force_phase(index % prepared.tone_count(), cosine[index], sine[index])
                .unwrap();
        }
    }
}

#[test]
fn rejected_harmonic_steps_leave_warm_mechanics_frame_filter_and_delay_unchanged() {
    for mechanical_refusal in [false, true] {
        let budget = ModalAcousticTimeBudget {
            maximum_abs_pressure_pa: if mechanical_refusal { 1e6 } else { 1e-4 },
            maximum_abs_velocity_m_sqrt_kg_per_s: if mechanical_refusal { 0.001 } else { 1e6 },
            ..ModalAcousticTimeBudget::audible_reference()
        };
        let modes = vec![mode(1200.0, 0.02), mode(1700.0, 0.04)];
        let mut observed = observer(
            modes.clone(),
            budget,
            if mechanical_refusal { 3 } else { 0 },
        );
        let mut reference = observer(modes, budget, if mechanical_refusal { 3 } else { 0 });
        let prepared = observed.prepare_harmonic(740.3).unwrap();
        let reference_plan = reference.prepare_harmonic(740.3).unwrap();
        for _ in 0..8 {
            observed
                .step_harmonic(&prepared, &[0.001, 0.001], &[0.0, 0.0])
                .unwrap();
            reference
                .step_harmonic(&reference_plan, &[0.001, 0.001], &[0.0, 0.0])
                .unwrap();
        }
        let before = (
            observed.model.states().to_vec(),
            observed.last_modal_frame().clone(),
            observed
                .filters
                .iter()
                .map(|f| f.state.clone())
                .collect::<Vec<_>>(),
            observed.delay.history.clone(),
            observed.delay.cursor,
        );
        assert!(
            observed
                .step_harmonic(&prepared, &[0.01, 1e6], &[0.0, 0.0])
                .is_err()
        );
        assert!(
            observed
                .step_harmonic(&prepared, &[f64::NAN, 0.0], &[0.0, 0.0])
                .is_err()
        );
        assert!(
            observed
                .step_harmonic(&reference_plan, &[0.0, 0.0], &[0.0, 0.0])
                .is_err()
        );
        let after = (
            observed.model.states().to_vec(),
            observed.last_modal_frame().clone(),
            observed
                .filters
                .iter()
                .map(|f| f.state.clone())
                .collect::<Vec<_>>(),
            observed.delay.history.clone(),
            observed.delay.cursor,
        );
        assert_eq!(before, after);
        assert_eq!(
            observed
                .step_harmonic(&prepared, &[0.001, 0.001], &[0.0, 0.0])
                .unwrap(),
            reference
                .step_harmonic(&reference_plan, &[0.001, 0.001], &[0.0, 0.0])
                .unwrap()
        );
        // A valid ramp following the refusal also retains the original warm history.
        observed.prepare_linear_forcing().unwrap();
        reference.prepare_linear_forcing().unwrap();
        assert_eq!(
            observed
                .step_linear(&[0.001, 0.001], &[0.002, 0.002])
                .unwrap(),
            reference
                .step_linear(&[0.001, 0.001], &[0.002, 0.002])
                .unwrap()
        );
    }
}

#[test]
fn dc_matches_held_observation_and_only_real_constant_feedthrough_is_admitted() {
    let physical_mode = mode(1200.0, 0.02);
    let mut continuous = observer(
        vec![physical_mode],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    let mut held = observer(
        vec![physical_mode],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    let dc = continuous.prepare_harmonic(0.0).unwrap();
    for _ in 0..64 {
        close(
            continuous.step_harmonic(&dc, &[0.4], &[-11.0]).unwrap(),
            held.step(&[0.4]).unwrap(),
        );
    }
    let mut gain = observer(
        vec![mode(0.0, 0.0)],
        ModalAcousticTimeBudget::audible_reference(),
        0,
    );
    gain.constant = Some(vec![C64::from_re(2.1)]);
    gain.filters.clear();
    let plan = gain.prepare_harmonic(600.0).unwrap();
    let (mut cosine, mut sine) = (0.4, -0.3);
    for _ in 0..32 {
        let pressure = gain.step_harmonic(&plan, &[cosine], &[sine]).unwrap();
        assert_eq!(
            pressure,
            2.1 * gain.model.states()[0].velocity_m_sqrt_kg_per_s
        );
        (cosine, sine) = plan.advance_force_phase(cosine, sine).unwrap();
    }
    gain.constant = Some(vec![C64::new(2.1, -0.37)]);
    assert!(gain.prepare_harmonic(600.0).is_err());
    assert!(
        gain.prepare_tones(&[600.0, 1300.0], ModalToneBudget::realtime_reference())
            .is_err()
    );
}
