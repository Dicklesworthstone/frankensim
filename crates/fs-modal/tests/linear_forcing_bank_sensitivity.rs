//! Physical-rate bank integration, immutable cache lifetime and transient fitting.

use fs_math::c64::C64;
use fs_modal::time::{
    LinearForceTransition, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn modes(frequencies: &[f64]) -> Vec<ModalAcousticMode> {
    frequencies
        .iter()
        .map(|&angular_frequency_rad_s| ModalAcousticMode {
            angular_frequency_rad_s,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        })
        .collect()
}

fn budget() -> ModalAcousticTimeBudget {
    ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 10.0,
        maximum_abs_velocity_m_sqrt_kg_per_s: 100.0,
        maximum_total_energy_j: 100.0,
        ..ModalAcousticTimeBudget::audible_reference()
    }
}

fn bank(frequencies: &[f64], rates: &[f64]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_damping_rates(20, modes(frequencies), rates, budget())
        .unwrap()
}

#[test]
fn g1_mixed_bank_derivatives_use_the_original_physical_rate_transitions() {
    let frequencies = [0.0, 1e-310, 2.0, 2.0, 2.0];
    let rates = [2.5, 3.0, 0.0, 4.0, 20.0];
    let initial = state(0.3, -0.2);
    let mut model = bank(&frequencies, &rates);
    model.restore_states(&[initial; 5]).unwrap();
    assert!(model.modes()[0].damping_ratio.is_infinite());
    assert!(model.modes()[1].damping_ratio.is_infinite());
    model.prepare_linear_forcing_sensitivities().unwrap();
    let before = model.states().to_vec();
    let mut expected_states = Vec::new();
    for i in 0..frequencies.len() {
        let scalar = LinearForceTransition::try_new_with_damping_rate(
            frequencies[i],
            rates[i],
            model.sample_period_s(),
        )
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
        let frame = model
            .advance_exact_linear_step_with_sensitivities(i, initial, 0.7, -0.4)
            .unwrap();
        assert_eq!(frame, scalar.advance(initial, 0.7, -0.4).unwrap());
        assert_eq!(
            frame.state,
            model
                .advance_exact_linear_step(i, initial, 0.7, -0.4)
                .unwrap()
        );
        assert_eq!(
            model
                .pullback_exact_linear_step(i, initial, 0.7, -0.4, 0.6, -0.3)
                .unwrap(),
            scalar.pullback(initial, 0.7, -0.4, 0.6, -0.3).unwrap()
        );
        expected_states.push(frame.state);
    }
    assert_eq!(model.states(), before);
    let mut workspace = ModalAcousticWorkspace::new(&model);
    model
        .step_linear_into(&[0.7; 5], &[-0.4; 5], &mut workspace)
        .unwrap();
    assert_eq!(model.states(), expected_states);
}

#[test]
fn g3_legacy_ratio_bank_preserves_the_original_prepared_arithmetic() {
    let mut original = modes(&[0.3, 2.7]);
    for (mode, ratio) in original.iter_mut().zip([0.17, 0.023]) {
        mode.damping_ratio = ratio;
    }
    let model = ModalAcousticTimeModel::try_new(20, original.clone(), budget()).unwrap();
    model.prepare_linear_forcing_sensitivities().unwrap();
    let initial = state(0.3, -0.2);
    for (i, mode) in original.iter().enumerate() {
        let expected = LinearForceTransition::try_new(
            mode.angular_frequency_rad_s,
            mode.damping_ratio,
            model.sample_period_s(),
        )
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
        assert_eq!(
            model
                .advance_exact_linear_step_with_sensitivities(i, initial, 0.7, -0.4)
                .unwrap(),
            expected.advance(initial, 0.7, -0.4).unwrap()
        );
    }
    // An unchanged legacy damping law retains the warmed sensitivity cache.
    let same = model.with_damping_ratios(&[0.17, 0.023]).unwrap();
    assert_eq!(
        model
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .unwrap(),
        same.pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .unwrap()
    );
}

#[test]
fn g0_preparation_is_shared_and_changed_mechanics_cannot_reuse_old_gradients() {
    let initial = state(-0.2, 0.3);
    let mut source = bank(&[0.0], &[2.5]);
    source.restore_states(&[initial]).unwrap();
    let clone = source.clone();
    assert!(
        source
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .is_err()
    );
    assert!(
        source
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .is_err()
    );
    clone.prepare_linear_forcing_sensitivities().unwrap();
    let before = source
        .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
        .unwrap();
    assert_eq!(
        before,
        clone
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .unwrap()
    );
    let same = source.with_damping_rates(&[2.5]).unwrap();
    assert_eq!(
        before,
        same.advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .unwrap()
    );

    let changed = source.with_damping_rates(&[5.0]).unwrap();
    // The diagnostic ratios are both infinity, while the physics has changed.
    assert_eq!(source.modes(), changed.modes());
    assert_eq!(changed.states(), source.states());
    assert!(
        changed
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .is_err()
    );
    assert!(
        changed
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .is_err()
    );
    changed.prepare_linear_forcing_sensitivities().unwrap();
    let after = changed
        .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
        .unwrap();
    assert_ne!(before.parameters, after.parameters);
    assert_eq!(
        before,
        source
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .unwrap()
    );
    assert_eq!(
        before,
        source
            .without_pressure_transfers()
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .unwrap()
    );

    let reset = source.with_damping_ratios(&[0.0]).unwrap();
    assert!(
        reset
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .is_err()
    );
    reset.prepare_linear_forcing_sensitivities().unwrap();
    let free = LinearForceTransition::try_new_with_damping_rate(0.0, 0.0, source.sample_period_s())
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
    assert_eq!(
        reset
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, 0.6, -0.3)
            .unwrap(),
        free.pullback(initial, 0.7, -0.4, 0.6, -0.3).unwrap()
    );
    assert!(
        source
            .advance_exact_linear_step_with_sensitivities(1, initial, 0.7, -0.4)
            .is_err()
    );
    assert!(
        source
            .pullback_exact_linear_step(1, initial, 0.7, -0.4, 0.6, -0.3)
            .is_err()
    );
    assert!(
        source
            .pullback_exact_linear_step(0, initial, 0.7, -0.4, f64::NAN, 0.0)
            .is_err()
    );
    source.prepare_linear_forcing_sensitivities().unwrap();
    assert_eq!(source.states(), &[initial]);
    assert_eq!(source.damping_rates_per_s(), &[2.5]);
    assert_eq!(
        before,
        source
            .advance_exact_linear_step_with_sensitivities(0, initial, 0.7, -0.4)
            .unwrap()
    );
}

const MASSES: [f64; 3] = [2.0, 0.5, 1.5];
const DASHPOT_SHAPES: [f64; 3] = [1.0, 0.4, 1.5];
const STEPS: usize = 24;

fn physical_bank(coefficient: f64) -> ModalAcousticTimeModel {
    // A shared dimensional coefficient C [kg/s] multiplies known dashpot
    // shapes. c_i = C*shape_i/m_i, with fixed physical masses and mode basis.
    let rates = std::array::from_fn::<_, 3, _>(|i| coefficient * DASHPOT_SHAPES[i] / MASSES[i]);
    let mut model = bank(&[0.0, 1e-310, 2.0], &rates);
    let q = [0.2, -0.1, 0.05];
    let v = [0.4, 0.25, 0.3];
    let initial =
        std::array::from_fn::<_, 3, _>(|i| state(q[i] * MASSES[i].sqrt(), v[i] * MASSES[i].sqrt()));
    model.restore_states(&initial).unwrap();
    model.prepare_linear_forcing_sensitivities().unwrap();
    model
}

fn forces(sample: usize) -> [f64; 3] {
    let force = 0.2 + 0.15 * (sample as f64 / 20.0);
    std::array::from_fn(|i| force * [1.0, -0.5, 0.3][i] / MASSES[i].sqrt())
}

fn record(model: &mut ModalAcousticTimeModel) -> Vec<[ModalAcousticState; 3]> {
    let mut history = vec![std::array::from_fn(|i| model.states()[i])];
    let mut workspace = ModalAcousticWorkspace::new(model);
    for sample in 0..STEPS {
        model
            .step_linear_into(&forces(sample), &forces(sample + 1), &mut workspace)
            .unwrap();
        history.push(std::array::from_fn(|i| model.states()[i]));
    }
    history
}

fn calibration(coefficient: f64, observed: &[[ModalAcousticState; 3]]) -> (f64, f64) {
    let mut model = physical_bank(coefficient);
    let history = record(&mut model);
    let mut loss = 0.0;
    for sample in 1..=STEPS {
        for (i, mass) in MASSES.iter().enumerate() {
            let q = history[sample][i].displacement_m_sqrt_kg
                - observed[sample][i].displacement_m_sqrt_kg;
            let v = history[sample][i].velocity_m_sqrt_kg_per_s
                - observed[sample][i].velocity_m_sqrt_kg_per_s;
            // Physical displacement/velocity observation weights; divide out
            // the fixed mass normalization before forming the measured loss.
            loss += 0.5 * (q * q + 0.25 * v * v) / mass;
        }
    }
    let mut qbar = [0.0; 3];
    let mut vbar = [0.0; 3];
    let mut gradient = 0.0;
    for sample in (0..STEPS).rev() {
        let start = forces(sample);
        let end = forces(sample + 1);
        for i in 0..3 {
            qbar[i] += (history[sample + 1][i].displacement_m_sqrt_kg
                - observed[sample + 1][i].displacement_m_sqrt_kg)
                / MASSES[i];
            vbar[i] += 0.25
                * (history[sample + 1][i].velocity_m_sqrt_kg_per_s
                    - observed[sample + 1][i].velocity_m_sqrt_kg_per_s)
                / MASSES[i];
            let pulled = model
                .pullback_exact_linear_step(
                    i,
                    history[sample][i],
                    start[i],
                    end[i],
                    qbar[i],
                    vbar[i],
                )
                .unwrap();
            gradient += pulled.damping_rate * DASHPOT_SHAPES[i] / MASSES[i];
            qbar[i] = pulled.inputs.initial_displacement;
            vbar[i] = pulled.inputs.initial_velocity;
        }
    }
    assert_eq!(model.states(), history.last().unwrap());
    (loss, gradient)
}

#[test]
fn g3_transient_bank_recovers_a_shared_dimensional_dashpot_through_critical_damping() {
    // At C=4, the elastic coordinate has c=4=2*omega. The other coordinates
    // retain positive drag despite undefined/overflowing diagnostic ratios.
    let truth = 4.0;
    let observed = record(&mut physical_bank(truth));
    let (value, gradient) = calibration(2.7, &observed);
    let h = 1e-4;
    let difference = (calibration(2.7 - 2.0 * h, &observed).0
        - 8.0 * calibration(2.7 - h, &observed).0
        + 8.0 * calibration(2.7 + h, &observed).0
        - calibration(2.7 + 2.0 * h, &observed).0)
        / (12.0 * h);
    assert!((gradient - difference).abs() < 2e-8 * difference.abs());
    assert!(value > 0.0);

    // A test-only scalar bracket locates the gradient root. Production users
    // can pass the same loss/gradient to the existing optimization engines.
    let (mut lower, mut upper) = (1.3, 6.8);
    assert!(calibration(lower, &observed).1 < 0.0);
    assert!(calibration(upper, &observed).1 > 0.0);
    for _ in 0..45 {
        let middle = 0.5 * (lower + upper);
        if calibration(middle, &observed).1 > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    let recovered = 0.5 * (lower + upper);
    assert!((recovered - truth).abs() < 2e-8, "recovered C={recovered}");
    assert!(calibration(recovered, &observed).0 < value * 1e-18);
}
