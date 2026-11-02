//! G1/G3 exact first-order-hold forcing for the shared scalar oscillator.

use fs_math::c64::C64;
use fs_modal::time::{
    LinearForceTransition, ModalAcousticMode, ModalAcousticState, advance_exact_zoh,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 8.0e-12 * scale.max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn g1_ramps_match_independent_high_precision_matrix_exponentials() {
    // Independent 100-decimal-digit Decimal matrix exponential of
    // h * [[0,1,0,0],[-omega^2,-2*zeta*omega,1,0],
    //      [0,0,0,1],[0,0,0,0]] applied to [q,v,F0,(F1-F0)/h].
    // Every input was lifted from its exact f64 value. The full 4x4 Taylor
    // exponential was scaled below entry norm 1/8, summed below 1e-98, and
    // squared back. No scalar gain or modal transition implementation enters
    // these constants. Cases include both sides of critical damping and a
    // strongly overdamped slow pole, as well as soft and high-Q modes.
    let cases = [
        (
            1.0e-6,
            0.05,
            1.0 / 48_000.0,
            100.0,
            300.0,
            0.0,
            0.0,
            3.616_898_148_145_887_6e-8,
            0.004_166_666_666_663_05,
        ),
        (
            1000.0,
            1.0e-12,
            1.0 / 48_000.0,
            24_200.0,
            -12_000.0,
            0.0002,
            -0.3,
            0.000_196_340_019_626_557_98,
            -0.177_040_761_195_166_73,
        ),
        (
            1000.0,
            0.999_999_995,
            0.002,
            3.0,
            11.0,
            0.0002,
            -0.3,
            3.947_346_362_220_115e-6,
            -0.010_345_539_922_932_839,
        ),
        (
            1000.0,
            1.0,
            0.002,
            3.0,
            11.0,
            0.0002,
            -0.3,
            3.947_346_982_656_295e-6,
            -0.010_345_540_023_080_949,
        ),
        (
            1000.0,
            1.000_000_005,
            0.002,
            3.0,
            11.0,
            0.0002,
            -0.3,
            3.947_347_603_092_472e-6,
            -0.010_345_540_123_229_057,
        ),
        (
            100.0,
            1.0e12,
            0.1,
            2.0e14,
            -1.0e14,
            0.0,
            0.0,
            0.024_999_999_999_877_502,
            -0.500_000_000_001_175,
        ),
        (
            17.0,
            0.15,
            0.3,
            9.0,
            -4.0,
            0.12,
            -0.6,
            0.007_831_563_946_827_692,
            0.404_752_711_637_670_1,
        ),
        (
            3.0,
            0.0,
            0.8,
            -0.4,
            1.2,
            0.3,
            0.2,
            -0.125_661_266_320_141_27,
            -0.459_369_870_668_605_96,
        ),
        (
            1.0e-160,
            0.0,
            1.0e-4,
            1.0,
            3.0,
            0.0,
            0.0,
            8.333_333_333_333_334e-9,
            0.0002,
        ),
    ];
    for (omega, zeta, dt, f0, f1, q, v, expected_q, expected_v) in cases {
        let step = LinearForceTransition::try_new(omega, zeta, dt).unwrap();
        let actual = step.advance(state(q, v), f0, f1).unwrap();
        close(actual.displacement_m_sqrt_kg, expected_q, expected_q.abs());
        close(
            actual.velocity_m_sqrt_kg_per_s,
            expected_v,
            expected_v.abs(),
        );
    }
}

#[test]
fn g1_free_coordinate_integrates_the_cubic_ramp_exactly() {
    let dt = 0.125;
    let initial = state(0.2, -0.4);
    let (f0, f1) = (3.0, -5.0);
    let step = LinearForceTransition::try_new(0.0, 0.0, dt).unwrap();
    let actual = step.advance(initial, f0, f1).unwrap();
    let q = initial.displacement_m_sqrt_kg
        + dt * initial.velocity_m_sqrt_kg_per_s
        + dt * dt * (f0 / 3.0 + f1 / 6.0);
    let v = initial.velocity_m_sqrt_kg_per_s + dt * (f0 + f1) / 2.0;
    close(actual.displacement_m_sqrt_kg, q, q.abs());
    close(actual.velocity_m_sqrt_kg_per_s, v, v.abs());
}

#[test]
fn g1_affine_motion_matches_viscous_and_elastic_force_balance() {
    // Prescribe q(t)=q0+v0*t, hence q''=0 and the exact required force is
    // F(t)=2*zeta*omega*v0+omega^2*q(t), itself an affine function.
    let (omega, dt) = (21.0, 0.3);
    let initial = state(0.03, -0.02);
    let expected_q = initial.displacement_m_sqrt_kg + dt * initial.velocity_m_sqrt_kg_per_s;
    for zeta in [0.0, 0.4, 1.0, 5.0] {
        let drag = 2.0 * zeta * omega * initial.velocity_m_sqrt_kg_per_s;
        let f0 = drag + omega * omega * initial.displacement_m_sqrt_kg;
        let f1 = drag + omega * omega * expected_q;
        let actual = LinearForceTransition::try_new(omega, zeta, dt)
            .unwrap()
            .advance(initial, f0, f1)
            .unwrap();
        close(actual.displacement_m_sqrt_kg, expected_q, 0.03);
        close(
            actual.velocity_m_sqrt_kg_per_s,
            initial.velocity_m_sqrt_kg_per_s,
            0.02,
        );
    }
}

#[test]
fn g3_linear_ramp_is_invariant_under_subdivision_and_replay() {
    let initial = state(0.23, -0.17);
    let (dt, f0, f1) = (0.4, 1.3, -0.8);
    for zeta in [0.0, 0.07, 0.999_999_999_9, 1.0, 1.000_000_000_1, 3.0, 1.0e5] {
        let whole = LinearForceTransition::try_new(9.2, zeta, dt).unwrap();
        let expected = whole.advance(initial, f0, f1).unwrap();
        let part = LinearForceTransition::try_new(9.2, zeta, dt / 17.0).unwrap();
        let mut actual = initial;
        for i in 0..17 {
            let left = f0 + (f1 - f0) * f64::from(i) / 17.0;
            let right = f0 + (f1 - f0) * f64::from(i + 1) / 17.0;
            actual = part.advance(actual, left, right).unwrap();
        }
        close(
            actual.displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            0.23,
        );
        close(
            actual.velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            1.0,
        );
        assert_eq!(whole.advance(initial, f0, f1).unwrap(), expected);
        assert_eq!(initial, state(0.23, -0.17));
    }
}

#[test]
fn g1_undamped_sampled_sinusoid_retains_resonant_growth() {
    // For piecewise-linear samples of F*sin(omega*t), the forcing's
    // resonant Fourier coefficient is F*sinc(pi/samples_per_cycle)^2.
    // Integrating the exact convolution over complete periods gives
    // q(N*T) = -F*N*pi/omega^2*sinc(...)^2 and v(N*T) = 0.
    let (omega, amplitude) = (7.0, 3.0);
    let samples_per_cycle = 16;
    let cycles = 20;
    let dt = core::f64::consts::TAU / (omega * f64::from(samples_per_cycle));
    let step = LinearForceTransition::try_new(omega, 0.0, dt).unwrap();
    let mut actual = state(0.0, 0.0);
    for i in 0..cycles * samples_per_cycle {
        let sample = |j| {
            amplitude
                * fs_math::det::sin(
                    core::f64::consts::TAU * f64::from(j % samples_per_cycle)
                        / f64::from(samples_per_cycle),
                )
        };
        actual = step.advance(actual, sample(i), sample(i + 1)).unwrap();
    }
    let half_angle = core::f64::consts::PI / f64::from(samples_per_cycle);
    let sinc = half_angle.sin() / half_angle;
    let expected_q =
        -amplitude * f64::from(cycles) * core::f64::consts::PI / (omega * omega) * sinc * sinc;
    close(actual.displacement_m_sqrt_kg, expected_q, expected_q.abs());
    close(
        actual.velocity_m_sqrt_kg_per_s,
        0.0,
        omega * expected_q.abs(),
    );
}

#[test]
fn g3_constant_endpoint_forces_agree_with_the_shared_held_force_path() {
    let initial = state(0.17, -0.29);
    for (omega, zeta) in [
        (0.0, 0.0),
        (11.0, 0.0),
        (11.0, 0.2),
        (11.0, 1.0),
        (11.0, 4.0),
    ] {
        let mode = ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        };
        for dt in [1.0e-8, 0.03, 2.0] {
            let actual = LinearForceTransition::try_new(omega, zeta, dt)
                .unwrap()
                .advance(initial, 1.4, 1.4)
                .unwrap();
            let expected = advance_exact_zoh(mode, initial, 1.4, dt);
            close(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                1.0,
            );
            close(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
                3.0,
            );
        }
    }
}

#[test]
fn g0_invalid_parameters_and_inputs_refuse_without_mutating_state() {
    for (omega, zeta, dt) in [
        (-1.0, 0.1, 0.1),
        (1.0, -0.1, 0.1),
        (1.0, 0.1, 0.0),
        (1.0, 0.1, -1.0),
        (0.0, 0.1, 0.1),
        (f64::INFINITY, 0.1, 0.1),
        (1.0, f64::NAN, 0.1),
        (1.0, 0.1, f64::INFINITY),
    ] {
        assert!(LinearForceTransition::try_new(omega, zeta, dt).is_err());
    }
    let step = LinearForceTransition::try_new(7.0, 0.02, 0.1).unwrap();
    let initial = state(0.2, -0.3);
    assert!(step.advance(initial, f64::NAN, 1.0).is_err());
    assert!(step.advance(initial, 1.0, f64::INFINITY).is_err());
    assert!(step.advance(state(f64::NAN, 0.0), 0.0, 0.0).is_err());
    assert!(step.advance(state(0.0, f64::INFINITY), 0.0, 0.0).is_err());
    assert_eq!(initial, state(0.2, -0.3));
}

#[test]
fn g3_pullback_is_the_transpose_in_every_damping_regime() {
    for (omega, zeta) in [(0.0, 0.0), (7.0, 0.0), (7.0, 0.2), (7.0, 1.0), (7.0, 4.0)] {
        for h in [1.0e-10, 0.03, 2.0] {
            let step = LinearForceTransition::try_new(omega, zeta, h).unwrap();
            let direction = [0.4, -0.7, 1.2, -0.8];
            let tangent = step
                .advance(
                    state(direction[0], direction[1]),
                    direction[2],
                    direction[3],
                )
                .unwrap();
            let gradient = step.pullback(-0.3, 0.9).unwrap();
            let forward =
                -0.3 * tangent.displacement_m_sqrt_kg + 0.9 * tangent.velocity_m_sqrt_kg_per_s;
            let reverse = gradient.initial_displacement * direction[0]
                + gradient.initial_velocity * direction[1]
                + gradient.force_start * direction[2]
                + gradient.force_end * direction[3];
            close(reverse, forward, 1.0);
            let seeded = step.pullback(0.0, 1.0).unwrap();
            let q_direction = step.advance(state(1.0, 0.0), 0.0, 0.0).unwrap();
            close(
                seeded.initial_displacement,
                q_direction.velocity_m_sqrt_kg_per_s,
                q_direction.velocity_m_sqrt_kg_per_s.abs(),
            );
        }
    }
}

#[test]
fn g3_trajectory_pullback_accumulates_shared_force_samples() {
    let step = LinearForceTransition::try_new(4.0, 0.3, 0.12).unwrap();
    let source = [0.2, -0.1, 1.0, -0.3, 0.8, 0.4]; // q0, v0, four force samples
    let evaluate = |x: &[f64; 6]| {
        let mut value = state(x[0], x[1]);
        for forces in x[2..].windows(2) {
            value = step.advance(value, forces[0], forces[1]).unwrap();
        }
        0.7 * value.displacement_m_sqrt_kg * value.displacement_m_sqrt_kg
            + 0.4 * value.velocity_m_sqrt_kg_per_s * value.velocity_m_sqrt_kg_per_s
    };
    let mut endpoint = state(source[0], source[1]);
    for forces in source[2..].windows(2) {
        endpoint = step.advance(endpoint, forces[0], forces[1]).unwrap();
    }
    let (mut q_bar, mut v_bar) = (
        1.4 * endpoint.displacement_m_sqrt_kg,
        0.8 * endpoint.velocity_m_sqrt_kg_per_s,
    );
    let mut derivatives = [0.0; 6];
    for interval in (0..3).rev() {
        let gradient = step.pullback(q_bar, v_bar).unwrap();
        derivatives[interval + 2] += gradient.force_start;
        derivatives[interval + 3] += gradient.force_end;
        (q_bar, v_bar) = (gradient.initial_displacement, gradient.initial_velocity);
    }
    derivatives[0] = q_bar;
    derivatives[1] = v_bar;
    for index in 0..6 {
        let (mut plus, mut minus) = (source, source);
        plus[index] += 1.0e-5;
        minus[index] -= 1.0e-5;
        let difference = (evaluate(&plus) - evaluate(&minus)) / 2.0e-5;
        close(derivatives[index], difference, 10.0);
    }
}

#[test]
fn g0_pullback_refuses_invalid_or_unrepresentable_cotangents() {
    let step = LinearForceTransition::try_new(7.0, 0.2, 0.1).unwrap();
    for seeds in [
        (f64::NAN, 0.0),
        (0.0, f64::INFINITY),
        (f64::NEG_INFINITY, 1.0),
    ] {
        assert!(step.pullback(seeds.0, seeds.1).is_err());
    }
    let long = LinearForceTransition::try_new(0.0, 0.0, 1.0e200).unwrap();
    assert!(long.pullback(1.0, 0.0).is_err());
    let zero = step.pullback(0.0, 0.0).unwrap();
    assert_eq!(
        [
            zero.initial_displacement,
            zero.initial_velocity,
            zero.force_start,
            zero.force_end
        ],
        [0.0; 4]
    );
    assert_eq!(
        step.pullback(0.4, -0.2).unwrap(),
        step.pullback(0.4, -0.2).unwrap()
    );
}
