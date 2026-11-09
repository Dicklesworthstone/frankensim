//! Independent variational references and full transient calibration gradients.

use fs_modal::time::{LinearForceTransition, ModalAcousticState};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= relative * expected.abs().max(1e-100),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

// G1: independent 85-digit Decimal matrix exponentials of the augmented
// variational ODE. In order [q,v,F,Fdot,q_k,v_k,q_c,v_c], its nonzero equations
// are q'=v; v'=-k*q-c*v+F; F'=Fdot; q_k'=v_k;
// v_k'=-q-k*q_k-c*v_k; q_c'=v_c; v_c'=-v-k*q_c-c*v_c.
// exp(h*A) used 80 Taylor terms at induced infinity norm <= 1/8, then ordinary
// dense matrix squaring. This oracle uses neither the production scalar moment
// recurrence nor finite differences. Constants are rounded once to binary64.
#[test]
fn g1_physical_coefficient_sensitivities_match_independent_variational_flow() {
    let cases = [
        (
            0.0,
            0.0,
            0.5,
            [
                2.41666666666666669e-1,
                -1.25e-1,
                -3.45833333333333340e-2,
                -1.33854166666666663e-1,
                1.61458333333333350e-2,
                5.83333333333333343e-2,
            ],
        ),
        (
            0.0,
            3.0,
            0.7,
            [
                2.79391964032218054e-1,
                -3.31758920966540791e-2,
                -3.84096388221598431e-2,
                -8.21959844116700311e-2,
                5.55837276303178425e-3,
                3.93291767868662011e-3,
            ],
        ),
        (
            2.0,
            0.0,
            0.7,
            [
                1.62714025399887767e-2,
                -6.06440197794925751e-1,
                -4.65232641839289640e-2,
                -7.61960621456850978e-2,
                7.16213973525839237e-2,
                2.20811797065799342e-1,
            ],
        ),
        (
            2.0,
            0.2,
            0.7,
            [
                3.00390274637473663e-2,
                -5.64645772707253579e-1,
                -4.47670910239980069e-2,
                -7.18118463194998280e-2,
                6.61429534517987139e-2,
                1.97621733036958042e-1,
            ],
        ),
        (
            2.0,
            4.0,
            0.7,
            [
                1.67621113312512698e-1,
                -2.26932614262147642e-1,
                -2.52203508041503009e-2,
                -3.30158114610107634e-2,
                1.87695509667265940e-2,
                3.57319000751387811e-2,
            ],
        ),
        (
            2.0,
            3.999999,
            0.7,
            [
                1.67621094542959204e-1,
                -2.26932649994053853e-1,
                -2.52203538202145153e-2,
                -3.30158160216332819e-2,
                1.87695560727067501e-2,
                3.57319123380307468e-2,
            ],
        ),
        (
            2.0,
            4.000001,
            0.7,
            [
                1.67621132082061114e-1,
                -2.26932578530253698e-1,
                -2.52203477880867770e-2,
                -3.30158069003897020e-2,
                1.87695458607483566e-2,
                3.57318878122529424e-2,
            ],
        ),
        (
            2.0,
            20.0,
            0.7,
            [
                2.59651462945895084e-1,
                -6.86553589808852754e-2,
                -8.59545089307815584e-3,
                -1.15450509919591722e-2,
                1.74520819295517879e-3,
                2.94647018193565586e-3,
            ],
        ),
        (
            10.0,
            1000.0,
            1.0,
            [
                2.71429107879122500e-1,
                -2.75445654099319946e-2,
                -2.71230561096823047e-4,
                -2.44358038992944799e-4,
                2.71447773705440073e-5,
                2.48342262969426916e-5,
            ],
        ),
        (
            1000.0,
            1e-8,
            0.005,
            [
                8.52896309634596228e-2,
                2.87619721113345975e2,
                7.18954875366226783e-7,
                -6.93859429985537016e-5,
                -3.57063225390303284e-4,
                -7.19145850690713795e-1,
            ],
        ),
    ];
    let initial = state(0.3, -0.2);
    for (omega, c, h, expected) in cases {
        let motion = LinearForceTransition::try_new_with_damping_rate(omega, c, h).unwrap();
        let prepared = motion.prepare_sensitivities().unwrap();
        let frame = prepared.advance(initial, 0.7, -0.4).unwrap();
        assert_eq!(frame.state, motion.advance(initial, 0.7, -0.4).unwrap());
        let derivative = frame.parameters;
        let actual = [
            frame.state.displacement_m_sqrt_kg,
            frame.state.velocity_m_sqrt_kg_per_s,
            derivative.displacement_by_stiffness,
            derivative.velocity_by_stiffness,
            derivative.displacement_by_damping_rate,
            derivative.velocity_by_damping_rate,
        ];
        for (value, reference) in actual.into_iter().zip(expected) {
            close(value, reference, 2e-11);
        }
    }
}

#[test]
fn g1_critical_sensitivities_retain_the_decaying_tail() {
    let motion = LinearForceTransition::try_new_with_damping_rate(1.0, 2.0, 40.0).unwrap();
    let frame = motion
        .prepare_sensitivities()
        .unwrap()
        .advance(state(1.0, 0.0), 0.0, 0.0)
        .unwrap();
    let d = frame.parameters;
    // Same independent augmented exponential as above, with q0=1, v0=F=0.
    for (actual, expected) in [
        (frame.state.displacement_m_sqrt_kg, 1.74182524466955141e-16),
        (
            frame.state.velocity_m_sqrt_kg_per_s,
            -1.69934170211663555e-16,
        ),
        (d.displacement_by_stiffness, -4.87144621273435562e-14),
        (d.velocity_by_stiffness, 4.51458445528986165e-14),
        (d.displacement_by_damping_rate, 4.53157787231102821e-14),
        (d.velocity_by_damping_rate, -4.19170953188770142e-14),
    ] {
        close(actual, expected, 2e-11);
    }
}

#[test]
fn g1_held_force_settling_keeps_damping_tails_and_exact_equilibrium() {
    for h in [20.0_f64, 40.0, 80.0] {
        let flow = LinearForceTransition::try_new_with_damping_rate(1.0, 2.0, h)
            .unwrap()
            .prepare_sensitivities()
            .unwrap();
        let driven = flow.advance(state(0.0, 0.0), 1.0, 1.0).unwrap().parameters;
        // Independent critical-limit derivatives at fixed k=1, c=2:
        // dP/dc = h^3*exp(-h)/6; dB/dc = (h^3/6-h^2/2)*exp(-h).
        // G=(1-P)/k, so the held-force displacement derivative is -dP/dc.
        let tail = (-h).exp();
        close(
            driven.displacement_by_damping_rate,
            -h * h * h * tail / 6.0,
            2e-11,
        );
        close(
            driven.velocity_by_damping_rate,
            h * h * (h / 6.0 - 0.5) * tail,
            2e-11,
        );
        let balanced = flow.advance(state(1.0, 0.0), 1.0, 1.0).unwrap().parameters;
        assert_eq!(balanced.displacement_by_damping_rate, 0.0);
        assert_eq!(balanced.velocity_by_damping_rate, 0.0);
    }
}

#[test]
fn g3_subdividing_the_same_ramp_preserves_coefficient_sensitivities() {
    let initial = state(0.3, -0.2);
    for c in [0.0, 0.2, 4.0, 20.0] {
        let full = LinearForceTransition::try_new_with_damping_rate(2.0, c, 0.7)
            .unwrap()
            .prepare_sensitivities()
            .unwrap();
        let half = LinearForceTransition::try_new_with_damping_rate(2.0, c, 0.35)
            .unwrap()
            .prepare_sensitivities()
            .unwrap();
        let first = half.advance(initial, 0.7, 0.15).unwrap();
        for (qbar, vbar) in [(1.0, 0.0), (0.0, 1.0), (0.7, -0.4)] {
            let whole = full.pullback(initial, 0.7, -0.4, qbar, vbar).unwrap();
            let second = half.pullback(first.state, 0.15, -0.4, qbar, vbar).unwrap();
            let first_bar = half
                .pullback(
                    initial,
                    0.7,
                    0.15,
                    second.inputs.initial_displacement,
                    second.inputs.initial_velocity,
                )
                .unwrap();
            close(
                first_bar.stiffness + second.stiffness,
                whole.stiffness,
                3e-12,
            );
            close(
                first_bar.damping_rate + second.damping_rate,
                whole.damping_rate,
                3e-12,
            );
        }
    }
}

fn objective(k: f64, c: f64, initial: ModalAcousticState, forces: &[f64]) -> f64 {
    let flow = LinearForceTransition::try_new_with_damping_rate(k.sqrt(), c, 0.13).unwrap();
    let mut current = initial;
    for force in forces.windows(2) {
        current = flow.advance(current, force[0], force[1]).unwrap();
    }
    0.6 * current.displacement_m_sqrt_kg - 0.3 * current.velocity_m_sqrt_kg_per_s
}

fn difference(function: impl Fn(f64) -> f64, value: f64) -> f64 {
    let epsilon = 1e-4;
    (function(value - 2.0 * epsilon) - 8.0 * function(value - epsilon)
        + 8.0 * function(value + epsilon)
        - function(value + 2.0 * epsilon))
        / (12.0 * epsilon)
}

#[test]
fn g3_transient_damping_calibration_reverse_sweep_matches_independent_perturbations() {
    let (k, c) = (4.0_f64, 0.7);
    let initial = state(0.2, -0.4);
    let forces = [0.4, -0.2, 0.7, 1.1, -0.8, 0.3];
    let original = LinearForceTransition::try_new_with_damping_rate(k.sqrt(), c, 0.13).unwrap();
    let flow = original.prepare_sensitivities().unwrap();
    let mut history = vec![initial];
    for force in forces.windows(2) {
        history.push(
            original
                .advance(*history.last().unwrap(), force[0], force[1])
                .unwrap(),
        );
    }
    let (mut qbar, mut vbar) = (0.6, -0.3);
    let (mut kbar, mut cbar) = (0.0, 0.0);
    let mut fbar = [0.0; 6];
    for interval in (0..forces.len() - 1).rev() {
        let gradient = flow
            .pullback(
                history[interval],
                forces[interval],
                forces[interval + 1],
                qbar,
                vbar,
            )
            .unwrap();
        kbar += gradient.stiffness;
        cbar += gradient.damping_rate;
        fbar[interval] += gradient.inputs.force_start;
        fbar[interval + 1] += gradient.inputs.force_end;
        qbar = gradient.inputs.initial_displacement;
        vbar = gradient.inputs.initial_velocity;
    }
    close(
        kbar,
        difference(|value| objective(value, c, initial, &forces), k),
        2e-8,
    );
    close(
        cbar,
        difference(|value| objective(k, value, initial, &forces), c),
        2e-8,
    );
    close(
        qbar,
        difference(|value| objective(k, c, state(value, -0.4), &forces), 0.2),
        2e-8,
    );
    close(
        vbar,
        difference(|value| objective(k, c, state(0.2, value), &forces), -0.4),
        2e-8,
    );
    for i in 0..forces.len() {
        close(
            fbar[i],
            difference(
                |value| {
                    let mut changed = forces;
                    changed[i] = value;
                    objective(k, c, initial, &changed)
                },
                forces[i],
            ),
            2e-8,
        );
    }
}

#[test]
fn g3_ratio_parameters_use_the_physical_rate_chain_rule() {
    let (omega, zeta, h) = (2.5, 0.2, 0.3);
    let initial = state(0.2, -0.4);
    let flow = LinearForceTransition::try_new(omega, zeta, h).unwrap();
    let gradient = flow
        .prepare_sensitivities()
        .unwrap()
        .pullback(initial, 0.7, -0.4, 0.6, -0.3)
        .unwrap();
    let value = |frequency, ratio| {
        let next = LinearForceTransition::try_new(frequency, ratio, h)
            .unwrap()
            .advance(initial, 0.7, -0.4)
            .unwrap();
        0.6 * next.displacement_m_sqrt_kg - 0.3 * next.velocity_m_sqrt_kg_per_s
    };
    close(
        2.0 * omega * gradient.stiffness + 2.0 * zeta * gradient.damping_rate,
        difference(|frequency| value(frequency, zeta), omega),
        2e-8,
    );
    close(
        2.0 * omega * gradient.damping_rate,
        difference(|ratio| value(omega, ratio), zeta),
        2e-8,
    );
    assert_eq!(gradient.inputs, flow.pullback(0.6, -0.3).unwrap());
}

#[test]
fn g0_sensitivity_refusals_are_nonmutating_and_replay_is_deterministic() {
    let original = LinearForceTransition::try_new_with_damping_rate(2.0, 0.7, 0.13).unwrap();
    let flow = original.prepare_sensitivities().unwrap();
    let initial = state(0.2, -0.4);
    let before = flow.advance(initial, 0.7, -0.4).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(flow.advance(state(bad, -0.4), 0.7, -0.4).is_err());
        assert!(flow.advance(initial, bad, -0.4).is_err());
        assert!(flow.pullback(initial, 0.7, -0.4, bad, 0.0).is_err());
    }
    let long = LinearForceTransition::try_new_with_damping_rate(0.0, 0.0, 1e155)
        .unwrap()
        .prepare_sensitivities()
        .unwrap();
    assert!(long.advance(state(1.0, 0.0), 0.0, 0.0).is_err());
    assert_eq!(before, flow.advance(initial, 0.7, -0.4).unwrap());
    let zero = flow.pullback(initial, 0.7, -0.4, 0.0, 0.0).unwrap();
    assert_eq!(zero.stiffness, 0.0);
    assert_eq!(zero.damping_rate, 0.0);
    assert_eq!(zero.inputs, original.pullback(0.0, 0.0).unwrap());
}
