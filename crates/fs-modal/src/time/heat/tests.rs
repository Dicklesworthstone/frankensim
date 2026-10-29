use super::super::*;

fn mode(omega: f64, zeta: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: zeta,
        pressure_per_modal_velocity: C64::ZERO,
    }
}
fn model(m: ModalAcousticMode, state: ModalAcousticState) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![m],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model.restore_states(&[state]).unwrap();
    model
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:.17e}, expected={expected:.17e}, relative={:.3e}",
        (actual - expected).abs() / expected.abs()
    );
}

#[test]
fn high_q_heat_survives_energy_subtraction_roundoff() {
    let omega = 3.0;
    let zeta = 1.0e-20;
    let dt = 0.05;
    let mut body = model(
        mode(omega, zeta),
        ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 1.0,
        },
    );
    let frame = body.step_duration(&[0.0], dt).unwrap();
    // To first order in zeta the undamped velocity is cos(omega*t).
    // The omitted relative correction is O(zeta*omega*dt), <2e-21 here.
    // det-ok: independent test-only analytic oracle, outside production evolution.
    let expected = 2.0 * zeta * omega * (0.5 * dt + (2.0 * omega * dt).sin() / (4.0 * omega));
    relative(frame.viscous_dissipation_j, expected, 2e-13);
    assert!(frame.viscous_dissipation_j < f64::EPSILON * frame.total_modal_energy_j);
    assert!(frame.viscous_dissipation_j > 1e10 * frame.dissipation_roundoff_tolerance_j);
    assert_eq!(
        frame.viscous_dissipation_j,
        frame.modal_viscous_losses[0].energy_j
    );
}

#[test]
fn soft_forced_mode_does_not_underflow_its_physical_heat_factor() {
    let mut body = model(mode(1e-200, 0.3), ModalAcousticState::default());
    let frame = body.step_duration(&[2.0], 1.0).unwrap();
    // q''=2 to working precision, so integral(c*v^2)=c*4/3.
    relative(frame.viscous_dissipation_j, 8e-201, 5e-13);
    relative(body.states()[0].velocity_m_sqrt_kg_per_s, 2.0, 1e-14);
    assert!(frame.viscous_dissipation_j > 0.0);
}

fn independent_rk4(
    m: ModalAcousticMode,
    initial: ModalAcousticState,
    force: f64,
    duration: f64,
) -> [f64; 4] {
    // Independent q/v/heat/work equations, with no exact-transition helpers.
    let c = 2.0 * m.damping_ratio * m.angular_frequency_rad_s;
    let k = m.angular_frequency_rad_s * m.angular_frequency_rad_s;
    let rate = |x: [f64; 4]| {
        [
            x[1],
            force - c * x[1] - k * x[0],
            c * x[1] * x[1],
            force * x[1],
        ]
    };
    let add = |x: [f64; 4], a: [f64; 4], h: f64| std::array::from_fn(|i| x[i] + h * a[i]);
    let mut x = [
        initial.displacement_m_sqrt_kg,
        initial.velocity_m_sqrt_kg_per_s,
        0.0,
        0.0,
    ];
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let a = rate(x);
        let b = rate(add(x, a, 0.5 * h));
        let c = rate(add(x, b, 0.5 * h));
        let d = rate(add(x, c, h));
        for i in 0..4 {
            x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    x
}

#[test]
fn forced_heat_and_work_agree_with_independent_velocity_space_integration() {
    for omega in [0.2_f64, 7.0, 350.0] {
        for zeta in [0.01_f64, 0.3, 1.0 - 1e-10, 1.0, 1.0 + 1e-10, 4.0, 100.0] {
            let m = mode(omega, zeta);
            let initial = ModalAcousticState {
                displacement_m_sqrt_kg: 0.1 / omega.max(1.0),
                velocity_m_sqrt_kg_per_s: -0.15,
            };
            let dt = 0.07_f64.min(0.4 / (omega * zeta.max(1.0)));
            for force in [-0.7, 0.0, 0.9] {
                let mut body = model(m, initial);
                let before = modal_energy(m, initial);
                let frame = body.step_duration(&[force], dt).unwrap();
                let expected = independent_rk4(m, initial, force, dt);
                relative(frame.viscous_dissipation_j, expected[2], 2e-10);
                assert!((frame.input_work_j - expected[3]).abs() < 1e-13);
                assert!((body.states()[0].velocity_m_sqrt_kg_per_s - expected[1]).abs() < 1e-12);
                assert!(
                    (frame.input_work_j
                        - (frame.total_modal_energy_j - before)
                        - frame.viscous_dissipation_j)
                        .abs()
                        < 1e-12
                );
            }
        }
    }
}

#[test]
fn long_oscillatory_and_stiff_intervals_compose_resolved_physical_heat() {
    for (zeta, dt) in [(0.01, 1000.0), (1.0, 40.0), (4.0, 200.0), (1e12, 1000.0)] {
        let m = mode(3.0, zeta);
        let initial = ModalAcousticState {
            displacement_m_sqrt_kg: 0.2,
            velocity_m_sqrt_kg_per_s: 0.3,
        };
        let mut whole = model(m, initial);
        let mut split = model(m, initial);
        let frame = whole.step_duration(&[0.0], dt).unwrap();
        let mut heat = 0.0;
        for _ in 0..16 {
            heat += split
                .step_duration(&[0.0], dt / 16.0)
                .unwrap()
                .viscous_dissipation_j;
        }
        relative(frame.viscous_dissipation_j, heat, 2e-11);
        // Independent overdamped roots give the long, stiff limit without
        // resolving the vanishingly short fast transient by time marching.
        if zeta > 1.0 {
            // det-ok: independent test-only analytic oracle, outside production evolution.
            let root = (zeta * zeta - 1.0).sqrt();
            let slow = -3.0 / (zeta + root);
            let fast = -3.0 * (zeta + root);
            let a = (initial.velocity_m_sqrt_kg_per_s - fast * initial.displacement_m_sqrt_kg)
                / (slow - fast);
            let b = (slow * initial.displacement_m_sqrt_kg - initial.velocity_m_sqrt_kg_per_s)
                / (slow - fast);
            // det-ok: independent test-only analytic oracle, outside production evolution.
            let q = a * (slow * dt).exp() + b * (fast * dt).exp();
            // det-ok: independent test-only analytic oracle, outside production evolution.
            let v = slow * a * (slow * dt).exp() + fast * b * (fast * dt).exp();
            let expected = modal_energy(m, initial) - 0.5 * (v * v + 9.0 * q * q);
            relative(frame.viscous_dissipation_j, expected, 2e-11);
        }
    }
}

#[test]
fn undamped_free_and_static_modes_do_not_turn_roundoff_into_heat() {
    let initial = ModalAcousticState {
        displacement_m_sqrt_kg: 0.2,
        velocity_m_sqrt_kg_per_s: 0.3,
    };
    let mut undamped = model(mode(17.0, 0.0), initial);
    for force in [0.0, 0.35, -0.9] {
        let frame = undamped.step_duration(&[force], 0.043).unwrap();
        assert_eq!(frame.viscous_dissipation_j, 0.0);
        assert_eq!(frame.modal_viscous_losses[0].roundoff_tolerance_j, 0.0);
    }
    let mut free = ModalAcousticTimeModel::try_free_mass(
        48_000,
        2.0,
        0.2,
        0.3,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    assert_eq!(
        free.step_duration(&[-1.0], 0.15)
            .unwrap()
            .viscous_dissipation_j,
        0.0
    );
    for zeta in [0.02, 1.0, 25.0] {
        let mut equilibrium = model(
            mode(2.0, zeta),
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.25,
                velocity_m_sqrt_kg_per_s: 0.0,
            },
        );
        let frame = equilibrium.step_duration(&[1.0], 0.2).unwrap();
        assert_eq!(frame.viscous_dissipation_j, 0.0);
    }
}

#[test]
fn cached_full_half_and_arbitrary_steps_route_the_same_small_heat() {
    let m = mode(37.0, 1e-18);
    let initial = ModalAcousticState {
        displacement_m_sqrt_kg: 0.1,
        velocity_m_sqrt_kg_per_s: 0.2,
    };
    let mut full = model(m, initial);
    let mut split = model(m, initial);
    let mut arbitrary = model(m, initial);
    let mut workspace = ModalAcousticWorkspace::new(&full);
    let whole = full
        .step_into(&[0.7], &mut workspace)
        .unwrap()
        .viscous_dissipation_j;
    let a = split
        .step_duration(&[0.7], 0.5 / 48_000.0)
        .unwrap()
        .viscous_dissipation_j;
    let b = split
        .step_duration(&[0.7], 0.5 / 48_000.0)
        .unwrap()
        .viscous_dissipation_j;
    let c = arbitrary
        .step_duration(&[0.7], 0.3 / 48_000.0)
        .unwrap()
        .viscous_dissipation_j;
    let d = arbitrary
        .step_duration(&[0.7], 0.7 / 48_000.0)
        .unwrap()
        .viscous_dissipation_j;
    relative(whole, a + b, 1e-12);
    relative(whole, c + d, 1e-12);
    let before = full.states().to_vec();
    assert!(full.step_into(&[f64::INFINITY], &mut workspace).is_err());
    assert_eq!(full.states(), before);
    let accepted = full.step_into(&[0.2], &mut workspace).unwrap();
    assert!(accepted.viscous_dissipation_j > 0.0);
}
