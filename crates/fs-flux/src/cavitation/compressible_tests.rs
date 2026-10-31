use super::*;

const RADIUS: f64 = 50.0e-6;
const AMBIENT: f64 = 101_325.0;

fn liquid(viscosity: f64, sound_speed: f64) -> CavitationState {
    CavitationState::try_new(293.15, 1_000.0, viscosity, 0.072, 2_338.0, sound_speed).unwrap()
}

fn limits() -> KellerMiksisLimits {
    KellerMiksisLimits {
        bubble: BubbleLimits {
            minimum_radius_m: RADIUS * 1e-3,
            maximum_radius_m: RADIUS * 10.0,
            maximum_wall_mach_number: 0.1,
            maximum_vapor_density_ratio: 0.01,
        },
        maximum_pressure_ratio: 0.01,
    }
}

fn model(viscosity: f64, sound_speed: f64) -> KellerMiksisBubble {
    let liquid = liquid(viscosity, sound_speed);
    KellerMiksisBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, RADIUS, AMBIENT, 1.4).unwrap(),
        limits(),
    )
    .unwrap()
}

fn numerics(duration: f64) -> BubbleNumerics {
    BubbleNumerics {
        relative_tolerance: 1e-12,
        radius_absolute_tolerance_m: 1e-18,
        velocity_absolute_tolerance_m_s: 1e-12,
        initial_step_s: duration / 1000.0,
        minimum_step_s: duration * 1e-12,
        maximum_step_s: duration / 100.0,
        maximum_attempts: 10_000,
    }
}

fn constant(pressure: f64) -> BubblePressureSample {
    BubblePressureSample {
        absolute_pressure_pa: pressure,
        time_derivative_pa_s: 0.0,
    }
}

#[test]
fn compressible_equilibrium_has_no_motion_or_radiation_source() {
    let model = model(0.001, 1500.0);
    let mut state = BubbleState::try_new(0.0, RADIUS, 0.0).unwrap();
    let step = model
        .advance(&mut state, 1e-4, &numerics(1e-4), &|_| constant(AMBIENT))
        .unwrap();
    assert_eq!(state.radius_m(), RADIUS);
    assert_eq!(state.wall_velocity_m_s(), 0.0);
    assert_eq!(step.observation.volume_acceleration_m3_s2, 0.0);
    assert_eq!(step.observation.viscous_dissipation_w, 0.0);
    assert_eq!(step.rhs_evaluations, 1 + 6 * step.accepted_steps);
}

// At a pressure-balanced motionless interface, only the prescribed pressure
// derivative acts. Differentiating viscous stress gives the extra 4 mu/R0 term
// in the denominator below; dropping it is strongly visible in the second case.
#[test]
fn analytic_drive_derivative_includes_viscous_acceleration_inertia() {
    let rate = 3e11;
    for viscosity in [0.0, 10.0] {
        let model = model(viscosity, 1500.0);
        let state = BubbleState::try_new(0.0, RADIUS, 0.0).unwrap();
        let observation = model
            .observe(
                &state,
                BubblePressureSample {
                    absolute_pressure_pa: AMBIENT,
                    time_derivative_pa_s: rate,
                },
            )
            .unwrap();
        let expected = -rate / (1_000.0 * 1500.0 + 4.0 * viscosity / RADIUS);
        assert!((observation.wall_acceleration_m_s2 / expected - 1.0).abs() < 3e-15);
        assert!(observation.volume_acceleration_m3_s2 < 0.0);
    }
}

#[test]
fn large_sound_speed_recovers_rp_acceleration_and_physical_trajectory() {
    let reference = model(0.001, 1500.0);
    let rp = RayleighPlessetBubble::try_new(reference.liquid(), reference.gas(), limits().bubble)
        .unwrap();
    let state = BubbleState::try_new(0.0, 0.8 * RADIUS, -0.4).unwrap();
    let pressure = BubblePressureSample {
        absolute_pressure_pa: 115_000.0,
        time_derivative_pa_s: 4e9,
    };
    let incompressible = rp
        .observe(&state, pressure.absolute_pressure_pa)
        .unwrap()
        .wall_acceleration_m_s2;
    let mut prior_error = f64::INFINITY;
    for c in [15_000.0, 150_000.0, 1_500_000.0] {
        let compressible = model(0.001, c)
            .observe(&state, pressure)
            .unwrap()
            .wall_acceleration_m_s2;
        let error = (compressible - incompressible).abs();
        assert!(error < 0.102 * prior_error);
        prior_error = error;
    }
    let almost_incompressible = model(0.001, 1.5e12);
    let mut a = BubbleState::try_new(0.0, 0.9 * RADIUS, -0.05).unwrap();
    let mut b = a;
    rp.advance(&mut a, 2e-6, &numerics(2e-6), &|_| AMBIENT)
        .unwrap();
    almost_incompressible
        .advance(&mut b, 2e-6, &numerics(2e-6), &|_| constant(AMBIENT))
        .unwrap();
    assert!((a.radius_m() - b.radius_m()).abs() / RADIUS < 1e-10);
    assert!((a.wall_velocity_m_s() - b.wall_velocity_m_s()).abs() < 1e-8);
}

// Closed-form linearized free motion, reconstructed from the physical inputs
// rather than obtained by integrating another version of the same ODE.
#[test]
fn small_amplitude_free_motion_has_physical_radiation_and_viscous_damping() {
    for (viscosity, c) in [(0.0, 1500.0), (0.02, 300.0)] {
        let model = model(viscosity, c);
        let pg0 = model.gas().reference_pressure_pa();
        let omega0_squared = (3.0 * 1.4 * pg0 - 2.0 * 0.072 / RADIUS) / (1000.0 * RADIUS * RADIUS);
        let a = 1.0 + 4.0 * viscosity / (1000.0 * c * RADIUS);
        let b = 4.0 * viscosity / (1000.0 * RADIUS * RADIUS) + omega0_squared * RADIUS / c;
        let beta = b / (2.0 * a);
        let omega_d = (omega0_squared / a - beta * beta).sqrt();
        let time = 7.3 * PI / omega_d;
        let amplitude = RADIUS * 1e-6;
        let decay = (-beta * time).exp();
        let angle = omega_d * time;
        let expected_x = amplitude * decay * (angle.cos() + beta / omega_d * angle.sin());
        let expected_v = -amplitude * decay * omega0_squared / (a * omega_d) * angle.sin();
        let mut state = BubbleState::try_new(0.0, RADIUS + amplitude, 0.0).unwrap();
        model
            .advance(&mut state, time, &numerics(time), &|_| constant(AMBIENT))
            .unwrap();
        assert!((state.radius_m() - RADIUS - expected_x).abs() / amplitude < 3e-5);
        assert!(
            (state.wall_velocity_m_s() - expected_v).abs() / (amplitude * omega0_squared.sqrt())
                < 3e-5
        );
        let coefficients = model.linearization().unwrap();
        assert!((coefficients.inertia_multiplier / a - 1.0).abs() < 3e-15);
        assert!((coefficients.damping_coefficient_per_s / b - 1.0).abs() < 3e-15);
        assert_eq!(coefficients.equilibrium_radius_m, RADIUS);
        assert_eq!(coefficients.equilibrium_absolute_pressure_pa, AMBIENT);
        assert_eq!(coefficients.liquid_density_kg_m3, 1000.0);
        assert_eq!(coefficients.pressure_derivative_time_s, RADIUS / c);
    }
}

// For exp(+i omega t), X/P = -(1+i omega R0/c) /
// [rho R0 (omega0^2-a omega^2+i b omega)]. The complex numerator is essential:
// a model that adds radiation damping but drops the drive derivative fails phase.
#[test]
fn forced_harmonic_amplitude_and_phase_match_independent_linear_response() {
    let viscosity = 0.001;
    let c = 300.0;
    let model = model(viscosity, c);
    let omega0_squared = (3.0 * 1.4 * model.gas().reference_pressure_pa() - 2.0 * 0.072 / RADIUS)
        / (1000.0 * RADIUS * RADIUS);
    let a = 1.0 + 4.0 * viscosity / (1000.0 * c * RADIUS);
    let b = 4.0 * viscosity / (1000.0 * RADIUS * RADIUS) + omega0_squared * RADIUS / c;
    let omega = 0.7 * omega0_squared.sqrt();
    let amplitude = 0.1;
    let nr = -amplitude / (1000.0 * RADIUS);
    let ni = nr * omega * RADIUS / c;
    let dr = omega0_squared - a * omega * omega;
    let di = b * omega;
    let denominator = dr * dr + di * di;
    let xr = (nr * dr + ni * di) / denominator;
    let xi = (ni * dr - nr * di) / denominator;
    let displacement_scale = xr.hypot(xi);
    let time = 0.37 * 2.0 * PI / omega;
    let mut state = BubbleState::try_new(0.0, RADIUS + xr, -omega * xi).unwrap();
    let mut controls = numerics(time);
    controls.relative_tolerance = 1e-13;
    controls.velocity_absolute_tolerance_m_s = 1e-13;
    model
        .advance(&mut state, time, &controls, &|t| BubblePressureSample {
            absolute_pressure_pa: AMBIENT + amplitude * (omega * t).cos(),
            time_derivative_pa_s: -amplitude * omega * (omega * t).sin(),
        })
        .unwrap();
    let angle = omega * time;
    let expected_x = xr * angle.cos() - xi * angle.sin();
    let expected_v = -omega * (xr * angle.sin() + xi * angle.cos());
    assert!((state.radius_m() - RADIUS - expected_x).abs() / displacement_scale < 5e-5);
    assert!((state.wall_velocity_m_s() - expected_v).abs() / (omega * displacement_scale) < 5e-5);
}

#[test]
fn compressible_model_refuses_invalid_forcing_and_large_pressure_at_zero_mach() {
    let model = model(0.001, 1500.0);
    let initial = BubbleState::try_new(0.0, RADIUS, 0.0).unwrap();
    let mut state = initial;
    assert!(matches!(
        model.advance(&mut state, 1e-6, &numerics(1e-6), &|_| constant(
            AMBIENT + 3e7
        )),
        Err(CavitationError::PressureRatioExceeded { .. })
    ));
    assert_eq!(state, initial);
    assert!(
        model
            .advance(&mut state, 1e-6, &numerics(1e-6), &|_| {
                BubblePressureSample {
                    absolute_pressure_pa: AMBIENT,
                    time_derivative_pa_s: f64::NAN,
                }
            })
            .is_err()
    );
    assert_eq!(state, initial);
    assert!(matches!(
        model.observe(
            &BubbleState::try_new(0.0, RADIUS, 151.0).unwrap(),
            constant(AMBIENT)
        ),
        Err(CavitationError::WallMachExceeded { .. })
    ));
    for bad_ratio in [0.0, 0.100001, f64::NAN] {
        let mut bad = limits();
        bad.maximum_pressure_ratio = bad_ratio;
        assert!(KellerMiksisBubble::try_new(model.liquid(), model.gas(), bad).is_err());
    }
    assert!(BubblePressureSample::try_new(AMBIENT, f64::INFINITY).is_err());
    let mut budget = numerics(1e-6);
    budget.maximum_attempts = 1;
    assert!(matches!(
        model.advance(&mut state, 1e-6, &budget, &|_| constant(AMBIENT)),
        Err(CavitationError::StepBudgetExceeded { .. })
    ));
    assert_eq!(state, initial);
}

#[test]
fn compressible_cancellation_and_retry_preserve_the_complete_physical_history() {
    let model = model(0.001, 1500.0);
    let initial = BubbleState::try_new(0.0, RADIUS * 1.01, -0.1).unwrap();
    let controls = numerics(2e-6);
    let mut expected = initial;
    let mut polls = 0;
    let accepted = model
        .advance_controlled(
            &mut expected,
            2e-6,
            &controls,
            &|_| constant(AMBIENT),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut state = initial;
        let mut calls = 0;
        assert_eq!(
            model.advance_controlled(
                &mut state,
                2e-6,
                &controls,
                &|_| constant(AMBIENT),
                &mut || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(CavitationError::Cancelled)
        );
        assert_eq!(state, initial);
        let replay = model
            .advance(&mut state, 2e-6, &controls, &|_| constant(AMBIENT))
            .unwrap();
        assert_eq!(state, expected);
        assert_eq!(replay, accepted);
    }
}

#[test]
fn linearization_refuses_unstable_or_out_of_domain_reference_equilibria() {
    let liquid = liquid(0.001, 1500.0);
    let vapor_only = BubbleGas::try_new(RADIUS, 0.0, 1.0).unwrap();
    let unstable = KellerMiksisBubble::try_new(liquid, vapor_only, limits()).unwrap();
    assert!(unstable.linearization().is_err());
    let out_of_domain = BubbleGas::try_new(RADIUS * 11.0, AMBIENT, 1.4).unwrap();
    let model = KellerMiksisBubble::try_new(liquid, out_of_domain, limits()).unwrap();
    assert!(matches!(
        model.linearization(),
        Err(CavitationError::RadiusOutsideDomain { .. })
    ));
}
