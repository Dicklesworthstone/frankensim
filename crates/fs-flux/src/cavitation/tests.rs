use super::*;

fn liquid(mu: f64, sigma: f64, vapor: f64) -> CavitationState {
    CavitationState::try_new(293.15, 1_000.0, mu, sigma, vapor, 1_500.0).unwrap()
}

fn limits(radius: f64) -> BubbleLimits {
    BubbleLimits {
        minimum_radius_m: radius * 1.0e-3,
        maximum_radius_m: radius * 10.0,
        maximum_wall_mach_number: 0.1,
        maximum_vapor_density_ratio: 0.01,
    }
}

fn numerics(time_scale: f64, radius_scale: f64) -> BubbleNumerics {
    BubbleNumerics {
        relative_tolerance: 1.0e-10,
        radius_absolute_tolerance_m: radius_scale * 1.0e-12,
        velocity_absolute_tolerance_m_s: 1.0e-10,
        initial_step_s: time_scale / 100.0,
        minimum_step_s: time_scale * 1.0e-10,
        maximum_step_s: time_scale / 20.0,
        maximum_attempts: 10_000,
    }
}

fn empty_bubble(radius: f64) -> RayleighPlessetBubble {
    RayleighPlessetBubble::try_new(
        liquid(0.0, 0.0, 0.0),
        BubbleGas::try_new(radius, 0.0, 1.0).unwrap(),
        limits(radius),
    )
    .unwrap()
}

// G0: equilibrium includes gas, vapor and capillary pressure. A missed vapor
// pressure or surface-tension factor produces motion in this fixture.
#[test]
fn equilibrium_remains_motionless_and_has_no_acoustic_source() {
    let radius = 50.0e-6;
    let liquid = liquid(0.001, 0.072, 2_338.0);
    let gas = BubbleGas::equilibrium(liquid, radius, 101_325.0, 1.4).unwrap();
    let model = RayleighPlessetBubble::try_new(liquid, gas, limits(radius)).unwrap();
    let mut state = BubbleState::try_new(0.0, radius, 0.0).unwrap();
    let report = model
        .advance(&mut state, 1.0e-4, &numerics(1.0e-4, radius), &|_| {
            101_325.0
        })
        .unwrap();
    assert_eq!(state.radius_m(), radius);
    assert_eq!(state.wall_velocity_m_s(), 0.0);
    assert_eq!(state.time_s(), 1.0e-4);
    assert_eq!(report.observation.volume_rate_m3_s, 0.0);
    assert_eq!(report.observation.volume_acceleration_m3_s2, 0.0);
    assert_eq!(report.observation.viscous_dissipation_w, 0.0);
    assert_eq!(report.rejected_steps, 0);
    assert_eq!(report.rhs_evaluations, 1 + 6 * report.accepted_steps);
}

// G1/G3: independent first integral and collapse-time quadrature of the
// Rayleigh equation. t/(R0 sqrt(rho/delta_p)) at R/R0 = 1/2 and 1/4 is obtained
// by integrating sqrt(3/2) sqrt(x^3/(1-x^3)) from R/R0 to 1, NOT by this ODE
// solver. Rayleigh total-collapse factor is 0.91468135650196245; the tests stop
// before its singular endpoint. See https://arxiv.org/abs/1205.4202, eqs. 3--8.
#[test]
fn rayleigh_collapse_recovers_independent_radius_velocity_and_energy() {
    for (radius, pressure) in [(1.0e-3, 100_000.0), (2.0e-4, 1_000.0)] {
        let model = empty_bubble(radius);
        let time_scale = radius * (1_000.0_f64 / pressure).sqrt();
        for (radius_ratio, time_ratio) in
            [(0.5_f64, 0.8254566798929558), (0.25, 0.8993172639099892)]
        {
            let mut state = BubbleState::try_new(0.0, radius, 0.0).unwrap();
            let step = model
                .advance(
                    &mut state,
                    time_ratio * time_scale,
                    &numerics(time_scale, radius),
                    &|_| pressure,
                )
                .unwrap();
            let expected_velocity =
                -((2.0 * pressure / 3_000.0) * (1.0 / radius_ratio.powi(3) - 1.0)).sqrt();
            assert!((state.radius_m() / radius - radius_ratio).abs() < 2.0e-8);
            assert!((state.wall_velocity_m_s() / expected_velocity - 1.0).abs() < 3.0e-8);
            let initial_energy = pressure * (4.0 * PI / 3.0) * radius.powi(3);
            let energy =
                step.observation.liquid_kinetic_energy_j + pressure * step.observation.volume_m3;
            assert!((energy / initial_energy - 1.0).abs() < 3.0e-9);
            assert!(step.maximum_normalized_local_error <= 1.0);
        }
    }
}

// G1: linearization about a gas/vapor equilibrium gives a damped oscillator
// with beta=2 mu/(rho R0^2) and omega0^2=(3 kappa pg0-2 sigma/R0)/(rho R0^2).
// Independent closed form checks phase, viscous damping, and capillary stiffness.
#[test]
fn small_amplitude_motion_recovers_damped_minnaert_limit() {
    let radius: f64 = 50.0e-6;
    let liquid = liquid(0.02, 0.072, 2_000.0);
    let gas = BubbleGas::equilibrium(liquid, radius, 101_325.0, 1.4).unwrap();
    let model = RayleighPlessetBubble::try_new(liquid, gas, limits(radius)).unwrap();
    let beta = 2.0 * liquid.dynamic_viscosity_pa_s() / (1_000.0 * radius * radius);
    let omega_squared = (3.0 * 1.4 * gas.reference_pressure_pa() - 2.0 * 0.072 / radius)
        / (1_000.0 * radius * radius);
    let omega_d = (omega_squared - beta * beta).sqrt();
    let time = PI / (2.0 * omega_d);
    let amplitude = radius * 1.0e-6;
    let mut state = BubbleState::try_new(0.0, radius + amplitude, 0.0).unwrap();
    let mut controls = numerics(time, radius);
    controls.relative_tolerance = 1.0e-12;
    controls.radius_absolute_tolerance_m = 1.0e-18;
    controls.velocity_absolute_tolerance_m_s = 1.0e-13;
    model
        .advance(&mut state, time, &controls, &|_| 101_325.0)
        .unwrap();
    let expected_displacement = amplitude * (-beta * time).exp() * beta / omega_d;
    let expected_velocity = -amplitude * (-beta * time).exp() * omega_squared / omega_d;
    assert!((state.radius_m() - radius - expected_displacement).abs() / amplitude < 3.0e-5);
    assert!((state.wall_velocity_m_s() / expected_velocity - 1.0).abs() < 3.0e-5);
}

#[test]
fn tensile_absolute_pressure_drives_growth_without_a_pressure_clamp() {
    let radius = 100.0e-6;
    let model = empty_bubble(radius);
    let mut state = BubbleState::try_new(0.0, radius, 0.0).unwrap();
    let observation = model.observe(&state, -10_000.0).unwrap();
    assert!((observation.wall_acceleration_m_s2 - 100_000.0).abs() < 1.0e-9);
    model
        .advance(&mut state, 1.0e-6, &numerics(1.0e-6, radius), &|_| {
            -10_000.0
        })
        .unwrap();
    assert!(state.radius_m() > radius);
    assert!(state.wall_velocity_m_s() > 0.0);
}

// G0: differential mechanical power balance of the exact RP equation,
// independent of the time integrator. No collapse energy is relabeled heat.
#[test]
fn observations_close_mechanical_power_and_geometric_source_identity() {
    let r0 = 100.0e-6;
    let liquid = liquid(0.001, 0.072, 2_338.0);
    let gas = BubbleGas::equilibrium(liquid, r0, 101_325.0, 1.4).unwrap();
    let model = RayleighPlessetBubble::try_new(liquid, gas, limits(r0)).unwrap();
    let state = BubbleState::try_new(0.0, 1.3 * r0, -0.5).unwrap();
    let obs = model.observe(&state, 80_000.0).unwrap();
    let r = state.radius_m();
    let v = state.wall_velocity_m_s();
    let a = obs.wall_acceleration_m_s2;
    let kinetic_rate = 2.0 * PI * 1_000.0 * (3.0 * r * r * v.powi(3) + 2.0 * r.powi(3) * v * a);
    let surface_rate = 8.0 * PI * 0.072 * r * v;
    let pressure_power = (2_338.0 + obs.gas_partial_pressure_pa - 80_000.0) * obs.volume_rate_m3_s;
    assert!(
        (kinetic_rate + surface_rate + obs.viscous_dissipation_w - pressure_power).abs()
            < 1.0e-12 * pressure_power.abs()
    );
    assert!(obs.viscous_dissipation_w > 0.0);
    // Differentiate V(R(t)) independently using symmetric nearby Taylor states.
    let h = 1.0e-8;
    let volume = |radius: f64| (4.0 * PI / 3.0) * radius.powi(3);
    let finite_difference = (volume(r + v * h + 0.5 * a * h * h) - 2.0 * volume(r)
        + volume(r - v * h + 0.5 * a * h * h))
        / (h * h);
    assert!((finite_difference / obs.volume_acceleration_m3_s2 - 1.0).abs() < 2.0e-6);
}

// G4/G5: cancellation at start, mid-integration and final publication never
// leaks a partially advanced state; retrying reproduces the accepted bits.
#[test]
fn cancellation_and_replay_are_transactional() {
    let radius = 1.0e-3;
    let model = empty_bubble(radius);
    let initial = BubbleState::try_new(0.0, radius, 0.0).unwrap();
    let controls = numerics(1.0e-4, radius);
    let end_time = 2.0e-5;
    let mut expected = initial;
    let mut total_polls = 0;
    let expected_report = model
        .advance_controlled(
            &mut expected,
            end_time,
            &controls,
            &|_| 100_000.0,
            &mut || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, total_polls / 2, total_polls - 1, total_polls] {
        let mut actual = initial;
        let mut polls = 0;
        assert_eq!(
            model.advance_controlled(
                &mut actual,
                end_time,
                &controls,
                &|_| 100_000.0,
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(CavitationError::Cancelled)
        );
        assert_eq!(actual, initial);
        let report = model
            .advance(&mut actual, end_time, &controls, &|_| 100_000.0)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(report, expected_report);
    }
}

#[test]
fn physical_boundaries_and_exhausted_work_never_publish_fake_rebounds() {
    let radius = 1.0e-3;
    let initial = BubbleState::try_new(0.0, radius, 0.0).unwrap();
    let controls = numerics(1.0e-4, radius);
    let base = empty_bubble(radius);
    let mut radius_limits = limits(radius);
    radius_limits.minimum_radius_m = 0.9 * radius;
    let radius_model =
        RayleighPlessetBubble::try_new(base.liquid(), base.gas(), radius_limits).unwrap();
    let mut state = initial;
    assert!(matches!(
        radius_model.advance(&mut state, 7.0e-5, &controls, &|_| 100_000.0),
        Err(CavitationError::RadiusOutsideDomain { .. })
    ));
    assert_eq!(state, initial);
    let mut mach_limits = limits(radius);
    mach_limits.maximum_wall_mach_number = 0.001;
    let mach_model =
        RayleighPlessetBubble::try_new(base.liquid(), base.gas(), mach_limits).unwrap();
    assert!(matches!(
        mach_model.advance(&mut state, 7.0e-5, &controls, &|_| 100_000.0),
        Err(CavitationError::WallMachExceeded { .. })
    ));
    assert_eq!(state, initial);
    let mut budget = controls;
    budget.maximum_attempts = 1;
    assert!(matches!(
        base.advance(&mut state, 7.0e-5, &budget, &|_| 100_000.0),
        Err(CavitationError::StepBudgetExceeded { attempts: 1, .. })
    ));
    assert_eq!(state, initial);
}

#[test]
fn invalid_drive_and_tolerances_refuse_before_publication() {
    let radius = 1.0e-3;
    let model = empty_bubble(radius);
    let initial = BubbleState::try_new(0.0, radius, 0.0).unwrap();
    let mut state = initial;
    let mut controls = numerics(1.0e-4, radius);
    assert!(matches!(
        model.advance(&mut state, 1.0e-5, &controls, &|_| f64::NAN),
        Err(CavitationError::InvalidInput("imposed absolute pressure"))
    ));
    assert_eq!(state, initial);
    controls.velocity_absolute_tolerance_m_s = 0.0;
    assert!(matches!(
        model.advance(&mut state, 1.0e-5, &controls, &|_| 100_000.0),
        Err(CavitationError::InvalidInput("bubble integration controls"))
    ));
    assert_eq!(state, initial);
    assert!(BubbleGas::try_new(radius, -1.0, 1.4).is_err());
    assert!(BubbleState::try_new(0.0, 0.0, 0.0).is_err());
    let mut invalid_limits = limits(radius);
    invalid_limits.maximum_wall_mach_number = 0.2;
    assert!(RayleighPlessetBubble::try_new(model.liquid(), model.gas(), invalid_limits).is_err());
}

#[test]
fn supplied_vapor_density_must_satisfy_the_interface_kinematic_approximation() {
    use fs_material::cavitation::SaturatedVaporState;
    let radius = 1.0e-3;
    let gas = BubbleGas::try_new(radius, 1.0e5, 1.4).unwrap();
    let make_liquid = |density| {
        liquid(0.001, 0.072, 2_338.0)
            .with_saturated_vapor(SaturatedVaporState::try_new(density, 2.0e6).unwrap())
            .unwrap()
    };
    assert!(RayleighPlessetBubble::try_new(make_liquid(10.0), gas, limits(radius)).is_ok());
    assert!(matches!(
        RayleighPlessetBubble::try_new(make_liquid(10.001), gas, limits(radius)),
        Err(CavitationError::VaporDensityRatioExceeded { .. })
    ));
    let mut too_permissive = limits(radius);
    too_permissive.maximum_vapor_density_ratio = 0.5;
    assert!(RayleighPlessetBubble::try_new(make_liquid(100.0), gas, too_permissive).is_err());
}
