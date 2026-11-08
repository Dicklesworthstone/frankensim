use super::*;
use crate::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubblePressureSample,
    KellerMiksisBubble, KellerMiksisLimits, RayleighPlessetBubble,
};
use fs_material::cavitation::CavitationState;

fn nucleus(
    radius: f64,
    sigma: f64,
    gas_pressure: f64,
    kappa: f64,
) -> RayleighPlessetBubble {
    RayleighPlessetBubble::try_new(
        CavitationState::try_new(
            293.15, 1_000.0, 0.001, sigma, 2_338.0, 1_500.0,
        )
        .unwrap(),
        BubbleGas::try_new(radius, gas_pressure, kappa).unwrap(),
        BubbleLimits {
            minimum_radius_m: radius * 0.1,
            maximum_radius_m: radius * 10.0,
            maximum_wall_mach_number: 0.1,
            maximum_vapor_density_ratio: 0.01,
        },
    )
    .unwrap()
}

fn controls() -> EquilibriumControls {
    EquilibriumControls {
        relative_radius_tolerance: 1.0e-11,
        absolute_radius_tolerance_m: 1.0e-18,
        absolute_pressure_tolerance_pa: 1.0e-8,
        maximum_iterations: 128,
    }
}

fn compressible(model: RayleighPlessetBubble) -> KellerMiksisBubble {
    KellerMiksisBubble::try_new(
        model.liquid(),
        model.gas(),
        KellerMiksisLimits {
            bubble: model.limits(),
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}

// G1: independently reduced isothermal Blake formula, without the numerical
// branch solver or its logarithmic critical-radius evaluation.
#[test]
fn isothermal_blake_threshold_recovers_closed_form_and_shared_static_limit() {
    let radius = 10.0e-6_f64;
    let sigma = 0.072;
    let gas_pressure = 100_000.0;
    let model = nucleus(radius, sigma, gas_pressure, 1.0);
    let expected_radius =
        (3.0 * gas_pressure * radius.powi(3) / (2.0 * sigma)).sqrt();
    let expected_pressure = 2_338.0 - 4.0 * sigma / (3.0 * expected_radius);
    let threshold = blake_threshold(&model).unwrap();
    assert!(
        (threshold.critical_radius_m / expected_radius - 1.0).abs() < 1.0e-13
    );
    assert!(
        (threshold.critical_absolute_pressure_pa - expected_pressure).abs()
            < 1.0e-9
    );
    let response =
        static_response(&model, threshold.critical_radius_m).unwrap();
    assert!(response.restoring_pressure_pa.abs() < 1.0e-8);
    assert!(
        (response.gas_partial_pressure_pa / threshold.gas_partial_pressure_pa
            - 1.0)
            .abs()
            < 1.0e-13
    );
    assert_eq!(blake_threshold(&compressible(model)).unwrap(), threshold);
}

// G1: at pg0 = sigma/R0, p_inf = p_v - sigma/R0, the isothermal
// equilibrium polynomial is x^3 - 2x^2 + 1 = (x-1)(x^2-x-1).
// The positive roots are exactly 1 and the golden ratio: two equilibria
// BELOW vapor pressure, with different mechanical stability.
#[test]
fn stable_and_unstable_roots_preserve_one_original_gas_inventory() {
    let radius = 10.0e-6;
    let model = nucleus(radius, 0.072, 7_200.0, 1.0);
    let pressure = 2_338.0 - 7_200.0;
    for (branch, expected_ratio) in [
        (EquilibriumBranch::Stable, 1.0),
        (EquilibriumBranch::Unstable, 0.5 * (1.0 + 5.0_f64.sqrt())),
    ] {
        let result =
            equilibrium_at_pressure(&model, pressure, branch, &controls())
                .unwrap();
        assert!(
            (result.response.radius_m / radius - expected_ratio).abs() < 1.0e-10
        );
        assert!(
            result.pressure_residual_pa.abs()
                <= controls().absolute_pressure_tolerance_pa
        );
        let expected_gas_pressure = 7_200.0 / expected_ratio.powi(3);
        assert!(
            (result.response.gas_partial_pressure_pa / expected_gas_pressure
                - 1.0)
                .abs()
                < 1.0e-10
        );
        let stable = result.response.restoring_pressure_pa > 0.0;
        assert_eq!(stable, branch == EquilibriumBranch::Stable);
        let accepted_state =
            BubbleState::try_new(0.0, result.response.radius_m, 0.0).unwrap();
        assert!(
            model
                .observe(&accepted_state, pressure)
                .unwrap()
                .wall_acceleration_m_s2
                .abs()
                < 1.0e-5
        );
        // Verify stability through the actual native equations, including the
        // compressible model's viscous inertia, independently of the slope.
        for perturbation in [-1.0e-5, 1.0e-5] {
            let state = BubbleState::try_new(
                0.0,
                result.response.radius_m * (1.0 + perturbation),
                0.0,
            )
            .unwrap();
            let acceleration =
                model.observe(&state, pressure).unwrap().wall_acceleration_m_s2;
            let compressible_acceleration = compressible(model)
                .observe(
                    &state,
                    BubblePressureSample::try_new(pressure, 0.0).unwrap(),
                )
                .unwrap()
                .wall_acceleration_m_s2;
            assert_eq!(acceleration * perturbation < 0.0, stable);
            assert_eq!(compressible_acceleration * perturbation < 0.0, stable);
        }
    }
}

// G1/G3: the fold is a pressure minimum for every admitted polytropic exponent.
// Once crossed, a nucleus initialized at Rc grows under the native ODE.
#[test]
fn crossing_blake_threshold_removes_equilibrium_and_drives_nucleus_growth() {
    for kappa in [0.5, 1.0, 1.4, 5.0 / 3.0] {
        let model = nucleus(10.0e-6, 0.072, 7_200.0, kappa);
        let threshold = blake_threshold(&model).unwrap();
        let radius = threshold.critical_radius_m;
        let critical = threshold.critical_absolute_pressure_pa;
        for ratio in [0.99, 1.01] {
            assert!(
                static_response(&model, radius * ratio)
                    .unwrap()
                    .equilibrium_absolute_pressure_pa
                    > critical
            );
        }
        let at_fold = equilibrium_at_pressure(
            &model,
            critical,
            EquilibriumBranch::Stable,
            &controls(),
        )
        .unwrap();
        assert_eq!(at_fold.response.radius_m, radius);
        let pressure = critical - 100.0;
        for branch in [EquilibriumBranch::Stable, EquilibriumBranch::Unstable] {
            assert!(matches!(
                equilibrium_at_pressure(&model, pressure, branch, &controls()),
                Err(InceptionError::EquilibriumOutsideDomain { .. })
            ));
        }
        let mut state = BubbleState::try_new(0.0, radius, 0.0).unwrap();
        let expected_acceleration = 100.0 / (1_000.0 * radius);
        assert!(
            (model.observe(&state, pressure).unwrap().wall_acceleration_m_s2
                / expected_acceleration
                - 1.0)
                .abs()
                < 1.0e-10
        );
        model
            .advance(
                &mut state,
                1.0e-6,
                &BubbleNumerics {
                    relative_tolerance: 1.0e-10,
                    radius_absolute_tolerance_m: 1.0e-18,
                    velocity_absolute_tolerance_m_s: 1.0e-10,
                    initial_step_s: 1.0e-8,
                    minimum_step_s: 1.0e-15,
                    maximum_step_s: 5.0e-8,
                    maximum_attempts: 10_000,
                },
                &|_| pressure,
            )
            .unwrap();
        assert!(state.radius_m() > radius);
        assert!(state.wall_velocity_m_s() > 0.0);
    }
}

#[test]
fn missing_roots_and_folds_do_not_become_pressure_or_radius_clamps() {
    let model = nucleus(10.0e-6, 0.072, 7_200.0, 1.0);
    assert!(matches!(
        equilibrium_at_pressure(
            &model,
            2_338.0,
            EquilibriumBranch::Unstable,
            &controls(),
        ),
        Err(InceptionError::EquilibriumOutsideDomain { .. })
    ));
    let far_endpoint =
        static_response(&model, model.limits().maximum_radius_m).unwrap();
    assert!(matches!(
        equilibrium_at_pressure(
            &model,
            far_endpoint.equilibrium_absolute_pressure_pa + 1.0,
            EquilibriumBranch::Unstable,
            &controls(),
        ),
        Err(InceptionError::EquilibriumOutsideDomain { .. })
    ));
    for (sigma, gas, kappa) in [
        (0.0, 7_200.0, 1.0),
        (0.072, 0.0, 1.0),
        (0.072, 7_200.0, 1.0 / 3.0),
    ] {
        let ideal = nucleus(10.0e-6, sigma, gas, kappa);
        assert_eq!(
            blake_threshold(&ideal),
            Err(InceptionError::NoFiniteBlakeMinimum)
        );
        assert!(static_response(&ideal, 10.0e-6).is_ok());
    }
    let outside = nucleus(10.0e-6, 0.072, 1.0e8, 1.0);
    assert!(matches!(
        blake_threshold(&outside),
        Err(InceptionError::Input(
            CavitationError::RadiusOutsideDomain { .. }
        ))
    ));
}

// The same analytic stable root R0 remains usable when the native model only
// admits a small neighborhood around it, excluding the larger Blake radius.
#[test]
fn stable_equilibrium_inside_native_window_does_not_require_admitting_fold() {
    let radius = 10.0e-6;
    let source = nucleus(radius, 0.072, 7_200.0, 1.0);
    let limits = BubbleLimits {
        minimum_radius_m: 0.9 * radius,
        maximum_radius_m: 1.1 * radius,
        ..source.limits()
    };
    let model =
        RayleighPlessetBubble::try_new(source.liquid(), source.gas(), limits)
            .unwrap();
    assert!(blake_threshold(&source).unwrap().critical_radius_m
        > limits.maximum_radius_m);
    assert!(matches!(
        blake_threshold(&model),
        Err(InceptionError::Input(
            CavitationError::RadiusOutsideDomain { .. }
        ))
    ));
    let result = equilibrium_at_pressure(
        &model,
        -4_862.0,
        EquilibriumBranch::Stable,
        &controls(),
    )
    .unwrap();
    assert!((result.response.radius_m / radius - 1.0).abs() < 1.0e-10);
    assert!(result.response.restoring_pressure_pa > 0.0);
    assert!(result.radius_bracket_m[0] >= limits.minimum_radius_m);
    assert!(result.radius_bracket_m[1] <= limits.maximum_radius_m);
    assert!(matches!(
        equilibrium_at_pressure(
            &model,
            -4_862.0,
            EquilibriumBranch::Unstable,
            &controls(),
        ),
        Err(InceptionError::EquilibriumOutsideDomain { .. })
    ));
    // Both physical endpoints must bound the requested pressure when the
    // native window ends before the fold.
    for pressure in [
        static_response(&model, limits.maximum_radius_m)
            .unwrap()
            .equilibrium_absolute_pressure_pa
            - 1.0,
        static_response(&model, limits.minimum_radius_m)
            .unwrap()
            .equilibrium_absolute_pressure_pa
            + 1.0,
    ] {
        assert!(matches!(
            equilibrium_at_pressure(
                &model,
                pressure,
                EquilibriumBranch::Stable,
                &controls(),
            ),
            Err(InceptionError::EquilibriumOutsideDomain { .. })
        ));
    }
}

// The golden-ratio unstable root remains usable when Rc is below the native
// window. The gas reference remains R0 even though R0 is outside this window.
#[test]
fn unstable_equilibrium_inside_native_window_does_not_require_admitting_fold() {
    let radius = 10.0e-6;
    let source = nucleus(radius, 0.072, 7_200.0, 1.0);
    let limits = BubbleLimits {
        minimum_radius_m: 1.4 * radius,
        maximum_radius_m: 1.8 * radius,
        ..source.limits()
    };
    let model =
        RayleighPlessetBubble::try_new(source.liquid(), source.gas(), limits)
            .unwrap();
    assert!(blake_threshold(&source).unwrap().critical_radius_m
        < limits.minimum_radius_m);
    assert!(matches!(
        blake_threshold(&model),
        Err(InceptionError::Input(
            CavitationError::RadiusOutsideDomain { .. }
        ))
    ));
    let result = equilibrium_at_pressure(
        &model,
        -4_862.0,
        EquilibriumBranch::Unstable,
        &controls(),
    )
    .unwrap();
    let expected_ratio = 0.5 * (1.0 + 5.0_f64.sqrt());
    assert!(
        (result.response.radius_m / radius - expected_ratio).abs() < 1.0e-10
    );
    assert!(result.response.restoring_pressure_pa < 0.0);
    assert!(result.radius_bracket_m[0] >= limits.minimum_radius_m);
    assert!(result.radius_bracket_m[1] <= limits.maximum_radius_m);
    assert!(matches!(
        equilibrium_at_pressure(
            &model,
            -4_862.0,
            EquilibriumBranch::Stable,
            &controls(),
        ),
        Err(InceptionError::EquilibriumOutsideDomain { .. })
    ));
    for pressure in [
        static_response(&model, limits.minimum_radius_m)
            .unwrap()
            .equilibrium_absolute_pressure_pa
            - 1.0,
        static_response(&model, limits.maximum_radius_m)
            .unwrap()
            .equilibrium_absolute_pressure_pa
            + 1.0,
    ] {
        assert!(matches!(
            equilibrium_at_pressure(
                &model,
                pressure,
                EquilibriumBranch::Unstable,
                &controls(),
            ),
            Err(InceptionError::EquilibriumOutsideDomain { .. })
        ));
    }
}

#[test]
fn branch_solve_has_bounded_work_cancellation_and_repeatability() {
    let model = nucleus(10.0e-6, 0.072, 7_200.0, 1.0);
    let pressure = -4_862.0;
    let branch = EquilibriumBranch::Unstable;
    let result =
        equilibrium_at_pressure(&model, pressure, branch, &controls()).unwrap();
    assert_eq!(
        equilibrium_at_pressure(&model, pressure, branch, &controls()).unwrap(),
        result
    );
    let mut one_iteration = controls();
    one_iteration.maximum_iterations = 1;
    assert_eq!(
        equilibrium_at_pressure(&model, pressure, branch, &one_iteration),
        Err(InceptionError::ToleranceNotMet { iterations: 1 })
    );
    let mut calls = 0;
    assert_eq!(
        equilibrium_at_pressure_controlled(
            &model,
            pressure,
            branch,
            &controls(),
            &mut || {
                calls += 1;
                calls == 5
            },
        ),
        Err(InceptionError::Input(CavitationError::Cancelled))
    );
    assert_eq!(calls, 5);
    let mut invalid = controls();
    invalid.absolute_pressure_tolerance_pa = f64::NAN;
    assert!(matches!(
        equilibrium_at_pressure(&model, pressure, branch, &invalid),
        Err(InceptionError::Input(CavitationError::InvalidInput(_)))
    ));
}

#[test]
fn tighter_radius_and_pressure_targets_converge_to_both_analytic_roots() {
    let radius = 10.0e-6;
    let model = nucleus(radius, 0.072, 7_200.0, 1.0);
    let coarse = EquilibriumControls {
        relative_radius_tolerance: 1.0e-3,
        absolute_radius_tolerance_m: 1.0e-12,
        absolute_pressure_tolerance_pa: 1.0,
        maximum_iterations: 128,
    };
    for (branch, expected_ratio) in [
        (EquilibriumBranch::Stable, 1.0),
        (EquilibriumBranch::Unstable, 0.5 * (1.0 + 5.0_f64.sqrt())),
    ] {
        let loose =
            equilibrium_at_pressure(&model, -4_862.0, branch, &coarse).unwrap();
        let tight =
            equilibrium_at_pressure(&model, -4_862.0, branch, &controls())
                .unwrap();
        let exact = radius * expected_ratio;
        assert!(
            (tight.response.radius_m - exact).abs()
                < (loose.response.radius_m - exact).abs() * 1.0e-3
        );
        assert!(tight.iterations > loose.iterations);
        let half_width =
            0.5 * (tight.radius_bracket_m[1] - tight.radius_bracket_m[0]);
        assert!(
            half_width
                <= controls().absolute_radius_tolerance_m
                    + controls().relative_radius_tolerance
                        * tight.response.radius_m
        );
    }
}

#[test]
fn coordinated_length_surface_tension_scaling_preserves_pressures() {
    let reference = nucleus(10.0e-6, 0.072, 7_200.0, 1.4);
    let threshold = blake_threshold(&reference).unwrap();
    let equilibrium = equilibrium_at_pressure(
        &reference,
        -4_862.0,
        EquilibriumBranch::Stable,
        &controls(),
    )
    .unwrap();
    for scale in [1.0e-3, 1.0e3] {
        let scaled =
            nucleus(10.0e-6 * scale, 0.072 * scale, 7_200.0, 1.4);
        let scaled_threshold = blake_threshold(&scaled).unwrap();
        assert!(
            (scaled_threshold.critical_radius_m
                / (threshold.critical_radius_m * scale)
                - 1.0)
                .abs()
                < 1.0e-12
        );
        assert!(
            (scaled_threshold.critical_absolute_pressure_pa
                - threshold.critical_absolute_pressure_pa)
                .abs()
                < 1.0e-8
        );
        let scaled_equilibrium = equilibrium_at_pressure(
            &scaled,
            -4_862.0,
            EquilibriumBranch::Stable,
            &controls(),
        )
        .unwrap();
        assert!(
            (scaled_equilibrium.response.radius_m
                / (equilibrium.response.radius_m * scale)
                - 1.0)
                .abs()
                < 1.0e-10
        );
    }
}
