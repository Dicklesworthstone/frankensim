//! Manufactured source cards; independent catenoid and original-ODE checks.
use super::*;
use fs_contact::capillary_bridge::axisymmetric::*;

fn numerical() -> AxisymmetricNumerics {
    AxisymmetricNumerics {
        quadrature_relative_tolerance: 2e-12,
        gap_relative_tolerance: 2e-12,
        maximum_quadrature_panels: 4096,
        maximum_root_iterations: 80,
    }
}
fn work_budget() -> AxisymmetricWorkBudget {
    AxisymmetricWorkBudget {
        equilibrium: numerical(),
        absolute_tolerance_j: 1e-18,
        relative_tolerance: 1e-6,
    }
}
fn envelope() -> AxisymmetricLimits {
    AxisymmetricLimits {
        minimum_gap_m: 1e-8,
        maximum_gap_m: 0.1,
        maximum_wetted_radius_m: 0.1,
        maximum_suction_pa: 1e8,
        maximum_abs_force_n: 10.0,
        maximum_abs_surface_energy_j: 1.0,
        maximum_abs_step_work_j: 1.0,
        gravity_m_s2: 0.0,
        maximum_bond_number: 0.01,
    }
}
fn branch() -> AxisymmetricBranch {
    AxisymmetricBranch {
        minimum_suction_parameter: 0.0,
        maximum_suction_parameter: 32.0,
        minimum_stability_margin: 0.01,
    }
}
fn make(
    theta: [f64; 2],
    sigma: f64,
    volume: f64,
    young: bool,
) -> SourceAxisymmetricCapillaryBridge {
    let (liquid, a) = sources(theta[0], sigma, 300.0, young);
    let (_, b) = sources(theta[1], sigma, 300.0, young);
    SourceAxisymmetricCapillaryBridge::new(
        liquid,
        [a[0].clone(), b[1].clone()],
        volume,
        branch(),
        envelope(),
    )
    .unwrap()
}

fn temperature_cards3d(
    young: bool,
    density_factor: f64,
) -> (MaterialCard, [InterfaceSystemCard; 2]) {
    let (material, lower) = source_cards(1.1, 0.05, young);
    let (_, upper) = source_cards(1.3, 0.05, young);
    let curves = |original: &ClaimSet| {
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims_ordered() {
            let factor = match original.key.name() {
                "density" => density_factor,
                SURFACE_TENSION_PROPERTY | EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY => 0.8,
                _ => 1.0,
            };
            let mut changed = original.clone();
            let PropertyValue::Scalar { value, dims } = original.value else {
                unreachable!()
            };
            changed.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(300.0, value), (320.0, factor * value)],
                dims,
            };
            changed.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(changed).unwrap();
        }
        claims
    };
    let material =
        MaterialCard::assemble(material.id().clone(), curves(material.claims()), vec![]).unwrap();
    let surfaces = [lower[0].clone(), upper[1].clone()].map(|card| {
        InterfaceSystemCard::assemble(
            card.surface_a().clone(),
            card.surface_b().clone(),
            card.context().clone(),
            curves(card.claims()),
            vec![],
        )
        .unwrap()
    });
    (material, surfaces)
}

#[test]
fn g3_axisymmetric_source_refresh_preserves_mass_gap_and_original_surface_claims() {
    let gate = CancelGate::new_clock_free();
    let hot = QueryPoint::new().with("T", 320.0).unwrap();
    let cold = QueryPoint::new().with("T", 300.0).unwrap();
    for young in [false, true] {
        let (card, surfaces) = temperature_cards3d(young, 0.95);
        let (liquid, wetting) = resolved_sources(&card, &surfaces, 300.0, young);
        let source =
            SourceAxisymmetricCapillaryBridge::new(liquid, wetting, VOLUME, branch(), envelope())
                .unwrap();
        let state = source.solve_at_suction(2.0, numerical(), &gate).unwrap();
        let before = *state.profile();
        assert!(
            source
                .refresh_at_query_point(
                    state,
                    &card,
                    surfaces.each_ref(),
                    &hot,
                    BridgeVolumeConstraint::FixedVolume,
                    work_budget(),
                    &gate
                )
                .is_err()
        );
        let refresh = source
            .refresh_at_query_point(
                state,
                &card,
                surfaces.each_ref(),
                &hot,
                BridgeVolumeConstraint::FixedMass,
                work_budget(),
                &gate,
            )
            .unwrap();
        let after = *refresh.state().profile();
        near(after.liquid_volume_m3, VOLUME / 0.95, 2e-15, 0.0);
        near(
            refresh.source().liquid().density_kg_m3() * after.liquid_volume_m3,
            source.liquid().density_kg_m3() * VOLUME,
            2e-15,
            0.0,
        );
        near(after.gap_m, before.gap_m, 3e-12, 0.0);
        assert_ne!(after.normal_force_n, before.normal_force_n);
        assert_eq!(
            refresh.parameter_work_j(),
            after.surface_free_energy_j - before.surface_free_energy_j
        );
        assert!(
            refresh.estimated_parameter_work_error_j()
                <= work_budget().absolute_tolerance_j
                    + work_budget().relative_tolerance * refresh.parameter_work_j().abs()
        );
        for i in 0..2 {
            assert_eq!(
                refresh.source().wetting()[i].resolved().query_point(),
                &[(String::from("T"), 320.0)]
            );
            assert_eq!(
                refresh.source().wetting()[i].angle_origin(),
                source.wetting()[i].angle_origin()
            );
        }
        let reverse = refresh
            .source()
            .refresh_at_query_point(
                refresh.state(),
                &card,
                surfaces.each_ref(),
                &cold,
                BridgeVolumeConstraint::FixedMass,
                work_budget(),
                &gate,
            )
            .unwrap();
        near(
            reverse.state().profile().liquid_volume_m3,
            VOLUME,
            3e-15,
            0.0,
        );
        near(
            reverse.parameter_work_j(),
            -refresh.parameter_work_j(),
            3e-10,
            1e-12,
        );
        assert_eq!(source.evaluate(state).unwrap(), before);
    }
}

#[test]
fn g0_g4_axisymmetric_refresh_refuses_unresolved_work_branch_exit_and_cancel() {
    let gate = CancelGate::new_clock_free();
    let hot = QueryPoint::new().with("T", 320.0).unwrap();
    for density_factor in [0.95, 1e6] {
        let (card, surfaces) = temperature_cards3d(true, density_factor);
        let (liquid, wetting) = resolved_sources(&card, &surfaces, 300.0, true);
        let source =
            SourceAxisymmetricCapillaryBridge::new(liquid, wetting, VOLUME, branch(), envelope())
                .unwrap();
        let state = source.solve_at_suction(2.0, numerical(), &gate).unwrap();
        let mut budget = work_budget();
        if density_factor < 1.0 {
            budget.absolute_tolerance_j = 1e-40;
            budget.relative_tolerance = 0.0;
        }
        let result = source.refresh_at_query_point(
            state,
            &card,
            surfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedMass,
            budget,
            &gate,
        );
        if density_factor < 1.0 {
            assert!(matches!(
                result,
                Err(AxisymmetricError::WorkUnresolved { .. })
            ));
        } else {
            assert!(matches!(
                result,
                Err(AxisymmetricError::GapOutsideBranch { .. })
            ));
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            source.refresh_at_query_point(
                state,
                &card,
                surfaces.each_ref(),
                &hot,
                BridgeVolumeConstraint::FixedMass,
                work_budget(),
                &cancelled
            ),
            Err(AxisymmetricError::Physical(BridgeError::Cancelled))
        ));
        assert_eq!(source.evaluate(state).unwrap(), *state.profile());
    }
}

#[test]
fn g1_axisymmetric_catenoid_matches_closed_form_and_actual_shape() {
    let gate = CancelGate::new_clock_free();
    for young in [false, true] {
        let theta: [f64; 2] = [1.1, 1.3];
        let model = make(theta, 0.05, VOLUME, young);
        let state = model.solve_at_suction(0.0, numerical(), &gate).unwrap();
        let p = state.profile();
        let u = theta.map(|t| t.cos().atanh());
        let integral: f64 = u.iter().map(|u| 0.5 * u + 0.25 * (2.0 * u).sinh()).sum();
        let neck = (VOLUME / (PI * integral)).cbrt();
        near(p.neck_radius_m, neck, 3e-12, 0.0);
        near(p.gap_m, neck * (u[0] + u[1]), 3e-12, 0.0);
        near(p.neck_height_m, neck * u[0], 3e-12, 0.0);
        near(p.liquid_gas_area_m2, 2.0 * VOLUME / neck, 3e-12, 0.0);
        assert_eq!(p.liquid_minus_gas_pressure_pa, 0.0);
        near(p.normal_force_n, -2.0 * PI * 0.05 * neck, 3e-12, 0.0);
        near(
            p.surface_free_energy_j,
            -0.5 * p.normal_force_n * p.gap_m,
            3e-12,
            0.0,
        );
        for i in 0..2 {
            near(p.contact_radii_m[i], neck * u[i].cosh(), 3e-12, 0.0);
        }
        for i in 0..=20 {
            let sample = model
                .meridian_point(state, f64::from(i) / 20.0, numerical(), &gate)
                .unwrap();
            let z = (sample.height_m - p.neck_height_m) / neck;
            near(sample.radius_m, neck * z.cosh(), 3e-12, 0.0);
            near(sample.radius_axial_derivative, z.sinh(), 3e-12, 1.0);
            near(
                sample.meridional_curvature_per_m + sample.azimuthal_curvature_per_m,
                0.0,
                1e-12,
                1.0,
            );
        }
        let solved = model.solve_at_gap(p.gap_m, numerical(), &gate).unwrap();
        near(solved.profile().gap_m, p.gap_m, 3e-12, 0.0);
    }
}

/// Direct RK4 integration of the original Young--Laplace graph ODE in z/a.
/// The production angle quadrature and algebraic radius are not reused.
fn original_ode(q: f64, height: f64) -> [f64; 4] {
    let derivative = |y: [f64; 4]| {
        let metric = 1.0 + y[1] * y[1];
        [
            y[1],
            metric / y[0] + q * metric * metric.sqrt(),
            y[0] * y[0],
            y[0] * metric.sqrt(),
        ]
    };
    let add = |x: [f64; 4], y: [f64; 4], factor: f64| std::array::from_fn(|i| x[i] + factor * y[i]);
    let step = height / 4096.0;
    let mut y = [1.0, 0.0, 0.0, 0.0];
    for _ in 0..4096 {
        let k1 = derivative(y);
        let k2 = derivative(add(y, k1, 0.5 * step));
        let k3 = derivative(add(y, k2, 0.5 * step));
        let k4 = derivative(add(y, k3, step));
        y = std::array::from_fn(|i| {
            y[i] + step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    y
}

#[test]
fn g1_axisymmetric_original_ode_volume_pressure_force_and_energy_agree() {
    let gate = CancelGate::new_clock_free();
    let theta = [1.1, 1.3];
    let model = make(theta, 0.05, VOLUME, true);
    for q in [0.3, 2.0, 12.0] {
        let state = model.solve_at_suction(q, numerical(), &gate).unwrap();
        let p = state.profile();
        let a = p.neck_radius_m;
        let lower = original_ode(q, p.neck_height_m / a);
        let upper = original_ode(q, (p.gap_m - p.neck_height_m) / a);
        near(PI * a * a * a * (lower[2] + upper[2]), VOLUME, 2e-10, 0.0);
        near(
            2.0 * PI * a * a * (lower[3] + upper[3]),
            p.liquid_gas_area_m2,
            2e-10,
            0.0,
        );
        for (i, oracle) in [lower, upper].into_iter().enumerate() {
            near(a * oracle[0], p.contact_radii_m[i], 2e-10, 0.0);
            near(oracle[1], theta[i].cos() / theta[i].sin(), 2e-10, 1.0);
            let r = p.contact_radii_m[i];
            let force =
                p.liquid_minus_gas_pressure_pa * PI * r * r - 2.0 * PI * r * 0.05 * theta[i].sin();
            near(force, p.normal_force_n, 3e-12, 0.0);
        }
        let energy = 0.05 * 2.0 * PI * a * a * (lower[3] + upper[3])
            + model.surface_laws()[0].wetted_surface_energy_difference_j_m2()
                * PI
                * (a * lower[0]).powi(2)
            + model.surface_laws()[1].wetted_surface_energy_difference_j_m2()
                * PI
                * (a * upper[0]).powi(2);
        near(energy, p.surface_free_energy_j, 2e-10, 1e-12);
        // Homogeneity at fixed Young angles: 2E=-F*h+3*Delta_p*V.
        near(
            2.0 * p.surface_free_energy_j,
            -p.normal_force_n * p.gap_m + 3.0 * p.liquid_minus_gas_pressure_pa * VOLUME,
            3e-11,
            1e-12,
        );
        let response = p.point_port([0.0, 0.0, 1.0], [0.02, 0.0, 0.0]).unwrap();
        assert_eq!(
            response.action_force_n,
            response.reaction_force_n.map(|x| -x)
        );
        assert_eq!(response.tangent_n_per_m, -p.force_gap_derivative_n_per_m);
        near(
            response.action_moment_n_m[1],
            -0.02 * p.normal_force_n,
            1e-14,
            0.0,
        );
    }
}

#[test]
fn g3_axisymmetric_energy_derivatives_and_work_conjugate_continuation() {
    let gate = CancelGate::new_clock_free();
    let model = make([1.1, 1.3], 0.05, VOLUME, true);
    let state = model.solve_at_suction(2.0, numerical(), &gate).unwrap();
    let p = state.profile();
    let delta = 1e-4 * p.gap_m;
    let lo = model
        .solve_at_gap(p.gap_m - delta, numerical(), &gate)
        .unwrap();
    let hi = model
        .solve_at_gap(p.gap_m + delta, numerical(), &gate)
        .unwrap();
    near(
        -(hi.profile().surface_free_energy_j - lo.profile().surface_free_energy_j) / (2.0 * delta),
        p.normal_force_n,
        3e-7,
        0.0,
    );
    near(
        (hi.profile().normal_force_n - lo.profile().normal_force_n) / (2.0 * delta),
        p.force_gap_derivative_n_per_m,
        3e-7,
        0.0,
    );
    for i in 0..2 {
        near(
            (hi.profile().contact_radii_m[i] - lo.profile().contact_radii_m[i]) / (2.0 * delta),
            p.contact_radius_gap_derivatives[i],
            3e-7,
            1.0,
        );
    }
    let step = model
        .continue_to_gap(state, 1.1 * p.gap_m, work_budget(), &gate)
        .unwrap();
    let end = step.state();
    let dh = end.profile().gap_m - p.gap_m;
    near(
        -step.mean_normal_force_n * dh,
        step.work_on_liquid_j,
        1e-14,
        0.0,
    );
    assert_eq!(
        step.work_on_liquid_j,
        end.profile().surface_free_energy_j - p.surface_free_energy_j
    );
    assert!(
        step.mean_normal_force_n > p.normal_force_n
            && step.mean_normal_force_n < end.profile().normal_force_n
    );
    assert!(step.estimated_work_error_j < 1e-6 * step.work_on_liquid_j);
    let plus = model
        .continue_to_gap(state, end.profile().gap_m + delta, work_budget(), &gate)
        .unwrap();
    let minus = model
        .continue_to_gap(state, end.profile().gap_m - delta, work_budget(), &gate)
        .unwrap();
    near(
        (plus.mean_normal_force_n - minus.mean_normal_force_n) / (2.0 * delta),
        step.mean_force_gap_derivative_n_per_m,
        2e-6,
        0.0,
    );
    let reverse = model
        .continue_to_gap(end, p.gap_m, work_budget(), &gate)
        .unwrap();
    near(reverse.work_on_liquid_j, -step.work_on_liquid_j, 3e-10, 0.0);
    let unchanged = model
        .continue_to_gap(state, p.gap_m, work_budget(), &gate)
        .unwrap();
    assert_eq!(unchanged.state(), state);
    assert_eq!(unchanged.work_on_liquid_j, 0.0);
    assert_eq!(unchanged.mean_normal_force_n, p.normal_force_n);
    let adjacent = f64::from_bits(p.gap_m.to_bits() + 1);
    assert!(matches!(
        model.continue_to_gap(state, adjacent, work_budget(), &gate),
        Err(AxisymmetricError::WorkUnresolved { .. })
    ));
}

#[test]
fn g3_axisymmetric_true_three_dimensional_scaling_and_plate_exchange() {
    let gate = CancelGate::new_clock_free();
    let model = make([1.1, 1.3], 0.05, VOLUME, false);
    let base = model.solve_at_suction(3.0, numerical(), &gate).unwrap();
    let p = base.profile();
    let scaled = make([1.1, 1.3], 0.05, 8.0 * VOLUME, false)
        .solve_at_suction(3.0, numerical(), &gate)
        .unwrap();
    let s = scaled.profile();
    near(s.gap_m, 2.0 * p.gap_m, 3e-13, 0.0);
    near(s.normal_force_n, 2.0 * p.normal_force_n, 3e-13, 0.0);
    near(
        s.surface_free_energy_j,
        4.0 * p.surface_free_energy_j,
        3e-13,
        0.0,
    );
    near(
        s.liquid_minus_gas_pressure_pa,
        0.5 * p.liquid_minus_gas_pressure_pa,
        3e-13,
        0.0,
    );
    let exchanged = make([1.3, 1.1], 0.05, VOLUME, false)
        .solve_at_suction(3.0, numerical(), &gate)
        .unwrap();
    let e = exchanged.profile();
    near(e.normal_force_n, p.normal_force_n, 2e-13, 0.0);
    near(e.surface_free_energy_j, p.surface_free_energy_j, 2e-13, 0.0);
    near(e.contact_radii_m[0], p.contact_radii_m[1], 2e-13, 0.0);
    near(e.neck_height_m, p.gap_m - p.neck_height_m, 2e-13, 0.0);
    let tension = make([1.1, 1.3], 0.1, VOLUME, false)
        .solve_at_suction(3.0, numerical(), &gate)
        .unwrap();
    near(
        tension.profile().normal_force_n,
        2.0 * p.normal_force_n,
        2e-13,
        0.0,
    );
    near(
        tension.profile().surface_free_energy_j,
        2.0 * p.surface_free_energy_j,
        2e-13,
        0.0,
    );
}

#[test]
fn g1_g3_tiny_axisymmetric_motion_resolves_force_work_and_tangent_without_energy_cancellation() {
    let gate = CancelGate::new_clock_free();
    let model = make([1.1, 1.3], 0.05, VOLUME, true);
    let previous = model.solve_at_suction(2.0, numerical(), &gate).unwrap();
    let p = previous.profile();
    for fraction in [-1e-7, -1e-9, 1e-9, 1e-7] {
        let step = model
            .continue_to_gap(previous, (1.0 + fraction) * p.gap_m, work_budget(), &gate)
            .unwrap();
        let end = step.state();
        let dh = end.profile().gap_m - p.gap_m;
        // Independent local Taylor limit of the original endpoint force. Both
        // opening and closing retain the correct finite-path half tangent.
        near(
            step.mean_normal_force_n,
            p.normal_force_n + 0.5 * p.force_gap_derivative_n_per_m * dh,
            3e-11,
            0.0,
        );
        near(
            step.mean_force_gap_derivative_n_per_m,
            0.5 * p.force_gap_derivative_n_per_m,
            3e-7,
            0.0,
        );
        near(
            -step.mean_normal_force_n * dh,
            step.work_on_liquid_j,
            2e-15,
            0.0,
        );
        assert_eq!(step.work_on_liquid_j.signum(), dh.signum());
        assert!(step.estimated_work_error_j < 1e-9 * step.work_on_liquid_j.abs());
        let subtraction_error = p.surface_free_energy_j.abs() * 32.0 * f64::EPSILON;
        assert!(step.estimated_work_error_j < subtraction_error);
        assert_eq!(
            step.surface_energy_balance_residual_j,
            end.profile().surface_free_energy_j - p.surface_free_energy_j - step.work_on_liquid_j
        );
        assert!(
            step.surface_energy_balance_residual_j.abs()
                <= previous.report().surface_energy_error_estimate_j
                    + end.report().surface_energy_error_estimate_j
                    + step.estimated_work_error_j
        );
        let reverse = model
            .continue_to_gap(end, p.gap_m, work_budget(), &gate)
            .unwrap();
        near(
            reverse.mean_normal_force_n,
            step.mean_normal_force_n,
            3e-11,
            0.0,
        );
        assert_eq!(
            previous.profile().liquid_volume_m3,
            end.profile().liquid_volume_m3
        );
    }
}

#[test]
fn g0_g4_axisymmetric_sources_stability_budgets_and_cancellation_are_atomic() {
    let gate = CancelGate::new_clock_free();
    let model = make([1.1, 1.3], 0.05, VOLUME, true);
    let state = model.solve_at_suction(2.0, numerical(), &gate).unwrap();
    let original = state;
    let addresses = model
        .wetting()
        .each_ref()
        .map(|w| w.resolved().properties().as_ptr());
    for surface in model.wetting() {
        assert_eq!(surface.angle_origin(), EquilibriumAngleOrigin::YoungDupre);
        assert_eq!(
            surface.resolved().query_point(),
            model.liquid().material().query_point()
        );
        assert!(
            surface
                .resolved()
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
    }
    let step = model
        .continue_to_gap(state, 1.05 * state.profile().gap_m, work_budget(), &gate)
        .unwrap();
    assert_eq!(
        addresses,
        model
            .wetting()
            .each_ref()
            .map(|w| w.resolved().properties().as_ptr())
    );
    assert_ne!(
        model.wetting()[0].resolved().identity(),
        model.wetting()[1].resolved().identity()
    );
    let other = make([1.1, 1.3], 0.05, 1.01 * VOLUME, true);
    assert!(matches!(
        other.continue_to_gap(state, state.profile().gap_m, work_budget(), &gate),
        Err(AxisymmetricError::StateMismatch)
    ));
    let mut tight = numerical();
    tight.maximum_quadrature_panels = 8;
    assert!(matches!(
        model.solve_at_suction(0.0, tight, &gate),
        Err(AxisymmetricError::QuadratureBudget)
    ));
    tight = numerical();
    tight.maximum_root_iterations = 1;
    assert!(matches!(
        model.solve_at_gap(state.profile().gap_m, tight, &gate),
        Err(AxisymmetricError::RootBudget)
    ));
    let low_angle = make([0.7, 0.7], 0.05, VOLUME, false);
    assert!(matches!(
        low_angle.solve_at_suction(0.0, numerical(), &gate),
        Err(AxisymmetricError::StabilityUnestablished { .. })
    ));
    let mut restricted = branch();
    restricted.minimum_suction_parameter = 16.0;
    restricted.maximum_suction_parameter = 64.0;
    let (liquid, wetting) = sources(0.2, 0.05, 300.0, false);
    let stable =
        SourceAxisymmetricCapillaryBridge::new(liquid, wetting, VOLUME, restricted, envelope())
            .unwrap();
    assert!(stable.solve_at_suction(32.0, numerical(), &gate).is_ok());
    let (liquid, wetting) = sources(1.1, 0.05, 300.0, false);
    let (_, hot) = sources(1.1, 0.05, 310.0, false);
    assert!(
        SourceAxisymmetricCapillaryBridge::new(
            liquid,
            [wetting[0].clone(), hot[1].clone()],
            VOLUME,
            branch(),
            envelope()
        )
        .is_err()
    );
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(matches!(
        model.continue_to_gap(state, step.state().profile().gap_m, work_budget(), &cancel),
        Err(AxisymmetricError::Physical(BridgeError::Cancelled))
    ));
    assert_eq!(state, original);
    assert_eq!(
        step,
        model
            .continue_to_gap(state, 1.05 * state.profile().gap_m, work_budget(), &gate)
            .unwrap()
    );
    let mut limits = envelope();
    limits.maximum_suction_pa = 1.0;
    let (liquid, wetting) = sources(1.1, 0.05, 300.0, false);
    let limited =
        SourceAxisymmetricCapillaryBridge::new(liquid, wetting, VOLUME, branch(), limits).unwrap();
    assert!(limited.solve_at_suction(2.0, numerical(), &gate).is_err());
    assert!(matches!(
        model.solve_at_gap(0.01, numerical(), &gate),
        Err(AxisymmetricError::GapOutsideBranch { .. })
    ));
}

#[test]
fn g1_mobile_contact_line_second_variation_includes_volume_and_wetting_terms() {
    let gate = CancelGate::new_clock_free();
    let model = make([1.1, 1.3], 0.05, VOLUME, true);
    let state = model.solve_at_suction(0.0, numerical(), &gate).unwrap();
    let p = state.profile();
    let h = p.gap_m;
    let a = p.neck_radius_m;
    // Exact catenoid and an asymmetric free-endpoint radial perturbation.
    let radius = |z: f64| a * ((z - p.neck_height_m) / a).cosh();
    let slope = |z: f64| ((z - p.neck_height_m) / a).sinh();
    let offset = simpson(|z| radius(z) * z / h, 0.0, h) / simpson(radius, 0.0, h);
    let u = |z: f64| z / h - offset;
    near(simpson(|z| radius(z) * u(z), 0.0, h), 0.0, 1e-12, a * h);
    let q = 2.0
        * PI
        * 0.05
        * simpson(
            |z| {
                let metric = 1.0 + slope(z).powi(2);
                radius(z) / (metric * metric.sqrt() * h * h)
                    - u(z).powi(2) / (radius(z) * metric.sqrt())
            },
            0.0,
            h,
        );
    let energy = |epsilon: f64| {
        let perturbed_volume = PI * simpson(|z| (radius(z) + epsilon * u(z)).powi(2), 0.0, h);
        // Preserve volume at every finite perturbation, not just to first order.
        let scale = (VOLUME / perturbed_volume).sqrt();
        let area = 2.0
            * PI
            * simpson(
                |z| {
                    scale
                        * (radius(z) + epsilon * u(z))
                        * (1.0 + (scale * (slope(z) + epsilon / h)).powi(2)).sqrt()
                },
                0.0,
                h,
            );
        0.05 * area
            + [0.0, h]
                .into_iter()
                .enumerate()
                .map(|(i, z)| {
                    model.surface_laws()[i].wetted_surface_energy_difference_j_m2()
                        * PI
                        * (scale * (radius(z) + epsilon * u(z))).powi(2)
                })
                .sum::<f64>()
    };
    let epsilon = 1e-4 * a;
    let observed = (energy(epsilon) - 2.0 * energy(0.0) + energy(-epsilon)) / (epsilon * epsilon);
    near(observed, q, 2e-5, 0.0);
    let bound =
        2.0 * PI * 0.05 * p.axisymmetric_stability_margin / a * simpson(|z| u(z).powi(2), 0.0, h);
    assert!(q > bound && bound > 0.0);
}
