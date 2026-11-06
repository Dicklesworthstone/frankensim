use super::*;

fn liquid(nu: f64, sigma: f64) -> SlenderJetParameters {
    SlenderJetParameters {
        density_kg_m3: 1.0,
        surface_tension_n_m: sigma,
        kinematic_viscosity_m2_s: nu,
    }
}

fn limits() -> SlenderJetLimits {
    SlenderJetLimits {
        minimum_radius_m: 0.1,
        maximum_abs_slope: 0.75,
        minimum_radius_cells: 2.0,
    }
}

fn budget() -> SlenderJetStepBudget {
    SlenderJetStepBudget {
        maximum_step_s: 0.01,
        cfl: 0.2,
        maximum_substeps: 100_000,
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

fn eigen_profile(n: usize, amplitude: f64, nu: f64) -> (f64, Vec<f64>, Vec<f64>, f64) {
    let k = 0.5;
    let length = 2.0 * core::f64::consts::PI / k;
    let dx = length / n as f64;
    let half_damping = 1.5 * nu * k * k;
    let drive = 0.5 * k * k * (1.0 - k * k);
    let growth = drive / ((half_damping * half_damping + drive).sqrt() + half_damping);
    let radius = (0..n)
        .map(|i| 1.0 + amplitude * fs_math::det::cos(k * i as f64 * dx))
        .collect();
    let velocity = (0..n)
        .map(|i| -2.0 * growth * amplitude / k * fs_math::det::sin(k * (i as f64 + 0.5) * dx))
        .collect();
    (length, radius, velocity, growth)
}

fn radius_mode(jet: &PeriodicSlenderJet) -> f64 {
    let n = jet.radius_squared_m2().len();
    let sum: f64 = jet
        .radii_m()
        .enumerate()
        .map(|(i, radius)| {
            (radius - 1.0) * fs_math::det::cos(2.0 * core::f64::consts::PI * i as f64 / n as f64)
        })
        .sum();
    2.0 * sum / n as f64
}

#[test]
fn g1_uniform_translating_cylinder_is_exactly_stationary_in_shape() {
    let gate = CancelGate::new_clock_free();
    for nu in [0.0, 0.1] {
        let mut jet = PeriodicSlenderJet::new(
            liquid(nu, 1.0),
            8.0,
            &[1.0; 64],
            &[0.375; 64],
            limits(),
            &gate,
        )
        .unwrap();
        let area = jet.radius_squared_m2().to_vec();
        let flux = jet.radius_squared_flux_m3_s().to_vec();
        let before = jet.diagnostics();
        let report = jet.advance(0.1, budget(), &gate).unwrap();
        assert_eq!(report.outcome, SlenderJetOutcome::Completed);
        assert_eq!(jet.radius_squared_m2(), area);
        assert_eq!(jet.radius_squared_flux_m3_s(), flux);
        assert_eq!(jet.diagnostics(), before);
        assert_eq!(jet.time_s(), 0.1);
        near(before.volume_m3, 8.0 * core::f64::consts::PI, 1e-15);
        near(before.mean_velocity_m_s, 0.375, 1e-15);
        near(
            before.axial_momentum_kg_m_s,
            before.volume_m3 * 0.375,
            1e-15,
        );
        assert_eq!(report.viscous_dissipation_j, 0.0);
        assert_eq!(report.upwind_dissipation_j, 0.0);
        assert_eq!(report.energy_balance_residual_j, 0.0);
    }
}

#[test]
fn g1_linear_rayleigh_growth_converges_to_independent_dispersion() {
    let gate = CancelGate::new_clock_free();
    let mut errors = Vec::new();
    for n in [64, 128] {
        let (length, radius, velocity, growth) = eigen_profile(n, 1e-5, 0.1);
        let mut jet = PeriodicSlenderJet::new(
            liquid(0.1, 1.0),
            length,
            &radius,
            &velocity,
            limits(),
            &gate,
        )
        .unwrap();
        let initial = radius_mode(&jet);
        let report = jet.advance(0.5, budget(), &gate).unwrap();
        assert_eq!(report.outcome, SlenderJetOutcome::Completed);
        let observed = fs_math::det::ln(radius_mode(&jet) / initial) / 0.5;
        errors.push((observed - growth).abs());
        assert!(
            (observed - growth).abs() < 0.002 * growth,
            "{observed} versus {growth}"
        );
        near(report.after.volume_m3, report.before.volume_m3, 3e-13);
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
}

#[test]
fn g1_nonlinear_transport_conserves_volume_and_axial_momentum() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, mut velocity, _) = eigen_profile(96, 0.08, 0.08);
    for value in &mut velocity {
        *value += 0.4;
    }
    let mut jet = PeriodicSlenderJet::new(
        liquid(0.08, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let before = jet.diagnostics();
    let initial_area = jet.radius_squared_m2().to_vec();
    let report = jet.advance(0.3, budget(), &gate).unwrap();
    assert_eq!(report.outcome, SlenderJetOutcome::Completed);
    assert_ne!(jet.radius_squared_m2(), initial_area);
    near(report.after.volume_m3, before.volume_m3, 3e-13);
    near(
        report.after.axial_momentum_kg_m_s,
        before.axial_momentum_kg_m_s,
        3e-13,
    );
    assert!(report.viscous_dissipation_j > 0.0);
    assert!(report.upwind_dissipation_j > 0.0);
    assert!(report.energy_balance_residual_j.is_finite());
}

#[test]
fn g1_viscosity_decays_velocity_and_energy_losses_are_separate() {
    let gate = CancelGate::new_clock_free();
    let n = 64;
    let length = 4.0 * core::f64::consts::PI;
    let initial_velocity: Vec<f64> = (0..n)
        .map(|i| {
            0.001 * fs_math::det::sin(2.0 * core::f64::consts::PI * (i as f64 + 0.5) / n as f64)
        })
        .collect();
    let mut jet = PeriodicSlenderJet::new(
        liquid(0.4, 0.0),
        length,
        &vec![1.0; n],
        &initial_velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let initial_energy = jet.diagnostics().kinetic_energy_j;
    let report = jet.advance(0.5, budget(), &gate).unwrap();
    assert_eq!(report.outcome, SlenderJetOutcome::Completed);
    let amplitude = 2.0 / n as f64
        * jet
            .velocities_m_s()
            .enumerate()
            .map(|(i, u)| {
                u * fs_math::det::sin(2.0 * core::f64::consts::PI * (i as f64 + 0.5) / n as f64)
            })
            .sum::<f64>();
    // Linearized Newtonian axial viscosity is 3*nu, not nu.
    near(
        amplitude,
        0.001 * fs_math::det::exp(-3.0 * 0.4 * 0.25 * 0.5),
        3e-4,
    );
    assert!(report.after.kinetic_energy_j < initial_energy);
    assert!(report.viscous_dissipation_j > 0.0);
    assert!(report.upwind_dissipation_j > 0.0);
    assert_eq!(report.after.surface_energy_j, 0.0);
    assert!(report.energy_balance_residual_j.abs() < 1e-6 * initial_energy);
    near(report.after.volume_m3, report.before.volume_m3, 3e-13);
}

#[test]
fn g0_domain_stops_retain_last_admitted_geometry() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.05, 0.05);
    let dx = length / radius.len() as f64;
    for (kind, authored) in [
        (
            SlenderJetLimitKind::MinimumRadius,
            SlenderJetLimits {
                minimum_radius_m: 0.94,
                ..limits()
            },
        ),
        (
            SlenderJetLimitKind::RadiusResolution,
            SlenderJetLimits {
                minimum_radius_cells: 0.94 / dx,
                ..limits()
            },
        ),
        (
            SlenderJetLimitKind::MaximumSlope,
            SlenderJetLimits {
                maximum_abs_slope: 0.028,
                ..limits()
            },
        ),
    ] {
        let mut jet = PeriodicSlenderJet::new(
            liquid(0.05, 1.0),
            length,
            &radius,
            &velocity,
            authored,
            &gate,
        )
        .unwrap();
        let report = jet.advance(1.0, budget(), &gate).unwrap();
        let SlenderJetOutcome::DomainLimit(violation) = report.outcome else {
            panic!("expected domain stop, got {:?}", report.outcome);
        };
        assert_eq!(violation.kind, kind);
        assert!(report.accepted_substeps > 0);
        assert!(report.end_time_s > 0.0 && report.end_time_s < 1.0);
        assert_eq!(jet.time_s(), report.end_time_s);
        assert_eq!(jet.diagnostics(), report.after);
        assert!(report.after.minimum_radius_m >= authored.minimum_radius_m);
        assert!(report.after.minimum_radius_m / dx >= authored.minimum_radius_cells);
        assert!(report.after.maximum_abs_slope <= authored.maximum_abs_slope);
        near(report.after.volume_m3, report.before.volume_m3, 3e-13);
    }
}

#[test]
fn g4_work_budget_and_cancellation_roll_back_the_whole_call() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.04, 0.1);
    let mut jet = PeriodicSlenderJet::new(
        liquid(0.1, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let mut independent = jet.clone();
    let initial_area = jet.radius_squared_m2().to_vec();
    let initial_flux = jet.radius_squared_flux_m3_s().to_vec();
    let initial_diagnostics = jet.diagnostics();
    assert_eq!(
        jet.advance(
            0.1,
            SlenderJetStepBudget {
                maximum_substeps: 1,
                ..budget()
            },
            &gate
        ),
        Err(NonlinearJetError::SubstepBudget)
    );
    assert_eq!(jet.radius_squared_m2(), initial_area);
    assert_eq!(jet.radius_squared_flux_m3_s(), initial_flux);
    assert_eq!(jet.diagnostics(), initial_diagnostics);
    assert_eq!(jet.time_s(), 0.0);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        jet.advance(0.1, budget(), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(jet.radius_squared_m2(), initial_area);
    assert_eq!(jet.radius_squared_flux_m3_s(), initial_flux);
    let retry = jet.advance(0.1, budget(), &gate).unwrap();
    let fresh = independent.advance(0.1, budget(), &gate).unwrap();
    assert_eq!(retry, fresh);
    assert_eq!(jet.radius_squared_m2(), independent.radius_squared_m2());
    assert_eq!(
        jet.radius_squared_flux_m3_s(),
        independent.radius_squared_flux_m3_s()
    );
}

#[test]
fn g3_si_similarity_preserves_the_dimensionless_trajectory() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.04, 0.1);
    let mut unit = PeriodicSlenderJet::new(
        liquid(0.1, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let scale: f64 = 0.001;
    let rho = 1000.0;
    let sigma = 0.072;
    let time_scale = (rho * scale * scale * scale / sigma).sqrt();
    let physical_radius: Vec<f64> = radius.iter().map(|h| scale * h).collect();
    let physical_velocity: Vec<f64> = velocity.iter().map(|u| scale / time_scale * u).collect();
    let mut physical = PeriodicSlenderJet::new(
        SlenderJetParameters {
            density_kg_m3: rho,
            surface_tension_n_m: sigma,
            kinematic_viscosity_m2_s: 0.1 * scale * scale / time_scale,
        },
        scale * length,
        &physical_radius,
        &physical_velocity,
        SlenderJetLimits {
            minimum_radius_m: 0.1 * scale,
            ..limits()
        },
        &gate,
    )
    .unwrap();
    let unit_report = unit.advance(0.15, budget(), &gate).unwrap();
    let physical_report = physical
        .advance(
            0.15 * time_scale,
            SlenderJetStepBudget {
                maximum_step_s: budget().maximum_step_s * time_scale,
                ..budget()
            },
            &gate,
        )
        .unwrap();
    assert_eq!(unit_report.outcome, SlenderJetOutcome::Completed);
    assert_eq!(physical_report.outcome, SlenderJetOutcome::Completed);
    for (a, b) in unit.radii_m().zip(physical.radii_m()) {
        near(b / scale, a, 2e-11);
    }
    for (a, b) in unit.velocities_m_s().zip(physical.velocities_m_s()) {
        assert!((b * time_scale / scale - a).abs() < 2e-11);
    }
    near(
        physical_report.after.volume_m3 / (scale * scale * scale),
        unit_report.after.volume_m3,
        5e-13,
    );
}

#[test]
fn g0_invalid_inputs_and_initial_domain_refuse() {
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        PeriodicSlenderJet::new(
            liquid(-0.1, 1.0),
            8.0,
            &[1.0; 32],
            &[0.0; 32],
            limits(),
            &gate
        ),
        Err(NonlinearJetError::InvalidInput(_))
    ));
    assert!(matches!(
        PeriodicSlenderJet::new(
            liquid(0.1, 1.0),
            8.0,
            &[1.0; 32],
            &[0.0; 31],
            limits(),
            &gate
        ),
        Err(NonlinearJetError::InvalidInput(_))
    ));
    assert!(matches!(
        PeriodicSlenderJet::new(
            liquid(0.1, 1.0),
            8.0,
            &[0.05; 32],
            &[0.0; 32],
            limits(),
            &gate
        ),
        Err(NonlinearJetError::InitialDomain(_))
    ));
    let mut jet = PeriodicSlenderJet::new(
        liquid(0.0, 0.0),
        8.0,
        &[1.0; 32],
        &[0.0; 32],
        limits(),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        jet.advance(f64::NAN, budget(), &gate),
        Err(NonlinearJetError::InvalidInput(_))
    ));
    assert!(matches!(
        jet.advance(
            0.1,
            SlenderJetStepBudget {
                cfl: 1.0,
                ..budget()
            },
            &gate
        ),
        Err(NonlinearJetError::InvalidInput(_))
    ));
    let report = jet.advance(0.0, budget(), &gate).unwrap();
    assert_eq!(report.accepted_substeps, 0);
    assert_eq!(report.before, report.after);
    assert_eq!(jet.time_s(), 0.0);
}
