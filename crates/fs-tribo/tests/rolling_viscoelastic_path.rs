//! G1/G3/G4: original exponential-rate integrals, material work and replay.

use fs_tribo::rolling_loss::viscoelastic::path::{ExponentialPathError, ExponentialRatePath};
use fs_tribo::rolling_loss::viscoelastic::{
    DeformingSphereSurface, QuasistaticSphereLimits, VISCOELASTIC_SPHERE_LAW_ID,
    ViscoelasticSphereParameters, ViscoelasticSphereRollingCard,
};
use fs_tribo::rolling_loss::{
    PatchCurvature, RollingKinematics, RollingLossApplicability, RollingLossChannel,
    RollingLossLaw, RollingLossState, RollingLossStateKind, RollingPatchReceipt,
    RollingWorkOwnership, ViscousContourCard,
};
use fs_tribo::{
    ApplicabilityRange, HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
};

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}; relative={relative:e}, absolute={absolute:e}"
    );
}
fn path(q0: f64, drive: f64, decay: f64, duration: f64) -> ExponentialRatePath {
    ExponentialRatePath::new(q0, drive, decay, duration, || false).unwrap()
}

#[test]
fn linear_limits_and_weak_decay_preserve_cancelled_endpoint_and_first_moment() {
    let p = path(2.0, -3.0, 0.0, 0.5);
    assert_eq!(p.initial_rate(), 2.0);
    assert_eq!(p.drive_rate_per_s(), -3.0);
    assert_eq!(p.decay_rate_per_s(), 0.0);
    assert_eq!(p.duration_s(), 0.5);
    close(p.endpoint_rate(), 0.5, 0.0, 0.0);
    close(p.rate_integral(), 0.625, 2e-15, 0.0);
    close(p.rate_squared_integral(), 0.875, 2e-15, 0.0);
    close(p.integrated_displacement(), 0.1875, 2e-15, 0.0);
    let zero = path(0.0, 0.0, 100.0, 3.0);
    assert_eq!(zero.endpoint_rate(), 0.0);
    assert_eq!(zero.rate_squared_integral(), 0.0);
    let empty = path(f64::MAX, f64::MAX, f64::MAX, 0.0);
    assert_eq!(empty.endpoint_rate(), f64::MAX);
    assert_eq!(empty.rate_integral(), 0.0);
    assert_eq!(empty.rate_squared_integral(), 0.0);
    assert_eq!(empty.integrated_displacement(), 0.0);
    for decay in [1e-12, 1e-20, 1e-100] {
        // q0+drive*T=0. The finite decay correction must survive subtraction
        // that would otherwise round both exp and phi to one.
        let endpoint = path(1.0, -1.0, decay, 1.0);
        close(endpoint.endpoint_rate(), -0.5 * decay, 2e-12, 0.0);
        // The zero-decay signed displacement is zero, but its first nonzero
        // term is -decay/6. Heat remains positive across the reversal.
        let reversal = path(1.0, -2.0, decay, 1.0);
        close(reversal.rate_integral(), -decay / 6.0, 2e-12, 0.0);
        close(reversal.rate_squared_integral(), 1.0 / 3.0, 2e-12, 0.0);
    }
}

// Independent composite Simpson integration of the continuous closed ODE
// solution. Production uses integrated analytic series/Gram and variance forms.
fn quadrature(q0: f64, drive: f64, decay: f64, duration: f64, n: usize) -> [f64; 3] {
    let h = duration / n as f64;
    let mut result = [0.0; 3];
    for i in 0..=n {
        let t = h * i as f64;
        let q = if decay == 0.0 {
            q0 + drive * t
        } else {
            q0 * (-decay * t).exp() + (drive / decay) * -(-decay * t).exp_m1()
        };
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        result[0] += weight * q;
        result[1] += weight * q * q;
        result[2] += weight * (duration - t) * q;
    }
    result.map(|x| x * h / 3.0)
}

#[test]
fn signed_moments_match_independent_time_integrals_through_decay_transition() {
    for z in [0.001, 0.5, 1.0, 1.000000001, 2.0, 12.0, 64.0] {
        let (q0, drive, duration) = (-0.3, 0.7, 0.6);
        let decay = z / duration;
        let p = path(q0, drive, decay, duration);
        let coarse = quadrature(q0, drive, decay, duration, 4096);
        let fine = quadrature(q0, drive, decay, duration, 8192);
        let actual = [
            p.rate_integral(),
            p.rate_squared_integral(),
            p.integrated_displacement(),
        ];
        for i in 0..3 {
            let extrapolated = (16.0 * fine[i] - coarse[i]) / 15.0;
            close(actual[i], extrapolated, 3e-12, 3e-16);
        }
        close(
            p.endpoint_rate() - q0,
            drive * duration - decay * p.rate_integral(),
            4e-13,
            3e-16,
        );
        close(
            p.rate_integral() - q0 * duration,
            0.5 * drive * duration * duration - decay * p.integrated_displacement(),
            4e-13,
            3e-16,
        );
        assert!(p.rate_squared_integral() > 0.0);
    }
    let p = path(-0.3, 0.7, 4.0, 0.6);
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        errors.push((quadrature(-0.3, 0.7, 4.0, 0.6, n)[1] - p.rate_squared_integral()).abs());
    }
    assert!(
        errors[0] / errors[1] > 14.0 && errors[1] / errors[2] > 14.0,
        "{errors:?}"
    );
}

#[test]
fn strong_decay_retains_short_transient_and_balanced_physical_scales() {
    // decay*time overflows, but all requested physical moments remain finite.
    let p = path(1e150, 0.0, 1e300, 1e300);
    assert_eq!(p.endpoint_rate(), 0.0);
    close(p.rate_integral(), 1e-150, 5e-15, 0.0);
    close(p.rate_squared_integral(), 0.5, 5e-15, 0.0);
    close(p.integrated_displacement(), 1e150, 5e-15, 0.0);
    let short = path(1e150, 0.0, 0.0, 1e-300);
    close(short.rate_integral(), 1e-150, 5e-15, 0.0);
    close(short.rate_squared_integral(), 1.0, 5e-15, 0.0);
    let long = path(1e-150, 0.0, 0.0, 1e150);
    close(long.rate_integral(), 1.0, 5e-15, 0.0);
    close(long.rate_squared_integral(), 1e-150, 5e-15, 0.0);
    close(long.integrated_displacement(), 0.5e150, 5e-15, 0.0);
    // Tiny final equilibrium and a large initial transient both contribute.
    let mixed = path(1e150, 1.0, 1e150, 1e150);
    close(mixed.endpoint_rate(), 1e-150, 5e-15, 0.0);
    close(mixed.rate_integral(), 2.0, 5e-15, 0.0);
    close(mixed.rate_squared_integral(), 0.5e150, 5e-15, 0.0);
}

#[test]
fn signed_scaling_and_split_continuation_preserve_all_path_moments() {
    let (q0, drive, decay, a, b) = (0.7, -0.4, 2.3, 0.6, 1.1);
    let all = path(q0, drive, decay, a + b);
    let first = path(q0, drive, decay, a);
    let second = path(first.endpoint_rate(), drive, decay, b);
    assert!(all.endpoint_rate() < 0.0);
    close(all.endpoint_rate(), second.endpoint_rate(), 2e-14, 1e-16);
    close(
        all.rate_integral(),
        first.rate_integral() + second.rate_integral(),
        2e-14,
        1e-16,
    );
    close(
        all.rate_squared_integral(),
        first.rate_squared_integral() + second.rate_squared_integral(),
        2e-14,
        1e-16,
    );
    close(
        all.integrated_displacement(),
        first.integrated_displacement()
            + b * first.rate_integral()
            + second.integrated_displacement(),
        2e-14,
        1e-16,
    );
    for scale in [-3.0, 0.125] {
        let transformed = path(scale * q0, scale * drive, decay, a + b);
        close(
            transformed.endpoint_rate(),
            scale * all.endpoint_rate(),
            2e-14,
            1e-16,
        );
        close(
            transformed.rate_integral(),
            scale * all.rate_integral(),
            2e-14,
            1e-16,
        );
        close(
            transformed.integrated_displacement(),
            scale * all.integrated_displacement(),
            2e-14,
            1e-16,
        );
        close(
            transformed.rate_squared_integral(),
            scale * scale * all.rate_squared_integral(),
            2e-14,
            1e-16,
        );
    }
}

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "sphere/plane",
        "rolling-history",
        "synthetic-material",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}
fn domain() -> RollingLossApplicability {
    RollingLossApplicability::new(
        ApplicabilityRange::new(280.0, 320.0).unwrap(),
        ApplicabilityRange::new(0.0, 100.0).unwrap(),
    )
    .unwrap()
}
fn fixture(tau: f64, indentation: f64) -> (RollingLossLaw, RollingPatchReceipt) {
    let young = 2e7;
    let nu = 0.3;
    let p = ViscoelasticSphereParameters {
        sphere_surface: DeformingSphereSurface::SurfaceA,
        radius_m: 0.02,
        young_modulus_pa: young,
        poisson_ratio: nu,
        density_kg_m3: 1100.0,
        shear_viscosity_pa_s: tau * young / (2.0 * (1.0 + nu)),
        bulk_viscosity_pa_s: tau * young / (3.0 * (1.0 - 2.0 * nu)),
    };
    let card = ViscoelasticSphereRollingCard::new(
        VISCOELASTIC_SPHERE_LAW_ID,
        "synthetic-viscosities",
        InputAuthority::SyntheticFixture,
        p,
        QuasistaticSphereLimits {
            maximum_indentation_radius_ratio: 0.03,
            maximum_speed_wave_ratio: 0.01,
            maximum_relaxation_ratio: 0.05,
        },
        domain(),
    )
    .unwrap();
    let normal = HertzSpherePlane {
        effective_radius: p.radius_m,
        reduced_modulus: young / (1.0 - nu * nu),
    }
    .response(&interface(), indentation)
    .unwrap();
    let patch = RollingPatchReceipt::new(
        "sphere-patch",
        "shared-Hertz",
        "synthetic-normal",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        core::f64::consts::PI * normal.contact_radius_m * normal.contact_radius_m,
        PatchCurvature::Principal {
            first_per_m: 50.0,
            second_per_m: 50.0,
        },
    )
    .unwrap();
    (RollingLossLaw::ViscoelasticSphere(card), patch)
}
fn work(interval: &str) -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "sphere-patch",
        interval,
        "lateral-rolling-angle",
        RollingLossChannel::RollingHysteresis,
    )
    .unwrap()
}

#[test]
fn original_material_path_keeps_positive_reversal_heat_with_zero_net_impulse() {
    let (law, patch) = fixture(2e-4, 2e-4);
    let angular = path(1.0, -2.0, 0.0, 1.0);
    let accepted = RollingLossState::zero();
    let result = law
        .advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.0,
            300.0,
            10.0,
            &work("whole"),
            &accepted,
            || false,
        )
        .unwrap();
    let coefficient = 2e-4 * patch.normal_load_n() * 0.02;
    assert_eq!(result.rolling_moment_impulse_nms, 0.0);
    assert_eq!(result.path.rate_integral(), 0.0);
    close(
        result.rolling.dissipation.total_heat_j,
        coefficient / 3.0,
        3e-15,
        0.0,
    );
    close(
        result.rolling.wrench.rolling_moment_nm,
        coefficient,
        3e-15,
        0.0,
    );
    close(
        result.rolling.generalized_work.endpoint_body_power_w,
        -coefficient,
        3e-15,
        0.0,
    );
    close(
        result.rolling.generalized_work.reconstructed_body_power_w,
        -coefficient / 3.0,
        3e-15,
        0.0,
    );
    assert_eq!(result.rolling.dissipation.contour_deformation_heat_j, 0.0);
    assert_eq!(result.rolling.stored_energy_j, 0.0);
    assert_eq!(
        result.rolling.state,
        RollingLossStateKind::RollingHysteresis
    );
    assert_eq!(result.rolling.next_state.accepted_steps(), 1);
    assert_eq!(
        result.rolling.next_state.cumulative_heat_j(),
        result.rolling.dissipation.total_heat_j
    );
    assert_eq!(accepted, RollingLossState::zero());
    let restored = law
        .restore_checkpoint(
            &patch,
            &interface(),
            &work("next"),
            &result.rolling.checkpoint,
        )
        .unwrap();
    assert_eq!(restored, result.rolling.next_state);
    let mean_rate = RollingKinematics::new(0.0, 0.0, 0.0, 300.0, 10.0, 1.0).unwrap();
    let wrong_mean_heat = law
        .advance(&patch, &interface(), mean_rate, &work("mean"), &accepted)
        .unwrap();
    assert_eq!(wrong_mean_heat.dissipation.total_heat_j, 0.0);
    assert!(result.rolling.dissipation.total_heat_j > 0.0);
}

#[test]
fn material_path_impulse_heat_and_checkpoint_match_independent_quadrature_and_split() {
    let (law, patch) = fixture(2e-4, 2e-4);
    let angular = path(0.8, -0.4, 1.7, 0.9);
    let all = law
        .advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.3,
            300.0,
            10.0,
            &work("whole"),
            &RollingLossState::zero(),
            || false,
        )
        .unwrap();
    let reference = quadrature(0.8, -0.4, 1.7, 0.9, 4096);
    let coefficient = 2e-4 * patch.normal_load_n() * 0.02;
    close(
        all.rolling_moment_impulse_nms,
        -coefficient * reference[0],
        1e-12,
        1e-20,
    );
    close(
        all.rolling.dissipation.total_heat_j,
        coefficient * reference[1],
        1e-12,
        1e-20,
    );
    let first_path = path(0.8, -0.4, 1.7, 0.35);
    let first = law
        .advance_exponential_sphere_path(
            &patch,
            &interface(),
            &first_path,
            0.3,
            300.0,
            10.0,
            &work("first"),
            &RollingLossState::zero(),
            || false,
        )
        .unwrap();
    let continuation = law
        .restore_checkpoint(
            &patch,
            &interface(),
            &work("second"),
            &first.rolling.checkpoint,
        )
        .unwrap();
    let second_path = path(first_path.endpoint_rate(), -0.4, 1.7, 0.55);
    let second = law
        .advance_exponential_sphere_path(
            &patch,
            &interface(),
            &second_path,
            0.3,
            300.0,
            10.0,
            &work("second"),
            &continuation,
            || false,
        )
        .unwrap();
    close(
        all.rolling_moment_impulse_nms,
        first.rolling_moment_impulse_nms + second.rolling_moment_impulse_nms,
        2e-14,
        1e-20,
    );
    close(
        all.rolling.dissipation.total_heat_j,
        second.rolling.next_state.cumulative_heat_j(),
        2e-14,
        1e-20,
    );
    assert_eq!(second.rolling.next_state.accepted_steps(), 2);
    assert!(
        law.restore_checkpoint(
            &patch,
            &InterfaceSystemRef::new(
                "other",
                "rolling-history",
                "synthetic-material",
                InputAuthority::SyntheticFixture,
                InterfaceMedium::Dry
            )
            .unwrap(),
            &work("third"),
            &second.rolling.checkpoint
        )
        .is_err()
    );
}

#[test]
fn source_limits_zero_loss_cancellation_and_failed_retry_are_transactional() {
    let (law, patch) = fixture(2e-4, 2e-4);
    let accepted = RollingLossState::zero();
    let angular = path(0.8, -0.4, 1.7, 0.9);
    let run = |temperature, frequency, path: &ExponentialRatePath| {
        law.advance_exponential_sphere_path(
            &patch,
            &interface(),
            path,
            0.0,
            temperature,
            frequency,
            &work("trial"),
            &accepted,
            || false,
        )
    };
    let good = run(300.0, 10.0, &angular).unwrap();
    assert!(matches!(
        run(400.0, 10.0, &angular),
        Err(ExponentialPathError::Constitutive(_))
    ));
    assert!(matches!(
        run(300.0, 101.0, &angular),
        Err(ExponentialPathError::Constitutive(_))
    ));
    // Final rate is admissible but the initial rate violates the original
    // quasistatic domain. The entire trajectory, not just its endpoint, refuses.
    let too_fast = path(5.0, 0.0, 10.0, 1.0);
    assert!(matches!(
        run(300.0, 10.0, &too_fast),
        Err(ExponentialPathError::Constitutive(_))
    ));
    assert_eq!(
        law.advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.0,
            300.0,
            10.0,
            &work("trial"),
            &accepted,
            || true
        ),
        Err(ExponentialPathError::Cancelled)
    );
    let mut count = 0;
    assert_eq!(
        law.advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.0,
            300.0,
            10.0,
            &work("trial"),
            &accepted,
            || {
                count += 1;
                count == 3
            }
        ),
        Err(ExponentialPathError::Cancelled)
    );
    let wrong = RollingWorkOwnership::new(
        "sphere-patch",
        "trial",
        "contour",
        RollingLossChannel::ContourDeformation,
    )
    .unwrap();
    assert!(
        law.advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.0,
            300.0,
            10.0,
            &wrong,
            &accepted,
            || false
        )
        .is_err()
    );
    let unsupported = RollingLossLaw::ViscousContour(
        ViscousContourCard::new(
            "viscous",
            "synthetic",
            InputAuthority::SyntheticFixture,
            0.1,
            domain(),
        )
        .unwrap(),
    );
    assert!(matches!(
        unsupported.advance_exponential_sphere_path(
            &patch,
            &interface(),
            &angular,
            0.0,
            300.0,
            10.0,
            &wrong,
            &accepted,
            || false
        ),
        Err(ExponentialPathError::Invalid(_))
    ));
    assert_eq!(accepted, RollingLossState::zero());
    assert_eq!(run(300.0, 10.0, &angular).unwrap(), good);
    for (tau, indentation) in [(0.0, 2e-4), (2e-4, 0.0)] {
        let (zero_law, zero_patch) = fixture(tau, indentation);
        let zero = zero_law
            .advance_exponential_sphere_path(
                &zero_patch,
                &interface(),
                &angular,
                0.0,
                300.0,
                10.0,
                &work("zero"),
                &accepted,
                || false,
            )
            .unwrap();
        assert_eq!(zero.rolling_moment_impulse_nms, 0.0);
        assert_eq!(zero.rolling.dissipation.total_heat_j, 0.0);
        assert_eq!(zero.rolling.state, RollingLossStateKind::Quiescent);
    }
    for (q, drive, decay, dt) in [
        (f64::NAN, 0.0, 1.0, 1.0),
        (1.0, f64::INFINITY, 1.0, 1.0),
        (1.0, 0.0, -1.0, 1.0),
        (1.0, 0.0, 1.0, -1.0),
    ] {
        assert!(ExponentialRatePath::new(q, drive, decay, dt, || false).is_err());
    }
    assert_eq!(
        ExponentialRatePath::new(1.0, 0.0, 0.5, 1.0, || true),
        Err(ExponentialPathError::Cancelled)
    );
    let mut count = 0;
    assert_eq!(
        ExponentialRatePath::new(1.0, 0.0, 0.5, 1.0, || {
            count += 1;
            count == 3
        }),
        Err(ExponentialPathError::Cancelled)
    );
    assert!(ExponentialRatePath::new(f64::MAX, 0.0, 0.0, 1.0, || false).is_err());
}
