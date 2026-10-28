//! Independent physical limits of the material-derived rolling deformation law.
//! All material/source data here are synthetic, not specimen calibration.

use fs_tribo::rolling_loss::viscoelastic::{
    DeformingSphereSurface, QuasistaticSphereLimits, VISCOELASTIC_SPHERE_LAW_ID,
    ViscoelasticSphereParameters, ViscoelasticSphereRollingCard,
};
use fs_tribo::rolling_loss::{
    PatchCurvature, RollingKinematics, RollingLossApplicability, RollingLossApplicabilityKind,
    RollingLossChannel, RollingLossLaw, RollingLossState, RollingLossStateKind,
    RollingPatchReceipt, RollingWorkOwnership,
};
use fs_tribo::{
    ApplicabilityRange, HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
};

fn close(actual: f64, expected: f64, relative: f64) {
    let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() <= relative * scale,
        "{actual:e} != {expected:e}"
    );
}

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "sphere-A/rigid-plane-B",
        "rolling-history-1",
        "fixture/material-1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}

fn parameters(nu: f64, shear_time: f64, bulk_time: f64) -> ViscoelasticSphereParameters {
    let young = 2.0e7;
    ViscoelasticSphereParameters {
        sphere_surface: DeformingSphereSurface::SurfaceA,
        radius_m: 0.02,
        young_modulus_pa: young,
        poisson_ratio: nu,
        density_kg_m3: 1100.0,
        shear_viscosity_pa_s: shear_time * (young / (2.0 * (1.0 + nu))),
        bulk_viscosity_pa_s: bulk_time * (young / (3.0 * (1.0 - 2.0 * nu))),
    }
}

fn limits() -> QuasistaticSphereLimits {
    QuasistaticSphereLimits {
        maximum_indentation_radius_ratio: 0.03,
        maximum_speed_wave_ratio: 0.01,
        maximum_relaxation_ratio: 0.05,
    }
}

fn domain() -> RollingLossApplicability {
    RollingLossApplicability::new(
        ApplicabilityRange::new(280.0, 320.0).unwrap(),
        ApplicabilityRange::new(0.0, 100.0).unwrap(),
    )
    .unwrap()
}

fn card(p: ViscoelasticSphereParameters) -> ViscoelasticSphereRollingCard {
    ViscoelasticSphereRollingCard::new(
        VISCOELASTIC_SPHERE_LAW_ID,
        "fixture/viscosities-1",
        InputAuthority::SyntheticFixture,
        p,
        limits(),
        domain(),
    )
    .unwrap()
}

fn patch(p: ViscoelasticSphereParameters, indentation: f64) -> RollingPatchReceipt {
    let normal = HertzSpherePlane {
        effective_radius: p.radius_m,
        reduced_modulus: p.young_modulus_pa / (1.0 - p.poisson_ratio * p.poisson_ratio),
    }
    .response(&interface(), indentation)
    .unwrap();
    RollingPatchReceipt::new(
        "sphere-patch",
        "shared-Hertz",
        "normal-source-1",
        InputAuthority::SyntheticFixture,
        normal.normal_force_n,
        core::f64::consts::PI * normal.contact_radius_m * normal.contact_radius_m,
        PatchCurvature::Principal {
            first_per_m: 1.0 / p.radius_m,
            second_per_m: 1.0 / p.radius_m,
        },
    )
    .unwrap()
}

fn motion(p: ViscoelasticSphereParameters, omega: f64, spin: f64, dt: f64) -> RollingKinematics {
    RollingKinematics::new(p.radius_m * omega, omega, spin, 300.0, 10.0, dt).unwrap()
}

fn ownership() -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "sphere-patch",
        "interval-1",
        "sphere-rolling-angle",
        RollingLossChannel::RollingHysteresis,
    )
    .unwrap()
}

// Independent frequency-domain material oracle: directly form the complex
// Kelvin–Voigt G(w), K(w), then E*(w)=4G(3K+G)/(3K+4G). No rolling coefficient
// or equivalent real-weight formula is reused here.
fn reduced_modulus_logarithmic_loss(p: ViscoelasticSphereParameters, omega: f64) -> f64 {
    let g = [
        p.young_modulus_pa / (2.0 * (1.0 + p.poisson_ratio)),
        omega * p.shear_viscosity_pa_s,
    ];
    let k = [
        p.young_modulus_pa / (3.0 * (1.0 - 2.0 * p.poisson_ratio)),
        omega * p.bulk_viscosity_pa_s,
    ];
    let l = [3.0 * k[0] + g[0], 3.0 * k[1] + g[1]];
    let d = [3.0 * k[0] + 4.0 * g[0], 3.0 * k[1] + 4.0 * g[1]];
    let n = [
        4.0 * (g[0] * l[0] - g[1] * l[1]),
        4.0 * (g[0] * l[1] + g[1] * l[0]),
    ];
    let den = d[0] * d[0] + d[1] * d[1];
    let e = [
        (n[0] * d[0] + n[1] * d[1]) / den,
        (n[1] * d[0] - n[0] * d[1]) / den,
    ];
    e[1] / (omega * e[0])
}

#[test]
fn coefficient_matches_independent_complex_modulus_and_proportional_viscosity() {
    for nu in [-0.8, -0.2, 0.0, 0.3, 0.499] {
        let p = parameters(nu, 4e-4, 1e-4);
        let exact = card(p).rolling_time_s();
        let coarse = reduced_modulus_logarithmic_loss(p, 10.0);
        let fine = reduced_modulus_logarithmic_loss(p, 5.0);
        let extrapolated = (4.0 * fine - coarse) / 3.0;
        close(exact, extrapolated, 1e-8);
        close(
            card(parameters(nu, 2e-4, 2e-4)).rolling_time_s(),
            2e-4,
            5e-15,
        );
        assert!(exact >= 1e-4 && exact <= 4e-4);
    }
}

#[test]
fn independent_viscous_channels_are_positive_linear_and_monotone() {
    for nu in [-0.3, 0.0, 0.3, 0.499] {
        let both = card(parameters(nu, 2e-4, 3e-4)).rolling_time_s();
        let shear = card(parameters(nu, 2e-4, 0.0)).rolling_time_s();
        let bulk = card(parameters(nu, 0.0, 3e-4)).rolling_time_s();
        assert!(shear > 0.0 && bulk > 0.0);
        close(both, shear + bulk, 5e-15);
        assert!(card(parameters(nu, 3e-4, 3e-4)).rolling_time_s() > both);
        assert!(card(parameters(nu, 2e-4, 4e-4)).rolling_time_s() > both);
    }
}

// Independently integrate the first moment of the advected Hertz stress
// correction -tau*V*d_x(p). r=a*sin(theta) removes the edge singularity; the
// full two-dimensional pressure gradient and area Jacobian remain explicit.
fn pressure_moment(load: f64, area: f64, tau: f64, speed: f64, n: usize) -> f64 {
    let a = (area / core::f64::consts::PI).sqrt();
    let p0 = 1.5 * load / area;
    let dtheta = core::f64::consts::FRAC_PI_2 / n as f64;
    let dphi = core::f64::consts::TAU / 64.0;
    let mut moment = 0.0;
    for i in 0..n {
        let theta = (i as f64 + 0.5) * dtheta;
        let r = a * theta.sin();
        let dr = a * theta.cos() * dtheta;
        for j in 0..64 {
            let phi = (j as f64 + 0.5) * dphi;
            let x = r * phi.cos();
            let pressure_gradient = -p0 * x / (a * a * theta.cos());
            let correction = -tau * speed * pressure_gradient;
            moment += x * correction * r * dr * dphi;
        }
    }
    moment
}

#[test]
fn actual_torque_and_heat_match_independent_pressure_moment_with_refinement() {
    let p = parameters(0.3, 2e-4, 2e-4);
    let law = RollingLossLaw::ViscoelasticSphere(card(p));
    let patch = patch(p, 2e-4);
    let kinematics = motion(p, 0.5, 0.0, 0.25);
    let step = law
        .advance(
            &patch,
            &interface(),
            kinematics,
            &ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let moment = pressure_moment(
            patch.normal_load_n(),
            patch.contact_area_m2(),
            2e-4,
            kinematics.contour_speed_mps(),
            n,
        );
        errors.push((moment + step.wrench.rolling_moment_nm).abs());
    }
    assert!(
        errors[1] < 0.1 * errors[0] && errors[2] < 0.1 * errors[1],
        "{errors:?}"
    );
    let moment = pressure_moment(
        patch.normal_load_n(),
        patch.contact_area_m2(),
        2e-4,
        kinematics.contour_speed_mps(),
        256,
    );
    close(-step.wrench.rolling_moment_nm, moment, 2e-10);
    close(step.dissipation.total_heat_j, moment * 0.5 * 0.25, 2e-10);
    assert_eq!(step.wrench.contour_force_n, 0.0);
    assert_eq!(step.wrench.spin_moment_nm, 0.0);
    assert_eq!(step.dissipation.contour_deformation_heat_j, 0.0);
    assert_eq!(
        step.dissipation.rolling_hysteresis_heat_j,
        step.dissipation.total_heat_j
    );
    assert_eq!(step.stored_energy_j, 0.0);
    assert_eq!(step.storage_change_j, 0.0);
    assert_eq!(
        step.applicability,
        RollingLossApplicabilityKind::FirstOrderViscoelasticSphere
    );
}

#[test]
fn zero_loss_release_rest_and_pure_spin_generate_no_rolling_heat() {
    let p = parameters(0.0, 2e-4, 2e-4);
    let c = card(p);
    let law = RollingLossLaw::ViscoelasticSphere(c.clone());
    for (normal, kinematics) in [
        (patch(p, 0.0), motion(p, 1e6, 0.0, 0.5)),
        (patch(p, 2e-4), motion(p, 0.0, 0.0, 0.5)),
        (patch(p, 2e-4), motion(p, 0.0, 500.0, 0.5)),
    ] {
        let step = law
            .advance(
                &normal,
                &interface(),
                kinematics,
                &ownership(),
                &RollingLossState::zero(),
            )
            .unwrap();
        assert_eq!(step.wrench.rolling_moment_nm, 0.0);
        assert_eq!(step.wrench.spin_moment_nm, 0.0);
        assert_eq!(step.dissipation.total_heat_j, 0.0);
        assert_eq!(step.state, RollingLossStateKind::Quiescent);
    }
    let elastic = card(parameters(0.0, 0.0, 0.0));
    assert_eq!(elastic.rolling_time_s(), 0.0);
    assert_eq!(
        elastic
            .response(&patch(p, 2e-4), &interface(), motion(p, 0.5, 0.0, 1.0))
            .unwrap()
            .dissipated_power_w,
        0.0
    );
}

#[test]
fn reversal_speed_load_and_geometry_scalings_follow_the_same_physics() {
    let p = parameters(0.3, 2e-4, 3e-4);
    let c = card(p);
    let normal = patch(p, 2e-4);
    let base = c
        .response(&normal, &interface(), motion(p, 0.5, 0.0, 0.5))
        .unwrap();
    let reverse = c
        .response(&normal, &interface(), motion(p, -0.5, -30.0, 0.5))
        .unwrap();
    assert_eq!(base.rolling_moment_nm, -reverse.rolling_moment_nm);
    assert_eq!(base.dissipated_power_w, reverse.dissipated_power_w);
    let faster = c
        .response(&normal, &interface(), motion(p, 1.0, 80.0, 0.5))
        .unwrap();
    close(
        faster.rolling_moment_nm,
        2.0 * base.rolling_moment_nm,
        5e-15,
    );
    close(
        faster.dissipated_power_w,
        4.0 * base.dissipated_power_w,
        5e-15,
    );
    let loaded = c
        .response(
            &patch(p, 2e-4 * 2.0_f64.powf(2.0 / 3.0)),
            &interface(),
            motion(p, 0.5, 0.0, 0.5),
        )
        .unwrap();
    close(
        loaded.rolling_moment_nm,
        2.0 * base.rolling_moment_nm,
        2e-14,
    );
    let scale = 3.5;
    let scaled_p = ViscoelasticSphereParameters {
        radius_m: p.radius_m * scale,
        ..p
    };
    let scaled = card(scaled_p)
        .response(
            &patch(scaled_p, 2e-4 * scale),
            &interface(),
            motion(scaled_p, 0.5, 0.0, 0.5),
        )
        .unwrap();
    close(
        scaled.rolling_moment_nm,
        scale * scale * scale * base.rolling_moment_nm,
        2e-14,
    );
    close(scaled.relaxation_ratio, base.relaxation_ratio, 2e-14);
}

#[test]
fn original_hertz_geometry_and_no_slip_mapping_are_required() {
    let p = parameters(0.3, 2e-4, 2e-4);
    let c = card(p);
    let normal = patch(p, 2e-4);
    let kinematics = motion(p, 0.5, 0.0, 0.5);
    for (load, area, curvature) in [
        (
            normal.normal_load_n() * 1.01,
            normal.contact_area_m2(),
            normal.curvature(),
        ),
        (
            normal.normal_load_n(),
            normal.contact_area_m2() * 1.01,
            normal.curvature(),
        ),
        (
            normal.normal_load_n(),
            normal.contact_area_m2(),
            PatchCurvature::Principal {
                first_per_m: 50.0,
                second_per_m: 25.0,
            },
        ),
        (
            normal.normal_load_n(),
            normal.contact_area_m2(),
            PatchCurvature::EquivalentRadiusApproximation {
                radius_m: p.radius_m,
                authority: InputAuthority::Estimated,
            },
        ),
    ] {
        let wrong = RollingPatchReceipt::new(
            "sphere-patch",
            "shared-Hertz",
            "normal-source-1",
            InputAuthority::SyntheticFixture,
            load,
            area,
            curvature,
        )
        .unwrap();
        assert!(c.response(&wrong, &interface(), kinematics).is_err());
    }
    let slip = RollingKinematics::new(0.011, 0.5, 0.0, 300.0, 10.0, 0.5).unwrap();
    assert!(c.response(&normal, &interface(), slip).is_err());
    let changed = card(ViscoelasticSphereParameters {
        young_modulus_pa: 2.1e7,
        ..p
    });
    assert!(changed.response(&normal, &interface(), kinematics).is_err());
    assert_eq!(
        c.response(&normal, &interface(), kinematics).unwrap(),
        c.response(&normal, &interface(), kinematics).unwrap()
    );
}

#[test]
fn approximation_and_source_domain_limits_refuse_before_any_heat_is_accepted() {
    let p = parameters(0.3, 2e-4, 2e-4);
    let c = card(p);
    let normal = patch(p, 2e-4);
    let kinematics = motion(p, 0.5, 0.0, 0.5);
    let admitted = c.response(&normal, &interface(), kinematics).unwrap();
    let tight = QuasistaticSphereLimits {
        maximum_indentation_radius_ratio: admitted.indentation_radius_ratio,
        maximum_speed_wave_ratio: admitted.speed_wave_ratio,
        maximum_relaxation_ratio: admitted.relaxation_ratio,
    };
    let bounded = ViscoelasticSphereRollingCard::new(
        VISCOELASTIC_SPHERE_LAW_ID,
        "fixture/tight",
        InputAuthority::SyntheticFixture,
        p,
        tight,
        domain(),
    )
    .unwrap();
    assert!(bounded.response(&normal, &interface(), kinematics).is_ok());
    for tighter in [
        QuasistaticSphereLimits {
            maximum_indentation_radius_ratio: 0.9 * tight.maximum_indentation_radius_ratio,
            ..tight
        },
        QuasistaticSphereLimits {
            maximum_speed_wave_ratio: 0.9 * tight.maximum_speed_wave_ratio,
            ..tight
        },
        QuasistaticSphereLimits {
            maximum_relaxation_ratio: 0.9 * tight.maximum_relaxation_ratio,
            ..tight
        },
    ] {
        let refused = ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            "fixture/tight",
            InputAuthority::SyntheticFixture,
            p,
            tighter,
            domain(),
        )
        .unwrap();
        assert!(refused.response(&normal, &interface(), kinematics).is_err());
    }
    let law = RollingLossLaw::ViscoelasticSphere(c);
    let accepted = law
        .advance(
            &normal,
            &interface(),
            kinematics,
            &ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    let saved = accepted.next_state.clone();
    for (temperature, frequency) in [(279.0, 10.0), (321.0, 10.0), (300.0, 101.0)] {
        let outside = RollingKinematics::new(0.01, 0.5, 0.0, temperature, frequency, 0.5).unwrap();
        assert!(
            law.advance(&normal, &interface(), outside, &ownership(), &saved)
                .is_err()
        );
    }
    assert_eq!(accepted.next_state, saved);
    let whole = law
        .advance(
            &normal,
            &interface(),
            motion(p, 0.5, 0.0, 1.0),
            &ownership(),
            &saved,
        )
        .unwrap();
    let half = law
        .advance(&normal, &interface(), kinematics, &ownership(), &saved)
        .unwrap();
    let halves = law
        .advance(
            &normal,
            &interface(),
            kinematics,
            &ownership(),
            &half.next_state,
        )
        .unwrap();
    close(
        whole.next_state.cumulative_heat_j(),
        halves.next_state.cumulative_heat_j(),
        2e-15,
    );
}

#[test]
fn checkpoints_retain_material_body_role_and_authority_and_bad_constants_refuse() {
    let p = parameters(0.3, 2e-4, 2e-4);
    let c = card(p);
    let law = RollingLossLaw::ViscoelasticSphere(c.clone());
    let normal = patch(p, 2e-4);
    let accepted = law
        .advance(
            &normal,
            &interface(),
            motion(p, 0.5, 0.0, 0.5),
            &ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(
        law.restore_checkpoint(&normal, &interface(), &ownership(), &accepted.checkpoint)
            .unwrap(),
        accepted.next_state
    );
    for changed in [
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            "another-source",
            InputAuthority::SyntheticFixture,
            p,
            limits(),
            domain(),
        )
        .unwrap(),
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            "fixture/viscosities-1",
            InputAuthority::Estimated,
            p,
            limits(),
            domain(),
        )
        .unwrap(),
        card(ViscoelasticSphereParameters {
            sphere_surface: DeformingSphereSurface::SurfaceB,
            ..p
        }),
        card(parameters(0.3, 3e-4, 2e-4)),
    ] {
        assert!(
            RollingLossLaw::ViscoelasticSphere(changed)
                .restore_checkpoint(&normal, &interface(), &ownership(), &accepted.checkpoint)
                .is_err()
        );
    }
    for malformed in [
        ViscoelasticSphereParameters {
            young_modulus_pa: 0.0,
            ..p
        },
        ViscoelasticSphereParameters {
            radius_m: -1.0,
            ..p
        },
        ViscoelasticSphereParameters {
            density_kg_m3: f64::INFINITY,
            ..p
        },
        ViscoelasticSphereParameters {
            poisson_ratio: 0.5,
            ..p
        },
        ViscoelasticSphereParameters {
            poisson_ratio: -1.0,
            ..p
        },
        ViscoelasticSphereParameters {
            shear_viscosity_pa_s: -1.0,
            ..p
        },
        ViscoelasticSphereParameters {
            bulk_viscosity_pa_s: f64::NAN,
            ..p
        },
    ] {
        assert!(
            ViscoelasticSphereRollingCard::new(
                VISCOELASTIC_SPHERE_LAW_ID,
                "fixture/bad",
                InputAuthority::SyntheticFixture,
                malformed,
                limits(),
                domain()
            )
            .is_err()
        );
    }
}
