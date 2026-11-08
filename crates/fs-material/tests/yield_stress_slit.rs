use fs_matdb::Provenance;
use fs_material::fluid_transport::rheology::yield_stress::{
    LaminarHerschelBulkleySlit, YieldStressSlitControls,
};
use fs_material::slurry::rheology::{
    HerschelBulkleyCard, HerschelBulkleyParameters, ShearResponse, SlurryConditions, SlurryDomain,
    SlurryError, SlurryTestContext,
};
use fs_qty::{Dimensionless, Frequency, Length, Pressure, Stress, Temperature, Velocity};

fn conditions() -> SlurryConditions {
    SlurryConditions {
        temperature: Temperature::new(298.15),
        solids_mass_fraction: Dimensionless::new(0.4),
    }
}

fn card_with_minimum(
    yield_stress: f64,
    reference_stress: f64,
    n: f64,
    minimum: f64,
) -> HerschelBulkleyCard {
    HerschelBulkleyCard::new(
        HerschelBulkleyParameters {
            yield_stress: Stress::new(yield_stress),
            reference_excess_stress: Stress::new(reference_stress),
            reference_shear_rate: Frequency::new(10.0),
            flow_index: Dimensionless::new(n),
        },
        SlurryDomain::new(
            [Frequency::new(minimum), Frequency::new(1000.0)],
            [Temperature::new(293.15), Temperature::new(303.15)],
            [Dimensionless::new(0.35), Dimensionless::new(0.45)],
        )
        .unwrap(),
        SlurryTestContext {
            specimen: "authored homogeneous suspension; not measured material data".into(),
            method: "independent exact slit and manufactured-profile regression".into(),
            provenance: Provenance {
                source: "authored arithmetic fixture using the original Herschel-Bulkley law"
                    .into(),
                license: "project test fixture".into(),
                artifact: None,
            },
        },
    )
    .unwrap()
}

fn card(y: f64, s: f64, n: f64) -> HerschelBulkleyCard {
    card_with_minimum(y, s, n, 0.0)
}

fn slit(card: &HerschelBulkleyCard, h: f64) -> LaminarHerschelBulkleySlit<'_> {
    LaminarHerschelBulkleySlit::new(
        card,
        conditions(),
        Length::new(0.03),
        Length::new(h),
        Length::new(0.5),
        Frequency::new(1000.0),
    )
    .unwrap()
}

fn controls() -> YieldStressSlitControls {
    YieldStressSlitControls {
        relative_wall_velocity_tolerance: 2.0e-13,
        maximum_iterations: 96,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-20),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn newtonian_combined_pressure_motion_profile_heat_and_tangent() {
    let source = card(0.0, 4.0, 1.0);
    let h: f64 = 0.002;
    let (w, length, mu) = (0.03, 0.5, 0.4);
    let owner = slit(&source, h);
    for dp in [-1000.0, 0.0, 1000.0] {
        let walls = [-0.01, 0.03];
        let relative = walls[1] - walls[0];
        let mean = 0.5 * (walls[0] + walls[1]);
        let profile = owner
            .profile_at_pressure(Pressure::new(dp), walls.map(Velocity::new), controls())
            .unwrap();
        let flow = profile.flow();
        let conductance = w * h.powi(3) / (12.0 * mu * length);
        close(
            flow.volume_flow_m3_s,
            w * h * mean + conductance * dp,
            3.0e-12,
        );
        close(
            flow.dissipation_w,
            conductance * dp * dp + w * length * mu * relative * relative / h,
            3.0e-12,
        );
        close(flow.midplane_stress_pa, mu * relative / h, 3.0e-12);
        let tangent = flow.tangent.unwrap();
        close(tangent.pressure_conductance_m3_s_pa, conductance, 3.0e-12);
        close(
            tangent.relative_wall_resistance_n_s_m,
            w * length * mu / h,
            3.0e-12,
        );
        assert!(tangent.flow_relative_wall_velocity_derivative_m2.abs() < 1.0e-16);
        for x in [0.0, 0.13, 0.5, 0.87, 1.0] {
            let sample = profile.sample(Length::new(h * x)).unwrap();
            let velocity =
                walls[0] + relative * x + dp * h * h * x * (1.0 - x) / (2.0 * mu * length);
            let rate = relative / h + dp * h * (1.0 - 2.0 * x) / (2.0 * mu * length);
            close(sample.velocity_m_s, velocity, 5.0e-12);
            close(sample.shear_rate_per_s, rate, 5.0e-12);
            close(sample.shear_stress_pa, mu * rate, 5.0e-12);
            close(sample.dissipation_w_m3, mu * rate * rate, 1.0e-11);
        }
    }
    let rest = owner
        .profile_at_pressure(Pressure::new(0.0), [Velocity::new(0.0); 2], controls())
        .unwrap();
    assert_eq!(rest.flow().dissipation_w, 0.0);
    close(
        rest.flow().tangent.unwrap().pressure_conductance_m3_s_pa,
        w * h.powi(3) / (12.0 * mu * length),
        3.0e-14,
    );
}

#[test]
fn exact_bingham_plug_flow_mobility_and_translating_static_stress_intervals() {
    let source = card(7.0, 12.0, 1.0);
    let (h, w, length, mu): (f64, f64, f64, f64) = (0.004, 0.03, 0.5, 1.2);
    let owner = slit(&source, h);
    let p = 20.0;
    let dp = 2.0 * length * p / h;
    let profile = owner
        .profile_at_pressure(Pressure::new(dp), [Velocity::new(0.0); 2], controls())
        .unwrap();
    let flow = profile.flow();
    let a: f64 = 7.0 / p;
    let conductance0 = w * h.powi(3) / (12.0 * mu * length);
    let expected_q = conductance0 * dp * (1.0 - 1.5 * a + 0.5 * a.powi(3));
    close(flow.volume_flow_m3_s, expected_q, 4.0e-14);
    close(flow.dissipation_w, dp * expected_q, 4.0e-14);
    let plug = flow.unyielded_interval_m.unwrap();
    close(plug[0], 0.5 * h * (1.0 - a), 2.0e-14);
    close(plug[1], 0.5 * h * (1.0 + a), 2.0e-14);
    let center = profile.sample(Length::new(0.5 * h)).unwrap();
    close(
        center.velocity_m_s,
        h * (p - 7.0).powi(2) / (4.0 * mu * p),
        3.0e-14,
    );
    assert!(center.unyielded);
    assert_eq!(center.shear_rate_per_s, 0.0);
    assert_eq!(center.dissipation_w_m3, 0.0);
    let tangent = flow.tangent.unwrap();
    close(
        tangent.pressure_conductance_m3_s_pa,
        conductance0 * (1.0 - a.powi(3)),
        4.0e-14,
    );
    close(
        tangent.relative_wall_resistance_n_s_m,
        w * length * mu / (h * (1.0 - a)),
        4.0e-14,
    );

    for p in [-7.0_f64, -3.0, 0.0, 3.0, 7.0] {
        let dp = 2.0 * length * p / h;
        let moving = owner
            .profile_at_pressure(Pressure::new(dp), [Velocity::new(0.2); 2], controls())
            .unwrap();
        let flow = moving.flow();
        close(flow.volume_flow_m3_s, w * h * 0.2, 2.0e-14);
        assert_eq!(flow.dissipation_w, 0.0);
        assert_eq!(flow.relative_volume_flow_m3_s, 0.0);
        assert_eq!(flow.unyielded_interval_m, Some([0.0, h]));
        assert_eq!(
            flow.static_midplane_stress_interval_pa,
            Some([-(7.0 - p.abs()), 7.0 - p.abs()])
        );
        assert!(flow.tangent.is_none());
        close(flow.pressure_work_w, dp * w * h * 0.2, 2.0e-14);
        close(flow.wall_work_w, -dp * w * h * 0.2, 2.0e-14);
        assert_eq!(flow.work_balance_residual_w, 0.0);
        assert_eq!(
            flow.wall_traction_bounds_pa[0],
            [-7.0 + p.abs() - p, 7.0 - p.abs() - p]
        );
        for x in [0.0, 0.25, 0.7, 1.0] {
            let sample = moving.sample(Length::new(h * x)).unwrap();
            close(sample.velocity_m_s, 0.2, 1.0e-14);
            assert_eq!(sample.shear_rate_per_s, 0.0);
            assert_eq!(sample.dissipation_w_m3, 0.0);
        }
    }
}

// An independent primitive of the signed inverse stress law. The native owner
// computes flux/heat through positive differential moments, not this formula.
fn primitive(stress: f64, y: f64, s: f64, m: f64) -> f64 {
    let excess = (stress.abs() - y).max(0.0);
    10.0 * excess * (excess / s).powf(m) / (m + 1.0)
}

#[test]
fn noninteger_mixed_flows_match_manufactured_profiles_heat_and_coupled_tangents() {
    let (h, w, length) = (0.004, 0.03, 0.5);
    let (y, s, p, c) = (4.0, 8.0, 15.0, 3.0);
    let dp = 2.0 * length * p / h;
    for n in [0.6, 1.5] {
        let m = 1.0 / n;
        let source = card(y, s, n);
        let owner = slit(&source, h);
        let relative = h * (primitive(c + p, y, s, m) - primitive(c - p, y, s, m)) / (2.0 * p);
        let walls = [-0.02, -0.02 + relative];
        let expected_velocity = |x: f64| {
            walls[0]
                + h * (primitive(c + p, y, s, m) - primitive(c + p - 2.0 * p * x, y, s, m))
                    / (2.0 * p)
        };
        let profile = owner
            .profile_at_pressure(Pressure::new(dp), walls.map(Velocity::new), controls())
            .unwrap();
        let flow = profile.flow();
        close(flow.midplane_stress_pa, c, 3.0e-12);
        for x in [0.0, 0.1, 0.35, 0.5, 0.6, 0.82, 1.0] {
            let sample = profile.sample(Length::new(h * x)).unwrap();
            close(sample.velocity_m_s, expected_velocity(x), 3.0e-11);
            assert!((sample.shear_stress_pa - (c + p * (1.0 - 2.0 * x))).abs() < 3.0e-11 * p);
            let tau = c + p * (1.0 - 2.0 * x);
            let gamma = tau.signum() * 10.0 * ((tau.abs() - y).max(0.0) / s).powf(m);
            close(sample.shear_rate_per_s, gamma, 3.0e-11);
        }
        let intervals = 16384;
        let mut velocity_integral = expected_velocity(0.0) + expected_velocity(1.0);
        for i in 1..intervals {
            velocity_integral += (if i % 2 == 0 { 2.0 } else { 4.0 })
                * expected_velocity(f64::from(i) / f64::from(intervals));
        }
        velocity_integral /= 3.0 * f64::from(intervals);
        close(flow.volume_flow_m3_s, w * h * velocity_integral, 2.0e-8);
        let heat: f64 = [c + p, c - p]
            .into_iter()
            .map(|stress| {
                let excess = stress.abs() - y;
                let ell = excess / (2.0 * p);
                w * length
                    * h
                    * ell
                    * 10.0
                    * (excess / s).powf(m)
                    * (y / (m + 1.0) + excess / (m + 2.0))
            })
            .sum();
        close(flow.dissipation_w, heat, 4.0e-12);
        // Use an absolute work scale here: the exact balance residual is near zero.
        assert!(
            (flow.pressure_work_w + flow.wall_work_w
                - flow.dissipation_w
                - flow.work_balance_residual_w)
                .abs()
                < 1.0e-13 * heat
        );
        let tangent = flow.tangent.unwrap();
        let pressure_step = dp * 1.0e-5;
        let at_pressure = |pressure: f64| {
            owner
                .profile_at_pressure(
                    Pressure::new(pressure),
                    walls.map(Velocity::new),
                    controls(),
                )
                .unwrap()
        };
        let numerical_pressure = (at_pressure(dp + pressure_step).flow().volume_flow_m3_s
            - at_pressure(dp - pressure_step).flow().volume_flow_m3_s)
            / (2.0 * pressure_step);
        close(
            tangent.pressure_conductance_m3_s_pa,
            numerical_pressure,
            3.0e-7,
        );
        let step = relative.abs() * 1.0e-5;
        let mean = 0.5 * (walls[0] + walls[1]);
        let at_relative = |relative: f64| {
            owner
                .profile_at_pressure(
                    Pressure::new(dp),
                    [
                        Velocity::new(mean - 0.5 * relative),
                        Velocity::new(mean + 0.5 * relative),
                    ],
                    controls(),
                )
                .unwrap()
        };
        let plus = at_relative(relative + step);
        let minus = at_relative(relative - step);
        close(
            tangent.flow_relative_wall_velocity_derivative_m2,
            (plus.flow().volume_flow_m3_s - minus.flow().volume_flow_m3_s) / (2.0 * step),
            3.0e-7,
        );
        close(
            tangent.relative_wall_resistance_n_s_m,
            w * length * (plus.flow().midplane_stress_pa - minus.flow().midplane_stress_pa)
                / (2.0 * step),
            3.0e-7,
        );
    }
}

#[test]
fn reversal_translation_ceiling_support_and_cancellation_are_explicit() {
    let source = card(3.0, 4.0, 0.6);
    let h = 0.004;
    let owner = slit(&source, h);
    let walls = [-0.01, 0.025];
    let original = owner
        .profile_at_pressure(Pressure::new(3000.0), walls.map(Velocity::new), controls())
        .unwrap();
    let reverse = owner
        .profile_at_pressure(
            Pressure::new(-3000.0),
            walls.map(|u| Velocity::new(-u)),
            controls(),
        )
        .unwrap();
    close(
        reverse.flow().volume_flow_m3_s,
        -original.flow().volume_flow_m3_s,
        3.0e-12,
    );
    close(
        reverse.flow().dissipation_w,
        original.flow().dissipation_w,
        3.0e-12,
    );
    let translated = owner
        .profile_at_pressure(
            Pressure::new(3000.0),
            walls.map(|u| Velocity::new(u + 0.1)),
            controls(),
        )
        .unwrap();
    close(
        translated.flow().volume_flow_m3_s,
        original.flow().volume_flow_m3_s + 0.03 * h * 0.1,
        3.0e-12,
    );
    close(
        translated.flow().dissipation_w,
        original.flow().dissipation_w,
        3.0e-12,
    );
    // A large common translation must not silently erase the small accepted
    // wall difference when rebuilding physical wall velocities.
    let large = 2.0_f64.powi(40);
    let adjacent = f64::from_bits(large.to_bits() + 1);
    match owner.profile_at_pressure(
        Pressure::new(3000.0),
        [Velocity::new(large), Velocity::new(adjacent)],
        controls(),
    ) {
        Ok(profile) => {
            let flow = profile.flow();
            let accepted = flow.accepted_wall_velocities_m_s;
            let represented = accepted[1] - accepted[0];
            let implied = adjacent - large + flow.relative_wall_velocity_residual_m_s;
            assert!(represented > 0.0);
            assert!(
                (represented - implied).abs()
                    <= controls().relative_wall_velocity_tolerance * (adjacent - large)
            );
        }
        Err(SlurryError::Invalid(
            "common wall motion hides the accepted relative slit velocity",
        )) => {}
        Err(error) => panic!("unexpected large-translation refusal: {error}"),
    }
    let before = source.clone();
    let mut checks = 0;
    let cancelled = owner.profile_at_pressure_with_check(
        Pressure::new(3000.0),
        walls.map(Velocity::new),
        controls(),
        || {
            checks += 1;
            if checks == 3 {
                Err(SlurryError::Invalid("authored cancellation"))
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(SlurryError::Invalid("authored cancellation"))
    ));
    assert_eq!(source, before);
    let retry = owner
        .profile_at_pressure(Pressure::new(3000.0), walls.map(Velocity::new), controls())
        .unwrap();
    assert_eq!(retry.flow(), original.flow());
    let bounded = owner.profile_at_pressure(
        Pressure::new(3000.0),
        walls.map(Velocity::new),
        YieldStressSlitControls {
            maximum_iterations: 1,
            ..controls()
        },
    );
    assert!(matches!(
        bounded,
        Err(SlurryError::NoConvergence { iterations: 1 })
    ));
    assert!(
        owner
            .profile_at_pressure(Pressure::new(1.0e9), walls.map(Velocity::new), controls())
            .is_err()
    );

    let h = 1.0 / 1024.0;
    let ceiling = 37.5;
    let limited = LaminarHerschelBulkleySlit::new(
        &source,
        conditions(),
        Length::new(0.03),
        Length::new(h),
        Length::new(0.5),
        Frequency::new(ceiling),
    )
    .unwrap();
    let endpoint = limited
        .profile_at_pressure(
            Pressure::new(0.0),
            [Velocity::new(0.0), Velocity::new(ceiling * h)],
            controls(),
        )
        .unwrap();
    assert_eq!(
        endpoint.sample(Length::new(0.0)).unwrap().shear_rate_per_s,
        ceiling
    );
    let ShearResponse::Yielded { stress, .. } = source
        .stress_at_shear_rate(Frequency::new(ceiling), conditions())
        .unwrap()
    else {
        panic!("positive endpoint must yield")
    };
    assert_eq!(endpoint.flow().midplane_stress_pa, stress.value());
    let above = f64::from_bits(ceiling.to_bits() + 1);
    assert!(
        limited
            .profile_at_pressure(
                Pressure::new(0.0),
                [Velocity::new(0.0), Velocity::new(above * h)],
                controls()
            )
            .is_err()
    );
    let positive_only = card_with_minimum(3.0, 4.0, 0.6, 0.1);
    assert!(
        LaminarHerschelBulkleySlit::new(
            &positive_only,
            conditions(),
            Length::new(0.03),
            Length::new(h),
            Length::new(0.5),
            Frequency::new(ceiling)
        )
        .is_err()
    );
    let invalid_conditions = SlurryConditions {
        temperature: Temperature::new(310.0),
        ..conditions()
    };
    assert!(
        LaminarHerschelBulkleySlit::new(
            &source,
            invalid_conditions,
            Length::new(0.03),
            Length::new(h),
            Length::new(0.5),
            Frequency::new(ceiling)
        )
        .is_err()
    );
}
