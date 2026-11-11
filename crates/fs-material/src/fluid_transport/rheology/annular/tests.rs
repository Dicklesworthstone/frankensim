//! Analytical and independent synthetic-law controls, not fluid calibration.
use super::*;
use crate::fluid_transport::rheology::slit::{CarreauSlitControls, LaminarCarreauSlit};
use crate::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law(constant: bool) -> ThermalCarreauYasuda {
    let curve = if constant {
        [3.0, 3.0, 0.6, 0.5, 1.0]
    } else {
        [12.0, 0.0, 0.6, 0.5, 1.0]
    };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(curve[0], curve[1], curve[2], curve[3], curve[4], [0.0, 1e9]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn tight() -> CarreauAnnularControls {
    CarreauAnnularControls {
        local: CarreauPipeControls {
            relative_constitutive_tolerance: 1e-13,
            relative_quadrature_tolerance: 2e-11,
            ..CarreauPipeControls::default()
        },
        relative_wall_velocity_tolerance: 1e-11,
        ..CarreauAnnularControls::default()
    }
}

fn annulus() -> LaminarCarreauAnnulus {
    LaminarCarreauAnnulus::new(law(false), 300.0, 0.01, 0.04, 0.5, 1e6).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let scale = actual.abs().max(expected.abs()).max(1e-20);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_newtonian_annular_profiles_forces_flux_heat_and_conductance() {
    for a in [0.01_f64, 0.0001] {
        let (b, length, mu) = (0.04_f64, 0.5, 3.0);
        let owner = LaminarCarreauAnnulus::new(law(true), 300.0, a, b, length, 1e6).unwrap();
        let log = (b / a).ln();
        let area = core::f64::consts::PI * (b * b - a * a);
        let conductance = core::f64::consts::PI / (8.0 * mu * length)
            * (b.powi(4) - a.powi(4) - (b * b - a * a).powi(2) / log);
        for pressure in [-2000.0, 0.0, 2000.0] {
            for walls in [[0.0, 0.0], [-0.2, 0.8], [1.0, 1.0]] {
                let profile = owner
                    .profile_at_pressure_and_walls(pressure, walls, tight())
                    .unwrap();
                let f = profile.flow();
                let relative = walls[1] - walls[0];
                let c = mu * relative / log + pressure * (b * b - a * a) / (4.0 * length * log);
                let expected_q = area * walls[0]
                    + core::f64::consts::PI * relative * (b * b - (b * b - a * a) / (2.0 * log))
                    + conductance * pressure;
                close(f.volumetric_flow_m3_s, expected_q, 5e-12);
                close(f.pressure_conductance_m3_s_pa, conductance, 5e-12);
                close(
                    f.dissipated_power_w,
                    conductance * pressure * pressure
                        + core::f64::consts::TAU * mu * length * relative * relative / log,
                    5e-12,
                );
                close(
                    f.wall_tractions_pa[0],
                    -(c / a - pressure * a / (2.0 * length)),
                    5e-12,
                );
                close(
                    f.wall_tractions_pa[1],
                    c / b - pressure * b / (2.0 * length),
                    5e-12,
                );
                let balance = pressure * area + f.wall_forces_n[0] + f.wall_forces_n[1];
                let force_scale =
                    (pressure * area).abs() + f.wall_forces_n[0].abs() + f.wall_forces_n[1].abs();
                assert!(balance.abs() <= 1e-12 * force_scale.max(1.0));
                assert!(f.work_residual_w.abs() <= 1e-11 * f.dissipated_power_w.max(1.0));
                for fraction in [0.0, 0.13, 0.5, 0.87, 1.0] {
                    let radius = a + fraction * (b - a);
                    let sample = profile.sample_gap_fraction(fraction).unwrap();
                    let expected = walls[0] + c / mu * (radius / a).ln()
                        - pressure * (radius * radius - a * a) / (4.0 * mu * length);
                    assert!(
                        (sample.axial_velocity_m_s - expected).abs()
                            <= 2e-9 * expected.abs().max(1.0)
                    );
                    close(
                        sample.radial_velocity_gradient_per_s,
                        (c / radius - pressure * radius / (2.0 * length)) / mu,
                        2e-10,
                    );
                }
            }
        }
    }
}

// Independent exact inverse for n=1/2, a=1, mu_inf=0. The production code
// instead uses the original source's bounded Newton constitutive inverse.
fn inverse(stress: f64) -> f64 {
    let z = stress.abs() * 0.6 / 12.0;
    ((z * z + z * (z * z + 4.0).sqrt()) / 1.2).copysign(stress)
}

fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let n = 8192;
    let mut sum = 0.0;
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * f(a + (b - a) * f64::from(i) / f64::from(n));
    }
    (b - a) * sum / (3.0 * f64::from(n))
}

fn oracle(pressure: f64, walls: [f64; 2]) -> (f64, f64, f64) {
    let (a, b, length) = (0.01, 0.04, 0.5);
    let (mut low, mut high) = (-100.0, 100.0);
    for _ in 0..80 {
        let c = 0.5 * (low + high);
        let relative = simpson(|r| inverse(c / r - pressure * r / (2.0 * length)), a, b);
        if relative < walls[1] - walls[0] {
            low = c;
        } else {
            high = c;
        }
    }
    let c = 0.5 * (low + high);
    let flow = core::f64::consts::PI * (b * b - a * a) * walls[0]
        + core::f64::consts::PI
            * simpson(
                |r| (b * b - r * r) * inverse(c / r - pressure * r / (2.0 * length)),
                a,
                b,
            );
    let heat = core::f64::consts::TAU
        * length
        * simpson(
            |r| {
                let stress = c / r - pressure * r / (2.0 * length);
                r * stress * inverse(stress)
            },
            a,
            b,
        );
    (c, flow, heat)
}

#[test]
fn g1_nonlinear_curved_profile_matches_independent_cylindrical_momentum() {
    let owner = annulus();
    for (pressure, walls) in [
        (1000.0, [-0.2, 0.6]),
        (-1000.0, [0.2, -0.6]),
        (0.0, [0.0, 1.0]),
        (3000.0, [0.0, 0.0]),
    ] {
        let profile = owner
            .profile_at_pressure_and_walls(pressure, walls, tight())
            .unwrap();
        let (c, flow, heat) = oracle(pressure, walls);
        close(profile.flow().volumetric_flow_m3_s, flow, 2e-8);
        close(profile.flow().dissipated_power_w, heat, 2e-8);
        assert!(profile.flow().work_residual_w.abs() <= 2e-8 * heat);
        for fraction in [0.1, 0.4, 0.8] {
            let radius = 0.01 + 0.03 * fraction;
            let sample = profile.sample_gap_fraction(fraction).unwrap();
            let stress = c / radius - pressure * radius;
            let velocity = walls[0] + simpson(|r| inverse(c / r - pressure * r), 0.01, radius);
            assert!((sample.axial_velocity_m_s - velocity).abs() <= 2e-8 * velocity.abs().max(1.0));
            close(sample.shear_stress_pa, stress, 2e-8);
            close(sample.dissipation_w_m3, stress * inverse(stress), 4e-8);
        }
    }
}

#[test]
fn g2_nonlinear_pressure_mobility_and_reciprocal_wall_coupling() {
    let owner = annulus();
    let walls = [-0.2, 0.6];
    let at = |p, walls| {
        owner
            .profile_at_pressure_and_walls(p, walls, tight())
            .unwrap()
    };
    let base = at(1000.0, walls);
    let dp = 0.02;
    let plus = at(1000.0 + dp, walls);
    let minus = at(1000.0 - dp, walls);
    let mobility =
        (plus.flow().volumetric_flow_m3_s - minus.flow().volumetric_flow_m3_s) / (2.0 * dp);
    close(base.flow().pressure_conductance_m3_s_pa, mobility, 3e-5);
    assert!(mobility > 0.0);
    let du = 1e-4;
    let upper = at(1000.0, [walls[0] - 0.5 * du, walls[1] + 0.5 * du]);
    let lower = at(1000.0, [walls[0] + 0.5 * du, walls[1] - 0.5 * du]);
    let cross =
        (upper.flow().volumetric_flow_m3_s - lower.flow().volumetric_flow_m3_s) / (2.0 * du);
    close(
        base.flow().flow_relative_wall_velocity_derivative_m2,
        cross,
        3e-5,
    );
    let force =
        |p: &CarreauAnnularProfile| 0.5 * (p.flow().wall_forces_n[1] - p.flow().wall_forces_n[0]);
    close((force(&plus) - force(&minus)) / (2.0 * dp), -cross, 3e-5);
    close(
        base.flow().relative_wall_resistance_n_s_m,
        (force(&upper) - force(&lower)) / (2.0 * du),
        3e-5,
    );
}

#[test]
fn g2_true_annular_geometry_converges_to_the_nonlinear_slit_limit() {
    let source = law(false);
    let mut errors = Vec::new();
    for gap in [1e-2, 1e-3, 1e-4] {
        let owner = LaminarCarreauAnnulus::new(source, 300.0, 1.0, 1.0 + gap, 0.5, 1e6).unwrap();
        let h = owner.gap_m();
        let width = core::f64::consts::TAU * (1.0 + 0.5 * h);
        let planar = LaminarCarreauSlit::new(source, 300.0, width, h, 0.5, 1e6).unwrap();
        let walls = [-5.0 * h, 15.0 * h];
        let pressure = 40.0 / h;
        let curved = owner
            .profile_at_pressure_and_walls(pressure, walls, tight())
            .unwrap();
        let flat = planar
            .profile_at_pressure_and_walls(
                pressure,
                walls,
                CarreauSlitControls {
                    local: tight().local,
                    relative_wall_velocity_tolerance: 1e-11,
                    ..CarreauSlitControls::default()
                },
            )
            .unwrap();
        let error = ((curved.flow().volumetric_flow_m3_s - flat.flow().volumetric_flow_m3_s)
            / flat.flow().volumetric_flow_m3_s)
            .abs();
        errors.push(error);
        close(
            curved.flow().dissipated_power_w,
            flat.flow().dissipated_power_w,
            2.0 * h,
        );
    }
    assert!(errors[2] < 2e-4 && errors[2] < 0.1 * errors[0] && errors[2] < 0.5 * errors[1]);

    // A squared-radius subtraction would lose this narrow-gap mobility.
    let narrow = LaminarCarreauAnnulus::new(law(true), 300.0, 1.0, 1.0 + 1e-10, 0.5, 1e6).unwrap();
    let p = narrow
        .profile_at_pressure_and_walls(200.0, [0.0; 2], tight())
        .unwrap();
    let h = narrow.gap_m();
    let width = core::f64::consts::TAU * (1.0 + 0.5 * h);
    close(
        p.flow().pressure_conductance_m3_s_pa,
        width * h.powi(3) / (12.0 * 3.0 * 0.5),
        5e-12,
    );
    close(
        p.sample_gap_fraction(0.5).unwrap().axial_velocity_m_s,
        200.0 * h * h / (8.0 * 3.0 * 0.5),
        2e-8,
    );
}

#[test]
fn g3_common_translation_and_original_thermal_rate_shift_preserve_the_profile() {
    let owner = annulus();
    let walls = [-0.25, 0.5];
    let base = owner
        .profile_at_pressure_and_walls(1000.0, walls, tight())
        .unwrap();
    let speed = 1e8;
    let translated = owner
        .profile_at_pressure_and_walls(1000.0, [walls[0] + speed, walls[1] + speed], tight())
        .unwrap();
    assert_eq!(translated.stress_offset_pa(), base.stress_offset_pa());
    assert_eq!(
        translated.flow().relative_flow_m3_s,
        base.flow().relative_flow_m3_s
    );
    assert_eq!(
        translated.flow().dissipated_power_w,
        base.flow().dissipated_power_w
    );
    assert_eq!(
        translated.flow().work_residual_w,
        base.flow().work_residual_w
    );
    close(
        translated.flow().volumetric_flow_m3_s - base.flow().volumetric_flow_m3_s,
        owner.area_m2() * speed,
        2e-14,
    );

    // VerticalAndRate gives tau_T(g)=tau_ref(a_T*g). Holding stress while
    // dividing both wall speeds by a_T must divide flow, mobility and heat by
    // that SAME original shift, without fitting an annular apparent viscosity.
    let source = *owner.law();
    let temperature = 330.0;
    let shift = source
        .evaluate(temperature, 0.0)
        .unwrap()
        .temperature()
        .factor();
    let heated = LaminarCarreauAnnulus::new(source, temperature, 0.01, 0.04, 0.5, 1e6)
        .unwrap()
        .profile_at_pressure_and_walls(1000.0, [walls[0] / shift, walls[1] / shift], tight())
        .unwrap();
    close(heated.stress_offset_pa(), base.stress_offset_pa(), 3e-9);
    close(
        heated.flow().volumetric_flow_m3_s,
        base.flow().volumetric_flow_m3_s / shift,
        3e-9,
    );
    close(
        heated.flow().dissipated_power_w,
        base.flow().dissipated_power_w / shift,
        3e-9,
    );
    close(
        heated.flow().pressure_conductance_m3_s_pa,
        base.flow().pressure_conductance_m3_s_pa / shift,
        3e-9,
    );
}

#[test]
fn g4_curved_wall_source_limits_and_cancellation_are_not_repaired() {
    let limited = LaminarCarreauAnnulus::new(law(true), 300.0, 0.001, 0.1, 0.5, 100.0).unwrap();
    // Mean gap shear is only about 10 /s, but the inner cylinder exceeds 100 /s.
    assert!(
        limited
            .profile_at_pressure_and_walls(0.0, [-0.5, 0.5], tight())
            .is_err()
    );
    let owner = annulus();
    assert!(
        owner
            .profile_at_pressure_and_walls(1e12, [0.0; 2], tight())
            .is_err()
    );
    assert!(
        owner
            .profile_at_pressure_and_walls(
                1000.0,
                [-0.2, 0.6],
                CarreauAnnularControls {
                    maximum_wall_stress_iterations: 1,
                    ..tight()
                }
            )
            .is_err()
    );
    for radii in [[0.0, 0.1], [0.1, 0.1], [0.2, 0.1], [f64::NAN, 0.1]] {
        assert!(
            LaminarCarreauAnnulus::new(law(false), 300.0, radii[0], radii[1], 0.5, 100.0).is_err()
        );
    }
    let original = owner;
    let failure = state_error("injected annular cancellation");
    let mut calls = 0;
    let refused = owner.profile_with_check(1000.0, [-0.2, 0.6], tight(), || {
        calls += 1;
        if calls == 80 {
            Err(failure.clone())
        } else {
            Ok(())
        }
    });
    assert_eq!(refused.unwrap_err(), failure);
    assert_eq!(calls, 80);
    assert_eq!(owner, original);
    let retry = owner
        .profile_at_pressure_and_walls(1000.0, [-0.2, 0.6], tight())
        .unwrap();
    assert_eq!(
        retry,
        original
            .profile_at_pressure_and_walls(1000.0, [-0.2, 0.6], tight())
            .unwrap()
    );
    assert!(retry.sample_gap_fraction(-0.1).is_err());
    assert!(retry.sample_gap_fraction(f64::NAN).is_err());
    assert!(retry.sample(0.0).is_err());
}

#[test]
fn g1_prescribed_newtonian_flow_recovers_pressure_with_axial_wall_motion() {
    let (a, b, length, mu) = (0.01_f64, 0.04_f64, 0.5, 3.0);
    let owner = LaminarCarreauAnnulus::new(law(true), 300.0, a, b, length, 1e6).unwrap();
    let log = (b / a).ln();
    let area = core::f64::consts::PI * (b * b - a * a);
    let conductance = core::f64::consts::PI / (8.0 * mu * length)
        * (b.powi(4) - a.powi(4) - (b * b - a * a).powi(2) / log);
    for pressure in [-1000.0, 0.0, 1000.0] {
        for walls in [[0.1, 0.9], [-0.3, -0.1], [0.25, 0.25]] {
            let flow = area * walls[0]
                + core::f64::consts::PI
                    * (walls[1] - walls[0])
                    * (b * b - (b * b - a * a) / (2.0 * log))
                + conductance * pressure;
            let inverse = owner
                .profile_at_flow_and_walls(flow, walls, tight(), 1e-11, 1)
                .unwrap();
            let result = inverse.profile().flow();
            assert!((result.pressure_drop_pa - pressure).abs() <= 2e-10 * pressure.abs().max(1.0));
            assert!(inverse.flow_iterations() <= 1);
            assert_eq!(inverse.requested_volumetric_flow_m3_s(), flow);
            assert!(inverse.flow_residual_m3_s().abs() <= 1e-11 * flow.abs().max(1e-12));
            assert!(result.work_residual_w.abs() <= 1e-11 * result.dissipated_power_w.max(1.0));
        }
    }
}

#[test]
fn g2_prescribed_nonlinear_flow_roundtrips_blockage_and_reverse_pumping() {
    let owner = annulus();
    for (pressure, walls) in [
        (1000.0, [-0.2, 0.6]),
        (-1000.0, [0.2, -0.6]),
        (0.0, [0.0, 1.0]),
        (3000.0, [0.0, 0.0]),
    ] {
        let forward = owner
            .profile_at_pressure_and_walls(pressure, walls, tight())
            .unwrap();
        let target = forward.flow().volumetric_flow_m3_s;
        let inverse = owner
            .profile_at_flow_and_walls(target, walls, tight(), 1e-10, 96)
            .unwrap();
        let recovered = inverse.profile();
        assert!(
            (recovered.flow().pressure_drop_pa - pressure).abs() <= 2e-7 * pressure.abs().max(1.0)
        );
        close(
            recovered.flow().dissipated_power_w,
            forward.flow().dissipated_power_w,
            2e-9,
        );
        let scale = (target - forward.flow().mean_wall_entrainment_m3_s)
            .abs()
            .max(0.5 * owner.area_m2() * (walls[1] - walls[0]).abs());
        assert!(inverse.flow_residual_m3_s().abs() <= 1e-10 * scale);
        for fraction in [0.1, 0.5, 0.9] {
            let actual = recovered.sample_gap_fraction(fraction).unwrap();
            let expected = forward.sample_gap_fraction(fraction).unwrap();
            assert!(
                (actual.axial_velocity_m_s - expected.axial_velocity_m_s).abs()
                    <= 2e-9 * expected.axial_velocity_m_s.abs().max(1.0)
            );
        }
    }
    let blocked = owner
        .profile_at_flow_and_walls(0.0, [0.0, 1.0], tight(), 1e-10, 96)
        .unwrap();
    let reversed = owner
        .profile_at_flow_and_walls(-0.001, [0.0, 1.0], tight(), 1e-10, 96)
        .unwrap();
    assert!(blocked.profile().flow().pressure_drop_pa < 0.0);
    assert!(reversed.profile().flow().pressure_drop_pa < blocked.profile().flow().pressure_drop_pa);
    assert!(reversed.profile().flow().volumetric_flow_m3_s < 0.0);
    for solution in [blocked, reversed] {
        let f = solution.profile().flow();
        let scale = 0.5 * owner.area_m2();
        assert!(solution.flow_residual_m3_s().abs() <= 1e-10 * scale);
        assert!(f.dissipated_power_w > 0.0);
        assert!(f.work_residual_w.abs() <= 2e-8 * f.dissipated_power_w);
        close(
            f.pressure_power_w,
            f.pressure_drop_pa * f.volumetric_flow_m3_s,
            2e-14,
        );
    }
}

#[test]
fn g3_prescribed_flow_admits_pressure_redistribution_beyond_zero_pressure_wall_limit() {
    // Equal positive wall stresses maximize the supported relative wall speed.
    // Their interior stress is tau_w*(a*b+r²)/((a+b)*r). This independent
    // physical-radius construction has no zero-pressure admissible profile.
    let (a, b, length) = (0.001, 0.1, 0.5);
    let wall_stress = 0.6 * 1200.0 / 61.0_f64.sqrt();
    let stress = |r: f64| wall_stress * (a * b + r * r) / ((a + b) * r);
    let relative = simpson(|r| inverse(stress(r)), a, b);
    let flow = core::f64::consts::PI * simpson(|r| (b * b - r * r) * inverse(stress(r)), a, b);
    let pressure = -2.0 * length * wall_stress / (a + b);
    let owner = LaminarCarreauAnnulus::new(law(false), 300.0, a, b, length, 100.0).unwrap();
    for sign in [-1.0, 1.0] {
        let walls = [0.0, sign * relative];
        assert!(
            owner
                .profile_at_pressure_and_walls(0.0, walls, tight())
                .is_err()
        );
        let solved = owner
            .profile_at_flow_and_walls(sign * flow, walls, tight(), 1e-10, 96)
            .unwrap();
        close(
            solved.profile().flow().pressure_drop_pa,
            sign * pressure,
            2e-7,
        );
        let accepted = solved.profile().flow();
        assert!(accepted.maximum_shear_rate_per_s <= 100.0);
        close(accepted.wall_tractions_pa[0], -sign * wall_stress, 2e-7);
        close(accepted.wall_tractions_pa[1], sign * wall_stress, 2e-7);
    }
}

#[test]
fn equal_wall_limit_reconstructs_both_walls_inside_every_neighbouring_ceiling() {
    // a=0.001, b=0.1: z=0.98, so 1-z=0.0198 turns every rounding of z into 49.5
    // roundings of the inner wall. The ceiling itself reconstructs one wall
    // above itself for a share of the 20000 consecutive doubles around
    // 1200/sqrt(61) Pa; every one is admitted with both walls inside, and none
    // is lowered by more than twice the stated bound.
    let owner = LaminarCarreauAnnulus::new(law(false), 300.0, 0.001, 0.1, 0.5, 100.0).unwrap();
    let z = 0.5 * (owner.gap_m / owner.mean_radius_m);
    let kappa = 4.0 * z / (1.0 - z) + 4.0;
    let mut ceiling = 1200.0 / 61.0_f64.sqrt();
    for _ in 0..10_000 {
        ceiling = ceiling.next_down();
    }
    let mut lowered = 0usize;
    for _ in 0..20_000 {
        let limit = owner.equal_wall_limit(ceiling, z).unwrap();
        for wall in owner.wall_stresses(limit, -z * limit) {
            assert!(wall.is_finite() && wall.abs() <= ceiling, "{ceiling:e}");
        }
        assert!(
            limit >= ceiling * (1.0 - 2.0 * kappa * f64::EPSILON),
            "{ceiling:e}"
        );
        lowered += usize::from(limit < ceiling);
        ceiling = ceiling.next_up();
    }
    assert!(
        lowered > 0,
        "the sweep must contain ceilings that need the bound"
    );
}

#[test]
fn equal_wall_limit_still_refuses_a_reconstruction_error_beyond_its_bound() {
    // z perturbed by 1e-12 relative is 4500 times the unit roundoff, where the
    // bound covers 4 roundings of z: the inner wall then sits 4.95e-11 of the
    // ceiling above it, 1100 times the lowering, so the verification refuses.
    let owner = LaminarCarreauAnnulus::new(law(false), 300.0, 0.001, 0.1, 0.5, 100.0).unwrap();
    let z = 0.5 * (owner.gap_m / owner.mean_radius_m);
    let ceiling = 1200.0 / 61.0_f64.sqrt();
    assert_eq!(
        owner
            .equal_wall_limit(ceiling, z * (1.0 - 1e-12))
            .unwrap_err(),
        state_error("Carreau annular equal-wall stress span is unrepresentable")
    );
}

#[test]
fn g4_prescribed_annular_flow_refuses_unsupported_targets_and_retries_after_cancellation() {
    let owner = annulus();
    for (tolerance, iterations) in [(0.0, 96), (f64::NAN, 96), (1e-9, 0), (1e-9, 129)] {
        assert!(
            owner
                .profile_at_flow_and_walls(0.001, [0.0, 0.0], tight(), tolerance, iterations,)
                .is_err()
        );
    }
    assert!(
        owner
            .profile_at_flow_and_walls(f64::NAN, [0.0, 0.0], tight(), 1e-9, 96)
            .is_err()
    );
    let limited = LaminarCarreauAnnulus::new(law(false), 300.0, 0.001, 0.1, 0.5, 100.0).unwrap();
    assert!(
        limited
            .profile_at_flow_and_walls(0.0, [0.0, 100.0 * limited.gap_m()], tight(), 1e-9, 96,)
            .is_err()
    );
    assert!(
        limited
            .profile_at_flow_and_walls(1.0, [0.0, 0.0], tight(), 1e-9, 96,)
            .is_err()
    );
    let target = owner
        .profile_at_pressure_and_walls(3000.0, [0.0, 0.0], tight())
        .unwrap()
        .flow()
        .volumetric_flow_m3_s;
    assert!(
        owner
            .profile_at_flow_and_walls(target, [0.0, 0.0], tight(), 1e-10, 1)
            .is_err()
    );

    let original = owner;
    let failure = state_error("injected prescribed-annular-flow cancellation");
    let mut calls = 0;
    let refused =
        owner.profile_at_flow_and_walls_with_check(0.002, [-0.2, 0.6], tight(), 1e-10, 96, || {
            calls += 1;
            if calls == 80 {
                Err(failure.clone())
            } else {
                Ok(())
            }
        });
    assert_eq!(refused.unwrap_err(), failure);
    assert_eq!(calls, 80);
    assert_eq!(owner, original);
    let retry = owner
        .profile_at_flow_and_walls(0.002, [-0.2, 0.6], tight(), 1e-10, 96)
        .unwrap();
    assert_eq!(
        retry,
        original
            .profile_at_flow_and_walls(0.002, [-0.2, 0.6], tight(), 1e-10, 96,)
            .unwrap()
    );
}
