//! Mathematical controls only; none of these synthetic laws calibrate a fluid.
use super::*;
use crate::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= relative * scale.max(1e-25),
        "actual {actual:e}, expected {expected:e}, tolerance {relative:e}"
    );
}

fn law(p: [f64; 5], support: [f64; 2], mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(p[0], p[1], p[2], p[3], p[4], support).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}

fn slit() -> LaminarCarreauSlit {
    LaminarCarreauSlit::new(
        law(
            [12.0, 0.0, 0.6, 0.5, 1.0],
            [0.0, 1e8],
            ViscosityShiftMode::Vertical,
        ),
        300.0,
        0.07,
        0.02,
        0.5,
        1e6,
    )
    .unwrap()
}

fn tight() -> CarreauSlitControls {
    CarreauSlitControls {
        local: CarreauPipeControls {
            relative_constitutive_tolerance: 1e-13,
            relative_quadrature_tolerance: 1e-10,
            ..CarreauPipeControls::default()
        },
        relative_wall_velocity_tolerance: 1e-11,
        ..CarreauSlitControls::default()
    }
}

#[test]
fn g1_newtonian_pressure_sliding_traction_heat_and_tangents_are_one_profile() {
    for p in [
        [3.0, 3.0, 0.6, 0.4, 0.6],
        [3.0, 0.1, 0.0, 0.4, 1.4],
        [3.0, 0.0, 0.7, 1.0, 0.3],
    ] {
        let source = law(p, [0.0, 1e8], ViscosityShiftMode::VerticalAndRate);
        let film = LaminarCarreauSlit::new(source, 330.0, 0.07, 0.02, 0.5, 1e6).unwrap();
        let mu = source
            .evaluate(330.0, 0.0)
            .unwrap()
            .dynamic_viscosity_pa_s();
        let (w, h, l) = (film.width_m(), film.gap_m(), film.length_m());
        assert_eq!(film.zero_shear_viscosity_pa_s(), mu);
        let rebound = film.with_gap(1.5 * h).unwrap();
        assert_eq!(rebound.gap_m(), 1.5 * h);
        assert_eq!(rebound.law(), film.law());
        assert_eq!(rebound.temperature_k(), film.temperature_k());
        assert_eq!(rebound.width_m(), film.width_m());
        assert_eq!(rebound.length_m(), film.length_m());
        assert_eq!(rebound.zero_shear_viscosity_pa_s(), mu);
        assert_eq!(
            rebound.maximum_laminar_shear_rate_per_s(),
            film.maximum_laminar_shear_rate_per_s()
        );
        assert!(film.with_gap(0.0).is_err());
        assert!(film.with_gap(f64::INFINITY).is_err());
        for dp in [-4000.0, 0.0, 4000.0] {
            for walls in [[-0.25, 0.75], [2.0, 2.0], [0.0, 0.0]] {
                let profile = film
                    .profile_at_pressure_and_walls(dp, walls, CarreauSlitControls::default())
                    .unwrap();
                let f = profile.flow();
                let relative = walls[1] - walls[0];
                let mean = 0.5 * (walls[0] + walls[1]);
                let c = mu * relative / h;
                close(profile.stress_offset_pa(), c, 5e-14);
                let s = dp * h / (2.0 * l);
                let g = w * h.powi(3) / (12.0 * mu * l);
                close(f.pressure_flow_m3_s, g * dp, 5e-14);
                close(f.volumetric_flow_m3_s, w * h * mean + g * dp, 5e-14);
                close(f.wall_tractions_pa[0], -c - s, 5e-14);
                close(f.wall_tractions_pa[1], c - s, 5e-14);
                close(
                    f.dissipated_power_w,
                    g * dp * dp + w * l * mu * relative * relative / h,
                    6e-14,
                );
                close(f.pressure_conductance_m3_s_pa, g, 5e-14);
                close(f.relative_wall_resistance_n_s_m, w * l * mu / h, 5e-14);
                assert_eq!(f.flow_relative_wall_velocity_derivative_m2, 0.0);
                assert_eq!(f.wall_stress_iterations, 0);
                for y in [0.0, 0.13 * h, 0.5 * h, 0.97 * h, h] {
                    let sample = profile.sample(y).unwrap();
                    let gamma = relative / h + dp * (0.5 * h - y) / (mu * l);
                    close(
                        sample.axial_velocity_m_s,
                        walls[0] + relative * y / h + dp * y * (h - y) / (2.0 * mu * l),
                        2e-12,
                    );
                    close(sample.velocity_gradient_per_s, gamma, 2e-12);
                    close(sample.dissipation_w_m3, mu * gamma * gamma, 5e-12);
                }
            }
        }
    }
}

#[test]
fn g1_pure_couette_uses_the_actual_source_rate_even_at_a_singular_scalar_slope() {
    let source = law(
        [12.0, 0.2, 0.6, 0.4, 0.3],
        [0.0, 1e8],
        ViscosityShiftMode::Vertical,
    );
    let film = LaminarCarreauSlit::new(source, 300.0, 0.07, 0.02, 0.5, 1e6).unwrap();
    for relative in [-0.75_f64, 0.0, 0.75] {
        let p = film
            .profile_at_pressure_and_walls(0.0, [-relative / 2.0, relative / 2.0], tight())
            .unwrap();
        let f = p.flow();
        let material = source
            .evaluate(300.0, relative.abs() / film.gap_m())
            .unwrap();
        let mu = material.dynamic_viscosity_pa_s();
        let tangent = material.stress_tangent_pa_s();
        close(f.wall_tractions_pa[1], mu * relative / film.gap_m(), 2e-14);
        close(
            f.dissipated_power_w,
            film.width_m() * film.length_m() * mu * relative.powi(2) / film.gap_m(),
            2e-14,
        );
        close(
            f.pressure_conductance_m3_s_pa,
            film.width_m() * film.gap_m().powi(3) / (12.0 * film.length_m() * tangent),
            2e-14,
        );
        assert_eq!(f.pressure_flow_m3_s, 0.0);
        assert_eq!(f.wall_stress_iterations, 0);
    }
}

// Independent exact inverse for a=1, n=1/2, mu_inf=0. The production code
// uses the canonical Carreau stress tangent and iterative inverse instead.
fn inverse(stress: f64) -> f64 {
    let s = stress.abs() * 0.6 / 12.0;
    ((s * s + s * (s * s + 4.0).sqrt()) / 1.2).copysign(stress)
}

fn simpson(f: impl Fn(f64) -> f64, end: f64) -> f64 {
    let n = 4096;
    let mut sum = 0.0;
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * f(end * f64::from(i) / f64::from(n));
    }
    end * sum / (3.0 * f64::from(n))
}

fn oracle(dp: f64, walls: [f64; 2]) -> (f64, f64, f64) {
    let film = slit();
    let s = dp * film.gap_m() / (2.0 * film.length_m());
    let target = (walls[1] - walls[0]) / film.gap_m();
    let (mut lo, mut hi) = (-1000.0, 1000.0);
    for _ in 0..70 {
        let c = 0.5 * (lo + hi);
        let mean = simpson(|t| inverse(c + s * (1.0 - 2.0 * t)), 1.0);
        if mean < target {
            lo = c;
        } else {
            hi = c;
        }
    }
    let c = if target == 0.0 { 0.0 } else { 0.5 * (lo + hi) };
    let excess = film.width_m()
        * film.gap_m().powi(2)
        * simpson(|t| (0.5 - t) * inverse(c + s * (1.0 - 2.0 * t)), 1.0);
    let heat = film.width_m()
        * film.length_m()
        * film.gap_m()
        * simpson(
            |t| {
                let tau = c + s * (1.0 - 2.0 * t);
                tau * inverse(tau)
            },
            1.0,
        );
    (c, excess, heat)
}

#[test]
fn g1_combined_nonlinear_flow_and_heat_match_independent_momentum_quadrature() {
    let film = slit();
    for (dp, walls) in [
        (4000.0, [-0.2, 0.5]),
        (-4000.0, [0.5, -0.2]),
        (5000.0, [0.0, 0.0]),
    ] {
        let p = film
            .profile_at_pressure_and_walls(dp, walls, tight())
            .unwrap();
        let f = p.flow();
        let (c, q, heat) = oracle(dp, walls);
        close(f.pressure_flow_m3_s, q, 3e-9);
        close(f.dissipated_power_w, heat, 3e-9);
        close(
            0.5 * (f.wall_tractions_pa[1] - f.wall_tractions_pa[0]),
            c,
            3e-9,
        );
        let requested_work = f.pressure_power_w + f.wall_power_w;
        let work_scale = f.pressure_power_w.abs() + f.wall_power_w.abs() + f.dissipated_power_w;
        assert!(
            (requested_work - f.dissipated_power_w - f.work_residual_w).abs()
                <= 32.0 * f64::EPSILON * work_scale
        );
        assert!(f.wall_velocity_residual_m_s.abs() <= 2e-10 * (walls[1] - walls[0]).abs());
        let s = dp * film.gap_m() / (2.0 * film.length_m());
        for t in [0.0, 0.13, 0.5, 0.87, 1.0] {
            let sample = p.sample(t * film.gap_m()).unwrap();
            let expected =
                walls[0] + film.gap_m() * simpson(|x| inverse(c + s * (1.0 - 2.0 * x)), t);
            assert!((sample.axial_velocity_m_s - expected).abs() <= 2e-8 * expected.abs().max(1.0));
            close(
                sample.velocity_gradient_per_s,
                inverse(c + s * (1.0 - 2.0 * t)),
                2e-9,
            );
        }
    }
}

#[test]
fn g1_pressure_and_sliding_are_not_separately_superposed() {
    let film = slit();
    let combined = film
        .profile_at_pressure_and_walls(5000.0, [-1.0, 1.0], tight())
        .unwrap();
    let pressure = film
        .profile_at_pressure_and_walls(5000.0, [0.0, 0.0], tight())
        .unwrap();
    let sliding = film
        .profile_at_pressure_and_walls(0.0, [-1.0, 1.0], tight())
        .unwrap();
    assert!(
        (combined.flow().pressure_flow_m3_s / pressure.flow().pressure_flow_m3_s - 1.0).abs() > 0.1
    );
    assert!(
        (combined.flow().dissipated_power_w
            / (pressure.flow().dissipated_power_w + sliding.flow().dissipated_power_w)
            - 1.0)
            .abs()
            > 0.01
    );
}

#[test]
fn g1_analytic_mobility_and_work_reciprocity_match_solved_directional_derivatives() {
    let film = slit();
    let dp = 4000.0;
    let delta = 0.7;
    let at = |p, d| {
        film.profile_at_pressure_and_walls(p, [-0.5 * d, 0.5 * d], tight())
            .unwrap()
    };
    let base = at(dp, delta);
    let dp_step = 0.4;
    let plus = at(dp + dp_step, delta);
    let minus = at(dp - dp_step, delta);
    close(
        base.flow().pressure_conductance_m3_s_pa,
        (plus.flow().pressure_flow_m3_s - minus.flow().pressure_flow_m3_s) / (2.0 * dp_step),
        2e-7,
    );
    let d_step = 7e-5;
    let plus_d = at(dp, delta + d_step);
    let minus_d = at(dp, delta - d_step);
    close(
        base.flow().flow_relative_wall_velocity_derivative_m2,
        (plus_d.flow().pressure_flow_m3_s - minus_d.flow().pressure_flow_m3_s) / (2.0 * d_step),
        2e-7,
    );
    let force = |p: &CarreauSlitProfile| {
        film.width_m()
            * film.length_m()
            * 0.5
            * (p.flow().wall_tractions_pa[1] - p.flow().wall_tractions_pa[0])
    };
    close(
        base.flow().relative_wall_resistance_n_s_m,
        (force(&plus_d) - force(&minus_d)) / (2.0 * d_step),
        2e-7,
    );
    close(
        -base.flow().flow_relative_wall_velocity_derivative_m2,
        (force(&plus) - force(&minus)) / (2.0 * dp_step),
        2e-7,
    );
    assert!(base.flow().pressure_conductance_m3_s_pa > 0.0);
    assert!(base.flow().relative_wall_resistance_n_s_m > 0.0);
}

#[test]
fn g3_common_translation_preserves_heat_and_resolves_small_pressure_flow() {
    let film = slit();
    let dp = 4000.0;
    let base = film
        .profile_at_pressure_and_walls(dp, [-0.25, 0.5], tight())
        .unwrap();
    let speed = 1048576.0;
    let moved = film
        .profile_at_pressure_and_walls(dp, [speed - 0.25, speed + 0.5], tight())
        .unwrap();
    assert_eq!(
        moved.flow().pressure_flow_m3_s,
        base.flow().pressure_flow_m3_s
    );
    assert_eq!(
        moved.flow().dissipated_power_w,
        base.flow().dissipated_power_w
    );
    assert_eq!(
        moved.flow().wall_tractions_pa,
        base.flow().wall_tractions_pa
    );
    close(
        moved.flow().volumetric_flow_m3_s - base.flow().volumetric_flow_m3_s,
        film.width_m() * film.gap_m() * speed,
        5e-15,
    );
    close(
        moved.flow().wall_power_w - base.flow().wall_power_w,
        -dp * film.width_m() * film.gap_m() * speed,
        5e-15,
    );
    let reverse = film
        .profile_at_pressure_and_walls(-dp, [0.25, -0.5], tight())
        .unwrap();
    close(
        reverse.flow().volumetric_flow_m3_s,
        -base.flow().volumetric_flow_m3_s,
        1e-12,
    );
    close(
        reverse.flow().dissipated_power_w,
        base.flow().dissipated_power_w,
        1e-12,
    );
    close(
        reverse.flow().pressure_conductance_m3_s_pa,
        base.flow().pressure_conductance_m3_s_pa,
        1e-12,
    );
}

#[test]
fn g1_both_thermal_shift_modes_retain_the_original_temperature_calibration() {
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let source = law([12.0, 0.2, 0.6, 0.45, 2.0], [0.0, 1e8], mode);
        let a = source.temperature_shift().evaluate(330.0).unwrap().factor();
        let heated = LaminarCarreauSlit::new(source, 330.0, 0.07, 0.02, 0.5, 1e6).unwrap();
        let time = if mode == ViscosityShiftMode::Vertical {
            0.6
        } else {
            0.6 * a
        };
        let shifted = LaminarCarreauSlit::new(
            law(
                [12.0 * a, 0.2 * a, time, 0.45, 2.0],
                [0.0, 1e8],
                ViscosityShiftMode::Vertical,
            ),
            300.0,
            0.07,
            0.02,
            0.5,
            1e6,
        )
        .unwrap();
        let p = heated
            .profile_at_pressure_and_walls(4000.0, [-0.2, 0.5], tight())
            .unwrap();
        let q = shifted
            .profile_at_pressure_and_walls(4000.0, [-0.2, 0.5], tight())
            .unwrap();
        close(
            p.flow().pressure_flow_m3_s,
            q.flow().pressure_flow_m3_s,
            3e-11,
        );
        close(
            p.flow().dissipated_power_w,
            q.flow().dissipated_power_w,
            3e-11,
        );
        close(
            p.flow().pressure_conductance_m3_s_pa,
            q.flow().pressure_conductance_m3_s_pa,
            3e-11,
        );
        assert_eq!(heated.law(), &source);
        assert_eq!(heated.temperature_k(), 330.0);
    }
}

#[test]
fn g4_source_support_invalid_controls_and_exhaustion_refuse_without_a_candidate() {
    let film = slit();
    let controls = CarreauSlitControls::default();
    assert!(LaminarCarreauSlit::new(*film.law(), 361.0, 0.07, 0.02, 0.5, 1e6).is_err());
    assert!(
        LaminarCarreauSlit::new(
            law(
                [12.0, 0.0, 0.6, 0.5, 1.0],
                [1.0, 1e8],
                ViscosityShiftMode::Vertical
            ),
            300.0,
            0.07,
            0.02,
            0.5,
            1e6
        )
        .is_err()
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LaminarCarreauSlit::new(*film.law(), 300.0, bad, 0.02, 0.5, 1e6).is_err());
    }
    assert!(
        film.profile_at_pressure_and_walls(1e20, [0.0, 0.0], controls)
            .is_err()
    );
    assert!(
        film.profile_at_pressure_and_walls(0.0, [0.0, 1e20], controls)
            .is_err()
    );
    assert!(
        film.profile_at_pressure_and_walls(f64::NAN, [0.0, 0.0], controls)
            .is_err()
    );
    for bad in [0.0, -1.0, 1.0, f64::NAN] {
        assert!(
            film.profile_at_pressure_and_walls(
                4000.0,
                [-0.2, 0.5],
                CarreauSlitControls {
                    relative_wall_velocity_tolerance: bad,
                    ..controls
                }
            )
            .is_err()
        );
    }
    assert!(
        film.profile_at_pressure_and_walls(
            4000.0,
            [-0.2, 0.5],
            CarreauSlitControls {
                maximum_wall_stress_iterations: 1,
                ..controls
            }
        )
        .is_err()
    );
    assert!(
        film.profile_at_pressure_and_walls(
            4000.0,
            [0.0, 0.0],
            CarreauSlitControls {
                local: CarreauPipeControls {
                    maximum_quadrature_refinements: 1,
                    ..controls.local
                },
                ..controls
            }
        )
        .is_err()
    );
    let p = film
        .profile_at_pressure_and_walls(4000.0, [-0.2, 0.5], controls)
        .unwrap();
    assert!(p.sample(-1e-10).is_err());
    assert!(p.sample(film.gap_m() + 1e-10).is_err());
    assert!(p.sample(f64::NAN).is_err());
}

#[test]
fn g4_late_checked_refusal_leaves_the_immutable_owner_exactly_retryable() {
    let film = slit();
    let initial = film;
    let mut checks = 0;
    let failure = state_error("injected caller cancellation");
    let result = film.profile_with_check(4000.0, [-0.2, 0.5], tight(), || {
        checks += 1;
        if checks == 100 {
            Err(failure.clone())
        } else {
            Ok(())
        }
    });
    assert_eq!(result.unwrap_err(), failure);
    assert_eq!(checks, 100);
    assert_eq!(film, initial);
    let retry = film
        .profile_at_pressure_and_walls(4000.0, [-0.2, 0.5], tight())
        .unwrap();
    let direct = initial
        .profile_at_pressure_and_walls(4000.0, [-0.2, 0.5], tight())
        .unwrap();
    assert_eq!(retry, direct);
    let mut sample_checks = 0;
    let result = retry.sample_with_check(0.37 * film.gap_m(), || {
        sample_checks += 1;
        if sample_checks == 7 {
            Err(failure.clone())
        } else {
            Ok(())
        }
    });
    assert_eq!(result.unwrap_err(), failure);
    assert_eq!(sample_checks, 7);
    assert_eq!(
        retry.sample(0.37 * film.gap_m()).unwrap(),
        direct.sample(0.37 * film.gap_m()).unwrap()
    );
}
