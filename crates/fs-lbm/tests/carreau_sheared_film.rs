//! G1/G3: prescribed surface traction drives one original nonlinear film.
//! Synthetic constitutive parameters; momentum balance and power integrals
//! supply independent references, including reversal and zero net flux.

use fs_lbm::film::carreau::{
    CarreauFilm, CarreauFilmControls, CarreauFilmError, CarreauFilmSpec, CarreauShearedFilm,
    CarreauShearedFilmSpec,
};
use fs_lbm::film::{FilmMaterial, ShearedFilm};
use fs_material::fluid_transport::rheology::slit::{CarreauSlitControls, LaminarCarreauSlit};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const CURVE: [f64; 5] = [0.08, 0.006, 0.4, 0.45, 2.0];

fn law(parameters: [f64; 5], support: [f64; 2]) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
            support,
        )
        .unwrap(),
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

fn controls() -> CarreauSlitControls {
    let mut value = CarreauSlitControls::default();
    value.local.relative_constitutive_tolerance = 1e-12;
    value.local.relative_quadrature_tolerance = 1e-10;
    // Traction determines stress directly; it must not need this wall-speed
    // iteration budget even for a strongly varying viscosity.
    value.maximum_wall_stress_iterations = 1;
    value
}

fn near(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= relative * scale.max(expected.abs()).max(1e-30),
        "{actual:e} != {expected:e}, scale={scale:e}, tolerance={relative:e}"
    );
}

fn film(spec: CarreauShearedFilmSpec) -> CarreauShearedFilm {
    CarreauShearedFilm::new(law(CURVE, [0.0, 1e6]), 300.0, spec, 1e4, controls()).unwrap()
}

#[test]
fn constant_source_recovers_signed_newtonian_profiles_and_zero_flux_heat() {
    let source = law([0.02, 0.001, 0.5, 1.0, 2.0], [0.0, 1e6]);
    let temperature = 330.0;
    let mu = source
        .evaluate(temperature, 0.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let (rho, h) = (1000.0, 0.001);
    for (g, traction) in [
        (1.7, 0.0),
        (1.7, 0.4),
        (1.7, -0.4),
        (1.7, -2.0 * rho * 1.7 * h / 3.0),
        (0.0, 0.4),
        (0.0, -0.4),
        (-1.7, 0.4),
        (0.0, 0.0),
    ] {
        let spec = CarreauShearedFilmSpec {
            density_kg_m3: rho,
            thickness_m: h,
            tangential_acceleration_m_s2: g,
            surface_shear_pa: traction,
        };
        let actual = CarreauShearedFilm::new(source, temperature, spec, 1e5, controls()).unwrap();
        let exact = ShearedFilm::new(
            FilmMaterial {
                density_kg_m3: rho,
                kinematic_viscosity_m2_s: mu / rho,
                surface_tension_n_m: 0.05,
            },
            h,
            g,
            traction,
        )
        .unwrap();
        let stress_scale = (rho * g * h).abs().max(traction.abs()).max(1e-12);
        let velocity_scale = stress_scale * h / mu;
        let flux_scale = velocity_scale * h;
        let power_scale = stress_scale * velocity_scale;
        near(
            actual.flow_per_width_m2_s(),
            exact.flow_per_width_m2_s(),
            flux_scale,
            2e-11,
        );
        near(
            actual.surface_velocity_m_s(),
            exact.surface_velocity_m_s(),
            velocity_scale,
            2e-11,
        );
        near(
            actual.wall_shear_pa(),
            exact.wall_shear_pa(),
            stress_scale,
            2e-12,
        );
        near(actual.surface_shear_residual_pa(), 0.0, stress_scale, 2e-12);
        for fraction in [0.0, 0.2, 0.7, 1.0] {
            let y = fraction * h;
            let sample = actual.sample(y).unwrap();
            near(
                sample.velocity_m_s,
                exact.velocity_m_s(y).unwrap(),
                velocity_scale,
                2e-11,
            );
            near(
                sample.shear_stress_pa,
                exact.shear_pa(y).unwrap(),
                stress_scale,
                2e-12,
            );
            near(sample.dynamic_viscosity_pa_s, mu, mu, 2e-12);
            near(sample.stress_tangent_pa_s, mu, mu, 2e-12);
            near(
                sample.velocity_curvature_per_m_s,
                -rho * g / mu,
                stress_scale / (h * mu),
                2e-12,
            );
        }
        let a = actual.power();
        let b = exact.power();
        near(a.body_work_w_m2, b.body_work_w_m2, power_scale, 2e-11);
        near(a.surface_work_w_m2, b.surface_work_w_m2, power_scale, 2e-11);
        near(
            a.viscous_dissipation_w_m2,
            b.viscous_dissipation_w_m2,
            power_scale,
            2e-11,
        );
        near(
            a.body_work_w_m2 + a.surface_work_w_m2,
            a.viscous_dissipation_w_m2,
            power_scale,
            2e-11,
        );
        near(actual.work_residual_w_m2(), 0.0, power_scale, 2e-12);
        near(
            actual.surface_velocity_traction_derivative_m_s_pa(),
            h / mu,
            h / mu,
            2e-12,
        );
        near(
            actual.flow_traction_derivative_m2_s_pa(),
            h * h / (2.0 * mu),
            h * h / mu,
            2e-12,
        );
        near(
            actual.flow_force_density_derivative_m3_s_pa(),
            h * h * h / (3.0 * mu),
            h * h * h / mu,
            2e-12,
        );
        if g == 1.7 && traction < -1.0 {
            near(actual.flow_per_width_m2_s(), 0.0, flux_scale, 2e-12);
            assert!(a.viscous_dissipation_w_m2 > 0.0);
        }
    }
}

// Independent platform-powf curve and bisection of mu(gamma)*gamma=|tau|.
// No production rheology, slit inverse, profile or quadrature is called here.
fn point_reference(stress: f64) -> (f64, f64, f64) {
    let [zero, infinity, time, n, a] = CURVE;
    let response = |rate: f64| {
        let x = (time * rate).powf(a);
        let excess = (zero - infinity) * (1.0 + x).powf((n - 1.0) / a);
        (
            infinity + excess,
            infinity + excess * (1.0 + (n - 1.0) * x / (1.0 + x)),
        )
    };
    let (mut low, mut high) = (0.0, stress.abs() / infinity);
    for _ in 0..80 {
        let rate = 0.5 * (low + high);
        if response(rate).0 * rate > stress.abs() {
            high = rate;
        } else {
            low = rate;
        }
    }
    let rate = 0.5 * (low + high);
    let (mu, tangent) = response(rate);
    (rate.copysign(stress), mu, tangent)
}

// Integrals [surface velocity, flux, heat, dUs/dtau, dq/dtau, dq/dF].
// q=int(h-y)*gamma dy and heat=int(tau*gamma)dy. The same independent
// differential viscosity gives the positive traction/body-force mobilities.
fn integral_reference(spec: CarreauShearedFilmSpec, upper: f64) -> [f64; 6] {
    let panels = 2048;
    let dy = upper / panels as f64;
    let force = spec.density_kg_m3 * spec.tangential_acceleration_m_s2;
    let mut total = [0.0; 6];
    for i in 0..=panels {
        let y = i as f64 * dy;
        let depth = spec.thickness_m - y;
        let stress = spec.surface_shear_pa + force * depth;
        let (rate, _, tangent) = point_reference(stress);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in total.iter_mut().zip([
            rate,
            depth * rate,
            stress * rate,
            1.0 / tangent,
            depth / tangent,
            depth * depth / tangent,
        ]) {
            *sum += weight * value;
        }
    }
    total.map(|value| value * dy / 3.0)
}

#[test]
fn nonlinear_counterflow_heat_and_mobilities_match_continuum_stress_integrals() {
    for traction in [0.7, -0.8, -2.0, -4.2] {
        let spec = CarreauShearedFilmSpec {
            density_kg_m3: 950.0,
            thickness_m: 0.002,
            tangential_acceleration_m_s2: 2.0,
            surface_shear_pa: traction,
        };
        let actual = film(spec);
        let exact = integral_reference(spec, spec.thickness_m);
        let force = spec.density_kg_m3 * spec.tangential_acceleration_m_s2;
        let stress_scale = traction
            .abs()
            .max((traction + force * spec.thickness_m).abs());
        let velocity_scale = stress_scale * spec.thickness_m / CURVE[1];
        near(
            actual.surface_velocity_m_s(),
            exact[0],
            velocity_scale,
            3e-8,
        );
        near(
            actual.flow_per_width_m2_s(),
            exact[1],
            velocity_scale * spec.thickness_m,
            3e-8,
        );
        near(
            actual.power().viscous_dissipation_w_m2,
            exact[2],
            exact[2],
            3e-8,
        );
        for (value, expected) in [
            actual.surface_velocity_traction_derivative_m_s_pa(),
            actual.flow_traction_derivative_m2_s_pa(),
            actual.flow_force_density_derivative_m3_s_pa(),
        ]
        .into_iter()
        .zip(&exact[3..])
        {
            near(value, *expected, expected.abs(), 3e-8);
        }
        for fraction in [0.0, 0.25, 0.6, 1.0] {
            let y = fraction * spec.thickness_m;
            let stress = traction + force * (spec.thickness_m - y);
            let (rate, mu, tangent) = point_reference(stress);
            let sample = actual.sample(y).unwrap();
            near(
                sample.velocity_m_s,
                integral_reference(spec, y)[0],
                velocity_scale,
                3e-8,
            );
            near(
                sample.velocity_gradient_per_s,
                rate,
                velocity_scale / spec.thickness_m,
                2e-10,
            );
            near(sample.dynamic_viscosity_pa_s, mu, mu, 2e-10);
            near(sample.stress_tangent_pa_s, tangent, tangent, 2e-10);
            near(
                sample.velocity_curvature_per_m_s,
                -force / tangent,
                force.abs() / tangent,
                2e-10,
            );
            near(
                sample.dissipation_w_m3,
                stress * rate,
                stress_scale * velocity_scale / spec.thickness_m,
                2e-10,
            );
            assert!(sample.dissipation_w_m3 >= 0.0);
        }
        let power = actual.power();
        near(
            power.body_work_w_m2 + power.surface_work_w_m2,
            power.viscous_dissipation_w_m2,
            power.viscous_dissipation_w_m2,
            2e-10,
        );
        near(
            actual.work_residual_w_m2(),
            0.0,
            power.viscous_dissipation_w_m2,
            2e-12,
        );
        if traction == -0.8 {
            assert!(power.surface_work_w_m2 < 0.0);
        }
        if traction == -2.0 {
            assert!(actual.sample(0.25 * spec.thickness_m).unwrap().velocity_m_s > 0.0);
            assert!(actual.surface_velocity_m_s() < 0.0);
        }
        if traction == -4.2 {
            assert!(power.body_work_w_m2 < 0.0);
        }
        assert!(exact[3] * exact[5] > exact[4] * exact[4]);
        if traction == -0.8 {
            let delta = 1e-4;
            let plus = film(CarreauShearedFilmSpec {
                surface_shear_pa: traction + delta,
                ..spec
            });
            let minus = film(CarreauShearedFilmSpec {
                surface_shear_pa: traction - delta,
                ..spec
            });
            near(
                (plus.flow_per_width_m2_s() - minus.flow_per_width_m2_s()) / (2.0 * delta),
                actual.flow_traction_derivative_m2_s_pa(),
                exact[4],
                2e-5,
            );
            let delta_force = 0.1;
            let plus = film(CarreauShearedFilmSpec {
                tangential_acceleration_m_s2: spec.tangential_acceleration_m_s2
                    + delta_force / spec.density_kg_m3,
                ..spec
            });
            let minus = film(CarreauShearedFilmSpec {
                tangential_acceleration_m_s2: spec.tangential_acceleration_m_s2
                    - delta_force / spec.density_kg_m3,
                ..spec
            });
            near(
                (plus.surface_velocity_m_s() - minus.surface_velocity_m_s()) / (2.0 * delta_force),
                actual.flow_traction_derivative_m2_s_pa(),
                exact[4],
                2e-5,
            );
            near(
                (plus.flow_per_width_m2_s() - minus.flow_per_width_m2_s()) / (2.0 * delta_force),
                actual.flow_force_density_derivative_m3_s_pa(),
                exact[5],
                2e-5,
            );
        }
    }
}

#[test]
fn zero_traction_continues_the_original_prescribed_flux_film() {
    let source = law(CURVE, [0.0, 1e6]);
    let base = CarreauFilm::at_flow(
        source,
        300.0,
        CarreauFilmSpec {
            density_kg_m3: 1000.0,
            surface_tension_n_m: 0.05,
            flow_per_width_m2_s: 1e-5,
            tangential_acceleration_m_s2: 3.0,
            normal_acceleration_m_s2: 3.0,
        },
        1e4,
        CarreauFilmControls {
            slit: controls(),
            ..CarreauFilmControls::default()
        },
    )
    .unwrap();
    let sheared = base.with_surface_shear(0.0, controls()).unwrap();
    assert_eq!(sheared.law(), base.law());
    assert_eq!(sheared.thickness_m(), base.thickness_m());
    near(
        sheared.flow_per_width_m2_s(),
        base.flow_per_width_m2_s(),
        base.flow_per_width_m2_s(),
        3e-8,
    );
    near(
        sheared.surface_velocity_m_s(),
        base.surface_velocity_m_s(),
        base.surface_velocity_m_s(),
        3e-8,
    );
    near(
        sheared.power().viscous_dissipation_w_m2,
        base.dissipation_w_m2(),
        base.dissipation_w_m2(),
        3e-8,
    );
    assert_eq!(sheared.surface_shear_residual_pa(), 0.0);
    let driven = base.with_surface_shear(0.2, controls()).unwrap();
    assert!(driven.flow_per_width_m2_s() > sheared.flow_per_width_m2_s());
    assert_eq!(driven.thickness_m(), base.thickness_m());
}

#[test]
fn pure_surface_shear_retains_thermal_yasuda_source_and_differential_mobility() {
    let source = law([0.08, 0.006, 0.4, 0.45, 0.7], [0.0, 1e6]);
    for traction in [-0.7, 0.7] {
        let spec = CarreauShearedFilmSpec {
            density_kg_m3: 950.0,
            thickness_m: 0.003,
            tangential_acceleration_m_s2: 0.0,
            surface_shear_pa: traction,
        };
        let actual = CarreauShearedFilm::new(source, 340.0, spec, 1e4, controls()).unwrap();
        let sample = actual.sample(0.4 * spec.thickness_m).unwrap();
        let state = source
            .evaluate(340.0, sample.velocity_gradient_per_s.abs())
            .unwrap();
        let rate = sample.velocity_gradient_per_s;
        let h = spec.thickness_m;
        assert_eq!(actual.law(), &source);
        assert_eq!(actual.temperature_k(), 340.0);
        near(
            sample.dynamic_viscosity_pa_s * rate,
            traction,
            traction.abs(),
            2e-11,
        );
        near(
            actual.surface_velocity_m_s(),
            h * rate,
            (h * rate).abs(),
            2e-11,
        );
        near(
            actual.flow_per_width_m2_s(),
            0.5 * h * h * rate,
            (h * h * rate).abs(),
            2e-11,
        );
        near(
            actual.power().viscous_dissipation_w_m2,
            traction * h * rate,
            traction * h * rate,
            2e-11,
        );
        near(
            actual.surface_velocity_traction_derivative_m_s_pa(),
            h / state.stress_tangent_pa_s(),
            h / state.stress_tangent_pa_s(),
            2e-11,
        );
        assert!(state.stress_tangent_pa_s() < state.dynamic_viscosity_pa_s());
        assert_eq!(actual.power().body_work_w_m2, 0.0);
        assert_eq!(sample.velocity_curvature_per_m_s, 0.0);
    }
}

#[test]
fn mixed_slit_traction_keeps_lower_wall_translation_and_original_heat() {
    let source = law(CURVE, [0.0, 1e6]);
    let (width, gap, length, pressure, traction) = (0.03, 0.002, 0.2, 380.0, -0.8);
    let slit = LaminarCarreauSlit::new(source, 300.0, width, gap, length, 1e4).unwrap();
    let stationary = slit
        .profile_at_pressure_and_upper_traction(pressure, 0.0, traction, controls())
        .unwrap();
    let speed = 0.37;
    let moving = slit
        .profile_at_pressure_and_upper_traction(pressure, speed, traction, controls())
        .unwrap();
    assert_eq!(moving.flow().wall_stress_iterations, 0);
    assert_eq!(stationary.flow().wall_stress_iterations, 0);
    near(
        moving.sample(0.0).unwrap().axial_velocity_m_s,
        speed,
        speed,
        2e-12,
    );
    near(
        moving.flow().wall_tractions_pa[1],
        traction,
        traction.abs(),
        2e-12,
    );
    near(
        moving.flow().volumetric_flow_m3_s - stationary.flow().volumetric_flow_m3_s,
        width * gap * speed,
        width * gap * speed,
        2e-11,
    );
    near(
        moving.flow().dissipated_power_w,
        stationary.flow().dissipated_power_w,
        stationary.flow().dissipated_power_w,
        2e-11,
    );
    let total_input = moving.flow().pressure_power_w + moving.flow().wall_power_w;
    near(
        total_input - moving.flow().dissipated_power_w,
        moving.flow().work_residual_w,
        moving.flow().dissipated_power_w,
        2e-10,
    );
    let mut velocity_controls = controls();
    velocity_controls.maximum_wall_stress_iterations = 96;
    let equivalent = slit
        .profile_at_pressure_and_walls(
            pressure,
            moving.flow().accepted_wall_velocities_m_s,
            velocity_controls,
        )
        .unwrap();
    near(
        equivalent.flow().wall_tractions_pa[1],
        traction,
        traction.abs(),
        3e-8,
    );
    near(
        equivalent.flow().volumetric_flow_m3_s,
        moving.flow().volumetric_flow_m3_s,
        moving.flow().volumetric_flow_m3_s.abs(),
        3e-8,
    );
}

#[test]
fn unsupported_sources_unrepresentable_traction_and_cancellation_refuse() {
    let source = law(CURVE, [0.0, 1e6]);
    let spec = CarreauShearedFilmSpec {
        density_kg_m3: 950.0,
        thickness_m: 0.002,
        tangential_acceleration_m_s2: 2.0,
        surface_shear_pa: -0.8,
    };
    assert!(CarreauShearedFilm::new(source, 270.0, spec, 1e4, controls()).is_err());
    let bounded = law(CURVE, [0.0, 1.0]);
    assert!(CarreauShearedFilm::new(bounded, 300.0, spec, 1.0, controls()).is_err());
    assert!(
        CarreauShearedFilm::new(
            source,
            300.0,
            CarreauShearedFilmSpec {
                surface_shear_pa: 1e-30,
                ..spec
            },
            1e4,
            controls()
        )
        .is_err()
    );
    for invalid in [
        CarreauShearedFilmSpec {
            density_kg_m3: 0.0,
            ..spec
        },
        CarreauShearedFilmSpec {
            thickness_m: -1.0,
            ..spec
        },
        CarreauShearedFilmSpec {
            surface_shear_pa: f64::NAN,
            ..spec
        },
    ] {
        assert!(CarreauShearedFilm::new(source, 300.0, invalid, 1e4, controls()).is_err());
    }
    let mut checks = 0;
    let cancelled = CarreauShearedFilm::with_cancel(source, 300.0, spec, 1e4, controls(), || {
        checks += 1;
        checks == 4
    });
    assert!(matches!(cancelled, Err(CarreauFilmError::Cancelled)));
    assert_eq!(checks, 4);
    let accepted = film(spec);
    assert!(matches!(
        accepted.sample_with_cancel(0.001, || true),
        Err(CarreauFilmError::Cancelled)
    ));
    assert!(accepted.sample(-0.001).is_err());
    assert!(accepted.sample(0.003).is_err());
}
