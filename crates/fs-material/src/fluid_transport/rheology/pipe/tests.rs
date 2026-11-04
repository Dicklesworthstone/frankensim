//! Synthetic mathematical controls; these do not calibrate a material or a
//! laminar-transition criterion. Radial energy integration is independent of
//! the production shear-coordinate volume-flow quadrature.

use super::*;
use crate::fluid_transport::rheology::CarreauYasuda;
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let scale = actual.abs().max(expected.abs());
    if scale != 0.0 {
        assert!(
            (actual / scale - expected / scale).abs() <= tolerance,
            "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
        );
    }
}

fn law(parameters: [f64; 5], support: [f64; 2], mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    let [zero, infinity, time, n, a] = parameters;
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(zero, infinity, time, n, a, support).unwrap(),
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

fn nonlinear() -> LaminarCarreauPipe {
    LaminarCarreauPipe::new(
        law(
            [12.0, 0.2, 0.6, 0.45, 2.0],
            [0.0, 1e8],
            ViscosityShiftMode::Vertical,
        ),
        300.0,
        0.01,
        0.5,
        1e6,
    )
    .unwrap()
}

#[test]
fn g0_newtonian_limits_recover_poiseuille_and_the_parabolic_thermal_profile() {
    for parameters in [
        [3.0, 0.2, 0.0, 0.4, 0.6],
        [3.0, 0.2, 0.6, 1.0, 0.6],
        [3.0, 3.0, 0.6, 0.4, 0.6],
    ] {
        for mode in [
            ViscosityShiftMode::Vertical,
            ViscosityShiftMode::VerticalAndRate,
        ] {
            let source = law(parameters, [0.0, 1e8], mode);
            let pipe = LaminarCarreauPipe::new(source, 320.0, 0.012, 0.7, 1e6).unwrap();
            let mu = source
                .evaluate(320.0, 0.0)
                .unwrap()
                .dynamic_viscosity_pa_s();
            for dp in [-1700.0, 1700.0] {
                let profile = pipe
                    .profile_at_pressure(dp, CarreauPipeControls::default())
                    .unwrap();
                assert_eq!(profile.pipe().law(), &source);
                assert_eq!(profile.pipe().temperature_k(), 320.0);
                close(
                    profile.flow().volumetric_flow_m3_s,
                    core::f64::consts::PI * pipe.radius_m().powi(4) * dp
                        / (8.0 * mu * pipe.length_m()),
                    2e-14,
                );
                close(profile.flow().pressure_drop_pa, dp, 2e-14);
                assert_eq!(profile.flow().quadrature_relative_change, 0.0);
                for fraction in [0.0, 0.2, 0.7, 1.0] {
                    let r = fraction * pipe.radius_m();
                    let sample = profile.sample(r).unwrap();
                    let gradient = -dp * r / (2.0 * mu * pipe.length_m());
                    close(
                        sample.axial_velocity_m_s,
                        dp * (pipe.radius_m().powi(2) - r * r) / (4.0 * mu * pipe.length_m()),
                        3e-14,
                    );
                    close(sample.radial_velocity_gradient_per_s, gradient, 3e-14);
                    close(sample.radial_shear_stress_pa, mu * gradient, 3e-14);
                    close(sample.dissipation_w_m3, mu * gradient * gradient, 4e-14);
                }
            }
        }
    }
}

#[test]
fn g1_radial_velocity_and_heat_integrate_to_the_same_accepted_pressure_work() {
    let pipe = nonlinear();
    let controls = CarreauPipeControls {
        relative_quadrature_tolerance: 1e-10,
        ..CarreauPipeControls::default()
    };
    let profile = pipe.profile_at_pressure(2000.0, controls).unwrap();
    let flow = profile.flow();
    // Integrate independently in physical radius, using sampled local u and
    // tau*gamma. The production Q uses a different variable and its stress tangent.
    let panels = 256;
    let mut volume_integral = 0.0;
    let mut heat_integral = 0.0;
    for i in 0..=panels {
        let x = f64::from(i) / f64::from(panels);
        let sample = profile.sample(x * pipe.radius_m()).unwrap();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        volume_integral += weight * sample.axial_velocity_m_s * x;
        heat_integral += weight * sample.dissipation_w_m3 * x;
        close(
            sample.radial_shear_stress_pa,
            -flow.pressure_drop_pa * x * pipe.radius_m() / (2.0 * pipe.length_m()),
            3e-14,
        );
        close(
            sample.radial_shear_stress_pa,
            sample.dynamic_viscosity_pa_s * sample.radial_velocity_gradient_per_s,
            2e-11,
        );
    }
    let area_factor =
        2.0 * core::f64::consts::PI * pipe.radius_m().powi(2) / (3.0 * f64::from(panels));
    close(
        area_factor * volume_integral,
        flow.volumetric_flow_m3_s,
        2e-7,
    );
    close(
        area_factor * pipe.length_m() * heat_integral,
        flow.dissipated_power_w,
        2e-7,
    );
    assert!(
        (flow.pressure_residual_pa / flow.requested_pressure_drop_pa).abs()
            <= controls.relative_constitutive_tolerance
    );
    // A bulk apparent-viscosity Poiseuille shortcut would overpredict this
    // flow: the higher-viscosity inner fluid is actually resolved.
    let wall_mu = profile
        .sample(pipe.radius_m())
        .unwrap()
        .dynamic_viscosity_pa_s;
    let wall_viscosity_shortcut =
        core::f64::consts::PI * pipe.radius_m().powi(4) * flow.pressure_drop_pa
            / (8.0 * wall_mu * pipe.length_m());
    assert!(flow.volumetric_flow_m3_s < 0.95 * wall_viscosity_shortcut);
    for x in [0.2, 0.6, 0.9] {
        let r = x * pipe.radius_m();
        let h = 1e-4 * pipe.radius_m();
        let derivative = (profile.sample(r + h).unwrap().axial_velocity_m_s
            - profile.sample(r - h).unwrap().axial_velocity_m_s)
            / (2.0 * h);
        close(
            derivative,
            profile.sample(r).unwrap().radial_velocity_gradient_per_s,
            2e-5,
        );
    }
}

#[test]
fn g1_pressure_reversal_and_the_rest_limit_preserve_passive_local_heating() {
    let pipe = nonlinear();
    let positive = pipe
        .profile_at_pressure(2000.0, CarreauPipeControls::default())
        .unwrap();
    let negative = pipe
        .profile_at_pressure(-2000.0, CarreauPipeControls::default())
        .unwrap();
    assert_eq!(
        positive.flow().volumetric_flow_m3_s,
        -negative.flow().volumetric_flow_m3_s
    );
    assert_eq!(
        positive.flow().dissipated_power_w,
        negative.flow().dissipated_power_w
    );
    for x in [0.0, 0.3, 0.8, 1.0] {
        let p = positive.sample(x * pipe.radius_m()).unwrap();
        let n = negative.sample(x * pipe.radius_m()).unwrap();
        assert_eq!(p.axial_velocity_m_s, -n.axial_velocity_m_s);
        assert_eq!(
            p.radial_velocity_gradient_per_s,
            -n.radial_velocity_gradient_per_s
        );
        assert_eq!(p.radial_shear_stress_pa, -n.radial_shear_stress_pa);
        assert_eq!(p.dynamic_viscosity_pa_s, n.dynamic_viscosity_pa_s);
        assert_eq!(p.dissipation_w_m3, n.dissipation_w_m3);
        assert!(p.dissipation_w_m3 >= 0.0);
    }
    let rest = pipe
        .profile_at_pressure(0.0, CarreauPipeControls::default())
        .unwrap();
    assert_eq!(rest.flow().volumetric_flow_m3_s, 0.0);
    assert_eq!(rest.flow().dissipated_power_w, 0.0);
    for x in [0.0, 0.5, 1.0] {
        let sample = rest.sample(x * pipe.radius_m()).unwrap();
        assert_eq!(sample.axial_velocity_m_s, 0.0);
        assert_eq!(sample.dissipation_w_m3, 0.0);
        assert_eq!(sample.radial_velocity_gradient_per_s, 0.0);
        assert_eq!(sample.dynamic_viscosity_pa_s, 12.0);
    }
}

#[test]
fn g1_low_and_high_rate_plateaus_approach_their_distinct_newtonian_flows() {
    let pipe = LaminarCarreauPipe::new(
        law(
            [8.0, 0.25, 0.7, 0.35, 2.0],
            [0.0, 1e15],
            ViscosityShiftMode::Vertical,
        ),
        300.0,
        0.02,
        1.0,
        1e15,
    )
    .unwrap();
    // These broad synthetic pressures check mathematical asymptotes only.
    for (pressure, plateau, tolerance) in [(1e-6, 8.0, 1e-9), (1e14, 0.25, 2e-6)] {
        let profile = pipe
            .profile_at_pressure(pressure, CarreauPipeControls::default())
            .unwrap();
        let expected = core::f64::consts::PI * pipe.radius_m().powi(4) * pressure
            / (8.0 * plateau * pipe.length_m());
        close(profile.flow().volumetric_flow_m3_s, expected, tolerance);
    }
}

#[test]
fn g0_singular_scalar_viscosity_slope_has_a_regular_pipe_core_and_near_wall_limit() {
    let source = law(
        [12.0, 0.2, 0.6, 0.45, 0.6],
        [0.0, 1e8],
        ViscosityShiftMode::Vertical,
    );
    assert_eq!(
        source
            .evaluate(300.0, 0.0)
            .unwrap()
            .shear_derivative_pa_s2(),
        None
    );
    let pipe = LaminarCarreauPipe::new(source, 300.0, 0.01, 0.5, 1e6).unwrap();
    let profile = pipe
        .profile_at_pressure(2000.0, CarreauPipeControls::default())
        .unwrap();
    let core = profile.sample(0.0).unwrap();
    assert!(core.axial_velocity_m_s > 0.0);
    assert_eq!(core.radial_velocity_gradient_per_s, 0.0);
    assert_eq!(core.dissipation_w_m3, 0.0);
    let near = pipe.radius_m().next_down();
    let sample = profile.sample(near).unwrap();
    assert!(sample.axial_velocity_m_s > 0.0);
    close(
        sample.axial_velocity_m_s,
        profile.flow().wall_shear_rate_per_s * (pipe.radius_m() - near),
        3e-10,
    );
    assert_eq!(
        profile.sample(pipe.radius_m()).unwrap().axial_velocity_m_s,
        0.0
    );
}

#[test]
fn g3_temperature_shift_and_the_entire_pipe_rate_domain_remain_source_owned() {
    let curve = CarreauYasuda::new(10.0, 0.1, 1.0, 0.5, 2.0, [0.0, 10.0]).unwrap();
    let shift = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 1.0,
            c2_k: 100.0,
        },
        300.0,
        [300.0, 400.0],
    )
    .unwrap();
    let vertical = ThermalCarreauYasuda::new(curve, shift, ViscosityShiftMode::Vertical);
    let combined = ThermalCarreauYasuda::new(curve, shift, ViscosityShiftMode::VerticalAndRate);
    assert!(LaminarCarreauPipe::new(vertical, 400.0, 1.0, 0.5, 20.0).is_err());
    let pipe = LaminarCarreauPipe::new(combined, 400.0, 1.0, 0.5, 20.0).unwrap();
    let ceiling_stress = combined
        .evaluate(400.0, 20.0)
        .unwrap()
        .dynamic_viscosity_pa_s()
        * 20.0;
    let profile = pipe
        .profile_at_pressure(0.5 * ceiling_stress, CarreauPipeControls::default())
        .unwrap();
    let wall = profile.sample(1.0).unwrap();
    close(
        wall.dynamic_viscosity_pa_s,
        combined
            .evaluate(400.0, wall.radial_velocity_gradient_per_s.abs())
            .unwrap()
            .dynamic_viscosity_pa_s(),
        2e-14,
    );
    assert!(
        pipe.profile_at_pressure(1.0001 * ceiling_stress, CarreauPipeControls::default())
            .is_err()
    );
    assert!(LaminarCarreauPipe::new(combined, 401.0, 1.0, 0.5, 20.0).is_err());
    let no_core = law(
        [10.0, 0.1, 1.0, 0.5, 2.0],
        [0.1, 100.0],
        ViscosityShiftMode::Vertical,
    );
    assert!(LaminarCarreauPipe::new(no_core, 300.0, 1.0, 0.5, 20.0).is_err());
}

#[test]
fn g0_invalid_controls_and_exhausted_inverse_or_quadrature_publish_no_profile() {
    let pipe = nonlinear();
    let inverse = pipe
        .profile_at_pressure(
            2000.0,
            CarreauPipeControls {
                maximum_constitutive_iterations: 1,
                ..CarreauPipeControls::default()
            },
        )
        .unwrap_err();
    assert!(
        inverse
            .to_string()
            .contains("constitutive inverse exhausted")
    );
    let quadrature = pipe
        .profile_at_pressure(
            2000.0,
            CarreauPipeControls {
                maximum_quadrature_refinements: 1,
                ..CarreauPipeControls::default()
            },
        )
        .unwrap_err();
    assert!(quadrature.to_string().contains("quadrature exhausted"));
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            pipe.profile_at_pressure(
                2000.0,
                CarreauPipeControls {
                    relative_constitutive_tolerance: bad,
                    ..CarreauPipeControls::default()
                }
            )
            .is_err()
        );
        assert!(LaminarCarreauPipe::new(*pipe.law(), 300.0, bad, 1.0, 10.0).is_err());
    }
    assert!(
        pipe.profile_at_pressure(f64::INFINITY, CarreauPipeControls::default())
            .is_err()
    );
    let profile = pipe
        .profile_at_pressure(2000.0, CarreauPipeControls::default())
        .unwrap();
    for bad in [-1.0, pipe.radius_m().next_up(), f64::INFINITY, f64::NAN] {
        assert!(profile.sample(bad).is_err());
    }
}

#[test]
fn g0_representable_flow_survives_extreme_intermediate_geometry_products() {
    for (radius, viscosity, ceiling, expected) in [
        (1e-150, 1e-200, 1e202, 0.25 * core::f64::consts::PI * 1e-250),
        (1e150, 1e200, 1e-198, 0.25 * core::f64::consts::PI * 1e250),
    ] {
        let source = law(
            [viscosity, viscosity, 1.0, 0.5, 2.0],
            [0.0, ceiling],
            ViscosityShiftMode::Vertical,
        );
        let pipe = LaminarCarreauPipe::new(source, 300.0, radius, radius, ceiling).unwrap();
        let profile = pipe
            .profile_at_pressure(2.0, CarreauPipeControls::default())
            .unwrap();
        close(profile.flow().volumetric_flow_m3_s, expected, 5e-12);
        close(profile.flow().dissipated_power_w, 2.0 * expected, 5e-12);
        assert!(profile.sample(0.5 * radius).unwrap().axial_velocity_m_s > 0.0);
    }
}
