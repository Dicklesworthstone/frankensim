use super::*;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

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
fn spec() -> CarreauFilmSpec {
    CarreauFilmSpec {
        density_kg_m3: 1000.0,
        surface_tension_n_m: 0.05,
        flow_per_width_m2_s: 1.0e-5,
        tangential_acceleration_m_s2: 3.0,
        normal_acceleration_m_s2: 3.0,
    }
}
fn controls() -> CarreauFilmControls {
    let mut controls = CarreauFilmControls::default();
    controls.slit.local.relative_constitutive_tolerance = 1.0e-12;
    controls.slit.local.relative_quadrature_tolerance = 1.0e-10;
    controls.relative_flow_tolerance = 1.0e-10;
    controls
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0e-25),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn constant_source_at_actual_temperature_recovers_the_native_nusselt_profile() {
    let source = law([0.02, 0.001, 0.5, 1.0, 2.0], [0.0, 1.0e5]);
    let temperature = 330.0;
    let mu = source
        .evaluate(temperature, 0.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let data = spec();
    let film = CarreauFilm::at_flow(source, temperature, data, 1.0e5, controls()).unwrap();
    let newtonian = NusseltFilm::new(
        FilmMaterial {
            density_kg_m3: data.density_kg_m3,
            kinematic_viscosity_m2_s: mu / data.density_kg_m3,
            surface_tension_n_m: data.surface_tension_n_m,
        },
        data.flow_per_width_m2_s,
        data.tangential_acceleration_m_s2,
        data.normal_acceleration_m_s2,
    )
    .unwrap();
    near(film.thickness_m(), newtonian.thickness_m(), 2.0e-12);
    near(
        film.surface_velocity_m_s(),
        newtonian.surface_velocity_m_s(),
        2.0e-12,
    );
    assert_eq!(film.law(), &source);
    assert_eq!(film.temperature_k(), temperature);
    for fraction in [0.0, 0.2, 0.7, 1.0] {
        let y = fraction * film.thickness_m();
        let sample = film.sample(y).unwrap();
        near(
            sample.velocity_m_s,
            newtonian.velocity_m_s(y).unwrap(),
            2.0e-11,
        );
        near(
            sample.shear_stress_pa,
            newtonian.shear_pa(y).unwrap(),
            2.0e-11,
        );
        assert_eq!(sample.dynamic_viscosity_pa_s, mu);
        assert_eq!(sample.stress_tangent_pa_s, mu);
    }
}

#[test]
fn shear_thinning_thickness_profile_and_heat_match_independent_stress_integration() {
    let source = law([0.02, 0.001, 0.5, 0.5, 2.0], [0.0, 1.0e5]);
    let data = spec();
    let film = CarreauFilm::at_flow(source, 300.0, data, 1.0e5, controls()).unwrap();
    // Independent scalar stress inversion + adaptive Gauss integration:
    // q(h)=integral_0^h d*gamma(rho*g*d) dd; u_s=integral_0^h gamma(rho*g*d) dd.
    near(film.thickness_m(), 0.000_301_515_480_792_150_53, 2.0e-8);
    near(
        film.surface_velocity_m_s(),
        0.046_293_761_532_898_13,
        2.0e-8,
    );
    near(
        film.sample(0.0).unwrap().velocity_gradient_per_s,
        380.447_727_008_235_63,
        2.0e-8,
    );
    near(
        film.flow_per_width_m2_s(),
        data.flow_per_width_m2_s,
        2.0e-10,
    );
    near(film.dissipation_w_m2(), 0.03, 2.0e-8);
    assert!(film.thickness_m() < 0.000_584_803_547_642_573_5);
    let wall = film.sample(0.0).unwrap();
    assert_eq!(wall.velocity_m_s, 0.0);
    assert!(wall.stress_tangent_pa_s < wall.dynamic_viscosity_pa_s);
    for fraction in [0.0, 0.5, 1.0] {
        let y = fraction * film.thickness_m();
        let sample = film.sample(y).unwrap();
        near(
            sample.shear_stress_pa,
            data.density_kg_m3 * data.tangential_acceleration_m_s2 * (film.thickness_m() - y),
            2.0e-12,
        );
        near(
            sample.stress_tangent_pa_s * sample.velocity_curvature_per_m_s,
            -data.density_kg_m3 * data.tangential_acceleration_m_s2,
            2.0e-12,
        );
    }
    let surface = film.sample(film.thickness_m()).unwrap();
    assert_eq!(surface.velocity_gradient_per_s, 0.0);
    assert_eq!(surface.shear_stress_pa, 0.0);
    assert_eq!(surface.dynamic_viscosity_pa_s, 0.02);
    assert_eq!(surface.stress_tangent_pa_s, 0.02);
}

#[test]
fn original_support_work_caps_and_cancellation_refuse_without_a_partial_film() {
    let source = law([0.02, 0.001, 0.5, 0.5, 2.0], [0.0, 1.0e5]);
    assert!(matches!(
        CarreauFilm::at_flow(source, 300.0, spec(), 5.0, controls()),
        Err(CarreauFilmError::OutsideRateSupport)
    ));
    let mut insufficient = controls();
    insufficient.maximum_thickness_iterations = 1;
    assert!(matches!(
        CarreauFilm::at_flow(source, 300.0, spec(), 1.0e5, insufficient),
        Err(CarreauFilmError::UnresolvedFlow { iterations: 1, .. })
    ));
    let unsupported = law([0.02, 0.001, 0.5, 0.5, 2.0], [1.0, 1.0e5]);
    assert!(matches!(
        CarreauFilm::at_flow(unsupported, 300.0, spec(), 1.0e5, controls()),
        Err(CarreauFilmError::Material(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        CarreauFilm::at_flow_with_cancel(source, 300.0, spec(), 1.0e5, controls(), || {
            polls += 1;
            polls == 12
        }),
        Err(CarreauFilmError::Cancelled)
    ));
}
