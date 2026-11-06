use super::*;
use crate::film::NusseltFilm;

const LIQUID: FilmMaterial = FilmMaterial {
    density_kg_m3: 2.0,
    kinematic_viscosity_m2_s: 0.5,
    surface_tension_n_m: 0.03,
};

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn signed_profiles_conserve_flux_momentum_and_mechanical_power() {
    let mu = LIQUID.density_kg_m3 * LIQUID.kinematic_viscosity_m2_s;
    for h in [0.125, 0.5] {
        for g in [-6.0, 0.0, 6.0] {
            for surface in [-8.0, -4.0, 0.0, 3.0] {
                let film = ShearedFilm::new(LIQUID, h, g, surface).unwrap();
                let mut flux = 0.0;
                let mut dissipation = 0.0;
                // Simpson is exact for both the quadratic velocity and tau².
                for (z, weight) in [(0.0, 1.0), (0.5, 4.0), (1.0, 1.0)] {
                    let y = z * h;
                    let reference_velocity = g * (h * y - y * y / 2.0)
                        / LIQUID.kinematic_viscosity_m2_s
                        + surface * y / mu;
                    let reference_stress = surface + LIQUID.density_kg_m3 * g * (h - y);
                    let u = film.velocity_m_s(y).unwrap();
                    let tau = film.shear_pa(y).unwrap();
                    near(u, reference_velocity, 1.0);
                    near(tau, reference_stress, 1.0);
                    flux += weight * u * h / 6.0;
                    dissipation += weight * tau * tau * h / (6.0 * mu);
                }
                near(film.flow_per_width_m2_s(), flux, 1.0);
                near(film.mean_velocity_m_s(), flux / h, 1.0);
                near(
                    film.wall_shear_pa() - film.surface_shear_pa(),
                    LIQUID.density_kg_m3 * g * h,
                    1.0,
                );
                assert_eq!(film.velocity_m_s(0.0).unwrap(), 0.0);
                assert_eq!(film.shear_pa(h).unwrap(), surface);
                let power = film.power();
                near(power.viscous_dissipation_w_m2, dissipation, 1.0);
                near(power.body_work_w_m2, LIQUID.density_kg_m3 * g * flux, 1.0);
                near(
                    power.surface_work_w_m2,
                    surface * film.velocity_m_s(h).unwrap(),
                    1.0,
                );
                near(
                    power.body_work_w_m2 + power.surface_work_w_m2,
                    dissipation,
                    1.0,
                );
                assert!(power.viscous_dissipation_w_m2 >= 0.0);
            }
        }
    }
}

#[test]
fn opposing_surface_stress_retains_counterflow_and_positive_heat_at_zero_net_flux() {
    let film = ShearedFilm::new(LIQUID, 0.5, 6.0, -4.0).unwrap();
    assert_eq!(film.flow_per_width_m2_s(), 0.0);
    assert!(film.velocity_m_s(0.125).unwrap() > 0.0);
    assert_eq!(film.surface_velocity_m_s(), -0.5);
    assert_eq!(film.wall_shear_pa(), 2.0);
    assert_eq!(film.power().body_work_w_m2, 0.0);
    near(film.power().surface_work_w_m2, 2.0, 1.0);
    near(film.power().viscous_dissipation_w_m2, 2.0, 1.0);

    let reversed = ShearedFilm::new(LIQUID, 0.5, 6.0, -6.0).unwrap();
    assert!(reversed.flow_per_width_m2_s() < 0.0);
    assert!(reversed.power().body_work_w_m2 < 0.0);
    assert!(reversed.power().surface_work_w_m2 > 0.0);
    assert!(reversed.power().viscous_dissipation_w_m2 > 0.0);
}

#[test]
fn zero_surface_shear_recovers_native_nusselt_and_horizontal_shear_gives_couette() {
    let original = NusseltFilm::new(LIQUID, 0.2, 6.0, 0.0).unwrap();
    let film = original.with_surface_shear(0.0).unwrap();
    assert_eq!(film.thickness_m(), original.thickness_m());
    near(film.flow_per_width_m2_s(), 0.2, 1.0);
    near(film.mean_velocity_m_s(), original.mean_velocity_m_s(), 1.0);
    near(
        film.surface_velocity_m_s(),
        original.surface_velocity_m_s(),
        1.0,
    );
    near(film.wall_shear_pa(), original.wall_shear_pa(), 1.0);
    for z in [0.0, 0.25, 0.5, 1.0] {
        let y = z * film.thickness_m();
        near(
            film.velocity_m_s(y).unwrap(),
            original.velocity_m_s(y).unwrap(),
            1.0,
        );
        near(
            film.shear_pa(y).unwrap(),
            original.shear_pa(y).unwrap(),
            1.0,
        );
    }
    let driven = original.with_surface_shear(3.0).unwrap();
    assert!(driven.flow_per_width_m2_s() > film.flow_per_width_m2_s());

    let couette = ShearedFilm::new(LIQUID, 0.5, 0.0, 2.0).unwrap();
    assert_eq!(couette.flow_per_width_m2_s(), 0.25);
    assert_eq!(couette.velocity_m_s(0.25).unwrap(), 0.5);
    assert_eq!(couette.surface_velocity_m_s(), 1.0);
    near(couette.power().viscous_dissipation_w_m2, 2.0, 1.0);
    let rest = ShearedFilm::new(LIQUID, 0.5, 0.0, 0.0).unwrap();
    assert_eq!(rest.power().viscous_dissipation_w_m2, 0.0);
}

#[test]
fn scaled_shear_power_does_not_square_away_representable_loss() {
    for (rho, nu, h, stress) in [
        (1e200, 1e-100, 1e-100, 1e200),
        (1e-200, 1e100, 1e100, 1e-200),
    ] {
        let liquid = FilmMaterial {
            density_kg_m3: rho,
            kinematic_viscosity_m2_s: nu,
            surface_tension_n_m: 0.0,
        };
        let film = ShearedFilm::new(liquid, h, 0.0, stress).unwrap();
        near(film.surface_velocity_m_s(), 1.0, 1.0);
        near(film.power().viscous_dissipation_w_m2 / stress, 1.0, 1.0);
        near(film.power().surface_work_w_m2 / stress, 1.0, 1.0);
    }
}

#[test]
fn invalid_geometry_sources_and_unrepresentable_coefficients_refuse() {
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(ShearedFilm::new(LIQUID, h, 1.0, 0.0).is_err());
    }
    assert!(ShearedFilm::new(LIQUID, 1.0, f64::NAN, 0.0).is_err());
    assert!(ShearedFilm::new(LIQUID, 1.0, 1.0, f64::INFINITY).is_err());
    assert!(
        ShearedFilm::new(
            FilmMaterial {
                kinematic_viscosity_m2_s: 0.0,
                ..LIQUID
            },
            1.0,
            1.0,
            0.0
        )
        .is_err()
    );
    assert_eq!(
        ShearedFilm::new(LIQUID, 1e300, 1.0, 0.0),
        Err(FilmError::NonRepresentable)
    );
    let film = ShearedFilm::new(LIQUID, 0.5, 1.0, 0.0).unwrap();
    for y in [-1.0, 0.6, f64::NAN, f64::INFINITY] {
        assert!(film.velocity_m_s(y).is_err());
        assert!(film.shear_pa(y).is_err());
    }
}
