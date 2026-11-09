//! Independent spherical geometry and isothermal d-squared refinement.
use super::*;
use crate::evaporation::{
    EvaporationProperties, InterfaceKinetics, KineticEvaporation, MOLAR_GAS_CONSTANT,
};

const D: f64 = 2e-5;
const RHO: f64 = 1000.0;
const CP: f64 = 4000.0;
const LATENT: f64 = 2.4e6;

struct Phase {
    expansion_per_k: f64,
}
impl EvaporationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        Ok(EvaporationProperties {
            saturation_pressure_pa: 4000.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: LATENT,
            liquid_density_kg_m3: RHO / (1.0 + self.expansion_per_k * (t - 300.0)),
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn droplet() -> SphericalDropletState {
    SphericalDropletState {
        mass_kg: (4.0 * core::f64::consts::PI / 3.0) * RHO * 1e-9,
        specific_enthalpy_j_kg: 0.0,
    }
}
fn temperature(h: f64) -> Result<f64, EvaporationError> {
    Ok(300.0 + h / CP)
}
fn film() -> SphericalGasFilm {
    SphericalGasFilm {
        temperature_k: 300.0,
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa: 1000.0,
        binary_diffusivity_m2_s: D,
    }
}
fn controls() -> SolverConfig {
    SolverConfig {
        specific_enthalpy_bracket_j_kg: [-80_000.0, 120_000.0],
        relative_tolerance: 1e-12,
        absolute_energy_tolerance_j: 1e-11,
        max_iterations: 128,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn sphere_geometry_and_stagnant_diffusion_follow_independent_volume_and_radius() {
    let sphere = SphericalDropletGeometry::from_mass(droplet().mass_kg, RHO).unwrap();
    near(sphere.radius_m(), 1e-3, 2e-17);
    near(sphere.diameter_m(), 2e-3, 4e-17);
    near(sphere.area_m2(), 4.0 * core::f64::consts::PI * 1e-6, 1e-18);
    near(
        sphere.volume_m3(),
        4.0 * core::f64::consts::PI / 3.0 * 1e-9,
        1e-22,
    );
    near(
        sphere.mass_transfer_coefficient_m_s(D).unwrap(),
        0.02,
        1e-14,
    );
    assert_eq!(sphere.mass_transfer_coefficient_m_s(0.0).unwrap(), 0.0);
    for invalid in [-1.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(SphericalDropletGeometry::from_mass(invalid, RHO).is_err());
        assert!(SphericalDropletGeometry::from_mass(droplet().mass_kg, invalid).is_err());
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(sphere.mass_transfer_coefficient_m_s(invalid).is_err());
    }
}

#[test]
fn supplied_latent_heat_keeps_real_calorics_isothermal_and_d_squared_converges() {
    let law = Phase {
        expansion_per_k: 0.0,
    };
    // Independent integrated binary diffusion, with fixed far-field coordinates.
    // Supplying its latent power holds the actual h-to-T chart at 300 K.
    let concentration =
        0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0) * (99_000.0_f64 / 96_000.0).ln();
    let horizon = 100.0;
    let exact_d_squared = 4e-6 - 8.0 * D * concentration * horizon / RHO;
    let errors: Vec<_> = [8, 16, 32]
        .into_iter()
        .map(|count| {
            let mut state = droplet();
            let dt = horizon / f64::from(count);
            for _ in 0..count {
                let sphere = state.geometry(&law, 300.0).unwrap();
                let independent_mass_rate =
                    4.0 * core::f64::consts::PI * sphere.radius_m() * D * concentration;
                let accepted = advance_droplet_with_heat(
                    &law,
                    film(),
                    state,
                    dt,
                    independent_mass_rate * dt * LATENT,
                    controls(),
                    temperature,
                    || Ok(()),
                )
                .unwrap();
                assert_eq!(accepted.phase.surface_temperature_k, 300.0);
                assert_eq!(accepted.droplet.specific_enthalpy_j_kg, 0.0);
                assert!(accepted.geometry[1].area_m2() < accepted.geometry[0].area_m2());
                assert!(accepted.phase.energy_residual_j.abs() <= 1e-11);
                near(
                    state.mass_kg - accepted.droplet.mass_kg,
                    accepted.phase.transferred_mass_kg,
                    0.0,
                );
                state = accepted.droplet;
            }
            let diameter = state.geometry(&law, 300.0).unwrap().diameter_m();
            (diameter * diameter - exact_d_squared).abs()
        })
        .collect();
    assert!(
        errors[0] > 1e-11,
        "inactive geometry refinement: {errors:?}"
    );
    assert!(
        errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn heat_changes_radius_through_original_density_even_when_diffusion_is_sealed() {
    let law = Phase {
        expansion_per_k: 0.001,
    };
    let initial = droplet();
    let sealed = SphericalGasFilm {
        binary_diffusivity_m2_s: 0.0,
        ..film()
    };
    let idle = advance_droplet_with_heat(
        &law,
        sealed,
        initial,
        1.0,
        0.0,
        controls(),
        temperature,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(idle.droplet, initial);
    assert_eq!(idle.geometry[0], idle.geometry[1]);
    let heat = initial.mass_kg * CP * 5.0;
    let warm = advance_droplet_with_heat(
        &law,
        sealed,
        initial,
        1.0,
        heat,
        controls(),
        temperature,
        || Ok(()),
    )
    .unwrap();
    near(warm.phase.surface_temperature_k, 305.0, 1e-8);
    assert_eq!(warm.droplet.mass_kg, initial.mass_kg);
    assert_eq!(warm.phase.transferred_mass_kg, 0.0);
    assert_eq!(warm.phase.recession_m, 0.0);
    assert!(warm.geometry[1].radius_m() > warm.geometry[0].radius_m());
    let density = RHO / (1.0 + 0.001 * (warm.phase.surface_temperature_k - 300.0));
    assert_eq!(warm.geometry[1].liquid_density_kg_m3(), density);
    near(
        warm.geometry[1].radius_m(),
        (3.0 * initial.mass_kg / (4.0 * core::f64::consts::PI * density)).cbrt(),
        2e-17,
    );
    near(
        initial.mass_kg * warm.droplet.specific_enthalpy_j_kg,
        heat,
        2e-11,
    );
}

#[test]
fn interface_sealing_inventory_refusal_and_cancellation_retain_the_original_seed() {
    let law = Phase {
        expansion_per_k: 0.0,
    };
    let initial = droplet();
    let sealed = KineticEvaporation::try_new(
        &law,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let idle = advance_droplet_with_heat(
        &sealed,
        film(),
        initial,
        1.0,
        0.0,
        controls(),
        temperature,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(idle.droplet, initial);
    let concentration =
        0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0) * (99_000.0_f64 / 96_000.0).ln();
    let sphere = initial.geometry(&law, 300.0).unwrap();
    let rate = 4.0 * core::f64::consts::PI * sphere.radius_m() * D * concentration;
    assert!(matches!(
        advance_droplet_with_heat(
            &law,
            film(),
            initial,
            2.0 * initial.mass_kg / rate,
            2.0 * initial.mass_kg * LATENT,
            controls(),
            temperature,
            || Ok(()),
        ),
        Err(EvaporationError::InventoryLimit { .. })
    ));
    assert!(matches!(
        advance_droplet_with_heat(
            &law,
            film(),
            initial,
            0.1,
            0.0,
            controls(),
            temperature,
            || Err(EvaporationError::Cancelled),
        ),
        Err(EvaporationError::Cancelled)
    ));
    assert_eq!(initial, droplet());
}
