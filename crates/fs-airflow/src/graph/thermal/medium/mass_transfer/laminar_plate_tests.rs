use super::*;
use fs_thermal::sublimation::{
    MOLAR_GAS_CONSTANT, SolverConfig, SublimationError, SublimationLaw, SublimationProperties,
    VaporReservoir,
    coating::{SupportedSurfaceState, advance_supported_surface_with_heat},
};

fn geometry() -> LaminarPlateGeometry {
    LaminarPlateGeometry {
        start_m: 0.0,
        end_m: 0.1,
        speed_m_s: 1.0,
    }
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= 1e-20 + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn analytic_plate_value_and_independent_diffusivity_speed_scaling() {
    let base = evaluate(geometry(), 4e-5, 2e-5).unwrap();
    near(base.reynolds(), 2500.0, 1e-15);
    near(base.schmidt(), 2.0, 1e-15);
    let expected_sh = 0.664 * 50.0 * 2.0_f64.powf(1.0 / 3.0);
    near(base.downstream_mean_sherwood(), expected_sh, 2e-13);
    near(
        base.mass_transfer_coefficient_m_s(),
        expected_sh * 2e-5 / 0.1,
        2e-13,
    );
    let more_diffusive = evaluate(geometry(), 4e-5, 4e-5).unwrap();
    near(
        more_diffusive.mass_transfer_coefficient_m_s() / base.mass_transfer_coefficient_m_s(),
        2.0_f64.powf(2.0 / 3.0),
        2e-13,
    );
    let faster = evaluate(
        LaminarPlateGeometry {
            speed_m_s: 4.0,
            ..geometry()
        },
        4e-5,
        2e-5,
    )
    .unwrap();
    near(
        faster.mass_transfer_coefficient_m_s(),
        2.0 * base.mass_transfer_coefficient_m_s(),
        2e-15,
    );
}

#[test]
fn refined_and_near_zero_width_strips_keep_the_physical_leading_edge() {
    let base = evaluate(geometry(), 4e-5, 2e-5).unwrap();
    let mut integrated = 0.0;
    for [start_m, end_m] in [[0.0, 0.001], [0.001, 0.037], [0.037, 0.1]] {
        let strip = evaluate(
            LaminarPlateGeometry {
                start_m,
                end_m,
                ..geometry()
            },
            4e-5,
            2e-5,
        )
        .unwrap();
        integrated += strip.mass_transfer_coefficient_m_s() * (end_m - start_m);
    }
    near(
        integrated,
        base.mass_transfer_coefficient_m_s() * 0.1,
        3e-13,
    );
    let thin = evaluate(
        LaminarPlateGeometry {
            start_m: 0.1_f64.next_down(),
            ..geometry()
        },
        4e-5,
        2e-5,
    )
    .unwrap();
    near(
        thin.mass_transfer_coefficient_m_s(),
        0.5 * base.mass_transfer_coefficient_m_s(),
        2e-15,
    );
}

#[test]
fn zero_flow_transition_low_schmidt_and_unrepresentable_coefficients_refuse() {
    assert!(
        evaluate(
            LaminarPlateGeometry {
                speed_m_s: 0.0,
                ..geometry()
            },
            4e-5,
            2e-5
        )
        .is_err()
    );
    assert!(
        evaluate(
            LaminarPlateGeometry {
                start_m: 0.1,
                ..geometry()
            },
            4e-5,
            2e-5
        )
        .is_err()
    );
    assert!(
        evaluate(
            LaminarPlateGeometry {
                end_m: f64::NAN,
                ..geometry()
            },
            4e-5,
            2e-5
        )
        .is_err()
    );
    assert!(matches!(
        evaluate(
            LaminarPlateGeometry {
                speed_m_s: 201.0,
                ..geometry()
            },
            4e-5,
            2e-5
        ),
        Err(LaminarPlateFilmError::OutsideCorrelation { .. })
    ));
    assert!(matches!(
        evaluate(geometry(), 4e-5, 1e-4),
        Err(LaminarPlateFilmError::OutsideCorrelation { .. })
    ));
    assert!(evaluate(geometry(), 0.0, 2e-5).is_err());
    assert!(evaluate(geometry(), 4e-5, 0.0).is_err());
    // Positive real transport may not silently underflow to a sealed film.
    assert!(
        evaluate(
            LaminarPlateGeometry {
                end_m: 1e308,
                speed_m_s: 1e-300,
                start_m: 0.0
            },
            1e308,
            1e-308
        )
        .is_err()
    );
    let base = evaluate(geometry(), 4e-5, 2e-5).unwrap();
    assert!(base.surface_kinetics(1.01).is_err());
    assert_eq!(
        base.surface_kinetics(0.0)
            .unwrap()
            .accommodation_coefficient,
        0.0
    );
}

struct Deposit;
impl SublimationLaw for Deposit {
    fn temperature_range_k(&self) -> [f64; 2] {
        [250.0, 350.0]
    }
    fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: 10.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}

#[test]
fn original_coating_owner_receives_the_film_and_conserves_bare_deposition_energy() {
    let film = evaluate(geometry(), 4e-5, 2e-5).unwrap();
    let state = SupportedSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 150000.0,
        deposit_mass_kg: 0.0,
        area_m2: 0.01,
    };
    let kinetic = (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * 300.0)).sqrt();
    let pressure_conductance =
        film.mass_transfer_coefficient_m_s() / (MOLAR_GAS_CONSTANT / 0.018 * 300.0);
    let expected_mass = -10.0 / (1.0 / kinetic + 1.0 / pressure_conductance) * state.area_m2 * 10.0;
    let cooling = expected_mass * 2.8e6;
    let accepted = advance_supported_surface_with_heat(
        &Deposit,
        film.surface_kinetics(1.0).unwrap(),
        state,
        VaporReservoir::Prescribed {
            partial_pressure_pa: 20.0,
            temperature_k: 300.0,
        },
        10.0,
        cooling,
        SolverConfig::default(),
        |h| Ok(h / 500.0),
        |t| Ok(2000.0 * t),
        || Ok(()),
    )
    .unwrap();
    near(accepted.phase.transferred_mass_kg, expected_mass, 2e-12);
    near(accepted.state.deposit_mass_kg, -expected_mass, 2e-12);
    near(accepted.phase.surface_temperature_k, 300.0, 2e-13);
    near(
        accepted.substrate_energy_change_j
            + accepted.deposit_energy_change_j
            + accepted.phase.vapor_energy_j,
        cooling,
        2e-10,
    );
    assert_eq!(accepted.state.substrate_mass_kg, state.substrate_mass_kg);
    assert!(accepted.phase.recession_m < 0.0);
}
