use super::*;
use crate::sublimation::{SublimationProperties, molecular::molecular_mass_flux_with_kinetics};

const M: f64 = 0.020;
const TS: f64 = 400.0;
const AREA: f64 = 0.02;
const KINETICS: SurfaceKinetics = SurfaceKinetics {
    accommodation_coefficient: 0.37,
    mass_transfer_coefficient_m_s: None,
};

// Manufactured frozen-rate law, not a sourced material correlation.
struct Law(f64);
impl SublimationLaw for Law {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 800.0]
    }
    fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: self.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: 200_000.0,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: M,
        })
    }
}
fn deposit(mass_kg: f64) -> SolidSurfaceState {
    SolidSurfaceState {
        mass_kg,
        specific_enthalpy_j_kg: -7.0e6,
        area_m2: AREA,
    }
}
fn finite_gas(mass: f64, tg: f64, volume: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: mass,
        volume_m3: volume,
        temperature_k: tg,
    }
}
fn vapor_mass(reservoir: VaporReservoir) -> f64 {
    match reservoir {
        VaporReservoir::Finite { vapor_mass_kg, .. } => vapor_mass_kg,
        _ => panic!("finite reservoir required"),
    }
}
fn rates(psat: f64, tg: f64, volume: f64) -> (f64, f64, f64) {
    // Independent dimensional Hertz-Knudsen and ideal-gas reduction.
    let outgoing = AREA
        * KINETICS.accommodation_coefficient
        * psat
        * (M / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * TS)).sqrt();
    let lambda = AREA
        * KINETICS.accommodation_coefficient
        * (M / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * tg)).sqrt()
        * (MOLAR_GAS_CONSTANT / M)
        * tg
        / volume;
    (outgoing, lambda, outgoing / lambda)
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn residuals(step: &MolecularTurnover, scale: f64) {
    for residual in [
        step.stream_mass_residual_kg,
        step.deposit_mass_residual_kg,
        step.vapor_mass_residual_kg,
    ] {
        assert!(residual.abs() <= 32.0 * f64::EPSILON * scale, "{step:?}");
    }
}

#[test]
fn g1_finite_vapor_matches_affine_solution_with_inventory_recycling() {
    let law = Law(1000.0);
    let tg = 500.0;
    let volume = 0.02;
    let (outgoing, lambda, equilibrium) = rates(law.0, tg, volume);
    for ratio in [0.2, 1.6] {
        let initial = ratio * equilibrium;
        let duration = 250.0 / lambda;
        let step = integrate_molecular_turnover(
            &law,
            TS,
            KINETICS,
            deposit(2.0 * equilibrium),
            finite_gas(initial, tg, volume),
            duration,
        )
        .unwrap();
        let net = (equilibrium - initial) * -(-lambda * duration).exp_m1();
        let expected_vapor = equilibrium + (initial - equilibrium) * (-lambda * duration).exp();
        close(step.transferred_mass_kg, net, 0.0, 2e-14);
        close(vapor_mass(step.reservoir), expected_vapor, 0.0, 2e-14);
        close(step.outgoing_mass_kg, outgoing * duration, 0.0, 2e-14);
        close(step.incoming_mass_kg, outgoing * duration - net, 0.0, 2e-14);
        assert!(step.outgoing_mass_kg > 100.0 * 2.0 * equilibrium);
        assert!(step.incoming_mass_kg > 100.0 * initial);
        assert_eq!(step.extinction_time_s, None);
        residuals(&step, outgoing * duration);
    }
}

#[test]
fn g1_extinction_time_and_gross_streams_compose_across_the_bare_boundary() {
    let law = Law(1000.0);
    let (outgoing, lambda, equilibrium) = rates(law.0, TS, 0.02);
    let initial_vapor = 0.2 * equilibrium;
    let initial_deposit = 0.3 * equilibrium;
    let initial = finite_gas(initial_vapor, TS, 0.02);
    let duration = 2.0 / lambda;
    let expected_event = -(1.0_f64 - 0.3 / 0.8).ln() / lambda;
    let whole = integrate_molecular_turnover(
        &law,
        TS,
        KINETICS,
        deposit(initial_deposit),
        initial,
        duration,
    )
    .unwrap();
    close(whole.extinction_time_s.unwrap(), expected_event, 0.0, 3e-14);
    assert_eq!(whole.deposit_mass_kg, 0.0);
    assert_eq!(whole.transferred_mass_kg, initial_deposit);
    close(
        vapor_mass(whole.reservoir),
        initial_vapor + initial_deposit,
        0.0,
        2e-14,
    );
    let expected_out = outgoing * expected_event + 0.5 * outgoing * (duration - expected_event);
    close(whole.outgoing_mass_kg, expected_out, 0.0, 3e-14);
    close(
        whole.incoming_mass_kg,
        expected_out - initial_deposit,
        0.0,
        3e-14,
    );
    for first_duration in [0.5 * expected_event, 1.5 * expected_event] {
        let first = integrate_molecular_turnover(
            &law,
            TS,
            KINETICS,
            deposit(initial_deposit),
            initial,
            first_duration,
        )
        .unwrap();
        let second = integrate_molecular_turnover(
            &law,
            TS,
            KINETICS,
            deposit(first.deposit_mass_kg),
            first.reservoir,
            duration - first_duration,
        )
        .unwrap();
        assert_eq!(second.deposit_mass_kg, 0.0);
        close(
            vapor_mass(second.reservoir),
            vapor_mass(whole.reservoir),
            0.0,
            3e-14,
        );
        close(
            first.outgoing_mass_kg + second.outgoing_mass_kg,
            whole.outgoing_mass_kg,
            0.0,
            3e-14,
        );
        close(
            first.incoming_mass_kg + second.incoming_mass_kg,
            whole.incoming_mass_kg,
            0.0,
            3e-14,
        );
    }
    residuals(&whole, whole.outgoing_mass_kg);
}

#[test]
fn g3_bare_birth_includes_reemission_and_undersaturation_cannot_evaporate_support() {
    let law = Law(1000.0);
    let (_, lambda, equilibrium) = rates(law.0, 500.0, 0.02);
    let duration = 0.3 / lambda;
    let initial = finite_gas(2.0 * equilibrium, 500.0, 0.02);
    let whole =
        integrate_molecular_turnover(&law, TS, KINETICS, deposit(0.0), initial, duration).unwrap();
    assert!(whole.deposit_mass_kg > 0.0);
    assert!(whole.outgoing_mass_kg > 0.0 && whole.incoming_mass_kg > whole.outgoing_mass_kg);
    assert_eq!(whole.extinction_time_s, None);
    let first =
        integrate_molecular_turnover(&law, TS, KINETICS, deposit(0.0), initial, duration / 2.0)
            .unwrap();
    let second = integrate_molecular_turnover(
        &law,
        TS,
        KINETICS,
        deposit(first.deposit_mass_kg),
        first.reservoir,
        duration / 2.0,
    )
    .unwrap();
    close(second.deposit_mass_kg, whole.deposit_mass_kg, 0.0, 2e-14);
    close(
        first.incoming_mass_kg + second.incoming_mass_kg,
        whole.incoming_mass_kg,
        0.0,
        2e-14,
    );
    let undersaturated = finite_gas(0.5 * equilibrium, 500.0, 0.02);
    let bare =
        integrate_molecular_turnover(&law, TS, KINETICS, deposit(0.0), undersaturated, duration)
            .unwrap();
    assert_eq!(bare.reservoir, undersaturated);
    assert_eq!(bare.deposit_mass_kg, 0.0);
    assert_eq!(bare.transferred_mass_kg, 0.0);
    assert_eq!(bare.incoming_mass_kg, bare.outgoing_mass_kg);
    assert!(bare.incoming_mass_kg > 0.0);
    assert_eq!(bare.extinction_time_s, Some(0.0));
    // Complete thermalization of this equal counterflow can exchange heat.
    assert!(bare.outgoing_mass_kg * (TS - 500.0) < 0.0);
}

#[test]
fn g1_prescribed_bath_retains_original_film_rates_and_exact_loss_limit() {
    let law = Law(1000.0);
    let kinetics = SurfaceKinetics {
        mass_transfer_coefficient_m_s: Some(0.03),
        ..KINETICS
    };
    for pressure in [300.0, 1800.0] {
        let reservoir = VaporReservoir::Prescribed {
            partial_pressure_pa: pressure,
            temperature_k: 500.0,
        };
        let flux = molecular_mass_flux_with_kinetics(&law, TS, kinetics, reservoir).unwrap();
        let initial = if flux.net_kg_m2_s > 0.0 { 1e-6 } else { 0.0 };
        let duration = if initial > 0.0 {
            2.0 * initial / (AREA * flux.net_kg_m2_s)
        } else {
            0.7
        };
        let step =
            integrate_molecular_turnover(&law, TS, kinetics, deposit(initial), reservoir, duration)
                .unwrap();
        assert_eq!(step.reservoir, reservoir);
        close(
            step.incoming_mass_kg,
            AREA * flux.incoming_kg_m2_s * duration,
            0.0,
            2e-14,
        );
        if initial > 0.0 {
            assert_eq!(step.transferred_mass_kg, initial);
            assert_eq!(step.deposit_mass_kg, 0.0);
            close(step.extinction_time_s.unwrap(), duration / 2.0, 0.0, 2e-14);
        } else {
            close(
                step.deposit_mass_kg,
                -AREA * flux.net_kg_m2_s * duration,
                0.0,
                2e-14,
            );
            assert!(step.outgoing_mass_kg > 0.0);
        }
        residuals(&step, step.incoming_mass_kg.max(step.outgoing_mass_kg));
    }
}

#[test]
fn g1_empty_vapor_retains_quadratic_return_and_tiny_extinction_time() {
    let law = Law(1000.0);
    let (outgoing, lambda, equilibrium) = rates(law.0, TS, 0.02);
    let duration = 1e-16 / lambda;
    let step = integrate_molecular_turnover(
        &law,
        TS,
        KINETICS,
        deposit(equilibrium),
        finite_gas(0.0, TS, 0.02),
        duration,
    )
    .unwrap();
    assert!(step.incoming_mass_kg > 0.0);
    close(
        step.incoming_mass_kg,
        outgoing * duration * (lambda * duration) / 2.0,
        0.0,
        2e-14,
    );
    close(step.transferred_mass_kg, outgoing * duration, 0.0, 2e-14);
    let tiny_deposit = equilibrium * 1e-17;
    let tiny_time = tiny_deposit / outgoing;
    let exhausted = integrate_molecular_turnover(
        &law,
        TS,
        KINETICS,
        deposit(tiny_deposit),
        finite_gas(0.0, TS, 0.02),
        2.0 * tiny_time,
    )
    .unwrap();
    assert_eq!(exhausted.deposit_mass_kg, 0.0);
    assert_eq!(exhausted.transferred_mass_kg, tiny_deposit);
    close(exhausted.extinction_time_s.unwrap(), tiny_time, 0.0, 2e-14);
    assert!(exhausted.incoming_mass_kg > 0.0);
}

#[test]
fn g0_exact_saturation_zero_time_and_sealed_interfaces_preserve_inventories() {
    let reservoir = finite_gas(0.001, TS, 0.02);
    let law = Law(reservoir.partial_pressure_pa(M).unwrap());
    for mass in [0.0, 0.003] {
        let step =
            integrate_molecular_turnover(&law, TS, KINETICS, deposit(mass), reservoir, 1000.0)
                .unwrap();
        assert_eq!(step.reservoir, reservoir);
        assert_eq!(step.deposit_mass_kg, mass);
        assert_eq!(step.transferred_mass_kg, 0.0);
        assert_eq!(step.incoming_mass_kg, step.outgoing_mass_kg);
        assert_eq!(step.extinction_time_s, None);
        residuals(&step, 0.0);
        for (kinetics, duration) in [
            (KINETICS, 0.0),
            (
                SurfaceKinetics {
                    accommodation_coefficient: 0.0,
                    ..KINETICS
                },
                1.0,
            ),
            (
                SurfaceKinetics {
                    mass_transfer_coefficient_m_s: Some(0.0),
                    ..KINETICS
                },
                1.0,
            ),
        ] {
            let step = integrate_molecular_turnover(
                &Law(0.0),
                TS,
                kinetics,
                deposit(mass),
                reservoir,
                duration,
            )
            .unwrap();
            assert_eq!(step.reservoir, reservoir);
            assert_eq!(step.deposit_mass_kg, mass);
            assert_eq!(step.incoming_mass_kg, 0.0);
            assert_eq!(step.outgoing_mass_kg, 0.0);
            assert_eq!(step.extinction_time_s, None);
        }
    }
}

#[test]
fn g0_representable_late_vapor_tail_and_encoding_residual_are_retained() {
    let volume = 1e300;
    let mass = 1e300;
    let (_, lambda, _) = rates(0.0, TS, volume);
    let step = integrate_molecular_turnover(
        &Law(0.0),
        TS,
        KINETICS,
        deposit(0.0),
        finite_gas(mass, TS, volume),
        800.0 / lambda,
    )
    .unwrap();
    let expected = (mass.ln() - 800.0).exp();
    assert!(vapor_mass(step.reservoir) > 0.0);
    close(vapor_mass(step.reservoir), expected, 0.0, 5e-13);
    // The transferred mass rounds to the original inventory, but the actual
    // vapor has not vanished. TwoSum preserves this smaller encoding defect.
    if step.transferred_mass_kg == -mass {
        assert_eq!(step.vapor_mass_residual_kg, vapor_mass(step.reservoir));
    }
    assert_eq!(step.extinction_time_s, None);
}

#[test]
fn g0_invalid_inputs_refuse_without_creating_a_candidate() {
    let law = Law(1000.0);
    let gas = finite_gas(0.001, TS, 0.02);
    for mass in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(integrate_molecular_turnover(&law, TS, KINETICS, deposit(mass), gas, 0.1).is_err());
    }
    for duration in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            integrate_molecular_turnover(&law, TS, KINETICS, deposit(0.01), gas, duration).is_err()
        );
    }
    for kinetics in [
        SurfaceKinetics {
            accommodation_coefficient: 1.1,
            ..KINETICS
        },
        SurfaceKinetics {
            mass_transfer_coefficient_m_s: Some(-1.0),
            ..KINETICS
        },
    ] {
        assert!(integrate_molecular_turnover(&law, TS, kinetics, deposit(0.01), gas, 0.1).is_err());
    }
    assert!(integrate_molecular_turnover(&law, 199.0, KINETICS, deposit(0.0), gas, 0.0).is_err());
    assert!(
        integrate_molecular_turnover(
            &law,
            TS,
            KINETICS,
            deposit(0.0),
            finite_gas(-1.0, TS, 0.02),
            0.1
        )
        .is_err()
    );
}
