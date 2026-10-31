//! Independent fixed-charge thermodynamics of the explicitly declared
//! D(E,T)=s(T)D_ref(E) profile. All source data are manufactured.
use super::*;

#[derive(Clone, Copy)]
struct ReferenceState {
    free: f64,
    entropy: f64,
    internal: f64,
    capacity: f64,
}

// Hand integrals of the three source segments, with A=gap=1. The source has
// (E,D_ref)=(0,0),(1,2),(3,3),(4,7), and s=1+alpha*(T-300).
// Fixed physical D means D_ref=|D|/s; it is not a fixed-field derivative.
fn reference_state(charge: f64, temperature: f64, alpha: f64) -> ReferenceState {
    let scale = 1.0 + alpha * (temperature - 300.0);
    let d = charge.abs() / scale;
    let (free_reference, coenergy_reference, differential_permittivity) = if d < 2.0 {
        (d * d / 4.0, d * d / 4.0, 2.0)
    } else if d < 3.0 {
        (d * d - 3.0 * d + 3.0, d * d - 3.0, 0.5)
    } else {
        (
            d * d / 8.0 + 9.0 * d / 4.0 - 39.0 / 8.0,
            d * d / 8.0 + 39.0 / 8.0,
            4.0,
        )
    };
    let free = scale * free_reference;
    let entropy = alpha * coenergy_reference;
    ReferenceState {
        free,
        entropy,
        internal: free + temperature * entropy,
        capacity: -temperature * alpha * alpha / scale * d * d / differential_permittivity,
    }
}

#[test]
fn fixed_charge_heating_reaches_independent_internal_energy_root_across_moving_charge_knots() {
    let body_capacity = 0.1;
    let initial_temperature = 300.0;
    let final_temperature = 320.0;
    let dt = 0.2;
    for alpha in [-0.001, 0.001] {
        // The second charge crosses the first physical q knot in opposite
        // directions for the two slopes; source E knots remain fixed.
        for magnitude in [2.5, 2.0 + 10.0 * alpha] {
            for polarity in [-1.0, 1.0] {
                let charge = polarity * magnitude;
                let before = reference_state(charge, initial_temperature, alpha);
                let after = reference_state(charge, final_temperature, alpha);
                let delta_u = after.internal - before.internal;
                let external_heat =
                    body_capacity * (final_temperature - initial_temperature) + delta_u;
                assert!(external_heat > 0.0);
                assert!(before.capacity < 0.0 && after.capacity < 0.0);
                assert!(body_capacity + before.capacity > 0.0);
                assert!(body_capacity + after.capacity > 0.0);
                let card = fixture_scaled(body_capacity, false, 0.0, [250.0, 500.0], alpha, 0.0);
                let caloric = caloric(&card);
                let thermal = thermal(&[&caloric], &[external_heat / dt]);
                let context = point(initial_temperature);
                let cells = [cell(0, &card, &caloric, &context)];
                let network = MaterialNetwork {
                    cells: &cells,
                    resistors: &[],
                    links: &[],
                };
                let (circuit, pole_law, instantaneous_law) =
                    initial_scaled(&card, charge, None, false);
                let original_state = circuit.state().to_vec();
                let original_thermal_energy = thermal.energy_j();
                let mut bound = binding(&pole_law, &[]);
                bound.maximum_instantaneous_charge_c = 3.0;
                let candidate = dielectric::step_candidate_with_reversible_storage(
                    &network,
                    &thermal,
                    dt,
                    &circuit,
                    &[],
                    bound,
                    &instantaneous_law,
                    controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
                let report = candidate.report();
                let source = &report.coupled.source;
                near(
                    candidate.thermal().temperature_k(DIELECTRIC).unwrap(),
                    final_temperature,
                    3e-8,
                );
                for (actual, expected) in source.dielectric_state.into_iter().zip([before, after]) {
                    near(actual.free_energy_j, expected.free, 2e-10);
                    near(actual.entropy_j_k, expected.entropy, 2e-12);
                    near(actual.internal_energy_j, expected.internal, 2e-10);
                    near(actual.heat_capacity_j_k, expected.capacity, 2e-12);
                }
                near(source.parameter_work_j, after.free - before.free, 2e-10);
                near(
                    source.thermal_offset_change_j,
                    final_temperature * after.entropy - initial_temperature * before.entropy,
                    2e-10,
                );
                near(source.dielectric_internal_energy_change_j, delta_u, 2e-10);
                near(source.reversible_heat_j, -delta_u, 2e-10);
                near(report.external_heat_j, external_heat, 2e-10);
                near(
                    candidate.thermal().energy_j() - original_thermal_energy,
                    body_capacity * (final_temperature - initial_temperature),
                    2e-10,
                );
                near(report.total_energy_residual_j, 0.0, 4e-10);
                assert_eq!(report.coupled.electrical.polarization_heat_j, 0.0);
                assert_eq!(
                    candidate.circuit().instantaneous_charge_c().to_bits(),
                    charge.to_bits()
                );
                assert_eq!(circuit.state(), original_state);
                assert_eq!(thermal.energy_j(), original_thermal_energy);
                assert_eq!(thermal.time_s(), 0.0);
            }
        }
    }
}

#[test]
fn positive_caloric_capacity_cannot_hide_negative_combined_fixed_charge_capacity() {
    let body_capacity = 0.002;
    for alpha in [-0.001, 0.001] {
        let charge = 2.5;
        let expected = reference_state(charge, 300.0, alpha);
        near(expected.capacity, -0.00375, 2e-18);
        assert!(body_capacity > 0.0 && body_capacity + expected.capacity < 0.0);
        let card = fixture_scaled(body_capacity, false, 0.0, [250.0, 500.0], alpha, 0.0);
        let caloric = caloric(&card);
        let thermal = thermal(&[&caloric], &[0.01]);
        let context = point(300.0);
        let cells = [cell(0, &card, &caloric, &context)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let (circuit, pole_law, instantaneous_law) = initial_scaled(&card, charge, None, false);
        let saved_state = circuit.state().to_vec();
        let saved_energy = thermal.energy_j();
        let mut bound = binding(&pole_law, &[]);
        bound.maximum_instantaneous_charge_c = 3.0;
        match dielectric::step_candidate_with_reversible_storage(
            &network,
            &thermal,
            0.1,
            &circuit,
            &[],
            bound,
            &instantaneous_law,
            controls(),
            &[0.0],
            || false,
        ) {
            Err(CircuitEnthalpyError::Source(dielectric::DebyeThermalError::Invalid(message))) => {
                assert!(
                    message.contains("capacit"),
                    "unexpected admission refusal: {message}"
                );
            }
            Err(error) => panic!("expected combined-capacity admission refusal, got {error:?}"),
            Ok(_) => panic!("negative total fixed-charge capacity cannot publish a thermal step"),
        }
        assert_eq!(circuit.state(), saved_state);
        assert_eq!(thermal.energy_j(), saved_energy);
        assert_eq!(thermal.time_s(), 0.0);
    }
}
