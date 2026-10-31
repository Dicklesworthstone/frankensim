//! Independent continuous closed-system oracle. It eliminates temperature with
//! conserved physical internal energy and integrates only the pole charge.
//! No native circuit, material response, or production thermal kernel is used
//! by this reference calculation.
use super::*;

const CALORIC_CAPACITY: f64 = 0.3;
const POLE_SLOPE: f64 = 0.0006;
const TOTAL_CHARGE: f64 = 4.0;
const END_TIME: f64 = 0.4;

fn reference_internal(q: f64, p: f64, t: f64, alpha: f64) -> f64 {
    let scale = 1.0 + alpha * (t - 300.0);
    let z = q / scale;
    let chi = 1.0 + POLE_SLOPE * (t - 300.0);
    let free = scale * density(z) + 0.5 * p * p / chi;
    let entropy = alpha * (z * voltage(z) - density(z)) + 0.5 * POLE_SLOPE * (p / chi).powi(2);
    free + t * entropy
}

fn reference_temperature(p: f64, alpha: f64, initial_energy: f64) -> f64 {
    let q = TOTAL_CHARGE - p;
    let (mut lower, mut upper) = (250.0, 500.0);
    let residual =
        |t| CALORIC_CAPACITY * (t - 300.0) + reference_internal(q, p, t, alpha) - initial_energy;
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..64 {
        let middle = f64::midpoint(lower, upper);
        if residual(middle) < 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    f64::midpoint(lower, upper)
}

fn reference_polarization(alpha: f64, steps: u32) -> f64 {
    let initial_energy = reference_internal(1.8, 2.2, 300.0, alpha);
    let rhs = |p| {
        let t = reference_temperature(p, alpha, initial_energy);
        let scale = 1.0 + alpha * (t - 300.0);
        let chi = 1.0 + POLE_SLOPE * (t - 300.0);
        // The authored rate curve is flat below 300 K and above 301 K.
        let tau = 0.2 + 2.0 * (t - 300.0).clamp(0.0, 1.0);
        (chi * voltage((TOTAL_CHARGE - p) / scale) - p) / tau
    };
    let dt = END_TIME / f64::from(steps);
    let mut p = 2.2;
    for _ in 0..steps {
        let a = rhs(p);
        let b = rhs(p + 0.5 * dt * a);
        let c = rhs(p + 0.5 * dt * b);
        let d = rhs(p + dt * c);
        p += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    p
}

#[test]
fn temperature_dependent_nonlinear_memory_refines_to_independent_closed_dynamics() {
    let mut final_temperatures = Vec::new();
    for alpha in [-0.001, 0.001] {
        let card = fixture_scaled(
            CALORIC_CAPACITY,
            true,
            0.0,
            [250.0, 500.0],
            alpha,
            POLE_SLOPE,
        );
        let curve = caloric(&card);
        let context = point(300.0);
        let cells = [cell(0, &card, &curve, &context)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let initial_energy = reference_internal(1.8, 2.2, 300.0, alpha);
        let reference = reference_polarization(alpha, 12000);
        near(reference, reference_polarization(alpha, 6000), 1e-8);
        let reference_t = reference_temperature(reference, alpha, initial_energy);
        final_temperatures.push(reference_t);
        let run = |steps: u32| {
            let mut heat = thermal(&[&curve], &[0.0]);
            let (mut circuit, pole_law, field_law) = initial_scaled(&card, 1.8, Some(2.2), false);
            let mut total_loss = 0.0;
            for step in 1..=steps {
                let mut bound = binding(&pole_law, &[3.0]);
                bound.maximum_instantaneous_charge_c = 3.5;
                let next = dielectric::step_candidate_with_reversible_storage(
                    &network,
                    &heat,
                    END_TIME * f64::from(step) / f64::from(steps),
                    &circuit,
                    &[],
                    bound,
                    &field_law,
                    controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
                near(next.report().total_energy_residual_j, 0.0, 4e-10);
                assert!(next.report().coupled.electrical.polarization_heat_j >= 0.0);
                total_loss += next.report().coupled.electrical.polarization_heat_j;
                let (thermal, electrical, _) = next.into_parts();
                heat = thermal;
                circuit = electrical;
            }
            let p = circuit.polarization_charges_c()[0];
            let q = circuit.instantaneous_charge_c();
            let t = heat.temperature_k(DIELECTRIC).unwrap();
            near(q + p, TOTAL_CHARGE, 2e-12);
            near(
                CALORIC_CAPACITY * (t - 300.0) + reference_internal(q, p, t, alpha),
                initial_energy,
                3e-8,
            );
            near(t, reference_temperature(p, alpha, initial_energy), 1e-7);
            let scale = 1.0 + alpha * (t - 300.0);
            assert!(
                q / scale > 2.0,
                "trajectory must cross the physical nonlinear knot"
            );
            near(circuit.terminal_voltage_v(), voltage(q / scale), 2e-12);
            assert!(total_loss > 0.0);
            assert!(circuit.thermal_reversible_storage().is_some());
            ((p - reference).abs(), (t - reference_t).abs())
        };
        let errors = [run(20), run(40), run(80)];
        for pair in errors.windows(2) {
            assert!(
                pair[1].0 < 0.7 * pair[0].0,
                "charge refinement errors {errors:?}"
            );
            assert!(
                pair[1].1 < 0.7 * pair[0].1,
                "thermal refinement errors {errors:?}"
            );
        }
        assert!(errors[2].0 < 0.005, "independent charge error {errors:?}");
    }
    // Changing the sign of the field's thermal response substantially changes
    // the actual relaxation trajectory, despite identical 300 K charge laws.
    assert!((final_temperatures[0] - final_temperatures[1]).abs() > 0.5);
}
