use super::*;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::tests::{model, point, close};

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}

#[test]
fn pressure_entropy_matches_constant_cp_and_nonlinear_integral() {
    for slope in [0.0, 0.002] {
        let model = model(vec![0.75, 0.25], slope);
        let branch = model.prepare_caloric_branch(controls(), || false).unwrap();
        let context = point(400.0, 100000.0);
        let r = model.evaluate(&context).unwrap().state().specific_gas_constant;
        for (t1, p1) in [(600.0_f64, 400000.0_f64), (300.0, 50000.0)] {
            let target = r * (3.25 * (t1 / 400.0).ln()
                + 0.75 * slope * (t1 - 400.0) - (p1 / 100000.0).ln());
            let result = branch.reconstruct_pressure_entropy_change(&context,
                Pressure::new(p1), target, 1e-9, || false).unwrap();
            close(result.state().state().temperature, t1);
            assert_eq!(result.state().state().pressure, p1);
            assert!(result.residual_j_kg_k().abs() <= 1e-9);
        }
        let reference = branch.reconstruct_pressure_entropy_change(&context,
            Pressure::new(100000.0), 0.0, 1e-9, || false).unwrap();
        assert_eq!(reference.state().state().temperature, 400.0);
        assert_eq!(reference.iterations(), 0);
    }
}

#[test]
fn pressure_entropy_keeps_typed_refusals_and_every_cancel_boundary() {
    let model = model(vec![0.75, 0.25], 0.002);
    let branch = model.prepare_caloric_branch(controls(), || false).unwrap();
    let context = point(400.0, 100000.0);
    for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(branch.reconstruct_pressure_entropy_change(&context,
            Pressure::new(pressure), 0.0, 1e-8, || false).is_err());
    }
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(branch.reconstruct_pressure_entropy_change(&context,
            Pressure::new(400000.0), 0.0, tolerance, || false).is_err());
    }
    assert!(branch.reconstruct_pressure_entropy_change(&context,
        Pressure::new(400000.0), f64::NAN, 1e-8, || false).is_err());
    let mut total = 0;
    let expected = branch.reconstruct_pressure_entropy_change(&context, Pressure::new(400000.0),
        0.0, 1e-8, || { total += 1; false }).unwrap();
    for stop in 1..=total {
        let mut polls = 0;
        assert!(matches!(branch.reconstruct_pressure_entropy_change(&context, Pressure::new(400000.0),
            0.0, 1e-8, || { polls += 1; polls == stop }), Err(GasEntropyError::Cancelled)));
    }
    assert_eq!(expected, branch.reconstruct_pressure_entropy_change(&context, Pressure::new(400000.0),
        0.0, 1e-8, || false).unwrap());
}

#[test]
fn tiny_pressure_change_uses_source_bounds_instead_of_a_full_temperature_search() {
    let source = model(vec![0.75, 0.25], 0.0);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let context = point(400.0, 100000.0);
    // Nonzero entropy at the reference is smaller than the authored tolerance.
    // It must not be accepted as an unchanged temperature merely for that reason.
    let pressure = 100000.0 * (1.0 + 1e-9);
    let result = branch.reconstruct_pressure_entropy_change(&context, Pressure::new(pressure),
        0.0, 1e-4, || false).unwrap();
    assert!(result.state().state().temperature > 400.0);
    assert!(result.iterations() < 16);
    let gas = source.evaluate(&context).unwrap();
    let expected = 400.0 * (pressure / 100000.0).powf(
        (gas.state().gamma - 1.0) / gas.state().gamma);
    // Physical tolerance is loose; the derived solver target resolves direction,
    // not an unrequested high-relative-accuracy heat/entropy certificate.
    assert!((result.state().state().temperature - expected).abs() < 2e-7);
}
