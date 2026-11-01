//! The representation floor of the caloric inversion. A temperature moves in
//! quanta q = capacity * (spacing of floats at T), so a request of r quanta can
//! only be answered by 0 or 1 quantum and the old acceptance `min(tol, 0.5 |request|)`
//! could never be met for 1/2 < r < 2/3. At a stagnated bracket the floor accepts the
//! better end when that end MOVED off the reference and its residual is within
//! min(q / 2, the caller's tolerance). A request the temperature cannot follow
//! (under half a quantum) is still refused, exactly as before, and the caller's
//! tolerance is never loosened. Everything that solved before returns the same bits.
use super::*;
use crate::thermochemical_gas::tests::{model, point};

const TEMPERATURES_K: [f64; 5] = [250.0, 300.0, 450.0, 1000.0, 1499.0];

fn controls(tolerance: f64) -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: tolerance, maximum_iterations: 128 }
}

/// The signed change of the model's own caloric law between `t` and its float
/// neighbour on the given side: one quantum of the temperature coordinate.
fn quantum(branch: &GasCaloricBranch<'_>, t: f64, up: bool, constraint: Constraint) -> f64 {
    let neighbour = if up { t.next_up() } else { t.next_down() };
    caloric_change(branch.model(), Temperature::new(t), neighbour, constraint).unwrap()
}

#[test]
fn a_request_between_half_and_two_thirds_of_a_quantum_returns_the_neighbouring_temperature() {
    for t in TEMPERATURES_K {
        let gas = model(vec![0.75, 0.25], 0.002);
        let context = point(t, 100000.0);
        let density = gas.evaluate(&context).unwrap().state().density;
        let branch = gas.prepare_caloric_branch(controls(1e-9), || false).unwrap();
        for up in [true, false] {
            let neighbour = if up { t.next_up() } else { t.next_down() };
            let q = quantum(&branch, t, up, Constraint::Density(density));
            let solved = branch.reconstruct_density_internal_energy_increment(
                &context, density, 0.63 * q, || false).unwrap();
            assert_eq!(solved.state().state().temperature.to_bits(), neighbour.to_bits(),
                "0.63 q at {t} K, up {up}");
            // The request lies 0.37 q from the neighbour and 0.63 q from the reference.
            assert!((solved.residual_j_kg() - 0.37 * q).abs() <= 1e-6 * q.abs(),
                "residual {:e} against {:e}", solved.residual_j_kg(), 0.37 * q);
            // The enthalpy (constant pressure) view of the same dead zone.
            let qh = quantum(&branch, t, up, Constraint::Pressure);
            let solved = branch.reconstruct_pressure_enthalpy_increment(
                &context, 0.63 * qh, || false).unwrap();
            assert_eq!(solved.state().state().temperature.to_bits(), neighbour.to_bits(),
                "0.63 qh at {t} K, up {up}");
        }
    }
}

#[test]
fn a_request_the_temperature_cannot_follow_is_refused_exactly_as_before() {
    // Under half a quantum the closest representable temperature is the reference
    // itself, and a nonzero request is never answered by "no change" (the contract
    // species/network heat_tests pins for a standalone increment).
    for t in TEMPERATURES_K {
        let gas = model(vec![0.75, 0.25], 0.002);
        let context = point(t, 100000.0);
        let density = gas.evaluate(&context).unwrap().state().density;
        let branch = gas.prepare_caloric_branch(controls(1e-9), || false).unwrap();
        for up in [true, false] {
            let constraint = Constraint::Density(density);
            let q = quantum(&branch, t, up, constraint);
            for fraction in [0.4, 0.1, 1e-3, 1e-6] {
                let request = fraction * q;
                assert!(branch.reconstruct_density_internal_energy_increment(
                    &context, density, request, || false).is_err(), "{fraction} q at {t} K");
                assert!(branch.sample_density_internal_energy_increment(
                    &context, density, request, || false).is_err(), "{fraction} q (numeric) at {t} K");
                let old = solve_temperature_with(&branch, &context, constraint, request,
                    &mut || false, false).err().unwrap();
                let new = solve_temperature_with(&branch, &context, constraint, request,
                    &mut || false, true).err().unwrap();
                assert_eq!(old, new, "{fraction} q at {t} K: the floor changed a refusal");
            }
            let qh = quantum(&branch, t, up, Constraint::Pressure);
            assert!(branch.reconstruct_pressure_enthalpy_increment(&context, 0.4 * qh, || false)
                .is_err(), "0.4 qh at {t} K");
        }
    }
}

#[test]
fn the_callers_tolerance_still_caps_the_floor() {
    // 0.63 q leaves 0.37 q against the neighbour. A caller that declared a tolerance
    // below that is refused (as before); one above it is answered. The old inner
    // tolerance min(tol, 0.5 |request|) = 0.315 q refuses both.
    for t in TEMPERATURES_K {
        let gas = model(vec![0.75, 0.25], 0.002);
        let context = point(t, 100000.0);
        let density = gas.evaluate(&context).unwrap().state().density;
        let probe = gas.prepare_caloric_branch(controls(1e-9), || false).unwrap();
        for up in [true, false] {
            let q = quantum(&probe, t, up, Constraint::Density(density)).abs();
            let request = if up { 0.63 * q } else { -0.63 * q };
            for (fraction, accepted) in [(0.30, false), (0.36, false), (0.38, true), (0.45, true)] {
                let branch = gas.prepare_caloric_branch(controls(fraction * q), || false).unwrap();
                let solved = branch.reconstruct_density_internal_energy_increment(
                    &context, density, request, || false);
                assert_eq!(solved.is_ok(), accepted, "tolerance {fraction} q at {t} K, up {up}");
            }
        }
    }
}

#[test]
fn the_solver_floor_is_consulted_only_at_stagnation_and_never_weakens_a_refusal() {
    // A step in the residual: no temperature has a residual below 1 in magnitude.
    let step = |t: f64| -> Result<root::Sample, ()> {
        Ok(root::Sample { residual: if t < 300.0 { -1.0 } else { 1.0 }, tangent: 1.0, pressure_pa: 1.0 })
    };
    let plain = root::solve(299.0, 301.0, 1e-3, 128, step, &mut || false);
    let below = root::solve_with_floor(299.0, 301.0, 1e-3, 128, step, &mut || false, |_, _| 0.5);
    for refused in [&plain, &below] {
        assert!(matches!(refused, Err(root::Error::NotConverged { .. })));
    }
    let (Err(root::Error::NotConverged { iterations: a, residual: ra }),
        Err(root::Error::NotConverged { iterations: b, residual: rb })) = (plain, below)
    else { unreachable!() };
    assert_eq!((a, ra.to_bits()), (b, rb.to_bits()));
    let admitted = root::solve_with_floor(299.0, 301.0, 1e-3, 128, step, &mut || false, |_, _| 1.0)
        .ok().unwrap();
    assert!(admitted.sample.residual.abs() <= 1.0);
    // Where the plain solver converges the floor is never looked at.
    let line = |t: f64| -> Result<root::Sample, ()> {
        Ok(root::Sample { residual: t - 300.5, tangent: 1.0, pressure_pa: 1.0 })
    };
    let plain = root::solve(299.0, 301.0, 1e-9, 128, line, &mut || false).ok().unwrap();
    let floored = root::solve_with_floor(299.0, 301.0, 1e-9, 128, line, &mut || false, |_, _| 1e9)
        .ok().unwrap();
    assert_eq!((plain.temperature_k.to_bits(), plain.iterations),
        (floored.temperature_k.to_bits(), floored.iterations));
}

#[test]
fn every_request_that_solved_before_the_floor_returns_the_same_bits() {
    let mut newly_accepted = 0;
    let mut identical = 0;
    for t in TEMPERATURES_K {
        let gas = model(vec![0.75, 0.25], 0.002);
        let context = point(t, 100000.0);
        let density = gas.evaluate(&context).unwrap().state().density;
        for tolerance in [1e-12, 1e-9, 1e-6] {
            let branch = gas.prepare_caloric_branch(controls(tolerance), || false).unwrap();
            for constraint in [Constraint::Density(density), Constraint::Pressure] {
                let up = quantum(&branch, t, true, constraint);
                let down = quantum(&branch, t, false, constraint);
                let mut requests = vec![];
                for fraction in [0.1, 0.3, 0.45, 0.5, 0.55, 0.6, 0.63, 0.66, 0.67, 0.7, 0.9, 1.0,
                    1.2, 1.5, 2.5, 7.5]
                {
                    requests.push(fraction * up);
                    requests.push(fraction * down);
                }
                for magnitude in [1e-12, 1e-11, 3e-11, 5e-11, 1e-10, 1e-9, 1e-6, 1e-3, 1.0, 1e3] {
                    requests.push(magnitude);
                    requests.push(-magnitude);
                }
                for request in requests {
                    let old = solve_temperature_with(&branch, &context, constraint, request,
                        &mut || false, false);
                    let new = solve_temperature_with(&branch, &context, constraint, request,
                        &mut || false, true);
                    match (old, new) {
                        (Ok(old), Ok(new)) => {
                            identical += 1;
                            assert_eq!(old.temperature_k.to_bits(), new.temperature_k.to_bits());
                            assert_eq!(old.sample.residual.to_bits(), new.sample.residual.to_bits());
                            assert_eq!(old.sample.tangent.to_bits(), new.sample.tangent.to_bits());
                            assert_eq!(old.sample.pressure_pa.to_bits(), new.sample.pressure_pa.to_bits());
                            assert_eq!(old.iterations, new.iterations);
                        }
                        (Err(GasEnergyError::NotConverged { .. }), Ok(new)) => {
                            newly_accepted += 1;
                            // Never the reference, never beyond the floor, never beyond the
                            // caller's tolerance.
                            assert_ne!(new.temperature_k.to_bits(), t.to_bits(),
                                "the reference answered {request:e} at {t} K");
                            let floor = representation_floor(new.temperature_k, new.sample.tangent)
                                .min(tolerance);
                            assert!(new.sample.residual.abs() <= floor,
                                "residual {:e} above floor {floor:e} at {t} K for {request:e}",
                                new.sample.residual);
                        }
                        (Err(old), Err(new)) => assert_eq!(old, new),
                        (old, new) => panic!("{} then {} at {t} K for {request:e}",
                            old.is_ok(), new.is_ok()),
                    }
                }
            }
        }
    }
    assert!(identical > 0 && newly_accepted > 0, "{identical} identical, {newly_accepted} new");
}
