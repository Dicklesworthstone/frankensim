//! Check the actual source-resolved nonlinear midpoint and its coordinate pullback.
use super::*;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
    }
}
fn parcel(x: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0 - x], slope),
        point(400.0, 100000.0), 0.001).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}");
}
fn geometry(q: [f64; 2]) -> (f64, [f64; 2]) {
    (0.001 + 0.01*q[0] - 0.002*q[1] + 0.015*q[0]*q[0]
        + 0.004*q[0]*q[1] + 0.003*q[1]*q[1],
     [0.01 + 0.03*q[0] + 0.004*q[1], -0.002 + 0.004*q[0] + 0.006*q[1]])
}
fn hessian_action(direction: [f64; 2]) -> [f64; 2] {
    [0.03*direction[0] + 0.004*direction[1],
     0.004*direction[0] + 0.006*direction[1]]
}

#[test]
fn nonlinear_curved_geometry_jacobian_matches_fresh_complete_material_states() {
    for (fraction, slope) in [(0.75, 0.0), (0.75, 0.002), (0.25, 0.001)] {
        let gas = parcel(fraction, slope);
        let q = [0.008, 0.015];
        let d = [0.4, -0.7];
        let (volume, gradient) = geometry(q);
        let trial = gas.midpoint_heat_volume(volume, 20.0, controls(), || false).unwrap();
        let load = trial.discrete_mechanical_load(volume, &gradient, 75000.0).unwrap();
        let tangent = load.tangent_action(&d, &hessian_action(d)).unwrap();
        let h = 1e-6;
        let evaluate = |sign: f64| {
            let (v, g) = geometry([q[0]+sign*h*d[0], q[1]+sign*h*d[1]]);
            let state = gas.midpoint_heat_volume(v, 20.0, controls(), || false).unwrap();
            state.discrete_mechanical_load(v, &g, 75000.0).unwrap()
                .generalized_forces().unwrap()
        };
        let (plus, minus) = (evaluate(1.0), evaluate(-1.0));
        for i in 0..2 {
            close(tangent[i], (plus[i]-minus[i])/(2.0*h), 2e-7);
        }
        let omitted_geometry = load.tangent_action(&d, &[0.0; 2]).unwrap();
        assert!((tangent[0]-omitted_geometry[0]).abs() > 100.0);
        let e = [-0.2, 0.9];
        let t_e = load.tangent_action(&e, &hessian_action(e)).unwrap();
        close(e[0]*tangent[0]+e[1]*tangent[1],
              d[0]*t_e[0]+d[1]*t_e[1], 1e-12);
        assert!(std::ptr::eq(load.candidate(), &trial));
    }
}

#[test]
fn nonlinear_heat_to_force_derivative_matches_fresh_source_resolutions() {
    for (fraction, slope) in [(0.75, 0.0), (0.75, 0.002), (0.25, 0.001)] {
        let gas = parcel(fraction, slope);
        let (volume, gradient) = geometry([0.008, 0.015]);
        let trial = gas.midpoint_heat_volume(volume, 20.0, controls(), || false).unwrap();
        let load = trial.discrete_mechanical_load(volume, &gradient, 75000.0).unwrap();
        let dq = 0.001;
        let eval = |heat| {
            let fresh = gas.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
            fresh.discrete_mechanical_load(volume, &gradient, 75000.0).unwrap()
                .generalized_forces().unwrap()
        };
        let plus = eval(20.0+dq);
        let minus = eval(20.0-dq);
        for (i, value) in load.heat_tangent().unwrap().iter().enumerate() {
            close(*value, (plus[i]-minus[i])/(2.0*dq), 2e-7);
        }
    }
}

#[test]
fn pressure_force_power_and_material_swap_share_one_candidate() {
    let (volume, gradient) = geometry([0.008, 0.015]);
    let mut old_force = None;
    for fraction in [0.75, 0.25] {
        let gas = parcel(fraction, 0.002);
        let candidate = gas.midpoint_heat_volume(volume, 30.0, controls(), || false).unwrap();
        let load = candidate.discrete_mechanical_load(volume, &gradient, 50000.0).unwrap();
        let forces = load.generalized_forces().unwrap();
        for i in 0..2 {
            close(forces[i], (candidate.midpoint().state().pressure - 50000.0)*gradient[i], 1e-12);
        }
        let rates = [0.3, -0.6];
        close(load.mechanical_power_w(&rates).unwrap(),
            forces[0]*rates[0]+forces[1]*rates[1], 1e-12);
        if let Some(previous) = old_force { assert_ne!(previous, forces[0]); }
        old_force = Some(forces[0]);
        assert_eq!(candidate.endpoint().mass_kg(), gas.mass_kg());
        assert!(candidate.step().energy_residual_j.abs() <= gas.mass_kg()*controls().energy_tolerance_j_kg);
    }
}

#[test]
fn bad_geometry_extreme_scalars_and_mismatched_vectors_fail_without_mutation() {
    let gas = parcel(0.75, 0.002);
    let original = gas.clone();
    let trial = gas.midpoint_heat_volume(0.001, 20.0, controls(), || false).unwrap();
    for pressure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(trial.discrete_mechanical_load(0.001, &[0.01], pressure).is_err());
    }
    for volume in [0.0, 0.000999, f64::NAN] {
        assert!(trial.discrete_mechanical_load(volume, &[0.01], 0.0).is_err());
    }
    for gradient in [vec![], vec![f64::NAN], vec![1.0; MAX_DISCRETE_GAS_COORDINATES+1]] {
        assert!(trial.discrete_mechanical_load(0.001, &gradient, 0.0).is_err());
    }
    let gradient = [0.01];
    let load = trial.discrete_mechanical_load(0.001, &gradient, 0.0).unwrap();
    assert!(load.tangent_action(&[], &[0.0]).is_err());
    assert!(load.tangent_action(&[0.0], &[]).is_err());
    assert!(load.tangent_action(&[f64::NAN], &[0.0]).is_err());
    assert!(load.mechanical_power_w(&[f64::INFINITY]).is_err());
    assert!(trial.discrete_mechanical_load(0.001, &[f64::MAX], 0.0)
        .unwrap().generalized_forces().is_err());
    assert_eq!(gas, original);
}
