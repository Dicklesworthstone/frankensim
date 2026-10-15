//! Real parcel/constitutive paths with synthetic gas and geometric test data.
use super::*;
use crate::thermochemical_gas::energy::entropy::GasEntropySolve;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
}
fn parcel(x: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope), point(400.0, 100000.0), 0.001).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance*b.abs().max(1.0), "{a} != {b}");
}
fn volume(q: [f64; 2]) -> (f64, [f64; 2]) {
    (0.001 + 0.01*q[0] - 0.002*q[1] + 0.015*q[0]*q[0]
        + 0.004*q[0]*q[1] + 0.003*q[1]*q[1],
        [0.01 + 0.03*q[0] + 0.004*q[1], -0.002 + 0.004*q[0] + 0.006*q[1]])
}
fn hessian_action(v: [f64; 2]) -> [f64; 2] {
    [0.03*v[0] + 0.004*v[1], 0.004*v[0] + 0.006*v[1]]
}

#[test]
fn curved_geometry_force_and_full_tangent_match_fresh_material_states() {
    let base = parcel(0.75, 0.002);
    let q = [0.008, 0.015];
    let d = [0.4, -0.7];
    let (v, g) = volume(q);
    let mut gas = base.clone();
    gas.set_volume_isentropic(v, controls(), || false).unwrap();
    let load = gas.mechanical_load(v, &g, 75000.0).unwrap();
    let action = load.tangent_action(&d, &hessian_action(d)).unwrap();
    let epsilon = 1e-6;
    let eval = |sign: f64| {
        let (vol, grad) = volume([q[0]+sign*epsilon*d[0], q[1]+sign*epsilon*d[1]]);
        let mut trial = base.clone();
        trial.set_volume_isentropic(vol, controls(), || false).unwrap();
        trial.mechanical_load(vol, &grad, 75000.0).unwrap().generalized_forces().unwrap()
    };
    let (plus, minus) = (eval(1.0), eval(-1.0));
    for i in 0..2 { near(action[i], (plus[i]-minus[i])/(2.0*epsilon), 2e-7); }
    let frozen_geometry = load.tangent_action(&d, &[0.0; 2]).unwrap();
    assert!((action[0]-frozen_geometry[0]).abs() > 100.0);
    let e = [-0.2, 0.9];
    let other = load.tangent_action(&e, &hessian_action(e)).unwrap();
    near(e[0]*action[0]+e[1]*action[1], d[0]*other[0]+d[1]*other[1], 1e-12);
    assert!(std::ptr::eq(load.state(), gas.state()));
}

#[test]
fn forces_are_the_energy_gradient_and_rate_pairing_is_power_conjugate() {
    let gas = parcel(0.75, 0.001);
    let gradient = [0.01]; // straight piston area; q in metres
    let load = gas.mechanical_load(0.001, &gradient, 30000.0).unwrap();
    let force = load.generalized_forces().unwrap()[0];
    near(force, 700.0, 1e-12);
    near(load.mechanical_power_w(&[-0.2]).unwrap(), -140.0, 1e-12);
    let h = 1e-6;
    let potential = |sign: f64| {
        let v = 0.001+gradient[0]*sign*h;
        let mut candidate = gas.clone();
        candidate.set_volume_isentropic(v, controls(), || false).unwrap();
        gas.mass_kg()*candidate.state().specific_internal_energy_j_kg() + 30000.0*v
    };
    near(force, -(potential(1.0)-potential(-1.0))/(2.0*h), 2e-7);
    let response = load.pressure_response();
    let medium = gas.state().state();
    near(response.bulk_modulus_pa, medium.density*medium.sound_speed.powi(2), 1e-12);
    near(load.tangent_action(&[1.0], &[0.0]).unwrap()[0],
        -response.bulk_modulus_pa/0.001*gradient[0].powi(2), 1e-12);
}

#[test]
fn coordinate_rescaling_and_permutation_preserve_force_work_and_tangent() {
    let gas = parcel(0.75, 0.001);
    let g = [0.01, -0.002]; let d = [0.4, -0.7];
    let h = hessian_action(d);
    let load = gas.mechanical_load(0.001, &g, 50000.0).unwrap();
    // q=[2*r1,-3*r0]; new derivative arrays follow the SAME chart transform.
    let gp = [-3.0*g[1], 2.0*g[0]];
    let dp = [-d[1]/3.0, d[0]/2.0];
    let hp = [-3.0*h[1], 2.0*h[0]];
    let mapped = gas.mechanical_load(0.001, &gp, 50000.0).unwrap();
    near(load.mechanical_power_w(&d).unwrap(), mapped.mechanical_power_w(&dp).unwrap(), 1e-12);
    let a = load.tangent_action(&d, &h).unwrap();
    let b = mapped.tangent_action(&dp, &hp).unwrap();
    near(b[0], -3.0*a[1], 1e-12); near(b[1], 2.0*a[0], 1e-12);
}

#[test]
fn composition_and_heat_change_mechanical_and_acoustic_response_without_geometry_changes() {
    let mut a = parcel(0.75, 0.001); let mut b = parcel(0.25, 0.001);
    for gas in [&mut a, &mut b] {
        gas.set_volume_isentropic(0.0005, controls(), || false).unwrap();
    }
    let g = [0.01];
    let load = a.mechanical_load(0.0005, &g, 100000.0).unwrap();
    let force = load.generalized_forces().unwrap()[0];
    let stiffness = load.tangent_action(&[1.0], &[0.0]).unwrap()[0];
    let sound = load.state().state().sound_speed;
    assert_ne!(force, b.mechanical_load(0.0005, &g, 100000.0).unwrap().generalized_forces().unwrap()[0]);
    a.add_heat(30.0, GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-6, maximum_iterations: 128 }, || false).unwrap();
    let hot = a.mechanical_load(0.0005, &g, 100000.0).unwrap();
    assert!(hot.generalized_forces().unwrap()[0] > force);
    assert!(hot.tangent_action(&[1.0], &[0.0]).unwrap()[0] < stiffness);
    assert!(hot.state().state().sound_speed > sound);
}

#[test]
fn stale_geometry_invalid_actions_and_unrepresentable_loads_refuse() {
    let gas = parcel(0.75, 0.001);
    for p in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(gas.mechanical_load(0.001, &[0.01], p).is_err());
    }
    for v in [0.002, f64::NAN, 0.0] { assert!(gas.mechanical_load(v, &[0.01], 0.0).is_err()); }
    for g in [vec![], vec![f64::NAN], vec![1.0; MAX_VOLUME_COORDINATES+1]] {
        assert!(gas.mechanical_load(0.001, &g, 0.0).is_err());
    }
    let load = gas.mechanical_load(0.001, &[0.01], 100000.0).unwrap();
    assert_eq!(load.generalized_forces().unwrap(), vec![0.0]);
    assert!(load.tangent_action(&[1.0], &[0.0]).unwrap()[0] < 0.0);
    assert!(load.tangent_action(&[], &[0.0]).is_err());
    assert!(load.tangent_action(&[1.0], &[f64::NAN]).is_err());
    assert!(load.mechanical_power_w(&[f64::INFINITY]).is_err());
    let huge = [f64::MAX];
    let load = gas.mechanical_load(0.001, &huge, 0.0).unwrap();
    assert!(load.generalized_forces().is_err());
    assert!(load.tangent_action(&[f64::MAX], &[0.0]).is_err());
    let g = vec![0.0; MAX_VOLUME_COORDINATES];
    assert_eq!(gas.mechanical_load(0.001, &g, 0.0).unwrap().generalized_forces().unwrap(), g);
}
