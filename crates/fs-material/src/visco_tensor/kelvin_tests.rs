//! Kelvin viscosity must survive the shared material's history/projection adapters.
use super::*;

fn law(kv: bool) -> TensorProny {
    let equilibrium = core::array::from_fn(|i| core::array::from_fn(|j|
        if i == j { 12.0 + i as f64 } else { 0.0 }));
    let mut b = [[0.0; 6]; 6];
    b[0] = [2.0, -1.0, 0.0, 3.0, 0.0, 0.0];
    b[1] = [0.0, 1.0, 0.5, 0.0, 2.0, -1.0];
    let branch = TensorPronyBranch::from_gram_rows(core::array::from_fn(|i|
        core::array::from_fn(|j| if i==j {0.5+i as f64/10.0} else {0.0})),0.08).unwrap();
    TensorProny::new(equilibrium, kv.then(|| gram(&b).unwrap()), vec![branch],
        (0.0, f64::MAX)).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-11*(1.0+a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn kelvin_ramps_retain_all_stress_components_and_independent_heat() {
    let material = law(true);
    let elastic = law(false);
    let viscosity = material.kv_pa_s().unwrap();
    let mut history = material.history();
    let mut plain = elastic.history();
    let mut previous = [0.0; 6];
    for next in [[0.01,-0.02,0.005,0.03,-0.01,0.02],
        [-0.03,0.01,0.02,-0.01,0.03,-0.01], [0.0;6]] {
        let dt = 0.03;
        let rate = core::array::from_fn(|i| (next[i]-previous[i])/dt);
        let mandel_rate = mandel(&rate);
        let viscous = tensor(&apply(viscosity, &mandel_rate));
        let heat = dt * dot(&mandel_rate, &apply(viscosity, &mandel_rate));
        let actual = history.advance(next, dt, || false).unwrap();
        let expected = plain.advance(next, dt, || false).unwrap();
        for i in 0..6 { close(actual.stress[i], expected.stress[i]+viscous[i]); }
        close(actual.stored_energy, expected.stored_energy);
        close(actual.work, expected.work+heat);
        close(actual.dissipated, expected.dissipated+heat);
        close(actual.balance_residual, 0.0);
        assert!(actual.dissipated > 0.0);
        previous = next;
    }
    let held = history.advance(previous, 0.03, || false).unwrap();
    assert!(held.dissipated>=0.0);
    close(history.work(), history.stored_energy()+history.dissipated());
}

#[test]
fn kelvin_consistent_tangent_and_instantaneous_jump_admission() {
    let material = law(true);
    let base = material.history();
    let target = [0.002,-0.001,0.0,0.003,0.001,-0.002];
    let tangent = material.tangent(0.04).unwrap();
    let h = 1e-7;
    for col in 0..6 {
        let mut plus = target; plus[col] += h;
        let mut minus = target; minus[col] -= h;
        let a = base.clone().advance(plus,0.04,||false).unwrap();
        let b = base.clone().advance(minus,0.04,||false).unwrap();
        for row in 0..6 { close((a.stress[row]-b.stress[row])/(2.0*h),tangent[row][col]); }
    }
    let mut history = base.clone();
    assert!(material.tangent(0.0).is_err());
    assert!(history.advance(target,0.0,||false).is_err());
    assert_eq!(history,base);
    assert!(history.advance(target,0.04,||true).is_err());
    assert_eq!(history,base);
    // A dt of f64::MIN_POSITIVE stays finite for this law and is admitted, as the
    // documented positive-duration ramp allows. The smallest subnormal step overflows
    // the Kelvin rate, and that refusal must publish no partial state.
    assert!(history.advance(target,f64::from_bits(1),||false).is_err());
    assert_eq!(history,base);
    let empty_jump = history.advance([0.0;6],0.0,||false).unwrap();
    assert_eq!(empty_jump.dissipated,0.0);
    let zero = TensorProny::new(*material.c_eq(),Some([[0.0;6];6]),vec![],
        (0.0,f64::MAX)).unwrap();
    assert!(zero.tangent(0.0).is_ok());
    assert_eq!(zero.history().advance(target,0.0,||false).unwrap().dissipated,0.0);
}

#[test]
fn projected_kelvin_force_power_and_volume_use_the_same_strain_map() {
    let material = law(true);
    let columns = [[1.0,0.3,0.0,-0.5,0.2,0.0],[-0.4,0.7,0.0,0.1,0.0,0.2]];
    let budget = TensorProjectionBudget {max_coordinates:2,max_scalar_terms:10000};
    let point = material.project_strain(&columns,0.2,budget,||false).unwrap();
    let twice = material.project_strain(&columns,0.4,budget,||false).unwrap();
    let v = [0.03,-0.02];
    let strain_rate = core::array::from_fn(|i| columns[0][i]*v[0]+columns[1][i]*v[1]);
    let stress = tensor(&apply(material.kv_pa_s().unwrap(),&mandel(&strain_rate)));
    let expected = point.conjugate_force(&stress).unwrap();
    let c = point.viscosity_matrix();
    assert_eq!(c[1].to_bits(),c[2].to_bits());
    for i in 0..2 {
        close(c[2*i]*v[0]+c[2*i+1]*v[1],expected[i]);
        for j in 0..2 {close(twice.viscosity_matrix()[2*i+j],2.0*c[2*i+j]);}
    }
    close(v[0]*expected[0]+v[1]*expected[1],
        0.2*crate::tensor::contract(&stress,&strain_rate));
    assert!(point.viscosity_matrix().iter().any(|x|*x!=0.0));
    let pure = law(false).project_strain(&columns,0.2,budget,||false).unwrap();
    assert!(pure.viscosity_matrix().iter().all(|x|*x==0.0));
}
