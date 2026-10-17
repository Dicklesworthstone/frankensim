use super::*;
use fs_material::visco_tensor::{TensorProny, TensorProjectionBudget};
use fs_math::det;

#[test]
fn projected_kelvin_damping_drives_the_exact_damped_oscillator_and_heat() {
    let gate = CancelGate::new_clock_free();
    let c = core::array::from_fn(|i|core::array::from_fn(|j|if i==j {18.0} else {0.0}));
    let v = core::array::from_fn(|i|core::array::from_fn(|j|if i==j {4.0} else {0.0}));
    let law = TensorProny::new(c,Some(v),vec![],(0.0,f64::MAX)).unwrap();
    let point = law.project_strain(&[[1.0,0.0,0.0,0.0,0.0,0.0]],1.0,
        TensorProjectionBudget {max_coordinates:1,max_scalar_terms:1000},||false).unwrap();
    let body = TensorMaterialBody::new(vec![2.0],vec![point],
        TensorMaterialBudget {max_state_coordinates:16,max_material_points:1},&gate).unwrap();
    assert_eq!(body.state_dimension(),2,"Kelvin loss must not manufacture relaxation states");
    let mut history = body.history(&[0.01],&[0.0],InitialMemory::Relaxed,&gate).unwrap();
    let initial = history.stored_energy_j();
    let wd = det::sqrt(8.0);
    for step in 1..=20 {
        let dt=0.02; let t=step as f64*dt;
        let rec = history.advance(&[0.0],dt,&gate).unwrap();
        let q = 0.01*det::exp(-t)*(det::cos(wd*t)+det::sin(wd*t)/wd);
        let velocity = -0.01*det::exp(-t)*(wd+1.0/wd)*det::sin(wd*t);
        assert!((history.displacement_m()[0]-q).abs()<2e-12);
        assert!((history.velocity_m_s()[0]-velocity).abs()<2e-12);
        assert!(rec.dissipated>=-1e-14);
        assert!((history.dissipated_heat_j()-(initial-9.0*q*q-velocity*velocity)).abs()<2e-12);
        assert!(history.balance_residual_j().abs()<2e-12);
    }
    assert!(history.dissipated_heat_j()>0.0);
}
