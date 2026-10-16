//! Force-driven anisotropic material-memory demonstration. All coefficients
//! below are declared SYNTHETIC test data, not a calibrated material/instrument.
//! Run: cargo run -p fs-couple --example tensor_relaxation -- 100 0.001
use fs_couple::tensor_material::{InitialMemory,TensorMaterialBody,TensorMaterialBudget};
use fs_exec::CancelGate;
use fs_material::visco_tensor::{TensorProny,TensorPronyBranch,TensorProjectionBudget};

fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>2 {return Err("usage: tensor_relaxation [steps [dt_seconds]]".into());}
    let steps=args.first().map(|s|s.parse::<usize>()).transpose()?.unwrap_or(100);
    let dt=args.get(1).map(|s|s.parse::<f64>()).transpose()?.unwrap_or(0.001);
    if steps==0 || steps>100_000 || !dt.is_finite() || dt<=0.0 {
        return Err("steps must be 1..=100000 and dt finite and positive".into());
    }
    let equilibrium=core::array::from_fn(|i|core::array::from_fn(|j|if i==j {12.0} else {0.0}));
    let mut factor=[[0.0;6];6];factor[0]=[2.0,-1.0,0.0,0.5,0.0,0.0];
    let material=TensorProny::new_mandel(equilibrium,
        vec![TensorPronyBranch::from_gram_rows(factor,0.05)?],0.1)?;
    let columns=[[1.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.2,0.0,0.0]];
    let contribution=material.project_strain(&columns,1.0,
        TensorProjectionBudget {max_coordinates:2,max_scalar_terms:1000},||false)?;
    let gate=CancelGate::new_clock_free();
    let body=TensorMaterialBody::new(vec![2.0,3.0],vec![contribution],
        TensorMaterialBudget {max_state_coordinates:32,max_material_points:1},&gate)?;
    let mut history=body.history(&[0.001,0.0],&[0.0,0.0],InitialMemory::Instantaneous,&gate)?;
    eprintln!("synthetic 3D strain-map material; {} retained states; exact held-force propagation; no radiation or measured-material claim",body.state_dimension());
    println!("time_s,u0_m,u1_m,stored_J,supplied_J,heat_J,balance_J");
    for tick in 0..steps {
        // An explicit 10 ms force pulse, resolved by the caller's held-force grid.
        let force=if tick as f64*dt<0.01 {[0.02,0.0]} else {[0.0,0.0]};
        history.advance(&force,dt,&gate)?;
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            history.time_s(),history.displacement_m()[0],history.displacement_m()[1],
            history.stored_energy_j(),history.supplied_work_j(),history.dissipated_heat_j(),history.balance_residual_j());
    }
    Ok(())
}
