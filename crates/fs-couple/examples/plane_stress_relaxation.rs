//! Declared synthetic coupon kinematics, not a measured material or instrument.
//! A prescribed parabolic bending shape supplies one fixed curvature column.
//! Two thickness stations integrate its quadratic strain energy exactly;
//! no branch or normal-strain response is fitted. This is not a mesh solution.
use fs_couple::tensor_material::{InitialMemory, PlaneStressPoint, TensorMaterialBody, TensorMaterialBudget};
use fs_exec::CancelGate;
use fs_material::{IsotropicElastic, visco_tensor::{TensorProny, TensorPronyBranch, TensorProjectionBudget}};
use fs_math::det;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let gate=CancelGate::new_clock_free();
    let (length,width,thickness,density)=(0.2_f64,0.05_f64,0.002_f64,1200.0_f64);
    let elastic=IsotropicElastic::new(2.0e8,0.3,0.01)?.stiffness();
    let mut factor=[[0.0;6];6];
    factor[0]=[3000.0,-1000.0,2000.0,0.0,0.0,0.0];
    factor[1]=[0.0,0.0,0.0,1500.0,0.0,0.0];
    let material=TensorProny::from_tensor_shear(elastic,
        vec![TensorPronyBranch::from_gram_rows(factor,0.02)?],0.01)?;
    let volume=length*width*thickness/2.0;
    let station=thickness/det::sqrt(12.0);
    let limits=TensorProjectionBudget{max_coordinates:4,max_scalar_terms:1_000_000};
    let points=[-station,station].into_iter().map(|z|
        PlaneStressPoint::new(&material,&[[8.0*z/(length*length),0.0,0.0]],volume,limits,&gate))
        .collect::<Result<Vec<_>,_>>()?;
    // Integral of (4*x*(L-x)/L^2)^2 over [0,L] is 8*L/15.
    let inertia=density*length*width*thickness*(8.0/15.0);
    let body=TensorMaterialBody::from_plane_stress(vec![inertia],points,
        TensorMaterialBudget{max_state_coordinates:64,max_material_points:2},&gate)?;
    let mut history=body.history(&[0.0],&[0.0],InitialMemory::Relaxed,&gate)?;
    println!("time_s,centre_displacement_m,normal_strain,normal_stress_pa,stored_j,work_j,heat_j,balance_j");
    for tick in 0..200 {
        history.advance(&[if tick<10 {0.002} else {0.0}],0.0005,&gate)?;
        println!("{},{},{},{},{},{},{},{}",history.time_s(),history.displacement_m()[0],
            history.plane_stress_strain(0)?[2],history.plane_stress_normal_stress_pa(0)?[0],
            history.stored_energy_j(),history.supplied_work_j(),history.dissipated_heat_j(),
            history.balance_residual_j());
    }
    Ok(())
}
