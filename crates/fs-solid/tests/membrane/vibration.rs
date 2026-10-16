//! Native tests of the original membrane -> modal ownership chain.
use super::*;
use fs_solid::membrane::vibration::{VibrationError, VibrationSettings};

fn modal_settings() -> VibrationSettings {
    VibrationSettings {max_free_dofs:16,relative_residual_tolerance:1e-8}
}

#[test]
fn loaded_pyramid_frequency_uses_the_actual_geometric_tangent_and_mass() {
    let p=patch(8000.0,8000.0,4000.0,0.2,0.0);
    let h=0.005_f64;let stretch=1.02_f64;let side=0.1_f64;
    let load=2.0*8000.0*(stretch.powi(2)-1.0)*h+4.0*8000.0*h.powi(3)/side.powi(2);
    let mut loads=vec![[0.0;3];5];loads[0][2]=-load;
    let result=p.solve_vibration(&loads,&prescribed(stretch),&[[0.0;3];5],
        settings(),modal_settings(),||false).unwrap();
    assert_eq!(result.free_dofs(),&[2]);assert_eq!(result.modes().len(),1);
    close(result.equilibrium().displacement_m[0][2],-h,1e-10,1e-10);
    let stiffness=2.0*8000.0*(stretch.powi(2)-1.0)+12.0*8000.0*h*h/(side*side);
    let mass=p.reference_mass_kg()/3.0;
    close(result.modes()[0].lambda,stiffness/mass,1e-4,1e-9);
    close(result.modes()[0].phi[0].abs(),1.0/mass.sqrt(),1e-12,1e-12);
    assert!(result.modes()[0].interval.0>0.0);
    // The centre shape integrates to one third of the stretched projected area.
    let projected_area=2.0*side*side*stretch*stretch;
    close(result.volume_coupling()[0].abs(),projected_area/(3.0*mass.sqrt()),1e-12,1e-11);
    assert_eq!(result.patch().sheets(),p.sheets());
}

#[test]
fn source_mass_stiffness_and_supports_change_frequencies_without_retuning() {
    let a=patch(8000.0,8000.0,4000.0,0.2,0.0);
    let stiff=patch(16000.0,16000.0,8000.0,0.2,0.0);
    let heavy=patch(8000.0,8000.0,4000.0,0.6,0.0);
    let zero=[[0.0;3];5];let bc=prescribed(1.02);
    let solve=|p:&MembranePatch|p.solve_vibration(&zero,&bc,&zero,settings(),modal_settings(),||false)
        .unwrap().modes()[0].lambda;
    let base=solve(&a);
    close(solve(&stiff),2.0*base,1e-6,1e-12);
    close(solve(&heavy),base/3.0,1e-6,1e-12);
    let tighter=a.solve_vibration(&zero,&prescribed(1.04),&zero,settings(),modal_settings(),||false).unwrap();
    close(tighter.modes()[0].lambda/base,(1.04_f64.powi(2)-1.0)/(1.02_f64.powi(2)-1.0),1e-10,1e-11);
}

#[test]
fn three_coordinate_modes_preserve_nodal_work_mass_and_oriented_volume_flux() {
    let p=patch(8000.0,2000.0,500.0,0.2,0.37);
    let bc:Vec<_>=prescribed(1.02).into_iter().filter(|(dof,_)|*dof>=3).collect();
    let mut loads=[[0.0;3];5];loads[0]=[0.25,-0.15,-0.12];
    let result=p.solve_vibration(&loads,&bc,&[[0.0;3];5],settings(),modal_settings(),||false).unwrap();
    assert_eq!(result.modes().len(),3);
    let qdot=[0.0004,-0.0007,0.0002];
    let velocity=result.reconstruct(&qdot).unwrap();
    assert!(velocity.iter().skip(1).flatten().all(|v|*v==0.0));
    let forces=[[0.7,-0.2,0.9],[0.2,0.1,0.3],[0.0,0.4,0.0],[-0.5,0.2,0.1],[0.1,0.2,-0.4]];
    let modal_force=result.project_nodal_force(&forces).unwrap();
    let dot=|a:&[f64],b:&[f64]|a.iter().zip(b).map(|(a,b)|a*b).sum::<f64>();
    let work:f64=forces.iter().zip(&velocity).map(|(f,v)|dot(f,v)).sum();
    close(work,dot(&modal_force,&qdot),1e-15,1e-12);
    let kinetic:f64=velocity.iter().zip(p.nodal_mass_kg()).map(|(v,m)|m*dot(v,v)).sum();
    close(kinetic,dot(&qdot,&qdot),1e-18,1e-11);
    // Independently integrate the linearized normal velocity on each actual
    // deformed triangle, including in-plane components of its tilted normal.
    let u=&result.equilibrium().displacement_m;
    let positions:Vec<[f64;3]>=p.mesh().nodes.iter().zip(u).map(|(x,u)|[x[0]+u[0],x[1]+u[1],u[2]]).collect();
    let mut flux=0.0;
    for t in &p.mesh().elems {
        let a:[f64;3]=std::array::from_fn(|i|positions[t[1]][i]-positions[t[0]][i]);
        let b:[f64;3]=std::array::from_fn(|i|positions[t[2]][i]-positions[t[0]][i]);
        let cross=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
        for &node in t {flux+=dot(&cross,&velocity[node])/6.0;}
    }
    close(flux,dot(result.volume_coupling(),&qdot),1e-16,1e-11);
    let pressure=3.7;
    let pressure_forces:Vec<_>=result.pressure_load_per_pa().iter().map(|a|a.map(|v|pressure*v)).collect();
    close(dot(&result.project_nodal_force(&pressure_forces).unwrap(),&qdot),pressure*flux,1e-16,1e-11);
    assert!(result.project_nodal_force(&[[0.0;3];4]).is_err());
    assert!(result.reconstruct(&[f64::NAN;3]).is_err());
}

#[test]
fn nonpositive_spectra_limits_and_cancellation_never_publish_pitched_modes() {
    let p=patch(8000.0,8000.0,4000.0,0.2,0.0);let zero=[[0.0;3];5];
    for stretch in [1.0,0.99] {
        assert!(matches!(p.solve_vibration(&zero,&prescribed(stretch),&zero,settings(),modal_settings(),||false),
            Err(VibrationError::UnresolvedMode {..})));
    }
    let mut cap=modal_settings();cap.max_free_dofs=0;
    assert!(matches!(p.solve_vibration(&zero,&prescribed(1.02),&zero,settings(),cap,||false),
        Err(VibrationError::Membrane(MembraneError::Budget(_)))));
    assert!(matches!(p.solve_vibration(&zero,&prescribed(1.02),&zero,settings(),modal_settings(),||true),
        Err(VibrationError::Membrane(MembraneError::Cancelled))));
    // Count the original equilibrium owner's callbacks, then request cancellation
    // on the first tangent-column boundary of the composed operation.
    let bc=prescribed(1.02);let mut equilibrium_polls=0;
    p.solve(&zero,&bc,&zero,settings(),||{equilibrium_polls+=1;false}).unwrap();
    let cancel_at=1+bc.len()+equilibrium_polls+1;let mut polls=0;
    assert!(matches!(p.solve_vibration(&zero,&bc,&zero,settings(),modal_settings(),||{polls+=1;polls==cancel_at}),
        Err(VibrationError::Membrane(MembraneError::Cancelled))));
    assert!(p.solve_vibration(&zero,&bc,&zero,settings(),modal_settings(),||false).is_ok());
}
