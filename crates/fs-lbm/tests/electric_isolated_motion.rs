//! Isolated material conduction through the production current/force/heat owners.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ElectricMotionInput, ElectricSetup, ElectricalContact, Electrode,
    ElectrodeDrive, MovingConductorSolution, Voxel, VoxelDomain,
};

fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} != {b}"); }
fn input<'a>(u: &'a [[f64;3]], b: &'a [[f64;3]]) -> ElectricMotionInput<'a> {
    ElectricMotionInput { velocity_m_s:u, magnetic_flux_density_t:b, max_cells:64, power_tolerance_w:1e-9 }
}
fn loop_field() -> [[f64;3];4] { [[0.0,0.0,1.0],[0.0,0.0,3.0],[0.0,0.0,1.0],[0.0,0.0,3.0]] }
fn square(rho: f64, contacts: &[ElectricalContact]) -> MovingConductorSolution {
    let domain=VoxelDomain::new(2,2,1,1.0).unwrap();
    ElectricSetup::default().solve_isolated_cell_motion(&domain,&[300.0;4],contacts,
        input(&[[1.0,0.0,0.0];4],&loop_field()),1e-10,
        |_,_,_| Ok(Some(rho)),&CancelGate::new()).unwrap()
}

#[test]
fn closed_loop_heats_and_brakes_without_a_physical_voltage_reference() {
    let result=square(1.0,&[]);
    // Four unit resistances in series; circulating EMF magnitude 2 V.
    // A numerical reference must NOT erase the first cell's resistance/force.
    assert!(result.electric.electrodes.is_empty());
    close(result.electric.delivered_w,0.0,0.0);
    for (p,expected) in result.electric.potential_v.iter().zip([0.0,0.5,-1.5,-2.0]) { close(*p,expected,1e-10); }
    for q in &result.electric.joule_w { close(*q,0.25,1e-10); }
    close(result.electric.total_joule_w,1.0,1e-10);
    close(result.mechanical_power_w,-1.0,1e-10);
    assert!(result.lorentz_force_density_n_m3[0].iter().any(|f| f.abs()>0.1));
    // Independent full continuity check at all four cells, including gauge 0.
    let mut residual=[0.0;4];
    for (a,b,emf) in [(0,1,0.0),(0,2,-1.0),(1,3,-3.0),(2,3,0.0)] {
        let current=result.electric.potential_v[a]-result.electric.potential_v[b]+emf;
        close(current.abs(),0.5,1e-10);
        residual[a]+=current;residual[b]-=current;
    }
    for r in residual { close(r,0.0,1e-10); }
    close(result.power_defect_w,0.0,1e-10);close(result.transfer_defect_w,0.0,1e-10);
}

#[test]
fn material_and_contact_data_change_circulation_without_changing_physics_code() {
    let base=square(1.0,&[]);let resistive=square(2.0,&[]);
    close(resistive.electric.total_joule_w,0.5*base.electric.total_joule_w,1e-10);
    for (a,b) in resistive.lorentz_force_density_n_m3.iter().flatten().zip(base.lorentz_force_density_n_m3.iter().flatten()) {
        close(*a,0.5*b,1e-10);
    }
    let contact=ElectricalContact {cells:[0,1],resistance_ohm_m2:3.0,heat_fraction_to_first:0.25};
    let joined=square(1.0,&[contact]);
    let reversed=square(1.0,&[ElectricalContact {cells:[1,0],heat_fraction_to_first:0.75,..contact}]);
    assert_eq!(joined,reversed);
    let i2=(2.0_f64/7.0).powi(2);
    for (q,r) in joined.electric.joule_w.iter().zip([1.75,3.25,1.0,1.0]) { close(*q,i2*r,1e-10); }
    close(joined.electric.total_joule_w,4.0/7.0,1e-10);
    close(joined.mechanical_power_w,-4.0/7.0,1e-10);
}

#[test]
fn uniform_translation_in_uniform_field_does_not_invent_closed_body_drag() {
    let domain=VoxelDomain::graded([vec![0.5,1.5],vec![0.7,1.3],vec![0.8,1.2]]).unwrap();
    let u=[[1.0,0.0,0.0];8];let b=[[0.0,0.0,2.0];8];
    let origin=domain.center(0,0,0);
    let result=ElectricSetup::default().solve_isolated_cell_motion(&domain,&[300.0;8],&[],input(&u,&b),1e-9,
        |c,_,_|Ok(Some(0.5+c as f64)),&CancelGate::new()).unwrap();
    for c in 0..8 {
        let [x,y,z]=domain.coords(c);let point=domain.center(x,y,z);
        close(result.electric.potential_v[c],-2.0*(point[1]-origin[1]),1e-9);
    }
    assert!(result.electric.total_joule_w < 1e-18);
    assert!(result.lorentz_force_density_n_m3.iter().flatten().all(|f| f.abs()<1e-9));
}

#[test]
fn disconnected_bodies_each_receive_a_gauge_and_keep_their_own_currents() {
    let domain=VoxelDomain::new(5,2,1,1.0).unwrap();let u=[[1.0,0.0,0.0];10];
    let b:Vec<_>=(0..10).map(|c|[0.0,0.0,if c%5==1 || c%5==4 {3.0}else{1.0}]).collect();
    let result=ElectricSetup::default().solve_isolated_cell_motion(&domain,&[300.0;10],&[],input(&u,&b),1e-9,
        |c,_,_|Ok(if c%5==2 {None}else{Some(if c%5<2 {1.0}else{2.0})}),&CancelGate::new()).unwrap();
    close(result.electric.potential_v[0],0.0,0.0);close(result.electric.potential_v[3],0.0,0.0);
    for c in [2,7] { assert!(result.electric.potential_v[c].is_nan());close(result.electric.joule_w[c],0.0,0.0); }
    for c in [0,1,5,6] { close(result.electric.joule_w[c],0.25,1e-10); }
    for c in [3,4,8,9] { close(result.electric.joule_w[c],0.125,1e-10); }
    close(result.electric.total_joule_w,1.5,1e-9);close(result.mechanical_power_w,-1.5,1e-9);
    let single=VoxelDomain::new(3,1,1,1.0).unwrap();
    let result=ElectricSetup::default().solve_isolated_cell_motion(&single,&[300.0;3],&[],
        input(&[[1.0,0.0,0.0];3],&[[0.0,0.0,2.0];3]),1e-10,
        |c,_,_|Ok(if c==1 {None}else{Some(1.0)}),&CancelGate::new()).unwrap();
    assert_eq!(result.electric.joule_w,vec![0.0;3]);close(result.electric.potential_v[0],0.0,0.0);close(result.electric.potential_v[2],0.0,0.0);
}

#[test]
fn mixed_mechanical_occupancy_does_not_select_electrical_materials() {
    let domain=VoxelDomain::from_fn(2,2,1,1.0,|p|if p[0]<1.0 {Voxel::Solid(27)}else{Voxel::Fluid}).unwrap();
    let before=domain.clone();let mut visited=Vec::new();
    let result=ElectricSetup::default().solve_isolated_cell_motion(&domain,&[300.0;4],&[],
        input(&[[1.0,0.0,0.0];4],&loop_field()),1e-10,|c,occupancy,t| {
            visited.push((c,occupancy,t));Ok(Some(1.0))
        },&CancelGate::new()).unwrap();
    assert_eq!(domain,before);assert_eq!(visited.len(),4);
    assert_eq!(visited[0],(0,Voxel::Solid(27),300.0));assert_eq!(visited[1],(1,Voxel::Fluid,300.0));
    assert_eq!(result,square(1.0,&[]));
}

#[test]
fn isolated_opt_in_preserves_old_refusals_and_rejects_invalid_inputs_without_queries() {
    let domain=VoxelDomain::new(2,2,1,1.0).unwrap();let setup=ElectricSetup::default();
    let u=[[1.0,0.0,0.0];4];let b=loop_field();let gate=CancelGate::new();
    assert!(setup.solve_with_cell_motion(&domain,&[300.0;4],&[],input(&u,&b),
        |_,_,_|panic!("legacy path still requires voltage reference"),&gate).is_err());
    let mut wired=setup.clone();wired.electrodes.push(Electrode {lo:[0,0,0],hi:[1,1,1],drive:ElectrodeDrive::Voltage(0.0)});
    assert!(wired.solve_isolated_cell_motion(&domain,&[300.0;4],&[],input(&u,&b),1e-10,
        |_,_,_|panic!("isolated path refuses electrode before queries"),&gate).is_err());
    for tolerance in [-1.0,f64::NAN,f64::INFINITY] {
        assert!(setup.solve_isolated_cell_motion(&domain,&[300.0;4],&[],input(&u,&b),tolerance,
            |_,_,_|panic!("invalid tolerance precedes queries"),&gate).is_err());
    }
    assert!(setup.solve_isolated_cell_motion(&domain,&[300.0,300.0,300.0,f64::NAN],&[],input(&u,&b),1e-10,
        |_,_,_|panic!("entire temperature field precedes queries"),&gate).is_err());
    assert!(setup.solve_isolated_cell_motion(&domain,&[300.0;4],&[],input(&u,&b[..3]),1e-10,
        |_,_,_|panic!("field shape precedes queries"),&gate).is_err());
    let refusal=ChtError::InvalidInput {field:"source",reason:"missing resistivity".into()};
    assert_eq!(setup.solve_isolated_cell_motion(&domain,&[300.0;4],&[],input(&u,&b),1e-10,
        |_,_,_|Err(refusal.clone()),&gate),Err(refusal));
    assert_eq!(setup.solve_isolated_cell_motion(&domain,&[300.0;4],&[],input(&u,&b),1e-10,
        |_,_,_| {gate.request();Ok(Some(1.0))},&gate),Err(ChtError::Cancelled));
}
