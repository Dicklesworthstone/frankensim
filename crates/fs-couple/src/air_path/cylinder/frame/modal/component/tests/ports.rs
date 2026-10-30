//! G1/G3: the same physical ports and one reciprocal two-body trajectory.
use super::*;
use crate::air_path::cylinder::frame::harmonic_ports::FrameHarmonicPortBudget;
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn port_budget() -> FrameHarmonicPortBudget {
    FrameHarmonicPortBudget { max_ports:4,max_projection_terms:4096,max_mobility_entries:16,
        power_absolute_tolerance_w:1e-10,power_relative_tolerance:1e-9 }
}

#[test]
fn original_point_moments_signed_drives_and_work_use_one_time_frequency_map() {
    let frame=frame(true,1000.0,2e6,500.0);
    let model=model(&frame,true,400.0,0.02);
    let gate=CancelGate::new_clock_free();
    let physical=model.harmonic_transverse_ports(&[(0,0.3),(1,0.7),(1,0.7)],port_budget(),&gate).unwrap();
    let (placement,network)=complete(&model,&initial(&model),false);
    let ports=placement.bind_ports(&physical,4096,&gate).unwrap();
    assert!(core::ptr::eq(ports.original(),&physical));
    let (mut q,mut v)=([0.0;3],[0.0;3]);
    ports.observe_into(network.components(),&mut q,&mut v,&gate).unwrap();
    let (mut u,mut nv)=([0.0;9],[0.0;9]);
    placement.reconstruct_into(network.components(),&mut u,&mut nv,&gate).unwrap();
    for i in 0..3 {
        near(q[i],physical.nodal_columns()[i].iter().zip(u).map(|(b,x)| b*x).sum());
        near(v[i],physical.nodal_columns()[i].iter().zip(nv).map(|(b,x)| b*x).sum());
    }
    assert_ne!(physical.nodal_columns()[0][2],0.0);
    assert_ne!(physical.nodal_columns()[0][5],0.0);
    let efforts=[0.01,-0.003,0.007];
    let mut force=vec![0.0;network.mode_count()];
    ports.accumulate_forces(network.components(),&efforts,&mut force,&gate).unwrap();
    let mut nodal=[0.0;9];
    for (column,effort) in physical.nodal_columns().iter().zip(efforts) {
        for (f,b) in nodal.iter_mut().zip(column) { *f+=b*effort; }
    }
    let native_power:f64=force.iter().zip(network.components().iter().flat_map(|c|c.states()))
        .map(|(f,s)| f*s.velocity_m_sqrt_kg_per_s).sum();
    near(native_power,efforts.iter().zip(v).map(|(f,v)| f*v).sum());
    near(native_power,nodal.iter().zip(nv).map(|(f,v)| f*v).sum());
    let expected=placement.project_nodal_force(&nodal,&gate).unwrap();
    for (a,b) in force.iter().zip(expected) { near(*a,b); }
    assert!(force[4..].iter().all(|f| *f==0.0));
    force.fill(0.0);
    ports.accumulate_forces(network.components(),&[0.0,1.0,-1.0],&mut force,&gate).unwrap();
    assert!(force.iter().all(|f| *f==0.0));
    assert_eq!(ports.attachment(1).unwrap(),ports.attachment(2).unwrap());
    assert!(ports.attachment(3).is_err());
}

fn monolithic(a:&MaterialFluidFrame<'_,FluidTransportState>,b:&MaterialFluidFrame<'_,FluidTransportState>,
    a_column:&[f64],b_column:&[f64],stiffness:f64,damping:f64)->PortHamiltonian {
    let n_a=a.system().state_dim();let n_b=b.system().state_dim();let n=n_a+n_b;
    let mut q=vec![0.0;n*n];let mut j=q.clone();let mut r=q.clone();let mut g=vec![0.0;n*2];
    for (offset,source) in [(0,a),(n_a,b)] {
        let ns=source.system().state_dim();let (js,rs,_)=source.system().structure();
        let qs=source.system().storage().exact_quadratic_matrix().unwrap();
        for row in 0..ns {for col in 0..ns {
            let out=(offset+row)*n+offset+col;let at=row*ns+col;
            q[out]=qs[at];j[out]=js[at];r[out]=rs[at];
        }}
    }
    let mut spring=vec![0.0;n];let mut resistance=vec![0.0;n];
    let ma=a.basis().modes().len();let mb=b.basis().modes().len();
    for (i,&value) in a_column.iter().enumerate() {
        spring[i]=value;resistance[ma+i]=value;g[(ma+i)*2]=value;
    }
    for (i,&value) in b_column.iter().enumerate() {
        spring[n_a+i]=-value;resistance[n_a+mb+i]=-value;g[(n_a+mb+i)*2+1]=value;
    }
    for row in 0..n {for col in 0..n {
        q[row*n+col]+=stiffness*spring[row]*spring[col];
        r[row*n+col]+=damping*resistance[row]*resistance[col];
    }}
    PortHamiltonian::new(n,2,j,r,g,Box::new(QuadraticStorage::new(q,n).unwrap())).unwrap()
}

#[test]
fn two_original_wet_bodies_exchange_reciprocal_force_in_one_complete_native_step() {
    let a_frame=frame(false,1000.0,2e6,500.0);
    let b_frame=frame(false,800.0,3e6,200.0);
    let a=model(&a_frame,false,400.0,0.02);
    let b=model(&b_frame,false,3.0,0.002);
    let gate=CancelGate::new_clock_free();
    let a_ports=a.harmonic_transverse_ports(&[(0,0.7)],port_budget(),&gate).unwrap();
    let b_ports=b.harmonic_transverse_ports(&[(0,0.3)],port_budget(),&gate).unwrap();
    let a_initial=initial(&a);let b_initial=vec![0.0;b.system().state_dim()];
    let (a_map,mut components,mut connections)=a.coupled_component(&a_initial,2000,
        ModalAcousticTimeBudget::audible_reference(),config(),&gate).unwrap().into_parts(0,0,&gate).unwrap();
    let (b_map,b_components,b_connections)=b.coupled_component(&b_initial,2000,
        ModalAcousticTimeBudget::audible_reference(),config(),&gate).unwrap()
        .into_parts(components.len(),connections.len(),&gate).unwrap();
    components.extend(b_components);connections.extend(b_connections);
    let a_port=a_map.bind_ports(&a_ports,4096,&gate).unwrap();
    let b_port=b_map.bind_ports(&b_ports,4096,&gate).unwrap();
    let joint=connections.len();
    connections.push(ModalConnection {left:a_port.attachment(0).unwrap(),right:b_port.attachment(0).unwrap(),
        stiffness_n_m:7.0,damping_n_s_m:0.03,rest_extension_m:0.0});
    let mut network=CoupledModalSystem::new(components,connections,config(),&gate).unwrap();
    let original=monolithic(&a,&b,&a_ports.modal_columns()[0],&b_ports.modal_columns()[0],7.0,0.03);
    let mut exact=PreparedExactZoh::new(&original,0.0005,
        LinearSplitBudget {max_modes:original.state_dim(),max_omega_dt:core::f64::consts::PI}).unwrap();
    let mut state=a_initial;state.extend(b_initial);
    let mut force=vec![0.0;network.mode_count()];
    let (mut aq,mut av,mut bq,mut bv)=([0.0],[0.0],[0.0],[0.0]);
    let mut received=false;
    for step in 0..24 {
        let effort=[if step<8 {0.01} else {0.0},0.0];
        force.fill(0.0);
        a_port.accumulate_forces(network.components(),&effort[..1],&mut force,&gate).unwrap();
        b_port.accumulate_forces(network.components(),&effort[1..],&mut force,&gate).unwrap();
        let expected=exact.step(&state,&effort).unwrap();
        let actual=network.step_under_gate(&force,&gate).unwrap().clone();
        let mut recovered=a_map.physical_state(network.components(),&gate).unwrap();
        recovered.extend(b_map.physical_state(network.components(),&gate).unwrap());
        for (x,y) in recovered.iter().zip(&expected.x) {near(*x,*y);}
        near(actual.external_work_j,expected.supplied);
        near(actual.total_dissipation_j(),expected.dissipated);
        near(network.total_energy_j().unwrap(),original.hamiltonian(&expected.x));
        a_port.observe_into(network.components(),&mut aq,&mut av,&gate).unwrap();
        b_port.observe_into(network.components(),&mut bq,&mut bv,&gate).unwrap();
        near(actual.connection_forces_n[joint],-7.0*(aq[0]-bq[0])-0.03*(av[0]-bv[0]));
        received|=bq[0].abs()>1e-12;
        state=expected.x;
    }
    assert!(received,"unforced second body must receive actual connector work");
    assert_eq!(network.samples_rendered(),24);
    assert_eq!(b_map.structural_component(),2);
    assert!(a_map.connections().fluid_regions.last().unwrap().end<=b_map.connections().elastic.start);
}

#[test]
fn physical_port_binding_refuses_foreign_sources_incomplete_vectors_and_cancellation() {
    let frame=frame(false,1000.0,2e6,500.0);
    let a=model(&frame,false,400.0,0.02);let b=model(&frame,false,400.0,0.02);
    let gate=CancelGate::new_clock_free();
    let pa=a.harmonic_transverse_ports(&[(0,0.7)],port_budget(),&gate).unwrap();
    let pb=b.harmonic_transverse_ports(&[(0,0.7)],port_budget(),&gate).unwrap();
    let (placement,network)=complete(&a,&initial(&a),false);
    assert!(placement.bind_ports(&pb,4096,&gate).is_err());
    assert!(placement.bind_ports(&pa,0,&gate).is_err());
    let port=placement.bind_ports(&pa,4096,&gate).unwrap();
    let saved=placement.physical_state(network.components(),&gate).unwrap();
    assert!(port.accumulate_forces(network.components(),&[1.0],&mut [],&gate).is_err());
    let mut force=vec![0.0;network.mode_count()];
    assert!(port.accumulate_forces(network.components(),&[],&mut force,&gate).is_err());
    assert!(port.accumulate_forces(network.components(),&[f64::NAN],&mut force,&gate).is_err());
    assert!(port.observe_into(network.components(),&mut [],&mut [],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(port.accumulate_forces(network.components(),&[1.0],&mut force,&cancelled).is_err());
    assert_eq!(placement.physical_state(network.components(),&gate).unwrap(),saved);
    assert_eq!(network.samples_rendered(),0);
}
