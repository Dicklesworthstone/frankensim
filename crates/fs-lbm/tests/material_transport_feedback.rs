//! Real kernels with synthetic local laws; no calibration authority is implied.
use fs_exec::CancelGate;
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};

fn material(k: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0,1000.0,k,mu).unwrap()
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.1).unwrap()
        .thermal_slab3(MaterialThermalSlab3 { nx:2,ny:1,nz:6,
            bottom_temperature_k:310.0,top_temperature_k:300.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16 }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.1).unwrap()
        .thermal_slab(MaterialThermalSlab { nx:2,ny:6,
            bottom_temperature_k:310.0,top_temperature_k:300.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16 }).unwrap()
}
fn limits(n: usize) -> HeatSourceLimits {
    HeatSourceLimits { max_cells:n,temperature_range_k:[290.0,330.0],balance_tolerance_j:1e-8 }
}
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a} != {b}"); }
fn power() -> Vec<f64> { let mut p=vec![0.0;16];p[4..12].fill(1e6);p }

#[test]
fn constant_material_is_the_existing_physical_step_in_both_dimensions() {
    let mut a=slab3();let mut b=slab3();let gate=CancelGate::new();
    let force=[[1.0,0.0,0.0];16];let p=power();
    let old=a.step_with_power_density(&p,&force,limits(16)).unwrap();
    let new=b.step_with_material_transport(&p,&force,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert_eq!(new.heat,old);assert_eq!(a.lattice.grid.f,b.lattice.grid.f);assert_eq!(a.lattice.gpop,b.lattice.gpop);
    assert_eq!(new.applied.len(),12);assert_eq!(new.endpoint.len(),12);
    for s in &new.applied {
        near(s.momentum_relaxation,0.8,1e-15);near(s.thermal_relaxation,0.9,1e-15);
    }
    let mut a=slab2();let mut b=slab2();let force=[[1.0,0.0];16];
    let old=a.step_with_power_density(&p,&force,0.01,limits(16)).unwrap();
    let new=b.step_with_material_transport(&p,&force,0.01,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert_eq!(new.heat,old);assert_eq!(a.lattice.grid.f,b.lattice.grid.f);assert_eq!(a.lattice.gpop,b.lattice.gpop);
    near(new.applied[0].thermal_relaxation,0.8,1e-15); // not the D3Q7 value
}

#[derive(Debug)]
struct SyntheticTicket { value:FluidTransportState, query_k:f64, condition:u32 }
impl FluidTransport for SyntheticTicket {
    fn transport_state(&self)->Result<FluidTransportState,MaterialError> { Ok(self.value) }
}

#[test]
fn actual_local_temperatures_and_concrete_source_payloads_survive_each_step() {
    let mut a=slab3();let mut control=slab3();let gate=CancelGate::new();
    a.lattice.grid.f[6][0]-=0.002;a.lattice.grid.f[6][1]+=0.001;a.lattice.grid.f[6][2]+=0.001;
    control.lattice.grid.f.clone_from(&a.lattice.grid.f);
    let p=power();let force=[[0.0;3];16];
    for _ in 0..2 {
        let old:Vec<_>=a.lattice.gpop.iter().map(|p|a.temperature_si(p.iter().sum()).unwrap()).collect();
        let r=a.step_with_material_transport(&p,&force,limits(16),&gate,|_,t,retained| {
            assert_eq!(*retained,material(1000.0,1.0));
            Ok(SyntheticTicket { value:material(1000.0*(1.0+(t-300.0)/10.0),1.0+(t-300.0)/10.0),query_k:t,condition:17 })
        }).unwrap();
        for s in &r.applied { assert_eq!(s.temperature_k,old[s.cell]);assert_eq!(s.material.query_k,s.temperature_k);assert_eq!(s.material.condition,17); }
        for s in &r.endpoint { assert_eq!(s.temperature_k,a.temperature_si(a.lattice.gpop[s.cell].iter().sum()).unwrap()); }
        assert!(r.applied[0].transport.thermal_conductivity_w_m_k()>r.applied[11].transport.thermal_conductivity_w_m_k());
        assert_eq!(a.lattice.grid.tau,vec![0.8;16]); // step-local coefficients, retained history
        control.step_with_power_density(&p,&force,limits(16)).unwrap();
    }
    assert_ne!(a.lattice.grid.f,control.lattice.grid.f);assert_ne!(a.lattice.gpop,control.lattice.gpop);
}

fn nonlinear_error(n:usize)->f64 {
    let dx=1.0/n as f64;let dt=80.0*dx*dx;
    let mut a=plan_scaling_from_material(&material(1000.0,1.0),dx,dt,0.0).unwrap()
        .thermal_slab3(MaterialThermalSlab3 {nx:1,ny:1,nz:n,bottom_temperature_k:310.0,
            top_temperature_k:300.0,gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:n+2}).unwrap();
    let gate=CancelGate::new();let p=vec![0.0;n+2];let forces=vec![[0.0;3];n+2];
    for _ in 0..30*n*n {
        a.step_with_material_transport(&p,&forces,limits(n+2),&gate,
            |_,t,_|Ok(material(1000.0*(1.0+(t-300.0)/10.0),1.0))).unwrap();
    }
    let mut error=0.0_f64;
    for i in 1..=n {
        let y=(i as f64-0.5)/n as f64;
        let exact=300.0+10.0*((4.0-3.0*y).sqrt()-1.0);
        error=error.max((a.temperature_si(a.lattice.temperature(0,0,i)).unwrap()-exact).abs());
    }
    error
}
#[test]
fn nonlinear_material_conduction_matches_a_kirchhoff_transform_reference() {
    let coarse=nonlinear_error(8);let fine=nonlinear_error(16);
    assert!(coarse<0.04);assert!(fine<0.3*coarse);
}

#[test]
fn density_capacity_and_endpoint_source_failures_leave_history_and_units_intact() {
    let mut a=slab3();let f=a.lattice.grid.f.clone();let g=a.lattice.gpop.clone();
    let tau=a.lattice.grid.tau.clone();let loads=a.lattice.grid.fext.clone();let gate=CancelGate::new();
    for changed in [
        FluidTransportState::try_new(999.0,1000.0,1000.0,1.0).unwrap(),
        FluidTransportState::try_new(1000.0,1001.0,1000.0,1.0).unwrap(),
    ] {
        let err=a.step_with_material_transport(&power(),&[[0.0;3];16],limits(16),&gate,|_,_,_|Ok(changed)).unwrap_err();
        assert!(matches!(err,MaterialTransportError::Invalid(_)));
    }
    let threshold=a.temperature_si(a.lattice.gpop[6].iter().sum()).unwrap()+0.01;
    let mut p=vec![0.0;16];p[6]=1e8;
    let err=a.step_with_material_transport(&p,&[[0.0;3];16],limits(16),&gate,|cell,t,_| {
        if cell==6 && t>threshold { Err(MaterialError::State { what: "outside the local source temperature support".into() }) }
        else { Ok(material(1000.0,1.0)) }
    }).unwrap_err();
    assert!(matches!(err,MaterialTransportError::State {cell:6,..}));
    assert_eq!(a.lattice.grid.f,f);assert_eq!(a.lattice.gpop,g);assert_eq!(a.lattice.grid.tau,tau);assert_eq!(a.lattice.grid.fext,loads);
    // Retry reaches the original accepted state, with no duplicated heating.
    let mut fresh=slab3();a.step_with_power_density(&power(),&[[0.0;3];16],limits(16)).unwrap();
    fresh.step_with_power_density(&power(),&[[0.0;3];16],limits(16)).unwrap();assert_eq!(a.lattice.gpop,fresh.lattice.gpop);
}

#[test]
fn cancelled_endpoint_and_invalid_inputs_never_publish_a_partial_step() {
    let mut a=slab3();let f=a.lattice.grid.f.clone();let g=a.lattice.gpop.clone();let gate=CancelGate::new();
    let mut calls=0;
    let result=a.step_with_material_transport(&power(),&[[0.0;3];16],limits(16),&gate,|_,_,_| {
        calls+=1;if calls==13 {gate.request();}Ok(material(1000.0,1.0))
    });
    assert!(matches!(result,Err(MaterialTransportError::Cancelled)));
    assert_eq!(a.lattice.grid.f,f);assert_eq!(a.lattice.gpop,g);
    for p in [vec![0.0;15],vec![1.0;16]] {
        assert!(a.step_with_material_transport(&p,&[[0.0;3];16],limits(16),&CancelGate::new(),
            |_,_,_|->Result<FluidTransportState,MaterialError>{panic!("preflight before material query")}).is_err());
    }
}

#[test]
fn depth_and_reversed_temperature_coordinates_preserve_physical_heating() {
    let gate=CancelGate::new();let mut a=slab2();let mut b=slab2();
    let r=a.step_with_material_transport(&power(),&[[0.0;2];16],0.01,limits(16),&gate,|_,t,_|Ok(material(1000.0+t,1.0))).unwrap();
    let s=b.step_with_material_transport(&power(),&[[0.0;2];16],0.02,limits(16),&gate,|_,t,_|Ok(material(1000.0+t,1.0))).unwrap();
    assert_eq!(a.lattice.gpop,b.lattice.gpop);near(s.heat.source_j,2.0*r.heat.source_j,1e-12);
    let mut inverted=plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.0).unwrap()
        .thermal_slab(MaterialThermalSlab {nx:2,ny:6,bottom_temperature_k:300.0,top_temperature_k:310.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16}).unwrap();
    let before=inverted.temperature_si(inverted.lattice.temperature(0,3)).unwrap();
    let r=inverted.step_with_material_transport(&power(),&[[0.0;2];16],0.01,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert!(r.heat.source_j>0.0);assert!(inverted.temperature_si(inverted.lattice.temperature(0,3)).unwrap()>before);
}
