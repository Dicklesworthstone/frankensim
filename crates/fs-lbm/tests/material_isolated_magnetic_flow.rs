//! Unwired conductive fluids use the existing current, flow and thermal steps.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::{MagneticMidpointError, MagneticMidpointOptions};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;
const DX:f64=0.01;
const DT:f64=0.01;
fn material() -> FluidTransportState {FluidTransportState::try_new(1000.0,1000.0,100.0,1.0).unwrap()}
fn slab3() -> MaterialSlab3 {
    let mut slab=plan_scaling_from_material(&material(),DX,DT,0.1).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx:2,ny:2,nz:2,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16,
    }).unwrap();
    for c in 4..12 {slab.lattice.grid.f[c]=fs_lbm::equilibrium3(1.0,[0.02,0.0,0.0]);}
    slab
}
fn slab2() -> MaterialSlab {
    let mut slab=plan_scaling_from_material(&material(),DX,DT,0.1).unwrap().thermal_slab(MaterialThermalSlab {
        nx:2,ny:2,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:8,
    }).unwrap();
    for c in 2..6 {slab.lattice.grid.f[c]=fs_lbm::equilibrium(1.0,0.02,0.0);}
    slab
}
fn field(n:usize) -> Vec<[f64;3]> {(0..n).map(|c|[0.0,0.0,if c%2==0 {20.0}else{60.0}]).collect()}
fn input(setup:&ElectricSetup) -> ElectricHeatInput<'_> {ElectricHeatInput {
    setup,contacts:&[],limits:ElectricHeatLimits {
        thermal:HeatSourceLimits {max_cells:16,temperature_range_k:[290.0,320.0],balance_tolerance_j:1e-9},
        coupling_tolerance_j:1e-10,
    },
}}
fn options() -> MagneticMidpointOptions {MagneticMidpointOptions {
    max_iterations:128,temperature_tolerance_k:1e-9,velocity_tolerance_m_s:1e-11,
    joule_tolerance_j:1e-10,impulse_tolerance_n_s:1e-12,work_tolerance_j:1e-10,relaxation:0.5,
}}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} != {b}");}

#[test]
fn unwired_three_dimensional_flow_receives_the_actual_heat_and_reaction() {
    let setup=ElectricSetup::default();let mut coupled=slab3();let mut direct=slab3();
    let external=[[0.0;3];16];let gate=CancelGate::new();
    let r=coupled.step_magnetothermal_isolated(input(&setup),&field(16),&external,1e-10,&gate,
        |_,_,_|Ok(Some(1.0))).unwrap();
    assert!(r.motion.electric.electrodes.is_empty());assert_eq!(r.electrode_energy_j,0.0);
    assert!(r.joule_energy_j>0.0 && r.lagged_mechanical_work_j<0.0);
    close(r.joule_energy_j,-r.lagged_mechanical_work_j,1e-10);
    // No dummy electrodes: all eight fluid cells, including the gauge, heat.
    assert!(r.motion.electric.joule_w[4..12].iter().all(|q|*q>0.0));
    let power:Vec<_>=r.motion.electric.joule_w.iter().map(|q|q/DX.powi(3)).collect();
    let heat=direct.step_with_power_density(&power,&r.motion.lorentz_force_density_n_m3,input(&setup).limits.thermal).unwrap();
    assert_eq!(heat,r.heat);assert_eq!(coupled.lattice.grid.f,direct.lattice.grid.f);assert_eq!(coupled.lattice.gpop,direct.lattice.gpop);
    assert!(r.coupling_defect_j.abs()<1e-9);
}

#[test]
fn two_dimensional_depth_scales_current_energy_not_temperature_or_velocity() {
    let setup=ElectricSetup::default();let mut a=slab2();let mut b=slab2();let gate=CancelGate::new();
    let r=a.step_magnetothermal_isolated(input(&setup),&field(8),&[[0.0;2];8],DX,1e-10,&gate,
        |_,_,_|Ok(Some(1.0))).unwrap();
    let s=b.step_magnetothermal_isolated(input(&setup),&field(8),&[[0.0;2];8],2.0*DX,1e-10,&gate,
        |_,_,_|Ok(Some(1.0))).unwrap();
    close(s.joule_energy_j,2.0*r.joule_energy_j,1e-12);
    for (a,b) in a.lattice.grid.f.iter().flatten().zip(b.lattice.grid.f.iter().flatten()) {close(*a,*b,1e-12);}
    for (a,b) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) {close(*a,*b,1e-12);}
}

#[test]
fn isolated_midpoint_advances_once_and_checks_actual_material_and_velocity_feedback() {
    let setup=ElectricSetup::default();let gate=CancelGate::new();let mut a=slab3();
    let old_t:Vec<_>=a.lattice.gpop.iter().map(|p|a.temperature_si(p.iter().sum()).unwrap()).collect();
    let old_u:Vec<_>=(0..16).map(|c|a.lattice.grid.moments3(c).1).collect();
    let mut direct=slab3();let o=options();
    let r=a.step_magnetothermal_midpoint_isolated(input(&setup),&field(16),&[[0.0;3];16],1e-10,o,&gate,
        |_,t,transport| {assert_eq!(*transport,material());Ok(Some(1.0+0.01*(t-300.0)))}).unwrap();
    assert!(r.iterations>1 && r.electrical_solves>r.iterations);
    assert_eq!(r.electrode_energy_j,0.0);assert!(r.joule_energy_j>0.0);
    assert!(r.source_residual.impulse_n_s<=o.impulse_tolerance_n_s);
    for c in 4..12 {
        close(r.evaluation_temperature_k[c],0.5*(old_t[c]+a.temperature_si(a.lattice.gpop[c].iter().sum()).unwrap()),o.temperature_tolerance_k);
        let end=a.lattice.grid.moments3(c).1;
        for k in 0..3 {close(r.evaluation_velocity_m_s[c][k],0.5*(old_u[c][k]+end[k]),o.velocity_tolerance_m_s);}
    }
    let power:Vec<_>=r.motion.electric.joule_w.iter().map(|q|q/DX.powi(3)).collect();
    let heat=direct.step_with_power_density(&power,&r.motion.lorentz_force_density_n_m3,input(&setup).limits.thermal).unwrap();
    assert_eq!(heat,r.heat);assert_eq!(a.lattice.grid.f,direct.lattice.grid.f);assert_eq!(a.lattice.gpop,direct.lattice.gpop);
    // The 2-D midpoint entry point uses the same machinery, not a dead API.
    let mut two=slab2();
    let r=two.step_magnetothermal_midpoint_isolated(input(&setup),&field(8),&[[0.0;2];8],DX,1e-10,o,&gate,
        |_,_,_|Ok(Some(1.0))).unwrap();
    assert!(r.joule_energy_j>0.0);assert!(r.coupling_defect_j.abs()<1e-9);
}

#[test]
fn isolated_flow_errors_never_publish_trial_populations_or_forces() {
    let setup=ElectricSetup::default();let mut a=slab3();let gate=CancelGate::new();
    let before_f=a.lattice.grid.f.clone();let before_g=a.lattice.gpop.clone();let before_force=a.lattice.grid.fext.clone();
    let mut wired=setup.clone();wired.electrodes.push(Electrode {lo:[0,0,1],hi:[1,1,2],drive:ElectrodeDrive::Voltage(0.0)});
    assert!(a.step_magnetothermal_isolated(input(&wired),&field(16),&[[0.0;3];16],1e-10,&gate,
        |_,_,_|panic!("no electrodes in this path")).is_err());
    assert!(a.step_magnetothermal_isolated(input(&setup),&field(16),&[[0.0;3];16],-1.0,&gate,
        |_,_,_|panic!("current tolerance before queries")).is_err());
    assert!(a.step_magnetothermal_isolated(input(&setup),&field(16),&[[0.0;3];16],1e-10,&gate,
        |_,_,_|Err(ChtError::InvalidInput {field:"source",reason:"missing resistivity".into()})).is_err());
    let e=a.step_magnetothermal_midpoint_isolated(input(&setup),&field(16),&[[0.0;3];16],1e-10,
        MagneticMidpointOptions {max_iterations:1,..options()},&gate,|_,_,_|Ok(Some(1.0))).unwrap_err();
    assert!(matches!(e,MagneticMidpointError::NonConvergence {..}));
    assert!(a.step_magnetothermal_isolated(input(&setup),&field(16),&[[0.0;3];16],1e-10,&gate,
        |_,_,_| {gate.request();Ok(Some(1.0))}).is_err());
    assert_eq!(a.lattice.grid.f,before_f);assert_eq!(a.lattice.gpop,before_g);assert_eq!(a.lattice.grid.fext,before_force);
}
