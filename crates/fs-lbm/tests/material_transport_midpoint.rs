//! Real current/local-transport composition; all coefficients are synthetic.
use std::cell::Cell;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::MagneticMidpointOptions;
use fs_lbm::material::electrothermal::magnetic::transport::ElectricalTransport;
use fs_lbm::material::electrothermal::magnetic::transport::midpoint::{MaterialMidpointError, MaterialMidpointOptions};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};

const DX: f64 = 0.01;
const N: usize = 16;
fn state(k: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, k, mu).unwrap()
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&state(100.0, 1.0), DX, 0.01, 0.05).unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: N,
        }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&state(100.0, 1.0), DX, 0.01, 0.05).unwrap()
        .thermal_slab(MaterialThermalSlab {
            nx: 2, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: N,
        }).unwrap()
}
fn setup(normal: usize, volts: f64) -> ElectricSetup {
    let (lo, hi, end_lo, end_hi) = if normal == 2 {
        ([0,0,1], [2,1,2], [0,0,6], [2,1,7])
    } else { ([0,1,0], [2,2,1], [0,6,0], [2,7,1]) };
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo, hi, drive: ElectrodeDrive::Voltage(volts) },
        Electrode { lo: end_lo, hi: end_hi, drive: ElectrodeDrive::Voltage(0.0) },
    ] }
}
fn input(setup: &ElectricSetup) -> ElectricHeatInput<'_> {
    ElectricHeatInput { setup, contacts: &[], limits: ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells: N, temperature_range_k: [290.0, 330.0], balance_tolerance_j: 1e-8 },
        coupling_tolerance_j: 1e-9,
    } }
}
fn options() -> MaterialMidpointOptions {
    MaterialMidpointOptions { coupling: MagneticMidpointOptions {
        max_iterations: 128, temperature_tolerance_k: 1e-8, velocity_tolerance_m_s: 1e-10,
        joule_tolerance_j: 1e-8, impulse_tolerance_n_s: 1e-10, work_tolerance_j: 1e-9,
        relaxation: 0.5,
    }, transport_relative_tolerance: 1e-8 }
}
fn law(_: usize, t: f64, _: &FluidTransportState) -> Result<ElectricalTransport<FluidTransportState>, MaterialError> {
    ElectricalTransport::try_new(state(100.0 * (1.0 + 0.02 * (t-300.0)),
        1.0 / (1.0 + 0.03 * (t-300.0))), Some(1.0 + 0.01 * (t-300.0)))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a} vs {b}, tolerance {tolerance}");
}

#[test]
fn constant_material_preserves_the_original_midpoint_fields_and_trial_count() {
    let mut a=slab3(); let mut b=slab3(); let drive=setup(2, 10.0); let gate=CancelGate::new();
    let field=[[0.0, 2.0, 0.0]; N]; let forces=[[1.0, 0.0, 0.0]; N];
    let old=a.step_magnetothermal_midpoint(input(&drive), &field, &forces, options().coupling,
        &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let new=b.step_magnetothermal_midpoint_with_material_transport(input(&drive), &field, &forces,
        options(), &gate, |_,_,_| ElectricalTransport::try_new(state(100.0,1.0),Some(1.0))).unwrap();
    assert_eq!(old.heat,new.step.heat); assert_eq!(old.motion.electric.joule_w,new.step.motion.electric.joule_w);
    assert_eq!(old.iterations,new.step.iterations); assert_eq!(old.electrical_solves,new.step.electrical_solves);
    assert_eq!(a.lattice.grid.f,b.lattice.grid.f); assert_eq!(a.lattice.gpop,b.lattice.gpop);
    assert_eq!(new.transport_residual.conductivity_relative,0.0);
    assert_eq!(new.transport_residual.viscosity_relative,0.0);
}

#[test]
fn applied_midpoint_objects_drive_the_actual_current_force_and_transport_once() {
    let mut slab=slab3(); let mut replay=slab3(); let drive=setup(2, 30.0); let gate=CancelGate::new();
    let contacts=[ElectricalContact { cells:[6,8],resistance_ohm_m2:0.01,heat_fraction_to_first:0.7 },
                  ElectricalContact { cells:[7,9],resistance_ohm_m2:0.01,heat_fraction_to_first:0.7 }];
    let request=ElectricHeatInput { contacts:&contacts,..input(&drive) };
    let old:Vec<_>=slab.lattice.gpop.iter().map(|p|slab.temperature_si(p.iter().sum()).unwrap()).collect();
    let forces=[[1.0,0.0,0.0];N];
    let r=slab.step_magnetothermal_midpoint_with_material_transport(request,&[[0.0,2.0,0.0];N],&forces,
        options(),&gate,law).unwrap();
    assert!(r.step.iterations>1);
    for (a,e) in r.applied.iter().zip(&r.endpoint) {
        assert_eq!(a.temperature_k,r.step.evaluation_temperature_k[a.cell]);
        assert_eq!(a.material,law(a.cell,a.temperature_k,&state(100.0,1.0)).unwrap());
        let final_t=slab.temperature_si(slab.lattice.gpop[a.cell].iter().sum()).unwrap();
        close(a.temperature_k,0.5*old[a.cell]+0.5*final_t,options().coupling.temperature_tolerance_k);
        assert_eq!(e.temperature_k,final_t);
    }
    let power:Vec<_>=r.step.motion.electric.joule_w.iter().map(|w|w/DX.powi(3)).collect();
    let total:Vec<[f64;3]>=r.step.motion.lorentz_force_density_n_m3.iter().zip(forces)
        .map(|(m,f)|core::array::from_fn(|a|m[a]+f[a])).collect();
    // Replay only one physical step with the accepted solved coefficient field.
    let h=replay.step_with_material_transport(&power,&total,request.limits.thermal,&gate,
        |cell,_,_| Ok(r.applied[cell-2].transport)).unwrap();
    assert_eq!(h.heat,r.step.heat); assert_eq!(replay.lattice.grid.f,slab.lattice.grid.f);
    assert_eq!(replay.lattice.gpop,slab.lattice.gpop);
    assert!(r.transport_residual.conductivity_relative<=options().transport_relative_tolerance);
    assert!(r.transport_residual.viscosity_relative<=options().transport_relative_tolerance);
}

#[test]
fn zero_electrical_residual_cannot_bypass_unconverged_material_transport() {
    let mut slab=slab3(); let drive=setup(2,0.0); let gate=CancelGate::new();
    // A non-affine initial temperature makes a zero-current transport trial
    // move temperature even before any nonequilibrium heat flux has developed.
    for value in &mut slab.lattice.gpop[6] { *value *= 1.2; }
    let before=slab.lattice.gpop.clone(); let old_f=slab.lattice.grid.f.clone();
    let loose=MaterialMidpointOptions { coupling:MagneticMidpointOptions {
        max_iterations:1,temperature_tolerance_k:100.0,velocity_tolerance_m_s:1.0,..options().coupling
    },transport_relative_tolerance:0.0 };
    let result=slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
        loose,&gate,law);
    match result {
        Err(MaterialMidpointError::NonConvergence {last_checked_transport:Some(r),..}) =>
            assert!(r.conductivity_relative>0.0 || r.viscosity_relative>0.0),
        other=>panic!("transport gate must reject the first coefficient trial: {other:?}"),
    }
    assert_eq!(slab.lattice.gpop,before); assert_eq!(slab.lattice.grid.f,old_f);
    let result=slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
        options(),&gate,law).unwrap();
    assert!(result.step.iterations>1); assert_eq!(result.step.source_residual.joule_j,0.0);
}

#[test]
fn non_clone_material_objects_are_extracted_once_and_retained_at_both_epochs() {
    #[derive(Debug)] struct Source { t:f64, reads:Cell<usize> }
    impl FluidTransport for Source {
        fn transport_state(&self)->Result<FluidTransportState,MaterialError> {
            assert_eq!(self.reads.replace(1),0,"re-extracted one source object");
            Ok(state(100.0,1.0))
        }
    }
    let mut slab=slab3();let drive=setup(2,10.0);let mut queries=0;
    let r=slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
        options(),&CancelGate::new(),|_,t,_| {
            queries+=1;ElectricalTransport::try_new(Source{t,reads:Cell::new(0)},Some(1.0))
        }).unwrap();
    assert_eq!(queries,(r.step.electrical_solves+1)*12);
    assert_eq!(r.applied.len(),12);assert_eq!(r.endpoint.len(),12);
    for sample in r.applied.iter().chain(&r.endpoint) {
        assert_eq!(sample.material.material().t,sample.temperature_k);
        assert_eq!(sample.material.material().reads.get(),1);
    }
}

#[test]
fn failed_endpoint_and_endpoint_cancellation_do_not_publish_the_converged_candidate() {
    let drive=setup(2,10.0);
    let loose=MaterialMidpointOptions { coupling:MagneticMidpointOptions {
        max_iterations:1,temperature_tolerance_k:100.0,velocity_tolerance_m_s:1.0,
        joule_tolerance_j:1.0,impulse_tolerance_n_s:1.0,work_tolerance_j:1.0,..options().coupling
    },transport_relative_tolerance:1.0 };
    for cancel in [false,true] {
        let mut slab=slab3();let gate=CancelGate::new();let mut calls=0;
        let f=slab.lattice.grid.f.clone();let g=slab.lattice.gpop.clone();
        let force=slab.lattice.grid.fext.clone();let tau=slab.lattice.grid.tau.clone();
        let result=slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
            loose,&gate,|_,_,_| {
                calls+=1;
                if calls==25 { // 12 applied + 12 recheck; first endpoint query
                    if cancel {gate.request();} else {return Err(MaterialError::State{what:"endpoint source refusal".into()});}
                }
                ElectricalTransport::try_new(state(100.0,1.0),Some(1.0))
            });
        assert!(result.is_err());assert!(calls>=25);
        assert_eq!(slab.lattice.grid.f,f);assert_eq!(slab.lattice.gpop,g);
        assert_eq!(slab.lattice.grid.fext,force);assert_eq!(slab.lattice.grid.tau,tau);
    }
}

#[test]
fn both_dimensions_support_isolated_midpoints_and_2d_depth_scales_extensive_heat() {
    let drive=setup(1,10.0);let gate=CancelGate::new();let mut a=slab2();let mut b=slab2();
    let ar=a.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0,0.0,2.0];N],&[[0.0;2];N],
        DX,options(),&gate,law).unwrap();
    let br=b.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0,0.0,2.0];N],&[[0.0;2];N],
        2.0*DX,options(),&gate,law).unwrap();
    close(br.step.joule_energy_j,2.0*ar.step.joule_energy_j,1e-9);
    for (x,y) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) {close(*x,*y,1e-9);}
    let isolated=ElectricSetup::default();let mut c=slab3();
    for cell in 2..14 {c.lattice.grid.f[cell]=fs_lbm::equilibrium3(1.0,[0.01,0.0,0.0]);}
    let varying:Vec<_>=(0..N).map(|cell|[0.0,if cell%2==0 {10.0}else{30.0},0.0]).collect();
    // By varies with x; divergence is zero, while motional voltage varies.
    let r=c.step_magnetothermal_midpoint_isolated_with_material_transport(input(&isolated),&varying,&[[0.0;3];N],
        1e-10,options(),&gate,law).unwrap();
    assert!(r.step.motion.electric.electrodes.is_empty());assert_eq!(r.step.electrode_energy_j,0.0);
    assert!(r.step.joule_energy_j>0.0);assert!(r.step.midpoint_mechanical_work_j<0.0);
    close(r.step.joule_energy_j,-r.step.midpoint_mechanical_work_j,1e-9);
    let r2=a.step_magnetothermal_midpoint_isolated_with_material_transport(input(&isolated),&[[0.0,0.0,2.0];N],&[[0.0;2];N],
        DX,1e-10,options(),&gate,law).unwrap();
    assert_eq!(r2.step.electrode_energy_j,0.0);
}

#[test]
fn invalid_controls_and_storage_changes_refuse_without_overwriting_history() {
    let drive=setup(2,0.0);let mut slab=slab3();let before=slab.lattice.gpop.clone();
    for tolerance in [-1.0,f64::NAN,f64::INFINITY,2.0] {
        let opts=MaterialMidpointOptions{transport_relative_tolerance:tolerance,..options()};
        let r=slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
            opts,&CancelGate::new(),|_,_,_|->Result<ElectricalTransport<FluidTransportState>,MaterialError>{panic!("validate controls first")});
        assert!(matches!(r,Err(MaterialMidpointError::InvalidTransportTolerance)));
    }
    assert!(slab.step_magnetothermal_midpoint_with_material_transport(input(&drive),&[[0.0;3];N],&[[0.0;3];N],
        options(),&CancelGate::new(),|_,_,_|ElectricalTransport::try_new(
            FluidTransportState::try_new(900.0,1000.0,100.0,1.0)?,Some(1.0))).is_err());
    assert_eq!(slab.lattice.gpop,before);
}
