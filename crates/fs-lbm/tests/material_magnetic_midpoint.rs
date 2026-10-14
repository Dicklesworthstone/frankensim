//! Synthetic G0/G1/G3/G4 regressions through the actual current and flow kernels.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricMotionInput, ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatError, ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::{MagneticMidpointError, MagneticMidpointOptions};
use fs_lbm::material::electrothermal::midpoint::MidpointOptions;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

const DX: f64 = 0.01;
const DT: f64 = 0.01;
fn transport() -> FluidTransportState { FluidTransportState::try_new(1000.0, 1000.0, 100.0, 1.0).unwrap() }
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&transport(), DX, DT, 0.1).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&transport(), DX, DT, 0.1).unwrap().thermal_slab(MaterialThermalSlab {
        nx: 2, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn electrodes(normal: usize, voltage: f64) -> ElectricSetup {
    let (a,b,c,d) = if normal == 2 { ([0,0,1],[2,1,2],[0,0,6],[2,1,7]) }
        else { ([0,1,0],[2,2,1],[0,6,0],[2,7,1]) };
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo:a, hi:b, drive:ElectrodeDrive::Voltage(voltage) },
        Electrode { lo:c, hi:d, drive:ElectrodeDrive::Voltage(0.0) },
    ] }
}
fn input(setup: &ElectricSetup) -> ElectricHeatInput<'_> { ElectricHeatInput {
    setup, contacts: &[], limits: ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells:16, temperature_range_k:[290.0,320.0], balance_tolerance_j:1e-9 },
        coupling_tolerance_j:1e-10,
    },
} }
fn options() -> MagneticMidpointOptions { MagneticMidpointOptions {
    max_iterations:128, temperature_tolerance_k:1e-9, velocity_tolerance_m_s:1e-11,
    joule_tolerance_j:1e-10, impulse_tolerance_n_s:1e-12, work_tolerance_j:1e-10, relaxation:0.5,
} }
fn law(_:usize, t:f64, _: &FluidTransportState) -> Result<Option<f64>,ChtError> { Ok(Some(1.0+0.02*(t-300.0))) }
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a} != {b}"); }

#[test]
fn zero_field_retains_the_existing_temperature_midpoint_equation() {
    let setup=electrodes(2,20.0); let gate=CancelGate::new();
    let mut a=slab3(); let mut b=slab3(); let forces=[[1.0,0.0,0.0];16];
    let o=options();
    let old=a.step_electrothermal_midpoint(input(&setup),&forces,MidpointOptions {
        max_iterations:o.max_iterations,temperature_tolerance_k:o.temperature_tolerance_k,
        source_tolerance_j:o.joule_tolerance_j,relaxation:o.relaxation,
    },&gate,law).unwrap();
    let new=b.step_magnetothermal_midpoint(input(&setup),&[[0.0;3];16],&forces,o,&gate,law).unwrap();
    for (x,y) in a.lattice.grid.f.iter().flatten().zip(b.lattice.grid.f.iter().flatten()) { near(*x,*y,1e-11); }
    for (x,y) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) { near(*x,*y,1e-11); }
    near(new.joule_energy_j,old.step.joule_energy_j,1e-10);
    near(new.midpoint_mechanical_work_j,0.0,0.0);
    near(new.source_residual.impulse_n_s,0.0,0.0);
}

#[test]
fn accepted_sources_close_at_actual_temperature_and_velocity_midpoints() {
    let mut slab=slab3(); let setup=electrodes(2,2.0); let gate=CancelGate::new();
    let t0:Vec<_>=slab.lattice.gpop.iter().map(|p|slab.temperature_si(p.iter().sum()).unwrap()).collect();
    let u0:Vec<_>=(0..16).map(|c|slab.lattice.grid.moments3(c).1).collect(); // dx/dt = 1 m/s
    let b=[[0.0,200.0,0.0];16]; let o=options();
    let result=slab.step_magnetothermal_midpoint(input(&setup),&b,&[[0.0;3];16],o,&gate,law).unwrap();
    let tm:Vec<_>=(0..16).map(|c|0.5*t0[c]+0.5*slab.temperature_si(slab.lattice.gpop[c].iter().sum()).unwrap()).collect();
    let um:Vec<[f64;3]>=(0..16).map(|c|core::array::from_fn(|a|0.5*u0[c][a]+0.5*slab.lattice.grid.moments3(c).1[a])).collect();
    for c in 0..16 {
        near(result.evaluation_temperature_k[c],tm[c],o.temperature_tolerance_k);
        for a in 0..3 { near(result.evaluation_velocity_m_s[c][a],um[c][a],o.velocity_tolerance_m_s); }
    }
    let domain=VoxelDomain::new(2,1,8,DX).unwrap();
    let checked=setup.solve_with_cell_motion(&domain,&tm,&[],ElectricMotionInput {
        velocity_m_s:&um,magnetic_flux_density_t:&b,max_cells:16,power_tolerance_w:1e-8,
    },|c,_,t| if !(2..14).contains(&c) {Ok(None)} else {law(c,t,&transport())},&gate).unwrap();
    let q:f64=result.motion.electric.joule_w.iter().zip(&checked.electric.joule_w).map(|(a,b)|(a-b).abs()*DT).sum();
    let impulse:f64=result.motion.lorentz_force_density_n_m3.iter().flatten()
        .zip(checked.lorentz_force_density_n_m3.iter().flatten()).map(|(a,b)|(a-b).abs()*DX.powi(3)*DT).sum();
    near(q,result.source_residual.joule_j,1e-14); near(impulse,result.source_residual.impulse_n_s,1e-14);
    assert!(q<=o.joule_tolerance_j && impulse<=o.impulse_tolerance_n_s);
    assert!(result.source_residual.mechanical_work_j<=o.work_tolerance_j);
    assert!(result.iterations>1 && result.electrical_solves>result.iterations);
    assert!(result.electrical_solves<=2*result.iterations);
    assert!(result.coupling_defect_j.abs()<1e-9);
}

#[test]
fn under_relaxed_braking_avoids_this_explicit_source_lag_reversal() {
    let setup=electrodes(2,0.0); let mut a=slab3(); let mut explicit=slab3(); let gate=CancelGate::new();
    for c in 2..14 {
        a.lattice.grid.f[c]=fs_lbm::equilibrium3(1.0,[0.02,0.0,0.0]);
        explicit.lattice.grid.f[c]=a.lattice.grid.f[c];
    }
    let b=[[0.0,350.0,0.0];16];
    explicit.step_magnetothermal(input(&setup),&b,&[[0.0;3];16],&gate,|_,_,_|Ok(Some(1.0))).unwrap();
    let r=a.step_magnetothermal_midpoint(input(&setup),&b,&[[0.0;3];16],options(),&gate,|_,_,_|Ok(Some(1.0))).unwrap();
    assert!(explicit.lattice.grid.moments3(6).1[0]<0.0);
    assert!(a.lattice.grid.moments3(6).1[0]>0.0);
    assert!(a.lattice.grid.moments3(6).1[0]<0.02);
    assert!(r.joule_energy_j>0.0 && r.midpoint_mechanical_work_j<0.0);
    near(r.electrode_energy_j,0.0,1e-12);
    near(r.joule_energy_j,-r.midpoint_mechanical_work_j,1e-10);
}

#[test]
fn force_and_work_rechecks_cannot_be_bypassed_by_loose_coordinate_tolerances() {
    let setup=electrodes(2,2.0); let mut slab=slab3();
    let before_f=slab.lattice.grid.f.clone(); let before_g=slab.lattice.gpop.clone();
    for gate_force in [true,false] {
        let o=MagneticMidpointOptions {max_iterations:1,temperature_tolerance_k:100.0,velocity_tolerance_m_s:1.0,
            joule_tolerance_j:1.0,impulse_tolerance_n_s:if gate_force {0.0}else{1.0},
            work_tolerance_j:if gate_force {1.0}else{0.0},..options()};
        let error=slab.step_magnetothermal_midpoint(input(&setup),&[[0.0,20.0,0.0];16],&[[0.0;3];16],o,
            &CancelGate::new(),|_,_,_|Ok(Some(1.0))).unwrap_err();
        let MagneticMidpointError::NonConvergence {source_residual:Some(r),..}=error else {panic!("expected independently rechecked source refusal")};
        assert!(if gate_force {r.impulse_n_s>0.0}else{r.mechanical_work_j>0.0});
        assert_eq!(slab.lattice.grid.f,before_f); assert_eq!(slab.lattice.gpop,before_g);
    }
}

#[test]
fn failed_and_cancelled_trials_preserve_all_accepted_state() {
    let setup=electrodes(2,2.0); let mut slab=slab3();
    let f=slab.lattice.grid.f.clone(); let g=slab.lattice.gpop.clone();
    let forces=slab.lattice.grid.fext.clone(); let tau=slab.lattice.grid.tau.clone(); let tg=slab.lattice.tau_g;
    let err=slab.step_magnetothermal_midpoint(input(&setup),&[[0.0,20.0,0.0];16],&[[0.0;3];16],
        MagneticMidpointOptions {max_iterations:1,..options()},&CancelGate::new(),law).unwrap_err();
    assert!(matches!(err,MagneticMidpointError::NonConvergence {..}));
    let gate=CancelGate::new();
    let err=slab.step_magnetothermal_midpoint(input(&setup),&[[0.0,20.0,0.0];16],&[[0.0;3];16],options(),&gate,
        |_,_,_| {gate.request();Ok(Some(1.0))}).unwrap_err();
    assert!(matches!(err,MagneticMidpointError::Step(ElectricHeatError::Cancelled)));
    assert_eq!(slab.lattice.grid.f,f);assert_eq!(slab.lattice.gpop,g);assert_eq!(slab.lattice.grid.fext,forces);
    assert_eq!(slab.lattice.grid.tau,tau);assert_eq!(slab.lattice.tau_g,tg);
}

#[test]
fn two_dimensional_depth_and_field_reversal_keep_the_same_coupled_equations() {
    let setup=electrodes(1,2.0); let mut a=slab2();let mut b=slab2();let mut reverse=slab2();let gate=CancelGate::new();
    let r=a.step_magnetothermal_midpoint(input(&setup),&[[0.0,0.0,20.0];16],&[[0.0;2];16],DX,options(),&gate,law).unwrap();
    let d=b.step_magnetothermal_midpoint(input(&setup),&[[0.0,0.0,20.0];16],&[[0.0;2];16],2.0*DX,options(),&gate,law).unwrap();
    let v=reverse.step_magnetothermal_midpoint(input(&setup),&[[0.0,0.0,-20.0];16],&[[0.0;2];16],DX,options(),&gate,law).unwrap();
    near(d.joule_energy_j,2.0*r.joule_energy_j,1e-10);near(v.joule_energy_j,r.joule_energy_j,1e-10);
    for c in 2..14 {
        near(a.lattice.grid.moments(c).u[0],b.lattice.grid.moments(c).u[0],1e-10);
        near(a.lattice.grid.moments(c).u[0],-reverse.lattice.grid.moments(c).u[0],1e-10);
    }
    assert!(a.step_magnetothermal_midpoint(input(&setup),&[[1.0,0.0,0.0];16],&[[0.0;2];16],DX,options(),&gate,
        |_,_,_|panic!("unsupported field before material queries")).is_err());
}

#[test]
fn material_only_resistivity_changes_drive_different_coupled_motion() {
    let setup=electrodes(2,0.2);let mut a=slab3();let mut b=slab3();let gate=CancelGate::new();
    let field=[[0.0,200.0,0.0];16];
    let r=a.step_magnetothermal_midpoint(input(&setup),&field,&[[0.0;3];16],options(),&gate,|_,_,_|Ok(Some(1.0))).unwrap();
    let s=b.step_magnetothermal_midpoint(input(&setup),&field,&[[0.0;3];16],options(),&gate,|_,_,_|Ok(Some(2.0))).unwrap();
    assert!(a.lattice.grid.moments3(6).1[0].abs()>b.lattice.grid.moments3(6).1[0].abs());
    assert!(r.joule_energy_j>s.joule_energy_j);
    assert!(r.source_residual.impulse_n_s<=options().impulse_tolerance_n_s);
    assert!(s.source_residual.impulse_n_s<=options().impulse_tolerance_n_s);
}

#[test]
fn invalid_controls_refuse_before_constitutive_queries() {
    let setup=electrodes(2,2.0);let mut slab=slab3();
    for o in [MagneticMidpointOptions {max_iterations:0,..options()},
        MagneticMidpointOptions {max_iterations:129,..options()},
        MagneticMidpointOptions {velocity_tolerance_m_s:f64::NAN,..options()},
        MagneticMidpointOptions {impulse_tolerance_n_s:-1.0,..options()},
        MagneticMidpointOptions {work_tolerance_j:f64::INFINITY,..options()},
        MagneticMidpointOptions {relaxation:0.0,..options()}] {
        let r=slab.step_magnetothermal_midpoint(input(&setup),&[[0.0,20.0,0.0];16],&[[0.0;3];16],o,
            &CancelGate::new(),|_,_,_|panic!("invalid controls before source query"));
        assert!(matches!(r,Err(MagneticMidpointError::InvalidOptions)));
    }
}
