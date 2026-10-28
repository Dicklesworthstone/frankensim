//! Unwired conductive fluids use the existing current, flow and thermal steps.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::{MagneticMidpointError, MagneticMidpointOptions, MagneticMidpointStep};
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

// Independent four-resistor loop: sum the prescribed motional voltages around
// the square. Each bulk half-cell contributes rho/(2*DX) ohm to each edge.
fn square_power(step: &MagneticMidpointStep, cells: [usize; 4], b: &[[f64; 3]], rho: f64) -> f64 {
    let [a, bcell, c, d] = cells;
    let mut emf = 0.0;
    for (p, q, axis, sign) in [(a,bcell,0,1.0),(bcell,d,1,1.0),(d,c,0,-1.0),(c,a,1,-1.0)] {
        let motional = |cell: usize| {
            let u = step.evaluation_velocity_m_s[cell];
            if axis == 0 { u[1] * b[cell][2] } else { -u[0] * b[cell][2] }
        };
        emf += sign * 0.5 * DX * (motional(p) + motional(q));
    }
    let resistance = 4.0 * rho / DX;
    let current = emf / resistance;
    for cell in cells {
        close(step.motion.electric.joule_w[cell], current.powi(2) * rho / DX, 1e-12);
    }
    emf.powi(2) / resistance
}

#[test]
fn insulating_slit_removes_magnetic_braking_and_heat_from_the_actual_midpoint_slab() {
    let setup = ElectricSetup::default();
    let gate = CancelGate::new();
    let b = field(8);
    let mut closed = slab2();
    let mut legacy = slab2();
    let mut open = slab2();
    let mut unforced = slab2();
    let r = closed.step_magnetothermal_midpoint_isolated_with_open_faces(
        input(&setup), &[], &b, &[[0.0;2];8], DX, 1e-10, options(), &gate,
        |_,_,_| Ok(Some(1.0)),
    ).unwrap();
    let old = legacy.step_magnetothermal_midpoint_isolated(
        input(&setup), &b, &[[0.0;2];8], DX, 1e-10, options(), &gate,
        |_,_,_| Ok(Some(1.0)),
    ).unwrap();
    assert_eq!(closed.lattice.grid.f, legacy.lattice.grid.f);
    assert_eq!(closed.lattice.gpop, legacy.lattice.gpop);
    assert_eq!(r.joule_energy_j, old.joule_energy_j);
    let expected = square_power(&r, [2,3,4,5], &b, 1.0);
    assert!(expected > 1e-9);
    close(r.joule_energy_j, expected * DT, 1e-12);
    close(r.midpoint_mechanical_work_j, -expected * DT, 1e-12);

    // Breaking one edge leaves a tree. Electrostatic charge separation cancels
    // ANY motional voltage on a tree, including the trial's changing velocity.
    let s = open.step_magnetothermal_midpoint_isolated_with_open_faces(
        input(&setup), &[[2,3]], &b, &[[0.0;2];8], DX, 1e-10, options(), &gate,
        |_,_,_| Ok(Some(1.0)),
    ).unwrap();
    let ordinary = unforced.step_with_power_density(
        &[0.0;8], &[[0.0;2];8], DX, input(&setup).limits.thermal,
    ).unwrap();
    assert!(s.iterations > 1); // Trial velocity still evolves by the real fluid step.
    assert!(s.joule_energy_j < 1e-20);
    assert!(s.motion.lorentz_force_density_n_m3.iter().flatten().all(|f| f.abs() < 1e-10));
    for (a,b) in open.lattice.grid.f.iter().flatten().zip(unforced.lattice.grid.f.iter().flatten()) { close(*a,*b,1e-14); }
    for (a,b) in open.lattice.gpop.iter().flatten().zip(unforced.lattice.gpop.iter().flatten()) { close(*a,*b,1e-14); }
    close(s.heat.stored_change_j, ordinary.stored_change_j, 1e-12);
    let speed = |slab: &MaterialSlab| (2..6).map(|c| slab.lattice.grid.moments(c).u[0]).sum::<f64>();
    assert!(speed(&closed) < speed(&open));
}

#[test]
fn isolated_laminations_keep_independent_currents_and_heat_during_three_dimensional_feedback() {
    let setup = ElectricSetup::default();
    let gate = CancelGate::new();
    let b = field(16);
    let mut slab = slab3();
    let openings = [[4,8],[5,9],[6,10],[7,11]];
    let result = slab.step_magnetothermal_midpoint_isolated_with_open_faces(
        input(&setup), &openings, &b, &[[0.0;3];16], 1e-10, options(), &gate,
        |c,_,_| Ok(Some(if c < 8 {1.0} else {2.0})),
    ).unwrap();
    assert!(result.iterations > 1);
    // Two electrically independent sheets, with all mechanical/thermal cells
    // still present. Gauge cells retain their resistance, heating and force.
    close(result.motion.electric.potential_v[4], 0.0, 0.0);
    close(result.motion.electric.potential_v[8], 0.0, 0.0);
    let power = square_power(&result, [4,5,6,7], &b, 1.0)
        + square_power(&result, [8,9,10,11], &b, 2.0);
    assert!(result.motion.electric.joule_w[4] > 0.0);
    assert!(result.motion.electric.joule_w[8] > 0.0);
    close(result.joule_energy_j, power * DT, 1e-12);
    close(result.midpoint_mechanical_work_j, -power * DT, 1e-12);
    close(result.heat.source_j, power * DT, 1e-12);
    assert!(result.source_residual.impulse_n_s <= options().impulse_tolerance_n_s);
    assert_eq!(slab.lattice.grid.f.len(), 16);
    assert_eq!(slab.lattice.gpop.len(), 16);
}

#[test]
fn open_midpoint_topology_refusals_preserve_the_complete_accepted_flow() {
    let setup = ElectricSetup::default();
    let mut slab = slab2();
    let before = (slab.lattice.grid.f.clone(), slab.lattice.gpop.clone(),
        slab.lattice.grid.fext.clone(), slab.lattice.grid.tau.clone());
    for openings in [vec![[2,5]], vec![[2,3],[3,2]]] {
        assert!(slab.step_magnetothermal_midpoint_isolated_with_open_faces(
            input(&setup), &openings, &field(8), &[[0.0;2];8], DX, 1e-10, options(),
            &CancelGate::new(), |_,_,_| Ok(Some(1.0)),
        ).is_err());
    }
    let gate = CancelGate::new();
    assert!(slab.step_magnetothermal_midpoint_isolated_with_open_faces(
        input(&setup), &[[2,3]], &field(8), &[[0.0;2];8], DX, 1e-10, options(), &gate,
        |_,_,_| { gate.request(); Ok(Some(1.0)) },
    ).is_err());
    assert_eq!(slab.lattice.grid.f, before.0);
    assert_eq!(slab.lattice.gpop, before.1);
    assert_eq!(slab.lattice.grid.fext, before.2);
    assert_eq!(slab.lattice.grid.tau, before.3);
}
