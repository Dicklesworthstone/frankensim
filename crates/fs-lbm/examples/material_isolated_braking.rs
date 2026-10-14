//! Electrode-free magnetic braking with explicit synthetic material/field data.
use std::error::Error;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{Conductor,ElectricSetup};
use fs_lbm::material::{MaterialThermalSlab3,plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput,ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::MagneticMidpointOptions;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

fn main() -> Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().collect();
    if !(2..=5).contains(&args.len()) {
        return Err("usage: material_isolated_braking MATERIAL.txt [B_low_T] [B_high_T] [steps:1..200]".into());
    }
    let low:f64=args.get(2).map_or(Ok(20.0),|v|v.parse())?;
    let high:f64=args.get(3).map_or(Ok(60.0),|v|v.parse())?;
    let steps:usize=args.get(4).map_or(Ok(20),|v|v.parse())?;
    if !low.is_finite() || !high.is_finite() || !(1..=200).contains(&steps) {return Err("finite field bounds and 1..200 steps required".into());}
    let text=std::fs::read_to_string(&args[1])?;
    let values:Vec<f64>=text.lines().flat_map(|l|l.split('#').next().unwrap_or("").split_whitespace())
        .map(str::parse).collect::<Result<_,_>>()?;
    if values.len()!=7 || !values.iter().all(|v|v.is_finite()) || values[4]<=0.0 || values[6]<=0.0 {
        return Err("expected rho cp k mu resistivity alpha Tref; positive resistivity and reference Kelvin".into());
    }
    let material=FluidTransportState::try_new(values[0],values[1],values[2],values[3])?;
    let law=Conductor::new(values[4]).with_temperature_coefficient(values[5],values[6]);
    let scale=plan_scaling_from_material(&material,0.01,0.01,0.1)?;
    let mut slab=scale.thermal_slab3(MaterialThermalSlab3 {
        nx:4,ny:2,nz:6,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:64,
    })?;
    // Initial motion is a physical initial condition, not an electrical drive.
    for c in 8..56 {slab.lattice.grid.f[c]=fs_lbm::equilibrium3(1.0,[scale.velocity_lattice(0.02)?,0.0,0.0]);}
    // Bz depends on x only, so this imposed field has zero divergence.
    let field:Vec<_>=(0..64).map(|c|[0.0,0.0,low+(high-low)*(c%4) as f64/3.0]).collect();
    let setup=ElectricSetup::default(); // no electrodes, including no fake ground
    let input=ElectricHeatInput {setup:&setup,contacts:&[],limits:ElectricHeatLimits {
        thermal:HeatSourceLimits {max_cells:64,temperature_range_k:[290.0,320.0],balance_tolerance_j:1e-9},
        coupling_tolerance_j:1e-10,
    }};
    let options=MagneticMidpointOptions {max_iterations:128,temperature_tolerance_k:1e-9,velocity_tolerance_m_s:1e-11,
        joule_tolerance_j:1e-10,impulse_tolerance_n_s:1e-12,work_tolerance_j:1e-10,relaxation:0.5};
    let gate=CancelGate::new();
    println!("# synthetic material; unwired fluid; imposed field; viscous wall losses also present; no total-energy claim");
    println!("time_s,mean_vx_m_s,mean_K,joule_W,magnetic_work_W,electrode_W,trials,coupling_defect_J");
    for i in 0..steps {
        let result=slab.step_magnetothermal_midpoint_isolated(input,&field,&[[0.0;3];64],1e-10,options,&gate,
            |_,t,_|Ok(Some(law.resistivity(t))))?;
        let mut u=0.0;let mut t=0.0;
        for c in 8..56 {
            u+=scale.velocity_si(slab.lattice.grid.moments3(c).1[0])?/48.0;
            t+=slab.temperature_si(slab.lattice.gpop[c].iter().sum())?/48.0;
        }
        println!("{:.6},{:.12e},{:.12},{:.12e},{:.12e},{:.12e},{},{:.6e}",(i+1) as f64*scale.timestep_s(),
            u,t,result.motion.electric.total_joule_w,result.motion.mechanical_power_w,result.motion.electric.delivered_w,
            result.iterations,result.coupling_defect_j);
    }
    Ok(())
}
