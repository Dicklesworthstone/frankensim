//! Synthetic material substitution through electric current, magnetic force and heat.
use std::error::Error;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{Conductor, ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if !(2..=4).contains(&args.len()) {
        return Err("usage: material_magnetic_flow MATERIAL.txt [B_y_tesla] [steps: 1..500]".into());
    }
    let b: f64 = args.get(2).map_or(Ok(2.0), |s| s.parse())?;
    let steps: usize = args.get(3).map_or(Ok(100), |s| s.parse())?;
    if !b.is_finite() || !(1..=500).contains(&steps) { return Err("finite B and 1..500 steps required".into()); }
    // Same seven-value material format as material_joule_flow, no name dispatch.
    let text = std::fs::read_to_string(&args[1])?;
    let values: Vec<f64> = text.lines().flat_map(|l| l.split('#').next().unwrap_or("").split_whitespace())
        .map(str::parse).collect::<Result<_,_>>()?;
    if values.len() != 7 || !values.iter().all(|x| x.is_finite()) || values[4] <= 0.0 || values[6] <= 0.0 {
        return Err("expected rho cp k mu resistivity alpha Tref; positive resistivity and Tref".into());
    }
    let material = FluidTransportState::try_new(values[0], values[1], values[2], values[3])?;
    let conductor = Conductor::new(values[4]).with_temperature_coefficient(values[5], values[6]);
    let scale = plan_scaling_from_material(&material, 0.01, 0.01, 0.1)?;
    let mut slab = scale.thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    })?;
    let setup = ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,1], hi: [2,1,2], drive: ElectrodeDrive::Voltage(2.0) },
        Electrode { lo: [0,0,6], hi: [2,1,7], drive: ElectrodeDrive::Voltage(0.0) },
    ] };
    let limits = ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells: 16, temperature_range_k: [290.0, 320.0], balance_tolerance_j: 1e-9 },
        coupling_tolerance_j: 1e-10,
    };
    let gate = CancelGate::new();
    println!("# synthetic coefficients; imposed B; explicit coupling; work is lagged, not kinetic-energy change");
    println!("time_s,current_A,joule_W,lorentz_W,mean_vx_m_s,mean_K,coupling_defect_J");
    for index in 0..steps {
        let r = slab.step_magnetothermal(ElectricHeatInput { setup: &setup, contacts: &[], limits },
            &[[0.0,b,0.0];16], &[[0.0;3];16], &gate,
            |_,t,_| Ok(Some(conductor.resistivity(t))))?;
        let mut mean_u = 0.0; let mut mean_t = 0.0;
        for c in 2..14 {
            mean_u += scale.velocity_si(slab.lattice.grid.moments3(c).1[0])? / 12.0;
            mean_t += slab.temperature_si(slab.lattice.gpop[c].iter().sum())? / 12.0;
        }
        println!("{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.6e}",
            (index+1) as f64 * scale.timestep_s(), r.motion.electric.electrodes[0].current_a,
            r.motion.electric.total_joule_w, r.motion.mechanical_power_w, mean_u, mean_t, r.coupling_defect_j);
    }
    Ok(())
}
