//! Data-only material substitution in a current-driven thermal/flow evolution.
//! Inputs are explicitly synthetic, not calibration data for a named liquid.
use std::error::Error;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{Conductor, ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::midpoint::MidpointOptions;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 2 || args.len() > 3 {
        return Err("usage: material_joule_flow MATERIAL.txt [steps: 2..500]".into());
    }
    let steps: usize = args.get(2).map_or(Ok(100), |s| s.parse())?;
    if !(2..=500).contains(&steps) { return Err("steps must be in 2..500".into()); }
    let text = std::fs::read_to_string(&args[1])?;
    let values: Vec<f64> = text.lines().flat_map(|line| line.split('#').next().unwrap_or("").split_whitespace())
        .map(str::parse).collect::<Result<_,_>>()?;
    if values.len() != 7 || !values.iter().all(|x| x.is_finite()) {
        return Err("expected seven finite SI values: rho cp k mu resistivity alpha Tref".into());
    }
    let material = FluidTransportState::try_new(values[0], values[1], values[2], values[3])?;
    if values[4] <= 0.0 || values[6] <= 0.0 { return Err("resistivity and reference Kelvin must be positive".into()); }
    let conductor = Conductor::new(values[4]).with_temperature_coefficient(values[5], values[6]);
    let dx = 0.01; let dt = 0.01;
    let scale = plan_scaling_from_material(&material, dx, dt, 0.05)?;
    let mut slab = scale.thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 9.81, expansion_per_k: 1e-4, max_cells: 16,
    })?;
    let mut setup = ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,1], hi: [2,1,2], drive: ElectrodeDrive::Current(0.5) },
        Electrode { lo: [0,0,6], hi: [2,1,7], drive: ElectrodeDrive::Voltage(0.0) },
    ] };
    let limits = ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells: 16, temperature_range_k: [290.0, 400.0], balance_tolerance_j: 1e-8 },
        coupling_tolerance_j: 1e-8,
    };
    let midpoint = MidpointOptions { max_iterations: 32, temperature_tolerance_k: 1e-8,
        source_tolerance_j: 1e-8, relaxation: 1.0 };
    let gate = CancelGate::new();
    let forces = [[0.0;3];16];
    let mut electrode_j = 0.0; let mut stored_j = 0.0; let mut wall_j = 0.0;
    println!("# synthetic material; ideal electrical constraints, Boussinesq sensible heat; no calibration claim");
    println!("time_s,current_A,power_W,mean_fluid_K,max_fluid_K,nonlinear_trials,combined_defect_J");
    for index in 0..steps {
        // Change only the drive. No population, temperature or clock reset.
        if index == steps/2 { setup.electrodes[0].drive = ElectrodeDrive::Current(0.0); }
        let result = slab.step_electrothermal_midpoint(ElectricHeatInput { setup: &setup, contacts: &[], limits },
            &forces, midpoint, &gate, |_,t,_| Ok(Some(conductor.resistivity(t))))?;
        let temperatures: Vec<_> = slab.lattice.gpop[2..14].iter()
            .map(|p| slab.temperature_si(p.iter().sum())).collect::<Result<_,_>>()?;
        let mean = temperatures.iter().sum::<f64>()/temperatures.len() as f64;
        let maximum = temperatures.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let record = &result.step;
        electrode_j += record.electrode_energy_j;
        stored_j += record.heat.stored_change_j;
        wall_j += record.heat.bottom_wall_j+record.heat.top_wall_j;
        println!("{:.6},{:.9},{:.9},{:.9},{:.9},{},{:.6e}", (index+1) as f64*dt,
            record.electric.electrodes[0].current_a, record.electric.delivered_w,
            mean, maximum, result.iterations, record.combined_defect_j);
    }
    eprintln!("electrode_J={electrode_j:.9} stored_change_J={stored_j:.9} wall_input_J={wall_j:.9} defect_J={:.6e}", stored_j-wall_j-electrode_j);
    Ok(())
}
