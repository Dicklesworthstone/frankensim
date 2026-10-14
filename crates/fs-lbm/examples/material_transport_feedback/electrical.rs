//! Generated electrical heat/force through the existing table-driven example.
use super::*;
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};

#[derive(Clone, Copy, Debug, Default)]
struct Record {
    time_s:f64, mean_k:f64, mean_vx:f64, current_a:f64, joule_w:f64,
    mechanical_w:f64, min_k:f64, max_k:f64, min_mu:f64, max_mu:f64,
    defect_j:f64,
}
#[derive(Default)]
struct Summary {
    last:Record, electrode_j:f64, joule_j:f64, mechanical_j:f64,
    stored_j:f64, wall_j:f64,
}
impl Summary {
    fn defect(&self)->f64 {self.stored_j-self.wall_j+self.mechanical_j-self.electrode_j}
}

fn simulate(table:&Table,steps:usize,voltage:f64,b:f64,mut emit:impl FnMut(Record))
    ->Result<Summary,Box<dyn Error>>
{
    if !(2..=2000).contains(&steps) || !voltage.is_finite() || !b.is_finite() {
        return Err("electrical mode requires 2..2000 steps and finite voltage/B".into());
    }
    // Missing resistance refuses rather than turning into an insulating fluid.
    let reference=table.electrical_at(305.0)?;
    let scale=plan_scaling_from_material(reference.material(),0.01,0.01,0.05)?;
    let mut slab=scale.thermal_slab3(MaterialThermalSlab3 {
        nx:2,ny:1,nz:8,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:20,
    })?;
    let mut setup=ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,1],hi:[2,1,2],drive:ElectrodeDrive::Voltage(voltage)},
        Electrode {lo:[0,0,8],hi:[2,1,9],drive:ElectrodeDrive::Voltage(0.0)},
    ]};
    let limits=ElectricHeatLimits {
        thermal:HeatSourceLimits {max_cells:20,temperature_range_k:[290.0,330.0],balance_tolerance_j:1e-8},
        coupling_tolerance_j:1e-10,
    };
    let field=[[0.0,b,0.0];20]; let force=[[2.0,0.0,0.0];20];
    let gate=CancelGate::new(); let mut totals=Summary::default();
    for index in 0..steps {
        // Zero voltage shorts the ideal electrodes; it does NOT open the circuit.
        // Motional currents and braking can still generate heat after switch-off.
        if index==steps/2 {setup.electrodes[0].drive=ElectrodeDrive::Voltage(0.0);}
        let result=slab.step_magnetothermal_with_material_transport(
            ElectricHeatInput {setup:&setup,contacts:&[],limits},&field,&force,&gate,
            |_,t,_|table.electrical_at(t),
        )?;
        let mean_k=result.endpoint.iter().map(|s|s.temperature_k).sum::<f64>()/16.0;
        let mean_vx=(2..18).map(|c|scale.velocity_si(slab.lattice.grid.moments3(c).1[0]))
            .collect::<Result<Vec<_>,_>>()?.iter().sum::<f64>()/16.0;
        let mut row=Record {time_s:(index+1) as f64*scale.timestep_s(),mean_k,mean_vx,
            current_a:result.step.motion.electric.electrodes[0].current_a,
            joule_w:result.step.motion.electric.total_joule_w,
            mechanical_w:result.step.motion.mechanical_power_w,
            min_k:f64::INFINITY,max_k:f64::NEG_INFINITY,min_mu:f64::INFINITY,max_mu:f64::NEG_INFINITY,
            defect_j:result.step.coupling_defect_j};
        for sample in &result.applied {
            let k=sample.transport.thermal_conductivity_w_m_k();
            let mu=sample.transport.dynamic_viscosity_pa_s();
            row.min_k=row.min_k.min(k); row.max_k=row.max_k.max(k);
            row.min_mu=row.min_mu.min(mu); row.max_mu=row.max_mu.max(mu);
        }
        let step=&result.step;
        totals.electrode_j+=step.electrode_energy_j; totals.joule_j+=step.joule_energy_j;
        totals.mechanical_j+=step.lagged_mechanical_work_j;
        totals.stored_j+=step.heat.stored_change_j;
        totals.wall_j+=step.heat.bottom_wall_j+step.heat.top_wall_j;
        totals.last=row;
        emit(row);
    }
    Ok(totals)
}

pub(super) fn run(table:&Table,steps:usize,voltage:f64,b:f64)->Result<(),Box<dyn Error>> {
    println!("# synthetic table; joint k(T)/mu(T)/resistivity(T); imposed B; explicit sources; work is not kinetic-energy change");
    println!("time_s,mean_K,mean_vx_m_s,current_A,joule_W,lorentz_W,min_k_W_m_K,max_k_W_m_K,min_mu_Pa_s,max_mu_Pa_s,coupling_defect_J");
    let totals=simulate(table,steps,voltage,b,|r| {
        println!("{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.6},{:.6},{:.6},{:.6},{:.6e}",
            r.time_s,r.mean_k,r.mean_vx,r.current_a,r.joule_w,r.mechanical_w,r.min_k,r.max_k,r.min_mu,r.max_mu,r.defect_j);
    })?;
    eprintln!("electrode_J={:.9} joule_J={:.9} lagged_mechanical_J={:.9} stored_J={:.9} wall_J={:.9} coupling_defect_J={:.6e}",
        totals.electrode_j,totals.joule_j,totals.mechanical_j,totals.stored_j,totals.wall_j,totals.defect());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn thermal_material_data_changes_subsequent_current_as_well_as_flow_and_heat() {
        let constant=Table::parse(include_str!("electrical_constant.txt")).unwrap();
        let variable=Table::parse(include_str!("electrical_temperature_dependent.txt")).unwrap();
        // Same rho/Cp and resistivity(T); only k and mu are changed in the data.
        for t in [300.0,305.0,310.0] {
            assert_eq!(constant.electrical_at(t).unwrap().resistivity_ohm_m(),variable.electrical_at(t).unwrap().resistivity_ohm_m());
        }
        let mut a=Vec::new(); let mut b=Vec::new();
        let x=simulate(&constant,40,10.0,0.5,|r|a.push(r)).unwrap();
        let y=simulate(&variable,40,10.0,0.5,|r|b.push(r)).unwrap();
        assert!((a[0].current_a-b[0].current_a).abs()<1e-12);
        assert!((x.last.current_a-y.last.current_a).abs()>1e-12);
        assert!((x.last.mean_vx-y.last.mean_vx).abs()>1e-8);
        assert!((x.last.mean_k-y.last.mean_k).abs()>1e-6);
        assert!(x.defect().abs()<1e-7 && y.defect().abs()<1e-7);
        assert!(a[20..].iter().any(|r|r.joule_w>0.0 && r.mechanical_w<0.0));
    }
    #[test]
    fn zero_magnetic_field_keeps_ohmic_heat_without_magnetic_work() {
        let table=Table::parse(include_str!("electrical_constant.txt")).unwrap();
        let summary=simulate(&table,20,2.0,0.0,|r|assert_eq!(r.mechanical_w,0.0)).unwrap();
        assert_eq!(summary.last.current_a,0.0);
        assert!(summary.joule_j>0.0); assert!(summary.last.mean_k>305.0);
        assert!(summary.defect().abs()<1e-7);
    }
}
