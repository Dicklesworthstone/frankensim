//! Example-only bounded table interpolation; no named-material calibration claim.
use std::{error::Error, io::Read};
use fs_exec::CancelGate;
use fs_lbm::{material::{MaterialThermalSlab3, plan_scaling_from_material}, thermal::sources::HeatSourceLimits};
use fs_material::{MaterialError, fluid_transport::FluidTransportState};

#[path = "material_transport_feedback/electrical.rs"]
mod electrical;
use fs_lbm::material::electrothermal::magnetic::transport::ElectricalTransport;

struct Table { rho:f64, cp:f64, knots:Vec<[f64;3]>, resistivity:Option<Vec<f64>> }
impl Table {
    fn parse(text:&str)->Result<Self,Box<dyn Error>> {
        let rows:Vec<Vec<f64>>=text.lines().filter_map(|l| {
            let line=l.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some(line)
        }).map(|l|l.split_whitespace().map(str::parse).collect()).collect::<Result<_,_>>()?;
        if rows.len()<3 || rows.len()>257 || rows[0].len()!=2 {
            return Err("expected rho cp followed by 2..256 rows of Kelvin conductivity viscosity".into());
        }
        let (rho,cp)=(rows[0][0],rows[0][1]);let mut knots=Vec::new();
        let columns=rows[1].len();
        if columns!=3 && columns!=4 {return Err("expected T k mu, optionally followed by resistivity".into());}
        let mut resistivity=(columns==4).then(Vec::new);
        for row in &rows[1..] {
            if row.len()!=columns || !row[0].is_finite() || row[0]<=0.0 {return Err("invalid Kelvin/k/mu row".into());}
            FluidTransportState::try_new(rho,cp,row[1],row[2])?;
            let knot=[row[0],row[1],row[2]];
            if knots.last().is_some_and(|p:&[f64;3]|p[0]>=knot[0]) {return Err("temperatures must strictly increase".into());}
            if let Some(values)=&mut resistivity {
                if !row[3].is_finite() || row[3]<=0.0 {return Err("resistivity must be finite and positive".into());}
                values.push(row[3]);
            }
            knots.push(knot);
        }
        Ok(Self {rho,cp,knots,resistivity})
    }
    fn at(&self,t:f64)->Result<FluidTransportState,MaterialError> {
        self.interpolate(t).map(|(transport,_)|transport)
    }
    fn electrical_at(&self,t:f64)->Result<ElectricalTransport<FluidTransportState>,MaterialError> {
        let (transport,resistivity)=self.interpolate(t)?;
        if resistivity.is_none() {return Err(MaterialError::State {what:"electrical mode requires a resistivity column; no fallback".into()});}
        ElectricalTransport::try_new(transport,resistivity)
    }
    fn interpolate(&self,t:f64)->Result<(FluidTransportState,Option<f64>),MaterialError> {
        if !t.is_finite() || t<self.knots[0][0] || t>self.knots[self.knots.len()-1][0] {
            return Err(MaterialError::State {what:format!("temperature {t} K outside example table support")});
        }
        let upper=self.knots.partition_point(|p|p[0]<t).max(1).min(self.knots.len()-1);
        let (a,b)=(self.knots[upper-1],self.knots[upper]);
        let f=(t-a[0])/(b[0]-a[0]);
        let transport=FluidTransportState::try_new(self.rho,self.cp,(1.0-f)*a[1]+f*b[1],(1.0-f)*a[2]+f*b[2])?;
        let resistivity=self.resistivity.as_ref().map(|r|(1.0-f)*r[upper-1]+f*r[upper]);
        Ok((transport,resistivity))
    }
}
fn midpoint_requested(args:&[String])->Result<bool,Box<dyn Error>> {
    if (2..=3).contains(&args.len()) { return Ok(false); }
    if args.get(3).is_some_and(|a|a=="--electrical") {
        if args.len()==6 { return Ok(false); }
        if args.len()==7 && args[6]=="--midpoint" { return Ok(true); }
    }
    Err("usage: material_transport_feedback TABLE.txt [steps: 2..2000] OR TABLE.txt STEPS --electrical VOLTAGE B_Y_TESLA [--midpoint]".into())
}

fn main()->Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().collect();
    let midpoint=midpoint_requested(&args)?;
    let steps:usize=args.get(2).map_or(Ok(200),|s|s.parse())?;
    if !(2..=2000).contains(&steps) {return Err("steps must be in 2..2000".into());}
    let mut text=String::new();
    std::fs::File::open(&args[1])?.take(65_537).read_to_string(&mut text)?;
    if text.len()>65_536 {return Err("table exceeds 64 KiB".into());}
    let table=Table::parse(&text)?;
    if args.len()>=6 {
        let (voltage,b)=(args[4].parse()?,args[5].parse()?);
        return if midpoint {electrical::run_midpoint(&table,steps,voltage,b)}
            else {electrical::run(&table,steps,voltage,b)};
    }
    let reference=table.at(305.0)?;
    let scale=plan_scaling_from_material(&reference,0.01,0.01,0.05)?;
    let mut slab=scale.thermal_slab3(MaterialThermalSlab3 {nx:2,ny:1,nz:8,
        bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:20})?;
    let limits=HeatSourceLimits {max_cells:20,temperature_range_k:[290.0,330.0],balance_tolerance_j:1e-8};
    let gate=CancelGate::new();let mut power=vec![0.0;20];power[2..18].fill(1e6);
    let force=[[2.0,0.0,0.0];20];
    let mut source_j=0.0;let mut storage_j=0.0;let mut walls_j=0.0;
    println!("# unqualified example table; fixed density/Cp; local k(T)/mu(T); explicit transport feedback");
    println!("time_s,mean_K,mean_vx_m_s,min_k_W_m_K,max_k_W_m_K,source_J,defect_J");
    for index in 0..steps {
        if index==steps/2 {power.fill(0.0);}
        let result=slab.step_with_material_transport(&power,&force,limits,&gate,|_,t,_|table.at(t))?;
        let mean_t=result.endpoint.iter().map(|s|s.temperature_k).sum::<f64>()/16.0;
        let mean_u=(2..18).map(|c|scale.velocity_si(slab.lattice.grid.moments3(c).1[0]))
            .collect::<Result<Vec<_>,_>>()?.iter().sum::<f64>()/16.0;
        let low=result.applied.iter().map(|s|s.transport.thermal_conductivity_w_m_k()).fold(f64::INFINITY,f64::min);
        let high=result.applied.iter().map(|s|s.transport.thermal_conductivity_w_m_k()).fold(f64::NEG_INFINITY,f64::max);
        let h=result.heat;source_j+=h.source_j;storage_j+=h.stored_change_j;walls_j+=h.bottom_wall_j+h.top_wall_j;
        println!("{:.6},{mean_t:.9},{mean_u:.9},{low:.6},{high:.6},{:.9},{:.6e}",
            (index+1) as f64*scale.timestep_s(),h.source_j,h.defect_j);
    }
    eprintln!("source_J={source_j:.9} storage_change_J={storage_j:.9} wall_input_J={walls_j:.9} defect_J={:.6e}",storage_j-source_j-walls_j);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_interpolation_is_bounded_and_changes_data_not_physics() {
        let t=Table::parse("1000 1000\n300 1000 1\n320 3000 0.25\n").unwrap();
        let s=t.at(310.0).unwrap();assert_eq!(s.thermal_conductivity_w_m_k(),2000.0);assert_eq!(s.dynamic_viscosity_pa_s(),0.625);
        assert!(t.at(299.0).is_err());assert!(t.at(321.0).is_err());
        assert!(Table::parse("1000 1000\n300 1000 1\n300 2000 1").is_err());
        assert!(Table::parse("1000 1000\n300 -1 1\n320 2000 1").is_err());
    }
    #[test]
    fn joint_table_uses_one_interpolation_and_refuses_missing_electrical_data() {
        let t=Table::parse("1000 1000\n300 100 1 2\n320 300 0.5 4\n").unwrap();
        let s=t.electrical_at(310.0).unwrap();
        assert_eq!(s.resistivity_ohm_m(),Some(3.0));
        assert_eq!(s.material().thermal_conductivity_w_m_k(),200.0);
        assert_eq!(s.material().dynamic_viscosity_pa_s(),0.75);
        assert!(t.electrical_at(321.0).is_err());
        let old=Table::parse("1000 1000\n300 100 1\n320 300 0.5\n").unwrap();
        assert!(old.at(310.0).is_ok()); assert!(old.electrical_at(310.0).is_err());
        for text in ["1000 1000\n300 100 1 2\n320 300 0.5\n", "1000 1000\n300 100 1 -2\n320 300 0.5 4\n"] {
            assert!(Table::parse(text).is_err());
        }
    }

    #[test]
    fn midpoint_flag_is_opt_in_and_preserves_original_command_forms() {
        let args=|tail:&str|std::iter::once("material_transport_feedback")
            .chain(tail.split_whitespace()).map(str::to_owned).collect::<Vec<_>>();
        for tail in ["table.txt", "table.txt 20", "table.txt 20 --electrical 10 0.5"] {
            assert!(!midpoint_requested(&args(tail)).unwrap());
        }
        assert!(midpoint_requested(&args("table.txt 20 --electrical 10 0.5 --midpoint")).unwrap());
        for tail in ["", "table.txt 20 --midpoint", "table.txt 20 --electrical 10",
            "table.txt 20 --electrical 10 0.5 --other", "table.txt 20 --electrical 10 0.5 --midpoint extra"] {
            assert!(midpoint_requested(&args(tail)).is_err());
        }
    }

}
