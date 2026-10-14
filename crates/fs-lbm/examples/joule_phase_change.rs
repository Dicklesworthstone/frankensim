//! Synthetic caloric/electrical data drive the existing current and heat owners.
use std::{error::Error,io::Read};
use fs_blake3::DomainHasher;
use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError,ElectricSetup,Electrode,ElectrodeDrive,VoxelDomain};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError,JouleEnthalpyInput,advance_joule_enthalpy};
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransport,advance_joule_enthalpy_with_material_transport};
use fs_material::phase::{EnthalpyPhaseKnot,EquilibriumEnthalpyPhaseCurve,EquilibriumPhaseState};

struct Table {
    curve:EquilibriumEnthalpyPhaseCurve,
    density:f64,
    conductivity:f64,
    electrical:Vec<[f64;2]>,
    local_conductivity:Option<Vec<f64>>,
}
impl Table {
    fn parse(text:&str)->Result<Self,Box<dyn Error>> {
        if text.len()>65_536 {return Err("example table exceeds 64 KiB".into());}
        let rows:Vec<Vec<f64>>=text.lines().filter_map(|line| {
            let line=line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some(line)
        }).map(|line|line.split_whitespace().map(str::parse).collect()).collect::<Result<_,_>>()?;
        if !(3..=257).contains(&rows.len()) || rows[0].len()!=2 {
            return Err("expected density conductivity, then 2..256 h T liquid_fraction resistivity rows".into());
        }
        let (density,conductivity)=(rows[0][0],rows[0][1]);
        if ![density,conductivity].iter().all(|v|v.is_finite()&&*v>0.0) {
            return Err("density and conductivity must be finite and positive".into());
        }
        let mut knots=Vec::new();let mut electrical=Vec::new();
        let columns=rows[1].len();
        if columns!=4 && columns!=5 {return Err("expected h T fraction resistivity [conductivity]".into());}
        let mut local_conductivity=(columns==5).then(Vec::new);
        for row in &rows[1..] {
            if row.len()!=columns || !row.iter().all(|v|v.is_finite()) || row[3]<=0.0 {
                return Err("invalid h/T/fraction/resistivity row".into());
            }
            knots.push(EnthalpyPhaseKnot {specific_enthalpy_j_kg:row[0],temperature_k:row[1],
                liquid_mass_fraction:row[2],bulk_density_kg_m3:density});
            electrical.push([row[0],row[3]]);
            if let Some(values)=&mut local_conductivity {
                if row[4]<0.0 {return Err("conductivity must be nonnegative".into());}
                values.push(row[4]);
            }
        }
        if knots[0].temperature_k!=300.0 {
            return Err("this example starts both specimen and ambient at an explicit 300 K; first knot must match".into());
        }
        // Identity of SYNTHETIC input bytes, not a claim of a sourced material card.
        let mut identity=DomainHasher::new("org.frankensim.example.joule-phase.v1");
        identity.update(text.as_bytes());
        let curve=EquilibriumEnthalpyPhaseCurve::try_new(identity.finalize(),knots)?;
        Ok(Self {curve,density,conductivity,electrical,local_conductivity})
    }
    fn resistivity(&self,h:f64)->Result<f64,ChtError> {
        self.interpolate(h).map(|(_,r,_)|r)
    }
    fn transport(&self,h:f64)->Result<EnthalpyTransport<EquilibriumPhaseState>,JouleEnthalpyError> {
        let (phase,r,k)=self.interpolate(h)?;
        EnthalpyTransport::try_new(phase,Some(r),k)
    }
    fn interpolate(&self,h:f64)->Result<(EquilibriumPhaseState,f64,f64),ChtError> {
        let phase=self.curve.state_at_specific_enthalpy(h).map_err(|e|ChtError::InvalidInput {
            field:"example.caloric",reason:e.to_string(),
        })?;
        let upper=self.electrical.partition_point(|p|p[0]<h).max(1).min(self.electrical.len()-1);
        let (a,b)=(self.electrical[upper-1],self.electrical[upper]);
        let f=(h-a[0])/(b[0]-a[0]);
        let k=self.local_conductivity.as_ref().map_or(self.conductivity,|values| {
            (1.0-f)*values[upper-1]+f*values[upper]
        });
        Ok((phase,(1.0-f)*a[1]+f*b[1],k))
    }
}

#[derive(Clone,Copy,Debug)]
struct Row {
    time:f64,temperature:f64,liquid_fraction:f64,enthalpy:f64,current:f64,
    stored:f64,electrode:f64,thermal_external:f64,residual:f64,solves:usize,conductivity:f64,
}

fn simulate(table:&Table,voltage:f64,duration:f64,max_dt:f64,mut emit:impl FnMut(Row))
    ->Result<Row,Box<dyn Error>> {
    if !voltage.is_finite() || !duration.is_finite() || !(0.0<duration&&duration<=600.0)
        || !max_dt.is_finite() || max_dt<=0.0 || max_dt>1.0 {
        return Err("finite voltage, duration in (0,600] s and max_dt in (0,1] s required".into());
    }
    // Three fixed 1 cm control volumes, with finite-mass thermal electrode cells.
    // This is a stationary contained caloric model, not a free-surface melt flow.
    let domain=VoxelDomain::new(3,1,1,0.01)?;let ids=[CellId(1),CellId(2),CellId(3)];
    let h0=table.electrical[0][0];let mut events=Vec::new();
    for c in 0..3 {
        let [x,y,z]=domain.coords(c);let centre=domain.center(x,y,z);let width=domain.widths(c);
        let mut start=centre;let mut end=centre;start[0]-=0.5*width[0];end[0]+=0.5*width[0];
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:ids[c],segment_id:c as u64,body_id:1,material_id:MaterialId(1),
            mass_kg:table.density*domain.volume(c),specific_h_j_kg:h0,
            geometry:CellGeometry {start_m:start,end_m:end,volume_m3:domain.volume(c),
                surface_area_m2:2.0*(0..3).map(|a|domain.face_area(c,a)).sum::<f64>()},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.1,reservoir:None,source_w:0.0},
        })});
    }
    for c in 0..2 {
        let area=domain.face_area(c,0);
        let length=0.5*domain.widths(c)[0]+0.5*domain.widths(c+1)[0];
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
            id:c as u64,a:ids[c],b:ids[c+1],area_m2:area,conductance_w_k:table.conductivity*area/length,
        })});
    }
    let config=StepConfig {max_dt_s:max_dt,energy_tolerance_j:1e-8};
    let mut heat=DepositionSimulation::new(events,&[(MaterialId(1),&table.curve as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:2,max_events:5,max_substeps:20_000})?
        .advance_to(0.0,config,||Ok(()))?.0;
    let initial_energy=heat.energy_j();let gate=CancelGate::new();
    let mut electric=ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive:ElectrodeDrive::Voltage(voltage)},
        Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ]};
    let (mut electrode_j,mut external_j,mut solves)=(0.0,0.0,0usize);
    let mut last=None;
    for index in 1..=60 {
        // Switch at an exact coupling endpoint, never inside a frozen-source step.
        if index==31 {electric.electrodes[0].drive=ElectrodeDrive::Voltage(0.0);}
        let target=if index==30 {0.5*duration} else {duration*(index as f64/60.0)};
        let input=JouleEnthalpyInput {domain:&domain,setup:&electric,contacts:&[],cell_ids:&ids,
            max_cells:3,max_current_solves:20_001,power_tolerance_w:1e-9,energy_tolerance_j:1e-8};
        let (next,report,conductivity)=if table.local_conductivity.is_some() {
            let (next,r)=advance_joule_enthalpy_with_material_transport(&heat,target,config,input,&[0,1],&gate,
                |_,cell,_|table.transport(cell.specific_h_j_kg),|patch,sides| {
                    // Explicit Cartesian half-cell series resistance, with no
                    // extra thermal interface resistance in this example.
                    let [a,b]=sides.map(|s|s.response.conductivity_w_m_k());
                    let g=if a==0.0 || b==0.0 {0.0} else {
                        let length=sides.map(|s|s.cell.geometry.end_m[0]-s.cell.geometry.start_m[0]);
                        patch.link.area_m2/(0.5*length[0]/a+0.5*length[1]/b)
                    };
                    Ok((g,()))
                })?;
            let k=r.endpoint.cells[1].response.conductivity_w_m_k();
            (next,r.joule,k)
        } else {
            let (next,r)=advance_joule_enthalpy(&heat,target,config,input,&gate,
                |_,cell,_|table.resistivity(cell.specific_h_j_kg).map(Some))?;
            (next,r,table.conductivity)
        };
        let state=table.curve.state_at_specific_enthalpy(next.cells()[&ids[1]].specific_h_j_kg)?;
        electrode_j+=report.electrode_energy_j;external_j+=report.thermal.external_heat_j;solves+=report.current_solves;
        let stored=next.energy_j()-initial_energy;
        let row=Row {time:next.time_s(),temperature:state.temperature_k(),liquid_fraction:state.liquid_mass_fraction(),
            enthalpy:next.cells()[&ids[1]].specific_h_j_kg,current:report.endpoint.electric.electrodes[0].current_a,
            stored,electrode:electrode_j,thermal_external:external_j,residual:stored-external_j,solves,conductivity};
        heat=next; // all actual electrical, enthalpy, endpoint and energy checks passed
        emit(row);last=Some(row);
    }
    Ok(last.expect("positive fixed observation count"))
}

fn main()->Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().collect();
    if !(2..=5).contains(&args.len()) {return Err("usage: joule_phase_change MATERIAL.txt [voltage=10] [duration_s=60] [max_dt_s=0.05]".into());}
    let mut text=String::new();std::fs::File::open(&args[1])?.take(65_537).read_to_string(&mut text)?;
    let table=Table::parse(&text)?;
    let value=|i:usize,default:f64|->Result<f64,Box<dyn Error>> {Ok(match args.get(i) {Some(s)=>s.parse()?,None=>default})};
    println!("# synthetic fixed-volume phase chart; explicit Ohmic source; no melt flow or density/pressure work");
    println!("time_s,centre_K,liquid_mass_fraction,centre_h_J_kg,endpoint_current_A,stored_change_J,electrode_work_J,total_external_heat_J,thermal_defect_J,current_solves{}",
        if table.local_conductivity.is_some() {",endpoint_centre_k_W_m_K"} else {""});
    simulate(&table,value(2,10.0)?,value(3,60.0)?,value(4,0.05)?,|r| {
        let suffix=if table.local_conductivity.is_some() {format!(",{:.9}",r.conductivity)} else {String::new()};
        println!("{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.6e},{}{}",
            r.time,r.temperature,r.liquid_fraction,r.enthalpy,r.current,r.stored,r.electrode,r.thermal_external,r.residual,r.solves,suffix);
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn same_temperature_does_not_erase_phase_dependent_electrical_response() {
        let table=Table::parse(include_str!("joule_phase_change/latent_20k.txt")).unwrap();
        assert_eq!(table.curve.state_at_specific_enthalpy(10000.0).unwrap().temperature_k(),310.0);
        assert_eq!(table.curve.state_at_specific_enthalpy(30000.0).unwrap().temperature_k(),310.0);
        assert_eq!(table.resistivity(10000.0).unwrap(),1.0);
        assert_eq!(table.resistivity(30000.0).unwrap(),2.0);
        assert!(table.resistivity(-1.0).is_err());
        assert!(Table::parse("1000 0.01\n0 300 0 1\n1 310 0 2").is_err());
    }
    #[test]
    fn data_only_latent_heat_changes_current_and_phase_with_power_off_memory() {
        let a=Table::parse(include_str!("joule_phase_change/latent_20k.txt")).unwrap();
        let b=Table::parse(include_str!("joule_phase_change/latent_40k.txt")).unwrap();
        let mut rows_a=Vec::new();let mut rows_b=Vec::new();
        let x=simulate(&a,10.0,60.0,0.05,|r|rows_a.push(r)).unwrap();
        let y=simulate(&b,10.0,60.0,0.05,|r|rows_b.push(r)).unwrap();
        assert!((x.temperature-310.0).abs()<1e-8 && (y.temperature-310.0).abs()<1e-8);
        assert!(x.liquid_fraction>y.liquid_fraction+0.1);
        assert!(rows_a[29].current<rows_b[29].current);
        assert!(x.stored>0.0 && y.stored>0.0 && x.current==0.0 && y.current==0.0);
        assert!(rows_a[30..].iter().all(|r|r.electrode==rows_a[29].electrode));
        assert!(x.residual.abs()<1e-6 && y.residual.abs()<1e-6);
    }
    #[test]
    fn joint_enthalpy_table_retains_plateau_transport_and_refuses_mixed_rows() {
        let t=Table::parse(include_str!("joule_phase_change/transport_phase_dependent.txt")).unwrap();
        let a=t.transport(10000.0).unwrap();let b=t.transport(20000.0).unwrap();
        assert_eq!(a.material().temperature_k(),b.material().temperature_k());
        assert_eq!(a.conductivity_w_m_k(),0.01);
        assert!((b.conductivity_w_m_k()-0.105).abs()<1e-12);
        assert_eq!(b.resistivity_ohm_m(),Some(1.5));
        assert!(t.transport(-1.0).is_err());
        for text in ["1000 0.01\n0 300 0 1 0.01\n30000 310 1 2",
            "1000 0.01\n0 300 0 1 -0.1\n30000 310 1 2 0.2"] {
            assert!(Table::parse(text).is_err());
        }
    }
    #[test]
    fn transport_only_data_changes_actual_phase_current_and_post_pulse_history() {
        let old=Table::parse(include_str!("joule_phase_change/latent_20k.txt")).unwrap();
        let a=Table::parse(include_str!("joule_phase_change/transport_constant.txt")).unwrap();
        let b=Table::parse(include_str!("joule_phase_change/transport_phase_dependent.txt")).unwrap();
        let legacy=simulate(&old,10.0,60.0,0.05,|_|{}).unwrap();
        let mut rows_a=Vec::new();let mut rows_b=Vec::new();
        let x=simulate(&a,10.0,60.0,0.05,|r|rows_a.push(r)).unwrap();
        let y=simulate(&b,10.0,60.0,0.05,|r|rows_b.push(r)).unwrap();
        assert!((legacy.enthalpy-x.enthalpy).abs()<1e-7);
        assert!((legacy.electrode-x.electrode).abs()<1e-9);
        assert!((x.liquid_fraction-y.liquid_fraction).abs()>1e-4);
        assert!((rows_a[29].current-rows_b[29].current).abs()>1e-7);
        assert!(y.conductivity>x.conductivity);
        assert!(y.enthalpy>0.0 && y.current==0.0);
        assert!(rows_b[30..].iter().all(|r|r.electrode==rows_b[29].electrode));
        assert!(x.residual.abs()<1e-6 && y.residual.abs()<1e-6);
    }

}
