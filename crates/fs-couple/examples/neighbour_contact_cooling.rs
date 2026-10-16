//! Two independently deformable roads sharing force and a declared heat patch.
//! Synthetic coefficients; point forces do not determine the supplied patch area.
use std::collections::BTreeMap;
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, ContactUpdate,
    DepositedCell, DepositionBudget, DepositionError, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::{bridge::{self, Bridge, neighbours, supports}, deposition};
use fs_solver::orthant::GramOrthantConfig;

type Result<T> = std::result::Result<T,Box<dyn std::error::Error>>;
const N: usize = 8;
const LENGTH: f64 = 0.02;
const AREA: f64 = 1e-6;
const PATCH_ID: u64 = 9000;
struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self,h: f64) -> std::result::Result<f64,DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static HEAT: Sensible = Sensible;

/// Explicit illustrative pressure-activated patch model, independent of mechanics.
/// Positive reaction above engage_force_n activates this GIVEN area and conductance.
#[derive(Clone,Copy)]
struct Patch { area_m2: f64, conductance_w_k: f64, engage_force_n: f64 }
#[derive(Clone)]
struct State<'m> { heat: DepositionSimulation<'static>, roads: BTreeMap<u64,Bridge<'m>> }
struct Step { thermal: EnergyAudit, mechanical: neighbours::Report }

fn beam_config() -> bridge::Config {
    bridge::Config { section: deposition::Config { force_tolerance_n:1e-7,
        moment_tolerance_nm:1e-9,energy_tolerance_j:1e-8,max_strain:0.03,max_rotation_rad:0.1 },
        maximum_deflection_m:0.001 }
}
fn mechanics_config() -> neighbours::Config {
    neighbours::Config { max_bodies:2,energy_tolerance_j:1e-8,
        contact:supports::Config { contact:GramOrthantConfig { max_variables:1,
            max_columns:2*N,max_setup_products:2*N,max_sweeps:4096,gradient_tolerance:1e-13 },
            gap_tolerance_m:1e-11,complementarity_tolerance_j:1e-13 } }
}
fn law(e: f64,branches: Vec<(f64,f64)>) -> Result<ThermoMaxwell> {
    Ok(ThermoMaxwell::new(e,branches,vec![
        ThermalKnot { temperature_k:300.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0 },
        ThermalKnot { temperature_k:500.0,modulus_scale:0.3,relaxation_scale:1.0,free_strain:0.002 },
    ],Limits { max_strain:0.03,max_strain_rate_s:10.0 })?)
}
fn thermal_id(body: u64,cell: u64) -> CellId { CellId(body*100+cell) }

impl<'m> State<'m> {
    fn new(upper: &'m ThermoMaxwell,lower: &'m ThermoMaxwell) -> Result<Self> {
        let dx=LENGTH/N as f64;
        let mut events=Vec::new(); let mut roads=BTreeMap::new();
        for (body,material,t) in [(1,upper,500.0),(2,lower,330.0)] {
            let parts: Vec<_>=(0..N).map(|i| (i as u64,dx)).collect();
            let road=Bridge::new(&parts,bridge::Support::SimplySupported,bridge::Budget {
                max_cells:N,max_layers_per_cell:1,max_branch_states:4*N })?
                .deposit_bonded_layer(7,deposition::LayerGeometry { area_m2:AREA,
                    centroid_z_m:0.0,thickness_m:0.001 },1200.0,material,
                    &parts.iter().map(|&(id,_)| (id,t)).collect(),&mut || false)?;
            roads.insert(body,road);
            for i in 0..N {
                let id=thermal_id(body,i as u64);
                let z=if body==1 { 0.001005 } else { 0.0 };
                events.push(TimedEvent { time_s:0.0,event:Event::Birth(DepositedCell {
                    id,segment_id:id.0,body_id:body,material_id:MaterialId(1),
                    mass_kg:1200.0*AREA*dx,specific_h_j_kg:2000.0*(t-300.0),
                    geometry:CellGeometry { start_m:[i as f64*dx,0.0,z],
                        end_m:[(i+1) as f64*dx,0.0,z],volume_m3:AREA*dx,
                        surface_area_m2:0.004*dx+2.0*AREA },
                    boundary:Boundary { ambient_k:300.0,htc_w_m2_k:20.0,reservoir:None,source_w:0.0 },
                }) });
            }
            for i in 1..N {
                events.push(TimedEvent { time_s:0.0,event:Event::Contact(ThermalLink {
                    id:body*1000+i as u64,a:thermal_id(body,(i-1) as u64),b:thermal_id(body,i as u64),
                    area_m2:AREA,conductance_w_k:0.2*AREA/dx,
                }) });
            }
        }
        let heat=DepositionSimulation::new(events,&[(MaterialId(1),&HEAT)],DepositionBudget {
            max_cells:2*N,max_links:2*N,max_events:4*N,max_substeps:10_000 })?
            .advance_to(0.0,StepConfig { max_dt_s:0.01,energy_tolerance_j:1e-8 },|| Ok(()))?.0;
        Ok(Self { heat,roads })
    }

    // Staggered heat -> solved relative contact -> next-interval thermal patch.
    // All three candidates stay local until the final cancellation checkpoint.
    fn advance(&self,target: f64,upper_extra: f64,patch: Patch,
        cancelled: &mut impl FnMut() -> bool) -> Result<(Self,Step)> {
        let checkpoint=|cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) }
        };
        checkpoint(cancelled)?;
        if !patch.area_m2.is_finite() || patch.area_m2<=0.0
            || !patch.conductance_w_k.is_finite() || patch.conductance_w_k<0.0
            || !patch.engage_force_n.is_finite() || patch.engage_force_n<0.0 {
            return Err(DepositionError::InvalidInput("invalid explicit heat patch").into());
        }
        let (heat,thermal)=self.heat.advance_to(target,
            StepConfig { max_dt_s:0.01,energy_tolerance_j:1e-8 },|| checkpoint(cancelled))?;
        let mut temperatures=BTreeMap::new();
        for (&body,road) in &self.roads {
            let mut field=BTreeMap::new();
            for cell in road.cells() {
                checkpoint(cancelled)?;
                field.insert((cell.id(),7),heat.temperature_k(thermal_id(body,cell.id()))?);
            }
            temperatures.insert(body,field);
        }
        let upper_load=[upper_extra;N]; let lower_load=[0.0;N];
        let inputs=[neighbours::Body { id:1,beam:&self.roads[&1],temperatures_k:&temperatures[&1],
            gravity_m_s2:9.81,additional_downward_n_m:&upper_load,config:beam_config() },
            neighbours::Body { id:2,beam:&self.roads[&2],temperatures_k:&temperatures[&2],
            gravity_m_s2:9.81,additional_downward_n_m:&lower_load,config:beam_config() }];
        let pair=neighbours::Contact { id:PATCH_ID,
            upper:neighbours::Station { body:1,reference_x_m:LENGTH/2.0 },
            lower:neighbours::Station { body:2,reference_x_m:LENGTH/2.0 },reference_gap_m:5e-6 };
        let solved=neighbours::advance_controlled(&inputs,&[pair],target,mechanics_config(),cancelled)?;
        let engaged=solved.report.reactions[0].force_on_upper_n>patch.engage_force_n;
        let replacement=engaged.then_some(ThermalLink { id:PATCH_ID,
            a:thermal_id(1,(N/2) as u64),b:thermal_id(2,(N/2) as u64),
            area_m2:patch.area_m2,conductance_w_k:patch.conductance_w_k });
        let heat=heat.replace_contacts(&[ContactUpdate { id:PATCH_ID,
            expected:heat.links().get(&PATCH_ID).copied(),replacement }],|| checkpoint(cancelled))?;
        checkpoint(cancelled)?;
        Ok((Self { heat,roads:solved.bodies },Step { thermal,mechanical:solved.report }))
    }
    fn mean_temperature(&self,body: u64) -> Result<f64> {
        let mut sum=0.0;
        for i in 0..N { sum+=self.heat.temperature_k(thermal_id(body,i as u64))?; }
        Ok(sum/N as f64)
    }
}

#[derive(Debug)]
struct ResultRow { lower_temperature_k: f64, maximum_contact_force_n: f64,
    open_after_lift: bool, energy_defect_j: f64 }
fn run(patch: Patch) -> Result<ResultRow> {
    let upper=law(1e7,vec![(9e7,0.2)])?; let lower=law(1e8,vec![])?;
    let mut state=State::new(&upper,&lower)?;
    let initial=state.heat.energy_j(); let mut external=0.0; let mut maximum=0.0_f64;
    for i in 1..=80 {
        let load=if i<=40 { 0.02 } else { -0.03 };
        let (next,step)=state.advance(i as f64*0.1,load,patch,&mut || false)?;
        external+=step.thermal.external_heat_j;
        maximum=maximum.max(step.mechanical.reactions[0].force_on_upper_n);
        assert!(step.mechanical.energy_residual_j.abs()<1e-8);
        state=next;
    }
    assert!((state.heat.mass_kg()-2.0*1200.0*AREA*LENGTH).abs()<1e-16);
    let defect=state.heat.energy_j()-initial-external;
    assert!(defect.abs()<1e-8);
    Ok(ResultRow { lower_temperature_k:state.mean_temperature(2)?,maximum_contact_force_n:maximum,
        open_after_lift:!state.heat.links().contains_key(&PATCH_ID),energy_defect_j:defect })
}
fn patch() -> Patch { Patch { area_m2:1e-6,conductance_w_k:0.003,engage_force_n:1e-12 } }
fn main() -> Result<()> {
    let conducting=run(patch())?;
    let insulating=run(Patch { conductance_w_k:0.0,..patch() })?;
    println!("Synthetic pressure-activated patch; both roads deform; endpoint-staggered coupling.");
    println!("conducting: {conducting:?}\ninsulating: {insulating:?}");
    assert!(conducting.maximum_contact_force_n>0.0 && conducting.open_after_lift);
    assert!(insulating.open_after_lift);
    assert!(conducting.lower_temperature_k>insulating.lower_temperature_k+0.1);
    assert!(conducting.energy_defect_j.abs()<1e-8 && insulating.energy_defect_j.abs()<1e-8);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn solved_neighbour_contact_exchanges_heat_and_lift_off_removes_it() {
        let conducting=run(patch()).unwrap();
        let insulating=run(Patch { conductance_w_k:0.0,..patch() }).unwrap();
        assert!(conducting.maximum_contact_force_n>0.0);
        assert!(conducting.open_after_lift && insulating.open_after_lift);
        assert!(conducting.lower_temperature_k>insulating.lower_temperature_k+0.1);
        assert!(conducting.energy_defect_j.abs()<1e-8 && insulating.energy_defect_j.abs()<1e-8);
    }
    #[test]
    fn late_cancellation_keeps_both_material_and_thermal_histories() {
        let upper=law(1e7,vec![(9e7,0.2)]).unwrap(); let lower=law(1e8,vec![]).unwrap();
        let state=State::new(&upper,&lower).unwrap();
        let before=state.clone(); let mut polls=0;
        state.advance(0.1,0.02,patch(),&mut || { polls+=1; false }).unwrap();
        let mut count=0;
        assert!(state.advance(0.1,0.02,patch(),&mut || { count+=1; count>=polls }).is_err());
        assert_eq!(state.heat.cells(),before.heat.cells());
        assert_eq!(state.heat.links(),before.heat.links());
        assert_eq!(state.roads,before.roads);
    }
}
