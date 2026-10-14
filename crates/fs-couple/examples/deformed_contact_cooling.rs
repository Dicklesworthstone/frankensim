//! Two freely contracting, flat roads with nonmatching meshes and thermal contact.
//! Reference-area ambient coefficients and projected-area contact coefficient
//! are explicitly synthetic. This demonstrates transfer, not pressure or adhesion.
use std::collections::{BTreeMap, BTreeSet};
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, ContactPatch,
    ContactPatchUpdate, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::{bridge::{self, Bridge, contact::pair}, deposition};

type Result<T> = std::result::Result<T,Box<dyn std::error::Error>>;
const LENGTH: f64 = 0.02;
const AREA: f64 = 1e-6;
struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self,h: f64) -> std::result::Result<f64,DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static HEAT: Sensible = Sensible;
fn thermal_id(body: u64,cell: u64) -> CellId { CellId(body*100+cell) }
fn thermal_config() -> StepConfig { StepConfig { max_dt_s: 0.01,energy_tolerance_j: 1e-8 } }
fn beam_config() -> bridge::Config {
    bridge::Config { section: deposition::Config { force_tolerance_n:1e-7,
        moment_tolerance_nm:1e-9,energy_tolerance_j:1e-8,max_strain:0.03,max_rotation_rad:0.1 },
        maximum_deflection_m:0.001 }
}
fn material() -> Result<ThermoMaxwell> {
    Ok(ThermoMaxwell::new(1e8,vec![],vec![
        ThermalKnot { temperature_k:300.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0 },
        ThermalKnot { temperature_k:550.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0025 },
    ],Limits { max_strain:0.03,max_strain_rate_s:10.0 })?)
}
#[derive(Clone)]
struct State<'a> {
    heat: DepositionSimulation<'static>,
    roads: BTreeMap<u64,Bridge<'a>>,
    patch_ids: BTreeSet<u64>,
    areas_m2: [f64;3], // lower reference, upper reference, common projected
}
impl<'a> State<'a> {
    fn new(law: &'a ThermoMaxwell, coefficient: f64) -> Result<Self> {
        let mut events=Vec::new(); let mut roads=BTreeMap::new();
        for (body,n,t,htc) in [(1,3,350.0,0.0),(2,5,500.0,100.0)] {
            let dx=LENGTH/n as f64;
            let parts: Vec<_>=(0..n).map(|i| (i as u64,dx)).collect();
            let initial=parts.iter().map(|&(id,_)| (id,t)).collect();
            let road=Bridge::new(&parts,bridge::Support::SimplySupported,bridge::Budget {
                max_cells:n,max_layers_per_cell:1,max_branch_states:0 })?
                .deposit_bonded_layer(7,deposition::LayerGeometry { area_m2:AREA,
                    centroid_z_m:0.0,thickness_m:0.001 },1200.0,law,&initial,&mut || false)?;
            roads.insert(body,road);
            for i in 0..n {
                let id=thermal_id(body,i as u64); let z=(body-1) as f64*0.001;
                events.push(TimedEvent { time_s:0.0,event:Event::Birth(DepositedCell {
                    id,segment_id:id.0,body_id:body,material_id:MaterialId(1),
                    mass_kg:1200.0*AREA*dx,specific_h_j_kg:2000.0*(t-300.0),
                    geometry:CellGeometry { start_m:[i as f64*dx,0.0,z],
                        end_m:[(i+1) as f64*dx,0.0,z],volume_m3:AREA*dx,
                        surface_area_m2:0.004*dx+2.0*AREA },
                    boundary:Boundary { ambient_k:300.0,htc_w_m2_k:htc,reservoir:None,source_w:0.0 },
                }) });
            }
            for i in 1..n {
                events.push(TimedEvent { time_s:0.0,event:Event::Contact(ThermalLink {
                    id:body*1000+i as u64,a:thermal_id(body,(i-1) as u64),b:thermal_id(body,i as u64),
                    area_m2:AREA,conductance_w_k:0.2*AREA/dx,
                }) });
            }
        }
        let heat=DepositionSimulation::new(events,&[(MaterialId(1),&HEAT)],DepositionBudget {
            max_cells:8,max_links:32,max_events:64,max_substeps:10000 })?
            .advance_to(0.0,thermal_config(),|| Ok(()))?.0;
        Self { heat,roads,patch_ids:BTreeSet::new(),areas_m2:[0.0;3] }.remap(coefficient)
    }

    // Same-time geometry transfer on a candidate, including patches that vanish
    // or split at the moving common partition. No full-field history is retained.
    fn remap(mut self, coefficient: f64) -> Result<Self> {
        if !coefficient.is_finite() || coefficient<0.0 {
            return Err(DepositionError::InvalidInput("invalid projected-area contact coefficient").into());
        }
        let footprint=self.roads[&1].neighbour_contact(&self.roads[&2],pair::Query {
            lower_origin_m:[0.0,0.0],upper_origin_m:[0.0,0.001],width_m:0.001,
            capture_gap_m:1e-10,penetration_tolerance_m:1e-11,
        },pair::Budget { max_cells:8,max_rows:16 },&mut || false)?;
        let mut desired=BTreeMap::new(); self.areas_m2=[0.0;3];
        for row in footprint.rows {
            if row.projected_area_m2==0.0 { continue; }
            // IDs are stable ordered cell pairs; meshes here have fewer than 100 cells.
            let id=1_000_000+100*row.lower_cell_id+row.upper_cell_id;
            desired.insert(id,ContactPatch { link:ThermalLink { id,
                a:thermal_id(1,row.lower_cell_id),b:thermal_id(2,row.upper_cell_id),
                area_m2:row.projected_area_m2,conductance_w_k:coefficient*row.projected_area_m2 },
                covered_area_m2:[row.lower_reference_area_m2,row.upper_reference_area_m2] });
            self.areas_m2[0]+=row.lower_reference_area_m2;
            self.areas_m2[1]+=row.upper_reference_area_m2;
            self.areas_m2[2]+=row.projected_area_m2;
        }
        let ids: BTreeSet<_>=self.patch_ids.iter().copied().chain(desired.keys().copied()).collect();
        let updates: Vec<_>=ids.into_iter().map(|id| ContactPatchUpdate {
            id,expected:self.heat.contact_patch(id),replacement:desired.get(&id).copied(),
        }).collect();
        self.heat=self.heat.replace_contact_patches(&updates,|| Ok(()))?;
        self.patch_ids=desired.into_keys().collect();
        Ok(self)
    }

    fn advance(&self,target: f64,coefficient: f64) -> Result<(Self,f64)> {
        let (heat,audit)=self.heat.advance_to(target,thermal_config(),|| Ok(()))?;
        let mut roads=BTreeMap::new();
        for (&body,road) in &self.roads {
            let mut temperatures=BTreeMap::new();
            for cell in road.cells() { temperatures.insert((cell.id(),7),heat.temperature_k(thermal_id(body,cell.id()))?); }
            let load=vec![0.0;road.cells().len()];
            let (updated,report)=road.advance_controlled(target,&temperatures,0.0,&load,beam_config(),&mut || false)?;
            assert!(report.energy_j.residual.abs()<1e-8);
            roads.insert(body,updated);
        }
        let next=Self { heat,roads,patch_ids:self.patch_ids.clone(),areas_m2:self.areas_m2 }.remap(coefficient)?;
        Ok((next,audit.external_heat_j))
    }
    fn mean(&self,body: u64) -> Result<f64> {
        let road=&self.roads[&body]; let mut temperature_mass=0.0; let mut mass=0.0;
        for cell in road.cells() {
            temperature_mass+=cell.mass_kg()*self.heat.temperature_k(thermal_id(body,cell.id()))?;
            mass+=cell.mass_kg();
        }
        Ok(temperature_mass/mass)
    }
}
#[derive(Debug)]
struct Outcome { lower_k:f64,upper_k:f64,areas_m2:[f64;3],energy_defect_j:f64 }
fn run(coefficient: f64) -> Result<Outcome> {
    let law=material()?; let mut state=State::new(&law,coefficient)?;
    let initial=state.heat.energy_j(); let mut external=0.0;
    for step in 1..=20 {
        let (next,q)=state.advance(step as f64*0.1,coefficient)?;
        external+=q; state=next;
    }
    let defect=state.heat.energy_j()-initial-external;
    assert!(defect.abs()<1e-8);
    assert!((state.areas_m2[0]-state.areas_m2[1]).abs()>1e-10);
    Ok(Outcome { lower_k:state.mean(1)?,upper_k:state.mean(2)?,
        areas_m2:state.areas_m2,energy_defect_j:defect })
}
fn main() -> Result<()> {
    let conducting=run(1000.0)?; let insulating=run(0.0)?;
    println!("Synthetic sliding flat faces: separate reference areas, projected-area conductance.");
    println!("conducting: {conducting:?}\ninsulating: {insulating:?}");
    assert!(conducting.lower_k>insulating.lower_k+1.0);
    assert!(conducting.upper_k<insulating.upper_k);
    assert!(conducting.areas_m2[0]<conducting.areas_m2[1]);
    assert!(conducting.energy_defect_j.abs()<1e-8);
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nonmatching_deformed_surfaces_drive_conservative_thermal_exchange() { main().unwrap(); }
    #[test]
    fn existing_frame_checkpoint_replays_geometry_heat_and_complete_side_coverage() {
        let law=material().unwrap(); let state=State::new(&law,1000.0).unwrap();
        let mid=state.advance(0.1,1000.0).unwrap().0;
        let a=mid.advance(0.2,1000.0).unwrap().0;
        let b=mid.clone().advance(0.2,1000.0).unwrap().0;
        assert_eq!(a.roads,b.roads); assert_eq!(a.heat.cells(),b.heat.cells());
        assert_eq!(a.patch_ids,b.patch_ids);
        for id in a.patch_ids { assert_eq!(a.heat.contact_patch(id),b.heat.contact_patch(id)); }
        assert!(mid.clone().remap(f64::NAN).is_err());
        assert_eq!(mid.heat.time_s(),0.1);
    }
}
