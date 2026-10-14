//! Explicit synthetic material and schedule; no calibrated printer parameters.
use std::collections::BTreeMap;
use fs_conduction::transient::deposition::*;
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::bridge::{self, Bridge, Support};
use fs_solid::bridge::drawing::Drawing;
use fs_solid::deposition::{Config as SectionConfig, LayerGeometry};
use super::driver::{Binding, Config, Error, GrowingRoad};

pub const AREA: f64 = 2.5e-7;
pub const DX: f64 = 0.0025;
pub const BODY: u64 = 77;
struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=500_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static HEAT: Sensible = Sensible;

pub fn material() -> ThermoMaxwell {
    ThermoMaxwell::new(3e7,vec![(7e7,0.5)],vec![
        ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 550.0, modulus_scale: 0.2, relaxation_scale: 0.1, free_strain: 0.0125 },
    ],Limits { max_strain: 0.03,max_strain_rate_s: 10.0 }).unwrap()
}
pub fn config() -> Config {
    Config { mechanics: bridge::Config {
        section: SectionConfig { force_tolerance_n: 1e-7,moment_tolerance_nm: 1e-9,
            energy_tolerance_j: 1e-9,max_strain: 0.03,max_rotation_rad: 0.1 },
        maximum_deflection_m: 0.002 }, gravity_m_s2: 9.81, drawing: Drawing::Slack,
        coordinate_tolerance_m: 1e-12, volume_tolerance_m3: 1e-18,
        mass_tolerance_kg: 1e-14, temperature_tolerance_k: 1e-10,
        max_bound_cells: 32,max_material_point_updates: 100_000 }
}
pub fn stepping(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt,energy_tolerance_j: 1e-9 } }
pub fn source_cell(id: u64, temperature_k: f64) -> DepositedCell {
    DepositedCell { id: CellId(id),segment_id: 1000+id,body_id: BODY,material_id: MaterialId(3),
        mass_kg: 1200.0*AREA*DX,specific_h_j_kg: (temperature_k-300.0)*2000.0,
        geometry: CellGeometry { start_m: [id as f64*DX,0.0,0.0],
            end_m: [(id+1) as f64*DX,0.0,0.0],volume_m3: AREA*DX,
            surface_area_m2: 0.002*DX+2.0*AREA },
        boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 25.0,reservoir: None,source_w: 0.0 } }
}
pub fn schedule() -> Vec<TimedEvent> {
    let mut out: Vec<_> = (0..4).map(|id| TimedEvent { time_s: 0.0,event: Event::Birth(source_cell(id,450.0)) }).collect();
    let link = |id| Event::Contact(ThermalLink { id,a: CellId(id-1),b: CellId(id),
        area_m2: AREA,conductance_w_k: 0.15*AREA/DX });
    out.extend((1..4).map(|id| TimedEvent { time_s: 0.0,event: link(id) }));
    for (id,time_s) in [(4,0.5),(5,0.5),(6,1.0),(7,1.5)] {
        out.push(TimedEvent { time_s,event: Event::Birth(source_cell(id,500.0)) });
        out.push(TimedEvent { time_s,event: link(id) });
    }
    out
}
pub fn initial<'a>(law: &'a ThermoMaxwell, events: Vec<TimedEvent>)
    -> Result<(DepositionSimulation<'static>,Bridge<'a>),Error> {
    let thermal = DepositionSimulation::new(events,&[(MaterialId(3),&HEAT)],DepositionBudget {
        max_cells: 32,max_links: 64,max_events: 100,max_substeps: 10_000 })?
        .advance_to(0.0,stepping(0.03125),|| Ok(()))?.0;
    let cells = (0..4).map(|id| (id,DX)).collect::<Vec<_>>();
    let t: BTreeMap<_,_> = cells.iter().map(|&(id,_)| Ok((id,thermal.temperature_k(CellId(id))?)))
        .collect::<Result<_,DepositionError>>()?;
    let beam = Bridge::new(&cells,Support::CantileverLeft,bridge::Budget {
        max_cells: 32,max_layers_per_cell: 1,max_branch_states: 64 })?
        .deposit_bonded_layer(9,LayerGeometry { area_m2: AREA,centroid_z_m: 0.0,thickness_m: 0.0005 },
            1200.0,law,&t,&mut || false)?;
    Ok((thermal,beam))
}
pub fn bindings(law: &ThermoMaxwell) -> Vec<Binding<'_>> {
    (4..8).map(|id| Binding { cell: CellId(id),thermal_material: MaterialId(3),mechanical_law: law }).collect()
}
pub fn build(law: &ThermoMaxwell) -> Result<GrowingRoad<'static,'_>,Error> {
    build_with_drawing(law, Drawing::Slack)
}
pub fn build_with_drawing(law: &ThermoMaxwell, drawing: Drawing)
    -> Result<GrowingRoad<'static,'_>,Error> {
    let (thermal,beam) = initial(law,schedule())?;
    let mut cfg=config();
    cfg.drawing=drawing;
    GrowingRoad::new(thermal,beam,BODY,bindings(law),cfg,|| Ok(()))
}
