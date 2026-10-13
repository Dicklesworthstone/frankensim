//! A hot second road reheats an existing road; then an explicit strength map
//! supplies a new cohesive pull test. All data are synthetic, not calibrated.
use std::error::Error;
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionError, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_conduction::transient::enthalpy::deposition::bonding::{BondBinding, BondBudget, BondedDeposition};
use fs_material::healing::{MobilityKnot, MobilityTable};
use fs_material::healing::strength::{StrengthKnot, StrengthMap};

struct Polymer;
impl ThermalLaw for Polymer {
    fn temperature_k(&self,h: f64) -> Result<f64,DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
fn birth(id: u64,time_s: f64,temperature: f64,htc: f64) -> TimedEvent {
    TimedEvent { time_s,event: Event::Birth(DepositedCell { id: CellId(id),segment_id: id,
        body_id: 1,material_id: MaterialId(1),mass_kg: 0.001,specific_h_j_kg: (temperature-300.0)*2000.0,
        geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],end_m: [0.01,0.0,id as f64*0.001],
            volume_m3: 1e-6,surface_area_m2: 0.0012 },
        boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: htc,source_w: 0.0,reservoir: None } }) }
}
fn main() -> Result<(),Box<dyn Error>> {
    let thermal_law = Polymer;
    let mobility = MobilityTable::new("synthetic-mobility-v1".into(),vec![
        MobilityKnot { temperature_k: 300.0,rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 400.0,rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 500.0,rate_per_s: 1.0 },
        MobilityKnot { temperature_k: 600.0,rate_per_s: 3.0 },
    ])?;
    let mapping = StrengthMap::new("synthetic-strength-not-measured".into(),"synthetic-mobility-v1".into(),vec![
        StrengthKnot { dose: 0.0,peak_traction_pa: 0.0,fracture_energy_j_m2: 0.0 },
        StrengthKnot { dose: 20.0,peak_traction_pa: 2e5,fracture_energy_j_m2: 100.0 },
    ],1e9,2e9,0.01)?;
    for (case,htc) in [("slow-cooling",10.0),("fast-cooling",100.0)] {
        let schedule = vec![birth(1,0.0,350.0,htc),birth(2,1.0,550.0,htc),TimedEvent {
            time_s: 1.0,event: Event::Contact(ThermalLink { id: 7,a: CellId(1),b: CellId(2),
                area_m2: 0.0002,conductance_w_k: 0.2 }) }];
        let state = BondedDeposition::new(schedule,&[(MaterialId(1),&thermal_law)],DepositionBudget {
            max_cells: 2,max_links: 1,max_events: 3,max_substeps: 10000 },vec![BondBinding {
                interface_id: 7,a: CellId(1),b: CellId(2),law: &mobility,temperature_weight_a: 0.5 }],
            BondBudget { max_interfaces: 1,max_updates: 10000 })?;
        let config = StepConfig { max_dt_s: 0.125,energy_tolerance_j: 1e-8 };
        let (before,_,_) = state.advance_to(1.0,config,|| Ok(()))?;
        let first_temperature = before.thermal().temperature_k(CellId(1))?;
        let (state,audit,updates) = before.advance_to(10.0,config,|| Ok(()))?;
        assert!(state.thermal().temperature_k(CellId(1))? > first_temperature);
        let handoff = state.freeze_bond(7,&mapping)?;
        let law = handoff.bond.law.ok_or("no developed bond under declared map")?;
        // Displacement-controlled fresh mechanical stage, NOT a reset of damage.
        let mut point = law.bonded_point();
        let mut work = 0.0;
        let mut peak_force: f64 = 0.0;
        for i in 0..=100 {
            let jump = if i == 100 { law.failure_separation_m() }
                else { law.failure_separation_m()*f64::from(i)/100.0 };
            let step = point.trial([jump,0.0,0.0],1e-8)?;
            work += handoff.area_m2*step.energy.work;
            peak_force = peak_force.max(handoff.area_m2*step.traction_pa[0]);
            point = step.point;
        }
        assert!(point.failed());
        assert!((work-handoff.area_m2*law.fracture_energy_j_m2()).abs() < 1e-8);
        println!("case={case}; dose={:.9}; area_m2={:.6e}; peak_pull_force_n={peak_force:.9}; fracture_work_j={work:.9}; thermal_defect_j={:.6e}; bond_updates={updates}; synthetic-only",
            handoff.bond.dose,handoff.area_m2,audit.residual_j);
    }
    Ok(())
}
