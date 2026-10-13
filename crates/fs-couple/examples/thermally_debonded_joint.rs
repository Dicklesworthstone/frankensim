//! Cooling-induced joint failure: the opening is solved, not prescribed.
//! Synthetic uniaxial coupon/restraint; no 3D peel or material calibration claim.
use std::collections::BTreeMap;
use std::error::Error;
use fs_conduction::transient::deposition::{Boundary,CellGeometry,CellId,DepositedCell,
    DepositionBudget,DepositionError,DepositionSimulation,Event,MaterialId,StepConfig,
    ThermalLaw,ThermalLink,TimedEvent};
use fs_conduction::transient::enthalpy::deposition::fracture::{ContactPolicy,FractureCooling,Interface};
use fs_material::cohesive::BilinearCohesive;
use fs_material::thermovisco::{Limits,ThermalKnot,ThermoMaxwell};
use fs_solid::cohesive::thermal_joint::{Budget,Config,Joint,Member};
use fs_solver::NewtonKrylovConfig;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self,h: f64) -> Result<f64,DepositionError> {
        if !(0.0..=440_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
fn main() -> Result<(),Box<dyn Error>> {
    let material=ThermoMaxwell::new(1e9,vec![(1e9,3.0)],vec![
        ThermalKnot { temperature_k: 280.0,modulus_scale: 1.0,relaxation_scale: 1.0,free_strain: -0.002 },
        ThermalKnot { temperature_k: 520.0,modulus_scale: 1.0,relaxation_scale: 1.0,free_strain: 0.022 },
    ],Limits { max_strain: 0.04,max_strain_rate_s: 1.0 })?;
    let cohesive=BilinearCohesive::new(1e10,2e5,10.0,2e10,0.01)?;
    let mut joint=Joint::new(0.0,vec![Member::new(1,0.01,1e-4,&material,500.0)?],7,1e-4,&cohesive,
        Budget { max_members: 1,max_branch_states: 1,max_krylov_basis_entries: 64 })?;
    let config=Config { equilibrium: fs_solid::cohesive::Config {
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-10,relative_tolerance: 1e-12,
            linear_restart: 4,forcing_maximum: 0.01,..NewtonKrylovConfig::default() },
        max_newton_iterations: 80,max_operator_evaluations: 2000,force_tolerance_n: 1e-7,
        local_energy_tolerance_j_m2: 1e-7,total_cohesive_energy_tolerance_j: 1e-9,
        maximum_displacement_m: 0.005 },extension_tolerance_m: 1e-11,energy_tolerance_j: 1e-7 };
    let heat=Sensible; let mut events=Vec::new();
    for (id,t) in [(1,500.0),(2,300.0)] {
        events.push(TimedEvent { time_s: 0.0,event: Event::Birth(DepositedCell {
            id: CellId(id),segment_id: id,body_id: id,material_id: MaterialId(1),mass_kg: 0.001,
            specific_h_j_kg: 2000.0*(t-300.0),geometry: CellGeometry {
                start_m: [0.0,0.0,id as f64*0.001],end_m: [0.01,0.0,id as f64*0.001],
                volume_m3: 1e-6,surface_area_m2: 0.001 },
            boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 100.0,reservoir: None,source_w: 0.0 },
        }) });
    }
    let link=ThermalLink { id: 7,a: CellId(1),b: CellId(2),area_m2: 1e-4,conductance_w_k: 0.02 };
    events.push(TimedEvent { time_s: 0.0,event: Event::Contact(link) });
    let thermal_config=StepConfig { max_dt_s: 0.01,energy_tolerance_j: 1e-7 };
    let thermal=DepositionSimulation::new(events,&[(MaterialId(1),&heat)],DepositionBudget {
        max_cells: 2,max_links: 1,max_events: 3,max_substeps: 1000,
    })?.advance_to(0.0,thermal_config,|| Ok(()))?.0;
    let mut coupled=FractureCooling::new(thermal,vec![Interface { reference: link,point: joint.bond(),closed: true,
        policy: ContactPolicy { open_gap_m: 1e-5,reclose_gap_m: 0.0,
            reclosed_conductance_w_k: 0.005,maximum_slip_m: 1e-4 } }],1,|| Ok(()))?;
    let initial_energy=coupled.thermal().energy_j(); let mut external=0.0;
    let mut opened=false; let mut fracture_work=0.0;
    println!("synthetic restrained axial coupon; columns=time_s,temperature_k,force_n,opening_m,damage,thermal_contact");
    for step in 1..=120 {
        let time=0.5*f64::from(step);
        let (heat_candidate,heat_report)=coupled.advance_to(time,thermal_config,|| Ok(()))?;
        let temperature=heat_candidate.thermal().temperature_k(CellId(1))?;
        let (joint_candidate,mechanics)=joint.advance_controlled(time,0.0,
            &BTreeMap::from([(1,temperature)]),config,&mut || false)?;
        let (candidate,contact)=heat_candidate.synchronize(time,&[(7,joint_candidate.bond())],
            config.equilibrium.local_energy_tolerance_j_m2,|| Ok(()))?;
        // Commit both only after thermal, mechanical and contact gates succeed.
        joint=joint_candidate; coupled=candidate;
        external+=heat_report.external_heat_j; fracture_work+=contact.fracture_dissipation_j;
        opened|=contact.opened>0;
        assert!(mechanics.energy_residual_j.abs()<config.energy_tolerance_j);
        if step%10==0 || contact.opened>0 {
            println!("{time},{temperature:.8},{:.8},{:.10},{:.8},{}",mechanics.force_n,
                mechanics.opening_m,joint.bond().damage(),coupled.interfaces()[&7].closed);
        }
    }
    assert!(opened && joint.bond().failed() && !coupled.interfaces()[&7].closed);
    assert!((fracture_work-1e-4*cohesive.fracture_energy_j_m2()).abs()<1e-9);
    assert!((coupled.thermal().energy_j()-initial_energy-external).abs()<1e-6);
    println!("fracture_j={fracture_work}; prescribed total extension=0 m; endpoint-staggered coupling, no fracture heat");
    Ok(())
}
