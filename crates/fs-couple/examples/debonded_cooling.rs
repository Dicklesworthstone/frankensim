//! A real cohesive equilibrium drives opening/reclosure in the thermal network.
//! Fixed geometry/material data and first-order staggered endpoints; synthetic only.
use std::error::Error;
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_conduction::transient::enthalpy::deposition::fracture::{ContactPolicy, FractureCooling, Interface};
use fs_material::cohesive::BilinearCohesive;
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

struct Heat;
impl ThermalLaw for Heat {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
struct Actuator;
impl LinearOp for Actuator {
    fn n(&self) -> usize { 1 }
    fn apply(&self, x: &[f64], y: &mut [f64]) { y[0] = 1e6*x[0]; }
}
fn birth(id: u64, temperature: f64) -> TimedEvent {
    TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
        mass_kg: 0.001, specific_h_j_kg: 2000.0*(temperature-300.0),
        geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],
            end_m: [0.01,0.0,id as f64*0.001], volume_m3: 1e-6,surface_area_m2: 0.001 },
        boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 20.0,source_w: 0.0,reservoir: None },
    }) }
}
fn main() -> Result<(), Box<dyn Error>> {
    let heat = Heat;
    let law = BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01)?;
    let link = ThermalLink { id: 7,a: CellId(1),b: CellId(2),area_m2: 0.0002,conductance_w_k: 0.2 };
    let cfg = StepConfig { max_dt_s: 0.03125,energy_tolerance_j: 1e-8 };
    let thermal = DepositionSimulation::new(vec![birth(1,500.0),birth(2,300.0),
        TimedEvent { time_s: 0.0,event: Event::Contact(link) }],&[(MaterialId(1),&heat)],
        DepositionBudget { max_cells: 2,max_links: 1,max_events: 3,max_substeps: 10000 })?
        .advance_to(0.0,cfg,|| Ok(()))?.0;
    let initial_energy = thermal.energy_j();
    let mut cooling = FractureCooling::new(thermal,vec![Interface {
        reference: link,point: law.bonded_point(),closed: true,policy: ContactPolicy {
            open_gap_m: 0.0001,reclose_gap_m: 0.0,reclosed_conductance_w_k: 0.05,maximum_slip_m: 0.002 },
    }],1,|| Ok(()))?;
    let actuator = Actuator;
    let solid = CohesiveSystem::new(&actuator,vec![CohesivePatch::new(7,link.area_m2,&law,
        vec![JumpTerm { dof: 0,local: [1.0,0.0,0.0] }])?],Budget {
            max_dofs: 1,max_patches: 1,max_jump_terms: 1,max_krylov_basis_entries: 64 })?;
    let mut mechanical = solid.bonded_state();
    let solve = Config { newton: NewtonKrylovConfig { linear_restart: 1,
        absolute_tolerance: 1e-9,relative_tolerance: 1e-13,forcing_maximum: 0.01,
        ..NewtonKrylovConfig::default() },max_newton_iterations: 80,max_operator_evaluations: 2000,
        force_tolerance_n: 1e-7,local_energy_tolerance_j_m2: 1e-8,
        total_cohesive_energy_tolerance_j: 1e-9,maximum_displacement_m: 0.005 };
    let (mut opened,mut reclosed,mut heat_j,mut fracture_j) = (0,0,0.0,0.0);
    println!("time_s,actuator_m,jump_m,reaction_n,damage,closed,hot_k,cold_k,exposed_m2");
    for i in 1..=16 {
        let t = f64::from(i)*0.125;
        let imposed = if i <= 8 { f64::from(i)*0.00015 }
            else { 0.0012-f64::from(i-8)*0.000175 };
        // Stage all three operations. Commit nothing if any later operation refuses.
        let (thermal_trial,audit) = cooling.advance_to(t,cfg,|| Ok(()))?;
        let proposal = mechanical.solve_controlled(&[1e6*imposed],&[[0.0;3]],solve,&mut || false)?;
        let point = proposal.state.points()[0];
        let (next,change) = thermal_trial.synchronize(t,&[(7,point)],1e-8,|| Ok(()))?;
        assert!((change.fracture_dissipation_j-proposal.cohesive_energy_j.dissipated).abs() < 1e-9);
        cooling = next;
        mechanical = proposal.state;
        heat_j += audit.external_heat_j;
        fracture_j += change.fracture_dissipation_j;
        opened += change.opened;
        reclosed += change.reclosed;
        let reaction = 1e6*(imposed-mechanical.displacements_m()[0]);
        assert!((reaction-link.area_m2*point.traction_pa()[0]).abs() < solve.force_tolerance_n);
        println!("{t},{imposed:.8e},{:.8e},{reaction:.8e},{},{},{},{},{}",point.jump_m()[0],point.damage(),
            cooling.interfaces()[&7].closed,cooling.thermal().temperature_k(CellId(1))?,
            cooling.thermal().temperature_k(CellId(2))?,cooling.thermal().exposed_area_m2(CellId(1))?);
    }
    assert_eq!(opened,1);
    assert_eq!(reclosed,1);
    assert!(mechanical.points()[0].failed());
    assert_eq!(cooling.thermal().links()[&7].conductance_w_k,0.05);
    assert!((cooling.thermal().energy_j()-initial_energy-heat_j).abs() < 1e-8);
    assert!((fracture_j-link.area_m2*law.fracture_energy_j_m2()).abs() < 1e-9);
    println!("fracture_j={fracture_j}; external_heat_j={heat_j}; fracture not injected as thermal heat; synthetic-only");
    Ok(())
}
