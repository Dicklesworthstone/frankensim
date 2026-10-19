//! Spatially moving absorbed heat, not a calibrated laser or manufacturing model.
//! cargo run -p fs-conduction --example scanned_melting -- 2 100
//! Arguments: absorbed power [W] and latent heat [J/kg]. Output: per-node CSV.
use std::error::Error;
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory,ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun,EnthalpyHistoryMaterials,EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::moving::{MovingHeatSource,MovingHeatBudget};
use fs_conduction::bc::history::march::moving::motion::ScanMotionLimit;
use fs_conduction::fixtures::box_grid;
use fs_conduction::transient::enthalpy::{EnthalpyBudget,EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh,ConductivityModel,ScalarField,ThermalBoundaryBuilder,COMPONENT_POWER_DIMS};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve,EnthalpyPhaseKnot};
use fs_qty::{QtyAny,Dims};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal,Interp};
use fs_solver::NewtonKrylovConfig;

fn run(cx:&Cx<'_>,power:f64,latent:f64)->Result<(),Box<dyn Error>> {
    let (complex,positions)=box_grid([8,2,2],[1.,0.2,0.2]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let history=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder(),&Environment::earth_lab(),Vec::new(),0.)?;
    let mut identity=[93;32];identity[..8].copy_from_slice(&latent.to_le_bytes());
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash(identity),
        [(0.,300.,0.),(100.,310.,0.),(100.+latent,310.,1.),(900.+latent,390.,1.)]
        .into_iter().map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let material=ConductivityModel::isotropic_declared(0.02)?;
    let materials=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    let assignment=vec![0;mesh.element_count()];
    let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&history,material:&material,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials {table:&materials,of_element:&assignment,phase_conductivity:None},
        &vec![80.;mesh.vertex_count()],EnthalpyHistoryBudget {spatial:EnthalpyBudget {
            max_vertices:mesh.vertex_count(),max_elements:mesh.element_count()},max_steps:128})?;
    let length=Dims([1,0,0,0,0,0]);
    let source=MovingHeatSource::new([TimeSignal::Ramp {t_start:0.,t_end:1.,
        from:QtyAny::new(0.125,length),to:QtyAny::new(0.875,length)},
        TimeSignal::Constant(QtyAny::new(0.1,length)),TimeSignal::Constant(QtyAny::new(0.1,length))],
        [0.24,0.16,0.16],TimeSignal::Table {times:vec![0.,1.],values:vec![power,0.],
            dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,MovingHeatBudget {
            max_vertices:mesh.vertex_count(),max_elements:mesh.element_count(),min_active_vertices:4,
            power_relative_tolerance:1e-10})?;
    let mut run=ScanningEnthalpyRun::new(cx,base,vec![source],1,None)?;
    run.set_motion_limit(cx,ScanMotionLimit {max_normalized_travel:0.05,max_chebyshev_coefficients:0})?;
    let config=EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-10,
        relative_tolerance:1e-12,linear_restart:40,max_linear_cycles:32,forcing_minimum:1e-12,
        forcing_maximum:1e-3,..Default::default()},max_newton_iterations:128,energy_tolerance_j:1e-7};
    println!("time_s,node,x_m,y_m,z_m,source_x_m,absorbed_W,temperature_K,enthalpy_J_kg,liquid_fraction,storage_change_J,energy_defect_J");
    while run.conduction().time_s()<1.2 {
        let interval=run.advance_event_step(cx,1.2,0.02,config)?;
        let step=interval.report.solution();let source=&interval.report.deposition[0];
        for (node,p) in mesh.positions().iter().enumerate() {
            println!("{},{node},{},{},{},{},{},{},{},{},{},{}",interval.end_time_s,p[0],p[1],p[2],
                source.center_m[0],source.absorbed_power_w,step.temperature[node],
                step.specific_enthalpy_j_kg[node],step.liquid_mass_fraction[node],
                step.stored_energy_change_j,step.energy_residual_j);
        }
    }
    Ok(())
}
fn main()->Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()!=2 {return Err("usage: scanned_melting ABSORBED_W LATENT_J_KG".into());}
    let power:f64=args[0].parse()?;let latent:f64=args[1].parse()?;
    if !power.is_finite() || power<=0. || !latent.is_finite() || latent<=0. {
        return Err("absorbed power and latent heat must be finite and positive".into());
    }
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|run(&Cx::new(&gate,arena,
        StreamKey {seed:93,kernel_id:893,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),power,latent))
}
